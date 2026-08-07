"""Per-molecule observation export (`count_bam_binned_observations`).

The counting path is untouched by the export, so binned↔legacy parity cannot detect a
broken export — parity compares `BaseCounts` only. The binding test here is therefore
load-bearing: for every variant the emitted rows must reconcile with the **fragment-level**
counts (`adf`/`rdf`/`dpf`).

It must be fragment-level, not read-level: a read excluded by the multi-allelic sibling
guard still contributes REF evidence to its fragment, so it counts in `rdf` and emits a REF
row even though read-level `rd` excluded it. Asserting against `rd`/`ad`/`dp` would fail
spuriously and mask real bugs.

Every `_rs`-level test below runs under **both** alignment backends. The reconciliation
invariant is only assertable at this level — `observe_molecules` returns rows without the
counts to reconcile against — and the two backends are not interchangeable: on real ACCESS
indels they disagree on ~5% of molecule calls. Pinned to one, the invariant would be proven
for a classifier half the callers never reach. `pairhmm` is the default at every layer; `sw`
remains selectable, so both are exercised here.
"""

from collections import Counter
from pathlib import Path

import pysam
import pytest
from helpers import build_bam, make_read
from pydantic import ValidationError

from gbcms._rs import Variant, count_bam_binned, count_bam_binned_observations
from gbcms.observations import ALLELE_ALT, ALLELE_N, ALLELE_OTHER, ALLELE_REF

REF_BASE = "A"
ALT_BASE = "T"
POS = 100  # 0-based
SEQ_LEN = 20


def _read(name, base_at_pos, start=90, flag=0):
    """A 20bp read spanning POS, carrying `base_at_pos` at the variant position."""
    seq = list(REF_BASE * SEQ_LEN)
    seq[POS - start] = base_at_pos
    return make_read(name, "".join(seq), start, ((0, SEQ_LEN),), flag=flag)


def _variant(ref=REF_BASE, alt=ALT_BASE):
    return Variant("chr1", POS, ref, alt, "SNP")


def _py_variant(ref=REF_BASE, alt=ALT_BASE):
    """The public-model twin of `_variant()`. `observe_molecules` takes these, not `_rs` ones."""
    from gbcms.models.core import Variant as PyVariant
    from gbcms.models.core import VariantType

    return PyVariant(chrom="chr1", pos=POS, ref=ref, alt=alt, variant_type=VariantType.SNP)


@pytest.fixture(params=["sw", "pairhmm"])
def backend(request):
    """Both alignment backends: `pairhmm` is the default everywhere, `sw` the explicit opt-in."""
    return request.param


def _kwargs(variants, backend, **overrides):
    """Standard counting args; `overrides` tweak individual ones per test.

    `backend` is positional and required so every test states which classifier it scored
    itself against rather than inheriting one — the same reason the Rust layer no longer
    has a `Default` impl for `AlignmentBackend`.
    """
    kwargs = {
        "decomposed": [None] * len(variants),
        "min_mapq": 20,
        "min_baseq": 20,
        "filter_duplicates": True,
        "filter_secondary": True,
        "filter_supplementary": True,
        "filter_qc_failed": False,
        "filter_improper_pair": False,
        "filter_indel": False,
        "threads": 1,
        "fragment_qual_threshold": 10,
        "sibling_variants": [[] for _ in variants],
        "alignment_backend": backend,
    }
    kwargs.update(overrides)
    return kwargs


def _observe(bam, variants, backend, **overrides):
    return count_bam_binned_observations(bam, variants, **_kwargs(variants, backend, **overrides))


def _count(bam, variants, backend, **overrides):
    return count_bam_binned(bam, variants, **_kwargs(variants, backend, **overrides))


def _assert_reconciles(counts, observations):
    """The load-bearing invariant: rows reconcile with FRAGMENT-level counts."""
    by_variant = {}
    for obs in observations:
        by_variant.setdefault(obs.variant_index, []).append(obs)
    for vi, c in enumerate(counts):
        rows = by_variant.get(vi, [])
        alleles = Counter(o.allele for o in rows)
        assert alleles[1] == c.adf, f"variant {vi}: #obs(ALT)={alleles[1]} != adf={c.adf}"
        assert alleles[0] == c.rdf, f"variant {vi}: #obs(REF)={alleles[0]} != rdf={c.rdf}"
        assert len(rows) == c.dpf, f"variant {vi}: rows={len(rows)} != dpf={c.dpf}"


# ── 1. counts unchanged + the reconciliation invariant ──────────────────────────────


def test_counts_are_identical_with_and_without_observations(tmp_path, backend):
    """`count_bam_binned` must be unaffected — same core, observations off."""
    bam = build_bam(
        tmp_path,
        [_read(f"alt{i}", ALT_BASE) for i in range(3)]
        + [_read(f"ref{i}", REF_BASE) for i in range(2)],
    )
    variants = [_variant()]
    counts_only = _count(bam, variants, backend)
    counts, observations = _observe(bam, variants, backend)

    fields = ("dp", "rd", "ad", "dpf", "rdf", "adf", "dp_fwd", "dp_rev")
    for f in fields:
        assert getattr(counts_only[0], f) == getattr(counts[0], f), f
    assert observations, "observations were requested but none were emitted"


def test_observations_reconcile_with_fragment_counts(tmp_path, backend):
    bam = build_bam(
        tmp_path,
        [_read(f"alt{i}", ALT_BASE) for i in range(4)]
        + [_read(f"ref{i}", REF_BASE) for i in range(3)],
    )
    counts, observations = _observe(bam, [_variant()], backend)
    _assert_reconciles(counts, observations)
    assert Counter(o.allele for o in observations) == {1: 4, 0: 3}


# ── 2. decomposed variants must not double-emit ─────────────────────────────────────


def test_decomposed_variant_emits_the_winning_forms_rows(tmp_path, backend):
    """The decomposed form WINS here, so the arbitration branch is actually executed.

    A decomposed variant runs the classifier twice and the higher-`ad` form wins. The
    observations must follow that same arbitration: pairing the winning counts with the
    *losing* rows would leave every count — and therefore parity — untouched while the
    export carried the wrong allele for every molecule. Asserting `used_decomposed` is
    what keeps this test honest; with an identical decomposed form the branch never runs
    and the test would pass vacuously.
    """
    bam = build_bam(
        tmp_path,
        [_read(f"alt{i}", ALT_BASE) for i in range(3)]
        + [_read(f"ref{i}", REF_BASE) for i in range(2)],
    )
    # original has NO support; the decomposed form carries all 3 ALT molecules
    variants = [Variant("chr1", POS, REF_BASE, "G", "SNP")]
    decomposed = [Variant("chr1", POS, REF_BASE, ALT_BASE, "SNP")]
    counts, observations = _observe(bam, variants, backend, decomposed=decomposed)
    assert counts[0].used_decomposed, "fixture failed to exercise the decomposed-wins branch"
    assert counts[0].adf == 3
    _assert_reconciles(counts, observations)
    # the winner's rows: 3 ALT (not the original form's 0)
    assert Counter(o.allele for o in observations)[1] == 3
    hashes = [o.molecule_hash for o in observations]
    assert len(hashes) == len(set(hashes)), "a molecule was emitted more than once"


# ── 3. multi-allelic siblings ───────────────────────────────────────────────────────


def test_sibling_variants_still_reconcile(tmp_path, backend):
    """With a sibling, reads excluded read-level still carry REF into their fragment."""
    bam = build_bam(
        tmp_path,
        [_read(f"alt{i}", ALT_BASE) for i in range(2)]
        + [_read(f"ref{i}", REF_BASE) for i in range(2)]
        + [_read("third", "G")],
    )
    variants = [_variant()]
    counts, observations = _observe(
        bam,
        variants,
        backend,
        sibling_variants=[[Variant("chr1", POS, REF_BASE, "G", "SNP")]],
    )
    _assert_reconciles(counts, observations)


# ── 4. determinism ──────────────────────────────────────────────────────────────────


def test_observation_order_is_deterministic(tmp_path, backend):
    """Rayon bin order and HashMap iteration are both nondeterministic; sorting erases both."""
    reads = [_read(f"alt{i}", ALT_BASE) for i in range(10)]
    reads += [_read(f"ref{i}", REF_BASE) for i in range(10)]
    bam = build_bam(tmp_path, reads)
    variants = [_variant()]

    _, first = _observe(bam, variants, backend, threads=4)
    for _ in range(3):
        _, again = _observe(bam, variants, backend, threads=4)
        assert [(o.variant_index, o.molecule_hash, o.allele) for o in first] == [
            (o.variant_index, o.molecule_hash, o.allele) for o in again
        ]
    keys = [(o.variant_index, o.molecule_hash) for o in first]
    assert keys == sorted(keys), "rows are not sorted by (variant_index, molecule_hash)"


# ── 5. molecule identity is shared across variants (what enables cross-locus linking) ─


def test_same_molecule_shares_one_hash_across_variants(tmp_path, backend):
    """A read spanning two variants yields one molecule_hash at both — the join key."""
    pos2 = POS + 5
    reads = []
    for i in range(3):
        seq = list(REF_BASE * SEQ_LEN)
        seq[POS - 90] = ALT_BASE
        seq[pos2 - 90] = ALT_BASE
        reads.append(make_read(f"m{i}", "".join(seq), 90, ((0, SEQ_LEN),)))
    bam = build_bam(tmp_path, reads)
    variants = [_variant(), Variant("chr1", pos2, REF_BASE, ALT_BASE, "SNP")]
    counts, observations = _observe(bam, variants, backend)
    _assert_reconciles(counts, observations)

    at0 = {o.molecule_hash for o in observations if o.variant_index == 0}
    at1 = {o.molecule_hash for o in observations if o.variant_index == 1}
    assert at0 and at0 == at1, "the same molecules must carry the same hash at both variants"


# ── 6. allele encoding ──────────────────────────────────────────────────────────────


def test_ref_and_alt_are_distinguished(tmp_path, backend):
    """REF is first-class, not 'absence of ALT'."""
    bam = build_bam(tmp_path, [_read("a", ALT_BASE), _read("r", REF_BASE)])
    _, observations = _observe(bam, [_variant()], backend)
    by_allele = {o.allele: o for o in observations}
    assert set(by_allele) == {0, 1}
    assert by_allele[1].best_qual > 0 and by_allele[0].best_qual > 0


def test_n_base_is_reported_as_N_not_OTHER(tmp_path, backend):
    """The N/OTHER split is the point: N marks an ambiguous base, OTHER a third allele.

    In consensus BAMs an N is a *strand-discordant* molecule, which is diagnostic — it must
    not be conflated with a third-allele observation. Collapsing the two (or swapping them)
    leaves every count identical, so only an explicit assertion catches it.
    """
    bam = build_bam(tmp_path, [_read("n", "N"), _read("third", "G"), _read("r", REF_BASE)])
    counts, observations = _observe(bam, [_variant()], backend)
    _assert_reconciles(counts, observations)
    alleles = Counter(o.allele for o in observations)
    assert alleles[ALLELE_N] == 1, f"expected one N row, got {dict(alleles)}"
    assert alleles[ALLELE_OTHER] == 1, f"expected one OTHER row, got {dict(alleles)}"
    assert alleles[ALLELE_REF] == 1
    # neither N nor OTHER has a called allele, so neither carries a quality
    for obs in observations:
        if obs.allele in (ALLELE_N, ALLELE_OTHER):
            assert obs.best_qual == 0


def test_paired_reads_collapse_to_one_molecule(tmp_path, backend):
    """Fragment-level, not read-level: R1+R2 of one template are ONE observation.

    This is the whole reason the export is per-molecule and the invariant is asserted
    against adf/rdf/dpf rather than ad/rd/dp — with unpaired fixtures the two are
    indistinguishable, so a read-level export would pass every other test in this file.
    """
    paired, mate = 0x1 | 0x2 | 0x40, 0x1 | 0x2 | 0x80
    bam = build_bam(
        tmp_path,
        [_read("frag", ALT_BASE, flag=paired), _read("frag", ALT_BASE, flag=mate)],
    )
    counts, observations = _observe(bam, [_variant()], backend)
    c = counts[0]
    assert c.ad == 2 and c.adf == 1, f"fixture not read-vs-fragment discriminating: {c.ad}/{c.adf}"
    assert len(observations) == 1, "R1+R2 of one template must yield ONE observation"
    assert observations[0].allele == ALLELE_ALT
    _assert_reconciles(counts, observations)


# ── 7. min_mapq ─────────────────────────────────────────────────────────────────────


def test_min_mapq_reports_the_measured_value_not_the_filter_threshold(tmp_path, backend):
    """Each molecule carries its own MAPQ, not the `--min-mapq` gate it happened to pass.

    The threshold is the reason this field has to exist: it *bounds* mapping confidence
    without measuring it. A consumer weighting evidence by error probability and forced to
    assume the threshold would charge every molecule ~1e-2, swamping the base-quality term
    and flattening exactly the per-molecule differences the weighting exists to express.
    """
    reads = [_read("hi", ALT_BASE), _read("mid", ALT_BASE), _read("lo", REF_BASE)]
    reads[0].mapping_quality = 60
    reads[1].mapping_quality = 42
    reads[2].mapping_quality = 25  # still above the min_mapq=20 gate below
    bam = build_bam(tmp_path, reads, filename="mapq.bam")

    counts, observations = _observe(bam, [_variant()], backend)
    _assert_reconciles(counts, observations)
    assert sorted(o.min_mapq for o in observations) == [25, 42, 60]


def test_min_mapq_takes_the_worst_read_of_a_fragment(tmp_path, backend):
    """R1+R2 are one molecule, and it is only as trustworthy as its worse-placed read.

    Taking the best (or last-seen) would let a cleanly-placed mate launder an ambiguously
    placed one — confidence manufactured from the pair rather than measured.
    """
    paired, mate = 0x1 | 0x2 | 0x40, 0x1 | 0x2 | 0x80
    r1 = _read("frag", ALT_BASE, flag=paired)
    r2 = _read("frag", ALT_BASE, flag=mate)
    r1.mapping_quality, r2.mapping_quality = 60, 23
    bam = build_bam(tmp_path, [r1, r2], filename="mapqpair.bam")

    counts, observations = _observe(bam, [_variant()], backend)
    assert counts[0].adf == 1, "fixture did not collapse to one fragment"
    assert len(observations) == 1
    assert observations[0].min_mapq == 23, "kept the better read's MAPQ instead of the worse"


def test_min_mapq_is_never_the_unavailable_sentinel_for_an_observed_molecule(tmp_path, backend):
    """255 means "unavailable"; a molecule built from real reads must report a real value.

    Guards the initializer: `min_mapq` starts at 255 so an unobserved fragment cannot
    masquerade as MAPQ 0 (a real value meaning multi-mapping). If the running-minimum update
    were ever skipped — gated behind is_ref/is_alt, say — rows would silently ship 255 and a
    consumer would degrade them as unknown rather than use the evidence.
    """
    bam = build_bam(
        tmp_path,
        [_read("a", ALT_BASE), _read("r", REF_BASE), _read("n", "N"), _read("o", "G")],
        filename="mapqsent.bam",
    )
    _, observations = _observe(bam, [_variant()], backend)
    assert observations, "fixture produced no rows"
    # includes the N and OTHER molecules — they have no called allele but are still placed
    assert all(o.min_mapq == 60 for o in observations), {o.allele: o.min_mapq for o in observations}


# ── 8. indels — the class where the two backends actually disagree ──────────────────


def test_deletion_rows_reconcile_under_both_backends(tmp_path, backend):
    """The reconciliation invariant on an INDEL, not a SNP.

    Every other `_rs`-level fixture in this file is a SNP, so until now the invariant was
    never asserted on the variant class that matters most downstream — indels dominate the
    reversion classes gbcms feeds, and indels are where the deletion-specific paths
    (CIGAR walking, the ambiguity window, decomposition) actually run.

    This particular deletion is unambiguous, and both backends agree on it (5 ALT / 4 REF) —
    a clean deletion is exactly what WFA triage resolves without reaching PairHMM. It is
    therefore not a divergence fixture, and is not claimed as one; the measured 5-of-104
    disagreement lives in repeat-context indels that synthetic reads do not reproduce
    honestly. What this pins down is narrower and still worth pinning: on an indel, under
    *both* classifiers, the exported rows reconcile with the counts beside them.

    The assertion is deliberately not a fixed allele mix. Which molecules a backend calls ALT
    is the classifier's business and the two are allowed to differ; the export's contract is
    only that it reports whatever that decision was.
    """
    from gbcms.core.kernel import CoordinateKernel

    ref = "A" * 100 + "TAG" + "A" * 100  # anchor T at 0-based 100, deleted A at 101

    def read(name, carries_del):
        if carries_del:  # 1bp deletion of the A at 101
            return make_read(name, ref[90:101] + ref[102:120], 90, ((0, 11), (2, 1), (0, 18)))
        return make_read(name, ref[90:120], 90, ((0, 30),))

    bam = build_bam(
        tmp_path,
        [read(f"d{i}", True) for i in range(5)] + [read(f"r{i}", False) for i in range(4)],
        filename="del.bam",
    )
    v = CoordinateKernel.vcf_to_internal("chr1", 101, "TA", "T")
    variants = [Variant(v.chrom, v.pos, v.ref, v.alt, v.variant_type.value)]

    counts, observations = _observe(bam, variants, backend)
    _assert_reconciles(counts, observations)
    assert counts[0].dpf == 9, "fixture lost molecules before classification"
    assert counts[0].adf > 0, f"{backend} found no ALT — fixture is not exercising the deletion"


# ── 9. the public wrapper ───────────────────────────────────────────────────────────


def test_settings_are_configurable_without_fabricating_a_config(tmp_path):
    """The motivating case: change one setting without inventing paths that mean nothing.

    `GbcmsDnaConfig` requires `variant_file`, `bam_files`, `reference_fasta` and `output` —
    a set with **zero** overlap with what this entry point reads. Before the individual
    arguments existed, adjusting a single filter meant constructing all four, so a caller
    had to supply an output directory and a variant path that are never touched. Paths that
    look load-bearing and are not is how the next reader gets misled.
    """
    import gbcms
    from gbcms.models.core import GbcmsDnaConfig, ReadFilters

    with pytest.raises(ValidationError):
        GbcmsDnaConfig()  # the four required-but-unread fields

    bam = build_bam(tmp_path, [_read("a", ALT_BASE)], filename="nocfg.bam")
    result = gbcms.observe_molecules(bam, [_py_variant()], filters=ReadFilters(duplicates=False))
    assert result.n_rows == 1


def test_filters_argument_actually_reaches_the_engine(tmp_path):
    """A forwarded setting must change the result, or it is decoration.

    Uses `duplicates` deliberately. `supplementary` would be the intuitive choice for an
    observation test, but supplementary reads are dropped unconditionally further in
    (msk-access/gbcms#82), so asserting a difference there would fail — and asserting
    *no* difference would pass whether or not the argument was plumbed at all.
    """
    import gbcms
    from gbcms.models.core import ReadFilters

    dup = _read("dup", ALT_BASE)
    dup.flag |= 0x400  # mark duplicate
    bam = build_bam(tmp_path, [_read("keep", ALT_BASE), dup], filename="dupfilt.bam")

    filtered = gbcms.observe_molecules(bam, [_py_variant()], filters=ReadFilters(duplicates=True))
    kept = gbcms.observe_molecules(bam, [_py_variant()], filters=ReadFilters(duplicates=False))
    assert filtered.n_rows == 1, "duplicate was not filtered"
    assert kept.n_rows == 2, "duplicates=False did not reach the engine"


def test_umi_tag_argument_decides_what_counts_as_one_molecule(tmp_path):
    """`umi_tag` is folded into `molecule_hash`, so it changes molecule identity itself."""
    import gbcms

    paired, mate = 0x1 | 0x2 | 0x40, 0x1 | 0x2 | 0x80
    r1 = _read("frag", ALT_BASE, flag=paired)
    r2 = _read("frag", ALT_BASE, flag=mate)
    r1.set_tag("MI", "famA")
    r2.set_tag("MI", "famB")  # same template, deliberately different families
    bam = build_bam(tmp_path, [r1, r2], filename="umitag.bam")

    by_qname = gbcms.observe_molecules(bam, [_py_variant()], umi_tag=None)
    by_umi = gbcms.observe_molecules(bam, [_py_variant()], umi_tag="MI")
    assert by_qname.n_rows == 1, "without a UMI tag the pair is one molecule"
    assert by_umi.n_rows == 2, "umi_tag did not reach the engine"


def _throwaway_config(tmp_path, bam, **over):
    """Build a `GbcmsDnaConfig` purely to have one — every required field is unread here.

    Not just four extra keyword arguments: `variant_file` and `reference_fasta` must name
    files that **exist on disk**, `bam_files` is a dict, and `output` is an `OutputConfig`.
    This helper is the cost the individual arguments remove, kept only so the override
    tests below have a config to override.
    """
    from gbcms.models.core import GbcmsDnaConfig, OutputConfig

    (tmp_path / "unused.txt").write_text("")
    (tmp_path / "unused.fa").write_text(">chr1\nA\n")
    return GbcmsDnaConfig(
        variant_file=tmp_path / "unused.txt",
        bam_files={"s": Path(bam)},
        reference_fasta=tmp_path / "unused.fa",
        output=OutputConfig(directory=tmp_path / "unused_out"),
        **over,
    )


def test_explicit_argument_overrides_config(tmp_path):
    """Both supplied → the individual argument wins, so a pipeline config can be adjusted."""
    import gbcms
    from gbcms.models.core import ReadFilters

    dup = _read("dup", ALT_BASE)
    dup.flag |= 0x400
    bam = build_bam(tmp_path, [_read("keep", ALT_BASE), dup], filename="override.bam")
    cfg = _throwaway_config(tmp_path, bam, filters=ReadFilters(duplicates=True))
    assert gbcms.observe_molecules(bam, [_py_variant()], config=cfg).n_rows == 1
    assert (
        gbcms.observe_molecules(
            bam, [_py_variant()], config=cfg, filters=ReadFilters(duplicates=False)
        ).n_rows
        == 2
    ), "the explicit argument did not override config"


def test_umi_tag_none_overrides_a_config_that_sets_one(tmp_path):
    """`None` is a real choice for `umi_tag`, not an absence — the one nullable argument.

    Every other argument treats `None` as "not supplied" and defers to `config`. If
    `umi_tag` did the same, a caller passing `umi_tag=None` to group by read pair would
    silently inherit the config's tag and group by UMI family instead — changing what
    counts as one molecule, with nothing to indicate it.
    """
    import gbcms

    paired, mate = 0x1 | 0x2 | 0x40, 0x1 | 0x2 | 0x80
    r1 = _read("frag", ALT_BASE, flag=paired)
    r2 = _read("frag", ALT_BASE, flag=mate)
    r1.set_tag("MI", "famA")
    r2.set_tag("MI", "famB")
    bam = build_bam(tmp_path, [r1, r2], filename="umi_none.bam")
    cfg = _throwaway_config(tmp_path, bam, umi_tag="MI")
    assert gbcms.observe_molecules(bam, [_py_variant()], config=cfg).n_rows == 2
    assert (
        gbcms.observe_molecules(bam, [_py_variant()], config=cfg, umi_tag=None).n_rows == 1
    ), "umi_tag=None was treated as 'not supplied' instead of an explicit override"


def test_public_wrapper_returns_observations(tmp_path):
    """`gbcms.observe_molecules` is the supported surface; `_rs` is internal."""
    import gbcms
    from gbcms.models.core import Variant as PyVariant
    from gbcms.models.core import VariantType

    bam = build_bam(
        tmp_path,
        [_read(f"alt{i}", ALT_BASE) for i in range(2)] + [_read("ref0", REF_BASE)],
    )
    result = gbcms.observe_molecules(
        bam,
        [
            PyVariant(
                chrom="chr1", pos=POS, ref=REF_BASE, alt=ALT_BASE, variant_type=VariantType.SNP
            )
        ],
    )
    assert result.n_rows == len(result.observations) == 3
    assert result.path is None
    assert Counter(o.allele for o in result.observations) == {1: 2, 0: 1}


def test_wrapper_forwards_decomposition_so_alt_is_not_lost(tmp_path):
    """Normalization's *decomposed* form must reach the counter, or ALT rows vanish.

    `prepare_variants` can rewrite a complex indel into a decomposed form that is what the
    reads actually carry. An earlier revision kept only `p.variant` and passed
    `decomposed=[None] * n`, so those molecules were exported as OTHER with **zero ALT
    rows** — silently, with PASS status. Regression guard for that.
    """
    import gbcms
    from gbcms.models.core import Variant as PyVariant
    from gbcms.models.core import VariantType

    ref = "A" * 100 + "CCCCCC" + "T" + "A" * 100  # CCCCCC at 0-based 100..105
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + ref + "\n")
    pysam.faidx(str(fasta))

    def read(name, carries_alt):
        seq = list(ref[90:120])
        if carries_alt:
            seq[105 - 90] = "T"  # the decomposed form: CCCCCC -> CCCCCT
        return make_read(name, "".join(seq), 90, ((0, 30),))

    bam = build_bam(
        tmp_path,
        [read(f"alt{i}", True) for i in range(6)] + [read(f"ref{i}", False) for i in range(4)],
        filename="decomp.bam",
    )
    variants = [
        PyVariant(chrom="chr1", pos=100, ref="CCCCCC", alt="T", variant_type=VariantType.DELETION)
    ]
    result = gbcms.observe_molecules(bam, variants, reference_fasta=str(fasta))
    alleles = Counter(o.allele for o in result.observations)
    assert alleles[ALLELE_ALT] == 6, f"ALT molecules lost: {dict(alleles)}"
    assert alleles[ALLELE_REF] == 4
    assert result.variant_status == ["PASS"]


def test_vcf_and_maf_representations_converge(tmp_path):
    """The same deletion, written VCF-style or MAF-style, must observe identically.

    VCF alleles are already anchored, so `is_maf=False` (the default) is correct and needs
    no reference lookup. MAF writes indels with '-' and the anchor base can only be resolved
    against the reference — so `is_maf=True` is required there. Getting that flag wrong
    silently drops the ALT molecules at PASS status, the same failure class as losing the
    decomposed form, which is why `is_maf` is an explicit parameter rather than a guess.
    """
    import gbcms
    from gbcms.core.kernel import CoordinateKernel

    ref = "A" * 100 + "TAG" + "A" * 100  # anchor T at 0-based 100, deleted A at 101
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + ref + "\n")
    pysam.faidx(str(fasta))

    def read(name, carries_del):
        if carries_del:  # 1bp deletion of the A at 101
            return make_read(name, ref[90:101] + ref[102:120], 90, ((0, 11), (2, 1), (0, 18)))
        return make_read(name, ref[90:120], 90, ((0, 30),))

    bam = build_bam(
        tmp_path,
        [read(f"d{i}", True) for i in range(5)] + [read(f"r{i}", False) for i in range(4)],
        filename="vcfmaf.bam",
    )

    # VCF: 1-based POS = anchor, REF = anchor+deleted, ALT = anchor
    v_vcf = CoordinateKernel.vcf_to_internal("chr1", 101, "TA", "T")
    # MAF: 1-based Start = first deleted base, REF = deleted base, ALT = '-'
    v_maf = CoordinateKernel.maf_to_internal("chr1", 102, 102, "A", "-")

    vcf = gbcms.observe_molecules(bam, [v_vcf], reference_fasta=str(fasta), is_maf=False)
    maf = gbcms.observe_molecules(bam, [v_maf], reference_fasta=str(fasta), is_maf=True)

    mix_vcf = Counter(o.allele for o in vcf.observations)
    mix_maf = Counter(o.allele for o in maf.observations)
    assert mix_vcf[ALLELE_ALT] == 5, f"VCF path lost ALT molecules: {dict(mix_vcf)}"
    assert mix_vcf == mix_maf, f"VCF {dict(mix_vcf)} != MAF {dict(mix_maf)}"


# ── 10. Parquet sink (at-scale path) ─────────────────────────────────────────────────


def test_parquet_sink_matches_the_in_memory_rows(tmp_path):
    """Writing to Parquet must carry exactly the rows the in-memory path returns.

    At panel/genome scale the FFI materialization — not the counting — is the bottleneck,
    so the rows are written from Rust and never become Python objects. The two paths must
    therefore agree exactly, or the at-scale path is silently a different answer.
    """
    import gbcms
    from gbcms.models.core import Variant as PyVariant
    from gbcms.models.core import VariantType

    bam = build_bam(
        tmp_path,
        [_read(f"alt{i}", ALT_BASE) for i in range(4)]
        + [_read(f"ref{i}", REF_BASE) for i in range(3)],
        filename="pq.bam",
    )
    variants = [
        PyVariant(chrom="chr1", pos=POS, ref=REF_BASE, alt=ALT_BASE, variant_type=VariantType.SNP)
    ]
    out = tmp_path / "obs.parquet"

    mem = gbcms.observe_molecules(bam, variants)
    written = gbcms.observe_molecules(bam, variants, observations_path=out)

    assert written.path == out and out.exists()
    assert written.observations == [], "rows must not cross the FFI boundary when written"
    assert written.n_rows == mem.n_rows == 7

    pq = pytest.importorskip("pyarrow.parquet")
    table = pq.read_table(out).to_pydict()
    # self-describing: the locus travels with the rows, so the file stands alone
    assert [f.name for f in pq.read_table(out).schema] == [
        "variant_index",
        "chrom",
        "pos",
        "ref",
        "alt",
        "molecule_hash",
        "allele",
        "best_qual",
        "min_mapq",
    ]
    assert set(table["chrom"]) == {"chr1"} and set(table["pos"]) == {POS}
    assert set(table["ref"]) == {REF_BASE} and set(table["alt"]) == {ALT_BASE}

    from_file = sorted(
        zip(
            table["variant_index"],
            table["molecule_hash"],
            table["allele"],
            table["best_qual"],
            table["min_mapq"],
            strict=True,
        )
    )
    from_mem = sorted(
        (o.variant_index, o.molecule_hash, o.allele, o.best_qual, o.min_mapq)
        for o in mem.observations
    )
    assert from_file == from_mem


def test_cli_flag_writes_observations_alongside_counts(tmp_path):
    """`--observations-parquet` must produce the side-file without disturbing counts.

    Every external consumer (Nextflow, ACCESS-Pipeline, CWL) reaches gbcms through the CLI,
    so a capability that exists only in the Python API is unreachable to them — which is
    exactly the audience the at-scale Parquet path was built for. Mirrors --mfsd-parquet.
    """
    pytest.importorskip("pyarrow.parquet")
    from typer.testing import CliRunner

    from gbcms.cli import app

    ref = "A" * 200
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + ref + "\n")
    pysam.faidx(str(fasta))

    bam = build_bam(
        tmp_path,
        [_read(f"alt{i}", ALT_BASE) for i in range(3)]
        + [_read(f"ref{i}", REF_BASE) for i in range(2)],
        filename="cli.bam",
    )
    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr1\t{POS + 1}\t.\t{REF_BASE}\t{ALT_BASE}\t.\tPASS\t.\n"
    )
    out = tmp_path / "out"

    args = ["dna", "-v", str(vcf), "-b", bam, "-f", str(fasta), "-o", str(out)]
    res_off = CliRunner().invoke(app, args)
    assert res_off.exit_code == 0, res_off.output
    assert not list(out.glob("*.observations.parquet")), "written without the flag"

    res_on = CliRunner().invoke(app, [*args, "--observations-parquet"])
    assert res_on.exit_code == 0, res_on.output
    written = list(out.glob("*.observations.parquet"))
    assert len(written) == 1, "flag did not produce the side-file"

    import pyarrow.parquet as pq

    table = pq.read_table(written[0])
    assert table.num_rows == 5  # one row per fragment
    assert [f.name for f in table.schema] == [
        "variant_index",
        "chrom",
        "pos",
        "ref",
        "alt",
        "molecule_hash",
        "allele",
        "best_qual",
        "min_mapq",
    ]
    # the counts output is still produced, unchanged by the flag
    assert list(out.glob("*.vcf")) or list(out.glob("*.maf"))
