"""Per-molecule observation export (`count_bam_binned_observations`).

The counting path is untouched by the export, so binned↔legacy parity cannot detect a
broken export — parity compares `BaseCounts` only. The binding test here is therefore
load-bearing: for every variant the emitted rows must reconcile with the **fragment-level**
counts (`adf`/`rdf`/`dpf`).

It must be fragment-level, not read-level: a read excluded by the multi-allelic sibling
guard still contributes REF evidence to its fragment, so it counts in `rdf` and emits a REF
row even though read-level `rd` excluded it. Asserting against `rd`/`ad`/`dp` would fail
spuriously and mask real bugs.
"""

from collections import Counter

import pysam
import pytest
from helpers import build_bam, make_read

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


def _kwargs(variants, **overrides):
    """Standard counting args; `overrides` tweak individual ones per test."""
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
    }
    kwargs.update(overrides)
    return kwargs


def _observe(bam, variants, **overrides):
    return count_bam_binned_observations(bam, variants, **_kwargs(variants, **overrides))


def _count(bam, variants, **overrides):
    return count_bam_binned(bam, variants, **_kwargs(variants, **overrides))


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


def test_counts_are_identical_with_and_without_observations(tmp_path):
    """`count_bam_binned` must be unaffected — same core, observations off."""
    bam = build_bam(
        tmp_path,
        [_read(f"alt{i}", ALT_BASE) for i in range(3)]
        + [_read(f"ref{i}", REF_BASE) for i in range(2)],
    )
    variants = [_variant()]
    counts_only = _count(bam, variants)
    counts, observations = _observe(bam, variants)

    fields = ("dp", "rd", "ad", "dpf", "rdf", "adf", "dp_fwd", "dp_rev")
    for f in fields:
        assert getattr(counts_only[0], f) == getattr(counts[0], f), f
    assert observations, "observations were requested but none were emitted"


def test_observations_reconcile_with_fragment_counts(tmp_path):
    bam = build_bam(
        tmp_path,
        [_read(f"alt{i}", ALT_BASE) for i in range(4)]
        + [_read(f"ref{i}", REF_BASE) for i in range(3)],
    )
    counts, observations = _observe(bam, [_variant()])
    _assert_reconciles(counts, observations)
    assert Counter(o.allele for o in observations) == {1: 4, 0: 3}


# ── 2. decomposed variants must not double-emit ─────────────────────────────────────


def test_decomposed_variant_emits_the_winning_forms_rows(tmp_path):
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
    counts, observations = count_bam_binned_observations(
        bam,
        variants,
        decomposed=decomposed,
        min_mapq=20,
        min_baseq=20,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
        fragment_qual_threshold=10,
        sibling_variants=[[]],
    )
    assert counts[0].used_decomposed, "fixture failed to exercise the decomposed-wins branch"
    assert counts[0].adf == 3
    _assert_reconciles(counts, observations)
    # the winner's rows: 3 ALT (not the original form's 0)
    assert Counter(o.allele for o in observations)[1] == 3
    hashes = [o.molecule_hash for o in observations]
    assert len(hashes) == len(set(hashes)), "a molecule was emitted more than once"


# ── 3. multi-allelic siblings ───────────────────────────────────────────────────────


def test_sibling_variants_still_reconcile(tmp_path):
    """With a sibling, reads excluded read-level still carry REF into their fragment."""
    bam = build_bam(
        tmp_path,
        [_read(f"alt{i}", ALT_BASE) for i in range(2)]
        + [_read(f"ref{i}", REF_BASE) for i in range(2)]
        + [_read("third", "G")],
    )
    variants = [_variant()]
    counts, observations = count_bam_binned_observations(
        bam,
        variants,
        decomposed=[None],
        min_mapq=20,
        min_baseq=20,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
        fragment_qual_threshold=10,
        sibling_variants=[[Variant("chr1", POS, REF_BASE, "G", "SNP")]],
    )
    _assert_reconciles(counts, observations)


# ── 4. determinism ──────────────────────────────────────────────────────────────────


def test_observation_order_is_deterministic(tmp_path):
    """Rayon bin order and HashMap iteration are both nondeterministic; sorting erases both."""
    reads = [_read(f"alt{i}", ALT_BASE) for i in range(10)]
    reads += [_read(f"ref{i}", REF_BASE) for i in range(10)]
    bam = build_bam(tmp_path, reads)
    variants = [_variant()]

    _, first = _observe(bam, variants, threads=4)
    for _ in range(3):
        _, again = _observe(bam, variants, threads=4)
        assert [(o.variant_index, o.molecule_hash, o.allele) for o in first] == [
            (o.variant_index, o.molecule_hash, o.allele) for o in again
        ]
    keys = [(o.variant_index, o.molecule_hash) for o in first]
    assert keys == sorted(keys), "rows are not sorted by (variant_index, molecule_hash)"


# ── 5. molecule identity is shared across variants (what enables cross-locus linking) ─


def test_same_molecule_shares_one_hash_across_variants(tmp_path):
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
    counts, observations = _observe(bam, variants)
    _assert_reconciles(counts, observations)

    at0 = {o.molecule_hash for o in observations if o.variant_index == 0}
    at1 = {o.molecule_hash for o in observations if o.variant_index == 1}
    assert at0 and at0 == at1, "the same molecules must carry the same hash at both variants"


# ── 6. allele encoding ──────────────────────────────────────────────────────────────


def test_ref_and_alt_are_distinguished(tmp_path):
    """REF is first-class, not 'absence of ALT'."""
    bam = build_bam(tmp_path, [_read("a", ALT_BASE), _read("r", REF_BASE)])
    _, observations = _observe(bam, [_variant()])
    by_allele = {o.allele: o for o in observations}
    assert set(by_allele) == {0, 1}
    assert by_allele[1].best_qual > 0 and by_allele[0].best_qual > 0


def test_n_base_is_reported_as_N_not_OTHER(tmp_path):
    """The N/OTHER split is the point: N marks an ambiguous base, OTHER a third allele.

    In consensus BAMs an N is a *strand-discordant* molecule, which is diagnostic — it must
    not be conflated with a third-allele observation. Collapsing the two (or swapping them)
    leaves every count identical, so only an explicit assertion catches it.
    """
    bam = build_bam(tmp_path, [_read("n", "N"), _read("third", "G"), _read("r", REF_BASE)])
    counts, observations = _observe(bam, [_variant()])
    _assert_reconciles(counts, observations)
    alleles = Counter(o.allele for o in observations)
    assert alleles[ALLELE_N] == 1, f"expected one N row, got {dict(alleles)}"
    assert alleles[ALLELE_OTHER] == 1, f"expected one OTHER row, got {dict(alleles)}"
    assert alleles[ALLELE_REF] == 1
    # neither N nor OTHER has a called allele, so neither carries a quality
    for obs in observations:
        if obs.allele in (ALLELE_N, ALLELE_OTHER):
            assert obs.best_qual == 0


def test_paired_reads_collapse_to_one_molecule(tmp_path):
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
    counts, observations = _observe(bam, [_variant()])
    c = counts[0]
    assert c.ad == 2 and c.adf == 1, f"fixture not read-vs-fragment discriminating: {c.ad}/{c.adf}"
    assert len(observations) == 1, "R1+R2 of one template must yield ONE observation"
    assert observations[0].allele == ALLELE_ALT
    _assert_reconciles(counts, observations)


# ── 7. the public wrapper ───────────────────────────────────────────────────────────


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


# ── 8. Parquet sink (at-scale path) ─────────────────────────────────────────────────


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
    ]
    assert set(table["chrom"]) == {"chr1"} and set(table["pos"]) == {POS}
    assert set(table["ref"]) == {REF_BASE} and set(table["alt"]) == {ALT_BASE}

    from_file = sorted(
        zip(
            table["variant_index"],
            table["molecule_hash"],
            table["allele"],
            table["best_qual"],
            strict=True,
        )
    )
    from_mem = sorted(
        (o.variant_index, o.molecule_hash, o.allele, o.best_qual) for o in mem.observations
    )
    assert from_file == from_mem
