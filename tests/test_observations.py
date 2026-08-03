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

from helpers import build_bam, make_read

from gbcms._rs import Variant, count_bam_binned, count_bam_binned_observations

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


def test_decomposed_variant_emits_one_set_of_rows(tmp_path):
    """A decomposed variant runs the classifier twice; only the winner's rows survive.

    Counts keep just the winning form, so a double-emit would leave parity green while
    the export carried two contradictory allele calls per molecule. The `dpf` leg of the
    invariant is what catches it.
    """
    bam = build_bam(
        tmp_path,
        [_read(f"alt{i}", ALT_BASE) for i in range(3)]
        + [_read(f"ref{i}", REF_BASE) for i in range(2)],
    )
    variants = [_variant()]
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
    _assert_reconciles(counts, observations)
    # one row per molecule, not two
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
