"""Tests for v4.2.0/v5.0.0 diagnostic flag computation (_compute_diagnostics).

Covers:
1.  ZERO_ALT flag fires when ad=0 on a PASS variant
2.  PARTIAL_DOMINANT flag fires when partial_alt > ad
3.  MNP_DISC_RATIO(n/m) flag always emitted for MNPs
4.  MNP_RESCUE_ELIGIBLE fires based on configurable threshold
5.  HIGH_N_FRACTION(f) flag fires when n_count/dp > 0.05
6.  Multi-flag combination (ZERO_ALT + MNP_DISC_RATIO + MNP_RESCUE_ELIGIBLE)
7.  FAIL variants never get diagnostics
8.  Clean PASS variant has empty gbcms_diagnostic
9.  Multi-value gbcms_status (PASS;MULTI_ALLELIC) still gets diagnostics
10. Parametric flag values are correctly formatted
11. PARTIAL_DOMINANT does NOT fire when partial_alt == ad
12. MNP_DISC_RATIO does NOT fire for SNPs or indels
13. MNP_RESCUE_ELIGIBLE respects threshold=0.50 (conservative mode)
"""

import types

from gbcms.pipeline import Pipeline

# ── Helpers ──────────────────────────────────────────────────────────────


def _mock_prepared_variant(
    *,
    chrom: str = "chr1",
    pos: int = 100,
    ref_allele: str = "A",
    alt_allele: str = "T",
    variant_type: str = "SNP",
    gbcms_status: str = "PASS",
    gbcms_status_reason: str = "",
):
    """Create a mock PreparedVariant with the minimum fields needed."""
    return types.SimpleNamespace(
        variant=types.SimpleNamespace(
            chrom=chrom,
            pos=pos,
            ref_allele=ref_allele,
            alt_allele=alt_allele,
            variant_type=variant_type,
        ),
        gbcms_status=gbcms_status,
        gbcms_status_reason=gbcms_status_reason,
        gbcms_diagnostic="",
        gbcms_rescue="",
    )


def _mock_counts(
    *,
    dp: int = 100,
    ad: int = 5,
    partial_alt: int = 0,
    any_alt: int = 5,
    n_count: int = 0,
):
    """Create a mock BaseCounts with diagnostic-relevant fields."""
    return types.SimpleNamespace(
        dp=dp,
        ad=ad,
        partial_alt=partial_alt,
        any_alt=any_alt,
        n_count=n_count,
    )


def _make_pipeline(rescue_mnp_threshold: float = 1.0):
    """Create a minimal Pipeline-like object for calling _compute_diagnostics.

    We only need self.config.rescue_mnp_threshold; all other config fields
    are irrelevant to diagnostic computation.
    """
    pipeline = object.__new__(Pipeline)
    pipeline.config = types.SimpleNamespace(rescue_mnp_threshold=rescue_mnp_threshold)
    return pipeline


# ── Test 1: ZERO_ALT ────────────────────────────────────────────────────


def test_zero_alt_fires_when_ad_is_zero():
    """ZERO_ALT should fire when ad=0 on a PASS variant."""
    pv = _mock_prepared_variant()
    counts = _mock_counts(ad=0, any_alt=0)
    _make_pipeline()._compute_diagnostics([pv], [counts])
    assert "ZERO_ALT" in pv.gbcms_diagnostic


# ── Test 2: PARTIAL_DOMINANT ─────────────────────────────────────────────


def test_partial_dominant_fires_when_partial_exceeds_ad():
    """PARTIAL_DOMINANT should fire when partial_alt > ad."""
    pv = _mock_prepared_variant()
    counts = _mock_counts(ad=2, partial_alt=5, any_alt=7)
    _make_pipeline()._compute_diagnostics([pv], [counts])
    assert "PARTIAL_DOMINANT" in pv.gbcms_diagnostic


# ── Test 3: MNP_DISC_RATIO (always emitted) ─────────────────────────────


def test_mnp_disc_ratio_fires_for_sparse_mnp():
    """MNP_DISC_RATIO should always fire for MNPs (here, 1/5 = sparse)."""
    pv = _mock_prepared_variant(ref_allele="AAAAA", alt_allele="AATAA", variant_type="DNP")
    counts = _mock_counts(ad=0)
    _make_pipeline()._compute_diagnostics([pv], [counts])
    assert "MNP_DISC_RATIO(1/5)" in pv.gbcms_diagnostic
    # With default threshold=1.0, should also be rescue eligible
    assert "MNP_RESCUE_ELIGIBLE" in pv.gbcms_diagnostic


def test_mnp_disc_ratio_fires_for_dense_mnp():
    """MNP_DISC_RATIO should fire even for fully dense DNPs (2/2)."""
    pv = _mock_prepared_variant(ref_allele="AT", alt_allele="GC", variant_type="DNP")
    counts = _mock_counts(ad=5)
    _make_pipeline()._compute_diagnostics([pv], [counts])
    assert "MNP_DISC_RATIO(2/2)" in pv.gbcms_diagnostic
    # With default threshold=1.0, even dense MNPs are rescue eligible
    assert "MNP_RESCUE_ELIGIBLE" in pv.gbcms_diagnostic


# ── Test 4: MNP_RESCUE_ELIGIBLE with conservative threshold ─────────────


def test_rescue_eligible_respects_threshold():
    """With threshold=0.50, dense DNPs (2/2=100%) should NOT be rescue eligible."""
    pv = _mock_prepared_variant(ref_allele="AT", alt_allele="GC", variant_type="DNP")
    counts = _mock_counts(ad=5)
    _make_pipeline(rescue_mnp_threshold=0.50)._compute_diagnostics([pv], [counts])
    assert "MNP_DISC_RATIO(2/2)" in pv.gbcms_diagnostic
    assert "MNP_RESCUE_ELIGIBLE" not in pv.gbcms_diagnostic


def test_rescue_eligible_sparse_with_conservative_threshold():
    """With threshold=0.50, sparse MNPs (1/5=20%) should still be rescue eligible."""
    pv = _mock_prepared_variant(ref_allele="AAAAA", alt_allele="AATAA", variant_type="DNP")
    counts = _mock_counts(ad=0)
    _make_pipeline(rescue_mnp_threshold=0.50)._compute_diagnostics([pv], [counts])
    assert "MNP_DISC_RATIO(1/5)" in pv.gbcms_diagnostic
    assert "MNP_RESCUE_ELIGIBLE" in pv.gbcms_diagnostic


# ── Test 5: HIGH_N_FRACTION ──────────────────────────────────────────────


def test_high_n_fraction_fires_above_threshold():
    """HIGH_N_FRACTION should fire when n_count/dp > 0.05."""
    pv = _mock_prepared_variant()
    counts = _mock_counts(dp=100, n_count=8, ad=5)
    _make_pipeline()._compute_diagnostics([pv], [counts])
    assert "HIGH_N_FRACTION(0.08)" in pv.gbcms_diagnostic


def test_high_n_fraction_does_not_fire_below_threshold():
    """HIGH_N_FRACTION should NOT fire when n_count/dp ≤ 0.05."""
    pv = _mock_prepared_variant()
    counts = _mock_counts(dp=100, n_count=3, ad=5)
    _make_pipeline()._compute_diagnostics([pv], [counts])
    assert "HIGH_N_FRACTION" not in pv.gbcms_diagnostic


# ── Test 6: Multi-flag combination ───────────────────────────────────────


def test_multi_flag_combination():
    """Multiple flags should combine with semicolons."""
    pv = _mock_prepared_variant(ref_allele="AAAAA", alt_allele="AATAA", variant_type="DNP")
    counts = _mock_counts(ad=0, partial_alt=0, any_alt=0)
    _make_pipeline()._compute_diagnostics([pv], [counts])
    assert "ZERO_ALT" in pv.gbcms_diagnostic
    assert "MNP_DISC_RATIO(1/5)" in pv.gbcms_diagnostic
    assert "MNP_RESCUE_ELIGIBLE" in pv.gbcms_diagnostic
    assert ";" in pv.gbcms_diagnostic


# ── Test 7: FAIL variants excluded ──────────────────────────────────────


def test_fail_variants_get_no_diagnostics():
    """FAIL variants should never receive diagnostic flags."""
    pv = _mock_prepared_variant(gbcms_status="FAIL_REF_MISMATCH")
    counts = _mock_counts(ad=0, n_count=50, dp=100)
    _make_pipeline()._compute_diagnostics([pv], [counts])
    assert pv.gbcms_diagnostic == ""


# ── Test 8: Clean PASS ──────────────────────────────────────────────────


def test_clean_pass_has_empty_diagnostic():
    """Clean PASS variant with normal counts has empty gbcms_diagnostic."""
    pv = _mock_prepared_variant()
    counts = _mock_counts(ad=10, partial_alt=0, n_count=1, dp=100)
    _make_pipeline()._compute_diagnostics([pv], [counts])
    assert pv.gbcms_diagnostic == ""


# ── Test 9: Multi-value status gets diagnostics ─────────────────────────


def test_multi_value_status_gets_diagnostics():
    """A PASS variant with reason tags (e.g. MULTI_ALLELIC) still gets diagnostics."""
    pv = _mock_prepared_variant(gbcms_status="PASS", gbcms_status_reason="MULTI_ALLELIC")
    counts = _mock_counts(ad=0)
    _make_pipeline()._compute_diagnostics([pv], [counts])
    assert "ZERO_ALT" in pv.gbcms_diagnostic


# ── Test 10: Parametric formatting ────────────────────────────────────────


def test_high_n_fraction_formatting():
    """HIGH_N_FRACTION should format fraction to 2 decimal places."""
    pv = _mock_prepared_variant()
    counts = _mock_counts(dp=200, n_count=15, ad=5)
    _make_pipeline()._compute_diagnostics([pv], [counts])
    # 15/200 = 0.075 → HIGH_N_FRACTION(0.07) with :.2f rounding
    assert "HIGH_N_FRACTION(0.07)" in pv.gbcms_diagnostic


# ── Test 11: PARTIAL_DOMINANT boundary ───────────────────────────────────


def test_partial_dominant_does_not_fire_when_equal():
    """PARTIAL_DOMINANT should NOT fire when partial_alt == ad."""
    pv = _mock_prepared_variant()
    counts = _mock_counts(ad=3, partial_alt=3, any_alt=6)
    _make_pipeline()._compute_diagnostics([pv], [counts])
    assert "PARTIAL_DOMINANT" not in pv.gbcms_diagnostic


# ── Test 12: MNP_DISC_RATIO not for SNPs/indels ────────────────────────


def test_mnp_disc_ratio_does_not_fire_for_snp():
    """MNP_DISC_RATIO should NOT fire for single-base SNPs."""
    pv = _mock_prepared_variant(ref_allele="A", alt_allele="T")
    counts = _mock_counts(ad=0)
    _make_pipeline()._compute_diagnostics([pv], [counts])
    assert "MNP_DISC_RATIO" not in pv.gbcms_diagnostic


def test_mnp_disc_ratio_does_not_fire_for_indel():
    """MNP_DISC_RATIO should NOT fire for indels (different ref/alt lengths)."""
    pv = _mock_prepared_variant(ref_allele="AAA", alt_allele="A", variant_type="DELETION")
    counts = _mock_counts(ad=0)
    _make_pipeline()._compute_diagnostics([pv], [counts])
    assert "MNP_DISC_RATIO" not in pv.gbcms_diagnostic


# ── Test 13: PARTIAL_DOMINANT with nearby evidence ───────────────────────


def test_partial_alt_indel_nonzero():
    """After v4.2.0, INDEL with nearby evidence should have partial_alt > 0."""
    pv = _mock_prepared_variant(ref_allele="AAA", alt_allele="A", variant_type="DELETION")
    counts = _mock_counts(ad=0, partial_alt=3, any_alt=3)
    _make_pipeline()._compute_diagnostics([pv], [counts])
    assert "PARTIAL_DOMINANT" in pv.gbcms_diagnostic
    assert "ZERO_ALT" in pv.gbcms_diagnostic


def test_non_discriminating_locus_flag():
    """NON_DISCRIMINATING_LOCUS fires when the Rust counts flag it — a sibling combo
    reconstructs REF, so REF/ALT are indistinguishable and reads tie to NEITHER."""
    pv = _mock_prepared_variant()
    counts = _mock_counts(ad=0, any_alt=0)
    counts.non_discriminating_locus = True
    _make_pipeline()._compute_diagnostics([pv], [counts])
    assert "NON_DISCRIMINATING_LOCUS" in pv.gbcms_diagnostic
