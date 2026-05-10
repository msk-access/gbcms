"""Tests for v4.2.0 diagnostic flag computation (_compute_diagnostics).

Covers:
1.  ZERO_ALT flag fires when ad=0 on a PASS variant
2.  PARTIAL_DOMINANT flag fires when partial_alt > ad
3.  MNP_SPARSE_DISC(n/m) flag fires for MNPs with ≤50% discriminating positions
4.  HIGH_N_FRACTION(f) flag fires when n_count/dp > 0.05
5.  Multi-flag combination (ZERO_ALT + MNP_SPARSE_DISC)
6.  FAIL variants never get diagnostics
7.  Clean PASS variant has empty gbcms_diagnostic
8.  Multi-value gbcms_status (PASS;MULTI_ALLELIC) still gets diagnostics
9.  Parametric flag values are correctly formatted
10. PARTIAL_DOMINANT does NOT fire when partial_alt == ad
11. MNP_SPARSE_DISC does NOT fire for SNPs or indels
"""

import types

import pytest

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


# ── Test 1: ZERO_ALT ────────────────────────────────────────────────────


def test_zero_alt_fires_when_ad_is_zero():
    """ZERO_ALT should fire when ad=0 on a PASS variant."""
    pv = _mock_prepared_variant()
    counts = _mock_counts(ad=0, any_alt=0)
    Pipeline._compute_diagnostics([pv], [counts])
    assert "ZERO_ALT" in pv.gbcms_diagnostic


# ── Test 2: PARTIAL_DOMINANT ─────────────────────────────────────────────


def test_partial_dominant_fires_when_partial_exceeds_ad():
    """PARTIAL_DOMINANT should fire when partial_alt > ad."""
    pv = _mock_prepared_variant()
    counts = _mock_counts(ad=2, partial_alt=5, any_alt=7)
    Pipeline._compute_diagnostics([pv], [counts])
    assert "PARTIAL_DOMINANT" in pv.gbcms_diagnostic


# ── Test 3: MNP_SPARSE_DISC ─────────────────────────────────────────────


def test_mnp_sparse_disc_fires_for_sparse_mnp():
    """MNP_SPARSE_DISC should fire for 5bp MNP with 1 discriminating position."""
    # AAAAA → AATAA → 1/5 = 0.20, which is ≤ 0.50
    pv = _mock_prepared_variant(ref_allele="AAAAA", alt_allele="AATAA", variant_type="DNP")
    counts = _mock_counts(ad=0)
    Pipeline._compute_diagnostics([pv], [counts])
    assert "MNP_SPARSE_DISC(1/5)" in pv.gbcms_diagnostic


def test_mnp_sparse_disc_does_not_fire_for_all_disc():
    """MNP_SPARSE_DISC should NOT fire when all positions are discriminating."""
    # AT → GC → 2/2 = 1.0, which is > 0.50
    pv = _mock_prepared_variant(ref_allele="AT", alt_allele="GC", variant_type="DNP")
    counts = _mock_counts(ad=5)
    Pipeline._compute_diagnostics([pv], [counts])
    assert "MNP_SPARSE_DISC" not in pv.gbcms_diagnostic


# ── Test 4: HIGH_N_FRACTION ──────────────────────────────────────────────


def test_high_n_fraction_fires_above_threshold():
    """HIGH_N_FRACTION should fire when n_count/dp > 0.05."""
    pv = _mock_prepared_variant()
    counts = _mock_counts(dp=100, n_count=8, ad=5)
    Pipeline._compute_diagnostics([pv], [counts])
    assert "HIGH_N_FRACTION(0.08)" in pv.gbcms_diagnostic


def test_high_n_fraction_does_not_fire_below_threshold():
    """HIGH_N_FRACTION should NOT fire when n_count/dp ≤ 0.05."""
    pv = _mock_prepared_variant()
    counts = _mock_counts(dp=100, n_count=3, ad=5)
    Pipeline._compute_diagnostics([pv], [counts])
    assert "HIGH_N_FRACTION" not in pv.gbcms_diagnostic


# ── Test 5: Multi-flag combination ───────────────────────────────────────


def test_multi_flag_combination():
    """Multiple flags should combine with semicolons."""
    pv = _mock_prepared_variant(ref_allele="AAAAA", alt_allele="AATAA", variant_type="DNP")
    counts = _mock_counts(ad=0, partial_alt=0, any_alt=0)
    Pipeline._compute_diagnostics([pv], [counts])
    assert "ZERO_ALT" in pv.gbcms_diagnostic
    assert "MNP_SPARSE_DISC(1/5)" in pv.gbcms_diagnostic
    assert ";" in pv.gbcms_diagnostic


# ── Test 6: FAIL variants excluded ──────────────────────────────────────


def test_fail_variants_get_no_diagnostics():
    """FAIL variants should never receive diagnostic flags."""
    pv = _mock_prepared_variant(gbcms_status="FAIL_REF_MISMATCH")
    counts = _mock_counts(ad=0, n_count=50, dp=100)
    Pipeline._compute_diagnostics([pv], [counts])
    assert pv.gbcms_diagnostic == ""


# ── Test 7: Clean PASS ──────────────────────────────────────────────────


def test_clean_pass_has_empty_diagnostic():
    """Clean PASS variant with normal counts has empty gbcms_diagnostic."""
    pv = _mock_prepared_variant()
    counts = _mock_counts(ad=10, partial_alt=0, n_count=1, dp=100)
    Pipeline._compute_diagnostics([pv], [counts])
    assert pv.gbcms_diagnostic == ""


# ── Test 8: Multi-value status gets diagnostics ─────────────────────────


def test_multi_value_status_gets_diagnostics():
    """PASS;MULTI_ALLELIC status should still receive diagnostic flags."""
    pv = _mock_prepared_variant(gbcms_status="PASS;MULTI_ALLELIC")
    counts = _mock_counts(ad=0)
    Pipeline._compute_diagnostics([pv], [counts])
    assert "ZERO_ALT" in pv.gbcms_diagnostic


# ── Test 9: Parametric formatting ────────────────────────────────────────


def test_high_n_fraction_formatting():
    """HIGH_N_FRACTION should format fraction to 2 decimal places."""
    pv = _mock_prepared_variant()
    counts = _mock_counts(dp=200, n_count=15, ad=5)
    Pipeline._compute_diagnostics([pv], [counts])
    # 15/200 = 0.075 → HIGH_N_FRACTION(0.07) with :.2f rounding
    assert "HIGH_N_FRACTION(0.07)" in pv.gbcms_diagnostic


# ── Test 10: PARTIAL_DOMINANT boundary ───────────────────────────────────


def test_partial_dominant_does_not_fire_when_equal():
    """PARTIAL_DOMINANT should NOT fire when partial_alt == ad."""
    pv = _mock_prepared_variant()
    counts = _mock_counts(ad=3, partial_alt=3, any_alt=6)
    Pipeline._compute_diagnostics([pv], [counts])
    assert "PARTIAL_DOMINANT" not in pv.gbcms_diagnostic


# ── Test 11: MNP_SPARSE_DISC not for SNPs/indels ────────────────────────


def test_mnp_sparse_disc_does_not_fire_for_snp():
    """MNP_SPARSE_DISC should NOT fire for single-base SNPs."""
    pv = _mock_prepared_variant(ref_allele="A", alt_allele="T")
    counts = _mock_counts(ad=0)
    Pipeline._compute_diagnostics([pv], [counts])
    assert "MNP_SPARSE_DISC" not in pv.gbcms_diagnostic


def test_mnp_sparse_disc_does_not_fire_for_indel():
    """MNP_SPARSE_DISC should NOT fire for indels (different ref/alt lengths)."""
    pv = _mock_prepared_variant(ref_allele="AAA", alt_allele="A", variant_type="DELETION")
    counts = _mock_counts(ad=0)
    Pipeline._compute_diagnostics([pv], [counts])
    assert "MNP_SPARSE_DISC" not in pv.gbcms_diagnostic


# ── Test 12: PARTIAL_DOMINANT with nearby evidence ───────────────────────


def test_partial_alt_indel_nonzero():
    """After v4.2.0, INDEL with nearby evidence should have partial_alt > 0."""
    pv = _mock_prepared_variant(ref_allele="AAA", alt_allele="A", variant_type="DELETION")
    counts = _mock_counts(ad=0, partial_alt=3, any_alt=3)
    Pipeline._compute_diagnostics([pv], [counts])
    assert "PARTIAL_DOMINANT" in pv.gbcms_diagnostic
    assert "ZERO_ALT" in pv.gbcms_diagnostic
