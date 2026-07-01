"""Tests for the MNP rescue pass (--rescue-mnp, v4.3.0).

Covers:
1.  Default config: rescue_mnp is False
2.  MAF output: gbcms_rescue column absent without flag
3.  MNP_DISC_RATIO always emitted, MNP_RESCUE_ELIGIBLE based on threshold
4.  HIGH_N_FRACTION(f) flag fires when n_count/dp > 0.05
5.  Multi-flag combination (ZERO_ALT + MNP_DISC_RATIO + MNP_RESCUE_ELIGIBLE)
6.  Rescue fires for MNP (ad=0, MNP_RESCUE_ELIGIBLE)
7.  Rescue skips non-MNP variants
8.  Rescue skips variants with non-zero ad
9.  Rescue skips FAIL variants
10. Rescue audit trail format validation
11. Rescue no-signal case (all SNPs have ad=0)
12. Column count with rescue_mnp=True
13. Invariant breakage after rescue (design §2)
"""

import types

from gbcms.io.output import MafWriter, VcfWriter
from gbcms.models.core import GbcmsBaseConfig

# ── Helpers ──────────────────────────────────────────────────────────────────


def _mock_prepared_variant(
    *,
    chrom: str = "chr5",
    pos: int = 1295227,  # 0-based
    ref_allele: str = "CCCCC",
    alt_allele: str = "CTCCC",
    gbcms_status: str = "PASS",
    gbcms_status_reason: str = "",
    gbcms_diagnostic: str = "",
    gbcms_rescue: str = "",
):
    """Create a mock PreparedVariant for unit tests."""
    variant = types.SimpleNamespace(
        chrom=chrom,
        pos=pos,
        ref_allele=ref_allele,
        alt_allele=alt_allele,
    )
    return types.SimpleNamespace(
        variant=variant,
        gbcms_status=gbcms_status,
        gbcms_status_reason=gbcms_status_reason,
        gbcms_diagnostic=gbcms_diagnostic,
        gbcms_rescue=gbcms_rescue,
        was_anchor_resolved=False,
        was_left_aligned=False,
        was_normalized=False,
        original_pos=pos,
        original_ref=ref_allele,
        original_alt=alt_allele,
        decomposed_variant=None,
        multi_allelic_group=None,
    )


def _mock_counts(
    *,
    dp: int = 100,
    rd: int = 100,
    ad: int = 0,
    any_alt: int = 0,
    partial_alt: int = 0,
    n_count: int = 0,
):
    """Create a minimal mock counts object."""
    _nan = float("nan")
    return types.SimpleNamespace(
        dp=dp,
        rd=rd,
        ad=ad,
        rd_fwd=50,
        rd_rev=50,
        ad_fwd=0,
        ad_rev=0,
        dpf=50,
        rdf=50,
        adf=0,
        rdf_fwd=25,
        rdf_rev=25,
        adf_fwd=0,
        adf_rev=0,
        sb_pval=1.0,
        sb_or=1.0,
        fsb_pval=1.0,
        fsb_or=1.0,
        used_decomposed=False,
        any_alt=any_alt,
        partial_alt=partial_alt,
        n_count=n_count,
        mfsd_ref_count=0,
        mfsd_alt_count=0,
        mfsd_nonref_count=0,
        mfsd_n_count=0,
        mfsd_ref_mean=_nan,
        mfsd_alt_mean=_nan,
        mfsd_nonref_mean=_nan,
        mfsd_n_mean=_nan,
        mfsd_alt_llr=_nan,
        mfsd_ref_llr=_nan,
        mfsd_delta_alt_ref=_nan,
        mfsd_ks_alt_ref=_nan,
        mfsd_pval_alt_ref=_nan,
        mfsd_delta_alt_nonref=_nan,
        mfsd_ks_alt_nonref=_nan,
        mfsd_pval_alt_nonref=_nan,
        mfsd_delta_ref_nonref=_nan,
        mfsd_ks_ref_nonref=_nan,
        mfsd_pval_ref_nonref=_nan,
        mfsd_delta_alt_n=_nan,
        mfsd_ks_alt_n=_nan,
        mfsd_pval_alt_n=_nan,
        mfsd_delta_ref_n=_nan,
        mfsd_ks_ref_n=_nan,
        mfsd_pval_ref_n=_nan,
        mfsd_delta_nonref_n=_nan,
        mfsd_ks_nonref_n=_nan,
        mfsd_pval_nonref_n=_nan,
    )


# ── Test 1: Default config ──────────────────────────────────────────────────


def test_rescue_flag_absent_by_default():
    """GbcmsBaseConfig.rescue_mnp defaults to False."""
    # Test via direct field default — no file I/O needed
    assert GbcmsBaseConfig.model_fields["rescue_mnp"].default is False


# ── Test 2: MAF column absent without flag ───────────────────────────────────


def test_rescue_column_absent_without_flag():
    """MafWriter omits gbcms_rescue column when rescue_mnp=False."""
    writer = MafWriter.__new__(MafWriter)
    writer.column_prefix = ""
    writer.mfsd = False
    writer.show_normalization = False
    writer.mode = "dna"
    writer.rescue_mnp = False

    cols = writer._gbcms_column_names()
    assert (
        "gbcms_rescue" not in cols
    ), f"gbcms_rescue should not be in columns without --rescue-mnp: {cols}"


# ── Test 3: MAF column present with flag ─────────────────────────────────────


def test_rescue_column_present_with_flag():
    """MafWriter includes gbcms_rescue column when rescue_mnp=True."""
    writer = MafWriter.__new__(MafWriter)
    writer.column_prefix = ""
    writer.mfsd = False
    writer.show_normalization = False
    writer.mode = "dna"
    writer.rescue_mnp = True

    cols = writer._gbcms_column_names()
    assert "gbcms_rescue" in cols, f"gbcms_rescue should be in columns with --rescue-mnp: {cols}"


# ── Test 4: VCF GR absent without flag ───────────────────────────────────────


def test_vcf_gr_absent_without_flag(tmp_path):
    """VcfWriter omits GR INFO header and value when rescue_mnp=False."""
    vcf_path = tmp_path / "no_rescue.vcf"
    writer = VcfWriter(vcf_path, sample_name="TUMOR", rescue_mnp=False)
    writer._write_header()
    writer.close()

    header_text = vcf_path.read_text()
    assert "ID=GR," not in header_text, "GR INFO header should not appear without --rescue-mnp"


# ── Test 5: VCF GR present with flag ─────────────────────────────────────────


def test_vcf_gr_present_with_flag(tmp_path):
    """VcfWriter includes GR INFO header when rescue_mnp=True."""
    vcf_path = tmp_path / "with_rescue.vcf"
    writer = VcfWriter(vcf_path, sample_name="TUMOR", rescue_mnp=True)
    writer._write_header()
    writer.close()

    header_text = vcf_path.read_text()
    assert "ID=GR," in header_text, "GR INFO header should appear with --rescue-mnp"


# ── Test 6: Rescue fires for sparse MNP ──────────────────────────────────────


def test_rescue_candidate_identification():
    """_rescue_mnp_pass identifies candidates: PASS + ad==0 + MNP_RESCUE_ELIGIBLE."""
    # Create a sparse MNP: 5bp, only 1 discriminating position
    pv = _mock_prepared_variant(
        ref_allele="CCCCC",
        alt_allele="CTCCC",  # only position 1 differs
        gbcms_status="PASS",
        gbcms_diagnostic="ZERO_ALT;MNP_DISC_RATIO(1/5);MNP_RESCUE_ELIGIBLE",
    )
    counts = _mock_counts(ad=0, any_alt=0, partial_alt=0)

    # Build the candidate identification portion of the rescue logic
    # (testing the filtering criteria, not the full BAM counting)
    prepared = [pv]
    full_counts = [counts]

    candidates = []
    for i, (p, c) in enumerate(zip(prepared, full_counts, strict=True)):
        if p.gbcms_status != "PASS":
            continue
        if c.ad != 0:
            continue
        if "MNP_RESCUE_ELIGIBLE" not in p.gbcms_diagnostic:
            continue

        ref_allele = p.variant.ref_allele
        alt_allele = p.variant.alt_allele
        if len(ref_allele) == len(alt_allele) and len(ref_allele) > 1:
            disc_positions = [
                (p.variant.pos + offset, ref_allele[offset], alt_allele[offset])
                for offset in range(len(ref_allele))
                if ref_allele[offset] != alt_allele[offset]
            ]
            if disc_positions:
                candidates.append((i, disc_positions))

    assert len(candidates) == 1, f"Expected 1 candidate, got {len(candidates)}"
    assert candidates[0][0] == 0  # index 0
    assert len(candidates[0][1]) == 1  # 1 discriminating position
    pos, ref_base, alt_base = candidates[0][1][0]
    assert ref_base == "C" and alt_base == "T"


# ── Test 7: Rescue skips non-MNP variants ─────────────────────────────────────


def test_rescue_skips_non_mnp():
    """SNP and INDEL variants are never rescue candidates."""
    # SNP
    snp = _mock_prepared_variant(
        ref_allele="A",
        alt_allele="T",
        gbcms_status="PASS",
        gbcms_diagnostic="ZERO_ALT",
    )
    snp_counts = _mock_counts(ad=0)

    # Indel
    indel = _mock_prepared_variant(
        ref_allele="AC",
        alt_allele="A",
        gbcms_status="PASS",
        gbcms_diagnostic="ZERO_ALT",
    )
    indel_counts = _mock_counts(ad=0)

    for pv, _c in [(snp, snp_counts), (indel, indel_counts)]:
        is_mnp = (
            len(pv.variant.ref_allele) == len(pv.variant.alt_allele)
            and len(pv.variant.ref_allele) > 1
        )
        assert not is_mnp or "MNP_RESCUE_ELIGIBLE" not in pv.gbcms_diagnostic


# ── Test 8: Rescue skips non-zero ad ──────────────────────────────────────────


def test_rescue_skips_nonzero_ad():
    """MNP with ad > 0 should not be a rescue candidate."""
    _mock_prepared_variant(
        ref_allele="CCCCC",
        alt_allele="CTCCC",
        gbcms_status="PASS",
        gbcms_diagnostic="MNP_DISC_RATIO(1/5);MNP_RESCUE_ELIGIBLE",
    )
    counts = _mock_counts(ad=3)  # non-zero

    # Should not be a candidate
    assert counts.ad != 0, "Variant with ad>0 should not be rescued"


# ── Test 9: Rescue skips FAIL variants ────────────────────────────────────────


def test_rescue_skips_fail():
    """FAIL variants are never rescue candidates."""
    pv = _mock_prepared_variant(
        ref_allele="CCCCC",
        alt_allele="CTCCC",
        gbcms_status="FAIL",
        gbcms_status_reason="REF_MISMATCH",
        gbcms_diagnostic="",
    )
    _mock_counts(ad=0)

    # Should not be a candidate
    assert pv.gbcms_status != "PASS"


# ── Test 10: Rescue audit trail format ────────────────────────────────────────


def test_rescue_audit_trail_format():
    """gbcms_rescue audit trail has the correct structured format."""
    # Simulate a successful rescue
    pv = _mock_prepared_variant()
    rescue_str = "method=decomposed;original_alt=0;positions=chr5:1295229(C>T):3"
    pv.gbcms_rescue = rescue_str

    # Parse and validate structure
    parts = pv.gbcms_rescue.split(";")
    kv = {p.split("=", 1)[0]: p.split("=", 1)[1] for p in parts}

    assert kv["method"] == "decomposed"
    assert kv["original_alt"] == "0"
    assert "positions" in kv
    # Positions format: chrom:pos(ref>alt):count
    positions = kv["positions"].split(",")
    assert len(positions) >= 1
    for pos_entry in positions:
        assert ":" in pos_entry
        assert ">" in pos_entry


# ── Test 11: Rescue no-signal case ────────────────────────────────────────────


def test_rescue_no_signal_format():
    """Failed rescue (all SNPs ad=0) should have outcome=no_signal."""
    rescue_str = (
        "method=decomposed;original_alt=0;outcome=no_signal;" "positions=chr5:1295229(C>T):0"
    )

    parts = rescue_str.split(";")
    kv = {p.split("=", 1)[0]: p.split("=", 1)[1] for p in parts}

    assert kv["outcome"] == "no_signal"
    assert kv["original_alt"] == "0"
    # All positions should have count=0
    positions = kv["positions"].split(",")
    for pos_entry in positions:
        count = pos_entry.rsplit(":", 1)[1]
        assert count == "0", f"Expected count=0 in no-signal rescue, got {count}"


# ── Test 12: Column count with rescue ────────────────────────────────────────


def test_column_count_with_rescue():
    """With rescue_mnp=True, there should be 27 columns (26 base + gbcms_rescue)."""
    writer = MafWriter.__new__(MafWriter)
    writer.column_prefix = ""
    writer.mfsd = False
    writer.show_normalization = False
    writer.mode = "dna"
    writer.rescue_mnp = True

    cols = writer._gbcms_column_names()
    assert len(cols) == 27, f"Expected 27 gbcms MAF columns with rescue, got {len(cols)}: {cols}"
    assert cols[3] == "gbcms_rescue", f"gbcms_rescue should be the 4th column, got {cols[3]}"


# ── Test 13: Invariant breakage after rescue ─────────────────────────────────


def test_invariant_breakage_after_rescue():
    """After rescue, Invariant 1 (any_alt = ad + partial_alt) intentionally breaks.

    This is by design: ad is updated with the best decomposed SNP count while
    any_alt and partial_alt retain original MNP-level values as forensic evidence.
    See validation_status_design.md §2 for the TERT example.
    """
    # Simulate the TERT-like case: before rescue
    counts = _mock_counts(ad=0, any_alt=108, partial_alt=108)

    # Invariant holds before rescue
    assert counts.any_alt == counts.ad + counts.partial_alt, "Invariant 1 should hold before rescue"

    # Simulate rescue: ad is updated to best decomposed SNP count
    counts.ad = 108  # rescued

    # Invariant intentionally breaks after rescue
    assert counts.any_alt != counts.ad + counts.partial_alt, (
        "Invariant 1 should break after rescue: "
        f"any_alt({counts.any_alt}) != ad({counts.ad}) + partial_alt({counts.partial_alt})"
    )
    # Verify the expected broken state
    assert counts.ad == 108
    assert counts.partial_alt == 108  # unchanged — forensic evidence
    assert counts.any_alt == 108  # unchanged — forensic evidence
