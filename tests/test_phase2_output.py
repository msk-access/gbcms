"""Tests for Phase 2: AAD/PAD/NAD (any_alt/partial_alt/n_count) output columns.

Covers:
1. MAF column presence (any_alt, partial_alt, n_count) in DNA mode
2. MAF column values from _populate_gbcms_counts
3. MAF column prefix support (t_any_alt, t_partial_alt, t_n_count)
4. VCF INFO header presence (AAD, PAD, NAD)
5. VCF FORMAT header presence (AAD, PAD, NAD)
6. VCF INFO data values
7. VCF FORMAT data values
8. _zero_counts() includes any_alt=0, partial_alt=0, n_count=0
9. Column count delta vs. Phase 0/1 (should be +3)
"""

import types

import pytest
from helpers import read_maf_output as _read_maf_output

from gbcms.io.output import MafWriter, VcfWriter
from gbcms.models.core import Variant, VariantType
from gbcms.pipeline import _zero_counts

# ── Helpers ───────────────────────────────────────────────────────────────

_nan = float("nan")


def _mock_counts(
    *,
    ad: int = 5,
    any_alt: int = 7,
    partial_alt: int = 2,
    n_count: int = 1,
):
    """Minimal mock BaseCounts with decomposed ALT counting fields.

    Default: ad=5, any_alt=7, partial_alt=2, n_count=1 → invariant holds (7=5+2).
    """
    return types.SimpleNamespace(
        dp=100,
        rd=95,
        ad=ad,
        dp_fwd=50,
        dp_rev=50,
        rd_fwd=48,
        rd_rev=47,
        ad_fwd=2,
        ad_rev=3,
        dpf=80,
        rdf=76,
        adf=4,
        rdf_fwd=38,
        rdf_rev=38,
        adf_fwd=2,
        adf_rev=2,
        sb_pval=0.5,
        sb_or=1.0,
        fsb_pval=0.5,
        fsb_or=1.0,
        # Phase 2 fields
        any_alt=any_alt,
        partial_alt=partial_alt,
        # Phase 2b: N-base diagnostic
        n_count=n_count,
        # mFSD (all NaN/zero — not under test)
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
        # RNA fields (zeroed)
        sense_depth=0,
        antisense_depth=0,
        sense_strand_alt_count=0,
        rna_editing_site_overlap=False,
        splice_spanning_count=0,
    )


@pytest.fixture
def mock_variant():
    return Variant(
        chrom="chr1",
        pos=100,
        ref="A",
        alt="T",
        variant_type=VariantType.SNP,
        original_id="rs42",
    )


# ── Test 1: MAF column presence ──────────────────────────────────────────


def test_maf_has_diagnostic_columns():
    """MafWriter._gbcms_column_names() includes any_alt, partial_alt, n_count."""
    writer = MafWriter.__new__(MafWriter)
    writer.column_prefix = ""
    writer.mfsd = False
    writer.show_normalization = False
    writer.mode = "dna"
    writer.rescue_mnp = False

    cols = writer._gbcms_column_names()
    assert "any_alt" in cols, f"'any_alt' missing from MAF columns: {cols}"
    assert "partial_alt" in cols, f"'partial_alt' missing from MAF columns: {cols}"
    assert "n_count" in cols, f"'n_count' missing from MAF columns: {cols}"


# ── Test 2: MAF column values ────────────────────────────────────────────


def test_maf_diagnostic_values(tmp_path, mock_variant):
    """MafWriter writes correct any_alt, partial_alt, and n_count values."""
    counts = _mock_counts(any_alt=12, partial_alt=3, n_count=4)
    maf_path = tmp_path / "test.maf"
    writer = MafWriter(maf_path)
    writer.write(mock_variant, counts)
    writer.close()

    row = next(_read_maf_output(maf_path))

    assert row["any_alt"] == "12", f"Expected any_alt=12, got {row['any_alt']}"
    assert row["partial_alt"] == "3", f"Expected partial_alt=3, got {row['partial_alt']}"
    assert row["n_count"] == "4", f"Expected n_count=4, got {row['n_count']}"


# ── Test 3: MAF column prefix ────────────────────────────────────────────


def test_maf_diagnostic_columns_respect_prefix():
    """MafWriter with column_prefix='t_' produces t_any_alt, t_partial_alt, t_n_count."""
    writer = MafWriter.__new__(MafWriter)
    writer.column_prefix = "t_"
    writer.mfsd = False
    writer.show_normalization = False
    writer.mode = "dna"
    writer.rescue_mnp = False

    cols = writer._gbcms_column_names()
    assert "t_any_alt" in cols, "'t_any_alt' missing from prefixed columns"
    assert "t_partial_alt" in cols, "'t_partial_alt' missing from prefixed columns"
    assert "t_n_count" in cols, "'t_n_count' missing from prefixed columns"
    # Unprefixed must NOT be present
    assert "any_alt" not in cols, "'any_alt' (unprefixed) should not appear with prefix='t_'"
    assert "n_count" not in cols, "'n_count' (unprefixed) should not appear with prefix='t_'"


# ── Test 4: VCF INFO header ──────────────────────────────────────────────


def test_vcf_has_diagnostic_info_headers(tmp_path):
    """VcfWriter header includes AAD, PAD, and NAD INFO lines."""
    vcf_path = tmp_path / "test.vcf"
    writer = VcfWriter(vcf_path, sample_name="TUMOR")
    writer._write_header()
    writer.close()

    content = vcf_path.read_text()
    assert "##INFO=<ID=AAD," in content, "AAD INFO header missing from VCF"
    assert "##INFO=<ID=PAD," in content, "PAD INFO header missing from VCF"
    assert "##INFO=<ID=NAD," in content, "NAD INFO header missing from VCF"


# ── Test 5: VCF FORMAT header ────────────────────────────────────────────


def test_vcf_has_diagnostic_format_headers(tmp_path):
    """VcfWriter header includes AAD, PAD, and NAD FORMAT lines."""
    vcf_path = tmp_path / "test.vcf"
    writer = VcfWriter(vcf_path, sample_name="TUMOR")
    writer._write_header()
    writer.close()

    content = vcf_path.read_text()
    assert "##FORMAT=<ID=AAD," in content, "AAD FORMAT header missing from VCF"
    assert "##FORMAT=<ID=PAD," in content, "PAD FORMAT header missing from VCF"
    assert "##FORMAT=<ID=NAD," in content, "NAD FORMAT header missing from VCF"


# ── Test 6: VCF INFO data ────────────────────────────────────────────────


def test_vcf_info_contains_diagnostic_values(tmp_path, mock_variant):
    """VCF data row INFO field includes AAD=7;PAD=2;NAD=1."""
    counts = _mock_counts(any_alt=7, partial_alt=2, n_count=1)
    vcf_path = tmp_path / "test.vcf"
    writer = VcfWriter(vcf_path, sample_name="TUMOR")
    writer.write(mock_variant, counts)
    writer.close()

    data_lines = [line for line in vcf_path.read_text().splitlines() if not line.startswith("#")]
    assert len(data_lines) == 1
    info_field = data_lines[0].split("\t")[7]
    assert "AAD=7" in info_field, f"AAD=7 not in INFO: {info_field}"
    assert "PAD=2" in info_field, f"PAD=2 not in INFO: {info_field}"
    assert "NAD=1" in info_field, f"NAD=1 not in INFO: {info_field}"


# ── Test 7: VCF FORMAT data ──────────────────────────────────────────────


def test_vcf_format_contains_diagnostic_values(tmp_path, mock_variant):
    """VCF FORMAT column and sample data include AAD, PAD, and NAD values."""
    counts = _mock_counts(any_alt=7, partial_alt=2, n_count=3)
    vcf_path = tmp_path / "test.vcf"
    writer = VcfWriter(vcf_path, sample_name="TUMOR")
    writer.write(mock_variant, counts)
    writer.close()

    data_lines = [line for line in vcf_path.read_text().splitlines() if not line.startswith("#")]
    cols = data_lines[0].split("\t")
    format_str = cols[8]
    sample_data = cols[9]

    # FORMAT must include AAD, PAD, and NAD
    format_tags = format_str.split(":")
    assert "AAD" in format_tags, f"AAD not in FORMAT: {format_str}"
    assert "PAD" in format_tags, f"PAD not in FORMAT: {format_str}"
    assert "NAD" in format_tags, f"NAD not in FORMAT: {format_str}"

    # Parse and verify values
    data_dict = dict(zip(format_str.split(":"), sample_data.split(":"), strict=True))
    assert data_dict["AAD"] == "7", f"Expected AAD=7, got {data_dict.get('AAD')}"
    assert data_dict["PAD"] == "2", f"Expected PAD=2, got {data_dict.get('PAD')}"
    assert data_dict["NAD"] == "3", f"Expected NAD=3, got {data_dict.get('NAD')}"


# ── Test 8: _zero_counts invariant ───────────────────────────────────────


def test_zero_counts_has_diagnostic_fields():
    """_zero_counts() must include any_alt=0, partial_alt=0, n_count=0."""
    z = _zero_counts()
    assert hasattr(z, "any_alt"), "_zero_counts() missing any_alt field"
    assert hasattr(z, "partial_alt"), "_zero_counts() missing partial_alt field"
    assert hasattr(z, "n_count"), "_zero_counts() missing n_count field"
    assert z.any_alt == 0, f"Expected any_alt=0, got {z.any_alt}"
    assert z.partial_alt == 0, f"Expected partial_alt=0, got {z.partial_alt}"
    assert z.n_count == 0, f"Expected n_count=0, got {z.n_count}"
    # Invariant: any_alt = ad + partial_alt
    assert z.any_alt == z.ad + z.partial_alt, (
        f"Invariant violated in _zero_counts: any_alt={z.any_alt} "
        f"!= ad={z.ad} + partial_alt={z.partial_alt}"
    )


# ── Test 9: Column count delta ───────────────────────────────────────────


def test_column_count_delta_is_three():
    """Phase 2+2b adds 3 columns (any_alt, partial_alt, n_count) to MAF output.

    Compare against the known Phase 0/1 baseline of 21 columns (DNA, no mfsd).
    """
    writer = MafWriter.__new__(MafWriter)
    writer.column_prefix = ""
    writer.mfsd = False
    writer.show_normalization = False
    writer.mode = "dna"
    writer.rescue_mnp = False

    cols = writer._gbcms_column_names()
    # Phase 0/1 baseline was 21 (validation_status + 4 core + 4 frag + 4 bias + 8 strand)
    # Phase 2+2b added 3 (any_alt + partial_alt + n_count) = 24
    # v4.2.0 added 1 (gbcms_diagnostic, replaces validation_status with gbcms_status) = 25
    # gbcms_rescue is conditional (only with --rescue-mnp), so excluded here
    assert len(cols) == 25, (
        f"Expected 25 gbcms MAF columns (21 baseline + 3 diagnostic + 1 new column), "
        f"got {len(cols)}: {cols}"
    )


# ── Test 10: Structural invariants ───────────────────────────────────────


def test_mock_counts_invariants():
    """Verify structural invariants hold across mock count scenarios.

    Invariants:
      1. any_alt = ad + partial_alt
      2. any_alt >= ad
      3. dp >= rd + ad + partial_alt + n_count  (depth decomposition)
    """
    # SNP-like: no partial, no N
    snp = _mock_counts(ad=10, any_alt=10, partial_alt=0, n_count=0)
    assert snp.any_alt == snp.ad + snp.partial_alt
    assert snp.any_alt >= snp.ad

    # MNP-like: partial evidence
    # dp=100, rd=88, ad=5, partial_alt=2, n_count=1 → 88+5+2+1=96 ≤ 100 (4 other/third-allele)
    mnp = _mock_counts(ad=5, any_alt=7, partial_alt=2, n_count=1)
    mnp.rd = 88  # Override default rd=95 to satisfy decomposition
    assert mnp.any_alt == mnp.ad + mnp.partial_alt
    assert mnp.any_alt >= mnp.ad
    assert mnp.dp >= mnp.rd + mnp.ad + mnp.partial_alt + mnp.n_count, (
        f"Depth decomposition violated: DP({mnp.dp}) < "
        f"RD({mnp.rd}) + AD({mnp.ad}) + partial_alt({mnp.partial_alt}) + n_count({mnp.n_count})"
    )

    # N-heavy (duplex masking hotspot)
    n_heavy = _mock_counts(ad=1, any_alt=1, partial_alt=0, n_count=30)
    assert n_heavy.any_alt == n_heavy.ad + n_heavy.partial_alt
    assert n_heavy.any_alt >= n_heavy.ad

    # Edge: all zero
    z = _zero_counts()
    assert z.any_alt == z.ad + z.partial_alt
    assert z.any_alt >= z.ad
    assert z.dp >= z.ad + z.partial_alt + z.n_count


# ── Test 11: VCF GS/GD pipe separator ───────────────────────────────────


def test_vcf_gs_gd_use_pipe_separator(tmp_path, mock_variant):
    """VCF INFO GS/GD values convert ';' to '|' to avoid VCF delimiter conflict.

    VCF uses ';' as the INFO field delimiter. Multi-value GS/GD/GR fields
    must use '|' as their internal separator (design §6).
    """
    counts = _mock_counts(any_alt=0, partial_alt=0, n_count=0)
    vcf_path = tmp_path / "test.vcf"
    writer = VcfWriter(vcf_path, sample_name="TUMOR", rescue_mnp=True)
    # Pass multi-value status and diagnostic with semicolons
    writer.write(
        mock_variant,
        counts,
        gbcms_status="PASS;WARN_REF_CORRECTED",
        gbcms_diagnostic="ZERO_ALT;PARTIAL_DOMINANT",
        gbcms_rescue="method=decomposed;original_alt=0",
    )
    writer.close()

    data_lines = [line for line in vcf_path.read_text().splitlines() if not line.startswith("#")]
    assert len(data_lines) == 1
    info_field = data_lines[0].split("\t")[7]

    # Split on ';' (VCF INFO delimiter) and find GS, GD, GR
    # Use maxsplit=1 because GR values contain '=' (e.g., method=decomposed|original_alt=0)
    info_parts = {
        kv.split("=", 1)[0]: kv.split("=", 1)[1] for kv in info_field.split(";") if "=" in kv
    }

    # GS should use pipe, not semicolon
    assert (
        info_parts["GS"] == "PASS|WARN_REF_CORRECTED"
    ), f"GS should use pipe separator, got: {info_parts['GS']}"
    # GD should use pipe, not semicolon
    assert (
        info_parts["GD"] == "ZERO_ALT|PARTIAL_DOMINANT"
    ), f"GD should use pipe separator, got: {info_parts['GD']}"
    # GR should use pipe, not semicolon
    assert (
        info_parts["GR"] == "method=decomposed|original_alt=0"
    ), f"GR should use pipe separator, got: {info_parts['GR']}"
