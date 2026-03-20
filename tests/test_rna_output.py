"""
Tests for RNA-specific output columns in MafWriter and VcfWriter.

Verifies that:
- RNA mode adds 5 RNA-specific columns to MAF output
- DNA mode does NOT have RNA columns
- RNA mode VCF has SEN/ANT/ASEN/RED/SPL in headers
- DNA mode VCF does NOT have RNA headers
"""

from gbcms.io.output import MafWriter, VcfWriter

# ── RNA vs DNA column names in MAF ────────────────────────────────────────

RNA_MAF_COLUMNS = [
    "rna_sense_depth",
    "rna_antisense_depth",
    "rna_alt_sense_count",
    "rna_editing_site",
    "rna_splice_spanning",
]


def test_rna_maf_has_rna_columns():
    """MafWriter in RNA mode includes the 5 RNA-specific column names."""
    writer = MafWriter.__new__(MafWriter)
    writer.column_prefix = ""
    writer.mfsd = False
    writer.show_normalization = False
    writer.mode = "rna"

    columns = writer._gbcms_column_names()
    for col in RNA_MAF_COLUMNS:
        assert col in columns, f"RNA column '{col}' missing from MAF output"


def test_dna_maf_lacks_rna_columns():
    """MafWriter in DNA mode does NOT include RNA-specific columns."""
    writer = MafWriter.__new__(MafWriter)
    writer.column_prefix = ""
    writer.mfsd = False
    writer.show_normalization = False
    writer.mode = "dna"

    columns = writer._gbcms_column_names()
    for col in RNA_MAF_COLUMNS:
        assert col not in columns, f"RNA column '{col}' should not be in DNA MAF"


# ── RNA vs DNA headers in VCF ─────────────────────────────────────────────


RNA_VCF_TAGS = ["SEN", "ANT", "ASEN", "RED", "SPL"]


def test_rna_vcf_has_rna_headers(tmp_path):
    """VcfWriter in RNA mode includes SEN/ANT/ASEN/RED/SPL INFO headers."""
    vcf_path = tmp_path / "output.vcf"
    writer = VcfWriter(vcf_path, sample_name="TEST", mode="rna")
    writer._write_header()
    writer.close()

    header_text = vcf_path.read_text()
    for tag in RNA_VCF_TAGS:
        assert f"ID={tag}," in header_text, f"RNA VCF header missing INFO tag '{tag}'"


def test_dna_vcf_lacks_rna_headers(tmp_path):
    """VcfWriter in DNA mode does NOT include RNA-specific INFO headers."""
    vcf_path = tmp_path / "output.vcf"
    writer = VcfWriter(vcf_path, sample_name="TEST", mode="dna")
    writer._write_header()
    writer.close()

    header_text = vcf_path.read_text()
    for tag in RNA_VCF_TAGS:
        # Check INFO lines only — skip the FORMAT/column header line
        info_lines = [line for line in header_text.splitlines() if line.startswith("##INFO")]
        for line in info_lines:
            assert f"ID={tag}," not in line, f"DNA VCF should not have RNA INFO tag '{tag}'"


def test_rna_maf_column_count():
    """RNA mode adds exactly 5 columns vs DNA mode."""
    dna_writer = MafWriter.__new__(MafWriter)
    dna_writer.column_prefix = ""
    dna_writer.mfsd = False
    dna_writer.show_normalization = False
    dna_writer.mode = "dna"

    rna_writer = MafWriter.__new__(MafWriter)
    rna_writer.column_prefix = ""
    rna_writer.mfsd = False
    rna_writer.show_normalization = False
    rna_writer.mode = "rna"

    dna_cols = dna_writer._gbcms_column_names()
    rna_cols = rna_writer._gbcms_column_names()

    diff = len(rna_cols) - len(dna_cols)
    assert diff == 5, f"Expected 5 extra RNA columns, got {diff}"


# ── Write round-trip tests — assert actual values, not just column names ─────
# These tests are the regression guard for the pipeline.py mode= bug:
# if mode= is ever dropped from writer construction, these will fail.

import csv
import types

import pytest


def _make_rna_counts() -> types.SimpleNamespace:
    """Minimal mock counts object with all fields required by MafWriter/VcfWriter in RNA mode."""
    _nan = float("nan")
    return types.SimpleNamespace(
        # Core read counts
        dp=20, rd=15, ad=5,
        rd_fwd=10, rd_rev=5, ad_fwd=2, ad_rev=3,
        # Fragment counts
        dpf=10, rdf=7, adf=3,
        rdf_fwd=4, rdf_rev=3, adf_fwd=1, adf_rev=2,
        # Strand bias stats
        sb_pval=0.05, sb_or=1.5, fsb_pval=0.1, fsb_or=1.2,
        # mFSD (all NaN/zero — not tested here)
        mfsd_ref_count=0, mfsd_alt_count=0, mfsd_nonref_count=0, mfsd_n_count=0,
        mfsd_ref_mean=_nan, mfsd_alt_mean=_nan, mfsd_nonref_mean=_nan, mfsd_n_mean=_nan,
        mfsd_alt_llr=_nan, mfsd_ref_llr=_nan,
        mfsd_delta_alt_ref=_nan, mfsd_ks_alt_ref=_nan, mfsd_pval_alt_ref=_nan,
        mfsd_delta_alt_nonref=_nan, mfsd_ks_alt_nonref=_nan, mfsd_pval_alt_nonref=_nan,
        mfsd_delta_ref_nonref=_nan, mfsd_ks_ref_nonref=_nan, mfsd_pval_ref_nonref=_nan,
        mfsd_delta_alt_n=_nan, mfsd_ks_alt_n=_nan, mfsd_pval_alt_n=_nan,
        mfsd_delta_ref_n=_nan, mfsd_ks_ref_n=_nan, mfsd_pval_ref_n=_nan,
        mfsd_delta_nonref_n=_nan, mfsd_ks_nonref_n=_nan, mfsd_pval_nonref_n=_nan,
        # RNA-specific fields under test
        sense_depth=12,
        antisense_depth=3,
        sense_strand_alt_count=4,
        rna_editing_site_overlap=False,
        splice_spanning_count=2,
    )


@pytest.fixture
def rna_variant():
    """A simple SNP variant with no metadata (VCF-origin)."""
    from gbcms.models.core import Variant, VariantType
    return Variant(chrom="chr1", pos=100, ref="A", alt="T",
                   variant_type=VariantType.SNP, original_id="rs999")


def test_rna_vcf_data_row_has_rna_info_values(tmp_path, rna_variant):
    """VcfWriter(mode='rna') writes correct SEN/ANT/ASEN/SPL values in the data row INFO."""
    counts = _make_rna_counts()
    vcf_path = tmp_path / "rna.vcf"
    writer = VcfWriter(vcf_path, sample_name="TUMOR", mode="rna")
    writer.write(rna_variant, counts)
    writer.close()

    data_lines = [l for l in vcf_path.read_text().splitlines() if not l.startswith("#")]
    assert len(data_lines) == 1, "Expected exactly one data row"
    info_field = data_lines[0].split("\t")[7]

    assert "SEN=12" in info_field, f"SEN=12 not found in INFO: {info_field}"
    assert "ANT=3" in info_field, f"ANT=3 not found in INFO: {info_field}"
    assert "ASEN=4" in info_field, f"ASEN=4 not found in INFO: {info_field}"
    assert "SPL=2" in info_field, f"SPL=2 not found in INFO: {info_field}"


def test_rna_vcf_editing_flag_present_when_overlap(tmp_path, rna_variant):
    """RED flag is present in INFO when rna_editing_site_overlap=True."""
    counts = _make_rna_counts()
    counts.rna_editing_site_overlap = True
    vcf_path = tmp_path / "rna_editing.vcf"
    writer = VcfWriter(vcf_path, sample_name="TUMOR", mode="rna")
    writer.write(rna_variant, counts)
    writer.close()

    data_lines = [l for l in vcf_path.read_text().splitlines() if not l.startswith("#")]
    info_field = data_lines[0].split("\t")[7]
    assert "RED" in info_field.split(";"), \
        f"Expected 'RED' flag in INFO when rna_editing_site_overlap=True: {info_field}"


def test_rna_vcf_editing_flag_absent_when_no_overlap(tmp_path, rna_variant):
    """RED flag is absent from INFO when rna_editing_site_overlap=False."""
    counts = _make_rna_counts()
    counts.rna_editing_site_overlap = False
    vcf_path = tmp_path / "rna_no_editing.vcf"
    writer = VcfWriter(vcf_path, sample_name="TUMOR", mode="rna")
    writer.write(rna_variant, counts)
    writer.close()

    data_lines = [l for l in vcf_path.read_text().splitlines() if not l.startswith("#")]
    info_field = data_lines[0].split("\t")[7]
    assert "RED" not in info_field.split(";"), \
        f"'RED' flag should be absent when rna_editing_site_overlap=False: {info_field}"


def test_rna_maf_data_row_has_rna_column_values(tmp_path, rna_variant):
    """MafWriter(mode='rna') writes correct rna_* column values in the data row."""
    counts = _make_rna_counts()
    maf_path = tmp_path / "rna.maf"
    writer = MafWriter(maf_path, mode="rna")
    writer.write(rna_variant, counts, sample_name="TUMOR")
    writer.close()

    with open(maf_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        row = next(reader)

    assert row["rna_sense_depth"] == "12", f"Expected rna_sense_depth=12, got {row['rna_sense_depth']}"
    assert row["rna_antisense_depth"] == "3", f"Expected rna_antisense_depth=3"
    assert row["rna_alt_sense_count"] == "4", f"Expected rna_alt_sense_count=4"
    assert row["rna_editing_site"] == "False", f"Expected rna_editing_site=False"
    assert row["rna_splice_spanning"] == "2", f"Expected rna_splice_spanning=2"
