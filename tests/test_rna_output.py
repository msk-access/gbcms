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
