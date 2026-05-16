"""Tests for v5.3.0 provenance and formatting features.

Covers plan steps:
  4d — Format helper unit tests (_fmt, _fmt_vcf, _fmt_sci, _fmt_vcf_sci)
  4e — VCF provenance header validation
  4f — MAF provenance comment line validation
  4g — VCF INFO strand bias NaN → '.' sentinel
"""

import math
from types import SimpleNamespace

import pytest

from gbcms.io.output import MafWriter, VcfWriter, _fmt, _fmt_sci, _fmt_vcf, _fmt_vcf_sci
from gbcms.models.core import Variant, VariantType

from helpers import read_maf_output

_nan = float("nan")


# ── Shared mock counts ───────────────────────────────────────────────────


def _mock_counts(**overrides):
    """Minimal mock counts with all required fields. Override any field via kwargs."""
    defaults = dict(
        dp=10, rd=8, ad=2, dp_fwd=5, dp_rev=5, rd_fwd=4, rd_rev=4,
        ad_fwd=1, ad_rev=1, dpf=5, rdf=4, adf=1,
        rdf_fwd=2, rdf_rev=2, adf_fwd=0, adf_rev=1,
        sb_pval=1.0, sb_or=0.0, fsb_pval=1.0, fsb_or=0.0,
        any_alt=2, partial_alt=0, n_count=0,
        mfsd_ref_count=0, mfsd_alt_count=0, mfsd_nonref_count=0, mfsd_n_count=0,
        mfsd_ref_mean=_nan, mfsd_alt_mean=_nan, mfsd_nonref_mean=_nan, mfsd_n_mean=_nan,
        mfsd_alt_llr=_nan, mfsd_ref_llr=_nan,
        mfsd_delta_alt_ref=_nan, mfsd_ks_alt_ref=_nan, mfsd_pval_alt_ref=_nan,
        mfsd_delta_alt_nonref=_nan, mfsd_ks_alt_nonref=_nan, mfsd_pval_alt_nonref=_nan,
        mfsd_delta_ref_nonref=_nan, mfsd_ks_ref_nonref=_nan, mfsd_pval_ref_nonref=_nan,
        mfsd_delta_alt_n=_nan, mfsd_ks_alt_n=_nan, mfsd_pval_alt_n=_nan,
        mfsd_delta_ref_n=_nan, mfsd_ks_ref_n=_nan, mfsd_pval_ref_n=_nan,
        mfsd_delta_nonref_n=_nan, mfsd_ks_nonref_n=_nan, mfsd_pval_nonref_n=_nan,
    )
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


@pytest.fixture
def snp_variant():
    return Variant(chrom="chr1", pos=100, ref="A", alt="T", variant_type=VariantType.SNP)


# ── 4d: Format helper unit tests ─────────────────────────────────────────


class TestFormatHelpers:
    """Unit tests for the NaN/Inf-guarded format helpers."""

    # -- _fmt (MAF, fixed notation) --

    def test_fmt_nan(self):
        assert _fmt(float("nan")) == "NA"

    def test_fmt_inf(self):
        assert _fmt(float("inf")) == "NA"

    def test_fmt_neg_inf(self):
        assert _fmt(float("-inf")) == "NA"

    def test_fmt_zero(self):
        assert _fmt(0.0) == "0.0000"

    def test_fmt_normal(self):
        assert _fmt(1.2345) == "1.2345"

    # -- _fmt_vcf (VCF, fixed notation) --

    def test_fmt_vcf_nan(self):
        assert _fmt_vcf(float("nan")) == "."

    def test_fmt_vcf_inf(self):
        assert _fmt_vcf(float("inf")) == "."

    def test_fmt_vcf_neg_inf(self):
        assert _fmt_vcf(float("-inf")) == "."

    def test_fmt_vcf_zero(self):
        assert _fmt_vcf(0.0) == "0.0000"

    def test_fmt_vcf_normal(self):
        assert _fmt_vcf(1.5) == "1.5000"

    # -- _fmt_sci (MAF, scientific notation) --

    def test_fmt_sci_nan(self):
        assert _fmt_sci(float("nan")) == "NA"

    def test_fmt_sci_inf(self):
        assert _fmt_sci(float("inf")) == "NA"

    def test_fmt_sci_zero(self):
        assert _fmt_sci(0.0) == "0.0000e+00"

    def test_fmt_sci_normal(self):
        assert _fmt_sci(0.05) == "5.0000e-02"

    # -- _fmt_vcf_sci (VCF, scientific notation) --

    def test_fmt_vcf_sci_nan(self):
        assert _fmt_vcf_sci(float("nan")) == "."

    def test_fmt_vcf_sci_inf(self):
        assert _fmt_vcf_sci(float("inf")) == "."

    def test_fmt_vcf_sci_zero(self):
        assert _fmt_vcf_sci(0.0) == "0.0000e+00"

    def test_fmt_vcf_sci_normal(self):
        assert _fmt_vcf_sci(0.24) == "2.4000e-01"


# ── 4e: VCF provenance header tests ─────────────────────────────────────


def test_vcf_header_contains_source_version(tmp_path):
    """VCF header must contain ##source=gbcms with version."""
    vcf_path = tmp_path / "test.vcf"
    writer = VcfWriter(vcf_path, sample_name="TEST")
    writer._write_header()
    writer.close()

    content = vcf_path.read_text()
    assert "##source=gbcms v" in content, "Missing ##source with version"


def test_vcf_header_contains_command_when_provided(tmp_path):
    """VCF header must include ##gbcms_command when command_line is provided."""
    cmd = "gbcms dna --bam tumor.bam --fasta ref.fa"
    vcf_path = tmp_path / "test.vcf"
    writer = VcfWriter(vcf_path, sample_name="TEST", command_line=cmd)
    writer._write_header()
    writer.close()

    content = vcf_path.read_text()
    assert f"##gbcms_command={cmd}" in content, "Missing ##gbcms_command header"


def test_vcf_header_contains_reference_when_provided(tmp_path):
    """VCF header must include ##reference when reference_fasta is provided."""
    vcf_path = tmp_path / "test.vcf"
    writer = VcfWriter(vcf_path, sample_name="TEST", reference_fasta="/path/to/ref.fa")
    writer._write_header()
    writer.close()

    content = vcf_path.read_text()
    assert "##reference=file:///path/to/ref.fa" in content, "Missing ##reference header"


def test_vcf_header_contains_contigs(tmp_path):
    """VCF header must include ##contig lines for all provided contigs."""
    contigs = [("chr1", 248956422), ("chr2", 242193529)]
    vcf_path = tmp_path / "test.vcf"
    writer = VcfWriter(vcf_path, sample_name="TEST", contigs=contigs)
    writer._write_header()
    writer.close()

    content = vcf_path.read_text()
    assert "##contig=<ID=chr1,length=248956422>" in content, "Missing chr1 contig"
    assert "##contig=<ID=chr2,length=242193529>" in content, "Missing chr2 contig"


def test_vcf_header_contains_filter_pass(tmp_path):
    """VCF header must include ##FILTER=<ID=PASS,...> per VCF 4.2 spec."""
    vcf_path = tmp_path / "test.vcf"
    writer = VcfWriter(vcf_path, sample_name="TEST")
    writer._write_header()
    writer.close()

    content = vcf_path.read_text()
    assert '##FILTER=<ID=PASS,Description="All filters passed">' in content, (
        "Missing FILTER PASS header"
    )


# ── 4f: MAF provenance comment line tests ────────────────────────────────


def test_maf_has_version_comment_line(tmp_path):
    """MAF output must start with #gbcms vX.Y.Z provenance comment."""
    maf_path = tmp_path / "test.maf"
    writer = MafWriter(maf_path)
    writer.close()

    with open(maf_path) as f:
        first_line = f.readline()

    assert first_line.startswith("#gbcms v"), (
        f"Expected first line to start with '#gbcms v', got: {first_line!r}"
    )


def test_maf_has_command_comment_line(tmp_path):
    """MAF output must include #command line when command_line is provided."""
    cmd = "gbcms dna --bam tumor.bam --fasta ref.fa"
    maf_path = tmp_path / "test.maf"
    writer = MafWriter(maf_path, command_line=cmd)
    writer.close()

    with open(maf_path) as f:
        lines = f.readlines()

    command_lines = [l for l in lines if l.startswith("#command")]
    assert len(command_lines) == 1, f"Expected 1 #command line, found {len(command_lines)}"
    assert cmd in command_lines[0], f"Command not found in #command line: {command_lines[0]!r}"


def test_maf_comment_lines_skipped_by_reader(tmp_path, snp_variant):
    """read_maf_output helper must skip comment lines and parse TSV correctly."""
    maf_path = tmp_path / "test.maf"
    counts = _mock_counts()

    writer = MafWriter(maf_path, command_line="gbcms dna --test")
    writer.write(snp_variant, counts)
    writer.close()

    # The shared helper must parse the file correctly
    reader = read_maf_output(maf_path)
    row = next(reader)

    assert "Chromosome" in row, f"Expected 'Chromosome' column, got keys: {list(row.keys())}"
    assert row["alt_count"] == "2", f"Expected alt_count=2, got {row['alt_count']!r}"


# ── 4g: VCF INFO NaN → '.' sentinel test ────────────────────────────────


def test_vcf_info_nan_strand_bias_renders_as_dot(tmp_path, snp_variant):
    """VCF INFO fields SB_PVAL/SB_OR must render NaN as '.' (not literal 'nan')."""
    counts = _mock_counts(
        ad=0, any_alt=0, sb_pval=1.0, sb_or=_nan, fsb_pval=1.0, fsb_or=_nan,
    )

    vcf_path = tmp_path / "nan_test.vcf"
    writer = VcfWriter(vcf_path, sample_name="TUMOR")
    writer.write(snp_variant, counts)
    writer.close()

    data_lines = [l for l in vcf_path.read_text().splitlines() if not l.startswith("#")]
    assert len(data_lines) == 1

    info_field = data_lines[0].split("\t")[7]
    assert "nan" not in info_field.lower(), f"Literal 'nan' in VCF INFO: {info_field}"
    assert "inf" not in info_field.lower(), f"Literal 'inf' in VCF INFO: {info_field}"
    assert "SB_OR=." in info_field, f"Expected SB_OR=. for NaN, got: {info_field}"
    assert "FSB_OR=." in info_field, f"Expected FSB_OR=. for NaN, got: {info_field}"


# ── RNA mode provenance tests ───────────────────────────────────────────


def _rna_mock_counts(**overrides):
    """Mock counts with RNA-specific fields."""
    defaults = dict(
        dp=10, rd=8, ad=2, dp_fwd=5, dp_rev=5, rd_fwd=4, rd_rev=4,
        ad_fwd=1, ad_rev=1, dpf=5, rdf=4, adf=1,
        rdf_fwd=2, rdf_rev=2, adf_fwd=0, adf_rev=1,
        sb_pval=1.0, sb_or=0.0, fsb_pval=1.0, fsb_or=0.0,
        any_alt=2, partial_alt=0, n_count=0,
        # RNA-specific
        sense_depth=6, antisense_depth=4,
        sense_strand_alt_count=1,
        rna_editing_site_overlap=0,
        splice_spanning_count=0,
    )
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


def test_vcf_header_rna_mode_contains_provenance(tmp_path):
    """VCF provenance headers are emitted in RNA mode, not just DNA."""
    cmd = "gbcms rna --bam tumor.bam --fasta ref.fa"
    vcf_path = tmp_path / "rna.vcf"
    writer = VcfWriter(
        vcf_path,
        sample_name="RNA_TUMOR",
        mode="rna",
        command_line=cmd,
        reference_fasta="/path/to/ref.fa",
        contigs=[("chr1", 248956422)],
    )
    writer._write_header()
    writer.close()

    content = vcf_path.read_text()
    assert "##source=gbcms v" in content, "Missing ##source in RNA mode"
    assert f"##gbcms_command={cmd}" in content, "Missing ##gbcms_command in RNA mode"
    assert "##reference=file:///path/to/ref.fa" in content, "Missing ##reference in RNA mode"
    assert "##contig=<ID=chr1,length=248956422>" in content, "Missing ##contig in RNA mode"
    assert '##FILTER=<ID=PASS,Description="All filters passed">' in content, (
        "Missing ##FILTER in RNA mode"
    )


def test_maf_rna_mode_has_provenance_comments(tmp_path, snp_variant):
    """MAF provenance comment lines are emitted in RNA mode."""
    cmd = "gbcms rna --bam tumor.bam --fasta ref.fa"
    maf_path = tmp_path / "rna.maf"
    counts = _rna_mock_counts()

    writer = MafWriter(maf_path, mode="rna", command_line=cmd)
    writer.write(snp_variant, counts)
    writer.close()

    with open(maf_path) as f:
        lines = f.readlines()

    # First line is #gbcms version
    assert lines[0].startswith("#gbcms v"), f"Expected '#gbcms v', got: {lines[0]!r}"
    # Second line is #command
    assert lines[1].startswith("#command "), f"Expected '#command', got: {lines[1]!r}"
    assert cmd in lines[1], f"Command not found in #command line: {lines[1]!r}"

    # Data row must still be parsable via the shared helper
    reader = read_maf_output(maf_path)
    row = next(reader)
    assert row["alt_count"] == "2"
    # RNA-specific column must be present
    assert "rna_sense_depth" in row, f"Missing RNA column, got: {list(row.keys())}"
    assert row["rna_sense_depth"] == "6"

