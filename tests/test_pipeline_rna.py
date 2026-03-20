"""
Pipeline-level integration test: RNA mode output columns.

Exercises Pipeline._write_output() directly (no Rust BAM engine call) to
verify that RNA-specific columns appear in both VCF and MAF output when
mode='rna'.  This is the regression guard for the pipeline.py bug where
mode= was not passed to VcfWriter/MafWriter, silently omitting RNA columns.

Does NOT require a real BAM or reference FASTA — only the output wiring
is under test here.  End-to-end counting accuracy is covered by test_accuracy.py.
"""

import csv
import types
from pathlib import Path

import pytest

from gbcms.io.output import MafWriter, VcfWriter
from gbcms.models.core import (
    GbcmsRnaConfig,
    OutputConfig,
    OutputFormat,
    Variant,
    VariantType,
)
from gbcms.pipeline import Pipeline, _zero_counts


# ── Shared fixtures ───────────────────────────────────────────────────────────


def _make_rna_counts() -> types.SimpleNamespace:
    """Full mock counts with realistic RNA field values."""
    zero = _zero_counts()
    # Override RNA-specific fields with non-zero values so their presence
    # in the output can be asserted.
    zero.dp = 30
    zero.rd = 20
    zero.ad = 10
    zero.rd_fwd = 12
    zero.rd_rev = 8
    zero.ad_fwd = 6
    zero.ad_rev = 4
    zero.rdf = 15
    zero.adf = 5
    zero.rdf_fwd = 8
    zero.rdf_rev = 7
    zero.adf_fwd = 3
    zero.adf_rev = 2
    zero.sense_depth = 18
    zero.antisense_depth = 5
    zero.sense_strand_alt_count = 9
    zero.rna_editing_site_overlap = False
    zero.splice_spanning_count = 3
    return zero


@pytest.fixture
def snp_variant() -> Variant:
    """A simple VCF-origin SNP variant (no metadata)."""
    return Variant(
        chrom="chr7", pos=149999, ref="C", alt="T",
        variant_type=VariantType.SNP, original_id=".",
    )


@pytest.fixture
def rna_pipeline(tmp_path) -> Pipeline:
    """A Pipeline instance configured for RNA mode.

    BAM files point to non-existent paths — _write_output() does not
    open BAM files, so this is safe for output-wiring tests.
    """
    dummy_vcf = tmp_path / "dummy.vcf"
    dummy_vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr7\t150000\t.\tC\tT\t.\t.\t.\n")
    dummy_bam = tmp_path / "sample.bam"
    dummy_bam.touch()  # must exist for Pydantic BAM validation
    dummy_fasta = tmp_path / "ref.fa"
    dummy_fasta.touch()
    dummy_fai = tmp_path / "ref.fa.fai"
    dummy_fai.touch()

    config = GbcmsRnaConfig(
        variant_file=dummy_vcf,
        bam_files={"rna_sample": dummy_bam},
        reference_fasta=dummy_fasta,
        output=OutputConfig(
            directory=tmp_path / "out",
            format=OutputFormat.VCF,
        ),
    )
    return Pipeline(config)


# ── VCF output tests ──────────────────────────────────────────────────────────


def test_rna_pipeline_vcf_has_rna_info_headers(rna_pipeline, snp_variant, tmp_path):
    """Pipeline._write_output() in RNA mode must include RNA ##INFO headers in VCF.

    This tests the full wiring from Pipeline._write_output() → VcfWriter(mode='rna')
    → _write_header().  If mode= is ever dropped from the VcfWriter constructor,
    the RNA ##INFO lines will be absent and this test will fail.
    """
    out_dir = tmp_path / "vcf_out"
    out_dir.mkdir()
    rna_pipeline.config.output.directory = out_dir
    rna_pipeline.config.output.format = OutputFormat.VCF

    counts = _make_rna_counts()
    rna_pipeline._write_output("rna_sample", [snp_variant], [counts], prepared=None)

    out_vcf = out_dir / "rna_sample.vcf"
    assert out_vcf.exists(), f"Expected output VCF at {out_vcf}"

    header_text = out_vcf.read_text()
    info_lines = [line for line in header_text.splitlines() if line.startswith("##INFO")]

    rna_tags = {"SEN", "ANT", "ASEN", "RED", "SPL"}
    for tag in rna_tags:
        assert any(f"ID={tag}," in line for line in info_lines), (
            f"RNA VCF INFO header 'ID={tag}' missing. "
            f"INFO lines found:\n" + "\n".join(info_lines)
        )


def test_rna_pipeline_vcf_data_row_has_rna_values(rna_pipeline, snp_variant, tmp_path):
    """Pipeline._write_output() in RNA mode writes correct RNA values into VCF data rows."""
    out_dir = tmp_path / "vcf_data_out"
    out_dir.mkdir()
    rna_pipeline.config.output.directory = out_dir
    rna_pipeline.config.output.format = OutputFormat.VCF

    counts = _make_rna_counts()
    rna_pipeline._write_output("rna_sample", [snp_variant], [counts], prepared=None)

    out_vcf = out_dir / "rna_sample.vcf"
    data_lines = [l for l in out_vcf.read_text().splitlines() if not l.startswith("#")]
    assert data_lines, "No data rows written to VCF"

    info_field = data_lines[0].split("\t")[7]
    assert "SEN=18" in info_field, f"SEN=18 not in INFO: {info_field}"
    assert "ANT=5" in info_field, f"ANT=5 not in INFO: {info_field}"
    assert "ASEN=9" in info_field, f"ASEN=9 not in INFO: {info_field}"
    assert "SPL=3" in info_field, f"SPL=3 not in INFO: {info_field}"
    # rna_editing_site_overlap=False → RED flag must be absent
    assert "RED" not in info_field.split(";"), \
        f"RED flag must be absent when rna_editing_site_overlap=False: {info_field}"


# ── MAF output tests ──────────────────────────────────────────────────────────


def test_rna_pipeline_maf_has_rna_column_headers(rna_pipeline, snp_variant, tmp_path):
    """Pipeline._write_output() in RNA mode must include rna_* columns in MAF header."""
    out_dir = tmp_path / "maf_out"
    out_dir.mkdir()
    rna_pipeline.config.output.directory = out_dir
    rna_pipeline.config.output.format = OutputFormat.MAF

    counts = _make_rna_counts()
    rna_pipeline._write_output("rna_sample", [snp_variant], [counts], prepared=None)

    out_maf = out_dir / "rna_sample.maf"
    assert out_maf.exists(), f"Expected output MAF at {out_maf}"

    with open(out_maf) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames or []

    rna_cols = {"rna_sense_depth", "rna_antisense_depth", "rna_alt_sense_count",
                "rna_editing_site", "rna_splice_spanning"}
    for col in rna_cols:
        assert col in fieldnames, (
            f"RNA MAF column '{col}' missing from header. "
            f"Columns found: {fieldnames}"
        )


def test_rna_pipeline_maf_data_row_has_rna_values(rna_pipeline, snp_variant, tmp_path):
    """Pipeline._write_output() in RNA mode writes correct rna_* values into MAF rows."""
    out_dir = tmp_path / "maf_data_out"
    out_dir.mkdir()
    rna_pipeline.config.output.directory = out_dir
    rna_pipeline.config.output.format = OutputFormat.MAF

    counts = _make_rna_counts()
    rna_pipeline._write_output("rna_sample", [snp_variant], [counts], prepared=None)

    out_maf = out_dir / "rna_sample.maf"
    with open(out_maf) as f:
        reader = csv.DictReader(f, delimiter="\t")
        row = next(reader)

    assert row["rna_sense_depth"] == "18", f"Expected 18, got {row['rna_sense_depth']!r}"
    assert row["rna_antisense_depth"] == "5", f"Expected 5, got {row['rna_antisense_depth']!r}"
    assert row["rna_alt_sense_count"] == "9", f"Expected 9, got {row['rna_alt_sense_count']!r}"
    assert row["rna_editing_site"] == "False", f"Expected False, got {row['rna_editing_site']!r}"
    assert row["rna_splice_spanning"] == "3", f"Expected 3, got {row['rna_splice_spanning']!r}"


def test_dna_pipeline_maf_lacks_rna_columns(tmp_path, snp_variant):
    """Pipeline._write_output() in DNA mode must NOT include rna_* columns in MAF."""
    from gbcms.models.core import GbcmsDnaConfig

    dummy_vcf = tmp_path / "dummy.vcf"
    dummy_vcf.write_text("#CHROM\tPOS\n")
    dummy_bam = tmp_path / "sample.bam"
    dummy_bam.touch()
    dummy_fasta = tmp_path / "ref.fa"
    dummy_fasta.touch()
    (tmp_path / "ref.fa.fai").touch()

    out_dir = tmp_path / "dna_maf_out"
    out_dir.mkdir()

    dna_config = GbcmsDnaConfig(
        variant_file=dummy_vcf,
        bam_files={"dna_sample": dummy_bam},
        reference_fasta=dummy_fasta,
        output=OutputConfig(directory=out_dir, format=OutputFormat.MAF),
    )
    pipeline = Pipeline(dna_config)
    pipeline._write_output("dna_sample", [snp_variant], [_zero_counts()], prepared=None)

    out_maf = out_dir / "dna_sample.maf"
    with open(out_maf) as f:
        fieldnames = csv.DictReader(f, delimiter="\t").fieldnames or []

    rna_cols = ["rna_sense_depth", "rna_antisense_depth", "rna_alt_sense_count",
                "rna_editing_site", "rna_splice_spanning"]
    for col in rna_cols:
        assert col not in fieldnames, (
            f"DNA MAF must NOT contain RNA column '{col}'. Found: {fieldnames}"
        )
