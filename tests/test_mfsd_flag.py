"""Tests for the --mfsd and --mfsd-parquet CLI flags and OutputConfig model validation.

Covers:
1. mfsd=False (default) — mFSD columns absent from MAF and VCF headers
2. mfsd=True — mFSD columns present in MAF and VCF headers
3. mfsd_parquet=True requires mfsd=True (model validator)
4. CLI cross-validation: --mfsd-parquet without --mfsd exits with code 1
5. MafWriter column count with/without mfsd
6. VcfWriter header INFO line count with/without mfsd
7. _zero_counts() does not contain ref_sizes or alt_sizes
"""

from pathlib import Path

import pytest
from pydantic import ValidationError
from typer.testing import CliRunner

from gbcms.cli import app
from gbcms.io.output import MafWriter, VcfWriter
from gbcms.models.core import OutputConfig, OutputFormat
from gbcms.pipeline import _zero_counts

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_nan = float("nan")


class _MockCounts:
    """Minimal mock of BaseCounts for output writer tests."""

    def __init__(self, *, mfsd: bool = False):
        # Standard fields
        self.dp = 100
        self.rd = 95
        self.ad = 5
        self.dp_fwd = 50
        self.rd_fwd = 48
        self.ad_fwd = 2
        self.dp_rev = 50
        self.rd_rev = 47
        self.ad_rev = 3
        self.dpf = 80
        self.rdf = 76
        self.adf = 4
        self.rdf_fwd = 38
        self.rdf_rev = 38
        self.adf_fwd = 2
        self.adf_rev = 2
        self.sb_pval = 0.5
        self.sb_or = 1.0
        self.fsb_pval = 0.5
        self.fsb_or = 1.0
        self.used_decomposed = False
        # Decomposed ALT counting
        self.any_alt = 5  # same as ad for SNP (no partial)
        self.partial_alt = 0
        self.n_count = 0
        # mFSD fields (always present on BaseCounts even without --mfsd)
        self.mfsd_ref_count = 10
        self.mfsd_alt_count = 3
        self.mfsd_nonref_count = 1
        self.mfsd_n_count = 0
        self.mfsd_ref_mean = 170.0
        self.mfsd_alt_mean = 142.0
        self.mfsd_nonref_mean = _nan
        self.mfsd_n_mean = _nan
        self.mfsd_alt_llr = 2.3
        self.mfsd_ref_llr = -0.5
        self.mfsd_delta_alt_ref = -28.0
        self.mfsd_ks_alt_ref = 0.45
        self.mfsd_pval_alt_ref = 0.12
        self.mfsd_qval_alt_ref = 0.30  # HI-11: BH-FDR q-value (>= p)
        self.mfsd_delta_alt_nonref = _nan
        self.mfsd_ks_alt_nonref = _nan
        self.mfsd_pval_alt_nonref = _nan
        self.mfsd_delta_ref_nonref = _nan
        self.mfsd_ks_ref_nonref = _nan
        self.mfsd_pval_ref_nonref = _nan
        self.mfsd_delta_alt_n = _nan
        self.mfsd_ks_alt_n = _nan
        self.mfsd_pval_alt_n = _nan
        self.mfsd_delta_ref_n = _nan
        self.mfsd_ks_ref_n = _nan
        self.mfsd_pval_ref_n = _nan
        self.mfsd_delta_nonref_n = _nan
        self.mfsd_ks_nonref_n = _nan
        self.mfsd_pval_nonref_n = _nan
        # Sub-nucleosomal / mono-nucleosomal fractions
        self.mfsd_sub_nuc_ref_frac = 0.15
        self.mfsd_sub_nuc_alt_frac = 0.25
        self.mfsd_sub_nuc_enrichment = 1.67
        self.mfsd_mono_nuc_ref_frac = 0.60
        self.mfsd_mono_nuc_alt_frac = 0.45


class _MockVariant:
    """Minimal mock of a Variant for output writer tests."""

    def __init__(self):
        self.chrom = "chr1"
        self.pos = 999  # 0-based
        self.ref = "A"
        self.alt = "T"
        self.original_id = "."  # VCF ID column
        self.metadata = None  # triggers VCF→MAF path in MafWriter


# ---------------------------------------------------------------------------
# Test 1: Default mFSD off — columns absent from MAF header
# ---------------------------------------------------------------------------


def test_maf_writer_no_mfsd_columns_by_default(tmp_path: Path):
    """Without --mfsd, no mFSD columns should appear in the MAF output header."""
    writer = MafWriter(tmp_path / "out.maf")
    cols = writer._gbcms_column_names()
    mfsd_cols = [c for c in cols if c.startswith("mfsd_")]
    assert mfsd_cols == [], f"Expected no mFSD columns without --mfsd, got: {mfsd_cols}"
    writer.close()


# ---------------------------------------------------------------------------
# Test 2: mfsd=True — 41 mFSD columns present in MAF header
# ---------------------------------------------------------------------------


def test_maf_writer_mfsd_columns_when_enabled(tmp_path: Path):
    """With mfsd=True, exactly 41 mFSD columns should appear in the MAF header.

    Breakdown:
    - 4 raw counts (ref/alt/nonref/n)
    - 2 LLR (alt, ref)
    - 4 mean sizes (ref/alt/nonref/n)
    - 18 pairwise KS (6 pairs x 3: delta, D-stat, p-value)
    - 1 FDR q-value (alt_ref, HI-11)
    - 6 derived metrics (error_rate, n_rate, size_ratio, quality_score,
      alt_confidence, ks_valid)
    - 5 nucleosomal fractions (sub_nuc_ref/alt/enrichment, mono_nuc_ref/alt)
    - 1 CH gene flag (ch_flag)
    Total = 41
    """
    writer = MafWriter(tmp_path / "out.maf", mfsd=True)
    cols = writer._gbcms_column_names()
    mfsd_cols = [c for c in cols if c.startswith("mfsd_")]
    assert (
        len(mfsd_cols) == 41  # HI-11 added mfsd_qval_alt_ref
    ), f"Expected 41 mFSD columns with mfsd=True, got {len(mfsd_cols)}: {mfsd_cols}"
    writer.close()


# ---------------------------------------------------------------------------
# Test 3: VcfWriter — no MFSD INFO lines by default
# ---------------------------------------------------------------------------


def test_vcf_writer_no_mfsd_info_by_default(tmp_path: Path):
    """Without mfsd=True, VCF header should not contain MFSD ##INFO lines."""
    path = tmp_path / "out.vcf"
    writer = VcfWriter(path, sample_name="TUMOR")
    writer._write_header()
    writer.close()
    content = path.read_text()
    assert "MFSD_" not in content, "Expected no MFSD INFO lines in VCF header without --mfsd"


# ---------------------------------------------------------------------------
# Test 4: VcfWriter — 13 MFSD INFO lines when mfsd=True
# ---------------------------------------------------------------------------


def test_vcf_writer_mfsd_info_when_enabled(tmp_path: Path):
    """With mfsd=True, VCF header should contain exactly 13 ##INFO=<ID=MFSD_...> lines.

    8 core (delta/KS/pval/qval/alt-llr/ref-llr/alt-count/ref-count; HI-11 added
    MFSD_QVAL_ALT_REF) + 5 sub/mono-nucleosomal fractions (ME-1, VCF↔MAF parity).
    """
    path = tmp_path / "out.vcf"
    writer = VcfWriter(path, sample_name="TUMOR", mfsd=True)
    writer._write_header()
    writer.close()
    content = path.read_text()
    mfsd_lines = [line_ for line_ in content.splitlines() if "##INFO=<ID=MFSD_" in line_]
    assert len(mfsd_lines) == 13, f"Expected 13 MFSD INFO header lines, got {len(mfsd_lines)}"


def test_vcf_writer_mfsd_subnuc_values_in_data_row(tmp_path: Path):
    """ME-1: the sub/mono-nucleosomal fractions reach the VCF data-row INFO (not just the
    header), matching the MAF columns — VCF↔MAF parity."""
    from gbcms.io.output import _fmt_vcf

    path = tmp_path / "out.vcf"
    writer = VcfWriter(path, sample_name="TUMOR", mfsd=True)
    writer.write(_MockVariant(), _MockCounts(mfsd=True))
    writer.close()

    info = [ln for ln in path.read_text().splitlines() if not ln.startswith("#")][0].split("\t")[7]
    assert f"MFSD_SUB_NUC_REF_FRAC={_fmt_vcf(0.15)}" in info, info
    assert f"MFSD_SUB_NUC_ALT_FRAC={_fmt_vcf(0.25)}" in info, info
    assert f"MFSD_SUB_NUC_ENRICHMENT={_fmt_vcf(1.67)}" in info, info
    assert f"MFSD_MONO_NUC_REF_FRAC={_fmt_vcf(0.60)}" in info, info
    assert f"MFSD_MONO_NUC_ALT_FRAC={_fmt_vcf(0.45)}" in info, info


# ---------------------------------------------------------------------------
# Test 5: OutputConfig model validator — mfsd_parquet requires mfsd
# ---------------------------------------------------------------------------


def test_output_config_mfsd_parquet_requires_mfsd():
    """OutputConfig should raise ValidationError if mfsd_parquet=True without mfsd=True."""
    with pytest.raises(ValidationError, match="mfsd_parquet"):
        OutputConfig(
            directory=Path("/tmp"),
            format=OutputFormat.MAF,
            mfsd=False,
            mfsd_parquet=True,
        )


def test_output_config_mfsd_parquet_valid_when_mfsd_true(tmp_path: Path):
    """OutputConfig should accept mfsd_parquet=True when mfsd=True."""
    config = OutputConfig(
        directory=tmp_path,
        format=OutputFormat.MAF,
        mfsd=True,
        mfsd_parquet=True,
    )
    assert config.mfsd is True
    assert config.mfsd_parquet is True


# ---------------------------------------------------------------------------
# Test 6: CLI cross-validation — --mfsd-parquet without --mfsd exits 1
# ---------------------------------------------------------------------------


def test_cli_mfsd_parquet_without_mfsd_exits_1(tmp_path: Path):
    """The CLI should exit with code 1 and log an error if --mfsd-parquet is used without --mfsd.

    We supply all required args except --mfsd so that Typer argument parsing
    completes and our cross-validation code is reached.
    """
    # Create dummy files so Typer parameter validation passes
    bam = tmp_path / "fake.bam"
    vcf = tmp_path / "fake.vcf"
    fasta = tmp_path / "fake.fa"
    bam.touch()
    vcf.touch()
    fasta.touch()

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "dna",
            "--bam",
            str(bam),
            "--variants",
            str(vcf),
            "--fasta",
            str(fasta),
            "--output-dir",
            str(tmp_path),
            "--mfsd-parquet",  # --mfsd deliberately omitted
        ],
    )
    # Must exit non-zero
    assert result.exit_code != 0, "Expected non-zero exit when --mfsd-parquet used without --mfsd"
    # Must include a helpful error message
    output = result.output or ""
    assert "mfsd" in output.lower(), f"Expected error mentioning mfsd in output, got: {output!r}"


# ---------------------------------------------------------------------------
# Test 7: _zero_counts() does not expose ref_sizes or alt_sizes
# ---------------------------------------------------------------------------


def test_zero_counts_no_size_arrays():
    """_zero_counts() must not expose ref_sizes or alt_sizes — those are Rust-internal."""
    counts = _zero_counts()
    assert not hasattr(
        counts, "ref_sizes"
    ), "_zero_counts() should not have ref_sizes (internal Rust field)"
    assert not hasattr(
        counts, "alt_sizes"
    ), "_zero_counts() should not have alt_sizes (internal Rust field)"
