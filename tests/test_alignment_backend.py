"""
Test alignment backend CLI flags and AlignmentConfig model.

Tests that:
1. Default backend is 'pairhmm' (normalized to 'hmm') with correct PairHMM parameter defaults
2. --alignment-backend hmm propagates to GbcmsConfig.alignment
3. Custom PairHMM parameters propagate correctly
4. Invalid backend name is rejected with clear error
"""

from unittest.mock import patch

import pytest
from typer.testing import CliRunner

from gbcms.cli import app
from gbcms.models.core import AlignmentConfig, GbcmsConfig

runner = CliRunner()


# ── Model-level tests ──


def test_alignment_config_defaults():
    """Default AlignmentConfig uses PairHMM (normalized to 'hmm') with documented defaults."""
    config = AlignmentConfig()
    assert config.backend == "pairhmm"  # default; normalized to 'hmm' only when explicitly passed
    assert config.hmm_llr_threshold == 2.3
    assert config.hmm_gap_open == 1e-4
    assert config.hmm_gap_extend == 0.1
    assert config.hmm_gap_open_repeat == 1e-2
    assert config.hmm_gap_extend_repeat == 0.5


def test_alignment_config_hmm():
    """AlignmentConfig accepts 'hmm' and 'pairhmm' backend values.

    'pairhmm' is the user-facing name (CLI --alignment-backend pairhmm),
    but the model normalizes it to 'hmm' for the Rust engine.
    """
    config_hmm = AlignmentConfig(backend="hmm")
    assert config_hmm.backend == "hmm"

    # 'pairhmm' → normalized to 'hmm' for Rust engine compatibility
    config_pairhmm = AlignmentConfig(backend="pairhmm")
    assert config_pairhmm.backend == "hmm"


def test_alignment_config_invalid_backend():
    """Invalid backend raises ValidationError with clear message."""
    with pytest.raises(Exception, match="Invalid alignment backend"):
        AlignmentConfig(backend="invalid")


def test_alignment_config_invalid_params():
    """Out-of-range gap probabilities are rejected."""
    with pytest.raises(ValueError):
        AlignmentConfig(hmm_gap_open=-0.1)  # negative

    with pytest.raises(ValueError):
        AlignmentConfig(hmm_gap_extend=1.5)  # > 1.0


# ── CLI-level tests ──


def _make_test_files(tmp_path):
    """Create dummy files for CLI invocation."""
    vcf = tmp_path / "test.vcf"
    vcf.touch()
    bam = tmp_path / "test.bam"
    bam.touch()
    fasta = tmp_path / "ref.fasta"
    fasta.touch()
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return vcf, bam, fasta, output_dir


def _base_args(vcf, bam, fasta, output_dir):
    """Build base CLI args for the 'run' command."""
    return [
        "dna",
        "-v",
        str(vcf),
        "-b",
        str(bam),
        "-f",
        str(fasta),
        "-o",
        str(output_dir),
    ]


@patch("gbcms.cli.Pipeline")
def test_cli_default_backend(mock_pipeline_cls, tmp_path):
    """Default invocation uses PairHMM backend."""
    vcf, bam, fasta, output_dir = _make_test_files(tmp_path)
    mock_pipeline_cls.return_value.run.return_value = {"samples_processed": 1, "failed_samples": []}

    result = runner.invoke(app, _base_args(vcf, bam, fasta, output_dir))
    assert result.exit_code == 0

    config = mock_pipeline_cls.call_args[0][0]
    assert isinstance(config, GbcmsConfig)
    assert config.alignment.backend == "hmm"  # CLI passes explicitly → validator normalizes
    assert config.alignment.hmm_llr_threshold == 2.3
    assert config.alignment.hmm_gap_open == 1e-4


@patch("gbcms.cli.Pipeline")
def test_cli_hmm_backend(mock_pipeline_cls, tmp_path):
    """--alignment-backend hmm propagates to config."""
    vcf, bam, fasta, output_dir = _make_test_files(tmp_path)
    mock_pipeline_cls.return_value.run.return_value = {"samples_processed": 1, "failed_samples": []}

    args = _base_args(vcf, bam, fasta, output_dir) + ["--alignment-backend", "hmm"]
    result = runner.invoke(app, args)
    assert result.exit_code == 0

    config = mock_pipeline_cls.call_args[0][0]
    assert config.alignment.backend == "hmm"


@patch("gbcms.cli.Pipeline")
def test_cli_custom_hmm_params(mock_pipeline_cls, tmp_path):
    """Custom PairHMM parameters propagate through CLI to config."""
    vcf, bam, fasta, output_dir = _make_test_files(tmp_path)
    mock_pipeline_cls.return_value.run.return_value = {"samples_processed": 1, "failed_samples": []}

    args = _base_args(vcf, bam, fasta, output_dir) + [
        "--alignment-backend",
        "hmm",
        "--llr-threshold",
        "3.0",
        "--gap-open-prob",
        "1e-3",
        "--gap-extend-prob",
        "0.2",
        "--repeat-gap-open-prob",
        "5e-2",
        "--repeat-gap-extend-prob",
        "0.6",
    ]
    result = runner.invoke(app, args)
    assert result.exit_code == 0

    config = mock_pipeline_cls.call_args[0][0]
    assert config.alignment.backend == "hmm"
    assert config.alignment.hmm_llr_threshold == 3.0
    assert config.alignment.hmm_gap_open == pytest.approx(1e-3)
    assert config.alignment.hmm_gap_extend == 0.2
    assert config.alignment.hmm_gap_open_repeat == pytest.approx(5e-2)
    assert config.alignment.hmm_gap_extend_repeat == 0.6


def test_cli_invalid_backend(tmp_path):
    """Invalid backend value exits with error."""
    vcf, bam, fasta, output_dir = _make_test_files(tmp_path)

    args = _base_args(vcf, bam, fasta, output_dir) + [
        "--alignment-backend",
        "invalid_backend",
    ]
    result = runner.invoke(app, args)
    # Should fail with validation error from AlignmentConfig
    assert result.exit_code != 0
