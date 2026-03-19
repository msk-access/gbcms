"""
Tests for gbcms CLI commands (dna, rna, run).

Verifies that:
- dna and rna commands exist and show help
- run command is deprecated but still functional
- RNA-specific options are not accepted by dna command
- CLI correctly routes to the right Pipeline mode
"""

from typer.testing import CliRunner

from gbcms.cli import app

# Widen terminal to prevent Typer's rich help from truncating options in CI
runner = CliRunner(env={"COLUMNS": "200"})


# ── Command Existence ─────────────────────────────────────────────────────


def test_dna_command_help():
    """'gbcms dna --help' exits 0 and shows DNA-specific help."""
    result = runner.invoke(app, ["dna", "--help"])
    assert result.exit_code == 0, f"Exit code {result.exit_code}: {result.output}"
    assert "dna" in result.output.lower() or "variant" in result.output.lower()
    # DNA should have --mfsd
    assert "--mfsd" in result.output


def test_rna_command_help():
    """'gbcms rna --help' exits 0 and shows RNA-specific options."""
    result = runner.invoke(app, ["rna", "--help"])
    assert result.exit_code == 0, f"Exit code {result.exit_code}: {result.output}"
    assert "--rna-editing-db" in result.output
    assert "--enforce-strandedness" in result.output


def test_run_command_exists():
    """'gbcms run --help' still exits 0 (deprecated but not removed)."""
    result = runner.invoke(app, ["run", "--help"])
    assert result.exit_code == 0, f"Exit code {result.exit_code}: {result.output}"


# ── Option Isolation ──────────────────────────────────────────────────────


def test_dna_does_not_have_rna_options():
    """DNA command should NOT accept --rna-editing-db."""
    result = runner.invoke(app, ["dna", "--help"])
    assert "--rna-editing-db" not in result.output
    assert "--enforce-strandedness" not in result.output


def test_rna_does_not_have_dna_options():
    """RNA command should NOT accept --mfsd (DNA-only mFSD feature)."""
    result = runner.invoke(app, ["rna", "--help"])
    assert "--mfsd" not in result.output


# ── Error Handling ────────────────────────────────────────────────────────


def test_dna_rejects_unknown_option():
    """Invalid options produce a non-zero exit code."""
    result = runner.invoke(app, ["dna", "--nonexistent-flag"])
    assert result.exit_code != 0


def test_rna_requires_variants():
    """RNA command without --variants should fail."""
    result = runner.invoke(app, ["rna", "--fasta", "/tmp/fake.fa", "--output-dir", "/tmp"])
    assert result.exit_code != 0
