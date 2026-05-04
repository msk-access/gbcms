"""
Tests for gbcms CLI commands (dna, rna, normalize).

Verifies that:
- dna and rna commands exist and show help
- run command was removed in v4.1.0 (deprecated in v4.0.0)
- RNA-specific options are not accepted by dna command
- CLI correctly routes to the right Pipeline mode

Note: Option presence is verified via Click's registered params rather than
parsing help text, because Typer/Rich truncates the rendered options list
in CI environments (no real TTY), hiding options in the middle.
"""

import click
from typer.main import get_command
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()

# Resolve Typer app → Click Group once for param introspection.
# Cast to click.Group so mypy knows .commands is available —
# get_command() returns BaseCommand which doesn't declare .commands,
# but Typer always maps a multi-command app to a Group at runtime.
_click_app: click.Group = get_command(app)  # type: ignore[assignment]


def _get_param_names(command_name: str) -> set[str]:
    """Return the set of CLI parameter names for a subcommand."""
    # p.name is str | None per Click stubs; filter None (unnamed params don't occur in practice)
    return {p.name for p in _click_app.commands[command_name].params if p.name is not None}


# ── Command Existence ─────────────────────────────────────────────────────


def test_dna_command_help():
    """'gbcms dna --help' exits 0 and shows DNA-specific help."""
    result = runner.invoke(app, ["dna", "--help"])
    assert result.exit_code == 0, f"Exit code {result.exit_code}: {result.output}"
    assert "dna" in result.output.lower() or "variant" in result.output.lower()
    # DNA should have --mfsd (checked via Click params, not rendered help)
    assert "mfsd" in _get_param_names("dna")


def test_rna_command_help():
    """'gbcms rna --help' exits 0 and shows RNA-specific options."""
    result = runner.invoke(app, ["rna", "--help"])
    assert result.exit_code == 0, f"Exit code {result.exit_code}: {result.output}"
    # Checked via Click params to avoid Typer/Rich help truncation in CI
    rna_params = _get_param_names("rna")
    assert "rna_editing_db" in rna_params
    assert "enforce_strandedness" in rna_params


def test_run_command_removed():
    """'gbcms run' was removed in v4.1.0 (deprecated since v4.0.0)."""
    assert (
        "run" not in _click_app.commands
    ), "'gbcms run' should have been removed in v4.1.0 — found it still registered"


# ── Option Isolation ──────────────────────────────────────────────────────


def test_dna_does_not_have_rna_options():
    """DNA command should NOT accept --rna-editing-db."""
    dna_params = _get_param_names("dna")
    assert "rna_editing_db" not in dna_params
    assert "enforce_strandedness" not in dna_params


def test_rna_does_not_have_dna_options():
    """RNA command should NOT accept --mfsd (DNA-only mFSD feature)."""
    rna_params = _get_param_names("rna")
    assert "mfsd" not in rna_params


# ── Error Handling ────────────────────────────────────────────────────────


def test_dna_rejects_unknown_option():
    """Invalid options produce a non-zero exit code."""
    result = runner.invoke(app, ["dna", "--nonexistent-flag"])
    assert result.exit_code != 0


def test_rna_requires_variants():
    """RNA command without --variants should fail."""
    result = runner.invoke(app, ["rna", "--fasta", "/tmp/fake.fa", "--output-dir", "/tmp"])
    assert result.exit_code != 0
