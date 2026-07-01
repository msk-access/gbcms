"""HI-1: sample failures (and empty/invalid runs) propagate to a non-zero exit code.

``Pipeline.run()`` catches per-sample errors and returns normally, so before this fix
the CLI exited ``0`` even when every sample failed — masking systematic failure as
success under Nextflow. Covered at two levels: the ``_exit_on_sample_failure()``
decision helper, and end-to-end through the ``dna`` command.
"""

import pytest
import typer
from typer.testing import CliRunner

from gbcms.cli import _exit_on_sample_failure, app

runner = CliRunner()


# ── Unit: the exit-code decision ──────────────────────────────────────────


def test_full_success_returns_without_exit():
    # processed > 0 and no failures → returns normally (the caller then exits 0).
    _exit_on_sample_failure({"samples_processed": 3, "failed_samples": []})


@pytest.mark.parametrize(
    "result",
    [
        {"samples_processed": 0, "failed_samples": [{"name": "s1", "error": "boom"}]},  # all failed
        {"samples_processed": 2, "failed_samples": [{"name": "s3", "error": "boom"}]},  # partial
        {"samples_processed": 0, "failed_samples": []},  # nothing processed (empty/invalid set)
        {},  # defensive: missing keys behave as "nothing processed"
    ],
)
def test_failure_and_empty_paths_exit_code_1(result):
    with pytest.raises(typer.Exit) as exc:
        _exit_on_sample_failure(result)
    assert exc.value.exit_code == 1


# ── End-to-end through the dna command ────────────────────────────────────


def test_dna_exits_nonzero_when_no_variants_processed(tmp_path, sample_bam, sample_fasta):
    """A run that processes zero samples (here: an empty variant set) exits non-zero —
    exercising the ``result → _exit_on_sample_failure`` wiring end-to-end."""
    empty_vcf = tmp_path / "empty.vcf"
    empty_vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    res = runner.invoke(
        app,
        ["dna", "-v", str(empty_vcf), "-b", sample_bam, "-f", sample_fasta, "-o", str(tmp_path)],
    )
    assert res.exit_code != 0, res.output


class _FakePipeline:
    """Stand-in for Pipeline that returns a fixed run() result — isolates the CLI
    exit-code wiring from the actual counting (which depends on fixture contigs)."""

    _result: dict = {}

    def __init__(self, config):
        pass

    def run(self):
        return self._result


def _patch_pipeline(monkeypatch, result):
    import gbcms.cli as cli

    _FakePipeline._result = result
    monkeypatch.setattr(cli, "Pipeline", _FakePipeline)


def test_dna_exits_zero_on_full_success(
    tmp_path, sample_bam, sample_fasta, sample_vcf, monkeypatch
):
    """A fully successful run exits 0 — the happy path is preserved."""
    _patch_pipeline(monkeypatch, {"samples_processed": 1, "failed_samples": []})
    res = runner.invoke(
        app,
        ["dna", "-v", sample_vcf, "-b", sample_bam, "-f", sample_fasta, "-o", str(tmp_path)],
    )
    assert res.exit_code == 0, res.output


def test_dna_exits_nonzero_when_a_sample_fails(
    tmp_path, sample_bam, sample_fasta, sample_vcf, monkeypatch
):
    """A run where a sample fails inside the engine exits non-zero (the core HI-1 case),
    exercising the failed_samples → exit-code path through the dna command."""
    _patch_pipeline(
        monkeypatch,
        {"samples_processed": 0, "failed_samples": [{"name": "s1", "error": "PyErr: panic"}]},
    )
    res = runner.invoke(
        app,
        ["dna", "-v", sample_vcf, "-b", sample_bam, "-f", sample_fasta, "-o", str(tmp_path)],
    )
    assert res.exit_code == 1, res.output
