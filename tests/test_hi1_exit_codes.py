"""HI-1: an actual sample *failure* propagates to a non-zero exit code.

``Pipeline.run()`` catches per-sample errors and returns normally, so before this fix
the CLI exited ``0`` even when every sample failed — masking systematic failure as
success under Nextflow.

Crucially, an **empty variant set is NOT a failure**: a sample can legitimately have no
variants called, and per-sample workflows must not fail that task. So the exit code keys
off ``failed_samples`` (a sample that raised), never off "zero samples processed".
"""

import pytest
import typer
from typer.testing import CliRunner

from gbcms.cli import _exit_on_sample_failure, app

runner = CliRunner()


# ── Unit: the exit-code decision ──────────────────────────────────────────


@pytest.mark.parametrize(
    "result",
    [
        {"samples_processed": 0, "failed_samples": [{"name": "s1", "error": "boom"}]},  # all failed
        {"samples_processed": 2, "failed_samples": [{"name": "s3", "error": "boom"}]},  # partial
    ],
)
def test_actual_sample_failures_exit_code_1(result):
    with pytest.raises(typer.Exit) as exc:
        _exit_on_sample_failure(result)
    assert exc.value.exit_code == 1


@pytest.mark.parametrize(
    "result",
    [
        {"samples_processed": 3, "failed_samples": []},  # normal success
        {"samples_processed": 0, "failed_samples": []},  # empty/rejected variant set — legitimate
        {},  # defensive: missing keys → no failures recorded
    ],
)
def test_success_and_empty_variant_set_do_not_exit(result):
    # No sample raised → must NOT fail the run, even when zero samples were processed
    # (a sample with nothing called is legitimate under per-sample workflows).
    _exit_on_sample_failure(result)  # returns without raising


# ── End-to-end through the dna command ────────────────────────────────────


def test_dna_exits_zero_on_empty_variant_file(tmp_path, sample_bam, sample_fasta):
    """A sample with no variants called PASSES (exit 0) — an empty variant set is
    legitimate, not a failure. This is the per-sample-workflow guarantee."""
    empty_vcf = tmp_path / "empty.vcf"
    empty_vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    res = runner.invoke(
        app,
        ["dna", "-v", str(empty_vcf), "-b", sample_bam, "-f", sample_fasta, "-o", str(tmp_path)],
    )
    assert res.exit_code == 0, res.output


def test_dna_all_rejected_writes_reasons_to_output(tmp_path, sample_bam, sample_fasta):
    """When every variant is rejected in prep (here: a contig the reference lacks →
    FAIL_FETCH_FAILED), the run still exits 0 AND writes each variant with its reason in
    the `gbcms_status` column — so failure reasons live in the OUTPUT, not just the log."""
    vcf = tmp_path / "absent_contig.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr_absent>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr_absent\t100\t.\tA\tT\t.\t.\t.\n"
        "chr_absent\t200\t.\tC\tG\t.\t.\t.\n"
    )
    out = tmp_path / "out"
    out.mkdir()
    res = runner.invoke(
        app,
        [
            "dna",
            "-v",
            str(vcf),
            "-b",
            sample_bam,
            "-f",
            sample_fasta,
            "-o",
            str(out),
            "--format",
            "maf",
        ],
    )
    assert res.exit_code == 0, res.output

    mafs = list(out.glob("*.maf"))
    assert len(mafs) == 1, f"expected one MAF written, got {mafs}"
    rows = [ln for ln in mafs[0].read_text().splitlines() if ln and not ln.startswith("#")]
    header, *data = rows
    cols = header.split("\t")
    status_idx = cols.index("gbcms_status")
    reason_idx = cols.index("gbcms_status_reason")
    assert data, "rejected variants must still be written as rows"
    assert all(
        r.split("\t")[status_idx] == "FAIL" for r in data
    ), "every rejected variant must carry the FAIL verdict in gbcms_status"
    assert all(
        r.split("\t")[reason_idx] for r in data
    ), "every rejected variant must carry a non-empty reason in gbcms_status_reason"


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
