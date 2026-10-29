"""Greptile review of #352 and #356 (operator, 2026-10-28).

1. A sample whose output cannot be read is not a clean sample. The QC summary names it, the
   sample count says "k of n", and the closing line is a WARNING ("1 of 3 samples not
   checked: ...") even when no row of the other samples needs a look. The exit code is
   unchanged. `gbcms summarize` does the same for an input it refuses or cannot read: it
   summarizes the rest, names the input, and still exits 1 as it did.
2. `multiqc_table` works out which flags fired once, so its work grows with the cohort, not
   with its square.
"""

import json
import logging
import pathlib

import pytest
from test_input_group4_contract import _files, _maf, _ref_reads
from test_qc_summarize import cli, make_maf
from test_qc_summary import _rows, maf_row, write_maf
from typer.testing import CliRunner

from gbcms.cli import app
from gbcms.qc_export import (
    MULTIQC_TABLE_NAME,
    MULTIQC_WARNING_NAME,
    multiqc_documents,
    multiqc_table,
    summarize_inputs,
)
from gbcms.qc_summary import (
    QC_FLAGS,
    QcSummary,
    Skipped,
    closing_line,
    log_qc_summary,
    render_summary,
    summarize_output,
)

runner = CliRunner()


def make_vcf_output(path, rows=2):
    """A VCF shaped as the writer's: its header declares the GS, GSR and GD INFO fields."""
    header = (
        "##fileformat=VCFv4.2\n##source=gbcms v0.0.0 (test)\n"
        '##INFO=<ID=GS,Number=1,Type=String,Description="gbcms verdict: PASS or FAIL">\n'
        '##INFO=<ID=GSR,Number=1,Type=String,Description="reasons">\n'
        '##INFO=<ID=GD,Number=1,Type=String,Description="diagnostics">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    )
    body = [
        f"chr1\t{100 + i}\t.\tA\tT\t.\t.\tGS=PASS;GSR=.;GD=ZERO_ALT\tGT:DP\t0/0:5"
        for i in range(rows)
    ]
    path.write_text(header + "".join(line + "\n" for line in body))
    return path


def clean(tmp_path, name="A"):
    return summarize_output(write_maf(tmp_path / f"{name}.maf", [maf_row(0), maf_row(1)]), name)


# ── the block and the closing line ──────────────────────────────────────────────


def test_the_block_names_the_samples_not_checked_and_counts_k_of_n(tmp_path):
    text = render_summary([clean(tmp_path)], skipped=[Skipped("B", "OSError: disk went away")])
    lines = text.splitlines()
    assert lines[0] == "QC summary for 1 of 2 samples: 0 of 2 rows need a look (0.0%)"
    assert "  Not checked: 1 of 2 samples" in lines
    assert "    B: OSError: disk went away" in lines


def test_the_closing_line_warns_when_a_sample_was_not_checked_though_nothing_else_needs_a_look(
    tmp_path, caplog
):
    skipped = [Skipped("B", "OSError: disk went away")]
    level, message = closing_line([clean(tmp_path)], skipped)
    assert level == logging.WARNING
    assert message.startswith("1 of 2 samples not checked: B (OSError: disk went away)")
    assert "no row needs a look in the 1 checked (2 rows)" in message
    with caplog.at_level(logging.INFO, logger="gbcms.qc_summary"):
        log_qc_summary([clean(tmp_path)], skipped)
    assert caplog.records[-1].levelno == logging.WARNING
    assert "No row needs a look" not in caplog.records[-1].getMessage()


def test_the_closing_line_keeps_both_facts_when_rows_need_a_look_too(tmp_path):
    bad = summarize_output(make_maf(tmp_path / "C.maf"), "C")
    level, message = closing_line([clean(tmp_path), bad], [Skipped("B", "gone")])
    assert level == logging.WARNING
    assert message.startswith("3 of 8 rows need a look in 1 of 2 samples (1 with a problem")
    assert "1 of 3 samples not checked: B (gone)" in message


def test_nothing_checked_is_still_a_warning():
    level, message = closing_line([], [Skipped("A", "gone"), Skipped("B", "gone too")])
    assert level == logging.WARNING
    assert "2 of 2 samples not checked: A (gone), B (gone too)" in message
    assert "nothing could be checked" in message
    assert render_summary([], skipped=[Skipped("A", "gone")]).splitlines()[0] == (
        "QC summary for 0 of 1 samples: no output could be read"
    )


def test_with_nothing_read_the_log_is_one_warning_line_and_no_block(caplog):
    with caplog.at_level(logging.INFO, logger="gbcms.qc_summary"):
        log_qc_summary([], [Skipped("S", "OSError: gone")])
    assert [r.levelno for r in caplog.records] == [logging.WARNING]
    assert caplog.records[0].getMessage().startswith("1 of 1 sample not checked: S (OSError: gone)")
    assert "QC summary above" not in caplog.records[0].getMessage()  # there is none


def test_many_samples_not_checked_are_capped_in_the_line_and_the_block():
    skipped = [Skipped(f"S{i}", "gone") for i in range(14)]
    level, message = closing_line([], skipped)
    assert "14 of 14 samples not checked: S0 (gone), S1 (gone), S2 (gone) and 11 more" in message
    block = render_summary([], skipped=skipped)
    assert "    S9: gone" in block and "    S10: gone" not in block
    assert "and 4 more" in block


def test_a_run_with_every_sample_checked_is_unchanged(tmp_path):
    """No sample skipped: the block and the closing line are what they were."""
    s = clean(tmp_path)
    assert render_summary([s], skipped=[]) == render_summary([s])
    assert closing_line([s], []) == closing_line([s])
    assert closing_line([s])[0] == logging.INFO


# ── through the run ─────────────────────────────────────────────────────────────


def two_sample_run(tmp_path, monkeypatch, unreadable="B"):
    """`gbcms dna` on two BAMs; the summary of `unreadable` cannot read its output."""
    import gbcms.pipeline as pipeline

    real = pipeline.summarize_output

    def flaky(path, sample, *args, **kwargs):
        if sample == unreadable:
            raise OSError("disk went away")
        return real(path, sample, *args, **kwargs)

    monkeypatch.setattr(pipeline, "summarize_output", flaky)
    fa, bam = _files(tmp_path, _ref_reads())
    maf = _maf(tmp_path, _rows(tmp_path)[:1])  # one row with no ALT read: informational only
    res = runner.invoke(
        app,
        ["dna", "-v", str(maf), "-f", str(fa), "-o", str(tmp_path / "out"), "--format", "maf"]
        + ["-b", f"A:{bam}", "-b", f"{unreadable}:{bam}"],
        env={"COLUMNS": "3000"},
    )
    return res, [ln for ln in res.output.splitlines() if ln.strip()]


def test_a_run_whose_second_output_cannot_be_read_does_not_end_as_clean(tmp_path, monkeypatch):
    res, lines = two_sample_run(tmp_path, monkeypatch)
    assert res.exit_code == 0, res.output  # the exit code is unchanged
    text = " ".join(res.output.split())
    assert "QC summary for 1 of 2 samples" in text
    assert "Not checked: 1 of 2 samples" in text and "B: OSError: disk went away" in text
    assert "WARNING" in lines[-1] and "1 of 2 samples not checked: B" in " ".join(lines[-1].split())
    assert "No row needs a look" not in res.output
    # The per-sample warning is still there, ahead of the summary.
    assert "QC summary skipped for sample B (OSError: disk went away)" in text
    # Both outputs are on disk, unchanged by any of this.
    assert (tmp_path / "out" / "A.maf").exists() and (tmp_path / "out" / "B.maf").exists()


def test_a_sample_whose_counting_failed_is_not_checked_either(tmp_path):
    fa, bam = _files(tmp_path, _ref_reads())
    garbage = tmp_path / "bad.bam"
    garbage.write_bytes(b"this is not a bam")
    (tmp_path / "bad.bam.bai").write_bytes(b"x")
    maf = _maf(tmp_path, _rows(tmp_path)[:1])
    res = runner.invoke(
        app,
        ["dna", "-v", str(maf), "-f", str(fa), "-o", str(tmp_path / "out"), "--format", "maf"]
        + ["-b", f"A:{bam}", "-b", f"B:{garbage}"],
        env={"COLUMNS": "3000"},
    )
    text = " ".join(res.output.split())
    assert res.exit_code == 1, res.output  # a failed sample already exits 1
    assert "QC summary for 1 of 2 samples" in text
    assert "B: counting failed" in text
    assert "No row needs a look" not in res.output


# ── gbcms summarize ─────────────────────────────────────────────────────────────


def with_one_unreadable(tmp_path):
    good = make_maf(tmp_path / "good.maf")
    other = tmp_path / "plain.maf"
    other.write_text("Chromosome\tStart_Position\nchr1\t5\n")
    return good, other


def test_summarize_checks_the_rest_names_the_input_and_still_exits_one(tmp_path):
    good, other = with_one_unreadable(tmp_path)
    res = cli(good, other)
    assert res.exit_code == 1, res.output
    text = " ".join(res.stdout.split())
    assert "QC summary for 1 of 2 samples: 3 of 6 rows need a look" in text
    assert "Not checked: 1 of 2 samples" in text and "plain: plain.maf has no gbcms_status" in text
    assert res.stdout.rstrip().splitlines()[-1].startswith("3 of 6 rows need a look")
    assert "1 of 2 samples not checked: plain" in res.stdout.rstrip().splitlines()[-1]
    assert "has no gbcms_status column" in " ".join(res.output.split())


def test_summarize_with_nothing_readable_names_every_input_and_exits_one(tmp_path):
    other = tmp_path / "plain.maf"
    other.write_text("Chromosome\tStart_Position\nchr1\t5\n")
    odd = tmp_path / "s.tsv"
    odd.write_text("x")
    res = cli(other, odd)
    assert res.exit_code == 1
    text = " ".join(res.stdout.split())
    assert "QC summary for 0 of 2 samples: no output could be read" in text
    assert "2 of 2 samples not checked" in text and "expected a .maf or .vcf" in text


def test_summarize_json_is_written_for_the_readable_inputs(tmp_path):
    good, other = with_one_unreadable(tmp_path)
    res = cli("--format", "json", good, other)
    assert res.exit_code == 1
    assert (tmp_path / "good.qc_summary.json").exists()
    assert not (tmp_path / "plain.qc_summary.json").exists()
    assert "not checked" in " ".join(res.output.split())


def test_summarize_a_repeated_sample_name_checks_the_first_only(tmp_path):
    (tmp_path / "d1").mkdir()
    (tmp_path / "d2").mkdir()
    a = make_maf(tmp_path / "d1" / "S.maf")
    b = make_maf(tmp_path / "d2" / "S.maf")
    res = cli(a, b)
    assert res.exit_code == 1
    text = " ".join(res.stdout.split())
    assert "QC summary for 1 of 2 samples" in text and "share the sample name S" in text


def test_multiqc_says_a_sample_was_not_checked_even_when_nothing_else_needs_a_look(tmp_path):
    ok = clean(tmp_path)
    skipped = [Skipped("B", "OSError: disk went away")]
    docs = multiqc_documents([ok], skipped)
    assert set(docs) == {MULTIQC_TABLE_NAME, MULTIQC_WARNING_NAME}
    warning = docs[MULTIQC_WARNING_NAME]
    assert warning["section_name"].startswith("WARNING")
    assert "1 of 2 samples not checked" in warning["description"]
    assert warning["data"]["B"] == {"not_checked": 1}
    assert "not_checked" in warning["headers"]
    columns = {c for row in warning["data"].values() for c in row}
    assert columns <= set(warning["headers"])


def test_multiqc_files_through_the_cli_when_an_input_is_unreadable(tmp_path):
    good, other = with_one_unreadable(tmp_path)
    res = cli("--format", "multiqc", "-o", tmp_path / "mq", good, other)
    assert res.exit_code == 1
    warning = json.loads((tmp_path / "mq" / MULTIQC_WARNING_NAME).read_text())
    assert set(warning["data"]) == {"good", "plain"}


def test_summarize_inputs_returns_what_it_could_not_read(tmp_path):
    good, other = with_one_unreadable(tmp_path)
    summaries, skipped = summarize_inputs([good, other])
    assert [s.sample for s in summaries] == ["good"] and summaries[0].source == str(good)
    assert [(s.sample, "gbcms_status" in s.reason) for s in skipped] == [("plain", True)]


# ── multiqc_table: linear in the cohort ─────────────────────────────────────────


class CountingCounter(dict):
    """A flag-count map that counts every lookup and fails fast past a bound."""

    lookups = 0
    limit = 0

    def __missing__(self, key):
        return 0

    def __getitem__(self, key):
        CountingCounter.lookups += 1
        if CountingCounter.lookups > CountingCounter.limit:
            raise AssertionError("more flag lookups than a pass or two over the cohort")
        return dict.__getitem__(self, key) if key in self else 0


def test_the_multiqc_table_makes_a_pass_or_two_over_the_cohort_not_one_per_sample():
    """A typical cohort: three flags fire, the other thirty-nine never do. A flag that never
    fires is not found by a scan that stops at the first sample with it, so a scan per
    sample and flag costs the cohort squared."""
    n = 3000
    fired_names = ["ZERO_ALT", "PARTIAL_DOMINANT", "OBSERVED_ALLELE"]
    CountingCounter.lookups = 0
    CountingCounter.limit = 3 * n * len(QC_FLAGS)
    summaries = []
    for i in range(n):
        s = QcSummary(sample=f"S{i}", rows=100, need_look=i % 3)
        s.flag_rows = CountingCounter({fired_names[i % 3]: 1 + i % 5})  # type: ignore[assignment]
        summaries.append(s)
    doc = multiqc_table(summaries)
    assert len(doc["data"]) == n
    assert CountingCounter.lookups <= CountingCounter.limit
    # A header and a value in every row for each flag that fired, and for no other.
    assert [h for h in doc["headers"] if h in {f.name for f in QC_FLAGS}] == [
        f.name for f in QC_FLAGS if f.name in fired_names
    ]
    assert all(f in row for row in doc["data"].values() for f in fired_names)


def test_the_multiqc_table_is_what_it_was(tmp_path):
    """The headers keep the table's order and the rows carry every fired flag."""
    a = summarize_output(make_maf(tmp_path / "a.maf"), "a")
    b = clean(tmp_path, "b")
    doc = multiqc_table([a, b])
    flag_headers = [h for h in doc["headers"] if h in {f.name for f in QC_FLAGS}]
    assert flag_headers == [f.name for f in QC_FLAGS if a.flag_rows[f.name]]
    assert list(doc["data"]["b"]) == list(doc["data"]["a"])
    assert doc["data"]["b"]["PARTIAL_DOMINANT"] == 0
    with pytest.raises(KeyError):
        doc["headers"]["COMBINED_NA"]


def test_the_docs_say_what_a_sample_not_checked_looks_like():
    root = pathlib.Path(__file__).resolve().parents[1]
    page = (root / "docs/reference/qc-flags.md").read_text()
    assert 'under "Not checked"' in page and "samples not checked" in page
    assert "**not checked**" in (root / "docs/cli/summarize.md").read_text()


# ── Greptile, third round (#356) ────────────────────────────────────────────────


def plain_vcf(path, header="##fileformat=VCFv4.2\n"):
    """An ordinary VCF: no gbcms INFO fields, as any caller writes."""
    path.write_text(
        header
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        + "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT:DP\t0/1:10\n"
    )
    return path


@pytest.mark.parametrize("content", ["plain", "empty", "header-only-no-gbcms"])
def test_an_ordinary_or_empty_vcf_is_not_a_clean_gbcms_output(tmp_path, content):
    """Without gbcms's own fields (the GS INFO the writer declares) every row would read as
    a PASS row with no flags: a VCF that was never checked would read as a clean one."""
    path = tmp_path / "v.vcf"
    if content == "plain":
        plain_vcf(path)
    elif content == "empty":
        path.write_text("")
    else:
        path.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    with pytest.raises(ValueError, match="not a gbcms VCF"):
        summarize_output(path)
    summaries, skipped = summarize_inputs([path])
    assert summaries == [] and [s.sample for s in skipped] == ["v"]
    assert "not a gbcms VCF" in skipped[0].reason


def test_summarize_names_an_ordinary_vcf_as_not_checked_beside_a_real_output(tmp_path):
    good = make_vcf_output(tmp_path / "good.vcf")
    other = plain_vcf(tmp_path / "other.vcf")
    res = cli(good, other)
    assert res.exit_code == 1, res.output
    text = " ".join(res.stdout.split())
    assert "QC summary for 1 of 2 samples" in text
    assert "other: other.vcf has no GS INFO field: not a gbcms VCF" in text


def test_a_gbcms_vcf_with_no_records_is_a_checked_output_of_zero_rows(tmp_path):
    path = make_vcf_output(tmp_path / "none.vcf", rows=0)
    s = summarize_output(path)
    assert (s.rows, s.need_look) == (0, 0)


def test_the_run_writes_the_field_the_summary_asks_for(tmp_path):
    """The check is the writer's own header line, not a copy of it."""
    fa, bam = _files(tmp_path, _ref_reads())
    res = runner.invoke(
        app,
        ["dna", "-v", str(_maf(tmp_path, _rows(tmp_path)[:1])), "-f", str(fa), "-o"]
        + [str(tmp_path / "out"), "--format", "vcf", "-b", f"S:{bam}"],
        env={"COLUMNS": "3000"},
    )
    assert res.exit_code == 0, res.output
    assert summarize_output(tmp_path / "out" / "S.vcf").rows == 1


def test_a_skipped_input_with_a_checked_samples_name_does_not_replace_its_row(tmp_path):
    flagged = summarize_output(make_maf(tmp_path / "S.maf"))
    warning = multiqc_documents([flagged], [Skipped("S", "two outputs share the sample name S")])[
        MULTIQC_WARNING_NAME
    ]
    row = warning["data"]["S"]
    assert row["need_look"] == 3 and row["not_checked"] == 0  # the checked sample, intact
    (other,) = (k for k in warning["data"] if k != "S")
    assert warning["data"][other] == {"not_checked": 1} and "S" in other
    columns = {c for r in warning["data"].values() for c in r}
    assert columns <= set(warning["headers"])


def test_skipped_inputs_that_share_a_name_each_get_a_row():
    docs = multiqc_documents(
        [], [Skipped("S", "a"), Skipped("S", "b"), Skipped("S", "c"), Skipped("T", "d")]
    )
    data = docs[MULTIQC_WARNING_NAME]["data"]
    assert len(data) == 4 and all(r == {"not_checked": 1} for r in data.values())


def test_the_same_through_the_cli_when_two_inputs_share_a_name(tmp_path):
    (tmp_path / "d1").mkdir()
    (tmp_path / "d2").mkdir()
    a = make_maf(tmp_path / "d1" / "S.maf")
    b = make_maf(tmp_path / "d2" / "S.maf")
    res = cli("--format", "multiqc", "-o", tmp_path / "mq", a, b)
    assert res.exit_code == 1
    data = json.loads((tmp_path / "mq" / MULTIQC_WARNING_NAME).read_text())["data"]
    assert data["S"]["need_look"] == 3 and len(data) == 2


def test_a_table_left_by_an_earlier_run_is_removed_when_nothing_could_be_read(tmp_path):
    out = tmp_path / "mq"
    assert cli("--format", "multiqc", "-o", out, make_maf(tmp_path / "good.maf")).exit_code == 0
    assert (out / MULTIQC_TABLE_NAME).exists()
    other = plain_vcf(tmp_path / "other.vcf")
    res = cli("--format", "multiqc", "-o", out, other)
    assert res.exit_code == 1
    # MultiQC must not read the earlier run's counts as this run's.
    assert sorted(p.name for p in out.iterdir()) == [MULTIQC_WARNING_NAME]
    assert "not checked" in json.loads((out / MULTIQC_WARNING_NAME).read_text())["description"]
