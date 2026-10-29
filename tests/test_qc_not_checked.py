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
