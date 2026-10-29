"""O10 #227 follow-ups: ``gbcms summarize``, the merge NA flag, and ``--qc_summary``.

Contracts (operator, 2026-10-28):
- ``gbcms summarize`` reads finished MAF or VCF outputs through the code the run's own
  summary uses (one path makes every format) and never writes into an input.
- ``--format text`` is the run's block; ``json`` writes a ``<sample>.qc_summary.json`` per
  input that round-trips; ``tsv`` one cohort table; ``multiqc`` custom content with the
  fields MultiQC needs.
- A merged MAF is read under each input's prefix, and a row whose combined count, total or
  VAF is NA is a "check the call" row (COMBINED_NA), built from what ``gbcms merge``
  already writes: no column is added.
- Nextflow: one parameter, ``qc_summary`` (default false), runs the command after the samples.
- ``gbcms dna`` and ``gbcms rna`` gain no option.
"""

import csv
import json
import pathlib
import re
import shlex

import pytest
from test_input_group4_contract import _files, _maf, _ref_reads
from test_qc_summary import (
    _dna,
    _rows,
    maf_row,
    vcf_row,
    write_maf,
    write_vcf,
)
from typer.main import get_command
from typer.testing import CliRunner

from gbcms.cli import app
from gbcms.merge import merge_mafs
from gbcms.models.core import MergeConfig
from gbcms.qc_export import (
    JSON_SCHEMA,
    MULTIQC_TABLE_NAME,
    MULTIQC_WARNING_NAME,
    SummaryFormat,
    multiqc_documents,
    render_tsv,
    summarize_paths,
    summary_from_dict,
    summary_to_dict,
    write_summaries,
)
from gbcms.qc_summary import (
    QC_FLAGS,
    Kind,
    QcSummary,
    Tier,
    closing_line,
    render_summary,
    summarize_output,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
runner = CliRunner()


def cli(*args, **kw):
    return runner.invoke(app, ["summarize", *map(str, args)], env={"COLUMNS": "3000"}, **kw)


# ── Synthetic outputs, one per format ────────────────────────────────────────


def make_maf(path, extra=()):
    rows = [
        maf_row(0, status="FAIL", reason="REF_MISMATCH", diag="REF_AT_OFFSET(2)"),
        maf_row(1, diag="PARTIAL_DOMINANT;COEXISTING_ALLELE(chr1:50:A>AT:5/2)"),
        maf_row(2, gene="GENE2", diag="OBSERVED_ALLELE(chr1:60:C>T:9/0)"),
        maf_row(3, reason="MULTI_ALLELIC"),
        maf_row(4, diag="ZERO_ALT"),
        maf_row(5),
        *extra,
    ]
    return write_maf(path, rows)


def make_vcf(path):
    return write_vcf(
        path,
        [
            vcf_row(0, gs="FAIL", gsr="REF_MISMATCH", gd="REF_AT_OFFSET(2)"),
            vcf_row(1, gd="PARTIAL_DOMINANT|COEXISTING_ALLELE(chr1:50:A>AT:5/2)"),
            vcf_row(2, gd="OBSERVED_ALLELE(chr1:60:C>T:9/0)"),
            vcf_row(3, gsr="MULTI_ALLELIC"),
            vcf_row(4, gd="ZERO_ALT"),
            vcf_row(5),
        ],
    )


def outputs(tmp_path, names=("a", "b")):
    """A MAF and a VCF in separate folders (their inputs stay untouched)."""
    (tmp_path / "in").mkdir(exist_ok=True)
    return [
        make_maf(tmp_path / "in" / f"{names[0]}.maf"),
        make_vcf(tmp_path / "in" / f"{names[1]}.vcf"),
    ]


# ── text ──────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("kind", ["maf", "vcf"])
def test_text_is_the_runs_block_and_its_closing_line(tmp_path, kind):
    path = (make_maf if kind == "maf" else make_vcf)(tmp_path / f"s.{kind}")
    res = cli(path)
    assert res.exit_code == 0, res.output
    block = render_summary([summarize_output(path)])
    assert res.stdout == block + "\n" + closing_line([summarize_output(path)])[1] + "\n"
    assert res.stdout.startswith("QC summary for sample s: 3 of 6 rows need a look (50.0%)")


def test_text_for_several_samples_is_the_cohort_block(tmp_path):
    res = cli(*outputs(tmp_path))
    assert res.exit_code == 0, res.output
    assert "QC summary for 2 samples: 6 of 12 rows need a look (50.0%)" in res.stdout
    assert "    REF_MISMATCH: 2 rows in 2 of 2 samples" in res.stdout
    assert res.stdout.rstrip().splitlines()[-1].startswith("6 of 12 rows need a look in 2 of 2")


def test_text_equals_the_block_the_run_logged_for_the_same_file(tmp_path):
    fa, bam = _files(tmp_path, _ref_reads())
    run, _ = _dna(tmp_path, _maf(tmp_path, _rows(tmp_path)), bam, fa)
    lines = run.output.splitlines()
    first = next(i for i, ln in enumerate(lines) if "QC summary for sample" in ln)
    offset = lines[first].index("QC summary for sample")
    last = next(i for i in range(first, len(lines)) if "What each flag means" in lines[i])
    logged = "\n".join(ln[offset:].rstrip() for ln in lines[first : last + 1])
    res = cli(tmp_path / "out" / "S.maf")
    assert res.exit_code == 0, res.output
    assert res.stdout.startswith(logged + "\n")


def test_the_exit_code_is_zero_whether_or_not_rows_need_a_look(tmp_path):
    clean = write_maf(tmp_path / "clean.maf", [maf_row(0), maf_row(1)])
    res = cli(clean)
    assert res.exit_code == 0 and "No row needs a look: 2 rows checked" in res.stdout
    assert cli(make_maf(tmp_path / "flagged.maf")).exit_code == 0


# ── json ──────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("kind", ["maf", "vcf"])
def test_json_is_a_file_per_input_beside_it(tmp_path, kind):
    path = (make_maf if kind == "maf" else make_vcf)(tmp_path / f"s.{kind}")
    res = cli("--format", "json", path)
    assert res.exit_code == 0, res.output
    data = json.loads((tmp_path / "s.qc_summary.json").read_text())
    assert data["schema"] == JSON_SCHEMA and data["sample"] == "s" and data["format"] == kind
    assert (data["rows"], data["need_look"], data["need_look_share"]) == (6, 3, 0.5)
    assert data["rows_by_tier"] == {"fix_input": 1, "check_call": 2, "informational": 3}
    # Every flag the table knows, those that did not fire too: a zero is evidence.
    assert [f["name"] for f in data["flags"]] == [f.name for f in QC_FLAGS]
    fired = {f["name"]: f["rows"] for f in data["flags"] if f["rows"]}
    assert fired == {
        "REF_MISMATCH": 1,
        "REF_AT_OFFSET": 1,
        "PARTIAL_DOMINANT": 1,
        "COEXISTING_ALLELE": 1,
        "OBSERVED_ALLELE": 1,
        "MULTI_ALLELIC": 1,
        "ZERO_ALT": 1,
    }
    one = next(f for f in data["flags"] if f["name"] == "OBSERVED_ALLELE")
    assert one["tier"] == "check_call" and one["family"] == "diagnostic"
    assert one["docs"].endswith("/reference/qc-flags/#diagnostics") and one["action"]
    # Every row that needs a look, input problems first, with the flags in full.
    assert data["rows_to_look_at_complete"] is True
    assert [r["flags"] for r in data["rows_to_look_at"]] == [
        ["REF_MISMATCH", "REF_AT_OFFSET(2)"],
        ["PARTIAL_DOMINANT"],
        ["OBSERVED_ALLELE(chr1:60:C>T:9/0)"],
    ]
    assert data["rows_to_look_at"][2]["gene"] == ("GENE2" if kind == "maf" else "")
    # The MAF names its build on the first line; this VCF has no ##source line.
    assert data["written_by"] == ("gbcms v0.0.0 (test)" if kind == "maf" else "")
    assert data["summarized_by"].startswith("gbcms v")


def test_json_lists_every_row_not_only_the_first_few(tmp_path):
    rows = [maf_row(i, diag="PARTIAL_DOMINANT") for i in range(40)]
    path = write_maf(tmp_path / "many.maf", rows)
    assert cli("--format", "json", path).exit_code == 0
    data = json.loads((tmp_path / "many.qc_summary.json").read_text())
    assert len(data["rows_to_look_at"]) == data["need_look"] == 40
    # The text block still names five.
    assert "Rows to look at (first 5 of 40" in cli(path).stdout


def test_json_goes_to_output_dir_when_given(tmp_path):
    paths = outputs(tmp_path)
    res = cli("--format", "json", "--output-dir", tmp_path / "qc" / "deep", *paths)
    assert res.exit_code == 0, res.output
    assert sorted(p.name for p in (tmp_path / "qc" / "deep").iterdir()) == [
        "a.qc_summary.json",
        "b.qc_summary.json",
    ]
    assert not list((tmp_path / "in").glob("*.json"))


@pytest.mark.parametrize("kind", ["maf", "vcf"])
def test_json_round_trips(tmp_path, kind):
    path = (make_maf if kind == "maf" else make_vcf)(tmp_path / f"s.{kind}")
    original = summarize_output(path, keep_rows=True)
    again = summary_from_dict(json.loads(json.dumps(summary_to_dict(original))))
    for name in (
        "sample",
        "fmt",
        "rows",
        "need_look",
        "info_only",
        "unknown_only",
        "written_by",
        "merged",
        "flag_rows",
        "tier_rows",
        "unlisted",
        "examples",
        "look_rows",
    ):
        assert getattr(again, name) == getattr(original, name), name
    assert render_summary([again]) == render_summary([original])


def test_json_of_a_file_with_unlisted_flags_round_trips(tmp_path):
    path = write_maf(
        tmp_path / "u.maf",
        [maf_row(0, diag="NEW_FLAG(3);ZERO_ALT"), maf_row(1, status="FAIL", reason="NEW_REASON")],
    )
    original = summarize_output(path, keep_rows=True)
    again = summary_from_dict(json.loads(json.dumps(summary_to_dict(original))))
    assert dict(again.unlisted) == {"NEW_FLAG": 1, "NEW_REASON": 1}
    assert render_summary([again]) == render_summary([original])


def test_a_json_of_another_schema_is_refused():
    with pytest.raises(ValueError, match="not a gbcms QC summary"):
        summary_from_dict({"schema": "other/1"})


# ── tsv ───────────────────────────────────────────────────────────────────────


def test_tsv_is_a_cohort_table_on_standard_output(tmp_path):
    res = cli("--format", "tsv", *outputs(tmp_path))
    assert res.exit_code == 0, res.output
    table = list(csv.DictReader(res.stdout.splitlines(), delimiter="\t"))
    assert [r["sample"] for r in table] == ["a", "b"]
    assert list(table[0])[:10] == [
        "sample",
        "format",
        "rows",
        "need_look",
        "need_look_share",
        "fix_input",
        "check_call",
        "informational",
        "only_informational",
        "only_unlisted",
    ]
    # A column per flag in the table's order, then the flags it does not know.
    assert list(table[0])[10:-1] == [f.name for f in QC_FLAGS] and list(table[0])[-1] == "unlisted"
    row = table[0]
    assert (row["format"], row["rows"], row["need_look"], row["need_look_share"]) == (
        "maf",
        "6",
        "3",
        "0.5000",
    )
    assert (row["fix_input"], row["check_call"], row["informational"]) == ("1", "2", "3")
    assert row["REF_MISMATCH"] == "1" and row["COMBINED_NA"] == "0" and row["unlisted"] == ""
    assert table[1]["format"] == "vcf"


def test_tsv_goes_to_a_file_with_output_dir(tmp_path):
    res = cli("--format", "tsv", "-o", tmp_path / "t", *outputs(tmp_path))
    assert res.exit_code == 0, res.output
    assert (tmp_path / "t" / "gbcms_qc_summary.tsv").read_text() == render_tsv(
        summarize_paths(outputs(tmp_path))
    )


def test_tsv_names_flags_the_table_does_not_know(tmp_path):
    path = write_maf(tmp_path / "u.maf", [maf_row(0, diag="NEW_FLAG;NEW_FLAG2"), maf_row(1)])
    (row,) = csv.DictReader(render_tsv([summarize_output(path)]).splitlines(), delimiter="\t")
    assert row["unlisted"] == "NEW_FLAG:1;NEW_FLAG2:1" and row["only_unlisted"] == "1"


# ── multiqc ───────────────────────────────────────────────────────────────────


def test_multiqc_writes_a_table_and_a_warning_when_rows_need_a_look(tmp_path):
    res = cli("--format", "multiqc", "-o", tmp_path / "mq", *outputs(tmp_path))
    assert res.exit_code == 0, res.output
    files = sorted(p.name for p in (tmp_path / "mq").iterdir())
    assert files == sorted([MULTIQC_TABLE_NAME, MULTIQC_WARNING_NAME])
    assert all(f.endswith("_mqc.json") for f in files)  # the name MultiQC finds a file by
    table = json.loads((tmp_path / "mq" / MULTIQC_TABLE_NAME).read_text())
    warning = json.loads((tmp_path / "mq" / MULTIQC_WARNING_NAME).read_text())
    assert table["id"] != warning["id"]
    assert warning["section_name"].startswith("WARNING")
    assert set(warning["data"]) == {"a", "b"}


def test_the_multiqc_table_has_the_fields_multiqc_needs(tmp_path):
    summaries = summarize_paths(outputs(tmp_path))
    for name, doc in multiqc_documents(summaries).items():
        assert doc["id"] and doc["section_name"] and doc["plot_type"] == "table", name
        assert doc["pconfig"]["id"] and doc["pconfig"]["title"], name
        # data: {sample: {column: number}}, every column described by a header.
        assert doc["data"] and all(isinstance(r, dict) and r for r in doc["data"].values())
        columns = {c for row in doc["data"].values() for c in row}
        assert columns <= set(doc["headers"]), (name, columns - set(doc["headers"]))
        assert all(isinstance(v, int | float) for r in doc["data"].values() for v in r.values())
        for header in doc["headers"].values():
            assert header["title"] and header["description"], name
            for rules in header.get("cond_formatting_rules", {}).values():
                assert all(set(rule) <= {"eq", "gt"} for rule in rules)
    table = multiqc_documents(summaries)[MULTIQC_TABLE_NAME]
    assert set(table["data"]) == {"a", "b"}
    assert table["data"]["a"]["need_look"] == 3 and table["data"]["a"]["fix_input"] == 1
    assert table["data"]["a"]["need_look_pct"] == 50.0
    # A column per flag that fired anywhere, hidden until switched on; none that never did.
    assert table["headers"]["PARTIAL_DOMINANT"]["hidden"] is True
    assert "COMBINED_NA" not in table["headers"] and "mfsd_ch_flag" not in table["headers"]
    assert table["section_href"].endswith("/reference/qc-flags/")
    assert "qc-flags" in table["helptext"]


def test_no_warning_section_when_no_row_needs_a_look(tmp_path):
    clean = write_maf(tmp_path / "clean.maf", [maf_row(0), maf_row(1, diag="ZERO_ALT")])
    docs = multiqc_documents(summarize_paths([clean]))
    assert list(docs) == [MULTIQC_TABLE_NAME]


def test_a_warning_left_by_an_earlier_run_is_removed(tmp_path):
    out = tmp_path / "mq"
    assert cli("--format", "multiqc", "-o", out, make_maf(tmp_path / "bad.maf")).exit_code == 0
    assert (out / MULTIQC_WARNING_NAME).exists()
    clean = write_maf(tmp_path / "ok.maf", [maf_row(0)])
    assert cli("--format", "multiqc", "-o", out, clean).exit_code == 0
    assert sorted(p.name for p in out.iterdir()) == [MULTIQC_TABLE_NAME]


def test_multiqc_defaults_to_the_current_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    path = make_maf(tmp_path / "in.maf")
    assert cli("--format", "multiqc", path).exit_code == 0
    assert (tmp_path / MULTIQC_TABLE_NAME).exists()


# ── the command never writes into an input ────────────────────────────────────


@pytest.mark.parametrize("fmt", ["text", "json", "tsv", "multiqc"])
def test_no_format_writes_into_an_input(tmp_path, fmt):
    paths = outputs(tmp_path)
    before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in paths}
    before_names = sorted(p.name for p in (tmp_path / "in").iterdir())
    res = cli("--format", fmt, "-o", tmp_path / "out", *paths)
    assert res.exit_code == 0, res.output
    assert {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in paths} == before
    assert sorted(p.name for p in (tmp_path / "in").iterdir()) == before_names


def test_json_beside_the_input_adds_only_the_json(tmp_path):
    paths = outputs(tmp_path)
    cli("--format", "json", *paths)
    assert sorted(p.name for p in (tmp_path / "in").iterdir()) == [
        "a.maf",
        "a.qc_summary.json",
        "b.qc_summary.json",
        "b.vcf",
    ]


def test_an_output_never_replaces_an_input(tmp_path):
    victim = tmp_path / "x.qc_summary.json"
    victim.write_text("an input")
    with pytest.raises(ValueError, match="is an input"):
        write_summaries([QcSummary(sample="x")], [victim], SummaryFormat.JSON)
    assert victim.read_text() == "an input"


# ── refusals ──────────────────────────────────────────────────────────────────


def test_two_outputs_with_one_sample_name_are_refused(tmp_path):
    (tmp_path / "d1").mkdir()
    (tmp_path / "d2").mkdir()
    a = make_maf(tmp_path / "d1" / "S.maf")
    b = make_maf(tmp_path / "d2" / "S.maf")
    res = cli(a, b)
    assert res.exit_code == 1 and "share the sample name S" in " ".join(res.output.split())


def test_a_file_that_is_not_a_gbcms_output_exits_one(tmp_path):
    other = tmp_path / "plain.maf"
    other.write_text("Chromosome\tStart_Position\nchr1\t5\n")
    res = cli(other)
    assert res.exit_code == 1 and "gbcms_status" in " ".join(res.output.split())


@pytest.mark.parametrize("name", ["s.maf.gz", "s.tsv", "s"])
def test_other_file_types_are_refused(tmp_path, name):
    path = tmp_path / name
    path.write_text("x")
    res = cli(path)
    assert res.exit_code == 1 and "expected a .maf or .vcf" in " ".join(res.output.split())


def test_a_missing_file_is_a_usage_error(tmp_path):
    assert cli(tmp_path / "nope.maf").exit_code == 2


def test_dna_and_rna_gain_no_option():
    """The summary is a command of its own, not a flag of the runs."""
    group = get_command(app)
    assert "summarize" in group.commands  # type: ignore[attr-defined]
    for mode in ("dna", "rna"):
        names = {p.name for p in group.commands[mode].params}  # type: ignore[attr-defined]
        assert not {n for n in names if "summar" in n}, (mode, names)


# ── a merged MAF ─────────────────────────────────────────────────────────────

_MERGE_COLS = [
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "gbcms_status",
    "gbcms_status_reason",
    "gbcms_diagnostic",
    "ref_count",
    "alt_count",
    "ref_count_forward",
    "ref_count_reverse",
    "alt_count_forward",
    "alt_count_reverse",
    "ref_count_fragment",
    "alt_count_fragment",
]


def _flavour_row(i, diag="", status="PASS", reason="", ref="10", alt="2", zero=False):
    strand = ("0", "0", "0", "0") if zero else ("5", "5", "1", "1")
    return dict(
        zip(
            _MERGE_COLS,
            [
                "G",
                "1",
                str(100 + i),
                str(100 + i),
                "A",
                "T",
                status,
                reason,
                diag,
                ref,
                alt,
                *strand,
                ref,
                alt,
            ],
            strict=True,
        )
    )


def _write_flavour(path, rows):
    lines = ["#gbcms v0 (test)", "\t".join(_MERGE_COLS)]
    lines += ["\t".join(r[c] for c in _MERGE_COLS) for r in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def merged_maf(tmp_path, legacy=False):
    zero = {"ref": "0", "alt": "0", "zero": True}
    duplex = _write_flavour(
        tmp_path / "duplex.maf",
        [
            _flavour_row(0),  # clean in both
            _flavour_row(1, diag="ZERO_ALT", **zero),  # all zero: a benign NA odds ratio
            _flavour_row(2, ref="NA"),  # a missing count: the combined cells are NA
            _flavour_row(3, diag="PARTIAL_DOMINANT"),  # a flag in one flavour only
            _flavour_row(4, status="FAIL", reason="REF_MISMATCH", **zero),
        ],
    )
    simplex = _write_flavour(
        tmp_path / "simplex.maf",
        [
            _flavour_row(0),
            _flavour_row(1, **zero),
            _flavour_row(2),
            _flavour_row(3, diag="OBSERVED_ALLELE(1:103:A>T:5/0)"),
            _flavour_row(4, status="FAIL", reason="REF_MISMATCH", **zero),
        ],
    )
    out = tmp_path / ("merged_legacy.maf" if legacy else "merged.maf")
    merge_mafs(
        MergeConfig(
            inputs={"duplex": duplex, "simplex": simplex},
            output=out,
            add_combined=True,
            legacy_naming=legacy,
        )
    )
    return out


@pytest.mark.parametrize("legacy", [False, True])
def test_a_merged_maf_is_read_under_each_inputs_prefix(tmp_path, legacy):
    path = merged_maf(tmp_path, legacy)
    s = summarize_output(path, keep_rows=True)
    assert s.merged is True and s.rows == 5
    # The union of the inputs' flags; one FAIL in either is a FAIL.
    assert dict(s.flag_rows) == {
        "ZERO_ALT": 1,
        "COMBINED_NA": 1,
        "PARTIAL_DOMINANT": 1,
        "OBSERVED_ALLELE": 1,
        "REF_MISMATCH": 1,
    }
    assert s.need_look == 3 and s.info_only == 1  # a row with two flags is one row
    by_start = {v.split()[0].split(":")[1]: flags for v, _, flags in s.look_rows}
    assert by_start["104"] == ["REF_MISMATCH"]
    assert by_start["102"] == ["COMBINED_NA"]
    assert by_start["103"] == ["PARTIAL_DOMINANT", "OBSERVED_ALLELE(1:103:A>T:5/0)"]
    assert "101" not in by_start and "100" not in by_start


def test_a_strand_bias_na_alone_is_not_the_flag(tmp_path):
    """The all-zero row has an NA strand-bias odds ratio with every count present."""
    path = merged_maf(tmp_path)
    lines = [ln for ln in path.read_text().splitlines() if not ln.startswith("#")]
    header = lines[0].split("\t")
    cells = next(ln for ln in lines if "\t101\t" in ln).split("\t")
    zero_row = dict(zip(header, cells, strict=True))
    assert zero_row["simplex_duplex_strand_bias_odds_ratio"] == "NA"
    assert zero_row["simplex_duplex_alt_count"] == "0"
    assert summarize_output(path).flag_rows["COMBINED_NA"] == 1  # the missing-count row only


def test_the_merge_flag_is_check_the_call_and_names_each_cause():
    flag = next(f for f in QC_FLAGS if f.name == "COMBINED_NA")
    assert flag.tier is Tier.CHECK_CALL and flag.kind is Kind.MERGE
    assert "missing" in flag.action and "rescue" in flag.action


def test_the_merge_writes_no_new_column_for_it(tmp_path):
    """The flag is built from the cells the merge already writes."""
    header = merged_maf(tmp_path).read_text().splitlines()[1].split("\t")
    assert not [c for c in header if "combined_na" in c.lower() or "na_flag" in c.lower()]


def test_a_merged_maf_through_the_cli_names_the_prefixed_fields(tmp_path):
    res = cli(merged_maf(tmp_path))
    assert res.exit_code == 0, res.output
    text = " ".join(res.stdout.split())
    assert "COMBINED_NA: 1 row" in text
    assert "a merged MAF has each under its input's prefix" in text
    assert "simplex_duplex_*" in text


def test_mixed_mnp_rescue_makes_the_combined_cells_na_and_the_flag(tmp_path):
    from helpers import read_maf_output

    cols = [*_MERGE_COLS, "gbcms_rescue"]
    rescued = (
        "method=decomposed;outcome=rescued;original_ref=9;original_alt=0;original_partial=4;"
        "original_confirmed=0;adopted=1:100(A>T)"
    )
    kept = rescued.replace("rescued", "no_improvement").split(";adopted")[0]
    paths = {}
    for flavour, audit in (("duplex", rescued), ("simplex", kept)):
        row = _flavour_row(0) | {"gbcms_rescue": audit}
        lines = ["#gbcms v0 (test)", "\t".join(cols), "\t".join(row[c] for c in cols)]
        paths[flavour] = tmp_path / f"{flavour}.maf"
        paths[flavour].write_text("\n".join(lines) + "\n")
    out = tmp_path / "mixed.maf"
    merge_mafs(MergeConfig(inputs=paths, output=out, add_combined=True))
    assert next(iter(read_maf_output(out)))["simplex_duplex_alt_count"] == "NA"
    s = summarize_output(out)
    assert s.flag_rows["COMBINED_NA"] == 1 and s.flag_rows["rescued"] == 1
    assert s.flag_rows["no_improvement"] == 1 and s.need_look == 1


# ── Nextflow: one parameter, a small module ──────────────────────────────────

NEXTFLOW = ROOT / "nextflow"


def test_qc_summary_is_one_parameter_off_by_default():
    config = (NEXTFLOW / "nextflow.config").read_text()
    assert re.search(r"^\s+qc_summary\s*=\s*false\b", config, re.M)
    assert not re.search(r"^\s+qc_\w+\s*=", config.replace("qc_summary", ""), re.M)


def test_the_module_runs_summarize_in_both_formats_over_the_collected_outputs():
    script = (NEXTFLOW / "modules/local/gbcms/qc_summary/main.nf").read_text()
    commands = re.findall(r"^\s+(gbcms summarize .*)$", script, re.M)
    assert [re.search(r"--format (\w+)", c).group(1) for c in commands] == ["multiqc", "json"]
    group = get_command(app).commands["summarize"]  # type: ignore[attr-defined]
    for command in commands:
        # The command parses against the CLI, with this file standing for the staged outputs.
        tokens = shlex.split(command.replace("${outputs}", str(pathlib.Path(__file__))))
        group.make_context("summarize", tokens[2:])
    assert "path outputs" in script and "emit: multiqc" in script and "emit: json" in script


def test_the_workflow_runs_it_only_with_the_parameter_after_the_samples():
    main = (NEXTFLOW / "main.nf").read_text()
    assert "include { QC_SUMMARY }" in main
    block = main[main.index("STEP 4") :]
    assert "if (asBool(params.qc_summary))" in block and "QC_SUMMARY(" in block
    assert ".collect()" in block, "one task after every sample, not one per sample"
    # Merged MAFs are summarized too (the NA flag shows only there).
    assert "GBCMS_DNA_WF.out.merged" in main
    dna = (NEXTFLOW / "workflows/dna.nf").read_text()
    assert "merged   = ch_merged" in dna and "ch_merged   = MERGE_COUNTS.out.merged_maf" in dna


def test_the_parameter_is_documented_with_its_default():
    docs = (ROOT / "docs/nextflow/parameters.md").read_text()
    assert re.search(r"^\| `--qc_summary` \| `false` \|", docs, re.M)
    assert "cli/summarize.md" in docs
    assert "`--qc_summary`" in (NEXTFLOW / "README.md").read_text()


# ── the docs ──────────────────────────────────────────────────────────────────


def test_the_summarize_examples_in_the_docs_parse():
    group = get_command(app).commands["summarize"]  # type: ignore[attr-defined]
    for p in group.params:
        if getattr(p.type, "exists", False):
            p.type.exists = False  # type: ignore[attr-defined]
    seen = 0
    for page in (ROOT / "docs").rglob("*.md"):
        if page.name == "CHANGELOG.md":
            continue
        in_bash = False
        for line in page.read_text().splitlines():
            if line.strip().startswith("```"):
                in_bash = line.strip().startswith("```bash") and not in_bash
                continue
            stripped = line.strip()
            if in_bash and stripped.startswith("gbcms summarize") and "<" not in stripped:
                tokens = shlex.split(stripped, comments=True)
                cut = next((i for i, t in enumerate(tokens) if t in (">", "|")), len(tokens))
                try:
                    group.make_context("summarize", tokens[2:cut])
                except Exception as e:  # `--help` exits 0: the flag exists
                    if type(e).__name__ != "Exit":
                        raise
                seen += 1
    assert seen >= 4


def test_the_summarize_page_is_in_the_nav_and_the_index():
    assert "cli/summarize.md" in (ROOT / "mkdocs.yml").read_text()
    assert "(summarize.md)" in (ROOT / "docs/cli/index.md").read_text()
    for flag in (f for f in QC_FLAGS if f.kind is Kind.MERGE):
        assert flag.name in (ROOT / "docs/cli/summarize.md").read_text()
