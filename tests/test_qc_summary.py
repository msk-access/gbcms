"""O10 #227: the end-of-run QC summary.

Contracts:
- The summary reads the output a run has written and changes nothing: counts,
  columns, rows and files are byte-identical with and without it, and a failure to
  read the output is a warning, never a failed run.
- One table (``QC_FLAGS``) holds each flag's tier and action; the QC Flags page
  repeats it. Every flag the page defines has an entry and the reverse, and the page's
  tables and example block are the code's own.
- A flag is in one tier: fix the input, check the call, informational. The first two
  are the rows that need a look. COEXISTING_ALLELE is informational and
  WARN_REF_CORRECTED is "check the call" (operator, 2026-10-05).
- The block ends with a WARNING when any row needs a look and an INFO line when
  none does; zeros are shown for the first two tiers.
"""

import glob
import logging
import pathlib
import re

import pytest
from helpers import read_maf_output
from rna_fixtures import mk_ref, write_bam, write_fasta, write_gtf
from rna_fixtures import write_vcf as write_variant_vcf
from test_asjd2_contract import _retention_setup
from test_input_group4_contract import REF, P, _files, _maf, _ref_reads, _row
from typer.testing import CliRunner

from gbcms.cli import app
from gbcms.pipeline import Pipeline
from gbcms.qc_summary import (
    QC_FLAGS,
    QC_FLAGS_URL,
    Kind,
    Tier,
    log_qc_summary,
    render_summary,
    summarize_output,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
PAGE = ROOT / "docs" / "reference" / "qc-flags.md"
runner = CliRunner()

# ── Writers for synthetic outputs ─────────────────────────────────────────────

_MAF_COLS = [
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "gbcms_status",
    "gbcms_status_reason",
    "gbcms_diagnostic",
]


def maf_row(n=0, gene="GENE1", status="PASS", reason="", diag="", **extra):
    return {
        "Hugo_Symbol": gene,
        "Chromosome": "chr1",
        "Start_Position": str(1000 + n),
        "Reference_Allele": "A",
        "Tumor_Seq_Allele2": "T",
        "gbcms_status": status,
        "gbcms_status_reason": reason,
        "gbcms_diagnostic": diag,
        **extra,
    }


def write_maf(path, rows, extra_cols=()):
    """A MAF as the writer lays it out: the provenance line, a header, tab-separated rows."""
    cols = [*_MAF_COLS, *extra_cols]
    lines = ["#gbcms v0.0.0 (test)", "\t".join(cols)]
    lines += ["\t".join(str(r.get(c, "")) for c in cols) for r in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def vcf_row(n=0, gs="PASS", gsr=".", gd=".", flags=(), **info):
    fields = [f"DP=10;GS={gs};GSR={gsr};GD={gd}"]
    fields += [f"{k}={v}" for k, v in info.items()] + list(flags)
    return f"chr1\t{1000 + n}\t.\tA\tT\t.\t.\t{';'.join(fields)}\tGT:DP:AD\t0/0:10:10,0"


def write_vcf(path, rows):
    head = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    path.write_text(head + "\n".join(rows) + "\n")
    return path


def summarize(tmp_path, rows, **kw):
    return summarize_output(
        write_maf(tmp_path / "o.maf", rows, kw.pop("extra_cols", ())), "S", **kw
    )


# ── The table and the QC Flags page share one definition ──────────────────────

_SECTION = re.compile(r"--8<-- \[start:([\w-]+)\](.*?)--8<-- \[end:\1\]", re.S)
# The mFSD report classes live in the HTML report, not in an output row.
_NOT_IN_ROWS = {"mfsd-classes"}


def _page_flags():
    """The flags the page defines: the first cell of each row of a family table."""
    flags = set()
    for name, body in _SECTION.findall(PAGE.read_text()):
        if name in _NOT_IN_ROWS:
            continue
        for line in body.splitlines():
            if line.startswith("| ") and not line.startswith("|:"):
                m = re.match(r"`([A-Za-z][\w-]*)", line.split("|")[1].strip())
                if m:
                    flags.add(m.group(1))
    return flags - {"Reason", "Flag", "Outcome", "Column", "Class"}


def _summary_section():
    text = PAGE.read_text()
    start = text.index("\n## End-of-run summary")
    nxt = re.search(r"\n## ", text[start + 1 :])
    return text[start : start + 1 + nxt.start()] if nxt else text[start:]


def _page_tiers():
    """{flag: (tier, summary text)} from the page's End-of-run summary tables."""
    titles = {t.title: t for t in Tier}
    tier, out = None, {}
    for line in _summary_section().splitlines():
        if line.startswith("### "):
            tier = titles.get(line[4:].strip())
        elif tier and line.startswith("| `"):
            name, text = (c.strip() for c in line.strip().strip("|").split("|", 1))
            assert name.strip("`") not in out, f"{name} listed twice"
            out[name.strip("`")] = (tier, text.replace("`", ""))
    return out


def test_every_flag_the_page_defines_has_a_summary_entry_and_the_reverse():
    entries = {f.name for f in QC_FLAGS}
    page = _page_flags()
    assert len(page) > 30
    assert sorted(page - entries) == [], "defined on the page, missing from QC_FLAGS"
    assert sorted(entries - page) == [], "in QC_FLAGS, not defined on the page"


def test_the_page_lists_each_flag_with_the_code_tier_and_action():
    page = _page_tiers()
    assert set(page) == {f.name for f in QC_FLAGS}
    wrong = [f.name for f in QC_FLAGS if page[f.name] != (f.tier, f.action)]
    assert wrong == [], f"the page and the code table differ for {wrong}"


def test_each_flag_is_in_the_family_the_page_defines_it_under():
    """A flag's family is the section of the page it is defined in."""
    family = {"status": Kind.REASON, "diagnostics": Kind.DIAGNOSTIC, "rescue": Kind.RESCUE}
    family |= {"asjd": Kind.ASJD, "qc-columns": Kind.COLUMN, "merge": Kind.MERGE}
    defined = {}
    for name, body in _SECTION.findall(PAGE.read_text()):
        if name in family:
            for line in body.splitlines():
                m = re.match(r"\| `([A-Za-z][\w-]*)", line)
                if m:
                    defined[m.group(1)] = family[name]
    wrong = [f.name for f in QC_FLAGS if defined[f.name] is not f.kind]
    assert wrong == []


def test_every_family_anchor_is_a_heading_of_the_page():
    slugs = {
        re.sub(r"[^\w\s-]", "", h.lower()).strip().replace(" ", "-")
        for h in re.findall(r"^## (.+)$", PAGE.read_text(), re.M)
    }
    assert {k.value for k in Kind} <= slugs


def test_the_documented_url_is_the_pages_address():
    site = re.search(r"^site_url:\s*(\S+)", (ROOT / "mkdocs.yml").read_text(), re.M).group(1)
    assert QC_FLAGS_URL == site.rstrip("/") + "/reference/qc-flags/"


def test_actions_fit_a_table_cell_and_a_log_line():
    for f in QC_FLAGS:
        assert f.action and f.action == f.action.strip(), f.name
        assert not set(f.action) & set("|[]`\n"), f.name
        # One line in the block for an informational gloss; two or three wrapped for the rest.
        assert len(f.action) <= (46 if f.tier is Tier.INFORMATIONAL else 175), f.name


def test_the_operator_decided_tiers():
    tier = {f.name: f.tier for f in QC_FLAGS}
    assert tier["COEXISTING_ALLELE"] is Tier.INFORMATIONAL
    assert tier["WARN_REF_CORRECTED"] is Tier.CHECK_CALL
    fail_reasons = (
        "REF_MISMATCH FETCH_FAILED EMPTY_ALLELE NON_SEQUENCE_ALLELE ALT_EQUALS_REF "
        "ALT_CONTAINS_N REF_AT_OFFSET"
    ).split()
    assert all(tier[r] is Tier.FIX_INPUT for r in fail_reasons)
    for name in "ZERO_ALT MULTI_ALLELIC TRACT_CLUSTER MNP_DISC_RATIO MNP_RESCUE_ELIGIBLE".split():
        assert tier[name] is Tier.INFORMATIONAL, name
    for name in "LOW_ALT_JUNC LOW_REF_JUNC".split():
        assert tier[name] is Tier.INFORMATIONAL, name
    check = (
        "OBSERVED_ALLELE PARTIAL_DOMINANT CLIP_CANDIDATES SW_FALLBACK RESCUED_COMPONENT "
        "SPLICE_SKIP_DOMINANT NON_DISCRIMINATING_LOCUS RETENTION_DOMINANT "
        "NOVEL_JUNC_AT_SPLICE_LOSS rna_editing_site"
    ).split()
    assert all(tier[c] is Tier.CHECK_CALL for c in check)


# ── The example on the page is the block the code renders ─────────────────────


def example_rows():
    """2,000 synthetic rows: two rejected, 38 to check, 107 with an informational flag only."""
    recipe = [  # (rows, status, reason, diagnostic)
        (1, "FAIL", "REF_MISMATCH", "REF_AT_OFFSET(2)"),
        (1, "FAIL", "REF_MISMATCH", ""),
        (24, "PASS", "", "PARTIAL_DOMINANT"),
        (6, "PASS", "", "PARTIAL_DOMINANT;COEXISTING_ALLELE(chr1:2000:A>AT:6/2)"),
        (7, "PASS", "", "OBSERVED_ALLELE(chr1:3000:C>T:9/0)"),
        (1, "PASS", "", "CLIP_CANDIDATES(3)"),
        (40, "PASS", "MULTI_ALLELIC", ""),
        (55, "PASS", "", "ZERO_ALT"),
        (12, "PASS", "", "COEXISTING_ALLELE(chr1:2000:A>AT:6/2)"),
    ]
    flags = [(status, reason, diag) for count, status, reason, diag in recipe for _ in range(count)]
    flags += [("PASS", "", "")] * (2000 - len(flags))
    return [
        maf_row(n, gene=f"GENE{n % 7 + 1}", status=status, reason=reason, diag=diag)
        for n, (status, reason, diag) in enumerate(flags)
    ]


def test_the_example_on_the_page_is_the_rendered_block(tmp_path):
    expected = render_summary([summarize(tmp_path, example_rows())])
    m = re.search(r"<!-- qc-summary-example -->\n```text\n(.*?)\n```", PAGE.read_text(), re.S)
    assert m, "the page has no <!-- qc-summary-example --> block"
    assert m.group(1) == expected, "the example is stale; the rendered block is:\n" + expected


# ── Tallying synthetic outputs ────────────────────────────────────────────────


def test_a_clean_output_needs_no_look_and_shows_the_zeros(tmp_path, caplog):
    s = summarize(tmp_path, [maf_row(i) for i in range(3)])
    assert (s.rows, s.need_look, s.info_only) == (3, 0, 0)
    text = render_summary([s])
    assert text.splitlines()[0] == "QC summary for sample S: 0 of 3 rows need a look (0.0%)"
    for tier in ("Fix the input", "Check the call", "Informational"):
        assert f"  {tier}: none" in text
    with caplog.at_level(logging.INFO, logger="gbcms.qc_summary"):
        log_qc_summary([s])
    last = caplog.records[-1]
    assert last.levelno == logging.INFO and "No row needs a look" in last.getMessage()
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]


def test_rows_are_counted_once_per_tier_and_once_in_need_a_look(tmp_path):
    s = summarize(
        tmp_path,
        [
            maf_row(0, status="FAIL", reason="REF_MISMATCH", diag="REF_AT_OFFSET(2)"),
            maf_row(1, reason="WARN_REF_CORRECTED|MULTI_ALLELIC"),
            maf_row(
                2, diag="PARTIAL_DOMINANT;COEXISTING_ALLELE(1:100:A>T:5/2);MNP_DISC_RATIO(2/5)"
            ),
            maf_row(3, diag="ZERO_ALT"),
            maf_row(4),
        ],
    )
    assert (s.rows, s.need_look, s.info_only) == (5, 3, 1)
    assert dict(s.flag_rows) == {
        "REF_MISMATCH": 1,
        "REF_AT_OFFSET": 1,
        "WARN_REF_CORRECTED": 1,
        "MULTI_ALLELIC": 1,
        "PARTIAL_DOMINANT": 1,
        "COEXISTING_ALLELE": 1,
        "MNP_DISC_RATIO": 1,
        "ZERO_ALT": 1,
    }
    assert dict(s.tier_rows) == {Tier.FIX_INPUT: 1, Tier.CHECK_CALL: 2, Tier.INFORMATIONAL: 3}
    assert not s.unlisted


def test_a_parametric_flag_counts_under_its_base_name(tmp_path):
    s = summarize(
        tmp_path,
        [maf_row(0, diag="HIGH_N_FRACTION(0.08);MNP_DISC_RATIO(1/3);SPLICE_SKIP_DOMINANT(12)")],
    )
    assert set(s.flag_rows) == {"HIGH_N_FRACTION", "MNP_DISC_RATIO", "SPLICE_SKIP_DOMINANT"}


def test_a_row_names_its_check_flags_in_full_and_not_its_informational_ones(tmp_path):
    s = summarize(
        tmp_path,
        [maf_row(0, gene="GENE9", diag="OBSERVED_ALLELE(1:100:A>T:5/0);ZERO_ALT")],
    )
    assert s.examples == [("chr1:1000 A>T (GENE9)", "OBSERVED_ALLELE(1:100:A>T:5/0)")]


def test_examples_put_input_problems_first_and_stop_at_five(tmp_path):
    rows = [maf_row(i, diag="PARTIAL_DOMINANT") for i in range(4)]
    rows += [maf_row(10 + i, status="FAIL", reason="EMPTY_ALLELE") for i in range(3)]
    s = summarize(tmp_path, rows)
    labels = [label for label, _ in s.examples]
    assert labels == [f"chr1:{1010 + i} A>T (GENE1)" for i in range(3)] + [
        "chr1:1000 A>T (GENE1)",
        "chr1:1001 A>T (GENE1)",
    ]
    assert s.need_look == 7
    text = render_summary([s])
    assert "Rows to look at (first 5 of 7, input problems first):" in text
    assert summarize(tmp_path, rows, max_examples=0).examples == []


def test_the_label_drops_an_unknown_gene_and_reads_allele1_when_allele2_is_the_ref(tmp_path):
    rows = [
        maf_row(0, gene="Unknown", diag="PARTIAL_DOMINANT"),
        maf_row(1, diag="PARTIAL_DOMINANT", Tumor_Seq_Allele2="A", Tumor_Seq_Allele1="G"),
    ]
    s = summarize(tmp_path, rows, extra_cols=["Tumor_Seq_Allele1"])
    assert [label for label, _ in s.examples] == ["chr1:1000 A>T", "chr1:1001 A>G (GENE1)"]


def test_a_flag_the_table_does_not_know_is_listed_apart_and_not_a_look(tmp_path):
    s = summarize(
        tmp_path,
        [
            maf_row(0, diag="NEW_FLAG(3);ZERO_ALT"),
            maf_row(1, diag="NEW_FLAG"),
            maf_row(2, status="FAIL", reason="NEW_REASON"),
        ],
    )
    assert dict(s.unlisted) == {"NEW_FLAG": 2, "NEW_REASON": 1}
    assert s.need_look == 1, "only the FAIL row; a FAIL verdict is an input problem"
    assert s.examples == [("chr1:1002 A>T (GENE1)", "NEW_REASON")]
    text = render_summary([s])
    assert "Not in this build's flag table" in text and "NEW_FLAG (2 rows)" in text
    assert s.unknown_only == 1
    assert (
        "  Rows with only informational flags: 1; with no flag: 0; "
        "with only flags this build does not know: 1"
    ) in text


def test_asjd_flags_and_qc_columns(tmp_path):
    cols = [
        "asjd_diagnostic",
        "rna_editing_site",
        "asjd_flag",
        "mfsd_ks_valid",
        "mfsd_alt_confidence",
        "mfsd_ch_flag",
    ]
    rows = [
        maf_row(0, asjd_diagnostic="LOW_ALT_JUNC;RETENTION_DOMINANT(65)", rna_editing_site="True"),
        maf_row(1, asjd_diagnostic="NOVEL_JUNC_AT_SPLICE_LOSS(5@100-200)", asjd_flag="True"),
        maf_row(2, mfsd_ks_valid="False", mfsd_alt_confidence="SPARSE", mfsd_ch_flag="True"),
        maf_row(3, mfsd_ks_valid="True", mfsd_alt_confidence="TESTABLE", mfsd_ch_flag="False"),
        maf_row(4, mfsd_alt_confidence="NONE", rna_editing_site="False", asjd_flag="False"),
    ]
    s = summarize(tmp_path, rows, extra_cols=cols)
    assert dict(s.flag_rows) == {
        "LOW_ALT_JUNC": 1,
        "RETENTION_DOMINANT": 1,
        "rna_editing_site": 1,
        "NOVEL_JUNC_AT_SPLICE_LOSS": 1,
        "asjd_flag": 1,
        "mfsd_ks_valid": 1,
        "mfsd_alt_confidence": 2,
        "mfsd_ch_flag": 1,
    }
    assert s.need_look == 2 and s.info_only == 2


def test_a_maf_without_the_qc_columns_has_no_column_flags(tmp_path):
    """DNA output has no RNA or mFSD columns; their absence is not a flag."""
    s = summarize(tmp_path, [maf_row(0), maf_row(1)])
    assert not s.flag_rows


_RESCUED = (
    "method=decomposed;outcome=rescued;original_ref=9;original_alt=0;original_partial=4;"
    "original_confirmed=0;adopted=chr1:100(A>T)"
)
_KEPT = _RESCUED.replace("rescued", "no_improvement").split(";adopted")[0]


def test_rescue_outcomes_are_read_from_the_audit_column(tmp_path):
    rows = [
        maf_row(0, diag="RESCUED_COMPONENT(chr1:100:A>T)", gbcms_rescue=_RESCUED),
        maf_row(1, gbcms_rescue=_KEPT),
        maf_row(2, gbcms_rescue=_KEPT.replace("no_improvement", "ref_validation_failed")),
        maf_row(3),
    ]
    s = summarize(tmp_path, rows, extra_cols=["gbcms_rescue"])
    assert dict(s.flag_rows) == {
        "RESCUED_COMPONENT": 1,
        "rescued": 1,
        "no_improvement": 1,
        "ref_validation_failed": 1,
    }
    assert s.need_look == 2 and s.info_only == 1
    # An outcome is named as the audit column writes it.
    text = render_summary([s])
    assert "    outcome=rescued: 1 row" in text and "outcome=no_improvement: 1 row (" in text
    assert ("chr1:1000 A>T (GENE1)", "RESCUED_COMPONENT(chr1:100:A>T), outcome=rescued") in (
        s.examples
    )


def test_every_rescue_outcome_the_audit_defines_has_an_entry():
    from gbcms import rescue_audit

    outcomes = {v for k, v in vars(rescue_audit).items() if k.startswith("OUTCOME_")}
    entries = {f.name for f in QC_FLAGS if f.kind is Kind.RESCUE}
    assert outcomes == entries


def test_an_outcome_the_table_does_not_know_is_listed_apart(tmp_path):
    s = summarize(
        tmp_path,
        [maf_row(0, gbcms_rescue=_KEPT.replace("no_improvement", "new_outcome"))],
        extra_cols=["gbcms_rescue"],
    )
    assert dict(s.unlisted) == {"new_outcome": 1} and s.need_look == 0


def test_the_reader_skips_provenance_lines_and_reads_quoted_cells(tmp_path):
    path = tmp_path / "q.maf"
    path.write_text(
        "#gbcms v0\n#command gbcms dna -v x\n"
        + "\t".join(_MAF_COLS)
        + '\nGENE1\tchr1\t5\tA\tT\tPASS\t\t"PARTIAL_DOMINANT;ZERO_ALT"\n'
    )
    s = summarize_output(path, "S")
    assert s.rows == 1 and set(s.flag_rows) == {"PARTIAL_DOMINANT", "ZERO_ALT"}


def test_an_output_that_is_not_gbcms_is_refused(tmp_path):
    path = tmp_path / "x.maf"
    path.write_text("Chromosome\tStart_Position\nchr1\t5\n")
    with pytest.raises(ValueError, match="gbcms_status"):
        summarize_output(path, "S")
    empty = tmp_path / "e.maf"
    empty.write_text("")
    with pytest.raises(ValueError):
        summarize_output(empty, "S")


def test_summarizing_changes_nothing_on_disk(tmp_path):
    path = write_maf(tmp_path / "o.maf", [maf_row(0, diag="PARTIAL_DOMINANT")])
    before = (path.read_bytes(), path.stat().st_mtime_ns)
    summarize_output(path, "S")
    assert (path.read_bytes(), path.stat().st_mtime_ns) == before
    assert sorted(p.name for p in tmp_path.iterdir()) == ["o.maf"]


# ── VCF output ────────────────────────────────────────────────────────────────


def test_vcf_flags_are_read_from_the_info_field(tmp_path):
    rows = [
        vcf_row(0, gs="FAIL", gsr="REF_MISMATCH", gd="REF_AT_OFFSET(1)"),
        vcf_row(1, gsr="WARN_REF_CORRECTED|MULTI_ALLELIC", gd="PARTIAL_DOMINANT|ZERO_ALT"),
        vcf_row(2, ASJDD="RETENTION_DOMINANT(65)|LOW_ALT_JUNC", flags=["RED", "ASJD"]),
        vcf_row(3, GR=_RESCUED.replace(";", "|")),
        vcf_row(4),
    ]
    path = write_vcf(tmp_path / "o.vcf", rows)
    s = summarize_output(path, "S")
    assert dict(s.flag_rows) == {
        "REF_MISMATCH": 1,
        "REF_AT_OFFSET": 1,
        "WARN_REF_CORRECTED": 1,
        "MULTI_ALLELIC": 1,
        "PARTIAL_DOMINANT": 1,
        "ZERO_ALT": 1,
        "RETENTION_DOMINANT": 1,
        "LOW_ALT_JUNC": 1,
        "rna_editing_site": 1,
        "asjd_flag": 1,
        "rescued": 1,
    }
    assert (s.fmt, s.rows, s.need_look) == ("vcf", 5, 4)
    assert [label for label, _ in s.examples][0] == "chr1:1000 A>T"
    assert "  Find them in the output (VCF INFO keys): ASJDD, GD, GR, GS, GSR, RED" in (
        render_summary([s])
    )


def test_a_short_vcf_line_is_skipped(tmp_path):
    path = write_vcf(tmp_path / "o.vcf", [vcf_row(0, gd="ZERO_ALT"), "chr1\t5\t.\tA"])
    assert summarize_output(path, "S").rows == 1


# ── Several samples ───────────────────────────────────────────────────────────


def test_several_samples_are_combined_and_the_worst_are_named(tmp_path):
    sums = []
    for i in range(12):
        rows = [maf_row(j, diag="PARTIAL_DOMINANT") for j in range(i % 4)] + [maf_row(99)]
        sums.append(summarize_output(write_maf(tmp_path / f"s{i}.maf", rows), f"S{i}"))
    text = render_summary(sums)
    lines = text.splitlines()
    need = sum(s.need_look for s in sums)
    assert lines[0] == f"QC summary for 12 samples: {need} of 30 rows need a look ({need / 30:.1%})"
    assert "    PARTIAL_DOMINANT: 18 rows in 9 of 12 samples" in lines
    assert "  Samples with rows to look at (first 9 of 9):" in lines
    assert "    S3: 3 of 4 rows" in lines
    assert "Rows to look at (first 5 of 18):" in text
    assert "    S1 chr1:1000 A>T (GENE1): PARTIAL_DOMINANT" in lines


def test_the_sample_list_stops_at_ten(tmp_path):
    sums = [
        summarize_output(
            write_maf(tmp_path / f"s{i}.maf", [maf_row(0, diag="SW_FALLBACK(2)")]), f"S{i}"
        )
        for i in range(13)
    ]
    assert "  Samples with rows to look at (first 10 of 13):" in render_summary(sums)


# ── The closing line ──────────────────────────────────────────────────────────


def test_the_run_ends_on_a_warning_when_rows_need_a_look(tmp_path, caplog):
    s = summarize(tmp_path, [maf_row(0, diag="PARTIAL_DOMINANT"), maf_row(1, status="FAIL")])
    with caplog.at_level(logging.INFO, logger="gbcms.qc_summary"):
        log_qc_summary([s])
    block, last = caplog.records[-2:]
    assert block.levelno == logging.INFO and block.getMessage().startswith(
        "QC summary for sample S"
    )
    assert last.levelno == logging.WARNING
    assert last.getMessage() == (
        f"2 of 2 rows need a look (1 with a problem in the input): see the QC summary above "
        f"and {QC_FLAGS_URL}"
    )
    # User text (a sample name, a gene) is not Rich markup.
    assert block.markup is False and last.markup is False  # type: ignore[attr-defined]


def test_the_closing_line_of_several_samples_names_how_many_have_rows(tmp_path, caplog):
    sums = [
        summarize_output(write_maf(tmp_path / "a.maf", [maf_row(0, diag="SW_FALLBACK(2)")]), "A"),
        summarize_output(write_maf(tmp_path / "b.maf", [maf_row(0), maf_row(1)]), "B"),
    ]
    with caplog.at_level(logging.INFO, logger="gbcms.qc_summary"):
        log_qc_summary(sums)
    last = caplog.records[-1]
    assert last.levelno == logging.WARNING
    assert last.getMessage().startswith(
        "1 of 3 rows needs a look in 1 of 2 samples (0 with a problem in the input): "
    )


def test_nothing_is_logged_for_a_run_that_wrote_no_sample(caplog):
    with caplog.at_level(logging.DEBUG, logger="gbcms.qc_summary"):
        log_qc_summary([])
    assert caplog.records == []


# ── Through the CLI ───────────────────────────────────────────────────────────


def _dna(tmp_path, maf, bam, fa, fmt="maf", name="out", sample="S"):
    res = runner.invoke(
        app,
        ["dna", "-v", str(maf), "-b", f"{sample}:{bam}", "-f", str(fa), "-o", str(tmp_path / name)]
        + ["--format", fmt],
        env={"COLUMNS": "3000"},
    )
    assert res.exit_code == 0, res.output
    return res, [ln for ln in res.output.splitlines() if ln.strip()]


def _rows(tmp_path):
    other = "A" if REF[P] != "A" else "C"
    wrong = "A" if REF[P + 19] != "A" else "C"
    return [
        _row(P + 1, REF[P], other),  # no read carries it: ZERO_ALT
        _row(P + 20, wrong, "T"),  # the REF is not the reference's: REF_MISMATCH
    ]


def test_a_run_ends_with_its_summary_and_a_warning(tmp_path):
    fa, bam = _files(tmp_path, _ref_reads())
    res, lines = _dna(tmp_path, _maf(tmp_path, _rows(tmp_path)), bam, fa)
    text = res.output
    assert "QC summary for sample S: 1 of 2 rows needs a look (50.0%)" in text
    assert "    REF_MISMATCH: 1 row" in text and "    ZERO_ALT: 1 row" in text
    assert "Fix the input: 1 row" in text and "Check the call: none" in text
    assert "WARNING" in lines[-1] and "1 of 2 rows needs a look" in lines[-1]
    assert text.index("Pipeline completed") < text.index("QC summary for sample")
    # The run's counts are the engine's: the invariants hold on every row written.
    for r in read_maf_output(tmp_path / "out" / "S.maf"):
        assert int(r["total_count"]) >= int(r["ref_count"]) + int(r["alt_count"])
        assert int(r["total_count_fragment"]) >= int(r["ref_count_fragment"]) + int(
            r["alt_count_fragment"]
        )
        assert int(r["ref_count"]) == int(r["ref_count_forward"]) + int(r["ref_count_reverse"])
        assert int(r["alt_count"]) == int(r["alt_count_forward"]) + int(r["alt_count_reverse"])


def test_a_run_with_nothing_to_look_at_ends_on_info(tmp_path):
    fa, bam = _files(tmp_path, _ref_reads())
    res, lines = _dna(tmp_path, _maf(tmp_path, _rows(tmp_path)[:1]), bam, fa, name="ok")
    assert "WARNING" not in lines[-1] and "No row needs a look: 1 row checked" in lines[-1]
    assert "  Fix the input: none" in res.output


def test_a_sample_name_with_brackets_is_logged_as_written(tmp_path):
    fa, bam = _files(tmp_path, _ref_reads())
    res, _ = _dna(tmp_path, _maf(tmp_path, _rows(tmp_path)), bam, fa, sample="S[x]")
    assert "QC summary for sample S[x]:" in res.output


def test_two_samples_get_one_combined_block(tmp_path):
    fa, bam = _files(tmp_path, _ref_reads())
    out = tmp_path / "multi"
    res = runner.invoke(
        app,
        ["dna", "-v", str(_maf(tmp_path, _rows(tmp_path))), "-f", str(fa), "-o", str(out)]
        + ["-b", f"A:{bam}", "-b", f"B:{bam}"],
        env={"COLUMNS": "3000"},
    )
    assert res.exit_code == 0, res.output
    assert res.output.count("QC summary for") == 1
    assert "QC summary for 2 samples: 2 of 4 rows need a look (50.0%)" in res.output
    assert "    REF_MISMATCH: 2 rows in 2 of 2 samples" in res.output


@pytest.mark.parametrize("fmt", ["maf", "vcf"])
def test_the_output_files_are_byte_identical_with_and_without_the_summary(
    tmp_path, monkeypatch, fmt
):
    fa, bam = _files(tmp_path, _ref_reads())
    maf = _maf(tmp_path, _rows(tmp_path))

    def snapshot():
        files = sorted(glob.glob(str(tmp_path / "out" / "*")))
        return {pathlib.Path(f).name: pathlib.Path(f).read_bytes() for f in files}

    with_summary, _ = _dna(tmp_path, maf, bam, fa, fmt)
    assert "QC summary for" in with_summary.output
    first = snapshot()
    for f in glob.glob(str(tmp_path / "out" / "*")):
        pathlib.Path(f).unlink()

    monkeypatch.setattr(Pipeline, "_tally_qc", lambda self, sample_name: None)
    without, _ = _dna(tmp_path, maf, bam, fa, fmt)
    assert "QC summary" not in without.output
    assert first and snapshot() == first


def test_a_summary_that_cannot_read_its_output_warns_and_the_run_succeeds(tmp_path, monkeypatch):
    fa, bam = _files(tmp_path, _ref_reads())
    maf = _maf(tmp_path, _rows(tmp_path))
    _dna(tmp_path, maf, bam, fa)

    def broken(*args, **kwargs):
        raise OSError("disk went away")

    monkeypatch.setattr("gbcms.pipeline.summarize_output", broken)
    res, lines = _dna(tmp_path, maf, bam, fa, name="broken")
    assert (
        "QC summary skipped for sample S (OSError: disk went away); the output is unaffected"
        in (" ".join(res.output.split()))
    )
    assert "QC summary for" not in res.output
    # The provenance lines carry the output directory; every row is as in the plain run.
    got = list(read_maf_output(tmp_path / "broken" / "S.maf"))
    assert got == list(read_maf_output(tmp_path / "out" / "S.maf"))
    assert sorted(p.name for p in (tmp_path / "broken").iterdir()) == ["S.maf"]


def test_an_rna_run_names_its_splice_flags_in_the_summary(tmp_path):
    """RNA mode, VCF output (the RNA default): ASJD flags are read from ASJDD."""
    ref = mk_ref()
    rows, reads = _retention_setup(ref)
    out = tmp_path / "rna"
    res = runner.invoke(
        app,
        ["rna", "-v", str(write_variant_vcf(tmp_path, rows)), "-f", str(write_fasta(tmp_path, ref))]
        + ["-b", f"S:{write_bam(tmp_path, ref, reads)}", "-o", str(out)]
        + ["--gtf", str(write_gtf(tmp_path))],
        env={"COLUMNS": "3000"},
    )
    assert res.exit_code == 0, res.output
    assert glob.glob(str(out / "*.vcf"))
    text = res.output
    assert "QC summary for sample S: 1 of 1 row needs a look (100.0%)" in text
    assert "    RETENTION_DOMINANT: 1 row" in text
    assert "VCF INFO keys): ASJDD" in text


def test_the_run_docs_point_at_the_summary():
    for page in ("cli/dna.md", "cli/rna.md", "resources/troubleshooting.md"):
        text = (ROOT / "docs" / page).read_text()
        assert "reference/qc-flags.md#end-of-run-summary" in text, page
