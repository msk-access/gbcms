"""The QC summary of finished outputs, in every format ``gbcms summarize`` writes.

One path makes all of them: :func:`summarize_paths` reads each output through
:mod:`gbcms.qc_summary` and every format is a function of the resulting
:class:`~gbcms.qc_summary.QcSummary` list.

* ``text``: the block the run logs (``render_summary``) and its closing line.
* ``json``: one ``<sample>.qc_summary.json`` per output, every flag with its zeros.
* ``tsv``: one cohort table, a row per sample.
* ``multiqc``: MultiQC custom content (``*_mqc.json``): a per-sample table coloured by
  tier and, only when a row needs a look, a "WARNING" section.

Nothing here writes into an input file.
"""

import csv
import io
import json
import logging
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from .io.atomic import atomic_output
from .io.output import build_identity
from .models.core import StrEnum
from .qc_summary import (
    MAX_EXAMPLES,
    QC_FLAGS,
    QC_FLAGS_URL,
    QcSummary,
    Skipped,
    Tier,
    _label,
    _not_checked_note,
    closing_line,
    render_summary,
    sample_name,
    summarize_output,
)

logger = logging.getLogger(__name__)

# The functions and the format enum. The module's constants (the file names, the schema and
# section ids) are imported by name: the QC flags docs test reads every upper-case string
# literal in the source as a flag the code emits, so listing them here would document them
# as flags.
__all__ = [
    "SummaryFormat",
    "multiqc_documents",
    "multiqc_table",
    "multiqc_warning",
    "render_tsv",
    "summarize_inputs",
    "summarize_paths",
    "summary_from_dict",
    "summary_to_dict",
    "text_report",
    "write_summaries",
]

JSON_SCHEMA = "gbcms.qc_summary/1"

JSON_SUFFIX = ".qc_summary.json"
TSV_NAME = "gbcms_qc_summary.tsv"
MULTIQC_TABLE_NAME = "gbcms_qc_summary_mqc.json"
MULTIQC_WARNING_NAME = "gbcms_qc_warning_mqc.json"

# MultiQC section ids: a pipeline orders its report by them (``custom_content: order``).
MULTIQC_TABLE_ID = "gbcms_qc_summary"
MULTIQC_WARNING_ID = "gbcms_qc_warning"


class SummaryFormat(StrEnum):
    """The formats of ``gbcms summarize --format``."""

    TEXT = "text"
    JSON = "json"
    TSV = "tsv"
    MULTIQC = "multiqc"


def summarize_inputs(
    paths: Sequence[Path | str], keep_rows: bool = False
) -> tuple[list[QcSummary], list[Skipped]]:
    """The summary of each output file that can be read, in the order given, and the
    inputs that could not be: refused (not a ``.maf`` or ``.vcf``), unreadable, not a
    gbcms output, or sharing a sample name (the file name without ``.maf`` or ``.vcf``)
    with an earlier input, which the cohort table and the per-sample files could not
    tell apart. A skipped input is named with its reason and is never dropped silently."""
    summaries: list[QcSummary] = []
    skipped: list[Skipped] = []
    taken: dict[str, Path] = {}
    for raw in paths:
        path = Path(raw)
        name = sample_name(path)
        if path.suffix not in (".maf", ".vcf"):
            reason = f"{path.name}: expected a .maf or .vcf output (compressed files are not read)"
            skipped.append(Skipped(name, reason))
        elif name in taken:
            reason = (
                f"two outputs share the sample name {name}: {path} is not checked, "
                f"{taken[name]} is (rename one, or summarize them separately)"
            )
            skipped.append(Skipped(name, reason))
        else:
            try:
                summaries.append(summarize_output(path, name, keep_rows=keep_rows))
            except (OSError, ValueError, csv.Error) as e:
                reason = str(e) if isinstance(e, ValueError) else f"{type(e).__name__}: {e}"
                skipped.append(Skipped(name, reason))
                continue
            taken[name] = path
    return summaries, skipped


def summarize_paths(paths: Sequence[Path | str], keep_rows: bool = False) -> list[QcSummary]:
    """:func:`summarize_inputs` for a caller that wants all or nothing: the first input
    that is skipped raises ``ValueError`` with its reason."""
    summaries, skipped = summarize_inputs(paths, keep_rows)
    if skipped:
        raise ValueError(skipped[0].reason)
    return summaries


# ── JSON ──────────────────────────────────────────────────────────────────────


def _tier_key(tier: Tier) -> str:
    return tier.name.lower()


def summary_to_dict(summary: QcSummary) -> dict[str, Any]:
    """One sample's summary as plain data: every flag with its zeros (a zero is
    evidence), the tier and action of each, the rows that need a look, and where the
    output came from."""
    rows = summary.rows
    look = summary.look_rows or [(label, "", [flags]) for label, flags in summary.examples]
    return {
        "schema": JSON_SCHEMA,
        "sample": summary.sample,
        "format": summary.fmt,
        "merged": summary.merged,
        "written_by": summary.written_by,
        "summarized_by": build_identity(),
        "docs": QC_FLAGS_URL,
        "rows": rows,
        "need_look": summary.need_look,
        "need_look_share": round(summary.need_look / rows, 6) if rows else 0.0,
        "rows_by_tier": {_tier_key(t): summary.tier_rows[t] for t in Tier},
        "rows_with_only_informational_flags": summary.info_only,
        "rows_with_only_unlisted_flags": summary.unknown_only,
        "rows_with_no_flag": rows - summary.need_look - summary.info_only - summary.unknown_only,
        "flags": [
            {
                "name": f.name,
                "family": f.kind.name.lower(),
                "tier": _tier_key(f.tier),
                "rows": summary.flag_rows[f.name],
                "share": round(summary.flag_rows[f.name] / rows, 6) if rows else 0.0,
                "action": f.action,
                "docs": f"{_url}#{f.kind.value}",
            }
            for f in QC_FLAGS
        ],
        "unlisted_flags": dict(summary.unlisted),
        "rows_to_look_at_complete": len(look) == summary.need_look,
        "rows_to_look_at": [
            {"variant": variant, "gene": gene, "flags": flags} for variant, gene, flags in look
        ],
    }


def summary_from_dict(data: dict[str, Any]) -> QcSummary:
    """The :class:`QcSummary` a :func:`summary_to_dict` file describes, so a written
    summary renders as the run's block does."""
    if data.get("schema") != JSON_SCHEMA:
        raise ValueError(f"not a gbcms QC summary (schema {data.get('schema')!r})")
    summary = QcSummary(
        sample=data["sample"],
        fmt=data["format"],
        rows=data["rows"],
        need_look=data["need_look"],
        info_only=data["rows_with_only_informational_flags"],
        unknown_only=data["rows_with_only_unlisted_flags"],
        written_by=data["written_by"],
        merged=data["merged"],
    )
    summary.flag_rows.update({f["name"]: f["rows"] for f in data["flags"] if f["rows"]})
    by_key = {_tier_key(t): t for t in Tier}
    summary.tier_rows.update({by_key[k]: n for k, n in data["rows_by_tier"].items() if n})
    summary.unlisted.update(data["unlisted_flags"])
    summary.look_rows = [(r["variant"], r["gene"], r["flags"]) for r in data["rows_to_look_at"]]
    summary.examples = [
        (_label(v, g), ", ".join(f)) for v, g, f in summary.look_rows[:MAX_EXAMPLES]
    ]
    return summary


# ── TSV ───────────────────────────────────────────────────────────────────────

_TSV_FIXED = [
    "sample",
    "format",
    "rows",
    "need_look",
    "need_look_share",
    *(_tier_key(t) for t in Tier),
    "only_informational",
    "only_unlisted",
]


def render_tsv(summaries: Sequence[QcSummary]) -> str:
    """The cohort table: a row per sample, the rows by tier, then a column per flag (the
    order of the flag table, zeros included) and the flags the table does not know."""
    out = io.StringIO()
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow([*_TSV_FIXED, *(f.name for f in QC_FLAGS), "unlisted"])
    for s in summaries:
        share = f"{s.need_look / s.rows:.4f}" if s.rows else "0.0000"
        writer.writerow(
            [
                s.sample,
                s.fmt,
                s.rows,
                s.need_look,
                share,
                *(s.tier_rows[t] for t in Tier),
                s.info_only,
                s.unknown_only,
                *(s.flag_rows[f.name] for f in QC_FLAGS),
                ";".join(f"{name}:{n}" for name, n in s.unlisted.most_common()),
            ]
        )
    return out.getvalue()


# ── MultiQC custom content ────────────────────────────────────────────────────

_url = QC_FLAGS_URL
_HELP = (
    "Every row of a gbcms output is put in a tier by its QC flags. **Fix the input**: the row "
    "was rejected or its input is wrong. **Check the call**: it was counted, but its count or "
    "VAF may not mean what it looks like. Both tiers *need a look*. **Informational**: a caveat "
    f"or a description. Each flag, and what to do about it, is defined at {_url}"
)


def _tier_headers(only_need: bool = False) -> dict[str, dict[str, Any]]:
    """Table columns for the tier counts; zero is green, a first-tier row red, a
    second-tier row amber."""
    headers: dict[str, dict[str, Any]] = {
        "need_look": {
            "title": "Need a look",
            "description": "Rows in the first two tiers (fix the input, check the call)",
            "cond_formatting_rules": {"pass": [{"eq": 0}], "warn": [{"gt": 0}]},
        },
        "fix_input": {
            "title": "Fix the input",
            "description": "Rows rejected, or whose input is wrong",
            "cond_formatting_rules": {"pass": [{"eq": 0}], "fail": [{"gt": 0}]},
        },
        "check_call": {
            "title": "Check the call",
            "description": "Rows counted, but whose count or VAF may not mean what it looks like",
            "cond_formatting_rules": {"pass": [{"eq": 0}], "warn": [{"gt": 0}]},
        },
    }
    if not only_need:
        headers["informational"] = {
            "title": "Informational",
            "description": "Rows with a caveat or a description; nothing to do",
        }
    for header in headers.values():
        header.update({"format": "{:,.0f}", "scale": False})
    return headers


def multiqc_table(summaries: Sequence[QcSummary]) -> dict[str, Any]:
    """The per-sample table: rows, the rows by tier, and a column per flag that fired
    (hidden until switched on in the table's column list)."""
    headers: dict[str, dict[str, Any]] = {
        "rows": {
            "title": "Rows",
            "description": "Rows in the output",
            "format": "{:,.0f}",
            "scale": False,
        },
        **_tier_headers(),
        "need_look_pct": {
            "title": "% need a look",
            "description": "Share of rows in the first two tiers",
            "suffix": "%",
            "format": "{:,.1f}",
            "scale": False,
        },
    }
    # Which flags fired anywhere, worked out once: a flag that never fires costs a scan of
    # the whole cohort to find out, so asking again for each sample would cost its square.
    fired = [f for f in QC_FLAGS if any(s.flag_rows[f.name] for s in summaries)]
    for f in fired:
        headers[f.name] = {
            "title": f.name,
            "description": f"{f.tier.title}: {f.action}",
            "format": "{:,.0f}",
            "scale": False,
            "hidden": True,
        }
    data: dict[str, dict[str, float]] = {}
    for s in summaries:
        data[s.sample] = {
            "rows": s.rows,
            "need_look": s.need_look,
            "fix_input": s.tier_rows[Tier.FIX_INPUT],
            "check_call": s.tier_rows[Tier.CHECK_CALL],
            "informational": s.tier_rows[Tier.INFORMATIONAL],
            "need_look_pct": round(100 * s.need_look / s.rows, 2) if s.rows else 0.0,
            **{f.name: s.flag_rows[f.name] for f in fired},
        }
    return {
        "id": MULTIQC_TABLE_ID,
        "section_name": "gbcms QC summary",
        "section_href": QC_FLAGS_URL,
        "description": (
            "Rows of each gbcms output by what they need: to fix in the input, to check, or "
            "nothing. Hover a column for what the flag means."
        ),
        "helptext": _HELP,
        "plot_type": "table",
        "pconfig": {
            "id": "gbcms_qc_summary_table",
            "title": "gbcms: QC summary",
            "namespace": "gbcms",
            "col1_header": "Sample",
        },
        "headers": headers,
        "data": data,
    }


def multiqc_warning(
    summaries: Sequence[QcSummary], skipped: Sequence[Skipped] = ()
) -> dict[str, Any] | None:
    """The section that says rows need a look or a sample was not checked, or None when
    neither: a report shows it only when there is something to read. MultiQC puts a
    custom-content section the report's ``custom_content: order`` does not name at the
    top."""
    flagged = [s for s in summaries if s.need_look]
    if not flagged and not skipped:
        return None
    rows = sum(s.rows for s in summaries)
    need = sum(s.need_look for s in summaries)
    headers = _tier_headers(only_need=True)
    data: dict[str, dict[str, float]] = {
        s.sample: {
            "need_look": s.need_look,
            "fix_input": s.tier_rows[Tier.FIX_INPUT],
            "check_call": s.tier_rows[Tier.CHECK_CALL],
        }
        for s in flagged
    }
    description = ""
    if flagged:
        description = (
            f"{need:,} of {rows:,} rows need a look, in {len(flagged)} of {len(summaries)} "
            "samples. Only the samples with such rows are listed."
        )
    if skipped:
        # A sample that was not checked is not a clean one.
        headers["not_checked"] = {
            "title": "Not checked",
            "description": "The output could not be read, or was refused: nothing was checked",
            "cond_formatting_rules": {"pass": [{"eq": 0}], "fail": [{"gt": 0}]},
            "format": "{:,.0f}",
            "scale": False,
        }
        for row in data.values():
            row["not_checked"] = 0
        for item in skipped:
            data[item.sample] = {"not_checked": 1}
        note = _not_checked_note(skipped, len(summaries) + len(skipped))
        description = (description + " " + note + ".").strip()
    return {
        "id": MULTIQC_WARNING_ID,
        "section_name": (
            "WARNING: gbcms rows need a look" if flagged else "WARNING: gbcms outputs not checked"
        ),
        "section_href": QC_FLAGS_URL,
        "description": description,
        "helptext": _HELP,
        "plot_type": "table",
        "pconfig": {
            "id": "gbcms_qc_warning_table",
            "title": "gbcms: rows that need a look" if flagged else "gbcms: outputs not checked",
            "namespace": "gbcms",
            "col1_header": "Sample",
        },
        "headers": headers,
        "data": data,
    }


def multiqc_documents(
    summaries: Sequence[QcSummary], skipped: Sequence[Skipped] = ()
) -> dict[str, dict[str, Any]]:
    """The MultiQC files to write: file name -> content (the warning only if needed)."""
    docs: dict[str, dict[str, Any]] = {}
    if summaries:
        docs[MULTIQC_TABLE_NAME] = multiqc_table(summaries)
    warning = multiqc_warning(summaries, skipped)
    if warning:
        docs[MULTIQC_WARNING_NAME] = warning
    return docs


# ── Writing ───────────────────────────────────────────────────────────────────


def _write_text(path: Path, text: str) -> Path:
    with atomic_output(path, "w", newline="") as fh:
        fh.write(text)
    logger.info("Wrote %s", path)
    return path


def _dump(document: dict[str, Any]) -> str:
    return json.dumps(document, indent=2, ensure_ascii=False) + "\n"


def write_summaries(
    summaries: Sequence[QcSummary],
    inputs: Sequence[Path],
    fmt: SummaryFormat,
    output_dir: Path | None = None,
    skipped: Sequence[Skipped] = (),
) -> list[Path]:
    """Write the files of a file format; return their paths. ``json`` puts a file
    beside each input unless ``output_dir`` is given; ``tsv`` and ``multiqc`` write one
    cohort file to ``output_dir`` (the current directory where none is given). ``inputs`` is
    every file the caller was given, read or not: an output never replaces one. ``skipped``
    reaches the MultiQC warning, so a report does not read as complete."""
    written: list[Path] = []
    taken = {Path(p).resolve() for p in inputs}

    def put(path: Path, text: str) -> None:
        if path.resolve() in taken:
            raise ValueError(f"{path} is an input; refusing to overwrite it")
        written.append(_write_text(path, text))

    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
    if fmt is SummaryFormat.JSON:
        for i, s in enumerate(summaries):
            folder = output_dir or (Path(s.source) if s.source else Path(inputs[i])).parent
            put(folder / (s.sample + JSON_SUFFIX), _dump(summary_to_dict(s)))
    elif fmt is SummaryFormat.TSV:
        put((output_dir or Path.cwd()) / TSV_NAME, render_tsv(summaries))
    elif fmt is SummaryFormat.MULTIQC:
        folder = output_dir or Path.cwd()
        documents = multiqc_documents(summaries, skipped)
        for name, document in documents.items():
            put(folder / name, _dump(document))
        stale = folder / MULTIQC_WARNING_NAME
        if MULTIQC_WARNING_NAME not in documents and stale.exists():
            # A warning left by an earlier run would still say rows need a look.
            stale.unlink()
            logger.info("Removed %s: no row needs a look now", stale)
    return written


def text_report(summaries: Sequence[QcSummary], skipped: Sequence[Skipped] = ()) -> str:
    """The ``text`` format: the run's block, then its closing line."""
    return render_summary(summaries, skipped=skipped) + "\n" + closing_line(summaries, skipped)[1]
