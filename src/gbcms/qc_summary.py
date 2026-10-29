"""End-of-run QC summary: which QC flags fired, on how many rows, which rows need
a look, and what to do about each.

The summary reads the output a run has already written (MAF or VCF) and tallies
the flags defined on the QC Flags page (``docs/reference/qc-flags.md``). It
changes no count, column or row, and writes no file: its product is one log block
and a closing line, a WARNING when any row needs a look and an INFO line when none
does. Reading the written output (not the run's in-memory results) keeps it
usable on any gbcms output, and a failure to read it is logged and never fails
the run.

Every flag has a tier:

* fix the input: the row was rejected, or its input is wrong;
* check the call: the row was counted, but its count or VAF may not mean what it
  looks like, so look at it;
* informational: a caveat or a description; nothing to do.

The first two tiers are the rows that "need a look".

``QC_FLAGS`` is the one table holding each flag's tier and one-line action. The QC
Flags page repeats it under "End-of-run summary" so a reader can see it, and
``tests/test_qc_summary.py`` fails when the page and the table differ, when a flag
the page defines has no entry here, and when an entry names no flag on the page.
Each flag's full definition stays on that page.
"""

import csv
import enum
import logging
import re
import textwrap
from collections import Counter
from collections.abc import Callable, Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

from .rescue_audit import rescue_outcome

logger = logging.getLogger(__name__)

# The classes and functions, as in io.output; the constants (``QC_FLAGS`` and the rest) are
# imported by name. A docs test reads every upper-case string literal in the source as a QC
# flag the code emits, so a constant's name listed here would have to be documented as one.
__all__ = [
    "Kind",
    "QcFlag",
    "QcSummary",
    "Tier",
    "log_qc_summary",
    "render_summary",
    "summarize_output",
]

# Where every flag is defined (the mkdocs site_url + the page); a test pins it.
QC_FLAGS_URL = "https://msk-access.github.io/gbcms/reference/qc-flags/"

# Rows named under "Rows to look at", and samples named in a cohort block.
MAX_EXAMPLES = 5
MAX_SAMPLES = 10

# Width an action is wrapped to, so a block reads the same on any console.
_WRAP = 88


class Tier(enum.IntEnum):
    """How a flagged row is read; the first two tiers are the rows that need a look."""

    FIX_INPUT = 1
    CHECK_CALL = 2
    INFORMATIONAL = 3

    @property
    def title(self) -> str:
        return _TIER_TITLES[self]


# The tier headings of the QC Flags page's "End-of-run summary" section.
_TIER_TITLES = {
    Tier.FIX_INPUT: "Fix the input",
    Tier.CHECK_CALL: "Check the call",
    Tier.INFORMATIONAL: "Informational",
}


class Kind(enum.Enum):
    """Where a flag is written; the value is its family's anchor on the QC Flags page."""

    REASON = "verdict-and-status-reasons"
    DIAGNOSTIC = "diagnostics"
    RESCUE = "mnp-rescue"
    ASJD = "asjd-splice-flags"
    COLUMN = "qc-columns"


@dataclass(frozen=True)
class QcFlag:
    """One flag the summary knows: its name as the QC Flags page spells it (a
    parametric flag by its base name), where it is written, its tier, and the
    one-line action the summary prints. No ``|`` and no square brackets in an
    action: it is a table cell on the page and a line in the log."""

    name: str
    kind: Kind
    tier: Tier
    action: str


# The flag table: tier | family | flag | the action the summary prints. A blank line
# separates the tiers. Written as text so each flag is one line, as on the QC Flags
# page; ``_parse_table`` checks every line. An action has no ``|`` and no square brackets.
_TABLE = """
fix | reason | REF_MISMATCH | The REF does not match the reference: check the build and coordinates; REF_AT_OFFSET says where it does match
fix | diagnostic | REF_AT_OFFSET | The REF matches the reference a few bases away: correct the input's coordinates (the row is not moved or counted)
fix | reason | FETCH_FAILED | The reference region could not be read: check contig naming and the FASTA
fix | reason | EMPTY_ALLELE | REF or ALT is empty: fix the input row
fix | reason | NON_SEQUENCE_ALLELE | An allele is not a base sequence: fix the input row
fix | reason | ALT_EQUALS_REF | ALT equals REF, so there is no change to count: fix the input row
fix | reason | START_BEFORE_CONTIG | The row starts before the contig's first base (Start_Position 0): fix the input row; MAF positions are 1-based
fix | reason | ALT_CONTAINS_N | The ALT contains an N base: fix the input row

check | reason | WARN_REF_CORRECTED | The input's REF was corrected to the reference: check the REF; the counts are for the corrected REF
check | reason | WARN_HOMOPOLYMER_DECOMP | The row reports a corrected allele's counts, not the given one: check which you want
check | diagnostic | OBSERVED_ALLELE | No read carries the given allele but others carry another one: the input is likely mis-described; the counts stay the given allele's
check | diagnostic | PARTIAL_DOMINANT | More reads partly match the ALT than carry it, often a coexisting allele the annotation misses: inspect the locus (--trace shows per-read lengths)
check | diagnostic | CLIP_CANDIDATES | An insertion with no ALT reads, but soft-clipped reads may hold its carriers: inspect in IGV
check | diagnostic | POSSIBLE_MINORITY_HAPLOTYPE | A change on a minority of nearby reads may be counted as ALT: weigh a low alt_count with that in mind (the counts are unchanged)
check | diagnostic | MNV_SHARED | Under the sequence judge the SNV's ALT reads include a co-listed MNV's carriers, so one molecule counts in both rows: read the two rows together (the counts are unchanged)
check | diagnostic | UNEXPLAINED_READS | Under the sequence judge some depth reads fit none of the row's candidate sequences: where they are a large share, foreign reads may count (the counts are unchanged)
check | diagnostic | SW_FALLBACK | Some reads were scored by a fallback aligner or left uncounted: the row's counts came partly from another scorer
check | diagnostic | SPLICE_SKIP_DOMINANT | More reads splice over the deleted span than confirm ALT: an alt_count of 0 may hide carriers written as junction reads
check | diagnostic | NON_DISCRIMINATING_LOCUS | A nearby germline combination makes REF and ALT the same sequence: a zero REF and ALT here is not an absence
check | diagnostic | RESCUED_COMPONENT | With --rescue-mnp the row reports a component SNV's counts, not the MNP's: read gbcms_rescue
check | rescue | rescued | The best component SNV beat the MNP, so its counts are written: read gbcms_rescue
check | rescue | ref_validation_failed | No component SNV survived preparation (an anomaly): the MNP's own counts stand; read the logged warning
check | asjd | NON_CANONICAL_MOTIF | The ALT junction's splice motif is not canonical: likely a mapping artifact
check | asjd | STRAND_DISCORDANT | The ALT junction has mixed transcript-strand support: likely an alignment artifact
check | asjd | RETENTION_DOMINANT | Most reads splice over the locus, so the counted reads are the intron-retaining minority: read vaf as the VAF within them, not as allelic balance
check | asjd | NOVEL_JUNC_AT_SPLICE_LOSS | Carriers of a splice-destroying variant may splice around the locus while alt_count is low: read asjd_diagnostic for the junction
check | column | rna_editing_site | The locus overlaps a known A-to-I editing site: an A>G or T>C call there is likely editing, not a somatic mutation

info | reason | MULTI_ALLELIC | another row has a different ALT here
info | reason | TRACT_CLUSTER | shares a repeat window with another indel
info | diagnostic | ZERO_ALT | no read carries the ALT
info | diagnostic | COEXISTING_ALLELE | another allele is more frequent: a VAF caveat
info | diagnostic | MNP_DISC_RATIO | the MNP's shape; nothing to do
info | diagnostic | MNP_RESCUE_ELIGIBLE | a rescue candidate (--rescue-mnp)
info | diagnostic | HIGH_N_FRACTION | over 5% of depth is N: counts are thinner
info | asjd | LOW_ALT_JUNC | fewer than 5 ALT junction fragments
info | asjd | LOW_REF_JUNC | fewer than 10 REF junction fragments
info | asjd | NOVEL_ALT_JUNC | the ALT uses an unannotated junction
info | asjd | MULTI_JUNCTION | ALT fragments use over two junctions
info | column | asjd_flag | ALT and REF junction usage differ
info | rescue | skipped_grouped | in a co-annotated group; rescue left it
info | rescue | haplotype_confirmed | a read shows the whole MNP haplotype
info | rescue | no_improvement | no component beat the MNP; its counts stand
info | column | mfsd_ks_valid | the KS test did not run
info | column | mfsd_alt_confidence | few or no ALT fragments in the size window
info | column | mfsd_ch_flag | in a clonal-hematopoiesis gene: a prior
"""

_TIERS = {"fix": Tier.FIX_INPUT, "check": Tier.CHECK_CALL, "info": Tier.INFORMATIONAL}
_KINDS = {
    "reason": Kind.REASON,
    "diagnostic": Kind.DIAGNOSTIC,
    "rescue": Kind.RESCUE,
    "asjd": Kind.ASJD,
    "column": Kind.COLUMN,
}


def _parse_table(text: str) -> tuple[QcFlag, ...]:
    flags = []
    for line in text.splitlines():
        if not line.strip():
            continue
        tier, kind, name, action = (cell.strip() for cell in line.split(" | ", 3))
        flags.append(QcFlag(name, _KINDS[kind], _TIERS[tier], action))
    return tuple(flags)


QC_FLAGS: tuple[QcFlag, ...] = _parse_table(_TABLE)

_BY_NAME: dict[str, QcFlag] = {f.name: f for f in QC_FLAGS}
if len(_BY_NAME) != len(QC_FLAGS):  # pragma: no cover - a table typo, caught on import
    raise RuntimeError("the QC flag table names a flag twice")

# Where each family is written: (MAF column, VCF INFO key). A column flag is its own
# MAF column; the VCF carries only some of them.
_FAMILY_FIELDS = {
    Kind.REASON: ("gbcms_status_reason", "GSR"),
    Kind.DIAGNOSTIC: ("gbcms_diagnostic", "GD"),
    Kind.RESCUE: ("gbcms_rescue", "GR"),
    Kind.ASJD: ("asjd_diagnostic", "ASJDD"),
}
_COLUMN_VCF = {"rna_editing_site": "RED", "asjd_flag": "ASJD"}

_STATUS_FAIL = "FAIL"


def _shown(flag: QcFlag) -> str:
    """How a flag is named in the log: a rescue outcome as the audit column writes it
    (``outcome=rescued``), every other flag by its name."""
    return f"outcome={flag.name}" if flag.kind is Kind.RESCUE else flag.name


def _field_of(flag: QcFlag, fmt: str) -> str:
    """The MAF column (``fmt`` ``maf``) or VCF INFO key a flag is written in."""
    if flag.kind is Kind.COLUMN:
        return flag.name if fmt == "maf" else _COLUMN_VCF.get(flag.name, flag.name)
    return _FAMILY_FIELDS[flag.kind][0 if fmt == "maf" else 1]


# ── Reading rows ──────────────────────────────────────────────────────────────


class _Row(NamedTuple):
    """The fields of one output row the summary reads. ``columns`` holds the QC
    columns that are set; ``label`` builds the row's ``chrom:pos ref>alt`` only
    when the row is shown."""

    status: str
    reasons: str
    diagnostic: str
    asjd_diagnostic: str
    rescue: str
    columns: tuple[str, ...]
    label: Callable[[], str]


def _maf_rows(path: Path) -> Iterator[_Row]:
    with open(path, newline="") as fh:
        reader = csv.reader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
        header = next(reader, None)
        if header is None or "gbcms_status" not in header:
            raise ValueError(f"{path.name} has no gbcms_status column: not a gbcms MAF")
        at = {name: i for i, name in enumerate(header)}

        def cell(row: list[str], name: str) -> str:
            i = at.get(name)
            return row[i] if i is not None and i < len(row) else ""

        for row in reader:
            columns = tuple(
                f.name
                for f in _COLUMN_FLAGS
                if f.name in at and _column_set(f.name, cell(row, f.name))
            )

            def label(row: list[str] = row) -> str:
                ref, alt = cell(row, "Reference_Allele"), cell(row, "Tumor_Seq_Allele2")
                if alt == ref:
                    alt = cell(row, "Tumor_Seq_Allele1") or alt
                gene = cell(row, "Hugo_Symbol")
                where = f"{cell(row, 'Chromosome')}:{cell(row, 'Start_Position')} {ref}>{alt}"
                return f"{where} ({gene})" if gene and gene != "Unknown" else where

            yield _Row(
                cell(row, "gbcms_status"),
                cell(row, "gbcms_status_reason"),
                cell(row, "gbcms_diagnostic"),
                cell(row, "asjd_diagnostic"),
                cell(row, "gbcms_rescue"),
                columns,
                label,
            )


def _vcf_rows(path: Path) -> Iterator[_Row]:
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            info = {}
            for part in f[7].split(";"):
                key, _, value = part.partition("=")
                info[key] = value
            columns = tuple(name for name, key in _COLUMN_VCF.items() if key in info)

            def label(f: list[str] = f) -> str:
                return f"{f[0]}:{f[1]} {f[3]}>{f[4]}"

            yield _Row(
                info.get("GS", ""),
                info.get("GSR", ""),
                info.get("GD", ""),
                info.get("ASJDD", ""),
                info.get("GR", ""),
                columns,
                label,
            )


_COLUMN_FLAGS = tuple(f for f in QC_FLAGS if f.kind is Kind.COLUMN)


def _column_set(name: str, value: str) -> bool:
    """Whether a QC column's value is the one the flag names (a boolean column set
    to ``True``, or the KS test not run)."""
    if name == "mfsd_ks_valid":
        return value == "False"
    if name == "mfsd_alt_confidence":
        return value in ("SPARSE", "NONE")
    return value == "True"


# ── Tallying ──────────────────────────────────────────────────────────────────

_SEPARATORS = re.compile(r"[;|]")


def _tokens(text: str) -> list[str]:
    """The flags of a ``;``- or ``|``-separated value (``.`` is the VCF's none)."""
    return [t for t in _SEPARATORS.split(text) if t and t != "."]


def _base(token: str) -> str:
    """``COEXISTING_ALLELE`` from ``COEXISTING_ALLELE(chr1:100:A>T:5/2)``."""
    return token.split("(", 1)[0]


@dataclass
class QcSummary:
    """The tally of one sample's output. ``flag_rows`` counts rows per flag;
    ``tier_rows`` counts rows carrying at least one flag of a tier; ``need_look``
    counts rows in the first two tiers (once each); ``unlisted`` counts flags the
    table does not know (another build's, or a newer one's)."""

    sample: str
    fmt: str = "maf"
    rows: int = 0
    need_look: int = 0
    info_only: int = 0
    unknown_only: int = 0
    flag_rows: Counter[str] = field(default_factory=Counter)
    tier_rows: Counter[Tier] = field(default_factory=Counter)
    unlisted: Counter[str] = field(default_factory=Counter)
    # (label, flags) of rows that need a look: input problems first, then the rest.
    examples: list[tuple[str, str]] = field(default_factory=list)


def _tally(
    rows: Iterable[_Row], sample: str, fmt: str = "maf", max_examples: int = MAX_EXAMPLES
) -> QcSummary:
    """Tally rows into a :class:`QcSummary`."""
    summary = QcSummary(sample=sample, fmt=fmt)
    fix_examples: list[tuple[str, str]] = []
    check_examples: list[tuple[str, str]] = []
    for row in rows:
        summary.rows += 1
        hits: dict[str, tuple[QcFlag, str]] = {}
        unlisted: dict[str, str] = {}
        for text in (row.reasons, row.diagnostic, row.asjd_diagnostic):
            for token in _tokens(text):
                flag = _BY_NAME.get(_base(token))
                if flag is None:
                    unlisted.setdefault(_base(token), token)
                else:
                    hits.setdefault(flag.name, (flag, token))
        outcome = rescue_outcome(row.rescue) if row.rescue else None
        if outcome in _BY_NAME:
            hits.setdefault(outcome, (_BY_NAME[outcome], _shown(_BY_NAME[outcome])))
        elif outcome:
            unlisted.setdefault(outcome, f"outcome={outcome}")
        for name in row.columns:
            hits.setdefault(name, (_BY_NAME[name], name))

        tiers = {flag.tier for flag, _ in hits.values()}
        failed = row.status == _STATUS_FAIL
        if failed:
            tiers.add(Tier.FIX_INPUT)
        for name in hits:
            summary.flag_rows[name] += 1
        for name in unlisted:
            summary.unlisted[name] += 1
        for tier in tiers:
            summary.tier_rows[tier] += 1

        if Tier.FIX_INPUT in tiers or Tier.CHECK_CALL in tiers:
            summary.need_look += 1
            shown = [t for flag, t in hits.values() if flag.tier is not Tier.INFORMATIONAL]
            if failed:
                shown += list(unlisted.values())
            target = fix_examples if Tier.FIX_INPUT in tiers else check_examples
            if len(target) < max_examples:
                target.append((row.label(), ", ".join(shown)))
        elif tiers:
            summary.info_only += 1
        elif unlisted:
            summary.unknown_only += 1
    summary.examples = (fix_examples + check_examples)[:max_examples]
    return summary


def summarize_output(
    path: Path | str, sample: str, fmt: str | None = None, max_examples: int = MAX_EXAMPLES
) -> QcSummary:
    """Tally a gbcms output file (``fmt`` is ``maf`` or ``vcf``; the extension by
    default). Read-only: the file is opened for reading and nothing is written."""
    path = Path(path)
    fmt = fmt or ("vcf" if path.suffix == ".vcf" else "maf")
    rows = _vcf_rows(path) if fmt == "vcf" else _maf_rows(path)
    return _tally(rows, sample, fmt, max_examples)


# ── Rendering ─────────────────────────────────────────────────────────────────


def _count(n: int, noun: str = "row") -> str:
    return f"{n:,} {noun}{'' if n == 1 else 's'}"


def _share(part: int, whole: int) -> str:
    return f"{100 * part / whole:.1f}%" if whole else "0%"


def _fill(text: str, first: int, rest: int) -> list[str]:
    """``text`` wrapped to the block's width, so it reads the same on any console."""
    return textwrap.wrap(
        text,
        _WRAP,
        initial_indent=" " * first,
        subsequent_indent=" " * rest,
        break_on_hyphens=False,
    )


def _tier_lines(tier: Tier, summaries: Sequence[QcSummary]) -> list[str]:
    """The lines of one tier: its row count, then each flag that fired, most rows
    first. A flag of the first two tiers carries its action; an informational one a
    short gloss on its own line."""
    many = len(summaries) > 1
    rows = sum(s.tier_rows[tier] for s in summaries)
    fired = [(f, sum(s.flag_rows[f.name] for s in summaries)) for f in QC_FLAGS if f.tier is tier]
    fired = sorted(((f, n) for f, n in fired if n), key=lambda fn: -fn[1])
    if not rows and not fired:
        return [f"  {tier.title}: none"]
    lines = [f"  {tier.title}: {_count(rows)}"]
    for f, n in fired:
        where = ""
        if many:
            k = sum(1 for s in summaries if s.flag_rows[f.name])
            where = f" in {k} of {len(summaries)} samples"
        if tier is Tier.INFORMATIONAL:
            lines += _fill(f"{_shown(f)}: {_count(n)}{where} ({f.action})", 4, 8)
        else:
            lines.append(f"    {_shown(f)}: {_count(n)}{where}")
            lines += _fill(f.action, 6, 6)
    return lines


def render_summary(summaries: Sequence[QcSummary], max_examples: int = MAX_EXAMPLES) -> str:
    """The summary block for one sample, or the cohort block for several."""
    many = len(summaries) > 1
    rows = sum(s.rows for s in summaries)
    need = sum(s.need_look for s in summaries)
    info_only = sum(s.info_only for s in summaries)
    unknown_only = sum(s.unknown_only for s in summaries)
    fmt = summaries[0].fmt

    who = f"{len(summaries)} samples" if many else f"sample {summaries[0].sample}"
    lines = [
        f"QC summary for {who}: {need:,} of {_count(rows)} {_need(need)} a look "
        f"({_share(need, rows)})",
        f"  Rows with only informational flags: {info_only:,}; with no flag: "
        f"{rows - need - info_only - unknown_only:,}"
        + (f"; with only flags this build does not know: {unknown_only:,}" if unknown_only else ""),
    ]
    for tier in Tier:
        lines += _tier_lines(tier, summaries)

    unlisted: Counter[str] = Counter()
    for s in summaries:
        unlisted.update(s.unlisted)
    if unlisted:
        named = ", ".join(f"{name} ({_count(n)})" for name, n in unlisted.most_common())
        lines += _fill("Not in this build's flag table (written by another build): " + named, 2, 4)

    if many and need:
        ranked = sorted((s for s in summaries if s.need_look), key=lambda s: -s.need_look)
        lines.append(
            f"  Samples with rows to look at (first {min(len(ranked), MAX_SAMPLES)} "
            f"of {len(ranked)}):"
        )
        for s in ranked[:MAX_SAMPLES]:
            lines.append(f"    {s.sample}: {s.need_look:,} of {_count(s.rows)}")

    shown = [(s.sample, label, flags) for s in summaries for label, flags in s.examples]
    if shown:
        order = "" if many else ", input problems first"
        lines.append(
            f"  Rows to look at (first {min(len(shown), max_examples)} of {need:,}{order}):"
        )
        for sample, label, shown_flags in shown[:max_examples]:
            prefix = f"{sample} " if many else ""
            lines.append(f"    {prefix}{label}: {shown_flags}")

    fired = {
        f
        for f in QC_FLAGS
        if f.tier is not Tier.INFORMATIONAL and any(s.flag_rows[f.name] for s in summaries)
    }
    fields = {_field_of(f, fmt) for f in fired}
    if any(s.tier_rows[Tier.FIX_INPUT] for s in summaries):
        fields.add("gbcms_status" if fmt == "maf" else "GS")
    if fields:
        kind = "MAF columns" if fmt == "maf" else "VCF INFO keys"
        lines.append(f"  Find them in the output ({kind}): {', '.join(sorted(fields))}")
    lines.append("  What each flag means and what to do: " + QC_FLAGS_URL)
    return "\n".join(lines)


def _need(n: int) -> str:
    return "needs" if n == 1 else "need"


def log_qc_summary(summaries: Sequence[QcSummary]) -> None:
    """Log the summary block (INFO), then, as the run's last line, a WARNING when
    any row needs a look and an INFO line when none does."""
    if not summaries:
        return
    # markup off: sample names and labels are user text, not Rich markup.
    logger.info("%s", render_summary(summaries), extra={"markup": False})
    rows = sum(s.rows for s in summaries)
    need = sum(s.need_look for s in summaries)
    if not need:
        logger.info("No row needs a look: %s checked (QC flags: %s)", _count(rows), QC_FLAGS_URL)
        return
    fix = sum(s.tier_rows[Tier.FIX_INPUT] for s in summaries)
    where = ""
    if len(summaries) > 1:
        k = sum(1 for s in summaries if s.need_look)
        where = f" in {k} of {len(summaries)} samples"
    logger.warning(
        "%s of %s %s a look%s (%s with a problem in the input): see the QC summary above " "and %s",
        f"{need:,}",
        _count(rows),
        _need(need),
        where,
        f"{fix:,}",
        QC_FLAGS_URL,
        extra={"markup": False},
    )
