"""
Multi-BAM MAF merge engine using Polars.

Merges per-BAM-type genotyped MAFs (e.g., duplex, simplex) into a single
output MAF with type-prefixed count columns.  Uses Polars lazy API for
WGS-scale performance.

Architecture:
    1. Scan each input MAF lazily via ``io.batch.scan_maf``
    2. Detect whether columns are already prefixed or need renaming
    3. Progressive outer join on the variant — contig, Start and alleles
       (plus the VCF record — vcf_pos / vcf_ref / vcf_alt — when every input
       carries it); End_Position is filled from the inputs that have each row
    4. Optionally compute additive ``simplex_duplex_*`` combined columns
    5. Materialize and write via ``io.batch.write_maf``

Design decisions:
    - All columns read as strings (MAF is text). Cast to numeric only
      for combined column arithmetic.
    - ``fillna("0")`` for count columns after outer join — missing variants
      get zero counts, not NULL (matches genotype_variants convention).
    - Logging at INFO for every operation (timing, row/col counts).
    - Every error includes file path and column context.
"""

import logging
import re
import time
from pathlib import Path

import polars as pl

from gbcms.core.kernel import CoordinateKernel
from gbcms.io.batch import scan_maf, write_maf
from gbcms.io.output import (
    _fmt,
    _fmt_sci,
    gbcms_column_basenames,
    gbcms_prefixed_basenames,
    provenance_line,
)
from gbcms.models.core import MergeConfig
from gbcms.rescue_audit import rescued_component

__all__ = ["merge_mafs"]

logger = logging.getLogger(__name__)


# ── Constants (single source of truth — canonical basenames from output.py) ──

# The MAF variant columns. Every input needs all but End_Position, which
# follows from Start and REF and is optional in a MAF.
VARIANT_KEY: list[str] = [
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
]

# Naming-independent contig key the joins use in place of Chromosome, so MAFs
# whose inputs name contigs differently (chr1 vs 1, chrM vs MT) still join.
_CONTIG_KEY = "_contig_key"
# A variant is its contig, Start and alleles. End_Position is not joined on:
# inputs that write it differently for one variant still join, and each row's
# End_Position is filled from the inputs that have the row.
JOIN_KEY: list[str] = [
    _CONTIG_KEY,
    *(k for k in VARIANT_KEY[1:] if k != "End_Position"),
]
# VCF-input MAFs also carry the VCF record each row came from. It is unique per
# record ALT, whereas two records can trim to one MAF record (TCT>TCG and T>G at
# the changed base), so the joins add it whenever every input has it. The MAF
# key stays in the join key: within one gbcms version it follows from the
# record (the pairing is the record's), and a full join fills only key
# columns, so a row only a later input has keeps its coordinates and alleles.
VCF_RECORD_KEY: list[str] = [*JOIN_KEY, "vcf_pos", "vcf_ref", "vcf_alt"]
# Prefixes of the other per-input join helpers (row numbers, each later
# input's own contig names and End_Position, a later-only row's annotations).
# Input columns with these names are rejected.
_HELPER_PREFIXES = ("_row_", "_chrom_", "_end_", "_ann_")


def _row_col(bam_type: str) -> str:
    """Per-input row-number column carried through the joins, so merged rows
    can be put back in the inputs' order (a full join guarantees none)."""
    return f"_row_{bam_type}"


def _reject_helper_columns(columns: list[str], bam_type: str, path: Path) -> None:
    """Refuse input columns that collide with merge's join helpers: they would
    be overwritten or shadow the helper and corrupt the join."""
    clash = [c for c in columns if c == _CONTIG_KEY or c.startswith(_HELPER_PREFIXES)]
    if clash:
        raise ValueError(
            f"{bam_type} MAF ({path}) has column(s) {clash} that gbcms merge reserves "
            "for its join helpers; rename or drop them"
        )


def _with_contig_key(lf: pl.LazyFrame, bam_type: str) -> pl.LazyFrame:
    """Add the join's contig key: :meth:`CoordinateKernel.contig_key` (the
    counting engine's rule), computed once per distinct contig name.

    Warns when this input names one contig more than one way (e.g. ``chrM``
    and ``MT``): the same variant under both names joins twice, so the merged
    output carries a duplicate row for it.
    """
    names = lf.select(pl.col("Chromosome").unique()).collect().to_series().drop_nulls()
    keys = {n: CoordinateKernel.contig_key(n) for n in names.to_list()}
    spellings: dict[str, list[str]] = {}
    for name, key in keys.items():
        spellings.setdefault(key, []).append(name)
    for aliases in spellings.values():
        if len(aliases) > 1:
            logger.warning(
                "  '%s' names one contig %d ways (%s): a variant listed under more than "
                "one of them joins once per name, so the merged output repeats it",
                bam_type,
                len(aliases),
                ", ".join(sorted(aliases)),
            )
    return lf.with_columns(
        pl.col("Chromosome")
        .replace_strict(keys, default=None, return_dtype=pl.String)
        .alias(_CONTIG_KEY)
    )


# gbcms count column basenames (without any prefix): the counts a row an input
# lacks gets as 0 for that input (see GBCMS_META_BASENAMES for its empty ones).
GBCMS_COUNT_BASENAMES: list[str] = [
    "ref_count",
    "alt_count",
    "any_alt",
    "partial_alt",
    "n_count",
    "total_count",
    "vaf",
    "ref_count_forward",
    "ref_count_reverse",
    "alt_count_forward",
    "alt_count_reverse",
    "ref_count_fragment",
    "alt_count_fragment",
    "total_count_fragment",
    "vaf_fragment",
    "ref_count_fragment_forward",
    "ref_count_fragment_reverse",
    "alt_count_fragment_forward",
    "alt_count_fragment_reverse",
]

# Non-count gbcms columns that also get type-prefixed.
# These are string/diagnostic columns (not numeric).
GBCMS_META_BASENAMES: list[str] = [
    "gbcms_status",
    "gbcms_status_reason",
    "gbcms_diagnostic",
    "gbcms_rescue",
    "strand_bias_p_value",
    "strand_bias_odds_ratio",
    "fragment_strand_bias_p_value",
    "fragment_strand_bias_odds_ratio",
]

# Every column gbcms writes, in any mode (mFSD, RNA, GTF, rescue,
# normalization), unprefixed: each is kept per input (``duplex_mfsd_ref_mean``).
# Taken from the writer so a new output column cannot fall back to "the first
# input's value" unnoticed.
ALL_GBCMS_BASENAMES: set[str] = set(gbcms_column_basenames())

# Additive count basenames for simplex+duplex combination.
# Duplex and simplex BAMs contain distinct consensus molecules — there is
# no double-counting — so ALL count levels are additive across BAM types.
#
# Column order follows output.py sections:
#   §2 read-level → §3 read strand → §4 fragment-level → §5 fragment strand
COMBINED_ADDITIVE_READ: list[str] = [
    "ref_count",
    "alt_count",
]
COMBINED_ADDITIVE_READ_STRAND: list[str] = [
    "ref_count_forward",
    "ref_count_reverse",
    "alt_count_forward",
    "alt_count_reverse",
]
COMBINED_ADDITIVE_FRAGMENT: list[str] = [
    "ref_count_fragment",
    "alt_count_fragment",
]
COMBINED_ADDITIVE_FRAGMENT_STRAND: list[str] = [
    "ref_count_fragment_forward",
    "ref_count_fragment_reverse",
    "alt_count_fragment_forward",
    "alt_count_fragment_reverse",
]
# Flat list of all additive basenames (for iteration).
COMBINED_ADDITIVE_ALL: list[str] = (
    COMBINED_ADDITIVE_READ
    + COMBINED_ADDITIVE_READ_STRAND
    + COMBINED_ADDITIVE_FRAGMENT
    + COMBINED_ADDITIVE_FRAGMENT_STRAND
)


def merge_mafs(config: MergeConfig) -> None:
    """Merge per-BAM-type genotyped MAFs into a single type-prefixed output.

    Performs an outer join on the variant (contig, Start and alleles; plus the
    VCF record — vcf_pos / vcf_ref / vcf_alt — when every input is
    VCF-derived), prefixes gbcms count columns with the BAM type label,
    optionally computes combined simplex_duplex columns, and writes the merged
    result.

    Args:
        config: Validated MergeConfig with inputs, output path, and options.

    Raises:
        FileNotFoundError: If any input MAF does not exist.
        ValueError: If column detection fails or variant key columns are missing.
    """
    t_start = time.perf_counter()
    logger.info(
        "Starting merge: %d inputs → %s",
        len(config.inputs),
        config.output,
    )

    # ── 1. Scan and rename ────────────────────────────────────────────────────
    frames: dict[str, pl.LazyFrame] = {}
    input_columns: dict[str, list[str]] = {}
    versions: dict[str, str | None] = {}
    for bam_type, path in config.inputs.items():
        versions[bam_type] = _version_line(path)
        logger.info("Scanning %s MAF: %s", bam_type, path)
        lf = scan_maf(path)

        # Validate variant key columns exist
        schema_names = lf.collect_schema().names()
        _validate_variant_key(schema_names, bam_type, path)
        _reject_helper_columns(schema_names, bam_type, path)

        # Detect and rename gbcms columns with type prefix
        rename_map = _build_rename_map(schema_names, bam_type)
        if rename_map:
            lf = lf.rename(rename_map)
            logger.info(
                "  Prefixed %d columns with '%s_'",
                len(rename_map),
                bam_type,
            )
        else:
            logger.info("  Columns already prefixed for '%s', using as-is", bam_type)

        frames[bam_type] = _with_contig_key(lf, bam_type).with_row_index(_row_col(bam_type))
        input_columns[bam_type] = schema_names

    _refuse_mixed_vcf_representations(input_columns, versions)
    _warn_mixed_versions(versions)

    # ── 2. Progressive outer join ─────────────────────────────────────────────
    join_key = _join_key(frames, input_columns)
    for bam_type, lf in frames.items():
        _warn_duplicate_keys(lf, join_key, bam_type)
    types = list(frames.keys())
    merged = frames[types[0]]

    for join_type in types[1:]:
        # Select only variant key + gbcms columns from the joining frame
        # to avoid duplicating annotation columns across inputs.
        # The joining frame's own Chromosome and End_Position are kept aside
        # (not join keys) so rows only it has keep them, and differences can
        # be reported after materialization.
        join_frame = frames[join_type]
        join_names = join_frame.collect_schema().names()
        count_cols = [c for c in join_names if _is_prefixed_gbcms_col(c, join_type)]
        aside = {"Chromosome": f"_chrom_{join_type}"}
        if "End_Position" in join_names and "End_Position" not in join_key:
            aside["End_Position"] = f"_end_{join_type}"
        merged = merged.join(
            join_frame.select([*join_key, *aside, _row_col(join_type), *count_cols]).rename(aside),
            on=join_key,
            how="full",
            coalesce=True,
        )
        logger.info("  Joined '%s' (%d count cols)", join_type, len(count_cols))

    # ── 3. A row an input lacks: its counts 0, its status columns empty ─────
    # A row the input has keeps its cells as written: a missing count stays
    # missing (it is not a zero), and the combined columns say NA for it.
    output_schema = set(merged.collect_schema().names())
    fills = []
    for t in types:
        absent = pl.col(_row_col(t)).is_null()
        for bases, value in ((GBCMS_COUNT_BASENAMES, "0"), (GBCMS_META_BASENAMES, "")):
            for base in bases:
                col = f"{t}_{base}"
                if col in output_schema:
                    fills.append(
                        pl.when(absent).then(pl.lit(value)).otherwise(pl.col(col)).alias(col)
                    )
    if fills:
        merged = merged.with_columns(fills)

    # ── 4. Combined simplex+duplex columns ────────────────────────────────────
    if config.add_combined and "simplex" in frames and "duplex" in frames:
        merged = _add_combined_columns(merged)
        logger.info("  Added simplex_duplex combined columns")

    # ── 5. Materialize and write ──────────────────────────────────────────────
    collected = _fill_later_only_rows(merged.collect(), frames, types, join_key)
    for t in types:
        lacking = collected.filter(pl.col(_row_col(t)).is_null()).height
        if lacking:
            logger.info(
                "  '%s' lacks %d of %d merged rows: its counts there are written as 0",
                t,
                lacking,
                collected.height,
            )
    result = _in_input_order(
        _resolve_end_positions(_resolve_contig_names(collected, types), types), types
    )
    logger.info(
        "Merged result: %d rows × %d columns",
        result.height,
        result.width,
    )

    if result.height == 0:
        logger.warning(
            "Merged output is empty (0 rows). Check that input MAFs "
            "contain data and share variant key columns."
        )

    if "duplex" in frames and "simplex" in frames:
        result = _mixed_rescue(result, combined=config.add_combined)

    # Legacy naming pass (rename {type}_{metric} → t_{metric}_{type})
    if config.legacy_naming:
        result = _apply_legacy_naming(result, types)
        logger.info("  Applied legacy t_{metric}_{type} naming")

    header = [provenance_line()]
    if config.command_line:
        header.append(f"#command {config.command_line}")
    header += [
        f"#input {t}: {versions[t] or 'no gbcms version line'} ({config.inputs[t]})" for t in types
    ]
    write_maf(result, config.output, header=header)
    elapsed = time.perf_counter() - t_start
    logger.info("Merge complete in %.1fs", elapsed)


# ── Helpers ──────────────────────────────────────────────────────────────────


def _version_line(path: Path) -> str | None:
    """The ``gbcms v...`` provenance an input MAF starts with (``#gbcms v6.5.0``,
    with a build commit since 6.6.0), or None when it has none."""
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                return None
            if line.startswith("#gbcms v"):
                return line[1:].strip()
    return None


_VERSION_LINE = re.compile(r"gbcms v(\S+)(?: \(([0-9a-f]+)\))?")


def _parse_version(line: str | None) -> tuple[str, str | None] | None:
    """(version, commit) from a ``gbcms v6.6.0.dev0 (9c371263)`` line."""
    m = _VERSION_LINE.match(line or "")
    return (m.group(1), m.group(2)) if m else None


def _release(version: str) -> tuple[int, ...]:
    """The numeric release of a version string (``6.6.0.dev0`` → (6, 6, 0))."""
    return tuple(int(x) for x in re.findall(r"\d+", version)[:3])


def _refuse_mixed_vcf_representations(
    input_columns: dict[str, list[str]], versions: dict[str, str | None]
) -> None:
    """Stop when one input is a VCF-input MAF from before 6.5.0 and another is
    not. Before 6.5.0 a VCF-input MAF carried ``vcf_pos`` (with ``vcf_id``,
    ``vcf_region``) but not ``vcf_ref``/``vcf_alt``, and many indel and MNP rows
    had other coordinates and alleles (6.5.0 follows vcf2maf), so their rows do
    not join: merged, the same variant appears twice, each half empty. Only an
    input whose version line is missing or older than 6.5.0 has that shape: a
    later MAF-input output can carry ``vcf_pos`` from vcf2maf."""

    def old_shape(t: str) -> bool:
        cols = set(input_columns[t])
        if "vcf_pos" not in cols or {"vcf_ref", "vcf_alt"} <= cols:
            return False
        parsed = _parse_version(versions.get(t))
        return parsed is None or _release(parsed[0]) < (6, 5, 0)

    old = sorted(t for t in input_columns if old_shape(t))
    if old and len(old) < len(input_columns):
        others = sorted(set(input_columns) - set(old))
        raise ValueError(
            f"Input(s) {old} are VCF-input MAFs from before gbcms 6.5.0 (vcf_pos without "
            f"vcf_ref/vcf_alt), and {others} are not: 6.5.0 changed the coordinates and "
            "alleles of VCF-input rows, so their rows would not join. Genotype every "
            "input with one gbcms version and merge again."
        )


def _warn_mixed_versions(versions: dict[str, str | None]) -> None:
    """Warn when the inputs come from different gbcms versions or builds (their
    counts may follow different rules, which their sums would hide), or when only
    some inputs say which version wrote them. Builds of one version differ only
    when both name a commit."""
    parsed = {t: _parse_version(v) for t, v in versions.items()}
    seen = list(parsed.values())

    def differ(a: tuple[str, str | None] | None, b: tuple[str, str | None] | None) -> bool:
        if a is None or b is None:
            return (a is None) != (b is None)
        return a[0] != b[0] or (a[1] is not None and b[1] is not None and a[1] != b[1])

    if any(differ(a, b) for i, a in enumerate(seen) for b in seen[i + 1 :]):
        logger.warning(
            "Inputs come from different gbcms versions (%s): counts made by different "
            "versions may follow different rules; merge outputs of one version",
            ", ".join(f"{t}: {versions[t] or 'no version line'}" for t in versions),
        )


def _join_key(frames: dict[str, pl.LazyFrame], input_columns: dict[str, list[str]]) -> list[str]:
    """The VCF record key when every input carries it (all VCF-derived), else
    the MAF variant key. End_Position joins too when every input has it and one
    lists a variant (contig, Start, alleles) more than once with different
    End_Position, so each such row pairs with its own counterpart rather than
    with every one of them."""
    key = JOIN_KEY
    if all("End_Position" in cols for cols in input_columns.values()):
        split = [t for t, lf in frames.items() if _end_splits_a_variant(lf)]
        if split:
            logger.info(
                "  Joining on End_Position too: %s list(s) a variant more than once with "
                "different End_Position",
                ", ".join(f"'{t}'" for t in split),
            )
            key = [*JOIN_KEY, "End_Position"]
    if all(set(VCF_RECORD_KEY[len(JOIN_KEY) :]) <= set(cols) for cols in input_columns.values()):
        logger.info("  Joining on the VCF record (vcf_pos, vcf_ref, vcf_alt): every input has it")
        return [*key, *VCF_RECORD_KEY[len(JOIN_KEY) :]]
    return key


def _end_splits_a_variant(lf: pl.LazyFrame) -> bool:
    """Whether an input lists one variant (the join key) with more than one
    End_Position."""
    split = lf.group_by(JOIN_KEY).agg(pl.col("End_Position").n_unique().alias("ends"))
    return bool(split.filter(pl.col("ends") > 1).select(pl.len()).collect().item())


def _warn_duplicate_keys(lf: pl.LazyFrame, join_key: list[str], bam_type: str) -> None:
    """Warn when an input lists one join key more than once: a full join pairs
    each such row with every matching row of the other inputs, so the merged
    output repeats them."""
    dups = lf.group_by(join_key).len().filter(pl.col("len") > 1).collect()
    if dups.height:
        logger.warning(
            "  '%s' has %d duplicate join key(s) (%d rows): each joins with every matching "
            "row of the other inputs, so the merged output repeats them",
            bam_type,
            dups.height,
            int(dups["len"].sum()),
        )


def _fill_later_only_rows(
    result: pl.DataFrame,
    frames: dict[str, pl.LazyFrame],
    types: list[str],
    join_key: list[str],
) -> pl.DataFrame:
    """Fill the first input's annotation columns of rows the first input lacks.

    Non-key columns come from the first input, so a row only a later input has
    would otherwise carry them empty (gene, sample barcode, classification...).
    Each such row takes them from the earliest later input that has the row and
    the column; the column set stays the first input's, and a row the first
    input has keeps its own values. Contig names and End_Position have their
    own rules (:func:`_resolve_contig_names`, :func:`_resolve_end_positions`).
    No work when every row is in the first input, as when every flavor was
    genotyped from one variant file.
    """
    first = types[0]
    lacks_first = pl.col(_row_col(first)).is_null()
    missing = result.filter(lacks_first)
    if missing.height == 0:
        return result
    first_names = frames[first].collect_schema().names()
    ann = [
        c
        for c in first_names
        if c not in join_key
        and c not in ("Chromosome", "End_Position")
        and not c.startswith("_")
        and not _is_prefixed_gbcms_col(c, first)
    ]
    filled = 0
    for t in types[1:]:
        names = set(frames[t].collect_schema().names())
        cols = [c for c in ann if c in names]
        if not cols:
            continue
        # By the later input's own row (not the join key): two rows sharing a key
        # (one variant listed for two samples) each keep their own annotations.
        row = _row_col(t)
        later = frames[t].select([row, *[pl.col(c).alias(f"_ann_{t}_{c}") for c in cols]])
        missing = missing.join(later.collect(), on=row, how="left")
        missing = missing.with_columns(
            [pl.coalesce(c, f"_ann_{t}_{c}").alias(c) for c in cols]
        ).drop([f"_ann_{t}_{c}" for c in cols])
        filled += len(cols)
    if filled:
        logger.info(
            "  %d row(s) only a later input has: their annotation columns are that input's",
            missing.height,
        )
    return pl.concat([result.filter(~lacks_first), missing.select(result.columns)])


def _resolve_contig_names(result: pl.DataFrame, types: list[str]) -> pl.DataFrame:
    """Settle Chromosome after naming-independent joins and drop the helpers.

    Each contig is written one way: the first input's name for it, or — for a
    contig the first input lacks — the name of the earliest input that has it.
    Rows the first input has keep its name as written. One INFO line per later
    input that names contigs differently from the output; its rows still
    joined, on the normalized contig.
    """
    later = [f"_chrom_{t}" for t in types[1:]]
    written: dict[str, str] = {}
    for col in ["Chromosome", *later]:
        pairs = result.select(_CONTIG_KEY, col).drop_nulls().unique(maintain_order=True)
        for key, name in pairs.iter_rows():
            written.setdefault(key, name)
    result = result.with_columns(
        pl.coalesce(
            "Chromosome",
            pl.col(_CONTIG_KEY).replace_strict(written, default=None, return_dtype=pl.String),
        ).alias("Chromosome")
    )
    for t, col in zip(types[1:], later, strict=True):
        differ = result.filter(pl.col(col).is_not_null() & (pl.col(col) != pl.col("Chromosome")))
        if differ.height:
            row = differ.row(0, named=True)
            logger.info(
                "  '%s' and the merged output name contigs differently ('%s' vs '%s', "
                "%d row(s)): joined on the normalized contig; the output names each "
                "contig one way",
                t,
                row[col],
                row["Chromosome"],
                differ.height,
            )
    return result.drop([_CONTIG_KEY, *later])


def _resolve_end_positions(result: pl.DataFrame, types: list[str]) -> pl.DataFrame:
    """Fill End_Position from the inputs that have each row and drop the helpers.

    End_Position is not a join key, so a row only a later input has would
    otherwise lose it. A row the first input has keeps its End_Position as
    written. One INFO line per later input that writes it differently for
    rows it shares with the output. No-op on the column set when the first
    input has no End_Position (the output's columns are the first input's).
    """
    later = [f"_end_{t}" for t in types[1:] if f"_end_{t}" in result.columns]
    if not later:
        return result
    if "End_Position" in result.columns:
        result = result.with_columns(pl.coalesce("End_Position", *later).alias("End_Position"))
        for col in later:
            differ = result.filter(
                pl.col(col).is_not_null() & (pl.col(col) != pl.col("End_Position"))
            )
            if differ.height:
                row = differ.row(0, named=True)
                logger.info(
                    "  '%s' writes End_Position differently from the merged output ('%s' vs "
                    "'%s' at %s:%s, %d row(s)): joined on Start and alleles; the output keeps "
                    "the earliest input's",
                    col.removeprefix("_end_"),
                    row[col],
                    row["End_Position"],
                    row["Chromosome"],
                    row["Start_Position"],
                    differ.height,
                )
    return result.drop(later)


def _in_input_order(result: pl.DataFrame, types: list[str]) -> pl.DataFrame:
    """Merged rows in the inputs' order, then the row numbers dropped.

    The first input's rows as it lists them, then rows only a later input has,
    in that input's order. A full join guarantees no order — without this the
    merged rows came out differently on every run.
    """
    rows = [_row_col(t) for t in types]
    return result.sort(rows, nulls_last=True).drop(rows)


def _mixed_rescue(result: pl.DataFrame, combined: bool) -> pl.DataFrame:
    """Rows whose duplex and simplex MNP rescue outcomes differ: warned, and
    their combined columns NA.

    A rescued row reports a component SNV's counts under the MNP's coordinates.
    When only one flavor was rescued — or the two adopted different components —
    the two flavors' counts describe different alleles in one row, so summing
    them in ``simplex_duplex_*`` would add two alleles: those cells are NA (a
    missing value, as for a missing count), and each row is named in the log.
    The per-flavor columns stay as they are. A row one flavor lacks is not mixed:
    nothing is summed across alleles there. No-op when rescue was run on neither
    flavor (no ``gbcms_rescue`` column). When it was run on only one, that is
    logged once and the other flavor is treated as reporting the MNP on every
    row, so each row rescued in the rescue-on flavor is mixed.
    """
    d, s = "duplex_gbcms_rescue", "simplex_gbcms_rescue"
    present = [col for col in (d, s) if col in result.columns]
    if not present:
        return result
    if len(present) == 1:
        ran, other = ("duplex", "simplex") if present[0] == d else ("simplex", "duplex")
        logger.warning(
            "MNP rescue was run on %s only (%s was genotyped without --rescue-mnp): rows "
            "rescued in %s report a component while %s reports the MNP",
            ran,
            other,
            ran,
            other,
        )
    statuses = [c for c in ("duplex_gbcms_status", "simplex_gbcms_status") if c in result.columns]
    key = [k for k in VARIANT_KEY if k in result.columns]
    mixed: list[int] = []
    for i, row in enumerate(result.select([*key, *present, *statuses]).iter_rows(named=True)):
        # A flavor that lacks the row has an empty status: nothing is mixed.
        if any(not row.get(c) for c in statuses):
            continue
        d_comp = rescued_component(row.get(d) or "")
        s_comp = rescued_component(row.get(s) or "")
        if d_comp == s_comp:
            continue
        mixed.append(i)
        logger.warning(
            "Mixed MNP rescue at %s:%s %s>%s — duplex %s, simplex %s%s",
            row["Chromosome"],
            row["Start_Position"],
            row["Reference_Allele"],
            row["Tumor_Seq_Allele2"],
            f"reports component {d_comp}" if d_comp else "reports the MNP",
            f"reports component {s_comp}" if s_comp else "reports the MNP",
            (
                "; its simplex_duplex_* columns are NA (they would add different alleles)"
                if combined
                else ""
            ),
        )
    if not mixed:
        return result
    logger.warning(
        "Mixed MNP rescue: %d row(s) where duplex and simplex rescue outcomes differ "
        "(see gbcms_rescue per flavor)%s",
        len(mixed),
        "; their simplex_duplex_* columns are NA" if combined else "",
    )
    combined_cols = [c for c in result.columns if c.startswith("simplex_duplex_")]
    if not combined_cols:
        return result
    is_mixed = pl.int_range(0, result.height).is_in(mixed)
    return result.with_columns(
        [
            pl.when(is_mixed).then(pl.lit("NA")).otherwise(pl.col(c).cast(pl.Utf8)).alias(c)
            for c in combined_cols
        ]
    )


def _validate_variant_key(
    columns: list[str],
    bam_type: str,
    path: Path,
) -> None:
    """Validate that all variant key columns exist in the MAF.

    Raises:
        ValueError: If any key column is missing, with actionable message.
    """
    missing = [k for k in VARIANT_KEY if k not in columns and k != "End_Position"]
    if missing:
        raise ValueError(
            f"Input MAF for '{bam_type}' ({path}) is missing variant key "
            f"columns: {missing}. Expected all of: "
            f"{[k for k in VARIANT_KEY if k != 'End_Position']} (End_Position is optional)"
        )


def _writer_prefix(columns: set[str]) -> str:
    """The ``--column-prefix`` an input's counts were written with: ``""``,
    ``duplex_`` as the pipeline runs it, ``t_`` for legacy names. Found from the
    core counts every gbcms output has (``ref_count`` and ``alt_count``); a name
    that is itself a gbcms column (``mfsd_ref_count``) is not a prefixed count."""
    if {"ref_count", "alt_count"} <= columns:
        return ""
    found = [
        col[: -len("ref_count")]
        for col in columns
        if col.endswith("ref_count")
        and col not in ALL_GBCMS_BASENAMES
        and f"{col[: -len('ref_count')]}alt_count" in columns
    ]
    return min(found, key=len) if found else ""


def _build_rename_map(columns: list[str], bam_type: str) -> dict[str, str]:
    """Rename every gbcms column of an input to ``{bam_type}_{basename}``.

    The writer's ``--column-prefix`` sits on the counts and normalization
    columns only; status, diagnostics, strand bias, mFSD and RNA columns are
    written unprefixed. Each column is found under the name the writer gave it,
    so an input written with ``--column-prefix duplex_`` (as the pipeline runs)
    or ``t_`` keeps all its gbcms columns per input, not only its counts.

    Args:
        columns: List of column names in the MAF.
        bam_type: BAM type label (e.g., "duplex", "simplex").

    Returns:
        Dict mapping original column name → prefixed column name.
    """
    target = f"{bam_type}_"
    present = set(columns)
    writer = _writer_prefix(present)
    prefixed = gbcms_prefixed_basenames()
    rename_map: dict[str, str] = {}
    for base in sorted(ALL_GBCMS_BASENAMES):
        src = f"{writer}{base}" if base in prefixed else base
        if src in present and src != f"{target}{base}" and f"{target}{base}" not in present:
            rename_map[src] = f"{target}{base}"

    if not rename_map and not any(f"{target}{b}" in present for b in ALL_GBCMS_BASENAMES):
        logger.warning(
            "No gbcms count columns found in '%s' MAF. " "Available columns: %s",
            bam_type,
            columns[:10],
        )

    return rename_map


def _strip_prefix(col: str) -> str:
    """Strip a type prefix (e.g., 'duplex_') from a column name.

    Handles the pattern ``{type}_{basename}`` by returning everything
    after the first underscore. Returns the original string if no
    underscore is present.
    """
    parts = col.split("_", 1)
    return parts[1] if len(parts) > 1 else col


def _is_prefixed_gbcms_col(col: str, bam_type: str) -> bool:
    """Check if a column is a gbcms column prefixed with the given BAM type.

    Args:
        col: Column name to check.
        bam_type: Expected BAM type prefix.

    Returns:
        True if column matches ``{bam_type}_{gbcms_basename}`` pattern.
    """
    prefix = f"{bam_type}_"
    if not col.startswith(prefix):
        return False
    basename = col[len(prefix) :]
    return basename in ALL_GBCMS_BASENAMES


def _add_combined_columns(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Add additive simplex_duplex_* combined columns.

    Computes the sum of simplex and duplex counts at ALL levels:
      - Read-level: ref_count, alt_count → total_count, vaf
      - Read-level strand: ref/alt_count_forward/reverse → strand_bias
      - Fragment-level: ref/alt_count_fragment → total_count_fragment, vaf_fragment
      - Fragment strand: ref/alt_count_fragment_forward/reverse → fragment_strand_bias

    Strand bias p-values and odds ratios are computed using the Rust
    ``fisher_exact_2x2`` function (same implementation as the per-BAM engine),
    ensuring numerical consistency.

    Column order matches the canonical layout in ``output.py``:
      §2 counts → §3 read strand + SB → §4 fragment counts → §5 fragment strand + FSB

    Args:
        lf: LazyFrame with simplex_* and duplex_* columns.

    Returns:
        LazyFrame with 20 additional simplex_duplex_* columns appended.
    """

    # ── Helper: cast + null-fill for a combined sum ───────────────────────
    def _num(col: str) -> pl.Expr:
        """String count column → Int64, tolerating float formatting.

        A pandas/R round-trip renders integer columns as '12.0'; a direct
        Int64 cast would null that, so cast through Float64 and round. A
        missing or non-numeric cell (empty, 'NA', 'nan', 'inf', text) is null:
        a missing count is not a zero, so every combined value built from it
        is NA (and counted in a warning per column below)."""
        return pl.col(col).cast(pl.Float64, strict=False).round(0).cast(pl.Int64, strict=False)

    def _sum(metric: str) -> pl.Expr:
        """Sum simplex_{metric} + duplex_{metric}, casting from string."""
        return (_num(f"simplex_{metric}") + _num(f"duplex_{metric}")).alias(
            f"simplex_duplex_{metric}"
        )

    def _total(ref_metric: str, alt_metric: str, total_name: str) -> pl.Expr:
        """Compute total = combined_ref + combined_alt."""
        return (
            pl.col(f"simplex_duplex_{ref_metric}") + pl.col(f"simplex_duplex_{alt_metric}")
        ).alias(f"simplex_duplex_{total_name}")

    def _vaf(alt_metric: str, total_name: str, vaf_name: str) -> pl.Expr:
        """Compute VAF = combined_alt / combined_total, 0/0 → 0.0, NA → NA."""
        alt = pl.col(f"simplex_duplex_{alt_metric}").cast(pl.Float64)
        total = pl.col(f"simplex_duplex_{total_name}").cast(pl.Float64)
        return (
            pl.when(total.is_null())
            .then(None)
            .when(total > 0)
            .then(alt / total)
            .otherwise(0.0)
            .alias(f"simplex_duplex_{vaf_name}")
        )

    # ── Determine which additive metrics exist in the merged output ─────
    schema_names = set(lf.collect_schema().names())

    def _has_both(metric: str) -> bool:
        """Check if both simplex_{metric} and duplex_{metric} are present."""
        return f"simplex_{metric}" in schema_names and f"duplex_{metric}" in schema_names

    available_additive = [m for m in COMBINED_ADDITIVE_ALL if _has_both(m)]
    if not available_additive:
        logger.warning(
            "No additive count columns found in merged output — "
            "skipping combined column computation"
        )
        return lf

    skipped = set(COMBINED_ADDITIVE_ALL) - set(available_additive)
    if skipped:
        logger.info(
            "  Skipping %d additive metrics not in input: %s",
            len(skipped),
            sorted(skipped),
        )

    # Materialize once: the bad-cell count below and the strand-bias pass would
    # otherwise each run the scan and join again.
    lf = lf.collect().lazy()

    # A count cell neither absent (a row the input lacks counts 0, filled
    # above) nor a finite number: the combined value is NA. One warning per
    # column, with the number of rows.
    def _bad(col: str) -> pl.Expr:
        v = pl.col(col).cast(pl.Float64, strict=False)
        return v.is_null() | ~v.is_finite()

    bad = (
        lf.select(
            [
                (_bad(f"simplex_{m}") | _bad(f"duplex_{m}")).sum().alias(m)
                for m in available_additive
            ]
        )
        .collect()
        .row(0, named=True)
    )
    for metric, n in bad.items():
        if n:
            logger.warning(
                "  %s: %d row(s) have a missing or non-numeric count in duplex or simplex; "
                "simplex_duplex_%s (and the totals and VAF built from it) is NA there",
                metric,
                n,
                metric,
            )

    # ── Phase 1: Additive sums (lazy, vectorized) ────────────────────────
    sum_exprs = [_sum(m) for m in available_additive]
    lf = lf.with_columns(sum_exprs)
    logger.info("  Computed %d additive simplex_duplex sums", len(sum_exprs))

    # ── Phase 2a: Derived totals (lazy, vectorized) ────────────────────────
    total_exprs = []
    vaf_exprs = []

    # §2: Read-level total + VAF (only if ref_count + alt_count were summed)
    if "ref_count" in available_additive and "alt_count" in available_additive:
        total_exprs.append(_total("ref_count", "alt_count", "total_count"))
        vaf_exprs.append(_vaf("alt_count", "total_count", "vaf"))

    # §4: Fragment-level total + VAF (only if fragment counts were summed)
    if "ref_count_fragment" in available_additive and "alt_count_fragment" in available_additive:
        total_exprs.append(
            _total("ref_count_fragment", "alt_count_fragment", "total_count_fragment")
        )
        vaf_exprs.append(_vaf("alt_count_fragment", "total_count_fragment", "vaf_fragment"))

    if total_exprs:
        lf = lf.with_columns(total_exprs)
    # ── Phase 2b: VAFs (separate pass — depends on totals from 2a) ────────
    if vaf_exprs:
        lf = lf.with_columns(vaf_exprs)
    if total_exprs or vaf_exprs:
        logger.info(
            "  Computed %d derived totals + %d VAFs",
            len(total_exprs),
            len(vaf_exprs),
        )

    # ── Phase 3: Strand bias via Rust Fisher exact test (eager, per-row) ──
    # Only compute when all 4 directional columns are available
    has_read_strand = all(m in available_additive for m in COMBINED_ADDITIVE_READ_STRAND)
    has_frag_strand = all(m in available_additive for m in COMBINED_ADDITIVE_FRAGMENT_STRAND)

    if has_read_strand or has_frag_strand:
        df = lf.collect()
        df = _compute_combined_strand_bias(
            df,
            compute_read_sb=has_read_strand,
            compute_fragment_sb=has_frag_strand,
        )
        lf = df.lazy()

    return _write_combined_as_text(lf)


def _write_combined_as_text(lf: pl.LazyFrame) -> pl.LazyFrame:
    """The combined counts and VAFs as the writers write theirs: integers, VAFs
    as ``f"{v:.4f}"``, NA where a value is missing. (Strand bias is formatted
    where it is computed.)"""
    names = lf.collect_schema().names()
    exprs = []
    for col in names:
        if not col.startswith("simplex_duplex_") or "strand_bias" in col:
            continue
        if col in ("simplex_duplex_vaf", "simplex_duplex_vaf_fragment"):
            text = pl.col(col).map_elements(lambda v: f"{v:.4f}", return_dtype=pl.Utf8)
            exprs.append(text.fill_null("NA").alias(col))
        else:
            exprs.append(pl.col(col).cast(pl.Utf8).fill_null("NA").alias(col))
    return lf.with_columns(exprs) if exprs else lf


def _compute_combined_strand_bias(
    df: pl.DataFrame,
    *,
    compute_read_sb: bool = True,
    compute_fragment_sb: bool = True,
) -> pl.DataFrame:
    """Compute Fisher strand bias on combined simplex+duplex counts.

    Runs the Rust ``fisher_exact_2x2`` on the 2×2 contingency table:

    .. code-block:: text

                    Forward                     Reverse
        Ref    simplex_duplex_ref_fwd    simplex_duplex_ref_rev
        Alt    simplex_duplex_alt_fwd    simplex_duplex_alt_rev

    Produces up to 4 columns (read-level SB + fragment-level FSB),
    depending on which directional columns are available.

    Args:
        df: DataFrame with simplex_duplex_*_forward/reverse columns.
        compute_read_sb: Whether to compute read-level strand bias.
        compute_fragment_sb: Whether to compute fragment-level strand bias.

    Returns:
        DataFrame with additional strand bias columns.
    """
    from gbcms._rs import fisher_exact_2x2

    # ── Read-level strand bias ────────────────────────────────────────────
    if compute_read_sb:
        sb_results = _apply_fisher(
            df,
            ref_fwd="simplex_duplex_ref_count_forward",
            ref_rev="simplex_duplex_ref_count_reverse",
            alt_fwd="simplex_duplex_alt_count_forward",
            alt_rev="simplex_duplex_alt_count_reverse",
            fisher_fn=fisher_exact_2x2,
        )
        df = df.with_columns(
            [
                pl.Series(
                    "simplex_duplex_strand_bias_p_value", [_fmt_sci(v) for v in sb_results[0]]
                ),
                pl.Series(
                    "simplex_duplex_strand_bias_odds_ratio", [_fmt(v) for v in sb_results[1]]
                ),
            ]
        )

    # ── Fragment-level strand bias ────────────────────────────────────────
    if compute_fragment_sb:
        fsb_results = _apply_fisher(
            df,
            ref_fwd="simplex_duplex_ref_count_fragment_forward",
            ref_rev="simplex_duplex_ref_count_fragment_reverse",
            alt_fwd="simplex_duplex_alt_count_fragment_forward",
            alt_rev="simplex_duplex_alt_count_fragment_reverse",
            fisher_fn=fisher_exact_2x2,
        )
        df = df.with_columns(
            [
                pl.Series(
                    "simplex_duplex_fragment_strand_bias_p_value",
                    [_fmt_sci(v) for v in fsb_results[0]],
                ),
                pl.Series(
                    "simplex_duplex_fragment_strand_bias_odds_ratio",
                    [_fmt(v) for v in fsb_results[1]],
                ),
            ]
        )
    # ── Sanitize NaN/Inf in strand bias columns ─────────────────────────────
    # Fisher exact test returns NaN for OR when alt_total ≤ 1.
    # Polars writes NaN as literal 'NaN' in CSV — convert to 'NA' for MAF.
    sb_cols = [c for c in df.columns if "strand_bias" in c and c.startswith("simplex_duplex_")]
    if sb_cols:
        df = df.with_columns(
            [
                pl.col(c).cast(pl.Utf8).str.replace("NaN", "NA").str.replace("inf", "NA")
                for c in sb_cols
            ]
        )
        logger.debug(
            "Sanitized %d combined strand bias columns (NaN/inf → NA)",
            len(sb_cols),
        )

    logger.debug(
        "Computed combined strand bias for %d variants (read-level + fragment-level)",
        df.height,
    )
    return df


def _apply_fisher(
    df: pl.DataFrame,
    *,
    ref_fwd: str,
    ref_rev: str,
    alt_fwd: str,
    alt_rev: str,
    fisher_fn,
) -> tuple[list[float], list[float]]:
    """Apply Fisher's exact test row-by-row on strand count columns.

    Args:
        df: DataFrame with the strand count columns.
        ref_fwd/ref_rev/alt_fwd/alt_rev: Column names for the 2×2 table.
        fisher_fn: Callable(a, b, c, d) → (p_value, odds_ratio).

    Returns:
        Tuple of (p_values_list, odds_ratios_list).
    """
    rf = df[ref_fwd].to_list()
    rr = df[ref_rev].to_list()
    af = df[alt_fwd].to_list()
    ar = df[alt_rev].to_list()

    p_values: list[float] = []
    odds_ratios: list[float] = []

    for i in range(df.height):
        # Values are Int64 from the additive sum phase
        # A missing combined count (NA) leaves no table to test: NA.
        if None in (rf[i], rr[i], af[i], ar[i]):
            p_values.append(float("nan"))
            odds_ratios.append(float("nan"))
            continue
        p, odds = fisher_fn(int(rf[i]), int(rr[i]), int(af[i]), int(ar[i]))
        p_values.append(p)
        odds_ratios.append(odds)

    return p_values, odds_ratios


def _apply_legacy_naming(
    df: pl.DataFrame,
    types: list[str],
) -> pl.DataFrame:
    """Rename ``{type}_{metric}`` → ``t_{metric}_{type}`` for genotype_variants compat.

    This reproduces the column naming convention used by the legacy
    ``create_duplex_simplex_dataframe.py`` in genotype_variants.

    Args:
        df: DataFrame with ``{type}_{metric}`` column names.
        types: List of BAM type labels used in the merge.

    Returns:
        DataFrame with renamed columns.
    """
    rename_map: dict[str, str] = {}
    for t in types:
        prefix = f"{t}_"
        for col in df.columns:
            if col.startswith(prefix):
                basename = col[len(prefix) :]
                if basename in ALL_GBCMS_BASENAMES:
                    rename_map[col] = f"t_{basename}_{t}"

    # Also rename combined columns if present
    for col in df.columns:
        if col.startswith("simplex_duplex_"):
            basename = col[len("simplex_duplex_") :]
            rename_map[col] = f"t_{basename}_simplex_duplex"

    if rename_map:
        df = df.rename(rename_map)
        logger.debug("Legacy rename: %d columns renamed", len(rename_map))

    return df
