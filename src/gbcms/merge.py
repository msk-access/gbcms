"""
Multi-BAM MAF merge engine using Polars.

Merges per-BAM-type genotyped MAFs (e.g., duplex, simplex) into a single
output MAF with type-prefixed count columns.  Uses Polars lazy API for
WGS-scale performance.

Architecture:
    1. Scan each input MAF lazily via ``io.batch.scan_maf``
    2. Detect whether columns are already prefixed or need renaming
    3. Progressive outer join on the 5-column variant key
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
import time
from pathlib import Path

import polars as pl

from gbcms.io.batch import scan_maf, write_maf
from gbcms.models.core import MergeConfig

__all__ = ["merge_mafs"]

logger = logging.getLogger(__name__)


# ── Constants (single source of truth — canonical basenames from output.py) ──

# 5-column variant key used for outer joins.
VARIANT_KEY: list[str] = [
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
]

# gbcms count column basenames (without any prefix).
# These columns contain numeric counts/rates and will be type-prefixed.
# Reference: MafWriter._gbcms_column_names() in io/output.py lines 190-219.
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
    "gbcms_diagnostic",
    "gbcms_rescue",
    "strand_bias_p_value",
    "strand_bias_odds_ratio",
    "fragment_strand_bias_p_value",
    "fragment_strand_bias_odds_ratio",
]

# Combined set of all gbcms basenames for detection.
ALL_GBCMS_BASENAMES: set[str] = set(GBCMS_COUNT_BASENAMES) | set(GBCMS_META_BASENAMES)

# Fragment-level metrics used for simplex+duplex additive combination.
# These are the only columns where addition is semantically meaningful
# (read-level counts are NOT additive across BAM types).
COMBINED_ADDITIVE: list[str] = [
    "ref_count_fragment",
    "alt_count_fragment",
]


def merge_mafs(config: MergeConfig) -> None:
    """Merge per-BAM-type genotyped MAFs into a single type-prefixed output.

    Performs an outer join on the 5-column variant key, prefixes gbcms count
    columns with the BAM type label, optionally computes combined
    simplex_duplex columns, and writes the merged result.

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
    for bam_type, path in config.inputs.items():
        logger.info("Scanning %s MAF: %s", bam_type, path)
        lf = scan_maf(path)

        # Validate variant key columns exist
        schema_names = lf.collect_schema().names()
        _validate_variant_key(schema_names, bam_type, path)

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

        frames[bam_type] = lf

    # ── 2. Progressive outer join ─────────────────────────────────────────────
    types = list(frames.keys())
    merged = frames[types[0]]

    for join_type in types[1:]:
        # Select only variant key + gbcms columns from the joining frame
        # to avoid duplicating annotation columns across inputs.
        join_frame = frames[join_type]
        join_cols = [
            c
            for c in join_frame.collect_schema().names()
            if c in VARIANT_KEY or _is_prefixed_gbcms_col(c, join_type)
        ]
        merged = merged.join(
            join_frame.select(join_cols),
            on=VARIANT_KEY,
            how="full",
            coalesce=True,
        )
        logger.info(
            "  Joined '%s' (%d count cols)",
            join_type,
            len(join_cols) - len(VARIANT_KEY),
        )

    # ── 3. Fill nulls → "0" for count columns, "" for meta columns ─────────
    output_schema = merged.collect_schema().names()
    count_cols_in_output = [
        c
        for c in output_schema
        if any(
            c == f"{t}_{base}" for t in types for base in GBCMS_COUNT_BASENAMES
        )
    ]
    meta_cols_in_output = [
        c
        for c in output_schema
        if any(
            c == f"{t}_{base}" for t in types for base in GBCMS_META_BASENAMES
        )
    ]
    if count_cols_in_output:
        merged = merged.with_columns(
            [pl.col(c).fill_null("0") for c in count_cols_in_output]
        )
        logger.info("  Filled nulls → '0' for %d count columns", len(count_cols_in_output))
    if meta_cols_in_output:
        merged = merged.with_columns(
            [pl.col(c).fill_null("") for c in meta_cols_in_output]
        )
        logger.info("  Filled nulls → '' for %d meta columns", len(meta_cols_in_output))

    # ── 4. Combined simplex+duplex columns ────────────────────────────────────
    if config.add_combined and "simplex" in frames and "duplex" in frames:
        merged = _add_combined_columns(merged)
        logger.info("  Added simplex_duplex combined columns")

    # ── 5. Materialize and write ──────────────────────────────────────────────
    result = merged.collect()
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

    # Log unmatched variant counts per type for monitoring
    # A variant is "unmatched" if its count columns are all "0" for that type
    for t in types:
        t_ref_col = f"{t}_ref_count"
        if t_ref_col in result.columns:
            n_unmatched = result.filter(pl.col(t_ref_col) == "0").height
            if n_unmatched > 0:
                logger.info(
                    "  %d/%d variants have no %s counts (filled with 0)",
                    n_unmatched,
                    result.height,
                    t,
                )

    # Legacy naming pass (rename {type}_{metric} → t_{metric}_{type})
    if config.legacy_naming:
        result = _apply_legacy_naming(result, types)
        logger.info("  Applied legacy t_{metric}_{type} naming")

    write_maf(result, config.output)
    elapsed = time.perf_counter() - t_start
    logger.info("Merge complete in %.1fs", elapsed)


# ── Helpers ──────────────────────────────────────────────────────────────────


def _validate_variant_key(
    columns: list[str],
    bam_type: str,
    path: Path,
) -> None:
    """Validate that all variant key columns exist in the MAF.

    Raises:
        ValueError: If any key column is missing, with actionable message.
    """
    missing = [k for k in VARIANT_KEY if k not in columns]
    if missing:
        raise ValueError(
            f"Input MAF for '{bam_type}' ({path}) is missing variant key "
            f"columns: {missing}. Expected all of: {VARIANT_KEY}"
        )


def _build_rename_map(columns: list[str], bam_type: str) -> dict[str, str]:
    """Build a column rename map: unprefixed gbcms cols → type-prefixed.

    Detects whether columns are already prefixed with the bam_type label.
    If already prefixed (e.g., ``duplex_ref_count``), returns empty dict
    to avoid double-prefixing.

    Args:
        columns: List of column names in the MAF.
        bam_type: BAM type label (e.g., "duplex", "simplex").

    Returns:
        Dict mapping original column name → prefixed column name.
        Empty dict if columns are already prefixed.
    """
    prefix = f"{bam_type}_"

    # Check if columns are already prefixed
    already_prefixed = any(c.startswith(prefix) for c in columns if _strip_prefix(c) in ALL_GBCMS_BASENAMES)
    if already_prefixed:
        logger.debug(
            "Columns in '%s' MAF already have '%s' prefix",
            bam_type,
            prefix,
        )
        return {}

    # Build rename map for unprefixed gbcms columns
    rename_map: dict[str, str] = {}
    for col in columns:
        if col in ALL_GBCMS_BASENAMES:
            rename_map[col] = f"{prefix}{col}"

    if not rename_map:
        logger.warning(
            "No gbcms count columns found in '%s' MAF. "
            "Available columns: %s",
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
    basename = col[len(prefix):]
    return basename in ALL_GBCMS_BASENAMES


def _add_combined_columns(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Add additive simplex_duplex_* combined columns.

    Computes the sum of simplex and duplex fragment-level counts for
    each metric in ``COMBINED_ADDITIVE``. Also computes a combined VAF.

    Only fragment-level metrics are summed — read-level counts are NOT
    additive across BAM types (a read can only come from one BAM).

    Args:
        lf: LazyFrame with simplex_* and duplex_* columns.

    Returns:
        LazyFrame with additional simplex_duplex_* columns appended.
    """
    exprs: list[pl.Expr] = []

    for metric in COMBINED_ADDITIVE:
        simplex_col = f"simplex_{metric}"
        duplex_col = f"duplex_{metric}"
        combined_col = f"simplex_duplex_{metric}"

        # Cast string → Int64 for arithmetic, fill null → 0
        exprs.append(
            (
                pl.col(simplex_col).cast(pl.Int64, strict=False).fill_null(0)
                + pl.col(duplex_col).cast(pl.Int64, strict=False).fill_null(0)
            ).alias(combined_col)
        )

    # Compute combined total and VAF
    exprs.append(
        (
            pl.col("simplex_ref_count_fragment").cast(pl.Int64, strict=False).fill_null(0)
            + pl.col("duplex_ref_count_fragment").cast(pl.Int64, strict=False).fill_null(0)
            + pl.col("simplex_alt_count_fragment").cast(pl.Int64, strict=False).fill_null(0)
            + pl.col("duplex_alt_count_fragment").cast(pl.Int64, strict=False).fill_null(0)
        ).alias("simplex_duplex_total_count_fragment")
    )

    # VAF = alt / (alt + ref), handle division by zero → 0.0
    combined_alt = (
        pl.col("simplex_alt_count_fragment").cast(pl.Int64, strict=False).fill_null(0)
        + pl.col("duplex_alt_count_fragment").cast(pl.Int64, strict=False).fill_null(0)
    )
    combined_total = (
        pl.col("simplex_ref_count_fragment").cast(pl.Int64, strict=False).fill_null(0)
        + pl.col("duplex_ref_count_fragment").cast(pl.Int64, strict=False).fill_null(0)
        + combined_alt
    )
    exprs.append(
        pl.when(combined_total > 0)
        .then(combined_alt.cast(pl.Float64) / combined_total.cast(pl.Float64))
        .otherwise(0.0)
        .alias("simplex_duplex_vaf_fragment")
    )

    return lf.with_columns(exprs)


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
                basename = col[len(prefix):]
                if basename in ALL_GBCMS_BASENAMES:
                    rename_map[col] = f"t_{basename}_{t}"

    # Also rename combined columns if present
    for col in df.columns:
        if col.startswith("simplex_duplex_"):
            basename = col[len("simplex_duplex_"):]
            rename_map[col] = f"t_{basename}_simplex_duplex"

    if rename_map:
        df = df.rename(rename_map)
        logger.debug("Legacy rename: %d columns renamed", len(rename_map))

    return df
