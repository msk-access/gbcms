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
    "gbcms_status_reason",
    "gbcms_diagnostic",
    "gbcms_rescue",
    "strand_bias_p_value",
    "strand_bias_odds_ratio",
    "fragment_strand_bias_p_value",
    "fragment_strand_bias_odds_ratio",
]

# Combined set of all gbcms basenames for detection.
ALL_GBCMS_BASENAMES: set[str] = set(GBCMS_COUNT_BASENAMES) | set(GBCMS_META_BASENAMES)

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
        if any(c == f"{t}_{base}" for t in types for base in GBCMS_COUNT_BASENAMES)
    ]
    meta_cols_in_output = [
        c
        for c in output_schema
        if any(c == f"{t}_{base}" for t in types for base in GBCMS_META_BASENAMES)
    ]
    if count_cols_in_output:
        merged = merged.with_columns([pl.col(c).fill_null("0") for c in count_cols_in_output])
        logger.info("  Filled nulls → '0' for %d count columns", len(count_cols_in_output))
    if meta_cols_in_output:
        merged = merged.with_columns([pl.col(c).fill_null("") for c in meta_cols_in_output])
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
    already_prefixed = any(
        c.startswith(prefix) for c in columns if _strip_prefix(c) in ALL_GBCMS_BASENAMES
    )
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
    def _sum(metric: str) -> pl.Expr:
        """Sum simplex_{metric} + duplex_{metric}, casting from string."""
        return (
            pl.col(f"simplex_{metric}").cast(pl.Int64, strict=False).fill_null(0)
            + pl.col(f"duplex_{metric}").cast(pl.Int64, strict=False).fill_null(0)
        ).alias(f"simplex_duplex_{metric}")

    def _total(ref_metric: str, alt_metric: str, total_name: str) -> pl.Expr:
        """Compute total = combined_ref + combined_alt."""
        return (
            pl.col(f"simplex_duplex_{ref_metric}") + pl.col(f"simplex_duplex_{alt_metric}")
        ).alias(f"simplex_duplex_{total_name}")

    def _vaf(alt_metric: str, total_name: str, vaf_name: str) -> pl.Expr:
        """Compute VAF = combined_alt / combined_total, 0/0 → 0.0."""
        alt = pl.col(f"simplex_duplex_{alt_metric}").cast(pl.Float64)
        total = pl.col(f"simplex_duplex_{total_name}").cast(pl.Float64)
        return (pl.when(total > 0).then(alt / total).otherwise(0.0)).alias(
            f"simplex_duplex_{vaf_name}"
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
        return df.lazy()

    return lf


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
                pl.Series("simplex_duplex_strand_bias_p_value", sb_results[0]),
                pl.Series("simplex_duplex_strand_bias_odds_ratio", sb_results[1]),
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
                pl.Series("simplex_duplex_fragment_strand_bias_p_value", fsb_results[0]),
                pl.Series("simplex_duplex_fragment_strand_bias_odds_ratio", fsb_results[1]),
            ]
        )
    # ── Sanitize NaN/Inf in strand bias columns ─────────────────────────────
    # Fisher exact test returns NaN for OR when alt_total ≤ 1 (issue #19).
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
        a = int(rf[i]) if rf[i] is not None else 0
        b = int(rr[i]) if rr[i] is not None else 0
        c = int(af[i]) if af[i] is not None else 0
        d = int(ar[i]) if ar[i] is not None else 0
        p, odds = fisher_fn(a, b, c, d)
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
