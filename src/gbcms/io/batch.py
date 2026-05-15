"""
Batch I/O helpers using Polars.

Provides functions for reading and writing complete files in batch mode.
Used by post-processing operations (merge, mFSD report) that need to load
entire files for joins, aggregation, or analysis.

NOT for the streaming core pipeline — that uses ``csv`` stdlib via
``io/input.py`` and ``io/output.py`` for memory-safe one-row-at-a-time
processing through the Rust counting engine.

Design rationale:
    - All MAF reads use ``infer_schema_length=0`` (all columns as Utf8)
      because MAF is a text format. Callers cast columns as needed.
    - Parquet reads infer schema from Parquet metadata (native types).
    - ``comment_prefix="#"`` matches the streaming MafReader convention
      in ``io/input.py`` (lines 99-105).
"""

import logging
from pathlib import Path

import polars as pl

__all__ = ["read_maf", "scan_maf", "read_parquet", "write_maf"]

logger = logging.getLogger(__name__)


def read_maf(path: Path, *, comment_prefix: str = "#") -> pl.DataFrame:
    """Read a MAF file into a Polars DataFrame, skipping comment lines.

    All columns are read as strings (``Utf8``) because MAF is a text format
    with mixed types. Callers should cast specific columns as needed for
    numeric operations (e.g., ``col.cast(pl.Int64)`` for count columns).

    Args:
        path: Path to the MAF file.
        comment_prefix: Lines starting with this string are skipped.
            Defaults to ``"#"`` (standard MAF comment convention).

    Returns:
        DataFrame with all columns as Utf8.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
        polars.exceptions.ComputeError: If file cannot be parsed as TSV.
    """
    if not path.exists():
        raise FileNotFoundError(f"MAF file not found: {path}")

    logger.debug("Reading MAF (batch): %s", path)
    df = pl.read_csv(
        path,
        separator="\t",
        comment_prefix=comment_prefix,
        infer_schema_length=0,  # All columns as strings
        truncate_ragged_lines=True,
    )
    logger.info("Loaded MAF: %s (%d rows × %d cols)", path.name, df.height, df.width)
    return df


def scan_maf(path: Path, *, comment_prefix: str = "#") -> pl.LazyFrame:
    """Lazy-scan a MAF file for deferred processing.

    Returns a ``LazyFrame`` that is not materialized until ``.collect()``
    is called. Use this for large-file joins where Polars can optimize
    the query plan before execution.

    Args:
        path: Path to the MAF file.
        comment_prefix: Lines starting with this string are skipped.

    Returns:
        LazyFrame for deferred execution.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
    """
    if not path.exists():
        raise FileNotFoundError(f"MAF file not found: {path}")

    logger.debug("Lazy-scanning MAF: %s", path)
    return pl.scan_csv(
        path,
        separator="\t",
        comment_prefix=comment_prefix,
        infer_schema_length=0,
        truncate_ragged_lines=True,
    )


def read_parquet(path: Path) -> pl.DataFrame:
    """Read a Parquet file into a Polars DataFrame.

    Schema is inferred from Parquet file metadata (native Arrow types).
    No string coercion — columns retain their original types (int, float,
    list, etc.).

    Args:
        path: Path to the Parquet file.

    Returns:
        DataFrame with schema inferred from Parquet metadata.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
    """
    if not path.exists():
        raise FileNotFoundError(f"Parquet file not found: {path}")

    logger.debug("Reading Parquet (batch): %s", path)
    df = pl.read_parquet(path)
    logger.info("Loaded Parquet: %s (%d rows × %d cols)", path.name, df.height, df.width)
    return df


def write_maf(df: pl.DataFrame, path: Path) -> None:
    """Write a Polars DataFrame as a tab-separated MAF file.

    Writes all columns as-is with tab separator. No comment header is
    added — the output is a plain TSV with a single header row.

    Args:
        df: DataFrame to write.
        path: Output file path.
    """
    df.write_csv(path, separator="\t")
    logger.info("Wrote MAF: %s (%d rows × %d cols)", path.name, df.height, df.width)
