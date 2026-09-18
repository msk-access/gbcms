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


def _validate_rectangular(path: Path, *, comment_prefix: str = "#") -> None:
    """Reject MAF files whose data rows do not match the header width.

    The polars parser raises for rows with MORE fields than the header but
    silently null-pads rows with FEWER — and a row that lost its trailing
    fields has lost exactly the gbcms count columns, which downstream
    null-filling would silently turn into zeros. One cheap text pass over
    the file catches both shapes with the offending line number.
    """
    header_width: int | None = None
    with open(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            if line.startswith(comment_prefix):
                continue
            width = line.rstrip("\n").count("\t") + 1
            if header_width is None:
                header_width = width
                continue
            if width != header_width:
                raise ValueError(
                    f"{path}: line {line_no} has {width} field(s) but the "
                    f"header has {header_width} — a ragged row would corrupt "
                    "the trailing count columns. Fix or remove the row."
                )


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

    _validate_rectangular(path, comment_prefix=comment_prefix)
    logger.debug("Reading MAF (batch): %s", path)
    df = pl.read_csv(
        path,
        separator="\t",
        comment_prefix=comment_prefix,
        infer_schema_length=0,  # All columns as strings
        # Over-length rows raise in the parser; under-length rows would be
        # silently null-padded by polars, so _validate_rectangular above
        # rejects both shapes before parsing — the gbcms count columns are
        # the LAST columns of the file, exactly what a lost tail would
        # silently zero.
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

    _validate_rectangular(path, comment_prefix=comment_prefix)
    logger.debug("Lazy-scanning MAF: %s", path)
    return pl.scan_csv(
        path,
        separator="\t",
        comment_prefix=comment_prefix,
        infer_schema_length=0,
        # Rectangularity is enforced by _validate_rectangular above — see
        # read_maf.
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
