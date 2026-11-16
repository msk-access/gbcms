"""
I/O module for gbcms.

Provides readers and writers for variant files (VCF, MAF format)
and batch helpers for post-processing operations.

Streaming I/O (input/output):
    Row-by-row processing through the Rust counting engine.
    Uses ``csv`` stdlib for zero-memory streaming.

Batch I/O (batch):
    Complete-file operations for merge, reports, and analytics.
    Uses Polars for efficient joins and aggregation.
"""

from typing import TYPE_CHECKING, Any

from .input import MafReader, VariantReader, VcfReader
from .output import MafWriter, OutputWriter, VcfWriter

if TYPE_CHECKING:
    from .batch import read_maf, read_parquet, scan_maf, write_maf

__all__ = [
    # Streaming I/O (core pipeline)
    "MafReader",
    "MafWriter",
    "OutputWriter",
    "VariantReader",
    "VcfReader",
    "VcfWriter",
    # Batch I/O (post-processing)
    "read_maf",
    "read_parquet",
    "scan_maf",
    "write_maf",
]


def __getattr__(name: str) -> Any:
    """The batch helpers are read on first use: importing them loads polars (about 0.14 s),
    which a counting run never needs."""
    if name in ("read_maf", "read_parquet", "scan_maf", "write_maf"):
        from . import batch

        return getattr(batch, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
