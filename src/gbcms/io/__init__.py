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

from .batch import read_maf, read_parquet, scan_maf, write_maf
from .input import MafReader, VariantReader, VcfReader
from .output import MafWriter, OutputWriter, VcfWriter

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
