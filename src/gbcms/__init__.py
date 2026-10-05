"""
gbcms (Get Base Counts Multi-Sample) - A tool for counting bases at variant positions.

This package provides a command-line interface and Python API for genotyping
variants in BAM files using a high-performance Rust counting engine.

Example usage:
    $ gbcms dna -v variants.vcf -b sample.bam -f reference.fa -o output/
"""

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _package_version

# The installed package's version, which maturin takes from rust/Cargo.toml: one
# source for the package, the extension and every output's provenance line.
try:
    __version__ = _package_version("gbcms")
except PackageNotFoundError:  # a source tree that was never installed
    __version__ = "unknown"

from .merge import merge_mafs
from .models.core import (
    GbcmsBaseConfig,
    GbcmsConfig,
    GbcmsDnaConfig,
    GbcmsRnaConfig,
    MergeConfig,
    OutputFormat,
    Variant,
    VariantType,
)
from .observations import ObservationResult, observe_molecules
from .pipeline import Pipeline

__all__ = [
    "__version__",
    "GbcmsBaseConfig",
    "GbcmsConfig",
    "GbcmsDnaConfig",
    "GbcmsRnaConfig",
    "MergeConfig",
    "OutputFormat",
    "Pipeline",
    "Variant",
    "VariantType",
    "ObservationResult",
    "merge_mafs",
    "observe_molecules",
]
