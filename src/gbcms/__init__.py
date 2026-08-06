"""
gbcms (Get Base Counts Multi-Sample) - A tool for counting bases at variant positions.

This package provides a command-line interface and Python API for genotyping
variants in BAM files using a high-performance Rust counting engine.

Example usage:
    $ gbcms dna -v variants.vcf -b sample.bam -f reference.fa -o output/
"""

__version__ = "6.1.0"

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
