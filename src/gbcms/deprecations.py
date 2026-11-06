"""Deprecation notices shared by the CLI, the Pipeline and the Python API."""

from __future__ import annotations

import logging
import warnings

__all__ = ["warn_rescue_homopolymer_deprecated"]

logger = logging.getLogger(__name__)

_rescue_homopolymer_notice = (
    "--rescue-homopolymer (rescue_homopolymer) is deprecated in 6.7.0 and removed in 6.8.0; "
    "it is ignored, and the row counts the given allele. It did nothing useful the "
    "diagnostics do not: OBSERVED_ALLELE and COEXISTING_ALLELE in gbcms_diagnostic name the "
    "allele the reads carry."
)


def warn_rescue_homopolymer_deprecated(stacklevel: int = 3) -> None:
    """Warn once per call that the homopolymer twin is deprecated and ignored."""
    logger.warning(_rescue_homopolymer_notice)
    warnings.warn(_rescue_homopolymer_notice, DeprecationWarning, stacklevel=stacklevel)
