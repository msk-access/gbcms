"""Report generation subpackage for gbcms.

Currently provides:
- mFSD per-variant HTML report (interactive fragment size distributions)
"""

from .mfsd_report import generate_mfsd_report

__all__ = ["generate_mfsd_report"]
