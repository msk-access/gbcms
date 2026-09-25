"""The MNP rescue audit contract, defined once for its writer and its readers.

The rescue pass (``pipeline``) writes a ``gbcms_rescue`` value for every rescue
candidate and ``gbcms merge`` reads it back from the per-flavor MAFs. Both sides
use the names and functions here, so renaming an outcome or field cannot leave
one side silently matching nothing.
"""

from typing import Any

# Diagnostic flag that makes an MNP a rescue candidate. Written by the
# per-row diagnostics and read by the rescue pass's candidate gate.
MNP_RESCUE_ELIGIBLE = "MNP_RESCUE_ELIGIBLE"

# Outcomes recorded as ``outcome=<name>`` in ``gbcms_rescue``.
OUTCOME_RESCUED = "rescued"
OUTCOME_SKIPPED_GROUPED = "skipped_grouped"
OUTCOME_HAPLOTYPE_CONFIRMED = "haplotype_confirmed"
OUTCOME_NO_IMPROVEMENT = "no_improvement"
OUTCOME_REF_VALIDATION_FAILED = "ref_validation_failed"

_FIELD_OUTCOME = "outcome"
_FIELD_ADOPTED = "adopted"


def format_rescue_audit(
    outcome: str,
    original: Any,
    positions: list[str] | None = None,
    adopted: str | None = None,
) -> str:
    """Build the ``gbcms_rescue`` value.

    Format: ``method=decomposed;outcome=<outcome>;original_ref=R;original_alt=A;
    original_partial=P;original_confirmed=C[;adopted=<label>]
    [;positions=<label>:<ad|ref_fail>+...]`` where a label is
    ``chrom:pos(REF>ALT)`` (1-based). Positions are joined with ``+``, never
    ``,``: the VCF writer emits this string as the Number=1 ``GR`` INFO value,
    which VCF parsers split at commas. ``original_*`` are the MNP's own counts —
    for a rescued row the only record of its evaluation, because the count
    columns then carry the adopted component's counts. ``original_confirmed``
    is the MNP's ``mnp_confirmed_alt``: reads that showed the whole haplotype.
    """
    parts = [
        "method=decomposed",
        f"{_FIELD_OUTCOME}={outcome}",
        f"original_ref={original.rd}",
        f"original_alt={original.ad}",
        f"original_partial={original.partial_alt}",
        f"original_confirmed={original.mnp_confirmed_alt}",
    ]
    if adopted is not None:
        parts.append(f"{_FIELD_ADOPTED}={adopted}")
    if positions:
        parts.append("positions=" + "+".join(positions))
    return ";".join(parts)


def rescued_component(audit: str) -> str | None:
    """The adopted component label of a ``gbcms_rescue`` value, or None if the
    row was not rescued (any other outcome, or an empty value)."""
    fields = dict(part.split("=", 1) for part in audit.split(";") if "=" in part)
    if fields.get(_FIELD_OUTCOME) != OUTCOME_RESCUED:
        return None
    return fields.get(_FIELD_ADOPTED)
