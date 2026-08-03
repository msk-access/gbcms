"""Per-molecule observation export — the supported public API.

`BaseCounts` aggregates every fragment covering a variant into totals. Some analyses need
the layer underneath: *which molecule* carried *which allele*. That is what this module
exposes — the same per-fragment call the counting engine already resolves, surfaced instead
of discarded.

Because a molecule observed at two variants in one call carries the same ``molecule_hash``,
a caller can link alleles **across loci** (read-backed phasing, allelic imbalance,
duplex-masking QC). gbcms performs no such linking itself: it reports what each molecule
showed, and the interpretation belongs to the caller.

This module is the **stable surface**. It delegates to ``gbcms._rs``, which is an internal
implementation detail and *not* covered by SemVer — depend on this module, not on ``_rs``.

Example::

    from gbcms import Variant, VariantType, observe_molecules

    variants = [Variant(chrom="17", pos=7676153, ref="T", alt="A",
                        variant_type=VariantType.SNP)]
    result = observe_molecules("sample.bam", variants)
    for obs in result.observations:
        ...  # obs.variant_index, obs.molecule_hash, obs.allele, obs.best_qual
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from gbcms._rs import Observation, count_bam_binned_observations, prepare_variants
from gbcms._rs import Variant as _RsVariant
from gbcms.models.core import QualityThresholds, ReadFilters, Variant

if TYPE_CHECKING:
    from gbcms.models.core import GbcmsDnaConfig

__all__ = [
    "ALLELE_ALT",
    "ALLELE_N",
    "ALLELE_OTHER",
    "ALLELE_REF",
    "ObservationResult",
    "observe_molecules",
]

# `Observation.allele` values. Mirrors the OBS_ALLELE_* constants in rust/src/types.rs.
ALLELE_REF = 0
ALLELE_ALT = 1
ALLELE_N = 2
ALLELE_OTHER = 3


@dataclass(frozen=True)
class ObservationResult:
    """Outcome of :func:`observe_molecules`.

    ``observations`` holds the rows when returned in memory, and is empty when they were
    streamed to ``path`` instead. ``n_rows`` is the true row count either way.
    """

    observations: list[Observation]
    path: Path | None
    n_rows: int


def observe_molecules(
    bam: str | Path,
    variants: Sequence[Variant],
    *,
    config: GbcmsDnaConfig | None = None,
    observations_path: str | Path | None = None,
) -> ObservationResult:
    """Observe the per-molecule allele at each variant.

    Args:
        bam: Indexed BAM/CRAM to read.
        variants: Variants to observe. ``Observation.variant_index`` indexes into this list.
        config: Counting configuration (filters, quality gates, UMI tag, threads). Library
            defaults are used when omitted. When ``config.reference_fasta`` is set, variants
            are normalized first (left-alignment + ref_context), which is what makes
            windowed indel detection work — mirroring :class:`~gbcms.pipeline.Pipeline`.
        observations_path: When given, rows are written to Parquet from Rust rather than
            materialized as Python objects — for panel- or genome-wide runs, where the FFI
            round-trip, not the counting, is the bottleneck.
            *(Reserved: the writer lands in a follow-up; passing a path currently raises.)*

    Returns:
        An :class:`ObservationResult`.

    Note:
        Molecules are grouped by ``hash_molecule(qname, umi)``, so reads must not be
        pre-filtered in a way that drops whole molecules. In particular do **not** enable
        ``filters.improper_pair`` for consensus BAMs whose reads carry no PROPER_PAIR flag
        (some UMI-collapsed pipelines emit none) — every read would be discarded.
    """
    if observations_path is not None:
        raise NotImplementedError(
            "Parquet streaming of observations is not implemented yet; omit "
            "`observations_path` to receive them in memory."
        )

    # Sub-model defaults, so a caller with no config never has to build a full
    # GbcmsDnaConfig (which requires variant/bam paths this entry point does not use).
    filters = config.filters if config is not None else ReadFilters()
    quality = config.quality if config is not None else QualityThresholds()
    threads = config.threads if config is not None else 1

    rs_variants = [_RsVariant(v.chrom, v.pos, v.ref, v.alt, v.variant_type.value) for v in variants]

    # Normalize when a reference is available: left-alignment + ref_context are what let the
    # engine recognise a shifted indel. Same call Pipeline makes before counting.
    fasta = getattr(config, "reference_fasta", None) if config is not None else None
    if fasta is not None:
        prepared = prepare_variants(
            rs_variants,
            str(fasta),
            quality.context_padding,
            False,
            threads,
            quality.adaptive_context,
        )
        rs_variants = [p.variant for p in prepared]

    _counts, observations = count_bam_binned_observations(
        str(bam),
        rs_variants,
        [None] * len(rs_variants),
        min_mapq=quality.min_mapping_quality,
        min_baseq=quality.min_base_quality,
        filter_duplicates=filters.duplicates,
        filter_secondary=filters.secondary,
        filter_supplementary=filters.supplementary,
        filter_qc_failed=filters.qc_failed,
        filter_improper_pair=filters.improper_pair,
        filter_indel=filters.indel,
        threads=threads,
        fragment_qual_threshold=quality.fragment_qual_threshold,
        apply_baq=config.apply_baq if config is not None else False,
        umi_tag=config.umi_tag if config is not None else None,
    )
    return ObservationResult(observations=observations, path=None, n_rows=len(observations))
