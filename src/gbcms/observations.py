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
    result = observe_molecules("sample.bam", variants, reference_fasta="hg19.fa")
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
from gbcms.models.core import AlignmentConfig, QualityThresholds, ReadFilters, Variant

if TYPE_CHECKING:
    from gbcms.models.core import GbcmsDnaConfig

__all__ = [
    "ALLELE_ALT",
    "ALLELE_N",
    "ALLELE_OTHER",
    "ALLELE_REF",
    "Observation",
    "ObservationResult",
    "observe_molecules",
]

# `Observation.allele` values. Mirrors the OBS_ALLELE_* constants in rust/src/types.rs.
ALLELE_REF = 0
ALLELE_ALT = 1
ALLELE_N = 2
ALLELE_OTHER = 3


def _parquet_row_count(path: Path) -> int:
    """Row count of a written Parquet file, from its footer metadata (no data read)."""
    try:
        import pyarrow.parquet as pq

        return int(pq.ParquetFile(path).metadata.num_rows)
    except ImportError:  # pragma: no cover - pyarrow is optional for readers
        return -1


@dataclass(frozen=True)
class ObservationResult:
    """Outcome of :func:`observe_molecules`.

    ``observations`` holds the rows when returned in memory, and is empty when they were
    streamed to ``path`` instead; ``n_rows`` is the true row count either way.

    ``variant_status`` is one entry per input variant, positionally aligned with the
    ``variants`` argument (and therefore with ``Observation.variant_index``). Anything other
    than ``"PASS"`` means normalization rejected that variant — e.g. ``REF_MISMATCH`` — and
    its rows should be discarded. It is ``None`` when no reference was supplied, since
    without one there is nothing to validate against.
    """

    observations: list[Observation]
    path: Path | None
    n_rows: int
    variant_status: list[str] | None = None


def observe_molecules(
    bam: str | Path,
    variants: Sequence[Variant],
    *,
    reference_fasta: str | Path | None = None,
    is_maf: bool = False,
    config: GbcmsDnaConfig | None = None,
    observations_path: str | Path | None = None,
) -> ObservationResult:
    """Observe the per-molecule allele at each variant.

    Args:
        bam: Indexed BAM/CRAM to read.
        variants: Variants to observe. ``Observation.variant_index`` indexes into this list,
            positionally — the list is never filtered or reordered, so the join key stays
            meaningful even when a variant fails validation.
        reference_fasta: Reference for normalization (left-alignment, ``ref_context``, and
            indel decomposition). **Strongly recommended for indels**: without it, a
            deletion that the aligner shifted, or one whose decomposed form carries the
            ALT support, is scored REF and its ALT molecules are exported as ``OTHER``.
        is_maf: Set when the variants came from a MAF, whose ``-`` alleles need anchor
            resolution during normalization. Ignored without ``reference_fasta``.
        config: Counting configuration (filters, quality gates, alignment backend, UMI tag,
            threads). Library defaults are used when omitted.
        observations_path: When given, rows are written to Parquet from Rust rather than
            materialized as Python objects — for panel- or genome-wide runs, where the FFI
            round-trip, not the counting, is the bottleneck. The returned
            ``observations`` list is then empty and ``path``/``n_rows`` describe the file,
            whose columns are ``variant_index, chrom, pos, ref, alt, molecule_hash, allele,
            best_qual`` (self-describing, so it stands alone once written).

    Returns:
        An :class:`ObservationResult`. Check ``variant_status`` before trusting rows for a
        given ``variant_index``.

    Note:
        Molecules are grouped by ``hash_molecule(qname, umi)``, so reads must not be
        pre-filtered in a way that drops whole molecules. In particular do **not** enable
        ``filters.improper_pair`` for consensus BAMs whose reads carry no PROPER_PAIR flag
        (some UMI-collapsed pipelines emit none) — every read would be discarded.
    """
    # Sub-model defaults, so a caller with no config never has to build a full
    # GbcmsDnaConfig (which requires variant/bam paths this entry point does not use).
    filters = config.filters if config is not None else ReadFilters()
    quality = config.quality if config is not None else QualityThresholds()
    align = config.alignment if config is not None else AlignmentConfig()
    threads = config.threads if config is not None else 1

    rs_variants = [_RsVariant(v.chrom, v.pos, v.ref, v.alt, v.variant_type.value) for v in variants]
    decomposed: list[_RsVariant | None] = [None] * len(rs_variants)
    status: list[str] | None = None

    # Normalize when a reference is available. Beyond left-alignment and ref_context, this
    # produces the *decomposed* form of complex indels — the form that often carries the ALT
    # support. Dropping it (as an earlier revision did) silently exported those molecules as
    # OTHER with zero ALT rows, so it is threaded through exactly as Pipeline does.
    # Variants are NOT filtered to PASS: `variant_index` is the caller's join key and must
    # stay positional. Failures are reported via `variant_status` instead.
    if reference_fasta is not None:
        prepared = prepare_variants(
            rs_variants,
            str(reference_fasta),
            quality.context_padding,
            is_maf,
            threads,
            quality.adaptive_context,
        )
        rs_variants = [p.variant for p in prepared]
        decomposed = [p.decomposed_variant for p in prepared]
        status = [p.gbcms_status for p in prepared]

    _counts, observations = count_bam_binned_observations(
        str(bam),
        rs_variants,
        decomposed,
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
        alignment_backend=align.backend,
        hmm_llr_threshold=align.hmm_llr_threshold,
        hmm_gap_open=align.hmm_gap_open,
        hmm_gap_extend=align.hmm_gap_extend,
        hmm_gap_open_repeat=align.hmm_gap_open_repeat,
        hmm_gap_extend_repeat=align.hmm_gap_extend_repeat,
        apply_baq=config.apply_baq if config is not None else False,
        umi_tag=config.umi_tag if config is not None else None,
        library_type=getattr(config, "library_type", "capture") if config else "capture",
        reference_fasta=str(reference_fasta) if reference_fasta is not None else None,
        observations_path=str(observations_path) if observations_path is not None else None,
    )
    if observations_path is not None:
        # Rows were written from Rust and never crossed the FFI boundary, so `observations`
        # comes back empty by design. Read the row count from the file rather than guessing.
        path = Path(observations_path)
        return ObservationResult(
            observations=[],
            path=path,
            n_rows=_parquet_row_count(path),
            variant_status=status,
        )
    return ObservationResult(
        observations=observations,
        path=None,
        n_rows=len(observations),
        variant_status=status,
    )
