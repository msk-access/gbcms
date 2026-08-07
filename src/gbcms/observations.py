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
        ...  # obs.variant_index, obs.molecule_hash, obs.allele, obs.best_qual, obs.min_mapq
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

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


class _Unset:
    """Sentinel distinguishing "argument not supplied" from an explicit ``None``.

    Needed only for ``umi_tag``, which is genuinely nullable: ``None`` is a *meaningful*
    value there (group by read pair, no UMI). Without this, a caller passing a ``config``
    that sets a UMI tag plus ``umi_tag=None`` to override it would silently inherit the
    config's tag and group by UMI family anyway — changing what counts as one molecule
    with nothing to indicate it.
    """

    def __repr__(self) -> str:  # pragma: no cover - debugging aid only
        return "<unset>"


_UNSET = _Unset()

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
    filters: ReadFilters | None = None,
    quality: QualityThresholds | None = None,
    alignment: AlignmentConfig | None = None,
    umi_tag: str | None = cast(Any, _UNSET),
    threads: int | None = None,
    apply_baq: bool | None = None,
    library_type: str | None = None,
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
        config: A full :class:`GbcmsDnaConfig` to take settings from. Convenient when you
            already have one; prefer the individual arguments below otherwise, since
            ``GbcmsDnaConfig`` requires ``variant_file``/``bam_files``/``reference_fasta``/
            ``output`` — none of which this entry point reads.
        filters: Read filters. **Set ``supplementary=False`` for cross-locus phasing**: a
            molecule spanning a large deletion reaches the far locus only through its
            supplementary alignment, and the default (filtered, matching counting) leaves
            that locus with no observation at all. Such records join their primary's
            fragment via the QNAME hash, so admitting them cannot double-count; they still
            never contribute to read-level ``dp``/``rd``/``ad``, which is why an admitted
            supplementary can raise ``dpf`` while leaving ``dp`` unchanged.
        quality: MAPQ/base-quality gates. Reads below them never reach fragment evidence,
            so they bound ``Observation.min_mapq`` from below.
        alignment: Alignment backend and PairHMM parameters. Defaults to ``pairhmm``.
        umi_tag: BAM tag holding the UMI (e.g. ``"MI"``). Folded into ``molecule_hash``, so
            it decides what counts as one molecule. Passing ``None`` is an **explicit**
            choice of qname-keyed (per-read-pair) grouping and overrides ``config`` — unlike
            the other arguments, omitting it and passing ``None`` are different things.
        threads: Worker threads. Defaults to 1.
        apply_baq: Apply BAQ before classification. Defaults to off.
        library_type: ``"capture"`` (default) or ``"amplicon"``. In amplicon mode the read
            number is folded into ``molecule_hash``, making the key per-read-end.
        observations_path: When given, rows are written to Parquet from Rust rather than
            materialized as Python objects — for panel- or genome-wide runs, where the FFI
            round-trip, not the counting, is the bottleneck. The returned
            ``observations`` list is then empty and ``path``/``n_rows`` describe the file,
            whose columns are ``variant_index, chrom, pos, ref, alt, molecule_hash, allele,
            best_qual, min_mapq`` (self-describing, so it stands alone once written).

    Returns:
        An :class:`ObservationResult`. Check ``variant_status`` before trusting rows for a
        given ``variant_index``.

    Note:
        An individual argument **overrides** the corresponding field of ``config`` when both
        are given, so a caller can reuse a pipeline config and adjust one setting. For every
        argument except ``umi_tag``, passing ``None`` means "not supplied" and defers to
        ``config``; ``umi_tag`` treats ``None`` as an explicit override, because ``None`` is
        a real grouping choice there rather than an absence.

    Note:
        Molecules are grouped by ``hash_molecule(qname, umi)``, so reads must not be
        pre-filtered in a way that drops whole molecules. In particular do **not** enable
        ``filters.improper_pair`` for consensus BAMs whose reads carry no PROPER_PAIR flag
        (some UMI-collapsed pipelines emit none) — every read would be discarded.

    Note:
        Filter and quality settings apply to the **counting**, not to the export layered on
        top: the returned rows and the counts they reconcile with (``adf``/``rdf``/``dpf``)
        come from one pass over one read set. That is why these default to the counting
        defaults rather than to something phasing-friendly — a different default here would
        make rows and counts disagree about which reads exist.
    """

    # Resolution order: explicit argument > `config` field > library default. `config` is
    # accepted for convenience, but its four required fields (variant_file, bam_files,
    # reference_fasta, output) are ones this entry point never reads — so requiring it just
    # to change a filter meant fabricating paths that look meaningful and are not.
    def _pick(explicit: object, attr: str, fallback: object, unset: object = None) -> object:
        """explicit argument > `config` field > library default.

        `unset` names the value that means "caller did not supply this". It is `None` for
        every parameter whose type excludes `None`, and `_UNSET` for `umi_tag`, where `None`
        is itself a meaningful choice (group by read pair rather than UMI family).
        """
        if explicit is not unset:
            return explicit
        if config is not None:
            return getattr(config, attr, fallback)
        return fallback

    filters = cast(ReadFilters, _pick(filters, "filters", ReadFilters()))
    quality = cast(QualityThresholds, _pick(quality, "quality", QualityThresholds()))
    align = cast(AlignmentConfig, _pick(alignment, "alignment", AlignmentConfig()))
    threads = cast(int, _pick(threads, "threads", 1))
    apply_baq = cast(bool, _pick(apply_baq, "apply_baq", False))
    umi_tag = cast("str | None", _pick(umi_tag, "umi_tag", None, unset=_UNSET))
    library_type = cast(str, _pick(library_type, "library_type", "capture"))

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
        apply_baq=apply_baq,
        umi_tag=umi_tag,
        library_type=library_type,
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
