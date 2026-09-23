"""
Pipeline Orchestrator: Manages the execution flow of gbcms.

This module handles:
1. Reading variants from input (VCF/MAF).
2. Preparing variants (MAF anchor, REF validation, left-alignment, ref_context)
   via the Rust ``prepare_variants()`` function.
3. Iterating over samples (BAM files).
4. Running the Rust-based counting engine for each sample.
5. Writing results to per-sample output files.
"""

import logging
import time
import types
from collections import Counter
from pathlib import Path
from typing import Any

from rich.console import Console
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeRemainingColumn,
)

from .core.kernel import CoordinateKernel
from .io.input import MafReader, VariantReader, VcfReader
from .io.output import MafWriter, VcfWriter
from .models.core import GbcmsBaseConfig, OutputFormat, Variant

_gbcms_rs = None


def _get_rs():
    """Lazy-import the Rust extension to avoid circular import with __init__.py."""
    global _gbcms_rs
    if _gbcms_rs is None:
        from gbcms import _rs

        _gbcms_rs = _rs
    return _gbcms_rs


logger = logging.getLogger(__name__)

__all__ = ["Pipeline"]


def _zero_counts():
    """Create a zero-count object mirroring BaseCounts for variants with no BAM coverage.

    All standard count fields default to 0; mFSD KS/LLR/delta/mean fields default to
    float('nan') since 0.0 would be scientifically misleading when the class is empty.
    Formatted as 'NA' in MAF output and '.' in VCF INFO via the _fmt/_fmt_vcf helpers.

    Note: ref_sizes/alt_sizes are NOT included — those are internal Rust fields written
    directly to Parquet by write_fsd_parquet() and never exposed to Python via PyO3.
    """
    _nan = float("nan")
    return types.SimpleNamespace(
        # Standard depth / allele counts
        dp=0,
        rd=0,
        ad=0,
        dp_fwd=0,
        rd_fwd=0,
        ad_fwd=0,
        dp_rev=0,
        rd_rev=0,
        ad_rev=0,
        dpf=0,
        rdf=0,
        adf=0,
        rdf_fwd=0,
        rdf_rev=0,
        adf_fwd=0,
        adf_rev=0,
        sb_pval=1.0,
        sb_or=1.0,
        fsb_pval=1.0,
        fsb_or=1.0,
        used_decomposed=False,
        # mFSD — raw counts
        mfsd_ref_count=0,
        mfsd_alt_count=0,
        mfsd_nonref_count=0,
        mfsd_n_count=0,
        # mFSD — mean sizes (NaN when class is empty)
        mfsd_ref_mean=_nan,
        mfsd_alt_mean=_nan,
        mfsd_nonref_mean=_nan,
        mfsd_n_mean=_nan,
        # mFSD — LLR (NaN when class is empty)
        mfsd_alt_llr=_nan,
        mfsd_ref_llr=_nan,
        # mFSD — pairwise KS triads (NaN when either class < MIN_FOR_KS=5)
        mfsd_delta_alt_ref=_nan,
        mfsd_ks_alt_ref=_nan,
        mfsd_pval_alt_ref=_nan,
        mfsd_delta_alt_nonref=_nan,
        mfsd_ks_alt_nonref=_nan,
        mfsd_pval_alt_nonref=_nan,
        mfsd_delta_ref_nonref=_nan,
        mfsd_ks_ref_nonref=_nan,
        mfsd_pval_ref_nonref=_nan,
        mfsd_delta_alt_n=_nan,
        mfsd_ks_alt_n=_nan,
        mfsd_pval_alt_n=_nan,
        mfsd_delta_ref_n=_nan,
        mfsd_ks_ref_n=_nan,
        mfsd_pval_ref_n=_nan,
        mfsd_delta_nonref_n=_nan,
        mfsd_ks_nonref_n=_nan,
        mfsd_pval_nonref_n=_nan,
        # Universal additions (both modes)
        mq0_count=0,
        alt_dist_end_median=_nan,
        ref_dist_end_median=_nan,
        singleton_alt_count=0,
        duplex_alt_count=0,
        # Decomposed ALT counting (invariant: any_alt = ad + partial_alt)
        any_alt=0,
        partial_alt=0,
        # N-base diagnostic (duplex masking QC)
        n_count=0,
        # RNA-specific (zeroed in DNA mode)
        sense_depth=0,
        antisense_depth=0,
        sense_strand_alt_count=0,
        antisense_strand_alt_count=0,
        rna_editing_site_overlap=False,
        splice_spanning_count=0,
        # mFSD — BH-corrected q-value and nucleosomal fractions (NaN like the
        # other continuous mFSD stats: 0.0 would read as a computed value)
        mfsd_qval_alt_ref=_nan,
        mfsd_sub_nuc_ref_frac=_nan,
        mfsd_sub_nuc_alt_frac=_nan,
        mfsd_sub_nuc_enrichment=_nan,
        mfsd_mono_nuc_ref_frac=_nan,
        mfsd_mono_nuc_alt_frac=_nan,
        # RNA+GTF annotation columns, mirroring what a counted variant the
        # annotator never touched would carry. ASJD p/q use the
        # not-applicable sentinel 1.0 (AsjdResult::empty): the writers render
        # values < 1.0 in scientific notation, so 0.0 would read as maximal
        # significance for a variant that was never counted.
        exon_boundary_dist=None,
        transcript_read_counts="",
        transcript_fragment_counts="",
        asjd_flag=False,
        asjd_pval=1.0,
        asjd_qval=1.0,
        asjd_ref_junction="",
        asjd_alt_junction="",
        asjd_ref_motif="",
        asjd_alt_motif="",
        asjd_ref_known=False,
        asjd_alt_known=False,
        asjd_n_ref_junc=0,
        asjd_n_alt_junc=0,
        asjd_n_ref_total=0,
        asjd_n_alt_total=0,
        asjd_diagnostic="",
    )


def read_variant_file(path: Path) -> list[Variant]:
    """Read raw variants from a ``.vcf``/``.vcf.gz``/``.vcf.bgz``/``.maf`` file.

    Format is selected by extension. This is the pre-normalization read, so no
    reference is required. Shared by the Pipeline (``_load_variants``) and the
    ``build-gtf-cache`` command, which only needs the variant chromosomes.
    """
    reader: VariantReader
    name_lower = path.name.lower()
    if (
        path.suffix.lower() == ".vcf"
        or name_lower.endswith(".vcf.gz")
        or name_lower.endswith(".vcf.bgz")
    ):
        reader = VcfReader(path)
    elif path.suffix.lower() == ".maf":
        reader = MafReader(path)
    else:
        raise ValueError(
            f"Unsupported variant file format: '{path.suffix}'. "
            "Expected .vcf, .vcf.gz, .vcf.bgz, or .maf."
        )
    variants = list(reader)
    if hasattr(reader, "close"):
        reader.close()
    return variants


def _is_mnp(ref_allele: str, alt_allele: str) -> bool:
    """Same-length multi-base substitution — what the engine dispatches to check_mnp."""
    return len(ref_allele) == len(alt_allele) > 1


def _mnp_discriminating_positions(variant: Any) -> list[tuple[int, str, str]]:
    """``(0-based pos, REF base, ALT base)`` for each MNP offset where REF and ALT differ."""
    return [
        (variant.pos + offset, ref_base, alt_base)
        for offset, (ref_base, alt_base) in enumerate(
            zip(variant.ref_allele, variant.alt_allele, strict=True)
        )
        if ref_base != alt_base
    ]


def _resolve_mnp_rescue(
    original_ad: int, component_ads: list[int | None]
) -> tuple[str, int | None]:
    """Decide a rescue candidate's outcome from its components' ALT counts.

    ``component_ads`` holds, per discriminating position, the synthetic SNV's
    ``ad`` — or ``None`` where that SNV failed REF validation and was not counted.

    Returns ``(outcome, index of the adopted component)``:
      - ``rescued``: the best component (highest ad; leftmost on ties) beats
        the haplotype's ``ad``.
      - ``no_improvement``: no counted component beats ``ad``. For consistent
        counts the candidate gate (``partial_alt > ad``) rules this out — every
        partial read matches ALT at an unmasked discriminating position, so the
        best component holds at least ``(partial_alt + ad) / 2`` reads — so the
        caller treats it as an anomaly.
      - ``ref_validation_failed``: no position could be counted.
    """
    counted = [(ad, idx) for idx, ad in enumerate(component_ads) if ad is not None]
    if not counted:
        return "ref_validation_failed", None
    best_ad, best_idx = max(counted, key=lambda pair: (pair[0], -pair[1]))
    if best_ad <= original_ad:
        return "no_improvement", None
    return "rescued", best_idx


def _format_rescue_audit(
    outcome: str,
    original: Any,
    positions: list[str] | None = None,
    adopted: str | None = None,
) -> str:
    """Build the ``gbcms_rescue`` value.

    Format: ``method=decomposed;outcome=<outcome>;original_ref=R;original_alt=A;
    original_partial=P[;adopted=<label>][;positions=<label>:<ad|ref_fail>,...]``
    where a label is ``chrom:pos(REF>ALT)`` (1-based). ``original_*`` are the
    MNP's own counts — for a rescued row the only record of its evaluation,
    because the count columns then carry the adopted component's counts.
    """
    parts = [
        "method=decomposed",
        f"outcome={outcome}",
        f"original_ref={original.rd}",
        f"original_alt={original.ad}",
        f"original_partial={original.partial_alt}",
    ]
    if adopted is not None:
        parts.append(f"adopted={adopted}")
    if positions:
        parts.append("positions=" + ",".join(positions))
    return ";".join(parts)


class Pipeline:
    """Main pipeline for processing BAM files and counting bases at variant positions."""

    def __init__(self, config: GbcmsBaseConfig):
        """
        Initialize the pipeline.

        Args:
            config: Configuration object with input/output paths and filter settings.
        """
        self.config = config
        self.console = Console()
        self._stats: dict[str, int | float] = {
            "samples_processed": 0,
            "total_variants": 0,
            "total_time": 0.0,
        }
        self._failed_samples: list[dict[str, str]] = []

    def run(self) -> dict:
        """
        Execute the pipeline.

        Returns:
            Dictionary with processing statistics.
        """
        start_time = time.perf_counter()
        logger.info("Starting gbcms pipeline")
        logger.info("Output directory: %s", self.config.output.directory)

        # Log all resolved parameters at DEBUG for full reproducibility (#19)
        logger.debug(
            "Parameters:\n"
            "  mode=%s\n"
            "  reference_fasta=%s\n"
            "  variant_file=%s\n"
            "  bam_files=%s\n"
            "  threads=%d\n"
            "  output_format=%s\n"
            "  output_suffix=%s\n"
            "  column_prefix=%s\n"
            "  min_mapq=%d\n"
            "  min_baseq=%d\n"
            "  fragment_qual_threshold=%d\n"
            "  context_padding=%d\n"
            "  adaptive_context=%s\n"
            "  alignment_backend=%s\n"
            "  apply_baq=%s\n"
            "  umi_tag=%s\n"
            "  show_normalization=%s\n"
            "  rescue_mnp=%s\n"
            "  mfsd=%s",
            self.config.mode,
            self.config.reference_fasta,
            self.config.variant_file,
            list(self.config.bam_files.keys()),
            self.config.threads,
            self.config.output.format.value,
            self.config.output.suffix or "(none)",
            self.config.output.column_prefix or "(none)",
            self.config.quality.min_mapping_quality,
            self.config.quality.min_base_quality,
            self.config.quality.fragment_qual_threshold,
            self.config.quality.context_padding,
            self.config.quality.adaptive_context,
            self.config.alignment.backend,
            self.config.apply_baq,
            self.config.umi_tag or "none",
            self.config.show_normalization,
            self.config.rescue_mnp,
            self.config.output.mfsd,
        )

        # 1. Load Variants (raw MAF/VCF coords)
        logger.debug("Loading variants from %s", self.config.variant_file)
        variants = self._load_variants()
        logger.info("Loaded %d variants", len(variants))

        if not variants:
            logger.error("No variants found. Exiting.")
            return self._stats

        # 2. Prepare variants: MAF anchor → validate REF → left-align → ref_context
        #    This replaces the old _validate_variants() + manual ref_context fetch.
        is_maf = self.config.variant_file.suffix.lower() == ".maf"
        rs_input = [
            _get_rs().Variant(
                v.chrom,
                v.pos,
                v.ref,
                v.alt,
                v.variant_type.value,
            )
            for v in variants
        ]
        prepared = _get_rs().prepare_variants(
            rs_input,
            str(self.config.reference_fasta),
            self.config.quality.context_padding,
            is_maf,
            self.config.threads,
            self.config.quality.adaptive_context,
        )

        # Split into valid (for counting) and all (for output)
        valid_indices = [i for i, p in enumerate(prepared) if p.gbcms_status == "PASS"]
        rs_variants = [prepared[i].variant for i in valid_indices]

        # Log validation results
        n_invalid = len(prepared) - len(valid_indices)
        logger.info(
            "Variant preparation: %d valid, %d rejected (%d total)",
            len(valid_indices),
            n_invalid,
            len(prepared),
        )
        invalid = [p for p in prepared if p.gbcms_status != "PASS"]
        for p in invalid[:5]:
            logger.warning(
                "Rejected variant: %s:%d %s>%s — %s (%s)",
                p.variant.chrom,
                p.original_pos + 1,
                p.original_ref,
                p.original_alt,
                p.gbcms_status,
                p.gbcms_status_reason or "no reason",
            )
        if len(invalid) > 5:
            logger.warning("... and %d more rejected variants", len(invalid) - 5)

        # Log variant type breakdown for transparency
        # MNPs (same-length multi-base substitutions) are classified as
        # COMPLEX by kernel.py but dispatched to check_mnp by the Rust
        # counting engine based on ref_len == alt_len.
        type_counts: dict[str, int] = {}
        mnp_count = 0
        for p in prepared:
            v = p.variant
            ref_len = len(v.ref_allele)
            alt_len = len(v.alt_allele)
            if ref_len == 1 and alt_len == 1:
                vtype = "SNP"
            elif _is_mnp(v.ref_allele, v.alt_allele):
                subtypes = {2: "DNP", 3: "TNP"}
                vtype = subtypes.get(ref_len, f"ONP({ref_len}bp)")
                mnp_count += 1
            elif ref_len > alt_len:
                vtype = "DEL"
            elif alt_len > ref_len:
                vtype = "INS"
            else:
                vtype = "COMPLEX"
            type_counts[vtype] = type_counts.get(vtype, 0) + 1
        type_str = ", ".join(f"{k}={v}" for k, v in sorted(type_counts.items()))
        logger.info("Variant types: %s", type_str)
        if mnp_count > 0:
            logger.info(
                "MNP counting: %d MNPs use selective discriminating-position "
                "quality gate (atomic block matching, no check_complex fallback)",
                mnp_count,
            )

        # Log normalization changes
        n_anchor = sum(1 for p in prepared if p.was_anchor_resolved)
        n_left = sum(1 for p in prepared if p.was_left_aligned)
        n_total = sum(1 for p in prepared if p.was_normalized)
        if n_total > 0:
            logger.info(
                "Normalized %d variants (%d anchor-resolved, %d left-aligned)",
                n_total,
                n_anchor,
                n_left,
            )

        if not rs_variants:
            # Every variant was rejected during preparation (verdict FAIL, e.g. a contig
            # mismatch → reason FETCH_FAILED, or EMPTY_ALLELE). This is NOT an empty variant
            # file (that returns earlier, before any variant exists) — the variants are real,
            # so we still fall through and write them per sample with their FAIL verdict + reason
            # in the `gbcms_status` column and zero counts, rather than silently emitting
            # no output. The run is not a failure (no sample raised); it exits 0.
            logger.warning(
                "No variants passed preparation (%d rejected); writing them with their "
                "FAIL_* status and zero counts so the reasons are in the output, not just the log.",
                len(prepared),
            )

        self._stats["total_variants"] = len(variants)
        self._stats["valid_variants"] = len(valid_indices)
        self._stats["mnp_variants"] = mnp_count

        # 3. Process Each Sample
        self.config.output.directory.mkdir(parents=True, exist_ok=True)
        samples = list(self.config.bam_files.items())

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TimeRemainingColumn(),
            console=self.console,
        ) as progress:
            task = progress.add_task("[cyan]Processing samples...", total=len(samples))

            for sample_name, bam_path in samples:
                progress.update(task, description=f"[cyan]Processing {sample_name}...")
                self._process_sample(
                    sample_name,
                    bam_path,
                    variants,
                    rs_variants,
                    prepared,
                    valid_indices,
                )
                progress.advance(task)

        # Calculate total time
        self._stats["total_time"] = time.perf_counter() - start_time

        # Log summary including failures
        if self._failed_samples:
            logger.error(
                "Pipeline completed with %d sample failure(s): %s",
                len(self._failed_samples),
                ", ".join(f"{s['name']} ({s['error']})" for s in self._failed_samples),
            )
        logger.info(
            "Pipeline completed: %d samples processed, %d failed, %.2fs",
            self._stats["samples_processed"],
            len(self._failed_samples),
            self._stats["total_time"],
        )

        # Include failed_samples in returned stats for callers
        return {**self._stats, "failed_samples": self._failed_samples}

    def _process_sample(
        self,
        sample_name: str,
        bam_path: Path,
        variants: list[Variant],
        rs_variants: list,
        prepared: list,
        valid_indices: list[int],
    ) -> None:
        """
        Process a single sample.

        Args:
            sample_name: Name of the sample.
            bam_path: Path to BAM file.
            variants: List of all input variants (for output).
            rs_variants: List of valid Rust variant objects (for counting).
            prepared: Full list of PreparedVariant objects.
            valid_indices: Indices of valid variants in the prepared list.
        """
        sample_start = time.perf_counter()
        logger.debug("Processing sample: %s (%s)", sample_name, bam_path)

        # Validate BAM Header
        if not self._validate_bam_header(bam_path, variants):
            logger.warning(
                "BAM %s may not contain variant chromosomes. Proceeding anyway.",
                sample_name,
            )

        try:
            # Run Rust Engine (only on valid variants)
            # Build decomposed variants list for dual-counting
            decomposed = [prepared[i].decomposed_variant for i in valid_indices]

            # Build sibling Variant objects for multi-allelic exclusion (Gap 1A)
            # For each variant in a multi-allelic group, collect the full Variant
            # objects of all OTHER variants in the same group. This allows the
            # Rust-side guard to run the complete classification pipeline
            # (check_allele_with_qual) for indels/complex/MNPs, not just SNPs.
            group_map: dict[int, list[int]] = {}
            for vi_pos, vi in enumerate(valid_indices):
                grp = prepared[vi].multi_allelic_group
                if grp is not None:
                    group_map.setdefault(grp, []).append(vi_pos)

            sibling_variants: list[list] = []
            for vi_pos, vi in enumerate(valid_indices):
                grp = prepared[vi].multi_allelic_group
                if grp is not None and grp in group_map:
                    siblings = [
                        prepared[valid_indices[j]].variant for j in group_map[grp] if j != vi_pos
                    ]
                    sibling_variants.append(siblings)
                else:
                    sibling_variants.append([])

            rust_start = time.perf_counter()
            if self.config.alignment.backend != "sw":
                logger.info("Using alignment backend: %s", self.config.alignment.backend)
            # One pass either way: when observations are requested, the same core also
            # writes the per-molecule rows straight to Parquet from Rust, so counts are
            # identical and the rows never cross the FFI boundary (flat memory at panel
            # scale). Mirrors how --mfsd-parquet delegates to the native writer.
            _obs_path: str | None = None
            if self.config.output.observations_parquet:
                _obs_path = str(
                    self.config.output.directory / f"{sample_name}.observations.parquet"
                )
            _count_fn = (
                _get_rs().count_bam_binned_observations if _obs_path else _get_rs().count_bam_binned
            )
            _result = _count_fn(
                str(bam_path),
                rs_variants,
                decomposed,
                sibling_variants=sibling_variants,
                **self._engine_kwargs(),
                **({"observations_path": _obs_path} if _obs_path else {}),
            )
            # The observations entry point returns (counts, rows); rows are empty because
            # they were written to Parquet. Counts are the same either way.
            counts_list = _result[0] if _obs_path else _result
            if _obs_path:
                logger.info("Wrote per-molecule observations → %s", _obs_path)
            rust_time = time.perf_counter() - rust_start
            logger.debug("Rust count_bam_binned completed in %.3fs", rust_time)

            # Append WARN_HOMOPOLYMER_DECOMP where the decomposed allele won.
            # Append (not overwrite) so a co-occurring WARN_REF_CORRECTED / MULTI_ALLELIC
            # survives; the verdict stays PASS. Reasons are '|'-separated.
            for idx, counts in zip(valid_indices, counts_list, strict=True):
                if counts.used_decomposed:
                    pv = prepared[idx]
                    if "WARN_HOMOPOLYMER_DECOMP" not in pv.gbcms_status_reason:
                        pv.gbcms_status_reason = (
                            f"{pv.gbcms_status_reason}|WARN_HOMOPOLYMER_DECOMP"
                            if pv.gbcms_status_reason
                            else "WARN_HOMOPOLYMER_DECOMP"
                        )

            # Merge counts back into full variant list
            # Valid variants get real counts; rejected variants get zero counts.
            full_counts = self._merge_counts(prepared, counts_list, valid_indices)

            # Post-counting: compute diagnostic flags (gbcms_diagnostic)
            self._compute_diagnostics(prepared, full_counts)

            # Post-counting: MNP rescue pass (optional, --rescue-mnp)
            if self.config.rescue_mnp:
                self._rescue_mnp_pass(prepared, full_counts, bam_path, sample_name)

            # Write Output (all variants, including rejected with zero counts)
            self._write_output(sample_name, variants, full_counts, prepared)
            self._stats["samples_processed"] += 1

            sample_time = time.perf_counter() - sample_start
            logger.debug("Sample %s completed in %.3fs", sample_name, sample_time)

        except Exception as e:
            # logger.exception captures the traceback; the stored message
            # names the exception class because str(e) alone can be a bare
            # dictionary key (KeyError) or even empty.
            logger.exception("Error processing sample %s", sample_name)
            self._failed_samples.append({"name": sample_name, "error": f"{type(e).__name__}: {e}"})

    def _engine_kwargs(self) -> dict[str, Any]:
        """Counting-engine settings shared by every ``count_bam_binned`` call for a sample.

        The main count and the MNP rescue's component re-count must classify
        reads identically, so both take their settings from here.
        """
        align_cfg = self.config.alignment
        return {
            "min_mapq": self.config.quality.min_mapping_quality,
            "min_baseq": self.config.quality.min_base_quality,
            "filter_duplicates": self.config.filters.duplicates,
            "filter_secondary": self.config.filters.secondary,
            "filter_supplementary": self.config.filters.supplementary,
            "filter_qc_failed": self.config.filters.qc_failed,
            "filter_improper_pair": self.config.filters.improper_pair,
            "filter_indel": self.config.filters.indel,
            "threads": self.config.threads,
            "fragment_qual_threshold": self.config.quality.fragment_qual_threshold,
            "alignment_backend": align_cfg.backend,
            "hmm_llr_threshold": align_cfg.hmm_llr_threshold,
            "hmm_gap_open": align_cfg.hmm_gap_open,
            "hmm_gap_extend": align_cfg.hmm_gap_extend,
            "hmm_gap_open_repeat": align_cfg.hmm_gap_open_repeat,
            "hmm_gap_extend_repeat": align_cfg.hmm_gap_extend_repeat,
            "apply_baq": self.config.apply_baq,
            "umi_tag": self.config.umi_tag,
            "mode": self.config.mode,
            "enforce_strandedness": getattr(self.config, "enforce_strandedness", False),
            "strandedness": getattr(self.config, "strandedness", "reverse"),
            "mfsd": self.config.output.mfsd,
            "rna_editing_db": (
                str(self.config.rna_editing_db)  # type: ignore[attr-defined]
                if getattr(self.config, "rna_editing_db", None)
                else None
            ),
            "gtf_path": (
                str(self.config.gtf)  # type: ignore[attr-defined]
                if getattr(self.config, "gtf", None)
                else None
            ),
            "gtf_cache_dir": (
                str(self.config.gtf_cache_dir)  # type: ignore[attr-defined]
                if getattr(self.config, "gtf_cache_dir", None)
                else None
            ),
            "reference_fasta": str(self.config.reference_fasta),
            "library_type": getattr(self.config, "library_type", "capture"),
        }

    @staticmethod
    def _merge_counts(
        prepared: list,
        counts_list: list,
        valid_indices: list[int],
    ) -> list:
        """Merge real counts for valid variants with zero counts for rejected ones.

        Returns a list with one BaseCounts per input variant (same order as prepared).
        """
        counts_by_idx: dict[int, object] = {}
        for offset, vi in enumerate(valid_indices):
            counts_by_idx[vi] = counts_list[offset]

        merged = []
        for i, _pv in enumerate(prepared):
            if i in counts_by_idx:
                merged.append(counts_by_idx[i])
            else:
                merged.append(_zero_counts())
        return merged

    def _compute_diagnostics(self, prepared: list, full_counts: list) -> None:
        """Populate gbcms_diagnostic (semicolon-separated) for every PASS variant.

        FAIL variants keep gbcms_diagnostic empty. Logs the per-sample flag
        distribution. Flag definitions: :meth:`_diagnostic_flags`.
        """
        flag_counts: dict[str, int] = {}

        for pv, counts in zip(prepared, full_counts, strict=True):
            if pv.gbcms_status != "PASS":
                continue
            flags = self._diagnostic_flags(pv.variant, counts)
            pv.gbcms_diagnostic = ";".join(flags)
            for flag_name in flags:
                # Parametric flags (e.g. MNP_DISC_RATIO(2/5)) count under their base name
                base_flag = flag_name.split("(")[0]
                flag_counts[base_flag] = flag_counts.get(base_flag, 0) + 1

        if flag_counts:
            summary = ", ".join(f"{flag}={count}" for flag, count in sorted(flag_counts.items()))
            logger.info("Diagnostic flags: %s", summary)
        else:
            logger.debug("No diagnostic flags triggered")

    def _diagnostic_flags(self, variant: Any, counts: Any) -> list[str]:
        """Diagnostic flags for one counted PASS variant, in output order.

        Flags:
            ZERO_ALT: ad == 0 and variant was successfully counted.
            PARTIAL_DOMINANT: partial_alt > ad (more structural evidence
                than confirmed ALT calls). For pure indels partial_alt
                includes wrong-length evidence — reads whose CIGAR proves an
                indel of a DIFFERENT length in the same tract (a distinct
                slippage allele; at >=50bp deletion loci, a real different
                large event). A dominant partial count therefore usually
                means the locus carries a coexisting allele the annotation
                does not describe; per-read lengths are visible with --trace.
            MNP_DISC_RATIO(n/m): for MNPs (ref_len == alt_len > 1),
                always emitted showing discriminating position ratio.
            MNP_RESCUE_ELIGIBLE: disc/len ≤ rescue_mnp_threshold.
            HIGH_N_FRACTION(f): n_count / dp > 0.05 (duplex masking hotspot).
            NON_DISCRIMINATING_LOCUS: a sibling combination reconstructs the
                reference haplotype, so REF and ALT are sequence-indistinguishable
                and reads tie to NEITHER (explains a zeroed RD/AD at a covered locus).
            SPLICE_SKIP_DOMINANT(n): deletion-type locus where more reads
                asserted splicing over the deleted span (CIGAR N, excluded from
                DP as no-observation) than confirmed ALT. RNA aligners represent
                large deletions as splices (STAR: ≥ alignIntronMin, default
                21bp), so AD=0 here may mean carriers exist as junction reads.
        """
        flags: list[str] = []

        # ZERO_ALT: no confirmed ALT reads despite successful counting
        if counts.ad == 0:
            flags.append("ZERO_ALT")

        # PARTIAL_DOMINANT: more structural/partial evidence than confirmed ALT
        if counts.partial_alt > counts.ad:
            flags.append("PARTIAL_DOMINANT")

        # MNP_DISC_RATIO: for MNPs, always emit discriminating position ratio
        # as a diagnostic signal. Additionally, mark rescue eligibility
        # based on the configurable --rescue-mnp-threshold.
        ref_allele = variant.ref_allele
        alt_allele = variant.alt_allele
        if _is_mnp(ref_allele, alt_allele):
            disc = len(_mnp_discriminating_positions(variant))
            flags.append(f"MNP_DISC_RATIO({disc}/{len(ref_allele)})")
            if disc / len(ref_allele) <= self.config.rescue_mnp_threshold:
                flags.append("MNP_RESCUE_ELIGIBLE")

        # HIGH_N_FRACTION: high rate of N-bases at discriminating positions
        if counts.dp > 0 and counts.n_count / counts.dp > 0.05:
            frac = counts.n_count / counts.dp
            flags.append(f"HIGH_N_FRACTION({frac:.2f})")

        # NON_DISCRIMINATING_LOCUS: a sibling combination reconstructs the
        # reference haplotype, so REF and ALT are sequence-indistinguishable and
        # reads tie to NEITHER — surfaces an otherwise-silent zeroed RD/AD.
        if getattr(counts, "non_discriminating_locus", False):
            flags.append("NON_DISCRIMINATING_LOCUS")

        # SPLICE_SKIP_DOMINANT: at a deletion-type locus, more reads
        # asserted splicing over the deleted span (CIGAR N — excluded
        # from DP as no-observation) than confirmed ALT. RNA aligners
        # write large deletions as splices (STAR: any deletion ≥
        # alignIntronMin, default 21bp, becomes N), so a zeroed AD here
        # can mean the carriers exist but are represented as junctions —
        # inspect the locus in IGV before trusting AD=0.
        excluded = getattr(counts, "splice_skip_excluded", 0)
        if len(ref_allele) > len(alt_allele) and excluded > counts.ad:
            flags.append(f"SPLICE_SKIP_DOMINANT({excluded})")

        return flags

    def _rescue_mnp_pass(
        self,
        prepared: list,
        full_counts: list,
        bam_path: Path,
        sample_name: str,
    ) -> None:
        """Report the best-supported component of MNPs whose haplotype is dominated.

        Rescue is for MNPs whose carriers hold only a component of the annotated
        haplotype (e.g. SNVs on different molecules annotated as one MNP): the
        engine correctly reports the haplotype as absent and those carriers as
        ``partial_alt``, and rescue reports the best-supported component instead.

        Candidates are PASS MNPs flagged ``MNP_RESCUE_ELIGIBLE`` with
        ``partial_alt > ad``. The gate is dominance, not ``ad == 0``: a component
        carrier whose other discriminating base is masked counts as full ALT, and
        one such read must not block rescue.

        Each candidate's discriminating positions are re-counted as synthetic
        SNVs with the sample's own settings. When the best component beats the
        MNP's ``ad``, its BaseCounts replace the row's wholesale — every count,
        fragment, strand, strand-bias, mFSD and RNA column then comes from one
        counting pass, so the counting invariants hold on the written row — and
        the row's diagnostics are recomputed from them. The MNP's own counts
        survive in ``gbcms_rescue`` (format: :func:`_format_rescue_audit`;
        outcomes: :func:`_resolve_mnp_rescue`).

        Grouped MNPs are skipped (``outcome=skipped_grouped``): their reads are
        exclusively assigned against co-annotated siblings, and a sibling-free
        component re-count would hand contested reads back.

        ``gbcms_rescue`` is reset for every variant first: ``prepared`` is shared
        by all samples of a run, so a value set for an earlier sample would
        otherwise be written on this sample's rows.
        """
        for pv in prepared:
            pv.gbcms_rescue = ""

        rescue_start = time.perf_counter()
        outcomes: Counter[str] = Counter()
        candidates: list[tuple[int, list[tuple[int, str, str]]]] = []
        for i, (pv, counts) in enumerate(zip(prepared, full_counts, strict=True)):
            if (
                pv.gbcms_status != "PASS"
                or "MNP_RESCUE_ELIGIBLE" not in pv.gbcms_diagnostic.split(";")
                or counts.partial_alt <= counts.ad
            ):
                continue
            v = pv.variant
            if pv.multi_allelic_group is not None:
                pv.gbcms_rescue = _format_rescue_audit("skipped_grouped", counts)
                outcomes["skipped_grouped"] += 1
                logger.debug(
                    "MNP rescue: %s:%d %s>%s skipped — co-annotated group %d owns its reads",
                    v.chrom,
                    v.pos + 1,
                    v.ref_allele,
                    v.alt_allele,
                    pv.multi_allelic_group,
                )
                continue
            candidates.append((i, _mnp_discriminating_positions(v)))

        component_counts = (
            self._count_mnp_components(prepared, candidates, bam_path) if candidates else []
        )
        for (i, positions), components in zip(candidates, component_counts, strict=True):
            pv, original = prepared[i], full_counts[i]
            v = pv.variant
            labels = [f"{v.chrom}:{pos + 1}({ref}>{alt})" for pos, ref, alt in positions]
            entries = [
                f"{label}:{'ref_fail' if c is None else c.ad}"
                for label, c in zip(labels, components, strict=True)
            ]
            outcome, best = _resolve_mnp_rescue(
                original.ad, [None if c is None else c.ad for c in components]
            )
            outcomes[outcome] += 1
            if best is None:
                pv.gbcms_rescue = _format_rescue_audit(outcome, original, entries)
                logger.warning(
                    "MNP rescue: %s:%d %s>%s not rescued (%s) — MNP ad=%d partial_alt=%d, "
                    "components %s; counts left as the MNP evaluation",
                    v.chrom,
                    v.pos + 1,
                    v.ref_allele,
                    v.alt_allele,
                    outcome,
                    original.ad,
                    original.partial_alt,
                    ",".join(entries),
                )
                continue
            adopted = components[best]
            assert adopted is not None, "_resolve_mnp_rescue adopts only counted components"
            full_counts[i] = adopted
            pv.gbcms_diagnostic = ";".join(self._diagnostic_flags(v, adopted))
            pv.gbcms_rescue = _format_rescue_audit(outcome, original, entries, labels[best])
            logger.debug(
                "MNP rescue: %s:%d %s>%s rescued — adopted %s (ad %d → %d, rd %d → %d)",
                v.chrom,
                v.pos + 1,
                v.ref_allele,
                v.alt_allele,
                labels[best],
                original.ad,
                adopted.ad,
                original.rd,
                adopted.rd,
            )

        if not outcomes:
            logger.debug("MNP rescue for %s: no candidates", sample_name)
            return
        logger.info(
            "MNP rescue for %s: %s (%.3fs)",
            sample_name,
            ", ".join(f"{name}={n}" for name, n in sorted(outcomes.items())),
            time.perf_counter() - rescue_start,
        )
        if outcomes["rescued"] and self.config.output.observations_parquet:
            logger.info(
                "Observations Parquet for %s records the MNP evaluation of %d rescued row(s); "
                "their written counts are the adopted component's (see gbcms_rescue)",
                sample_name,
                outcomes["rescued"],
            )

    def _count_mnp_components(
        self,
        prepared: list,
        candidates: list[tuple[int, list[tuple[int, str, str]]]],
        bam_path: Path,
    ) -> list[list[Any | None]]:
        """Count every candidate's discriminating positions as synthetic SNVs.

        One batched prepare + count call. Returns, per candidate, one entry per
        position: the SNV's BaseCounts, or ``None`` where the synthetic SNV
        failed preparation. The MNP itself passed REF validation, so a failing
        component is an inconsistency — each one is logged as a warning.
        """
        rs = _get_rs()
        snvs = [
            rs.Variant(prepared[i].variant.chrom, pos, ref, alt, "SNP")
            for i, positions in candidates
            for pos, ref, alt in positions
        ]
        snv_prepared = rs.prepare_variants(
            snvs,
            str(self.config.reference_fasta),
            self.config.quality.context_padding,
            False,  # is_maf: synthetic 0-based coordinates, no MAF anchor handling
            self.config.threads,
            self.config.quality.adaptive_context,
        )
        valid = [sp for sp in snv_prepared if sp.gbcms_status == "PASS"]
        for sp in snv_prepared:
            if sp.gbcms_status != "PASS":
                logger.warning(
                    "MNP rescue: synthetic SNV %s:%d %s>%s failed preparation (%s) — "
                    "reported as ref_fail",
                    sp.variant.chrom,
                    sp.variant.pos + 1,
                    sp.variant.ref_allele,
                    sp.variant.alt_allele,
                    sp.gbcms_status_reason or sp.gbcms_status,
                )

        counted = iter(
            rs.count_bam_binned(
                str(bam_path),
                [sp.variant for sp in valid],
                [None] * len(valid),
                sibling_variants=[[] for _ in valid],
                **self._engine_kwargs(),
            )
            if valid
            else []
        )
        flat = [next(counted) if sp.gbcms_status == "PASS" else None for sp in snv_prepared]

        per_candidate: list[list[Any | None]] = []
        offset = 0
        for _, positions in candidates:
            per_candidate.append(flat[offset : offset + len(positions)])
            offset += len(positions)
        return per_candidate

    def _load_variants(self) -> list[Variant]:
        """Load variants based on file extension.

        Delegates to the module-level :func:`read_variant_file`. The CLI pre-checks
        the extension at parse time (before Pydantic); an unsupported extension here
        means Pipeline was called programmatically — :func:`read_variant_file` raises
        ValueError as a defensive backstop.
        """
        return read_variant_file(self.config.variant_file)

    def _validate_bam_header(self, bam_path: Path, variants: list[Variant]) -> bool:
        """Check if BAM/CRAM header contains chromosomes from variants.

        Uses pysam auto-detect mode (no explicit format flag) so both BAM and
        CRAM files are handled transparently.  For CRAM files, ``reference_filename``
        is passed for correct header decoding.
        """
        try:
            import pysam

            # Auto-detect format (BAM/CRAM) — reference_filename is required for CRAM
            # decoding but harmless for BAM.
            with pysam.AlignmentFile(
                str(bam_path),
                reference_filename=str(self.config.reference_fasta),
            ) as bam:
                bam_chroms = set(bam.references)

            norm_bam_chroms = {CoordinateKernel.normalize_chromosome(c) for c in bam_chroms}

            if variants:
                v = variants[0]
                norm_v_chrom = CoordinateKernel.normalize_chromosome(v.chrom)
                if norm_v_chrom not in norm_bam_chroms:
                    return False
            return True
        except Exception as e:
            logger.warning("Could not validate BAM/CRAM header: %s", e)
            return True

    def _load_contigs_from_fai(self) -> list[tuple[str, int]]:
        """Load contig names and lengths from the FASTA index (.fai).

        Returns a list of (name, length) tuples for VCF ##contig headers.
        Falls back to an empty list if the FAI is missing or malformed —
        the pipeline should not fail just because of missing contig headers.
        """
        fai_path = Path(str(self.config.reference_fasta) + ".fai")
        if not fai_path.exists():
            logger.warning(
                "FASTA index not found: %s — VCF ##contig headers will be omitted",
                fai_path,
            )
            return []

        contigs: list[tuple[str, int]] = []
        try:
            with open(fai_path) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        contigs.append((parts[0], int(parts[1])))
            logger.debug("Loaded %d contigs from %s", len(contigs), fai_path.name)
        except (ValueError, OSError) as e:
            logger.warning(
                "Failed to parse FASTA index %s: %s — VCF ##contig headers will be omitted",
                fai_path,
                e,
            )
            return []
        return contigs

    def _write_output(
        self,
        sample_name: str,
        variants: list[Variant],
        counts_list: list,
        prepared: list | None = None,
    ) -> None:
        """Write results to output file.

        Args:
            sample_name: Name of the sample.
            variants: Original input variants.
            counts_list: Merged counts (one per variant, including zero-count stubs).
            prepared: PreparedVariant objects for normalization/validation info.
        """
        ext = "vcf" if self.config.output.format == OutputFormat.VCF else "maf"
        suffix = self.config.output.suffix
        output_path = self.config.output.directory / f"{sample_name}{suffix}.{ext}"

        writer: VcfWriter | MafWriter
        if self.config.output.format == OutputFormat.VCF:
            # Load contigs from FAI for VCF ##contig headers (lazy — only for VCF output)
            contigs = self._load_contigs_from_fai()

            # mode= is required so RNA-specific INFO/FORMAT headers and data fields
            # (SEN, ANT, ASEN, RED, SPL) are included when self.config.mode == "rna".
            # Without it, VcfWriter defaults to mode="dna" and RNA columns are silently absent.
            writer = VcfWriter(
                output_path,
                sample_name=sample_name,
                show_normalization=self.config.show_normalization,
                mfsd=self.config.output.mfsd,
                mode=self.config.mode,
                rescue_mnp=self.config.rescue_mnp,
                has_gtf=bool(getattr(self.config, "gtf", None)),
                command_line=self.config.command_line,
                reference_fasta=str(self.config.reference_fasta),
                contigs=contigs,
            )
        else:
            # mode= is required so RNA-specific MAF columns (rna_sense_depth, etc.)
            # are appended when self.config.mode == "rna".
            # Without it, MafWriter defaults to mode="dna" and RNA columns are silently absent.
            writer = MafWriter(
                output_path,
                column_prefix=self.config.output.column_prefix,
                preserve_barcode=self.config.output.preserve_barcode,
                show_normalization=self.config.show_normalization,
                mfsd=self.config.output.mfsd,
                mode=self.config.mode,
                rescue_mnp=self.config.rescue_mnp,
                has_gtf=bool(getattr(self.config, "gtf", None)),
                command_line=self.config.command_line,
            )
        logger.debug(
            "Writer initialised: format=%s, mode=%s, sample=%s, path=%s",
            self.config.output.format.value,
            self.config.mode,
            sample_name,
            output_path,
        )

        for i, (v, counts) in enumerate(zip(variants, counts_list, strict=True)):
            pv = prepared[i] if prepared else None

            # Build norm_variant only when normalization display is enabled
            norm_v = None
            if pv and pv.was_normalized:
                norm_v = Variant(
                    chrom=pv.variant.chrom,
                    pos=pv.variant.pos,
                    ref=pv.variant.ref_allele,
                    alt=pv.variant.alt_allele,
                    variant_type=v.variant_type,
                )

            writer.write(
                v,
                counts,
                sample_name=sample_name,
                gbcms_status=pv.gbcms_status if pv else "PASS",
                gbcms_status_reason=pv.gbcms_status_reason if pv else "",
                gbcms_diagnostic=pv.gbcms_diagnostic if pv else "",
                gbcms_rescue=pv.gbcms_rescue if pv else "",
                norm_variant=norm_v,
            )

        writer.close()
        logger.debug("Results written to %s", output_path)

        # Write companion mFSD Parquet when --mfsd-parquet is enabled.
        # Delegates to the native Rust writer (no pyarrow dep). Rejected
        # variants carry Python zero-count stubs that PyO3 cannot cast as
        # BaseCounts, and they have no fragment sizes to write anyway — pass
        # only the counted (real BaseCounts) rows and say how many were left
        # out, so a row-count difference vs the MAF/VCF is explained.
        if self.config.output.mfsd_parquet:
            fsd_path = output_path.with_suffix("").with_suffix(".fsd.parquet")
            counted = [
                (v, c)
                for v, c in zip(variants, counts_list, strict=True)
                if not isinstance(c, types.SimpleNamespace)
            ]
            excluded = len(variants) - len(counted)
            _get_rs().write_fsd_parquet(
                str(fsd_path),
                [v.chrom for v, _ in counted],
                [v.pos + 1 for v, _ in counted],  # 1-based MAF/VCF convention
                [v.ref for v, _ in counted],
                [v.alt for v, _ in counted],
                [c for _, c in counted],
            )
            if excluded:
                logger.info(
                    "mFSD Parquet written: %s (%d variants; %d rejected "
                    "variant(s) excluded — no fragment data exists for them)",
                    fsd_path,
                    len(counted),
                    excluded,
                )
            else:
                logger.info(
                    "mFSD Parquet written: %s (%d variants)",
                    fsd_path,
                    len(counted),
                )

            # Generate mFSD HTML report when --mfsd-report is enabled.
            # Runs after parquet write since it reads the parquet file.
            if self.config.output.mfsd_report:
                try:
                    from .report import generate_mfsd_report

                    report_path = output_path.with_suffix("").with_suffix(".mfsd_report.html")
                    generate_mfsd_report(
                        parquet_path=fsd_path,
                        maf_path=output_path,
                        output_path=report_path,
                        min_alt=self.config.output.mfsd_report_min_alt,
                        max_variants=self.config.output.mfsd_report_max_variants,
                        sample_name=sample_name,
                    )
                except Exception:
                    logger.exception(
                        "mFSD report generation failed (non-fatal); " "main output is unaffected."
                    )
        elif self.config.output.mfsd:
            logger.debug(
                "mFSD analysis enabled but --mfsd-parquet not set; "
                "raw fragment size arrays not written to disk."
            )
