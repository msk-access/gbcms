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
from pathlib import Path

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
        # RNA-specific (zeroed in DNA mode)
        sense_depth=0,
        antisense_depth=0,
        sense_strand_alt_count=0,
        antisense_strand_alt_count=0,
        rna_editing_site_overlap=False,
        splice_spanning_count=0,
    )


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
        valid_indices = [
            i for i, p in enumerate(prepared) if p.validation_status.startswith("PASS")
        ]
        rs_variants = [prepared[i].variant for i in valid_indices]

        # Log validation results
        n_invalid = len(prepared) - len(valid_indices)
        logger.info(
            "Variant preparation: %d valid, %d rejected (%d total)",
            len(valid_indices),
            n_invalid,
            len(prepared),
        )
        invalid = [p for p in prepared if not p.validation_status.startswith("PASS")]
        for p in invalid[:5]:
            logger.warning(
                "Rejected variant: %s:%d %s>%s — %s",
                p.variant.chrom,
                p.original_pos + 1,
                p.original_ref,
                p.original_alt,
                p.validation_status,
            )
        if len(invalid) > 5:
            logger.warning("... and %d more rejected variants", len(invalid) - 5)

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
            logger.error("No valid variants remaining after preparation. Exiting.")
            return self._stats

        self._stats["total_variants"] = len(variants)
        self._stats["valid_variants"] = len(valid_indices)

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
            align_cfg = self.config.alignment
            if align_cfg.backend != "sw":
                logger.info("Using alignment backend: %s", align_cfg.backend)
            counts_list = _get_rs().count_bam_binned(
                str(bam_path),
                rs_variants,
                decomposed,
                min_mapq=self.config.quality.min_mapping_quality,
                min_baseq=self.config.quality.min_base_quality,
                filter_duplicates=self.config.filters.duplicates,
                filter_secondary=self.config.filters.secondary,
                filter_supplementary=self.config.filters.supplementary,
                filter_qc_failed=self.config.filters.qc_failed,
                filter_improper_pair=self.config.filters.improper_pair,
                filter_indel=self.config.filters.indel,
                threads=self.config.threads,
                fragment_qual_threshold=self.config.quality.fragment_qual_threshold,
                sibling_variants=sibling_variants,
                alignment_backend=align_cfg.backend,
                hmm_llr_threshold=align_cfg.hmm_llr_threshold,
                hmm_gap_open=align_cfg.hmm_gap_open,
                hmm_gap_extend=align_cfg.hmm_gap_extend,
                hmm_gap_open_repeat=align_cfg.hmm_gap_open_repeat,
                hmm_gap_extend_repeat=align_cfg.hmm_gap_extend_repeat,
                apply_baq=self.config.apply_baq,
                umi_tag=self.config.umi_tag,
                mode=self.config.mode,
                enforce_strandedness=getattr(self.config, "enforce_strandedness", False),
                rna_editing_db=(
                    str(self.config.rna_editing_db)  # type: ignore[attr-defined]
                    if getattr(self.config, "rna_editing_db", None)
                    else None
                ),
            )
            rust_time = time.perf_counter() - rust_start
            logger.debug("Rust count_bam_binned completed in %.3fs", rust_time)

            # Update validation_status for variants where decomposed allele won
            for idx, counts in zip(valid_indices, counts_list, strict=True):
                if counts.used_decomposed:
                    prepared[idx].validation_status = "PASS_WARN_HOMOPOLYMER_DECOMP"

            # Merge counts back into full variant list
            # Valid variants get real counts; rejected variants get zero counts.
            full_counts = self._merge_counts(prepared, counts_list, valid_indices)

            # Write Output (all variants, including rejected with zero counts)
            self._write_output(sample_name, variants, full_counts, prepared)
            self._stats["samples_processed"] += 1

            sample_time = time.perf_counter() - sample_start
            logger.debug("Sample %s completed in %.3fs", sample_name, sample_time)

        except Exception as e:
            logger.error("Error processing sample %s: %s", sample_name, e)
            self._failed_samples.append({"name": sample_name, "error": str(e)})

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

    def _load_variants(self) -> list[Variant]:
        """Load variants based on file extension."""
        path = self.config.variant_file
        reader: VariantReader

        # Note: The CLI pre-checks the extension at parse time (before Pydantic), so
        # reaching this branch with an unsupported extension means Pipeline was called
        # programmatically rather than via the CLI.  We raise ValueError here as a
        # defensive backstop.
        #
        # Accepted compressed formats:
        #   .vcf.gz  — gzip/bgzip (standard tabix format)
        #   .vcf.bgz — explicit bgzip extension used by some pipelines
        # Both are block-gzip compatible and handled identically by pysam.
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

    def _validate_bam_header(self, bam_path: Path, variants: list[Variant]) -> bool:
        """Check if BAM header contains chromosomes from variants."""
        try:
            import pysam

            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                bam_chroms = set(bam.references)

            norm_bam_chroms = {CoordinateKernel.normalize_chromosome(c) for c in bam_chroms}

            if variants:
                v = variants[0]
                norm_v_chrom = CoordinateKernel.normalize_chromosome(v.chrom)
                if norm_v_chrom not in norm_bam_chroms:
                    return False
            return True
        except Exception as e:
            logger.warning("Could not validate BAM header: %s", e)
            return True

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
            # mode= is required so RNA-specific INFO/FORMAT headers and data fields
            # (SEN, ANT, ASEN, RED, SPL) are included when self.config.mode == "rna".
            # Without it, VcfWriter defaults to mode="dna" and RNA columns are silently absent.
            writer = VcfWriter(
                output_path,
                sample_name=sample_name,
                show_normalization=self.config.show_normalization,
                mfsd=self.config.output.mfsd,
                mode=self.config.mode,
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
                validation_status=pv.validation_status if pv else "PASS",
                norm_variant=norm_v,
            )

        writer.close()
        logger.debug("Results written to %s", output_path)

        # Write companion mFSD Parquet when --mfsd-parquet is enabled.
        # Delegates to the native Rust writer (no pyarrow dep).
        if self.config.output.mfsd_parquet:
            fsd_path = output_path.with_suffix("").with_suffix(".fsd.parquet")
            _get_rs().write_fsd_parquet(
                str(fsd_path),
                [v.chrom for v in variants],
                [v.pos + 1 for v in variants],  # 1-based MAF/VCF convention
                [v.ref for v in variants],
                [v.alt for v in variants],
                counts_list,
            )
            logger.info(
                "mFSD Parquet written: %s (%d variants)",
                fsd_path,
                len(variants),
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
