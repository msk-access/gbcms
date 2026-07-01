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
        valid_indices = [i for i, p in enumerate(prepared) if p.gbcms_status.startswith("PASS")]
        rs_variants = [prepared[i].variant for i in valid_indices]

        # Log validation results
        n_invalid = len(prepared) - len(valid_indices)
        logger.info(
            "Variant preparation: %d valid, %d rejected (%d total)",
            len(valid_indices),
            n_invalid,
            len(prepared),
        )
        invalid = [p for p in prepared if not p.gbcms_status.startswith("PASS")]
        for p in invalid[:5]:
            logger.warning(
                "Rejected variant: %s:%d %s>%s — %s",
                p.variant.chrom,
                p.original_pos + 1,
                p.original_ref,
                p.original_alt,
                p.gbcms_status,
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
            elif ref_len == alt_len and ref_len > 1:
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
            # Every variant was rejected during preparation (e.g. a contig mismatch →
            # FAIL_FETCH_FAILED, or FAIL_EMPTY_ALLELE). This is NOT an empty variant file
            # (that returns earlier, before any variant exists) — the variants are real,
            # so we still fall through and write them per sample with their FAIL_* reason
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
                strandedness=getattr(self.config, "strandedness", "reverse"),
                mfsd=self.config.output.mfsd,
                rna_editing_db=(
                    str(self.config.rna_editing_db)  # type: ignore[attr-defined]
                    if getattr(self.config, "rna_editing_db", None)
                    else None
                ),
                gtf_path=(
                    str(self.config.gtf)  # type: ignore[attr-defined]
                    if getattr(self.config, "gtf", None)
                    else None
                ),
                gtf_cache_dir=(
                    str(self.config.gtf_cache_dir)  # type: ignore[attr-defined]
                    if getattr(self.config, "gtf_cache_dir", None)
                    else None
                ),
                reference_fasta=str(self.config.reference_fasta),
                library_type=getattr(self.config, "library_type", "capture"),
            )
            rust_time = time.perf_counter() - rust_start
            logger.debug("Rust count_bam_binned completed in %.3fs", rust_time)

            # Update gbcms_status for variants where decomposed allele won
            for idx, counts in zip(valid_indices, counts_list, strict=True):
                if counts.used_decomposed:
                    prepared[idx].gbcms_status = "PASS;WARN_HOMOPOLYMER_DECOMP"

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

    def _compute_diagnostics(self, prepared: list, full_counts: list) -> None:
        """Compute post-counting diagnostic flags and populate gbcms_diagnostic.

        Diagnostic flags are semicolon-separated and stored in each
        PreparedVariant's gbcms_diagnostic field. Only PASS variants receive
        diagnostics; FAIL variants keep gbcms_diagnostic empty.

        Flags (per design §3):
            ZERO_ALT: ad == 0 and variant was successfully counted.
            PARTIAL_DOMINANT: partial_alt > ad (more structural evidence
                than confirmed ALT calls).
            MNP_DISC_RATIO(n/m): for MNPs (ref_len == alt_len > 1),
                always emitted showing discriminating position ratio.
            MNP_RESCUE_ELIGIBLE: disc/len ≤ rescue_mnp_threshold.
            HIGH_N_FRACTION(f): n_count / dp > 0.05 (duplex masking hotspot).
            NON_DISCRIMINATING_LOCUS: a sibling combination reconstructs the
                reference haplotype, so REF and ALT are sequence-indistinguishable
                and reads tie to NEITHER (explains a zeroed RD/AD at a covered locus).
        """
        flag_counts: dict[str, int] = {}

        for pv, counts in zip(prepared, full_counts, strict=True):
            # Only PASS variants get diagnostics; FAIL variants are not counted
            if not pv.gbcms_status.startswith("PASS"):
                continue

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
            ref_allele = pv.variant.ref_allele
            alt_allele = pv.variant.alt_allele
            ref_len = len(ref_allele)
            alt_len = len(alt_allele)
            if ref_len == alt_len and ref_len > 1:
                # Count positions where ref != alt (discriminating positions)
                disc = sum(1 for r, a in zip(ref_allele, alt_allele, strict=False) if r != a)
                ratio = disc / ref_len if ref_len > 0 else 0.0
                flags.append(f"MNP_DISC_RATIO({disc}/{ref_len})")
                if ratio <= self.config.rescue_mnp_threshold:
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

            # Populate the diagnostic field
            diagnostic = ";".join(flags)
            pv.gbcms_diagnostic = diagnostic

            # Track flag counts for logging
            for flag_name in flags:
                # Normalize parametric flags for counting
                base_flag = flag_name.split("(")[0]
                flag_counts[base_flag] = flag_counts.get(base_flag, 0) + 1

        # Log diagnostic summary
        if flag_counts:
            summary = ", ".join(f"{flag}={count}" for flag, count in sorted(flag_counts.items()))
            logger.info("Diagnostic flags: %s", summary)
        else:
            logger.debug("No diagnostic flags triggered")

    def _rescue_mnp_pass(
        self,
        prepared: list,
        full_counts: list,
        bam_path: Path,
        sample_name: str,
    ) -> None:
        """MNP rescue pass: decompose MNPs into individual SNPs and re-count.

        For each PASS variant where:
          - ad == 0 (no confirmed ALT reads)
          - variant is MNP (ref_len == alt_len > 1)
          - MNP_RESCUE_ELIGIBLE is flagged (disc/len ≤ rescue_mnp_threshold)

        The method creates synthetic SNP variants at each discriminating position,
        prepares and counts them via the Rust engine, and takes the best (highest
        alt_count) as the rescued value.

        After rescue, Invariant 1 (any_alt = ad + partial_alt) intentionally breaks.
        The original any_alt and partial_alt are forensic evidence from the original
        MNP check; only ad is updated with the rescued value.

        Design reference: validation_status_design.md §2, §5.

        Args:
            prepared: Full list of PreparedVariant objects.
            full_counts: Merged counts (one per variant).
            bam_path: Path to BAM file for re-counting.
            sample_name: Sample name for logging.
        """
        rescue_start = time.perf_counter()
        rs = _get_rs()

        # 1. Identify rescue candidates
        candidates: list[tuple[int, list[tuple[int, str, str]]]] = []
        for i, (pv, counts) in enumerate(zip(prepared, full_counts, strict=True)):
            if not pv.gbcms_status.startswith("PASS"):
                continue
            if counts.ad != 0:
                continue
            if "MNP_RESCUE_ELIGIBLE" not in pv.gbcms_diagnostic:
                continue

            ref_allele = pv.variant.ref_allele
            alt_allele = pv.variant.alt_allele
            ref_len = len(ref_allele)
            alt_len = len(alt_allele)
            if ref_len != alt_len or ref_len <= 1:
                continue

            # Extract discriminating positions (where ref != alt)
            disc_positions: list[tuple[int, str, str]] = []
            for offset in range(ref_len):
                if ref_allele[offset] != alt_allele[offset]:
                    abs_pos = pv.variant.pos + offset
                    disc_positions.append((abs_pos, ref_allele[offset], alt_allele[offset]))

            if disc_positions:
                candidates.append((i, disc_positions))

        if not candidates:
            logger.debug("MNP rescue: no candidates found for %s", sample_name)
            return

        logger.info("MNP rescue: %d candidate(s) for %s", len(candidates), sample_name)

        # 2. Build synthetic SNP variants for all candidates (batched)
        snp_variants: list = []
        snp_map: list[tuple[int, int]] = []  # (candidate_idx, disc_idx)

        for cand_idx, (pv_idx, disc_positions) in enumerate(candidates):
            pv = prepared[pv_idx]
            for disc_idx, (abs_pos, ref_base, alt_base) in enumerate(disc_positions):
                snp_v = rs.Variant(pv.variant.chrom, abs_pos, ref_base, alt_base, "SNP")
                snp_variants.append(snp_v)
                snp_map.append((cand_idx, disc_idx))

        # 3. Prepare synthetic SNPs (REF validation + ref_context)
        snp_prepared = rs.prepare_variants(
            snp_variants,
            str(self.config.reference_fasta),
            self.config.quality.context_padding,
            False,  # is_maf=False — these are synthetic
            self.config.threads,
            self.config.quality.adaptive_context,
        )

        # Filter to only valid SNPs
        valid_snp_indices = [
            j for j, sp in enumerate(snp_prepared) if sp.gbcms_status.startswith("PASS")
        ]
        valid_snp_variants = [snp_prepared[j].variant for j in valid_snp_indices]

        if not valid_snp_variants:
            logger.warning(
                "MNP rescue: all synthetic SNPs failed REF validation for %s",
                sample_name,
            )
            # Mark all candidates as failed rescue
            for _cand_idx, (pv_idx, _disc_positions) in enumerate(candidates):
                pv = prepared[pv_idx]
                pv.gbcms_rescue = "method=decomposed;original_alt=0;outcome=ref_validation_failed"
            return

        # 4. Count synthetic SNPs against the BAM
        # No decomposed variants or siblings for simple SNPs
        snp_decomposed = [None] * len(valid_snp_variants)
        snp_siblings: list[list] = [[] for _ in valid_snp_variants]

        align_cfg = self.config.alignment
        snp_counts = rs.count_bam_binned(
            str(bam_path),
            valid_snp_variants,
            snp_decomposed,
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
            sibling_variants=snp_siblings,
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
            strandedness=getattr(self.config, "strandedness", "reverse"),
            mfsd=self.config.output.mfsd,
            rna_editing_db=(
                str(self.config.rna_editing_db)  # type: ignore[attr-defined]
                if getattr(self.config, "rna_editing_db", None)
                else None
            ),
            gtf_path=(
                str(self.config.gtf)  # type: ignore[attr-defined]
                if getattr(self.config, "gtf", None)
                else None
            ),
            gtf_cache_dir=(
                str(self.config.gtf_cache_dir)  # type: ignore[attr-defined]
                if getattr(self.config, "gtf_cache_dir", None)
                else None
            ),
            reference_fasta=str(self.config.reference_fasta),
            library_type=getattr(self.config, "library_type", "capture"),
        )

        # 5. Map counts back to valid SNP indices
        # Build a full-index → count map for valid SNPs
        snp_count_by_idx: dict[int, Any] = {}
        for offset, j in enumerate(valid_snp_indices):
            snp_count_by_idx[j] = snp_counts[offset]

        # 6. For each candidate, find the best rescue position
        rescued_count = 0
        attempted_count = len(candidates)

        for cand_idx, (pv_idx, disc_positions) in enumerate(candidates):
            pv = prepared[pv_idx]
            counts = full_counts[pv_idx]

            best_alt = 0
            positions_str_parts: list[str] = []

            for disc_idx, (abs_pos, ref_base, alt_base) in enumerate(disc_positions):
                # Find the global SNP index for this disc position
                global_snp_idx = sum(len(candidates[c][1]) for c in range(cand_idx)) + disc_idx

                snp_alt = 0
                if global_snp_idx in snp_count_by_idx:
                    snp_alt = snp_count_by_idx[global_snp_idx].ad

                positions_str_parts.append(
                    f"{pv.variant.chrom}:{abs_pos + 1}({ref_base}>{alt_base}):{snp_alt}"
                )

                if snp_alt > best_alt:
                    best_alt = snp_alt

            positions_str = ",".join(positions_str_parts)

            if best_alt > 0:
                # Successful rescue: replace counts with a copy carrying the
                # best decomposed SNP alt_count.  BaseCounts is a frozen PyO3
                # struct (#[pyo3(get)] only), so we use copy-on-write via
                # with_ad() rather than direct field mutation.
                full_counts[pv_idx] = counts.with_ad(best_alt)
                pv.gbcms_rescue = f"method=decomposed;original_alt=0;positions={positions_str}"
                rescued_count += 1
                logger.debug(
                    "MNP rescue: %s:%d %s>%s → rescued alt=%d via decomposed SNPs",
                    pv.variant.chrom,
                    pv.variant.pos + 1,
                    pv.variant.ref_allele,
                    pv.variant.alt_allele,
                    best_alt,
                )
            else:
                # Failed rescue: no signal at any disc position
                pv.gbcms_rescue = (
                    f"method=decomposed;original_alt=0;outcome=no_signal;"
                    f"positions={positions_str}"
                )
                logger.debug(
                    "MNP rescue: %s:%d %s>%s → no signal at any disc position",
                    pv.variant.chrom,
                    pv.variant.pos + 1,
                    pv.variant.ref_allele,
                    pv.variant.alt_allele,
                )

        rescue_time = time.perf_counter() - rescue_start
        failed_count = attempted_count - rescued_count
        logger.info(
            "MNP rescue: %d/%d rescued, %d failed (%.3fs) for %s",
            rescued_count,
            attempted_count,
            failed_count,
            rescue_time,
            sample_name,
        )

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
                gbcms_diagnostic=pv.gbcms_diagnostic if pv else "",
                gbcms_rescue=pv.gbcms_rescue if pv else "",
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
