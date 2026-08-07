"""
CLI Entry Point: Exposes the gbcms functionality via command line.

Validation order (enforced for every option):
  1. Parse-time  — Typer Enum / min=/max= for constrained choices and numeric ranges.
  2. Pre-model   — Explicit checks in the command body before Pydantic construction
                   (file extensions, cross-option semantics, charset validation).
  3. Model-time  — Pydantic field constraints and validators in models/core.py.
  4. No silent skips — Missing inputs fail-fast unless the caller opts-out explicitly
                       (e.g. --lenient-bam).
"""

import logging
import os
import re
import sys
from pathlib import Path

import typer

from . import __version__
from .models.core import (
    AlignmentConfig,
    GbcmsDnaConfig,
    GbcmsRnaConfig,
    MergeConfig,
    OutputConfig,
    OutputFormat,
    QualityThresholds,
    ReadFilters,
    StrEnum,  # canonical backport (Python ≤ 3.10 compatible), defined in models.core
)
from .pipeline import Pipeline
from .utils import setup_logging

__all__ = ["app", "dna", "rna", "normalize", "merge"]

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constrained-choice enums (parse-time validation via Typer)
# ---------------------------------------------------------------------------


class AlignmentBackend(StrEnum):
    """CLI-exposed alignment backend options."""

    SW = "sw"
    HMM = "hmm"
    PAIRHMM = "pairhmm"


# Valid variant file extensions (checked before Pydantic config construction).
# .vcf.gz  — bgzip/gzip-compressed VCF (tabix-indexable, most common HPC format)
# .vcf.bgz — BGZip-compressed VCF, alternative extension used by some pipelines
# .vcf     — uncompressed VCF
# .maf     — Mutation Annotation Format
_VALID_VARIANT_EXTENSIONS: frozenset[str] = frozenset({".vcf", ".maf"})
_COMPRESSED_VCF_SUFFIXES: tuple[str, ...] = (".vcf.gz", ".vcf.bgz")

# Column-prefix charset: only letters, digits, underscores
_COLUMN_PREFIX_RE = re.compile(r"^[A-Za-z0-9_]*$")

app = typer.Typer(help="gbcms: Get Base Counts Multi-Sample")


def _log_command() -> str:
    """Log the full CLI command for reproducibility and return the command string.

    Logs the reconstructed command at INFO level so users and pipeline logs
    capture exactly what was invoked. Returns the string for embedding in
    VCF/MAF provenance headers.
    """
    command_line = " ".join(sys.argv)
    logger.info("Command: %s", command_line)
    return command_line


def _is_compressed_vcf(path: Path) -> bool:
    """Return True if *path* has a compressed-VCF suffix (.vcf.gz or .vcf.bgz).

    BGZip (.vcf.bgz) and gzip (.vcf.gz) are both block-gzip compatible and
    are handled identically by pysam.  We accept both so users can pass
    tabix-indexed .bgz files without renaming.
    """
    name_lower = path.name.lower()
    return any(name_lower.endswith(suffix) for suffix in _COMPRESSED_VCF_SUFFIXES)


def _exit_on_sample_failure(result: dict) -> None:
    """Propagate per-sample *failures* to the process exit code (HI-1).

    ``Pipeline.run()`` catches per-sample errors, records them in ``failed_samples``,
    and returns normally, so a run where a BAM failed (e.g. a Rust panic surfaced as
    ``PyErr``) would otherwise exit ``0`` and read as success to an orchestrator like
    Nextflow. This exits **non-zero** only when a sample actually failed.

    An **empty variant set is NOT a failure**: a sample can legitimately have no variants
    called, and per-sample workflows must not fail that task. Those runs process zero
    samples but record no ``failed_samples`` (the pipeline logs "No variants found …"),
    and they exit ``0`` here.

    Must be called **outside** the command's ``try/except Exception`` block: ``typer.Exit``
    subclasses ``RuntimeError``, so raising it inside that block would be swallowed and
    re-logged as a generic "Pipeline failed".
    """
    failed = result.get("failed_samples", [])
    if not failed:
        # Full success, OR a legitimately empty/rejected variant set — both exit 0.
        return

    processed = int(result.get("samples_processed", 0))
    if processed > 0:
        logger.error(
            "%d of %d sample(s) failed — exiting with code 1.",
            len(failed),
            processed + len(failed),
        )
    else:
        logger.error("All %d sample(s) failed — exiting with code 1.", len(failed))
    raise typer.Exit(code=1)


def version_callback(value: bool) -> None:
    """Print version and exit."""
    if value:
        typer.echo(f"gbcms {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool | None = typer.Option(
        None, "--version", callback=version_callback, is_eager=True, help="Show version and exit."
    ),
) -> None:
    """
    gbcms: Get Base Counts Multi-Sample
    """
    pass


@app.command()
def dna(
    # Input options
    variant_file: Path = typer.Option(
        ..., "--variants", "-v", help="Path to VCF or MAF file containing variants"
    ),
    bam_files: list[Path] | None = typer.Option(
        None, "--bam", "-b", help="Path to BAM/CRAM file(s). Can be specified multiple times."
    ),
    bam_list: Path | None = typer.Option(
        None, "--bam-list", "-L", help="File containing list of BAM paths (one per line)"
    ),
    reference: Path = typer.Option(..., "--fasta", "-f", help="Path to reference FASTA file"),
    # Output options
    output_dir: Path = typer.Option(
        ..., "--output-dir", "-o", help="Directory to write output files"
    ),
    output_format: OutputFormat = typer.Option(
        OutputFormat.VCF, "--format", help="Output format (vcf or maf)"
    ),
    output_suffix: str = typer.Option(
        "", "--suffix", "-S", help="Suffix to append to output filename (e.g. '.genotyped')"
    ),
    column_prefix: str = typer.Option(
        "",
        "--column-prefix",
        help=(
            "Prefix for gbcms count columns in MAF output. "
            "Default: no prefix (e.g., 'ref_count'). "
            "Use 't_' for legacy compatibility (e.g., 't_ref_count')."
        ),
    ),
    preserve_barcode: bool = typer.Option(
        False,
        "--preserve-barcode",
        help=(
            "Preserve original Tumor_Sample_Barcode from input MAF "
            "instead of overriding with BAM sample name. "
            "Only applies to MAF→MAF output."
        ),
    ),
    # mFSD options (DNA-only)
    mfsd: bool = typer.Option(
        False,
        "--mfsd",
        help=(
            "Enable Mutant Fragment Size Distribution (mFSD) analysis. "
            "Adds 41 mFSD columns (KS test, LLR, mean sizes, pairwise "
            "comparisons, derived metrics) to MAF output and 13 MFSD INFO "
            "fields to VCF. See docs/reference/counting-metrics.md#mfsd."
        ),
    ),
    mfsd_parquet: bool = typer.Option(
        False,
        "--mfsd-parquet",
        help=(
            "Write a companion <sample>.fsd.parquet with per-variant raw "
            "fragment size arrays (ref_sizes, alt_sizes). Enables downstream "
            "mFSD visualizations. Requires --mfsd."
        ),
    ),
    observations_parquet: bool = typer.Option(
        False,
        "--observations-parquet",
        help=(
            "Write a companion <sample>.observations.parquet with the per-molecule allele "
            "call at each variant (which molecule carried REF/ALT/N/OTHER). Enables "
            "read-backed phasing and allelic-imbalance analyses downstream. Counts are "
            "unchanged; rows are written from Rust, so memory is flat at panel scale."
        ),
    ),
    mfsd_report: bool = typer.Option(
        False,
        "--mfsd-report",
        help=(
            "Generate an interactive HTML report with per-variant fragment "
            "size distributions and CH-vs-ctDNA fragment origin signals. "
            "Implies --mfsd and --mfsd-parquet. Output: "
            "<sample>.mfsd_report.html alongside the main output."
        ),
    ),
    mfsd_report_min_alt: int = typer.Option(
        3,
        "--mfsd-report-min-alt",
        help="Minimum ALT fragment count for a variant to appear in the mFSD report.",
    ),
    mfsd_report_max_variants: int = typer.Option(
        20,
        "--mfsd-report-max-variants",
        help=(
            "Maximum number of variants to include in the mFSD report. "
            "Variants are ranked by ALT fragment count (descending). "
            "Use -1 for no limit."
        ),
    ),
    # BAQ (both modes)
    apply_baq: bool = typer.Option(
        False,
        "--apply-baq/--no-baq",
        help="Apply heuristic BAQ quality downgrade near indels. Off by default "
        "(BQSR/fgbio consensus already recalibrates).",
    ),
    # UMI (both modes)
    umi_tag: str | None = typer.Option(
        None,
        "--umi-tag",
        help="BAM tag for UMI barcode (e.g. 'RX'). Enables UMI-aware fragment grouping.",
    ),
    # Quality thresholds
    min_mapq: int = typer.Option(20, "--min-mapq", help="Minimum mapping quality"),
    min_baseq: int = typer.Option(20, "--min-baseq", help="Minimum base quality"),
    fragment_qual_threshold: int = typer.Option(
        10,
        "--fragment-qual-threshold",
        help=(
            "Quality difference threshold for fragment consensus. "
            "When R1 and R2 disagree, the higher-quality allele wins only if "
            "the difference exceeds this threshold; otherwise the fragment is discarded."
        ),
    ),
    context_padding: int = typer.Option(
        5,
        "--context-padding",
        min=1,
        max=50,
        help=(
            "Minimum flanking reference bases around indel/complex variants for "
            "haplotype construction and SW alignment. Range 1–50 enforced at "
            "parse time. Auto-increased in repeat regions when --adaptive-context is enabled."
        ),
    ),
    adaptive_context: bool = typer.Option(
        True,
        "--adaptive-context/--no-adaptive-context",
        help="Dynamically increase context padding in tandem repeat regions.",
    ),
    # Read filters
    filter_duplicates: bool = typer.Option(True, help="Filter duplicate reads"),
    filter_secondary: bool = typer.Option(
        True,
        help="Exclude secondary alignments from the read cache. They never count toward read-level depth regardless; turning this off admits them to FRAGMENT-level evidence (dpf/rdf/adf), which is what cross-locus phasing needs.",
    ),
    filter_supplementary: bool = typer.Option(
        True,
        help="Exclude supplementary alignments from the read cache. They never count toward read-level depth regardless; turning this off admits them to FRAGMENT-level evidence (dpf/rdf/adf), which is what cross-locus phasing needs.",
    ),
    filter_qc_failed: bool = typer.Option(True, help="Filter reads failing QC"),
    filter_improper_pair: bool = typer.Option(False, help="Filter improperly paired reads"),
    filter_indel: bool = typer.Option(False, help="Filter reads containing indels"),
    # Normalization
    show_normalization: bool = typer.Option(
        False,
        "--show-normalization",
        help="Add norm_* columns showing left-aligned coordinates in output.",
    ),
    # MNP rescue
    rescue_mnp: bool = typer.Option(
        False,
        "--rescue-mnp",
        help=(
            "Enable MNP rescue pass for multi-base substitutions. "
            "When alt_count=0, decomposes the MNP into individual SNPs "
            "and re-counts using the best discriminating position. "
            "Populates gbcms_rescue with a structured audit trail."
        ),
    ),
    rescue_mnp_threshold: float = typer.Option(
        1.0,
        "--rescue-mnp-threshold",
        min=0.0,
        max=1.0,
        help=(
            "Maximum disc/len ratio for MNP rescue eligibility (0.0–1.0). "
            "1.0 = rescue ALL MNPs (default, C++ compatible). "
            "0.5 = only rescue sparse MNPs (≤50%% discriminating positions). "
            "0.0 = disable rescue eligibility (diagnostics still emitted). "
            "Only relevant when --rescue-mnp is enabled."
        ),
    ),
    # Performance
    threads: int = typer.Option(
        1,
        "--threads",
        "-t",
        min=1,
        help="Total worker-thread budget for this sample. gbcms keeps all parallelism "
        "within this budget (Nextflow passes the task's allocated cores).",
    ),
    verbose: bool = typer.Option(False, "--verbose", "-V", help="Enable verbose debug logging"),
    trace: bool = typer.Option(
        False,
        "--trace",
        "-T",
        help="Enable per-read Rust trace logging (slow). Implies --verbose. "
        "Shows detailed per-read classification diagnostics from the counting engine.",
    ),
    # BAM robustness
    lenient_bam: bool = typer.Option(
        False,
        "--lenient-bam",
        help=(
            "Skip missing BAM files and continue with the remaining samples. "
            "Default (off): any missing BAM causes an immediate exit. "
            "Use this flag when running on batch lists where some files may be absent."
        ),
    ),
    # Alignment backend (advanced)
    alignment_backend: AlignmentBackend = typer.Option(
        AlignmentBackend.PAIRHMM,
        "--alignment-backend",
        help=(
            "Alignment backend for Phase 3 classification: "
            "'pairhmm' (WFA2 + PairHMM, default) or 'sw' (Smith-Waterman). "
            "Invalid values are rejected at parse time."
        ),
    ),
    hmm_llr_threshold: float = typer.Option(
        2.3,
        "--llr-threshold",
        help="PairHMM log-likelihood ratio threshold for confident calls (default: ln(10) ≈ 2.3).",
    ),
    hmm_gap_open: float = typer.Option(
        1e-4,
        "--gap-open-prob",
        help="PairHMM gap-open probability for non-repeat regions.",
    ),
    hmm_gap_extend: float = typer.Option(
        0.1,
        "--gap-extend-prob",
        help="PairHMM gap-extend probability for non-repeat regions.",
    ),
    hmm_gap_open_repeat: float = typer.Option(
        1e-2,
        "--repeat-gap-open-prob",
        help="PairHMM gap-open probability for tandem repeat regions.",
    ),
    hmm_gap_extend_repeat: float = typer.Option(
        0.5,
        "--repeat-gap-extend-prob",
        help="PairHMM gap-extend probability for tandem repeat regions.",
    ),
):
    """
    Count alleles in cfDNA/somatic DNA BAMs at known variant sites.
    """
    # ── 1. Logging (must be first so all subsequent checks log correctly) ──────
    setup_logging(verbose=verbose, trace=trace)
    command_line = _log_command()
    logger.info("Running gbcms v%s in DNA mode", __version__)
    if rescue_mnp:
        logger.info(
            "MNP rescue pass enabled (--rescue-mnp, threshold=%.2f)",
            rescue_mnp_threshold,
        )

    # ── 2. Pre-model validation (semantic + cross-option checks) ───────────────

    # GAP 12: Reject unsupported variant file extensions before any I/O.
    _is_vcf_gz = _is_compressed_vcf(variant_file)
    _ext = variant_file.suffix.lower()
    if not _is_vcf_gz and _ext not in _VALID_VARIANT_EXTENSIONS:
        logger.error(
            "Unsupported variant file extension '%s'. "
            "Expected .vcf, .vcf.gz, .vcf.bgz, or .maf. Got: %s",
            _ext,
            variant_file,
        )
        raise typer.Exit(code=1)

    # GAP 10: Validate --column-prefix charset (letters, digits, underscores only).
    if column_prefix and not _COLUMN_PREFIX_RE.match(column_prefix):
        logger.error(
            "Invalid --column-prefix '%s': only letters, digits, and underscores are allowed. "
            "Whitespace and special characters would produce malformed column names.",
            column_prefix,
        )
        raise typer.Exit(code=1)

    # GAP 9: Warn when --preserve-barcode is used with non-MAF input (it is a no-op).
    if preserve_barcode and not _is_vcf_gz and _ext != ".maf":
        logger.warning(
            "--preserve-barcode has no effect when the variant file is not a MAF "
            "(got '%s'). The BAM sample name will be used in all output rows.",
            variant_file.suffix,
        )

    # --mfsd-report implies --mfsd and --mfsd-parquet (user convenience; no need
    # to specify all three flags). Log the auto-enable so it's not a surprise.
    if mfsd_report:
        if not mfsd:
            mfsd = True
            logger.info("--mfsd-report implies --mfsd; auto-enabled.")
        if not mfsd_parquet:
            mfsd_parquet = True
            logger.info("--mfsd-report implies --mfsd-parquet; auto-enabled.")

    # Validate --mfsd-parquet requires --mfsd (CLI-level check matches model-level validator).
    if mfsd_parquet and not mfsd:
        logger.error(
            "--mfsd-parquet requires --mfsd to be enabled. "
            "Re-run with both --mfsd and --mfsd-parquet."
        )
        raise typer.Exit(code=1)

    # GAP 8: Advisory warning when threads exceeds available CPUs.
    cpu_count = os.cpu_count() or 1
    if threads > cpu_count:
        logger.warning(
            "--threads %d exceeds os.cpu_count() (%d). "
            "Performance may degrade due to CPU oversubscription.",
            threads,
            cpu_count,
        )

    # ── 3. Parse BAM inputs (fail-fast by default; --lenient-bam opts out) ─────
    bams_dict = _parse_bam_inputs(bam_files, bam_list, lenient=lenient_bam)

    if not bams_dict:
        logger.error(
            "No BAM files to process. Provide at least one via --bam <path> "
            "or --bam-list <file>. Use --lenient-bam to allow partial BAM lists."
        )
        raise typer.Exit(code=1)

    logger.info("Found %d BAM file(s) to process", len(bams_dict))
    logger.info(
        "Config: min_mapq=%d, apply_baq=%s, alignment_backend=%s, umi_tag=%s",
        min_mapq,
        apply_baq,
        alignment_backend.value,
        umi_tag or "none",
    )

    try:
        # Build nested config objects
        output_config = OutputConfig(
            directory=output_dir,
            format=output_format,
            suffix=output_suffix,
            column_prefix=column_prefix,
            preserve_barcode=preserve_barcode,
            mfsd=mfsd,
            mfsd_parquet=mfsd_parquet,
            observations_parquet=observations_parquet,
            mfsd_report=mfsd_report,
            mfsd_report_min_alt=mfsd_report_min_alt,
            mfsd_report_max_variants=mfsd_report_max_variants,
        )

        quality_config = QualityThresholds(
            min_mapping_quality=min_mapq,
            min_base_quality=min_baseq,
            fragment_qual_threshold=fragment_qual_threshold,
            context_padding=context_padding,
            adaptive_context=adaptive_context,
        )

        filter_config = ReadFilters(
            duplicates=filter_duplicates,
            secondary=filter_secondary,
            supplementary=filter_supplementary,
            qc_failed=filter_qc_failed,
            improper_pair=filter_improper_pair,
            indel=filter_indel,
        )

        # Pass .value so AlignmentConfig receives a plain str, not the enum wrapper.
        # This is required because AlignmentConfig.validate_backend operates on str.
        alignment_config = AlignmentConfig(
            backend=alignment_backend.value,
            hmm_llr_threshold=hmm_llr_threshold,
            hmm_gap_open=hmm_gap_open,
            hmm_gap_extend=hmm_gap_extend,
            hmm_gap_open_repeat=hmm_gap_open_repeat,
            hmm_gap_extend_repeat=hmm_gap_extend_repeat,
        )

        config = GbcmsDnaConfig(
            variant_file=variant_file,
            bam_files=bams_dict,
            reference_fasta=reference,
            output=output_config,
            quality=quality_config,
            filters=filter_config,
            threads=threads,
            command_line=command_line,
            alignment=alignment_config,
            show_normalization=show_normalization,
            apply_baq=apply_baq,
            umi_tag=umi_tag,
            rescue_mnp=rescue_mnp,
            rescue_mnp_threshold=rescue_mnp_threshold,
        )

        result = Pipeline(config).run()

    except Exception as e:
        logger.exception("Pipeline failed: %s", e)
        raise typer.Exit(code=1) from e

    # HI-1: exit non-zero if any sample failed (or none were processed). Outside the
    # try so typer.Exit isn't caught by `except Exception` above.
    _exit_on_sample_failure(result)


@app.command()
def rna(
    # Input options (shared with dna)
    variant_file: Path = typer.Option(
        ..., "--variants", "-v", help="Path to VCF or MAF file containing variants"
    ),
    bam_files: list[Path] | None = typer.Option(
        None, "--bam", "-b", help="Path to BAM/CRAM file(s). Can be specified multiple times."
    ),
    bam_list: Path | None = typer.Option(
        None, "--bam-list", "-L", help="File containing list of BAM paths (one per line)"
    ),
    reference: Path = typer.Option(..., "--fasta", "-f", help="Path to reference FASTA file"),
    # Output options (shared with dna, no mFSD)
    output_dir: Path = typer.Option(
        ..., "--output-dir", "-o", help="Directory to write output files"
    ),
    output_format: OutputFormat = typer.Option(
        OutputFormat.VCF, "--format", help="Output format (vcf or maf)"
    ),
    output_suffix: str = typer.Option(
        "", "--suffix", "-S", help="Suffix to append to output filename (e.g. '.genotyped')"
    ),
    column_prefix: str = typer.Option(
        "",
        "--column-prefix",
        help="Prefix for gbcms count columns in MAF output.",
    ),
    preserve_barcode: bool = typer.Option(
        False,
        "--preserve-barcode",
        help="Preserve original Tumor_Sample_Barcode from input MAF.",
    ),
    # BAQ (on by default for RNA — upstream pipelines typically lack BQSR)
    apply_baq: bool = typer.Option(
        True,
        "--apply-baq/--no-baq",
        help="Apply heuristic BAQ quality downgrade near indels and splice "
        "junctions. On by default for RNA (upstream pipelines typically "
        "lack BQSR). Disable with --no-baq if BAMs are pre-calibrated.",
    ),
    # UMI (shared)
    umi_tag: str | None = typer.Option(
        None,
        "--umi-tag",
        help="BAM tag for UMI barcode (e.g. 'RX'). Enables UMI-aware fragment grouping.",
    ),
    # RNA-specific options
    enforce_strandedness: bool = typer.Option(
        True,
        "--enforce-strandedness/--no-strandedness",
        help=(
            "Filter reads by strand orientation relative to gene strand. "
            "Disable for unstranded RNA-seq libraries (or pass --strandedness unstranded)."
        ),
    ),
    strandedness: str = typer.Option(
        "reverse",
        "--strandedness",
        help=(
            "RNA library strand protocol: 'reverse' (default; dUTP/fr-firststrand, "
            "featureCounts -s 2 — the FORTE default), 'forward' (fr-secondstrand, -s 1), "
            "or 'unstranded' (-s 0). Sets the read->transcript-strand fold used by both "
            "--enforce-strandedness and ASJD strand-discordance. 'unstranded' disables both."
        ),
    ),
    rna_editing_db: Path | None = typer.Option(
        None,
        "--rna-editing-db",
        help="Path to REDIportal TABLE1 file (tab-delimited) of known A-to-I RNA editing sites.",
    ),
    gtf: Path | None = typer.Option(
        None,
        "--gtf",
        help=(
            "Path to GTF annotation file (Ensembl/GENCODE). Enables exon "
            "boundary distance calculation and BAQ suppression at annotated "
            "splice junctions. Only chromosomes with variants are loaded."
        ),
    ),
    gtf_cache_dir: Path | None = typer.Option(
        None,
        "--gtf-cache-dir",
        help=(
            "Directory for caching the parsed GTF index. On first use the parsed "
            "annotation is written here; later runs over the same GTF and variant "
            "set reuse it, skipping the multi-second GTF text parse. Point every "
            "sample in a cohort at one shared directory to parse the GTF only once."
        ),
    ),
    # P5: Library type flag
    library_type: str = typer.Option(
        "capture",
        "--library-type",
        help=(
            "RNA library type: 'capture' (default, IDT xGen-style) or "
            "'amplicon'. In amplicon mode, R1/R2 pairs are treated as "
            "independent observations (no fragment consensus) and "
            "strandedness filtering is automatically disabled."
        ),
    ),
    # Quality thresholds (different defaults for RNA)
    min_mapq: int = typer.Option(
        1,
        "--min-mapq",
        help="Minimum mapping quality (default: 1 for STAR-aligned reads).",
    ),
    min_baseq: int = typer.Option(20, "--min-baseq", help="Minimum base quality"),
    fragment_qual_threshold: int = typer.Option(
        10,
        "--fragment-qual-threshold",
        help="Quality difference threshold for fragment consensus.",
    ),
    context_padding: int = typer.Option(
        5,
        "--context-padding",
        min=1,
        max=50,
        help="Minimum flanking reference bases around indel/complex variants.",
    ),
    adaptive_context: bool = typer.Option(
        True,
        "--adaptive-context/--no-adaptive-context",
        help="Dynamically increase context padding in tandem repeat regions.",
    ),
    # Read filters (shared)
    filter_duplicates: bool = typer.Option(True, help="Filter duplicate reads"),
    filter_secondary: bool = typer.Option(
        True,
        help="Exclude secondary alignments from the read cache. They never count toward read-level depth regardless; turning this off admits them to FRAGMENT-level evidence (dpf/rdf/adf), which is what cross-locus phasing needs.",
    ),
    filter_supplementary: bool = typer.Option(
        True,
        help="Exclude supplementary alignments from the read cache. They never count toward read-level depth regardless; turning this off admits them to FRAGMENT-level evidence (dpf/rdf/adf), which is what cross-locus phasing needs.",
    ),
    filter_qc_failed: bool = typer.Option(True, help="Filter reads failing QC"),
    filter_improper_pair: bool = typer.Option(False, help="Filter improperly paired reads"),
    filter_indel: bool = typer.Option(False, help="Filter reads containing indels"),
    # Normalization
    observations_parquet: bool = typer.Option(
        False,
        "--observations-parquet",
        help=(
            "Write a companion <sample>.observations.parquet with the per-molecule allele "
            "call at each variant (which molecule carried REF/ALT/N/OTHER). Enables "
            "read-backed phasing and allelic-imbalance analyses downstream. Counts are "
            "unchanged; rows are written from Rust, so memory is flat at panel scale."
        ),
    ),
    show_normalization: bool = typer.Option(
        False,
        "--show-normalization",
        help="Add norm_* columns showing left-aligned coordinates in output.",
    ),
    # MNP rescue
    rescue_mnp: bool = typer.Option(
        False,
        "--rescue-mnp",
        help=(
            "Enable MNP rescue pass for multi-base substitutions. "
            "When alt_count=0, decomposes the MNP into individual SNPs "
            "and re-counts using the best discriminating position. "
            "Populates gbcms_rescue with a structured audit trail."
        ),
    ),
    rescue_mnp_threshold: float = typer.Option(
        1.0,
        "--rescue-mnp-threshold",
        min=0.0,
        max=1.0,
        help=(
            "Maximum disc/len ratio for MNP rescue eligibility (0.0–1.0). "
            "1.0 = rescue ALL MNPs (default, C++ compatible). "
            "0.5 = only rescue sparse MNPs (≤50%% discriminating positions). "
            "0.0 = disable rescue eligibility (diagnostics still emitted). "
            "Only relevant when --rescue-mnp is enabled."
        ),
    ),
    # Performance
    threads: int = typer.Option(
        1,
        "--threads",
        "-t",
        min=1,
        help="Total worker-thread budget for this sample. gbcms keeps all parallelism "
        "within this budget (Nextflow passes the task's allocated cores).",
    ),
    verbose: bool = typer.Option(False, "--verbose", "-V", help="Enable verbose debug logging"),
    trace: bool = typer.Option(
        False,
        "--trace",
        "-T",
        help="Enable per-read Rust trace logging (slow). Implies --verbose.",
    ),
    lenient_bam: bool = typer.Option(
        False,
        "--lenient-bam",
        help="Skip missing BAM files and continue with the remaining samples.",
    ),
    # Alignment backend (RNA defaults to pairhmm = WFA2 + marginalized PairHMM)
    alignment_backend: AlignmentBackend = typer.Option(
        AlignmentBackend.PAIRHMM,
        "--alignment-backend",
        help=(
            "Alignment backend for Phase 3 classification: "
            "'sw' (Smith-Waterman) or 'pairhmm' (WFA2 + PairHMM, default for RNA)."
        ),
    ),
    hmm_llr_threshold: float = typer.Option(
        2.3,
        "--llr-threshold",
        help="PairHMM log-likelihood ratio threshold for confident calls.",
    ),
    hmm_gap_open: float = typer.Option(
        5e-3,
        "--gap-open-prob",
        help="PairHMM gap-open probability (default: 5e-3 for RT stutter tolerance).",
    ),
    hmm_gap_extend: float = typer.Option(
        0.25,
        "--gap-extend-prob",
        help="PairHMM gap-extend probability (default: 0.25 for RT stutter tolerance).",
    ),
    hmm_gap_open_repeat: float = typer.Option(
        1e-2,
        "--repeat-gap-open-prob",
        help="PairHMM gap-open probability for tandem repeat regions.",
    ),
    hmm_gap_extend_repeat: float = typer.Option(
        0.5,
        "--repeat-gap-extend-prob",
        help="PairHMM gap-extend probability for tandem repeat regions.",
    ),
):
    """
    Count alleles in RNA-seq BAMs with transcriptome-aware filtering.

    RNA mode includes:
    - NH:i:1 MAPQ rescue for STAR-aligned reads
    - dUTP strandedness filtering (--enforce-strandedness)
    - A-to-I RNA editing site flagging (--rna-editing-db)
    - GTF-informed splice boundary annotation (--gtf)
    - RT-aware PairHMM gap penalties
    """
    # ── 1. Logging ──
    setup_logging(verbose=verbose, trace=trace)
    command_line = _log_command()
    logger.info("Running gbcms v%s in RNA mode", __version__)
    if rescue_mnp:
        logger.info(
            "MNP rescue pass enabled (--rescue-mnp, threshold=%.2f)",
            rescue_mnp_threshold,
        )

    # ── 2. Pre-model validation ──
    _is_vcf_gz = _is_compressed_vcf(variant_file)
    _ext = variant_file.suffix.lower()
    if not _is_vcf_gz and _ext not in _VALID_VARIANT_EXTENSIONS:
        logger.error(
            "Unsupported variant file extension '%s'. "
            "Expected .vcf, .vcf.gz, .vcf.bgz, or .maf. Got: %s",
            _ext,
            variant_file,
        )
        raise typer.Exit(code=1)

    if column_prefix and not _COLUMN_PREFIX_RE.match(column_prefix):
        logger.error(
            "Invalid --column-prefix '%s': only letters, digits, and underscores are allowed.",
            column_prefix,
        )
        raise typer.Exit(code=1)

    if preserve_barcode and not _is_vcf_gz and _ext != ".maf":
        logger.warning(
            "--preserve-barcode has no effect when the variant file is not a MAF.",
        )

    cpu_count = os.cpu_count() or 1
    if threads > cpu_count:
        logger.warning(
            "--threads %d exceeds os.cpu_count() (%d).",
            threads,
            cpu_count,
        )

    # ── 3. Parse BAM inputs ──
    bams_dict = _parse_bam_inputs(bam_files, bam_list, lenient=lenient_bam)
    if not bams_dict:
        logger.error("No BAM files to process.")
        raise typer.Exit(code=1)

    logger.info("Found %d BAM file(s) to process", len(bams_dict))

    # Normalize the strand protocol (the model validates the value; we normalize here
    # for the interaction checks below).
    strandedness = strandedness.lower().strip()

    # Amplicon libraries are not stranded — treat them as unstranded so reads are not
    # folded under a protocol that does not apply.
    if library_type == "amplicon" and strandedness != "unstranded":
        logger.info(
            "--library-type=amplicon: forcing --strandedness=unstranded "
            "(amplicon libraries are not strand-specific)"
        )
        strandedness = "unstranded"

    # An unstranded protocol has no transcript strand, so strand enforcement is a no-op
    # — disable it explicitly rather than letting it silently pass every read.
    if strandedness == "unstranded" and enforce_strandedness:
        enforce_strandedness = False
        logger.info(
            "--strandedness=unstranded: auto-disabled --enforce-strandedness "
            "(no transcript strand to filter against)"
        )

    logger.info(
        "Config: min_mapq=%d, apply_baq=%s, alignment_backend=%s, "
        "enforce_strandedness=%s, strandedness=%s, library_type=%s, umi_tag=%s",
        min_mapq,
        apply_baq,
        alignment_backend.value,
        enforce_strandedness,
        strandedness,
        library_type,
        umi_tag or "none",
    )
    if rna_editing_db:
        logger.info("RNA editing database: %s", rna_editing_db)
    if gtf:
        logger.info("GTF annotation: %s", gtf)

    try:
        output_config = OutputConfig(
            directory=output_dir,
            format=output_format,
            suffix=output_suffix,
            column_prefix=column_prefix,
            preserve_barcode=preserve_barcode,
            mfsd=False,  # mFSD not applicable to RNA
            mfsd_parquet=False,
            observations_parquet=observations_parquet,
        )

        quality_config = QualityThresholds(
            min_mapping_quality=min_mapq,
            min_base_quality=min_baseq,
            fragment_qual_threshold=fragment_qual_threshold,
            context_padding=context_padding,
            adaptive_context=adaptive_context,
        )

        filter_config = ReadFilters(
            duplicates=filter_duplicates,
            secondary=filter_secondary,
            supplementary=filter_supplementary,
            qc_failed=filter_qc_failed,
            improper_pair=filter_improper_pair,
            indel=filter_indel,
        )

        alignment_config = AlignmentConfig(
            backend=alignment_backend.value,
            hmm_llr_threshold=hmm_llr_threshold,
            hmm_gap_open=hmm_gap_open,
            hmm_gap_extend=hmm_gap_extend,
            hmm_gap_open_repeat=hmm_gap_open_repeat,
            hmm_gap_extend_repeat=hmm_gap_extend_repeat,
        )

        config = GbcmsRnaConfig(
            variant_file=variant_file,
            bam_files=bams_dict,
            reference_fasta=reference,
            output=output_config,
            quality=quality_config,
            filters=filter_config,
            threads=threads,
            command_line=command_line,
            alignment=alignment_config,
            show_normalization=show_normalization,
            apply_baq=apply_baq,
            umi_tag=umi_tag,
            enforce_strandedness=enforce_strandedness,
            strandedness=strandedness,
            rna_editing_db=rna_editing_db,
            gtf=gtf,
            gtf_cache_dir=gtf_cache_dir,
            library_type=library_type,
            rescue_mnp=rescue_mnp,
            rescue_mnp_threshold=rescue_mnp_threshold,
        )

        result = Pipeline(config).run()

    except Exception as e:
        logger.exception("Pipeline failed: %s", e)
        raise typer.Exit(code=1) from e

    # HI-1: exit non-zero if any sample failed (or none were processed). Outside the
    # try so typer.Exit isn't caught by `except Exception` above.
    _exit_on_sample_failure(result)


@app.command("build-gtf-cache")
def build_gtf_cache(
    gtf: Path = typer.Option(
        ...,
        "--gtf",
        "-g",
        exists=True,
        help="Path to the GTF annotation file (Ensembl/GENCODE).",
    ),
    variants: Path = typer.Option(
        ...,
        "--variants",
        "-v",
        exists=True,
        help=(
            "Variant file (VCF/MAF) for the cohort. Only its chromosome set is used. "
            "It MUST be the same variant file the per-sample 'gbcms rna' runs use, so "
            "the cache key lines up and those runs reuse this entry."
        ),
    ),
    gtf_cache_dir: Path = typer.Option(
        ...,
        "--gtf-cache-dir",
        help=(
            "Shared directory to write the cache into (created if missing). Point every "
            "per-sample 'gbcms rna --gtf-cache-dir' at this same directory."
        ),
    ),
    verbose: bool = typer.Option(False, "--verbose", "-V", help="Enable verbose debug logging"),
):
    """
    Pre-build the GTF index cache so a cohort parses the GTF only once.

    Parses the GTF for the chromosomes covered by --variants and writes the
    serialized index into --gtf-cache-dir. Run this ONCE before fanning out the
    per-sample 'gbcms rna' jobs (all pointed at the same --gtf-cache-dir): each then
    loads the prebuilt index in ~0.05s instead of re-parsing the GTF (~9s).

    Why a separate step: when many samples launch concurrently they all cold-miss
    and each re-parses the GTF, so the cache alone saves nothing until a later wave.
    Building it up front lets every sample start warm.
    """
    from gbcms import _rs
    from gbcms.pipeline import read_variant_file

    setup_logging(verbose=verbose, trace=False)

    # Extension pre-check (mirrors the dna/rna/normalize commands).
    if (
        not _is_compressed_vcf(variants)
        and variants.suffix.lower() not in _VALID_VARIANT_EXTENSIONS
    ):
        logger.error(
            "Unsupported variant file extension '%s'. Expected .vcf, .vcf.gz, .vcf.bgz, or .maf.",
            variants.suffix,
        )
        raise typer.Exit(code=1)

    chroms = [v.chrom for v in read_variant_file(variants)]
    if not chroms:
        logger.error("No variants found in %s — nothing to scope the GTF cache to.", variants)
        raise typer.Exit(code=1)

    logger.info("Building GTF index cache for %d variants -> %s", len(chroms), gtf_cache_dir)
    n_exons = _rs.build_gtf_cache(str(gtf), chroms, str(gtf_cache_dir))
    logger.info(
        "GTF index cache ready in %s (%d exons across %d chromosomes). Per-sample runs "
        "using --gtf-cache-dir %s will now skip the GTF parse.",
        gtf_cache_dir,
        n_exons,
        len(set(chroms)),
        gtf_cache_dir,
    )


@app.command()
def normalize(
    variant_file: Path = typer.Option(
        ..., "--variants", "-v", help="Path to VCF or MAF file containing variants"
    ),
    reference: Path = typer.Option(..., "--fasta", "-f", help="Path to reference FASTA file"),
    output: Path = typer.Option(
        ..., "--output", "-o", help="Output file path (TSV with normalization results)"
    ),
    threads: int = typer.Option(
        1, "--threads", "-t", min=1, help="Total worker-thread budget for this run."
    ),
    verbose: bool = typer.Option(False, "--verbose", "-V", help="Enable verbose debug logging"),
    trace: bool = typer.Option(
        False,
        "--trace",
        "-T",
        help="Enable per-read Rust trace logging (slow). Implies --verbose.",
    ),
):
    """
    Normalize variants (left-align + validate REF) without counting.

    Reads variants from a VCF or MAF file, applies MAF anchor resolution,
    REF validation, and bcftools-style left-alignment, then writes results
    to a TSV file showing original and normalized coordinates.
    """
    from .normalize import normalize_variants

    setup_logging(verbose=verbose, trace=trace)

    # Apply the same file extension pre-check as the 'run' command.
    _is_vcf_gz = _is_compressed_vcf(variant_file)
    _ext = variant_file.suffix.lower()
    if not _is_vcf_gz and _ext not in _VALID_VARIANT_EXTENSIONS:
        logger.error(
            "Unsupported variant file extension '%s'. "
            "Expected .vcf, .vcf.gz, .vcf.bgz, or .maf. Got: %s",
            _ext,
            variant_file,
        )
        raise typer.Exit(code=1)

    normalize_variants(
        variant_file=variant_file,
        reference=reference,
        output=output,
        threads=threads,
    )


@app.command()
def merge(
    inputs: list[str] = typer.Option(
        ...,
        "--input",
        "-i",
        help=(
            "Input MAF in 'type:path' format (e.g., duplex:sample1-duplex.maf). "
            "Repeatable. At least 2 required."
        ),
    ),
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output merged MAF path.",
    ),
    add_combined: bool = typer.Option(
        True,
        "--add-combined/--no-combined",
        help=(
            "Compute additive simplex_duplex_* fragment columns "
            "when both duplex and simplex inputs are present."
        ),
    ),
    legacy_naming: bool = typer.Option(
        False,
        "--legacy-naming",
        help="Use t_{metric}_{type} naming (genotype_variants compatible).",
    ),
    verbose: bool = typer.Option(False, "--verbose", "-V", help="Enable debug logging"),
):
    """
    Merge per-BAM-type genotyped MAFs into a single type-prefixed output.

    Performs an outer join on the 5-column variant key (Chromosome,
    Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2),
    prefixes all gbcms count columns with the BAM type label, and
    optionally computes additive simplex+duplex combined columns.

    Example::

        gbcms merge \\
            --input duplex:sample1-duplex.maf \\
            --input simplex:sample1-simplex.maf \\
            --output sample1-merged.maf
    """
    from pydantic import ValidationError

    from .merge import merge_mafs

    setup_logging(verbose=verbose, trace=False)
    _log_command()
    logger.info("gbcms merge v%s", __version__)

    # ── Pre-model: parse type:path pairs ──────────────────────────────────
    parsed: dict[str, Path] = {}
    for inp in inputs:
        if ":" not in inp:
            logger.error(
                "Invalid --input format '%s'. " "Expected type:path (e.g., duplex:sample.maf)",
                inp,
            )
            raise typer.Exit(code=1)
        label, path_str = inp.split(":", 1)
        label = label.strip().lower()
        if not label:
            logger.error("Empty BAM type label in --input '%s'", inp)
            raise typer.Exit(code=1)
        if not _COLUMN_PREFIX_RE.match(label):
            logger.error(
                "Invalid BAM type label '%s' — only letters, digits, " "underscores allowed.",
                label,
            )
            raise typer.Exit(code=1)
        if label in parsed:
            logger.error(
                "Duplicate BAM type label '%s'. Each --input must have " "a unique type label.",
                label,
            )
            raise typer.Exit(code=1)
        parsed[label] = Path(path_str)

    # ── Model-time: Pydantic validation ───────────────────────────────────
    try:
        config = MergeConfig(
            inputs=parsed,
            output=output,
            add_combined=add_combined,
            legacy_naming=legacy_naming,
        )
    except ValidationError as e:
        logger.error("Configuration error: %s", e)
        raise typer.Exit(code=1) from e

    # ── Execute ───────────────────────────────────────────────────────────
    try:
        merge_mafs(config)
    except Exception as e:
        logger.exception("Merge failed: %s", e)
        raise typer.Exit(code=1) from e

    logger.info("Done.")


def _parse_bam_inputs(
    bam_files: list[Path] | None,
    bam_list: Path | None,
    *,
    lenient: bool = False,
) -> dict[str, Path]:
    """
    Parse BAM inputs from direct arguments and/or a BAM list file.

    Validation behaviour:
    - **Fail-fast (default)**: If any BAM path does not exist, all missing paths
      are logged at ERROR level and ``typer.Exit(code=1)`` is raised.
    - **Lenient mode** (``lenient=True``, enabled via ``--lenient-bam``): Missing
      paths are logged as errors but skipped; the run continues with the
      remaining samples.
    - **BAM list file not found**: Always fails immediately regardless of lenient
      mode.  The list file itself is a required input, not an optional sample.

    Args:
        bam_files: List of BAM paths (optionally with ``sample_id:path`` format).
        bam_list: Path to a file containing BAM paths (one per line,
            optionally ``sample_name<whitespace>path``).
        lenient: When True, skip missing BAM files instead of exiting.

    Returns:
        Dictionary mapping sample names to resolved BAM ``Path`` objects.

    Raises:
        typer.Exit: If any BAM file or the list file itself is missing and
            ``lenient`` is False.
    """
    bams_dict: dict[str, Path] = {}

    # ── 1. Process direct --bam arguments ────────────────────────────────────
    if bam_files:
        missing: list[str] = []
        for bam_arg in bam_files:
            sample_name, bam_path = _parse_bam_arg(bam_arg)

            if not bam_path.exists():
                logger.error("BAM file not found: %s", bam_path)
                missing.append(str(bam_path))
                continue

            logger.debug("Registered BAM sample '%s': %s", sample_name, bam_path)
            bams_dict[sample_name] = bam_path

        if missing and not lenient:
            logger.error(
                "%d BAM file(s) not found. "
                "Add --lenient-bam to skip missing files and continue with the rest.",
                len(missing),
            )
            raise typer.Exit(code=1)

    # ── 2. Process --bam-list file ────────────────────────────────────────────
    if bam_list:
        # The list file itself is always required — lenient mode does not apply here.
        if not bam_list.exists():
            logger.error(
                "BAM list file not found: %s. "
                "Note: --lenient-bam does not apply to the list file itself.",
                bam_list,
            )
            raise typer.Exit(code=1)

        logger.debug("Reading BAM list from: %s", bam_list)
        try:
            with open(bam_list) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue  # skip blanks and comment lines

                    parts = line.split()
                    if len(parts) >= 2:
                        sample_name = parts[0]
                        bam_path = Path(parts[1])
                    else:
                        bam_path = Path(parts[0])
                        sample_name = bam_path.stem

                    if not bam_path.exists():
                        # Upgraded from WARNING to ERROR — a missing BAM in the
                        # list is always unexpected, whether in lenient mode or not.
                        logger.error(
                            "BAM file from list not found: %s (sample '%s')",
                            bam_path,
                            sample_name,
                        )
                        continue

                    logger.debug("Registered BAM sample from list '%s': %s", sample_name, bam_path)
                    bams_dict[sample_name] = bam_path

        except OSError as e:
            logger.error("Error reading BAM list file %s: %s", bam_list, e)

    return bams_dict


def _parse_bam_arg(bam_arg: Path) -> tuple[str, Path]:
    """
    Parse a BAM argument that may be in sample_id:path format.

    Args:
        bam_arg: Path object (may contain sample_id:path as string).

    Returns:
        Tuple of (sample_name, bam_path).
    """
    bam_str = str(bam_arg)
    if ":" in bam_str:
        parts = bam_str.split(":", 1)
        return parts[0], Path(parts[1])
    return bam_arg.stem, bam_arg


if __name__ == "__main__":
    app()
