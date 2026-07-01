"""
Core data models for gbcms.

This module defines the data models for variants, configuration, and nested
config groups (filters, quality thresholds, output settings).

Config hierarchy:
    GbcmsBaseConfig        — shared fields for all counting modes
    ├── GbcmsDnaConfig     — DNA/cfDNA-specific defaults (mFSD, sw backend)
    └── GbcmsRnaConfig     — RNA-specific defaults (strandedness, editing, pairhmm)
"""

import sys
from pathlib import Path

if sys.version_info >= (3, 11):
    from enum import StrEnum
else:
    from enum import Enum

    class StrEnum(str, Enum):
        """Backport of StrEnum for Python 3.10."""

        pass


from pydantic import BaseModel, Field, field_validator, model_validator

__all__ = [
    "VariantType",
    "Variant",
    "OutputFormat",
    "ReadFilters",
    "QualityThresholds",
    "AlignmentConfig",
    "OutputConfig",
    "GbcmsBaseConfig",
    "GbcmsDnaConfig",
    "GbcmsRnaConfig",
    "GbcmsConfig",  # deprecated alias for GbcmsDnaConfig
    "MergeConfig",
]


class VariantType(StrEnum):
    """Type of genomic variant."""

    SNP = "SNP"
    INSERTION = "INSERTION"
    DELETION = "DELETION"
    COMPLEX = "COMPLEX"


class Variant(BaseModel):
    """Normalized representation of a genomic variant."""

    chrom: str
    pos: int = Field(ge=0, description="0-based position of the variant")
    ref: str
    alt: str
    variant_type: VariantType

    # Original input metadata (optional)
    original_id: str | None = None
    metadata: dict[str, str] = Field(
        default_factory=dict, description="Original input metadata/columns"
    )


class OutputFormat(StrEnum):
    """Supported output formats for gbcms."""

    VCF = "vcf"
    MAF = "maf"


# =============================================================================
# Nested Configuration Models
# =============================================================================


class ReadFilters(BaseModel):
    """
    Filters for read selection during BAM processing.

    These flags control which reads are excluded from counting.
    When True, reads with the corresponding flag are filtered out.
    """

    duplicates: bool = Field(default=True, description="Filter duplicate reads")
    secondary: bool = Field(default=True, description="Filter secondary alignments")
    supplementary: bool = Field(default=True, description="Filter supplementary alignments")
    qc_failed: bool = Field(default=True, description="Filter reads failing QC")
    improper_pair: bool = Field(default=False, description="Filter improperly paired reads")
    indel: bool = Field(default=False, description="Filter reads containing indels")


class QualityThresholds(BaseModel):
    """Quality score thresholds for filtering reads and bases.

    Defaults are chosen to match the CLI defaults so that programmatic callers
    and CLI users observe identical behaviour without explicit overrides.
    """

    min_mapping_quality: int = Field(
        default=20, ge=0, description="Minimum mapping quality (MAPQ). Default 20 matches CLI."
    )
    min_base_quality: int = Field(
        default=20,
        ge=0,
        description=(
            "Minimum base quality (BQ). Default 20 matches CLI `--min-baseq` default. "
            "Previously this was 0 (maximally permissive), which diverged from the CLI."
        ),
    )
    fragment_qual_threshold: int = Field(
        default=10,
        ge=0,
        le=93,
        description=(
            "Quality difference threshold for fragment consensus. "
            "When R1 and R2 disagree, the allele with higher base quality wins "
            "only if the difference exceeds this threshold. Ambiguous fragments "
            "(within threshold) are discarded to preserve VAF accuracy."
        ),
    )
    context_padding: int = Field(
        default=5,
        ge=1,
        le=50,
        description=(
            "Number of flanking reference bases fetched around indel/complex "
            "variants for haplotype construction and Smith-Waterman alignment. "
            "Larger values improve sensitivity for shifted indels in repeat "
            "regions at minimal computational cost."
        ),
    )
    adaptive_context: bool = Field(
        default=True,
        description=(
            "Dynamically increase context padding in tandem repeat regions. "
            "When enabled, the effective padding is max(context_padding, "
            "repeat_span/2 + 3), capped at 50bp."
        ),
    )


class AlignmentConfig(BaseModel):
    """Alignment backend configuration for Phase 3 read classification.

    Controls which algorithm is used for haplotype-level alignment when
    variant-type-specific CIGAR matching is inconclusive. PairHMM (hmm)
    uses per-base quality-aware emission probabilities and is the default;
    Smith-Waterman (sw) is a simpler alternative.
    """

    backend: str = Field(
        default="pairhmm",
        description="Alignment backend: 'sw' (Smith-Waterman) or 'pairhmm' (PairHMM, default).",
    )
    hmm_llr_threshold: float = Field(
        default=2.3,
        gt=0.0,
        description=(
            "PairHMM log-likelihood ratio threshold for confident calls. "
            "ln(10) ≈ 2.3 = 10:1 odds ratio (default)."
        ),
    )
    hmm_gap_open: float = Field(
        default=1e-4,
        gt=0.0,
        lt=1.0,
        description="PairHMM gap-open probability for non-repeat regions.",
    )
    hmm_gap_extend: float = Field(
        default=0.1,
        gt=0.0,
        lt=1.0,
        description="PairHMM gap-extend probability for non-repeat regions.",
    )
    hmm_gap_open_repeat: float = Field(
        default=1e-2,
        gt=0.0,
        lt=1.0,
        description="PairHMM gap-open probability for tandem repeat regions.",
    )
    hmm_gap_extend_repeat: float = Field(
        default=0.5,
        gt=0.0,
        lt=1.0,
        description="PairHMM gap-extend probability for tandem repeat regions.",
    )

    @field_validator("backend")
    @classmethod
    def validate_backend(cls, v: str) -> str:
        """Validate backend is a supported value."""
        v = v.lower().strip()
        if v not in ("sw", "hmm", "pairhmm"):
            raise ValueError(f"Invalid alignment backend '{v}'. Must be 'sw' or 'pairhmm'.")
        # Normalize 'pairhmm' to 'hmm' for the Rust engine.
        if v == "pairhmm":
            v = "hmm"
        return v


class OutputConfig(BaseModel):
    """Output configuration settings."""

    directory: Path = Field(description="Directory to write output files")
    format: OutputFormat = Field(default=OutputFormat.VCF, description="Output format (vcf or maf)")
    suffix: str = Field(default="", description="Suffix to append to output filename")
    column_prefix: str = Field(
        default="",
        description=(
            "Prefix for gbcms count columns in MAF output. "
            "Default: no prefix (e.g., 'ref_count'). "
            "Use 't_' for legacy compatibility (e.g., 't_ref_count')."
        ),
    )
    preserve_barcode: bool = Field(
        default=False,
        description=(
            "When True, preserve the original Tumor_Sample_Barcode from "
            "input MAF instead of overriding with the BAM sample name. "
            "Only applies to MAF→MAF output; VCF→MAF always uses BAM name."
        ),
    )
    mfsd: bool = Field(
        default=False,
        description=(
            "Enable Mutant Fragment Size Distribution (mFSD) analysis. "
            "Adds 34 mFSD columns (KS test, LLR, mean sizes, pairwise comparisons, "
            "derived metrics) to MAF output and 7 MFSD INFO fields to VCF. "
            "Required when mfsd_parquet=True."
        ),
    )
    mfsd_parquet: bool = Field(
        default=False,
        description=(
            "Write a companion <sample>.fsd.parquet file with per-variant raw "
            "fragment size arrays (REF and ALT insert sizes). "
            "Enables downstream mFSD visualizations (e.g. density plots). "
            "Requires mfsd=True."
        ),
    )
    mfsd_report: bool = Field(
        default=False,
        description=(
            "Generate an interactive HTML report with per-variant fragment size "
            "distributions. Implies mfsd=True and mfsd_parquet=True. Output: "
            "<sample>.mfsd_report.html alongside the main output."
        ),
    )
    mfsd_report_min_alt: int = Field(
        default=3,
        ge=1,
        description=(
            "Minimum ALT fragment count for a variant to appear in the mFSD report. "
            "Variants with fewer ALT fragments are excluded from the report."
        ),
    )
    mfsd_report_max_variants: int = Field(
        default=20,
        description=(
            "Maximum number of variants to include in the mFSD report. "
            "Variants are ranked by ALT fragment count (descending). "
            "Use -1 for no limit."
        ),
    )

    @field_validator("directory")
    @classmethod
    def validate_output_dir(cls, v: Path) -> Path:
        """Ensure output path is not a file."""
        if v.exists() and v.is_file():
            raise ValueError(f"Output path must be a directory, not a file: {v}")
        return v

    @model_validator(mode="after")
    def validate_mfsd_dependencies(self) -> "OutputConfig":
        """Enforce mFSD dependency chain: mfsd_report → mfsd_parquet → mfsd.

        This is validated at model construction so both CLI and programmatic
        callers get the same fail-fast behaviour.
        """
        if self.mfsd_report and not self.mfsd_parquet:
            raise ValueError(
                "mfsd_report=True requires mfsd_parquet=True. "
                "Enable mFSD Parquet export before requesting report generation."
            )
        if self.mfsd_parquet and not self.mfsd:
            raise ValueError(
                "mfsd_parquet=True requires mfsd=True. "
                "Enable mFSD analysis before requesting Parquet export."
            )
        return self


class GbcmsBaseConfig(BaseModel):
    """
    Base configuration shared between DNA and RNA modes.

    Groups related settings into nested models for cleaner organization.
    Subclassed by GbcmsDnaConfig and GbcmsRnaConfig which override defaults.
    """

    # Input files
    variant_file: Path
    bam_files: dict[str, Path]  # sample_name -> bam_path
    reference_fasta: Path

    # Mode identifier — set by subclasses, used by Rust engine for branching
    mode: str = "dna"

    # Nested configuration groups
    output: OutputConfig
    filters: ReadFilters = Field(default_factory=ReadFilters)
    quality: QualityThresholds = Field(default_factory=QualityThresholds)

    # BAQ quality downgrade (both modes, off by default in base/DNA —
    # upstream BQSR/consensus already recalibrates base qualities;
    # applying BAQ on top may double-penalize. RNA overrides to True.)
    apply_baq: bool = Field(
        default=False,
        description=(
            "Apply heuristic BAQ (Base Alignment Quality) downgrade near "
            "indels and splice junctions. Subtracts 20 from base qualities "
            "within 5bp of alignment indels or CIGAR N (splice junctions) to "
            "reduce false positive variant calls from alignment artifacts. "
            "Off by default for DNA (upstream BQSR/consensus already "
            "recalibrates). Overridden to True for RNA mode."
        ),
    )

    # UMI-aware fragment grouping (both modes)
    umi_tag: str | None = Field(
        default=None,
        description=(
            "BAM tag containing the UMI barcode (e.g. 'RX'). When set, QNAME+UMI "
            "are hashed together for fragment grouping, treating reads with different "
            "UMIs as distinct molecules."
        ),
    )

    # Performance
    threads: int = Field(default=1, ge=1, description="Number of threads")

    # Provenance — full CLI command for VCF/MAF headers and log reproducibility
    command_line: str = Field(
        default="",
        description=(
            "Full reconstructed CLI command (sys.argv joined). Embedded in "
            "VCF ##gbcms_command and MAF #command provenance headers."
        ),
    )

    # Alignment backend
    alignment: AlignmentConfig = Field(default_factory=AlignmentConfig)

    # Advanced
    show_normalization: bool = Field(
        default=False,
        description="Add normalization columns showing left-aligned coordinates to output.",
    )

    # MNP rescue pass (v4.3.0 — design §2)
    rescue_mnp: bool = Field(
        default=False,
        description=(
            "Enable MNP rescue pass for multi-base substitutions. "
            "When ad=0, decomposes the MNP into individual SNP positions "
            "and re-counts using the best discriminating position. "
            "Controlled by --rescue-mnp-threshold for eligibility gating. "
            "Populates gbcms_rescue with a structured audit trail. "
            "Intentionally breaks Invariant 1 (any_alt = ad + partial_alt) "
            "for rescued variants."
        ),
    )
    rescue_mnp_threshold: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        description=(
            "Maximum discriminating-position ratio (disc/len) for MNP rescue "
            "eligibility. 1.0 = rescue ALL MNPs (C++ gbcms permissive behavior). "
            "0.5 = only rescue sparse MNPs (≤50%% discriminating positions). "
            "0.0 = disable rescue eligibility (diagnostics still emitted). "
            "Only relevant when --rescue-mnp is enabled."
        ),
    )

    @field_validator("variant_file", "reference_fasta")
    @classmethod
    def validate_file_exists(cls, v: Path) -> Path:
        """Validate that input files exist."""
        if not v.exists():
            raise ValueError(f"File not found: {v}")
        return v

    @model_validator(mode="after")
    def validate_bams(self) -> "GbcmsBaseConfig":
        """Validate that all BAM files exist."""
        for name, path in self.bam_files.items():
            if not path.exists():
                raise ValueError(f"BAM file for sample '{name}' not found: {path}")
        return self


class GbcmsDnaConfig(GbcmsBaseConfig):
    """DNA/cfDNA-specific configuration.

    Inherits all shared fields from GbcmsBaseConfig. Overrides:
    - mode = "dna"
    - Default MAPQ = 20
    - Default alignment backend = pairhmm
    - mFSD available
    """

    mode: str = "dna"


class GbcmsRnaConfig(GbcmsBaseConfig):
    """RNA-seq-specific configuration.

    Inherits all shared fields from GbcmsBaseConfig. Overrides:
    - mode = "rna"
    - Default MAPQ = 1 (STAR assigns low MAPQ to novel splice junctions)
    - Default alignment backend = pairhmm (WFA2 + marginalized PairHMM)
    - Gap penalties relaxed for RT stutter tolerance
    - Strandedness filtering enabled by default
    - RNA editing database support
    - BAQ enabled by default (RNA BAMs typically lack BQSR; the
      splice-proximity penalty reduces false positives at exon
      boundaries — see docs/reference/rna-splice-handling.md)
    """

    mode: str = "rna"

    # RNA BAMs typically do not go through GATK BQSR or consensus
    # calling (e.g., fgbio), so raw sequencer BQ values are preserved.
    # Enabling BAQ applies a -20 BQ penalty within 5bp of indels AND
    # splice junctions (CIGAR N/RefSkip), mimicking the effect of GATK
    # SplitNCigarReads' overhang clipping without modifying the BAM.
    apply_baq: bool = Field(
        default=True,
        description=(
            "Apply heuristic BAQ quality downgrade near indels and splice "
            "junctions. On by default for RNA because upstream pipelines "
            "typically do not run BQSR or consensus calling. "
            "Disable with --no-baq if BAMs have already been quality-adjusted."
        ),
    )

    # RNA-specific quality overrides
    quality: QualityThresholds = Field(
        default_factory=lambda: QualityThresholds(min_mapping_quality=1),
    )

    # RNA-specific alignment defaults (WFA2 + marginalized PairHMM)
    alignment: AlignmentConfig = Field(
        default_factory=lambda: AlignmentConfig(
            backend="pairhmm",
            hmm_gap_open=5e-3,
            hmm_gap_extend=0.25,
        ),
    )

    # dUTP sense-strand filtering for stranded RNA-seq libraries
    enforce_strandedness: bool = Field(
        default=True,
        description=(
            "Filter reads by strand orientation relative to gene_strand. "
            "In dUTP stranded libraries, R1 is antisense and R2 is sense. "
            "Antisense reads at variant sites are noise, not signal."
        ),
    )

    # REDIportal A-to-I RNA editing site database
    rna_editing_db: Path | None = Field(
        default=None,
        description=(
            "Path to REDIportal VCF of known RNA editing sites. "
            "Variants overlapping these sites are flagged (not filtered) "
            "to distinguish genuine mutations from A-to-I editing."
        ),
    )

    # GTF annotation for splice-site-aware counting
    gtf: Path | None = Field(
        default=None,
        description=(
            "Path to GTF annotation file (Ensembl/GENCODE). When provided, "
            "enables exon boundary distance calculation and BAQ suppression "
            "at annotated splice junctions. Adds exon_boundary_dist column "
            "to output. Only chromosomes with variants are loaded "
            "(variant-guided streaming)."
        ),
    )

    gtf_cache_dir: Path | None = Field(
        default=None,
        description=(
            "Directory for caching the parsed GTF index (M5a). When set, the parsed "
            "annotation intermediate is persisted here and reused across runs over "
            "the same GTF and variant set, skipping the GTF text parse. Intended as "
            "a directory shared across a Nextflow cohort so the GTF is parsed once. "
            "Caching is best-effort: any cache error falls back to a normal parse."
        ),
    )

    # P5: Library type flag — controls fragment consensus behavior
    library_type: str = Field(
        default="capture",
        description=(
            "RNA library type: 'capture' (default, IDT xGen-style) or "
            "'amplicon'. In amplicon mode, R1/R2 read pairs are treated "
            "as independent observations (no fragment consensus), and "
            "strandedness filtering is automatically disabled."
        ),
    )

    @field_validator("library_type")
    @classmethod
    def validate_library_type(cls, v: str) -> str:
        """Validate library type is supported."""
        v = v.lower().strip()
        if v not in ("capture", "amplicon"):
            raise ValueError(f"Invalid library_type '{v}'. Must be 'capture' or 'amplicon'.")
        return v

    # RNA library strand protocol — controls the read→transcript-strand fold
    strandedness: str = Field(
        default="reverse",
        description=(
            "RNA library strand protocol: 'reverse' (default; dUTP / fr-firststrand, "
            "featureCounts -s 2 — the FORTE pipeline default), 'forward' "
            "(fr-secondstrand, -s 1), or 'unstranded' (-s 0). Controls how a read's "
            "orientation folds to its transcript strand for both strandedness filtering "
            "and ASJD strand-discordance detection. 'unstranded' disables both."
        ),
    )

    @field_validator("strandedness")
    @classmethod
    def validate_strandedness(cls, v: str) -> str:
        """Validate the strand protocol is supported."""
        v = v.lower().strip()
        if v not in ("reverse", "forward", "unstranded"):
            raise ValueError(
                f"Invalid strandedness '{v}'. Must be 'reverse', 'forward', or 'unstranded'."
            )
        return v

    @field_validator("rna_editing_db")
    @classmethod
    def validate_editing_db(cls, v: Path | None) -> Path | None:
        """Validate RNA editing database file exists."""
        if v is not None and not v.exists():
            raise ValueError(f"RNA editing database not found: {v}")
        return v

    @field_validator("gtf")
    @classmethod
    def validate_gtf(cls, v: Path | None) -> Path | None:
        """Validate GTF annotation file exists and has correct extension."""
        if v is not None:
            if not v.exists():
                raise ValueError(f"GTF annotation file not found: {v}")
            if not v.name.endswith((".gtf", ".gtf.gz")):
                raise ValueError(f"GTF file must have .gtf or .gtf.gz extension, got: {v.name}")
        return v

    @model_validator(mode="after")
    def validate_amplicon_strandedness(self) -> "GbcmsRnaConfig":
        """P5: Auto-disable strandedness for amplicon libraries.

        Amplicon libraries are not strand-specific, so enforcing dUTP
        strandedness filtering would incorrectly discard ~50% of reads.
        This model-level guard ensures both CLI and programmatic callers
        get correct behavior without needing to know the dependency.
        """
        if self.library_type == "amplicon" and self.enforce_strandedness:
            import logging

            logging.getLogger("gbcms.models").warning(
                "library_type='amplicon' is incompatible with enforce_strandedness=True. "
                "Auto-disabling strandedness filtering (amplicon libraries are not strand-specific)."
            )
            self.enforce_strandedness = False
        return self


# Deprecated alias for backward compatibility — use GbcmsDnaConfig directly.
GbcmsConfig = GbcmsDnaConfig


class MergeConfig(BaseModel):
    """Configuration for ``gbcms merge`` subcommand.

    Merges multiple per-BAM-type genotyped MAFs into a single
    type-prefixed output MAF.  Primary use case: duplex + simplex
    ACCESS cfDNA genotyping.

    Attributes:
        inputs: Mapping of BAM type label to MAF path.
            Example: ``{"duplex": Path("s1-duplex.maf"), "simplex": Path("s1-simplex.maf")}``.
        output: Path for the merged output MAF.
        add_combined: When True and both ``duplex`` and ``simplex`` inputs
            are present, compute additive ``simplex_duplex_*`` columns.
        legacy_naming: When True, use ``t_{metric}_{type}`` column naming
            for backward compatibility with genotype_variants.
    """

    inputs: dict[str, Path] = Field(
        description=(
            "Mapping of BAM type label (e.g., 'duplex', 'simplex') to "
            "genotyped MAF file path. At least 2 required."
        ),
    )
    output: Path = Field(description="Output merged MAF file path.")
    add_combined: bool = Field(
        default=True,
        description=(
            "Compute additive simplex_duplex_* columns when both duplex "
            "and simplex inputs are present."
        ),
    )
    legacy_naming: bool = Field(
        default=False,
        description=(
            "Use t_{metric}_{type} column naming for backward compatibility "
            "with genotype_variants merge output."
        ),
    )

    @field_validator("inputs")
    @classmethod
    def validate_inputs(cls, v: dict[str, Path]) -> dict[str, Path]:
        """Require at least 2 inputs; validate all paths exist."""
        if len(v) < 2:
            raise ValueError(
                f"At least 2 input MAFs required for merge, got {len(v)}. "
                "Provide multiple --input type:path arguments."
            )
        for label, path in v.items():
            if not path.exists():
                raise ValueError(f"Input MAF for '{label}' not found: {path}")
        return v
