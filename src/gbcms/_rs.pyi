"""Type stubs for the gbcms._rs Rust extension module.

Kept in sync with rust/src/types.rs and rust/src/lib.rs.
"""

class Variant:
    chrom: str
    pos: int
    ref_allele: str
    alt_allele: str
    variant_type: str
    ref_context: str | None
    ref_context_start: int
    repeat_span: int
    gene_strand: str | None
    def __init__(
        self,
        chrom: str,
        pos: int,
        ref_allele: str,
        alt_allele: str,
        variant_type: str,
        ref_context: str | None = None,
        ref_context_start: int = 0,
        repeat_span: int = 0,
        gene_strand: str | None = None,
    ) -> None: ...

class BaseCounts:
    # Basic counts
    dp: int
    rd: int
    ad: int
    # Strand-specific counts
    dp_fwd: int
    rd_fwd: int
    ad_fwd: int
    dp_rev: int
    rd_rev: int
    ad_rev: int
    # Fragment counts
    dpf: int
    rdf: int
    adf: int
    rdf_fwd: int
    rdf_rev: int
    adf_fwd: int
    adf_rev: int
    # Stats
    sb_pval: float
    sb_or: float
    fsb_pval: float
    fsb_or: float
    used_decomposed: bool
    # mFSD counts
    mfsd_ref_count: int
    mfsd_alt_count: int
    mfsd_nonref_count: int
    mfsd_n_count: int
    # mFSD means
    mfsd_ref_mean: float
    mfsd_alt_mean: float
    mfsd_nonref_mean: float
    mfsd_n_mean: float
    # mFSD LLR
    mfsd_alt_llr: float
    mfsd_ref_llr: float
    # mFSD pairwise KS triads
    mfsd_delta_alt_ref: float
    mfsd_ks_alt_ref: float
    mfsd_pval_alt_ref: float
    mfsd_qval_alt_ref: float
    mfsd_delta_alt_nonref: float
    mfsd_ks_alt_nonref: float
    mfsd_pval_alt_nonref: float
    mfsd_delta_ref_nonref: float
    mfsd_ks_ref_nonref: float
    mfsd_pval_ref_nonref: float
    mfsd_delta_alt_n: float
    mfsd_ks_alt_n: float
    mfsd_pval_alt_n: float
    mfsd_delta_ref_n: float
    mfsd_ks_ref_n: float
    mfsd_pval_ref_n: float
    mfsd_delta_nonref_n: float
    mfsd_ks_nonref_n: float
    mfsd_pval_nonref_n: float
    # Sub-nucleosomal / mono-nucleosomal fractions
    mfsd_sub_nuc_ref_frac: float
    mfsd_sub_nuc_alt_frac: float
    mfsd_sub_nuc_enrichment: float
    mfsd_mono_nuc_ref_frac: float
    mfsd_mono_nuc_alt_frac: float
    # Universal additions (both modes)
    mq0_count: int
    alt_dist_end_median: float
    ref_dist_end_median: float
    singleton_alt_count: int
    duplex_alt_count: int
    # Decomposed ALT counting (diagnostic, all variant types)
    # Invariant: any_alt = ad + partial_alt
    any_alt: int
    partial_alt: int
    # N-base diagnostic: reads with N at ≥1 discriminating position (NAD in VCF).
    # Tracks duplex masking burden for QC. Follows bam-readcount N:count model.
    n_count: int
    # RNA-specific (zeroed in DNA mode)
    sense_depth: int
    antisense_depth: int
    sense_strand_alt_count: int
    antisense_strand_alt_count: int
    rna_editing_site_overlap: bool
    splice_spanning_count: int
    # GTF-informed annotation (None when no GTF)
    exon_boundary_dist: int | None
    # P4b: Per-transcript counts (empty string when no GTF or no overlap)
    transcript_read_counts: str
    transcript_fragment_counts: str
    # P4c: ASJD fields (defaults: flag=False, pval/qval=0.0, strings="", bools=False, ints=0)
    asjd_flag: bool
    asjd_pval: float
    asjd_qval: float
    asjd_ref_junction: str
    asjd_alt_junction: str
    asjd_ref_motif: str
    asjd_alt_motif: str
    asjd_ref_known: bool
    asjd_alt_known: bool
    asjd_n_ref_junc: int
    asjd_n_alt_junc: int
    asjd_n_ref_total: int
    asjd_n_alt_total: int
    asjd_diagnostic: str
    # NON_DISCRIMINATING_LOCUS marker: a sibling combo reconstructs REF, so REF/ALT
    # are sequence-indistinguishable (PairHMM backend) and reads tie to NEITHER.
    non_discriminating_locus: bool
    # Copy-on-write method for MNP rescue pass (BaseCounts is frozen from Python)
    def with_ad(self, new_ad: int) -> BaseCounts:
        """Return a copy with `ad` replaced by `new_ad`."""
        ...

class PreparedVariant:
    variant: Variant
    original_pos: int
    original_ref: str
    original_alt: str
    gbcms_status: str
    # Post-counting diagnostic flags (set by pipeline._compute_diagnostics).
    # Semicolon-separated. Empty string when no diagnostics.
    # Examples: "ZERO_ALT", "PARTIAL_DOMINANT;MNP_DISC_RATIO(2/5);MNP_RESCUE_ELIGIBLE".
    gbcms_diagnostic: str
    # Rescue audit trail (set by pipeline._rescue_mnp_pass).
    # Semicolon-separated key=value pairs. Empty string when no rescue attempted.
    # Only populated when --rescue-mnp is enabled.
    gbcms_rescue: str
    was_anchor_resolved: bool
    was_left_aligned: bool
    @property
    def was_normalized(self) -> bool: ...
    decomposed_variant: Variant | None
    multi_allelic_group: int | None

def count_bam(
    bam_path: str,
    variants: list[Variant],
    decomposed: list[Variant | None],
    min_mapq: int = 20,
    min_baseq: int = 20,
    filter_duplicates: bool = True,
    filter_secondary: bool = True,
    filter_supplementary: bool = True,
    filter_qc_failed: bool = True,
    filter_improper_pair: bool = False,
    filter_indel: bool = False,
    threads: int = 1,
    fragment_qual_threshold: int = 10,
    sibling_variants: list[list[Variant]] | None = None,
    alignment_backend: str = "sw",
    hmm_llr_threshold: float = 2.3,
    hmm_gap_open: float = 1e-4,
    hmm_gap_extend: float = 0.1,
    hmm_gap_open_repeat: float = 1e-2,
    hmm_gap_extend_repeat: float = 0.5,
    mode: str = "dna",
    enforce_strandedness: bool = False,
    strandedness: str = "reverse",
) -> list[BaseCounts]: ...
def count_bam_binned(
    bam_path: str,
    variants: list[Variant],
    decomposed: list[Variant | None],
    min_mapq: int = 20,
    min_baseq: int = 20,
    filter_duplicates: bool = True,
    filter_secondary: bool = True,
    filter_supplementary: bool = True,
    filter_qc_failed: bool = True,
    filter_improper_pair: bool = False,
    filter_indel: bool = False,
    threads: int = 1,
    fragment_qual_threshold: int = 10,
    sibling_variants: list[list[Variant]] | None = None,
    alignment_backend: str = "sw",
    hmm_llr_threshold: float = 2.3,
    hmm_gap_open: float = 1e-4,
    hmm_gap_extend: float = 0.1,
    hmm_gap_open_repeat: float = 1e-2,
    hmm_gap_extend_repeat: float = 0.5,
    apply_baq: bool = False,
    umi_tag: str | None = None,
    mode: str = "dna",
    enforce_strandedness: bool = False,
    strandedness: str = "reverse",
    mfsd: bool = False,
    rna_editing_db: str | None = None,
    gtf_path: str | None = None,
    gtf_cache_dir: str | None = None,
    reference_fasta: str | None = None,
    library_type: str = "capture",
) -> list[BaseCounts]: ...
def build_gtf_cache(
    gtf_path: str,
    variant_chroms: list[str],
    cache_dir: str,
) -> int: ...
def prepare_variants(
    variants: list[Variant],
    reference_fasta: str,
    context_padding: int,
    is_maf: bool,
    threads: int,
    adaptive_context: bool,
) -> list[PreparedVariant]: ...
def write_fsd_parquet(
    output_path: str,
    chroms: list[str],
    positions: list[int],
    refs: list[str],
    alts: list[str],
    counts: list[BaseCounts],
) -> None: ...
def fisher_exact_2x2(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    """Compute Fisher's exact test on a 2×2 contingency table.

    Args:
        a: Top-left cell (e.g., ref_forward)
        b: Top-right cell (e.g., ref_reverse)
        c: Bottom-left cell (e.g., alt_forward)
        d: Bottom-right cell (e.g., alt_reverse)

    Returns:
        (p_value, odds_ratio) tuple. p_value is the two-sided Fisher exact
        probability; odds_ratio is ad/bc (inf when bc=0).
    """
    ...
