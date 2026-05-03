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

class PreparedVariant:
    variant: Variant
    original_pos: int
    original_ref: str
    original_alt: str
    validation_status: str
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
    rna_editing_db: str | None = None,
) -> list[BaseCounts]: ...
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
