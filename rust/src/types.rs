use pyo3::prelude::*;

#[pyclass]
#[derive(Debug, Clone)]
pub struct Variant {
    #[pyo3(get, set)]
    pub chrom: String,
    #[pyo3(get, set)]
    pub pos: i64, // 0-based
    #[pyo3(get, set)]
    pub ref_allele: String,
    #[pyo3(get, set)]
    pub alt_allele: String,
    #[pyo3(get, set)]
    pub variant_type: String, // "SNP", "INSERTION", "DELETION", "COMPLEX"
    /// Reference sequence around the variant for windowed indel detection.
    /// Covers [ref_context_start, ref_context_start + len) in genomic coords.
    /// Used by Safeguard 3 to verify shifted indels are biologically valid.
    #[pyo3(get, set)]
    pub ref_context: Option<String>,
    /// Genomic start position (0-based) of the ref_context string.
    #[pyo3(get, set)]
    pub ref_context_start: i64,
    /// Span of the tandem repeat region surrounding the variant (0 if not in a repeat).
    /// Used to dynamically tune Smith-Waterman gap penalties: repeat_span >= 10
    /// triggers gap_extend = 0 to absorb polymerase slippage noise.
    #[pyo3(get, set)]
    pub repeat_span: usize,

    /// Transcript strand ('+' or '-') for RNA strandedness filtering.
    /// None in DNA mode — zero cost, no branching impact.
    #[pyo3(get, set)]
    pub gene_strand: Option<char>,
}

#[pymethods]
impl Variant {
    #[new]
    #[pyo3(signature = (chrom, pos, ref_allele, alt_allele, variant_type, ref_context=None, ref_context_start=0, repeat_span=0, gene_strand=None))]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        chrom: String,
        pos: i64,
        ref_allele: String,
        alt_allele: String,
        variant_type: String,
        ref_context: Option<String>,
        ref_context_start: i64,
        repeat_span: usize,
        gene_strand: Option<char>,
    ) -> Self {
        Variant {
            chrom,
            pos,
            ref_allele,
            alt_allele,
            variant_type,
            ref_context,
            ref_context_start,
            repeat_span,
            gene_strand,
        }
    }
}

#[pyclass]
#[derive(Debug, Clone, Default)]
pub struct BaseCounts {
    // Basic counts
    #[pyo3(get)]
    pub dp: u32,
    #[pyo3(get)]
    pub rd: u32,
    #[pyo3(get)]
    pub ad: u32,

    // Strand-specific counts
    #[pyo3(get)]
    pub dp_fwd: u32,
    #[pyo3(get)]
    pub rd_fwd: u32,
    #[pyo3(get)]
    pub ad_fwd: u32,
    #[pyo3(get)]
    pub dp_rev: u32,
    #[pyo3(get)]
    pub rd_rev: u32,
    #[pyo3(get)]
    pub ad_rev: u32,

    // Fragment counts (Majority Rule)
    #[pyo3(get)]
    pub dpf: u32,
    #[pyo3(get)]
    pub rdf: u32,
    #[pyo3(get)]
    pub adf: u32,

    // Fragment strand counts
    #[pyo3(get)]
    pub rdf_fwd: u32,
    #[pyo3(get)]
    pub rdf_rev: u32,
    #[pyo3(get)]
    pub adf_fwd: u32,
    #[pyo3(get)]
    pub adf_rev: u32,

    // Stats
    #[pyo3(get)]
    pub sb_pval: f64,
    #[pyo3(get)]
    pub sb_or: f64,
    #[pyo3(get)]
    pub fsb_pval: f64,
    #[pyo3(get)]
    pub fsb_or: f64,

    // Homopolymer decomposition
    /// True if the decomposed (corrected) allele was used because it
    /// produced a higher alt_count than the original allele.
    #[pyo3(get)]
    pub used_decomposed: bool,

    // ── mFSD: Fragment Size Distribution counts ───────────────────────────────
    // Number of fragments in each Krewlyzer size class (valid cfDNA range only;
    // 50–1000 bp with a non-zero TLEN). Gate downstream stats on mfsd_ks_valid.
    /// Fragments classified as REF and with a valid insert size.
    #[pyo3(get)]
    pub mfsd_ref_count: u32,
    /// Fragments classified as ALT and with a valid insert size.
    #[pyo3(get)]
    pub mfsd_alt_count: u32,
    /// Fragments that are neither REF nor ALT, nor an N base (third allele).
    #[pyo3(get)]
    pub mfsd_nonref_count: u32,
    /// Fragments where the base at the variant position was called 'N'.
    #[pyo3(get)]
    pub mfsd_n_count: u32,

    // ── mFSD: Mean fragment sizes ─────────────────────────────────────────────
    /// Mean fragment size (bp) for REF-classified fragments. 0.0 when empty.
    #[pyo3(get)]
    pub mfsd_ref_mean: f64,
    /// Mean fragment size (bp) for ALT-classified fragments. 0.0 when empty.
    #[pyo3(get)]
    pub mfsd_alt_mean: f64,
    /// Mean fragment size (bp) for NonREF-classified fragments. 0.0 when empty.
    #[pyo3(get)]
    pub mfsd_nonref_mean: f64,
    /// Mean fragment size (bp) for N-classified fragments. 0.0 when empty.
    #[pyo3(get)]
    pub mfsd_n_mean: f64,

    // ── mFSD: Log-Likelihood Ratios ───────────────────────────────────────────
    // LLR = Σ log(P_tumor(size) / P_healthy(size)) over all fragments in class.
    // Positive = tumor-like (short fragments); negative = healthy-like (long).
    /// LLR for ALT-classified fragments.
    #[pyo3(get)]
    pub mfsd_alt_llr: f64,
    /// LLR for REF-classified fragments.
    #[pyo3(get)]
    pub mfsd_ref_llr: f64,

    // ── mFSD: Pairwise KS comparisons (6 pairs × 3 values = 18 fields) ───────
    // Each triad: (delta = mean_A - mean_B, KS D-statistic, KS p-value).
    // NaN when either class has fewer than mfsd::MIN_FOR_KS (5) fragments.
    // Check mfsd_ks_valid (Python-derived) before interpreting these values.

    /// ALT vs REF: mean(ALT) − mean(REF)
    #[pyo3(get)]
    pub mfsd_delta_alt_ref: f64,
    /// ALT vs REF: KS D-statistic
    #[pyo3(get)]
    pub mfsd_ks_alt_ref: f64,
    /// ALT vs REF: KS p-value
    #[pyo3(get)]
    pub mfsd_pval_alt_ref: f64,
    /// ALT vs REF: Benjamini-Hochberg FDR q-value for `mfsd_pval_alt_ref`,
    /// corrected across all variants with a valid alt-vs-REF KS test in the
    /// sample. The report classifies TUMOR-LIKE/CH-LIKE on this q-value,
    /// not the raw p-value. NaN when the KS test was invalid (too few fragments);
    /// equals the p-value until the post-counting BH pass runs.
    #[pyo3(get)]
    pub mfsd_qval_alt_ref: f64,

    /// ALT vs NonREF: mean(ALT) − mean(NonREF)
    #[pyo3(get)]
    pub mfsd_delta_alt_nonref: f64,
    /// ALT vs NonREF: KS D-statistic
    #[pyo3(get)]
    pub mfsd_ks_alt_nonref: f64,
    /// ALT vs NonREF: KS p-value
    #[pyo3(get)]
    pub mfsd_pval_alt_nonref: f64,

    /// REF vs NonREF: mean(REF) − mean(NonREF)
    #[pyo3(get)]
    pub mfsd_delta_ref_nonref: f64,
    /// REF vs NonREF: KS D-statistic
    #[pyo3(get)]
    pub mfsd_ks_ref_nonref: f64,
    /// REF vs NonREF: KS p-value
    #[pyo3(get)]
    pub mfsd_pval_ref_nonref: f64,

    /// ALT vs N: mean(ALT) − mean(N)
    #[pyo3(get)]
    pub mfsd_delta_alt_n: f64,
    /// ALT vs N: KS D-statistic
    #[pyo3(get)]
    pub mfsd_ks_alt_n: f64,
    /// ALT vs N: KS p-value
    #[pyo3(get)]
    pub mfsd_pval_alt_n: f64,

    /// REF vs N: mean(REF) − mean(N)
    #[pyo3(get)]
    pub mfsd_delta_ref_n: f64,
    /// REF vs N: KS D-statistic
    #[pyo3(get)]
    pub mfsd_ks_ref_n: f64,
    /// REF vs N: KS p-value
    #[pyo3(get)]
    pub mfsd_pval_ref_n: f64,

    /// NonREF vs N: mean(NonREF) − mean(N)
    #[pyo3(get)]
    pub mfsd_delta_nonref_n: f64,
    /// NonREF vs N: KS D-statistic
    #[pyo3(get)]
    pub mfsd_ks_nonref_n: f64,
    /// NonREF vs N: KS p-value
    #[pyo3(get)]
    pub mfsd_pval_nonref_n: f64,

    // ── mFSD: Nucleosomal fraction fields ────────────────────────────────────
    // Sub-nucleosomal (<150bp) and mono-nucleosomal (150–200bp) fractions for
    // REF and ALT fragment populations. These enable CH-vs-ctDNA differentiation:
    // ctDNA tends to show sub-nucleosomal enrichment; CH mirrors REF distribution.
    // NaN when the denominator (class count) is zero.

    /// Fraction of REF fragments < 150bp (sub-nucleosomal).
    #[pyo3(get)]
    pub mfsd_sub_nuc_ref_frac: f64,
    /// Fraction of ALT fragments < 150bp (sub-nucleosomal).
    #[pyo3(get)]
    pub mfsd_sub_nuc_alt_frac: f64,
    /// Sub-nucleosomal enrichment ratio: ALT frac / REF frac.
    /// Values > 1.0 suggest ALT fragments are enriched in short sizes (ctDNA-like).
    /// Values ≈ 1.0 suggest ALT mirrors REF distribution (CH-like).
    #[pyo3(get)]
    pub mfsd_sub_nuc_enrichment: f64,
    /// Fraction of REF fragments in 150–200bp range (mono-nucleosomal).
    #[pyo3(get)]
    pub mfsd_mono_nuc_ref_frac: f64,
    /// Fraction of ALT fragments in 150–200bp range (mono-nucleosomal).
    #[pyo3(get)]
    pub mfsd_mono_nuc_alt_frac: f64,

    // ── mFSD: Raw size arrays (for --mfsd-parquet export) ────────────────────
    // Populated in all runs but only copied to disk when --mfsd-parquet is set.
    // NOT exported via PyO3 — written directly to Parquet by write_fsd_parquet()
    // in parquet_writer.rs, avoiding an FFI round-trip and the pyarrow dependency.
    /// Raw REF fragment sizes (bp). Internal only; use write_fsd_parquet() to persist.
    pub ref_sizes: Vec<u32>,
    /// Raw ALT fragment sizes (bp). Internal only; use write_fsd_parquet() to persist.
    pub alt_sizes: Vec<u32>,

    // ── Universal additions (both DNA and RNA modes) ───────────────────────
    /// Reads with MAPQ == 0 at this locus.
    #[pyo3(get)]
    pub mq0_count: u32,
    /// Median distance of ALT-supporting bases to read end.
    #[pyo3(get)]
    pub alt_dist_end_median: f64,
    /// Median distance of REF-supporting bases to read end.
    #[pyo3(get)]
    pub ref_dist_end_median: f64,
    /// ALT reads from singleton UMI families (no mate confirmation).
    #[pyo3(get)]
    pub singleton_alt_count: u32,
    /// ALT reads from duplex UMI families (both strands confirmed).
    #[pyo3(get)]
    pub duplex_alt_count: u32,

    // ── Decomposed ALT counting (diagnostic, all variant types) ──────────
    // Enables DMP-compatible "any evidence of ALT" counting alongside the
    // strict block-match AD. Invariant: any_alt = ad + partial_alt.
    //
    // For SNPs/indels: any_alt == ad, partial_alt == 0 (no partial concept).
    // For MNPs: any_alt >= ad when reads match ALT at some but not all
    //           discriminating positions.
    // For complex: any_alt >= ad from masked comparison partial matches.

    /// Reads with ANY evidence of ALT at ≥1 discriminating position.
    /// Relaxed counting: includes both full ALT matches (ad) and partial
    /// matches where some positions match ALT. DMP-compatible metric.
    #[pyo3(get)]
    pub any_alt: u32,

    /// Reads with PARTIAL ALT match only (some but not all discriminating
    /// positions match ALT). any_alt = ad + partial_alt.
    #[pyo3(get)]
    pub partial_alt: u32,

    /// Reads with N base at ≥1 discriminating position (NAD in VCF).
    /// N bases arise from duplex collapsing (fgbio masks disagreeing bases)
    /// or sequencer failure. These reads are uninformative (neither REF nor
    /// ALT) but tracking them separately from true third-allele reads enables
    /// downstream QC: n_count/DP ratio flags duplex masking hotspots.
    /// Follows bam-readcount's explicit N:count separation model.
    ///
    /// Depth decomposition: DP = RD + AD + partial_alt + n_count + other
    /// where "other" = third-allele + low-BQ reads (implicit).
    #[pyo3(get)]
    pub n_count: u32,

    // ── RNA-specific (zeroed in DNA mode via Default) ─────────────────────
    /// Total reads on the transcript sense strand (SEN in VCF).
    #[pyo3(get)]
    pub sense_depth: u32,
    /// Total reads on the antisense strand (ANT in VCF).
    #[pyo3(get)]
    pub antisense_depth: u32,
    /// ALT reads on the transcript sense strand (ASEN in VCF).
    #[pyo3(get)]
    pub sense_strand_alt_count: u32,
    /// ALT reads on the antisense strand.
    #[pyo3(get)]
    pub antisense_strand_alt_count: u32,
    /// True if locus overlaps a known A-to-I RNA editing site (RED in VCF).
    #[pyo3(get)]
    pub rna_editing_site_overlap: bool,
    /// Reads supporting ALT that span a splice junction — CIGAR N (SPL in VCF).
    #[pyo3(get)]
    pub splice_spanning_count: u32,

    // ── GTF-informed annotation (None when no GTF provided) ──────────────
    /// Distance (bp) to nearest annotated exon boundary (EBD in VCF).
    /// None when no GTF is provided. 0 = at boundary. Used for BAQ suppression
    /// at splice sites. Only populated in RNA mode with `--gtf`.
    #[pyo3(get)]
    pub exon_boundary_dist: Option<i32>,

    // ── Per-transcript counts (empty when no GTF or no overlap) ─────
    /// Per-transcript read-level counts (TXRC in VCF).
    /// Format: "ENST...:AD,RD,DP;ENST...:AD,RD,DP"
    /// Empty string when no GTF, no overlapping transcripts, or DNA mode.
    #[pyo3(get)]
    pub transcript_read_counts: String,

    /// Per-transcript fragment-level counts (TXFC in VCF).
    /// Format: "ENST...:ADF,RDF,DPF;ENST...:ADF,RDF,DPF"
    /// Empty string when no GTF, no overlapping transcripts, or DNA mode.
    #[pyo3(get)]
    pub transcript_fragment_counts: String,

    // ── Allele-Specific Junction Divergence (ASJD) ──────────────────
    // All fields default to false/0.0/"" via #[derive(Default)], producing
    // clean output when no GTF is provided or in DNA mode.

    /// True when Fisher's exact test p < 0.05 for junction divergence.
    #[pyo3(get)]
    pub asjd_flag: bool,
    /// Raw Fisher exact p-value (1.0 when no divergence or insufficient data).
    #[pyo3(get)]
    pub asjd_pval: f64,
    /// Benjamini-Hochberg FDR-corrected q-value (set post-counting).
    #[pyo3(get)]
    pub asjd_qval: f64,

    /// Dominant REF junction as "start-end" (empty when no junction reads).
    #[pyo3(get)]
    pub asjd_ref_junction: String,
    /// Dominant ALT junction as "start-end".
    #[pyo3(get)]
    pub asjd_alt_junction: String,

    /// Splice motif for REF junction (GT-AG, GC-AG, AT-AC, or OTHER).
    #[pyo3(get)]
    pub asjd_ref_motif: String,
    /// Splice motif for ALT junction.
    #[pyo3(get)]
    pub asjd_alt_motif: String,

    /// Whether the REF dominant junction matches a GTF-annotated intron.
    #[pyo3(get)]
    pub asjd_ref_known: bool,
    /// Whether the ALT dominant junction matches a GTF-annotated intron.
    #[pyo3(get)]
    pub asjd_alt_known: bool,

    /// REF reads supporting the dominant junction.
    #[pyo3(get)]
    pub asjd_n_ref_junc: u32,
    /// ALT reads supporting the dominant junction.
    #[pyo3(get)]
    pub asjd_n_alt_junc: u32,
    /// Total REF reads with any splice junction.
    #[pyo3(get)]
    pub asjd_n_ref_total: u32,
    /// Total ALT reads with any splice junction.
    #[pyo3(get)]
    pub asjd_n_alt_total: u32,

    /// Semicolon-separated diagnostic flags (e.g., "LOW_ALT_JUNC;NOVEL_ALT_JUNC").
    #[pyo3(get)]
    pub asjd_diagnostic: String,

    /// True when the pangenomic haplotype matrix has a REF-class haplotype byte-
    /// identical to an ALT-class one (a sibling combination reconstructs the
    /// reference). REF and ALT are then sequence-indistinguishable, so reads tie to
    /// NEITHER — surfaced as the `NON_DISCRIMINATING_LOCUS` diagnostic so the zeroed
    /// RD/AD is explained rather than silent. PairHMM backend only.
    #[pyo3(get)]
    pub non_discriminating_locus: bool,
}

#[pymethods]
impl BaseCounts {
    /// Return a copy with `ad` replaced by `new_ad`.
    ///
    /// Used by the Python MNP rescue pass (`--rescue-mnp`) which needs to
    /// update `ad` after decomposing a sparse MNP into individual SNPs.
    /// Preserves immutability of the original struct from the Python side —
    /// all fields are `#[pyo3(get)]` only, so this copy-on-write method is
    /// the only way to produce a modified `BaseCounts` from Python.
    fn with_ad(&self, new_ad: u32) -> BaseCounts {
        BaseCounts {
            ad: new_ad,
            ..self.clone()
        }
    }
}
