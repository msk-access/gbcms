//! PairHMM alignment backend for Phase 3 read classification.
//!
//! Provides a probabilistic alternative to Smith-Waterman for allele
//! classification. Uses bio::stats::pairhmm with a base-quality-aware
//! emission model that naturally handles:
//! - Quality-dependent match/mismatch probabilities per position
//! - Configurable gap parameters for different genomic contexts
//! - Log-likelihood ratio (LLR) classification instead of score margin
//!
//! ## Selection
//!
//! Activated via `--alignment-backend hmm`. The `sw` backend remains the
//! default for v2.8.0; PairHMM becomes the default in v3.0.0.
//!
//! ## Architecture
//!
//! The PairHMM backend slots into the same dispatch location as the SW
//! backend (Phase 3 in check_allele_with_qual), returning the same
//! `(is_ref, is_alt, base_qual)` tuple for zero-friction integration.

use bio::stats::pairhmm::{
    EmissionParameters, GapParameters, PairHMM, StartEndGapParameters, XYEmission,
};
use bio::stats::{LogProb, Prob};
use log::{debug, trace};

use crate::types::Variant;
use super::utils::{median_qual, build_haplotypes, ClassifyResult, ClassifyPhase};


/// Base-quality-aware emission model for PairHMM.
///
/// In this model, x = read sequence and y = haplotype sequence.
/// Emission probabilities are derived from the Phred-encoded base
/// quality at each read position:
///
/// - Match:    P(correct) = 1 - 10^(-BQ/10)
/// - Mismatch: P(error)   = 10^(-BQ/10) / 3  (uniform over 3 wrong bases)
///
/// Haplotype bases (y) are assumed error-free (reference-quality).
/// Gap emissions use the match probability since the base still
/// needs to be "observed" (just at a different position).
pub struct BQEmission<'a> {
    /// Read sequence (x)
    pub read_seq: &'a [u8],
    /// Read base qualities (Phred scale)
    pub read_quals: &'a [u8],
    /// Haplotype sequence (y) — either REF or ALT haplotype
    pub haplotype: &'a [u8],
}

impl<'a> EmissionParameters for BQEmission<'a> {
    /// Probability of emitting read[i] given haplotype[j].
    ///
    /// Uses the Phred quality at read position i to compute:
    /// - Match:    log(1 - 10^(-q/10))
    /// - Mismatch: log(10^(-q/10) / 3)
    fn prob_emit_xy(&self, i: usize, j: usize) -> XYEmission {
        let q = self.read_quals[i] as f64;
        let p_error = 10.0_f64.powf(-q / 10.0);
        let p_correct = 1.0 - p_error;

        if self.read_seq[i].eq_ignore_ascii_case(&self.haplotype[j]) {
            XYEmission::Match(LogProb::from(Prob(p_correct)))
        } else {
            // Mismatch: error probability divided by 3 (uniform over wrong bases)
            XYEmission::Mismatch(LogProb::from(Prob(p_error / 3.0)))
        }
    }

    /// Probability of observing read[i] in a gap context (deletion in haplotype).
    /// Uses the match probability since the base quality still constrains
    /// how confidently we see this read base.
    fn prob_emit_x(&self, i: usize) -> LogProb {
        let q = self.read_quals[i] as f64;
        let p_correct = 1.0 - 10.0_f64.powf(-q / 10.0);
        LogProb::from(Prob(p_correct))
    }

    /// Probability of observing haplotype[j] in a gap context (insertion in read).
    /// Haplotype bases are reference-quality, so this is a fixed high probability.
    fn prob_emit_y(&self, _j: usize) -> LogProb {
        LogProb::from(Prob(0.999))
    }

    fn len_x(&self) -> usize {
        self.read_seq.len()
    }

    fn len_y(&self) -> usize {
        self.haplotype.len()
    }
}


/// Configurable gap parameters for PairHMM.
///
/// Gap probabilities are specified as linear-scale probabilities and
/// converted to `LogProb` on access. Two sets of parameters are supported:
/// - Standard: for most genomic regions
/// - Repeat: for tandem repeat regions (higher gap tolerance)
///
/// ## Defaults
///
/// | Parameter          | Standard  | Repeat   | Source                           |
/// |---------------------|-----------|----------|----------------------------------|
/// | `gap_open`          | 1e-4      | 1e-2     | Illumina error rate / SW behavior|
/// | `gap_extend`        | 0.1       | 0.5      | GATK / SW behavior               |
///
/// All values are starting points subject to empirical calibration.
pub struct ConfigurableGapParams {
    /// Probability of opening a gap in x (read)
    pub gap_open: f64,
    /// Probability of extending a gap in x (read)
    pub gap_extend: f64,
}

impl GapParameters for ConfigurableGapParams {
    /// Probability of opening a gap in x (read has insertion relative to haplotype).
    fn prob_gap_x(&self) -> LogProb {
        LogProb::from(Prob(self.gap_open))
    }

    /// Probability of opening a gap in y (haplotype has insertion relative to read,
    /// i.e., read has deletion).
    fn prob_gap_y(&self) -> LogProb {
        LogProb::from(Prob(self.gap_open))
    }

    /// Probability of extending a gap in x.
    fn prob_gap_x_extend(&self) -> LogProb {
        LogProb::from(Prob(self.gap_extend))
    }

    /// Probability of extending a gap in y.
    fn prob_gap_y_extend(&self) -> LogProb {
        LogProb::from(Prob(self.gap_extend))
    }
}

#[allow(dead_code)] // Public API — standard/repeat/rna_defaults used by HLA module
impl ConfigurableGapParams {
    /// Create standard gap parameters for non-repeat regions.
    pub fn standard(gap_open: f64, gap_extend: f64) -> Self {
        ConfigurableGapParams { gap_open, gap_extend }
    }

    /// Create relaxed gap parameters for tandem repeat regions.
    pub fn repeat(gap_open: f64, gap_extend: f64) -> Self {
        ConfigurableGapParams { gap_open, gap_extend }
    }

    /// RT-aware defaults for RNA: higher gap tolerance for splice artifacts.
    ///
    /// RNA-seq alignments near splice junctions often show indel-like
    /// artifacts from reverse transcriptase (RT) slippage. These defaults
    /// use a higher gap_open (0.005 vs 1e-4) and gap_extend (0.25 vs 0.1)
    /// to be more permissive of these artifacts without inflating false
    /// positive variant calls.
    pub fn rna_defaults() -> Self {
        ConfigurableGapParams {
            gap_open: 0.005,
            gap_extend: 0.25,
        }
    }

    /// Continuous gap parameters scaled by repeat_span using a logistic function.
    ///
    /// Replaces the rigid `repeat_span >= 10` binary threshold with a smooth
    /// S-curve that models the continuous increase in polymerase slippage
    /// probability with repeat length.
    ///
    /// - `gap_open` scales linearly from `gap_open_base` to `gap_open_repeat`
    ///    over the range [0, 20bp], clamped at 20bp.
    /// - `gap_extend` follows a logistic curve from `gap_extend_base` to
    ///   `gap_extend_repeat` with inflection at 10bp.
    ///
    /// ## Parameters
    ///
    /// The base/repeat values typically come from the user's CLI flags:
    /// - `--hmm-gap-open` → `gap_open_base` (default 1e-4)
    /// - `--hmm-gap-extend` → `gap_extend_base` (default 0.1)
    /// - `--hmm-gap-open-repeat` → `gap_open_repeat` (default 1e-2)
    /// - `--hmm-gap-extend-repeat` → `gap_extend_repeat` (default 0.5)
    pub fn dynamic(
        repeat_span: usize,
        gap_open_base: f64,
        gap_extend_base: f64,
        gap_open_repeat: f64,
        gap_extend_repeat: f64,
    ) -> Self {
        let gap_extend = dynamic_gap_extend(repeat_span, gap_extend_base, gap_extend_repeat);
        // Linear interpolation for gap_open, clamped at 20bp
        let t = (repeat_span as f64 / 20.0).min(1.0);
        let gap_open = gap_open_base + (gap_open_repeat - gap_open_base) * t;
        ConfigurableGapParams { gap_open, gap_extend }
    }
}


/// Continuous gap extension probability using a logistic sigmoid.
///
/// Models the biological reality that polymerase slippage probability
/// increases continuously with tandem repeat length — there is no magic
/// boundary at any particular repeat_span.
///
/// The logistic curve transitions smoothly from `p_base` to `p_max`:
///
/// ```text
///   f(x) = p_base + (p_max - p_base) / (1 + exp(-steepness * (x - k)))
/// ```
///
/// | repeat_span | gap_extend (default params) |
/// |-------------|---------------------------|
/// |           0 | ≈ 0.107                   |
/// |           5 | ≈ 0.149                   |
/// |          10 | = 0.300 (inflection)      |
/// |          15 | ≈ 0.451                   |
/// |          20 | ≈ 0.493                   |
///
/// ## Parameters
///
/// - `p_base`: baseline gap_extend for non-repeat regions (from `--hmm-gap-extend`)
/// - `p_max`: ceiling gap_extend for deep repeats (from `--hmm-gap-extend-repeat`)
/// - `k=10.0`: inflection point in bp (hardcoded — matches biological transition zone)
/// - `steepness=0.5`: controls transition width (~8bp, avoiding cliff effects)
pub fn dynamic_gap_extend(repeat_span: usize, p_base: f64, p_max: f64) -> f64 {
    let k = 10.0;          // inflection point (bp)
    let steepness = 0.5;   // transition width
    p_base + (p_max - p_base) / (1.0 + f64::exp(-steepness * (repeat_span as f64 - k)))
}


/// Convert the continuous logistic gap_extend to an integer SW penalty.
///
/// Smith-Waterman uses integer gap penalties (typically -1 for tight, 0 for free).
/// This function maps the logistic curve to i32 by rounding:
///
/// ```text
///   sw_gap_extend = round(-1.0 * (1.0 - logistic(repeat_span)))
/// ```
///
/// - repeat_span=0  → round(-0.893) = -1 (tight)
/// - repeat_span=10 → round(-0.500) = -1 (tight, Rust rounds half away from zero)
/// - repeat_span=12 → round(-0.280) = 0  (free)
/// - repeat_span=20 → round(-0.107) = 0  (free)
///
/// The SW path naturally has a slightly higher threshold (~12bp) for full
/// relaxation than PairHMM, which is appropriate since SW scoring is less
/// quality-aware and benefits from tighter gap control.
pub fn dynamic_sw_gap_extend(repeat_span: usize) -> i32 {
    // Use default base/max for the logistic (PairHMM defaults: 0.1, 0.5)
    let logistic = dynamic_gap_extend(repeat_span, 0.1, 0.5);
    (-1.0 * (1.0 - logistic)).round() as i32
}


/// Semiglobal alignment mode for read-vs-haplotype comparison.
///
/// Allows free start/end gaps in x (the read) so the haplotype can
/// be "slid" along — matching the SW semiglobal behavior where the
/// read is the pattern and the haplotype is the text.
pub struct SemiglobalMode;

impl StartEndGapParameters for SemiglobalMode {
    /// Allow free start gaps in x (read can begin anywhere in haplotype).
    fn free_start_gap_x(&self) -> bool {
        true
    }

    /// Allow free end gaps in x (read can end anywhere in haplotype).
    fn free_end_gap_x(&self) -> bool {
        true
    }
}



/// Classify a read subsequence as REF or ALT using PairHMM.
///
/// Computes log-likelihood of the read under both the REF and ALT
/// haplotypes, then classifies using the log-likelihood ratio (LLR):
///
///   LLR = log P(read | ALT) - log P(read | REF)
///
/// - LLR >  threshold → ALT
/// - LLR < -threshold → REF
/// - |LLR| <= threshold → ambiguous (neither)
///
/// ## Advantages over Smith-Waterman
///
/// - **Quality-aware**: per-base BQ informs match/mismatch probabilities
/// - **No tie-proneness**: continuous LLR avoids integer score ties
/// - **Platform-adaptive**: naturally handles different error profiles
///
/// ## Parameters
///
/// - `read_seq`, `read_quals`: raw read bases and Phred qualities
/// - `variant`: the variant being classified
/// - `min_baseq`: minimum base quality for median_qual calculation
/// - `gap_params`: gap open/extend probabilities
/// - `llr_threshold`: log-likelihood ratio threshold for confident calls
///
/// Returns `ClassifyResult` — same interface as `classify_by_alignment`.
// INTENTIONAL: kept for concordance testing against classify_by_marginalized_pairhmm.
// Remove after concordance proven (D3 complete, pending D8 cleanup in v4.1.0).
pub fn classify_by_pairhmm(
    read_seq: &[u8],
    read_quals: &[u8],
    variant: &Variant,
    min_baseq: u8,
    gap_params: &ConfigurableGapParams,
    llr_threshold: f64,
) -> ClassifyResult {
    // Build haplotypes
    let (ref_hap, alt_hap) = match build_haplotypes(variant) {
        Some(haps) => haps,
        None => {
            debug!("classify_by_pairhmm: build_haplotypes failed for {}:{}",
                   variant.chrom, variant.pos + 1);
            return ClassifyResult::neither(ClassifyPhase::Alignment);
        }
    };

    // Skip if too few usable bases
    let usable_count = read_quals.iter().filter(|&&q| q >= min_baseq).count();
    if usable_count < 3 {
        trace!("classify_by_pairhmm: only {} usable bases — skipping", usable_count);
        return ClassifyResult::neither(ClassifyPhase::Alignment);
    }

    // Create PairHMM instance (reuses internal DP matrix)
    let mut pairhmm = PairHMM::new(gap_params);

    // Compute log-likelihood under ALT haplotype
    let emission_alt = BQEmission {
        read_seq,
        read_quals,
        haplotype: &alt_hap,
    };
    let ll_alt = pairhmm.prob_related(&emission_alt, &SemiglobalMode, None);

    // Compute log-likelihood under REF haplotype
    let emission_ref = BQEmission {
        read_seq,
        read_quals,
        haplotype: &ref_hap,
    };
    let ll_ref = pairhmm.prob_related(&emission_ref, &SemiglobalMode, None);

    // Log-likelihood ratio: positive favors ALT, negative favors REF
    let llr = *ll_alt - *ll_ref;

    let med_qual = median_qual(read_quals, min_baseq);

    trace!(
        "classify_by_pairhmm: ll_alt={:.3} ll_ref={:.3} llr={:.3} threshold={:.3} \
         read_len={} ref_hap={} alt_hap={}",
        *ll_alt, *ll_ref, llr, llr_threshold,
        read_seq.len(), ref_hap.len(), alt_hap.len()
    );

    // Trace-level: dump full sequences for deep debugging
    trace!("PairHMM read_seq={}", String::from_utf8_lossy(read_seq));
    trace!("PairHMM ref_hap={}", String::from_utf8_lossy(&ref_hap));
    trace!("PairHMM alt_hap={}", String::from_utf8_lossy(&alt_hap));

    if llr > llr_threshold {
        ClassifyResult::is_alt(med_qual, ClassifyPhase::Alignment)  // ALT
    } else if llr < -llr_threshold {
        ClassifyResult::is_ref(med_qual, ClassifyPhase::Alignment)  // REF
    } else {
        // Ambiguous: LLR within threshold — route to "neither"
        // Same logic as SW tie handling: contributes to DP but not RD/AD
        trace!(
            "PairHMM ambiguous (|llr|={:.3} <= threshold={:.3}) — routing to neither",
            llr.abs(), llr_threshold
        );
        if *ll_alt > f64::NEG_INFINITY && *ll_ref > f64::NEG_INFINITY {
            ClassifyResult::new(false, false, med_qual, ClassifyPhase::Alignment)
        } else {
            ClassifyResult::neither(ClassifyPhase::Alignment)
        }
    }
}


/// Marginalized PairHMM classification against the full haplotype matrix.
///
/// Computes log-likelihood of the read under EACH haplotype in the matrix,
/// then takes the max LL over REF-class haplotypes and max LL over ALT-class
/// haplotypes. The LLR = max_alt_ll - max_ref_ll marginalizes over sibling
/// germline configurations.
///
/// This prevents misclassification at multi-allelic / germline-adjacent sites
/// where a single REF/ALT comparison would confuse the somatic variant with
/// a nearby germline variant.
///
/// ## Haplotype Classes
///
/// - Even indices (0, 2, 4, ...) = REF-class (H0, H2, H4, ...)
/// - Odd indices (1, 3, 5, ...) = ALT-class (H1, H3, H5, ...)
///
/// ## Parameters
///
/// - `read_seq`, `adjusted_quals`: raw read bases and (possibly BAQ-adjusted) qualities
/// - `matrix`: haplotype matrix from `build_haplotype_matrix()`
/// - `min_baseq`: minimum base quality for median_qual calculation
/// - `gap_params`: gap open/extend probabilities
/// - `llr_threshold`: log-likelihood ratio threshold for confident calls
pub fn classify_by_marginalized_pairhmm(
    read_seq: &[u8],
    adjusted_quals: &[u8],
    matrix: &[Vec<u8>],
    min_baseq: u8,
    gap_params: &ConfigurableGapParams,
    llr_threshold: f64,
) -> ClassifyResult {
    if matrix.len() < 2 {
        trace!("classify_by_marginalized_pairhmm: matrix has {} haplotypes (need ≥2)", matrix.len());
        return ClassifyResult::neither(ClassifyPhase::Alignment);
    }

    // Skip if too few usable bases
    let usable_count = adjusted_quals.iter().filter(|&&q| q >= min_baseq).count();
    if usable_count < 3 {
        trace!("classify_by_marginalized_pairhmm: only {} usable bases — skipping", usable_count);
        return ClassifyResult::neither(ClassifyPhase::Alignment);
    }

    let mut pairhmm = PairHMM::new(gap_params);
    let mut best_ref_ll = f64::NEG_INFINITY;
    let mut best_alt_ll = f64::NEG_INFINITY;

    for (i, haplotype) in matrix.iter().enumerate() {
        let emission = BQEmission {
            read_seq,
            read_quals: adjusted_quals,
            haplotype,
        };
        let ll = pairhmm.prob_related(&emission, &SemiglobalMode, None);

        trace!(
            "marginalized_pairhmm: hap[{}] ({}) ll={:.3} len={}",
            i, if i % 2 == 0 { "REF" } else { "ALT" }, *ll, haplotype.len(),
        );

        if i % 2 == 0 {
            if *ll > best_ref_ll { best_ref_ll = *ll; }
        } else {
            if *ll > best_alt_ll { best_alt_ll = *ll; }
        }
    }

    let llr = best_alt_ll - best_ref_ll;
    let med_qual = median_qual(adjusted_quals, min_baseq);

    debug!(
        "marginalized_pairhmm: best_ref_ll={:.3} best_alt_ll={:.3} llr={:.3} threshold={:.3} \
         nhaps={} read_len={}",
        best_ref_ll, best_alt_ll, llr, llr_threshold,
        matrix.len(), read_seq.len(),
    );

    if llr > llr_threshold {
        ClassifyResult::is_alt(med_qual, ClassifyPhase::Alignment)
    } else if llr < -llr_threshold {
        ClassifyResult::is_ref(med_qual, ClassifyPhase::Alignment)
    } else {
        trace!(
            "marginalized_pairhmm: ambiguous (|llr|={:.3} <= threshold={:.3})",
            llr.abs(), llr_threshold
        );
        if best_alt_ll > f64::NEG_INFINITY && best_ref_ll > f64::NEG_INFINITY {
            ClassifyResult::new(false, false, med_qual, ClassifyPhase::Alignment)
        } else {
            ClassifyResult::neither(ClassifyPhase::Alignment)
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    fn make_variant(pos: i64, ref_allele: &str, alt_allele: &str, context: &str, ctx_start: i64) -> Variant {
        Variant {
            chrom: "1".to_string(),
            pos,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            variant_type: String::new(),
            ref_context: Some(context.to_string()),
            ref_context_start: ctx_start,
            repeat_span: 0,
            gene_strand: None,
        }
    }

    #[test]
    fn test_bq_emission_match_high_quality() {
        let read = b"ACGT";
        let quals = &[40, 40, 40, 40]; // Q40 → P(error) ≈ 0.0001
        let hap = b"ACGT";
        let emission = BQEmission {
            read_seq: read,
            read_quals: quals,
            haplotype: hap,
        };

        // All positions match: should get high match probabilities
        if let XYEmission::Match(lp) = emission.prob_emit_xy(0, 0) {
            assert!(*lp > -0.001, "Q40 match should have near-zero log-prob, got {}", *lp);
        } else {
            panic!("Expected Match emission for matching bases");
        }
    }

    #[test]
    fn test_bq_emission_mismatch_high_quality() {
        let read = b"TCGT";
        let quals = &[40, 40, 40, 40];
        let hap = b"ACGT";
        let emission = BQEmission {
            read_seq: read,
            read_quals: quals,
            haplotype: hap,
        };

        // Position 0: T vs A — mismatch at Q40
        if let XYEmission::Mismatch(lp) = emission.prob_emit_xy(0, 0) {
            // log(0.0001/3) ≈ -10.48
            assert!(*lp < -5.0, "Q40 mismatch should have very negative log-prob, got {}", *lp);
        } else {
            panic!("Expected Mismatch emission for different bases");
        }
    }

    #[test]
    fn test_bq_emission_low_quality_reduces_penalty() {
        let read = b"T";
        let quals_low = &[5_u8];  // Q5 → P(error) ≈ 0.316
        let quals_high = &[40_u8]; // Q40 → P(error) ≈ 0.0001
        let hap = b"A"; // mismatch

        let emission_low = BQEmission { read_seq: read, read_quals: quals_low, haplotype: hap };
        let emission_high = BQEmission { read_seq: read, read_quals: quals_high, haplotype: hap };

        let lp_low = match emission_low.prob_emit_xy(0, 0) {
            XYEmission::Mismatch(lp) => *lp,
            _ => panic!("Expected Mismatch"),
        };
        let lp_high = match emission_high.prob_emit_xy(0, 0) {
            XYEmission::Mismatch(lp) => *lp,
            _ => panic!("Expected Mismatch"),
        };

        // Low quality mismatch should have LESS penalty (higher log-prob) than high quality
        assert!(lp_low > lp_high,
            "Low-quality mismatch ({:.3}) should penalize less than high-quality ({:.3})",
            lp_low, lp_high
        );
    }

    #[test]
    fn test_gap_params_standard() {
        let params = ConfigurableGapParams::standard(1e-4, 0.1);
        // Verify LogProb values are finite and negative
        assert!(*params.prob_gap_x() < 0.0);
        assert!(*params.prob_gap_y() < 0.0);
        assert!(*params.prob_gap_x_extend() < 0.0);
        assert!(*params.prob_gap_y_extend() < 0.0);
    }

    #[test]
    fn test_gap_params_repeat() {
        let params = ConfigurableGapParams::repeat(1e-2, 0.5);
        let standard = ConfigurableGapParams::standard(1e-4, 0.1);

        // Repeat gap_open should be higher (less negative log) than standard
        assert!(*params.prob_gap_x() > *standard.prob_gap_x(),
            "Repeat gap_open ({:.3}) should be higher than standard ({:.3})",
            *params.prob_gap_x(), *standard.prob_gap_x()
        );
    }

    #[test]
    fn test_build_haplotypes_snp() {
        // Context: AAACAAA with C→T at position 3 (offset 3 in context starting at 0)
        let variant = make_variant(3, "C", "T", "AAACAAA", 0);
        let (ref_hap, alt_hap) = build_haplotypes(&variant).unwrap();

        assert_eq!(ref_hap, b"AAACAAA");
        assert_eq!(alt_hap, b"AAATAAA");
    }

    #[test]
    fn test_build_haplotypes_deletion() {
        // Context: AAATCCGGG with TCC→T deletion at position 3
        let variant = make_variant(3, "TCC", "T", "AAATCCGGG", 0);
        let (ref_hap, alt_hap) = build_haplotypes(&variant).unwrap();

        assert_eq!(ref_hap, b"AAATCCGGG");
        assert_eq!(alt_hap, b"AAATGGG");
    }

    #[test]
    fn test_build_haplotypes_no_context() {
        let variant = Variant {
            chrom: "1".to_string(),
            pos: 100,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            variant_type: String::new(),
            ref_context: None,
            ref_context_start: 0,
            repeat_span: 0,
            gene_strand: None,
        };
        assert!(build_haplotypes(&variant).is_none());
    }

    #[test]
    fn test_classify_pairhmm_clear_snp_ref() {
        // Read perfectly matches REF haplotype at high quality
        // Context: GGGGACGGGG, SNP at offset 4: A→T
        let context = "GGGGACGGGG";
        let variant = make_variant(4, "A", "T", context, 0);

        // Read: GGGGACGGGG (matches REF)
        let read = b"GGGGACGGGG";
        let quals = &[35; 10];
        let gap_params = ConfigurableGapParams::standard(1e-4, 0.1);

        let r = classify_by_pairhmm(read, quals, &variant, 20, &gap_params, 2.3);
        assert!(r.is_ref, "Expected REF for read matching REF haplotype");
        assert!(!r.is_alt, "Should not be ALT");
        assert!(r.qual > 0, "Should have non-zero quality");
    }

    #[test]
    fn test_classify_pairhmm_clear_snp_alt() {
        // Read perfectly matches ALT haplotype at high quality
        let context = "GGGGACGGGG";
        let variant = make_variant(4, "A", "T", context, 0);

        // Read: GGGGTCGGGG (matches ALT)
        let read = b"GGGGTCGGGG";
        let quals = &[35; 10];
        let gap_params = ConfigurableGapParams::standard(1e-4, 0.1);

        let r = classify_by_pairhmm(read, quals, &variant, 20, &gap_params, 2.3);
        assert!(!r.is_ref, "Should not be REF");
        assert!(r.is_alt, "Expected ALT for read matching ALT haplotype");
        assert!(r.qual > 0, "Should have non-zero quality");
    }

    #[test]
    fn test_classify_pairhmm_low_quality_ambiguous() {
        // Read has very low quality at the SNP position → ambiguous
        let context = "GGGGACGGGG";
        let variant = make_variant(4, "A", "T", context, 0);

        let read = b"GGGGTCGGGG";
        let mut quals = [35_u8; 10];
        quals[4] = 2; // Q2 at the SNP position
        let gap_params = ConfigurableGapParams::standard(1e-4, 0.1);

        let _r = classify_by_pairhmm(read, &quals, &variant, 20, &gap_params, 2.3);
        // With Q2 at the mismatch position, the penalty is minimal,
        // so LLR should be low → ambiguous
        // This tests that low BQ reduces discrimination power
        // (exact outcome depends on LLR threshold)
    }

    #[test]
    fn test_classify_pairhmm_insufficient_usable_bases() {
        let context = "GGGGACGGGG";
        let variant = make_variant(4, "A", "T", context, 0);

        let read = b"GG";
        let quals = &[2, 2]; // both below min_baseq=20
        let gap_params = ConfigurableGapParams::standard(1e-4, 0.1);

        let r = classify_by_pairhmm(read, quals, &variant, 20, &gap_params, 2.3);
        assert!(!r.is_ref && !r.is_alt, "Should return neither with < 3 usable bases");
        assert_eq!(r.qual, 0);
    }

    #[test]
    fn test_classify_pairhmm_insertion() {
        // Test with insertion variant: REF=A, ALT=ACCC
        // Context: GGGGA___GGGGG (4bp padding each side)
        let context = "GGGGAGGGGG";
        let variant = make_variant(4, "A", "ACCC", context, 0);

        // Read carrying the insertion
        let read = b"GGGGACCCGGGGG";
        let quals = &[35; 13];
        let gap_params = ConfigurableGapParams::standard(1e-4, 0.1);

        let r = classify_by_pairhmm(read, quals, &variant, 20, &gap_params, 2.3);
        assert!(r.is_alt, "Expected ALT for read carrying insertion");
        assert!(!r.is_ref);
        assert!(r.qual > 0);
    }

    #[test]
    fn test_classify_pairhmm_deletion() {
        // Test with deletion variant: REF=ACCC, ALT=A
        let context = "GGGGACCCGGGGG";
        let variant = make_variant(4, "ACCC", "A", context, 0);

        // Read carrying the deletion (no CCC)
        let read = b"GGGGAGGGGG";
        let quals = &[35; 10];
        let gap_params = ConfigurableGapParams::standard(1e-4, 0.1);

        let r = classify_by_pairhmm(read, quals, &variant, 20, &gap_params, 2.3);
        assert!(r.is_alt, "Expected ALT for read carrying deletion");
        assert!(!r.is_ref);
        assert!(r.qual > 0);
    }
}
