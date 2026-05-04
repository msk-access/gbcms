//! Mutant Fragment Size Distribution (mFSD) statistics.
//!
//! Provides distributional statistics for comparing fragment size profiles
//! across REF, ALT, NonREF, and N fragment classes. Ported from Krewlyzer's
//! `mfsd.rs` implementation.
//!
//! ## Overview
//! Fragment size distributions carry biological signal in cfDNA:
//! - Healthy cfDNA peaks near 167 bp (mono-nucleosome protection)
//! - Tumor-derived cfDNA is enriched in shorter fragments (~120–145 bp)
//! - The KS test detects distributional shifts between allele classes
//! - The LLR scores each fragment class relative to a Gaussian tumor/healthy model
//!
//! All functions operate on raw, unweighted fragment size slices. GC correction
//! is not applied here — GC bias affects count depth, not fragment length, so
//! distributional tests are already unbiased on observed sizes.
//!
//! ## Usage
//! ```ignore
//! let physical = mfsd::calc_physical_insert_size(&record);
//! let (ks_d, ks_p) = mfsd::ks_test(&alt_sizes, &ref_sizes);
//! let llr = mfsd::calc_llr(&alt_sizes);
//! let mean = mfsd::calc_mean(&alt_sizes);
//! ```

use log::trace;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;

/// Minimum number of fragments required in each class for the KS test.
/// Below this threshold, `ks_test` returns `(f64::NAN, 1.0)` to signal
/// insufficient data rather than a spurious result.
pub const MIN_FOR_KS: usize = 5;

// ── Model Parameters ──────────────────────────────────────────────────────────

/// Gaussian model parameters for cfDNA fragment size classification.
///
/// Healthy cfDNA populates the mono-nucleosome window (~167 bp ± 30 bp).
/// Tumor-derived cfDNA is enriched in shorter, sub-nucleosomal fragments
/// (~145 bp ± 35 bp). These defaults match the MSK-ACCESS cohort as used
/// in the Krewlyzer mFSD implementation.
///
/// Per-study calibration may improve accuracy for other sequencing protocols.
pub struct LlrModelParams {
    /// Mean fragment size for healthy cfDNA (bp). Default: 167.0
    pub healthy_mu: f64,
    /// Std dev for healthy cfDNA distribution (bp). Default: 30.0
    pub healthy_sigma: f64,
    /// Mean fragment size for tumor-derived cfDNA (bp). Default: 145.0
    pub tumor_mu: f64,
    /// Std dev for tumor-derived cfDNA distribution (bp). Default: 35.0
    pub tumor_sigma: f64,
}

impl LlrModelParams {
    /// Human cfDNA defaults calibrated to the MSK-ACCESS cohort.
    pub fn human() -> Self {
        Self {
            healthy_mu: 167.0,
            healthy_sigma: 30.0,
            tumor_mu: 145.0,
            tumor_sigma: 35.0,
        }
    }
}

impl Default for LlrModelParams {
    fn default() -> Self {
        Self::human()
    }
}

// ── Core Statistical Functions ────────────────────────────────────────────────

/// Arithmetic mean of a fragment size slice.
///
/// Returns `0.0` for an empty slice (caller should check count before using).
///
/// # Arguments
/// * `v` – Slice of fragment sizes in base pairs.
pub fn calc_mean(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.iter().sum::<f64>() / v.len() as f64
}

/// Fraction of fragment sizes falling within `[lo, hi)`.
///
/// Returns `NaN` for an empty slice. Used for sub-nucleosomal (<150bp)
/// and mono-nucleosomal (150–200bp) fraction computation.
///
/// # Arguments
/// * `v`  – Slice of fragment sizes in base pairs.
/// * `lo` – Inclusive lower bound (bp).
/// * `hi` – Exclusive upper bound (bp). Use `f64::INFINITY` for open-ended.
pub fn calc_fraction_in_range(v: &[f64], lo: f64, hi: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    let count = v.iter().filter(|&&x| x >= lo && x < hi).count();
    count as f64 / v.len() as f64
}

/// Gaussian probability density function (unnormalised for LLR use).
///
/// # Arguments
/// * `x` – Fragment size (bp)
/// * `mu` – Distribution mean (bp)
/// * `sigma` – Distribution standard deviation (bp)
#[inline]
fn gaussian_pdf(x: f64, mu: f64, sigma: f64) -> f64 {
    let z = (x - mu) / sigma;
    (-0.5 * z * z).exp() / (sigma * std::f64::consts::TAU.sqrt())
}

/// Log-Likelihood Ratio for a fragment size slice vs. the human cfDNA model.
///
/// For each fragment, computes `log(P_tumor(size) / P_healthy(size))` and sums
/// the results. Positive totals indicate tumor-like fragment length enrichment;
/// negative totals indicate healthy-like (long) fragment enrichment.
///
/// Returns `0.0` for an empty slice.
///
/// # Arguments
/// * `lengths` – Fragment sizes in base pairs.
pub fn calc_llr(lengths: &[f64]) -> f64 {
    calc_llr_with_params(lengths, &LlrModelParams::human())
}

/// Log-Likelihood Ratio with caller-supplied model parameters.
///
/// Internal version used for testing alternative models.
///
/// # Arguments
/// * `lengths` – Fragment sizes in base pairs.
/// * `params` – Model parameters (healthy/tumor mu and sigma).
pub fn calc_llr_with_params(lengths: &[f64], params: &LlrModelParams) -> f64 {
    if lengths.is_empty() {
        return 0.0;
    }
    lengths.iter().map(|&x| {
        let p_tumor   = gaussian_pdf(x, params.tumor_mu,   params.tumor_sigma);
        let p_healthy = gaussian_pdf(x, params.healthy_mu, params.healthy_sigma);
        // Guard: clamp denominator to avoid log(0) for extreme sizes
        let ratio = if p_healthy < f64::EPSILON {
            f64::INFINITY
        } else {
            p_tumor / p_healthy
        };
        ratio.ln()
    }).sum()
}

// ── Physical Fragment Sizing ─────────────────────────────────────────────────

/// Compute the physical fragment insert size from CIGAR, correcting TLEN for indels.
///
/// BAM TLEN measures the reference span between the outermost aligned bases of a
/// read pair. For fragments carrying indels, TLEN ≠ physical fragment length:
/// - A deletion makes TLEN *longer* than the physical DNA fragment
/// - An insertion makes TLEN *shorter* than the physical DNA fragment
///
/// **Formula:** `physical_size = |TLEN| - sum(D_ops) + sum(I_ops)`
///
/// Validated on 6 real MSK-ACCESS duplex BAMs (EGFR 15bp del, MET 35bp del,
/// ERBB2 12bp ins, KIT 6bp ins, KRAS G12D SNP, TP53 DNP). All corrections
/// match expected indel sizes exactly; REF fragments and SNPs are unaffected.
///
/// For MSK-ACCESS cfDNA (~167bp fragments, ~150bp reads), R1 and R2 overlap
/// by ~133bp, so both reads carry the same indel CIGAR. The caller (`observe()`)
/// stores `min(R1, R2)` as a defensive measure for non-cfDNA contexts where
/// only one read may span the indel.
///
/// # Arguments
/// * `record` – BAM record (must have valid CIGAR and TLEN)
///
/// # Returns
/// Physical insert size in bp (always ≥ 0). Returns `0` if TLEN is 0 (unpaired).
pub fn calc_physical_insert_size(record: &Record) -> i32 {
    let raw_tlen = record.insert_size();
    if raw_tlen == 0 {
        return 0; // Unpaired/unmapped mate — no correction possible
    }

    let abs_tlen = raw_tlen.unsigned_abs() as i32;
    let mut del_bp: i32 = 0;
    let mut ins_bp: i32 = 0;

    for op in record.cigar().iter() {
        match op {
            Cigar::Del(n) => del_bp += *n as i32,
            Cigar::Ins(n) => ins_bp += *n as i32,
            _ => {}
        }
    }

    let physical = abs_tlen - del_bp + ins_bp;

    // Guard: pathological CIGAR where correction overshoots (should never happen
    // with well-formed BAMs, but avoid returning nonsensical negative sizes)
    if physical <= 0 {
        trace!(
            "calc_physical_insert_size: pathological correction — TLEN={} D={} I={} → {} (clamped to |TLEN|)",
            raw_tlen, del_bp, ins_bp, physical
        );
        return abs_tlen; // Fall back to raw |TLEN| rather than a nonsensical value
    }

    trace!(
        "calc_physical_insert_size: TLEN={} D={} I={} → physical={}",
        raw_tlen, del_bp, ins_bp, physical
    );

    physical
}

// ── Two-Sample Kolmogorov-Smirnov Test ───────────────────────────────────────

/// Two-sample Kolmogorov-Smirnov test.
///
/// Computes the KS D-statistic (maximum absolute difference between empirical
/// CDFs) and an approximate p-value using the Kolmogorov distribution series.
///
/// Returns `(f64::NAN, 1.0)` if either slice has fewer than [`MIN_FOR_KS`]
/// fragments — callers should check `mfsd_ks_valid` before interpreting results.
///
/// # Arguments
/// * `a` – Fragment sizes for the first class (e.g., ALT fragments).
/// * `b` – Fragment sizes for the second class (e.g., REF fragments).
///
/// # Returns
/// `(d_statistic, p_value)`
pub fn ks_test(a: &[f64], b: &[f64]) -> (f64, f64) {
    if a.len() < MIN_FOR_KS || b.len() < MIN_FOR_KS {
        return (f64::NAN, 1.0);
    }

    // Sort copies for CDF walks
    let mut a_sorted = a.to_vec();
    let mut b_sorted = b.to_vec();
    a_sorted.sort_unstable_by(|x, y| x.partial_cmp(y).unwrap());
    b_sorted.sort_unstable_by(|x, y| x.partial_cmp(y).unwrap());

    let n = a_sorted.len() as f64;
    let m = b_sorted.len() as f64;

    // Walk merged sorted values computing CDF difference at each step
    let mut d: f64 = 0.0;
    let mut i = 0usize;
    let mut j = 0usize;

    // Merge-walk to track CDF of each sample
    while i < a_sorted.len() && j < b_sorted.len() {
        let val = if a_sorted[i] <= b_sorted[j] {
            a_sorted[i]
        } else {
            b_sorted[j]
        };
        // Advance all entries equal to `val` in both arrays
        while i < a_sorted.len() && a_sorted[i] <= val { i += 1; }
        while j < b_sorted.len() && b_sorted[j] <= val { j += 1; }

        let cdf_a = i as f64 / n;
        let cdf_b = j as f64 / m;
        d = d.max((cdf_a - cdf_b).abs());
    }

    let p = ks_p_value(d, n, m);
    (d, p)
}

/// Approximate p-value for a KS D-statistic using the Kolmogorov distribution.
///
/// Uses the series approximation: Q_KS(λ) = 2 Σ (-1)^(k-1) exp(-2k²λ²)
/// where λ = D * sqrt(n*m / (n+m)).
///
/// The series converges rapidly; we truncate at 100 terms (negligible error).
fn ks_p_value(d: f64, n: f64, m: f64) -> f64 {
    let lambda = d * (n * m / (n + m)).sqrt();
    if lambda < f64::EPSILON {
        return 1.0;
    }
    let mut sum = 0.0;
    for k in 1..=100i64 {
        let term = (-2.0 * (k as f64 * lambda).powi(2)).exp();
        let signed = if k % 2 == 0 { -1.0 } else { 1.0 } * term;
        sum += signed;
        // Early convergence: term negligible
        if term < 1e-12 { break; }
    }
    (2.0 * sum).clamp(0.0, 1.0)
}

// ── Unit Tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ─── calc_mean ────────────────────────────────────────────────────────────

    #[test]
    fn test_calc_mean_empty() {
        assert_eq!(calc_mean(&[]), 0.0);
    }

    #[test]
    fn test_calc_mean_single() {
        assert_eq!(calc_mean(&[100.0]), 100.0);
    }

    #[test]
    fn test_calc_mean_values() {
        let v: Vec<f64> = (1..=5).map(|x| x as f64).collect();
        assert!((calc_mean(&v) - 3.0).abs() < 1e-10);
    }

    // ─── calc_llr ─────────────────────────────────────────────────────────────

    #[test]
    fn test_llr_empty() {
        assert_eq!(calc_llr(&[]), 0.0);
    }

    #[test]
    fn test_llr_tumor_like() {
        // Short fragments (145 bp) are more likely under the tumor model
        let sizes: Vec<f64> = vec![145.0; 20];
        assert!(calc_llr(&sizes) > 0.0, "LLR should be positive for tumor-like sizes");
    }

    #[test]
    fn test_llr_healthy_like() {
        // Long fragments (167 bp) are more likely under the healthy model
        let sizes: Vec<f64> = vec![167.0; 20];
        assert!(calc_llr(&sizes) < 0.0, "LLR should be negative for healthy-like sizes");
    }

    // ─── ks_test ──────────────────────────────────────────────────────────────

    #[test]
    fn test_ks_insufficient_n() {
        let a: Vec<f64> = vec![100.0; 3]; // below MIN_FOR_KS
        let b: Vec<f64> = vec![200.0; 10];
        let (d, p) = ks_test(&a, &b);
        assert!(d.is_nan(), "D should be NaN for insufficient n");
        assert_eq!(p, 1.0, "p-value should be 1.0 for insufficient n");
    }

    #[test]
    fn test_ks_identical_distributions() {
        let v: Vec<f64> = (100..=200).map(|x| x as f64).collect();
        let (d, p) = ks_test(&v, &v);
        assert!(d.abs() < 1e-10, "D should be ~0 for identical distributions");
        // p is not tested here as it depends on the approximation for D=0
        let _ = p;
    }

    #[test]
    fn test_ks_distinct_distributions() {
        // Two non-overlapping distributions — D should be 1.0
        let a: Vec<f64> = (100..=150).map(|x| x as f64).collect();
        let b: Vec<f64> = (200..=250).map(|x| x as f64).collect();
        let (d, p) = ks_test(&a, &b);
        assert!((d - 1.0).abs() < 1e-10, "D should be ~1.0 for non-overlapping distributions, got {d}");
        assert!(p < 0.05, "p should be < 0.05 for well-separated distributions, got {p}");
    }

    // ─── calc_physical_insert_size ────────────────────────────────────────────

    use rust_htslib::bam::record::CigarString;

    /// Helper: build a minimal BAM record with given CIGAR and TLEN.
    fn mock_record(cigar: &CigarString, tlen: i64) -> Record {
        // Compute query-consuming length from CIGAR (M, I, S, X, = consume query)
        let seq_len: u32 = cigar.0.iter().map(|op| match op {
            Cigar::Match(n) | Cigar::Ins(n) | Cigar::SoftClip(n)
            | Cigar::Equal(n) | Cigar::Diff(n) => *n,
            _ => 0,
        }).sum();

        let seq: Vec<u8> = vec![b'A'; seq_len as usize];
        let qual: Vec<u8> = vec![255u8; seq_len as usize];

        let mut rec = Record::new();
        rec.set(b"r1", Some(cigar), &seq, &qual);
        rec.set_insert_size(tlen);
        rec
    }

    #[test]
    fn test_physical_size_zero_tlen() {
        let cigar = CigarString(vec![Cigar::Match(150)]);
        let rec = mock_record(&cigar, 0);
        assert_eq!(calc_physical_insert_size(&rec), 0, "TLEN=0 should return 0");
    }

    #[test]
    fn test_physical_size_no_indels() {
        // 150M, TLEN=167 → physical = 167 (no correction)
        let cigar = CigarString(vec![Cigar::Match(150)]);
        let rec = mock_record(&cigar, 167);
        assert_eq!(calc_physical_insert_size(&rec), 167);
    }

    #[test]
    fn test_physical_size_negative_tlen_no_indels() {
        // 150M, TLEN=-167 → physical = 167 (abs value)
        let cigar = CigarString(vec![Cigar::Match(150)]);
        let rec = mock_record(&cigar, -167);
        assert_eq!(calc_physical_insert_size(&rec), 167);
    }

    #[test]
    fn test_physical_size_with_deletion() {
        // 89M15D1M, TLEN=182 → physical = 182 - 15 = 167
        let cigar = CigarString(vec![Cigar::Match(89), Cigar::Del(15), Cigar::Match(1)]);
        let rec = mock_record(&cigar, 182);
        assert_eq!(calc_physical_insert_size(&rec), 167,
            "15bp deletion should subtract 15 from TLEN");
    }

    #[test]
    fn test_physical_size_with_large_deletion() {
        // 80M35D11M, TLEN=340 → physical = 340 - 35 = 305
        let cigar = CigarString(vec![Cigar::Match(80), Cigar::Del(35), Cigar::Match(11)]);
        let rec = mock_record(&cigar, 340);
        assert_eq!(calc_physical_insert_size(&rec), 305,
            "35bp deletion should subtract 35 from TLEN");
    }

    #[test]
    fn test_physical_size_with_insertion() {
        // 62M12I17M, TLEN=140 → physical = 140 + 12 = 152
        let cigar = CigarString(vec![Cigar::Match(62), Cigar::Ins(12), Cigar::Match(17)]);
        let rec = mock_record(&cigar, 140);
        assert_eq!(calc_physical_insert_size(&rec), 152,
            "12bp insertion should add 12 to TLEN");
    }

    #[test]
    fn test_physical_size_with_small_insertion() {
        // 83M6I2M, TLEN=119 → physical = 119 + 6 = 125
        let cigar = CigarString(vec![Cigar::Match(83), Cigar::Ins(6), Cigar::Match(2)]);
        let rec = mock_record(&cigar, 119);
        assert_eq!(calc_physical_insert_size(&rec), 125,
            "6bp insertion should add 6 to TLEN");
    }

    #[test]
    fn test_physical_size_combined_del_and_ins() {
        // 50M10D20M5I10M, TLEN=200 → physical = 200 - 10 + 5 = 195
        let cigar = CigarString(vec![
            Cigar::Match(50), Cigar::Del(10), Cigar::Match(20),
            Cigar::Ins(5), Cigar::Match(10),
        ]);
        let rec = mock_record(&cigar, 200);
        assert_eq!(calc_physical_insert_size(&rec), 195,
            "Combined D=10 I=5: 200 - 10 + 5 = 195");
    }

    #[test]
    fn test_physical_size_softclips_ignored() {
        // 5S85M15D1M5S, TLEN=182 → physical = 182 - 15 = 167
        // Soft-clips should NOT be included in the correction
        let cigar = CigarString(vec![
            Cigar::SoftClip(5), Cigar::Match(85), Cigar::Del(15),
            Cigar::Match(1), Cigar::SoftClip(5),
        ]);
        let rec = mock_record(&cigar, 182);
        assert_eq!(calc_physical_insert_size(&rec), 167,
            "Soft-clips should be ignored in physical sizing");
    }
}
