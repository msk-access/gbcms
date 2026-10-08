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
//! - The LLR scores each fragment class relative to a Gaussian tumor/healthy model,
//!   reported as the mean per fragment (comparable across classes of any size)
//!
//! All functions operate on raw, unweighted fragment size slices. GC correction
//! is not applied here — GC bias affects count depth, not fragment length, so
//! distributional tests are already unbiased on observed sizes.
//!
//! ## Usage
//! ```ignore
//! let physical = mfsd::calc_physical_insert_size(&record);
//! let (ks_d, ks_p) = mfsd::ks_test(&alt_sizes, &ref_sizes);
//! let llr = mfsd::calc_llr_mean(&alt_sizes);
//! let mean = mfsd::calc_mean(&alt_sizes);
//! ```

use log::trace;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;

/// Minimum number of fragments required in each class for the KS test.
/// Below this threshold, `ks_test` returns `(f64::NAN, 1.0)` to signal
/// insufficient data rather than a spurious result.
pub const MIN_FOR_KS: usize = 5;

/// Largest lattice (n·m cells) for which the KS p-value is computed exactly. Every
/// class pair a targeted cfDNA panel produces fits (a few thousand fragments a
/// side); beyond it both classes are large and the corrected asymptotic series is
/// accurate (within a few percent once the smaller class has a few hundred).
const KS_EXACT_MAX_CELLS: u64 = 10_000_000;

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

// (gaussian_pdf removed — the LLR now uses the closed-form Gaussian
// log-ratio directly, so the density function and its underflow are gone.)

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
    lengths
        .iter()
        .map(|&x| {
            // Closed-form Gaussian log-ratio: ln(P_tumor(x)) − ln(P_healthy(x)).
            // This is finite for every finite x. The previous code formed the
            // pdf ratio p_tumor/p_healthy and took its log, which underflowed to
            // ±∞ in the tails (p_healthy < f64::EPSILON → +∞; a symmetric p_tumor
            // underflow → ln(0) = −∞), so a single tail fragment pinned the whole
            // sum to ±Infinity. The 0.5·ln(2π) normalisers cancel in the ratio,
            // leaving 0.5·(z_h² − z_t²) + ln(σ_h/σ_t).
            let z_t = (x - params.tumor_mu) / params.tumor_sigma;
            let z_h = (x - params.healthy_mu) / params.healthy_sigma;
            0.5 * (z_h * z_h - z_t * z_t) + (params.healthy_sigma / params.tumor_sigma).ln()
        })
        .sum()
}

/// Mean per-fragment log-likelihood ratio: [`calc_llr`] divided by the number of
/// fragments, so it does not grow with depth (the sum grows about 22x from n < 20
/// to n >= 100 on real cfDNA, while the mean stays flat). NaN for an empty class:
/// with no fragments there is no ratio to report.
pub fn calc_llr_mean(lengths: &[f64]) -> f64 {
    if lengths.is_empty() {
        return f64::NAN;
    }
    calc_llr(lengths) / lengths.len() as f64
}

// ── Physical Fragment Sizing ─────────────────────────────────────────────────

/// Compute the physical fragment insert size from CIGAR, correcting TLEN for indels.
///
/// BAM TLEN measures the reference span between the outermost aligned bases of a
/// read pair. For fragments carrying indels, TLEN ≠ physical fragment length:
/// - A deletion makes TLEN *longer* than the physical DNA fragment
/// - An insertion makes TLEN *shorter* than the physical DNA fragment
///
/// **Formula:** `physical_size = |TLEN| - sum(D_ops) - sum(N_ops) + sum(I_ops)`
///
/// `N` (RefSkip) ops are introns in spliced RNA reads: TLEN spans the genomic
/// distance *including* the intron, but the physical (mature-mRNA) fragment does
/// not, so introns are discounted like deletions. DNA reads carry no `N` ops, so
/// the term is zero there and the cfDNA correction is unchanged.
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
    let mut skip_bp: i32 = 0;

    for op in record.cigar().iter() {
        match op {
            Cigar::Del(n) => del_bp += *n as i32,
            Cigar::Ins(n) => ins_bp += *n as i32,
            Cigar::RefSkip(n) => skip_bp += *n as i32,
            _ => {}
        }
    }

    let physical = abs_tlen - del_bp - skip_bp + ins_bp;

    // Guard: pathological CIGAR where correction overshoots (should never happen
    // with well-formed BAMs, but avoid returning nonsensical negative sizes)
    if physical <= 0 {
        trace!(
            "calc_physical_insert_size: pathological correction — TLEN={} D={} N={} I={} → {} (clamped to |TLEN|)",
            raw_tlen, del_bp, skip_bp, ins_bp, physical
        );
        return abs_tlen; // Fall back to raw |TLEN| rather than a nonsensical value
    }

    trace!(
        "calc_physical_insert_size: TLEN={} D={} N={} I={} → physical={}",
        raw_tlen, del_bp, skip_bp, ins_bp, physical
    );

    physical
}

// ── Two-Sample Kolmogorov-Smirnov Test ───────────────────────────────────────

/// Two-sample Kolmogorov-Smirnov test.
///
/// Computes the KS D-statistic (maximum absolute difference between empirical
/// CDFs, evaluated after each block of tied values) and its p-value: exact up to
/// `KS_EXACT_MAX_CELLS` lattice cells, else the finite-sample-corrected asymptotic
/// series (see [`ks_p_value`]).
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
    let (n, m) = (a_sorted.len(), b_sorted.len());

    // Merge-walk the two ECDFs, taking the gap after each block of tied values. The
    // gap is kept in integers scaled by n·m (|i·m − j·n|), so the exact p-value can
    // compare lattice points against the observed deviation without float error.
    let mut dev: i64 = 0;
    let (mut i, mut j) = (0usize, 0usize);
    while i < n && j < m {
        let val = if a_sorted[i] <= b_sorted[j] { a_sorted[i] } else { b_sorted[j] };
        // Advance all entries equal to `val` in both arrays
        while i < n && a_sorted[i] <= val { i += 1; }
        while j < m && b_sorted[j] <= val { j += 1; }
        dev = dev.max((i as i64 * m as i64 - j as i64 * n as i64).abs());
    }

    let d = dev as f64 / (n as f64 * m as f64);
    (d, ks_p_value(dev, d, n, m))
}

/// Two-sample KS p-value: P(D ≥ d) under the null, for an observed deviation `dev`
/// = D·n·m (an integer; `d` is the same value as a fraction).
///
/// Exact whenever the lattice has at most `KS_EXACT_MAX_CELLS` cells. A cfDNA
/// variant usually pairs a few ALT fragments with thousands of REF fragments; the
/// asymptotic Kolmogorov series overstates p there (about 1.7x at 5 ALT fragments
/// near p = 0.05, and 45x for a strong shift), and the old switch on `n·m` (R's
/// `ks.test` rule) sent most real ALT-vs-REF pairs to it. Beyond the limit the
/// lattice is large (for realistic shapes, both classes are), and the series with
/// Stephens' correction is used.
///
/// The exact p-value treats the sizes as continuous; with tied (integer) sizes it
/// is conservative (measured on real cfDNA: 3.7–4.7% of null draws below 0.05).
fn ks_p_value(dev: i64, d: f64, n: usize, m: usize) -> f64 {
    if (n as u64) * (m as u64) <= KS_EXACT_MAX_CELLS {
        ks_p_value_exact(dev, n, m)
    } else {
        ks_p_value_asymptotic(d, n as f64, m as f64)
    }
}

/// Exact two-sample KS p-value from the monotone lattice paths (Hodges 1957).
///
/// A path from (0,0) to (n,m) takes one step per sorted observation; under the null
/// every path is equally likely. A lattice point is inside the band when its
/// deviation |i·m − j·n| is below the observed `dev` (integers, so a point exactly
/// at the observed deviation counts as reaching D: `P(D ≥ d)`). `v(i, j)` is the
/// share of the paths reaching (i, j) that have already left the band — 1 at a point
/// outside it, and inside it `v(i,j) = v(i−1,j)·i/(i+j) + v(i,j−1)·j/(i+j)`, since
/// those are the shares of paths arriving from each neighbour. The p-value is
/// `v(n, m)`, computed directly (no `1 − u` cancellation, so a tiny p keeps its
/// digits), with every value in [0, 1] (no path count or binomial can overflow).
/// Only the cells inside the band are visited; the rest stay at 1.
fn ks_p_value_exact(dev: i64, n: usize, m: usize) -> f64 {
    if dev <= 0 {
        return 1.0;
    }
    // Rows over the larger class, columns over the smaller: the band is symmetric in
    // (n, m), and the working row is the smaller one.
    let (n, m) = if n >= m { (n, m) } else { (m, n) };
    let (ni, mi) = (n as i64, m as i64);
    // Columns j inside the band on row i: |i·m − j·n| < dev, i.e.
    // (i·m − dev)/n < j < (i·m + dev)/n.
    let band = |i: usize| -> (usize, usize) {
        let c = i as i64 * mi;
        let lo = (c - dev).div_euclid(ni) + 1;
        let hi = (c + dev - 1).div_euclid(ni);
        (lo.max(0) as usize, hi.min(mi) as usize)
    };

    let mut row = vec![1.0f64; m + 1]; // outside the band: every path has left it
    let (lo0, hi0) = band(0);
    for v in row.iter_mut().take(hi0 + 1).skip(lo0) {
        *v = 0.0; // row 0 inside the band: reached only along the edge, never outside
    }
    let mut prev_lo = lo0;
    for i in 1..=n {
        let (lo, hi) = band(i);
        // Cells that left the band on the left since the previous row read as 1.
        for v in row.iter_mut().take(lo.min(m + 1)).skip(prev_lo) {
            *v = 1.0;
        }
        if lo <= hi {
            let fi = i as f64;
            for j in lo..=hi {
                let up = row[j]; // v(i−1, j): 1 if it was outside the previous band
                row[j] = if j == 0 {
                    up
                } else {
                    let k = fi + j as f64;
                    (up * fi + row[j - 1] * j as f64) / k
                };
            }
        }
        prev_lo = lo;
    }
    row[m].clamp(0.0, 1.0)
}

/// Asymptotic KS p-value via the Kolmogorov series Q_KS(λ) = 2 Σ (−1)^(k−1) e^(−2k²λ²)
/// with Stephens' finite-sample correction, λ = (√Nₑ + 0.12 + 0.11/√Nₑ)·D,
/// Nₑ = n·m/(n+m). Used only above `KS_EXACT_MAX_CELLS`.
fn ks_p_value_asymptotic(d: f64, n: f64, m: f64) -> f64 {
    let sqrt_ne = (n * m / (n + m)).sqrt();
    let lambda = (sqrt_ne + 0.12 + 0.11 / sqrt_ne) * d;
    if lambda < f64::EPSILON {
        return 1.0;
    }
    let mut sum = 0.0;
    for k in 1..=100i64 {
        let term = (-2.0 * (k as f64 * lambda).powi(2)).exp();
        sum += if k % 2 == 0 { -term } else { term };
        if term < 1e-12 {
            break;
        }
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

    // ─── LLR finiteness + exact KS ────────────────────────────────────────────

    #[test]
    fn test_llr_finite_in_tails() {
        // Extreme fragment sizes used to underflow the pdf ratio to ±∞ and pin the
        // whole sum to ±Infinity. The closed-form log-ratio is always finite.
        let llr = calc_llr(&[10.0, 300.0, 600.0, 1000.0]);
        assert!(llr.is_finite(), "LLR must be finite for tail fragments, got {llr}");
    }

    // Reference values: SciPy ks_2samp(method='exact'). Baked in as literals so the
    // test has no scipy dependency (the n=5 disjoint case is also self-evident: 2/C(10,5)).
    #[test]
    fn test_ks_exact_disjoint_n5() {
        let a = vec![100.0, 110.0, 120.0, 130.0, 140.0];
        let b = vec![160.0, 170.0, 180.0, 190.0, 200.0];
        let (d, p) = ks_test(&a, &b);
        assert!((d - 1.0).abs() < 1e-9, "D={d}");
        assert!((p - 2.0 / 252.0).abs() < 1e-6, "exact p={p}, want {}", 2.0 / 252.0);
    }

    #[test]
    fn test_ks_exact_overlap_n5() {
        let a = vec![100.0, 105.0, 110.0, 115.0, 120.0];
        let b = vec![112.0, 118.0, 122.0, 128.0, 132.0];
        let (d, p) = ks_test(&a, &b);
        assert!((d - 0.6).abs() < 1e-9, "D={d}");
        assert!((p - 0.357143).abs() < 1e-5, "exact p={p}, want 0.357143");
    }

    #[test]
    fn test_ks_exact_n8() {
        let a = vec![140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0];
        let b = vec![150.0, 152.0, 154.0, 156.0, 158.0, 160.0, 162.0, 164.0];
        let (d, p) = ks_test(&a, &b);
        assert!((d - 0.375).abs() < 1e-9, "D={d}");
        assert!((p - 0.660140).abs() < 1e-5, "exact p={p}, want 0.660140");
    }

    // ─── mean LLR + exact KS at any depth ─────────────────────────────────────

    #[test]
    fn test_llr_mean_is_the_sum_over_n() {
        let sizes = [118.0, 126.0, 133.0, 141.0, 152.0];
        assert!((calc_llr_mean(&sizes) - calc_llr(&sizes) / 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_llr_mean_of_an_empty_class_is_nan() {
        assert!(calc_llr_mean(&[]).is_nan());
    }

    // Few ALT fragments against a deep REF class: 5 x 2,400 = 12,000 lattice cells,
    // past R's exact limit. Reference: the share-of-paths recursion in exact Python
    // arithmetic (tests/test_statistics_group6_contract.py::_ks_exact); the
    // asymptotic series gave 0.0346 here.
    #[test]
    fn test_ks_exact_few_alt_deep_ref() {
        let b: Vec<f64> = (0..2400).map(|k| (140 + (k * 37) % 80) as f64).collect();
        let a = vec![143.0, 150.0, 156.0, 161.0, 168.0];
        let (d, p) = ks_test(&a, &b);
        assert!((d - 0.6375).abs() < 1e-12, "D={d}");
        assert!((p - 0.017188813774322798).abs() < 1e-9, "exact p={p}");
    }

    #[test]
    fn test_ks_exact_is_symmetric_in_the_classes() {
        let a: Vec<f64> = (0..7).map(|k| (120 + 9 * k) as f64).collect();
        let b: Vec<f64> = (0..300).map(|k| (130 + (k * 13) % 70) as f64).collect();
        let (d1, p1) = ks_test(&a, &b);
        let (d2, p2) = ks_test(&b, &a);
        assert!((d1 - d2).abs() < 1e-12 && (p1 - p2).abs() < 1e-12, "{p1} vs {p2}");
    }

    #[test]
    fn test_ks_large_classes_use_the_corrected_series() {
        // 4,000 x 4,000 cells is past the exact limit: the p-value is the Kolmogorov
        // series at Stephens' corrected lambda (the uncorrected lambda gives another value).
        let a: Vec<f64> = (0..4000).map(|k| (120 + (k * 7) % 90) as f64).collect();
        let b: Vec<f64> = (0..4000).map(|k| (123 + (k * 11) % 90) as f64).collect();
        let (d, p) = ks_test(&a, &b);
        let series = |lambda: f64| -> f64 {
            let s: f64 = (1..=100i64)
                .map(|k| { let t = (-2.0 * (k as f64 * lambda).powi(2)).exp(); if k % 2 == 0 { -t } else { t } })
                .sum();
            (2.0 * s).clamp(0.0, 1.0)
        };
        let sqrt_ne = (4000.0f64 * 4000.0 / 8000.0).sqrt();
        let corrected = series((sqrt_ne + 0.12 + 0.11 / sqrt_ne) * d);
        let uncorrected = series(sqrt_ne * d);
        assert!(p > 1e-4 && p < 0.9, "pick data with a mid-range p: {p}");
        assert!((p - corrected).abs() < 1e-12, "p={p} corrected={corrected}");
        assert!((corrected - uncorrected).abs() > 1e-6);
    }

    /// Sizes from a 32-bit LCG (Numerical Recipes constants), reproducible in tests.
    fn lcg_sizes(seed: u32, k: usize, lo: u32, span: u32) -> Vec<f64> {
        let mut x = seed;
        (0..k)
            .map(|_| {
                x = x.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                (lo + (x >> 16) % span) as f64
            })
            .collect()
    }

    /// Independent reference: P(D ≥ d_obs) over every placement of the a-values among
    /// the pooled sorted values (tie-free data; the continuous null), by enumeration.
    fn brute_force_p(a: &[f64], b: &[f64]) -> f64 {
        let (n, m) = (a.len(), b.len());
        let (d_obs, _) = ks_test(a, b);
        let mut pooled: Vec<f64> = a.iter().chain(b.iter()).copied().collect();
        pooled.sort_by(f64::total_cmp);
        let (mut hit, mut total) = (0u64, 0u64);
        for mask in 0u32..(1 << (n + m)) {
            if mask.count_ones() as usize != n {
                continue;
            }
            total += 1;
            let (mut i, mut j, mut d) = (0usize, 0usize, 0f64);
            for k in 0..n + m {
                if mask >> k & 1 == 1 { i += 1 } else { j += 1 }
                d = d.max((i as f64 / n as f64 - j as f64 / m as f64).abs());
            }
            let _ = &pooled;
            if d >= d_obs - 1e-12 {
                hit += 1;
            }
        }
        hit as f64 / total as f64
    }

    #[test]
    fn test_ks_exact_matches_brute_force_enumeration() {
        for (n, m) in [(5usize, 5usize), (5, 7), (6, 5), (5, 9)] {
            for seed in 1..6u32 {
                // distinct values (no ties), so enumeration and the lattice agree exactly
                let mut a: Vec<f64> = lcg_sizes(seed, n, 100, 1000);
                let mut b: Vec<f64> = lcg_sizes(seed + 77, m, 100, 1000);
                for (k, v) in a.iter_mut().enumerate() { *v += 0.001 * k as f64; }
                for (k, v) in b.iter_mut().enumerate() { *v += 0.0005 + 0.001 * k as f64; }
                let (_, p) = ks_test(&a, &b);
                let want = brute_force_p(&a, &b);
                assert!((p - want).abs() < 1e-12, "n={n} m={m} seed={seed}: p={p} want={want}");
            }
        }
    }

    /// Reference p-value over the full lattice: the observed deviation from the ECDF
    /// walk in integers, every cell visited, p = 1 − (share of paths kept inside).
    fn full_grid_p(a: &[f64], b: &[f64]) -> f64 {
        let (mut sa, mut sb) = (a.to_vec(), b.to_vec());
        sa.sort_by(f64::total_cmp);
        sb.sort_by(f64::total_cmp);
        let (n, m) = (sa.len() as i64, sb.len() as i64);
        let (mut i, mut j, mut dev) = (0usize, 0usize, 0i64);
        while i < sa.len() && j < sb.len() {
            let v = sa[i].min(sb[j]);
            while i < sa.len() && sa[i] <= v { i += 1; }
            while j < sb.len() && sb[j] <= v { j += 1; }
            dev = dev.max((i as i64 * m - j as i64 * n).abs());
        }
        let inside = |i: usize, j: usize| (i as i64 * m - j as i64 * n).abs() < dev;
        let mut row = vec![0f64; sb.len() + 1];
        row[0] = 1.0;
        for jj in 1..=sb.len() {
            row[jj] = if inside(0, jj) { row[jj - 1] } else { 0.0 };
        }
        for ii in 1..=sa.len() {
            row[0] = if inside(ii, 0) { row[0] } else { 0.0 };
            for jj in 1..=sb.len() {
                row[jj] = if inside(ii, jj) {
                    let k = (ii + jj) as f64;
                    row[jj] * ii as f64 / k + row[jj - 1] * jj as f64 / k
                } else {
                    0.0
                };
            }
        }
        (1.0 - row[sb.len()]).clamp(0.0, 1.0)
    }

    // The observed deviation must count as reaching D (P(D >= d), not P(D > d)) at any
    // lattice size. Here (2,583 x 3,871, just under the exact cap) a float band
    // d*n*m - 1e-9 lands above the integer deviation; the reference compares integers.
    #[test]
    fn test_ks_exact_counts_the_observed_point_in_integers() {
        let a = lcg_sizes(3369, 2583, 100, 161);
        let b = lcg_sizes(103_369, 3871, 104, 161);
        let (_, p) = ks_test(&a, &b);
        let want = full_grid_p(&a, &b);
        assert!((p - want).abs() <= 1e-12, "p={p} want={want}");
    }

    // The banded walk visits only cells inside the band; it must equal the full lattice
    // for any shape and either class order (rows run over the larger class).
    #[test]
    fn test_ks_banded_walk_equals_the_full_lattice() {
        for seed in 1..40u32 {
            let n = 5 + (seed as usize * 37) % 300;
            let m = 5 + (seed as usize * 53) % 400;
            let a = lcg_sizes(seed, n, 120, 40 + seed % 60);
            let b = lcg_sizes(seed + 500, m, 118 + seed % 9, 50);
            let (_, p) = ks_test(&a, &b);
            let want = full_grid_p(&a, &b);
            assert!((p - want).abs() <= 1e-11, "seed={seed} n={n} m={m}: p={p} want={want}");
        }
    }

    // Disjoint classes: P(D >= 1) = 2 / C(n+m, n) (all ALT below, or all above, all REF).
    // About 3e-15 here; computing p as 1 - u loses it to cancellation.
    #[test]
    fn test_ks_tiny_p_is_accurate() {
        let a: Vec<f64> = (0..5).map(|k| 100.0 + k as f64).collect();
        let b: Vec<f64> = (0..2400).map(|k| 200.0 + (k % 300) as f64).collect();
        let (d, p) = ks_test(&a, &b);
        assert!((d - 1.0).abs() < 1e-12);
        let comb: f64 = (0..5).map(|k| (2405 - k) as f64 / (k + 1) as f64).product();
        let want = 2.0 / comb;
        assert!(((p - want) / want).abs() < 1e-6, "p={p} want={want}");
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
    fn test_physical_size_spliced_read_discounts_intron() {
        // Spliced RNA read 50M2000N50M, TLEN spans the intron (100 + 2000 = 2100)
        // → physical mature-mRNA fragment = 2100 - 2000 = 100.
        let cigar = CigarString(vec![Cigar::Match(50), Cigar::RefSkip(2000), Cigar::Match(50)]);
        let rec = mock_record(&cigar, 2100);
        assert_eq!(calc_physical_insert_size(&rec), 100,
            "2000bp intron (N) should be discounted, not inflate the fragment size");
    }

    #[test]
    fn test_physical_size_spliced_with_indels() {
        // 40M500N5D40M3I, intron + 5bp del + 3bp ins; TLEN=588
        // → 588 - 500 (N) - 5 (D) + 3 (I) = 86.
        let cigar = CigarString(vec![
            Cigar::Match(40), Cigar::RefSkip(500), Cigar::Del(5),
            Cigar::Match(40), Cigar::Ins(3),
        ]);
        let rec = mock_record(&cigar, 588);
        assert_eq!(calc_physical_insert_size(&rec), 86,
            "intron, deletion and insertion corrections combine");
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
