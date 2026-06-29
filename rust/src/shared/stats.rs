//! Statistical tests shared across analysis modes.
//!
//! Provides:
//! - [`fisher_exact_2x2`]: Fisher's exact test for any 2×2 contingency table
//!   (strand bias, ASJD junction comparison).
//! - [`fisher_exact_2x2_py`]: Python-callable wrapper (exposed via `gbcms._rs`).
//! - [`benjamini_hochberg`]: Benjamini-Hochberg FDR correction for multiple
//!   hypothesis testing (ASJD q-values).
//!
//! Used by:
//! - `counting/engine.rs` — Fisher strand bias per variant
//! - `counting/engine.rs` — ASJD junction comparison
//! - `counting/engine.rs` — BH correction across variants (ASJD)
//! - `merge.py` — combined strand bias on merged simplex+duplex counts

use pyo3::prelude::*;
use statrs::distribution::{Discrete, Hypergeometric};

/// Calculate Fisher's Exact Test for a 2×2 contingency table.
///
/// Returns (p-value, odds ratio).
///
/// Generic 2×2 table layout:
/// ```text
///     [[a, b],
///      [c, d]]
/// ```
///
/// **Strand bias usage** (existing):
/// ```text
///     [[ref_fwd, ref_rev],
///      [alt_fwd, alt_rev]]
/// ```
///
/// **ASJD usage**:
/// ```text
///     [[ref_junction_A, ref_junction_B],
///      [alt_junction_A, alt_junction_B]]
/// ```
///
/// Renamed from `fisher_strand_bias` in v5.0.0 for generality.
pub fn fisher_exact_2x2(a: u32, b: u32, c: u32, d: u32) -> (f64, f64) {
    let a = a as u64;
    let b = b as u64;
    let c = c as u64;
    let d = d as u64;

    let n = a + b + c + d;
    if n == 0 {
        return (1.0, f64::NAN);
    }

    // Guard: strand bias is undefined/underpowered with ≤1 observation in
    // either row. For the primary use-case (strand bias), row 2 = ALT counts.
    // With 0–1 ALT reads there is no statistical power to detect asymmetry;
    // the Hypergeometric distribution becomes degenerate (K ≈ N), causing
    // floating-point underflow that produces p ≈ 0 (see GitHub issue #19).
    //
    // Returning NaN for OR signals "undefined" — downstream writers format
    // NaN as '.' in VCF (spec-compliant missing value).
    let row2_total = c + d;
    if row2_total <= 1 {
        return (1.0, f64::NAN);
    }

    // Calculate Odds Ratio: (a*d) / (b*c)
    // When either cell in the denominator is zero, the OR is mathematically
    // undefined (0/0 or x/0). We return NaN rather than Infinity because:
    //   1. VCF Float fields do not support 'inf' (spec violation)
    //   2. NaN → '.' in output (standard missing value sentinel)
    let numerator = (a as f64) * (d as f64);
    let denominator = (b as f64) * (c as f64);
    let odds_ratio = if denominator == 0.0 {
        f64::NAN
    } else {
        numerator / denominator
    };

    // Fisher's Exact Test (Two-sided) using Hypergeometric distribution
    // We want the probability of observing a table as extreme or more extreme than the current one,
    // given fixed marginals.
    // Hypergeometric(N, K, n) where:
    // N = total population size (a+b+c+d)
    // K = number of successes in population (a+b) (Row 1 sum)
    // n = sample size (a+c) (Col 1 sum)
    // k = number of successes in sample (a)

    let row1_sum = a + b;
    let col1_sum = a + c;

    // If any marginal is 0, p-value is 1.0
    if row1_sum == 0 || col1_sum == 0 || row1_sum == n || col1_sum == n {
        return (1.0, odds_ratio);
    }

    let dist = match Hypergeometric::new(n, row1_sum, col1_sum) {
        Ok(d) => d,
        Err(_) => return (1.0, odds_ratio), // Should not happen with checks above
    };

    let p_observed = dist.pmf(a);
    let mut p_value = 0.0;

    // Sum probabilities of all tables with p <= p_observed
    // Range of possible values for cell 'a' is [max(0, row1_sum + col1_sum - n), min(row1_sum, col1_sum)]
    let min_a = (row1_sum + col1_sum).saturating_sub(n);
    let max_a = if row1_sum < col1_sum {
        row1_sum
    } else {
        col1_sum
    };

    for k in min_a..=max_a {
        let p = dist.pmf(k);
        if p <= p_observed + 1e-10 {
            // Add epsilon for float comparison
            p_value += p;
        }
    }

    // Cap at 1.0
    if p_value > 1.0 {
        p_value = 1.0;
    }

    (p_value, odds_ratio)
}

/// Backward-compatible alias for existing callers.
///
/// Delegates to [`fisher_exact_2x2`]. Existing code can continue to use this
/// name without modification.
#[inline]
pub fn fisher_strand_bias(ref_fwd: u32, ref_rev: u32, alt_fwd: u32, alt_rev: u32) -> (f64, f64) {
    fisher_exact_2x2(ref_fwd, ref_rev, alt_fwd, alt_rev)
}

/// Python-accessible Fisher's exact test for a 2×2 contingency table.
///
/// Returns ``(p_value, odds_ratio)``.
///
/// This is a thin wrapper around [`fisher_exact_2x2`] exposed to Python
/// via PyO3. Used by the merge engine (``merge.py``) to compute combined
/// strand bias statistics on summed simplex+duplex forward/reverse counts.
///
/// # Example (Python)
///
/// ```python
/// from gbcms._rs import fisher_exact_2x2
/// p_value, odds_ratio = fisher_exact_2x2(10, 12, 8, 15)
/// ```
#[pyfunction]
#[pyo3(name = "fisher_exact_2x2")]
pub fn fisher_exact_2x2_py(a: u32, b: u32, c: u32, d: u32) -> (f64, f64) {
    fisher_exact_2x2(a, b, c, d)
}

/// Benjamini-Hochberg FDR correction for multiple hypothesis testing.
///
/// Given a slice of raw p-values, returns a Vec of adjusted q-values.
/// The q-values are monotonically non-decreasing and capped at 1.0.
///
/// # Algorithm
///
/// 1. Sort p-values in descending order (by index).
/// 2. Walk from largest to smallest: q_i = min(p_i * n / rank, q_{i+1}).
/// 3. Cap all values at 1.0.
///
/// This matches R's `p.adjust(method = "BH")` exactly.
///
/// # Parameters
///
/// - `pvalues`: raw p-values (must be in [0, 1]).
///
/// # Returns
///
/// Vec of adjusted q-values, same length and order as input.
///
/// # Example
///
/// ```
/// # use _rs::shared::stats::benjamini_hochberg;
/// let pvals = vec![0.01, 0.04, 0.03, 0.10, 0.50];
/// let qvals = benjamini_hochberg(&pvals);
/// // Matches R: p.adjust(c(0.01, 0.04, 0.03, 0.10, 0.50), method="BH")
/// // → [0.05, 0.0667, 0.05, 0.125, 0.50]
/// ```
pub fn benjamini_hochberg(pvalues: &[f64]) -> Vec<f64> {
    let n = pvalues.len();
    if n == 0 {
        return vec![];
    }

    // Sort indices by descending p-value
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| {
        pvalues[b]
            .partial_cmp(&pvalues[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut qvalues = vec![0.0; n];
    let mut cummin = f64::INFINITY;

    for (rank_from_end, &idx) in indices.iter().enumerate() {
        // rank_from_end: 0 = largest p-value, n-1 = smallest
        let rank = n - rank_from_end; // 1-based rank from smallest
        let adjusted = (pvalues[idx] * n as f64) / rank as f64;
        cummin = cummin.min(adjusted).min(1.0);
        qvalues[idx] = cummin;
    }

    qvalues
}

/// Benjamini-Hochberg FDR over a *sub-family* selected by index.
///
/// Corrects only the p-values at `valid_indices`, so the family size `n` equals
/// `valid_indices.len()` — not the length of the full result vector. Each corrected
/// entry is returned as `(index, q-value)`, ready to scatter back into the caller's
/// per-variant records.
///
/// This is the multiplicity guard for the mFSD and ASJD q-values: variants with no
/// real test (too few fragments → placeholder p = 1.0, or no junction reads) are
/// excluded by the caller's predicate and never enter the family, so they cannot
/// pad `n` and over-correct the variants that *were* tested. `pvalue_at(i)` reads the
/// raw p-value for index `i` from the caller's records without allocating a full
/// p-value vector.
///
/// # Example
///
/// ```
/// # use _rs::shared::stats::benjamini_hochberg_family;
/// // Indices 1 and 3 are no-test padding; only 0, 2, 4 form the family (n = 3).
/// let pvals = [0.01, 1.0, 0.04, 1.0, 0.50];
/// let corrected = benjamini_hochberg_family(|i| pvals[i], &[0, 2, 4]);
/// assert_eq!(corrected.len(), 3);
/// ```
pub fn benjamini_hochberg_family(
    pvalue_at: impl Fn(usize) -> f64,
    valid_indices: &[usize],
) -> Vec<(usize, f64)> {
    let sub: Vec<f64> = valid_indices.iter().map(|&i| pvalue_at(i)).collect();
    benjamini_hochberg(&sub)
        .into_iter()
        .zip(valid_indices.iter().copied())
        .map(|(q, i)| (i, q))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── fisher_exact_2x2 tests ──

    #[test]
    fn test_fisher_no_bias() {
        // Balanced table: no strand bias
        let (p, _or) = fisher_exact_2x2(10, 10, 10, 10);
        assert!(p > 0.9, "Balanced table should have p ~1.0, got {}", p);
    }

    #[test]
    fn test_fisher_strong_bias() {
        // Extreme bias: all ref on fwd, all alt on rev
        let (p, _or) = fisher_exact_2x2(20, 0, 0, 20);
        assert!(p < 0.001, "Extreme bias should have p < 0.001, got {}", p);
    }

    #[test]
    fn test_fisher_empty_table() {
        let (p, or) = fisher_exact_2x2(0, 0, 0, 0);
        assert_eq!(p, 1.0);
        assert!(or.is_nan(), "Empty table OR should be NaN, got {}", or);
    }

    #[test]
    fn test_fisher_single_alt_read_runx1() {
        // Regression: RUNX1 duplex BAM variant with 1 ALT read (GitHub #19).
        // Previously returned (p≈0.0, OR=inf) due to degenerate Hypergeometric.
        // With ≤1 ALT read, strand bias is undefined — must return (1.0, NaN).
        let (p, or) = fisher_exact_2x2(1852, 1484, 0, 1);
        assert_eq!(p, 1.0, "1 ALT read should give p=1.0, got {}", p);
        assert!(or.is_nan(), "1 ALT read should give OR=NaN, got {}", or);
    }

    #[test]
    fn test_fisher_zero_alt_reads() {
        // No ALT reads at all — strand bias is undefined.
        let (p, or) = fisher_exact_2x2(100, 100, 0, 0);
        assert_eq!(p, 1.0, "0 ALT reads should give p=1.0, got {}", p);
        assert!(or.is_nan(), "0 ALT reads should give OR=NaN, got {}", or);
    }

    #[test]
    fn test_fisher_single_alt_one_sided() {
        // 1 ALT read on forward only — still ≤1 total ALT, so underpowered.
        let (p, or) = fisher_exact_2x2(50, 50, 1, 0);
        assert_eq!(p, 1.0, "1 ALT read (one-sided) should give p=1.0, got {}", p);
        assert!(or.is_nan(), "1 ALT read (one-sided) should give OR=NaN, got {}", or);
    }

    #[test]
    fn test_fisher_two_alt_reads_computes_normally() {
        // Boundary: exactly 2 ALT reads — should compute normally (not early return).
        // With 2 ALT reads (1 fwd, 1 rev), balanced → p should be high (no bias).
        let (p, or) = fisher_exact_2x2(50, 50, 1, 1);
        assert!(p > 0.5, "2 balanced ALT reads should have p > 0.5, got {}", p);
        assert!(!or.is_nan(), "2 ALT reads should produce a real OR, got NaN");
    }

    #[test]
    fn test_fisher_backward_compat() {
        // Verify the alias produces identical results
        let (p1, or1) = fisher_exact_2x2(5, 15, 12, 3);
        let (p2, or2) = fisher_strand_bias(5, 15, 12, 3);
        assert_eq!(p1, p2);
        assert_eq!(or1, or2);
    }

    // ── benjamini_hochberg tests ──

    #[test]
    fn test_bh_empty() {
        let q = benjamini_hochberg(&[]);
        assert!(q.is_empty());
    }

    #[test]
    fn test_bh_single() {
        let q = benjamini_hochberg(&[0.05]);
        assert_eq!(q, vec![0.05]);
    }

    #[test]
    fn test_bh_matches_r() {
        // Standard BH procedure on p = [0.01, 0.04, 0.03, 0.10, 0.50]:
        // Ascending ranks: p[0]=0.01 (rank 1), p[2]=0.03 (rank 2), p[1]=0.04 (rank 3),
        //                  p[3]=0.10 (rank 4), p[4]=0.50 (rank 5)
        // Adjusted: [0.05, 0.0667, 0.075, 0.125, 0.50]
        // After cummin (descending): [0.05, 0.0667, 0.0667, 0.125, 0.50]
        let pvals = vec![0.01, 0.04, 0.03, 0.10, 0.50];
        let qvals = benjamini_hochberg(&pvals);

        let expected = [0.05, 1.0 / 15.0, 1.0 / 15.0, 0.125, 0.50];
        for (i, (&q, &e)) in qvals.iter().zip(expected.iter()).enumerate() {
            assert!(
                (q - e).abs() < 1e-10,
                "q[{}] = {}, expected {}",
                i,
                q,
                e
            );
        }
    }

    #[test]
    fn test_bh_preserves_order() {
        // All same p-values → all same q-values
        let pvals = vec![0.05; 10];
        let qvals = benjamini_hochberg(&pvals);
        for q in &qvals {
            assert!((q - 0.05).abs() < 1e-10, "q = {}", q);
        }
    }

    #[test]
    fn test_bh_caps_at_one() {
        // Large p-values should be capped at 1.0
        let pvals = vec![0.9, 0.95, 0.99];
        let qvals = benjamini_hochberg(&pvals);
        for q in &qvals {
            assert!(*q <= 1.0, "q should be <= 1.0, got {}", q);
        }
    }

    #[test]
    fn test_bh_monotonic() {
        // Sorted p-values should produce monotonically non-decreasing q-values
        let pvals = vec![0.001, 0.01, 0.05, 0.1, 0.5];
        let qvals = benjamini_hochberg(&pvals);
        for i in 1..qvals.len() {
            assert!(
                qvals[i] >= qvals[i - 1] - 1e-10,
                "q[{}]={} < q[{}]={}",
                i,
                qvals[i],
                i - 1,
                qvals[i - 1]
            );
        }
    }

    #[test]
    fn test_bh_family_excludes_padding() {
        // Indices 1 and 3 are "no-test" padding (placeholder p = 1.0). Only the
        // genuine tests (indices 0, 2, 4) may form the family, so n = 3 — the
        // padding must NOT inflate it to 5.
        let pvals = [0.01, 1.0, 0.04, 1.0, 0.50];
        let valid = [0usize, 2, 4];
        let corrected = benjamini_hochberg_family(|i| pvals[i], &valid);

        // q-values must match BH over ONLY the 3 real p-values (n = 3),
        // not BH over all 5.
        let expected = benjamini_hochberg(&[0.01, 0.04, 0.50]);
        let inflated = benjamini_hochberg(&pvals);
        assert_eq!(corrected.len(), 3);
        for ((idx, q), (&vidx, &eq)) in corrected.iter().zip(valid.iter().zip(expected.iter())) {
            assert_eq!(*idx, vidx);
            assert!((q - eq).abs() < 1e-12, "index {idx}: {q} vs n=3 BH {eq}");
            // And strictly tighter than the padded n=5 correction would give.
            assert!(*q <= inflated[vidx] + 1e-12);
        }
        // The padding indices never appear in the corrected output.
        assert!(corrected.iter().all(|(i, _)| *i != 1 && *i != 3));
    }

    #[test]
    fn test_bh_family_empty_is_noop() {
        // No valid tests → nothing to correct, no panic.
        let corrected = benjamini_hochberg_family(|_| 0.5, &[]);
        assert!(corrected.is_empty());
    }
}
