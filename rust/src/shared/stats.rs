//! Statistical tests shared across analysis modes.
//!
//! Currently provides Fisher's exact test for strand bias (used by counting).
//! Will also host Beta-Binomial LOH test (used by HLA) in a future release.
//!
//! Used by:
//! - `counting/engine.rs` — Fisher strand bias per variant
//! - `hla/aggregate.rs` (future) — Beta-Binomial LOH p-value

use statrs::distribution::{Discrete, Hypergeometric};

/// Calculate Fisher's Exact Test for strand bias.
///
/// Returns (p-value, odds ratio).
///
/// Contingency table:
/// [[ref_fwd, ref_rev],
///  [alt_fwd, alt_rev]]
pub fn fisher_strand_bias(ref_fwd: u32, ref_rev: u32, alt_fwd: u32, alt_rev: u32) -> (f64, f64) {
    let a = ref_fwd as u64;
    let b = ref_rev as u64;
    let c = alt_fwd as u64;
    let d = alt_rev as u64;

    let n = a + b + c + d;
    if n == 0 {
        return (1.0, 0.0);
    }

    // Calculate Odds Ratio: (a*d) / (b*c)
    // Add small epsilon to avoid division by zero if needed, but standard is 0/Inf
    let numerator = (a as f64) * (d as f64);
    let denominator = (b as f64) * (c as f64);
    let odds_ratio = if denominator == 0.0 {
        if numerator == 0.0 { 0.0 } else { f64::INFINITY }
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
