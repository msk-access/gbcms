//! Shared utilities used by both `counting` and `hla` modules.
//!
//! This module provides reusable infrastructure that is not tied to any
//! specific analysis mode. Domain-specific logic stays in its respective
//! module (`counting/` for variant counting, `hla/` for HLA typing).
//!
//! ## Submodules
//!
//! - [`fragment`] — Fragment evidence tracking, QNAME hashing, insert size
//! - [`stats`] — Statistical tests (Fisher's exact, strand bias)
//! - [`bam_utils`] — BAM utility functions (quality helpers, position lookup)
//! - [`filters`] — Configurable BAM read filters (duplicate, secondary, MAPQ)
//! - [`baq`] — Base Alignment Quality heuristic (Li 2011)
//! - [`contig`] — Cross-source chromosome-name normalization

pub mod fragment;
pub mod stats;
pub mod bam_utils;
pub mod filters;
pub mod baq;
pub mod contig;

/// Resolve the worker-thread budget for one gbcms invocation.
///
/// `--threads` is the **total** thread budget for a single sample/process. Multi-sample
/// parallelism is provided by Nextflow, which runs gbcms as one of many concurrent
/// processes each pinned to `task.cpus` — so every parallel section here must stay
/// *within* this budget rather than grabbing all node cores. All rayon pools are sized
/// from this value; any future htslib decode threads must **subdivide** it, never add
/// to it.
///
/// Guards a rayon foot-gun: `ThreadPoolBuilder::num_threads(0)` is interpreted as "use
/// all logical cores", which would silently oversubscribe a small SLURM allocation. A
/// request of 0 is coerced to 1 with a warning (never silently).
pub fn resolve_thread_budget(requested: usize) -> usize {
    if requested == 0 {
        log::warn!("threads budget of 0 is invalid (rayon would grab all cores); using 1 instead");
        1
    } else {
        requested
    }
}

#[cfg(test)]
mod tests {
    use super::resolve_thread_budget;

    #[test]
    fn test_resolve_thread_budget_guards_zero() {
        assert_eq!(resolve_thread_budget(0), 1, "0 must clamp to 1, not rayon's all-cores");
        assert_eq!(resolve_thread_budget(1), 1);
        assert_eq!(resolve_thread_budget(8), 8);
    }
}
