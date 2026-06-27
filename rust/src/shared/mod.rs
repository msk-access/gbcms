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
