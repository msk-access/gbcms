//! Variant normalization: left-alignment, MAF anchor resolution, REF validation.
//!
//! Consolidates all FASTA-dependent variant preparation into a single pass:
//! 1. MAF→VCF anchor base fetch (if `is_maf`)
//! 2. REF allele validation against reference
//! 3. bcftools-style left-alignment (`realign_left`)
//! 4. `ref_context` fetch for Smith-Waterman haplotype alignment
//!
//! Uses rayon `par_iter().map_init()` with thread-local FASTA readers,
//! matching the pattern in `counting::engine::count_bam()`.
//!
//! ## Submodules
//!
//! - [`types`] — `PreparedVariant` output struct
//! - [`decomp`] — Homopolymer decomposition detection
//! - [`left_align`] — bcftools `realign_left()` algorithm
//! - [`fasta`] — FASTA I/O, REF validation, MAF anchor resolution
//! - [`repeat`] — Tandem repeat detection and adaptive padding
//! - [`engine`] — `prepare_variants` orchestration + `prepare_single_variant`

mod types;
mod decomp;
mod left_align;
pub(crate) mod fasta;
// pub(crate): counting::variant_checks reuses find_tandem_repeat for the
// insertion truncation-containment rule (low-complexity gate).
pub(crate) mod repeat;
mod engine;

// Re-exports for lib.rs
pub use types::PreparedVariant;
pub use engine::prepare_variants;
