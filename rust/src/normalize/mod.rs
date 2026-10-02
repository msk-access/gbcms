//! Variant normalization: left-alignment, MAF anchor resolution, REF validation.
//!
//! Consolidates all FASTA-dependent variant preparation into a single pass:
//! 1. Malformed-allele rejection (EMPTY_ALLELE, ALT_EQUALS_REF)
//! 2. MAF→VCF anchor base fetch (if `is_maf`)
//! 3. REF allele validation against reference, and ALT N-base rejection
//! 4. bcftools-style left-alignment (`realign_left`)
//! 5. `ref_context` fetch for Phase-3 haplotype alignment (adaptively padded)
//! 6. Homopolymer twin (only with `rescue_homopolymer`)
//! 7. `repeat_span`, the shift region and the event reference (`event_ref`)
//!
//! then groups co-annotated variants (MULTI_ALLELIC / TRACT_CLUSTER).
//!
//! Uses rayon `par_iter().map_init()` with thread-local FASTA readers.
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
