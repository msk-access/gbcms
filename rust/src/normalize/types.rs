//! `PreparedVariant` — output type from variant preparation.

use pyo3::prelude::*;
use crate::types::Variant;

/// Result of variant preparation: normalized coords, ref_context, and validation info.
///
/// Every input variant produces exactly one `PreparedVariant`, even if validation
/// fails — this ensures the output always has the same row count as input.
#[pyclass]
#[derive(Debug, Clone)]
pub struct PreparedVariant {
    /// Ready-to-count variant (normalized coords + ref_context populated).
    /// For invalid variants, contains best-effort coords with no ref_context.
    #[pyo3(get)]
    pub variant: Variant,

    /// gbcms normalization status. Semicolon-separated multi-value.
    /// First token is always PASS or FAIL_*. Replaces old `validation_status`.
    /// Examples: "PASS", "PASS;WARN_REF_CORRECTED", "FAIL_REF_MISMATCH".
    #[pyo3(get, set)]
    pub gbcms_status: String,

    /// Post-counting diagnostic flags. Semicolon-separated.
    /// Empty string = no diagnostics. Set by Python pipeline after counting.
    /// Examples: "ZERO_ALT", "PARTIAL_DOMINANT;MNP_SPARSE_DISC(2/5)".
    #[pyo3(get, set)]
    pub gbcms_diagnostic: String,

    /// Rescue audit trail. Semicolon-separated key=value pairs.
    /// Contains original_alt=N when rescue was attempted.
    /// Empty string = no rescue. Only populated with --rescue-mnp (P1).
    #[pyo3(get, set)]
    pub gbcms_rescue: String,

    /// True if MAF anchor resolution (Step 1) changed pos/ref/alt.
    /// Only set for MAF input with dash alleles or different-length non-dash indels.
    #[pyo3(get)]
    pub was_anchor_resolved: bool,

    /// True if left-alignment (Step 3) shifted the variant's coordinates.
    #[pyo3(get)]
    pub was_left_aligned: bool,

    /// Original 0-based position before any transformation.
    #[pyo3(get)]
    pub original_pos: i64,

    /// Original REF allele before any transformation.
    #[pyo3(get)]
    pub original_ref: String,

    /// Original ALT allele before any transformation.
    #[pyo3(get)]
    pub original_alt: String,

    /// Corrected variant for homopolymer decomposition dual-counting.
    /// When a complex variant spans a homopolymer and appears to be a
    /// miscollapsed D(n)+SNV event, this holds the corrected allele
    /// (e.g., CCCCCC→CCCCT instead of CCCCCC→T).
    /// `None` for normal variants where no decomposition is detected.
    #[pyo3(get)]
    pub decomposed_variant: Option<Variant>,

    /// Group ID for overlapping multi-allelic variants at the same locus.
    /// `None` for isolated variants, `Some(id)` when multiple variants share
    /// overlapping genomic footprints (same chrom, overlapping REF spans).
    #[pyo3(get)]
    pub multi_allelic_group: Option<u32>,
}

#[pymethods]
impl PreparedVariant {
    /// Combined normalization flag: True if any transformation changed pos/ref/alt.
    ///
    /// Backward-compatible replacement for the old `was_normalized` field.
    /// Returns `was_anchor_resolved || was_left_aligned`.
    #[getter]
    pub fn was_normalized(&self) -> bool {
        self.was_anchor_resolved || self.was_left_aligned
    }
}
