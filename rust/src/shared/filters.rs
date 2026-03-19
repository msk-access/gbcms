//! Configurable BAM read filters shared across analysis modes.
//!
//! Provides a `ReadFilter` struct that encapsulates all universal BAM flag
//! checks (duplicates, secondary, supplementary, QC-failed, improper pair,
//! indel CIGAR). Mode-specific filtering (RNA NH rescue, MAPQ=0 tracking)
//! remains in the respective module's engine.
//!
//! Used by:
//! - `counting/engine.rs` — Phase 0 universal read filtering
//! - `hla/extract.rs` (future) — MHC region read extraction filtering

use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;
use log::trace;

/// Configurable BAM read filter for universal flag checks.
///
/// This struct handles filters that are the same across all analysis modes.
/// Mode-specific behavior (RNA NH:i tag rescue, MAPQ=0 tracking for MQ0
/// annotation) stays in the respective engine code.
pub struct ReadFilter {
    pub filter_duplicates: bool,
    pub filter_secondary: bool,
    pub filter_supplementary: bool,
    pub filter_qc_failed: bool,
    pub filter_improper_pair: bool,
    pub filter_indel: bool,
}

/// Tracks how many reads were rejected by each filter category.
///
/// Used for debug logging to diagnose unexpected filtering behavior.
#[derive(Debug, Default)]
pub struct FilterCounts {
    pub duplicates: u64,
    pub secondary: u64,
    pub supplementary: u64,
    pub qc_failed: u64,
    pub improper_pair: u64,
    pub indel: u64,
}

impl FilterCounts {
    /// Total number of reads rejected across all filter categories.
    pub fn total(&self) -> u64 {
        self.duplicates + self.secondary + self.supplementary
            + self.qc_failed + self.improper_pair + self.indel
    }
}

impl ReadFilter {
    /// Create a new `ReadFilter` with all filters enabled.
    pub fn all_enabled() -> Self {
        ReadFilter {
            filter_duplicates: true,
            filter_secondary: true,
            filter_supplementary: true,
            filter_qc_failed: true,
            filter_improper_pair: true,
            filter_indel: false, // Off by default — only for specific modes
        }
    }

    /// Check if a BAM record passes all enabled universal filters.
    ///
    /// Returns `true` if the record passes, `false` if it should be filtered out.
    /// When a record is filtered, the corresponding `FilterCounts` counter is
    /// incremented for debug logging.
    ///
    /// **Note:** MAPQ filtering is NOT handled here — it has mode-specific
    /// behavior (RNA NH rescue, MAPQ=0 tracking) and stays in the engine.
    pub fn passes(&self, record: &Record, counts: &mut FilterCounts) -> bool {
        if self.filter_duplicates && record.is_duplicate() {
            counts.duplicates += 1;
            trace!("ReadFilter: rejected duplicate");
            return false;
        }
        if self.filter_secondary && record.is_secondary() {
            counts.secondary += 1;
            trace!("ReadFilter: rejected secondary alignment");
            return false;
        }
        if self.filter_supplementary && record.is_supplementary() {
            counts.supplementary += 1;
            trace!("ReadFilter: rejected supplementary alignment");
            return false;
        }
        if self.filter_qc_failed && record.is_quality_check_failed() {
            counts.qc_failed += 1;
            trace!("ReadFilter: rejected QC-failed read");
            return false;
        }
        if self.filter_improper_pair && !record.is_proper_pair() {
            counts.improper_pair += 1;
            trace!("ReadFilter: rejected improper pair");
            return false;
        }
        if self.filter_indel {
            let has_indel = record.cigar().iter()
                .any(|op| matches!(op, Cigar::Ins(_) | Cigar::Del(_)));
            if has_indel {
                counts.indel += 1;
                trace!("ReadFilter: rejected read with CIGAR indel");
                return false;
            }
        }
        true
    }
}
