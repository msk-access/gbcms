//! Fragment-level evidence tracking for read pair consensus.
//!
//! Provides `FragmentEvidence` for quality-weighted allele consensus across
//! R1/R2 reads of a fragment, and `hash_qname` for memory-efficient
//! fragment tracking using u64 keys.

use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;

/// Evidence accumulated for a single fragment (read pair) at a variant site.
/// Tracks the best base quality seen for each allele across both reads,
/// enabling quality-weighted consensus to resolve R1-vs-R2 conflicts.
///
/// Orientation tracking is per-allele: the strand direction stored for each
/// allele is the orientation of the read that provided the best-quality
/// evidence for that allele. This ensures Fisher's Strand Bias (FSB)
/// reflects the actual strand of the evidence, not just R1's strand.
///
/// ## mFSD fields
/// `insert_size` and `has_n_base` feed the mFSD engine during fragment
/// resolution, building the four class size vectors
/// (`ref_sizes`, `alt_sizes`, `nonref_sizes`, `n_sizes`).
#[derive(Debug, Clone)]
pub struct FragmentEvidence {
    /// Best base quality seen supporting REF across reads in this fragment
    pub best_ref_qual: u8,
    /// Best base quality seen supporting ALT across reads in this fragment
    pub best_alt_qual: u8,
    /// Orientation of the read providing best REF evidence
    best_ref_orientation: Option<bool>,
    /// Orientation of the read providing best ALT evidence
    best_alt_orientation: Option<bool>,
    /// Orientation of read 1 (fallback for orientation when allele-specific is unavailable)
    read1_orientation: Option<bool>,
    /// Orientation of read 2 (fallback)
    read2_orientation: Option<bool>,

    // ── mFSD Fragment Size Distribution fields ────────────────────────────
    /// Absolute insert size (|TLEN|) in base pairs, from the first non-zero TLEN
    /// seen across reads of this fragment. TLEN=0 (unpaired/unmapped mate) is
    /// ignored. Only sizes in the cfDNA range (50–1000 bp) are later aggregated.
    pub insert_size: Option<i32>,
    /// True if the base at the variant position was 'N' (ambiguous) on any read.
    /// Sticky: once set, never cleared. Used to split "neither-REF-nor-ALT"
    /// fragments into the N class vs. the NonREF class for mFSD analysis.
    pub has_n_base: bool,
}

impl FragmentEvidence {
    pub fn new() -> Self {
        FragmentEvidence {
            best_ref_qual: 0,
            best_alt_qual: 0,
            best_ref_orientation: None,
            best_alt_orientation: None,
            read1_orientation: None,
            read2_orientation: None,
            insert_size: None,
            has_n_base: false,
        }
    }

    /// Record a read's allele call, orientation, and fragment size into this
    /// fragment's evidence.
    ///
    /// Tracks per-allele orientation: when a new best-quality observation is
    /// recorded for REF or ALT, the orientation of THAT read is stored.
    /// This couples the strand direction to the winning evidence, not just R1.
    ///
    /// ## mFSD tracking
    /// - `tlen`: signed TLEN from the BAM record. Absolute value stored once
    ///   (first non-zero value wins). TLEN=0 (unpaired/unmapped mate) is skipped.
    /// - `is_n_base`: set `true` when the base at the variant position is 'N'.
    ///   Sticky across reads of the pair — once set, not cleared.
    #[allow(clippy::too_many_arguments)] // 5 existing + tlen + is_n_base: unavoidable
    pub fn observe(
        &mut self,
        is_ref: bool,
        is_alt: bool,
        base_qual: u8,
        is_read1: bool,
        is_forward: bool,
        tlen: i32,
        is_n_base: bool,
    ) {
        if is_ref && base_qual > self.best_ref_qual {
            self.best_ref_qual = base_qual;
            self.best_ref_orientation = Some(is_forward);
        }
        if is_alt && base_qual > self.best_alt_qual {
            self.best_alt_qual = base_qual;
            self.best_alt_orientation = Some(is_forward);
        }
        // Track R1/R2 orientation as fallback
        if is_read1 {
            self.read1_orientation = Some(is_forward);
        } else {
            self.read2_orientation = Some(is_forward);
        }

        // mFSD: capture insert size once (first non-zero value from either read)
        if self.insert_size.is_none() && tlen != 0 {
            self.insert_size = Some(tlen.abs());
        }
        // mFSD: sticky N flag — once a read sees 'N' at this position, it stays
        if is_n_base {
            self.has_n_base = true;
        }
    }

    /// Resolve this fragment's allele call using quality-weighted consensus.
    /// Returns (is_ref, is_alt) — exactly one will be true.
    pub fn resolve(&self, qual_diff_threshold: u8) -> (bool, bool) {
        let has_ref = self.best_ref_qual > 0;
        let has_alt = self.best_alt_qual > 0;

        match (has_ref, has_alt) {
            (true, false) => (true, false),   // Only REF evidence
            (false, true) => (false, true),   // Only ALT evidence
            (true, true) => {
                // Conflict: both alleles seen across reads in this fragment.
                // Use quality-weighted consensus: higher quality wins.
                // If quality difference is within threshold, discard the fragment
                // to avoid biasing VAF in either direction — critical for low-VAF
                // cfDNA detection where every fragment matters.
                if self.best_ref_qual > self.best_alt_qual + qual_diff_threshold {
                    (true, false)  // REF wins by quality margin
                } else if self.best_alt_qual > self.best_ref_qual + qual_diff_threshold {
                    (false, true)  // ALT wins by quality margin
                } else {
                    // Within threshold — ambiguous, discard to preserve VAF accuracy
                    (false, false)
                }
            }
            (false, false) => (false, false),  // Should not happen (filtered earlier)
        }
    }

    /// Get orientation for REF allele. Uses the strand of the read that provided
    /// the best REF evidence, falling back to R1 > R2 if not available.
    pub fn ref_orientation(&self) -> Option<bool> {
        self.best_ref_orientation
            .or(self.read1_orientation)
            .or(self.read2_orientation)
    }

    /// Get orientation for ALT allele. Uses the strand of the read that provided
    /// the best ALT evidence, falling back to R1 > R2 if not available.
    pub fn alt_orientation(&self) -> Option<bool> {
        self.best_alt_orientation
            .or(self.read1_orientation)
            .or(self.read2_orientation)
    }
}

/// Hash a QNAME to u64 for memory-efficient fragment tracking.
/// Using DefaultHasher for speed — collision probability is negligible
/// for typical variant-level read counts (~1000 fragments).
#[inline]
pub fn hash_qname(qname: &[u8]) -> u64 {
    let mut hasher = DefaultHasher::new();
    qname.hash(&mut hasher);
    hasher.finish()
}

/// Hash QNAME + optional UMI barcode for UMI-aware fragment grouping.
///
/// When `umi` is `Some(tag_bytes)`, the UMI is appended to the QNAME hash
/// with a separator byte (0xFF, which cannot appear in ASCII QNAME/UMI),
/// so reads with different UMIs are treated as distinct molecules even if
/// they share the same QNAME. This is critical for libraries where PCR
/// duplicates share QNAMEs but have distinct UMI barcodes.
///
/// When `umi` is `None`, this is identical to `hash_qname()`.
#[inline]
pub fn hash_molecule(qname: &[u8], umi: Option<&[u8]>) -> u64 {
    let mut hasher = DefaultHasher::new();
    qname.hash(&mut hasher);
    if let Some(umi_bytes) = umi {
        // Separator byte prevents QNAME="AB" + UMI="CD" from colliding
        // with QNAME="ABC" + UMI="D". 0xFF is invalid in ASCII BAM fields.
        0xFFu8.hash(&mut hasher);
        umi_bytes.hash(&mut hasher);
    }
    hasher.finish()
}
