//! Fragment-level evidence tracking for read pair consensus.
//!
//! Provides `FragmentEvidence` for quality-weighted allele consensus across
//! R1/R2 reads of a fragment, and `hash_qname` / `hash_molecule` for
//! memory-efficient fragment tracking using u64 keys.
//!
//! Used by:
//! - `counting/engine.rs` — variant-level fragment dedup and consensus
//! - `hla/extract.rs` (future) — HLA read extraction dedup

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
    /// Physical fragment insert size (CIGAR-corrected) in base pairs.
    /// Updated via min() across both reads of the pair to keep the most
    /// corrected value (critical for deletions where only one read may span
    /// the indel in non-cfDNA contexts). TLEN=0 (unpaired/unmapped mate)
    /// is ignored. Only sizes in the cfDNA range (50–1000 bp) are later
    /// aggregated into mFSD class vectors.
    pub insert_size: Option<i32>,
    /// True if the base at the variant position was 'N' (ambiguous) on any read.
    /// Sticky: once set, never cleared. Used to split "neither-REF-nor-ALT"
    /// fragments into the N class vs. the NonREF class for mFSD analysis.
    pub has_n_base: bool,

    // ── INDEL structural evidence ───────────────────────────────────────────────────
    /// True if any read in this fragment had structural CIGAR evidence
    /// (I/D op) confirming ALT. Sticky across reads — once set, never
    /// cleared. Used by `resolve()` to prioritize structural evidence
    /// over base-quality comparisons in R1-vs-R2 INDEL conflicts.
    ///
    /// Set by `observe()` when `is_structural=true && is_alt=true`.
    /// Follows the same sticky-flag pattern as `has_n_base`.
    pub has_structural_alt: bool,
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
            has_structural_alt: false,
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
    /// - `tlen`: CIGAR-corrected physical insert size (`|TLEN| - D + I`).
    ///   Updated via `min()` across both reads to keep the most corrected value.
    ///   TLEN=0 (unpaired/unmapped mate) is skipped.
    /// - `is_n_base`: set `true` when the base at the variant position is 'N'.
    ///   Sticky across reads of the pair — once set, not cleared.
    /// ## Structural INDEL tracking
    /// - `is_structural`: set `true` when the read's classification came from
    ///   a direct CIGAR I/D op match (via `ClassifyResult::is_alt_structural`).
    ///   Sticky across reads of the pair — once set, not cleared. Enables
    ///   `resolve()` to prioritize structural evidence over BQ comparisons.
    #[allow(clippy::too_many_arguments)]
    pub fn observe(
        &mut self,
        is_ref: bool,
        is_alt: bool,
        base_qual: u8,
        is_read1: bool,
        is_forward: bool,
        tlen: i32,
        is_n_base: bool,
        is_structural: bool,
    ) {
        // Guard: a read cannot be both REF and ALT simultaneously.
        // ClassifyResult constructors enforce mutual exclusivity, so this
        // should never fire. Compiles out in --release builds; in production,
        // both best_*_qual update and resolve() falls through to quality
        // comparison — suboptimal but not catastrophic.
        debug_assert!(
            !(is_ref && is_alt),
            "FragmentEvidence::observe() called with both is_ref=true and is_alt=true"
        );

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

        // Structural CIGAR evidence: sticky flag (same pattern as has_n_base).
        // Once a read in this fragment shows a matching I/D op for ALT,
        // the flag stays set even if the other read has non-structural evidence.
        if is_structural && is_alt {
            self.has_structural_alt = true;
        }

        // mFSD: capture physical insert size — keep the MOST corrected value
        // across both reads. For deletions, the read carrying the D op gives a
        // smaller (more accurate) value; min() picks it. For insertions and
        // SNPs, both reads agree, so min() is a no-op.
        if tlen != 0 {
            let abs_physical = tlen.abs();
            match self.insert_size {
                None => self.insert_size = Some(abs_physical),
                Some(existing) => {
                    self.insert_size = Some(existing.min(abs_physical));
                }
            }
        }
        // mFSD: sticky N flag — once a read sees 'N' at this position, it stays
        if is_n_base {
            self.has_n_base = true;
        }
    }

    /// Resolve this fragment's allele call using quality-weighted consensus.
    ///
    /// Returns `(is_ref, is_alt)`:
    /// - `(true, false)` = REF wins
    /// - `(false, true)` = ALT wins
    /// - `(false, false)` = ambiguous, fragment discarded (counted in DPF
    ///   but not RDF/ADF)
    ///
    /// ## Structural priority (INDEL-aware)
    ///
    /// When a fragment has structural ALT evidence (CIGAR I/D op matching
    /// the target variant), ALT wins unconditionally regardless of quality.
    ///
    /// **Why**: For INDELs, both REF and ALT reads report anchor base quality
    /// (the base before the insertion/deletion). This quality measures "how
    /// confident is the anchor base call", NOT "how confident is the INDEL
    /// detection." Comparing anchor BQs to resolve an INDEL conflict is
    /// semantically meaningless — the CIGAR I/D op is the only signal that
    /// discriminates between the two alleles.
    ///
    /// This applies to both insertions and deletions:
    /// - INS: REF read's M-block is absence-of-insertion (zero-width in ref space)
    /// - DEL: REF read's M-block is the aligner's default when it didn't detect
    ///   the deletion (validated on DNMT3A duplex: all 7 conflict fragments were
    ///   genuine D-op evidence with MAPQ=60 and correct length)
    ///
    /// ## Known limitations
    ///
    /// Variants classified through Phase 3 alignment (complex variants,
    /// wrong-length INDELs) have `is_structural=false` and continue to use
    /// quality-weighted consensus. This is intentional: Phase 3 classifications
    /// are probabilistic, and quality arbitration is appropriate for them.
    pub fn resolve(&self, qual_diff_threshold: u8) -> (bool, bool) {
        let has_ref = self.best_ref_qual > 0;
        let has_alt = self.best_alt_qual > 0;

        match (has_ref, has_alt) {
            (true, false) => (true, false),   // Only REF evidence
            (false, true) => (false, true),   // Only ALT evidence
            (true, true) => {
                // Structural CIGAR evidence (I/D op) takes priority over
                // base-quality comparison. When a read has a matching I/D op,
                // the other read's "REF" is the aligner's default M-block —
                // absence of INDEL detection, not counter-evidence.
                if self.has_structural_alt {
                    return (false, true);
                }
                // Non-structural conflict (SNPs, Phase 3 returns):
                // quality-weighted consensus with threshold-based discard.
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
