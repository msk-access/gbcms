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

    // ── Mapping confidence ──────────────────────────────────────────────────────────
    /// Worst (minimum) MAPQ among the reads that contributed evidence to this fragment.
    ///
    /// Minimum rather than best, deliberately: a fragment is only as trustworthy as its
    /// least confidently placed read, and a consumer weighting evidence by error
    /// probability needs the pessimistic bound, not the flattering one.
    ///
    /// Initialized to `u8::MAX` (255), which the SAM spec defines as "mapping quality is
    /// unavailable" — so a fragment that somehow recorded no read reports *unknown* rather
    /// than the maximally-confident 0-error reading that a `0` initializer would imply.
    /// (0 is a real, meaningful MAPQ: multi-mapping. It must never arise by default.)
    pub min_mapq: u8,
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
            min_mapq: u8::MAX,
        }
    }

    /// Record a read's allele call, orientation, and fragment size into this
    /// fragment's evidence.
    ///
    /// Tracks per-allele orientation: when a new best-quality observation is
    /// recorded for REF or ALT, the orientation of THAT read is stored.
    /// This couples the strand direction to the winning evidence, not just R1.
    ///
    /// Tie-break (LO-6): the best-quality update uses a strict `>`, so on an
    /// exact base-quality tie the *first-observed* mate's orientation is kept.
    /// Reads arrive in coordinate order from a sorted BAM, so this is fully
    /// deterministic per input — but it follows iteration order, not a fixed
    /// R1/R2 rule. For a symmetric equal-quality overlap the FSB strand can
    /// therefore reflect whichever mate the fetch yielded first. This is a
    /// documented, minor source of FSB noise; an R1-preferring tie-break was
    /// considered and not adopted (it would perturb FSB p-values with no
    /// diagnostic benefit).
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
    /// - `mapq`: this read's MAPQ. Tracked as a running minimum in `min_mapq`.
    ///   Appended last on purpose: every other trailing parameter is a `bool`, so a
    ///   mis-ordered call fails to compile rather than silently swapping two `u8`s
    ///   (which is what inserting it next to `base_qual` would have risked).
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
        mapq: u8,
    ) {
        // Unconditional, and before any allele branching: a read that is neither REF nor
        // ALT still counts toward DPF, so its mapping confidence still describes the
        // fragment. Gating this on is_ref/is_alt would leave those fragments reporting 255.
        self.min_mapq = self.min_mapq.min(mapq);

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
///
/// Note (LO-4): the fragment map is keyed on this u64 with no stored QNAME, so a
/// birthday collision would silently merge two fragments. At realistic per-locus
/// depth the probability is ~5e-10 — not worth the memory of storing keys.
/// Revisit only if per-locus fragment counts grow orders of magnitude, then
/// store the QNAME for a collision tiebreak.
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

// ── Unit Tests ──────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Typical high-confidence MAPQ. These tests exercise allele resolution, not mapping
    /// confidence, so one uniform value keeps `min_mapq` out of the way; the tracking
    /// itself is covered by its own tests below.
    const TEST_MAPQ: u8 = 60;

    /// Helper: create a FragmentEvidence with specific REF/ALT observations.
    /// Simulates observe() calls with given qualities and structural flag.
    fn evidence_with(
        ref_qual: u8,
        alt_qual: u8,
        structural_alt: bool,
    ) -> FragmentEvidence {
        let mut ev = FragmentEvidence::new();
        if ref_qual > 0 {
            ev.observe(true, false, ref_qual, true, true, 200, false, false, TEST_MAPQ);
        }
        if alt_qual > 0 {
            ev.observe(false, true, alt_qual, false, false, 200, false, structural_alt, TEST_MAPQ);
        }
        ev
    }

    // ── min_mapq tracking ────────────────────────────────────────────

    #[test]
    fn min_mapq_keeps_the_worst_read_not_the_best() {
        // A fragment is only as trustworthy as its least confidently placed read. Taking
        // the max (or last) would let one well-placed mate launder a badly-placed one.
        let mut ev = FragmentEvidence::new();
        ev.observe(true, false, 30, true, true, 200, false, false, 60);
        ev.observe(false, true, 30, false, false, 200, false, false, 11);
        assert_eq!(ev.min_mapq, 11);
        // order must not matter
        let mut rev = FragmentEvidence::new();
        rev.observe(false, true, 30, false, false, 200, false, false, 11);
        rev.observe(true, false, 30, true, true, 200, false, false, 60);
        assert_eq!(rev.min_mapq, 11);
    }

    #[test]
    fn min_mapq_is_tracked_for_neither_ref_nor_alt_reads() {
        // A read that is neither REF nor ALT still counts toward DPF, so it still describes
        // the fragment's mapping confidence. Gating the update on is_ref/is_alt would leave
        // those fragments reporting the "unavailable" sentinel instead of what was measured.
        let mut ev = FragmentEvidence::new();
        ev.observe(false, false, 0, true, true, 200, false, false, 7);
        assert_eq!(ev.min_mapq, 7);
    }

    #[test]
    fn unobserved_fragment_reports_unavailable_not_zero() {
        // 0 is a REAL MAPQ meaning multi-mapping. If it were the initializer, a fragment
        // that recorded nothing would be indistinguishable from one placed ambiguously —
        // and a consumer weighting by error probability would silently discard good
        // evidence. 255 is the SAM spec's "unavailable".
        assert_eq!(FragmentEvidence::new().min_mapq, u8::MAX);
        assert_ne!(FragmentEvidence::new().min_mapq, 0);
    }

    #[test]
    fn mapq_zero_is_recorded_faithfully() {
        // The flip side: a genuine MAPQ 0 must survive as 0, not be treated as "missing".
        let mut ev = FragmentEvidence::new();
        ev.observe(true, false, 30, true, true, 200, false, false, 0);
        assert_eq!(ev.min_mapq, 0);
    }

    // ── resolve() structural priority tests ───────────────────────────

    #[test]
    fn resolve_structural_alt_wins_despite_equal_quality() {
        // Core fix: structural ALT wins even when REF and ALT have identical BQ.
        // This is the exact scenario that was causing fragment discard before
        // the INDEL consensus fix (both reads report anchor BQ = 79).
        let ev = evidence_with(79, 79, true);
        assert_eq!(ev.resolve(10), (false, true), "structural ALT should win");
    }

    #[test]
    fn resolve_structural_alt_wins_despite_higher_ref_quality() {
        // Structural ALT wins even when REF has significantly higher BQ.
        // The quality comparison is meaningless for INDEL evidence.
        let ev = evidence_with(90, 30, true);
        assert_eq!(ev.resolve(10), (false, true), "structural ALT should win");
    }

    #[test]
    fn resolve_non_structural_ref_wins_by_quality() {
        // Without structural flag, normal quality consensus applies.
        // REF quality exceeds ALT + threshold → REF wins.
        let ev = evidence_with(90, 30, false);
        assert_eq!(ev.resolve(10), (true, false), "REF should win by quality");
    }

    #[test]
    fn resolve_non_structural_alt_wins_by_quality() {
        // Without structural flag, ALT wins when ALT quality exceeds REF + threshold.
        let ev = evidence_with(30, 90, false);
        assert_eq!(ev.resolve(10), (false, true), "ALT should win by quality");
    }

    #[test]
    fn resolve_non_structural_tie_discards() {
        // Without structural flag, equal quality within threshold → discard.
        // This is the SNP conflict behavior (unchanged by this fix).
        let ev = evidence_with(30, 30, false);
        assert_eq!(ev.resolve(10), (false, false), "tie should discard");
    }

    #[test]
    fn resolve_ref_only_returns_ref() {
        // Only REF evidence seen → REF wins. No conflict to resolve.
        let ev = evidence_with(50, 0, false);
        assert_eq!(ev.resolve(10), (true, false), "REF-only should return REF");
    }

    #[test]
    fn resolve_alt_only_returns_alt() {
        // Only ALT evidence seen → ALT wins. No conflict to resolve.
        let ev = evidence_with(0, 50, false);
        assert_eq!(ev.resolve(10), (false, true), "ALT-only should return ALT");
    }

    #[test]
    fn resolve_no_evidence_returns_neither() {
        // No evidence at all → neither. Should not normally happen
        // (filtered upstream), but the function handles it gracefully.
        let ev = FragmentEvidence::new();
        assert_eq!(ev.resolve(10), (false, false), "no evidence → neither");
    }

    // ── observe() structural flag tests ──────────────────────────────

    #[test]
    fn observe_structural_flag_is_sticky() {
        // Once has_structural_alt is set, it should persist even if
        // a subsequent non-structural observation is made.
        let mut ev = FragmentEvidence::new();
        // First read: structural ALT
        ev.observe(false, true, 30, true, true, 200, false, true, TEST_MAPQ);
        assert!(ev.has_structural_alt, "should be set after structural ALT");
        // Second read: non-structural REF
        ev.observe(true, false, 90, false, false, 200, false, false, TEST_MAPQ);
        assert!(ev.has_structural_alt, "should remain set (sticky)");
    }

    #[test]
    fn observe_structural_ref_does_not_set_flag() {
        // is_structural=true on a REF observation should NOT set the flag.
        // This guards against a hypothetical bug where someone passes
        // is_structural=true with is_ref=true (shouldn't happen, but
        // defensive programming).
        let mut ev = FragmentEvidence::new();
        ev.observe(true, false, 50, true, true, 200, false, true, TEST_MAPQ);
        assert!(!ev.has_structural_alt, "REF obs should not set structural ALT flag");
    }

    #[test]
    fn observe_non_structural_alt_does_not_set_flag() {
        // Non-structural ALT (e.g., Phase 3 alignment) should not set the flag.
        let mut ev = FragmentEvidence::new();
        ev.observe(false, true, 50, true, true, 200, false, false, TEST_MAPQ);
        assert!(!ev.has_structural_alt, "non-structural ALT should not set flag");
    }
}
