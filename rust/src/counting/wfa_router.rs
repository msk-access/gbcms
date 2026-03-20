//! WFA (Wavefront Alignment) fast-path router for read triage.
//!
//! Uses edit-distance WFA to quickly evaluate a read against the full
//! haplotype matrix. Clear-cut reads (perfect match to one allele,
//! high distance to others) are classified directly without PairHMM.
//! Ambiguous reads fall through to the full PairHMM pipeline.
//!
//! ## Triage Logic
//!
//! 1. Compute edit distance of read vs each haplotype in the matrix
//! 2. Find best REF-class score (H0, H2, H4, ...) and best ALT-class score (H1, H3, H5, ...)
//! 3. If one is 0 and the other > 0 → definitive classification
//! 4. If both > OFF_TARGET_THRESHOLD → off-target read, classify as NEITHER
//! 5. Otherwise → ambiguous, return None for PairHMM fallback
//!
//! ## Performance
//!
//! WFA edit-distance is O(s²) where s is the alignment score (edit distance).
//! For reads with ≤2 mismatches vs the correct haplotype, this is extremely fast.
//! This should resolve 70-80% of reads without touching PairHMM.

use log::{debug, trace};
use wfa2lib_rs::aligner::EditAligner;
use wfa2lib_rs::penalties::WavefrontPenalties;

use super::utils::{ClassifyResult, ClassifyPhase};


/// Score threshold above which both alleles are considered off-target.
/// If the best edit distance to both REF-class and ALT-class haplotypes
/// exceeds this, the read is too divergent to classify reliably.
const OFF_TARGET_THRESHOLD: i32 = 20;


/// WFA fast-path triage against the pangenomic haplotype matrix.
///
/// Computes edit distance of `read_seq` against every haplotype in `matrix`.
/// Haplotypes at even indices (0, 2, 4, ...) are REF-class,
/// odd indices (1, 3, 5, ...) are ALT-class.
///
/// Returns `Some(ClassifyResult)` for definitive calls:
///   - Perfect match (score=0) to one class, >0 to other → REF or ALT
///   - Both >OFF_TARGET_THRESHOLD → NEITHER (off-target)
///
/// Returns `None` for ambiguous reads that need PairHMM:
///   - Tied scores between REF and ALT classes
///   - Small non-zero scores where BQ-aware classification matters
///
/// ## Parameters
///
/// - `read_seq`: raw read bases
/// - `matrix`: haplotype matrix from `build_haplotype_matrix()`
/// - `med_qual`: pre-computed median quality for the ClassifyResult
pub fn wfa_fast_path(
    read_seq: &[u8],
    matrix: &[Vec<u8>],
    med_qual: u8,
) -> Option<ClassifyResult> {
    if matrix.len() < 2 {
        trace!("wfa_fast_path: matrix has {} haplotypes (need ≥2) — skipping", matrix.len());
        return None;
    }

    let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());

    // Compute edit distances: even indices = REF-class, odd = ALT-class
    let mut best_ref_score = i32::MAX;
    let mut best_alt_score = i32::MAX;

    for (i, haplotype) in matrix.iter().enumerate() {
        let score = aligner.align_end2end(read_seq, haplotype);

        trace!(
            "wfa_fast_path: hap[{}] ({}) score={}",
            i, if i % 2 == 0 { "REF" } else { "ALT" }, score,
        );

        if i % 2 == 0 {
            // REF-class haplotype
            best_ref_score = best_ref_score.min(score);
        } else {
            // ALT-class haplotype
            best_alt_score = best_alt_score.min(score);
        }
    }

    debug!(
        "wfa_fast_path: best_ref={} best_alt={} read_len={}",
        best_ref_score, best_alt_score, read_seq.len(),
    );

    // Triage decision
    if best_ref_score == 0 && best_alt_score > 0 {
        // Perfect REF match, ALT has mismatches → definitive REF
        Some(ClassifyResult::is_ref(med_qual, ClassifyPhase::Alignment))
    } else if best_alt_score == 0 && best_ref_score > 0 {
        // Perfect ALT match, REF has mismatches → definitive ALT
        Some(ClassifyResult::is_alt(med_qual, ClassifyPhase::Alignment))
    } else if best_ref_score > OFF_TARGET_THRESHOLD && best_alt_score > OFF_TARGET_THRESHOLD {
        // Both scores very high → off-target read
        debug!(
            "wfa_fast_path: off-target (ref={} alt={} > threshold={})",
            best_ref_score, best_alt_score, OFF_TARGET_THRESHOLD,
        );
        Some(ClassifyResult::neither(ClassifyPhase::Alignment))
    } else {
        // Ambiguous: tied or close scores → need BQ-aware PairHMM
        None
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_perfect_ref_match() {
        let read = b"AACGT";
        let matrix = vec![
            b"AACGT".to_vec(),  // H0: REF (perfect match)
            b"AATGT".to_vec(),  // H1: ALT (1 mismatch)
        ];
        let result = wfa_fast_path(read, &matrix, 30);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(r.is_ref);
        assert!(!r.is_alt);
    }

    #[test]
    fn test_perfect_alt_match() {
        let read = b"AATGT";
        let matrix = vec![
            b"AACGT".to_vec(),  // H0: REF (1 mismatch)
            b"AATGT".to_vec(),  // H1: ALT (perfect match)
        ];
        let result = wfa_fast_path(read, &matrix, 30);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(!r.is_ref);
        assert!(r.is_alt);
    }

    #[test]
    fn test_ambiguous_returns_none() {
        // Read matches both equally (1 mismatch each)
        let read = b"AAXGT";
        let matrix = vec![
            b"AACGT".to_vec(),  // H0: REF (1 mismatch at pos 2)
            b"AATGT".to_vec(),  // H1: ALT (1 mismatch at pos 2)
        ];
        let result = wfa_fast_path(read, &matrix, 30);
        assert!(result.is_none());
    }

    #[test]
    fn test_empty_matrix() {
        let result = wfa_fast_path(b"ACGT", &[], 30);
        assert!(result.is_none());
    }

    #[test]
    fn test_with_siblings_ref_class_wins() {
        // H0: REF, H1: ALT, H2: REF+sib, H3: ALT+sib
        // Read matches H2 perfectly (REF + sibling)
        let read = b"AAYGT";
        let matrix = vec![
            b"AACGT".to_vec(),  // H0: REF
            b"AATGT".to_vec(),  // H1: ALT
            b"AAYGT".to_vec(),  // H2: REF+sib (perfect match)
            b"AAXGT".to_vec(),  // H3: ALT+sib
        ];
        let result = wfa_fast_path(read, &matrix, 30);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(r.is_ref);
    }
}
