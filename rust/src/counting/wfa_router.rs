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

use super::utils::{ClassifyResult, ClassifyPhase, MIN_USABLE_BASES};


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
/// ## Base-quality gate (CR-2)
///
/// A *definitive* REF/ALT call requires an exact (edit-distance 0) match on the
/// **BQ-masked** read: sub-`min_baseq` bases are masked to `N` first, so the
/// discriminating base — and every base in the window — must clear the same
/// base-quality gate the SW/PairHMM fallback applies. A low-BQ base becomes `N`,
/// breaks the exact match, and the read defers to the BQ-aware PairHMM instead of
/// getting a confident call the fallback would reject. This mirrors the low-BQ→`N`
/// masking in `alignment.rs`. The **off-target** decision stays on the *raw* read,
/// because sequence divergence is BQ-independent and masking must not inflate a
/// low-quality read into a false "off-target".
///
/// ## Parameters
///
/// - `read_seq`: raw read bases
/// - `read_quals`: per-base qualities, aligned 1:1 with `read_seq`
/// - `min_baseq`: bases below this are masked to `N` for the definitive-call triage
/// - `matrix`: haplotype matrix from `build_haplotype_matrix()`
/// - `med_qual`: pre-computed median quality for the ClassifyResult
pub fn wfa_fast_path(
    read_seq: &[u8],
    read_quals: &[u8],
    min_baseq: u8,
    matrix: &[Vec<u8>],
    med_qual: u8,
) -> Option<ClassifyResult> {
    if matrix.len() < 2 {
        trace!("wfa_fast_path: matrix has {} haplotypes (need ≥2) — skipping", matrix.len());
        return None;
    }

    // Share the usable-base gate with the SW/PairHMM backends (the "one quality
    // contract"): too few high-BQ bases means there is not enough signal for a
    // definitive call, so defer to the equally-gated fallback rather than risk a
    // call driven by one or two bases. In practice a definitive call already
    // requires an all-high-BQ exact match, so this only guards pathological
    // short/over-masked reads — but it keeps every backend's floor identical.
    let usable_count = read_quals.iter().filter(|&&q| q >= min_baseq).count();
    if usable_count < MIN_USABLE_BASES {
        trace!("wfa_fast_path: only {} usable bases — deferring to fallback", usable_count);
        return None;
    }

    // Mask sub-min_baseq bases to N for the definitive-call triage. A missing
    // qual (length mismatch) is treated as BQ 0 → masked, which is conservative.
    let masked: Vec<u8> = read_seq
        .iter()
        .enumerate()
        .map(|(i, &b)| {
            if read_quals.get(i).copied().unwrap_or(0) < min_baseq {
                b'N'
            } else {
                b
            }
        })
        .collect();

    // Common case (all bases ≥ min_baseq, e.g. duplex consensus): the masked read
    // equals the raw read, so the second alignment per haplotype is redundant.
    let needs_mask = masked.as_slice() != read_seq;

    let mut aligner = EditAligner::new(WavefrontPenalties::new_edit());

    // Best edit distance per class, on the RAW read (off-target decision) and the
    // BQ-MASKED read (definitive calls). Even indices = REF-class, odd = ALT-class.
    let (mut best_ref_raw, mut best_alt_raw) = (i32::MAX, i32::MAX);
    let (mut best_ref_masked, mut best_alt_masked) = (i32::MAX, i32::MAX);

    for (i, haplotype) in matrix.iter().enumerate() {
        let raw = aligner.align_end2end(read_seq, haplotype);
        let msk = if needs_mask {
            aligner.align_end2end(&masked, haplotype)
        } else {
            raw
        };

        trace!(
            "wfa_fast_path: hap[{}] ({}) raw={} masked={}",
            i, if i % 2 == 0 { "REF" } else { "ALT" }, raw, msk,
        );

        if i % 2 == 0 {
            best_ref_raw = best_ref_raw.min(raw);
            best_ref_masked = best_ref_masked.min(msk);
        } else {
            best_alt_raw = best_alt_raw.min(raw);
            best_alt_masked = best_alt_masked.min(msk);
        }
    }

    debug!(
        "wfa_fast_path: raw(ref={} alt={}) masked(ref={} alt={}) read_len={}",
        best_ref_raw, best_alt_raw, best_ref_masked, best_alt_masked, read_seq.len(),
    );

    // Definitive calls require an exact match on the BQ-MASKED read (CR-2): the
    // discriminating base must be high-BQ, else masking breaks the exact match.
    if best_ref_masked == 0 && best_alt_masked > 0 {
        Some(ClassifyResult::is_ref(med_qual, ClassifyPhase::Alignment))
    } else if best_alt_masked == 0 && best_ref_masked > 0 {
        Some(ClassifyResult::is_alt(med_qual, ClassifyPhase::Alignment))
    } else if best_ref_raw > OFF_TARGET_THRESHOLD && best_alt_raw > OFF_TARGET_THRESHOLD {
        // Off-target uses the RAW read — divergence is BQ-independent.
        debug!(
            "wfa_fast_path: off-target (raw ref={} alt={} > threshold={})",
            best_ref_raw, best_alt_raw, OFF_TARGET_THRESHOLD,
        );
        Some(ClassifyResult::neither(ClassifyPhase::Alignment))
    } else {
        // Ambiguous, or a low-BQ discriminating base → need BQ-aware PairHMM.
        None
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    const MINQ: u8 = 20;

    /// All-high-BQ quals of the given length (no masking).
    fn hq(n: usize) -> Vec<u8> {
        vec![30u8; n]
    }

    #[test]
    fn test_perfect_ref_match() {
        let read = b"AACGT";
        let matrix = vec![
            b"AACGT".to_vec(), // H0: REF (perfect match)
            b"AATGT".to_vec(), // H1: ALT (1 mismatch)
        ];
        let result = wfa_fast_path(read, &hq(read.len()), MINQ, &matrix, 30);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(r.is_ref);
        assert!(!r.is_alt);
    }

    #[test]
    fn test_perfect_alt_match() {
        let read = b"AATGT";
        let matrix = vec![
            b"AACGT".to_vec(), // H0: REF (1 mismatch)
            b"AATGT".to_vec(), // H1: ALT (perfect match)
        ];
        let result = wfa_fast_path(read, &hq(read.len()), MINQ, &matrix, 30);
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
            b"AACGT".to_vec(), // H0: REF (1 mismatch at pos 2)
            b"AATGT".to_vec(), // H1: ALT (1 mismatch at pos 2)
        ];
        let result = wfa_fast_path(read, &hq(read.len()), MINQ, &matrix, 30);
        assert!(result.is_none());
    }

    #[test]
    fn test_insufficient_usable_bases_defers() {
        // The read perfectly matches REF, but only 2 bases clear min_baseq, so the
        // shared usable-base gate defers to the fallback instead of calling REF.
        let read = b"AACGT";
        let matrix = vec![
            b"AACGT".to_vec(), // H0: REF — would be a perfect match
            b"AATGT".to_vec(), // H1: ALT
        ];
        let quals = vec![30u8, 30, 5, 5, 5]; // only 2 usable (< MIN_USABLE_BASES)
        assert!(
            wfa_fast_path(read, &quals, MINQ, &matrix, 30).is_none(),
            "fewer than MIN_USABLE_BASES usable bases must defer to the fallback",
        );
    }

    #[test]
    fn test_empty_matrix() {
        let result = wfa_fast_path(b"ACGT", &hq(4), MINQ, &[], 30);
        assert!(result.is_none());
    }

    #[test]
    fn test_with_siblings_ref_class_wins() {
        // H0: REF, H1: ALT, H2: REF+sib, H3: ALT+sib
        // Read matches H2 perfectly (REF + sibling)
        let read = b"AAYGT";
        let matrix = vec![
            b"AACGT".to_vec(), // H0: REF
            b"AATGT".to_vec(), // H1: ALT
            b"AAYGT".to_vec(), // H2: REF+sib (perfect match)
            b"AAXGT".to_vec(), // H3: ALT+sib
        ];
        let result = wfa_fast_path(read, &hq(read.len()), MINQ, &matrix, 30);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(r.is_ref);
    }

    // ── CR-2: base-quality gate on the discriminating base ──

    #[test]
    fn test_low_bq_discriminating_base_defers_to_pairhmm() {
        // Read matches ALT exactly, but the DISCRIMINATING base (pos 2) is Q2.
        // Without the gate this is a confident ALT; with it, the base masks to N,
        // breaks the exact match, and the read defers (None) to BQ-aware PairHMM.
        let read = b"AATGT";
        let quals = vec![30, 30, 2, 30, 30]; // pos 2 (the T vs C site) is low-BQ
        let matrix = vec![
            b"AACGT".to_vec(), // H0: REF
            b"AATGT".to_vec(), // H1: ALT (read matches this raw)
        ];
        let result = wfa_fast_path(read, &quals, MINQ, &matrix, 30);
        assert!(
            result.is_none(),
            "low-BQ discriminating base must defer to PairHMM, not call ALT",
        );
    }

    #[test]
    fn test_high_bq_discriminating_base_still_calls_alt() {
        // Same read/matrix as above but the discriminating base is high-BQ:
        // the fast path still makes the definitive ALT call.
        let read = b"AATGT";
        let quals = hq(read.len());
        let matrix = vec![b"AACGT".to_vec(), b"AATGT".to_vec()];
        let result = wfa_fast_path(read, &quals, MINQ, &matrix, 30);
        assert!(result.is_some());
        assert!(result.unwrap().is_alt);
    }
}
