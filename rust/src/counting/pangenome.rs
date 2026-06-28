//! Pangenomic haplotype matrix construction.
//!
//! Builds the multi-haplotype evaluation matrix for read classification.
//! Instead of a simple REF/ALT pair, the pangenomic approach constructs
//! all plausible local haplotypes by combining the test variant with
//! nearby germline (sibling) variants using power-set combinatorics:
//!
//!   H0:       pure REF context (baseline)
//!   H1:       test ALT only
//!   H2..H2n:  every combination of sibling ALTs ± test ALT
//!
//! ## Right-to-Left Application Algorithm
//!
//! Mutations are applied from the rightmost genomic position to the
//! leftmost. This eliminates coordinate drift: each splice only shifts
//! bases to its right, preserving exact offsets for all remaining
//! variants to its left. This makes underflow/overflow mathematically
//! impossible by construction.
//!
//! ## Even/Odd Parity Contract
//!
//! Downstream consumers (`wfa_fast_path`, `classify_by_marginalized_pairhmm`)
//! classify haplotypes by matrix index:
//!   - Even indices (0, 2, 4, ...) = REF-class
//!   - Odd indices  (1, 3, 5, ...) = ALT-class
//!
//! This module guarantees that parity by always emitting (H_ref, H_alt)
//! pairs, with H0/H1 unconditionally placed at positions 0 and 1.
//!
//! ## Usage
//!
//! ```ignore
//! let matrix = build_haplotype_matrix(variant, siblings);
//! // matrix[0] = REF, matrix[1] = ALT, matrix[2..] = sibling combos
//! // Feed to wfa_fast_path() then classify_by_marginalized_pairhmm()
//! ```

use log::{debug, trace};
use crate::types::Variant;

/// Maximum haplotype length to prevent PairHMM O(N²) memory bloat
/// from stacked insertions. 400bp covers the largest adaptive context
/// window (ref_len + 2 × 50bp padding) with headroom.
const MAX_HAP_LEN: usize = 400;

/// Maximum siblings included in power-set combinatorics.
/// 2^4 = 16 combinations × 2 (ref+alt) = 32 haplotypes max.
/// Siblings are sorted by proximity; the closest 4 matter most.
const MAX_SIBS: usize = 4;


/// Splice an allele into a sequence, returning `None` if bounds are invalid.
///
/// Replaces `seq[offset..offset+ref_len]` with `alt_allele`.
/// Unlike the previous `splice_allele`, this function is bounds-checked:
/// it returns `None` instead of panicking if `offset` or `offset + ref_len`
/// exceeds `seq.len()`.
#[inline]
fn safe_splice_allele(
    seq: &[u8],
    offset: usize,
    ref_len: usize,
    alt_allele: &[u8],
) -> Option<Vec<u8>> {
    // Strict boundary check before any slicing occurs.
    // Catches: offset > seq.len(), integer overflow in offset + ref_len,
    // and ref_len extending past end of sequence.
    if offset > seq.len() || offset + ref_len > seq.len() {
        return None;
    }

    Some(
        seq[..offset]
            .iter()
            .chain(alt_allele.iter())
            .chain(seq[offset + ref_len..].iter())
            .copied()
            .collect(),
    )
}


/// Apply a set of variants to a reference context using right-to-left ordering.
///
/// Sorts variants by genomic position, validates no overlapping REF spans,
/// then applies from rightmost to leftmost. Right-to-left application
/// eliminates coordinate drift — each splice only shifts bases to its right,
/// preserving the exact offsets of all remaining variants to its left.
///
/// # Arguments
///
/// * `ref_ctx` — The reference context byte string to mutate
/// * `ctx_start` — Genomic position (0-based) of the first byte in `ref_ctx`
/// * `variants` — Tuples of (Variant, alt_allele_bytes) to apply.
///   The caller specifies allele bytes independently so the same Variant
///   can be applied with its REF or ALT allele without cloning/mutation.
/// * `max_hap_len` — Maximum allowed output length (memory bloat guard)
///
/// # Returns
///
/// * `Some(Vec<u8>)` — The mutated haplotype sequence
/// * `None` if:
///   - Any variant pair has overlapping REF spans (mutually exclusive on same strand)
///   - Any variant falls outside the context bounds (`pos < ctx_start`)
///   - The resulting haplotype exceeds `max_hap_len` (stacked insertion guard)
///   - A splice operation fails bounds checks
fn apply_variants_to_context(
    ref_ctx: &[u8],
    ctx_start: i64,
    variants: &mut [(&Variant, &[u8])],
    max_hap_len: usize,
) -> Option<Vec<u8>> {
    if variants.is_empty() {
        return Some(ref_ctx.to_vec());
    }

    // Sort left-to-right by genomic position for overlap detection
    variants.sort_by_key(|(v, _)| v.pos);

    // Overlap check: adjacent variants must not have overlapping REF spans.
    // Two variants physically overlap if left_end > right.pos (exclusive end).
    // Adjacent (left_end == right.pos) is allowed — different codons can co-exist.
    for i in 0..variants.len() - 1 {
        let (left, _) = &variants[i];
        let (right, _) = &variants[i + 1];
        let left_end = left.pos + left.ref_allele.len() as i64;
        if left_end > right.pos {
            trace!(
                "apply_variants_to_context: overlap — {}:{} (ref_end={}) overlaps {}:{}",
                left.chrom, left.pos + 1, left_end + 1,
                right.chrom, right.pos + 1,
            );
            return None;
        }
    }

    let mut current_seq = ref_ctx.to_vec();

    // Right-to-left application: process rightmost variant first.
    // Each splice shifts bases to the right, preserving left offsets.
    for (var, alt_bytes) in variants.iter().rev() {
        let offset_i64 = var.pos - ctx_start;

        // Guard: variant position before context start (signed check)
        if offset_i64 < 0 {
            trace!(
                "apply_variants_to_context: {}:{} before context start {} — aborting",
                var.chrom, var.pos + 1, ctx_start + 1,
            );
            return None;
        }

        current_seq = safe_splice_allele(
            &current_seq,
            offset_i64 as usize,
            var.ref_allele.len(),
            alt_bytes,
        )?;

        // Memory bloat guard: reject if stacked insertions blow past limit
        if current_seq.len() > max_hap_len {
            trace!(
                "apply_variants_to_context: haplotype length {} exceeds MAX_HAP_LEN {} — aborting",
                current_seq.len(), max_hap_len,
            );
            return None;
        }
    }

    Some(current_seq)
}


/// Build the pangenomic haplotype matrix for read evaluation.
///
/// Returns `Vec<Vec<u8>>` where:
///   - `[0]` = H0: pure REF haplotype (ref_context unchanged)
///   - `[1]` = H1: ALT allele spliced into REF context
///   - `[2..2n]` = power-set combinations of sibling ALTs with/without test ALT
///
/// Even indices are REF-class, odd indices are ALT-class. This parity
/// is enforced for compatibility with `wfa_fast_path` and
/// `classify_by_marginalized_pairhmm`.
///
/// Returns `None` if:
///   - The variant has no `ref_context`
///   - The test variant doesn't fit within its context window
///   - H1 construction fails (test ALT outside bounds)
///
/// ## Performance
///
/// Matrix construction is O(2^k × k × L) where k = min(siblings, 4)
/// and L = ref_context length. For typical calls (k ≤ 4, L ≈ 40bp),
/// this is ~5KB of byte ops — negligible vs PairHMM alignment cost.
pub fn build_haplotype_matrix(
    variant: &Variant,
    siblings: &[Variant],
) -> Option<Vec<Vec<u8>>> {
    let ref_context = variant.ref_context.as_ref()?.as_bytes();
    let ctx_start = variant.ref_context_start;

    // Validate test variant fits within its own context window.
    // Uses signed arithmetic to prevent underflow.
    let test_offset_i64 = variant.pos - ctx_start;
    if test_offset_i64 < 0
        || (test_offset_i64 as usize) + variant.ref_allele.len() > ref_context.len()
    {
        trace!(
            "build_haplotype_matrix: test variant {}:{} outside context [{}, +{})",
            variant.chrom, variant.pos + 1, ctx_start + 1, ref_context.len(),
        );
        return None;
    }

    // === Filter siblings: same chrom, within context bounds ===
    let mut valid_siblings: Vec<&Variant> = siblings
        .iter()
        .filter(|sib| {
            if sib.chrom != variant.chrom {
                trace!(
                    "build_haplotype_matrix: skipping sibling on different chrom: {} vs {}",
                    sib.chrom, variant.chrom,
                );
                return false;
            }
            let sib_offset = sib.pos - ctx_start;
            if sib_offset < 0 {
                trace!(
                    "build_haplotype_matrix: sibling at {}:{} before context start {} — skipping",
                    sib.chrom, sib.pos + 1, ctx_start + 1,
                );
                return false;
            }
            let sib_off_usize = sib_offset as usize;
            if sib_off_usize + sib.ref_allele.len() > ref_context.len() {
                trace!(
                    "build_haplotype_matrix: sibling at {}:{} extends past context end — skipping",
                    sib.chrom, sib.pos + 1,
                );
                return false;
            }
            true
        })
        .collect();

    // Sort by proximity to test variant (closest siblings matter most)
    valid_siblings.sort_by_key(|sib| (sib.pos - variant.pos).unsigned_abs());

    let max_sibs = valid_siblings.len().min(MAX_SIBS);

    // === H0 (pure REF) and H1 (test ALT only) — always at positions 0 and 1 ===
    let h0 = ref_context.to_vec();

    let mut test_only: Vec<(&Variant, &[u8])> = vec![
        (variant, variant.alt_allele.as_bytes()),
    ];
    let h1 = match apply_variants_to_context(ref_context, ctx_start, &mut test_only, MAX_HAP_LEN) {
        Some(h) => h,
        None => {
            trace!(
                "build_haplotype_matrix: H1 construction failed for {}:{}",
                variant.chrom, variant.pos + 1,
            );
            return None;
        }
    };

    let mut matrix: Vec<Vec<u8>> = vec![h0, h1];

    // === Power-set sibling combinatorics (masks 1..2^max_sibs) ===
    // mask=0 is H0/H1 (already built above).
    // For each mask, build:
    //   H_ref_sib (even): sibling combo only (test locus stays REF)
    //   H_alt_sib (odd):  sibling combo + test ALT
    for mask in 1u32..(1u32 << max_sibs) {
        // Build the sibling subset for this mask
        let sib_subset: Vec<&Variant> = (0..max_sibs)
            .filter(|i| (mask & (1 << i)) != 0)
            .map(|i| valid_siblings[i])
            .collect();

        // H_ref_sib: siblings only (test locus stays REF)
        let mut ref_vars: Vec<(&Variant, &[u8])> = sib_subset
            .iter()
            .map(|s| (*s, s.alt_allele.as_bytes()))
            .collect();

        let h_ref_sib = match apply_variants_to_context(
            ref_context, ctx_start, &mut ref_vars, MAX_HAP_LEN,
        ) {
            Some(h) => h,
            None => {
                // This sibling combo is invalid (overlap/bounds) — skip entirely.
                // Not pushing anything preserves matrix parity.
                trace!(
                    "build_haplotype_matrix: H_ref_sib failed for mask={:#06b} — skipping combo",
                    mask,
                );
                continue;
            }
        };

        // H_alt_sib: test ALT + siblings
        let mut alt_vars: Vec<(&Variant, &[u8])> = vec![
            (variant, variant.alt_allele.as_bytes()),
        ];
        alt_vars.extend(sib_subset.iter().map(|s| (*s, s.alt_allele.as_bytes())));

        let h_alt_sib = match apply_variants_to_context(
            ref_context, ctx_start, &mut alt_vars, MAX_HAP_LEN,
        ) {
            Some(h) => h,
            None => {
                // H_ref_sib succeeded but H_alt_sib failed (e.g., test ALT
                // overlaps a sibling). Don't push H_ref_sib alone — that
                // would break even/odd parity. Skip the entire combo.
                trace!(
                    "build_haplotype_matrix: H_alt_sib failed for mask={:#06b} — skipping combo",
                    mask,
                );
                continue;
            }
        };

        // Both succeeded — push the pair (REF-class even, ALT-class odd)
        if !matrix.contains(&h_ref_sib) || !matrix.contains(&h_alt_sib) {
            matrix.push(h_ref_sib);
            matrix.push(h_alt_sib);
        } else {
            // Both haplotypes are duplicates of existing entries.
            // Skip to avoid redundant PairHMM evaluations.
            trace!(
                "build_haplotype_matrix: duplicate haplotype pair for mask={:#06b} — skipping",
                mask,
            );
        }
    }

    // Minimum viable matrix: need at least H0 + H1
    if matrix.len() < 2 {
        return None;
    }

    debug!(
        "build_haplotype_matrix: {} haplotypes for {}:{} ({} valid siblings of {} total)",
        matrix.len(), variant.chrom, variant.pos + 1,
        valid_siblings.len().min(MAX_SIBS), siblings.len(),
    );

    Some(matrix)
}

/// Detect a non-discriminating locus: a REF-class haplotype (even index) that is
/// byte-identical to an ALT-class one (odd index). This arises when the test ALT,
/// combined with a sibling germline combination, reconstructs a reference sequence
/// (e.g. a homopolymer deletion cancelled by a nearby insertion of the same base).
/// When it occurs no read can tell REF from ALT — every read aligns equally to both
/// classes, the LLR is 0, and the locus would otherwise return all-NEITHER with no
/// signal. Callers surface this as a `NON_DISCRIMINATING_LOCUS` diagnostic.
pub fn has_ref_alt_collision(matrix: &[Vec<u8>]) -> bool {
    let ref_haps: std::collections::HashSet<&Vec<u8>> = matrix.iter().step_by(2).collect();
    matrix.iter().skip(1).step_by(2).any(|h| ref_haps.contains(h))
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_has_ref_alt_collision_homopolymer_cancel() {
        // 5-A homopolymer TG[AAAAA]GT (positions 102..=106), context starts at 100.
        // Test variant: delete one A (102 "AA">"A"). Sibling: insert one A (105 "A">"AA").
        // The deletion and insertion cancel, so test-ALT + sibling reconstructs REF:
        //   H0 (REF) = TGAAAAAGT, H3 (ALT = test del + sib ins) = TGAAAAAGT == H0.
        let ctx = "TGAAAAAGT";
        let test = make_variant(102, "AA", "A", ctx, 100);
        let sib = make_variant(105, "A", "AA", ctx, 100);
        let matrix = build_haplotype_matrix(&test, &[sib]).unwrap();
        assert!(
            has_ref_alt_collision(&matrix),
            "test-ALT + sibling reconstructs a REF-class haplotype — collision expected",
        );
    }

    #[test]
    fn test_no_collision_for_distinct_alleles() {
        // A normal multi-allelic locus where the two ALTs are distinct from REF and
        // from each other: no REF-class haplotype equals an ALT-class one.
        let ctx = "TGCATGCAT";
        let test = make_variant(102, "C", "G", ctx, 100); // SNP C>G
        let sib = make_variant(105, "G", "T", ctx, 100); // sibling SNP G>T
        let matrix = build_haplotype_matrix(&test, &[sib]).unwrap();
        assert!(
            !has_ref_alt_collision(&matrix),
            "distinct alleles must not collide",
        );
    }

    fn make_variant(pos: i64, ref_a: &str, alt_a: &str, ctx: &str, ctx_start: i64) -> Variant {
        Variant {
            chrom: "1".to_string(),
            pos,
            ref_allele: ref_a.to_string(),
            alt_allele: alt_a.to_string(),
            variant_type: String::new(),
            ref_context: Some(ctx.to_string()),
            ref_context_start: ctx_start,
            repeat_span: 0,
            gene_strand: None,
        }
    }

    // ── Existing tests (updated for new API) ──

    #[test]
    fn test_no_siblings_produces_two_haplotypes() {
        // Context: AAAAACGGGGG (pos 5 = C→T)
        let v = make_variant(5, "C", "T", "AAAAACGGGGG", 0);
        let matrix = build_haplotype_matrix(&v, &[]).unwrap();
        assert_eq!(matrix.len(), 2);
        assert_eq!(matrix[0], b"AAAAACGGGGG"); // H0: REF
        assert_eq!(matrix[1], b"AAAAATGGGGG"); // H1: ALT
    }

    #[test]
    fn test_one_sibling_produces_four_haplotypes() {
        // Context: AAAAACGTTTT (pos 5 = C→T, sibling at pos 7 T→A)
        let v = make_variant(5, "C", "T", "AAAAACGTTTT", 0);
        let sib = make_variant(7, "T", "A", "AAAAACGTTTT", 0);
        let matrix = build_haplotype_matrix(&v, &[sib]).unwrap();
        assert_eq!(matrix.len(), 4);
        assert_eq!(matrix[0], b"AAAAACGTTTT"); // H0: REF
        assert_eq!(matrix[1], b"AAAAATGTTTT"); // H1: test ALT only
        assert_eq!(matrix[2], b"AAAAACGATTT"); // H2: sibling ALT only (REF-class)
        assert_eq!(matrix[3], b"AAAAATGATTT"); // H3: both (ALT-class)
    }

    #[test]
    fn test_no_ref_context_returns_none() {
        let v = Variant {
            chrom: "1".to_string(),
            pos: 5,
            ref_allele: "C".to_string(),
            alt_allele: "T".to_string(),
            variant_type: String::new(),
            ref_context: None,
            ref_context_start: 0,
            repeat_span: 0,
            gene_strand: None,
        };
        assert!(build_haplotype_matrix(&v, &[]).is_none());
    }

    #[test]
    fn test_overlapping_sibling_is_skipped() {
        // Test variant at pos 5 ref_len=1, sibling at pos 5 (overlaps test)
        // Overlap: test REF span [5,6) overlaps sib REF span [5,6)
        // The sibling combo will fail in apply_variants_to_context overlap check
        let v = make_variant(5, "C", "T", "AAAAACGTTTT", 0);
        let sib = make_variant(5, "C", "A", "AAAAACGTTTT", 0);
        let matrix = build_haplotype_matrix(&v, &[sib]).unwrap();
        // Overlapping sibling combo skipped: only H0 + H1
        assert_eq!(matrix.len(), 2);
    }

    #[test]
    fn test_safe_splice_allele_insertion() {
        let seq = b"AACGG";
        let result = safe_splice_allele(seq, 2, 1, b"TTT").unwrap(); // replace C with TTT
        assert_eq!(result, b"AATTTGG");
    }

    #[test]
    fn test_safe_splice_allele_deletion() {
        let seq = b"AACCCGG";
        let result = safe_splice_allele(seq, 2, 3, b"T").unwrap(); // replace CCC with T
        assert_eq!(result, b"AATGG");
    }

    #[test]
    fn test_upstream_sibling() {
        // Context: XXXXXATTTTT (test at pos 5 A→G, sibling at pos 2 X→Y)
        let v = make_variant(5, "A", "G", "XXXXXATTTTT", 0);
        let sib = make_variant(2, "X", "Y", "XXXXXATTTTT", 0);
        let matrix = build_haplotype_matrix(&v, &[sib]).unwrap();
        assert_eq!(matrix.len(), 4);
        assert_eq!(matrix[0], b"XXXXXATTTTT"); // H0: REF
        assert_eq!(matrix[1], b"XXXXXGTTTTT"); // H1: test ALT
        assert_eq!(matrix[2], b"XXYXXATTTTT"); // H2: sib ALT only (REF-class)
        assert_eq!(matrix[3], b"XXYXXGTTTTT"); // H3: both (ALT-class)
    }

    // ── New regression tests for underflow prevention ──

    #[test]
    fn test_upstream_sib_large_deletion_no_panic() {
        // Upstream sibling with large deletion (sib_ref_len >> sib_alt_len).
        // Old code: delta = 1 - 10 = -9, new_offset = (5 + (-9)) as usize = usize::MAX → PANIC
        // New code: overlap detected (sib REF span [2,12) overlaps test at [5,6)) → skipped
        let v = make_variant(5, "C", "T", "AAAAACGGGGGGGGGG", 0); // 16bp context
        let sib = make_variant(2, "AAACGGGGGG", "X", "AAAAACGGGGGGGGGG", 0); // 10bp del at pos 2
        let matrix = build_haplotype_matrix(&v, &[sib]).unwrap();
        // Sibling's REF span [2, 12) overlaps test variant at pos 5 → combo skipped
        assert_eq!(matrix.len(), 2); // Only H0 + H1
    }

    #[test]
    fn test_downstream_sib_large_test_deletion_no_panic() {
        // Test variant is a large deletion, downstream sibling is close.
        // Old code: delta = 2 - 20 = -18, new_sib_offset = (25 + (-18)) as usize → could underflow
        // if sib_offset < |delta|.
        // New code: overlap detected (test REF span [5,25) overlaps sib at [10,...)) → skipped
        let ctx = "AAAAACCCCCCCCCCCCCCCCCCCCCTTTTTTTTTT"; // 35bp context
        let v = make_variant(5, "CCCCCCCCCCCCCCCCCCCC", "CT", ctx, 0); // 20bp del
        let sib = make_variant(10, "C", "G", ctx, 0); // SNP at pos 10 (inside del region)
        let matrix = build_haplotype_matrix(&v, &[sib]).unwrap();
        // Sib at pos 10 is inside test's REF span [5, 25) → overlap → skipped
        assert_eq!(matrix.len(), 2);
    }

    #[test]
    fn test_sib_before_context_start_filtered() {
        // Sibling position is before context_start.
        // Old code: sib.pos - ctx_start = 3 - 5 = -2 → (as usize) = usize::MAX → PANIC
        // New code: signed check filters it out safely.
        let v = make_variant(10, "C", "T", "AAAAACGGGGG", 5); // ctx_start=5
        let sib = make_variant(3, "A", "G", "AAAAACGGGGG", 5); // sib at pos 3, before ctx_start=5
        let matrix = build_haplotype_matrix(&v, &[sib]).unwrap();
        // Sibling filtered out (pos 3 < ctx_start 5) → only H0 + H1
        assert_eq!(matrix.len(), 2);
    }

    #[test]
    fn test_memory_bloat_guard() {
        // 4 siblings each inserting 100bp. Total haplotype would be
        // 20 + 4*100 = 420bp, exceeding MAX_HAP_LEN=400.
        let ctx = "AAAAACCCCCTTTTTTTTTT"; // 20bp context
        let v = make_variant(5, "C", "T", ctx, 0);
        let big_insert = "X".repeat(100);
        let sib1 = make_variant(7, "C", &big_insert, ctx, 0);
        let sib2 = make_variant(10, "T", &big_insert, ctx, 0);
        let sib3 = make_variant(13, "T", &big_insert, ctx, 0);
        let sib4 = make_variant(16, "T", &big_insert, ctx, 0);
        let matrix = build_haplotype_matrix(&v, &[sib1, sib2, sib3, sib4]).unwrap();
        // H0 + H1 always present. Combos with many insertions should be rejected
        // by the MAX_HAP_LEN guard, but individual sibling combos (1 sib = 120bp)
        // are under the limit.
        assert!(matrix.len() >= 2);
        // No haplotype should exceed MAX_HAP_LEN
        for h in &matrix {
            assert!(h.len() <= MAX_HAP_LEN, "Haplotype length {} > {}", h.len(), MAX_HAP_LEN);
        }
    }

    // ── safe_splice_allele boundary tests ──

    #[test]
    fn test_safe_splice_allele_out_of_bounds() {
        let seq = b"ACGT";
        // offset beyond seq length
        assert!(safe_splice_allele(seq, 10, 1, b"X").is_none());
        // offset + ref_len beyond seq length
        assert!(safe_splice_allele(seq, 3, 5, b"X").is_none());
        // exact boundary: offset=4, ref_len=0 → valid (empty splice at end)
        assert!(safe_splice_allele(seq, 4, 0, b"X").is_some());
    }

    // ── apply_variants_to_context tests ──

    #[test]
    fn test_apply_variants_empty() {
        let ctx = b"ACGTACGT";
        let result = apply_variants_to_context(ctx, 0, &mut [], 400);
        assert_eq!(result.unwrap(), ctx.to_vec());
    }

    #[test]
    fn test_apply_variants_right_to_left_correctness() {
        // Two non-overlapping SNPs: pos 2 (C→X) and pos 5 (C→Y)
        // Right-to-left: apply pos 5 first, then pos 2 — both offsets correct.
        let ctx = b"AACGTACGT"; // 9bp, ctx_start=0
        let v1 = make_variant(2, "C", "X", "AACGTACGT", 0);
        let v2 = make_variant(5, "A", "Y", "AACGTACGT", 0);
        let mut vars: Vec<(&Variant, &[u8])> = vec![
            (&v2, b"Y"),
            (&v1, b"X"),
        ];
        let result = apply_variants_to_context(ctx, 0, &mut vars, 400).unwrap();
        assert_eq!(result, b"AAXGTYCGT");
    }

    #[test]
    fn test_apply_variants_overlap_rejected() {
        // Two variants whose REF spans overlap: pos 2 ref_len=3 and pos 4 ref_len=1
        // left_end = 2 + 3 = 5 > 4 = right.pos → overlap
        let ctx = b"AACGTACGT";
        let v1 = make_variant(2, "CGT", "X", "AACGTACGT", 0);
        let v2 = make_variant(4, "A", "Y", "AACGTACGT", 0);
        let mut vars: Vec<(&Variant, &[u8])> = vec![
            (&v1, b"X"),
            (&v2, b"Y"),
        ];
        assert!(apply_variants_to_context(ctx, 0, &mut vars, 400).is_none());
    }
}
