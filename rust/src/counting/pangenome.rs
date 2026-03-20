//! Pangenomic haplotype matrix construction.
//!
//! Builds the multi-haplotype evaluation matrix for read classification.
//! Instead of a simple REF/ALT pair, the pangenomic approach constructs
//! all plausible local haplotypes by combining the test variant with
//! nearby germline (sibling) variants:
//!
//!   H0:  REF context (baseline)
//!   H1:  ALT injected into REF context
//!   H2…Hn:   REF + each sibling variant applied
//!   Hn+1…H2n: ALT + each sibling variant applied
//!
//! This prevents misclassification at multi-allelic or germline-adjacent
//! sites where a single REF/ALT comparison would confuse the somatic
//! variant with a nearby germline variant.
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


/// Build the pangenomic haplotype matrix for read evaluation.
///
/// Returns `Vec<Vec<u8>>` where:
///   - `[0]` = H0: pure REF haplotype (ref_context unchanged)
///   - `[1]` = H1: ALT allele spliced into REF context
///   - `[2..n+1]` = H2…Hn: each sibling's ALT spliced into REF context (without test ALT)
///   - `[n+1..2n+1]` = Hn+1…H2n: test ALT + each sibling's ALT both spliced in
///
/// Sibling variants are expected to be on the same chromosome within the
/// ref_context window. Siblings outside the context window are skipped
/// with a `trace!` log.
///
/// Returns `None` if the variant has no `ref_context` or the offset is invalid.
///
/// ## Performance
///
/// Matrix construction is O(k × L) where k = number of siblings + 1
/// and L = ref_context length. For typical variant calls (k < 10, L ≈ 40bp)
/// this is negligible vs the alignment cost.
pub fn build_haplotype_matrix(
    variant: &Variant,
    siblings: &[Variant],
) -> Option<Vec<Vec<u8>>> {
    let ref_context = variant.ref_context.as_ref()?.as_bytes();
    let ctx_start = variant.ref_context_start;
    let ctx_len = ref_context.len();

    // Offset of the test variant within the ref_context
    let offset = (variant.pos - ctx_start) as usize;
    let ref_len = variant.ref_allele.len();

    if offset + ref_len > ctx_len {
        trace!(
            "build_haplotype_matrix: offset {} + ref_len {} exceeds context len {}",
            offset, ref_len, ctx_len,
        );
        return None;
    }

    let left_ctx = &ref_context[..offset];
    let right_ctx = &ref_context[offset + ref_len..];

    // H0: pure REF haplotype (ref_context itself)
    let h0_ref: Vec<u8> = ref_context.to_vec();

    // H1: ALT spliced into REF context
    let h1_alt: Vec<u8> = left_ctx
        .iter()
        .chain(variant.alt_allele.as_bytes())
        .chain(right_ctx.iter())
        .copied()
        .collect();

    let mut matrix = vec![h0_ref, h1_alt];

    // For each sibling, create:
    //   H_ref_sib: sibling ALT spliced into REF context (test locus stays REF)
    //   H_alt_sib: both test ALT and sibling ALT spliced in
    for sib in siblings {
        // Check sibling is on the same chromosome
        if sib.chrom != variant.chrom {
            trace!(
                "build_haplotype_matrix: skipping sibling on different chrom: {} vs {}",
                sib.chrom, variant.chrom,
            );
            continue;
        }

        let sib_offset = (sib.pos - ctx_start) as usize;
        let sib_ref_len = sib.ref_allele.len();

        // Check sibling falls within the ref_context window
        if sib_offset + sib_ref_len > ctx_len {
            trace!(
                "build_haplotype_matrix: sibling at pos {} outside context window [{}, {})",
                sib.pos, ctx_start, ctx_start + ctx_len as i64,
            );
            continue;
        }

        // Skip siblings that overlap the test variant position
        // (would produce ambiguous haplotypes)
        if sib_offset < offset + ref_len && sib_offset + sib_ref_len > offset {
            trace!(
                "build_haplotype_matrix: sibling at pos {} overlaps test variant at pos {} — skipping",
                sib.pos, variant.pos,
            );
            continue;
        }

        // H_ref_sib: apply sibling ALT to REF context (test locus stays REF)
        let h_ref_sib = splice_allele(ref_context, sib_offset, sib_ref_len, sib.alt_allele.as_bytes());
        matrix.push(h_ref_sib);

        // H_alt_sib: apply BOTH test ALT and sibling ALT
        // Apply the one with the lower offset first to keep positions valid
        let h_alt_sib = if sib_offset < offset {
            // Sibling is upstream: apply sibling first, then test ALT
            let tmp = splice_allele(ref_context, sib_offset, sib_ref_len, sib.alt_allele.as_bytes());
            // Adjust test offset: if sibling ALT is shorter/longer, shift accordingly
            let delta = sib.alt_allele.len() as i64 - sib_ref_len as i64;
            let new_offset = (offset as i64 + delta) as usize;
            splice_allele(&tmp, new_offset, ref_len, variant.alt_allele.as_bytes())
        } else {
            // Sibling is downstream: apply test ALT first, then sibling
            let tmp = splice_allele(ref_context, offset, ref_len, variant.alt_allele.as_bytes());
            // Adjust sibling offset for any length change from test ALT
            let delta = variant.alt_allele.len() as i64 - ref_len as i64;
            let new_sib_offset = (sib_offset as i64 + delta) as usize;
            splice_allele(&tmp, new_sib_offset, sib_ref_len, sib.alt_allele.as_bytes())
        };
        matrix.push(h_alt_sib);
    }

    debug!(
        "build_haplotype_matrix: {} haplotypes for {}:{} ({} siblings)",
        matrix.len(), variant.chrom, variant.pos + 1, siblings.len(),
    );

    Some(matrix)
}


/// Splice an allele into a sequence at the given offset.
///
/// Replaces `seq[offset..offset+ref_len]` with `alt_allele`.
/// This is a purely functional helper — doesn't modify the input.
#[inline]
fn splice_allele(seq: &[u8], offset: usize, ref_len: usize, alt_allele: &[u8]) -> Vec<u8> {
    seq[..offset]
        .iter()
        .chain(alt_allele.iter())
        .chain(seq[offset + ref_len..].iter())
        .copied()
        .collect()
}


#[cfg(test)]
mod tests {
    use super::*;

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
        // Context: AAAAACGTTTT (pos 5 = C→T, sibling at pos 7 G→A)
        let v = make_variant(5, "C", "T", "AAAAACGTTTT", 0);
        let sib = make_variant(7, "T", "A", "AAAAACGTTTT", 0);
        let matrix = build_haplotype_matrix(&v, &[sib]).unwrap();
        assert_eq!(matrix.len(), 4);
        assert_eq!(matrix[0], b"AAAAACGTTTT"); // H0: REF
        assert_eq!(matrix[1], b"AAAAATGTTTT"); // H1: test ALT only
        assert_eq!(matrix[2], b"AAAAACGATTT"); // H2: sibling ALT only (REF at test locus)
        assert_eq!(matrix[3], b"AAAAATGATTT"); // H3: both test ALT + sibling ALT
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
        // Test variant at pos 5 ref_len=1, sibling at pos 5 (overlaps)
        let v = make_variant(5, "C", "T", "AAAAACGTTTT", 0);
        let sib = make_variant(5, "C", "A", "AAAAACGTTTT", 0);
        let matrix = build_haplotype_matrix(&v, &[sib]).unwrap();
        // Overlapping sibling should be skipped: only H0 + H1
        assert_eq!(matrix.len(), 2);
    }

    #[test]
    fn test_splice_allele_insertion() {
        let seq = b"AACGG";
        let result = splice_allele(seq, 2, 1, b"TTT"); // replace C with TTT
        assert_eq!(result, b"AATTTGG");
    }

    #[test]
    fn test_splice_allele_deletion() {
        let seq = b"AACCCGG";
        let result = splice_allele(seq, 2, 3, b"T"); // replace CCC with T
        assert_eq!(result, b"AATGG");
    }

    #[test]
    fn test_upstream_sibling() {
        // Context: XXXXXATTTTT (test at pos 5 A→G, sibling at pos 2 X→Y)
        // After sibling: XXYXXATTTTT → sibling replaces pos 2
        let v = make_variant(5, "A", "G", "XXXXXATTTTT", 0);
        let sib = make_variant(2, "X", "Y", "XXXXXATTTTT", 0);
        let matrix = build_haplotype_matrix(&v, &[sib]).unwrap();
        assert_eq!(matrix.len(), 4);
        assert_eq!(matrix[0], b"XXXXXATTTTT"); // H0: REF
        assert_eq!(matrix[1], b"XXXXXGTTTTT"); // H1: test ALT
        assert_eq!(matrix[2], b"XXYXXATTTTT"); // H2: sib ALT only
        assert_eq!(matrix[3], b"XXYXXGTTTTT"); // H3: both
    }
}
