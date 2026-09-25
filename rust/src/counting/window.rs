//! Discrimination windows: the reference stretch a read must span to tell a
//! variant's alleles apart.
//!
//! An indel in a repeat can sit anywhere along its shift-equivalence region
//! (the stretch it slides over without changing the haplotype). A read that
//! starts or ends inside that region matches both alleles: its bases read the
//! same whether the event is there or not, so the aligner places no gap and
//! the read says nothing about the event. GATK's AD likewise counts only
//! informative reads.
//!
//! Two uses:
//! - A REF call on a pure indel stands only when the read is informative:
//!   it covers a flank of the region and runs one base past the first base
//!   where REF and ALT differ reading inward from that flank
//!   ([`read_is_informative`]). Otherwise the read counts toward depth only.
//! - A co-annotated sibling's carriers are excluded from a row's REF only when
//!   the sibling's change lies inside that row's discrimination window (the
//!   region plus one base on each side, [`discrimination_window`]). A carrier
//!   whose change lies elsewhere shows the reference across every base that
//!   could distinguish this row's alleles, so it is REF here.

use rust_htslib::bam::record::{Cigar, Record};

use crate::normalize::repeat::first_change_offset;
use crate::types::Variant;

/// Reference interval `[lo, hi)` (0-based, half-open) holding every base at
/// which the variant's haplotype can differ from REF, over every equivalent
/// placement of the event.
///
/// - Pure deletion: the deleted bases slid both ways over their equivalence
///   region (a homopolymer deletion covers the whole tract).
/// - Pure insertion: the boundaries the insertion can slide over. It is empty
///   (`lo == hi`) in unique sequence, where the insert sits between `lo - 1`
///   and `lo`.
/// - Anything else (SNV, MNP, delins): the span after the shared leading bases.
///
/// Prep stores a pure indel's region on the variant (`shift_region`),
/// measured over a fetch sized to the event. Without it the slide runs over
/// `ref_context` and stops at its edges, so a region that runs past the
/// context is cut there: never wider than the true region, only narrower.
pub(crate) fn change_interval(v: &Variant) -> (i64, i64) {
    if let Some(region) = v.shift_region {
        return region;
    }
    let ctx = v.ref_context.as_deref().map_or(&[][..], str::as_bytes);
    shift_region_over(v.pos, &v.ref_allele, &v.alt_allele, |g| {
        let i = g - v.ref_context_start;
        (i >= 0 && (i as usize) < ctx.len()).then(|| ctx[i as usize].to_ascii_uppercase())
    })
}

/// Whether the alleles are a pure insertion or deletion: one allele is the
/// other plus bases after their shared prefix.
pub(crate) fn is_pure_indel(ref_allele: &str, alt_allele: &str) -> bool {
    let (r, a) = (ref_allele.len() as i64, alt_allele.len() as i64);
    let p = first_change_offset(ref_allele, alt_allele);
    (a < r && p == a) || (r < a && p == r)
}

/// [`change_interval`] for an event at `pos`, sliding over the reference
/// bases `base` returns (uppercase; None past its known extent).
pub(crate) fn shift_region_over(
    pos: i64,
    ref_allele: &str,
    alt_allele: &str,
    base: impl Fn(i64) -> Option<u8>,
) -> (i64, i64) {
    let r = ref_allele.as_bytes();
    let a = alt_allele.as_bytes();
    let p = first_change_offset(ref_allele, alt_allele) as usize;
    let start = pos + p as i64;

    if a.len() < r.len() && p == a.len() {
        // Pure deletion of [start, start + len).
        let len = (r.len() - p) as i64;
        let mut left = 0;
        while matches!((base(start - 1 - left), base(start + len - 1 - left)), (Some(x), Some(y)) if x == y) {
            left += 1;
        }
        let mut right = 0;
        while matches!((base(start + right), base(start + len + right)), (Some(x), Some(y)) if x == y) {
            right += 1;
        }
        return (start - left, start + len + right);
    }
    if r.len() < a.len() && p == r.len() {
        // Pure insertion of `ins` at the boundary before `start`.
        let ins: Vec<u8> = a[p..].iter().map(u8::to_ascii_uppercase).collect();
        let n = ins.len() as i64;
        let mut left = 0;
        while base(start - 1 - left) == Some(ins[(n - 1 - left % n) as usize]) {
            left += 1;
        }
        let mut right = 0;
        while base(start + right) == Some(ins[(right % n) as usize]) {
            right += 1;
        }
        return (start - left, start + right);
    }
    (start.min(pos + r.len() as i64), pos + r.len() as i64)
}

/// Reference interval `[lo, hi)` a read must span to discriminate the
/// variant's alleles: the change interval plus one base on each side for a
/// length-changing variant, the substituted bases themselves otherwise.
pub(crate) fn discrimination_window(v: &Variant) -> (i64, i64) {
    let (lo, hi) = change_interval(v);
    if v.ref_allele.len() == v.alt_allele.len() {
        (lo, hi)
    } else {
        (lo - 1, hi + 1)
    }
}

/// The two reference intervals `[lo, hi)`, one per side, of which a read must
/// span at least one to tell a length-changing variant's alleles apart.
///
/// Reading inward from a flank of the change interval, REF and ALT agree
/// until one base (the first discriminating base), because the repeat looks
/// the same with or without the event up to there. Each window runs from the
/// flank base to one base past that first discriminating base. The extra base
/// is the margin an aligner needs: an ALT read ending on the discriminating
/// base keeps one terminal mismatch (cheaper than a clip) and would look like
/// REF; one more base makes the mismatch pair cost more than the gap or clip.
///
/// - Deletion of `L` bases over `[lo, hi)`: first discriminating bases are
///   `hi - L` (from the left) and `lo + L - 1` (from the right). A homopolymer
///   deletion needs its whole tract and both flanks; a deletion longer than a
///   read needs only one junction.
/// - Insertion over boundaries `lo..=hi`: `hi` from the left, `lo - 1` from the
///   right, so both sides need the whole stretch.
///
/// `None` for anything but a pure indel: a substitution-bearing event has no
/// shift region, and its first base already discriminates.
pub(crate) fn informative_windows(v: &Variant) -> Option<[(i64, i64); 2]> {
    if !is_pure_indel(&v.ref_allele, &v.alt_allele) {
        return None;
    }
    let (lo, hi) = change_interval(v);
    let (r, a) = (v.ref_allele.len() as i64, v.alt_allele.len() as i64);
    if a < r {
        let len = r - a;
        Some([(lo - 1, hi - len + 2), (lo + len - 2, hi + 1)])
    } else {
        Some([(lo - 1, hi + 2), (lo - 2, hi + 1)])
    }
}

/// Whether the read is informative for a pure indel: its aligned reference
/// extent (M/=/X/D/N; clips excluded) spans either [`informative_windows`]
/// interval. Always true for other variants, which this rule does not cover.
pub(crate) fn read_is_informative(record: &Record, v: &Variant) -> bool {
    let Some(windows) = informative_windows(v) else {
        return true;
    };
    let start = record.pos();
    let mut end = start;
    for op in record.cigar().iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len)
            | Cigar::Del(len) | Cigar::RefSkip(len) => end += *len as i64,
            _ => {}
        }
    }
    windows.iter().any(|&(lo, hi)| start <= lo && end >= hi)
}

/// The siblings whose change lies inside `variant`'s discrimination window:
/// the only ones whose carriers are not REF testimony for `variant`.
pub(crate) fn siblings_in_window<'a>(variant: &Variant, siblings: &'a [Variant]) -> Vec<&'a Variant> {
    let (w_lo, w_hi) = discrimination_window(variant);
    siblings
        .iter()
        .filter(|s| {
            let (c_lo, c_hi) = change_interval(s);
            // Half-open overlap; an empty insertion interval [b, b) meets the
            // window when both of its flanking bases lie inside it.
            c_lo < w_hi && w_lo < c_hi
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A variant with a reference context starting at 0.
    fn var(ctx: &str, pos: i64, r: &str, a: &str) -> Variant {
        Variant::new(
            "1".into(), pos, r.into(), a.into(), "X".into(),
            Some(ctx.into()), 0, 0, None, None,
        )
    }

    //                0123456789012345
    const HOMO: &str = "TGCAAAAAAGTCCTGA"; // A-run at 3..9, G at 9

    #[test]
    fn homopolymer_deletion_covers_the_tract() {
        let v = var(HOMO, 2, "CA", "C");
        assert_eq!(change_interval(&v), (3, 9));
        assert_eq!(discrimination_window(&v), (2, 10));
    }

    #[test]
    fn homopolymer_deletion_not_left_aligned_still_covers_the_tract() {
        let v = var(HOMO, 5, "AA", "A");
        assert_eq!(change_interval(&v), (3, 9));
    }

    #[test]
    fn homopolymer_insertion_covers_the_tract() {
        let v = var(HOMO, 2, "C", "CA");
        assert_eq!(change_interval(&v), (3, 9));
        assert_eq!(discrimination_window(&v), (2, 10));
    }

    #[test]
    fn str_unit_deletion_slides_over_partial_copies() {
        // CACACAC at 3..10: three full CA copies and a trailing partial C.
        let v = var("TTGCACACACGTT", 2, "GCA", "G");
        assert_eq!(change_interval(&v), (3, 10));
    }

    #[test]
    fn single_base_deletion_inside_a_dinucleotide_repeat_stays_local() {
        // Deleting one A of ACACAC cannot slide: its neighbours are C.
        let v = var("TTGACACACGTT", 2, "GA", "G");
        assert_eq!(change_interval(&v), (3, 4));
    }

    #[test]
    fn unique_deletion_is_its_span() {
        let v = var("ACGTTGCAGTCA", 3, "TTG", "T");
        assert_eq!(change_interval(&v), (4, 6));
        assert_eq!(discrimination_window(&v), (3, 7));
    }

    #[test]
    fn unique_insertion_is_empty_and_its_window_is_the_two_flanks() {
        let v = var("ACGTAGCA", 3, "T", "TC");
        assert_eq!(change_interval(&v), (4, 4));
        assert_eq!(discrimination_window(&v), (3, 5));
    }

    #[test]
    fn snv_and_mnp_windows_are_their_bases() {
        assert_eq!(discrimination_window(&var(HOMO, 5, "A", "G")), (5, 6));
        assert_eq!(discrimination_window(&var(HOMO, 5, "AA", "GT")), (5, 7));
    }

    #[test]
    fn delins_is_its_trimmed_span_plus_one() {
        let v = var(HOMO, 9, "GTC", "GAAAT");
        assert_eq!(change_interval(&v), (10, 12));
        assert_eq!(discrimination_window(&v), (9, 13));
    }

    #[test]
    fn slide_stops_at_the_context_edge() {
        // Context ends inside the tract: the region is cut there.
        let v = var("TGCAAAAA", 2, "CA", "C");
        assert_eq!(change_interval(&v), (3, 8));
    }

    #[test]
    fn without_context_the_region_is_the_event_itself() {
        let v = Variant::new("1".into(), 2, "CA".into(), "C".into(), "DELETION".into(), None, 0, 0, None, None);
        assert_eq!(change_interval(&v), (3, 4));
    }

    #[test]
    fn homopolymer_deletion_needs_the_tract_and_both_flanks() {
        // A-run 3..9: X=C at 2, Y=G at 9.
        let v = var(HOMO, 2, "CA", "C");
        assert_eq!(informative_windows(&v), Some([(2, 10), (2, 10)]));
    }

    #[test]
    fn long_unique_deletion_needs_one_junction() {
        let mut ctx = String::from("ACGT");
        for i in 0..30 {
            ctx.push(b"ACGGTTCA"[i % 8] as char);
        }
        let v = var(&ctx, 3, &ctx[3..24], &ctx[3..4]); // 20bp deleted at 4..24
        let (lo, hi) = change_interval(&v);
        assert_eq!(informative_windows(&v), Some([(lo - 1, lo + 2), (hi - 2, hi + 1)]));
    }

    #[test]
    fn homopolymer_insertion_needs_the_tract_both_flanks_and_one_more() {
        let v = var(HOMO, 2, "C", "CA");
        assert_eq!(informative_windows(&v), Some([(2, 11), (1, 10)]));
    }

    #[test]
    fn substitution_bearing_events_have_no_informative_windows() {
        assert_eq!(informative_windows(&var(HOMO, 2, "C", "GA")), None);
        assert_eq!(informative_windows(&var(HOMO, 5, "A", "G")), None);
    }

    #[test]
    fn a_stored_shift_region_wins_over_the_context_slide() {
        // The context would cut this tract at 8; prep measured it to 12.
        let mut v = var("TGCAAAAA", 2, "CA", "C");
        assert_eq!(change_interval(&v), (3, 8));
        v.shift_region = Some((3, 12));
        assert_eq!(change_interval(&v), (3, 12));
        assert_eq!(informative_windows(&v), Some([(2, 13), (2, 13)]));
    }

    #[test]
    fn tandem_duplication_slides_over_the_duplicated_segment() {
        // Inserting a copy of ACGTTG right after itself: it can sit anywhere
        // along the segment, plus the next base that repeats the pattern.
        let ctx = "TTTACGTTGAGGG";
        let v = var(ctx, 2, "T", "TACGTTG");
        assert_eq!(change_interval(&v), (3, 10));
    }

    #[test]
    fn siblings_filter_by_window() {
        let ctx = "TGCAAAAAAGTCCTGACTTTTTTTGA"; // A-run 3..9, T-run 17..24
        let a = var(ctx, 2, "CA", "C"); // window 2..10
        let inside = var(ctx, 5, "A", "G");
        let twin = var(ctx, 2, "C", "CA");
        let outside = var(ctx, 16, "CT", "C");
        let left_of_window = var(ctx, 1, "G", "T");
        let sibs = vec![inside, twin, outside, left_of_window];
        let got: Vec<i64> = siblings_in_window(&a, &sibs).iter().map(|s| s.pos).collect();
        assert_eq!(got, vec![5, 2]);
    }
}
