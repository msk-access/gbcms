//! Heuristic Base Alignment Quality (BAQ) adjustment (Li 2011).
//!
//! Downgrades base qualities near CIGAR insertions and deletions to
//! reduce false positive variant calls caused by alignment artifacts.
//! Walks the CIGAR string and subtracts `BAQ_PENALTY` from base
//! qualities within `BAQ_RADIUS` bp of any Ins/Del operation.
//!
//! A splice junction (CIGAR N, `RefSkip`) lowers no base's quality. A splice
//! moves no base, so it says nothing about the quality of the bases beside it,
//! and a base is judged by its own quality. GATK's SplitNCigarReads changes no
//! quality at a splice, and htslib's BAQ skips a read with a reference skip whole.
//!
//! Used by:
//! - `counting/engine.rs` — Phase 0 quality adjustment before classification

use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;

/// Number of base-pairs on either side of an insertion or deletion to penalize.
const BAQ_RADIUS: usize = 5;
/// Quality penalty subtracted from BQ near an insertion or deletion (clamped to 0).
const BAQ_PENALTY: u8 = 20;

/// Apply heuristic BAQ quality downgrade to a read's base qualities.
///
/// Walks the CIGAR string, finds Ins/Del operations, and subtracts
/// `BAQ_PENALTY` (20) from base qualities within `BAQ_RADIUS` (5bp) of
/// the boundary on the read. Qualities are clamped to 0 (never negative).
/// A splice junction (RefSkip) lowers nothing.
///
/// Returns `None` if the read has no insertion or deletion (caller
/// should use original quals), a spliced read without one included.
/// Returns `Some(adjusted_quals)` with the modified quality vector otherwise.
///
/// This function is intentionally allocation-free for reads without
/// indels (the common case), only allocating the Vec when adjustment is needed.
///
/// `spare` is the genomic span `[lo, hi]` of the variant being counted: a read's
/// insertion or deletion touching it IS that variant's evidence, not an
/// alignment artifact beside it, so its bases are not penalized (at Q37 the
/// penalty would put every ALT read of a small indel below min BQ, while REF
/// reads, with no indel, keep full quality). Indels elsewhere in the read are
/// penalized as before.
pub fn apply_heuristic_baq(record: &Record, spare: Option<(i64, i64)>) -> Option<Vec<u8>> {
    let cigar = record.cigar();

    // Quick scan: does this read have any insertion or deletion? A splice
    // junction (RefSkip, CIGAR N) is not one: it lowers no base's quality.
    let needs_adjustment = cigar.iter().any(|op| matches!(op, Cigar::Ins(_) | Cigar::Del(_)));
    if !needs_adjustment {
        return None;
    }

    // Trace-level logging for per-read BAQ diagnostics (--trace flag).
    if log::log_enabled!(log::Level::Trace) {
        let n_indel = cigar
            .iter()
            .filter(|op| matches!(op, Cigar::Ins(_) | Cigar::Del(_)))
            .count();
        log::trace!(
            "BAQ: adjusting read {} ({} indels)",
            std::str::from_utf8(record.qname()).unwrap_or("?"),
            n_indel,
        );
    }

    let quals = record.qual().to_vec();
    let mut adjusted = quals;
    let read_len = adjusted.len();

    // Walk CIGAR to find read-coordinate positions of indel boundaries.
    // `read_pos` tracks the current position on the read (query) sequence and
    // `ref_pos` on the reference, to tell the variant's own indel from others.
    let mut read_pos: usize = 0;
    let mut ref_pos = record.pos();
    let mut changed = false;
    let own = |start: i64, end: i64| spare.is_some_and(|(lo, hi)| start <= hi && end >= lo);
    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                read_pos += *len as usize;
                ref_pos += *len as i64;
            }
            // The variant's own insertion: its bases are the evidence.
            Cigar::Ins(len) if own(ref_pos, ref_pos) => read_pos += *len as usize,
            // The variant's own deletion.
            Cigar::Del(len) if own(ref_pos, ref_pos + *len as i64) => ref_pos += *len as i64,
            Cigar::Ins(len) => {
                // Insertion: penalize BAQ_RADIUS bases before the insertion start
                // and after the insertion end on the read.
                let ins_start = read_pos;
                let ins_end = read_pos + *len as usize;

                let pen_start = ins_start.saturating_sub(BAQ_RADIUS);
                let pen_end = (ins_end + BAQ_RADIUS).min(read_len);

                for q in adjusted[pen_start..pen_end].iter_mut() {
                    *q = q.saturating_sub(BAQ_PENALTY);
                }
                changed = true;
                read_pos = ins_end;
            }
            // A splice junction consumes reference only and lowers nothing.
            Cigar::RefSkip(len) => ref_pos += *len as i64,
            Cigar::Del(len) => {
                // Deletion: penalize BAQ_RADIUS bases on either side of the
                // boundary on the read. It consumes reference but not query,
                // so read_pos stays.
                let pen_start = read_pos.saturating_sub(BAQ_RADIUS);
                let pen_end = (read_pos + BAQ_RADIUS).min(read_len);

                for q in adjusted[pen_start..pen_end].iter_mut() {
                    *q = q.saturating_sub(BAQ_PENALTY);
                }
                changed = true;
                ref_pos += *len as i64;
            }
            Cigar::SoftClip(len) => {
                read_pos += *len as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {
                // No read or reference consumption
            }
        }
    }

    changed.then_some(adjusted)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;

    fn record(seq: &[u8], pos: i64, cigar: Vec<Cigar>, q: u8) -> Record {
        let mut rec = Record::new();
        rec.set(b"r", Some(&CigarString(cigar)), seq, &vec![q; seq.len()]);
        rec.set_pos(pos);
        rec
    }

    /// A spliced read, 20 bases each side of a 100-base intron at 1000.
    fn spliced() -> Record {
        record(&[b'A'; 40], 980, vec![Cigar::Match(20), Cigar::RefSkip(100), Cigar::Match(20)], 37)
    }

    #[test]
    fn a_splice_lowers_no_base_quality() {
        assert_eq!(apply_heuristic_baq(&spliced(), None), None, "a spliced read without an indel is left alone");
    }

    #[test]
    fn an_indel_keeps_its_penalty_on_a_spliced_read() {
        // a 1-base deletion 10 bases into the second exon, after the junction
        let rec = record(
            &[b'A'; 40],
            980,
            vec![Cigar::Match(20), Cigar::RefSkip(100), Cigar::Match(10), Cigar::Del(1), Cigar::Match(10)],
            37,
        );
        let q = apply_heuristic_baq(&rec, None).unwrap();
        assert_eq!(&q[..25], &[37; 25], "the junction's bases and the bases before the deletion's radius keep their quality");
        assert_eq!(&q[25..35], &[17; 10], "the deletion's five bases each side are lowered");
        assert_eq!(&q[35..], &[37; 5]);
    }

    #[test]
    fn a_splice_beside_a_deletion_does_not_move_the_deletions_bases() {
        // the splice comes first: the deletion's read position must not shift
        let with_splice = record(
            &[b'A'; 40],
            980,
            vec![Cigar::Match(20), Cigar::RefSkip(100), Cigar::Match(10), Cigar::Del(1), Cigar::Match(10)],
            37,
        );
        let without = record(&[b'A'; 40], 980, vec![Cigar::Match(30), Cigar::Del(1), Cigar::Match(10)], 37);
        assert_eq!(apply_heuristic_baq(&with_splice, None), apply_heuristic_baq(&without, None));
    }

    #[test]
    fn a_variants_own_deletion_is_spared_on_a_spliced_read() {
        // the deletion sits at reference 1112..1113 on the spliced read (after the 100-base intron)
        let rec = record(
            &[b'A'; 40],
            980,
            vec![Cigar::Match(20), Cigar::RefSkip(100), Cigar::Match(10), Cigar::Del(1), Cigar::Match(10)],
            37,
        );
        assert_eq!(apply_heuristic_baq(&rec, Some((1109, 1113))), None);
    }
}
