//! Heuristic Base Alignment Quality (BAQ) adjustment (Li 2011).
//!
//! Downgrades base qualities near CIGAR indels to reduce false positive
//! variant calls caused by alignment artifacts. Walks the CIGAR string
//! and subtracts `BAQ_PENALTY` from base qualities within `BAQ_RADIUS`
//! bp of any Ins/Del operation.
//!
//! Used by:
//! - `counting/engine.rs` — Phase 0 quality adjustment before classification
//! - `hla/router.rs` (future) — quality adjustment before PairHMM scoring

use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;

/// Number of base-pairs on either side of an indel to penalize.
const BAQ_RADIUS: usize = 5;
/// Quality penalty subtracted from BQ near indels (clamped to 0).
const BAQ_PENALTY: u8 = 20;

/// Apply heuristic BAQ quality downgrade to a read's base qualities.
///
/// Walks the CIGAR string, finds Ins/Del operations, and subtracts
/// `BAQ_PENALTY` (20) from base qualities within `BAQ_RADIUS` (5bp) of
/// the indel boundary on the read. Qualities are clamped to 0 (never negative).
///
/// Returns `None` if the read has no indels (caller should use original quals).
/// Returns `Some(adjusted_quals)` with the modified quality vector otherwise.
///
/// This function is intentionally allocation-free for reads without indels
/// (the common case), only allocating the Vec when adjustment is needed.
pub fn apply_heuristic_baq(record: &Record) -> Option<Vec<u8>> {
    let cigar = record.cigar();

    // Quick scan: does this read have any indels?
    let has_indel = cigar.iter().any(|op| matches!(op, Cigar::Ins(_) | Cigar::Del(_)));
    if !has_indel {
        return None;
    }

    let quals = record.qual().to_vec();
    let mut adjusted = quals;
    let read_len = adjusted.len();

    // Walk CIGAR to find read-coordinate positions of indel boundaries.
    // `read_pos` tracks the current position on the read (query) sequence.
    let mut read_pos: usize = 0;
    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                read_pos += *len as usize;
            }
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
                read_pos = ins_end;
            }
            Cigar::Del(_) | Cigar::RefSkip(_) => {
                // Deletion/skip: penalize BAQ_RADIUS bases on either side
                // of the deletion boundary on the read.
                // Deletion consumes reference but not query — read_pos stays.
                let pen_start = read_pos.saturating_sub(BAQ_RADIUS);
                let pen_end = (read_pos + BAQ_RADIUS).min(read_len);

                for q in adjusted[pen_start..pen_end].iter_mut() {
                    *q = q.saturating_sub(BAQ_PENALTY);
                }
                // Del/RefSkip don't advance read_pos (no query consumption)
            }
            Cigar::SoftClip(len) => {
                read_pos += *len as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {
                // No read or reference consumption
            }
        }
    }

    Some(adjusted)
}
