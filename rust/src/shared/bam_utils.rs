//! BAM utility functions shared across analysis modes.
//!
//! Provides CIGAR-aware position lookup and quality computation helpers
//! that are independent of any variant-specific logic.
//!
//! Used throughout `counting` (re-exported by `counting::utils`).

use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;

/// Find the read index corresponding to a genomic position.
///
/// Walks the CIGAR string to translate a reference-coordinate position
/// into the corresponding index in the read's query sequence. Returns
/// `None` if the position falls in a deletion or is not covered by the read.
pub fn find_read_pos(record: &Record, target_pos: i64) -> Option<usize> {
    let cigar = record.cigar();
    let mut ref_pos = record.pos();
    let mut read_pos = 0;

    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                if target_pos >= ref_pos && target_pos < ref_pos + *len as i64 {
                    return Some(read_pos + (target_pos - ref_pos) as usize);
                }
                ref_pos += *len as i64;
                read_pos += *len as usize;
            }
            Cigar::Ins(len) => {
                read_pos += *len as usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                if target_pos >= ref_pos && target_pos < ref_pos + *len as i64 {
                    return None; // Position is deleted
                }
                ref_pos += *len as i64;
            }
            Cigar::SoftClip(len) => {
                read_pos += *len as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    None
}

/// End (exclusive, 0-based) of the read's aligned reference extent: its start
/// plus every reference-consuming op (M, =, X, D, N); clips and insertions add
/// nothing.
pub fn ref_end(record: &Record) -> i64 {
    record.cigar().end_pos()
}

/// Lengths of the read's leading and trailing soft clips, behind any hard clip.
/// A CIGAR of one (non-hard-clip) op has no trailing clip.
pub fn soft_clips(record: &Record) -> (u32, u32) {
    let cigar = record.cigar();
    let mut ops = cigar.iter().filter(|op| !matches!(op, Cigar::HardClip(_)));
    let clip = |op: Option<&Cigar>| match op {
        Some(Cigar::SoftClip(n)) => *n,
        _ => 0,
    };
    let lead = clip(ops.next());
    (lead, clip(ops.next_back()))
}

/// Compute the median quality of bases that pass the minimum threshold.
///
/// Follows the GATK `BaseQuality` annotation standard (median rather than min)
/// to prevent a single low-quality outlier from penalizing an entire read's
/// contribution to fragment consensus.
///
/// Returns 0 if no qualifying bases.
///
/// Even-count convention: returns the *upper* of the two middle values
/// (`filtered[len / 2]`), a slight upward bias, rather than averaging them. This
/// is deliberate — `med_qual` feeds classification thresholds in `variant_checks`,
/// `alignment`, and `pairhmm`, so switching to an exact average could flip
/// borderline REF/ALT calls for no diagnostic gain. Kept as-is intentionally.
#[inline]
pub fn median_qual(quals: &[u8], min_baseq: u8) -> u8 {
    let mut filtered: Vec<u8> = quals.iter()
        .copied()
        .filter(|&q| q >= min_baseq)
        .collect();
    if filtered.is_empty() { return 0; }
    filtered.sort_unstable();
    filtered[filtered.len() / 2]
}
