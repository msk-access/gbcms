//! BAM utility functions shared across analysis modes.
//!
//! Provides CIGAR-aware position lookup and quality computation helpers
//! that are independent of any variant-specific logic.
//!
//! Used throughout `counting` (re-exported by `counting::utils`).

use rust_htslib::bam::record::{Cigar, CigarString};
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

/// The read's query positions [lo, hi) that lie inside its fragment, when the
/// fragment is well defined: paired, mate mapped on the same contig in the
/// opposite orientation, TLEN set. TLEN runs from the forward read's 5' end to
/// the reverse read's 5' end (as BWA-MEM writes it), so past the mate's 5' end a
/// read runs into adapter (read-through). A forward read keeps every base up to
/// the last aligned one before its fragment end, a reverse read every base from
/// the first aligned one at or after its fragment start; an insertion beside the
/// boundary lies outside. Past an aligned end the boundary falls inside the soft
/// clip, counted base for base. None otherwise.
pub(crate) fn fragment_query_span(record: &Record) -> Option<(usize, usize)> {
    if !record.is_paired()
        || record.is_unmapped()
        || record.is_mate_unmapped()
        || record.insert_size() == 0
        || record.tid() != record.mtid()
        || record.is_reverse() == record.is_mate_reverse()
    {
        return None;
    }
    let len = record.seq_len();
    let tlen = record.insert_size().abs();
    let forward = !record.is_reverse();
    let (start, end) = (record.pos(), ref_end(record));
    let (frag_start, frag_end) = if forward { (start, start + tlen) } else { (end - tlen, end) };
    let (mut q, mut r) = (0usize, start);
    let (mut first_q, mut after_q) = (None, 0usize);
    let (mut hi_inside, mut lo_inside) = (0usize, None);
    for op in record.cigar().iter() {
        match op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                let n = *n as i64;
                first_q.get_or_insert(q);
                let before_end = (frag_end - r).clamp(0, n) as usize;
                if before_end > 0 {
                    hi_inside = q + before_end;
                }
                let before_start = (frag_start - r).clamp(0, n);
                if lo_inside.is_none() && before_start < n {
                    lo_inside = Some(q + before_start as usize);
                }
                q += n as usize;
                r += n;
                after_q = q;
            }
            Cigar::Ins(n) | Cigar::SoftClip(n) => q += *n as usize,
            Cigar::Del(n) | Cigar::RefSkip(n) => r += *n as i64,
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    let first_q = first_q?;
    if forward {
        let hi = if frag_end >= end { after_q + (frag_end - end) as usize } else { hi_inside };
        Some((0, hi.min(len)))
    } else {
        let lo = if frag_start <= start {
            first_q.saturating_sub((start - frag_start) as usize)
        } else {
            lo_inside.unwrap_or(after_q)
        };
        Some((lo, len))
    }
}

/// A read ends at its fragment end: bases past it (read-through into adapter, an
/// insert shorter than the read) are neither the read's bases nor its reach.
/// Returns the read with those bases soft-clipped and their qualities zeroed,
/// its start moved past any aligned bases it loses there, its tags kept; the read
/// itself when it has none or no well-defined fragment; None when none of its
/// aligned bases lies inside its fragment.
pub(crate) fn clip_to_fragment(mut record: Record) -> Option<Record> {
    let Some((lo, hi)) = fragment_query_span(&record) else {
        return Some(record);
    };
    let len = record.seq_len();
    if lo == 0 && hi == len {
        return Some(record);
    }
    let (mut lead_hard, mut tail_hard) = (None, None);
    let mut mid: Vec<Cigar> = Vec::new();
    let (mut q, mut r) = (0usize, record.pos());
    // The first kept aligned base (query, reference); where the last one ends
    // (query, ops kept up to it).
    let mut first: Option<(usize, i64)> = None;
    let mut last = (0usize, 0usize);
    for op in record.cigar().iter() {
        let n = op.len() as usize;
        let (a, b) = (q.max(lo), (q + n).min(hi)); // the op's query bases inside
        match op {
            Cigar::HardClip(_) if q == 0 => lead_hard = Some(*op),
            Cigar::HardClip(_) => tail_hard = Some(*op),
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
                if a < b {
                    first.get_or_insert((a, r + (a - q) as i64));
                    let k = (b - a) as u32;
                    mid.push(match op {
                        Cigar::Equal(_) => Cigar::Equal(k),
                        Cigar::Diff(_) => Cigar::Diff(k),
                        _ => Cigar::Match(k),
                    });
                    last = (b, mid.len());
                }
                q += n;
                r += n as i64;
            }
            // An insertion before the first kept aligned base joins the soft clip.
            Cigar::Ins(_) => {
                if a < b && first.is_some() {
                    mid.push(Cigar::Ins((b - a) as u32));
                }
                q += n;
            }
            Cigar::SoftClip(_) => q += n,
            Cigar::Del(_) | Cigar::RefSkip(_) => {
                if first.is_some() {
                    mid.push(*op);
                }
                r += n as i64;
            }
            Cigar::Pad(_) => {}
        }
    }
    let (first_q, pos) = first?;
    mid.truncate(last.1); // no insertion, deletion or skip after the last kept base
    let mut ops = Vec::with_capacity(mid.len() + 4);
    ops.extend(lead_hard);
    if first_q > 0 {
        ops.push(Cigar::SoftClip(first_q as u32));
    }
    ops.extend(mid);
    if last.0 < len {
        ops.push(Cigar::SoftClip((len - last.0) as u32));
    }
    ops.extend(tail_hard);
    let mut quals = record.qual().to_vec();
    quals[..lo].fill(0);
    quals[hi..].fill(0);
    let (qname, seq) = (record.qname().to_vec(), record.seq().as_bytes());
    record.set(&qname, Some(&CigarString(ops)), &seq, &quals);
    record.set_pos(pos);
    Some(record)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::Aux;

    /// A read of `cigar` at `pos` in a pair with fragment length `tlen`, carrying
    /// an NH tag; forward unless `reverse`.
    fn read(cigar: Vec<Cigar>, pos: i64, tlen: i64, reverse: bool) -> Record {
        let cigar = CigarString(cigar);
        let len = cigar.iter().filter(|op| matches!(op, Cigar::Match(_) | Cigar::Ins(_) | Cigar::SoftClip(_))).map(|op| op.len() as usize).sum();
        let mut rec = Record::new();
        let seq: Vec<u8> = (0..len).map(|i| b"ACGT"[i % 4]).collect();
        rec.set(b"r", Some(&cigar), &seq, &vec![30; len]);
        rec.set_pos(pos);
        rec.set_tid(0);
        rec.set_mtid(0);
        rec.set_flags(if reverse { 0x1 | 0x10 } else { 0x1 | 0x20 });
        rec.set_insert_size(if reverse { -tlen } else { tlen });
        rec.push_aux(b"NH", Aux::U8(1)).unwrap();
        rec
    }

    fn cigar(rec: &Record) -> String {
        rec.cigar().to_string()
    }

    #[test]
    fn a_forward_read_past_its_fragment_end_is_clipped_there() {
        let out = clip_to_fragment(read(vec![Cigar::Match(100)], 1000, 50, false)).unwrap();
        assert_eq!((out.pos(), cigar(&out)), (1000, "50M50S".into()));
        assert!(out.qual()[..50].iter().all(|&q| q == 30));
        assert!(out.qual()[50..].iter().all(|&q| q == 0));
        assert_eq!(out.aux(b"NH").unwrap(), Aux::U8(1));
    }

    #[test]
    fn a_reverse_read_before_its_fragment_start_is_clipped_and_moved() {
        let out = clip_to_fragment(read(vec![Cigar::Match(100)], 1000, 40, true)).unwrap();
        assert_eq!((out.pos(), cigar(&out)), (1060, "60S40M".into()));
        assert!(out.qual()[..60].iter().all(|&q| q == 0));
        assert!(out.qual()[60..].iter().all(|&q| q == 30));
    }

    #[test]
    fn a_read_inside_its_fragment_is_unchanged() {
        let out = clip_to_fragment(read(vec![Cigar::Match(100)], 1000, 300, false)).unwrap();
        assert_eq!((out.pos(), cigar(&out)), (1000, "100M".into()));
        let unpaired = {
            let mut r = read(vec![Cigar::Match(100)], 1000, 50, false);
            r.set_flags(0);
            r
        };
        assert_eq!(cigar(&clip_to_fragment(unpaired).unwrap()), "100M");
    }

    #[test]
    fn a_boundary_inside_a_deletion_ends_at_the_last_aligned_base_before_it() {
        // 30M 10D 70M from 1000: the fragment ends at 1035, inside the deletion.
        let out = clip_to_fragment(read(vec![Cigar::Match(30), Cigar::Del(10), Cigar::Match(70)], 1000, 35, false)).unwrap();
        assert_eq!(cigar(&out), "30M70S");
    }

    #[test]
    fn an_insertion_beside_the_boundary_lies_outside() {
        // 40M 5I 55M: the fragment ends after the 40th aligned base.
        let out = clip_to_fragment(read(vec![Cigar::Match(40), Cigar::Ins(5), Cigar::Match(55)], 1000, 40, false)).unwrap();
        assert_eq!(cigar(&out), "40M60S");
        assert!(out.qual()[40..].iter().all(|&q| q == 0));
    }

    #[test]
    fn a_boundary_inside_the_trailing_clip_masks_only_past_it() {
        // 80M 20S from 1000 with the fragment ending at 1090: ten clipped bases
        // are the molecule's, the last ten adapter.
        let out = clip_to_fragment(read(vec![Cigar::Match(80), Cigar::SoftClip(20)], 1000, 90, false)).unwrap();
        assert_eq!(cigar(&out), "80M20S");
        assert!(out.qual()[..90].iter().all(|&q| q == 30));
        assert!(out.qual()[90..].iter().all(|&q| q == 0));
    }

    #[test]
    fn hard_clips_stay_outermost() {
        let out = clip_to_fragment(read(vec![Cigar::HardClip(5), Cigar::Match(100), Cigar::HardClip(3)], 1000, 50, false)).unwrap();
        assert_eq!(cigar(&out), "5H50M50S3H");
    }
}
