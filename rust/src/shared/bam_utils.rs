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

/// The read's aligned blocks `[start, end)` on the reference: its reference
/// extent split at every splice (N). M, =, X and D extend a block; clips and
/// insertions add nothing. An unspliced read is one block, `[pos, ref_end)`.
pub fn aligned_blocks(record: &Record) -> Vec<(i64, i64)> {
    let mut blocks = Vec::with_capacity(2);
    let (mut start, mut pos) = (record.pos(), record.pos());
    for op in record.cigar().iter() {
        match op {
            Cigar::RefSkip(n) => {
                blocks.push((start, pos));
                pos += *n as i64;
                start = pos;
            }
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Del(n) => pos += *n as i64,
            _ => {}
        }
    }
    blocks.push((start, pos));
    blocks
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
/// fragment is well defined: paired, mate mapped on the same contig, the pair
/// facing inward (forward-reverse). TLEN runs from the forward read's 5' end to
/// the reverse read's 5' end (as BWA-MEM, samtools fixmate and Picard write it),
/// positive on the forward read, so past the mate's 5' end a read runs into
/// adapter (read-through). An outward-facing pair (the forward read's TLEN
/// negative, as at a tandem-duplication junction) defines no fragment. A forward read keeps every base up to
/// the last aligned one before its fragment end, a reverse read every base from
/// the first aligned one at or after its fragment start; an insertion beside the
/// boundary lies outside. Past an aligned end the boundary falls inside the soft
/// clip, counted base for base. None otherwise.
pub(crate) fn fragment_query_span(record: &Record) -> Option<(usize, usize)> {
    if !record.is_paired()
        || record.is_unmapped()
        || record.is_mate_unmapped()
        || record.tid() != record.mtid()
        || record.is_reverse() == record.is_mate_reverse()
        || (record.insert_size() > 0) == record.is_reverse()
        || record.insert_size() == 0
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

/// Most aligned bases past the fragment end that can still be adapter. An aligner
/// extends a read into adapter by a base or two of chance alignment (BWA-MEM
/// keeps an end mismatch over a clip): on MSK data one aligned base past the
/// boundary mismatches the reference 75% of the time, two 30%, ten or more 0.7%.
const ADAPTER_ALIGNED_MAX: usize = 2;

/// Whether the read's bases outside [lo, hi) look like adapter: soft-clipped, or
/// at most `ADAPTER_ALIGNED_MAX` of them aligned, and none inserted. More aligned
/// bases past the boundary, or an insertion there, are the molecule's: TLEN, a
/// reference distance, leaves out a molecule's inserted bases (a tandem
/// duplication) and a mate's clipped 5' bases.
fn adapter_like(record: &Record, lo: usize, hi: usize) -> bool {
    let (mut q, mut aligned) = (0usize, 0usize);
    for op in record.cigar().iter() {
        let n = op.len() as usize;
        let outside = n - ((q + n).min(hi).saturating_sub(q.max(lo)));
        match op {
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => aligned += outside,
            Cigar::Ins(_) if outside > 0 => return false,
            _ => {}
        }
        if matches!(op, Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) | Cigar::Ins(_) | Cigar::SoftClip(_)) {
            q += n;
        }
    }
    aligned <= ADAPTER_ALIGNED_MAX
}

/// A read ends at its fragment end: bases past it (read-through into adapter, an
/// insert shorter than the read) are neither the read's bases nor its reach.
/// Returns the read with those bases hard-clipped (removed from its sequence and
/// qualities, as if trimmed: a masked base would still fill a window as a
/// match), its start moved past any aligned bases it loses there, its tags kept;
/// the read itself when it has none, no well-defined fragment, or bases past the
/// boundary that are the molecule's rather than adapter (see `adapter_like`);
/// None when none of its aligned bases lies inside its fragment.
pub(crate) fn clip_to_fragment(mut record: Record) -> Option<Record> {
    let Some((lo, hi)) = fragment_query_span(&record) else {
        return Some(record);
    };
    let len = record.seq_len();
    if (lo == 0 && hi == len) || !adapter_like(&record, lo, hi) {
        return Some(record);
    }
    let (mut lead_hard, mut tail_hard) = (0u32, 0u32);
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
            Cigar::HardClip(n) if q == 0 => lead_hard = *n,
            Cigar::HardClip(n) => tail_hard = *n,
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
    // Bases outside [lo, hi) join the hard clips; soft clips inside it stay.
    let (lead_hard, tail_hard) = (lead_hard + lo as u32, tail_hard + (len - hi) as u32);
    let mut ops = Vec::with_capacity(mid.len() + 4);
    ops.extend((lead_hard > 0).then_some(Cigar::HardClip(lead_hard)));
    ops.extend((first_q > lo).then(|| Cigar::SoftClip((first_q - lo) as u32)));
    ops.extend(mid);
    ops.extend((last.0 < hi).then(|| Cigar::SoftClip((hi - last.0) as u32)));
    ops.extend((tail_hard > 0).then_some(Cigar::HardClip(tail_hard)));
    let quals = record.qual()[lo..hi].to_vec();
    let (qname, seq) = (record.qname().to_vec(), record.seq().as_bytes()[lo..hi].to_vec());
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
        // 51M 49S, the fragment 50 bases: one adapter base aligned, the rest clipped.
        let out = clip_to_fragment(read(vec![Cigar::Match(51), Cigar::SoftClip(49)], 1000, 50, false)).unwrap();
        assert_eq!((out.pos(), cigar(&out)), (1000, "50M50H".into()));
        assert_eq!((out.seq_len(), out.qual().len()), (50, 50));
        assert_eq!(out.seq().as_bytes(), (0..50).map(|i| b"ACGT"[i % 4]).collect::<Vec<u8>>());
        assert_eq!(out.aux(b"NH").unwrap(), Aux::U8(1));
    }

    #[test]
    fn a_reverse_read_before_its_fragment_start_is_clipped_and_moved() {
        // 58S 42M from 1058, the fragment [1060, 1100): two adapter bases aligned.
        let out = clip_to_fragment(read(vec![Cigar::SoftClip(58), Cigar::Match(42)], 1058, 40, true)).unwrap();
        assert_eq!((out.pos(), cigar(&out)), (1060, "60H40M".into()));
        assert_eq!(out.seq().as_bytes(), (60..100).map(|i| b"ACGT"[i % 4]).collect::<Vec<u8>>());
    }

    #[test]
    fn a_read_aligning_on_past_its_fragment_end_keeps_its_bases() {
        // Three or more aligned bases past TLEN's end are the molecule's (TLEN
        // understates it when the mate's 5' end is clipped); two are adapter.
        for (ops, kept) in [
            (vec![Cigar::Match(100)], true),
            (vec![Cigar::Match(53), Cigar::SoftClip(47)], true),
            (vec![Cigar::Match(52), Cigar::SoftClip(48)], false),
        ] {
            let before = cigar(&read(ops.clone(), 1000, 50, false));
            let out = clip_to_fragment(read(ops, 1000, 50, false)).unwrap();
            assert_eq!(cigar(&out) == before, kept, "{before} -> {}", cigar(&out));
        }
    }

    #[test]
    fn an_insertion_past_the_boundary_is_the_molecules() {
        // 40M 5I 55M with the fragment ending after the 40th aligned base: TLEN, a
        // reference distance, leaves out the molecule's inserted bases.
        let out = clip_to_fragment(read(vec![Cigar::Match(40), Cigar::Ins(5), Cigar::SoftClip(55)], 1000, 40, false)).unwrap();
        assert_eq!(cigar(&out), "40M5I55S");
    }

    #[test]
    fn a_read_inside_its_fragment_is_unchanged() {
        let out = clip_to_fragment(read(vec![Cigar::Match(100)], 1000, 300, false)).unwrap();
        assert_eq!((out.pos(), cigar(&out)), (1000, "100M".into()));
        let unpaired = {
            let mut r = read(vec![Cigar::Match(51), Cigar::SoftClip(49)], 1000, 50, false);
            r.set_flags(0);
            r
        };
        assert_eq!(cigar(&clip_to_fragment(unpaired).unwrap()), "51M49S");
    }

    #[test]
    fn an_outward_facing_pair_defines_no_fragment() {
        // The forward read's TLEN is negative (BWA, 5' to 5'): its mate lies
        // behind it, so nothing it reads is past a fragment end.
        let fwd = read(vec![Cigar::Match(51), Cigar::SoftClip(49)], 1000, -13, false);
        assert_eq!(fragment_query_span(&fwd), None);
        assert_eq!(cigar(&clip_to_fragment(fwd).unwrap()), "51M49S");
        let rev = read(vec![Cigar::Match(100)], 1000, -13, true); // TLEN +13 on the reverse read
        assert_eq!(fragment_query_span(&rev), None);
    }

    #[test]
    fn a_boundary_inside_a_deletion_ends_at_the_last_aligned_base_before_it() {
        // 30M 10D 2M 68S from 1000: the fragment ends at 1035, inside the deletion.
        let ops = vec![Cigar::Match(30), Cigar::Del(10), Cigar::Match(2), Cigar::SoftClip(68)];
        let out = clip_to_fragment(read(ops, 1000, 35, false)).unwrap();
        assert_eq!(cigar(&out), "30M70H");
    }

    #[test]
    fn a_boundary_inside_the_trailing_clip_keeps_the_molecules_clipped_bases() {
        // 80M 20S from 1000 with the fragment ending at 1090: ten clipped bases
        // are the molecule's, the last ten adapter.
        let out = clip_to_fragment(read(vec![Cigar::Match(80), Cigar::SoftClip(20)], 1000, 90, false)).unwrap();
        assert_eq!((cigar(&out), out.seq_len()), ("80M10S10H".into(), 90));
    }

    #[test]
    fn hard_clips_stay_outermost() {
        let ops = vec![Cigar::HardClip(5), Cigar::Match(51), Cigar::SoftClip(49), Cigar::HardClip(3)];
        assert_eq!(cigar(&clip_to_fragment(read(ops, 1000, 50, false)).unwrap()), "5H50M53H");
        let ops = vec![Cigar::HardClip(5), Cigar::SoftClip(58), Cigar::Match(42)];
        let out = clip_to_fragment(read(ops, 1058, 40, true)).unwrap();
        assert_eq!((out.pos(), cigar(&out)), (1060, "65H40M".into()));
    }

    #[test]
    fn a_clipped_read_lies_inside_its_fragment() {
        // Clipping twice changes nothing: every base left is the molecule's.
        let shapes = [
            (vec![Cigar::Match(51), Cigar::SoftClip(49)], 1000, 50, false),
            (vec![Cigar::SoftClip(58), Cigar::Match(42)], 1058, 40, true),
            (vec![Cigar::Match(30), Cigar::Del(10), Cigar::Match(2), Cigar::SoftClip(68)], 1000, 35, false),
            (vec![Cigar::SoftClip(60), Cigar::Match(1), Cigar::Del(10), Cigar::Match(39)], 1000, 39, true),
            (vec![Cigar::SoftClip(10), Cigar::Match(90)], 1000, 88, true),
            (vec![Cigar::Match(80), Cigar::SoftClip(20)], 1000, 90, false),
        ];
        for (ops, pos, tlen, reverse) in shapes {
            let once = clip_to_fragment(read(ops, pos, tlen, reverse)).unwrap();
            assert_eq!(fragment_query_span(&once), Some((0, once.seq_len())), "{}", cigar(&once));
        }
    }
}
