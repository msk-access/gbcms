//! BAM utility functions shared across analysis modes.
//!
//! Provides CIGAR-aware position lookup and quality computation helpers
//! that are independent of any variant-specific logic, and a bin's read cache
//! ([`ReadCache`]) with the runs of it each rule reads.
//!
//! Used throughout `counting` (re-exported by `counting::utils`).

use std::borrow::Cow;

use rust_htslib::bam::record::{Aux, Cigar, CigarString, CigarStringView};
use rust_htslib::bam::Record;

/// The read's CIGAR: the list the record keeps (`Record::cache_cigar`), borrowed, where
/// it keeps one; else a list decoded now, as `Record::cigar` decodes it. Every record of
/// a bin's read cache keeps its list ([`ReadCache::new`]), so the rules reading it decode
/// no CIGAR again.
pub fn cigar_of(record: &Record) -> Cow<'_, CigarStringView> {
    match record.cigar_cached() {
        Some(c) => {
            // A kept list is made at the record's final place and length; a record is not
            // moved or rebuilt after (`Record::set` drops the list it kept).
            debug_assert!(c.pos() == record.pos() && c.len() == record.cigar_len());
            Cow::Borrowed(c)
        }
        None => Cow::Owned(record.cigar()),
    }
}

/// Find the read index corresponding to a genomic position.
///
/// Walks the CIGAR string to translate a reference-coordinate position
/// into the corresponding index in the read's query sequence. Returns
/// `None` if the position falls in a deletion or is not covered by the read.
pub fn find_read_pos(record: &Record, target_pos: i64) -> Option<usize> {
    let cigar = cigar_of(record);
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
/// nothing. Read from the record's packed CIGAR, with no list decoded: the same sum as
/// `Record::cigar().end_pos()`. A CIGAR holding an op code past 8 is read by
/// `Record::cigar()` itself, which stops the run there, as it always has.
pub fn ref_end(record: &Record) -> i64 {
    packed_end(record.pos(), record.raw_cigar()).unwrap_or_else(|| record.cigar().end_pos())
}

/// The aligned end of a read at `pos` with the packed CIGAR `words`: `pos` plus every
/// reference-consuming op. None where a word holds an op code past 8 (rust-htslib
/// decodes none of those; nothing is decoded or stopped here).
fn packed_end(pos: i64, words: &[u32]) -> Option<i64> {
    let mut end = pos;
    for &op in words {
        match op & 0xf {
            // M, D, N, = and X take reference bases
            0 | 2 | 3 | 7 | 8 => end += i64::from(op >> 4),
            // I, S, H and P take none
            1 | 4 | 5 | 6 => {}
            _ => return None,
        }
    }
    Some(end)
}

/// The bounds a read cache keeps for a record at `pos` with the packed CIGAR `words`: its
/// aligned span `[pos, ref_end)`; for a CIGAR holding an op code past 8, which the cache
/// does not decode, the whole line [`UNREAD`].
fn cache_bounds(pos: i64, words: &[u32]) -> (i64, i64) {
    packed_end(pos, words).map_or(UNREAD, |end| (pos, end))
}

/// The bounds of a record whose CIGAR the cache cannot decode: they hold every interval, so
/// the record lies in every run and every rule meets it in cache order, as walking the whole
/// cache, and decodes it (and stops) where that walk would.
const UNREAD: (i64, i64) = (i64::MIN, i64::MAX);

/// Whether a read starting at `start` and ending (exclusive) at `end` lies outside the
/// interval `[lo, hi)`: it starts at or after `hi`, or ends at or before `lo`. For a rule
/// that reads both values first; its walk reads only what can meet the interval.
pub(crate) fn outside(start: i64, end: i64, (lo, hi): (i64, i64)) -> bool {
    start >= hi || end <= lo
}

/// [`outside`] for a record, its end read only where it starts before `hi`.
pub(crate) fn lies_outside(record: &Record, (lo, hi): (i64, i64)) -> bool {
    record.pos() >= hi || ref_end(record) <= lo
}

/// The records of a read cache that a rule reads for the genome interval `[lo, hi)`:
/// one run of them, in cache order, holding every record whose aligned span
/// `[pos, ref_end)` meets the interval (`pos < hi` and `ref_end > lo`). It may hold
/// others too, so each rule keeps its own test on every record of the run; a record
/// outside the run fails that test.
pub(crate) trait Reads {
    fn meeting(&self, lo: i64, hi: i64) -> &[Record];
}

/// A plain list of records is read whole, every record in its order.
impl Reads for [Record] {
    fn meeting(&self, _lo: i64, _hi: i64) -> &[Record] {
        self
    }
}

/// A plain list of records is read whole, every record in its order.
impl Reads for Vec<Record> {
    fn meeting(&self, _lo: i64, _hi: i64) -> &[Record] {
        self
    }
}

/// A bin's read cache: its records in the order its fetch returned them, each keeping
/// its decoded CIGAR, with two running bounds that give a rule the run of records that
/// can meet its interval without a look at the others.
///
/// - `reach[i]` is the largest aligned end among `records[..=i]`, so it never falls;
/// - `floor[i]` is the smallest start among `records[i..]`, so it never falls either.
///
/// A record before the first `i` with `reach[i] > lo` ends at or before `lo`; a record
/// from the first `i` (from there) with `floor[i] >= hi` on starts at or after `hi`.
/// Neither meets `[lo, hi)`. Both bounds are monotone however the records are ordered,
/// so the run is right for any order; the order sets only how long the run is.
///
/// A record whose CIGAR holds an op code past 8 is not decoded here: it keeps no list
/// and has the bounds [`UNREAD`], and every run holds it (`unread`, the first and last
/// such index, widens a run for the intervals at the line's ends too). Each rule then
/// meets it in cache order and decodes it (`Record::cigar`, which stops the run) exactly
/// where it would walking the whole cache.
pub(crate) struct ReadCache {
    records: Vec<Record>,
    reach: Vec<i64>,
    floor: Vec<i64>,
    unread: Option<(usize, usize)>,
}

impl ReadCache {
    /// The cache of `records`, in their order. A record keeps its CIGAR decoded from
    /// here on (once, at its final place: a record keeping a list made before it was
    /// moved or rebuilt decodes it again), unless its CIGAR holds an op code past 8.
    pub(crate) fn new(mut records: Vec<Record>) -> Self {
        let mut bounds = Vec::with_capacity(records.len());
        let mut unread: Option<(usize, usize)> = None;
        for (i, record) in records.iter_mut().enumerate() {
            let span = cache_bounds(record.pos(), record.raw_cigar());
            if span == UNREAD {
                unread = Some(unread.map_or((i, i), |(first, _)| (first, i)));
            } else {
                let kept = record.cigar_cached().is_some_and(|c| c.pos() == record.pos() && c.len() == record.cigar_len());
                if !kept {
                    record.cache_cigar();
                }
            }
            bounds.push(span);
        }
        let mut reach = Vec::with_capacity(records.len());
        let mut most = i64::MIN;
        for &(_, end) in &bounds {
            most = most.max(end);
            reach.push(most);
        }
        let mut floor = vec![0; records.len()];
        let mut least = i64::MAX;
        for (i, &(start, _)) in bounds.iter().enumerate().rev() {
            least = least.min(start);
            floor[i] = least;
        }
        Self { records, reach, floor, unread }
    }

    /// The run `[a, b)` of cache indices [`Reads::meeting`] reads for `[lo, hi)`: from the
    /// bounds, widened to hold every record the cache did not decode.
    pub(crate) fn run(&self, lo: i64, hi: i64) -> std::ops::Range<usize> {
        let a = self.reach.partition_point(|&end| end <= lo);
        let b = a + self.floor[a..].partition_point(|&start| start < hi);
        match self.unread {
            Some((first, last)) => a.min(first)..b.max(last + 1),
            None => a..b,
        }
    }
}

impl Reads for ReadCache {
    fn meeting(&self, lo: i64, hi: i64) -> &[Record] {
        &self.records[self.run(lo, hi)]
    }
}

/// The read's aligned blocks `[start, end)` on the reference: its reference
/// extent split at every splice (N). M, =, X and D extend a block; clips and
/// insertions add nothing. An unspliced read is one block, `[pos, ref_end)`.
pub fn aligned_blocks(record: &Record) -> Vec<(i64, i64)> {
    let mut blocks = Vec::with_capacity(2);
    let (mut start, mut pos) = (record.pos(), record.pos());
    for op in cigar_of(record).iter() {
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
    let cigar = cigar_of(record);
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
    if !fragment_defined(record) {
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
    for op in cigar_of(record).iter() {
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

/// Whether the read's fragment is well defined ([`fragment_query_span`]): paired, mate
/// mapped on the same contig, the pair facing inward, TLEN not 0.
fn fragment_defined(record: &Record) -> bool {
    record.is_paired()
        && !record.is_unmapped()
        && !record.is_mate_unmapped()
        && record.tid() == record.mtid()
        && record.is_reverse() != record.is_mate_reverse()
        && (record.insert_size() > 0) != record.is_reverse()
        && record.insert_size() != 0
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
    for op in cigar_of(record).iter() {
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
    // A read with a fragment has its CIGAR decoded once, here, where the span's walk
    // first reads it; the span, the adapter test and the clip borrow that list, and the
    // read keeps it (a clip rebuilding the read drops it).
    if fragment_defined(&record) && record.cigar_cached().is_none() {
        record.cache_cigar();
    }
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
    for op in cigar_of(&record).iter() {
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

/// An integer tag of any width, or None where it is absent or not an integer.
pub(crate) fn int_tag(record: &Record, tag: &[u8]) -> Option<i64> {
    match record.aux(tag) {
        Ok(Aux::I8(v)) => Some(i64::from(v)),
        Ok(Aux::U8(v)) => Some(i64::from(v)),
        Ok(Aux::I16(v)) => Some(i64::from(v)),
        Ok(Aux::U16(v)) => Some(i64::from(v)),
        Ok(Aux::I32(v)) => Some(i64::from(v)),
        Ok(Aux::U32(v)) => Some(i64::from(v)),
        _ => None,
    }
}

/// The aligner's own scores for the read's present placement: its alignment score (AS) and
/// the best score it found elsewhere (XS), each None where the tag is absent or not an
/// integer (a character XS, as STAR and HISAT2 write for the strand, is no score).
///
/// Both are None where the tags do not describe the present placement (RJ-29):
/// - the read carries a `YO` string beginning `N/A`: ABRA2 writes `YO:Z:N/A` (or
///   `YO:Z:N/A:<start>`) for a read it placed after BWA left it unmapped, and the tags stay
///   BWA's; or
/// - `AS` and `XS` are both present as integers and both 0: BWA's placeholder pair for an
///   unmapped read, which a realigner that writes no tag of its own leaves on the read.
///
/// Callers then fall back on the mapping quality, as for a read with no XS.
pub(crate) fn placement_scores(record: &Record) -> (Option<i64>, Option<i64>) {
    if let Ok(Aux::String(yo)) = record.aux(b"YO") {
        if yo.starts_with("N/A") {
            return (None, None);
        }
    }
    let (best, elsewhere) = (int_tag(record, b"AS"), int_tag(record, b"XS"));
    if best == Some(0) && elsewhere == Some(0) {
        return (None, None);
    }
    (best, elsewhere)
}

/// Seeded random reads for the tests of the read cache and of the rules that read it.
#[cfg(test)]
pub(crate) mod oracle {
    use super::*;

    /// The repo's deterministic LCG (as the binning property test).
    pub(crate) struct Lcg(u64);

    impl Lcg {
        pub(crate) fn new(seed: u64) -> Self {
            Self(seed)
        }

        /// A number in `[0, n)` (0 where `n` is 0).
        pub(crate) fn below(&mut self, n: u64) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (self.0 >> 33) % n.max(1)
        }

        /// True `k` times in `of`.
        pub(crate) fn chance(&mut self, k: u64, of: u64) -> bool {
            self.below(of) < k
        }
    }

    /// A random CIGAR: M I D N S H P = X, lengths 1 to 60 (a splice sometimes up to 50 kb),
    /// soft and hard clips at either end, sometimes no reference op at all.
    pub(crate) fn random_cigar(rng: &mut Lcg) -> Vec<Cigar> {
        let mut ops = Vec::new();
        if rng.chance(1, 6) {
            ops.push(Cigar::HardClip(1 + rng.below(20) as u32));
        }
        if rng.chance(1, 3) {
            ops.push(Cigar::SoftClip(1 + rng.below(60) as u32));
        }
        // mostly one aligned block, sometimes a few ops of any kind
        let middle = if rng.chance(1, 2) { 1 } else { rng.below(6) };
        for k in 0..middle {
            let len = 1 + rng.below(60) as u32;
            ops.push(if middle == 1 && k == 0 {
                Cigar::Match(30 + rng.below(121) as u32)
            } else {
                match rng.below(10) {
                    0..=2 => Cigar::Match(len),
                    3 => Cigar::Ins(len),
                    4 => Cigar::Del(len),
                    5 => Cigar::RefSkip(if rng.chance(1, 2) { 1 + rng.below(50_000) as u32 } else { len }),
                    6 => Cigar::Equal(len),
                    7 => Cigar::Diff(len),
                    8 => Cigar::Pad(len),
                    _ => Cigar::SoftClip(len),
                }
            });
        }
        if rng.chance(1, 3) {
            ops.push(Cigar::SoftClip(1 + rng.below(60) as u32));
        }
        if rng.chance(1, 6) {
            ops.push(Cigar::HardClip(1 + rng.below(20) as u32));
        }
        ops
    }

    /// The query bases a CIGAR holds (M, I, S, =, X).
    pub(crate) fn query_len(cigar: &[Cigar]) -> usize {
        cigar
            .iter()
            .filter(|op| matches!(op, Cigar::Match(_) | Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::Equal(_) | Cigar::Diff(_)))
            .map(|op| op.len() as usize)
            .sum()
    }

    /// A record named `name` of `cigar` at `pos`, its bases random (some N), qualities 2 to 41.
    pub(crate) fn record(rng: &mut Lcg, name: &[u8], cigar: Vec<Cigar>, pos: i64) -> Record {
        let len = query_len(&cigar);
        let seq: Vec<u8> = (0..len).map(|_| b"ACGTN"[if rng.chance(1, 50) { 4 } else { rng.below(4) as usize }]).collect();
        let quals: Vec<u8> = (0..len).map(|_| 2 + rng.below(40) as u8).collect();
        let mut rec = Record::new();
        rec.set(name, Some(&CigarString(cigar)), &seq, &quals);
        rec.set_pos(pos);
        rec.set_tid(0);
        rec.set_mapq(60);
        rec
    }

    /// Makes `rec` one read of a pair on contig 0: forward or reverse, its mate the other way
    /// round (mostly), TLEN random about its aligned span, sometimes 0 or the wrong sign.
    pub(crate) fn pair(rng: &mut Lcg, rec: &mut Record) {
        let reverse = rng.chance(1, 2);
        let mate_reverse = if rng.chance(1, 10) { reverse } else { !reverse };
        let mut flags = 0x1 | if rng.chance(1, 2) { 0x40 } else { 0x80 };
        if reverse {
            flags |= 0x10;
        }
        if mate_reverse {
            flags |= 0x20;
        }
        if rng.chance(1, 20) {
            flags |= 0x8;
        }
        rec.set_flags(flags);
        rec.set_mtid(if rng.chance(1, 20) { 1 } else { 0 });
        let span = ref_end(rec) - rec.pos();
        let tlen = match rng.below(10) {
            0 => 0,
            1..=4 => (span - rng.below(12) as i64).max(1),
            _ => span + rng.below(300) as i64,
        };
        let sign = if rng.chance(1, 10) { reverse } else { !reverse };
        rec.set_insert_size(if sign { tlen } else { -tlen });
    }

    /// A random cluster of `n` reads over about `[0, span)`: starts drawn with ties and
    /// kept in order, most reads one of a pair, each passed through `clip_to_fragment`
    /// as the bin's fetch passes it (a read the clip drops is left out); shuffled where
    /// `shuffle`.
    pub(crate) fn cluster(rng: &mut Lcg, n: usize, span: i64, shuffle: bool) -> Vec<Record> {
        let mut starts: Vec<i64> = (0..n).map(|_| rng.below(span.max(1) as u64) as i64).collect();
        starts.sort_unstable();
        let mut out = Vec::with_capacity(n);
        for (i, &pos) in starts.iter().enumerate() {
            let name = format!("r{}", i / 2);
            let cigar = random_cigar(rng);
            let mut rec = record(rng, name.as_bytes(), cigar, pos);
            if rng.chance(4, 5) {
                pair(rng, &mut rec);
            }
            if let Some(rec) = clip_to_fragment(rec) {
                out.push(rec);
            }
        }
        if shuffle {
            for i in (1..out.len()).rev() {
                out.swap(i, rng.below(i as u64 + 1) as usize);
            }
        }
        out
    }

    /// A random interval near `[0, span)`: an ordinary one, empty, reversed, one base, one
    /// outside every read, or the whole line either way round.
    pub(crate) fn interval(rng: &mut Lcg, span: i64) -> (i64, i64) {
        let at = rng.below(span.max(1) as u64 + 400) as i64 - 200;
        match rng.below(10) {
            0 => (at, at),
            1 => (at, at - 1 - rng.below(300) as i64),
            2 => (at, at + 1),
            3 => (-100_000, -99_000),
            4 => (span + 200_000, span + 201_000),
            5 => (i64::MIN, i64::MAX),
            6 => (i64::MAX, i64::MIN),
            _ => (at, at + 1 + rng.below(2_000) as i64),
        }
    }

    /// A random reference of `n` bases, with a short repeat now and then.
    pub(crate) fn reference(rng: &mut Lcg, n: usize) -> Vec<u8> {
        let mut out = Vec::with_capacity(n);
        while out.len() < n {
            if rng.chance(1, 12) {
                let unit: Vec<u8> = (0..1 + rng.below(3)).map(|_| b"ACGT"[rng.below(4) as usize]).collect();
                for _ in 0..2 + rng.below(5) {
                    out.extend_from_slice(&unit);
                }
            } else {
                out.push(b"ACGT"[rng.below(4) as usize]);
            }
        }
        out.truncate(n);
        out
    }

    /// `n` reads over `reference` (genome position 0 at its first base) at the change of
    /// `r` to `a` at `pos`: each spells REF or the ALT haplotype from a random start, is
    /// written as the aligner may write it (REF reads as M; ALT reads with the change as
    /// its own ops, or aligned up to it and soft-clipped past it, or clipped before it and
    /// aligned after), with a recurring change beside the event on some, an error now and
    /// then, a low quality here and there, both mates of some fragments, and a few reads of
    /// random shape. Fragments are named; reads of one fragment share the name.
    pub(crate) fn allele_reads(rng: &mut Lcg, reference: &[u8], pos: usize, r: &[u8], a: &[u8], n: usize) -> Vec<Record> {
        let change = Change::new(rng, reference, pos, r, a);
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let name = format!("f{}", if rng.chance(1, 4) && i > 0 { i - 1 } else { i });
            if rng.chance(1, 15) {
                let cigar = random_cigar(rng);
                let start = rng.below(reference.len() as u64) as i64;
                out.push(record(rng, name.as_bytes(), cigar, start));
                continue;
            }
            let alt = rng.chance(1, 2);
            let seq_of: &[u8] = if alt { &change.hap } else { reference };
            let len = (40 + rng.below(111) as usize).min(seq_of.len());
            let h0 = rng.below((seq_of.len() - len) as u64 + 1) as usize;
            out.push(change.read(rng, &name, alt, h0, len));
        }
        out
    }

    /// Reads of `len` bases over `reference` at the change of `r` to `a` at `pos`, one starting
    /// at every `step`-th base of the haplotype it spells (REF or the ALT, at random), written
    /// as [`allele_reads`] writes them, each its own fragment, in order of their start. All of
    /// one length, the largest end so far is about each read's own end, so a run holds little
    /// besides the reads that meet its interval: reads that end one base into an interval or
    /// start one base before its end, wherever it lies, are cut from a run one base too narrow.
    /// A quarter of them have their last (or first) aligned bases soft-clipped, so an aligned
    /// block ends (or starts) at any base, its bases read on past it; a quarter of the ALT reads
    /// are written without their gap ([`Change::read_as`]).
    pub(crate) fn tiled_reads(rng: &mut Lcg, reference: &[u8], pos: usize, r: &[u8], a: &[u8], len: usize, step: usize) -> Vec<Record> {
        let change = Change::new(rng, reference, pos, r, a);
        let mut out = Vec::new();
        for (i, h0) in (0..reference.len()).step_by(step.max(1)).enumerate() {
            let alt = rng.chance(1, 2);
            let seq_len = if alt { change.hap.len() } else { reference.len() };
            if h0 + len <= seq_len {
                let gapless = (alt && rng.chance(1, 4)).then(|| rng.chance(1, 2));
                let rec = change.read_as(rng, &format!("t{i}"), alt, h0, len, gapless);
                let k = 1 + rng.below(len as u64 / 2) as u32;
                let tail = rng.chance(1, 2);
                out.push(if rng.chance(1, 4) { clipped_in(&rec, k, tail).unwrap_or(rec) } else { rec });
            }
        }
        out.sort_by_key(Record::pos);
        out
    }

    /// `rec` with its last (`tail`) or first `k` aligned bases soft-clipped, where its CIGAR ends
    /// (starts) with an aligned block of more than `k` bases: its bases stay, its aligned end
    /// (start) moves `k` bases in. None otherwise.
    fn clipped_in(rec: &Record, k: u32, tail: bool) -> Option<Record> {
        let mut ops: Vec<Cigar> = rec.cigar().iter().copied().collect();
        let at = if tail { ops.len().checked_sub(1)? } else { 0 };
        let Cigar::Match(n) = *ops.get(at)? else { return None };
        if n <= k {
            return None;
        }
        let mut start = rec.pos();
        if tail {
            ops[at] = Cigar::Match(n - k);
            ops.push(Cigar::SoftClip(k));
        } else {
            ops[0] = Cigar::Match(n - k);
            ops.insert(0, Cigar::SoftClip(k));
            start += i64::from(k);
        }
        let mut out = Record::new();
        out.set(rec.qname(), Some(&CigarString(ops)), &rec.seq().as_bytes(), rec.qual());
        out.set_pos(start);
        out.set_tid(rec.tid());
        out.set_mapq(rec.mapq());
        out.set_flags(rec.flags());
        Some(out)
    }

    /// A change of `r` to `a` at `pos` on `reference`, with the ALT haplotype, and a base beside
    /// it that some reads of both alleles change.
    pub(crate) struct Change<'a> {
        reference: &'a [u8],
        pos: usize,
        r: &'a [u8],
        a: &'a [u8],
        pub(crate) hap: Vec<u8>,
        neighbour: usize,
    }

    impl<'a> Change<'a> {
        pub(crate) fn new(rng: &mut Lcg, reference: &'a [u8], pos: usize, r: &'a [u8], a: &'a [u8]) -> Self {
            let hap: Vec<u8> = [&reference[..pos], a, &reference[pos + r.len()..]].concat();
            let neighbour = pos + r.len() + 1 + rng.below(4) as usize;
            Self { reference, pos, r, a, hap, neighbour }
        }

        /// The read named `name` of `len` bases of REF (or of the ALT haplotype, where `alt`)
        /// from its base `h0`, written as the aligner may write it.
        pub(crate) fn read(&self, rng: &mut Lcg, name: &str, alt: bool, h0: usize, len: usize) -> Record {
            self.read_as(rng, name, alt, h0, len, None)
        }

        /// [`Change::read`]; where `gapless` is given, an ALT read is written as one aligned
        /// block, its gap left out (as an aligner writes a read ending in a repeat), placed from
        /// its first base (`Some(false)`: the ALT's bases past the change shifted against the
        /// reference) or from its last (`Some(true)`: those before it shifted).
        pub(crate) fn read_as(&self, rng: &mut Lcg, name: &str, alt: bool, h0: usize, len: usize, gapless: Option<bool>) -> Record {
            let (reference, pos, r, a, hap, neighbour) = (self.reference, self.pos, self.r, self.a, &self.hap, self.neighbour);
            let d = a.len() as i64 - r.len() as i64;
            // the bases REF and the ALT share before the change
            let lead = r.iter().zip(a).take_while(|(x, y)| x == y).count();
            let seq_of: &[u8] = if alt { hap } else { reference };
            let mut seq = seq_of[h0..h0 + len].to_vec();
            // a change beside the event on some reads of both alleles (more of the ALT's)
            let on = if alt { rng.chance(2, 5) } else { rng.chance(1, 5) };
            let nb = if alt { (neighbour as i64 + d) as usize } else { neighbour };
            if on && nb >= h0 && nb < h0 + len {
                seq[nb - h0] = if seq[nb - h0] == b'A' { b'C' } else { b'A' };
            }
            for b in seq.iter_mut() {
                if rng.chance(1, 100) {
                    *b = b"ACGT"[rng.below(4) as usize];
                }
            }
            let at = pos + lead; // the first base the alleles may differ at, in both
            let alt_end = pos + a.len(); // past the change, on the haplotype
            let (start, cigar) = if !alt || d == 0 || gapless.is_some() || h0 + len <= at || h0 >= alt_end {
                // aligned without a gap (REF, an MNP, or an ALT read off the change)
                // a gapless ALT read is placed from its first base, or from its last (the shift
                // then shows before the change)
                let from_right = alt && gapless == Some(true);
                let start = if alt && (h0 >= alt_end || from_right) { h0 as i64 - d } else { h0 as i64 };
                (start, vec![Cigar::Match(len as u32)])
            } else if h0 > at {
                // starting inside the inserted bases: clipped up to the right flank
                let clip = (alt_end - h0).min(len);
                let mut ops = vec![Cigar::SoftClip(clip as u32)];
                ops.extend((len > clip).then(|| Cigar::Match((len - clip) as u32)));
                ((pos + r.len()) as i64, ops)
            } else {
                let before = at - h0;
                let ins = ((a.len() - lead) as u32).min((len - before) as u32);
                let del = (r.len() - lead) as u32;
                let after = len - before - ins as usize;
                let mut ops = Vec::new();
                match rng.below(4) {
                    // the change as its own ops
                    0 | 1 => {
                        ops.extend((before > 0).then_some(Cigar::Match(before as u32)));
                        ops.extend((ins > 0).then_some(Cigar::Ins(ins)));
                        if after > 0 || before > 0 {
                            ops.extend((del > 0).then_some(Cigar::Del(del)));
                        }
                        ops.extend((after > 0).then_some(Cigar::Match(after as u32)));
                        if ops.is_empty() {
                            ops.push(Cigar::SoftClip(len as u32));
                        }
                        (h0 as i64, ops)
                    }
                    // aligned up to the change, clipped past it
                    2 => {
                        ops.extend((before > 0).then_some(Cigar::Match(before as u32)));
                        ops.push(Cigar::SoftClip((len - before) as u32));
                        (h0 as i64, ops)
                    }
                    // clipped before the change, aligned after it
                    _ => {
                        let clip = before + ins as usize;
                        ops.push(Cigar::SoftClip(clip as u32));
                        ops.extend((len > clip).then(|| Cigar::Match((len - clip) as u32)));
                        ((pos + r.len()) as i64, ops)
                    }
                }
            };
            let mut rec = Record::new();
            let quals: Vec<u8> = (0..len).map(|_| if rng.chance(1, 30) { 5 } else { 30 + rng.below(10) as u8 }).collect();
            rec.set(name.as_bytes(), Some(&CigarString(cigar)), &seq, &quals);
            rec.set_pos(start);
            rec.set_tid(0);
            rec.set_mapq(if rng.chance(1, 20) { 0 } else { 60 });
            rec.set_flags(if rng.chance(1, 2) { 0x10 } else { 0 });
            rec
        }
    }

    /// Writes an op code past 8 (9 to 15, which rust-htslib cannot decode) into one word of
    /// the record's packed CIGAR, through its data bytes (name and its NUL padding, CIGAR,
    /// bases, qualities: these records carry no tags), with no `unsafe`. False where the
    /// record has no CIGAR.
    pub(crate) fn unknown_op(rng: &mut Lcg, record: &mut Record) -> bool {
        let mut words = record.raw_cigar().to_vec();
        if words.is_empty() {
            return false;
        }
        let k = rng.below(words.len() as u64) as usize;
        words[k] = (words[k] & !0xf) | (9 + rng.below(7) as u32);
        let (qname, seq, qual) = (record.qname().to_vec(), record.seq().encoded.to_vec(), record.qual().to_vec());
        let mut data = qname.clone();
        let qlen = qname.len() + 1;
        data.resize(qlen + if qlen % 4 == 0 { 0 } else { 4 - qlen % 4 }, 0);
        for w in &words {
            data.extend(w.to_ne_bytes());
        }
        data.extend(&seq);
        data.extend(&qual);
        record.set_data(&data);
        assert_eq!(
            (record.qname(), record.raw_cigar(), record.seq().encoded, record.qual()),
            (&qname[..], &words[..], &seq[..], &qual[..]),
        );
        true
    }

    /// Whether the record's CIGAR holds an op code past 8.
    pub(crate) fn unreadable(r: &Record) -> bool {
        packed_end(r.pos(), r.raw_cigar()).is_none()
    }

    /// Where a run stopped: the panic's message, and the file and line it was raised at.
    pub(crate) type Stop = (String, String, u32);

    thread_local! {
        /// Whether this thread catches its stops ([`stopped`]), and where the last one was raised.
        static CATCHING: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
        static RAISED: std::cell::RefCell<Option<(String, u32)>> = const { std::cell::RefCell::new(None) };
    }

    /// Runs `f`, catching a stop with the place it was raised (the panic hook records it, on
    /// this thread only; every other panic goes to the hook before it, as without this).
    pub(crate) fn stopped<T>(f: impl FnOnce() -> T) -> Result<T, Stop> {
        static HOOK: std::sync::Once = std::sync::Once::new();
        HOOK.call_once(|| {
            let before = std::panic::take_hook();
            std::panic::set_hook(Box::new(move |info| {
                if CATCHING.with(std::cell::Cell::get) {
                    let at = info.location().map(|l| (l.file().to_string(), l.line()));
                    RAISED.with(|r| *r.borrow_mut() = at);
                } else {
                    before(info);
                }
            }));
        });
        CATCHING.with(|c| c.set(true));
        let out = std::panic::catch_unwind(std::panic::AssertUnwindSafe(f));
        CATCHING.with(|c| c.set(false));
        out.map_err(|e| {
            let (file, line) = RAISED.with(|r| r.borrow_mut().take()).unwrap_or_default();
            (panic_text(&e), file, line)
        })
    }

    /// A walk of a rule: its output, or where it stopped, and the index of each record its
    /// `counted` was asked about, in order.
    pub(crate) type Walk<T> = (Result<T, Stop>, Vec<usize>);

    /// Runs `rule` over `whole` (every record, in order, as a whole-cache walk reads them)
    /// and over `cache` (the same records), each time with a `counted` that notes the index
    /// of every record it is asked about and admits those of MAPQ above 0. (`whole` is a
    /// list, not a slice: a slice cannot be passed as `&dyn Reads`.)
    #[allow(clippy::ptr_arg)]
    pub(crate) fn both_walks<T>(
        whole: &Vec<Record>,
        cache: &ReadCache,
        rule: impl Fn(&dyn Reads, &dyn Fn(&Record) -> bool) -> T,
    ) -> [Walk<T>; 2] {
        let walk = |reads: &dyn Reads, records: &[Record]| {
            let asked = std::cell::RefCell::new(Vec::new());
            let counted = |r: &Record| {
                asked.borrow_mut().push(records.iter().position(|x| std::ptr::eq(x, r)).expect("a record of the walk"));
                r.mapq() > 0
            };
            let out = stopped(|| rule(reads, &counted));
            (out, asked.into_inner())
        };
        [walk(whole, whole), walk(cache, &cache.records)]
    }

    /// Asserts the two walks of [`both_walks`] agree: the same output, or the same stop (its
    /// message, file and line, rust-htslib's); the records the cache's walk asked about are,
    /// in order, among those the whole walk asked about; and where they stopped, both stopped
    /// at the same record (the last asked about). Whether they stopped.
    pub(crate) fn agree<T: PartialEq + std::fmt::Debug>([(whole, whole_asked), (run, run_asked)]: [Walk<T>; 2]) -> bool {
        assert_eq!(run, whole);
        if let Err((_, file, _)) = &whole {
            assert!(file.ends_with("bam/record.rs"), "stopped in {file}, not by rust-htslib");
        }
        let mut rest = whole_asked.iter();
        assert!(run_asked.iter().all(|i| rest.any(|j| j == i)), "the cache's walk asked about a record the whole walk did not, or out of order");
        if whole.is_err() {
            assert_eq!(run_asked.last(), whole_asked.last(), "stopped at another record");
        }
        whole.is_err()
    }

    /// The message a panic carried.
    pub(crate) fn panic_text(e: &Box<dyn std::any::Any + Send>) -> String {
        e.downcast_ref::<&str>().map(|s| s.to_string()).or_else(|| e.downcast_ref::<String>().cloned()).unwrap_or_default()
    }

    /// Whether the record's aligned span `[pos, ref_end)` meets `[lo, hi)`.
    pub(crate) fn meets(r: &Record, (lo, hi): (i64, i64)) -> bool {
        r.pos() < hi && ref_end(r) > lo
    }

    /// For a rule admitting the records `admits` takes, which every one meeting `interval`:
    /// the cache indices it admits walking the run of `cache`, against those it admits
    /// walking the whole cache, and the indices outside the run it would admit (none).
    pub(crate) fn admitted(
        cache: &ReadCache,
        interval: (i64, i64),
        admits: &dyn Fn(&Record) -> bool,
    ) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
        let run = cache.run(interval.0, interval.1);
        let whole: Vec<usize> = (0..cache.records.len()).filter(|&i| admits(&cache.records[i])).collect();
        let in_run: Vec<usize> = run.clone().filter(|&i| admits(&cache.records[i])).collect();
        let outside: Vec<usize> = whole.iter().copied().filter(|i| !run.contains(i)).collect();
        assert_eq!(cache.meeting(interval.0, interval.1).len(), run.len());
        (in_run, whole, outside)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::Aux;

    use super::oracle::{cluster, interval, meets, query_len, random_cigar, record, Lcg};

    /// Develop's fragment span, adapter test and clip, each decoding the CIGAR afresh at
    /// every read of it: the reference the one-decode clip is checked against.
    mod develop {
        use super::super::*;

        pub(super) fn develop_fragment_query_span(record: &Record) -> Option<(usize, usize)> {
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
            let (start, end) = (record.pos(), record.cigar().end_pos());
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


        fn develop_adapter_like(record: &Record, lo: usize, hi: usize) -> bool {
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

        pub(super) fn develop_clip_to_fragment(mut record: Record) -> Option<Record> {
            let Some((lo, hi)) = develop_fragment_query_span(&record) else {
                return Some(record);
            };
            let len = record.seq_len();
            if (lo == 0 && hi == len) || !develop_adapter_like(&record, lo, hi) {
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
    }

    #[test]
    fn the_end_read_from_the_packed_cigar_is_the_decoded_lists_end() {
        // Every op code 0 to 8 (M I D N S H P = X), lengths 0 to 2^20, any order.
        let mut rng = Lcg::new(0x5eed_0001);
        for _ in 0..5_000 {
            let ops: Vec<Cigar> = (0..rng.below(12))
                .map(|_| {
                    let len = rng.below(1 << 20) as u32 * u32::from(!rng.chance(1, 20));
                    match rng.below(9) {
                        0 => Cigar::Match(len),
                        1 => Cigar::Ins(len),
                        2 => Cigar::Del(len),
                        3 => Cigar::RefSkip(len),
                        4 => Cigar::SoftClip(len),
                        5 => Cigar::HardClip(len),
                        6 => Cigar::Pad(len),
                        7 => Cigar::Equal(len),
                        _ => Cigar::Diff(len),
                    }
                })
                .collect();
            let mut rec = Record::new();
            rec.set(b"r", Some(&CigarString(ops.clone())), &[], &[]);
            rec.set_pos(rng.below(1 << 31) as i64);
            assert_eq!(ref_end(&rec), rec.cigar().end_pos(), "{}", CigarString(ops.clone()));
            rec.cache_cigar();
            assert_eq!(ref_end(&rec), cigar_of(&rec).end_pos());
        }
    }

    #[test]
    fn a_kept_cigar_reads_as_the_list_decoded_afresh() {
        let mut rng = Lcg::new(0x5eed_0002);
        for trial in 0..200 {
            for kept in cluster(&mut rng, 40, 3_000, trial % 2 == 1) {
                // a clone keeps no list: it decodes afresh
                let fresh = kept.clone();
                assert!(fresh.cigar_cached().is_none());
                let mut kept = kept;
                if kept.cigar_cached().is_none() {
                    kept.cache_cigar();
                }
                assert!(matches!(cigar_of(&kept), Cow::Borrowed(_)));
                assert!(matches!(cigar_of(&fresh), Cow::Owned(_)));
                assert_eq!(*cigar_of(&kept), fresh.cigar());
                assert_eq!(*cigar_of(&fresh), fresh.cigar());
                assert_eq!(ref_end(&kept), fresh.cigar().end_pos());
                assert_eq!(aligned_blocks(&kept), aligned_blocks(&fresh));
                assert_eq!(soft_clips(&kept), soft_clips(&fresh));
                assert_eq!(fragment_query_span(&kept), develop::develop_fragment_query_span(&fresh));
                for _ in 0..10 {
                    let at = kept.pos() - 5 + rng.below(200) as i64;
                    assert_eq!(find_read_pos(&kept, at), find_read_pos(&fresh, at));
                }
            }
        }
    }

    #[test]
    fn the_clip_reading_one_decoded_cigar_gives_the_record_develops_clip_gives() {
        // Pairs of random shapes, TLENs and strands: the record, its start, CIGAR, bases and
        // qualities are develop's, and a list it keeps is its own CIGAR at its own start.
        let mut rng = Lcg::new(0x5eed_0003);
        let (mut clipped, mut dropped) = (0, 0);
        for i in 0..20_000 {
            let (cigar, pos) = (random_cigar(&mut rng), rng.below(5_000) as i64);
            let mut rec = record(&mut rng, format!("r{i}").as_bytes(), cigar, pos);
            super::oracle::pair(&mut rng, &mut rec);
            let before = rec.clone();
            let (got, want) = (clip_to_fragment(rec.clone()), develop::develop_clip_to_fragment(rec));
            match (got, want) {
                (Some(got), Some(want)) => {
                    clipped += usize::from(got.cigar().to_string() != before.cigar().to_string());
                    assert!(got == want, "{} at {}", before.cigar(), before.pos());
                    assert_eq!((got.pos(), got.cigar().to_string()), (want.pos(), want.cigar().to_string()));
                    assert_eq!((got.seq().as_bytes(), got.qual()), (want.seq().as_bytes(), want.qual()));
                    if let Some(c) = got.cigar_cached() {
                        assert_eq!(*c, want.cigar());
                    }
                }
                (None, None) => dropped += 1,
                (got, want) => panic!("{} at {}: {:?} against {:?}", before.cigar(), before.pos(), got.is_some(), want.is_some()),
            }
        }
        // the cases reach the clip, and the (rare) drop
        assert!(clipped > 1_000 && dropped > 0, "clipped {clipped}, dropped {dropped}");
    }

    #[test]
    fn the_clip_stops_on_an_unknown_op_code_where_develops_clip_stops() {
        // Pairs and single reads whose CIGAR holds an op code past 8: develop's clip decodes
        // it in its fragment span (`ref_end`, `Record::cigar`), this one where it keeps the
        // list (`Record::cache_cigar`); both stop on the same record with the same message
        // at the same rust-htslib line, and a read with no fragment passes both undecoded.
        let mut rng = Lcg::new(0x5eed_0008);
        let (mut stops, mut passed) = (0, 0);
        for i in 0..3_000 {
            let (cigar, pos) = (random_cigar(&mut rng), rng.below(5_000) as i64);
            let mut rec = record(&mut rng, format!("r{i}").as_bytes(), cigar, pos);
            if rng.chance(4, 5) {
                super::oracle::pair(&mut rng, &mut rec);
            }
            if !super::oracle::unknown_op(&mut rng, &mut rec) {
                continue;
            }
            let got = super::oracle::stopped(|| clip_to_fragment(rec.clone()));
            let want = super::oracle::stopped(|| develop::develop_clip_to_fragment(rec.clone()));
            match (got, want) {
                (Ok(got), Ok(want)) => {
                    assert!(got == want, "read {i}");
                    assert!(got.as_ref().is_some_and(|r| r.cigar_cached().is_none()), "read {i}: not decoded");
                    passed += 1;
                }
                (Err(got), Err(want)) => {
                    assert_eq!(got, want, "read {i}");
                    assert!(got.1.ends_with("bam/record.rs") && got.0 == "Unexpected cigar operation", "{got:?}");
                    stops += 1;
                }
                (got, want) => panic!("read {i}: {:?} against {:?}", got.is_ok(), want.is_ok()),
            }
        }
        assert!(stops > 500 && passed > 500, "{stops} stopped, {passed} passed");
    }

    #[test]
    fn a_run_holds_every_record_meeting_its_interval_in_any_order() {
        // Random clusters (clipped pairs, splices up to 50 kb, records with no reference
        // op), sorted and shuffled; random intervals, empty and reversed ones included: the
        // records meeting the interval are the run's, index for index, and none outside it.
        let mut rng = Lcg::new(0x5eed_0004);
        let (mut cases, mut met) = (0usize, 0usize);
        for trial in 0..400 {
            let n = rng.below(120) as usize;
            let records = cluster(&mut rng, n, 4_000, trial % 2 == 1);
            let cache = ReadCache::new(records);
            assert!(cache.records.iter().all(|r| r.cigar_cached().is_some_and(|c| *c == r.cigar())));
            for _ in 0..30 {
                let iv = interval(&mut rng, 4_000);
                let (in_run, whole, outside) = super::oracle::admitted(&cache, iv, &|r| meets(r, iv));
                assert_eq!(in_run, whole, "interval {iv:?}");
                assert!(outside.is_empty());
                cases += 1;
                met += whole.len();
            }
        }
        assert!(cases == 12_000 && met > 50_000, "{cases} cases, {met} records met");
    }

    #[test]
    fn a_run_is_cut_by_the_running_bounds() {
        // Three reads [0, 10), [5, 100), [20, 30) in that order: the long read raises the
        // running end, so [40, 50) reads from it on, and [0, 1) stops before it (it starts
        // at 5); nothing starts at or past 50.
        let mut rng = Lcg::new(1);
        let reads = [(0, 10), (5, 95), (20, 10)].map(|(pos, len)| record(&mut rng, b"r", vec![Cigar::Match(len)], pos));
        let cache = ReadCache::new(reads.into_iter().collect());
        assert_eq!(cache.run(40, 50), 1..3);
        assert_eq!(cache.run(0, 1), 0..1);
        assert_eq!(cache.run(100, 200), 3..3);
        assert_eq!(cache.run(-10, 0), 0..0);
        assert_eq!(query_len(&[Cigar::Match(3), Cigar::Del(2), Cigar::SoftClip(1)]), 4);
    }

    #[test]
    fn the_cache_bounds_are_the_span_or_the_whole_line() {
        // Packed words with every op code 0 to 15: an op code past 8 anywhere gives the whole
        // line; op codes 0 to 8 give [pos, end), the end the decoded list gives.
        let mut rng = Lcg::new(0x5eed_0005);
        let (mut unread, mut read) = (0, 0);
        for _ in 0..5_000 {
            let pos = rng.below(1 << 31) as i64;
            let codes = if rng.chance(1, 3) { 16 } else { 9 };
            let words: Vec<u32> = (0..rng.below(10)).map(|_| (rng.below(1 << 20) as u32) << 4 | rng.below(codes) as u32).collect();
            let bounds = cache_bounds(pos, &words);
            if words.iter().any(|w| w & 0xf > 8) {
                assert_eq!(bounds, UNREAD, "{words:?}");
                unread += 1;
            } else {
                let ops: Vec<Cigar> = words
                    .iter()
                    .map(|w| {
                        let n = w >> 4;
                        [Cigar::Match(n), Cigar::Ins(n), Cigar::Del(n), Cigar::RefSkip(n), Cigar::SoftClip(n), Cigar::HardClip(n), Cigar::Pad(n), Cigar::Equal(n), Cigar::Diff(n)][(w & 0xf) as usize]
                    })
                    .collect();
                let mut rec = Record::new();
                rec.set(b"r", Some(&CigarString(ops)), &[], &[]);
                rec.set_pos(pos);
                assert_eq!(rec.raw_cigar(), &words[..]);
                assert_eq!(bounds, (pos, rec.cigar().end_pos()), "{words:?}");
                read += 1;
            }
        }
        assert!(unread > 500 && read > 500, "{unread} unread, {read} read");
    }

    #[test]
    fn the_end_of_a_cigar_with_an_unknown_op_code_stops_where_the_decoder_does() {
        // ref_end hands such a CIGAR to Record::cigar: the same stop, with the same message, at
        // the same rust-htslib line.
        let mut rng = Lcg::new(0x5eed_0006);
        for i in 0..200 {
            let cigar = random_cigar(&mut rng);
            let mut rec = record(&mut rng, format!("r{i}").as_bytes(), cigar, 100);
            if !super::oracle::unknown_op(&mut rng, &mut rec) {
                continue;
            }
            let ours = super::oracle::stopped(|| ref_end(&rec)).expect_err("a stop");
            let theirs = super::oracle::stopped(|| rec.cigar().end_pos()).expect_err("a stop");
            assert_eq!(ours, theirs);
            assert!(ours.0 == "Unexpected cigar operation" && ours.1.ends_with("bam/record.rs"), "{ours:?}");
        }
    }

    #[test]
    fn a_run_holds_every_record_the_cache_cannot_decode() {
        // Random clusters with some records' CIGARs holding an op code past 8, sorted and
        // shuffled: the cache is built without decoding them (no stop, no kept list); every
        // run, for any interval (the line's ends included), holds each of them; the records it
        // can read keep their lists, and those meeting an interval are in its run.
        let mut rng = Lcg::new(0x5eed_0007);
        let (mut held, mut met) = (0usize, 0usize);
        for trial in 0..300 {
            let n = rng.below(80) as usize;
            let mut records = cluster(&mut rng, n, 4_000, trial % 2 == 1);
            for rec in records.iter_mut() {
                if rng.chance(1, 15) {
                    super::oracle::unknown_op(&mut rng, rec);
                }
            }
            let unread: Vec<usize> = (0..records.len()).filter(|&i| super::oracle::unreadable(&records[i])).collect();
            let cache = ReadCache::new(records);
            for (i, r) in cache.records.iter().enumerate() {
                assert_eq!(r.cigar_cached().is_some(), !unread.contains(&i));
            }
            for k in 0..30 {
                let iv = match k {
                    0 => (i64::MAX, i64::MAX),
                    1 => (i64::MIN, i64::MIN),
                    _ => interval(&mut rng, 4_000),
                };
                let run = cache.run(iv.0, iv.1);
                assert!(unread.iter().all(|i| run.contains(i)), "{iv:?}: {run:?} misses one of {unread:?}");
                held += unread.len();
                let meeting = |i: &usize| !unread.contains(i) && meets(&cache.records[*i], iv);
                let whole: Vec<usize> = (0..cache.records.len()).filter(meeting).collect();
                assert!(whole.iter().all(|i| run.contains(i)), "{iv:?}");
                met += whole.len();
            }
        }
        assert!(held > 2_000 && met > 10_000, "{held} held, {met} met");
    }

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

    /// A read carrying the given tags, for the placement tests.
    fn tagged(tags: Vec<(&[u8], Aux)>) -> Record {
        let mut rec = read(vec![Cigar::Match(10)], 100, 0, false);
        for (t, v) in tags {
            rec.push_aux(t, v).unwrap();
        }
        rec
    }

    #[test]
    fn placement_scores_none_for_a_realigner_placed_read() {
        let r = tagged(vec![(b"AS", Aux::I32(0)), (b"XS", Aux::I32(0)), (b"YO", Aux::String("N/A"))]);
        assert_eq!(placement_scores(&r), (None, None));
        let r = tagged(vec![(b"AS", Aux::I32(52)), (b"XS", Aux::I32(20)), (b"YO", Aux::String("N/A"))]);
        assert_eq!(placement_scores(&r), (None, None));
    }

    #[test]
    fn placement_scores_none_for_a_realigner_placed_read_with_a_start() {
        let r = tagged(vec![(b"AS", Aux::I32(0)), (b"XS", Aux::I32(0)), (b"YO", Aux::String("N/A:1234"))]);
        assert_eq!(placement_scores(&r), (None, None));
    }

    #[test]
    fn placement_scores_none_for_the_unmapped_placeholder_pair() {
        let r = tagged(vec![(b"AS", Aux::U8(0)), (b"XS", Aux::I16(0))]);
        assert_eq!(placement_scores(&r), (None, None));
    }

    #[test]
    fn placement_scores_keep_a_zero_score_without_an_xs() {
        let r = tagged(vec![(b"AS", Aux::I32(0))]);
        assert_eq!(placement_scores(&r), (Some(0), None));
    }

    #[test]
    fn placement_scores_keep_equal_nonzero_scores() {
        let r = tagged(vec![(b"AS", Aux::I32(52)), (b"XS", Aux::I32(52))]);
        assert_eq!(placement_scores(&r), (Some(52), Some(52)));
    }

    #[test]
    fn placement_scores_take_a_character_xs_as_no_score() {
        let r = tagged(vec![(b"AS", Aux::I32(43)), (b"XS", Aux::Char(b'+'))]);
        assert_eq!(placement_scores(&r), (Some(43), None));
    }

    #[test]
    fn placement_scores_keep_the_tags_for_a_yo_not_beginning_n_a() {
        let r = tagged(vec![(b"AS", Aux::I32(52)), (b"XS", Aux::I32(52)), (b"YO", Aux::String("5M1D"))]);
        assert_eq!(placement_scores(&r), (Some(52), Some(52)));
    }

    #[test]
    fn placement_scores_keep_a_zero_score_with_a_score_elsewhere() {
        // AS 0 with XS above 0 is not the placeholder pair: the aligner doubts the read
        let r = tagged(vec![(b"AS", Aux::I32(0)), (b"XS", Aux::I32(30))]);
        assert_eq!(placement_scores(&r), (Some(0), Some(30)));
    }

    #[test]
    fn placement_scores_keep_a_score_with_none_elsewhere() {
        // XS 0 with AS above 0 is not the placeholder pair: the aligner places the read here
        let r = tagged(vec![(b"AS", Aux::I32(52)), (b"XS", Aux::I32(0))]);
        assert_eq!(placement_scores(&r), (Some(52), Some(0)));
    }

    #[test]
    fn placement_scores_ignore_a_yo_that_is_not_a_string() {
        let r = tagged(vec![(b"AS", Aux::I32(52)), (b"XS", Aux::I32(52)), (b"YO", Aux::I32(0))]);
        assert_eq!(placement_scores(&r), (Some(52), Some(52)));
    }
}
