//! Exact-carrier classification for complex variants.
//!
//! A complex variant (a delins, a deletion whose anchor also changes, or an MNP
//! read with an indel in or right beside its block) is REF or ALT for a read only when the read's own
//! bases carry that allele across the whole event, with [`FLANK`] reference bases
//! on each side. The read's bases include soft clips; bases below min BQ (and N)
//! match anything, the one quality rule every backend shares. Nothing else is
//! tolerated: a read one confident base off the given allele carries a different
//! allele, and a read that ends inside the event cannot show either.
//!
//! - **The event** is every base where REF and ALT can differ over all equivalent
//!   placements. That is the union of the minimal difference trimmed from the
//!   left first and from the right first, grown through any tandem repeat
//!   touching either end on either allele (unit 1-6, or the length change):
//!   an aligner may place an indel anywhere in such a repeat, including a run
//!   the ALT's own bases continue.
//! - **Windows are read where they sit.** Each window is compared with the read's
//!   bases starting at the query position of the window's first reference base
//!   (left anchor), and ending at its last (right anchor); the closer reading
//!   counts. The anchors are flank bases outside the event, aligned the same way
//!   on either allele; within the event the read's bases are compared as a
//!   contiguous stretch, so where the aligner put an indel (or a clip) does not
//!   matter, and one it placed past one anchor leaves the other in place. A copy
//!   of a window elsewhere in the read is never a match.
//! - **Equal windows.** The shorter allele's window is padded with reference
//!   flank until both are the same length, so neither allele is favoured by where
//!   reads start.
//! - **Long events.** When the windows exceed [`LONG_EVENT`] bases no read can
//!   hold them whole. Both alleles are then judged by equal-length junction
//!   windows at each end, reading inward: a flank through one base past the first
//!   base where they differ, and through the shorter allele when that fits in
//!   [`LONG_EVENT`] bases (the right junction from just before its first
//!   difference, so growth on the left favours neither allele). A read must match
//!   the same allele at every junction it holds.
//! - **Outcomes.** A read decides only where it holds both the REF and the ALT
//!   window (padding can let one complete where the other cannot). A read
//!   matching both (through masked bases) is neither. A read holding no pair is
//!   depth only, like a pure-indel read ending inside its tract: no allele, no
//!   partial evidence, no mFSD class. A read holding the windows and matching
//!   neither carries another allele: partial evidence when closer to ALT.
//! - **Spliced reads.** A read's splice in a window's flank or padding cuts both
//!   alleles' windows at the exon edge (the read's next bases come from the next
//!   exon), anchored at the junction. A splice through the bases where the
//!   alleles differ leaves the read depth only.

use rust_htslib::bam::record::{Cigar, Record};

use super::utils::{ClassifyPhase, ClassifyResult};
use crate::types::Variant;

/// Reference bases required on each side of the event.
const FLANK: usize = 2;
/// Window length beyond which the whole allele cannot be expected in one read.
const LONG_EVENT: usize = 50;
/// Longest tandem-repeat unit the event is grown through (besides the length change).
const MAX_UNIT: usize = 6;
/// Masked base: matches any haplotype base.
const WILD: u8 = b'N';

/// A haplotype string over a reference span, read from the read at one or both
/// of its reference ends.
struct Window {
    seq: Vec<u8>,
    /// Reference position of the window's first base, when it can anchor on it.
    left: Option<i64>,
    /// Reference position one past the window's last base, when it can anchor there.
    right: Option<i64>,
    /// Bases of the window before the event (growth, flank and padding) and
    /// after it, counted from where the alleles differ.
    before: usize,
    after: usize,
}

/// The REF and ALT windows, paired: whole-event windows (one pair) or junction
/// windows (a left pair and a right pair).
struct Windows {
    pairs: Vec<(Window, Window)>,
    /// Where the alleles differ, genomic [start, end): the trims' union, before
    /// repeat growth. A spliced read never holds growth past its junction.
    event: (i64, i64),
}

impl Windows {
    /// The windows as a read spliced at `skips` sees them. A splice in a window's
    /// flank or padding ends the window at the exon edge on that side: the read's
    /// next bases come from the next exon, and it is anchored at the junction.
    /// Both alleles' windows are cut alike.
    fn cut_at(&self, skips: &[(i64, i64)]) -> Windows {
        let (ev_lo, ev_hi) = self.event;
        let exon_start = skips.iter().filter(|s| s.1 <= ev_lo).map(|s| s.1).max();
        let exon_end = skips.iter().filter(|s| s.0 >= ev_hi).map(|s| s.0).min();
        let cut = |w: &Window| -> Window {
            let (mut seq, mut left, mut right, mut before, mut after) = (w.seq.clone(), w.left, w.right, w.before, w.after);
            if let Some(n1) = exon_start.filter(|&n1| n1 > ev_lo - before as i64) {
                let drop = (n1 - (ev_lo - before as i64)) as usize;
                seq.drain(..drop);
                before -= drop;
                left = left.map(|_| n1);
            }
            if let Some(n0) = exon_end.filter(|&n0| n0 < ev_hi + after as i64) {
                let drop = (ev_hi + after as i64 - n0) as usize;
                seq.truncate(seq.len() - drop);
                after -= drop;
                right = right.map(|_| n0);
            }
            Window { seq, left, right, before, after }
        };
        Windows { pairs: self.pairs.iter().map(|(r, a)| (cut(r), cut(a))).collect(), event: self.event }
    }

    /// The reference span the windows cover.
    fn span(&self) -> (i64, i64) {
        self.pairs.iter().fold((i64::MAX, i64::MIN), |(lo, hi), (r, a)| {
            let starts = [r.left, a.left].into_iter().flatten();
            let ends = [r.right, a.right].into_iter().flatten();
            (starts.fold(lo, i64::min), ends.fold(hi, i64::max))
        })
    }
}

/// A window read from the read: mismatching bases, masked bases (below min BQ,
/// or N; they match anything), whether an N sat among them, and the query span.
struct Reading {
    mismatches: usize,
    masked: usize,
    had_n: bool,
    span: (usize, usize),
}

/// Classify a read at a complex variant by the exact-carrier rule, or None when
/// the rule cannot judge this variant: it needs both alleles non-empty and a
/// reference (prep's `event_ref`, else the `ref_context`) holding the event, its
/// flank and padding. Prep fetches enough for that, so None means an unprepared
/// variant or an event at a contig end; the caller then uses the previous
/// classifier.
pub(crate) fn classify(record: &Record, variant: &Variant, quals: &[u8], min_baseq: u8) -> Option<ClassifyResult> {
    let win = windows(variant)?;
    let seq = record.seq().as_bytes();
    if seq.is_empty() || quals.len() < seq.len() {
        // A record stored without its bases (SEQ '*') shows nothing here.
        return Some(ClassifyResult::no_coverage(ClassifyPhase::MaskedCompare));
    }
    let skips = ref_skips(record, win.span());
    let (ev_lo, ev_hi) = win.event;
    if skips.iter().any(|&(n0, n1)| n0 < ev_hi && n1 > ev_lo) {
        // Spliced through the event: the read skips bases where the alleles
        // differ, so it shows neither (depth only).
        let mut r = ClassifyResult::neither(ClassifyPhase::CigarRecon);
        r.ref_uninformative = true;
        return Some(r);
    }
    let cut;
    let win = if skips.is_empty() {
        &win
    } else {
        cut = win.cut_at(&skips);
        &cut
    };

    // Per pair the read holds: (matches REF, matches ALT).
    let mut held: Vec<(bool, bool)> = Vec::with_capacity(win.pairs.len());
    let (mut mm_ref, mut mm_alt, mut alt_masked, mut had_n) = (0usize, 0usize, 0usize, false);
    let mut qual_bases: Vec<u8> = Vec::new();
    let mut read_spans: Vec<(usize, usize)> = Vec::new();
    for (rw, aw) in &win.pairs {
        let (Some(r), Some(a)) = (read_window(record, &seq, quals, min_baseq, rw), read_window(record, &seq, quals, min_baseq, aw))
        else {
            continue;
        };
        mm_ref += r.mismatches;
        mm_alt += a.mismatches;
        alt_masked += a.masked;
        had_n |= r.had_n || a.had_n;
        for (lo, hi) in [r.span, a.span] {
            qual_bases.extend_from_slice(&quals[lo..hi]);
            read_spans.push((lo, hi));
        }
        held.push((r.mismatches == 0, a.mismatches == 0));
    }
    let qual = median(&qual_bases, min_baseq);
    let every = |want: (bool, bool)| !held.is_empty() && held.iter().all(|&h| h == want);
    let mut result = if every((false, true)) {
        ClassifyResult::is_alt(qual, ClassifyPhase::MaskedCompare)
    } else if every((true, false)) {
        ClassifyResult::is_ref(qual, ClassifyPhase::MaskedCompare)
    } else if !held.is_empty() && mm_alt < mm_ref {
        // Another allele, closer to ALT: partial evidence.
        ClassifyResult::neither_with_nearby(qual, ClassifyPhase::MaskedCompare)
    } else {
        ClassifyResult::neither(ClassifyPhase::MaskedCompare)
    };
    // A read holding no pair cannot show the allele: depth only (as a pure-indel
    // read that ends inside its tract), and no mFSD class.
    result.ref_uninformative = held.is_empty();
    result.has_n_base = had_n;
    // An MNP ALT read with every window base read unmasked shows the whole
    // haplotype, as a fully read block does on the base-by-base path.
    result.mnp_confirmed =
        result.is_alt && alt_masked == 0 && variant.ref_allele.len() == variant.alt_allele.len();
    // A read decided from its own bases may count although its aligned span stops
    // short of the variant (its allele lies in soft-clipped bases), when every base
    // read lies inside its fragment: past the fragment end a clip is adapter.
    result.clip_admissible = (result.is_ref || result.is_alt)
        && fragment_query_span(record).is_some_and(|(lo, hi)| read_spans.iter().all(|&(a, b)| a >= lo && b <= hi));
    Some(result)
}

/// The read's query positions [lo, hi) that lie inside its fragment, when the
/// fragment is well defined: paired, mate mapped on the same contig in the
/// opposite orientation, TLEN set. Past the mate's 5' end the read runs into
/// adapter (read-through). None otherwise.
fn fragment_query_span(record: &Record) -> Option<(usize, usize)> {
    if !record.is_paired()
        || record.is_mate_unmapped()
        || record.insert_size() == 0
        || record.tid() != record.mtid()
        || record.is_reverse() == record.is_mate_reverse()
    {
        return None;
    }
    let len = record.seq_len() as i64;
    let tlen = record.insert_size().abs();
    let (mut ref_end, mut first_q, mut clip_tail) = (record.pos(), None, 0i64);
    let mut q = 0i64;
    for op in record.cigar().iter() {
        match op {
            Cigar::SoftClip(n) => {
                if first_q.is_none() {
                    q += *n as i64;
                } else {
                    clip_tail += *n as i64;
                }
            }
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                first_q.get_or_insert(q);
                q += *n as i64;
                ref_end += *n as i64;
            }
            Cigar::Ins(n) => q += *n as i64,
            Cigar::Del(n) | Cigar::RefSkip(n) => ref_end += *n as i64,
            _ => {}
        }
    }
    let first_q = first_q?;
    let after_q = len - clip_tail; // one past the last aligned base
    if !record.is_reverse() {
        // Forward: the fragment starts with the read and runs TLEN bases.
        let frag_end = record.pos() + tlen;
        let hi = if frag_end >= ref_end {
            after_q + (frag_end - ref_end)
        } else {
            qpos(record, frag_end - 1).map_or(after_q, |q| q as i64 + 1)
        };
        Some((0, hi.clamp(0, len) as usize))
    } else {
        // Reverse: the fragment ends with the read and starts TLEN bases earlier.
        let frag_start = ref_end - tlen;
        let lo = if frag_start <= record.pos() {
            first_q - (record.pos() - frag_start)
        } else {
            qpos(record, frag_start).map_or(first_q, |q| q as i64)
        };
        Some((lo.clamp(0, len) as usize, len as usize))
    }
}

/// Whether the read has an insertion or deletion in the variant's block or right
/// beside it. An aligner may write an MNP that way (a block shifted by one base
/// as an insertion before it and a deletion after), so such a read is judged by
/// its own bases, like any MNP read with an indel in the block. An indel further
/// off (a germline one nearby) leaves the block's bases where IGV shows them.
pub(crate) fn indel_at_block(record: &Record, variant: &Variant) -> bool {
    let (start, end) = (variant.pos, variant.pos + variant.ref_allele.len() as i64);
    let mut pos = record.pos();
    for op in record.cigar().iter() {
        match op {
            // An insertion sits between reference bases pos-1 and pos.
            Cigar::Ins(_) if pos >= start && pos <= end => return true,
            Cigar::Del(len) => {
                let d_end = pos + *len as i64;
                if pos <= end && d_end >= start {
                    return true;
                }
                pos = d_end;
            }
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) | Cigar::RefSkip(len) => pos += *len as i64,
            _ => {}
        }
    }
    false
}

/// Whether a reference fetched from `start` is too short to judge the variant:
/// (on the left, on the right), for prep to fetch more. None when it holds the
/// event with its flank and padding, or the variant is not judged by this rule.
pub(crate) fn reference_short(start: i64, reference: &str, v: &Variant) -> Option<(bool, bool)> {
    let (r, a) = (v.ref_allele.as_bytes(), v.alt_allele.as_bytes());
    let del_snv = a.len() == 1 && r.len() > 1 && !a[0].eq_ignore_ascii_case(&r[0]);
    if !(r.len() > 1 && a.len() > 1 || del_snv) {
        return None; // SNVs and pure indels have their own classifiers
    }
    let reference = upper(reference);
    let ev = event(start, &reference, v)?;
    let (need_l, need_r) = ev.room_needed();
    let short = (ev.lo < need_l, ev.hi + need_r > reference.len());
    (short.0 || short.1).then_some(short)
}

fn upper(s: &str) -> Vec<u8> {
    s.bytes().map(|b| b.to_ascii_uppercase()).collect()
}

/// The variant's event on a reference: its haplotype, where REF and ALT first and
/// last differ, and the event grown through repeats (REF offsets).
struct Event {
    hap: Vec<u8>,
    d: i64,
    /// First base where REF and ALT differ, reading from the left.
    first: usize,
    /// One past the last base where they differ, reading from the right.
    last: usize,
    /// Where the alleles can differ over every placement (the trims' union)...
    core: (usize, usize),
    /// ...and that grown through repeats touching it.
    lo: usize,
    hi: usize,
}

impl Event {
    /// Whether the REF or ALT window (event plus flank) exceeds [`LONG_EVENT`].
    fn long(&self) -> bool {
        let ref_len = self.hi - self.lo + 2 * FLANK;
        (ref_len as i64 + self.d.max(0)) as usize > LONG_EVENT
    }

    /// Reference bases needed left and right of the event: the flank, plus the
    /// shorter window's padding for whole windows.
    fn room_needed(&self) -> (usize, usize) {
        if self.long() {
            return (FLANK, FLANK);
        }
        let pad = self.d.unsigned_abs() as usize;
        (FLANK + pad / 2, FLANK + pad - pad / 2)
    }
}

/// The event on `reference` (starting at `start`). None with an empty allele or a
/// REF that does not match the reference there.
fn event(start: i64, reference: &[u8], v: &Variant) -> Option<Event> {
    let (r_al, a_al) = (upper(&v.ref_allele), upper(&v.alt_allele));
    if r_al.is_empty() || a_al.is_empty() {
        return None;
    }
    let p = usize::try_from(v.pos - start).ok()?;
    if p + r_al.len() > reference.len() || reference[p..p + r_al.len()] != r_al[..] {
        return None;
    }
    let mut hap = reference[..p].to_vec();
    hap.extend_from_slice(&a_al);
    hap.extend_from_slice(&reference[p + r_al.len()..]);
    let d = hap.len() as i64 - reference.len() as i64;

    let (l_lo, l_hi) = trimmed(reference, &hap, true);
    let (r_lo, r_hi) = trimmed(reference, &hap, false);
    let (b_lo, b_hi) = (l_lo.min(r_lo).min(p), l_hi.max(r_hi).max(p + r_al.len()));
    // Grow through repeats on either allele. Left of the event the alleles share
    // coordinates; right of it the haplotype's are shifted by d.
    let unit = MAX_UNIT.max(d.unsigned_abs() as usize);
    let lo = repeat_start(reference, b_lo, unit).min(repeat_start(&hap, b_lo, unit));
    let hap_hi = repeat_end(&hap, (b_hi as i64 + d) as usize, unit) as i64 - d;
    let hi = repeat_end(reference, b_hi, unit).max(hap_hi as usize);
    Some(Event { hap, d, first: l_lo, last: r_hi, core: (b_lo, b_hi), lo, hi })
}

/// REF and ALT windows for the variant. None without a reference that holds the
/// event and its flank, or with an empty allele.
fn windows(v: &Variant) -> Option<Windows> {
    let (start, reference) = match (&v.event_ref, &v.ref_context) {
        (Some((s, seq)), _) => (*s, upper(seq)),
        (None, Some(ctx)) => (v.ref_context_start, upper(ctx)),
        (None, None) => return None,
    };
    let ev = event(start, &reference, v)?;
    if ev.lo < FLANK || ev.hi + FLANK > reference.len() {
        return None; // the reference does not hold the event's flank
    }
    let (hap, d) = (&ev.hap, ev.d);
    let lo = ev.lo - FLANK;
    let hi = ev.hi + FLANK;
    let alt_hi = (hi as i64 + d) as usize;
    let (ref_win, alt_win) = (&reference[lo..hi], &hap[lo..alt_hi]);
    let g = |off: usize| start + off as i64;

    if ev.long() {
        // Equal-length junction windows, reading inward from each end: through one
        // base past the first difference, and through the shorter allele when it
        // fits, so a short ALT is read base by base. The left reads it from the
        // left flank; the right from one base before its first difference, not
        // back through growth on the left, which the shorter allele's reads
        // starting inside it could not hold while the longer allele's could.
        let short = ref_win.len().min(alt_win.len());
        let fits = short <= LONG_EVENT;
        let short_end = (hi as i64 + d.min(0)) as usize;
        let reach = |to_diff: usize, whole: usize| (to_diff + 2).max(FLANK + 2).max(whole).min(short);
        let j_left = reach(ev.first - lo, if fits { short } else { 0 });
        let j_right = reach(hi - ev.last, if fits { short_end + 1 - ev.first } else { 0 });
        // `end` is the allele's event end in its own offsets (REF ev.hi, ALT ev.hi + d).
        let (c_lo, c_hi) = ev.core;
        let left = |seq: &[u8], end: usize| Window {
            seq: seq[..j_left].to_vec(),
            left: Some(g(lo)),
            right: None,
            before: c_lo - lo,
            after: (lo + j_left).saturating_sub(end),
        };
        let right = |seq: &[u8]| Window {
            seq: seq[seq.len() - j_right..].to_vec(),
            left: None,
            right: Some(g(hi)),
            before: c_lo.saturating_sub(lo + seq.len() - j_right),
            after: hi - c_hi,
        };
        let alt_end = (c_hi as i64 + d) as usize;
        return Some(Windows {
            pairs: vec![(left(ref_win, c_hi), left(alt_win, alt_end)), (right(ref_win), right(alt_win))],
            event: (g(c_lo), g(c_hi)),
        });
    }

    // Pad the shorter window with reference flank so both have equal length; at a
    // contig end the padding moves to the other side.
    let (short_is_ref, pad) = if ref_win.len() < alt_win.len() {
        (true, alt_win.len() - ref_win.len())
    } else {
        (false, ref_win.len() - alt_win.len())
    };
    let (mut pl, mut pr) = (pad / 2, pad - pad / 2);
    if pl > lo {
        pr += pl - lo;
        pl = lo;
    }
    if pr > reference.len() - hi {
        pl = (pl + pr - (reference.len() - hi)).min(lo);
        pr = reference.len() - hi;
    }
    let (ref_w, alt_w, span) = if short_is_ref {
        (reference[lo - pl..hi + pr].to_vec(), alt_win.to_vec(), [(lo - pl, hi + pr), (lo, hi)])
    } else {
        let alt_hi_ext = (alt_hi + pr).min(hap.len());
        (ref_win.to_vec(), hap[lo - pl..alt_hi_ext].to_vec(), [(lo, hi), (lo - pl, hi + pr)])
    };
    let (c_lo, c_hi) = ev.core;
    let whole = |seq: Vec<u8>, (a, b): (usize, usize)| Window {
        seq,
        left: Some(g(a)),
        right: Some(g(b)),
        before: c_lo - a,
        after: b - c_hi,
    };
    Some(Windows { pairs: vec![(whole(ref_w, span[0]), whole(alt_w, span[1]))], event: (g(c_lo), g(c_hi)) })
}

/// The minimal differing interval of `a` (REF offsets) against `b`, trimming the
/// shared prefix first (`prefix_first`) or the shared suffix first.
fn trimmed(a: &[u8], b: &[u8], prefix_first: bool) -> (usize, usize) {
    let limit = a.len().min(b.len());
    let prefix = |cap: usize| a.iter().zip(b).take(cap).take_while(|(x, y)| x == y).count();
    let suffix = |cap: usize| a.iter().rev().zip(b.iter().rev()).take(cap).take_while(|(x, y)| x == y).count();
    let (pre, suf) = if prefix_first {
        let pre = prefix(limit);
        (pre, suffix(limit - pre))
    } else {
        let suf = suffix(limit);
        (prefix(limit - suf), suf)
    };
    (pre, a.len() - suf)
}

/// Start of the longest tandem repeat (unit 1..=`max_unit`, at least one whole
/// unit) that runs left from `b`, continuing the event's own bases.
fn repeat_start(r: &[u8], b: usize, max_unit: usize) -> usize {
    let mut best = b;
    for k in 1..=max_unit {
        let mut lo = b;
        while lo >= 1 && lo - 1 + k < r.len() && r[lo - 1] == r[lo - 1 + k] {
            lo -= 1;
        }
        if b - lo >= k {
            best = best.min(lo);
        }
    }
    best
}

/// End (exclusive) of the longest tandem repeat that runs right from `b`.
fn repeat_end(r: &[u8], b: usize, max_unit: usize) -> usize {
    let mut best = b;
    for k in 1..=max_unit {
        let mut hi = b;
        while hi < r.len() && hi >= k && r[hi] == r[hi - k] {
            hi += 1;
        }
        if hi - b >= k {
            best = best.max(hi);
        }
    }
    best
}

/// The query position aligned to reference position `g` (an M/=/X base), if any.
fn qpos(record: &Record, g: i64) -> Option<usize> {
    let (mut ref_pos, mut read_pos) = (record.pos(), 0usize);
    for op in record.cigar().iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len = *len as i64;
                if g >= ref_pos && g < ref_pos + len {
                    return Some(read_pos + (g - ref_pos) as usize);
                }
                ref_pos += len;
                read_pos += len as usize;
            }
            Cigar::Ins(len) | Cigar::SoftClip(len) => read_pos += *len as usize,
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                if g >= ref_pos && g < ref_pos + *len as i64 {
                    return None;
                }
                ref_pos += *len as i64;
            }
            _ => {}
        }
    }
    None
}

/// Read `w` from the read at the window's own position: from the query position
/// of its first reference base, and back from its last. For a read aligned with
/// the event between the anchors the two readings are the same bases; the closer
/// one counts. None when the read holds the window from neither anchor.
fn read_window(record: &Record, seq: &[u8], quals: &[u8], min_baseq: u8, w: &Window) -> Option<Reading> {
    let n = w.seq.len();
    let from_left = w.left.and_then(|g| qpos(record, g)).filter(|&q| q + n <= seq.len()).map(|q| (q, q + n));
    let to_right = w
        .right
        .and_then(|g| qpos(record, g - 1))
        .map(|q| q + 1)
        .filter(|&e| e >= n && e <= seq.len())
        .map(|e| (e - n, e));
    [from_left, to_right]
        .into_iter()
        .flatten()
        .map(|(a, b)| score(&seq[a..b], &quals[a..b], min_baseq, &w.seq, (a, b)))
        .min_by_key(|r| (r.mismatches, r.masked))
}

/// Mismatches of read bases against haplotype bases position by position, masked
/// bases (below `min_baseq`, or N) matching anything.
fn score(bases: &[u8], quals: &[u8], min_baseq: u8, hap: &[u8], span: (usize, usize)) -> Reading {
    let (mut mismatches, mut masked, mut had_n) = (0, 0, false);
    for ((&b, &q), &h) in bases.iter().zip(quals).zip(hap) {
        let b = b.to_ascii_uppercase();
        had_n |= b == WILD;
        if b == WILD || q < min_baseq {
            masked += 1;
        } else if b != h {
            mismatches += 1;
        }
    }
    Reading { mismatches, masked, had_n, span }
}

/// The read's splices (CIGAR N, genomic [start, end)) that overlap `[lo, hi)`.
fn ref_skips(record: &Record, (lo, hi): (i64, i64)) -> Vec<(i64, i64)> {
    let mut skips = Vec::new();
    let mut pos = record.pos();
    for op in record.cigar().iter() {
        match op {
            Cigar::RefSkip(len) => {
                let end = pos + *len as i64;
                if *len > 0 && pos < hi && end > lo {
                    skips.push((pos, end));
                }
                pos = end;
            }
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) | Cigar::Del(len) => pos += *len as i64,
            _ => {}
        }
    }
    skips
}

/// Median of the qualities at or above `min_baseq` (0 when none).
fn median(quals: &[u8], min_baseq: u8) -> u8 {
    let mut q: Vec<u8> = quals.iter().copied().filter(|&x| x >= min_baseq).collect();
    if q.is_empty() {
        return 0;
    }
    q.sort_unstable();
    q[q.len() / 2]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn var(ctx: &str, pos: i64, r: &str, a: &str) -> Variant {
        Variant::new(
            "1".into(), pos, r.into(), a.into(), "COMPLEX".into(),
            None, 0, 0, None, None, Some((0, ctx.into())),
        )
    }

    fn only_pair(w: &Windows) -> (&Window, &Window) {
        assert_eq!(w.pairs.len(), 1);
        (&w.pairs[0].0, &w.pairs[0].1)
    }

    #[test]
    fn windows_are_the_event_with_two_flank_bases_and_equal_length() {
        // GGAC [TA] GTCGTT -> GGAC [GCC] GTCGTT: the event plus two flank bases on
        // each side; the REF window takes one more flank base to match ALT's length.
        let w = windows(&var("GGACTAGTCGTT", 4, "TA", "GCC")).unwrap();
        let (r, a) = only_pair(&w);
        assert_eq!(a.seq, b"ACGCCGT".to_vec());
        assert_eq!(r.seq, b"ACTAGTC".to_vec());
        assert_eq!((r.left, r.right), (Some(2), Some(9)));
        assert_eq!((a.left, a.right), (Some(2), Some(8)));
    }

    #[test]
    fn the_event_grows_through_runs_on_either_side() {
        // CAC>A inside CCCC A CCCC: a C deleted from each run, placeable anywhere
        // in it; the windows are anchored outside both runs.
        let w = windows(&var("TTGCCCCACCCCTTA", 6, "CAC", "A")).unwrap();
        let (r, a) = only_pair(&w);
        assert_eq!(r.left, Some(1)); // two flank bases before the first C run (3..)
        assert_eq!(r.seq, b"TGCCCCACCCCTT".to_vec());
        assert_eq!(a.seq.len(), r.seq.len());
    }

    #[test]
    fn the_event_grows_through_a_run_the_alt_continues() {
        // GTC AAA [TA] CGT -> GTC AAA [ACC] CGT: the ALT's A continues the run, so
        // an inserted A can sit anywhere in it; the left anchor moves before it.
        let w = windows(&var("GCGTCAAATACGTGG", 8, "TA", "ACC")).unwrap();
        let (r, _) = only_pair(&w);
        assert_eq!(r.left, Some(3)); // two flank bases before the run (5..)
    }

    #[test]
    fn repeat_growth_needs_a_whole_unit() {
        assert_eq!(repeat_start(b"GACTACAGG", 4, MAX_UNIT), 4);
        assert_eq!(repeat_start(b"GCCCCA", 4, MAX_UNIT), 1);
        assert_eq!(repeat_end(b"ACCCCG", 2, MAX_UNIT), 5);
    }

    #[test]
    fn trimming_from_either_side() {
        assert_eq!(trimmed(b"GAAAAC", b"GAAAC", true), (4, 5));
        assert_eq!(trimmed(b"GAAAAC", b"GAAAC", false), (1, 2));
    }

    #[test]
    fn long_events_read_the_short_allele_at_both_junctions() {
        let mut ctx = String::from("ACGTTGCA");
        let body: String = (0..70).map(|i| b"ACGGTTCA"[i % 8] as char).collect();
        ctx.push_str(&body);
        ctx.push_str("TGCATTGC");
        let w = windows(&var(&ctx, 7, &ctx[7..76], "TG")).unwrap();
        assert_eq!(w.pairs.len(), 2);
        let alt_whole = &w.pairs[0].1.seq; // the left pair holds the whole ALT window
        for (r, a) in &w.pairs {
            assert_eq!(r.seq.len(), a.seq.len());
            assert_ne!(r.seq, a.seq);
        }
        // The right pair reads the ALT from one base before its first difference.
        let right_alt = &w.pairs[1].1.seq;
        assert!(alt_whole.ends_with(right_alt));
        assert!(right_alt.windows(2).any(|x| x == b"TG"));
    }

    #[test]
    fn a_short_reference_is_reported_on_its_side() {
        // The event ends a run reaching the reference's left end: prep must fetch
        // more on the left only.
        let v = var("AAAAAAAAAATTCGTGGTCAGCTTGACCA", 9, "ATT", "GC");
        assert_eq!(reference_short(0, "AAAAAAAAAATTCGTGGTCAGCTTGACCA", &v), Some((true, false)));
        let v = var("GCGTCGAAAAAAAAAATTCGTGGTCAGCTTGACCA", 15, "ATT", "GC");
        assert_eq!(reference_short(0, "GCGTCGAAAAAAAAAATTCGTGGTCAGCTTGACCA", &v), None);
    }
}
