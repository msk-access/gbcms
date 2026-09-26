//! Exact-carrier classification for complex variants.
//!
//! A complex variant (a delins, a deletion whose anchor also changes, or an MNP
//! read carrying an indel) is REF or ALT for a read only when the read's own
//! bases carry that allele across the whole event, with [`FLANK`] reference bases
//! on each side. The read's bases include soft clips; bases below min BQ (and N)
//! match anything, the one quality rule every backend shares. Nothing else is
//! tolerated: a read one confident base off the given allele carries a different
//! allele, and a read that ends inside the event cannot show either.
//!
//! - **The event** is every base where REF and ALT can differ over all equivalent
//!   placements. That is the union of the minimal difference trimmed from the
//!   left first and from the right first, grown through any tandem repeat (unit
//!   1-6) touching either end, since an aligner may place an indel anywhere in
//!   such a repeat.
//! - **Windows are read where they sit.** Each window is compared with the read's
//!   bases starting at the query position of the window's first reference base
//!   (left anchor), or ending at its last (right anchor). The anchors are flank
//!   bases outside the event, aligned the same way on either allele; within the
//!   event the read's bases are compared as a contiguous stretch, so where the
//!   aligner put an indel (or a clip) does not matter. A copy of a window
//!   elsewhere in the read is never a match.
//! - **Equal windows.** The shorter allele's window is padded with reference
//!   flank until both are the same length, so neither allele is favoured by where
//!   reads start.
//! - **Long events.** When the windows exceed [`LONG_EVENT`] bases no read can
//!   hold them whole. Both alleles are then judged by equal-length junction
//!   windows at each end: a flank through one base past the first base where they
//!   differ, reading inward. A read covering both ends must agree at both.
//! - **Outcomes.** A read matching both windows (through masked bases) is neither.
//!   A read that cannot complete both windows of a pair is uninformative for it:
//!   depth only, like a pure-indel read ending inside its tract. A read that completes the windows
//!   and matches neither carries another allele. Either is partial evidence when
//!   closer to ALT than to REF.

use rust_htslib::bam::record::{Cigar, Record};

use super::utils::{ClassifyPhase, ClassifyResult};
use crate::types::Variant;

/// Reference bases required on each side of the event.
const FLANK: usize = 2;
/// Window length beyond which the whole allele cannot be expected in one read.
const LONG_EVENT: usize = 50;
/// Longest tandem-repeat unit the event is grown through.
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
}

/// The REF and ALT windows, paired: whole-event windows (one pair) or junction
/// windows (a left pair and a right pair).
struct Windows {
    pairs: Vec<(Window, Window)>,
}

/// How a window reads in the read.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Read {
    Match,
    Mismatch,
    /// The read does not reach across the window from any anchor.
    Incomplete,
}

/// A window compared with the read: outcome, mismatching bases over the part
/// the read covers, and whether an N sat in the compared bases.
struct Compared {
    read: Read,
    mismatches: usize,
    had_n: bool,
}

/// Classify a read at a complex variant by the exact-carrier rule, or None when
/// the rule cannot judge this variant: it needs both alleles non-empty and a
/// reference (prep's `event_ref`, else the `ref_context`) holding the event and
/// its flank. The caller then uses the previous classifier, whose `SW_FALLBACK`
/// flag keeps that malformed-input case visible.
pub(crate) fn classify(record: &Record, variant: &Variant, quals: &[u8], min_baseq: u8) -> Option<ClassifyResult> {
    let win = windows(variant)?;
    let (lo, hi) = win.pairs.iter().fold((i64::MAX, i64::MIN), |(lo, hi), (r, a)| {
        let starts = [r.left, a.left].into_iter().flatten();
        let ends = [r.right, a.right].into_iter().flatten();
        (starts.fold(lo, i64::min), ends.fold(hi, i64::max))
    });
    if splices_over(record, (lo, hi)) {
        return Some(ClassifyResult::neither(ClassifyPhase::CigarRecon));
    }

    let mut calls: Vec<Option<bool>> = Vec::new(); // per pair: Some(true)=ALT, Some(false)=REF, None=both
    let (mut complete, mut incomplete, mut had_n) = (false, false, false);
    let (mut mm_ref, mut mm_alt) = (0usize, 0usize);
    let mut qual_bases: Vec<u8> = Vec::new();
    for (rw, aw) in &win.pairs {
        let r = compare(record, quals, min_baseq, rw, &mut qual_bases);
        let a = compare(record, quals, min_baseq, aw, &mut qual_bases);
        had_n |= r.had_n || a.had_n;
        mm_ref += r.mismatches;
        mm_alt += a.mismatches;
        // A pair decides only when the read holds both windows: padding can let one
        // allele's window complete where the other's cannot (a read ending in the
        // padding), and a read that cannot show both alleles cannot choose.
        match (r.read, a.read) {
            (Read::Incomplete, _) | (_, Read::Incomplete) => incomplete = true,
            (Read::Match, Read::Match) => {
                complete = true;
                calls.push(None);
            }
            (Read::Mismatch, Read::Match) => {
                complete = true;
                calls.push(Some(true));
            }
            (Read::Match, Read::Mismatch) => {
                complete = true;
                calls.push(Some(false));
            }
            (Read::Mismatch, Read::Mismatch) => complete = true,
        }
    }
    let qual = median(&qual_bases, min_baseq);
    let decided: Vec<bool> = calls.iter().flatten().copied().collect();
    let agreed = !decided.is_empty() && !calls.contains(&None) && decided.iter().all(|&x| x == decided[0]);
    let mut result = if agreed && decided[0] {
        ClassifyResult::is_alt(qual, ClassifyPhase::MaskedCompare)
    } else if agreed {
        ClassifyResult::is_ref(qual, ClassifyPhase::MaskedCompare)
    } else if mm_alt < mm_ref {
        // Another allele, or too little of the read: closer to ALT is partial evidence.
        ClassifyResult::neither_with_nearby(qual, ClassifyPhase::MaskedCompare)
    } else {
        ClassifyResult::neither(ClassifyPhase::MaskedCompare)
    };
    if !result.is_ref && !result.is_alt && incomplete && !complete {
        // The read cannot show the allele: depth only (as a pure-indel read that
        // ends inside its tract), and no mFSD class.
        result.ref_uninformative = true;
    }
    result.has_n_base = had_n;
    Some(result)
}

/// REF and ALT windows for the variant. None without a reference that holds the
/// event and its flank, or with an empty allele.
fn windows(v: &Variant) -> Option<Windows> {
    let upper = |s: &str| -> Vec<u8> { s.bytes().map(|b| b.to_ascii_uppercase()).collect() };
    let (start, reference) = match (&v.event_ref, &v.ref_context) {
        (Some((s, seq)), _) => (*s, upper(seq)),
        (None, Some(ctx)) => (v.ref_context_start, upper(ctx)),
        (None, None) => return None,
    };
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

    // The event over every equivalent placement, widened to the given alleles and
    // grown through repeats touching either end.
    let (l_lo, l_hi) = trimmed(&reference, &hap, true);
    let (r_lo, r_hi) = trimmed(&reference, &hap, false);
    let e_lo = repeat_start(&reference, l_lo.min(r_lo).min(p));
    let e_hi = repeat_end(&reference, l_hi.max(r_hi).max(p + r_al.len()));
    if e_lo < FLANK || e_hi + FLANK > reference.len() {
        return None; // the reference does not hold the event's flank
    }
    let lo = e_lo - FLANK;
    let hi = e_hi + FLANK;
    let alt_hi = (hi as i64 + d) as usize;
    let (ref_win, alt_win) = (&reference[lo..hi], &hap[lo..alt_hi]);
    let g = |off: usize| start + off as i64;

    if ref_win.len().max(alt_win.len()) > LONG_EVENT {
        // Equal-length junction windows: flank through one base past the first
        // difference, reading inward from each end.
        let j_left = (l_lo.min(r_lo) + 2).saturating_sub(lo).max(FLANK + 2).min(ref_win.len()).min(alt_win.len());
        let j_right = (hi + 2).saturating_sub(l_hi.max(r_hi)).max(FLANK + 2).min(ref_win.len()).min(alt_win.len());
        let left = |seq: &[u8]| Window { seq: seq[..j_left].to_vec(), left: Some(g(lo)), right: None };
        let right = |seq: &[u8]| Window { seq: seq[seq.len() - j_right..].to_vec(), left: None, right: Some(g(hi)) };
        return Some(Windows { pairs: vec![(left(ref_win), left(alt_win)), (right(ref_win), right(alt_win))] });
    }

    // Pad the shorter window with reference flank so both have equal length.
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
    let whole = |seq: Vec<u8>, (a, b): (usize, usize)| Window { seq, left: Some(g(a)), right: Some(g(b)) };
    Some(Windows { pairs: vec![(whole(ref_w, span[0]), whole(alt_w, span[1]))] })
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

/// Start of the longest tandem repeat (unit 1..=MAX_UNIT, at least one whole
/// unit) that runs left from `b`, continuing the event's own bases.
fn repeat_start(r: &[u8], b: usize) -> usize {
    let mut best = b;
    for k in 1..=MAX_UNIT {
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
fn repeat_end(r: &[u8], b: usize) -> usize {
    let mut best = b;
    for k in 1..=MAX_UNIT {
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

/// Compare `w` with the read's bases at the window's own position: from the
/// query position of its first reference base, else back from its last.
fn compare(record: &Record, quals: &[u8], min_baseq: u8, w: &Window, qual_bases: &mut Vec<u8>) -> Compared {
    let seq = record.seq().as_bytes();
    let n = w.seq.len();
    let from_left = w.left.and_then(|g| qpos(record, g));
    let to_right = w.right.and_then(|g| qpos(record, g - 1)).map(|q| q + 1);
    let whole = from_left
        .filter(|&q| q + n <= seq.len())
        .map(|q| (q, q + n))
        .or_else(|| to_right.filter(|&e| e >= n).map(|e| (e - n, e)));
    if let Some((a, b)) = whole {
        let (mismatches, had_n) = score(&seq[a..b], &quals[a..b], min_baseq, &w.seq, qual_bases);
        let read = if mismatches == 0 { Read::Match } else { Read::Mismatch };
        return Compared { read, mismatches, had_n };
    }
    // The read ends inside the window: score what it covers from an anchor.
    let (mismatches, had_n) = if let Some(q) = from_left {
        let b = seq.len().min(q + n);
        score(&seq[q..b], &quals[q..b], min_baseq, &w.seq[..b - q], qual_bases)
    } else if let Some(e) = to_right {
        let covered = e.min(n);
        score(&seq[e - covered..e], &quals[e - covered..e], min_baseq, &w.seq[n - covered..], qual_bases)
    } else {
        (0, false)
    };
    Compared { read: Read::Incomplete, mismatches, had_n }
}

/// Mismatches of read bases against haplotype bases position by position, masked
/// bases matching anything; whether an N was among them. Collects the qualities
/// of the compared bases for the fragment-consensus quality.
fn score(bases: &[u8], quals: &[u8], min_baseq: u8, hap: &[u8], qual_bases: &mut Vec<u8>) -> (usize, bool) {
    let mut had_n = false;
    let mut mismatches = 0;
    for ((&b, &q), &h) in bases.iter().zip(quals).zip(hap) {
        let b = b.to_ascii_uppercase();
        had_n |= b == WILD;
        qual_bases.push(q);
        if b != WILD && q >= min_baseq && b != h {
            mismatches += 1;
        }
    }
    (mismatches, had_n)
}

/// Whether a splice N of the read overlaps `[lo, hi)`.
fn splices_over(record: &Record, (lo, hi): (i64, i64)) -> bool {
    let mut pos = record.pos();
    for op in record.cigar().iter() {
        match op {
            Cigar::RefSkip(len) => {
                let end = pos + *len as i64;
                if *len > 0 && pos < hi && end > lo {
                    return true;
                }
                pos = end;
            }
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) | Cigar::Del(len) => pos += *len as i64,
            _ => {}
        }
    }
    false
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
        // GGAC [TA] CAGGTT -> GGAC [GCC] CAGGTT: the event plus two flank bases on
        // each side; the REF window takes one more flank base to match ALT's length.
        let w = windows(&var("GGACTACAGGTT", 4, "TA", "GCC")).unwrap();
        let (r, a) = only_pair(&w);
        assert_eq!(a.seq, b"ACGCCCA".to_vec());
        assert_eq!(r.seq, b"ACTACAG".to_vec());
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
    fn repeat_growth_needs_a_whole_unit() {
        assert_eq!(repeat_start(b"GACTACAGG", 4), 4);
        assert_eq!(repeat_start(b"GCCCCA", 4), 1);
        assert_eq!(repeat_end(b"ACCCCG", 2), 5);
    }

    #[test]
    fn trimming_from_either_side() {
        assert_eq!(trimmed(b"GAAAAC", b"GAAAC", true), (4, 5));
        assert_eq!(trimmed(b"GAAAAC", b"GAAAC", false), (1, 2));
    }

    #[test]
    fn long_events_get_equal_junction_windows() {
        let mut ctx = String::from("ACGTTGCA");
        let body: String = (0..70).map(|i| b"ACGGTTCA"[i % 8] as char).collect();
        ctx.push_str(&body);
        ctx.push_str("TGCATTGC");
        let w = windows(&var(&ctx, 7, &ctx[7..76], "AC")).unwrap();
        assert_eq!(w.pairs.len(), 2);
        for (r, a) in &w.pairs {
            assert_eq!(r.seq.len(), a.seq.len());
            assert!(r.seq.len() >= FLANK + 2);
        }
    }
}
