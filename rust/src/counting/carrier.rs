//! Exact-carrier classification for complex variants.
//!
//! A complex variant (a delins, a deletion whose anchor also changes, or an MNP
//! read carrying an indel) is REF or ALT for a read only when the read's own
//! bases carry that allele across the whole event, with [`FLANK`] reference bases
//! on each side. The read's bases include soft clips; bases below min BQ are
//! wildcards (the one quality rule every backend shares). Nothing else is
//! tolerated: a read one confident base off the given allele carries a
//! different allele, and a read that ends inside the event cannot show either.
//!
//! - **The event** is every base where REF and ALT can differ over all equivalent
//!   placements: the union of the minimal differing interval trimmed from the
//!   left first and from the right first, so a delins in a repeat is judged over
//!   its whole ambiguity.
//! - **Equal windows.** The shorter allele's window is padded with reference
//!   flank until both are the same length, so neither allele is favoured by where
//!   reads start.
//! - **Long events.** When the windows exceed [`LONG_EVENT`] bases no read can
//!   hold them whole. Both alleles are then judged by equal-length junction
//!   windows at each end: a flank through one base past the first base where
//!   they differ, reading inward.
//! - **Outcomes.** A read that carries both windows (only through masked bases)
//!   is neither. A read carrying neither is partial evidence when it is closer to
//!   ALT than to REF.

use rust_htslib::bam::record::{Cigar, Record};

use super::utils::{ClassifyPhase, ClassifyResult};
use crate::types::Variant;

/// Reference bases required on each side of the event.
const FLANK: usize = 2;
/// Window length beyond which the whole allele cannot be expected in one read.
const LONG_EVENT: usize = 50;
/// Masked base: matches any haplotype base.
const WILD: u8 = b'N';

/// The REF and ALT windows a read is tested against, as haplotype strings, and
/// the genomic stretch they cover (for gathering the read's bases).
struct Windows {
    reference: Vec<Vec<u8>>,
    alternate: Vec<Vec<u8>>,
    genomic: (i64, i64),
}

/// Whether the exact-carrier rule can judge this variant: both alleles are
/// non-empty and a reference (prep's `event_ref`, else the `ref_context`)
/// covers the event. When it cannot (a variant built without prep, or a
/// failed reference fetch), the caller uses the previous classifier, whose
/// `SW_FALLBACK` flag keeps that malformed-input case visible.
pub(crate) fn can_judge(variant: &Variant) -> bool {
    (variant.event_ref.is_some() || variant.ref_context.is_some()) && windows(variant).is_some()
}

/// Classify a read at a complex variant by the exact-carrier rule. Call only
/// when [`can_judge`] holds.
pub(crate) fn check_complex_exact(
    record: &Record,
    variant: &Variant,
    quals: &[u8],
    min_baseq: u8,
) -> ClassifyResult {
    let Some(win) = windows(variant) else {
        return ClassifyResult::neither(ClassifyPhase::MaskedCompare);
    };
    if splices_over(record, win.genomic) {
        return ClassifyResult::neither(ClassifyPhase::CigarRecon);
    }
    let longest = win.reference.iter().chain(&win.alternate).map(Vec::len).max().unwrap_or(0);
    let (seq, qs, had_n) = local_bases(record, quals, min_baseq, win.genomic, longest as i64 + 10);
    if seq.is_empty() {
        return ClassifyResult::neither(ClassifyPhase::MaskedCompare);
    }
    let qual = median(&qs, min_baseq);
    let alt = win.alternate.iter().any(|w| contains(&seq, w));
    let reference = win.reference.iter().any(|w| contains(&seq, w));
    let mut result = match (reference, alt) {
        (false, true) => ClassifyResult::is_alt(qual, ClassifyPhase::MaskedCompare),
        (true, false) => ClassifyResult::is_ref(qual, ClassifyPhase::MaskedCompare),
        (true, true) => ClassifyResult::neither(ClassifyPhase::MaskedCompare),
        (false, false) => {
            let d_alt = win.alternate.iter().map(|w| fit_distance(w, &seq)).min().unwrap_or(usize::MAX);
            let d_ref = win.reference.iter().map(|w| fit_distance(w, &seq)).min().unwrap_or(usize::MAX);
            if d_alt < d_ref {
                ClassifyResult::neither_with_nearby(qual, ClassifyPhase::MaskedCompare)
            } else {
                ClassifyResult::neither(ClassifyPhase::MaskedCompare)
            }
        }
    };
    result.has_n_base = had_n && !result.is_ref && !result.is_alt;
    result
}

/// REF and ALT windows for the variant, from prep's `event_ref` (else the
/// `ref_context`). None without a reference that holds the event and its flank,
/// or with an empty allele.
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

    // The event over every equivalent placement, widened to the given alleles.
    let (l_lo, l_hi) = trimmed(&reference, &hap, true);
    let (r_lo, r_hi) = trimmed(&reference, &hap, false);
    let e_lo = l_lo.min(r_lo).min(p);
    let e_hi = l_hi.max(r_hi).max(p + r_al.len());
    if e_lo < FLANK || e_hi + FLANK > reference.len() {
        return None; // the reference does not hold the event's flank
    }
    let lo = e_lo - FLANK;
    let hi = e_hi + FLANK;
    let alt_hi = (hi as i64 + d) as usize;
    let (ref_win, alt_win) = (&reference[lo..hi], &hap[lo..alt_hi]);

    if ref_win.len().max(alt_win.len()) > LONG_EVENT {
        // Junction windows of equal length on both alleles.
        let j_left = (l_lo + 2).min(ref_win.len() + lo).saturating_sub(lo).min(alt_win.len());
        let j_right = (hi - r_hi.min(hi) + 2).min(ref_win.len()).min(alt_win.len());
        return Some(Windows {
            reference: vec![ref_win[..j_left].to_vec(), ref_win[ref_win.len() - j_right..].to_vec()],
            alternate: vec![alt_win[..j_left].to_vec(), alt_win[alt_win.len() - j_right..].to_vec()],
            genomic: (start + lo as i64, start + hi as i64),
        });
    }

    // Pad the shorter window with reference flank so both have equal length.
    let (short_is_ref, pad) = if ref_win.len() < alt_win.len() {
        (true, alt_win.len() - ref_win.len())
    } else {
        (false, ref_win.len() - alt_win.len())
    };
    let (mut pl, mut pr) = (pad / 2, pad - pad / 2);
    let left_room = lo;
    let right_room = reference.len() - hi;
    if pl > left_room {
        pr += pl - left_room;
        pl = left_room;
    }
    if pr > right_room {
        pl = (pl + pr - right_room).min(left_room);
        pr = right_room;
    }
    let (ref_w, alt_w) = if short_is_ref {
        (reference[lo - pl..hi + pr].to_vec(), alt_win.to_vec())
    } else {
        let alt_lo = lo - pl;
        let alt_hi_ext = (alt_hi + pr).min(hap.len());
        (ref_win.to_vec(), hap[alt_lo..alt_hi_ext].to_vec())
    };
    let g_lo = start + (lo - pl) as i64;
    let g_hi = start + (hi + pr) as i64;
    Some(Windows { reference: vec![ref_w], alternate: vec![alt_w], genomic: (g_lo, g_hi) })
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

/// The read's bases around `[lo, hi)` extended by `margin` on each side: aligned
/// bases there, insertions among them, and soft clips whose boundary lies there.
/// Bases below `min_baseq` (and N) become [`WILD`]. Returns the bases, their
/// qualities, and whether an N base was seen.
fn local_bases(
    record: &Record,
    quals: &[u8],
    min_baseq: u8,
    (lo, hi): (i64, i64),
    margin: i64,
) -> (Vec<u8>, Vec<u8>, bool) {
    let (g_lo, g_hi) = (lo - margin, hi + margin);
    let seq = record.seq().as_bytes();
    let (mut q_lo, mut q_hi): (Option<usize>, Option<usize>) = (None, None);
    let (mut ref_pos, mut read_pos) = (record.pos(), 0usize);
    let cigar = record.cigar();
    let n_ops = cigar.len();
    for (i, op) in cigar.iter().enumerate() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len = *len as i64;
                let a = ref_pos.max(g_lo);
                let b = (ref_pos + len).min(g_hi);
                if a < b {
                    let qa = read_pos + (a - ref_pos) as usize;
                    let qb = read_pos + (b - ref_pos) as usize;
                    q_lo = Some(q_lo.map_or(qa, |q| q.min(qa)));
                    q_hi = Some(q_hi.map_or(qb, |q| q.max(qb)));
                }
                ref_pos += len;
                read_pos += len as usize;
            }
            Cigar::Ins(len) => read_pos += *len as usize,
            Cigar::Del(len) | Cigar::RefSkip(len) => ref_pos += *len as i64,
            Cigar::SoftClip(len) => {
                let len = *len as usize;
                if (g_lo..=g_hi).contains(&ref_pos) {
                    let leading = i == 0 || (i == 1 && matches!(cigar.first(), Some(Cigar::HardClip(_))));
                    if leading {
                        q_lo = Some(0);
                        q_hi = Some(q_hi.map_or(len, |q| q.max(len)));
                    } else if i + 1 == n_ops || matches!(cigar.get(i + 1), Some(Cigar::HardClip(_))) {
                        q_lo = Some(q_lo.map_or(read_pos, |q| q.min(read_pos)));
                        q_hi = Some(read_pos + len);
                    }
                }
                read_pos += len;
            }
            _ => {}
        }
    }
    let (Some(a), Some(b)) = (q_lo, q_hi) else {
        return (Vec::new(), Vec::new(), false);
    };
    let b = b.min(seq.len());
    let mut had_n = false;
    let bases = (a..b)
        .map(|i| {
            let base = seq[i].to_ascii_uppercase();
            if base == b'N' {
                had_n = true;
            }
            if base == b'N' || quals.get(i).copied().unwrap_or(0) < min_baseq { WILD } else { base }
        })
        .collect();
    (bases, quals[a..b].to_vec(), had_n)
}

/// Whether `text` contains `pat` exactly, [`WILD`] bases of `text` matching anything.
fn contains(text: &[u8], pat: &[u8]) -> bool {
    !pat.is_empty()
        && text.len() >= pat.len()
        && text.windows(pat.len()).any(|w| w.iter().zip(pat).all(|(t, p)| t == p || *t == WILD))
}

/// Edit distance of `pat` against its best-matching stretch of `text` (fitting
/// alignment), [`WILD`] bases of `text` matching anything.
fn fit_distance(pat: &[u8], text: &[u8]) -> usize {
    let mut prev: Vec<usize> = vec![0; text.len() + 1];
    for (i, &p) in pat.iter().enumerate() {
        let mut cur = vec![i + 1; text.len() + 1];
        for (j, &t) in text.iter().enumerate() {
            let sub = prev[j] + usize::from(!(t == p || t == WILD));
            cur[j + 1] = sub.min(prev[j + 1] + 1).min(cur[j] + 1);
        }
        prev = cur;
    }
    prev.into_iter().min().unwrap_or(pat.len())
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

    #[test]
    fn windows_are_the_event_with_two_flank_bases_and_equal_length() {
        // GGAC [TA] CAGGTT -> GGAC [GCC] CAGGTT: the event plus two flank bases on
        // each side; the REF window takes one more flank base to match ALT's length.
        let w = windows(&var("GGACTACAGGTT", 4, "TA", "GCC")).unwrap();
        assert_eq!(w.alternate, vec![b"ACGCCCA".to_vec()]);
        assert_eq!(w.reference, vec![b"ACTACAG".to_vec()]);
    }

    #[test]
    fn an_event_in_a_repeat_is_judged_over_its_whole_ambiguity() {
        // TTTT -> TTA inside a T run: the differing interval moves with trimming
        // order; the window must cover both placements.
        let w = windows(&var("GCATTTTTGC", 3, "TTTT", "TTA")).unwrap();
        for h in w.reference.iter().chain(&w.alternate) {
            assert!(h.starts_with(b"CA"), "{:?}", String::from_utf8_lossy(h));
        }
    }

    #[test]
    fn containment_masks_only_wild_bases() {
        assert!(contains(b"AAGTNCAA", b"GTGCA"));
        assert!(!contains(b"AAGTACAA", b"GTGCA"));
    }

    #[test]
    fn fitting_distance_is_free_at_the_text_ends() {
        assert_eq!(fit_distance(b"GCC", b"AAGCCTT"), 0);
        assert_eq!(fit_distance(b"GCC", b"AAGCATT"), 1);
        assert_eq!(fit_distance(b"GCCA", b"TTGC"), 2);
    }

    #[test]
    fn trimming_from_either_side() {
        // REF AAAA vs ALT AAA: the differing base floats across the run.
        assert_eq!(trimmed(b"GAAAAC", b"GAAAC", true), (4, 5));
        assert_eq!(trimmed(b"GAAAAC", b"GAAAC", false), (1, 2));
    }
}
