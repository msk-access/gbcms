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
//! Three uses:
//! - A REF call on a pure indel stands only when the read is informative:
//!   it covers a flank of the region and runs one base past the first base
//!   where REF and ALT differ reading inward from that flank
//!   ([`read_is_informative`]). Otherwise the read counts toward depth only.
//! - An ALT call stands when the read spans the informative windows read on the
//!   ALT haplotype ([`alt_read_is_informative`]), or else when its own bases tell
//!   the alleles apart ([`alt_bases_discriminate`]): its gap alone is placement.
//! - A co-annotated sibling's carriers are excluded from a row's REF only when
//!   the sibling's change lies inside that row's discrimination window (the
//!   region plus one base on each side, [`discrimination_window`]). A carrier
//!   whose change lies elsewhere shows the reference across every base that
//!   could distinguish this row's alleles, so it is REF here.

use rust_htslib::bam::record::{Cigar, Record};

use crate::normalize::repeat::first_change_offset;
use crate::shared::bam_utils::{find_read_pos, ref_end};
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

/// What a row's alleles are, by their lengths and whether the ALT keeps the REF's
/// first (anchor) base, compared case-insensitively: one classification for the
/// counting dispatcher, splice triage, the exact-carrier rule's scope and prep's
/// variant-type label.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AlleleKind {
    /// One base each.
    Snv,
    /// Equal lengths, more than one base.
    Mnp,
    /// A one-base REF whose ALT keeps it and adds bases.
    Insertion,
    /// A one-base ALT that keeps the REF's first base.
    Deletion,
    /// Anything else: a delins, or an insertion or deletion whose anchor changes.
    Complex,
}

/// The row's allele kind; None when either allele is empty.
pub(crate) fn allele_kind(ref_allele: &str, alt_allele: &str) -> Option<AlleleKind> {
    let (r, a) = (ref_allele.as_bytes(), alt_allele.as_bytes());
    let anchor_kept = r.first()?.eq_ignore_ascii_case(a.first()?);
    Some(match (r.len(), a.len()) {
        (1, 1) => AlleleKind::Snv,
        (m, n) if m == n => AlleleKind::Mnp,
        (1, _) if anchor_kept => AlleleKind::Insertion,
        (_, 1) if anchor_kept => AlleleKind::Deletion,
        _ => AlleleKind::Complex,
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

/// How far past its anchor a variant's reads and windowed scans reach: 5 bases,
/// or its repeat span plus 2 when wider, so a read reaches past the tract.
pub(crate) fn scan_pad(variant: &Variant) -> i64 {
    pad_for_repeat_span(variant.repeat_span as i64)
}

/// `scan_pad` for a repeat span (a bin pads by its widest member's).
pub(crate) fn pad_for_repeat_span(repeat_span: i64) -> i64 {
    std::cmp::max(5, repeat_span + 2)
}

/// The variant's read window `[start, end)`: its REF span padded by `scan_pad`
/// on each side. The cached reads overlapping it are the variant's reads; each
/// bin's fetch holds it for every member (`build_genomic_bins`).
pub(crate) fn read_window(variant: &Variant) -> (i64, i64) {
    let pad = scan_pad(variant);
    ((variant.pos - pad).max(0), variant.pos + variant.ref_allele.len() as i64 + pad)
}

/// The windowed indel scan's reference range `[start, end]` (inclusive): `window`
/// bases each side of the anchor, widened to reach every placement of the variant
/// in its shift region `[lo, hi)`: every junction `lo..=hi` of an insertion, every
/// start `lo..=hi - len` of a deletion of `len` bases (a start inside the deleted
/// span past that is another deletion, not this one written elsewhere).
pub(crate) fn scan_window(variant: &Variant, window: i64) -> (i64, i64) {
    let (mut start, mut end) = (variant.pos - window, variant.pos + window);
    if let Some((lo, hi)) = variant.shift_region {
        let del_len = (variant.ref_allele.len() as i64 - variant.alt_allele.len() as i64).max(0);
        start = start.min(lo);
        end = end.max(hi - del_len);
    }
    (start.max(0), end)
}

/// Reference bases `[lo, hi)`, upper-case, from the variant's prepared reference:
/// its event reference, else its `ref_context`. None when neither holds them.
pub(crate) fn reference_span(v: &Variant, lo: i64, hi: i64) -> Option<Vec<u8>> {
    let slice = |start: i64, seq: &str| {
        let (a, b) = (lo - start, hi - start);
        (a >= 0 && a <= b && b as usize <= seq.len())
            .then(|| seq.as_bytes()[a as usize..b as usize].to_ascii_uppercase())
    };
    v.event_ref
        .as_ref()
        .and_then(|(start, seq)| slice(*start, seq))
        .or_else(|| v.ref_context.as_ref().and_then(|ctx| slice(v.ref_context_start, ctx)))
}

/// Reference bases from `lo` towards `hi`, upper-case, from the variant's prepared
/// reference (its event reference, else its `ref_context`), stopping where that
/// reference ends (a contig end). None when it does not hold `lo`.
pub(crate) fn reference_from(v: &Variant, lo: i64, hi: i64) -> Option<Vec<u8>> {
    let slice = |start: i64, seq: &str| {
        let a = lo - start;
        let b = (hi - start).min(seq.len() as i64);
        (a >= 0 && a < b).then(|| seq.as_bytes()[a as usize..b as usize].to_ascii_uppercase())
    };
    v.event_ref
        .as_ref()
        .and_then(|(start, seq)| slice(*start, seq))
        .or_else(|| v.ref_context.as_ref().and_then(|ctx| slice(v.ref_context_start, ctx)))
}

/// Reference bases from `lo` up to `hi`, upper-case, starting where the variant's
/// prepared reference starts if that is after `lo` (a contig start), with the
/// position of the first base. None when it does not hold the base before `hi`.
pub(crate) fn reference_to(v: &Variant, lo: i64, hi: i64) -> Option<(i64, Vec<u8>)> {
    let slice = |start: i64, seq: &str| {
        let a = (lo - start).max(0);
        let b = hi - start;
        (a < b && b <= seq.len() as i64)
            .then(|| (start + a, seq.as_bytes()[a as usize..b as usize].to_ascii_uppercase()))
    };
    v.event_ref
        .as_ref()
        .and_then(|(start, seq)| slice(*start, seq))
        .or_else(|| v.ref_context.as_ref().and_then(|ctx| slice(v.ref_context_start, ctx)))
}

/// Whether a read's own bases tell a pure indel's ALT from its REF, read where
/// they sit. Rightwards from the read's aligned base just left of the
/// discrimination window, or leftwards from its aligned base just right of it,
/// the read must read, unmasked, a base where the two alleles differ (the first
/// such base, or a later one when that is masked), every base up to there fitting
/// the ALT; a base below `min_baseq`, or N, fits anything. Unlike the REF rule
/// ([`read_is_informative`]) there is no margin base past it: that margin guards
/// CIGAR-only REF calls against a hidden terminal mismatch, while this reads the
/// deciding base. Soft-clipped bases are not read. Used for an ALT read that
/// spans neither ALT-side window: the CIGAR's gap alone is placement, and a
/// carrier ending inside the repeat holds only bases the alleles share, while a
/// truncated long insertion's carrier holds the inserted bases. A read starting
/// (or ending) on a flank base reads from it. Each reading uses the reference it
/// can hold, stopping at a contig end. A read with no aligned base on either side
/// of the window has nothing to read from: false. None when the read has a side
/// to read from but no prepared reference holds it (an unprepared variant), or
/// its deciding base lies past a contig edge (REF has no base there, but a
/// circular contig's continues at its start): it cannot be judged.
pub(crate) fn alt_bases_discriminate(record: &Record, v: &Variant, quals: &[u8], min_baseq: u8) -> Option<bool> {
    if !is_pure_indel(&v.ref_allele, &v.alt_allele) {
        return Some(true);
    }
    let (lo, hi) = change_interval(v);
    let (dlo, dhi) = (lo - 1, hi + 1);
    let (r, a) = (v.ref_allele.as_bytes(), v.alt_allele.as_bytes());
    let n = (a.len() as i64 - r.len() as i64).abs();
    let ins: Vec<u8> = if a.len() > r.len() { a[r.len()..].to_ascii_uppercase() } else { Vec::new() };
    let j0 = v.pos + r.len().min(a.len()) as i64;
    // Two bases past the window on the reading's far side are enough for either
    // allele to show the first differing base and one more (the deletion's
    // removed bases are taken from inside the window).
    let ext = 2;
    // The ALT haplotype of reference [from, to): the insert added, or the deleted
    // bases removed, at j0. None when the stretch stops short of the event (it
    // reaches a contig end inside it).
    let alt_of = |from: i64, seq: &[u8]| -> Option<Vec<u8>> {
        let k = usize::try_from(j0 - from).ok()?;
        if ins.is_empty() {
            let rest = k + n as usize;
            (rest <= seq.len()).then(|| [&seq[..k], &seq[rest..]].concat())
        } else {
            (k <= seq.len()).then(|| [&seq[..k], ins.as_slice(), &seq[k..]].concat())
        }
    };
    // An insertion's ALT stretch is n bases longer than the reference it is built
    // from: the REF stretch reaches as far, where the prepared reference holds it;
    // past the REF stretch (a contig end) a base decides nothing.
    let nn = if ins.is_empty() { 0 } else { n };
    let right = reference_from(v, dlo, dhi + ext + nn);
    let left = reference_to(v, dlo - ext - nn, dhi);
    let seq = record.seq().as_bytes();
    // A flank base read where the read starts (or ends) on it: unmasked and the
    // reference's. Masked, it would fit anything, and for a deletion sliding
    // through a repeat the flank is what tells the read from REF placed further
    // along the run.
    let flank_read = |q: usize, flank: Option<&u8>| -> bool {
        let b = seq[q].to_ascii_uppercase();
        b != b'N' && quals.get(q).is_some_and(|&x| x >= min_baseq) && flank == Some(&b)
    };
    // Where reading starts: rightwards from the read's aligned base just left of the
    // window, or from the flank base itself when the read starts on it and reads it
    // (the windows start at the flank too; the flank is a base both alleles share).
    // Leftwards likewise, from one past the read's last base to read.
    let ql = find_read_pos(record, dlo - 1).map(|q| q + 1).or_else(|| {
        find_read_pos(record, dlo).filter(|&q| flank_read(q, right.as_ref().and_then(|r| r.first())))
    });
    let qr = find_read_pos(record, dhi).or_else(|| {
        find_read_pos(record, dhi - 1)
            .filter(|&q| flank_read(q, left.as_ref().and_then(|(_, l)| l.last())))
            .map(|q| q + 1)
    });
    if ql.is_none() && qr.is_none() {
        return Some(false);
    }
    // A stretch cut short by the prepared reference's end (a contig edge): a read
    // still undecided where it ends has its deciding base past the edge.
    let right_cut = right.as_ref().is_some_and(|r| (r.len() as i64) < dhi + ext + nn - dlo);
    let left_cut = left.as_ref().is_some_and(|(_, l)| (l.len() as i64) < dhi - (dlo - ext - nn));
    let cig: Vec<Cigar> = record.cigar().iter().copied().collect();
    // Soft clips at either end (behind any hard clip) are not read.
    let clip = |ops: &mut dyn Iterator<Item = &Cigar>| -> usize {
        ops.take_while(|c| matches!(c, Cigar::SoftClip(_) | Cigar::HardClip(_)))
            .map(|c| if let Cigar::SoftClip(k) = c { *k as usize } else { 0 })
            .sum()
    };
    let lead = clip(&mut cig.iter());
    let trail = seq.len().saturating_sub(clip(&mut cig.iter().rev()));
    // The read's bases at query offsets `idx`, judged against `refh` and `alth`:
    // every unmasked base fits the ALT, up to a base at or past the first
    // differing one where the alleles differ, read unmasked. No margin base past
    // it: the REF rule's margin protects a CIGAR-only REF call from a hidden
    // terminal mismatch, while this reads the deciding base itself, on a read the
    // CIGAR already calls ALT. Reading stops where either stretch ends: Some(true)
    // ALT, Some(false) not, None when the read still had bases past the stretch's end.
    let holds = |idx: &mut dyn Iterator<Item = usize>, refh: &[u8], alth: &[u8]| -> Option<bool> {
        // The first base where the alleles differ: past the shorter stretch when
        // they agree over it (a REF stretch cut at a contig edge).
        let n = alth.len().min(refh.len());
        let d = refh.iter().zip(alth).position(|(x, y)| x != y).unwrap_or(n);
        for (i, q) in idx.enumerate() {
            if i >= n {
                return None;
            }
            let b = seq[q].to_ascii_uppercase();
            let masked = b == b'N' || quals.get(q).is_none_or(|&x| x < min_baseq);
            if masked {
                continue;
            }
            if b != alth[i] {
                return Some(false);
            }
            if i >= d && refh[i] != alth[i] {
                return Some(true);
            }
        }
        Some(false)
    };
    // A side the read can be read from but no reference holds, or whose deciding
    // base lies past a contig edge, leaves the read unjudged.
    let mut unjudged = false;
    if let Some(q) = ql {
        match right.as_deref().and_then(|r| alt_of(dlo, r).map(|a| (r, a))) {
            Some((right, alth)) => match holds(&mut (q..trail), right, &alth) {
                Some(true) => return Some(true),
                None if right_cut => unjudged = true,
                _ => {}
            },
            None => unjudged = true,
        }
    }
    if let Some(q) = qr {
        match left.as_ref().and_then(|(from, l)| alt_of(*from, l).map(|a| (l, a))) {
            Some((left, alth)) => {
                let refr: Vec<u8> = left.iter().rev().copied().collect();
                let altr: Vec<u8> = alth.iter().rev().copied().collect();
                match holds(&mut (lead..q).rev(), &refr, &altr) {
                    Some(true) => return Some(true),
                    None if left_cut => unjudged = true,
                    _ => {}
                }
            }
            None => unjudged = true,
        }
    }
    if unjudged { None } else { Some(false) }
}

/// Whether the read's bases between its nearest aligned base left of a pure indel's
/// discrimination window and its nearest aligned base right of it (its own gap may
/// cover a flank base) are exactly the ALT over that stretch (bases below
/// `min_baseq`, or N, fit anything; at least one base read), and wherever the
/// haplotype its own alignment proposes differs from the ALT the read shows the
/// ALT's base unmasked (else it fits another allele as well). Judges a read whose
/// CIGAR writes the event somewhere it gives another haplotype: compensating
/// mismatches can make its bases the ALT nonetheless.
pub(crate) fn read_spells_alt(record: &Record, v: &Variant, quals: &[u8], min_baseq: u8) -> bool {
    read_spells_alt_checked(record, v, quals, min_baseq, true)
}

/// `read_spells_alt` without the check that the read says, unmasked, which allele
/// where its own alignment proposes another haplotype: for a read whose gap covers
/// the anchor, where the alignment's placement is a tie-break, not evidence.
pub(crate) fn read_bases_fit_alt(record: &Record, v: &Variant, quals: &[u8], min_baseq: u8) -> bool {
    read_spells_alt_checked(record, v, quals, min_baseq, false)
}

fn read_spells_alt_checked(record: &Record, v: &Variant, quals: &[u8], min_baseq: u8, check_claimed: bool) -> bool {
    if !is_pure_indel(&v.ref_allele, &v.alt_allele) {
        return false;
    }
    let (lo, hi) = change_interval(v);
    let (dlo, dhi) = (lo - 1, hi + 1);
    let (r, a) = (v.ref_allele.as_bytes(), v.alt_allele.as_bytes());
    // The read's nearest aligned bases outside the window: its own (misplaced)
    // gap may cover a flank base.
    let reach = (r.len() as i64 - a.len() as i64).abs() + 5;
    let Some(left) = (dlo - 1 - reach..dlo).rev().find(|&p| find_read_pos(record, p).is_some()) else {
        return false;
    };
    let Some(right) = (dhi..dhi + reach).find(|&p| find_read_pos(record, p).is_some()) else {
        return false;
    };
    let (Some(ql), Some(qr)) = (find_read_pos(record, left), find_read_pos(record, right)) else {
        return false;
    };
    let Some(refw) = reference_span(v, left + 1, right) else {
        return false;
    };
    let j = (v.pos + r.len().min(a.len()) as i64 - (left + 1)) as usize;
    let altw: Vec<u8> = if a.len() > r.len() {
        [&refw[..j], &a[r.len()..].to_ascii_uppercase(), &refw[j..]].concat()
    } else {
        [&refw[..j], &refw[j + r.len() - a.len()..]].concat()
    };
    if qr <= ql || qr - ql - 1 != altw.len() {
        return false;
    }
    // The haplotype the read's own alignment proposes between the anchors: the
    // reference where it aligns, its bases where it inserts.
    let seq = record.seq().as_bytes();
    let mut claimed: Vec<u8> = Vec::with_capacity(altw.len());
    let (mut rp, mut qp) = (record.pos(), 0usize);
    for op in record.cigar().iter() {
        match op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                for k in 0..*n as usize {
                    if qp + k > ql && qp + k < qr {
                        claimed.push(refw[(rp + k as i64 - (left + 1)) as usize]);
                    }
                }
                rp += *n as i64;
                qp += *n as usize;
            }
            Cigar::Ins(n) | Cigar::SoftClip(n) => {
                for k in 0..*n as usize {
                    if qp + k > ql && qp + k < qr {
                        claimed.push(seq[qp + k].to_ascii_uppercase());
                    }
                }
                qp += *n as usize;
            }
            Cigar::Del(n) | Cigar::RefSkip(n) => rp += *n as i64,
            _ => {}
        }
    }
    if claimed.len() != altw.len() {
        return false;
    }
    // Every unmasked base fits the ALT, and wherever the alignment's haplotype and
    // the ALT disagree the read must say which, with an unmasked base: a masked one
    // there fits both alleles (an N where "deleted a C" and "deleted a G" differ).
    let mut read_any = false;
    for (i, q) in (ql + 1..qr).enumerate() {
        let b = seq[q].to_ascii_uppercase();
        let masked = b == b'N' || quals.get(q).is_none_or(|&x| x < min_baseq);
        if masked {
            if check_claimed && claimed[i] != altw[i] {
                return false;
            }
            continue;
        }
        if b != altw[i] {
            return false;
        }
        read_any = true;
    }
    read_any
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
    let (start, end) = (record.pos(), ref_end(record));
    windows.iter().any(|&(lo, hi)| start <= lo && end >= hi)
}

/// Whether an ALT read's CIGAR alone settles a pure indel: its aligned reference
/// extent spans one of the informative windows as seen from the ALT haplotype,
/// `[lo - 1, hi + 2)` or `[lo - 2, hi + 1)`. An insertion's are its REF windows (its
/// carrier's extent leaves out the inserted bases, so spanning them it reads more,
/// not less). A deletion's carrier's extent counts its own gap, so spanning the
/// REF windows it may hold up to `L - 1` fewer bases than it needs to show where
/// the alleles differ: read on the ALT haplotype, a deletion's windows are an
/// insertion's. A window that starts or ends on the read's own end base also
/// needs that flank base read (unmasked, the reference's). Always true for other
/// variants.
pub(crate) fn alt_read_is_informative(record: &Record, v: &Variant, quals: &[u8], min_baseq: u8) -> bool {
    if !is_pure_indel(&v.ref_allele, &v.alt_allele) {
        return true;
    }
    let (lo, hi) = change_interval(v);
    let (start, end) = (record.pos(), ref_end(record));
    // A window that starts (or ends) on the read's own first (or last) base needs
    // that flank base read: for a deletion sliding through a repeat, it is all that
    // tells the read from REF placed one base along the run. Unverifiable without a
    // prepared reference, where the extent stands as before.
    let flank_read = |pos: i64| -> bool {
        let Some(q) = find_read_pos(record, pos) else { return false };
        let Some(r) = reference_span(v, pos, pos + 1) else { return true };
        let b = record.seq()[q].to_ascii_uppercase();
        b != b'N' && quals.get(q).is_some_and(|&x| x >= min_baseq) && b == r[0]
    };
    // [lo - 1, hi + 2): from the left flank; [lo - 2, hi + 1): to the right flank
    // (`end` is exclusive, so a read whose last base is the flank ends at hi + 1).
    let from_left = start < lo && end >= hi + 2 && (start < lo - 1 || flank_read(lo - 1));
    let to_right = start < lo - 1 && end > hi && (end > hi + 1 || flank_read(hi));
    from_left || to_right
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

    #[test]
    fn allele_kind_reads_lengths_and_the_anchor_case_insensitively() {
        use AlleleKind::*;
        assert_eq!(allele_kind("A", "T"), Some(Snv));
        assert_eq!(allele_kind("AC", "GT"), Some(Mnp));
        assert_eq!(allele_kind("A", "ACC"), Some(Insertion));
        assert_eq!(allele_kind("a", "ACC"), Some(Insertion));
        assert_eq!(allele_kind("ACC", "a"), Some(Deletion));
        assert_eq!(allele_kind("C", "TA"), Some(Complex)); // the anchor changes
        assert_eq!(allele_kind("GC", "T"), Some(Complex));
        assert_eq!(allele_kind("GCA", "TT"), Some(Complex));
        assert_eq!(allele_kind("", "A"), None);
        assert_eq!(allele_kind("A", ""), None);
    }

    /// A variant with a reference context starting at 0.
    fn var(ctx: &str, pos: i64, r: &str, a: &str) -> Variant {
        Variant::new(
            "1".into(), pos, r.into(), a.into(), "X".into(),
            Some(ctx.into()), 0, 0, None, None, None, None,
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
        let v = Variant::new("1".into(), 2, "CA".into(), "C".into(), "DELETION".into(), None, 0, 0, None, None, None, None);
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
    fn the_scan_reaches_every_placement_and_no_further() {
        // A 2bp insertion over a 20-junction region: every junction to 30.
        let mut ins = var(HOMO, 9, "G", "GAC");
        ins.shift_region = Some((10, 30));
        assert_eq!(scan_window(&ins, 5), (4, 30));
        // A 4bp deletion sliding over [10, 30): starts up to 26, not into its
        // own deleted span past that.
        let mut del = var(HOMO, 9, "GACGT", "G");
        del.shift_region = Some((10, 30));
        assert_eq!(scan_window(&del, 5), (4, 26));
        // Inside the repeat-span reach, the region changes nothing.
        del.shift_region = Some((10, 16));
        assert_eq!(scan_window(&del, 5), (4, 14));
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
