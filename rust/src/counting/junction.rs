//! The junction judge for pure deletions of [`LONG_DELETION`] bases or more (RJ-24).
//!
//! A read is REF or ALT for such a deletion by its own bases at the POS breakpoint,
//! however the aligner wrote it (a D op, a soft clip, a split read):
//!
//! - **One window per allele, of equal length, read from the left flank.** [`K`]
//!   flank bases ending at the anchor, the microhomology (the bases the deletion can
//!   slide right through, which both alleles share), then [`K`] bases of the
//!   allele's own continuation: the deleted bases for REF, the right flank for ALT.
//!   A REF and an ALT molecule hold their windows from the same read starts, so
//!   neither allele is favoured by where reads start.
//! - **Nothing past the window is read**, for either allele: checking REF through
//!   every deleted base it covers would expose it to errors the ALT reading never
//!   meets (REJ-20261009-003).
//!
//! - **Placement.** A read is read from the query position of its aligned first
//!   window base, so both windows start where the read's flank does. A read whose
//!   left flank is not aligned there (in its soft clip, under its own deletion, or
//!   aligned over the deleted bases' other copy) is placed from the right flank,
//!   which both alleles share, as the ALT's window sits there (RJ-13); REF's window
//!   ends inside the deleted bases, which no ALT molecule holds. Its confident bases
//!   (at or above `--min-baseq`, not N) must favour the left flank over the deleted
//!   bases REF has before the right flank by one base read at `--min-baseq` (the
//!   base next to the junction always differs: the deletion is placed leftmost). It
//!   is ALT where it fits, else another allele; aligned over the other copy, only
//!   where it fits. Hard-clipped bases are absent: such a record holds no window on
//!   its own.
//! - **Masked bases.** Bases below `--min-baseq` (and N) fit anything, the one
//!   quality rule every backend shares; so do the window offsets the neighbouring
//!   change guard masks. An ALT call also needs the read's bases to favour ALT over
//!   REF by one minimum-quality base's evidence (RJ-16).
//! - **Outcomes.** A read holding the windows is REF or ALT when it fits exactly one;
//!   fitting neither, it carries another allele (partial evidence when closer to
//!   ALT). A read that does not hold them, or whose windows cross a splice, is depth
//!   only. A decided read whose windows read soft-clipped bases inside its fragment
//!   counts toward depth in DNA (clip admission). A read that fits REF but whose own
//!   alignment deletes or inserts [`MIN_OWN_INDEL`] bases or more inside the deleted
//!   span carries another allele (C26): partial evidence.

use log::trace;
use rust_htslib::bam::record::{Cigar, Record};

use crate::shared::bam_utils::fragment_query_span;
use crate::types::Variant;

use super::carrier::{
    aligned_query_range, confident_weight, one_base_evidence, readable_range, score, Guard, Placement, ReadRules,
    Reading, GUARD_MIN_READS, GUARD_MIN_SHARE, GUARD_REF_SHARE, WILD,
};
use super::observed::{canonical, ObservedAllele};
use super::utils::{find_read_pos, median_qual, ClassifyPhase, ClassifyResult};
use super::variant_checks::deleted_anchor_query;
use super::window::{allele_kind, AlleleKind};

/// Deleted bases from which a pure deletion is judged at its junction.
pub(crate) const LONG_DELETION: usize = 50;
/// Clean bases each window holds before the junction (the flank) and past the
/// microhomology (the allele's own continuation).
const K: usize = 4;
/// A read's own indel at least this long inside the deleted span makes it another
/// allele (C26); shorter ones in a deleted repeat are stutter or alignment noise, as
/// the CIGAR checks treat 1-4 bp deletions in their scan window.
const MIN_OWN_INDEL: u32 = 5;

/// A long deletion's junction windows on the prepared reference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Junction {
    /// Genome position of the first deleted base at the leftmost placement.
    lo: i64,
    /// Deleted bases.
    len: usize,
    /// Bases the deletion slides right through (the microhomology).
    mh: usize,
    /// REF's window: K flank bases, the microhomology, K deleted bases.
    pub(crate) ref_w: Vec<u8>,
    /// ALT's window: the same flank and microhomology, K right-flank bases.
    pub(crate) alt_w: Vec<u8>,
    /// REF read back from ALT's window's last base: the K deleted bases before the
    /// right flank where ALT has the left flank, then the bases ALT's window ends with.
    ref_back: Vec<u8>,
}

impl Junction {
    /// Genome position of both windows' first base (the flank's first).
    pub(crate) fn first(&self) -> i64 {
        self.lo - K as i64
    }

    /// Genome position of an allele's window base at `offset` (0-based).
    pub(crate) fn genome(&self, alt: bool, offset: usize) -> i64 {
        if alt && offset >= K {
            self.lo + self.len as i64 + (offset - K) as i64
        } else {
            self.first() + offset as i64
        }
    }

    /// Genome position of an allele's window's last base.
    pub(crate) fn last(&self, alt: bool) -> i64 {
        self.genome(alt, self.ref_w.len() - 1)
    }

    /// The bases where the alleles differ over every placement, genomic [lo, hi):
    /// the deleted bases and the microhomology.
    fn span(&self) -> (i64, i64) {
        (self.lo, self.lo + (self.len + self.mh) as i64)
    }
}

/// Whether the junction judge reads `v`: a pure deletion of [`LONG_DELETION`] bases
/// or more, prepared against a reference (prep's `event_ref`).
pub(crate) fn judges(v: &Variant) -> bool {
    v.event_ref.is_some()
        && allele_kind(&v.ref_allele, &v.alt_allele) == Some(AlleleKind::Deletion)
        && v.ref_allele.len() - v.alt_allele.len() >= LONG_DELETION
}

/// The deletion's junction windows on its prepared reference. None when it is not
/// judged here, or the reference does not hold the windows.
pub(crate) fn junction(v: &Variant) -> Option<Junction> {
    if !judges(v) {
        return None;
    }
    let (start, seq) = v.event_ref.as_ref()?;
    let reference: Vec<u8> = seq.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let at = |g: i64| usize::try_from(g - start).ok().filter(|&i| i < reference.len());
    let len = v.ref_allele.len() - v.alt_allele.len();
    // The deleted bases follow the anchor; slide them to their leftmost placement,
    // then count how far they slide right (prep left-aligns, so the left slide is
    // a guard, not the rule).
    let mut lo = at(v.pos + v.alt_allele.len() as i64)?;
    while lo > 0 && reference.get(lo - 1) == reference.get(lo - 1 + len) {
        lo -= 1;
    }
    let mut mh = 0;
    while lo + mh + len < reference.len() && reference[lo + mh] == reference[lo + mh + len] {
        mh += 1;
    }
    let reach = mh + K;
    if lo < K || lo + len + reach > reference.len() {
        return None;
    }
    let flank = &reference[lo - K..lo];
    let ref_w = [flank, &reference[lo..lo + reach]].concat();
    let alt_w = [flank, &reference[lo + len..lo + len + reach]].concat();
    let ref_back = reference[lo + len - K..lo + len + reach].to_vec();
    Some(Junction { lo: start + lo as i64, len, mh, ref_w, alt_w, ref_back })
}

/// The genome span REF's window reads, [its first base, its last + 1): the flank, the
/// microhomology and K deleted bases. Nothing past it is read, for either allele.
/// None when the junction judge does not judge `v`.
pub(crate) fn ref_window_span(v: &Variant) -> Option<(i64, i64)> {
    let j = junction(v)?;
    Some((j.first(), j.last(false) + 1))
}

/// A read's call at a long deletion (RJ-24), or None when the junction judge does
/// not judge the variant (the caller then uses the CIGAR checks).
pub(crate) fn classify(
    record: &Record,
    v: &Variant,
    quals: &[u8],
    min_baseq: u8,
    rules: &ReadRules,
) -> Option<ClassifyResult> {
    let j = junction(v)?;
    let seq = record.seq().as_bytes();
    if seq.is_empty() || quals.len() < seq.len() {
        // A record stored without its bases (SEQ '*') shows nothing here.
        return Some(ClassifyResult::no_coverage(ClassifyPhase::MaskedCompare));
    }
    let Some(readable) = readable_range(record, rules) else {
        return Some(ClassifyResult::no_coverage(ClassifyPhase::MaskedCompare));
    };
    // Offsets the guard masks (a change recurring beside the event, a separate
    // event) fit any base, in each allele's window at its own genome position.
    let masked = rules.guard.and_then(|g| g.get(v)).map_or(&[][..], |g| g.masked.as_slice());
    let window = |alt: bool, w: &[u8]| -> Vec<u8> {
        w.iter()
            .enumerate()
            .map(|(o, &b)| if masked.contains(&j.genome(alt, o)) { WILD } else { b })
            .collect()
    };
    let (ref_w, alt_w) = (window(false, &j.ref_w), window(true, &j.alt_w));
    let n = ref_w.len();
    let read_from = |s: i64| -> Option<(Reading, Reading)> {
        let s = usize::try_from(s).ok().filter(|&s| s >= readable.0 && s + n <= readable.1)?;
        let (bases, qs, span) = (&seq[s..s + n], &quals[s..s + n], (s, s + n));
        Some((score(bases, qs, min_baseq, &ref_w, span), score(bases, qs, min_baseq, &alt_w, span)))
    };
    // The windows' reading, and how it was placed.
    let reading = match find_read_pos(record, j.first()) {
        Some(q) => read_from(q as i64).map(|p| (p, Placement::Flank)),
        None => {
            // The left flank not aligned where the windows start: placed from the
            // right flank as the ALT's window sits there, when the read's confident
            // bases favour the left flank over the deleted bases REF has there by one
            // base read at --min-baseq, and by more than nothing (at --min-baseq 0
            // that base weighs nothing). A masked flank offset weighs nothing for
            // ALT; the deleted bases it is weighed against were never masked.
            find_read_pos(record, j.last(true)).and_then(|q| {
                let s = q as i64 - (n as i64 - 1);
                let evidence = confident_weight(&seq, quals, min_baseq, &j.ref_back, s, readable)
                    - confident_weight(&seq, quals, min_baseq, &alt_w, s, readable);
                if evidence <= 0.0 || evidence < one_base_evidence(min_baseq) {
                    trace!("junction: placed from the right flank, its bases do not show the left flank ({:.2}) → depth only", evidence);
                    return None;
                }
                let (r, a) = read_from(s)?;
                // Held: the left flank in the read's clip, or under its own deletion;
                // else aligned over the deleted bases' other copy.
                let held = aligned_query_range(record).is_some_and(|(first, _)| s < first as i64)
                    || deleted_anchor_query(record, j.first()).is_some();
                let placement = match (a.mismatches == 0, held) {
                    (true, true) => Placement::OtherFlank,
                    (true, false) => Placement::ThroughAligned,
                    (false, true) => Placement::Unfit { alt: true },
                    (false, false) => return None,
                };
                Some(((r, a), placement))
            })
        }
    };
    let Some(((r, a), placement)) = reading else {
        // The read does not hold the windows: depth only, no allele, no mFSD class.
        let mut out = ClassifyResult::neither(ClassifyPhase::MaskedCompare);
        out.uninformative = true;
        return Some(out);
    };
    let (lo, hi) = r.span;
    if crosses_splice(record, lo, hi) {
        // Spliced inside the windows: the read skips bases the windows compare.
        let mut out = ClassifyResult::neither(ClassifyPhase::CigarRecon);
        out.uninformative = true;
        return Some(out);
    }
    let qual = median_qual(&quals[lo..hi], min_baseq);
    // The read's evidence for ALT over REF (log10), every base weighed by its quality.
    let evidence = r.mismatch_weight - a.mismatch_weight;
    let (fits_ref, fits_alt) = (r.mismatches == 0, a.mismatches == 0);
    let mut out = if placement == (Placement::Unfit { alt: true }) {
        // Placed as the ALT, which its bases do not fit: another allele.
        ClassifyResult::neither_with_nearby(qual, ClassifyPhase::MaskedCompare)
    } else if fits_alt && !fits_ref && evidence < one_base_evidence(min_baseq) {
        ClassifyResult::neither(ClassifyPhase::MaskedCompare)
    } else if fits_alt && !fits_ref {
        ClassifyResult::is_alt(qual, ClassifyPhase::MaskedCompare)
    } else if fits_ref && !fits_alt {
        ClassifyResult::is_ref(qual, ClassifyPhase::MaskedCompare)
    } else if a.mismatches < r.mismatches {
        // Another allele, closer to ALT: partial evidence.
        ClassifyResult::neither_with_nearby(qual, ClassifyPhase::MaskedCompare)
    } else {
        ClassifyResult::neither(ClassifyPhase::MaskedCompare)
    };
    if out.is_ref && own_indel_inside(record, j.span()) {
        // The read's own alignment deletes or inserts bases inside the deleted span:
        // another allele (C26), whatever its window shows.
        out = ClassifyResult::neither_with_nearby(qual, ClassifyPhase::MaskedCompare);
    }
    out.has_n_base = r.had_n || a.had_n;
    // A decided read whose windows read soft-clipped bases counts toward depth (in
    // DNA) when every base read lies inside its fragment: past it a clip is adapter.
    let clipped = aligned_query_range(record).is_some_and(|(first, after)| lo < first || hi > after);
    out.clip_admissible = (out.is_ref || out.is_alt)
        && clipped
        && fragment_query_span(record).is_some_and(|(f0, f1)| lo >= f0 && hi <= f1);
    // Decided from aligned bases placed from the other side: the read starts past
    // the variant position as aligned, while the same molecule of the other allele
    // would span it, so it counts as that one would (as RJ-23's right reading).
    out.read_from_right = (out.is_ref || out.is_alt) && placement != Placement::Flank && !clipped;
    Some(out)
}

/// Bases past the windows that tell a read's allele where the [`K`] bases past the
/// microhomology do not (a run among them), and that the bases after a run must
/// rejoin its allele's continuation over.
const REJOIN: usize = 12;
/// Bases past the windows the guard reads: a run running on past them, and the
/// [`REJOIN`] after it.
const PAST: usize = 2 * REJOIN;

/// The neighbouring-change guard for a long deletion (the carrier rule's, RJ-23,
/// over the junction windows): from the reads `counted` admits, once per fragment,
/// each read holding the windows from its aligned first window base is taken as the
/// allele whose continuation it follows (its [`K`] bases past the microhomology, or,
/// where it reads them, the [`REJOIN`] bases past the windows, fit that allele to
/// within one base and the other by two or more; the two disagreeing, neither), and
/// its confident changes are tallied per window offset, on the base before the
/// windows, and on [`PAST`] bases on along its allele past them (up to a splice). A
/// change recurs on at least [`GUARD_MIN_READS`] fragments and [`GUARD_MIN_SHARE`] of
/// that allele's reads. Masked at its offset in both windows are, in the flank and
/// microhomology both alleles share, the changes the REF reads carry (a germline SNP
/// or MNP, background on both haplotypes, runs included, on [`GUARD_MIN_READS`] REF
/// reads or more); each isolated change of the ALT reads' own (no own change at
/// either neighbouring base); past the junction each isolated change on the REF
/// reads; and each run past the junction, on either allele's reads, that is a
/// phased neighbour ([`neighbour`]: the bases after it rejoin the allele's
/// continuation, and no other breakpoint gives them). Any other run of the ALT reads'
/// own changes is an indel at the junction, and any other run past it on reads
/// fitting REF another breakpoint: another allele, not masked. Masks never take every
/// offset where REF and ALT differ: those stay read, and a read changed there is
/// neither allele. On the ALT reads alone (under [`GUARD_REF_SHARE`] of REF's)
/// isolated recurring changes and phased runs also make the larger allele the
/// diagnostic names (the deletion with the changes, through a run's last base past
/// the windows), and beside any other run none is named. None when nothing is masked
/// or named, or the variant is not judged here.
pub(crate) fn guard_for(
    read_cache: &[Record],
    v: &Variant,
    counted: &dyn Fn(&Record) -> bool,
    min_baseq: u8,
) -> Option<Guard> {
    let j = junction(v)?;
    let (start, seq) = v.event_ref.as_ref()?;
    let reference: Vec<u8> = seq.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let at = |g: i64| usize::try_from(g - start).ok().filter(|&i| i < reference.len());
    let base = |g: i64| at(g).map(|i| reference[i]);
    let n = j.ref_w.len();
    let cont = K + j.mh;
    // Each allele's bases by edged offset e: the base before the windows (e 0), the
    // window (window offset e - 1), then PAST bases on along the allele (the deleted
    // bases for REF, the right flank for ALT); None past the prepared reference.
    let edged = |alt: bool, w: &[u8]| -> Vec<Option<u8>> {
        let past = (1..=PAST as i64).map(|k| base(j.last(alt) + k));
        [base(j.first() - 1)].into_iter().chain(w.iter().map(|&b| Some(b))).chain(past).collect()
    };
    let windows = [edged(false, &j.ref_w), edged(true, &j.alt_w)];
    let m = n + 1 + PAST;
    // Per allele (REF 0, ALT 1) and edged offset: fragments reading it, changing it,
    // and the bases the changes show (A C G T).
    let (mut cover, mut change) = ([vec![0u32; m], vec![0u32; m]], [vec![0u32; m], vec![0u32; m]]);
    let mut shown = [vec![[0u32; 4]; m], vec![[0u32; 4]; m]];
    let mut fragments: std::collections::HashSet<&[u8]> = std::collections::HashSet::new();
    let mut scanned = 0u32;
    for record in read_cache {
        if !counted(record) {
            continue;
        }
        let bases = record.seq().as_bytes();
        let quals = record.qual();
        let Some(q) = find_read_pos(record, j.first()).filter(|&q| q + n <= bases.len() && quals.len() >= q + n)
        else {
            continue;
        };
        // Past the windows a read is read up to its next splice: the next exon is not
        // its allele's continuation.
        let stop = splice_after(record, q).unwrap_or(bases.len());
        // The read's base at edged offset e, and whether it is read confidently.
        let read = |e: usize| -> Option<(u8, bool)> {
            let i = (q + e).checked_sub(1).filter(|&i| e <= n || i < stop)?;
            let b = bases.get(i)?.to_ascii_uppercase();
            Some((b, b != WILD && *quals.get(i)? >= min_baseq))
        };
        let off = |w: &[Option<u8>], e: usize| matches!((read(e), w[e]), (Some((b, true)), Some(x)) if b != x);
        let told = |es: std::ops::RangeInclusive<usize>| {
            let misses = |w: &[Option<u8>]| es.clone().filter(|&e| off(w, e)).count();
            match (misses(&windows[0]), misses(&windows[1])) {
                (0..=1, 2..) => Some(0),
                (2.., 0..=1) => Some(1),
                _ => None,
            }
        };
        let past = n + 1..=n + REJOIN;
        let read_on = past.clone().all(|e| read(e).is_some()).then(|| told(past)).flatten();
        let allele = match (told(cont + 1..=n), read_on) {
            (Some(a), Some(b)) if a != b => continue,
            (Some(a), _) | (None, Some(a)) => a,
            (None, None) => continue,
        };
        if !fragments.insert(record.qname()) {
            continue;
        }
        scanned += 1;
        for e in 0..m {
            let (Some((b, _)), Some(_)) = (read(e), windows[allele][e]) else {
                continue;
            };
            cover[allele][e] += 1;
            if off(&windows[allele], e) {
                change[allele][e] += 1;
                if let Some(k) = b"ACGT".iter().position(|&x| x == b) {
                    shown[allele][e][k] += 1;
                }
            }
        }
    }
    let share = |a: usize, e: usize| if cover[a][e] == 0 { 0.0 } else { change[a][e] as f64 / cover[a][e] as f64 };
    let recurs = |a: usize, e: usize| change[a][e] >= GUARD_MIN_READS && share(a, e) >= GUARD_MIN_SHARE;
    let on_ref = |e: usize| share(0, e) >= GUARD_REF_SHARE;
    // Background, a change the REF reads carry too (a germline SNP or MNP), only
    // where both windows read the same genome bases (the flank and the
    // microhomology): recurring on the REF reads, or on the ALT reads and on as many
    // REF reads as a recurring change needs (one REF read's error is no background).
    let background = |e: usize| {
        e <= cont && (recurs(0, e) || (recurs(1, e) && change[0][e] >= GUARD_MIN_READS && on_ref(e)))
    };
    // The ALT reads' own changes, isolated when no own change sits beside them; past
    // the junction a REF reads' change, when none recurs beside it on the REF reads.
    let own = |e: usize| recurs(1, e) && !background(e);
    let isolated = |e: usize| !own(e - 1) && !own(e + 1);
    let ref_isolated = |e: usize| recurs(0, e) && !recurs(0, e - 1) && !recurs(0, e + 1);
    // Each allele's runs reaching the windows (edged [first, last], two or more own
    // changes side by side): the ALT reads' own, and past the junction the REF
    // reads'. Those that are phased neighbours are masked; any other is another
    // allele.
    let own_of = |a: usize, e: usize| if a == 1 { own(e) } else { e > cont && recurs(0, e) };
    let runs_of = |a: usize| -> Vec<(usize, usize)> {
        let mut out: Vec<(usize, usize)> = Vec::new();
        for e in (0..m).filter(|&e| own_of(a, e)) {
            match out.last_mut() {
                Some(r) if r.1 + 1 == e => r.1 = e,
                _ => out.push((e, e)),
            }
        }
        out.retain(|&(e0, e1)| e1 > e0 && e0 <= n);
        out
    };
    let phased = |a: usize, run: (usize, usize)| {
        let most = |e: usize| b"ACGT"[(0..4).max_by_key(|&k| shown[a][e][k]).unwrap_or(0)];
        neighbour(&j, &base, &windows[a], a == 1, run, &|e| cover[a][e], &|e| recurs(a, e), &most)
    };
    let runs = [runs_of(0), runs_of(1)];
    let neighbours: [Vec<(usize, usize)>; 2] = [0, 1].map(|a| runs[a].iter().copied().filter(|&r| phased(a, r)).collect());
    let in_neighbour = |a: usize, e: usize| neighbours[a].iter().any(|&(e0, e1)| (e0..=e1).contains(&e));
    let mut offsets: Vec<usize> = (1..=n)
        .filter(|&e| {
            background(e)
                || (own(e) && isolated(e))
                || (e > cont && ref_isolated(e))
                || in_neighbour(0, e)
                || in_neighbour(1, e)
        })
        .collect();
    let tells = |e: &usize| j.ref_w[e - 1] != j.alt_w[e - 1];
    if (1..=n).filter(tells).all(|e| offsets.contains(&e)) {
        // Masking them all would leave the windows identical: they stay read.
        offsets.retain(|e| !tells(e));
    }
    let alt_only: Vec<usize> =
        (1..m).filter(|&e| recurs(1, e) && !on_ref(e) && (e <= n || in_neighbour(1, e))).collect();
    let mut masked = Vec::new();
    for &e in &offsets {
        for p in [j.genome(false, e - 1), j.genome(true, e - 1)] {
            if !masked.contains(&p) {
                masked.push(p);
            }
        }
    }
    // The larger allele: the ALT haplotype from the windows' first base through its
    // farthest change, and at least through the junction (a change in the flank is
    // carried with the deletion), each changed base as most ALT reads show it,
    // against the genome over the same span (the deleted bases included). Beside a
    // run that is no phased neighbour, none: the tallies cannot spell another allele
    // (an MNP at the junction, or an indel shifting the bases past it); the observed
    // allele names it from the reads' alignments.
    let run = runs[1].len() > neighbours[1].len();
    let larger = alt_only.iter().max().filter(|_| !run).and_then(|&far| {
        let far = (far - 1).max(K);
        let alt: Vec<u8> = (0..=far)
            .map(|o| match alt_only.contains(&(o + 1)) {
                true => Some(b"ACGT"[(0..4).max_by_key(|&k| shown[1][o + 1][k]).unwrap_or(0)]),
                false => windows[1][o + 1],
            })
            .collect::<Option<_>>()?;
        let (g0, g1) = (j.first(), j.genome(true, far));
        let ref_part = reference.get(at(g0)?..=at(g1)?)?.to_vec();
        let carriers = alt_only.iter().map(|&e| change[1][e]).min().unwrap_or(0);
        let given = alt_only.iter().map(|&e| cover[1][e] - change[1][e]).max().unwrap_or(0);
        let text = |x: Vec<u8>| String::from_utf8_lossy(&x).into_owned();
        canonical(g0, ref_part, alt, &base).map(|(pos, r, a)| ObservedAllele {
            pos,
            ref_allele: text(r),
            alt_allele: text(a),
            carriers,
            given_carriers: given,
        })
    });
    if masked.is_empty() && larger.is_none() {
        return None;
    }
    trace!(
        "junction guard at {}:{} {}>{}: {} window offsets masked, {} on the ALT reads alone, {} phased runs ({} fragments read)",
        v.chrom, v.pos + 1, v.ref_allele, v.alt_allele, offsets.len(), alt_only.len(),
        neighbours[0].len() + neighbours[1].len(), scanned,
    );
    Some(Guard::new(v, masked, larger, scanned))
}

/// Whether a run of an allele's own changes (edged offsets [e0, e1], `run`) is a phased
/// neighbour of the deletion, not another breakpoint (#285). It lies past the junction
/// (`windows` is the allele's edged bases, `cover` and `recurs` its reads' tallies,
/// `most` the base most of them show at a changed offset), and:
/// - **its bases rejoin the allele's continuation** after it: each of the [`REJOIN`]
///   bases after its last is read by [`GUARD_MIN_READS`] of the allele's fragments or
///   more, and no change recurs on any;
/// - **at the given breakpoint only**: the allele as its reads show it (the run's
///   bases as most show them, each other recurring change not compared), from the
///   windows' first base through those [`REJOIN`] bases, is within one base of no
///   other breakpoint's: neither the reference with one deletion at another start or
///   of another length, nor the allele with one more insertion or deletion from the
///   junction through the run, within the reference held. A run at the junction whose
///   bases are the deleted bases' own is the deletion shifted right; one before a
///   repeat that a deletion of another length gives whole is that deletion.
#[allow(clippy::too_many_arguments)]
fn neighbour(
    j: &Junction,
    base: &dyn Fn(i64) -> Option<u8>,
    windows: &[Option<u8>],
    alt: bool,
    (e0, e1): (usize, usize),
    cover: &dyn Fn(usize) -> u32,
    recurs: &dyn Fn(usize) -> bool,
    most: &dyn Fn(usize) -> u8,
) -> bool {
    let cont = K + j.mh;
    let after = e1 + 1..=e1 + REJOIN;
    if e0 <= cont || *after.end() >= windows.len() || !after.clone().all(|e| cover(e) >= GUARD_MIN_READS && !recurs(e)) {
        return false;
    }
    // The allele's base at window offset o, on along it past the windows.
    let own = |o: usize| base(j.genome(alt, o));
    let span = e1 + REJOIN;
    let shows: Vec<Option<u8>> = (0..span)
        .map(|o| match o + 1 {
            e if (e0..=e1).contains(&e) => Some(most(e)),
            e if recurs(e) => None,
            _ => own(o),
        })
        .collect();
    // Whether another breakpoint's bases (Some(None): an inserted base, any; None: off
    // the reference) differ from the reads' at one base or none. The given allele
    // differs at every base of the run.
    let gives = |other: &dyn Fn(usize) -> Option<Option<u8>>| -> bool {
        let mut differ = 0;
        for (o, &shown) in shows.iter().enumerate() {
            let Some(x) = other(o) else { return false };
            if let (Some(s), Some(x)) = (shown, x) {
                differ += usize::from(s != x);
                if differ > 1 {
                    return false;
                }
            }
        }
        true
    };
    let first = j.first();
    let held = (0..).take_while(|&d| base(first + (span + d) as i64).is_some()).count();
    // The reference with one deletion of d bases at genome x (none: d 0).
    let deletion = (first..=first + span as i64).any(|x| {
        (0..held).any(|d| {
            gives(&|o: usize| {
                let g = first + o as i64;
                base(if g < x { g } else { g + d as i64 }).map(Some)
            })
        })
    });
    // The allele with one more deletion of k bases, or insertion of k inside the run,
    // at window offset p from the junction through the run.
    let indel = (cont..e1).any(|p| {
        (1..held).any(|k| gives(&|o: usize| own(if o < p { o } else { o + k }).map(Some)))
            || (1..=e1 - p).any(|k| {
                gives(&|o: usize| match o {
                    o if o < p => own(o).map(Some),
                    o if o < p + k => Some(None),
                    o => own(o - k).map(Some),
                })
            })
    });
    !deletion && !indel
}

/// Query index of the read's first splice (N) after query index `q`, if any.
fn splice_after(record: &Record, q: usize) -> Option<usize> {
    let mut i = 0usize;
    for op in record.cigar().iter() {
        match op {
            Cigar::RefSkip(_) if i > q => return Some(i),
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) | Cigar::Ins(l) | Cigar::SoftClip(l) => i += *l as usize,
            _ => {}
        }
    }
    None
}

/// Whether the read's own alignment carries a deletion or an insertion of
/// [`MIN_OWN_INDEL`] bases or more starting inside the genomic span `[lo, hi)`.
fn own_indel_inside(record: &Record, (lo, hi): (i64, i64)) -> bool {
    let mut rp = record.pos();
    for op in record.cigar().iter() {
        match op {
            Cigar::Del(l) | Cigar::Ins(l) if *l >= MIN_OWN_INDEL && rp >= lo && rp < hi => return true,
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) | Cigar::Del(l) | Cigar::RefSkip(l) => rp += *l as i64,
            _ => {}
        }
    }
    false
}

/// Whether one of the read's splices (N) falls between query positions `lo` and
/// `hi`: the read's bases there are not one contiguous stretch of the genome.
fn crosses_splice(record: &Record, lo: usize, hi: usize) -> bool {
    let mut q = 0usize;
    for op in record.cigar().iter() {
        match op {
            Cigar::RefSkip(_) if q > lo && q < hi => return true,
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) | Cigar::Ins(l) | Cigar::SoftClip(l) => {
                q += *l as usize;
            }
            _ => {}
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;

    /// A deletion of `len` bases after the anchor at `pos`, prepared on `ctx`
    /// (which starts at genome position 0).
    fn del(ctx: &str, pos: i64, len: usize) -> Variant {
        let p = pos as usize;
        Variant::new(
            "1".into(), pos, ctx[p..p + 1 + len].into(), ctx[p..p + 1].into(), "DEL".into(),
            None, 0, 0, None, None, Some((0, ctx.into())), None,
        )
    }

    /// `n` C/G/T bases from a fixed seed (no A, so planted A's stand out).
    fn cgt(n: usize, seed: u32) -> String {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                ['C', 'G', 'T'][(x >> 16) as usize % 3]
            })
            .collect()
    }

    #[test]
    fn only_prepared_pure_deletions_of_fifty_bases_or_more_are_judged() {
        let ctx = cgt(300, 7);
        assert!(judges(&del(&ctx, 99, 50)));
        assert!(!judges(&del(&ctx, 99, 49)));
        let mut unprepared = del(&ctx, 99, 60);
        unprepared.event_ref = None;
        assert!(!judges(&unprepared));
    }

    #[test]
    fn the_windows_hold_the_flank_and_each_allele_s_continuation() {
        // 60 bases deleted after the anchor at 99 (deleted [100, 160)), no slide.
        let mut c: Vec<u8> = cgt(300, 11).into_bytes();
        c[100] = b'A'; // the first deleted base differs from the right flank's first
        c[99] = b'A'; // and the anchor from the last deleted base: no slide either way
        c[159] = b'C';
        c[160] = b'G';
        let ctx = String::from_utf8(c).unwrap();
        let j = junction(&del(&ctx, 99, 60)).unwrap();
        assert_eq!((j.lo, j.mh), (100, 0));
        assert_eq!(j.ref_w, &ctx.as_bytes()[96..104]);
        assert_eq!(j.alt_w, [&ctx[96..100], &ctx[160..164]].concat().as_bytes());
        assert_eq!((j.first(), j.last(false), j.last(true)), (96, 103, 163));
        assert_eq!((j.genome(true, 3), j.genome(true, 4)), (99, 160));
    }

    #[test]
    fn the_microhomology_is_read_by_both_windows() {
        // The deleted bases start with the right flank's first three: the deletion
        // slides three bases right, and both windows hold them before their own K.
        let mut c: Vec<u8> = cgt(300, 13).into_bytes();
        let (lo, len) = (100, 60);
        for i in 0..3 {
            c[lo + i] = c[lo + len + i];
        }
        c[lo + 3] = if c[lo + len + 3] == b'A' { b'C' } else { b'A' };
        c[lo - 1] = if c[lo + len - 1] == b'A' { b'C' } else { b'A' };
        let ctx = String::from_utf8(c).unwrap();
        let j = junction(&del(&ctx, lo as i64 - 1, len)).unwrap();
        assert_eq!(j.mh, 3);
        assert_eq!(j.ref_w.len(), 2 * K + 3);
        assert_eq!(&j.ref_w[K..K + 3], &j.alt_w[K..K + 3]);
        assert_ne!(j.ref_w[K + 3], j.alt_w[K + 3]);
    }

    /// `n` reads of `seq` aligned at `pos` by `cigar`, each its own fragment.
    fn reads(seq: &[u8], cigar: Vec<Cigar>, pos: i64, n: usize, tag: &str) -> Vec<Record> {
        (0..n)
            .map(|i| {
                let mut r = Record::new();
                let name = format!("{tag}{i}");
                r.set(name.as_bytes(), Some(&CigarString(cigar.clone())), seq, &vec![30; seq.len()]);
                r.set_pos(pos);
                r
            })
            .collect()
    }

    /// A deletion of 60 bases after the anchor at 99 (deleted [100, 160)), no slide;
    /// four REF reads from 40, and the ALT haplotype from 40 (the junction at 60).
    fn sixty() -> (String, Variant, Junction, Vec<u8>, Vec<Record>) {
        let mut c: Vec<u8> = cgt(300, 19).into_bytes();
        c[100] = b'A';
        c[99] = b'A';
        let ctx = String::from_utf8(c).unwrap();
        let v = del(&ctx, 99, 60);
        let j = junction(&v).unwrap();
        let hap = [&ctx.as_bytes()[40..100], &ctx.as_bytes()[160..200]].concat();
        let refs = reads(&ctx.as_bytes()[40..140], vec![Cigar::Match(100)], 40, 4, "r");
        (ctx, v, j, hap, refs)
    }

    /// `seq` with a third base (A, or C where it holds A) at each index in `at`.
    fn changed(seq: &[u8], at: &[usize]) -> Vec<u8> {
        let mut s = seq.to_vec();
        for &i in at {
            s[i] = if s[i] == b'A' { b'C' } else { b'A' };
        }
        s
    }

    /// Changes at `at` and on two of every three bases for 12 past it: bases that keep
    /// mismatching after a run, as a shifted breakpoint's do.
    fn keeps(at: usize) -> Vec<usize> {
        (0..14).filter(|k| k % 3 != 2).map(|k| at + k).collect()
    }

    #[test]
    fn the_guard_masks_an_isolated_change_and_a_run_that_rejoins_the_flank() {
        let (_ctx, v, j, hap, refs) = sixty();
        let d = vec![Cigar::Match(60), Cigar::Del(60), Cigar::Match(40)];
        let guard = |at: &[usize]| {
            let all = [reads(&changed(&hap, at), d.clone(), 40, 4, "c"), refs.clone()].concat();
            guard_for(&all, &v, &|_| true, 20)
        };
        // one change two bases past the junction on every carrier: masked, named
        let g = guard(&[61]).expect("a recurring change");
        assert_eq!(g.masked, vec![j.genome(false, K + 1), j.genome(true, K + 1)]);
        assert!(g.larger.is_some());
        // one in the flank: masked once (both windows read that base)
        assert_eq!(guard(&[58]).expect("a recurring change").masked, vec![j.genome(false, 2)]);
        // two adjacent changes in the flank: a run (another allele), not masked
        assert!(guard(&[57, 58]).is_none());
        // a change on the windows' last base is masked; so is a run from it on whose
        // bases rejoin the flank after it (a phased DNP), named through its last base
        // past the window; one whose bases keep mismatching is another breakpoint
        let o = j.ref_w.len() - 1;
        let last = 60 + o - K;
        let at_last = vec![j.genome(false, o), j.genome(true, o)];
        assert!(guard(&[last]).is_some_and(|g| g.masked == at_last));
        let g = guard(&[last, last + 1]).expect("a run that rejoins the flank");
        assert_eq!(g.masked, at_last);
        let l = g.larger.expect("the deletion with the run");
        assert_eq!((l.ref_allele.len() - l.alt_allele.len(), l.carriers), (60, 4), "{l:?}");
        assert_eq!(l.pos + l.ref_allele.len() as i64, j.genome(true, o + 1) + 1, "{l:?}");
        assert!(guard(&keeps(last)).is_none());
    }

    #[test]
    fn a_run_past_the_junction_is_another_allele_only_where_the_bases_after_it_keep_mismatching() {
        // Reads fitting REF, with bases changed from the windows' last one on. Where
        // the bases after the run rejoin REF's (a germline MNP in the deleted bases)
        // it is background on the REF reads, masked; where they keep mismatching it
        // is another breakpoint, another allele, not masked. One REF change there
        // alone is a separate event, masked.
        let (ctx, v, j, _hap, _refs) = sixty();
        let o = j.ref_w.len() - 1;
        let last = (j.genome(false, o) - 40) as usize;
        let at_last = vec![j.genome(false, o), j.genome(true, o)];
        let refs = |at: &[usize]| reads(&changed(&ctx.as_bytes()[40..140], at), vec![Cigar::Match(100)], 40, 4, "r");
        let g = guard_for(&refs(&[last, last + 1]), &v, &|_| true, 20).expect("a run that rejoins REF");
        assert_eq!((g.masked, g.larger), (at_last.clone(), None));
        assert!(guard_for(&refs(&keeps(last)), &v, &|_| true, 20).is_none());
        let g = guard_for(&refs(&[last]), &v, &|_| true, 20).expect("an isolated REF change");
        assert_eq!(g.masked, at_last);
    }

    #[test]
    fn a_run_inside_the_window_is_seen_from_the_bases_past_it() {
        // A phased DNP on the first two continuation bases misses two of the four: the
        // reads are each allele's by the bases past the windows, and the run, which
        // rejoins, is masked in both windows on either allele's reads.
        let (ctx, v, j, hap, _refs) = sixty();
        let inside = vec![j.genome(false, K + 1), j.genome(true, K + 1), j.genome(false, K + 2), j.genome(true, K + 2)];
        let d = vec![Cigar::Match(60), Cigar::Del(60), Cigar::Match(40)];
        let g = guard_for(&reads(&changed(&hap, &[61, 62]), d, 40, 4, "c"), &v, &|_| true, 20).expect("carriers");
        assert_eq!(g.masked, inside);
        assert!(g.larger.is_some());
        let refs = reads(&changed(&ctx.as_bytes()[40..140], &[61, 62]), vec![Cigar::Match(100)], 40, 4, "r");
        assert_eq!(guard_for(&refs, &v, &|_| true, 20).expect("REF reads").masked, inside);
    }

    #[test]
    fn a_run_another_junction_gives_is_another_breakpoint() {
        // The right flank GGGG C, then twenty A. Carriers of a deletion one base longer
        // differ from the ALT on the window's last base and the base past it, then
        // match the A run; so do carriers of the given deletion with TT there. Only the
        // TT is the given deletion's neighbour: the longer deletion gives the other's
        // bases whole. Carriers of the deletion two bases right (a run at the junction
        // continuing the deleted bases, then the flank) are another breakpoint too.
        let mut c: Vec<u8> = cgt(300, 19).into_bytes();
        c[160..185].copy_from_slice(b"GGGGCAAAAAAAAAAAAAAAAAAAA");
        c[99] = b'A'; // the anchor differs from the last deleted base: no slide left
        c[100] = b'T'; // nor right
        c[101] = b'C';
        let ctx = String::from_utf8(c).unwrap();
        let (b, v) = (ctx.as_bytes(), del(&ctx, 99, 60));
        let j = junction(&v).unwrap();
        let o = j.ref_w.len() - 1;
        let at = |cigar: Vec<Cigar>, seq: Vec<u8>| guard_for(&reads(&seq, cigar, 40, 4, "c"), &v, &|_| true, 20);
        let longer = [&b[40..100], &b[161..201]].concat();
        assert!(at(vec![Cigar::Match(60), Cigar::Del(61), Cigar::Match(40)], longer).is_none());
        let hap = [&b[40..100], &b[160..200]].concat();
        let mut tt = hap.clone();
        tt[63..65].copy_from_slice(b"TT");
        let g = at(vec![Cigar::Match(60), Cigar::Del(60), Cigar::Match(40)], tt).expect("a phased TT");
        assert_eq!(g.masked, vec![j.genome(false, o), j.genome(true, o)]);
        let right = [&b[40..102], &b[162..200]].concat();
        assert!(at(vec![Cigar::Match(62), Cigar::Del(60), Cigar::Match(38)], right).is_none());
    }

    #[test]
    fn masks_never_leave_the_windows_identical() {
        // The continuations differ at their first two bases only; the REF haplotype
        // carries a third base at the first, the ALT haplotype at the second.
        let mut c: Vec<u8> = cgt(300, 23).into_bytes();
        c[100..104].copy_from_slice(b"CGTT");
        c[160..164].copy_from_slice(b"GCTT");
        c[99] = if c[159] == b'A' { b'C' } else { b'A' };
        let ctx = String::from_utf8(c).unwrap();
        let v = del(&ctx, 99, 60);
        let j = junction(&v).unwrap();
        let hap = [&ctx.as_bytes()[40..100], &ctx.as_bytes()[160..200]].concat();
        let refs = reads(&changed(&ctx.as_bytes()[40..140], &[60]), vec![Cigar::Match(100)], 40, 4, "r");
        let alts = reads(&changed(&hap, &[61]), vec![Cigar::Match(60), Cigar::Del(60), Cigar::Match(40)], 40, 4, "c");
        let g = guard_for(&[refs, alts].concat(), &v, &|_| true, 20).expect("both changes recur");
        let tells: Vec<i64> = (K..j.ref_w.len())
            .filter(|&o| j.ref_w[o] != j.alt_w[o])
            .flat_map(|o| [j.genome(false, o), j.genome(true, o)])
            .collect();
        assert!(tells.iter().any(|p| !g.masked.contains(p)), "masked {:?}", g.masked);
    }

    #[test]
    fn the_larger_allele_from_a_flank_change_carries_the_deletion() {
        // Carriers written as a clip with one change in the flank, on them alone.
        let (_ctx, v, _j, hap, refs) = sixty();
        let clipped = reads(&changed(&hap, &[58]), vec![Cigar::Match(60), Cigar::SoftClip(40)], 40, 4, "c");
        let g = guard_for(&[clipped, refs].concat(), &v, &|_| true, 20).expect("a recurring change");
        let l = g.larger.expect("a larger allele");
        assert_eq!(l.ref_allele.len() - l.alt_allele.len(), 60, "{l:?}");
    }

    #[test]
    fn a_reference_too_short_for_the_windows_gives_none() {
        let ctx = cgt(162, 17);
        assert!(junction(&del(&ctx, 99, 60)).is_none()); // 2 right-flank bases only
    }
}
