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
//! - **Nothing past the window is read to call a read**, for either allele: checking
//!   REF through every deleted base it covers would expose it to errors the ALT
//!   reading never meets (REJ-20261009-003). The neighbouring-change guard reads past
//!   the windows only to tell a fragment's allele and whether a run of changes
//!   rejoins its allele's continuation ([`guard_for`]).
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
//! - **Another unit count beside a run** (#297, RJ-24 amended 2026-10-30). Where the bases
//!   past the deletion repeat with the period of another allele at its breakpoint (one more
//!   or one fewer unit of a run whose deleted copy differs from it, so the windows stop at
//!   the copy's differing base), that allele reads like the ALT through the windows and on
//!   to its first difference ([`Hidden`]). Where it is present at the site (listed, or
//!   spelled in full on [`GUARD_MIN_READS`] fragments and [`GUARD_REF_SHARE`] of those
//!   spelling the ALT there, as RJ-25 finds another allele: [`is_found`]), an ALT call stands
//!   only where the read's confident bases where the two differ tell the ALT ([`tells`]);
//!   otherwise it is the ALT or that allele: neither, with partial evidence (RJ-25). A listed
//!   count is held wherever it first differs; a read that spells an unlisted one in full
//!   ([`spells_in_full`]) is never ALT, whatever the floor (RJ-6). The windows do not grow:
//!   REF, told from both at the junction, and the sibling claim keep theirs. Past what a read
//!   can spell an unlisted count nothing is held, and the row is named (UNIT_COUNT_UNSPANNED).

use log::trace;
use rust_htslib::bam::record::{Aux, Cigar, Record};

use crate::shared::bam_utils::fragment_query_span;
use crate::types::Variant;

use super::carrier::{
    aligned_query_range, confident_weight, one_base_evidence, readable_range, score, Guard, OtherAllele, Placement,
    is_found, one_error_forgiven, ReadRules, Reading, GUARD_MIN_READS, GUARD_MIN_SHARE, GUARD_REF_SHARE, SHIFT_REACH, WILD,
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

/// Reference bases prep holds past a deletion judged here, beyond its usual margin: as
/// many as it deletes, so the guard can test a deletion up to twice its length as
/// another breakpoint (one more unit of a repeat the deletion removes whole units of,
/// #285). Zero for any other variant.
pub(crate) fn reference_past(v: &Variant) -> i64 {
    let deleted = v.ref_allele.len().saturating_sub(v.alt_allele.len());
    let judged = allele_kind(&v.ref_allele, &v.alt_allele) == Some(AlleleKind::Deletion) && deleted >= LONG_DELETION;
    if judged {
        deleted as i64
    } else {
        0
    }
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

/// Another allele at the breakpoint past what a read can spell: its length change (deleted
/// bases negative), the bases a read must hold from the windows' first base to spell it in
/// full ([`Hidden::need`]), and the longest read at the site (UNIT_COUNT_UNSPANNED).
pub(crate) type Unspanned = (i64, usize, usize);

/// Bases past an allele's first difference from the ALT that its bases are held for: a
/// whole read's reach (2x250 and 2x300 runs included).
const HIDDEN_REACH: usize = super::carrier::WHOLE_READ_REACH as usize;

/// Another allele at a long deletion's breakpoint that reads like the ALT through the
/// junction windows and on (#297): the haplotype that, from the deletion's leftmost first
/// base on, continues where the reference does `shift` bases further (more bases deleted)
/// or nearer (fewer; before it, bases inserted again), or a listed row. Where the bases
/// past the deletion repeat with the period of that difference (a tandem run of the
/// deleted unit, the deleted copy differing from it), the two read alike until the run's
/// copies end.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Hidden {
    /// Its length change from REF (deleted bases negative), as the diagnostic names it.
    pub(crate) change: i64,
    /// Its bases from the windows' first base, held [`HIDDEN_REACH`] past `first_diff`
    /// (fewer where the reference held ends).
    pub(crate) bases: Vec<u8>,
    /// The first offset (from the windows' first base) where its bases and the ALT's
    /// differ; past the windows. A read must hold `first_diff + K` bases from there to
    /// spell it in full.
    pub(crate) first_diff: usize,
    /// The input lists it (a co-annotated row); otherwise the reads spell it.
    pub(crate) listed: bool,
    /// As a row (0-based POS of its anchor, REF, ALT; upper case), for the sequence judge.
    pub(crate) row: (i64, Vec<u8>, Vec<u8>),
    /// Where present at the site (`present_hidden`): the fragments placed here that spell it
    /// in full, those that tell the ALT from it, and those fitting the ALT's windows whose ALT
    /// call it takes (none of their reads tells the given count from it); the diagnostic names
    /// it, with the last two, where it takes one.
    pub(crate) spelled: u32,
    pub(crate) told: u32,
    pub(crate) held: u32,
}

impl Hidden {
    /// A listed count whose REF runs past the reference held: the bases that follow its last
    /// deleted base are not held, so its haplotype cannot be built and no read can be shown to
    /// tell it from the ALT. Its `bases` are empty, and its `first_diff` is the reference held
    /// (at least a whole read), past any read's reach (RJ-24).
    pub(crate) fn unresolved(&self) -> bool {
        self.listed && self.bases.is_empty()
    }

    /// Bases a read must hold from the windows' first base to spell it in full: [`PAST`] bases
    /// on from its first difference, as the guard reads an allele on past the windows before it
    /// places a fragment. The same bases the reach, the spelling test and the diagnostic use (a
    /// count that first differs sooner than that before the longest read's end is in no read's
    /// reach to spell).
    pub(crate) fn need(&self) -> usize {
        self.first_diff + PAST
    }
}

/// The ALT's bases from the windows' first base over the reference held for reads (prep's
/// `whole_ref`, else the prepared reference), and every other deletion at the deletion's
/// breakpoint (another unit count of a run beside it, more units deleted or fewer), or row
/// listed among `siblings`, whose bases read like the ALT's through the windows ([`Hidden`]),
/// nearest length first, but for those that are the ALT with a listed row applied over every
/// base a read of `reach` bases holds from the windows' first base. Each differs from REF
/// inside the windows (it reads like the ALT there), so REF is told from it at the junction.
/// An insertion at the breakpoint (the bases before it again: a tandem duplication, a run
/// expanded) is not looked for: a REF read inside the copy fits both, the general case of a
/// read ending inside a repeat (RJ-1), which the junction windows do not hold. None where
/// the deletion is not judged here.
pub(crate) fn hidden_alleles(v: &Variant, siblings: &[Variant], reach: usize) -> Option<(Vec<u8>, Vec<Hidden>)> {
    let j = junction(v)?;
    let (ws, seq) = v.whole_ref.as_ref().or(v.event_ref.as_ref())?;
    let w: Vec<u8> = seq.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let lo = usize::try_from(j.lo - ws).ok()?;
    let g0 = lo.checked_sub(K)?;
    let n = j.ref_w.len();
    let len = j.len;
    if lo + len > w.len() {
        return None;
    }
    // A haplotype that leaves the reference before the deletion's first base and, from it
    // on, reads it from `from` (the ALT: `lo + len`), `upto` bases from the windows' first.
    let from_shift = |from: usize, upto: usize| -> Vec<u8> {
        w[g0..lo].iter().chain(&w[from..]).take(upto).copied().collect()
    };
    let alt = from_shift(lo + len, usize::MAX);
    let telling = |b: u8| matches!(b, b'A' | b'C' | b'G' | b'T');
    let first_diff = |other: &[u8]| {
        (0..).find(|&k| match (alt.get(k), other.get(k)) {
            (Some(&a), Some(&b)) => a != b && telling(a) && telling(b),
            _ => true,
        })
    };
    let mut out: Vec<Hidden> = Vec::new();
    // Another number of bases deleted at the breakpoint: fewer than the given, or more, up to
    // twice it (as RJ-24's other breakpoints).
    for from in lo + 1..=(lo + 2 * len).min(w.len()) {
        if from == lo + len {
            continue;
        }
        // how far it reads like the ALT past the junction (most shifts differ at once); a
        // reference N tells nothing, so the two read on past it; where that runs to the end of
        // the reference held, its first difference lies further still, past any read's reach
        let agree = (0..)
            .take_while(|&i| {
                lo + len + i < w.len()
                    && from + i < w.len()
                    && (w[lo + len + i] == w[from + i] || !telling(w[lo + len + i]) || !telling(w[from + i]))
            })
            .count();
        let d = K + agree;
        if d < n {
            continue; // told inside the windows
        }
        let change = lo as i64 - from as i64;
        if d + PAST > reach {
            // Past what the longest read can spell: nothing is held for it, and only its length
            // change and where it first differs are used (the diagnostic names the nearest).
            // Its bases are kept only up to the read's end, and only where the listed-neighbour
            // test reads them (`spells_neighbour`: its first difference [`K`] or more before that
            // end); its row waits until a listed row is found to be it. Building every shift's
            // bases and row at once held memory growing with the square of the deletion (a long
            // deletion beside a long simple repeat).
            let bases = if d + K < reach { from_shift(from, reach) } else { Vec::new() };
            out.push(Hidden { change, bases, first_diff: d, listed: false, row: (*ws + lo as i64 - 1, Vec::new(), Vec::new()), spelled: 0, told: 0, held: 0 });
            continue;
        }
        let other = from_shift(from, d + HIDDEN_REACH);
        out.push(Hidden { change, bases: other, first_diff: d, listed: false, row: row_at(&w, *ws, lo, from), spelled: 0, told: 0, held: 0 });
    }
    // The listed rows: each pure indel's haplotype over the same bases, read from the
    // windows' first base where it keeps the bases before them.
    for sib in siblings {
        if !super::window::is_pure_indel(&sib.ref_allele, &sib.alt_allele) {
            continue;
        }
        let (r, a) = (sib.ref_allele.to_ascii_uppercase().into_bytes(), sib.alt_allele.to_ascii_uppercase().into_bytes());
        // its REF may run past the reference held (a longer count of the run); the part held must match
        let inside = |p: usize| (w.len() - p).min(r.len());
        let Some(p) = usize::try_from(sib.pos - ws).ok().filter(|&p| p < w.len() && w[p..p + inside(p)] == r[..inside(p)]) else {
            continue;
        };
        if p + 1 < g0 {
            continue; // it changes bases before the windows: their flank tells it
        }
        if p + r.len() > w.len() {
            // The bases after its last deleted base are not held, so where it first differs from
            // the ALT cannot be read. Where the bases held already tell it inside the windows, a
            // read of the ALT is told from it there, as for any allele. Otherwise it stays a
            // candidate (the input is taken as correct, rule X) that no read can be shown to tell:
            // past reach, held, and named (RJ-24).
            let head: Vec<u8> = w[..p].iter().chain(&a).copied().collect();
            let told_inside = head.get(g0..).is_some_and(|known| {
                (0..known.len().min(alt.len()).min(n)).any(|k| known[k] != alt[k] && telling(known[k]) && telling(alt[k]))
            });
            if !told_inside {
                let change = a.len() as i64 - r.len() as i64;
                out.push(Hidden {
                    change,
                    bases: Vec::new(),
                    first_diff: alt.len().max(reach).max(n),
                    listed: true,
                    row: (sib.pos, r, a),
                    spelled: 0,
                    told: 0,
                    held: 0,
                });
            }
            continue;
        }
        let hap: Vec<u8> = w[..p].iter().chain(&a).chain(&w[p + r.len()..]).copied().collect();
        let Some(other) = hap.get(g0..) else { continue };
        let Some(d) = first_diff(other).filter(|&d| d >= n) else { continue };
        let change = a.len() as i64 - r.len() as i64;
        let bases = other[..(d + HIDDEN_REACH).min(other.len())].to_vec();
        // a generated allele past reach (kept small, above) is the haplotype when it first differs
        // here and its bases over the same stretch are
        let from_of = |h: &Hidden| usize::try_from(lo as i64 - h.change).ok();
        let past_reach = |h: &Hidden| !h.listed && h.need() > reach;
        let same = |h: &Hidden| {
            h.bases == bases
                || (h.change == change && h.first_diff == d)
                || (past_reach(h) && h.first_diff == d && from_of(h).is_some_and(|from| from_shift(from, d + HIDDEN_REACH) == bases))
        };
        match out.iter_mut().find(|h| same(h)) {
            Some(h) => {
                if past_reach(h) {
                    if let Some(from) = from_of(h) {
                        h.bases = from_shift(from, h.first_diff + HIDDEN_REACH);
                        h.row = row_at(&w, *ws, lo, from);
                    }
                }
                h.listed = true;
            }
            None => out.push(Hidden { change, bases, first_diff: d, listed: true, row: (sib.pos, r, a), spelled: 0, told: 0, held: 0 }),
        }
    }
    // An allele whose bases, past where it and the ALT first differ, are the ALT's with a
    // listed row applied past the junction (a phased neighbour the input names: a 1-bp
    // insertion that shifts the carriers' bases past a run onto one more unit deleted) is the
    // given deletion with that neighbour: the input is taken as correct (RJ-24), so it is
    // never another allele. Compared over the bases a read of `reach` bases holds: past them
    // no read can tell the two apart.
    let neighbours = with_listed_rows(v, &j, &alt);
    let spells_neighbour = |h: &Hidden| {
        neighbours.iter().any(|x| {
            let to = h.bases.len().min(x.len()).min(reach);
            to > h.first_diff + K && h.bases[h.first_diff..to] == x[h.first_diff..to]
        })
    };
    out.retain(|h| h.listed || !spells_neighbour(h));
    let given = -(len as i64);
    out.sort_by_key(|h| ((h.change - given).abs(), h.change));
    let mut alt = alt;
    let far = out.iter().filter(|h| !h.unresolved()).map(|h| h.first_diff + HIDDEN_REACH).max().unwrap_or(0);
    alt.truncate(far.max(n));
    // A substitution another input row lists on the ALT's bases is that row's allele, a
    // neighbour the carriers may hold (RJ-24): its offset tells no allele from another, so it
    // is read as no base by either (the offsets are the ALT's, past the junction, which every
    // such allele shares up to its first difference).
    for (p, r, _) in &v.listed_neighbours {
        for i in 0..r.len() as i64 {
            let g = p + i;
            let k = if g < j.lo { g - j.first() } else { g - (j.lo + len as i64) + K as i64 };
            if let Ok(k) = usize::try_from(k) {
                if k >= n {
                    for bases in std::iter::once(&mut alt).chain(out.iter_mut().map(|h| &mut h.bases)) {
                        if let Some(b) = bases.get_mut(k) {
                            *b = WILD;
                        }
                    }
                }
            }
        }
    }
    Some((alt, out))
}

/// The diagnostic's unit counts among `hidden` (nearest length first): on each side of the
/// given deletion (`given`, its length change), more bases deleted and fewer, the nearest
/// allele that reads like the ALT past the junction for the bases it shifts the reference by
/// and [`REJOIN`] more (a copy of that length at least twice over: a unit count of a run, not
/// a shift a few bases read alike by chance, nor a far count a repeat whose copies differ
/// lets read alike).
fn nearest_counts(hidden: &[Hidden], given: i64) -> Vec<&Hidden> {
    let unit = |h: &Hidden| (h.change - given).unsigned_abs() as usize;
    // (an unresolved count is named wherever it is nearest: nothing says it reads unlike the ALT)
    let count = |h: &&Hidden| h.unresolved() || h.first_diff.saturating_sub(K) >= unit(h) + REJOIN;
    let more = hidden.iter().filter(|h| h.change < given).filter(count).min_by_key(|h| unit(h));
    let fewer = hidden.iter().filter(|h| h.change > given).filter(count).min_by_key(|h| unit(h));
    let mut out: Vec<&Hidden> = more.into_iter().chain(fewer).collect();
    out.sort_by_key(|h| (unit(h), h.change));
    out
}

/// The ALT's bases from the windows' first base (`alt`) with each listed row past the junction
/// applied, one at a time (`Variant::listed_neighbours` and `listed_indels`): the given
/// deletion with a neighbour the input names.
pub(crate) fn with_listed_rows(v: &Variant, j: &Junction, alt: &[u8]) -> Vec<Vec<u8>> {
    v.listed_neighbours
        .iter()
        .chain(&v.listed_indels)
        .filter_map(|(p, r, a)| {
            let k = usize::try_from(p - (j.lo + j.len as i64) + K as i64).ok().filter(|&k| k >= K)?;
            let (r, a) = (r.to_ascii_uppercase().into_bytes(), a.to_ascii_uppercase().into_bytes());
            let held = alt.get(k..k + r.len())?;
            let fits = held.iter().zip(&r).all(|(&x, &y)| x == y || x == WILD);
            fits.then(|| [&alt[..k], &a[..], &alt[k + r.len()..]].concat())
        })
        .collect()
}

/// The read's confident mismatches against `hap` (bases from the windows' first base, read
/// from query position `s`) over offsets `from` up to query position `end`; None where it reads
/// none of them.
fn misses_over(seq: &[u8], quals: &[u8], min_baseq: u8, s: usize, hap: &[u8], from: usize, end: usize) -> Option<usize> {
    let (mut seen, mut misses) = (0usize, 0usize);
    for (k, &y) in hap.iter().enumerate().skip(from) {
        let q = s + k;
        if q >= end || q >= seq.len() || q >= quals.len() {
            break;
        }
        let (base, bq) = (seq[q].to_ascii_uppercase(), quals[q]);
        if base == WILD || bq < min_baseq || !matches!(y, b'A' | b'C' | b'G' | b'T') {
            continue;
        }
        seen += 1;
        misses += usize::from(base != y);
    }
    (seen > 0).then_some(misses)
}

/// Whether the read's confident bases show a listed row's change (`x`: the ALT's bases `alt`
/// with that row applied, [`with_listed_rows`]): at one base or more of offsets `from` up to
/// query position `end` where `x` and the ALT differ, it reads `x`'s base. A row whose change
/// no base the read holds shows (it lies past the read, or before `from` where the read cannot
/// tell), or whose base the read reads as the ALT's or not confidently, is not carried by it.
#[allow(clippy::too_many_arguments)]
fn shows_change(seq: &[u8], quals: &[u8], min_baseq: u8, s: usize, alt: &[u8], x: &[u8], from: usize, end: usize) -> bool {
    (from..x.len()).take_while(|&k| s + k < end && s + k < seq.len() && s + k < quals.len()).any(|k| {
        let (base, bq) = (seq[s + k].to_ascii_uppercase(), quals[s + k]);
        matches!(x[k], b'A' | b'C' | b'G' | b'T') && alt.get(k).is_some_and(|&y| y != x[k]) && base != WILD && bq >= min_baseq && base == x[k]
    })
}

/// Whether the read is the given deletion with a listed neighbour (`neighbours`, from
/// [`with_listed_rows`], over the ALT's bases `alt`): its confident bases show that row's change
/// ([`shows_change`]) and fit the deletion with it over offsets `from` up to query position `end`
/// as well as allele `b`, and but for one base. It is that, the input taken as correct
/// (RJ-24), never `b`; a read that shows no listed row's change is not (RJ-28: a row's count
/// does not depend on an input row the read does not carry).
#[allow(clippy::too_many_arguments)]
fn a_listed_neighbour(
    seq: &[u8],
    quals: &[u8],
    min_baseq: u8,
    s: usize,
    b: &[u8],
    (alt, neighbours): (&[u8], &[Vec<u8>]),
    from: usize,
    end: usize,
) -> bool {
    let own = misses_over(seq, quals, min_baseq, s, b, from, end).unwrap_or(0);
    neighbours.iter().any(|x| {
        misses_over(seq, quals, min_baseq, s, x, from, end).is_some_and(|m| m <= 1 && m <= own)
            && shows_change(seq, quals, min_baseq, s, alt, x, from, end)
    })
}

/// The row of the haplotype that, from the deletion's first base `lo` (offset in `w`,
/// which starts at genome position `ws`), reads the reference from `from` (past `lo`): a
/// deletion of `from - lo` bases after the anchor.
fn row_at(w: &[u8], ws: i64, lo: usize, from: usize) -> (i64, Vec<u8>, Vec<u8>) {
    (ws + lo as i64 - 1, w[lo - 1..from].to_vec(), vec![w[lo - 1]])
}

/// Whether the aligner places the read here: not elsewhere at least as well (XS below AS),
/// and, where it wrote no XS, at MAPQ above 0 (RJ-29's placement test). A read it doubts may
/// belong to another copy, so it spells no allele here.
fn placed_here(record: &Record) -> bool {
    let tag = |t: &[u8]| match record.aux(t) {
        Ok(Aux::I8(v)) => Some(i64::from(v)),
        Ok(Aux::U8(v)) => Some(i64::from(v)),
        Ok(Aux::I16(v)) => Some(i64::from(v)),
        Ok(Aux::U16(v)) => Some(i64::from(v)),
        Ok(Aux::I32(v)) => Some(i64::from(v)),
        Ok(Aux::U32(v)) => Some(i64::from(v)),
        _ => None,
    };
    match (tag(b"AS"), tag(b"XS")) {
        (Some(a), Some(x)) => x < a,
        _ => record.mapq() > 0,
    }
}

/// Which of two alleles a read's confident bases tell (true: `a`), its bases read from query
/// position `s` against each allele's bases from the windows' first base, from offset
/// `from` up to query position `end` (exclusive): only the bases where the two differ weigh.
/// A base fitting neither costs both alike, an error (RJ-26); a base below `--min-baseq` or
/// N weighs nothing. An allele is told where the read shows the other's base nowhere (or at
/// one base where it shows this one's at three or more, a one-off error, RJ-27) and its
/// bases favour it by one base read at `--min-baseq` (RJ-16). None where neither is told: the
/// read stops short of their first difference, or its bases there are unread.
#[allow(clippy::too_many_arguments)]
pub(crate) fn tells(seq: &[u8], quals: &[u8], min_baseq: u8, s: i64, (a, b): (&[u8], &[u8]), from: usize, end: usize) -> Option<bool> {
    let (mut miss_a, mut miss_b, mut wa, mut wb) = (0usize, 0usize, 0.0f64, 0.0f64);
    for k in from..a.len().min(b.len()) {
        let Some(q) = usize::try_from(s + k as i64).ok() else { continue };
        if q >= end || q >= seq.len() || q >= quals.len() {
            break;
        }
        let (x, y) = (a[k], b[k]);
        if x == y || !matches!(x, b'A' | b'C' | b'G' | b'T') || !matches!(y, b'A' | b'C' | b'G' | b'T') {
            continue;
        }
        let (base, bq) = (seq[q].to_ascii_uppercase(), quals[q]);
        if base == WILD || bq < min_baseq {
            continue;
        }
        if base == x {
            miss_b += 1;
            wb += one_base_evidence(bq);
        } else if base == y {
            miss_a += 1;
            wa += one_base_evidence(bq);
        }
    }
    let need = one_base_evidence(min_baseq);
    let favours = |own: usize, other: usize, w_own: f64, w_other: f64| {
        other >= 1 && (own == 0 || one_error_forgiven(own, other)) && w_other - w_own >= need && w_other - w_own > 0.0
    };
    if favours(miss_a, miss_b, wa, wb) {
        Some(true)
    } else if favours(miss_b, miss_a, wb, wa) {
        Some(false)
    } else {
        None
    }
}

/// Whether a read placed from query position `s` spells allele `b` in full against the ALT
/// `a` and REF `r` (each allele's bases from the windows' first base), from offset `from` (the
/// windows' end) up to query position `end`: at most one of its confident bases there
/// mismatches `b`, none where `a` and `b` differ shows `a`'s base, two or more there show
/// `b`'s, and REF read the same way misses two more (a REF molecule of the sample's own
/// repeat copies, read past windows its changes let fit the ALT, is no such allele: RJ-25's
/// found allele is no worse than REF).
fn spells(seq: &[u8], quals: &[u8], min_baseq: u8, s: usize, (a, b, r): (&[u8], &[u8], &[u8]), from: usize, end: usize) -> bool {
    let acgt = |x: u8| matches!(x, b'A' | b'C' | b'G' | b'T');
    let (mut misses, mut shown, mut ref_misses) = (0, 0, 0);
    for (k, &y) in b.iter().enumerate().skip(from) {
        let q = s + k;
        if q >= end || q >= seq.len() || q >= quals.len() {
            break;
        }
        let (base, bq) = (seq[q].to_ascii_uppercase(), quals[q]);
        if base == WILD || bq < min_baseq || !acgt(y) {
            continue;
        }
        ref_misses += usize::from(r.get(k).is_some_and(|&x| acgt(x) && x != base));
        let differ = a.get(k).is_some_and(|&x| x != y && acgt(x));
        if base == y {
            shown += usize::from(differ);
            continue;
        }
        if differ && a.get(k) == Some(&base) {
            return false;
        }
        misses += 1;
        if misses > 1 {
            return false;
        }
    }
    shown >= 2 && ref_misses >= misses + 2
}

/// Whether a read placed from query position `q` and read up to `end` spells allele `h` in
/// full (`told`: what [`tells`] found for it): it holds [`PAST`] bases on from the allele's first
/// difference from the ALT, its bases where the two differ show the allele's and not the ALT's,
/// past the windows they fit the allele but for one ([`spells`]; `bases`: the ALT's, REF's and
/// the ALT with each listed row applied, from the windows' first base), and it is not the
/// deletion with a listed neighbour. The one test of a read that shows another unit count:
/// [`present_hidden`] counts the fragments passing it, and [`classify`] holds the read.
#[allow(clippy::too_many_arguments)]
fn spells_in_full(
    seq: &[u8],
    quals: &[u8],
    min_baseq: u8,
    (q, end): (usize, usize),
    told: Option<bool>,
    (alt, reference, neighbours): (&[u8], &[u8], &[Vec<u8>]),
    n: usize,
    h: &Hidden,
) -> bool {
    told == Some(false)
        && q + h.need() <= end
        && spells(seq, quals, min_baseq, q, (alt, &h.bases, reference), n, end)
        && !a_listed_neighbour(seq, quals, min_baseq, q, &h.bases, (alt, neighbours), n, end)
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
    let guard = rules.guard.and_then(|g| g.get(v));
    let masked = guard.map_or(&[][..], |g| g.masked.as_slice());
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
    if out.is_alt {
        // Another allele at the breakpoint present at the site that reads like the ALT past
        // the windows (#297): the call stands only where the read's bases where the two
        // differ tell the ALT; else it is the ALT or that allele, partial evidence (RJ-25).
        if let Some(g) = guard.filter(|g| !g.unit_others.is_empty() || !g.unit_absent.is_empty()) {
            let end = splice_after(record, hi).unwrap_or(seq.len()).min(readable.1);
            let untold = g.unit_others.iter().find(|h| {
                tells(&seq, quals, min_baseq, lo as i64, (&g.unit_alt, &h.bases), h.first_diff, end) != Some(true)
                    && !a_listed_neighbour(&seq, quals, min_baseq, lo, &h.bases, (&g.unit_alt, &g.unit_neighbours), h.first_diff, end)
            });
            if let Some(h) = untold {
                trace!(
                    "junction: an ALT read not told from another allele at the breakpoint ({:+} bp, first differing {} bases from the windows' first base) → another allele",
                    h.change, h.first_diff,
                );
                out = ClassifyResult::neither_with_nearby(qual, ClassifyPhase::MaskedCompare);
            } else if let Some(h) = g.unit_absent.iter().find(|h| {
                let told = tells(&seq, quals, min_baseq, lo as i64, (&g.unit_alt, &h.bases), h.first_diff, end);
                spells_in_full(&seq, quals, min_baseq, (lo, end), told, (&g.unit_alt, &g.unit_ref, &g.unit_neighbours), n, h)
            }) {
                // A read whose own bases spell another unit count is no ALT, whatever the
                // presence floor: the floor decides only the reads that cannot tell (RJ-6).
                // Partial evidence where the aligner places it here alone, else depth only
                // (RJ-29: a read it doubts may belong to another copy). The side the window
                // was placed from decides nothing: a read placed from the right flank is ALT
                // where it fits, else another allele (RJ-24).
                let here = placed_here(record);
                trace!(
                    "junction: an ALT read that spells another allele at the breakpoint ({:+} bp) → {}",
                    h.change, if here { "partial evidence" } else { "depth only" },
                );
                out = if here {
                    ClassifyResult::neither_with_nearby(qual, ClassifyPhase::MaskedCompare)
                } else {
                    let mut depth = ClassifyResult::neither(ClassifyPhase::MaskedCompare);
                    depth.uninformative = true;
                    depth
                };
            }
        }
    }
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

/// A fragment's base at a guard offset: the base, its quality, and whether its mates
/// disagree there at equal quality (then it is read, but not confidently).
type FragmentBase = Option<(u8, u8, bool)>;

/// Bases past the windows that tell a read's allele (the [`K`] bases past the
/// microhomology may hold a run, whose bases may be the other allele's), and that the
/// bases after a run must rejoin its allele's continuation over.
const REJOIN: usize = 12;
/// Bases past the windows the guard reads: a run running on past them, and the
/// [`REJOIN`] after it.
const PAST: usize = 2 * REJOIN;

/// How the guard reads a record (F6): its base qualities as the classifiers read them (heuristic
/// BAQ where it applies, `--apply-baq`) and the query range they may read (in RNA a soft clip
/// reaching an exon edge is left out, [`readable_range`]); None where the record shows nothing
/// here. The presence test of another unit count reads the reads as the hold does: one quality
/// contract (invariant 2).
pub(crate) type ReadAs<'a> = &'a dyn for<'r> Fn(&'r Record) -> Readable<'r>;

/// A record's qualities and the query range they may read ([`ReadAs`]).
pub(crate) type Readable<'r> = Option<(std::borrow::Cow<'r, [u8]>, (usize, usize))>;

/// A record as stored: its qualities, every base readable (DNA, no BAQ).
#[cfg(test)]
fn as_stored(record: &Record) -> Readable<'_> {
    Some((std::borrow::Cow::Borrowed(record.qual()), (0, record.seq_len())))
}

/// [`guard_read_as`] with the records as stored (tests).
#[cfg(test)]
pub(crate) fn guard_for(
    read_cache: &[Record],
    v: &Variant,
    counted: &dyn Fn(&Record) -> bool,
    min_baseq: u8,
    siblings: &[Variant],
) -> Option<Guard> {
    guard_read_as(read_cache, v, counted, min_baseq, siblings, &as_stored)
}

/// The neighbouring-change guard for a long deletion (the carrier rule's, RJ-23,
/// over the junction windows): each fragment of the reads `counted` admits that holds
/// the windows from a read's aligned first window base (its mates' bases merged, the
/// better-quality one's where both read a base) is taken as the allele whose
/// continuation it follows (as its bases past the windows place it, [`placed_past`];
/// else as its [`K`] bases past the microhomology tell it: they fit that allele to
/// within one base and the other by two or more), and its confident changes are
/// tallied per
/// window offset, on the base before the windows, and on [`PAST`] bases on along its
/// allele past them (up to a splice). A change recurs on at least [`GUARD_MIN_READS`]
/// fragments and [`GUARD_MIN_SHARE`] of that allele's reads. Masked at its offset in both windows are, in the flank and
/// microhomology both alleles share, the changes the REF reads carry (a germline SNP
/// or MNP, background on both haplotypes, runs included, on [`GUARD_MIN_READS`] REF
/// reads or more); each isolated change of the ALT reads' own (no own change at
/// either neighbouring base); past the junction each isolated change on the REF
/// reads; and each run past the junction, on either allele's reads, whose bases rejoin
/// the allele's continuation after it ([`rejoins`]): a phased neighbour. Any other run
/// of the ALT reads' own changes is an indel at the junction, and any other run past it
/// on reads fitting REF another breakpoint: another allele, not masked. An allele's own
/// changes are taken together first: where another breakpoint gives the reads' bases
/// exactly ([`another_breakpoint`]), none of them is masked, isolated or run. A change
/// on the ALT reads that another input row lists (`Variant::listed_neighbours`) is that
/// row's allele, not one of their own: masked, never another breakpoint's, not named. Masks never take every
/// offset where REF and ALT differ: those stay read, and a read changed there is
/// neither allele. On the ALT reads alone (under [`GUARD_REF_SHARE`] of REF's)
/// isolated recurring changes and phased runs also make the larger allele the
/// diagnostic names (the deletion with the changes, through a run's last base past
/// the windows), and beside any other run none is named. None when nothing is masked
/// or named, or the variant is not judged here.
pub(crate) fn guard_read_as(
    read_cache: &[Record],
    v: &Variant,
    counted: &dyn Fn(&Record) -> bool,
    min_baseq: u8,
    siblings: &[Variant],
    read_as: ReadAs,
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
    // Each fragment's bases by edged offset, from the reads holding the windows from
    // their aligned first window base. Where both mates read a base, the
    // better-quality one's stands.
    let mut order: Vec<&[u8]> = Vec::new();
    let mut fragments: std::collections::HashMap<&[u8], Vec<FragmentBase>> = std::collections::HashMap::new();
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
        let merged = fragments.entry(record.qname()).or_insert_with(|| {
            order.push(record.qname());
            vec![None; m]
        });
        for (e, slot) in merged.iter_mut().enumerate() {
            let Some(i) = (q + e).checked_sub(1).filter(|&i| e <= n || i < stop) else { continue };
            let (Some(&b), Some(&bq)) = (bases.get(i), quals.get(i)) else { continue };
            let b = b.to_ascii_uppercase();
            *slot = match *slot {
                Some((x, xq, tie)) if x == b => Some((x, xq.max(bq), tie && bq <= xq)),
                Some((x, xq, tie)) if xq >= bq => Some((x, xq, tie || xq == bq)),
                _ => Some((b, bq, false)),
            };
        }
    }
    let mut scanned = 0u32;
    for name in order {
        let fragment = &fragments[name];
        // The fragment's base at edged offset e, and whether it is read confidently.
        let read = |e: usize| fragment[e].map(|(b, bq, tie)| (b, !tie && b != WILD && bq >= min_baseq));
        let off = |w: &[Option<u8>], e: usize| matches!((read(e), w[e]), (Some((b, true)), Some(x)) if b != x);
        // The allele whose continuation it follows: as the bases past the windows place
        // it, else as its K bases past the microhomology tell it (a change among those
        // may be the other allele's base there).
        let told = |es: std::ops::RangeInclusive<usize>| {
            let misses = |w: &[Option<u8>]| es.clone().filter(|&e| off(w, e)).count();
            match (misses(&windows[0]), misses(&windows[1])) {
                (0..=1, 2..) => Some(0),
                (2.., 0..=1) => Some(1),
                _ => None,
            }
        };
        let Some(allele) = placed_past(&read, &windows, n).or_else(|| told(cont + 1..=n)) else {
            continue;
        };
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
    let most = |a: usize, e: usize| b"ACGT"[(0..4).max_by_key(|&k| shown[a][e][k]).unwrap_or(0)];
    // A change on the ALT reads that another input row lists (a substitution giving that
    // base there) is that row's allele: the input is taken as correct, so the reads are
    // the deletion with a listed neighbour. It is masked, never evidence of another
    // breakpoint, and not named (#285).
    let listed = |e: usize| {
        e >= 1 && recurs(1, e) && {
            let (g, b) = (j.genome(true, e - 1), most(1, e));
            v.listed_neighbours.iter().any(|(p, r, a)| {
                let i = usize::try_from(g - p).ok();
                i.and_then(|i| Some((r.as_bytes().get(i)?, a.as_bytes().get(i)?)))
                    .is_some_and(|(&x, &y)| y.to_ascii_uppercase() == b && x.to_ascii_uppercase() != b)
            })
        }
    };
    // The ALT reads' own changes, isolated when no own change sits beside them; past
    // the junction a REF reads' change, when none recurs beside it on the REF reads.
    let own = |e: usize| recurs(1, e) && !background(e) && !listed(e);
    let isolated = |e: usize| !own(e - 1) && !own(e + 1);
    let ref_isolated = |e: usize| recurs(0, e) && !recurs(0, e - 1) && !recurs(0, e + 1);
    // Each allele's own changes: the ALT reads' (not background, not listed), and past
    // the junction the REF reads'. Those in the windows, with any running on past them,
    // taken together are another breakpoint's where one gives the reads' bases exactly:
    // then none is masked, and the reads carrying them are another allele (#296).
    let own_of = |a: usize, e: usize| if a == 1 { own(e) } else { e > cont && recurs(0, e) };
    let elsewhere: [bool; 2] = [0, 1].map(|a| {
        // The changes in the windows, and on past them each within REJOIN bases of the
        // last one taken.
        let mut changed: Vec<usize> = (1..=n).filter(|&e| own_of(a, e)).collect();
        while let Some(&last) = changed.last() {
            match (last + 1..(last + 1 + REJOIN).min(m)).find(|&e| own_of(a, e)) {
                Some(e) => changed.push(e),
                None => break,
            }
        }
        let compared = |e: usize| e < m && cover[a][e] >= GUARD_MIN_READS && (changed.contains(&e) || !recurs(a, e));
        another_breakpoint(&j, &base, a == 1, &changed, &compared, &|e| most(a, e))
    });
    // Each allele's runs reaching the windows (edged [first, last], two or more own
    // changes side by side). Those whose bases rejoin the allele's continuation are
    // phased neighbours, masked; any other is another allele.
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
    let runs = [runs_of(0), runs_of(1)];
    let neighbours: [Vec<(usize, usize)>; 2] = [0, 1].map(|a| {
        let phased = |&(e0, e1): &(usize, usize)| rejoins(&j, (e0, e1), m, &|e| cover[a][e], &|e| recurs(a, e));
        runs[a].iter().copied().filter(|r| !elsewhere[a] && phased(r)).collect()
    });
    let in_neighbour = |a: usize, e: usize| neighbours[a].iter().any(|&(e0, e1)| (e0..=e1).contains(&e));
    let mut offsets: Vec<usize> = (1..=n)
        .filter(|&e| {
            background(e)
                || listed(e)
                || (own(e) && isolated(e) && !elsewhere[1])
                || (e > cont && ref_isolated(e) && !elsewhere[0])
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
        (1..m).filter(|&e| recurs(1, e) && !on_ref(e) && !listed(e) && (e <= n || in_neighbour(1, e))).collect();
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
    // run that is no phased neighbour, or changes another breakpoint gives, none: the
    // tallies cannot spell another allele (an MNP at the junction, or an indel shifting
    // the bases past it); the observed allele names it from the reads' alignments.
    let run = runs[1].len() > neighbours[1].len() || elsewhere[1];
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
    let recurring = |a: usize| {
        (0..m).filter(|&e| recurs(a, e)).map(|e| format!("{}:{}/{}", e, change[a][e], cover[a][e])).collect::<Vec<_>>().join(" ")
    };
    trace!(
        "junction guard tallies at {}:{}: {} fragments; recurring changes (edged offset: changed/read) REF [{}] ALT [{}]; another breakpoint REF/ALT {}/{}",
        v.chrom, v.pos + 1, scanned, recurring(0), recurring(1), elsewhere[0], elsewhere[1],
    );
    // Other alleles at the breakpoint that read like the ALT past the windows (#297): those
    // present at the site within a read's reach hold an ALT call to their first difference;
    // one past the longest read at the site names the row.
    let found = present_hidden(read_cache, v, counted, min_baseq, &j, &masked, siblings, read_as);
    if masked.is_empty() && larger.is_none() && found.others.is_empty() && found.absent.is_empty() && found.unspanned.is_none() {
        return None;
    }
    trace!(
        "junction guard at {}:{} {}>{}: {} window offsets masked, {} on the ALT reads alone, {} phased runs, another breakpoint on REF/ALT reads {}/{} ({} fragments read)",
        v.chrom, v.pos + 1, v.ref_allele, v.alt_allele, offsets.len(), alt_only.len(),
        neighbours[0].len() + neighbours[1].len(), elsewhere[0], elsewhere[1], scanned,
    );
    let mut guard = Guard::new(v, masked, larger, scanned);
    // the found ones are site states for the sequence judge (the listed ones are its rows)
    guard.judge_others = found
        .others
        .iter()
        .filter(|h| !h.listed)
        .map(|h| OtherAllele { pos: h.row.0, ref_allele: h.row.1.clone(), alt_allele: h.row.2.clone() })
        .collect();
    (guard.unit_alt, guard.unit_others, guard.unit_absent, guard.unit_ref, guard.unit_neighbours, guard.unit_unspanned) =
        (found.alt, found.others, found.absent, found.reference, found.neighbours, found.unspanned);
    Some(guard)
}

/// The other alleles at a long deletion's breakpoint that read like the ALT past the
/// windows ([`hidden_alleles`], `siblings` the listed rows) and are present at the site
/// within the longest read's reach, with the ALT's bases they are told from, and the
/// nearest one past that reach (its change, the bases a read must hold, the longest read).
/// Present: listed, or spelled in full (from the windows' first base through [`K`] bases
/// from its first difference, as the windows hold [`K`] bases of each allele's
/// continuation) on [`GUARD_MIN_READS`] fragments or more and on [`GUARD_REF_SHARE`] of
/// the fragments spelling the ALT over the same bases (RJ-25's floor). Fragments are read
/// from the reads `counted` admits that are placed from the left flank (their aligned first
/// window base, as RJ-25 places them) and fit the ALT's windows (the guard's `masked`
/// positions fitting any base) but not REF's, each up to its next splice; mates telling the
/// two apart differently spell neither.
#[allow(clippy::too_many_arguments)]
fn present_hidden(
    read_cache: &[Record],
    v: &Variant,
    counted: &dyn Fn(&Record) -> bool,
    min_baseq: u8,
    j: &Junction,
    masked: &[i64],
    siblings: &[Variant],
    read_as: ReadAs,
) -> PresentHidden {
    let (first, last) = (j.first(), j.last(true));
    let at_site = |r: &&Record| counted(r) && r.pos() <= last && crate::shared::bam_utils::ref_end(r) > first;
    let longest = read_cache.iter().filter(at_site).map(Record::seq_len).max().unwrap_or(0);
    let Some((alt, hidden)) = hidden_alleles(v, siblings, longest).filter(|(_, h)| !h.is_empty()) else {
        return PresentHidden::default();
    };
    // (no read at the site reads nothing, past any reach: nothing to name)
    let unspanned = nearest_counts(&hidden, -(j.len as i64))
        .into_iter()
        .find(|h| longest > 0 && h.need() > longest)
        .map(|h| (h.change, h.need(), longest));
    // a listed count is held wherever it first differs: the reads that cannot reach the
    // difference cannot tell it from the ALT (invariant 7: the input is taken as correct)
    let reach: Vec<Hidden> = hidden.into_iter().filter(|h| h.listed || h.need() <= longest).collect();
    if reach.is_empty() {
        return PresentHidden { unspanned, ..PresentHidden::default() };
    }
    let window = |is_alt: bool, w: &[u8]| -> Vec<u8> {
        w.iter().enumerate().map(|(o, &b)| if masked.contains(&j.genome(is_alt, o)) { WILD } else { b }).collect()
    };
    let (ref_w, alt_w) = (window(false, &j.ref_w), window(true, &j.alt_w));
    let n = ref_w.len();
    let neighbours = with_listed_rows(v, j, &alt);
    // REF's bases from the windows' first base, as far as the ALT's are held
    let reference: Vec<u8> = v
        .whole_ref
        .as_ref()
        .or(v.event_ref.as_ref())
        .and_then(|(ws, seq)| {
            let g0 = usize::try_from(first - ws).ok()?;
            Some(seq.bytes().skip(g0).take(alt.len()).map(|b| b.to_ascii_uppercase()).collect())
        })
        .unwrap_or_default();
    // per fragment, per allele in reach (`Tally`)
    let mut fragments: std::collections::HashMap<&[u8], Vec<Tally>> = std::collections::HashMap::new();
    for record in read_cache.iter().filter(at_site) {
        // the qualities and the bases it may read, as `classify` reads them (F6)
        let Some((quals, readable)) = read_as(record) else { continue };
        let (seq, quals) = (record.seq().as_bytes(), &quals[..]);
        // placed from the left flank, as RJ-25 places a fragment that spells an allele; else
        // from the right one, as `classify` reads it, to tally the calls the hold takes only
        let left = find_read_pos(record, first);
        let Some(q) = left
            .or_else(|| find_read_pos(record, last).and_then(|q| q.checked_sub(n - 1)))
            .filter(|&q| q >= readable.0 && q + n <= readable.1 && q + n <= seq.len() && q + n <= quals.len())
        else {
            continue;
        };
        let (r, a) = (score(&seq[q..q + n], &quals[q..q + n], min_baseq, &ref_w, (q, q + n)), score(&seq[q..q + n], &quals[q..q + n], min_baseq, &alt_w, (q, q + n)));
        if a.mismatches > 0 || r.mismatches == 0 {
            continue;
        }
        let end = splice_after(record, q).unwrap_or(seq.len()).min(readable.1).min(quals.len());
        // a read the aligner places elsewhere as well spells no allele here (RJ-29), nor does
        // one read from the right flank
        let here = placed_here(record) && left.is_some();
        let calls = fragments.entry(record.qname()).or_insert_with(|| vec![Tally::default(); reach.len()]);
        for (h, slot) in reach.iter().zip(calls.iter_mut()) {
            let told_here = tells(&seq, quals, min_baseq, q as i64, (&alt, &h.bases), h.first_diff, end);
            // a read the hold would take from the ALT, as `classify` holds one
            slot.held |= told_here != Some(true)
                && !a_listed_neighbour(&seq, quals, min_baseq, q, &h.bases, (&alt, &neighbours), h.first_diff, end);
            // (a mate that tells the ALT is recorded here, before the length check below: ending
            // short of the allele's first difference plus 24 bases, it still disagrees with a
            // mate that spells the allele, `Tally::disagrees`)
            slot.alt |= told_here == Some(true);
            // read in full: [`PAST`] bases on from the first difference, as the guard reads an
            // allele on past the windows before it places a fragment
            if q + h.first_diff + PAST > end {
                continue;
            }
            let Some(told) = told_here else { continue };
            // Spelled in full: past the windows its bases fit the allele, but for one, wherever
            // they read (not only where it and the ALT differ), showing its base at two or more
            // of the bases where the two differ and the ALT's at none. A read of a third allele
            // (in a repeat whose copies differ, the sample's own copy changes) is not taken for
            // it.
            if !told {
                // placed here alone, else it is shown but not tallied (RJ-29): the hold still takes it
                let shown = spells_in_full(&seq, quals, min_baseq, (q, end), told_here, (&alt, &reference, &neighbours), n, h);
                slot.shown |= shown;
                if !(here && shown) {
                    continue;
                }
            }
            trace!(
                "junction: read {} {} the allele {:+} bp at the breakpoint",
                String::from_utf8_lossy(record.qname()), if told { "tells the ALT over" } else { "spells" }, h.change,
            );
            match slot.told {
                Some(x) if x != told => slot.conflict = true,
                _ => slot.told = Some(told),
            }
        }
    }
    let (mut present, mut absent): (Vec<Hidden>, Vec<Hidden>) = (Vec::new(), Vec::new());
    for (k, h) in reach.into_iter().enumerate() {
        let count = |want: bool| fragments.values().filter(|c| !c[k].disagrees() && c[k].told == Some(want)).count();
        let (spelled, alt_told) = (count(false), count(true));
        // fragments whose ALT call it would take: a read fitting the ALT's windows that
        // tells the given count from it nowhere, its mate no better
        let held = fragments.values().filter(|c| c[k].held && !c[k].alt).count();
        let found = is_found(spelled, alt_told);
        if found || h.listed {
            trace!(
                "junction: another allele at {}:{} ({:+} bp, reading like the ALT for {} bases from the windows' first base) {}: {} fragments spell it, {} the ALT over the same bases",
                v.chrom, v.pos + 1, h.change, h.first_diff, if h.listed { "listed" } else { "found" }, spelled, alt_told,
            );
            present.push(Hidden { spelled: spelled as u32, told: alt_told as u32, held: held as u32, ..h });
        } else if fragments.values().any(|c| c[k].shown) {
            // not at the site by the floor, yet a read spells it: that read is no ALT
            absent.push(h);
        }
    }
    PresentHidden { alt, others: present, absent, neighbours, reference, unspanned }
}

/// One fragment's reading against one other allele at the breakpoint ([`present_hidden`]):
/// read in full, the ALT told (`Some(true)`) or the allele spelled (`Some(false)`), mates
/// disagreeing there (`conflict`); and, wherever it reads, a read of it the hold takes from the
/// ALT (`held`) or one telling the ALT (`alt`), and whether any read of it spells the allele in
/// full (`shown`).
#[derive(Clone, Default)]
struct Tally {
    told: Option<bool>,
    conflict: bool,
    held: bool,
    alt: bool,
    /// A read of the fragment spells the allele in full, placed here or not (`spells_in_full`).
    shown: bool,
}

impl Tally {
    /// The fragment's mates disagree: they read in full and differ (`conflict`), or one spelled
    /// the allele in full and a mate of any length confidently told the ALT over it. A mate that
    /// ends before it holds the allele in full still tells the ALT where it reads a base the two
    /// differ by; the fragment then spells neither, and counts nowhere
    /// (RJ-25).
    fn disagrees(&self) -> bool {
        self.conflict || (self.told == Some(false) && self.alt)
    }
}

/// [`present_hidden`]'s finding: the ALT's bases and the present other alleles' (from the
/// windows' first base), those within a read's reach that are not present (a read that spells
/// one in full is no ALT: partial evidence, RJ-6), the ALT with each listed row applied, REF's
/// bases (for [`spells_in_full`]) and the allele past the longest read's reach.
#[derive(Default)]
struct PresentHidden {
    alt: Vec<u8>,
    others: Vec<Hidden>,
    absent: Vec<Hidden>,
    neighbours: Vec<Vec<u8>>,
    reference: Vec<u8>,
    unspanned: Option<Unspanned>,
}

/// The allele a fragment's bases past the windows place it as (`read`: its base at an
/// edged offset and whether it is read confidently; `windows`: each allele's edged
/// bases; `n`: the windows' length), read on as RJ-27 reads a pure indel's: the
/// [`REJOIN`] bases past them fit that allele to within one base and the other by two
/// or more, and fit each allele shifted by up to [`SHIFT_REACH`] bases either way worse
/// by two or more than they fit it (in a repeat a read of another length of either
/// allele reads the same there); and the [`PAST`] bases it reads on along the allele
/// fit it to within one base all told (a read of another allele whose bases past the
/// windows happen to read the same goes on to differ). None where the fragment does not
/// read all [`PAST`] of them (one ending sooner may end before it differs), or they
/// place it nowhere unambiguously.
fn placed_past(read: &dyn Fn(usize) -> Option<(u8, bool)>, windows: &[Vec<Option<u8>>; 2], n: usize) -> Option<usize> {
    let past = n + 1..=n + REJOIN;
    if !(n + 1..=n + PAST).all(|e| read(e).is_some()) {
        return None;
    }
    let misses = |a: usize, s: i64, es: std::ops::RangeInclusive<usize>| {
        es.filter(|&e| {
            let at = usize::try_from(e as i64 + s).ok().and_then(|i| windows[a].get(i).copied().flatten());
            matches!((read(e), at), (Some((b, true)), Some(x)) if b != x)
        })
        .count()
    };
    let allele = match (misses(0, 0, past.clone()), misses(1, 0, past.clone())) {
        (0..=1, 2..) => 0,
        (2.., 0..=1) => 1,
        _ => return None,
    };
    let fit = misses(allele, 0, past.clone());
    let shifted = (1..=SHIFT_REACH).flat_map(|k| [k, -k]).any(|s| (0..2).any(|x| misses(x, s, past.clone()) < fit + 2));
    let read_on = misses(allele, 0, n + 1..=n + PAST) <= 1;
    (!shifted && read_on).then_some(allele)
}

/// Whether a run of an allele's own changes (edged offsets [e0, e1]) lies past the
/// junction and its bases rejoin the allele's continuation after it (#285): each of the
/// [`REJOIN`] bases after its last is read by [`GUARD_MIN_READS`] of the allele's
/// fragments or more (`cover`), and no change recurs on any (`recurs`); `m` is the
/// guard's edged reach.
fn rejoins(j: &Junction, (e0, e1): (usize, usize), m: usize, cover: &dyn Fn(usize) -> u32, recurs: &dyn Fn(usize) -> bool) -> bool {
    let after = e1 + 1..=e1 + REJOIN;
    e0 > K + j.mh && *after.end() < m && after.clone().all(|e| cover(e) >= GUARD_MIN_READS && !recurs(e))
}

/// Whether another breakpoint gives an allele's bases as its reads show them, exactly
/// (#285, #296): the allele from the windows' first base through [`REJOIN`] bases past
/// its last change, each change (`changed`, edged offsets, ascending) as most of its
/// reads show it (`most`), each base `compared` does not take (another recurring change,
/// or one too few fragments read) fitting any; the [`REJOIN`] bases after the last change
/// must all be compared, or no breakpoint is established. The other breakpoints are the reference
/// with one deletion at another start or of another length, and the allele with one more
/// insertion or deletion from the junction through its last change, within the reference
/// held (as many bases past the deletion as it deletes: [`reference_past`]). A run at the
/// junction whose bases are the deleted bases' own is the deletion shifted right; a
/// change or run before a repeat that a deletion of another length gives whole (another
/// unit count of a VNTR) is that deletion. The reads' bases are a consensus, so another
/// breakpoint must give them exactly: one that differs at a base is the given deletion's
/// with a neighbour as much as its own.
fn another_breakpoint(
    j: &Junction,
    base: &dyn Fn(i64) -> Option<u8>,
    alt: bool,
    changed: &[usize],
    compared: &dyn Fn(usize) -> bool,
    most: &dyn Fn(usize) -> u8,
) -> bool {
    let Some(&last) = changed.last() else { return false };
    let cont = K + j.mh;
    // The allele's base at window offset o, on along it past the windows.
    let own = |o: usize| base(j.genome(alt, o));
    let span = last + REJOIN;
    let shows: Vec<Option<u8>> = (0..span)
        .map(|o| match o + 1 {
            e if changed.contains(&e) => Some(most(e)),
            e if !compared(e) => None,
            _ => own(o),
        })
        .collect();
    if shows[last..].iter().any(Option::is_none) {
        // The bases after the last change not all read: no breakpoint is established.
        return false;
    }
    // Whether bases `at` gives (None: off the reference) are the reads' over window
    // offsets `os`, each base read (an offset not compared fits).
    let matches = |at: &dyn Fn(usize) -> Option<u8>, os: std::ops::Range<usize>| {
        os.into_iter().all(|o| match (at(o), shows[o]) {
            (None, _) => false,
            (Some(x), Some(s)) => x == s,
            (Some(_), None) => true,
        })
    };
    // Every alternative leaves the bases before its breakpoint as they stand and shifts
    // the rest, the bases after the last change among them: the shifts that give those
    // first, then each breakpoint up to where the unshifted bases stop being the reads'.
    let tail = last..span;
    let reach = |at: &dyn Fn(usize) -> Option<u8>| (0..).take_while(|&d| at(span - 1 + d).is_some()).count();
    let first = j.first();
    let genome = |o: usize| base(first + o as i64);
    let unshifted = (0..span).find(|&o| !matches(&genome, o..o + 1)).unwrap_or(span);
    let deletion = (0..reach(&genome))
        .filter(|&d| matches(&|o| genome(o + d), tail.clone()))
        .any(|d| (0..=unshifted.min(last)).any(|x| matches(&|o| genome(o + d), x..last)));
    // The allele with one more deletion, or an insertion among its changes, at an
    // offset from the junction through its first change.
    let kept = (0..span).find(|&o| !matches(&own, o..o + 1)).unwrap_or(span).min(last - 1);
    let more_deleted = (1..reach(&own))
        .filter(|&k| matches(&|o| own(o + k), tail.clone()))
        .any(|k| (cont..=kept).any(|p| matches(&|o| own(o + k), p..last)));
    let inserted = (cont..=kept).any(|p| (1..=last - p).any(|k| matches(&|o| own(o - k), p + k..span)));
    deletion || more_deleted || inserted
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

    /// REF = flank, u1, `m` copies of u, flank (#297): u is A then 59 C/G/T bases, u1 is u
    /// with an A at its third base, the flanks C/G/T only (so neither runs into u); the
    /// row deletes u1 after the anchor at 99.
    fn beside_run(m: usize) -> (String, Variant) {
        let mut left: Vec<u8> = cgt(100, 3).into_bytes();
        let u = format!("A{}", cgt(59, 5));
        let mut u1 = u.clone().into_bytes();
        u1[2] = b'A';
        if left[99] == u.as_bytes()[59] {
            left[99] = if u.as_bytes()[59] == b'C' { b'G' } else { b'C' };
        }
        let ctx = format!("{}{}{}{}", String::from_utf8(left).unwrap(), String::from_utf8(u1).unwrap(), u.repeat(m), cgt(400, 7));
        let v = del(&ctx, 99, 60);
        (ctx, v)
    }

    #[test]
    fn one_more_unit_deleted_reads_like_the_alt_until_the_run_s_last_copy() {
        // The ALT continues u u, one more unit deleted u then the flank: they first differ
        // 60 bases past the junction, K + 60 from the windows' first base.
        let (ctx, v) = beside_run(2);
        let j = junction(&v).unwrap();
        assert_eq!((j.lo, j.mh, j.ref_w.len()), (100, 2, 2 * K + 2));
        let (alt, hidden) = hidden_alleles(&v, &[], 151).unwrap();
        let two = hidden.iter().find(|h| h.change == -120).expect("one more unit deleted");
        assert_eq!((two.first_diff, two.need(), two.listed), (K + 60, K + 60 + PAST, false));
        assert_eq!(two.row, (99, ctx.as_bytes()[99..220].to_vec(), vec![ctx.as_bytes()[99]]));
        assert_eq!(&alt[..K], &ctx.as_bytes()[96..100]);
        assert_eq!(&alt[K..K + 120], &ctx.as_bytes()[160..280], "the ALT reads the run on");
        assert_eq!(&two.bases[K..K + 60], &ctx.as_bytes()[220..280]);
        assert_ne!(alt[two.first_diff], two.bases[two.first_diff]);
        assert_eq!(alt[..two.first_diff], two.bases[..two.first_diff]);
        // Every other allele found differs past the windows (none inside them is listed), and
        // each deletes more bases or fewer: no insertion at the breakpoint is looked for.
        assert!(hidden.iter().all(|h| h.first_diff >= j.ref_w.len() && h.change < 0));
        // A reference N inside the run tells nothing: the first difference stays put.
        let mut with_n = ctx.clone().into_bytes();
        with_n[190] = b'N';
        let with_n = del(&String::from_utf8(with_n).unwrap(), 99, 60);
        let (_, hidden) = hidden_alleles(&with_n, &[], 151).unwrap();
        assert_eq!(hidden.iter().find(|h| h.change == -120).map(|h| h.first_diff), Some(K + 60));
        // With four copies, one more unit deleted reads alike for three of them; a deletion
        // of more than twice the given bases is not looked for (as RJ-24's other
        // breakpoints).
        let (_, hidden) = hidden_alleles(&beside_run(4).1, &[], 151).unwrap();
        assert_eq!(hidden.iter().find(|h| h.change == -120).map(|h| h.first_diff), Some(K + 180));
        assert!(hidden.iter().all(|h| h.change >= -120));
    }

    #[test]
    fn a_listed_unit_count_is_the_allele_found_at_the_breakpoint() {
        // The two-unit deletion listed (left-aligned as prep writes it) is the allele the
        // breakpoint gives, marked listed; a listed SNV, or a deletion before the windows,
        // is not one.
        let (ctx, v) = beside_run(2);
        let two = del(&ctx, 99, 120);
        let snv = Variant::new("1".into(), 230, ctx[230..231].into(), "A".into(), "SNP".into(), None, 0, 0, None, None, None, None);
        let before = del(&ctx, 40, 50);
        let (_, hidden) = hidden_alleles(&v, &[two, snv, before], 151).unwrap();
        let listed: Vec<i64> = hidden.iter().filter(|h| h.listed).map(|h| h.change).collect();
        assert_eq!(listed, vec![-120]);
    }

    #[test]
    fn a_listed_count_whose_ref_runs_past_the_reference_held_is_an_unresolved_candidate() {
        // The reference held ends 300 bases past the junction; the listed 600-bp deletion at the
        // same anchor ends at 700. The bases after its last deleted base are not held: it cannot
        // be built, so it stays a candidate past any read's reach, named where it is the nearest.
        let (ctx, v) = beside_run(12);
        let held = Variant { event_ref: Some((0, ctx[..400].into())), ..v.clone() };
        let (_, hidden) = hidden_alleles(&held, &[del(&ctx, 99, 600)], 151).unwrap();
        let unresolved: Vec<&Hidden> = hidden.iter().filter(|h| h.unresolved()).collect();
        assert_eq!(unresolved.len(), 1);
        let h = unresolved[0];
        assert_eq!((h.change, h.listed, h.bases.len()), (-600, true, 0));
        assert_eq!(h.row, (99, ctx.as_bytes()[99..700].to_vec(), vec![ctx.as_bytes()[99]]));
        assert!(h.need() > 151, "past any read's reach");
        assert_eq!(nearest_counts(&hidden, -60).iter().filter(|h| h.unresolved()).count(), 0, "the generated one deleting a unit more is nearer");
        let alone: Vec<Hidden> = hidden.iter().filter(|h| h.listed).cloned().collect();
        assert_eq!(nearest_counts(&alone, -60).iter().map(|h| h.change).collect::<Vec<_>>(), vec![-600]);
        // a listed deletion from inside the run, also running past the reference held, differs from
        // the ALT inside the windows by the bases held: told there, nothing to hold
        let (_, hidden) = hidden_alleles(&held, &[del(&ctx, 219, 600)], 151).unwrap();
        assert!(hidden.iter().all(|h| !h.listed), "told inside the windows");
        // a listed row whose held part does not match the reference is no row of this site
        let mut wrong = del(&ctx, 99, 600);
        wrong.ref_allele = format!("{}{}", &ctx[99..300], "A".repeat(400));
        let (_, hidden) = hidden_alleles(&held, &[wrong], 151).unwrap();
        assert!(hidden.iter().all(|h| !h.listed));
    }

    #[test]
    fn an_allele_that_is_the_alt_with_a_listed_neighbour_is_not_another_allele() {
        // The ALT continues T: 12-base units for 40 + 12 bases, then an A, then 13-base units.
        // One more unit deleted reads T from its 12th base: it first differs at T's 40th (the
        // A) and on is the ALT with that A inserted after its 39th base. Listed as an input row,
        // that insertion is the carriers' neighbour (RJ-24): no other allele.
        let mut left: Vec<u8> = cgt(100, 23).into_bytes();
        let deleted = format!("A{}", cgt(59, 29));
        let u = cgt(12, 31).into_bytes();
        let (p, at) = (12usize, 40usize);
        let mut t: Vec<u8> = u.clone();
        while t.len() < 200 {
            let i = t.len();
            t.push(if i < at + p { t[i - p] } else if i == at + p { b'A' } else { t[i - p - 1] });
        }
        if left[99] == deleted.as_bytes()[59] {
            left[99] = if deleted.as_bytes()[59] == b'C' { b'G' } else { b'C' };
        }
        let ctx = format!("{}{}{}{}", String::from_utf8(left).unwrap(), deleted, String::from_utf8(t.clone()).unwrap(), cgt(300, 37));
        let v = del(&ctx, 99, 60);
        let (_, hidden) = hidden_alleles(&v, &[], 151).unwrap();
        let more = hidden.iter().find(|h| h.change == -72).expect("one more unit deleted");
        assert_eq!(more.first_diff, K + at);
        let g = 160 + at as i64 - 1; // the anchor of the insertion, on the ALT's continuation
        let row = (g, (t[at - 1] as char).to_string(), format!("{}A", t[at - 1] as char));
        let listed = Variant { listed_indels: vec![row], ..v.clone() };
        let (_, hidden) = hidden_alleles(&listed, &[], 151).unwrap();
        assert!(hidden.iter().all(|h| h.change != -72), "the ALT with its listed neighbour");
    }

    #[test]
    fn the_diagnostic_weighs_the_nearest_unit_count_on_each_side() {
        // Given a 60-bp deletion: a 2-bp shift reading alike for 6 bases by chance, one unit
        // (-120 and 0) each side, and a far count reading alike for 160 bases.
        let h = |change: i64, agree: usize| Hidden { change, bases: Vec::new(), first_diff: K + agree, listed: false, row: (0, Vec::new(), Vec::new()), spelled: 0, told: 0, held: 0 };
        let all = [h(-62, 6), h(-120, 100), h(-58, 5), h(-30, 180), h(-240, 160)];
        let near: Vec<i64> = nearest_counts(&all, -60).iter().map(|h| h.change).collect();
        assert_eq!(near, vec![-30, -120], "fewer by 30 (180 >= 30 + 12) and more by 60, not the chance shifts or the far count");
        let short = [h(-120, 60)];
        assert!(nearest_counts(&short, -60).is_empty(), "60 bases alike for a 60-base shift: no copy twice over");
    }

    #[test]
    fn mates_disagree_when_one_spells_the_allele_and_a_mate_of_any_length_tells_the_alt() {
        let t = |told: Option<bool>, conflict: bool, alt: bool| Tally { told, conflict, alt, ..Tally::default() };
        assert!(t(Some(false), false, true).disagrees(), "a spelling mate and a mate telling the ALT, however short");
        assert!(t(Some(true), true, true).disagrees(), "two full mates that differ");
        assert!(!t(Some(false), false, false).disagrees(), "the allele spelled, no mate tells the ALT");
        assert!(!t(Some(true), false, true).disagrees(), "the ALT told, no mate spells the allele");
        assert!(!t(None, false, true).disagrees(), "a mate telling the ALT alone");
    }

    /// A deletion of `len` bases after the anchor at 99, then a simple repeat (AC)n of `len`
    /// bases and a random flank: every shift by an even number of bases into the repeat reads
    /// like the ALT up to the repeat's end.
    fn before_a_simple_repeat(len: usize) -> (String, Variant) {
        let (mut left, deleted) = (cgt(100, 3).into_bytes(), cgt(len, 5));
        if left[99] == deleted.as_bytes()[len - 1] {
            left[99] = if deleted.as_bytes()[len - 1] == b'C' { b'G' } else { b'C' };
        }
        let ctx = format!("{}{}{}{}", String::from_utf8(left).unwrap(), deleted, "AC".repeat(len / 2), cgt(400, 7));
        let v = del(&ctx, 99, len);
        (ctx, v)
    }

    #[test]
    fn an_allele_past_a_reads_reach_keeps_a_small_record() {
        // Every even shift into a long repeat is another allele reading like the ALT to the
        // repeat's end: about len/2 of them, each first differing thousands of bases on, far past
        // any read. Building each one's bases and row held memory growing with the square of the
        // deletion (a 24-kb deletion held some 580 MB); only the alleles a read could spell, and
        // the few whose bases the listed-neighbour test reads up to the read's end, keep any.
        let reach = 151;
        for (len, check) in [(2_000usize, true), (24_000, false)] {
            let (ctx, v) = before_a_simple_repeat(len);
            let j = junction(&v).unwrap();
            let started = std::time::Instant::now();
            let (alt, hidden) = hidden_alleles(&v, &[], reach).unwrap();
            let kept: usize = hidden.iter().map(|h| h.bases.len() + h.row.1.len() + h.row.2.len()).sum();
            // the alleles a read could spell keep a row of about the deletion's length (at most
            // one for each base the read holds); none of the others keeps anything
            let bound = reach * 2 * len + 1_000_000;
            assert!(kept < bound, "{len}: {kept} bytes kept for {} alleles (bound {bound})", hidden.len());
            assert!(started.elapsed().as_secs() < 30, "{len}: {:?}", started.elapsed());
            assert!(hidden.len() >= (len - 4) / 2, "{len}: {} alleles", hidden.len());
            if !check {
                continue;
            }
            // each against a plain reading of the reference: where it first differs, its length
            // change, and, where a read could spell it, its bases and row
            let (b, lo) = (ctx.as_bytes(), j.lo as usize);
            for h in &hidden {
                let from = (lo as i64 - h.change) as usize;
                let agree = (0..).take_while(|&i| lo + len + i < b.len() && from + i < b.len() && b[lo + len + i] == b[from + i]).count();
                assert_eq!(h.first_diff, K + agree, "{:+}", h.change);
                if h.need() <= reach {
                    let want: Vec<u8> = b[lo - K..lo].iter().chain(&b[from..]).take(h.first_diff + HIDDEN_REACH).copied().collect();
                    assert_eq!(h.bases, want, "{:+}", h.change);
                    assert_eq!(h.row, (lo as i64 - 1, b[lo - 1..from].to_vec(), vec![b[lo - 1]]), "{:+}", h.change);
                }
            }
            assert!(alt.len() > len, "the ALT is held as far as the reference is");
            // a listed row that is one of them (the same haplotype) brings its bases and row back
            let (_, hidden) = hidden_alleles(&v, &[del(&ctx, 99, len + 2)], reach).unwrap();
            let h = hidden.iter().find(|h| h.change == -(len as i64 + 2)).expect("the listed count");
            let from = lo + len + 2;
            let want: Vec<u8> = b[lo - K..lo].iter().chain(&b[from..]).take(h.first_diff + HIDDEN_REACH).copied().collect();
            assert!(h.listed && !h.unresolved());
            assert_eq!((&h.bases, &h.row), (&want, &(lo as i64 - 1, b[lo - 1..from].to_vec(), vec![b[lo - 1]])));
            assert_eq!(hidden.len(), hidden.iter().filter(|x| !x.listed).count() + 1, "merged, not added");
        }
    }

    #[test]
    fn a_read_tells_two_alleles_only_where_they_differ() {
        // The alleles agree on their first 8 bases, then differ at 8 and 9 (and agree at 10).
        let (a, b) = (b"CCCCGGGGAAT".as_slice(), b"CCCCGGGGTTT".as_slice());
        let tell = |read: &[u8], quals: &[u8], end: usize| tells(read, quals, 20, 0, (a, b), 8, end);
        let q = [30u8; 11];
        assert_eq!(tell(b"CCCCGGGGAAT", &q, 11), Some(true));
        assert_eq!(tell(b"CCCCGGGGTTT", &q, 11), Some(false));
        // an error where they agree costs both alike; a third base where they differ too
        assert_eq!(tell(b"CACCGGGGAAT", &q, 11), Some(true));
        assert_eq!(tell(b"CCCCGGGGGAT", &q, 11), Some(true));
        // ending before they differ, or reading their differing bases as N or below
        // --min-baseq, tells nothing
        assert_eq!(tell(b"CCCCGGGGAAT", &q, 8), None);
        assert_eq!(tell(b"CCCCGGGGNNT", &q, 11), None);
        let mut low = q;
        low[8] = 10;
        low[9] = 10;
        assert_eq!(tell(b"CCCCGGGGAAT", &low, 11), None);
        // one base of each: neither is told
        assert_eq!(tell(b"CCCCGGGGATT", &q, 11), None);
    }

    #[test]
    fn a_read_shows_a_listed_rows_change_only_where_it_reads_the_rows_base() {
        // The ALT's bases, the ALT with a listed SNV at offset 8 (C for G), and with a listed
        // 1-bp insertion after offset 10 (an A): a read is the deletion with that neighbour
        // only where its confident bases, past `from`, read the row's change.
        let alt = b"CCCCGGGGCCTTTGGAAC".as_slice();
        let mut snv = alt.to_vec();
        snv[8] = b'G';
        let ins = [&alt[..10], b"TA".as_slice(), &alt[11..]].concat();
        let q = [30u8; 19];
        let shows = |read: &[u8], quals: &[u8], x: &[u8], from: usize, end: usize| shows_change(read, quals, 20, 0, alt, x, from, end);
        assert!(shows(&snv, &q, &snv, 8, 18), "a read with the SNV's base shows it");
        assert!(!shows(alt, &q, &snv, 8, 18), "the ALT's base there: no");
        assert!(!shows(&snv, &q, &snv, 9, 18), "the change lies before `from`: the read cannot show it there");
        assert!(!shows(&snv, &q, &snv, 8, 8), "the read ends before it");
        let mut low = q;
        low[8] = 10;
        assert!(!shows(&snv, &low, &snv, 8, 18), "below --min-baseq: not confident");
        assert!(shows(&ins, &q, &ins, 8, 19), "a read with the insertion shifts onto the row's bases");
        assert!(!shows(alt, &q, &ins, 8, 18), "the ALT's bases over the same offsets: no");
        // a deletion row 270 bases past the read changes nothing the read holds: x is the ALT
        assert!(!shows(alt, &q, alt, 8, 18), "no difference to show");
    }

    #[test]
    fn a_read_spells_an_allele_only_where_it_fits_it_and_not_ref() {
        // From offset 4: the ALT a, the allele b (differing at 4, 5 and 9), REF r.
        let (a, b, r) = (b"CCCCAAGGGAGGG".as_slice(), b"CCCCTTGGGTGGG".as_slice(), b"CCCCTTCCCTCCC".as_slice());
        let q = [30u8; 13];
        let spell = |read: &[u8]| spells(read, &q, 20, 0, (a, b, r), 4, 13);
        assert!(spell(b"CCCCTTGGGTGGG"), "b's bases where it and a differ, and r worse by six");
        assert!(spell(b"CCCCTTGGCTGGG"), "one other change aside");
        assert!(!spell(b"CCCCTTGCCTGGG"), "two other changes");
        assert!(!spell(b"CCCCTAGGGTGGG"), "a's base where they differ");
        assert!(!spell(b"CCCCTCGGGNGGG"), "b's base at one telling base only");
        // a read REF fits as well (the sample's own copy, read past windows it fits by chance)
        let r2 = b"CCCCTTGGGTGGC".as_slice();
        assert!(!spells(b"CCCCTTGGGTGGG", &q, 20, 0, (a, b, r2), 4, 13), "REF misses one only");
    }

    #[test]
    fn one_error_against_three_telling_bases_is_forgiven() {
        // Differing at 4..8: a read showing b's base at one and a's at three is a, as RJ-27
        // reads on past the windows; at one against two it tells neither.
        let (a, b) = (b"CCCCAAAAG".as_slice(), b"CCCCTTTTG".as_slice());
        let q = [30u8; 9];
        assert_eq!(tells(b"CCCCTAAAG", &q, 20, 0, (a, b), 4, 9), Some(true));
        assert_eq!(tells(b"CCCCTAANG", &q, 20, 0, (a, b), 4, 9), None);
    }

    #[test]
    fn a_read_spells_another_count_in_full_only_where_it_holds_the_bases_and_tells_it() {
        // The ALT and another count agree for 6 bases, then differ at 6, 9 and 12; REF is
        // other bases throughout. The allele is held [`PAST`] bases past its first difference.
        let common = cgt(60, 3).into_bytes();
        let (mut alt, mut other) = (common.clone(), common.clone());
        for k in [6, 9, 12] {
            alt[k] = b'A';
            other[k] = b'T';
        }
        let reference = cgt(60, 9).into_bytes();
        let h = Hidden { change: -60, bases: other.clone(), first_diff: 6, listed: false, row: (0, Vec::new(), Vec::new()), spelled: 0, told: 0, held: 0 };
        let q = [30u8; 60];
        let told = |read: &[u8], end: usize| tells(read, &q, 20, 0, (&alt, &other), 6, end);
        let spelt = |read: &[u8], end: usize| spells_in_full(read, &q, 20, (0, end), told(read, end), (&alt, &reference, &[]), 6, &h);
        assert!(spelt(&other, 60), "its bases are the allele's, held to 24 past the first difference");
        assert!(!spelt(&alt, 60), "the ALT's bases tell the ALT");
        assert!(!spelt(&other, 29), "it ends before 24 bases past the first difference");
        assert!(spelt(&other, 30));
        // the deletion with a listed neighbour (here the ALT with a row giving the allele's bases)
        // is never another allele
        let neighbours = vec![other.clone()];
        assert!(!spells_in_full(&other, &q, 20, (0, 60), told(&other, 60), (&alt, &reference, &neighbours), 6, &h));
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
            guard_for(&all, &v, &|_| true, 20, &[])
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
        let g = guard_for(&refs(&[last, last + 1]), &v, &|_| true, 20, &[]).expect("a run that rejoins REF");
        assert_eq!((g.masked, g.larger), (at_last.clone(), None));
        assert!(guard_for(&refs(&keeps(last)), &v, &|_| true, 20, &[]).is_none());
        let g = guard_for(&refs(&[last]), &v, &|_| true, 20, &[]).expect("an isolated REF change");
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
        let g = guard_for(&reads(&changed(&hap, &[61, 62]), d, 40, 4, "c"), &v, &|_| true, 20, &[]).expect("carriers");
        assert_eq!(g.masked, inside);
        assert!(g.larger.is_some());
        let refs = reads(&changed(&ctx.as_bytes()[40..140], &[61, 62]), vec![Cigar::Match(100)], 40, 4, "r");
        assert_eq!(guard_for(&refs, &v, &|_| true, 20, &[]).expect("REF reads").masked, inside);
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
        let at = |cigar: Vec<Cigar>, seq: Vec<u8>| guard_for(&reads(&seq, cigar, 40, 4, "c"), &v, &|_| true, 20, &[]);
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
        let g = guard_for(&[refs, alts].concat(), &v, &|_| true, 20, &[]).expect("both changes recur");
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
        let g = guard_for(&[clipped, refs].concat(), &v, &|_| true, 20, &[]).expect("a recurring change");
        let l = g.larger.expect("a larger allele");
        assert_eq!(l.ref_allele.len() - l.alt_allele.len(), 60, "{l:?}");
    }

    #[test]
    fn the_bases_past_the_windows_place_a_read_only_unambiguously() {
        // Edged windows of 8 bases with 24 past them; a read reading ALT's bases.
        let seq = |seed: u32| -> Vec<Option<u8>> { cgt(33, seed).bytes().map(Some).collect() };
        let (r, alt) = (seq(3), seq(5));
        let reads = |w: &[Option<u8>]| {
            let bases: Vec<u8> = w.iter().map(|b| b.unwrap()).collect();
            move |e: usize| Some((bases[e], true))
        };
        let windows = [r.clone(), alt.clone()];
        assert_eq!(placed_past(&reads(&alt), &windows, 8), Some(1));
        assert_eq!(placed_past(&reads(&r), &windows, 8), Some(0));
        // ALT's bases past the windows are REF's one base on (a repeat): a read of REF
        // with one base deleted reads them too, so they place no read
        let mut shifted = alt.clone();
        shifted[9..32].copy_from_slice(&r[10..33]);
        assert_eq!(placed_past(&reads(&shifted), &[r.clone(), shifted.clone()], 8), None);
        // a read reading ALT's twelve bases past the windows, then two others of its own
        let mut other = alt.clone();
        for e in [23, 26] {
            other[e] = Some(if alt[e] == Some(b'A') { b'C' } else { b'A' });
        }
        assert_eq!(placed_past(&reads(&other), &[r, alt], 8), None);
    }

    #[test]
    fn a_reference_too_short_for_the_windows_gives_none() {
        let ctx = cgt(162, 17);
        assert!(junction(&del(&ctx, 99, 60)).is_none()); // 2 right-flank bases only
    }
}
