//! Exact-carrier classification for complex variants.
//!
//! A complex variant (a delins, a deletion whose anchor also changes, an
//! insertion whose ALT changes the anchor such as C>TA, or an MNP read with an
//! indel or clip at its block) is REF or ALT for a read only when the read's own
//! bases carry that allele across the whole event, with [`FLANK`] reference bases
//! on each side. The read's bases include soft clips, except an RNA read's clip
//! that reaches an exon edge or a junction end (STAR soft-clips a junction
//! overhang it cannot splice, so such a clip holds the next exon's bases). Bases
//! below min BQ (and N) match anything, the one quality rule every backend shares. Nothing else is
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
//! - **One pair of equal-length windows at one flank.** Both alleles are read
//!   inward from the same flank: the shorter allele's whole window (its event with
//!   [`FLANK`] bases each side), and as many of the longer allele's bases. A read
//!   holding them is also read on inward, as far as it reaches, against the rest of
//!   the longer allele's window: a read whose later bases contradict an allele (one A
//!   more than the ALT's run) is not that allele. Windows of equal length at one
//!   anchor give a REF and an ALT molecule the same number of read starts that can
//!   judge them, so neither allele is favoured by where reads start, with no padding
//!   (which cost reads) and no second junction (which lets the longer allele count at
//!   both ends). The flank is the left (POS) one; an RNA variant whose left flank an
//!   exon edge cuts reads from the right. When both alleles are longer than
//!   [`LONG_EVENT`] bases the windows end one base past the first difference
//!   instead, read from the flank nearer it (an event ending a long run is out of a
//!   read's reach from the other). A read judged from the right counts though it
//!   starts past the variant position, as the other allele's molecule starting at
//!   the same place would. Past the shorter allele's event its flank may repeat the
//!   longer allele's next bases: both windows run on through them to [`FLANK`]
//!   bases from where the alleles differ, so a read carrying only a substitution of
//!   the event's bases is neither allele.
//! - **Placed from the other flank.** A read whose reading flank lies in its soft
//!   clip, under its own deletion, or aligned at the event's other copy is placed
//!   from the other flank as the allele its bases favour would sit there: each
//!   allele's window placed from that flank at its own length, its confident bases
//!   (at or above min BQ, not N) must favour one by at least one base read at min BQ.
//!   It counts where that allele fits, else carries another allele; one placed
//!   through aligned bases counts only where the allele fits. A read aligned straight
//!   across the windows that fits the shorter allele from its reading flank, while
//!   read whole from the other flank its bases are the longer allele's, had that
//!   allele's other copy written over the reading flank: the longer allele where it
//!   holds that allele's window from there and reaches past it, else depth only.
//! - **A change beside the event.** A change that recurs at a fixed distance from
//!   the event is a separate event, as another indel outside the window is, masked
//!   in both windows ([`guard_for`]) so a read carrying the given allele counts
//!   whatever it shows there: one the REF fragments carry too (a germline SNP or
//!   MNP), and an isolated one on the ALT fragments. A run of changes on the ALT
//!   fragments alone is an indel or MNP beside the event, another allele. The ALT
//!   fragments' own isolated changes make their haplotype larger than the given
//!   allele; the diagnostic names that allele as any observed allele is (beside a
//!   run, the observed allele names what the reads' alignments carry).
//! - **Outcomes.** A read decides only where it holds both the REF and the ALT
//!   window. A read matching both (through masked bases) is neither. A read
//!   holding no pair is depth only, like a pure-indel read ending inside its
//!   tract: no allele, no partial evidence, no mFSD class. A read holding the
//!   windows and matching neither carries another allele: partial evidence when
//!   closer to ALT.
//! - **Spliced reads.** A read spliced in a window's flank is judged on the
//!   windows built (event, growth and flank) over the reference spliced at its
//!   own junctions: its bases past a splice are the next exon's, so the alleles'
//!   haplotypes continue there, and a length change that pushes the exon's bases
//!   past the junction shows against them. A junction that starts inside those
//!   bases and runs past them (or ends inside them) splices the haplotypes at the
//!   event's edge, so the bases decide however the aligner wrote the gap; a splice
//!   through the bases where the alleles differ otherwise leaves the read depth
//!   only.

use std::collections::HashMap;
use std::rc::Rc;

use log::trace;
use rust_htslib::bam::record::{Cigar, Record};

use super::observed::{canonical, ObservedAllele};
use super::rna;
use super::variant_checks::deleted_anchor_query;
use super::window::AlleleKind;
use super::utils::{find_read_pos, median_qual, soft_clips, ClassifyPhase, ClassifyResult};
use crate::normalize::fasta::CachedFasta;
use crate::shared::bam_utils::fragment_query_span;
use crate::types::Variant;

/// How the rule reads a read past its aligned blocks, by library.
#[derive(Clone, Copy)]
pub(crate) struct ReadRules<'a> {
    /// RNA: the exon edges and junction ends near the variant (annotated intron
    /// boundaries and the junctions its reads splice at). A soft clip whose
    /// reference reach crosses one, or whose aligned bases end within
    /// [`EDGE_SLACK`] of one, may hold the next exon's bases and is not read; any
    /// other clip is the read's own bases. None (DNA): every clip is read.
    pub clip_edges: Option<&'a ClipEdges<'a>>,
    /// The reference, for the far side of a spliced read's junctions (RNA).
    pub reference: Option<&'a CachedFasta>,
    /// Spliced windows already built for this variant's reads (RNA).
    pub spliced: Option<&'a SplicedCache>,
    /// The variant's neighbouring-change guard, found from its reads on first use.
    pub guard: Option<&'a CarrierGuard<'a>>,
    /// The substitutions the variant's reads share, found from its reads on first
    /// use: no error, so never tolerated. None tolerates every one.
    pub tolerance: Option<&'a Tolerance<'a>>,
}

impl ReadRules<'static> {
    /// DNA, no guard: clipped bases are read; no read is spliced.
    pub(crate) const DNA: ReadRules<'static> = ReadRules { clip_edges: None, reference: None, spliced: None, guard: None, tolerance: None };
}

/// The exon edges and junction ends near a variant, found on first use: only the
/// exact-carrier rule reads them, so an SNV row never pays for them. `find` gives the
/// exon starts (the ends of introns, where the bases before are intronic) and the
/// exon ends (the first bases of introns, where the bases from on are intronic).
pub(crate) struct ClipEdges<'a> {
    edges: std::cell::OnceCell<Edges>,
    find: Box<dyn Fn() -> (Vec<i64>, Vec<i64>) + 'a>,
}

/// Every edge, and the exon starts and ends among them, each ascending.
struct Edges {
    all: Vec<i64>,
    starts: Vec<i64>,
    ends: Vec<i64>,
}

impl<'a> ClipEdges<'a> {
    pub(crate) fn new(find: impl Fn() -> (Vec<i64>, Vec<i64>) + 'a) -> Self {
        Self { edges: std::cell::OnceCell::new(), find: Box::new(find) }
    }

    fn edges(&self) -> &Edges {
        self.edges.get_or_init(|| {
            let (mut starts, mut ends) = (self.find)();
            let mut all = [starts.as_slice(), ends.as_slice()].concat();
            for v in [&mut all, &mut starts, &mut ends] {
                v.sort_unstable();
                v.dedup();
            }
            Edges { all, starts, ends }
        })
    }

    fn get(&self) -> &[i64] {
        &self.edges().all
    }

    /// Whether an exon edge cuts the flank of `flank` bases before `lo`: an exon start
    /// in it or at `lo` (the bases before it intronic), or an exon end inside it.
    fn cut_left(&self, lo: i64, flank: i64) -> bool {
        let e = self.edges();
        starts_cut_left(e.starts.iter().copied(), lo, flank) || e.ends.iter().any(|&n| n > lo - flank && n < lo)
    }

    /// Whether an exon edge lies past `hi` within `flank` bases (the right flank's test
    /// as it stood before the edges were typed: an event filling a short exon, both
    /// flanks cut, keeps reading from the right).
    fn cut_right(&self, hi: i64, flank: i64) -> bool {
        self.edges().all.iter().any(|&x| x > hi && x <= hi + flank)
    }
}

/// How far an RNA read's aligned bases may run past an exon edge before a clip:
/// STAR extends an overhang shorter than its novel-junction minimum
/// (`alignSJoverhangMin`, 5) into the intron with mismatches rather than splice it.
const EDGE_SLACK: i64 = 5;

/// Spliced windows built for a variant's reads, keyed by the variant and the
/// splice segments: reads sharing a junction share their windows.
#[derive(Default)]
pub(crate) struct SplicedCache(std::cell::RefCell<HashMap<SplicedKey, Option<Rc<SplicedWindows>>>>);

/// The variant (position, REF, ALT) and the read's splice segments.
type SplicedKey = (i64, String, String, Vec<(i64, i64, i64)>);

/// Windows built over a spliced reference, with the map back to the genome.
type SplicedWindows = (Windows, SpliceMap);

/// Query positions `[first, after)` the rule may read: the whole read in DNA; in
/// RNA, a soft clip is left out when its reach crosses an exon edge or junction
/// end in `edges` (or the aligned bases end within [`EDGE_SLACK`] of one).
pub(crate) fn readable_range(record: &Record, rules: &ReadRules) -> Option<(usize, usize)> {
    let len = record.seq_len();
    let Some(edges) = rules.clip_edges.map(ClipEdges::get) else {
        return Some((0, len));
    };
    let (first, after) = aligned_query_range(record)?;
    let (lead, tail) = (first as i64, (len - after) as i64);
    let near = |lo: i64, hi: i64| {
        let i = edges.partition_point(|&b| b < lo);
        edges.get(i).is_some_and(|&b| b <= hi)
    };
    let (start, end) = (record.pos(), crate::shared::bam_utils::ref_end(record));
    let first = if lead > 0 && !near(start - lead, start + EDGE_SLACK) { 0 } else { first };
    let after = if tail > 0 && !near(end - EDGE_SLACK, end + tail) { len } else { after };
    Some((first, after))
}

/// Reference bases required on each side of the event.
const FLANK: usize = 2;
/// Reference bases prep holds each side of an event this rule judges, for reading a
/// read whole (`Variant::whole_ref`): a whole read's reach (2x250 and 2x300 runs
/// included), so a read is read whole over all its bases, never the part held.
pub(crate) const WHOLE_READ_REACH: i64 = 300;

/// Reference bases a pure insertion's or deletion's windows hold each side of its
/// tract: one base past it decides, as the base past a repeat does for the field's
/// read-level judges (RJ-3; operator, 2026-10-15). A read holds the shorter allele's
/// window (RJ-23), so a read ending inside the tract is depth only (RJ-1), and the
/// longer allele's read decides from its tract's last base as the shorter allele's
/// does from the base past its own.
const PURE_FLANK: usize = 1;
/// Window length beyond which the whole allele cannot be expected in one read.
const LONG_EVENT: usize = 50;
/// How far from the event a neighbouring change on the ALT reads is looked for.
const GUARD_REACH: usize = 25;
/// Reference bases prep fetches on each side of the event: the flank, and the
/// reach the neighbouring-change guard may grow it by.
const ROOM: usize = FLANK + GUARD_REACH;
/// Fragments, and share of the fragments reading that base, a neighbouring change
/// must recur on before it is taken as part of the reads' haplotype.
pub(crate) const GUARD_MIN_READS: u32 = 3;
pub(crate) const GUARD_MIN_SHARE: f64 = 0.2;
/// Fragments placed apart from the event that must read a change's base before its
/// absence there shows it is not background: at this many, a haplotype carrying it on
/// a fifth of them shows on three on average.
const GUARD_APART_READS: u32 = 15;
/// Share of the REF fragments a change at the ALT's shifted base may show on as
/// errors: on more, the fragments placed apart from the event are read for it.
const GUARD_APART_ERRORS: f64 = 0.01;
/// Share of the fragments placed apart from the event a change must show on, with
/// [`GUARD_MIN_READS`] of them, to be a haplotype's of their own: above the errors at
/// one base (random ones give a particular base on a third of 1% at Q20; a systematic
/// one, oxoG or a sequence context, may give a few percent at depth).
const GUARD_APART_SHARE: f64 = 0.03;
/// Chance below which the carriers' pair orientations, against the reads placed apart's
/// own balance, are no haplotype's (read in both as the library is) but an artifact's.
const GUARD_ORIENTATION_ALPHA: f64 = 1e-3;
/// Share of the REF reads on which a change is no longer the ALT's alone.
pub(crate) const GUARD_REF_SHARE: f64 = 0.05;
/// Longest tandem-repeat unit the event is grown through (besides the length change).
const MAX_UNIT: usize = 6;
/// Masked base: matches any haplotype base.
pub(crate) const WILD: u8 = b'N';
/// The bases a masked window position takes, one bit per base (A 1, C 2, G 4, T 8).
pub(crate) type Bases = u8;
/// Every base: a change beside the event is a separate event, so a read carrying the
/// allele counts whatever it shows there.
const ANY: Bases = 0b1111;
/// IUPAC codes by the bases they take: a window base masked to a set is written as its
/// code (N, [`WILD`], for any base), so a read base outside the set still mismatches.
const IUPAC: &[u8; 16] = b"-ACMGRSVTWYHKDBN";

/// A base's bit (none for N or anything else).
fn bit(b: u8) -> Bases {
    match b {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        _ => 0,
    }
}

/// The bases a window base takes: its own, or every base of its IUPAC code.
fn takes(h: u8) -> Bases {
    IUPAC.iter().position(|&c| c == h).map_or(0, |i| i as Bases)
}

/// Whether read base `b` fits window base `h`.
fn fits(b: u8, h: u8) -> bool {
    bit(b) & takes(h) != 0
}

/// Base `x`'s class at a position masked to `set` (0: unmasked): the whole set where
/// it holds `x`, else `x` alone. Two bases tell their alleles apart only where their
/// classes differ: a background change is no difference, a third base still is.
fn class(set: Bases, x: u8) -> Bases {
    if set & bit(x) != 0 {
        set
    } else {
        bit(x)
    }
}
/// Bases past the other flank a read's own bases must continue the reference over (one
/// error aside) for its alignment there to anchor a placement from that flank.
const OTHER_FLANK_CHECK: usize = 12;
/// Bases past the event a read placed from the left flank must fit as its allele's
/// (and, for a REF molecule, read): a read of the other length, shifted by the length
/// change, does not fit that many by chance.
const PLACED_REACH: usize = 8;
/// Changes beside the event a placed molecule may carry and still be its allele's:
/// two phased SNPs, or a germline DNP.
const PLACED_CHANGES: usize = 2;

/// One allele's window, its event (grown through repeats) with the flank on each
/// side. A read must hold `hold` of its bases from the reading flank, the same
/// number for both alleles, and is read on, as far as it reaches, through the rest.
#[derive(Clone)]
struct Window {
    /// The allele's whole window.
    allele: Vec<u8>,
    /// Bases a read must hold, counted from the window's left end (its right end
    /// when `tail`).
    hold: usize,
    /// Read from the right flank: the bases held are the window's last `hold`.
    tail: bool,
    /// Reference position of the window's first base (the left flank's), when a
    /// read can anchor on it.
    left: Option<i64>,
    /// Reference position one past the window's last base (the right flank's).
    right: Option<i64>,
    /// Bases of the window before the event (growth and flank) and after it,
    /// counted from where the alleles differ.
    before: usize,
    after: usize,
    /// The reference's bases past the window's right end, which both alleles continue
    /// with; none for windows read from the right or cut at an exon edge there.
    past: Vec<u8>,
    /// The reference's bases beyond the window's reading-side end (before its left end,
    /// or past its right end when `tail`), as many as the alleles' lengths differ, which
    /// both alleles share: the shorter allele extended by them is read at the longer
    /// one's length from the other flank. Cut at an exon edge there.
    outer: Vec<u8>,
    /// The event grown through repeats, from the first base the input names (a pure
    /// indel's anchor), as window positions [start, end): the rest is flank, beside
    /// the event.
    event_at: (usize, usize),
}

impl Window {
    /// Whether a read reads the window from its right flank: windows laid out from
    /// the right, or windows cut at an exon edge on the left (no flank base left
    /// there; the event's own first base is no anchor, as an aligner may write it as
    /// an insertion), while a right flank remains.
    fn reads_from_right(&self) -> bool {
        if self.tail { self.after > 0 } else { self.before == 0 && self.after > 0 }
    }
}

/// The REF and ALT windows, paired.
#[derive(Clone)]
struct Windows {
    pair: (Window, Window),
    /// Where the alleles differ, genomic [start, end): the trims' union, before
    /// repeat growth. A read spliced through it shows neither allele; a read
    /// spliced beside it is judged on windows rebuilt around it over the reference
    /// spliced along the read's junctions.
    event: (i64, i64),
    /// The given allele's own bases placed once, genomic [start, end) (for a pure
    /// indel its leftmost placement): right of them the ALT's bases are shifted by the
    /// length change.
    own: (i64, i64),
    /// Each other allele of the ALT's length at the locus (RJ-25), as the ALT's window
    /// over the same span: its own bases where it differs, the masks kept; and each of
    /// REF's length, as REF's window.
    others: Vec<Vec<u8>>,
    ref_others: Vec<Vec<u8>>,
    /// The bases held before masks left none telling the alleles apart and both
    /// windows read on to one that does (RJ-23); None when they did not.
    masked_hold: Option<usize>,
}

impl Windows {
    /// The windows as a read spliced at `skips` sees them. A splice in a window's
    /// flank ends the window at the exon edge on that side: the read's next bases
    /// come from the next exon, and it is anchored at the junction. Both alleles'
    /// windows are cut alike, and the bases held shrink with them. A splice past a
    /// window's right end ends its reading on there.
    fn cut_at(&self, skips: &[(i64, i64)]) -> Windows {
        let (ev_lo, ev_hi) = self.event;
        let exon_start = skips.iter().filter(|s| s.1 <= ev_lo).map(|s| s.1).max();
        let exon_end = skips.iter().filter(|s| s.0 >= ev_hi).map(|s| s.0).min();
        let cut = |w: &Window| -> (Window, usize, usize) {
            let (mut allele, mut left, mut right, mut before, mut after) =
                (w.allele.clone(), w.left, w.right, w.before, w.after);
            let (mut dl, mut dr) = (0, 0);
            if let Some(n1) = exon_start.filter(|&n1| n1 > ev_lo - before as i64) {
                dl = (n1 - (ev_lo - before as i64)) as usize;
                allele.drain(..dl);
                before -= dl;
                left = left.map(|_| n1);
            }
            if let Some(n0) = exon_end.filter(|&n0| n0 < ev_hi + after as i64) {
                dr = (ev_hi + after as i64 - n0) as usize;
                allele.truncate(allele.len() - dr);
                after -= dr;
                right = right.map(|_| n0);
            }
            let past = match exon_end.zip(right) {
                _ if dr > 0 => Vec::new(),
                Some((n0, r)) if n0 >= r => w.past[..usize::try_from(n0 - r).map_or(0, |k| k.min(w.past.len()))].to_vec(),
                _ => w.past.clone(),
            };
            let outer = match (w.tail, exon_start.zip(left), exon_end.zip(right)) {
                (false, _, _) if dl > 0 => Vec::new(),
                (true, _, _) if dr > 0 => Vec::new(),
                (false, Some((n1, l)), _) if n1 <= l => w.outer[w.outer.len() - usize::try_from(l - n1).map_or(0, |k| k.min(w.outer.len()))..].to_vec(),
                (true, _, Some((n0, r))) if n0 >= r => w.outer[..usize::try_from(n0 - r).map_or(0, |k| k.min(w.outer.len()))].to_vec(),
                _ => w.outer.clone(),
            };
            let event_at = (w.event_at.0.saturating_sub(dl).min(allele.len()), w.event_at.1.saturating_sub(dl).min(allele.len()));
            (Window { allele, hold: w.hold, tail: w.tail, left, right, before, after, past, outer, event_at }, dl, dr)
        };
        let ((mut r, dl, dr), (mut a, adl, adr)) = (cut(&self.pair.0), cut(&self.pair.1));
        let hold = r.hold.saturating_sub(if r.tail { dr } else { dl }).min(r.allele.len()).min(a.allele.len());
        r.hold = hold;
        a.hold = hold;
        // the other alleles share their allele's span: cut as it is
        let others = self.others.iter().map(|w| w[adl..w.len() - adr].to_vec()).collect();
        let ref_others = self.ref_others.iter().map(|w| w[dl..w.len() - dr].to_vec()).collect();
        Windows { pair: (r, a), event: self.event, own: self.own, others, ref_others, masked_hold: None }
    }

    /// The windows read on past their right end by the reference's next bases, both
    /// alleles continuing with them (the bases a read must hold unchanged). None for
    /// windows read from the right, or with no reference past them.
    fn read_on(&self) -> Option<Windows> {
        let (r, a) = &self.pair;
        let more = Some(r.past.as_slice()).filter(|m| !r.tail && !m.is_empty())?;
        let on = |w: &[u8]| [w, more].concat();
        let window = |w: &Window| Window {
            allele: on(&w.allele),
            hold: w.hold,
            tail: w.tail,
            left: w.left,
            right: w.right.map(|g| g + more.len() as i64),
            before: w.before,
            after: w.after + more.len(),
            past: Vec::new(),
            outer: w.outer.clone(),
            event_at: w.event_at,
        };
        Some(Windows {
            pair: (window(r), window(a)),
            event: self.event,
            own: self.own,
            others: self.others.iter().map(|o| on(o)).collect(),
            ref_others: self.ref_others.iter().map(|o| on(o)).collect(),
            masked_hold: self.masked_hold,
        })
    }

    /// The windows with a read's separate event (RJ-8) on both alleles: `ins` inserted
    /// before genome base `at`, or `del` reference bases deleted from it. The read is
    /// judged as either allele carrying it, holding as many more (fewer) bases. None
    /// unless it lies inside the windows past their anchors, clear of the given
    /// allele's own bases.
    fn with_separate(&self, at: i64, ins: &[u8], del: usize) -> Option<Windows> {
        let (r, a) = &self.pair;
        let (left, right) = (r.left?, r.right?);
        let end = at + del as i64;
        let clear = if del > 0 { end <= self.own.0 || at >= self.own.1 } else { at <= self.own.0 || at >= self.own.1 };
        // a deletion running past the windows' right end cuts them there (read from the
        // left; read from the right that end is the anchor)
        let past_end = end > right;
        if at <= left || at >= right || (past_end && r.tail) || !clear {
            return None;
        }
        // right of the given allele's own bases (past an insertion's junction)
        let right_of = if del > 0 || self.own.0 < self.own.1 { at >= self.own.1 } else { at > self.own.0 };
        let n = ins.len() as i64 - del as i64;
        let d = a.allele.len() as i64 - r.allele.len() as i64;
        let i = (at - left) as usize;
        let ia = if right_of { (i as i64 + d) as usize } else { i };
        let edit = |w: &[u8], k: usize| -> Vec<u8> {
            let mut out = w[..k.min(w.len())].to_vec();
            out.extend_from_slice(ins);
            out.extend_from_slice(&w[(k + del).min(w.len())..]);
            out
        };
        // Inside the bases held from the reading flank, a deletion holds as many fewer
        // (the same reference bases); an insertion's bases count as held (the event is
        // held before or after them, with its flank base).
        let held = |w: &Window, k: usize| if w.tail { k + del > w.allele.len() - w.hold } else { k < w.hold };
        let hold = if del > 0 && (held(r, i) || held(a, ia)) { r.hold.saturating_sub(del).max(1) } else { r.hold };
        // never more than either edited window holds
        let edited_len = |w: &Window, k: usize| k.min(w.allele.len()) + ins.len() + w.allele.len().saturating_sub(k + del);
        let hold = hold.min(edited_len(r, i)).min(edited_len(a, ia));
        let window = |w: &Window, k: usize| Window {
            allele: edit(&w.allele, k),
            hold,
            tail: w.tail,
            left: w.left,
            right: w.right,
            before: if right_of { w.before } else { usize::try_from(w.before as i64 + n).ok().unwrap_or(0) },
            after: if right_of { usize::try_from(w.after as i64 + n).ok().unwrap_or(0) } else { w.after },
            past: if past_end { Vec::new() } else { w.past.clone() },
            outer: w.outer.clone(),
            event_at: if k <= w.event_at.0 {
                let shift = |x: usize| usize::try_from(x as i64 + n).unwrap_or(0);
                (shift(w.event_at.0), shift(w.event_at.1))
            } else {
                w.event_at
            },
        };
        Some(Windows {
            pair: (window(r, i), window(a, ia)),
            event: self.event,
            own: self.own,
            others: self.others.iter().map(|o| edit(o, ia)).collect(),
            ref_others: self.ref_others.iter().map(|o| edit(o, i)).collect(),
            masked_hold: None,
        })
    }

    /// The reference span the windows cover.
    fn span(&self) -> (i64, i64) {
        let (r, a) = &self.pair;
        let lo = [r.left, a.left].into_iter().flatten().fold(i64::MAX, i64::min);
        let hi = [r.right, a.right].into_iter().flatten().fold(i64::MIN, i64::max);
        (lo, hi)
    }
}

/// A window read from the read: mismatching bases, masked bases (below min BQ,
/// or N; they match anything), whether an N sat among them, and the query span.
#[derive(Clone)]
pub(crate) struct Reading {
    pub(crate) mismatches: usize,
    pub(crate) masked: usize,
    pub(crate) had_n: bool,
    pub(crate) span: (usize, usize),
    /// Query positions of every base read: the held ones and the window's rest.
    pub(crate) read: (usize, usize),
    /// The read's evidence against this allele (log10): the weight of every
    /// base that mismatches it, by its quality, masked bases included.
    pub(crate) mismatch_weight: f64,
    /// Query position of the window's first base (it may lie before the read).
    pub(crate) start: i64,
    /// Query positions of the confident mismatches.
    pub(crate) misses: Vec<usize>,
}

/// A substitution a reading tolerates: the allele read (0 REF, 1 ALT), the reference
/// position its window is anchored at (its left end, or one past its right end when
/// read from the right), the base's offset from there, the read's base, and whether it
/// lies in the event (grown through repeats, or a base the input names) rather than its
/// flank.
pub(crate) type Tolerated = (u8, i64, i64, u8, bool);

/// The substitution `own`, an allele's reading of window `w` (`side`: 0 REF, 1 ALT),
/// tolerates beside the other allele's reading `other` of `wo`: its one confident
/// mismatch, which the other allele does not fit either (both mismatch the same read
/// base, costing both alike: it decides nothing, an error; operator, 2026-10-15). Past
/// the other allele's window (the shorter allele's ends first) its base is the
/// reference's next, which both alleles continue with. None with no mismatch, two or
/// more, or one telling the alleles apart.
fn tolerated(seq: &[u8], own: &Reading, other: &Reading, w: &Window, wo: &Window, side: u8) -> Option<Tolerated> {
    let &[q] = own.misses.as_slice() else { return None };
    let i = q as i64 - own.start;
    let b = seq[q].to_ascii_uppercase();
    let shared = match usize::try_from(q as i64 - other.start) {
        Ok(j) if j < wo.allele.len() => other.misses.contains(&q),
        Ok(j) if !wo.tail => wo.past.get(j - wo.allele.len()).is_some_and(|&h| !fits(b, h)),
        _ => false,
    };
    let in_event = usize::try_from(i).is_ok_and(|i| i >= w.event_at.0 && i < w.event_at.1);
    let anchor = if w.tail { w.right.map(|r| (r, i - w.allele.len() as i64)) } else { w.left.map(|l| (l, i)) };
    anchor.filter(|_| shared).map(|(at, off)| (side, at, off, b, in_event))
}

/// The substitution `own`'s one confident mismatch makes, keyed as [`tolerated`] keys
/// one, whether or not the other allele's reading shares it: read on past the windows,
/// [`judge`] forgives it there where the other reading misses three or more, the read
/// writes the allele's length, and it recurs on no other fragment.
fn forgiven(seq: &[u8], own: &Reading, w: &Window, side: u8) -> Option<Tolerated> {
    let &[q] = own.misses.as_slice() else { return None };
    let i = q as i64 - own.start;
    let in_event = usize::try_from(i).is_ok_and(|i| i >= w.event_at.0 && i < w.event_at.1);
    let anchor = if w.tail { w.right.map(|r| (r, i - w.allele.len() as i64)) } else { w.left.map(|l| (l, i)) };
    anchor.map(|(at, off)| (side, at, off, seq[q].to_ascii_uppercase(), in_event))
}

/// Fragments sharing a tolerated substitution in the event (grown through repeats):
/// two or more and [`GUARD_REF_SHARE`] of the fragments whose reads fit that allele make
/// it no error but another allele at the locus (RJ-25), which the scan from the left
/// flank may not place (reads starting inside a run). An error is a one-off. In the
/// flank a shared substitution is a phased neighbour (masked and named from three
/// fragments, RJ-23), never another allele: tolerated there.
pub(crate) const ERROR_RECURS_AT: usize = 2;

/// Calls back with each of a variant's reads its counts read and the qualities they
/// are read at.
pub(crate) type EachRead<'a> = Box<dyn Fn(&mut dyn FnMut(&Record, &[u8])) + 'a>;

/// A variant's recurring substitutions ([`ERROR_RECURS_AT`]), found from its reads on
/// first use: only a read that would tolerate one pays for them.
pub(crate) struct Tolerance<'a> {
    key: (i64, String, String),
    found: std::cell::OnceCell<std::collections::HashSet<Tolerated>>,
    reads: EachRead<'a>,
}

impl<'a> Tolerance<'a> {
    pub(crate) fn new(v: &Variant, reads: impl Fn(&mut dyn FnMut(&Record, &[u8])) + 'a) -> Self {
        let key = (v.pos, v.ref_allele.clone(), v.alt_allele.clone());
        Self { key, found: std::cell::OnceCell::new(), reads: Box::new(reads) }
    }

    /// Whether `k` recurs among `v`'s reads, judged under `rules` tolerating every
    /// substitution (false for another variant: a decomposed form keeps its own).
    fn recurs(&self, v: &Variant, min_baseq: u8, rules: &ReadRules, k: &Tolerated) -> bool {
        if (self.key.0, self.key.1.as_str(), self.key.2.as_str()) != (v.pos, v.ref_allele.as_str(), v.alt_allele.as_str()) {
            return false;
        }
        self.found
            .get_or_init(|| {
                use std::collections::HashSet;
                let all = ReadRules { tolerance: None, ..*rules };
                // per allele, the fragments whose reads fit it alone; per substitution,
                // those tolerating it
                let mut fit: [HashSet<Vec<u8>>; 2] = Default::default();
                let mut shared: HashMap<Tolerated, HashSet<Vec<u8>>> = HashMap::new();
                (self.reads)(&mut |record, quals| {
                    let Some(j) = judge(record, v, quals, min_baseq, &all) else { return };
                    match j.held {
                        Some((true, false)) => fit[0].insert(record.qname().to_vec()),
                        Some((false, true)) => fit[1].insert(record.qname().to_vec()),
                        _ => false,
                    };
                    if let Some(k) = j.tolerated {
                        shared.entry(k).or_default().insert(record.qname().to_vec());
                    }
                });
                shared
                    .into_iter()
                    .filter(|(k, frags)| {
                        frags.len() >= ERROR_RECURS_AT && frags.len() as f64 >= GUARD_REF_SHARE * fit[k.0 as usize].len() as f64
                    })
                    .map(|(k, _)| k)
                    .collect()
            })
            .contains(k)
    }
}

/// Classify a read at a complex variant by the exact-carrier rule, or None when
/// the rule cannot judge this variant: it needs both alleles non-empty and a
/// reference (prep's `event_ref`, else the `ref_context`) holding the REF allele
/// and the event with its flank. None whenever `windows` is None: no prepared
/// reference, a reference that does not hold the REF allele, or one that does not
/// hold the event with its flank (a contig end, or a repeat past prep's 16 kb
/// fetch cap). The caller then uses the previous classifier.
pub(crate) fn classify(
    record: &Record,
    variant: &Variant,
    quals: &[u8],
    min_baseq: u8,
    rules: &ReadRules,
) -> Option<ClassifyResult> {
    judge(record, variant, quals, min_baseq, rules).map(|j| j.result)
}

/// A read judged by the exact-carrier rule: its call, whether its readings fit REF
/// and the ALT (None: it holds no pair), and the substitution its fit tolerates.
struct Judged {
    result: ClassifyResult,
    held: Option<(bool, bool)>,
    tolerated: Option<Tolerated>,
}

impl From<ClassifyResult> for Judged {
    fn from(result: ClassifyResult) -> Self {
        Judged { result, held: None, tolerated: None }
    }
}

/// A read's REF and ALT readings with their placement, when it holds the pair, and
/// whether each fits its allele.
type Readings = (Option<(Reading, Reading, Placement)>, Option<(bool, bool)>);

/// [`classify`], with what the call rests on.
fn judge(record: &Record, variant: &Variant, quals: &[u8], min_baseq: u8, rules: &ReadRules) -> Option<Judged> {
    let layout = Layout::for_variant(variant, rules);
    let win = windows(variant, &layout)?;
    let seq = record.seq().as_bytes();
    if seq.is_empty() || quals.len() < seq.len() {
        // A record stored without its bases (SEQ '*') shows nothing here.
        return Some(ClassifyResult::no_coverage(ClassifyPhase::MaskedCompare).into());
    }
    // A pure indel's read whose first base is its anchoring flank base reads from that
    // flank only when it reads that base, the reference's (RJ-3): unread or another
    // base, it is masked and the read is placed from its other flank.
    let rw0 = &win.pair.0;
    let pure = matches!(super::window::allele_kind(&variant.ref_allele, &variant.alt_allele), Some(AlleleKind::Insertion | AlleleKind::Deletion));
    let unanchored = (pure && !rw0.tail)
        .then(|| rw0.left.and_then(|g| find_read_pos(record, g)))
        .flatten()
        .filter(|&q| readable_range(record, rules).is_some_and(|(first, _)| q == first))
        .filter(|&q| !(quals[q] >= min_baseq && fits(seq[q].to_ascii_uppercase(), rw0.allele[0])));
    let masked_quals;
    let quals = match unanchored {
        Some(q) => {
            masked_quals = [&quals[..q], &[0], &quals[q + 1..]].concat();
            &masked_quals[..]
        }
        None => quals,
    };
    let span = win.span();
    let guard = rules.guard.and_then(|g| g.get(variant));
    let own = own_event(record, &seq, quals, min_baseq, variant, &win, guard);
    let skips = rna::splice_junctions_in(record, span);
    // A junction past the windows that their reading on reaches: there the read's bases
    // are the next exon's, not the genome's continuation.
    let reach = prepared_span(variant).map_or_else(|| skips.clone(), |s| rna::splice_junctions_in(record, s));
    let (ev_lo, ev_hi) = win.event;
    // A junction that starts inside the event and runs past it (or ends inside it,
    // coming from before) leaves the read showing the event's first (last) bases
    // and then the next exon: judged on the haplotypes spliced at the event's edge.
    let enters = |&(n0, n1): &(i64, i64)| {
        (n0 > ev_lo && n0 < ev_hi && n1 >= ev_hi) || (n1 > ev_lo && n1 < ev_hi && n0 <= ev_lo)
    };
    if skips.iter().any(|j| j.0 < ev_hi && j.1 > ev_lo && (rules.reference.is_none() || !enters(j))) {
        // Spliced through the event: the read skips bases where the alleles
        // differ, so it shows neither (depth only).
        let mut r = ClassifyResult::neither(ClassifyPhase::CigarRecon);
        r.uninformative = true;
        return Some(r.into());
    }
    // Such a read is also REF spliced at its own junction (an alternative donor or
    // acceptor, or a splice site the event straddles): it can count only ALT, and
    // only when its bases favour ALT over that reading too.
    let entering = skips.iter().any(enters);
    if splice_beside_deletion(record, span) && matches!(super::window::allele_kind(&variant.ref_allele, &variant.alt_allele), Some(AlleleKind::Deletion)) {
        // A deletion written right beside the read's splice has no aligned flank: its
        // bases equal REF spliced at an acceptor further on (RJ-17): depth only.
        trace!("carrier: a deletion written beside its splice has no aligned flank → depth only");
        let mut r = ClassifyResult::neither(ClassifyPhase::CigarRecon);
        r.uninformative = true;
        return Some(r.into());
    }
    let spliced: Rc<SplicedWindows>;
    let cut;
    let separate;
    let at_edge;
    let identity = SpliceMap::default();
    let mut map = &identity;
    let base = if reach.is_empty() {
        &win
    } else if let Some(sw) = rules.reference.and_then(|r| {
        spliced_windows(variant, &rna::extract_splice_junctions(record), win.event, r, rules.spliced, &layout)
    }) {
        // Spliced inside the windows or past them, within their reading's reach: judged
        // on the reference spliced along its junctions, read on along its own mRNA.
        spliced = sw;
        map = &spliced.1;
        &spliced.0
    } else if entering {
        // No spliced windows to judge a junction entering the event on: depth only.
        trace!("carrier: spliced windows unavailable for a junction entering the event → depth only");
        let mut r = ClassifyResult::neither(ClassifyPhase::CigarRecon);
        r.uninformative = true;
        return Some(r.into());
    } else {
        // The far exon cannot be read, or the spliced reference no longer holds
        // the REF allele (its shared bases reach across the junction): judge the
        // read on the windows cut at its exon edges.
        trace!("carrier: spliced windows unavailable → windows cut at the exon edge");
        cut = win.cut_at(&reach);
        &cut
    };
    // A read with no junction inside the windows holds them where the genome does (one
    // past them only continues them along its next exon), and its own indel was placed
    // on them.
    let held_as_genome = skips.is_empty() && base.pair.0.left == win.pair.0.left;
    let win = if held_as_genome {
        // An unspliced read starting at an exon start that cuts the left flank, or
        // within the overhang an aligner writes into the intron as mismatches
        // (EDGE_SLACK), reads the windows from the right up to that edge, as a spliced
        // read does: on its mRNA the bases before the edge are the previous exon's, so
        // it reads from inside a run that crosses the junction. A pre-mRNA read through
        // the intron reads its bases as they are.
        let win = match layout.cut_left.filter(|&e| record.pos() >= e - EDGE_SLACK) {
            Some(e) => {
                at_edge = base.cut_at(&[(e, e)]);
                &at_edge
            }
            None => base,
        };
        // Another indel outside the discrimination window is a separate event (RJ-8):
        // the read is judged as either allele carrying it, so it keeps its call.
        match &own {
            OwnEvent::Separate { at, ins, del } => match win.with_separate(*at, ins, *del) {
                Some(w) => {
                    trace!("carrier: its own indel at {} lies outside the window → judged as either allele carrying it", at);
                    separate = w;
                    &separate
                }
                None => win,
            },
            _ => win,
        }
    } else {
        base
    };
    let Some(readable) = readable_range(record, rules) else {
        return Some(ClassifyResult::no_coverage(ClassifyPhase::MaskedCompare).into());
    };

    // An allele's reading fits it exactly, or but for one substitution the other
    // reading shares (an error: it decides nothing) that the variant's reads do not
    // share (then no error). Never an indel's footprint: a read whose own alignment
    // writes another length change inside the windows fits exactly or not.
    // An allele tolerates an error only where the read writes its length (none for REF,
    // the ALT's change for the ALT), or a separate event judged on both.
    let d = variant.alt_allele.len() as i64 - variant.ref_allele.len() as i64;
    let own_change: i64 = if matches!(own, OwnEvent::Given { .. }) { d } else { own_indels(record, span).iter().map(|o| o.1).sum() };
    let separate = matches!(own, OwnEvent::Separate { .. });
    // where the read's bases put the windows' first base, when they place its gap
    let left = match own {
        OwnEvent::Given { left } if held_as_genome => Some(left),
        _ => None,
    };
    let tolerate = |side: u8| {
        separate || own_change == if side == 0 { 0 } else { d }
    };
    let kept = |own: &Reading, other: &Reading, w: &Window, wo: &Window, side: u8| {
        tolerated(&seq, own, other, w, wo, side)
            .filter(|k| tolerate(side) && !(k.4 && rules.tolerance.is_some_and(|t| t.recurs(variant, min_baseq, rules, k))))
    };
    let fit = |own: &Reading, other: &Reading, w: &Window, wo: &Window, side: u8| {
        own.mismatches == 0 || kept(own, other, w, wo, side).is_some()
    };

    // The read's REF and ALT readings, when it holds the pair, and whether they fit
    // (matches REF, matches ALT).
    let read_bases = |w: &Windows, readable: (usize, usize), seq: &[u8]| {
        let placed = read_pair(record, seq, quals, min_baseq, &w.pair.0, &w.pair.1, map, readable, &fit, left, unanchored.is_some());
        let (rw, aw) = &w.pair;
        // An error both readings show alike decides nothing: where they miss exactly the
        // same bases (no confidently read base tells the alleles apart) and the read's own
        // bases place no length change, the allele tolerating it for the length the read
        // is written with (its alignment, not its bases) fits no better than the other.
        let held = placed.as_ref().map(|(r, a, _)| match (fit(r, a, rw, aw, 0), fit(a, r, aw, rw, 1)) {
            (f, g) if f != g && matches!(own, OwnEvent::None) && r.mismatches > 0 && r.misses == a.misses => (true, true),
            fits => fits,
        });
        (placed, held)
    };
    let read = |w: &Windows, readable: (usize, usize)| read_bases(w, readable, &seq);
    // The read's bases with the errors both readings show alike unread (N): read on past
    // the windows, such an error decides nothing either, nor spends the one error an
    // allele tolerates.
    let shared_unread = |(placed, _): &Readings| -> Option<Vec<u8>> {
        let (r, _, _) = placed.as_ref().filter(|(r, a, _)| matches!(own, OwnEvent::None) && r.mismatches > 0 && r.misses == a.misses)?;
        let mut bases = seq.clone();
        r.misses.iter().for_each(|&q| bases[q] = WILD);
        Some(bases)
    };
    // Readings that fit both alleles, or miss both alike (masked bases, or the shorter
    // allele's window ending before the read's bases show its length): the read is read
    // on past the windows, both alleles continuing with the reference, as the windows
    // read on where the guard's masks hide the alleles' difference.
    let undecided = |(placed, held): &Readings| match held {
        Some((true, true)) => true,
        Some((false, false)) => placed.as_ref().is_some_and(|(r, a, _)| r.mismatches == a.mismatches),
        _ => false,
    };
    // A read placed from the reading flank, fitting the ALT of a length change, whose
    // own bases place no length change there (it writes none, or one its bases do not
    // place: a substitution inside the windows can spell the shorter allele's bases) or
    // another one is read on too: ALT stands only where its further bases fit the ALT
    // better than REF, as a carrier's do however the aligner wrote its gap.
    let unplaced = |(placed, held): &Readings| {
        d != 0
            && matches!(own, OwnEvent::None | OwnEvent::Inside)
            && *held == Some((false, true))
            && placed.as_ref().is_some_and(|(_, _, p)| *p == Placement::Flank)
    };
    // A pure indel's fit resting on fewer bases than tell the alleles from the others
    // nearest them is confirmed past the windows, where another allele's bases, shifted
    // by its other length, show (one base past the tract may be its extra base, the
    // flank's twin): the ALT where the read's own bases place no change or another one,
    // REF where it writes none, and either allele's fit on a tolerated mismatch. From
    // its first mismatch, or the window's end, the read must not fit the allele shifted
    // better than as it stands (`shifted_better`); a phased change beside the event
    // mismatches at its own bases under any shift and keeps the call.
    let confirms = |held: Option<(bool, bool)>, placed: &Option<(Reading, Reading, Placement)>| {
        let Some((r, a, Placement::Flank)) = placed else { return false };
        pure && match held {
            Some((false, true)) => a.mismatches > 0 || matches!(own, OwnEvent::None | OwnEvent::Inside),
            Some((true, false)) => r.mismatches > 0 || matches!(own, OwnEvent::None),
            _ => false,
        }
    };
    let shifted = |w: &Windows, placed: &Option<(Reading, Reading, Placement)>, alt: bool, end: usize| {
        let Some((r, a, _)) = placed else { return false };
        let (reading, win) = if alt { (a, &w.pair.1) } else { (r, &w.pair.0) };
        let from = reading.misses.first().map_or(reading.start + win.allele.len() as i64, |&q| q as i64);
        let (Ok(q), Ok(k)) = (usize::try_from(from), usize::try_from(from - reading.start)) else { return false };
        let expected: Vec<u8> = win.allele.iter().chain(&win.past).skip(k).copied().collect();
        !win.tail && shifted_better(&seq, quals, min_baseq, (q, end), &expected)
    };
    let longer;
    let mut past_error = None;
    let mut judged = read(win, readable);
    // A read short of the base the masks made both windows read on to is still decided
    // by its own bases where, held as far as before, they fit one allele alone (a REF
    // read showing the reference's base where the masks take either): the masks leave
    // no base telling the alleles apart for every read, not for every read's bases.
    let short;
    let win = match win.masked_hold.filter(|_| judged.1.is_none()) {
        Some(h) => {
            let mut w = win.clone();
            (w.pair.0.hold, w.pair.1.hold, w.masked_hold) = (h, h, None);
            let s = read(&w, readable);
            // an exact fit: one error must not decide what the masks left undecided
            let exact = s.0.as_ref().is_some_and(|(r, a, _)| match s.1 {
                Some((true, false)) => r.mismatches == 0,
                Some((false, true)) => a.mismatches == 0,
                _ => false,
            });
            if exact {
                trace!("carrier: short of the masks' read-on, its own bases fit one allele → judged there");
                judged = s;
                short = w;
                &short
            } else {
                win
            }
        }
        None => win,
    };
    // Read on past the windows only as far as the read's next own indel past them (one
    // written on their edge is the read's own event): its bases beyond are shifted by
    // that separate event. Past them as the read is judged: on its spliced reference an
    // indel written after the junction the windows run through is at the event.
    let past_end = map.to_genome(win.span().1 - 1) + 1;
    let on = (readable.0, next_indel_query(record, past_end.max(span.1 + 1)).map_or(readable.1, |q| q.min(readable.1)));
    let win = if undecided(&judged) {
        match win.read_on() {
            Some(w) => {
                longer = w;
                trace!("carrier: its readings do not tell the alleles apart → read on past the windows");
                judged = match shared_unread(&judged) {
                    Some(bases) => {
                        let mut j = read_bases(&longer, on, &bases);
                        // the bases set unread are the read's own, no N and not masked
                        if let Some((r, a, _)) = j.0.as_mut() {
                            for x in [r, a] {
                                let (q0, q1) = x.read;
                                let own = seq.get(q0..q1).unwrap_or_default();
                                x.had_n = own.iter().any(|&b| b.to_ascii_uppercase() == WILD);
                                x.masked = x.masked.saturating_sub((q0..q1).filter(|&q| bases[q] == WILD && seq[q].to_ascii_uppercase() != WILD).count());
                            }
                        }
                        j
                    }
                    None => read(&longer, on),
                };
                // Past the windows the alleles differ at most bases (the other allele is
                // shifted by its length): one error there the other allele's reading does not
                // share is forgiven where that reading misses three or more, the read writes
                // the allele's length (as within them), and it is an error, a one-off
                // (recurring, it is another allele's base, [`ERROR_RECURS_AT`]). Within the
                // windows only a shared error is tolerated: there one base may be all that
                // tells another allele of the same length.
                if let (Some((r, a, _)), Some((false, false))) = (&judged.0, judged.1) {
                    let forgive = |own: &Reading, other: &Reading, w: &Window, held: usize, side: u8| {
                        let &[q] = own.misses.as_slice() else { return None };
                        let past = usize::try_from(q as i64 - own.start).is_ok_and(|i| i >= held);
                        (past && tolerate(side) && other.mismatches >= 3)
                            .then(|| forgiven(&seq, own, w, side))
                            .flatten()
                            .filter(|k| !rules.tolerance.is_some_and(|t| t.recurs(variant, min_baseq, rules, k)))
                    };
                    if let Some(k) = forgive(r, a, &longer.pair.0, win.pair.0.allele.len(), 0) {
                        (judged.1, past_error) = (Some((true, false)), Some(k));
                    } else if let Some(k) = forgive(a, r, &longer.pair.1, win.pair.1.allele.len(), 1) {
                        (judged.1, past_error) = (Some((false, true)), Some(k));
                    }
                }
                &longer
            }
            None => win,
        }
    } else {
        if unplaced(&judged) {
            // only to confirm the call: the further bases may hold errors past the windows
            let better = win.read_on().map(|w| read(&w, on).0).is_none_or(|ext| {
                ext.as_ref().is_some_and(|(r, a, _)| a.mismatch_weight < r.mismatch_weight)
            });
            if !better {
                trace!("carrier: its bases place no length change, or another, and past the windows fit the ALT no better than REF → neither");
                judged.1 = judged.1.map(|_| (false, false));
            }
        }
        if confirms(judged.1, &judged.0) && shifted(win, &judged.0, judged.1 == Some((false, true)), on.1) {
            trace!("carrier: past its fit its bases fit the allele shifted by another length → neither");
            judged.1 = judged.1.map(|_| (false, false));
        }
        win
    };
    let (rw, aw) = &win.pair;
    let (placed, mut held) = judged;
    if matches!(own, OwnEvent::Given { .. }) {
        // Its bases place a change of the ALT's length inside the window: never REF,
        // however little of it the REF window reaches.
        held = held.map(|(_, alt)| (false, alt));
    }
    // A read the aligner placed a unit or a copy off in a repeat: ungapped, its bases fit
    // an allele over the windows where it was written and miss the reference beside
    // them, or its own indels absorb the units' differences inside the windows. Wherever
    // its alignment is not clean (indels the windows do not place, or, writing none, two
    // confident changes from the reference) it is read whole, without its alignment
    // ([`whole_fit`]): the allele its bases fit there decides, and a read fitting the
    // longer allele only from inside its event is depth only. A read whose own bases place
    // its indel outside the windows (a separate event, RJ-8), or the ALT's change inside
    // them, keeps its judgment; so does a read whose bases read whole decide nothing, and
    // a spliced read (its bases are not one stretch of the genome).
    let spliced = record.cigar().iter().any(|op| matches!(op, Cigar::RefSkip(_)));
    let built = if spliced || !matches!(own, OwnEvent::Inside | OwnEvent::None) {
        None
    } else if let Some(cg) = rules.guard {
        Some(cg.whole(variant, || whole_haps(variant, &layout)))
    } else {
        Some(std::sync::Arc::new(whole_haps(variant, &layout)))
    };
    let haps = built.as_deref().and_then(Option::as_ref);
    let unclean = haps.filter(|w| matches!(own, OwnEvent::Inside) || aligned_changes(record, &seq, quals, min_baseq, w) >= 2);
    let whole = unclean.and_then(|w| whole_fit(w, &seq, quals, min_baseq, readable));
    match whole {
        Some(Some(alt)) => {
            if held != Some((!alt, alt)) {
                trace!("carrier: read whole, without its own alignment, its bases fit {} → that allele", if alt { "the ALT" } else { "REF" });
            }
            held = Some((!alt, alt));
        }
        Some(None) => {
            // as far as its bases go it starts inside the longer allele's event, showing
            // only the other junction (RJ-23): depth only
            trace!("carrier: read whole, without its own alignment, it fits the longer allele only from inside its event → depth only");
            held = None;
        }
        None => {}
    }
    // a verdict read whole stands for the read's own indels too
    let own = if whole.is_some() { OwnEvent::None } else { own };

    let placement = placed.as_ref().map(|(_, _, p)| *p);
    let reading = placed.map(|(r, a, _)| (r, a));
    let leaned = reading.as_ref().and_then(|(r, a)| match held {
        Some((true, false)) if r.mismatches > 0 => kept(r, a, rw, aw, 0),
        Some((false, true)) if a.mismatches > 0 => kept(a, r, aw, rw, 1),
        _ => None,
    });
    let (mm_ref, mm_alt) = reading.as_ref().map_or((0, 0), |(r, a)| (r.mismatches, a.mismatches));
    let alt_masked = reading.as_ref().map_or(0, |(_, a)| a.masked);
    let had_n = reading.as_ref().is_some_and(|(r, a)| r.had_n || a.had_n);
    // The read's evidence for ALT over REF (log10), every base weighed by its
    // quality: a low-quality base counts for little instead of fitting anything.
    // Bases matching both alleles cancel: the evidence is what mismatches REF less
    // what mismatches ALT, so the two readings need not cover the same bases. A
    // tolerated error costs both alike: past the other allele's window (the shorter
    // allele's ends first) the reference's continuation mismatches it too.
    let past_other = |own: &Reading, other: &Reading, wo: &Window| {
        own.misses
            .first()
            .filter(|&&q| usize::try_from(q as i64 - other.start).is_ok_and(|j| j >= wo.allele.len()))
            .map_or(0.0, |&q| one_base_evidence(quals[q]))
    };
    let evidence = reading.as_ref().map_or(0.0, |(r, a)| {
        let e = r.mismatch_weight - a.mismatch_weight;
        match (leaned.is_some(), held) {
            (true, Some((false, true))) => e + past_other(a, r, rw),
            (true, Some((true, false))) => e - past_other(r, a, aw),
            _ => e,
        }
    });
    let alt_weight = reading.as_ref().map_or(0.0, |(_, a)| a.mismatch_weight);
    let alt_spans: Vec<(usize, usize)> = reading.iter().map(|(_, a)| a.span).collect();
    let read_spans: Vec<(usize, usize)> = reading.iter().flat_map(|(r, a)| [r.span, a.span]).collect();
    let qual_bases: Vec<u8> = read_spans.iter().flat_map(|&(lo, hi)| quals[lo..hi].iter().copied()).collect();
    // The call's quality is its held bases'; where none is read at min BQ (masked at a
    // junction, the read decided reading on), every base its readings read: a call
    // the fragment layer takes as no evidence would leave a decided read uncounted
    // there.
    let qual = match median_qual(&qual_bases, min_baseq) {
        0 => {
            let read: Vec<u8> = reading.iter().flat_map(|(r, a)| [r.read, a.read]).flat_map(|(lo, hi)| quals[lo..hi].iter().copied()).collect();
            median_qual(&read, min_baseq)
        }
        q => q,
    };
    let every = |want: (bool, bool)| held == Some(want);
    let mut result = if let Some(Placement::Unfit { alt }) = placement {
        // Placed as the allele its bases favour, which they do not fit: another
        // allele, partial evidence when that is the ALT.
        if alt {
            ClassifyResult::neither_with_nearby(qual, ClassifyPhase::MaskedCompare)
        } else {
            ClassifyResult::neither(ClassifyPhase::MaskedCompare)
        }
    } else if every((false, true)) && whole.is_none() && evidence < one_base_evidence(min_baseq) {
        // Its clearly read bases fit the ALT, but weighed by quality its bases do not
        // favour ALT over REF by as much as one base read at the minimum quality: a
        // REF molecule read poorly around the event, with one error fitting the ALT.
        trace!("carrier: fits the ALT but its evidence {:.2} is below one base's → neither", evidence);
        ClassifyResult::neither(ClassifyPhase::MaskedCompare)
    } else if every((false, true))
        && whole.is_none()
        && !reading.as_ref().is_some_and(|(_, a)| tells_others(&seq, quals, min_baseq, aw, &win.others, a.start, readable))
    {
        // Another allele of the ALT's length is at the locus and the read reads no
        // base telling the two apart (RJ-25): it is that allele or the ALT, partial
        // evidence.
        trace!("carrier: fits the ALT but reads no base telling it from another allele of its length → neither, partial");
        ClassifyResult::neither_with_nearby(qual, ClassifyPhase::MaskedCompare)
    } else if every((false, true)) {
        ClassifyResult::is_alt(qual, ClassifyPhase::MaskedCompare)
    } else if every((true, false)) && matches!(own, OwnEvent::Inside) {
        // Its bases fit REF over the windows, but its own alignment changes the length
        // inside them with an indel no background recurs as, and its bases place that
        // indel inside the windows: another indel inside the window (RJ-7) the windows
        // do not reach (an insertion whose first bases repeat the flank) - another
        // allele, partial evidence. Bases decide: a read whose bases still fit REF
        // (ending inside a run) stays REF, and one whose indel its bases place only
        // outside the windows is a separate event (RJ-8).
        trace!("carrier: fits REF but its bases place another indel inside the windows → neither, partial");
        ClassifyResult::neither_with_nearby(qual, ClassifyPhase::MaskedCompare)
    } else if every((true, false))
        && whole.is_none()
        && !reading.as_ref().is_some_and(|(r, _)| tells_others(&seq, quals, min_baseq, rw, &win.ref_others, r.start, readable))
    {
        // REF alike (RJ-25): another allele of REF's length is at the locus and the read
        // reads no base telling the two apart: that allele or REF, no evidence for ALT.
        trace!("carrier: fits REF but reads no base telling it from another allele of its length → neither");
        ClassifyResult { ref_untold: true, ..ClassifyResult::neither(ClassifyPhase::MaskedCompare) }
    } else if every((true, false)) {
        ClassifyResult::is_ref(qual, ClassifyPhase::MaskedCompare)
    } else if held.is_some() && mm_alt < mm_ref {
        // Another allele, closer to ALT: partial evidence.
        ClassifyResult::neither_with_nearby(qual, ClassifyPhase::MaskedCompare)
    } else {
        ClassifyResult::neither(ClassifyPhase::MaskedCompare)
    };
    // A read holding no pair cannot show the allele: depth only (as a pure-indel
    // read that ends inside its tract), and no mFSD class.
    result.uninformative = held.is_none();
    let given_other = matches!(own, OwnEvent::Given { .. }) && held.is_some();
    if !result.is_ref && !result.is_alt && !result.has_nearby_evidence && (matches!(own, OwnEvent::Inside) || given_other) {
        // Another indel inside the window, placed there by the read's own bases (RJ-7;
        // a deleted anchor, RJ-15), whether or not the read holds the windows; or one
        // of the ALT's length a read holding them does not spell as the ALT (RJ-21):
        // another allele, partial evidence.
        trace!("carrier: its bases place another indel inside the windows → partial evidence");
        result.has_nearby_evidence = true;
    }
    if entering && !result.uninformative {
        // Against REF spliced where the read splices (its own alignment), over the
        // bases its ALT reading read: ALT needs one minimum-quality base's evidence
        // there too; anything else is depth only.
        let own = rules
            .reference
            .and_then(|r| own_alignment_weight(record, &seq, quals, min_baseq, &alt_spans, &variant.chrom, r));
        let alt_stands = result.is_alt && own.is_some_and(|w| w - alt_weight >= one_base_evidence(min_baseq));
        if !alt_stands {
            trace!("carrier: junction entering the event, bases fit REF spliced there ({:?}) → depth only", own);
            result = ClassifyResult::neither(ClassifyPhase::MaskedCompare);
            result.uninformative = true;
        }
    }
    result.has_n_base = had_n;
    // An MNP ALT read with every window base read unmasked shows the whole
    // haplotype, as a fully read block does on the base-by-base path.
    result.mnp_confirmed =
        result.is_alt && alt_masked == 0 && variant.ref_allele.len() == variant.alt_allele.len();
    // A read decided from its own bases may count although its aligned span stops
    // short of the variant, when its allele lies in soft-clipped bases (a window it
    // was decided on reads them) and every base read lies inside its fragment: past
    // the fragment end a clip is adapter. An unclipped read that misses the variant
    // position (one starting inside a long deletion) is not admitted: no ALT read
    // can start there.
    result.clip_admissible = (result.is_ref || result.is_alt)
        && aligned_query_range(record).is_some_and(|(first, after)| read_spans.iter().any(|&(a, b)| a < first || b > after))
        && fragment_query_span(record).is_some_and(|(lo, hi)| read_spans.iter().all(|&(a, b)| a >= lo && b <= hi));
    // Windows read from the right flank are held from there by a REF and an ALT
    // molecule alike, so one starting past the variant position (inside the
    // longer allele) is judged as the other allele's molecule starting there is:
    // it counts without overlapping the variant position.
    // A read placed from the other flank through aligned bases starts past the
    // variant position as aligned, while the same molecule of the other allele would
    // span it: it counts as that one would.
    result.read_from_right =
        (result.is_ref || result.is_alt) && (rw.reads_from_right() || placement == Some(Placement::ThroughAligned));
    Some(Judged { result, held, tolerated: leaned.or(past_error) })
}

/// Whether the read, with an allele's window `w` starting at query position `start`,
/// reads for every other allele's window a base telling that allele from it: read
/// inside `[first, after)` at or above `min_baseq`, fitting the allele's window base
/// and not the other's (a masked base takes the bases of its code). True with no other
/// allele.
fn tells_others(
    seq: &[u8],
    quals: &[u8],
    min_baseq: u8,
    w: &Window,
    others: &[Vec<u8>],
    start: i64,
    (first, after): (usize, usize),
) -> bool {
    // the read's confident base at window position i
    let read = |i: usize| {
        let q = usize::try_from(start + i as i64).ok().filter(|&q| q >= first && q < after)?;
        (*quals.get(q)? >= min_baseq).then(|| seq.get(q).map(|b| b.to_ascii_uppercase()))?
    };
    others.iter().all(|other| {
        other.iter().zip(&w.allele).enumerate().any(|(i, (&o, &a))| read(i).is_some_and(|b| fits(b, a) && !fits(b, o)))
    })
}

/// Whether a read, as aligned, reads for every other allele of `guard` beside REF or
/// the ALT of `v` (`side`) a base that tells that allele from it: a position where
/// their haplotypes differ past the guard's masks (a background base no difference, a
/// third base still one), read as the allele's base, or its background, at or above
/// `min_baseq` (RJ-25). The MNP check's rule for every read it calls REF or ALT
/// (it compares only the bases where REF and ALT differ, so a read spelling another
/// allele where the ALT keeps the reference's base would pass it); the block is
/// aligned without an indel, so the haplotypes and the read share coordinates. True
/// with no other allele it can place.
pub(crate) fn tells_others_aligned(record: &Record, quals: &[u8], min_baseq: u8, v: &Variant, guard: &Guard, side: Side) -> bool {
    let others = if side == Side::Ref { &guard.ref_others } else { &guard.others };
    if others.is_empty() {
        return true;
    }
    let Some((start, reference)) = prepared_reference(v) else { return true };
    let own = if side == Side::Ref {
        Some(reference.clone())
    } else {
        usize::try_from(v.pos - start).ok().and_then(|p| haplotype(&reference, p, &upper(&v.ref_allele), &upper(&v.alt_allele)))
    };
    let Some(own) = own else { return true };
    let seq = record.seq().as_bytes();
    // the bases the guard masks at reference offset i to (0: none)
    let set_at = |i: usize| {
        let p = start + i as i64;
        if guard.masked.contains(&p) {
            guard.background.iter().find(|&&(q, _)| q == p).map_or(ANY, |&(_, s)| s)
        } else {
            0
        }
    };
    let reads_own = |i: usize, set: Bases| {
        find_read_pos(record, start + i as i64).is_some_and(|q| {
            seq.get(q).is_some_and(|&b| {
                let b = b.to_ascii_uppercase();
                bit(b) != 0 && class(set, b) == class(set, own[i])
            }) && quals.get(q).is_some_and(|&b| b >= min_baseq)
        })
    };
    others.iter().all(|o| {
        let Some(h) = o.on(start, &reference).filter(|h| h.len() == own.len()) else {
            trace!("carrier: another allele at {}:{} does not sit on the reference → not read", v.chrom, o.pos + 1);
            return true;
        };
        let (a, b) = o.span();
        let lo = usize::try_from(a.min(v.pos) - start).unwrap_or(0);
        let hi = usize::try_from(b.max(v.pos + v.ref_allele.len() as i64) - start).unwrap_or(0).min(own.len());
        (lo..hi).any(|i| {
            let set = set_at(i);
            class(set, h[i]) != class(set, own[i]) && reads_own(i, set)
        })
    })
}

/// The read's evidence against REF as it is aligned (spliced at its own
/// junctions), over the query `spans`: the weight of every base aligned to a
/// different reference base (by its quality, N excepted), of every inserted base,
/// and of one minimum-quality base per deletion. None when the reference cannot
/// be read there.
fn own_alignment_weight(
    record: &Record,
    seq: &[u8],
    quals: &[u8],
    min_baseq: u8,
    spans: &[(usize, usize)],
    chrom: &str,
    reference: &CachedFasta,
) -> Option<f64> {
    let in_span = |q: usize| spans.iter().any(|&(a, b)| q >= a && q < b);
    let (mut qp, mut rp, mut weight) = (0usize, record.pos(), 0.0f64);
    for op in record.cigar().iter() {
        match op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                let n = *n as usize;
                let (lo, hi) = (qp.max(spans.iter().map(|s| s.0).min()?), (qp + n).min(spans.iter().map(|s| s.1).max()?));
                if lo < hi {
                    let bases = reference.bases(chrom, rp + (lo - qp) as i64, rp + (hi - qp) as i64)?;
                    for (k, q) in (lo..hi).enumerate() {
                        let b = seq[q].to_ascii_uppercase();
                        if in_span(q) && b != WILD && b != bases[k] {
                            weight += one_base_evidence(quals[q]);
                        }
                    }
                }
                qp += n;
                rp += n as i64;
            }
            Cigar::Ins(n) => {
                for q in qp..qp + *n as usize {
                    if in_span(q) && seq[q].to_ascii_uppercase() != WILD {
                        weight += one_base_evidence(quals[q]);
                    }
                }
                qp += *n as usize;
            }
            Cigar::Del(n) => {
                if in_span(qp) && qp > 0 && in_span(qp - 1) {
                    weight += one_base_evidence(min_baseq);
                }
                rp += *n as i64;
            }
            Cigar::RefSkip(n) => rp += *n as i64,
            Cigar::SoftClip(n) => qp += *n as usize,
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    Some(weight)
}

/// Positions of a reference spliced along a read's junctions, mapped back to the
/// genome. Spliced space is a run of genome segments, the read's exons: the event
/// itself, then on each side the genome up to the read's next junction, the exon
/// after it, and so on.
#[derive(Default)]
struct SpliceMap {
    /// (spliced start, genome start, length), in spliced order.
    segments: Vec<(i64, i64, i64)>,
}

impl SpliceMap {
    fn to_genome(&self, g: i64) -> i64 {
        self.segments
            .iter()
            .find(|&&(s, _, len)| g >= s && g < s + len)
            .map_or(g, |&(s, gs, _)| gs + (g - s))
    }

    /// The spliced-space position of genome position `g`, when a segment holds it.
    fn to_spliced(&self, g: i64) -> Option<i64> {
        self.segments.iter().find(|&&(_, gs, len)| g >= gs && g < gs + len).map(|&(s, gs, _)| s + (g - gs))
    }
}

/// The windows a spliced read is judged on: built by [`windows`] over the
/// variant's reference spliced along the read's junctions on each side of the
/// event (every junction the windows reach, so a short exon between two is read
/// through; a junction entering the event joins the next exon at the event's
/// edge), the far exons read from `reference`, with the map from spliced space
/// back to the genome. Reads with the same splice segments share them through
/// `cache`. None when a far exon cannot be read or the spliced reference does not
/// hold the REF allele and the event's flank.
fn spliced_windows(
    v: &Variant,
    junctions: &[(i64, i64)],
    (ev_lo, ev_hi): (i64, i64),
    reference: &CachedFasta,
    cache: Option<&SplicedCache>,
    layout: &Layout,
) -> Option<Rc<SplicedWindows>> {
    let (start, context) = prepared_reference(v)?;
    let end = start + context.len() as i64;
    let segments = splice_segments((start, end), junctions, (ev_lo, ev_hi));
    let key = cache.map(|_| (v.pos, v.ref_allele.clone(), v.alt_allele.clone(), segments.clone()));
    if let (Some(c), Some(k)) = (cache, key.as_ref()) {
        if let Some(hit) = c.0.borrow().get(k) {
            return hit.clone();
        }
    }
    let built = (|| {
        let mut spliced = Vec::with_capacity(context.len());
        for &(_, gs, len) in &segments {
            if gs >= start && gs + len <= end {
                spliced.extend_from_slice(&context[(gs - start) as usize..(gs - start + len) as usize]);
            } else {
                spliced.extend(reference.bases(&v.chrom, gs, gs + len)?);
            }
        }
        let mut on_spliced = v.clone();
        // A pure indel's anchor is no part of its change: on the spliced reference it
        // is whatever base precedes the change there (an intronic anchor spliced out).
        let pure = matches!(super::window::allele_kind(&v.ref_allele, &v.alt_allele), Some(AlleleKind::Insertion | AlleleKind::Deletion));
        if let Some(&b) = usize::try_from(v.pos - start).ok().and_then(|p| spliced.get(p)).filter(|_| pure) {
            let anchored = |a: &str| format!("{}{}", b as char, &a[1..]);
            on_spliced.ref_allele = anchored(&v.ref_allele);
            on_spliced.alt_allele = anchored(&v.alt_allele);
        }
        on_spliced.event_ref = Some((start, String::from_utf8(spliced).ok()?));
        on_spliced.ref_context = None;
        // masked positions are genome ones: each kept where a segment holds it
        let map = SpliceMap { segments };
        let masked = layout.masked.iter().filter_map(|&p| map.to_spliced(p)).collect();
        let background = layout.background.iter().filter_map(|&(p, b)| map.to_spliced(p).map(|q| (q, b))).collect();
        // other alleles are genome ones too: each kept where its whole span is where
        // it sits in the genome (no junction of the read crosses it)
        let same = |p: i64| map.to_spliced(p) == Some(p);
        let in_place = |alleles: &[OtherAllele]| -> Vec<OtherAllele> {
            alleles
                .iter()
                .filter(|o| {
                    let (a, b) = o.span();
                    same(a) && same(b - 1)
                })
                .cloned()
                .collect()
        };
        // An exon edge cutting the left flank that the read splices at is none on its
        // spliced reference: the read holds the exon it comes from there, so its windows
        // read from the left (and on past the event, where bases at the junction are
        // too poor to decide).
        let own_edge = grown_event(v).is_some_and(|(lo, _)| starts_cut_left(junctions.iter().map(|j| j.1), lo, FLANK as i64));
        let spliced_layout = Layout {
            masked,
            background,
            others: in_place(&layout.others),
            ref_others: in_place(&layout.ref_others),
            from_right: layout.from_right && !own_edge,
            cut_left: layout.cut_left.filter(|_| !own_edge),
        };
        Some(Rc::new((windows(&on_spliced, &spliced_layout)?, map)))
    })();
    if let (Some(c), Some(k)) = (cache, key) {
        c.0.borrow_mut().insert(k, built.clone());
    }
    built
}

/// The segments `(spliced start, genome start, length)` of a reference spliced
/// along a read's junctions over `[start, end)`: the event `[ev_lo, ev_hi)` as it
/// is, then on each side the genome up to the next junction, the exon after it,
/// and so on. A junction starting inside the event (ending inside it) joins the
/// next (previous) exon at the event's edge.
fn splice_segments((start, end): (i64, i64), junctions: &[(i64, i64)], (ev_lo, ev_hi): (i64, i64)) -> Vec<(i64, i64, i64)> {
    // Right of the event: the genome from ev_hi, jumping each junction past it.
    let mut right = Vec::new();
    let into_right = junctions.iter().find(|j| j.0 > ev_lo && j.0 < ev_hi && j.1 >= ev_hi).map(|j| j.1);
    let (mut sp, mut g) = (ev_hi, into_right.unwrap_or(ev_hi));
    let mut after: Vec<(i64, i64)> = junctions.iter().filter(|j| j.0 >= ev_hi).copied().collect();
    after.sort_unstable();
    for (n0, n1) in after {
        if sp >= end {
            break;
        }
        if n0 < g {
            continue;
        }
        let len = (n0 - g).min(end - sp);
        if len > 0 {
            right.push((sp, g, len));
        }
        sp += n0 - g;
        g = n1;
    }
    if sp < end {
        right.push((sp, g, end - sp));
    }
    // Left of the event: the genome back from ev_lo, jumping each junction before it.
    let mut left = Vec::new();
    let into_left = junctions.iter().find(|j| j.1 > ev_lo && j.1 < ev_hi && j.0 <= ev_lo).map(|j| j.0);
    let (mut sp, mut g) = (ev_lo, into_left.unwrap_or(ev_lo));
    let mut before: Vec<(i64, i64)> = junctions.iter().filter(|j| j.1 <= ev_lo).copied().collect();
    before.sort_unstable_by(|a, b| b.cmp(a));
    for (n0, n1) in before {
        if sp <= start {
            break;
        }
        if n1 > g {
            continue;
        }
        let len = (g - n1).min(sp - start);
        if len > 0 {
            left.push((sp - len, g - len, len));
        }
        sp -= g - n1;
        g = n0;
    }
    if sp > start {
        left.push((start, g - (sp - start), sp - start));
    }
    left.reverse();
    left.into_iter().chain(std::iter::once((ev_lo, ev_lo, ev_hi - ev_lo))).chain(right).collect()
}

/// Query positions [first, after) of the read's aligned bases: its soft clips lie
/// before `first` and from `after` on. None for a read with no aligned base.
pub(crate) fn aligned_query_range(record: &Record) -> Option<(usize, usize)> {
    let (lead, tail) = soft_clips(record);
    let (lead, tail) = (lead as usize, tail as usize);
    let len = record.seq_len();
    (lead + tail < len).then_some((lead, len - tail))
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

/// Whether this rule judges the variant: every substitution-bearing event that is
/// not an SNV — an MNP (for reads with an indel or clip at the block), a delins,
/// and a deletion or insertion whose ALT also changes the anchor base (GC>T,
/// C>TA). SNVs and pure indels have their own classifiers. The dispatcher sends
/// exactly these here, and prep widens the reference for exactly these.
pub(crate) fn judges(v: &Variant) -> bool {
    // Every pure indel is judged by its bases, a deletion of 50+ bp by the junction
    // judge at its POS breakpoint (RJ-24).
    match super::window::allele_kind(&v.ref_allele, &v.alt_allele) {
        Some(AlleleKind::Mnp | AlleleKind::Complex | AlleleKind::Insertion) => true,
        Some(AlleleKind::Deletion) => v.ref_allele.len() - v.alt_allele.len() < super::junction::LONG_DELETION,
        _ => false,
    }
}

/// Whether a reference fetched from `start` is too short to judge the variant:
/// (on the left, on the right), for prep to fetch more. None when it holds the
/// event with [`ROOM`] bases each side, or the variant is not judged by this rule.
pub(crate) fn reference_short(start: i64, reference: &str, v: &Variant) -> Option<(bool, bool)> {
    if !judges(v) {
        return None;
    }
    let reference = upper(reference);
    let ev = event(start, &reference, v)?;
    let mut short = (ev.lo < ROOM, ev.hi + ROOM > reference.len());
    if !short.0 && !short.1 {
        // The windows run through bases the shorter allele's flank shares with the
        // longer allele: those must lie inside the fetched reference too.
        let (lo, hi) = (ev.lo - FLANK, ev.hi + FLANK);
        let alt_hi = (hi as i64 + ev.d) as usize;
        short = (
            shared_past(&reference, &ev, (lo, hi, alt_hi), true).is_none(),
            shared_past(&reference, &ev, (lo, hi, alt_hi), false).is_none(),
        );
    }
    if !short.0 && !short.1 {
        // Whatever the guard masks, the windows read on to a base telling the alleles
        // apart: past every base it can mask in them, the reference holds one, on the
        // side they read on to from either flank.
        let reach = (start + ev.lo.saturating_sub(GUARD_REACH) as i64, start + (ev.hi + GUARD_REACH) as i64);
        for from_right in [false, true] {
            let plain = Layout { from_right, ..Layout::default() };
            let Some(w) = windows_on(v, start, reference.clone(), &plain) else { continue };
            let every = Layout { from_right, masked: (reach.0..reach.1).collect(), ..Layout::default() };
            if windows_on(v, start, reference.clone(), &every).is_none() {
                if w.pair.0.tail {
                    short.0 = true;
                } else {
                    short.1 = true;
                }
            }
        }
    }
    (short.0 || short.1).then_some(short)
}

/// The variant's event grown through repeats touching it on either allele,
/// genomic [start, end): every base an aligner may place its change on. None
/// without a reference that holds it, or with an empty allele.
pub(crate) fn grown_event(v: &Variant) -> Option<(i64, i64)> {
    let (start, reference) = prepared_reference(v)?;
    let ev = event(start, &reference, v)?;
    Some((start + ev.lo as i64, start + ev.hi as i64))
}

/// Why the rule cannot judge a variant (`classify` gives None for every read),
/// for its once-per-variant warning.
pub(crate) fn unjudged_reason(v: &Variant) -> &'static str {
    match prepared_reference(v) {
        None => "no prepared reference (the variant was not prepared against a FASTA, or failed prep)",
        Some((start, reference)) if event(start, &reference, v).is_none() => {
            "the prepared reference does not hold its REF allele"
        }
        Some(_) if windows(v, &Layout::default()).is_some() => {
            "the recurring changes masked beside it leave no base telling the alleles apart \
             within the prepared reference"
        }
        Some(_) => {
            "the prepared reference does not hold the event with its flank (a contig end, or a \
             repeat past the fetch cap)"
        }
    }
}

/// The reference the rule reads, upper case, and its genomic start: the event's
/// widened reference when prep fetched one, else the reference context.
fn prepared_reference(v: &Variant) -> Option<(i64, Vec<u8>)> {
    prepared(v).map(|(s, seq)| (s, upper(seq)))
}

/// The reference a whole read is read over: prep's whole-read reach, else the prepared
/// reference.
fn whole_reference(v: &Variant) -> Option<(i64, Vec<u8>)> {
    v.whole_ref.as_ref().map(|(s, seq)| (*s, upper(seq))).or_else(|| prepared_reference(v))
}

/// The genome span `[start, end)` of the variant's prepared reference: as far as its
/// windows and their reading on reach.
fn prepared_span(v: &Variant) -> Option<(i64, i64)> {
    prepared(v).map(|(s, seq)| (s, s + seq.len() as i64))
}

fn prepared(v: &Variant) -> Option<(i64, &str)> {
    match (&v.event_ref, &v.ref_context) {
        (Some((s, seq)), _) => Some((*s, seq.as_str())),
        (None, Some(ctx)) => Some((v.ref_context_start, ctx.as_str())),
        (None, None) => None,
    }
}

fn upper(s: &str) -> Vec<u8> {
    s.bytes().map(|b| b.to_ascii_uppercase()).collect()
}

/// `reference` with `r` at offset `p` replaced by `a`. None when the reference does
/// not hold `r` there.
fn haplotype(reference: &[u8], p: usize, r: &[u8], a: &[u8]) -> Option<Vec<u8>> {
    if p + r.len() > reference.len() || reference[p..p + r.len()] != r[..] {
        return None;
    }
    let mut hap = Vec::with_capacity(reference.len() + a.len());
    hap.extend_from_slice(&reference[..p]);
    hap.extend_from_slice(a);
    hap.extend_from_slice(&reference[p + r.len()..]);
    Some(hap)
}

/// Another allele of the ALT's or REF's length at the variant's locus (RJ-25), on the
/// reference: a co-annotated sibling, or an allele the reads spell in full.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct OtherAllele {
    /// Genomic (0-based) position of its first base.
    pub pos: i64,
    /// Its REF and ALT bases, upper case.
    pub ref_allele: Vec<u8>,
    pub alt_allele: Vec<u8>,
}

impl OtherAllele {
    /// Its haplotype over a reference starting at genome position `start`. None when
    /// the reference does not hold its REF bases (a spliced reference whose junction
    /// crosses it, or one too short).
    fn on(&self, start: i64, reference: &[u8]) -> Option<Vec<u8>> {
        haplotype(reference, usize::try_from(self.pos - start).ok()?, &self.ref_allele, &self.alt_allele)
    }

    /// Its genome span `[pos, pos + REF length)`.
    fn span(&self) -> (i64, i64) {
        (self.pos, self.pos + self.ref_allele.len() as i64)
    }
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
    /// ...and the given allele's own bases placed once: for a pure insertion or
    /// deletion its leftmost placement (deleted bases, or none where bases are
    /// inserted), so the rest of the run it slides through is bases it keeps, where a
    /// change both alleles' molecules carry is masked at its place on each; for any
    /// other event the core.
    own: (usize, usize),
    /// ...and that grown through repeats touching it.
    lo: usize,
    hi: usize,
    /// A pure deletion: its own bases are REF's alone (the ALT has none there).
    deletes: bool,
    /// A pure insertion or deletion.
    pure: bool,
}

impl Event {
    /// The flank bases its windows hold each side.
    fn flank(&self) -> usize {
        if self.pure { PURE_FLANK } else { FLANK }
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
    let hap = haplotype(reference, p, &r_al, &a_al)?;
    let d = hap.len() as i64 - reference.len() as i64;

    let (l_lo, l_hi) = trimmed(reference, &hap, true);
    let (r_lo, r_hi) = trimmed(reference, &hap, false);
    // A pure insertion's or deletion's anchor is a base both alleles keep: its event is
    // where they differ over every placement (the trims' union, the tract). Any other
    // event holds the bases its alleles name.
    let pure = matches!(super::window::allele_kind(&v.ref_allele, &v.alt_allele), Some(AlleleKind::Insertion | AlleleKind::Deletion));
    let (b_lo, b_hi) =
        if pure { (l_lo.min(r_lo), l_hi.max(r_hi)) } else { (l_lo.min(r_lo).min(p), l_hi.max(r_hi).max(p + r_al.len())) };
    // Grow through repeats: a complex event's on either allele (another representation
    // of it may sit anywhere in one); a pure indel's on the allele holding its bases
    // (the ALT for an insertion, REF for a deletion), the repeat its change lengthens
    // or shortens, whose other lengths a read must span it to tell apart. A run on the
    // other allele only is beside the change. Left of the event the alleles share
    // coordinates; right of it the haplotype's are shifted by d.
    let unit = MAX_UNIT.max(d.unsigned_abs() as usize);
    let (on_ref, on_hap) = if pure { (d < 0, d > 0) } else { (true, true) };
    let hap_hi = (repeat_end(&hap, (b_hi as i64 + d) as usize, unit) as i64 - d) as usize;
    let lo = [on_ref.then(|| repeat_start(reference, b_lo, unit)), on_hap.then(|| repeat_start(&hap, b_lo, unit))]
        .into_iter()
        .flatten()
        .fold(b_lo, usize::min);
    let hi = [on_ref.then(|| repeat_end(reference, b_hi, unit)), on_hap.then_some(hap_hi)].into_iter().flatten().fold(b_hi, usize::max);
    let own = if pure { (r_lo, r_hi) } else { (b_lo, b_hi) };
    let deletes = pure && d < 0;
    Some(Event { hap, d, first: l_lo, last: r_hi, core: (b_lo, b_hi), own, lo, hi, deletes, pure })
}

/// How a variant's windows are laid out beyond the event itself: which flank they
/// are read from, the genome positions the neighbouring-change guard masks, and the
/// other alleles of the ALT's and of REF's length at the locus.
#[derive(Clone, Default)]
struct Layout {
    from_right: bool,
    /// The exon edge cutting the left flank when the windows read from the right: an
    /// unspliced read's windows end there, as a spliced read's do (the bases before it
    /// are the intron's, or on the mRNA the previous exon's).
    cut_left: Option<i64>,
    masked: Vec<i64>,
    /// The masked positions that take only some bases (background inside the event's
    /// span: the reference's base and the change's), with those bases.
    background: Vec<(i64, Bases)>,
    others: Vec<OtherAllele>,
    ref_others: Vec<OtherAllele>,
}

impl Layout {
    /// The layout for `v` under `rules`: its guard (when it was found for this very
    /// variant), and the flank it is read from (the left one, unless an RNA exon
    /// edge cuts it and not the right).
    fn for_variant(v: &Variant, rules: &ReadRules) -> Layout {
        let mut layout = rules
            .guard
            .and_then(|g| g.get(v))
            .map_or_else(Layout::default, |g| Layout {
                from_right: false,
                cut_left: None,
                masked: g.masked.clone(),
                background: g.background.clone(),
                others: g.others.clone(),
                ref_others: g.ref_others.clone(),
            });
        if let (Some(clip), Some((lo, hi))) = (rules.clip_edges, grown_event(v)) {
            let flank = FLANK as i64;
            if clip.cut_left(lo, flank) && !clip.cut_right(hi, flank) {
                trace!("carrier: an exon edge cuts the left flank of {}:{} → windows read from the right", v.chrom, v.pos + 1);
                layout.from_right = true;
                // an exon start there: the bases before it are intronic
                layout.cut_left = clip.edges().starts.iter().copied().filter(|&e| e > lo - flank && e <= lo).max();
            }
        }
        layout
    }
}

/// Whether one of the exon starts `starts` cuts the flank of `flank` bases before an
/// event whose first base is `lo`: the bases before it are intronic.
fn starts_cut_left(mut starts: impl Iterator<Item = i64>, lo: i64, flank: i64) -> bool {
    starts.any(|s| s > lo - flank && s <= lo)
}

/// REF and ALT windows for the variant, laid out by `layout`. None without a
/// reference that holds the event and its flank, or with an empty allele.
fn windows(v: &Variant, layout: &Layout) -> Option<Windows> {
    let (start, reference) = prepared_reference(v)?;
    windows_on(v, start, reference, layout)
}

/// [`windows`] on `reference`, fetched from `start`.
fn windows_on(v: &Variant, start: i64, reference: Vec<u8>, layout: &Layout) -> Option<Windows> {
    let ev = event(start, &reference, v)?;
    if ev.lo < FLANK || ev.hi + FLANK > reference.len() {
        return None; // the reference does not hold the event's flank
    }
    let d = ev.d;
    let flank = ev.flank();
    let lo = ev.lo - flank;
    let hi = ev.hi + flank;
    let alt_hi = (hi as i64 + d) as usize;
    let g = |off: usize| start + off as i64;
    // Windows of one length for both alleles: the shorter allele's whole window,
    // read from the POS (left) flank, or the right one when an RNA exon edge cuts
    // it; or, when both alleles are long, through one base past the first
    // difference from the flank nearer it (an event ending a long run is out of a
    // read's reach from the other).
    let short = (hi - lo).min(alt_hi - lo);
    let long = short > LONG_EVENT;
    let j_left = (ev.first - lo + 2).max(FLANK + 2).min(short);
    let j_right = (hi - ev.last + 2).max(FLANK + 2).min(short);
    let tail = layout.from_right || (long && j_right < j_left);
    // The shorter allele's flank may repeat the longer allele's next bases: both
    // windows run on through them to FLANK bases from where the alleles differ.
    let ext = if long { 0 } else { shared_past(&reference, &ev, (lo, hi, alt_hi), tail)? };
    let (lo, hi, alt_hi) = if tail { (lo - ext, hi, alt_hi) } else { (lo, hi + ext, alt_hi + ext) };
    let j = if !long { short + ext } else if tail { j_right } else { j_left };
    let (c_lo, c_hi) = ev.core;
    // Each other allele of the ALT's length and of REF's, on this reference (one it
    // does not hold, or of another length, is not read here).
    let other_haps: Vec<Vec<u8>> =
        layout.others.iter().filter_map(|o| o.on(start, &reference)).filter(|h| h.len() == ev.hap.len()).collect();
    let ref_other_haps: Vec<Vec<u8>> =
        layout.ref_others.iter().filter_map(|o| o.on(start, &reference)).filter(|h| h.len() == reference.len()).collect();
    // A masked base matches any read base where a change recurs beside the event, and
    // only the reference's base or the change's where the REF molecules carry one inside
    // the repeat it was grown through, or at a base an MNP keeps (a third base there is
    // another allele). Left of the given allele the alleles share coordinates; right of
    // it the ALT's are shifted by the length change. The given allele's own bases are
    // never masked (a base an MNP keeps is no change of its own).
    let (mut reference, mut hap) = (reference, ev.hap.clone());
    let mask = |b: u8, p: i64| IUPAC[(bit(b) | layout.background.iter().find(|&&(q, _)| q == p).map_or(ANY, |&(_, s)| s)) as usize];
    let (o_lo, o_hi) = ev.own;
    for &p in &layout.masked {
        let off = p - start;
        let (r, h) = if off >= lo as i64 && off < o_lo as i64 {
            (off as usize, off as usize)
        } else if off >= o_hi as i64 && off < hi as i64 {
            (off as usize, (off + d) as usize)
        } else if d == 0 && off >= c_lo as i64 && off < c_hi as i64 && reference[off as usize] == hap[off as usize] {
            (off as usize, off as usize)
        } else if ev.deletes && off >= o_lo as i64 && off < o_hi as i64 {
            // a base the deletion deletes: REF's alone
            reference[off as usize] = mask(reference[off as usize], p);
            continue;
        } else {
            continue;
        };
        reference[r] = mask(reference[r], p);
        hap[h] = mask(hap[h], p);
    }
    // Masks never take every base where the alleles differ: where they leave none in
    // the bases a read holds (a germline change that joins the run the event slides
    // through to the next, read to its end only), both windows read on to the flank bases
    // past the first base where the masked alleles take no base in common.
    let (lo, hi, alt_hi, j, masked_hold) = {
        let at = |t: usize| -> Option<(u8, u8)> {
            if tail {
                Some((*reference.get(hi.checked_sub(1 + t)?)?, *hap.get(alt_hi.checked_sub(1 + t)?)?))
            } else {
                Some((*reference.get(lo + t)?, *hap.get(lo + t)?))
            }
        };
        let tells = |t: usize| at(t).is_some_and(|(r, a)| takes(r) & takes(a) == 0);
        if (0..j).any(tells) {
            (lo, hi, alt_hi, j, None)
        } else {
            let first = (j..).take_while(|&t| at(t).is_some()).find(|&t| tells(t))?;
            let e = first + flank - j;
            let held = if tail { lo >= e } else { hi + e <= reference.len() && alt_hi + e <= hap.len() };
            if !held {
                return None;
            }
            trace!("carrier: masks leave no base telling the alleles apart in the windows → both read on {} bases", e);
            if tail { (lo - e, hi, alt_hi, j + e, Some(j)) } else { (lo, hi + e, alt_hi + e, j + e, Some(j)) }
        }
    };
    let (ref_win, alt_win) = (&reference[lo..hi], &hap[lo..alt_hi]);
    let n = ref_win.len().abs_diff(alt_win.len());
    let outer = if tail { reference[hi..(hi + n).min(reference.len())].to_vec() } else { reference[lo.saturating_sub(n)..lo].to_vec() };
    let kept = |haps: &[Vec<u8>], win: &[u8], (a, b): (usize, usize)| -> Vec<Vec<u8>> {
        // an allele takes the window's mask where its base is one the mask takes
        let masked = |o: u8, w: u8| if takes(w).count_ones() > 1 && fits(o, w) { w } else { o };
        haps.iter().map(|h| h[a..b].iter().zip(win).map(|(&o, &w)| masked(o, w)).collect()).collect()
    };
    let (others, ref_others) = (kept(&other_haps, alt_win, (lo, alt_hi)), kept(&ref_other_haps, ref_win, (lo, hi)));
    // Both windows are anchored at the window's ends: the left flank's first base
    // and one past the right flank's last, shared by the two alleles.
    // the event, from the first base the input names (a pure indel's anchor)
    let named = usize::try_from(v.pos - start).map_or(ev.lo, |p| p.max(lo));
    let window = |allele: &[u8], ev_end: usize| Window {
        allele: allele.to_vec(),
        hold: j,
        tail,
        left: Some(g(lo)),
        right: Some(g(hi)),
        before: c_lo - lo,
        after: hi - c_hi,
        past: if tail { Vec::new() } else { reference[hi..].to_vec() },
        outer: outer.clone(),
        event_at: (ev.lo.min(named) - lo, ev_end - lo),
    };
    let pair = (window(ref_win, ev.hi), window(alt_win, (ev.hi as i64 + d) as usize));
    Some(Windows { pair, event: (g(c_lo), g(c_hi)), own: (g(o_lo), g(o_hi)), others, ref_others, masked_hold })
}

/// Bases both windows must run on past the shorter allele's window (`lo`, `hi`,
/// `alt_hi`: its REF and ALT ends before the run), reading inward from the left
/// flank, or the right one when `tail`. Past the shorter allele's event its flank
/// may repeat the longer allele's next bases; the alleles are told apart only where
/// they differ, so the windows run through the shared bases to [`FLANK`] bases from
/// there, as a long deletion's run through its microhomology. Zero when the flank
/// differs at once, or for an event that keeps the length (its alleles share every
/// base past it). None when the windows run off the reference (through the shared
/// bases, or the flank past them).
fn shared_past(reference: &[u8], ev: &Event, (lo, hi, alt_hi): (usize, usize, usize), tail: bool) -> Option<usize> {
    if ev.d == 0 {
        return Some(0);
    }
    let short = (hi - lo).min(alt_hi - lo);
    // Index, from the reading flank, one past the shorter allele's event.
    let alt_ev_hi = (ev.hi as i64 + ev.d) as usize;
    let end = if tail { (hi - ev.lo).min(alt_hi - ev.lo) } else { (ev.hi - lo).min(alt_ev_hi - lo) };
    let bases = |t: usize| -> Option<(u8, u8)> {
        if tail {
            Some((*reference.get(hi.checked_sub(1 + t)?)?, *ev.hap.get(alt_hi.checked_sub(1 + t)?)?))
        } else {
            Some((*reference.get(lo + t)?, *ev.hap.get(lo + t)?))
        }
    };
    let mut first_diff = end;
    loop {
        let (r, a) = bases(first_diff)?;
        if r != a {
            break;
        }
        first_diff += 1;
    }
    let ext = (first_diff + ev.flank()).saturating_sub(short);
    let held = if tail { lo >= ext } else { hi + ext <= reference.len() && alt_hi + ext <= ev.hap.len() };
    held.then_some(ext)
}

/// The neighbouring-change guard of one variant: the genome positions of the
/// changes recurring beside the event, masked (a separate event, part of neither
/// allele), and the larger allele the ALT reads carry (for the diagnostic), with
/// the variant it was found for.
#[derive(Debug, Clone)]
pub(crate) struct Guard {
    key: (i64, String, String),
    pub masked: Vec<i64>,
    /// The masked positions that take only some bases (background inside the event's
    /// span), with those bases; every other masked position takes any.
    pub background: Vec<(i64, Bases)>,
    /// The reads' larger allele (0-based POS, left-aligned minimal VCF form): its
    /// carriers are the fewest ALT reads showing any one change, its given carriers
    /// the most showing the reference base at one.
    pub larger: Option<ObservedAllele>,
    /// Fragments read: those showing one allele's core exactly between the flank
    /// bases.
    pub scanned: u32,
    /// The other alleles of the ALT's length at the locus (RJ-25), and of REF's.
    pub others: Vec<OtherAllele>,
    pub ref_others: Vec<OtherAllele>,
    /// The indels the reads' own alignments write inside the windows' span that recur
    /// as background (on [`GUARD_MIN_READS`] fragments and [`GUARD_MIN_SHARE`] of
    /// those spanning the windows: a germline indel beside the event): its genome
    /// position and length (insertion +, deletion -), as [`own_indels`] gives them.
    pub indels: Vec<(i64, i64)>,
    /// A change at the ALT's shifted base the reads placed apart from the event show
    /// on too few of them, or in one pair orientation, to be a haplotype's, left
    /// unmasked: a minority haplotype carrying it would count ALT where its reads read
    /// as the ALT. Its genome position, the reference's base and the change's, and the
    /// carriers among the reads placed apart (the most carried, where several are).
    /// Diagnostic only (POSSIBLE_MINORITY_HAPLOTYPE).
    pub minority: Option<(i64, u8, u8, u32, u32)>,
}

impl Guard {
    /// The guard found for `v`.
    pub(crate) fn new(v: &Variant, masked: Vec<i64>, larger: Option<ObservedAllele>, scanned: u32) -> Self {
        Self {
            key: (v.pos, v.ref_allele.clone(), v.alt_allele.clone()),
            masked,
            background: Vec::new(),
            larger,
            scanned,
            others: Vec::new(),
            ref_others: Vec::new(),
            indels: Vec::new(),
            minority: None,
        }
    }
}

/// The insertions and deletions a read's own alignment writes inside `[lo, hi)`
/// (genome): an insertion between two bases of the span at the position after it, a
/// deletion overlapping it at its first deleted base, each with its length
/// (insertion +, deletion -).
fn own_indels(record: &Record, (lo, hi): (i64, i64)) -> Vec<(i64, i64)> {
    let (mut rp, mut out) = (record.pos(), Vec::new());
    for op in record.cigar().iter() {
        match op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::RefSkip(n) => rp += *n as i64,
            Cigar::Ins(n) => {
                if rp > lo && rp < hi {
                    out.push((rp, *n as i64));
                }
            }
            Cigar::Del(n) => {
                if rp < hi && rp + *n as i64 > lo {
                    out.push((rp, -(*n as i64)));
                }
                rp += *n as i64;
            }
            Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    out
}

/// A variant's [`Guard`], found from its reads on first use: only the exact-carrier
/// rule and the junction judge read it, so an SNV or a short pure indel never pays
/// for it.
pub(crate) struct CarrierGuard<'a> {
    found: std::cell::OnceCell<Option<Guard>>,
    find: Box<dyn Fn() -> Option<Guard> + 'a>,
    /// Per variant judged under these rules (the row's, and its co-annotated siblings
    /// judged for their claims), its alleles over the whole-read reach.
    whole: std::cell::RefCell<Vec<std::sync::Arc<Option<WholeHaps>>>>,
}

impl<'a> CarrierGuard<'a> {
    pub(crate) fn new(find: impl Fn() -> Option<Guard> + 'a) -> Self {
        Self { found: std::cell::OnceCell::new(), find: Box::new(find), whole: std::cell::RefCell::new(Vec::new()) }
    }

    /// The variant's alleles over the whole-read reach, built by `build` on its first
    /// use (they depend on the variant and its guard alone).
    fn whole(&self, v: &Variant, build: impl FnOnce() -> Option<WholeHaps>) -> std::sync::Arc<Option<WholeHaps>> {
        let of_v = |w: &WholeHaps| w.key.0 == v.pos && w.key.1 == v.ref_allele && w.key.2 == v.alt_allele;
        if let Some(w) = self.whole.borrow().iter().find(|w| w.as_ref().as_ref().is_some_and(of_v)) {
            return w.clone();
        }
        let w = std::sync::Arc::new(build());
        if w.is_some() {
            self.whole.borrow_mut().push(w.clone());
        }
        w
    }

    /// The guard when it was found for `v` (the dual-counted decomposed form of a
    /// variant shares its rules but not its guard).
    pub(crate) fn get(&self, v: &Variant) -> Option<&Guard> {
        self.found
            .get_or_init(|| (self.find)())
            .as_ref()
            .filter(|g| g.key == (v.pos, v.ref_allele.clone(), v.alt_allele.clone()))
    }
}

/// What the reads show placed by their bases from the event's left flank, as the
/// windows read them, not where the aligner wrote the gap (RJ-25): per fragment,
/// the spellings of another allele of REF's or the ALT's length, and the bases of the
/// REF and the ALT molecules over the event and [`FLANK`] bases each side.
struct Placed<'a> {
    /// Per allele's window (REF's `[lo, hi)` 0, the ALT's `[lo, alt_hi)` 1: the event
    /// with FLANK bases each side), the spellings that differ from that allele inside
    /// the event, every base read, by the fragments showing each.
    spellings: [HashMap<Vec<u8>, std::collections::HashSet<&'a [u8]>>; 2],
    /// Per reference offset in `[lo, hi)`: REF molecules reading each base (A C G T),
    /// and reading the base at all; the same for ALT molecules per haplotype offset in
    /// `[lo, alt_hi)`.
    ref_bases: Vec<[u32; 4]>,
    ref_cover: Vec<u32>,
    alt_bases: Vec<[u32; 4]>,
    alt_cover: Vec<u32>,
    /// Per allele's window, every spelling read in full, by the fragments showing
    /// each: a spelling must recur above the errors of the allele it is nearest, and a
    /// co-annotated sibling be spelled at all, counted from these once the masks are
    /// known.
    full: [HashMap<Vec<u8>, std::collections::HashSet<&'a [u8]>>; 2],
}

/// The [`Placed`] tallies of the variant's own reads in `read_cache` (those
/// overlapping its read window, `window`, as its counts read them: a bin's reads for
/// other rows never enter) at `ev`, over windows `[lo, hi)` (REF) and `[lo, alt_hi)`
/// (ALT), `lo` being FLANK bases before the event. A read is placed at its aligned
/// base on `lo`, and read through [`PLACED_REACH`] bases past the window. It is a
/// molecule of an allele when its confident bases (masked ones fit anything) differ
/// from that allele at no more than [`PLACED_CHANGES`] bases outside the event
/// (phased SNPs, an error), and, for an event changing the length, mismatch it at
/// fewer bases over the window and past it than the other allele placed the same way:
/// a molecule of the other length, shifted by the change, does not, so a carrier
/// written with its gap outside the event is no REF molecule, and where the bases past
/// the event repeat themselves (a run) the event's own bases tell the two apart.
/// Inside such an event a change fits (a germline SNP in the repeat it grew through),
/// and is tallied. For an event that keeps the length the alleles differ only inside
/// it, so each molecule must show its own allele's bases where the other has others.
fn placed_reads<'a>(
    read_cache: &'a [Record],
    counted: &dyn Fn(&Record) -> bool,
    min_baseq: u8,
    (start, reference): (i64, &[u8]),
    ev: &Event,
    (w_start, w_end): (i64, i64),
) -> Option<Placed<'a>> {
    let lo = ev.lo.checked_sub(FLANK)?;
    let (hi, alt_hi) = (ev.hi + FLANK, (ev.hi as i64 + ev.d) as usize + FLANK);
    if hi > reference.len() || alt_hi > ev.hap.len() {
        return None;
    }
    let mut placed = Placed {
        spellings: [HashMap::new(), HashMap::new()],
        ref_bases: vec![[0; 4]; hi - lo],
        ref_cover: vec![0; hi - lo],
        alt_bases: vec![[0; 4]; alt_hi - lo],
        alt_cover: vec![0; alt_hi - lo],
        full: [HashMap::new(), HashMap::new()],
    };
    let (mut ref_frags, mut alt_frags) = (std::collections::HashSet::new(), std::collections::HashSet::new());
    for record in read_cache {
        if !counted(record) || record.pos() >= w_end || crate::shared::bam_utils::ref_end(record) <= w_start {
            continue;
        }
        let (seq, quals) = (record.seq().as_bytes(), record.qual());
        let Some(q0) = find_read_pos(record, start + lo as i64) else { continue };
        if seq.len() > quals.len() {
            continue;
        }
        let read = |i: usize| seq.get(q0 + i).map(|&b| (b.to_ascii_uppercase(), quals[q0 + i] >= min_baseq));
        // Whether the read is a molecule of the allele whose haplotype is `allele`,
        // its window ending at `end` (the event's end plus FLANK): see above.
        let molecule = |allele: &[u8], end: usize, other: &[u8]| -> bool {
            let past = end + PLACED_REACH;
            if past > allele.len() {
                return false;
            }
            let (mut beside, mut misfits, mut other_misfits) = (0usize, 0usize, 0usize);
            for (i, &a) in allele.iter().enumerate().take(past).skip(lo) {
                match read(i - lo) {
                    Some((b, true)) if b != WILD => {
                        let miss = usize::from(b != a);
                        misfits += miss;
                        if i < ev.lo || i >= end - FLANK {
                            beside += miss;
                        }
                        other_misfits += usize::from(other.get(i).is_some_and(|&o| o != b));
                    }
                    Some(_) => {}
                    None => return false,
                }
            }
            let own = if ev.d != 0 {
                // A tie (a change joining the run to the next one makes both alleles
                // fit as far): the read's further bases, as far as both reach, decide.
                let mut i = past;
                while misfits == other_misfits && i < allele.len() && i < other.len() {
                    match read(i - lo) {
                        Some((b, true)) if b != WILD => {
                            misfits += usize::from(b != allele[i]);
                            other_misfits += usize::from(b != other[i]);
                        }
                        Some(_) => {}
                        None => break,
                    }
                    i += 1;
                }
                misfits < other_misfits
            } else {
                (ev.lo..end - FLANK)
                    .filter(|&i| allele[i] != other[i])
                    .all(|i| read(i - lo).is_some_and(|(b, ok)| ok && b == allele[i]))
            };
            beside <= PLACED_CHANGES && own
        };
        let qname = record.qname();
        if molecule(reference, hi, &ev.hap) && ref_frags.insert(qname) {
            tally(&mut placed.ref_bases, &mut placed.ref_cover, hi - lo, &read);
        }
        if molecule(&ev.hap, alt_hi, reference) && alt_frags.insert(qname) {
            tally(&mut placed.alt_bases, &mut placed.alt_cover, alt_hi - lo, &read);
        }
        // Another allele of either allele's length (REF's 0, the ALT's 1): every base of
        // that allele's window is read; for an event changing the length, past it the
        // read fits that allele placed this way at all but PLACED_CHANGES bases at most,
        // and no worse than the other allele placed the same way (an event keeping the
        // length places both alike); and the bases differ from that allele's inside the
        // event (beside it too, maybe: an allele's own trims may leave a base the other
        // allele changes in the flank) and are not the other allele's own.
        for (side, own, counter, end) in [(0usize, reference, &ev.hap[..], hi), (1, &ev.hap[..], reference, alt_hi)] {
            let n = end - lo;
            let bases: Option<Vec<u8>> = (0..n).map(|i| read(i).filter(|&(b, ok)| ok && b != WILD).map(|(b, _)| b)).collect();
            let Some(bases) = bases else { continue };
            placed.full[side].entry(bases.clone()).or_default().insert(qname);
            let past_fits = ev.d == 0 || {
                let past = end..end + PLACED_REACH;
                let misfits = |h: &[u8]| {
                    past.clone().filter(|&i| read(i - lo).is_some_and(|(b, ok)| ok && b != WILD && h.get(i) != Some(&b))).count()
                };
                past.end <= own.len()
                    && past.clone().all(|i| read(i - lo).is_some())
                    && misfits(own) <= PLACED_CHANGES
                    && misfits(own) <= misfits(counter)
            };
            let other = (FLANK..n - FLANK).any(|i| bases[i] != own[lo + i]) && counter.get(lo..end) != Some(&bases[..]);
            if past_fits && other {
                placed.spellings[side].entry(bases).or_default().insert(qname);
            }
        }
    }
    Some(placed)
}

/// The fragments placed apart from the event that read `off`, and those showing `base`
/// there: among the bases where the reference and the ALT shifted onto reads written
/// without its gap differ (`tells`, each with that ALT's base, all on one side of `off`),
/// those they hold aligned with `off` (no length change between) show the reference's
/// base more often than the shifted ALT's, so no placement of the event shifts their
/// base at `off` (an ALT read placed through the event shows its own at them). One read
/// per fragment; only reads overlapping the variant's guard window (`window`), which
/// every bin's fetch holds, so the tally is the same in any bin geometry. With them,
/// those of the F1R2 pair orientation (read 1 forward, read 2 reverse): of the
/// carriers, and of all.
fn read_apart(
    read_cache: &[Record],
    counted: &dyn Fn(&Record) -> bool,
    min_baseq: u8,
    (start, reference): (i64, &[u8]),
    (off, tells, base): (usize, &[(usize, u8)], u8),
    (w_start, w_end): (i64, i64),
) -> (u32, u32, u32, u32) {
    let mut fragments = std::collections::HashSet::new();
    let (mut n, mut c, mut f1r2, mut n_f1r2) = (0, 0, 0, 0);
    let own = |r: &&Record| counted(r) && r.pos() < w_end && crate::shared::bam_utils::ref_end(r) > w_start;
    for record in read_cache.iter().filter(own) {
        let Some((lo, hi, q_lo)) = ungapped_block(record, start + off as i64) else { continue };
        let (seq, quals) = (record.seq().as_bytes(), record.qual());
        let q = |o: usize| q_lo + (start + o as i64 - lo) as usize;
        let read = |o: usize| quals.get(q(o)).is_some_and(|&x| x >= min_baseq);
        if quals.len() < seq.len() || !read(off) {
            continue;
        }
        let held = tells.iter().filter(|&&(o, _)| start + o as i64 >= lo && (start + o as i64) < hi && read(o));
        let (on_ref, on_alt) = held.fold((0u32, 0u32), |(r, a), &(o, alt)| {
            let b = seq[q(o)].to_ascii_uppercase();
            (r + u32::from(b == reference[o]), a + u32::from(b == alt))
        });
        if on_ref <= on_alt || !fragments.insert(record.qname()) {
            continue;
        }
        n += 1;
        let (carries, oriented) = (seq[q(off)].to_ascii_uppercase() == base, record.is_reverse() == record.is_last_in_template());
        c += u32::from(carries);
        f1r2 += u32::from(carries && oriented);
        n_f1r2 += u32::from(oriented);
    }
    (n, c, f1r2, n_f1r2)
}

/// The genome span `[lo, hi)` a read aligns around genome position `at` with no length
/// change (no insertion, deletion, splice or clip inside it), and the query position of
/// its first base. None where the read does not align `at`.
fn ungapped_block(record: &Record, at: i64) -> Option<(i64, i64, usize)> {
    let (mut g, mut q, mut block) = (record.pos(), 0usize, None::<(i64, usize)>);
    let holds = |block: Option<(i64, usize)>, g: i64| block.filter(|&(s, _)| at >= s && at < g).map(|(s, sq)| (s, g, sq));
    for op in record.cigar().iter() {
        match *op {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                block = block.or(Some((g, q)));
                g += l as i64;
                q += l as usize;
                continue;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => continue,
            _ => {}
        }
        if let Some(found) = holds(block, g) {
            return Some(found);
        }
        block = None;
        match *op {
            Cigar::Ins(l) | Cigar::SoftClip(l) => q += l as usize,
            Cigar::Del(l) | Cigar::RefSkip(l) => g += l as i64,
            _ => {}
        }
    }
    holds(block, g)
}

/// What the fragments placed apart from the event tell of a change at a base where the
/// ALT's own base would show shifted ([`read_apart`]: `n` read the base, `c` show the
/// change; `f1r2` of the carriers and `n_f1r2` of all in the F1R2 pair orientation):
/// carried on a haplotype of their own (`Some(true)`) where [`GUARD_MIN_READS`] and
/// [`GUARD_APART_SHARE`] of them show it, in the pair orientations the reads are read in
/// (a haplotype's molecules are read in either; an oxoG or FFPE artifact sits in one: where
/// a fifth of all the reads placed apart lie in each orientation, carriers split one way
/// with a chance under [`GUARD_ORIENTATION_ALPHA`] against that balance are an artifact's;
/// an amplicon library, read in one orientation, weighs none); else the ALT's own base
/// (`Some(false)`) where at least [`GUARD_APART_READS`] read it; None where too few read it
/// to tell.
fn carried_apart(n: u32, c: u32, f1r2: u32, n_f1r2: u32) -> Option<bool> {
    let both = n_f1r2.min(n - n_f1r2) * 5 >= n.max(1);
    let p = n_f1r2 as f64 / n.max(1) as f64;
    let one_sided = both && binomial_cdf(f1r2, c, p).min(binomial_cdf(c - f1r2, c, 1.0 - p)) < GUARD_ORIENTATION_ALPHA;
    if c >= GUARD_MIN_READS && share(c, n) >= GUARD_APART_SHARE && !one_sided {
        Some(true)
    } else if n >= GUARD_APART_READS {
        Some(false)
    } else {
        None
    }
}

/// The chance of `x` or fewer successes in `n` trials at chance `p` each (summed in logs:
/// no term overflows at any depth).
fn binomial_cdf(x: u32, n: u32, p: f64) -> f64 {
    if x >= n || p <= 0.0 {
        return 1.0;
    }
    if p >= 1.0 {
        return 0.0;
    }
    let (lp, lq) = (p.ln(), (1.0 - p).ln());
    let (mut ln_c, mut sum) = (0.0f64, 0.0f64);
    for k in 0..=x {
        if k > 0 {
            ln_c += ((n - k + 1) as f64).ln() - (k as f64).ln();
        }
        sum += (ln_c + k as f64 * lp + (n - k) as f64 * lq).exp();
    }
    sum.min(1.0)
}

/// A share, 0 of nothing.
fn share(n: u32, of: u32) -> f64 {
    if of == 0 { 0.0 } else { n as f64 / of as f64 }
}

/// The reference offsets inside the repeat the event was grown through, off the
/// given allele's own bases, masked as beside the event (RJ-23, 2026-10-14), each with
/// the bases it takes: background takes the reference's base and the change's only (a
/// third base there is another allele); a phased neighbour any base, as beside the
/// event. A change
/// the REF molecules carry is background: recurring on them, or recurring on the ALT
/// molecules and on [`GUARD_MIN_READS`] and [`GUARD_REF_SHARE`] of them. A change the
/// ALT molecules alone carry there is none for an event changing the length (inside
/// the repeat it may be the event's own sequence); for an MNP, which has nothing to
/// place, the run is beside it, so an isolated one is its phased neighbour, masked as
/// beside any event (a run of them is another allele). For an MNP, the bases it keeps
/// inside its core too: they sit on both alleles as the grown part's do, so background
/// there is masked (a germline SNP in cis or in trans with it), while a change the ALT
/// molecules alone carry there is another allele. Read from the molecules placed by
/// their bases, so a carrier however aligned feeds it.
fn grown_background(placed: &Placed, reference: &[u8], ev: &Event) -> Vec<(usize, Bases)> {
    let lo = ev.lo - FLANK;
    // the given allele's own bases placed once (a pure indel's leftmost placement, the
    // rest of its run kept; any other event's core)
    let (c_lo, c_hi) = ev.own;
    let recurs = |n: u32, of: u32| n >= GUARD_MIN_READS && share(n, of) >= GUARD_MIN_SHARE;
    let grown = |off: usize| (ev.lo..c_lo).contains(&off) || (c_hi..ev.hi).contains(&off);
    // the bases a pure deletion deletes: REF's alone, so a change the REF molecules carry
    // there is background on REF (the ALT has no base there to carry it)
    let deleted = |off: usize| ev.deletes && (c_lo..c_hi).contains(&off);
    // the ALT's haplotype offset: left of the given allele the same, right of it
    // shifted by the length change
    let on_alt = |off: usize| if off < c_lo { off } else { (off as i64 + ev.d) as usize };
    // the background bases at an offset (none: no background there)
    let background = |off: usize| -> Bases {
        let h = on_alt(off);
        (0..4)
            .filter(|&k| {
                let (n, of) = (placed.ref_bases[off - lo][k], placed.ref_cover[off - lo]);
                let (an, aof) = (placed.alt_bases[h - lo][k], placed.alt_cover[h - lo]);
                b"ACGT"[k] != reference[off]
                    && (recurs(n, of) || (recurs(an, aof) && n >= GUARD_MIN_READS && share(n, of) >= GUARD_REF_SHARE))
            })
            .fold(0, |set, k| set | bit(b"ACGT"[k]))
    };
    let own = |off: usize| {
        ev.d == 0
            && grown(off)
            && background(off) == 0
            && (0..4).any(|k| b"ACGT"[k] != reference[off] && recurs(placed.alt_bases[off - lo][k], placed.alt_cover[off - lo]))
    };
    let isolated = |off: usize| !(own(off + 1) || off > 0 && own(off - 1));
    let kept = |off: usize| ev.d == 0 && reference[off] == ev.hap[off];
    // background takes the reference's base and its own; a phased neighbour any base
    let masked = |off: usize, b: Bases| (off, b | bit(reference[off]));
    (ev.lo..c_lo)
        .chain(c_hi..ev.hi)
        .filter_map(|off| match background(off) {
            0 => (own(off) && isolated(off)).then_some((off, ANY)),
            b => Some(masked(off, b)),
        })
        .chain((c_lo..c_hi).filter(|&off| kept(off)).filter_map(|off| match background(off) {
            0 => None,
            b => Some(masked(off, b)),
        }))
        .chain((c_lo..c_hi).filter(|&off| deleted(off)).filter_map(|off| {
            let b = (0..4)
                .filter(|&k| b"ACGT"[k] != reference[off] && recurs(placed.ref_bases[off - lo][k], placed.ref_cover[off - lo]))
                .fold(0, |set, k| set | bit(b"ACGT"[k]));
            (b != 0).then(|| masked(off, b))
        }))
        .collect()
}

/// Which allele a set of other alleles stands beside: REF (its window `[lo, hi)`) or
/// the ALT (`[lo, alt_hi)`); for an event keeping the length the two windows are one.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Side {
    Ref,
    Alt,
}

/// The other alleles of one allele's length at the locus (RJ-25), beside REF or the
/// ALT (`side`): the spellings read in full on [`GUARD_MIN_READS`] fragments or more
/// and on [`GUARD_REF_SHARE`] of the fragments spelling the allele they are nearest,
/// REF's own bases or the ALT's (the larger, when as near): an allele recurs above the
/// errors of the molecules it could be an error of, not as a recurring error does.
/// Save one whose every difference from the allele is a change from the other allele's
/// base that the other allele's molecules carry too (on GUARD_MIN_READS and
/// GUARD_REF_SHARE of those reading it, where the allele's base sits on the other's
/// haplotype: left-aligned, or right-aligned past the length change): a separate event
/// on both haplotypes. The other allele's own base is no such change: an allele lacking
/// part of the event, or carrying part of it, is another allele. And the co-annotated siblings that spell one (the allele's
/// haplotype save inside the event and [`FLANK`] bases each side, differing from it
/// inside the event) where the fragments spelling them in full are at least one and
/// GUARD_REF_SHARE of those spelling the allele they are nearest: a component of the
/// event the reads never show apart from it, but for an error, is part of it. None differs from the
/// allele only at `masked` positions (genome): that is the allele with a separate
/// event, which no read could tell from it. For an event keeping the length (an MNP)
/// none differs only beside the given allele's own bases: with nothing to place, that
/// is the allele with a phased neighbour (named, as beside any event).
fn other_alleles(
    placed: &Placed,
    (start, reference): (i64, &[u8]),
    ev: &Event,
    siblings: &[Variant],
    masks: &[(i64, Bases)],
    side: Side,
) -> Vec<OtherAllele> {
    let lo = ev.lo - FLANK;
    let (hi, alt_hi) = (ev.hi + FLANK, (ev.hi as i64 + ev.d) as usize + FLANK);
    let (c_lo, c_hi) = ev.own;
    // The masked positions on each allele's coordinates: REF's as the genome's; the
    // ALT's left of the given allele's own bases the same, right of them shifted by the
    // length change (for an MNP, the same throughout).
    let masked_ref: Vec<(usize, Bases)> =
        masks.iter().filter_map(|&(p, b)| usize::try_from(p - start).ok().map(|off| (off, b))).collect();
    let masked_hap: Vec<(usize, Bases)> = masked_ref
        .iter()
        .filter_map(|&(off, b)| match off {
            off if off < c_lo || ev.d == 0 => Some((off, b)),
            off if off >= c_hi => usize::try_from(off as i64 + ev.d).ok().map(|h| (h, b)),
            _ => None,
        })
        .collect();
    // Whether two bases at an offset are the same past the masks there (a background
    // base no difference, a third base still one).
    let same = |at: &[(usize, Bases)], i: usize, x: u8, y: u8| {
        x == y || at.iter().find(|m| m.0 == i).is_some_and(|&(_, b)| class(b, x) == class(b, y))
    };
    // this allele's haplotype, window end, masks and spellings; the other allele's
    let (own, counter): (&[u8], &[u8]) = if side == Side::Ref { (reference, &ev.hap) } else { (&ev.hap, reference) };
    let (masked_own, masked_counter) = if side == Side::Ref { (&masked_ref, &masked_hap) } else { (&masked_hap, &masked_ref) };
    let end = if side == Side::Ref { hi } else { alt_hi };
    let full = &placed.full[if side == Side::Ref { 0 } else { 1 }];
    let (own_spelling, counter_spelling) = (&own[lo..end], counter.get(lo..end));
    // How far a spelling is from an allele's, at the bases that can tell: not where the
    // guard masks (on that allele's coordinates), nor, for an MNP, beside its own bases
    // (a phased neighbour there is no difference).
    let apart = |b: &[u8], allele: &[u8], masked_at: &[(usize, Bases)]| {
        b.iter()
            .zip(allele)
            .enumerate()
            .filter(|&(k, (&x, &y))| {
                let i = lo + k;
                !same(masked_at, i, x, y) && (ev.d != 0 || (c_lo..c_hi).contains(&i))
            })
            .count()
    };
    let tellable = |h: &[u8]| {
        h.iter().zip(own).enumerate().any(|(i, (&x, &y))| !same(masked_own, i, x, y))
            && (ev.d != 0 || (c_lo..c_hi).any(|i| !same(masked_own, i, h[i], own[i])))
    };
    // the fragments spelling an allele's bases over the window, so measured
    let spellers = |allele: Option<&[u8]>, masked_at: &[(usize, Bases)]| -> u32 {
        let Some(allele) = allele else { return 0 };
        let mut frags = std::collections::HashSet::new();
        for (b, f) in full {
            if apart(b, allele, masked_at) == 0 {
                frags.extend(f.iter().copied());
            }
        }
        frags.len() as u32
    };
    let (own_n, counter_n) = (spellers(Some(own_spelling), masked_own), spellers(counter_spelling, masked_counter));
    // the fragments spelling the allele a spelling is nearest
    let source = |b: &[u8]| {
        let to_own = apart(b, own_spelling, masked_own);
        match counter_spelling.map(|r| apart(b, r, masked_counter)) {
            Some(to_counter) if to_counter < to_own => counter_n,
            Some(to_counter) if to_counter == to_own => own_n.max(counter_n),
            _ => own_n,
        }
    };
    // A change from the other allele's base that its molecules carry too, at a base the
    // event keeps: this allele's base sits on the other allele's haplotype (left-aligned,
    // or right-aligned past the length change) and is the same base there. Beside the
    // ALT, a non-reference change the REF molecules carry; beside REF, a change from the
    // ALT's base the ALT molecules carry. At a base the event changes the other allele's
    // molecules may be a third allele's, so a change there is no separate event.
    let (counter_bases, counter_cover) =
        if side == Side::Ref { (&placed.alt_bases, &placed.alt_cover) } else { (&placed.ref_bases, &placed.ref_cover) };
    let (counter_end, shift) = if side == Side::Ref { (ev.hi as i64 + ev.d, ev.d) } else { (ev.hi as i64, -ev.d) };
    let counter_carries = |i: usize, b: u8| {
        let left = ((i as i64) < counter_end).then_some(i);
        let right = usize::try_from(i as i64 + shift).ok().filter(|&j| i >= ev.lo && j >= ev.lo);
        let k = b"ACGT".iter().position(|&x| x == b);
        [left, right].into_iter().flatten().any(|j| {
            own[i] == counter[j]
                && k.is_some_and(|k| {
                    let (n, of) = (counter_bases[j - lo][k], counter_cover[j - lo]);
                    b != counter[j] && n >= GUARD_MIN_READS && share(n, of) >= GUARD_REF_SHARE
                })
        })
    };
    let separate =
        |b: &[u8]| b.iter().enumerate().all(|(k, &x)| same(masked_own, lo + k, x, own[lo + k]) || counter_carries(lo + k, x));
    let mut others: Vec<OtherAllele> = placed.spellings[if side == Side::Ref { 0 } else { 1 }]
        .iter()
        .filter(|(b, f)| {
            let n = f.len() as u32;
            counter_spelling.is_none_or(|r| apart(b, r, masked_counter) > 0)
                && n >= GUARD_MIN_READS
                && n as f64 >= GUARD_REF_SHARE * source(b) as f64
                && !separate(b)
        })
        // the spelling in place of REF's window bases (either allele's window starts at lo)
        .map(|(b, _)| OtherAllele { pos: start + lo as i64, ref_allele: reference[lo..hi].to_vec(), alt_allele: b.clone() })
        .collect();
    others.retain(|o| o.on(start, reference).is_some_and(|h| tellable(&h)));
    let mut haps: Vec<Vec<u8>> = others.iter().filter_map(|o| o.on(start, reference)).collect();
    for sib in siblings {
        let o = OtherAllele { pos: sib.pos, ref_allele: upper(&sib.ref_allele), alt_allele: upper(&sib.alt_allele) };
        let Some(h) = o.on(start, reference) else { continue };
        let spells = h.len() == own.len()
            && h[..lo] == own[..lo]
            && h[end..] == own[end..]
            && (ev.lo..end - FLANK).any(|i| h[i] != own[i]);
        // Spelled in full past the masks by at least one fragment and by GUARD_REF_SHARE
        // of those spelling the allele it is nearest: the annotation lowers the
        // GUARD_MIN_READS floor, not the share, so one error read at depth shows nothing.
        let seen = || {
            let n = spellers(Some(&h[lo..end]), masked_own);
            n > 0 && n as f64 >= GUARD_REF_SHARE * source(&h[lo..end]) as f64
        };
        if spells && tellable(&h) && !separate(&h[lo..end]) && !haps.contains(&h) && seen() {
            haps.push(h);
            others.push(o);
        }
    }
    others.sort_by(|a, b| (a.pos, &a.alt_allele).cmp(&(b.pos, &b.alt_allele)));
    others
}

/// Add a placed read's confident bases over a window's first `n` bases to `bases`
/// (A C G T) and `cover`.
fn tally(bases: &mut [[u32; 4]], cover: &mut [u32], n: usize, read: &dyn Fn(usize) -> Option<(u8, bool)>) {
    for i in 0..n {
        if let Some((b, true)) = read(i) {
            if let Some(k) = b"ACGT".iter().position(|&x| x == b) {
                cover[i] += 1;
                bases[i][k] += 1;
            }
        }
    }
}

/// The guard for a complex variant, from the reads `counted` admits: each read that
/// holds both of the event's flank bases aligned and shows one allele's core
/// exactly between them (masked bases fit) is read outward up to [`GUARD_REACH`]
/// bases on each side (as far as the reference holds), and its confident changes
/// from the reference are tallied per distance, once per fragment. A change on at
/// least [`GUARD_MIN_READS`] ALT fragments and [`GUARD_MIN_SHARE`] of those reading
/// that base, and on under [`GUARD_REF_SHARE`] of the REF fragments reading it, is
/// the ALT reads' alone: with the given allele, these make the larger allele named
/// in the diagnostic, unless the ALT fragments' own changes include a run (the
/// tallies cannot spell another allele: an MNP, or an indel shifting the bases past
/// it). A change beside the event is a separate event, so a read carrying the given
/// allele counts whatever it shows there: masked are the changes the REF fragments
/// carry too (a germline SNP or MNP, background on both haplotypes, runs included:
/// recurring on the REF fragments, or on the ALT fragments and on at least
/// [`GUARD_MIN_READS`] and [`GUARD_REF_SHARE`] of the REF fragments) and each
/// isolated change of the ALT fragments' own. A run of their own changes (a
/// neighbouring distance recurring there, the bases past the window included) is an
/// indel or MNP beside the event, another allele: not masked. The variant's own reads,
/// placed by their bases ([`placed_reads`]), add the background inside the repeat the
/// event was grown through and at the bases an MNP keeps ([`grown_background`]; the
/// scan lets it fit) and the other
/// alleles of the ALT's and of REF's length at the locus ([`other_alleles`], RJ-25).
/// None when nothing is masked, named or found, or the variant is not judged by this
/// rule.
pub(crate) fn guard_for(
    read_cache: &[Record],
    v: &Variant,
    counted: &dyn Fn(&Record) -> bool,
    min_baseq: u8,
    siblings: &[Variant],
) -> Option<Guard> {
    if !judges(v) {
        return None;
    }
    let pure = matches!(super::window::allele_kind(&v.ref_allele, &v.alt_allele), Some(AlleleKind::Insertion | AlleleKind::Deletion));
    let (start, reference) = prepared_reference(v)?;
    let ev = event(start, &reference, v)?;
    // The scan beside the event starts past the repeat it was grown through (for an
    // MNP too, so a read written with its block shifted is still read, and its reach
    // still meets the windows' flank); inside the repeat, the molecules placed by their
    // bases tell what both alleles carry and, for an MNP, its phased neighbours.
    if ev.lo == 0 || ev.hi >= reference.len() {
        return None;
    }
    let reach_l = GUARD_REACH.min(ev.lo);
    let reach_r = GUARD_REACH.min(reference.len() - ev.hi);
    let alt_core = &ev.hap[ev.lo..(ev.hi as i64 + ev.d) as usize];
    let ref_core = &reference[ev.lo..ev.hi];
    let g = |off: usize| start + off as i64;
    // The variant's own reads placed by their bases (RJ-25): the background inside the
    // repeat the event was grown through, which the scan beside the event lets fit
    // too, so a read carrying it is still read on past the event.
    let placed = placed_reads(read_cache, counted, min_baseq, (start, &reference), &ev, super::window::read_window(v));
    let grown = placed.as_ref().map_or_else(Vec::new, |p| grown_background(p, &reference, &ev));
    let c_lo = ev.own.0;
    let on_alt = |off: usize| if off < c_lo { off } else { (off as i64 + ev.d) as usize };
    // the bases each core position takes in the scan's fit (0: its own only)
    let ref_wild: Vec<Bases> = (ev.lo..ev.hi).map(|off| grown.iter().find(|m| m.0 == off).map_or(0, |m| m.1)).collect();
    let alt_wild: Vec<Bases> =
        (ev.lo..ev.lo + alt_core.len()).map(|h| grown.iter().find(|m| on_alt(m.0) == h).map_or(0, |m| m.1)).collect();
    // Per allele (REF 0, ALT 1) and side (left 0, right 1), by distance from the
    // event (index 0 is the flank base beside it): reads reading the base, and reads
    // showing a confident change there.
    let mut cover = vec![vec![vec![0u32; GUARD_REACH]; 2]; 2];
    let mut scanned = 0u32;
    let mut change = cover.clone();
    // the bases the ALT reads show where they change, per side and distance (A C G T)
    let mut alt_bases = vec![vec![[0u32; 4]; GUARD_REACH]; 2];
    // and the bases the REF fragments show where they change
    let mut ref_bases = alt_bases.clone();
    // one tally per fragment: overlapping mates are one molecule's evidence
    let mut fragments: std::collections::HashSet<&[u8]> = std::collections::HashSet::new();
    for record in read_cache {
        if !counted(record) {
            continue;
        }
        let seq = record.seq().as_bytes();
        let quals = record.qual();
        if seq.is_empty() || quals.len() < seq.len() {
            continue;
        }
        let (Some(ql), Some(qr)) = (find_read_pos(record, g(ev.lo - 1)), find_read_pos(record, g(ev.hi))) else {
            continue;
        };
        if qr <= ql {
            continue;
        }
        // a read writing an indel beside the event shifts its bases there: no molecule
        // whose changes the tally can place
        let (e_lo, e_hi) = (g(ev.lo), g(ev.hi));
        let beside = own_indels(record, (e_lo - reach_l as i64, e_hi + reach_r as i64))
            .into_iter()
            .any(|(p, n)| if n > 0 { p < e_lo || p > e_hi } else { p - n <= e_lo || p >= e_hi });
        if beside {
            continue;
        }
        let between = &seq[ql + 1..qr];
        let fits = |core: &[u8], wild: &[Bases]| {
            core.len() == between.len()
                && between.iter().zip(core).enumerate().all(|(i, (&b, &c))| {
                    let b = b.to_ascii_uppercase();
                    b == c || b == WILD || quals[ql + 1 + i] < min_baseq || class(wild[i], b) == class(wild[i], c)
                })
        };
        let allele = match (fits(ref_core, &ref_wild), fits(alt_core, &alt_wild)) {
            (true, false) => 0,
            (false, true) => 1,
            _ => continue,
        };
        if !fragments.insert(record.qname()) {
            continue;
        }
        scanned += 1;
        for (side, reach) in [(0usize, reach_l), (1usize, reach_r)] {
            for i in 0..reach {
                let off = if side == 0 { ev.lo - 1 - i } else { ev.hi + i };
                let Some(q) = find_read_pos(record, g(off)) else { break };
                cover[allele][side][i] += 1;
                let b = seq[q].to_ascii_uppercase();
                if b != WILD && quals[q] >= min_baseq && b != reference[off] {
                    change[allele][side][i] += 1;
                    if let Some(k) = b"ACGT".iter().position(|&x| x == b) {
                        let tally = if allele == 1 { &mut alt_bases } else { &mut ref_bases };
                        tally[side][i][k] += 1;
                    }
                }
            }
        }
    }
    // Whether a change recurs on an allele's fragments (REF 0, ALT 1) at a side and
    // distance; past the reach, nothing is known to recur.
    let recurs = |allele: usize, side: usize, i: usize| {
        let (n, of) = (change[allele][side].get(i).copied().unwrap_or(0), cover[allele][side].get(i).copied().unwrap_or(0));
        n >= GUARD_MIN_READS && share(n, of) >= GUARD_MIN_SHARE
    };
    let ref_changes = |side: usize, i: usize| change[0][side].get(i).copied().unwrap_or(0);
    let on_ref =
        |side: usize, i: usize| share(ref_changes(side, i), cover[0][side].get(i).copied().unwrap_or(0)) >= GUARD_REF_SHARE;
    // Background, a change the REF fragments carry too (a germline SNP or MNP):
    // recurring there, or on the ALT fragments and on as many REF fragments as a
    // recurring change needs (one REF read's error is no background).
    let background = |side: usize, i: usize| {
        recurs(0, side, i) || (recurs(1, side, i) && ref_changes(side, i) >= GUARD_MIN_READS && on_ref(side, i))
    };
    // The ALT fragments' own changes; isolated when no own change sits at either
    // neighbouring distance (the event itself is no neighbour).
    let own = |side: usize, i: usize| recurs(1, side, i) && !background(side, i);
    let isolated = |side: usize, i: usize| !(own(side, i + 1) || i > 0 && own(side, i - 1));
    // per side, the distances of the changes the ALT reads alone carry; and whether
    // their own changes include a run
    // The bases the input names (a pure indel's anchor beside its tract) are the given
    // allele's own: never masked, as the event's never are.
    let named = usize::try_from(v.pos - start).map_or(0..0, |p| p..p + v.ref_allele.len());
    // A change giving the reads the ALT's own base there, shifted onto them where the
    // aligner wrote no length change (the base an ALT read placed from the event's far
    // flank shows; a germline change may give the same base), is told by the fragments
    // placed apart from the event ([`read_apart`]: by the bases past it, away from the
    // event, where the reference and that shifted ALT differ). Carried by enough of them
    // ([`carried_apart`]), a haplotype holds it, read as the ALT over a short read:
    // masked as background, whatever its share of the REF fragments. Otherwise, read there
    // by enough of them, a background change is the ALT's own base on reads placed
    // through the event: not masked, and spelling no other allele. With too few of them,
    // it is masked where it is background.
    // (read over the whole-read reach: the telling bases may lie past the prepared
    // reference, a repeat's unit or more away)
    let wide = whole_reference(v).and_then(|(ws, wr)| event(ws, &wr, v).map(|we| (ws, wr, we)));
    // (the verdict, the reads placed apart, the carriers, the reference's base and the
    // change's)
    let apart = |side: usize, i: usize, off: usize| -> Option<(Option<bool>, u32, u32, u8, u8)> {
        let (ws, wr, we) = wide.as_ref()?;
        let off = usize::try_from(g(off) - ws).ok()?;
        let alt_at = |o: usize| {
            if side == 1 { we.hap.get(o).copied() } else { usize::try_from(o as i64 + we.d).ok().and_then(|h| we.hap.get(h).copied()) }
        };
        let base = alt_at(off).filter(|&b| wr.get(off).is_some_and(|&r| b != r))?;
        let k = b"ACGT".iter().position(|&x| x == base)?;
        // the base most REF fragments changing there show, background or on more of
        // them than errors would (few REF fragments may read it: they span the event)
        let (shown, of) = (ref_bases[side][i][k], cover[0][side][i]);
        let most = (0..4).max_by_key(|&k| ref_bases[side][i][k]) == Some(k);
        if !most || !(background(side, i) || (shown > 0 && share(shown, of) > GUARD_APART_ERRORS)) {
            return None;
        }
        // the bases past it, away from the event, where the reference and that shifted
        // ALT differ, with the ALT's base there
        let tell = |o: usize| alt_at(o).filter(|&b| b != wr[o]).map(|b| (o, b));
        let tells: Vec<(usize, u8)> =
            if side == 1 { (off + 1..wr.len()).filter_map(tell).collect() } else { (0..off).filter_map(tell).collect() };
        if tells.is_empty() {
            return None;
        }
        let (n, c, f1r2, n_f1r2) = read_apart(read_cache, counted, min_baseq, (*ws, wr), (off, &tells, base), super::window::guard_window(v));
        Some((carried_apart(n, c, f1r2, n_f1r2), n, c, wr[off], base))
    };
    let mut minority: Option<(i64, u8, u8, u32, u32)> = None;
    let mut shifted = Vec::new();
    let (mut masked, mut alt_only, mut run) = (Vec::new(), [Vec::new(), Vec::new()], false);
    for (side, reach) in [(0usize, reach_l), (1usize, reach_r)] {
        for i in 0..reach {
            let off = if side == 0 { ev.lo - 1 - i } else { ev.hi + i };
            if named.contains(&off) {
                continue;
            }
            if recurs(1, side, i) && !on_ref(side, i) {
                alt_only[side].push(i);
            }
            let found = apart(side, i, off);
            let told = found.and_then(|f| f.0);
            if background(side, i) && told == Some(false) {
                shifted.push(g(off));
            } else if background(side, i) || told == Some(true) || (own(side, i) && isolated(side, i)) {
                masked.push(g(off));
            }
            // carriers too few (or one-sided) for a haplotype, the change left unmasked:
            // from two of them, above the errors' share (one is an error's)
            let named = |f: &(Option<bool>, u32, u32, u8, u8)| f.2 >= 2 && share(f.2, f.1) > GUARD_APART_ERRORS;
            if let Some((_, n, c, r, b)) = found.filter(|f| named(f) && f.0 != Some(true) && !masked.contains(&g(off))) {
                if minority.is_none_or(|m| c > m.3) {
                    minority = Some((g(off), r, b, c, n));
                }
            }
            run |= own(side, i) && !isolated(side, i);
        }
    }
    // The larger allele: the given one with every change the ALT reads alone carry,
    // each base as most of them show it; its carriers are the fewest showing any one.
    // Beside a run, none: the tallies cannot spell another allele (an MNP, or an
    // indel shifting the bases past it); the observed allele names it from the
    // reads' own alignments.
    let larger = if run || alt_only.iter().all(Vec::is_empty) {
        None
    } else {
        let far = |side: usize| alt_only[side].iter().max().copied();
        let lo = far(0).map_or(ev.lo, |i| ev.lo - 1 - i);
        let hi = far(1).map_or(ev.hi, |i| ev.hi + i + 1);
        let (mut carriers, mut given) = (u32::MAX, 0u32);
        for side in 0..2 {
            for &i in &alt_only[side] {
                carriers = carriers.min(change[1][side][i]);
                given = given.max(cover[1][side][i] - change[1][side][i]);
            }
        }
        let base_at = |side: usize, i: usize, off: usize| {
            if !alt_only[side].contains(&i) {
                return reference[off];
            }
            let k = (0..4).max_by_key(|&k| alt_bases[side][i][k]).unwrap_or(0);
            b"ACGT"[k]
        };
        let mut alt = Vec::with_capacity(hi - lo + alt_core.len());
        alt.extend((lo..ev.lo).map(|off| base_at(0, ev.lo - 1 - off, off)));
        alt.extend_from_slice(alt_core);
        alt.extend((ev.hi..hi).map(|off| base_at(1, off - ev.hi, off)));
        // named like every observed allele: left-aligned, minimal VCF form
        let base = |p: i64| usize::try_from(p - start).ok().and_then(|o| reference.get(o).copied());
        let text = |x: Vec<u8>| String::from_utf8_lossy(&x).into_owned();
        canonical(g(lo), reference[lo..hi].to_vec(), alt, &base).map(|(pos, r, a)| ObservedAllele {
            pos,
            ref_allele: text(r),
            alt_allele: text(a),
            carriers,
            given_carriers: given,
        })
    };
    // The background inside the grown event joins the masks; the other alleles of each
    // allele's length are those a read can tell from it past every mask.
    let mut background = Vec::new();
    for &(off, b) in &grown {
        let p = g(off);
        if !masked.contains(&p) {
            masked.push(p);
            if b != ANY {
                background.push((p, b));
            }
        }
    }
    let masks: Vec<(i64, Bases)> =
        masked.iter().map(|&p| (p, background.iter().find(|&&(q, _)| q == p).map_or(ANY, |&(_, b)| b))).collect();
    // The ALT's own bases shifted onto reads placed through the event spell no other
    // allele of REF's length: those reads are the ALT's.
    let apart: Vec<(i64, Bases)> = masks.iter().copied().chain(shifted.iter().map(|&p| (p, ANY))).collect();
    let others_of = |side| {
        placed.as_ref().map_or_else(Vec::new, |p| other_alleles(p, (start, &reference), &ev, siblings, &apart, side))
    };
    let (others, ref_others) = (others_of(Side::Alt), others_of(Side::Ref));
    // The indels the reads' own alignments write inside the windows' span, the event's
    // own aside, per fragment, against the fragments spanning the windows: those
    // recurring are background (a germline indel beside the event), every other one
    // another allele.
    let mut indels = Vec::new();
    if let Some(span) = windows(v, &Layout::default()).map(|w| w.span()) {
        let mut carrying: HashMap<(i64, i64), std::collections::HashSet<&[u8]>> = HashMap::new();
        let mut spanning: std::collections::HashSet<&[u8]> = std::collections::HashSet::new();
        for record in read_cache.iter().filter(|r| counted(r)) {
            if record.pos() <= span.0 && crate::shared::bam_utils::ref_end(record) >= span.1 {
                spanning.insert(record.qname());
            }
            // a read whose indels there change the length by the event's own change
            // writes the event itself, however the aligner placed it: no indel beside it
            let ops = own_indels(record, span);
            if ops.iter().map(|o| o.1).sum::<i64>() == ev.d {
                continue;
            }
            // Inside the span a length change slides through another indel is another
            // length of that run, an allele at the locus, never background: for a pure
            // indel its discrimination window (the tract with a base each side), for any
            // other event the run it was grown through.
            let (e_lo, e_hi) = if pure { super::window::discrimination_window(v) } else { (g(ev.lo), g(ev.hi)) };
            let beside = |&&(p, n): &&(i64, i64)| {
                ev.d == 0 || if n > 0 { p <= e_lo || p >= e_hi } else { p - n <= e_lo || p >= e_hi }
            };
            for &o in ops.iter().filter(beside) {
                carrying.entry(o).or_default().insert(record.qname());
            }
        }
        let of = spanning.len() as u32;
        indels = carrying
            .into_iter()
            .filter(|(_, f)| f.len() as u32 >= GUARD_MIN_READS && share(f.len() as u32, of) >= GUARD_MIN_SHARE)
            .map(|(o, _)| o)
            .collect();
        indels.sort_unstable();
    }
    if masked.is_empty() && larger.is_none() && others.is_empty() && ref_others.is_empty() && indels.is_empty() && minority.is_none() {
        return None;
    }
    trace!(
        "carrier guard at {}:{} {}>{}: {} recurring changes beside the event masked ({} inside the repeat it grew \
         through), {} the ALT's own bases shifted onto reads placed through it (not masked), {} on the ALT reads \
         alone, {} other alleles of its length, {} of REF's, {} recurring indels inside the windows",
        v.chrom, v.pos + 1, v.ref_allele, v.alt_allele, masked.len(), grown.len(), shifted.len(),
        alt_only.iter().map(Vec::len).sum::<usize>(), others.len(), ref_others.len(), indels.len(),
    );
    let mut guard = Guard::new(v, masked, larger, scanned);
    guard.background = background;
    guard.others = others;
    guard.ref_others = ref_others;
    guard.indels = indels;
    guard.minority = minority;
    Some(guard)
}

/// Whether a read's own indels (`ops`) change the length with one the guard does not
/// hold as background (a germline indel beside the event). An insertion and a deletion
/// that cancel (a block written shifted, an insertion before it and a deletion after)
/// change no length.
fn foreign(ops: &[(i64, i64)], guard: Option<&Guard>) -> bool {
    ops.iter().map(|o| o.1).sum::<i64>() != 0 && ops.iter().any(|o| !guard.is_some_and(|g| g.indels.contains(o)))
}

/// How far beside the windows a read's own indel is placed by its bases: an insert
/// written outside them that can sit inside the window is never a separate event
/// (RJ-22, up to 64 junctions each way).
const PLACE_REACH: i64 = 64;

/// Whether a read's own alignment writes a deletion directly beside one of its splices
/// (an `N` op next to a `D` op) inside `span` (genome).
fn splice_beside_deletion(record: &Record, (lo, hi): (i64, i64)) -> bool {
    let (mut rp, mut prev): (i64, Option<char>) = (record.pos(), None);
    for op in record.cigar().iter() {
        let (kind, len) = match op {
            Cigar::RefSkip(n) => ('N', *n as i64),
            Cigar::Del(n) => ('D', *n as i64),
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => ('M', *n as i64),
            Cigar::Ins(_) => ('I', 0),
            _ => continue,
        };
        let pair = matches!((prev, kind), (Some('N'), 'D') | (Some('D'), 'N'));
        if pair && rp < hi && rp + len.max(1) > lo {
            return true;
        }
        rp += len;
        prev = Some(kind);
    }
    false
}

/// The query position where a read's own first indel at or past genome position
/// `from` starts (the first base its alignment shifts there), if any.
fn next_indel_query(record: &Record, from: i64) -> Option<usize> {
    let (mut rp, mut qp) = (record.pos(), 0usize);
    for op in record.cigar().iter() {
        match op {
            Cigar::Ins(n) => {
                if rp >= from {
                    return Some(qp);
                }
                qp += *n as usize;
            }
            Cigar::Del(n) => {
                if rp >= from {
                    return Some(qp);
                }
                rp += *n as i64;
            }
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                rp += *n as i64;
                qp += *n as usize;
            }
            Cigar::RefSkip(n) => rp += *n as i64,
            Cigar::SoftClip(n) => qp += *n as usize,
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    None
}

/// What a read's own alignment writes inside the windows besides the given allele, by
/// the read's bases (RJ-7, RJ-8).
enum OwnEvent {
    /// Nothing the bases show: no length change, bases that show no indel (a read
    /// ending in a run), or a change of the ALT's length its bases place nowhere.
    None,
    /// A change the windows judge as written: of the ALT's length, its bases unread
    /// from the windows' first base (spliced, or under its own deletion), or background
    /// (a germline indel beside the event).
    AsWritten,
    /// The ALT's own length inside the window, the read's bases spelling the ALT (one
    /// error aside): the windows judge it from the query position its bases put their
    /// first base at (`left`), however the aligner wrote the gap.
    Given { left: i64 },
    /// Another event inside the discrimination window (RJ-7), or one no placement fits:
    /// another allele.
    Inside,
    /// Another event outside the discrimination window (RJ-8): `ins` inserted before
    /// genome base `at`, or `del` reference bases deleted from it.
    Separate { at: i64, ins: Vec<u8>, del: usize },
}

/// A read's own event inside the windows ([`OwnEvent`]): its own alignment changes the
/// length inside them, and its bases, from the windows' first base to its end, its next
/// splice or its next indel, place that change where they fit the reference with it best (the fewest
/// mismatches, one at most), whatever the aligner wrote. Inside the window (RJ-7, for
/// a pure indel its tract with a base each side, for any other event the windows) when
/// one best placement is, outside it (RJ-8) when none is. Of the given allele's own
/// length, or an indel the guard holds as background (a germline indel beside the
/// event), inside the window or placed nowhere, the windows judge it (the ALT, or
/// another allele of its length); any other change placed nowhere is another allele. A read whose bases fit the reference as it is shows no
/// indel (one ending inside a run).
fn own_event(record: &Record, seq: &[u8], quals: &[u8], min_baseq: u8, v: &Variant, win: &Windows, guard: Option<&Guard>) -> OwnEvent {
    let span = win.span();
    // The read's indels inside the windows; with none there, the one nearest them
    // within reach, placed by its bases alone (others further off lie past it: each a
    // separate event, never summed into one change).
    let ops = match own_indels(record, span) {
        inside if !inside.is_empty() => inside,
        _ => own_indels(record, (span.0 - PLACE_REACH, span.1 + PLACE_REACH))
            .into_iter()
            .min_by_key(|&(p, n)| if p < span.0 { span.0 - p + n.min(0) } else { p - span.1 })
            .into_iter()
            .collect(),
    };
    let n: i64 = ops.iter().map(|o| o.1).sum();
    if n == 0 {
        return OwnEvent::None;
    }
    // an indel recurring as background (a germline indel beside the event) is judged
    // by the windows inside the window, a separate event outside it
    let background = !foreign(&ops, guard);
    let Some((start, reference)) = prepared_reference(v) else { return OwnEvent::Inside };
    // read from the base before the read's first indel there, or the windows' first
    let first = ops.iter().map(|o| o.0 - 1).min().unwrap_or(span.0).min(span.0).max(start);
    // A read whose bases cannot be read from there (spliced, or deleted) places
    // nothing: the given allele's own length is the windows' to judge.
    let own_length = n == v.alt_allele.len() as i64 - v.ref_allele.len() as i64;
    let unplaced = if own_length || background { OwnEvent::AsWritten } else { OwnEvent::Inside };
    let Some((x0, q0)) = [first, span.0].into_iter().find_map(|x| find_read_pos(record, x).map(|q| (x, q))) else {
        return unplaced;
    };
    let Ok(r0) = usize::try_from(x0 - start) else { return unplaced };
    // the read's bases from the windows' first base to its end, its next splice, or its
    // next indel past the ones placed (a separate event shifts the bases beyond it, as
    // reading on past the windows stops there)
    let (mut qp, mut stop) = (0usize, seq.len());
    for op in record.cigar().iter() {
        match op {
            Cigar::RefSkip(_) if qp > q0 => {
                stop = qp;
                break;
            }
            Cigar::Match(k) | Cigar::Equal(k) | Cigar::Diff(k) | Cigar::Ins(k) | Cigar::SoftClip(k) => qp += *k as usize,
            _ => {}
        }
    }
    let placed_end = ops.iter().map(|&(p, n)| p + (-n).max(0)).max().unwrap_or(span.1).max(span.1);
    if let Some(q) = next_indel_query(record, placed_end + 1) {
        stop = stop.min(q);
    }
    let bases: Vec<(u8, bool)> = (q0..stop.min(seq.len()).min(quals.len()))
        .map(|q| (seq[q].to_ascii_uppercase(), quals[q] >= min_baseq && seq[q].to_ascii_uppercase() != WILD))
        .collect();
    let refr = &reference[r0.min(reference.len())..];
    let m = n.unsigned_abs() as usize;
    // Every placement is scored over the same read bases: as far as the reference
    // reaches with the length change anywhere (it ends a few bases past the windows).
    let reach = bases.len().min(if n > 0 { refr.len() + m } else { refr.len().saturating_sub(m) });
    // the read's confident mismatches against REF with the length change at read
    // offset p (its inserted bases the read's own), counted up to `cap`
    let misses_at = |p: usize, cap: usize| {
        let mut errors = 0;
        for (k, &(b, ok)) in bases[..reach].iter().enumerate() {
            let h = match (n > 0, k < p, k < p + m) {
                (_, true, _) => refr.get(k),
                (true, false, true) => Some(&b),
                (true, false, false) => refr.get(k - m),
                (false, false, _) => refr.get(k + m),
            };
            // past the reference held here nothing is compared (no placement reaches it)
            let Some(&h) = h else { return cap + 1 };
            if ok && !fits(b, h) {
                errors += 1;
                if errors > cap {
                    break;
                }
            }
        }
        errors
    };
    // an insertion sits before a reference base it holds (or at its end), a deletion
    // removes bases it holds
    let places = reach.min(if n > 0 { refr.len() + 1 } else { refr.len().saturating_sub(m) + 1 });
    let errors: Vec<usize> = (0..places).map(|p| misses_at(p, 1)).collect();
    // Bases decide: a read whose bases fit the reference as it is (one error) shows no
    // indel, however its alignment writes one (a read ending inside a run, or one base
    // past it: RJ-1). Where the read runs on past the reference held here, one mismatch
    // the change resolves shows it: the read's next bases would show more.
    let plain = bases[..reach.min(refr.len())].iter().zip(refr).filter(|&(&(b, ok), &h)| ok && !fits(b, h)).count();
    let runs_on = bases.len() > refr.len();
    if plain == 0 || plain == 1 && !(runs_on && errors.contains(&0)) {
        return OwnEvent::None;
    }
    // The best placements: the fewest mismatches, one at most (an error) for the bases
    // to place the change. Inside the discrimination window when one of them is (RJ-22:
    // an insert that can sit there is never a separate event); a separate event at the
    // leftmost when none is.
    // Of the given allele's own length, a read no placement fits is the windows' to
    // judge (the ALT, or another allele of its length), as one placed inside is, its
    // bases placing no change; background is judged as written.
    let Some(&fewest) = errors.iter().min().filter(|&&e| e <= 1) else {
        return match (background, own_length) {
            (true, _) => OwnEvent::AsWritten,
            (false, true) => OwnEvent::None,
            (false, false) => OwnEvent::Inside,
        };
    };
    // RJ-7 and RJ-8's window: a pure indel's discrimination window (its tract with a
    // base each side); any other event's windows
    let (w_lo, w_hi) = match super::window::allele_kind(&v.ref_allele, &v.alt_allele) {
        Some(AlleleKind::Insertion | AlleleKind::Deletion) => super::window::discrimination_window(v),
        _ => span,
    };
    let inside_at = |p: usize| {
        let at = x0 + p as i64;
        if n > 0 { w_lo < at && at < w_hi } else { at < w_hi && at + m as i64 > w_lo }
    };
    // An insert with no readable base does not slide (RJ-20): it sits where written.
    let mut inserted = Vec::new();
    let mut qp = 0usize;
    for op in record.cigar().iter() {
        match op {
            Cigar::Ins(k) => {
                inserted.extend(qp..qp + *k as usize);
                qp += *k as usize;
            }
            Cigar::Match(k) | Cigar::Equal(k) | Cigar::Diff(k) | Cigar::SoftClip(k) => qp += *k as usize,
            _ => {}
        }
    }
    let readable = |q: &usize| quals.get(*q).is_some_and(|&x| x >= min_baseq) && seq.get(*q).is_some_and(|b| b.to_ascii_uppercase() != WILD);
    let unreadable = n > 0 && !inserted.is_empty() && !inserted.iter().any(readable);
    let written = ops.iter().find(|o| o.1 > 0).and_then(|o| usize::try_from(o.0 - x0).ok()).filter(|&w| w < bases.len());
    let mut best: Box<dyn Iterator<Item = usize>> = match written.filter(|_| unreadable) {
        Some(w) => Box::new(std::iter::once(w)),
        None => Box::new((0..errors.len()).filter(|&p| errors[p] == fewest)),
    };
    // Equivalent placements spell one haplotype: a separate event is applied to the
    // windows at one inside them where there is one, else at the leftmost.
    let (w_left, w_right) = (win.pair.0.left.unwrap_or(i64::MIN), win.pair.0.right.unwrap_or(i64::MAX));
    let p = written.filter(|_| unreadable).unwrap_or_else(|| {
        let fewest_at = (0..errors.len()).filter(|&p| errors[p] == fewest);
        fewest_at.clone().find(|&p| x0 + p as i64 > w_left && x0 + (p as i64) < w_right).or(fewest_at.min()).unwrap_or(0)
    });
    let (at, inside) = (x0 + p as i64, best.any(inside_at));
    // Of the ALT's own length inside the window, the read's bases spell the ALT (one
    // error aside; the windows judge it) or another allele of its length (RJ-21:
    // another allele, partial evidence).
    let spells_alt = || {
        let Some(ev) = event(start, &reference, v) else { return true };
        let hap = ev.hap.get(r0..).unwrap_or_default();
        bases[..reach].iter().zip(hap).filter(|&(&(b, ok), &h)| ok && !fits(b, h)).count() <= 1
    };
    let left = win.pair.0.left.map_or(q0 as i64, |l| q0 as i64 + (l - x0));
    // Inside a pure indel's window a recurring insert is another length of its run, an
    // allele at the locus, never background; inside a complex event's windows a
    // germline indel beside the event keeps its reads' call (RJ-7).
    let pure = matches!(super::window::allele_kind(&v.ref_allele, &v.alt_allele), Some(AlleleKind::Insertion | AlleleKind::Deletion));
    let background = background && !(inside && pure);
    match (inside, background || own_length && spells_alt()) {
        (true, true) if background => OwnEvent::AsWritten,
        (true, true) => OwnEvent::Given { left },
        (true, false) => OwnEvent::Inside,
        (false, _) => OwnEvent::Separate {
            at,
            ins: if n > 0 { bases.get(p..(p + m).min(bases.len())).unwrap_or_default().iter().map(|&(b, _)| b).collect() } else { Vec::new() },
            del: if n < 0 { m } else { 0 },
        },
    }
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

/// How a read's windows were placed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Placement {
    /// From the reading flank's end base, aligned.
    Flank,
    /// From the other flank, the reading flank in the read's clip or under its own
    /// deletion.
    OtherFlank,
    /// From the other flank through aligned bases: the aligner wrote the reading
    /// flank at the event's other copy.
    ThroughAligned,
    /// From the other flank, the reading flank in the clip or under its own deletion,
    /// as the allele its bases favour sits (the ALT when `alt`), which they do not
    /// fit: another allele, closer to that one.
    Unfit { alt: bool },
}

/// The REF and ALT windows read from the read, and how they were placed. A read
/// that aligns the reading flank's end base (the left flank's first base, or the
/// right flank's last) reads both windows from it: the windows share it, so the
/// read is placed truly whatever its own gap. Otherwise it is placed from the other
/// flank's aligned end base as the allele its bases favour would sit there: each
/// allele's window placed from that anchor at its own length, the read's confident
/// bases (at or above `min_baseq`, not N) must mismatch one allele's by at least one
/// base read at `min_baseq` more than the other's. Bases below min BQ weigh nothing,
/// so a REF molecule behind a poor clip shows no placement, whatever its clip's
/// bases happen to fit. With the flank in the read's soft clip (the aligner clipped
/// the event and the flank) or under its own deletion (it used a microhomology at
/// the event's far end), the favoured allele counts where it fits, else the read
/// carries another allele; with the flank aligned elsewhere (an aligner writes a
/// carrier's flank over the event's other copy, as mismatches), only where it fits:
/// a REF read there shows the event's bases, not the flank. A read that merely
/// starts (or ends) inside the windows holds them from neither flank: depth only,
/// for either allele alike. None when the read cannot be placed or does not hold
/// both windows' held bases. Only the query bases `[first, after)` are read, and
/// windows built in spliced space are placed through `map`.
#[allow(clippy::too_many_arguments)]
fn read_pair(
    record: &Record,
    seq: &[u8],
    quals: &[u8],
    min_baseq: u8,
    rw: &Window,
    aw: &Window,
    map: &SpliceMap,
    readable: (usize, usize),
    fit: &dyn Fn(&Reading, &Reading, &Window, &Window, u8) -> bool,
    left: Option<i64>,
    unanchored: bool,
) -> Option<(Reading, Reading, Placement)> {
    // a read that does not read its left flank base is not placed from it
    let at_left = if unanchored {
        None
    } else {
        left.or_else(|| rw.left.and_then(|g| find_read_pos(record, map.to_genome(g))).map(|q| q as i64))
    };
    let at_right = rw.right.and_then(|g| find_read_pos(record, map.to_genome(g - 1))).map(|q| q as i64 + 1);
    let tail = rw.reads_from_right();
    let read = |w: &Window, s: i64| read_at(seq, quals, min_baseq, w, s, tail, readable);
    // Both windows from one placement of the reading flank's end base: their shared
    // start (left), or their shared end (right; each starts its own length before it).
    let from = |anchor: i64| {
        let start = |w: &Window| if tail { anchor - w.allele.len() as i64 } else { anchor };
        Some((read(rw, start(rw))?, read(aw, start(aw))?))
    };
    // The allele confident bases favour by `evidence` (log10, ALT over REF): by at least
    // `bases` bases read at min BQ, and more than none (at --min-baseq 0 a base weighs
    // nothing); a tie favours neither.
    let one = one_base_evidence(min_baseq);
    let favoured = |evidence: f64, bases: f64| {
        if evidence > 0.0 && evidence >= bases * one {
            Some(true)
        } else if evidence < 0.0 && evidence <= -bases * one {
            Some(false)
        } else {
            None
        }
    };
    // The read placed from the other flank's `anchor` as the allele `alt` names sits
    // there: the reading flank's end base where that allele puts it.
    let as_placed = |anchor: i64, alt: bool| -> Option<(Reading, Reading, Placement)> {
        let w = if alt { aw } else { rw };
        let s = if tail { anchor + w.allele.len() as i64 } else { anchor - w.allele.len() as i64 };
        let (r, a) = from(s)?;
        let (first, after) = aligned_query_range(record)?;
        let flank_end = if tail { rw.right.map(|g| g - 1) } else { rw.left };
        let deleted = flank_end.is_some_and(|g| deleted_anchor_query(record, map.to_genome(g)).is_some());
        // Held: the flank's end base lies in the clip on the reading side, or under the
        // read's own deletion; else it is aligned elsewhere.
        let held = deleted || if tail { s > after as i64 } else { s < first as i64 };
        let fits = if alt { fit(&a, &r, aw, rw, 1) } else { fit(&r, &a, rw, aw, 0) };
        let placement = match (fits, held) {
            (true, true) => Placement::OtherFlank,
            (true, false) => Placement::ThroughAligned,
            (false, true) => Placement::Unfit { alt },
            (false, false) => return None,
        };
        Some((r, a, placement))
    };
    let other = if tail { at_left } else { at_right };
    if let Some(s) = if tail { at_right } else { at_left } {
        let (r, a) = from(s)?;
        // Confirmed from the other flank, read whole. A fit of the shorter allele rests on
        // its few bases at the reading flank (one past the tract decides, RJ-27): where
        // the aligner wrote the longer allele's other copy over the reading flank (an
        // insert, or a deleted sequence, whose end resembles the flank, aligned onto it
        // as mismatches, the rest clipped), those bases fit the shorter allele and the
        // read's bases before them are the longer allele's. Read from the other flank at
        // the longer allele's length (the shorter extended on the reading side by the
        // reference beyond it, which both share), its confident bases then favour the
        // longer allele by two bases, as a shifted read must show. A read holding the
        // longer allele's window from there and reaching past it, as one read from the
        // reading flank does, is that allele; one starting inside the event shows only
        // the other junction, which the shorter allele's molecules starting there have no
        // counterpart of: depth only, never the shorter allele (RJ-23's one reading
        // flank, operator 2026-10-16). Only a read aligned straight across the windows:
        // a read whose own indel there placed them holds the other flank shifted by it.
        // A fit of the longer allele reads its whole event and stands.
        let fitted = (fit(&r, &a, rw, aw, 0), fit(&a, &r, aw, rw, 1));
        let longer_alt = aw.allele.len() > rw.allele.len();
        let (fits_longer, fits_shorter) = if longer_alt { (fitted.1, fitted.0) } else { (fitted.0, fitted.1) };
        let straight = left.is_none()
            && rw.left.zip(rw.right).is_some_and(|(lo, hi)| own_indels(record, (map.to_genome(lo), map.to_genome(hi - 1) + 1)).is_empty());
        // ... and holding the other flank where its own bases say: past it they continue
        // the reference as aligned (one error aside) and fit it shifted by the alleles'
        // length difference worse, by two bases. A carrier the aligner wrote without its
        // gap holds the other flank shifted; past it a repeat whose period divides the
        // length change fits either way and anchors nothing.
        let shift = rw.allele.len().abs_diff(aw.allele.len()) as i64;
        let anchored = !tail
            && other.and_then(|q| usize::try_from(q).ok()).is_some_and(|q| {
                // confident bases compared and missed over the first OTHER_FLANK_CHECK
                // reference bases past the other flank, the read's bases from `q + read`
                // against the reference's from `past`
                let misses = |read: usize, past: usize| {
                    (0..OTHER_FLANK_CHECK)
                        .filter_map(|i| Some((q + read + i, *rw.past.get(past + i)?)))
                        .filter(|&(x, _)| x < readable.1 && quals[x] >= min_baseq && seq[x].to_ascii_uppercase() != WILD)
                        .fold((0, 0), |(n, m), (x, h)| (n + 1, m + usize::from(!fits(seq[x].to_ascii_uppercase(), h))))
                };
                let (seen, missed) = misses(0, 0);
                // the read's bases there shifted by the length change, either way, fit
                // worse; a shift the read does not reach that far shows nothing
                let d = usize::try_from(shift).unwrap_or(0);
                let worse = |(n, m): (usize, usize)| n < OTHER_FLANK_CHECK / 2 || m >= missed + 2;
                seen >= OTHER_FLANK_CHECK / 2 && missed <= 1 && worse(misses(0, d)) && worse(misses(d, 0))
            });
        let longest = rw.allele.len().max(aw.allele.len());
        let whole = |w: &Window| -> Option<Vec<u8>> {
            let n = longest - w.allele.len();
            let outer = if tail { w.outer.get(..n)? } else { w.outer.get(w.outer.len().checked_sub(n)?..)? };
            Some(if tail { [w.allele.as_slice(), outer].concat() } else { [outer, w.allele.as_slice()].concat() })
        };
        let longer_seen = straight
            && anchored
            && rw.allele.len() != aw.allele.len()
            && !fits_longer
            && other.zip(whole(rw)).zip(whole(aw)).is_some_and(|((anchor, rh), ah)| {
                let start = if tail { anchor } else { anchor - longest as i64 };
                let weight = |h: &[u8]| confident_weight(seq, quals, min_baseq, h, start, readable);
                favoured(weight(&rh) - weight(&ah), 2.0) == Some(longer_alt)
            });
        if longer_seen {
            match other.and_then(|anchor| as_placed(anchor, longer_alt)) {
                Some(placed) if matches!(placed.2, Placement::OtherFlank | Placement::ThroughAligned) => {
                    trace!("carrier: read whole from the other flank its bases are the longer allele's → placed from there");
                    return Some(placed);
                }
                _ if fits_shorter => {
                    trace!("carrier: it fits the shorter allele, but read whole from the other flank its bases are the longer allele's, whose window it does not hold → depth only");
                    return None;
                }
                _ => {}
            }
        }
        return Some((r, a, Placement::Flank));
    }
    // Placed from the other flank: each allele's window starts there (right reading)
    // or ends there (left reading), at its own length.
    let anchor = other?;
    let start = |w: &Window| if tail { anchor } else { anchor - w.allele.len() as i64 };
    let weight = |w: &Window| confident_weight(seq, quals, min_baseq, &w.allele, start(w), readable);
    let evidence = weight(rw) - weight(aw);
    let Some(alt) = favoured(evidence, 1.0) else {
        trace!("carrier: placed from the other flank, its confident bases favour neither allele ({:.2}) → depth only", evidence);
        return None;
    };
    as_placed(anchor, alt)
}

/// A variant's alleles over the whole-read reach, for reading a read whole
/// ([`whole_fit`]). They depend on the variant and its guard alone, so they are built
/// once per variant ([`CarrierGuard::whole`]).
pub(crate) struct WholeHaps {
    key: (i64, String, String),
    start: i64,
    /// The reference as fetched: a read's changes from it ([`aligned_changes`]).
    plain: Vec<u8>,
    /// REF and the ALT, the guard's masks laid as the windows lay them.
    reference: Vec<u8>,
    hap: Vec<u8>,
    /// Each with its event's bases open (any base): the other allele a read is weighed
    /// against must be shown by its length or its bases past the event.
    open_ref: Vec<u8>,
    open_hap: Vec<u8>,
    /// The windows: REF's `[lo, hi)`, the ALT's `[lo, alt_hi)`.
    lo: usize,
    hi: usize,
    alt_hi: usize,
    /// The other alleles of the ALT's length and of REF's at the locus, built on the
    /// reference as fetched, each taking an allele's mask where its base is one the
    /// mask takes (as the windows' other alleles do).
    others: Vec<Vec<u8>>,
    ref_others: Vec<Vec<u8>>,
}

/// The variant's alleles over the whole-read reach (prep's `whole_ref`, else the
/// prepared reference), laid out by `layout`'s masks and other alleles. None without a
/// reference holding the event and its flank.
fn whole_haps(v: &Variant, layout: &Layout) -> Option<WholeHaps> {
    let (start, plain) = whole_reference(v)?;
    let ev = event(start, &plain, v)?;
    let flank = ev.flank();
    let (lo, hi) = (ev.lo.checked_sub(flank)?, ev.hi + flank);
    let alt_hi = usize::try_from(hi as i64 + ev.d).ok()?;
    if hi > plain.len() || alt_hi > ev.hap.len() {
        return None;
    }
    // The guard's masks, laid as the windows lay them (`windows_on`): a masked base
    // takes the bases its mask takes (any, beside the event; the reference's and the
    // change's, inside the repeat it grew through), left of the given allele on both
    // alleles at its place, right of it on the ALT shifted by the length change; at a
    // base an MNP keeps on both; at a base a deletion deletes on REF alone; never on
    // the given allele's own bases.
    let (mut reference, mut hap) = (plain.clone(), ev.hap.clone());
    let mask = |b: u8, p: i64| IUPAC[(bit(b) | layout.background.iter().find(|&&(q, _)| q == p).map_or(ANY, |&(_, s)| s)) as usize];
    let ((o_lo, o_hi), (c_lo, c_hi)) = (ev.own, ev.core);
    for &p in &layout.masked {
        let Some(off) = usize::try_from(p - start).ok().filter(|&off| off < plain.len()) else { continue };
        let (r, h) = if off < o_lo {
            (true, Some(off))
        } else if off >= o_hi {
            (true, usize::try_from(off as i64 + ev.d).ok())
        } else if ev.d == 0 && off >= c_lo && off < c_hi && plain[off] == ev.hap[off] {
            (true, Some(off))
        } else {
            (ev.deletes, None)
        };
        if r {
            reference[off] = mask(reference[off], p);
        }
        if let Some(h) = h.filter(|&h| h < hap.len()) {
            hap[h] = mask(hap[h], p);
        }
    }
    let open = |h: &[u8], (e0, e1): (usize, usize)| -> Vec<u8> {
        h.iter().enumerate().map(|(i, &b)| if i >= e0 && i < e1 { WILD } else { b }).collect()
    };
    let alt_ev_hi = usize::try_from(ev.hi as i64 + ev.d).ok()?;
    let (open_ref, open_hap) = (open(&reference, (ev.lo, ev.hi)), open(&hap, (ev.lo, alt_ev_hi)));
    let kept = |list: &[OtherAllele], allele: &[u8]| -> Vec<Vec<u8>> {
        list.iter()
            .filter_map(|o| o.on(start, &plain))
            .filter(|h| h.len() == allele.len())
            .map(|h| h.iter().zip(allele).map(|(&o, &w)| if takes(w).count_ones() > 1 && fits(o, w) { w } else { o }).collect())
            .collect()
    };
    let (others, ref_others) = (kept(&layout.others, &hap), kept(&layout.ref_others, &reference));
    let key = (v.pos, v.ref_allele.clone(), v.alt_allele.clone());
    Some(WholeHaps { key, start, plain, reference, hap, open_ref, open_hap, lo, hi, alt_hi, others, ref_others })
}

/// The read's confident changes from the reference over its aligned bases inside the
/// whole-read reach (its own indels each count as one): where it is two or more, the
/// aligner's placement is no clean one.
fn aligned_changes(record: &Record, seq: &[u8], quals: &[u8], min_baseq: u8, w: &WholeHaps) -> usize {
    let (start, reference) = (w.start, &w.plain);
    let end = start + reference.len() as i64;
    let (mut changes, mut g, mut q) = (0usize, record.pos(), 0usize);
    for op in record.cigar().iter() {
        match *op {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                for i in 0..l as usize {
                    let at = g + i as i64;
                    if at >= start && at < end {
                        let (b, r) = (seq[q + i].to_ascii_uppercase(), reference[(at - start) as usize]);
                        changes += usize::from(quals.get(q + i).is_some_and(|&x| x >= min_baseq) && b != WILD && b != r);
                    }
                }
                g += l as i64;
                q += l as usize;
            }
            Cigar::Ins(l) => {
                changes += usize::from(g >= start && g < end);
                q += l as usize;
            }
            Cigar::Del(l) | Cigar::RefSkip(l) => {
                changes += usize::from(matches!(*op, Cigar::Del(_)) && g >= start && g < end);
                g += l as i64;
            }
            Cigar::SoftClip(l) => q += l as usize,
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    changes
}

/// Mismatches counted per placement when a read is read whole: only none, one, or
/// more decide ([`whole_fit`]).
const WHOLE_CAP: usize = 2;

/// The allele a read's bases fit read whole (`readable`, without its own alignment's
/// indels): placed ungapped on each allele's haplotype over the whole-read reach, the
/// whole read on it, at every placement holding that allele's windows (event and flank),
/// its confident mismatches (at or above `min_baseq`, not N; the guard's masks take
/// their bases, as in the windows). An allele fits
/// at a placement reading its windows exactly (the flank bases anchor it there) and
/// missing at most one base past them (an error, or a phased change); ALT (`Some(Some(true))`)
/// or REF where one allele fits so and the other, and every other allele of either
/// length at the locus, misses by more wherever it holds as much of its windows as the
/// shorter allele's (its event's bases open: a read must show the other's length or
/// its bases past the event). `Some(None)`, depth only, where it fits the shorter allele's
/// windows exactly and the longer allele only from inside its event (RJ-23); None
/// otherwise.
fn whole_fit(w: &WholeHaps, seq: &[u8], quals: &[u8], min_baseq: u8, (first, after): (usize, usize)) -> Option<Option<bool>> {
    let (lo, hi, alt_hi) = (w.lo, w.hi, w.alt_hi);
    let (read, q) = (seq.get(first..after)?, quals.get(first..after)?);
    let n = read.len();
    let miss = |hap: &[u8], o: usize, i: usize| {
        let b = read[i].to_ascii_uppercase();
        q[i] >= min_baseq && b != WILD && !fits(b, hap[o + i])
    };
    // over the placements holding [w0, w1): the fewest mismatches (up to WHOLE_CAP), and
    // the fewest where none lies inside [e0, e1) (the windows read exactly)
    let best = |hap: &[u8], (w0, w1): (usize, usize), (e0, e1): (usize, usize)| -> Option<(usize, Option<usize>)> {
        // the whole read on the held reference: no base of it goes unread
        let (from, to) = ((w1 as i64 - n as i64).max(0), (w0 as i64).min(hap.len() as i64 - n as i64));
        if from > to {
            return None;
        }
        let (mut any, mut clean) = (WHOLE_CAP, None::<usize>);
        for o in from as usize..=to as usize {
            let inside = (e0.saturating_sub(o).min(n)..e1.saturating_sub(o).min(n)).any(|i| miss(hap, o, i));
            let mut m = 0;
            for i in 0..n {
                if miss(hap, o, i) {
                    m += 1;
                    if m == WHOLE_CAP {
                        break;
                    }
                }
            }
            any = any.min(m);
            if !inside {
                clean = Some(clean.map_or(m, |c| c.min(m)));
            }
            if any == 0 && clean == Some(0) {
                break;
            }
        }
        Some((any, clean))
    };
    let (_, r_clean) = best(&w.reference, (lo, hi), (lo, hi))?;
    let (_, a_clean) = best(&w.hap, (lo, alt_hi), (lo, alt_hi))?;
    // The other allele is weighed wherever the read holds as much of its windows as the
    // shorter allele's, from either flank (as a read holding the windows must): a read
    // too short to hold the longer allele's windows does not rule it out.
    let short = (hi - lo).min(alt_hi - lo);
    let held_short = |hap: &[u8], w1: usize| -> Option<usize> {
        let left = best(hap, (lo, lo + short), (lo, lo + short)).map(|x| x.0);
        let right = best(hap, (w1 - short, w1), (w1 - short, w1)).map(|x| x.0);
        left.into_iter().chain(right).min()
    };
    let r_any = held_short(&w.open_ref, hi).unwrap_or(usize::MAX);
    let a_any = held_short(&w.open_hap, alt_hi).unwrap_or(usize::MAX);
    // the other alleles of each length at the locus, each over its own allele's windows
    let others = |list: &[Vec<u8>], win: (usize, usize)| -> usize {
        list.iter().filter_map(|h| best(h, win, win)).map(|x| x.0).min().unwrap_or(usize::MAX)
    };
    let rest = others(&w.others, (lo, alt_hi)).min(others(&w.ref_others, (lo, hi)));
    let fits = |clean: Option<usize>, other: usize| clean.is_some_and(|c| c <= 1 && other > c && rest > c);
    // A read fitting the shorter allele's whole windows exactly as written, and the
    // longer allele exactly only where it starts inside the longer allele's event (it
    // holds not its whole windows), shows only the other junction: depth only, never
    // the shorter allele (RJ-23).
    let (shorter_whole, longer_whole, longer_part) = if hi - lo <= alt_hi - lo {
        (r_clean, a_clean, held_short(&w.hap, alt_hi))
    } else {
        (a_clean, r_clean, held_short(&w.reference, hi))
    };
    if fits(a_clean, r_any) {
        Some(Some(true))
    } else if fits(r_clean, a_any) {
        Some(Some(false))
    } else if shorter_whole == Some(0) && longer_whole.is_none() && longer_part == Some(0) {
        Some(None)
    } else {
        None
    }
}

/// The weight (log10) of the read's confident bases (at or above `min_baseq`, not
/// N) that mismatch `allele` placed at query position `s`, each by its quality, over
/// the readable bases `[first, after)`. A masked window base mismatches nothing.
pub(crate) fn confident_weight(seq: &[u8], quals: &[u8], min_baseq: u8, allele: &[u8], s: i64, (first, after): (usize, usize)) -> f64 {
    allele
        .iter()
        .enumerate()
        .filter_map(|(i, &h)| {
            let q = usize::try_from(s + i as i64).ok().filter(|&q| q >= first && q < after)?;
            let (b, bq) = (seq.get(q)?.to_ascii_uppercase(), *quals.get(q)?);
            (b != WILD && bq >= min_baseq && !fits(b, h)).then(|| one_base_evidence(bq))
        })
        .sum()
}

/// One allele's window read with the window starting at query position `s` (which
/// may lie before the read): its held bases (its last `hold` when `tail`), which must lie in `[first, after)`, then
/// on through the rest of the window as far as the read reaches, whose mismatches
/// count too (the quality and span stay the held bases').
#[allow(clippy::too_many_arguments)]
fn read_at(
    seq: &[u8],
    quals: &[u8],
    min_baseq: u8,
    w: &Window,
    s: i64,
    tail: bool,
    (first, after): (usize, usize),
) -> Option<Reading> {
    let (n, hold) = (w.allele.len() as i64, w.hold as i64);
    let off = if tail { n - hold } else { 0 };
    let (a, b) = (s + off, s + off + hold);
    if a < first as i64 || b > after as i64 {
        return None;
    }
    let (a, b, o) = (a as usize, b as usize, off as usize);
    let mut r = score(&seq[a..b], &quals[a..b], min_baseq, &w.allele[o..o + w.hold], (a, b));
    r.start = s;
    let (lo, hi) = if tail { (s.max(first as i64) as usize, a) } else { (b, (s + n).min(after as i64) as usize) };
    if lo < hi {
        let h = (lo as i64 - s) as usize;
        let more = score(&seq[lo..hi], &quals[lo..hi], min_baseq, &w.allele[h..h + (hi - lo)], (lo, hi));
        r.mismatches += more.mismatches;
        r.mismatch_weight += more.mismatch_weight;
        r.misses.extend(more.misses);
        r.read = (r.read.0.min(lo), r.read.1.max(hi));
    }
    Some(r)
}

/// Bases a shift past the windows may take ([`shifted_better`]): another allele's other
/// length, a repeat unit or a few.
pub(crate) const SHIFT_REACH: i64 = 8;
/// Read bases a shift is scored over.
const SHIFT_SPAN: usize = 12;

/// Whether the read's confident bases in `[from, to)` fit `expected` (an allele's bases
/// there, then the reference past the windows) shifted by up to [`SHIFT_REACH`] bases
/// better, by two, than as they stand, where as they stand two or more mismatch: the
/// read is another allele, whose other length shifts its further bases. A change beside
/// the event mismatches at its own bases under any shift.
fn shifted_better(seq: &[u8], quals: &[u8], min_baseq: u8, (from, to): (usize, usize), expected: &[u8]) -> bool {
    let mismatches = |s: i64| {
        let (q0, e0) = if s >= 0 { (from + s as usize, 0) } else { (from, s.unsigned_abs() as usize) };
        let (mut seen, mut miss) = (0usize, 0usize);
        for (i, &h) in expected.iter().skip(e0).take(SHIFT_SPAN).enumerate() {
            let q = q0 + i;
            let (Some(&b), Some(&bq)) = (seq.get(q).filter(|_| q < to), quals.get(q)) else { break };
            let b = b.to_ascii_uppercase();
            if b == WILD || bq < min_baseq {
                continue;
            }
            seen += 1;
            miss += usize::from(!fits(b, h));
        }
        (miss, seen)
    };
    let (standing, _) = mismatches(0);
    standing >= 2
        && (1..=SHIFT_REACH).flat_map(|k| [k, -k]).any(|s| {
            let (miss, seen) = mismatches(s);
            seen >= 3 && miss + 2 <= standing
        })
}

/// A base's error probability from its quality, kept away from 0 and from a
/// coin flip among four bases.
fn error_probability(q: u8) -> f64 {
    10f64.powf(-(q.min(60) as f64) / 10.0).clamp(1e-6, 0.75)
}

/// The evidence one base read at quality `q` gives for the allele it matches over
/// one it does not (log10): the least an ALT call must show.
pub(crate) fn one_base_evidence(q: u8) -> f64 {
    let e = error_probability(q);
    ((1.0 - e) / (e / 3.0)).log10()
}


/// Mismatches of read bases against haplotype bases position by position, masked read
/// bases (below `min_baseq`, or N) matching anything, and a masked haplotype base the
/// bases its code takes.
pub(crate) fn score(bases: &[u8], quals: &[u8], min_baseq: u8, hap: &[u8], span: (usize, usize)) -> Reading {
    let (mut mismatches, mut masked, mut had_n, mut mismatch_weight) = (0, 0, false, 0.0);
    let mut misses = Vec::new();
    for (i, ((&b, &q), &h)) in bases.iter().zip(quals).zip(hap).enumerate() {
        let b = b.to_ascii_uppercase();
        had_n |= b == WILD;
        if b != WILD && !fits(b, h) {
            mismatch_weight += one_base_evidence(q);
        }
        if b == WILD || q < min_baseq {
            masked += 1;
        } else if !fits(b, h) {
            mismatches += 1;
            misses.push(span.0 + i);
        }
    }
    Reading { mismatches, masked, had_n, span, read: span, mismatch_weight, start: span.0 as i64, misses }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;

    /// An unpaired read of `cigar` at `pos`.
    fn clipped(cigar: Vec<Cigar>, pos: i64) -> Record {
        let cigar = CigarString(cigar);
        let len: usize = cigar
            .iter()
            .filter(|op| matches!(op, Cigar::Match(_) | Cigar::Ins(_) | Cigar::SoftClip(_)))
            .map(|op| op.len() as usize)
            .sum();
        let mut rec = Record::new();
        let seq: Vec<u8> = (0..len).map(|i| b"ACGT"[i % 4]).collect();
        rec.set(b"r", Some(&cigar), &seq, &vec![30; len]);
        rec.set_pos(pos);
        rec
    }

    fn readable(rec: &Record, edges: Option<Vec<i64>>) -> Option<(usize, usize)> {
        let found = edges.map(|e| ClipEdges::new(move || (e.clone(), Vec::new())));
        readable_range(rec, &ReadRules { clip_edges: found.as_ref(), reference: None, spliced: None, guard: None, tolerance: None })
    }

    #[test]
    fn an_rna_clip_is_read_unless_it_reaches_an_exon_edge() {
        // 90M10S at 200: aligned [200, 290), the clip reaching [290, 300).
        let rec = clipped(vec![Cigar::Match(90), Cigar::SoftClip(10)], 200);
        assert_eq!(readable(&rec, None), Some((0, 100)), "DNA reads every clip");
        assert_eq!(readable(&rec, Some(vec![])), Some((0, 100)), "no edge near");
        assert_eq!(readable(&rec, Some(vec![295])), Some((0, 90)), "an edge inside the clip's reach");
        assert_eq!(readable(&rec, Some(vec![290])), Some((0, 90)), "the clip starts at the edge");
        assert_eq!(readable(&rec, Some(vec![286])), Some((0, 90)), "aligned four bases past an edge");
        assert_eq!(readable(&rec, Some(vec![285])), Some((0, 90)), "aligned five bases past an edge");
        assert_eq!(readable(&rec, Some(vec![284])), Some((0, 100)), "six bases past: the read's own clip");
        assert_eq!(readable(&rec, Some(vec![301])), Some((0, 100)), "past the clip's reach");
        // 8S92M at 508: a leading clip reaching [500, 508), an acceptor at 500.
        let lead = clipped(vec![Cigar::SoftClip(8), Cigar::Match(92)], 508);
        assert_eq!(readable(&lead, Some(vec![500])), Some((8, 100)), "acceptor inside the leading clip's reach");
        assert_eq!(readable(&lead, Some(vec![513])), Some((8, 100)), "aligned five bases into the exon");
        assert_eq!(readable(&lead, Some(vec![514])), Some((0, 100)), "six bases in: the read's own clip");
        // Each clip on its own.
        let both = clipped(vec![Cigar::SoftClip(5), Cigar::Match(90), Cigar::SoftClip(5)], 205);
        assert_eq!(readable(&both, Some(vec![297])), Some((0, 95)), "only the trailing clip reaches the edge");
    }

    fn var(ctx: &str, pos: i64, r: &str, a: &str) -> Variant {
        Variant::new(
            "1".into(), pos, r.into(), a.into(), "COMPLEX".into(),
            None, 0, 0, None, None, Some((0, ctx.into())), None,
        )
    }

    fn only_pair(w: &Windows) -> (&Window, &Window) {
        (&w.pair.0, &w.pair.1)
    }

    #[test]
    fn a_pure_indel_grows_through_the_repeat_holding_its_bases_only() {
        let grown = |ctx: &str, pos: usize, r: &str, a: &str| {
            let ev = event(0, ctx.as_bytes(), &var(ctx, pos as i64, r, a)).unwrap();
            (ev.lo, ev.hi)
        };
        // An insert of other bases after the first G of a GGG run on REF: the run is
        // beside the change, the event its junction alone.
        let ctx = "CTCATCTC".to_owned() + "GGG" + "TCTCACTC";
        assert_eq!(grown(&ctx, 8, "G", "GATCCACT"), (9, 9));
        // Two units deleted between two G runs: the G run the ALT joins is beside it.
        let ctx = "CTCATCAG".to_owned() + "TC" + "GGGGCTCATC";
        assert_eq!(grown(&ctx, 7, "GTC", "G"), (8, 10));
        // An insert of GAAA and GAAG units where a (GAAA)n run meets a (GAAG)n run:
        // the event is the repeat it lengthens, both runs on the ALT.
        let ctx = "CTCTTCTC".to_owned() + &"GAAA".repeat(5) + &"GAAG".repeat(5) + "CTTCTCTC";
        let (lo, hi) = grown(&ctx, 8 + 20 - 1, "A", "AGAAAGAAG");
        assert_eq!((lo, hi), (8, 48), "the whole (GAAA)5(GAAG)5 repeat");
    }

    #[test]
    fn windows_are_the_event_with_two_flank_bases_held_to_one_length() {
        // GGAC [TA] GTCGTT -> GGAC [GCC] GTCGTT: each allele's window is the event
        // with two flank bases on each side; a read holds as many bases of each as
        // the shorter (REF) window has, from the left flank, both anchored on the
        // same flank bases.
        let w = windows(&var("GGACTAGTCGTT", 4, "TA", "GCC"), &Layout::default()).unwrap();
        let (r, a) = only_pair(&w);
        assert_eq!(a.allele, b"ACGCCGT".to_vec());
        assert_eq!(r.allele, b"ACTAGT".to_vec());
        assert_eq!((r.hold, a.hold, r.tail), (6, 6, false));
        assert_eq!((r.left, r.right), (Some(2), Some(8)));
        assert_eq!((a.left, a.right), (Some(2), Some(8)));
    }

    #[test]
    fn the_windows_run_through_bases_the_alleles_share_past_the_event() {
        // GGAC [TTTA] TTGACGTCAG -> GGAC [G] TTGACGTCAG: past the one-base ALT its flank
        // TT repeats REF's next two bases, and the alleles differ at the third (G vs
        // A). Both windows run on to two bases from there: a read with only REF's
        // first base substituted shows A where the ALT's flank has G.
        let w = windows(&var("GGACTTTATTGACGTCAG", 4, "TTTA", "G"), &Layout::default()).unwrap();
        let (r, a) = only_pair(&w);
        assert_eq!(a.allele, b"ACGTTGA".to_vec());
        assert_eq!(r.allele, b"ACTTTATTGA".to_vec());
        assert_eq!((r.hold, a.hold), (7, 7));
        assert_eq!((r.right, a.right), (Some(12), Some(12)));
        // No shared bases (the flank differs at once): the windows stay the event and
        // its two flank bases.
        let w = windows(&var("GGACTTTACCGACGTCAG", 4, "TTTA", "G"), &Layout::default()).unwrap();
        assert_eq!(only_pair(&w).1.allele, b"ACGCC".to_vec());
    }

    #[test]
    fn the_event_grows_through_runs_on_either_side() {
        // CAC>A inside CCCC A CCCC: a C deleted from each run, placeable anywhere
        // in it; the windows are anchored outside both runs.
        let w = windows(&var("TTGCCCCACCCCTTA", 6, "CAC", "A"), &Layout::default()).unwrap();
        let (r, a) = only_pair(&w);
        assert_eq!(r.left, Some(1)); // two flank bases before the first C run (3..)
        assert_eq!(r.allele, b"TGCCCCACCCCTT".to_vec());
        assert_eq!((r.hold, a.hold), (a.allele.len(), a.allele.len()));
    }

    #[test]
    fn the_event_grows_through_a_run_the_alt_continues() {
        // GTC AAA [TA] CGT -> GTC AAA [ACC] CGT: the ALT's A continues the run, so
        // an inserted A can sit anywhere in it; the left anchor moves before it.
        let w = windows(&var("GCGTCAAATACGTGG", 8, "TA", "ACC"), &Layout::default()).unwrap();
        let (r, _) = only_pair(&w);
        assert_eq!(r.left, Some(3)); // two flank bases before the run (5..)
    }

    #[test]
    fn one_base_evidence_moves_with_the_minimum_quality() {
        // A Q20 base: (0.99 / (0.01 / 3)), about 2.47 log10; Q30 about 3.48; Q0 none.
        assert!((one_base_evidence(20) - 2.47).abs() < 0.01);
        assert!((one_base_evidence(30) - 3.48).abs() < 0.01);
        assert!(one_base_evidence(0).abs() < 1e-9);
    }

    #[test]
    fn low_quality_mismatches_weigh_little_and_n_nothing() {
        // One mismatch (G vs T): at Q30 it weighs one Q30 base, at Q5 far less; the
        // N and the matching bases weigh nothing.
        let full = score(b"ACGN", &[30, 30, 30, 30], 20, b"ACTA", (0, 4)).mismatch_weight;
        let low = score(b"ACGN", &[30, 30, 5, 30], 20, b"ACTA", (0, 4)).mismatch_weight;
        assert!((full - one_base_evidence(30)).abs() < 1e-9);
        assert!(low < 1.0, "{low}");
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
    fn a_long_allele_against_a_short_one_holds_the_short_window_whole() {
        // A 69-base REF replaced by TG: one pair, the whole short ALT window held and
        // as many REF bases from the same flank, the rest of the REF read on.
        let mut ctx = String::from("ACGTTGCA");
        let body: String = (0..70).map(|i| b"ACGGTTCA"[i % 8] as char).collect();
        ctx.push_str(&body);
        ctx.push_str("TGCATTGC");
        let w = windows(&var(&ctx, 7, &ctx[7..76], "TG"), &Layout::default()).unwrap();
        let (r, a) = only_pair(&w);
        assert_eq!((a.hold, r.hold), (a.allele.len(), a.allele.len()));
        assert!(r.allele.len() > LONG_EVENT);
        assert_ne!(&r.allele[..r.hold], &a.allele[..]);
        // Read from the right flank (an exon edge cuts the left), the REF's last
        // bases are held instead.
        let from_right = Layout { from_right: true, ..Layout::default() };
        let w = windows(&var(&ctx, 7, &ctx[7..76], "TG"), &from_right).unwrap();
        let (r, a) = only_pair(&w);
        assert!(r.tail && a.tail);
        assert_eq!(&r.allele[r.allele.len() - 2..], &a.allele[a.allele.len() - 2..]);
    }

    #[test]
    fn both_long_alleles_are_read_from_the_flank_nearer_the_first_difference() {
        // 70 bases replaced by 60 others (A/T against a C/G/T reference: they differ
        // from the first base): past 50 bases neither window fits a read, so both
        // hold one base past the first difference from the flank nearer it.
        let mut x = 11u32;
        let mut cgt = |n: usize| -> String {
            (0..n)
                .map(|_| {
                    x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                    ['C', 'G', 'T'][(x >> 16) as usize % 3]
                })
                .collect()
        };
        let ctx = format!("{}{}{}", cgt(40), cgt(70), cgt(40));
        let alt: String = (0..60).map(|i| if i % 3 == 0 { 'A' } else { 'T' }).collect();
        let v = var(&ctx, 40, &ctx[40..110], &alt);
        let w = windows(&v, &Layout::default()).unwrap();
        let (r, a) = only_pair(&w);
        let (start, reference) = prepared_reference(&v).unwrap();
        let ev = event(start, &reference, &v).unwrap();
        let from_left = ev.first - (ev.lo - FLANK) + 2;
        let from_right = (ev.hi + FLANK) - ev.last + 2;
        assert_eq!((r.tail, a.tail), (from_right < from_left, from_right < from_left));
        let held = from_left.min(from_right);
        assert!(held < LONG_EVENT);
        assert_eq!((r.hold, a.hold), (held, held));
    }

    #[test]
    fn a_short_reference_is_reported_on_its_side() {
        // Prep must hold the event with the flank and the guard's reach on each side:
        // an event ending a run that reaches the reference's left end is short on the
        // left only; with room on both sides, not short.
        let right: String = "TTCGTGGTCAGCTTGACCAGCATCGTACGTGATCGG".into();
        let ctx = format!("AAAAAAAAAA{right}");
        let v = var(&ctx, 9, "ATT", "GC");
        assert_eq!(reference_short(0, &ctx, &v), Some((true, false)));
        let ctx = format!("GCGTCGTAGCATCGATGCATCGTAGCTAGCATGCTAAAAAAAAAA{right}");
        let v = var(&ctx, 44, "ATT", "GC");
        assert_eq!(reference_short(0, &ctx, &v), None);
    }

    #[test]
    fn a_reference_ending_inside_the_flank_past_the_shared_bases_is_short() {
        // A 31-to-1 delins whose right flank repeats the deleted bases' first 26 and
        // then differs at P: the windows run to FLANK bases past P, so a reference
        // ending one base past P is short on the right, and one holding them is not.
        let (a, del, run) = (30usize, 30usize, 26usize);
        let mut x = 11u32;
        let mut ctx: Vec<u8> = (0..a + 1 + del + run + 2)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                b"CGT"[(x >> 16) as usize % 3]
            })
            .collect();
        (ctx[a - 1], ctx[a]) = (b'G', b'C');
        for u in 0..run {
            ctx[a + 1 + del + u] = ctx[a + 1 + u];
        }
        let p = a + 1 + del + run;
        if ctx[p] == ctx[p - del] {
            ctx[p] = if ctx[p] == b'G' { b'T' } else { b'G' };
        }
        let ctx = String::from_utf8(ctx).unwrap();
        let v = var(&ctx, a as i64, &ctx[a..a + 1 + del], "A");
        let ev = event(0, ctx.as_bytes(), &v).unwrap();
        assert_eq!((ev.lo, ev.hi), (a, a + 1 + del), "no repeat grows the event");
        assert_eq!(reference_short(0, &ctx[..p + 1], &v), Some((false, true)));
        assert!(windows(&var(&ctx[..p + 1], a as i64, &ctx[a..a + 1 + del], "A"), &Layout::default()).is_none());
        assert_eq!(reference_short(0, &ctx[..p + 2], &v), None);
        assert!(windows(&var(&ctx[..p + 2], a as i64, &ctx[a..a + 1 + del], "A"), &Layout::default()).is_some());
    }

    /// A read of `seq` aligned at `pos` by `cigar`, at Q30.
    fn aligned(seq: &[u8], cigar: Vec<Cigar>, pos: i64) -> Record {
        let mut rec = Record::new();
        rec.set(b"r", Some(&CigarString(cigar)), seq, &vec![30; seq.len()]);
        rec.set_pos(pos);
        rec
    }

    /// `n` copies of `rec`, each its own fragment (named `tag` and a number).
    #[test]
    fn a_reads_own_indels_inside_the_windows_are_another_unless_they_recur() {
        let read = |cigar: Vec<Cigar>| aligned(&[b'A'; 40], cigar, 100);
        let span = (110, 120);
        let ins = read(vec![Cigar::Match(15), Cigar::Ins(5), Cigar::Match(20)]);
        assert_eq!(own_indels(&ins, span), vec![(115, 5)], "an insertion inside");
        let del = read(vec![Cigar::Match(18), Cigar::Del(4), Cigar::Match(22)]);
        assert_eq!(own_indels(&del, span), vec![(118, -4)], "a deletion overlapping the end");
        let before = read(vec![Cigar::Match(10), Cigar::Ins(3), Cigar::Match(27)]);
        assert!(own_indels(&before, span).is_empty(), "an insertion at the span's first base is outside");
        let shifted = read(vec![Cigar::Match(12), Cigar::Ins(1), Cigar::Match(2), Cigar::Del(1), Cigar::Match(25)]);
        assert!(!foreign(&own_indels(&shifted, span), None), "an insertion and a deletion that cancel");
        assert!(foreign(&own_indels(&ins, span), None), "no background: another indel");
        let mut g = Guard::new(&var(&"A".repeat(200), 115, "A", "G"), vec![], None, 0);
        g.indels = vec![(115, 5)];
        assert!(!foreign(&own_indels(&ins, span), Some(&g)), "background: a separate event");
    }

    fn fragments(rec: &Record, n: usize, tag: &str) -> Vec<Record> {
        (0..n)
            .map(|i| {
                let mut r = rec.clone();
                r.set_qname(format!("{tag}{i}").as_bytes());
                r
            })
            .collect()
    }

    #[test]
    fn the_guard_reads_only_the_flank_the_reference_holds() {
        // A delins three bases from the reference's start (a contig start, where
        // prep cannot fetch the guard's reach): the left side reads three bases.
        let ctx = format!("GAG{}{}", "A".repeat(4), "GAGGAGAGGAAGAGGGAGAAGGAGAGAAGGAGAGGA");
        let v = var(&ctx, 3, "AAAA", "TC");
        let hap = format!("GAGTC{}", &ctx[7..]);
        let mut seq = hap.as_bytes().to_vec();
        seq[5 + 4] = b'C';
        let carrier = aligned(&seq, vec![Cigar::Match(3), Cigar::Ins(2), Cigar::Del(4), Cigar::Match(36)], 0);
        let reads = fragments(&carrier, 4, "a");
        let g = guard_for(&reads, &v, &|_| true, 20, &[]).expect("the change past the event");
        // the C four bases past the deleted A's (genome 11), named with the ALT
        assert_eq!(g.masked, vec![11]);
        assert!(g.larger.is_some());
        // Two fragments, each read by two overlapping mates: under the fragments needed.
        let mates = [fragments(&carrier, 2, "m"), fragments(&carrier, 2, "m")].concat();
        assert!(guard_for(&mates, &v, &|_| true, 20, &[]).is_none(), "mates are one fragment");
    }

    #[test]
    fn the_guard_masks_background_and_isolated_changes_and_names_the_alt_reads_own() {
        // 30 A/G bases, a 28-base delins to TTCTCT, 40 A/G bases. Four carriers show
        // a C three and six bases past the event (grown through the runs at its
        // edges); with `snp`, four REF reads show the C six bases past it too.
        let mut x = 7u32;
        let ctx: String = (0..98)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                if x >> 16 & 1 == 0 { 'A' } else { 'G' }
            })
            .collect();
        let v = var(&ctx, 30, &ctx[30..58], "TTCTCT");
        let (start, reference) = prepared_reference(&v).unwrap();
        let hi = event(start, &reference, &v).unwrap().hi;
        assert!(start == 0 && hi + 7 < ctx.len());
        let with_c = |seq: &str, at: &[usize]| -> Vec<u8> {
            let mut b = seq.as_bytes().to_vec();
            at.iter().for_each(|&i| b[i] = b'C');
            b
        };
        // past the event, a carrier's base for genome position g sits at g - 22
        let carrier_hap = format!("{}TTCTCT{}", &ctx[..30], &ctx[58..]);
        let carrier = aligned(
            &with_c(&carrier_hap, &[hi + 3 - 22, hi + 6 - 22]),
            vec![Cigar::Match(30), Cigar::Ins(6), Cigar::Del(28), Cigar::Match(40)],
            0,
        );
        // A run on the ALT reads alone: changes on two neighbouring bases (an indel
        // or MNP beside the event), another allele, neither masked nor named.
        let run = aligned(
            &with_c(&carrier_hap, &[hi + 3 - 22, hi + 4 - 22]),
            vec![Cigar::Match(30), Cigar::Ins(6), Cigar::Del(28), Cigar::Match(40)],
            0,
        );
        let plain = aligned(ctx.as_bytes(), vec![Cigar::Match(98)], 0);
        let reads = [fragments(&run, 4, "a"), fragments(&plain, 4, "r")].concat();
        assert!(guard_for(&reads, &v, &|_| true, 20, &[]).is_none(), "a run is not masked");
        // With an isolated change further on, that change is masked, and no larger
        // allele is named: it would leave out the run every ALT read carries.
        let both = aligned(
            &with_c(&carrier_hap, &[hi + 3 - 22, hi + 4 - 22, hi + 9 - 22]),
            vec![Cigar::Match(30), Cigar::Ins(6), Cigar::Del(28), Cigar::Match(40)],
            0,
        );
        let reads = [fragments(&both, 4, "a"), fragments(&plain, 4, "r")].concat();
        let g = guard_for(&reads, &v, &|_| true, 20, &[]).expect("an isolated change");
        assert_eq!(g.masked, vec![(hi + 9) as i64]);
        assert!(g.larger.is_none());
        // One REF read's error at a base of the run is no background (the REF
        // fragments must show it as often as a recurring change): still another
        // allele, not masked.
        let erred = aligned(&with_c(&ctx, &[hi + 3]), vec![Cigar::Match(98)], 0);
        let reads = [fragments(&run, 4, "a"), fragments(&plain, 5, "r"), fragments(&erred, 1, "e")].concat();
        assert!(guard_for(&reads, &v, &|_| true, 20, &[]).is_none(), "one REF read's error");
        // A change of the ALT reads' own beside a germline SNP (on the REF reads too)
        // is isolated among their own changes: both masked, the own one named.
        let beside = aligned(
            &with_c(&carrier_hap, &[hi + 3 - 22, hi + 4 - 22]),
            vec![Cigar::Match(30), Cigar::Ins(6), Cigar::Del(28), Cigar::Match(40)],
            0,
        );
        let germline = aligned(&with_c(&ctx, &[hi + 4]), vec![Cigar::Match(98)], 0);
        let reads = [fragments(&beside, 4, "a"), fragments(&germline, 4, "r")].concat();
        let g = guard_for(&reads, &v, &|_| true, 20, &[]).expect("an own change and a germline SNP");
        assert_eq!(g.masked, vec![(hi + 3) as i64, (hi + 4) as i64]);
        assert!(g.larger.is_some());
        // The same run on the REF reads too is background (a germline MNP): masked,
        // and not the ALT reads' own.
        let ref_run = aligned(&with_c(&ctx, &[hi + 3, hi + 4]), vec![Cigar::Match(98)], 0);
        let reads = [fragments(&run, 4, "a"), fragments(&ref_run, 4, "r")].concat();
        let g = guard_for(&reads, &v, &|_| true, 20, &[]).expect("a germline run");
        assert_eq!(g.masked, vec![(hi + 3) as i64, (hi + 4) as i64]);
        assert!(g.larger.is_none());
        for snp in [false, true] {
            let shared: &[usize] = if snp { &[hi + 6] } else { &[] };
            let ref_read = aligned(&with_c(&ctx, shared), vec![Cigar::Match(98)], 0);
            let reads = [fragments(&carrier, 4, "a"), fragments(&ref_read, 4, "r")].concat();
            let g = guard_for(&reads, &v, &|_| true, 20, &[]).expect("a recurring change");
            let l = g.larger.clone().expect("a larger allele");
            let (pos, r, a, n, m) = (l.pos, l.ref_allele, l.alt_allele, l.carriers, l.given_carriers);
            // left-aligned and minimal: the grown event's shared bases are trimmed
            // both changes recur beside the event: masked, a separate event
            assert_eq!(g.masked, vec![(hi + 3) as i64, (hi + 6) as i64]);
            if snp {
                // the shared C is not the ALT reads' alone: the larger allele stops
                // at the other one
                assert_eq!((pos, r.as_str(), n, m), (30, &ctx[30..hi + 4], 4, 0));
                assert_eq!(a, format!("TTCTCT{}C", &ctx[58..hi + 3]));
            } else {
                assert_eq!((pos, r.as_str(), n, m), (30, &ctx[30..hi + 7], 4, 0));
                assert_eq!(a, format!("TTCTCT{}C{}C", &ctx[58..hi + 3], &ctx[hi + 4..hi + 6]));
            }
        }
    }

    /// `n` bases of A and G, from a seeded generator.
    fn ag(n: usize, seed: u32) -> String {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                if x >> 16 & 1 == 0 { 'A' } else { 'G' }
            })
            .collect()
    }

    /// A read of `seq` aligned at 0 by `cigar`, at Q30 save Q5 at query positions `low`.
    fn aligned_low(seq: &[u8], cigar: Vec<Cigar>, low: &[usize]) -> Record {
        let mut quals = vec![30; seq.len()];
        low.iter().for_each(|&i| quals[i] = 5);
        let mut rec = Record::new();
        rec.set(b"r", Some(&CigarString(cigar)), seq, &quals);
        rec.set_pos(0);
        rec
    }

    fn spelled(o: &OtherAllele, ctx: &str) -> String {
        String::from_utf8(o.on(0, ctx.as_bytes()).expect("on the reference")).unwrap()
    }

    #[test]
    fn a_haplotype_needs_the_reference_to_hold_the_allele() {
        assert_eq!(haplotype(b"ACGTAC", 2, b"GT", b"T").as_deref(), Some(&b"ACTAC"[..]));
        assert_eq!(haplotype(b"ACGTAC", 2, b"GA", b"T"), None, "another REF");
        assert_eq!(haplotype(b"ACGTAC", 5, b"CA", b"T"), None, "past the end");
    }

    #[test]
    fn the_guard_finds_the_other_alleles_of_the_alts_length() {
        // A 28-base delins to TTCTCT in A/G sequence; TTCGCT is another allele.
        let ctx = ag(98, 7);
        let v = var(&ctx, 30, &ctx[30..58], "TTCTCT");
        let hap = |alt: &str| format!("{}{alt}{}", &ctx[..30], &ctx[58..]);
        let carrier =
            |alt: &str| aligned(hap(alt).as_bytes(), vec![Cigar::Match(30), Cigar::Ins(6), Cigar::Del(28), Cigar::Match(40)], 0);
        let plain = aligned(ctx.as_bytes(), vec![Cigar::Match(98)], 0);
        let base = [fragments(&carrier("TTCTCT"), 4, "a"), fragments(&plain, 4, "r")].concat();
        let others = |reads: &[Record], sibs: &[Variant]| {
            guard_for(reads, &v, &|_| true, 20, sibs).map_or_else(Vec::new, |g| g.others.iter().map(|o| spelled(o, &ctx)).collect())
        };
        let found = others(&[base.clone(), fragments(&carrier("TTCGCT"), 3, "o")].concat(), &[]);
        assert_eq!(found, vec![hap("TTCGCT")], "spelled in full on three fragments");
        let two = others(&[base.clone(), fragments(&carrier("TTCGCT"), 2, "o")].concat(), &[]);
        assert!(two.is_empty(), "two fragments are not looked for");
        let sib = var(&ctx, 30, &ctx[30..58], "TTCGCT");
        assert!(others(&base, std::slice::from_ref(&sib)).is_empty(), "co-annotated, never spelled: part of the event");
        let once = [base.clone(), fragments(&carrier("TTCGCT"), 1, "o")].concat();
        assert_eq!(others(&once, std::slice::from_ref(&sib)), vec![hap("TTCGCT")], "co-annotated, spelled once");
        let longer = var(&ctx, 30, &ctx[30..58], "TTCGCTC");
        assert!(others(&base, &[longer, v.clone()]).is_empty(), "another length, or the ALT itself");
    }

    #[test]
    fn an_allele_lacking_part_of_the_event_is_another_and_a_shared_change_is_not() {
        // GC>AT beside the germline allele GT: GT differs from the ALT only by the
        // reference's own base, so the REF molecules carrying it do not excuse it.
        let ctx = format!("{}CGCA{}", ag(40, 11), ag(40, 12));
        let p = 41;
        let v = var(&ctx, p as i64, "GC", "AT");
        let block = |b: &str| format!("{}{b}{}", &ctx[..p], &ctx[p + 2..]);
        let read = |b: &str| aligned(block(b).as_bytes(), vec![Cigar::Match(ctx.len() as u32)], 0);
        let reads = [fragments(&read("GC"), 4, "r"), fragments(&read("AT"), 4, "a"), fragments(&read("GT"), 4, "g")].concat();
        let g = guard_for(&reads, &v, &|_| true, 20, &[]).expect("the germline allele");
        assert_eq!(g.others.iter().map(|o| spelled(o, &ctx)).collect::<Vec<_>>(), vec![block("GT")]);
        // CCCAC>ACCAG with a germline C>T at its middle base on the REF and ALT reads
        // alike: background at a base the MNP keeps, masked (RJ-23, 2026-10-14), so
        // neither the ALT nor REF with it is another allele.
        let ctx = format!("{}GCCCACA{}", ag(40, 13), ag(40, 14));
        let v = var(&ctx, p as i64, "CCCAC", "ACCAG");
        let block = |b: &str| format!("{}{b}{}", &ctx[..p], &ctx[p + 5..]);
        let read = |b: &str| aligned(block(b).as_bytes(), vec![Cigar::Match(ctx.len() as u32)], 0);
        let reads = [
            fragments(&read("CCCAC"), 4, "r"),
            fragments(&read("CCTAC"), 4, "s"),
            fragments(&read("ACTAG"), 4, "a"),
            fragments(&read("ACCAG"), 4, "b"),
        ]
        .concat();
        let g = guard_for(&reads, &v, &|_| true, 20, &[]).expect("the background at the kept base");
        assert!(g.others.is_empty(), "no other allele");
        assert!(g.masked.contains(&((p + 2) as i64)), "the middle base masked");
    }

    #[test]
    fn a_carrier_aligned_with_its_gap_before_the_event_is_no_ref_molecule() {
        // C>TA after a TT: an aligner may write the carrier as a T inserted before the
        // run and C>A. Placed by its bases it is no REF molecule: the REF tally holds
        // no A where REF has the C.
        let left = ag(40, 21);
        let ctx = format!("{left}GTTCAAAAAAG{}", ag(40, 22));
        let p = left.len() + 3;
        let v = var(&ctx, p as i64, "C", "TA");
        let hap = format!("{}TA{}", &ctx[..p], &ctx[p + 1..]);
        let carrier = aligned(hap.as_bytes(), vec![Cigar::Match((p - 2) as u32), Cigar::Ins(1), Cigar::Match((ctx.len() - p + 2) as u32)], 0);
        let plain = aligned(ctx.as_bytes(), vec![Cigar::Match(ctx.len() as u32)], 0);
        let reads = [fragments(&carrier, 4, "a"), fragments(&plain, 4, "r")].concat();
        let (start, reference) = prepared_reference(&v).unwrap();
        let ev = event(start, &reference, &v).unwrap();
        let placed = placed_reads(&reads, &|_| true, 20, (start, &reference), &ev, super::super::window::read_window(&v))
            .expect("the windows fit");
        let at = p - (ev.lo - FLANK);
        assert_eq!(placed.ref_cover[at], 4, "the four REF reads");
        assert_eq!(placed.ref_bases[at][0], 0, "no carrier's A");
    }

    /// G then CACACA then T: a G>TCA grown through the repeat, with an A>G at the
    /// repeat's fourth base on the REF reads, the carriers, or both.
    fn ca_run(on_ref: bool, on_alt: bool) -> (String, usize, Option<Guard>) {
        let left = ag(40, 31);
        let ctx = format!("{left}CGCACACATC{}", ag(40, 32));
        let p = left.len() + 1;
        let v = var(&ctx, p as i64, "G", "TCA");
        let snp = |s: &str, at: usize| format!("{}G{}", &s[..at], &s[at + 1..]);
        let hap = format!("{}TCA{}", &ctx[..p], &ctx[p + 1..]);
        let carrier = |s: &str| aligned(s.as_bytes(), vec![Cigar::Match((p + 1) as u32), Cigar::Ins(2), Cigar::Match((ctx.len() - p - 1) as u32)], 0);
        let plain = |s: &str| aligned(s.as_bytes(), vec![Cigar::Match(ctx.len() as u32)], 0);
        let mut reads = [fragments(&plain(&ctx), 4, "r"), fragments(&carrier(&hap), 4, "a")].concat();
        reads.extend(fragments(&plain(&if on_ref { snp(&ctx, p + 4) } else { ctx.clone() }), 4, "s"));
        reads.extend(fragments(&carrier(&if on_alt { snp(&hap, p + 6) } else { hap.clone() }), 4, "b"));
        let g = guard_for(&reads, &v, &|_| true, 20, &[]);
        (ctx, p, g)
    }

    #[test]
    fn a_change_the_ref_molecules_carry_inside_the_grown_event_is_masked() {
        for (on_ref, on_alt, masked) in [(true, false, true), (true, true, true), (false, true, false)] {
            let (_, p, g) = ca_run(on_ref, on_alt);
            let got = g.is_some_and(|g| g.masked.contains(&((p + 4) as i64)));
            assert_eq!(got, masked, "on REF {on_ref}, on ALT {on_alt}");
        }
    }

    /// GCA>ACT (the C kept) with four fragments each of REF and the ALT, and four more of
    /// each with a T at the kept base (`on_ref`, `on_alt`) or without.
    fn kept_base(on_ref: bool, on_alt: bool) -> (usize, Option<Guard>) {
        let left = ag(40, 31);
        let ctx = format!("{left}GCA{}", ag(40, 32));
        let p = left.len();
        let v = var(&ctx, p as i64, "GCA", "ACT");
        let t = |s: &str| format!("{}T{}", &s[..p + 1], &s[p + 2..]);
        let hap = format!("{}ACT{}", &ctx[..p], &ctx[p + 3..]);
        let plain = |s: &str| aligned(s.as_bytes(), vec![Cigar::Match(ctx.len() as u32)], 0);
        let mut reads = [fragments(&plain(&ctx), 4, "r"), fragments(&plain(&hap), 4, "a")].concat();
        reads.extend(fragments(&plain(&if on_ref { t(&ctx) } else { ctx.clone() }), 4, "s"));
        reads.extend(fragments(&plain(&if on_alt { t(&hap) } else { hap.clone() }), 4, "b"));
        (p, guard_for(&reads, &v, &|_| true, 20, &[]))
    }

    #[test]
    fn a_change_the_ref_molecules_carry_at_a_base_an_mnp_keeps_is_masked() {
        for (on_ref, on_alt, masked) in [(true, false, true), (true, true, true), (false, true, false)] {
            let (p, g) = kept_base(on_ref, on_alt);
            let got = g.is_some_and(|g| {
                let at = (p + 1) as i64;
                g.masked.contains(&at) && g.background.contains(&(at, bit(b'C') | bit(b'T')))
            });
            assert_eq!(got, masked, "on REF {on_ref}, on ALT {on_alt}");
        }
    }

    #[test]
    fn a_masked_window_base_takes_only_the_bases_of_its_code() {
        // the C at the window's second base masked to C or T (Y); N takes any base
        let q = [30u8; 4];
        assert_eq!(score(b"ATGT", &q, 20, b"AYGT", (0, 4)).mismatches, 0, "the background base");
        assert_eq!(score(b"ACGT", &q, 20, b"AYGT", (0, 4)).mismatches, 0, "the reference's base");
        assert_eq!(score(b"AGGT", &q, 20, b"AYGT", (0, 4)).mismatches, 1, "a third base");
        assert_eq!(score(b"AGGT", &q, 20, b"ANGT", (0, 4)).mismatches, 0, "any base beside the event");
        assert!(confident_weight(b"AGGT", &q, 20, b"AYGT", 0, (0, 4)) > 0.0);
        // a read tells the allele from another only by a base fitting the one, not the other
        let aw = Window { allele: b"AYGT".to_vec(), hold: 4, tail: false, left: None, right: None, before: 0, after: 0, past: Vec::new(), outer: Vec::new(), event_at: (0, 0) };
        assert!(tells_others(b"ACGT", &q, 20, &aw, &[b"AGGT".to_vec()], 0, (0, 4)), "a third base is another allele's");
        assert!(!tells_others(b"ATGT", &q, 20, &aw, &[b"AYGT".to_vec()], 0, (0, 4)), "background tells nothing");
    }

    #[test]
    fn masks_inside_the_grown_event_sit_where_each_allele_has_the_base() {
        let left = ag(40, 31);
        let ctx = format!("{left}CGCACACATC{}", ag(40, 32));
        let p = left.len() + 1;
        let v = var(&ctx, p as i64, "G", "TCA");
        let (start, reference) = prepared_reference(&v).unwrap();
        let lo = event(start, &reference, &v).unwrap().lo - FLANK;
        let other = OtherAllele { pos: p as i64, ref_allele: b"G".to_vec(), alt_allele: b"ACA".to_vec() };
        let layout = Layout { masked: vec![(p + 4) as i64], others: vec![other], ..Layout::default() };
        let w = windows(&v, &layout).unwrap();
        let (r, a) = only_pair(&w);
        // right of the given allele: the REF's base at its place, the ALT's two later
        assert_eq!(r.allele[p + 4 - lo], WILD);
        assert_eq!(a.allele[p + 6 - lo], WILD);
        assert_eq!(a.allele.iter().filter(|&&b| b == WILD).count(), 1);
        // the other allele's window: the ALT's, save its own first base
        let o = &w.others[0];
        let differ: Vec<usize> = (0..o.len()).filter(|&i| o[i] != a.allele[i]).collect();
        assert_eq!(differ, vec![p - lo]);
        assert_eq!(o[p + 6 - lo], WILD, "the masks kept");
    }

    #[test]
    fn a_read_tells_the_alt_from_another_allele_only_by_a_base_it_reads() {
        let aw = Window { allele: b"AACGT".to_vec(), hold: 5, tail: false, left: None, right: None, before: 0, after: 0, past: Vec::new(), outer: Vec::new(), event_at: (0, 0) };
        let other = vec![b"AACTT".to_vec()];
        let (seq, good) = (b"AACGT", [30u8; 5]);
        assert!(tells_others(seq, &good, 20, &aw, &other, 0, (0, 5)), "the G read");
        let mut low = good;
        low[3] = 5;
        assert!(!tells_others(seq, &low, 20, &aw, &other, 0, (0, 5)), "the G below min BQ");
        assert!(!tells_others(seq, &good, 20, &aw, &other, 0, (0, 3)), "the G not readable");
        assert!(!tells_others(seq, &good, 20, &aw, &[b"AACNT".to_vec()], 0, (0, 5)), "a masked difference tells nothing");
        assert!(tells_others(seq, &low, 20, &aw, &[], 0, (0, 5)), "no other allele");
    }

    #[test]
    fn an_mnp_read_tells_the_alt_from_another_allele_by_a_base_it_reads() {
        let ctx = format!("{}CGCA{}", ag(40, 11), ag(40, 12));
        let p = 41;
        let v = var(&ctx, p as i64, "GC", "AT");
        let mut guard = Guard::new(&v, Vec::new(), None, 0);
        let read = |low: &[usize]| {
            let seq = format!("{}AT{}", &ctx[..p], &ctx[p + 2..]);
            aligned_low(seq.as_bytes(), vec![Cigar::Match(ctx.len() as u32)], low)
        };
        let quals = |r: &Record| r.qual().to_vec();
        let tells = |r: &Record, g: &Guard| tells_others_aligned(r, &quals(r), 20, &v, g, Side::Alt);
        assert!(tells(&read(&[p]), &guard), "no other allele");
        guard.others = vec![OtherAllele { pos: p as i64, ref_allele: b"GC".to_vec(), alt_allele: b"GT".to_vec() }];
        assert!(!tells(&read(&[p]), &guard), "?T: only the base GT shares");
        assert!(tells(&read(&[p + 1]), &guard), "A?: the A tells");
        assert!(tells(&read(&[]), &guard), "both read");
    }

    #[test]
    fn an_allele_the_reads_spell_must_recur_above_the_errors_of_the_one_it_is_nearest() {
        // GC>AT: GT is one base from REF and from the ALT, so it must reach 5% of the
        // larger (REF here); AA is one base from the ALT and two from REF, so 5% of the
        // ALT's fragments.
        let ctx = format!("{}CGCA{}", ag(40, 11), ag(40, 12));
        let p = 41;
        let v = var(&ctx, p as i64, "GC", "AT");
        let block = |b: &str| format!("{}{b}{}", &ctx[..p], &ctx[p + 2..]);
        let read = |b: &str| aligned(block(b).as_bytes(), vec![Cigar::Match(ctx.len() as u32)], 0);
        for (other, n_ref, found) in [("GT", 20, true), ("GT", 100, false), ("AA", 100, true)] {
            let reads = [fragments(&read("GC"), n_ref, "r"), fragments(&read("AT"), 4, "a"), fragments(&read(other), 3, "o")].concat();
            let others = guard_for(&reads, &v, &|_| true, 20, &[]).map_or(0, |g| g.others.len());
            assert_eq!(others == 1, found, "{other} beside {n_ref} REF reads");
        }
    }

    #[test]
    fn an_allele_differing_only_where_the_guard_masks_is_not_another() {
        // The CA repeat's germline SNP on half the REF reads is masked; a sibling merging
        // the event with it differs from the ALT only there.
        let left = ag(40, 31);
        let ctx = format!("{left}CGCACACATC{}", ag(40, 32));
        let p = left.len() + 1;
        let v = var(&ctx, p as i64, "G", "TCA");
        let hap = format!("{}TCA{}", &ctx[..p], &ctx[p + 1..]);
        let snp = format!("{}G{}", &ctx[..p + 4], &ctx[p + 5..]);
        let carrier = aligned(hap.as_bytes(), vec![Cigar::Match((p + 1) as u32), Cigar::Ins(2), Cigar::Match((ctx.len() - p - 1) as u32)], 0);
        let plain = |s: &str| aligned(s.as_bytes(), vec![Cigar::Match(ctx.len() as u32)], 0);
        let reads = [fragments(&plain(&ctx), 4, "r"), fragments(&plain(&snp), 4, "s"), fragments(&carrier, 8, "a")].concat();
        let merged = var(&ctx, p as i64, &ctx[p..p + 5], "TCACACG");
        let g = guard_for(&reads, &v, &|_| true, 20, std::slice::from_ref(&merged)).expect("the background");
        assert!(g.masked.contains(&((p + 4) as i64)));
        assert!(g.others.is_empty(), "only a masked difference");
        // a sibling with another base inside the event as well is another allele, once a
        // fragment spells it
        let other = var(&ctx, p as i64, &ctx[p..p + 5], "ACACACG");
        let spelled = format!("{}ACACACG{}", &ctx[..p], &ctx[p + 5..]);
        let one = aligned(spelled.as_bytes(), vec![Cigar::Match((p + 1) as u32), Cigar::Ins(2), Cigar::Match((ctx.len() - p - 1) as u32)], 0);
        let reads = [reads, fragments(&one, 1, "o")].concat();
        let g = guard_for(&reads, &v, &|_| true, 20, std::slice::from_ref(&other)).expect("the background");
        assert_eq!(g.others.len(), 1);
    }

    #[test]
    fn a_read_is_read_whole_over_the_whole_read_reach_not_the_prepared_reference() {
        // A 10-bp unit inserted in a two-unit repeat; prep's reference holds 60 bases each
        // side of the event, as production's does, the whole-read reach the contig.
        let (left, unit, right) = (ag(200, 51), "CTTGACTCGT", ag(200, 52));
        let ctx = format!("{left}{unit}{unit}{right}");
        let alt_hap = format!("{left}{unit}{unit}{unit}{right}");
        let mut v = var(&ctx, 219, "T", &format!("T{unit}"));
        v.event_ref = Some((140, ctx[140..280].to_owned()));
        let wide = Variant { whole_ref: Some((0, ctx.clone())), ..v.clone() };
        let ev = event(0, ctx.as_bytes(), &wide).unwrap();
        let alt_hi = ev.hi + 10 + ev.flank();
        let fit = |v: &Variant, seq: &[u8]| {
            whole_haps(v, &Layout::default()).and_then(|w| whole_fit(&w, seq, &vec![30; seq.len()], 20, (0, seq.len())))
        };
        // 150 bases from 120: the left flank, the three units and the right flank
        let alt_read = alt_hap.as_bytes()[120..270].to_vec();
        let ref_read = ctx.as_bytes()[120..270].to_vec();
        assert_eq!(fit(&v, &alt_read), None, "longer than the prepared reference: never read in part");
        assert_eq!(fit(&wide, &alt_read), Some(Some(true)));
        assert_eq!(fit(&wide, &ref_read), Some(Some(false)));
        // an error past the windows is tolerated, two are not, nor one on a window's base
        let with = |at: &[usize]| {
            let mut s = alt_read.clone();
            at.iter().for_each(|&i| s[i] = b'T');
            s
        };
        assert_eq!(fit(&wide, &with(&[10])), Some(Some(true)), "one error past the windows");
        assert_eq!(fit(&wide, &with(&[10, 140])), None, "two errors");
        assert_eq!(fit(&wide, &with(&[alt_hi - 1 - 120])), None, "the window's flank base misread");
    }

    #[test]
    fn reads_placed_apart_carry_a_change_beyond_an_error_or_leave_it_the_alts_own() {
        assert_eq!(carried_apart(14, 0, 0, 7), None, "too few read it");
        assert_eq!(carried_apart(14, 2, 1, 7), None, "too few, two errors");
        assert_eq!(carried_apart(15, 0, 0, 7), Some(false), "enough read it, none carries it");
        assert_eq!(carried_apart(73, 1, 0, 36), Some(false), "one is an error");
        assert_eq!(carried_apart(73, 2, 1, 36), Some(false), "two are errors");
        assert_eq!(carried_apart(10, 3, 1, 5), Some(true), "three carry it, however few read it");
        assert_eq!(carried_apart(73, 3, 1, 36), Some(true), "a haplotype on a twentieth of them");
        assert_eq!(carried_apart(247, 3, 1, 123), Some(false), "three errors at depth");
        assert_eq!(carried_apart(2500, 60, 30, 1250), Some(false), "a systematic error on 2.4% of them");
        assert_eq!(carried_apart(2500, 75, 37, 1250), Some(true), "a haplotype on 3% of them");
        assert_eq!(carried_apart(2500, 75, 75, 1250), Some(false), "all in one pair orientation: an artifact");
        assert_eq!(carried_apart(2500, 75, 15, 1250), Some(false), "a fifth one way: too unlikely for a haplotype");
        assert_eq!(carried_apart(2500, 75, 30, 1250), Some(true), "two fifths one way: a haplotype's chance");
        assert_eq!(carried_apart(78, 8, 8, 39), Some(true), "eight carriers all one way: a haplotype's chance (0.4%)");
        assert_eq!(carried_apart(200, 12, 12, 100), Some(false), "twelve all one way: an artifact (0.02%)");
        assert_eq!(carried_apart(2500, 75, 75, 2500), Some(true), "an amplicon library, read in one orientation: unweighed");
        assert_eq!(carried_apart(2500, 75, 75, 2400), Some(true), "under a fifth of all the reads the other way");
        assert!((binomial_cdf(0, 10, 0.5) - 0.000977).abs() < 1e-6);
        assert!((binomial_cdf(5, 10, 0.5) - 0.623047).abs() < 1e-6);
        assert_eq!(binomial_cdf(10, 10, 0.3), 1.0);
        assert!(binomial_cdf(3, 5000, 0.5) < 1e-300, "deep: no term overflows");
    }

    #[test]
    fn a_reads_ungapped_block_ends_at_its_indels_splices_and_clips() {
        let read = |cigar: Vec<Cigar>| aligned(&[b'A'; 40], cigar, 100);
        let ins = read(vec![Cigar::Match(10), Cigar::Ins(2), Cigar::Match(10), Cigar::SoftClip(18)]);
        assert_eq!(ungapped_block(&ins, 105), Some((100, 110, 0)));
        assert_eq!(ungapped_block(&ins, 112), Some((110, 120, 12)));
        assert_eq!(ungapped_block(&ins, 99), None);
        assert_eq!(ungapped_block(&ins, 120), None, "past the aligned bases");
        let del = read(vec![Cigar::SoftClip(5), Cigar::Match(10), Cigar::Del(2), Cigar::Match(25)]);
        assert_eq!(ungapped_block(&del, 105), Some((100, 110, 5)));
        assert_eq!(ungapped_block(&del, 110), None, "a deleted base");
        assert_eq!(ungapped_block(&del, 112), Some((112, 137, 15)));
        let spliced = read(vec![Cigar::Match(10), Cigar::RefSkip(5), Cigar::Match(30)]);
        assert_eq!(ungapped_block(&spliced, 116), Some((115, 145, 10)));
        let split = read(vec![Cigar::Match(10), Cigar::Equal(10), Cigar::Diff(20)]);
        assert_eq!(ungapped_block(&split, 125), Some((100, 140, 0)), "match ops of one kind or another, no gap");
    }

    #[test]
    fn the_whole_read_alleles_are_kept_per_variant() {
        let ctx = format!("{}CGCACACATC{}", ag(40, 31), ag(40, 32));
        let (v, sibling) = (var(&ctx, 41, "G", "TCA"), var(&ctx, 41, "G", "TA"));
        let cg = CarrierGuard::new(|| None);
        let a = cg.whole(&v, || whole_haps(&v, &Layout::default()));
        let b = cg.whole(&sibling, || whole_haps(&sibling, &Layout::default()));
        assert!(a.as_ref().as_ref().is_some_and(|x| x.key.2 == "TCA"));
        assert!(b.as_ref().as_ref().is_some_and(|x| x.key.2 == "TA"), "a sibling judged under the row's rules: its own");
        let again = cg.whole(&v, || panic!("built once per variant"));
        assert!(std::sync::Arc::ptr_eq(&a, &again));
    }

    #[test]
    fn a_read_read_whole_must_tell_the_alt_from_another_allele_whose_span_holds_a_mask() {
        // G>TCA with another allele of its length, GCACA>ACACACA, differing from it only
        // at the T; a carrier reading the T below the minimum quality is either. The guard
        // masks a change inside the other allele's REF span: it is still weighed.
        let left = ag(40, 31);
        let ctx = format!("{left}CGCACACATC{}", ag(40, 32));
        let p = left.len() + 1;
        let v = var(&ctx, p as i64, "G", "TCA");
        let other = OtherAllele { pos: p as i64, ref_allele: ctx.as_bytes()[p..p + 5].to_vec(), alt_allele: b"ACACACA".to_vec() };
        let hap = format!("{}TCA{}", &ctx[..p], &ctx[p + 1..]);
        let read = &hap.as_bytes()[5..hap.len() - 5];
        let fit = |masked: Vec<i64>, t_qual: u8| {
            let mut quals = vec![30u8; read.len()];
            quals[p - 5] = t_qual;
            let layout = Layout { masked, others: vec![other.clone()], ..Layout::default() };
            whole_haps(&v, &layout).and_then(|w| whole_fit(&w, read, &quals, 20, (0, read.len())))
        };
        assert_eq!(fit(vec![], 5), None, "the T unread: the ALT or the other allele");
        assert_eq!(fit(vec![(p + 4) as i64], 5), None, "a mask inside the other allele's REF span");
        assert_eq!(fit(vec![(p + 4) as i64], 30), Some(Some(true)), "the T read: the ALT");
    }
}
