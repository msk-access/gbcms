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
//!   through aligned bases counts only where the allele fits.
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
}

impl ReadRules<'static> {
    /// DNA, no guard: clipped bases are read; no read is spliced.
    pub(crate) const DNA: ReadRules<'static> = ReadRules { clip_edges: None, reference: None, spliced: None, guard: None };
}

/// The exon edges and junction ends near a variant, found on first use: only the
/// exact-carrier rule reads them, so an SNV or pure-indel row never pays for them.
pub(crate) struct ClipEdges<'a> {
    edges: std::cell::OnceCell<Vec<i64>>,
    find: Box<dyn Fn() -> Vec<i64> + 'a>,
}

impl<'a> ClipEdges<'a> {
    pub(crate) fn new(find: impl Fn() -> Vec<i64> + 'a) -> Self {
        Self { edges: std::cell::OnceCell::new(), find: Box::new(find) }
    }

    fn get(&self) -> &[i64] {
        self.edges.get_or_init(|| (self.find)())
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
/// Share of the REF reads on which a change is no longer the ALT's alone.
pub(crate) const GUARD_REF_SHARE: f64 = 0.05;
/// Longest tandem-repeat unit the event is grown through (besides the length change).
const MAX_UNIT: usize = 6;
/// Masked base: matches any haplotype base.
pub(crate) const WILD: u8 = b'N';
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
struct Windows {
    pair: (Window, Window),
    /// Where the alleles differ, genomic [start, end): the trims' union, before
    /// repeat growth. A read spliced through it shows neither allele; a read
    /// spliced beside it is judged on windows rebuilt around it over the reference
    /// spliced along the read's junctions.
    event: (i64, i64),
    /// Each other allele of the ALT's length at the locus (RJ-25), as the ALT's window
    /// over the same span: its own bases where it differs, the masks kept.
    others: Vec<Vec<u8>>,
}

impl Windows {
    /// The windows as a read spliced at `skips` sees them. A splice in a window's
    /// flank ends the window at the exon edge on that side: the read's next bases
    /// come from the next exon, and it is anchored at the junction. Both alleles'
    /// windows are cut alike, and the bases held shrink with them.
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
            (Window { allele, hold: w.hold, tail: w.tail, left, right, before, after }, dl, dr)
        };
        let ((mut r, dl, dr), (mut a, adl, adr)) = (cut(&self.pair.0), cut(&self.pair.1));
        let hold = r.hold.saturating_sub(if r.tail { dr } else { dl }).min(r.allele.len()).min(a.allele.len());
        r.hold = hold;
        a.hold = hold;
        // the other alleles share the ALT's span: cut as it is
        let others = self.others.iter().map(|w| w[adl..w.len() - adr].to_vec()).collect();
        Windows { pair: (r, a), event: self.event, others }
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
    /// The read's evidence against this allele (log10): the weight of every
    /// base that mismatches it, by its quality, masked bases included.
    pub(crate) mismatch_weight: f64,
    /// Query position of the window's first base (it may lie before the read).
    pub(crate) start: i64,
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
    let layout = Layout::for_variant(variant, rules);
    let win = windows(variant, &layout)?;
    let seq = record.seq().as_bytes();
    if seq.is_empty() || quals.len() < seq.len() {
        // A record stored without its bases (SEQ '*') shows nothing here.
        return Some(ClassifyResult::no_coverage(ClassifyPhase::MaskedCompare));
    }
    let skips = rna::splice_junctions_in(record, win.span());
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
        return Some(r);
    }
    // Such a read is also REF spliced at its own junction (an alternative donor or
    // acceptor, or a splice site the event straddles): it can count only ALT, and
    // only when its bases favour ALT over that reading too.
    let entering = skips.iter().any(enters);
    let spliced: Rc<SplicedWindows>;
    let cut;
    let identity = SpliceMap::default();
    let mut map = &identity;
    let win = if skips.is_empty() {
        &win
    } else if let Some(sw) = rules.reference.and_then(|r| {
        spliced_windows(variant, &rna::extract_splice_junctions(record), win.event, r, rules.spliced, &layout)
    }) {
        spliced = sw;
        map = &spliced.1;
        &spliced.0
    } else if entering {
        // No spliced windows to judge a junction entering the event on: depth only.
        trace!("carrier: spliced windows unavailable for a junction entering the event → depth only");
        let mut r = ClassifyResult::neither(ClassifyPhase::CigarRecon);
        r.uninformative = true;
        return Some(r);
    } else {
        // The far exon cannot be read, or the spliced reference no longer holds
        // the REF allele (its shared bases reach across the junction): judge the
        // read on the windows cut at its exon edges.
        trace!("carrier: spliced windows unavailable → windows cut at the exon edge");
        cut = win.cut_at(&skips);
        &cut
    };
    let Some(readable) = readable_range(record, rules) else {
        return Some(ClassifyResult::no_coverage(ClassifyPhase::MaskedCompare));
    };

    // The read's REF and ALT readings, when it holds the pair.
    let (rw, aw) = &win.pair;
    let placed = read_pair(record, &seq, quals, min_baseq, rw, aw, map, readable);
    let placement = placed.as_ref().map(|(_, _, p)| *p);
    let reading = placed.map(|(r, a, _)| (r, a));
    // (matches REF, matches ALT)
    let held = reading.as_ref().map(|(r, a)| (r.mismatches == 0, a.mismatches == 0));
    let (mm_ref, mm_alt) = reading.as_ref().map_or((0, 0), |(r, a)| (r.mismatches, a.mismatches));
    let alt_masked = reading.as_ref().map_or(0, |(_, a)| a.masked);
    let had_n = reading.as_ref().is_some_and(|(r, a)| r.had_n || a.had_n);
    // The read's evidence for ALT over REF (log10), every base weighed by its
    // quality: a low-quality base counts for little instead of fitting anything.
    // Bases matching both alleles cancel: the evidence is what mismatches REF less
    // what mismatches ALT, so the two readings need not cover the same bases.
    let evidence = reading.as_ref().map_or(0.0, |(r, a)| r.mismatch_weight - a.mismatch_weight);
    let alt_weight = reading.as_ref().map_or(0.0, |(_, a)| a.mismatch_weight);
    let alt_spans: Vec<(usize, usize)> = reading.iter().map(|(_, a)| a.span).collect();
    let read_spans: Vec<(usize, usize)> = reading.iter().flat_map(|(r, a)| [r.span, a.span]).collect();
    let qual_bases: Vec<u8> = read_spans.iter().flat_map(|&(lo, hi)| quals[lo..hi].iter().copied()).collect();
    let qual = median_qual(&qual_bases, min_baseq);
    let every = |want: (bool, bool)| held == Some(want);
    let mut result = if let Some(Placement::Unfit { alt }) = placement {
        // Placed as the allele its bases favour, which they do not fit: another
        // allele, partial evidence when that is the ALT.
        if alt {
            ClassifyResult::neither_with_nearby(qual, ClassifyPhase::MaskedCompare)
        } else {
            ClassifyResult::neither(ClassifyPhase::MaskedCompare)
        }
    } else if every((false, true)) && evidence < one_base_evidence(min_baseq) {
        // Its clearly read bases fit the ALT, but weighed by quality its bases do not
        // favour ALT over REF by as much as one base read at the minimum quality: a
        // REF molecule read poorly around the event, with one error fitting the ALT.
        trace!("carrier: fits the ALT but its evidence {:.2} is below one base's → neither", evidence);
        ClassifyResult::neither(ClassifyPhase::MaskedCompare)
    } else if every((false, true))
        && !reading.as_ref().is_some_and(|(_, a)| tells_others(&seq, quals, min_baseq, aw, &win.others, a.start, readable))
    {
        // Another allele of the ALT's length is at the locus and the read reads no
        // base telling the two apart (RJ-25): it is that allele or the ALT, partial
        // evidence.
        trace!("carrier: fits the ALT but reads no base telling it from another allele of its length → neither, partial");
        ClassifyResult::neither_with_nearby(qual, ClassifyPhase::MaskedCompare)
    } else if every((false, true)) {
        ClassifyResult::is_alt(qual, ClassifyPhase::MaskedCompare)
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
    Some(result)
}

/// Whether the read, with the ALT's window starting at query position `start`, reads
/// for every other allele's window a base telling the ALT from it: a position where
/// the two differ (neither masked), read inside `[first, after)` as the ALT's base at
/// or above `min_baseq`. True with no other allele.
fn tells_others(
    seq: &[u8],
    quals: &[u8],
    min_baseq: u8,
    aw: &Window,
    others: &[Vec<u8>],
    start: i64,
    (first, after): (usize, usize),
) -> bool {
    let reads_alt = |i: usize| {
        let q = start + i as i64;
        q >= first as i64
            && q < after as i64
            && seq.get(q as usize).is_some_and(|&b| b.to_ascii_uppercase() == aw.allele[i])
            && quals.get(q as usize).is_some_and(|&b| b >= min_baseq)
    };
    others.iter().all(|other| {
        other.iter().zip(&aw.allele).enumerate().any(|(i, (&o, &a))| o != a && o != WILD && a != WILD && reads_alt(i))
    })
}

/// Whether a read, as aligned, reads for every other allele of `guard` a base that
/// tells the ALT of `v` from it: a position where their haplotypes differ and the
/// guard masks nothing, read as the ALT's base at or above `min_baseq` (RJ-25). The
/// MNP check's rule for every read it calls ALT (it compares only the bases where REF
/// and ALT differ, so a read spelling another allele where the ALT keeps the
/// reference's base would pass it); the block is aligned without an indel, so the
/// haplotypes and the read share coordinates. True with no other allele it can place.
pub(crate) fn tells_others_aligned(record: &Record, quals: &[u8], min_baseq: u8, v: &Variant, guard: &Guard) -> bool {
    if guard.others.is_empty() {
        return true;
    }
    let Some((start, reference)) = prepared_reference(v) else { return true };
    let alt = usize::try_from(v.pos - start)
        .ok()
        .and_then(|p| haplotype(&reference, p, &upper(&v.ref_allele), &upper(&v.alt_allele)));
    let Some(alt) = alt else { return true };
    let seq = record.seq().as_bytes();
    let reads_alt = |i: usize| {
        find_read_pos(record, start + i as i64).is_some_and(|q| {
            seq.get(q).is_some_and(|&b| b.to_ascii_uppercase() == alt[i]) && quals.get(q).is_some_and(|&b| b >= min_baseq)
        })
    };
    guard.others.iter().all(|o| {
        let Some(h) = o.on(start, &reference).filter(|h| h.len() == alt.len()) else {
            trace!("carrier: another allele at {}:{} does not sit on the reference → not read", v.chrom, o.pos + 1);
            return true;
        };
        let (a, b) = o.span();
        let lo = usize::try_from(a.min(v.pos) - start).unwrap_or(0);
        let hi = usize::try_from(b.max(v.pos + v.ref_allele.len() as i64) - start).unwrap_or(0).min(alt.len());
        (lo..hi).any(|i| h[i] != alt[i] && !guard.masked.contains(&(start + i as i64)) && reads_alt(i))
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
        on_spliced.event_ref = Some((start, String::from_utf8(spliced).ok()?));
        on_spliced.ref_context = None;
        // masked positions are genome ones: each kept where a segment holds it
        let map = SpliceMap { segments };
        let masked = layout.masked.iter().filter_map(|&p| map.to_spliced(p)).collect();
        // other alleles are genome ones too: each kept where its whole span is where
        // it sits in the genome (no junction of the read crosses it)
        let same = |p: i64| map.to_spliced(p) == Some(p);
        let others = layout
            .others
            .iter()
            .filter(|o| {
                let (a, b) = o.span();
                same(a) && same(b - 1)
            })
            .cloned()
            .collect();
        let spliced_layout = Layout { masked, others, ..layout.clone() };
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
    matches!(
        super::window::allele_kind(&v.ref_allele, &v.alt_allele),
        Some(AlleleKind::Mnp | AlleleKind::Complex)
    )
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
        Some(_) => {
            "the prepared reference does not hold the event with its flank (a contig end, or a \
             repeat past the fetch cap)"
        }
    }
}

/// The reference the rule reads, upper case, and its genomic start: the event's
/// widened reference when prep fetched one, else the reference context.
fn prepared_reference(v: &Variant) -> Option<(i64, Vec<u8>)> {
    match (&v.event_ref, &v.ref_context) {
        (Some((s, seq)), _) => Some((*s, upper(seq))),
        (None, Some(ctx)) => Some((v.ref_context_start, upper(ctx))),
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

/// Another allele of the ALT's length at the variant's locus (RJ-25), on the
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
    /// ...and that grown through repeats touching it.
    lo: usize,
    hi: usize,
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
    let (b_lo, b_hi) = (l_lo.min(r_lo).min(p), l_hi.max(r_hi).max(p + r_al.len()));
    // Grow through repeats on either allele. Left of the event the alleles share
    // coordinates; right of it the haplotype's are shifted by d.
    let unit = MAX_UNIT.max(d.unsigned_abs() as usize);
    let lo = repeat_start(reference, b_lo, unit).min(repeat_start(&hap, b_lo, unit));
    let hap_hi = repeat_end(&hap, (b_hi as i64 + d) as usize, unit) as i64 - d;
    let hi = repeat_end(reference, b_hi, unit).max(hap_hi as usize);
    Some(Event { hap, d, first: l_lo, last: r_hi, core: (b_lo, b_hi), lo, hi })
}

/// How a variant's windows are laid out beyond the event itself: which flank they
/// are read from, the genome positions the neighbouring-change guard masks, and the
/// other alleles of the ALT's length at the locus.
#[derive(Clone, Default)]
struct Layout {
    from_right: bool,
    masked: Vec<i64>,
    others: Vec<OtherAllele>,
}

impl Layout {
    /// The layout for `v` under `rules`: its guard (when it was found for this very
    /// variant), and the flank it is read from (the left one, unless an RNA exon
    /// edge cuts it and not the right).
    fn for_variant(v: &Variant, rules: &ReadRules) -> Layout {
        let mut layout = rules
            .guard
            .and_then(|g| g.get(v))
            .map_or_else(Layout::default, |g| Layout { from_right: false, masked: g.masked.clone(), others: g.others.clone() });
        if let (Some(edges), Some((lo, hi))) = (rules.clip_edges.map(ClipEdges::get), grown_event(v)) {
            let cuts = |a: i64, b: i64| edges.iter().any(|&e| e > a && e <= b);
            let flank = FLANK as i64;
            if cuts(lo - flank, lo) && !cuts(hi, hi + flank) {
                trace!("carrier: an exon edge cuts the left flank of {}:{} → windows read from the right", v.chrom, v.pos + 1);
                layout.from_right = true;
            }
        }
        layout
    }
}

/// REF and ALT windows for the variant, laid out by `layout`. None without a
/// reference that holds the event and its flank, or with an empty allele.
fn windows(v: &Variant, layout: &Layout) -> Option<Windows> {
    let (start, reference) = prepared_reference(v)?;
    let ev = event(start, &reference, v)?;
    if ev.lo < FLANK || ev.hi + FLANK > reference.len() {
        return None; // the reference does not hold the event's flank
    }
    let d = ev.d;
    let lo = ev.lo - FLANK;
    let hi = ev.hi + FLANK;
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
    // Each other allele of the ALT's length, on this reference (one it does not hold,
    // or of another length, is not read here).
    let other_haps: Vec<Vec<u8>> =
        layout.others.iter().filter_map(|o| o.on(start, &reference)).filter(|h| h.len() == ev.hap.len()).collect();
    // A masked base (a change recurring beside the event, or one the REF molecules
    // carry inside the repeat it was grown through) matches any read base. Left of
    // the given allele the alleles share coordinates; right of it the ALT's are
    // shifted by the length change. The given allele's own bases are never masked.
    let (mut reference, mut hap) = (reference, ev.hap.clone());
    for &p in &layout.masked {
        let off = p - start;
        if off >= lo as i64 && off < c_lo as i64 {
            reference[off as usize] = WILD;
            hap[off as usize] = WILD;
        } else if off >= c_hi as i64 && off < hi as i64 {
            reference[off as usize] = WILD;
            hap[(off + d) as usize] = WILD;
        }
    }
    let (ref_win, alt_win) = (&reference[lo..hi], &hap[lo..alt_hi]);
    let others = other_haps
        .iter()
        .map(|h| h[lo..alt_hi].iter().zip(alt_win).map(|(&o, &a)| if a == WILD { WILD } else { o }).collect())
        .collect();
    // Both windows are anchored at the window's ends: the left flank's first base
    // and one past the right flank's last, shared by the two alleles.
    let window = |allele: &[u8]| Window {
        allele: allele.to_vec(),
        hold: j,
        tail,
        left: Some(g(lo)),
        right: Some(g(hi)),
        before: c_lo - lo,
        after: hi - c_hi,
    };
    let pair = (window(ref_win), window(alt_win));
    Some(Windows { pair, event: (g(c_lo), g(c_hi)), others })
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
    let ext = (first_diff + FLANK).saturating_sub(short);
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
    /// The reads' larger allele (0-based POS, left-aligned minimal VCF form): its
    /// carriers are the fewest ALT reads showing any one change, its given carriers
    /// the most showing the reference base at one.
    pub larger: Option<ObservedAllele>,
    /// Fragments read: those showing one allele's core exactly between the flank
    /// bases.
    pub scanned: u32,
    /// The other alleles of the ALT's length at the locus (RJ-25).
    pub others: Vec<OtherAllele>,
}

impl Guard {
    /// The guard found for `v`.
    pub(crate) fn new(v: &Variant, masked: Vec<i64>, larger: Option<ObservedAllele>, scanned: u32) -> Self {
        Self { key: (v.pos, v.ref_allele.clone(), v.alt_allele.clone()), masked, larger, scanned, others: Vec::new() }
    }
}

/// A variant's [`Guard`], found from its reads on first use: only the exact-carrier
/// rule and the junction judge read it, so an SNV or a short pure indel never pays
/// for it.
pub(crate) struct CarrierGuard<'a> {
    found: std::cell::OnceCell<Option<Guard>>,
    find: Box<dyn Fn() -> Option<Guard> + 'a>,
}

impl<'a> CarrierGuard<'a> {
    pub(crate) fn new(find: impl Fn() -> Option<Guard> + 'a) -> Self {
        Self { found: std::cell::OnceCell::new(), find: Box::new(find) }
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
/// the spellings of another allele of the ALT's length, and the bases of the REF and
/// the ALT molecules over the event and [`FLANK`] bases each side.
struct Placed<'a> {
    /// Spellings of the ALT haplotype's `[lo, alt_hi)` (the event with FLANK bases
    /// each side) that differ from the ALT inside the event, every base read, by the
    /// fragments showing each.
    spellings: HashMap<Vec<u8>, std::collections::HashSet<&'a [u8]>>,
    /// Per reference offset in `[lo, hi)`: REF molecules reading each base (A C G T),
    /// and reading the base at all; the same for ALT molecules per haplotype offset in
    /// `[lo, alt_hi)`.
    ref_bases: Vec<[u32; 4]>,
    ref_cover: Vec<u32>,
    alt_bases: Vec<[u32; 4]>,
    alt_cover: Vec<u32>,
    /// Every spelling of `[lo, alt_hi)` read in full, by the fragments showing each:
    /// a spelling must recur above the errors of the allele it is nearest, counted
    /// from these once the masks are known.
    full: HashMap<Vec<u8>, std::collections::HashSet<&'a [u8]>>,
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
        spellings: HashMap::new(),
        ref_bases: vec![[0; 4]; hi - lo],
        ref_cover: vec![0; hi - lo],
        alt_bases: vec![[0; 4]; alt_hi - lo],
        alt_cover: vec![0; alt_hi - lo],
        full: HashMap::new(),
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
        // Another allele of the ALT's length: every base of the window is read; for an
        // event changing the length, past it the read fits the ALT placed this way at
        // all but PLACED_CHANGES bases at most, and no worse than REF placed the same
        // way (an event keeping the length places both alike); and the bases differ
        // from the ALT's inside the event (beside it too, maybe: the ALT's own trims may
        // leave a base the other allele changes in the flank) and are not REF's own.
        let n = alt_hi - lo;
        let bases: Option<Vec<u8>> = (0..n).map(|i| read(i).filter(|&(b, ok)| ok && b != WILD).map(|(b, _)| b)).collect();
        let Some(bases) = bases else { continue };
        placed.full.entry(bases.clone()).or_default().insert(qname);
        let refs = reference.get(lo..lo + n);
        let past_fits = ev.d == 0 || {
            let past = alt_hi..alt_hi + PLACED_REACH;
            let misfits = |h: &[u8]| {
                past.clone().filter(|&i| read(i - lo).is_some_and(|(b, ok)| ok && b != WILD && h.get(i) != Some(&b))).count()
            };
            past.end <= ev.hap.len()
                && past.clone().all(|i| read(i - lo).is_some())
                && misfits(&ev.hap) <= PLACED_CHANGES
                && misfits(&ev.hap) <= misfits(reference)
        };
        let other = (FLANK..n - FLANK).any(|i| bases[i] != ev.hap[lo + i]) && refs != Some(&bases[..]);
        if past_fits && other {
            placed.spellings.entry(bases).or_default().insert(qname);
        }
    }
    Some(placed)
}

/// A share, 0 of nothing.
fn share(n: u32, of: u32) -> f64 {
    if of == 0 { 0.0 } else { n as f64 / of as f64 }
}

/// The reference offsets inside the repeat the event was grown through, off the
/// given allele's own bases, masked as beside the event (RJ-23, 2026-10-14). A change
/// the REF molecules carry is background: recurring on them, or recurring on the ALT
/// molecules and on [`GUARD_MIN_READS`] and [`GUARD_REF_SHARE`] of them. A change the
/// ALT molecules alone carry there is none for an event changing the length (inside
/// the repeat it may be the event's own sequence); for an MNP, which has nothing to
/// place, the run is beside it, so an isolated one is its phased neighbour, masked as
/// beside any event (a run of them is another allele). Read from the molecules placed
/// by their bases, so a carrier however aligned feeds it.
fn grown_background(placed: &Placed, reference: &[u8], ev: &Event) -> Vec<usize> {
    let lo = ev.lo - FLANK;
    let (c_lo, c_hi) = ev.core;
    let recurs = |n: u32, of: u32| n >= GUARD_MIN_READS && share(n, of) >= GUARD_MIN_SHARE;
    let grown = |off: usize| (ev.lo..c_lo).contains(&off) || (c_hi..ev.hi).contains(&off);
    // the ALT's haplotype offset: left of the given allele the same, right of it
    // shifted by the length change
    let on_alt = |off: usize| if off < c_lo { off } else { (off as i64 + ev.d) as usize };
    let background = |off: usize| {
        let h = on_alt(off);
        (0..4).any(|k| {
            let (n, of) = (placed.ref_bases[off - lo][k], placed.ref_cover[off - lo]);
            let (an, aof) = (placed.alt_bases[h - lo][k], placed.alt_cover[h - lo]);
            b"ACGT"[k] != reference[off]
                && (recurs(n, of) || (recurs(an, aof) && n >= GUARD_MIN_READS && share(n, of) >= GUARD_REF_SHARE))
        })
    };
    let own = |off: usize| {
        ev.d == 0
            && grown(off)
            && !background(off)
            && (0..4).any(|k| b"ACGT"[k] != reference[off] && recurs(placed.alt_bases[off - lo][k], placed.alt_cover[off - lo]))
    };
    let isolated = |off: usize| !(own(off + 1) || off > 0 && own(off - 1));
    (ev.lo..c_lo).chain(c_hi..ev.hi).filter(|&off| background(off) || (own(off) && isolated(off))).collect()
}

/// The other alleles of the ALT's length at the locus (RJ-25): the spellings read in
/// full on [`GUARD_MIN_READS`] fragments or more and on [`GUARD_REF_SHARE`] of the
/// fragments spelling the allele they are nearest, REF's own bases or the ALT's (the
/// larger, when as near): an allele recurs above the errors of the molecules it could
/// be an error of, not as a recurring error does. Save one whose every difference
/// from the ALT is a change from the reference the REF molecules carry too (on
/// GUARD_MIN_READS and GUARD_REF_SHARE of those reading it, where the ALT's base sits
/// on the reference: left-aligned, or right-aligned past the length change): a
/// separate event on both haplotypes. The reference's own base is no such change: an
/// allele lacking part of the event is another allele. And the co-annotated siblings
/// that spell one: the ALT's haplotype save inside the event and [`FLANK`] bases each
/// side, differing from it inside the event. None differs from the ALT only at
/// `masked` positions (genome): that is the ALT with a separate event, which no read
/// could tell from it. For an event keeping the length (an MNP) none differs from the
/// ALT only beside the given allele's own bases: with nothing to place, that is the
/// ALT with a phased neighbour (named, as beside any event), not another allele.
fn other_alleles(
    placed: &Placed,
    (start, reference): (i64, &[u8]),
    ev: &Event,
    siblings: &[Variant],
    masked: &[i64],
) -> Vec<OtherAllele> {
    let lo = ev.lo - FLANK;
    let (hi, alt_hi) = (ev.hi + FLANK, (ev.hi as i64 + ev.d) as usize + FLANK);
    // the masked positions on the ALT's haplotype: left of the given allele the same,
    // right of it shifted by the length change
    let (c_lo, c_hi) = ev.core;
    let masked_hap: Vec<usize> = masked
        .iter()
        .filter_map(|&p| usize::try_from(p - start).ok())
        .filter_map(|off| match off {
            off if off < c_lo => Some(off),
            off if off >= c_hi => usize::try_from(off as i64 + ev.d).ok(),
            _ => None,
        })
        .collect();
    let tellable = |h: &[u8]| {
        h.iter().zip(&ev.hap).enumerate().any(|(i, (x, y))| x != y && !masked_hap.contains(&i))
            && (ev.d != 0 || (c_lo..c_hi).any(|i| h[i] != ev.hap[i]))
    };
    let (ref_spelling, alt_spelling) = (reference.get(lo..alt_hi), &ev.hap[lo..alt_hi]);
    // How far a spelling is from REF's or the ALT's, at the bases that can tell: not
    // where the guard masks (on that allele's coordinates), nor, for an MNP, beside its
    // own bases (a phased neighbour there is no difference).
    let masked_ref: Vec<usize> = masked.iter().filter_map(|&p| usize::try_from(p - start).ok()).collect();
    let apart = |b: &[u8], allele: &[u8], masked_at: &[usize]| {
        b.iter()
            .zip(allele)
            .enumerate()
            .filter(|&(k, (x, y))| {
                let i = lo + k;
                x != y && !masked_at.contains(&i) && (ev.d != 0 || (c_lo..c_hi).contains(&i))
            })
            .count()
    };
    // the fragments spelling REF's own bases and the ALT's, so measured
    let spellers = |allele: Option<&[u8]>, masked_at: &[usize]| -> u32 {
        let Some(allele) = allele else { return 0 };
        let mut frags = std::collections::HashSet::new();
        for (b, f) in &placed.full {
            if apart(b, allele, masked_at) == 0 {
                frags.extend(f.iter().copied());
            }
        }
        frags.len() as u32
    };
    let (ref_n, alt_n) = (spellers(ref_spelling, &masked_ref), spellers(Some(alt_spelling), &masked_hap));
    // the fragments spelling the allele a spelling is nearest
    let source = |b: &[u8]| {
        let to_alt = apart(b, alt_spelling, &masked_hap);
        match ref_spelling.map(|r| apart(b, r, &masked_ref)) {
            Some(to_ref) if to_ref < to_alt => ref_n,
            Some(to_ref) if to_ref == to_alt => ref_n.max(alt_n),
            _ => alt_n,
        }
    };
    let ref_carries = |i: usize, b: u8| {
        let left = (i < ev.hi).then_some(i);
        let right = usize::try_from(i as i64 - ev.d).ok().filter(|&j| i >= ev.lo && j >= ev.lo);
        let k = b"ACGT".iter().position(|&x| x == b);
        [left, right].into_iter().flatten().any(|j| {
            k.is_some_and(|k| {
                let (n, of) = (placed.ref_bases[j - lo][k], placed.ref_cover[j - lo]);
                b != reference[j] && n >= GUARD_MIN_READS && share(n, of) >= GUARD_REF_SHARE
            })
        })
    };
    let mut others: Vec<OtherAllele> = placed
        .spellings
        .iter()
        .filter(|(b, f)| {
            let n = f.len() as u32;
            n >= GUARD_MIN_READS
                && n as f64 >= GUARD_REF_SHARE * source(b) as f64
                && b.iter().enumerate().any(|(k, &x)| x != ev.hap[lo + k] && !ref_carries(lo + k, x))
        })
        .map(|(b, _)| OtherAllele { pos: start + lo as i64, ref_allele: reference[lo..hi].to_vec(), alt_allele: b.clone() })
        .collect();
    others.retain(|o| o.on(start, reference).is_some_and(|h| tellable(&h)));
    let mut haps: Vec<Vec<u8>> = others.iter().filter_map(|o| o.on(start, reference)).collect();
    for sib in siblings {
        let o = OtherAllele { pos: sib.pos, ref_allele: upper(&sib.ref_allele), alt_allele: upper(&sib.alt_allele) };
        let Some(h) = o.on(start, reference) else { continue };
        let spells = h.len() == ev.hap.len()
            && h[..lo] == ev.hap[..lo]
            && h[alt_hi..] == ev.hap[alt_hi..]
            && (ev.lo..alt_hi - FLANK).any(|i| h[i] != ev.hap[i]);
        if spells && tellable(&h) && !haps.contains(&h) {
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
/// indel or MNP beside the event, another allele: not masked. None when nothing is masked or
/// named, or the variant is not judged by this rule.
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
    let (start, reference) = prepared_reference(v)?;
    let ev = event(start, &reference, v)?;
    // The scan beside the event starts where the event ends: for one changing the
    // length, past the repeat it was grown through (the change may sit anywhere in
    // it); for an MNP, at its own bases (it has nothing to place, so the run beside it
    // is beside it, and a change there a phased neighbour or background).
    let (s_lo, s_hi) = if ev.d == 0 { ev.core } else { (ev.lo, ev.hi) };
    if s_lo == 0 || s_hi >= reference.len() || ev.lo == 0 {
        return None;
    }
    let reach_l = GUARD_REACH.min(s_lo);
    let reach_r = GUARD_REACH.min(reference.len() - s_hi);
    let alt_core = &ev.hap[s_lo..(s_hi as i64 + ev.d) as usize];
    let ref_core = &reference[s_lo..s_hi];
    let g = |off: usize| start + off as i64;
    // The variant's own reads placed by their bases (RJ-25): the background inside the
    // repeat the event was grown through, which the scan beside the event lets fit
    // too, so a read carrying it is still read on past the event.
    let placed = placed_reads(read_cache, counted, min_baseq, (start, &reference), &ev, super::window::read_window(v));
    let grown = placed.as_ref().map_or_else(Vec::new, |p| grown_background(p, &reference, &ev));
    let c_lo = ev.core.0;
    let on_alt = |off: usize| if off < c_lo { off } else { (off as i64 + ev.d) as usize };
    let ref_wild: Vec<bool> = (s_lo..s_hi).map(|off| grown.contains(&off)).collect();
    let alt_wild: Vec<bool> = (s_lo..s_lo + alt_core.len()).map(|h| grown.iter().any(|&off| on_alt(off) == h)).collect();
    // Per allele (REF 0, ALT 1) and side (left 0, right 1), by distance from the
    // event (index 0 is the flank base beside it): reads reading the base, and reads
    // showing a confident change there.
    let mut cover = vec![vec![vec![0u32; GUARD_REACH]; 2]; 2];
    let mut scanned = 0u32;
    let mut change = cover.clone();
    // the bases the ALT reads show where they change, per side and distance (A C G T)
    let mut alt_bases = vec![vec![[0u32; 4]; GUARD_REACH]; 2];
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
        let (Some(ql), Some(qr)) = (find_read_pos(record, g(s_lo - 1)), find_read_pos(record, g(s_hi))) else {
            continue;
        };
        if qr <= ql {
            continue;
        }
        let between = &seq[ql + 1..qr];
        let fits = |core: &[u8], wild: &[bool]| {
            core.len() == between.len()
                && between.iter().zip(core).enumerate().all(|(i, (&b, &c))| {
                    let b = b.to_ascii_uppercase();
                    b == c || b == WILD || quals[ql + 1 + i] < min_baseq || wild[i]
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
                let off = if side == 0 { s_lo - 1 - i } else { s_hi + i };
                let Some(q) = find_read_pos(record, g(off)) else { break };
                cover[allele][side][i] += 1;
                let b = seq[q].to_ascii_uppercase();
                if b != WILD && quals[q] >= min_baseq && b != reference[off] {
                    change[allele][side][i] += 1;
                    if let (1, Some(k)) = (allele, b"ACGT".iter().position(|&x| x == b)) {
                        alt_bases[side][i][k] += 1;
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
    let (mut masked, mut alt_only, mut run) = (Vec::new(), [Vec::new(), Vec::new()], false);
    for (side, reach) in [(0usize, reach_l), (1usize, reach_r)] {
        for i in 0..reach {
            let off = if side == 0 { s_lo - 1 - i } else { s_hi + i };
            if recurs(1, side, i) && !on_ref(side, i) {
                alt_only[side].push(i);
            }
            if background(side, i) || (own(side, i) && isolated(side, i)) {
                masked.push(g(off));
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
        let lo = far(0).map_or(s_lo, |i| s_lo - 1 - i);
        let hi = far(1).map_or(s_hi, |i| s_hi + i + 1);
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
        alt.extend((lo..s_lo).map(|off| base_at(0, s_lo - 1 - off, off)));
        alt.extend_from_slice(alt_core);
        alt.extend((s_hi..hi).map(|off| base_at(1, off - s_hi, off)));
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
    // The background inside the grown event joins the masks; the other alleles of the
    // ALT's length are those a read can tell from it past every mask.
    for p in grown.iter().map(|&off| g(off)) {
        if !masked.contains(&p) {
            masked.push(p);
        }
    }
    let others =
        placed.as_ref().map_or_else(Vec::new, |p| other_alleles(p, (start, &reference), &ev, siblings, &masked));
    if masked.is_empty() && larger.is_none() && others.is_empty() {
        return None;
    }
    trace!(
        "carrier guard at {}:{} {}>{}: {} recurring changes beside the event masked ({} inside the repeat it grew \
         through), {} on the ALT reads alone, {} other alleles of its length",
        v.chrom, v.pos + 1, v.ref_allele, v.alt_allele, masked.len(), grown.len(),
        alt_only.iter().map(Vec::len).sum::<usize>(), others.len(),
    );
    let mut guard = Guard::new(v, masked, larger, scanned);
    guard.others = others;
    Some(guard)
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
) -> Option<(Reading, Reading, Placement)> {
    let at_left = rw.left.and_then(|g| find_read_pos(record, map.to_genome(g))).map(|q| q as i64);
    let at_right = rw.right.and_then(|g| find_read_pos(record, map.to_genome(g - 1))).map(|q| q as i64 + 1);
    let tail = rw.reads_from_right();
    let read = |w: &Window, s: i64| read_at(seq, quals, min_baseq, w, s, tail, readable);
    // Both windows from one placement of the reading flank's end base: their shared
    // start (left), or their shared end (right; each starts its own length before it).
    let from = |anchor: i64| {
        let start = |w: &Window| if tail { anchor - w.allele.len() as i64 } else { anchor };
        Some((read(rw, start(rw))?, read(aw, start(aw))?))
    };
    match (tail, at_left, at_right) {
        (false, Some(s), _) => return from(s).map(|(r, a)| (r, a, Placement::Flank)),
        (true, _, Some(e)) => return from(e).map(|(r, a)| (r, a, Placement::Flank)),
        _ => {}
    }
    // Placed from the other flank: each allele's window starts there (right reading)
    // or ends there (left reading), at its own length.
    let anchor = if tail { at_left? } else { at_right? };
    let start = |w: &Window| if tail { anchor } else { anchor - w.allele.len() as i64 };
    let weight = |w: &Window| confident_weight(seq, quals, min_baseq, &w.allele, start(w), readable);
    let evidence = weight(rw) - weight(aw);
    // At least one base read at min BQ, and more than none (at --min-baseq 0 that
    // base weighs nothing): a tie favours neither.
    let one = one_base_evidence(min_baseq);
    let (w, alt) = if evidence > 0.0 && evidence >= one {
        (aw, true)
    } else if evidence < 0.0 && evidence <= -one {
        (rw, false)
    } else {
        trace!("carrier: placed from the other flank, its confident bases favour neither allele ({:.2}) → depth only", evidence);
        return None;
    };
    // The reading flank's end base where the favoured allele puts it.
    let s = if tail { anchor + w.allele.len() as i64 } else { start(w) };
    let (r, a) = from(s)?;
    let (first, after) = aligned_query_range(record)?;
    let flank_end = if tail { rw.right.map(|g| g - 1) } else { rw.left };
    let deleted = flank_end.is_some_and(|g| deleted_anchor_query(record, map.to_genome(g)).is_some());
    // Held: the flank's end base lies in the clip on the reading side, or under the
    // read's own deletion; else it is aligned elsewhere.
    let held = deleted || if tail { s > after as i64 } else { s < first as i64 };
    let fits = if alt { a.mismatches == 0 } else { r.mismatches == 0 };
    let placement = match (fits, held) {
        (true, true) => Placement::OtherFlank,
        (true, false) => Placement::ThroughAligned,
        (false, true) => Placement::Unfit { alt },
        (false, false) => return None,
    };
    Some((r, a, placement))
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
            (b != WILD && h != WILD && bq >= min_baseq && b != h).then(|| one_base_evidence(bq))
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
    }
    Some(r)
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


/// Mismatches of read bases against haplotype bases position by position, masked
/// bases (below `min_baseq`, or N) matching anything.
pub(crate) fn score(bases: &[u8], quals: &[u8], min_baseq: u8, hap: &[u8], span: (usize, usize)) -> Reading {
    let (mut mismatches, mut masked, mut had_n, mut mismatch_weight) = (0, 0, false, 0.0);
    for ((&b, &q), &h) in bases.iter().zip(quals).zip(hap) {
        let b = b.to_ascii_uppercase();
        had_n |= b == WILD;
        if b != WILD && h != WILD && b != h {
            mismatch_weight += one_base_evidence(q);
        }
        if b == WILD || q < min_baseq {
            masked += 1;
        } else if b != h && h != WILD {
            mismatches += 1;
        }
    }
    Reading { mismatches, masked, had_n, span, mismatch_weight, start: span.0 as i64 }
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
        let found = edges.map(|e| ClipEdges::new(move || e.clone()));
        readable_range(rec, &ReadRules { clip_edges: found.as_ref(), reference: None, spliced: None, guard: None })
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
        assert_eq!(others(&base, std::slice::from_ref(&sib)), vec![hap("TTCGCT")], "co-annotated, with no read");
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
        // alike: the ALT with it is no other allele (a separate event on both
        // haplotypes); the REF with it is one (it lacks the event).
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
        let g = guard_for(&reads, &v, &|_| true, 20, &[]).expect("the REF haplotype with the SNP");
        assert_eq!(g.others.iter().map(|o| spelled(o, &ctx)).collect::<Vec<_>>(), vec![block("CCTAC")]);
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

    #[test]
    fn masks_inside_the_grown_event_sit_where_each_allele_has_the_base() {
        let left = ag(40, 31);
        let ctx = format!("{left}CGCACACATC{}", ag(40, 32));
        let p = left.len() + 1;
        let v = var(&ctx, p as i64, "G", "TCA");
        let (start, reference) = prepared_reference(&v).unwrap();
        let lo = event(start, &reference, &v).unwrap().lo - FLANK;
        let other = OtherAllele { pos: p as i64, ref_allele: b"G".to_vec(), alt_allele: b"ACA".to_vec() };
        let layout = Layout { from_right: false, masked: vec![(p + 4) as i64], others: vec![other] };
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
        let aw = Window { allele: b"AACGT".to_vec(), hold: 5, tail: false, left: None, right: None, before: 0, after: 0 };
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
        let tells = |r: &Record, g: &Guard| tells_others_aligned(r, &quals(r), 20, &v, g);
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
        // a sibling with another base inside the event as well is another allele
        let other = var(&ctx, p as i64, &ctx[p..p + 5], "ACACACG");
        let g = guard_for(&reads, &v, &|_| true, 20, std::slice::from_ref(&other)).expect("the background");
        assert_eq!(g.others.len(), 1);
    }
}
