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
//!   both ends). When both alleles are longer than [`LONG_EVENT`] bases the windows
//!   end one base past the first difference instead. The flank is the left one; an
//!   RNA variant whose left flank an exon edge cuts reads from the right.
//! - **A neighbouring change on the ALT reads.** A change that recurs at a fixed
//!   distance from the event on the reads carrying the ALT core, and not on the REF
//!   reads, makes their haplotype larger than the given allele: both windows grow
//!   past it ([`guard_for`]), so a read counts ALT only when it shows the reference
//!   base there. One that recurs on REF reads too (a germline SNP) is masked in both
//!   windows instead.
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

use super::observed::canonical;
use super::rna;
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
fn readable_range(record: &Record, rules: &ReadRules) -> Option<(usize, usize)> {
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
/// Reads, and share of the reads reading that base, a neighbouring change must
/// recur on before it is taken as part of the reads' haplotype.
const GUARD_MIN_READS: u32 = 3;
const GUARD_MIN_SHARE: f64 = 0.2;
/// Share of the REF reads on which a change is no longer the ALT's alone.
const GUARD_REF_SHARE: f64 = 0.05;
/// Longest tandem-repeat unit the event is grown through (besides the length change).
const MAX_UNIT: usize = 6;
/// Masked base: matches any haplotype base.
const WILD: u8 = b'N';

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

/// The REF and ALT windows, paired.
struct Windows {
    pairs: Vec<(Window, Window)>,
    /// Where the alleles differ, genomic [start, end): the trims' union, before
    /// repeat growth. A read spliced through it shows neither allele; a read
    /// spliced beside it is judged on windows rebuilt around it over the reference
    /// spliced along the read's junctions.
    event: (i64, i64),
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
        let pairs = self
            .pairs
            .iter()
            .map(|(r, a)| {
                let ((mut r, dl, dr), (mut a, _, _)) = (cut(r), cut(a));
                let hold = r.hold.saturating_sub(if r.tail { dr } else { dl }).min(r.allele.len()).min(a.allele.len());
                r.hold = hold;
                a.hold = hold;
                (r, a)
            })
            .collect();
        Windows { pairs, event: self.event }
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
#[derive(Clone)]
struct Reading {
    mismatches: usize,
    masked: usize,
    had_n: bool,
    span: (usize, usize),
    /// The read's evidence against this allele (log10): the weight of every
    /// base that mismatches it, by its quality, masked bases included.
    mismatch_weight: f64,
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

    // Per pair the read holds: (matches REF, matches ALT).
    let mut held: Vec<(bool, bool)> = Vec::with_capacity(win.pairs.len());
    let (mut mm_ref, mut mm_alt, mut alt_masked, mut had_n) = (0usize, 0usize, 0usize, false);
    // The read's evidence for ALT over REF (log10), every base weighed by its
    // quality: a low-quality base counts for little instead of fitting anything.
    let mut evidence = 0.0f64;
    let mut alt_weight = 0.0f64;
    let mut alt_spans: Vec<(usize, usize)> = Vec::new();
    let mut qual_bases: Vec<u8> = Vec::new();
    let mut read_spans: Vec<(usize, usize)> = Vec::new();
    for (rw, aw) in &win.pairs {
        let Some((r, a)) = read_pair(record, &seq, quals, min_baseq, rw, aw, map, readable) else {
            continue;
        };
        mm_ref += r.mismatches;
        mm_alt += a.mismatches;
        alt_masked += a.masked;
        // Bases matching both alleles cancel: the evidence is what mismatches REF
        // less what mismatches ALT, so the two readings need not cover the same
        // bases.
        evidence += r.mismatch_weight - a.mismatch_weight;
        alt_weight += a.mismatch_weight;
        alt_spans.push(a.span);
        had_n |= r.had_n || a.had_n;
        for (lo, hi) in [r.span, a.span] {
            qual_bases.extend_from_slice(&quals[lo..hi]);
            read_spans.push((lo, hi));
        }
        held.push((r.mismatches == 0, a.mismatches == 0));
    }
    let qual = median_qual(&qual_bases, min_baseq);
    let every = |want: (bool, bool)| !held.is_empty() && held.iter().all(|&h| h == want);
    let mut result = if every((false, true)) && evidence < one_base_evidence(min_baseq) {
        // Its clearly read bases fit the ALT, but weighed by quality its bases do not
        // favour ALT over REF by as much as one base read at the minimum quality: a
        // REF molecule read poorly around the event, with one error fitting the ALT.
        trace!("carrier: fits the ALT but its evidence {:.2} is below one base's → neither", evidence);
        ClassifyResult::neither(ClassifyPhase::MaskedCompare)
    } else if every((false, true)) {
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
    result.uninformative = held.is_empty();
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
    Some(result)
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
        // genome positions are not spliced-space ones: masks stay in genome space
        let spliced_layout = Layout { masked: Vec::new(), ..layout.clone() };
        Some(Rc::new((windows(&on_spliced, &spliced_layout)?, SpliceMap { segments })))
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
fn aligned_query_range(record: &Record) -> Option<(usize, usize)> {
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
/// event with its flank and padding, or the variant is not judged by this rule.
pub(crate) fn reference_short(start: i64, reference: &str, v: &Variant) -> Option<(bool, bool)> {
    if !judges(v) {
        return None;
    }
    let reference = upper(reference);
    let ev = event(start, &reference, v)?;
    let (need_l, need_r) = ev.room_needed();
    let short = (ev.lo < need_l, ev.hi + need_r > reference.len());
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
    /// Reference bases wanted left and right of the event: the flank, and the reach
    /// the neighbouring-change guard may grow it by.
    fn room_needed(&self) -> (usize, usize) {
        (FLANK + GUARD_REACH, FLANK + GUARD_REACH)
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

/// How a variant's windows are laid out beyond the event itself: which flank they
/// are read from, how far the neighbouring-change guard grows each flank, and the
/// genome positions it masks.
#[derive(Clone, Default)]
struct Layout {
    from_right: bool,
    grow_left: usize,
    grow_right: usize,
    masked: Vec<i64>,
}

impl Layout {
    /// The layout for `v` under `rules`: its guard (when it was found for this very
    /// variant), and the flank it is read from (the left one, unless an RNA exon
    /// edge cuts it and not the right).
    fn for_variant(v: &Variant, rules: &ReadRules) -> Layout {
        let mut layout = rules.guard.and_then(|g| g.get(v)).map_or_else(Layout::default, |g| Layout {
            from_right: false,
            grow_left: g.grow_left,
            grow_right: g.grow_right,
            masked: g.masked.clone(),
        });
        if let (Some(edges), Some((lo, hi))) = (rules.clip_edges.map(ClipEdges::get), grown_event(v)) {
            let cuts = |a: i64, b: i64| edges.iter().any(|&e| e > a && e <= b);
            let (fl, fr) = ((FLANK + layout.grow_left) as i64, (FLANK + layout.grow_right) as i64);
            if cuts(lo - fl, lo) && !cuts(hi, hi + fr) {
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
    // The flank on each side: FLANK bases, grown past a neighbouring change on the
    // ALT reads, as far as the prepared reference holds (the guard never grows past it).
    let fl = (FLANK + layout.grow_left).min(ev.lo);
    let fr = (FLANK + layout.grow_right).min(reference.len().saturating_sub(ev.hi));
    if fl < FLANK || fr < FLANK {
        return None; // the reference does not hold the event's flank
    }
    let d = ev.d;
    let lo = ev.lo - fl;
    let hi = ev.hi + fr;
    let alt_hi = (hi as i64 + d) as usize;
    let g = |off: usize| start + off as i64;
    // A masked flank base (a change the REF reads share) matches any read base.
    let (mut reference, mut hap) = (reference, ev.hap.clone());
    for &p in &layout.masked {
        let off = p - start;
        if off >= lo as i64 && off < ev.lo as i64 {
            reference[off as usize] = WILD;
            hap[off as usize] = WILD;
        } else if off >= ev.hi as i64 && off < hi as i64 {
            reference[off as usize] = WILD;
            hap[(off + d) as usize] = WILD;
        }
    }
    let (ref_win, alt_win) = (&reference[lo..hi], &hap[lo..alt_hi]);
    let (c_lo, c_hi) = ev.core;
    // Windows of one length for both alleles: the shorter allele's whole window or,
    // when both alleles are long, through one base past the first difference read
    // from the flank, from the flank nearer it (a structural choice, the same for
    // every read and both alleles).
    let short = ref_win.len().min(alt_win.len());
    let j_left = (ev.first - lo + 2).max(fl + 2).min(short);
    let j_right = (hi - ev.last + 2).max(fr + 2).min(short);
    let tail = layout.from_right || (short > LONG_EVENT && j_right < j_left);
    let j = if short <= LONG_EVENT { short } else if tail { j_right } else { j_left };
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
    Some(Windows { pairs: vec![pair], event: (g(c_lo), g(c_hi)) })
}

/// The neighbouring-change guard of one variant: how far each flank grows, the
/// genome positions masked, and the larger allele the ALT reads carry (for the
/// diagnostic), with the variant it was found for.
#[derive(Debug, Clone)]
pub(crate) struct Guard {
    key: (i64, String, String),
    pub grow_left: usize,
    pub grow_right: usize,
    pub masked: Vec<i64>,
    /// 0-based POS, REF, ALT (left-aligned, minimal VCF form) of the reads' larger
    /// allele; reads carrying it; reads carrying the ALT core with the reference
    /// base where it differs.
    pub larger: Option<(i64, String, String, u32, u32)>,
}

/// A variant's [`Guard`], found from its reads on first use: only the exact-carrier
/// rule reads it, so an SNV or a pure indel never pays for it.
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

/// The guard for a complex variant, from the reads `counted` admits: each read that
/// holds both of the event's flank bases aligned and shows one allele's core
/// exactly between them (masked bases fit) is read outward up to [`GUARD_REACH`]
/// bases on each side, and its confident changes from the reference are tallied per
/// distance. A change on at least [`GUARD_MIN_READS`] ALT reads and
/// [`GUARD_MIN_SHARE`] of those reading that base, and on under [`GUARD_REF_SHARE`]
/// of the REF reads reading it, grows that flank past it; one that recurs on REF
/// reads too (or on ALT and REF alike) is masked. None when nothing recurs, or the
/// variant is not judged by this rule.
pub(crate) fn guard_for(
    read_cache: &[Record],
    v: &Variant,
    counted: &dyn Fn(&Record) -> bool,
    min_baseq: u8,
) -> Option<Guard> {
    if !judges(v) {
        return None;
    }
    let (start, reference) = prepared_reference(v)?;
    let ev = event(start, &reference, v)?;
    if ev.lo == 0 || ev.hi >= reference.len() {
        return None;
    }
    let reach_l = GUARD_REACH.min(ev.lo);
    let reach_r = GUARD_REACH.min(reference.len() - ev.hi);
    let alt_core = &ev.hap[ev.lo..(ev.hi as i64 + ev.d) as usize];
    let ref_core = &reference[ev.lo..ev.hi];
    let g = |off: usize| start + off as i64;
    // Per allele (REF 0, ALT 1) and side (left 0, right 1), by distance from the
    // event (index 0 is the flank base beside it): reads reading the base, and reads
    // showing a confident change there.
    let mut cover = vec![vec![vec![0u32; GUARD_REACH]; 2]; 2];
    let mut change = cover.clone();
    // the bases the ALT reads show where they change, per side and distance (A C G T)
    let mut alt_bases = vec![vec![[0u32; 4]; GUARD_REACH]; 2];
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
        let between = &seq[ql + 1..qr];
        let fits = |core: &[u8]| {
            core.len() == between.len()
                && between.iter().zip(core).enumerate().all(|(i, (&b, &c))| {
                    let b = b.to_ascii_uppercase();
                    b == c || b == WILD || quals[ql + 1 + i] < min_baseq
                })
        };
        let allele = match (fits(ref_core), fits(alt_core)) {
            (true, false) => 0,
            (false, true) => 1,
            _ => continue,
        };
        for (side, reach) in [(0usize, reach_l), (1usize, reach_r)] {
            for i in 0..reach {
                let off = if side == 0 { ev.lo - 1 - i } else { ev.hi + i };
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
    let share = |n: u32, of: u32| if of == 0 { 0.0 } else { n as f64 / of as f64 };
    // per side, the distances of the changes the flank grows past
    let (mut grow, mut masked, mut grown) = ([0usize; 2], Vec::new(), [Vec::new(), Vec::new()]);
    for (side, reach) in [(0usize, reach_l), (1usize, reach_r)] {
        for i in 0..reach {
            let (alt_n, alt_of) = (change[1][side][i], cover[1][side][i]);
            let (ref_n, ref_of) = (change[0][side][i], cover[0][side][i]);
            let on_alt = alt_n >= GUARD_MIN_READS && share(alt_n, alt_of) >= GUARD_MIN_SHARE;
            let on_ref = share(ref_n, ref_of) >= GUARD_REF_SHARE;
            let recurs_on_ref = ref_n >= GUARD_MIN_READS && share(ref_n, ref_of) >= GUARD_MIN_SHARE;
            let off = if side == 0 { ev.lo - 1 - i } else { ev.hi + i };
            if on_alt && !on_ref {
                // the flank must reach the base: FLANK plus the growth covers i + 1 bases
                grow[side] = grow[side].max((i + 1).saturating_sub(FLANK));
                grown[side].push(i);
            } else if (on_alt && on_ref) || recurs_on_ref {
                masked.push(g(off));
            }
        }
    }
    if grow == [0, 0] && masked.is_empty() {
        return None;
    }
    // The larger allele: the given one with every grown change, each base as most of
    // the ALT reads show it; its carriers are the fewest showing any one change.
    let larger = if grow == [0, 0] {
        None
    } else {
        let far = |side: usize| grown[side].iter().max().copied();
        let lo = far(0).map_or(ev.lo, |i| ev.lo - 1 - i);
        let hi = far(1).map_or(ev.hi, |i| ev.hi + i + 1);
        let (mut carriers, mut given) = (u32::MAX, 0u32);
        for side in 0..2 {
            for &i in &grown[side] {
                carriers = carriers.min(change[1][side][i]);
                given = given.max(cover[1][side][i] - change[1][side][i]);
            }
        }
        let base_at = |side: usize, i: usize, off: usize| {
            if !grown[side].contains(&i) {
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
        canonical(g(lo), reference[lo..hi].to_vec(), alt, &base)
            .map(|(pos, r, a)| (pos, text(r), text(a), carriers, given))
    };
    trace!(
        "carrier guard at {}:{} {}>{}: flanks grown by {:?} past a change on the ALT reads only, {} positions masked",
        v.chrom, v.pos + 1, v.ref_allele, v.alt_allele, grow, masked.len(),
    );
    Some(Guard {
        key: (v.pos, v.ref_allele.clone(), v.alt_allele.clone()),
        grow_left: grow[0],
        grow_right: grow[1],
        masked,
        larger,
    })
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

/// The REF and ALT windows read from the read. A read that aligns the reading
/// flank's end base (the left flank's first base, or the right flank's last) reads
/// both windows from it: the windows share it, so the read is placed truly whatever
/// its own gap. A read whose soft clip holds that flank (the aligner clipped the
/// event and the flank) is placed from the other flank as each allele would sit
/// there, at its own length; a placement counts only when it falls in the clip,
/// and both windows are read from it. A read that merely starts (or ends) inside
/// the windows holds them from neither flank: depth only, for either allele alike.
/// Of the placements, one where its allele fits counts; with none, the closer
/// reading. None when the read cannot be placed or does not hold both windows' held
/// bases. Only the query bases `[first, after)` are read, and windows built in
/// spliced space are placed through `map`.
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
) -> Option<(Reading, Reading)> {
    let at_left = rw.left.and_then(|g| find_read_pos(record, map.to_genome(g))).map(|q| q as i64);
    let at_right = rw.right.and_then(|g| find_read_pos(record, map.to_genome(g - 1))).map(|q| q as i64 + 1);
    // A window is read from a flank only while it keeps a flank base there: windows
    // cut at an exon edge on the reading side are read from the other one (the
    // event's own first base is no anchor; an aligner may write it as an insertion).
    let tail = if rw.tail { rw.after > 0 } else { rw.before == 0 && rw.after > 0 };
    let read = |w: &Window, s: i64| read_at(seq, quals, min_baseq, w, s, tail, readable);
    // Both windows from one placement of the reading flank's end base: their shared
    // start (left), or their shared end (right; each starts its own length before it).
    let from = |anchor: i64| {
        let start = |w: &Window| if tail { anchor - w.allele.len() as i64 } else { anchor };
        Some((read(rw, start(rw))?, read(aw, start(aw))?))
    };
    match (tail, at_left, at_right) {
        (false, Some(s), _) => return from(s),
        (true, _, Some(e)) => return from(e),
        _ => {}
    }
    // The reading flank in the read's clip: placed from the other flank, as each
    // allele would sit, when that falls inside the clip.
    let (first, after) = aligned_query_range(record)?;
    let candidates: Vec<(i64, bool)> = [(rw, false), (aw, true)]
        .iter()
        .filter_map(|&(w, is_alt)| {
            let n = w.allele.len() as i64;
            let s = if tail { at_left.map(|s| s + n)? } else { at_right.map(|e| e - n)? };
            // the flank's end base must lie in the clip on the reading side
            let in_clip = if tail { s > after as i64 } else { s < first as i64 };
            in_clip.then_some((s, is_alt))
        })
        .collect();
    let pairs: Vec<((Reading, Reading), bool)> =
        candidates.iter().filter_map(|&(anchor, is_alt)| Some((from(anchor)?, is_alt))).collect();
    // A placement where its own allele fits counts; else the closer reading.
    pairs
        .iter()
        .find(|((r, a), is_alt)| if *is_alt { a.mismatches == 0 } else { r.mismatches == 0 })
        .or_else(|| pairs.iter().min_by_key(|((r, a), _)| r.mismatches.min(a.mismatches)))
        .map(|(p, _)| (p.0.clone(), p.1.clone()))
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
fn one_base_evidence(q: u8) -> f64 {
    let e = error_probability(q);
    ((1.0 - e) / (e / 3.0)).log10()
}


/// Mismatches of read bases against haplotype bases position by position, masked
/// bases (below `min_baseq`, or N) matching anything.
fn score(bases: &[u8], quals: &[u8], min_baseq: u8, hap: &[u8], span: (usize, usize)) -> Reading {
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
    Reading { mismatches, masked, had_n, span, mismatch_weight }
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
        assert_eq!(w.pairs.len(), 1);
        (&w.pairs[0].0, &w.pairs[0].1)
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
    fn long_events_are_read_from_one_flank_too() {
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

    /// A read of `seq` aligned at `pos` by `cigar`, at Q30.
    fn aligned(seq: &[u8], cigar: Vec<Cigar>, pos: i64) -> Record {
        let mut rec = Record::new();
        rec.set(b"r", Some(&CigarString(cigar)), seq, &vec![30; seq.len()]);
        rec.set_pos(pos);
        rec
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
        let reads: Vec<Record> = (0..4).map(|_| carrier.clone()).collect();
        let g = guard_for(&reads, &v, &|_| true, 20).expect("the change past the event");
        assert_eq!((g.grow_left, g.grow_right), (0, 3));
    }

    #[test]
    fn the_guard_names_every_change_the_alt_reads_carry_and_masks_a_shared_one() {
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
        for snp in [false, true] {
            let shared: &[usize] = if snp { &[hi + 6] } else { &[] };
            let ref_read = aligned(&with_c(&ctx, shared), vec![Cigar::Match(98)], 0);
            let reads: Vec<Record> = (0..4).map(|_| carrier.clone()).chain((0..4).map(|_| ref_read.clone())).collect();
            let g = guard_for(&reads, &v, &|_| true, 20).expect("a recurring change");
            let (pos, r, a, n, m) = g.larger.clone().expect("a larger allele");
            // left-aligned and minimal: the grown event's shared bases are trimmed
            if snp {
                // the shared C is masked; the larger allele stops at the other one
                assert_eq!(g.masked, vec![(hi + 6) as i64]);
                assert_eq!((pos, r.as_str(), n, m), (30, &ctx[30..hi + 4], 4, 0));
                assert_eq!(a, format!("TTCTCT{}C", &ctx[58..hi + 3]));
            } else {
                assert!(g.masked.is_empty());
                assert_eq!((pos, r.as_str(), n, m), (30, &ctx[30..hi + 7], 4, 0));
                assert_eq!(a, format!("TTCTCT{}C{}C", &ctx[58..hi + 3], &ctx[hi + 4..hi + 6]));
            }
        }
    }
}
