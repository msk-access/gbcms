//! The haplotypes a row X's reads are weighed against, built once per variant from its
//! prepared rows, its reads and its guard:
//!
//! - **The window:** the reference over X and every other listed row within a read's
//!   reach of it, [`READ_REACH`](super::READ_REACH) bases each side.
//! - **The site:** the states X's site can hold: REF, X, every listed row at the site,
//!   and, for an SNV, every other base; for an MNP, an indel or a delins, the other
//!   alleles at the locus the guard found in the reads. Two rows that both change the
//!   length are at one site when any of their equivalent placements meet (an insertion
//!   left-aligned one base before another's anchor is still at its place); otherwise
//!   when their REF spans meet (an insertion's span is its anchor base).
//! - **The contexts:** the other listed rows within the reads' reach, with the guard's
//!   recurring background changes (they explain reads but are never a site state). A
//!   listed length change whose run X lies inside is also placed at the run's other end:
//!   on the reference the two placements are one allele, beside X two molecules.
//! - **The graphs:** one per state, every context row a bubble a read's path may take or
//!   skip ([`super::graph`]), so every combination of them is weighed at once; but for a
//!   context row at the state's place, and but for another allele that spells X's or
//!   REF's sequence (the reads' spelling of X beside a listed neighbour is X's). Every
//!   context lies with REF and with X, so only a read's bases at X's site can tell those
//!   apart.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use rust_htslib::bam::record::Cigar;

use super::graph::{self, Graph, GraphBuilder, GraphRead, ShortReads, Workspace};
use super::score::{self, Alignment, Band, Placed, BAND_SLACK, DELTA, WILD};
use crate::counting::carrier::{GUARD_MIN_READS, GUARD_REF_SHARE};


/// A listed allele the haplotypes apply: its 0-based genome position and its REF and
/// ALT bases, upper case.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct Row {
    pub(crate) pos: i64,
    pub(crate) ref_al: Vec<u8>,
    pub(crate) alt_al: Vec<u8>,
}

impl Row {
    pub(crate) fn new(pos: i64, ref_al: &[u8], alt_al: &[u8]) -> Row {
        Row { pos, ref_al: ref_al.to_ascii_uppercase(), alt_al: alt_al.to_ascii_uppercase() }
    }

    fn end(&self) -> i64 {
        self.pos + self.ref_al.len() as i64
    }

    fn changes_length(&self) -> bool {
        self.ref_al.len() != self.alt_al.len()
    }
}

/// What a haplotype carries at X's site.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Class {
    /// X itself.
    Alt = 0,
    /// The reference's bases.
    Ref = 1,
    /// Another allele at the site.
    Other = 2,
}

/// A read as the judge reads it: its readable bases (masked to [`WILD`] below
/// `--min-baseq` or where not A, C, G or T) with their qualities, each base's reference
/// position (a soft clip's extrapolated, an inserted base's at its anchor), and its
/// aligned span.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ReadView {
    pub(crate) seq: Vec<u8>,
    pub(crate) quals: Vec<u8>,
    refpos: Vec<i64>,
    ref_start: i64,
    ref_end: i64,
    /// The insertions its own alignment writes: the genome position after the anchor,
    /// and the inserted bases (upper case, as read).
    pub(crate) inserts: Vec<(i64, Vec<u8>)>,
    /// Per indel its own alignment writes, the bases (indices into `seq`) on either side
    /// of it: the last before it and the first after it (either may lie outside `seq`).
    indels: Vec<(isize, isize)>,
    /// Per base of `seq`, whether its alignment places it on a reference base.
    aligned: Vec<bool>,
    /// The read's name (its fragment), empty where not given.
    pub(crate) name: Vec<u8>,
    /// Its aligner's own evidence on where it belongs, for the "fits nothing" rule.
    pub(crate) tags: AlnTags,
}

/// What a read's aligner wrote about its placement: MAPQ, the AS and XS tags (the best
/// alignment's score and the best score elsewhere, BWA-MEM's scale) and NM. None where
/// the record carries no such tag.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub(crate) struct AlnTags {
    pub(crate) mapq: u8,
    pub(crate) best: Option<i64>,
    pub(crate) elsewhere: Option<i64>,
    pub(crate) nm: Option<i64>,
}

impl AlnTags {
    /// Whether the aligner places the read elsewhere at least as well as here (XS at
    /// least AS); where it wrote no XS, whether it gave the read MAPQ 0.
    pub(crate) fn fits_elsewhere(&self) -> bool {
        match (self.best, self.elsewhere) {
            (Some(a), Some(x)) => x >= a,
            _ => self.mapq == 0,
        }
    }
}

/// Some of a view's fields (views equal in all are equal in these), for the judge's maps of
/// views alike (`judge::Haplotypes`), which compare views whole.
impl std::hash::Hash for ReadView {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.seq.hash(state);
        self.quals.hash(state);
        self.ref_start.hash(state);
        self.ref_end.hash(state);
        self.tags.hash(state);
    }
}

impl ReadView {
    /// The readable bases `[first, after)` of a read aligned at `ref_start` with `cigar`.
    /// None when none is readable.
    pub(crate) fn new(
        seq: &[u8],
        quals: &[u8],
        ref_start: i64,
        cigar: &[Cigar],
        (first, after): (usize, usize),
        min_baseq: u8,
    ) -> Option<ReadView> {
        let mut refpos = Vec::with_capacity(seq.len());
        let mut aligned = Vec::with_capacity(seq.len());
        let (mut inserts, mut indels) = (Vec::new(), Vec::new());
        let mut r = ref_start;
        // the read's bases either side of an indel, as indices into the readable bases
        let side = |q: usize| q as isize - first as isize;
        for op in cigar {
            match *op {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    for _ in 0..n {
                        refpos.push(r);
                        aligned.push(true);
                        r += 1;
                    }
                }
                Cigar::Ins(n) => {
                    let q = refpos.len();
                    if let Some(bases) = seq.get(q..q + n as usize) {
                        inserts.push((r, bases.to_ascii_uppercase()));
                    }
                    indels.push((side(q) - 1, side(q + n as usize)));
                    refpos.extend(std::iter::repeat_n((r - 1).max(ref_start), n as usize));
                    aligned.extend(std::iter::repeat_n(false, n as usize));
                }
                Cigar::SoftClip(n) => {
                    let n = n as i64;
                    if refpos.is_empty() {
                        refpos.extend((0..n).map(|t| r - (n - t)));
                    } else {
                        refpos.extend((0..n).map(|t| r + t));
                    }
                    aligned.extend(std::iter::repeat_n(false, n as usize));
                }
                Cigar::Del(n) => {
                    let q = refpos.len();
                    indels.push((side(q) - 1, side(q)));
                    r += n as i64
                }
                Cigar::RefSkip(n) => r += n as i64,
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }
        // A CIGAR that disagrees with the bases' length: pad or cut, never fail.
        refpos.resize(seq.len(), r);
        aligned.resize(seq.len(), false);
        let after = after.min(seq.len()).min(quals.len());
        if first >= after {
            return None;
        }
        let masked = seq[first..after]
            .iter()
            .zip(&quals[first..after])
            .map(|(&b, &q)| {
                let b = b.to_ascii_uppercase();
                if q < min_baseq || !matches!(b, b'A' | b'C' | b'G' | b'T') {
                    WILD
                } else {
                    b
                }
            })
            .collect();
        Some(ReadView {
            seq: masked,
            quals: quals[first..after].to_vec(),
            refpos: refpos[first..after].to_vec(),
            ref_start,
            ref_end: r.max(ref_start + 1),
            inserts,
            indels,
            aligned: aligned[first..after].to_vec(),
            name: Vec::new(),
            tags: AlnTags::default(),
        })
    }

    /// Measurement only: the readable bases its own alignment leaves unplaced at
    /// each end (its soft clips there), leading and trailing.
    pub(crate) fn clipped_ends(&self) -> (usize, usize) {
        let lead = self.aligned.iter().take_while(|&&a| !a).count();
        let tail = if lead == self.aligned.len() { 0 } else { self.aligned.iter().rev().take_while(|&&a| !a).count() };
        (lead, tail)
    }

    /// The view with the read's name (its fragment).
    pub(crate) fn named(mut self, name: &[u8]) -> Self {
        self.name = name.to_vec();
        self
    }

    /// The view with its aligner's placement evidence.
    pub(crate) fn tagged(mut self, tags: AlnTags) -> Self {
        self.tags = tags;
        self
    }

    /// Its bases read confidently (not masked).
    pub(crate) fn confident(&self) -> usize {
        self.seq.iter().filter(|&&b| b != WILD).count()
    }

    /// The confident mismatches its base qualities lead one to expect: the sum of each
    /// confident base's error probability.
    pub(crate) fn expected_errors(&self) -> f64 {
        self.seq.iter().zip(&self.quals).filter(|(&b, _)| b != WILD).map(|(_, &q)| 10f64.powf(-f64::from(q) / 10.0)).sum()
    }

    /// The allele the read's own bases spell at each indel its alignment writes: from the
    /// last base before it to the first after it that start [`SPELL_STOP`] bases reading
    /// the reference's (read confidently, aligned), the read's bases between them, as a
    /// row anchored on the base before. A mismatch the aligner wrote beside the gap is part of the allele (a
    /// delins written as a deletion and a substitution is that delins). None for an indel
    /// whose bases there are not all read. A base the guard masks (`masked`, genome
    /// positions sorted: a germline change beside the event) is no part of the allele.
    pub(crate) fn spelled_indels(&self, base: &dyn Fn(i64) -> Option<u8>, masked: &[i64]) -> Vec<Row> {
        let n = self.seq.len() as isize;
        let fits = |i: isize| (0..n).contains(&i) && self.aligned[i as usize] && self.seq[i as usize] != WILD;
        let matches = |i: isize| {
            let p = self.refpos[i as usize];
            base(p) == Some(self.seq[i as usize]) || masked.binary_search(&p).is_ok()
        };
        // The allele ends at SPELL_STOP read bases in a row that read the reference's (a
        // base of the allele may match the reference by chance).
        let stops = |i: isize, step: isize| (0..SPELL_STOP as isize).all(|k| fits(i + k * step) && matches(i + k * step));
        self.indels
            .iter()
            .filter_map(|&(l, r)| {
                let (mut l, mut r) = (l, r);
                while fits(l) && !stops(l, -1) {
                    l -= 1;
                }
                while fits(r) && !stops(r, 1) {
                    r += 1;
                }
                if !fits(l) || !fits(r) {
                    return None;
                }
                let bases = &self.seq[(l + 1) as usize..r as usize];
                if bases.contains(&WILD) {
                    return None;
                }
                let (a, b) = (self.refpos[l as usize], self.refpos[r as usize]);
                let reference: Option<Vec<u8>> = (a..b).map(base).collect();
                let reference = reference?;
                Some(Row::new(a, &reference, &[&reference[..1], bases].concat()))
            })
            .collect()
    }

    pub(crate) fn len(&self) -> usize {
        self.seq.len()
    }

    /// The view without its qualities, its aligner's tags and its name: its bases (masked)
    /// and placement alone, which a read's programmes and alignments are fixed by.
    pub(crate) fn placed_alone(&self) -> ReadView {
        ReadView { quals: Vec::new(), tags: AlnTags::default(), name: Vec::new(), ..self.clone() }
    }

    /// The genome position its own alignment places base `i` at; None for a base it leaves
    /// unplaced (inserted or clipped).
    pub(crate) fn own_position(&self, i: usize) -> Option<i64> {
        (i < self.seq.len() && self.aligned[i]).then(|| self.refpos[i])
    }

    /// The base its own alignment places at genome position `g`, read confidently; None
    /// where it places none there or reads it below `--min-baseq`.
    fn aligned_base(&self, g: i64) -> Option<u8> {
        (0..self.seq.len()).find(|&i| self.aligned[i] && self.refpos[i] == g).map(|i| self.seq[i]).filter(|&b| b != WILD)
    }
}

/// The window with a set of rows applied that can be applied together.
#[derive(Debug, Clone)]
pub(crate) struct Hap {
    /// The rows applied (indices into the site's rows), ascending.
    rows: Vec<usize>,
    pub(crate) class: Class,
    seq: Vec<u8>,
    /// Per applied row, ascending: its window-relative REF span and haplotype span.
    segs: Vec<(i64, i64, i64, i64)>,
    /// X's allele's haplotype span, where X is applied.
    x_seg: Option<(i64, i64)>,
}

impl Hap {
    /// Measurement only: the haplotype's base at column `c` (a masked base is WILD).
    pub(crate) fn base(&self, c: usize) -> u8 {
        self.seq.get(c).copied().unwrap_or(WILD)
    }

    fn build(site: &Site, rows: &[usize]) -> Hap {
        let mut order: Vec<usize> = rows.to_vec();
        order.sort_by_key(|&r| site.spans[r].0);
        let mut seq = site.reference.clone();
        // A masked base takes any base where no applied row covers it.
        for &m in &site.masked {
            if !order.iter().any(|&r| site.spans[r].0 <= m && m < site.spans[r].1) {
                seq[m as usize] = WILD;
            }
        }
        for &r in order.iter().rev() {
            let (s, e) = site.spans[r];
            seq.splice(s as usize..e as usize, site.row_bases(r));
        }
        let mut segs = Vec::with_capacity(order.len());
        let mut x_seg = None;
        let mut delta = 0i64;
        for &r in &order {
            let (s, e) = site.spans[r];
            let hs = s + delta;
            let he = hs + site.rows[r].alt_al.len() as i64;
            delta += site.rows[r].alt_al.len() as i64 - (e - s);
            segs.push((s, e, hs, he));
            if r == 0 {
                x_seg = Some((hs, he));
            }
        }
        let class = if rows.contains(&0) || rows.iter().any(|&r| site.holds_x[r]) {
            Class::Alt
        } else if rows.iter().any(|&r| site.is_state(r)) {
            Class::Other
        } else {
            Class::Ref
        };
        let mut rows = rows.to_vec();
        rows.sort_unstable();
        Hap { rows, class, seq, segs, x_seg }
    }

    /// Window coordinate → haplotype coordinate. One inside an applied row maps to the
    /// row's haplotype start (for a start) or end (for an end), so a trimmed haplotype
    /// always holds an applied allele whole.
    fn to_hap(&self, p: i64, is_end: bool) -> i64 {
        let mut delta = 0;
        for &(rs, re, hs, he) in &self.segs {
            if p >= re {
                delta = he - re;
                continue;
            }
            if p > rs {
                return if is_end { he } else { hs };
            }
            break;
        }
        p + delta
    }

    /// The window's length.
    fn window_len(&self) -> i64 {
        self.seq.len() as i64 - self.segs.iter().map(|&(rs, re, hs, he)| (he - hs) - (re - rs)).sum::<i64>()
    }

    /// The part a read can reach: its aligned span `[s, e)` (window coordinates) mapped
    /// into the haplotype, widened by `l` haplotype bases each side.
    fn reach(&self, s: i64, e: i64, l: i64) -> (usize, usize) {
        let (n, wn) = (self.seq.len() as i64, self.window_len());
        let a = (self.to_hap(s.clamp(0, wn), false) - l).clamp(0, n) as usize;
        let b = (self.to_hap(e.clamp(0, wn), true) + l).clamp(0, n) as usize;
        (a, b.max(a))
    }
}

/// X's haplotypes over its window.
#[derive(Debug, Clone)]
pub(crate) struct Site {
    /// Genome start of the window, and its reference.
    start: i64,
    reference: Vec<u8>,
    /// X (index 0), the other site states, then the context rows.
    rows: Vec<Row>,
    n_states: usize,
    /// Per row, whether it is a state holding X's change at X's site (a substitution X:
    /// an allele of the same length over X's bases with X's bases there, such as an MNV
    /// listed beside its component SNV). Its haplotypes are X's ALT (decided 2026-10-25).
    holds_x: Vec<bool>,
    /// Window-relative REF spans of `rows`, and the spans their equivalent placements
    /// cover (a length change's; any other row's REF span).
    spans: Vec<(i64, i64)>,
    places: Vec<(i64, i64)>,
    /// The most any haplotype's rows change the window's length by.
    max_shift: i64,
    /// Window positions every haplotype takes any base at (the guard's masks), sorted.
    masked: Vec<i64>,
    /// X's site, `[s, e)` window coordinates: the bases it changes and, for a length
    /// change, every place it can sit and the repeat it changes, a base each side.
    site_window: (i64, i64),
    /// The graph scorer's graphs, one per state of X's site: REF's bases, X, then the
    /// other states.
    pub(crate) graphs: Vec<StateGraph>,
    /// What a covered telling base read confidently as neither allele does to a REF or ALT
    /// call ([`neither_rule`]).
    pub(crate) neither: NeitherRule,
    /// Which states a read reaching them only through an indel of its own loses
    /// ([`own_gap_rule`]).
    pub(crate) own_gap: OwnGapRule,
    /// What a confident base of the other allele, or a third, at X's telling bases does to
    /// a call ([`ref_base_rule`]).
    pub(crate) ref_base: RefBaseRule,
    /// At a substitution of two bases or more or a delins X: per genome position over the
    /// window and base, the reads whose own alignments read that base there confidently.
    tally: HashMap<(i64, u8), u32>,
    /// The guard's background, genome coordinates: its recurring single-base changes, and
    /// the positions it masks (any base). A base it explains is no allele's third base.
    background_bases: HashSet<(i64, u8)>,
    background_masks: Vec<i64>,
    /// Where a pure indel X's placements touch a contig edge (its repeat running to it, or X
    /// at the contig's first or last base: the base that tells its lengths apart lies past
    /// the edge, where a circular contig runs on at its start): that edge, genome
    /// coordinates, as (first base, one past the last).
    edge: (Option<i64>, Option<i64>),
}

/// Measurement only, TEMPORARY: what a covered telling base read confidently as neither
/// of the row's alleles there does to a REF or ALT call on a substitution of two bases or
/// more or a delins (`GBCMS_JUDGE_NEITHER`): `margin` (the default, decided 2026-10-24: the
/// margin alone decides, so a read one confident base from REF there is REF, and an ALT
/// call needing every changed base read leaves it depth only), or `other` (such a base
/// makes the read another allele, partial evidence, as an MNV read whose block is neither
/// allele counts under the rules and in GetBaseCounts). An SNV's third base is a state of
/// its site under either. Read when a site is built.
///
/// `other:N` moves such reads only at rows whose allele is under N bases ([`Site::x_size`]);
/// `recur:K` only where a neither base it reads is read the same at the same place by K
/// other reads of the row, by their own alignments ([`Site::recurrence`]); `window` only where
/// a neither base lies where both alleles' reads are read alike: a base a substitution changes
/// (not one it keeps, where the alleles agree), or within [`JUNCTION_WINDOW`] telling bases of
/// either end of a length change's (its junctions, as RJ-24 reads a long deletion; deeper in,
/// REF's molecules hold bases X's never do, so an error there would weigh on REF alone).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum NeitherRule {
    Margin,
    Other,
    Below(usize),
    Recur(u32),
    Window,
}

pub(crate) fn neither_rule() -> NeitherRule {
    let v = std::env::var("GBCMS_JUDGE_NEITHER").unwrap_or_default();
    match v.split_once(':') {
        None if v == "other" => NeitherRule::Other,
        None if v == "window" => NeitherRule::Window,
        Some(("other", n)) => n.parse().map_or(NeitherRule::Margin, NeitherRule::Below),
        Some(("recur", k)) => k.parse().map_or(NeitherRule::Margin, NeitherRule::Recur),
        _ => NeitherRule::Margin,
    }
}

/// Measurement only, TEMPORARY: with `GBCMS_JUDGE_NEITHER_TRACE=1`, every REF or ALT call on
/// a substitution of two bases or more or a delins traces what the read shows at X's telling
/// bases that is neither allele's (each such base's quality, place and recurrence), and its
/// mismatches elsewhere on the called haplotype ("judge neither" lines).
/// Measurement only, TEMPORARY: which states a read loses that it reaches only through an
/// indel of its own at the state's change (`GBCMS_JUDGE_OWNGAP`): `tie`, the default, another
/// allele tied with REF where X is beaten (the decided tie rule); `all`, any state (REF, X or
/// another allele) within the margin of the read's best, before the margin decides, for every
/// pair of states and for wins as for ties; `others`, as `all` for another allele's states
/// only, REF and X staying for the margin to weigh against each other. Read when a site is
/// built.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum OwnGapRule {
    Tie,
    All,
    Others,
}

pub(crate) fn own_gap_rule() -> OwnGapRule {
    match std::env::var("GBCMS_JUDGE_OWNGAP").as_deref() {
        Ok("all") => OwnGapRule::All,
        Ok("others") => OwnGapRule::Others,
        _ => OwnGapRule::Tie,
    }
}

/// Measurement only, TEMPORARY: what a confident base at one of X's telling bases that is
/// not the called allele's does to a REF or ALT call at a substitution of two bases or more
/// or a delins (`GBCMS_JUDGE_REFBASE`). `margin`, the default: the margin alone decides a
/// delins call, and an MNP's ALT call needs every base X changes read as X's (a third base
/// there leaves it depth only). `partial`: the other allele's base there (under the read's
/// own placement: at a delins, the other allele's base at the same telling index counted
/// from either end) is that allele's evidence, not an error, so the read is another
/// allele (RJ-25, RJ-28); a third base, neither allele's, is an error the margin forgives
/// (RJ-26), at an MNP as at a delins. Read when a site is built.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RefBaseRule {
    Margin,
    Partial,
}

pub(crate) fn ref_base_rule() -> RefBaseRule {
    match std::env::var("GBCMS_JUDGE_REFBASE").as_deref() {
        Ok("partial") => RefBaseRule::Partial,
        _ => RefBaseRule::Margin,
    }
}

/// Measurement only: whether to trace each background item the respelling check weighs
/// (`GBCMS_JUDGE_RESPELL_TRACE=1`): its kind, direction, state and fate.
pub(crate) fn respell_trace() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GBCMS_JUDGE_RESPELL_TRACE").is_ok_and(|v| v == "1"))
}

pub(crate) fn neither_trace() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GBCMS_JUDGE_NEITHER_TRACE").is_ok_and(|v| v == "1"))
}

/// One state of X's site and its graph: the window with the state's row applied (none for
/// REF's bases) and every context row in reach that can lie with it as a bubble.
#[derive(Debug, Clone)]
pub(crate) struct StateGraph {
    pub(crate) class: Class,
    /// The state's row (X is 0), None for REF's bases.
    pub(crate) row: Option<usize>,
    pub(crate) graph: Graph,
    /// An earlier state whose graph is laid out alike, and the first base where their bases
    /// differ ([`Graph::layout_until`]): a substitution's site states differ at its bases
    /// alone, so a read's programme takes the columns before them from that state's
    /// ([`graph::decided_states`]). None for a state laid out as no earlier one.
    pub(crate) share: Option<(usize, usize)>,
}

/// The rows a site is built from, besides X: the input's other rows within reach (a
/// state where at X's site, else a context), the other alleles the reads show at the
/// locus (states), and the background changes recurring there (contexts only).
#[derive(Debug, Default, Clone)]
pub(crate) struct Candidates {
    pub(crate) listed: Vec<Row>,
    pub(crate) discovered: Vec<Row>,
    pub(crate) background: Vec<Row>,
    /// Genome positions where a change recurs beside the event that may take any base
    /// (the guard's masks): every haplotype takes any base there, unless a row it applies
    /// covers it, so a read's base there tells no haplotype from another.
    pub(crate) masked: Vec<i64>,
    /// Indels the reads' own alignments write recurringly (any length), each with the
    /// reads (indices into the site's reads) that write it: states where they lie at X's
    /// site and those reads' bases tell them from X and REF, else left out.
    pub(crate) found: Vec<(Row, Vec<usize>)>,
    /// The contig's edges inside the window, genome coordinates: its first base where the
    /// window would run past its start, one past its last where it would run past its end.
    pub(crate) edges: (Option<i64>, Option<i64>),
}

/// Measurement only, TEMPORARY: which of the judge's other-allele options are on
/// (`GBCMS_JUDGE_OTHER`, comma-separated: `runs` the other lengths of X's repeat as
/// states, `found` the indels the reads write at X's site as states (`foundK`: written by
/// K reads, default [`FOUND_MIN_READS`]; `found3` the guard's floor), `told` or `toldK`
/// for which of them stay states ([`FoundRule`]), `fit` the absolute-fit check, `edges` the check
/// passing over indels inside a long event's own bases; `none`; `both` runs, found and
/// fit). Default `found,fit,told,edges`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct OtherOption {
    pub(crate) runs: bool,
    pub(crate) found: bool,
    pub(crate) fit: bool,
    /// Which found alleles stay states ([`Site::drop_graph_spellings`]).
    pub(crate) rule: FoundRule,
    /// The fit check passes over an indel wholly inside a long event's own bases
    /// ([`Site::long_interior`]).
    pub(crate) edges: bool,
    /// The reads that must write an indel (and [`GUARD_REF_SHARE`] of the site's) for it to
    /// be found (`foundK`; default [`FOUND_MIN_READS`]).
    pub(crate) found_min: u32,
}

/// Reads whose own alignments write an indel (with [`GUARD_REF_SHARE`] of the site's
/// reads) that make it a candidate state: two reads writing one allele make it a
/// candidate, as two reads showing an allele do for the community's callers (GATK's
/// assembly-graph pruning, FreeBayes's `--min-alternate-count`); whether it is a state is
/// the reads' telling it from X and REF ([`FoundRule`]), not the aligner's writing it. A
/// sibling allele whose carriers mostly end inside a repeat, or soft-clip its end, is
/// written whole by few reads but told apart by many.
pub(crate) const FOUND_MIN_READS: u32 = 2;

/// Measurement only, TEMPORARY: which found alleles stay site states. A read ending inside
/// X's own bases (a long insertion's first bases, one copy of a deleted duplication, a
/// repeat X lengthens) writes the part of X it holds as another indel; that spelling is
/// X's, never another allele, while a sibling allele its few reads spell partially is
/// another allele all the same. A read *tells* a found allele where its best haplotype
/// holding it beats every haplotype carrying X and every one holding REF's bases by
/// [`DELTA`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum FoundRule {
    /// Every found allele at X's site is a state (no token).
    All,
    /// `told`: the found floor of reads ([`GUARD_MIN_READS`] and [`GUARD_REF_SHARE`] of the
    /// site's) tell it; `toldK`: K reads.
    Told(usize),
}

pub(crate) fn other_option() -> OtherOption {
    static OPT: std::sync::OnceLock<OtherOption> = std::sync::OnceLock::new();
    *OPT.get_or_init(|| parse_other(&std::env::var("GBCMS_JUDGE_OTHER").unwrap_or_else(|_| "found,fit,told,edges".into())))
}

fn parse_other(v: &str) -> OtherOption {
    let v = match v {
        "both" => "runs,found,fit".to_string(),
        "states" => "runs,found".to_string(),
        other => other.to_string(),
    };
    let tokens: Vec<&str> = v.split(',').map(str::trim).collect();
    let on = |k: &str| tokens.contains(&k);
    let numbered = |stem: &str| {
        tokens.iter().find_map(|t| t.strip_prefix(stem).and_then(|n| if n.is_empty() { Some(0) } else { n.parse().ok() }))
    };
    let rule = numbered("told").map_or(FoundRule::All, FoundRule::Told);
    let found = numbered("found");
    let found_min = found.filter(|&k| k > 0).map_or(FOUND_MIN_READS, |k| k as u32);
    OtherOption { runs: on("runs"), found: found.is_some(), fit: on("fit"), rule, edges: on("edges"), found_min }
}

/// An event's own bases from which an indel inside them, away from its ends, is a change
/// on its molecule rather than another allele of it (the length from which the rules
/// judge a deletion at its junction and expect no read to hold an allele whole).
const LONG_EVENT: i64 = crate::counting::junction::LONG_DELETION as i64;
/// Bases inside a long event's ends where an indel still meets its junction.
const EDGE: i64 = 1;
/// A long deletion's deleted bases, past its microhomology, that a read over its anchor
/// reads to tell REF from it (RJ-24's REF window: four clean bases of the deleted bases).
pub(crate) const JUNCTION_WINDOW: i64 = 4;

/// Read bases in a row reading the reference that end an allele the read spells.
const SPELL_STOP: usize = 3;

/// Masks this many bases in a row or more are a shifted molecule's, not applied.
const MASK_RUN: usize = 3;

/// The longest repeat the most units may be added to or taken from: unit length.
const RUN_MAX_UNIT: usize = 6;
/// Units added to or taken from a run for the other lengths of it a site holds.
const RUN_STEPS: [i64; 4] = [-2, -1, 1, 2];

/// The longest tandem run of unit `p` in `h` with a unit starting or ending at the
/// boundary `edge` (before base `edge`): `[lo, hi)`, at least that one unit long.
fn run_through(h: &[u8], edge: usize, p: usize) -> (usize, usize) {
    let grow = |mut lo: usize, mut hi: usize| {
        while lo >= 1 && h[lo - 1] == h[lo - 1 + p] {
            lo -= 1;
        }
        while hi < h.len() && h[hi] == h[hi - p] {
            hi += 1;
        }
        (lo, hi)
    };
    let right = (edge + p <= h.len()).then(|| grow(edge, edge + p));
    let left = (edge >= p).then(|| grow(edge - p, edge));
    [right, left].into_iter().flatten().max_by_key(|&(lo, hi)| hi - lo).unwrap_or((edge, edge))
}

/// The repeat X's length change lengthens or shortens: the longest run (unit `p` up to
/// [`RUN_MAX_UNIT`], dividing the change, two copies or more) on the longer allele that
/// holds as many of the bases where the alleles differ as the change adds, and the run
/// of that unit through the same place on the shorter allele. Per allele (REF, then X),
/// its run `[lo, hi)` in that allele's coordinates, with the unit length. None when X
/// keeps the length or lies in no such repeat.
fn x_run(x: &Row, start: i64, w: &[u8]) -> Option<([(usize, usize); 2], usize)> {
    let d = x.alt_al.len() as i64 - x.ref_al.len() as i64;
    if d == 0 {
        return None;
    }
    let (xs, xe) = ((x.pos - start) as usize, (x.end() - start) as usize);
    let (on_ref, on_alt) = telling(x);
    let hx: Vec<u8> = [&w[..xs], &x.alt_al[..], &w[xe..]].concat();
    let alleles = [(w, xs + on_ref.start, xs + on_ref.end), (&hx[..], xs + on_alt.start, xs + on_alt.end)];
    let (long, short) = if d > 0 { (1, 0) } else { (0, 1) };
    let k = d.unsigned_abs() as usize;
    let (h, a, b) = alleles[long];
    let mut best: Option<(usize, usize, usize)> = None;
    for p in (1..=RUN_MAX_UNIT).filter(|p| k.is_multiple_of(*p)) {
        for edge in [a, b] {
            let (lo, hi) = run_through(h, edge, p);
            let held = hi.min(b).saturating_sub(lo.max(a));
            if hi - lo >= 2 * p && hi - lo >= 3 && held >= k && best.is_none_or(|(l, e, _)| hi - lo > e - l) {
                best = Some((lo, hi, p));
            }
        }
    }
    let (lo, hi, p) = best?;
    let (hs, sa, sb) = alleles[short];
    let on_short = [sa, sb].into_iter().map(|edge| run_through(hs, edge, p)).max_by_key(|&(l, e)| e - l)?;
    let mut runs = [(0, 0); 2];
    runs[long] = (lo, hi);
    runs[short] = on_short;
    Some((runs, p))
}

/// The other lengths of the repeat X's change lengthens or shortens ([`x_run`]), on REF
/// and on X: each allele with one or two units of it more or fewer, as rows over the
/// window. Another length of the repeat is another allele (partial evidence), never
/// REF or X.
fn run_alleles(x: &Row, start: i64, w: &[u8]) -> Vec<Row> {
    let Some((runs, p)) = x_run(x, start, w) else {
        return Vec::new();
    };
    let (xs, xe) = ((x.pos - start) as usize, (x.end() - start) as usize);
    let hx: Vec<u8> = [&w[..xs], &x.alt_al[..], &w[xe..]].concat();
    let mut out = Vec::new();
    for (h, (lo, hi)) in [(w, runs[0]), (&hx[..], runs[1])] {
        for k in RUN_STEPS {
            let alt: Vec<u8> = if k > 0 {
                if hi - lo < p {
                    continue;
                }
                let unit = &h[lo..lo + p];
                [&h[..lo], &unit.repeat(k as usize)[..], &h[lo..]].concat()
            } else {
                let cut = p * k.unsigned_abs() as usize;
                if hi - lo < cut {
                    continue;
                }
                [&h[..lo], &h[lo + cut..]].concat()
            };
            if let Some(r) = row_between(w, &alt, start) {
                out.push(r);
            }
        }
    }
    out
}

/// The longest unit an imperfect repeat's period is looked for up to.
const APPROX_MAX_UNIT: usize = 50;

/// The repeat, perfect or not, holding the bases `[a, b)` of `h`: the longest stretch
/// around them where each base mostly repeats the one a unit before it (a unit of up to
/// [`APPROX_MAX_UNIT`] bases; a unit of 5 or more may differ at one base in a unit's
/// length, as the units of a VNTR do), the smallest unit reaching it; at least two units
/// past the bases. `[lo, hi)`, or None.
fn approx_repeat(h: &[u8], a: usize, b: usize) -> Option<(usize, usize)> {
    let k = b - a;
    let mut best: Option<(usize, usize)> = None;
    for p in 1..=APPROX_MAX_UNIT.min(h.len() / 2) {
        let slack = usize::from(p >= 5);
        // right: base i repeats base i - p
        let mut hi = b.max(p);
        let mut misses: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
        let mut last_good = hi;
        while hi < h.len() {
            if h[hi] != h[hi - p] {
                misses.push_back(hi);
            }
            while misses.front().is_some_and(|&m| m + p <= hi) {
                misses.pop_front();
            }
            if misses.len() > slack {
                break;
            }
            if h[hi] == h[hi - p] {
                last_good = hi + 1;
            }
            hi += 1;
        }
        let hi = last_good;
        // left: base j repeats base j + p
        let mut lo = a;
        misses.clear();
        let mut first_good = lo;
        while lo >= 1 && lo - 1 + p < h.len() {
            let j = lo - 1;
            if h[j] != h[j + p] {
                misses.push_back(j);
            }
            while misses.front().is_some_and(|&m| m >= j + p) {
                misses.pop_front();
            }
            if misses.len() > slack {
                break;
            }
            if h[j] == h[j + p] {
                first_good = j;
            }
            lo -= 1;
        }
        let lo = first_good;
        // the changed bases must follow the period too
        let inside = (a.max(p)..b).filter(|&i| h[i] == h[i - p]).count();
        let checked = b.saturating_sub(a.max(p));
        if checked > 0 && inside * 10 < checked * 9 {
            continue;
        }
        if hi - lo >= k + 2 * p && best.is_none_or(|(l, e)| hi - lo > e - l) {
            best = Some((lo, hi));
        }
    }
    best
}

/// X's site in window coordinates `[s, e)`, where another indel in a read makes it
/// another allele (RJ-7; #290: another length of the repeat X lengthens, wherever it is
/// written in it): for a length change its discrimination window (every place it can
/// sit, `place` in genome coordinates, with a base each side), the repeat it lengthens
/// or shortens ([`x_run`]) and the repeat, perfect or not, its changed bases lie in
/// ([`approx_repeat`]), on REF, with a base each side; else its REF span.
fn x_site(x: &Row, start: i64, w: &[u8], place: (i64, i64)) -> (i64, i64) {
    let (xs, xe) = (x.pos - start, x.end() - start);
    if !x.changes_length() {
        return (xs, xe);
    }
    let (mut s, mut e) = (place.0 - start, place.1 - start + 1);
    if let Some((runs, _)) = x_run(x, start, w) {
        s = s.min(runs[0].0 as i64 - 1);
        e = e.max(runs[0].1 as i64 + 1);
    }
    // the repeat on the longer allele, in REF coordinates
    let d = x.alt_al.len() as i64 - x.ref_al.len() as i64;
    let (on_ref, on_alt) = telling(x);
    let (xs_u, xe_u) = (xs as usize, xe as usize);
    // REF's repeat at the bases X changes (for an insertion, at its junction), and the
    // longer allele's, in REF coordinates.
    let mut found = vec![approx_repeat(w, xs_u + on_ref.start, xs_u + on_ref.end).map(|(lo, hi)| (lo as i64, hi as i64))];
    if d > 0 {
        let hx: Vec<u8> = [&w[..xs_u], &x.alt_al[..], &w[xe_u..]].concat();
        let (a, b) = (xs_u + on_alt.start, xs_u + on_alt.end);
        found.push(approx_repeat(&hx, a, b).map(|(lo, hi)| {
            let back = |p: usize| if p <= a { p as i64 } else if p >= b { p as i64 - d } else { a as i64 };
            (back(lo), back(hi))
        }));
    }
    for (lo, hi) in found.into_iter().flatten() {
        s = s.min(lo - 1);
        e = e.max(hi + 1);
    }
    (s.max(0), e.min(w.len() as i64))
}

/// Where `r` (a same-length change over every base a substitution `x` changes) reads
/// X's or REF's bases at X's changed bases: its other changes, one row per base; else
/// None.
fn split_beside(x: &Row, r: &Row) -> Option<Vec<Row>> {
    if x.changes_length() || r.changes_length() || r.pos > x.pos || r.end() < x.end() {
        return None;
    }
    let at = |p: i64| (p - r.pos) as usize;
    let (on_ref, _) = telling(x);
    let span = (x.pos + on_ref.start as i64)..(x.pos + on_ref.end as i64);
    let reads_x = span.clone().all(|p| r.alt_al[at(p)] == x.alt_al[(p - x.pos) as usize]);
    let reads_ref = span.clone().all(|p| r.alt_al[at(p)] == r.ref_al[at(p)]);
    if !(reads_x || reads_ref) {
        return None;
    }
    Some(
        (r.pos..r.end())
            .filter(|p| !span.contains(p) && r.alt_al[at(*p)] != r.ref_al[at(*p)])
            .map(|p| Row::new(p, &r.ref_al[at(p)..=at(p)], &r.alt_al[at(p)..=at(p)]))
            .collect(),
    )
}

/// The row turning `w` (starting at genome `start`) into `h`: their differing bases,
/// anchored on the base before them when one side is empty. None when they are equal.
fn row_between(w: &[u8], h: &[u8], start: i64) -> Option<Row> {
    let short = w.len().min(h.len());
    let pre = w.iter().zip(h).take_while(|(a, b)| a == b).count().min(short);
    if pre == short && w.len() == h.len() {
        return None;
    }
    let suf = w.iter().rev().zip(h.iter().rev()).take(short - pre).take_while(|(a, b)| a == b).count();
    let (mut s, we, he) = (pre, w.len() - suf, h.len() - suf);
    if s == we || s == he {
        s = s.checked_sub(1)?;
    }
    Some(Row::new(start + s as i64, &w[s..we], &h[s..he]))
}

/// Which of the guard's background (its recurring changes, any-base masks and recurring
/// indels beside the event) and of the found alleles the judge keeps as contexts
/// (`GBCMS_JUDGE_BACKGROUND`): `respell`, the default (decided 2026-10-27), all but the items
/// that respell the site's own candidate sequences ([`respellings`]), alone or joined into
/// one change at X's site, either way, for any state and any X, where no read needs them, the
/// masks kept; and, for measurement only, `asis` (the narrow rule it generalises: at a
/// substitution X, a recurring indel with which REF reads as X, left out where no read needs
/// it, the rest kept) and `off` (none of the background, the found alleles kept). Read when
/// a site is built.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum BackgroundRule {
    AsIs,
    Respell,
    Off,
}

pub(crate) fn background_rule() -> BackgroundRule {
    static RULE: std::sync::OnceLock<BackgroundRule> = std::sync::OnceLock::new();
    *RULE.get_or_init(|| {
        let v = std::env::var("GBCMS_JUDGE_BACKGROUND").unwrap_or_default();
        parse_background(&v).unwrap_or_else(|| {
            log::warn!("GBCMS_JUDGE_BACKGROUND={v:?} is not respell, asis or off; the default, respell, is used");
            BackgroundRule::Respell
        })
    })
}

fn parse_background(v: &str) -> Option<BackgroundRule> {
    match v {
        "respell" | "" => Some(BackgroundRule::Respell),
        "asis" => Some(BackgroundRule::AsIs),
        "off" => Some(BackgroundRule::Off),
        _ => None,
    }
}

/// Context rows nearest X's site, every set of them weighed on each side, when asking
/// whether a background item respells a candidate sequence ([`respellings`]).
const RESPELL_MAX_ROWS: usize = 5;

/// The stretch of the window a read over X's site can reach (genome `[lo, hi)`), spelled
/// with rows applied.
struct Stretch<'a> {
    reference: &'a [u8],
    start: i64,
    lo: i64,
    hi: i64,
}

impl Stretch<'_> {
    /// The stretch with `rows` applied (none meeting another, each inside it) and `wild`
    /// taking any base where no applied row covers it; None otherwise.
    fn spell(&self, rows: &[&Row], wild: &[i64]) -> Option<Vec<u8>> {
        let mut order: Vec<&Row> = rows.to_vec();
        order.sort_by_key(|r| (r.pos, r.end()));
        if order.iter().any(|r| r.pos < self.lo || r.end() > self.hi) || order.windows(2).any(|w| w[0].end() > w[1].pos) {
            return None;
        }
        let at = |p: i64| (p - self.start) as usize;
        let mut seq = self.reference[at(self.lo)..at(self.hi)].to_vec();
        for &m in wild {
            if (self.lo..self.hi).contains(&m) && !order.iter().any(|r| r.pos <= m && m < r.end()) {
                seq[(m - self.lo) as usize] = WILD;
            }
        }
        for r in order.iter().rev() {
            seq.splice((r.pos - self.lo) as usize..(r.end() - self.lo) as usize, r.alt_al.iter().copied());
        }
        Some(seq)
    }

    /// Genome position `p` in the stretch spelled with `rows`: past each row ending at or
    /// before it by the row's length change (one inside a row at the row's start).
    fn at(&self, rows: &[&Row], p: i64) -> i64 {
        let mut q = p - self.lo;
        for r in rows {
            if r.end() <= p {
                q += r.alt_al.len() as i64 - r.ref_al.len() as i64;
            } else if r.pos < p {
                q -= p - r.pos;
            }
        }
        q
    }
}

/// Spellings of the stretch: the rows applied to each, and its bases.
type Spellings<'r> = Vec<(Vec<&'r Row>, Vec<u8>)>;

/// How far two spellings agree from their left (or, `right`, their right) end, a
/// [`WILD`] base on either side matching any.
fn agree(a: &[u8], b: &[u8], right: bool) -> usize {
    let eq = |(x, y): (&u8, &u8)| x == y || *x == WILD || *y == WILD;
    if right {
        a.iter().rev().zip(b.iter().rev()).take_while(|&p| eq(p)).count()
    } else {
        a.iter().zip(b.iter()).take_while(|&p| eq(p)).count()
    }
}

/// Which of the guard's background items respell the site's own candidate sequences:
/// per change row in `rows` and per any-base mask in `masks` (genome positions), whether
/// with it one of REF and another candidate sequence (X, or another allele at X's site,
/// `states` with X first) spells the other over a read reaching X's site, where without
/// it the two differ there: 2 where it turns one into the other whole (a recurring
/// indel with which X spells REF's haplotype in a repeat), 1 where only over a read (a
/// mask on a duplication's one telling base; also a real change beside a repeat, which
/// only the reads can tell apart: [`Site::build_backed`] asks them), 0 where not. Such an
/// item is the aligner's respelling of the row's allele, not a change beside the event:
/// as a context it would tie the reads that tell the candidates apart.
///
/// With `one_way`, only REF with the item spelling a candidate counts (not a candidate with
/// it spelling REF); [`respell_directions`] weighs either direction alone.
///
/// Reads are taken as placed, of length `l`: one reading in from the left flank tells two
/// spellings apart once it passes the first base where they differ, so an item that moves
/// that base on, where it lies at X's site or within a read of it, unmakes the telling
/// (and from the right flank alike). Contexts are the other rows (`others`, the
/// [`RESPELL_MAX_ROWS`] nearest X's site), every set of them that can lie together, on
/// either side, so a read's other differences are explained as the judge would.
#[allow(clippy::too_many_arguments)]
fn respellings<'r>(
    reference: &[u8],
    start: i64,
    site: (i64, i64),
    l: i64,
    states: &'r [Row],
    others: &'r [Row],
    rows: &'r [Row],
    masks: &[i64],
) -> (Vec<u8>, Vec<u8>) {
    respell_levels(reference, start, site, l, states, others, rows, masks, false)
}

/// [`respellings`], with `one_way` as it says.
#[allow(clippy::too_many_arguments)]
fn respell_levels<'r>(
    reference: &[u8],
    start: i64,
    site: (i64, i64),
    l: i64,
    states: &'r [Row],
    others: &'r [Row],
    rows: &'r [Row],
    masks: &[i64],
    one_way: bool,
) -> (Vec<u8>, Vec<u8>) {
    respell_directions(reference, start, site, l, states, others, rows, masks, (true, !one_way))
}

/// Measurement only: [`respell_levels`] weighing REF with the item spelling a candidate
/// (`dirs.0`) and a candidate with it spelling REF (`dirs.1`), each where asked.
#[allow(clippy::too_many_arguments)]
fn respell_directions<'r>(
    reference: &[u8],
    start: i64,
    site: (i64, i64),
    l: i64,
    states: &'r [Row],
    others: &'r [Row],
    rows: &'r [Row],
    masks: &[i64],
    dirs: (bool, bool),
) -> (Vec<u8>, Vec<u8>) {
    let st = Stretch { reference, start, lo: (site.0 - l).max(start), hi: (site.1 + l).min(start + reference.len() as i64) };
    let inside = |r: &Row| r.pos >= st.lo && r.end() <= st.hi;
    let meet = |a: &Row, b: &Row| a.pos < b.end() && b.pos < a.end();
    let near = |r: &Row| (r.pos - site.1).max(site.0 - r.end()).max(0);
    let (mut out_rows, mut out_masks) = (vec![0u8; rows.len()], vec![0u8; masks.len()]);
    if !states.iter().all(&inside) {
        return (out_rows, out_masks);
    }
    // The context pool: the other rows and the background, each once, nearest first.
    let mut all: Vec<&'r Row> = Vec::new();
    for r in others.iter().chain(rows).filter(|r| inside(r)) {
        if !all.contains(&r) {
            all.push(r);
        }
    }
    all.sort_by_key(|r| (near(r), r.pos));
    // Every set of the pool (the item weighed left out) that can lie together.
    let sets_without = |item: Option<&Row>| -> Vec<Vec<&'r Row>> {
        let pool: Vec<&'r Row> = all.iter().copied().filter(|r| item.is_none_or(|b| *r != b)).take(RESPELL_MAX_ROWS).collect();
        (0u32..1 << pool.len())
            .map(|m| (0..pool.len()).filter(|i| m >> i & 1 == 1).map(|i| pool[i]).collect::<Vec<_>>())
            .filter(|c| c.iter().enumerate().all(|(i, a)| c[i + 1..].iter().all(|b| !meet(a, b))))
            .collect()
    };
    // A side's spellings: (its rows, the stretch spelled) for each set it can lie with.
    let side = |sets: &[Vec<&'r Row>], state: Option<&'r Row>, item: Option<&'r Row>, wild: &[i64]| -> Spellings<'r> {
        sets.iter()
            .filter_map(|c| {
                let rows: Vec<&Row> = state.into_iter().chain(item).chain(c.iter().copied()).collect();
                st.spell(&rows, wild).map(|h| (rows, h))
            })
            .collect()
    };
    // How far `with` reads into one of `target`'s spellings past where `without` does,
    // from a first difference at X's site or within a read of it (`target_m`: the target
    // spelled with the item, where the item is a mask on every spelling): 2 where a `with`
    // spelling is the target whole (the item turns one candidate into the other), 1
    // where it reads further over a read only, 0 where it does not.
    let unmakes = |target: &Spellings, target_m: Option<&Spellings>, without: &Spellings, with: &Spellings| -> u8 {
        let mut level = 0u8;
        for (k, (t_rows, t)) in target.iter().enumerate() {
            let tm = target_m.map_or(t, |m| &m[k].1);
            let (ts, te) = (st.at(t_rows, site.0), st.at(t_rows, site.1));
            let n = t.len() as i64;
            for right in [false, true] {
                let (Some(before), Some(after)) = (
                    without.iter().map(|(_, h)| agree(h, t, right) as i64).max(),
                    with.iter().map(|(_, h)| agree(h, tm, right) as i64).max(),
                ) else {
                    continue;
                };
                // the first difference reading in from that side, as a position from it
                let (s0, e0) = if right { (n - te, n - ts) } else { (ts, te) };
                if after > before && before + 1 >= s0 && before < e0 + l - 1 {
                    let whole = with.iter().any(|(_, h)| h.len() == tm.len() && agree(h, tm, false) == h.len());
                    level = level.max(if whole { 2 } else { 1 });
                }
            }
        }
        level
    };
    for state in states {
        for (i, b) in rows.iter().enumerate() {
            if out_rows[i] == 2 || !inside(b) {
                continue;
            }
            let sets = sets_without(Some(b));
            let (refs, alts) = (side(&sets, None, None, &[]), side(&sets, Some(state), None, &[]));
            // REF with the item spells the candidate, or the candidate with it spells REF.
            let with_b: Vec<Vec<&'r Row>> = sets.iter().filter(|c| c.iter().all(|r| !meet(r, b))).cloned().collect();
            let ref_b = side(&with_b, None, Some(b), &[]);
            let alt_b = if meet(b, state) { Vec::new() } else { side(&with_b, Some(state), Some(b), &[]) };
            let back = if dirs.1 { unmakes(&refs, None, &alts, &alt_b) } else { 0 };
            let forth = if dirs.0 { unmakes(&alts, None, &refs, &ref_b) } else { 0 };
            out_rows[i] = out_rows[i].max(forth).max(back);
        }
        let sets = sets_without(None);
        let (refs, alts) = (side(&sets, None, None, &[]), side(&sets, Some(state), None, &[]));
        for (i, &m) in masks.iter().enumerate() {
            if out_masks[i] == 2 || !(st.lo..st.hi).contains(&m) {
                continue;
            }
            // the mask on every spelling: REF's against the candidate's, both masked
            let (refs_m, alts_m) = (side(&sets, None, None, &[m]), side(&sets, Some(state), None, &[m]));
            if alts_m.len() == alts.len() {
                out_masks[i] = out_masks[i].max(unmakes(&alts, Some(&alts_m), &refs, &refs_m));
            }
        }
    }
    (out_rows, out_masks)
}

/// The context rows a site is built without (the respelling check's).
#[derive(Debug, Clone, Default)]
struct Dropped {
    rows: Vec<Row>,
}

/// Whether a read's aligned span reaches `r`.
fn reaches(read: &ReadView, r: &Row) -> bool {
    read.ref_start < r.end() && r.pos < read.ref_end
}

/// A row's type as a trace names it: SNV, MNV, INS, DEL or DELINS.
fn x_type(r: &Row) -> &'static str {
    let (on_ref, on_alt) = telling(r);
    match (r.ref_al.len() == r.alt_al.len(), on_ref.is_empty(), on_alt.is_empty()) {
        (true, _, _) if on_ref.len() <= 1 => "SNV",
        (true, _, _) => "MNV",
        (false, true, _) => "INS",
        (false, _, true) => "DEL",
        _ => "DELINS",
    }
}

/// A row as a trace names it.
fn named(r: &Row) -> String {
    format!("{}:{}>{}", r.pos + 1, String::from_utf8_lossy(&r.ref_al), String::from_utf8_lossy(&r.alt_al))
}

/// Reads weighed when asking whether the reads need an item ([`needed_by`]), the first
/// reaching it.
const ATTEST_MAX_READS: usize = 400;

/// Whether the reads reaching `item` need it: on [`FOUND_MIN_READS`] of them and
/// [`GUARD_REF_SHARE`], a read the site's candidate sequences call REF or ALT with it
/// (`kept`) gets another call without it (`without`). An item no read needs for its call
/// (a germline deletion beside the row does: its reads keep their call only with it)
/// explains none of their differences.
fn needed_by(kept: &Site, without: &Site, reads: &[ReadView], item: &Row) -> bool {
    let near: Vec<&ReadView> = reads.iter().filter(|r| reaches(r, item)).take(ATTEST_MAX_READS).collect();
    let mut ws = Workspace::default();
    let mut call = |site: &Site, read: &ReadView| {
        let gr = site.graph_read(read);
        let mut best: [Option<i32>; 3] = [None; 3];
        for g in &site.graphs {
            let s = graph::decided(&ShortReads, &gr, &g.graph, &mut ws).0;
            let b = &mut best[g.class as usize];
            *b = Some(b.map_or(s, |v| v.max(s)));
        }
        super::decide::decide(best)
    };
    let need = near
        .iter()
        .filter(|r| {
            let with = call(kept, r);
            matches!(with, super::decide::Call::Alt | super::decide::Call::Ref) && call(without, r) != with
        })
        .count();
    need as u32 >= FOUND_MIN_READS && need as f64 >= GUARD_REF_SHARE * near.len() as f64
}

/// Context items (the guard's background, the found contexts) weighed together when
/// asking whether a set of them spells one change at X's site ([`joined_at_site`]).
const JOINED_MAX_ROWS: usize = 6;

/// Per item of `items` (context rows none of which meets X alone), whether it belongs to a
/// set of them that, applied to REF together, spells one change (the minimal row between
/// REF and that spelling) that `meets_x`, the set's every item needed for it: a carrier's
/// event the aligner wrote as two indels on either side of X (chr16:23024630 TA>AC, HG002:
/// a 7 bp deletion and an inserted C spelling X with a 6 bp deletion beside it). Such a
/// set is an allele at X's site, not context: with REF it would spell the carrier.
/// Weighed over the [`JOINED_MAX_ROWS`] items nearest X's `site` within a read (`l`) of it.
fn joined_at_site(reference: &[u8], start: i64, site: (i64, i64), l: i64, items: &[Row], meets_x: impl Fn(&Row) -> bool) -> Vec<bool> {
    let mut out = vec![false; items.len()];
    let near = |r: &Row| (r.pos - site.1).max(site.0 - r.end()).max(0);
    let mut idx: Vec<usize> = (0..items.len()).filter(|&i| near(&items[i]) <= l).collect();
    idx.sort_by_key(|&i| (near(&items[i]), items[i].pos));
    idx.dedup_by(|a, b| items[*a] == items[*b]);
    idx.truncate(JOINED_MAX_ROWS);
    if idx.len() < 2 {
        return out;
    }
    let whole = Stretch { reference, start, lo: start, hi: start + reference.len() as i64 };
    let n = idx.len();
    // Per set (bitmask over `idx`), whether REF with it spells one change meeting X.
    let mut at_site = vec![false; 1 << n];
    for m in 1usize..1 << n {
        let set: Vec<&Row> = (0..n).filter(|k| m >> k & 1 == 1).map(|k| &items[idx[k]]).collect();
        if set.len() < 2 {
            continue;
        }
        let Some(spelled) = whole.spell(&set, &[]) else { continue };
        at_site[m] = row_between(reference, &spelled, start).is_some_and(|r| meets_x(&r));
        // every item needed: no set of one fewer does it already
        let minimal = (0..n).filter(|k| m >> k & 1 == 1).all(|k| !at_site[m & !(1 << k)]);
        if at_site[m] && minimal {
            for k in (0..n).filter(|k| m >> k & 1 == 1) {
                out[idx[k]] = true;
            }
        }
    }
    out
}

/// The bases a row's alleles differ at, as offsets into REF and into ALT: past the
/// bases they share at the start, then at the end. Empty on one side for a pure
/// insertion or deletion.
fn telling(r: &Row) -> (std::ops::Range<usize>, std::ops::Range<usize>) {
    let (a, b) = (&r.ref_al, &r.alt_al);
    let short = a.len().min(b.len());
    let p = a.iter().zip(b.iter()).take_while(|(x, y)| x == y).count().min(short);
    let s = a.iter().rev().zip(b.iter().rev()).take_while(|(x, y)| x == y).count().min(short - p);
    (p..a.len() - s, p..b.len() - s)
}

/// The bases a substitution or delins `x` keeps among those it changes: pairs of offsets
/// (into REF, into ALT) inside its telling bases where its two alleles read alike, aligned
/// from the left and, for a delins, from the right. A change there that the REF molecules
/// carry too is a separate event on both haplotypes (RJ-25), never another allele. Empty
/// for a pure insertion or deletion.
fn kept_bases(x: &Row) -> Vec<(usize, usize)> {
    if pure_indel(x) {
        return Vec::new();
    }
    let (on_ref, on_alt) = telling(x);
    let n = on_ref.len().min(on_alt.len());
    let mut out: Vec<(usize, usize)> = (0..n)
        .map(|j| (on_ref.start + j, on_alt.start + j))
        .filter(|&(i, k)| x.ref_al[i] == x.alt_al[k])
        .collect();
    if x.changes_length() {
        for j in 0..n {
            let (i, k) = (on_ref.end - 1 - j, on_alt.end - 1 - j);
            if x.ref_al[i] == x.alt_al[k] && !out.iter().any(|&(a, _)| a == i) {
                out.push((i, k));
            }
        }
    }
    out
}

/// Whether `r` changes the length alone: its alleles, the bases they share at either end set
/// aside, one empty and the other not. A pure insertion or deletion anchored on the base before
/// it, or, at a contig's first base, on the base after it (the base-after form, which
/// [`pure_indel`] calls complex).
fn length_change_only(r: &Row) -> bool {
    let (on_ref, on_alt) = telling(r);
    on_ref.is_empty() != on_alt.is_empty()
}

fn pure_indel(r: &Row) -> bool {
    super::super::window::is_pure_indel(&String::from_utf8_lossy(&r.ref_al), &String::from_utf8_lossy(&r.alt_al))
}

impl Site {
    fn is_state(&self, r: usize) -> bool {
        r >= 1 && r <= self.n_states
    }

    /// Whether two rows are at one place: never on one haplotype.
    fn conflict(&self, a: usize, b: usize) -> bool {
        let meet = |(s, e): (i64, i64), (t, f): (i64, i64)| s < f && t < e;
        if self.rows[a].changes_length() && self.rows[b].changes_length() {
            meet(self.places[a], self.places[b])
        } else {
            meet(self.spans[a], self.spans[b])
        }
    }

    /// Whether X substitutes bases without changing the length (an SNV or an MNP).
    pub(crate) fn is_substitution(&self) -> bool {
        !self.rows[0].changes_length()
    }

    /// X's anchor (its first REF base) in `hap`: X's own first base where X is applied.
    pub(crate) fn anchor_in(&self, hap: &Hap) -> i64 {
        hap.x_seg.map_or_else(|| self.column_in(hap, self.spans[0].0), |(hs, _)| hs)
    }

    /// The column of window position `p` in `hap`: inside an applied row that keeps the
    /// length (a state holding X's change), the base at its offset there; else as
    /// [`Hap::to_hap`] maps it.
    fn column_in(&self, hap: &Hap, p: i64) -> i64 {
        hap.segs
            .iter()
            .find(|&&(rs, re, hs, he)| rs <= p && p < re && he - hs == re - rs)
            .map_or_else(|| hap.to_hap(p, false), |&(rs, _, hs, _)| hs + (p - rs))
    }

    /// The haplotype positions, with their bases, where `hap` tells X from REF at X's
    /// site: X's own differing bases (ALT) or REF's (REF). None for a pure insertion or
    /// deletion (its tract is read by the scores) and for another allele's haplotype. A
    /// base a substitution keeps tells nothing (REF's and X's read alike there): a change
    /// there is a separate event or another allele, the scores weighing it on both sides.
    pub(crate) fn telling_in(&self, hap: &Hap) -> Vec<(i64, u8)> {
        let x = &self.rows[0];
        if x.changes_length() && super::super::window::is_pure_indel(
            &String::from_utf8_lossy(&x.ref_al),
            &String::from_utf8_lossy(&x.alt_al),
        ) {
            return Vec::new();
        }
        let (on_ref, on_alt) = telling(x);
        let same = !x.changes_length();
        let on_ref = on_ref.filter(|&k| !same || x.ref_al[k] != x.alt_al[k]);
        let on_alt = on_alt.filter(|&k| !same || x.ref_al[k] != x.alt_al[k]);
        match (hap.class, hap.x_seg) {
            (Class::Alt, Some((hs, _))) => on_alt.map(|k| (hs + k as i64, x.alt_al[k])).collect(),
            // a state holding X's change: X's bases at their places in it
            (Class::Alt, None) if !x.changes_length() => {
                on_alt.map(|k| (self.column_in(hap, self.spans[0].0 + k as i64), x.alt_al[k])).collect()
            }
            (Class::Ref, _) => on_ref.map(|k| (hap.to_hap(self.spans[0].0 + k as i64, false), x.ref_al[k])).collect(),
            _ => Vec::new(),
        }
    }

    /// For a substitution X, the haplotype positions of X's telling bases in `hap` (a
    /// haplotype carrying X or holding REF's bases there), each with REF's base and X's
    /// there; empty for a length change and on another allele's haplotype.
    pub(crate) fn telling_pairs(&self, hap: &Hap) -> Vec<(i64, u8, u8)> {
        let x = &self.rows[0];
        if x.changes_length() || hap.class == Class::Other {
            return Vec::new();
        }
        let (on_ref, _) = telling(x);
        on_ref
            .map(|k| {
                let p = self.spans[0].0 + k as i64;
                let col = match (hap.class, hap.x_seg) {
                    (Class::Alt, Some((hs, _))) => hs + k as i64,
                    (Class::Alt, None) => self.column_in(hap, p),
                    _ => hap.to_hap(p, false),
                };
                (col, x.ref_al[k], x.alt_al[k])
            })
            .collect()
    }

    /// For a delins X, the haplotype positions of the telling bases of the allele `hap`
    /// carries (X's on a haplotype carrying X, REF's on REF's), each with that allele's base
    /// and the other allele's bases at the same telling index counted from the left and from
    /// the right (where a read's own placement of the length change may put them). Bases X
    /// keeps (its two alleles alike there from one end) are left out: they tell nothing.
    /// Empty for a substitution, a pure indel and another allele's haplotype.
    pub(crate) fn telling_with_other(&self, hap: &Hap) -> Vec<(i64, u8, [Option<u8>; 2])> {
        let x = &self.rows[0];
        if !x.changes_length() || pure_indel(x) {
            return Vec::new();
        }
        let (on_ref, on_alt) = telling(x);
        let kept = kept_bases(x);
        // the other allele's bases at telling index `k` of a range `own`, over `other`
        let at = |own: &std::ops::Range<usize>, k: usize, other: &std::ops::Range<usize>, bases: &[u8]| -> [Option<u8>; 2] {
            let (j, from_end) = (k - own.start, own.end - 1 - k);
            [(j < other.len()).then(|| bases[other.start + j]), (from_end < other.len()).then(|| bases[other.end - 1 - from_end])]
        };
        match (hap.class, hap.x_seg) {
            (Class::Alt, Some((hs, _))) => on_alt
                .clone()
                .filter(|&k| !kept.iter().any(|&(_, a)| a == k))
                .map(|k| (hs + k as i64, x.alt_al[k], at(&on_alt, k, &on_ref, &x.ref_al)))
                .collect(),
            (Class::Ref, _) => on_ref
                .clone()
                .filter(|&i| !kept.iter().any(|&(r, _)| r == i))
                .map(|i| (hap.to_hap(self.spans[0].0 + i as i64, false), x.ref_al[i], at(&on_ref, i, &on_alt, &x.alt_al)))
                .collect(),
            _ => Vec::new(),
        }
    }

    /// X's first base, 0-based.
    pub(crate) fn x_pos(&self) -> i64 {
        self.rows[0].pos
    }

    /// The read's alignment to state `k`'s haplotype `hap` (the one its best path spells,
    /// [`Site::hap_on_path`]) read off that path (`steps`, [`graph::kept_steps`]): each read
    /// base's graph base mapped to its column in `hap` (a reference base by its coordinate past
    /// the rows applied before it, an allele's by its place in the row's span on `hap`);
    /// columns are `hap`'s own (the trim starts at 0). The read's call is examined on it: the
    /// alignment its score came from, not a second one.
    pub(crate) fn alignment_on_path(&self, k: usize, steps: &graph::PathSteps, hap: &Hap, read: &[u8]) -> (Alignment, usize) {
        let g = &self.graphs[k].graph;
        let column = |v: usize| -> usize {
            let seg = g.segment_of(v);
            let c = match seg.row {
                None => hap.to_hap(g.coord[v], false),
                Some(r) => {
                    let (s, e) = self.spans[r];
                    let hs = hap.segs.iter().find(|q| q.0 == s && q.1 == e).map_or(0, |q| q.2);
                    hs + g.offset_in_segment(v) as i64
                }
            };
            debug_assert_eq!(hap.base(c as usize), g.bases[v], "a path's base is its haplotype's");
            c as usize
        };
        let placed: Vec<Placed> = steps
            .steps
            .iter()
            .map(|s| match *s {
                graph::Step::At(v) => Placed::At(column(v as usize)),
                graph::Step::Inserted => Placed::Inserted,
                graph::Step::LeftOver => Placed::LeftOver,
            })
            .collect();
        let costs = (0..read.len())
            .map(|i| match placed[i] {
                Placed::At(c) => score::substitution_cost(read[i], hap.base(c)),
                Placed::Inserted if i > 0 && placed[i - 1] == Placed::Inserted => -score::GAP_EXTEND,
                Placed::Inserted => -(score::GAP_OPEN + score::GAP_EXTEND),
                Placed::LeftOver => score::overhang_cost(read[i]),
            })
            .collect();
        (Alignment { score: steps.score, placed, skipped_after: steps.skipped_after.clone(), costs }, 0)
    }

    /// X's size, as the class bins measure it: the longer of its two alleles once the bases
    /// they share at either end are set aside (its telling bases on each allele).
    pub(crate) fn x_size(&self) -> usize {
        let (on_ref, on_alt) = telling(&self.rows[0]);
        on_ref.len().max(on_alt.len())
    }

    /// How many reads besides the one asking read `base` at genome position `pos` by their
    /// own alignments (at a substitution of two bases or more or a delins X).
    pub(crate) fn recurrence(&self, pos: i64, base: u8) -> u32 {
        self.tally.get(&(pos, base)).map_or(0, |n| n.saturating_sub(1))
    }

    /// Whether the guard's background explains `base` at genome position `pos`: one of its
    /// recurring changes there, or any base where it masks.
    pub(crate) fn background_explains(&self, pos: i64, base: u8) -> bool {
        self.background_bases.contains(&(pos, base)) || self.background_masks.binary_search(&pos).is_ok()
    }

    /// The bases X changes (REF's and X's differ), for a substitution X; 0 for a length
    /// change.
    pub(crate) fn changed_bases(&self) -> usize {
        let x = &self.rows[0];
        if x.changes_length() {
            return 0;
        }
        x.ref_al.iter().zip(&x.alt_al).filter(|(a, b)| a != b).count()
    }

    /// The states holding X's change (`holds_x`), in order.
    pub(crate) fn holders(&self) -> Vec<&Row> {
        (1..=self.n_states).filter(|&r| self.holds_x[r]).map(|r| &self.rows[r]).collect()
    }

    /// Which of [`Site::holders`] (a bit each, in order, the first 32) `read` spells whole
    /// on its alignment `al` (columns offset by `off`) to `hap`: every base such a state
    /// changes, X's own included, read confidently (at or above `--min-baseq`) as the
    /// state's. An ALT call does not vouch for X's bases: an MNP read with a changed base
    /// below the floor counts ALT, and spells no longer MNV there.
    pub(crate) fn spelled_holders(&self, read: &ReadView, hap: &Hap, al: &Alignment, off: usize) -> u32 {
        if !self.holds_x.contains(&true) {
            return 0;
        }
        let at: HashMap<i64, u8> = al
            .placed
            .iter()
            .enumerate()
            .filter_map(|(i, p)| if let Placed::At(c) = *p { Some(((c + off) as i64, read.seq[i])) } else { None })
            .collect();
        let spells = |r: usize| {
            let row = &self.rows[r];
            if row.changes_length() {
                // X with a change at a base it keeps (a delins): no listed MNV
                return false;
            }
            let (s, _) = self.spans[r];
            (0..row.ref_al.len()).filter(|&i| row.ref_al[i] != row.alt_al[i]).all(|i| {
                let col = self.column_in(hap, s + i as i64);
                at.get(&col).is_some_and(|&b| b != WILD && b == row.alt_al[i])
            })
        };
        (1..=self.n_states)
            .filter(|&r| self.holds_x[r])
            .take(32)
            .enumerate()
            .filter(|&(_, r)| spells(r))
            .fold(0u32, |m, (k, _)| m | 1 << k)
    }

    /// X's site in `hap`'s columns, `[a, b)`: the bases X changes and, for a length
    /// change, every place it can sit and the repeat it lengthens or shortens, with a base
    /// each side.
    /// X's site in genome coordinates, `[s, e)`.
    pub(crate) fn site_span(&self) -> (i64, i64) {
        (self.site_window.0 + self.start, self.site_window.1 + self.start)
    }

    /// The read's best fit anywhere on the window's reference (every row's REF, no
    /// context), whole read, ends free: its fit to the local reference neighbourhood.
    pub(crate) fn neighbourhood_fit(&self, read: &ReadView) -> i32 {
        super::score::full_score(&read.seq, &self.reference)
    }

    /// The columns of `hap` inside a long event's own bases, more than [`EDGE`] from
    /// either end: X's inserted (or substituted) bases on a haplotype carrying X, its
    /// deleted (or substituted) REF bases on one holding REF's, where they number
    /// [`LONG_EVENT`] or more (an event a read is not expected to hold whole, told from
    /// REF at its junctions). None for a shorter event, and on another allele's haplotype.
    fn long_interior(&self, hap: &Hap) -> Option<(i64, i64)> {
        let (on_ref, on_alt) = telling(&self.rows[0]);
        let (c0, c1) = match (hap.class, hap.x_seg) {
            (Class::Alt, Some((hs, _))) => (hs + on_alt.start as i64, hs + on_alt.end as i64),
            (Class::Ref, _) => {
                let xs = self.spans[0].0;
                (hap.to_hap(xs + on_ref.start as i64, false), hap.to_hap(xs + on_ref.end as i64, true))
            }
            _ => return None,
        };
        (c1 - c0 >= LONG_EVENT).then_some((c0 + EDGE, c1 - EDGE))
    }

    fn site_columns(&self, hap: &Hap) -> (i64, i64) {
        let (s, e) = self.site_window;
        (hap.to_hap(s, false), hap.to_hap(e, true))
    }

    /// Whether an alignment of a read to `hap` (its columns offset by `off`) inserts or
    /// deletes bases at X's site: at any place the indel can slide to (an indel in a
    /// repeat is one allele wherever it is written in it).
    pub(crate) fn gap_at_site(&self, hap: &Hap, read: &[u8], al: &Alignment, off: usize) -> bool {
        // a gap lying wholly inside a long event's own bases is a change on that molecule,
        // never another allele of X
        let interior = if other_option().edges { self.long_interior(hap) } else { None };
        let inside = |lo: i64, hi: i64| interior.is_some_and(|(i0, i1)| i0 <= lo && hi <= i1);
        gap_between(hap, read, al, off, self.site_columns(hap), inside, false)
    }

    /// Measurement only (`GBCMS_JUDGE_OWNGAP=all`): whether an alignment of a read to the
    /// state `hap` is (its columns offset by `off`) inserts or deletes bases of its own at
    /// that state's change, wherever the indel can slide to, read bases inserted before its
    /// first aligned base included. REF's and X's change is X's site (the bases a
    /// substitution changes; a length change's run and a base each side, the discrimination
    /// window); another allele's is the bases it changes, and for one changing the length a
    /// base each side too ([`Site::gap_at_state`]). An indel beside a substitution's bases is
    /// a separate event (RJ-8), never a reason the read does not carry it. A gap wholly
    /// inside a long event's own bases is a change on that molecule, as the fit check has it.
    pub(crate) fn own_gap_at(&self, hap: &Hap, read: &[u8], al: &Alignment, off: usize) -> bool {
        if hap.class == Class::Other || (hap.class == Class::Alt && hap.x_seg.is_none()) {
            let Some(&r) = hap.rows.iter().find(|&&r| self.is_state(r)) else {
                return false;
            };
            if self.rows[r].changes_length() {
                return self.gap_at_state(hap, read, al, off);
            }
            let Some(&(_, _, hs, he)) = hap.segs.iter().find(|&&(rs, re, _, _)| (rs, re) == self.spans[r]) else {
                return false;
            };
            return gap_between(hap, read, al, off, (hs, he), |_, _| false, true);
        }
        let interior = if other_option().edges { self.long_interior(hap) } else { None };
        let inside = |lo: i64, hi: i64| interior.is_some_and(|(i0, i1)| i0 <= lo && hi <= i1);
        gap_between(hap, read, al, off, self.site_columns(hap), inside, true)
    }

    /// Whether an alignment of a read to another allele's haplotype `hap` (its columns
    /// offset by `off`) inserts or deletes bases of its own where that allele changes the
    /// window, a base each side, wherever the indel can slide to, read bases inserted before
    /// its first aligned base included: the read writes another indel across the allele's
    /// change rather than carrying it.
    pub(crate) fn gap_at_state(&self, hap: &Hap, read: &[u8], al: &Alignment, off: usize) -> bool {
        let Some(&r) = hap.rows.iter().find(|&&r| self.is_state(r)) else {
            return false;
        };
        let Some(&(_, _, hs, he)) = hap.segs.iter().find(|&&(rs, re, _, _)| (rs, re) == self.spans[r]) else {
            return false;
        };
        gap_between(hap, read, al, off, (hs - 1, he + 1), |_, _| false, true)
    }
}

/// Whether an alignment of `read` to `hap` (its columns offset by `off`) inserts or deletes
/// bases meeting the columns `[a, b)`, at any place the indel can slide to, but for one lying
/// wholly where `inside` says; with `leading`, read bases inserted before the alignment's
/// first column too.
#[allow(clippy::too_many_arguments)]
fn gap_between(
    hap: &Hap,
    read: &[u8],
    al: &Alignment,
    off: usize,
    (a, b): (i64, i64),
    inside: impl Fn(i64, i64) -> bool,
    leading: bool,
) -> bool {
    let h = &hap.seq;
    let m = h.len() as i64;
    let at = |c: i64| (0..m).contains(&c).then(|| h[c as usize]);
    let mut last: Option<i64> = None;
    let mut i = 0;
    while i < al.placed.len() {
        match al.placed[i] {
            Placed::At(c) => {
                let c = (c + off) as i64;
                let k = i64::from(al.skipped_after[i]);
                if k > 0 {
                    // hap columns [lo, lo + k) skipped, slid each way
                    let (mut lo, mut hi) = (c + 1, c + 1);
                    while at(lo - 1).is_some() && at(lo - 1) == at(lo - 1 + k) {
                        lo -= 1;
                    }
                    while at(hi + k).is_some() && at(hi) == at(hi + k) {
                        hi += 1;
                    }
                    if lo < b && hi + k > a && !inside(lo, hi + k) {
                        return true;
                    }
                }
                last = Some(c);
                i += 1;
            }
            Placed::Inserted => {
                // the read bases [i, j) inserted before column `last` + 1
                let j = (i..al.placed.len()).find(|&t| al.placed[t] != Placed::Inserted).unwrap_or(al.placed.len());
                let first = || match al.placed.get(j) {
                    Some(&Placed::At(c)) if leading => Some((c + off) as i64 - 1),
                    _ => None,
                };
                if let Some(c) = last.or_else(first) {
                    let ins = &read[i..j];
                    let n = ins.len();
                    // slide left: the base before the boundary equals the insert's last
                    // base (rotated as it slides); right: the base after it its first
                    let mut lo = c + 1;
                    let mut rot = n;
                    while lo > 0 && at(lo - 1) == Some(ins[(rot - 1) % n]) {
                        lo -= 1;
                        rot = if rot == 1 { n } else { rot - 1 };
                    }
                    let mut hi = c + 1;
                    let mut r = 0;
                    while at(hi) == Some(ins[r % n]) {
                        hi += 1;
                        r += 1;
                    }
                    if a < hi && lo < b && !inside(lo, hi) {
                        return true;
                    }
                }
                i = j;
            }
            Placed::LeftOver => i += 1,
        }
    }
    false
}

impl Site {
    /// Whether `read` reaches where another state's haplotype `other` differs from REF's:
    /// on REF's haplotype with `other`'s contexts (its rows but the state's), its aligned
    /// bases reach past their common prefix and start before their common suffix. In a
    /// repeat the two differ only past its end (reading in from either side), so a read
    /// ending inside it reaches nothing that tells them apart.
    pub(crate) fn reaches_state(&self, read: &ReadView, other: &Hap) -> bool {
        let rows: Vec<usize> = other.rows.iter().copied().filter(|&r| !self.is_state(r)).collect();
        let reference = Hap::build(self, &rows);
        let (a, b) = (&reference.seq, &other.seq);
        if a == b {
            return false;
        }
        let prefix = a.iter().zip(b).take_while(|(x, y)| x == y).count();
        let suffix = a.iter().rev().zip(b.iter().rev()).take_while(|(x, y)| x == y).count();
        let (al, off) = Scorer::new(read).align(self, &reference);
        al.covers().is_some_and(|(f, l)| l + off >= prefix && f + off + suffix < a.len())
    }

    /// The haplotypes of `target` over the window `(start, reference)` (upper case),
    /// from `candidates` and the variant's `reads`. None when the window does not hold
    /// the target's REF.
    pub(crate) fn build(target: Row, candidates: &Candidates, window: (i64, Vec<u8>), reads: &[ReadView]) -> Option<Site> {
        Site::build_backed(target, candidates, window, reads, other_option().rule, background_rule())
    }

    /// [`Site::build`] with the found alleles kept as states by `rule` and the guard's
    /// background kept by `bg_rule`.
    pub(crate) fn build_backed(
        target: Row,
        candidates: &Candidates,
        window: (i64, Vec<u8>),
        reads: &[ReadView],
        rule: FoundRule,
        bg_rule: BackgroundRule,
    ) -> Option<Site> {
        Site::build_core(target, candidates, window, reads, rule, bg_rule, None)
    }

    /// [`Site::build_backed`], the context items in `dropped` left out where given (the
    /// respelling check's own builds).
    #[allow(clippy::too_many_arguments)]
    fn build_core(
        target: Row,
        candidates: &Candidates,
        (start, reference): (i64, Vec<u8>),
        reads: &[ReadView],
        rule: FoundRule,
        bg_rule: BackgroundRule,
        dropped: Option<&Dropped>,
    ) -> Option<Site> {
        let wn = reference.len() as i64;
        let holds = |r: &Row| {
            let (s, e) = (r.pos - start, r.end() - start);
            s >= 0 && e <= wn && reference[s as usize..e as usize] == r.ref_al[..] && r.ref_al != r.alt_al
        };
        if !holds(&target) || target.ref_al.is_empty() {
            return None;
        }
        let base = |g: i64| usize::try_from(g - start).ok().and_then(|i| reference.get(i).copied());
        let pure = |r: &Row| {
            super::super::window::is_pure_indel(&String::from_utf8_lossy(&r.ref_al), &String::from_utf8_lossy(&r.alt_al))
        };
        // A pure length change's run: the boundaries an insertion slides over, or the
        // bases a deletion slides through (genome, [lo, hi)).
        let run = |r: &Row| -> (i64, i64) {
            let (ref_s, alt_s) = (String::from_utf8_lossy(&r.ref_al), String::from_utf8_lossy(&r.alt_al));
            super::super::window::shift_region_over(r.pos, &ref_s, &alt_s, base)
        };
        // A length change's equivalent placements: from the anchor before the first
        // place it can sit to the end of the last (the slide runs over the window).
        let place = |r: &Row| -> (i64, i64) {
            if !r.changes_length() {
                return (r.pos, r.end());
            }
            let (lo, hi) = run(r);
            ((lo - 1).min(r.pos), hi.max(r.end()))
        };
        let meets = |a: &Row, b: &Row| {
            let m = |(s, e): (i64, i64), (t, f): (i64, i64)| s < f && t < e;
            if a.changes_length() && b.changes_length() {
                m(place(a), place(b))
            } else {
                m((a.pos, a.end()), (b.pos, b.end()))
            }
        };
        // The bases a long pure deletion X removes wherever it sits (genome, [lo, hi)), past
        // the window its junctions are read by at each end. A substitution there changes
        // bases only REF's molecules hold, never X's, and away from the junctions it is no
        // allele at X's site (as RJ-24 has a co-annotated sibling claim REF reads only
        // through REF's window): a context on REF's side, lying with REF and with any state
        // keeping those bases, never with X. A shorter deletion's bases are all its window:
        // a row there stays another allele, as RJ-5 has it.
        let deleted_core: Option<(i64, i64)> =
            (pure(&target) && target.ref_al.len() >= target.alt_al.len() + LONG_EVENT as usize)
                .then(|| {
                    let (lo, hi) = run(&target);
                    let k = (target.ref_al.len() - target.alt_al.len()) as i64;
                    (hi - k + JUNCTION_WINDOW, lo + k - JUNCTION_WINDOW)
                })
                .filter(|(a, b)| a < b);
        let in_deleted = |r: &Row| {
            let (on_ref, _) = telling(r);
            let (s, e) = (r.pos + on_ref.start as i64, r.pos + on_ref.end as i64);
            !r.changes_length() && s < e && deleted_core.is_some_and(|(a, b)| a <= s && e <= b)
        };
        // The site's states and the contexts' rows, each once.
        let mut seen: BTreeSet<Row> = BTreeSet::from([target.clone()]);
        let (mut states, mut contexts): (Vec<Row>, Vec<Row>) = (Vec::new(), Vec::new());
        let mut take = |r: &Row, state: bool| -> bool {
            if holds(r) && seen.insert(r.clone()) {
                if state {
                    states.push(r.clone());
                } else {
                    contexts.push(r.clone());
                }
                return true;
            }
            false
        };
        // A pure indel anchored on a substitution X's last base changes no base X names:
        // anchored on the base after it instead (the same haplotype), it is beside X, so it
        // can lie with X and is weighed with every state (a germline indel written right
        // after the block).
        let beside = |r: &Row| -> Row {
            let xe = target.end();
            if target.changes_length() || r.pos != xe - 1 || !pure(r) {
                return r.clone();
            }
            let at = |g: i64| usize::try_from(g - start).ok().filter(|&i| i < reference.len());
            let moved = if r.alt_al.len() > r.ref_al.len() {
                at(xe).map(|i| Row::new(xe, &reference[i..=i], &[&r.alt_al[1..], &reference[i..=i]].concat()))
            } else {
                let k = (r.ref_al.len() - 1) as i64;
                at(xe).zip(at(xe + k)).map(|(i, j)| Row::new(xe, &reference[i..=j], &reference[j..=j]))
            };
            moved.unwrap_or_else(|| r.clone())
        };
        for r in candidates.listed.iter().map(beside) {
            let state = meets(&target, &r) && !in_deleted(&r);
            take(&r, state);
        }
        if target.ref_al.len() == 1 && target.alt_al.len() == 1 {
            // An SNV's site holds every base (a third base is another allele).
            for &b in b"ACGT" {
                take(&Row::new(target.pos, &target.ref_al, &[b]), true);
            }
        }
        // The contexts the reads' own alleles supply (the found source), each once.
        let mut found_contexts: Vec<Row> = Vec::new();
        // The found indels clear of X's site taken as contexts: each stays only where the
        // reads need it (below).
        let mut beside_found: Vec<Row> = Vec::new();
        // An allele the guard found, spelled over a window: a state where the bases it
        // changes meet X's site, else a change beside X (RJ-8: another indel outside the
        // window is a separate event), a context.
        for r in &candidates.discovered {
            let (s0, e0) = ((r.pos - start) as usize, (r.end() - start) as usize);
            if s0 > e0 || e0 > reference.len() {
                continue;
            }
            let spelled = [&reference[..s0], &r.alt_al[..], &reference[e0..]].concat();
            match row_between(&reference, &spelled, start) {
                Some(minimal) if !meets(&target, &minimal) || in_deleted(&minimal) => {
                    if take(&minimal, false) {
                        found_contexts.push(minimal);
                    }
                }
                _ => {
                    take(r, true);
                }
            }
        }
        // Another length of X's repeat, and an indel the reads write recurringly at X's
        // site: other alleles, which explain their reads as neither REF nor X.
        if other_option().runs {
            for r in run_alleles(&target, start, &reference) {
                take(&r, true);
            }
        }
        // The found states, each with the reads that spell it.
        let mut found_states: Vec<(Row, &[usize])> = Vec::new();
        if other_option().found {
            // X's site: where another indel makes a read another allele (its discrimination
            // window, the repeat it changes, perfect or not: #290).
            let (site_s, site_e) = {
                let (s0, e0) = x_site(&target, start, &reference, place(&target));
                (s0 + start, e0 + start)
            };
            for (r, by) in candidates.found.iter().map(|(r, by)| (beside(r), &by[..])) {
                // An indel the reads write recurringly clear of X's site is a separate event
                // (RJ-8): a context explaining the reads that carry it, wherever their
                // aligner wrote it (one written inside the site, its bases placing it
                // beside, fits a state with it there; RJ-27). One meeting the site but not
                // X stays out, as before (another length of the repeat X changes is another
                // allele wherever written, #290).
                if !meets(&target, &r) {
                    let (ps, pe) = place(&r);
                    if (pe <= site_s || ps >= site_e) && take(&r, false) {
                        beside_found.push(r.clone());
                        found_contexts.push(r);
                    }
                    continue;
                }
                // An allele of the same length over the bases a substitution X changes,
                // which reads X's (or REF's) bases there, is X (or REF) with a change
                // beside it, a phased neighbour: each change beside X is a context, never
                // a state.
                if let Some(beside_x) = split_beside(&target, &r) {
                    for c in beside_x {
                        if !meets(&target, &c) && take(&c, false) {
                            found_contexts.push(c);
                        }
                    }
                    continue;
                }
                if in_deleted(&r) {
                    take(&r, false);
                    continue;
                }
                if take(&r, true) {
                    found_states.push((r, by));
                }
            }
        }
        // A change at a base X keeps (REF's and X's alleles read alike there, inside the
        // bases X changes) that the REF molecules carry too is a separate event on both
        // haplotypes, not another allele (RJ-25; operator, 2026-10-14). It explains a read's
        // base there on REF and on X alike: a context lying with REF (and any state keeping
        // that base), and X with it, a state of X's class. A listed, found or discovered
        // allele that is REF or X but for such changes is that, never another allele. The
        // REF molecules carry it where the guard's background holds it, or where it recurs
        // ([`GUARD_MIN_READS`] fragments and [`GUARD_REF_SHARE`] of those reading the base)
        // on reads showing REF's bases where X changes them.
        let kept = kept_bases(&target);
        let spell = |r: &Row| -> Vec<u8> {
            let (s, e) = ((r.pos - start) as usize, (r.end() - start) as usize);
            [&reference[..s], &r.alt_al[..], &reference[e..]].concat()
        };
        let mut x_kept: BTreeSet<Vec<u8>> = BTreeSet::new();
        if !kept.is_empty() {
            let (on_ref, _) = telling(&target);
            let x0 = target.pos - start;
            let x_hap = spell(&target);
            // (REF offset, ALT offset, base) of each change at a kept base a candidate shows
            let mut changes: BTreeSet<(usize, usize, u8)> = BTreeSet::new();
            for r in candidates.background.iter().chain(&candidates.listed) {
                if !r.changes_length() && r.ref_al.len() == 1 {
                    if let Some(&(i, k)) = kept.iter().find(|&&(i, _)| target.pos + i as i64 == r.pos) {
                        changes.insert((i, k, r.alt_al[0]));
                    }
                }
            }
            for s in &states {
                let h = spell(s);
                let diff = |a: &[u8], b: &[u8]| -> Vec<usize> { (0..a.len()).filter(|&c| a[c] != b[c]).collect() };
                if h.len() == reference.len() {
                    let d = diff(&h, &reference);
                    let at: Vec<Option<&(usize, usize)>> = d.iter().map(|&c| kept.iter().find(|&&(i, _)| x0 + i as i64 == c as i64)).collect();
                    if !d.is_empty() && at.iter().all(Option::is_some) {
                        for (&c, k) in d.iter().zip(at) {
                            let &(i, k) = k.expect("checked");
                            changes.insert((i, k, h[c]));
                        }
                    }
                }
                if h.len() == x_hap.len() {
                    let d = diff(&h, &x_hap);
                    let at: Vec<Option<&(usize, usize)>> = d.iter().map(|&c| kept.iter().find(|&&(_, k)| x0 + k as i64 == c as i64)).collect();
                    if !d.is_empty() && at.iter().all(Option::is_some) {
                        for (&c, k) in d.iter().zip(at) {
                            let &(i, k) = k.expect("checked");
                            changes.insert((i, k, h[c]));
                        }
                    }
                }
            }
            // Fragments reading base `b` at genome `p` among those whose own alignments read
            // `want` there (each position its base): (carrying it, reading p).
            let carried = |p: i64, b: u8, want: &[(i64, u8)]| -> (usize, usize) {
                let (mut reading, mut carrying): (BTreeSet<Vec<u8>>, BTreeSet<Vec<u8>>) = Default::default();
                for (n, r) in reads.iter().enumerate() {
                    if !want.iter().all(|&(g, w)| r.aligned_base(g) == Some(w)) {
                        continue;
                    }
                    let Some(c) = r.aligned_base(p) else { continue };
                    let key = if r.name.is_empty() { n.to_string().into_bytes() } else { r.name.clone() };
                    if c == b {
                        carrying.insert(key.clone());
                    }
                    reading.insert(key);
                }
                (carrying.len(), reading.len())
            };
            let recurs = |(c, r): (usize, usize)| c as u32 >= GUARD_MIN_READS && c as f64 >= GUARD_REF_SHARE * r as f64;
            // REF's bases where X changes them, and (a substitution) X's
            let on_ref_side: Vec<(i64, u8)> = on_ref
                .clone()
                .filter(|i| !kept.iter().any(|&(a, _)| a == *i))
                .filter_map(|i| base(target.pos + i as i64).map(|b| (target.pos + i as i64, b)))
                .collect();
            let on_alt_side: Vec<(i64, u8)> = if target.changes_length() {
                Vec::new()
            } else {
                (0..target.ref_al.len()).filter(|&i| target.ref_al[i] != target.alt_al[i]).map(|i| (target.pos + i as i64, target.alt_al[i])).collect()
            };
            let with = |set: &[(usize, usize, u8)]| -> Row {
                let mut alt = target.alt_al.clone();
                for &(_, k, b) in set {
                    alt[k] = b;
                }
                Row::new(target.pos, &target.ref_al, &alt)
            };
            let ref_with = |set: &[(usize, usize, u8)]| -> Vec<u8> {
                let mut h = reference.clone();
                for &(i, _, b) in set {
                    h[(x0 + i as i64) as usize] = b;
                }
                h
            };
            let state_spells: BTreeSet<Vec<u8>> = states.iter().map(spell).collect();
            let listed_spells: BTreeSet<Vec<u8>> = candidates.listed.iter().filter(|r| holds(r)).map(spell).collect();
            // A change on both haplotypes (the REF molecules and the ALT's, each on the
            // floor), or the guard's background there (a base it masks as REF's or the
            // change's) unless the input lists REF with it (a co-annotated sibling, RJ-5).
            let mut both: BTreeSet<(usize, usize, u8)> = BTreeSet::new();
            let changes: Vec<(usize, usize, u8)> = changes
                .into_iter()
                .filter(|&(i, k, b)| {
                    let p = target.pos + i as i64;
                    if b == target.ref_al[i] {
                        return false;
                    }
                    let alt_side = (!on_alt_side.is_empty() && recurs(carried(p, b, &on_alt_side)))
                        || state_spells.contains(&spell(&with(&[(i, k, b)])));
                    if alt_side && recurs(carried(p, b, &on_ref_side)) {
                        both.insert((i, k, b));
                        return true;
                    }
                    let bg = candidates.background.iter().any(|r| !r.changes_length() && r.ref_al.len() == 1 && r.pos == p && r.alt_al[0] == b);
                    bg && !listed_spells.contains(&ref_with(&[(i, k, b)]))
                })
                .collect();
            let mut sets: Vec<Vec<(usize, usize, u8)>> = changes.iter().map(|&c| vec![c]).collect();
            if changes.len() > 1 {
                sets.push(changes.clone());
            }
            for set in &sets {
                let row = with(set);
                x_kept.insert(spell(&row));
                if holds(&row) && seen.insert(row.clone()) {
                    states.push(row);
                }
            }
            if !changes.is_empty() {
                // REF but for such changes: REF with a separate event, a context's work (a
                // listed one stays a state unless the change is on both haplotypes).
                let ref_kept: BTreeSet<Vec<u8>> = sets.iter().map(|set| ref_with(set)).collect();
                let ref_both: BTreeSet<Vec<u8>> =
                    sets.iter().filter(|set| set.iter().all(|c| both.contains(c))).map(|set| ref_with(set)).collect();
                states.retain(|s| {
                    let h = spell(s);
                    !(ref_both.contains(&h) || (ref_kept.contains(&h) && !listed_spells.contains(&h)))
                });
                found_states.retain(|(r, _)| !x_kept.contains(&spell(r)) && !ref_kept.contains(&spell(r)));
                for &(i, _, b) in &changes {
                    let c = Row::new(target.pos + i as i64, &target.ref_al[i..=i], &[b]);
                    seen.insert(c.clone());
                    if !contexts.contains(&c) {
                        contexts.push(c);
                    }
                }
                log::debug!(
                    "haplotype judge {}: {} change(s) at bases it keeps that the REF molecules carry: separate events",
                    named(&target),
                    changes.len(),
                );
            }
        }
        // A co-annotated sibling of X's length at X's site is another allele only where the
        // fragments spelling it in full (past the masks) are at least one and
        // [`GUARD_REF_SHARE`] of those spelling the allele it is nearest, REF or X (the
        // larger when as near): a component the reads never show apart from the event is
        // part of it, and at depth one read spelling it is an error's share (RJ-25;
        // operator, 2026-10-14). Read on the reads' own alignments over the sibling's and
        // X's bases; a sibling changing the length, and one holding X's change (RJ-28), are
        // left as they are.
        {
            let listed_states: Vec<Row> = candidates
                .listed
                .iter()
                .map(beside)
                .filter(|r| meets(&target, r) && !in_deleted(r) && !r.changes_length() && states.contains(r))
                .collect();
            let x_same = !target.changes_length();
            for row in listed_states {
                let holds = x_same
                    && row.pos <= target.pos
                    && row.end() >= target.end()
                    && (0..target.ref_al.len()).all(|k| {
                        target.ref_al[k] == target.alt_al[k] || row.alt_al[(target.pos + k as i64 - row.pos) as usize] == target.alt_al[k]
                    });
                if holds {
                    continue;
                }
                let (u0, u1) = (row.pos.min(target.pos), row.end().max(target.end()));
                let on = |allele: &Row, g: i64| -> Option<u8> {
                    (allele.pos <= g && g < allele.end()).then(|| allele.alt_al[(g - allele.pos) as usize])
                };
                let ref_at = |g: i64| base(g).unwrap_or(WILD);
                let sib_at = |g: i64| on(&row, g).unwrap_or_else(|| ref_at(g));
                let x_at = |g: i64| if x_same { on(&target, g).unwrap_or_else(|| ref_at(g)) } else { WILD };
                // the positions telling the three apart must be read confidently
                let telling: Vec<i64> = (u0..u1).filter(|&g| sib_at(g) != ref_at(g) || (x_same && x_at(g) != ref_at(g))).collect();
                let spells = |r: &ReadView, at: &dyn Fn(i64) -> u8| -> bool {
                    (u0..u1).all(|g| match r.aligned_base(g) {
                        Some(b) => b == at(g),
                        None => !telling.contains(&g) && (0..r.seq.len()).any(|i| r.aligned[i] && r.refpos[i] == g),
                    })
                };
                let mut by_sib: BTreeSet<Vec<u8>> = BTreeSet::new();
                let mut by_ref: BTreeSet<Vec<u8>> = BTreeSet::new();
                let mut by_x: BTreeSet<Vec<u8>> = BTreeSet::new();
                let mut reading = 0usize;
                for (n, r) in reads.iter().enumerate() {
                    let key = if r.name.is_empty() { n.to_string().into_bytes() } else { r.name.clone() };
                    reading += usize::from((u0..u1).all(|g| (0..r.seq.len()).any(|i| r.aligned[i] && r.refpos[i] == g)));
                    if spells(r, &sib_at) {
                        by_sib.insert(key);
                    } else if spells(r, &ref_at) {
                        by_ref.insert(key);
                    } else if x_same && spells(r, &x_at) {
                        by_x.insert(key);
                    }
                }
                let dist = |at: &dyn Fn(i64) -> u8| (u0..u1).filter(|&g| sib_at(g) != at(g)).count();
                let (dr, dx) = (dist(&ref_at), if x_same { dist(&x_at) } else { usize::MAX });
                let near = match dr.cmp(&dx) {
                    std::cmp::Ordering::Less => by_ref.len(),
                    std::cmp::Ordering::Greater => by_x.len(),
                    std::cmp::Ordering::Equal => by_ref.len().max(by_x.len()),
                };
                // with no read over the sibling's and X's bases nothing shows it absent
                let seen_apart = reading == 0 || (!by_sib.is_empty() && by_sib.len() as f64 >= GUARD_REF_SHARE * near as f64);
                log::debug!(
                    "haplotype judge {}: the co-annotated {} spelled on {} fragment(s), its nearest allele on {} → {}",
                    named(&target),
                    named(&row),
                    by_sib.len(),
                    near,
                    if seen_apart { "another allele, a state" } else { "part of the event, not a state" },
                );
                if !seen_apart {
                    states.retain(|r| *r != row);
                }
            }
        }
        // A background change at X's place cannot share a haplotype with X: crossed with
        // REF alone it would favour REF, so it is left to the reads' bases there. Under the
        // default (`respell`) an item of the guard's background or of the found contexts that
        // respells the site's own candidate sequences ([`respellings`]) is left out where no
        // read needs it; the guard's masks stay.
        let background: Vec<Row> = match bg_rule {
            BackgroundRule::Off => Vec::new(),
            _ => candidates.background.iter().map(beside).filter(|r| !meets(&target, r)).collect(),
        };
        let nb = background.len();
        let items: Vec<Row> = background.iter().chain(&found_contexts).cloned().collect();
        // A whole respelling is left out; one over a read only where no reads need it: an
        // item the reads attest (a real change beside a repeat, a germline SNP whose haplotype
        // reads as X over a short read) explains them.
        let settle = |whole: Dropped, partial: Vec<&Row>| -> Dropped {
            if partial.is_empty() {
                return whole;
            }
            let build = |d: &Dropped| {
                Site::build_core(target.clone(), candidates, (start, reference.clone()), reads, rule, bg_rule, Some(d))
            };
            let Some(kept) = build(&whole) else { return whole };
            let mut out = whole.clone();
            for item in partial {
                let mut d = whole.clone();
                d.rows.push(item.clone());
                let attested = build(&d).is_some_and(|without| needed_by(&kept, &without, reads, item));
                log::debug!(
                    "haplotype judge {}: {} respells a candidate sequence over a read; {}",
                    named(&target),
                    named(item),
                    if attested { "the reads need it, kept" } else { "no reads need it, left out" },
                );
                if !attested {
                    out.rows.push(item.clone());
                }
            }
            out
        };
        let drop: Dropped = match dropped {
            Some(d) => d.clone(),
            None if bg_rule == BackgroundRule::Respell && !(items.is_empty() && candidates.masked.is_empty()) => {
                let l = reads.iter().map(|r| r.len() as i64).max().unwrap_or(0);
                let (s0, e0) = x_site(&target, start, &reference, place(&target));
                let all_states: Vec<Row> = std::iter::once(target.clone()).chain(states.iter().cloned()).collect();
                let listed: Vec<Row> = contexts.iter().filter(|r| !found_contexts.contains(r)).cloned().collect();
                let site = (s0 + start, e0 + start);
                // The guard's masks stay whatever they spell: a masked germline SNP whose
                // haplotype reads as X over a short read never counts ALT (decided), and its
                // reads cannot be told from a duplication's carriers at its one telling base.
                let (mut rows, _) = respellings(&reference, start, site, l, &all_states, &listed, &items, &[]);
                // Items that together, with REF, spell one change meeting X (a carrier's
                // event the aligner wrote as indels on either side of X): an allele at X's
                // site, never contexts.
                for (i, joined) in joined_at_site(&reference, start, site, l, &items, |r| meets(&target, r)).into_iter().enumerate() {
                    if joined {
                        rows[i] = 2;
                    }
                }
                let whole = Dropped { rows: items.iter().zip(&rows).filter(|(_, &k)| k == 2).map(|(r, _)| r.clone()).collect() };
                // An indel the reads write beside X's site is a separate event only where the
                // reads need it for their calls, as a respelling over a read is: one the event's
                // own reads write (a long insertion's copy of its flank, read into the
                // reference) is explained by X already, and as a context it would tie them.
                let settled = settle(
                    whole,
                    items.iter().zip(&rows).filter(|(r, &k)| k == 1 || (k == 0 && beside_found.contains(r))).map(|(r, _)| r).collect(),
                );
                if respell_trace() {
                    let (forth, _) = respell_directions(&reference, start, site, l, &all_states, &listed, &items, &[], (true, false));
                    let (back, _) = respell_directions(&reference, start, site, l, &all_states, &listed, &items, &[], (false, true));
                    let (on_x, _) = respellings(&reference, start, site, l, std::slice::from_ref(&target), &listed, &items, &[]);
                    let joined = joined_at_site(&reference, start, site, l, &items, |r| meets(&target, r));
                    for (i, item) in items.iter().enumerate().filter(|&(i, _)| rows[i] > 0) {
                        let kind = if i >= nb { "found" } else if item.changes_length() { "indel" } else { "substitution" };
                        let dir = match (forth[i] > 0, back[i] > 0) {
                            (true, true) => "both",
                            (true, false) => "ref_to_state",
                            (false, true) => "state_to_ref",
                            (false, false) => "joined_only",
                        };
                        log::trace!(
                            "judge respell {} xtype={} rule=respell item={} kind={kind} level={} joined={} dir={dir} state={} fate={}",
                            named(&target), x_type(&target), named(item), if rows[i] == 2 { "whole" } else { "read" },
                            u8::from(joined[i]), if on_x[i] > 0 { "x" } else { "other" },
                            if settled.rows.contains(item) { "left_out" } else { "kept" },
                        );
                    }
                }
                settled
            }
            // Under `asis` (the narrow rule `respell` generalises, kept for measurement), at a
            // substitution X, a recurring indel of the guard's with which REF reads as X over
            // a read reaching X's site is X's change respelled (a repeat X lengthens, written
            // as more units of it beside X: AGCA>TGCT between a GCT run and a GCA run, REF
            // with two GCT units more reads as X into the GCA run), never a context: left out
            // where no reads need it for their call. A substitution's carriers hold its bases
            // at X; an indel there is a shift of the repeat beside it. The masks, the
            // substitutions and the found alleles stay contexts, and an indel X's site keeps
            // its background as it is.
            None if bg_rule == BackgroundRule::AsIs && !target.changes_length() && background.iter().any(|r| r.changes_length()) => {
                let l = reads.iter().map(|r| r.len() as i64).max().unwrap_or(0);
                let (s0, e0) = x_site(&target, start, &reference, place(&target));
                let site = (s0 + start, e0 + start);
                let indels: Vec<Row> = background.iter().filter(|r| r.changes_length()).cloned().collect();
                let others: Vec<Row> = contexts
                    .iter()
                    .filter(|r| !found_contexts.contains(r))
                    .chain(items.iter().filter(|r| !r.changes_length() || found_contexts.contains(r)))
                    .cloned()
                    .collect();
                let (rows, _) = respell_levels(&reference, start, site, l, std::slice::from_ref(&target), &others, &indels, &[], true);
                let whole = Dropped { rows: indels.iter().zip(&rows).filter(|(_, &k)| k == 2).map(|(r, _)| r.clone()).collect() };
                let settled = settle(whole, indels.iter().zip(&rows).filter(|(_, &k)| k == 1).map(|(r, _)| r).collect());
                if respell_trace() {
                    for (i, item) in indels.iter().enumerate().filter(|&(i, _)| rows[i] > 0) {
                        log::trace!(
                            "judge respell {} xtype={} rule=asis item={} kind=indel level={} joined=0 dir=ref_to_state state=x fate={}",
                            named(&target), x_type(&target), named(item), if rows[i] == 2 { "whole" } else { "read" },
                            if settled.rows.contains(item) { "left_out" } else { "kept" },
                        );
                    }
                }
                settled
            }
            None => Dropped::default(),
        };
        let respelt: Vec<bool> = items.iter().map(|r| drop.rows.contains(r)).collect();
        let left_out = |r: &Row, whose: &str| {
            log::debug!(
                "haplotype judge {}:{}>{}: {whose} {}:{}>{} respells a candidate sequence; left out",
                target.pos + 1, String::from_utf8_lossy(&target.ref_al), String::from_utf8_lossy(&target.alt_al),
                r.pos + 1, String::from_utf8_lossy(&r.ref_al), String::from_utf8_lossy(&r.alt_al),
            );
        };
        for (r, &respells) in found_contexts.iter().zip(&respelt[nb..]) {
            if respells {
                left_out(r, "the found allele");
                contexts.retain(|c| c != r);
            }
        }
        for (r, &respells) in background.into_iter().zip(&respelt[..nb]) {
            if respells {
                left_out(&r, "the guard's");
            } else if holds(&r) && seen.insert(r.clone()) {
                contexts.push(r);
            }
        }
        let kept_masks: Vec<i64> = match bg_rule {
            BackgroundRule::Off => Vec::new(),
            _ => candidates.masked.clone(),
        };
        // A state of a substitution X's site that keeps the length may also change bases
        // beside X's (an allele spelled over a window, a listed MNV reaching past X): each
        // such change is a separate event, so it is also a context of every state, and a
        // read carrying X (or REF) with it is X's (or REF's) beside a phased neighbour,
        // never left tied between X and that state.
        if !target.changes_length() {
            let beside_x: Vec<Row> = states
                .iter()
                .filter(|r| !r.changes_length())
                .flat_map(|r| {
                    (0..r.ref_al.len())
                        .filter(|&i| r.ref_al[i] != r.alt_al[i])
                        .map(|i| Row::new(r.pos + i as i64, &r.ref_al[i..=i], &r.alt_al[i..=i]))
                        .filter(|c| c.pos < target.pos || c.pos >= target.end())
                        .collect::<Vec<_>>()
                })
                .collect();
            for c in beside_x {
                if holds(&c) && seen.insert(c.clone()) {
                    contexts.push(c);
                }
            }
        }
        states.sort();
        contexts.sort();
        // A row spelling REF, X, or an earlier row of its kind over the window is that
        // one. (A state spelling a context row's haplotype stays a row: its haplotypes
        // that spell X's or REF's are dropped below, never the listed context.)
        let spelled = |r: &Row| -> Vec<u8> {
            let (s, e) = ((r.pos - start) as usize, (r.end() - start) as usize);
            [&reference[..s], &r.alt_al[..], &reference[e..]].concat()
        };
        let given: BTreeSet<Vec<u8>> = BTreeSet::from([reference.clone(), spelled(&target)]);
        let distinct = |rows: Vec<Row>| -> Vec<Row> {
            let mut seen = given.clone();
            rows.into_iter().filter(|r| seen.insert(spelled(r))).collect()
        };
        let (states, mut contexts) = (distinct(states), distinct(contexts));
        // A pure length change whose run X lies inside sits on the reference at either end
        // of it alike, but beside X they are two molecules: the other end's placement joins
        // the contexts (at one place with the row, so never with it), where it can lie with X.
        let (xs, xe) = (target.pos, target.end());
        let mut others_side: Vec<Row> = Vec::new();
        for r in contexts.iter().filter(|r| pure(r)) {
            let (lo, hi) = run(r);
            if !(xs < hi && lo < xe) {
                continue;
            }
            let k = r.ref_al.len().abs_diff(r.alt_al.len()) as i64;
            let hap = spelled(r);
            let at = |g: i64| usize::try_from(g - start).ok();
            let ends = if r.alt_al.len() > r.ref_al.len() {
                // inserted after the anchor before boundary b: the haplotype's bases there
                [lo, hi].map(|b| {
                    let (a, b0) = (at(b - 1)?, at(b)?);
                    let inserted = hap.get(b0..b0 + k as usize)?;
                    Some(Row::new(b - 1, &reference[a..a + 1], &[&reference[a..a + 1], inserted].concat()))
                })
            } else {
                // the deleted bases [d, d + k) with the anchor before them
                [lo, hi - k].map(|d| {
                    let (a, e) = (at(d - 1)?, at(d + k)?);
                    Some(Row::new(d - 1, &reference[a..e], &reference[a..a + 1]))
                })
            };
            for placed in ends.into_iter().flatten() {
                let beside_x = placed.end() <= xs || placed.pos >= xe;
                if placed != *r && beside_x && holds(&placed) && !others_side.contains(&placed) && !contexts.contains(&placed) {
                    others_side.push(placed);
                }
            }
        }
        contexts.extend(others_side);
        contexts.sort();
        let n_states = states.len();
        let rows: Vec<Row> = std::iter::once(target).chain(states).chain(contexts).collect();
        // A state holding X's change at X's site: for a substitution X, an allele of the
        // same length over X's bases that reads X's bases where X changes them; and X with
        // a change both haplotypes carry at a base it keeps (above).
        let holds_x: Vec<bool> = {
            let x = &rows[0];
            let (on_ref, _) = telling(x);
            rows.iter()
                .enumerate()
                .map(|(r, row)| {
                    (1..=n_states).contains(&r)
                        && ((!x.changes_length()
                            && !row.changes_length()
                            && row.pos <= x.pos
                            && row.end() >= x.end()
                            && on_ref.clone().all(|k| row.alt_al[(x.pos + k as i64 - row.pos) as usize] == x.alt_al[k]))
                            || x_kept.contains(&spelled(row)))
                })
                .collect()
        };
        let spans: Vec<(i64, i64)> = rows.iter().map(|r| (r.pos - start, r.end() - start)).collect();
        let places = rows.iter().map(place).map(|(s, e)| (s - start, e - start)).collect();
        let max_shift = rows.iter().map(|r| r.ref_al.len().abs_diff(r.alt_al.len()) as i64).sum();
        // A mask where a row changes a base is that row's base: the row explains the reads
        // carrying it (a context with every state, or a state of the site). A row spelled
        // over a wider window (an allele the guard found) changes only some of its bases.
        // A base a row of the same length keeps (its two alleles alike there) it does not
        // change: a mask there stays, on that row's haplotypes as on the rest.
        let changed: Vec<(i64, i64)> = rows
            .iter()
            .map(|r| {
                let (on_ref, _) = telling(r);
                (r.pos - start + on_ref.start as i64, r.pos - start + on_ref.end as i64)
            })
            .collect();
        let mut masked: Vec<i64> = kept_masks
            .iter()
            .map(|&p| p - start)
            .filter(|&m| {
                m >= 0
                    && m < wn
                    && rows.iter().zip(&changed).all(|(r, &(s, e))| {
                        !(s <= m && m < e)
                            || (!r.changes_length() && r.ref_al[(m - (r.pos - start)) as usize] == r.alt_al[(m - (r.pos - start)) as usize])
                    })
            })
            .collect();
        masked.sort_unstable();
        masked.dedup();
        // A run of masks (three bases or more in a row) is no germline change beside the
        // event but the shifted bases of reads of another length written without their
        // gap: any base there would make those reads fit X and REF alike. Left unmasked,
        // the reads' bases there are weighed.
        let run_at = |i: usize| {
            let (mut lo, mut hi) = (i, i);
            while lo > 0 && masked[lo - 1] + 1 == masked[lo] {
                lo -= 1;
            }
            while hi + 1 < masked.len() && masked[hi] + 1 == masked[hi + 1] {
                hi += 1;
            }
            hi - lo + 1
        };
        let masked: Vec<i64> = (0..masked.len()).filter(|&i| run_at(i) < MASK_RUN).map(|i| masked[i]).collect();
        let site_window = x_site(&rows[0], start, &reference, place(&rows[0]));
        // Per position and base over the window, the reads whose own alignments read it there:
        // a neither base's recurrence, and a mismatch's elsewhere (measurement only).
        let tally = {
            let x = &rows[0];
            let (on_ref, on_alt) = telling(x);
            let changes = x.ref_al.iter().zip(&x.alt_al).filter(|(a, b)| a != b).count();
            let applicable = if x.changes_length() { !on_ref.is_empty() && !on_alt.is_empty() } else { changes >= 2 };
            let mut t: HashMap<(i64, u8), u32> = HashMap::new();
            if applicable {
                let (lo, hi) = (start, start + reference.len() as i64);
                for r in reads {
                    for i in 0..r.seq.len() {
                        if r.aligned[i] && r.seq[i] != WILD && (lo..hi).contains(&r.refpos[i]) {
                            *t.entry((r.refpos[i], r.seq[i])).or_insert(0) += 1;
                        }
                    }
                }
            }
            t
        };
        let (background_bases, background_masks) = if bg_rule == BackgroundRule::Off {
            (HashSet::new(), Vec::new())
        } else {
            let bases: HashSet<(i64, u8)> = candidates
                .background
                .iter()
                .filter(|r| !r.changes_length())
                .flat_map(|r| {
                    let changed = r.ref_al.iter().zip(&r.alt_al).enumerate().filter(|(_, (a, b))| a != b);
                    changed.map(move |(k, (_, &b))| (r.pos + k as i64, b))
                })
                .collect();
            let mut masks = candidates.masked.clone();
            masks.sort_unstable();
            (bases, masks)
        };
        // X's placements, its anchor included, touching a contig edge: its repeat running to
        // the edge, or X itself at the contig's first or last base (a MAF insertion at Start
        // 0, the base-after form at position 1, anchored on the base after it).
        let edge = if length_change_only(&rows[0]) {
            let (lo, hi) = place(&rows[0]);
            (candidates.edges.0.filter(|&first| lo <= first), candidates.edges.1.filter(|&end| hi >= end))
        } else {
            (None, None)
        };
        let mut site = Site {
            start,
            reference,
            rows,
            n_states,
            holds_x,
            spans,
            places,
            max_shift,
            masked,
            site_window,
            graphs: Vec::new(),
            neither: neither_rule(),
            own_gap: own_gap_rule(),
            ref_base: ref_base_rule(),
            tally,
            background_bases,
            background_masks,
            edge,
        };

        // The context rows a read over X can reach: REF span meeting X's ± the longest read
        // (a long deletion starting there runs past it, yet a read over X can carry its
        // junction).
        let l = reads.iter().map(|r| r.len() as i64).max().unwrap_or(0);
        let (xs, xe) = site.spans[0];
        let rest: Vec<usize> = (1 + site.n_states..site.rows.len())
            .filter(|&r| site.spans[r].1 > xs - l && site.spans[r].0 < xe + l)
            .collect();
        // Every context at once: one graph per state, no enumeration.
        site.graphs = site.state_graphs(&rest);
        site.drop_spelled_states(&rest);
        site.drop_graph_spellings(&found_states, reads, rule);
        site.link_layouts();
        Some(site)
    }

}

/// A zone's alternatives (the sets of its rows at one place that can lie together) up to
/// which it is laid out as parallel copies; past it the zone's rows are plain bubbles, so
/// two of them at one place may share a path (a haplotype more, never one fewer).
const ZONE_MAX_ALTS: usize = 64;
/// Rows inside an other state's span up to which every set of them is tried when asking
/// whether the state spells X or REF beside them.
const SPELL_MAX_ROWS: usize = 6;

/// The graph scorer's side of a site.
impl Site {
    /// The graph per state of X's site: REF's bases, X, then each other state; each with
    /// every row of `rest` (the context rows in the reads' reach) that can lie with the
    /// state as a bubble.
    fn state_graphs(&self, rest: &[usize]) -> Vec<StateGraph> {
        let mut states: Vec<(Class, Option<usize>)> = vec![(Class::Ref, None), (Class::Alt, Some(0))];
        states.extend((1..=self.n_states).map(|r| (if self.holds_x[r] { Class::Alt } else { Class::Other }, Some(r))));
        states
            .into_iter()
            .map(|(class, row)| {
                let optional: Vec<usize> = rest.iter().copied().filter(|&c| row.is_none_or(|s| !self.conflict(s, c))).collect();
                StateGraph { class, row, graph: self.graph_of(row, &optional), share: None }
            })
            .collect()
    }

    /// Row `r`'s bases as its haplotypes spell them: its ALT, a base it keeps (a row of
    /// one length, its alleles alike there) taking any base where the guard masks it, as
    /// REF's haplotypes take it there.
    fn row_bases(&self, r: usize) -> Vec<u8> {
        let row = &self.rows[r];
        let mut b = row.alt_al.clone();
        if !row.changes_length() {
            let s = self.spans[r].0;
            for (i, c) in b.iter_mut().enumerate() {
                if row.ref_al[i] == *c && self.masked.binary_search(&(s + i as i64)).is_ok() {
                    *c = WILD;
                }
            }
        }
        b
    }

    fn spans_meet(&self, a: usize, b: usize) -> bool {
        let ((s, e), (t, f)) = (self.spans[a], self.spans[b]);
        s < f && t < e
    }

    /// The window's graph with `forced` applied and `optional` as bubbles. A path cannot
    /// take two rows whose REF spans meet; two rows at one place whose spans do not meet
    /// (two length changes in one repeat) are kept apart by a zone: the window over them
    /// laid out once per set of them that can lie together, each copy with the zone's
    /// other rows as bubbles.
    pub(crate) fn graph_of(&self, forced: Option<usize>, optional: &[usize]) -> Graph {
        let mut all: Vec<usize> = forced.into_iter().chain(optional.iter().copied()).collect();
        all.sort_by_key(|&r| (self.spans[r], r));
        all.dedup();
        // rows at one place whose spans do not meet, joined
        let mut parent: Vec<usize> = (0..all.len()).collect();
        fn find(p: &mut [usize], x: usize) -> usize {
            let mut r = x;
            while p[r] != r {
                r = p[r];
            }
            p[x] = r;
            r
        }
        let mut paired = vec![false; all.len()];
        for i in 0..all.len() {
            for j in i + 1..all.len() {
                let (a, b) = (all[i], all[j]);
                if self.conflict(a, b) && !self.spans_meet(a, b) {
                    let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                    parent[ri] = rj;
                    paired[i] = true;
                    paired[j] = true;
                }
            }
        }
        // zones: [lo, hi) and the rows at one place in it
        let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for i in 0..all.len() {
            if paired[i] {
                let r = find(&mut parent, i);
                groups.entry(r).or_default().push(all[i]);
            }
        }
        let mut zones: Vec<((i64, i64), Vec<usize>)> = groups
            .into_values()
            .map(|rows| {
                let lo = rows.iter().map(|&r| self.spans[r].0).min().unwrap_or(0);
                let hi = rows.iter().map(|&r| self.spans[r].1).max().unwrap_or(0);
                ((lo, hi), rows)
            })
            .collect();
        // merged where they meet, widened over every row reaching across an end
        loop {
            zones.sort_by_key(|z| z.0);
            let mut merged: Vec<((i64, i64), Vec<usize>)> = Vec::new();
            for z in zones {
                match merged.last_mut() {
                    Some(m) if z.0 .0 < m.0 .1 => {
                        m.0 .1 = m.0 .1.max(z.0 .1);
                        m.1.extend(z.1);
                    }
                    _ => merged.push(z),
                }
            }
            let mut grew = false;
            for z in merged.iter_mut() {
                for &r in &all {
                    let (s, e) = self.spans[r];
                    if s < z.0 .1 && e > z.0 .0 && (s < z.0 .0 || e > z.0 .1) {
                        z.0 = (z.0 .0.min(s), z.0 .1.max(e));
                        grew = true;
                    }
                }
            }
            zones = merged;
            if !grew {
                break;
            }
        }
        // each zone's alternatives; a zone with too many is left to plain bubbles
        // per zone, its span and the sets of its rows at one place laid out side by side
        type Zone = ((i64, i64), Vec<Vec<usize>>);
        let laid: Vec<Zone> = zones
            .into_iter()
            .filter_map(|(span, at_one)| {
                let alts = self.lying_together(&at_one, forced, ZONE_MAX_ALTS);
                if alts.is_none() {
                    log::debug!("haplotype judge: a zone of {} rows at one place has over {} sets; its rows are bubbles", at_one.len(), ZONE_MAX_ALTS);
                }
                alts.map(|a| (span, a))
            })
            .collect();
        let wn = self.reference.len() as i64;
        let inside = |r: usize, (lo, hi): (i64, i64)| self.spans[r].0 >= lo && self.spans[r].1 <= hi;
        let mut gb = GraphBuilder::default();
        let mut prev: Vec<usize> = Vec::new();
        let mut pos = 0i64;
        fn join(gb: &mut GraphBuilder, prev: &mut Vec<usize>, entry: Vec<usize>, exit: Vec<usize>) {
            for &a in prev.iter() {
                for &b in &entry {
                    gb.edge(a, b);
                }
            }
            *prev = exit;
        }
        for ((lo, hi), alts) in &laid {
            if pos < *lo {
                let f: Vec<usize> = forced.into_iter().filter(|&r| inside(r, (pos, *lo))).collect();
                let o: Vec<usize> = all.iter().copied().filter(|&r| Some(r) != forced && inside(r, (pos, *lo))).collect();
                let (entry, exit) = self.region(&mut gb, (pos, *lo), &f, &o);
                join(&mut gb, &mut prev, entry, exit);
            }
            let zone_rows: Vec<usize> = all.iter().copied().filter(|&r| inside(r, (*lo, *hi))).collect();
            let (mut entries, mut exits) = (Vec::new(), Vec::new());
            for set in alts {
                let mut f: Vec<usize> = set.clone();
                if let Some(x) = forced.filter(|&x| inside(x, (*lo, *hi)) && !set.contains(&x)) {
                    f.push(x);
                }
                // the zone's other rows: those at one place with none of the set's
                let o: Vec<usize> = zone_rows
                    .iter()
                    .copied()
                    .filter(|&r| !f.contains(&r) && !alts.iter().flatten().any(|&c| c == r))
                    .filter(|&r| f.iter().all(|&c| !self.conflict(c, r)))
                    .collect();
                let (entry, exit) = self.region(&mut gb, (*lo, *hi), &f, &o);
                entries.extend(entry);
                exits.extend(exit);
            }
            join(&mut gb, &mut prev, entries, exits);
            pos = *hi;
        }
        if pos < wn {
            let f: Vec<usize> = forced.into_iter().filter(|&r| inside(r, (pos, wn))).collect();
            let o: Vec<usize> = all.iter().copied().filter(|&r| Some(r) != forced && inside(r, (pos, wn))).collect();
            let (entry, exit) = self.region(&mut gb, (pos, wn), &f, &o);
            join(&mut gb, &mut prev, entry, exit);
        }
        gb.finish()
    }

    /// The sets of `rows` (rows at one place) that can lie together, `forced` in each
    /// where it is one of them; None past `cap`.
    fn lying_together(&self, rows: &[usize], forced: Option<usize>, cap: usize) -> Option<Vec<Vec<usize>>> {
        let mut rows: Vec<usize> = rows.to_vec();
        rows.sort_unstable();
        rows.dedup();
        let must = forced.filter(|f| rows.contains(f));
        let mut out: Vec<Vec<usize>> = Vec::new();
        fn rec(site: &Site, rows: &[usize], i: usize, chosen: &mut Vec<usize>, must: Option<usize>, out: &mut Vec<Vec<usize>>, cap: usize) -> bool {
            if out.len() > cap {
                return false;
            }
            if i == rows.len() {
                if must.is_none_or(|m| chosen.contains(&m)) {
                    out.push(chosen.clone());
                }
                return out.len() <= cap;
            }
            let r = rows[i];
            if chosen.iter().all(|&c| !site.conflict(c, r)) {
                chosen.push(r);
                let ok = rec(site, rows, i + 1, chosen, must, out, cap);
                chosen.pop();
                if !ok {
                    return false;
                }
            }
            if Some(r) == must {
                return true;
            }
            rec(site, rows, i + 1, chosen, must, out, cap)
        }
        rec(self, &rows, 0, &mut Vec::new(), must, &mut out, cap).then_some(out)
    }

    /// The window `[a, b)` as reference segments with `forced` applied (the reference cut
    /// under them) and `optional` as bubbles (those meeting a forced row left out). Its
    /// entry segments (starting at `a`) and exit segments (ending at `b`).
    fn region(&self, gb: &mut GraphBuilder, (a, b): (i64, i64), forced: &[usize], optional: &[usize]) -> (Vec<usize>, Vec<usize>) {
        let rows: Vec<usize> = forced
            .iter()
            .copied()
            .chain(optional.iter().copied().filter(|&o| forced.iter().all(|&f| !self.spans_meet(f, o) && f != o)))
            .collect();
        let mut cuts: BTreeSet<i64> = BTreeSet::from([a, b]);
        for &r in &rows {
            cuts.insert(self.spans[r].0);
            cuts.insert(self.spans[r].1);
        }
        let cuts: Vec<i64> = cuts.into_iter().collect();
        let covered = |p: i64| forced.iter().any(|&f| self.spans[f].0 <= p && p < self.spans[f].1);
        let (mut starting, mut ending): (BTreeMap<i64, Vec<usize>>, BTreeMap<i64, Vec<usize>>) = Default::default();
        for w in cuts.windows(2) {
            let (p, q) = (w[0], w[1]);
            if covered(p) {
                continue;
            }
            let bases: Vec<u8> = (p..q)
                .map(|c| if self.masked.binary_search(&c).is_ok() { WILD } else { self.reference[c as usize] })
                .collect();
            let id = gb.segment(None, (p, q), bases);
            starting.entry(p).or_default().push(id);
            ending.entry(q).or_default().push(id);
        }
        for &r in &rows {
            let (s, e) = self.spans[r];
            if self.rows[r].alt_al.is_empty() {
                continue;
            }
            let id = gb.segment(Some(r), (s, e), self.row_bases(r));
            starting.entry(s).or_default().push(id);
            ending.entry(e).or_default().push(id);
        }
        for (p, ins) in &ending {
            if *p == b {
                continue;
            }
            if let Some(outs) = starting.get(p) {
                for &x in ins {
                    for &y in outs {
                        gb.edge(x, y);
                    }
                }
            }
        }
        (starting.remove(&a).unwrap_or_default(), ending.remove(&b).unwrap_or_default())
    }

    /// The haplotype the graph of `state` spells along a path taking `rows`: the state's
    /// row with them.
    pub(crate) fn hap_on_path(&self, state: usize, rows: &[usize]) -> Hap {
        let mut all: Vec<usize> = self.graphs[state].row.into_iter().chain(rows.iter().copied()).collect();
        all.sort_unstable();
        all.dedup();
        Hap::build(self, &all)
    }

    /// The read as the graph scorer reads it: its aligned positions in window coordinates,
    /// made non-decreasing.
    pub(crate) fn graph_read<'r>(&self, read: &'r ReadView) -> GraphRead<'r> {
        let mut rp: Vec<i64> = Vec::with_capacity(read.refpos.len());
        let mut last = i64::MIN;
        for &p in &read.refpos {
            last = last.max(p - self.start);
            rp.push(last);
        }
        let first = read.refpos.iter().copied().min().unwrap_or(read.ref_start).min(read.ref_start);
        let end = read.refpos.iter().copied().max().unwrap_or(read.ref_end).max(read.ref_end);
        GraphRead { seq: &read.seq, rp, span: (first - self.start, end - self.start) }
    }

    /// Per state graph, the first earlier one laid out alike that shares from none itself,
    /// and where their bases first differ ([`StateGraph::share`]).
    fn link_layouts(&mut self) {
        for k in 0..self.graphs.len() {
            let share = (0..k).filter(|&s| self.graphs[s].share.is_none()).find_map(|s| self.graphs[s].graph.layout_until(&self.graphs[k].graph).map(|p| (s, p)));
            self.graphs[k].share = share;
        }
    }

    /// An other state whose allele spells X's haplotype (or REF's) with rows of the reads'
    /// reach inside its place (an allele the reads spell over X with a listed neighbour
    /// beside it) is that haplotype: X's or REF's, never another allele's. Its graph is
    /// dropped.
    fn drop_spelled_states(&mut self, rest: &[usize]) {
        let mut dropped: Vec<usize> = Vec::new();
        for (k, st) in self.graphs.iter().enumerate() {
            let Some(o) = st.row.filter(|_| st.class == Class::Other) else { continue };
            let (ps, pe) = (self.spans[o].0.min(self.places[o].0), self.spans[o].1.max(self.places[o].1));
            let inside: Vec<usize> = rest.iter().copied().filter(|&r| self.spans[r].0 >= ps && self.spans[r].1 <= pe).take(SPELL_MAX_ROWS).collect();
            let own = Hap::build(self, &[o]).seq;
            let sets = self.lying_together(&inside, None, 1 << SPELL_MAX_ROWS).unwrap_or_else(|| vec![Vec::new()]);
            let spells = sets.iter().any(|c| {
                let with_x = c.iter().all(|&r| !self.conflict(0, r)) && {
                    let mut h: Vec<usize> = std::iter::once(0).chain(c.iter().copied()).collect();
                    h.sort_unstable();
                    Hap::build(self, &h).seq == own
                };
                with_x || Hap::build(self, c).seq == own
            });
            if spells {
                dropped.push(k);
            }
        }
        for k in dropped.into_iter().rev() {
            self.graphs.remove(k);
        }
    }

    /// A found allele stays a state only where the reads tell it from X and from REF by
    /// [`DELTA`], the state's graph beating REF's and X's, on the found floor of reads
    /// (`told`: [`GUARD_MIN_READS`] and [`GUARD_REF_SHARE`] of the site's; `toldK`: K reads),
    /// those that write it first, then the rest.
    fn drop_graph_spellings(&mut self, found: &[(Row, &[usize])], reads: &[ReadView], rule: FoundRule) {
        let floor = match rule {
            FoundRule::All => return,
            FoundRule::Told(k) if k > 0 => k as f64,
            FoundRule::Told(_) => (GUARD_REF_SHARE * reads.len() as f64).max(f64::from(GUARD_MIN_READS)),
        };
        let mut ws = Workspace::default();
        let mut dropped: Vec<usize> = Vec::new();
        let given: Vec<usize> = (0..self.graphs.len()).filter(|&k| self.graphs[k].class != Class::Other).collect();
        for (row, by) in found {
            let Some(o) = self.rows.iter().position(|r| r == row).filter(|&o| self.is_state(o)) else {
                continue;
            };
            let Some(k) = self.graphs.iter().position(|g| g.row == Some(o)) else { continue };
            let mut told = 0usize;
            let rest = (0..reads.len()).filter(|i| by.binary_search(i).is_err());
            for read in by.iter().copied().chain(rest).filter_map(|i| reads.get(i)) {
                let gr = self.graph_read(read);
                let own = graph::score(&ShortReads, &gr, &self.graphs[k].graph, &mut ws);
                let tells = given.iter().all(|&g| graph::score(&ShortReads, &gr, &self.graphs[g].graph, &mut ws) <= own - DELTA);
                if tells {
                    told += 1;
                    if told as f64 >= floor {
                        break;
                    }
                }
            }
            let keep = told as f64 >= floor;
            log::debug!(
                "haplotype judge: the found allele {}:{}>{} at {}:{}>{} (written by {} of {} reads, told on {told}) → {}",
                row.pos + 1,
                String::from_utf8_lossy(&row.ref_al),
                String::from_utf8_lossy(&row.alt_al),
                self.rows[0].pos + 1,
                String::from_utf8_lossy(&self.rows[0].ref_al),
                String::from_utf8_lossy(&self.rows[0].alt_al),
                by.len(),
                reads.len(),
                if keep { "another allele, a state" } else { "X's or REF's spelling, not a state" },
            );
            if !keep {
                dropped.push(k);
            }
        }
        dropped.sort_unstable();
        dropped.dedup();
        let gone: Vec<usize> = dropped.iter().filter_map(|&k| self.graphs[k].row).collect();
        for k in dropped.into_iter().rev() {
            self.graphs.remove(k);
        }
        if !gone.is_empty() {
            self.max_shift = (0..self.rows.len())
                .filter(|r| !gone.contains(r))
                .map(|r| self.rows[r].ref_al.len().abs_diff(self.rows[r].alt_al.len()) as i64)
                .sum();
        }
    }

    /// Whether the genome span `[s, e)` of a read reaches a contig edge X's placements touch
    /// ([`Site::edge`]; RJ-3).
    pub(crate) fn edge_reached(&self, s: i64, e: i64) -> bool {
        self.edge.0.is_some_and(|first| s <= first) || self.edge.1.is_some_and(|end| e >= end)
    }
}

/// Aligns one read to a haplotype of the site, trimmed to the read's reach and banded
/// around its own alignment mapped into the haplotype (the alignments a call is examined
/// on).
pub(crate) struct Scorer<'r> {
    read: &'r ReadView,
    /// The read's reference positions, first and last (soft clips extrapolated).
    span: (i64, i64),
}

impl<'r> Scorer<'r> {
    pub(crate) fn new(read: &'r ReadView) -> Self {
        let first = read.refpos.iter().copied().min().unwrap_or(read.ref_start);
        let last = read.refpos.iter().copied().max().unwrap_or(read.ref_end);
        Scorer { read, span: (first.min(read.ref_start), last.max(read.ref_end)) }
    }

    /// The rows `hap` applies whose REF span meets the read's reference span widened by
    /// `by` (window coordinates).
    fn rows_within(&self, site: &Site, hap: &Hap, by: i64) -> Vec<usize> {
        let (lo, hi) = (self.span.0 - site.start - by, self.span.1 - site.start + by);
        hap.rows.iter().copied().filter(|&r| site.spans[r].0 < hi && site.spans[r].1 > lo).collect()
    }

    /// The haplotype trimmed to the read's reach (and where the trim starts), and the
    /// read's band in it: per row of the programme the column its own alignment puts it
    /// at, the band's half-width the slack plus the length changes the rows within the
    /// band (`near`) carry.
    fn prepare<'h>(&self, site: &Site, hap: &'h Hap, near: &[usize]) -> (&'h [u8], usize, Band) {
        self.prepare_with(site, hap, near, BAND_SLACK)
    }

    /// [`Scorer::prepare`] with `slack` columns each side of the diagonal.
    fn prepare_with<'h>(&self, site: &Site, hap: &'h Hap, near: &[usize], slack: i64) -> (&'h [u8], usize, Band) {
        let read = self.read;
        let l = read.len() as i64;
        let (a, b) = hap.reach(read.ref_start - site.start, read.ref_end - site.start, l);
        let wn = hap.window_len();
        let at = |p: i64| hap.to_hap((p - site.start).clamp(0, wn), false) - a as i64;
        let mut centers = Vec::with_capacity(read.len() + 1);
        centers.push(at(read.refpos[0]));
        centers.extend(read.refpos.iter().map(|&p| at(p) + 1));
        let shifts: i64 = near.iter().map(|&r| site.rows[r].ref_al.len().abs_diff(site.rows[r].alt_al.len()) as i64).sum();
        (&hap.seq[a..b], a, Band::around(&centers, slack + shifts, b - a))
    }

    /// The rows within the band's reach: its slack and every length change the site's
    /// rows could carry, twice (the band's width and the read's mapped span).
    fn near(&self, site: &Site, hap: &Hap) -> Vec<usize> {
        self.rows_within(site, hap, BAND_SLACK + 2 * site.max_shift + 2)
    }

    /// The read's alignment in `hap`, and where its trimmed haplotype starts (add it to an
    /// aligned column for the column in `hap`).
    pub(crate) fn align(&self, site: &Site, hap: &Hap) -> (Alignment, usize) {
        let near = self.near(site, hap);
        let (seq, a, band) = self.prepare(site, hap, &near);
        (score::align(&self.read.seq, seq, Some(&band)), a)
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// Deterministic pseudo-random bases.
    pub(crate) fn random_seq(n: usize, mut state: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(state >> 33) as usize % 4]
            })
            .collect()
    }

    /// A 600-base contig with `motif` written at 0-based `at`.
    pub(crate) fn contig(seed: u64, at: usize, motif: &[u8]) -> Vec<u8> {
        let mut s = random_seq(600, seed);
        s[at..at + motif.len()].copy_from_slice(motif);
        s
    }

    /// A row at 1-based `pos`, as a VCF writes it.
    pub(crate) fn row(pos: i64, r: &str, a: &str) -> Row {
        Row::new(pos - 1, r.as_bytes(), a.as_bytes())
    }

    /// The site of `target` among `listed` rows on `contig`, with reads of `l` bases.
    pub(crate) fn site(contig: &[u8], target: &Row, listed: &[Row], l: usize) -> Site {
        let others: Vec<Row> = listed.iter().filter(|r| *r != target).cloned().collect();
        let read = ReadView::new(&vec![b'A'; l], &vec![30; l], 0, &[Cigar::Match(l as u32)], (0, l), 20).unwrap();
        let candidates = Candidates { listed: others, ..Candidates::default() };
        Site::build(target.clone(), &candidates, (0, contig.to_vec()), &[read]).unwrap()
    }

    #[test]
    fn the_kept_other_allele_tokens_parse() {
        let o = parse_other("found,fit,told,edges");
        assert!(o.found && o.fit && o.edges && !o.runs);
        assert_eq!((o.rule, o.found_min), (FoundRule::Told(0), FOUND_MIN_READS));
        let o = parse_other("found3,fit,told2,edges");
        assert_eq!((o.found, o.rule, o.found_min), (true, FoundRule::Told(2), 3));
        let o = parse_other("none");
        assert_eq!((o.found, o.fit, o.rule), (false, false, FoundRule::All));
    }

    #[test]
    fn every_context_source_adds_its_rows() {
        let c = random_seq(600, 31);
        let other = |p: usize| if c[p] == b'A' { b"C".to_vec() } else { b"A".to_vec() };
        let x = Row::new(300, &c[300..301], &other(300));
        // A listed deletion over X (a state) and a listed SNV beside it (a context).
        let over_x = Row::new(298, &c[298..303], &c[298..299]);
        let beside = Row::new(330, &c[330..331], &other(330));
        // The guard's: a recurring change, a mask, and another allele spelled over a window
        // that changes one base beside X only.
        let background = Row::new(270, &c[270..271], &other(270));
        let mut spelled = c[340..350].to_vec();
        spelled[5] = other(345)[0];
        let discovered = Row::new(340, &c[340..350], &spelled);
        let candidates = Candidates {
            listed: vec![over_x, beside],
            discovered: vec![discovered],
            background: vec![background],
            masked: vec![250],
            found: Vec::new(),
            edges: (None, None),
        };
        let read = ReadView::new(&[b'A'; 150], &[30; 150], 230, &[Cigar::Match(150)], (0, 150), 20).unwrap();
        // Before the respelling check (the default) leaves out the items no read needs: the
        // narrow background rule keeps every one here. The site's states (the two other
        // bases and the deletion); the contexts (the listed SNV beside X, the guard's
        // recurring change and its allele changing one base beside X); the guard's mask.
        let (window, reads) = ((0, c.clone()), std::slice::from_ref(&read));
        let s = Site::build_backed(x.clone(), &candidates, window, reads, FoundRule::All, BackgroundRule::AsIs)
            .unwrap();
        assert_eq!((s.n_states, s.rows.len() - 1 - s.n_states, s.masked.len()), (3, 3, 1));
    }

    /// The chr5:44868216 A>AT shape (#246 arm, HG002): an insertion of T in an A run
    /// beside a T, with homozygous SNVs at the bases before and two after it and a listed
    /// SNV A>T on the next base. 0-based: A x12 at 290..302, T 302, A 303, A 304, T 305,
    /// A x13 at 306..319, T 319.
    pub(crate) fn respelt_insertion() -> (Vec<u8>, Row, Vec<Row>, Row) {
        let c = contig(43, 290, b"AAAAAAAAAAAATAATAAAAAAAAAAAAAT");
        let x = Row::new(303, b"A", b"AT");
        let listed = vec![Row::new(302, b"T", b"A"), Row::new(304, b"A", b"T"), Row::new(305, b"T", b"A")];
        // The guard's recurring 1 bp T deletion at 305: X with it spells the SNV haplotype.
        let deletion = Row::new(304, b"AT", b"A");
        (c, x, listed, deletion)
    }

    /// The FLT3-ITD shape (gate set): a 36 bp insertion after 0-based 300 copying the next
    /// 36 bases but one (offset 15), whose base the guard masks: REF with that base masked
    /// spells the insertion's first 36 bases.
    pub(crate) fn respelt_duplication() -> (Vec<u8>, Row, i64) {
        let c = random_seq(600, 61);
        let mut copy = c[301..337].to_vec();
        copy[15] = if copy[15] == b'A' { b'C' } else { b'A' };
        let x = Row::new(300, &c[300..301], &[&c[300..301], &copy[..]].concat());
        (c, x, 316)
    }

    /// The chr16:23024630 TA>AC shape (HG002, the clip audit): its carriers hold X with a
    /// 6 bp deletion beside it, which the reads write as a 7 bp deletion and a C inserted
    /// after X. 0-based, from 290: ATTGTAGTTT, X's TA at 299, ACTTGCATTT.
    pub(crate) fn split_carrier() -> (Vec<u8>, Row, Row, [Row; 2]) {
        let c = contig(47, 290, b"ATTGTAGTTTAACTTGCATTT");
        let x = Row::new(299, b"TA", b"AC");
        let carrier = Row::new(293, b"GTAGTTTA", b"AC");
        (c, x, carrier, [Row::new(290, b"ATTGTAGT", b"A"), Row::new(300, b"A", b"AC")])
    }

    #[test]
    fn items_spelling_one_change_at_xs_site_together_are_found() {
        let (c, x, _, [deletion, raw]) = split_carrier();
        // the insertion as the site takes it, anchored past X's last base
        let insertion = Row::new(301, b"A", &[b"C", &c[301..302]].concat());
        assert_eq!((raw.pos, &raw.alt_al[..]), (300, &b"AC"[..]));
        assert_eq!(insertion.alt_al, b"CA");
        let far = Row::new(420, &c[420..421], if c[420] == b'A' { b"C" } else { b"A" });
        let meets_x = |r: &Row| r.pos < x.end() && x.pos < r.end();
        let got = joined_at_site(&c, 0, (x.pos, x.end()), 100, &[deletion.clone(), insertion.clone(), far.clone()], meets_x);
        assert_eq!(got, vec![true, true, false], "the deletion and the insertion spell X's carrier together");
        // alone, neither does
        assert_eq!(joined_at_site(&c, 0, (x.pos, x.end()), 100, &[deletion, far], meets_x), vec![false, false]);
    }

    #[test]
    fn the_default_background_rule_is_respell() {
        // Decided 2026-10-27: `respell` is the default; `asis` and `off` stay for measurement,
        // and an unknown value falls back to the default.
        assert_eq!(parse_background(""), Some(BackgroundRule::Respell));
        assert_eq!(parse_background("respell"), Some(BackgroundRule::Respell));
        assert_eq!(parse_background("asis"), Some(BackgroundRule::AsIs));
        assert_eq!(parse_background("off"), Some(BackgroundRule::Off));
        assert_eq!(parse_background("other"), None);
    }

    #[test]
    fn a_background_item_respelling_a_candidate_sequence_is_found() {
        // A recurring deletion with which X spells REF's SNV haplotype; a mask on the one
        // base telling a duplication from the bases it copies.
        let (c, x, listed, deletion) = respelt_insertion();
        let far = Row::new(420, &c[420..421], if c[420] == b'A' { b"C" } else { b"A" });
        let near_snv = Row::new(298, &c[298..299], b"G");
        let site = |x: &Row| (x.pos - 12, x.end() + 20);
        let (rows, masks) = respellings(&c, 0, site(&x), 100, std::slice::from_ref(&x), &listed, &[deletion, far, near_snv], &[]);
        assert_eq!(rows, vec![2, 0, 0], "only the deletion respells, turning X into REF's SNV haplotype whole");
        assert!(masks.is_empty());
        let (c, x, telling) = respelt_duplication();
        let site = (x.pos - 1, x.pos + 17);
        let (_, masks) = respellings(&c, 0, site, 100, std::slice::from_ref(&x), &[], &[], &[telling, 250, 330]);
        assert_eq!(masks, vec![1, 0, 0], "only the mask on the telling base, over a read");
    }

    /// An indel the reads write beside X's site is a context only where the reads need it
    /// for their calls (as a respelling over a read is): a one-base deletion 40 bases past
    /// an SNV, on six REF reads, costs REF and X alike and changes no read's call, so it
    /// explains none of their differences and is left out. One the event's own reads write
    /// (a long insertion's copy of its flank, read into the reference: HG002 chr16:34590143,
    /// a 1,795-bp insertion whose reads write a 1-bp insertion four bases before it) would,
    /// as a context, tie the carriers it does not explain.
    #[test]
    fn a_found_indel_beside_the_site_no_read_needs_is_no_context() {
        let c = random_seq(600, 31);
        let x = Row::new(299, &[c[299]], &[if c[299] == b'A' { b'C' } else { b'A' }]);
        let del = Row::new(338, &c[338..340], &c[338..339]);
        let alt: Vec<u8> = [&c[..299], &x.alt_al[..], &c[300..]].concat();
        let read = |h: &[u8], s: usize, with_del: bool| {
            if with_del {
                let seq: Vec<u8> = [&h[s..339], &h[340..s + 101]].concat();
                let a = 339 - s;
                ReadView::new(&seq, &[30; 100], s as i64, &[Cigar::Match(a as u32), Cigar::Del(1), Cigar::Match((100 - a) as u32)], (0, 100), 20).unwrap()
            } else {
                ReadView::new(&h[s..s + 100], &[30; 100], s as i64, &[Cigar::Match(100)], (0, 100), 20).unwrap()
            }
        };
        let mut reads: Vec<ReadView> = (0..10).map(|i| read(&c, 250 + i, i < 6)).collect();
        reads.extend((0..10).map(|i| read(&alt, 250 + i, false)));
        let candidates = Candidates { found: vec![(del.clone(), (0..6).collect())], ..Candidates::default() };
        let s = Site::build(x.clone(), &candidates, (0, c.clone()), &reads).unwrap();
        assert!(!s.rows.contains(&del), "no read needs the deletion: no context");
    }
}

/// Tests of the call examined on the traced path (`graph::kept_steps`, `Site::alignment_on_path`).
#[cfg(test)]
mod path_tests {
    use super::super::decide::judge;
    use super::tests::random_seq;
    use super::*;

    /// `c` with `rows` applied, and per base its reference offset (None inserted).
    fn applied(c: &[u8], rows: &[&Row]) -> (Vec<u8>, Vec<Option<i64>>) {
        let mut rows: Vec<&Row> = rows.to_vec();
        rows.sort_by_key(|r| r.pos);
        let (mut h, mut o, mut at) = (Vec::new(), Vec::new(), 0i64);
        for r in rows {
            h.extend_from_slice(&c[at as usize..r.pos as usize]);
            o.extend((at..r.pos).map(Some));
            let shared = r.ref_al.len().min(r.alt_al.len());
            h.extend_from_slice(&r.alt_al);
            o.extend((0..r.alt_al.len()).map(|k| (k < shared).then_some(r.pos + k as i64)));
            at = r.end();
        }
        h.extend_from_slice(&c[at as usize..]);
        o.extend((at..c.len() as i64).map(Some));
        (h, o)
    }

    /// A read of `len` bases of a haplotype from offset `s`, its CIGAR from the bases'
    /// origins (`gapless`: written with no gap at its first base's origin, as an aligner
    /// without the indel would), its quality 2 at `masked` offsets.
    fn read_on(h: &[u8], o: &[Option<i64>], s: usize, len: usize, gapless: bool, masked: &[usize]) -> Option<ReadView> {
        let s = (s..h.len()).find(|&i| o[i].is_some())?;
        let mut e = (s + len).min(h.len());
        while e > s && o[e - 1].is_none() {
            e -= 1;
        }
        let seq = &h[s..e];
        let quals: Vec<u8> = (0..seq.len()).map(|i| if masked.contains(&i) { 2 } else { 30 }).collect();
        let start = o[s]?;
        let cigar = if gapless {
            vec![Cigar::Match(seq.len() as u32)]
        } else {
            let mut ops: Vec<Cigar> = Vec::new();
            let mut last: Option<i64> = None;
            let push = |ops: &mut Vec<Cigar>, op: Cigar| match (ops.last_mut(), op) {
                (Some(Cigar::Match(n)), Cigar::Match(1)) => *n += 1,
                (Some(Cigar::Ins(n)), Cigar::Ins(1)) => *n += 1,
                _ => ops.push(op),
            };
            for &x in &o[s..e] {
                match x {
                    None => push(&mut ops, Cigar::Ins(1)),
                    Some(p) => {
                        if let Some(l) = last.filter(|&l| p > l + 1) {
                            ops.push(Cigar::Del((p - l - 1) as u32));
                        }
                        push(&mut ops, Cigar::Match(1));
                        last = Some(p);
                    }
                }
            }
            ops
        };
        ReadView::new(seq, &quals, start, &cigar, (0, seq.len()), 20)
    }

    fn roll(state: &mut u64, n: u64) -> u64 {
        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (*state >> 33) % n
    }

    /// A random row at `p`: an SNV, an MNP, an insertion or a deletion.
    fn random_row(c: &[u8], p: usize, st: &mut u64) -> Row {
        let other = |b: u8, k: u64| *b"ACGT".iter().filter(|&&x| x != b).nth(k as usize % 3).unwrap();
        match roll(st, 4) {
            0 => Row::new(p as i64, &c[p..=p], &[other(c[p], roll(st, 3))]),
            1 => {
                let n = 2 + roll(st, 2) as usize;
                let alt: Vec<u8> = c[p..p + n].iter().map(|&b| other(b, roll(st, 3))).collect();
                Row::new(p as i64, &c[p..p + n], &alt)
            }
            2 => {
                let n = 1 + roll(st, 6) as usize;
                let ins: Vec<u8> = (0..n).map(|_| b"ACGT"[roll(st, 4) as usize]).collect();
                Row::new(p as i64, &c[p..=p], &[&c[p..=p], &ins[..]].concat())
            }
            _ => {
                let n = 1 + roll(st, 8) as usize;
                Row::new(p as i64, &c[p..=p + n], &c[p..=p])
            }
        }
    }

    /// An alignment's columns as haplotype columns (its trim's start added), its costs and skips.
    fn laid(a: &(Alignment, usize)) -> (Vec<Placed>, Vec<i32>, Vec<u32>) {
        let cols = a.0.placed.iter().map(|p| if let Placed::At(c) = *p { Placed::At(c + a.1) } else { *p }).collect();
        (cols, a.0.costs.clone(), a.0.skipped_after.clone())
    }

    /// An alignment's score from its own costs and deletions.
    fn summed(a: &Alignment) -> i32 {
        let gaps: i32 = a.skipped_after.iter().filter(|&&k| k > 0).map(|&k| -(score::GAP_OPEN + score::GAP_EXTEND * k as i32)).sum();
        a.costs.iter().sum::<i32>() + gaps
    }

    #[test]
    fn a_reads_call_is_examined_on_the_alignment_its_score_came_from() {
        // The call (X's telling bases, the anchor, the deciding quality) is read off the read's
        // alignment to the haplotype a state's best path spells. That alignment is the path
        // itself: its bases are the haplotype's, its costs sum to the state's score. Random
        // sites of every kind of row, reads with errors, masks and indels of their own; and a
        // 4-base deletion's site in an 18-base repeat beside a listed 14-base deletion, whose
        // read the linear aligner aligns below its state's score (-20 against -17: its band
        // misses the path the graph finds), so a second alignment would examine the call on
        // other bases than its score's.
        let mut st = 83u64;
        let mut seen = 0usize;
        for case in 0..48u64 {
            let c = random_seq(700, 7300 + case);
            let x = match case % 3 {
                0 => Row::new(350, &c[350..351], &[*b"ACGT".iter().find(|&&b| b != c[350]).unwrap()]),
                1 => Row::new(350, &c[350..351], &[&c[350..351], b"TGA"].concat()),
                _ => Row::new(350, &c[350..356], &c[350..351]),
            };
            let mut rows: Vec<Row> = Vec::new();
            for _ in 0..roll(&mut st, 3) {
                let p = if roll(&mut st, 2) == 0 { 300 + roll(&mut st, 40) } else { 362 + roll(&mut st, 40) } as usize;
                let r = random_row(&c, p, &mut st);
                if !rows.contains(&r) && (r.end() < 345 || r.pos > 360) {
                    rows.push(r);
                }
            }
            let mut reads = Vec::new();
            for k in 0..10 {
                let mut on: Vec<&Row> = rows.iter().filter(|_| roll(&mut st, 2) == 0).take(1).collect();
                if k % 2 == 0 && on.iter().all(|o| o.end() <= x.pos || x.end() <= o.pos) {
                    on.push(&x);
                }
                let (mut h, mut o) = applied(&c, &on);
                if roll(&mut st, 3) == 0 {
                    let at = 335 + roll(&mut st, 30) as usize;
                    h.insert(at, h[at]);
                    o.insert(at, None);
                }
                let s = 250 + roll(&mut st, 90) as usize;
                for _ in 0..roll(&mut st, 3) {
                    let i = (s + roll(&mut st, 110) as usize).min(h.len() - 1);
                    h[i] = b"ACGT"[roll(&mut st, 4) as usize];
                }
                let masked: Vec<usize> = (0..roll(&mut st, 3)).map(|_| roll(&mut st, 110) as usize).collect();
                if let Some(r) = read_on(&h, &o, s, 100 + roll(&mut st, 40) as usize, k % 3 == 0, &masked) {
                    reads.push(r);
                }
            }
            let cand = Candidates { listed: rows.clone(), ..Candidates::default() };
            let site = Site::build(x.clone(), &cand, (0, c.clone()), &reads).unwrap();
            for read in &reads {
                for (k, (score, ex, _)) in super::super::decide::examined_alignments(&site, read).into_iter().enumerate() {
                    assert_eq!((ex.0.score, summed(&ex.0)), (score, score), "case {case} state {k}: the examined alignment is the state's path");
                    seen += 1;
                }
            }
        }
        assert!(seen > 1000, "{seen} alignments");
        // the deletion beside a listed deletion, in a repeat
        let mut c = random_seq(700, 31167);
        let unit = random_seq(18, 167);
        for i in 0..60 {
            c[320 + i] = unit[i % 18];
        }
        let x = Row::new(350, &c[350..355], &c[350..351]);
        let listed = Row::new(301, &c[301..316], &c[301..302]);
        let seq = b"TNCANNCTNCGNNTANCANNCTNCGNNTANCCNNCANACNNCGNCTNNCCNCCNNACNACNNTCNTGNNGCNTANNTCNCTNNGCNTGNNCTNTTNNGGNAANNAGNCCNNCGNCCNNGANTCNNAGNACNNGTNG";
        let quals: Vec<u8> = (0..seq.len()).map(|i| [30, 2, 22, 25, 2, 2, 21][i % 7]).collect();
        let cigar = [Cigar::Match(17), Cigar::Del(4), Cigar::Match(1), Cigar::Ins(18), Cigar::Match(83), Cigar::SoftClip(17)];
        let read = ReadView::new(seq, &quals, 334, &cigar, (0, seq.len()), 20).unwrap();
        let cand = Candidates { listed: vec![listed], ..Candidates::default() };
        let site = Site::build(x.clone(), &cand, (0, c.clone()), std::slice::from_ref(&read)).unwrap();
        let both = super::super::decide::examined_alignments(&site, &read);
        let (score, ex, lin) = &both[0];
        assert_eq!((*score, ex.0.score, lin.0.score), (-17, -17, -20), "REF's state: its score, the examined alignment's, the linear aligner's");
    }

    #[test]
    fn on_tied_placements_in_a_run_the_traced_path_is_the_aligners_own_alignment() {
        // A base more in a G run (or two, or a T more in the T run beside it) can be inserted at
        // any of the run's bases at one cost. The graph's traceback and the linear aligner take
        // the same one (both read back from the last base: an aligned base before a deletion
        // before an insertion), on every state's haplotype: these calls are as they were.
        let c = super::tests::contig(5, 344, b"GGGGGGTTTTTT");
        let x = Row::new(350, b"T", b"C");
        for (name, base, extra, at) in [("X, a G more", b'C', b"G".as_slice(), 344usize), ("X, two Gs more", b'C', b"GG", 344), ("REF, a G more", b'T', b"G", 344), ("X, a T more", b'C', b"T", 351)] {
            let mut h = c.clone();
            h[350] = base;
            let seq: Vec<u8> = [&h[290..at], extra, &h[at..410]].concat();
            let cigar = [Cigar::Match((at - 290) as u32), Cigar::Ins(extra.len() as u32), Cigar::Match((410 - at) as u32)];
            let read = ReadView::new(&seq, &vec![30; seq.len()], 290, &cigar, (0, seq.len()), 20).unwrap();
            let site = Site::build(x.clone(), &Candidates::default(), (0, c.clone()), std::slice::from_ref(&read)).unwrap();
            for (k, (_, ex, lin)) in super::super::decide::examined_alignments(&site, &read).into_iter().enumerate() {
                assert_eq!(laid(&ex), laid(&lin), "{name}, state {k}");
            }
        }
    }

    #[test]
    fn a_deletion_tied_across_repeat_units_is_called_on_the_paths_placement() {
        // A period-17 repeat ends 13 bases before X (C>G); the read carries X and one unit
        // deleted, its aligner writing the deletion at one unit's end (71M17D63M). Every
        // unit's end gives the same score. The call's quality is the lowest of the read's
        // confident bases where the winning alignment (X's) and the best losing one (REF's)
        // differ: X, and the bases either side of each alignment's deletion. On the traced
        // paths X's own quality, 25, is the lowest; on the linear aligner's placements (a
        // second alignment, the deletion placed apart on the two haplotypes another way) a
        // base of quality 21 beside a deletion was.
        let mut c = random_seq(700, 31072);
        let unit = random_seq(17, 72);
        for i in 0..60 {
            c[320 + i] = unit[i % 17];
        }
        let x = Row::new(350, &c[350..351], b"G");
        assert_eq!(&c[350..351], b"C");
        let seq = b"AAAGAAGTTGGACTTCGATGCACTCACTCCGGAGATCCCGTTTCCCCGCAGATCCCGTTTCCCGGCAGATCCCGTTTTAAAGCGACTACGGGTGGCCTCTGCTTTGACAAACCTGCACCGAGAACCAGCACCAC";
        let quals: Vec<u8> = (0..seq.len()).map(|i| [25, 30, 21, 30, 30, 30, 25][i % 7]).collect();
        let read = ReadView::new(seq, &quals, 287, &[Cigar::Match(71), Cigar::Del(17), Cigar::Match(63)], (0, seq.len()), 20).unwrap();
        let at_x = read.refpos.iter().position(|&p| p == 350).unwrap();
        let site = Site::build(x.clone(), &Candidates::default(), (0, c.clone()), std::slice::from_ref(&read)).unwrap();
        let v = judge(&site, &read, 20, true);
        assert_eq!((v.call, v.qual, read.quals[at_x]), (super::super::decide::Call::Alt, 25, 25), "ALT at X's own quality");
        let both = super::super::decide::examined_alignments(&site, &read);
        let deleted_after = |a: &Alignment| a.skipped_after.iter().position(|&k| k > 0);
        let path = (deleted_after(&both[0].1 .0), deleted_after(&both[1].1 .0));
        let lin = (deleted_after(&both[0].2 .0), deleted_after(&both[1].2 .0));
        assert!(path != lin, "the deletion's places on REF's and X's haplotypes: paths {path:?}, aligner {lin:?}");
    }
}
