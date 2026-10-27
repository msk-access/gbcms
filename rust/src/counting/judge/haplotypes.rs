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
//! - **The contexts:** sets of the other listed rows within the reads' reach, with the
//!   guard's recurring background changes (they explain reads but are never a site
//!   state). A listed length change whose run X lies inside is also placed at the run's
//!   other end: on the reference the two placements are one allele, beside X two
//!   molecules. With at most [`ALL_SUBSETS`] such rows every set of them that can be
//!   applied together, within [`ALL_SUBSETS_CAP`] haplotypes; with more, the empty set,
//!   then each row, each pair reads show together and each set reads show together, those
//!   explaining more reads first, within [`MAX_HAPS`].
//! - **The haplotypes:** every context crossed with every state, but for a state at a
//!   context row's place, and but for another allele's haplotype that spells one of X's
//!   or REF's (the reads' spelling of X beside a listed neighbour is X's). Every context
//!   appears with REF and with X, so only a read's bases at X's site can tell those
//!   apart.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use rust_htslib::bam::record::Cigar;

use super::graph::{self, Graph, GraphBuilder, GraphRead, ShortReads, Workspace};
use super::score::{self, Alignment, Band, Placed, BAND_SLACK, DELTA, WILD};
use crate::counting::carrier::{GUARD_MIN_READS, GUARD_REF_SHARE};

/// Rows within the reads' reach up to which every applicable set of them is a context.
pub(crate) const ALL_SUBSETS: usize = 4;
/// Haplotypes every applicable set of up to [`ALL_SUBSETS`] rows may fill.
pub(crate) const ALL_SUBSETS_CAP: usize = 64;
/// Haplotypes the contexts of more rows may fill.
pub(crate) const MAX_HAPS: usize = 64;
/// Rows a read's evidence adds to those its bases show on the reference, at most.
const MAX_STEPS: usize = 16;
/// Haplotypes a read's evidence search may score, at most (a site with a few listed rows
/// needs a few dozen).
const MAX_EVALS: usize = 128;
/// Reads whose evidence is searched at one site, at most (spread evenly over its reads).
const MAX_SEARCHED: usize = 400;
/// How far from where a length change can sit a substitution's bases move with it.
const MOVE_NEAR: i64 = 32;
/// At a dense cluster: the distinct sets of rows reads show that are intersected for the
/// sets reads show together (those most reads show)...
const MAX_SHOWN_SETS: usize = 128;
/// ...the sets collected from them...
const MAX_TOGETHER: usize = 1024;
/// ...the molecules the reads assemble into...
const MAX_MOLECULES: usize = 64;
/// ...and the candidate contexts ranked by the reads they explain (those most reads
/// show). Each bound is far above what a site with a few listed rows reaches.
const MAX_RANKED: usize = 512;


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
#[derive(Debug, Clone)]
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
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
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

    /// The genome position its own alignment places base `i` at; None for a base it leaves
    /// unplaced (inserted or clipped).
    pub(crate) fn own_position(&self, i: usize) -> Option<i64> {
        (i < self.seq.len() && self.aligned[i]).then(|| self.refpos[i])
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
            seq.splice(s as usize..e as usize, site.rows[r].alt_al.iter().copied());
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
    pub(crate) haps: Vec<Hap>,
    /// Contexts kept, and whether any candidate context was left out by the cap.
    pub(crate) n_contexts: usize,
    pub(crate) cap_exceeded: bool,
    /// How reads are scored: against `haps` (enumerated) or against `graphs`.
    pub(crate) scorer: ScorerKind,
    /// The graph scorer's graphs, one per state of X's site: REF's bases, X, then the
    /// other states.
    pub(crate) graphs: Vec<StateGraph>,
    /// The coloured k-mer index triage reads first, where it is on and the site allows it.
    pub(crate) triage: Option<super::triage::Index>,
    pub(crate) triage_mode: TriageMode,
    /// What a covered telling base read confidently as neither allele does to a REF or ALT
    /// call ([`neither_rule`]).
    pub(crate) neither: NeitherRule,
    /// At a substitution of two bases or more or a delins X: per genome position over the
    /// window and base, the reads whose own alignments read that base there confidently.
    tally: HashMap<(i64, u8), u32>,
}

/// Which scorer weighs a read (`GBCMS_JUDGE_SCORER`): `graph`, the judge's scorer and the
/// default, one programme per site state over a sequence graph holding every context;
/// `enumerated`, the contexts enumerated and each haplotype scored. TEMPORARY: the
/// enumeration stays behind the switch for one round, as the reference the graph is
/// checked against (identical scores and calls where it is exact, four context rows or
/// fewer). Read when a site is built, so a process may switch it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ScorerKind {
    Graph,
    Enumerated,
}

pub(crate) fn scorer_kind() -> ScorerKind {
    match std::env::var("GBCMS_JUDGE_SCORER").as_deref() {
        Ok("enumerated") => ScorerKind::Enumerated,
        _ => ScorerKind::Graph,
    }
}

/// TEMPORARY, kept as a measured negative result: coloured k-mer triage in front of the
/// graph scorer (`GBCMS_JUDGE_TRIAGE`): `off`, the default; `on`, a read triage decides is
/// not scored; `check`, triage and the scorer both run, every triage-decided read is traced
/// with both calls, and the scorer's call counts. As built it disagreed with the full
/// programme on 0.1-0.6% of the reads it decided; confirmed exactly
/// (`GBCMS_JUDGE_TRIAGE_EXACT`) it disagrees on none but saves only 2-18% of the time, so
/// it stays off.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum TriageMode {
    Off,
    On,
    Check,
}

pub(crate) fn triage_mode() -> TriageMode {
    match std::env::var("GBCMS_JUDGE_TRIAGE").as_deref() {
        Ok("on") => TriageMode::On,
        Ok("check") => TriageMode::Check,
        _ => TriageMode::Off,
    }
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
/// other reads of the row, by their own alignments ([`Site::recurrence`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum NeitherRule {
    Margin,
    Other,
    Below(usize),
    Recur(u32),
}

pub(crate) fn neither_rule() -> NeitherRule {
    let v = std::env::var("GBCMS_JUDGE_NEITHER").unwrap_or_default();
    match v.split_once(':') {
        None if v == "other" => NeitherRule::Other,
        Some(("other", n)) => n.parse().map_or(NeitherRule::Margin, NeitherRule::Below),
        Some(("recur", k)) => k.parse().map_or(NeitherRule::Margin, NeitherRule::Recur),
        _ => NeitherRule::Margin,
    }
}

/// Measurement only, TEMPORARY: with `GBCMS_JUDGE_NEITHER_TRACE=1`, every REF or ALT call on
/// a substitution of two bases or more or a delins traces what the read shows at X's telling
/// bases that is neither allele's (each such base's quality, place and recurrence), and its
/// mismatches elsewhere on the called haplotype ("judge neither" lines).
pub(crate) fn neither_trace() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GBCMS_JUDGE_NEITHER_TRACE").is_ok_and(|v| v == "1"))
}

/// Measurement only, TEMPORARY: where an allele whose equivalent placements end exactly at
/// an insertion X's point stands (`GBCMS_JUDGE_ENDS`): `site`, the default, at X's site (its
/// placements meet X's), another allele there; or `beside`, a change beside X, a context, so
/// a read carrying it is read at X by its bases there. Read when a site is built.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum EndsRule {
    Site,
    Beside,
}

pub(crate) fn ends_rule() -> EndsRule {
    match std::env::var("GBCMS_JUDGE_ENDS").as_deref() {
        Ok("beside") => EndsRule::Beside,
        _ => EndsRule::Site,
    }
}

/// One state of X's site and its graph: the window with the state's row applied (none for
/// REF's bases) and every context row in reach that can lie with it as a bubble.
#[derive(Debug, Clone)]
pub(crate) struct StateGraph {
    pub(crate) class: Class,
    /// The state's row (X is 0), None for REF's bases.
    pub(crate) row: Option<usize>,
    pub(crate) graph: Graph,
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
}

/// Measurement only, TEMPORARY: which of the judge's other-allele options are on
/// (`GBCMS_JUDGE_OTHER`, comma-separated: `runs` the other lengths of X's repeat as
/// states, `found` the indels the reads write at X's site as states (`foundK`: written by
/// K reads, default [`FOUND_MIN_READS`]; `found3` the guard's floor), a [`FoundRule`]
/// token for which of them stay states, `fit` the absolute-fit check, `edges` the check
/// passing over indels inside a long event's own bases; `none`; `both` runs, found and
/// fit). Default `found,fit,told,edges`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct OtherOption {
    pub(crate) runs: bool,
    pub(crate) found: bool,
    pub(crate) fit: bool,
    /// Which found alleles stay states ([`Site::drop_spellings`]).
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
    /// `spell`: dropped where no writer (a read whose alignment writes it) tells it from X,
    /// or none from REF, by score; `spellK`: kept where K writers tell it.
    Spell(usize),
    /// `exact`, `exactK`: as `spell`, a writer telling it from X (REF) only where the bases
    /// it shows of the allele's haplotype (error-free) lie in no haplotype carrying X
    /// (holding REF's bases) as well.
    Exact(usize),
    /// `fragK`: K fragments tell it as `exact` does (a writer or its mate among the site's
    /// reads) and neither mate is told from it by X or REF.
    Frag(usize),
    /// `cleanK` (`clean`: K = 1): the found floor tells it, or K reads show it whole
    /// ([`Pick::Clean`]): reads whose own bases are that allele's, X's and REF's haplotypes
    /// each missing them by a confident base wherever placed. Noise rarely makes a carrier
    /// of X (or a REF read) spell another found allele with no confident edit and no copy
    /// in X's or REF's haplotypes, while one such read is a read carrying another allele.
    Clean(usize),
    /// `share`: one read or more tells it, and they are at least [`GUARD_REF_SHARE`] of the
    /// reads telling it from X either way (those spanning where they differ).
    Share,
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
    let rule = if let Some(k) = numbered("told") {
        FoundRule::Told(k)
    } else if let Some(k) = numbered("spell") {
        FoundRule::Spell(k)
    } else if let Some(k) = numbered("exact") {
        FoundRule::Exact(k)
    } else if let Some(k) = numbered("frag") {
        FoundRule::Frag(k)
    } else if let Some(k) = numbered("clean") {
        FoundRule::Clean(k.max(1))
    } else if on("share") {
        FoundRule::Share
    } else {
        FoundRule::All
    };
    let found = numbered("found");
    let found_min = found.filter(|&k| k > 0).map_or(FOUND_MIN_READS, |k| k as u32);
    OtherOption { runs: on("runs"), found: found.is_some(), fit: on("fit"), rule, edges: on("edges"), found_min }
}

/// Measurement only, TEMPORARY: which sources of context explain a read's *other*
/// differences, so they don't pass for X's allele (`GBCMS_JUDGE_CONTEXTS`). The site's
/// states are the same under every level: REF, X, the listed rows over X's site, an SNV's
/// other bases, and the guard's and the found alleles at X's site. Each level adds one
/// source to the one before it:
/// - `seq`: no context; the candidate sequences differ only at X's site.
/// - `listed`: the listed rows in reach beside X's site, as bubbles.
/// - `nofound`: and the guard's recurring background: its changes, its any-base masks and
///   its recurring indels.
/// - `all` (the default): and the alleles the reads spell away from X's site: the guard's
///   other alleles that change bases beside the site only, and a found allele's changes
///   beside a substitution X.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum ContextSources {
    Seq,
    Listed,
    NoFound,
    All,
}

impl ContextSources {
    fn listed(self) -> bool {
        self >= ContextSources::Listed
    }

    fn background(self) -> bool {
        self >= ContextSources::NoFound
    }

    fn found(self) -> bool {
        self == ContextSources::All
    }
}

pub(crate) fn context_sources() -> ContextSources {
    static SOURCES: std::sync::OnceLock<ContextSources> = std::sync::OnceLock::new();
    *SOURCES.get_or_init(|| {
        let v = std::env::var("GBCMS_JUDGE_CONTEXTS").unwrap_or_default();
        parse_contexts(&v).unwrap_or_else(|| {
            log::warn!("GBCMS_JUDGE_CONTEXTS={v:?} is not seq, listed, nofound or all; every context source is used");
            ContextSources::All
        })
    })
}

fn parse_contexts(v: &str) -> Option<ContextSources> {
    match v {
        "seq" => Some(ContextSources::Seq),
        "listed" => Some(ContextSources::Listed),
        "nofound" => Some(ContextSources::NoFound),
        "all" | "" => Some(ContextSources::All),
        _ => None,
    }
}

/// An event's own bases from which an indel inside them, away from its ends, is a change
/// on its molecule rather than another allele of it (the length from which the rules
/// judge a deletion at its junction and expect no read to hold an allele whole).
const LONG_EVENT: i64 = crate::counting::junction::LONG_DELETION as i64;
/// Bases inside a long event's ends where an indel still meets its junction.
const EDGE: i64 = 1;
/// A long deletion's deleted bases, past its microhomology, that a read over its anchor
/// reads to tell REF from it (RJ-24's REF window: four clean bases of the deleted bases).
const JUNCTION_WINDOW: i64 = 4;

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

/// Measurement only, TEMPORARY: which of the guard's background (its recurring changes,
/// any-base masks and recurring indels beside the event) the judge keeps as contexts
/// (`GBCMS_JUDGE_BACKGROUND`): `asis` (all of it, the default), `respell` (all but the
/// items that respell the site's own candidate sequences, [`respellings`]) or `off`
/// (none, the found alleles kept). Read when a site is built.
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
            log::warn!("GBCMS_JUDGE_BACKGROUND={v:?} is not asis, respell or off; the background is kept as it is");
            BackgroundRule::AsIs
        })
    })
}

fn parse_background(v: &str) -> Option<BackgroundRule> {
    match v {
        "asis" | "" => Some(BackgroundRule::AsIs),
        "respell" => Some(BackgroundRule::Respell),
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
/// it spelling REF).
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
            let back = if one_way { 0 } else { unmakes(&refs, None, &alts, &alt_b) };
            out_rows[i] = out_rows[i].max(unmakes(&alts, None, &refs, &ref_b)).max(back);
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

    fn independent(&self, set: &[usize]) -> bool {
        (0..set.len()).all(|i| (i + 1..set.len()).all(|j| !self.conflict(set[i], set[j])))
    }

    fn distance_to_target(&self, r: usize) -> i64 {
        let (s, e) = self.spans[r];
        let (xs, xe) = self.spans[0];
        if s < xe && xs < e {
            0
        } else if e <= xs {
            xs - e + 1
        } else {
            s - xe + 1
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
    /// deletion (its tract is read by the scores) and for another allele's haplotype.
    pub(crate) fn telling_in(&self, hap: &Hap) -> Vec<(i64, u8)> {
        let x = &self.rows[0];
        if x.changes_length() && super::super::window::is_pure_indel(
            &String::from_utf8_lossy(&x.ref_al),
            &String::from_utf8_lossy(&x.alt_al),
        ) {
            return Vec::new();
        }
        let (on_ref, on_alt) = telling(x);
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

    /// X's first base, 0-based.
    pub(crate) fn x_pos(&self) -> i64 {
        self.rows[0].pos
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
    /// changes beside X's read confidently as the state's, not only X's part of it.
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
        let (xs, xe) = self.spans[0];
        let spells = |r: usize| {
            let row = &self.rows[r];
            let (s, _) = self.spans[r];
            (0..row.ref_al.len()).filter(|&i| row.ref_al[i] != row.alt_al[i] && !(xs..xe).contains(&(s + i as i64))).all(|i| {
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
        Site::build_with(scorer_kind(), target, candidates, window, reads)
    }

    /// [`Site::build`] with the scorer named.
    pub(crate) fn build_with(
        scorer: ScorerKind,
        target: Row,
        candidates: &Candidates,
        window: (i64, Vec<u8>),
        reads: &[ReadView],
    ) -> Option<Site> {
        Site::build_ruled(scorer, target, candidates, window, reads, other_option().rule)
    }

    /// [`Site::build`] with the scorer named, the found alleles kept as states by `rule`.
    fn build_ruled(
        scorer: ScorerKind,
        target: Row,
        candidates: &Candidates,
        window: (i64, Vec<u8>),
        reads: &[ReadView],
        rule: FoundRule,
    ) -> Option<Site> {
        Site::build_sourced(scorer, target, candidates, window, reads, rule, context_sources())
    }

    /// [`Site::build_ruled`] with the context sources named.
    fn build_sourced(
        scorer: ScorerKind,
        target: Row,
        candidates: &Candidates,
        window: (i64, Vec<u8>),
        reads: &[ReadView],
        rule: FoundRule,
        sources: ContextSources,
    ) -> Option<Site> {
        Site::build_backed(scorer, target, candidates, window, reads, rule, sources, background_rule())
    }

    /// [`Site::build_sourced`] with the guard's background kept by `bg_rule`.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn build_backed(
        scorer: ScorerKind,
        target: Row,
        candidates: &Candidates,
        window: (i64, Vec<u8>),
        reads: &[ReadView],
        rule: FoundRule,
        sources: ContextSources,
        bg_rule: BackgroundRule,
    ) -> Option<Site> {
        Site::build_core(scorer, target, candidates, window, reads, rule, sources, bg_rule, None, ends_rule())
    }

    /// [`Site::build`] with an allele ending at an insertion X's point placed by `ends`.
    #[cfg(test)]
    pub(crate) fn build_ended(target: Row, candidates: &Candidates, window: (i64, Vec<u8>), reads: &[ReadView], ends: EndsRule) -> Option<Site> {
        let (rule, sources, bg) = (other_option().rule, context_sources(), background_rule());
        Site::build_core(scorer_kind(), target, candidates, window, reads, rule, sources, bg, None, ends)
    }

    /// [`Site::build_backed`], the context items in `dropped` left out where given (the
    /// respelling check's own builds).
    #[allow(clippy::too_many_arguments)]
    fn build_core(
        scorer: ScorerKind,
        target: Row,
        candidates: &Candidates,
        (start, reference): (i64, Vec<u8>),
        reads: &[ReadView],
        rule: FoundRule,
        sources: ContextSources,
        bg_rule: BackgroundRule,
        dropped: Option<&Dropped>,
        ends: EndsRule,
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
        // Under `GBCMS_JUDGE_ENDS=beside`: an allele changing the length whose equivalent
        // placements end exactly at an insertion X's point is beside X, a context.
        let ends_at_x = |r: &Row| {
            ends == EndsRule::Beside
                && pure(&target)
                && target.alt_al.len() > target.ref_al.len()
                && r.changes_length()
                && *r != target
                && place(r).1 == target.end()
        };
        for r in candidates.listed.iter().map(beside) {
            let state = meets(&target, &r) && !in_deleted(&r) && !ends_at_x(&r);
            if state || sources.listed() {
                take(&r, state);
            }
        }
        if target.ref_al.len() == 1 && target.alt_al.len() == 1 {
            // An SNV's site holds every base (a third base is another allele).
            for &b in b"ACGT" {
                take(&Row::new(target.pos, &target.ref_al, &[b]), true);
            }
        }
        // The contexts the reads' own alleles supply (the found source), each once.
        let mut found_contexts: Vec<Row> = Vec::new();
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
                    if sources.found() && take(&minimal, false) {
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
            for (r, by) in candidates.found.iter().map(|(r, by)| (beside(r), &by[..])).filter(|(r, _)| meets(&target, r)) {
                if ends_at_x(&r) {
                    if sources.found() && take(&r, false) {
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
                        if !meets(&target, &c) && sources.found() && take(&c, false) {
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
        // A background change at X's place cannot share a haplotype with X: crossed with
        // REF alone it would favour REF, so it is left to the reads' bases there. Under
        // `GBCMS_JUDGE_BACKGROUND=respell` an item of the guard's background or of the found
        // contexts that respells the site's own candidate sequences ([`respellings`]) is
        // left out, change or mask.
        let bg_rule = if sources.background() { bg_rule } else { BackgroundRule::Off };
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
            if partial.is_empty() || scorer != ScorerKind::Graph {
                return whole;
            }
            let build = |d: &Dropped| {
                Site::build_core(scorer, target.clone(), candidates, (start, reference.clone()), reads, rule, sources, bg_rule, Some(d), ends)
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
                settle(whole, items.iter().zip(&rows).filter(|(_, &k)| k == 1).map(|(r, _)| r).collect())
            }
            // Under the default, at a substitution X, a recurring indel of the guard's with
            // which REF reads as X over a read reaching X's site is X's change respelled (a
            // repeat X lengthens, written as more units of it beside X: AGCA>TGCT between a
            // GCT run and a GCA run, REF with two GCT units more reads as X into the GCA
            // run), never a context: left out where no reads need it for their call. A
            // substitution's carriers hold its bases at X; an indel there is a shift of the
            // repeat beside it. The masks, the substitutions and the found alleles stay
            // contexts, and an indel X's site keeps its background as it is.
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
                settle(whole, indels.iter().zip(&rows).filter(|(_, &k)| k == 1).map(|(r, _)| r).collect())
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
        // same length over X's bases that reads X's bases where X changes them.
        let holds_x: Vec<bool> = {
            let x = &rows[0];
            let (on_ref, _) = telling(x);
            rows.iter()
                .enumerate()
                .map(|(r, row)| {
                    (1..=n_states).contains(&r)
                        && !x.changes_length()
                        && !row.changes_length()
                        && row.pos <= x.pos
                        && row.end() >= x.end()
                        && on_ref.clone().all(|k| row.alt_al[(x.pos + k as i64 - row.pos) as usize] == x.alt_al[k])
                })
                .collect()
        };
        let spans: Vec<(i64, i64)> = rows.iter().map(|r| (r.pos - start, r.end() - start)).collect();
        let places = rows.iter().map(place).map(|(s, e)| (s - start, e - start)).collect();
        let max_shift = rows.iter().map(|r| r.ref_al.len().abs_diff(r.alt_al.len()) as i64).sum();
        // A mask where a row changes a base is that row's base: the row explains the reads
        // carrying it (a context with every state, or a state of the site). A row spelled
        // over a wider window (an allele the guard found) changes only some of its bases.
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
            .filter(|&m| m >= 0 && m < wn && !changed.iter().any(|&(s, e)| s <= m && m < e))
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
            haps: Vec::new(),
            n_contexts: 0,
            cap_exceeded: false,
            scorer,
            graphs: Vec::new(),
            triage: None,
            triage_mode: triage_mode(),
            neither: neither_rule(),
            tally,
        };

        // The context rows a read over X can reach: REF span meeting X's ± the longest read
        // (a long deletion starting there runs past it, yet a read over X can carry its
        // junction).
        let l = reads.iter().map(|r| r.len() as i64).max().unwrap_or(0);
        let (xs, xe) = site.spans[0];
        let rest: Vec<usize> = (1 + site.n_states..site.rows.len())
            .filter(|&r| site.spans[r].1 > xs - l && site.spans[r].0 < xe + l)
            .collect();
        if scorer == ScorerKind::Graph {
            // Every context at once: one graph per state, no enumeration.
            site.graphs = site.state_graphs(&rest);
            site.drop_spelled_states(&rest);
            site.drop_graph_spellings(&found_states, reads, rule);
            if site.triage_mode != TriageMode::Off {
                site.triage = super::triage::Index::build(&site, l);
            }
            return Some(site);
        }
        let (ordered, cap) = if rest.len() <= ALL_SUBSETS {
            (site.subsets(&rest), ALL_SUBSETS_CAP)
        } else {
            let evidence = site.evidence(reads, &rest);
            (site.ranked_contexts(&rest, &evidence), MAX_HAPS)
        };
        let (sets, n_contexts, cap_exceeded) = site.close_over_site(ordered, cap);
        let haps: Vec<Hap> = sets.iter().map(|rows| Hap::build(&site, rows)).collect();
        // Another allele that spells X or REF in one of the contexts is that haplotype (an
        // allele the reads spell over X with a listed neighbour beside it): X's or REF's,
        // never another allele's. A haplotype its class already holds is scored once.
        let given: BTreeSet<&[u8]> = haps.iter().filter(|h| h.class != Class::Other).map(|h| &h.seq[..]).collect();
        let mut held: BTreeSet<(usize, &[u8])> = BTreeSet::new();
        let keep: Vec<bool> = haps
            .iter()
            .map(|h| (h.class != Class::Other || !given.contains(&h.seq[..])) && held.insert((h.class as usize, &h.seq[..])))
            .collect();
        site.haps = haps.into_iter().zip(keep).filter_map(|(h, k)| k.then_some(h)).collect();
        site.n_contexts = n_contexts;
        site.cap_exceeded = cap_exceeded;
        site.drop_spellings(&found_states, reads, rule);
        Some(site)
    }

    /// Every set of `rows` that can be applied together, the empty one first, then by
    /// size and rows.
    fn subsets(&self, rows: &[usize]) -> Vec<Vec<usize>> {
        let mut sets: Vec<Vec<usize>> = (0u32..1 << rows.len())
            .map(|m| (0..rows.len()).filter(|i| m >> i & 1 == 1).map(|i| rows[i]).collect::<Vec<_>>())
            .filter(|h| self.independent(h))
            .collect();
        sets.sort_by(|a, b| a.len().cmp(&b.len()).then(a.cmp(b)));
        sets
    }

    /// Per read, the context rows it shows and those it refutes, among `rows` it
    /// reaches. A substitution is read off the read's alignment to the reference: shown
    /// where every base it changes is read confidently as its ALT's, refuted where every
    /// one is read as REF's. The length changes are chosen together, one at a time on
    /// the haplotype of those already shown (the substitutions it shows first): each
    /// time the one that fits the read best, while it fits `DELTA` better than without
    /// it; a length change that then fits it `DELTA` worse is refuted. So two rows that
    /// each explain part of one molecule's indel (two insertions in one repeat) are not
    /// both shown by a read carrying one of them, and a read carrying several changes
    /// shows the set its molecule carries.
    fn evidence(&self, reads: &[ReadView], rows: &[usize]) -> Vec<Evidence> {
        let reference = Hap::build(self, &[]);
        let mut built: HashMap<Vec<usize>, Hap> = HashMap::new();
        let wn = self.reference.len() as i64;
        // At a deep site the search runs on reads spread evenly over it, at most
        // MAX_SEARCHED; the rest keep what their bases show on the reference.
        let stride = reads.len().div_ceil(MAX_SEARCHED).max(1);
        reads
            .iter()
            .enumerate()
            .map(|(i, read)| {
                let mut scorer = Scorer::new(read);
                let l = read.len() as i64;
                let (a, b) = ((read.ref_start - self.start - l).clamp(0, wn), (read.ref_end - self.start + l).clamp(0, wn));
                // the reference column each read base is aligned to
                let (on_ref, off) = scorer.align(self, &reference);
                let mut at: HashMap<i64, u8> = HashMap::new();
                for (i, p) in on_ref.placed.iter().enumerate() {
                    if let Placed::At(c) = *p {
                        at.insert((c + off) as i64, read.seq[i]);
                    }
                }
                let mut e = Evidence::default();
                let mut lengths = Vec::new();
                for &r in rows {
                    let (rs, re) = self.spans[r];
                    if rs < a || re > b {
                        continue;
                    }
                    let row = &self.rows[r];
                    if row.changes_length() {
                        lengths.push(r);
                        continue;
                    }
                    let reads: Vec<Option<u8>> = (0..row.ref_al.len())
                        .filter(|&k| row.ref_al[k] != row.alt_al[k])
                        .map(|k| {
                            at.get(&(rs + k as i64)).copied().filter(|&b| b != WILD).map(|b| {
                                if b == row.alt_al[k] {
                                    1
                                } else if b == row.ref_al[k] {
                                    0
                                } else {
                                    2
                                }
                            })
                        })
                        .collect();
                    if reads.iter().all(|v| *v == Some(1)) {
                        e.shown.push(r);
                    } else if reads.iter().all(|v| *v == Some(0)) {
                        e.refuted.push(r);
                    }
                }
                if !lengths.is_empty() && i % stride == 0 {
                    // Seeded with the substitutions shown that can be applied together, the
                    // rows the read's molecule carries are added one at a time: each time the
                    // one that fits the read best, while it fits DELTA better. A length change
                    // moves the read's bases, so once one is chosen the substitutions are
                    // weighed by the fit too (their bases read off the reference may sit
                    // elsewhere on the molecule).
                    let mut chosen: Vec<usize> = Vec::new();
                    for &r in &e.shown {
                        if chosen.iter().all(|&c| !self.conflict(c, r)) {
                            chosen.push(r);
                        }
                    }
                    // A set's fit with REF at X's site, or with X (where X can lie with it):
                    // the molecule a read carries is chosen whichever allele it holds there,
                    // each allele's search apart (a read carrying one molecule's indels can
                    // fit another molecule with X almost as well, step by step).
                    let mut score_of = |set: &[usize], with_x: bool, scorer: &mut Scorer| -> i32 {
                        let mut key = set.to_vec();
                        if with_x {
                            if key.iter().any(|&c| self.conflict(0, c)) {
                                return i32::MIN / 4;
                            }
                            key.push(0);
                        }
                        key.sort_unstable();
                        let hap = built.entry(key.clone()).or_insert_with(|| Hap::build(self, &key));
                        scorer.score(self, hap)
                    };
                    let reachable: Vec<usize> =
                        rows.iter().copied().filter(|&r| self.spans[r].0 >= a && self.spans[r].1 <= b).collect();
                    // From the seed, and from the seed with each length change first (one
                    // molecule's indels can each fit a read less well alone than another
                    // molecule's one does), the best set reached.
                    let seed = chosen.clone();
                    // A substitution near a length change (within MOVE_NEAR of where it can
                    // sit) is read anew with it: the change moves the read's bases there.
                    let near = |r: usize, set: &[usize]| {
                        let (s, e) = self.spans[r];
                        set.iter().filter(|&&c| self.rows[c].changes_length()).any(|&c| {
                            let (ps, pe) = self.places[c];
                            s < pe + MOVE_NEAR && ps - MOVE_NEAR < e
                        })
                    };
                    let mut evals = 0usize;
                    let mut score_counted = |set: &[usize], with_x: bool, scorer: &mut Scorer, evals: &mut usize| {
                        *evals += 1;
                        score_of(set, with_x, scorer)
                    };
                    // the seed fitting the read whole with REF at X: nothing fits it better
                    let seed_fit = score_counted(&seed, false, &mut scorer, &mut evals);
                    let starts: Vec<Option<usize>> = if seed_fit == 0 {
                        Vec::new()
                    } else {
                        std::iter::once(None).chain(lengths.iter().map(|&r| Some(r))).collect()
                    };
                    let mut reached: Option<(i32, Vec<usize>)> = Some((seed_fit, seed.clone()));
                    let mut reached_x = false;
                    'starts: for (with_x, first) in [false, true].into_iter().flat_map(|x| starts.iter().map(move |&f| (x, f))) {
                        // the seed but for the rows at the first change's place
                        let mut set: Vec<usize> =
                            seed.iter().copied().filter(|&c| first.is_none_or(|f| !self.conflict(c, f))).collect();
                        set.extend(first);
                        let mut current = score_counted(&set, with_x, &mut scorer, &mut evals);
                        for _ in 0..MAX_STEPS {
                            if evals >= MAX_EVALS {
                                break 'starts;
                            }
                            // the best move: a row added (a length change, or a substitution
                            // near one chosen), or such a substitution dropped
                            let mut best: Option<(i32, Vec<usize>)> = None;
                            for &r in reachable.iter().filter(|&&r| {
                                !set.contains(&r)
                                    && (self.rows[r].changes_length() || near(r, &set))
                                    && set.iter().all(|&c| !self.conflict(c, r))
                            }) {
                                let next = [&set[..], &[r]].concat();
                                let s = score_counted(&next, with_x, &mut scorer, &mut evals);
                                if s - current >= DELTA && best.as_ref().is_none_or(|(b, _)| s > *b) {
                                    best = Some((s, next));
                                }
                            }
                            for &r in set.iter().filter(|&&r| !self.rows[r].changes_length() && near(r, &set)) {
                                let next: Vec<usize> = set.iter().copied().filter(|&c| c != r).collect();
                                let s = score_counted(&next, with_x, &mut scorer, &mut evals);
                                if s - current >= DELTA && best.as_ref().is_none_or(|(b, _)| s > *b) {
                                    best = Some((s, next));
                                }
                            }
                            let Some((s, next)) = best else { break };
                            set = next;
                            current = s;
                        }
                        if reached.as_ref().is_none_or(|(s, r)| current > *s || (current == *s && set.len() < r.len())) {
                            reached = Some((current, set));
                            reached_x = with_x;
                        }
                    }
                    let (current, chosen) = reached.unwrap_or((seed_fit, seed));
                    // A substitution near a chosen length change has the evidence the fit
                    // gives it; every other keeps what its bases read off the reference show.
                    e.shown.retain(|&r| !near(r, &chosen) || chosen.contains(&r));
                    e.refuted.retain(|&r| !near(r, &chosen));
                    for &r in &reachable {
                        if chosen.contains(&r) {
                            if !e.shown.contains(&r) {
                                e.shown.push(r);
                            }
                        } else if (self.rows[r].changes_length() || near(r, &chosen))
                            && !e.refuted.contains(&r)
                            && chosen.iter().all(|&c| !self.conflict(c, r))
                            && evals < 2 * MAX_EVALS
                            && current - score_counted(&[&chosen[..], &[r]].concat(), reached_x, &mut scorer, &mut evals) >= DELTA
                        {
                            e.refuted.push(r);
                        }
                    }
                }
                e.shown.sort_unstable();
                e.refuted.sort_unstable();
                e
            })
            .collect()
    }

    /// Contexts from more than [`ALL_SUBSETS`] rows, in the order they are kept while
    /// there is room: each row, each pair two reads or more show together, each set reads
    /// show together and each molecule the reads assemble into, ranked by the reads they
    /// explain (a context missing a row the reads carry with another leaves them fitting
    /// no haplotype, and the state that best stands in for the row then wins them), then
    /// by fewer rows and nearer X. At a dense cluster (rows on nearly every base, deep
    /// reads) the enumeration is bounded: the [`MAX_RANKED`] candidates most reads show
    /// are ranked (every candidate where there are no more).
    fn ranked_contexts(&self, rest: &[usize], evidence: &[Evidence]) -> Vec<Vec<usize>> {
        let conflicts = self.conflicts_among(rest);
        // Per row, the reads showing it (a bitset): a set two reads show together holds
        // only rows two reads or more show.
        let words = evidence.len().div_ceil(64);
        let mut shown_by: HashMap<usize, Vec<u64>> = HashMap::new();
        for (i, e) in evidence.iter().enumerate() {
            for &r in &e.shown {
                shown_by.entry(r).or_insert_with(|| vec![0u64; words])[i / 64] |= 1 << (i % 64);
            }
        }
        let count = |bits: &[u64]| bits.iter().map(|w| w.count_ones() as usize).sum::<usize>();
        let support: HashMap<usize, usize> = shown_by.iter().map(|(&r, b)| (r, count(b))).collect();
        // (reads showing it, the set)
        let mut candidates: BTreeMap<Vec<usize>, usize> = BTreeMap::new();
        let add = |set: Vec<usize>, n: usize, candidates: &mut BTreeMap<Vec<usize>, usize>| {
            let e = candidates.entry(set).or_insert(0);
            *e = (*e).max(n);
        };
        for &r in rest {
            add(vec![r], support.get(&r).copied().unwrap_or(0), &mut candidates);
        }
        let paired: Vec<usize> = rest.iter().copied().filter(|r| support.get(r).is_some_and(|&n| n >= 2)).collect();
        for (i, &a) in paired.iter().enumerate() {
            for &b in &paired[i + 1..] {
                if conflicts[&a].contains(&b) {
                    continue;
                }
                let together: usize = shown_by[&a].iter().zip(&shown_by[&b]).map(|(x, y)| (x & y).count_ones() as usize).sum();
                if together >= 2 {
                    add(vec![a, b], together, &mut candidates);
                }
            }
        }
        for (set, n) in self.shown_together(evidence, &support) {
            add(set, n, &mut candidates);
        }
        for (set, n) in self.molecules(evidence, &conflicts) {
            add(set, n, &mut candidates);
        }
        let mut pool: Vec<(Vec<usize>, usize)> = candidates.into_iter().collect();
        if pool.len() > MAX_RANKED {
            pool.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.len().cmp(&b.0.len())).then(a.0.cmp(&b.0)));
            pool.truncate(MAX_RANKED);
        }
        self.rank(evidence, pool.into_iter().map(|(set, _)| set).collect(), &conflicts)
    }

    /// Per row of `rows`, the rows of `rows` at its place.
    fn conflicts_among(&self, rows: &[usize]) -> HashMap<usize, Vec<usize>> {
        rows.iter()
            .map(|&a| (a, rows.iter().copied().filter(|&b| b != a && self.conflict(a, b)).collect()))
            .collect()
    }

    /// The molecules the reads assemble into: reads whose evidence agrees (none shows a
    /// row the others refute or that sits at the place of one they show, nor refutes one
    /// they show) join, the rows they show together one context, so a haplotype whose
    /// rows no one read spans whole is still a context. Larger reads' sets first; each
    /// read joins the first molecule it agrees with (a read agreeing with none once
    /// [`MAX_MOLECULES`] are open joins none). Those of two reads or more, with their reads.
    fn molecules(&self, evidence: &[Evidence], conflicts: &HashMap<usize, Vec<usize>>) -> Vec<(Vec<usize>, usize)> {
        let mut order: Vec<&Evidence> = evidence.iter().filter(|e| !e.shown.is_empty()).collect();
        order.sort_by(|a, b| b.shown.len().cmp(&a.shown.len()).then(a.shown.cmp(&b.shown)).then(a.refuted.cmp(&b.refuted)));
        // per molecule: the rows shown, those refuted, and its reads
        let mut found: Vec<(BTreeSet<usize>, BTreeSet<usize>, usize)> = Vec::new();
        for e in order {
            let agrees = |(shown, refuted, _): &(BTreeSet<usize>, BTreeSet<usize>, usize)| {
                e.shown.iter().all(|r| {
                    !refuted.contains(r) && conflicts.get(r).is_none_or(|cs| cs.iter().all(|c| !shown.contains(c)))
                }) && e.refuted.iter().all(|r| !shown.contains(r))
            };
            let open = found.len() < MAX_MOLECULES;
            match found.iter_mut().find(|m| agrees(m)) {
                Some(m) => {
                    m.0.extend(e.shown.iter().copied());
                    m.1.extend(e.refuted.iter().copied());
                    m.2 += 1;
                }
                None if open => found.push((e.shown.iter().copied().collect(), e.refuted.iter().copied().collect(), 1)),
                None => {}
            }
        }
        // A read that agrees with several molecules (one showing only rows they share,
        // germline rows on both haplotypes) carries its rows to each of them.
        let mut found: Vec<_> = found.into_iter().filter(|m| m.2 >= 2).collect();
        for e in evidence.iter().filter(|e| !e.shown.is_empty()) {
            for m in found.iter_mut() {
                let agrees = e.shown.iter().all(|r| {
                    !m.1.contains(r) && conflicts.get(r).is_none_or(|cs| cs.iter().all(|c| !m.0.contains(c)))
                }) && e.refuted.iter().all(|r| !m.0.contains(r));
                if agrees {
                    m.0.extend(e.shown.iter().copied());
                }
            }
        }
        found.into_iter().map(|m| (m.0.into_iter().collect(), m.2)).collect()
    }

    /// Every maximal set of rows that can be applied together shown together by two
    /// reads or more, with the reads of the two sets it came from. Bounded at a dense
    /// cluster: only rows two reads or more show take part (no other can be in two reads'
    /// sets); the [`MAX_SHOWN_SETS`] distinct sets most reads show are intersected, and
    /// the collection stops at [`MAX_TOGETHER`] sets.
    fn shown_together(&self, evidence: &[Evidence], support: &HashMap<usize, usize>) -> Vec<(Vec<usize>, usize)> {
        let mut distinct: HashMap<Vec<usize>, usize> = HashMap::new();
        for e in evidence {
            let set: Vec<usize> = e.shown.iter().copied().filter(|r| support.get(r).is_some_and(|&n| n >= 2)).collect();
            if !set.is_empty() {
                *distinct.entry(set).or_insert(0) += 1;
            }
        }
        let mut distinct: Vec<(Vec<usize>, usize)> = distinct.into_iter().collect();
        distinct.sort_by(|a, b| b.1.cmp(&a.1).then(b.0.len().cmp(&a.0.len())).then(a.0.cmp(&b.0)));
        distinct.truncate(MAX_SHOWN_SETS);
        let mut together: BTreeMap<Vec<usize>, usize> = BTreeMap::new();
        let mut seen: std::collections::HashSet<Vec<usize>> = std::collections::HashSet::new();
        'pairs: for a in 0..distinct.len() {
            for b in a..distinct.len() {
                if a == b && distinct[a].1 < 2 {
                    continue;
                }
                let inter: Vec<usize> =
                    distinct[a].0.iter().copied().filter(|r| distinct[b].0.binary_search(r).is_ok()).collect();
                // Many read pairs share an intersection: enumerate each one once.
                if !inter.is_empty() && seen.insert(inter.clone()) {
                    let mut sets = BTreeSet::new();
                    self.maximal_independent_sets(&inter, &mut sets);
                    let n = if a == b { distinct[a].1 } else { distinct[a].1 + distinct[b].1 };
                    for s in sets {
                        let e = together.entry(s).or_insert(0);
                        *e = (*e).max(n);
                    }
                    if together.len() >= MAX_TOGETHER {
                        break 'pairs;
                    }
                }
            }
        }
        // the maximal ones: none inside a longer one
        let mut by_len: Vec<(Vec<usize>, usize)> = together.into_iter().collect();
        by_len.sort_by(|a, b| b.0.len().cmp(&a.0.len()).then(a.0.cmp(&b.0)));
        let mut kept: Vec<(Vec<usize>, usize)> = Vec::new();
        for (t, n) in by_len {
            if !kept.iter().any(|(u, _)| u.len() > t.len() && is_subset(&t, u)) {
                kept.push((t, n));
            }
        }
        kept
    }

    /// The maximal sets of `set` (sorted) that can be applied together, bounded.
    fn maximal_independent_sets(&self, set: &[usize], out: &mut BTreeSet<Vec<usize>>) {
        if self.independent(set) {
            out.insert(set.to_vec());
            return;
        }
        const MAX_SETS: usize = 64;
        let mut budget = 20_000usize;
        let mut found: Vec<Vec<usize>> = Vec::new();
        fn rec(
            site: &Site,
            set: &[usize],
            i: usize,
            chosen: &mut Vec<usize>,
            found: &mut Vec<Vec<usize>>,
            budget: &mut usize,
        ) {
            if *budget == 0 || found.len() >= MAX_SETS {
                return;
            }
            *budget -= 1;
            if i == set.len() {
                // maximal: every row left out is at a chosen one's place
                if set.iter().all(|r| chosen.contains(r) || chosen.iter().any(|&c| site.conflict(c, *r))) {
                    found.push(chosen.clone());
                }
                return;
            }
            let r = set[i];
            if chosen.iter().all(|&c| !site.conflict(c, r)) {
                chosen.push(r);
                rec(site, set, i + 1, chosen, found, budget);
                chosen.pop();
            }
            rec(site, set, i + 1, chosen, found, budget);
        }
        rec(self, set, 0, &mut Vec::new(), &mut found, &mut budget);
        out.extend(found);
    }

    /// A context explains a read when it holds none of the rows the read refutes, at
    /// least one it shows, and accounts for every one it shows (holding it, or a row at
    /// its place: rows at one place can all beat the reference in one read). `h` sorted.
    fn explains(&self, h: &[usize], e: &Evidence, conflicts: &HashMap<usize, Vec<usize>>) -> bool {
        let holds = |r: &usize| h.binary_search(r).is_ok();
        !e.refuted.iter().any(holds)
            && e.shown.iter().any(holds)
            && e.shown.iter().all(|s| holds(s) || conflicts.get(s).is_some_and(|cs| cs.iter().any(holds)))
    }

    /// Candidates ranked: reads explained (more first), fewer rows, nearer X, positions.
    fn rank(&self, evidence: &[Evidence], cands: Vec<Vec<usize>>, conflicts: &HashMap<usize, Vec<usize>>) -> Vec<Vec<usize>> {
        let mut keyed: Vec<(usize, usize, i64, Vec<usize>)> = cands
            .into_iter()
            .map(|h| {
                let support = evidence.iter().filter(|e| self.explains(&h, e, conflicts)).count();
                let d = h.iter().map(|&r| self.distance_to_target(r)).min().unwrap_or(0);
                (support, h.len(), d, h)
            })
            .collect();
        keyed.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)).then(a.3.cmp(&b.3)));
        keyed.into_iter().map(|k| k.3).collect()
    }

    /// The bases a read shows of `hap`: the haplotype's bases its alignment there covers,
    /// from the first base it places to the last, any base where the read's own is
    /// unreadable. What the read would show, were it from `hap` without an error.
    fn shown(&self, read: &ReadView, hap: &Hap) -> Vec<u8> {
        let (al, off) = Scorer::new(read).align(self, hap);
        let Some((a, b)) = al.covers() else {
            return Vec::new();
        };
        let mut shown = hap.seq[off + a..=off + b].to_vec();
        for (i, p) in al.placed.iter().enumerate() {
            if let Placed::At(c) = *p {
                if read.seq[i] == WILD {
                    shown[c - a] = WILD;
                }
            }
        }
        shown
    }

    /// One read's best scores on the haplotypes `holds` marks (a found allele's), on those
    /// carrying X and on those holding REF's bases, with the first haplotype reaching the
    /// allele's best.
    fn telling(&self, read: &ReadView, holds: &[bool]) -> Telling {
        let mut scorer = Scorer::new(read);
        let none = i32::MIN / 4;
        let (mut f, mut h, mut x, mut r) = (none, 0usize, none, none);
        for (k, (hap, &held)) in self.haps.iter().zip(holds).enumerate() {
            let s = scorer.score(self, hap);
            if held && s > f {
                (f, h) = (s, k);
            }
            match hap.class {
                Class::Alt => x = x.max(s),
                Class::Ref => r = r.max(s),
                Class::Other => {}
            }
        }
        Telling { f, x, r, h, lies: None }
    }

    /// Whether the bases `read` shows of haplotype `h` ([`Site::shown`]) lie whole in a
    /// haplotype carrying X, and in one holding REF's bases: wherever placed there (a copy
    /// of a repeat the scorer's band, which follows the read's own alignment, misses).
    fn lies_in_x_or_ref(&self, read: &ReadView, h: usize) -> (bool, bool) {
        let shown = self.shown(read, &self.haps[h]);
        let lies = |class: Class| self.haps.iter().filter(|hap| hap.class == class).any(|hap| holds_bases(&hap.seq, &shown));
        (lies(Class::Alt), lies(Class::Ref))
    }

    /// A found allele (an indel the reads' own alignments write at X's site) stays a state
    /// only where `rule` ([`FoundRule`]) finds it another allele rather than X's or REF's
    /// spelling. A read ending inside X's own bases (a long insertion's first bases, one
    /// copy of a deleted duplication, a repeat X lengthens) writes the part of X it holds
    /// as an indel that reads, over its span, as X does: such an allele is X's (or REF's)
    /// spelling, never another allele, and its haplotypes are dropped. A sibling allele
    /// whose few reads spell it whole (most of its carriers ending in a repeat or clipping
    /// its end) is told apart by many: another allele, a state.
    fn drop_spellings(&mut self, found: &[(Row, &[usize])], reads: &[ReadView], rule: FoundRule) {
        // Measurement only, TEMPORARY: with GBCMS_JUDGE_FOUND_PROBE=1 every rule's decision
        // is logged per found allele (the rule in force still decides).
        let probe = std::env::var("GBCMS_JUDGE_FOUND_PROBE").is_ok_and(|v| v == "1");
        if rule == FoundRule::All && !probe {
            return;
        }
        let mut dropped: Vec<usize> = Vec::new();
        for (row, by) in found {
            let Some(o) = self.rows.iter().position(|r| r == row).filter(|&o| self.is_state(o)) else {
                continue;
            };
            let holds: Vec<bool> = self.haps.iter().map(|h| h.rows.binary_search(&o).is_ok()).collect();
            if !holds.contains(&true) {
                continue;
            }
            let mut tally = Tally::new(self, reads, holds, by);
            let keep = tally.decide(rule);
            if probe {
                let rules = [
                    ("told", FoundRule::Told(0)),
                    ("told1", FoundRule::Told(1)),
                    ("spell", FoundRule::Spell(0)),
                    ("spell1", FoundRule::Spell(1)),
                    ("exact1", FoundRule::Exact(1)),
                    ("frag1", FoundRule::Frag(1)),
                    ("share", FoundRule::Share),
                    ("clean", FoundRule::Clean(1)),
                ];
                let verdicts: Vec<String> =
                    rules.iter().map(|&(name, r)| format!("{name}={}", u8::from(tally.decide(r)))).collect();
                let n = |tally: &mut Tally, pick: Pick| (0..reads.len()).filter(|&i| tally.picks(i, pick)).count();
                let (tellers, clean, spanning) = (n(&mut tally, Pick::Tells), n(&mut tally, Pick::Clean), n(&mut tally, Pick::SpansX));
                let w: Vec<Telling> = by.iter().map(|&i| tally.get(i)).collect();
                log::debug!(
                    "haplotype judge: found-allele probe {}:{}>{} at {}:{}>{} writers {} of {}: tellers {tellers} (writers {}), \
                     clean {clean}, spanning X {spanning}, writers X-over {} REF-over {}; {}",
                    row.pos + 1,
                    String::from_utf8_lossy(&row.ref_al),
                    String::from_utf8_lossy(&row.alt_al),
                    self.rows[0].pos + 1,
                    String::from_utf8_lossy(&self.rows[0].ref_al),
                    String::from_utf8_lossy(&self.rows[0].alt_al),
                    by.len(),
                    reads.len(),
                    w.iter().filter(|t| t.tells()).count(),
                    w.iter().filter(|t| t.x - t.f >= DELTA).count(),
                    w.iter().filter(|t| t.r - t.f >= DELTA).count(),
                    verdicts.join(" "),
                );
            }
            log::debug!(
                "haplotype judge: the found allele {}:{}>{} (written by {} of {} reads) → {}",
                row.pos + 1,
                String::from_utf8_lossy(&row.ref_al),
                String::from_utf8_lossy(&row.alt_al),
                by.len(),
                reads.len(),
                if keep { "another allele, a state" } else { "X's or REF's spelling, not a state" },
            );
            if !keep {
                dropped.push(o);
            }
        }
        if dropped.is_empty() {
            return;
        }
        self.haps.retain(|h| !h.rows.iter().any(|r| dropped.contains(r)));
        self.max_shift = (0..self.rows.len())
            .filter(|r| !dropped.contains(r))
            .map(|r| self.rows[r].ref_al.len().abs_diff(self.rows[r].alt_al.len()) as i64)
            .sum();
    }

    /// The haplotypes: the empty context, then the `ordered` contexts while the
    /// haplotypes stay within `cap`, each crossed with every state (REF, X, each other
    /// state) but for a state at one of the context's rows' place. Returns the row sets,
    /// the contexts kept, and whether one was left out.
    fn close_over_site(&self, ordered: Vec<Vec<usize>>, cap: usize) -> (Vec<Vec<usize>>, usize, bool) {
        let mut states: Vec<Vec<usize>> = vec![vec![], vec![0]];
        states.extend((1..=self.n_states).map(|r| vec![r]));
        let cross = |ctx: &Vec<usize>| -> Vec<Vec<usize>> {
            states
                .iter()
                .filter(|st| st.iter().all(|&s| ctx.iter().all(|&c| !self.conflict(s, c))))
                .map(|st| {
                    let mut h: Vec<usize> = ctx.iter().chain(st).copied().collect();
                    h.sort_unstable();
                    h
                })
                .collect()
        };
        let mut seen: BTreeSet<Vec<usize>> = BTreeSet::from([vec![]]);
        let mut haps: Vec<Vec<usize>> = cross(&vec![]);
        let (mut n_contexts, mut left_out) = (1usize, false);
        for c in ordered {
            if !seen.insert(c.clone()) {
                continue;
            }
            let add = cross(&c);
            if haps.len() + add.len() > cap {
                left_out = true;
                continue;
            }
            haps.extend(add);
            n_contexts += 1;
        }
        (haps, n_contexts, left_out)
    }
}

/// One read's best scores at a found allele ([`Site::telling`]): on the haplotypes holding
/// it (`h` the first reaching it), on those carrying X and on those holding REF's bases;
/// once read, whether the bases it shows of haplotype `h` lie whole in one carrying X, and
/// in one holding REF's ([`Site::lies_in_x_or_ref`]).
#[derive(Clone, Copy, Debug)]
struct Telling {
    f: i32,
    x: i32,
    r: i32,
    h: usize,
    lies: Option<(bool, bool)>,
}

impl Telling {
    /// The allele's best haplotype beats every one carrying X and every one holding REF's
    /// bases by [`DELTA`].
    fn tells(&self) -> bool {
        self.f - self.x >= DELTA && self.f - self.r >= DELTA
    }
}

/// Which reads a [`Tally`] counts.
#[derive(Clone, Copy, Debug)]
enum Pick {
    /// Those that tell the allele ([`Telling::tells`]).
    Tells,
    /// Those that tell it and fit it whole (no confident base of theirs edited on its
    /// haplotype), the bases they show of it lying in no haplotype carrying X or holding
    /// REF's bases wherever placed: reads whose own bases are that allele's, X and REF each
    /// missing them by a confident base.
    Clean,
    /// Those that tell it from X either way (span where they differ).
    SpansX,
}

/// The reads' tellings at one found allele, each read on first use.
struct Tally<'a> {
    site: &'a Site,
    reads: &'a [ReadView],
    holds: Vec<bool>,
    /// The reads writing it (indices into `reads`, sorted).
    by: &'a [usize],
    /// The writers, then the rest.
    order: Vec<usize>,
    cache: Vec<Option<Telling>>,
}

impl<'a> Tally<'a> {
    fn new(site: &'a Site, reads: &'a [ReadView], holds: Vec<bool>, by: &'a [usize]) -> Self {
        let order = by.iter().copied().chain((0..reads.len()).filter(|i| by.binary_search(i).is_err())).collect();
        Tally { site, reads, holds, by, order, cache: vec![None; reads.len()] }
    }

    fn get(&mut self, i: usize) -> Telling {
        if let Some(t) = self.cache[i] {
            return t;
        }
        let t = self.site.telling(&self.reads[i], &self.holds);
        self.cache[i] = Some(t);
        t
    }

    /// Read `i`'s bases on the allele's best haplotype: lying in one carrying X, in one
    /// holding REF's bases.
    fn lies(&mut self, i: usize) -> (bool, bool) {
        let t = self.get(i);
        if let Some(l) = t.lies {
            return l;
        }
        let l = self.site.lies_in_x_or_ref(&self.reads[i], t.h);
        self.cache[i] = Some(Telling { lies: Some(l), ..t });
        l
    }

    fn picks(&mut self, i: usize, pick: Pick) -> bool {
        let t = self.get(i);
        match pick {
            Pick::Tells => t.tells(),
            Pick::Clean => t.tells() && t.f >= 0 && self.lies(i) == (false, false),
            Pick::SpansX => (t.f - t.x).abs() >= DELTA,
        }
    }

    /// The reads `pick` takes, writers first, counted up to `enough`.
    fn count(&mut self, enough: f64, pick: Pick) -> usize {
        let mut n = 0usize;
        for k in 0..self.order.len() {
            if self.picks(self.order[k], pick) {
                n += 1;
                if n as f64 >= enough {
                    break;
                }
            }
        }
        n
    }

    /// Whether the writers tell it from X and from REF (`exact`: their bases there lying
    /// in no haplotype carrying X, or holding REF's): some writer from each (`k` 0), or
    /// `k` writers from both.
    fn spell(&mut self, k: usize, exact: bool) -> bool {
        let (mut tx, mut tr, mut both) = (0usize, 0usize, 0usize);
        for j in 0..self.by.len() {
            let i = self.by[j];
            let t = self.get(i);
            let (mut x, mut r) = (t.f - t.x >= DELTA, t.f - t.r >= DELTA);
            if exact && (x || r) {
                let (in_x, in_ref) = self.lies(i);
                (x, r) = (x && !in_x, r && !in_ref);
            }
            tx += usize::from(x);
            tr += usize::from(r);
            both += usize::from(x && r);
        }
        if k == 0 {
            tx > 0 && tr > 0
        } else {
            both >= k
        }
    }

    /// Whether `k` writers' fragments (a writer and its mate among the site's reads) tell
    /// it as `exact` does, neither mate told from it by X or REF.
    fn fragments(&mut self, k: usize) -> bool {
        let mut by_name: BTreeMap<Vec<u8>, Vec<usize>> = BTreeMap::new();
        for &i in self.by {
            by_name.entry(self.reads[i].name.clone()).or_default().push(i);
        }
        for (i, read) in self.reads.iter().enumerate() {
            if let Some(v) = by_name.get_mut(&read.name).filter(|_| !read.name.is_empty()) {
                if !v.contains(&i) {
                    v.push(i);
                }
            }
        }
        let mut n = 0usize;
        for members in by_name.into_values() {
            let t: Vec<Telling> = members.iter().map(|&i| self.get(i)).collect();
            let refuted = t.iter().any(|t| t.x - t.f >= DELTA || t.r - t.f >= DELTA);
            if !refuted && members.iter().any(|&i| self.get(i).tells() && self.lies(i) == (false, false)) {
                n += 1;
            }
        }
        n >= k.max(1)
    }

    fn decide(&mut self, rule: FoundRule) -> bool {
        let floor = (GUARD_REF_SHARE * self.reads.len() as f64).max(f64::from(GUARD_MIN_READS));
        match rule {
            FoundRule::All => true,
            FoundRule::Told(0) => self.count(floor, Pick::Tells) as f64 >= floor,
            FoundRule::Told(k) => self.count(k as f64, Pick::Tells) >= k,
            FoundRule::Spell(k) => self.spell(k, false),
            FoundRule::Exact(k) => self.spell(k, true),
            FoundRule::Frag(k) => self.fragments(k),
            FoundRule::Clean(k) => self.count(floor, Pick::Tells) as f64 >= floor || self.count(k as f64, Pick::Clean) >= k,
            FoundRule::Share => {
                let n = self.count(f64::INFINITY, Pick::Tells);
                let spanning = self.count(f64::INFINITY, Pick::SpansX);
                n >= 1 && n as f64 >= GUARD_REF_SHARE * spanning as f64
            }
        }
    }
}

/// Whether `needle` lies whole in `hay`, a base [`WILD`] on either side matching any.
fn holds_bases(hay: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() {
        return true;
    }
    hay.windows(needle.len()).any(|w| w.iter().zip(needle).all(|(&h, &n)| h == n || h == WILD || n == WILD))
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
                StateGraph { class, row, graph: self.graph_of(row, &optional) }
            })
            .collect()
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
            let id = gb.segment(Some(r), (s, e), self.rows[r].alt_al.clone());
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

    /// An other state whose allele spells X's haplotype (or REF's) with rows of the reads'
    /// reach inside its place (an allele the reads spell over X with a listed neighbour
    /// beside it) is that haplotype: X's or REF's, never another allele's. Its graph is
    /// dropped, as the enumeration drops each haplotype of it spelling one of theirs.
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

    /// [`Site::drop_spellings`] on the graphs: a found allele stays a state only where the
    /// reads tell it from X and from REF by [`DELTA`], the state's graph beating REF's and
    /// X's, on the found floor of reads (`told`: [`GUARD_MIN_READS`] and
    /// [`GUARD_REF_SHARE`] of the site's; `toldK`: K reads), those that write it first,
    /// then the rest. The other rules read the enumerated haplotypes' bases and are
    /// weighed under the enumerated scorer only; here they fall back to `told`, with a
    /// warning.
    fn drop_graph_spellings(&mut self, found: &[(Row, &[usize])], reads: &[ReadView], rule: FoundRule) {
        let floor = match rule {
            FoundRule::All => return,
            FoundRule::Told(k) if k > 0 => k as f64,
            FoundRule::Told(_) => (GUARD_REF_SHARE * reads.len() as f64).max(f64::from(GUARD_MIN_READS)),
            other => {
                static WARNED: std::sync::Once = std::sync::Once::new();
                WARNED.call_once(|| {
                    log::warn!(
                        "sequence judge: the found-allele rule {other:?} is weighed under the enumerated scorer only \
                         (GBCMS_JUDGE_SCORER=enumerated); the graph scorer applies told"
                    )
                });
                (GUARD_REF_SHARE * reads.len() as f64).max(f64::from(GUARD_MIN_READS))
            }
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

    /// X's site in window coordinates, `[s, e)`.
    pub(crate) fn site_window(&self) -> (i64, i64) {
        self.site_window
    }

    /// The reference window: its genome start and bases.
    pub(crate) fn window(&self) -> (i64, &[u8]) {
        (self.start, &self.reference)
    }

    /// X's REF span in window coordinates.
    pub(crate) fn x_span(&self) -> (i64, i64) {
        self.spans[0]
    }
}

/// One read's evidence on the context rows (sorted row indices).
#[derive(Default)]
struct Evidence {
    shown: Vec<usize>,
    refuted: Vec<usize>,
}

/// The half-width (past its length changes) of the wider band tried before the full
/// programme.
const WIDE_SLACK: i64 = 4 * BAND_SLACK;

/// The rows within a read's band reach in `hap` ([`Scorer::near`]).
fn near_rows(scorer: &Scorer, site: &Site, hap: &Hap) -> Vec<usize> {
    scorer.near(site, hap)
}

/// The k-mer length of [`kmer_bound`].
const BOUND_K: usize = 8;

/// An upper bound on the score of `read` (masked bases [`WILD`]) anywhere in `hap`: the
/// read's k-mers holding no masked base that occur nowhere in `hap` are each broken by
/// some edit of any alignment; one edit breaks at most `BOUND_K` of them and costs at
/// least a mismatch (a gap or a base left over more per k-mer). 0 when `hap` holds a
/// masked base (any base fits there).
fn kmer_bound(read: &[u8], hap: &[u8]) -> i32 {
    let code = |b: u8| match b {
        b'A' => Some(0u32),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    };
    if read.len() < BOUND_K || hap.iter().any(|&b| code(b).is_none()) {
        return 0;
    }
    let mask = (1u32 << (2 * BOUND_K)) - 1;
    let mut present = vec![0u64; 1 << (2 * BOUND_K - 6)];
    let mut v = 0u32;
    for (i, &b) in hap.iter().enumerate() {
        v = ((v << 2) | code(b).unwrap_or(0)) & mask;
        if i + 1 >= BOUND_K {
            present[(v >> 6) as usize] |= 1 << (v & 63);
        }
    }
    let (mut missing, mut v, mut run) = (0usize, 0u32, 0usize);
    for &b in read {
        match code(b) {
            Some(c) => {
                v = ((v << 2) | c) & mask;
                run += 1;
            }
            None => run = 0,
        }
        if run >= BOUND_K && present[(v >> 6) as usize] & (1 << (v & 63)) == 0 {
            missing += 1;
        }
    }
    score::MISMATCH * missing.div_ceil(BOUND_K) as i32
}

fn is_subset(a: &[usize], b: &[usize]) -> bool {
    a.iter().all(|x| b.binary_search(x).is_ok())
}

/// Scores one read against the site's haplotypes, each trimmed to the read's reach and
/// banded around its own alignment mapped into the haplotype. Haplotypes applying the
/// same rows within the read's trimmed reach read the same there and are scored once;
/// those applying the same rows within its band share the banded programme.
pub(crate) struct Scorer<'r> {
    read: &'r ReadView,
    /// The read's reference positions, first and last (soft clips extrapolated).
    span: (i64, i64),
    finals: HashMap<Vec<usize>, i32>,
    banded: HashMap<Vec<usize>, score::Run>,
}

impl<'r> Scorer<'r> {
    pub(crate) fn new(read: &'r ReadView) -> Self {
        let first = read.refpos.iter().copied().min().unwrap_or(read.ref_start);
        let last = read.refpos.iter().copied().max().unwrap_or(read.ref_end);
        Scorer { read, span: (first.min(read.ref_start), last.max(read.ref_end)), finals: HashMap::new(), banded: HashMap::new() }
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

    pub(crate) fn score(&mut self, site: &Site, hap: &Hap) -> i32 {
        match self.score_or_bound(site, hap, None) {
            Ok(s) => s,
            Err(_) => unreachable!("no bound asked"),
        }
    }

    /// The read's score on `hap`; or, where the banded programme cannot vouch for it
    /// and the full one would be needed, `Err(bound)` when `beaten` is given and an upper
    /// bound on the score is already at most `beaten` (the score cannot matter there).
    /// The bound counts the read's confident k-mers found nowhere in the haplotype's
    /// trim: each edit of an alignment breaks at most [`BOUND_K`] of them and costs at
    /// least a mismatch, so the score is at most a mismatch per `BOUND_K` missing.
    pub(crate) fn score_or_bound(&mut self, site: &Site, hap: &Hap, beaten: Option<i32>) -> Result<i32, i32> {
        let wide = self.rows_within(site, hap, self.read.len() as i64 + site.max_shift);
        if let Some(&s) = self.finals.get(&wide) {
            return Ok(s);
        }
        let near = self.near(site, hap);
        let (seq, _, band) = self.prepare(site, hap, &near);
        let read = &self.read.seq;
        // A band the haplotype's ends cut is not the same relative to another trim.
        let banded = if band.cut {
            score::banded(read, seq, &band)
        } else {
            *self.banded.entry(near).or_insert_with(|| score::banded(read, seq, &band))
        };
        let s = if score::stands(&banded) && (banded.score == 0 || !score::copy_off_band(read, seq, &band)) {
            banded.score
        } else {
            if let Some(b) = beaten {
                let bound = kmer_bound(read, seq);
                if bound <= b {
                    return Err(bound);
                }
            }
            // A band WIDE_SLACK wide decides as the narrow one does, by its own edge and
            // floor (an indel long enough to leave it costs more than its score); else the
            // full programme.
            let (wseq, _, wide) = self.prepare_with(site, hap, &near_rows(self, site, hap), WIDE_SLACK);
            let w = score::banded(read, wseq, &wide);
            if score::stands_in(&w, WIDE_SLACK) && (w.score == 0 || !score::copy_off_band(read, wseq, &wide)) {
                w.score
            } else {
                score::full_score(read, seq)
            }
        };
        self.finals.insert(wide, s);
        Ok(s)
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

    /// The site of `target` among `listed` rows on `contig`, with reads of `l` bases, scored
    /// as `GBCMS_JUDGE_SCORER` says (the graph by default).
    pub(crate) fn site(contig: &[u8], target: &Row, listed: &[Row], l: usize) -> Site {
        site_with(scorer_kind(), contig, target, listed, l)
    }

    /// [`site`] with the scorer named.
    pub(crate) fn site_with(scorer: ScorerKind, contig: &[u8], target: &Row, listed: &[Row], l: usize) -> Site {
        let others: Vec<Row> = listed.iter().filter(|r| *r != target).cloned().collect();
        let read = ReadView::new(&vec![b'A'; l], &vec![30; l], 0, &[Cigar::Match(l as u32)], (0, l), 20).unwrap();
        let candidates = Candidates { listed: others, ..Candidates::default() };
        Site::build_with(scorer, target.clone(), &candidates, (0, contig.to_vec()), &[read]).unwrap()
    }

    /// [`site`] under the enumerated scorer: the tests of its contexts.
    pub(crate) fn enumerated_site(contig: &[u8], target: &Row, listed: &[Row], l: usize) -> Site {
        site_with(ScorerKind::Enumerated, contig, target, listed, l)
    }

    #[test]
    fn up_to_four_rows_every_subset_meets_every_state() {
        let c = random_seq(600, 9);
        let base = |p: i64| (c[p as usize - 1] as char).to_string();
        let other = |p: i64| if c[p as usize - 1] == b'A' { "C".to_string() } else { "A".to_string() };
        let x = row(300, &base(300), &other(300));
        let near: Vec<Row> = [250, 270, 330, 360].iter().map(|&p| row(p, &base(p), &other(p))).collect();
        let s = enumerated_site(&c, &x, &near, 100);
        // An SNV's site: REF, X and the two other bases; 16 contexts of four rows.
        assert_eq!(s.n_states, 2);
        assert_eq!(s.n_contexts, 16);
        assert_eq!(s.haps.len(), 16 * 4);
        let alt = s.haps.iter().filter(|h| h.class == Class::Alt).count();
        let reference = s.haps.iter().filter(|h| h.class == Class::Ref).count();
        assert_eq!((alt, reference), (16, 16), "every context with X and with REF");
        // A row out of the reads' reach is no context.
        let far = row(450, &base(450), &other(450));
        let s = enumerated_site(&c, &x, &[near[0].clone(), far], 100);
        assert_eq!(s.n_contexts, 2);
    }

    #[test]
    fn past_four_rows_the_contexts_are_each_row_and_the_pairs_reads_show_together() {
        let c = random_seq(600, 13);
        let base = |p: i64| (c[p as usize - 1] as char).to_string();
        let other = |p: i64| if c[p as usize - 1] == b'A' { "C".to_string() } else { "A".to_string() };
        let x = row(300, &base(300), &other(300));
        let near: Vec<Row> = [260, 270, 280, 320, 330, 340].iter().map(|&p| row(p, &base(p), &other(p))).collect();
        // Three reads show the rows at 260 and 270 together, two the one at 330 alone.
        let mut both = c.clone();
        for r in &near[..2] {
            both[r.pos as usize] = r.alt_al[0];
        }
        let mut lone = c.clone();
        lone[near[4].pos as usize] = near[4].alt_al[0];
        let view = |h: &[u8], s: usize| ReadView::new(&h[s..s + 100], &[30; 100], s as i64, &[Cigar::Match(100)], (0, 100), 20).unwrap();
        let reads: Vec<ReadView> =
            (0..3).map(|i| view(&both, 230 + i)).chain((0..2).map(|i| view(&lone, 300 + i))).chain([view(&c, 250)]).collect();
        let candidates = Candidates { listed: near, ..Candidates::default() };
        let s = Site::build_with(ScorerKind::Enumerated, x, &candidates, (0, c.clone()), &reads).unwrap();
        // {}, six singletons and the pair: eight contexts, each with REF, X and two bases.
        assert_eq!((s.n_contexts, s.haps.len(), s.cap_exceeded), (8, 32, false));
    }

    #[test]
    fn the_scorers_caches_never_change_a_score() {
        // Rows of every kind near X; reads over them carrying some. Each haplotype scored
        // through the caches, in order, equals it scored afresh.
        let c = random_seq(700, 23);
        let other = |p: usize| if c[p] == b'A' { b"C".to_vec() } else { b"A".to_vec() };
        let x = Row::new(350, &c[350..351], &other(350));
        let rows = vec![
            Row::new(300, &c[300..301], &other(300)),
            Row::new(320, &c[320..326], &c[320..321]),
            Row::new(362, &c[362..363], &[&c[362..363], b"GATTACA"].concat()),
            Row::new(380, &c[380..382], &[other(380)[0], other(381)[0]]),
            Row::new(410, &c[410..411], &other(410)),
            Row::new(430, &c[430..431], &other(430)),
        ];
        let carrier = |with: &[&Row], s: usize| {
            let mut all: Vec<&Row> = with.to_vec();
            all.sort_by_key(|r| std::cmp::Reverse(r.pos));
            let mut h = c.clone();
            for r in all {
                h.splice(r.pos as usize..r.pos as usize + r.ref_al.len(), r.alt_al.iter().copied());
            }
            ReadView::new(&h[s..s + 120], &[30; 120], s as i64, &[Cigar::Match(120)], (0, 120), 20).unwrap()
        };
        let reads = vec![
            carrier(&[&x, &rows[1]], 260),
            carrier(&[&x, &rows[2], &rows[3]], 290),
            carrier(&[&rows[0], &rows[4]], 290),
            carrier(&[], 330),
            carrier(&[&x, &rows[5]], 340),
        ];
        let candidates = Candidates { listed: rows, ..Candidates::default() };
        let s = Site::build_with(ScorerKind::Enumerated, x, &candidates, (0, c.clone()), &reads).unwrap();
        assert!(s.haps.len() > 8);
        for read in &reads {
            let mut cached = Scorer::new(read);
            for hap in &s.haps {
                let fresh = Scorer::new(read);
                let (seq, _, band) = fresh.prepare(&s, hap, &fresh.near(&s, hap));
                assert_eq!(cached.score(&s, hap), score::score(&read.seq, seq, Some(&band)));
            }
        }
    }

    #[test]
    fn a_state_at_a_context_rows_place_is_skipped() {
        // X: an SNV at 300; a listed deletion over 300 is a state; an SNV at 310 inside
        // a listed 20-base deletion is a context its deletion cannot share a haplotype with.
        let c = random_seq(600, 11);
        let base = |p: i64, n: usize| String::from_utf8(c[p as usize - 1..p as usize - 1 + n].to_vec()).unwrap();
        let alt_of = |p: i64| if c[p as usize - 1] == b'A' { "C".to_string() } else { "A".to_string() };
        let x = row(300, &base(300, 1), &alt_of(300));
        let del_at_x = row(298, &base(298, 5), &base(298, 1));
        let ctx_snv = row(340, &base(340, 1), &alt_of(340));
        let ctx_del = row(335, &base(335, 10), &base(335, 1));
        let s = enumerated_site(&c, &x, &[del_at_x.clone(), ctx_snv, ctx_del], 100);
        assert_eq!(s.n_states, 3, "the two other bases and the deletion over X");
        // Contexts: {}, {snv}, {del} (the two context rows meet): 3 contexts × 5 states.
        assert_eq!(s.n_contexts, 3);
        assert_eq!(s.haps.len(), 15);
        // Two insertions left-aligned a base apart are one site.
        let c = contig(23, 299, b"AC");
        let ins1 = row(300, "A", "AGT");
        let ins2 = row(300, "A", "ATTA");
        let s = enumerated_site(&c, &ins1, &[ins2], 100);
        assert_eq!((s.n_states, s.n_contexts, s.haps.len()), (1, 1, 3));
    }

    #[test]
    fn the_other_allele_tokens_parse() {
        let o = parse_other("found,fit,told,edges");
        assert!(o.found && o.fit && o.edges && !o.runs);
        assert_eq!((o.rule, o.found_min), (FoundRule::Told(0), FOUND_MIN_READS));
        let o = parse_other("found3,fit,clean,edges");
        assert_eq!((o.found, o.rule, o.found_min), (true, FoundRule::Clean(1), 3));
        let o = parse_other("none");
        assert_eq!((o.found, o.fit, o.rule), (false, false, FoundRule::All));
    }

    #[test]
    fn each_context_source_adds_only_its_own_rows() {
        use ContextSources::*;
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
        };
        let read = ReadView::new(&[b'A'; 150], &[30; 150], 230, &[Cigar::Match(150)], (0, 150), 20).unwrap();
        let shape = |sources| {
            let s = Site::build_sourced(ScorerKind::Graph, x.clone(), &candidates, (0, c.clone()), std::slice::from_ref(&read), FoundRule::All, sources)
                .unwrap();
            (s.n_states, s.rows.len() - 1 - s.n_states, s.masked.len())
        };
        // The site's states (the two other bases and the deletion) under every level.
        assert_eq!(shape(Seq), (3, 0, 0));
        assert_eq!(shape(Listed), (3, 1, 0));
        assert_eq!(shape(NoFound), (3, 2, 1));
        assert_eq!(shape(All), (3, 3, 1));
        assert_eq!(parse_contexts(""), Some(All));
        assert_eq!(parse_contexts("nofound"), Some(NoFound));
        assert_eq!(parse_contexts("Seq"), None);
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

    #[test]
    fn a_needle_lies_in_a_haystack_any_base_matching_wild() {
        assert!(holds_bases(b"ACGTACGT", b"GTAC"));
        assert!(holds_bases(b"ACGTACGT", b"GNAC"));
        assert!(holds_bases(b"ACNTACGT", b"GTAC"), "a wild base in the haystack takes any");
        assert!(!holds_bases(b"ACGTACGT", b"GTTC"));
        assert!(holds_bases(b"ACGT", b""));
    }

    /// (AT)10 T23 between unique flanks; X makes it (AT)11 T28; a found allele (AT)12 T23
    /// (the reference with +ATAT), written by reads that fit none of them whole.
    fn sibling_case() -> (Vec<u8>, Row, Row) {
        let mut c = random_seq(700, 41);
        let run: Vec<u8> = [&b"AT".repeat(10)[..], &b"T".repeat(23)[..]].concat();
        c[300] = b'C';
        c[301..301 + run.len()].copy_from_slice(&run);
        c[301 + run.len()] = b'G';
        let anchor = 301 + 19; // the run's last T
        let x = Row::new(anchor, b"T", b"TATTTTTT");
        let found = Row::new(anchor, b"T", b"TATAT");
        (c, x, found)
    }

    #[test]
    fn a_found_allele_one_read_shows_whole_stays_a_state_only_under_clean() {
        let (c, x, found) = sibling_case();
        let apply = |r: &Row| [&c[..r.pos as usize], &r.alt_al[..], &c[r.end() as usize..]].concat();
        let view = |h: &[u8], s: usize, n: usize| {
            ReadView::new(&h[s..s + n], &vec![30; n], s as i64, &[Cigar::Match(n as u32)], (0, n), 20).unwrap()
        };
        // the sibling (AT)12 T28, read from the left flank to 19 of its T's: the found
        // allele's haplotype whole; X misses it by a base, REF by two
        let sib = apply(&Row::new(x.pos, b"T", b"TATATTTTTT"));
        let teller = view(&sib, 301 + 24 + 19 - 100, 100);
        let reads = vec![teller, view(&apply(&x), 240, 100), view(&c, 245, 100)];
        let candidates = Candidates { found: vec![(found.clone(), vec![])], ..Candidates::default() };
        let states = |rule: FoundRule| {
            let s = Site::build_ruled(ScorerKind::Enumerated, x.clone(), &candidates, (0, c.clone()), &reads, rule).unwrap();
            s.haps.iter().filter(|h| h.class == Class::Other).count()
        };
        assert!(states(FoundRule::All) > 0, "the found allele is a state before any rule");
        assert_eq!(states(FoundRule::Told(0)), 0, "one read is below the found floor");
        assert!(states(FoundRule::Clean(1)) > 0, "one read shows it whole: another allele");
        // the same read with a confident error elsewhere tells it, yet not whole: dropped
        let mut noisy = sib.clone();
        let at = 301 + 24 + 19 - 100 + 10;
        noisy[at] = if noisy[at] == b'A' { b'C' } else { b'A' };
        let reads = vec![view(&noisy, 301 + 24 + 19 - 100, 100), view(&apply(&x), 240, 100), view(&c, 245, 100)];
        let s = Site::build_ruled(ScorerKind::Enumerated, x.clone(), &candidates, (0, c.clone()), &reads, FoundRule::Clean(1)).unwrap();
        assert_eq!(s.haps.iter().filter(|h| h.class == Class::Other).count(), 0);
    }
}

#[cfg(test)]
mod graph_tests {
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

    fn class_best(s: &Site, read: &ReadView) -> [Option<i32>; 3] {
        let mut sc = Scorer::new(read);
        let mut e = [None; 3];
        for h in &s.haps {
            let v = sc.score(s, h);
            let b = &mut e[h.class as usize];
            if b.is_none_or(|x| v > x) {
                *b = Some(v);
            }
        }
        e
    }

    fn graph_best(s: &Site, read: &ReadView) -> [Option<i32>; 3] {
        let gr = s.graph_read(read);
        let mut ws = Workspace::default();
        let mut best = [None; 3];
        for st in &s.graphs {
            let v = graph::score(&ShortReads, &gr, &st.graph, &mut ws);
            assert_eq!(v, graph::full_score(&ShortReads, &read.seq, &st.graph), "the band's score is the full programme's");
            let b = &mut best[st.class as usize];
            if b.is_none_or(|x| v > x) {
                *b = Some(v);
            }
        }
        best
    }

    #[test]
    fn at_four_rows_or_fewer_the_graph_scores_and_calls_as_the_enumeration() {
        // Sites of an SNV, an insertion or a deletion with one to four rows in reach of
        // every kind (two length changes in one repeat among them, rows at one place the
        // enumeration never applies together), and reads of random haplotypes with errors
        // and masked bases, written with their gaps or without. Where the enumeration
        // builds every combination, the graph's best per class is the enumeration's, and
        // so is every call.
        let mut st = 47u64;
        let mut compared = 0usize;
        for case in 0..96u64 {
            let mut c = random_seq(700, 9000 + case);
            if case % 3 == 0 {
                // a CA repeat beside X: two length changes in it are at one place
                c[372..392].copy_from_slice(b"CACACACACACACACACACA");
            }
            let x = match case % 3 {
                0 => Row::new(350, &c[350..351], &[*b"ACGT".iter().find(|&&b| b != c[350]).unwrap()]),
                1 => Row::new(350, &c[350..351], &[&c[350..351], b"TGA"].concat()),
                _ => Row::new(350, &c[350..356], &c[350..351]),
            };
            let n_rows = 1 + roll(&mut st, 4) as usize;
            let mut rows: Vec<Row> = Vec::new();
            if case % 3 == 0 {
                rows.push(Row::new(371, &c[371..372], &[&c[371..372], b"CA"].concat()));
                rows.push(Row::new(385, &c[385..388], &c[385..386]));
            }
            while rows.len() < n_rows {
                let side = if roll(&mut st, 2) == 0 { 290 + roll(&mut st, 45) } else { 375 + roll(&mut st, 45) };
                let r = random_row(&c, side as usize, &mut st);
                if !rows.contains(&r) && (r.end() < 345 || r.pos > 360) {
                    rows.push(r);
                }
            }
            // reads: haplotypes of random sets of the rows (that can lie together), with X or REF
            let mut reads: Vec<ReadView> = Vec::new();
            for k in 0..14 {
                let mut on: Vec<&Row> = Vec::new();
                for r in &rows {
                    if roll(&mut st, 2) == 0 && on.iter().all(|o| o.end() <= r.pos || r.end() <= o.pos) {
                        on.push(r);
                    }
                }
                if k % 2 == 0 && on.iter().all(|o| o.end() <= x.pos || x.end() <= o.pos) {
                    on.push(&x);
                }
                let (mut h, o) = applied(&c, &on);
                let s = 210 + roll(&mut st, 140) as usize;
                for _ in 0..roll(&mut st, 3) {
                    let i = (s + roll(&mut st, 120) as usize).min(h.len() - 1);
                    h[i] = b"ACGT"[roll(&mut st, 4) as usize];
                }
                let masked: Vec<usize> = (0..roll(&mut st, 3)).map(|_| roll(&mut st, 120) as usize).collect();
                if let Some(r) = read_on(&h, &o, s, 100 + roll(&mut st, 52) as usize, k % 3 == 0, &masked) {
                    reads.push(r);
                }
            }
            let cand = Candidates { listed: rows.clone(), ..Candidates::default() };
            let e = Site::build_with(ScorerKind::Enumerated, x.clone(), &cand, (0, c.clone()), &reads).unwrap();
            let g = Site::build_with(ScorerKind::Graph, x.clone(), &cand, (0, c.clone()), &reads).unwrap();
            if e.cap_exceeded {
                continue;
            }
            for (i, read) in reads.iter().enumerate() {
                let want = class_best(&e, read);
                assert_eq!(graph_best(&g, read), want, "case {case} read {i}: best per class (ALT, REF, OTHER)");
                assert_eq!(judge(&g, read, 20, true).call, judge(&e, read, 20, true).call, "case {case} read {i}: the call");
                compared += 1;
            }
        }
        assert!(compared > 800, "{compared} reads compared");
    }
}

