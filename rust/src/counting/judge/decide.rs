//! A read's call from its scores against the site's haplotypes, and the call mapped onto
//! the engine's `ClassifyResult`.
//!
//! With H_alt the haplotypes carrying X at its site, H_ref those with the reference's
//! bases there and H_other those with another allele: ALT where the best of H_alt beats
//! the best of the rest by [`DELTA`] (one confident base), REF and OTHER likewise, else
//! depth only. No priors; a tie is never credited.

use log::trace;
use rust_htslib::bam::record::{Aux, Cigar, Record};


use super::graph::{self, ShortReads, Workspace};
use super::haplotypes::{neither_trace, other_option, AlnTags, Class, Hap, NeitherRule, OwnGapRule, ReadView, RefBaseRule, Scorer, ScorerKind, Site};
use super::score::{Alignment, Placed, DELTA, WILD};
use crate::counting::carrier::{readable_range, ReadRules};
use crate::counting::utils::{soft_clips, ClassifyPhase, ClassifyResult};
use crate::counting::window::{allele_kind, AlleleKind};
use crate::shared::bam_utils::{find_read_pos, fragment_query_span, ref_end};
use crate::types::Variant;

/// The phase a judged call reports: the one downstream treats as decided by the read's
/// bases (as the junction judge's calls are). `Structural` would exempt an ALT call from
/// the sibling contest, `Alignment` would demote a pure indel's there; the judge weighs
/// the site's alleles itself, so neither special case is wanted.
pub(crate) const PHASE: ClassifyPhase = ClassifyPhase::MaskedCompare;

/// A read's call at X.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Call {
    Alt,
    Ref,
    Other,
    Depth,
}

/// A read's call, with the quality of the bases that decided it (REF and ALT calls),
/// its best score per class (ALT, REF, OTHER), and, where asked, whether its alignment
/// to the haplotype it was called on holds X's anchor there.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Verdict {
    pub(crate) call: Call,
    pub(crate) qual: u8,
    pub(crate) best: [Option<i32>; 3],
    pub(crate) holds_anchor: bool,
    /// The read bases (`ReadView` indices, `[first, after)`) its alignment to the
    /// haplotype it was called on places, where asked for the anchor.
    pub(crate) placed: Option<(usize, usize)>,
    /// An ALT call on a substitution that read every base X changes confidently, as X's.
    pub(crate) confirmed: bool,
    /// For an ALT call, the states holding X's change the read spells whole
    /// ([`Site::spelled_holders`]).
    pub(crate) spells: u32,
    /// Measurement only, the "fits nothing" decision package: the site states the full
    /// programme decided (the band could not), whether the rule sent the read to depth,
    /// and its fit to the local reference neighbourhood where that was computed.
    pub(crate) fits: FitsTrace,
    /// The read's best score over every candidate sequence is below [`FITS_NOTHING`]: it
    /// fits none of them (an upper bound under the enumerated scorer), whether or not the
    /// fit rule made it depth only. Diagnostic (UNEXPLAINED_READS).
    pub(crate) fits_none: bool,
    /// Measurement only, TEMPORARY: where a substitution's telling bases lie in the read's
    /// soft clip (asked for the anchor), the clip's bases around them ([`ClipBases`]).
    pub(crate) clip: Option<ClipBases>,
    /// Measurement only: what the read shows at X's telling bases that is neither allele's
    /// ([`NeitherTrace`]), where traced or the switch is on.
    pub(crate) neither: NeitherTrace,
    /// Measurement only: the states the general own-gap rule dropped ([`OwnGapTrace`]).
    pub(crate) owngap: OwnGapTrace,
}

/// Measurement only (`GBCMS_JUDGE_OWNGAP=all`): the states (ALT 0, REF 1, OTHER 2) a read
/// reached only through an indel of its own at their change, dropped before the margin
/// decided; the call the margin gave before; and whether more than one state stood within
/// the margin of the best (a tie) or one beat the rest (a win).
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct OwnGapTrace {
    pub(crate) dropped: [bool; 3],
    pub(crate) before: Option<Call>,
    pub(crate) tie: bool,
}

/// Measurement only, TEMPORARY: a REF or ALT call's read at X's telling bases (on the called
/// haplotype): how many there are and are read confidently, the neither bases among them (up
/// to four: offset among the telling bases, base, quality, how many other reads of the row
/// read it the same at the same place, and that place relative to X by the read's own
/// alignment, `i32::MIN` where it leaves it unplaced), and the read's confident bases and
/// mismatches on that haplotype elsewhere.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct NeitherTrace {
    pub(crate) on: bool,
    pub(crate) call: u8,
    pub(crate) size: u16,
    pub(crate) telling: u16,
    pub(crate) conf: u16,
    pub(crate) n: u8,
    pub(crate) bases: [(u16, u8, u8, u16, i32); 4],
    pub(crate) max_recur: u32,
    pub(crate) bg_conf: u32,
    pub(crate) bg_mm: u32,
    /// Of those mismatches, the ones no other read of the row reads the same at that place.
    pub(crate) bg_single: u32,
    pub(crate) moved: bool,
    /// The neither bases lying where both alleles' reads are read alike (`window`).
    pub(crate) window: u8,
}

/// Measurement only: what the "fits nothing" rule saw of one read (see [`Verdict::fits`]).
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct FitsTrace {
    pub(crate) full: u8,
    pub(crate) states: u8,
    pub(crate) fired: bool,
    pub(crate) neighbourhood: Option<i32>,
}

/// Measurement only, TEMPORARY (the clip-admission measurement): a substitution's
/// telling bases read in a soft clip, and the clip's other bases on the called
/// haplotype: `beyond` past them towards the read's end, `between` from the aligned
/// bases to them; each (confident and matching, confident and not, unconfident).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ClipBases {
    pub(crate) tail: bool,
    pub(crate) beyond: (u16, u16, u16),
    pub(crate) between: (u16, u16, u16),
    /// Telling bases inside the clip, of all of them.
    pub(crate) in_clip: (u16, u16),
}

/// Measurement only, TEMPORARY: which reads off the anchor a substitution row admits
/// (`GBCMS_JUDGE_CLIPSUB`): `admit` (RJ-13, the default); `off` (none, by a clip or
/// not, as the rules and the sign-out counter admit none); `anchored:N` (a clip only
/// where it holds N confident bases past X's telling bases, matching the called
/// haplotype); `clean:N` (as anchored, with no confident mismatch among the clip's
/// other bases).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ClipSub {
    Admit,
    Off,
    Anchored(u16),
    Clean(u16),
}

pub(crate) fn clip_sub() -> ClipSub {
    static OPT: std::sync::OnceLock<ClipSub> = std::sync::OnceLock::new();
    *OPT.get_or_init(|| {
        let v = std::env::var("GBCMS_JUDGE_CLIPSUB").unwrap_or_default();
        let n = |s: &str| s.parse::<u16>().unwrap_or(0);
        match v.split_once(':') {
            _ if v == "off" => ClipSub::Off,
            Some(("anchored", k)) => ClipSub::Anchored(n(k)),
            Some(("clean", k)) => ClipSub::Clean(n(k)),
            _ => ClipSub::Admit,
        }
    })
}

/// Measurement only: [`ClipBases`] for a read whose telling bases (`idx`, read
/// indices) lie in its soft clip, on the alignment `wa` to `hap` (its columns from `off`);
/// None where none does.
fn clip_bases(read: &ReadView, hap: &Hap, wa: &Alignment, off: usize, idx: &[usize]) -> Option<ClipBases> {
    let (lo, hi) = (*idx.iter().min()?, *idx.iter().max()?);
    let n = read.seq.len();
    let (lead, tail) = read.clipped_ends();
    let tally = |r: std::ops::Range<usize>| {
        let mut t = (0u16, 0u16, 0u16);
        for i in r.filter(|i| !idx.contains(i)) {
            let b = read.seq[i];
            match wa.placed[i] {
                _ if b == WILD => t.2 += 1,
                Placed::At(c) if hap.base(c + off) == WILD => t.2 += 1,
                Placed::At(c) if hap.base(c + off) == b => t.0 += 1,
                _ => t.1 += 1,
            }
        }
        t
    };
    let inside = |f: &dyn Fn(usize) -> bool| idx.iter().filter(|&&i| f(i)).count() as u16;
    if lo < lead {
        Some(ClipBases {
            tail: false,
            beyond: tally(0..lo),
            between: tally(hi + 1..lead.max(hi + 1)),
            in_clip: (inside(&|i| i < lead), idx.len() as u16),
        })
    } else if tail > 0 && hi >= n - tail {
        Some(ClipBases {
            tail: true,
            beyond: tally(hi + 1..n),
            between: tally((n - tail).min(lo)..lo),
            in_clip: (inside(&|i| i >= n - tail), idx.len() as u16),
        })
    } else {
        None
    }
}

/// The fit rule: which test, if any, makes a read that fits none of the site's candidate
/// sequences depth only on the graph scorer (decided 2026-10-25: `xs`). The other options
/// stay behind `GBCMS_JUDGE_FITS` for measurement, TEMPORARY. The enumerated scorer keeps
/// applying [`FITS_NOTHING`] to upper bounds on every haplotype's score. Whatever the rule,
/// such a read is named in the row's `UNEXPLAINED_READS` diagnostic.
/// - `xs` (the default): the aligner places the read elsewhere at least as well (XS at
///   least AS; MAPQ 0 where it wrote no XS) and its best score is below [`FITS_NOTHING`]:
///   a read the aligner doubts and no candidate sequence explains. A read that fits none
///   but that the aligner places here alone is judged by its margin, as any other.
/// - `exact`: the read's best score over every state below [`FITS_NOTHING`].
/// - `off`: no fit test; the margin alone decides.
/// - `scaled[:R]`: below R confident mismatches per 100 confident read bases ([`SCALED_RATE`]
///   where R is not given), so the budget follows the read's length.
/// - `relative[:M]`: below the read's fit to the local reference neighbourhood (anywhere
///   on the window's reference) by more than M (default one confident base): the read fits
///   somewhere nearby better than in any state here.
/// - `xspure`: the aligner places the read elsewhere at least as well, whatever its score;
///   decided from the tags before any programme runs.
/// - `xsscaled[:R]`: depth only where both `xs` and `scaled:R` would make it so: the aligner
///   doubts the read and its best score is below both [`FITS_NOTHING`] and R confident
///   mismatches per 100 confident read bases. Every read either rule judges is judged.
///
/// With `GBCMS_JUDGE_FITS_MEASURE=1`, every read whose best score is below a confident
/// base also gets its neighbourhood fit computed, for the trace.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum FitsRule {
    Exact,
    Off,
    Scaled(f64),
    Relative(i32),
    Xs,
    XsPure,
    XsScaled(f64),
}

/// The `scaled` option's rate where none is given: confident mismatches per 100 confident
/// read bases.
pub(crate) const SCALED_RATE: f64 = 10.0;

pub(crate) fn fits_rule() -> FitsRule {
    static RULE: std::sync::OnceLock<FitsRule> = std::sync::OnceLock::new();
    *RULE.get_or_init(|| parse_fits(&std::env::var("GBCMS_JUDGE_FITS").unwrap_or_default()))
}

/// The fit rule a `GBCMS_JUDGE_FITS` value names; `xs` where it names none.
fn parse_fits(v: &str) -> FitsRule {
    let (name, arg) = v.split_once(':').map_or((v, None), |(n, a)| (n, Some(a)));
    match name {
        "off" => FitsRule::Off,
        "scaled" => FitsRule::Scaled(arg.and_then(|a| a.parse().ok()).unwrap_or(SCALED_RATE)),
        "relative" => FitsRule::Relative(arg.and_then(|a| a.parse().ok()).unwrap_or(DELTA)),
        "exact" => FitsRule::Exact,
        "xspure" => FitsRule::XsPure,
        "xsscaled" => FitsRule::XsScaled(arg.and_then(|a| a.parse().ok()).unwrap_or(SCALED_RATE)),
        _ => FitsRule::Xs,
    }
}

fn fits_measure() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GBCMS_JUDGE_FITS_MEASURE").is_ok_and(|v| v == "1"))
}

/// A read whose best score on every candidate sequence is below this (ten confident
/// mismatches) fits none of them: depth only where the fit rule says so ([`FitsRule`]),
/// and named in the row's `UNEXPLAINED_READS` diagnostic.
pub(crate) const FITS_NOTHING: i32 = 10 * crate::counting::judge::score::MISMATCH;

/// Whether `rule` makes a read whose best score over the site's states is `top` depth
/// only; the neighbourhood fit, where the rule or the measurement computed it.
fn fits_nothing(rule: FitsRule, site: &Site, read: &ReadView, top: i32) -> (bool, Option<i32>) {
    let neighbourhood = || site.neighbourhood_fit(read);
    let measured = (fits_measure() && top <= -DELTA).then(neighbourhood);
    let fired = match rule {
        FitsRule::Exact => top < FITS_NOTHING,
        FitsRule::Off => false,
        FitsRule::Scaled(rate) => top < scaled_floor(rate, read),
        FitsRule::Relative(margin) => {
            // a neighbourhood fit is at most 0, so a read within the margin of 0 is kept
            if top >= -margin {
                false
            } else {
                let nb = measured.unwrap_or_else(neighbourhood);
                return (top < nb - margin, Some(nb));
            }
        }
        FitsRule::Xs => top < FITS_NOTHING && read.tags.fits_elsewhere(),
        FitsRule::XsPure => read.tags.fits_elsewhere(),
        FitsRule::XsScaled(rate) => top < FITS_NOTHING && top < scaled_floor(rate, read) && read.tags.fits_elsewhere(),
    };
    (fired, measured)
}

/// The `scaled` option's floor for `read`: `rate` confident mismatches per 100 of its
/// confident bases, rounded up.
fn scaled_floor(rate: f64, read: &ReadView) -> i32 {
    let budget = (rate * read.confident() as f64 / 100.0).ceil() as i32;
    budget * crate::counting::judge::score::MISMATCH
}

/// The class whose best beats every other class's best by [`DELTA`]; else, where
/// another allele and X (or REF) both beat the third class by as much, partial
/// evidence (the read is that allele or X, never credited: RJ-25); else depth only.
pub(crate) fn decide(best: [Option<i32>; 3]) -> Call {
    let v = best.map(|b| b.map_or(i64::MIN / 4, i64::from));
    let calls = [Call::Alt, Call::Ref, Call::Other];
    let delta = i64::from(DELTA);
    for k in 0..3 {
        let rest = (0..3).filter(|&o| o != k).map(|o| v[o]).max().unwrap_or(i64::MIN / 4);
        if best[k].is_some() && v[k] - rest >= delta {
            return calls[k];
        }
    }
    // ALT 0, REF 1, OTHER 2: another allele tied with one of REF and X, the other beaten
    if tied_with(best).is_some() {
        return Call::Other;
    }
    Call::Depth
}

/// Where no class beats both others by [`DELTA`]: the class (ALT 0, REF 1) another allele
/// ties with, the remaining one beaten by both.
fn tied_with(best: [Option<i32>; 3]) -> Option<usize> {
    let v = best.map(|b| b.map_or(i64::MIN / 4, i64::from));
    let delta = i64::from(DELTA);
    [(0, 1), (1, 0)]
        .into_iter()
        .find(|&(tied, beaten)| {
            best[2].is_some() && best[tied].is_some() && (v[2] - v[tied]).abs() < delta && v[2].min(v[tied]) - v[beaten] >= delta
        })
        .map(|(tied, _)| tied)
}

/// [`judge_kept`] with nothing kept from another read.
#[cfg(test)]
pub(crate) fn judge(site: &Site, read: &ReadView, min_baseq: u8, ask_anchor: bool) -> Verdict {
    judge_kept(site, read, min_baseq, ask_anchor, None).0
}

/// Judge one read against the site's haplotypes. A REF or ALT call stands only where
/// the read's alignment to the winning haplotype reads the bases that tell X from REF
/// at X's site: an ALT call (X a substitution or a delins) every one of X's own,
/// confidently and as X's; a REF call on a substitution every one of REF's, covered (as
/// the size-specific rules require of an MNP's block). A read that shows part of them is
/// depth only. With `ask_anchor` (a read the genome places
/// off X's anchor), the verdict says whether its alignment to the winning haplotype
/// holds the anchor there.
///
/// `kept`: what an earlier read alike in its bases (masked) and placement left (its
/// programmes' scores, bands and paths, and the alignments built since, all fixed by those
/// alone); under the graph scorer the read's [`Scored`] is given back for the next such
/// read. The call is the read's own: the aligner's tags, its qualities and whether its
/// anchor is asked for are weighed afresh.
pub(crate) fn judge_kept(site: &Site, read: &ReadView, min_baseq: u8, ask_anchor: bool, kept: Option<std::rc::Rc<Scored>>) -> (Verdict, Option<std::rc::Rc<Scored>>) {
    match site.scorer {
        ScorerKind::Enumerated => (judge_enumerated(site, read, min_baseq, ask_anchor), None),
        ScorerKind::Graph => judge_graph(site, read, min_baseq, ask_anchor, kept),
    }
}

/// A read's programmes over a site's state graphs, fixed by its bases (masked) and
/// placement: each state's best score (per class, the best and the first state reaching
/// it), the band that decided each, the rows each state's best path takes, and the
/// haplotypes those spell and the read's alignments to them, built as asked for and kept
/// for the reads alike.
pub(crate) struct Scored {
    best: [Option<(i32, usize)>; 3],
    halves: Vec<Option<i64>>,
    rows: Vec<Vec<usize>>,
    haps: std::cell::RefCell<Vec<(usize, std::rc::Rc<Hap>)>>,
    aligns: std::cell::RefCell<Vec<(usize, Aligned)>>,
}

/// [`judge`] against the enumerated haplotypes.
fn judge_enumerated(site: &Site, read: &ReadView, min_baseq: u8, ask_anchor: bool) -> Verdict {
    let mut scorer = Scorer::new(read);
    // Per class, the best score and the first haplotype reaching it. A haplotype whose
    // banded score cannot stand is scored in full only where an upper bound on its score
    // could still raise its class's best (or reach it first), so a read fitting none of
    // the haplotypes well (a dense cluster) does not run the full programme on each.
    let mut best: [Option<(i32, usize)>; 3] = [None; 3];
    let mut pending = Vec::new();
    let mut top = i32::MIN;
    for (h, hap) in site.haps.iter().enumerate() {
        match scorer.score_or_bound(site, hap, Some(i32::MAX)) {
            Ok(s) => {
                top = top.max(s);
                raise(&mut best, hap.class as usize, s, h)
            }
            Err(bound) => {
                top = top.max(bound);
                pending.push(h)
            }
        }
    }
    // A read that fits none of the haplotypes within FITS_NOTHING (an upper bound on its
    // best score below it: a paralog's read, or one from a cluster of changes the site's
    // haplotypes do not hold) is depth only: a margin between two misfits is no evidence.
    if top < FITS_NOTHING {
        trace!("haplotype judge: the read fits no haplotype (best at most {top}) → depth only");
        return Verdict { call: Call::Depth, qual: 0, best: [None; 3], holds_anchor: false, placed: None, confirmed: false, spells: 0, clip: None, fits: FitsTrace::default(), fits_none: true, neither: NeitherTrace::default(), owngap: OwnGapTrace::default() };
    }
    for h in pending {
        let hap = &site.haps[h];
        let class = hap.class as usize;
        // a bound at most the class's best cannot raise it, nor reach it first past it
        let beaten = best[class].map(|(b, bh)| if h < bh { b - 1 } else { b });
        let s = match scorer.score_or_bound(site, hap, beaten) {
            Ok(s) => s,
            Err(_) => continue,
        };
        raise(&mut best, class, s, h);
    }
    let aligns = std::cell::RefCell::new(Vec::new());
    let ex = Examined::new(site, &scorer, |h| std::rc::Rc::new(site.haps[h].clone()), &aligns);
    conclude(site, read, &ex, best, min_baseq, ask_anchor)
}

/// Per class, the best score and the first candidate (haplotype or state) reaching it.
fn raise(best: &mut [Option<(i32, usize)>; 3], class: usize, s: i32, h: usize) {
    let slot = &mut best[class];
    if slot.is_none_or(|(b, bh)| s > b || (s == b && h < bh)) {
        *slot = Some((s, h));
    }
}

/// [`judge`] on the site's graphs: per state, the read's best path over the window with
/// every context row in reach as a bubble. The haplotype a call is examined on (the fit
/// check, the bases that tell X from REF, the deciding quality, the anchor) is the one the
/// state's best path spells.
fn judge_graph(site: &Site, read: &ReadView, min_baseq: u8, ask_anchor: bool, kept: Option<std::rc::Rc<Scored>>) -> (Verdict, Option<std::rc::Rc<Scored>>) {
    let rule = fits_rule();
    // The aligner's own tags decide before any programme runs.
    if rule == FitsRule::XsPure && read.tags.fits_elsewhere() {
        let fits = FitsTrace { states: site.graphs.len() as u8, fired: true, ..FitsTrace::default() };
        let verdict = Verdict { call: Call::Depth, qual: 0, best: [None; 3], holds_anchor: false, placed: None, confirmed: false, spells: 0, clip: None, fits, fits_none: false, neither: NeitherTrace::default(), owngap: OwnGapTrace::default() };
        return (verdict, kept);
    }
    let gr = site.graph_read(read);
    let scored = kept.unwrap_or_else(|| std::rc::Rc::new(graph::with_workspaces(site.graphs.len(), |wss| score_states(site, &gr, wss))));
    let (best, halves) = (scored.best, &scored.halves);
    let top = best.iter().flatten().map(|b| b.0).max().unwrap_or(i32::MIN);
    let (fired, neighbourhood) = fits_nothing(rule, site, read, top);
    let fits = FitsTrace { full: halves.iter().filter(|h| h.is_none()).count() as u8, states: halves.len() as u8, fired, neighbourhood };
    let verdict = if fired {
        trace!("haplotype judge: the read fits no state's graph (best {top}, {rule:?}) → depth only");
        Verdict { call: Call::Depth, qual: 0, best: best.map(|b| b.map(|(s, _)| s)), holds_anchor: false, placed: None, confirmed: false, spells: 0, clip: None, fits, fits_none: false, neither: NeitherTrace::default(), owngap: OwnGapTrace::default() }
    } else {
        let scorer = Scorer::new(read);
        let ex = Examined::new(
            site,
            &scorer,
            |k| {
                if let Some((_, h)) = scored.haps.borrow().iter().find(|(j, _)| *j == k) {
                    return h.clone();
                }
                let h = std::rc::Rc::new(site.hap_on_path(k, &scored.rows[k]));
                scored.haps.borrow_mut().push((k, h.clone()));
                h
            },
            &scored.aligns,
        );
        conclude(site, read, &ex, best, min_baseq, ask_anchor)
    };
    (Verdict { fits, fits_none: top < FITS_NOTHING, ..verdict }, Some(scored))
}

/// The read's programmes over the site's states, a workspace each (`wss`): each state's
/// score, the band that decided it, and its best path traced back from the values its
/// scoring left there ([`Scored`]).
fn score_states(site: &Site, gr: &graph::GraphRead, wss: &mut [Workspace]) -> Scored {
    let mut best: [Option<(i32, usize)>; 3] = [None; 3];
    let graphs: Vec<&graph::Graph> = site.graphs.iter().map(|st| &st.graph).collect();
    let shares: Vec<Option<(usize, usize)>> = site.graphs.iter().map(|st| st.share).collect();
    let decided = graph::decided_states(&ShortReads, gr, &graphs, &shares, wss);
    let mut halves = Vec::with_capacity(site.graphs.len());
    for (k, (st, &(s, half))) in site.graphs.iter().zip(&decided).enumerate() {
        halves.push(half);
        raise(&mut best, st.class as usize, s, k);
    }
    let rows = (0..site.graphs.len())
        .map(|k| {
            let p = graph::kept_path(&ShortReads, gr, &site.graphs[k].graph, halves[k], &wss[k]);
            if best.iter().flatten().any(|&(s, h)| h == k && s != p.score) {
                trace!("haplotype judge: state {k}'s traced path scores {} against its score", p.score);
            }
            p.rows
        })
        .collect();
    Scored { best, halves, rows, haps: Default::default(), aligns: Default::default() }
}

/// The candidates a read's call is examined on, each built where first asked for and kept
/// for the call: its haplotype (from `build`: under the graph scorer the haplotype the
/// state's best path spells, from a traced programme) and the read's alignment to it. The
/// margin's checks and the call ask for one candidate's more than once (the fit check and
/// the winning haplotype are mostly one); both are fixed by the candidate, so each is
/// computed once.
pub(crate) struct Examined<'h, 'a, F: Fn(usize) -> std::rc::Rc<Hap>> {
    site: &'h Site,
    scorer: &'a Scorer<'a>,
    build: F,
    haps: std::cell::RefCell<Vec<(usize, std::rc::Rc<Hap>)>>,
    /// The read's alignments to the candidates: the call's own, or those kept for the reads
    /// alike in bases and placement ([`Scored`]), which the alignment is fixed by.
    aligns: &'a std::cell::RefCell<Vec<(usize, Aligned)>>,
}

/// A read's alignment to a candidate's haplotype, and where the haplotype's trim starts.
type Aligned = std::rc::Rc<(Alignment, usize)>;

impl<'h, 'a, F: Fn(usize) -> std::rc::Rc<Hap>> Examined<'h, 'a, F> {
    pub(crate) fn new(site: &'h Site, scorer: &'a Scorer<'a>, build: F, aligns: &'a std::cell::RefCell<Vec<(usize, Aligned)>>) -> Self {
        Examined { site, scorer, build, haps: Default::default(), aligns }
    }

    /// Candidate `k`'s haplotype.
    fn hap(&self, k: usize) -> std::rc::Rc<Hap> {
        if let Some((_, h)) = self.haps.borrow().iter().find(|(j, _)| *j == k) {
            return h.clone();
        }
        let h = (self.build)(k);
        self.haps.borrow_mut().push((k, h.clone()));
        h
    }

    /// The read's alignment to candidate `k`'s haplotype, and where its trim starts.
    fn align(&self, k: usize) -> Aligned {
        if let Some((_, a)) = self.aligns.borrow().iter().find(|(j, _)| *j == k) {
            return a.clone();
        }
        let hap = self.hap(k);
        let a = std::rc::Rc::new(self.scorer.align(self.site, &hap));
        self.aligns.borrow_mut().push((k, a.clone()));
        a
    }
}

/// A read's call from its best score per class (`best`: score and candidate), its
/// candidates examined through `ex` (built where needed): the margin decides; the fit
/// check, X's telling bases, the anchor and the deciding quality are read off the
/// alignments to the winning and best losing haplotypes.
fn conclude<'h, F: Fn(usize) -> std::rc::Rc<Hap>>(
    site: &'h Site,
    read: &ReadView,
    ex: &Examined<'h, '_, F>,
    best: [Option<(i32, usize)>; 3],
    min_baseq: u8,
    ask_anchor: bool,
) -> Verdict {
    if site.own_gap == OwnGapRule::Tie {
        return conclude_margin(site, read, ex, best, min_baseq, ask_anchor);
    }
    // Measurement only (`GBCMS_JUDGE_OWNGAP=all` or `others`): a state the read reaches only
    // through an indel of its own at that state's change is not the read's, whatever it ties
    // with or beats: each state within the margin of the best (the best first, then any the
    // margin could weigh against it; under `others` another allele's only) is dropped where
    // its alignment writes such an indel, and the margin decides among the rest.
    // The fit check (RJ-7: the best haplotype's alignment leaving an indel of the read's own
    // at X's site is another allele) stands first, on every state: a read it sends to
    // another allele carries no state there, whichever the margin would weigh after.
    if other_option().fit {
        let overall = best.iter().flatten().copied().max_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)));
        if let Some((_, h)) = overall.filter(|&(s, _)| s < 0) {
            let (candidate, al) = (ex.hap(h), ex.align(h));
            if site.gap_at_site(&candidate, &read.seq, &al.0, al.1) {
                return conclude_margin(site, read, ex, best, min_baseq, ask_anchor);
            }
        }
    }
    let score = |c: usize, dropped: &[bool; 3]| best[c].filter(|_| !dropped[c]).map(|(s, _)| s);
    let mut og = OwnGapTrace { before: Some(decide(best.map(|b| b.map(|(s, _)| s)))), ..OwnGapTrace::default() };
    let top0 = best.iter().flatten().map(|b| b.0).max();
    og.tie = top0.is_some_and(|t| best.iter().flatten().filter(|b| t - b.0 < DELTA).count() > 1);
    let mut checked = [false; 3];
    while let Some(top) = (0..3).filter_map(|c| score(c, &og.dropped)).max() {
        let mut changed = false;
        for c in 0..3 {
            let Some((s, h)) = best[c].filter(|_| !og.dropped[c] && !checked[c]) else { continue };
            // under `others`, REF and X stay for the margin to weigh against each other
            if top - s >= DELTA || (site.own_gap == OwnGapRule::Others && c != Class::Other as usize) {
                continue;
            }
            checked[c] = true;
            let (state, al) = (ex.hap(h), ex.align(h));
            if site.own_gap_at(&state, &read.seq, &al.0, al.1) {
                og.dropped[c] = true;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    if !og.dropped.iter().any(|&d| d) {
        return conclude_margin(site, read, ex, best, min_baseq, ask_anchor);
    }
    let mut kept = best;
    for (k, &d) in kept.iter_mut().zip(&og.dropped) {
        if d {
            *k = None;
        }
    }
    if kept.iter().all(Option::is_none) {
        trace!("haplotype judge: every state is reached only through an indel of the read's own → another allele");
        let scores = best.map(|b| b.map(|(s, _)| s));
        return Verdict { call: Call::Other, qual: 0, best: scores, holds_anchor: false, placed: None, confirmed: false, spells: 0, clip: None, fits: FitsTrace::default(), fits_none: false, neither: NeitherTrace::default(), owngap: og };
    }
    trace!("haplotype judge: states reached only through an indel of the read's own dropped ({:?})", og.dropped);
    Verdict { owngap: og, ..conclude_margin(site, read, ex, kept, min_baseq, ask_anchor) }
}

/// [`conclude`] by the margin alone.
fn conclude_margin<'h, F: Fn(usize) -> std::rc::Rc<Hap>>(
    site: &'h Site,
    read: &ReadView,
    ex: &Examined<'h, '_, F>,
    best: [Option<(i32, usize)>; 3],
    min_baseq: u8,
    ask_anchor: bool,
) -> Verdict {
    let scores = best.map(|b| b.map(|(s, _)| s));
    let mut call = decide(scores);
    // A tie of REF with another allele, X beaten, is evidence of that allele only where
    // the read reaches the bases that tell it from REF: short of them (an allele inside a
    // long deletion's bases, past where a read over X's anchor ends) the read shows REF
    // where REF and X differ, and nothing of the other allele.
    let mut ref_by_reach = false;
    // Nor where the read writes an indel of its own across that allele's change (an
    // unlisted deletion beside a listed one, each fitting it badly): it does not carry it,
    // and reads REF where REF and X differ.
    if call == Call::Other && tied_with(scores) == Some(Class::Ref as usize) {
        if let Some((_, k)) = best[Class::Other as usize] {
            let other = ex.hap(k);
            if !site.reaches_state(read, &other) {
                trace!("haplotype judge: REF ties another allele the read never reaches → REF");
                call = Call::Ref;
                ref_by_reach = true;
            } else {
                let al = ex.align(k);
                if site.gap_at_state(&other, &read.seq, &al.0, al.1) {
                    trace!("haplotype judge: REF ties another allele the read writes its own indel across → REF");
                    call = Call::Ref;
                    ref_by_reach = true;
                }
            }
        }
    }
    // The haplotype the read fits best of all, where its alignment there leaves an indel
    // unexplained inside X's site: the read carries another allele there, never REF or X.
    let overall = best.iter().flatten().copied().max_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)));
    let mut fit_other: Option<usize> = None;
    if other_option().fit {
        if let Some((s, h)) = overall.filter(|&(s, _)| s < 0) {
            let (candidate, al) = (ex.hap(h), ex.align(h));
            if site.gap_at_site(&candidate, &read.seq, &al.0, al.1) {
                trace!("haplotype judge: the best haplotype ({s}) leaves an indel at X's site unexplained → another allele");
                call = Call::Other;
                fit_other = Some(h);
            }
        }
    }
    // A substitution X's read whose confident bases at X's telling bases carry some of X's
    // changes and REF's bases at others carries another allele (RJ-28: at an MNV's row a
    // read carrying only some of its components), wherever the margin leaves it: REF by a
    // base (one change of three), X by a base (two of three) or a tie of the two (one of
    // two). A base below `--min-baseq` there weighs nothing either way.
    let mut mixed_on: Option<usize> = None;
    if call == Call::Depth && site.changed_bases() >= 2 {
        if let (Some((a, _)), Some((r, h))) = (best[Class::Alt as usize], best[Class::Ref as usize]) {
            if (a - r).abs() < DELTA {
                let (candidate, al) = (ex.hap(h), ex.align(h));
                if read_at_x(site, read, &candidate, &al.0, al.1).mixed() {
                    trace!("haplotype judge: REF and X tie on a read carrying some of X's changes → another allele");
                    call = Call::Other;
                    mixed_on = Some(h);
                }
            }
        }
    }
    let won = match call {
        Call::Alt => Some(Class::Alt),
        Call::Ref => Some(Class::Ref),
        Call::Other => Some(Class::Other),
        Call::Depth => None,
    };
    let Some(won) = won else {
        return Verdict { call, qual: 0, best: scores, holds_anchor: false, placed: None, confirmed: false, spells: 0, clip: None, fits: FitsTrace::default(), fits_none: false, neither: NeitherTrace::default(), owngap: OwnGapTrace::default() };
    };
    if call == Call::Other && !ask_anchor {
        return Verdict { call, qual: 0, best: scores, holds_anchor: false, placed: None, confirmed: false, spells: 0, clip: None, fits: FitsTrace::default(), fits_none: false, neither: NeitherTrace::default(), owngap: OwnGapTrace::default() };
    }
    let w = fit_other.or(mixed_on).unwrap_or_else(|| best[won as usize].map_or(0, |(_, h)| h));
    let (hap_w, w_al) = (ex.hap(w), ex.align(w));
    let hap_w: &Hap = &hap_w;
    let (wa, off) = (&w_al.0, w_al.1);
    // Measurement only (GBCMS_JUDGE_NEITHER): a covered telling base read confidently as
    // neither allele makes a REF or ALT call another allele, at every size (`other`), below a
    // size (`other:N`), or where other reads of the row read it too (`recur:K`).
    let mut nt = NeitherTrace::default();
    if matches!(call, Call::Alt | Call::Ref) && (site.neither != NeitherRule::Margin || neither_trace()) {
        if let Some(t) = neither_bases(site, read, hap_w, wa, off) {
            nt = NeitherTrace { on: neither_trace(), call: if call == Call::Alt { 1 } else { 0 }, ..t };
            nt.moved = nt.n > 0
                && match site.neither {
                    NeitherRule::Margin => false,
                    NeitherRule::Other => true,
                    NeitherRule::Below(n) => site.x_size() < n,
                    NeitherRule::Recur(k) => nt.max_recur >= k,
                    NeitherRule::Window => nt.window > 0,
                };
        }
    }
    if nt.moved {
        trace!("haplotype judge: a {call:?} call on a read with a telling base read as neither allele → another allele");
        call = Call::Other;
        if !ask_anchor {
            return Verdict { call, qual: 0, best: scores, holds_anchor: false, placed: None, confirmed: false, spells: 0, clip: None, fits: FitsTrace::default(), fits_none: false, neither: nt, owngap: OwnGapTrace::default() };
        }
    }
    if matches!(call, Call::Alt | Call::Ref) && site.changed_bases() >= 2 && read_at_x(site, read, hap_w, wa, off).mixed() {
        trace!("haplotype judge: a {call:?} call on a read carrying some of X's changes → another allele");
        call = Call::Other;
        if !ask_anchor {
            return Verdict { call, qual: 0, best: scores, holds_anchor: false, placed: None, confirmed: false, spells: 0, clip: None, fits: FitsTrace::default(), fits_none: false, neither: NeitherTrace::default(), owngap: OwnGapTrace::default() };
        }
    }
    // Measurement only (GBCMS_JUDGE_REFBASE=partial): at a delins, a REF or ALT call on a
    // read showing the other allele's base at one of the called allele's telling bases is
    // another allele (RJ-25: an allele lacking part of the event), as at an MNP (RJ-28).
    if site.ref_base == RefBaseRule::Partial
        && matches!(call, Call::Alt | Call::Ref)
        && !site.is_substitution()
        && other_allele_base(site, read, hap_w, wa, off)
    {
        trace!("haplotype judge: a {call:?} call on a read showing the other allele's base at a telling base → another allele");
        call = Call::Other;
        if !ask_anchor {
            return Verdict { call, qual: 0, best: scores, holds_anchor: false, placed: None, confirmed: false, spells: 0, clip: None, fits: FitsTrace::default(), fits_none: false, neither: nt, owngap: OwnGapTrace::default() };
        }
    }
    let column = |p: &Placed| if let Placed::At(c) = *p { Some((c + off) as i64) } else { None };
    let holds_anchor = wa.covers().is_some_and(|(a, b)| {
        let anchor = site.anchor_in(hap_w);
        ((a + off) as i64) <= anchor && anchor <= (b + off) as i64
    });
    let at = |p: &Placed| matches!(p, Placed::At(_));
    let placed = wa.placed.iter().position(at).zip(wa.placed.iter().rposition(at)).map(|(a, b)| (a, b + 1));
    if call == Call::Other {
        return Verdict { call, qual: 0, best: scores, holds_anchor, placed, confirmed: false, spells: 0, clip: None, fits: FitsTrace::default(), fits_none: false, neither: nt, owngap: OwnGapTrace::default() };
    }
    // The bases that tell X from REF at its site, covered by the winning alignment: X's
    // own for an ALT call, none of them read confidently as another base; REF's for a REF
    // call on a substitution (the rules read an MNP's whole block); never on a delins (the
    // rules read equal-length windows there). A base below `--min-baseq` there weighs
    // nothing either way, as an MNP read with a differing base below it counts ALT under
    // the rules (RJ-25); the margin already rests the call on its confident bases.
    let need = site.is_substitution();
    let telling = if need { site.telling_in(hap_w) } else { Vec::new() };
    // Measurement only (GBCMS_JUDGE_REFBASE=partial): REF's base at each of them, where a
    // third base (neither REF's nor X's) is an error the margin forgives (RJ-26).
    let ref_at: Vec<(i64, u8)> = if site.ref_base == RefBaseRule::Partial && site.is_substitution() {
        site.telling_pairs(hap_w).into_iter().map(|(c, rb, _)| (c, rb)).collect()
    } else {
        Vec::new()
    };
    let mut confirmed = call == Call::Alt;
    let mut reads_all = true;
    let mut telling_at = Vec::new();
    for (col, base) in telling {
        match wa.placed.iter().position(|p| column(p) == Some(col)) {
            None => reads_all = false,
            Some(i) if call == Call::Alt => {
                telling_at.push(i);
                let b = read.seq[i];
                confirmed &= b == base;
                let third = b != base && ref_at.iter().any(|&(c, rb)| c == col && rb != b);
                reads_all &= b == WILD || b == base || third;
            }
            Some(i) => telling_at.push(i),
        }
    }
    let clip = if ask_anchor && site.is_substitution() && reads_all { clip_bases(read, hap_w, wa, off, &telling_at) } else { None };
    if !reads_all {
        call = Call::Depth;
        return Verdict { call, qual: 0, best: scores, holds_anchor: false, placed: None, confirmed: false, spells: 0, clip: None, fits: FitsTrace::default(), fits_none: false, neither: nt, owngap: OwnGapTrace::default() };
    }
    // The best losing haplotype: the best of the other classes, the first on a tie (X's,
    // for a REF call the other allele tied: the bases that beat X decided it).
    let floor = min_baseq.max(1);
    let lost = (0..3)
        .filter(|&k| k != won as usize && !(ref_by_reach && k == Class::Other as usize))
        .filter_map(|k| best[k])
        .max_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)))
        .map(|(_, h)| h);
    let qual = lost.map_or(floor, |l| {
        let la = ex.align(l);
        deciding_quality(read, wa, &la.0).map_or(floor, |q| q.max(floor))
    });
    let spells = if call == Call::Alt { site.spelled_holders(read, hap_w, wa, off) } else { 0 };
    Verdict { call, qual, best: scores, holds_anchor, placed, confirmed, spells, clip, fits: FitsTrace::default(), fits_none: false, neither: nt, owngap: OwnGapTrace::default() }
}

/// What a read's confident bases show at a substitution X's telling bases, on its
/// alignment `al` (columns offset by `off`) to `hap` (carrying X or holding REF's bases):
/// at a base X changes, X's base, REF's or a third; where the two alleles agree (a base
/// inside an MNV that it keeps), a base other than theirs is a third.
#[derive(Debug, Clone, Copy, Default)]
struct AtX {
    x: bool,
    reference: bool,
    third: bool,
}

impl AtX {
    /// Some of X's changes read and REF's bases at others: another allele.
    fn mixed(self) -> bool {
        self.x && self.reference
    }
}

fn read_at_x(site: &Site, read: &ReadView, hap: &Hap, al: &Alignment, off: usize) -> AtX {
    let mut out = AtX::default();
    for (col, rb, ab) in site.telling_pairs(hap) {
        let Some(i) = al.placed.iter().position(|p| matches!(*p, Placed::At(c) if (c + off) as i64 == col)) else {
            continue;
        };
        let b = read.seq[i];
        if b == WILD {
            continue;
        }
        if rb != ab && b == ab {
            out.x = true;
        } else if b == rb {
            out.reference |= rb != ab;
        } else {
            out.third = true;
        }
    }
    out
}

/// Measurement only (`GBCMS_JUDGE_REFBASE=partial`): whether `read`, on its alignment `al`
/// to `hap` (carrying a delins X or holding REF's bases there; columns offset by `off`),
/// reads confidently, at one of that allele's telling bases, the other allele's base at the
/// same telling index from either end: that allele's evidence, not an error. A base the
/// haplotype holds there (a context it applies) or the guard's background explains is
/// neither.
fn other_allele_base(site: &Site, read: &ReadView, hap: &Hap, al: &Alignment, off: usize) -> bool {
    site.telling_with_other(hap).into_iter().any(|(col, _, other)| {
        let Some(i) = al.placed.iter().position(|p| matches!(*p, Placed::At(c) if (c + off) as i64 == col)) else {
            return false;
        };
        let (got, held) = (read.seq[i], hap.base(col as usize));
        got != WILD
            && held != WILD
            && got != held
            && !read.own_position(i).is_some_and(|p| site.background_explains(p, got))
            && other.contains(&Some(got))
    })
}

/// Measurement only: what `read`'s alignment `al` to `hap` (a haplotype carrying X or holding
/// REF's bases there, its columns offset by `off`) reads at X's telling bases that is neither
/// allele's: on a substitution of two bases or more, a confident base neither REF's nor X's
/// at a base X changes, or other than theirs where they agree; on a delins, a confident base
/// other than the called allele's; never the haplotype's own base there nor one the guard's
/// background explains (a masked base takes any). None at any other row.
fn neither_bases(site: &Site, read: &ReadView, hap: &Hap, al: &Alignment, off: usize) -> Option<NeitherTrace> {
    let telling: Vec<(i64, Option<u8>, u8)> = if site.changed_bases() >= 2 {
        site.telling_pairs(hap).into_iter().map(|(col, rb, ab)| (col, Some(rb), ab)).collect()
    } else if !site.is_substitution() {
        site.telling_in(hap).into_iter().map(|(col, b)| (col, None, b)).collect()
    } else {
        return None;
    };
    if telling.is_empty() {
        return None;
    }
    let at = |col: i64| al.placed.iter().position(|p| matches!(*p, Placed::At(c) if (c + off) as i64 == col));
    let x_pos = site.x_pos();
    let mut t = NeitherTrace { size: site.x_size() as u16, telling: telling.len() as u16, ..NeitherTrace::default() };
    for (k, &(col, rb, b)) in telling.iter().enumerate() {
        let Some(i) = at(col) else { continue };
        let (got, held) = (read.seq[i], hap.base(col as usize));
        if got == WILD || held == WILD {
            continue;
        }
        t.conf += 1;
        let pos = read.own_position(i);
        // A base the read's haplotype holds there (a context it applies), or the guard's
        // background explains, is no third allele's.
        let neither = got != held
            && !pos.is_some_and(|p| site.background_explains(p, got))
            && match rb {
                Some(rb) if rb != b => got != rb && got != b,
                Some(rb) => got != rb,
                None => got != b,
            };
        if neither {
            // a base the alleles both read: one a substitution changes, or a length change's
            // telling base within the junction window of either end
            let jw = super::haplotypes::JUNCTION_WINDOW as usize;
            let alike = match rb {
                Some(rb) => rb != b,
                None => k < jw || k + jw >= telling.len(),
            };
            t.window = t.window.saturating_add(u8::from(alike));
            let recur = pos.map_or(0, |p| site.recurrence(p, got));
            t.max_recur = t.max_recur.max(recur);
            if (t.n as usize) < t.bases.len() {
                let off_x = pos.map_or(i32::MIN, |p| (p - x_pos) as i32);
                t.bases[t.n as usize] = (k as u16, got, read.quals[i], recur.min(u32::from(u16::MAX)) as u16, off_x);
            }
            t.n = t.n.saturating_add(1);
        }
    }
    let cols: Vec<i64> = telling.iter().map(|t| t.0).collect();
    for (i, p) in al.placed.iter().enumerate() {
        if let Placed::At(c) = *p {
            let col = (c + off) as i64;
            let (got, want) = (read.seq[i], hap.base(c + off));
            if !cols.contains(&col) && got != WILD && want != WILD {
                t.bg_conf += 1;
                if got != want {
                    t.bg_mm += 1;
                    t.bg_single += u32::from(read.own_position(i).is_some_and(|p| site.recurrence(p, got) == 0));
                }
            }
        }
    }
    Some(t)
}

/// What one read base costs in an alignment, and the haplotype bases skipped on each
/// side of it.
fn base_costs(a: &Alignment) -> Vec<(i32, u32, u32)> {
    (0..a.costs.len())
        .map(|i| {
            let before = if i > 0 { a.skipped_after[i - 1] } else { 0 };
            (a.costs[i], before, a.skipped_after[i])
        })
        .collect()
}

/// The lowest quality among the read's confident bases where the winning and the best
/// losing haplotype's alignments differ: the bases that decided the call. None when no
/// confident base differs.
fn deciding_quality(read: &ReadView, won: &Alignment, lost: &Alignment) -> Option<u8> {
    let (w, l) = (base_costs(won), base_costs(lost));
    (0..read.seq.len()).filter(|&i| w[i] != l[i] && read.seq[i] != WILD).map(|i| read.quals[i]).min()
}

/// The read as the judge reads it, or None: a record whose CIGAR splices (judged depth
/// only for now), or one with no readable base.
pub(crate) fn read_view(record: &Record, quals: &[u8], min_baseq: u8, rules: &ReadRules) -> Option<ReadView> {
    let cigar = record.cigar();
    // TODO(step 6 of the judge's build): judge a spliced read on the reference spliced
    // along its junctions (`carrier::splice_segments`), as the exact-carrier rule does;
    // until then it is depth only.
    if cigar.iter().any(|op| matches!(op, Cigar::RefSkip(_))) {
        return None;
    }
    let seq = record.seq().as_bytes();
    let readable = readable_range(record, rules)?;
    let tag = |t: &[u8]| match record.aux(t) {
        Ok(Aux::I8(v)) => Some(i64::from(v)),
        Ok(Aux::U8(v)) => Some(i64::from(v)),
        Ok(Aux::I16(v)) => Some(i64::from(v)),
        Ok(Aux::U16(v)) => Some(i64::from(v)),
        Ok(Aux::I32(v)) => Some(i64::from(v)),
        Ok(Aux::U32(v)) => Some(i64::from(v)),
        _ => None,
    };
    let tags = AlnTags { mapq: record.mapq(), best: tag(b"AS"), elsewhere: tag(b"XS"), nm: tag(b"NM") };
    ReadView::new(&seq, quals, record.pos(), &cigar, readable, min_baseq).map(|v| v.named(record.qname()).tagged(tags))
}

/// The haplotype judge's call for one read at `variant` (after the splice-skip triage),
/// as the engine counts it:
/// - ALT / REF: `is_alt` / `is_ref`, with `qual` the lowest quality among the bases
///   that decided it (at least `--min-baseq`), so mates resolve by the bases a call
///   rests on;
/// - OTHER: neither, partial evidence (`partial_match_count = 1`);
/// - depth only: neither, `uninformative`.
///
/// `has_n_base` where the read shows an N over the row's REF span, and an MNP's ALT call
/// `mnp_confirmed` where it read every base the alleles differ at, as the rules set them.
///
/// Admission keeps DP's meaning: a read the genome places off the row's anchor counts
/// only where its alignment to the haplotype it was called on holds the anchor there
/// (an insertion's carrier whose clip holds inserted bases alone does not, nor a read
/// starting inside a deletion's bases): `clip_admissible` where its soft clip is read
/// (the engine admits it in DNA), else `read_from_right`. The verdict is kept per read,
/// so the per-transcript counts and ASJD reuse the main counts' call.
pub(crate) fn classify(record: &Record, variant: &Variant, quals: &[u8], min_baseq: u8, rules: &ReadRules) -> ClassifyResult {
    let depth = || ClassifyResult { uninformative: true, ..ClassifyResult::neither(PHASE) };
    if record.seq_len() == 0 || quals.len() < record.seq_len() {
        // A record stored without its bases shows nothing here.
        return ClassifyResult::no_coverage(PHASE);
    }
    let Some(haplotypes) = rules.haplotypes else {
        return depth();
    };
    let Some(site) = haplotypes.get(variant, rules) else {
        return depth();
    };
    let key = (
        record.qname().to_vec(),
        record.flags(),
        record.pos(),
        [record.raw_cigar().iter().flat_map(|op| op.to_le_bytes()).collect::<Vec<u8>>(), record.seq().encoded.to_vec()].concat(),
    );
    if let Some(kept) = haplotypes.verdict(&key) {
        return kept;
    }
    let (start, end) = (record.pos(), ref_end(record));
    let overlaps_anchor = start <= variant.pos && end > variant.pos;
    let out = match read_view(record, quals, min_baseq, rules) {
        None => depth(),
        Some(mut read) => {
            // The judge never reads a read's name: a read alike in all else (bases, qualities,
            // placement, the aligner's tags) to one this row judged gets that read's verdict.
            read.name.clear();
            // the per-read trace lines' arguments are built only where they are written
            let tracing = log::log_enabled!(log::Level::Trace);
            let view = (read, !overlaps_anchor, min_baseq);
            let clock = std::time::Instant::now();
            let judged = haplotypes.judged(&view);
            let verdict = judged.unwrap_or_else(|| {
                // what a read alike in bases (masked) and placement left: its programmes
                let placed = view.0.placed_alone();
                let kept = haplotypes.scored(&placed);
                let fresh = kept.is_none();
                let (verdict, scored) = judge_kept(site, &view.0, min_baseq, !overlaps_anchor, kept);
                if let Some(s) = scored.filter(|_| fresh) {
                    haplotypes.keep_scored(placed, s);
                }
                verdict
            });
            let ns = clock.elapsed().as_nanos();
            let read = &view.0;
            // X's placements touch a contig edge (its repeat runs to it, or X sits at the
            // contig's first or last base) and the read reaches it: an ALT call there
            // rests on REF holding no base past the edge, which a circular contig would
            // hold (RJ-3), so the read is depth only, counted and warned. A REF call rests
            // on REF's own bases inside the contig and stands. The judge's verdict is kept
            // as it is for the reads alike (`Haplotypes::judged`); the edge is weighed on
            // every read, alike reads alike (one placement).
            let edge_alt = verdict.call == Call::Alt && site.edge_reached(start, end);
            let call = if edge_alt { Call::Depth } else { verdict.call };
            let mut out = match call {
                Call::Alt => ClassifyResult::is_alt(verdict.qual, PHASE),
                Call::Ref => ClassifyResult::is_ref(verdict.qual, PHASE),
                Call::Other => ClassifyResult::neither_with_partial(PHASE, 1, false),
                Call::Depth => depth(),
            };
            // An MNP's ALT call that read every base where the alleles differ, confidently
            // and as the ALT's (MNP rescue keeps such a row's haplotype).
            out.mnp_confirmed = out.is_alt && verdict.confirmed;
            // A read no candidate sequence explains, counted or not (UNEXPLAINED_READS).
            out.fits_none = verdict.fits_none;
            out.alt_unjudged = edge_alt;
            if out.is_alt && verdict.spells != 0 {
                out.spells_listed = spelled_listed(site, variant, verdict.spells);
            }
            // Admission off the anchor, by the alignment to the called haplotype. A clip
            // admits the read only where every base that alignment places lies inside the
            // read's fragment (RJ-13): past the fragment end a clip is adapter, and a read
            // with no fragment (unpaired, or its pair not facing inward) has none to bound
            // it.
            if !overlaps_anchor && call != Call::Depth && verdict.holds_anchor {
                // Measurement only, TEMPORARY: GBCMS_JUDGE_CLIPSUB at a substitution row.
                let sub = site.is_substitution();
                let option = clip_sub();
                let allowed = !sub
                    || match option {
                        ClipSub::Admit => true,
                        ClipSub::Off => false,
                        ClipSub::Anchored(k) => verdict.clip.is_none_or(|c| c.beyond.0 >= k),
                        ClipSub::Clean(k) => verdict.clip.is_none_or(|c| c.beyond.0 >= k && c.beyond.1 == 0 && c.between.1 == 0),
                    };
                let (lead, tail) = soft_clips(record);
                let mut inside = false;
                if lead + tail > 0 {
                    let first = readable_range(record, rules).map_or(0, |r| r.0);
                    inside = verdict.placed.zip(fragment_query_span(record)).is_some_and(|((a, b), (lo, hi))| {
                        first + a >= lo && first + b <= hi
                    });
                    out.clip_admissible = inside && allowed;
                } else {
                    out.read_from_right = allowed;
                }
                if sub && tracing {
                    let c = verdict.clip;
                    trace!(
                        "judge clipsub {}:{} {}>{} read={} mate={} call={:?} soft=({},{}) clip={} inside={} admitted={} \
                         beyond={:?} between={:?} in_clip={:?} option={:?}",
                        variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele,
                        String::from_utf8_lossy(record.qname()), if record.is_first_in_template() { 1 } else { 2 },
                        call, lead, tail,
                        c.map_or("none", |c| if c.tail { "tail" } else { "lead" }),
                        inside, out.clip_admissible || out.read_from_right,
                        c.map(|c| c.beyond), c.map(|c| c.between), c.map(|c| c.in_clip), option,
                    );
                }
            }
            if tracing {
                trace!(
                    "haplotype judge {}:{} {}>{} read={}: {:?} (best alt={:?} ref={:?} other={:?}) qual={} \
                     holds the anchor={} over {} haplotypes",
                    variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele,
                    String::from_utf8_lossy(record.qname()), call, verdict.best[0], verdict.best[1],
                    verdict.best[2], out.qual, overlaps_anchor || verdict.holds_anchor, site.haps.len(),
                );
            }
            // Measurement only, the "fits nothing" decision package: per judged read, its best
            // score per class (kept where the rule made it depth), its confident bases and the
            // errors its qualities lead one to expect, the states the full programme decided,
            // the time the judge took, the aligner's tags and the neighbourhood fit.
            let nt = verdict.neither;
            if nt.on {
                let bases: Vec<String> = nt.bases[..usize::from(nt.n).min(4)]
                    .iter()
                    .map(|&(k, b, q, r, o)| format!("{k}:{}:{q}:{r}:{}", b as char, if o == i32::MIN { "-".to_string() } else { o.to_string() }))
                    .collect();
                trace!(
                    "judge neither {}:{} {}>{} read={} mate={} side={} size={} telling={} conf={} n={} bases={} bg={}/{}/{} moved={} final={:?}",
                    variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele,
                    String::from_utf8_lossy(record.qname()), if record.is_first_in_template() { 1 } else { 2 },
                    if nt.call == 1 { "ALT" } else { "REF" }, nt.size, nt.telling, nt.conf, nt.n,
                    if bases.is_empty() { "-".to_string() } else { bases.join(";") }, nt.bg_single, nt.bg_mm, nt.bg_conf, nt.moved, call,
                );
            }
            let og = verdict.owngap;
            if let Some(before) = og.before.filter(|_| og.dropped.iter().any(|&d| d)) {
                let names = ["ALT", "REF", "OTHER"];
                let dropped: Vec<&str> = (0..3).filter(|&c| og.dropped[c]).map(|c| names[c]).collect();
                trace!(
                    "judge owngap {}:{} {}>{} read={} mate={} before={:?} after={:?} dropped={} tie={}",
                    variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele,
                    String::from_utf8_lossy(record.qname()), if record.is_first_in_template() { 1 } else { 2 },
                    before, call, dropped.join(","), og.tie,
                );
            }
            let f = verdict.fits;
            if tracing {
                trace!(
                    "judge fits {}:{} {}>{} read={} mate={} call={:?} alt={} ref={} other={} conf={} len={} experr={:.2} \
                     full={}/{} ns={} mapq={} as={} xs={} nm={} nb={} fired={} rule={:?}",
                    variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele,
                    String::from_utf8_lossy(record.qname()), if record.is_first_in_template() { 1 } else { 2 },
                    call, opt(verdict.best[0]), opt(verdict.best[1]), opt(verdict.best[2]), read.confident(), read.seq.len(),
                    read.expected_errors(), f.full, f.states, ns, read.tags.mapq, opt(read.tags.best), opt(read.tags.elsewhere),
                    opt(read.tags.nm), opt(f.neighbourhood), f.fired, fits_rule(),
                );
            }
            if judged.is_none() {
                haplotypes.keep_judged(view, verdict);
            }
            out
        }
    };
    let out = mark(out, record, variant);
    haplotypes.keep_verdict(key, out);
    out
}

/// The co-listed rows holding `variant`'s change (`super::listed_holders`) among the site's
/// holders a read spells (`spells`, a bit per [`Site::holders`]), as bits in their order.
fn spelled_listed(site: &Site, variant: &Variant, spells: u32) -> u32 {
    let listed = super::listed_holders(variant);
    site.holders()
        .into_iter()
        .take(32)
        .enumerate()
        .filter(|&(k, _)| spells >> k & 1 == 1)
        .filter_map(|(_, h)| {
            listed.iter().position(|&(p, r, a)| *p == h.pos && r.as_bytes().eq_ignore_ascii_case(&h.ref_al) && a.as_bytes().eq_ignore_ascii_case(&h.alt_al))
        })
        .filter(|&j| j < 32)
        .fold(0u32, |m, j| m | 1 << j)
}

/// A value for the trace, `-` where absent.
fn opt<T: std::fmt::Display>(v: Option<T>) -> String {
    v.map_or_else(|| "-".to_string(), |x| x.to_string())
}

/// The flags a call carries beside its allele: an N the read shows over the row's REF
/// span (the duplex-masking diagnostic); an MNP ALT call that read every base where the
/// alleles differ, confidently and as the ALT's (MNP rescue keeps such a row's
/// haplotype); and, measurement only, the CIGAR's exact indel.
fn mark(mut out: ClassifyResult, record: &Record, variant: &Variant) -> ClassifyResult {
    let kind = allele_kind(&variant.ref_allele, &variant.alt_allele);
    let ref_len = variant.ref_allele.len() as i64;
    let long_deletion = kind == Some(AlleleKind::Deletion)
        && variant.ref_allele.len() - variant.alt_allele.len() >= crate::counting::junction::LONG_DELETION;
    let span = if long_deletion { (variant.pos, variant.pos + 1) } else { (variant.pos, variant.pos + ref_len) };
    let seq = record.seq();
    out.has_n_base = (span.0..span.1).any(|p| find_read_pos(record, p).is_some_and(|q| seq[q] == b'N'));
    // An MNP's ALT call that read every base where the alleles differ, confidently and as
    // the ALT's, on the alignment it was called on (however the CIGAR wrote them).
    out.mnp_confirmed &= out.is_alt && kind == Some(AlleleKind::Mnp);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::counting::judge::haplotypes::tests::{random_seq, row, site};
    use crate::counting::judge::{Haplotypes, ReadJudge, Row};
    use rust_htslib::bam::record::CigarString;

    fn record(seq: &[u8], pos: i64, cigar: Vec<Cigar>) -> Record {
        let mut rec = Record::new();
        rec.set(b"r", Some(&CigarString(cigar)), seq, &vec![30; seq.len()]);
        rec.set_pos(pos);
        rec
    }

    /// `rec` as the forward read of a pair facing inward on contig 0, its mate 150 bp
    /// on: a defined fragment (TLEN 300) that holds its leading clip.
    fn in_pair(mut rec: Record) -> Record {
        rec.set_tid(0);
        rec.set_mtid(0);
        rec.set_mpos(rec.pos() + 150);
        rec.set_insert_size(300);
        rec.set_flags(0x1 | 0x2 | 0x20 | 0x40);
        rec
    }

    fn variant(r: &Row) -> Variant {
        let text = |b: &[u8]| String::from_utf8(b.to_vec()).unwrap();
        Variant::new("1".into(), r.pos, text(&r.ref_al), text(&r.alt_al), "X".into(), None, 0, 0, None, None, None, None)
    }

    fn applied(contig: &[u8], r: &Row) -> Vec<u8> {
        [&contig[..r.pos as usize], &r.alt_al[..], &contig[r.pos as usize + r.ref_al.len()..]].concat()
    }

    /// `f` with `x`'s haplotypes on `contig` in the rules.
    fn with_rules(contig: &[u8], x: &Row, f: impl FnOnce(&Variant, &ReadRules)) {
        let s = site(contig, x, &[], 100);
        let v = variant(x);
        let haps = Haplotypes::new(&v, move |_| Some(s.clone()));
        let rules = ReadRules { judge: ReadJudge::Haplotype, haplotypes: Some(&haps), ..ReadRules::DNA };
        f(&v, &rules);
    }

    #[test]
    fn a_call_maps_onto_the_engines_result_and_admission() {
        let c = random_seq(600, 17);
        let x = row(301, &(c[300] as char).to_string(), if c[300] == b'A' { "C" } else { "A" });
        let h = applied(&c, &x);
        with_rules(&c, &x, |v, rules| {
            let call = |rec: &Record| classify(rec, v, rec.qual(), 20, rules);
            let alt = call(&record(&h[250..350], 250, vec![Cigar::Match(100)]));
            assert!(alt.is_alt && alt.qual == 30 && alt.phase == PHASE && !alt.uninformative);
            assert!(!alt.read_from_right && !alt.clip_admissible && !alt.is_structural);
            let reference = call(&record(&c[251..351], 251, vec![Cigar::Match(100)]));
            assert!(reference.is_ref && !reference.is_alt);
            // X in the leading clip of a read paired inward, aligned from the base after
            // it: its alignment to X's haplotype holds the anchor, by its clipped bases.
            let clipped = call(&in_pair(record(&h[281..381], 301, vec![Cigar::SoftClip(20), Cigar::Match(80)])));
            assert!(clipped.is_alt && clipped.clip_admissible && !clipped.read_from_right);
            // Unpaired, the same read has no fragment to bound its clip: not admitted.
            let unpaired = call(&record(&h[281..381], 301, vec![Cigar::SoftClip(20), Cigar::Match(80)]));
            assert!(!unpaired.clip_admissible, "{unpaired:?}");
            // A read starting right after X holds no anchor on any haplotype.
            let after = call(&record(&c[301..401], 301, vec![Cigar::Match(100)]));
            assert!(!after.clip_admissible && !after.read_from_right, "{after:?}");
            // A read ending before X: depth only.
            let short = call(&record(&c[190..290], 190, vec![Cigar::Match(100)]));
            assert!(!short.is_ref && !short.is_alt && short.uninformative);
            // A third base at the SNV: another allele.
            let mut third = h[250..350].to_vec();
            third[50] = *b"GT".iter().find(|&&b| b != x.ref_al[0] && b != x.alt_al[0]).unwrap();
            let other = call(&record(&third, 250, vec![Cigar::Match(100)]));
            assert!(!other.is_ref && !other.is_alt && other.partial_match_count == 1 && !other.uninformative);
        });
    }

    #[test]
    fn a_deletion_admits_reads_off_its_anchor_only_where_they_hold_it() {
        let c = random_seq(600, 19);
        for len in [10usize, 60] {
            let x = Row::new(250, &c[250..251 + len], &c[250..251]);
            let h = applied(&c, &x);
            with_rules(&c, &x, |v, rules| {
                // A REF read starting inside the deleted bases: REF's, but no anchor held.
                let inside = record(&c[255..355], 255, vec![Cigar::Match(100)]);
                let r = classify(&inside, v, inside.qual(), 20, rules);
                assert!(r.is_ref && !r.read_from_right && !r.clip_admissible, "{len}bp: inside the deleted bases");
                // A carrier aligned after the junction, the anchor and its flank clipped.
                let clipped = in_pair(record(&h[230..330], 251 + len as i64, vec![Cigar::SoftClip(21), Cigar::Match(79)]));
                let r = classify(&clipped, v, clipped.qual(), 20, rules);
                assert!(r.is_alt && r.clip_admissible, "{len}bp: the clip holds the anchor");
            });
        }
    }

    #[test]
    fn a_read_fitting_no_candidate_sequence_is_depth_only_only_where_its_aligner_doubts_it() {
        use crate::counting::judge::haplotypes::AlnTags;
        // The fit rule (xs, the default): an SNV's carrier with twelve confident mismatches
        // away from X fits none of the site's candidate sequences within ten. Placed
        // elsewhere as well by its aligner (XS at least AS) it is depth only; placed here
        // alone, its bases at X decide it. Either way it is a read no candidate sequence
        // explains (UNEXPLAINED_READS); a clean carrier is not.
        assert_eq!(parse_fits(""), FitsRule::Xs);
        assert_eq!(parse_fits("xs"), FitsRule::Xs);
        assert_eq!(parse_fits("exact"), FitsRule::Exact);
        let c = random_seq(600, 23);
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let x = row(301, &(c[300] as char).to_string(), &(other(c[300]) as char).to_string());
        let h = applied(&c, &x);
        let s = site(&c, &x, &[], 100);
        let mut noisy = h[250..350].to_vec();
        for i in (0..100).step_by(8).take(12) {
            noisy[i] = other(noisy[i]);
        }
        let view = |seq: &[u8], tags: AlnTags| {
            ReadView::new(seq, &[30; 100], 250, &[Cigar::Match(100)], (0, 100), 20).unwrap().tagged(tags)
        };
        let doubted = AlnTags { mapq: 0, best: Some(52), elsewhere: Some(52), nm: Some(13) };
        let sure = AlnTags { mapq: 60, best: Some(52), elsewhere: Some(20), nm: Some(13) };
        let got = |seq: &[u8], tags| {
            let v = judge(&s, &view(seq, tags), 20, false);
            (v.call, v.fits_none)
        };
        assert_eq!(got(&noisy, doubted), (Call::Depth, true), "the aligner doubts it");
        assert_eq!(got(&noisy, sure), (Call::Alt, true), "placed here alone: its bases at X");
        assert_eq!(got(&h[250..350], sure), (Call::Alt, false), "a clean carrier");
        // No XS written: MAPQ 0 stands for the aligner's doubt.
        let no_xs = AlnTags { mapq: 0, best: Some(52), elsewhere: None, nm: Some(13) };
        assert_eq!(got(&noisy, no_xs), (Call::Depth, true));
    }

    #[test]
    fn the_xs_and_scaled_option_is_depth_only_where_both_rules_would_make_it_so() {
        use crate::counting::judge::haplotypes::AlnTags;
        // GBCMS_JUDGE_FITS=xsscaled:R (measurement only): a read is depth only where the
        // aligner doubts it (xs) and its misfit passes both ten confident mismatches and R
        // per 100 confident bases (scaled). Twelve mismatches on a 100-base read pass ten
        // but not 14 per 100: judged by its bases. Sixteen pass both: depth only where
        // doubted, judged where the aligner places it here alone.
        assert_eq!(parse_fits("xsscaled:14"), FitsRule::XsScaled(14.0));
        assert_eq!(parse_fits("xsscaled"), FitsRule::XsScaled(SCALED_RATE));
        let c = random_seq(600, 29);
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let x = row(301, &(c[300] as char).to_string(), &(other(c[300]) as char).to_string());
        let h = applied(&c, &x);
        let s = site(&c, &x, &[], 100);
        let noisy = |n: usize| {
            let mut r = h[250..350].to_vec();
            for i in (0..100).step_by(6).filter(|&i| i != 50).take(n) {
                r[i] = other(r[i]);
            }
            r
        };
        let view = |seq: &[u8], tags: AlnTags| ReadView::new(seq, &[30; 100], 250, &[Cigar::Match(100)], (0, 100), 20).unwrap().tagged(tags);
        let doubted = AlnTags { mapq: 0, best: Some(40), elsewhere: Some(40), nm: Some(16) };
        let sure = AlnTags { mapq: 60, best: Some(40), elsewhere: Some(10), nm: Some(16) };
        let fired = |rule: FitsRule, seq: &[u8], tags: AlnTags| {
            let read = view(seq, tags);
            let top = judge(&s, &read, 20, false).best.iter().flatten().copied().max().unwrap();
            fits_nothing(rule, &s, &read, top).0
        };
        let both = FitsRule::XsScaled(14.0);
        assert!(fired(FitsRule::Xs, &noisy(12), doubted), "xs: past ten mismatches, doubted");
        assert!(!fired(both, &noisy(12), doubted), "within 14 per 100: judged");
        assert!(fired(both, &noisy(16), doubted), "past both, doubted: depth only");
        assert!(!fired(both, &noisy(16), sure), "placed here alone: judged");
    }
}
