//! The sequence judge (`--read-judge sequence`, #304; `haplotype` is its former name):
//! one read judgment for every small-variant row. It resolves which local sequence a
//! read's own bases spell at the row: REF, the row's ALT, or another allele at its site.
//! The read is scored whole against the row's candidate sequences (every state of the
//! row's site, with every context of the other listed rows, the guard's background and
//! the found alleles within read reach) and counts ALT or REF only by one confident base;
//! otherwise it is depth only. Contexts only explain the read's other differences: no
//! linkage, no phasing. Design: `docs/proposals/sequence-judge.md`. Internal names
//! (haplotypes, `ReadJudge::Haplotype`) stay until the module is restructured.
//!
//! - [`graph`]: the scorer, a read against each state's sequence graph holding every
//!   context at once (the default).
//! - [`score`]: a read against one haplotype (the enumerated scorer behind
//!   `GBCMS_JUDGE_SCORER=enumerated`, and the alignments a call is examined on).
//! - [`haplotypes`]: the window, the site's states and the contexts, built once per row.
//! - [`decide`]: a read's call, mapped onto `ClassifyResult`.
//! - [`triage`]: coloured k-mer triage in front of the graph scorer (off).

mod decide;
mod graph;
mod haplotypes;
mod score;
mod triage;

use std::cell::{OnceCell, RefCell};
use std::collections::HashMap;

use log::{debug, warn};
use rust_htslib::bam::record::Record;

use super::carrier::{bases_in, Guard, ReadRules, GUARD_REF_SHARE};
use super::utils::ClassifyResult;
use crate::normalize::fasta::CachedFasta;
use crate::types::Variant;

pub(crate) use decide::classify;
pub(crate) use haplotypes::{Candidates, ReadView, Row, Site};

/// How far a read reaches each side of a row's REF span: rows within it can share a read
/// with the row (2x250 and 2x300 runs included), so they shape its haplotypes.
pub(crate) const READ_REACH: i64 = super::carrier::WHOLE_READ_REACH;

/// Which read judgment counts a row's reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ReadJudge {
    /// The size-specific rules (RJ-1 to RJ-27): the default.
    Rules,
    /// The sequence judge (`sequence`; `haplotype`, its former name).
    Haplotype,
}

/// Parse the `read_judge` token strictly, as `alignment_backend` is: a typo must fail
/// loudly rather than count with the other judgment.
pub(crate) fn parse_read_judge(token: &str) -> Result<ReadJudge, String> {
    match token {
        "rules" => Ok(ReadJudge::Rules),
        "sequence" | "haplotype" => Ok(ReadJudge::Haplotype),
        other => Err(format!("unknown read_judge {other:?}; expected \"rules\" or \"sequence\"")),
    }
}

/// The input's rows within reach of a substitution `v` that hold its change at its site
/// (`Variant::reach_rows`): rows of the same length over `v`'s bases reading its ALT bases
/// where it changes them, such as an MNV listed beside its component SNV. Their carriers
/// are `v`'s ALT reads (decided 2026-10-25); the diagnostic MNV_SHARED says how many of
/// `v`'s ALT reads spell each whole. Empty for an indel or delins `v`.
pub(crate) fn listed_holders(v: &Variant) -> Vec<&(i64, String, String)> {
    let (r0, a0) = (v.ref_allele.as_bytes(), v.alt_allele.as_bytes());
    if r0.len() != a0.len() {
        return Vec::new();
    }
    let end = v.pos + r0.len() as i64;
    v.reach_rows
        .iter()
        .filter(|(p, r, a)| {
            let (r, a) = (r.as_bytes(), a.as_bytes());
            r.len() == a.len()
                && *p <= v.pos
                && *p + r.len() as i64 >= end
                && (r, a) != (r0, a0)
                && (0..r0.len()).filter(|&k| !r0[k].eq_ignore_ascii_case(&a0[k])).all(|k| {
                    let i = (v.pos - *p) as usize + k;
                    a[i].eq_ignore_ascii_case(&a0[k])
                })
        })
        .collect()
}

/// A record's identity at one row: its name, flags, position, and CIGAR and bases.
pub(crate) type ReadKey = (Vec<u8>, u16, i64, Vec<u8>);

/// A row's haplotypes, built on first use from its reads (the main counts, the
/// per-transcript counts and ASJD read the same ones, whichever asks first), and the
/// call each read got, which they share.
pub(crate) struct Haplotypes<'a> {
    key: (i64, String, String),
    site: OnceCell<Option<Site>>,
    #[allow(clippy::type_complexity)]
    build: Box<dyn Fn(&ReadRules) -> Option<Site> + 'a>,
    verdicts: RefCell<HashMap<ReadKey, ClassifyResult>>,
    /// Each read view this row's site judged (its name cleared), whether it was asked for
    /// the anchor, and `--min-baseq`, with its verdict: the judge is a function of these
    /// alone, so a read alike to one judged gets its verdict (deep duplex reads are often
    /// alike: one start, one end, binned qualities).
    judged: RefCell<HashMap<JudgedKey, decide::Verdict, std::hash::BuildHasherDefault<ViewHasher>>>,
    /// Each read's bases (masked) and placement this row's site scored (a view without its
    /// qualities, tags and name), with what its programmes gave: a read alike there has
    /// the same scores, paths and alignments, its call its own.
    scored: RefCell<HashMap<ReadView, std::rc::Rc<decide::Scored>, std::hash::BuildHasherDefault<ViewHasher>>>,
}

/// The judged views' hasher: keys are compared whole on a hit, so the hash need only
/// spread them; a multiply-rotate over each word, cheaper than the default's per read.
#[derive(Default)]
pub(crate) struct ViewHasher(u64);

impl ViewHasher {
    #[inline]
    fn add(&mut self, w: u64) {
        self.0 = (self.0.rotate_left(5) ^ w).wrapping_mul(0x517c_c1b7_2722_0a95);
    }
}

impl std::hash::Hasher for ViewHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        let (chunks, rest) = bytes.as_chunks::<8>();
        for c in chunks {
            self.add(u64::from_le_bytes(*c));
        }
        if !rest.is_empty() {
            let mut w = [0u8; 8];
            w[..rest.len()].copy_from_slice(rest);
            self.add(u64::from_le_bytes(w) ^ (rest.len() as u64) << 56);
        }
    }
    fn write_u8(&mut self, i: u8) {
        self.add(u64::from(i));
    }
    fn write_u32(&mut self, i: u32) {
        self.add(u64::from(i));
    }
    fn write_u64(&mut self, i: u64) {
        self.add(i);
    }
    fn write_usize(&mut self, i: usize) {
        self.add(i as u64);
    }
    fn write_i64(&mut self, i: i64) {
        self.add(i as u64);
    }
    fn write_isize(&mut self, i: isize) {
        self.add(i as u64);
    }
}

/// A judged read view, whether its anchor was asked for, and `--min-baseq`.
pub(crate) type JudgedKey = (ReadView, bool, u8);

impl<'a> Haplotypes<'a> {
    pub(crate) fn new(v: &Variant, build: impl Fn(&ReadRules) -> Option<Site> + 'a) -> Self {
        Self {
            key: (v.pos, v.ref_allele.clone(), v.alt_allele.clone()),
            site: OnceCell::new(),
            build: Box::new(build),
            verdicts: RefCell::new(HashMap::new()),
            judged: RefCell::new(HashMap::default()),
            scored: RefCell::new(HashMap::default()),
        }
    }

    /// The verdict a read alike to `view`'s got at this row ([`Haplotypes::judged`]'s key).
    pub(crate) fn judged(&self, view: &JudgedKey) -> Option<decide::Verdict> {
        self.judged.borrow().get(view).copied()
    }

    pub(crate) fn keep_judged(&self, view: JudgedKey, verdict: decide::Verdict) {
        self.judged.borrow_mut().insert(view, verdict);
    }

    /// What the programmes gave a read alike to `placed` (bases and placement) at this row.
    pub(crate) fn scored(&self, placed: &ReadView) -> Option<std::rc::Rc<decide::Scored>> {
        self.scored.borrow().get(placed).cloned()
    }

    pub(crate) fn keep_scored(&self, placed: ReadView, scored: std::rc::Rc<decide::Scored>) {
        self.scored.borrow_mut().insert(placed, scored);
    }

    /// The call a read already got at this row.
    pub(crate) fn verdict(&self, key: &ReadKey) -> Option<ClassifyResult> {
        self.verdicts.borrow().get(key).copied()
    }

    pub(crate) fn keep_verdict(&self, key: ReadKey, result: ClassifyResult) {
        self.verdicts.borrow_mut().insert(key, result);
    }

    /// The haplotypes when they were found for `v` (the dual-counted decomposed form of
    /// a variant has its own). None when its window cannot be read.
    pub(crate) fn get(&self, v: &Variant, rules: &ReadRules) -> Option<&Site> {
        let (pos, r, a) = &self.key;
        if *pos != v.pos || *r != v.ref_allele || *a != v.alt_allele {
            return None;
        }
        self.site.get_or_init(|| (self.build)(rules)).as_ref()
    }
}

/// The variant's reads as the judge reads them, from `(record, qualities)` pairs.
pub(crate) fn read_views<'r>(
    reads: impl Iterator<Item = (&'r Record, std::borrow::Cow<'r, [u8]>)>,
    min_baseq: u8,
    rules: &ReadRules,
) -> Vec<ReadView> {
    reads.filter_map(|(record, quals)| decide::read_view(record, &quals, min_baseq, rules)).collect()
}

/// `v`'s haplotypes: its window read from `reference`, its rows within reach
/// (`Variant::reach_rows`), and from its guard the other alleles the reads show at the
/// locus (site states) and the background changes recurring there (contexts). None, with
/// a warning, when the reference does not hold the window or the row's REF.
pub(crate) fn site_for(v: &Variant, guard: Option<&Guard>, reference: &CachedFasta, reads: &[ReadView]) -> Option<Site> {
    let target = Row::new(v.pos, v.ref_allele.as_bytes(), v.alt_allele.as_bytes());
    let listed: Vec<Row> = v.reach_rows.iter().map(|(p, r, a)| Row::new(*p, r.as_bytes(), a.as_bytes())).collect();
    let discovered: Vec<Row> = guard
        .map(|g| g.others.iter().chain(&g.ref_others).map(|o| Row::new(o.pos, &o.ref_allele, &o.alt_allele)).collect())
        .unwrap_or_default();
    let ends = || std::iter::once(&target).chain(&listed).chain(&discovered);
    let want_lo = ends().map(|r| r.pos).min().unwrap_or(v.pos) - READ_REACH;
    let lo = want_lo.max(0);
    let hi = ends().map(|r| r.pos + r.ref_al.len() as i64).max().unwrap_or(v.pos + 1) + READ_REACH;
    let Some(window) = reference.window(&v.chrom, lo, hi) else {
        warn!(
            "sequence judge: the reference does not hold {}:{}-{} around {}:{} {}>{}; its reads are depth only",
            v.chrom, lo + 1, hi, v.chrom, v.pos + 1, v.ref_allele, v.alt_allele,
        );
        return None;
    };
    let base = |p: i64| usize::try_from(p - lo).ok().and_then(|i| window.get(i).copied());
    // The guard's recurring changes beside the event explain reads, never credit them:
    // a masked base taking only some bases is a context row per change; one taking any
    // base is any base in every haplotype; a recurring indel (a germline indel beside
    // the event) is a context row.
    let mut background: Vec<Row> = guard
        .map(|g| {
            g.background
                .iter()
                .filter_map(|&(p, set)| base(p).map(|r| (p, r, set)))
                .flat_map(|(p, r, set)| bases_in(set).filter(move |&b| b != r).map(move |b| Row::new(p, &[r], &[b])))
                .collect()
        })
        .unwrap_or_default();
    let masked: Vec<i64> = guard
        .map(|g| g.masked.iter().copied().filter(|p| !g.background.iter().any(|(q, _)| q == p)).collect())
        .unwrap_or_default();
    if let Some(g) = guard {
        background.extend(g.indels.iter().filter_map(|&(p, n)| indel_row(p, n, &base, reads)));
    }
    let mut guarded: Vec<i64> = guard.map(|g| g.masked.clone()).unwrap_or_default();
    guarded.sort_unstable();
    let found = found_indels(&base, reads, &guarded);
    // A contig edge inside the window: past it the reference is not known (a circular
    // contig runs on at its start), so a call resting on no base there is withdrawn (RJ-3).
    let right_edge = (window.len() as i64) < hi - lo;
    let edges = ((want_lo < 0).then_some(lo), right_edge.then_some(lo + window.len() as i64));
    let candidates = Candidates { listed, discovered, background, masked, found, edges };
    let site = Site::build(target, &candidates, (lo, window), reads);
    match &site {
        Some(s) if s.scorer == haplotypes::ScorerKind::Graph => debug!(
            "haplotype judge {}:{} {}>{}: {} state graphs ({} bases){} from {} reads; site {}-{}",
            v.chrom, v.pos + 1, v.ref_allele, v.alt_allele, s.graphs.len(),
            s.graphs.iter().map(|g| g.graph.len()).sum::<usize>(),
            if s.triage.is_some() { ", triage on" } else { "" }, reads.len(), s.site_span().0 + 1, s.site_span().1,
        ),
        Some(s) => debug!(
            "haplotype judge {}:{} {}>{}: {} haplotypes over {} contexts{} from {} reads; site {}-{}",
            v.chrom, v.pos + 1, v.ref_allele, v.alt_allele, s.haps.len(), s.n_contexts,
            if s.cap_exceeded { " (the base contexts alone pass the cap)" } else { "" }, reads.len(),
            s.site_span().0 + 1, s.site_span().1,
        ),
        None => warn!(
            "sequence judge: the reference does not hold the REF of {}:{} {}>{}; its reads are depth only",
            v.chrom, v.pos + 1, v.ref_allele, v.alt_allele,
        ),
    }
    site
}

/// The alleles the reads' own bases spell at the indels their alignments write
/// (`ReadView::spelled_indels`), on [`haplotypes::FOUND_MIN_READS`] reads or more and
/// [`GUARD_REF_SHARE`] of them, each with the reads (indices into `reads`) spelling it:
/// alleles the reads show, whatever their length (a state where one lies at X's site and
/// the reads tell it apart from X and REF, `haplotypes::FoundRule`). X's own carriers,
/// however the aligner wrote them, spell X, or, ending inside it, the part of X they hold.
fn found_indels(base: &dyn Fn(i64) -> Option<u8>, reads: &[ReadView], masked: &[i64]) -> Vec<(Row, Vec<usize>)> {
    let mut tally: std::collections::BTreeMap<Row, Vec<usize>> = std::collections::BTreeMap::new();
    for (i, r) in reads.iter().enumerate() {
        let mut once = std::collections::BTreeSet::new();
        for row in r.spelled_indels(base, masked) {
            if once.insert(row.clone()) {
                tally.entry(row).or_default().push(i);
            }
        }
    }
    let floor = (GUARD_REF_SHARE * reads.len() as f64).max(f64::from(haplotypes::other_option().found_min));
    tally.into_iter().filter(|(_, by)| by.len() as f64 >= floor).collect()
}

/// A recurring indel the reads' own alignments write (`Guard::indels`: the genome
/// position after its anchor, its length, + inserted, - deleted) as a row: a deletion's
/// bases from the reference, an insertion's as the reads writing it there show them most
/// (the first such sequence on a tie).
fn indel_row(p: i64, n: i64, base: &dyn Fn(i64) -> Option<u8>, reads: &[ReadView]) -> Option<Row> {
    let anchor = base(p - 1)?;
    if n < 0 {
        let deleted: Option<Vec<u8>> = (p..p - n).map(base).collect();
        return Some(Row::new(p - 1, &[&[anchor][..], &deleted?].concat(), &[anchor]));
    }
    let mut seen: Vec<(&[u8], usize)> = Vec::new();
    for bases in reads.iter().flat_map(|r| &r.inserts).filter(|(q, b)| *q == p && b.len() as i64 == n).map(|(_, b)| &b[..]) {
        match seen.iter_mut().find(|(b, _)| *b == bases) {
            Some(s) => s.1 += 1,
            None => seen.push((bases, 1)),
        }
    }
    let most = seen.iter().map(|s| s.1).max()?;
    let bases = seen.iter().find(|s| s.1 == most)?.0;
    Some(Row::new(p - 1, &[anchor], &[&[anchor][..], bases].concat()))
}

#[cfg(test)]
mod tests {
    use super::decide::{judge, Call};
    use super::haplotypes::tests::{contig, random_seq, row, site};
    use super::haplotypes::{ReadView, ScorerKind};
    use super::*;
    use rust_htslib::bam::record::Cigar;

    #[test]
    fn the_read_judge_token_is_parsed_strictly() {
        assert_eq!(parse_read_judge("rules").unwrap(), ReadJudge::Rules);
        assert_eq!(parse_read_judge("sequence").unwrap(), ReadJudge::Haplotype);
        assert_eq!(parse_read_judge("haplotype").unwrap(), ReadJudge::Haplotype, "the former name");
        for bad in ["Sequence", "sequences", "Haplotype", "haplotypes", "", "rule"] {
            assert!(parse_read_judge(bad).is_err(), "{bad:?} must be rejected");
        }
    }

    const Q: u8 = 37;
    const MIN_BQ: u8 = 20;

    /// `contig` with `rows` applied, right to left.
    fn haplotype(contig: &[u8], rows: &[&Row]) -> Vec<u8> {
        let mut rows: Vec<&Row> = rows.to_vec();
        rows.sort_by_key(|r| std::cmp::Reverse(r.pos));
        let mut s = contig.to_vec();
        for r in rows {
            s.splice(r.pos as usize..r.pos as usize + r.ref_al.len(), r.alt_al.iter().copied());
        }
        s
    }

    /// A CIGAR string such as "75M6I70M".
    fn cigar(s: &str) -> Vec<Cigar> {
        let mut out = Vec::new();
        let mut n = 0u32;
        for b in s.bytes() {
            if b.is_ascii_digit() {
                n = n * 10 + u32::from(b - b'0');
            } else {
                out.push(match b {
                    b'M' => Cigar::Match(n),
                    b'I' => Cigar::Ins(n),
                    b'D' => Cigar::Del(n),
                    b'S' => Cigar::SoftClip(n),
                    _ => panic!("op {b}"),
                });
                n = 0;
            }
        }
        out
    }

    fn read_q(seq: &[u8], quals: Vec<u8>, ref_start: i64, cig: &str) -> ReadView {
        ReadView::new(seq, &quals, ref_start, &cigar(cig), (0, seq.len()), MIN_BQ).unwrap()
    }

    /// `len` bases of `hap` from its offset `from`, aligned at `ref_start` with `cig`.
    fn read(hap: &[u8], from: usize, len: usize, ref_start: i64, cig: &str) -> ReadView {
        read_q(&hap[from..from + len], vec![Q; len], ref_start, cig)
    }

    fn calls(contig: &[u8], rows: &[Row], target: &Row, reads: &[ReadView]) -> Vec<Call> {
        let l = reads.iter().map(ReadView::len).max().unwrap();
        let s = site(contig, target, rows, l);
        reads.iter().map(|r| judge(&s, r, MIN_BQ, false).call).collect()
    }

    /// An insertion A>ACGTTAG anchored at 0-based 299, and an SNV A>C at 300 whose ALT
    /// is the first inserted base: the shape an insertion carrier written through (no I
    /// in its CIGAR) takes.
    fn case_246() -> (Vec<u8>, Row, Row) {
        (contig(11, 299, b"AA"), row(300, "A", "ACGTTAG"), row(301, "A", "C"))
    }

    #[test]
    fn an_insertion_carrier_in_trans_is_alt_for_the_insertion_only() {
        use Call::*;
        let (c, ins, snv) = case_246();
        let rows = vec![ins.clone(), snv.clone()];
        let (h_ins, h_snv) = (haplotype(&c, &[&ins]), haplotype(&c, &[&snv]));
        let reads = vec![
            read(&h_ins, 220, 151, 220, "151M"), // written through
            read(&h_ins, 230, 151, 230, "151M"),
            read(&h_ins, 225, 151, 225, "75M6I70M"),
            read(&h_snv, 221, 151, 221, "151M"),
            read(&h_snv, 231, 151, 231, "151M"),
            read(&c, 222, 151, 222, "151M"),
            read(&c, 232, 151, 232, "151M"),
        ];
        assert_eq!(calls(&c, &rows, &ins, &reads), vec![Alt, Alt, Alt, Ref, Ref, Ref, Ref], "insertion row");
        assert_eq!(calls(&c, &rows, &snv, &reads), vec![Ref, Ref, Ref, Alt, Alt, Ref, Ref], "SNV row");
    }

    #[test]
    fn a_reads_name_never_changes_its_verdict() {
        // Reads alike but for their names (duplex reads of one start and end) share one
        // verdict per row (Haplotypes::judged): the judge reads no name.
        let (c, ins, snv) = case_246();
        let rows = vec![ins.clone(), snv.clone()];
        let (h_ins, h_snv) = (haplotype(&c, &[&ins]), haplotype(&c, &[&snv]));
        let mut quals = vec![Q; 151];
        quals[80] = 2;
        let reads = [
            read(&h_ins, 220, 151, 220, "151M"),
            read(&h_ins, 225, 151, 225, "75M6I70M"),
            read(&h_snv, 221, 151, 221, "151M"),
            read(&c, 232, 151, 232, "151M"),
            read_q(&h_snv[231..382], quals, 231, "151M"),
            read(&h_ins, 290, 60, 290, "10S50M"),
        ];
        for target in [&ins, &snv] {
            let s = site(&c, target, &rows, 151);
            for (k, r) in reads.iter().enumerate() {
                for ask in [false, true] {
                    let (a, b) = (judge(&s, &r.clone().named(b"one"), MIN_BQ, ask), judge(&s, &r.clone().named(b"two"), MIN_BQ, ask));
                    assert_eq!(format!("{a:?}"), format!("{b:?}"), "read {k}, anchor asked {ask}");
                }
            }
        }
    }

    /// Each read judged afresh, against the same read named otherwise and the same bases and
    /// placement with other qualities (none below the floor) judged on the first's kept
    /// scoring: the fast paths' verdicts are the slow path's.
    fn agrees_with_afresh(s: &Site, reads: &[ReadView], what: &str) {
        for (k, r) in reads.iter().enumerate() {
            for ask in [false, true] {
                let afresh = judge(s, r, MIN_BQ, ask);
                let renamed = judge(s, &r.clone().named(b"another"), MIN_BQ, ask);
                assert_eq!(format!("{afresh:?}"), format!("{renamed:?}"), "{what}: read {k}, anchor asked {ask}, renamed");
                let (_, kept) = super::decide::judge_kept(s, r, MIN_BQ, ask, None);
                let kept = kept.expect("the graph scorer keeps its scoring");
                let mut other = r.clone();
                for (i, q) in other.quals.iter_mut().enumerate() {
                    *q = if i % 3 == 0 { 25 } else { 41 };
                }
                assert_eq!(other.placed_alone(), r.placed_alone(), "{what}: read {k} alike");
                let (on_kept, _) = super::decide::judge_kept(s, &other, MIN_BQ, ask, Some(kept));
                let other_afresh = judge(s, &other, MIN_BQ, ask);
                assert_eq!(format!("{other_afresh:?}"), format!("{on_kept:?}"), "{what}: read {k}, anchor asked {ask}, kept scoring");
            }
        }
    }

    #[test]
    fn the_fast_paths_agree_with_judging_afresh_on_the_fix_rounds_sites() {
        use Call::*;
        let len = 120;
        // A change at a base an MNP keeps, on both haplotypes (a separate event): CAG>TAC
        // keeping its A, a germline A>T there on half the REF reads and half the carriers.
        let c = contig(21, 300, b"CAG");
        let x = Row::new(300, b"CAG", b"TAC");
        let h = haplotype(&c, &[&x]);
        let germ = |hap: &[u8], s: usize, on: bool| {
            let mut b = hap[s..s + len].to_vec();
            if on {
                b[301 - s] = b'T';
            }
            read_q(&b, vec![Q; len], s as i64, "120M")
        };
        let mut reads: Vec<ReadView> = (0..20).map(|i| germ(&c, 230 + i, i % 2 == 1)).collect();
        reads.extend((0..10).map(|i| germ(&h, 235 + i, i % 2 == 1)));
        let s = Site::build(x.clone(), &Candidates { listed: vec![Row::new(301, b"A", b"T")], ..Candidates::default() }, (0, c.clone()), &reads).unwrap();
        let got: Vec<Call> = reads.iter().map(|r| judge(&s, r, MIN_BQ, false).call).collect();
        assert!(got[..20].iter().all(|&g| g == Ref) && got[20..].iter().all(|&g| g == Alt), "kept-base change: {got:?}");
        agrees_with_afresh(&s, &reads, "kept-base change");

        // An indel the reads write beside the site, a context where they need it: an SNV at
        // 300 and a G inserted after 339 on ten REF reads.
        let c = random_seq(600, 22);
        let alt = if c[300] == b'A' { b'C' } else { b'A' };
        let x = Row::new(300, &c[300..301], &[alt]);
        let ins = Row::new(339, &c[339..340], &[c[339], b'G']);
        let (h, hi) = (haplotype(&c, &[&x]), haplotype(&c, &[&ins]));
        let mut reads: Vec<ReadView> = (0..10).map(|i| read(&c, 240 + i, len, 240 + i as i64, "120M")).collect();
        reads.extend((0..10).map(|i| {
            let s = 250 + i;
            read(&hi, s, len, s as i64, &format!("{}M1I{}M", 340 - s, len - (340 - s) - 1))
        }));
        reads.extend((0..10).map(|i| read(&h, 245 + i, len, 245 + i as i64, "120M")));
        let found = vec![(ins.clone(), (10..20).collect())];
        let s = Site::build(x.clone(), &Candidates { found, ..Candidates::default() }, (0, c.clone()), &reads).unwrap();
        let got: Vec<Call> = reads.iter().map(|r| judge(&s, r, MIN_BQ, false).call).collect();
        assert!(got[..20].iter().all(|&g| g == Ref) && got[20..].iter().all(|&g| g == Alt), "beside-site indel: {got:?}");
        agrees_with_afresh(&s, &reads, "beside-site indel");

        // The REF-base option at a delins (CGC>T): two carriers showing REF's C where the ALT
        // has its T are another allele under `partial`.
        let mut c = random_seq(600, 23);
        c[299] = b'A';
        c[300..303].copy_from_slice(b"CGC");
        c[303] = b'A';
        let x = Row::new(300, b"CGC", b"T");
        let h = haplotype(&c, &[&x]);
        let carrier = |s: usize, odd: Option<u8>| {
            let mut b = h[s..s + len].to_vec();
            if let Some(o) = odd {
                b[300 - s] = o;
            }
            read_q(&b, vec![Q; len], s as i64, &format!("{}M1I3D{}M", 300 - s, len - (300 - s) - 1))
        };
        let mut reads: Vec<ReadView> = (0..20).map(|i| read(&c, 240 + i, len, 240 + i as i64, "120M")).collect();
        reads.extend((0..10).map(|i| carrier(245 + i, None)));
        reads.extend((0..2).map(|i| carrier(250 + i, Some(b'C'))));
        let mut s = Site::build(x.clone(), &Candidates::default(), (0, c.clone()), &reads).unwrap();
        s.ref_base = haplotypes::RefBaseRule::Partial;
        let got: Vec<Call> = reads.iter().map(|r| judge(&s, r, MIN_BQ, false).call).collect();
        assert!(got[30..].iter().all(|&g| g == Other), "REF-base option: {got:?}");
        agrees_with_afresh(&s, &reads, "REF-base option");
    }

    #[test]
    fn a_read_alike_in_bases_and_placement_is_called_on_the_kept_scoring_as_afresh() {
        // Reads alike in their bases (masked) and placement share their programmes' scores,
        // paths and alignments (decide::Scored); qualities above the floor, the aligner's
        // tags and the anchor asked are each read's own, and its call is as if judged afresh.
        let (c, ins, snv) = case_246();
        let rows = vec![ins.clone(), snv.clone()];
        let (h_ins, h_snv) = (haplotype(&c, &[&ins]), haplotype(&c, &[&snv]));
        let reads = [
            read(&h_ins, 220, 151, 220, "151M"),
            read(&h_ins, 225, 151, 225, "75M6I70M"),
            read(&h_snv, 221, 151, 221, "151M"),
            read(&c, 232, 151, 232, "151M"),
            read(&h_ins, 290, 60, 290, "10S50M"),
        ];
        let tags = [
            haplotypes::AlnTags { mapq: 60, best: Some(151), elsewhere: Some(40), nm: Some(0) },
            haplotypes::AlnTags { mapq: 0, best: Some(140), elsewhere: Some(140), nm: Some(2) },
        ];
        for target in [&ins, &snv] {
            let s = site(&c, target, &rows, 151);
            for (k, r) in reads.iter().enumerate() {
                for ask in [false, true] {
                    let (_, kept) = decide::judge_kept(&s, r, MIN_BQ, ask, None);
                    let kept = kept.expect("the graph scorer keeps its scoring");
                    for (t, &tag) in tags.iter().enumerate() {
                        // the same bases and placement, other qualities (none below the floor)
                        let mut other = r.clone().tagged(tag);
                        for (i, q) in other.quals.iter_mut().enumerate() {
                            *q = if i % 3 == 0 { 25 } else { 41 };
                        }
                        assert_eq!(other.placed_alone(), r.placed_alone(), "read {k}: alike");
                        let afresh = judge(&s, &other, MIN_BQ, ask);
                        let (on_kept, _) = decide::judge_kept(&s, &other, MIN_BQ, ask, Some(kept.clone()));
                        assert_eq!(format!("{afresh:?}"), format!("{on_kept:?}"), "read {k}, tags {t}, anchor asked {ask}");
                    }
                }
            }
        }
    }

    #[test]
    fn a_carrier_of_both_in_cis_is_alt_at_both_rows() {
        use Call::*;
        let (c, ins, snv) = case_246();
        let rows = vec![ins.clone(), snv.clone()];
        let both = haplotype(&c, &[&ins, &snv]);
        let reads = vec![
            read(&both, 220, 151, 220, "151M"),
            read(&both, 230, 151, 230, "70M6I75M"),
            read(&c, 222, 151, 222, "151M"),
            read(&c, 232, 151, 232, "151M"),
        ];
        assert_eq!(calls(&c, &rows, &ins, &reads), vec![Alt, Alt, Ref, Ref], "insertion row");
        assert_eq!(calls(&c, &rows, &snv, &reads), vec![Alt, Alt, Ref, Ref], "SNV row");
    }

    #[test]
    fn a_compound_het_counts_each_rows_own_carriers_and_the_others_as_other() {
        use Call::*;
        let c = contig(23, 299, b"AC");
        let (ins1, ins2) = (row(300, "A", "AGT"), row(300, "A", "ATTA"));
        let rows = vec![ins1.clone(), ins2.clone()];
        let (h1, h2) = (haplotype(&c, &[&ins1]), haplotype(&c, &[&ins2]));
        let reads = vec![
            read(&h1, 220, 151, 220, "80M2I69M"),
            read(&h1, 230, 151, 230, "151M"),
            read(&h2, 221, 151, 221, "79M3I69M"),
            read(&h2, 231, 151, 231, "151M"),
            read(&c, 225, 151, 225, "151M"),
        ];
        assert_eq!(calls(&c, &rows, &ins1, &reads), vec![Alt, Alt, Other, Other, Ref], "ins1");
        assert_eq!(calls(&c, &rows, &ins2, &reads), vec![Other, Other, Alt, Alt, Ref], "ins2");
    }

    #[test]
    fn a_read_one_repeat_unit_off_is_partial_evidence() {
        // Another length of the repeat is another allele: partial evidence, never REF or ALT.
        // G anchor at 299, (CA)x8 at 300..316, then a T.
        let c = contig(37, 297, b"TTGCACACACACACACACAT");
        let ins = row(300, "G", "GCA");
        let two_units = row(300, "G", "GCACA"); // not listed: the read's own allele
        let (h_off, h_alt) = (haplotype(&c, &[&two_units]), haplotype(&c, &[&ins]));
        let reads = vec![
            read(&h_off, 220, 151, 220, "80M4I67M"),
            read(&h_alt, 222, 151, 222, "78M2I71M"),
            read(&c, 224, 151, 224, "151M"),
        ];
        assert_eq!(calls(&c, std::slice::from_ref(&ins), &ins, &reads), vec![Call::Other, Call::Alt, Call::Ref]);
    }

    #[test]
    fn a_read_ending_before_the_telling_base_is_depth_only() {
        // A anchor at 299, T x6 at 300..306, then G (the telling base) at 306.
        let c = contig(41, 298, b"CATTTTTTG");
        let ins = row(300, "A", "AT");
        let h_alt = haplotype(&c, &[&ins]);
        let reads = vec![
            read(&c, 152, 151, 152, "151M"), // last base at 302, inside the run
            read(&h_alt, 153, 151, 153, "151M"),
            read(&c, 230, 151, 230, "151M"),
        ];
        assert_eq!(calls(&c, std::slice::from_ref(&ins), &ins, &reads), vec![Call::Depth, Call::Depth, Call::Ref]);
    }

    #[test]
    fn a_masked_telling_base_is_depth_only() {
        use Call::*;
        let c = contig(53, 299, b"GAT");
        let snv = row(301, "A", "C");
        let h_alt = haplotype(&c, &[&snv]);
        let mut with_n = h_alt[225..376].to_vec();
        with_n[300 - 225] = b'N';
        let q2: Vec<u8> = (0..151).map(|i| if i == 300 - 226 { 2 } else { Q }).collect();
        let reads = vec![
            read_q(&with_n, vec![Q; 151], 225, "151M"),
            read_q(&h_alt[226..377], q2, 226, "151M"),
            read(&h_alt, 227, 151, 227, "151M"),
            read(&c, 228, 151, 228, "151M"),
        ];
        assert_eq!(calls(&c, std::slice::from_ref(&snv), &snv, &reads), vec![Depth, Depth, Alt, Ref]);
    }

    #[test]
    fn a_linked_variant_never_credits_a_read_blind_at_x() {
        use Call::*;
        let c = contig(61, 299, b"GAT");
        let snv = row(301, "A", "C"); // X at 0-based 300
        let linked = row(251, &(c[250] as char).to_string(), if c[250] == b'G' { "T" } else { "G" });
        let rows = vec![linked.clone(), snv.clone()];
        let h = haplotype(&c, &[&linked, &snv]);
        let mut blind = h[200..351].to_vec();
        blind[300 - 200] = b'N';
        let reads = vec![
            read_q(&blind, vec![Q; 151], 200, "151M"),
            read(&h, 210, 151, 210, "151M"),
            read(&c, 220, 151, 220, "151M"),
        ];
        assert_eq!(calls(&c, &rows, &snv, &reads), vec![Depth, Alt, Ref]);
        // The linked row itself is read by all three.
        assert_eq!(calls(&c, &rows, &linked, &reads), vec![Alt, Alt, Ref]);
    }

    #[test]
    fn an_allele_the_reads_spell_with_a_listed_neighbour_is_that_neighbours() {
        use Call::*;
        // The guard's other alleles at the insertion, spelled over its window: the
        // insertion with the SNV beside it (the cis carriers), and REF with the SNV (the
        // SNV's carriers). Each is a listed row's haplotype, never another allele.
        let (c, ins, snv) = case_246();
        let (lo, hi) = (296usize, 304usize);
        let window = |h: &[u8], extra: usize| Row::new(lo as i64, &c[lo..hi], &h[lo..hi + extra]);
        let (both, h_snv) = (haplotype(&c, &[&ins, &snv]), haplotype(&c, &[&snv]));
        let candidates = Candidates {
            listed: vec![snv.clone()],
            discovered: vec![window(&both, 6), window(&h_snv, 0)],
            ..Candidates::default()
        };
        let reads = [read(&both, 220, 151, 220, "151M"), read(&h_snv, 221, 151, 221, "151M"), read(&c, 222, 151, 222, "151M")];
        let s = Site::build_with(ScorerKind::Enumerated, ins.clone(), &candidates, (0, c.clone()), &reads).unwrap();
        let got: Vec<Call> = reads.iter().map(|r| judge(&s, r, MIN_BQ, false).call).collect();
        assert_eq!(got, vec![Alt, Ref, Ref]);
        assert_eq!(s.haps.len(), 4, "REF, X, and each with the SNV: the spelled alleles are theirs");
    }

    #[test]
    fn an_snv_in_a_repeat_is_read_beside_a_listed_expansion_on_either_side_of_it() {
        // G anchor at 299, (CA)x8 at 300..316, then T. X: unit 5's A. A listed one-unit
        // expansion (or contraction) sits left-aligned at the anchor; a molecule may carry
        // it after X's unit instead, the same row on the reference, another sequence with X.
        let c = contig(37, 299, b"GCACACACACACACACAT");
        let x = row(310, "A", "G");
        for (ctx, right_of_x) in
            [(row(300, "G", "GCA"), Row::new(315, b"A", b"ACA")), (row(300, "GCA", "G"), Row::new(313, b"ACA", b"A"))]
        {
            let s = site(&c, &x, std::slice::from_ref(&ctx), 151);
            for placed in [&ctx, &right_of_x] {
                let h = haplotype(&c, &[placed, &x]);
                let r = read(&h, 230, 151, 230, "151M");
                assert_eq!(judge(&s, &r, MIN_BQ, false).call, Call::Alt, "X beside {placed:?}: its bases say G at X");
            }
        }
    }

    #[test]
    #[ignore = "open: a junk clip next to X is absorbed through an insertion (a junk base \
                aligned, the rest inserted at one a base, X's base with it), which no end-gap \
                charge reaches; see the step 3 fix-round report"]
    fn a_junk_end_never_absorbs_x() {
        // A 40-base junk clip before the bases carrying X, two bases into them: the junk
        // should cost as much as mismatches would, so it cannot take X's base along.
        let c = random_seq(600, 3);
        let x = Row::new(301, &c[301..302], if c[301] == b'A' { b"C" } else { b"A" });
        let s = site(&c, &x, &[], 140);
        let h = haplotype(&c, &[&x]);
        let seq = [&random_seq(40, 77)[..], &h[299..399]].concat();
        let r = read_q(&seq, vec![Q; 140], 299, "40S100M");
        assert_eq!(judge(&s, &r, MIN_BQ, false).call, Call::Alt);
    }

    #[test]
    fn an_mnp_alt_call_needs_every_changed_base_read() {
        use Call::*;
        let c = random_seq(600, 5);
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let alt: Vec<u8> = c[300..303].iter().map(|&b| other(b)).collect();
        let mnp = Row::new(300, &c[300..303], &alt);
        let h = haplotype(&c, &[&mnp]);
        let reads = vec![
            read(&h, 250, 100, 250, "100M"), // the whole MNP
            read(&h, 201, 100, 201, "100M"), // ends on its first base
            read(&c, 201, 100, 201, "100M"), // REF, ends on its first base
            read(&c, 250, 100, 250, "100M"),
        ];
        assert_eq!(calls(&c, std::slice::from_ref(&mnp), &mnp, &reads), vec![Alt, Depth, Depth, Ref]);
    }

    #[test]
    fn many_rows_in_reach_never_take_more_than_the_cap() {
        let c = random_seq(900, 21);
        let other = |b: u8| if b == b'A' { "C" } else { "A" };
        let rows: Vec<Row> = (300..600).step_by(3).map(|p| Row::new(p, &c[p as usize..p as usize + 1], other(c[p as usize]).as_bytes())).collect();
        let x = rows[50].clone();
        let reads: Vec<ReadView> = (0..30).map(|i| read(&c, 380 + i, 150, 380 + i as i64, "150M")).collect();
        let candidates = Candidates { listed: rows.into_iter().filter(|r| *r != x).collect(), ..Candidates::default() };
        let s = Site::build_with(ScorerKind::Enumerated, x, &candidates, (0, c.clone()), &reads).unwrap();
        assert!(s.haps.len() <= haplotypes::MAX_HAPS, "{} haplotypes", s.haps.len());
    }

    #[test]
    fn a_call_rests_on_the_quality_of_the_bases_that_decided_it() {
        let c = random_seq(600, 5);
        let alt = if c[299] == b'A' { "C" } else { "A" };
        let snv = row(300, &(c[299] as char).to_string(), alt);
        let h_alt = haplotype(&c, &[&snv]);
        let mut quals = vec![Q; 151];
        quals[299 - 220] = 25;
        quals[10] = 21; // a base both haplotypes read alike
        let carrier = read_q(&h_alt[220..371], quals, 220, "151M");
        let s = site(&c, &snv, &[], 151);
        let v = judge(&s, &carrier, MIN_BQ, false);
        assert_eq!((v.call, v.qual), (Call::Alt, 25));
    }

    /// A 113-base deletion after the anchor at 0-based 250 on a random contig.
    fn long_deletion() -> (Vec<u8>, Row) {
        let c = random_seq(700, 83);
        let del = Row::new(250, &c[250..364], &c[250..251]);
        (c, del)
    }

    #[test]
    fn a_substitution_inside_a_deletions_bases_is_a_context_on_refs_side() {
        use Call::*;
        // A listed SNV 86 bases into the deleted bases: REF's bases can carry it, X's
        // never hold them. Reads over the anchor short of it are REF, whatever it is.
        let (c, del) = long_deletion();
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let snv = Row::new(336, &c[336..337], &[other(c[336])]);
        let rows = vec![del.clone(), snv.clone()];
        let (h_del, h_snv) = (haplotype(&c, &[&del]), haplotype(&c, &[&snv]));
        let reads = vec![
            read(&c, 180, 100, 180, "100M"),     // REF over the anchor, short of the SNV
            read(&c, 200, 100, 200, "100M"),
            read(&h_snv, 240, 100, 240, "100M"), // REF's bases with the SNV, over the anchor
            read(&h_del, 180, 100, 180, "71M113D29M"),
        ];
        assert_eq!(calls(&c, &rows, &del, &reads), vec![Ref, Ref, Ref, Alt], "the deletion's row");
        assert_eq!(calls(&c, &rows, &snv, &reads[2..]), vec![Alt, Other], "the SNV's row: the deletion is another allele there");
    }

    #[test]
    fn a_substitution_in_a_deletions_window_stays_another_allele() {
        use Call::*;
        // A 20-base deletion's bases are all its window (RJ-5), and a long deletion's first
        // bases its junction's (RJ-24): a listed SNV there is another allele, so a read over
        // the anchor carrying it is partial evidence; one short of it shows REF.
        let c = random_seq(700, 89);
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let short = Row::new(250, &c[250..271], &c[250..251]);
        let (_, long) = long_deletion();
        let c_long = long_deletion().0;
        for (c, del, at, cig) in [(&c, &short, 262usize, "71M20D29M"), (&c_long, &long, 253, "71M113D29M")] {
            let snv = Row::new(at as i64, &c[at..=at], &[other(c[at])]);
            let rows = vec![del.clone(), snv.clone()];
            let (h_del, h_snv) = (haplotype(c, &[del]), haplotype(c, &[&snv]));
            let reads = vec![
                read(c, at - 100, 100, (at - 100) as i64, "100M"), // over the anchor, ends before the SNV
                read(&h_snv, 200, 100, 200, "100M"),
                read(&h_del, 180, 100, 180, cig),
            ];
            assert_eq!(calls(c, &rows, del, &reads), vec![Ref, Other, Alt], "SNV at {at}");
        }
    }

    #[test]
    fn a_long_deletion_starting_within_a_reads_reach_is_a_context() {
        use Call::*;
        // X: a 1-base deletion in an A run at 300; a listed 113-base deletion from 16 bases
        // on, its bases running past any read's reach of X. A read over X carrying the long
        // deletion fits REF (or X) with it, and is judged by its bases at X.
        let c = contig(97, 299, b"GAAAAAC");
        let x = Row::new(299, b"GA", b"G");
        let long = Row::new(316, &c[316..430], &c[316..317]);
        let rows = vec![x.clone(), long.clone()];
        let (h_ref, h_alt) = (haplotype(&c, &[&long]), haplotype(&c, &[&x, &long]));
        let reads = vec![
            read(&h_ref, 250, 100, 250, "67M113D33M"),
            read(&h_alt, 250, 100, 250, "50M1D16M113D34M"),
            read(&c, 251, 100, 251, "100M"),
        ];
        assert_eq!(calls(&c, &rows, &x, &reads), vec![Ref, Alt, Ref]);
    }

    #[test]
    fn a_tie_with_another_allele_the_read_never_reaches_is_no_evidence_of_it() {
        use Call::*;
        // A listed 60-base deletion from 80 bases into X's deleted bases to past X's end:
        // another allele at X's site. A REF read over X's anchor ending before it starts
        // reads nothing that tells REF from it: REF, X beaten. A read carrying it is it.
        let (c, del) = long_deletion();
        let sib = Row::new(330, &c[330..391], &c[330..331]);
        let rows = vec![del.clone(), sib.clone()];
        let (h_del, h_sib) = (haplotype(&c, &[&del]), haplotype(&c, &[&sib]));
        let reads = vec![
            read(&c, 180, 100, 180, "100M"),
            read(&c, 220, 100, 220, "100M"),
            read(&h_sib, 240, 100, 240, "91M60D9M"),
            read(&h_del, 180, 100, 180, "71M113D29M"),
        ];
        assert_eq!(calls(&c, &rows, &del, &reads), vec![Ref, Ref, Other, Alt]);
    }

    #[test]
    fn an_mnp_carrier_with_a_changed_base_below_the_quality_floor_is_alt() {
        use Call::*;
        // Every changed base covered, one below --min-baseq: the confident one decides,
        // as an MNP read with a differing base below the floor counts ALT under the rules
        // (RJ-25). A read short of a changed base stays depth only.
        let c = random_seq(600, 5);
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let alt: Vec<u8> = c[300..302].iter().map(|&b| other(b)).collect();
        let mnp = Row::new(300, &c[300..302], &alt);
        let h = haplotype(&c, &[&mnp]);
        let low = |i: usize| -> Vec<u8> { (0..100).map(|k| if k == i { 2 } else { Q }).collect() };
        let reads = [
            read_q(&h[250..350], low(51), 250, "100M"), // the second changed base low
            read_q(&h[250..350], low(50), 250, "100M"), // the first
            read_q(&c[250..350], low(51), 250, "100M"), // REF, the second low
            read_q(&h[250..350], vec![Q; 100], 250, "100M"),
            read(&h, 201, 100, 201, "100M"), // ends on the first changed base
        ];
        let s = site(&c, &mnp, &[], 100);
        let got: Vec<(Call, bool)> = reads.iter().map(|r| judge(&s, r, MIN_BQ, false)).map(|v| (v.call, v.confirmed)).collect();
        assert_eq!(got, vec![(Alt, false), (Alt, false), (Ref, false), (Alt, true), (Depth, false)]);
    }

    #[test]
    fn a_change_another_allele_carries_beside_x_is_a_context_of_every_state() {
        use Call::*;
        // X: an MNP at 300-301. An allele the reads show over 298-305 (the guard's), with
        // X's second base, REF's first and a change at 304 beside X: another allele at X's
        // site. X's carriers carrying the change at 304 too are X's; the allele's own
        // carriers are another allele.
        let c = random_seq(600, 29);
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let x = Row::new(300, &c[300..302], &[other(c[300]), other(c[301])]);
        let mut spelled = c[298..306].to_vec();
        spelled[3] = other(c[301]);
        spelled[6] = other(c[304]);
        let allele = Row::new(298, &c[298..306], &spelled);
        let beside = Row::new(304, &c[304..305], &[other(c[304])]);
        let candidates = Candidates { discovered: vec![allele.clone()], ..Candidates::default() };
        let (h_x, h_allele) = (haplotype(&c, &[&x, &beside]), haplotype(&c, &[&allele]));
        let reads = [read(&h_x, 250, 100, 250, "100M"), read(&h_allele, 251, 100, 251, "100M"), read(&c, 252, 100, 252, "100M")];
        let s = Site::build(x, &candidates, (0, c.clone()), &reads).unwrap();
        let got: Vec<Call> = reads.iter().map(|r| judge(&s, r, MIN_BQ, false).call).collect();
        assert_eq!(got, vec![Alt, Other, Ref]);
    }

    #[test]
    fn an_snvs_carriers_include_a_listed_mnv_holding_its_change() {
        use Call::*;
        // An SNV at 300 and a listed MNV over 300-301 with the SNV's base at 300. At the
        // SNV's row a read carrying the MNV carries its change: ALT, and it spells the MNV
        // whole. At the MNV's row a read carrying the SNV alone is another allele.
        let c = random_seq(600, 31);
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let snv = Row::new(300, &c[300..301], &[other(c[300])]);
        let mnv = Row::new(300, &c[300..302], &[other(c[300]), other(c[301])]);
        let rows = vec![snv.clone(), mnv.clone()];
        let (h_snv, h_mnv) = (haplotype(&c, &[&snv]), haplotype(&c, &[&mnv]));
        let reads = vec![
            read(&h_mnv, 250, 100, 250, "100M"),
            read(&h_snv, 251, 100, 251, "100M"),
            read(&c, 252, 100, 252, "100M"),
            read(&h_mnv, 201, 100, 201, "100M"), // ends on the SNV's base: the SNV or the MNV
        ];
        assert_eq!(calls(&c, &rows, &snv, &reads), vec![Alt, Alt, Ref, Alt], "the SNV's row");
        assert_eq!(calls(&c, &rows, &mnv, &reads), vec![Alt, Other, Ref, Other], "the MNV's row");
        let s = site(&c, &snv, &rows, 100);
        let spells: Vec<u32> = reads.iter().map(|r| judge(&s, r, MIN_BQ, false).spells).collect();
        assert_eq!(spells, vec![1, 0, 0, 0], "only the read showing the MNV's second base spells it");
    }

    /// Calls on `x`'s site with the guard's background kept by `bg`.
    fn calls_backed(contig: &[u8], x: &Row, candidates: &Candidates, reads: &[ReadView], bg: haplotypes::BackgroundRule) -> Vec<Call> {
        let s = Site::build_backed(
            ScorerKind::Graph,
            x.clone(),
            candidates,
            (0, contig.to_vec()),
            reads,
            haplotypes::other_option().rule,
            haplotypes::ContextSources::All,
            bg,
        )
        .unwrap();
        reads.iter().map(|r| judge(&s, r, MIN_BQ, false).call).collect()
    }

    #[test]
    fn a_background_change_respelling_the_rows_allele_is_left_out() {
        use haplotypes::BackgroundRule::{AsIs, Respell};
        use Call::*;
        // chr5:44868216 A>AT (#246 arm): with the guard's recurring T deletion two bases on
        // as a context, X spells the SNV haplotype (REF at X's row) exactly, so its reads
        // tie between X and REF. Left out, they are REF; the insertion's carriers stay ALT.
        let (c, x, listed, deletion) = haplotypes::tests::respelt_insertion();
        let h_snv = haplotype(&c, &[&listed[0], &listed[1], &listed[2]]);
        let h_ins = haplotype(&c, &[&listed[0], &x, &listed[2]]);
        let reads = vec![
            read(&h_snv, 250, 100, 250, "100M"),
            read(&h_snv, 256, 100, 256, "100M"),
            read(&h_ins, 250, 100, 250, "54M1I45M"),
            read(&h_ins, 257, 100, 257, "47M1I52M"),
        ];
        let candidates = Candidates { listed, background: vec![deletion], ..Candidates::default() };
        assert_eq!(calls_backed(&c, &x, &candidates, &reads, AsIs), vec![Depth, Depth, Alt, Alt]);
        assert_eq!(calls_backed(&c, &x, &candidates, &reads, Respell), vec![Ref, Ref, Alt, Alt]);
    }

    #[test]
    fn background_indels_spelling_a_carrier_together_are_left_out() {
        use haplotypes::BackgroundRule::{AsIs, Respell};
        use Call::*;
        // chr16:23024630 TA>AC (HG002): the carriers hold X with a 6 bp deletion beside
        // it, written as the guard's recurring 7 bp deletion and an inserted C. REF with
        // both spells a carrier exactly, so it scores REF. Left out, the carriers are X's
        // (X with a phased deletion beside it); REF's reads stay REF.
        let (c, x, carrier, background) = haplotypes::tests::split_carrier();
        let h = haplotype(&c, &[&carrier]);
        let reads = vec![
            read(&h, 240, 100, 240, "53M6D47M"),
            read(&h, 246, 100, 246, "47M6D53M"),
            read(&c, 241, 100, 241, "100M"),
        ];
        let candidates = Candidates { background: background.to_vec(), ..Candidates::default() };
        assert_eq!(calls_backed(&c, &x, &candidates, &reads, AsIs), vec![Ref, Ref, Ref]);
        assert_eq!(calls_backed(&c, &x, &candidates, &reads, Respell), vec![Alt, Alt, Ref]);
    }

    #[test]
    fn a_mask_on_a_duplications_one_telling_base_stays() {
        use haplotypes::BackgroundRule::{AsIs, Respell};
        use Call::*;
        // A 36 bp duplication differing from the bases it copies at one base, which the
        // guard masks (FLT3 ITD, gate set): a carrier ending 25 bases into it reads as REF
        // with the base masked, so it ties. The refinement keeps the mask: its reads cannot
        // be told from a germline SNP's whose haplotype reads as the insertion over a short
        // read, which never counts ALT (decided).
        let (c, x, telling) = haplotypes::tests::respelt_duplication();
        let h = haplotype(&c, &[&x]);
        let reads = vec![read(&h, 226, 100, 226, "100M"), read(&c, 226, 100, 226, "100M")];
        let candidates = Candidates { masked: vec![telling], ..Candidates::default() };
        assert_eq!(calls_backed(&c, &x, &candidates, &reads, AsIs), vec![Depth, Ref]);
        assert_eq!(calls_backed(&c, &x, &candidates, &reads, Respell), vec![Depth, Ref]);
    }

    #[test]
    fn a_germline_change_the_reads_need_stays_though_it_respells_over_a_read() {
        use haplotypes::BackgroundRule::{AsIs, Respell};
        use Call::*;
        // #278's shape: X deletes one G of GGG before eight Cs; every molecule carries a
        // germline C at the run's last G (the guard's background). With REF it reads as X
        // over a read ending in the Cs, but the REF molecules need it: without it they are
        // depth only. Kept, they are REF.
        let c = contig(91, 296, b"ACAGCGGGCCCCCCCCAGGG");
        let x = Row::new(300, b"CG", b"C");
        let germ = Row::new(303, b"G", b"C");
        let h = haplotype(&c, &[&germ]);
        let reads: Vec<ReadView> = (0..20).map(|i| read(&h, 240 + i, 100, 240 + i as i64, "100M")).collect();
        let candidates = Candidates { background: vec![germ], ..Candidates::default() };
        let want = vec![Ref; 20];
        assert_eq!(calls_backed(&c, &x, &candidates, &reads, AsIs), want);
        assert_eq!(calls_backed(&c, &x, &candidates, &reads, Respell), want);
    }

    #[test]
    fn a_read_carrying_some_of_an_mnvs_changes_but_not_all_is_another_allele() {
        use Call::*;
        // A read whose bases carry some of a substitution's changes and REF's bases at the
        // others is another allele (RJ-28), whichever of REF and X its margin favours: one
        // change of three (REF by a base), two of three (X by a base), one of two (a tie).
        // A changed base below --min-baseq weighs nothing, as before.
        let c = random_seq(600, 5);
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let alt: Vec<u8> = c[300..303].iter().map(|&b| other(b)).collect();
        let carrying = |k: &[usize]| {
            let mut h = c.clone();
            for &i in k {
                h[300 + i] = alt[i];
            }
            h
        };
        let tnp = Row::new(300, &c[300..303], &alt);
        let low = |i: usize| -> Vec<u8> { (0..100).map(|k| if k == i { 2 } else { Q }).collect() };
        let reads = vec![
            read(&carrying(&[2]), 250, 100, 250, "100M"),
            read(&carrying(&[0]), 251, 100, 251, "100M"),
            read(&carrying(&[0, 1]), 252, 100, 252, "100M"),
            read(&carrying(&[0, 1, 2]), 253, 100, 253, "100M"),
            read(&c, 254, 100, 254, "100M"),
            read_q(&carrying(&[0])[250..350], low(50), 250, "100M"), // the change unread: REF
        ];
        assert_eq!(calls(&c, std::slice::from_ref(&tnp), &tnp, &reads), vec![Other, Other, Other, Alt, Ref, Ref], "TNP");
        let dnp = Row::new(300, &c[300..302], &alt[..2]);
        let reads = vec![
            read(&carrying(&[1]), 250, 100, 250, "100M"),
            read(&carrying(&[0]), 251, 100, 251, "100M"),
            read(&carrying(&[0, 1]), 252, 100, 252, "100M"),
            read(&c, 253, 100, 253, "100M"),
        ];
        assert_eq!(calls(&c, std::slice::from_ref(&dnp), &dnp, &reads), vec![Other, Other, Alt, Ref], "DNP");
    }

    #[test]
    #[ignore = "open: X with a listed change on its inserted copy is no candidate, so such a read \
                is REF by a base; as a candidate it cost 23 true-REF reads at HG002's 1 kb+ \
                insertions (truth v2db) and 23 on the gate set (oracle v3)"]
    fn a_listed_snv_inside_a_tandem_duplications_copied_bases_may_lie_on_either_copy() {
        use Call::*;
        // X: a 26 bp tandem duplication after 0-based 299, copying 300..326; a listed SNV 13
        // bases into the copied bases (0-based 313). X's molecules may carry it on either
        // copy: a read ending 20 bases past the anchor with the SNV's base reads as REF with
        // it and as X with it on the inserted copy alike, depth only. Reads past the copy
        // tell them apart.
        let c = random_seq(600, 71);
        let unit = c[300..326].to_vec();
        let x = Row::new(299, &c[299..300], &[&c[299..300], &unit[..]].concat());
        let alt = if c[313] == b'C' { b'G' } else { b'C' };
        let snv = Row::new(313, &c[313..314], &[alt]);
        let mut on_copy = unit.clone();
        on_copy[13] = alt;
        let x_copy = Row::new(299, &c[299..300], &[&c[299..300], &on_copy[..]].concat());
        let (h_snv, h_copy, h_orig) = (haplotype(&c, &[&snv]), haplotype(&c, &[&x_copy]), haplotype(&c, &[&x, &snv]));
        let reads = vec![
            read(&h_snv, 220, 100, 220, "100M"),
            read(&h_copy, 220, 100, 220, "100M"),
            read(&h_snv, 250, 100, 250, "100M"),
            read(&h_copy, 250, 100, 250, "50M26I24M"),
            read(&h_orig, 250, 100, 250, "50M26I24M"),
            read(&c, 250, 100, 250, "100M"),
        ];
        assert_eq!(calls(&c, &[x.clone(), snv], &x, &reads), vec![Depth, Depth, Ref, Alt, Alt, Ref]);
    }

    #[test]
    fn a_guard_indel_with_which_ref_reads_as_x_is_left_out_where_no_read_needs_it() {
        use haplotypes::BackgroundRule::AsIs;
        use Call::*;
        // AGCA>TGCT between a (GCT)3 run and a (GCA)20 run: X's ALT is (GCT)5 (GCA)20. The
        // guard's recurring indels there are two GCT units more in the first run and one GCA
        // unit more before X. REF with the GCT units is (GCT)5 (GCA)22: X's ALT over every
        // read ending inside the GCA run, so as a context it ties X's carriers with REF. No
        // read needs it for its call: it is left out, and X's carriers are X's. The GCA unit
        // spells no candidate and stays.
        let mut motif = b"GCGCCAGTTTGCTGCTGCTGCAGCA".to_vec();
        motif.extend(b"GCA".repeat(20));
        motif.extend(b"AGAGACTAGCCCC");
        let c = contig(53, 260, &motif);
        let x = row(282, "AGCA", "TGCT");
        let (units, unit) = (Row::new(268, b"T", b"TTGCTGC"), Row::new(278, b"T", b"TGCA"));
        let h_x = haplotype(&c, &[&x]);
        let mut reads: Vec<ReadView> = (226..236).map(|s| read(&h_x, s, 100, s as i64, "100M")).collect();
        reads.extend((228..238).map(|s| read(&c, s, 100, s as i64, "100M")));
        let candidates = Candidates { background: vec![units, unit], ..Candidates::default() };
        assert_eq!(calls_backed(&c, &x, &candidates, &reads, AsIs), [vec![Alt; 10], vec![Ref; 10]].concat());
    }

    #[test]
    fn a_read_tying_ref_with_another_allele_it_writes_its_own_indel_across_is_ref() {
        use Call::*;
        // X: a T inserted in a TT run (after 0-based 300). A listed 33 bp deletion from
        // 0-based 261 whose bases repeat their first seven after them: placed right, it
        // ends at X's anchor, so it is another allele at X's site. A read REF at X with an
        // unlisted 11 bp deletion ending 9 bases before X, starting inside the listed deletion's
        // bases, fits REF and that allele about as badly (each wants an indel of its own):
        // REF where X is beaten, never evidence of the deletion. A read carrying the
        // deletion itself is that allele.
        let c = contig(59, 261, b"CAGAGAACAAAAGAAAGAAAGATTTTCTTACTTCAGAGAATTCTTTGCCACGTATTTCTAGC");
        let x = row(301, "A", "AT");
        let del = Row::new(260, &c[260..294], &c[260..261]);
        let rows = vec![x.clone(), del.clone()];
        let own = [&c[275..281], &c[292..386]].concat();
        let h_del = haplotype(&c, &[&del]);
        let reads = vec![
            read_q(&own, vec![Q; 100], 275, "6M11D94M"),
            read(&c, 240, 100, 240, "100M"),
            read(&h_del, 200, 100, 200, "61M33D39M"),
        ];
        assert_eq!(calls(&c, &rows, &x, &reads), vec![Ref, Ref, Other]);
    }

    #[test]
    fn a_covered_telling_base_read_as_neither_allele_is_another_allele_under_the_switch() {
        use haplotypes::NeitherRule::{Margin, Other as Neither};
        use Call::*;
        // CCG>TTA-like: a read with a confident third base at one of X's changed bases and
        // REF's at the others is REF by a base under the margin (the default), and another
        // allele when the switch says such a base makes it one. X's carriers with a third
        // base are depth only under the margin. A delins likewise at REF's telling bases.
        let c = random_seq(600, 5);
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let alt: Vec<u8> = c[300..303].iter().map(|&b| other(b)).collect();
        let third = |b: u8, a: u8| *b"ACGT".iter().find(|&&t| t != b && t != a).unwrap();
        let with = |h: &[u8], i: usize| {
            let mut h = h.to_vec();
            h[300 + i] = third(c[300 + i], alt[i]);
            h
        };
        let calls_under = |x: &Row, reads: &[ReadView], rule| {
            let mut s = site(&c, x, &[], 100);
            s.neither = rule;
            reads.iter().map(|r| judge(&s, r, MIN_BQ, false).call).collect::<Vec<_>>()
        };
        let tnp = Row::new(300, &c[300..303], &alt);
        let h_x = haplotype(&c, &[&tnp]);
        let reads = vec![
            read(&with(&c, 1), 250, 100, 250, "100M"),
            read(&with(&h_x, 1), 251, 100, 251, "100M"),
            read(&c, 252, 100, 252, "100M"),
            read(&h_x, 253, 100, 253, "100M"),
        ];
        assert_eq!(calls_under(&tnp, &reads, Margin), vec![Ref, Depth, Ref, Alt], "the margin");
        assert_eq!(calls_under(&tnp, &reads, Neither), vec![Other, Other, Ref, Alt], "the switch");
        let delins = Row::new(300, &c[300..304], b"GG");
        let h = haplotype(&c, &[&delins]);
        let reads = vec![read(&with(&c, 1), 250, 100, 250, "100M"), read(&c, 251, 100, 251, "100M"), read(&h, 252, 100, 252, "100M")];
        assert_eq!(calls_under(&delins, &reads, Margin), vec![Ref, Ref, Alt], "delins, the margin");
        assert_eq!(calls_under(&delins, &reads, Neither), vec![Other, Ref, Alt], "delins, the switch");
    }

    #[test]
    fn the_neither_switch_below_a_size_leaves_longer_alleles_to_the_margin() {
        use haplotypes::NeitherRule::{Below, Other as Neither};
        use Call::*;
        // A 4-base delins and a 24-base one (their size: the longer allele once the bases
        // both share at either end are set aside), each read by a REF read with a third base
        // at a REF telling base. Below 20 bases the switch makes the short one's read another
        // allele and leaves the long one's REF by the margin; at every size both move.
        let c = random_seq(600, 5);
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let mut reread = c.clone();
        reread[302] = other(c[302]);
        let reads = [read(&reread, 250, 100, 250, "100M")];
        for (x, below_20, every) in [(Row::new(300, &c[300..304], b"GG"), Other, Other), (Row::new(300, &c[300..324], b"GG"), Ref, Other)] {
            let call = |rule| {
                let mut s = site(&c, &x, &[], 100);
                s.neither = rule;
                judge(&s, &reads[0], MIN_BQ, false).call
            };
            assert_eq!(call(Below(20)), below_20, "{} bases, other:20", x.ref_al.len());
            assert_eq!(call(Neither), every, "{} bases, other", x.ref_al.len());
        }
    }

    /// The calls of `reads` at `target` among `rows` on `contig` under an own-gap rule.
    fn calls_own_gap(contig: &[u8], rows: &[Row], target: &Row, reads: &[ReadView], rule: haplotypes::OwnGapRule) -> Vec<Call> {
        let l = reads.iter().map(ReadView::len).max().unwrap();
        let mut s = site(contig, target, rows, l);
        s.own_gap = rule;
        reads.iter().map(|r| judge(&s, r, MIN_BQ, false).call).collect()
    }

    #[test]
    fn under_the_general_own_gap_rule_a_read_beside_another_allele_by_its_own_indel_keeps_the_margins_call() {
        use haplotypes::OwnGapRule::{All, Others, Tie};
        use Call::*;
        // X a T inserted in a TT run (after 0-based 300), beside a listed 33 bp deletion
        // whose right placement ends at X's anchor (another allele at X's site). A read
        // with an unlisted 11 bp deletion inside the listed deletion's bases reaches that
        // allele only through its own gap at its change, and REF and X without one at
        // X's site. REF at X, it ties REF with that allele: REF under either rule (the
        // decided tie rule). Carrying X's T, X beats both: ALT under either. Alike where
        // the read begins two bases before its own gap (a read end), and where the bases
        // beside the gap are below --min-baseq.
        let c = contig(59, 261, b"CAGAGAACAAAAGAAAGAAAGATTTTCTTACTTCAGAGAATTCTTTGCCACGTATTTCTAGC");
        let x = row(301, "A", "AT");
        let del = Row::new(260, &c[260..294], &c[260..261]);
        let rows = vec![x.clone(), del.clone()];
        let hx = haplotype(&c, &[&x]);
        let low = |at: &[usize]| (0..100).map(|i| if at.contains(&i) { 5 } else { Q }).collect::<Vec<u8>>();
        let reads = vec![
            read_q(&[&c[275..281], &c[292..386]].concat(), vec![Q; 100], 275, "6M11D94M"),
            read_q(&[&c[279..281], &c[292..390]].concat(), vec![Q; 100], 279, "2M11D98M"),
            read_q(&[&c[275..281], &c[292..386]].concat(), low(&[4, 5, 6, 7]), 275, "6M11D94M"),
            read_q(&[&hx[275..281], &hx[292..386]].concat(), vec![Q; 100], 275, "6M11D9M1I84M"),
            read_q(&[&hx[279..281], &hx[292..390]].concat(), vec![Q; 100], 279, "2M11D9M1I88M"),
            read_q(&[&hx[275..281], &hx[292..386]].concat(), low(&[4, 5, 6, 7]), 275, "6M11D9M1I84M"),
        ];
        let want = vec![Ref, Ref, Ref, Alt, Alt, Alt];
        assert_eq!(calls_own_gap(&c, &rows, &x, &reads, All), want, "all");
        assert_eq!(calls_own_gap(&c, &rows, &x, &reads, Others), want, "others");
        assert_eq!(calls_own_gap(&c, &rows, &x, &reads, Tie), want, "tie");
    }

    #[test]
    fn under_the_general_own_gap_rule_a_read_writing_its_own_indel_at_x_carries_no_state() {
        use haplotypes::OwnGapRule::{All, Tie};
        use Call::*;
        // Unique sequence: an SNV whose base the read deletes, an MNV inside whose block it
        // inserts a base, a delins beside whose junction it deletes one: every state is
        // reached only through the read's own gap at X's change, so another allele. A read
        // deleting a base four bases from the SNV, outside its change, is a separate event
        // (RJ-8): REF. In a run: X one T inserted, a listed sibling of five; a read of every
        // other length (one to seven) is another allele (RJ-6), X's length ALT; alike for
        // deletions in a G run (X one G, the sibling three). Under either rule.
        let c = random_seq(600, 5);
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let snv = Row::new(300, &c[300..301], &[other(c[300])]);
        let mnv = Row::new(300, &c[300..303], &c[300..303].iter().map(|&b| other(b)).collect::<Vec<_>>());
        let delins = Row::new(300, &c[300..304], b"GG");
        let hd = haplotype(&c, &[&delins]);
        let cases: Vec<(Row, ReadView, Call, &str)> = vec![
            (snv.clone(), read_q(&[&c[250..300], &c[301..351]].concat(), vec![Q; 100], 250, "50M1D50M"), Other, "SNV, its base deleted"),
            (snv.clone(), read_q(&[&c[250..304], &c[305..351]].concat(), vec![Q; 100], 250, "54M1D46M"), Ref, "SNV, a base four away deleted"),
            (mnv.clone(), read_q(&[&c[250..302], b"T".as_slice(), &c[302..349]].concat(), vec![Q; 100], 250, "52M1I47M"), Other, "MNV, a base inserted in its block"),
            (delins.clone(), read_q(&[&hd[250..302], &hd[303..351]].concat(), vec![Q; 100], 250, "52M1D48M"), Other, "delins, a base deleted at its junction"),
        ];
        for rule in [All, Tie] {
            for (x, r, want, what) in &cases {
                assert_eq!(calls_own_gap(&c, std::slice::from_ref(x), x, std::slice::from_ref(r), rule), vec![*want], "{what}, {rule:?}");
            }
        }
        let a = contig(61, 299, b"GCTTA");
        let x = row(301, "C", "CT");
        let sib = row(301, "C", "CTTTTT");
        let reads: Vec<ReadView> = (1..=7usize)
            .map(|n| {
                let h = [&a[..301], "T".repeat(n).as_bytes(), &a[301..]].concat();
                read_q(&h[220..320], vec![Q; 100], 220, &format!("81M{n}I{}M", 19 - n))
            })
            .collect();
        let want = vec![Alt, Other, Other, Other, Other, Other, Other];
        for rule in [All, Tie] {
            assert_eq!(calls_own_gap(&a, &[x.clone(), sib.clone()], &x, &reads, rule), want, "T run, {rule:?}");
        }
        let d = contig(62, 295, b"AGGGGGGGGC");
        let (xd, sd) = (row(296, "AG", "A"), row(296, "AGGG", "A"));
        let reads: Vec<ReadView> = (1..=5usize)
            .map(|n| read_q(&[&d[220..297], &d[297 + n..320 + n]].concat(), vec![Q; 100], 220, &format!("77M{n}D23M")))
            .collect();
        let want = vec![Alt, Other, Other, Other, Other];
        for rule in [All, Tie] {
            assert_eq!(calls_own_gap(&d, &[xd.clone(), sd.clone()], &xd, &reads, rule), want, "G run, {rule:?}");
        }
    }

    #[test]
    fn the_window_switch_moves_a_neither_base_only_where_both_alleles_reads_are_read_alike() {
        use haplotypes::NeitherRule::{Other as Neither, Window};
        use Call::*;
        // An MNV keeping its middle base: a REF read with a third base at the base it keeps
        // (the alleles agree there) stays REF under the window form, one at a base it changes
        // moves. A 24-base delins: a REF read with a third base two bases in (its junction
        // window) moves, one twelve bases in (bases X's molecules never hold) stays REF.
        // Every one moves under `other`.
        let c = random_seq(600, 5);
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let third = |b: u8, a: u8| *b"ACGT".iter().find(|&&t| t != b && t != a).unwrap();
        let reread = |p: usize, b: u8| {
            let mut h = c.clone();
            h[p] = b;
            read(&h, 250, 100, 250, "100M")
        };
        let mnv = Row::new(300, &c[300..303], &[other(c[300]), c[301], other(c[302])]);
        let delins = Row::new(300, &c[300..324], b"GG");
        let cases = [
            (&mnv, reread(301, other(c[301])), Ref, "MNV, a kept base"),
            (&mnv, reread(302, third(c[302], other(c[302]))), Other, "MNV, a changed base"),
            (&delins, reread(302, other(c[302])), Other, "delins, its junction window"),
            (&delins, reread(312, other(c[312])), Ref, "delins, twelve bases in"),
        ];
        for (x, r, under_window, what) in cases {
            let call = |rule| {
                let mut s = site(&c, x, &[], 100);
                s.neither = rule;
                judge(&s, &r, MIN_BQ, false).call
            };
            assert_eq!(call(Window), under_window, "{what}, window");
            assert_eq!(call(Neither), Other, "{what}, other");
        }
    }

    #[test]
    fn the_recurrence_switch_moves_a_neither_base_only_where_other_reads_show_it() {
        use haplotypes::NeitherRule::Recur;
        use Call::*;
        // A TNP: three REF reads read the same third base at its second changed base, one
        // reads another third base at its third (alone). With one other read required, the
        // three are another allele and the lone one stays REF; with three required, all stay
        // REF by the margin.
        let c = random_seq(600, 5);
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let alt: Vec<u8> = c[300..303].iter().map(|&b| other(b)).collect();
        let third = |b: u8, a: u8| *b"ACGT".iter().find(|&&t| t != b && t != a).unwrap();
        let tnp = Row::new(300, &c[300..303], &alt);
        let with = |i: usize| {
            let mut h = c.clone();
            h[300 + i] = third(c[300 + i], alt[i]);
            h
        };
        let mut reads: Vec<ReadView> = (0..3).map(|k| read(&with(1), 250 + k, 100, 250 + k as i64, "100M")).collect();
        reads.push(read(&with(2), 260, 100, 260, "100M"));
        reads.extend((0..4).map(|k| read(&c, 255 + k, 100, 255 + k as i64, "100M")));
        let calls_under = |rule| {
            let mut s = Site::build(tnp.clone(), &Candidates::default(), (0, c.clone()), &reads).unwrap();
            s.neither = rule;
            reads.iter().map(|r| judge(&s, r, MIN_BQ, false).call).collect::<Vec<_>>()
        };
        assert_eq!(calls_under(Recur(1)), [vec![Other; 3], vec![Ref; 5]].concat(), "one other read");
        assert_eq!(calls_under(Recur(3)), vec![Ref; 8], "three other reads");
    }

    #[test]
    fn an_allele_ending_at_an_insertions_point_is_beside_it_under_the_measurement_option() {
        use haplotypes::EndsRule::{Beside, Site as AtSite};
        use Call::*;
        // Item 6's shape: a listed 33 bp deletion whose right placement ends exactly where X
        // inserts its T. At X's site by its placements (the default), its carriers are
        // another allele; beside X (the option), they are REF at X's row.
        let c = contig(59, 261, b"CAGAGAACAAAAGAAAGAAAGATTTTCTTACTTCAGAGAATTCTTTGCCACGTATTTCTAGC");
        let x = row(301, "A", "AT");
        let del = Row::new(260, &c[260..294], &c[260..261]);
        let h_del = haplotype(&c, &[&del]);
        let reads = vec![read(&h_del, 200, 100, 200, "61M33D39M"), read(&c, 240, 100, 240, "100M")];
        let candidates = Candidates { listed: vec![del], ..Candidates::default() };
        let calls_under = |ends| {
            let s = Site::build_ended(x.clone(), &candidates, (0, c.clone()), &reads, ends).unwrap();
            reads.iter().map(|r| judge(&s, r, MIN_BQ, false).call).collect::<Vec<_>>()
        };
        assert_eq!(calls_under(AtSite), vec![Other, Ref], "at X's site");
        assert_eq!(calls_under(Beside), vec![Ref, Ref], "beside X");
    }

    #[test]
    fn an_allele_inserting_at_an_insertions_point_stays_at_its_site_under_the_beside_option() {
        use haplotypes::EndsRule::{Beside, Site as AtSite};
        use Call::*;
        // Another insertion at X's anchor, and X's insertion behind a substituted anchor (a
        // delins): each ends at X's point but inserts there, so it is at X's site under the
        // option too, and its carriers stay another allele (two insertions at one anchor
        // count their own carriers; a carrier behind an anchor substitution is another
        // allele). Only an allele that ends there and inserts nothing is beside X.
        let c = contig(59, 261, b"CAGAGAACAAAAGAAAGAAAGATTTTCTTACTTCAGAGAATTCTTTGCCACGTATTTCTAGC");
        let x = row(301, "A", "AT");
        for sibling in [row(301, "A", "AG"), row(301, "A", "CT")] {
            let h = haplotype(&c, &[&sibling]);
            let reads = vec![read(&h, 240, 100, 240, "61M1I38M"), read(&c, 240, 100, 240, "100M")];
            let candidates = Candidates { listed: vec![sibling.clone()], ..Candidates::default() };
            for ends in [AtSite, Beside] {
                let s = Site::build_ended(x.clone(), &candidates, (0, c.clone()), &reads, ends).unwrap();
                let calls: Vec<Call> = reads.iter().map(|r| judge(&s, r, MIN_BQ, false).call).collect();
                assert_eq!(calls, vec![Other, Ref], "{}>{}, {ends:?}", String::from_utf8_lossy(&sibling.ref_al), String::from_utf8_lossy(&sibling.alt_al));
            }
        }
    }

    #[test]
    fn an_allele_with_a_placement_clear_of_an_indels_bases_is_beside_it_under_the_apart_option() {
        use haplotypes::EndsRule::{Apart, Beside, Site as AtSite};
        use Call::*;
        // The general form of the beside option: the listed 33 bp deletion ending at X's
        // point at its rightmost placement is clear of X's bases at its leftmost, so beside
        // X; another insertion at X's anchor and X's behind a substituted anchor touch X's
        // anchor wherever written, so at X's site; and a delins removing bases that replaces
        // X's anchor (one placement, ending at X's point) is at X's site, where the narrower
        // option takes it as beside.
        let c = contig(59, 261, b"CAGAGAACAAAAGAAAGAAAGATTTTCTTACTTCAGAGAATTCTTTGCCACGTATTTCTAGC");
        let x = row(301, "A", "AT");
        let del = Row::new(260, &c[260..294], &c[260..261]);
        let short = Row::new(297, &c[297..301], b"G");
        let cases: [(Row, &str, [Call; 3]); 4] = [
            (del.clone(), "61M33D39M", [Other, Ref, Ref]),
            (row(301, "A", "AG"), "61M1I38M", [Other, Other, Other]),
            (row(301, "A", "CT"), "61M1I38M", [Other, Other, Other]),
            (short.clone(), "57M3D43M", [Other, Ref, Other]),
        ];
        for (allele, cig, want) in cases {
            let h = haplotype(&c, &[&allele]);
            let from = if allele == del { 200 } else { 240 };
            let reads = vec![read(&h, from, 100, from as i64, cig), read(&c, 240, 100, 240, "100M")];
            let candidates = Candidates { listed: vec![allele.clone()], ..Candidates::default() };
            for (ends, w) in [(AtSite, want[0]), (Beside, want[1]), (Apart, want[2])] {
                let s = Site::build_ended(x.clone(), &candidates, (0, c.clone()), &reads, ends).unwrap();
                let calls: Vec<Call> = reads.iter().map(|r| judge(&s, r, MIN_BQ, false).call).collect();
                let name = format!("{}>{}", allele.ref_al.len(), allele.alt_al.len());
                assert_eq!(calls, vec![w, Ref], "{name}, {ends:?}");
            }
        }
    }

    #[test]
    fn the_neither_switch_leaves_a_base_the_guards_background_explains_to_the_margin() {
        use haplotypes::NeitherRule::{Other as Neither, Recur};
        use Call::*;
        // A germline change in the guard's background at a base an MNV keeps, and one inside
        // a delins's deleted bases: REF reads carrying it read background there, no third
        // allele's base (the background at a base the MNV keeps is masked, decided
        // 2026-10-14), so the switch leaves them REF. A REF read with a base the background
        // does not explain there moves.
        let c = random_seq(600, 5);
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let third = |b: u8, a: u8| *b"ACGT".iter().find(|&&t| t != b && t != a).unwrap();
        let mnv = Row::new(300, &c[300..303], &[other(c[300]), c[301], other(c[302])]);
        let delins = Row::new(300, &c[300..304], b"GG");
        for (x, g) in [(mnv, 301usize), (delins, 302)] {
            let germ = Row::new(g as i64, &c[g..=g], &[other(c[g])]);
            let h = haplotype(&c, &[&germ]);
            let mut reads: Vec<ReadView> = (0..6).map(|k| read(&h, 250 + k, 100, 250 + k as i64, "100M")).collect();
            let mut lone = c.clone();
            lone[g] = third(c[g], other(c[g]));
            reads.push(read(&lone, 260, 100, 260, "100M"));
            let candidates = Candidates { background: vec![germ], ..Candidates::default() };
            for rule in [Neither, Recur(0)] {
                let mut s = Site::build(x.clone(), &candidates, (0, c.clone()), &reads).unwrap();
                s.neither = rule;
                let calls: Vec<Call> = reads.iter().map(|r| judge(&s, r, MIN_BQ, false).call).collect();
                assert_eq!(calls, [vec![Ref; 6], vec![Other]].concat(), "{} bases, {rule:?}", x.ref_al.len());
            }
        }
    }

    #[test]
    fn two_reads_writing_an_indel_find_it_within_the_share_of_the_sites_reads() {
        let c = random_seq(600, 17);
        let base = |p: i64| c.get(p as usize).copied();
        let ins = Row::new(299, &c[299..300], &[&c[299..300], b"GG"].concat());
        let h = haplotype(&c, &[&ins]);
        let writer = |s: usize| read(&h, s, 150, s as i64, &format!("{}M2I{}M", 300 - s, 148 - (300 - s)));
        let plain = |s: usize| read(&c, s, 150, s as i64, "150M");
        let found = |reads: &[ReadView]| found_indels(&base, reads, &[]).into_iter().map(|(r, by)| (r, by.len())).collect::<Vec<_>>();
        // one read writing it is no candidate; two are, unless below the share of the site's
        let one = [writer(220), plain(230)];
        assert!(found(&one).is_empty());
        let two = [writer(220), writer(230), plain(240)];
        let got = found(&two);
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].1, 2);
        let deep: Vec<ReadView> = [writer(220), writer(230)].into_iter().chain((0..60).map(|i| plain(200 + i))).collect();
        assert!(found(&deep).is_empty(), "2 of 62 is below the share");
    }
}
