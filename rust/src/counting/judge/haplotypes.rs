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

use super::score::{self, Alignment, Band, Placed, BAND_SLACK, DELTA, WILD};

/// Rows within the reads' reach up to which every applicable set of them is a context.
pub(crate) const ALL_SUBSETS: usize = 4;
/// Haplotypes every applicable set of up to [`ALL_SUBSETS`] rows may fill.
pub(crate) const ALL_SUBSETS_CAP: usize = 64;
/// Haplotypes the contexts of more rows may fill.
pub(crate) const MAX_HAPS: usize = 64;


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
        let mut r = ref_start;
        for op in cigar {
            match *op {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    for _ in 0..n {
                        refpos.push(r);
                        r += 1;
                    }
                }
                Cigar::Ins(n) => refpos.extend(std::iter::repeat_n((r - 1).max(ref_start), n as usize)),
                Cigar::SoftClip(n) => {
                    let n = n as i64;
                    if refpos.is_empty() {
                        refpos.extend((0..n).map(|t| r - (n - t)));
                    } else {
                        refpos.extend((0..n).map(|t| r + t));
                    }
                }
                Cigar::Del(n) | Cigar::RefSkip(n) => r += n as i64,
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }
        // A CIGAR that disagrees with the bases' length: pad or cut, never fail.
        refpos.resize(seq.len(), r);
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
        })
    }

    pub(crate) fn len(&self) -> usize {
        self.seq.len()
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
    fn build(site: &Site, rows: &[usize]) -> Hap {
        let mut order: Vec<usize> = rows.to_vec();
        order.sort_by_key(|&r| site.spans[r].0);
        let mut seq = site.reference.clone();
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
        let class = if rows.contains(&0) {
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
    /// Window-relative REF spans of `rows`, and the spans their equivalent placements
    /// cover (a length change's; any other row's REF span).
    spans: Vec<(i64, i64)>,
    places: Vec<(i64, i64)>,
    /// The most any haplotype's rows change the window's length by.
    max_shift: i64,
    pub(crate) haps: Vec<Hap>,
    /// Contexts kept, and whether any candidate context was left out by the cap.
    pub(crate) n_contexts: usize,
    pub(crate) cap_exceeded: bool,
}

/// The rows a site is built from, besides X: the input's other rows within reach (a
/// state where at X's site, else a context), the other alleles the reads show at the
/// locus (states), and the background changes recurring there (contexts only).
#[derive(Debug, Default, Clone)]
pub(crate) struct Candidates {
    pub(crate) listed: Vec<Row>,
    pub(crate) discovered: Vec<Row>,
    pub(crate) background: Vec<Row>,
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
        hap.x_seg.map_or_else(|| hap.to_hap(self.spans[0].0, false), |(hs, _)| hs)
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
            (Class::Ref, _) => on_ref.map(|k| (hap.to_hap(self.spans[0].0 + k as i64, false), x.ref_al[k])).collect(),
            _ => Vec::new(),
        }
    }

    /// The haplotypes of `target` over the window `(start, reference)` (upper case),
    /// from `candidates` and the variant's `reads`. None when the window does not hold
    /// the target's REF.
    pub(crate) fn build(
        target: Row,
        candidates: &Candidates,
        (start, reference): (i64, Vec<u8>),
        reads: &[ReadView],
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
        // The site's states and the contexts' rows, each once.
        let mut seen: BTreeSet<Row> = BTreeSet::from([target.clone()]);
        let (mut states, mut contexts): (Vec<Row>, Vec<Row>) = (Vec::new(), Vec::new());
        let mut take = |r: &Row, state: bool| {
            if holds(r) && seen.insert(r.clone()) {
                if state {
                    states.push(r.clone());
                } else {
                    contexts.push(r.clone());
                }
            }
        };
        for r in &candidates.listed {
            take(r, meets(&target, r));
        }
        if target.ref_al.len() == 1 && target.alt_al.len() == 1 {
            // An SNV's site holds every base (a third base is another allele).
            for &b in b"ACGT" {
                take(&Row::new(target.pos, &target.ref_al, &[b]), true);
            }
        }
        for r in &candidates.discovered {
            take(r, true);
        }
        // A background change at X's place cannot share a haplotype with X: crossed with
        // REF alone it would favour REF, so it is left to the reads' bases there.
        for r in candidates.background.iter().filter(|r| !meets(&target, r)) {
            take(r, false);
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
        let spans = rows.iter().map(|r| (r.pos - start, r.end() - start)).collect();
        let places = rows.iter().map(place).map(|(s, e)| (s - start, e - start)).collect();
        let max_shift = rows.iter().map(|r| r.ref_al.len().abs_diff(r.alt_al.len()) as i64).sum();
        let mut site = Site {
            start,
            reference,
            rows,
            n_states,
            spans,
            places,
            max_shift,
            haps: Vec::new(),
            n_contexts: 0,
            cap_exceeded: false,
        };

        // The context rows a read over X can reach: REF span within X's ± the longest read.
        let l = reads.iter().map(|r| r.len() as i64).max().unwrap_or(0);
        let (xs, xe) = site.spans[0];
        let rest: Vec<usize> = (1 + site.n_states..site.rows.len())
            .filter(|&r| site.spans[r].0 >= xs - l && site.spans[r].1 <= xe + l)
            .collect();
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
    /// one is read as REF's. A length change is scored: `DELTA` better with the row than
    /// on the reference shows it, as much worse refutes it.
    fn evidence(&self, reads: &[ReadView], rows: &[usize]) -> Vec<Evidence> {
        let reference = Hap::build(self, &[]);
        let singles: Vec<Option<Hap>> =
            rows.iter().map(|&r| self.rows[r].changes_length().then(|| Hap::build(self, &[r]))).collect();
        let wn = self.reference.len() as i64;
        reads
            .iter()
            .map(|read| {
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
                for (&r, single) in rows.iter().zip(&singles) {
                    let (rs, re) = self.spans[r];
                    if rs < a || re > b {
                        continue;
                    }
                    match single {
                        None => {
                            let row = &self.rows[r];
                            let reads: Vec<Option<u8>> = (0..row.ref_al.len())
                                .filter(|&k| row.ref_al[k] != row.alt_al[k])
                                .map(|k| at.get(&(rs + k as i64)).copied().filter(|&b| b != WILD).map(|b| {
                                    if b == row.alt_al[k] {
                                        1
                                    } else if b == row.ref_al[k] {
                                        0
                                    } else {
                                        2
                                    }
                                }))
                                .collect();
                            if reads.iter().all(|v| *v == Some(1)) {
                                e.shown.push(r);
                            } else if reads.iter().all(|v| *v == Some(0)) {
                                e.refuted.push(r);
                            }
                        }
                        Some(single) => {
                            let with_r = scorer.score(self, single);
                            let without = scorer.score(self, &reference);
                            if with_r - without >= DELTA {
                                e.shown.push(r);
                            } else if without - with_r >= DELTA {
                                e.refuted.push(r);
                            }
                        }
                    }
                }
                e
            })
            .collect()
    }

    /// Contexts from more than [`ALL_SUBSETS`] rows, in the order they are kept while
    /// there is room: each row, each pair two reads or more show together, and each set
    /// reads show together, ranked by the reads they explain (a context missing a row the
    /// reads carry with another leaves them fitting no haplotype, and the state that
    /// best stands in for the row then wins them), then by fewer rows and nearer X.
    fn ranked_contexts(&self, rest: &[usize], evidence: &[Evidence]) -> Vec<Vec<usize>> {
        let mut candidates: BTreeSet<Vec<usize>> = rest.iter().map(|&r| vec![r]).collect();
        for (i, &a) in rest.iter().enumerate() {
            for &b in &rest[i + 1..] {
                if self.conflict(a, b) {
                    continue;
                }
                let together = evidence
                    .iter()
                    .filter(|e| e.shown.binary_search(&a).is_ok() && e.shown.binary_search(&b).is_ok())
                    .count();
                if together >= 2 {
                    candidates.insert(vec![a, b]);
                }
            }
        }
        candidates.extend(self.shown_together(evidence));
        self.rank(evidence, candidates.into_iter().collect())
    }

    /// Every maximal set of rows that can be applied together shown together by two
    /// reads or more.
    fn shown_together(&self, evidence: &[Evidence]) -> BTreeSet<Vec<usize>> {
        let mut distinct: BTreeMap<&Vec<usize>, usize> = BTreeMap::new();
        for e in evidence.iter().filter(|e| !e.shown.is_empty()) {
            *distinct.entry(&e.shown).or_insert(0) += 1;
        }
        let distinct: Vec<(&Vec<usize>, usize)> = distinct.into_iter().collect();
        let mut together: BTreeSet<Vec<usize>> = BTreeSet::new();
        let mut seen: BTreeSet<Vec<usize>> = BTreeSet::new();
        for a in 0..distinct.len() {
            for b in a..distinct.len() {
                if a == b && distinct[a].1 < 2 {
                    continue;
                }
                let inter: Vec<usize> =
                    distinct[a].0.iter().copied().filter(|r| distinct[b].0.binary_search(r).is_ok()).collect();
                // Many read pairs share an intersection: enumerate each one once.
                if !inter.is_empty() && seen.insert(inter.clone()) {
                    self.maximal_independent_sets(&inter, &mut together);
                }
            }
        }
        together
            .iter()
            .filter(|t| !together.iter().any(|u| u.len() > t.len() && is_subset(t, u)))
            .cloned()
            .collect()
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
    /// its place: rows at one place can all beat the reference in one read).
    fn explains(&self, h: &[usize], e: &Evidence) -> bool {
        !h.iter().any(|r| e.refuted.binary_search(r).is_ok())
            && h.iter().any(|r| e.shown.binary_search(r).is_ok())
            && e.shown.iter().all(|&s| h.iter().any(|&r| r == s || self.conflict(r, s)))
    }

    /// Candidates ranked: reads explained (more first), fewer rows, nearer X, positions.
    fn rank(&self, evidence: &[Evidence], cands: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
        let mut keyed: Vec<(usize, usize, i64, Vec<usize>)> = cands
            .into_iter()
            .map(|h| {
                let support = evidence.iter().filter(|e| self.explains(&h, e)).count();
                let d = h.iter().map(|&r| self.distance_to_target(r)).min().unwrap_or(0);
                (support, h.len(), d, h)
            })
            .collect();
        keyed.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)).then(a.3.cmp(&b.3)));
        keyed.into_iter().map(|k| k.3).collect()
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

/// One read's evidence on the context rows (sorted row indices).
#[derive(Default)]
struct Evidence {
    shown: Vec<usize>,
    refuted: Vec<usize>,
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
        let read = self.read;
        let l = read.len() as i64;
        let (a, b) = hap.reach(read.ref_start - site.start, read.ref_end - site.start, l);
        let wn = hap.window_len();
        let at = |p: i64| hap.to_hap((p - site.start).clamp(0, wn), false) - a as i64;
        let mut centers = Vec::with_capacity(read.len() + 1);
        centers.push(at(read.refpos[0]));
        centers.extend(read.refpos.iter().map(|&p| at(p) + 1));
        let shifts: i64 = near.iter().map(|&r| site.rows[r].ref_al.len().abs_diff(site.rows[r].alt_al.len()) as i64).sum();
        (&hap.seq[a..b], a, Band::around(&centers, BAND_SLACK + shifts, b - a))
    }

    /// The rows within the band's reach: its slack and every length change the site's
    /// rows could carry, twice (the band's width and the read's mapped span).
    fn near(&self, site: &Site, hap: &Hap) -> Vec<usize> {
        self.rows_within(site, hap, BAND_SLACK + 2 * site.max_shift + 2)
    }

    pub(crate) fn score(&mut self, site: &Site, hap: &Hap) -> i32 {
        let wide = self.rows_within(site, hap, self.read.len() as i64 + site.max_shift);
        if let Some(&s) = self.finals.get(&wide) {
            return s;
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
            score::full_score(read, seq)
        };
        self.finals.insert(wide, s);
        s
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
    fn up_to_four_rows_every_subset_meets_every_state() {
        let c = random_seq(600, 9);
        let base = |p: i64| (c[p as usize - 1] as char).to_string();
        let other = |p: i64| if c[p as usize - 1] == b'A' { "C".to_string() } else { "A".to_string() };
        let x = row(300, &base(300), &other(300));
        let near: Vec<Row> = [250, 270, 330, 360].iter().map(|&p| row(p, &base(p), &other(p))).collect();
        let s = site(&c, &x, &near, 100);
        // An SNV's site: REF, X and the two other bases; 16 contexts of four rows.
        assert_eq!(s.n_states, 2);
        assert_eq!(s.n_contexts, 16);
        assert_eq!(s.haps.len(), 16 * 4);
        let alt = s.haps.iter().filter(|h| h.class == Class::Alt).count();
        let reference = s.haps.iter().filter(|h| h.class == Class::Ref).count();
        assert_eq!((alt, reference), (16, 16), "every context with X and with REF");
        // A row out of the reads' reach is no context.
        let far = row(450, &base(450), &other(450));
        let s = site(&c, &x, &[near[0].clone(), far], 100);
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
        let s = Site::build(x, &candidates, (0, c.clone()), &reads).unwrap();
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
        let s = Site::build(x, &candidates, (0, c.clone()), &reads).unwrap();
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
        let s = site(&c, &x, &[del_at_x.clone(), ctx_snv, ctx_del], 100);
        assert_eq!(s.n_states, 3, "the two other bases and the deletion over X");
        // Contexts: {}, {snv}, {del} (the two context rows meet): 3 contexts × 5 states.
        assert_eq!(s.n_contexts, 3);
        assert_eq!(s.haps.len(), 15);
        // Two insertions left-aligned a base apart are one site.
        let c = contig(23, 299, b"AC");
        let ins1 = row(300, "A", "AGT");
        let ins2 = row(300, "A", "ATTA");
        let s = site(&c, &ins1, &[ins2], 100);
        assert_eq!((s.n_states, s.n_contexts, s.haps.len()), (1, 1, 3));
    }
}
