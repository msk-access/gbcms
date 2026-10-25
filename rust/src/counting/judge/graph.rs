//! The graph scorer (the judge's default; `GBCMS_JUDGE_SCORER=enumerated` selects the
//! enumeration it replaces): per state of X's site, one dynamic programme of the read over
//! a small local sequence graph, in place of enumerating the contexts of the other rows and
//! scoring the read against each haplotype.
//!
//! - **The graph** is segments joined by edges: reference segments, and a segment per
//!   allele a row spells, its edges leaving the reference before the row's REF span and
//!   rejoining it after (a bubble the read may take or skip). A masked base is a wildcard
//!   in the reference segment. A path through the graph spells one haplotype of the window,
//!   so the best path is the best of every combination of the rows at once. Segments and
//!   edges are the whole model: a splice (RNA) or a breakend junction (an SV) is one more
//!   edge, joining two reference segments that are not adjacent.
//! - **The programme** is the linear scorer's (`score::run`) over a DAG: the read consumed
//!   whole, the graph's ends free, affine gaps, a read base left over at either end
//!   costing as the error model says; costs come from one [`ErrorModel`] so a long-read
//!   model can stand in for [`ShortReads`].
//! - **The band:** each graph base holds a reference coordinate (an allele's bases spread
//!   over the REF span they replace); read base i meets the graph bases whose coordinate
//!   lies within the band's half-width of its own aligned position. As in the linear
//!   scorer, the band's result stands only where its best path stays off the band's edge,
//!   scores above what leaving the band costs and (short of exact) no copy of the read's
//!   bases lies on the graph off the band; else a wider band, then every base in the
//!   read's reach.

use super::score::{self, BAND_SLACK, COPY_K, GAP_EXTEND, GAP_OPEN, WILD};

/// The costs a read is scored with: what a read base costs against a graph base, what a
/// read base left over at an end of the alignment costs, and a gap's opening and
/// extension. One interface for every read technology's error model.
pub(crate) trait ErrorModel {
    fn substitution(&self, read_base: u8, graph_base: u8) -> i32;
    fn overhang(&self, read_base: u8) -> i32;
    /// A gap's first base: opening plus one extension.
    fn gap_open(&self) -> i32;
    fn gap_extend(&self) -> i32;
}

/// Short reads: a confident mismatch 4, a gap of k bases 6 + k, a masked base (below
/// `--min-baseq`, or N) nothing; a base left over at an end as a mismatch.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct ShortReads;

impl ErrorModel for ShortReads {
    #[inline]
    fn substitution(&self, read_base: u8, graph_base: u8) -> i32 {
        score::substitution_cost(read_base, graph_base)
    }
    #[inline]
    fn overhang(&self, read_base: u8) -> i32 {
        score::overhang_cost(read_base)
    }
    #[inline]
    fn gap_open(&self) -> i32 {
        GAP_OPEN + GAP_EXTEND
    }
    #[inline]
    fn gap_extend(&self) -> i32 {
        GAP_EXTEND
    }
}

/// Unreachable.
const NEG: i32 = i32::MIN / 4;

/// One segment of the graph: the reference's bases over a span, or the allele a row
/// spells in place of the REF bases over its span.
#[derive(Debug, Clone)]
pub(crate) struct Segment {
    /// The row (an index into the site's rows) whose allele it spells; None for the
    /// reference.
    pub(crate) row: Option<usize>,
    /// The reference span it stands for, window coordinates `[s, e)`.
    pub(crate) span: (i64, i64),
    /// Its bases, indices into the graph's bases `[first, after)`.
    first: u32,
    after: u32,
}

/// A segment before it is laid out: its row, its reference span and its bases.
type Pending = (Option<usize>, (i64, i64), Vec<u8>);

/// Segments and the edges between them, before the bases are laid out.
#[derive(Debug, Default)]
pub(crate) struct GraphBuilder {
    segments: Vec<Pending>,
    edges: Vec<(usize, usize)>,
}

impl GraphBuilder {
    /// A segment: the bases spelled over the reference span `span` (window coordinates),
    /// for `row` (None: the reference). Its id.
    pub(crate) fn segment(&mut self, row: Option<usize>, span: (i64, i64), bases: Vec<u8>) -> usize {
        self.segments.push((row, span, bases));
        self.segments.len() - 1
    }

    /// An edge from the end of segment `from` to the start of `to`. Edges run forward: `to`
    /// starts where `from` ends or later (a splice or a deletion junction jumps ahead).
    pub(crate) fn edge(&mut self, from: usize, to: usize) {
        self.edges.push((from, to));
    }

    /// The graph: segments in order of their spans (a topological order, as every edge
    /// runs forward), their bases laid out, each base's predecessors (the reference's
    /// first, so a traceback prefers it on a tie).
    pub(crate) fn finish(self) -> Graph {
        let mut order: Vec<usize> = (0..self.segments.len()).filter(|&s| !self.segments[s].2.is_empty()).collect();
        // an empty segment (an allele of no bases) joins its neighbours directly
        let empty: Vec<bool> = self.segments.iter().map(|s| s.2.is_empty()).collect();
        order.sort_by_key(|&s| (self.segments[s].1 .0, self.segments[s].1 .1, self.segments[s].0.is_some(), s));
        let mut at = vec![u32::MAX; self.segments.len()];
        let mut segments = Vec::with_capacity(order.len());
        let (mut bases, mut coord, mut seg_of) = (Vec::new(), Vec::new(), Vec::new());
        for (k, &s) in order.iter().enumerate() {
            let (row, span, b) = &self.segments[s];
            at[s] = k as u32;
            let first = bases.len() as u32;
            let width = span.1 - span.0;
            for (j, &x) in b.iter().enumerate() {
                bases.push(x);
                // the allele's bases spread over the REF span they replace
                coord.push(span.0 + (j as i64 * width) / b.len() as i64);
                seg_of.push(k as u32);
            }
            segments.push(Segment { row: *row, span: *span, first, after: bases.len() as u32 });
        }
        // Per segment, its predecessors (empty ones resolved through), the reference first.
        let mut into: Vec<Vec<usize>> = vec![Vec::new(); self.segments.len()];
        for &(a, b) in &self.edges {
            debug_assert!(self.segments[b].1 .0 >= self.segments[a].1 .1 || empty[a] || empty[b], "edges run forward");
            into[b].push(a);
        }
        fn ends(s: usize, into: &[Vec<usize>], empty: &[bool], out: &mut Vec<usize>, depth: usize) {
            if !empty[s] {
                out.push(s);
            } else if depth < 64 {
                for &p in &into[s] {
                    ends(p, into, empty, out, depth + 1);
                }
            }
        }
        let n = bases.len();
        let mut pred_off = Vec::with_capacity(n + 1);
        let mut preds: Vec<u32> = Vec::new();
        pred_off.push(0u32);
        for &s in &order {
            let k = at[s] as usize;
            let seg = &segments[k];
            let mut from: Vec<usize> = Vec::new();
            for &p in &into[s] {
                ends(p, &into, &empty, &mut from, 0);
            }
            from.sort_by_key(|&p| (segments[at[p] as usize].row.is_some(), at[p]));
            from.dedup();
            for j in seg.first..seg.after {
                if j == seg.first {
                    preds.extend(from.iter().map(|&p| segments[at[p] as usize].after - 1));
                } else {
                    preds.push(j - 1);
                }
                pred_off.push(preds.len() as u32);
            }
        }
        // successors, for the band's edge
        let mut succ_count = vec![0u32; n + 1];
        for &p in &preds {
            succ_count[p as usize + 1] += 1;
        }
        for v in 0..n {
            succ_count[v + 1] += succ_count[v];
        }
        let mut succs = vec![0u32; preds.len()];
        let mut fill = succ_count.clone();
        for v in 0..n {
            for &p in &preds[pred_off[v] as usize..pred_off[v + 1] as usize] {
                succs[fill[p as usize] as usize] = v as u32;
                fill[p as usize] += 1;
            }
        }
        // each row's length change once, for the band's half-width
        let mut shifts: Vec<(i64, i64, i64)> = Vec::new();
        let mut seen = std::collections::BTreeSet::new();
        for seg in &segments {
            if let Some(r) = seg.row {
                if seen.insert(r) {
                    let d = ((seg.after - seg.first) as i64 - (seg.span.1 - seg.span.0)).abs();
                    if d > 0 {
                        shifts.push((seg.span.0, seg.span.1, d));
                    }
                }
            }
        }
        let mut end_max = Vec::with_capacity(segments.len());
        for seg in &segments {
            let last: i64 = end_max.last().copied().unwrap_or(i64::MIN);
            end_max.push(last.max(seg.span.1));
        }
        Graph { segments, bases, coord, seg_of, pred_off, preds, succ_off: succ_count, succs, shifts, end_max }
    }
}

/// A local sequence graph laid out for the programme: its bases in a topological order.
#[derive(Debug, Clone, Default)]
pub(crate) struct Graph {
    pub(crate) segments: Vec<Segment>,
    pub(crate) bases: Vec<u8>,
    /// Per base, the reference coordinate it stands at (window coordinates).
    pub(crate) coord: Vec<i64>,
    seg_of: Vec<u32>,
    pred_off: Vec<u32>,
    preds: Vec<u32>,
    succ_off: Vec<u32>,
    succs: Vec<u32>,
    /// Per row the graph holds (once), its REF span and the bases it changes the length by.
    shifts: Vec<(i64, i64, i64)>,
    /// Per segment (in order), the furthest end of it and those before it.
    end_max: Vec<i64>,
}

impl Graph {
    pub(crate) fn len(&self) -> usize {
        self.bases.len()
    }

    #[inline]
    pub(crate) fn preds(&self, v: usize) -> &[u32] {
        &self.preds[self.pred_off[v] as usize..self.pred_off[v + 1] as usize]
    }

    #[inline]
    fn succs(&self, v: usize) -> &[u32] {
        &self.succs[self.succ_off[v] as usize..self.succ_off[v + 1] as usize]
    }

    /// The segment base `v` lies in.
    pub(crate) fn segment_of(&self, v: usize) -> &Segment {
        &self.segments[self.seg_of[v] as usize]
    }

    /// The bases `[v0, v1)` of the segments meeting the reference span `trim`.
    fn range(&self, trim: (i64, i64)) -> (usize, usize) {
        let k0 = self.end_max.partition_point(|&e| e <= trim.0);
        let k1 = self.segments.partition_point(|s| s.span.0 < trim.1);
        if k0 >= k1 {
            return (0, 0);
        }
        (self.segments[k0].first as usize, self.segments[k1 - 1].after as usize)
    }

    /// Whether base `v` lies in the read's reach `trim`: a reference base by its
    /// coordinate, an allele's with its whole segment where it meets the reach.
    #[inline]
    fn in_reach(&self, v: usize, trim: (i64, i64)) -> bool {
        let seg = self.segment_of(v);
        if seg.row.is_some() {
            seg.span.0 < trim.1 && seg.span.1 > trim.0
        } else {
            trim.0 <= self.coord[v] && self.coord[v] < trim.1
        }
    }

    /// The length changes of the rows within `by` of the reference span `[lo, hi]`.
    pub(crate) fn shifts_near(&self, lo: i64, hi: i64, by: i64) -> i64 {
        self.shifts.iter().filter(|&&(s, e, _)| s < hi + by && e > lo - by).map(|s| s.2).sum()
    }
}

/// Per graph base in the read's reach, the programme's rows (read bases consumed,
/// `1..=n`) it holds: `[lo, hi]`, empty where `lo > hi`; and where each base's cells start
/// in the arrays. Bases `[v0, v1)` cover the reach; a base outside holds no row.
#[derive(Debug, Clone)]
pub(crate) struct GraphBand {
    v0: usize,
    v1: usize,
    lo: Vec<u32>,
    hi: Vec<u32>,
    off: Vec<u32>,
    cells: usize,
    trim: (i64, i64),
}

impl GraphBand {
    /// The band `half` coordinates each side of each read base's aligned position (`rp`,
    /// window coordinates, non-decreasing), over the bases in the read's reach `trim`;
    /// None for every row of every base in reach (the full programme).
    pub(crate) fn new(g: &Graph, rp: &[i64], half: Option<i64>, trim: (i64, i64)) -> GraphBand {
        let n = rp.len() as u32;
        let (v0, v1) = g.range(trim);
        let len = v1 - v0;
        let (mut lo, mut hi, mut off) = (Vec::with_capacity(len), Vec::with_capacity(len), Vec::with_capacity(len + 1));
        let mut cells = 0usize;
        for v in v0..v1 {
            let c = g.coord[v];
            let (a, b) = if !g.in_reach(v, trim) || n == 0 {
                (1, 0)
            } else if let Some(h) = half {
                (1 + rp.partition_point(|&p| p < c - h) as u32, rp.partition_point(|&p| p <= c + h) as u32)
            } else {
                (1, n)
            };
            off.push(cells as u32);
            if a <= b {
                cells += (b - a + 1) as usize;
            }
            lo.push(a);
            hi.push(b);
        }
        off.push(cells as u32);
        GraphBand { v0, v1, lo, hi, off, cells, trim }
    }

    #[inline]
    fn holds(&self, v: usize, i: u32) -> bool {
        v >= self.v0 && v < self.v1 && self.lo[v - self.v0] <= i && i <= self.hi[v - self.v0]
    }

    #[inline]
    fn cell(&self, v: usize, i: u32) -> Option<usize> {
        self.holds(v, i).then(|| (self.off[v - self.v0] + i - self.lo[v - self.v0]) as usize)
    }

    /// Whether base `v` is in the read's reach at all: a base out of reach is no neighbour.
    #[inline]
    fn reached(&self, g: &Graph, v: usize) -> bool {
        v >= self.v0 && v < self.v1 && g.in_reach(v, self.trim)
    }
}

/// The best alignment of a read in a graph: its score, where it ends (row, base), and
/// whether its path (or one as good) touched the band's edge.
#[derive(Debug, Clone, Copy)]
pub(crate) struct GraphRun {
    pub(crate) score: i32,
    end: Option<(u32, u32)>,
    clipped: bool,
}

/// Traceback bits per cell: 0-1 where the cell's best path came from (0 the diagonal, 1 a
/// deletion, 2 an insertion), 2 a deletion extended, 3 an insertion extended, 4 the
/// diagonal started there (the read bases before it left over).
const FROM_E: u8 = 1;
const FROM_F: u8 = 2;
const E_EXT: u8 = 4;
const F_EXT: u8 = 8;
const D_START: u8 = 16;

/// Reusable arrays for the programme.
#[derive(Debug, Default)]
pub(crate) struct Workspace {
    a: Vec<i32>,
    e: Vec<i32>,
    touched: Vec<u8>,
    dir: Vec<u8>,
    dpred: Vec<u8>,
    epred: Vec<u8>,
    lead: Vec<i32>,
}

/// The programme over `band`. `TRACK`: follow whether paths touch the band's edge.
/// `TRACE`: keep the directions for [`traceback`].
///
/// Per cell (row i, base v): D, a path whose last move aligned read base i - 1 to v
/// (after an aligned path at a predecessor, or starting there after the read bases before
/// it left over); E, a deletion (v skipped after an aligned base); F, an insertion (read
/// base i - 1 after an aligned base at v); A the best of them. A path ends on an aligned
/// base, the read bases after it left over.
fn run<M: ErrorModel, const TRACK: bool, const TRACE: bool>(
    m: &M,
    read: &[u8],
    g: &Graph,
    band: &GraphBand,
    ws: &mut Workspace,
) -> GraphRun {
    let n = read.len() as u32;
    let (open, ext) = (m.gap_open(), m.gap_extend());
    ws.lead.clear();
    ws.lead.push(0);
    for &b in read {
        let last = *ws.lead.last().unwrap();
        ws.lead.push(last + m.overhang(b));
    }
    let lead = &ws.lead;
    let cells = band.cells;
    ws.a.clear();
    ws.a.resize(cells, NEG);
    ws.e.clear();
    ws.e.resize(cells, NEG);
    if TRACK {
        ws.touched.clear();
        ws.touched.resize(cells, 0);
    }
    if TRACE {
        for v in [&mut ws.dir, &mut ws.dpred, &mut ws.epred] {
            v.clear();
            v.resize(cells, 0);
        }
    }
    let total_lead = lead[n as usize];
    let (mut score, mut end, mut clipped): (i32, Option<(u32, u32)>, bool) = (total_lead, None, false);
    for v in band.v0..band.v1 {
        let (lo, hi) = (band.lo[v - band.v0], band.hi[v - band.v0]);
        if lo > hi {
            continue;
        }
        let base = g.bases[v];
        let preds = g.preds(v);
        let (mut f, mut tf_prev) = (NEG, 0u8);
        for i in lo..=hi {
            let idx = (band.off[v - band.v0] + i - lo) as usize;
            // D: from a predecessor's aligned path at row i - 1, or starting here
            let (mut diag, mut dk, mut dt) = (NEG, 0usize, 0u8);
            if i >= 2 {
                for (k, &u) in preds.iter().enumerate() {
                    if let Some(x) = band.cell(u as usize, i - 1) {
                        let val = ws.a[x];
                        if val > diag {
                            (diag, dk) = (val, k);
                            if TRACK {
                                dt = ws.touched[x] & 1;
                            }
                        } else if TRACK && val == diag && val > NEG / 2 {
                            dt |= ws.touched[x] & 1;
                        }
                    }
                }
            }
            let start = lead[i as usize - 1];
            let sub = m.substitution(read[i as usize - 1], base);
            let (d, from_a) = if diag >= start { (diag + sub, true) } else { (start + sub, false) };
            // E: v skipped after a predecessor's aligned base (or deletion) at row i
            let (mut en, mut ek, mut e_ext, mut te) = (NEG, 0usize, false, 0u8);
            for (k, &u) in preds.iter().enumerate() {
                if let Some(x) = band.cell(u as usize, i) {
                    let (eo, ee) = (ws.a[x] - open, ws.e[x] - ext);
                    let val = eo.max(ee);
                    let tt = if TRACK {
                        let t = ws.touched[x];
                        (u8::from(eo == val) & t) | (u8::from(ee == val) & (t >> 2) & 1)
                    } else {
                        0
                    };
                    if val > en {
                        (en, ek, e_ext, te) = (val, k, ee > eo, tt);
                    } else if TRACK && val == en {
                        te |= tt;
                    }
                }
            }
            // F: read base i - 1 inserted after an aligned base at v
            let (fo, fe) = if i > lo { (ws.a[idx - 1] - open, f - ext) } else { (NEG, NEG) };
            let fnew = fo.max(fe);
            let a = d.max(en).max(fnew);
            let mut td = 0u8;
            if TRACK {
                let edge = u8::from(edge_cell(g, band, v, i, lo, hi, n));
                let tf = edge | if i > lo { (u8::from(fo == fnew) & ws.touched[idx - 1] & 1) | (u8::from(fe == fnew) & tf_prev) } else { 0 };
                te |= edge;
                td = edge | if from_a { dt } else { 0 };
                let ta = (u8::from(d == a) & td) | (u8::from(en == a) & te) | (u8::from(fnew == a) & tf);
                ws.touched[idx] = (ta & 1) | ((tf & 1) << 1) | ((te & 1) << 2);
                tf_prev = tf & 1;
            }
            if TRACE {
                let src = if d == a {
                    0
                } else if en == a {
                    FROM_E
                } else {
                    FROM_F
                };
                ws.dir[idx] = src
                    | if e_ext { E_EXT } else { 0 }
                    | if fe > fo { F_EXT } else { 0 }
                    | if from_a { 0 } else { D_START };
                ws.dpred[idx] = dk as u8;
                ws.epred[idx] = ek as u8;
            }
            ws.a[idx] = a;
            ws.e[idx] = en;
            f = fnew;
            if d > NEG / 2 {
                // the best end, the first (row, then base) on a tie, as the linear scorer's;
                // a tie with the read left over whole keeps that (no base aligned)
                let total = d + (total_lead - lead[i as usize]);
                if total > score {
                    (score, end, clipped) = (total, Some((i, v as u32)), td != 0);
                } else if total == score {
                    clipped |= td != 0;
                    if end.is_some_and(|e| (i, v as u32) < e) {
                        end = Some((i, v as u32));
                    }
                }
            }
        }
    }
    GraphRun { score, end, clipped }
}

/// Whether a cell has a neighbour in the read's reach that the band leaves out: a path
/// through it may leave the band, or enter it from outside.
#[inline]
fn edge_cell(g: &Graph, band: &GraphBand, v: usize, i: u32, lo: u32, hi: u32, n: u32) -> bool {
    if (i == lo && lo > 1) || (i == hi && hi < n) {
        return true;
    }
    let out = |u: usize, r: u32| !band.holds(u, r) && band.reached(g, u);
    g.preds(v).iter().any(|&u| out(u as usize, i) || (i >= 2 && out(u as usize, i - 1)))
        || g.succs(v).iter().any(|&w| out(w as usize, i) || (i < n && out(w as usize, i + 1)))
}

/// The lowest score the band answers for by its edge alone (as the linear scorer's).
fn floor_for(slack: i64) -> i32 {
    -(GAP_OPEN + GAP_EXTEND * (slack as i32 + 1))
}

/// Whether one of the read's exact k-mers (at its start, middle and end, unmasked) lies
/// along a path of the graph ending where the band does not hold it: a copy of the read's
/// bases elsewhere, which an alignment reaches with no indel. Shift-And over the DAG.
fn copy_off_band(read: &[u8], g: &Graph, band: &GraphBand) -> bool {
    let n = read.len();
    if n < COPY_K {
        return false;
    }
    let mut offsets = vec![0, (n - COPY_K) / 2, n - COPY_K];
    offsets.dedup();
    let top = 1u32 << (COPY_K - 1);
    let mut state = vec![0u32; band.v1 - band.v0];
    offsets.into_iter().any(|o| {
        let kmer = &read[o..o + COPY_K];
        if kmer.contains(&WILD) {
            return false;
        }
        let mask = |b: u8| -> u32 {
            if b == WILD {
                return u32::MAX;
            }
            kmer.iter().enumerate().filter(|(_, &k)| k == b).fold(0u32, |m, (j, _)| m | (1 << j))
        };
        let masks = [b'A', b'C', b'G', b'T'].map(mask);
        let pick = |b: u8| match b {
            b'A' => masks[0],
            b'C' => masks[1],
            b'G' => masks[2],
            b'T' => masks[3],
            _ => u32::MAX,
        };
        let row = (o + COPY_K) as u32;
        for v in band.v0..band.v1 {
            if !band.reached(g, v) {
                state[v - band.v0] = 0;
                continue;
            }
            let from =
                g.preds(v).iter().filter(|&&u| band.reached(g, u as usize)).fold(0u32, |s, &u| s | state[u as usize - band.v0]);
            let s = ((from << 1) | 1) & pick(g.bases[v]);
            state[v - band.v0] = s;
            if s & top != 0 && !band.holds(v, row) {
                return true;
            }
        }
        false
    })
}

/// A read as the graph scorer reads it: its bases (masked to [`WILD`]) and each base's
/// aligned position in window coordinates, non-decreasing.
pub(crate) struct GraphRead<'r> {
    pub(crate) seq: &'r [u8],
    pub(crate) rp: Vec<i64>,
    /// Its reference span (soft clips extrapolated), window coordinates.
    pub(crate) span: (i64, i64),
}

/// Half-width of the wider band tried before the full programme (past the length changes).
const WIDE_SLACK: i64 = 4 * BAND_SLACK;

/// The read's reach in the graph (its reference span widened by its length and by the
/// length changes of the rows near it: a deletion's far side lies that much further along
/// the reference) and the half-widths of the band and the wider band.
fn bands(g: &Graph, read: &GraphRead) -> ((i64, i64), [i64; 2]) {
    let l = read.seq.len() as i64;
    // the rows the read can pass: those over its span, widened by what they shift it
    let mut shifts = 0;
    loop {
        let more = g.shifts_near(read.span.0, read.span.1, BAND_SLACK + 2 + 2 * shifts);
        if more == shifts {
            break;
        }
        shifts = more;
    }
    let trim = (read.span.0 - l - shifts, read.span.1 + l + shifts);
    (trim, [BAND_SLACK + shifts, WIDE_SLACK + shifts])
}

/// The read's best score in `g`: in the band where it stands, else the wider band where
/// it stands, else the full programme over the read's reach.
pub(crate) fn score<M: ErrorModel>(m: &M, read: &GraphRead, g: &Graph, ws: &mut Workspace) -> i32 {
    decided(m, read, g, ws).0
}

/// [`score`], with the band that decided it (its half-width; None for the full
/// programme), which [`path`] traces in.
pub(crate) fn decided<M: ErrorModel>(m: &M, read: &GraphRead, g: &Graph, ws: &mut Workspace) -> (i32, Option<i64>) {
    if read.seq.is_empty() {
        return (0, None);
    }
    let (trim, halves) = bands(g, read);
    for (half, slack) in halves.into_iter().zip([BAND_SLACK, WIDE_SLACK]) {
        let band = GraphBand::new(g, &read.rp, Some(half), trim);
        let r = run::<M, true, false>(m, read.seq, g, &band, ws);
        if !r.clipped && r.score >= floor_for(slack) && (r.score == 0 || !copy_off_band(read.seq, g, &band)) {
            return (r.score, Some(half));
        }
    }
    let full = GraphBand::new(g, &read.rp, None, trim);
    (run::<M, false, false>(m, read.seq, g, &full, ws).score, None)
}

/// The best path of the read in a graph: its score and the rows whose alleles it takes.
#[derive(Debug, Clone)]
pub(crate) struct GraphPath {
    pub(crate) score: i32,
    pub(crate) rows: Vec<usize>,
}

/// The read's best path in `g`, traced back in the band [`decided`] chose (`half`).
pub(crate) fn path<M: ErrorModel>(m: &M, read: &GraphRead, g: &Graph, half: Option<i64>, ws: &mut Workspace) -> GraphPath {
    if read.seq.is_empty() {
        return GraphPath { score: 0, rows: Vec::new() };
    }
    let (trim, _) = bands(g, read);
    let band = GraphBand::new(g, &read.rp, half, trim);
    let r = run::<M, false, true>(m, read.seq, g, &band, ws);
    traceback(g, &band, ws, r)
}

fn traceback(g: &Graph, band: &GraphBand, ws: &Workspace, r: GraphRun) -> GraphPath {
    let mut visited: Vec<u32> = Vec::new();
    let Some((mut i, mut v)) = r.end else {
        return GraphPath { score: r.score, rows: Vec::new() };
    };
    // states: 0 D, 1 A, 2 E, 3 F
    let mut state = 0u8;
    while i > 0 {
        let Some(idx) = band.cell(v as usize, i) else { break };
        let dir = ws.dir[idx];
        match state {
            0 => {
                visited.push(v);
                if dir & D_START != 0 || i == 1 {
                    break;
                }
                v = g.preds(v as usize)[ws.dpred[idx] as usize];
                i -= 1;
                state = 1;
            }
            1 => {
                state = match dir & 3 {
                    0 => 0,
                    FROM_E => 2,
                    _ => 3,
                }
            }
            2 => {
                visited.push(v);
                state = if dir & E_EXT != 0 { 2 } else { 1 };
                v = g.preds(v as usize)[ws.epred[idx] as usize];
            }
            _ => {
                state = if dir & F_EXT != 0 { 3 } else { 1 };
                i -= 1;
            }
        }
    }
    let mut rows: Vec<usize> = visited.iter().filter_map(|&v| g.segment_of(v as usize).row).collect();
    rows.sort_unstable();
    rows.dedup();
    GraphPath { score: r.score, rows }
}

/// The full programme over every base of the graph, scores only (for tests).
#[cfg(test)]
pub(crate) fn full_score<M: ErrorModel>(m: &M, read: &[u8], g: &Graph) -> i32 {
    let rp: Vec<i64> = vec![0; read.len()];
    let full = GraphBand::new(g, &rp, None, (i64::MIN / 4, i64::MAX / 4));
    run::<M, false, false>(m, read, g, &full, &mut Workspace::default()).score
}

#[cfg(test)]
mod tests {
    use super::*;

    fn linear(seq: &[u8]) -> Graph {
        let mut b = GraphBuilder::default();
        b.segment(None, (0, seq.len() as i64), seq.to_vec());
        b.finish()
    }

    #[test]
    fn a_linear_graph_scores_as_the_linear_scorer() {
        for (read, hap) in [
            (&b"GGATCC"[..], &b"AAAGGATCCTTT"[..]),
            (b"GGTTCC", b"AAAGGATCCTTT"),
            (b"CTAGGATCCTGA", b"TTCTAGGTCCTGATT"),
            (b"CTAGGTCCTGA", b"TTCTAGGATTCCTGATT"),
            (b"TTTGG", b"AAATTT"),
            (b"NNNNNNAAAT", b"AAATTT"),
        ] {
            assert_eq!(full_score(&ShortReads, read, &linear(hap)), score::full_score(read, hap), "{read:?} in {hap:?}");
        }
    }

    #[test]
    fn a_bubble_scores_the_best_of_its_two_paths() {
        // AAAC[G|T]TTTA: a read with G or with T fits exactly; one with A mismatches once.
        let mut b = GraphBuilder::default();
        let l = b.segment(None, (0, 4), b"AAAC".to_vec());
        let g = b.segment(None, (4, 5), b"G".to_vec());
        let t = b.segment(Some(1), (4, 5), b"T".to_vec());
        let r = b.segment(None, (5, 9), b"TTTA".to_vec());
        for s in [g, t] {
            b.edge(l, s);
            b.edge(s, r);
        }
        let gr = b.finish();
        assert_eq!(full_score(&ShortReads, b"AACGTTT", &gr), 0);
        assert_eq!(full_score(&ShortReads, b"AACTTTT", &gr), 0);
        assert_eq!(full_score(&ShortReads, b"AACATTT", &gr), score::MISMATCH);
    }
}
