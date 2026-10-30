//! The judge's scorer: per state of X's site, one dynamic programme of the read over a
//! small local sequence graph, in place of enumerating the contexts of the other rows and
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
    /// Where a read base costs nothing against its own base or a masked one ([`WILD`]) on
    /// either side and one cost against any other: that cost (so a column's substitutions
    /// are compared bases alone). None where costs depend on the bases otherwise.
    fn uniform_mismatch(&self) -> Option<i32> {
        None
    }
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
    #[inline]
    fn uniform_mismatch(&self) -> Option<i32> {
        Some(score::MISMATCH)
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
        let (mut shifts, mut grows): (Vec<(i64, i64, i64)>, Vec<bool>) = (Vec::new(), Vec::new());
        let mut seen = std::collections::BTreeSet::new();
        for seg in &segments {
            if let Some(r) = seg.row {
                if seen.insert(r) {
                    let d = (seg.after - seg.first) as i64 - (seg.span.1 - seg.span.0);
                    if d != 0 {
                        shifts.push((seg.span.0, seg.span.1, d.abs()));
                        grows.push(d > 0);
                    }
                }
            }
        }
        let mut end_max = Vec::with_capacity(segments.len());
        for seg in &segments {
            let last: i64 = end_max.last().copied().unwrap_or(i64::MIN);
            end_max.push(last.max(seg.span.1));
        }
        // every predecessor laid out before its base (edges run forward)
        let topological = (0..n).all(|v| preds[pred_off[v] as usize..pred_off[v + 1] as usize].iter().all(|&u| (u as usize) < v));
        Graph { segments, bases, coord, seg_of, pred_off, preds, succ_off: succ_count, succs, shifts, grows, end_max, topological }
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
    /// Per entry of `shifts`, whether the row lengthens the haplotype (an insertion).
    grows: Vec<bool>,
    /// Per segment (in order), the furthest end of it and those before it.
    end_max: Vec<i64>,
    /// Whether every base's predecessors lie before it: [`run`] fills a column from the
    /// columns before it, so a graph that is not (none the builder makes) is scored by
    /// [`run_reference`].
    topological: bool,
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

    /// Base `v`'s place in its segment (0 for its first base).
    pub(crate) fn offset_in_segment(&self, v: usize) -> usize {
        v - self.segment_of(v).first as usize
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

    /// Where `other`'s bases first differ from this graph's (its length where none does),
    /// where the two are laid out alike: the same segments (spans, bases' places, whether
    /// each spells a row), coordinates, edges and length changes, so a programme over them
    /// differs only from that base on. None where they are not laid out alike (as REF's
    /// graph, not cut at X's span, and X's).
    pub(crate) fn layout_until(&self, other: &Graph) -> Option<usize> {
        let segments = self.segments.len() == other.segments.len()
            && self.segments.iter().zip(&other.segments).all(|(a, b)| {
                (a.span, a.first, a.after, a.row.is_some()) == (b.span, b.first, b.after, b.row.is_some())
            });
        let alike = segments
            && self.topological
            && other.topological
            && self.coord == other.coord
            && self.seg_of == other.seg_of
            && self.pred_off == other.pred_off
            && self.preds == other.preds
            && self.succ_off == other.succ_off
            && self.succs == other.succs
            && self.shifts == other.shifts
            && self.grows == other.grows
            && self.end_max == other.end_max;
        alike.then(|| self.bases.iter().zip(&other.bases).position(|(a, b)| a != b).unwrap_or(self.bases.len()))
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
    /// The bytes it holds (its arrays' capacities and itself), for the work space's count.
    fn heap_bytes(&self) -> usize {
        (self.lo.capacity() + self.hi.capacity() + self.off.capacity()) * 4 + std::mem::size_of::<GraphBand>()
    }

    /// The band `half` coordinates each side of each read base's aligned position (`rp`,
    /// window coordinates, non-decreasing), over the bases in the read's reach `trim`;
    /// None for every row of every base in reach (the full programme).
    pub(crate) fn new(g: &Graph, rp: &[i64], half: Option<i64>, trim: (i64, i64)) -> GraphBand {
        let n = rp.len() as u32;
        let (v0, v1) = g.range(trim);
        let len = v1 - v0;
        let (mut lo, mut hi, mut off) = (Vec::with_capacity(len), Vec::with_capacity(len), Vec::with_capacity(len + 1));
        let mut cells = 0usize;
        // The rows a base holds are the read bases whose positions lie within `half` of its
        // coordinate: two partition points of `rp`. A segment's coordinates do not fall, so
        // within one the points only move on (found by search at its first base).
        let (mut seg, mut p_lo, mut p_hi) = (u32::MAX, 0usize, 0usize);
        for v in v0..v1 {
            let c = g.coord[v];
            let (a, b) = if !g.in_reach(v, trim) || n == 0 {
                (1, 0)
            } else if let Some(h) = half {
                if g.seg_of[v] != seg {
                    seg = g.seg_of[v];
                    p_lo = rp.partition_point(|&p| p < c - h);
                    p_hi = rp.partition_point(|&p| p <= c + h);
                } else {
                    while p_lo < rp.len() && rp[p_lo] < c - h {
                        p_lo += 1;
                    }
                    while p_hi < rp.len() && rp[p_hi] <= c + h {
                        p_hi += 1;
                    }
                }
                (1 + p_lo as u32, p_hi as u32)
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

    /// [`GraphBand::new`] as first written, two searches per base (for tests).
    #[cfg(test)]
    fn new_reference(g: &Graph, rp: &[i64], half: Option<i64>, trim: (i64, i64)) -> GraphBand {
        let n = rp.len() as u32;
        let (v0, v1) = g.range(trim);
        let (mut lo, mut hi, mut off) = (Vec::new(), Vec::new(), Vec::new());
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

thread_local! {
    /// The thread's workspaces, one per state ([`with_workspaces`]); each grows to the largest
    /// band a read has needed of it, never to the full programme ahead of need.
    static POOL: std::cell::RefCell<Vec<Workspace>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// The bytes this thread's pooled workspaces hold (their capacities, and the bands they keep
/// for the paths, a band states laid out alike share counted once), for the row caches' budget.
pub(crate) fn pool_bytes() -> usize {
    POOL.with(|p| {
        let pool = p.borrow();
        let mut seen: Vec<*const GraphBand> = Vec::new();
        let mut bands = 0;
        for (b, _) in pool.iter().filter_map(|ws| ws.last.as_ref()) {
            if !seen.contains(&std::rc::Rc::as_ptr(b)) {
                seen.push(std::rc::Rc::as_ptr(b));
                bands += b.heap_bytes();
            }
        }
        pool.iter().map(Workspace::heap_bytes).sum::<usize>() + bands
    })
}

/// Frees this thread's pooled workspaces (they grow again as reads need them).
pub(crate) fn free_pool() {
    POOL.with(|p| p.borrow_mut().clear());
}

impl Workspace {
    /// The bytes it holds on the heap (its capacities).
    fn heap_bytes(&self) -> usize {
        let words = self.a.capacity() + self.e.capacity() + self.f.capacity() + self.lead.capacity() + self.col.diag.capacity() + self.col.en.capacity() + self.marks.capacity();
        let bytes = self.touched.capacity() + self.dir.capacity() + self.dpred.capacity() + self.epred.capacity();
        words * 4 + bytes + (self.ends.capacity() + self.interiors.capacity() + self.snap.1.capacity()) * 8 + self.stack.capacity() * 12
    }
}

/// Reusable arrays for the programme: per cell of the band, A, E and F (the insertion), the
/// read's leftover costs, one column's view of its predecessors, the cells the best ends
/// lie at, and the edge search's marks and stack; the reference programme's edge bits and
/// directions.
#[derive(Debug, Default)]
pub(crate) struct Workspace {
    a: Vec<i32>,
    e: Vec<i32>,
    f: Vec<i32>,
    lead: Vec<i32>,
    col: Column,
    ends: Vec<(u32, u32)>,
    /// Per cell, the edge search that last marked it (above the low four bits) and the
    /// kinds it marked (D, E, F, A), so no search clears the marks of the one before.
    marks: Vec<u32>,
    mark_gen: u32,
    stack: Vec<(u32, u32, u8)>,
    interiors: Vec<(u32, u32)>,
    /// The band that decided the last score ([`decided`]) and its run, whose values the
    /// arrays still hold: the path is traced back from them ([`kept_path`]).
    last: Option<(std::rc::Rc<GraphBand>, GraphRun)>,
    /// The programme's best total and ends as they stood before a base's column, kept for
    /// the states laid out alike that start from it ([`decided_states`]).
    snap: (i32, Vec<(u32, u32)>),
    touched: Vec<u8>,
    dir: Vec<u8>,
    dpred: Vec<u8>,
    epred: Vec<u8>,
}

/// One column's view of its predecessors, per row (`0..=n`): the best aligned path one row
/// up (the diagonal) and the best deletion at the row.
#[derive(Debug, Default)]
struct Column {
    diag: Vec<i32>,
    en: Vec<i32>,
}

/// The rows `[i0, i1]` of base `v`'s column (rows `lo..=hi` of `n`) holding no cell with a
/// neighbour in the read's reach that the band leaves out: the cells off the band's edge,
/// one interval (empty where `i0 > i1`). A row is on the edge where it is the column's first
/// (past row 1) or last (short of row n), or where a predecessor in reach lacks the row or
/// the row above it, or a successor in reach lacks the row or the row below it.
#[inline]
fn interior(g: &Graph, band: &GraphBand, v: usize, lo: u32, hi: u32, n: u32) -> (u32, u32) {
    let mut i0 = if lo > 1 { lo + 1 } else { lo };
    let mut i1 = if hi < n { hi - 1 } else { hi };
    for &u in g.preds(v) {
        let u = u as usize;
        if band.reached(g, u) {
            let (ulo, uhi) = (band.lo[u - band.v0], band.hi[u - band.v0]);
            // the row and the row above it (none above row 1)
            i0 = i0.max(if ulo <= 1 { ulo } else { ulo + 1 });
            i1 = i1.min(uhi);
        }
    }
    for &w in g.succs(v) {
        let w = w as usize;
        if band.reached(g, w) {
            let (wlo, whi) = (band.lo[w - band.v0], band.hi[w - band.v0]);
            // the row and the row below it (none below row n)
            i0 = i0.max(wlo);
            i1 = i1.min(if whi >= n { whi } else { whi.saturating_sub(1) });
        }
    }
    (i0, i1)
}

/// The programme over `band`. `TRACK`: whether a best path (or one as good) touches the
/// band's edge (`GraphRun::clipped`).
///
/// Per cell (row i, base v): D, a path whose last move aligned read base i - 1 to v
/// (after an aligned path at a predecessor, or starting there after the read bases before
/// it left over); E, a deletion (v skipped after an aligned base); F, an insertion (read
/// base i - 1 after an aligned base at v); A the best of them. A path ends on an aligned
/// base, the read bases after it left over.
///
/// Column by column (each base's rows in one pass): the predecessors' cells are read a
/// predecessor at a time over the rows it holds, into the column's diagonal and deletion,
/// then the column's rows are filled in order (the insertion runs down the column). A, E
/// and F are kept per cell, so the edge (`TRACK`, [`touches_edge`]) and the traceback
/// ([`traceback`]) are read back from them for the cells they need ([`sources`]) instead of
/// being carried through every cell. Every value, the best end and the edge flag are
/// [`run_reference`]'s, cell for cell: a maximum keeps the first predecessor reaching it,
/// the edge flag is the OR over every tie, and the best end the first by row, then base.
fn run<M: ErrorModel, const TRACK: bool>(m: &M, read: &[u8], g: &Graph, band: &GraphBand, ws: &mut Workspace) -> GraphRun {
    run_from::<M, TRACK>(m, read, g, band, ws, None, None)
}

/// [`run`], its columns before base `start.1` taken from `start.0`'s run over the same
/// band of a graph laid out alike whose bases are the same before it ([`Graph::layout_until`];
/// its best total and ends there from its snapshot), and with `snap_at` its own best total
/// and ends kept as they stand before that base's column. Every cell, the score and its end
/// are [`run`]'s.
fn run_from<M: ErrorModel, const TRACK: bool>(
    m: &M,
    read: &[u8],
    g: &Graph,
    band: &GraphBand,
    ws: &mut Workspace,
    start: Option<(&Workspace, usize)>,
    snap_at: Option<usize>,
) -> GraphRun {
    let n = read.len() as u32;
    let (open, ext) = (m.gap_open(), m.gap_extend());
    ws.lead.clear();
    ws.lead.push(0);
    for &b in read {
        let last = *ws.lead.last().unwrap();
        ws.lead.push(last + m.overhang(b));
    }
    let cells = band.cells;
    // Every cell of the band is written before it is read (a predecessor's column comes
    // first, a row's cell before the next row's), so the arrays need no clearing.
    for v in [&mut ws.a, &mut ws.e, &mut ws.f] {
        if v.len() < cells {
            v.resize(cells, NEG);
        }
    }
    let rows = n as usize + 1;
    for v in [&mut ws.col.diag, &mut ws.col.en] {
        if v.len() < rows {
            v.resize(rows, NEG);
        }
    }
    ws.ends.clear();
    let total_lead = ws.lead[n as usize];
    let mut score = total_lead;
    let mut from = band.v0;
    if let Some((src, q)) = start {
        let q = q.clamp(band.v0, band.v1);
        let cells = band.off[q - band.v0] as usize;
        ws.a[..cells].copy_from_slice(&src.a[..cells]);
        ws.e[..cells].copy_from_slice(&src.e[..cells]);
        ws.f[..cells].copy_from_slice(&src.f[..cells]);
        score = src.snap.0;
        ws.ends.extend_from_slice(&src.snap.1);
        from = q;
    }
    let snap_at = snap_at.map(|q| q.clamp(band.v0, band.v1));
    {
        let Workspace { a, e, f, lead, col, ends, snap, .. } = &mut *ws;
        let lead = &lead[..];
        for v in from..band.v1 {
            if snap_at == Some(v) {
                *snap = (score, ends.clone());
            }
            let (lo, hi) = (band.lo[v - band.v0], band.hi[v - band.v0]);
            if lo > hi {
                continue;
            }
            let len = (hi - lo + 1) as usize;
            let lo_u = lo as usize;
            let off = band.off[v - band.v0] as usize;
            let base = g.bases[v];
            // The predecessors' diagonal (row i from row i - 1, from row 2 on) and deletion,
            // over this column's rows (index j: row lo + j).
            let diag = &mut col.diag[lo_u..lo_u + len];
            let en = &mut col.en[lo_u..lo_u + len];
            diag.fill(NEG);
            en.fill(NEG);
            for &u in g.preds(v) {
                let u = u as usize;
                if u < band.v0 || u >= band.v1 {
                    continue;
                }
                let (ulo, uhi) = (band.lo[u - band.v0], band.hi[u - band.v0]);
                if ulo > uhi {
                    continue;
                }
                let uoff = band.off[u - band.v0] as usize;
                // diagonal: rows i in [max(lo, ulo + 1, 2), min(hi, uhi + 1)] read row i - 1
                let (d0, d1) = (lo.max(ulo + 1).max(2), hi.min(uhi + 1));
                if d0 <= d1 {
                    let cnt = (d1 - d0 + 1) as usize;
                    let (j0, x0) = ((d0 - lo) as usize, uoff + (d0 - 1 - ulo) as usize);
                    for (dst, &val) in diag[j0..j0 + cnt].iter_mut().zip(&a[x0..x0 + cnt]) {
                        *dst = (*dst).max(val);
                    }
                }
                // deletion: rows i in [max(lo, ulo), min(hi, uhi)] read row i
                let (e0, e1) = (lo.max(ulo), hi.min(uhi));
                if e0 <= e1 {
                    let cnt = (e1 - e0 + 1) as usize;
                    let (j0, x0) = ((e0 - lo) as usize, uoff + (e0 - ulo) as usize);
                    for ((dst, &sa), &se) in en[j0..j0 + cnt].iter_mut().zip(&a[x0..x0 + cnt]).zip(&e[x0..x0 + cnt]) {
                        *dst = (*dst).max((sa - open).max(se - ext));
                    }
                }
            }
            // The column's rows: D and the best of D and E first, row by row; then the
            // insertion down the column, as the best opening above with its extensions:
            // F at row t is max over k < t of (D-or-E at k) - open - (t - 1 - k) * ext. An
            // insertion is never opened from one (open is at least ext), so the best of A
            // above, F included, is the best of D and E above: the same F, row for row, as
            // extending it a row at a time.
            let a_col = &mut a[off..off + len];
            let e_col = &mut e[off..off + len];
            let f_col = &mut f[off..off + len];
            let (lead_up, lead_at, rd) = (&lead[lo_u - 1..lo_u - 1 + len], &lead[lo_u..lo_u + len], &read[lo_u - 1..lo_u - 1 + len]);
            let (diag, en) = (&*diag, &*en);
            // D (kept in F's place for now), the best of D and E (in A's), and the column's
            // best total
            let mut col_best = NEG;
            let rows = diag.iter().zip(lead_up).zip(lead_at).zip(rd).zip(en.iter());
            let cells = a_col.iter_mut().zip(e_col.iter_mut()).zip(f_col.iter_mut());
            // the substitution as compared bases where the model allows (no branch per row)
            let mismatch = m.uniform_mismatch().map(|c| if base == WILD { 0 } else { c });
            for (((((&dg, &lu), &la), &r), &enj), ((ao, eo), fo)) in rows.zip(cells) {
                let sub = match mismatch {
                    Some(c) => {
                        if r == base || r == WILD {
                            0
                        } else {
                            c
                        }
                    }
                    None => m.substitution(r, base),
                };
                let d = dg.max(lu) + sub;
                *ao = d.max(enj);
                *eo = enj;
                *fo = d;
                col_best = col_best.max(d - la);
            }
            // the best ends: every cell scoring the best total (a tie with the read left over
            // whole keeps that, no base aligned)
            if col_best + total_lead >= score {
                for j in 0..len {
                    let total = f_col[j] + (total_lead - lead_at[j]);
                    if total >= score {
                        if total > score {
                            score = total;
                            ends.clear();
                        }
                        ends.push((lo + j as u32, v as u32));
                    }
                }
            }
            if open >= ext {
                let mut best = NEG;
                for j in 0..len {
                    let a0 = a_col[j];
                    let fnew = if j > 0 { best - open - (j as i32 - 1) * ext } else { NEG };
                    a_col[j] = a0.max(fnew);
                    f_col[j] = fnew;
                    best = best.max(a0 + j as i32 * ext);
                }
            } else {
                let (mut fv, mut a_prev) = (NEG, NEG);
                for j in 0..len {
                    let fnew = if j > 0 { (a_prev - open).max(fv - ext) } else { NEG };
                    let av = a_col[j].max(fnew);
                    a_col[j] = av;
                    f_col[j] = fnew;
                    a_prev = av;
                    fv = fnew;
                }
            }
        }
    }
    if snap_at == Some(band.v1) {
        ws.snap = (score, ws.ends.clone());
    }
    // The first end by row, then base; none where the read left over whole scores as well.
    let end = if score > total_lead { ws.ends.iter().copied().min() } else { None };
    let clipped = TRACK && touches_edge(m, read, g, band, ws);
    GraphRun { score, end, clipped }
}

/// A cell's sources, read back from the programme's arrays: its diagonal (the best of its
/// predecessors' A one row up, [`NEG`] where none holds that row), D and whether D extends
/// an aligned path (else it starts there), E, the insertion's opening and extension and F,
/// and A.
#[derive(Debug, Clone, Copy)]
struct Sources {
    diag: i32,
    d: i32,
    from_a: bool,
    en: i32,
    fo: i32,
    fe: i32,
    f: i32,
    a: i32,
}

#[inline]
fn sources<M: ErrorModel>(m: &M, read: &[u8], g: &Graph, band: &GraphBand, ws: &Workspace, v: usize, i: u32) -> Sources {
    let (open, ext) = (m.gap_open(), m.gap_extend());
    let lo = band.lo[v - band.v0];
    let idx = (band.off[v - band.v0] + i - lo) as usize;
    let mut diag = NEG;
    if i >= 2 {
        for &u in g.preds(v) {
            if let Some(x) = band.cell(u as usize, i - 1) {
                diag = diag.max(ws.a[x]);
            }
        }
    }
    let start = ws.lead[i as usize - 1];
    let sub = m.substitution(read[i as usize - 1], g.bases[v]);
    let (d, from_a) = if diag >= start { (diag + sub, true) } else { (start + sub, false) };
    let (fo, fe) = if i > lo { (ws.a[idx - 1] - open, ws.f[idx - 1] - ext) } else { (NEG, NEG) };
    Sources { diag, d, from_a, en: ws.e[idx], fo, fe, f: ws.f[idx], a: ws.a[idx] }
}

/// Whether a best end's path (or one as good) touches the band's edge: from each cell
/// scoring the best total, back through every source that ties (D through the
/// predecessors' A one row up when it extends an aligned path; E through each tying
/// predecessor's A by an opening or its E by an extension; F through the cell above's A by
/// an opening or its F by an extension; A through whichever of D, E and F it equals), to
/// any cell on the edge ([`interior`]). The reference programme carries the same flag
/// through every cell instead.
fn touches_edge<M: ErrorModel>(m: &M, read: &[u8], g: &Graph, band: &GraphBand, ws: &mut Workspace) -> bool {
    let n = read.len() as u32;
    let (open, ext) = (m.gap_open(), m.gap_extend());
    const D: u8 = 1;
    const E: u8 = 2;
    const F: u8 = 4;
    const A: u8 = 8;
    let mut marks = std::mem::take(&mut ws.marks);
    let mut stack = std::mem::take(&mut ws.stack);
    let mut interiors = std::mem::take(&mut ws.interiors);
    ws.mark_gen += 1;
    if ws.mark_gen >= 1 << 27 {
        marks.fill(0);
        ws.mark_gen = 1;
    }
    if marks.len() < band.cells {
        marks.resize(band.cells, 0);
    }
    let gen = ws.mark_gen;
    stack.clear();
    stack.extend(ws.ends.iter().map(|&(i, v)| (i, v, D)));
    // per column, its cells off the edge, once asked (u32::MAX: not yet)
    interiors.clear();
    interiors.resize(band.v1 - band.v0, (u32::MAX, 0));
    let mut found = false;
    while let Some((i, v, kind)) = stack.pop() {
        let vu = v as usize;
        let (lo, hi) = (band.lo[vu - band.v0], band.hi[vu - band.v0]);
        let idx = (band.off[vu - band.v0] + i - lo) as usize;
        let held = if marks[idx] >> 4 == gen { marks[idx] & 15 } else { 0 };
        if held & u32::from(kind) != 0 {
            continue;
        }
        marks[idx] = gen << 4 | held | u32::from(kind);
        let slot = &mut interiors[vu - band.v0];
        if slot.0 == u32::MAX {
            *slot = interior(g, band, vu, lo, hi, n);
        }
        let (i0, i1) = *slot;
        if i < i0 || i > i1 {
            found = true;
            break;
        }
        let c = sources(m, read, g, band, ws, vu, i);
        match kind {
            D => {
                if c.from_a {
                    for &u in g.preds(vu) {
                        if band.cell(u as usize, i - 1).is_some_and(|x| ws.a[x] == c.diag) && c.diag > NEG / 2 {
                            stack.push((i - 1, u, A));
                        }
                    }
                }
            }
            E => {
                for &u in g.preds(vu) {
                    if let Some(x) = band.cell(u as usize, i) {
                        let (eo, ee) = (ws.a[x] - open, ws.e[x] - ext);
                        let val = eo.max(ee);
                        if val == c.en {
                            if eo == val {
                                stack.push((i, u, A));
                            }
                            if ee == val {
                                stack.push((i, u, E));
                            }
                        }
                    }
                }
            }
            F => {
                if i > lo {
                    if c.fo == c.f {
                        stack.push((i - 1, v, A));
                    }
                    if c.fe == c.f {
                        stack.push((i - 1, v, F));
                    }
                }
            }
            _ => {
                if c.d == c.a {
                    stack.push((i, v, D));
                }
                if c.en == c.a {
                    stack.push((i, v, E));
                }
                if c.f == c.a {
                    stack.push((i, v, F));
                }
            }
        }
    }
    ws.marks = marks;
    ws.stack = stack;
    ws.interiors = interiors;
    found
}

/// The programme as first written, cell by cell: the reference [`run`] is checked against
/// (every cell, every bit, the score and its end), and the programme of a graph whose
/// layout is not topological. `TRACK`: follow whether paths touch the band's edge.
/// `TRACE`: keep the directions for [`traceback`].
///
/// Per cell (row i, base v): D, a path whose last move aligned read base i - 1 to v
/// (after an aligned path at a predecessor, or starting there after the read bases before
/// it left over); E, a deletion (v skipped after an aligned base); F, an insertion (read
/// base i - 1 after an aligned base at v); A the best of them. A path ends on an aligned
/// base, the read bases after it left over.
fn run_reference<M: ErrorModel, const TRACK: bool, const TRACE: bool>(
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

/// The lowest score a band of `slack` answers for by its edge alone, for `read` (as the linear
/// scorer's: [`score::stand_floor`], the read's free shift along a repeat counted).
fn floor_for(slack: i64, read: &GraphRead) -> i32 {
    score::stand_floor(slack, read.free_shift)
}

/// Whether one of the read's exact k-mers (at its start, middle and end, unmasked) lies
/// along a path of the graph ending where the band does not hold it: a copy of the read's
/// bases elsewhere, which an alignment reaches with no indel. Shift-And over the DAG.
fn copy_off_band(read: &[u8], g: &Graph, band: &GraphBand) -> bool {
    if g.topological {
        copy_off_band_lanes(read, g, band)
    } else {
        copy_off_band_reference(read, g, band)
    }
}

/// [`copy_off_band`] with the three k-mers run at once, one 21-bit lane of a word each
/// (Shift-And: a lane's bit j says the k-mer's first j + 1 bases end at the base; bit 20
/// takes the shift out of bit 19 and is cleared, so no lane reaches the next). A base's
/// lanes are its predecessors' in reach ORed, shifted a base, each lane's first bit set,
/// and kept where its base matches; a lane's last bit at a base whose row past that k-mer
/// the band does not hold is a copy off the band. The same bases and rows as one k-mer at
/// a time.
fn copy_off_band_lanes(read: &[u8], g: &Graph, band: &GraphBand) -> bool {
    const LANE: u32 = 21;
    let n = read.len();
    if n < COPY_K {
        return false;
    }
    let mut offsets = vec![0, (n - COPY_K) / 2, n - COPY_K];
    offsets.dedup();
    // per lane: its k-mer's row past its end, and its bits (none for a k-mer with a masked base)
    let (mut firsts, mut tops, mut pick, mut rows) = (0u64, [0u64; 3], [0u64; 4], [u32::MAX; 3]);
    let mut any = false;
    for (l, &o) in offsets.iter().enumerate() {
        let kmer = &read[o..o + COPY_K];
        if kmer.contains(&WILD) {
            continue;
        }
        any = true;
        let shift = LANE * l as u32;
        firsts |= 1u64 << shift;
        tops[l] = 1u64 << (shift + COPY_K as u32 - 1);
        rows[l] = (o + COPY_K) as u32;
        for (j, &b) in kmer.iter().enumerate() {
            let c = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                _ => 3,
            };
            pick[c] |= 1u64 << (shift + j as u32);
        }
    }
    if !any {
        return false;
    }
    // a masked graph base matches every base of every lane
    let lanes = firsts * ((1u64 << COPY_K) - 1);
    let top_any = tops[0] | tops[1] | tops[2];
    let mut state = vec![0u64; band.v1 - band.v0];
    let reach: Vec<bool> = (band.v0..band.v1).map(|v| band.reached(g, v)).collect();
    for v in band.v0..band.v1 {
        let k = v - band.v0;
        if !reach[k] {
            state[k] = 0;
            continue;
        }
        let mut from = 0u64;
        for &u in g.preds(v) {
            let u = u as usize;
            if u >= band.v0 && u < band.v1 && reach[u - band.v0] {
                from |= state[u - band.v0];
            }
        }
        let mask = match g.bases[v] {
            b'A' => pick[0],
            b'C' => pick[1],
            b'G' => pick[2],
            b'T' => pick[3],
            _ => lanes,
        };
        let st = ((from << 1) | firsts) & mask & lanes;
        state[k] = st;
        if st & top_any != 0 && (0..3).any(|l| st & tops[l] != 0 && !band.holds(v, rows[l])) {
            return true;
        }
    }
    false
}

/// [`copy_off_band`] one k-mer at a time, as first written (and for a graph not laid out
/// in order, whose predecessors' states it reads as they stand).
fn copy_off_band_reference(read: &[u8], g: &Graph, band: &GraphBand) -> bool {
    let n = read.len();
    if n < COPY_K {
        return false;
    }
    let mut offsets = vec![0, (n - COPY_K) / 2, n - COPY_K];
    offsets.dedup();
    let top = 1u32 << (COPY_K - 1);
    let mut state = vec![0u32; band.v1 - band.v0];
    // which bases are in the read's reach, once for the three k-mers
    let reach: Vec<bool> = (band.v0..band.v1).map(|v| band.reached(g, v)).collect();
    let reached = |u: usize| u >= band.v0 && u < band.v1 && reach[u - band.v0];
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
        let masks = (*b"ACGT").map(mask);
        let pick = |b: u8| match b {
            b'A' => masks[0],
            b'C' => masks[1],
            b'G' => masks[2],
            b'T' => masks[3],
            _ => u32::MAX,
        };
        let row = (o + COPY_K) as u32;
        for v in band.v0..band.v1 {
            if !reach[v - band.v0] {
                state[v - band.v0] = 0;
                continue;
            }
            let from = g.preds(v).iter().filter(|&&u| reached(u as usize)).fold(0u32, |s, &u| s | state[u as usize - band.v0]);
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
    /// The longest tandem-repeat period in its reach: how far its own placement may sit off its
    /// best alignment at no cost, which lowers what a band answers for ([`score::stand_floor`]).
    pub(crate) free_shift: i64,
    /// Measurement only (`GBCMS_JUDGE_BAND=full`, or `longinsall` where a state skips its bands):
    /// no band; the full programme over its reach.
    pub(crate) full: bool,
    /// Measurement only (`GBCMS_JUDGE_BAND=longfull`, `longins` or `longinsall`): what a long
    /// listed length change inside the read's reach but outside the band's widening window does
    /// ([`LongRule`]).
    pub(crate) long: LongRule,
}

/// Measurement only (#369): a listed long insertion inside a read's reach (its span and its
/// length each side) but outside the window its band widens for (the slack and two bases each
/// side) can hold the read's best path, which the band then misses: a read whose body lies in a
/// copy of its own sequence inside the insertion (an expanded repeat, a duplication). All of an
/// insertion's bases stand at its anchor's coordinate, so a near one widens the band already, and
/// a path to a far deletion pays the gap to it and gains nothing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) enum LongRule {
    /// As before: the band ignores it.
    #[default]
    Off,
    /// `longfull` (option 2 as first measured): a state with a length change past the band's
    /// window longer than its slack, an insertion or a deletion, goes straight to the full programme.
    Full,
    /// `longins` (option 2 on insertions): a state with such an insertion goes straight to the full
    /// programme.
    Ins,
    /// `longinsall` (the whole read): where any state of the site has one, every state of the
    /// read runs the full programme ([`GraphRead::full`], set by the site).
    InsAll,
}

/// Half-width of the wider band tried before the full programme (past the length changes).
const WIDE_SLACK: i64 = 4 * BAND_SLACK;

/// The read's reach in the graph (its reference span widened by its length and by the
/// length changes of the rows near it: a deletion's far side lies that much further along
/// the reference) and the half-widths of the band and the wider band.
fn bands(g: &Graph, read: &GraphRead) -> ((i64, i64), [i64; 2]) {
    let l = read.seq.len() as i64;
    let shifts = near_shift(g, read);
    let trim = (read.span.0 - l - shifts, read.span.1 + l + shifts);
    (trim, [BAND_SLACK + shifts, WIDE_SLACK + shifts])
}

/// The length changes the read's band widens for: the rows the read can pass, those over its
/// span, widened by what they shift it (the slack and two bases, and twice the shift, each side).
fn near_shift(g: &Graph, read: &GraphRead) -> i64 {
    let mut shifts = 0;
    loop {
        let more = g.shifts_near(read.span.0, read.span.1, BAND_SLACK + 2 + 2 * shifts);
        if more == shifts {
            return shifts;
        }
        shifts = more;
    }
}

/// Measurement only (#369): the rows' length changes in the read's reach ([`bands`]) but
/// outside the window its band widens for ([`near_shift`]): per row its REF span and its change
/// (positive: an insertion; negative: a deletion).
pub(crate) fn far_changes(g: &Graph, read: &GraphRead) -> Vec<(i64, i64, i64)> {
    if g.shifts.is_empty() {
        return Vec::new();
    }
    let l = read.seq.len() as i64;
    let shifts = near_shift(g, read);
    let by = BAND_SLACK + 2 + 2 * shifts;
    let trim = (read.span.0 - l - shifts, read.span.1 + l + shifts);
    let near = |s: i64, e: i64| s < read.span.1 + by && e > read.span.0 - by;
    g.shifts
        .iter()
        .zip(&g.grows)
        .filter(|(&(s, e, _), _)| s < trim.1 && e > trim.0 && !near(s, e))
        .map(|(&(s, e, d), &grows)| (s, e, if grows { d } else { -d }))
        .collect()
}

/// Measurement only (#369): whether the read goes straight to the full programme in `g` under
/// its [`LongRule`]: a length change past its band's window longer than the band's slack (an
/// insertion; under [`LongRule::Full`] a deletion too).
pub(crate) fn skips_bands(g: &Graph, read: &GraphRead) -> bool {
    if read.long == LongRule::Off || g.shifts.is_empty() {
        return false;
    }
    let deletions = read.long == LongRule::Full;
    far_changes(g, read).iter().any(|&(_, _, d)| d > BAND_SLACK || (deletions && -d > BAND_SLACK))
}

/// The read's best score in `g`: in the band where it stands, else the wider band where
/// it stands, else the full programme over the read's reach.
pub(crate) fn score<M: ErrorModel>(m: &M, read: &GraphRead, g: &Graph, ws: &mut Workspace) -> i32 {
    decided(m, read, g, ws).0
}

/// The banded programme, its edge flag searched for only where the band could stand on it
/// (the score at least `floor`): below it the band falls to the next whatever the flag.
fn run_tracked<M: ErrorModel>(m: &M, read: &[u8], g: &Graph, band: &GraphBand, ws: &mut Workspace, floor: i32) -> GraphRun {
    let r = run::<M, false>(m, read, g, band, ws);
    if r.score < floor {
        return r;
    }
    GraphRun { clipped: touches_edge(m, read, g, band, ws), ..r }
}

/// [`decided`] for each of a site's state graphs at once (`graphs[k]`, its workspace
/// `wss[k]`, which keeps its deciding run for [`kept_path`]): every state's narrow band,
/// then the wider band of those it did not decide, then the full programme of the rest. A
/// graph laid out as an earlier one (`shares[k] = Some((s, p))`: state `s`'s graph, its
/// bases the same before base `p`, [`Graph::layout_until`]) takes the columns before the
/// first such base from state `s`'s run over the same band at that step, where `s` ran it
/// (the same band: bands are fixed by the layout). Each state's score and band are its
/// [`decided`]'s.
pub(crate) fn decided_states<M: ErrorModel>(m: &M, read: &GraphRead, graphs: &[&Graph], shares: &[Option<(usize, usize)>], wss: &mut [Workspace]) -> Vec<(i32, Option<i64>)> {
    let k = graphs.len();
    let mut out: Vec<Option<(i32, Option<i64>)>> = vec![None; k];
    for ws in wss.iter_mut() {
        ws.last = None;
    }
    if read.seq.is_empty() {
        return vec![(0, None); k];
    }
    let planned: Vec<((i64, i64), [i64; 2])> = graphs.iter().map(|g| bands(g, read)).collect();
    for step in (if read.full { 2 } else { 0 })..3 {
        // a state that skips the bands (measurement only, option 2) runs the full programme alone
        let runs: Vec<bool> = out.iter().enumerate().map(|(i, o)| o.is_none() && !(step < 2 && skips_bands(graphs[i], read))).collect();
        // this step's band per state; a state laid out as the one it shares from has its band
        let mut stepped: Vec<Option<std::rc::Rc<GraphBand>>> = vec![None; k];
        // a state shares at this step from one that runs it too; that state keeps its total
        // and ends before the first base any of them takes over
        let from: Vec<Option<usize>> = (0..k).map(|i| shares[i].filter(|&(s, _)| runs[i] && runs[s]).map(|(s, _)| s)).collect();
        let mut snap_at: Vec<Option<usize>> = vec![None; k];
        for i in 0..k {
            if let (Some(s), Some((_, p))) = (from[i], shares[i]) {
                snap_at[s] = Some(snap_at[s].map_or(p, |q: usize| q.min(p)));
            }
        }
        for i in 0..k {
            if !runs[i] {
                continue;
            }
            let g = graphs[i];
            let (trim, halves) = planned[i];
            let half = (step < 2).then(|| halves[step]);
            let band = match from[i].and_then(|s| stepped[s].clone()) {
                Some(b) => b,
                None => std::rc::Rc::new(GraphBand::new(g, &read.rp, half, trim)),
            };
            stepped[i] = Some(band.clone());
            let (head, rest) = wss.split_at_mut(i);
            let ws = &mut rest[0];
            let start = from[i].map(|s| (&head[s], snap_at[s].unwrap_or(0)));
            let r = if !g.topological {
                if step < 2 {
                    run_reference::<M, true, false>(m, read.seq, g, &band, ws)
                } else {
                    run_reference::<M, false, false>(m, read.seq, g, &band, ws)
                }
            } else {
                run_from::<M, false>(m, read.seq, g, &band, ws, start, snap_at[i])
            };
            let stands = match step {
                0 | 1 => {
                    let floor = floor_for([BAND_SLACK, WIDE_SLACK][step], read);
                    let clipped = if g.topological { r.score >= floor && touches_edge(m, read.seq, g, &band, ws) } else { r.clipped };
                    !clipped && r.score >= floor && (r.score == 0 || !copy_off_band(read.seq, g, &band))
                }
                _ => true,
            };
            if stands {
                out[i] = Some((r.score, half));
                if g.topological {
                    ws.last = Some((band, r));
                }
            }
        }
    }
    out.into_iter().map(|o| o.expect("the full programme decides every state")).collect()
}

/// [`score`], with the band that decided it (its half-width; None for the full
/// programme), which [`path`] traces in.
pub(crate) fn decided<M: ErrorModel>(m: &M, read: &GraphRead, g: &Graph, ws: &mut Workspace) -> (i32, Option<i64>) {
    ws.last = None;
    if read.seq.is_empty() {
        return (0, None);
    }
    let (trim, halves) = bands(g, read);
    let steps = if read.full || skips_bands(g, read) { 0 } else { 2 };
    for (half, slack) in halves.into_iter().zip([BAND_SLACK, WIDE_SLACK]).take(steps) {
        let band = GraphBand::new(g, &read.rp, Some(half), trim);
        let floor = floor_for(slack, read);
        let r = if g.topological { run_tracked(m, read.seq, g, &band, ws, floor) } else { run_reference::<M, true, false>(m, read.seq, g, &band, ws) };
        if !r.clipped && r.score >= floor && (r.score == 0 || !copy_off_band(read.seq, g, &band)) {
            if g.topological {
                ws.last = Some((std::rc::Rc::new(band), r));
            }
            return (r.score, Some(half));
        }
    }
    let full = GraphBand::new(g, &read.rp, None, trim);
    let r = if g.topological { run::<M, false>(m, read.seq, g, &full, ws) } else { run_reference::<M, false, false>(m, read.seq, g, &full, ws) };
    if g.topological {
        ws.last = Some((std::rc::Rc::new(full), r));
    }
    (r.score, None)
}

/// Tests only: the full programme over the read's reach (the last of [`decided`]'s steps,
/// what its bands stand for).
#[cfg(test)]
pub(crate) fn reach_score<M: ErrorModel>(m: &M, read: &GraphRead, g: &Graph) -> i32 {
    if read.seq.is_empty() {
        return 0;
    }
    let (trim, _) = bands(g, read);
    let full = GraphBand::new(g, &read.rp, None, trim);
    let mut ws = Workspace::default();
    if g.topological {
        run::<M, false>(m, read.seq, g, &full, &mut ws).score
    } else {
        run_reference::<M, false, false>(m, read.seq, g, &full, &mut ws).score
    }
}

/// The best path of the read in a graph: its score and the rows whose alleles it takes.
#[derive(Debug, Clone)]
pub(crate) struct GraphPath {
    pub(crate) score: i32,
    pub(crate) rows: Vec<usize>,
}

/// The read's best path in `g` traced back from the values `ws` still holds from
/// [`decided`] on this read and graph (its band and run, the band `half`): the programme
/// [`path`] would run again holds the same values. Where `ws` holds none (an empty read, a
/// graph not laid out in order), [`path`] runs it.
pub(crate) fn kept_path<M: ErrorModel>(m: &M, read: &GraphRead, g: &Graph, half: Option<i64>, ws: &Workspace) -> GraphPath {
    match &ws.last {
        Some((band, r)) => traceback(m, read.seq, g, band, ws, *r),
        None => path(m, read, g, half, &mut Workspace::default()),
    }
}

/// Workspaces for `k` programmes (one per state of a site), reused across reads on this
/// thread: their arrays only grow.
pub(crate) fn with_workspaces<R>(k: usize, f: impl FnOnce(&mut [Workspace]) -> R) -> R {
    // taken out for the call (a call inside it starts a pool of its own)
    let mut pool = POOL.with(|p| std::mem::take(&mut *p.borrow_mut()));
    if pool.len() < k {
        pool.resize_with(k, Workspace::default);
    }
    let out = f(&mut pool[..k]);
    POOL.with(|p| *p.borrow_mut() = pool);
    out
}

/// The read's best path in `g`, traced back in the band [`decided`] chose (`half`).
pub(crate) fn path<M: ErrorModel>(m: &M, read: &GraphRead, g: &Graph, half: Option<i64>, ws: &mut Workspace) -> GraphPath {
    if read.seq.is_empty() {
        return GraphPath { score: 0, rows: Vec::new() };
    }
    let (trim, _) = bands(g, read);
    let band = GraphBand::new(g, &read.rp, half, trim);
    if !g.topological {
        let r = run_reference::<M, false, true>(m, read.seq, g, &band, ws);
        return traceback_reference(g, &band, ws, r);
    }
    let r = run::<M, false>(m, read.seq, g, &band, ws);
    traceback(m, read.seq, g, &band, ws, r)
}

/// The best path back from its end, each cell's move read back from the programme's arrays
/// ([`sources`]) as the reference programme's directions hold it: D extends an aligned path
/// from the first predecessor reaching the diagonal, else starts there; A is D where they
/// are equal, else E, else F; E extends a deletion where the first predecessor reaching it
/// does so by more than by opening one; F extends an insertion where that beats opening one.
fn traceback<M: ErrorModel>(m: &M, read: &[u8], g: &Graph, band: &GraphBand, ws: &Workspace, r: GraphRun) -> GraphPath {
    let (open, ext) = (m.gap_open(), m.gap_extend());
    let mut visited: Vec<u32> = Vec::new();
    let Some((mut i, mut v)) = r.end else {
        return GraphPath { score: r.score, rows: Vec::new() };
    };
    // states: 0 D, 1 A, 2 E, 3 F
    let mut state = 0u8;
    while i > 0 {
        if band.cell(v as usize, i).is_none() {
            break;
        }
        let c = sources(m, read, g, band, ws, v as usize, i);
        match state {
            0 => {
                visited.push(v);
                if !c.from_a || i == 1 {
                    break;
                }
                // the first predecessor reaching the diagonal
                let first = g.preds(v as usize).iter().copied().find(|&u| band.cell(u as usize, i - 1).is_some_and(|x| ws.a[x] == c.diag));
                let Some(u) = first else { break };
                v = u;
                i -= 1;
                state = 1;
            }
            1 => {
                state = if c.d == c.a {
                    0
                } else if c.en == c.a {
                    2
                } else {
                    3
                }
            }
            2 => {
                visited.push(v);
                // the first predecessor reaching the deletion, and whether it extends one
                let first = g.preds(v as usize).iter().copied().find_map(|u| {
                    let x = band.cell(u as usize, i)?;
                    let (eo, ee) = (ws.a[x] - open, ws.e[x] - ext);
                    (eo.max(ee) == c.en).then_some((u, ee > eo))
                });
                let Some((u, extends)) = first else { break };
                state = if extends { 2 } else { 1 };
                v = u;
            }
            _ => {
                state = if c.fe > c.fo { 3 } else { 1 };
                i -= 1;
            }
        }
    }
    let mut rows: Vec<usize> = visited.iter().filter_map(|&v| g.segment_of(v as usize).row).collect();
    rows.sort_unstable();
    rows.dedup();
    GraphPath { score: r.score, rows }
}

/// What a read's best path does with each of its bases: aligned to a graph base, inserted, or
/// left over at an end; and the graph bases it skips right after each (a deletion). The
/// read's call is examined on this alignment ([`Site::alignment_on_path`]).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PathSteps {
    pub(crate) score: i32,
    pub(crate) steps: Vec<Step>,
    pub(crate) skipped_after: Vec<u32>,
}

/// One read base's step on a path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Step {
    At(u32),
    Inserted,
    LeftOver,
}

/// [`kept_path`] as steps per read base ([`PathSteps`]): the same traceback, read base by
/// read base. None where `ws` keeps no run (an empty read, a graph not laid out in order).
pub(crate) fn kept_steps<M: ErrorModel>(m: &M, read: &GraphRead, g: &Graph, ws: &Workspace) -> Option<PathSteps> {
    let (band, r) = ws.last.as_ref()?;
    let (open, ext) = (m.gap_open(), m.gap_extend());
    let n = read.seq.len();
    let mut steps = vec![Step::LeftOver; n];
    let mut skipped_after = vec![0u32; n];
    if let Some((mut i, mut v)) = r.end {
        // states: 0 D, 1 A, 2 E, 3 F (as [`traceback`])
        let mut state = 0u8;
        while i > 0 {
            if band.cell(v as usize, i).is_none() {
                break;
            }
            let c = sources(m, read.seq, g, band, ws, v as usize, i);
            match state {
                0 => {
                    steps[i as usize - 1] = Step::At(v);
                    if !c.from_a || i == 1 {
                        break;
                    }
                    let first = g.preds(v as usize).iter().copied().find(|&u| band.cell(u as usize, i - 1).is_some_and(|x| ws.a[x] == c.diag));
                    let Some(u) = first else { break };
                    v = u;
                    i -= 1;
                    state = 1;
                }
                1 => {
                    state = if c.d == c.a {
                        0
                    } else if c.en == c.a {
                        2
                    } else {
                        3
                    }
                }
                2 => {
                    // graph base v skipped after read base i - 1
                    skipped_after[i as usize - 1] += 1;
                    let first = g.preds(v as usize).iter().copied().find_map(|u| {
                        let x = band.cell(u as usize, i)?;
                        let (eo, ee) = (ws.a[x] - open, ws.e[x] - ext);
                        (eo.max(ee) == c.en).then_some((u, ee > eo))
                    });
                    let Some((u, extends)) = first else { break };
                    state = if extends { 2 } else { 1 };
                    v = u;
                }
                _ => {
                    steps[i as usize - 1] = Step::Inserted;
                    state = if c.fe > c.fo { 3 } else { 1 };
                    i -= 1;
                }
            }
        }
    }
    Some(PathSteps { score: r.score, steps, skipped_after })
}

/// The traceback over the reference programme's directions.
fn traceback_reference(g: &Graph, band: &GraphBand, ws: &Workspace, r: GraphRun) -> GraphPath {
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
    run::<M, false>(m, read, g, &full, &mut Workspace::default()).score
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_threads_work_space_counts_the_bands_its_workspaces_keep() {
        // A read's programmes leave each state's deciding band in its workspace (`last`), its
        // path traced back from it; the bands stay in the thread's pool with the arrays, and
        // the work space the budget counts holds them, a band states laid out alike share once.
        use super::super::haplotypes::tests::{random_seq, site};
        use super::super::haplotypes::{ReadView, Row};
        use rust_htslib::bam::record::Cigar;
        free_pool();
        let c = random_seq(600, 51);
        let x = Row::new(300, &c[300..301], if c[300] == b'A' { b"C" } else { b"A" });
        let s = site(&c, &x, &[], 100);
        let read = ReadView::new(&c[250..350], &[30; 100], 250, &[Cigar::Match(100)], (0, 100), 20).unwrap();
        let gr = s.graph_read(&read);
        let graphs: Vec<&Graph> = s.graphs.iter().map(|g| &g.graph).collect();
        let shares: Vec<Option<(usize, usize)>> = s.graphs.iter().map(|g| g.share).collect();
        with_workspaces(graphs.len(), |wss| decided_states(&ShortReads, &gr, &graphs, &shares, wss));
        let (arrays, bands) = POOL.with(|p| {
            let pool = p.borrow();
            let mut seen: Vec<*const GraphBand> = Vec::new();
            let mut bands = 0;
            for ws in pool.iter() {
                if let Some((b, _)) = &ws.last {
                    if !seen.contains(&std::rc::Rc::as_ptr(b)) {
                        seen.push(std::rc::Rc::as_ptr(b));
                        bands += (b.lo.capacity() + b.hi.capacity() + b.off.capacity()) * 4 + std::mem::size_of::<GraphBand>();
                    }
                }
            }
            (pool.iter().map(Workspace::heap_bytes).sum::<usize>(), bands)
        });
        assert!(bands > 0, "the workspaces keep their bands");
        assert_eq!(pool_bytes(), arrays + bands, "the work space counted, against the arrays and bands the pool holds");
    }

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

    /// Deterministic pseudo-random numbers (xorshift), `0..n`.
    fn roll(state: &mut u64, n: u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state % n.max(1)
    }

    fn base(st: &mut u64) -> u8 {
        if roll(st, 20) == 0 {
            WILD
        } else {
            b"ACGT"[roll(st, 4) as usize]
        }
    }

    /// A random sequence graph as a site lays one out: a reference of up to 90 bases (some
    /// masked) cut where the rows' spans start and end, a segment per row's allele (an SNV,
    /// an MNV, an insertion, a deletion to an anchor, an empty allele, rows overlapping one
    /// another), each joined to the segments ending where it starts and starting where it
    /// ends; a repeat now and then, so copies of a read's bases lie off its band.
    fn random_graph(st: &mut u64) -> Graph {
        random_graph_of(st, 90, 7)
    }

    /// [`random_graph`] up to `max_len` bases and `max_rows` rows.
    fn random_graph_of(st: &mut u64, max_len: u64, max_rows: u64) -> Graph {
        let len = 4 + roll(st, max_len - 3) as i64;
        let unit: Vec<u8> = (0..1 + roll(st, 4)).map(|_| base(st)).collect();
        let repeat = roll(st, 3) == 0;
        let reference: Vec<u8> = (0..len).map(|i| if repeat { unit[i as usize % unit.len()] } else { base(st) }).collect();
        let mut rows: Vec<((i64, i64), Vec<u8>)> = Vec::new();
        for _ in 0..roll(st, max_rows) {
            let s = roll(st, len as u64) as i64;
            let e = (s + 1 + roll(st, 4) as i64).min(len);
            let alt: Vec<u8> = (0..roll(st, 7)).map(|_| base(st)).collect();
            rows.push(((s, e), alt));
        }
        let mut cuts: Vec<i64> = vec![0, len];
        for ((s, e), _) in &rows {
            cuts.push(*s);
            cuts.push(*e);
        }
        cuts.sort_unstable();
        cuts.dedup();
        let mut b = GraphBuilder::default();
        let mut starting: Vec<(i64, usize)> = Vec::new();
        let mut ending: Vec<(i64, usize)> = Vec::new();
        for w in cuts.windows(2) {
            let id = b.segment(None, (w[0], w[1]), reference[w[0] as usize..w[1] as usize].to_vec());
            starting.push((w[0], id));
            ending.push((w[1], id));
        }
        for (r, ((s, e), alt)) in rows.iter().enumerate() {
            let id = b.segment(Some(r), (*s, *e), alt.clone());
            starting.push((*s, id));
            ending.push((*e, id));
        }
        for &(p, x) in &ending {
            for &(q, y) in &starting {
                if p == q && p < len && x != y {
                    b.edge(x, y);
                }
            }
        }
        b.finish()
    }

    /// A random read over `g`: a path's bases from a random start, mutated (substitutions,
    /// indels, masked bases, a junk end), and each base's aligned position (the path's
    /// coordinates, shifted at random, made non-decreasing).
    fn random_read(g: &Graph, st: &mut u64) -> (Vec<u8>, Vec<i64>) {
        random_read_of(g, st, 60)
    }

    /// [`random_read`] of up to `max_take` path bases.
    fn random_read_of(g: &Graph, st: &mut u64, max_take: u64) -> (Vec<u8>, Vec<i64>) {
        let mut path: Vec<usize> = Vec::new();
        let starts: Vec<usize> = (0..g.len()).filter(|&v| g.preds(v).is_empty()).collect();
        let mut v = if starts.is_empty() || roll(st, 4) == 0 { roll(st, g.len() as u64) as usize } else { starts[roll(st, starts.len() as u64) as usize] };
        loop {
            path.push(v);
            let next = g.succs(v);
            if next.is_empty() || path.len() > 2000 {
                break;
            }
            v = next[roll(st, next.len() as u64) as usize] as usize;
        }
        let from = roll(st, path.len() as u64) as usize;
        let take = (1 + roll(st, max_take) as usize).min(path.len() - from);
        let (mut seq, mut rp): (Vec<u8>, Vec<i64>) = (Vec::new(), Vec::new());
        let shift = roll(st, 9) as i64 - 4;
        for &v in &path[from..from + take] {
            match roll(st, 30) {
                0 => {}
                1 => {
                    seq.push(base(st));
                    rp.push(g.coord[v] + shift);
                    seq.push(g.bases[v]);
                    rp.push(g.coord[v] + shift);
                }
                2 | 3 => {
                    seq.push(b"ACGT"[roll(st, 4) as usize]);
                    rp.push(g.coord[v] + shift);
                }
                4 => {
                    seq.push(WILD);
                    rp.push(g.coord[v] + shift);
                }
                _ => {
                    seq.push(if g.bases[v] == WILD { b"ACGT"[roll(st, 4) as usize] } else { g.bases[v] });
                    rp.push(g.coord[v] + shift);
                }
            }
        }
        if roll(st, 5) == 0 {
            for _ in 0..1 + roll(st, 8) {
                seq.push(b"ACGT"[roll(st, 4) as usize]);
                rp.push(rp.last().copied().unwrap_or(0) + 1);
            }
        }
        if seq.is_empty() {
            seq.push(b'A');
            rp.push(0);
        }
        if roll(st, 6) == 0 {
            let at = roll(st, 100) as i64 - 10;
            rp = (0..seq.len() as i64).map(|i| at + i).collect();
        }
        let mut last = i64::MIN;
        for p in rp.iter_mut() {
            last = last.max(*p);
            *p = last;
        }
        (seq, rp)
    }

    /// Both programmes over one band, every output compared: the score, its end and the edge
    /// flag (searched for, and carried through every cell), A and E per cell, and the best
    /// path traced back (read back from the values, and from the directions).
    fn same_run(read: &[u8], g: &Graph, band: &GraphBand, ws: &mut Workspace, case: u64) {
        let mut w0 = Workspace::default();
        let r0 = run_reference::<ShortReads, true, true>(&ShortReads, read, g, band, &mut w0);
        let r1 = run::<ShortReads, true>(&ShortReads, read, g, band, ws);
        let c = band.cells;
        assert_eq!((r0.score, r0.end, r0.clipped), (r1.score, r1.end, r1.clipped), "case {case}");
        assert_eq!(w0.a[..c], ws.a[..c], "A, case {case}");
        assert_eq!(w0.e[..c], ws.e[..c], "E, case {case}");
        let (p0, p1) = (traceback_reference(g, band, &w0, r0), traceback(&ShortReads, read, g, band, ws, r1));
        assert_eq!((p0.score, &p0.rows), (p1.score, &p1.rows), "path, case {case}");
        // the plain programme: the same values and end, no flag
        let r2 = run::<ShortReads, false>(&ShortReads, read, g, band, ws);
        assert_eq!((r0.score, r0.end, false), (r2.score, r2.end, r2.clipped), "plain, case {case}");
    }

    /// [`decided`] and [`path`] on the reference programme.
    fn decided_reference(read: &GraphRead, g: &Graph) -> (i32, Option<i64>, GraphPath) {
        let mut ws = Workspace::default();
        let (trim, halves) = bands(g, read);
        let steps = if skips_bands(g, read) { 0 } else { 2 };
        for (half, slack) in halves.into_iter().zip([BAND_SLACK, WIDE_SLACK]).take(steps) {
            let band = GraphBand::new(g, &read.rp, Some(half), trim);
            let r = run_reference::<ShortReads, true, true>(&ShortReads, read.seq, g, &band, &mut ws);
            if !r.clipped && r.score >= floor_for(slack, read) && (r.score == 0 || !copy_off_band(read.seq, g, &band)) {
                return (r.score, Some(half), traceback_reference(g, &band, &ws, r));
            }
        }
        let full = GraphBand::new(g, &read.rp, None, trim);
        let r = run_reference::<ShortReads, false, true>(&ShortReads, read.seq, g, &full, &mut ws);
        (r.score, None, traceback_reference(g, &full, &ws, r))
    }

    #[test]
    fn the_column_programme_is_the_cell_programme_on_random_graphs_reads_and_bands() {
        let mut st = 0x9e37_79b9_7f4a_7c15u64;
        // one workspace throughout: what an earlier read left in it never reaches a score
        let mut ws = Workspace::default();
        let (mut clipped, mut copies) = (0, 0);
        for case in 0..6000u64 {
            let g = random_graph(&mut st);
            assert!(g.topological, "case {case}: the builder lays a graph out in order");
            let (read, rp) = random_read(&g, &mut st);
            let half = match roll(&mut st, 5) {
                0 => None,
                _ => Some(roll(&mut st, 24) as i64),
            };
            let (lo, hi) = (rp[0], *rp.last().unwrap());
            let trim = match roll(&mut st, 4) {
                0 => (i64::MIN / 4, i64::MAX / 4),
                _ => (lo - roll(&mut st, 70) as i64, hi + 1 + roll(&mut st, 70) as i64),
            };
            let band = GraphBand::new(&g, &rp, half, trim);
            let b0 = GraphBand::new_reference(&g, &rp, half, trim);
            assert_eq!((&b0.lo, &b0.hi, &b0.off, b0.cells, b0.v0, b0.v1), (&band.lo, &band.hi, &band.off, band.cells, band.v0, band.v1), "band, case {case}");
            // the edge rows, cell by cell
            let n = read.len() as u32;
            for v in band.v0..band.v1 {
                let (lo, hi) = (band.lo[v - band.v0], band.hi[v - band.v0]);
                if lo > hi {
                    continue;
                }
                let (i0, i1) = interior(&g, &band, v, lo, hi, n);
                for i in lo..=hi {
                    assert_eq!(edge_cell(&g, &band, v, i, lo, hi, n), i < i0 || i > i1, "case {case} base {v} row {i}");
                }
            }
            same_run(&read, &g, &band, &mut ws, case);
            clipped += usize::from(run::<ShortReads, true>(&ShortReads, &read, &g, &band, &mut ws).clipped);
            let copy = copy_off_band_reference(&read, &g, &band);
            assert_eq!(copy_off_band_lanes(&read, &g, &band), copy, "copy off the band, case {case}");
            copies += usize::from(copy);
            // and the scores, bands and paths the judge reads
            let gr = GraphRead { seq: &read, rp: rp.clone(), span: (lo, hi + 1), free_shift: 0, full: false, long: LongRule::Off };
            let (s0, h0, p0) = decided_reference(&gr, &g);
            let (s1, h1) = decided(&ShortReads, &gr, &g, &mut ws);
            assert_eq!((s0, h0), (s1, h1), "decided, case {case}");
            let kept = kept_path(&ShortReads, &gr, &g, h1, &ws);
            assert_eq!((p0.score, &p0.rows), (kept.score, &kept.rows), "kept path, case {case}");
            let p1 = path(&ShortReads, &gr, &g, h1, &mut ws);
            assert_eq!((p0.score, &p0.rows), (p1.score, &p1.rows), "path, case {case}");
        }
        assert!(clipped > 600, "the cases reach the band's edge ({clipped})");
        assert!(copies > 30, "the cases find copies off the band ({copies})");
    }

    #[test]
    fn the_column_programme_is_the_cell_programme_on_windows_and_reads_of_a_sites_size() {
        // windows of up to 700 bases with up to 40 rows, reads of up to 160 bases, bands from
        // the narrow to the full: as the judge's sites and reads are
        let mut st = 0x5151_2626_3737_4848u64;
        let mut ws = Workspace::default();
        for case in 0..150u64 {
            let g = random_graph_of(&mut st, 700, 40);
            let (read, rp) = random_read_of(&g, &mut st, 160);
            let (lo, hi) = (rp[0], *rp.last().unwrap());
            for half in [Some(16), Some(40 + roll(&mut st, 40) as i64), None] {
                let band = GraphBand::new(&g, &rp, half, (lo - 170, hi + 170));
                same_run(&read, &g, &band, &mut ws, case);
            }
            let gr = GraphRead { seq: &read, rp: rp.clone(), span: (lo, hi + 1), free_shift: 0, full: false, long: LongRule::Off };
            let (s0, h0, p0) = decided_reference(&gr, &g);
            let (s1, h1) = decided(&ShortReads, &gr, &g, &mut ws);
            assert_eq!((s0, h0), (s1, h1), "decided, case {case}");
            let kept = kept_path(&ShortReads, &gr, &g, h1, &ws);
            assert_eq!((p0.score, &p0.rows), (kept.score, &kept.rows), "kept path, case {case}");
            let p1 = path(&ShortReads, &gr, &g, h1, &mut ws);
            assert_eq!((p0.score, &p0.rows), (p1.score, &p1.rows), "path, case {case}");
        }
    }

    #[test]
    fn the_copy_check_by_lanes_is_the_check_a_k_mer_at_a_time_in_repeats() {
        // repeats with bubbles, reads from them placed a few bases off: copies off the band
        let mut st = 0x7777_1234_abcd_0001u64;
        let mut copies = 0;
        for case in 0..8000u64 {
            let g = random_graph_of(&mut st, 300, 12);
            let (mut read, rp) = random_read_of(&g, &mut st, 120);
            if roll(&mut st, 3) == 0 && read.len() > 25 {
                let k = roll(&mut st, read.len() as u64) as usize;
                read[k] = WILD;
            }
            let (lo, hi) = (rp[0], *rp.last().unwrap());
            let shift = roll(&mut st, 40) as i64 - 20;
            let rp: Vec<i64> = rp.iter().map(|p| p + shift).collect();
            let half = Some(roll(&mut st, 20) as i64);
            let band = GraphBand::new(&g, &rp, half, (lo - 150, hi + 150));
            let copy = copy_off_band_reference(&read, &g, &band);
            assert_eq!(copy_off_band_lanes(&read, &g, &band), copy, "case {case}");
            copies += usize::from(copy);
        }
        assert!(copies > 300, "the cases find copies off the band ({copies})");
    }

    #[test]
    fn states_laid_out_alike_take_their_first_columns_from_one_another_and_decide_alike() {
        // a substitution's site: graphs laid out alike whose bases differ from some base on
        // (one base, or a base and the bases after it), beside a graph laid out otherwise
        let mut st = 0x2468_ace0_1357_9bdfu64;
        let mut shared = 0;
        for case in 0..1500u64 {
            let g0 = random_graph_of(&mut st, 200, 6);
            let mut graphs: Vec<Graph> = vec![g0.clone()];
            for _ in 0..1 + roll(&mut st, 3) {
                let mut g = g0.clone();
                let p = roll(&mut st, g.len() as u64) as usize;
                for v in p..(p + 1 + roll(&mut st, 3) as usize).min(g.len()) {
                    g.bases[v] = b"ACGT"[roll(&mut st, 4) as usize];
                }
                graphs.push(g);
            }
            if roll(&mut st, 2) == 0 {
                graphs.insert(roll(&mut st, graphs.len() as u64 + 1) as usize, random_graph_of(&mut st, 200, 6));
            }
            let mut shares: Vec<Option<(usize, usize)>> = Vec::new();
            for k in 0..graphs.len() {
                let share = (0..k).filter(|&s| shares[s].is_none()).find_map(|s| graphs[s].layout_until(&graphs[k]).map(|p| (s, p)));
                shares.push(share);
            }
            shared += shares.iter().flatten().count();
            let base = &graphs[shares.iter().position(|s| s.is_some()).map_or(0, |k| shares[k].unwrap().0)];
            let (read, rp) = random_read_of(base, &mut st, 100);
            let (lo, hi) = (rp[0], *rp.last().unwrap());
            let gr = GraphRead { seq: &read, rp: rp.clone(), span: (lo, hi + 1), free_shift: 0, full: false, long: LongRule::Off };
            let refs: Vec<&Graph> = graphs.iter().collect();
            let mut wss: Vec<Workspace> = (0..graphs.len()).map(|_| Workspace::default()).collect();
            let all = decided_states(&ShortReads, &gr, &refs, &shares, &mut wss);
            for (k, g) in graphs.iter().enumerate() {
                let mut ws = Workspace::default();
                let one = decided(&ShortReads, &gr, g, &mut ws);
                assert_eq!(all[k], one, "case {case} state {k}: score and band");
                let (p0, p1) = (kept_path(&ShortReads, &gr, g, one.1, &ws), kept_path(&ShortReads, &gr, g, all[k].1, &wss[k]));
                assert_eq!((p0.score, p0.rows), (p1.score, p1.rows), "case {case} state {k}: path");
            }
        }
        assert!(shared > 1500, "the cases share layouts ({shared})");
    }

    #[test]
    fn a_graph_laid_out_out_of_order_is_scored_by_the_reference_programme() {
        // A base whose predecessor lies after it (none the builder makes from a site): the
        // column programme needs every predecessor first, so the reference one runs.
        let mut b = GraphBuilder::default();
        let x = b.segment(None, (0, 3), b"ACG".to_vec());
        let y = b.segment(None, (3, 6), b"TTA".to_vec());
        b.edge(x, y);
        let mut g = b.finish();
        assert!(g.topological);
        g.pred_off = vec![0, 1, 2, 3, 4, 5, 6];
        g.preds = vec![4, 0, 1, 2, 3, 4];
        g.topological = false;
        let gr = GraphRead { seq: b"ACGTA", rp: (0..5).collect(), span: (0, 5), free_shift: 0, full: false, long: LongRule::Off };
        let mut ws = Workspace::default();
        let (s0, h0, p0) = decided_reference(&gr, &g);
        let (s1, h1) = decided(&ShortReads, &gr, &g, &mut ws);
        assert_eq!((s0, h0), (s1, h1));
        let p1 = path(&ShortReads, &gr, &g, h1, &mut ws);
        assert_eq!((p0.score, p0.rows), (p1.score, p1.rows));
    }
}
