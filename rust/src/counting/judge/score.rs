//! The judge's scorer: a read against one haplotype, semi-global and affine. The read is
//! consumed whole (every readable base), the haplotype's ends are free, so a read may lie
//! anywhere along it. Match 0, a confident mismatch −4, a gap of k bases −(6 + k) (the
//! BWA-MEM scale); a read base below `--min-baseq`, or N, costs nothing against any
//! haplotype base: the one quality gate every backend shares. Read bases left over at
//! either end of the alignment (an overhang past the haplotype's free end, a junk clip)
//! cost a mismatch each, masked ones nothing: an end gap is no cheaper than mismatching
//! the bases, so it cannot take a telling base along cheaply.
//!
//! The dynamic programme runs in a band around the read's own alignment mapped into the
//! haplotype (see `Band`). The full programme decides instead where the band may have
//! cut a better alignment off: where the best banded one reaches the band's edge; where
//! it scores below what an indel long enough to leave the band costs (an alignment that
//! starts near the read's own diagonal needs one to leave it); and, where it is not
//! exact, where one of the read's exact k-mers lies on the haplotype off the band (a
//! copy of its bases elsewhere, which no indel leads to). A traceback variant ([`align`])
//! gives the alignment itself, for the bases that decide a call and the haplotype
//! positions the read covers.

/// A confident mismatch.
pub(crate) const MISMATCH: i32 = -4;
/// A gap of k bases costs `GAP_OPEN + k * GAP_EXTEND`.
pub(crate) const GAP_OPEN: i32 = 6;
pub(crate) const GAP_EXTEND: i32 = 1;
/// The margin a call needs: one confident base.
pub(crate) const DELTA: i32 = 4;
/// Haplotype columns the band holds each side of the read's diagonal, past the length
/// changes the haplotype carries there.
pub(crate) const BAND_SLACK: i64 = 16;
/// The length of the read's k-mers looked for off the band.
pub(crate) const COPY_K: usize = 20;

/// Unreachable.
const NEG: i32 = i32::MIN / 4;

/// The masked read base: matches any haplotype base.
pub(crate) const WILD: u8 = b'N';

/// A read base against a haplotype base.
#[inline]
pub(crate) fn substitution_cost(read_base: u8, hap_base: u8) -> i32 {
    if read_base == WILD || hap_base == WILD || read_base == hap_base {
        0
    } else {
        MISMATCH
    }
}

/// A read base left over at an end of the alignment: a mismatch, or nothing when masked.
#[inline]
pub(crate) fn overhang_cost(read_base: u8) -> i32 {
    if read_base == WILD {
        0
    } else {
        MISMATCH
    }
}

/// Per programme row (read bases consumed, `0..=n`), the haplotype columns `[lo, hi]`
/// the band holds. `lo` rises by at most one column a row and `hi` never falls, so an
/// alignment can leave the band only through a row's first or last column.
#[derive(Debug, Clone)]
pub(crate) struct Band {
    lo: Vec<usize>,
    hi: Vec<usize>,
    /// Whether the haplotype's ends cut the band somewhere (its edge there is the
    /// haplotype's own, so the band is not the same relative to another trim of it).
    pub(crate) cut: bool,
}

impl Band {
    /// The band around `centers` (per row, the column the read's own alignment puts it
    /// at), `half` columns each side, on a haplotype of `m` bases.
    pub(crate) fn around(centers: &[i64], half: i64, m: usize) -> Band {
        let mut cut = false;
        let mut clamp = |c: i64| {
            cut |= c < 0 || c > m as i64;
            c.clamp(0, m as i64) as usize
        };
        let mut lo: Vec<usize> = Vec::with_capacity(centers.len());
        let mut hi: Vec<usize> = Vec::with_capacity(centers.len());
        for (i, &c) in centers.iter().enumerate() {
            let (a, b) = (clamp(c - half), clamp(c + half));
            if i == 0 {
                lo.push(a);
                hi.push(b);
            } else {
                lo.push(a.min(lo[i - 1] + 1));
                hi.push(b.max(hi[i - 1]));
            }
        }
        Band { lo, hi, cut }
    }

    /// Every column, every row: the full programme.
    fn full(n: usize, m: usize) -> Band {
        Band { lo: vec![0; n + 1], hi: vec![m; n + 1], cut: true }
    }

    fn holds(&self, i: usize, j: usize) -> bool {
        j >= self.lo[i] && j <= self.hi[i]
    }

    /// Where each row's cells start in a traceback held over the band alone (row `i`'s
    /// column `j` at `starts[i] + j - lo[i]`), and how many cells it holds.
    fn trace_layout(&self) -> (Vec<usize>, usize) {
        let mut starts = Vec::with_capacity(self.lo.len());
        let mut cells = 0usize;
        for (&lo, &hi) in self.lo.iter().zip(&self.hi) {
            starts.push(cells);
            cells += (hi + 1).saturating_sub(lo);
        }
        (starts, cells)
    }
}

/// The best alignment's score, where it ends (the row and column of its last aligned
/// read base; None for a read left over whole), and whether its path (or one as good)
/// touched the band's edge.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Run {
    pub(crate) score: i32,
    end: Option<(usize, usize)>,
    clipped: bool,
}

/// Traceback directions per cell: bits 0-1 where the cell's best (aligned) path came
/// from (0 the diagonal, 1 a deletion, 2 an insertion), bit 2 a deletion extended, bit 3
/// an insertion extended, bit 4 the diagonal started there (the read bases before it
/// left over).
const FROM_E: u8 = 1;
const FROM_F: u8 = 2;
const E_EXT: u8 = 4;
const F_EXT: u8 = 8;
const D_START: u8 = 16;

/// The programme over `band`. `TRACK`: follow whether paths touch the band's edge.
/// `TRACE`: fill `trace`, the directions of the band's cells laid out by `starts`
/// ([`Band::trace_layout`]).
///
/// Per cell: D, a path whose last move aligned read base i - 1 to haplotype base j - 1
/// (from an aligned path, or starting there after the read bases before it left over);
/// E and F, a deletion or an insertion after an aligned base; A, the best of them. A
/// path ends on an aligned base, the read bases after it left over. The previous row's
/// cells the next one reads outside its band hold no path.
fn run<const TRACK: bool, const TRACE: bool>(read: &[u8], hap: &[u8], band: &Band, trace: &mut [u8], starts: &[usize]) -> Run {
    let (n, m) = (read.len(), hap.len());
    let open = GAP_OPEN + GAP_EXTEND;
    let width = m + 1;
    // Left over before read base i (prefix) and from it on (suffix).
    let mut lead = vec![0i32; n + 1];
    for i in 0..n {
        lead[i + 1] = lead[i] + overhang_cost(read[i]);
    }
    // A and F of the previous and current rows; per cell, bit 0 A touched the edge,
    // bit 1 F did.
    let (mut a_prev, mut f_prev, mut t_prev) = (vec![NEG; width], vec![NEG; width], vec![0u8; width]);
    let (mut a_cur, mut f_cur, mut t_cur) = (vec![NEG; width], vec![NEG; width], vec![0u8; width]);
    // The whole read left over: no base aligned.
    let (mut score, mut end, mut clipped) = (lead[n], None, false);
    for i in 1..=n {
        let rb = read[i - 1];
        let (lo, hi) = (band.lo[i], band.hi[i]);
        let start = lead[i - 1];
        let (mut e, mut te) = (NEG, 0u8);
        let (mut row_best, mut row_end, mut row_clipped) = (NEG, 0usize, 0u8);
        // This row's and the previous row's cells over the band (and the one before it),
        // so the loop reads them without bounds checks.
        let base = lo.saturating_sub(1);
        let (ap, fp, tp) = (&a_prev[base..=hi], &f_prev[base..=hi], &t_prev[base..=hi]);
        let (ac, fc, tc) = (&mut a_cur[base..=hi], &mut f_cur[base..=hi], &mut t_cur[base..=hi]);
        let hs = &hap[base..hi];
        for j in lo..=hi {
            let k = j - base; // this cell, in the slices; k - 1 the cell before it
            let edge = u8::from(TRACK && ((j == lo && lo > 0) || (j == hi && hi < m)));
            // F: the read base inserted after an aligned base.
            let (fo, fe) = (ap[k] - open, fp[k] - GAP_EXTEND);
            let f = fo.max(fe);
            // E: a haplotype base skipped after an aligned base.
            let first = j == lo;
            let (eo, ee) = if first { (NEG, NEG) } else { (ac[k - 1] - open, e - GAP_EXTEND) };
            let e_new = eo.max(ee);
            // D: the read base against haplotype base j - 1, after an aligned path or
            // starting there.
            let (d, from_a) = if j > 0 {
                let before = ap[k - 1];
                let sub = substitution_cost(rb, hs[k - 1]);
                if before >= start {
                    (before + sub, true)
                } else {
                    (start + sub, false)
                }
            } else {
                (NEG, false)
            };
            let a = d.max(e_new).max(f);
            if TRACK {
                let tf = edge | (u8::from(fo == f) & tp[k]) | (u8::from(fe == f) & (tp[k] >> 1));
                let teo = if first { 0 } else { tc[k - 1] & 1 };
                te = edge | (u8::from(eo == e_new) & teo) | (u8::from(ee == e_new) & te);
                let td = edge | if from_a { tp[k - 1] & 1 } else { 0 };
                let ta = (u8::from(d == a) & td) | (u8::from(e_new == a) & te) | (u8::from(f == a) & tf);
                tc[k] = (ta & 1) | ((tf & 1) << 1);
                if d > row_best {
                    (row_best, row_end, row_clipped) = (d, j, td);
                } else if d == row_best {
                    row_clipped |= td;
                }
            } else if d > row_best {
                (row_best, row_end) = (d, j);
            }
            e = e_new;
            ac[k] = a;
            fc[k] = f;
            if TRACE {
                let src = if d == a {
                    0
                } else if e_new == a {
                    FROM_E
                } else {
                    FROM_F
                };
                trace[starts[i] + j - lo] = src
                    | if ee > eo { E_EXT } else { 0 }
                    | if fe > fo { F_EXT } else { 0 }
                    | if from_a { 0 } else { D_START };
            }
        }
        if row_best > NEG / 2 {
            let total = row_best + (lead[n] - lead[i]);
            if total > score {
                (score, end, clipped) = (total, Some((i, row_end)), row_clipped != 0);
            } else if total == score {
                clipped |= row_clipped != 0;
            }
        }
        // The cells the next row reads outside this row's band hold no path.
        if i < n {
            let (nlo, nhi) = (band.lo[i + 1], band.hi[i + 1]);
            for j in nlo.saturating_sub(1)..lo.min(nhi + 1) {
                (a_cur[j], f_cur[j], t_cur[j]) = (NEG, NEG, 0);
            }
            for j in (hi + 1).max(nlo.saturating_sub(1))..=nhi.min(m) {
                (a_cur[j], f_cur[j], t_cur[j]) = (NEG, NEG, 0);
            }
        }
        std::mem::swap(&mut a_prev, &mut a_cur);
        std::mem::swap(&mut f_prev, &mut f_cur);
        std::mem::swap(&mut t_prev, &mut t_cur);
    }
    Run { score, end, clipped }
}

/// A banded programme's values, per cell of the band (row `i`'s column `j` at
/// `starts[i] + j - lo[i]`): A, E and F; and per row its best D and the first column
/// reaching it (rows `1..=n`; the read left over before them, no row 0 is computed).
#[derive(Debug, Default)]
struct Values {
    starts: Vec<usize>,
    cells: usize,
    a: Vec<i32>,
    e: Vec<i32>,
    f: Vec<i32>,
    lead: Vec<i32>,
    rows: Vec<(i32, usize)>,
}

impl Values {
    /// Cell `(i, j)`'s index where the band holds it (row 0 holds none).
    #[inline]
    fn at(&self, band: &Band, i: usize, j: usize) -> Option<usize> {
        (i > 0 && band.holds(i, j)).then(|| self.starts[i] + j - band.lo[i])
    }
}

/// [`run`] keeping every cell's values ([`Values`]) instead of directions or edge bits, and
/// no edge flag (`Run::clipped` false): the edge ([`touches_edge`]) and the traceback
/// ([`traceback_values`]) are read back from the values where needed. Every value, the
/// score and its end are [`run`]'s.
fn run_values(read: &[u8], hap: &[u8], band: &Band, vals: &mut Values) -> Run {
    let (n, m) = (read.len(), hap.len());
    let open = GAP_OPEN + GAP_EXTEND;
    let (starts, cells) = band.trace_layout();
    vals.starts = starts;
    vals.cells = cells;
    // every cell of rows 1..=n is written before it is read (row 0 is never read), so the
    // arrays only grow
    for v in [&mut vals.a, &mut vals.e, &mut vals.f] {
        if v.len() < cells {
            v.resize(cells, NEG);
        }
    }
    vals.lead.clear();
    vals.lead.push(0);
    for &b in read {
        let last = *vals.lead.last().unwrap();
        vals.lead.push(last + overhang_cost(b));
    }
    vals.rows.clear();
    vals.rows.push((NEG, 0));
    let mut pa: Vec<i32> = Vec::new();
    let mut pf: Vec<i32> = Vec::new();
    let mut score = vals.lead[n];
    let mut end = None;
    for i in 1..=n {
        let rb = read[i - 1];
        let (lo, hi) = (band.lo[i], band.hi[i]);
        let start = vals.lead[i - 1];
        // the previous row's A and F over this row's columns and the one before them
        // (none off its band, nor in row 0)
        let base = lo.saturating_sub(1);
        let w = hi + 1 - base;
        pa.clear();
        pa.resize(w, NEG);
        pf.clear();
        pf.resize(w, NEG);
        if i > 1 {
            let (plo, phi) = (band.lo[i - 1], band.hi[i - 1]);
            let (c0, c1) = (plo.max(base), phi.min(hi));
            if c0 <= c1 {
                let x0 = vals.starts[i - 1] + c0 - plo;
                let k0 = c0 - base;
                pa[k0..=k0 + (c1 - c0)].copy_from_slice(&vals.a[x0..=x0 + (c1 - c0)]);
                pf[k0..=k0 + (c1 - c0)].copy_from_slice(&vals.f[x0..=x0 + (c1 - c0)]);
            }
        }
        let x0 = vals.starts[i];
        let len = hi + 1 - lo;
        let (ac, ec, fc) = (&mut vals.a[x0..x0 + len], &mut vals.e[x0..x0 + len], &mut vals.f[x0..x0 + len]);
        let (mut row_best, mut row_end) = (NEG, 0usize);
        // Column 0 (a row whose band starts there) aligns no base: D is none.
        let skip = usize::from(lo == 0);
        if skip == 1 {
            let f = (pa[0] - open).max(pf[0] - GAP_EXTEND);
            (ac[0], fc[0]) = (NEG.max(f), f);
        }
        // From column max(lo, 1): F from A and F above (pa[k], pf[k]), D from A on the
        // diagonal (pa[k - 1]) against haplotype base j - 1, and the best of the two (A
        // before E, kept in A's place), over equal-length slices.
        let j0 = lo + skip;
        let cnt = hi + 1 - j0;
        let k0 = j0 - base;
        let (up_a, up_f, diag) = (&pa[k0..k0 + cnt], &pf[k0..k0 + cnt], &pa[k0 - 1..k0 - 1 + cnt]);
        let hs = &hap[j0 - 1..j0 - 1 + cnt];
        {
            let (ac, fc) = (&mut ac[skip..skip + cnt], &mut fc[skip..skip + cnt]);
            let cells = up_a.iter().zip(up_f).zip(diag).zip(hs).zip(ac.iter_mut().zip(fc.iter_mut()));
            for (t, ((((&ua, &uf), &dg), &h), (ao, fo))) in cells.enumerate() {
                // F: the read base inserted after an aligned base.
                let f = (ua - open).max(uf - GAP_EXTEND);
                // D: the read base against haplotype base j - 1, after an aligned path or
                // starting there (a masked base costs nothing against any).
                let sub = if rb == WILD || h == WILD || h == rb { 0 } else { MISMATCH };
                let d = dg.max(start) + sub;
                *ao = d.max(f);
                *fo = f;
                if d > row_best {
                    (row_best, row_end) = (d, j0 + t);
                }
            }
        }
        // E along the row: a haplotype base skipped after an aligned base (none in the
        // band's first column), as the best opening to its left with its extensions. A
        // deletion is never opened from one (open is at least the extension), so the best
        // A to the left, E included, is the best of D and F there: E at cell u is
        // max over k < u of (D-or-F at k) - open - (u - 1 - k) * extend.
        let mut best = NEG;
        for (u, (ao, eo)) in ac.iter_mut().zip(ec.iter_mut()).enumerate() {
            let a0 = *ao;
            let e = if u > 0 { best - open - (u as i32 - 1) * GAP_EXTEND } else { NEG };
            *ao = a0.max(e);
            *eo = e;
            best = best.max(a0 + u as i32 * GAP_EXTEND);
        }
        vals.rows.push((row_best, row_end));
        if row_best > NEG / 2 {
            let total = row_best + (vals.lead[n] - vals.lead[i]);
            if total > score {
                (score, end) = (total, Some((i, row_end)));
            }
        }
    }
    let _ = m;
    Run { score, end, clipped: false }
}

/// A cell's sources read back from [`Values`]: D and whether it extends an aligned path,
/// E's opening and extension, F's opening and extension, and A, E and F.
#[derive(Debug, Clone, Copy)]
struct CellSources {
    d: i32,
    from_a: bool,
    eo: i32,
    ee: i32,
    fo: i32,
    fe: i32,
    a: i32,
    e: i32,
    f: i32,
}

#[inline]
fn cell_sources(read: &[u8], hap: &[u8], band: &Band, vals: &Values, i: usize, j: usize, x: usize) -> CellSources {
    let open = GAP_OPEN + GAP_EXTEND;
    let up = |jj: usize| vals.at(band, i - 1, jj);
    let (d, from_a) = if j > 0 {
        let before = up(j - 1).map_or(NEG, |y| vals.a[y]);
        let start = vals.lead[i - 1];
        let sub = substitution_cost(read[i - 1], hap[j - 1]);
        if before >= start {
            (before + sub, true)
        } else {
            (start + sub, false)
        }
    } else {
        (NEG, false)
    };
    let (eo, ee) = if j == band.lo[i] { (NEG, NEG) } else { (vals.a[x - 1] - open, vals.e[x - 1] - GAP_EXTEND) };
    let (pa, pf) = up(j).map_or((NEG, NEG), |y| (vals.a[y], vals.f[y]));
    CellSources { d, from_a, eo, ee, fo: pa - open, fe: pf - GAP_EXTEND, a: vals.a[x], e: vals.e[x], f: vals.f[x] }
}

/// Whether a best end's path (or one as good) touches the band's edge, as [`run`]'s
/// `TRACK` flag: from each cell of each row scoring the best total whose D is the row's
/// best, back through every source that ties (D through the cell before it on the
/// diagonal when it extends an aligned path; E through the cell to its left, A by an
/// opening or E by an extension; F through the cell above, A by an opening or F by an
/// extension; A through whichever of D, E and F it equals), to a row's first column (past
/// column 0) or last (short of the haplotype's end).
fn touches_edge(read: &[u8], hap: &[u8], band: &Band, vals: &Values, score: i32) -> bool {
    let (n, m) = (read.len(), hap.len());
    const D: u8 = 1;
    const E: u8 = 2;
    const F: u8 = 4;
    const A: u8 = 8;
    let mut marks = vec![0u8; vals.cells];
    let mut stack: Vec<(usize, usize, u8)> = Vec::new();
    for i in 1..=n {
        let (row_best, _) = vals.rows[i];
        if row_best > NEG / 2 && row_best + (vals.lead[n] - vals.lead[i]) == score {
            for j in band.lo[i]..=band.hi[i] {
                let x = vals.starts[i] + j - band.lo[i];
                if cell_sources(read, hap, band, vals, i, j, x).d == row_best {
                    stack.push((i, j, D));
                }
            }
        }
    }
    while let Some((i, j, kind)) = stack.pop() {
        let Some(x) = vals.at(band, i, j) else { continue };
        if marks[x] & kind != 0 {
            continue;
        }
        marks[x] |= kind;
        let (lo, hi) = (band.lo[i], band.hi[i]);
        if (j == lo && lo > 0) || (j == hi && hi < m) {
            return true;
        }
        let c = cell_sources(read, hap, band, vals, i, j, x);
        match kind {
            D => {
                if c.from_a {
                    stack.push((i - 1, j - 1, A));
                }
            }
            E => {
                if j > lo {
                    if c.eo == c.e {
                        stack.push((i, j - 1, A));
                    }
                    if c.ee == c.e {
                        stack.push((i, j - 1, E));
                    }
                }
            }
            F => {
                if c.fo == c.f {
                    stack.push((i - 1, j, A));
                }
                if c.fe == c.f {
                    stack.push((i - 1, j, F));
                }
            }
            _ => {
                if c.d == c.a {
                    stack.push((i, j, D));
                }
                if c.e == c.a {
                    stack.push((i, j, E));
                }
                if c.f == c.a {
                    stack.push((i, j, F));
                }
            }
        }
    }
    false
}

/// The lowest score the band answers for by its edge alone: an alignment that starts
/// near the read's own diagonal needs an indel longer than the band's slack to leave
/// it, so one scoring below that may have been beaten off it.
const BAND_FLOOR: i32 = -(GAP_OPEN + GAP_EXTEND * (BAND_SLACK as i32 + 1));

/// Whether one of the read's exact k-mers (at its start, middle and end, unmasked) lies
/// on the haplotype off the band at its row: a copy of the read's bases elsewhere along
/// the haplotype, which an alignment reaches with no indel.
pub(crate) fn copy_off_band(read: &[u8], hap: &[u8], band: &Band) -> bool {
    let (n, m) = (read.len(), hap.len());
    if n < COPY_K || m < COPY_K {
        return false;
    }
    let mut offsets = vec![0, (n - COPY_K) / 2, n - COPY_K];
    offsets.dedup();
    offsets.into_iter().any(|o| {
        let kmer = &read[o..o + COPY_K];
        !kmer.contains(&WILD)
            && hap.windows(COPY_K).enumerate().any(|(q, w)| w == kmer && !band.holds(o + 1, q + 1))
    })
}

/// The banded programme, scores only.
pub(crate) fn banded(read: &[u8], hap: &[u8], band: &Band) -> Run {
    run::<true, false>(read, hap, band, &mut [], &[])
}

/// Whether a banded run stands by itself: its best path stays off the band's edge and
/// scores at least [`BAND_FLOOR`]. One that stands short of exact also needs no copy of
/// the read off the band ([`copy_off_band`]).
pub(crate) fn stands(r: &Run) -> bool {
    !r.clipped && r.score >= BAND_FLOOR
}

/// [`stands`] for a band `slack` columns each side of the read's diagonal (past the
/// length changes it carries).
pub(crate) fn stands_in(r: &Run, slack: i64) -> bool {
    !r.clipped && r.score >= -(GAP_OPEN + GAP_EXTEND * (slack as i32 + 1))
}

/// The full programme, scores only.
pub(crate) fn full_score(read: &[u8], hap: &[u8]) -> i32 {
    if read.is_empty() {
        return 0;
    }
    run::<false, false>(read, hap, &Band::full(read.len(), hap.len()), &mut [], &[]).score
}

/// The best score of `read` (masked bases [`WILD`]) consumed whole inside `hap`, ends
/// free: in `band` where it stands (and is exact, or no copy of the read lies off it),
/// else over every column. The rule `haplotypes::Scorer` applies with its caches.
#[cfg(test)]
pub(crate) fn score(read: &[u8], hap: &[u8], band: Option<&Band>) -> i32 {
    if read.is_empty() {
        return 0;
    }
    if let Some(b) = band {
        let r = banded(read, hap, b);
        if stands(&r) && (r.score == 0 || !copy_off_band(read, hap, b)) {
            return r.score;
        }
    }
    full_score(read, hap)
}

/// What an alignment does with one read base.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Placed {
    /// Aligned to the haplotype base at this column (of the haplotype scored).
    At(usize),
    /// Inserted between aligned bases.
    Inserted,
    /// Left over at an end of the alignment.
    LeftOver,
}

/// One alignment of a read inside a haplotype, per read base.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Alignment {
    pub(crate) score: i32,
    pub(crate) placed: Vec<Placed>,
    /// Haplotype bases skipped right after each read base (a deletion in the read).
    pub(crate) skipped_after: Vec<u32>,
    /// What each read base costs there: its substitution, its share of an insertion (the
    /// opening on the first), or its cost left over.
    pub(crate) costs: Vec<i32>,
}

impl Alignment {
    /// The first and last haplotype columns the read's aligned bases cover.
    pub(crate) fn covers(&self) -> Option<(usize, usize)> {
        let at = |p: &Placed| if let Placed::At(c) = *p { Some(c) } else { None };
        Some((self.placed.iter().find_map(at)?, self.placed.iter().rev().find_map(at)?))
    }
}

/// [`score`]'s alignment, traced back: the same band rule, so the same score.
pub(crate) fn align(read: &[u8], hap: &[u8], band: Option<&Band>) -> Alignment {
    let (n, m) = (read.len(), hap.len());
    let mut placed = vec![Placed::LeftOver; n];
    let mut skipped_after = vec![0u32; n];
    let costs_of = |placed: &[Placed]| -> Vec<i32> {
        (0..n)
            .map(|i| match placed[i] {
                Placed::At(c) => substitution_cost(read[i], hap[c]),
                Placed::Inserted if i > 0 && placed[i - 1] == Placed::Inserted => -GAP_EXTEND,
                Placed::Inserted => -(GAP_OPEN + GAP_EXTEND),
                Placed::LeftOver => overhang_cost(read[i]),
            })
            .collect()
    };
    if n == 0 {
        return Alignment { score: 0, placed, skipped_after, costs: Vec::new() };
    }
    // The band's values are kept (the edge and the traceback are read back from them); the
    // full programme's directions (a cell off the band reads as none, as it holds none).
    if let Some(b) = band {
        thread_local! {
            static VALUES: std::cell::RefCell<Values> = std::cell::RefCell::new(Values::default());
        }
        // reused across calls on this thread (taken out for the call)
        let mut vals = VALUES.with(|v| std::mem::take(&mut *v.borrow_mut()));
        let r = run_values(read, hap, b, &mut vals);
        let stands = r.score >= BAND_FLOOR && !touches_edge(read, hap, b, &vals, r.score) && (r.score == 0 || !copy_off_band(read, hap, b));
        let out = stands.then(|| traceback_values(read, hap, b, &vals, r));
        VALUES.with(|v| *v.borrow_mut() = vals);
        if let Some(a) = out {
            return a;
        }
    }
    let full = Band::full(n, m);
    let (starts, cells) = full.trace_layout();
    let mut trace = vec![0u8; cells];
    let r = run::<false, true>(read, hap, &full, &mut trace, &starts);
    let b = &full;
    let dir_at = |i: usize, j: usize| if b.holds(i, j) { trace[starts[i] + j - b.lo[i]] } else { 0 };
    let Some((mut i, mut j)) = r.end else {
        let costs = costs_of(&placed);
        return Alignment { score: r.score, placed, skipped_after, costs };
    };
    // States: 0 D (the read base i - 1 aligned), 1 A, 2 E (deletion), 3 F (insertion).
    let mut state = 0u8;
    while i > 0 {
        let dir = dir_at(i, j);
        match state {
            0 => {
                placed[i - 1] = Placed::At(j - 1);
                if dir & D_START != 0 {
                    break;
                }
                i -= 1;
                j -= 1;
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
                // Haplotype base j - 1 skipped, after read base i - 1.
                skipped_after[i - 1] += 1;
                state = if dir & E_EXT != 0 { 2 } else { 1 };
                j -= 1;
            }
            _ => {
                // Read base i - 1 inserted.
                placed[i - 1] = Placed::Inserted;
                state = if dir & F_EXT != 0 { 3 } else { 1 };
                i -= 1;
            }
        }
    }
    let costs = costs_of(&placed);
    Alignment { score: r.score, placed, skipped_after, costs }
}

/// [`align`]'s traceback over a band's [`Values`]: each cell's move read back as [`run`]'s
/// directions hold it (D extends an aligned path or starts there; A is D where they are
/// equal, else E, else F; E extends a deletion where that beats opening one, F an
/// insertion likewise).
fn traceback_values(read: &[u8], hap: &[u8], band: &Band, vals: &Values, r: Run) -> Alignment {
    let n = read.len();
    let mut placed = vec![Placed::LeftOver; n];
    let mut skipped_after = vec![0u32; n];
    if let Some((mut i, mut j)) = r.end {
        // States: 0 D (the read base i - 1 aligned), 1 A, 2 E (deletion), 3 F (insertion).
        let mut state = 0u8;
        while i > 0 {
            // a cell off the band holds no direction: D, from an aligned path
            let c = vals.at(band, i, j).map(|x| cell_sources(read, hap, band, vals, i, j, x));
            match state {
                0 => {
                    placed[i - 1] = Placed::At(j - 1);
                    if c.is_some_and(|c| !c.from_a) {
                        break;
                    }
                    i -= 1;
                    j -= 1;
                    state = 1;
                }
                1 => {
                    state = match c {
                        Some(c) if c.d != c.a && c.e == c.a => 2,
                        Some(c) if c.d != c.a => 3,
                        _ => 0,
                    }
                }
                2 => {
                    skipped_after[i - 1] += 1;
                    state = if c.is_some_and(|c| c.ee > c.eo) { 2 } else { 1 };
                    j -= 1;
                }
                _ => {
                    placed[i - 1] = Placed::Inserted;
                    state = if c.is_some_and(|c| c.fe > c.fo) { 3 } else { 1 };
                    i -= 1;
                }
            }
        }
    }
    let costs = (0..n)
        .map(|i| match placed[i] {
            Placed::At(c) => substitution_cost(read[i], hap[c]),
            Placed::Inserted if i > 0 && placed[i - 1] == Placed::Inserted => -GAP_EXTEND,
            Placed::Inserted => -(GAP_OPEN + GAP_EXTEND),
            Placed::LeftOver => overhang_cost(read[i]),
        })
        .collect();
    Alignment { score: r.score, placed, skipped_after, costs }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The full programme, unbanded.
    fn full(read: &[u8], hap: &[u8]) -> i32 {
        score(read, hap, None)
    }

    #[test]
    fn an_exact_substring_scores_zero_wherever_it_lies() {
        assert_eq!(full(b"GGATCC", b"AAAGGATCCTTT"), 0);
        assert_eq!(full(b"AAAG", b"AAAGGATCCTTT"), 0, "at the haplotype's start: its ends are free");
        assert_eq!(full(b"CCTTT", b"AAAGGATCCTTT"), 0, "at its end");
    }

    #[test]
    fn mismatches_and_gaps_are_affine() {
        assert_eq!(full(b"GGTTCC", b"AAAGGATCCTTT"), MISMATCH);
        // one inserted base in the read: 6 + 1
        assert_eq!(full(b"CTAGGATCCTGA", b"TTCTAGGTCCTGATT"), -7);
        // two deleted haplotype bases: 6 + 2, not twice 6 + 1
        assert_eq!(full(b"CTAGGTCCTGA", b"TTCTAGGATTCCTGATT"), -8);
        // the read is consumed whole: bases past the haplotype's end are left over
        assert_eq!(full(b"TTTGG", b"AAATTT"), 2 * MISMATCH);
    }

    #[test]
    fn a_masked_base_costs_nothing_against_any_base() {
        assert_eq!(full(b"GGNTCC", b"AAAGGATCCTTT"), 0);
        assert_eq!(full(b"GGNNCC", b"AAAGGATCCTTT"), 0);
    }

    /// Deterministic pseudo-random bases.
    fn bases(n: usize, state: &mut u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(*state >> 33) as usize % 4]
            })
            .collect()
    }

    fn roll(state: &mut u64, n: u64) -> u64 {
        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (*state >> 33) % n
    }

    #[test]
    fn the_banded_score_is_the_full_programmes_on_random_cases() {
        let mut s = 304u64;
        for case in 0..400 {
            let hap = bases(300, &mut s);
            let start = 80 + roll(&mut s, 60) as usize;
            let mut read = hap[start..start + 100].to_vec();
            // Substitutions, an indel of up to 30 bases, masked bases, a junk end.
            for _ in 0..roll(&mut s, 5) {
                let k = roll(&mut s, 100) as usize;
                read[k] = b"ACGT"[roll(&mut s, 4) as usize];
            }
            let at = 10 + roll(&mut s, 80) as usize;
            let k = roll(&mut s, 31) as usize;
            match roll(&mut s, 3) {
                0 => read.splice(at..at, bases(k, &mut s)).for_each(drop),
                1 => read.drain(at..(at + k).min(read.len())).for_each(drop),
                _ => {}
            }
            if roll(&mut s, 4) == 0 {
                let k = roll(&mut s, 25) as usize;
                let junk = bases(k, &mut s);
                read.splice(0..k.min(read.len()), junk).for_each(drop);
            }
            if roll(&mut s, 3) == 0 {
                let k = roll(&mut s, read.len() as u64) as usize;
                read[k] = WILD;
            }
            // The read's own diagonal, as its CIGAR would give it (no indel written):
            // sometimes off by a few columns, as where a haplotype's own length change
            // shifts it, which the band's width holds (as it holds the change).
            let shift = roll(&mut s, 9) as i64 - 4;
            let centers: Vec<i64> =
                (0..=read.len()).map(|i| start as i64 + i as i64 + shift).collect();
            let band = Band::around(&centers, BAND_SLACK + shift.abs(), hap.len());
            assert_eq!(
                score(&read, &hap, Some(&band)),
                full(&read, &hap),
                "case {case}: banded and full programmes disagree"
            );
            let a = align(&read, &hap, Some(&band));
            assert_eq!(a.score, full(&read, &hap), "case {case}: the traceback's score");
            assert_eq!(a.placed.len(), read.len());
        }
    }

    /// [`align`] as first written: the directions in an `(n + 1) * (m + 1)` matrix, every
    /// cell off the band reading as none.
    fn align_reference(read: &[u8], hap: &[u8], band: Option<&Band>) -> Alignment {
        let (n, m) = (read.len(), hap.len());
        let mut placed = vec![Placed::LeftOver; n];
        let mut skipped_after = vec![0u32; n];
        let costs_of = |placed: &[Placed]| -> Vec<i32> {
            (0..n)
                .map(|i| match placed[i] {
                    Placed::At(c) => substitution_cost(read[i], hap[c]),
                    Placed::Inserted if i > 0 && placed[i - 1] == Placed::Inserted => -GAP_EXTEND,
                    Placed::Inserted => -(GAP_OPEN + GAP_EXTEND),
                    Placed::LeftOver => overhang_cost(read[i]),
                })
                .collect()
        };
        if n == 0 {
            return Alignment { score: 0, placed, skipped_after, costs: Vec::new() };
        }
        let width = m + 1;
        // the matrix's own layout: row i's column j at i * width + j
        let matrix = |b: &Band| -> Vec<usize> { (0..=n).map(|i| i * width + b.lo[i]).collect() };
        let mut trace = vec![0u8; (n + 1) * (m + 1)];
        let mut r = None;
        if let Some(b) = band {
            let banded = run::<true, true>(read, hap, b, &mut trace, &matrix(b));
            if stands(&banded) && (banded.score == 0 || !copy_off_band(read, hap, b)) {
                r = Some(banded);
            }
        }
        let full = Band::full(n, m);
        let r = r.unwrap_or_else(|| run::<false, true>(read, hap, &full, &mut trace, &matrix(&full)));
        let Some((mut i, mut j)) = r.end else {
            let costs = costs_of(&placed);
            return Alignment { score: r.score, placed, skipped_after, costs };
        };
        let mut state = 0u8;
        while i > 0 {
            let dir = trace[i * width + j];
            match state {
                0 => {
                    placed[i - 1] = Placed::At(j - 1);
                    if dir & D_START != 0 {
                        break;
                    }
                    i -= 1;
                    j -= 1;
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
                    skipped_after[i - 1] += 1;
                    state = if dir & E_EXT != 0 { 2 } else { 1 };
                    j -= 1;
                }
                _ => {
                    placed[i - 1] = Placed::Inserted;
                    state = if dir & F_EXT != 0 { 3 } else { 1 };
                    i -= 1;
                }
            }
        }
        let costs = costs_of(&placed);
        Alignment { score: r.score, placed, skipped_after, costs }
    }

    #[test]
    fn the_traceback_kept_over_the_band_is_the_whole_matrixs() {
        let mut s = 6262u64;
        let mut clipped = 0;
        for case in 0..3000 {
            let m = 20 + roll(&mut s, 280) as usize;
            let unit = 1 + roll(&mut s, 4) as usize;
            let motif = bases(unit, &mut s);
            let hap: Vec<u8> = if roll(&mut s, 3) == 0 { (0..m).map(|i| motif[i % unit]).collect() } else { bases(m, &mut s) };
            let len = (1 + roll(&mut s, 120) as usize).min(m);
            let start = roll(&mut s, (m - len + 1) as u64) as usize;
            let mut read = hap[start..start + len].to_vec();
            for _ in 0..roll(&mut s, 6) {
                let k = roll(&mut s, read.len() as u64) as usize;
                read[k] = if roll(&mut s, 4) == 0 { WILD } else { b"ACGT"[roll(&mut s, 4) as usize] };
            }
            if read.len() > 4 {
                let at = roll(&mut s, read.len() as u64) as usize;
                let k = roll(&mut s, 12) as usize;
                match roll(&mut s, 3) {
                    0 => read.splice(at..at, bases(k, &mut s)).for_each(drop),
                    1 => read.drain(at..(at + k).min(read.len())).for_each(drop),
                    _ => {}
                }
            }
            if read.is_empty() {
                read.push(b'A');
            }
            let shift = roll(&mut s, 13) as i64 - 6;
            let centers: Vec<i64> = (0..=read.len()).map(|i| start as i64 + i as i64 + shift).collect();
            let band = Band::around(&centers, roll(&mut s, 20) as i64, hap.len());
            for b in [Some(&band), None] {
                assert_eq!(align(&read, &hap, b), align_reference(&read, &hap, b), "case {case} banded={}", b.is_some());
            }
            // the values' score, end and edge flag against the edge bits carried through
            // every cell
            let tracked = banded(&read, &hap, &band);
            let mut vals = Values::default();
            let r = run_values(&read, &hap, &band, &mut vals);
            assert_eq!((r.score, r.end), (tracked.score, tracked.end), "case {case}: score and end");
            assert_eq!(touches_edge(&read, &hap, &band, &vals, r.score), tracked.clipped, "case {case}: the edge");
            clipped += usize::from(tracked.clipped);
        }
        assert!(clipped > 300, "the cases reach the band's edge ({clipped})");
    }

    #[test]
    fn a_band_that_would_cut_the_optimum_off_falls_back_to_the_full_programme() {
        let mut s = 7u64;
        let hap = bases(400, &mut s);
        // The read carries a 40-base deletion its band (16 columns each side) cannot follow.
        let read: Vec<u8> = [&hap[100..150], &hap[190..240]].concat();
        let centers: Vec<i64> = (0..=read.len()).map(|i| 100 + i as i64).collect();
        let band = Band::around(&centers, BAND_SLACK, hap.len());
        assert_eq!(score(&read, &hap, Some(&band)), -(GAP_OPEN + 40 * GAP_EXTEND));
    }

    #[test]
    fn a_copy_of_the_read_off_its_band_is_found() {
        // A near-copy of [100, 200) at [230, 330), one base apart; the read is the second
        // copy exactly while its own alignment puts it at the first.
        let mut s = 99u64;
        let mut hap = bases(400, &mut s);
        let seg: Vec<u8> = hap[100..200].to_vec();
        hap.splice(230..330, seg);
        hap[150] = if hap[280] == b'A' { b'C' } else { b'A' };
        let read = hap[230..330].to_vec();
        let centers: Vec<i64> = (0..=read.len()).map(|i| 100 + i as i64).collect();
        let band = Band::around(&centers, BAND_SLACK, hap.len());
        assert_eq!(score(&read, &hap, Some(&band)), 0);
    }

    #[test]
    fn read_bases_past_the_haplotypes_end_cost_a_mismatch_each() {
        // An overhang is no cheaper than mismatching the bases (one base past the end
        // costs a mismatch, not a gap opening), masked ones aside.
        assert_eq!(full(b"TTTG", b"AAATTT"), MISMATCH);
        assert_eq!(full(b"TTTGG", b"AAATTT"), 2 * MISMATCH);
        assert_eq!(full(b"NNNNNNAAAT", b"AAATTT"), 0);
        assert_eq!(full(b"AAATNNNN", b"AAATTT"), 0);
    }

    #[test]
    fn the_traceback_gives_each_read_base_its_haplotype_base() {
        // A 2-base deletion in the read after its 7th base.
        let a = align(b"CTAGGTCCTGA", b"TTCTAGGATTCCTGATT", None);
        assert_eq!(a.score, -8);
        assert_eq!(a.skipped_after.iter().sum::<u32>(), 2);
        assert!(a.placed.iter().all(|p| matches!(p, Placed::At(_))));
        // An inserted base.
        let a = align(b"CTAGGATCCTGA", b"TTCTAGGTCCTGATT", None);
        assert_eq!(a.score, -7);
        assert_eq!(a.placed.iter().filter(|p| **p == Placed::Inserted).count(), 1);
    }
}
