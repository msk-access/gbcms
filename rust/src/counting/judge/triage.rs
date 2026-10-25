//! Coloured k-mer triage in front of the graph scorer (`GBCMS_JUDGE_TRIAGE`): every k-mer
//! along the paths of each state's graph over the reads' reach, coloured by the states
//! whose graph holds it. A read whose confident k-mers all lie in exactly one state's
//! graph, with every k-mer over X's site known, is called for that state without the
//! programme.
//!
//! **k:** per site, the smallest odd k from [`K_MIN`] whose k-mers are each unique in the
//! reference over the reads' reach (so a k-mer names one place there, and a repeat no k
//! spans leaves the site to the scorer), and shorter than the reads (a read holds several).
//!
//! What triage needs of a read, so its call is the scorer's: the read reaches past X's
//! site on both sides and reads every base within k of it confidently (a masked base
//! there leaves the call to the margin); every k-mer holding a base of the site is known
//! (an unknown one is an allele no state holds: the fit check's case); at most
//! [`MAX_UNKNOWN`] k-mers elsewhere are unknown (a read that fits no state well is the
//! scorer's to call depth only); and the colours of its known k-mers meet in one state.

use std::collections::HashMap;

use super::decide::Call;
use super::graph::{Graph, GraphRead};
use super::haplotypes::{Class, ReadView, Site};
use super::score::WILD;

/// The shortest k tried, and the longest (two bits a base in 128).
const K_MIN: usize = 15;
const K_MAX: usize = 63;
/// The k-mers ending at one graph base followed at most (a dense cluster's paths
/// multiply), and k-mers kept per site: past either, the site is left to the scorer.
const MAX_PER_BASE: usize = 256;
const MAX_KMERS: usize = 1 << 20;

/// A read's unknown k-mers away from X's site at most, in multiples of k (two isolated
/// errors' worth). Measurement only, TEMPORARY: `GBCMS_JUDGE_TRIAGE_UNKNOWN` sets it (0: a
/// read triage decides has no unknown k-mer anywhere).
const MAX_UNKNOWN: usize = 2;

/// Measurement only, TEMPORARY: `GBCMS_JUDGE_TRIAGE_EXACT=on` confirms a triage call only
/// where the whole read is an exact path of the named state's graph (a masked base matching
/// any), so the call is the scorer's: that state scores 0, every other lacks a k-mer over
/// the site (one confident edit or more), and neither the fit check nor the fits-nothing
/// rule applies to a read scoring 0.
fn exact_wanted() -> bool {
    std::env::var("GBCMS_JUDGE_TRIAGE_EXACT").is_ok_and(|v| v == "on")
}

/// Whether some path of `g` spells `read` exactly, a masked base on either side matching
/// any: Shift-And over the DAG, one bit per read base.
fn spelled_exactly(g: &Graph, read: &[u8]) -> bool {
    let n = read.len();
    if n == 0 {
        return false;
    }
    let words = n.div_ceil(64);
    let mut masks = vec![vec![0u64; words]; 4];
    for (j, &b) in read.iter().enumerate() {
        for (c, m) in masks.iter_mut().enumerate() {
            if b == WILD || code(b) == Some(c as u128) {
                m[j / 64] |= 1 << (j % 64);
            }
        }
    }
    let all = vec![u64::MAX; words];
    let top = (n - 1) / 64;
    let top_bit = 1u64 << ((n - 1) % 64);
    let mut state: Vec<u64> = vec![0; g.len() * words];
    for v in 0..g.len() {
        let mut from = vec![0u64; words];
        for &u in g.preds(v) {
            let u = u as usize;
            for w in 0..words {
                from[w] |= state[u * words + w];
            }
        }
        let mask = match code(g.bases[v]) {
            Some(c) => &masks[c as usize],
            None => &all,
        };
        let mut carry = 1u64;
        for w in 0..words {
            let shifted = (from[w] << 1) | carry;
            carry = from[w] >> 63;
            state[v * words + w] = shifted & mask[w];
        }
        if state[v * words + top] & top_bit != 0 {
            return true;
        }
    }
    false
}

fn max_unknown() -> usize {
    std::env::var("GBCMS_JUDGE_TRIAGE_UNKNOWN").ok().and_then(|v| v.parse().ok()).unwrap_or(MAX_UNKNOWN)
}

/// A site's coloured k-mer index.
#[derive(Debug, Clone)]
pub(crate) struct Index {
    k: usize,
    colours: HashMap<u128, u64>,
    classes: Vec<Class>,
    /// Unknown k-mers away from the site a decided read may hold.
    max_unknown: usize,
    /// Whether a call is confirmed by the read spelled exactly along the state's graph.
    exact: bool,
}

/// Triage's call for a read, with the lowest quality among the read's bases over X's site
/// and the unknown k-mers it holds away from the site.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Triaged {
    pub(crate) call: Call,
    pub(crate) qual: u8,
    pub(crate) unknown: u16,
}

#[inline]
fn code(b: u8) -> Option<u128> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

impl Index {
    /// The index of `site`'s graphs over the reach of reads of up to `l` bases around X;
    /// None where no k fits or the paths are too many.
    pub(crate) fn build(site: &Site, l: i64) -> Option<Index> {
        let graphs = &site.graphs;
        if graphs.is_empty() || graphs.len() > 64 {
            return None;
        }
        let (xs, xe) = site.x_span();
        let (_, reference) = site.window();
        // Every base the programme may align a read over X's site to: the read reaches a
        // read's length past X, its programme a read's length past the read, each widened by
        // the length changes of the rows there (a path through a deletion reads the
        // reference that much further on).
        let mut shifts = 0;
        loop {
            let more = graphs.iter().map(|g| g.graph.shifts_near(xs, xe, 2 * l + 2 * shifts)).max().unwrap_or(0);
            if more == shifts {
                break;
            }
            shifts = more;
        }
        let reach = ((xs - 2 * l - shifts).max(0), (xe + 2 * l + shifts).min(reference.len() as i64));
        let k = choose_k(&reference[reach.0 as usize..reach.1 as usize], l)?;
        let mut colours: HashMap<u128, u64> = HashMap::new();
        for (s, st) in graphs.iter().enumerate() {
            for c in kmers_of(&st.graph, k, reach)? {
                *colours.entry(c).or_insert(0) |= 1 << s;
            }
            if colours.len() > MAX_KMERS {
                return None;
            }
        }
        Some(Index { k, colours, classes: graphs.iter().map(|g| g.class).collect(), max_unknown: max_unknown() * k, exact: exact_wanted() })
    }

    /// Triage's call for `read` (its graph view `gr`), or None for the scorer to call.
    pub(crate) fn decide(&self, site: &Site, read: &ReadView, gr: &GraphRead, min_baseq: u8) -> Option<Triaged> {
        let k = self.k;
        let n = read.seq.len();
        if n < k + 2 {
            return None;
        }
        let (ss, se) = site.site_window();
        let rp = &gr.rp;
        // the read's bases over X's site, and bases past it each side
        let first = rp.partition_point(|&p| p < ss);
        let after = rp.partition_point(|&p| p < se);
        if first == 0 || after >= n || first >= after {
            return None;
        }
        let (t0, t1) = (first, after - 1);
        if read.seq[t0.saturating_sub(k - 1)..(t1 + k).min(n)].contains(&WILD) {
            return None;
        }
        let mask: u128 = if 2 * k == 128 { u128::MAX } else { (1u128 << (2 * k)) - 1 };
        let (mut inter, mut unknown) = (u64::MAX, 0usize);
        let (mut c, mut run) = (0u128, 0usize);
        for (t, &b) in read.seq.iter().enumerate() {
            match code(b) {
                Some(x) => {
                    c = ((c << 2) | x) & mask;
                    run += 1;
                }
                None => run = 0,
            }
            if run < k {
                continue;
            }
            let s = t + 1 - k; // the k-mer's first base
            let over_site = s <= t1 && t >= t0;
            match self.colours.get(&c) {
                None if over_site => return None,
                None => unknown += 1,
                Some(&col) => inter &= col,
            }
        }
        if unknown > self.max_unknown {
            return None;
        }
        let live = if self.classes.len() == 64 { u64::MAX } else { (1u64 << self.classes.len()) - 1 };
        let inter = inter & live;
        if inter.count_ones() != 1 {
            return None;
        }
        let state = inter.trailing_zeros() as usize;
        if self.exact && !spelled_exactly(&site.graphs[state].graph, &read.seq) {
            return None;
        }
        let call = match self.classes[state] {
            Class::Alt => Call::Alt,
            Class::Ref => Call::Ref,
            Class::Other => Call::Other,
        };
        let floor = min_baseq.max(1);
        let qual = if call == Call::Other { 0 } else { read.quals[t0..=t1].iter().copied().min().unwrap_or(floor).max(floor) };
        Some(Triaged { call, qual, unknown: unknown.min(u16::MAX as usize) as u16 })
    }
}

/// The smallest odd k from [`K_MIN`], shorter than the reads, whose k-mers are each unique
/// in `region` (k-mers holding a base other than A, C, G or T aside).
fn choose_k(region: &[u8], l: i64) -> Option<usize> {
    let top = K_MAX.min(usize::try_from(l).ok()?.saturating_sub(1));
    (K_MIN..=top).step_by(2).find(|&k| {
        let mask: u128 = if 2 * k == 128 { u128::MAX } else { (1u128 << (2 * k)) - 1 };
        let mut seen = std::collections::HashSet::new();
        let (mut c, mut run) = (0u128, 0usize);
        for &b in region {
            match code(b) {
                Some(x) => {
                    c = ((c << 2) | x) & mask;
                    run += 1;
                }
                None => run = 0,
            }
            if run >= k && !seen.insert(c) {
                return false;
            }
        }
        true
    })
}

/// Every k-mer along a path of `g` over the reference span `reach` (window coordinates);
/// a wildcard base stands for each base. None past the bounds.
fn kmers_of(g: &Graph, k: usize, reach: (i64, i64)) -> Option<Vec<u128>> {
    let mask: u128 = if 2 * k == 128 { u128::MAX } else { (1u128 << (2 * k)) - 1 };
    let n = g.len();
    let inside = |v: usize| reach.0 <= g.coord[v] && g.coord[v] < reach.1;
    // per base, the distinct (length up to k, code) suffixes of the paths ending there
    let mut ends: Vec<Vec<(u8, u128)>> = vec![Vec::new(); n];
    let mut out: Vec<u128> = Vec::new();
    for v in 0..n {
        if !inside(v) {
            continue;
        }
        let b = g.bases[v];
        let bases: Vec<u128> = match code(b) {
            Some(x) => vec![x],
            None if b == WILD => vec![0, 1, 2, 3],
            None => Vec::new(),
        };
        let preds: Vec<usize> = g.preds(v).iter().map(|&u| u as usize).filter(|&u| inside(u)).collect();
        let mut cur: Vec<(u8, u128)> = Vec::new();
        for &x in &bases {
            if preds.is_empty() {
                cur.push((1, x));
            }
            for &u in &preds {
                for &(len, c) in &ends[u] {
                    let nl = (len as usize + 1).min(k) as u8;
                    cur.push((nl, ((c << 2) | x) & mask));
                }
            }
        }
        cur.sort_unstable();
        cur.dedup();
        if cur.len() > MAX_PER_BASE {
            return None;
        }
        out.extend(cur.iter().filter(|e| e.0 as usize == k).map(|e| e.1));
        if out.len() > MAX_KMERS {
            return None;
        }
        ends[v] = cur;
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::super::decide::judge;
    use super::super::haplotypes::tests::random_seq;
    use super::super::haplotypes::{Candidates, Row, ScorerKind};
    use super::*;
    use rust_htslib::bam::record::Cigar;

    fn view(h: &[u8], s: usize, len: usize, quals: &[u8]) -> ReadView {
        ReadView::new(&h[s..s + len], quals, s as i64, &[Cigar::Match(len as u32)], (0, len), 20).unwrap()
    }

    #[test]
    fn every_read_triage_decides_gets_the_scorers_call() {
        // An SNV X with a listed SNV beside it; reads of REF, of X, of the neighbour, of a
        // third base at X, some with an error or a masked base: triage decides the clean
        // ones, always as the scorer does, and leaves a read masked at X to it.
        let c = random_seq(800, 31);
        let other = |b: u8, k: usize| *b"ACGT".iter().filter(|&&x| x != b).nth(k).unwrap();
        let x = Row::new(400, &c[400..401], &[other(c[400], 0)]);
        let near = Row::new(430, &c[430..431], &[other(c[430], 1)]);
        let third = Row::new(400, &c[400..401], &[other(c[400], 2)]);
        let apply = |rows: &[&Row]| {
            let mut h = c.clone();
            for r in rows {
                h[r.pos as usize] = r.alt_al[0];
            }
            h
        };
        let haps = [apply(&[]), apply(&[&x]), apply(&[&near]), apply(&[&x, &near]), apply(&[&third])];
        let mut reads = Vec::new();
        for (k, h) in haps.iter().enumerate() {
            for s in [300usize, 320, 340, 360, 380] {
                let mut q = vec![30u8; 120];
                if s == 360 {
                    q[400 - s] = 2; // X's base masked
                }
                let mut h = h.clone();
                if s == 340 && k % 2 == 0 {
                    h[s + 5] = other(h[s + 5], 0); // an error far from X
                }
                reads.push(view(&h, s, 120, &q));
            }
        }
        let cand = Candidates { listed: vec![near], ..Candidates::default() };
        let site = Site::build_with(ScorerKind::Graph, x, &cand, (0, c.clone()), &reads).unwrap();
        let ix = Index::build(&site, 120).expect("a unique k fits a random window");
        let (mut decided, mut masked_left) = (0, true);
        for (i, r) in reads.iter().enumerate() {
            let gr = site.graph_read(r);
            if let Some(t) = ix.decide(&site, r, &gr, 20) {
                decided += 1;
                assert_eq!(t.call, judge(&site, r, 20, false).call, "read {i}");
                masked_left &= i % 5 != 3;
            }
        }
        assert!(masked_left, "a read masked at X is the scorer's");
        assert!(decided >= 15, "{decided} of {} decided", reads.len());
    }

    #[test]
    fn a_read_is_spelled_exactly_only_along_one_path() {
        use super::super::graph::GraphBuilder;
        // AAAC[G|T]TTTA[C|-]GG: reads along a path are spelled; one mixing both alleles of
        // a bubble, or with a mismatch, is not; a masked read base matches any.
        let mut b = GraphBuilder::default();
        let l = b.segment(None, (0, 4), b"AAAC".to_vec());
        let g = b.segment(None, (4, 5), b"G".to_vec());
        let t = b.segment(Some(1), (4, 5), b"T".to_vec());
        let m = b.segment(None, (5, 9), b"TTTA".to_vec());
        let c = b.segment(None, (9, 10), b"C".to_vec());
        let r = b.segment(None, (10, 12), b"GG".to_vec());
        for s in [g, t] {
            b.edge(l, s);
            b.edge(s, m);
        }
        b.edge(m, c);
        b.edge(c, r);
        b.edge(m, r);
        let gr = b.finish();
        for (read, want) in [
            (&b"AACGTTTACGG"[..], true),
            (b"AACTTTTAGG", true),
            (b"ACGTTTAGG", true),
            (b"AACNTTTACG", true),
            (b"AACGTTTAAGG", false),
            (b"AACGTTTACCG", false),
        ] {
            assert_eq!(spelled_exactly(&gr, read), want, "{}", String::from_utf8_lossy(read));
        }
    }

    #[test]
    fn no_k_fits_a_repeat_longer_than_every_k() {
        let mut c = random_seq(800, 33);
        let unit = b"CAG";
        for i in 0..40 {
            c[350 + 3 * i..353 + 3 * i].copy_from_slice(unit);
        }
        assert!(choose_k(&c[200..600], 150).is_none());
        assert!(choose_k(&random_seq(400, 35), 150).is_some_and(|k| (K_MIN..=K_MAX).contains(&k)));
    }
}
