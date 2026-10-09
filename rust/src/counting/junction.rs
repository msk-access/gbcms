//! The junction judge for pure deletions of [`LONG_DELETION`] bases or more (RJ-24).
//!
//! A read is REF or ALT for such a deletion by its own bases at the POS breakpoint,
//! however the aligner wrote it (a D op, a soft clip, a split read):
//!
//! - **One window per allele, of equal length, read from the left flank.** [`K`]
//!   flank bases ending at the anchor, the microhomology (the bases the deletion can
//!   slide right through, which both alleles share), then [`K`] bases of the
//!   allele's own continuation: the deleted bases for REF, the right flank for ALT.
//!   A REF and an ALT molecule hold their windows from the same read starts, so
//!   neither allele is favoured by where reads start.
//! - **Nothing past the window is read**, for either allele: checking REF through
//!   every deleted base it covers would expose it to errors the ALT reading never
//!   meets (REJ-20261009-003).
//!
//! - **Placement.** A read is read from the query position of its aligned first
//!   window base, so both windows start where the read's flank does. A read whose
//!   soft clip holds the flank, or whose own deletion covers it, is placed from the
//!   other side at each allele's length (RJ-13); a placement counts where its own
//!   allele fits, else the closer one does. Hard-clipped bases are absent: such a
//!   record holds no window on its own.
//! - **Masked bases.** Bases below `--min-baseq` (and N) fit anything, the one
//!   quality rule every backend shares; so do the window offsets the neighbouring
//!   change guard masks. An ALT call also needs the read's bases to favour ALT over
//!   REF by one minimum-quality base's evidence (RJ-16).
//! - **Outcomes.** A read holding the windows is REF or ALT when it fits exactly one;
//!   fitting neither, it carries another allele (partial evidence when closer to
//!   ALT). A read that does not hold them, or whose windows cross a splice, is depth
//!   only. A decided read whose windows read soft-clipped bases inside its fragment
//!   counts toward depth in DNA (clip admission). A read that fits REF but whose own
//!   alignment deletes or inserts [`MIN_OWN_INDEL`] bases or more inside the deleted
//!   span carries another allele (C26): partial evidence.

use log::trace;
use rust_htslib::bam::record::{Cigar, Record};

use crate::shared::bam_utils::fragment_query_span;
use crate::types::Variant;

use super::carrier::{
    aligned_query_range, one_base_evidence, readable_range, score, Guard, ReadRules, Reading, GUARD_MIN_READS,
    GUARD_MIN_SHARE, GUARD_REF_SHARE, WILD,
};
use super::observed::{canonical, ObservedAllele};
use super::utils::{find_read_pos, median_qual, ClassifyPhase, ClassifyResult};
use super::variant_checks::deleted_anchor_query;
use super::window::{allele_kind, AlleleKind};

/// Deleted bases from which a pure deletion is judged at its junction.
pub(crate) const LONG_DELETION: usize = 50;
/// Clean bases each window holds before the junction (the flank) and past the
/// microhomology (the allele's own continuation).
const K: usize = 4;
/// A read's own indel at least this long inside the deleted span makes it another
/// allele (C26); shorter ones in a deleted repeat are stutter or alignment noise, as
/// the CIGAR checks treat 1-4 bp deletions in their scan window.
const MIN_OWN_INDEL: u32 = 5;

/// A long deletion's junction windows on the prepared reference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Junction {
    /// Genome position of the first deleted base at the leftmost placement.
    lo: i64,
    /// Deleted bases.
    len: usize,
    /// Bases the deletion slides right through (the microhomology).
    mh: usize,
    /// REF's window: K flank bases, the microhomology, K deleted bases.
    pub(crate) ref_w: Vec<u8>,
    /// ALT's window: the same flank and microhomology, K right-flank bases.
    pub(crate) alt_w: Vec<u8>,
}

impl Junction {
    /// Genome position of both windows' first base (the flank's first).
    pub(crate) fn first(&self) -> i64 {
        self.lo - K as i64
    }

    /// Genome position of an allele's window base at `offset` (0-based).
    pub(crate) fn genome(&self, alt: bool, offset: usize) -> i64 {
        if alt && offset >= K {
            self.lo + self.len as i64 + (offset - K) as i64
        } else {
            self.first() + offset as i64
        }
    }

    /// Genome position of an allele's window's last base.
    pub(crate) fn last(&self, alt: bool) -> i64 {
        self.genome(alt, self.ref_w.len() - 1)
    }

    /// The bases where the alleles differ over every placement, genomic [lo, hi):
    /// the deleted bases and the microhomology.
    fn span(&self) -> (i64, i64) {
        (self.lo, self.lo + (self.len + self.mh) as i64)
    }
}

/// Whether the junction judge reads `v`: a pure deletion of [`LONG_DELETION`] bases
/// or more, prepared against a reference (prep's `event_ref`).
pub(crate) fn judges(v: &Variant) -> bool {
    v.event_ref.is_some()
        && allele_kind(&v.ref_allele, &v.alt_allele) == Some(AlleleKind::Deletion)
        && v.ref_allele.len() - v.alt_allele.len() >= LONG_DELETION
}

/// The deletion's junction windows on its prepared reference. None when it is not
/// judged here, or the reference does not hold the windows.
pub(crate) fn junction(v: &Variant) -> Option<Junction> {
    if !judges(v) {
        return None;
    }
    let (start, seq) = v.event_ref.as_ref()?;
    let reference: Vec<u8> = seq.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let at = |g: i64| usize::try_from(g - start).ok().filter(|&i| i < reference.len());
    let len = v.ref_allele.len() - v.alt_allele.len();
    // The deleted bases follow the anchor; slide them to their leftmost placement,
    // then count how far they slide right (prep left-aligns, so the left slide is
    // a guard, not the rule).
    let mut lo = at(v.pos + v.alt_allele.len() as i64)?;
    while lo > 0 && reference.get(lo - 1) == reference.get(lo - 1 + len) {
        lo -= 1;
    }
    let mut mh = 0;
    while lo + mh + len < reference.len() && reference[lo + mh] == reference[lo + mh + len] {
        mh += 1;
    }
    let reach = mh + K;
    if lo < K || lo + len + reach > reference.len() {
        return None;
    }
    let flank = &reference[lo - K..lo];
    let ref_w = [flank, &reference[lo..lo + reach]].concat();
    let alt_w = [flank, &reference[lo + len..lo + len + reach]].concat();
    Some(Junction { lo: start + lo as i64, len, mh, ref_w, alt_w })
}

/// A read's call at a long deletion (RJ-24), or None when the junction judge does
/// not judge the variant (the caller then uses the CIGAR checks).
pub(crate) fn classify(
    record: &Record,
    v: &Variant,
    quals: &[u8],
    min_baseq: u8,
    rules: &ReadRules,
) -> Option<ClassifyResult> {
    let j = junction(v)?;
    let seq = record.seq().as_bytes();
    if seq.is_empty() || quals.len() < seq.len() {
        // A record stored without its bases (SEQ '*') shows nothing here.
        return Some(ClassifyResult::no_coverage(ClassifyPhase::MaskedCompare));
    }
    let Some(readable) = readable_range(record, rules) else {
        return Some(ClassifyResult::no_coverage(ClassifyPhase::MaskedCompare));
    };
    // Offsets the guard masks (a change recurring beside the event, a separate
    // event) fit any base, in each allele's window at its own genome position.
    let masked = rules.guard.and_then(|g| g.get(v)).map_or(&[][..], |g| g.masked.as_slice());
    let window = |alt: bool, w: &[u8]| -> Vec<u8> {
        w.iter()
            .enumerate()
            .map(|(o, &b)| if masked.contains(&j.genome(alt, o)) { WILD } else { b })
            .collect()
    };
    let (ref_w, alt_w) = (window(false, &j.ref_w), window(true, &j.alt_w));
    let n = ref_w.len();
    let read_from = |s: i64| -> Option<(Reading, Reading)> {
        let s = usize::try_from(s).ok().filter(|&s| s >= readable.0 && s + n <= readable.1)?;
        let (bases, qs, span) = (&seq[s..s + n], &quals[s..s + n], (s, s + n));
        Some((score(bases, qs, min_baseq, &ref_w, span), score(bases, qs, min_baseq, &alt_w, span)))
    };
    let reading = match find_read_pos(record, j.first()) {
        Some(q) => read_from(q as i64),
        None => {
            // The flank in the read's clip, or under its own deletion: placed from
            // the other side as each allele would sit there.
            let (first, _) = aligned_query_range(record)?;
            let deleted = deleted_anchor_query(record, j.first()).is_some();
            let placed: Vec<((Reading, Reading), bool)> = [false, true]
                .into_iter()
                .filter_map(|alt| {
                    let s = find_read_pos(record, j.last(alt))? as i64 - (n as i64 - 1);
                    (s < first as i64 || deleted).then(|| read_from(s)).flatten().map(|p| (p, alt))
                })
                .collect();
            placed
                .iter()
                .find(|((r, a), alt)| if *alt { a.mismatches == 0 } else { r.mismatches == 0 })
                .or_else(|| placed.iter().min_by_key(|((r, a), _)| r.mismatches.min(a.mismatches)))
                .map(|(p, _)| p.clone())
        }
    };
    let Some((r, a)) = reading else {
        // The read does not hold the windows: depth only, no allele, no mFSD class.
        let mut out = ClassifyResult::neither(ClassifyPhase::MaskedCompare);
        out.uninformative = true;
        return Some(out);
    };
    let (lo, hi) = r.span;
    if crosses_splice(record, lo, hi) {
        // Spliced inside the windows: the read skips bases the windows compare.
        let mut out = ClassifyResult::neither(ClassifyPhase::CigarRecon);
        out.uninformative = true;
        return Some(out);
    }
    let qual = median_qual(&quals[lo..hi], min_baseq);
    // The read's evidence for ALT over REF (log10), every base weighed by its quality.
    let evidence = r.mismatch_weight - a.mismatch_weight;
    let (fits_ref, fits_alt) = (r.mismatches == 0, a.mismatches == 0);
    let mut out = if fits_alt && !fits_ref && evidence < one_base_evidence(min_baseq) {
        ClassifyResult::neither(ClassifyPhase::MaskedCompare)
    } else if fits_alt && !fits_ref {
        ClassifyResult::is_alt(qual, ClassifyPhase::MaskedCompare)
    } else if fits_ref && !fits_alt {
        ClassifyResult::is_ref(qual, ClassifyPhase::MaskedCompare)
    } else if a.mismatches < r.mismatches {
        // Another allele, closer to ALT: partial evidence.
        ClassifyResult::neither_with_nearby(qual, ClassifyPhase::MaskedCompare)
    } else {
        ClassifyResult::neither(ClassifyPhase::MaskedCompare)
    };
    if out.is_ref && own_indel_inside(record, j.span()) {
        // The read's own alignment deletes or inserts bases inside the deleted span:
        // another allele (C26), whatever its window shows.
        out = ClassifyResult::neither_with_nearby(qual, ClassifyPhase::MaskedCompare);
    }
    out.has_n_base = r.had_n || a.had_n;
    // A decided read whose windows read soft-clipped bases counts toward depth (in
    // DNA) when every base read lies inside its fragment: past it a clip is adapter.
    out.clip_admissible = (out.is_ref || out.is_alt)
        && aligned_query_range(record).is_some_and(|(first, after)| lo < first || hi > after)
        && fragment_query_span(record).is_some_and(|(f0, f1)| lo >= f0 && hi <= f1);
    out.junction_judged = true;
    Some(out)
}

/// The neighbouring-change guard for a long deletion (the carrier rule's, RJ-23,
/// over the junction windows): from the reads `counted` admits, once per fragment,
/// each read holding the windows from its aligned first window base is taken as the
/// allele whose continuation it follows (its [`K`] bases past the microhomology fit
/// that allele to within one base and the other by two or more), and its confident
/// changes are tallied per window offset. A change on at least [`GUARD_MIN_READS`]
/// fragments and [`GUARD_MIN_SHARE`] of that allele's reads, with neither
/// neighbouring offset recurring too (a run of changes is an indel at the junction:
/// another allele, not masked), is masked at its offset in both windows. On the ALT
/// reads alone (under [`GUARD_REF_SHARE`] of REF's) it also makes the larger allele
/// the diagnostic names. None when nothing recurs, or the variant is not judged
/// here.
pub(crate) fn guard_for(
    read_cache: &[Record],
    v: &Variant,
    counted: &dyn Fn(&Record) -> bool,
    min_baseq: u8,
) -> Option<Guard> {
    let j = junction(v)?;
    let (start, seq) = v.event_ref.as_ref()?;
    let reference: Vec<u8> = seq.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let n = j.ref_w.len();
    let cont = K + j.mh;
    let windows = [&j.ref_w, &j.alt_w];
    // Per allele (REF 0, ALT 1) and offset: fragments reading it, and changing it;
    // the bases ALT's changes show (A C G T).
    let (mut cover, mut change) = ([vec![0u32; n], vec![0u32; n]], [vec![0u32; n], vec![0u32; n]]);
    let mut shown = vec![[0u32; 4]; n];
    let mut fragments: std::collections::HashSet<&[u8]> = std::collections::HashSet::new();
    let mut scanned = 0u32;
    for record in read_cache {
        if !counted(record) {
            continue;
        }
        let bases = record.seq().as_bytes();
        let quals = record.qual();
        let Some(q) = find_read_pos(record, j.first()).filter(|&q| q + n <= bases.len() && quals.len() >= q + n)
        else {
            continue;
        };
        let read: Vec<u8> = bases[q..q + n].iter().map(u8::to_ascii_uppercase).collect();
        let confident = |o: usize| read[o] != WILD && quals[q + o] >= min_baseq;
        let off = |w: &[u8], o: usize| confident(o) && read[o] != w[o];
        let misses = |w: &[u8]| (cont..n).filter(|&o| off(w, o)).count();
        let allele = match (misses(&j.ref_w), misses(&j.alt_w)) {
            (0..=1, 2..) => 0,
            (2.., 0..=1) => 1,
            _ => continue,
        };
        if !fragments.insert(record.qname()) {
            continue;
        }
        scanned += 1;
        for o in 0..n {
            cover[allele][o] += 1;
            if off(windows[allele], o) {
                change[allele][o] += 1;
                if let (1, Some(k)) = (allele, b"ACGT".iter().position(|&b| b == read[o])) {
                    shown[o][k] += 1;
                }
            }
        }
    }
    let share = |a: usize, o: usize| if cover[a][o] == 0 { 0.0 } else { change[a][o] as f64 / cover[a][o] as f64 };
    let recurs = |a: usize, o: usize| change[a][o] >= GUARD_MIN_READS && share(a, o) >= GUARD_MIN_SHARE;
    let neighbour_recurs = |a: usize, o: usize| (o > 0 && recurs(a, o - 1)) || (o + 1 < n && recurs(a, o + 1));
    let isolated = |a: usize, o: usize| !neighbour_recurs(a, o);
    let (mut masked, mut alt_only) = (Vec::new(), Vec::new());
    for o in 0..n {
        let hit = [0, 1].map(|a| recurs(a, o) && isolated(a, o));
        if hit[0] || hit[1] {
            for p in [j.genome(false, o), j.genome(true, o)] {
                if !masked.contains(&p) {
                    masked.push(p);
                }
            }
        }
        if hit[1] && share(0, o) < GUARD_REF_SHARE {
            alt_only.push(o);
        }
    }
    if masked.is_empty() {
        return None;
    }
    // The larger allele: the ALT haplotype from the windows' first base through its
    // farthest change, each changed base as most ALT reads show it, against the
    // genome over the same span (the deleted bases included).
    let larger = alt_only.iter().max().and_then(|&far| {
        let alt: Vec<u8> = (0..=far)
            .map(|o| match alt_only.contains(&o) {
                true => b"ACGT"[(0..4).max_by_key(|&k| shown[o][k]).unwrap_or(0)],
                false => j.alt_w[o],
            })
            .collect();
        let (g0, g1) = (j.first(), j.genome(true, far));
        let at = |g: i64| usize::try_from(g - start).ok().filter(|&i| i < reference.len());
        let ref_part = reference.get(at(g0)?..=at(g1)?)?.to_vec();
        let base = |g: i64| at(g).map(|i| reference[i]);
        let carriers = alt_only.iter().map(|&o| change[1][o]).min().unwrap_or(0);
        let given = alt_only.iter().map(|&o| cover[1][o] - change[1][o]).max().unwrap_or(0);
        let text = |x: Vec<u8>| String::from_utf8_lossy(&x).into_owned();
        canonical(g0, ref_part, alt, &base).map(|(pos, r, a)| ObservedAllele {
            pos,
            ref_allele: text(r),
            alt_allele: text(a),
            carriers,
            given_carriers: given,
        })
    });
    trace!(
        "junction guard at {}:{} {}>{}: {} window offsets masked, {} on the ALT reads alone ({} fragments read)",
        v.chrom, v.pos + 1, v.ref_allele, v.alt_allele, masked.len() / 2, alt_only.len(), scanned,
    );
    Some(Guard::new(v, masked, larger, scanned))
}

/// Whether the read's own alignment carries a deletion or an insertion of
/// [`MIN_OWN_INDEL`] bases or more starting inside the genomic span `[lo, hi)`.
fn own_indel_inside(record: &Record, (lo, hi): (i64, i64)) -> bool {
    let mut rp = record.pos();
    for op in record.cigar().iter() {
        match op {
            Cigar::Del(l) | Cigar::Ins(l) if *l >= MIN_OWN_INDEL && rp >= lo && rp < hi => return true,
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) | Cigar::Del(l) | Cigar::RefSkip(l) => rp += *l as i64,
            _ => {}
        }
    }
    false
}

/// Whether one of the read's splices (N) falls between query positions `lo` and
/// `hi`: the read's bases there are not one contiguous stretch of the genome.
fn crosses_splice(record: &Record, lo: usize, hi: usize) -> bool {
    let mut q = 0usize;
    for op in record.cigar().iter() {
        match op {
            Cigar::RefSkip(_) if q > lo && q < hi => return true,
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) | Cigar::Ins(l) | Cigar::SoftClip(l) => {
                q += *l as usize;
            }
            _ => {}
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;

    /// A deletion of `len` bases after the anchor at `pos`, prepared on `ctx`
    /// (which starts at genome position 0).
    fn del(ctx: &str, pos: i64, len: usize) -> Variant {
        let p = pos as usize;
        Variant::new(
            "1".into(), pos, ctx[p..p + 1 + len].into(), ctx[p..p + 1].into(), "DEL".into(),
            None, 0, 0, None, None, Some((0, ctx.into())), None,
        )
    }

    /// `n` C/G/T bases from a fixed seed (no A, so planted A's stand out).
    fn cgt(n: usize, seed: u32) -> String {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                ['C', 'G', 'T'][(x >> 16) as usize % 3]
            })
            .collect()
    }

    #[test]
    fn only_prepared_pure_deletions_of_fifty_bases_or_more_are_judged() {
        let ctx = cgt(300, 7);
        assert!(judges(&del(&ctx, 99, 50)));
        assert!(!judges(&del(&ctx, 99, 49)));
        let mut unprepared = del(&ctx, 99, 60);
        unprepared.event_ref = None;
        assert!(!judges(&unprepared));
    }

    #[test]
    fn the_windows_hold_the_flank_and_each_allele_s_continuation() {
        // 60 bases deleted after the anchor at 99 (deleted [100, 160)), no slide.
        let mut c: Vec<u8> = cgt(300, 11).into_bytes();
        c[100] = b'A'; // the first deleted base differs from the right flank's first
        c[99] = b'A'; // and the anchor from the last deleted base: no slide either way
        c[159] = b'C';
        c[160] = b'G';
        let ctx = String::from_utf8(c).unwrap();
        let j = junction(&del(&ctx, 99, 60)).unwrap();
        assert_eq!((j.lo, j.mh), (100, 0));
        assert_eq!(j.ref_w, &ctx.as_bytes()[96..104]);
        assert_eq!(j.alt_w, [&ctx[96..100], &ctx[160..164]].concat().as_bytes());
        assert_eq!((j.first(), j.last(false), j.last(true)), (96, 103, 163));
        assert_eq!((j.genome(true, 3), j.genome(true, 4)), (99, 160));
    }

    #[test]
    fn the_microhomology_is_read_by_both_windows() {
        // The deleted bases start with the right flank's first three: the deletion
        // slides three bases right, and both windows hold them before their own K.
        let mut c: Vec<u8> = cgt(300, 13).into_bytes();
        let (lo, len) = (100, 60);
        for i in 0..3 {
            c[lo + i] = c[lo + len + i];
        }
        c[lo + 3] = if c[lo + len + 3] == b'A' { b'C' } else { b'A' };
        c[lo - 1] = if c[lo + len - 1] == b'A' { b'C' } else { b'A' };
        let ctx = String::from_utf8(c).unwrap();
        let j = junction(&del(&ctx, lo as i64 - 1, len)).unwrap();
        assert_eq!(j.mh, 3);
        assert_eq!(j.ref_w.len(), 2 * K + 3);
        assert_eq!(&j.ref_w[K..K + 3], &j.alt_w[K..K + 3]);
        assert_ne!(j.ref_w[K + 3], j.alt_w[K + 3]);
    }

    /// `n` reads of `seq` aligned at `pos` by `cigar`, each its own fragment.
    fn reads(seq: &[u8], cigar: Vec<Cigar>, pos: i64, n: usize, tag: &str) -> Vec<Record> {
        (0..n)
            .map(|i| {
                let mut r = Record::new();
                let name = format!("{tag}{i}");
                r.set(name.as_bytes(), Some(&CigarString(cigar.clone())), seq, &vec![30; seq.len()]);
                r.set_pos(pos);
                r
            })
            .collect()
    }

    #[test]
    fn the_guard_masks_an_isolated_recurring_change_and_not_a_run() {
        // 60 bases deleted after the anchor at 99; carriers from 40 read the flank,
        // the junction and 40 right-flank bases.
        let mut c: Vec<u8> = cgt(300, 19).into_bytes();
        c[100] = b'A';
        c[99] = b'A';
        let ctx = String::from_utf8(c).unwrap();
        let v = del(&ctx, 99, 60);
        let j = junction(&v).unwrap();
        let hap = [&ctx.as_bytes()[40..100], &ctx.as_bytes()[160..200]].concat();
        let cigar = vec![Cigar::Match(60), Cigar::Del(60), Cigar::Match(40)];
        let refs = reads(&ctx.as_bytes()[40..140], vec![Cigar::Match(100)], 40, 4, "r");
        // one change two bases past the junction on every carrier: masked, named
        let mut one = hap.clone();
        one[61] = if one[61] == b'A' { b'C' } else { b'A' };
        let all = [reads(&one, cigar.clone(), 40, 4, "c"), refs.clone()].concat();
        let g = guard_for(&all, &v, &|_| true, 20).expect("a recurring change");
        assert_eq!(g.masked, vec![j.genome(false, K + 1), j.genome(true, K + 1)]);
        assert!(g.larger.is_some());
        // two inserted bases at the junction shift every base after it: a run of
        // changes (another allele), not masked
        let ins = [&hap[..60], b"TT", &hap[60..98]].concat();
        let cig = vec![Cigar::Match(60), Cigar::Del(60), Cigar::Ins(2), Cigar::Match(38)];
        let all = [reads(&ins, cig, 40, 4, "t"), refs].concat();
        assert!(guard_for(&all, &v, &|_| true, 20).is_none());
    }

    #[test]
    fn a_reference_too_short_for_the_windows_gives_none() {
        let ctx = cgt(162, 17);
        assert!(junction(&del(&ctx, 99, 60)).is_none()); // 2 right-flank bases only
    }
}
