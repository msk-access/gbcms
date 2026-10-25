//! A read's call from its scores against the site's haplotypes, and the call mapped onto
//! the engine's `ClassifyResult`.
//!
//! With H_alt the haplotypes carrying X at its site, H_ref those with the reference's
//! bases there and H_other those with another allele: ALT where the best of H_alt beats
//! the best of the rest by [`DELTA`] (one confident base), REF and OTHER likewise, else
//! depth only. No priors; a tie is never credited.

use log::trace;
use rust_htslib::bam::record::{Cigar, Record};

use super::haplotypes::{other_option, Class, ReadView, Scorer, Site};
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
}

/// Measurement only, TEMPORARY: what a delins call needs read (`GBCMS_JUDGE_DELINS`):
/// `asym`, an ALT call X's telling bases and a REF call none; `sym`, each call its own
/// allele's; `margin`, neither (the margin alone). Default `margin`.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum DelinsOption {
    Asym,
    Sym,
    Margin,
}

pub(crate) fn delins_option() -> DelinsOption {
    static OPT: std::sync::OnceLock<DelinsOption> = std::sync::OnceLock::new();
    *OPT.get_or_init(|| match std::env::var("GBCMS_JUDGE_DELINS").as_deref() {
        Ok("asym") => DelinsOption::Asym,
        Ok("sym") => DelinsOption::Sym,
        _ => DelinsOption::Margin,
    })
}

/// A read whose best score on every haplotype is below this (ten confident mismatches)
/// fits none of them: depth only.
pub(crate) const FITS_NOTHING: i32 = 10 * crate::counting::judge::score::MISMATCH;

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
    for (tied, beaten) in [(0, 1), (1, 0)] {
        if best[2].is_some() && best[tied].is_some() && v[2].min(v[tied]) - v[beaten] >= delta {
            return Call::Other;
        }
    }
    Call::Depth
}

/// Judge one read against the site's haplotypes. A REF or ALT call stands only where
/// the read's alignment to the winning haplotype reads the bases that tell X from REF
/// at X's site: an ALT call (X a substitution or a delins) every one of X's own,
/// confidently and as X's; a REF call on a substitution every one of REF's, covered (as
/// the size-specific rules require of an MNP's block). A read that shows part of them is
/// depth only. With `ask_anchor` (a read the genome places
/// off X's anchor), the verdict says whether its alignment to the winning haplotype
/// holds the anchor there.
pub(crate) fn judge(site: &Site, read: &ReadView, min_baseq: u8, ask_anchor: bool) -> Verdict {
    let mut scorer = Scorer::new(read);
    // Per class, the best score and the first haplotype reaching it. A haplotype whose
    // banded score cannot stand is scored in full only where an upper bound on its score
    // could still raise its class's best (or reach it first), so a read fitting none of
    // the haplotypes well (a dense cluster) does not run the full programme on each.
    let mut best: [Option<(i32, usize)>; 3] = [None; 3];
    let raise = |best: &mut [Option<(i32, usize)>; 3], class: usize, s: i32, h: usize| {
        let slot = &mut best[class];
        if slot.is_none_or(|(b, bh)| s > b || (s == b && h < bh)) {
            *slot = Some((s, h));
        }
    };
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
        return Verdict { call: Call::Depth, qual: 0, best: [None; 3], holds_anchor: false, placed: None };
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
    let scores = best.map(|b| b.map(|(s, _)| s));
    let mut call = decide(scores);
    // The haplotype the read fits best of all, where its alignment there leaves an indel
    // unexplained inside X's site: the read carries another allele there, never REF or X.
    let overall = best.iter().flatten().copied().max_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)));
    let mut fit_other = None;
    if other_option().fit {
        if let Some((s, h)) = overall.filter(|&(s, _)| s < 0) {
            let (al, off) = scorer.align(site, &site.haps[h]);
            if site.gap_at_site(&site.haps[h], &read.seq, &al, off) {
                trace!("haplotype judge: the best haplotype ({s}) leaves an indel at X's site unexplained → another allele");
                call = Call::Other;
                fit_other = Some(h);
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
        return Verdict { call, qual: 0, best: scores, holds_anchor: false, placed: None };
    };
    let w = fit_other.unwrap_or_else(|| best[won as usize].map_or(0, |(_, h)| h));
    let hap = &site.haps[w];
    if call == Call::Other && !ask_anchor {
        return Verdict { call, qual: 0, best: scores, holds_anchor: false, placed: None };
    }
    let (wa, off) = scorer.align(site, hap);
    let column = |p: &Placed| if let Placed::At(c) = *p { Some((c + off) as i64) } else { None };
    let holds_anchor = wa.covers().is_some_and(|(a, b)| {
        let anchor = site.anchor_in(hap);
        ((a + off) as i64) <= anchor && anchor <= (b + off) as i64
    });
    let at = |p: &Placed| matches!(p, Placed::At(_));
    let placed = wa.placed.iter().position(at).zip(wa.placed.iter().rposition(at)).map(|(a, b)| (a, b + 1));
    if call == Call::Other {
        return Verdict { call, qual: 0, best: scores, holds_anchor, placed };
    }
    // The bases that tell X from REF at its site, read by the winning alignment: X's own
    // for an ALT call; REF's for a REF call on a substitution (the rules read an MNP's
    // whole block), never on a delins (the rules read equal-length windows there).
    let need = site.is_substitution()
        || match delins_option() {
            DelinsOption::Asym => call == Call::Alt,
            DelinsOption::Sym => true,
            DelinsOption::Margin => false,
        };
    let telling = if need { site.telling_in(hap) } else { Vec::new() };
    let reads_all = telling.into_iter().all(|(col, base)| {
        wa.placed.iter().enumerate().any(|(i, p)| {
            column(p) == Some(col) && (call == Call::Ref || (read.seq[i] != WILD && read.seq[i] == base))
        })
    });
    if !reads_all {
        call = Call::Depth;
        return Verdict { call, qual: 0, best: scores, holds_anchor: false, placed: None };
    }
    // The best losing haplotype: the best of the other classes, the first on a tie.
    let floor = min_baseq.max(1);
    let lost = (0..3)
        .filter(|&k| k != won as usize)
        .filter_map(|k| best[k])
        .max_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)))
        .map(|(_, h)| h);
    let qual = lost.map_or(floor, |l| {
        let (la, _) = scorer.align(site, &site.haps[l]);
        deciding_quality(read, &wa, &la).map_or(floor, |q| q.max(floor))
    });
    Verdict { call, qual, best: scores, holds_anchor, placed }
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
    ReadView::new(&seq, quals, record.pos(), &cigar, readable, min_baseq).map(|v| v.named(record.qname()))
}

/// Measurement only, TEMPORARY, to be removed before the merge (design §10): with
/// `GBCMS_JUDGE_STRUCTURAL=keep` in the environment, an ALT call whose CIGAR writes the
/// row's exact indel at its position keeps `is_structural`, so the arms can measure both
/// ways of resolving mates. Read at each such call, so a process may switch it.
fn keep_structural() -> bool {
    std::env::var("GBCMS_JUDGE_STRUCTURAL").is_ok_and(|v| v == "keep")
}

/// Whether the read's CIGAR writes `v`'s pure insertion or deletion exactly: the op,
/// its length (and an insertion's bases) right after the anchor.
fn writes_exact_op(record: &Record, v: &Variant) -> bool {
    let kind = allele_kind(&v.ref_allele, &v.alt_allele);
    let (mut rp, mut qp) = (record.pos(), 0usize);
    let after_anchor = v.pos + 1;
    for op in record.cigar().iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                rp += n as i64;
                qp += n as usize;
            }
            Cigar::Ins(n) => {
                if kind == Some(AlleleKind::Insertion) && rp == after_anchor && n as usize + 1 == v.alt_allele.len() {
                    let seq = record.seq().as_bytes();
                    return seq.get(qp..qp + n as usize).is_some_and(|b| b.eq_ignore_ascii_case(&v.alt_allele.as_bytes()[1..]));
                }
                qp += n as usize;
            }
            Cigar::Del(n) => {
                if kind == Some(AlleleKind::Deletion) && rp == after_anchor && n as usize + 1 == v.ref_allele.len() {
                    return true;
                }
                rp += n as i64;
            }
            Cigar::RefSkip(n) => rp += n as i64,
            Cigar::SoftClip(n) => qp += n as usize,
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    false
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
        Some(read) => {
            let verdict = judge(site, &read, min_baseq, !overlaps_anchor);
            let mut out = match verdict.call {
                Call::Alt => ClassifyResult::is_alt(verdict.qual, PHASE),
                Call::Ref => ClassifyResult::is_ref(verdict.qual, PHASE),
                Call::Other => ClassifyResult::neither_with_partial(PHASE, 1, false),
                Call::Depth => depth(),
            };
            // Admission off the anchor, by the alignment to the called haplotype. A clip
            // admits the read only where every base that alignment places lies inside the
            // read's fragment (RJ-13): past the fragment end a clip is adapter, and a read
            // with no fragment (unpaired, or its pair not facing inward) has none to bound
            // it.
            if !overlaps_anchor && verdict.call != Call::Depth && verdict.holds_anchor {
                let (lead, tail) = soft_clips(record);
                if lead + tail > 0 {
                    let first = readable_range(record, rules).map_or(0, |r| r.0);
                    out.clip_admissible = verdict.placed.zip(fragment_query_span(record)).is_some_and(|((a, b), (lo, hi))| {
                        first + a >= lo && first + b <= hi
                    });
                } else {
                    out.read_from_right = true;
                }
            }
            trace!(
                "haplotype judge {}:{} {}>{} read={}: {:?} (best alt={:?} ref={:?} other={:?}) qual={} \
                 holds the anchor={} over {} haplotypes",
                variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele,
                String::from_utf8_lossy(record.qname()), verdict.call, verdict.best[0], verdict.best[1],
                verdict.best[2], out.qual, overlaps_anchor || verdict.holds_anchor, site.haps.len(),
            );
            out
        }
    };
    let out = mark(out, record, variant);
    haplotypes.keep_verdict(key, out);
    out
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
    // An MNP's ALT call already read every base where the alleles differ, confidently
    // and as the ALT's, on the alignment it was called on (however the CIGAR wrote them).
    out.mnp_confirmed = out.is_alt && kind == Some(AlleleKind::Mnp);
    if out.is_alt && matches!(kind, Some(AlleleKind::Insertion | AlleleKind::Deletion)) && keep_structural() {
        out.is_structural = writes_exact_op(record, variant);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::counting::judge::haplotypes::tests::{contig, random_seq, row, site};
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
    fn the_cigar_writes_the_exact_indel_or_not() {
        let c = contig(71, 0, b"");
        let ins = row(300, &(c[299] as char).to_string(), &format!("{}CG", c[299] as char));
        let v = variant(&ins);
        let seq = |ins: &[u8]| [&c[250..300], ins, &c[300..348]].concat();
        let exact = record(&seq(b"CG"), 250, vec![Cigar::Match(50), Cigar::Ins(2), Cigar::Match(48)]);
        let other_bases = record(&seq(b"TT"), 250, vec![Cigar::Match(50), Cigar::Ins(2), Cigar::Match(48)]);
        let elsewhere = record(&seq(b"CG"), 250, vec![Cigar::Match(49), Cigar::Ins(2), Cigar::Match(49)]);
        assert!(writes_exact_op(&exact, &v));
        assert!(!writes_exact_op(&other_bases, &v));
        assert!(!writes_exact_op(&elsewhere, &v));
        let del = Row::new(299, &c[299..302], &c[299..300]);
        let v = variant(&del);
        let rec = record(&[&c[250..300], &c[302..352]].concat(), 250, vec![Cigar::Match(50), Cigar::Del(2), Cigar::Match(50)]);
        assert!(writes_exact_op(&rec, &v));
    }
}
