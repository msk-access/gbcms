//! The allele the reads carry at a row, when it is not the given one.
//!
//! gbcms counts the given allele. When the input is mis-described (a curator's
//! edit, a caller merging events), those counts are honest but low, so the row
//! must say what the reads show instead: `OBSERVED_ALLELE` in `gbcms_diagnostic`.
//! Nothing here changes a count.
//!
//! Each spanning read is rebuilt over the event's core: its change interval
//! (`window::change_interval`, the shift region for pure indels) plus one base on
//! each side. A read must cover the core plus [`FLANK`] bases to be read, so it
//! is anchored. The flank itself is not compared, so a germline SNP beside the
//! event is not an allele of this row. Reads with a base below min BQ in the core
//! are skipped: the quality contract the classifiers use. Identical core
//! sequences are counted; the most frequent one that is neither REF, the given
//! ALT, nor a co-annotated sibling's ALT (already an input row) is named when
//! enough reads carry it exactly.

use std::collections::HashMap;

use rust_htslib::bam::record::{Cigar, Record};

use super::variant_checks::reconstruct_span;
use super::window;
use crate::types::Variant;

/// Reference bases a read must cover beyond the core on each side.
const FLANK: i64 = 5;
/// Exact carriers needed before an allele is named.
const MIN_CARRIERS: u32 = 3;
/// Share of the scanned (spanning) reads the observed allele must reach.
const MIN_FRACTION: f64 = 0.05;

/// A named allele: 0-based POS, trimmed VCF-style REF/ALT, and read counts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ObservedAllele {
    pub pos: i64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub carriers: u32,
    pub given_carriers: u32,
}

/// An allele in canonical form: 0-based POS, REF, ALT (left-aligned, minimal).
type Allele = (i64, Vec<u8>, Vec<u8>);

/// The most frequent allele the reads carry at the event that is neither REF,
/// the given allele, nor a co-annotated sibling's, if at least [`MIN_CARRIERS`]
/// reads carry it exactly, more than carry the given allele, and at least
/// [`MIN_FRACTION`] of the scanned reads. None without prep's `event_ref`.
pub(crate) fn observed_allele(
    read_cache: &[Record],
    variant: &Variant,
    siblings: &[Variant],
    min_mapq: u8,
    min_baseq: u8,
) -> Option<ObservedAllele> {
    let (ev_start, ev_seq) = variant.event_ref.as_ref()?;
    let ev_start = *ev_start;
    let refseq: Vec<u8> = ev_seq.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let ev_end = ev_start + refseq.len() as i64;
    let base = |g: i64| -> Option<u8> { (g >= ev_start && g < ev_end).then(|| refseq[(g - ev_start) as usize]) };

    let (c_lo, c_hi) = window::change_interval(variant);
    let core_lo = (c_lo - 1).min(variant.pos);
    let core_hi = (c_hi + 1).max(variant.pos + variant.ref_allele.len() as i64);
    if core_lo - FLANK - 1 < ev_start || core_hi + FLANK > ev_end {
        return None;
    }

    let upper = |s: &str| -> Vec<u8> { s.bytes().map(|b| b.to_ascii_uppercase()).collect() };
    let given = canonical(variant.pos, upper(&variant.ref_allele), upper(&variant.alt_allele), &base)?;
    let known: Vec<Allele> = siblings
        .iter()
        .filter_map(|s| canonical(s.pos, upper(&s.ref_allele), upper(&s.alt_allele), &base))
        .collect();

    let mut seen: HashMap<Allele, u32> = HashMap::new();
    let mut scanned: u32 = 0;
    for record in read_cache {
        if record.is_secondary() || record.is_supplementary() || record.mapq() < min_mapq {
            continue;
        }
        let read_end = window::ref_end(record);
        if record.pos() > core_lo - FLANK || read_end < core_hi + FLANK {
            continue;
        }
        // Widen the comparison over any indel of the read that could sit on the
        // core: its shift region (every equivalent placement) touches it. A longer
        // event is then read whole, and one indel aligned at different places in a
        // repeat is read the same way.
        let (mut lo, mut hi) = (core_lo, core_hi);
        let (mut ref_pos, mut read_pos) = (record.pos(), 0usize);
        let seq = record.seq();
        for op in record.cigar().iter() {
            match op {
                Cigar::Del(len) => {
                    let d_end = ref_pos + *len as i64;
                    let anchor = base(ref_pos - 1);
                    let deleted: Option<Vec<u8>> = (ref_pos..d_end).map(&base).collect();
                    if let (Some(anc), Some(del)) = (anchor, deleted) {
                        let mut r = vec![anc];
                        r.extend(del);
                        let (r0, r1) = window::shift_region_over(
                            ref_pos - 1,
                            &String::from_utf8_lossy(&r),
                            &String::from_utf8_lossy(&[anc]),
                            base,
                        );
                        if r0 <= hi && r1 >= lo {
                            lo = lo.min(r0);
                            hi = hi.max(r1);
                        }
                    }
                    ref_pos = d_end;
                }
                Cigar::Ins(len) => {
                    let ins: Vec<u8> = (read_pos..read_pos + *len as usize)
                        .map(|i| seq[i].to_ascii_uppercase())
                        .collect();
                    if let Some(anc) = base(ref_pos - 1) {
                        let mut a = vec![anc];
                        a.extend(ins);
                        let (r0, r1) = window::shift_region_over(
                            ref_pos - 1,
                            &String::from_utf8_lossy(&[anc]),
                            &String::from_utf8_lossy(&a),
                            base,
                        );
                        if r0 <= hi && r1 >= lo {
                            lo = lo.min(r0);
                            hi = hi.max(r1);
                        }
                    }
                    read_pos += *len as usize;
                }
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    ref_pos += *len as i64;
                    read_pos += *len as usize;
                }
                Cigar::RefSkip(len) => ref_pos += *len as i64,
                Cigar::SoftClip(len) => read_pos += *len as usize,
                _ => {}
            }
        }
        // One reference base before the stretch anchors the allele (VCF form).
        let lo = lo - 1;
        if lo < ev_start || hi > ev_end || record.pos() > lo || read_end < hi {
            continue;
        }
        let recon = reconstruct_span(record, record.qual(), lo, hi);
        if recon.splice_skip || recon.quals.iter().any(|&q| q < min_baseq) {
            continue;
        }
        scanned += 1;
        let ref_part: Vec<u8> = refseq[(lo - ev_start) as usize..(hi - ev_start) as usize].to_vec();
        let read_part: Vec<u8> = recon.seq.iter().map(u8::to_ascii_uppercase).collect();
        if let Some(allele) = canonical(lo, ref_part, read_part, &base) {
            *seen.entry(allele).or_insert(0) += 1;
        }
    }

    let given_n = seen.get(&given).copied().unwrap_or(0);
    // Most carriers first; ties broken by allele so the choice is stable.
    let (best, &n) = seen
        .iter()
        .filter(|(a, _)| **a != given && !known.contains(*a))
        .max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0)))?;
    if n < MIN_CARRIERS || n <= given_n || (n as f64) < MIN_FRACTION * scanned as f64 {
        return None;
    }
    Some(ObservedAllele {
        pos: best.0,
        ref_allele: String::from_utf8_lossy(&best.1).into_owned(),
        alt_allele: String::from_utf8_lossy(&best.2).into_owned(),
        carriers: n,
        given_carriers: given_n,
    })
}

/// Canonical (left-aligned, minimal VCF) form of `REF>ALT` at 0-based `pos`,
/// or None when the alleles are equal (no change) or left-alignment runs off
/// the known reference. The standard algorithm: while the last bases agree,
/// drop them; when an allele empties, extend both one reference base left;
/// then drop shared leading bases, keeping one.
fn canonical(pos: i64, r: Vec<u8>, a: Vec<u8>, base: &dyn Fn(i64) -> Option<u8>) -> Option<Allele> {
    if r == a {
        return None;
    }
    let (mut pos, mut r, mut a) = (pos, r, a);
    loop {
        if !r.is_empty() && !a.is_empty() && r.last() == a.last() {
            r.pop();
            a.pop();
        } else if r.is_empty() || a.is_empty() {
            let b = base(pos - 1)?;
            r.insert(0, b);
            a.insert(0, b);
            pos -= 1;
        } else {
            break;
        }
    }
    while r.len() > 1 && a.len() > 1 && r[0] == a[0] {
        r.remove(0);
        a.remove(0);
        pos += 1;
    }
    Some((pos, r, a))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn at(seq: &'static [u8], start: i64) -> impl Fn(i64) -> Option<u8> {
        move |g| (g >= start && ((g - start) as usize) < seq.len()).then(|| seq[(g - start) as usize])
    }

    #[test]
    fn a_substitution_is_trimmed_to_its_base() {
        let b = at(b"ACAGT", 10);
        assert_eq!(canonical(10, b"ACA".to_vec(), b"AGA".to_vec(), &b), Some((11, b"C".to_vec(), b"G".to_vec())));
    }

    #[test]
    fn a_homopolymer_deletion_left_aligns_to_the_run_start() {
        // G CCCCCC T: deleting the last C is the same allele as deleting the first.
        let b = at(b"AGCCCCCCTA", 0);
        assert_eq!(canonical(6, b"CC".to_vec(), b"C".to_vec(), &b), Some((1, b"GC".to_vec(), b"G".to_vec())));
    }

    #[test]
    fn a_delins_keeps_its_own_bases() {
        // G CCCCCC T -> G CCCCT T: CC>T at the run's end.
        let b = at(b"AGCCCCCCTA", 0);
        assert_eq!(
            canonical(1, b"GCCCCCCT".to_vec(), b"GCCCCTT".to_vec(), &b),
            Some((6, b"CC".to_vec(), b"T".to_vec()))
        );
    }

    #[test]
    fn equal_alleles_are_no_change() {
        let b = at(b"ACGT", 0);
        assert_eq!(canonical(0, b"ACG".to_vec(), b"ACG".to_vec(), &b), None);
    }
}
