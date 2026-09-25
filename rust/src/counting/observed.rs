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

use rust_htslib::bam::Record;

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

/// The most frequent non-REF, non-given allele the reads carry over the
/// event's core, if at least [`MIN_CARRIERS`] reads carry it exactly, more
/// than carry the given ALT, and at least [`MIN_FRACTION`] of the scanned
/// reads. None without prep's `event_ref`.
pub(crate) fn observed_allele(
    read_cache: &[Record],
    variant: &Variant,
    siblings: &[Variant],
    min_mapq: u8,
    min_baseq: u8,
) -> Option<ObservedAllele> {
    // The core and its reference bases, as prep fetched them (every variant type).
    let (core_lo, core_seq) = variant.event_ref.as_ref()?;
    let core_lo = *core_lo;
    let ref_core: Vec<u8> = core_seq.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let core_hi = core_lo + ref_core.len() as i64;
    if variant.pos < core_lo || variant.pos + variant.ref_allele.len() as i64 > core_hi {
        return None;
    }
    let alt_core = apply_to_core(&ref_core, core_lo, variant)?;
    // Alleles already in the input: co-annotated siblings that fall inside the core.
    let known: Vec<Vec<u8>> = siblings.iter().filter_map(|s| apply_to_core(&ref_core, core_lo, s)).collect();

    let (a_lo, a_hi) = (core_lo - FLANK, core_hi + FLANK);
    let mut seen: HashMap<Vec<u8>, u32> = HashMap::new();
    let mut scanned: u32 = 0;
    for record in read_cache {
        if record.is_secondary() || record.is_supplementary() || record.mapq() < min_mapq {
            continue;
        }
        if record.pos() > a_lo || window::ref_end(record) < a_hi {
            continue;
        }
        let recon = reconstruct_span(record, record.qual(), core_lo, core_hi);
        if recon.splice_skip || recon.quals.iter().any(|&q| q < min_baseq) {
            continue;
        }
        scanned += 1;
        let key: Vec<u8> = recon.seq.iter().map(u8::to_ascii_uppercase).collect();
        *seen.entry(key).or_insert(0) += 1;
    }

    let given = seen.get(&alt_core).copied().unwrap_or(0);
    // Most carriers first; ties broken by sequence so the choice is stable.
    let (best, &n) = seen
        .iter()
        .filter(|(s, _)| **s != ref_core && **s != alt_core && !known.contains(*s))
        .max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0)))?;
    if n < MIN_CARRIERS || n <= given || (n as f64) < MIN_FRACTION * scanned as f64 {
        return None;
    }
    let (pos, r, a) = trim_to_vcf(core_lo, &ref_core, best);
    Some(ObservedAllele {
        pos,
        ref_allele: String::from_utf8_lossy(&r).into_owned(),
        alt_allele: String::from_utf8_lossy(&a).into_owned(),
        carriers: n,
        given_carriers: given,
    })
}

/// The core with `v`'s ALT applied, or None when `v`'s REF does not lie inside
/// the core or does not match it.
fn apply_to_core(ref_core: &[u8], core_lo: i64, v: &Variant) -> Option<Vec<u8>> {
    let start = v.pos - core_lo;
    let end = start + v.ref_allele.len() as i64;
    if start < 0 || end as usize > ref_core.len() {
        return None;
    }
    let (start, end) = (start as usize, end as usize);
    if !ref_core[start..end].eq_ignore_ascii_case(v.ref_allele.as_bytes()) {
        return None;
    }
    let mut out = ref_core[..start].to_vec();
    out.extend(v.alt_allele.bytes().map(|b| b.to_ascii_uppercase()));
    out.extend_from_slice(&ref_core[end..]);
    Some(out)
}

/// Trim shared trailing then leading bases, keeping one base in each allele
/// (VCF form). Returns the 0-based position of the first kept base.
fn trim_to_vcf(start: i64, r: &[u8], a: &[u8]) -> (i64, Vec<u8>, Vec<u8>) {
    let (mut r, mut a) = (r.to_vec(), a.to_vec());
    while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
        r.pop();
        a.pop();
    }
    let mut k = 0;
    while k + 1 < r.len() && k + 1 < a.len() && r[k] == a[k] {
        k += 1;
    }
    (start + k as i64, r[k..].to_vec(), a[k..].to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trimming_keeps_one_anchor_base() {
        // X CCCCCC T -> X CCCCT T: the change is CC>T at offset 5.
        assert_eq!(trim_to_vcf(10, b"GCCCCCCT", b"GCCCCTT"), (15, b"CC".to_vec(), b"T".to_vec()));
        // A pure deletion keeps its anchor.
        assert_eq!(trim_to_vcf(0, b"ACGT", b"AT"), (0, b"ACG".to_vec(), b"A".to_vec()));
        // A substitution.
        assert_eq!(trim_to_vcf(4, b"ACA", b"AGA"), (5, b"C".to_vec(), b"G".to_vec()));
    }
}
