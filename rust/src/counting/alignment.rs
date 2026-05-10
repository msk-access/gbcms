//! Smith-Waterman alignment backend for Phase 3 read classification.
//!
//! Contains `classify_by_alignment` (dual-haplotype semiglobal SW with
//! local fallback for complex variants), and read extraction helpers
//! used by Phase 0 structural anomaly detection.
//!
//! See also `pairhmm.rs` for the probabilistic PairHMM alternative
//! backend, selectable via `--alignment-backend hmm`.

use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;
use bio::alignment::pairwise::Aligner;
use log::{debug, trace};

use crate::types::Variant;
use super::utils::{median_qual, build_haplotypes, ClassifyResult, ClassifyPhase};

/// Extract contiguous **raw** read bases spanning a genomic window `[win_start, win_end)`.
///
/// Unlike `extract_read_subsequence` which concatenates CIGAR-projected bases
/// (breaking complex variants represented as DEL+INS), this function finds the
/// read position range overlapping the window and returns the **original read
/// sequence** — preserving the true underlying allele for SW alignment.
///
/// ## Why this matters for complex variants
///
/// For EPHA7 (REF=TCC, ALT=CT), BWA represents ALT reads as `91M2D6M4I`.
/// CIGAR-projected extraction produces `GGAAACTCTCCAAAA` which matches
/// neither REF (`GGAAATCCACTCC`) nor ALT (`GGAAACTACTCC`).
/// Raw extraction produces `GGAAACTCTCCAAAA` from the same read positions,
/// which SW alignment can correctly classify by scoring against both haplotypes.
///
/// ## Algorithm
///
/// 1. Walk CIGAR to find the first read_pos where ref enters the window
/// 2. Walk CIGAR to find the last read_pos where ref is still in the window
///    (including insertion bases at window boundaries)
/// 3. Return `seq[first_read_pos..=last_read_pos]`
///
/// Returns `(bases, quals)` — empty if the read doesn't overlap the window.
///
/// `variant_pos` and `variant_ref_len` restrict soft-clip inclusion to clips
/// adjacent to the variant site, preventing adapter/chimeric garbage from
/// large leading clips from polluting Phase 3 alignment.
pub fn extract_raw_read_window(
    record: &Record,
    quals: &[u8],
    win_start: i64,
    win_end: i64,
    variant_pos: i64,
    variant_ref_len: usize,
) -> Option<(Vec<u8>, Vec<u8>)> {
    let cigar = record.cigar();
    let mut ref_pos = record.pos();
    let mut read_pos: usize = 0;
    let seq = record.seq();
    // NOTE: quals is passed from the caller — either raw record.qual() or BAQ-adjusted.

    let mut first_read_pos: Option<usize> = None;
    let mut last_read_pos: Option<usize> = None;

    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len_i64 = *len as i64;
                let len_usize = *len as usize;

                // Check overlap with [win_start, win_end)
                let overlap_start = std::cmp::max(ref_pos, win_start);
                let overlap_end = std::cmp::min(ref_pos + len_i64, win_end);
                if overlap_start < overlap_end {
                    let offset = (overlap_start - ref_pos) as usize;
                    let count = (overlap_end - overlap_start) as usize;
                    let rp_start = read_pos + offset;
                    let rp_end = rp_start + count - 1;
                    if first_read_pos.is_none() {
                        first_read_pos = Some(rp_start);
                    }
                    last_read_pos = Some(rp_end);
                }
                ref_pos += len_i64;
                read_pos += len_usize;
            }
            Cigar::Ins(len) => {
                let len_usize = *len as usize;
                // Include insertion bases if they fall within the window boundary.
                // An insertion at ref_pos X appears between ref X-1 and X.
                if ref_pos >= win_start && ref_pos <= win_end {
                    let rp_end = read_pos + len_usize - 1;
                    if first_read_pos.is_none() {
                        first_read_pos = Some(read_pos);
                    }
                    last_read_pos = Some(rp_end);
                }
                read_pos += len_usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            }
            Cigar::SoftClip(len) => {
                let len_usize = *len as usize;
                // Only include soft clips adjacent to the variant position
                // itself (within ±1bp), not the broad fetching window.
                let var_end = variant_pos + variant_ref_len as i64;
                if ref_pos >= variant_pos - 1 && ref_pos <= var_end + 1 {
                    let rp_end = read_pos + len_usize - 1;
                    if first_read_pos.is_none() {
                        first_read_pos = Some(read_pos);
                    }
                    last_read_pos = Some(rp_end);
                }
                read_pos += len_usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    // Extract contiguous raw read bases from first to last overlapping position.
    // Returns None instead of empty Vecs to avoid unnecessary heap allocation.
    match (first_read_pos, last_read_pos) {
        (Some(first), Some(last)) if first <= last && last < seq.len() => {
            let window_len = last - first + 1;
            let mut bases = Vec::with_capacity(window_len);
            for i in first..=last {
                bases.push(seq[i]);
            }
            let base_quals = quals[first..=last].to_vec();
            Some((bases, base_quals))
        }
        _ => None,
    }
}

/// Quick pre-filter before Smith-Waterman alignment (indelpost pattern).
///
/// Returns `true` only if the read shows evidence that could change the
/// allele call via SW realignment:
/// - Soft-clipping within or near the variant window
/// - CIGAR contains indels within `[win_start, win_end)`
///
/// Clean M-only reads that span the variant are skipped (the vast majority
/// at any locus), eliminating ~80-90% of unnecessary SW calls.
pub fn is_worth_realignment(record: &Record, win_start: i64, win_end: i64) -> bool {
    let cigar = record.cigar();
    let mut ref_pos = record.pos();

    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                ref_pos += *len as i64;
            }
            Cigar::Ins(_) => {
                // Insertion at this ref_pos — check if within window
                if ref_pos >= win_start && ref_pos <= win_end {
                    return true;
                }
            }
            Cigar::Del(len) => {
                let del_end = ref_pos + *len as i64;
                // Deletion overlapping window
                if ref_pos < win_end && del_end > win_start {
                    return true;
                }
                ref_pos = del_end;
            }
            Cigar::SoftClip(_) => {
                // Soft-clip near the variant window suggests misalignment
                if ref_pos >= win_start - 5 && ref_pos <= win_end + 5 {
                    return true;
                }
            }
            Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            }
            _ => {}
        }
    }
    false
}

/// Classify a read subsequence as REF or ALT using dual-haplotype
/// Smith-Waterman alignment (inspired by indelpost).
///
/// Builds ref/alt haplotypes from `variant.ref_context`, then aligns
/// `read_seq` against both using the provided reusable `Aligner` instances.
///
/// **Alignment mode**: Semiglobal — haplotype (query) must be fully aligned,
/// read (text) has free overhangs. This answers "does this read *contain*
/// this haplotype?" — the standard approach for pattern-in-text problems.
///
/// **Haplotype trimming**: When a haplotype exceeds the read length (common
/// with large `context_padding` and complex CIGARs), the haplotype is
/// symmetrically trimmed from flanks to fit. This ensures semiglobal alignment
/// penalizes only allele differences, not excess length. A fixed score margin
/// of 2 works reliably because shared flanking bases contribute equally to
/// both allele scores — trimming preserves this invariant.
///
/// **Memory optimization**: Aligners are created once per variant in
/// `count_single_variant()` and reused for all reads, avoiding repeated
/// O(n×m) DP matrix allocation (indelpost pattern).
///
/// Returns `(is_ref, is_alt, base_qual)` where `base_qual` is the median
/// quality across the read subsequence bases (GATK standard).
///
/// # Alignment semantics
/// Uses semiglobal(read, haplotype): the read is the pattern (globally consumed,
/// no free overhangs) and the haplotype is the text (free leading/trailing gaps).
/// This allows the aligner to find the best-scoring window within a longer
/// haplotype that matches the read — critical for large deletions where the
/// REF haplotype (e.g. 501bp for a 500bp deletion) far exceeds the read length
/// (~150bp), and for off-center reads that don't symmetrically overlap the variant.
///
/// See also [`classify_by_marginalized_pairhmm`](super::pairhmm::classify_by_marginalized_pairhmm) for a
/// probabilistic alternative using BQ-aware emission probabilities.
pub fn classify_by_alignment<F: Fn(u8, u8) -> i32>(
    read_seq: &[u8],
    read_quals: &[u8],
    variant: &Variant,
    min_baseq: u8,
    alt_aligner: &mut Aligner<F>,
    ref_aligner: &mut Aligner<F>,
) -> ClassifyResult {
    let (ref_hap, alt_hap) = match build_haplotypes(variant) {
        Some(haps) => haps,
        None => {
            debug!("classify_by_alignment: build_haplotypes failed for {}:{}",
                   variant.chrom, variant.pos + 1);
            return ClassifyResult::neither(ClassifyPhase::Alignment);
        }
    };

    // Mask low-quality bases as N so they don't bias scoring.
    let masked_seq: Vec<u8> = read_seq
        .iter()
        .zip(read_quals.iter())
        .map(|(&b, &q)| if q >= min_baseq { b } else { b'N' })
        .collect();

    // Skip if too few usable bases.
    let read_len = read_seq.len();
    let usable_count = read_quals.iter().filter(|&&q| q >= min_baseq).count();
    if usable_count < 3 {
        trace!("classify_by_alignment: only {} usable bases — skipping", usable_count);
        return ClassifyResult::neither(ClassifyPhase::Alignment);
    }

    // Use the provided reusable aligners (created once per variant).
    trace!(
        "classify_by_alignment: read_len={} ref_hap_len={} alt_hap_len={} usable={}",
        read_len, ref_hap.len(), alt_hap.len(), usable_count
    );
    trace!(
        "classify_by_alignment seqs: read={} alt_hap={} ref_hap={}",
        String::from_utf8_lossy(&masked_seq),
        String::from_utf8_lossy(&alt_hap),
        String::from_utf8_lossy(&ref_hap)
    );

    // Semiglobal alignment: read = query (fully consumed, every base must
    // participate), haplotype = text (free leading/trailing overhangs).
    // The aligner slides the read along the haplotype to find the best window.
    // This correctly handles:
    //  - Large deletions: 150bp read vs 501bp REF haplotype
    //  - Off-center reads: no symmetric trimming needed
    //  - Short ALT haplotypes: read may overhang, but all read bases align
    let alt_aln = alt_aligner.semiglobal(&masked_seq, &alt_hap);
    let ref_aln = ref_aligner.semiglobal(&masked_seq, &ref_hap);

    let med_qual = median_qual(read_quals, min_baseq);

    // Score margin: require at least 2 points difference to call.
    // Uses >= (not >) because with larger context_padding, the SW aligner
    // can find shifted alignments that reduce the score difference from 3
    // to exactly 2. Strict > would incorrectly classify these as ambiguous.
    let margin = 2;
    let mut is_alt = alt_aln.score >= ref_aln.score + margin;
    let mut is_ref = ref_aln.score >= alt_aln.score + margin;

    trace!(
        "classify_by_alignment: alt_score={} ref_score={} margin={} is_ref={} is_alt={} \
         read_len={} ref_hap={} alt_hap={}",
        alt_aln.score, ref_aln.score, margin, is_ref, is_alt,
        read_len, ref_hap.len(), alt_hap.len()
    );

    // Dual-trigger local fallback for indel/complex variants.
    //
    // When the MAF/VCF definition is incomplete (e.g., TCC→CT missing an
    // adjacent SNV), the ALT haplotype has a "frameshifted flank": bases in
    // the right context that don't match the biological read. Semiglobal
    // alignment forces gap penalties through this invalid flank, producing
    // confident but WRONG calls (e.g., EPHA7: alt=-4 ref=-2 → REF).
    //
    // Local alignment soft-clips the bad flank, finding the best matching
    // substring without penalizing overhangs on either side.
    //
    // Two triggers detect low-confidence semiglobal results:
    //  - Borderline: score difference is within margin+1 (barely decisive)
    //  - Poor quality: best score is < half the "perfect" score (both
    //    haplotypes are fighting severe penalties)
    //
    // Only applied to **complex** variants (both alleles > 1bp), where the
    // MAF definition is most likely to be incomplete. Pure insertions and
    // deletions are well-handled by semiglobal alignment.
    let ref_allele_len = variant.ref_allele.len();
    let alt_allele_len = variant.alt_allele.len();
    let is_complex = ref_allele_len > 1 && alt_allele_len > 1
        && ref_allele_len != alt_allele_len;

    if is_complex {
        let diff = (alt_aln.score - ref_aln.score).abs();
        let max_score = std::cmp::max(alt_aln.score, ref_aln.score);
        let perfect_score = read_len as i32;
        let is_borderline = diff <= margin + 1;
        let is_poor = max_score < perfect_score / 2;

        if is_borderline || is_poor {
            trace!(
                "Phase 3 low-confidence (diff={}, max_score={}/{}, borderline={}, poor={}) \
                 → local fallback",
                diff, max_score, perfect_score, is_borderline, is_poor
            );
            let alt_local = alt_aligner.local(&masked_seq, &alt_hap);
            let ref_local = ref_aligner.local(&masked_seq, &ref_hap);

            is_alt = alt_local.score >= ref_local.score + margin;
            is_ref = ref_local.score >= alt_local.score + margin;

            trace!(
                "Phase 3 local result: alt_local={} ref_local={} is_alt={} is_ref={}",
                alt_local.score, ref_local.score, is_alt, is_ref
            );
        }
    }

    if is_alt {
        ClassifyResult::is_alt(med_qual, ClassifyPhase::Alignment)
    } else if is_ref {
        let mut r = ClassifyResult::is_ref(med_qual, ClassifyPhase::Alignment);
        // Propagate nearby evidence: read was classified as REF, but if ALT
        // alignment score was close (within 2× margin), there is structural
        // ALT support worth tracking for partial_alt/PARTIAL_DOMINANT diagnostic.
        if alt_aln.score > 0 && alt_aln.score + 2 * margin >= ref_aln.score {
            r.has_nearby_evidence = true;
            trace!(
                "Phase 3 SW: REF wins but ALT score close (alt={}, ref={}, margin={}) → has_nearby_evidence=true",
                alt_aln.score, ref_aln.score, margin
            );
        }
        r
    } else {
        // Tie or ambiguous: the read cannot confidently distinguish REF from
        // ALT.  Route to "neither" so this read still contributes to physical
        // DP (via the anchor overlap gate) but does NOT inflate RD or AD.
        // This preserves an unbiased VAF = AD/(RD+AD) — critical for low-VAF
        // cfDNA detection where even a handful of misrouted reads can push
        // a variant below the Limit of Detection.
        let max_score = std::cmp::max(alt_aln.score, ref_aln.score);
        if max_score >= (read_len as i32) / 2 {
            trace!("Ambiguous tie (alt={}, ref={}) — routing to neither to preserve unbiased VAF", alt_aln.score, ref_aln.score);
            let mut r = ClassifyResult::new(false, false, med_qual, ClassifyPhase::Alignment);
            // Propagate nearby evidence: tie means ALT had non-trivial
            // alignment support — worth tracking for partial_alt.
            if alt_aln.score > 0 {
                r.has_nearby_evidence = true;
                trace!(
                    "Phase 3 SW: tie with ALT score {} > 0 → has_nearby_evidence=true",
                    alt_aln.score
                );
            }
            r
        } else {
            ClassifyResult::neither(ClassifyPhase::Alignment)
        }
    }
}
