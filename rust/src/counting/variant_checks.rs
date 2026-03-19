//! Per-variant-type classification functions.
//!
//! Contains the variant-type-specific allele checkers:
//! - `check_snp` — single nucleotide polymorphism
//! - `check_mnp` — multi-nucleotide polymorphism (with `MnpResult`)
//! - `check_insertion` — insertion with windowed CIGAR scan
//! - `check_deletion` — deletion with windowed CIGAR scan
//! - `check_complex` — complex variant (indel + substitution) with
//!   phased Masked Comparison → Levenshtein → SW alignment pipeline

use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;
use bio::alignment::distance::levenshtein;
use bio::alignment::pairwise::Aligner;
use log::{debug, trace, warn};

use crate::types::Variant;
use super::alignment::{classify_by_alignment, extract_raw_read_window, is_worth_realignment};
// INTENTIONAL: classify_by_pairhmm retained for concordance testing against
// the new marginalized pipeline. Remove after concordance proven.
use super::pairhmm::{classify_by_pairhmm, classify_by_marginalized_pairhmm, ConfigurableGapParams};
use super::pangenome::build_haplotype_matrix;
use super::wfa_router::wfa_fast_path;
use super::utils::{find_read_pos, masked_dual_compare, masked_single_compare, median_qual, ClassifyResult, ClassifyPhase};
use super::AlignmentBackend;


/// Backend-aware Phase 3 classification.
///
/// Routes to Smith-Waterman (`classify_by_alignment`) or the pangenomic
/// WFA+PairHMM pipeline based on the active backend. Called from all
/// Phase 3 fallback sites in variant_checks (check_complex, check_insertion,
/// check_deletion).
///
/// For PairHMM backend:
/// 1. Build pangenomic haplotype matrix (H0/H1/H2..H2n) from variant + siblings
/// 2. Try WFA fast-path triage (edit distance; resolves ~70-80% of reads)
/// 3. Ambiguous reads fall through to marginalized PairHMM (BQ-aware)
/// 4. If haplotype matrix construction fails, falls back to SW
fn phase3_classify<F: Fn(u8, u8) -> i32>(
    record: &Record,
    variant: &Variant,
    siblings: &[Variant],
    quals: &[u8],
    min_baseq: u8,
    alt_aligner: &mut Aligner<F>,
    ref_aligner: &mut Aligner<F>,
    backend: &AlignmentBackend,
) -> ClassifyResult {
    match backend {
        AlignmentBackend::SmithWaterman => {
            // Full Phases 0-2.5 then Phase 3 SW
            check_complex(record, variant, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend)
        }
        AlignmentBackend::PairHMM {
            llr_threshold,
            gap_open,
            gap_extend,
            gap_open_repeat,
            gap_extend_repeat,
        } => {
            // ── PANGENOMIC WFA + MARGINALIZED PairHMM PIPELINE ──
            //
            // Step 1: Build multi-haplotype evaluation matrix.
            //   H0=REF, H1=ALT, H2..H2n=sibling germline combos
            if let Some(matrix) = build_haplotype_matrix(variant, siblings) {
                if let Some(ref _ctx) = variant.ref_context {
                    let win_start = variant.ref_context_start;
                    let win_end = win_start + variant.ref_context.as_ref().unwrap().len() as i64;

                    if let Some((sub_seq, sub_quals)) = extract_raw_read_window(
                        record, quals, win_start, win_end, variant.pos, variant.ref_allele.len()
                    ) {
                        if sub_seq.len() >= 3 {
                            // Step 2: WFA fast-path triage (edit distance).
                            // Clear-cut reads classified without PairHMM.
                            let med_qual = median_qual(&sub_quals, min_baseq);
                            if let Some(result) = wfa_fast_path(&sub_seq, &matrix, med_qual) {
                                trace!(
                                    "phase3: WFA resolved → is_ref={} is_alt={} qual={} at {}:{}",
                                    result.is_ref, result.is_alt, result.qual,
                                    variant.chrom, variant.pos + 1,
                                );
                                return result;
                            }

                            // Step 3: Ambiguous reads → marginalized PairHMM.
                            // Evaluates read against ALL haplotypes in the matrix,
                            // taking max LL over REF/ALT classes.
                            trace!(
                                "phase3: WFA ambiguous → escalating to marginalized PairHMM at {}:{}",
                                variant.chrom, variant.pos + 1,
                            );
                            let gap_params = ConfigurableGapParams::dynamic(
                                variant.repeat_span as usize,
                                *gap_open, *gap_extend,
                                *gap_open_repeat, *gap_extend_repeat,
                            );
                            return classify_by_marginalized_pairhmm(
                                &sub_seq, &sub_quals, &matrix,
                                min_baseq, &gap_params, *llr_threshold,
                            );
                        } else {
                            trace!(
                                "phase3: sub_seq too short ({} < 3) → SW fallback at {}:{}",
                                sub_seq.len(), variant.chrom, variant.pos + 1,
                            );
                        }
                    } else {
                        trace!(
                            "phase3: read window extraction failed → SW fallback at {}:{}",
                            variant.chrom, variant.pos + 1,
                        );
                    }
                } else {
                    debug!(
                        "phase3: variant {}:{} has no ref_context → SW fallback",
                        variant.chrom, variant.pos + 1,
                    );
                }
            } else {
                debug!(
                    "phase3: matrix construction failed for {}:{} → SW fallback",
                    variant.chrom, variant.pos + 1,
                );
            }
            // Fallback: if matrix construction or extraction fails, use SW
            trace!(
                "phase3: using check_complex (SW) fallback for {}:{}",
                variant.chrom, variant.pos + 1,
            );
            check_complex(record, variant, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend)
        }
    }
}


/// Result from MNP classification — distinguishes failure reason so the
/// caller can decide whether Phase 3 SW fallback is appropriate.
///
/// This replaces the previous `(bool, bool, u8)` return which conflated
/// quality failures, third alleles, and structural issues into a single
/// `(false, false, 0)` that always fell to Phase 3 SW (causing ~99%
/// MNP ALT loss via ties in haplotype-similar SW alignments).
pub enum MnpResult {
    /// All bases match REF, with median quality across the block.
    Ref(u8),
    /// All bases match ALT, with median quality across the block.
    Alt(u8),
    /// min(BQ) across the MNP block < threshold — skip read entirely.
    /// Matches C++ GBCMS behavior. Phase 3 fallback NOT appropriate:
    /// SW would count the read toward DP while C++ wouldn't.
    LowQuality,
    /// All bases pass quality, but the read sequence matches neither
    /// REF nor ALT. Phase 3 fallback NOT appropriate: if bases clearly
    /// don't match either allele, SW will likely generate a tie.
    ThirdAllele,
    /// Structural issue prevents string comparison:
    ///   - Read doesn't cover the entire MNP region
    ///   - Position not found in CIGAR walk
    ///   - Indel within the MNP block (contiguity check failed)
    ///
    /// Phase 3 fallback IS appropriate: the read may carry the variant
    /// through a complex alignment (e.g., complex variant annotated as MNP).
    Structural,
}


/// Returns `ClassifyResult` for SNP variants. Always Phase 0 (Structural).
/// Quality is the base quality at the variant position.
pub fn check_snp(record: &Record, variant: &Variant, quals: &[u8], min_baseq: u8) -> ClassifyResult {
    let read_pos = match find_read_pos(record, variant.pos) {
        Some(p) => p,
        None => return ClassifyResult::neither(ClassifyPhase::Structural),
    };

    let qual = quals[read_pos];
    if qual < min_baseq {
        return ClassifyResult::neither(ClassifyPhase::Structural);
    }

    let base = record.seq()[read_pos] as char;
    let ref_char = variant
        .ref_allele
        .chars()
        .next()
        .unwrap()
        .to_ascii_uppercase();
    let alt_char = variant
        .alt_allele
        .chars()
        .next()
        .unwrap()
        .to_ascii_uppercase();
    let base_upper = base.to_ascii_uppercase();

    let is_ref = base_upper == ref_char;
    let is_alt = base_upper == alt_char;
    // Return quality for fragment consensus scoring
    let base_qual = if is_ref || is_alt { qual } else { 0 };
    ClassifyResult::new(is_ref, is_alt, base_qual, ClassifyPhase::Structural)
}

/// Classify a read for an MNP variant using min-BQ-across-block strategy.
///
/// Uses the minimum base quality across the entire MNP block as the
/// quality gate (matching C++ GBCMS `baseCountDNP` behavior), rather
/// than per-base hard rejection which caused catastrophic ALT loss.
///
/// The contiguity check is performed FIRST (fail-fast for structural
/// issues) before quality and sequence comparison.
pub fn check_mnp(record: &Record, variant: &Variant, quals: &[u8], min_baseq: u8) -> MnpResult {
    let len = variant.ref_allele.len();

    // ── Step 1: Find read position of the first MNP base ──
    let start_read_pos = match find_read_pos(record, variant.pos) {
        Some(p) => p,
        None => return MnpResult::Structural,
    };

    // ── Step 2: Check full coverage ──
    if start_read_pos + len > record.seq().len() {
        return MnpResult::Structural;
    }

    // ── Step 3: Contiguity check FIRST (fail-fast for structural issues) ──
    // Verify no indels within the MNP block before doing quality/sequence.
    // This catches complex variants misannotated as MNPs.
    let end_read_pos = match find_read_pos(record, variant.pos + len as i64 - 1) {
        Some(p) => p,
        None => return MnpResult::Structural,
    };
    if end_read_pos - start_read_pos != len - 1 {
        return MnpResult::Structural; // Indel within MNP block
    }

    // ── Step 4: Min-BQ-across-block quality gate (C++ strategy) ──
    // NOTE: quals is passed from the caller — either raw record.qual() or BAQ-adjusted.
    let mut min_qual = u8::MAX;
    for i in 0..len {
        min_qual = min_qual.min(quals[start_read_pos + i]);
    }
    if min_qual < min_baseq {
        return MnpResult::LowQuality;
    }

    // ── Step 5: Direct string comparison (all bases passed quality) ──
    let seq = record.seq();
    let seq_bytes = seq.as_bytes();
    let ref_bytes = variant.ref_allele.as_bytes();
    let alt_bytes = variant.alt_allele.as_bytes();

    let mut matches_ref = true;
    let mut matches_alt = true;
    let mut mnp_quals: Vec<u8> = Vec::with_capacity(len);

    for i in 0..len {
        let pos = start_read_pos + i;
        mnp_quals.push(quals[pos]);

        let base = seq_bytes[pos].to_ascii_uppercase();
        if base != ref_bytes[i].to_ascii_uppercase() {
            matches_ref = false;
        }
        if base != alt_bytes[i].to_ascii_uppercase() {
            matches_alt = false;
        }
    }

    if matches_ref {
        MnpResult::Ref(median_qual(&mnp_quals, min_baseq))
    } else if matches_alt {
        MnpResult::Alt(median_qual(&mnp_quals, min_baseq))
    } else {
        MnpResult::ThirdAllele
    }
}


/// Check if a read supports a complex variant (indel + substitution).
///
/// Uses **haplotype reconstruction**: walks the CIGAR to rebuild what the read
/// shows for the genomic region covered by REF, then compares the reconstructed
/// sequence to both REF and ALT using **quality-aware masked comparison**.
///
/// ## Masked Comparison ("Reliable Intersection")
///
/// Instead of requiring exact byte-for-byte match, bases with quality below
/// `min_baseq` are **masked out** — they cannot vote for either allele. Only
/// "reliable" (high-quality) bases participate in the comparison.
///
/// Three cases based on reconstructed sequence length:
/// - **Case A** (`recon == alt == ref` length): simultaneous REF/ALT check with
///   ambiguity detection. If reliable bases match *both*, read is discarded.
/// - **Case B** (`recon == alt` length only): masked comparison against ALT only.
/// - **Case C** (`recon == ref` length only): masked comparison against REF only.
///
/// Returns `ClassifyResult` where base_qual is the median quality
/// across the reconstructed haplotype bases, used for fragment consensus.
pub fn check_complex<F: Fn(u8, u8) -> i32>(
    record: &Record,
    variant: &Variant,
    _siblings: &[Variant],
    quals: &[u8],
    min_baseq: u8,
    alt_aligner: &mut Aligner<F>,
    ref_aligner: &mut Aligner<F>,
    backend: &AlignmentBackend,
) -> ClassifyResult {
    let start_pos = variant.pos;
    let end_pos = variant.pos + variant.ref_allele.len() as i64; // exclusive

    let cigar = record.cigar();
    let mut ref_pos = record.pos();
    let mut read_pos: usize = 0;
    let seq = record.seq();
    // NOTE: quals is passed from the caller — either raw record.qual() or BAQ-adjusted.

    // Pre-allocate reconstruction buffers for performance
    let capacity = seq.len();
    let mut reconstructed_seq: Vec<u8> = Vec::with_capacity(capacity);
    let mut quals_per_base: Vec<u8> = Vec::with_capacity(capacity);

    trace!(
        "check_complex start: pos={} ref={} alt={}",
        start_pos, variant.ref_allele, variant.alt_allele
    );

    // --- Phase 0: Structural Anomaly Fast-Track ---
    // If the read has Soft-Clips (S) or explicit indels (I/D) within the window,
    // Phase 1's CIGAR-projected reconstruction will produce a severely truncated
    // or garbaged sequence. Phase 2 (masked comparison) then artificially matches
    // this truncated string perfectly to REF, erroneously rejecting complex ALT reads!
    // To prevent this false REF classification, we IMMEDIATELY route all
    // structurally anomalous reads (is_worth_realignment) to Phase 3 SW, which natively
    // extracts the raw unstructured bases and aligns against the full haplotype geometries.
    if let Some(ref ctx) = variant.ref_context {
        let win_start = variant.ref_context_start;
        let win_end = win_start + ctx.len() as i64;

        if is_worth_realignment(record, win_start, win_end) {
            if let Some((sub_seq, sub_quals)) = extract_raw_read_window(
                record, quals, win_start, win_end, variant.pos, variant.ref_allele.len()
            ) {
                if sub_seq.len() >= 3 {
                    trace!(
                        "check_complex: bypassing Phase 1/2 due to soft-clips/indels, Phase 3 extracted {} bases",
                        sub_seq.len()
                    );
                    return match backend {
                        AlignmentBackend::SmithWaterman => classify_by_alignment(
                            &sub_seq, &sub_quals, variant, min_baseq,
                            alt_aligner, ref_aligner,
                        ),
                        AlignmentBackend::PairHMM {
                            llr_threshold, gap_open, gap_extend,
                            gap_open_repeat, gap_extend_repeat,
                        } => {
                            let gap_params = ConfigurableGapParams::dynamic(
                                variant.repeat_span as usize,
                                *gap_open, *gap_extend,
                                *gap_open_repeat, *gap_extend_repeat,
                            );
                            classify_by_pairhmm(
                                &sub_seq, &sub_quals, variant, min_baseq,
                                &gap_params, *llr_threshold,
                            )
                        }
                    };
                }
            }
        }
    }

    // --- Phase 1: Haplotype Reconstruction ---
    // Walk the CIGAR to reconstruct what the read shows for [start_pos, end_pos).
    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len_i64 = *len as i64;
                let len_usize = *len as usize;

                // Intersection of [ref_pos, ref_pos + len) and [start_pos, end_pos)
                let overlap_start = std::cmp::max(ref_pos, start_pos);
                let overlap_end = std::cmp::min(ref_pos + len_i64, end_pos);

                if overlap_start < overlap_end {
                    let offset_in_op = (overlap_start - ref_pos) as usize;
                    let overlap_len = (overlap_end - overlap_start) as usize;
                    let current_read_pos = read_pos + offset_in_op;

                    for i in 0..overlap_len {
                        let p = current_read_pos + i;
                        if p >= seq.len() {
                            break;
                        }
                        reconstructed_seq.push(seq[p]);
                        quals_per_base.push(quals[p]);
                    }
                }
                ref_pos += len_i64;
                read_pos += len_usize;
            }
            Cigar::Ins(len) => {
                let len_usize = *len as usize;
                if ref_pos >= start_pos && ref_pos <= end_pos {
                    for i in 0..len_usize {
                        let p = read_pos + i;
                        if p >= seq.len() {
                            break;
                        }
                        reconstructed_seq.push(seq[p]);
                        quals_per_base.push(quals[p]);
                    }
                }
                read_pos += len_usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            }
            Cigar::SoftClip(len) => {
                let len_usize = *len as usize;
                // P1-2: Include soft-clipped bases that overlap the variant window.
                // Soft clips don't consume reference, so ref_pos is unchanged.
                // This recovers evidence from reads where the aligner clipped
                // the variant-supporting bases (inspired by VarDict's approach).
                if ref_pos >= start_pos && ref_pos < end_pos {
                    for i in 0..len_usize {
                        let p = read_pos + i;
                        if p >= seq.len() { break; }
                        reconstructed_seq.push(seq[p]);
                        quals_per_base.push(quals[p]);
                    }
                }
                read_pos += len_usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    let reconstructed_str = String::from_utf8_lossy(&reconstructed_seq);
    trace!(
        "Reconstructed: '{}' (len={})",
        reconstructed_str,
        reconstructed_seq.len()
    );

    // Median quality across the reconstructed haplotype for fragment consensus.
    // This is the quality we return for whichever allele matches.
    let med_haplotype_qual = median_qual(&quals_per_base, min_baseq);

    // --- Phase 2: Quality-Aware Masked Comparison ---
    // Mask out low-quality bases. Only reliable bases (qual >= min_baseq) vote.
    let alt_bytes = variant.alt_allele.as_bytes();
    let ref_bytes = variant.ref_allele.as_bytes();
    let recon_len = reconstructed_seq.len();
    let matches_alt_len = recon_len == alt_bytes.len();
    let matches_ref_len = recon_len == ref_bytes.len();

    // Guard: pathologically short reconstruction for large-REF variants.
    // When REF is much longer than ALT, a truncated reconstruction can
    // trivially match ALT length — reads that don't fully span the REF
    // region (soft-clipped, partial coverage) produce short reconstructions
    // that coincidentally match alt_len. This causes massive overcounting.
    //
    // Two tiers:
    // 1. Very large REF (>50bp or >1/3 read length): if recon < 10% of REF,
    //    skip Phase 2 entirely (original guard for 1kb+ deletions).
    // 2. REF significantly longer than ALT (>2x): if recon matches alt_len
    //    but not ref_len, the reconstruction is likely truncated, not ALT
    //    evidence. Skip Phase 2 for this case — let Phase 3 (SW/HMM) decide
    //    with full haplotype context.
    //
    // Example: ARID1A 1:27024008 REF=42bp ALT=5bp. Reads with 5bp recon
    // would match alt_len=5 in Case B → false ALT. Guard catches ref_len=42
    // > 2*5=10 and skips to Phase 3.
    let ref_len = ref_bytes.len();
    let alt_len = alt_bytes.len();
    let read_len = seq.len();
    let large_ref_threshold = std::cmp::max(50, read_len / 3);

    let skip_phase2 = if ref_len > large_ref_threshold && recon_len > 0 && recon_len < ref_len / 10 {
        // Tier 1: Massive REF (e.g., 1kb deletion), tiny recon
        trace!(
            "check_complex: recon_len={} is <10% of ref_len={} — \
             skipping Phase 2 (unreliable direct comparison)",
            recon_len, ref_len
        );
        true
    } else if ref_len > 2 * alt_len && matches_alt_len && !matches_ref_len {
        // Tier 2: REF >> ALT and recon matches only ALT length.
        // Reconstruction is likely truncated, not true ALT evidence.
        trace!(
            "check_complex: ref_len={} > 2*alt_len={}, recon_len={} matches alt_len \
             but not ref_len — skipping Phase 2 (likely truncated recon)",
            ref_len, alt_len, recon_len
        );
        true
    } else {
        false
    };

    if skip_phase2 {
        // Fall through to Phase 2.5 / Phase 3
    } else if matches_alt_len && matches_ref_len {
        // Case A: Equal-length REF and ALT — need simultaneous check + ambiguity detection
        let (mismatches_alt, mismatches_ref, reliable_count) =
            masked_dual_compare(&reconstructed_seq, &quals_per_base, alt_bytes, ref_bytes, min_baseq);

        trace!(
            "Case A: reliable={} mm_alt={} mm_ref={}",
            reliable_count, mismatches_alt, mismatches_ref
        );

        // Step 1: No reliable data → discard (MUST come first)
        if reliable_count == 0 {
            trace!("No reliable bases — discarding");
            return ClassifyResult::neither(ClassifyPhase::MaskedCompare);
        }

        // Step 2: Ambiguity — reliable bases match both alleles → discard
        if mismatches_alt == 0 && mismatches_ref == 0 {
            trace!("Ambiguous: reliable bases match both REF and ALT — discarding");
            return ClassifyResult::neither(ClassifyPhase::MaskedCompare);
        }

        // Step 3: Unambiguous match
        if mismatches_alt == 0 {
            trace!("Matches ALT on {} reliable bases, med_qual={}", reliable_count, med_haplotype_qual);
            return ClassifyResult::is_alt(med_haplotype_qual, ClassifyPhase::MaskedCompare);
        }
        if mismatches_ref == 0 {
            trace!("Matches REF on {} reliable bases, med_qual={}", reliable_count, med_haplotype_qual);
            return ClassifyResult::is_ref(med_haplotype_qual, ClassifyPhase::MaskedCompare);
        }

        // Step 4: Neither matches on reliable bases
        trace!("No match: mm_alt={} mm_ref={}", mismatches_alt, mismatches_ref);
    } else if matches_alt_len {
        // Case B: Only ALT length matches (e.g., DelIns) — no ambiguity possible
        let (mismatches, reliable_count) =
            masked_single_compare(&reconstructed_seq, &quals_per_base, alt_bytes, min_baseq);

        trace!(
            "Case B (ALT-only): reliable={} mismatches={}",
            reliable_count, mismatches
        );

        if reliable_count > 0 && mismatches == 0 {
            trace!("Matches ALT on {} reliable bases, med_qual={}", reliable_count, med_haplotype_qual);
            return ClassifyResult::is_alt(med_haplotype_qual, ClassifyPhase::MaskedCompare);
        }
    } else if matches_ref_len {
        // Case C: Only REF length matches — no ambiguity possible
        let (mismatches, reliable_count) =
            masked_single_compare(&reconstructed_seq, &quals_per_base, ref_bytes, min_baseq);

        trace!(
            "Case C (REF-only): reliable={} mismatches={}",
            reliable_count, mismatches
        );

        if reliable_count > 0 && mismatches == 0 {
            trace!("Matches REF on {} reliable bases, med_qual={}", reliable_count, med_haplotype_qual);
            return ClassifyResult::is_ref(med_haplotype_qual, ClassifyPhase::MaskedCompare);
        }
    } else {
        trace!(
            "Length mismatch: recon={} alt={} ref={}",
            recon_len,
            alt_bytes.len(),
            ref_bytes.len()
        );

        // Phase 2.5: Fuzzy edit distance fallback.
        // When reconstruction length doesn't match REF or ALT exactly
        // (e.g., incomplete MAF definition drops an adjacent SNV),
        // Levenshtein distance can still discriminate the closest allele.
        // Requires >1 edit margin for safety on very short strings.
        // Skip for large variants (>50bp) — O(n×m) is wasteful when Phase 3
        // SW handles them correctly with affine gap penalties.
        //
        // Also skip when REF >> ALT (>2x): Levenshtein is structurally
        // biased toward the shorter allele. A 20bp reconstruction has
        // d_alt ≈ 15 to a 5bp ALT, but d_ref ≈ 22-37 to a 42bp REF,
        // causing massive false ALT overcounting. Phase 3's full
        // haplotype alignment handles this correctly.
        if recon_len >= 2 && ref_bytes.len() <= 50 && alt_bytes.len() <= 50
            && ref_len <= 2 * alt_len
        {
            let d_ref = levenshtein(&reconstructed_seq, ref_bytes);
            let d_alt = levenshtein(&reconstructed_seq, alt_bytes);
            trace!(
                "Phase 2.5: edit_dist to_ref={} to_alt={} recon_len={}",
                d_ref, d_alt, recon_len
            );
            if d_alt + 1 < d_ref {
                trace!("Phase 2.5 → ALT (edit distance margin)");
                return ClassifyResult::is_alt(med_haplotype_qual, ClassifyPhase::Levenshtein);
            } else if d_ref + 1 < d_alt {
                trace!("Phase 2.5 → REF (edit distance margin)");
                return ClassifyResult::is_ref(med_haplotype_qual, ClassifyPhase::Levenshtein);
            }
            // else: ambiguous, fall through to Phase 3
        }
    }

    // --- Phase 3: Alignment-based fallback (indelpost approach) ---
    // When narrow-window reconstruction fails (FM1: D-truncation, FM2: adjacent I),
    // expand to the full ref_context window and use dual-haplotype SW alignment.
    //
    // CRITICAL: Use raw read window extraction (not CIGAR-projected) to preserve
    // the true biological sequence. For complex variants (e.g. EPHA7 REF=TCC ALT=CT),
    // BWA represents ALT reads as DEL+INS CIGARs. CIGAR-projected extraction
    // produces a hybrid sequence matching neither REF nor ALT haplotype.
    // Raw extraction gives the contiguous read bases that SW can correctly align.
    //
    // Pre-filter (indelpost pattern): only attempt SW for reads showing
    // evidence of carrying the variant (soft-clips, indels near window).
    // This eliminates ~80-90% of clean REF reads from expensive alignment.
    if let Some(ref ctx) = variant.ref_context {
        let win_start = variant.ref_context_start;
        let win_end = win_start + ctx.len() as i64;

        let is_mnp = variant.ref_allele.len() == variant.alt_allele.len() && variant.ref_allele.len() > 1;

        if !is_mnp && !is_worth_realignment(record, win_start, win_end) {
            trace!(
                "Phase 3 skipped: read has clean CIGAR over [{}, {})",
                win_start, win_end
            );
            return ClassifyResult::neither(ClassifyPhase::Alignment);
        }

        if let Some((sub_seq, sub_quals)) = extract_raw_read_window(
            record, quals, win_start, win_end, variant.pos, variant.ref_allele.len()
        ) {
            if sub_seq.len() >= 3 {
                trace!(
                    "Phase 3 fallback: extracted {} raw bases over [{}, {})",
                    sub_seq.len(), win_start, win_end
                );
                return match backend {
                    AlignmentBackend::SmithWaterman => classify_by_alignment(
                        &sub_seq, &sub_quals, variant, min_baseq,
                        alt_aligner, ref_aligner,
                    ),
                    AlignmentBackend::PairHMM {
                        llr_threshold, gap_open, gap_extend,
                        gap_open_repeat, gap_extend_repeat,
                    } => {
                        let gap_params = ConfigurableGapParams::dynamic(
                            variant.repeat_span as usize,
                            *gap_open, *gap_extend,
                            *gap_open_repeat, *gap_extend_repeat,
                        );
                        classify_by_pairhmm(
                            &sub_seq, &sub_quals, variant, min_baseq,
                            &gap_params, *llr_threshold,
                        )
                    }
                };
            }
        }
    }

    ClassifyResult::neither(ClassifyPhase::Alignment)
}


/// Check if a read supports an insertion variant.
///
/// Returns (is_ref, is_alt, base_qual) where base_qual is the quality of the
/// anchor base, used for fragment-level consensus scoring.
///
/// Uses a single CIGAR walk with three detection strategies:
/// 1. **Backward boundary check:** When anchor falls at the start of an M block
///    and the previous CIGAR op was a matching insertion (fixes off-by-one at
///    M/I/M boundaries where the aligner splits the match block).
/// 2. **Strict match (fast path):** Insertion immediately after the anchor base.
///    Returns ALT immediately if length + sequence match.
/// 3. **Windowed scan (fallback):** Any insertion within ±5bp of the anchor,
///    validated by three safeguards:
///    - S1: Inserted sequence matches expected ALT bases (quality-masked)
///    - S2: Closest match wins (minimum |shift_pos - anchor_pos|)
///    - S3: Reference base at shifted anchor matches original anchor base
///      (via variant.ref_context)
/// 4. **Phase 3 haplotype fallback:** When a length-matching insertion exists
///    nearby but fails the sequence check (e.g., same biological event
///    represented differently by caller vs aligner), falls back to
///    check_complex for Smith-Waterman haplotype comparison.
pub fn check_insertion<F: Fn(u8, u8) -> i32>(
    record: &Record,
    variant: &Variant,
    siblings: &[Variant],
    quals: &[u8],
    min_baseq: u8,
    alt_aligner: &mut Aligner<F>,
    ref_aligner: &mut Aligner<F>,
    backend: &AlignmentBackend,
) -> ClassifyResult {
    let cigar_view = record.cigar();
    // NOTE: quals is passed from the caller — either raw record.qual() or BAQ-adjusted.
    let mut ref_pos = record.pos();
    let mut read_pos: usize = 0;

    let anchor_pos = variant.pos;
    let expected_ins_len = variant.alt_allele.len() - 1; // VCF ALT includes anchor
    let expected_ins_seq = &variant.alt_allele.as_bytes()[1..]; // ALT without anchor
    let original_anchor_base = variant.ref_allele.as_bytes()[0].to_ascii_uppercase();

    // Windowed scan parameters — scales with repeat_span for MSI regions
    let window: i64 = std::cmp::max(5, variant.repeat_span as i64 + 2);
    let window_start = (anchor_pos - window).max(0);
    let window_end = anchor_pos + window;

    // State tracked across the CIGAR walk
    let mut found_ref_coverage = false;
    let mut anchor_read_pos: Option<usize> = None; // read position of anchor base
    let mut best_windowed_match: Option<u64> = None; // distance of best windowed match
    let mut has_nearby_length_match = false; // length-matching Ins found but seq check failed

    for (i, op) in cigar_view.iter().enumerate() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len_i64 = *len as i64;
                let block_end = ref_pos + len_i64;

                // Track anchor read position whenever we encounter it
                if anchor_pos >= ref_pos && anchor_pos < block_end {
                    let offset = (anchor_pos - ref_pos) as usize;
                    anchor_read_pos = Some(read_pos + offset);
                }

                // --- P0-2: Backward boundary check ---
                // When anchor falls at block_end of prior M block, CIGAR geometry
                // places it at ref_pos of THIS block (after the Ins was consumed).
                // Check backward: was the previous op a matching insertion?
                if anchor_pos == ref_pos && i > 0 {
                    if let Some(Cigar::Ins(ins_len)) = cigar_view.get(i - 1) {
                        let ins_len_usize = *ins_len as usize;
                        if ins_len_usize == expected_ins_len {
                            let ins_read_start = read_pos - ins_len_usize;
                            if ins_read_start + ins_len_usize <= record.seq().len() {
                                let ins_seq = &record.seq().as_bytes()
                                    [ins_read_start..ins_read_start + ins_len_usize];
                                // P1-1: Quality-aware fuzzy match
                                let ins_quals = &quals[ins_read_start..ins_read_start + ins_len_usize];
                                let (mismatches, reliable) = masked_single_compare(
                                    ins_seq, ins_quals, expected_ins_seq, min_baseq
                                );
                                if reliable > 0 && mismatches == 0 {
                                    let qual = if read_pos < quals.len() { quals[read_pos] } else { 0 };
                                    trace!(
                                        "check_insertion: backward boundary match at pos {}, qual={}",
                                        anchor_pos, qual
                                    );
                                    return ClassifyResult::is_alt(qual, ClassifyPhase::Structural); // ALT — backward match
                                }
                            }
                        }
                    }
                }

                // --- Strict fast path: anchor at end of this block ---
                if anchor_pos >= ref_pos && anchor_pos < block_end {
                    if anchor_pos == block_end - 1 {
                        // Anchor is the last base of this match block.
                        // Check if next op is an insertion with matching length + sequence.
                        if let Some(Cigar::Ins(ins_len)) = cigar_view.get(i + 1) {
                            let ins_len_usize = *ins_len as usize;
                            if ins_len_usize == expected_ins_len {
                                let ins_start = read_pos + *len as usize;
                                if ins_start + ins_len_usize <= record.seq().len() {
                                    let ins_seq = &record.seq().as_bytes()
                                        [ins_start..ins_start + ins_len_usize];
                                    // P1-1: Quality-aware fuzzy match
                                    let ins_quals = &quals[ins_start..ins_start + ins_len_usize];
                                    let (mismatches, reliable) = masked_single_compare(
                                        ins_seq, ins_quals, expected_ins_seq, min_baseq
                                    );
                                    if reliable > 0 && mismatches == 0 {
                                        let arp = anchor_read_pos.unwrap_or(0);
                                        let qual = if arp < quals.len() { quals[arp] } else { 0 };
                                        trace!(
                                            "check_insertion: strict match at pos {}, anchor_qual={}",
                                            anchor_pos, qual
                                        );
                                        return ClassifyResult::is_alt(qual, ClassifyPhase::Structural); // ALT — strict match
                                    }
                                }
                            }
                        }
                        // Anchor at end but no matching insertion → REF coverage
                        found_ref_coverage = true;
                    } else {
                        // Anchor in middle of match block → read covers anchor without insertion
                        found_ref_coverage = true;
                    }
                }

                // --- Windowed scan: check if any Ins after this block is within window ---
                if let Some(Cigar::Ins(ins_len)) = cigar_view.get(i + 1) {
                    let ins_ref_pos = block_end; // genomic position where insertion occurs
                    if ins_ref_pos >= window_start && ins_ref_pos <= window_end
                        && ins_ref_pos != anchor_pos + 1 // skip strict position (already handled)
                    {
                        let ins_len_usize = *ins_len as usize;
                        // Safeguard 1: length must match
                        if ins_len_usize == expected_ins_len {
                            let ins_start = read_pos + *len as usize;
                            if ins_start + ins_len_usize <= record.seq().len() {
                                let ins_seq = &record.seq().as_bytes()
                                    [ins_start..ins_start + ins_len_usize];
                                // Quality-aware fuzzy match for inserted bases
                                let ins_quals = &quals[ins_start..ins_start + ins_len_usize];
                                let (mismatches, reliable) = masked_single_compare(
                                    ins_seq, ins_quals, expected_ins_seq, min_baseq
                                );
                                if reliable > 0 && mismatches == 0 {
                                    // Safeguard 3: verify anchor base at shifted position
                                    let shifted_anchor_pos = ins_ref_pos - 1;
                                    let anchor_ok = match &variant.ref_context {
                                        Some(ctx) => {
                                            let ctx_offset = (shifted_anchor_pos
                                                - variant.ref_context_start)
                                                as usize;
                                            if ctx_offset < ctx.len() {
                                                ctx.as_bytes()[ctx_offset].to_ascii_uppercase()
                                                    == original_anchor_base
                                            } else {
                                                trace!(
                                                    "ref_context offset {} out of bounds (len={}), rejecting",
                                                    ctx_offset, ctx.len()
                                                );
                                                false
                                            }
                                        }
                                        None => {
                                            warn!("ref_context is None for variant at {}:{} — S3 cannot validate shifted insertion",
                                                  variant.chrom, variant.pos + 1);
                                            false
                                        },
                                    };

                                    if anchor_ok {
                                        // Safeguard 2: track closest match
                                        let distance =
                                            (ins_ref_pos - (anchor_pos + 1)).unsigned_abs();
                                        if best_windowed_match
                                            .is_none_or(|prev| distance < prev)
                                        {
                                            best_windowed_match = Some(distance);
                                        }
                                    } else {
                                        trace!(
                                            "check_insertion: S3 reject at shifted pos {} \
                                             (anchor base mismatch)",
                                            shifted_anchor_pos
                                        );
                                    }
                                } else {
                                    // Length matches but sequence differs — the caller
                                    // and aligner may represent the same event
                                    // differently (e.g., shifted insertion in a repeat).
                                    // Track this so Phase 3 SW can arbitrate.
                                    has_nearby_length_match = true;
                                    trace!(
                                        "check_insertion: windowed I({}) at pos {} seq \
                                         mismatch (mismatches={}, reliable={}), \
                                         flagging for Phase 3 fallback",
                                        ins_len_usize, ins_ref_pos,
                                        mismatches, reliable
                                    );
                                }
                            }
                        }
                    }
                }

                ref_pos = block_end;
                read_pos += *len as usize;
            }
            Cigar::Ins(len) => {
                read_pos += *len as usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            }
            Cigar::SoftClip(len) => {
                read_pos += *len as usize;
            }
            _ => {}
        }
    }

    // Anchor quality for fragment consensus (used for both ALT and REF returns)
    let anchor_qual = anchor_read_pos
        .filter(|&p| p < quals.len())
        .map(|p| quals[p])
        .unwrap_or(0);

    // Evaluate results after full CIGAR walk
    if best_windowed_match.is_some() {
        trace!(
            "check_insertion: windowed match for variant at pos {}, anchor_qual={}",
            anchor_pos, anchor_qual
        );
        return ClassifyResult::is_alt(anchor_qual, ClassifyPhase::CigarRecon); // ALT — windowed match
    }

    // Phase 3 haplotype fallback: when a length-matching insertion exists nearby
    // but the sequence check failed, the caller and aligner may represent the
    // same biological event differently (e.g., FLT4 A→GGAT placed at different
    // positions with different inserted bases). Suppress found_ref_coverage to
    // let check_complex → Smith-Waterman arbitrate using full haplotype comparison.
    if has_nearby_length_match && found_ref_coverage {
        trace!(
            "check_insertion: nearby I({}) with seq mismatch at pos {}, \
             falling back to check_complex for Phase 3 SW",
            expected_ins_len, anchor_pos
        );
        return phase3_classify(record, variant, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend);
    }

    if found_ref_coverage {
        // CIGAR is definitive for pure insertions: if the read's Match op
        // covers the anchor and no matching I op was found (strict or
        // windowed), the read is REF. Soft-clipping elsewhere cannot
        // represent a missing insertion at the anchor position.
        //
        // Note: reads with soft-clip AT the anchor have
        // found_ref_coverage=false and are correctly routed to Phase 3
        // via the !found_ref_coverage path above.
        return ClassifyResult::is_ref(anchor_qual, ClassifyPhase::Structural);
    }
    ClassifyResult::neither(ClassifyPhase::Structural) // Read does not cover the variant region
}

/// Check if a read supports a deletion variant.
///
/// Returns (is_ref, is_alt, base_qual) where base_qual is the quality of the
/// anchor base, used for fragment-level consensus scoring.
///
/// Uses the same single-walk strategy as check_insertion:
/// 1. **Strict match (fast path):** Deletion immediately after anchor, length matches.
/// 2. **Windowed scan (fallback):** Any deletion within ±5bp, validated by:
///    - S1: Deletion length matches expected
///    - S2: Closest match wins
///    - S3: Reference bases at shifted position match expected deleted sequence
///      (via variant.ref_context)
/// 3. **Haplotype fallback:** When CIGAR geometry doesn't match (e.g. different
///    breakpoint placement or wrong deletion length), delegates to `check_complex`
///    for quality-aware haplotype comparison.
pub fn check_deletion<F: Fn(u8, u8) -> i32>(
    record: &Record,
    variant: &Variant,
    siblings: &[Variant],
    quals: &[u8],
    min_baseq: u8,
    alt_aligner: &mut Aligner<F>,
    ref_aligner: &mut Aligner<F>,
    backend: &AlignmentBackend,
) -> ClassifyResult {
    let cigar_view = record.cigar();
    // NOTE: quals is passed from the caller — either raw record.qual() or BAQ-adjusted.
    let mut ref_pos = record.pos();
    let mut read_pos: usize = 0;

    let anchor_pos = variant.pos;
    let expected_del_len = variant.ref_allele.len() - 1; // REF without anchor
    // The expected deleted bases (REF without the anchor base)
    let expected_del_seq = &variant.ref_allele.as_bytes()[1..];

    // Windowed scan parameters — scales with repeat_span for MSI regions
    let window: i64 = std::cmp::max(5, variant.repeat_span as i64 + 2);
    let window_start = (anchor_pos - window).max(0);
    let window_end = anchor_pos + window;

    let mut found_ref_coverage = false;
    let mut anchor_read_pos: Option<usize> = None; // read position of anchor base
    let mut best_windowed_match: Option<u64> = None;


    for (i, op) in cigar_view.iter().enumerate() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len_i64 = *len as i64;
                let block_end = ref_pos + len_i64;

                // Track anchor read position whenever we encounter it
                if anchor_pos >= ref_pos && anchor_pos < block_end {
                    let offset = (anchor_pos - ref_pos) as usize;
                    anchor_read_pos = Some(read_pos + offset);
                }

                // --- Strict fast path ---
                if anchor_pos >= ref_pos && anchor_pos < block_end {
                    if anchor_pos == block_end - 1 {
                        // Anchor at end of match. Check if next op is a deletion.
                        if let Some(Cigar::Del(del_len)) = cigar_view.get(i + 1) {
                            if *del_len as usize == expected_del_len {
                                let arp = anchor_read_pos.unwrap_or(0);
                                let qual = if arp < quals.len() { quals[arp] } else { 0 };
                                trace!(
                                    "check_deletion: strict match at pos {}, anchor_qual={}",
                                    anchor_pos, qual
                                );
                                return ClassifyResult::is_alt(qual, ClassifyPhase::Structural); // ALT — strict match
                            } else {
                                // P0-3: D found at anchor but wrong length.
                                // Use SV-caller-style reciprocal overlap matching:
                                // aligners often report slightly different breakpoints
                                // for the same large deletion, producing different
                                // CIGAR D lengths. If both start at the same anchor
                                // and share ≥50% reciprocal overlap, treat as the
                                // same biological event.
                                // Precedent: SURVIVOR uses ≥50% overlap, BEDTools
                                // uses configurable reciprocal overlap for SV matching.
                                let found_del_len = *del_len as usize;
                                let min_del = expected_del_len.min(found_del_len);
                                let max_del = expected_del_len.max(found_del_len);
                                let reciprocal_overlap =
                                    min_del as f64 / max_del as f64;

                                if expected_del_len >= 50 && reciprocal_overlap >= 0.5 {
                                    // Large deletion with significant overlap → ALT
                                    let arp = anchor_read_pos.unwrap_or(0);
                                    let qual = if arp < quals.len() {
                                        quals[arp]
                                    } else {
                                        0
                                    };
                                    trace!(
                                        "check_deletion: tolerant match D({}) ≈ D({}) \
                                         at pos {} (reciprocal_overlap={:.2}, \
                                         threshold=0.50, anchor_qual={})",
                                        found_del_len,
                                        expected_del_len,
                                        anchor_pos,
                                        reciprocal_overlap,
                                        qual
                                    );
                                    return ClassifyResult::is_alt(qual, ClassifyPhase::Structural); // ALT — tolerant match
                                }

                                // Small deletion or low overlap: fall back to
                                // check_complex for haplotype-based comparison.
                                trace!(
                                    "check_deletion: D({}) at anchor {} but expected D({}), \
                                     reciprocal_overlap={:.2} (below 0.50 or del<50bp), \
                                     falling back to check_complex",
                                    found_del_len,
                                    anchor_pos,
                                    expected_del_len,
                                    reciprocal_overlap
                                );
                                return phase3_classify(record, variant, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend);
                            }
                        }
                        found_ref_coverage = true;
                    } else {
                        // Anchor in middle of match → REF coverage
                        found_ref_coverage = true;
                    }
                }

                // --- Windowed scan: check for Del after this block within window ---
                if let Some(Cigar::Del(del_len)) = cigar_view.get(i + 1) {
                    let del_ref_pos = block_end; // genomic position where deletion starts
                    if del_ref_pos >= window_start && del_ref_pos <= window_end
                        && del_ref_pos != anchor_pos + 1 // skip strict position
                    {
                        let del_len_usize = *del_len as usize;

                        // Safeguard 1: deletion length check.
                        // Accept exact matches, OR for large deletions (≥50bp),
                        // accept reciprocal overlap ≥50% (same logic as Fix 1
                        // on the strict path). Aligners often report slightly
                        // different breakpoints for the same biological event.
                        let length_ok = if del_len_usize == expected_del_len {
                            true
                        } else if expected_del_len >= 50 {
                            let min_del = expected_del_len.min(del_len_usize);
                            let max_del = expected_del_len.max(del_len_usize);
                            let overlap = min_del as f64 / max_del as f64;
                            if overlap >= 0.5 {
                                trace!(
                                    "check_deletion: windowed tolerant match \
                                     D({}) ≈ D({}) at pos {} (overlap={:.2})",
                                    del_len_usize, expected_del_len,
                                    del_ref_pos, overlap
                                );

                                true
                            } else {
                                false
                            }
                        } else {
                            false
                        };

                        if length_ok {
                            // Safeguard 3: verify the deleted reference bases match
                            // (only for exact-length matches; skip for tolerant
                            // matches where lengths differ)
                            let del_ok = if del_len_usize != expected_del_len {
                                true // tolerant match — length differs, skip seq check
                            } else {
                                match &variant.ref_context {
                                    Some(ctx) => {
                                        let ctx_bytes = ctx.as_bytes();
                                        let ctx_offset_i64 =
                                            del_ref_pos - variant.ref_context_start;
                                        if ctx_offset_i64 < 0 {
                                            trace!(
                                                "ref_context offset negative ({}), rejecting shifted deletion",
                                                ctx_offset_i64
                                            );
                                            false
                                        } else {
                                            let ctx_offset = ctx_offset_i64 as usize;
                                            if ctx_offset + del_len_usize <= ctx_bytes.len() {
                                                let ref_at_shift = &ctx_bytes
                                                    [ctx_offset..ctx_offset + del_len_usize];
                                                let ref_quals = vec![u8::MAX; del_len_usize];
                                                let (mismatches, reliable) = masked_single_compare(
                                                    ref_at_shift, &ref_quals, expected_del_seq, 0
                                                );
                                                reliable > 0 && mismatches == 0
                                            } else {
                                                trace!(
                                                    "ref_context offset {} out of bounds (len={}), rejecting",
                                                    ctx_offset, ctx_bytes.len()
                                                );
                                                false
                                            }
                                        }
                                    }
                                    None => {
                                        warn!("ref_context is None for variant at {}:{} — S3 cannot validate shifted deletion",
                                              variant.chrom, variant.pos + 1);
                                        false
                                    },
                                }
                            };

                            if del_ok {
                                // Safeguard 2: track closest match
                                let distance =
                                    (del_ref_pos - (anchor_pos + 1)).unsigned_abs();
                                if best_windowed_match
                                    .is_none_or(|prev| distance < prev)
                                {
                                    best_windowed_match = Some(distance);
                                }
                            } else {
                                trace!(
                                    "check_deletion: S3 reject at shifted pos {} \
                                     (deleted bases mismatch)",
                                    del_ref_pos
                                );
                            }
                        }
                    }
                }

                ref_pos = block_end;
                read_pos += *len as usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            }
            Cigar::Ins(len) => {
                read_pos += *len as usize;
            }
            Cigar::SoftClip(len) => {
                read_pos += *len as usize;
            }
            _ => {}
        }
    }

    // Anchor quality for fragment consensus (used for both ALT and REF returns)
    let anchor_qual = anchor_read_pos
        .filter(|&p| p < quals.len())
        .map(|p| quals[p])
        .unwrap_or(0);

    // Evaluate results after full CIGAR walk
    if best_windowed_match.is_some() {
        trace!(
            "check_deletion: windowed match for variant at pos {}, anchor_qual={}",
            anchor_pos, anchor_qual
        );
        return ClassifyResult::is_alt(anchor_qual, ClassifyPhase::CigarRecon); // ALT — windowed match
    }

    // Note: an earlier version had an "interior REF guard" here that classified
    // any read starting inside a large deletion span as REF. This was REMOVED
    // because it massively overcounted: for a 1023bp deletion (TP53), it claimed
    // ~4,000 reads mapping anywhere in the 1kb interior as REF evidence, when
    // IGV shows only ~770 reads at the anchor. Reads that don't cover the anchor
    // have no information about whether the deletion is present and must not be
    // counted. With the SW swap fix (Issue #1), Phase 3 correctly handles large
    // deletions: the read (pattern) slides along the longer haplotype (text)
    // without gap penalty, so false ALT calls no longer occur.

    // P0-3: Haplotype fallback — when strict/windowed CIGAR matching found no
    // deletion match and the read doesn't cover the anchor, try check_complex
    // which reconstructs the read's haplotype and does quality-aware comparison.
    // Only fall back when NOT found_ref_coverage to avoid false positives on
    // reads that genuinely show REF at this position.
    //
    // CRITICAL: Only attempt Phase 3 if the read actually overlaps the anchor
    // position (variant.pos). Reads that map entirely inside a large deletion
    // span (e.g., a 1023bp TP53 deletion) have no information about the variant
    // and must not be counted. Without this guard, the SW aligner would classify
    // interior reads as REF, massively inflating ref_count.
    if !found_ref_coverage && best_windowed_match.is_none() {
        // Compute the read's reference span end from CIGAR
        let read_ref_end = {
            let mut rend = record.pos();
            for op in record.cigar().iter() {
                match op {
                    Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len)
                    | Cigar::Del(len) | Cigar::RefSkip(len) => {
                        rend += *len as i64;
                    }
                    _ => {}
                }
            }
            rend
        };

        // Read must span the anchor position to have variant information.
        // anchor_pos is the 0-based position of the base before the deletion.
        if record.pos() <= anchor_pos && read_ref_end > anchor_pos {
            trace!(
                "check_deletion: no CIGAR match at pos {}, falling back to check_complex",
                anchor_pos
            );
            return phase3_classify(record, variant, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend);
        }
        // Otherwise: read doesn't overlap the anchor → no variant info
    }

    if found_ref_coverage {
        // CIGAR is definitive for pure deletions: if the read's Match op
        // covers the anchor and no matching D op was found (strict or
        // windowed), the read is REF. Soft-clipping elsewhere cannot
        // represent a missing deletion at the anchor position.
        //
        // Note: reads with soft-clip AT the anchor have
        // found_ref_coverage=false and are correctly routed to Phase 3
        // via the !found_ref_coverage path above.
        return ClassifyResult::is_ref(anchor_qual, ClassifyPhase::Structural);
    }
    ClassifyResult::neither(ClassifyPhase::Structural) // Read does not cover the variant region
}
