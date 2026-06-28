//! Per-variant-type classification functions.
//!
//! Contains the variant-type-specific allele checkers:
//! - `check_snp` — single nucleotide polymorphism
//! - `check_mnp` — multi-nucleotide polymorphism (with `MnpResult`)
//! - `check_insertion` — insertion with windowed CIGAR scan
//! - `check_deletion` — deletion with windowed CIGAR scan
//! - `check_complex` — complex variant (indel + substitution) with
//!   phased Masked Comparison → Levenshtein → SW alignment pipeline
//!
//! ## N-base handling
//!
//! All checkers detect N bases (from duplex masking or sequencer failure)
//! and signal `has_n_base=true` via `ClassifyResult`:
//! - **SNP**: N base → `ClassifyResult::neither_n()` (uninformative)
//! - **MNP**: N at discriminating position → masked (doesn't vote), but
//!   `had_n_base` flag propagated through `MnpResult`
//! - **Complex**: N in reconstructed haplotype → detected via scan,
//!   propagated through all 10 return paths
//!
//! The engine uses `has_n_base` to increment `n_count` for duplex masking QC.

use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;
use bio::alignment::distance::levenshtein;
use bio::alignment::pairwise::Aligner;
use log::{debug, trace, warn};

use crate::types::Variant;
use super::alignment::{classify_by_alignment, extract_raw_read_window, is_worth_realignment};
use super::pairhmm::{classify_by_marginalized_pairhmm, ConfigurableGapParams};
use super::pangenome::build_haplotype_matrix;
use super::wfa_router::wfa_fast_path;
use super::utils::{find_read_pos, masked_dual_compare, masked_single_compare, median_qual, ClassifyResult, ClassifyPhase};
use super::AlignmentBackend;


/// Pangenomic WFA + marginalized PairHMM classification pipeline.
///
/// Standalone helper that takes pre-extracted read data and attempts
/// classification against the full sibling-aware haplotype matrix.
///
/// ## Pipeline
///
/// 1. Build pangenomic haplotype matrix (H0/H1/H2..H2n) from variant + siblings
/// 2. WFA edit-distance triage — resolves ~70-80% of reads instantly
/// 3. Ambiguous reads → marginalized PairHMM (BQ-aware probabilistic)
///
/// ## Returns
///
/// - `Some(ClassifyResult)` if classification succeeded
/// - `None` if the pipeline cannot produce a result (matrix build failure,
///   no ref_context, etc.) — caller should fall back to SW
///
/// ## Usage
///
/// Called from:
/// - `phase3_classify()` — for insertion/deletion fallbacks
/// - `check_complex()` — Phase 0 (structural bypass) and Phase 3 (fallback)
#[allow(clippy::too_many_arguments)]
fn pangenomic_classify(
    sub_seq: &[u8],
    sub_quals: &[u8],
    variant: &Variant,
    siblings: &[Variant],
    min_baseq: u8,
    gap_open: f64,
    gap_extend: f64,
    gap_open_repeat: f64,
    gap_extend_repeat: f64,
    llr_threshold: f64,
) -> Option<ClassifyResult> {
    // Step 1: Build multi-haplotype evaluation matrix.
    //   H0=REF, H1=ALT, H2..H2n=sibling germline combos
    let matrix = match build_haplotype_matrix(variant, siblings) {
        Some(m) => m,
        None => {
            debug!(
                "pangenomic_classify: matrix construction failed for {}:{} (siblings={})",
                variant.chrom, variant.pos + 1, siblings.len(),
            );
            return None;
        }
    };

    // Step 2: WFA fast-path triage (edit distance).
    // Clear-cut reads classified without PairHMM.
    let med_qual = median_qual(sub_quals, min_baseq);
    if let Some(result) = wfa_fast_path(sub_seq, sub_quals, min_baseq, &matrix, med_qual) {
        trace!(
            "pangenomic_classify: WFA resolved → is_ref={} is_alt={} qual={} at {}:{}",
            result.is_ref, result.is_alt, result.qual,
            variant.chrom, variant.pos + 1,
        );
        return Some(result);
    }

    // Step 3: Ambiguous reads → marginalized PairHMM.
    // Evaluates read against ALL haplotypes in the matrix,
    // taking max LL over REF/ALT classes.
    trace!(
        "pangenomic_classify: WFA ambiguous → escalating to marginalized PairHMM at {}:{} (haplotypes={})",
        variant.chrom, variant.pos + 1, matrix.len(),
    );
    let gap_params = ConfigurableGapParams::dynamic(
        variant.repeat_span,
        gap_open, gap_extend,
        gap_open_repeat, gap_extend_repeat,
    );
    Some(classify_by_marginalized_pairhmm(
        sub_seq, sub_quals, &matrix,
        min_baseq, &gap_params, llr_threshold,
    ))
}


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
#[allow(clippy::too_many_arguments)]
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
            // Try pangenomic WFA+PairHMM pipeline via pangenomic_classify().
            // Handles: matrix build → WFA triage → marginalized PairHMM.
            if let Some(ref _ctx) = variant.ref_context {
                let win_start = variant.ref_context_start;
                let win_end = win_start + variant.ref_context.as_ref().unwrap().len() as i64;

                if let Some((sub_seq, sub_quals)) = extract_raw_read_window(
                    record, quals, win_start, win_end, variant.pos, variant.ref_allele.len()
                ) {
                    if sub_seq.len() >= 3 {
                        if let Some(result) = pangenomic_classify(
                            &sub_seq, &sub_quals, variant, siblings,
                            min_baseq, *gap_open, *gap_extend,
                            *gap_open_repeat, *gap_extend_repeat, *llr_threshold,
                        ) {
                            return result;
                        }
                        // pangenomic_classify returned None (matrix build failed) —
                        // fall through to SW fallback below.
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
            // Fallback: if pangenomic pipeline fails at any stage, use SW via check_complex
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
#[derive(Debug)]
pub enum MnpResult {
    /// All unmasked bases match REF. (quality, had_n_at_any_position)
    Ref(u8, bool),
    /// All unmasked bases match ALT. (quality, had_n_at_any_position)
    Alt(u8, bool),
    /// All discriminating positions masked (BQ < threshold or N).
    /// Carries (positions_matching_alt, had_n_at_any_position).
    /// Used for `partial_alt` counting when positions_matching_alt > 0.
    /// Matches C++ GBCMS behavior: read not counted for AD.
    /// Phase 3 fallback NOT appropriate: SW would count toward DP while C++ wouldn't.
    LowQuality(u8, bool),
    /// Unmasked bases match neither REF nor ALT (mixed or third-allele).
    /// Carries (positions_matching_alt, had_n_at_any_position).
    /// Phase 3 fallback NOT appropriate: if bases clearly
    /// don't match either allele, SW will likely generate a tie.
    ThirdAllele(u8, bool),
    /// Structural issue prevents string comparison:
    ///   - Read doesn't cover the entire MNP region
    ///   - Position not found in CIGAR walk
    ///   - Indel within the MNP block (contiguity check failed)
    ///
    /// Phase 3 fallback IS appropriate: the read may carry the variant
    /// through a complex alignment (e.g., complex variant annotated as MNP).
    /// No had_n field: structural issues route to check_complex which
    /// independently detects N in the reconstructed haplotype and propagates
    /// has_n_base through its own Phase 2 return paths.
    Structural,
}


/// Returns `ClassifyResult` for SNP variants. Always Phase 0 (Structural).
///
/// Quality is the base quality at the variant position.
/// N bases are treated as uninformative (neither REF nor ALT), matching GATK's
/// approach: the base quality gate (BQ < `min_baseq`) catches most N bases from
/// duplex collapsing (fgbio assigns BQ ≈ 2 to masked positions), but this
/// explicit guard handles raw BAMs where N may have arbitrary BQ.
pub fn check_snp(record: &Record, variant: &Variant, quals: &[u8], min_baseq: u8) -> ClassifyResult {
    let read_pos = match find_read_pos(record, variant.pos) {
        Some(p) => p,
        None => return ClassifyResult::neither(ClassifyPhase::Structural),
    };

    let qual = quals[read_pos];
    if qual < min_baseq {
        return ClassifyResult::neither(ClassifyPhase::Structural);
    }

    let base = record.seq()[read_pos];

    // N base = uninformative (duplex masking, sequencer failure).
    // Treat as masked regardless of BQ — do not count as REF, ALT, or third allele.
    // Defense-in-depth: fgbio assigns BQ ≈ 2 to N bases, but raw BAMs may not.
    // Use neither_n() to signal has_n_base=true → engine increments n_count.
    if base == b'N' || base == b'n' {
        trace!("SNP N guard: base={} at read_pos={}, returning neither_n", base as char, read_pos);
        return ClassifyResult::neither_n(ClassifyPhase::Structural);
    }

    let ref_char = variant.ref_allele.as_bytes()[0].to_ascii_uppercase();
    let alt_char = variant.alt_allele.as_bytes()[0].to_ascii_uppercase();
    let base_upper = base.to_ascii_uppercase();

    let is_ref = base_upper == ref_char;
    let is_alt = base_upper == alt_char;
    // Return quality for fragment consensus scoring
    let base_qual = if is_ref || is_alt { qual } else { 0 };
    ClassifyResult::new(is_ref, is_alt, base_qual, ClassifyPhase::Structural)
}

/// Classify a read for an MNP variant using selective quality gating.
///
/// Quality is gated only at **discriminating positions** (where REF ≠ ALT).
/// Uninformative positions (where REF == ALT) don't affect allele
/// classification and should not cause valid reads to be dropped.
///
/// This is an improvement over C++ GBCMS `baseCountDNP`, which gates on
/// min(BQ) across ALL positions. For GC-rich regions (e.g., TERT promoter)
/// or long ONPs with few discriminating positions, the old strategy dropped
/// reads due to low quality at positions that don't matter for classification.
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

    // ── Step 4+5: Masked per-position evaluation ──
    // Each discriminating position (REF ≠ ALT) is independently evaluated:
    //   - Masked: BQ < min_baseq OR N base → cannot vote (uninformative)
    //   - Unmasked: votes REF if matches REF, ALT if matches ALT, other if neither
    //
    // This replaces the old aggregate min-BQ gate that dropped the entire read
    // when ANY discriminating position had low quality. The masked approach
    // recovers reads where only some positions are low-quality (e.g., duplex N
    // at one MNP position), matching the GATK approach.
    //
    // N bases at discriminating positions are treated as BQ=0 regardless of
    // their reported quality (defense-in-depth for raw BAMs; duplex BAMs
    // already assign BQ ≈ 2 to N bases via fgbio).
    let ref_bytes = variant.ref_allele.as_bytes();
    let alt_bytes = variant.alt_allele.as_bytes();
    let seq = record.seq();
    let seq_bytes = seq.as_bytes();

    let mut n_discriminating: usize = 0;
    let mut n_masked: usize = 0;
    let mut n_unmasked_match_alt: usize = 0;
    let mut n_unmasked_match_ref: usize = 0;
    let mut n_unmasked_match_neither: usize = 0;
    // Total UNMASKED discriminating positions matching ALT, for partial counting.
    // N bases and low-BQ positions are excluded — only high-confidence evidence counts.
    let mut positions_matching_alt: usize = 0;
    // Track whether ANY discriminating position had an N base (for n_count).
    let mut had_n_base = false;
    let mut mnp_quals: Vec<u8> = Vec::with_capacity(len);

    for i in 0..len {
        let pos = start_read_pos + i;
        mnp_quals.push(quals[pos]);
        let base = seq_bytes[pos];
        let base_upper = base.to_ascii_uppercase();
        let ref_upper = ref_bytes[i].to_ascii_uppercase();
        let alt_upper = alt_bytes[i].to_ascii_uppercase();

        // Skip non-discriminating positions (REF == ALT)
        if ref_upper == alt_upper {
            continue;
        }
        n_discriminating += 1;

        // Mask condition: BQ too low OR base is N (uninformative)
        // N bases do NOT contribute to positions_matching_alt — they are
        // uninformative and should not inflate partial counting.
        let is_n = base == b'N' || base == b'n';
        if quals[pos] < min_baseq || is_n {
            n_masked += 1;
            if is_n {
                had_n_base = true;
                trace!("MNP position {} masked: N base (BQ={})", i, quals[pos]);
            }
            continue; // Cannot vote
        }

        // Unmasked position: votes based on allele match
        if base_upper == alt_upper {
            n_unmasked_match_alt += 1;
            // Track unmasked positions matching ALT for partial counting.
            // Only unmasked positions count — N bases are excluded.
            positions_matching_alt += 1;
        } else if base_upper == ref_upper {
            n_unmasked_match_ref += 1;
        } else {
            n_unmasked_match_neither += 1;
        }
    }

    // Safety: an MNP with zero discriminating positions is degenerate
    // (REF == ALT). Treat as ThirdAllele to surface the issue.
    if n_discriminating == 0 {
        warn!(
            "MNP with zero discriminating positions: {}>{} at {}:{}",
            variant.ref_allele, variant.alt_allele, variant.chrom, variant.pos + 1
        );
        return MnpResult::ThirdAllele(0, false);
    }

    let n_unmasked = n_discriminating - n_masked;

    // All discriminating positions masked → LowQuality (contributes to DP only).
    // No unmasked positions can vote, so we can't classify this read.
    if n_unmasked == 0 {
        trace!(
            "MNP all {} discriminating positions masked ({} N, {} low-BQ)",
            n_discriminating, n_masked, n_discriminating
        );
        return MnpResult::LowQuality(positions_matching_alt as u8, had_n_base);
    }

    // Classification based on unmasked votes:
    let med_qual = median_qual(&mnp_quals, min_baseq);

    if n_unmasked_match_ref == n_unmasked && n_unmasked_match_alt == 0 {
        // All unmasked discriminating positions match REF
        MnpResult::Ref(med_qual, had_n_base)
    } else if n_unmasked_match_alt == n_unmasked && n_unmasked_match_ref == 0 {
        // All unmasked discriminating positions match ALT
        MnpResult::Alt(med_qual, had_n_base)
    } else {
        // Mixed or neither — log per-position breakdown for diagnostics
        if log::log_enabled!(log::Level::Trace) {
            let mut pos_info = Vec::with_capacity(n_discriminating);
            for i in 0..len {
                let ref_upper = ref_bytes[i].to_ascii_uppercase();
                let alt_upper = alt_bytes[i].to_ascii_uppercase();
                if ref_upper == alt_upper { continue; }

                let pos = start_read_pos + i;
                let base = seq_bytes[pos];
                let base_upper = base.to_ascii_uppercase();
                let is_masked = quals[pos] < min_baseq || base == b'N' || base == b'n';

                let label = if is_masked {
                    "MASKED"
                } else if base_upper == alt_upper {
                    "ALT"
                } else if base_upper == ref_upper {
                    "REF"
                } else {
                    "other"
                };
                pos_info.push(format!("pos{}:{}={}", i, base_upper as char, label));
            }
            trace!(
                "MNP ThirdAllele {}>{}: {} (unmasked: {} ref, {} alt, {} other; {} masked)",
                variant.ref_allele, variant.alt_allele, pos_info.join(", "),
                n_unmasked_match_ref, n_unmasked_match_alt,
                n_unmasked_match_neither, n_masked
            );
        }
        MnpResult::ThirdAllele(positions_matching_alt as u8, had_n_base)
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
#[allow(clippy::too_many_arguments)]
pub fn check_complex<F: Fn(u8, u8) -> i32>(
    record: &Record,
    variant: &Variant,
    siblings: &[Variant],
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
    // structurally anomalous reads (is_worth_realignment) to alignment-based
    // classification, which extracts raw bases and aligns against full haplotypes.
    if let Some(ref ctx) = variant.ref_context {
        let win_start = variant.ref_context_start;
        let win_end = win_start + ctx.len() as i64;

        if is_worth_realignment(record, win_start, win_end) {
            if let Some((sub_seq, sub_quals)) = extract_raw_read_window(
                record, quals, win_start, win_end, variant.pos, variant.ref_allele.len()
            ) {
                if sub_seq.len() >= 3 {
                    trace!(
                        "check_complex: Phase 0 bypass (soft-clips/indels), extracted {} bases",
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
                            // Pangenomic WFA → marginalized PairHMM pipeline.
                            // Falls back to SW if matrix construction fails.
                            pangenomic_classify(
                                &sub_seq, &sub_quals, variant, siblings,
                                min_baseq, *gap_open, *gap_extend,
                                *gap_open_repeat, *gap_extend_repeat, *llr_threshold,
                            ).unwrap_or_else(|| {
                                trace!(
                                    "check_complex: Phase 0 pangenomic failed, SW fallback at {}:{}",
                                    variant.chrom, variant.pos + 1,
                                );
                                classify_by_alignment(
                                    &sub_seq, &sub_quals, variant, min_baseq,
                                    alt_aligner, ref_aligner,
                                )
                            })
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
                // Inclusive of end_pos (deliberately, unlike the Match/SoftClip
                // branches): an insertion at the exclusive REF end is a *trailing
                // insertion* that belongs to the ALT haplotype (e.g. REF=AB, ALT=ABC).
                // Capturing it lets the reconstruction reach alt_len and match the
                // ALT; excluding it would reconstruct only the REF span and
                // misclassify such reads as REF. A read whose reconstruction equals
                // the ALT genuinely IS ALT evidence — Phase 2 still requires the
                // inserted bases to match the ALT exactly, so this cannot manufacture
                // a false ALT. (Soft-clips differ: clipped bases are unaligned and
                // uncertain, so that branch stays exclusive.)
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

    // N-base detection: scan the reconstructed haplotype for N bases.
    // N bases in the reconstructed sequence are masked by masked_dual_compare /
    // masked_single_compare (they don't affect classification), but we need to
    // propagate has_n_base=true so the engine can increment n_count for
    // duplex masking QC. Matches the MNP/SNP N-tracking pattern.
    let had_n = reconstructed_seq.iter().any(|&b| b == b'N' || b == b'n');
    if had_n {
        trace!(
            "check_complex: N base detected in reconstructed haplotype ({} bases)",
            reconstructed_seq.len()
        );
    }

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
            trace!("No reliable bases — discarding (had_n={})", had_n);
            let mut r = ClassifyResult::neither(ClassifyPhase::MaskedCompare);
            r.has_n_base = had_n;
            return r;
        }

        // Step 2: Ambiguity — reliable bases match both alleles → discard
        if mismatches_alt == 0 && mismatches_ref == 0 {
            trace!("Ambiguous: reliable bases match both REF and ALT — discarding (had_n={})", had_n);
            let mut r = ClassifyResult::neither(ClassifyPhase::MaskedCompare);
            r.has_n_base = had_n;
            return r;
        }

        // Step 3: Unambiguous match
        if mismatches_alt == 0 {
            trace!("Matches ALT on {} reliable bases, med_qual={}, had_n={}", reliable_count, med_haplotype_qual, had_n);
            let mut r = ClassifyResult::is_alt(med_haplotype_qual, ClassifyPhase::MaskedCompare);
            r.has_n_base = had_n;
            return r;
        }
        if mismatches_ref == 0 {
            trace!("Matches REF on {} reliable bases, med_qual={}, had_n={}", reliable_count, med_haplotype_qual, had_n);
            let mut r = ClassifyResult::is_ref(med_haplotype_qual, ClassifyPhase::MaskedCompare);
            r.has_n_base = had_n;
            return r;
        }

        // Step 4: Neither matches cleanly on reliable bases.
        // Report partial ALT evidence if some reliable bases matched ALT.
        let partial_alt_bases = reliable_count.saturating_sub(mismatches_alt);
        if partial_alt_bases > 0 {
            trace!(
                "Case A: {} of {} reliable bases match ALT (partial evidence)",
                partial_alt_bases, reliable_count
            );
            return ClassifyResult::neither_with_partial(
                ClassifyPhase::MaskedCompare,
                partial_alt_bases as u8,
                had_n,
            );
        }
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
            trace!("Matches ALT on {} reliable bases, med_qual={}, had_n={}", reliable_count, med_haplotype_qual, had_n);
            let mut r = ClassifyResult::is_alt(med_haplotype_qual, ClassifyPhase::MaskedCompare);
            r.has_n_base = had_n;
            return r;
        }
        // Case B fall-through: partial ALT evidence if some bases matched
        if reliable_count > 0 && mismatches > 0 {
            let partial_alt_bases = reliable_count.saturating_sub(mismatches);
            if partial_alt_bases > 0 {
                trace!(
                    "Case B: {} of {} reliable bases match ALT (partial evidence)",
                    partial_alt_bases, reliable_count
                );
                return ClassifyResult::neither_with_partial(
                    ClassifyPhase::MaskedCompare,
                    partial_alt_bases as u8,
                    had_n,
                );
            }
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
            trace!("Matches REF on {} reliable bases, med_qual={}, had_n={}", reliable_count, med_haplotype_qual, had_n);
            let mut r = ClassifyResult::is_ref(med_haplotype_qual, ClassifyPhase::MaskedCompare);
            r.has_n_base = had_n;
            return r;
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
                trace!("Phase 2.5 → ALT (edit distance margin, had_n={})", had_n);
                let mut r = ClassifyResult::is_alt(med_haplotype_qual, ClassifyPhase::Levenshtein);
                r.has_n_base = had_n;
                return r;
            } else if d_ref + 1 < d_alt {
                trace!("Phase 2.5 → REF (edit distance margin, had_n={})", had_n);
                let mut r = ClassifyResult::is_ref(med_haplotype_qual, ClassifyPhase::Levenshtein);
                r.has_n_base = had_n;
                // Propagate nearby evidence: Levenshtein classified as REF, but
                // if ALT edit distance is close (within 3 edits of REF), the read
                // has partial ALT evidence worth tracking for diagnostics.
                if d_alt <= d_ref + 3 {
                    r.has_nearby_evidence = true;
                    trace!(
                        "Phase 2.5: ALT edit distance {} close to REF {} → has_nearby_evidence=true",
                        d_alt, d_ref
                    );
                }
                return r;
            }
            // else: ambiguous, fall through to Phase 3
        }
    }

    // --- Phase 3: Alignment-based fallback (indelpost approach) ---
    // When narrow-window reconstruction fails (FM1: D-truncation, FM2: adjacent I),
    // expand to the full ref_context window and use alignment-based classification.
    //
    // CRITICAL: Use raw read window extraction (not CIGAR-projected) to preserve
    // the true biological sequence. For complex variants (e.g. EPHA7 REF=TCC ALT=CT),
    // BWA represents ALT reads as DEL+INS CIGARs. CIGAR-projected extraction
    // produces a hybrid sequence matching neither REF nor ALT haplotype.
    // Raw extraction gives the contiguous read bases that alignment can correctly classify.
    //
    // For PairHMM backend: uses pangenomic WFA → marginalized PairHMM pipeline,
    // which evaluates the read against ALL sibling haplotypes, not just H0/H1.
    //
    // Pre-filter (indelpost pattern): only attempt alignment for reads showing
    // evidence of carrying the variant (soft-clips, indels near window).
    // This eliminates ~80-90% of clean REF reads from expensive alignment.
    if let Some(ref ctx) = variant.ref_context {
        let win_start = variant.ref_context_start;
        let win_end = win_start + ctx.len() as i64;

        let is_mnp = variant.ref_allele.len() == variant.alt_allele.len() && variant.ref_allele.len() > 1;

        if !is_mnp && !is_worth_realignment(record, win_start, win_end) {
            // Clean CIGAR over the window: no indels or soft-clips near the
            // variant, so expensive realignment is not needed.
            //
            // For deletion-direction complex alleles (ref_len > alt_len, e.g.
            // a 100bp deletion with 7bp replacement like NF2), REF reads have
            // no CIGAR deletion and land here. Returning `neither` for them
            // causes ref=0 even when 56 true REF reads are visible in IGV.
            //
            // Fix: if the variant removes bases (ref_len > alt_len) and the
            // read's M-blocks cover the anchor, classify the read as REF.
            // This mirrors the `found_ref_coverage → REF` logic in
            // check_deletion, but applied here for complex alleles that arrive
            // via the `else` routing (both ref_len > 1 and alt_len > 1).
            //
            // Precondition: is_worth_realignment is false, so the read has a
            // clean CIGAR with no nearby indels — we can trust M coverage.
            if ref_bytes.len() > alt_bytes.len() {
                let anchor_pos = variant.pos;
                let mut rpos = record.pos();
                let mut anchor_qual: Option<u8> = None;
                'cigar_walk: for op in record.cigar().iter() {
                    match op {
                        rust_htslib::bam::record::Cigar::Match(len)
                        | rust_htslib::bam::record::Cigar::Equal(len)
                        | rust_htslib::bam::record::Cigar::Diff(len) => {
                            let block_end = rpos + *len as i64;
                            if anchor_pos >= rpos && anchor_pos < block_end {
                                // Compute read position for this anchor
                                let qp = {
                                    let mut q_off = 0usize;
                                    let mut r_off = record.pos();
                                    for op2 in record.cigar().iter() {
                                        match op2 {
                                            rust_htslib::bam::record::Cigar::Match(l)
                                            | rust_htslib::bam::record::Cigar::Equal(l)
                                            | rust_htslib::bam::record::Cigar::Diff(l) => {
                                                if anchor_pos >= r_off && anchor_pos < r_off + *l as i64 {
                                                    q_off += (anchor_pos - r_off) as usize;
                                                    break;
                                                }
                                                r_off += *l as i64;
                                                q_off += *l as usize;
                                            }
                                            rust_htslib::bam::record::Cigar::Del(l)
                                            | rust_htslib::bam::record::Cigar::RefSkip(l) => {
                                                r_off += *l as i64;
                                            }
                                            rust_htslib::bam::record::Cigar::Ins(l)
                                            | rust_htslib::bam::record::Cigar::SoftClip(l)
                                            | rust_htslib::bam::record::Cigar::HardClip(l) => {
                                                q_off += *l as usize;
                                            }
                                            _ => {}
                                        }
                                    }
                                    q_off
                                };
                                let qual_val = if qp < quals.len() { quals[qp] } else { 0 };
                                if qual_val >= min_baseq {
                                    anchor_qual = Some(qual_val);
                                }
                                break 'cigar_walk;
                            }
                            rpos = block_end;
                        }
                        rust_htslib::bam::record::Cigar::Del(len)
                        | rust_htslib::bam::record::Cigar::RefSkip(len) => {
                            rpos += *len as i64;
                        }
                        _ => {}
                    }
                }
                if let Some(aq) = anchor_qual {
                    trace!(
                        "check_complex: clean M covers anchor {} for del-direction complex \
                         allele (ref_len={} > alt_len={}) → REF",
                        anchor_pos,
                        ref_bytes.len(),
                        alt_bytes.len(),
                    );
                    return ClassifyResult::is_ref(aq, ClassifyPhase::Alignment);
                }
            }

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
                        // Pangenomic WFA → marginalized PairHMM pipeline.
                        // Falls back to SW if matrix construction fails.
                        pangenomic_classify(
                            &sub_seq, &sub_quals, variant, siblings,
                            min_baseq, *gap_open, *gap_extend,
                            *gap_open_repeat, *gap_extend_repeat, *llr_threshold,
                        ).unwrap_or_else(|| {
                            trace!(
                                "check_complex: Phase 3 pangenomic failed, SW fallback at {}:{}",
                                variant.chrom, variant.pos + 1,
                            );
                            classify_by_alignment(
                                &sub_seq, &sub_quals, variant, min_baseq,
                                alt_aligner, ref_aligner,
                            )
                        })
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
#[allow(clippy::too_many_arguments)]
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

    // VCF/MAF left-anchored invariant — REF and ALT share a leading anchor
    // base, so both are non-empty for a well-formed insertion. Defend it: an empty
    // ALT/REF (malformed or non-left-anchored record) would underflow `len() - 1`
    // and panic the `[1..]` / `[0]` slices below, surfacing as an opaque PyErr.
    // debug! not warn!: this runs per-read, so a single malformed variant would
    // otherwise log once per overlapping read. Loud once-per-variant surfacing
    // belongs at prep time (tracked follow-up); here we just avoid the panic.
    if variant.alt_allele.is_empty() || variant.ref_allele.is_empty() {
        debug!(
            "check_insertion: empty REF/ALT at {}:{} (ref={:?} alt={:?}) — classifying neither",
            variant.chrom,
            variant.pos + 1,
            variant.ref_allele,
            variant.alt_allele,
        );
        return ClassifyResult::neither(ClassifyPhase::Structural);
    }
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
    let mut has_nearby_length_match = false; // nearby Ins needs Phase 3: wrong-seq or wrong-length

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
                                        "check_insertion: backward boundary match at pos {}, qual={} (structural)",
                                        anchor_pos, qual
                                    );
                                    return ClassifyResult::is_alt_structural(qual, ClassifyPhase::Structural); // ALT — backward match
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
                                            "check_insertion: strict match at pos {}, anchor_qual={} (structural)",
                                            anchor_pos, qual
                                        );
                                        return ClassifyResult::is_alt_structural(qual, ClassifyPhase::Structural); // ALT — strict match
                                    }
                                }
                            } else {
                                // Wrong-length insertion at anchor: I(n) where n ≠ expected.
                                // Mirrors check_deletion's strict wrong-length handling:
                                // route to Phase 3 (SW/PairHMM) for haplotype-level
                                // arbitration. The read has an insertion at the exact
                                // anchor but of a different length — structural evidence
                                // of a third allele (e.g., PAX5 I(1) when expecting I(2)).
                                //
                                // Unlike deletions, insertions are point events in reference
                                // space, so reciprocal overlap is not applicable — different
                                // lengths genuinely indicate different alleles, not alignment
                                // breakpoint ambiguity. Always fall back to Phase 3.
                                let found_ins_len = ins_len_usize;
                                trace!(
                                    "check_insertion: I({}) at anchor {} but expected I({}), \
                                     falling back to phase3_classify",
                                    found_ins_len, anchor_pos, expected_ins_len
                                );
                                let mut result = phase3_classify(
                                    record, variant, siblings, quals, min_baseq,
                                    alt_aligner, ref_aligner, backend,
                                );
                                // Propagate nearby evidence: the read has an insertion at
                                // the variant anchor (just wrong length), which is structural
                                // evidence worth tracking even if Phase 3 returns REF.
                                // Consumed by engine to increment partial_alt, enabling the
                                // PARTIAL_DOMINANT diagnostic flag.
                                if !result.is_alt {
                                    result.has_nearby_evidence = true;
                                    trace!(
                                        "check_insertion: Phase 3 did not confirm ALT, but I({}) \
                                         at anchor exists → has_nearby_evidence=true",
                                        found_ins_len
                                    );
                                }
                                return result;
                            }
                        }
                        // Anchor at end but no insertion at all → REF coverage
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
                        } else {
                            // Different-length insertion in window: I(n) where
                            // n ≠ expected. Flag for Phase 3 fallback so the
                            // post-walk handler can route to
                            // haplotype alignment and set has_nearby_evidence.
                            //
                            // Note: has_nearby_length_match is reused here despite
                            // the name implying "same-length" — it means "needs
                            // Phase 3 arbitration because an insertion exists nearby".
                            // Both wrong-sequence and wrong-length cases get the same
                            // post-walk treatment: Phase 3 + has_nearby_evidence.
                            //
                            // This also covers backward boundary insertions
                            // (anchor_pos == ref_pos): the same insertion is at
                            // block_end of the previous M block, which the windowed
                            // scan processes on the prior loop iteration.
                            has_nearby_length_match = true;
                            trace!(
                                "check_insertion: windowed I({}) at pos {} (expected I({})), \
                                 different-length → flagging for Phase 3 fallback",
                                ins_len_usize, ins_ref_pos, expected_ins_len
                            );
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
        return ClassifyResult::is_alt_structural(anchor_qual, ClassifyPhase::CigarRecon); // ALT — windowed INS match (structural)
    }

    // Phase 3 haplotype fallback: when a nearby insertion exists but doesn't match
    // the expected variant — either wrong sequence (same length, S3 failed) or wrong
    // length (different allele). In both cases, the aligner placed an insertion near
    // the anchor that warrants haplotype-level arbitration. Route to Phase 3
    // (Smith-Waterman/PairHMM) for full comparison, and propagate has_nearby_evidence
    // if Phase 3 doesn't confirm ALT, so the engine counts it as partial_alt.
    if has_nearby_length_match && found_ref_coverage {
        trace!(
            "check_insertion: nearby insertion evidence at pos {} \
             (wrong-seq or wrong-length), falling back to phase3_classify",
            anchor_pos
        );
        let mut result = phase3_classify(record, variant, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend);
        // Propagate nearby evidence: Phase 3 may return is_ref or neither,
        // but the CIGAR proved an insertion exists nearby. Mark it so the
        // engine can count this read as partial_alt evidence, enabling the
        // PARTIAL_DOMINANT diagnostic flag.
        if !result.is_alt {
            result.has_nearby_evidence = true;
            trace!(
                "check_insertion: Phase 3 did not confirm ALT, but nearby insertion \
                 exists at pos {} → has_nearby_evidence=true",
                anchor_pos
            );
        }
        return result;
    }

    // P0-3: Haplotype fallback — when strict/windowed CIGAR matching found no
    // insertion match and the read doesn't cover the anchor, try Phase 3
    // (check_complex → Smith-Waterman/PairHMM) for haplotype-level comparison.
    // Only fall back when NOT found_ref_coverage to avoid false positives on
    // reads that genuinely show REF at this position.
    //
    // Mirrors check_deletion's !found_ref_coverage path.
    //
    // CRITICAL: Only attempt Phase 3 if the read actually overlaps the anchor
    // position (variant.pos). Reads that don't span the anchor have no
    // information about the variant and must not be counted.
    if !found_ref_coverage && best_windowed_match.is_none() {
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

        if record.pos() <= anchor_pos && read_ref_end > anchor_pos {
            trace!(
                "check_insertion: no CIGAR evidence at pos {}, read spans anchor \
                 ({}..{}), falling back to phase3_classify",
                anchor_pos, record.pos(), read_ref_end
            );
            return phase3_classify(
                record, variant, siblings, quals, min_baseq,
                alt_aligner, ref_aligner, backend,
            );
        }
        // Otherwise: read doesn't overlap the anchor → no variant info
    }

    if found_ref_coverage {
        // CIGAR is definitive for pure insertions: if the read's Match op
        // covers the anchor and no matching I op was found (strict or
        // windowed), the read is REF. Soft-clipping elsewhere cannot
        // represent a missing insertion at the anchor position.
        //
        // Note: reads with soft-clip AT the anchor have
        // found_ref_coverage=false and are correctly routed to Phase 3
        // via the !found_ref_coverage haplotype fallback above.
        return ClassifyResult::is_ref(anchor_qual, ClassifyPhase::Structural);
    }
    ClassifyResult::neither(ClassifyPhase::Structural) // Read does not cover the variant region
}

/// Verify that the reference bases at an observed deletion position match the
/// variant's expected deleted bases over `compare_len` positions.
///
/// Used by `check_deletion`'s Safeguard 3 for exact-length matches AND for
/// tolerant (different-length) matches over the *shared* span. Comparing the
/// reference at the observed breakpoint against `expected_del_seq` rejects an
/// unrelated SV that merely shares ≥50% length overlap with the target deletion,
/// while still accepting the same deletion reported with a slightly different
/// breakpoint length (its overlapping prefix still matches). Both operands are
/// reference-derived, so a true match is exact. Returns `true` only on a reliable,
/// zero-mismatch concordance; a `None` ref_context or out-of-bounds offset returns
/// `false` (the caller routes those to the Phase 3 haplotype fallback).
fn verify_deleted_bases(
    variant: &Variant,
    del_ref_pos: i64,
    expected_del_seq: &[u8],
    compare_len: usize,
) -> bool {
    let ctx = match &variant.ref_context {
        Some(c) => c.as_bytes(),
        None => {
            warn!(
                "check_deletion: ref_context is None at {}:{} — cannot validate deleted bases",
                variant.chrom,
                variant.pos + 1,
            );
            return false;
        }
    };
    let ctx_offset_i64 = del_ref_pos - variant.ref_context_start;
    if ctx_offset_i64 < 0 {
        trace!(
            "verify_deleted_bases: negative ref_context offset ({}), rejecting",
            ctx_offset_i64
        );
        return false;
    }
    let ctx_offset = ctx_offset_i64 as usize;
    if compare_len == 0
        || compare_len > expected_del_seq.len()
        || ctx_offset + compare_len > ctx.len()
    {
        trace!(
            "verify_deleted_bases: out of bounds (offset={} len={} ctx={} exp={}), rejecting",
            ctx_offset,
            compare_len,
            ctx.len(),
            expected_del_seq.len(),
        );
        return false;
    }
    let ref_at = &ctx[ctx_offset..ctx_offset + compare_len];
    let ref_quals = vec![u8::MAX; compare_len];
    let (mismatches, reliable) =
        masked_single_compare(ref_at, &ref_quals, &expected_del_seq[..compare_len], 0);
    reliable > 0 && mismatches == 0
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
#[allow(clippy::too_many_arguments)]
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

    // Left-anchored invariant — an empty REF/ALT (malformed or
    // non-left-anchored record) would underflow `len() - 1` and panic the
    // `[1..]` slice below. Defend it and classify as neither. debug! not warn!:
    // runs per-read; loud once-per-variant surfacing belongs at prep (follow-up).
    if variant.ref_allele.is_empty() || variant.alt_allele.is_empty() {
        debug!(
            "check_deletion: empty REF/ALT at {}:{} (ref={:?} alt={:?}) — classifying neither",
            variant.chrom,
            variant.pos + 1,
            variant.ref_allele,
            variant.alt_allele,
        );
        return ClassifyResult::neither(ClassifyPhase::Structural);
    }
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
    // Tracks a windowed Del with matching length but failed S3 sequence check.
    // When set alongside found_ref_coverage, the read carries a same-length deletion
    // placed at a different position by the aligner (often due to left-alignment
    // in a repeat context) — Phase 3 SW can arbitrate correctly.
    let mut has_nearby_length_match = false; // nearby Del needs Phase 3: wrong-seq or wrong-length


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
                                    "check_deletion: strict match at pos {}, anchor_qual={} (structural)",
                                    anchor_pos, qual
                                );
                                return ClassifyResult::is_alt_structural(qual, ClassifyPhase::Structural); // ALT — strict match
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
                                    // A ≥50% length overlap alone is not enough —
                                    // verify the OVERLAPPING deleted bases match
                                    // expected_del_seq, so an unrelated SV at the anchor
                                    // is not accepted. The D begins at anchor_pos + 1.
                                    let compare_len = found_del_len.min(expected_del_len);
                                    if verify_deleted_bases(
                                        variant,
                                        anchor_pos + 1,
                                        expected_del_seq,
                                        compare_len,
                                    ) {
                                        let arp = anchor_read_pos.unwrap_or(0);
                                        let qual = if arp < quals.len() { quals[arp] } else { 0 };
                                        trace!(
                                            "check_deletion: tolerant match D({}) ≈ D({}) at pos {} \
                                             (reciprocal_overlap={:.2}, seq-verified over {}bp, \
                                             anchor_qual={}) (structural)",
                                            found_del_len,
                                            expected_del_len,
                                            anchor_pos,
                                            reciprocal_overlap,
                                            compare_len,
                                            qual,
                                        );
                                        return ClassifyResult::is_alt_structural(qual, ClassifyPhase::Structural); // ALT — tolerant, seq-verified
                                    }
                                    // Overlap met but the deleted bases don't match the
                                    // target deletion → likely a different SV. Fall
                                    // through to check_complex for haplotype arbitration.
                                    trace!(
                                        "check_deletion: tolerant structural D({}) at anchor {} REJECTED \
                                         — overlapping deleted bases mismatch expected over {}bp → check_complex",
                                        found_del_len,
                                        anchor_pos,
                                        compare_len,
                                    );
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
                                // Large deletion (≥50bp) but low reciprocal overlap
                                // (<50%): significantly different deletion events.
                                // Flag for Phase 3 haplotype arbitration.
                                has_nearby_length_match = true;
                                trace!(
                                    "check_deletion: windowed D({}) at pos {} (expected D({})), \
                                     low reciprocal overlap ({:.2} < 0.50) → flagging for Phase 3 fallback",
                                    del_len_usize, del_ref_pos, expected_del_len, overlap
                                );
                                false
                            }
                        } else {
                            // Small wrong-length deletion in window (< 50bp,
                            // reciprocal overlap n/a). Flag for Phase 3 fallback
                            // for dels ≥ 5bp — short (1-4bp) wrong-length
                            // deletions are almost certainly spurious noise in
                            // homopolymer/STR regions; CIGAR is definitive for
                            // those. Mirrors insertion windowed fix (Step 1.3).
                            if del_len_usize >= 5 {
                                has_nearby_length_match = true;
                                trace!(
                                    "check_deletion: windowed D({}) at pos {} (expected D({})), \
                                     different-length (≥5bp) → flagging for Phase 3 fallback",
                                    del_len_usize, del_ref_pos, expected_del_len
                                );
                            } else {
                                trace!(
                                    "check_deletion: windowed D({}) at pos {} (expected D({})), \
                                     different-length (<5bp) → CIGAR definitive, not flagging",
                                    del_len_usize, del_ref_pos, expected_del_len
                                );
                            }
                            false
                        };

                        if length_ok {
                            // Safeguard 3: verify the deleted reference bases match
                            // expected_del_seq. Exact-length matches compare the full
                            // span; tolerant (different-length) matches compare
                            // the OVERLAPPING span instead of skipping verification, so
                            // an unrelated SV sharing only ≥50% length overlap is not
                            // accepted as ALT. A genuinely shifted/different deletion
                            // fails here and is routed to the Phase 3 fallback below.
                            let compare_len = del_len_usize.min(expected_del_len);
                            let del_ok =
                                verify_deleted_bases(variant, del_ref_pos, expected_del_seq, compare_len);

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
                                // S3 failed: deleted bases at the shifted position
                                // don't match expected_del_seq. This happens when
                                // left-alignment moves the anchor further left than
                                // the aligner's CIGAR Del position, so the reference
                                // slice at del_ref_pos encodes different bases than
                                // expected_del_seq (e.g. TP53 GACCGTGCAAGT→- where
                                // left-alignment shifts anchor 3bp left of the actual
                                // D(12) in reads). Track for Phase 3 fallback.
                                //
                                // Only flag for Phase 3 when del_len >= 5bp.
                                // Short (1-4bp) same-length deletions that fail S3
                                // are almost certainly spurious/unrelated noise —
                                // CIGAR remains definitive for them. Longer deletions
                                // are more susceptible to BWA left-alignment shifting
                                // the anchor multiple positions away from the CIGAR D.
                                if del_len_usize >= 5 {
                                    has_nearby_length_match = true;
                                    trace!(
                                        "check_deletion: S3 reject at shifted pos {} \
                                         (deleted bases mismatch, del_len={} >= 5), \
                                         flagging for Phase 3 fallback",
                                        del_ref_pos, del_len_usize
                                    );
                                } else {
                                    trace!(
                                        "check_deletion: S3 reject at shifted pos {} \
                                         (deleted bases mismatch, del_len={} < 5, \
                                         CIGAR definitive — not flagging Phase 3)",
                                        del_ref_pos, del_len_usize
                                    );
                                }
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
        return ClassifyResult::is_alt_structural(anchor_qual, ClassifyPhase::CigarRecon); // ALT — windowed DEL match (structural)
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

    // Phase 3 haplotype fallback: when a nearby deletion exists but doesn't match
    // the expected variant — either wrong sequence (same length, S3 failed), wrong
    // length (different allele ≥5bp), or large deletion with low reciprocal overlap
    // (≥50bp, <50%). In all cases, the aligner placed a deletion near the anchor
    // that warrants haplotype-level arbitration. Route to Phase 3 for full
    // comparison, and propagate has_nearby_evidence if Phase 3 doesn't confirm ALT.
    if has_nearby_length_match && found_ref_coverage {
        trace!(
            "check_deletion: nearby deletion evidence at pos {} \
             (wrong-seq or wrong-length), falling back to phase3_classify",
            anchor_pos
        );
        let mut result = phase3_classify(record, variant, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend);
        // Propagate nearby evidence: Phase 3 may return is_ref or neither,
        // but the CIGAR proved a deletion exists nearby. Mark it so the
        // engine can count this read as partial_alt evidence, enabling the
        // PARTIAL_DOMINANT diagnostic flag.
        if !result.is_alt {
            result.has_nearby_evidence = true;
            trace!(
                "check_deletion: Phase 3 did not confirm ALT, but nearby deletion \
                 exists at pos {} → has_nearby_evidence=true",
                anchor_pos
            );
        }
        return result;
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

#[cfg(test)]
mod tests {
    use super::*;

    /// A deletion Variant carrying a ref_context window for the seq-check tests.
    fn deletion_with_context(ref_context: &str, ref_context_start: i64, ref_allele: &str) -> Variant {
        Variant {
            chrom: "1".to_string(),
            pos: ref_context_start,
            ref_allele: ref_allele.to_string(),
            alt_allele: ref_allele.get(..1).unwrap_or("N").to_string(),
            variant_type: "DELETION".to_string(),
            ref_context: Some(ref_context.to_string()),
            ref_context_start,
            repeat_span: 0,
            gene_strand: None,
        }
    }

    // ── verify_deleted_bases — partial seq-check for tolerant deletions ──
    // ref_context "NNNCGTACGTCCCC" at genomic 100 → pos 103 = "CGTACGT...".

    #[test]
    fn test_verify_deleted_bases_exact_match() {
        let v = deletion_with_context("NNNCGTACGTCCCC", 100, "ACGTACGT");
        assert!(verify_deleted_bases(&v, 103, b"CGTACGT", 7));
    }

    #[test]
    fn test_verify_deleted_bases_tolerant_same_event_matches() {
        // Same start, observed deletion longer than expected → compare the SHARED
        // prefix only; the same biological event still verifies (sensitivity kept).
        let v = deletion_with_context("NNNCGTACGTCCCC", 100, "ACGTACGT");
        assert!(verify_deleted_bases(&v, 103, b"CGTACGT", 7)); // shared span = 7
    }

    #[test]
    fn test_verify_deleted_bases_shifted_different_sv_rejected() {
        // A deletion at a SHIFTED position (105) whose reference bases differ
        // from the target's expected_del_seq is rejected — a different SV that
        // merely shares length overlap must not be counted as ALT.
        let v = deletion_with_context("NNNCGTACGTCCCC", 100, "ACGTACGT");
        // ref at 105 = "TACGTCC" != expected "CGTACGT"
        assert!(!verify_deleted_bases(&v, 105, b"CGTACGT", 7));
    }

    #[test]
    fn test_verify_deleted_bases_none_context_rejects() {
        let mut v = deletion_with_context("NNNCGTACGTCCCC", 100, "ACGTACGT");
        v.ref_context = None;
        assert!(!verify_deleted_bases(&v, 103, b"CGTACGT", 7));
    }

    #[test]
    fn test_verify_deleted_bases_out_of_bounds_rejects() {
        let v = deletion_with_context("NNNCGT", 100, "ACGT");
        assert!(!verify_deleted_bases(&v, 103, b"CGTACGT", 7)); // past ref_context end
        assert!(!verify_deleted_bases(&v, 50, b"CGT", 3)); // negative offset
    }
}
