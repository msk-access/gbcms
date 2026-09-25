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

use crate::normalize::repeat::find_tandem_repeat;
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
/// Routes to `check_complex` (SW backend) or the pangenomic WFA+PairHMM
/// pipeline (PairHMM backend). Called from the Phase 3 fallback sites in
/// check_insertion and check_deletion.
///
/// For PairHMM backend:
/// 1. Build pangenomic haplotype matrix (H0/H1/H2..H2n) from variant + siblings
/// 2. Try WFA fast-path triage (edit distance; resolves ~70-80% of reads)
/// 3. Ambiguous reads fall through to marginalized PairHMM (BQ-aware)
/// 4. If the pangenomic pipeline cannot run (no ref_context, window
///    extraction failure, short window, matrix build failure), falls back to
///    `check_complex`'s full pipeline — where Smith-Waterman itself runs only
///    if the matrix cannot be built there either
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
                        // fall through to the check_complex fallback below.
                    } else {
                        trace!(
                            "phase3: sub_seq too short ({} < 3) → check_complex fallback at {}:{}",
                            sub_seq.len(), variant.chrom, variant.pos + 1,
                        );
                    }
                } else {
                    trace!(
                        "phase3: read window extraction failed → check_complex fallback at {}:{}",
                        variant.chrom, variant.pos + 1,
                    );
                }
            } else {
                debug!(
                    "phase3: variant {}:{} has no ref_context → check_complex fallback",
                    variant.chrom, variant.pos + 1,
                );
            }
            // Fallback when the pangenomic pipeline cannot run here: check_complex's
            // FULL pipeline (structural bypass → CIGAR reconstruction → masked
            // comparison → Levenshtein → its own Phase 3). Under this PairHMM
            // backend, check_complex's Phase 3 retries the pangenomic pipeline;
            // actual Smith-Waterman runs only if the haplotype matrix cannot be
            // built there either.
            trace!(
                "phase3: falling back to check_complex full pipeline for {}:{}",
                variant.chrom, variant.pos + 1,
            );
            check_complex(record, variant, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend)
        }
    }
}


/// How a read's CIGAR observes a reference interval `[start, end)`.
///
/// - `aligned`: ≥1 M/=/X base inside the interval — a direct observation.
/// - `deleted`: ≥1 D op overlapping it — the aligner asserts those
///   reference bases are absent from the molecule (deletion evidence).
/// - `skipped`: ≥1 N op overlapping it — the aligner asserts the interval
///   is spliced-out reference; the read observes nothing there.
struct SpanObservation {
    aligned: bool,
    deleted: bool,
    skipped: bool,
}

fn observe_read_span(record: &Record, start: i64, end: i64) -> SpanObservation {
    let mut obs = SpanObservation { aligned: false, deleted: false, skipped: false };
    let mut ref_pos = record.pos();
    for op in record.cigar().iter() {
        // Zero-length ops (BAM-representable) have an empty reference
        // interval [x, x), which would still pass the half-open overlap test
        // at interior points — they observe nothing and must not set flags.
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let block_end = ref_pos + *len as i64;
                if *len > 0 && ref_pos < end && block_end > start {
                    obs.aligned = true;
                }
                ref_pos = block_end;
            }
            Cigar::Del(len) => {
                let block_end = ref_pos + *len as i64;
                if *len > 0 && ref_pos < end && block_end > start {
                    obs.deleted = true;
                }
                ref_pos = block_end;
            }
            Cigar::RefSkip(len) => {
                let block_end = ref_pos + *len as i64;
                if *len > 0 && ref_pos < end && block_end > start {
                    obs.skipped = true;
                }
                ref_pos = block_end;
            }
            // I/S/H/P consume no reference — nothing to observe in ref space.
            _ => {}
        }
        // Ref-consuming ops are ordered, so nothing later can reach the interval.
        if ref_pos >= end {
            break;
        }
    }
    obs
}

/// Splice-skip triage: decide, before per-type classification, whether a
/// spliced read observes this variant's discriminating positions at all.
///
/// A CIGAR `N` (RefSkip) asserts spliced-out reference. Unlike `D`, it is
/// not a claim that the molecule lacks those bases — it is a claim that the
/// read observes nothing there (samtools pileup reports zero coverage inside
/// an N gap). A read whose N spans every discriminating position therefore
/// cannot distinguish REF from ALT and must not testify — or count depth —
/// in either direction. Without this gate such reads fall through the
/// checkers' anchor-based fast paths as definitive REF, so the REF/ALT
/// ledger flips on the aligner's arbitrary D-vs-N representation choice.
///
/// Returns `Some(no_coverage)` when the read has no observation:
/// - no aligned (M/=/X) base on any discriminating position, AND
/// - no D op overlapping them (a D there is deletion evidence and proceeds
///   to normal classification), AND
/// - at least one discriminating position inside an N op.
///
/// Discriminating positions by variant geometry (0-based, half-open):
/// - pure deletion (anchor preserved): the deleted span
///   `[pos+1, pos+ref_len)` — the retained anchor base alone cannot
///   distinguish the alleles.
/// - pure insertion (single-base REF): the two bases flanking the insertion
///   junction, `[pos, pos+2)`. A read spliced over both flanks has nothing
///   to say about bases inserted between them. A read whose M ends exactly
///   at the junction keeps its aligned anchor and classifies normally.
/// - SNV / MNP / delins / anchor-substituting Del+SNV: every REF position,
///   `[pos, pos+ref_len)`.
///
/// Reads without any N op return `None` immediately, so DNA-mode
/// classification is untouched.
pub fn splice_skip_triage(record: &Record, variant: &Variant) -> Option<ClassifyResult> {
    if !super::rna::has_splice_junction(record) {
        return None;
    }
    let ref_len = variant.ref_allele.len() as i64;
    let alt_len = variant.alt_allele.len() as i64;
    if ref_len == 0 || alt_len == 0 {
        return None; // malformed alleles — the checkers' own guards handle these
    }
    let anchor_preserved = variant
        .ref_allele
        .as_bytes()
        .first()
        .map(|b| b.to_ascii_uppercase())
        == variant
            .alt_allele
            .as_bytes()
            .first()
            .map(|b| b.to_ascii_uppercase());
    let (span_start, span_end) = if ref_len == 1 && alt_len > 1 {
        (variant.pos, variant.pos + 2)
    } else if ref_len > 1 && alt_len == 1 && anchor_preserved {
        (variant.pos + 1, variant.pos + ref_len)
    } else {
        (variant.pos, variant.pos + ref_len)
    };
    let obs = observe_read_span(record, span_start, span_end);
    if obs.skipped && !obs.aligned && !obs.deleted {
        // Shifted-representation escape: the aligner can place the event's
        // I/D op just outside the annotated span (repeat-tract shift) while
        // extending the N over the span itself. Such a read still carries
        // indel evidence — the windowed scans (including the post-N
        // inspection) exist precisely to judge it, with their own sequence
        // safeguards. Exclusion is only for reads with no indel op anywhere
        // in the scan window: a clean junction read observes nothing.
        let window: i64 = std::cmp::max(5, variant.repeat_span as i64 + 2);
        if has_indel_in_window(record, variant.pos - window, variant.pos + window) {
            trace!(
                "splice_skip_triage: N spans [{}, {}) at {}:{} but the read has an \
                 indel op inside the scan window — deferring to classification",
                span_start,
                span_end,
                variant.chrom,
                variant.pos + 1
            );
            return None;
        }
        trace!(
            "splice_skip_triage: N spans [{}, {}) at {}:{} with no aligned base, \
             no D, and no indel op in the scan window — read observes nothing \
             here → no coverage",
            span_start,
            span_end,
            variant.chrom,
            variant.pos + 1
        );
        return Some(ClassifyResult::no_coverage(ClassifyPhase::Structural));
    }
    None
}

/// Whether the read carries any I/D op whose reference start lies inside
/// `[wstart, wend]` — the same candidate-position bounds the windowed scans
/// use. Used by `splice_skip_triage` to keep shifted-representation carriers
/// in classification.
fn has_indel_in_window(record: &Record, wstart: i64, wend: i64) -> bool {
    let mut ref_pos = record.pos();
    for op in record.cigar().iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) | Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            }
            Cigar::Del(len) => {
                if *len > 0 && ref_pos >= wstart && ref_pos <= wend {
                    return true;
                }
                ref_pos += *len as i64;
            }
            Cigar::Ins(len) if *len > 0 && ref_pos >= wstart && ref_pos <= wend => {
                return true;
            }
            _ => {}
        }
        if ref_pos > wend {
            return false;
        }
    }
    false
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
    /// All unmasked bases match ALT. (quality, had_n_at_any_position,
    /// confirmed) — `confirmed` when no discriminating base was masked, i.e.
    /// every one was read and matched ALT: the read itself shows the whole
    /// haplotype rather than inferring it from the unmasked subset.
    Alt(u8, bool, bool),
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
        MnpResult::Alt(med_qual, had_n_base, n_masked == 0)
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


/// A read's CIGAR-projected reconstruction over a genomic span.
pub(crate) struct SpanRecon {
    /// The bases the read shows for [start_pos, end_pos), with insertions at
    /// or inside the span included (trailing insertions deliberately so —
    /// see the Ins branch) and deletions consuming reference silently.
    pub seq: Vec<u8>,
    /// Per-base qualities aligned with `seq`.
    pub quals: Vec<u8>,
    /// A CIGAR `N` (RefSkip) overlapped the span: the reconstruction would
    /// stitch exon arms together and masquerade as deletion evidence, so
    /// callers must refuse to string-compare when this is set.
    pub splice_skip: bool,
}

/// Walk the CIGAR and reconstruct what the read shows for [start_pos,
/// end_pos). Shared by `check_complex` Phase 1 and the engine's
/// multi-allelic AD-claiming contest (span-explanation cost).
pub(crate) fn reconstruct_span(
    record: &Record,
    quals: &[u8],
    start_pos: i64,
    end_pos: i64,
) -> SpanRecon {
    let cigar = record.cigar();
    let mut ref_pos = record.pos();
    let mut read_pos: usize = 0;
    let seq = record.seq();
    let mut reconstructed_seq: Vec<u8> = Vec::with_capacity(seq.len());
    let mut quals_per_base: Vec<u8> = Vec::with_capacity(seq.len());
    let mut splice_skip_in_window = false;
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
            Cigar::Del(len) => {
                ref_pos += *len as i64;
            }
            Cigar::RefSkip(len) => {
                let block_end = ref_pos + *len as i64;
                // Zero-length N: empty interval, observes/skips nothing.
                if *len > 0 && ref_pos < end_pos && block_end > start_pos {
                    splice_skip_in_window = true;
                }
                ref_pos = block_end;
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
    SpanRecon { seq: reconstructed_seq, quals: quals_per_base, splice_skip: splice_skip_in_window }
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

    // NOTE: quals is passed from the caller — either raw record.qual() or BAQ-adjusted.
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
            // A structurally anomalous read whose splice N overlaps the
            // context window cannot be alignment-classified (extraction
            // refuses to stitch across the splice) — and letting it fall
            // through to Phase 1/2 string comparison is exactly the false-REF
            // path this Phase-0 bypass exists to prevent: with its indel
            // shifted outside the variant span, the reconstruction over the
            // span is clean REF sequence and Phase 2 absorbs an ALT carrier
            // into rd. No clean evidence exists for such a read → neither.
            // (Splice-aware Phase-3 scoring would need a spliced haplotype
            // with an explicit genomic→spliced coordinate map and
            // junction-compatible extraction — tracked in issue #94; until
            // real-data measurement justifies that machinery, this stays
            // conservative.)
            if super::rna::has_splice_junction(record)
                && observe_read_span(record, win_start, win_end).skipped
            {
                trace!(
                    "check_complex: structurally anomalous read has a splice N inside \
                     the context window [{}, {}) at {}:{} — unscoreable across the \
                     splice, refusing string comparison → neither",
                    win_start, win_end, variant.chrom, variant.pos + 1
                );
                return ClassifyResult::neither(ClassifyPhase::Alignment);
            }
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
                                let mut r = classify_by_alignment(
                                    &sub_seq, &sub_quals, variant, min_baseq,
                                    alt_aligner, ref_aligner,
                                );
                                r.sw_fallback = true;
                                r
                            })
                        }
                    };
                }
            }
        }
    }

    // --- Phase 1: Haplotype Reconstruction ---
    // Walk the CIGAR to reconstruct what the read shows for [start_pos,
    // end_pos) (reconstruct_span — shared with the engine's multi-allelic
    // AD-claiming contest).
    let recon = reconstruct_span(record, quals, start_pos, end_pos);
    if recon.splice_skip {
        // The read asserts splicing over part of [start_pos, end_pos): the
        // reconstruction joined bases across the N gap, so comparing it
        // against the alleles would read exon-stitching as deletion evidence.
        // No clean string comparison exists for such a read — classify
        // neither. (Reads whose N covers EVERY discriminating position never
        // get here: splice_skip_triage already excluded them from coverage.)
        trace!(
            "check_complex: splice N overlaps variant span [{}, {}) at {}:{} — \
             reconstruction crosses the gap, no string-comparison evidence → neither",
            start_pos,
            end_pos,
            variant.chrom,
            variant.pos + 1
        );
        return ClassifyResult::neither(ClassifyPhase::CigarRecon);
    }
    let reconstructed_seq = recon.seq;
    let quals_per_base = recon.quals;

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
    let read_len = record.seq_len();
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
                            let mut r = classify_by_alignment(
                                &sub_seq, &sub_quals, variant, min_baseq,
                                alt_aligner, ref_aligner,
                            );
                            r.sw_fallback = true;
                            r
                        })
                    }
                };
            }
        }
    }

    // A read reaching here without a reference context was never evaluated by
    // the requested scorer. For a length-changing variant (indel or delins —
    // the variants prep fetches a context for) that means prep's context fetch
    // failed (it logs once), so no haplotype matrix — and no SW either — can be
    // built. Under PairHMM this is the same malformed-input condition the SW
    // fallback covers; mark it so the row carries SW_FALLBACK instead of
    // silently losing the read. MNPs carry no context by design (no Phase 3),
    // so their structural reads ending here are expected, not a fallback.
    let mut r = ClassifyResult::neither(ClassifyPhase::Alignment);
    if variant.ref_context.is_none()
        && variant.ref_allele.len() != variant.alt_allele.len()
        && matches!(backend, AlignmentBackend::PairHMM { .. })
    {
        r.sw_fallback = true;
    }
    r
}


/// Outcome of resolving an I op found at the exact expected junction
/// (`anchor_pos + 1`).
enum AnchorInsertionOutcome {
    /// Definitive resolution — the walk returns it immediately.
    Classified(ClassifyResult),
    /// Same-length I whose inserted bases cannot be verified (every base
    /// below `min_baseq`, or the op runs past the read end): the caller
    /// flags `has_shifted_same_length` and the walk continues to post-walk
    /// Phase-3 arbitration.
    Unverifiable,
}

/// Resolve a CIGAR I op sitting at the exact expected insertion junction
/// (`anchor_pos + 1`) against a pure-insertion variant.
///
/// Shared by the M-arm strict fast path and the post-N inspection — an op
/// directly after a splice N is a first-class candidate at the same genomic
/// position, and before the post-N inspection existed it was structurally
/// invisible. `ins_start` is the read index of the first inserted base;
/// `qual` is the fragment-consensus quality the caller attributes to this
/// evidence (the anchor base when it is aligned, or the first inserted base
/// when the anchor itself is spliced out).
fn resolve_anchor_insertion_candidate(
    record: &Record,
    variant: &Variant,
    found_ins_len: usize,
    ins_start: usize,
    qual: u8,
    quals: &[u8],
    min_baseq: u8,
) -> AnchorInsertionOutcome {
    let anchor_pos = variant.pos;
    let expected_ins_len = variant.alt_allele.len() - 1; // VCF ALT includes anchor
    let expected_ins_seq = &variant.alt_allele.as_bytes()[1..];

    if found_ins_len == expected_ins_len {
        if ins_start + found_ins_len <= record.seq().len() {
            let ins_seq = &record.seq().as_bytes()[ins_start..ins_start + found_ins_len];
            // Quality-aware fuzzy match for the inserted bases
            let ins_quals = &quals[ins_start..ins_start + found_ins_len];
            let (mismatches, reliable) =
                masked_single_compare(ins_seq, ins_quals, expected_ins_seq, min_baseq);
            if reliable > 0 && mismatches == 0 {
                trace!(
                    "check_insertion: I({}) match at expected junction after pos {}, \
                     qual={} (structural)",
                    found_ins_len, anchor_pos, qual
                );
                return AnchorInsertionOutcome::Classified(ClassifyResult::is_alt_structural(
                    qual,
                    ClassifyPhase::Structural,
                ));
            }
            if reliable > 0 {
                // Confident mismatch on ≥min_baseq bases: the read carries a
                // same-length insertion of DIFFERENT bases at the exact
                // junction — a third allele. Like the wrong-length rule: not
                // REF (rd must not absorb it), not the queried ALT, and
                // Phase 3 must not arbitrate — alignment scoring promotes a
                // wrong-sequence insert to ALT because it still beats the
                // gapped REF alignment. Partial evidence.
                trace!(
                    "check_insertion: I({}) at junction after {} matches expected \
                     length but bases mismatch (mismatches={}, reliable={}) — \
                     distinct allele → neither + partial evidence",
                    found_ins_len, anchor_pos, mismatches, reliable
                );
                return AnchorInsertionOutcome::Classified(ClassifyResult::neither_with_nearby(
                    qual,
                    ClassifyPhase::Structural,
                ));
            }
            // Every inserted base is below min_baseq: no confident read of
            // the inserted sequence, and a definitive call on
            // quality-rejected bases would break the cross-backend quality
            // contract. Post-walk Phase-3 arbitration (BQ-aware; windowed
            // matches still win first) propagates partial evidence on
            // non-ALT.
            trace!(
                "check_insertion: I({}) at junction after {} matches expected length \
                 but all inserted bases are below min_baseq → Phase-3 arbitration",
                found_ins_len, anchor_pos
            );
            return AnchorInsertionOutcome::Unverifiable;
        }
        // Insert runs past the read end (truncated record): the inserted
        // bases cannot be verified — same Phase-3 arbitration as the
        // low-quality case.
        trace!(
            "check_insertion: I({}) at junction after {} extends past the read end \
             → Phase-3 arbitration",
            found_ins_len, anchor_pos
        );
        return AnchorInsertionOutcome::Unverifiable;
    }

    // Wrong-length insertion at the exact junction: I(n), n ≠ expected.
    // Insertions are point events in reference space, so this is one of two
    // things, resolved in order below: a truncated copy of the expected
    // insert (same event), or a distinct allele in the same tract (+A vs
    // +AA slippage ladder) — partial evidence, never REF or ALT.

    // Truncation containment: sequencing loses bases from long inserted
    // sequences, so reads carry shorter I ops whose bases match a slice of
    // the expected insert. Same event → ALT (structural, like the
    // exact-length match).
    if ins_start + found_ins_len <= record.seq().len() {
        let ins_seq = &record.seq().as_bytes()[ins_start..ins_start + found_ins_len];
        if insert_truncation_match(ins_seq, expected_ins_seq) {
            trace!(
                "check_insertion: I({}) at junction after {} is a truncation of \
                 expected I({}) (≥90% substring identity) (structural)",
                found_ins_len, anchor_pos, expected_ins_len
            );
            return AnchorInsertionOutcome::Classified(ClassifyResult::is_alt_structural(
                qual,
                ClassifyPhase::Structural,
            ));
        }
    }

    // A wrong-length I that is not a truncation is a DIFFERENT allele: not
    // REF — rd must not absorb it — and not the queried ALT. Phase 3 must
    // not arbitrate: its haplotype window is length-blind inside repeat
    // tracts (and, with the default narrow context padding, can promote
    // not-the-event reads to definitive calls). Split representations of
    // the expected insert reach ALT through the containment check above —
    // a real split's anchor piece is a prefix of the expected insert.
    trace!(
        "check_insertion: I({}) at junction after {} vs expected I({}) — \
         wrong-length insertion → neither + partial evidence",
        found_ins_len, anchor_pos, expected_ins_len
    );
    AnchorInsertionOutcome::Classified(ClassifyResult::neither_with_nearby(
        qual,
        ClassifyPhase::Structural,
    ))
}

/// Evaluate a windowed (non-junction) I candidate at `ins_ref_pos`.
///
/// Shared by the M-arm windowed scan and the post-N inspection; mutates the
/// walk's resolution state. Candidates outside `[window_start, window_end]`
/// or at the exact expected junction (owned by
/// `resolve_anchor_insertion_candidate`) are ignored.
#[allow(clippy::too_many_arguments)]
fn scan_windowed_insertion_candidate(
    record: &Record,
    variant: &Variant,
    ins_ref_pos: i64,
    ins_len_usize: usize,
    ins_start: usize,
    quals: &[u8],
    min_baseq: u8,
    window_start: i64,
    window_end: i64,
    best_windowed_match: &mut Option<u64>,
    has_shifted_same_length: &mut bool,
    has_wrong_length_nearby: &mut bool,
) {
    let anchor_pos = variant.pos;
    if ins_ref_pos < window_start || ins_ref_pos > window_end || ins_ref_pos == anchor_pos + 1 {
        return;
    }
    let expected_ins_len = variant.alt_allele.len() - 1;
    let expected_ins_seq = &variant.alt_allele.as_bytes()[1..];
    let original_anchor_base = variant.ref_allele.as_bytes()[0].to_ascii_uppercase();

    // Safeguard 1: length must match
    if ins_len_usize == expected_ins_len {
        if ins_start + ins_len_usize <= record.seq().len() {
            let ins_seq = &record.seq().as_bytes()[ins_start..ins_start + ins_len_usize];
            // Quality-aware fuzzy match for inserted bases
            let ins_quals = &quals[ins_start..ins_start + ins_len_usize];
            let (mismatches, reliable) =
                masked_single_compare(ins_seq, ins_quals, expected_ins_seq, min_baseq);
            if reliable > 0 && mismatches == 0 {
                // Safeguard 3: verify anchor base at shifted position
                let shifted_anchor_pos = ins_ref_pos - 1;
                let anchor_ok = match &variant.ref_context {
                    Some(ctx) => {
                        let ctx_offset =
                            (shifted_anchor_pos - variant.ref_context_start) as usize;
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
                        warn!(
                            "ref_context is None for variant at {}:{} — S3 cannot validate shifted insertion",
                            variant.chrom, variant.pos + 1
                        );
                        false
                    }
                };

                if anchor_ok {
                    // Safeguard 2: track closest match
                    let distance = (ins_ref_pos - (anchor_pos + 1)).unsigned_abs();
                    if best_windowed_match.is_none_or(|prev| distance < prev) {
                        *best_windowed_match = Some(distance);
                    }
                } else {
                    trace!(
                        "check_insertion: S3 reject at shifted pos {} (anchor base mismatch)",
                        shifted_anchor_pos
                    );
                }
            } else {
                // Length matches but sequence differs — the caller and
                // aligner may represent the same event differently (e.g.,
                // shifted insertion in a repeat). Track this so Phase 3 can
                // arbitrate.
                *has_shifted_same_length = true;
                trace!(
                    "check_insertion: windowed I({}) at pos {} seq mismatch \
                     (mismatches={}, reliable={}), flagging for Phase 3 fallback",
                    ins_len_usize, ins_ref_pos, mismatches, reliable
                );
            }
        }
    } else {
        // Different-length insertion in window: I(n) where n ≠ expected — a
        // distinct-allele candidate, resolved by the post-walk handler.
        //
        // This also covers backward boundary insertions (anchor_pos ==
        // ref_pos): the same insertion is at block_end of the previous M
        // block, which the windowed scan processes on the prior loop
        // iteration.
        *has_wrong_length_nearby = true;
        trace!(
            "check_insertion: windowed I({}) at pos {} (expected I({})), \
             wrong length → distinct-allele candidate",
            ins_len_usize, ins_ref_pos, expected_ins_len
        );
    }
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
/// 4. **Wrong-length rule:** An I op at the anchor with a different length is
///    either a truncation of the expected insert (≥90% substring identity,
///    both sequences non-low-complexity → ALT; this also admits a real split
///    representation's anchor-prefix piece) or a distinct allele in the same
///    tract (neither + partial evidence, never Phase 3). Applies to every
///    single-base-REF variant the dispatcher sends here, including
///    anchor-substituting Ins+SNV (A>CCC).
/// 5. **Phase 3 haplotype fallback:** When a length-matching insertion exists
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

    // Windowed scan parameters — scales with repeat_span for MSI regions
    let window: i64 = std::cmp::max(5, variant.repeat_span as i64 + 2);
    let window_start = (anchor_pos - window).max(0);
    let window_end = anchor_pos + window;

    // State tracked across the CIGAR walk
    let mut found_ref_coverage = false;
    let mut anchor_read_pos: Option<usize> = None; // read position of anchor base
    let mut best_windowed_match: Option<u64> = None; // distance of best windowed match
    // Windowed I with the RIGHT length whose bases failed the sequence check:
    // the caller and aligner may represent the same event differently (e.g.
    // a shifted insertion in a repeat) — Phase 3 arbitrates after the walk.
    let mut has_shifted_same_length = false;
    // Windowed I with the WRONG length: a distinct-allele candidate,
    // resolved after the walk. Flagged at ANY length — unlike the deletion
    // side's ≥5bp noise gate — because wrong-length insertions must never be
    // silently absorbed into REF (the engine-level windowed-I contract).
    let mut has_wrong_length_nearby = false;

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
                        // Anchor is the last base of this match block: an I op
                        // directly after it sits at the exact expected junction.
                        if let Some(Cigar::Ins(ins_len)) = cigar_view.get(i + 1) {
                            let ins_start = read_pos + *len as usize;
                            let arp = anchor_read_pos.unwrap_or(0);
                            let qual = if arp < quals.len() { quals[arp] } else { 0 };
                            match resolve_anchor_insertion_candidate(
                                record, variant, *ins_len as usize, ins_start, qual,
                                quals, min_baseq,
                            ) {
                                AnchorInsertionOutcome::Classified(result) => return result,
                                AnchorInsertionOutcome::Unverifiable => {
                                    has_shifted_same_length = true;
                                }
                            }
                        }
                        // Anchor at the block end: either no I op followed (plain
                        // REF coverage) or an unverified same-length I was flagged
                        // above — found_ref_coverage still gates the post-walk
                        // Phase-3 branch for that flag.
                        found_ref_coverage = true;
                    } else {
                        // Anchor in middle of match block → read covers anchor without insertion
                        found_ref_coverage = true;
                    }
                }

                // --- Windowed scan: check if any Ins after this block is within window ---
                if let Some(Cigar::Ins(ins_len)) = cigar_view.get(i + 1) {
                    scan_windowed_insertion_candidate(
                        record, variant,
                        block_end, // genomic position where the insertion occurs
                        *ins_len as usize,
                        read_pos + *len as usize,
                        quals, min_baseq, window_start, window_end,
                        &mut best_windowed_match,
                        &mut has_shifted_same_length,
                        &mut has_wrong_length_nearby,
                    );
                }

                ref_pos = block_end;
                read_pos += *len as usize;
            }
            Cigar::Ins(len) => {
                read_pos += *len as usize;
            }
            Cigar::Del(len) => {
                ref_pos += *len as i64;
            }
            Cigar::RefSkip(len) => {
                let n_end = ref_pos + *len as i64;
                // --- Post-N inspection: an I op directly after a splice N sits
                // at genomic position n_end and gets the same inspection every
                // M arm gives its following op. Without this, an insertion at
                // an exon boundary reached only through the splice (M-N-I-M)
                // is structurally invisible: no M arm precedes the I, so
                // carriers fell through to Phase 3 or REF.
                if let Some(Cigar::Ins(ins_len)) = cigar_view.get(i + 1) {
                    let ins_len_usize = *ins_len as usize;
                    if n_end == anchor_pos + 1 {
                        // I at the exact expected junction with the anchor base
                        // spliced out: attribute the evidence to the first
                        // inserted base (there is no anchor base to read).
                        let qual = if read_pos < quals.len() { quals[read_pos] } else { 0 };
                        match resolve_anchor_insertion_candidate(
                            record, variant, ins_len_usize, read_pos, qual,
                            quals, min_baseq,
                        ) {
                            AnchorInsertionOutcome::Classified(result) => return result,
                            AnchorInsertionOutcome::Unverifiable => {
                                // The M-arm twin flags has_shifted_same_length
                                // and lets post-walk Phase 3 arbitrate with
                                // partial propagation — but Phase 3 cannot
                                // score across this read's splice (extraction
                                // refuses N-crossing windows), so that route
                                // ends in a bare neither and silently drops
                                // the structural evidence. A same-length I at
                                // the exact junction whose bases cannot be
                                // verified is partial evidence: say so
                                // directly.
                                return ClassifyResult::neither_with_nearby(
                                    qual,
                                    ClassifyPhase::Structural,
                                );
                            }
                        }
                    } else {
                        scan_windowed_insertion_candidate(
                            record, variant, n_end, ins_len_usize,
                            read_pos, // N consumes no read bases
                            quals, min_baseq, window_start, window_end,
                            &mut best_windowed_match,
                            &mut has_shifted_same_length,
                            &mut has_wrong_length_nearby,
                        );
                    }
                }
                ref_pos = n_end;
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

    // Phase 3 haplotype fallback: a same-length insertion exists nearby but its
    // sequence failed the check — the caller and aligner may represent the same
    // event differently (shifted insertion in a repeat). Phase 3
    // (Smith-Waterman/PairHMM) arbitrates; propagate has_nearby_evidence if it
    // doesn't confirm ALT, so the engine counts the read as partial_alt.
    if has_shifted_same_length && found_ref_coverage {
        trace!(
            "check_insertion: nearby same-length insertion (seq mismatch) at pos {}, \
             falling back to phase3_classify",
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

    // Wrong-length insertion in the window (no same-length candidate). Note:
    // NOT the strict-path resolution — the truncation-containment check does
    // not run here (a windowed I's bases are not extracted during the scan);
    // shifted truncations of a long insert therefore land in partial, not ALT.
    // A wrong-length I inside a repeat tract is a distinct slippage allele →
    // neither + partial evidence (Phase 3 must not arbitrate: it is
    // length-blind inside repeat tracts). In unique context the
    // anchor-covering M is definitive REF and the stray I is alignment
    // noise — keep rd, but carry the partial evidence so the read is not
    // silently absorbed (matches the pre-rule Phase-3 outcome).
    if has_wrong_length_nearby && found_ref_coverage {
        if variant.repeat_span >= 2 {
            trace!(
                "check_insertion: lone wrong-length I in repeat tract at pos {} → \
                 neither + partial evidence",
                anchor_pos
            );
            return ClassifyResult::neither_with_nearby(anchor_qual, ClassifyPhase::Structural);
        }
        trace!(
            "check_insertion: lone wrong-length I in window at pos {} (unique \
             context) → REF at anchor + partial evidence",
            anchor_pos
        );
        let mut result = ClassifyResult::is_ref(anchor_qual, ClassifyPhase::Structural);
        result.has_nearby_evidence = true;
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
/// Used by `check_deletion`'s Safeguard 3 for windowed matches (exact-length
/// and in-band) and as the context-integrity guard of the strict-path band.
/// For a SHIFTED position, comparing the reference at the observed breakpoint
/// against `expected_del_seq` rejects a deletion of different bases; at the
/// anchor itself both operands are reference-derived at the same coordinates,
/// so the comparison can only fail on a missing/short ref_context. Returns
/// `true` only on a reliable, zero-mismatch concordance; a `None` ref_context
/// or out-of-bounds offset returns `false`.
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

/// CIGAR indel summary near a variant, over I/D ops whose reference placement
/// overlaps `[region_start, region_end]` (0-based, inclusive; insertions are
/// point events at the boundary they precede).
///
/// The span fields are placement-AWARE: `deleted_in_span` clips each D op to
/// the expected deleted interval `[span_start, span_end)`, so a read that
/// deletes elsewhere in the region cannot masquerade as the expected event
/// (a D(2)-at-anchor + M(53) + D(49) read nets ~50 deleted bases but deletes
/// only 2 of an expected 50-base span — the M(53) proves the event is absent).
struct RegionIndelSummary {
    /// Deleted bases falling INSIDE `[span_start, span_end)` (clipped overlap).
    deleted_in_span: i64,
    /// Deleted bases from region ops falling OUTSIDE the span, plus all
    /// inserted bases in the region — reference-length changes the expected
    /// deletion does not explain. Small values (breakpoint wobble, split-gap
    /// slop) are tolerable; large values mean a different event.
    changed_outside_span: i64,
}

/// Walk the read's CIGAR once and summarize its indel ops near the variant.
/// Feeds `check_deletion`'s placement-aware large-deletion band.
fn summarize_indels_in_region(
    record: &Record,
    region_start: i64,
    region_end: i64,
    span_start: i64,
    span_end: i64,
) -> RegionIndelSummary {
    let mut ref_pos = record.pos();
    let mut summary = RegionIndelSummary {
        deleted_in_span: 0,
        changed_outside_span: 0,
    };
    for op in record.cigar().iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) | Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            }
            Cigar::Del(len) => {
                let del_len = *len as i64;
                if ref_pos <= region_end && ref_pos + del_len > region_start {
                    let in_span =
                        (ref_pos + del_len).min(span_end) - ref_pos.max(span_start);
                    let in_span = in_span.max(0);
                    summary.deleted_in_span += in_span;
                    summary.changed_outside_span += del_len - in_span;
                }
                ref_pos += del_len;
            }
            Cigar::Ins(len) if ref_pos >= region_start && ref_pos <= region_end => {
                summary.changed_outside_span += *len as i64;
            }
            _ => {}
        }
    }
    summary
}

/// True when an observed inserted sequence reads as a *truncation* of the
/// expected insert: sequencing loses bases from long insertions, so reads
/// carry shorter I ops whose bases still match a contiguous slice of the
/// expected insert (a validated long-insertion cohort shows this smear —
/// e.g. I(11)/I(14) prefixes alongside the full-length I(19) population).
/// Sign-out counts those as the same event; so must we.
///
/// Four gates keep this from admitting distinct alleles:
/// - the observed insert must be ≥4bp and strictly shorter than expected —
///   1-3bp fragments match almost any sequence;
/// - neither the expected insert nor the OBSERVED insert may be
///   low-complexity (longest tandem repeat covering ≥80% of its length):
///   inside a repeat tract every wrong-length insert matches a substring
///   trivially, and a repeat-only observed insert matching a repeat
///   subregion of a complex expected insert is slippage, not truncation —
///   in both cases different lengths are distinct alleles;
/// - the best-matching window of the expected insert must be ≥90% identical
///   to the observed bases (tolerates isolated sequencing errors).
fn insert_truncation_match(observed: &[u8], expected: &[u8]) -> bool {
    const MIN_OBSERVED_LEN: usize = 4;
    if observed.len() < MIN_OBSERVED_LEN || observed.len() >= expected.len() {
        return false;
    }
    let low_complexity = |seq: &[u8]| {
        let max_span = (0..seq.len())
            .map(|i| find_tandem_repeat(seq, i).1)
            .max()
            .unwrap_or(1);
        max_span as f64 >= 0.8 * seq.len() as f64
    };
    if low_complexity(expected) || low_complexity(observed) {
        return false;
    }
    let mut best_matches = 0usize;
    for start in 0..=(expected.len() - observed.len()) {
        let matches = observed
            .iter()
            .zip(&expected[start..start + observed.len()])
            .filter(|(o, e)| o.eq_ignore_ascii_case(e))
            .count();
        best_matches = best_matches.max(matches);
    }
    best_matches as f64 >= 0.9 * observed.len() as f64
}

/// Resolve a CIGAR D op that starts exactly at the expected deletion start
/// (`anchor_pos + 1`) against a pure-deletion variant.
///
/// Shared by the M-arm strict fast path and the post-N inspection — an op
/// directly after a splice N is a first-class candidate at the same genomic
/// position, and before the post-N inspection existed it was structurally
/// invisible. Every branch is definitive: exact length → ALT; in the
/// large-deletion band with verified context → ALT; otherwise a distinct
/// allele in the same tract → neither + partial evidence. `qual` is the
/// fragment-consensus quality the caller attributes to this evidence (the
/// anchor base when it is aligned, or the first aligned base after the
/// deletion when the anchor itself is spliced out).
fn resolve_anchor_deletion_candidate(
    record: &Record,
    variant: &Variant,
    found_del_len: usize,
    qual: u8,
    window_start: i64,
    window: i64,
) -> ClassifyResult {
    let anchor_pos = variant.pos;
    let expected_del_len = variant.ref_allele.len() - 1; // REF without anchor
    let expected_del_seq = &variant.ref_allele.as_bytes()[1..];

    if found_del_len == expected_del_len {
        trace!(
            "check_deletion: D({}) match at expected span after pos {}, qual={} (structural)",
            found_del_len, anchor_pos, qual
        );
        return ClassifyResult::is_alt_structural(qual, ClassifyPhase::Structural);
    }

    // D found at the expected start but with the WRONG length. For a PURE
    // deletion (the dispatcher guarantees purity — delins go to
    // check_complex) this is one of two things, resolved in order below:
    // the same large event with breakpoint wobble or a split representation
    // (net band), or a distinct allele in the same tract — partial
    // evidence, never REF or ALT.
    let span_start = anchor_pos + 1;
    let span_end = span_start + expected_del_len as i64;
    let region_end = anchor_pos + expected_del_len as i64 + window;
    let indels = summarize_indels_in_region(record, window_start, region_end, span_start, span_end);

    // Large-deletion band: validated pure large deletions align as a single
    // exact-length D at the anchor (zero wobble across 12-539bp events), so
    // the same event allows at most minimal slop. Placement-aware, both
    // ways: the read must DELETE essentially the whole expected span (≤3
    // retained bases — covers split representations like D(60)+2M+D(40)),
    // and must not delete/insert more than 3bp beyond it (breakpoint
    // wobble). A read that nets the right total while matching reference
    // across the span (D(2) at anchor + M(53) + D(49) downstream) fails the
    // first test — its M ops prove the event is absent. This replaces the
    // earlier SV-caller-style ≥50% reciprocal-overlap rule
    // (SURVIVOR/BEDTools precedent): at a SHARED anchor that admitted any
    // deletion sharing half the length — a D(60) passed for a 100bp variant.
    if expected_del_len >= 50
        && expected_del_len as i64 - indels.deleted_in_span <= 3
        && indels.changed_outside_span <= 3
    {
        // Context-integrity guard only: both operands are reference-derived
        // at the same coordinates, so this can fail solely on a
        // missing/short ref_context — never on sequence. It keeps a variant
        // with broken context out of the structural-ALT fast path.
        let compare_len = found_del_len.min(expected_del_len);
        if verify_deleted_bases(variant, anchor_pos + 1, expected_del_seq, compare_len) {
            trace!(
                "check_deletion: band match at pos {} — deleted {}/{} span bases, \
                 {}bp outside (D({}) at expected start, qual={}) (structural)",
                anchor_pos,
                indels.deleted_in_span,
                expected_del_len,
                indels.changed_outside_span,
                found_del_len,
                qual,
            );
            return ClassifyResult::is_alt_structural(qual, ClassifyPhase::Structural);
        }
        trace!(
            "check_deletion: in-band D({}) at expected start after {} REJECTED — \
             ref_context missing or too short for validation",
            found_del_len,
            anchor_pos,
        );
    }

    // Outside the band, a wrong-length D at the expected start is definitive
    // evidence of a DIFFERENT allele in the same tract: coexisting
    // distinct-length populations are distinct slippage alleles (D1 vs D2 in
    // a homopolymer, the GGC-repeat ladder), and the placement-aware band
    // above already recognized every same-event composition (lone wobble and
    // split representations alike). Not REF, not the queried ALT — partial
    // evidence. Phase 3 must not arbitrate: its haplotype window is
    // length-blind inside repeat tracts (and, with narrow context padding,
    // can promote not-the-event reads to definitive calls).
    trace!(
        "check_deletion: D({}) at expected start after {} vs expected D({}) — \
         wrong-length pure deletion → neither + partial evidence",
        found_del_len, anchor_pos, expected_del_len
    );
    ClassifyResult::neither_with_nearby(qual, ClassifyPhase::Structural)
}

/// Evaluate a windowed (non-anchor) D candidate at `del_ref_pos`.
///
/// Shared by the M-arm windowed scan and the post-N inspection; mutates the
/// walk's resolution state. Candidates outside `[window_start, window_end]`
/// or at the exact expected start (owned by
/// `resolve_anchor_deletion_candidate`) are ignored.
#[allow(clippy::too_many_arguments)]
fn scan_windowed_deletion_candidate(
    variant: &Variant,
    del_ref_pos: i64,
    del_len_usize: usize,
    window_start: i64,
    window_end: i64,
    best_windowed_match: &mut Option<u64>,
    has_shifted_same_length: &mut bool,
    has_wrong_length_nearby: &mut bool,
) {
    let anchor_pos = variant.pos;
    if del_ref_pos < window_start || del_ref_pos > window_end || del_ref_pos == anchor_pos + 1 {
        return;
    }
    let expected_del_len = variant.ref_allele.len() - 1;
    let expected_del_seq = &variant.ref_allele.as_bytes()[1..];

    // Safeguard 1: deletion length check.
    // Accept exact matches, OR for large deletions (≥50bp) a ≤3bp length
    // difference — the same breakpoint band as the anchor-candidate path
    // (validated pure large deletions show zero length wobble, so anything
    // beyond minimal wobble is a different event, not the same one
    // re-reported).
    let len_diff = expected_del_len.abs_diff(del_len_usize);
    let length_ok = if del_len_usize == expected_del_len {
        true
    } else if expected_del_len >= 50 && len_diff <= 3 {
        trace!(
            "check_deletion: windowed band match D({}) ≈ D({}) at pos {} (|Δlen|={})",
            del_len_usize, expected_del_len, del_ref_pos, len_diff
        );
        true
    } else if del_len_usize >= 5 {
        // Wrong-length D in the window: a distinct-allele candidate,
        // resolved after the walk (partial, or Phase 3 when the read looks
        // like a split representation).
        *has_wrong_length_nearby = true;
        trace!(
            "check_deletion: windowed D({}) at pos {} (expected D({})), \
             wrong length → distinct-allele candidate",
            del_len_usize, del_ref_pos, expected_del_len
        );
        false
    } else {
        // Short (1-4bp) wrong-length deletions in the window are spurious
        // alignment noise in homopolymer/STR regions; CIGAR is definitive
        // for those.
        trace!(
            "check_deletion: windowed D({}) at pos {} (expected D({})), \
             different-length (<5bp) → CIGAR definitive, not flagging",
            del_len_usize, del_ref_pos, expected_del_len
        );
        false
    };

    if length_ok {
        // Safeguard 3: verify the deleted reference bases match
        // expected_del_seq. Exact-length matches compare the full span;
        // in-band (different-length) matches compare the OVERLAPPING span
        // instead of skipping verification, so an unrelated deletion of
        // different bases is not accepted as ALT. A genuinely
        // shifted/different deletion fails here and is routed to the
        // Phase 3 fallback after the walk.
        let compare_len = del_len_usize.min(expected_del_len);
        let del_ok = verify_deleted_bases(variant, del_ref_pos, expected_del_seq, compare_len);

        if del_ok {
            // Safeguard 2: track closest match
            let distance = (del_ref_pos - (anchor_pos + 1)).unsigned_abs();
            if best_windowed_match.is_none_or(|prev| distance < prev) {
                *best_windowed_match = Some(distance);
            }
        } else {
            // S3 failed: deleted bases at the shifted position don't match
            // expected_del_seq. This happens when left-alignment moves the
            // anchor further left than the aligner's CIGAR Del position, so
            // the reference slice at del_ref_pos encodes different bases
            // than expected_del_seq (e.g. TP53 GACCGTGCAAGT→- where
            // left-alignment shifts anchor 3bp left of the actual D(12) in
            // reads). Track for Phase 3 fallback.
            //
            // Only flag for Phase 3 when del_len >= 5bp. Short (1-4bp)
            // same-length deletions that fail S3 are almost certainly
            // spurious/unrelated noise — CIGAR remains definitive for them.
            // Longer deletions are more susceptible to BWA left-alignment
            // shifting the anchor multiple positions away from the CIGAR D.
            if del_len_usize >= 5 {
                *has_shifted_same_length = true;
                trace!(
                    "check_deletion: S3 reject at shifted pos {} (deleted bases \
                     mismatch, del_len={} >= 5), flagging for Phase 3 fallback",
                    del_ref_pos, del_len_usize
                );
            } else {
                trace!(
                    "check_deletion: S3 reject at shifted pos {} (deleted bases \
                     mismatch, del_len={} < 5, CIGAR definitive — not flagging Phase 3)",
                    del_ref_pos, del_len_usize
                );
            }
        }
    }
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
/// 3. **Wrong-length rule (pure deletions):** A D op at the anchor with a
///    different length is the same large event only when the read deletes
///    essentially the whole expected span (≤3 retained bases, ≤3 changed
///    outside — ≥50bp deletions; covers breakpoint wobble and split
///    representations → ALT); otherwise it is a distinct allele in the same
///    tract (neither + partial evidence, never Phase 3).
/// 4. **Haplotype fallback:** When CIGAR geometry doesn't match (e.g. a
///    soft-clip at the anchor, or a shifted same-length D whose bases fail
///    S3), delegates to `check_complex` for quality-aware haplotype
///    comparison.
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

    // Left-anchored invariant (LO-8): `variant.pos` MUST be the retained
    // reference anchor base immediately *before* the deletion — never a
    // coordinate inside the deleted span. The candidate helpers take the
    // deleted bases as `ref_allele[1..]` (VCF/MAF deletions carry the anchor as
    // ref_allele[0]); prep-time left-alignment/normalization guarantees this
    // shape. An empty REF/ALT (malformed or non-left-anchored record) would
    // underflow their `len() - 1` and panic the `[1..]` slice. Defend it and
    // classify as neither. debug! not warn!: runs per-read; loud
    // once-per-variant surfacing belongs at prep (follow-up).
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
    // Deleted-span bounds for span-aligned REF testimony (0-based, half-open)
    let span_start = anchor_pos + 1;
    let span_end = anchor_pos + variant.ref_allele.len() as i64;

    // Windowed scan parameters — scales with repeat_span for MSI regions
    let window: i64 = std::cmp::max(5, variant.repeat_span as i64 + 2);
    let window_start = (anchor_pos - window).max(0);
    let window_end = anchor_pos + window;

    let mut found_ref_coverage = false;
    let mut anchor_read_pos: Option<usize> = None; // read position of anchor base
    let mut best_windowed_match: Option<u64> = None;
    // Windowed Del with the RIGHT length (or within the large-deletion band)
    // whose bases at the shifted position failed the S3 sequence check: BWA
    // left-alignment can move the anchor several bases away from the CIGAR D
    // in repeat context — Phase 3 haplotype comparison arbitrates after the walk.
    let mut has_shifted_same_length = false;
    // Windowed Del with the WRONG length (≥5bp, outside the band): a
    // distinct-allele candidate, resolved after the walk (partial, or
    // Phase 3 for split representations).
    let mut has_wrong_length_nearby = false;
    // Span-aligned REF testimony state: whether the anchor base sits inside
    // a splice N of THIS read, how many deleted-span positions its M ops
    // cover, and the read index of the first covered span base (quality
    // attribution when there is no anchor base to read).
    let mut anchor_in_refskip = false;
    let mut span_aligned: i64 = 0;
    let mut span_first_read_pos: Option<usize> = None;

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

                // Track aligned coverage of the deleted span (for span-aligned
                // REF testimony when the anchor itself is spliced out).
                let ov_start = ref_pos.max(span_start);
                let ov_end = block_end.min(span_end);
                if ov_start < ov_end {
                    span_aligned += ov_end - ov_start;
                    if span_first_read_pos.is_none() {
                        span_first_read_pos =
                            Some(read_pos + (ov_start - ref_pos) as usize);
                    }
                }

                // --- Strict fast path ---
                if anchor_pos >= ref_pos && anchor_pos < block_end {
                    if anchor_pos == block_end - 1 {
                        // Anchor at end of match: a D op directly after it
                        // starts at the exact expected deletion position.
                        if let Some(Cigar::Del(del_len)) = cigar_view.get(i + 1) {
                            let arp = anchor_read_pos.unwrap_or(0);
                            let qual = if arp < quals.len() { quals[arp] } else { 0 };
                            // Every branch is definitive (ALT, band-ALT, or
                            // distinct allele → partial), so return directly.
                            return resolve_anchor_deletion_candidate(
                                record, variant, *del_len as usize, qual,
                                window_start, window,
                            );
                        }
                        found_ref_coverage = true;
                    } else {
                        // Anchor in middle of match → REF coverage
                        found_ref_coverage = true;
                    }
                }

                // --- Windowed scan: check for Del after this block within window ---
                if let Some(Cigar::Del(del_len)) = cigar_view.get(i + 1) {
                    scan_windowed_deletion_candidate(
                        variant,
                        block_end, // genomic position where the deletion starts
                        *del_len as usize,
                        window_start, window_end,
                        &mut best_windowed_match,
                        &mut has_shifted_same_length,
                        &mut has_wrong_length_nearby,
                    );
                }

                ref_pos = block_end;
                read_pos += *len as usize;
            }
            Cigar::Del(len) => {
                ref_pos += *len as i64;
            }
            Cigar::RefSkip(len) => {
                let n_end = ref_pos + *len as i64;
                // The anchor base inside this read's splice N: the read
                // cannot show it (span-aligned REF testimony may still apply
                // after the walk).
                if *len > 0 && anchor_pos >= ref_pos && anchor_pos < n_end {
                    anchor_in_refskip = true;
                }
                // --- Post-N inspection: a D op directly after a splice N
                // starts at genomic position n_end and gets the same
                // inspection every M arm gives its following op. Without
                // this, a deletion at an exon boundary reached only through
                // the splice (M-N-D-M) is structurally invisible: no M arm
                // precedes the D, so carriers fell through to Phase 3 or
                // neither.
                if let Some(Cigar::Del(del_len)) = cigar_view.get(i + 1) {
                    let del_len_usize = *del_len as usize;
                    if n_end == anchor_pos + 1 {
                        // D at the exact expected span with the anchor base
                        // spliced out: attribute the evidence to the first
                        // aligned base after the deletion (there is no anchor
                        // base to read).
                        let qual = if read_pos < quals.len() { quals[read_pos] } else { 0 };
                        return resolve_anchor_deletion_candidate(
                            record, variant, del_len_usize, qual, window_start, window,
                        );
                    }
                    scan_windowed_deletion_candidate(
                        variant, n_end, del_len_usize, window_start, window_end,
                        &mut best_windowed_match,
                        &mut has_shifted_same_length,
                        &mut has_wrong_length_nearby,
                    );
                }
                ref_pos = n_end;
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

    // Span-aligned REF testimony: the anchor base is inside this read's
    // splice N (asserted splicing — the read cannot show it), but every
    // deleted-span base is covered by aligned M ops: the annotated deletion
    // is demonstrably absent from this read. The span, not the anchor, is
    // the discriminating fact for a pure deletion, so this is REF — the
    // same coverage-based standard the anchor fast path applies (presence,
    // not per-base identity; base mismatches inside the span are SNV-level
    // signal, not deletion evidence). Quality is attributed to the first
    // aligned span base (there is no anchor base to read). Guards: the FULL
    // span must be aligned (partial coverage cannot rule the deletion out),
    // and any competing indel candidate on this read (shifted same-length,
    // wrong-length nearby) keeps its existing arbitration path. Typical
    // population: exon-anchored deletion annotations left-aligned to the
    // last intronic base — at a validated FORTE acceptor locus, 823 of the
    // 824 anchor-spliced junction reads convert (the residual covers only
    // part of the span and stays neither).
    if !found_ref_coverage
        && anchor_in_refskip
        && !has_shifted_same_length
        && !has_wrong_length_nearby
        && span_aligned == span_end - span_start
    {
        let qual = span_first_read_pos
            .filter(|&p| p < quals.len())
            .map(|p| quals[p])
            .unwrap_or(0);
        trace!(
            "check_deletion: anchor {} spliced out but deleted span [{}, {}) fully \
             aligned — span-aligned REF testimony, qual={}",
            anchor_pos, span_start, span_end, qual
        );
        // Structural: the coverage is the evidence. The carried qual is the
        // first span base — right after the junction, where BAQ can zero it
        // — and fragment consensus must not let a zero erase the
        // observation (has_structural_ref mirrors the ALT-side rule).
        return ClassifyResult::is_ref_structural(qual, ClassifyPhase::Structural);
    }

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

    // Phase 3 haplotype fallback: a same-length (or in-band) deletion exists
    // nearby but its shifted-position bases failed the sequence check — BWA
    // left-alignment can shift the anchor several bases away from the CIGAR D
    // in repeat context. Phase 3 rebuilds the haplotype and arbitrates;
    // propagate has_nearby_evidence if it doesn't confirm ALT.
    if has_shifted_same_length && found_ref_coverage {
        trace!(
            "check_deletion: nearby same-length deletion (seq mismatch) at pos {}, \
             falling back to phase3_classify",
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

    // Wrong-length deletion in the window (no same-length candidate; flagged
    // ops are ≥5bp — smaller ones were dropped as noise during the scan).
    // The read shows REF at the anchor junction, so the expected event is
    // absent; the wrong-length D nearby is a distinct slippage allele inside
    // a repeat tract → neither + partial evidence (Phase 3 must not
    // arbitrate: it is length-blind inside repeat tracts). In unique context
    // the anchor-covering M is definitive REF — keep rd, but carry the
    // partial evidence so the nearby deletion is not silently hidden
    // (matches the pre-rule Phase-3 outcome).
    if has_wrong_length_nearby && found_ref_coverage {
        if variant.repeat_span >= 2 {
            trace!(
                "check_deletion: lone wrong-length D in repeat tract at pos {} → \
                 neither + partial evidence",
                anchor_pos
            );
            return ClassifyResult::neither_with_nearby(anchor_qual, ClassifyPhase::Structural);
        }
        trace!(
            "check_deletion: lone wrong-length D in window at pos {} (unique \
             context) → REF at anchor + partial evidence",
            anchor_pos
        );
        let mut result = ClassifyResult::is_ref(anchor_qual, ClassifyPhase::Structural);
        result.has_nearby_evidence = true;
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

    // ── insert_truncation_match — wrong-length insertion containment rule ──
    // Expected insert must be complex; observed must be ≥4bp, shorter, and
    // ≥90% identical to some window of the expected insert.

    const COMPLEX_INS: &[u8] = b"CTTAGTCACCTTCGTGGCA"; // 19bp, no dominant repeat

    #[test]
    fn test_truncation_prefix_matches() {
        assert!(insert_truncation_match(&COMPLEX_INS[..11], COMPLEX_INS));
    }

    #[test]
    fn test_truncation_interior_slice_matches() {
        assert!(insert_truncation_match(&COMPLEX_INS[5..13], COMPLEX_INS));
    }

    #[test]
    fn test_truncation_tolerates_one_error_in_eleven() {
        // 10/11 identity = 0.909 ≥ 0.90
        let mut obs = COMPLEX_INS[..11].to_vec();
        obs[6] = b'T'; // was 'C'
        assert!(insert_truncation_match(&obs, COMPLEX_INS));
    }

    #[test]
    fn test_truncation_rejects_two_errors_in_eleven() {
        // 9/11 identity = 0.818 < 0.90
        let mut obs = COMPLEX_INS[..11].to_vec();
        obs[6] = b'T';
        obs[7] = b'G'; // was 'A'
        assert!(!insert_truncation_match(&obs, COMPLEX_INS));
    }

    #[test]
    fn test_truncation_rejects_unrelated_sequence() {
        assert!(!insert_truncation_match(b"GGGGGGGGGGG", COMPLEX_INS));
    }

    #[test]
    fn test_truncation_rejects_low_complexity_expected() {
        // In a slippage tract every wrong-length insert matches trivially —
        // those are distinct alleles, so the rule must not fire.
        assert!(!insert_truncation_match(b"AAAA", b"AAAAAAAAAA"));
        assert!(!insert_truncation_match(b"GGCGGC", b"GGCGGCGGCGGC"));
    }

    #[test]
    fn test_truncation_rejects_low_complexity_observed() {
        // A repeat-only observed insert matching a repeat SUBREGION of a
        // complex expected insert is slippage, not truncation: the expected
        // insert passes the whole-sequence gate (9bp poly-A of 14bp < 80%),
        // so the observed-side gate must reject it.
        assert!(!insert_truncation_match(b"AAAAAA", b"ACGTAAAAAAAAAC"));
        // A complex slice of the same expected insert still matches.
        assert!(insert_truncation_match(b"ACGTAA", b"ACGTAAAAAAAAAC"));
    }

    #[test]
    fn test_truncation_rejects_short_and_non_shorter_observed() {
        assert!(!insert_truncation_match(&COMPLEX_INS[..3], COMPLEX_INS)); // <4bp
        assert!(!insert_truncation_match(COMPLEX_INS, COMPLEX_INS)); // equal length
        let longer = [COMPLEX_INS, b"TT"].concat();
        assert!(!insert_truncation_match(&longer, COMPLEX_INS)); // longer than expected
    }
}
