//! Counting engine: BAM/CRAM read classification → allele counts.
//!
//! Entry point: `count_bam_binned()` groups variants into genomic bins,
//! issues one `bam.fetch()` per bin, and classifies each read against all
//! variants in the bin via type-specific dispatchers (check_snp, check_mnp,
//! check_complex, check_insertion, check_deletion).
//!
//! ## Invariants maintained during counting
//!
//! - `any_alt = ad + partial_alt` (decomposed ALT counting)
//! - `any_alt >= ad` (partial_alt is non-negative)
//! - `DP >= RD + AD + partial_alt + n_count` (depth decomposition)
//! - N-base reads increment `n_count` but do NOT contribute to RD, AD,
//!   any_alt, or partial_alt (uninformative signal).
//!
//! ## Output: `BaseCounts`
//!
//! Each variant produces a `BaseCounts` struct with:
//! - Core: dp, rd, ad (read-level) + dpf, rdf, adf (fragment-level)
//! - Strand: dp_fwd/rev, rd_fwd/rev, ad_fwd/rev, rdf_fwd/rev, adf_fwd/rev
//! - Bias: sb_pval, sb_or, fsb_pval, fsb_or
//! - Diagnostic: any_alt, partial_alt, n_count (Phase 2/2b)
//! - RNA: sense_depth, antisense_depth, sense_strand_alt_count, splice_spanning_count
//! - Annotation: exon_boundary_dist (GTF-informed, RNA mode only)

use pyo3::prelude::*;
#[cfg(test)]
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{self, Read, Record};
use std::collections::{HashMap, HashSet};

use crate::annotation::AnnotationIndex;
use crate::shared::stats::fisher_strand_bias;
use crate::types::{
    BaseCounts, Observation, Variant, OBS_ALLELE_ALT, OBS_ALLELE_N, OBS_ALLELE_OTHER,
    OBS_ALLELE_REF,
};

use rayon::prelude::*;

/// What one genomic bin produces: `(vi, counts)` pairs, the per-molecule rows for
/// those variants (empty unless observations were requested), and how many fetched
/// records the read filter dropped for having no bases (warned once per BAM).
/// Per bin: (variant index, counts) pairs, observations, and the identities of the
/// records dropped for having no bases / no base qualities. Overlapping bins fetch
/// a record more than once, so the identities are merged as sets: each record is
/// counted once in the per-BAM warning.
type BinOutput = (
    Vec<(usize, BaseCounts)>,
    Vec<Observation>,
    (HashSet<u64>, HashSet<u64>),
);

use anyhow::{Context, Result};
use log::{debug, info, trace, warn};
use bio::alignment::pairwise::Aligner;

use super::fragment::{FragmentEvidence, MoleculeClass, hash_qname, hash_molecule};
use super::alignment::{SW_GAP_EXTEND, SW_GAP_OPEN};
use super::variant_checks::{check_snp, check_mnp, check_complex, check_insertion, check_deletion, splice_skip_triage, reconstruct_span, MnpResult};
use bio::alignment::distance::levenshtein;
use super::utils::{find_read_pos, ref_end, soft_clips, ClassifyResult, ClassifyPhase};
use super::mfsd;
use super::rna;
use super::window::{self, AlleleKind};
use super::observed;
use super::carrier;
use crate::shared::baq::apply_heuristic_baq;


// BAQ (Base Alignment Quality) heuristic now lives in shared::baq.
// See crate::shared::baq::apply_heuristic_baq for implementation details.


/// Compute the median of a u32 vector. Returns 0.0 for empty input.
///
/// Uses `sort_unstable()` for optimal performance (no allocation,
/// not a stable sort — fine for numeric values).
fn compute_median_u32(v: &mut [u32]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_unstable();
    let mid = v.len() / 2;
    if v.len().is_multiple_of(2) {
        (v[mid - 1] as f64 + v[mid] as f64) / 2.0
    } else {
        v[mid] as f64
    }
}


// ── Genomic Binning ──────────────────────────────────────────────────────
// Groups co-located variants into ~10kb bins so the binned engine can do
// a single bam.fetch() per bin instead of per variant. This reduces I/O
// overhead significantly when variants are clustered (e.g., MAF files with
// thousands of variants on the same gene).
//
// A bin spans at least BIN_WINDOW from its first variant; each variant added
// extends its end to that variant's REF span plus half a window, so close variants
// chain it further (a floor, never capped: every member's full read window must be
// fetched). It closes at a variant past its end or at BIN_MAX_VARIANTS. The C++
// GBCMS used --max_block_size=200 and --max_block_dist=10000.

/// Default bin window size in base pairs.
const BIN_WINDOW: i64 = 10_000;

/// Maximum number of variants per bin.
///
/// When exceeded, the bin is split to prevent O(V × R) blowup in the
/// shared-read classification loop. Matches the original C++ GBCMS default.
/// Not a CLI flag: an internal performance constant that does not affect output
/// (the binning-invariance tests vary it, with the window, through the
/// `bin_window` / `bin_max_variants` test arguments).
const BIN_MAX_VARIANTS: usize = 200;

/// The bin geometry for one call: the production constants unless a test passes
/// its own. Bin geometry is performance only, so tests vary it to check that
/// counts do not depend on it; a window or cap below 1 is rejected.
fn bin_geometry(bin_window: Option<i64>, bin_max_variants: Option<i64>) -> PyResult<(i64, usize)> {
    let window = bin_window.unwrap_or(BIN_WINDOW);
    let cap = bin_max_variants.unwrap_or(BIN_MAX_VARIANTS as i64);
    if window < 1 || cap < 1 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "bin_window and bin_max_variants must be at least 1 (got {window}, {cap})"
        )));
    }
    Ok((window, cap as usize))
}

/// A genomic region containing one or more co-located variants.
///
/// The engine fetches reads once for the entire bin, then classifies
/// each read against all variants in the bin. Padding ensures that
/// reads overlapping bin boundaries are not missed.
#[derive(Debug)]
struct GenomicBin {
    /// BAM target ID (tid) for this chromosome.
    tid: u32,
    /// 0-based start coordinate of the bin (includes padding).
    start: i64,
    /// 0-based end coordinate of the bin (includes padding).
    end: i64,
    /// Indices into the original variant array.
    variant_indices: Vec<usize>,
}

/// Resolve a variant chromosome name to a BAM target id, tolerating contig-naming
/// differences between the variant source and the BAM (e.g. `chr1` vs `1`, `chrM` vs
/// `MT`) via `normalize_contig`. The exact name is tried first so the common (matching)
/// case stays O(1); the normalized scan of the header only runs on a miss. Returns
/// `None` only when no BAM contig matches even after normalization — i.e. the variant's
/// chromosome genuinely isn't in the BAM (a real mismatch worth surfacing, not silencing).
fn resolve_tid(bam_header: &bam::HeaderView, chrom: &str) -> Option<u32> {
    if let Some(t) = bam_header.tid(chrom.as_bytes()) {
        return Some(t);
    }
    let norm = crate::shared::contig::normalize_contig(chrom);
    (0..bam_header.target_count()).find(|&tid| {
        crate::shared::contig::normalize_contig(&String::from_utf8_lossy(bam_header.tid2name(tid)))
            == norm
    })
}

/// Build genomic bins from a list of variants.
///
/// Variants are grouped by chromosome and position into bins of at least
/// `BIN_WINDOW` (a floor that close variants chain past; see the module comment).
/// The algorithm is O(n) — a single pass over sorted variants. Bins are padded by
/// the widest member's read window, `max(5, repeat_span + 2)` (`window::pad_for_repeat_span`).
///
/// # Arguments
/// * `variants` — Variants (need not be sorted; sorting is internal)
/// * `bam_header` — BAM header for chromosome → tid lookup
/// * `window` — Bin window size in bp (production: `BIN_WINDOW`)
/// * `max_variants` — Variants per bin before a split (production: `BIN_MAX_VARIANTS`)
fn build_genomic_bins(
    variants: &[Variant],
    bam_header: &bam::HeaderView,
    window: i64,
    max_variants: usize,
) -> Vec<GenomicBin> {
    if variants.is_empty() {
        return Vec::new();
    }

    // Build index sorted by (chrom, pos) — we sort indices, not variants,
    // to preserve the original variant order for output.
    let mut sorted_indices: Vec<usize> = (0..variants.len()).collect();
    sorted_indices.sort_by(|&a, &b| {
        variants[a].chrom.cmp(&variants[b].chrom)
            .then(variants[a].pos.cmp(&variants[b].pos))
    });

    // A variant's reads can extend to `pos + ref_allele.len()` — the right
    // breakpoint of a deletion. The bin's fetch end must cover that for EVERY
    // variant in the bin, including the anchor; half a window of slack matches
    // the per-variant extension below. (Seeding the end at only
    // `bin_start + window` under-fetched a bin anchored by a deletion whose ref
    // span exceeds the window, dropping reads aligned past it.)
    let span_end = |idx: usize| -> i64 {
        (variants[idx].pos + variants[idx].ref_allele.len() as i64).saturating_add(window / 2)
    };

    let mut bins: Vec<GenomicBin> = Vec::new();
    let mut i = 0;

    while i < sorted_indices.len() {
        let first_idx = sorted_indices[i];
        let chrom = &variants[first_idx].chrom;
        let tid = match resolve_tid(bam_header, chrom) {
            Some(t) => t,
            None => {
                // Genuinely absent from the BAM (even after contig-name normalization):
                // warn ONCE per chromosome (variants are chrom-sorted) so a naming
                // mismatch surfaces instead of silently producing zero counts.
                warn!(
                    "Variant chromosome '{}' not found in the BAM header (checked with \
                     contig-name normalization); its variants will get zero counts. Likely a \
                     contig-naming mismatch between the variant file and the BAM.",
                    chrom
                );
                while i < sorted_indices.len() && &variants[sorted_indices[i]].chrom == chrom {
                    i += 1;
                }
                continue;
            }
        };

        let bin_start = variants[first_idx].pos;
        // Cover at least one window, but also the anchor variant's full ref span.
        let mut bin_end = bin_start.saturating_add(window).max(span_end(first_idx));
        let mut indices = vec![first_idx];
        let mut max_repeat_span: i64 = variants[first_idx].repeat_span as i64;

        // Extend bin while next variant is on same chrom and within window
        let mut j = i + 1;
        while j < sorted_indices.len() {
            let jdx = sorted_indices[j];
            if variants[jdx].chrom != *chrom {
                break;
            }
            if variants[jdx].pos >= bin_end {
                break;
            }
            // Enforce max variants per bin — split if exceeded
            if indices.len() >= max_variants {
                debug!(
                    "Bin split: {} variants reached the per-bin cap {} at {}:{}",
                    indices.len(), max_variants,
                    chrom, variants[jdx].pos + 1,
                );
                break;
            }
            indices.push(jdx);
            max_repeat_span = max_repeat_span.max(variants[jdx].repeat_span as i64);
            // Extend bin end to cover this variant's full ref span.
            bin_end = bin_end.max(span_end(jdx));
            j += 1;
        }

        // Pad by the widest member's read window (`window::read_window`), the
        // filter each variant applies to the cached reads.
        let padding = window::pad_for_repeat_span(max_repeat_span);
        bins.push(GenomicBin {
            tid,
            start: (bin_start - padding).max(0),
            end: bin_end.saturating_add(padding),
            variant_indices: indices,
        });

        i = j;
    }

    info!(
        "Built {} genomic bins from {} variants (window={}bp)",
        bins.len(),
        variants.len(),
        window,
    );

    bins
}


/// Alignment backend selection for Phase 3 fallback classification.
///
/// Controls which algorithm is used when variant-type-specific checkers
/// (SNP, Ins, Del, MNP, Complex) need to fall back to haplotype-level
/// alignment for ambiguous reads.
///
/// Selectable via `--alignment-backend` CLI flag; `pairhmm` is the default at every layer.
///
/// There is deliberately **no `Default` impl**. The backend materially changes indel calls
/// (measured: 5 of 104 molecule calls differ on real ACCESS deletions), so every entry point
/// states its choice rather than inheriting one silently. A `#[default]` here previously
/// still named `SmithWaterman` long after `pairhmm` became the real default everywhere else.
#[derive(Clone, Debug, PartialEq)]
pub enum AlignmentBackend {
    /// Smith-Waterman with affine gap penalties. The default through v2.8.0; kept as an
    /// explicit opt-in for backend comparison and for callers wanting the cheaper classifier.
    /// Score-margin classification: (alt_score - ref_score) > 0 → ALT.
    SmithWaterman,
    /// PairHMM with BQ-aware emissions — **the default since v3.0.0**. WFA edit-distance
    /// triage (`wfa_router`) runs in front of it and resolves clear-cut reads directly;
    /// only ambiguous reads reach the full model. WFA is not a separate backend.
    /// LLR classification: log P(read|ALT) - log P(read|REF) > threshold → ALT.
    PairHMM {
        /// Log-likelihood ratio threshold for confident calls (default: 2.3 ≈ 10:1 odds).
        llr_threshold: f64,
        /// Gap-open probability (linear scale, default: 1e-4).
        gap_open: f64,
        /// Gap-extend probability (linear scale, default: 0.1).
        gap_extend: f64,
        /// Gap-open probability for repeat regions (linear scale, default: 1e-2).
        gap_open_repeat: f64,
        /// Gap-extend probability for repeat regions (linear scale, default: 0.5).
        gap_extend_repeat: f64,
    },
}

impl AlignmentBackend {
    /// PairHMM backend with the CLI's default parameters, for tests. In
    /// production the Python layer passes the parameters to `count_bam_binned()`.
    #[cfg(test)]
    pub fn pairhmm_default() -> Self {
        AlignmentBackend::PairHMM {
            llr_threshold: 2.3,
            gap_open: 1e-4,
            gap_extend: 0.1,
            gap_open_repeat: 1e-2,
            gap_extend_repeat: 0.5,
        }
    }
}

/// Parse the `alignment_backend` token, rejecting anything unrecognized.
///
/// Deliberately strict. This previously fell back to Smith-Waterman on any unmatched string,
/// so a typo or wrong casing (`"PairHMM"`, `"smith-waterman"`, `"pairhm"`) silently computed
/// with the other classifier and returned plausible-looking counts. The backends genuinely
/// disagree on ambiguous indels, so that failure was invisible and wrong rather than loud.
/// Mirrors the `strandedness` parse in the same call path, which has always erred loudly.
///
/// Accepts exactly the tokens the CLI enum can emit: `sw`, `hmm`, `pairhmm`.
fn parse_alignment_backend(
    token: &str,
    llr_threshold: f64,
    gap_open: f64,
    gap_extend: f64,
    gap_open_repeat: f64,
    gap_extend_repeat: f64,
) -> PyResult<AlignmentBackend> {
    match token {
        "hmm" | "pairhmm" => Ok(AlignmentBackend::PairHMM {
            llr_threshold,
            gap_open,
            gap_extend,
            gap_open_repeat,
            gap_extend_repeat,
        }),
        "sw" => Ok(AlignmentBackend::SmithWaterman),
        other => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "unknown alignment_backend {other:?}; expected \"pairhmm\" (alias \"hmm\") or \"sw\""
        ))),
    }
}


/// Bin-centric parallel BAM counting: groups variants into ~10kb genomic bins,
/// fetches reads once per bin, then classifies each read against every variant in
/// the bin. Bin geometry is performance only; counts never depend on it (the
/// binning-invariance tests vary it).
///
/// Shared implementation behind `count_bam_binned` (counts only) and
/// `count_bam_binned_observations` (counts + per-molecule rows). Both public entry
/// points delegate here, so the export cannot drift from the counts. `emit_obs=false`
/// allocates nothing and returns an empty observation Vec (invariant 3: output-aware,
/// no compute-then-discard).
#[allow(clippy::too_many_arguments)]
fn count_bam_binned_core(
    py: Python<'_>,
    emit_obs: bool,
    observations_path: Option<&str>,
    bam_path: String,
    mut variants: Vec<Variant>,
    mut decomposed: Vec<Option<Variant>>,
    min_mapq: u8,
    min_baseq: u8,
    filter_duplicates: bool,
    filter_secondary: bool,
    filter_supplementary: bool,
    filter_qc_failed: bool,
    filter_improper_pair: bool,
    filter_indel: bool,
    threads: usize,
    fragment_qual_threshold: u8,
    sibling_variants: Vec<Vec<Variant>>,
    alignment_backend: &str,
    hmm_llr_threshold: f64,
    hmm_gap_open: f64,
    hmm_gap_extend: f64,
    hmm_gap_open_repeat: f64,
    hmm_gap_extend_repeat: f64,
    apply_baq: bool,
    umi_tag: Option<&str>,
    mode: &str,
    enforce_strandedness: bool,
    strandedness: &str,
    mfsd: bool,
    rna_editing_db: Option<&str>,
    gtf_path: Option<&str>,
    reference_fasta: Option<&str>,
    library_type: &str,
    bin_window: Option<i64>,
    bin_max_variants: Option<i64>,
    warn_per_bam: bool,
) -> PyResult<(Vec<BaseCounts>, Vec<Observation>)> {
    let (window, max_variants) = bin_geometry(bin_window, bin_max_variants)?;
    // The decomposed and sibling lists run parallel to the variants: a shorter
    // one is padded (no twin, no siblings), a longer one is a caller error, never
    // silently cut.
    for (name, len) in [("decomposed", decomposed.len()), ("sibling_variants", sibling_variants.len())] {
        if len > variants.len() {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "{name} has {len} entries for {} variants: it must run parallel to them",
                variants.len(),
            )));
        }
    }
    decomposed.resize_with(variants.len(), || None);
    let mut sibling_variants = sibling_variants;
    sibling_variants.resize_with(variants.len(), Vec::new);
    let backend = parse_alignment_backend(
        alignment_backend,
        hmm_llr_threshold,
        hmm_gap_open,
        hmm_gap_extend,
        hmm_gap_open_repeat,
        hmm_gap_extend_repeat,
    )?;

    // ── Load RNA editing site database (if provided) ──
    // Loaded ONCE at init, then shared across all bins/threads via Arc.
    // DB-only strategy: flag = True only when the variant matches a REDIportal edit.
    // No DB = no editing flags (no pattern-matching guessing). Each entry carries the
    // catalogued (chrom, 0-based pos, ref base, edited base) so the flag can require
    // the variant's substitution to match, not merely its position.
    let editing_sites: Option<HashSet<(String, i64, u8, u8)>> = match rna_editing_db {
        Some(path) => {
            let sites = rna::build_rna_editing_set(path)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
                    format!("Failed to load RNA editing DB: {}", e)
                ))?;
            info!("Loaded {} RNA editing sites from {}", sites.len(), path);
            Some(sites)
        }
        None => None,
    };
    // Wrap in Arc for thread-safe sharing across rayon workers
    let editing_sites = std::sync::Arc::new(editing_sites);

    // ── Build GTF annotation index (if GTF provided, RNA mode only) ──
    // Loaded ONCE at init, then shared across all bins/threads via Arc.
    // Only builds for chromosomes that have variants (variant-guided streaming).
    let annotation: Option<std::sync::Arc<AnnotationIndex>> = match (mode, gtf_path) {
        ("rna", Some(path)) => {
            let variant_chroms: HashSet<String> = variants.iter()
                .map(|v| crate::shared::contig::normalize_contig(&v.chrom))
                .collect();
            let annot = crate::annotation::parse_gtf(path, &variant_chroms).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
                format!("Failed to load GTF annotation: {}", e)
            ))?;
            info!(
                "GTF annotation: built index from {} — {} exons, {} transcripts, {} chromosomes",
                path, annot.n_exons(), annot.n_transcripts(), annot.n_chromosomes(),
            );
            Some(std::sync::Arc::new(annot))
        }
        ("rna", None) => {
            debug!("GTF annotation: no GTF provided, annotation features disabled");
            None
        }
        _ => None,  // DNA mode: no annotation, no log noise
    };

    // Resolve each variant's gene strand from the annotation, so RNA strandedness
    // enforcement and sense/antisense partitioning have a strand to act on. The
    // engine reads `variant.gene_strand` per read; without this it stays None,
    // every read is treated as sense, and `enforce_strandedness` is a silent no-op.
    // Only fill when absent, so a strand supplied upstream is never overwritten.
    if let Some(ref annot) = annotation {
        let mut resolved = 0usize;
        for v in variants.iter_mut() {
            if v.gene_strand.is_none() {
                let key = crate::shared::contig::normalize_contig(&v.chrom);
                if let Some(strand) = annot.strand_at(&key, v.pos) {
                    v.gene_strand = Some(strand);
                    resolved += 1;
                }
            }
        }
        debug!(
            "Resolved gene strand from the annotation for {}/{} variants",
            resolved, variants.len(),
        );
    }
    // A decomposed twin shares its original's contig and position, so it takes
    // the same strand (resolved above or supplied upstream); without it the twin
    // counts antisense reads the original excludes wherever it wins the dual-count.
    for (v, twin) in variants.iter().zip(decomposed.iter_mut()) {
        if let Some(twin) = twin {
            if twin.gene_strand.is_none() {
                twin.gene_strand = v.gene_strand;
            }
        }
    }
    // Strandedness enforcement needs a gene strand per variant; variants left
    // without one (no GTF, an intergenic locus, genes of both strands over it, or
    // a GTF/variant contig mismatch) pass every read as sense. Say so once, naming
    // them, instead of enforcing nothing while the run banner claims
    // enforce_strandedness=true.
    if enforce_strandedness {
        let unresolved: Vec<&Variant> = variants.iter().filter(|v| v.gene_strand.is_none()).collect();
        if !unresolved.is_empty() {
            const NAMED: usize = 10;
            let named: Vec<String> = unresolved
                .iter()
                .take(NAMED)
                .map(|v| format!("{}:{} {}>{}", v.chrom, v.pos + 1, v.ref_allele, v.alt_allele))
                .collect();
            warn!(
                "--enforce-strandedness: {}/{} variants have no gene strand ({}) — \
                 strandedness is NOT enforced for them and antisense reads are counted: {}{}",
                unresolved.len(),
                variants.len(),
                if annotation.is_some() {
                    "no transcript spans the locus, transcripts of both strands do, or contig naming mismatch"
                } else {
                    "no --gtf provided"
                },
                named.join(", "),
                if unresolved.len() > NAMED { format!(" and {} more", unresolved.len() - NAMED) } else { String::new() },
            );
        }
    }

    // Store FASTA path for thread-local readers (used by ASJD motif classification)
    let fasta_path_owned: Option<String> = reference_fasta.map(|p| p.to_string());

    // Convert library_type string to boolean for amplicon mode.
    // In amplicon mode, R1/R2 hash to separate "fragments" (no consensus).
    let amplicon_mode = library_type == "amplicon";

    // Parse the RNA library strand protocol (loud error on an unknown token). Amplicon
    // libraries are not stranded, so force Unstranded there — reads must not be folded
    // under a protocol that doesn't apply (matches the CLI auto-disabling
    // --enforce-strandedness for amplicon).
    let strandedness = rna::Strandedness::from_protocol(strandedness)
        .map_err(PyErr::new::<pyo3::exceptions::PyValueError, _>)?;
    let strandedness = if amplicon_mode {
        rna::Strandedness::Unstranded
    } else {
        strandedness
    };

    if mode == "rna" {
        info!(
            "count_bam_binned: {} variants, apply_baq={}, umi_tag={:?}, backend={:?}, \
             rna_editing_db={}, gtf={}, library_type={}, strandedness={:?}, threads={}",
            variants.len(), apply_baq, umi_tag, alignment_backend,
            rna_editing_db.unwrap_or("none"),
            gtf_path.unwrap_or("none"),
            library_type,
            strandedness,
            threads,
        );
    } else {
        info!(
            "count_bam_binned: {} variants, apply_baq={}, umi_tag={:?}, backend={:?}, mfsd={}, threads={}",
            variants.len(), apply_baq, umi_tag, alignment_backend, mfsd, threads,
        );
    }
    warn_degraded_variants(&variants);

    // Build genomic bins using a temporary BAM/CRAM reader for header access
    let mut header_reader = bam::IndexedReader::from_path(&bam_path).map_err(|e| {
        pyo3::exceptions::PyRuntimeError::new_err(format!("Failed to open BAM/CRAM: {}", e))
    })?;
    // Set CRAM reference for header decoding (safe no-op for BAM)
    if let Some(ref fasta) = fasta_path_owned {
        header_reader.set_reference(fasta).map_err(|e| {
            pyo3::exceptions::PyRuntimeError::new_err(
                format!("Failed to set CRAM reference: {}", e)
            )
        })?;
    }
    let bins = build_genomic_bins(&variants, header_reader.header(), window, max_variants);

    let n = variants.len();

    // Kept for the post-run UMI-tag check (bam_path moves into the workers).
    let bam_label = bam_path.clone();

    // Parse UMI tag for thread-local use
    let umi_tag_owned: Option<[u8; 2]> = umi_tag.and_then(|tag| {
        let bytes = tag.as_bytes();
        if bytes.len() == 2 {
            Some([bytes[0], bytes[1]])
        } else {
            log::warn!("UMI tag '{}' is not 2 characters, ignoring", tag);
            None
        }
    });

    // Configure thread pool. `--threads` is the TOTAL budget for this process; all
    // bin-level parallelism stays within it (see shared::resolve_thread_budget).
    let threads = crate::shared::resolve_thread_budget(threads);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
            format!("Failed to build thread pool: {}", e)
        ))?;

    // Result array: one BaseCounts per variant, initialized to default
    let mut all_counts: Vec<BaseCounts> = (0..n).map(|_| BaseCounts::default()).collect();
    // With mFSD on, every row starts as a variant with no fragments: NaN means, LLR
    // and KS statistic. A row the bins never reach (its contig is absent from the
    // BAM) keeps that state, so it reads as "no test" and stays out of the BH family
    // below; the 0.0 defaults would pass for a test with p = 0 and deflate every
    // real variant's q-value. Counted rows overwrite it.
    if mfsd {
        for (c, v) in all_counts.iter_mut().zip(variants.iter()) {
            compute_mfsd_stats(c, Vec::new(), Vec::new(), Vec::new(), Vec::new(), v);
        }
    }

    // The Parquet rows echo (chrom, pos, ref, alt) so the file is self-describing — a bare
    // variant_index means nothing once the data outlives the call. `variants` is moved into
    // the rayon closure below, so keep the loci here. One small clone per call, not per row,
    // and only when a path was actually requested.
    let obs_loci: Option<Vec<Variant>> = if emit_obs && observations_path.is_some() {
        Some(variants.clone())
    } else {
        None
    };

    // Process bins in parallel, each bin does one bam.fetch()
    #[allow(clippy::type_complexity)]
    let bin_results: Result<BinOutput, anyhow::Error> = py.detach(move || {
        pool.install(|| {
            bins.par_iter()
                .map_init(
                    || {
                        let bam_reader = (|| -> Result<bam::IndexedReader, anyhow::Error> {
                            let mut reader = bam::IndexedReader::from_path(&bam_path).map_err(|e| {
                                anyhow::anyhow!("Failed to open BAM/CRAM: {}", e)
                            })?;
                            // CRAM files require a reference FASTA for decoding.
                            // set_reference() is a safe no-op for BAM files.
                            if let Some(ref fasta) = fasta_path_owned {
                                reader.set_reference(fasta).map_err(|e| {
                                    anyhow::anyhow!("Failed to set CRAM reference: {}", e)
                                })?;
                            }
                            Ok(reader)
                        })();
                        // Thread-local FASTA reader for splice motif classification
                        // (opened whenever a FASTA path is provided)
                        let fasta_reader: Option<bio::io::fasta::IndexedReader<std::fs::File>> =
                            fasta_path_owned.as_ref().and_then(|path| {
                                bio::io::fasta::IndexedReader::from_file(path).ok()
                            });
                        // RNA: a cached reader for the far exon of a spliced read's
                        // junctions (the exact-carrier rule reads across them). A
                        // worker that cannot open it fails the run: judging its bins
                        // differently from the others' would be silent.
                        let far_reference: Result<Option<crate::normalize::fasta::CachedFasta>, String> =
                            match (mode == "rna", fasta_path_owned.as_deref()) {
                                (true, Some(path)) => crate::normalize::fasta::CachedFasta::open(path)
                                    .map(Some)
                                    .ok_or_else(|| format!("cannot open the reference FASTA {path}")),
                                _ => Ok(None),
                            };
                        (bam_reader, fasta_reader, far_reference)
                    },
                    |(bam_result, fasta_reader, far_reference), bin| {
                        let bam = match bam_result {
                            Ok(b) => b,
                            Err(e) => return Err(anyhow::anyhow!("BAM init failed: {}", e)),
                        };
                        let far_reference = match far_reference {
                            Ok(r) => r.as_ref(),
                            Err(e) => return Err(anyhow::anyhow!("reference init failed: {}", e)),
                        };

                        debug!(
                            "Processing bin tid={} {}-{} ({} variants)",
                            bin.tid, bin.start, bin.end,
                            bin.variant_indices.len(),
                        );

                        // ── Shared-read optimization ──
                        // Single bam.fetch() per bin; reads shared across all
                        // variants via count_bin_shared.
                        #[allow(clippy::needless_question_mark)]
                        Ok(count_bin_shared(
                            bam,
                            bin,
                            &variants,
                            &decomposed,
                            &sibling_variants,
                            min_mapq,
                            min_baseq,
                            filter_duplicates,
                            filter_secondary,
                            filter_supplementary,
                            filter_qc_failed,
                            filter_improper_pair,
                            filter_indel,
                            fragment_qual_threshold,
                            &backend,
                            apply_baq,
                            umi_tag_owned,
                            mode,
                            enforce_strandedness,
                            strandedness,
                            mfsd,
                            &editing_sites,
                            &annotation,
                            fasta_reader,
                            far_reference,
                            amplicon_mode,
                            emit_obs,
                        )?)
                    },
                )
                // Each bin owns its Vecs; they are concatenated here. A shared `&mut`
                // sink cannot cross the rayon closure, and a Mutex would serialize the
                // hot deep-coverage bins — so per-bin ownership + merge is the shape.
                // (`bins.par_iter()` is indexed, so this reduce splits contiguously and
                // joins adjacent results: bin ORDER is preserved. The nondeterminism the
                // sort below fixes comes from `HashMap` iteration *within* a variant, not
                // from here.)
                .try_reduce(
                    || (Vec::new(), Vec::new(), (HashSet::new(), HashSet::new())),
                    |mut acc: BinOutput, batch| {
                        acc.0.extend(batch.0);
                        acc.1.extend(batch.1);
                        acc.2 .0.extend(batch.2 .0);
                        acc.2 .1.extend(batch.2 .1);
                        Ok(acc)
                    },
                )
        })
    });

    // Scatter results back to variant-order array
    match bin_results {
        Ok((pairs, mut observations, (no_bases, no_quals))) => {
            for (vi, counts) in pairs {
                all_counts[vi] = counts;
            }

            // Which rule decided the pass's depth reads, once per BAM (logs only).
            let mut t = crate::types::DecisionTally::default();
            for c in &all_counts {
                t.add(&c.decisions);
            }
            let sw: u64 = all_counts.iter().map(|c| c.sw_fallback_reads as u64).sum();
            info!(
                "{}: depth reads by deciding rule over {} variant(s): withdrawn as uninformative \
                 {} REF, {} ALT ({} unjudged); ALT kept by its bases {}; exact-carrier rule {} \
                 judged, {} fell back; sibling guards {} REF excluded, {} ALT claimed; clip \
                 admissions {}; SW fallback {}",
                bam_label, all_counts.len(), t.ref_withdrawn, t.alt_withdrawn, t.alt_unjudged,
                t.alt_by_bases, t.carrier_judged, t.carrier_fallback, t.sibling_ref_excluded,
                t.sibling_alt_claimed, t.clip_admitted, sw,
            );

            // Records stored without bases (SEQ '*') show no allele, so the read
            // filter drops them. Say so once per BAM: a BAM stripped of its sequences
            // would otherwise count nothing with no word above DEBUG. Each record is
            // counted once, however many overlapping bins fetched it. A second pass
            // over the same BAM (the --rescue-mnp recount) runs with warn_per_bam off,
            // so the main pass's warnings stand for the BAM.
            if warn_per_bam && !no_bases.is_empty() {
                warn!(
                    "{}: skipped {} record(s) stored without bases (SEQ '*'); they show no \
                     allele and are not counted",
                    bam_label, no_bases.len(),
                );
            }
            if warn_per_bam && !no_quals.is_empty() {
                warn!(
                    "{}: skipped {} record(s) stored without base qualities (QUAL '*'); their \
                     bases have no stated quality and are not counted",
                    bam_label, no_quals.len(),
                );
            }

            // A requested UMI tag that no processed read carries means fragment
            // grouping silently fell back to QNAME for this BAM. Counts are
            // unaffected (QNAME grouping is the no-UMI behaviour), and mixed
            // tagged/untagged workflows are legitimate — so warn, don't fail.
            if let (Some(tag), Some(tag_bytes)) = (umi_tag, umi_tag_owned) {
                let tagged: u64 = all_counts.iter().map(|c| c.umi_tagged_reads as u64).sum();
                let depth: u64 = all_counts.iter().map(|c| c.dp as u64).sum();
                if warn_per_bam && tagged == 0 && depth > 0 {
                    warn!(
                        "--umi-tag {}: no processed read in {} carries the {} tag; fragment \
                         grouping fell back to read names (QNAME) for this BAM",
                        tag, bam_label, String::from_utf8_lossy(&tag_bytes),
                    );
                }
            }

            // Deterministic observation order. `HashMap<u64, FragmentEvidence>` iteration
            // inside each variant is nondeterministic (RandomState), and it does not affect
            // counts — those are order-invariant sums — so comparing counts never shows it.
            // Sorting by the join key makes the emitted rows reproducible run-to-run
            // (verified on real BAMs across 1/4/8 threads).
            if emit_obs {
                observations.sort_by_key(|o| (o.variant_index, o.molecule_hash));
            }

            // At-scale path: stream the rows to Parquet from here rather than returning
            // them. A panel- or genome-wide run yields 10^6–10^7 observations, where the
            // FFI materialization (10^7 PyObjects, roughly doubled peak RSS) — not the
            // counting — is the bottleneck. Same reasoning as `write_fsd_parquet`.
            // The returned Vec is then emptied, so the caller never pays that cost.
            if let (Some(path), Some(loci)) = (observations_path, obs_loci.as_ref()) {
                crate::counting::parquet_writer::write_observations_parquet(
                    path,
                    &observations,
                    loci,
                )?;
                observations.clear();
                observations.shrink_to_fit();
            }

            // ── BH-FDR correction for ASJD p-values ──
            // Requires all p-values simultaneously, so must run after all bins
            // are processed. Only runs in RNA mode with annotation present.
            // Correct ONLY over variants with a real ASJD test (junction
            // reads present). Non-junction variants keep the default asjd_pval and
            // would pad the family — inflating n and distorting the FDR — the same
            // bug fixed for the mFSD alt-vs-REF correction. Filtering on junction-read
            // presence also subsumes the old `p < 1.0` data guard.
            if mode == "rna" {
                let asjd_valid: Vec<usize> = all_counts
                    .iter()
                    .enumerate()
                    .filter(|(_, c)| c.asjd_n_alt_junc + c.asjd_n_ref_junc > 0)
                    .map(|(i, _)| i)
                    .collect();
                let corrected = crate::shared::stats::benjamini_hochberg_family(
                    |i| all_counts[i].asjd_pval,
                    &asjd_valid,
                );
                for (i, q) in &corrected {
                    all_counts[*i].asjd_qval = *q;
                }
                if !corrected.is_empty() {
                    debug!(
                        "ASJD BH-FDR: corrected {} ASJD p-values (over real tests)",
                        corrected.len(),
                    );
                }
            }

            // ── BH-FDR correction for the mFSD alt-vs-REF KS p-values ──
            // This p-value drives the report's LEANS-SOMATIC class, so correct
            // it for multiplicity across the sample. BH runs ONLY over variants whose
            // KS test actually RAN — a variant with too few fragments returns
            // (D = NaN, p = 1.0) from ks_test, so the *D-statistic* (not the p-value)
            // marks a real test. Filtering on the D-stat avoids padding the family
            // with no-test variants (which would inflate n and over-correct);
            // a genuine D=0 test legitimately keeps its p=1.0.
            //
            // When `--mfsd` is off, mFSD was never computed so the KS fields sit at
            // their BaseCounts default (0.0, not NaN). The `mfsd &&` guard keeps the family
            // empty in that case — otherwise every variant would be treated as a real test
            // and BH-corrected over default p-values.
            let mfsd_valid: Vec<usize> = all_counts
                .iter()
                .enumerate()
                .filter(|(_, c)| mfsd && !c.mfsd_ks_alt_ref.is_nan())
                .map(|(i, _)| i)
                .collect();
            let corrected = crate::shared::stats::benjamini_hochberg_family(
                |i| all_counts[i].mfsd_pval_alt_ref,
                &mfsd_valid,
            );
            for (i, q) in &corrected {
                all_counts[*i].mfsd_qval_alt_ref = *q;
            }
            if !corrected.is_empty() {
                debug!(
                    "mFSD BH-FDR: corrected {} mFSD alt-vs-REF p-values",
                    corrected.len(),
                );
            }

            Ok((all_counts, observations))
        }
        Err(e) => Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e))),
    }
}

/// Bin-centric parallel BAM counting with BAQ and UMI support.
///
/// Returns `list[BaseCounts]`, one per variant in input order. Delegates to
/// `count_bam_binned_core` with observations off. For per-molecule rows use
/// `count_bam_binned_observations`. `bin_window` / `bin_max_variants` are test
/// arguments (production passes neither): counts must not depend on them.
#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (bam_path, variants, decomposed, min_mapq, min_baseq, filter_duplicates, filter_secondary, filter_supplementary, filter_qc_failed, filter_improper_pair, filter_indel, threads, fragment_qual_threshold=10, sibling_variants=Vec::new(), alignment_backend="pairhmm", hmm_llr_threshold=2.3, hmm_gap_open=1e-4, hmm_gap_extend=0.1, hmm_gap_open_repeat=1e-2, hmm_gap_extend_repeat=0.5, apply_baq=false, umi_tag=None, mode="dna", enforce_strandedness=false, strandedness="reverse", mfsd=false, rna_editing_db=None, gtf_path=None, reference_fasta=None, library_type="capture", bin_window=None, bin_max_variants=None, warn_per_bam=true))]
pub fn count_bam_binned(
    py: Python<'_>,
    bam_path: String,
    variants: Vec<Variant>,
    decomposed: Vec<Option<Variant>>,
    min_mapq: u8,
    min_baseq: u8,
    filter_duplicates: bool,
    filter_secondary: bool,
    filter_supplementary: bool,
    filter_qc_failed: bool,
    filter_improper_pair: bool,
    filter_indel: bool,
    threads: usize,
    fragment_qual_threshold: u8,
    sibling_variants: Vec<Vec<Variant>>,
    alignment_backend: &str,
    hmm_llr_threshold: f64,
    hmm_gap_open: f64,
    hmm_gap_extend: f64,
    hmm_gap_open_repeat: f64,
    hmm_gap_extend_repeat: f64,
    apply_baq: bool,
    umi_tag: Option<&str>,
    mode: &str,
    enforce_strandedness: bool,
    strandedness: &str,
    mfsd: bool,
    rna_editing_db: Option<&str>,
    gtf_path: Option<&str>,
    reference_fasta: Option<&str>,
    library_type: &str,
    bin_window: Option<i64>,
    bin_max_variants: Option<i64>,
    warn_per_bam: bool,
) -> PyResult<Vec<BaseCounts>> {
    let (counts, _observations) = count_bam_binned_core(
        py, false, None, bam_path, variants, decomposed, min_mapq, min_baseq, filter_duplicates,
        filter_secondary, filter_supplementary, filter_qc_failed, filter_improper_pair,
        filter_indel, threads, fragment_qual_threshold, sibling_variants, alignment_backend,
        hmm_llr_threshold, hmm_gap_open, hmm_gap_extend, hmm_gap_open_repeat,
        hmm_gap_extend_repeat, apply_baq, umi_tag, mode, enforce_strandedness, strandedness,
        mfsd, rna_editing_db, gtf_path, reference_fasta, library_type,
        bin_window, bin_max_variants, warn_per_bam,
    )?;
    Ok(counts)
}

/// Same counting as `count_bam_binned`, plus the per-molecule allele calls it normally
/// aggregates away — returned as `(list[BaseCounts], list[Observation])` from one pass.
///
/// One `Observation` per fragment per variant, carrying the molecule identity, the resolved
/// allele, and its best supporting base quality. The same molecule seen at two variants in
/// this call shares a `molecule_hash`, which is what lets a caller link alleles across loci;
/// gbcms does no such linking itself. Rows are sorted by `(variant_index, molecule_hash)`,
/// so output is deterministic despite parallel bin processing.
///
/// Counts are byte-identical to `count_bam_binned` — same core, same classifier.
#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (bam_path, variants, decomposed, min_mapq, min_baseq, filter_duplicates, filter_secondary, filter_supplementary, filter_qc_failed, filter_improper_pair, filter_indel, threads, fragment_qual_threshold=10, sibling_variants=Vec::new(), alignment_backend="pairhmm", hmm_llr_threshold=2.3, hmm_gap_open=1e-4, hmm_gap_extend=0.1, hmm_gap_open_repeat=1e-2, hmm_gap_extend_repeat=0.5, apply_baq=false, umi_tag=None, mode="dna", enforce_strandedness=false, strandedness="reverse", mfsd=false, rna_editing_db=None, gtf_path=None, reference_fasta=None, library_type="capture", observations_path=None, bin_window=None, bin_max_variants=None))]
pub fn count_bam_binned_observations(
    py: Python<'_>,
    bam_path: String,
    variants: Vec<Variant>,
    decomposed: Vec<Option<Variant>>,
    min_mapq: u8,
    min_baseq: u8,
    filter_duplicates: bool,
    filter_secondary: bool,
    filter_supplementary: bool,
    filter_qc_failed: bool,
    filter_improper_pair: bool,
    filter_indel: bool,
    threads: usize,
    fragment_qual_threshold: u8,
    sibling_variants: Vec<Vec<Variant>>,
    alignment_backend: &str,
    hmm_llr_threshold: f64,
    hmm_gap_open: f64,
    hmm_gap_extend: f64,
    hmm_gap_open_repeat: f64,
    hmm_gap_extend_repeat: f64,
    apply_baq: bool,
    umi_tag: Option<&str>,
    mode: &str,
    enforce_strandedness: bool,
    strandedness: &str,
    mfsd: bool,
    rna_editing_db: Option<&str>,
    gtf_path: Option<&str>,
    reference_fasta: Option<&str>,
    library_type: &str,
    observations_path: Option<&str>,
    bin_window: Option<i64>,
    bin_max_variants: Option<i64>,
) -> PyResult<(Vec<BaseCounts>, Vec<Observation>)> {
    count_bam_binned_core(
        py, true, observations_path, bam_path, variants, decomposed, min_mapq, min_baseq, filter_duplicates,
        filter_secondary, filter_supplementary, filter_qc_failed, filter_improper_pair,
        filter_indel, threads, fragment_qual_threshold, sibling_variants, alignment_backend,
        hmm_llr_threshold, hmm_gap_open, hmm_gap_extend, hmm_gap_open_repeat,
        hmm_gap_extend_repeat, apply_baq, umi_tag, mode, enforce_strandedness, strandedness,
        mfsd, rna_editing_db, gtf_path, reference_fasta, library_type,
        bin_window, bin_max_variants, true,
    )
}

// ── Shared-Read Bin Processing ─────────────────────────────────────────────
//
// A port of the original C++ GBCMS block-processing architecture.
// Instead of calling bam.fetch() per variant, we fetch once for the entire
// bin and share the read buffer across all variants. This eliminates
// redundant I/O for co-located variants.
//
// Architecture:
//   Phase 0: Single bam.fetch(bin.start..bin.end) → apply universal filters
//            → store passing reads in Vec<Record> (the read cache)
//   Phase 1: For each variant, iterate the cached reads and classify
//            (allele check, fragment tracking, QC metrics)
//
// The original C++ GBCMS used:
//   my_bam_reader.SetRegion(block_start, block_end);
//   vector<BamAlignment> bam_vec;       // ← our Vec<Record>
//   while(GetNextAlignment(al)) { bam_vec.push_back(al); }
//   for (variant in block) { baseCountSNP(variant, bam_vec, ...); }

/// Process all variants in a genomic bin using a shared read cache.
///
/// Fetches reads once for the entire bin region, applies universal filters,
/// then classifies each read against each variant in the bin (the shared-read
/// optimization).
///
/// # Arguments
///
/// * `bam` — Indexed BAM reader (mutable for fetch)
/// * `bin` — Genomic bin with variant indices
/// * `variants` — All variants (bin.variant_indices indexes into this)
/// * `decomposed` — Parallel array of decomposed variants for dual-counting
/// * `sibling_variants` — Per-variant sibling arrays for multi-allelic guard
/// * All filter/config params — as passed to `count_bam_binned`
///
/// # Returns
///
/// Vec of (variant_index, BaseCounts) pairs for scatter-back into results.
///
/// # Errors
///
/// Returns error if BAM fetch fails or record reading fails.
/// No silent failures — all I/O errors are propagated.
#[allow(clippy::too_many_arguments)]
fn count_bin_shared(
    bam: &mut bam::IndexedReader,
    bin: &GenomicBin,
    variants: &[Variant],
    decomposed: &[Option<Variant>],
    sibling_variants: &[Vec<Variant>],
    min_mapq: u8,
    min_baseq: u8,
    filter_duplicates: bool,
    filter_secondary: bool,
    filter_supplementary: bool,
    filter_qc_failed: bool,
    filter_improper_pair: bool,
    filter_indel: bool,
    fragment_qual_threshold: u8,
    backend: &AlignmentBackend,
    apply_baq: bool,
    umi_tag: Option<[u8; 2]>,
    mode: &str,
    enforce_strandedness: bool,
    strandedness: rna::Strandedness,
    mfsd: bool,
    editing_sites: &Option<HashSet<(String, i64, u8, u8)>>,
    annotation: &Option<std::sync::Arc<AnnotationIndex>>,
    fasta_reader: &mut Option<bio::io::fasta::IndexedReader<std::fs::File>>,
    far_reference: Option<&crate::normalize::fasta::CachedFasta>,
    amplicon_mode: bool,
    emit_obs: bool,
) -> Result<BinOutput> {

    // ══════════════════════════════════════════════════════════════════════
    // PHASE 0: Single fetch + universal filtering → read cache
    // ══════════════════════════════════════════════════════════════════════
    //
    // Universal filters are applied ONCE here. Per-variant filters
    // (strandedness, anchor overlap) are applied in Phase 1 because they
    // depend on variant-specific properties (gene_strand, variant.pos).

    bam.fetch((bin.tid, bin.start, bin.end))
        .context("count_bin_shared: failed to fetch bin region")?;

    let mut read_cache: Vec<Record> = Vec::new();
    let mut mapq_filtered: u64 = 0;

    // Construct ReadFilter from the boolean params passed in.
    let read_filter = crate::shared::filters::ReadFilter {
        filter_duplicates,
        filter_secondary,
        filter_supplementary,
        filter_qc_failed,
        filter_improper_pair,
        filter_indel,
    };
    let mut filter_counts = crate::shared::filters::FilterCounts::default();
    let (mut no_bases_seen, mut no_quals_seen) = (HashSet::new(), HashSet::new());

    for result in bam.records() {
        let record = result.context("count_bin_shared: error reading BAM record")?;

        // Universal flag filters (delegated to shared::filters::ReadFilter).
        let (no_bases_before, no_quals_before) = (filter_counts.no_bases, filter_counts.no_quals);
        if !read_filter.passes(&record, &mut filter_counts) {
            if filter_counts.no_bases != no_bases_before {
                no_bases_seen.insert(record_identity(&record));
            } else if filter_counts.no_quals != no_quals_before {
                no_quals_seen.insert(record_identity(&record));
            }
            continue;
        }

        // MAPQ filter: DNA mode uses simple threshold, RNA mode uses NH rescue.
        // NOTE: strandedness is NOT filtered here — gene_strand varies per variant.
        // NOTE: MAPQ=0 reads are KEPT in the cache for MQ0 tracking in Phase 1.
        //       They are counted per-variant (only for overlapping reads) then
        //       skipped for classification, matching GATK MappingQualityZero
        //       behavior where MQ0 is tracked BEFORE the MAPQ filter.
        if mode == "rna" {
            // In RNA mode, MAPQ=0 reads are rescued by NH:i:1 tag.
            // is_valid_rna_alignment handles this: returns true for MAPQ≥threshold
            // OR for MAPQ<threshold with NH:i:1. So MAPQ=0+NH:i:1 stays in cache.
            // Only filter reads that truly fail RNA validation (MAPQ<threshold AND NH>1).
            if record.mapq() > 0 && !rna::is_valid_rna_alignment(&record, min_mapq) {
                mapq_filtered += 1;
                continue;
            }
            // MAPQ=0 reads: keep in cache for MQ0 tracking; will be handled in Phase 1
            if record.mapq() == 0 && !rna::is_valid_rna_alignment(&record, min_mapq) {
                // Still keep in cache — MQ0 tracking needs them
                // Phase 1 will count them and then skip classification
            }
        } else if record.mapq() < min_mapq && record.mapq() > 0 {
            // DNA mode: filter reads with 0 < MAPQ < min_mapq.
            // MAPQ=0 reads are KEPT for MQ0 tracking in Phase 1.
            mapq_filtered += 1;
            continue;
        }

        // A read ends at its fragment end: read-through bases past it are adapter,
        // hard-clipped here, as if trimmed, so that no rule sees them as bases or reach.
        let Some(record) = crate::shared::bam_utils::clip_to_fragment(record) else {
            trace!("clip_to_fragment: no aligned base inside the fragment");
            continue;
        };
        read_cache.push(record);
    }

    debug!(
        "Bin tid={} {}-{}: {} reads cached ({} filtered: dup={} sec={} supp={} qc={} pair={} indel={} no_bases={} unmapped={} no_quals={} mapq={})",
        bin.tid, bin.start, bin.end, read_cache.len(),
        filter_counts.total() + mapq_filtered,
        filter_counts.duplicates, filter_counts.secondary, filter_counts.supplementary,
        filter_counts.qc_failed, filter_counts.improper_pair, filter_counts.indel,
        filter_counts.no_bases, filter_counts.unmapped, filter_counts.no_quals, mapq_filtered,
    );

    // ══════════════════════════════════════════════════════════════════════
    // PHASE 1: Per-variant classification from cached reads
    // ══════════════════════════════════════════════════════════════════════
    //
    // For each variant (and its decomposed twin), iterate the read cache and
    // classify each read. Per-variant filters (strandedness, anchor overlap) are
    // applied here.

    let mut results = Vec::with_capacity(bin.variant_indices.len());
    let mut bin_observations: Vec<Observation> = Vec::new();

    for &vi in &bin.variant_indices {
        let variant = &variants[vi];
        let siblings = &sibling_variants[vi];
        // How the exact-carrier rule reads past a read's aligned blocks: in RNA a
        // soft clip reaching an exon edge or junction end is not evidence, and a
        // spliced read continues into its next exon.
        let edges = carrier::ClipEdges::new(|| clip_edges(variant, &read_cache, annotation));
        let spliced = carrier::SplicedCache::default();
        // A change recurring beside a complex event is masked in the exact-carrier
        // windows, and one the ALT reads alone carry names their larger allele; found
        // from the reads the counts read, on first use, for the main and
        // per-transcript counts alike. The decomposed form is another allele: it
        // finds its own.
        let cache = &read_cache;
        let guard_of = |v| {
            carrier::CarrierGuard::new(move || {
                let counted = |r: &Record| counted_read(r, v, mode, enforce_strandedness, strandedness, min_mapq);
                carrier::guard_for(cache, v, &counted, min_baseq)
            })
        };
        let guard = guard_of(variant);
        let rules = if mode == "rna" {
            carrier::ReadRules {
                clip_edges: Some(&edges),
                reference: far_reference,
                spliced: Some(&spliced),
                guard: Some(&guard),
            }
        } else {
            carrier::ReadRules { guard: Some(&guard), ..carrier::ReadRules::DNA }
        };

        let (counts_orig, obs_orig) = count_variant_from_cache(
            &read_cache, variant, siblings,
            min_mapq, min_baseq,
            filter_improper_pair, filter_indel,
            fragment_qual_threshold, backend,
            apply_baq, umi_tag, mode, enforce_strandedness, strandedness, mfsd,
            editing_sites, annotation, amplicon_mode,
            emit_obs, vi as u32, &rules,
        )?;

        // Dual-count for decomposed variants: run the same classification
        // against the decomposed variant form and take the higher ALT count.
        //
        // Observations follow the SAME arbitration as the counts. A decomposed variant
        // runs the classifier twice, so keeping both sets would emit two contradictory
        // allele calls per molecule (including the losing allele form) — and the counts
        // would stay right while the export was wrong. Only the winner's rows survive.
        let (mut final_counts, final_obs) = if let Some(ref decomp) = decomposed[vi] {
            let decomp_guard = guard_of(decomp);
            let (counts_decomp, obs_decomp) = count_variant_from_cache(
                &read_cache, decomp, siblings,
                min_mapq, min_baseq,
                filter_improper_pair, filter_indel,
                fragment_qual_threshold, backend,
                apply_baq, umi_tag, mode, enforce_strandedness, strandedness, mfsd,
                editing_sites, annotation, amplicon_mode,
                emit_obs, vi as u32, &carrier::ReadRules { guard: Some(&decomp_guard), ..rules },
            )?;

            if counts_decomp.ad > counts_orig.ad {
                (
                    BaseCounts {
                        used_decomposed: true,
                        ..counts_decomp
                    },
                    obs_decomp,
                )
            } else {
                (counts_orig, obs_orig)
            }
        } else {
            (counts_orig, obs_orig)
        };
        if emit_obs {
            bin_observations.extend(final_obs);
        }

        // ── Non-discriminating-locus detection (PairHMM backend) ──
        // When a sibling combination reconstructs the reference haplotype, REF and
        // ALT are sequence-indistinguishable and every read ties to NEITHER. Detect it
        // once per variant (the matrix is read-independent) and flag it, so the zeroed
        // RD/AD is explained by NON_DISCRIMINATING_LOCUS rather than left silent.
        if matches!(backend, AlignmentBackend::PairHMM { .. }) {
            if let Some(matrix) =
                crate::counting::pangenome::build_haplotype_matrix(variant, siblings)
            {
                if crate::counting::pangenome::has_ref_alt_collision(&matrix) {
                    final_counts.non_discriminating_locus = true;
                    warn!(
                        "non-discriminating locus at {}:{}: a sibling combination reconstructs \
                         the reference haplotype — REF and ALT are sequence-indistinguishable, \
                         reads tie to NEITHER",
                        variant.chrom, variant.pos + 1,
                    );
                }
            }
        }

        // ── Per-transcript counting (RNA + GTF only) ──
        // For each overlapping transcript, count reads whose splice junctions
        // are compatible with that transcript's intron structure. Reuses the
        // same read cache — no additional BAM I/O.
        if let Some(ref annot) = *annotation {
            // The main counts' BAQ rule on the main counts' own distance (a
            // decomposed form keeps the variant's contig, position and REF).
            let use_baq = baq_applies(apply_baq, final_counts.exon_boundary_dist);
            let (read_cts, frag_cts) = count_per_transcript(
                &read_cache, variant, siblings, annot,
                min_mapq, min_baseq, fragment_qual_threshold,
                backend, use_baq, umi_tag, enforce_strandedness, strandedness,
                amplicon_mode, &rules,
            );
            final_counts.transcript_read_counts = read_cts;
            final_counts.transcript_fragment_counts = frag_cts;

            // ── Allele-Specific Junction Divergence (ASJD) ──
            // Compare splice junction usage between REF and ALT reads.
            // asjd_qval initialized to asjd_pval; corrected post-counting
            // via benjamini_hochberg() in count_bam_binned().
            let asjd = detect_asjd(
                &read_cache, variant, siblings, annot,
                min_mapq, min_baseq, backend, use_baq, enforce_strandedness, strandedness,
                fasta_reader, &rules,
            );
            final_counts.asjd_flag = asjd.flag;
            final_counts.asjd_pval = asjd.pval;
            final_counts.asjd_qval = asjd.pval; // Placeholder — BH-corrected post-counting
            final_counts.asjd_ref_junction = asjd.ref_junction;
            final_counts.asjd_alt_junction = asjd.alt_junction;
            final_counts.asjd_ref_motif = asjd.ref_motif;
            final_counts.asjd_alt_motif = asjd.alt_motif;
            final_counts.asjd_ref_known = asjd.ref_known;
            final_counts.asjd_alt_known = asjd.alt_known;
            final_counts.asjd_n_ref_junc = asjd.n_ref_junc;
            final_counts.asjd_n_alt_junc = asjd.n_alt_junc;
            final_counts.asjd_n_ref_total = asjd.n_ref_total;
            final_counts.asjd_n_alt_total = asjd.n_alt_total;
            final_counts.asjd_diagnostic = asjd.diagnostic;
        }

        results.push((vi, final_counts));
    }

    Ok((results, bin_observations, (no_bases_seen, no_quals_seen)))
}


/// Compute the mFSD fragment-size statistics for one variant and store them on
/// `counts`: the four class counts, means, alt/ref LLR, the six pairwise KS triads
/// (delta/D/p), the sub-/mono-nucleosomal fractions, and the raw size arrays for
/// `--mfsd-parquet`. Consumes the size vectors. Called only when mFSD output is
/// requested (the engine is output-aware — see the `mfsd` gate in
/// `count_variant_from_cache`).
fn compute_mfsd_stats(
    counts: &mut BaseCounts,
    mut ref_sizes: Vec<f64>,
    mut alt_sizes: Vec<f64>,
    mut nonref_sizes: Vec<f64>,
    mut n_sizes: Vec<f64>,
    variant: &Variant,
) {
    // The sizes arrive in the fragments' hash order; floating-point sums depend on
    // order, so sort once and every statistic (and the size arrays written to
    // --mfsd-parquet) is the same bit for bit run to run.
    for sizes in [&mut ref_sizes, &mut alt_sizes, &mut nonref_sizes, &mut n_sizes] {
        sizes.sort_by(f64::total_cmp);
    }
    counts.mfsd_ref_count    = ref_sizes.len()    as u32;
    counts.mfsd_alt_count    = alt_sizes.len()    as u32;
    counts.mfsd_nonref_count = nonref_sizes.len() as u32;
    counts.mfsd_n_count      = n_sizes.len()      as u32;

    // An empty class has no mean size (NaN, written NA), not a mean of 0 bp.
    let mean = |v: &[f64]| if v.is_empty() { f64::NAN } else { mfsd::calc_mean(v) };
    counts.mfsd_ref_mean    = mean(&ref_sizes);
    counts.mfsd_alt_mean    = mean(&alt_sizes);
    counts.mfsd_nonref_mean = mean(&nonref_sizes);
    counts.mfsd_n_mean      = mean(&n_sizes);

    // Mean per fragment, so the value does not grow with depth; n is mfsd_*_count.
    counts.mfsd_alt_llr = mfsd::calc_llr_mean(&alt_sizes);
    counts.mfsd_ref_llr = mfsd::calc_llr_mean(&ref_sizes);

    // KS helper: pairwise delta + D-statistic + p-value.
    // delta = mean(a) - mean(b); ks_test returns (NaN, 1.0) when either class < MIN_FOR_KS.
    let ks_pair = |a: &[f64], b: &[f64]| -> (f64, f64, f64) {
        let (d, p) = mfsd::ks_test(a, b);
        let delta = if a.is_empty() || b.is_empty() {
            f64::NAN
        } else {
            mfsd::calc_mean(a) - mfsd::calc_mean(b)
        };
        (delta, d, p)
    };

    (counts.mfsd_delta_alt_ref,    counts.mfsd_ks_alt_ref,    counts.mfsd_pval_alt_ref)    = ks_pair(&alt_sizes, &ref_sizes);
    counts.mfsd_qval_alt_ref = counts.mfsd_pval_alt_ref; // placeholder; BH-corrected post-counting
    (counts.mfsd_delta_alt_nonref, counts.mfsd_ks_alt_nonref, counts.mfsd_pval_alt_nonref) = ks_pair(&alt_sizes, &nonref_sizes);
    (counts.mfsd_delta_ref_nonref, counts.mfsd_ks_ref_nonref, counts.mfsd_pval_ref_nonref) = ks_pair(&ref_sizes, &nonref_sizes);
    (counts.mfsd_delta_alt_n,      counts.mfsd_ks_alt_n,      counts.mfsd_pval_alt_n)      = ks_pair(&alt_sizes, &n_sizes);
    (counts.mfsd_delta_ref_n,      counts.mfsd_ks_ref_n,      counts.mfsd_pval_ref_n)      = ks_pair(&ref_sizes, &n_sizes);
    (counts.mfsd_delta_nonref_n,   counts.mfsd_ks_nonref_n,   counts.mfsd_pval_nonref_n)   = ks_pair(&nonref_sizes, &n_sizes);

    // Sub-nucleosomal (<150bp): ctDNA enrichment indicator. Mono-nucleosomal
    // (150–200bp): dominant cfDNA peak. Computed before the size arrays are consumed.
    counts.mfsd_sub_nuc_ref_frac = mfsd::calc_fraction_in_range(&ref_sizes, 0.0, 150.0);
    counts.mfsd_sub_nuc_alt_frac = mfsd::calc_fraction_in_range(&alt_sizes, 0.0, 150.0);
    counts.mfsd_sub_nuc_enrichment = if counts.mfsd_sub_nuc_ref_frac > 0.0 {
        counts.mfsd_sub_nuc_alt_frac / counts.mfsd_sub_nuc_ref_frac
    } else {
        f64::NAN
    };
    counts.mfsd_mono_nuc_ref_frac = mfsd::calc_fraction_in_range(&ref_sizes, 150.0, 200.0);
    counts.mfsd_mono_nuc_alt_frac = mfsd::calc_fraction_in_range(&alt_sizes, 150.0, 200.0);

    // Raw size arrays for --mfsd-parquet export (held on BaseCounts across all variants).
    counts.ref_sizes = ref_sizes.into_iter().map(|v| v as u32).collect();
    counts.alt_sizes = alt_sizes.into_iter().map(|v| v as u32).collect();

    debug!(
        "mFSD {}:{} {}>{}: ref={} alt={} nonref={} n={} delta={:.1} ks_d={:.3} ks_p={:.3e} alt_llr={:.2} sizing=physical",
        variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele,
        counts.mfsd_ref_count, counts.mfsd_alt_count,
        counts.mfsd_nonref_count, counts.mfsd_n_count,
        counts.mfsd_delta_alt_ref,
        counts.mfsd_ks_alt_ref,
        counts.mfsd_pval_alt_ref,
        counts.mfsd_alt_llr,
    );
}

/// Classify and count reads from a pre-fetched cache for a single variant.
///
/// Operates on the bin's `&[Record]` cache rather than a per-variant
/// `bam.fetch()`. Universal filters (dup, secondary, supp, QC, MAPQ) have
/// already been applied in Phase 0; only per-variant filters remain:
/// - Strandedness (gene_strand varies per variant)
/// - Anchor overlap (variant.pos varies per variant)
/// - Read window overlap (variant.pos ± window_pad)
///
/// # Performance
///
/// Each call creates its own per-variant state (SW aligners, fragment map,
/// dist vectors). The read cache is borrowed immutably — no cloning.
#[allow(clippy::too_many_arguments)]
fn count_variant_from_cache(
    read_cache: &[Record],
    variant: &Variant,
    sibling_variants: &[Variant],
    min_mapq: u8,        // Used in Phase 1 for MAPQ skip after MQ0 tracking
    min_baseq: u8,
    _filter_improper_pair: bool, // Already filtered in Phase 0
    _filter_indel: bool,         // Already filtered in Phase 0
    fragment_qual_threshold: u8,
    backend: &AlignmentBackend,
    apply_baq: bool,
    umi_tag: Option<[u8; 2]>,
    mode: &str,
    enforce_strandedness: bool,
    strandedness: rna::Strandedness,
    mfsd: bool,
    editing_sites: &Option<HashSet<(String, i64, u8, u8)>>,
    annotation: &Option<std::sync::Arc<AnnotationIndex>>,
    amplicon_mode: bool,
    // ── Observation export (additive; counting is untouched) ──
    // `emit_obs` is the output-aware gate (invariant 3): when false, no per-molecule
    // rows are built or allocated. `variant_index` is stamped onto each row so the
    // caller can join back to its input `variants` list.
    emit_obs: bool,
    variant_index: u32,
    rules: &carrier::ReadRules,
) -> Result<(BaseCounts, Vec<Observation>)> {

    let mut counts = BaseCounts::default();

    // ── Compute exon boundary distance (GTF-informed) ──
    // Set once per variant, not per read: the least distance from any base of the
    // variant's REF span (or the span it was given, `boundary_span`), 0 when a
    // boundary lies inside it, so an MNP reaching into an exon's last bases is as
    // close as its nearest base. Used for BAQ suppression and as an output column.
    // None when no GTF is provided.
    let (span_first, span_last) = variant.boundary_span.unwrap_or_else(|| {
        (variant.pos, variant.pos + variant.ref_allele.len().max(1) as i64 - 1)
    });
    let exon_boundary_dist: Option<i32> = annotation.as_ref().and_then(|annot| {
        annot.nearest_splice_distance(&variant.chrom, span_first, span_last)
    });
    counts.exon_boundary_dist = exon_boundary_dist;

    // The allele the reads carry when it is not the given one: diagnostic only
    // (OBSERVED_ALLELE); no count below depends on it. It reads the reads the
    // counts read, so its n/m compare with alt_count and ref_count.
    let counted = |record: &Record| counted_read(record, variant, mode, enforce_strandedness, strandedness, min_mapq);
    if let Some(o) = observed::observed_allele(read_cache, variant, sibling_variants, &counted, min_baseq) {
        counts.observed_pos = o.pos + 1;
        counts.observed_ref = o.ref_allele;
        counts.observed_alt = o.alt_allele;
        counts.observed_reads = o.carriers;
        counts.observed_given_reads = o.given_carriers;
    } else if let Some(o) = rules.guard.and_then(|g| g.get(variant)).and_then(|g| {
        g.larger.as_ref().and_then(|l| observed::larger_allele(l, g.scanned, variant, sibling_variants))
    }) {
        // The ALT reads carry the given allele with a change beside it: the counts
        // stay the given allele's; name the larger allele they carry.
        counts.observed_pos = o.pos + 1;
        counts.observed_ref = o.ref_allele;
        counts.observed_alt = o.alt_allele;
        counts.observed_reads = o.carriers;
        counts.observed_given_reads = o.given_carriers;
    }
    let use_baq = baq_applies(apply_baq, exon_boundary_dist);
    // The variant's own indel span that BAQ leaves alone, once per variant.
    let baq_spare = if use_baq { baq_own_span(variant) } else { None };
    if apply_baq && !use_baq {
        debug!(
            "BAQ skipped at {}:{} {}>{}: bases {}-{} are {}bp from an annotated exon boundary",
            variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele,
            span_first + 1, span_last + 1, exon_boundary_dist.unwrap_or_default(),
        );
    }

    // Fragment tracking: QNAME hash -> FragmentEvidence
    let mut fragments: HashMap<u64, FragmentEvidence> = HashMap::new();
    let qual_diff_threshold: u8 = fragment_qual_threshold;

    // Distance-to-read-end tracking for QC metrics
    let mut alt_dists: Vec<u32> = Vec::with_capacity(500);
    let mut ref_dists: Vec<u32> = Vec::with_capacity(500);

    // Create SW aligners ONCE per variant (indelpost pattern), scored by `sw_score`.
    let mut alt_aligner = Aligner::new(SW_GAP_OPEN, SW_GAP_EXTEND, &sw_score);
    let mut ref_aligner = Aligner::new(SW_GAP_OPEN, SW_GAP_EXTEND, &sw_score);
    // Siblings whose change lies inside this row's discrimination window:
    // the only ones whose carriers the REF guard excludes.
    let ref_guard_siblings = window::siblings_in_window(variant, sibling_variants);

    // Per-phase classification counters
    let mut phase_counts = [0u32; 5];

    // The variant's read window: which cached reads overlap this variant (the bin's
    // fetch holds it for every member; see `build_genomic_bins`)
    let (v_start, v_end) = window::read_window(variant);

    // ── ref_context is ALWAYS genomic: consensus introns are never spliced out.
    // Splicing it without a coordinate map would shrink the context while
    // ref_context_start stays genomic, so every `pos - ref_context_start`
    // indexer right of a snipped intron — the deleted-bases check, the
    // large-deletion band's context guard, haplotype offsets — would read the
    // wrong bases, and pre-mRNA / intron-retention reads, whose bases genuinely
    // include intron sequence, would be scored against a haplotype missing
    // those bases. Junction-spanning reads never reach Phase 3 scoring anyway:
    // the splice-aware evidence rule refuses to extract or string-compare
    // across an N. Splice-aware Phase-3 scoring, if real-data measurement shows
    // it is needed, requires an explicit genomic→spliced coordinate map with
    // junction-compatible extraction and translated variant/sibling offsets.
    let mut reads_considered = 0u32;

    for record in read_cache {
        // ── Read window filter: the bin's cache holds the reads of every member's
        // read window; each variant sees only the cached reads that overlap its own.
        let r_start = record.pos();
        let r_end = ref_end(record);
        if r_start >= v_end || r_end <= v_start {
            continue; // Read doesn't overlap the variant's read window
        }

        // Supplementary/secondary alignments share a QNAME with their primary, so they are
        // never *first-class* read-level observations: counting one toward DP/RD/AD would
        // report two reads where one physical read exists. That holds regardless of the
        // filter flags, and is what the CLI documents ("secondary/supplementary never count
        // toward read-level depth regardless").
        //
        // They ARE admitted to fragment evidence when the caller opts out via
        // `filter_supplementary` / `filter_secondary` (the flags gate cache admission
        // above). That cannot double-count: `hash_molecule` keys on QNAME, so a primary and
        // its supplementary collapse into ONE `FragmentEvidence`. What it does fix is a
        // locus covered *only* by a supplementary segment, which previously reported
        // `dpf=0` — not a filtered read but a wrong answer, and the reason a molecule
        // spanning a large deletion was invisible to cross-locus phasing.
        //
        // Previously an unconditional skip sat here, discarding these records before
        // fragment evidence too, which made both flags unable to change any output at all.
        let first_class = !(record.is_supplementary() || record.is_secondary());
        if first_class {
            reads_considered += 1;
        }

        // ── MQ0 TRACKING: Count MAPQ=0 reads BEFORE any MAPQ-based skip
        // AND before the strandedness filter — an antisense MAPQ-0 read is
        // still a physical read at the locus.
        // Mirrors GATK's MappingQualityZero annotation — a high MQ0 count
        // is a locus-level red flag for regions with high homology or
        // pseudogenes, even when those reads are filtered for classification.
        // Read-level, so first-class records only (see `first_class` above).
        if first_class && record.mapq() == 0 {
            counts.mq0_count += 1;
        }

        // ── RNA STRANDEDNESS FILTER: per-variant because gene_strand differs.
        // An antisense read under enforcement counts nowhere, but it is classified
        // below as a sense read would be and tallied in rna_antisense_depth when it
        // is a first-class REF or ALT read over the anchor, then dropped (as
        // mq0_count counts reads the MAPQ skip drops). The column then means the
        // same with and without enforcement.
        let antisense_excluded = mode == "rna"
            && enforce_strandedness
            && !rna::is_sense_strand(record, variant.gene_strand, strandedness);
        // Only a first-class read over the anchor can be tallied (RNA admits no read
        // by its soft clip), so any other excluded read stops here, unclassified.
        if antisense_excluded && !(first_class && r_start <= variant.pos && r_end > variant.pos) {
            continue;
        }

        // ── MAPQ SKIP (Phase 1): MAPQ=0 reads were kept in the cache
        // specifically for MQ0 tracking above. Now skip them for
        // classification — they should not contribute to DP/RD/AD/DPF.
        if !mapping_admits(record, mode, min_mapq) {
            continue;
        }

        // ── HEURISTIC BAQ: resolve adjusted qualities for classification.
        // When BAQ is enabled, bases near indels and splice junctions
        // (CIGAR N) are downgraded. Default: off for DNA (upstream BQSR),
        // on for RNA (no upstream BQ recalibration). Skipped at exon edges:
        // `use_baq` is `baq_applies`, resolved once per variant above.
        let quals = effective_quals(record, use_baq, baq_spare);
        let effective_quals: &[u8] = &quals;

        // ── Allele classification
        let result = check_allele_with_qual(
            record, variant, sibling_variants, effective_quals, min_baseq,
            &mut alt_aligner, &mut ref_aligner, backend, rules,
        );
        // ── SPLICE-SKIP EXCLUSION: covers_locus=false means the read's
        // CIGAR N spans every discriminating position — it observes nothing
        // here, so it contributes to neither DP nor fragment depth, and it
        // is not a classification (kept out of phase_counts). At skipped
        // positions this matches samtools pileup's zero coverage; at an
        // anchor-preserved deletion it is deliberately stricter than pileup
        // at POS (the anchor base may be aligned) — an unobservant read in
        // DP would only deflate VAF. Must run before the anchor-overlap
        // gate below, which is span-based (ref_end includes N) and
        // would otherwise admit these reads.
        if !result.covers_locus {
            if !antisense_excluded {
                counts.splice_skip_excluded += 1;
            }
            continue;
        }

        let base_qual = result.qual;

        // ── MULTI-ALLELIC AD-CLAIMING GUARD: an ALT match contested and won
        // by a co-annotated sibling is the sibling's molecule. Downgrade
        // it here — before distance tracking, fragment evidence, and
        // read-level counting — so AD and ADF exclude it consistently. The
        // read still counts toward DP/DPF and is recorded as
        // partial_alt/any_alt below.
        let claimed_by_sibling = sibling_claims_alt(
            record, variant, &result, sibling_variants, effective_quals, min_baseq,
        );
        let is_alt = result.is_alt && !claimed_by_sibling;
        // ── MULTI-ALLELIC REF GUARD: a REF read that is ALT for a sibling
        // whose change lies inside this row's discrimination window is the
        // sibling's molecule. Excluded here, beside the AD guard, so RD and
        // RDF drop the same molecules; recorded as partial_alt/any_alt below.
        let ref_claimed_by_sibling = sibling_claims_ref(
            record, variant, &result, &ref_guard_siblings, effective_quals, min_baseq,
            &mut alt_aligner, &mut ref_aligner, backend, rules,
        );
        let is_ref = result.is_ref && !ref_claimed_by_sibling;
        // One line per classified read, named, so a count can be traced back to
        // the reads behind it (read-level validation against the BAM).
        trace!(
            "read call {}:{} {}>{} read={} mate={} ref={} alt={} phase={:?} partial={} nearby={} \
             sibling_claimed={} ref_sibling_claimed={} uninformative={} antisense_excluded={} \
             rule={}",
            variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele,
            read_name(record),
            if record.is_first_in_template() { 1 } else { 2 },
            is_ref, is_alt, result.phase, result.partial_match_count,
            result.has_nearby_evidence, claimed_by_sibling, ref_claimed_by_sibling,
            result.uninformative, antisense_excluded, deciding_rule(&result),
        );

        // An antisense read excluded by strandedness (first-class, over the anchor:
        // checked above): tallied as the sense reads are at the end of this loop,
        // when REF or ALT, then dropped before any count.
        if antisense_excluded {
            if is_ref || is_alt {
                counts.antisense_depth += 1;
                if is_alt {
                    counts.antisense_strand_alt_count += 1;
                }
            }
            continue;
        }
        phase_counts[result.phase as usize] += 1;
        tally_clip_candidate(&mut counts, record, variant, first_class);

        // ── DISTANCE TO READ END: Track how close the variant-supporting
        // base is to the nearest end of the read. Bases near read ends
        // have higher error rates and misalignment probability.
        // Stored per-allele for median computation after the loop.
        if first_class && (is_ref || is_alt) {
            if let Some(read_idx) = find_read_pos(record, variant.pos) {
                let read_len = record.seq_len();
                let dist = std::cmp::min(read_idx, read_len.saturating_sub(1 + read_idx)) as u32;
                if is_alt { alt_dists.push(dist); }
                if is_ref { ref_dists.push(dist); }
            }
        }

        // ── ANCHOR OVERLAP CHECK (strict): DP, RD, and AD are all defined
        // exclusively as depth at the variant anchor position (VCF POS).
        // This matches samtools pileup, GATK FORMAT/DP, and VarDict conventions.
        //
        // Reads that are in the classification window (±window_pad) but do NOT
        // overlap the anchor are used solely for haplotype evidence during
        // allele classification above. They must NOT contribute to DP/RD/AD,
        // because:
        //   1. Their bases are not at the locus being reported.
        //   2. Including them inflates DP above the true pileup depth.
        //   3. It makes VAF (AD/DP) inconsistent with standard tools.
        //
        // REF+ALT ≤ DP is guaranteed because RD and AD are strict subsets
        // of the anchor-overlap read set counted in DP here.
        let overlaps_anchor = r_start <= variant.pos && r_end > variant.pos;
        // A read whose allele lies in its soft-clipped bases (an aligner clips an
        // ALT read near its end where the REF reads beside it align in full): the
        // exact-carrier rule decided it from its own bases inside its fragment,
        // so it counts like a read aligned over the anchor. DNA only: an RNA
        // read's clip may hold the next exon's bases.
        let clip_admitted = !overlaps_anchor && mode != "rna" && result.clip_admissible;
        // Decided from windows read from the right flank: a REF and an ALT molecule
        // hold them from the same starts, so one starting past the variant position
        // counts as the other allele's would.
        let right_admitted = !overlaps_anchor && result.read_from_right;
        if !overlaps_anchor && !clip_admitted && !right_admitted {
            continue;
        }
        if clip_admitted || right_admitted {
            trace!(
                "read {} admitted at {}:{} {}>{} by {} (ref={} alt={})",
                String::from_utf8_lossy(record.qname()), variant.chrom, variant.pos + 1,
                variant.ref_allele, variant.alt_allele,
                if right_admitted { "windows read from the right flank" } else { "its soft-clipped bases" },
                is_ref, is_alt,
            );
        }

        // ── TOTAL DEPTH: all anchor-overlapping reads count toward DP,
        // regardless of allele classification (REF, ALT, or other/ambiguous).
        // This ensures DP reflects true physical coverage at the locus.
        // First-class records only: a supplementary segment sharing a QNAME with its
        // primary is the same physical read, and counting both would report depth 2 where
        // one read exists. It still reaches fragment evidence below.
        let is_reverse = record.is_reverse();
        if first_class {
            counts.dp += 1;
            if is_reverse {
                counts.dp_rev += 1;
            } else {
                counts.dp_fwd += 1;
            }
            // Counted with DP, not at classification: SW_FALLBACK claims the
            // row's counts came partly from a different scorer, which is only
            // true for reads that contribute to those counts.
            if result.sw_fallback {
                counts.sw_fallback_reads += 1;
            }
            // Which rule decided the read (logs only).
            let d = &mut counts.decisions;
            d.ref_withdrawn += u64::from(result.ref_withdrawn);
            d.alt_withdrawn += u64::from(result.alt_withdrawn);
            d.alt_unjudged += u64::from(result.alt_unjudged);
            d.alt_by_bases += u64::from(result.alt_by_bases);
            d.carrier_judged += u64::from(result.carrier_judged);
            d.carrier_fallback += u64::from(result.carrier_fallback);
            d.sibling_ref_excluded += u64::from(result.is_ref && ref_claimed_by_sibling);
            d.sibling_alt_claimed += u64::from(result.is_alt && claimed_by_sibling);
            d.clip_admitted += u64::from(clip_admitted);
        }

        // ── FRAGMENT TRACKING: track ALL fragments for DPF.
        // FragmentEvidence::observe() correctly handles (false, false) —
        // it skips updating best_ref_qual/best_alt_qual but still tracks
        // the fragment for DPF in the downstream resolution loop.
        // UMI-aware fragment grouping: when umi_tag is set, reads with
        // different UMIs are treated as distinct molecules. The UMI is
        // extracted from the BAM aux tag (e.g., RX:Z:ACGT).
        let (mol_hash, umi_tagged) = molecule_key(record, umi_tag, amplicon_mode);
        if umi_tagged {
            counts.umi_tagged_reads += 1;
        }
        let is_read1 = record.is_first_in_template();
        let is_forward = !is_reverse;

        let evidence = fragments.entry(mol_hash).or_insert_with(FragmentEvidence::new);

        // mFSD: compute physical fragment size from CIGAR, correcting TLEN for indels.
        // Formula: physical = |TLEN| - D + I (validated on real MSK-ACCESS BAMs).
        // observe() stores min(R1, R2) to handle cases where only one read
        // spans the indel (defensive for WGS/WES; no-op for cfDNA overlap).
        // is_n_base: uses the explicit has_n_base flag from variant classification
        // (set by check_snp/check_mnp/check_complex when N detected at a
        // discriminating position) rather than the previous heuristic
        // (base_qual==0 && !is_ref && !is_alt) which could mis-classify
        // true third-allele reads with qual=0 as N-class fragments.
        let tlen = mfsd::calc_physical_insert_size(record);
        let informative = !result.uninformative || ref_claimed_by_sibling;
        evidence.observe(is_ref, is_alt, base_qual, is_read1, is_forward, tlen, result.has_n_base, result.is_structural, record.mapq(), informative);

        // Secondary/supplementary records end here: they are fragment evidence
        // only. Every counter below (n_count, any_alt/partial_alt, RD/AD,
        // sense/antisense depth, splice-spanning) is read-level and defined over
        // the same first-class read set as DP — counting these records there
        // broke DP >= RD+AD whenever the secondary or supplementary filter was
        // disabled.
        if !first_class {
            continue;
        }

        // ── ALLELE-SPECIFIC COUNTS: only REF/ALT reads contribute to RD/AD.
        // DP and DPF are already recorded above.
        //
        // N-base counting (diagnostic):
        // Reads with N at ≥1 discriminating position are counted separately
        // for duplex masking QC. This is independent of allele classification —
        // a read classified as ALT via masked MNP evaluation can still have
        // had N at one masked position.
        if result.has_n_base {
            counts.n_count += 1;
            trace!("n_count++: read has N at discriminating position (total={})", counts.n_count);
        }
        //
        // Decomposed counting (any_alt / partial_alt):
        // - Full ALT match: ad++, any_alt++ (invariant: any_alt = ad + partial_alt)
        // - Partial ALT match (some discriminating positions match ALT): any_alt++, partial_alt++
        // - Nearby evidence (right-length INDEL, close alignment score): any_alt++, partial_alt++
        // - Neither/REF with no evidence: no any_alt/partial_alt change
        //
        // Note: has_nearby_evidence propagates structural evidence from variant checkers
        // and alignment backends. This captures reads with right-length INDELs but wrong
        // sequences (e.g., PAX5 A>CCC) that were previously lost as silent REF calls.
        if !is_ref && !is_alt {
            // Check for partial ALT evidence before skipping. A read whose
            // ALT or REF call was claimed by a sibling is partial evidence
            // for this row: the molecule carries a variant in this tract but
            // belongs to the sibling's representation.
            let sibling_claimed = claimed_by_sibling || ref_claimed_by_sibling;
            if result.partial_match_count > 0 || result.has_nearby_evidence || sibling_claimed {
                counts.any_alt += 1;
                counts.partial_alt += 1;
                trace!("partial_alt++: partial_match={} nearby_evidence={} sibling_claimed={} (any_alt={}, partial_alt={})",
                    result.partial_match_count, result.has_nearby_evidence, sibling_claimed,
                    counts.any_alt, counts.partial_alt);
            }
            continue;
        }

        // is_ref with nearby evidence: the read is classified as REF, but
        // the checker found structural evidence of the variant (e.g., right-length
        // INDEL with wrong sequence). Count as partial_alt to enable the
        // PARTIAL_DOMINANT diagnostic flag. The read still counts as rd++.
        if is_ref && result.has_nearby_evidence {
            counts.any_alt += 1;
            counts.partial_alt += 1;
            trace!("partial_alt++ (nearby evidence on REF read): any_alt={}, partial_alt={}",
                counts.any_alt, counts.partial_alt);
        }

        if is_ref {
            counts.rd += 1;
            if is_reverse { counts.rd_rev += 1; } else { counts.rd_fwd += 1; }
        } else if is_alt {
            counts.ad += 1;
            counts.any_alt += 1; // Full ALT → counts toward any_alt
            if result.mnp_confirmed { counts.mnp_confirmed_alt += 1; }
            if is_reverse { counts.ad_rev += 1; } else { counts.ad_fwd += 1; }

            // ── RNA-SPECIFIC ALT TRACKING ──
            // Splice-spanning count: ALT reads that cross a splice junction
            if mode == "rna" && rna::has_splice_junction(record) {
                counts.splice_spanning_count += 1;
            }
        }

        // ── RNA SENSE/ANTISENSE DEPTH: track strand-specific depth.
        // Uses the same dUTP logic as is_sense_strand to classify reads.
        if mode == "rna" {
            if rna::is_sense_strand(record, variant.gene_strand, strandedness) {
                counts.sense_depth += 1;
                if is_alt { counts.sense_strand_alt_count += 1; }
            } else {
                counts.antisense_depth += 1;
                if is_alt { counts.antisense_strand_alt_count += 1; }
            }
        }
    }

    // ── QC MEDIAN COMPUTATION: compute median distance-to-end for REF/ALT.
    counts.alt_dist_end_median = compute_median_u32(&mut alt_dists);
    counts.ref_dist_end_median = compute_median_u32(&mut ref_dists);

    // ── RNA editing site flag (DB-only) ────────────────────────────
    // Flag is True ONLY when a REDIportal database is provided AND this
    // variant's position is found in the database. No DB = no flagging.
    // This avoids false positives from pattern-matching-only heuristics
    // (e.g. every A→G SNP being flagged as a potential editing site).
    if mode == "rna" {
        counts.rna_editing_site_overlap = editing_sites.as_ref().is_some_and(|sites| {
            // Flag only when the variant's substitution matches the catalogued edit
            // (genomic-forward Ref>Ed, e.g. A>G on '+' or T>C on '-'). Position alone
            // would over-flag unrelated substitutions sitting on an editing coordinate.
            let found = variant.ref_allele.len() == 1
                && variant.alt_allele.len() == 1
                && {
                    let chrom = crate::shared::contig::normalize_contig(&variant.chrom);
                    sites.contains(&(
                        chrom,
                        variant.pos,
                        variant.ref_allele.as_bytes()[0].to_ascii_uppercase(),
                        variant.alt_allele.as_bytes()[0].to_ascii_uppercase(),
                    ))
                };
            if found {
                trace!(
                    "editing site match at {}:{} {}>{}",
                    variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele,
                );
            }
            found
        });
    }

    // ── FRAGMENT RESOLUTION: Resolve fragment-level counts using quality-weighted
    // consensus. Each fragment contributes exactly ONE allele call (REF xor ALT),
    // preventing the double-counting bug where R1=REF + R2=ALT inflated both
    // rdf and adf.
    //
    // Strand bias uses allele-specific orientation: the strand of the read
    // that provided the best evidence for the winning allele, not just R1.
    // Example: R1=Fwd/REF(Q10) + R2=Rev/ALT(Q30) → ALT wins, counted as
    // adf_rev (not adf_fwd).
    //
    // mFSD size vectors: one per Krewlyzer fragment class.
    // Populated below during resolution. Only sizes in the cfDNA-valid range
    // (50–1000 bp) with a known TLEN are added. GC correction is not applied —
    // GC bias affects count depth, not fragment length, so these raw sizes are
    // already unbiased samples of the true size distribution.
    //
    // mFSD is output-aware. When `--mfsd` is off the size vectors are never
    // reserved or filled and the stats block below is skipped — sparing the per-variant
    // `counts.ref_sizes`/`alt_sizes` arrays that are otherwise held on every BaseCounts
    // for the whole run (the dominant mFSD memory cost under Nextflow fan-out).
    let mfsd_cap = if mfsd { fragments.len() } else { 0 };
    let mut ref_sizes:    Vec<f64> = Vec::with_capacity(mfsd_cap);
    let mut alt_sizes:    Vec<f64> = Vec::with_capacity(mfsd_cap);
    let mut nonref_sizes: Vec<f64> = Vec::with_capacity(mfsd_cap);
    let mut n_sizes:      Vec<f64> = Vec::with_capacity(mfsd_cap);

    // Observation export: one row per fragment, mirroring the mFSD gate above —
    // capacity 0 (no allocation) when the caller did not ask for observations.
    let mut observations: Vec<Observation> =
        Vec::with_capacity(if emit_obs { fragments.len() } else { 0 });

    for (molecule_hash, evidence) in fragments.iter() {
        let (frag_ref, frag_alt) = evidence.resolve(qual_diff_threshold);
        // One class per molecule, from the same resolved call that feeds rdf/adf/dpf
        // below, read by both the export and mFSD, so neither can diverge from the counts.
        let class = evidence.class(frag_ref, frag_alt);

        // Per-molecule export: every molecule is a row (the rows reconcile with DPF), so a
        // molecule with no readable allele is OTHER here.
        if emit_obs {
            let (allele, best_qual) = match class {
                MoleculeClass::Ref => (OBS_ALLELE_REF, evidence.best_ref_qual),
                MoleculeClass::Alt => (OBS_ALLELE_ALT, evidence.best_alt_qual),
                MoleculeClass::N => (OBS_ALLELE_N, 0),
                MoleculeClass::Other | MoleculeClass::Unread => (OBS_ALLELE_OTHER, 0),
            };
            observations.push(Observation {
                variant_index,
                molecule_hash: *molecule_hash,
                allele,
                best_qual,
                min_mapq: evidence.min_mapq,
            });
        }

        // Count every fragment in dpf regardless of consensus outcome.
        // Discarded fragments (ambiguous R1-vs-R2 within quality threshold)
        // are still real molecules — tracking them in dpf makes the gap
        // dpf - (rdf + adf) a useful quality metric for the locus.
        counts.dpf += 1;

        if frag_ref {
            counts.rdf += 1;
            // Use REF-specific orientation (strand of best REF evidence)
            if let Some(ori) = evidence.ref_orientation() {
                if ori { counts.rdf_fwd += 1; } else { counts.rdf_rev += 1; }
            }
        } else if frag_alt {
            counts.adf += 1;
            // Use ALT-specific orientation (strand of best ALT evidence)
            if let Some(ori) = evidence.alt_orientation() {
                if ori { counts.adf_fwd += 1; } else { counts.adf_rev += 1; }
            }
        }

        // mFSD: classify fragment into size class vectors (only when --mfsd is on;
        // an output-aware gate — leaves the size vectors empty otherwise)
        if let Some(sz) = evidence.insert_size {
            if mfsd && (50..=1000).contains(&sz) {
                let sz_f = sz as f64;
                match class {
                    MoleculeClass::Ref => ref_sizes.push(sz_f),
                    MoleculeClass::Alt => alt_sizes.push(sz_f),
                    MoleculeClass::N => n_sizes.push(sz_f),
                    MoleculeClass::Other => nonref_sizes.push(sz_f),
                    // No readable allele, so no size class.
                    MoleculeClass::Unread => {}
                }
            }
        }
    }

    // ── STRAND BIAS
    let (sb_pval, sb_or) =
        fisher_strand_bias(counts.rd_fwd, counts.rd_rev, counts.ad_fwd, counts.ad_rev);
    counts.sb_pval = sb_pval;
    counts.sb_or = sb_or;

    let (fsb_pval, fsb_or) = fisher_strand_bias(
        counts.rdf_fwd, counts.rdf_rev, counts.adf_fwd, counts.adf_rev,
    );
    counts.fsb_pval = fsb_pval;
    counts.fsb_or = fsb_or;

    // ── mFSD Statistics — only when requested (output-aware gate). When off,
    // the size vectors are empty and all mFSD fields stay at their BaseCounts default;
    // the writers omit the mFSD columns and the post-counting BH-FDR pass skips them.
    if mfsd {
        compute_mfsd_stats(&mut counts, ref_sizes, alt_sizes, nonref_sizes, n_sizes, variant);
    }

    warn_sw_fallback(variant, counts.sw_fallback_reads);
    let d = counts.decisions;
    warn_unjudged(variant, d.alt_unjudged, d.carrier_fallback);

    // Per-phase classification breakdown, and which rule decided the depth reads
    debug!(
        "Phase stats {}:{} {}>{}: P0={} P1={} P2={} P2.5={} P3={} splice_skip_excluded={} \
         sw_fallback={} clip_candidates={} | withdrawn ref={} alt={} (unjudged {}) \
         alt_by_bases={} carrier judged={} fallback={} sibling ref_excluded={} alt_claimed={} \
         clip_admitted={} ({} backend, {} reads overlapping the window)",
        variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele,
        phase_counts[0], phase_counts[1], phase_counts[2], phase_counts[3], phase_counts[4],
        counts.splice_skip_excluded, counts.sw_fallback_reads, counts.clip_candidates,
        d.ref_withdrawn, d.alt_withdrawn, d.alt_unjudged, d.alt_by_bases, d.carrier_judged,
        d.carrier_fallback, d.sibling_ref_excluded, d.sibling_alt_claimed, d.clip_admitted,
        match backend {
            AlignmentBackend::SmithWaterman => "SW",
            AlignmentBackend::PairHMM { .. } => "HMM",
        },
        reads_considered,
    );

    Ok((counts, observations))
}


/// Build the REF and ALT haplotypes of `v` restricted to the genomic
/// window `[w_lo, w_hi)`. Requires the variant's (genomic) ref_context to
/// cover the window and the window to contain the full event span —
/// otherwise the comparison would be lopsided and the caller keeps the
/// read's classification.
fn window_haplotypes(v: &Variant, w_lo: i64, w_hi: i64) -> Option<(Vec<u8>, Vec<u8>)> {
    let ctx = v.ref_context.as_ref()?.as_bytes();
    let cs = v.ref_context_start;
    let ce = cs + ctx.len() as i64;
    let span_end = v.pos + v.ref_allele.len() as i64;
    if w_lo < cs || w_hi > ce || w_lo > v.pos || w_hi < span_end {
        return None;
    }
    let left = &ctx[(w_lo - cs) as usize..(v.pos - cs) as usize];
    let right = &ctx[(span_end - cs) as usize..(w_hi - cs) as usize];
    let mut ref_hap = Vec::with_capacity(left.len() + v.ref_allele.len() + right.len());
    ref_hap.extend_from_slice(left);
    ref_hap.extend_from_slice(v.ref_allele.as_bytes());
    ref_hap.extend_from_slice(right);
    let mut alt_hap = Vec::with_capacity(left.len() + v.alt_allele.len() + right.len());
    alt_hap.extend_from_slice(left);
    alt_hap.extend_from_slice(v.alt_allele.as_bytes());
    alt_hap.extend_from_slice(right);
    Some((ref_hap, alt_hap))
}

/// Mask sub-threshold and N bases to a sentinel byte that matches no
/// haplotype base: the contest then charges them equally against every
/// candidate instead of letting sequencing noise coincidentally vote for
/// one — the same bases the BQ-aware classification masked (cross-backend
/// quality contract).
fn mask_low_qual(seq: &mut [u8], quals: &[u8], min_baseq: u8) {
    for (b, &q) in seq.iter_mut().zip(quals.iter()) {
        if q < min_baseq || *b == b'N' || *b == b'n' {
            *b = 0;
        }
    }
}

/// Whether a read classified REF for `variant` (or whose REF call was
/// withdrawn as uninformative) is ALT for a co-annotated sibling whose change
/// lies inside `variant`'s discrimination window
/// (`window_siblings`, from `window::siblings_in_window`). Such a read
/// carries a different allele where this row's REF is read, so it is not REF
/// testimony here. A carrier of a sibling elsewhere in the group shows the
/// reference across every base that could tell this row's alleles apart, and
/// stays REF (as IGV shows it, and as GATK counts REF at a site unless an
/// event overlaps it). Callers apply this before fragment evidence, so REF
/// reads and REF fragments exclude the same molecules.
#[allow(clippy::too_many_arguments)]
fn sibling_claims_ref<F: Fn(u8, u8) -> i32>(
    record: &Record,
    variant: &Variant,
    result: &ClassifyResult,
    window_siblings: &[&Variant],
    quals: &[u8],
    min_baseq: u8,
    alt_aligner: &mut Aligner<F>,
    ref_aligner: &mut Aligner<F>,
    backend: &AlignmentBackend,
    rules: &carrier::ReadRules,
) -> bool {
    if !result.is_ref && !result.uninformative {
        return false;
    }
    for sib in window_siblings {
        let sib_result = check_allele_with_qual(
            record, sib, &[], quals, min_baseq, alt_aligner, ref_aligner, backend, rules,
        );
        if sib_result.is_alt {
            trace!(
                "Multi-allelic guard: read={} is ALT for sibling {}>{} at {}:{} inside the \
                 window of {}>{} — excluded from its REF",
                read_name(record), sib.ref_allele, sib.alt_allele, sib.chrom, sib.pos + 1,
                variant.ref_allele, variant.alt_allele,
            );
            return true;
        }
    }
    false
}

/// Multi-allelic AD-claiming guard: decide whether a representation-tolerant
/// ALT classification really belongs to this row at a co-annotated locus.
///
/// Representation-tolerant matching (windowed matches at shifted
/// placements, BQ-masked comparison, PairHMM alignment) lets one physical
/// event satisfy several co-annotated rows in the same tract — and lets
/// unannotated same-tract ladder events be absorbed as ALT — so without
/// this guard the per-locus AD sum exceeds the number of distinct ALT
/// molecules (measured 2-5x at a real hypermutation cluster; sign-out uses
/// exclusive assignment).
///
/// Anchor-exact evidence (Phase 0: a structural CIGAR op at the annotated
/// left-aligned position, or a direct SNP base observation) is never
/// contested — a molecule genuinely carrying two anchor-exact ops counts
/// full AD on both rows. Everything else faces three demotion tests, each
/// decisive for a distinct failure mode measured on signed-out data:
///
/// 1. **REF test** (foreign events): over the read-covered slice of this
///    variant's genomic context, the read's reconstruction must explain
///    strictly better with the row's ALT haplotype than its REF haplotype
///    (Levenshtein). A read carrying only an event in the flanks ties or
///    favors REF → demote.
/// 2. **Probabilistic pure indel** (unannotated ladders): an
///    Alignment-phase (Phase 3) ALT on a *pure* indel row carries no
///    matching structural op at any position — in a contested tract that
///    is ambiguity (PairHMM absorbs D11 into a D14 row because 3 edits
///    beat 11), so demote. Complex/MNP rows are exempt: the exact-carrier
///    rule judges their carriers by their own bases.
/// 3. **Strictly-better sibling** (same-tract competitors): over a window
///    covering both spans, a sibling whose ALT haplotype explains the
///    read at strictly lower cost claims it. Equal cost means equivalent
///    representations of the same event (e.g. a delins double-annotated
///    as an insertion) — both rows keep the read rather than zeroing one.
///
/// Demoted reads surface as partial evidence (any_alt/partial_alt) and are
/// excluded from AD and ADF (the caller downgrades before fragment
/// evidence). Reads that do not fully span the event, or splice-poisoned
/// windows, keep their classification. Isolated variants are untouched.
fn sibling_claims_alt(
    record: &Record,
    variant: &Variant,
    result: &ClassifyResult,
    sibling_variants: &[Variant],
    quals: &[u8],
    min_baseq: u8,
) -> bool {
    if !result.is_alt || sibling_variants.is_empty() {
        return false;
    }
    if result.phase == ClassifyPhase::Structural {
        return false;
    }

    let read_lo = record.pos();
    let read_hi = ref_end(record);
    let own_ctx_end = variant.ref_context_start
        + variant.ref_context.as_ref().map_or(0, |c| c.len() as i64);
    let w_lo = read_lo.max(variant.ref_context_start);
    let w_hi = read_hi.min(own_ctx_end);

    // Test 1: the read's window must favor ALT strictly over REF.
    let (ref_hap, alt_hap) = match window_haplotypes(variant, w_lo, w_hi) {
        Some(h) => h,
        None => return false,
    };
    let mut recon = reconstruct_span(record, quals, w_lo, w_hi);
    if recon.splice_skip || recon.seq.is_empty() {
        return false;
    }
    let recon_quals = std::mem::take(&mut recon.quals);
    mask_low_qual(&mut recon.seq, &recon_quals, min_baseq);
    let cost_alt = levenshtein(&recon.seq, &alt_hap);
    let cost_ref = levenshtein(&recon.seq, &ref_hap);
    if cost_alt >= cost_ref {
        trace!(
            "AD-claiming guard: read={} ALT for {}>{} at {}:{} not favored over REF \
             (ALT cost {} vs REF cost {}) — partial_alt, not ad",
            read_name(record), variant.ref_allele, variant.alt_allele,
            variant.chrom, variant.pos + 1, cost_alt, cost_ref,
        );
        return true;
    }

    // Test 2: probabilistic call on a pure indel row. Known tradeoff: a true
    // carrier whose only evidence is soft-clipped (no I/D op anywhere) is
    // also demoted here — measured to be a rare sensitivity tail (clip
    // survey: 12/12 ITD loci I-op-dominant), it stays visible as
    // partial_alt, and clip rescue is tracked separately (CLIP_CANDIDATES).
    // Pure = anchor-preserved
    // deletion/insertion; a delins that merely has a 1-base side (e.g.
    // CAG>T) is complex and exempt (the exact-carrier rule judges its carriers by
    // their own bases).
    if result.phase == ClassifyPhase::Alignment && window::is_pure_indel(&variant.ref_allele, &variant.alt_allele) {
        trace!(
            "AD-claiming guard: read={} alignment-phase ALT on pure indel {}>{} at {}:{} \
             in a co-annotated cluster (no structural op) — partial_alt, not ad",
            read_name(record), variant.ref_allele, variant.alt_allele,
            variant.chrom, variant.pos + 1,
        );
        return true;
    }

    // Test 3: a sibling that explains the read strictly better claims it. One
    // that explains it exactly, as this row's ALT does (the two alleles differ
    // only at bases the read has masked), leaves it ambiguous: it is either
    // allele, so it is neither row's AD.
    for sib in sibling_variants {
        let sib_ctx_end = sib.ref_context_start
            + sib.ref_context.as_ref().map_or(0, |c| c.len() as i64);
        let s_lo = w_lo.max(sib.ref_context_start);
        let s_hi = w_hi.min(sib_ctx_end);
        let own_hap_s = match window_haplotypes(variant, s_lo, s_hi) {
            Some((_, a)) => a,
            None => continue,
        };
        let sib_hap_s = match window_haplotypes(sib, s_lo, s_hi) {
            Some((_, a)) => a,
            None => continue,
        };
        let recon_s = if (s_lo, s_hi) == (w_lo, w_hi) {
            recon.seq.clone()
        } else {
            let mut r = reconstruct_span(record, quals, s_lo, s_hi);
            if r.splice_skip || r.seq.is_empty() {
                continue;
            }
            let rq = std::mem::take(&mut r.quals);
            mask_low_qual(&mut r.seq, &rq, min_baseq);
            r.seq
        };
        let own_c = levenshtein(&recon_s, &own_hap_s);
        let sib_c = levenshtein(&recon_s, &sib_hap_s);
        if sib_c < own_c || (equal_but_masked(&recon_s, &own_hap_s) && equal_but_masked(&recon_s, &sib_hap_s)) {
            trace!(
                "AD-claiming guard: read={} ALT for {}>{} at {}:{} (cost {}) claimed by \
                 sibling {}>{} at {}:{} (cost {}) — partial_alt, not ad",
                read_name(record), variant.ref_allele, variant.alt_allele,
                variant.chrom, variant.pos + 1, own_c,
                sib.ref_allele, sib.alt_allele,
                sib.chrom, sib.pos + 1, sib_c,
            );
            return true;
        }
    }
    false
}

/// Whether a reconstructed read equals a haplotype base for base, its masked
/// bases (0, from `mask_low_qual`) matching anything.
fn equal_but_masked(read: &[u8], hap: &[u8]) -> bool {
    read.len() == hap.len() && read.iter().zip(hap).all(|(&r, &h)| r == 0 || r.eq_ignore_ascii_case(&h))
}

/// Minimum soft-clip length for an insertion clip candidate: shorter clips
/// are routine adapter/quality trimming, not unaligned inserted sequence.
const CLIP_CANDIDATE_MIN_LEN: u32 = 8;

/// Slack (bp) added to the insert length when bounding where a
/// clip-represented tandem-duplication carrier's clip boundary can land
/// relative to the anchor.
const CLIP_REACH_SLACK: i64 = 10;

/// Insertion loci only: count a first-class read carrying a clip candidate.
/// Shared by the read loops and called right after
/// classification — deliberately before the anchor-overlap gate, because a
/// clip-represented tandem-duplication carrier can align entirely past the
/// anchor (its clip covers the inserted copy) and still be the evidence the
/// CLIP_CANDIDATES flag points at.
fn tally_clip_candidate(
    counts: &mut BaseCounts,
    record: &Record,
    variant: &Variant,
    first_class: bool,
) {
    let ins_len = variant.alt_allele.len() as i64 - variant.ref_allele.len() as i64;
    if first_class && ins_len > 0 {
        let reach = ins_len + CLIP_REACH_SLACK;
        if has_clip_boundary_in(record, variant.pos - reach, variant.pos + reach) {
            counts.clip_candidates += 1;
        }
    }
}

/// Whether the read has a soft clip of at least `CLIP_CANDIDATE_MIN_LEN`
/// whose boundary (the reference position where the clipped bases would
/// continue the alignment) lies in `[lo, hi]`. Hard clips at the read ends
/// are skipped when locating the soft clip.
fn has_clip_boundary_in(record: &Record, lo: i64, hi: i64) -> bool {
    let (lead, tail) = soft_clips(record);
    (lead >= CLIP_CANDIDATE_MIN_LEN && (lo..=hi).contains(&record.pos()))
        || (tail >= CLIP_CANDIDATE_MIN_LEN && (lo..=hi).contains(&ref_end(record)))
}

/// Once per counting pass: rows the rules can only judge degraded. Prep rejects
/// or prepares every such row, so these come from callers passing unprepared
/// variants to the engine directly.
fn warn_degraded_variants(variants: &[Variant]) {
    let empty = |v: &&Variant| v.ref_allele.is_empty() || v.alt_allele.is_empty();
    let unprepared = variants
        .iter()
        .filter(|v| !empty(v) && v.ref_allele.len() != v.alt_allele.len() && v.ref_context.is_none())
        .count();
    if unprepared > 0 {
        warn!(
            "{} insertion, deletion or complex variant(s) carry no prepared reference context \
             (not prepared against a FASTA, or prep failed): the rules that read the reference \
             around the event run degraded for them",
            unprepared,
        );
    }
    let degenerate = variants
        .iter()
        .filter(|v| v.ref_allele.len() > 1 && v.ref_allele.eq_ignore_ascii_case(&v.alt_allele))
        .count();
    if degenerate > 0 {
        warn!(
            "{} MNP row(s) have REF equal to ALT: no read can show either allele, so they count \
             depth only",
            degenerate,
        );
    }
    let empties = variants.iter().filter(empty).count();
    if empties > 0 {
        warn!(
            "{} row(s) have an empty allele: no read can show it, so they count depth only",
            empties,
        );
    }
}

/// Once per variant: reads a rule could not judge for want of reference around
/// the event, so their loss is never silent.
fn warn_unjudged(variant: &Variant, alt_unjudged: u64, carrier_fallback: u64) {
    if alt_unjudged > 0 {
        warn!(
            "{}:{} {}>{}: {} ALT read(s) span neither informative window and could not be judged \
             by their bases: no prepared reference holds the event (the variant was not \
             prepared against a FASTA, or failed prep), or the base that would decide lies \
             past the prepared reference's end (a contig edge, or a short fetch) — counted \
             as depth only",
            variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele, alt_unjudged,
        );
    }
    if carrier_fallback > 0 {
        warn!(
            "{}:{} {}>{}: {} read(s) could not be judged by the exact-carrier rule: {} — \
             classified by the previous complex classifier",
            variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele, carrier_fallback,
            carrier::unjudged_reason(variant),
        );
    }
}

/// One WARN per variant when depth reads could not be evaluated by the
/// PairHMM backend's pangenomic haplotype matrix. The Smith-Waterman fallback
/// is kept (operator decision) but must not be silent: it only fires when the
/// variant's reference context is missing or does not contain it, i.e.
/// upstream input was malformed. Without any context no scorer can run, so
/// those reads end NEITHER — the message says which outcome applied.
fn warn_sw_fallback(variant: &Variant, n: u32) {
    if n > 0 {
        let outcome = if variant.ref_context.is_none() {
            "no scorer can run without a reference context, so they were left NEITHER"
        } else {
            "they were scored by the Smith-Waterman fallback instead (NEITHER where SW \
             could not run either)"
        };
        warn!(
            "{}:{} {}>{}: {} read(s) could not be evaluated by the pangenomic haplotype \
             matrix ({}); {} — flagged SW_FALLBACK({})",
            variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele,
            n, super::pangenome::matrix_failure_reason(variant), outcome, n,
        );
    }
}

/// The exon edges and junction ends near a variant, ascending: the annotated
/// intron boundaries and the junctions the reads over its window splice at (an
/// unannotated junction the data show is an edge too).
fn clip_edges(variant: &Variant, read_cache: &[Record], annotation: &Option<std::sync::Arc<AnnotationIndex>>) -> Vec<i64> {
    // Annotated boundaries this far from the read window cannot meet a clip.
    const REACH: i64 = 1000;
    let (v_start, v_end) = window::read_window(variant);
    let mut edges = annotation.as_ref().map_or_else(Vec::new, |a| {
        let chrom = crate::shared::contig::normalize_contig(&variant.chrom);
        a.intron_boundaries_in(&chrom, v_start - REACH, v_end + REACH)
    });
    for record in read_cache {
        if record.pos() >= v_end || ref_end(record) <= v_start {
            continue;
        }
        for (n0, n1) in rna::extract_splice_junctions(record) {
            edges.push(n0);
            edges.push(n1);
        }
    }
    edges.sort_unstable();
    edges.dedup();
    edges
}

/// Whether the mapping rule admits a read to the counts: MAPQ at least
/// `min_mapq`, or in RNA a unique mapper (`NH:i:1`) below it
/// (`rna::is_valid_rna_alignment`).
fn mapping_admits(record: &Record, mode: &str, min_mapq: u8) -> bool {
    if mode == "rna" {
        rna::is_valid_rna_alignment(record, min_mapq)
    } else {
        record.mapq() >= min_mapq
    }
}

/// Whether a read is one the read-level counts read: a first-class record (not
/// secondary or supplementary) the mapping rule admits and, under RNA
/// strandedness enforcement, a sense read.
fn counted_read(
    record: &Record,
    variant: &Variant,
    mode: &str,
    enforce_strandedness: bool,
    strandedness: rna::Strandedness,
    min_mapq: u8,
) -> bool {
    let first_class = !(record.is_secondary() || record.is_supplementary());
    let antisense_excluded =
        mode == "rna" && enforce_strandedness && !rna::is_sense_strand(record, variant.gene_strand, strandedness);
    first_class && !antisense_excluded && mapping_admits(record, mode, min_mapq)
}

/// A REF call on an insertion or deletion stands only when the read is
/// informative (`window::read_is_informative`). A read that starts or ends
/// inside the repeat tract matches both alleles (the aligner places no gap
/// either way), so it counts toward depth and fragment depth but is neither
/// REF nor ALT. The anchor quality and any nearby-indel evidence are kept.
fn ref_needs_the_window(record: &Record, variant: &Variant, mut result: ClassifyResult) -> ClassifyResult {
    if result.is_ref && !window::read_is_informative(record, variant) {
        trace!(
            "{}:{} {}>{} read={}: no aligned block of {:?} spans an informative window {:?} (0-based) \
             — uninformative, not REF",
            variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele, read_name(record),
            crate::shared::bam_utils::aligned_blocks(record), window::informative_windows(variant),
        );
        result.is_ref = false;
        result.is_structural = false;
        result.uninformative = true;
        result.ref_withdrawn = true;
    }
    result
}

/// An ALT call on an insertion or deletion from a read that spans neither of the
/// informative windows read on the ALT haplotype (`window::alt_read_is_informative`)
/// stands only when the read's own bases tell the alleles apart
/// (`window::alt_bases_discriminate`). The CIGAR's gap alone is placement: a
/// carrier ending inside the repeat holds only bases both alleles share, so it
/// counts toward depth and fragment depth but is neither REF nor ALT (the REF
/// side's rule, on the ALT side). A read whose bases discriminate (a truncated
/// long insertion's carrier, whose insert consumes the read's span) keeps its ALT.
fn alt_needs_the_window(
    record: &Record,
    variant: &Variant,
    quals: &[u8],
    min_baseq: u8,
    mut result: ClassifyResult,
) -> ClassifyResult {
    if !result.is_alt {
        return result;
    }
    let kept = |why: &str| {
        trace!(
            "{}:{} {}>{} read={}: ALT kept — {}",
            variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele, read_name(record), why,
        );
    };
    if window::alt_read_is_informative(record, variant, quals, min_baseq) {
        kept("it spans an ALT-side window");
        return result;
    }
    let judged = window::alt_bases_discriminate(record, variant, quals, min_baseq);
    if judged == Some(true) {
        kept("its own bases tell the alleles apart");
        result.alt_by_bases = true;
        return result;
    }
    let why = match judged {
        Some(_) => "its bases fit both alleles",
        None if variant.event_ref.is_none() && variant.ref_context.is_none() => {
            "no prepared reference to read its bases against"
        }
        None => "the prepared reference ends before the base that would decide",
    };
    trace!(
        "{}:{} {}>{} read={}: no aligned block of {:?} spans an ALT-side window around {:?} (0-based) \
         and {} — uninformative, not ALT",
        variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele, read_name(record),
        crate::shared::bam_utils::aligned_blocks(record), window::change_interval(variant), why,
    );
    result.is_alt = false;
    result.is_structural = false;
    result.has_nearby_evidence = false;
    result.partial_match_count = 0;
    result.uninformative = true;
    result.alt_withdrawn = true;
    result.alt_unjudged = judged.is_none();
    result
}


/// Check if a read supports the reference or alternate allele.
/// Returns `ClassifyResult` containing (is_ref, is_alt, base_quality, phase)
/// where base_quality is the quality score at the variant position
/// (used for fragment consensus) and phase indicates which classification
/// stage resolved the read.
///
/// Each variant-type handler returns quality directly from its own CIGAR
/// walk, ensuring correct quality extraction even for reads carrying
/// indels at the variant position.
///
/// The `alt_aligner` and `ref_aligner` are reusable SW aligners created
/// once per variant in the read loop and threaded through to avoid per-read
/// allocation (indelpost pattern).
#[allow(clippy::too_many_arguments)]
fn check_allele_with_qual<F: Fn(u8, u8) -> i32>(
    record: &Record,
    variant: &Variant,
    siblings: &[Variant],
    quals: &[u8],
    min_baseq: u8,
    alt_aligner: &mut Aligner<F>,
    ref_aligner: &mut Aligner<F>,
    backend: &AlignmentBackend,
    rules: &carrier::ReadRules,
) -> ClassifyResult {
    // Splice-skip triage first: a read whose CIGAR N spans every
    // discriminating position observes nothing at this locus and must not
    // testify for either allele (or count depth — the engine skips
    // covers_locus=false reads entirely). Checked before dispatch so no
    // checker's anchor-based fast path can read an asserted splice as REF.
    // Reads without N ops return None immediately, so DNA classification
    // pays one CIGAR-flag scan and nothing else.
    if let Some(result) = splice_skip_triage(record, variant) {
        return result;
    }

    // Dispatch on the alleles (`window::allele_kind`), never on the variant_type
    // label, which callers write inconsistently (e.g. "COMPLEX" for what is a pure
    // deletion after normalization).
    trace!(
        "check_allele {}:{} {}>{} read={}",
        variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele, read_name(record),
    );
    let Some(kind) = window::allele_kind(&variant.ref_allele, &variant.alt_allele) else {
        // A row with an empty allele shows no allele (prep rejects such rows; the
        // counting pass warns once about any passed to it directly).
        return ClassifyResult::neither(ClassifyPhase::Structural);
    };

    if kind == AlleleKind::Snv {
        // SNP: single base substitution — no Phase 3 needed
        check_snp(record, variant, quals, min_baseq)
    } else if kind == AlleleKind::Mnp && carrier::indel_at_block(record, variant) {
        // An MNP read with an indel in or right beside the block (an aligner may
        // write a shifted block as an insertion before it and a deletion after):
        // it counts only if its own bases carry the whole allele (exact-carrier rule).
        let route = "an MNP read with an indel at the block";
        classify_complex(record, variant, route, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend, rules)
    } else if kind == AlleleKind::Mnp {
        // MNP: selective discriminating-position quality gate with no Phase 3 fallback.
        match check_mnp(record, variant, quals, min_baseq) {
            MnpResult::Ref(q, had_n) => {
                let mut r = ClassifyResult::is_ref(q, ClassifyPhase::MaskedCompare);
                r.has_n_base = had_n;
                r
            }
            MnpResult::Alt(q, had_n, confirmed) => {
                let mut r = ClassifyResult::is_alt(q, ClassifyPhase::MaskedCompare);
                r.has_n_base = had_n;
                r.mnp_confirmed = confirmed;
                r
            }
            MnpResult::LowQuality(partial, had_n) => {
                // After masked per-position evaluation, LowQuality means ALL
                // discriminating positions were masked (BQ < min_baseq or N).
                // Do NOT route to check_complex: PairHMM/SW is designed for
                // indel realignment, not MNP classification, and is biased
                // toward REF for multi-base substitutions.
                // C++ GBCMS (baseCountDNP) has no fallback — reads are simply
                // not counted. Match that behavior.
                // Fragment impact: observe(false, false) → DPF++ but not
                // RDF/ADF. If mate read provides evidence, mate's call wins.
                // `partial` carries positions_matching_alt for partial_alt counting.
                trace!(
                    "MNP LowQuality: all discriminating positions masked, partial_alt_positions={}, had_n={}",
                    partial, had_n
                );
                ClassifyResult::neither_with_partial(ClassifyPhase::MaskedCompare, partial, had_n)
            }
            MnpResult::ThirdAllele(partial, had_n) => {
                // Unmasked positions show mixed REF/ALT or third-allele bases.
                // `partial` carries positions_matching_alt for partial_alt counting.
                if partial > 0 {
                    trace!(
                        "MNP ThirdAllele with {} positions matching ALT (partial evidence), had_n={}",
                        partial, had_n
                    );
                }
                ClassifyResult::neither_with_partial(ClassifyPhase::MaskedCompare, partial, had_n)
            }
            MnpResult::Structural => {
                // The read carries an indel or clip at the MNP: it counts only if
                // its own bases carry the whole allele (exact-carrier rule).
                let route = "an MNP read with an indel or clip at the block";
                classify_complex(record, variant, route, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend, rules)
            }
        }
    } else if kind == AlleleKind::Insertion {
        // Pure insertion: CIGAR-based fast paths, then backend-aware Phase 3 fallback
        let result = check_insertion(record, variant, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend);
        ref_needs_the_window(record, variant, alt_needs_the_window(record, variant, quals, min_baseq, result))
    } else if kind == AlleleKind::Deletion {
        // Pure deletion: CIGAR-based fast paths, then backend-aware Phase 3 fallback
        let result = check_deletion(record, variant, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend);
        ref_needs_the_window(record, variant, alt_needs_the_window(record, variant, quals, min_baseq, result))
    } else {
        // Complex, judged by its whole allele with the exact-carrier rule
        // (`classify_complex`), reaching check_complex only when the rule cannot
        // judge the variant:
        // - a one-base REF whose ALT changes it (C>TA, A>CCC) is a delins, not an
        //   insertion of its tail: the insertion check compares only the inserted
        //   bases, so reads that keep the anchor would count ALT;
        // - a deletion whose anchor also changes (GC>T): check_deletion judges the
        //   gap and the deleted bases (ref_allele[1..]) and never reads the anchor,
        //   so it cannot tell such a carrier from a pure-deletion carrier;
        // - a delins with both alleles longer than one base.
        let route = match (variant.ref_allele.len(), variant.alt_allele.len()) {
            (1, _) => "an insertion whose anchor changes",
            (_, 1) => "a deletion whose anchor changes",
            _ => "a delins",
        };
        classify_complex(record, variant, route, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend, rules)
    }
}

/// A complex variant (delins, a deletion whose anchor also changes, an insertion
/// whose ALT changes the anchor such as C>TA, or an MNP read carrying an indel or
/// clip at the block) counts a read only when the read's own bases carry the
/// whole allele: the exact-carrier rule (`carrier`). A variant the rule cannot
/// judge (no prepared reference, or one that does not hold the REF allele or the
/// event with its flank) goes to the previous classifier (`check_complex`); such
/// reads are counted and warned once per variant. `route` names why the read
/// came here, for the trace.
#[allow(clippy::too_many_arguments)]
fn classify_complex<F: Fn(u8, u8) -> i32>(
    record: &Record,
    variant: &Variant,
    route: &str,
    siblings: &[Variant],
    quals: &[u8],
    min_baseq: u8,
    alt_aligner: &mut Aligner<F>,
    ref_aligner: &mut Aligner<F>,
    backend: &AlignmentBackend,
    rules: &carrier::ReadRules,
) -> ClassifyResult {
    if let Some(mut result) = carrier::classify(record, variant, quals, min_baseq, rules) {
        result.carrier_judged = true;
        trace!(
            "{}:{} {}>{} read={}: exact-carrier rule ({}): ref={} alt={} nearby={} uninformative={} \
             clip_admissible={} read_from_right={} mnp_confirmed={}",
            variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele, read_name(record),
            route, result.is_ref, result.is_alt, result.has_nearby_evidence, result.uninformative,
            result.clip_admissible, result.read_from_right, result.mnp_confirmed,
        );
        return result;
    }
    trace!(
        "{}:{} {}>{} read={}: the exact-carrier rule cannot judge ({}) → previous complex classifier",
        variant.chrom, variant.pos + 1, variant.ref_allele, variant.alt_allele, read_name(record), route,
    );
    let mut result = check_complex(record, variant, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend);
    result.carrier_fallback = true;
    result
}

/// Smith-Waterman scoring for the read classifiers: an N in the read or the
/// haplotype is neutral (0), evidence for neither a match nor a mismatch (as in
/// GATK), so duplex-masked bases decide nothing.
fn sw_score(a: u8, b: u8) -> i32 {
    if a == b'N' || b == b'N' { 0 } else if a == b { 1 } else { -1 }
}

/// The read's base qualities as the classifiers read them: heuristic BAQ when
/// `use_baq` (sparing the variant's own evidence, `baq_spare`), else as stored.
fn effective_quals(record: &Record, use_baq: bool, baq_spare: Option<(i64, i64)>) -> std::borrow::Cow<'_, [u8]> {
    match use_baq.then(|| apply_heuristic_baq(record, baq_spare)).flatten() {
        Some(adjusted) => std::borrow::Cow::Owned(adjusted),
        None => std::borrow::Cow::Borrowed(record.qual()),
    }
}

/// One record's identity across the bins that fetch it: name, flags, contig and
/// position (a read and its mate, or its supplementary pieces, differ in flags or
/// position).
fn record_identity(record: &Record) -> u64 {
    let mut h = std::collections::hash_map::DefaultHasher::new();
    std::hash::Hash::hash(&(record.qname(), record.flags(), record.tid(), record.pos()), &mut h);
    std::hash::Hasher::finish(&h)
}

/// The read's molecule key, and whether it carried the UMI: a hash of its QNAME,
/// with its UMI when `umi_tag` is set (reads with different UMIs are different
/// molecules). In amplicon mode R1 and R2 key apart, so each read is its own
/// observation (no fragment consensus).
fn molecule_key(record: &Record, umi_tag: Option<[u8; 2]>, amplicon_mode: bool) -> (u64, bool) {
    let umi = umi_tag.and_then(|tag| match record.aux(&tag) {
        Ok(rust_htslib::bam::record::Aux::String(s)) => Some(s.as_bytes()),
        _ => None,
    });
    let mut key = if umi_tag.is_some() {
        hash_molecule(record.qname(), umi)
    } else {
        hash_qname(record.qname())
    };
    if amplicon_mode {
        key ^= if record.is_first_in_template() { 0x1 } else { 0x2 };
    }
    (key, umi.is_some())
}

/// The rule that decided a read's call, for its `read call` trace line.
fn deciding_rule(r: &ClassifyResult) -> &'static str {
    if r.alt_unjudged {
        "alt_unjudged"
    } else if r.alt_withdrawn {
        "alt_withdrawn"
    } else if r.ref_withdrawn {
        "ref_withdrawn"
    } else if r.alt_by_bases {
        "alt_by_bases"
    } else if r.carrier_fallback {
        "carrier_fallback"
    } else if r.carrier_judged {
        "exact_carrier"
    } else if r.sw_fallback {
        "sw_fallback"
    } else {
        "checks"
    }
}

/// A read's name for trace lines (lines from bins counted in parallel interleave,
/// so each per-read line names its read).
fn read_name(record: &Record) -> std::borrow::Cow<'_, str> {
    String::from_utf8_lossy(record.qname())
}


// ══════════════════════════════════════════════════════════════════════════════
// PER-TRANSCRIPT COUNTING
// ══════════════════════════════════════════════════════════════════════════════

/// Per-transcript read and fragment counts for a single variant.
///
/// For each transcript whose exons overlap the variant, this function:
/// 1. Filters the read cache by splice-junction compatibility
/// 2. Re-invokes allele classification on compatible reads
/// 3. Applies fragment consensus per-transcript
/// 4. Formats results as `|`-separated strings
///
/// Returns `(transcript_read_counts, transcript_fragment_counts)`.
///
/// Both are empty strings when:
/// - No transcripts overlap the variant position
/// - The variant is intronic
/// - The chromosome has no annotation
///
/// # Format
///
/// Read:     `"ENST...:AD,RD,DP|ENST...:AD,RD,DP"`
/// Fragment: `"ENST...:ADF,RDF,DPF|ENST...:ADF,RDF,DPF"`
/// Transcripts are separated by `|` (VCF-INFO-safe); fields within an entry by `:`/`,`.
///
/// # Performance
///
/// Re-invokes allele classification per transcript × per read, but this is
/// acceptable given:
/// - Typical gene loci have 1–3 overlapping transcripts
/// - RNA-seq depth is modest (50–200×)
/// - No additional BAM I/O — reuses the existing read cache
///
/// `use_baq`: whether heuristic BAQ applies at this variant, as resolved by
/// `baq_applies` (the main counts' rule).
#[allow(clippy::too_many_arguments)]
fn count_per_transcript(
    read_cache: &[Record],
    variant: &Variant,
    sibling_variants: &[Variant],
    annotation: &AnnotationIndex,
    min_mapq: u8,
    min_baseq: u8,
    fragment_qual_threshold: u8,
    backend: &AlignmentBackend,
    use_baq: bool,
    umi_tag: Option<[u8; 2]>,
    enforce_strandedness: bool,
    strandedness: rna::Strandedness,
    amplicon_mode: bool,
    rules: &carrier::ReadRules,
) -> (String, String) {
    let baq_spare = if use_baq { baq_own_span(variant) } else { None };
    let ref_guard_siblings = window::siblings_in_window(variant, sibling_variants);
    // Step 1: Find overlapping transcripts
    let chrom = crate::shared::contig::normalize_contig(&variant.chrom);
    let transcript_ids = annotation.overlapping_transcripts(&chrom, variant.pos);

    if transcript_ids.is_empty() {
        return (String::new(), String::new());
    }

    trace!(
        "per-transcript: {} overlapping transcripts at {}:{} ({})",
        transcript_ids.len(), variant.chrom, variant.pos + 1,
        transcript_ids.join(", "),
    );

    // Step 2: The variant's read window (as in the main counts)
    let (v_start, v_end) = window::read_window(variant);

    // Step 3: Create aligners (same pattern as count_variant_from_cache)

    // Step 4: Per-transcript counting
    let mut read_entries: Vec<String> = Vec::with_capacity(transcript_ids.len());
    let mut frag_entries: Vec<String> = Vec::with_capacity(transcript_ids.len());

    for tx_id in &transcript_ids {
        let tx_introns = match annotation.get_transcript_introns(tx_id) {
            Some(ti) => ti,
            None => {
                // Should not happen if overlapping_transcripts returned this ID,
                // but guard against index inconsistency.
                debug!(
                    "per-transcript: transcript {} has no intron data, skipping",
                    tx_id,
                );
                continue;
            }
        };

        // Per-transcript counters
        let mut tx_ad: u32 = 0;
        let mut tx_rd: u32 = 0;
        let mut tx_dp: u32 = 0;
        let mut tx_fragments: HashMap<u64, FragmentEvidence> = HashMap::new();

        // Fresh aligners per transcript to avoid cross-contamination
        let mut alt_aligner = Aligner::new(SW_GAP_OPEN, SW_GAP_EXTEND, &sw_score);
        let mut ref_aligner = Aligner::new(SW_GAP_OPEN, SW_GAP_EXTEND, &sw_score);

        for record in read_cache {
            // ── Overlap check: does this read overlap the variant window?
            let r_start = record.pos();
            let r_end = ref_end(record);
            if r_start >= v_end || r_end <= v_start {
                continue;
            }

            // ── Standard filters (same as count_variant_from_cache)
            if record.mapq() < min_mapq && !super::rna::is_valid_rna_alignment(record, min_mapq) {
                continue;
            }

            // ── Strandedness filter
            if enforce_strandedness && !super::rna::is_sense_strand(record, variant.gene_strand, strandedness) {
                continue;
            }

            // ── Splice-junction compatibility check
            let observed_junctions = super::rna::extract_splice_junctions(record);
            if !annotation.is_read_compatible(&observed_junctions, tx_introns, 5) {
                continue; // Incompatible junctions → skip for this transcript
            }

            // ── BAQ under the main counts' rule (`use_baq`), so the spliced
            // reads at an exon edge count here as they do in the main counts.
            let quals = effective_quals(record, use_baq, baq_spare);
            let effective_quals: &[u8] = &quals;

            // ── Allele classification
            let result = check_allele_with_qual(
                record, variant, sibling_variants, effective_quals, min_baseq,
                &mut alt_aligner, &mut ref_aligner, backend, rules,
            );

            // ── Splice-skip exclusion (same as main counting): a read whose
            // N spans every discriminating position observes nothing here —
            // no tx_dp, no fragment evidence.
            if !result.covers_locus {
                continue;
            }

            // ── Anchor overlap check, as in the main counts: a read decided from
            // windows read from the right flank counts without overlapping the
            // variant position (clip admission is DNA only, and per-transcript
            // counts run only in RNA mode)
            let overlaps_anchor = r_start <= variant.pos && r_end > variant.pos;
            if !(overlaps_anchor || result.read_from_right) {
                continue;
            }

            // ── Count reads (first-class only, same rule as the main engine:
            // secondary/supplementary records feed fragment evidence but never
            // read-level tx_dp/tx_rd/tx_ad)
            let first_class = !(record.is_supplementary() || record.is_secondary());
            if first_class {
                tx_dp += 1;
            }

            // ── Fragment tracking
            let (mol_hash, _) = molecule_key(record, umi_tag, amplicon_mode);
            let is_read1 = record.is_first_in_template();
            let is_forward = !record.is_reverse();
            let tlen = mfsd::calc_physical_insert_size(record);

            // ── Multi-allelic guards (the main engine's read-level rules): an
            // ALT match won by a sibling is excluded from tx_ad and from
            // ALT fragment evidence; a REF-classified read that is ALT for a
            // sibling whose change lies inside this variant's discrimination
            // window is excluded from tx_rd and from REF fragment evidence.
            // Both still count tx_dp.
            let claimed_by_sibling = sibling_claims_alt(
                record, variant, &result, sibling_variants, effective_quals, min_baseq,
            );
            let is_alt = result.is_alt && !claimed_by_sibling;
            let is_ref = result.is_ref
                && !sibling_claims_ref(
                    record, variant, &result, &ref_guard_siblings, effective_quals, min_baseq,
                    &mut alt_aligner, &mut ref_aligner, backend, rules,
                );

            // Only the resolved call reads this evidence (`resolve`), so whether the
            // read was informative is not needed here.
            let evidence = tx_fragments.entry(mol_hash).or_insert_with(FragmentEvidence::new);
            evidence.observe(is_ref, is_alt, result.qual, is_read1, is_forward, tlen, result.has_n_base, result.is_structural, record.mapq(), !result.uninformative);

            if first_class {
                if is_ref {
                    tx_rd += 1;
                } else if is_alt {
                    tx_ad += 1;
                }
            }
        }

        // ── Fragment resolution for this transcript
        let mut tx_adf: u32 = 0;
        let mut tx_rdf: u32 = 0;
        let mut tx_dpf: u32 = 0;
        let qual_diff = fragment_qual_threshold;

        for evidence in tx_fragments.values() {
            let (frag_ref, frag_alt) = evidence.resolve(qual_diff);
            tx_dpf += 1;
            if frag_ref { tx_rdf += 1; }
            else if frag_alt { tx_adf += 1; }
        }

        // Format: "ENST...:AD,RD,DP"
        read_entries.push(format!("{}:{},{},{}", tx_id, tx_ad, tx_rd, tx_dp));
        frag_entries.push(format!("{}:{},{},{}", tx_id, tx_adf, tx_rdf, tx_dpf));

        trace!(
            "per-transcript: {} → read AD={} RD={} DP={}, frag ADF={} RDF={} DPF={}",
            tx_id, tx_ad, tx_rd, tx_dp, tx_adf, tx_rdf, tx_dpf,
        );
    }

    if read_entries.is_empty() {
        return (String::new(), String::new());
    }

    // Separate transcripts with '|' (not ';'): ';' is the VCF INFO field separator, so a
    // ';'-joined value would corrupt VCF INFO parsing. Joining with '|' here makes MAF,
    // VCF, and the documented `ENST:AD,RD,DP|…` header all agree. (Within an entry,
    // ':'/',' are the field separators.)
    (read_entries.join("|"), frag_entries.join("|"))
}


// ══════════════════════════════════════════════════════════════════════════════
// ALLELE-SPECIFIC JUNCTION DIVERGENCE (ASJD)
// ══════════════════════════════════════════════════════════════════════════════

/// Result of ASJD detection for a single variant.
#[derive(Debug)]
struct AsjdResult {
    flag: bool,
    pval: f64,
    ref_junction: String,
    alt_junction: String,
    ref_motif: String,
    alt_motif: String,
    ref_known: bool,
    alt_known: bool,
    n_ref_junc: u32,
    n_alt_junc: u32,
    n_ref_total: u32,
    n_alt_total: u32,
    diagnostic: String,
}

impl AsjdResult {
    /// Empty result when ASJD detection is not applicable (no annotation,
    /// no junction reads, etc.).
    fn empty() -> Self {
        Self {
            flag: false,
            pval: 1.0,
            ref_junction: String::new(),
            alt_junction: String::new(),
            ref_motif: String::new(),
            alt_motif: String::new(),
            ref_known: false,
            alt_known: false,
            n_ref_junc: 0,
            n_alt_junc: 0,
            n_ref_total: 0,
            n_alt_total: 0,
            diagnostic: String::new(),
        }
    }
}

/// Per-junction forward/reverse read counts for strand discordance detection.
///
/// In stranded (dUTP) RNA-seq libraries, reads spanning a real splice junction
/// should originate from a single transcript strand. If the dominant ALT junction
/// has substantial support from **both transcript strands** (minority fraction
/// ≥ 30%), it indicates the junction may be an alignment artifact (e.g., DNA
/// contamination or mismapping) rather than a genuine RNA splice event.
///
/// Reads are binned by their dUTP-folded *transcript* strand
/// ([`super::rna::read_transcript_strand`]), not raw genomic orientation —
/// otherwise the two mates of a normal FR pair land on opposite genomic strands
/// and a genuine junction looks fully mixed, firing the flag spuriously.
///
/// The `STRAND_DISCORDANT` diagnostic flag fires when:
/// `min(plus, minus) / (plus + minus) >= 0.30`
#[derive(Default, Debug)]
struct JunctionStrandCounts {
    /// Reads whose transcript strand resolves to '+'.
    plus: u32,
    /// Reads whose transcript strand resolves to '-'.
    minus: u32,
}

impl JunctionStrandCounts {
    /// Total read count across both transcript strands.
    fn total(&self) -> u32 {
        self.plus + self.minus
    }

    /// Minority strand fraction: 0.0 = perfectly stranded, 0.5 = fully mixed.
    fn minority_strand_fraction(&self) -> f64 {
        let total = self.total();
        if total == 0 {
            return 0.0;
        }
        let minority = std::cmp::min(self.plus, self.minus);
        minority as f64 / total as f64
    }
}

/// Per-fragment junction-strand tally for ASJD.
///
/// A fragment (identified by its QNAME hash) votes once toward `n_total` and once per
/// junction it spans, regardless of how many of its mates cross. Counting both mates
/// of an overlapping pair would inflate junction totals (~1.38x on real RNA) and fire
/// spurious `STRAND_DISCORDANT` at low depth. Mates always fold to the same transcript
/// strand (verified 0/319k disagreements on real RNA), so the first vote per
/// (fragment, junction) wins unambiguously.
#[derive(Default)]
struct JunctionTally {
    /// Per-junction transcript-strand counts, deduped to one vote per fragment.
    counts: HashMap<(i64, i64), JunctionStrandCounts>,
    /// Distinct fragments contributing at least one junction-spanning read.
    n_total: u32,
    frag_seen: std::collections::HashSet<u64>,
    junc_seen: std::collections::HashSet<(i64, i64, u64)>,
    /// Mates skipped because their fragment already voted at that junction (monitoring).
    collapsed: u32,
}

impl JunctionTally {
    /// Record one junction-spanning read of a fragment.
    ///
    /// - `qhash`: fragment identity (QNAME hash).
    /// - `junctions`: the junctions this read spans.
    /// - `tx_minus`: folded transcript strand (`true` = '-'; unstranded reads pass
    ///   `false` so they count toward totals without skewing the strand split).
    fn add(&mut self, qhash: u64, junctions: &[(i64, i64)], tx_minus: bool) {
        if self.frag_seen.insert(qhash) {
            self.n_total += 1;
        }
        for &j in junctions {
            if self.junc_seen.insert((j.0, j.1, qhash)) {
                let entry = self.counts.entry(j).or_default();
                if tx_minus {
                    entry.minus += 1;
                } else {
                    entry.plus += 1;
                }
            } else {
                self.collapsed += 1;
            }
        }
    }
}

/// Coordinate tolerance (bp) for matching an observed junction endpoint to an
/// annotated one — shared by ASJD's known-junction test and the ASJD-2
/// anchoring test so the two cannot disagree about what "annotated" means.
const JUNCTION_TOLERANCE: i32 = 5;

/// Heuristic BAQ is skipped at variants whose REF span comes within this
/// distance (bp) of an annotated exon boundary: its CIGAR-N penalty would land
/// on exactly the reads that splice there, which are the evidence at an exon
/// edge.
const BAQ_BOUNDARY_SUPPRESS_BP: i32 = 5;

/// The span whose read indels BAQ spares for `v`: the variant's own event (grown
/// through repeats, or its discrimination window without a reference) plus two
/// flank bases. None for an SNV, whose reads' indels are never its evidence.
fn baq_own_span(v: &Variant) -> Option<(i64, i64)> {
    if v.ref_allele.len() == 1 && v.alt_allele.len() == 1 {
        return None;
    }
    let (lo, hi) = carrier::grown_event(v).unwrap_or_else(|| window::discrimination_window(v));
    Some((lo - 2, hi + 2))
}

/// Whether heuristic BAQ applies at a variant. The one rule every view that
/// classifies alleles uses (main counts, per-transcript counts, ASJD), so a
/// transcript's counts decompose the main counts: BAQ was requested, and no base
/// of the variant's REF span is within `BAQ_BOUNDARY_SUPPRESS_BP` of an annotated
/// exon boundary (`exon_boundary_dist`; None without a GTF).
fn baq_applies(apply_baq: bool, exon_boundary_dist: Option<i32>) -> bool {
    apply_baq && !matches!(exon_boundary_dist, Some(d) if d <= BAQ_BOUNDARY_SUPPRESS_BP)
}

/// Minimum fragments of REF-side junction evidence ASJD speaks on (below it:
/// `LOW_REF_JUNC`). Also the floor for `RETENTION_DOMINANT`, whose count is
/// the spliced (overwhelmingly wild-type) population — REF-side evidence.
const ASJD_MIN_REF_JUNC: u32 = 10;

/// Minimum fragments of ALT-side junction evidence ASJD speaks on (below it:
/// `LOW_ALT_JUNC`). Also the floor for `NOVEL_JUNC_AT_SPLICE_LOSS`, whose
/// count is the mutant allele's splicing outcome — ALT-side evidence.
const ASJD_MIN_ALT_JUNC: u32 = 5;

/// ASJD-2 splice-disruption markers for `asjd_diagnostic`.
///
/// ASJD's tallies only see allele-classified reads, but at a splice-site
/// variant the informative reads are often the ones the splice-aware
/// evidence rule excludes: reads whose CIGAR N spans the variant observe
/// nothing there. Both markers read that excluded population, and both are
/// population comparisons — no tuned rates:
///
/// - `RETENTION_DOMINANT(n)`: spliced-over fragments (`n`) outnumber the
///   allele-classified ones, and the classified fragments are mostly
///   junction-free. The reads that genotype this locus are then the
///   intron-retaining minority, so `vaf` is the VAF *within that
///   population*, not allelic balance (an allele-specific retention reads as
///   a very high `vaf`).
/// - `NOVEL_JUNC_AT_SPLICE_LOSS(n@start-end)`: the excluded population's
///   top unannotated junction — anchored to an annotated splice site (an
///   exon-skip or alternative-site event in the gene's splice graph, not
///   aligner noise) and not the deletion itself written as a splice — is
///   carried by more fragments than confirm ALT: the mutant allele's
///   splicing outcome is visible here while `ad` is not. Coordinates follow
///   the `asjd_*_junction` convention (0-based, half-open intron).
///
/// Both require the variant's REF span to reach within two bases of an
/// annotated intron boundary of a transcript on the variant's gene strand
/// (the canonical splice dinucleotide and the adjacent exonic bases;
/// transcript termini and antisense genes' sites do not count), and both
/// speak only above ASJD's own junction-evidence floors
/// (`ASJD_MIN_REF_JUNC` for the spliced population, `ASJD_MIN_ALT_JUNC` for
/// the novel junction). Returns the markers to append, possibly empty.
#[allow(clippy::too_many_arguments)]
fn splice_disruption_markers(
    annotation: &AnnotationIndex,
    chrom: &str,
    variant: &Variant,
    window_pad: i64,
    excluded: &JunctionTally,
    classified_frags: &std::collections::HashSet<u64>,
    classified_with_junc: usize,
    n_alt_frags: usize,
) -> Vec<String> {
    let span_start = variant.pos;
    let span_end = variant.pos + variant.ref_allele.len() as i64;
    // Span [s, e) intersects a boundary B's window [B-2, B+2) iff B lies in [s-1, e+1].
    if !annotation.intron_boundary_in_range(chrom, span_start - 1, span_end + 1, variant.gene_strand) {
        return Vec::new();
    }

    let mut markers = Vec::new();
    // A fragment with one classified mate is classified, not excluded.
    let n_excluded = excluded
        .frag_seen
        .iter()
        .filter(|qh| !classified_frags.contains(qh))
        .count();
    let n_classified = classified_frags.len();
    let junction_free = n_classified.saturating_sub(classified_with_junc);

    if n_excluded >= ASJD_MIN_REF_JUNC as usize
        && n_excluded > n_classified
        && junction_free > classified_with_junc
    {
        debug!(
            "ASJD-2 RETENTION_DOMINANT at {}:{}: {} spliced-over fragments vs {} classified \
             ({} junction-free)",
            variant.chrom, variant.pos + 1, n_excluded, n_classified, junction_free,
        );
        markers.push(format!("RETENTION_DOMINANT({})", n_excluded));
    }

    let del_len = variant.ref_allele.len().saturating_sub(variant.alt_allele.len()) as i64;
    let tol = JUNCTION_TOLERANCE as i64;
    let anchored = |x: i64| {
        annotation.intron_boundary_in_range(chrom, x - tol, x + tol, variant.gene_strand)
    };
    let top_novel = excluded
        .counts
        .iter()
        .filter(|(&(a, b), _)| {
            // The deletion itself, written by the aligner as a splice of the
            // same length at the locus, is the ALT carriers — not a splicing
            // consequence (SPLICE_SKIP_DOMINANT covers that representation).
            let is_deletion_as_splice = del_len > 0
                && b - a == del_len
                && (a - (span_start + 1)).abs() <= window_pad;
            !is_deletion_as_splice
                && !annotation.is_junction_known(chrom, a, b, JUNCTION_TOLERANCE)
                && (anchored(a) || anchored(b))
        })
        .map(|(&j, sc)| (sc.total(), j))
        // Highest depth first; ties broken by leftmost coordinates (deterministic).
        .max_by(|x, y| x.0.cmp(&y.0).then_with(|| y.1.cmp(&x.1)));
    if let Some((n, (a, b))) = top_novel {
        if n >= ASJD_MIN_ALT_JUNC && n as usize > n_alt_frags {
            debug!(
                "ASJD-2 NOVEL_JUNC_AT_SPLICE_LOSS at {}:{}: junction {}-{} on {} excluded \
                 fragments vs {} ALT fragments",
                variant.chrom, variant.pos + 1, a, b, n, n_alt_frags,
            );
            markers.push(format!("NOVEL_JUNC_AT_SPLICE_LOSS({}@{}-{})", n, a, b));
        } else {
            trace!(
                "ASJD-2: top anchored novel junction {}-{} ({} fragments) at {}:{} is below \
                 the ALT-evidence floor ({}) or does not exceed {} ALT fragments — no marker",
                a, b, n, variant.chrom, variant.pos + 1, ASJD_MIN_ALT_JUNC, n_alt_frags,
            );
        }
    }
    markers
}

/// Classify the splice motif at a junction by reading donor/acceptor dinucleotides
/// from the reference FASTA.
///
/// Splice junctions are defined by intron boundaries:
/// - **Donor** (5' end): 2bp at `junction_start` (first 2 bases of intron)
/// - **Acceptor** (3' end): 2bp at `junction_end - 2` (last 2 bases of intron)
///
/// | Donor | Acceptor | Motif  | Spliceosome |
/// |-------|----------|--------|-------------|
/// | GT    | AG       | GT-AG  | U2 major    |
/// | GC    | AG       | GC-AG  | U2 minor    |
/// | AT    | AC       | AT-AC  | U12         |
/// | other | other    | OTHER  | Non-canonical |
///
/// On the minus strand the donor/acceptor roles swap and each dinucleotide is
/// reverse-complemented, so a canonical minus-strand intron (genomic `CT..AC`) is
/// recognized rather than misread as `OTHER`. With strand unknown a canonical motif
/// in either orientation is accepted rather than guessing non-canonical.
/// Returns `"UNKNOWN"` if the FASTA reader is unavailable or a fetch fails.
fn classify_splice_motif(
    fasta_reader: &mut Option<bio::io::fasta::IndexedReader<std::fs::File>>,
    chrom: &str,
    junction_start: i64,
    junction_end: i64,
    gene_strand: Option<char>,
) -> String {
    let reader = match fasta_reader.as_mut() {
        Some(r) => r,
        None => return "UNKNOWN".to_string(),
    };

    // Fetch the two genomic dinucleotides bracketing the intron (forward strand):
    // `left` at the intron start, `right` at the intron end.
    let left = match crate::normalize::fasta::fetch_region(
        reader, chrom, junction_start as u64, (junction_start + 2) as u64,
    ) {
        Ok(bases) if bases.len() == 2 => {
            [bases[0].to_ascii_uppercase(), bases[1].to_ascii_uppercase()]
        }
        _ => return "UNKNOWN".to_string(),
    };
    let right = match crate::normalize::fasta::fetch_region(
        reader, chrom, (junction_end - 2) as u64, junction_end as u64,
    ) {
        Ok(bases) if bases.len() == 2 => {
            [bases[0].to_ascii_uppercase(), bases[1].to_ascii_uppercase()]
        }
        _ => return "UNKNOWN".to_string(),
    };

    // Orient donor/acceptor by gene strand. On '+', donor = left, acceptor = right.
    // On '-', the intron reads in reverse: donor = revcomp(right), acceptor =
    // revcomp(left). With strand unknown, accept a canonical motif in either reading.
    match gene_strand {
        Some('-') => motif_label(revcomp2(right), revcomp2(left)),
        Some('+') => motif_label(left, right),
        _ => canonical_motif(left, right)
            .or_else(|| canonical_motif(revcomp2(right), revcomp2(left)))
            .unwrap_or("OTHER")
            .to_string(),
    }
}

/// Reverse-complement a 2-base motif (donor/acceptor dinucleotide).
fn revcomp2(d: [u8; 2]) -> [u8; 2] {
    let complement = |b: u8| match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    };
    [complement(d[1]), complement(d[0])]
}

/// Canonical spliceosomal motif name for a donor/acceptor pair, or `None` if
/// non-canonical (GT-AG = U2 major, GC-AG = U2 minor, AT-AC = U12).
fn canonical_motif(donor: [u8; 2], acceptor: [u8; 2]) -> Option<&'static str> {
    match (donor, acceptor) {
        ([b'G', b'T'], [b'A', b'G']) => Some("GT-AG"),
        ([b'G', b'C'], [b'A', b'G']) => Some("GC-AG"),
        ([b'A', b'T'], [b'A', b'C']) => Some("AT-AC"),
        _ => None,
    }
}

/// Two junctions are the same splice event when both ends agree within
/// `JUNCTION_TOLERANCE`: aligners place a junction a few bases apart when the
/// sequence at the boundary is ambiguous.
fn same_junction(a: (i64, i64), b: (i64, i64)) -> bool {
    let tol = JUNCTION_TOLERANCE as i64;
    (a.0 - b.0).abs() <= tol && (a.1 - b.1).abs() <= tol
}

/// The junctions tied for the most fragments in one ASJD partition, leftmost
/// first. Empty for an empty partition.
fn top_junctions(counts: &HashMap<(i64, i64), JunctionStrandCounts>) -> Vec<(i64, i64)> {
    let max = counts.values().map(|sc| sc.total()).max().unwrap_or(0);
    let mut top: Vec<(i64, i64)> =
        counts.iter().filter(|(_, sc)| sc.total() == max).map(|(j, _)| *j).collect();
    top.sort_unstable();
    top
}

/// Of one partition's tied junctions (`top`, leftmost first, non-empty), the
/// one the other partition's fragments use most; the leftmost on a further tie.
fn most_supported(
    top: &[(i64, i64)],
    other: &HashMap<(i64, i64), JunctionStrandCounts>,
) -> (i64, i64) {
    let support = |j: &(i64, i64)| other.get(j).map_or(0, |sc| sc.total());
    // `top` is sorted, so the first maximum is the leftmost.
    let best = top.iter().map(support).max().unwrap_or(0);
    *top.iter().find(|j| support(j) == best).expect("top junctions are non-empty")
}

/// Canonical motif name, or `"OTHER"` when non-canonical.
fn motif_label(donor: [u8; 2], acceptor: [u8; 2]) -> String {
    canonical_motif(donor, acceptor).unwrap_or("OTHER").to_string()
}

/// Detect allele-specific junction divergence at a variant site.
///
/// Partitions reads from the cache into REF- and ALT-classified sets,
/// collects splice junctions from each partition, and tests whether
/// the junction distributions differ significantly (Fisher's exact test).
///
/// # Algorithm
///
/// 1. Re-classify reads (same as count_variant_from_cache) to partition into REF/ALT
/// 2. Extract CIGAR N ops from each partition → per-allele junction multiset with strand info
/// 3. Find the dominant junction in each partition
/// 4. If dominant junctions differ, run Fisher's exact 2x2 test
/// 5. Classify splice motifs (FASTA lookup), check GTF annotation, build diagnostic flags
/// 6. Check strand discordance on the ALT dominant junction (dUTP artifact detection)
///
/// # Parameters
///
/// - `read_cache`: All reads in the genomic bin (shared across variants)
/// - `variant`: The variant being analyzed
/// - `sibling_variants`: Multi-allelic sibling variants at the same locus
/// - `annotation`: The GTF annotation index
/// - `min_mapq`, `min_baseq`, etc.: Standard counting parameters
/// - `use_baq`: whether heuristic BAQ applies at this variant, as resolved by
///   `baq_applies` (the main counts' rule)
#[allow(clippy::too_many_arguments)]
fn detect_asjd(
    read_cache: &[Record],
    variant: &Variant,
    sibling_variants: &[Variant],
    annotation: &AnnotationIndex,
    min_mapq: u8,
    min_baseq: u8,
    backend: &AlignmentBackend,
    use_baq: bool,
    enforce_strandedness: bool,
    strandedness: rna::Strandedness,
    fasta_reader: &mut Option<bio::io::fasta::IndexedReader<std::fs::File>>,
    rules: &carrier::ReadRules,
) -> AsjdResult {
    let baq_spare = if use_baq { baq_own_span(variant) } else { None };
    // Bind the normalized key, then borrow it as &str so the junction/motif lookups
    // below are unchanged.
    let chrom_key = crate::shared::contig::normalize_contig(&variant.chrom);
    let chrom = chrom_key.as_str();
    let window_pad = window::scan_pad(variant);
    let (v_start, v_end) = window::read_window(variant);

    // Create aligners
    let mut alt_aligner = Aligner::new(SW_GAP_OPEN, SW_GAP_EXTEND, &sw_score);
    let mut ref_aligner = Aligner::new(SW_GAP_OPEN, SW_GAP_EXTEND, &sw_score);

    // Step 1: Partition reads into REF/ALT and collect junctions with strand info.
    // Tracks per-junction transcript-strand counts for STRAND_DISCORDANT detection:
    // a genuine splice junction is supported from a single transcript strand, so
    // mixed-strand evidence suggests an alignment artifact.
    // ASJD junction evidence is counted per FRAGMENT, not per mate. On real RNA, 35.6%
    // of fragment×junction incidences have both mates spanning the same junction;
    // counting both inflates junction totals ~1.38x on average and fires spurious
    // STRAND_DISCORDANT at low depth. `JunctionTally` dedups each fragment to one vote
    // per allele-total and per junction (mates always fold to the same transcript
    // strand — verified 0/319k disagreements on real RNA — so first-seen wins).
    let mut ref_tally = JunctionTally::default();
    let mut alt_tally = JunctionTally::default();
    // Splice-disruption populations (ASJD-2): fragments whose CIGAR N spans
    // the variant (excluded from counting as no-observation) with the
    // junctions they carry over the locus, and the allele-classified
    // fragments split by whether they carry any junction.
    let mut excluded_tally = JunctionTally::default();
    let mut classified_frags: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut alt_frags: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let var_span = (variant.pos, variant.pos + variant.ref_allele.len() as i64);

    for record in read_cache {
        let r_start = record.pos();
        let r_end = ref_end(record);
        if r_start >= v_end || r_end <= v_start {
            continue;
        }

        // Standard filters
        if record.mapq() < min_mapq && !super::rna::is_valid_rna_alignment(record, min_mapq) {
            continue;
        }
        if enforce_strandedness && !super::rna::is_sense_strand(record, variant.gene_strand, strandedness) {
            continue;
        }

        // BAQ under the main counts' rule (`use_baq`): at an exon edge the
        // spliced reads are exactly ASJD's evidence.
        let quals = effective_quals(record, use_baq, baq_spare);
        let effective_quals: &[u8] = &quals;

        // Classify allele
        let result = check_allele_with_qual(
            record, variant, sibling_variants, effective_quals, min_baseq,
            &mut alt_aligner, &mut ref_aligner, backend, rules,
        );

        let qh = crate::shared::fragment::hash_qname(record.qname());
        if result.is_ref || result.is_alt {
            classified_frags.insert(qh);
            if result.is_alt {
                alt_frags.insert(qh);
            }
        }

        // Only interested in reads with splice junctions
        let junctions = super::rna::extract_splice_junctions(record);
        if junctions.is_empty() {
            continue;
        }

        // Transcript-strand-folded under the library protocol, so both mates of a
        // normal FR pair agree rather than splitting a genuine junction across both
        // genomic strands. For an unstranded library there is no transcript strand
        // (None) — such reads fall to `plus` so they still count toward the junction
        // total, while the strand split (used only for discordance below) is left
        // inert and the STRAND_DISCORDANT emission is gated off.
        let tx_minus = super::rna::read_transcript_strand(record, strandedness) == Some('-');

        // Dedup by fragment (QNAME hash): vote once per allele-total and per junction.
        if result.is_ref {
            ref_tally.add(qh, &junctions, tx_minus);
        } else if result.is_alt {
            alt_tally.add(qh, &junctions, tx_minus);
        } else if !result.covers_locus {
            // Spliced over the locus: keep only the junction(s) spanning the
            // variant — the ones that explain why the read saw nothing here.
            let over_locus: Vec<(i64, i64)> = junctions
                .iter()
                .copied()
                .filter(|&(a, b)| a < var_span.1 && var_span.0 < b)
                .collect();
            if !over_locus.is_empty() {
                excluded_tally.add(qh, &over_locus, tx_minus);
            }
        }
    }

    // Unpack the deduped tallies into the names the rest of the routine reads.
    let n_ref_total = ref_tally.n_total;
    let n_alt_total = alt_tally.n_total;
    let collapsed_mates = ref_tally.collapsed + alt_tally.collapsed;
    let ref_junction_counts = ref_tally.counts;
    let alt_junction_counts = alt_tally.counts;

    if collapsed_mates > 0 {
        debug!(
            "ASJD QNAME-dedup at {}:{}: collapsed {} duplicate junction-spanning mates",
            variant.chrom, variant.pos + 1, collapsed_mates,
        );
    }

    // Step 2: Check minimum evidence thresholds
    let mut diag_flags: Vec<String> = Vec::new();

    if n_ref_total < ASJD_MIN_REF_JUNC {
        diag_flags.push("LOW_REF_JUNC".to_string());
    }
    if n_alt_total < ASJD_MIN_ALT_JUNC {
        diag_flags.push("LOW_ALT_JUNC".to_string());
    }

    // Step 2b: splice-disruption markers (ASJD-2). Computed before the
    // no-junction early return below: they explain exactly the loci where
    // the classified reads carry no junctions (LOW_*_JUNC).
    let classified_with_junc = ref_tally
        .frag_seen
        .union(&alt_tally.frag_seen)
        .filter(|qh| classified_frags.contains(qh))
        .count();
    diag_flags.extend(splice_disruption_markers(
        annotation,
        chrom,
        variant,
        window_pad,
        &excluded_tally,
        &classified_frags,
        classified_with_junc,
        alt_frags.len(),
    ));

    // If either partition has no junction reads, no divergence can be detected
    if ref_junction_counts.is_empty() || alt_junction_counts.is_empty() {
        return AsjdResult {
            diagnostic: diag_flags.join(";"),
            n_ref_total,
            n_alt_total,
            ..AsjdResult::empty()
        };
    }

    // Step 3: Find dominant junction in each partition (by total fragments across
    // both strands). A tie is read the least divergent way, from the reads, never
    // from hash order (which varies run to run):
    // - when a REF top junction and an ALT top junction are the same splice
    //   event (`same_junction`, an exact match preferred), each partition
    //   reports its own of that pair — a tie is not a divergence;
    // - otherwise each partition reports the tied junction the other partition
    //   supports most, then the leftmost.
    let ref_top = top_junctions(&ref_junction_counts);
    let alt_top = top_junctions(&alt_junction_counts);
    let pairs = || ref_top.iter().flat_map(|r| alt_top.iter().map(move |a| (*r, *a)));
    let (ref_dom_junc, alt_dom_junc) = pairs()
        .find(|(r, a)| r == a)
        .or_else(|| pairs().find(|&(r, a)| same_junction(r, a)))
        .unwrap_or_else(|| {
            (
                most_supported(&ref_top, &alt_junction_counts),
                most_supported(&alt_top, &ref_junction_counts),
            )
        });
    let ref_dom_strand_info = &ref_junction_counts[&ref_dom_junc];
    let n_ref_junc = ref_dom_strand_info.total();
    let alt_dom_strand_info = &alt_junction_counts[&alt_dom_junc];
    let n_alt_junc = alt_dom_strand_info.total();

    // Check for multi-junction in ALT reads
    let alt_distinct_junctions = alt_junction_counts.len();
    if alt_distinct_junctions > 2 {
        diag_flags.push("MULTI_JUNCTION".to_string());
    }

    // Step 4: Compare dominant junctions
    let is_same = same_junction(ref_dom_junc, alt_dom_junc);
    let (pval, flag) = if is_same {
        // Same dominant junction → no divergence
        (1.0, false)
    } else {
        // Different dominant junctions → Fisher's exact test
        // 2x2 table: [REF on ref_junc, REF on alt_junc] vs [ALT on ref_junc, ALT on alt_junc]
        let ref_on_ref_junc = ref_junction_counts.get(&ref_dom_junc).map_or(0, |sc| sc.total());
        let ref_on_alt_junc = ref_junction_counts.get(&alt_dom_junc).map_or(0, |sc| sc.total());
        let alt_on_ref_junc = alt_junction_counts.get(&ref_dom_junc).map_or(0, |sc| sc.total());
        let alt_on_alt_junc = alt_junction_counts.get(&alt_dom_junc).map_or(0, |sc| sc.total());

        let (p, _or) = crate::shared::stats::fisher_exact_2x2(
            ref_on_ref_junc, ref_on_alt_junc,
            alt_on_ref_junc, alt_on_alt_junc,
        );

        (p, p < 0.05 && n_alt_junc >= ASJD_MIN_ALT_JUNC && n_ref_junc >= ASJD_MIN_REF_JUNC)
    };

    // Step 5: Classify splice motifs and GTF annotation
    let ref_known = annotation.is_junction_known(chrom, ref_dom_junc.0, ref_dom_junc.1, JUNCTION_TOLERANCE);
    let alt_known = annotation.is_junction_known(chrom, alt_dom_junc.0, alt_dom_junc.1, JUNCTION_TOLERANCE);

    if !alt_known && !is_same {
        diag_flags.push("NOVEL_ALT_JUNC".to_string());
    }

    let ref_junction_str = format!("{}-{}", ref_dom_junc.0, ref_dom_junc.1);
    let alt_junction_str = format!("{}-{}", alt_dom_junc.0, alt_dom_junc.1);

    // Step 5b: Classify splice motifs from reference sequence
    // Splice junctions have a donor (5') and acceptor (3') dinucleotide:
    //   Donor:   2bp at junction_start (intron start)
    //   Acceptor: 2bp at junction_end - 2 (intron end)
    // Canonical motifs: GT-AG (U2), GC-AG (U2 minor), AT-AC (U12)
    let ref_motif = classify_splice_motif(
        fasta_reader, chrom, ref_dom_junc.0, ref_dom_junc.1, variant.gene_strand,
    );
    let alt_motif = classify_splice_motif(
        fasta_reader, chrom, alt_dom_junc.0, alt_dom_junc.1, variant.gene_strand,
    );

    // NON_CANONICAL_MOTIF diagnostic flag
    if !is_same && alt_motif == "OTHER" {
        diag_flags.push("NON_CANONICAL_MOTIF".to_string());
    }

    // Step 5c: Strand discordance detection on the dominant ALT junction.
    // In a stranded library, reads from a genuine splice event should be predominantly
    // on one transcript strand. A minority strand fraction ≥ 30% (with ≥5 total reads
    // to avoid noise at low depth) indicates the junction may be an alignment artifact.
    // Undefined for an unstranded library (no transcript strand), so gate it off there.
    let alt_minority_frac = alt_dom_strand_info.minority_strand_fraction();
    if strandedness != rna::Strandedness::Unstranded
        && !is_same
        && n_alt_junc >= ASJD_MIN_ALT_JUNC
        && alt_minority_frac >= 0.30
    {
        diag_flags.push("STRAND_DISCORDANT".to_string());
        debug!(
            "STRAND_DISCORDANT: {}:{} alt_junc={}-{} tx+={} tx-={} minority_frac={:.2}",
            variant.chrom, variant.pos + 1,
            alt_dom_junc.0, alt_dom_junc.1,
            alt_dom_strand_info.plus, alt_dom_strand_info.minus,
            alt_minority_frac,
        );
    }

    trace!(
        "ASJD: {}:{} flag={} pval={:.4e} ref={}({}) alt={}({}) ref_n={}/{} alt_n={}/{} alt_strand={}tx+/{}tx-",
        variant.chrom, variant.pos + 1, flag, pval,
        ref_junction_str, ref_motif, alt_junction_str, alt_motif,
        n_ref_junc, n_ref_total, n_alt_junc, n_alt_total,
        alt_dom_strand_info.plus, alt_dom_strand_info.minus,
    );

    AsjdResult {
        flag,
        pval,
        ref_junction: ref_junction_str,
        alt_junction: alt_junction_str,
        ref_motif,
        alt_motif,
        ref_known,
        alt_known,
        n_ref_junc,
        n_alt_junc,
        n_ref_total,
        n_alt_total,
        diagnostic: diag_flags.join(";"),
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;
    use std::ffi::CString;

    #[test]
    fn test_baq_applies_threshold() {
        // Skipped within BAQ_BOUNDARY_SUPPRESS_BP of an annotated exon edge,
        // applied beyond it and without annotation, never when not requested.
        assert!(!baq_applies(true, Some(0)));
        assert!(!baq_applies(true, Some(BAQ_BOUNDARY_SUPPRESS_BP)));
        assert!(baq_applies(true, Some(BAQ_BOUNDARY_SUPPRESS_BP + 1)));
        assert!(baq_applies(true, None));
        assert!(!baq_applies(false, Some(50)));
        assert!(!baq_applies(false, None));
    }

    #[test]
    fn test_same_junction_tolerance() {
        let j = (300, 500);
        assert!(same_junction(j, j));
        assert!(same_junction(j, (300 + JUNCTION_TOLERANCE as i64, 500 - JUNCTION_TOLERANCE as i64)));
        assert!(!same_junction(j, (300, 500 + JUNCTION_TOLERANCE as i64 + 1)));
    }

    #[test]
    fn test_splice_motif_orientation_by_strand() {
        // A canonical GT-AG intron on the MINUS strand appears on the forward genome
        // as left=CT (revcomp of acceptor AG) and right=AC (revcomp of donor GT).
        let left_minus = *b"CT";
        let right_minus = *b"AC";
        // Read forward (plus orientation) it is non-canonical...
        assert_eq!(motif_label(left_minus, right_minus), "OTHER");
        // ...but minus orientation (revcomp + donor/acceptor swap) recovers GT-AG.
        assert_eq!(motif_label(revcomp2(right_minus), revcomp2(left_minus)), "GT-AG");

        // Plus-strand canonical intron: left=GT donor, right=AG acceptor.
        assert_eq!(motif_label(*b"GT", *b"AG"), "GT-AG");

        // Building-block sanity.
        assert_eq!(revcomp2(*b"AC"), *b"GT");
        assert_eq!(revcomp2(*b"CT"), *b"AG");
        assert_eq!(canonical_motif(*b"AT", *b"AC"), Some("AT-AC"));
        assert_eq!(canonical_motif(*b"CC", *b"GG"), None);
    }

    /// Build a synthetic BAM record for testing.
    ///
    /// Creates a minimal Record with the given sequence, qualities, CIGAR,
    /// and mapping position. All other fields are set to sensible defaults
    /// (unmapped=false, mapq=60, proper_pair=true, forward strand).
    fn build_record(seq: &[u8], qual: &[u8], cigar: &CigarString, pos: i64) -> Record {
        let mut record = Record::new();
        let qname = CString::new("test_read").unwrap();
        record.set(
            qname.as_bytes(), // qname
            Some(cigar),       // cigar
            seq,               // seq
            qual,              // qual
        );
        record.set_pos(pos);
        record.set_tid(0);
        record.set_mapq(60);
        // Set as mapped, proper pair, first in template
        record.set_flags(0x01 | 0x02 | 0x40); // paired + proper + read1
        record
    }

    /// Helper: build a Variant with the given position, ref, and alt alleles.
    fn build_variant(pos: i64, ref_allele: &str, alt_allele: &str) -> Variant {
        Variant {
            chrom: "1".to_string(),
            pos,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            variant_type: String::new(),
            ref_context: None,
            ref_context_start: 0,
            repeat_span: 0,
            gene_strand: None,
            shift_region: None,
            event_ref: None,
            boundary_span: None,
        }
    }

    /// Build a minimal single-contig HeaderView for binning tests.
    fn build_header_with_contig(name: &str, len: i64) -> bam::HeaderView {
        let mut header = bam::Header::new();
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", name);
        sq.push_tag(b"LN", len);
        header.push_record(&sq);
        bam::HeaderView::from_header(&header)
    }

    // ── contig-naming reconciliation (chr1 vs 1) ──

    #[test]
    fn test_resolve_tid_reconciles_contig_naming() {
        // BAM uses UCSC "chr1"; a variant may arrive as "1" (or vice-versa). resolve_tid
        // must reconcile both so reads are found — the fix for the silent zero-count when
        // naming differs — and return None only for a genuine miss (which then warns).
        let ucsc = build_header_with_contig("chr1", 1_000);
        assert_eq!(resolve_tid(&ucsc, "chr1"), Some(0)); // exact
        assert_eq!(resolve_tid(&ucsc, "1"), Some(0)); // "1" reconciled to "chr1"
        assert_eq!(resolve_tid(&ucsc, "chr2"), None); // genuine miss

        let b37 = build_header_with_contig("1", 1_000);
        assert_eq!(resolve_tid(&b37, "1"), Some(0)); // exact
        assert_eq!(resolve_tid(&b37, "chr1"), Some(0)); // "chr1" reconciled to "1"
        assert_eq!(resolve_tid(&b37, "chrX"), None); // genuine miss
    }

    // ── bin fetch-end must cover the anchor variant's full ref span ──

    #[test]
    fn test_bin_covers_anchor_deletion_ref_span() {
        // A bin anchored by a deletion whose ref span exceeds the window must
        // still fetch the deletion's right breakpoint. Before the anchor-aware
        // seed, bin.end stopped at bin_start + window and dropped reads aligned
        // past it (clip-admitted carriers of a long event).
        let header = build_header_with_contig("1", 1_000_000);
        let window = 10_000_i64;

        // Anchor is a 12kb deletion (ref span > window); a SNP sits within it.
        let big_ref = "A".repeat(12_000);
        let variants = vec![
            build_variant(1_000, &big_ref, "A"),
            build_variant(2_000, "C", "T"),
        ];

        let bins = build_genomic_bins(&variants, &header, window, BIN_MAX_VARIANTS);
        assert_eq!(bins.len(), 1, "both variants belong in one bin");
        let bin = &bins[0];

        // End-coverage invariant for EVERY variant, anchor included.
        for (k, v) in variants.iter().enumerate() {
            let ref_end = v.pos + v.ref_allele.len() as i64; // exclusive ref end
            assert!(
                bin.end >= ref_end,
                "bin.end ({}) must cover variant {} ref end ({})",
                bin.end,
                k,
                ref_end,
            );
            assert!(bin.start <= v.pos, "bin.start must cover variant {} pos", k);
        }
        // The anchor deletion's right breakpoint is at 1000 + 12000 = 13000.
        assert!(
            bin.end >= 13_000,
            "anchor deletion right breakpoint (13000) must be fetched, got bin.end={}",
            bin.end,
        );
    }

    #[test]
    fn every_variant_lands_in_one_bin_whose_fetch_holds_its_window() {
        // Bin geometry is performance only: under any window and cap, each variant
        // belongs to exactly one bin, and that bin's fetch holds the variant's whole
        // read window (`window::read_window`, the filter each
        // variant applies to the bin's cached reads). Random clusters of SNVs and
        // indels up to 300bp, with repeat spans, over a deterministic LCG.
        let header = build_header_with_contig("1", 1_000_000);
        let mut state: u64 = 0x5eed;
        let mut next = |n: u64| -> u64 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) % n
        };
        for _trial in 0..200 {
            let count = 1 + next(40) as usize;
            let mut variants = Vec::with_capacity(count);
            for _ in 0..count {
                let pos = 1 + next(30_000) as i64;
                let ref_len = if next(3) == 0 { 1 + next(300) as usize } else { 1 };
                let mut v = build_variant(pos, &"A".repeat(ref_len), "A");
                v.repeat_span = next(40) as usize;
                variants.push(v);
            }
            for &window in &[1_i64, 7, 100, 10_000] {
                for &cap in &[1_usize, 2, 3, 200] {
                    let bins = build_genomic_bins(&variants, &header, window, cap);
                    let mut seen = vec![0_usize; variants.len()];
                    for bin in &bins {
                        assert!(bin.variant_indices.len() <= cap);
                        for &k in &bin.variant_indices {
                            seen[k] += 1;
                            let v = &variants[k];
                            let (lo, hi) = window::read_window(v);
                            assert!(
                                bin.start <= lo && bin.end >= hi,
                                "window {window} cap {cap}: bin [{}, {}) misses variant at {} \
                                 (read window [{lo}, {hi}))",
                                bin.start, bin.end, v.pos,
                            );
                        }
                    }
                    assert!(seen.iter().all(|&n| n == 1), "window {window} cap {cap}: {seen:?}");
                }
            }
        }
    }

    // ── ASJD per-fragment junction dedup ──

    #[test]
    fn test_junction_tally_dedup_removes_spurious_discordance() {
        // A junction supported by 3 single-mate '+' fragments and 1 fragment whose
        // BOTH mates span it on '-' (the overlapping-mate case, common on real RNA).
        let mut t = JunctionTally::default();
        let j = (100i64, 200i64);
        t.add(1, &[j], false); // frag 1, '+'
        t.add(2, &[j], false); // frag 2, '+'
        t.add(3, &[j], false); // frag 3, '+'
        t.add(4, &[j], true); // frag 4, '-' (mate A)
        t.add(4, &[j], true); // frag 4, '-' (mate B) → must collapse

        let c = &t.counts[&j];
        assert_eq!((c.plus, c.minus), (3, 1), "fragment 4 votes once, not twice");
        assert_eq!(t.n_total, 4, "four distinct fragments");
        assert_eq!(t.collapsed, 1, "frag 4's second mate collapsed");
        // Deduped minority fraction 1/4 = 0.25 is below the 0.30 STRAND_DISCORDANT
        // threshold; the raw per-mate count (+3/-2 = 0.40) would have fired it.
        assert!(
            c.minority_strand_fraction() < 0.30,
            "deduped fraction must clear the discordance threshold"
        );
    }

    #[test]
    fn test_junction_tally_fragment_spans_two_junctions() {
        // One fragment spanning two distinct junctions votes once per junction and
        // counts once toward the allele total (no collapse — different junctions).
        let mut t = JunctionTally::default();
        let (j1, j2) = ((100i64, 200i64), (300i64, 400i64));
        t.add(1, &[j1, j2], false);
        assert_eq!(t.n_total, 1, "one fragment");
        assert_eq!(t.counts[&j1].plus, 1);
        assert_eq!(t.counts[&j2].plus, 1);
        assert_eq!(t.collapsed, 0, "distinct junctions are not collapses");
    }

    // ── Phase 0 N-base defense tests ──

    #[test]
    fn test_check_snp_n_base_high_bq() {
        // N base with BQ=30 should be classified as neither, not third allele.
        // Defense-in-depth: fgbio assigns BQ ≈ 2 to N, but raw BAMs may not.
        // Without the explicit N guard in check_snp, N would pass the BQ gate
        // (30 >= 20) and fall through to the REF/ALT comparison as "neither"
        // with qual=0 — correct result but wrong reason (no explicit N handling).
        let seq = b"AANTTCC";
        let qual = &[30, 30, 30, 30, 30, 30, 30]; // N at pos 2 has BQ=30
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "A", "T");
        let result = check_snp(&record, &variant, record.qual(), 20);

        assert!(!result.is_ref, "N should not be classified as REF");
        assert!(!result.is_alt, "N should not be classified as ALT");
        assert_eq!(result.qual, 0, "N should have qual=0 (uninformative)");
        assert!(result.has_n_base, "SNP N base should set has_n_base=true");
    }

    #[test]
    fn test_check_snp_n_base_low_bq_also_filtered() {
        // N base with BQ=2 (typical duplex) — caught by BQ gate first.
        // Verifies the BQ gate is the first line of defense.
        let seq = b"AANTTCC";
        let qual = &[30, 30, 2, 30, 30, 30, 30]; // N at pos 2 has BQ=2
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "A", "T");
        let result = check_snp(&record, &variant, record.qual(), 20);

        assert!(!result.is_ref, "N with low BQ should not be REF");
        assert!(!result.is_alt, "N with low BQ should not be ALT");
        // Low BQ N is caught by the BQ gate first (before the N check),
        // so has_n_base is NOT set — it never reaches the N guard.
        // This is correct: the read is filtered for quality reasons,
        // not N-masking reasons. The has_n_base flag is for reads that
        // *pass* BQ but carry N (defense-in-depth).
    }

    // ── MNP check_mnp tests ──

    #[test]
    fn test_mnp_high_quality_matches_ref() {
        // 2bp MNP: REF=AT, ALT=CG. Read has AT (matches REF).
        // All qualities >= 20.
        let seq = b"AAATGCC";
        let qual = &[30, 30, 30, 35, 30, 30, 30]; // all high quality
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        // Variant at pos 102-103 (0-based), REF=AT, ALT=CG
        let variant = build_variant(102, "AT", "CG");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        match result {
            MnpResult::Ref(q, _) => assert!(q >= 20, "Quality should be >= 20, got {}", q),
            other => panic!("Expected MnpResult::Ref, got {:?}", format_mnp_result(&other)),
        }
    }

    #[test]
    fn test_mnp_high_quality_matches_alt() {
        // 2bp MNP: REF=AT, ALT=CG. Read has CG (matches ALT).
        let seq = b"AACGGCC";
        let qual = &[30, 30, 35, 38, 30, 30, 30];
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "AT", "CG");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        match result {
            MnpResult::Alt(q, _, confirmed) => {
                assert!(q >= 20, "Quality should be >= 20, got {}", q);
                assert!(confirmed, "both discriminating bases read as ALT → confirmed");
            }
            other => panic!("Expected MnpResult::Alt, got {:?}", format_mnp_result(&other)),
        }
    }

    #[test]
    fn test_mnp_one_low_quality_base_classified_by_unmasked() {
        // 2bp MNP: REF=AT, ALT=CG. Read has CG (matches ALT).
        // Pos 0 (C): Q=35 ≥ 20 → unmasked, matches ALT.
        // Pos 1 (G): Q=8 < 20 → masked, cannot vote.
        // OLD behavior: min(35, 8) = 8 < 20 → LowQuality (entire read dropped).
        // NEW behavior: 1 unmasked position matches ALT → Alt (read recovered).
        // This is the key improvement: masked per-position evaluation
        // recovers reads where only some discriminating positions are low-quality.
        let seq = b"AACGGCC";
        let qual = &[30, 30, 35, 8, 30, 30, 30]; // pos 3: Q=8 < 20 → masked
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "AT", "CG");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        match result {
            MnpResult::Alt(q, _, _) => assert!(q > 0, "Expected Alt with quality > 0, got {}", q),
            other => panic!(
                "Expected Alt for DNP with one masked position (recovered by masked eval), got {:?}",
                format_mnp_result(&other)
            ),
        }
    }

    #[test]
    fn test_mnp_all_low_quality() {
        // Both bases below threshold.
        let seq = b"AACGGCC";
        let qual = &[30, 30, 5, 8, 30, 30, 30]; // pos 2: Q=5, pos 3: Q=8
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "AT", "CG");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        assert!(
            matches!(result, MnpResult::LowQuality(..)),
            "Expected LowQuality when all bases below threshold"
        );
    }

    #[test]
    fn test_mnp_third_allele() {
        // 2bp MNP: REF=AT, ALT=CG. Read has TT (matches neither).
        // All qualities pass threshold.
        let seq = b"AATTGCC";
        let qual = &[30, 30, 35, 38, 30, 30, 30];
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "AT", "CG");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        assert!(
            matches!(result, MnpResult::ThirdAllele(..)),
            "Expected ThirdAllele when bases match neither REF nor ALT"
        );
    }

    #[test]
    fn test_mnp_indel_in_block_structural() {
        // 3bp MNP: REF=ATG, ALT=CGC. Read has CIGAR 1M 1I 1D 2M,
        // creating a non-contiguous mapping across the MNP block.
        // The contiguity check should detect this and return Structural.
        let seq = b"AACXGCC"; // 7 bases, but CIGAR rearranges alignment
        let qual = &[30, 30, 35, 38, 35, 30, 30];
        // CIGAR: 2M 1I 1D 3M = consumes 2+0+1+3=6 ref, 2+1+0+3=6 read
        let cigar = CigarString(vec![
            Cigar::Match(2),
            Cigar::Ins(1),
            Cigar::Del(1),
            Cigar::Match(3),
        ]);
        let record = build_record(seq, qual, &cigar, 100);

        // Variant spans pos 102-104 (3bp MNP)
        // With the indel, positions won't be contiguous in read space
        let variant = build_variant(102, "ATG", "CGC");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        assert!(
            matches!(result, MnpResult::Structural),
            "Expected Structural when indel exists in MNP block"
        );
    }

    #[test]
    fn test_mnp_partial_coverage_structural() {
        // 3bp MNP but read only covers first 2 positions.
        // Read: 3 bases starting at pos 101, so covers 101-103.
        // Variant: pos 102-104 → pos 104 not covered.
        let seq = b"ACG";
        let qual = &[30, 35, 38];
        let cigar = CigarString(vec![Cigar::Match(3)]);
        let record = build_record(seq, qual, &cigar, 101);

        let variant = build_variant(102, "ATG", "CGC");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        assert!(
            matches!(result, MnpResult::Structural),
            "Expected Structural when read doesn't fully cover MNP"
        );
    }

    #[test]
    fn test_mnp_position_not_found_structural() {
        // Read starts after the variant position.
        let seq = b"ATCGATCG";
        let qual = &[30; 8];
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let record = build_record(seq, qual, &cigar, 200); // far beyond variant

        let variant = build_variant(102, "AT", "CG");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        assert!(
            matches!(result, MnpResult::Structural),
            "Expected Structural when variant position not in read"
        );
    }

    #[test]
    fn test_mnp_5bp_tert_pattern() {
        // 5bp MNP mimicking TERT: GAGGG→AAGGA
        // Read carries the ALT allele AAGGA at high quality.
        let seq = b"CCAAGGATTT";
        let qual = &[30, 30, 35, 38, 32, 36, 34, 30, 30, 30];
        let cigar = CigarString(vec![Cigar::Match(10)]);
        let record = build_record(seq, qual, &cigar, 100);

        // Variant at pos 102-106
        let variant = build_variant(102, "GAGGG", "AAGGA");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        match result {
            MnpResult::Alt(q, _, _) => assert!(q >= 20, "Quality should be >= 20, got {}", q),
            other => panic!("Expected MnpResult::Alt for TERT-like 5bp MNP, got {:?}",
                           format_mnp_result(&other)),
        }
    }

    #[test]
    fn test_mnp_min_bq_boundary() {
        // Test the exact boundary: min(BQ) == min_baseq should PASS.
        let seq = b"AACGGCC";
        let qual = &[30, 30, 20, 30, 30, 30, 30]; // min BQ = 20 == threshold
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "AT", "CG");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        match result {
            MnpResult::Alt(q, _, _) => assert!(q >= 20, "Should pass at exact threshold boundary"),
            other => panic!("Expected Alt at exact BQ threshold, got {:?}",
                           format_mnp_result(&other)),
        }
    }

    /// Helper to format MnpResult for panic messages
    fn format_mnp_result(result: &MnpResult) -> String {
        match result {
            MnpResult::Ref(q, n) => format!("Ref(q={}, had_n={})", q, n),
            MnpResult::Alt(q, n, c) => format!("Alt(q={}, had_n={}, confirmed={})", q, n, c),
            MnpResult::LowQuality(p, n) => format!("LowQuality(partial={}, had_n={})", p, n),
            MnpResult::ThirdAllele(p, n) => format!("ThirdAllele(partial={}, had_n={})", p, n),
            MnpResult::Structural => "Structural".to_string(),
        }
    }

    // ── Selective discriminating-position quality gate regression tests ──

    #[test]
    fn test_mnp_low_qual_non_discriminating_passes() {
        // TERT-like ONP: GAGGG→AAGGA (5bp, positions 0 and 4 discriminating).
        // Positions 1-3 are non-discriminating (REF==ALT: A=A, G=G, G=G).
        // Read carries AAGGA (ALT). Low quality at non-discriminating pos 2.
        // Old behavior: min(BQ) across ALL = 5 < 20 → LowQuality (WRONG).
        // New behavior: min(BQ) across discriminating (pos 0, 4) = 32 ≥ 20 → ALT (CORRECT).
        let seq = b"CCAAGGACC";
        let qual = &[30, 30, 35, 38, 5, 36, 32, 30, 30]; // pos 4 (offset +2): Q=5 (non-discriminating)
        let cigar = CigarString(vec![Cigar::Match(9)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "GAGGG", "AAGGA");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        match result {
            MnpResult::Alt(q, _, _) => assert!(q >= 20,
                "Should classify as ALT when non-discriminating bases are low quality, got q={}", q),
            other => panic!(
                "Expected Alt for TERT ONP with low-qual non-discriminating pos, got {:?}",
                format_mnp_result(&other)),
        }
    }

    #[test]
    fn test_mnp_low_qual_one_discriminating_classified_by_unmasked() {
        // All-discriminating DNP: GG→AA (both positions differ).
        // Read carries AA (ALT). Pos 0: Q=8 < 20 → masked. Pos 1: Q=38 ≥ 20 → unmasked.
        // OLD behavior: min(discriminating BQ) = 8 < 20 → LowQuality.
        // NEW behavior: 1 unmasked position matches ALT → Alt (read recovered).
        let seq = b"CCAAGCC";
        let qual = &[30, 30, 8, 38, 30, 30, 30]; // pos 2 (discriminating): Q=8
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "GG", "AA");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        match result {
            MnpResult::Alt(q, _, _) => assert!(q > 0, "Expected Alt with quality > 0, got {}", q),
            other => panic!(
                "Expected Alt for DNP with one masked discriminating position, got {:?}",
                format_mnp_result(&other)
            ),
        }
    }

    #[test]
    fn test_mnp_third_allele_partial_alt_match() {
        // TERT-like ONP: GAGGG→AAGGA. Read carries AAGGG (only pos 0 mutated).
        // This is the typical misannotated compound SNP pattern.
        // All discriminating positions (0 and 4) have high quality.
        let seq = b"CCAAGGGTTT";
        let qual = &[30, 30, 35, 38, 32, 36, 34, 30, 30, 30];
        let cigar = CigarString(vec![Cigar::Match(10)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "GAGGG", "AAGGA");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        assert!(
            matches!(result, MnpResult::ThirdAllele(..)),
            "Expected ThirdAllele for partial ALT match (only pos 0 mutated)"
        );
    }

    // ── Phase 1: Masked per-position evaluation tests ──

    #[test]
    fn test_mnp_partial_match_count_in_third_allele() {
        // TERT-like ONP: GAGGG→AAGGA (discriminating positions: 0 and 4).
        // Read carries AAGGG: pos 0 matches ALT (G→A), pos 4 matches REF (G, not A).
        // Expected: ThirdAllele with positions_matching_alt=1
        let seq = b"CCAAGGGTTT";
        let qual = &[30, 30, 35, 38, 32, 36, 34, 30, 30, 30];
        let cigar = CigarString(vec![Cigar::Match(10)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "GAGGG", "AAGGA");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        match result {
            MnpResult::ThirdAllele(partial, _) => assert_eq!(partial, 1,
                "Expected 1 position matching ALT (pos 0), got {}", partial),
            other => panic!("Expected ThirdAllele(1), got {:?}", format_mnp_result(&other)),
        }
    }

    #[test]
    fn test_mnp_all_masked_returns_low_quality_with_zero_partial() {
        // DNP: GG→AA. Read has AA (would match ALT) but both positions have Q < 20.
        // All discriminating positions masked → LowQuality.
        // positions_matching_alt only counts UNMASKED positions,
        // so when ALL positions are masked, partial=0 (no reliable evidence).
        let seq = b"CCAAGCC";
        let qual = &[30, 30, 5, 8, 30, 30, 30]; // both discriminating positions low-Q
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "GG", "AA");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        match result {
            MnpResult::LowQuality(partial, _) => assert_eq!(partial, 0,
                "All positions masked → no reliable ALT evidence, expected 0, got {}", partial),
            other => panic!("Expected LowQuality(0), got {:?}", format_mnp_result(&other)),
        }
    }

    #[test]
    fn test_mnp_n_base_at_one_discriminating_position_recovers_alt() {
        // TERT-like ONP: GAGGG→AAGGA (discriminating positions: 0 and 4).
        // Read carries NAGGA: pos 0 is N (masked), pos 4 matches ALT (A).
        // 1 unmasked position matches ALT → Alt (read recovered).
        let seq = b"CCNAGGATTTT";
        let qual = &[30, 30, 30, 38, 32, 36, 34, 30, 30, 30, 30];
        let cigar = CigarString(vec![Cigar::Match(11)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "GAGGG", "AAGGA");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        match result {
            MnpResult::Alt(q, had_n, confirmed) => {
                assert!(q > 0,
                    "Expected Alt when N masks one position but other unmasked matches ALT, got q={}", q);
                assert!(had_n,
                    "Expected had_n=true when N base present at discriminating position");
                assert!(!confirmed, "an N-masked discriminating base → not confirmed");
            }
            other => panic!(
                "Expected Alt for ONP with N at one discriminating position, got {:?}",
                format_mnp_result(&other)),
        }
    }

    #[test]
    fn test_mnp_n_base_at_all_discriminating_returns_low_quality() {
        // DNP: GG→AA. Read has NA (N at pos 0 masks it, pos 1 has low Q).
        // Both discriminating positions masked → LowQuality.
        let seq = b"CCNAGCC";
        let qual = &[30, 30, 30, 5, 30, 30, 30]; // pos 3: Q=5 < 20
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "GG", "AA");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        assert!(
            matches!(result, MnpResult::LowQuality(..)),
            "Expected LowQuality when all discriminating positions masked (N + low-Q)"
        );
    }

    #[test]
    fn test_mnp_all_n_high_bq_both_masked() {
        // Phase 0.6 plan: "DNP CC→TT: read carries NN with BQ=[30, 30].
        // Both positions masked → LowQuality (contributes to DP only)."
        //
        // This is distinct from test_mnp_n_base_at_all_discriminating_returns_low_quality
        // (which uses N + low-Q). Here BOTH positions are N with HIGH BQ, verifying
        // that the N guard masks independently of the BQ gate. If the N guard were
        // removed, this test would fail (high-BQ N would pass the BQ gate and
        // fall through as ThirdAllele), while the N+low-Q test would still pass.
        let seq = b"CCNNGCC";
        let qual = &[30, 30, 30, 30, 30, 30, 30]; // ALL high BQ — N masking is sole defense
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "CC", "TT");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        match result {
            MnpResult::LowQuality(partial, had_n) => {
                assert_eq!(partial, 0,
                    "All-N should have 0 partial ALT matches (no reliable evidence)");
                assert!(had_n,
                    "All-N with high BQ should report had_n=true for n_count tracking");
            }
            other => panic!(
                "Expected LowQuality for NN read with high BQ, got {:?}. \
                 If this returned ThirdAllele, the N guard may be missing.",
                format_mnp_result(&other)),
        }
    }

    #[test]
    fn test_classify_result_partial_match_count() {
        // Verify ClassifyResult::neither_with_partial carries partial count and has_n
        let r = ClassifyResult::neither_with_partial(ClassifyPhase::MaskedCompare, 3, false);
        assert!(!r.is_ref, "neither_with_partial should not be REF");
        assert!(!r.is_alt, "neither_with_partial should not be ALT");
        assert_eq!(r.partial_match_count, 3, "partial_match_count should be 3");
        assert!(!r.has_n_base, "has_n_base should be false when has_n=false");

        // Verify has_n=true is propagated
        let r_n = ClassifyResult::neither_with_partial(ClassifyPhase::MaskedCompare, 2, true);
        assert!(r_n.has_n_base, "has_n_base should be true when has_n=true");
        assert_eq!(r_n.partial_match_count, 2, "partial_match_count should be 2");

        // Verify standard constructors have partial_match_count == 0
        let r2 = ClassifyResult::is_alt(30, ClassifyPhase::Structural);
        assert_eq!(r2.partial_match_count, 0, "is_alt should have partial=0");
        let r3 = ClassifyResult::neither(ClassifyPhase::Structural);
        assert_eq!(r3.partial_match_count, 0, "neither should have partial=0");
    }

    #[test]
    fn test_base_counts_any_alt_invariant() {
        // Verify structural invariants of BaseCounts:
        //   1. any_alt = ad + partial_alt
        //   2. any_alt >= ad (partial_alt is non-negative)
        //   3. DP >= RD + AD + partial_alt + n_count (decomposition)
        //
        // These must hold regardless of variant type or classification path.
        use crate::types::BaseCounts;

        // ── Case 1: MNP-like scenario (mixed full + partial ALT) ──
        let mut counts = BaseCounts { dp: 20, ..Default::default() };

        // Simulate: 3 full ALT reads + 2 partial reads + 5 REF + 2 N + 8 other
        for _ in 0..3 {
            counts.ad += 1;
            counts.any_alt += 1; // Full ALT → any_alt++
        }
        for _ in 0..2 {
            counts.any_alt += 1;
            counts.partial_alt += 1; // Partial → any_alt++, partial_alt++
        }
        counts.rd = 5;
        counts.n_count = 2;

        // Invariant 1: any_alt = ad + partial_alt
        assert_eq!(
            counts.any_alt,
            counts.ad + counts.partial_alt,
            "Invariant 1 violated: any_alt({}) != ad({}) + partial_alt({})",
            counts.any_alt, counts.ad, counts.partial_alt
        );
        // Invariant 2: any_alt >= ad
        assert!(
            counts.any_alt >= counts.ad,
            "Invariant 2 violated: any_alt({}) < ad({})",
            counts.any_alt, counts.ad
        );
        // Invariant 3: DP >= RD + AD + partial_alt + n_count
        let decomposed = counts.rd + counts.ad + counts.partial_alt + counts.n_count;
        assert!(
            counts.dp >= decomposed,
            "Invariant 3 violated: DP({}) < RD({}) + AD({}) + partial_alt({}) + n_count({}) = {}",
            counts.dp, counts.rd, counts.ad, counts.partial_alt, counts.n_count, decomposed
        );
        assert_eq!(counts.any_alt, 5);
        assert_eq!(counts.ad, 3);
        assert_eq!(counts.partial_alt, 2);

        // ── Case 2: SNP/Indel (no partial concept) ──
        let snp_counts = BaseCounts {
            dp: 100, rd: 80, ad: 20, any_alt: 20,
            partial_alt: 0, n_count: 0, ..Default::default()
        };
        assert_eq!(snp_counts.any_alt, snp_counts.ad + snp_counts.partial_alt);
        assert!(snp_counts.any_alt >= snp_counts.ad);
        assert!(snp_counts.dp >= snp_counts.rd + snp_counts.ad + snp_counts.partial_alt + snp_counts.n_count);

        // ── Case 3: High N-count (duplex masking hotspot) ──
        let n_heavy = BaseCounts {
            dp: 50, rd: 10, ad: 2, any_alt: 2,
            partial_alt: 0, n_count: 30, ..Default::default()
        };
        assert_eq!(n_heavy.any_alt, n_heavy.ad + n_heavy.partial_alt);
        assert!(n_heavy.dp >= n_heavy.rd + n_heavy.ad + n_heavy.partial_alt + n_heavy.n_count);
        // n_count/DP ratio for QC
        let n_ratio = n_heavy.n_count as f64 / n_heavy.dp as f64;
        assert!(n_ratio > 0.5, "N-heavy site should have high n_count/DP ratio: {}", n_ratio);
    }

    // ── SW vs PairHMM concordance tests ──
    //
    // These tests send identical synthetic reads through both backends
    // and assert that unambiguous cases produce the same REF/ALT decision.
    // Uses check_allele_with_qual (the real dispatch entry point).

    fn build_variant_with_context(
        pos: i64, ref_allele: &str, alt_allele: &str,
        ref_context: &str, ctx_start: i64,
    ) -> Variant {
        Variant {
            chrom: "1".to_string(),
            pos,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            variant_type: String::new(),
            ref_context: Some(ref_context.to_string()),
            ref_context_start: ctx_start,
            repeat_span: 0,
            gene_strand: None,
            shift_region: None,
            event_ref: None,
            boundary_span: None,
        }
    }

    /// Run check_allele_with_qual through both SW and HMM backends for concordance testing.
    fn run_both_backends(
        record: &Record, variant: &Variant, min_baseq: u8,
    ) -> (ClassifyResult, ClassifyResult) {
        // We need separate aligner instances because check_allele_with_qual takes &mut
        let scoring_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

        let mut alt_a1 = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE)
                .yclip(0),
        );
        let mut ref_a1 = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE)
                .yclip(0),
        );

        let quals = record.qual();
        let sw_result = check_allele_with_qual(
            record, variant, &[], quals, min_baseq,
            &mut alt_a1, &mut ref_a1,
            &AlignmentBackend::SmithWaterman, &carrier::ReadRules::DNA,
        );
        let hmm_result = check_allele_with_qual(
            record, variant, &[], quals, min_baseq,
            &mut alt_a1, &mut ref_a1,
            &AlignmentBackend::pairhmm_default(), &carrier::ReadRules::DNA,
        );

        (sw_result, hmm_result)
    }

    #[test]
    fn test_concordance_snp_ref() {
        // Read matches REF haplotype at high quality — both backends should say REF
        let variant = build_variant_with_context(5, "C", "T", "GGGGGCGGGGG", 0);
        let cigar = CigarString(vec![rust_htslib::bam::record::Cigar::Match(11)]);
        let record = build_record(b"GGGGGCGGGGG", &[35_u8; 11], &cigar, 0);

        let (sw, hmm) = run_both_backends(&record, &variant, 20);
        assert_eq!((sw.is_ref, sw.is_alt), (true, false), "SW should classify as REF");
        assert_eq!((hmm.is_ref, hmm.is_alt), (true, false), "PairHMM should classify as REF");
    }

    #[test]
    fn test_concordance_snp_alt() {
        // Read matches ALT haplotype at high quality — both backends should say ALT
        let variant = build_variant_with_context(5, "C", "T", "GGGGGCGGGGG", 0);
        let cigar = CigarString(vec![rust_htslib::bam::record::Cigar::Match(11)]);
        let record = build_record(b"GGGGGTGGGGG", &[35_u8; 11], &cigar, 0);

        let (sw, hmm) = run_both_backends(&record, &variant, 20);
        assert_eq!((sw.is_ref, sw.is_alt), (false, true), "SW should classify as ALT");
        assert_eq!((hmm.is_ref, hmm.is_alt), (false, true), "PairHMM should classify as ALT");
    }

    #[test]
    fn test_concordance_insertion_alt() {
        // Read carries a 3bp insertion (A→ACCC) — both backends should agree ALT
        let variant = build_variant_with_context(4, "A", "ACCC", "GGGGAGGGGG", 0);
        let cigar = CigarString(vec![
            rust_htslib::bam::record::Cigar::Match(5),
            rust_htslib::bam::record::Cigar::Ins(3),
            rust_htslib::bam::record::Cigar::Match(5),
        ]);
        let record = build_record(b"GGGGACCCGGGGG", &[35_u8; 13], &cigar, 0);

        let (sw, hmm) = run_both_backends(&record, &variant, 20);
        assert!(sw.is_alt, "SW should classify insertion as ALT");
        assert!(hmm.is_alt, "PairHMM should classify insertion as ALT");
    }

    #[test]
    fn test_concordance_deletion_alt() {
        // Read carries a 3bp deletion (ACCC→A) — both backends should agree ALT
        let variant = build_variant_with_context(4, "ACCC", "A", "GGGGACCCGGGGG", 0);
        let cigar = CigarString(vec![
            rust_htslib::bam::record::Cigar::Match(5),
            rust_htslib::bam::record::Cigar::Del(3),
            rust_htslib::bam::record::Cigar::Match(5),
        ]);
        let record = build_record(b"GGGGAGGGGG", &[35_u8; 10], &cigar, 0);

        let (sw, hmm) = run_both_backends(&record, &variant, 20);
        assert!(sw.is_alt, "SW should classify deletion as ALT");
        assert!(hmm.is_alt, "PairHMM should classify deletion as ALT");
    }

    #[test]
    fn test_concordance_complex_delins() {
        // Complex variant: TC→GA (2bp substitution, same length) at pos 5
        // This is simpler than a DelIns — both haplotypes are the same length,
        // so both backends can classify reliably.
        // REF context: GGGGG TC GGGGG (13bp, variant at offset 5)
        // ALT read:    GGGGG GA GGGGG (matches ALT)
        let variant = build_variant_with_context(5, "TC", "GA", "GGGGGTCGGGGG", 0);
        let cigar = CigarString(vec![rust_htslib::bam::record::Cigar::Match(12)]);
        let record = build_record(b"GGGGGGAGGGGG", &[35_u8; 12], &cigar, 0);

        let (sw, hmm) = run_both_backends(&record, &variant, 20);
        // Both backends should agree on ALT for this unambiguous 2bp substitution
        assert!(sw.is_alt, "SW should classify 2bp sub as ALT, got is_ref={} is_alt={}", sw.is_ref, sw.is_alt);
        assert!(hmm.is_alt, "PairHMM should classify 2bp sub as ALT, got is_ref={} is_alt={}", hmm.is_ref, hmm.is_alt);
    }

    // ── n_count accumulation integration tests ──
    //
    // Verify that has_n_base flows through the classification dispatch
    // and would correctly drive counts.n_count += 1 in the engine.
    // These are component-level integration tests (dispatch → ClassifyResult)
    // since full engine tests require BAM file I/O infrastructure.

    #[test]
    fn test_n_count_snp_dispatch_sets_has_n_base() {
        // SNP A→T: read carries N at variant position with BQ=30.
        // check_allele_with_qual dispatches to check_snp → neither_n() → has_n_base=true.
        // The engine would then do: counts.n_count += 1.
        let seq = b"GGGGNGGGGG";
        let qual = &[35, 35, 35, 35, 30, 35, 35, 35, 35, 35]; // N at pos 4, BQ=30
        let cigar = CigarString(vec![Cigar::Match(10)]);
        let record = build_record(seq, qual, &cigar, 0);

        let variant = build_variant_with_context(4, "A", "T", "GGGGAGGGGG", 0);
        let scoring_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut alt_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let mut ref_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );

        let result = check_allele_with_qual(
            &record, &variant, &[], record.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman, &carrier::ReadRules::DNA,
        );

        assert!(!result.is_ref, "N at SNP should not be REF");
        assert!(!result.is_alt, "N at SNP should not be ALT");
        assert!(result.has_n_base, "N at SNP should set has_n_base=true for n_count accumulation");
    }

    #[test]
    fn test_n_count_mnp_dispatch_with_partial_n() {
        // MNP AT→CG: read has N at first discriminating pos, G at second.
        // check_allele_with_qual dispatches to check_mnp → Alt(qual, had_n=true).
        // The engine would then do: ad++, any_alt++, AND counts.n_count += 1.
        let seq = b"GGNGGGGGGG";
        //          pos: 0 1 2 3 4 5 6 7 8 9
        // Variant at pos 2-3: REF=AT, ALT=CG
        // Read has N at pos 2 (masked), G at pos 3 (matches ALT)
        let qual = &[35, 35, 30, 35, 35, 35, 35, 35, 35, 35]; // N at pos 2, BQ=30
        let cigar = CigarString(vec![Cigar::Match(10)]);
        let record = build_record(seq, qual, &cigar, 0);

        let variant = build_variant(2, "AT", "CG");
        let result = check_mnp(&record, &variant, record.qual(), 20);

        // N at first position is masked, G at second matches ALT → classified as ALT
        match result {
            MnpResult::Alt(_, had_n, confirmed) => {
                assert!(had_n, "MNP Alt with N at one position should report had_n=true");
                assert!(!confirmed, "an N-masked discriminating base → not confirmed");
            }
            other => panic!("Expected MnpResult::Alt, got {:?}", other),
        }
    }

    #[test]
    fn test_n_count_mnp_dispatch_propagates_through_check_allele() {
        // Same scenario as above but through check_allele_with_qual dispatch.
        // Verifies the MnpResult → ClassifyResult → has_n_base chain.
        let seq = b"GGNGGGGGGG";
        let qual = &[35, 35, 30, 35, 35, 35, 35, 35, 35, 35];
        let cigar = CigarString(vec![Cigar::Match(10)]);
        let record = build_record(seq, qual, &cigar, 0);

        let variant = build_variant_with_context(2, "AT", "CG", "GGATGGGGGG", 0);
        let scoring_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut alt_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let mut ref_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );

        let result = check_allele_with_qual(
            &record, &variant, &[], record.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman, &carrier::ReadRules::DNA,
        );

        assert!(result.is_alt, "MNP with N-masked + ALT-matching should be ALT");
        assert!(result.has_n_base, "MNP ALT with N at one position should propagate has_n_base=true");
        assert!(!result.mnp_confirmed, "an N-masked discriminating base → not confirmed");

        // Same read with the N replaced by the ALT base: every discriminating
        // base read → confirmed through the dispatch.
        let full = build_record(b"GGCGGGGGGG", qual, &cigar, 0);
        let result = check_allele_with_qual(
            &full, &variant, &[], full.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman, &carrier::ReadRules::DNA,
        );
        assert!(result.is_alt && result.mnp_confirmed, "fully read MNP ALT → confirmed");

        // An SNV ALT call through the same dispatch is never MNP-confirmed.
        let snv = build_variant_with_context(2, "A", "C", "GGATGGGGGG", 0);
        let result = check_allele_with_qual(
            &full, &snv, &[], full.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman, &carrier::ReadRules::DNA,
        );
        assert!(result.is_alt && !result.mnp_confirmed, "SNV ALT → never MNP-confirmed");
    }

    // ── check_complex N base propagation tests ──
    //
    // Verify that check_complex detects N in the reconstructed haplotype
    // and propagates has_n_base through all Phase 2 return paths.

    #[test]
    fn test_complex_n_in_reconstructed_haplotype_alt_match() {
        // Complex variant (equal-length, routes to Case A masked_dual_compare):
        // REF=TC, ALT=GA at pos 5. ref_len=2, alt_len=2 → same length.
        // HOWEVER: equal-length REF/ALT with len>1 routes to check_mnp, not check_complex.
        // To force check_complex, we need ref_len != alt_len.
        //
        // Strategy: REF=TCC (3bp) ALT=GAG (3bp) at pos 5.
        // WAIT: that's same-length → MNP. We need different lengths.
        //
        // Strategy: use check_complex directly (it's pub) with a same-length
        // complex variant. The engine dispatches MNPs to check_mnp, but
        // check_complex itself handles equal-length via Case A.
        // REF=TC, ALT=GA at pos 5, read has "GN" at pos 5-6.
        // Reconstruction: seq[5..7] = "GN" (2bp = ref_len = alt_len → Case A).
        // masked_dual_compare: N masked, G matches ALT[0] but not REF[0].
        // → mismatches_alt=0 (G matches ALT[0]), mismatches_ref=1 (G ≠ REF[0]='T')
        // → ALT match on 1 reliable base.
        let seq = b"GGGGGGNGGGG";
        //          pos: 0 1 2 3 4 5 6 7 8 9 10
        // Variant at pos 5: REF=TC (2bp), ALT=GA (2bp)
        // Reconstruction from pos 5-6: "GN" → Case A (both same length)
        let qual = &[35, 35, 35, 35, 35, 35, 30, 35, 35, 35, 35];
        let cigar = CigarString(vec![Cigar::Match(11)]);
        let record = build_record(seq, qual, &cigar, 0);

        // Call check_complex directly (bypassing MNP dispatch).
        // This tests the N detection and propagation in check_complex itself.
        let variant = build_variant_with_context(5, "TC", "GA", "GGGGGTCGGGGG", 0);
        let scoring_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut alt_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let mut ref_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );

        let result = check_complex(
            &record, &variant, &[], record.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman,
        );

        // The reconstruction "GN" has N at position 1. N is masked by
        // masked_dual_compare. Only position 0 (G) is reliable:
        // G matches ALT[0]='G' → mismatches_alt=0 → ALT.
        // has_n_base should be true because N was detected in the reconstruction.
        assert!(result.is_alt,
            "Complex with N-masked + ALT-matching reliable base should be ALT, \
             got is_ref={}, is_alt={}", result.is_ref, result.is_alt);
        assert!(result.has_n_base,
            "check_complex with N in reconstructed haplotype should set has_n_base=true");
    }

    #[test]
    fn test_complex_anchor_quality_skips_a_leading_hard_clip() {
        // A REF read over a long deletion-direction complex allele (100 REF bases
        // replaced by one), ending 8 bases past the anchor: too little of the REF
        // span for the masked comparison, so the clean-coverage path calls it REF
        // from the anchor's quality. The read is hard-clipped by 5 bases at its
        // start; hard-clipped bases are not in SEQ, so that quality is the one at
        // the anchor's own query position (Q30), not 5 bases on (Q5).
        let ref_allele = format!("G{}", "ACGT".repeat(24) + "ACG");
        let context = format!("{}{}{}", "A".repeat(10), ref_allele, "T".repeat(20));
        let seq = &context.as_bytes()[..18];
        let mut qual = [30u8; 18];
        qual[15] = 5;
        let cigar = CigarString(vec![Cigar::HardClip(5), Cigar::Match(18)]);
        let record = build_record(seq, &qual, &cigar, 0);
        let variant = build_variant_with_context(10, &ref_allele, "A", &context, 0);
        let scoring_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut alt_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let mut ref_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let result = check_complex(
            &record, &variant, &[], record.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman,
        );
        assert!(
            result.is_ref && result.qual == 30,
            "expected REF at the anchor's Q30, got is_ref={} qual={}", result.is_ref, result.qual
        );
    }

    #[test]
    fn test_complex_no_n_in_reconstructed_haplotype() {
        // Complex variant with no N bases — has_n_base should be false.
        // REF=TCC, ALT=GA at pos 5.
        // Read has "GA" at the variant position (matches ALT), no N bases.
        let seq = b"GGGGGGAGGGGG";
        let qual = &[35_u8; 12];
        let cigar = CigarString(vec![Cigar::Match(12)]);
        let record = build_record(seq, qual, &cigar, 0);

        let variant = build_variant_with_context(5, "TCC", "GA", "GGGGGTCCGGGGG", 0);
        let scoring_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut alt_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let mut ref_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );

        let result = check_complex(
            &record, &variant, &[], record.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman,
        );

        assert!(!result.has_n_base,
            "check_complex with no N in reconstructed haplotype should have has_n_base=false");
    }

    // ── INDEL wrong-length / nearby-evidence tests ──
    //
    // These tests verify the wrong-length INDEL handling added to fix
    // PAX5-class discordances (originally Phase-3 fallbacks; now the
    // wrong-length rule resolves lone ops directly as neither + partial
    // evidence, with Phase 3 kept for shifted same-length candidates and in-band
    // large deletions whose bases differ). Each test documents which code path in
    // check_insertion/check_deletion it exercises.

    // Aligner construction is inlined in each test below because:
    // 1. Rust's impl Trait creates distinct opaque types per call site
    // 2. The scoring closure must outlive the Aligner that borrows it
    // Pattern: let scoring_fn = ...; let mut alt_a = Aligner::with_capacity_and_scoring(...);
    // This matches the existing check_complex tests in this file.

    #[test]
    fn test_insertion_correct_length_seq_match() {
        // Regression test: read has I(2) matching ALT="ACC" at anchor.
        // Expected path: strict fast path → exact match → ALT.
        //
        // Geometry:
        //   Ref:  ...GGGGA----GGGGG...   (A at pos 10, no insertion in ref)
        //   Read: ...GGGGACCGGGGG        (I(2) = "CC" after anchor A)
        //   CIGAR: 5M 2I 5M
        //
        // Variant: pos=14, REF=A, ALT=ACC (2bp insertion after anchor)
        // ref_context covers positions 10-19: "GGGGAGGGGG"
        let seq = b"GGGGACCGGGGG";
        let qual = &[35_u8; 12];
        let cigar = CigarString(vec![Cigar::Match(5), Cigar::Ins(2), Cigar::Match(5)]);
        let record = build_record(seq, qual, &cigar, 10);

        // Insertion after pos 14 (anchor = last base of 5M block = 10+5-1 = 14)
        let variant = build_variant_with_context(14, "A", "ACC", "GGGGAGGGGG", 10);
        let scoring_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut alt_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let mut ref_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );

        let result = check_insertion(
            &record, &variant, &[], record.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman,
        );

        assert!(result.is_alt, "Correct-length insertion with matching sequence should be ALT");
        assert!(!result.is_ref, "Should not be REF");
        assert!(!result.has_nearby_evidence, "Exact match should not set has_nearby_evidence");
    }

    #[test]
    fn test_insertion_no_insertion_at_anchor() {
        // Regression test: read has only M blocks covering the anchor → REF.
        // Expected path: strict fast path → anchor in middle of M block →
        //   found_ref_coverage = true → REF.
        //
        // Geometry:
        //   Ref:  ...GGGGAGGGGG...
        //   Read: ...GGGGAGGGGG     (no insertion, matches ref)
        //   CIGAR: 10M
        //
        // Variant: pos=14, REF=A, ALT=ACC
        let seq = b"GGGGAGGGGG";
        let qual = &[35_u8; 10];
        let cigar = CigarString(vec![Cigar::Match(10)]);
        let record = build_record(seq, qual, &cigar, 10);

        let variant = build_variant_with_context(14, "A", "ACC", "GGGGAGGGGG", 10);
        let scoring_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut alt_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let mut ref_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );

        let result = check_insertion(
            &record, &variant, &[], record.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman,
        );

        assert!(result.is_ref, "Read with no insertion at anchor should be REF");
        assert!(!result.is_alt, "Should not be ALT");
        assert!(!result.has_nearby_evidence, "Clean REF should not set has_nearby_evidence");
    }

    #[test]
    fn test_insertion_wrong_length_at_anchor() {
        // PAX5-class test: read has I(1) at strict anchor but expected I(2).
        // Expected path: strict wrong-length rule — no truncation (observed
        // 1bp is below the containment floor), no split ops → lone wrong-
        // length I → neither + has_nearby_evidence, because I(1) at the
        // anchor is structural evidence of a third allele.
        //
        // Geometry:
        //   Ref:  ...GGGGA---GGGGG...   (anchor A at pos 14)
        //   Read: ...GGGGACGGGGG         (I(1) = "C" instead of expected "CC")
        //   CIGAR: 5M 1I 5M
        //
        // Variant: pos=14, REF=A, ALT=ACC (expected 2bp insertion)
        let seq = b"GGGGACGGGGG";
        let qual = &[35_u8; 11];
        let cigar = CigarString(vec![Cigar::Match(5), Cigar::Ins(1), Cigar::Match(5)]);
        let record = build_record(seq, qual, &cigar, 10);

        let variant = build_variant_with_context(14, "A", "ACC", "GGGGAGGGGG", 10);
        let scoring_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut alt_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let mut ref_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );

        let result = check_insertion(
            &record, &variant, &[], record.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman,
        );

        // Phase 3 may classify as REF (wrong allele) or ALT (if SW finds
        // partial support). Either way, has_nearby_evidence must be true
        // because there IS an insertion at the anchor position.
        assert!(
            result.has_nearby_evidence || result.is_alt,
            "Wrong-length I(1) at anchor for expected I(2) must either be ALT \
             or have has_nearby_evidence=true (for partial_alt counting). \
             Got is_ref={}, is_alt={}, has_nearby_evidence={}",
            result.is_ref, result.is_alt, result.has_nearby_evidence
        );
    }

    #[test]
    fn test_insertion_same_length_wrong_sequence() {
        // Read has I(2) at anchor but with wrong bases ("TT" vs expected "CC")
        // on confident (≥min_baseq) bases: a same-length THIRD allele.
        //
        // Geometry:
        //   Ref:  ...GGGGA---GGGGG...
        //   Read: ...GGGGATTGGGGG       (I(2) = "TT" instead of "CC")
        //   CIGAR: 5M 2I 5M
        //
        // Expected path: strict fast path → length matches → confident seq
        // mismatch → neither + has_nearby_evidence. Not REF (rd must not
        // absorb a read that provably carries an insertion), not ALT (the
        // bases are provably not the queried insert), and not Phase 3 —
        // alignment scoring would promote the wrong-sequence insert to ALT
        // because it still beats the gapped REF alignment.
        let seq = b"GGGGATTGGGGG";
        let qual = &[35_u8; 12];
        let cigar = CigarString(vec![Cigar::Match(5), Cigar::Ins(2), Cigar::Match(5)]);
        let record = build_record(seq, qual, &cigar, 10);

        let variant = build_variant_with_context(14, "A", "ACC", "GGGGAGGGGG", 10);
        let scoring_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut alt_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let mut ref_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );

        let result = check_insertion(
            &record, &variant, &[], record.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman,
        );

        // A confident wrong-sequence insertion is a distinct allele:
        // neither REF nor ALT, with partial evidence for any_alt/partial_alt.
        assert!(!result.is_alt,
            "Same-length insertion with wrong sequence must not be ALT");
        assert!(!result.is_ref,
            "Same-length insertion with wrong sequence must not be absorbed into REF");
        assert!(result.has_nearby_evidence,
            "Same-length wrong-sequence insertion must carry partial evidence");
    }

    #[test]
    fn test_insertion_wrong_length_in_window() {
        // Read has I(1) at a windowed position (2bp from anchor) for expected I(2),
        // outside the discrimination window in unique context: a separate event,
        // so the read's bases across the window are REF (operator, 2026-10-01).
        //
        // Geometry:
        //   Ref:    ...GGGGGAGGGGG...   (anchor A at pos 15)
        //   Read:   ...GGGCGGAGGGGG     (I(1)="C" at ref pos 13, 2bp before anchor)
        //   CIGAR: 3M 1I 2M ...
        //   Anchor at pos 15 is in the second M block (pos 13..15), so
        //   found_ref_coverage = true. Insertion at pos 13 is in window
        //   [15-5, 15+5] = [10, 20]. Step 1.3's else clause sets
        //   has_nearby_length_match = true.
        //
        let seq = b"GGGCGGAGGGGG";
        let qual = &[35_u8; 12];
        // 3M at pos 10 → covers 10,11,12
        // 1I (1bp insertion at ref pos 13)
        // 8M at pos 13 → covers 13..21 (includes anchor at 15)
        let cigar = CigarString(vec![Cigar::Match(3), Cigar::Ins(1), Cigar::Match(8)]);
        let record = build_record(seq, qual, &cigar, 10);

        // Variant: 2bp insertion expected after pos 15
        let variant = build_variant_with_context(15, "A", "ACC", "GGGGGAGGGGG", 10);
        let scoring_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut alt_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let mut ref_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );

        let result = check_insertion(
            &record, &variant, &[], record.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman,
        );

        assert!(
            result.is_ref && !result.has_nearby_evidence,
            "An I(1) outside the window is a separate event: REF, no partial evidence. \
             Got is_ref={}, is_alt={}, has_nearby_evidence={}",
            result.is_ref, result.is_alt, result.has_nearby_evidence
        );
    }

    #[test]
    fn test_insertion_no_ref_coverage_spans_anchor() {
        // Read has no insertion and is soft-clipped at the anchor position,
        // so found_ref_coverage is never set. The !found_ref_coverage path
        // (Step 1.4) should invoke Phase 3 because the read spans the anchor.
        //
        // Geometry:
        //   Read starts at pos 10, CIGAR: 4M 3S
        //   Covers ref positions 10-13 (4M), then 3bp soft-clipped
        //   Anchor at pos 14 — read's ref end is 14, so ref_end (14) > anchor_pos (14)?
        //   No: 14 > 14 is false. Let's adjust.
        //   Use: 5M 3S at pos 10 → covers 10-14, ref_end=15
        //   Anchor at 14 → anchor is the last M base → BUT that means anchor_pos == block_end - 1
        //   which IS the strict path. Since no I follows (next op is S), it sets
        //   found_ref_coverage = true. Not what we want.
        //
        //   Better: Use a CIGAR that doesn't cover the anchor in M blocks.
        //   3M 2I 3S at pos 12 → M covers 12,13,14 → actually this covers anchor.
        //
        //   Simplest: read starts at pos 10, CIGAR: 3M 1D 1M 3S
        //   M covers 10-12 (3bp), D skips 13, M covers 14, then S.
        //   At the M(1) block starting at ref_pos=14: anchor_pos=14 is within [14,15).
        //   anchor_pos == block_end - 1 (14 == 14). Check next op → S(3).
        //   S is not I, so falls to line 1175: found_ref_coverage = true. Still triggers.
        //
        //   The !found_ref_coverage path fires when the read has NO M block covering
        //   the anchor. This happens with unusual CIGAR geometry like:
        //   M blocks end before the anchor, but ref_pos reaches past anchor via D ops.
        //
        //   Use: 3M 1I 2D at pos 10 → M covers 10-12, I (point), D skips 13-14.
        //   The M block ends at ref_pos=13 (block_end=13), anchor at 14 → not covered.
        //   After I: ref_pos still 13. After D(2): ref_pos=15.
        //   found_ref_coverage stays false. ref_end = 10 + 3(M) + 2(D) = 15.
        //   15 > 14 → spans anchor → Phase 3.
        //   BUT: this read doesn't have another M after the D, so it's truncated.
        //   We need more sequence. Add more M at the end:
        //   3M 1I 2D 3M at pos 10 → M covers 10-12, D skips 13-14, M covers 15-17.
        //   Anchor at 14: not in any M block → found_ref_coverage stays false.
        //   ref_end = 10 + 3 + 2 + 3 = 18. Spans anchor. Phase 3 invoked.
        let seq = b"GGGCGGG"; // 3M + 1I + 3M = 7 read bases
        let qual = &[35_u8; 7];
        let cigar = CigarString(vec![
            Cigar::Match(3), Cigar::Ins(1), Cigar::Del(2), Cigar::Match(3)
        ]);
        let record = build_record(seq, qual, &cigar, 10);

        // Anchor at 14: not covered by any M block
        let variant = build_variant_with_context(14, "A", "ACC", "GGGGGAGGGGG", 10);
        let scoring_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut alt_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let mut ref_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );

        let result = check_insertion(
            &record, &variant, &[], record.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman,
        );

        // Phase 3 is invoked (no CIGAR evidence, read spans anchor).
        // Should NOT be silently classified as neither without Phase 3.
        // Result depends on Phase 3, but should not be a silent drop.
        assert!(
            result.is_ref || result.is_alt || result.has_nearby_evidence,
            "Read spanning anchor without M coverage should invoke Phase 3, \
             not silently return neither. Got is_ref={}, is_alt={}, has_nearby_evidence={}",
            result.is_ref, result.is_alt, result.has_nearby_evidence
        );
    }

    #[test]
    fn test_deletion_wrong_length_windowed_small() {
        // Read has D(6) at a windowed position (NOT the strict anchor+1 position)
        // when expected D(10). Both are ≥5bp so Step 2.1 should flag
        // has_nearby_length_match.
        //
        // Geometry:
        //   Anchor at pos 15 (A), expected D(10) starts at pos 16 (strict).
        //   Read CIGAR: 7M 6D 7M at pos 10.
        //   7M covers pos 10-16 (includes anchor at 15). block_end=17.
        //   D(6) at ref pos 17 (del_ref_pos = block_end = 17).
        //   Strict position: anchor_pos + 1 = 16 ≠ 17 → NOT skipped.
        //   Window: [15-5, 15+5] = [10, 20]. 17 is in [10, 20] → windowed.
        //   del_len=6 ≠ expected_del_len=10, expected < 50 → Step 2.1.
        //   6 ≥ 5 → has_nearby_length_match = true.
        //
        // Post-walk: has_nearby_length_match && found_ref_coverage → Phase 3.
        let seq = b"GGGGGAGXXXXXXX"; // 7M + 7M = 14 read bases
        let qual = &[35_u8; 14];
        let cigar = CigarString(vec![Cigar::Match(7), Cigar::Del(6), Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 10);

        // Variant: 10bp deletion after anchor at pos 15
        // REF = "AXXXXXXXXXX" (anchor + 10 deleted), ALT = "A"
        // expected_del_len = 10. Strict position = anchor_pos + 1 = 16.
        let variant = build_variant_with_context(
            15, "AXXXXXXXXXX", "A",
            "GGGGGAXXXXXXXXXXGGGGG", 10,
        );
        let scoring_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut alt_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let mut ref_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(SW_GAP_OPEN, SW_GAP_EXTEND, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );

        let result = check_deletion(
            &record, &variant, &[], record.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman,
        );

        // D(6) at pos 17 is in window, 6 ≥ 5bp → has_nearby_length_match → Phase 3.
        // Phase 3 may return REF or ALT. If REF, has_nearby_evidence must be set.
        assert!(
            result.has_nearby_evidence || result.is_alt,
            "Wrong-length D(6) in window for expected D(10) (both ≥5bp) must \
             either be ALT or have has_nearby_evidence=true. \
             Got is_ref={}, is_alt={}, has_nearby_evidence={}",
            result.is_ref, result.is_alt, result.has_nearby_evidence
        );
    }
}
