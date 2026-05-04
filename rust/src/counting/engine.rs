//! Counting engine: BAM read classification → allele counts.
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

use pyo3::prelude::*;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{self, Read, Record};
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};

use crate::shared::stats::fisher_strand_bias;
use crate::types::{BaseCounts, Variant};

use rayon::prelude::*;

use anyhow::{Context, Result};
use log::{debug, info, trace};
use bio::alignment::pairwise::Aligner;

use super::fragment::{FragmentEvidence, hash_qname, hash_molecule};
use super::pairhmm::dynamic_sw_gap_extend;
use super::variant_checks::{check_snp, check_mnp, check_complex, check_insertion, check_deletion, MnpResult};
use super::utils::{find_read_pos, ClassifyResult, ClassifyPhase};
use super::mfsd;
use super::rna;
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
// Bins are split when either the genomic distance or variant count exceeds
// the limits. This matches the original C++ GBCMS architecture which used
// --max_block_size=200 and --max_block_dist=10000.

/// Default bin window size in base pairs.
const BIN_WINDOW: i64 = 10_000;

/// Maximum number of variants per bin.
///
/// When exceeded, the bin is split to prevent O(V × R) blowup in the
/// shared-read classification loop. Matches the original C++ GBCMS default.
/// Not a CLI flag — this is an internal performance constant that does not
/// affect output (validated by D1 parity tests).
const BIN_MAX_VARIANTS: usize = 200;

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

/// Build genomic bins from a list of variants.
///
/// Variants are grouped by chromosome and position into `BIN_WINDOW`-sized
/// bins. The algorithm is O(n) — a single pass over sorted variants.
/// Bins are padded by the maximum `repeat_span` + 5bp to ensure reads
/// at bin edges are captured.
///
/// # Arguments
/// * `variants` — Variants (need not be sorted; sorting is internal)
/// * `bam_header` — BAM header for chromosome → tid lookup
/// * `window` — Bin window size in bp (default: 10,000)
fn build_genomic_bins(
    variants: &[Variant],
    bam_header: &bam::HeaderView,
    window: i64,
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

    let mut bins: Vec<GenomicBin> = Vec::new();
    let mut i = 0;

    while i < sorted_indices.len() {
        let first_idx = sorted_indices[i];
        let chrom = &variants[first_idx].chrom;
        let tid = match bam_header.tid(chrom.as_bytes()) {
            Some(t) => t,
            None => {
                // Skip variants on chromosomes not in BAM
                debug!("Skipping variant on chromosome not in BAM: {}", chrom);
                i += 1;
                continue;
            }
        };

        let bin_start = variants[first_idx].pos;
        let mut bin_end = bin_start + window;
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
            if indices.len() >= BIN_MAX_VARIANTS {
                debug!(
                    "Bin split: {} variants reached BIN_MAX_VARIANTS={} at {}:{}",
                    indices.len(), BIN_MAX_VARIANTS,
                    chrom, variants[jdx].pos + 1,
                );
                break;
            }
            indices.push(jdx);
            max_repeat_span = max_repeat_span.max(variants[jdx].repeat_span as i64);
            // Extend bin end to cover this variant's ref span
            let var_end = variants[jdx].pos + variants[jdx].ref_allele.len() as i64;
            bin_end = bin_end.max(var_end + window / 2);
            j += 1;
        }

        // Pad bin by max repeat_span + 5bp to capture shifted reads
        // (same logic as count_single_variant's window_pad)
        let padding = std::cmp::max(5, max_repeat_span + 2);
        bins.push(GenomicBin {
            tid,
            start: (bin_start - padding).max(0),
            end: bin_end + padding,
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
/// Selectable via `--alignment-backend` CLI flag.
#[derive(Clone, Debug, Default, PartialEq)]
pub enum AlignmentBackend {
    /// Smith-Waterman with affine gap penalties (default for v2.8.0).
    /// Score-margin classification: (alt_score - ref_score) > 0 → ALT.
    #[default]
    SmithWaterman,
    /// PairHMM with BQ-aware emissions (planned default for v3.0.0).
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
    /// Create PairHMM backend with default parameters.
    ///
    /// Convenience constructor for tests and downstream consumers.
    /// In production, the Python CLI passes params directly to the
    /// `PairHMM { ... }` variant via `count_bam()`.
    #[allow(dead_code)]
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


/// Count bases for a list of variants in a BAM file.
///
/// When `decomposed` is provided (same length as `variants`), variants with
/// a `Some(decomposed_variant)` are counted twice — once with the original
/// allele and once with the corrected allele. The result with the higher
/// `ad` (alt_count) is returned, with `used_decomposed` set accordingly.
///
/// // INTENTIONAL: This per-variant codepath is retained alongside count_bam_binned()
/// // for parity testing. Both codepaths must produce identical BaseCounts for the
/// // same inputs. Once parity is confirmed across the 22-BAM regression suite,
/// // pipeline.py can switch to the binned codepath.
#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (bam_path, variants, decomposed, min_mapq, min_baseq, filter_duplicates, filter_secondary, filter_supplementary, filter_qc_failed, filter_improper_pair, filter_indel, threads, fragment_qual_threshold=10, sibling_variants=Vec::new(), alignment_backend="sw", hmm_llr_threshold=2.3, hmm_gap_open=1e-4, hmm_gap_extend=0.1, hmm_gap_open_repeat=1e-2, hmm_gap_extend_repeat=0.5, mode="dna", enforce_strandedness=false))]
pub fn count_bam(
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
    mode: &str,
    enforce_strandedness: bool,
) -> PyResult<Vec<BaseCounts>> {
    // Parse alignment backend from string
    let backend = match alignment_backend {
        "hmm" | "pairhmm" => AlignmentBackend::PairHMM {
            llr_threshold: hmm_llr_threshold,
            gap_open: hmm_gap_open,
            gap_extend: hmm_gap_extend,
            gap_open_repeat: hmm_gap_open_repeat,
            gap_extend_repeat: hmm_gap_extend_repeat,
        },
        _ => AlignmentBackend::SmithWaterman,
    };

    // We cannot share a single IndexedReader across threads because it's not Sync.
    // Instead, we use rayon's map_init to initialize a reader for each thread.
    // This is efficient because map_init reuses the thread-local state (the reader)
    // for multiple items processed by that thread.

    // Configure thread pool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Failed to build thread pool: {}", e)))?;

    // Pad sibling_variants to match variants length (handles default empty case)
    let mut sibling_variants = sibling_variants;
    let n = variants.len();
    sibling_variants.resize_with(n, Vec::new);

    // Zip variants with their decomposed counterparts and sibling alts for parallel iteration
    let paired: Vec<_> = variants.into_iter()
        .zip(decomposed)
        .zip(sibling_variants)
        .map(|((v, d), s)| (v, d, s))
        .collect();

    // Release GIL for parallel execution
    #[allow(deprecated)]
    let results: Result<Vec<BaseCounts>, anyhow::Error> = py.allow_threads(move || {
        pool.install(|| {
            paired
                .par_iter()
                .map_init(
                    || {
                        // Initialize thread-local BAM reader
                        bam::IndexedReader::from_path(&bam_path).map_err(|e| {
                            anyhow::anyhow!("Failed to open BAM: {}", e)
                        })
                    },
                    |bam_result, (variant, decomp_opt, siblings)| {
                        // Get the reader or return error if initialization failed
                        let bam = match bam_result {
                            Ok(b) => b,
                            Err(e) => return Err(anyhow::anyhow!("BAM init failed: {}", e)),
                        };

                        let counts_orig = count_single_variant(
                            bam,
                            variant,
                            siblings,
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
                            false,  // apply_baq: legacy codepath
                            None,   // umi_tag: legacy codepath
                            mode,
                            enforce_strandedness,
                        )?;

                        // Dual-count: if a decomposed variant exists, count it too
                        // and return whichever has the higher alt_count.
                        if let Some(decomp) = decomp_opt {
                            let counts_decomp = count_single_variant(
                                bam,
                                decomp,
                                siblings,
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
                                false,  // apply_baq: legacy codepath
                                None,   // umi_tag: legacy codepath
                                mode,
                                enforce_strandedness,
                            )?;

                            if counts_decomp.ad > counts_orig.ad {
                                // Sanity: both hypotheses count the same reads at the
                                // same locus, so DP should be nearly identical. A large
                                // divergence indicates a counting bug.
                                if (counts_decomp.dp as i64 - counts_orig.dp as i64).abs() > 2 {
                                    log::warn!(
                                        "DP mismatch in dual-counting: decomp={} orig={} at {}:{} {}→{}",
                                        counts_decomp.dp, counts_orig.dp,
                                        variant.chrom, variant.pos + 1,
                                        variant.ref_allele, decomp.alt_allele
                                    );
                                }
                                debug!(
                                    "Homopolymer decomp: corrected allele wins \
                                     (ad={} vs orig ad={}, dp_decomp={}, dp_orig={}) \
                                     for {}:{} {}→{}",
                                    counts_decomp.ad, counts_orig.ad,
                                    counts_decomp.dp, counts_orig.dp,
                                    variant.chrom, variant.pos + 1,
                                    variant.ref_allele, decomp.alt_allele,
                                );
                                return Ok(BaseCounts {
                                    used_decomposed: true,
                                    ..counts_decomp
                                });
                            }
                        }

                        Ok(counts_orig)
                    },
                )
                .collect()
        })
    });

    // Map anyhow::Error back to PyErr
    match results {
        Ok(r) => Ok(r),
        Err(e) => Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e))),
    }
}


/// Bin-centric parallel BAM counting with BAQ and UMI support.
///
/// Groups variants into ~10kb genomic bins, fetches reads once per bin,
/// then classifies each read against all variants in the bin. This reduces
/// I/O overhead vs per-variant `bam.fetch()` in `count_bam()`, especially
/// for MAF files with clustered variants.
///
/// New features over `count_bam()`:
/// - `apply_baq`: Heuristic BAQ quality downgrade near indels (both modes).
/// - `umi_tag`: UMI-aware fragment grouping via `hash_molecule()`.
///
/// // INTENTIONAL: This is the new codepath. count_bam() is retained for
/// // parity testing until the 22-BAM regression confirms identical counts.
#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (bam_path, variants, decomposed, min_mapq, min_baseq, filter_duplicates, filter_secondary, filter_supplementary, filter_qc_failed, filter_improper_pair, filter_indel, threads, fragment_qual_threshold=10, sibling_variants=Vec::new(), alignment_backend="sw", hmm_llr_threshold=2.3, hmm_gap_open=1e-4, hmm_gap_extend=0.1, hmm_gap_open_repeat=1e-2, hmm_gap_extend_repeat=0.5, apply_baq=false, umi_tag=None, mode="dna", enforce_strandedness=false, rna_editing_db=None))]
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
    rna_editing_db: Option<&str>,
) -> PyResult<Vec<BaseCounts>> {
    // Parse alignment backend from string
    let backend = match alignment_backend {
        "hmm" | "pairhmm" => AlignmentBackend::PairHMM {
            llr_threshold: hmm_llr_threshold,
            gap_open: hmm_gap_open,
            gap_extend: hmm_gap_extend,
            gap_open_repeat: hmm_gap_open_repeat,
            gap_extend_repeat: hmm_gap_extend_repeat,
        },
        _ => AlignmentBackend::SmithWaterman,
    };

    // ── D7: Load RNA editing site database (if provided) ──
    // Loaded ONCE at init, then shared across all bins/threads via Arc.
    // DB-only strategy: flag = True only when site is in REDIportal.
    // No DB = no editing flags (no pattern-matching guessing).
    let editing_sites: Option<HashSet<(String, i64)>> = match rna_editing_db {
        Some(path) => {
            let sites = rna::build_rna_editing_set(path)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(
                    format!("Failed to load RNA editing DB: {}", e)
                ))?;
            info!("D7: Loaded {} RNA editing sites from {}", sites.len(), path);
            Some(sites)
        }
        None => None,
    };
    // Wrap in Arc for thread-safe sharing across rayon workers
    let editing_sites = std::sync::Arc::new(editing_sites);

    if mode == "rna" {
        info!(
            "count_bam_binned: {} variants, apply_baq={}, umi_tag={:?}, backend={:?}, rna_editing_db={}",
            variants.len(), apply_baq, umi_tag, alignment_backend,
            rna_editing_db.unwrap_or("none"),
        );
    } else {
        info!(
            "count_bam_binned: {} variants, apply_baq={}, umi_tag={:?}, backend={:?}",
            variants.len(), apply_baq, umi_tag, alignment_backend,
        );
    }

    // Build genomic bins using a temporary BAM reader for header access
    let header_reader = bam::IndexedReader::from_path(&bam_path).map_err(|e| {
        pyo3::exceptions::PyRuntimeError::new_err(format!("Failed to open BAM: {}", e))
    })?;
    let bins = build_genomic_bins(&variants, header_reader.header(), BIN_WINDOW);

    // Pre-pad sibling_variants
    let mut sibling_variants = sibling_variants;
    let n = variants.len();
    sibling_variants.resize_with(n, Vec::new);

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

    // Configure thread pool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
            format!("Failed to build thread pool: {}", e)
        ))?;

    // Result array: one BaseCounts per variant, initialized to default
    let mut all_counts: Vec<BaseCounts> = (0..n).map(|_| BaseCounts::default()).collect();

    // Process bins in parallel, each bin does one bam.fetch()
    #[allow(deprecated)]
    let bin_results: Result<Vec<(usize, BaseCounts)>, anyhow::Error> = py.allow_threads(move || {
        pool.install(|| {
            bins.par_iter()
                .map_init(
                    || {
                        bam::IndexedReader::from_path(&bam_path).map_err(|e| {
                            anyhow::anyhow!("Failed to open BAM: {}", e)
                        })
                    },
                    |bam_result, bin| {
                        let bam = match bam_result {
                            Ok(b) => b,
                            Err(e) => return Err(anyhow::anyhow!("BAM init failed: {}", e)),
                        };

                        debug!(
                            "Processing bin tid={} {}-{} ({} variants)",
                            bin.tid, bin.start, bin.end,
                            bin.variant_indices.len(),
                        );

                        // ── D10: Shared-read optimization ──
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
                            &editing_sites,
                        )?)
                    },
                )
                .try_reduce(
                    Vec::new,
                    |mut acc, batch| {
                        acc.extend(batch);
                        Ok(acc)
                    },
                )
        })
    });

    // Scatter results back to variant-order array
    match bin_results {
        Ok(pairs) => {
            for (vi, counts) in pairs {
                all_counts[vi] = counts;
            }
            Ok(all_counts)
        }
        Err(e) => Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e))),
    }
}

/// Compute the reference-consumed end position of an aligned record.
///
/// Walks the CIGAR string summing reference-consuming operations (M/=/X/D/N).
/// This is equivalent to `record.cigar().end_pos()` but avoids the costly
/// CigarString allocation that `end_pos()` performs internally.
///
/// Used in both `count_bin_shared` and `count_single_variant` to determine
/// anchor overlap. Extracted as a helper to eliminate the duplicated inline
/// CIGAR walk that previously existed in `count_single_variant`.
#[inline]
fn read_ref_end(record: &Record) -> i64 {
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
}


// ── Shared-Read Bin Processing ─────────────────────────────────────────────
//
// D10: Port of the original C++ GBCMS block-processing architecture.
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
/// then classifies each read against each variant in the bin. This is the
/// D10 shared-read optimization.
///
/// # Arguments
///
/// * `bam` — Indexed BAM reader (mutable for fetch)
/// * `bin` — Genomic bin with variant indices
/// * `variants` — All variants (bin.variant_indices indexes into this)
/// * `decomposed` — Parallel array of decomposed variants for dual-counting
/// * `sibling_variants` — Per-variant sibling arrays for multi-allelic guard
/// * All filter/config params — same as count_single_variant
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
    editing_sites: &Option<HashSet<(String, i64)>>,
) -> Result<Vec<(usize, BaseCounts)>> {

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

    for result in bam.records() {
        let record = result.context("count_bin_shared: error reading BAM record")?;

        // Universal flag filters (delegated to shared::filters::ReadFilter).
        if !read_filter.passes(&record, &mut filter_counts) {
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

        read_cache.push(record);
    }

    debug!(
        "Bin tid={} {}-{}: {} reads cached ({} filtered: dup={} sec={} supp={} qc={} pair={} indel={} mapq={})",
        bin.tid, bin.start, bin.end, read_cache.len(),
        filter_counts.total() + mapq_filtered,
        filter_counts.duplicates, filter_counts.secondary, filter_counts.supplementary,
        filter_counts.qc_failed, filter_counts.improper_pair, filter_counts.indel,
        mapq_filtered,
    );

    // ══════════════════════════════════════════════════════════════════════
    // PHASE 1: Per-variant classification from cached reads
    // ══════════════════════════════════════════════════════════════════════
    //
    // For each variant (and its decomposed twin), iterate the read cache
    // and perform the same classification logic as count_single_variant.
    // Per-variant filters (strandedness, anchor overlap) are applied here.

    let mut results = Vec::with_capacity(bin.variant_indices.len());

    for &vi in &bin.variant_indices {
        let variant = &variants[vi];
        let siblings = &sibling_variants[vi];

        let counts_orig = count_variant_from_cache(
            &read_cache, variant, siblings,
            min_mapq, min_baseq,
            filter_improper_pair, filter_indel,
            fragment_qual_threshold, backend,
            apply_baq, umi_tag, mode, enforce_strandedness,
            editing_sites,
        )?;

        // Dual-count for decomposed variants: run the same classification
        // against the decomposed variant form and take the higher ALT count.
        let final_counts = if let Some(ref decomp) = decomposed[vi] {
            let counts_decomp = count_variant_from_cache(
                &read_cache, decomp, siblings,
                min_mapq, min_baseq,
                filter_improper_pair, filter_indel,
                fragment_qual_threshold, backend,
                apply_baq, umi_tag, mode, enforce_strandedness,
                editing_sites,
            )?;

            if counts_decomp.ad > counts_orig.ad {
                BaseCounts {
                    used_decomposed: true,
                    ..counts_decomp
                }
            } else {
                counts_orig
            }
        } else {
            counts_orig
        };

        results.push((vi, final_counts));
    }

    Ok(results)
}


/// Classify and count reads from a pre-fetched cache for a single variant.
///
/// This is the Phase 1 workhorse of D10. It performs the same logic as
/// `count_single_variant` but operates on `&[Record]` instead of doing a
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
    editing_sites: &Option<HashSet<(String, i64)>>,
) -> Result<BaseCounts> {

    let mut counts = BaseCounts::default();

    // Fragment tracking: QNAME hash -> FragmentEvidence
    let mut fragments: HashMap<u64, FragmentEvidence> = HashMap::new();
    let qual_diff_threshold: u8 = fragment_qual_threshold;

    // Distance-to-read-end tracking for QC metrics
    let mut alt_dists: Vec<u32> = Vec::with_capacity(500);
    let mut ref_dists: Vec<u32> = Vec::with_capacity(500);

    // Create SW aligners ONCE per variant (indelpost pattern).
    // Score function: N in read or haplotype → 0 (neutral, uninformative).
    // Matches GATK approach: N contributes no evidence for match or mismatch.
    // Handles N-in-read (duplex masking / sequencer failure) and N-in-haplotype (rare).
    let score_fn = |a: u8, b: u8| -> i32 {
        if a == b'N' || b == b'N' { 0 } else if a == b { 1 } else { -1 }
    };
    let gap_open: i32 = -5;
    let gap_extend: i32 = dynamic_sw_gap_extend(variant.repeat_span);
    let mut alt_aligner = Aligner::new(gap_open, gap_extend, &score_fn);
    let mut ref_aligner = Aligner::new(gap_open, gap_extend, &score_fn);

    // Per-phase classification counters
    let mut phase_counts = [0u32; 5];

    // Compute the variant fetch window — same logic as count_single_variant
    // used for determining which cached reads overlap this variant
    let window_pad: i64 = std::cmp::max(5, variant.repeat_span as i64 + 2);
    let v_start = (variant.pos - window_pad).max(0);
    let v_end = variant.pos + (variant.ref_allele.len() as i64) + window_pad;

    // ── D6: RNA CONSENSUS SPLICING ──────────────────────────────────────
    // For RNA mode, snip consensus introns from the variant's ref_context
    // so that haplotype alignment in Phase 3 uses a mature-mRNA-like
    // reference instead of genomic (intron-containing) sequence.
    //
    // The read cache provides local reads — only reads overlapping this
    // variant's window contribute their CIGAR N ops for intron discovery.
    // apply_consensus_splicing uses >50% consensus threshold to filter
    // alignment artifacts and rare alternative splicing.
    //
    // Uses Cow to avoid cloning Variant when no introns are found
    // (common case: >80% of variants have no overlapping introns).
    let effective_variant: Cow<'_, Variant>;
    if mode == "rna" {
        if let Some(ref ctx) = variant.ref_context {
            let local_reads: Vec<&Record> = read_cache.iter()
                .filter(|r| r.pos() < v_end && read_ref_end(r) > v_start)
                .collect();

            let spliced_ctx = rna::apply_consensus_splicing(
                ctx.as_bytes(), &local_reads, variant.ref_context_start,
            );

            if spliced_ctx.len() != ctx.len() {
                // Introns were snipped — create a modified variant copy
                let mut v2 = variant.clone();
                v2.ref_context = Some(String::from_utf8_lossy(&spliced_ctx).into_owned());
                debug!(
                    "D6 splice: {}:{} ref_context {} → {} bases ({} intron bases removed)",
                    variant.chrom, variant.pos + 1,
                    ctx.len(), spliced_ctx.len(),
                    ctx.len() - spliced_ctx.len(),
                );
                effective_variant = Cow::Owned(v2);
            } else {
                effective_variant = Cow::Borrowed(variant);
            }
        } else {
            effective_variant = Cow::Borrowed(variant);
        }
    } else {
        effective_variant = Cow::Borrowed(variant);
    }
    // Shadow `variant` with the potentially-spliced version.
    // All downstream code (Phase 3 classification, haplotype matrix
    // construction) automatically uses the mature mRNA ref_context.
    let variant = effective_variant.as_ref();

    let mut reads_considered = 0u32;

    for record in read_cache {
        // ── Per-variant overlap check: does this read overlap the variant window?
        // This is the key filter that ensures each variant only sees reads
        // that would have been fetched by a per-variant bam.fetch().
        let r_start = record.pos();
        let r_end = read_ref_end(record);
        if r_start >= v_end || r_end <= v_start {
            continue; // Read doesn't overlap variant fetch window
        }
        reads_considered += 1;

        // ── RNA STRANDEDNESS FILTER: per-variant because gene_strand differs
        if mode == "rna" && enforce_strandedness && !rna::is_sense_strand(record, variant.gene_strand) {
            continue;
        }

        // ── MQ0 TRACKING: Count MAPQ=0 reads BEFORE any MAPQ-based skip.
        // Mirrors GATK's MappingQualityZero annotation — a high MQ0 count
        // is a locus-level red flag for regions with high homology or
        // pseudogenes, even when those reads are filtered for classification.
        if record.mapq() == 0 {
            counts.mq0_count += 1;
        }

        // ── MAPQ SKIP (Phase 1): MAPQ=0 reads were kept in the cache
        // specifically for MQ0 tracking above. Now skip them for
        // classification — they should not contribute to DP/RD/AD/DPF.
        if mode == "rna" {
            if !rna::is_valid_rna_alignment(record, min_mapq) {
                continue;
            }
        } else if record.mapq() < min_mapq {
            continue;
        }

        // ── HEURISTIC BAQ: resolve adjusted qualities for classification.
        // When BAQ is enabled, bases near indels and splice junctions
        // (CIGAR N) are downgraded. Default: off for DNA (upstream BQSR),
        // on for RNA (no upstream BQ recalibration).
        let baq_adjusted = if apply_baq {
            apply_heuristic_baq(record)
        } else {
            None
        };
        let effective_quals: &[u8] = match &baq_adjusted {
            Some(adj) => adj,
            None => record.qual(),
        };

        // ── Allele classification
        let result = check_allele_with_qual(
            record, variant, sibling_variants, effective_quals, min_baseq,
            &mut alt_aligner, &mut ref_aligner, backend,
        );
        let is_ref = result.is_ref;
        let is_alt = result.is_alt;
        let base_qual = result.qual;
        phase_counts[result.phase as usize] += 1;

        // ── DISTANCE TO READ END: Track how close the variant-supporting
        // base is to the nearest end of the read. Bases near read ends
        // have higher error rates and misalignment probability.
        // Stored per-allele for median computation after the loop.
        if is_ref || is_alt {
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
        if !overlaps_anchor {
            continue;
        }

        // ── TOTAL DEPTH: all anchor-overlapping reads count toward DP,
        // regardless of allele classification (REF, ALT, or other/ambiguous).
        // This ensures DP reflects true physical coverage at the locus.
        counts.dp += 1;
        let is_reverse = record.is_reverse();
        if is_reverse {
            counts.dp_rev += 1;
        } else {
            counts.dp_fwd += 1;
        }

        // ── FRAGMENT TRACKING: track ALL fragments for DPF.
        // FragmentEvidence::observe() correctly handles (false, false) —
        // it skips updating best_ref_qual/best_alt_qual but still tracks
        // the fragment for DPF in the downstream resolution loop.
        // UMI-aware fragment grouping: when umi_tag is set, reads with
        // different UMIs are treated as distinct molecules. The UMI is
        // extracted from the BAM aux tag (e.g., RX:Z:ACGT).
        let mol_hash = if let Some(tag) = umi_tag {
            let umi_bytes = record.aux(&tag)
                .ok()
                .and_then(|aux| match aux {
                    rust_htslib::bam::record::Aux::String(s) => Some(s.as_bytes()),
                    _ => None,
                });
            hash_molecule(record.qname(), umi_bytes)
        } else {
            hash_qname(record.qname())
        };
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
        evidence.observe(is_ref, is_alt, base_qual, is_read1, is_forward, tlen, result.has_n_base);

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
        // - Neither/REF: no any_alt/partial_alt change
        if !is_ref && !is_alt {
            // Check for partial ALT evidence before skipping
            if result.partial_match_count > 0 {
                counts.any_alt += 1;
                counts.partial_alt += 1;
                trace!("partial_alt++: {} positions matched ALT (any_alt={}, partial_alt={})",
                    result.partial_match_count, counts.any_alt, counts.partial_alt);
            }
            continue;
        }

        if is_ref {
            // Multi-allelic guard: if this read is classified as ALT for any
            // sibling variant at this locus, don't count it as REF for this
            // variant. This handles overlapping indels/complex variants where
            // a read carrying one variant's ALT could be miscounted as REF
            // for another variant at the same locus.
            if !sibling_variants.is_empty() {
                let mut is_sibling_alt = false;
                for sib in sibling_variants {
                    let sib_result = check_allele_with_qual(
                        record, sib, &[], effective_quals, min_baseq,
                        &mut alt_aligner, &mut ref_aligner, backend,
                    );
                    if sib_result.is_alt {
                        trace!(
                            "Multi-allelic guard: read is ALT for sibling {}>{} at {}:{}, \
                             excluding from REF for {}>{}",
                            sib.ref_allele, sib.alt_allele,
                            variant.chrom, variant.pos + 1,
                            variant.ref_allele, variant.alt_allele,
                        );
                        is_sibling_alt = true;
                        break;
                    }
                }
                if is_sibling_alt {
                    continue;
                }
            }
            counts.rd += 1;
            if is_reverse { counts.rd_rev += 1; } else { counts.rd_fwd += 1; }
        } else if is_alt {
            counts.ad += 1;
            counts.any_alt += 1; // Full ALT → counts toward any_alt
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
            if rna::is_sense_strand(record, variant.gene_strand) {
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

    // ── D7: RNA EDITING SITE FLAG (DB-only) ────────────────────────────
    // Flag is True ONLY when a REDIportal database is provided AND this
    // variant's position is found in the database. No DB = no flagging.
    // This avoids false positives from pattern-matching-only heuristics
    // (e.g. every A→G SNP being flagged as a potential editing site).
    if mode == "rna" {
        counts.rna_editing_site_overlap = editing_sites.as_ref().is_some_and(|sites| {
            let chrom = variant.chrom.trim_start_matches("chr");
            let found = sites.contains(&(chrom.to_string(), variant.pos));
            if found {
                trace!(
                    "D7: editing site flagged at {}:{}",
                    variant.chrom, variant.pos + 1,
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
    let mut ref_sizes:    Vec<f64> = Vec::with_capacity(fragments.len());
    let mut alt_sizes:    Vec<f64> = Vec::with_capacity(fragments.len());
    let mut nonref_sizes: Vec<f64> = Vec::with_capacity(fragments.len());
    let mut n_sizes:      Vec<f64> = Vec::with_capacity(fragments.len());

    for evidence in fragments.values() {
        let (frag_ref, frag_alt) = evidence.resolve(qual_diff_threshold);

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

        // mFSD: classify fragment into size class vectors
        if let Some(sz) = evidence.insert_size {
            if (50..=1000).contains(&sz) {
                let sz_f = sz as f64;
                if frag_ref {
                    ref_sizes.push(sz_f);
                } else if frag_alt {
                    alt_sizes.push(sz_f);
                } else if evidence.has_n_base {
                    n_sizes.push(sz_f);
                } else {
                    nonref_sizes.push(sz_f);
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

    // ── mFSD Statistics
    counts.mfsd_ref_count    = ref_sizes.len()    as u32;
    counts.mfsd_alt_count    = alt_sizes.len()    as u32;
    counts.mfsd_nonref_count = nonref_sizes.len() as u32;
    counts.mfsd_n_count      = n_sizes.len()      as u32;

    counts.mfsd_ref_mean    = mfsd::calc_mean(&ref_sizes);
    counts.mfsd_alt_mean    = mfsd::calc_mean(&alt_sizes);
    counts.mfsd_nonref_mean = mfsd::calc_mean(&nonref_sizes);
    counts.mfsd_n_mean      = mfsd::calc_mean(&n_sizes);

    counts.mfsd_alt_llr = mfsd::calc_llr(&alt_sizes);
    counts.mfsd_ref_llr = mfsd::calc_llr(&ref_sizes);

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
    (counts.mfsd_delta_alt_nonref, counts.mfsd_ks_alt_nonref, counts.mfsd_pval_alt_nonref) = ks_pair(&alt_sizes, &nonref_sizes);
    (counts.mfsd_delta_ref_nonref, counts.mfsd_ks_ref_nonref, counts.mfsd_pval_ref_nonref) = ks_pair(&ref_sizes, &nonref_sizes);
    (counts.mfsd_delta_alt_n,      counts.mfsd_ks_alt_n,      counts.mfsd_pval_alt_n)      = ks_pair(&alt_sizes, &n_sizes);
    (counts.mfsd_delta_ref_n,      counts.mfsd_ks_ref_n,      counts.mfsd_pval_ref_n)      = ks_pair(&ref_sizes, &n_sizes);
    (counts.mfsd_delta_nonref_n,   counts.mfsd_ks_nonref_n,   counts.mfsd_pval_nonref_n)   = ks_pair(&nonref_sizes, &n_sizes);

    // ── mFSD: Sub-nucleosomal / mono-nucleosomal fractions ──────────────────
    // Computed before ref_sizes/alt_sizes are consumed by into_iter().
    // Sub-nucleosomal (<150bp): ctDNA enrichment indicator.
    // Mono-nucleosomal (150–200bp): dominant cfDNA peak.
    counts.mfsd_sub_nuc_ref_frac = mfsd::calc_fraction_in_range(&ref_sizes, 0.0, 150.0);
    counts.mfsd_sub_nuc_alt_frac = mfsd::calc_fraction_in_range(&alt_sizes, 0.0, 150.0);
    counts.mfsd_sub_nuc_enrichment = if counts.mfsd_sub_nuc_ref_frac > 0.0 {
        counts.mfsd_sub_nuc_alt_frac / counts.mfsd_sub_nuc_ref_frac
    } else {
        f64::NAN
    };
    counts.mfsd_mono_nuc_ref_frac = mfsd::calc_fraction_in_range(&ref_sizes, 150.0, 200.0);
    counts.mfsd_mono_nuc_alt_frac = mfsd::calc_fraction_in_range(&alt_sizes, 150.0, 200.0);

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

    // Log per-phase classification breakdown + reads considered
    debug!(
        "Phase stats {}:{} {}→{}: P0={} P1={} P2={} P2.5={} P3={} ({} backend, {} reads from cache)",
        variant.chrom, variant.pos, variant.ref_allele, variant.alt_allele,
        phase_counts[0], phase_counts[1], phase_counts[2], phase_counts[3], phase_counts[4],
        match backend {
            AlignmentBackend::SmithWaterman => "SW",
            AlignmentBackend::PairHMM { .. } => "HMM",
        },
        reads_considered,
    );

    Ok(counts)
}



// INTENTIONAL: count_single_variant is retained for the legacy count_bam API.
// count_bam_binned uses count_bin_shared/count_variant_from_cache instead.
// Do NOT remove until count_bam itself is removed (D8b cleanup).
//
// NOTE: Consensus splicing (D6) is NOT applied in this legacy path because
// it requires buffered reads (two-pass). D6 is only available via
// count_bam_binned which has the D10 read cache.
#[allow(clippy::too_many_arguments)]
fn count_single_variant(
    bam: &mut bam::IndexedReader,
    variant: &Variant,
    sibling_variants: &[Variant],
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
) -> Result<BaseCounts> {
    let tid = bam.header().tid(variant.chrom.as_bytes()).ok_or_else(|| {
        anyhow::anyhow!("Chromosome not found in BAM: {}", variant.chrom)
    })?;

    // Fetch region around the variant. For windowed indel detection,
    // we expand the window so that reads with shifted indels are also retrieved.
    // The window scales with repeat_span to capture indels that aligners
    // shift beyond 5bp in long homopolymers/microsatellites.
    let window_pad: i64 = std::cmp::max(5, variant.repeat_span as i64 + 2);
    let start = (variant.pos - window_pad).max(0);
    let end = variant.pos + (variant.ref_allele.len() as i64) + window_pad;

    bam.fetch((tid, start, end)).context("Failed to fetch region")?;

    let mut counts = BaseCounts::default();

    // Fragment tracking: QNAME hash -> FragmentEvidence
    // Using u64 hash keys instead of String for memory efficiency.
    let mut fragments: HashMap<u64, FragmentEvidence> = HashMap::new();

    // Quality threshold for fragment consensus tiebreaking.
    // When R1 and R2 disagree, the allele with higher quality wins
    // only if the quality difference exceeds this threshold.
    // Configurable via --fragment-qual-threshold (default: 10).
    let qual_diff_threshold: u8 = fragment_qual_threshold;

    // Distance-to-read-end tracking for QC metrics.
    // Pre-allocated with conservative capacity to avoid resizing.
    let mut alt_dists: Vec<u32> = Vec::with_capacity(500);
    let mut ref_dists: Vec<u32> = Vec::with_capacity(500);

    // Create SW aligners ONCE per variant, not per read (indelpost pattern).
    // bio::alignment::pairwise::Aligner reuses internal DP buffers on
    // subsequent calls, avoiding repeated O(n×m) heap allocation.
    // Score function: N in read or haplotype → 0 (neutral, uninformative).
    // Matches GATK approach: N contributes no evidence for match or mismatch.
    // Handles N-in-read (duplex masking / sequencer failure) and N-in-haplotype (rare).
    let score_fn = |a: u8, b: u8| -> i32 {
        if a == b'N' || b == b'N' { 0 } else if a == b { 1 } else { -1 }
    };
    // ALT + REF: Same affine gap penalties for fair comparison.
    // Continuous gap_extend: uses logistic sigmoid to smoothly transition
    // from tight (-1) to free (0) as repeat_span increases.
    // Replaces the previous rigid `repeat_span >= 10` binary threshold
    // to prevent boundary artifacts at the transition point.
    let gap_open: i32 = -5;
    let gap_extend: i32 = dynamic_sw_gap_extend(variant.repeat_span);
    let mut alt_aligner = Aligner::new(gap_open, gap_extend, &score_fn);
    let mut ref_aligner = Aligner::new(gap_open, gap_extend, &score_fn);

    // Per-phase classification counters
    // Indices: 0=Structural, 1=CigarRecon, 2=MaskedCompare, 3=Levenshtein, 4=Alignment
    let mut phase_counts = [0u32; 5];

    // Construct ReadFilter from the boolean params.
    let read_filter = crate::shared::filters::ReadFilter {
        filter_duplicates,
        filter_secondary,
        filter_supplementary,
        filter_qc_failed,
        filter_improper_pair,
        filter_indel,
    };
    let mut filter_counts = crate::shared::filters::FilterCounts::default();

    for result in bam.records() {
        let record = result.context("Error reading BAM record")?;

        // Universal flag filters (delegated to shared::filters::ReadFilter).
        if !read_filter.passes(&record, &mut filter_counts) {
            continue;
        }

        // ── MQ0 TRACKING: Count MAPQ=0 reads BEFORE the MAPQ filter.
        // Mirrors GATK's MappingQualityZero annotation — a high MQ0
        // count is a locus-level red flag for regions with high homology
        // or pseudogenes, even when those reads are filtered out.
        if record.mapq() == 0 {
            counts.mq0_count += 1;
        }

        // ── RNA MAPQ FILTER: In RNA mode, use NH:i:1 rescue logic.
        // STAR assigns low MAPQ to reads at novel splice junctions despite
        // unique mapping — NH:i:1 rescue recovers these informative reads.
        // In DNA mode, use standard MAPQ filter.
        if mode == "rna" {
            if !rna::is_valid_rna_alignment(&record, min_mapq) {
                continue;
            }
        } else if record.mapq() < min_mapq {
            continue;
        }

        // ── RNA STRANDEDNESS FILTER: In RNA mode with strandedness enforced,
        // reject reads on the wrong strand relative to the gene annotation.
        // This prevents antisense artifacts from inflating variant counts.
        if mode == "rna" && enforce_strandedness && !rna::is_sense_strand(&record, variant.gene_strand) {
            continue;
        }

        // ── HEURISTIC BAQ: When enabled, downgrade base qualities near
        // alignment indels and splice junctions before allele classification.
        // This affects which bases pass the min_baseq gate in
        // check_snp/check_mnp and the quality used for fragment consensus
        // tiebreaking.
        //
        // BAQ is applied lazily: apply_heuristic_baq() returns None for
        // reads without indels or splice junctions (zero allocation for
        // the common case).
        let baq_adjusted = if apply_baq {
            apply_heuristic_baq(&record)
        } else {
            None
        };
        let effective_quals: &[u8] = match &baq_adjusted {
            Some(adj) => adj,
            None => record.qual(),
        };

        // Determine allele status and base quality at variant position
        let result = check_allele_with_qual(
            &record, variant, sibling_variants, effective_quals, min_baseq, &mut alt_aligner, &mut ref_aligner, backend,
        );
        let is_ref = result.is_ref;
        let is_alt = result.is_alt;
        let base_qual = result.qual;
        phase_counts[result.phase as usize] += 1;

        // ── DISTANCE TO READ END: Track how close the variant-supporting
        // base is to the nearest end of the read. Bases near read ends
        // have higher error rates and misalignment probability.
        // Stored per-allele for median computation after the loop.
        if is_ref || is_alt {
            if let Some(read_idx) = find_read_pos(&record, variant.pos) {
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
        let read_start = record.pos();
        let read_end = read_ref_end(&record);
        let overlaps_anchor = read_start <= variant.pos && read_end > variant.pos;
        if !overlaps_anchor {
            continue;
        }

        // ── TOTAL DEPTH: all anchor-overlapping reads count toward DP,
        // regardless of allele classification (REF, ALT, or other/ambiguous).
        // This ensures DP reflects true physical coverage at the locus.
        counts.dp += 1;
        let is_reverse = record.is_reverse();

        if is_reverse {
            counts.dp_rev += 1;
        } else {
            counts.dp_fwd += 1;
        }

        // ── FRAGMENT TRACKING: track ALL fragments for DPF.
        // FragmentEvidence::observe() correctly handles (false, false) —
        // it skips updating best_ref_qual/best_alt_qual but still tracks
        // the fragment for DPF in the downstream resolution loop.
        // UMI-aware fragment grouping: when umi_tag is set, reads with
        // different UMIs are treated as distinct molecules. The UMI is
        // extracted from the BAM aux tag (e.g., RX:Z:ACGT).
        let mol_hash = if let Some(tag) = umi_tag {
            let umi_bytes = record.aux(&tag)
                .ok()
                .and_then(|aux| match aux {
                    rust_htslib::bam::record::Aux::String(s) => Some(s.as_bytes()),
                    _ => None,
                });
            hash_molecule(record.qname(), umi_bytes)
        } else {
            hash_qname(record.qname())
        };
        let is_read1 = record.is_first_in_template();
        let is_forward = !is_reverse;

        let evidence = fragments.entry(mol_hash).or_insert_with(FragmentEvidence::new);

        // mFSD: compute physical fragment size from CIGAR, correcting TLEN for indels.
        // Formula: physical = |TLEN| - D + I (validated on real MSK-ACCESS BAMs).
        // observe() stores min(R1, R2) for defensive correctness.
        // is_n_base: a fragment is in the N class when the base at the variant
        // position was 'N' — proxied by base_qual==0 with no REF or ALT call.
        let tlen = mfsd::calc_physical_insert_size(&record);
        let is_n_base = base_qual == 0 && !is_ref && !is_alt;

        evidence.observe(is_ref, is_alt, base_qual, is_read1, is_forward, tlen, is_n_base);

        // ── ALLELE-SPECIFIC COUNTS: only REF/ALT reads contribute to RD/AD.
        // DP and DPF are already recorded above.
        //
        // Decomposed counting (any_alt / partial_alt):
        // - Full ALT match: ad++, any_alt++ (invariant: any_alt = ad + partial_alt)
        // - Partial ALT match (some discriminating positions match ALT): any_alt++, partial_alt++
        // - Neither/REF: no any_alt/partial_alt change
        if !is_ref && !is_alt {
            // Check for partial ALT evidence before skipping
            if result.partial_match_count > 0 {
                counts.any_alt += 1;
                counts.partial_alt += 1;
            }
            continue;
        }

        if is_ref {
            // Multi-allelic guard: if this read is classified as ALT for any
            // sibling variant at this locus, don't count it as REF for this
            // variant. This handles overlapping indels/complex variants where
            // a read carrying one variant's ALT could be miscounted as REF
            // for another variant at the same locus.
            if !sibling_variants.is_empty() {
                let mut is_sibling_alt = false;
                for sib in sibling_variants {
                    let sib_result = check_allele_with_qual(
                        &record, sib, &[], effective_quals, min_baseq,
                        &mut alt_aligner, &mut ref_aligner, backend,
                    );
                    let sib_alt = sib_result.is_alt;
                    if sib_alt {
                        trace!(
                            "Multi-allelic guard: read is ALT for sibling {}>{} at {}:{}, \
                             excluding from REF for {}>{}",
                            sib.ref_allele, sib.alt_allele,
                            variant.chrom, variant.pos + 1,
                            variant.ref_allele, variant.alt_allele,
                        );
                        is_sibling_alt = true;
                        break;
                    }
                }
                if is_sibling_alt {
                    continue; // Skip REF counting — this read belongs to a sibling
                }
            }
            counts.rd += 1;
            if is_reverse {
                counts.rd_rev += 1;
            } else {
                counts.rd_fwd += 1;
            }
        } else if is_alt {
            counts.ad += 1;
            counts.any_alt += 1; // Full ALT → counts toward any_alt
            if is_reverse {
                counts.ad_rev += 1;
            } else {
                counts.ad_fwd += 1;
            }

            // ── RNA-SPECIFIC ALT TRACKING ──
            if mode == "rna" {
                // Splice-spanning count: ALT reads that cross a splice junction
                if rna::has_splice_junction(&record) {
                    counts.splice_spanning_count += 1;
                }
            }
        }

        // ── RNA SENSE/ANTISENSE DEPTH: track strand-specific depth.
        // Uses the same dUTP logic as is_sense_strand to classify reads.
        if mode == "rna" {
            if rna::is_sense_strand(&record, variant.gene_strand) {
                counts.sense_depth += 1;
                if is_alt {
                    counts.sense_strand_alt_count += 1;
                }
            } else {
                counts.antisense_depth += 1;
                if is_alt {
                    counts.antisense_strand_alt_count += 1;
                }
            }
        }
    }

    // ── QC MEDIAN COMPUTATION: compute median distance-to-end for REF/ALT.
    counts.alt_dist_end_median = compute_median_u32(&mut alt_dists);
    counts.ref_dist_end_median = compute_median_u32(&mut ref_dists);

    // ── RNA EDITING SITE FLAG (legacy path) ──
    // D7 DB-only editing flagging is only available via count_bam_binned.
    // Legacy path does not receive the editing_sites HashSet, so
    // rna_editing_site_overlap stays false (no pattern-matching guessing).

    // Resolve fragment-level counts using quality-weighted consensus.
    // Each fragment contributes exactly ONE allele call (REF xor ALT),
    // preventing the double-counting bug where R1=REF + R2=ALT
    // inflated both rdf and adf.
    //
    // Strand bias uses allele-specific orientation: the strand of the read
    // that provided the best evidence for the winning allele, not just R1.
    // Example: R1=Fwd/REF(Q10) + R2=Rev/ALT(Q30) → ALT wins, counted as
    // adf_rev (not adf_fwd).
    // ── mFSD size vectors: one per Krewlyzer fragment class ─────────────────
    // Populated below during resolution. Only sizes in the cfDNA-valid range
    // (50–1000 bp) with a known TLEN are added. GC correction is not applied —
    // GC bias affects count depth, not fragment length, so these raw sizes are
    // already unbiased samples of the true size distribution.
    let mut ref_sizes:    Vec<f64> = Vec::with_capacity(fragments.len());
    let mut alt_sizes:    Vec<f64> = Vec::with_capacity(fragments.len());
    let mut nonref_sizes: Vec<f64> = Vec::with_capacity(fragments.len());
    let mut n_sizes:      Vec<f64> = Vec::with_capacity(fragments.len());

    for evidence in fragments.values() {
        let (frag_ref, frag_alt) = evidence.resolve(qual_diff_threshold);

        // Count every fragment in dpf regardless of consensus outcome.
        // Discarded fragments (ambiguous R1-vs-R2 within quality threshold)
        // are still real molecules — tracking them in dpf makes the gap
        // dpf - (rdf + adf) a useful quality metric for the locus.
        counts.dpf += 1;

        if frag_ref {
            counts.rdf += 1;
            // Use REF-specific orientation (strand of best REF evidence)
            if let Some(ori) = evidence.ref_orientation() {
                if ori {
                    counts.rdf_fwd += 1;
                } else {
                    counts.rdf_rev += 1;
                }
            }
        } else if frag_alt {
            counts.adf += 1;
            // Use ALT-specific orientation (strand of best ALT evidence)
            if let Some(ori) = evidence.alt_orientation() {
                if ori {
                    counts.adf_fwd += 1;
                } else {
                    counts.adf_rev += 1;
                }
            }
        }

        // mFSD: classify fragment into one of four size class vectors.
        // Only fragments with a known, in-range insert size contribute.
        if let Some(sz) = evidence.insert_size {
            if (50..=1000).contains(&sz) {
                let sz_f = sz as f64;
                if frag_ref {
                    ref_sizes.push(sz_f);
                } else if frag_alt {
                    alt_sizes.push(sz_f);
                } else if evidence.has_n_base {
                    // N class: ambiguous base at variant position
                    n_sizes.push(sz_f);
                } else {
                    // NonREF class: definite non-ref, non-alt, non-N base
                    nonref_sizes.push(sz_f);
                }
            }
        }
    }

    // Calculate stats
    let (sb_pval, sb_or) =
        fisher_strand_bias(counts.rd_fwd, counts.rd_rev, counts.ad_fwd, counts.ad_rev);
    counts.sb_pval = sb_pval;
    counts.sb_or = sb_or;

    let (fsb_pval, fsb_or) = fisher_strand_bias(
        counts.rdf_fwd,
        counts.rdf_rev,
        counts.adf_fwd,
        counts.adf_rev,
    );
    counts.fsb_pval = fsb_pval;
    counts.fsb_or = fsb_or;

    // ── mFSD Statistics ──────────────────────────────────────────────────────
    // All distributional tests (KS, LLR, delta) use raw unweighted size arrays.
    // Counts populated here match mfsd_*_count fields on BaseCounts.

    counts.mfsd_ref_count    = ref_sizes.len()    as u32;
    counts.mfsd_alt_count    = alt_sizes.len()    as u32;
    counts.mfsd_nonref_count = nonref_sizes.len() as u32;
    counts.mfsd_n_count      = n_sizes.len()      as u32;

    // Means (0.0 for empty — callers should gate on mfsd_*_count)
    counts.mfsd_ref_mean    = mfsd::calc_mean(&ref_sizes);
    counts.mfsd_alt_mean    = mfsd::calc_mean(&alt_sizes);
    counts.mfsd_nonref_mean = mfsd::calc_mean(&nonref_sizes);
    counts.mfsd_n_mean      = mfsd::calc_mean(&n_sizes);

    // LLR: sum(log P_tumor / P_healthy) per fragment; positive = tumor-like
    counts.mfsd_alt_llr = mfsd::calc_llr(&alt_sizes);
    counts.mfsd_ref_llr = mfsd::calc_llr(&ref_sizes);

    // KS helper closure: pairwise delta + D-statistic + p-value
    // delta = mean(a) - mean(b);  (NaN, 1.0) when either class < MIN_FOR_KS
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
    (counts.mfsd_delta_alt_nonref, counts.mfsd_ks_alt_nonref, counts.mfsd_pval_alt_nonref) = ks_pair(&alt_sizes, &nonref_sizes);
    (counts.mfsd_delta_ref_nonref, counts.mfsd_ks_ref_nonref, counts.mfsd_pval_ref_nonref) = ks_pair(&ref_sizes, &nonref_sizes);
    (counts.mfsd_delta_alt_n,      counts.mfsd_ks_alt_n,      counts.mfsd_pval_alt_n)      = ks_pair(&alt_sizes, &n_sizes);
    (counts.mfsd_delta_ref_n,      counts.mfsd_ks_ref_n,      counts.mfsd_pval_ref_n)      = ks_pair(&ref_sizes, &n_sizes);
    (counts.mfsd_delta_nonref_n,   counts.mfsd_ks_nonref_n,   counts.mfsd_pval_nonref_n)   = ks_pair(&nonref_sizes, &n_sizes);

    // ── mFSD: Sub-nucleosomal / mono-nucleosomal fractions ──────────────────
    // Computed before ref_sizes/alt_sizes are consumed by into_iter().
    // Sub-nucleosomal (<150bp): ctDNA enrichment indicator.
    // Mono-nucleosomal (150–200bp): dominant cfDNA peak.
    counts.mfsd_sub_nuc_ref_frac = mfsd::calc_fraction_in_range(&ref_sizes, 0.0, 150.0);
    counts.mfsd_sub_nuc_alt_frac = mfsd::calc_fraction_in_range(&alt_sizes, 0.0, 150.0);
    counts.mfsd_sub_nuc_enrichment = if counts.mfsd_sub_nuc_ref_frac > 0.0 {
        counts.mfsd_sub_nuc_alt_frac / counts.mfsd_sub_nuc_ref_frac
    } else {
        f64::NAN
    };
    counts.mfsd_mono_nuc_ref_frac = mfsd::calc_fraction_in_range(&ref_sizes, 150.0, 200.0);
    counts.mfsd_mono_nuc_alt_frac = mfsd::calc_fraction_in_range(&alt_sizes, 150.0, 200.0);

    // Store raw size arrays for --mfsd-parquet export
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

    // Log per-phase classification breakdown
    debug!(
        "Phase stats {}:{} {}→{}: P0(Structural)={} P1(CigarRecon)={} P2(Masked)={} P2.5(Lev)={} P3(Align)={} ({} backend)",
        variant.chrom, variant.pos, variant.ref_allele, variant.alt_allele,
        phase_counts[0], phase_counts[1], phase_counts[2], phase_counts[3], phase_counts[4],
        match backend {
            AlignmentBackend::SmithWaterman => "SW",
            AlignmentBackend::PairHMM { .. } => "HMM",
        }
    );

    Ok(counts)
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
/// once per variant in `count_single_variant()` and threaded through to
/// avoid per-read allocation (indelpost pattern).
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
) -> ClassifyResult {
    // Dispatch based on allele lengths rather than the variant_type string.
    // This is more robust than relying on upstream type labels, which can be
    // inconsistent (e.g., a caller emitting "COMPLEX" for what is really a
    // pure deletion after normalization).
    let ref_len = variant.ref_allele.len();
    let alt_len = variant.alt_allele.len();
    trace!(
        "check_allele ref_len={} alt_len={} pos={} ref={} alt={}",
        ref_len, alt_len, variant.pos, variant.ref_allele, variant.alt_allele
    );

    if ref_len == 1 && alt_len == 1 {
        // SNP: single base substitution — no Phase 3 needed
        check_snp(record, variant, quals, min_baseq)
    } else if ref_len == alt_len {
        // MNP: selective discriminating-position quality gate with no Phase 3 fallback.
        match check_mnp(record, variant, quals, min_baseq) {
            MnpResult::Ref(q, had_n) => {
                let mut r = ClassifyResult::is_ref(q, ClassifyPhase::MaskedCompare);
                r.has_n_base = had_n;
                r
            }
            MnpResult::Alt(q, had_n) => {
                let mut r = ClassifyResult::is_alt(q, ClassifyPhase::MaskedCompare);
                r.has_n_base = had_n;
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
                trace!("MNP structural issue, falling back to Phase 3");
                check_complex(record, variant, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend)
            }
        }
    } else if ref_len == 1 {
        // Pure insertion: CIGAR-based fast paths, then backend-aware Phase 3 fallback
        check_insertion(record, variant, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend)
    } else if alt_len == 1 {
        // Distinguish pure deletion (anchor base preserved) from complex Del+SNV
        // (anchor base also substituted, e.g. GC→T where G is both the anchor
        // AND changes to T).
        //
        // Pure deletion example:  GC → G  (anchor G kept, C deleted)
        // Complex Del+SNV:        GC → T  (C deleted AND G→T at anchor position)
        //
        // check_deletion's S3 safeguard validates shifted deletions by comparing
        // the deleted reference bases against `expected_del_seq` (ref_allele[1..]).
        // For complex Del+SNV the anchor mismatch causes a cascade: reads whose
        // deletion left-shifts away from the anchor pass the S3 check at the
        // shifted position only if the reference base there matches expected_del_seq;
        // when it doesn't S3 rejects the windowed match, `found_ref_coverage` is
        // set to true (the anchor M-block still covers the anchor), and the read
        // is definitively classified as REF — hiding the true ALT reads.
        //
        // Routing complex Del+SNV to check_complex lets Phase 3 (PairHMM/SW)
        // align the read against the full REF and ALT haplotype contexts where
        // both the deletion and the anchor substitution are captured correctly.
        let anchor_preserved = variant
            .alt_allele
            .as_bytes()
            .first()
            .map(|b| b.to_ascii_uppercase())
            == variant
                .ref_allele
                .as_bytes()
                .first()
                .map(|b| b.to_ascii_uppercase());

        if anchor_preserved {
            // Pure deletion: CIGAR-based fast paths, then backend-aware Phase 3 fallback
            check_deletion(record, variant, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend)
        } else {
            // Complex Del+SNV: anchor base also substituted — route to Phase 3
            trace!(
                "Complex Del+SNV at {}:{} (ref[0]={} ≠ alt[0]={}): routing to check_complex",
                variant.chrom,
                variant.pos + 1,
                variant.ref_allele.chars().next().unwrap_or('?'),
                variant.alt_allele.chars().next().unwrap_or('?'),
            );
            check_complex(record, variant, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend)
        }
    } else {
        // Complex: ref_len != alt_len, both > 1 (e.g., DelIns)
        check_complex(record, variant, siblings, quals, min_baseq, alt_aligner, ref_aligner, backend)
    }
}




#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;
    use std::ffi::CString;

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
        }
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
            MnpResult::Alt(q, _) => assert!(q >= 20, "Quality should be >= 20, got {}", q),
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
            MnpResult::Alt(q, _) => assert!(q > 0, "Expected Alt with quality > 0, got {}", q),
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
            MnpResult::Alt(q, _) => assert!(q >= 20, "Quality should be >= 20, got {}", q),
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
            MnpResult::Alt(q, _) => assert!(q >= 20, "Should pass at exact threshold boundary"),
            other => panic!("Expected Alt at exact BQ threshold, got {:?}",
                           format_mnp_result(&other)),
        }
    }

    /// Helper to format MnpResult for panic messages
    fn format_mnp_result(result: &MnpResult) -> String {
        match result {
            MnpResult::Ref(q, n) => format!("Ref(q={}, had_n={})", q, n),
            MnpResult::Alt(q, n) => format!("Alt(q={}, had_n={})", q, n),
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
            MnpResult::Alt(q, _) => assert!(q >= 20,
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
            MnpResult::Alt(q, _) => assert!(q > 0, "Expected Alt with quality > 0, got {}", q),
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
        // After G1 fix: positions_matching_alt only counts UNMASKED positions,
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
            MnpResult::Alt(q, had_n) => {
                assert!(q > 0,
                    "Expected Alt when N masks one position but other unmasked matches ALT, got q={}", q);
                assert!(had_n,
                    "Expected had_n=true when N base present at discriminating position");
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
        }
    }

    /// Run check_allele_with_qual through both SW and HMM backends for concordance testing.
    fn run_both_backends(
        record: &Record, variant: &Variant, min_baseq: u8,
    ) -> (ClassifyResult, ClassifyResult) {
        // We need separate aligner instances because check_allele_with_qual takes &mut
        let scoring_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

        let mut alt_a1 = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(-5, -1, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE)
                .yclip(0),
        );
        let mut ref_a1 = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(-5, -1, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE)
                .yclip(0),
        );

        let quals = record.qual();
        let sw_result = check_allele_with_qual(
            record, variant, &[], quals, min_baseq,
            &mut alt_a1, &mut ref_a1,
            &AlignmentBackend::SmithWaterman,
        );
        let hmm_result = check_allele_with_qual(
            record, variant, &[], quals, min_baseq,
            &mut alt_a1, &mut ref_a1,
            &AlignmentBackend::pairhmm_default(),
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

    // ── G5: n_count accumulation integration tests ──
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
            200, 200, bio::alignment::pairwise::Scoring::new(-5, -1, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let mut ref_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(-5, -1, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );

        let result = check_allele_with_qual(
            &record, &variant, &[], record.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman,
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
            MnpResult::Alt(_, had_n) => {
                assert!(had_n, "MNP Alt with N at one position should report had_n=true");
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
            200, 200, bio::alignment::pairwise::Scoring::new(-5, -1, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let mut ref_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(-5, -1, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );

        let result = check_allele_with_qual(
            &record, &variant, &[], record.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman,
        );

        assert!(result.is_alt, "MNP with N-masked + ALT-matching should be ALT");
        assert!(result.has_n_base, "MNP ALT with N at one position should propagate has_n_base=true");
    }

    // ── G4: check_complex N base propagation tests ──
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
            200, 200, bio::alignment::pairwise::Scoring::new(-5, -1, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let mut ref_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(-5, -1, &scoring_fn)
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
            200, 200, bio::alignment::pairwise::Scoring::new(-5, -1, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );
        let mut ref_a = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(-5, -1, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE).yclip(0),
        );

        let result = check_complex(
            &record, &variant, &[], record.qual(), 20,
            &mut alt_a, &mut ref_a, &AlignmentBackend::SmithWaterman,
        );

        assert!(!result.has_n_base,
            "check_complex with no N in reconstructed haplotype should have has_n_base=false");
    }
}
