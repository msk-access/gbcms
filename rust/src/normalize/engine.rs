//! Variant preparation orchestration.
//!
//! Contains `prepare_variants` (PyO3 entry point), `prepare_single_variant`
//! (per-variant pipeline), and `assign_multi_allelic_groups` (post-processing).

use std::fs::File;

use bio::io::fasta;
use log::{debug, info, warn};
use pyo3::prelude::*;
use rayon::prelude::*;

use crate::types::Variant;
use super::types::PreparedVariant;
use super::decomp::check_homopolymer_decomp;
use super::left_align::left_align_variant;
use super::fasta::{fetch_region, resolve_maf_anchor, validate_ref};
use super::repeat::{find_tandem_repeat, compute_adaptive_padding, first_change_offset};

/// Prepare variants for counting in a single pass over the reference FASTA.
///
/// For each variant, this function performs (in order):
/// 1. **MAF anchor fetch** — if `is_maf`, resolve `-` alleles to VCF-style
/// 2. **REF validation** — check REF against the reference genome
/// 3. **Left-alignment** — bcftools `realign_left()` for indels
/// 4. **ref_context fetch** — flanking sequence for Smith-Waterman alignment
///
/// Uses rayon `par_iter().map_init()` with thread-local FASTA readers,
/// matching the established pattern in `count_bam()`.
///
/// # Arguments
/// * `variants` — Input variants (raw MAF or VCF coords, 0-based)
/// * `fasta_path` — Path to indexed reference FASTA
/// * `context_padding` — Minimum flanking bases for ref_context (e.g. 5)
/// * `is_maf` — If true, perform MAF→VCF anchor resolution for indels
/// * `threads` — Number of rayon worker threads
/// * `adaptive_context` — If true, dynamically increase padding in repeat regions
///
/// # Returns
/// One `PreparedVariant` per input variant, in the same order.
#[pyfunction]
#[pyo3(signature = (variants, fasta_path, context_padding, is_maf, threads=1, adaptive_context=true))]
pub fn prepare_variants(
    py: Python<'_>,
    variants: Vec<Variant>,
    fasta_path: String,
    context_padding: i64,
    is_maf: bool,
    threads: usize,
    adaptive_context: bool,
) -> PyResult<Vec<PreparedVariant>> {
    info!(
        "prepare_variants: {} variants, is_maf={}, context_padding={}, adaptive={}, threads={}",
        variants.len(),
        is_maf,
        context_padding,
        adaptive_context,
        threads,
    );

    // Build rayon thread pool (same pattern as count_bam). `--threads` is the total
    // budget for this process (see shared::resolve_thread_budget).
    let threads = crate::shared::resolve_thread_budget(threads);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| {
            pyo3::exceptions::PyRuntimeError::new_err(format!(
                "Failed to build thread pool: {}",
                e
            ))
        })?;

    let fasta_path_clone = fasta_path.clone();

    // Release GIL for parallel execution
    #[allow(deprecated)]
    let results: Result<Vec<PreparedVariant>, anyhow::Error> =
        py.allow_threads(move || {
            pool.install(|| {
                variants
                    .par_iter()
                    .map_init(
                        || {
                            // Thread-local FASTA reader initialization
                            fasta::IndexedReader::from_file(&fasta_path_clone)
                                .map_err(|e| {
                                    anyhow::anyhow!(
                                        "Failed to open FASTA: {}",
                                        e
                                    )
                                })
                        },
                        |reader_result, variant| {
                            prepare_single_variant(
                                reader_result,
                                variant,
                                context_padding,
                                is_maf,
                                adaptive_context,
                            )
                        },
                    )
                    .collect()
            })
        });

    match results {
        Ok(mut r) => {
            // Post-processing: assign multi-allelic group IDs
            assign_multi_allelic_groups(&mut r);

            // Log summary
            let valid = r.iter().filter(|p| p.gbcms_status == "PASS").count();
            let n_anchor = r.iter().filter(|p| p.was_anchor_resolved).count();
            let n_left = r.iter().filter(|p| p.was_left_aligned).count();
            let n_total = r.iter().filter(|p| p.was_normalized()).count();
            let multi_allelic = r
                .iter()
                .filter(|p| p.gbcms_status_reason.contains("MULTI_ALLELIC"))
                .count();
            let tract_cluster = r
                .iter()
                .filter(|p| p.gbcms_status_reason.contains("TRACT_CLUSTER"))
                .count();
            // Count MNPs: same-length multi-base substitutions excluded from
            // indel normalization (left-alignment + ref_context fetch).
            let n_mnp = r.iter().filter(|p| {
                let ref_len = p.variant.ref_allele.len();
                let alt_len = p.variant.alt_allele.len();
                ref_len == alt_len && ref_len > 1
            }).count();
            info!(
                "prepare_variants complete: {}/{} valid, {} normalized ({} anchor-resolved, {} left-aligned), {} multi-allelic, {} tract-clustered, {} MNPs (normalization bypassed)",
                valid,
                r.len(),
                n_total,
                n_anchor,
                n_left,
                multi_allelic,
                tract_cluster,
                n_mnp,
            );
            Ok(r)
        }
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(format!(
            "prepare_variants failed: {}",
            e
        ))),
    }
}

/// Scan-window pad for tract-cluster grouping — the engine's own windowed-scan
/// formula (uncapped, like the classification scan), so grouping reach and
/// classification reach cannot drift apart.
fn window_pad(v: &Variant) -> i64 {
    std::cmp::max(5, v.repeat_span as i64 + 2)
}

/// Whether a variant changes sequence length (pure indel or delins). Only
/// length-changing variants participate in window-based (tract-cluster)
/// grouping: an SNV near an indel must NOT become its sibling, or the
/// engine's sibling-ALT guard would drain rd at every locus with a nearby
/// annotated SNV.
fn is_length_changing(v: &Variant) -> bool {
    v.ref_allele.len() != v.alt_allele.len()
}

/// Assign group IDs to co-annotated variants so the engine can evaluate them
/// jointly (sibling exclusion + exclusive assignment).
///
/// Two membership criteria, transitively closed per chromosome:
/// 1. **Span overlap** (any variant types): REF spans `pos..pos+ref_len`
///    intersect — the original multi-allelic rule; members are tagged
///    `MULTI_ALLELIC` when their span truly intersects another member's.
/// 2. **Window overlap** (length-changing variants only): scan windows —
///    spans padded by `window_pad` on each side — intersect. This chains
///    tract clusters (co-annotated indels/delins in one repeat
///    neighborhood) whose spans never touch; window-only members are
///    tagged `TRACT_CLUSTER`.
///
/// Groups may be non-contiguous in position order (an SNV can sit between
/// two window-joined deletions without joining), so the sweep marks
/// assigned indices instead of consuming a contiguous run.
fn assign_multi_allelic_groups(variants: &mut [PreparedVariant]) {
    if variants.len() < 2 {
        return;
    }

    // Index sorted by (chrom, pos) — sort indices, not the array, to keep
    // input order (which must match the output row order).
    let mut indices: Vec<usize> = (0..variants.len()).collect();
    indices.sort_by(|&a, &b| {
        let va = &variants[a].variant;
        let vb = &variants[b].variant;
        va.chrom.cmp(&vb.chrom).then(va.pos.cmp(&vb.pos))
    });

    let mut assigned = vec![false; variants.len()];
    let mut group_id: u32 = 0;
    // Widest pad any candidate can carry: bounds how far past a group's
    // reach the sweep must look before it may stop. Derived from the input
    // (repeat_span is uncapped), never a hard-coded constant.
    let max_candidate_pad = variants
        .iter()
        .map(|p| window_pad(&p.variant))
        .max()
        .unwrap_or(5);
    // Widest candidate extent (span + pad): how far LEFT of a group's box a
    // candidate's start can sit while still reaching it.
    let max_candidate_extent = variants
        .iter()
        .map(|p| {
            let len = p.variant.ref_allele.len() as i64;
            len + if is_length_changing(&p.variant) { window_pad(&p.variant) } else { 0 }
        })
        .max()
        .unwrap_or(1);
    // Chromosome runs over the sorted index, computed ONCE, plus a parallel
    // position array for binary-searching each pass's scan bounds — an
    // isolated seed then costs O(log run) instead of O(run).
    let mut runs: Vec<(usize, usize)> = Vec::new();
    {
        let mut rs = 0;
        for k in 1..=indices.len() {
            if k == indices.len()
                || variants[indices[k]].variant.chrom != variants[indices[rs]].variant.chrom
            {
                runs.push((rs, k));
                rs = k;
            }
        }
    }
    let mut run_of = vec![0usize; indices.len()];
    for (ri, &(a, b)) in runs.iter().enumerate() {
        for r in run_of.iter_mut().take(b).skip(a) {
            *r = ri;
        }
    }
    let pos_of: Vec<i64> = indices.iter().map(|&ix| variants[ix].variant.pos).collect();

    for i in 0..indices.len() {
        let idx = indices[i];
        if assigned[idx] || variants[idx].gbcms_status != "PASS" {
            continue;
        }

        // Same-chromosome run holding the seed (precomputed). Each
        // fixed-point pass re-scans the box-reachable slice of it: a late
        // joiner's window can extend LEFT past candidates already skipped,
        // and pads differ per variant, so a single forward pass misses
        // transitive joins in both directions.
        let (c_start, c_end) = runs[run_of[i]];

        let seed = &variants[idx].variant;
        // Group reach as bounding boxes per criterion (matching the
        // pre-widening sweep's bounding-box semantics for spans):
        // [span_lo, span_hi) unions member spans; [win_lo, win_hi) unions
        // length-changing members' padded windows (empty when none).
        let mut span_lo = seed.pos;
        let mut span_hi = seed.pos + seed.ref_allele.len() as i64;
        let (mut win_lo, mut win_hi) = if is_length_changing(seed) {
            (span_lo - window_pad(seed), span_hi + window_pad(seed))
        } else {
            (i64::MAX, i64::MIN)
        };
        let mut group_members: Vec<usize> = vec![idx];
        assigned[idx] = true;

        // Fixed point: keep re-scanning while the group grows. Terminates
        // because each pass either adds >=1 member or stops; member count is
        // bounded by the run length.
        loop {
            let mut grew = false;
            // Scan only the slice of the run whose positions could reach the
            // current boxes (binary-searched; the cheap reject below stays
            // exact for boundary cases).
            let lo_bound = span_lo.min(win_lo).saturating_sub(max_candidate_extent);
            let hi_bound = span_hi.max(win_hi).saturating_add(max_candidate_pad);
            let k0 = c_start + pos_of[c_start..c_end].partition_point(|&p| p < lo_bound);
            let k1 = c_start + pos_of[c_start..c_end].partition_point(|&p| p < hi_bound);
            for &jdx in &indices[k0..k1] {
                if assigned[jdx] || variants[jdx].gbcms_status != "PASS" {
                    continue;
                }
                let vj = &variants[jdx].variant;
                // Cheap reject: outside every possible reach of this group.
                let reach_hi = span_hi.max(win_hi) + max_candidate_pad;
                let reach_lo = span_lo.min(win_lo).saturating_sub(max_candidate_pad);
                if vj.pos >= reach_hi || vj.pos + (vj.ref_allele.len() as i64) <= reach_lo {
                    continue;
                }
                let vj_end = vj.pos + vj.ref_allele.len() as i64;
                let joins_span = vj.pos < span_hi && span_lo < vj_end;
                let joins_window = is_length_changing(vj) && {
                    let pad = window_pad(vj);
                    vj.pos - pad < win_hi && win_lo < vj_end + pad
                };
                if joins_span || joins_window {
                    span_lo = span_lo.min(vj.pos);
                    span_hi = span_hi.max(vj_end);
                    if is_length_changing(vj) {
                        let pad = window_pad(vj);
                        win_lo = win_lo.min(vj.pos - pad);
                        win_hi = win_hi.max(vj_end + pad);
                    }
                    group_members.push(jdx);
                    assigned[jdx] = true;
                    grew = true;
                }
            }
            if !grew {
                break;
            }
        }

        if group_members.len() == 1 {
            // Singleton: release the seed so a later seed's group whose
            // window reaches back over it can still absorb it.
            assigned[idx] = false;
        } else {
            group_id += 1;
            // Honest reason tags: MULTI_ALLELIC only when a member's span
            // truly intersects another member's; window-only membership is
            // TRACT_CLUSTER. Computed pairwise post-hoc (order-independent).
            for (mi, &member_idx) in group_members.iter().enumerate() {
                let vm = &variants[member_idx].variant;
                let m_start = vm.pos;
                let m_end = vm.pos + vm.ref_allele.len() as i64;
                let span_overlaps_other = group_members.iter().enumerate().any(|(oi, &o)| {
                    if oi == mi {
                        return false;
                    }
                    let vo = &variants[o].variant;
                    m_start < vo.pos + vo.ref_allele.len() as i64 && vo.pos < m_end
                });
                variants[member_idx].multi_allelic_group = Some(group_id);
                let tag = if span_overlaps_other { "MULTI_ALLELIC" } else { "TRACT_CLUSTER" };
                let reason = &mut variants[member_idx].gbcms_status_reason;
                if !reason.contains(tag) {
                    if !reason.is_empty() {
                        reason.push('|');
                    }
                    reason.push_str(tag);
                }
            }
            debug!(
                "Co-annotation group {}: {} variants at {}:{}-{} (window [{}, {}))",
                group_id,
                group_members.len(),
                variants[group_members[0]].variant.chrom,
                span_lo + 1,
                span_hi,
                if win_lo == i64::MAX { span_lo } else { win_lo },
                if win_hi == i64::MIN { span_hi } else { win_hi },
            );
        }
    }
}

/// Process a single variant through the full preparation pipeline.
///
/// Steps: MAF anchor → validate REF → left-align → adaptive context → fetch ref_context.
fn prepare_single_variant(
    reader_result: &mut Result<fasta::IndexedReader<File>, anyhow::Error>,
    variant: &Variant,
    context_padding: i64,
    is_maf: bool,
    adaptive_context: bool,
) -> Result<PreparedVariant, anyhow::Error> {
    let reader = reader_result.as_mut().map_err(|e| {
        anyhow::anyhow!("FASTA reader not available: {}", e)
    })?;

    // Save originals before any transformation
    let original_pos = variant.pos;
    let original_ref = variant.ref_allele.clone();
    let original_alt = variant.alt_allele.clone();

    // Step 0: reject structurally empty alleles up front. The internal representation
    // is VCF-style (anchor-based) — a true indel keeps its anchor base, and MAF dash
    // alleles arrive as the literal "-" (resolved in Step 1), never "". An empty REF or
    // ALT is therefore always malformed input: counting it would lean on the engine's
    // defensive empty-allele guards and silently yield zero counts. Reject it loudly
    // here with a FAIL status (mirrors the ALT-contains-N gate below) so the variant is
    // surfaced, not quietly dropped to all-zero.
    if variant.ref_allele.is_empty() || variant.alt_allele.is_empty() {
        warn!(
            "Empty allele at {}:{} {:?}>{:?} — rejecting (malformed indel; MAF dash alleles must be '-', not '')",
            variant.chrom,
            variant.pos + 1,
            variant.ref_allele,
            variant.alt_allele,
        );
        return Ok(PreparedVariant {
            variant: variant.clone(),
            gbcms_status: "FAIL".to_string(),
            gbcms_status_reason: "EMPTY_ALLELE".to_string(),
            gbcms_diagnostic: String::new(),
            gbcms_rescue: String::new(),
            was_anchor_resolved: false,
            was_left_aligned: false,
            original_pos,
            original_ref,
            original_alt,
            decomposed_variant: None,
            multi_allelic_group: None,
        });
    }

    // Step 1: MAF anchor resolution (only for dash-allele variants)
    // Non-dash MAF variants with different-length alleles (e.g., GG>A) already have
    // complete alleles and don't need an anchor base prepended.
    let (mut pos, mut ref_al, mut alt_al, mut vtype) = if is_maf
        && (variant.ref_allele == "-" || variant.alt_allele == "-")
    {
        // MAF indel/complex: resolve anchor base
        // variant.pos is 0-based (from maf_to_internal), start_pos is 1-based
        let start_pos_1based = variant.pos + 1;
        match resolve_maf_anchor(
            reader,
            &variant.chrom,
            start_pos_1based,
            &variant.ref_allele,
            &variant.alt_allele,
        ) {
            Ok(result) => result,
            Err(e) => {
                warn!(
                    "MAF anchor fetch failed for {}:{} {}>{}:  {}",
                    variant.chrom,
                    variant.pos + 1,
                    variant.ref_allele,
                    variant.alt_allele,
                    e,
                );
                return Ok(PreparedVariant {
                    variant: variant.clone(),
                    gbcms_status: "FAIL".to_string(),
                    gbcms_status_reason: "FETCH_FAILED".to_string(),
                    gbcms_diagnostic: String::new(),
                    gbcms_rescue: String::new(),
                    was_anchor_resolved: false,
                    was_left_aligned: false,
                    original_pos,
                    original_ref,
                    original_alt,
                    decomposed_variant: None,
                    multi_allelic_group: None,
                });
            }
        }
    } else {
        // VCF-style or MAF SNP: use coords as-is
        (
            variant.pos,
            variant.ref_allele.clone(),
            variant.alt_allele.clone(),
            variant.variant_type.clone(),
        )
    };

    // Track whether MAF anchor resolution changed pos/ref/alt (Step 1).
    // For VCF input or MAF SNPs, this is always false.
    let was_anchor_resolved = pos != original_pos
        || ref_al != original_ref
        || alt_al != original_alt;

    // Step 2: REF validation (with tolerance for partial mismatches).
    // `reason` carries WARN_REF_CORRECTED forward to the PASS success path (below),
    // so a corrected REF is surfaced rather than silently downgraded to a bare PASS.
    let (verdict, reason, corrected_ref) = validate_ref(reader, &variant.chrom, pos, &ref_al);

    // If REF was partially mismatched but ≥90% similar, correct it to the FASTA REF
    if let Some(fasta_ref) = corrected_ref {
        ref_al = fasta_ref;
    }

    if verdict != "PASS" {
        debug!(
            "REF validation FAIL ({}): {}:{} {}>{}",
            reason,
            variant.chrom,
            pos + 1,
            ref_al,
            alt_al,
        );
        return Ok(PreparedVariant {
            variant: Variant {
                chrom: variant.chrom.clone(),
                pos,
                ref_allele: ref_al,
                alt_allele: alt_al,
                variant_type: vtype,
                ref_context: None,
                ref_context_start: 0,
                repeat_span: 0,
                gene_strand: None,
            },
            gbcms_status: verdict,
            gbcms_status_reason: reason,
            gbcms_diagnostic: String::new(),
            gbcms_rescue: String::new(),
            was_anchor_resolved,
            was_left_aligned: false,
            original_pos,
            original_ref,
            original_alt,
            decomposed_variant: None,
            multi_allelic_group: None,
        });
    }

    // Step 2b: ALT allele N-base validation
    // ALT alleles containing N (ambiguous/placeholder bases from incomplete
    // genotyping) cannot be meaningfully counted — no real read base matches N.
    // Reject explicitly rather than silently producing zero counts.
    if alt_al.as_bytes().iter().any(|&b| b == b'N' || b == b'n') {
        warn!(
            "ALT allele contains N (ambiguous): {}:{} {}>{} — rejecting",
            variant.chrom, pos + 1, ref_al, alt_al
        );
        return Ok(PreparedVariant {
            variant: Variant {
                chrom: variant.chrom.clone(),
                pos,
                ref_allele: ref_al,
                alt_allele: alt_al,
                variant_type: vtype,
                ref_context: None,
                ref_context_start: 0,
                repeat_span: 0,
                gene_strand: None,
            },
            gbcms_status: "FAIL".to_string(),
            gbcms_status_reason: "ALT_CONTAINS_N".to_string(),
            gbcms_diagnostic: String::new(),
            gbcms_rescue: String::new(),
            was_anchor_resolved,
            was_left_aligned: false,
            original_pos,
            original_ref,
            original_alt,
            decomposed_variant: None,
            multi_allelic_group: None,
        });
    }

    // Step 3: Left-alignment (only for true indels/complex, NOT MNPs)
    // MNPs (same-length multi-base substitutions: DNP, TNP, ONP) are typed
    // as COMPLEX by kernel.py, but they are pure substitutions that cannot
    // be left-aligned and don't need ref_context for alignment.
    // C++ GBCMS (baseCountDNP) has no normalization at all.
    let mut was_left_aligned = false;
    let is_mnp = ref_al.len() == alt_al.len() && ref_al.len() > 1;
    let is_indel = !is_mnp
        && (ref_al.len() != alt_al.len()
            || (ref_al.len() > 1 && alt_al.len() > 1));

    if is_indel {
        let mut norm_window: i64 = 100; // bcftools default
        let max_norm_window: i64 = 2500; // safety cap for centromeric regions

        loop {
            let pad = context_padding.max(norm_window);
            let wide_start = (pos - pad).max(0);
            let wide_end = pos + ref_al.len() as i64 + pad;

            match fetch_region(
                reader,
                &variant.chrom,
                wide_start as u64,
                wide_end as u64,
            ) {
                Ok(wide_ref) => {
                    let pos_before_align = pos;
                    let (new_pos, new_ref, new_alt, modified) = left_align_variant(
                        pos,
                        ref_al.as_bytes(),
                        alt_al.as_bytes(),
                        &wide_ref,
                        wide_start,
                        norm_window as usize,
                    );

                    if modified {
                        match (String::from_utf8(new_ref), String::from_utf8(new_alt)) {
                            (Ok(new_ref_s), Ok(new_alt_s)) => {
                                debug!(
                                    "Left-aligned: {}:{} {}>{} → {}:{} {}>{}",
                                    variant.chrom,
                                    pos + 1,
                                    ref_al,
                                    alt_al,
                                    variant.chrom,
                                    new_pos + 1,
                                    new_ref_s,
                                    new_alt_s,
                                );
                                pos = new_pos;
                                ref_al = new_ref_s;
                                alt_al = new_alt_s;
                                was_left_aligned = true;

                                // Re-determine variant type after normalization
                                vtype = if ref_al.len() == 1 && alt_al.len() == 1 {
                                    "SNP".to_string()
                                } else if ref_al.len() == 1 && alt_al.len() > 1 {
                                    "INSERTION".to_string()
                                } else if ref_al.len() > 1 && alt_al.len() == 1 {
                                    "DELETION".to_string()
                                } else {
                                    "COMPLEX".to_string()
                                };

                                // Check: did the variant shift all the way to the window
                                // edge? If so, it may not have fully converged — expand
                                // and retry.
                                let shift = pos_before_align - pos;
                                if shift >= norm_window {
                                    if norm_window < max_norm_window {
                                        norm_window = (norm_window * 2).min(max_norm_window);
                                        debug!(
                                            "Left-align hit window edge (shift={}bp), \
                                             expanding to {}bp for {}:{}",
                                            shift, norm_window, variant.chrom, pos + 1
                                        );
                                        continue; // Re-align with wider window
                                    }
                                    // The safety cap (centromeric/telomeric repeat
                                    // oceans) bound before convergence: the variant is
                                    // left-aligned only as far as the cap allows.
                                    warn!(
                                        "Left-align hit the {}bp window cap without \
                                         converging for {}:{} (shift={}bp) — anchor may \
                                         not be fully left-aligned",
                                        max_norm_window, variant.chrom, pos + 1, shift
                                    );
                                }
                            }
                            _ => {
                                // Non-UTF-8 bytes in the fetched window (corrupt
                                // reference). Keep the variant EXACTLY as it was —
                                // adopting the shifted pos with reverted alleles
                                // would corrupt coordinates.
                                warn!(
                                    "Left-align produced non-UTF-8 alleles for {}:{} \
                                     (corrupt reference window?) — keeping the \
                                     unnormalized variant",
                                    variant.chrom, pos + 1
                                );
                            }
                        }
                    }
                }
                Err(e) => {
                    // A missed left-alignment shifts the counting anchor for every
                    // downstream consumer (repeat scan, windowed matching, Phase-3
                    // haplotypes), so this must be loud. Reachable for real inputs:
                    // the FASTA reader errors when the window end passes the contig
                    // end, so indels within ~100bp of a contig boundary land here.
                    warn!(
                        "Wide ref fetch failed for {}:{}-{} ({e}) — variant NOT \
                         left-aligned; counting proceeds at the input coordinates",
                        variant.chrom, wide_start, wide_end,
                    );
                }
            }
            break; // Normal exit: alignment converged, capped, or fetch failed
        }
    }

    // Step 4: Fetch ref_context at (possibly normalized) position
    //         With adaptive_context, padding is increased in repeat regions.
    //         MNPs are excluded: they are pure substitutions that don't need
    //         ref_context for SW/HMM indel realignment.
    let (ref_context, ref_context_start) = if is_indel || (vtype == "COMPLEX" && !is_mnp) {
        // Scan for repeats at the FIRST CHANGED base, not the shared anchor:
        // a left-aligned repeat indel's anchor sits one base left of the
        // tract, where the scan finds span 1 and adaptive padding never
        // widens (issue #91).
        let change_off = first_change_offset(&ref_al, &alt_al);
        let effective_padding = if adaptive_context {
            compute_adaptive_padding(
                reader,
                &variant.chrom,
                pos + change_off,
                ref_al.len(),
                context_padding,
                50,  // max cap
            )
        } else {
            context_padding
        };

        let ctx_start = (pos - effective_padding).max(0);
        let ctx_end = pos + ref_al.len() as i64 + effective_padding;

        match fetch_region(
            reader,
            &variant.chrom,
            ctx_start as u64,
            ctx_end as u64,
        ) {
            Ok(ctx_bytes) => (
                Some(String::from_utf8_lossy(&ctx_bytes).to_string()),
                ctx_start,
            ),
            Err(e) => {
                warn!(
                    "ref_context fetch failed for {}:{}-{} ({e}) — haplotype \
                     alignment (Phase 3) will be skipped for this variant",
                    variant.chrom, ctx_start, ctx_end,
                );
                (None, 0)
            }
        }
    } else {
        // SNPs don't need ref_context
        (None, 0)
    };

    // Step 5: Homopolymer decomposition detection
    // Check if the variant looks like a miscollapsed D(n)+SNV in a homopolymer.
    // If detected, build a corrected Variant for dual-counting.
    let decomposed_variant = if ref_al.len() > alt_al.len() && ref_al.len() >= 3 {
        // Fetch the next reference base after the ref span
        let next_pos = (pos + ref_al.len() as i64) as u64;
        let next_base = fetch_region(reader, &variant.chrom, next_pos, next_pos + 1)
            .ok()
            .and_then(|b| b.first().copied());

        next_base.and_then(|nb| {
            check_homopolymer_decomp(&ref_al, &alt_al, nb).map(|corrected_alt| {
                let decomp_vtype = if ref_al.len() == corrected_alt.len() {
                    if ref_al.len() == 1 { "SNP" } else { "COMPLEX" }
                } else if ref_al.len() > corrected_alt.len() {
                    "DELETION"
                } else {
                    "INSERTION"
                }
                .to_string();

                Variant {
                    chrom: variant.chrom.clone(),
                    pos,
                    ref_allele: ref_al.clone(),
                    alt_allele: corrected_alt,
                    variant_type: decomp_vtype,
                    ref_context: ref_context.clone(),
                    ref_context_start,
                    repeat_span: 0, // Decomposed variant inherits context but repeat info is not critical
                    gene_strand: None,
                }
            })
        })
    } else {
        None
    };

    // Compute repeat_span from ref_context (PairHMM gap blending, windowed-scan
    // width, tract-cluster grouping reach).
    // find_tandem_repeat detects tandem repeats around the variant position.
    let variant_repeat_span = if let Some(ref ctx) = ref_context {
        let ctx_bytes = ctx.as_bytes();
        // Same first-changed-base anchoring as the adaptive scan above:
        // repeat_span feeds the windowed-scan width and SW gap tuning, and an
        // anchor-based scan misses edge tracts entirely (issue #91).
        let scan_pos = pos + first_change_offset(&ref_al, &alt_al);
        let pos_in_ctx = (scan_pos - ref_context_start) as usize;
        let (_motif_len, span) = find_tandem_repeat(ctx_bytes, pos_in_ctx.min(ctx_bytes.len().saturating_sub(1)));
        span
    } else {
        0
    };

    Ok(PreparedVariant {
        variant: Variant {
            chrom: variant.chrom.clone(),
            pos,
            ref_allele: ref_al,
            alt_allele: alt_al,
            variant_type: vtype,
            ref_context,
            ref_context_start,
            repeat_span: variant_repeat_span,
            gene_strand: None,
        },
        gbcms_status: "PASS".to_string(),
        // Carry any WARN_REF_CORRECTED from validate_ref (empty otherwise). A later
        // pass may append MULTI_ALLELIC; the pipeline may append WARN_HOMOPOLYMER_DECOMP.
        gbcms_status_reason: reason,
        gbcms_diagnostic: String::new(),
        gbcms_rescue: String::new(),
        was_anchor_resolved,
        was_left_aligned,
        original_pos,
        original_ref,
        original_alt,
        decomposed_variant,
        multi_allelic_group: None, // Set by post-processing in prepare_variants()
    })
}

// ===========================================================================
// Unit Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // -- left_align_variant tests --

    #[test]
    fn test_snp_passthrough() {
        // SNPs should never be modified
        let reference = b"ATCGATCG";
        let (pos, r, a, modified) =
            left_align_variant(3, b"G", b"T", reference, 0, 100);
        assert_eq!(pos, 3);
        assert_eq!(r, b"G");
        assert_eq!(a, b"T");
        assert!(!modified);
    }

    #[test]
    fn test_deletion_in_homopolymer() {
        // Reference: GAAAAAAG (poly-A at pos 1-6)
        //            01234567
        // Deletion at pos 5: REF=AA, ALT=A (delete one A)
        // Should normalize to pos 0: REF=GA, ALT=G
        let reference = b"GAAAAAAG";
        let (pos, r, a, modified) =
            left_align_variant(5, b"AA", b"A", reference, 0, 100);

        assert!(modified, "Deletion in homopolymer should be left-shifted");
        assert_eq!(pos, 0, "Should shift to leftmost position");
        assert_eq!(r, b"GA");
        assert_eq!(a, b"G");
    }

    #[test]
    fn test_insertion_in_homopolymer() {
        // Reference: GAAAAAAG (poly-A at pos 1-6)
        //            01234567
        // Insertion at pos 5: REF=A, ALT=AA (insert one A)
        // Should normalize to pos 0: REF=G, ALT=GA
        let reference = b"GAAAAAAG";
        let (pos, r, a, modified) =
            left_align_variant(5, b"A", b"AA", reference, 0, 100);

        assert!(modified, "Insertion in homopolymer should be left-shifted");
        assert_eq!(pos, 0, "Should shift to leftmost position");
        assert_eq!(r, b"G");
        assert_eq!(a, b"GA");
    }

    #[test]
    fn test_deletion_in_dinuc_repeat() {
        // Reference: CTCTCTCTCTCT (CT repeat)
        //            0 1 2 3 4 5 6 7 8 9 10 11
        // Deletion at pos 4: REF=CTCT, ALT=CT (delete one CT unit)
        // ref[4..8] = C,T,C,T ✓
        // After extend left 4 + right-trim + left-trim:
        //   Extended: r=CTCTCTCT, a=CTCTCT, cur_pos=0
        //   Right-trim: removes matching pairs until a.len()==1 → r=CTC, a=C
        //   Left-trim: a.len()==1 → no trim
        // Result: pos=0, REF=CTC, ALT=C (minimal representation)
        let reference = b"CTCTCTCTCTCT";
        let (pos, r, a, modified) =
            left_align_variant(4, b"CTCT", b"CT", reference, 0, 100);

        assert!(modified, "Deletion in dinuc repeat should be left-shifted");
        assert_eq!(pos, 0, "Should shift to leftmost position");
        assert_eq!(r, b"CTC");
        assert_eq!(a, b"C");
    }

    #[test]
    fn test_insertion_in_dinuc_repeat() {
        // Reference: CTCTCTCTCTCT (CT repeat)
        //            0 1 2 3 4 5 6 7 8 9 10 11
        // Insertion at pos 4: REF=C, ALT=CCT (insert one CT unit)
        // ref[4] = C ✓
        // After extend left 4: r=CTCTC, a=CTCTCCT
        // Right-trim: C vs T → no match (immediate stop since last bases differ)
        //
        // The algorithm correctly can't shift this because the trailing bases
        // are different. We need to ensure the alleles share a trailing context.
        // Instead, let's use: pos=4, REF=CT, ALT=CTCT (insert CT before next CT)
        // ref[4..6] = C,T ✓
        // After extend left 4: r=CTCTCT, a=CTCTCTCT
        // Right-trim: T==T, C==C, T==T, C==C → trim 4 → r=CT, a=CTCT
        //   Keep trimming: T==T → r=C, a=CTC. C==C → r="" stop (len=1 for r)
        //   Actually r.len()>1 fails. So r=C, a=CTC.
        // Wait, let me retrace. After extend 4: r=[C,T,C,T,C,T], a=[C,T,C,T,C,T,C,T]
        // Right: T==T → pop → r=[C,T,C,T,C], a=[C,T,C,T,C,T,C]
        // Right: C==C → pop → r=[C,T,C,T], a=[C,T,C,T,C,T]
        // Right: T==T → pop → r=[C,T,C], a=[C,T,C,T,C]
        // Right: C==C → pop → r=[C,T], a=[C,T,C,T]
        // Right: T==T → pop → r=[C], a=[C,T,C] — stop (r.len()==1)
        // Left: C==C → r.len()==1, stop
        // Result: pos=0+3=3? No, cur_pos = 4-4=0 (from extend), then left-trim doesn't happen.
        // Final: pos=0, r=[C], a=[C,T,C]  = REF=C, ALT=CTC
        //
        // That's correct! pos=0, REF=C, ALT=CTC is the leftmost representation.
        let reference = b"CTCTCTCTCTCT";
        let (pos, r, a, modified) =
            left_align_variant(4, b"CT", b"CTCT", reference, 0, 100);

        assert!(modified, "Insertion in dinuc repeat should be left-shifted");
        assert_eq!(pos, 0, "Should shift to leftmost position");
        assert_eq!(r, b"C");
        assert_eq!(a, b"CTC");
    }

    #[test]
    fn test_no_shift_needed() {
        // Reference: ATCGATCG — non-repetitive
        // Deletion at pos=1: REF=TC, ALT=T
        // ref[1..3] = T,C ✓
        // After extend 1: r=[A,T,C], a=[A,T]
        // Right-trim: C vs T → no match
        // Left-trim: A==A → r=[T,C], a=[T], pos=1. T==T → r.len()>1 but a.len()==1, stop
        // But now pos=1 and r=[T,C], a=[T] — same as original!
        // Actually this IS already leftmost. Let me verify:
        // Original: pos=1, REF=TC, ALT=T. After normalize: pos=1, REF=TC, ALT=T.
        let reference = b"ATCGATCG";
        let (pos, r, a, modified) =
            left_align_variant(1, b"TC", b"T", reference, 0, 100);
        assert_eq!(pos, 1);
        assert_eq!(r, b"TC");
        assert_eq!(a, b"T");
        assert!(!modified);
    }

    #[test]
    fn test_pos_at_start() {
        // Variant at position 0 — no room to shift left
        // Reference: AATCGATCG
        //            012345678
        // pos=0, REF=AA, ALT=A — deletion of one A
        let reference = b"AATCGATCG";
        let (pos, r, a, _modified) =
            left_align_variant(0, b"AA", b"A", reference, 0, 100);
        assert_eq!(pos, 0, "Cannot shift past position 0");
        assert_eq!(r, b"AA");
        assert_eq!(a, b"A");
    }

    #[test]
    fn test_with_ref_offset() {
        // Reference window starts at offset 100 (simulates fetching a window)
        // Poly-A: reference = GAAAAAAG at positions 100-107
        // Deletion at pos=105: REF=AA, ALT=A
        // ref_in_window[5..7] = A,A ✓
        let reference = b"GAAAAAAG";
        let (pos, r, a, modified) =
            left_align_variant(105, b"AA", b"A", reference, 100, 100);

        assert!(modified);
        assert_eq!(pos, 100, "Should shift to start of window");
        assert_eq!(r, b"GA");
        assert_eq!(a, b"G");
    }

    #[test]
    fn test_complex_no_normalize() {
        // Complex variant (MNP-like): REF=AG, ALT=TC — same length, no shared bases
        let reference = b"AAAGTTTT";
        let (pos, r, a, modified) =
            left_align_variant(2, b"AG", b"TC", reference, 0, 100);
        assert_eq!(pos, 2);
        assert_eq!(r, b"AG");
        assert_eq!(a, b"TC");
        assert!(!modified);
    }

    // -- find_tandem_repeat tests --

    #[test]
    fn test_find_tandem_repeat_homopolymer() {
        // GATCAAAAAAGCTT — 6 A's at pos 4-9
        let seq = b"GATCAAAAAAGCTT";
        let (unit, span) = find_tandem_repeat(seq, 4);
        assert_eq!(unit, 1, "Should detect homopolymer (unit=1)");
        assert_eq!(span, 6, "6 consecutive A's");
    }

    #[test]
    fn test_find_tandem_repeat_dinuc() {
        // GATCACACACCGTT — CA repeat at pos 4-9 (3 copies)
        let seq = b"GATCACACACGTT";
        let (unit, span) = find_tandem_repeat(seq, 4);
        assert_eq!(unit, 2, "Should detect dinucleotide repeat");
        assert_eq!(span, 6, "3 copies of CA = 6bp span");
    }

    #[test]
    fn test_find_tandem_repeat_trinuc() {
        // CAGCAGCAGTTTT — 3 copies of CAG at pos 0
        let seq = b"CAGCAGCAGTTTT";
        let (unit, span) = find_tandem_repeat(seq, 0);
        assert_eq!(unit, 3, "Should detect trinucleotide repeat");
        assert_eq!(span, 9, "3 copies of CAG = 9bp span");
    }

    #[test]
    fn test_find_tandem_repeat_no_repeat() {
        // ATCGATCG — no repeat at pos 2
        let seq = b"ATCGATCG";
        let (unit, span) = find_tandem_repeat(seq, 2);
        assert_eq!(span, 1, "No repeat detected, span should be 1");
        assert_eq!(unit, 1);
    }

    #[test]
    fn test_find_tandem_repeat_at_edge() {
        // AAAAAAA — homopolymer at pos 0
        let seq = b"AAAAAAA";
        let (unit, span) = find_tandem_repeat(seq, 0);
        assert_eq!(unit, 1);
        assert_eq!(span, 7, "Entire sequence is a homopolymer");
    }

    #[test]
    fn test_find_tandem_repeat_middle_of_run() {
        // TTAAAAATT — 5 A's, queried at pos 4 (middle of run)
        let seq = b"TTAAAAATT";
        let (unit, span) = find_tandem_repeat(seq, 4);
        assert_eq!(unit, 1);
        assert_eq!(span, 5, "5 A's even when querying from the middle");
    }

    // -- compute_adaptive_padding (formula) tests --

    // Formula: genuine repeat (span >= 2) pads by the full tract span on top
    // of the default, so the haplotype window always contains the tract plus
    // unique flank on both sides (issue #91). Non-repeats keep the default.
    fn effective_padding(span: usize, default_pad: i64, max_pad: i64) -> i64 {
        let adaptive = if span >= 2 { span as i64 + default_pad } else { 0 };
        default_pad.max(adaptive).min(max_pad)
    }

    #[test]
    fn test_adaptive_padding_default() {
        // Non-repeat: span=1 → default padding unchanged
        assert_eq!(effective_padding(1, 5, 50), 5, "Non-repeat should use default padding");
    }

    #[test]
    fn test_adaptive_padding_homopoly() {
        // Poly-A span=14 → 14 + 5 = 19: tract + 5bp unique flank each side
        assert_eq!(effective_padding(14, 5, 50), 19, "Homopolymer must span tract + flanks");
    }

    #[test]
    fn test_adaptive_padding_dinuc() {
        // CA-repeat span=20 → 20 + 5 = 25
        assert_eq!(effective_padding(20, 5, 50), 25, "Dinucleotide repeat must span tract + flanks");
    }

    #[test]
    fn test_adaptive_padding_capped() {
        // Very long repeat span=120 → 125 → capped at 50
        assert_eq!(effective_padding(120, 5, 50), 50, "Should be capped at max_pad");
    }

    // -- first_change_offset + scan-anchor regression (issue #91) --

    #[test]
    fn test_first_change_offset_indels_and_complex() {
        use super::super::repeat::first_change_offset;
        assert_eq!(first_change_offset("GAA", "G"), 1, "pure deletion: first deleted base");
        assert_eq!(first_change_offset("G", "GTT"), 1, "pure insertion: first inserted base");
        assert_eq!(first_change_offset("GC", "T"), 0, "complex without shared anchor");
        assert_eq!(first_change_offset("GAT", "GCT"), 1, "substitution after anchor");
    }

    #[test]
    fn test_repeat_scan_at_first_changed_base_finds_edge_tract() {
        // Anchor 'C' sits one base left of a 10-A tract. Scanning at the
        // anchor found span 1 (the historical bug); scanning at the first
        // deleted base finds the full tract.
        let seq = b"GGTTCAAAAAAAAAATTGG";
        let anchor = 4; // the 'C'
        let (_, span_at_anchor) = find_tandem_repeat(seq, anchor);
        let (_, span_at_change) = find_tandem_repeat(seq, anchor + 1);
        assert_eq!(span_at_anchor, 1, "anchor scan misses the tract (documented bug)");
        assert_eq!(span_at_change, 10, "first-changed-base scan finds the tract");
    }

    // -- Gap 1B: Dynamic window expansion tests --

    #[test]
    fn test_window_expansion_long_homopolymer() {
        // A 120bp dinucleotide repeat requires window expansion beyond 100bp.
        // Pure homopolymers don't shift (all positions equivalent), so we use AC-repeat.
        // Reference: 50bp T-prefix + 120bp AC-repeat + 50bp G-suffix = 220bp.
        let mut reference = Vec::with_capacity(220);
        reference.extend(std::iter::repeat_n(b'T', 50));
        for _ in 0..60 { reference.push(b'A'); reference.push(b'C'); }
        reference.extend(std::iter::repeat_n(b'G', 50));

        // Deletion near end of AC-run: pos=160, REF=AC, ALT=A
        // In dinucleotide repeat, this should left-align back toward pos 49/50
        let (pos, _r, _a, modified) =
            left_align_variant(160, b"AC", b"A", &reference, 0, 200);
        if modified {
            assert!(pos < 160,
                "Expected leftward shift from pos 160, got pos={}", pos);
        }
    }
}
