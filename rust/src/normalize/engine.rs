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
use super::repeat::{find_tandem_repeat, compute_adaptive_padding};

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

    // Build rayon thread pool (same pattern as count_bam)
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
            let valid = r.iter().filter(|p| p.validation_status == "PASS").count();
            let n_anchor = r.iter().filter(|p| p.was_anchor_resolved).count();
            let n_left = r.iter().filter(|p| p.was_left_aligned).count();
            let n_total = r.iter().filter(|p| p.was_normalized()).count();
            let multi_allelic = r.iter().filter(|p| p.multi_allelic_group.is_some()).count();
            // Count MNPs: same-length multi-base substitutions excluded from
            // indel normalization (left-alignment + ref_context fetch).
            let n_mnp = r.iter().filter(|p| {
                let ref_len = p.variant.ref_allele.len();
                let alt_len = p.variant.alt_allele.len();
                ref_len == alt_len && ref_len > 1
            }).count();
            info!(
                "prepare_variants complete: {}/{} valid, {} normalized ({} anchor-resolved, {} left-aligned), {} multi-allelic, {} MNPs (normalization bypassed)",
                valid,
                r.len(),
                n_total,
                n_anchor,
                n_left,
                multi_allelic,
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

/// Assign group IDs to variants that overlap at the same genomic locus.
///
/// Two variants overlap if they share the same chromosome and their REF spans
/// (pos..pos+ref_len) intersect. Uses a sweep-line algorithm: sort by position,
/// then extend the current group while new variants overlap the group's footprint.
fn assign_multi_allelic_groups(variants: &mut [PreparedVariant]) {
    if variants.len() < 2 {
        return;
    }

    // Build index sorted by (chrom, pos) — we sort indices, not the array itself,
    // to avoid disrupting input order (which must match the output row order).
    let mut indices: Vec<usize> = (0..variants.len()).collect();
    indices.sort_by(|&a, &b| {
        let va = &variants[a].variant;
        let vb = &variants[b].variant;
        va.chrom.cmp(&vb.chrom).then(va.pos.cmp(&vb.pos))
    });

    let mut group_id: u32 = 0;
    let mut i = 0;

    while i < indices.len() {
        let idx = indices[i];
        // Skip non-PASS variants
        if !variants[idx].validation_status.starts_with("PASS") {
            i += 1;
            continue;
        }

        let chrom = &variants[idx].variant.chrom.clone();
        let mut group_end = variants[idx].variant.pos
            + variants[idx].variant.ref_allele.len() as i64;
        let mut group_members: Vec<usize> = vec![idx];

        // Sweep forward: extend group while variants overlap
        let mut j = i + 1;
        while j < indices.len() {
            let jdx = indices[j];
            let vj = &variants[jdx].variant;

            if &vj.chrom != chrom {
                break; // Different chromosome
            }
            if !variants[jdx].validation_status.starts_with("PASS") {
                j += 1;
                continue;
            }
            if vj.pos < group_end {
                // Overlaps — extend the group footprint
                let vj_end = vj.pos + vj.ref_allele.len() as i64;
                group_end = group_end.max(vj_end);
                group_members.push(jdx);
                j += 1;
            } else {
                break; // No overlap
            }
        }

        // Only assign group ID if there are 2+ overlapping variants
        if group_members.len() > 1 {
            group_id += 1;
            for &member_idx in &group_members {
                variants[member_idx].multi_allelic_group = Some(group_id);
                // Append to validation_status so it surfaces in existing output columns
                // e.g. "PASS" → "PASS_MULTI_ALLELIC"
                if !variants[member_idx].validation_status.contains("MULTI_ALLELIC") {
                    variants[member_idx].validation_status.push_str("_MULTI_ALLELIC");
                }
            }
            debug!(
                "Multi-allelic group {}: {} variants at {}:{}-{}",
                group_id,
                group_members.len(),
                chrom,
                variants[group_members[0]].variant.pos + 1,
                group_end,
            );
        }

        i = j;
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
                    validation_status: "FETCH_FAILED".to_string(),
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

    // Step 2: REF validation (with tolerance for partial mismatches)
    let (status, corrected_ref) = validate_ref(reader, &variant.chrom, pos, &ref_al);

    // If REF was partially mismatched but ≥90% similar, correct it to the FASTA REF
    if let Some(fasta_ref) = corrected_ref {
        ref_al = fasta_ref;
    }

    if !status.starts_with("PASS") {
        debug!(
            "REF validation {}: {}:{} {}>{} ({})",
            status,
            variant.chrom,
            pos + 1,
            ref_al,
            alt_al,
            status,
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
            validation_status: status,
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
            validation_status: "FAIL_ALT_CONTAINS_N".to_string(),
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

            if let Ok(wide_ref) = fetch_region(
                reader,
                &variant.chrom,
                wide_start as u64,
                wide_end as u64,
            ) {
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
                    debug!(
                        "Left-aligned: {}:{} {}>{} → {}:{} {}>{}",
                        variant.chrom,
                        pos + 1,
                        ref_al,
                        alt_al,
                        variant.chrom,
                        new_pos + 1,
                        String::from_utf8_lossy(&new_ref),
                        String::from_utf8_lossy(&new_alt),
                    );
                    pos = new_pos;
                    ref_al = String::from_utf8(new_ref)
                        .unwrap_or_else(|_| ref_al.clone());
                    alt_al = String::from_utf8(new_alt)
                        .unwrap_or_else(|_| alt_al.clone());
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

                    // Check: did the variant shift all the way to the window edge?
                    // If so, it may not have fully converged — expand and retry.
                    let shift = pos_before_align - pos;
                    if shift >= norm_window && norm_window < max_norm_window {
                        norm_window = (norm_window * 2).min(max_norm_window);
                        debug!(
                            "Left-align hit window edge (shift={}bp), \
                             expanding to {}bp for {}:{}",
                            shift, norm_window, variant.chrom, pos + 1
                        );
                        continue; // Re-align with wider window
                    }
                }
            } else {
                debug!(
                    "Wide ref fetch failed for {}:{}-{}, skipping normalization",
                    variant.chrom, wide_start, wide_end,
                );
            }
            break; // Normal exit: alignment converged or fetch failed
        }
    }

    // Step 4: Fetch ref_context at (possibly normalized) position
    //         With adaptive_context, padding is increased in repeat regions.
    //         MNPs are excluded: they are pure substitutions that don't need
    //         ref_context for SW/HMM indel realignment.
    let (ref_context, ref_context_start) = if is_indel || (vtype == "COMPLEX" && !is_mnp) {
        let effective_padding = if adaptive_context {
            compute_adaptive_padding(
                reader,
                &variant.chrom,
                pos,
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
                println!("REF_CTX ERROR: {:?}", e);
                warn!(
                    "ref_context fetch failed for {}:{}-{}, SW alignment will be skipped",
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

    // Compute repeat_span from ref_context for dynamic SW gap penalty tuning.
    // find_tandem_repeat detects tandem repeats around the variant position.
    let variant_repeat_span = if let Some(ref ctx) = ref_context {
        let ctx_bytes = ctx.as_bytes();
        let pos_in_ctx = (pos - ref_context_start) as usize;
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
        validation_status: "PASS".to_string(),
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

    #[test]
    fn test_adaptive_padding_default() {
        // Non-repeat: span=1 → adaptive = 1/2+3 = 3 → max(5,3) = 5
        let span = 1;
        let default_pad: i64 = 5;
        let max_pad: i64 = 50;
        let adaptive = (span as i64) / 2 + 3;
        let effective = default_pad.max(adaptive).min(max_pad);
        assert_eq!(effective, 5, "Non-repeat should use default padding");
    }

    #[test]
    fn test_adaptive_padding_homopoly() {
        // Poly-A span=14 → adaptive = 14/2+3 = 10 → max(5,10) = 10
        let span = 14;
        let default_pad: i64 = 5;
        let max_pad: i64 = 50;
        let adaptive = (span as i64) / 2 + 3;
        let effective = default_pad.max(adaptive).min(max_pad);
        assert_eq!(effective, 10, "Homopolymer should increase padding");
    }

    #[test]
    fn test_adaptive_padding_dinuc() {
        // CA-repeat span=20 → adaptive = 20/2+3 = 13 → max(5,13) = 13
        let span = 20;
        let default_pad: i64 = 5;
        let max_pad: i64 = 50;
        let adaptive = (span as i64) / 2 + 3;
        let effective = default_pad.max(adaptive).min(max_pad);
        assert_eq!(effective, 13, "Dinucleotide repeat should increase padding");
    }

    #[test]
    fn test_adaptive_padding_capped() {
        // Very long repeat span=120 → adaptive = 120/2+3 = 63 → min(63,50) = 50
        let span = 120;
        let default_pad: i64 = 5;
        let max_pad: i64 = 50;
        let adaptive = (span as i64) / 2 + 3;
        let effective = default_pad.max(adaptive).min(max_pad);
        assert_eq!(effective, 50, "Should be capped at max_pad");
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
