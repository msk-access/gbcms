//! Tandem repeat detection and adaptive context padding.

use std::fs::File;

use bio::io::fasta;
use log::debug;

use super::fasta::fetch_region;

/// Detect the longest tandem repeat (motif 1–6bp) touching a given position.
///
/// Scans the provided sequence `seq` for tandem repeats whose footprint
/// includes `pos_in_seq`.  Returns `(motif_length, total_span)` for the
/// longest repeat found, or `(1, 1)` if no repeat is detected.
///
/// # Algorithm
///
/// For each candidate motif length *k* (1..=6):
/// 1. Extract the motif starting at `pos_in_seq`
/// 2. Extend left while preceding *k* bases match the motif
/// 3. Extend right while following *k* bases match the motif
/// 4. If the span covers ≥ 2 full copies of the motif, record it
///
/// The longest span (breaking ties in favour of larger motif) wins.
/// Offset of the first base that differs between REF and ALT alleles.
///
/// For anchored VCF-style indels (`GAA>G`, `G>GTT`) this is 1 — the first
/// deleted/inserted base — which is where a repeat scan must start. Scanning
/// at the shared anchor base instead (offset 0) sits one base LEFT of the
/// tract for every left-aligned repeat indel and finds span 1 (issue #91).
/// For complex alleles whose first bases already differ, the offset is 0.
pub(crate) fn first_change_offset(ref_al: &str, alt_al: &str) -> i64 {
    let r = ref_al.as_bytes();
    let a = alt_al.as_bytes();
    let shared = r.len().min(a.len());
    for i in 0..shared {
        if !r[i].eq_ignore_ascii_case(&a[i]) {
            return i as i64;
        }
    }
    shared as i64
}

pub(crate) fn find_tandem_repeat(seq: &[u8], pos_in_seq: usize) -> (usize, usize) {
    let len = seq.len();
    let mut best_unit = 1;
    let mut best_span = 1;

    for k in 1..=6usize {
        if pos_in_seq + k > len {
            break;
        }

        let motif = &seq[pos_in_seq..pos_in_seq + k];

        // Extend left: how many full motif copies precede pos_in_seq?
        let mut left_copies = 0;
        let mut cursor = pos_in_seq;
        while cursor >= k {
            if &seq[cursor - k..cursor] == motif {
                left_copies += 1;
                cursor -= k;
            } else {
                break;
            }
        }

        // Extend right: how many full motif copies follow pos_in_seq?
        let mut right_copies = 0;
        let mut cursor = pos_in_seq + k;
        while cursor + k <= len {
            if &seq[cursor..cursor + k] == motif {
                right_copies += 1;
                cursor += k;
            } else {
                break;
            }
        }

        let total_copies = left_copies + 1 + right_copies; // 1 = the copy at pos_in_seq
        let span = total_copies * k;

        // Must have ≥2 copies to count as a repeat
        if total_copies >= 2 && span > best_span {
            best_unit = k;
            best_span = span;
        }
    }

    (best_unit, best_span)
}

/// Compute adaptive context padding based on nearby tandem repeats.
///
/// Strategy:
/// 1. Fetch a wide scan window (±`scan_radius` bp) around the variant
/// 2. Run `find_tandem_repeat()` at the variant position within the window
/// 3. `effective = max(default_pad, repeat_span + default_pad)` for genuine
///    repeats (span >= 2), capped at `max_pad` — the window must contain the
///    whole tract plus unique flank on both sides (issue #91)
/// 4. Falls back to `default_pad` on FASTA fetch failure
///
/// # Arguments
/// * `reader` — FASTA reader (mutable)
/// * `chrom` — Chromosome name
/// * `pos` — 0-based variant position
/// * `ref_len` — Length of the REF allele
/// * `default_pad` — Minimum padding (from `--context-padding`)
/// * `max_pad` — Hard cap on adaptive padding
pub(crate) fn compute_adaptive_padding(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    pos: i64,
    ref_len: usize,
    default_pad: i64,
    max_pad: i64,
) -> i64 {
    let scan_radius: i64 = 30;
    let scan_start = (pos - scan_radius).max(0);
    let scan_end = pos + ref_len as i64 + scan_radius;

    let scan_seq = match fetch_region(reader, chrom, scan_start as u64, scan_end as u64) {
        Ok(s) => s,
        Err(_) => {
            // warn (not debug): a silent fallback here narrows the Phase-3
            // haplotype window in exactly the repeat contexts where a narrow
            // window mis-arbitrates wrong-length indels (issue #91).
            log::warn!(
                "Adaptive scan fetch failed for {}:{}-{}, using default padding {}",
                chrom, scan_start, scan_end, default_pad
            );
            return default_pad;
        }
    };

    // Position within the scan window
    let pos_in_scan = (pos - scan_start) as usize;

    let (motif_len, repeat_span) = find_tandem_repeat(&scan_seq, pos_in_scan);
    // A Phase-3 haplotype window can only distinguish tract lengths when it
    // contains the ENTIRE repeat tract plus unique flank on both sides, so a
    // genuine repeat (span >= 2) pads by the full span on top of the default.
    // The previous `span/2 + 3` could never cover an edge-anchored tract
    // (issue #91). Non-repeat positions keep the default unchanged.
    let adaptive = if repeat_span >= 2 {
        repeat_span as i64 + default_pad
    } else {
        0
    };
    let effective = default_pad.max(adaptive).min(max_pad);

    if effective > default_pad {
        debug!(
            "Adaptive padding for {}:{}: repeat motif={}bp span={}bp → padding {} (default {})",
            chrom, pos + 1, motif_len, repeat_span, effective, default_pad
        );
    }

    effective
}
