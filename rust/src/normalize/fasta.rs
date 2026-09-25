//! FASTA I/O, REF validation, and MAF anchor resolution.

use std::fs::File;

use bio::io::fasta;
use log::warn;

/// Fetch a region from FASTA under any name the contig goes by: as given, with
/// or without a `chr` prefix, and — for the mitochondrion — each of its
/// spellings, so the FASTA reconciles contigs as the BAM side does
/// (`normalize_contig`: `chrM` ~ `M` ~ `MT` ~ `chrMT`).
pub(crate) fn fetch_region(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    start: u64,
    end: u64,
) -> anyhow::Result<Vec<u8>> {
    let mut names = vec![chrom.to_string(), format!("chr{}", chrom)];
    if let Some(stripped) = chrom.strip_prefix("chr") {
        names.push(stripped.to_string());
    }
    if crate::shared::contig::normalize_contig(chrom) == "MT" {
        names.extend(
            crate::shared::contig::MITO_SPELLINGS
                .iter()
                .filter(|m| !names.iter().any(|n| n == *m))
                .map(|m| m.to_string())
                .collect::<Vec<_>>(),
        );
    }

    let mut buf = Vec::new();
    for name in &names {
        if reader.fetch(name, start, end).is_ok() {
            buf.clear();
            reader.read(&mut buf)?;
            if !buf.is_empty() {
                return Ok(buf);
            }
        }
    }

    anyhow::bail!("FASTA fetch failed for {}:{}-{}", chrom, start, end)
}

/// Fetch a single base from the reference, delegating to `fetch_region`.
pub(crate) fn fetch_single_base(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    pos_0based: i64,
) -> anyhow::Result<u8> {
    let pos = pos_0based as u64;
    let buf = fetch_region(reader, chrom, pos, pos + 1)?;
    buf.first()
        .map(|b| b.to_ascii_uppercase())
        .ok_or_else(|| anyhow::anyhow!(
            "FASTA fetch returned empty for {}:{}", chrom, pos_0based
        ))
}

/// Validate that the REF allele matches the reference genome.
///
/// Case-insensitive comparison with tolerance for partial mismatches.
/// Tries both the given chromosome name and with/without "chr" prefix.
///
/// # Returns
/// Tuple of `(verdict, reason, Option<fasta_ref>)`, where `verdict` is `"PASS"` or
/// `"FAIL"` and `reason` is the (possibly empty) reason tag:
/// - `("PASS", "", None)` — exact match
/// - `("PASS", "WARN_REF_CORRECTED", Some(fasta_ref))` — ≥90% match, corrected
/// - `("FAIL", "REF_MISMATCH", None)` — <90% match, rejected
/// - `("FAIL", "FETCH_FAILED", None)` — region could not be fetched
pub(crate) fn validate_ref(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    pos_0based: i64,
    ref_allele: &str,
) -> (String, String, Option<String>) {
    let ref_len = ref_allele.len() as u64;
    let pos = pos_0based as u64;

    let fetch_result = fetch_region(reader, chrom, pos, pos + ref_len);
    match fetch_result {
        Ok(ref_seq) => {
            // Fast path: exact match
            if ref_seq.eq_ignore_ascii_case(ref_allele.as_bytes()) {
                ("PASS".to_string(), String::new(), None)
            } else {
                // Compute similarity: count matching bases (case-insensitive)
                let max_len = ref_seq.len().max(ref_allele.len());
                if max_len == 0 {
                    return ("FAIL".to_string(), "REF_MISMATCH".to_string(), None);
                }
                let matches = ref_seq
                    .iter()
                    .zip(ref_allele.as_bytes())
                    .filter(|(a, b)| a.eq_ignore_ascii_case(b))
                    .count();
                let similarity = matches as f64 / max_len as f64;

                if similarity >= 0.90 {
                    let fasta_ref = String::from_utf8_lossy(&ref_seq).to_uppercase();
                    warn!(
                        "REF partially mismatched at {}:{} — {}/{} bases match ({:.1}%), \
                         correcting to FASTA REF",
                        chrom,
                        pos + 1,
                        matches,
                        max_len,
                        similarity * 100.0,
                    );
                    (
                        "PASS".to_string(),
                        "WARN_REF_CORRECTED".to_string(),
                        Some(fasta_ref),
                    )
                } else {
                    ("FAIL".to_string(), "REF_MISMATCH".to_string(), None)
                }
            }
        }
        Err(_) => ("FAIL".to_string(), "FETCH_FAILED".to_string(), None),
    }
}

/// Convert a MAF `-` allele to VCF style by fetching the anchor base.
///
/// - `ref = "-"` (insertion): MAF Start is the base before the insertion; that
///   base is the anchor, prepended to ALT (and it is the REF).
/// - `alt = "-"` (deletion): MAF Start is the first deleted base; the base
///   before it is the anchor, prepended to REF (and it is the ALT).
///
/// Called only for rows with a `-` allele; sequence alleles are used as
/// written. The prepared label is derived from the result's alleles.
///
/// # Returns
/// `(pos_0based, vcf_ref, vcf_alt)` or error if the FASTA fetch fails.
pub(crate) fn resolve_maf_anchor(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    start_pos: i64,
    ref_allele: &str,
    alt_allele: &str,
) -> anyhow::Result<(i64, String, String)> {
    let is_insertion = ref_allele == "-";
    // Insertion: anchor at Start (1-based) -> Start - 1 (0-based).
    // Deletion: anchor one base before Start -> Start - 2 (0-based).
    let anchor_pos_0based = if is_insertion { start_pos - 1 } else { start_pos - 2 };

    // Fetch anchor base, trying both chrom names
    let anchor_base = fetch_single_base(reader, chrom, anchor_pos_0based)?;
    let anchor = (anchor_base as char).to_uppercase().to_string();

    let (vcf_ref, vcf_alt) = if is_insertion {
        (anchor.clone(), format!("{anchor}{alt_allele}"))
    } else {
        (format!("{anchor}{ref_allele}"), anchor)
    };
    Ok((anchor_pos_0based, vcf_ref, vcf_alt))
}
