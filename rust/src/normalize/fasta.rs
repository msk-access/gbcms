//! FASTA I/O, REF validation, and MAF anchor resolution.

use std::fs::File;

use bio::io::fasta;
use log::warn;

/// Fetch a region from FASTA, trying chrom variants (with/without chr prefix).
pub(crate) fn fetch_region(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    start: u64,
    end: u64,
) -> anyhow::Result<Vec<u8>> {
    let mut buf = Vec::new();

    // Try as-is
    if reader.fetch(chrom, start, end).is_ok() {
        buf.clear();
        reader.read(&mut buf)?;
        if !buf.is_empty() {
            return Ok(buf);
        }
    }

    // Try with chr prefix
    let chr_name = format!("chr{}", chrom);
    if reader.fetch(&chr_name, start, end).is_ok() {
        buf.clear();
        reader.read(&mut buf)?;
        if !buf.is_empty() {
            return Ok(buf);
        }
    }

    // Try stripping chr prefix
    if let Some(stripped) = chrom.strip_prefix("chr") {
        if reader.fetch(stripped, start, end).is_ok() {
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
/// Tuple of `(status, Option<fasta_ref>)`:
/// - `("PASS", None)` — exact match
/// - `("PASS_WARN_REF_CORRECTED", Some(fasta_ref))` — ≥90% match, corrected
/// - `("REF_MISMATCH", None)` — <90% match, rejected
/// - `("FETCH_FAILED", None)` — region could not be fetched
pub(crate) fn validate_ref(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    pos_0based: i64,
    ref_allele: &str,
) -> (String, Option<String>) {
    let ref_len = ref_allele.len() as u64;
    let pos = pos_0based as u64;

    let fetch_result = fetch_region(reader, chrom, pos, pos + ref_len);
    match fetch_result {
        Ok(ref_seq) => {
            // Fast path: exact match
            if ref_seq.eq_ignore_ascii_case(ref_allele.as_bytes()) {
                ("PASS".to_string(), None)
            } else {
                // Compute similarity: count matching bases (case-insensitive)
                let max_len = ref_seq.len().max(ref_allele.len());
                if max_len == 0 {
                    return ("FAIL_REF_MISMATCH".to_string(), None);
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
                    ("PASS;WARN_REF_CORRECTED".to_string(), Some(fasta_ref))
                } else {
                    ("FAIL_REF_MISMATCH".to_string(), None)
                }
            }
        }
        Err(_) => ("FAIL_FETCH_FAILED".to_string(), None),
    }
}

/// Convert a MAF-style indel to VCF-style by fetching the anchor base.
///
/// Handles three cases based on the `is_maf` flag:
/// - `ref = "-"` (insertion): Anchor is at `start_pos` (1-based), prepend to ALT
/// - `alt = "-"` (deletion): Anchor is at `start_pos - 1` (1-based), prepend to REF
/// - Complex (different-length, non-dash): Anchor at `start_pos - 1`, prepend to both
///
/// # Returns
/// `(pos_0based, vcf_ref, vcf_alt, variant_type)` or error if FASTA fetch fails.
pub(crate) fn resolve_maf_anchor(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    start_pos: i64,
    ref_allele: &str,
    alt_allele: &str,
) -> anyhow::Result<(i64, String, String, String)> {
    let is_insertion = ref_allele == "-";
    let is_deletion = alt_allele == "-";
    let is_snp = ref_allele.len() == 1
        && alt_allele.len() == 1
        && !is_insertion
        && !is_deletion;

    if is_snp {
        // SNPs: MAF Start_Position == VCF POS (1-based), convert to 0-based
        return Ok((
            start_pos - 1,
            ref_allele.to_string(),
            alt_allele.to_string(),
            "SNP".to_string(),
        ));
    }

    // Determine anchor position (0-based) and variant type
    let (anchor_pos_0based, vtype) = if is_insertion {
        // MAF insertion: Start_Position is the base BEFORE the insertion
        // Anchor is at start_pos (1-based) → start_pos - 1 (0-based)
        (start_pos - 1, "INSERTION")
    } else if is_deletion {
        // MAF deletion: Start_Position is the FIRST DELETED base
        // Anchor is one base before → start_pos - 2 (0-based)
        (start_pos - 2, "DELETION")
    } else {
        // Complex: different-length, non-dash alleles
        // Anchor is one base before start → start_pos - 2 (0-based)
        (start_pos - 2, "COMPLEX")
    };

    // Fetch anchor base, trying both chrom names
    let anchor_base = fetch_single_base(reader, chrom, anchor_pos_0based)?;
    let anchor_upper = (anchor_base as char).to_uppercase().to_string();

    // Build VCF-style alleles
    let (vcf_ref, vcf_alt) = if is_insertion {
        (anchor_upper.clone(), format!("{}{}", anchor_upper, alt_allele))
    } else if is_deletion {
        (format!("{}{}", anchor_upper, ref_allele), anchor_upper.clone())
    } else {
        // Complex
        (
            format!("{}{}", anchor_upper, ref_allele),
            format!("{}{}", anchor_upper, alt_allele),
        )
    };

    // VCF POS (0-based) = anchor position
    Ok((anchor_pos_0based, vcf_ref, vcf_alt, vtype.to_string()))
}
