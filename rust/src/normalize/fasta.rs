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
    let mut buf = Vec::new();
    for name in &contig_names(chrom) {
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

/// Fetch a window of reference around a locus: the part of `[start, end)` that
/// lies on the contig. The FASTA reader rejects a window that passes the contig
/// end, and prep pads its windows on both sides, so an indel near the end lost
/// every window. Callers index from `start`, so only the end is clamped: a
/// shorter result means the window reached the contig end. Exact fetches (a REF
/// allele, a MAF anchor) use [`fetch_region`] and fail instead.
pub(crate) fn fetch_window(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    start: u64,
    end: u64,
) -> anyhow::Result<Vec<u8>> {
    // Inside the contig, as nearly every window is: the plain fetch, at no extra cost.
    if let Ok(seq) = fetch_region(reader, chrom, start, end) {
        return Ok(seq);
    }
    // Past the contig end: clamp to the length of the record fetch_region reads
    // (the first of the contig's names the index holds), and read that record,
    // so a FASTA holding two of its names cannot mix their lengths.
    let (name, len) = resolve_contig(reader, chrom)
        .ok_or_else(|| anyhow::anyhow!("FASTA has no contig {}", chrom))?;
    let end = end.min(len);
    if start >= end {
        anyhow::bail!("FASTA window {}:{}-{} lies past the contig end ({})", chrom, start, end, len);
    }
    let mut buf = Vec::new();
    reader.fetch(&name, start, end)?;
    reader.read(&mut buf)?;
    Ok(buf)
}

/// The FASTA record a contig is read from, with its length: the first of its
/// names, in the order [`fetch_region`] tries them, that the index holds.
fn resolve_contig(reader: &fasta::IndexedReader<File>, chrom: &str) -> Option<(String, u64)> {
    let records = reader.index.sequences();
    contig_names(chrom)
        .into_iter()
        .find_map(|n| records.iter().find(|r| r.name == n).map(|r| (n, r.len)))
}

/// The names a contig may go by in the FASTA: as given, with or without a `chr`
/// prefix, and each mitochondrial spelling.
fn contig_names(chrom: &str) -> Vec<String> {
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
    names
}

/// Fetch a single base from the reference, delegating to `fetch_region`.
pub(crate) fn fetch_single_base(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    pos_0based: i64,
) -> anyhow::Result<u8> {
    if pos_0based < 0 {
        // E.g. the anchor of a MAF deletion at Start 1: there is no base there.
        anyhow::bail!("no reference base before the start of {} (position {})", chrom, pos_0based);
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fetch_single_base_before_contig_start_is_an_error() {
        // A MAF deletion at Start 1 asks for the anchor at 0-based -1: an
        // error (FETCH_FAILED), never a u64 wrap.
        let dir = std::env::temp_dir();
        let fa = dir.join(format!("gbcms-fetch-neg-{}.fa", std::process::id()));
        std::fs::write(&fa, ">1\nACGT\n").unwrap();
        std::fs::write(fa.with_extension("fa.fai"), "1\t4\t3\t4\t5\n").unwrap();
        let mut reader = fasta::IndexedReader::from_file(&fa).unwrap();
        assert!(fetch_single_base(&mut reader, "1", -1).is_err());
        assert_eq!(fetch_single_base(&mut reader, "1", 0).unwrap(), b'A');
        let _ = std::fs::remove_file(fa.with_extension("fa.fai"));
        let _ = std::fs::remove_file(&fa);
    }

    #[test]
    fn test_a_window_is_clamped_to_the_contig_and_an_exact_fetch_is_not() {
        // chr-named FASTA, looked up as "1": the clamp finds the contig under
        // any name, as the fetch does.
        let dir = std::env::temp_dir();
        let fa = dir.join(format!("gbcms-fetch-window-{}.fa", std::process::id()));
        std::fs::write(&fa, ">chr1\nACGTACGTAC\n").unwrap();
        std::fs::write(fa.with_extension("fa.fai"), "chr1\t10\t6\t10\t11\n").unwrap();
        let mut reader = fasta::IndexedReader::from_file(&fa).unwrap();
        assert_eq!(fetch_window(&mut reader, "1", 6, 50).unwrap(), b"GTAC");
        assert_eq!(fetch_window(&mut reader, "1", 2, 6).unwrap(), b"GTAC");
        assert!(fetch_region(&mut reader, "1", 6, 50).is_err(), "an exact fetch stays exact");
        assert!(fetch_window(&mut reader, "1", 10, 20).is_err(), "nothing on the contig");
        let _ = std::fs::remove_file(fa.with_extension("fa.fai"));
        let _ = std::fs::remove_file(&fa);
    }

    #[test]
    fn test_a_window_is_clamped_to_the_record_it_reads() {
        // Two names of one contig at different lengths, the shorter listed first:
        // a window on "1" reads and clamps to "1", never to "chr1".
        let dir = std::env::temp_dir();
        let fa = dir.join(format!("gbcms-fetch-alias-{}.fa", std::process::id()));
        std::fs::write(&fa, ">chr1\nAAAA\n>1\nCCGGTTACGT\n").unwrap();
        std::fs::write(fa.with_extension("fa.fai"), "chr1\t4\t6\t4\t5\n1\t10\t14\t10\t11\n").unwrap();
        let mut reader = fasta::IndexedReader::from_file(&fa).unwrap();
        assert_eq!(fetch_window(&mut reader, "1", 6, 50).unwrap(), b"ACGT");
        assert_eq!(fetch_window(&mut reader, "chr1", 1, 50).unwrap(), b"AAA");
        let _ = std::fs::remove_file(fa.with_extension("fa.fai"));
        let _ = std::fs::remove_file(&fa);
    }
}
