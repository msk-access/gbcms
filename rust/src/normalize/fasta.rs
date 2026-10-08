//! FASTA I/O, REF validation, and MAF anchor resolution.

use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;

use bio::io::fasta;
use log::warn;

/// A FASTA reader that remembers the regions it has fetched: many reads ask for
/// the same few places (a spliced read's far exon), so each is read once. One per
/// worker thread (not `Sync`); the cache is cleared when it grows past
/// [`CachedFasta::MAX_REGIONS`].
pub(crate) struct CachedFasta {
    reader: RefCell<fasta::IndexedReader<File>>,
    cache: RefCell<FetchedRegions>,
}

/// Regions fetched: (contig, start, end) → the bases, or None when not held whole.
type FetchedRegions = HashMap<(String, i64, i64), Option<Vec<u8>>>;

impl CachedFasta {
    const MAX_REGIONS: usize = 100_000;

    /// Open an indexed FASTA; None when it cannot be opened.
    pub(crate) fn open(path: &str) -> Option<Self> {
        let reader = fasta::IndexedReader::from_file(&path).ok()?;
        Some(Self { reader: RefCell::new(reader), cache: RefCell::new(HashMap::new()) })
    }

    /// The reference bases `[start, end)` of `chrom` (any naming, as
    /// [`fetch_region`]), uppercased; None when the region is not held whole.
    pub(crate) fn bases(&self, chrom: &str, start: i64, end: i64) -> Option<Vec<u8>> {
        if start < 0 || end <= start {
            return None;
        }
        let key = (chrom.to_string(), start, end);
        if let Some(hit) = self.cache.borrow().get(&key) {
            return hit.clone();
        }
        let got = fetch_region(&mut self.reader.borrow_mut(), chrom, start as u64, end as u64)
            .ok()
            .filter(|b| b.len() == (end - start) as usize)
            .map(|b| b.iter().map(u8::to_ascii_uppercase).collect::<Vec<u8>>());
        let mut cache = self.cache.borrow_mut();
        if cache.len() >= Self::MAX_REGIONS {
            cache.clear();
        }
        cache.insert(key, got.clone());
        got
    }
}

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
/// allele, a MAF anchor) use [`fetch_region`] and fail instead. Upper-case: a
/// soft-masked FASTA's lower case is the same reference, and every reader of a
/// window (left-alignment, the shift region, the event reference, `ref_context`
/// and the aligners scoring against it) compares bases case-sensitively.
pub(crate) fn fetch_window(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    start: u64,
    end: u64,
) -> anyhow::Result<Vec<u8>> {
    // Inside the contig, as nearly every window is: the plain fetch, at no extra cost.
    if let Ok(mut seq) = fetch_region(reader, chrom, start, end) {
        seq.make_ascii_uppercase();
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
    buf.make_ascii_uppercase();
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

/// How far from the given position [`ref_offsets`] looks. Measured on 1.13M
/// sign-out rows: of the `REF_MISMATCH` rows whose REF (3+ bases) matches the
/// reference exactly within 10 bases, 91% do so within 3, and further out a
/// short REF matches by chance more often.
const REF_OFFSET_REACH: i64 = 3;
/// The shortest REF [`ref_offsets`] places: a 2-base REF matches the reference
/// within 3 bases by chance about one time in three.
const REF_OFFSET_MIN_LEN: usize = 3;

/// Offsets from `pos_0based`, within ±[`REF_OFFSET_REACH`], at which a REF
/// that failed validation matches the reference exactly (case-insensitive),
/// nearest first and the left one first at equal distance. Empty for a REF
/// shorter than [`REF_OFFSET_MIN_LEN`] (a MAF `-` included). It explains a
/// `REF_MISMATCH` row, typically a REF written at the wrong coordinate; the row
/// is neither moved nor counted (the input is taken as given).
pub(crate) fn ref_offsets(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    pos_0based: i64,
    ref_allele: &str,
) -> Vec<i64> {
    if ref_allele.len() < REF_OFFSET_MIN_LEN {
        return Vec::new();
    }
    let len = ref_allele.len() as i64;
    let mut offsets: Vec<i64> = (1..=REF_OFFSET_REACH).flat_map(|d| [-d, d]).collect();
    offsets.retain(|&o| {
        let start = pos_0based + o;
        start >= 0
            && fetch_region(reader, chrom, start as u64, (start + len) as u64)
                .is_ok_and(|seq| seq.eq_ignore_ascii_case(ref_allele.as_bytes()))
    });
    offsets
}

/// Convert a MAF `-` allele to VCF style by fetching the anchor base.
///
/// - `ref = "-"` (insertion): MAF Start is the base before the insertion; that
///   base is the anchor, prepended to ALT (and it is the REF).
/// - `alt = "-"` (deletion): MAF Start is the first deleted base; the base
///   before it is the anchor, prepended to REF (and it is the ALT). A deletion
///   at Start 1 has no base before it: the base after it is appended instead
///   (the VCF spec's form at position 1, as the VCF writer writes it), so it is
///   counted as the same event given as VCF is.
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
    if !is_insertion && start_pos == 1 {
        let after = fetch_single_base(reader, chrom, ref_allele.len() as i64)?;
        let after = (after as char).to_uppercase().to_string();
        return Ok((0, format!("{ref_allele}{after}"), after));
    }
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
        // A base before the contig start (0-based -1) is an error, never a u64
        // wrap. (A MAF deletion at Start 1 takes the base after it instead.)
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
