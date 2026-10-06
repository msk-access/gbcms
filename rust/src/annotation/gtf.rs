//! GTF file parser for building [`AnnotationIndex`].
//!
//! Reads plain or gzip/BGZF-compressed GTF (detected by the gzip magic bytes, not
//! the extension) and loads the `exon` records of the variant chromosomes. Each line's
//! feature column is read first, so the other ~half of a GTF costs one tab scan;
//! exon lines are then checked column by column ([`super::gtf_line`], noodles-gtf's
//! grammar). On top: variant-guided streaming, chromosome normalization, GENCODE
//! version stripping, and intron derivation.
//!
//! # Supported GTF formats
//!
//! | Format | Chromosome style | Example attribute |
//! |--------|-----------------|-------------------|
//! | Ensembl v75+ | `1`, `X`, `MT` | `gene_id "ENSG00000141510"` |
//! | GENCODE v19+ | `chr1`, `chrX`, `chrM` | `gene_id "ENSG00000141510.11"` |
//!
//! Chromosome normalization: names are canonicalized via
//! `shared::contig::normalize_contig` (strips `chr`, folds `M`/`MT` aliases), so
//! `chr1`/`1` and `chrM`/`MT` reconcile across BAM, GTF, variants and editing DBs.

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};

use log::{debug, info, warn};

use super::gtf_line::{columns, LineError};
use super::{build_exon_trees, AnnotationIndex, ExonRecord, TranscriptIntrons};

/// Parse a GTF file and build an [`AnnotationIndex`].
///
/// Only loads exon records from chromosomes present in `variant_chroms`
/// (variant-guided streaming), typically reducing memory by 40-60% for targeted
/// panels. Malformed exon lines on those chromosomes are skipped, counted and
/// warned about once (with the first one's line number and reason).
///
/// # Parameters
///
/// - `gtf_path`: path to the GTF file, plain or gzip/BGZF-compressed.
/// - `variant_chroms`: set of normalized chromosome names (no "chr" prefix)
///   that have variants. Only these chromosomes are loaded.
///
/// # Errors
///
/// Returns `Err` if the file cannot be opened or read (including a corrupt gzip
/// stream).
pub fn parse_gtf(
    gtf_path: &str,
    variant_chroms: &HashSet<String>,
) -> anyhow::Result<AnnotationIndex> {
    info!("Loading GTF annotation from: {}", gtf_path);
    debug!(
        "Variant-guided filter: loading {} chromosomes: {:?}",
        variant_chroms.len(),
        variant_chroms,
    );

    let mut reader = open_gtf(gtf_path)?;
    let read_err = |e: std::io::Error| anyhow::anyhow!("Failed to read GTF file '{}': {}", gtf_path, e);

    let mut exons: Vec<ExonRecord> = Vec::new();
    let mut chrom_map: HashMap<String, u32> = HashMap::new();
    let mut next_chrom_id: u32 = 0;

    // Track exon coordinates per transcript for intron derivation
    let mut transcript_exons: HashMap<String, Vec<(i32, i32)>> = HashMap::new();
    // Track transcript → chrom_id for TranscriptIntrons construction
    let mut transcript_chrom_ids: HashMap<String, u32> = HashMap::new();

    let mut total_lines = 0u64;
    let mut skipped_non_exon = 0u64;
    let mut skipped_chrom = 0u64;
    let mut unreadable = 0u64;
    let mut rejected = 0u64;
    let mut first_rejected: Option<(u64, LineError)> = None;

    // The last seqname seen and its normalized form / filter verdict: GTFs group
    // a chromosome's lines, so this normalizes once per chromosome, not per line.
    let mut last_seqname = String::new();
    let mut last_chrom = String::new();
    let mut last_kept = false;

    let mut buf: Vec<u8> = Vec::with_capacity(4096);
    loop {
        buf.clear();
        if reader.read_until(b'\n', &mut buf).map_err(read_err)? == 0 {
            break;
        }
        total_lines += 1;
        let bytes = buf.strip_suffix(b"\n").unwrap_or(&buf);
        let bytes = bytes.strip_suffix(b"\r").unwrap_or(bytes);

        // Skip comment lines (GTF header)
        if bytes.is_empty() || bytes[0] == b'#' {
            continue;
        }
        let Ok(line) = std::str::from_utf8(bytes) else {
            unreadable += 1;
            debug!("GTF line {} is not UTF-8 text", total_lines);
            continue;
        };

        let mut keep_chrom = |seqname: &str| {
            if seqname != last_seqname {
                last_seqname.clear();
                last_seqname.push_str(seqname);
                // Canonicalize the chromosome so chr1/1 and chrM/MT reconcile across sources.
                last_chrom = crate::shared::contig::normalize_contig(seqname);
                last_kept = variant_chroms.contains(&last_chrom);
            }
            last_kept
        };
        let mut reject = |e: LineError| {
            rejected += 1;
            first_rejected.get_or_insert((total_lines, e));
            debug!("GTF line {} rejected: {}", total_lines, e);
        };

        let cols = match columns(line) {
            Ok(c) => c,
            Err(e) => {
                // Too few columns: a malformed exon line if its third column says
                // exon and it would have loaded; otherwise not a record at all.
                let mut it = line.split('\t');
                let (seqname, feature) = (it.next().unwrap_or(""), it.nth(1));
                if feature == Some("exon") && keep_chrom(seqname) {
                    reject(e);
                } else {
                    unreadable += 1;
                }
                continue;
            }
        };

        // Only parse exon records
        if cols.feature != "exon" {
            skipped_non_exon += 1;
            continue;
        }

        // Variant-guided filter: skip chromosomes without variants
        if !keep_chrom(cols.seqname) {
            skipped_chrom += 1;
            continue;
        }

        let exon = match cols.exon() {
            Ok(x) => x,
            Err(e) => {
                reject(e);
                continue;
            }
        };
        let transcript_id = match exon.transcript_id {
            Some(id) if !id.is_empty() => strip_gencode_version(id).to_string(),
            _ => {
                reject(LineError::MissingTranscriptId);
                continue;
            }
        };
        // Coordinates are 0-based half-open; the strand keeps unstranded ('.')
        // distinct from '+', so it propagates downstream as "no strand"
        // (gene_strand = None, no enforcement) rather than a false plus-strand call
        // that would mis-orient strandedness and splice-motif checks.
        let (start, end, strand) = (exon.start, exon.end, exon.strand);

        // Assign chromosome numeric ID
        let chrom_id = match chrom_map.get(&last_chrom) {
            Some(&id) => id,
            None => {
                let id = next_chrom_id;
                next_chrom_id += 1;
                chrom_map.insert(last_chrom.clone(), id);
                id
            }
        };

        // Store exon with chrom_id for tree construction
        exons.push(ExonRecord { transcript_id: transcript_id.clone(), chrom_id, start, end, strand });

        // Track for intron derivation
        transcript_exons
            .entry(transcript_id.clone())
            .or_default()
            .push((start, end));
        // Track chrom_id per transcript (all exons of a transcript share the same chromosome)
        transcript_chrom_ids.entry(transcript_id).or_insert(chrom_id);
    }

    if let Some((line_no, why)) = first_rejected {
        warn!(
            "GTF parser: {} exon line{} on the variant chromosomes rejected as malformed and \
             not loaded (first at line {}: {}) in {}",
            rejected,
            if rejected == 1 { "" } else { "s" },
            line_no,
            why,
            gtf_path,
        );
    }

    if exons.is_empty() {
        // Annotation is inert either way (splice distance, per-transcript counts and
        // strand all become no-ops), so make the *reason* loud and actionable. An
        // exon record that reached the chromosome filter either loaded, was rejected,
        // or bumped `skipped_chrom`; so `skipped_chrom > 0` means exons existed but
        // matched no variant chromosome, whereas `== 0` means the file had no
        // loadable `exon` records at all (likely the wrong file or feature column).
        // Common spellings (chr prefix, M/MT) are already normalized on both sides
        // before comparison, so a residual mismatch means the GTF genuinely lacks
        // those contigs or uses spellings the normalizer does not cover (e.g.
        // accession-style names).
        if skipped_chrom > 0 {
            warn!(
                "GTF parser: exon records exist but none on the variant chromosomes {:?} in {} \
                 ({} exon rows skipped by the chromosome filter; names shown are normalized). \
                 The GTF lacks these contigs, or uses spellings the normalizer does not cover \
                 (e.g. accession-style names). RNA annotation will be inert.",
                variant_chroms, gtf_path, skipped_chrom,
            );
        } else {
            warn!(
                "GTF parser: no 'exon' feature records found in {} ({} lines read, {} non-exon, \
                 {} unreadable, {} rejected) — wrong file or feature column? RNA annotation will \
                 be inert.",
                gtf_path, total_lines, skipped_non_exon, unreadable, rejected,
            );
        }
    }

    // Partial coverage is as silent a failure as an empty index, per chromosome:
    // a variant chromosome with zero loaded exons gets no chrom_map entry, so
    // splice distance, per-transcript counts, ASJD and strand resolution are all
    // inert for its variants while working normally elsewhere. Name the gaps.
    if !exons.is_empty() {
        warn_uncovered_variant_chroms(&chrom_map, variant_chroms);
    }

    info!(
        "GTF parser: {} lines read, {} exons loaded, {} skipped (non-exon: {}, chrom-filter: {}, \
         unreadable: {}, rejected: {})",
        total_lines,
        exons.len(),
        skipped_non_exon + skipped_chrom + unreadable + rejected,
        skipped_non_exon,
        skipped_chrom,
        unreadable,
        rejected,
    );

    let n_unstranded = exons.iter().filter(|e| e.strand == '.').count();
    if n_unstranded > 0 {
        debug!(
            "GTF parser: {} loaded exon records are unstranded ('.'); variants over \
             them will not have strandedness enforced or splice motifs oriented",
            n_unstranded,
        );
    }

    // ── Build sorted splice_sites per chromosome ─────────────────────────────

    let mut splice_sites: HashMap<u32, Vec<i32>> = HashMap::new();

    for exon in &exons {
        let sites = splice_sites.entry(exon.chrom_id).or_default();
        sites.push(exon.start);
        sites.push(exon.end);
    }

    // Sort and deduplicate each chromosome's splice sites
    for sites in splice_sites.values_mut() {
        sites.sort_unstable();
        sites.dedup();
    }

    debug!(
        "Splice sites: {} chromosomes, {} total boundary positions",
        splice_sites.len(),
        splice_sites.values().map(|v| v.len()).sum::<usize>(),
    );

    // ── Derive introns per transcript ────────────────────────────────────────

    let mut transcript_introns: HashMap<String, TranscriptIntrons> = HashMap::new();

    for (tx_id, mut exon_coords) in transcript_exons {
        // Sort exons by start position
        exon_coords.sort_by_key(|&(start, _)| start);
        // Deduplicate (some GTFs have duplicate exon entries)
        exon_coords.dedup();

        // Derive introns from consecutive exon pairs
        let mut introns: Vec<(i32, i32)> = Vec::new();
        for window in exon_coords.windows(2) {
            let intron_start = window[0].1; // end of previous exon
            let intron_end = window[1].0; // start of next exon
            if intron_end > intron_start {
                introns.push((intron_start, intron_end));
            }
        }

        if !introns.is_empty() {
            transcript_introns.insert(
                tx_id.clone(),
                TranscriptIntrons {
                    transcript_id: tx_id.clone(),
                    chrom_id: *transcript_chrom_ids.get(&tx_id).unwrap_or(&0),
                    introns,
                },
            );
        }
    }

    debug!(
        "Transcript introns: {} transcripts with intron data",
        transcript_introns.len(),
    );

    let exon_trees = build_exon_trees(&exons);
    Ok(AnnotationIndex::new(exon_trees, exons, splice_sites, transcript_introns, chrom_map))
}

/// Open a GTF for line reading, decompressing it when it starts with the gzip
/// magic bytes (plain gzip and BGZF alike: BGZF is a series of gzip members).
fn open_gtf(gtf_path: &str) -> anyhow::Result<Box<dyn BufRead>> {
    const BUF: usize = 1 << 20;
    let open_err = |e: std::io::Error| anyhow::anyhow!("Failed to open GTF file '{}': {}", gtf_path, e);
    let mut file = std::fs::File::open(gtf_path).map_err(open_err)?;
    let mut magic = [0u8; 2];
    let mut got = 0;
    while got < 2 {
        match file.read(&mut magic[got..]).map_err(open_err)? {
            0 => break,
            n => got += n,
        }
    }
    file.seek(SeekFrom::Start(0)).map_err(open_err)?;
    Ok(if got == 2 && magic == [0x1f, 0x8b] {
        debug!("GTF '{}' is gzip-compressed", gtf_path);
        let gz = flate2::read::MultiGzDecoder::new(BufReader::with_capacity(BUF, file));
        Box::new(BufReader::with_capacity(BUF, gz))
    } else {
        Box::new(BufReader::with_capacity(BUF, file))
    })
}

/// Warn once for every variant chromosome that has no loaded exons.
///
/// A variant chromosome absent from `chrom_map` makes splice distance,
/// per-transcript counts, ASJD and strand resolution silently inert for its
/// variants while annotation works normally elsewhere. Both sides of the
/// comparison are already contig-normalized (chr prefix stripped, M/MT folded),
/// so a gap means the GTF genuinely lacks the contig or spells it in a form the
/// normalizer does not cover (e.g. accession-style names).
fn warn_uncovered_variant_chroms(chrom_map: &HashMap<String, u32>, variant_chroms: &HashSet<String>) {
    let mut uncovered: Vec<&String> = variant_chroms
        .iter()
        .filter(|c| !chrom_map.contains_key(*c))
        .collect();
    if !uncovered.is_empty() {
        uncovered.sort();
        warn!(
            "GTF annotation: no exons loaded for variant chromosome(s) {:?} \
             (names shown are normalized) — splice distance, per-transcript \
             counts, ASJD and gene strand are inert for variants there. The \
             GTF lacks these contigs, or uses spellings the normalizer does \
             not cover (e.g. accession-style names).",
            uncovered,
        );
    }
}

// ─── GENCODE Version Stripping ───────────────────────────────────────────────

/// Strip GENCODE version suffix from an identifier.
///
/// GENCODE IDs have version suffixes (e.g., `ENST00000269305.8`).
/// Ensembl IDs do not. This function strips the suffix only if the
/// part after the last dot is purely numeric.
///
/// # Examples
///
/// ```text
/// "ENST00000269305.8" → "ENST00000269305"
/// "ENSG00000141510"   → "ENSG00000141510" (unchanged)
/// "gene.name.1"       → "gene.name"       (numeric after last dot)
/// ```
fn strip_gencode_version(id: &str) -> &str {
    if let Some(dot_pos) = id.rfind('.') {
        // Only strip if the part after the dot is numeric (version suffix)
        if id[dot_pos + 1..].chars().all(|c| c.is_ascii_digit()) {
            return &id[..dot_pos];
        }
    }
    id
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strip_gencode_version() {
        assert_eq!(strip_gencode_version("ENST00000269305.8"), "ENST00000269305");
        assert_eq!(strip_gencode_version("ENSG00000141510.11"), "ENSG00000141510");
        assert_eq!(strip_gencode_version("ENSG00000141510"), "ENSG00000141510");
        assert_eq!(strip_gencode_version("gene.name"), "gene.name"); // non-numeric after dot
    }

    #[test]
    fn test_parse_gtf_mini() {
        use std::io::Write;

        // Create a temporary GTF file
        let dir = std::env::temp_dir();
        let gtf_path = dir.join("test_mini.gtf");
        let mut f = std::fs::File::create(&gtf_path).unwrap();
        writeln!(f, "#!genome-build GRCh37").unwrap();
        writeln!(
            f,
            "1\tensembl\texon\t101\t200\t.\t+\t.\tgene_id \"ENSG001\"; transcript_id \"ENST001\";"
        )
        .unwrap();
        writeln!(
            f,
            "1\tensembl\texon\t301\t400\t.\t+\t.\tgene_id \"ENSG001\"; transcript_id \"ENST001\";"
        )
        .unwrap();
        writeln!(
            f,
            "2\tensembl\texon\t501\t600\t.\t-\t.\tgene_id \"ENSG002\"; transcript_id \"ENST002\";"
        )
        .unwrap();

        // Only load chrom "1"
        let mut variant_chroms = HashSet::new();
        variant_chroms.insert("1".to_string());

        let idx = parse_gtf(gtf_path.to_str().unwrap(), &variant_chroms).unwrap();

        // Should load 2 exons on chrom 1, skip chrom 2
        assert_eq!(idx.n_exons(), 2);
        assert_eq!(idx.n_chromosomes(), 1);

        // Should have 1 transcript with 1 intron [200, 300)
        assert_eq!(idx.n_transcripts(), 1);
        let tx = idx.get_transcript_introns("ENST001").unwrap();
        assert_eq!(tx.introns, vec![(200, 300)]);

        // Splice distance
        assert_eq!(idx.nearest_splice_distance("1", 100, 100), Some(0)); // exon start (0-based: 101-1=100)
        assert_eq!(idx.nearest_splice_distance("1", 200, 200), Some(0)); // exon end
        assert_eq!(idx.nearest_splice_distance("2", 550, 550), None); // filtered out

        // Cleanup
        std::fs::remove_file(&gtf_path).ok();
    }

    #[test]
    fn test_parse_gtf_gencode_chr_strip() {
        use std::io::Write;

        let dir = std::env::temp_dir();
        let gtf_path = dir.join("test_gencode.gtf");
        let mut f = std::fs::File::create(&gtf_path).unwrap();
        writeln!(
            f,
            "chr1\tGENCODE\texon\t101\t200\t.\t+\t.\tgene_id \"ENSG001.5\"; transcript_id \"ENST001.3\";"
        ).unwrap();

        let mut variant_chroms = HashSet::new();
        variant_chroms.insert("1".to_string()); // normalized, no "chr"

        let idx = parse_gtf(gtf_path.to_str().unwrap(), &variant_chroms).unwrap();
        assert_eq!(idx.n_exons(), 1);

        // Version suffix should be stripped
        let txs = idx.overlapping_transcripts("1", 150);
        assert_eq!(txs, vec!["ENST001"]);

        std::fs::remove_file(&gtf_path).ok();
    }

    #[test]
    fn test_parse_gtf_unstranded_strand() {
        // A GTF '.' strand must parse to an unstranded exon, so strand_at returns
        // None (no enforcement) rather than being coerced to '+'. The strand_at(None)
        // contract itself is covered in mod.rs; this guards the parser's strand
        // mapping end-to-end.
        use std::io::Write;
        let dir = std::env::temp_dir();
        let gtf_path = dir.join("test_unstranded.gtf");
        let mut f = std::fs::File::create(&gtf_path).unwrap();
        writeln!(f, "1\tensembl\texon\t101\t200\t.\t+\t.\tgene_id \"GP\"; transcript_id \"TP\";")
            .unwrap();
        writeln!(f, "1\tensembl\texon\t301\t400\t.\t.\t.\tgene_id \"GU\"; transcript_id \"TU\";")
            .unwrap();

        let mut variant_chroms = HashSet::new();
        variant_chroms.insert("1".to_string());
        let idx = parse_gtf(gtf_path.to_str().unwrap(), &variant_chroms).unwrap();

        assert_eq!(idx.n_exons(), 2);
        assert_eq!(idx.strand_at("1", 150), Some('+'), "stranded exon keeps '+'");
        assert_eq!(idx.strand_at("1", 350), None, "unstranded '.' exon → no enforcement");

        std::fs::remove_file(&gtf_path).ok();
    }

    #[test]
    fn test_parse_gtf_empty_index_no_panic() {
        // An empty index must be returned (Ok, not an error/panic) for both
        // distinguished causes — (a) exons exist but none on a variant chromosome
        // (contig-naming mismatch), and (b) the file carries no exon rows at all.
        use std::io::Write;
        let dir = std::env::temp_dir();

        // (a) naming mismatch: exon on "1", but the variant set asks for "7".
        let p1 = dir.join("test_empty_mismatch.gtf");
        let mut f1 = std::fs::File::create(&p1).unwrap();
        writeln!(f1, "1\tensembl\texon\t101\t200\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";")
            .unwrap();
        let mut chroms_a = HashSet::new();
        chroms_a.insert("7".to_string());
        let idx_a = parse_gtf(p1.to_str().unwrap(), &chroms_a).unwrap();
        assert_eq!(idx_a.n_exons(), 0, "no exon on the requested chromosome → empty");
        assert_eq!(idx_a.strand_at("7", 150), None);

        // (b) no exon rows at all (only a non-exon 'gene' feature).
        let p2 = dir.join("test_empty_noexon.gtf");
        let mut f2 = std::fs::File::create(&p2).unwrap();
        writeln!(f2, "1\tensembl\tgene\t101\t200\t.\t+\t.\tgene_id \"G\";").unwrap();
        let mut chroms_b = HashSet::new();
        chroms_b.insert("1".to_string());
        let idx_b = parse_gtf(p2.to_str().unwrap(), &chroms_b).unwrap();
        assert_eq!(idx_b.n_exons(), 0, "file has no exon rows → empty");

        std::fs::remove_file(&p1).ok();
        std::fs::remove_file(&p2).ok();
    }
}
