//! GTF file parser for building [`AnnotationIndex`].
//!
//! Reads plain or gzip/BGZF-compressed GTF (detected by the gzip magic bytes, not
//! the extension) and loads the `exon` records of the variant chromosomes. Each line's
//! feature and chromosome are read from its bytes first, so the lines it skips cost a
//! scan to their third tab; the exon lines it keeps are then checked column by column
//! ([`super::gtf_line`], noodles-gtf's grammar). On top: variant-guided streaming, chromosome normalization, GENCODE
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
use super::{build_exon_trees, AnnotationIndex, ExonRecord, TranscriptIntrons, TranscriptKey};

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
    let started = std::time::Instant::now();
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

    // Exon coordinates per transcript, for intron derivation. A transcript is an ID on
    // one chromosome and strand: an ID the GTF reuses elsewhere is another transcript,
    // and its exons must not be sorted in with this one's (the consecutive pairs would
    // make introns that exist nowhere and miss the real ones).
    let mut transcript_exons: HashMap<TranscriptKey, Vec<(i32, i32)>> = HashMap::new();

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

        // Feature and chromosome from the raw bytes: only the exon lines of the
        // variant chromosomes are decoded and checked.
        let mut fields = bytes.splitn(4, |&b| b == b'\t');
        let (seqname, feature) = (fields.next().unwrap_or_default(), fields.nth(1));
        match feature {
            Some(b"exon") => {}
            Some(_) => {
                skipped_non_exon += 1;
                continue;
            }
            None => {
                unreadable += 1; // fewer than three columns: not a record
                continue;
            }
        }
        if seqname != last_seqname.as_bytes() {
            last_seqname = String::from_utf8_lossy(seqname).into_owned();
            // Canonicalize the chromosome so chr1/1 and chrM/MT reconcile across sources.
            last_chrom = crate::shared::contig::normalize_contig(&last_seqname);
            last_kept = std::str::from_utf8(seqname).is_ok() && variant_chroms.contains(&last_chrom);
        }
        // Variant-guided filter: skip chromosomes without variants
        if !last_kept {
            skipped_chrom += 1;
            continue;
        }

        let mut reject = |e: LineError| {
            rejected += 1;
            first_rejected.get_or_insert((total_lines, e));
            debug!("GTF line {} rejected: {}", total_lines, e);
        };
        let Ok(line) = std::str::from_utf8(bytes) else {
            reject(LineError::NotUtf8);
            continue;
        };
        let cols = match columns(line) {
            Ok(c) => c,
            Err(e) => {
                reject(e);
                continue;
            }
        };

        let exon = match cols.exon() {
            Ok(x) => x,
            Err(e) => {
                reject(e);
                continue;
            }
        };
        // Checked after the version strip, so "." and ".1" are as empty as "".
        let transcript_id = match exon.transcript_id.map(strip_gencode_version) {
            Some(id) if !id.is_empty() => id.to_string(),
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
            .entry(TranscriptKey { id: transcript_id, chrom_id, strand })
            .or_default()
            .push((start, end));
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
        if rejected > 0 {
            warn!(
                "GTF parser: every exon line on the variant chromosomes {:?} in {} was rejected \
                 as malformed ({} lines; names shown are normalized). RNA annotation will be \
                 inert.",
                variant_chroms, gtf_path, rejected,
            );
        } else if skipped_chrom > 0 {
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
                 {} unreadable) — wrong file or feature column? RNA annotation will be inert.",
                gtf_path, total_lines, skipped_non_exon, unreadable,
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

    warn_reused_transcript_ids(transcript_exons.keys(), gtf_path);

    // ── Derive introns per transcript ────────────────────────────────────────

    let mut transcript_introns: HashMap<TranscriptKey, TranscriptIntrons> = HashMap::new();

    for (key, mut exon_coords) in transcript_exons {
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
                key.clone(),
                TranscriptIntrons {
                    transcript_id: key.id,
                    chrom_id: key.chrom_id,
                    strand: key.strand,
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
    let index = AnnotationIndex::new(exon_trees, exons, splice_sites, transcript_introns, chrom_map);
    info!("GTF annotation loaded in {:.1}s", started.elapsed().as_secs_f64());
    Ok(index)
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

/// Warn once when the GTF uses a transcript ID at more than one chromosome or strand
/// (among the loaded chromosomes). Each place is read as its own transcript; the
/// warning names how many IDs are reused and the first of them, so a GTF that merges
/// distinct transcripts under one ID (version-stripped RefSeq IDs, say) is noticed.
fn warn_reused_transcript_ids<'a>(keys: impl Iterator<Item = &'a TranscriptKey>, gtf_path: &str) {
    let mut places: HashMap<&str, u32> = HashMap::new();
    for key in keys {
        *places.entry(key.id.as_str()).or_default() += 1;
    }
    let mut reused: Vec<&str> = places.iter().filter(|(_, &n)| n > 1).map(|(&id, _)| id).collect();
    if reused.is_empty() {
        return;
    }
    reused.sort_unstable();
    warn!(
        "GTF parser: {} transcript ID{} used on more than one chromosome or strand in {} \
         (first: {}); each chromosome and strand is read as its own transcript, with its own \
         introns.",
        reused.len(),
        if reused.len() == 1 { "" } else { "s" },
        gtf_path,
        reused[0],
    );
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
        let tx = idx.spliced_transcripts_at("1", 150)[0];
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

    /// Parse GTF text for the given variant chromosomes, through a temp file.
    fn parse_text(tag: &str, text: &str, chroms: &[&str]) -> AnnotationIndex {
        use std::io::Write;
        let path = std::env::temp_dir().join(format!("gbcms-gtf-{tag}-{}.gtf", std::process::id()));
        std::fs::File::create(&path).unwrap().write_all(text.as_bytes()).unwrap();
        let chroms: HashSet<String> = chroms.iter().map(|c| c.to_string()).collect();
        let idx = parse_gtf(path.to_str().unwrap(), &chroms).unwrap();
        std::fs::remove_file(&path).ok();
        idx
    }

    /// One exon line of transcript `id` (1-based closed `start`..`end`).
    fn exon_line(chrom: &str, start: i32, end: i32, strand: char, id: &str) -> String {
        format!(
            "{chrom}\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id \"G\"; transcript_id \"{id}\";\n"
        )
    }

    #[test]
    fn an_id_reused_on_another_chromosome_keeps_each_places_introns() {
        // TD: exons [100,200) [300,400) on chrom 1 (intron [200,300)) and [1000,1100)
        // [1200,1300) on chrom 2 (intron [1100,1200)). Joined by ID, the exons of both
        // sort together: an intron [400,1000) that exists nowhere, chrom 2's real intron
        // missing from chrom 2, and chrom 1's intron lost.
        let gtf = [
            exon_line("1", 101, 200, '+', "TD"),
            exon_line("1", 301, 400, '+', "TD"),
            exon_line("2", 1001, 1100, '+', "TD"),
            exon_line("2", 1201, 1300, '+', "TD"),
        ]
        .concat();
        let idx = parse_text("reuse-chrom", &gtf, &["1", "2"]);
        assert_eq!(idx.n_transcripts(), 2, "one transcript per place");
        assert!(idx.is_junction_known("1", 200, 300, 0));
        assert!(idx.is_junction_known("2", 1100, 1200, 0));
        assert!(!idx.is_junction_known("1", 1100, 1200, 0), "chrom 2's intron is not on chrom 1");
        assert!(!idx.is_junction_known("2", 200, 300, 0), "chrom 1's intron is not on chrom 2");
        for chrom in ["1", "2"] {
            assert!(!idx.is_junction_known(chrom, 400, 1000, 0), "the joined exons' intron, chrom {chrom}");
        }
        assert!(idx.intron_boundary_in_range("1", 200, 300, None));
        assert!(!idx.intron_boundary_in_range("1", 1100, 1200, None));
        assert!(idx.intron_boundary_in_range("2", 1100, 1200, None));
        assert!(!idx.intron_boundary_in_range("2", 200, 300, None));
    }

    #[test]
    fn an_id_reused_on_the_other_strand_of_a_chromosome_is_another_transcript() {
        let gtf = [
            exon_line("1", 101, 200, '+', "TD"),
            exon_line("1", 301, 400, '+', "TD"),
            exon_line("1", 1001, 1100, '-', "TD"),
            exon_line("1", 1201, 1300, '-', "TD"),
        ]
        .concat();
        let idx = parse_text("reuse-strand", &gtf, &["1"]);
        assert_eq!(idx.n_transcripts(), 2, "one transcript per strand");
        assert!(idx.is_junction_known("1", 200, 300, 0));
        assert!(idx.is_junction_known("1", 1100, 1200, 0));
        assert!(!idx.is_junction_known("1", 400, 1000, 0), "the joined exons' intron");
        // Each intron's boundaries belong to its own strand.
        assert!(idx.intron_boundary_in_range("1", 200, 200, Some('+')));
        assert!(!idx.intron_boundary_in_range("1", 200, 200, Some('-')));
        assert!(idx.intron_boundary_in_range("1", 1100, 1100, Some('-')));
        assert!(!idx.intron_boundary_in_range("1", 1100, 1100, Some('+')));
    }

    #[test]
    fn copies_of_a_transcript_at_one_place_stay_one_transcript() {
        // The same exons listed twice, and a version-stripped twin (T.1 and T.2 become T)
        // on the same chromosome and strand, are one transcript: exact duplicates drop out.
        let gtf = [
            exon_line("1", 101, 200, '+', "T.1"),
            exon_line("1", 301, 400, '+', "T.1"),
            exon_line("1", 101, 200, '+', "T.2"),
            exon_line("1", 301, 400, '+', "T.2"),
        ]
        .concat();
        let idx = parse_text("same-place", &gtf, &["1"]);
        assert_eq!(idx.n_transcripts(), 1);
        assert!(idx.is_junction_known("1", 200, 300, 0));
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
