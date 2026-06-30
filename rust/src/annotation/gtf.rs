//! GTF file parser for building [`AnnotationIndex`].
//!
//! Uses the [`noodles_gtf`] crate for standards-compliant record parsing,
//! with our business logic layered on top: variant-guided streaming,
//! chromosome normalization, GENCODE version stripping, and intron derivation.
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
//! for consistent matching with BAM contigs (same approach as `rna.rs::build_rna_editing_set`).

use std::collections::{HashMap, HashSet};
use std::io::BufRead;

use log::{debug, info, warn};
use noodles_gtf as gtf;

use super::cache::{GtfIndexBundle, CACHE_FORMAT_VERSION};
use super::{AnnotationIndex, ExonRecord, TranscriptIntrons};

/// Parse a GTF file and build an [`AnnotationIndex`].
///
/// Uses `noodles_gtf::Reader` for standards-compliant GTF parsing. Only loads
/// exon records from chromosomes present in `variant_chroms` (variant-guided
/// streaming), typically reducing memory by 40-60% for targeted panels.
///
/// # Parameters
///
/// - `gtf_path`: path to the GTF file (plain text, not gzipped).
/// - `variant_chroms`: set of normalized chromosome names (no "chr" prefix)
///   that have variants. Only these chromosomes are loaded.
///
/// # Errors
///
/// Returns `Err` if the file cannot be opened or contains unparseable records.
pub fn parse_gtf(
    gtf_path: &str,
    variant_chroms: &HashSet<String>,
) -> anyhow::Result<AnnotationIndex> {
    Ok(parse_gtf_to_bundle(gtf_path, variant_chroms)?.into_index())
}

/// Parse a GTF into the serializable [`GtfIndexBundle`] — the M5a cache payload:
/// everything an [`AnnotationIndex`] needs *except* the COITrees, which are rebuilt
/// from the exon records by [`GtfIndexBundle::into_index`]. Splitting the parse out
/// here lets the cache layer persist/restore the bundle without touching the
/// arch-specific trees. This is the function that does the ~8.7s text parse.
pub(crate) fn parse_gtf_to_bundle(
    gtf_path: &str,
    variant_chroms: &HashSet<String>,
) -> anyhow::Result<GtfIndexBundle> {
    info!("Loading GTF annotation from: {}", gtf_path);
    debug!(
        "Variant-guided filter: loading {} chromosomes: {:?}",
        variant_chroms.len(),
        variant_chroms,
    );

    let file = std::fs::File::open(gtf_path)
        .map_err(|e| anyhow::anyhow!("Failed to open GTF file '{}': {}", gtf_path, e))?;
    let buf_reader = std::io::BufReader::new(file);

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
    let mut skipped_parse = 0u64;

    // Use noodles-gtf line-by-line parsing with our filtering logic.
    // We read raw lines and parse via Record::from_str to handle comments
    // and maintain variant-guided streaming (skip chroms without variants).
    for line_result in buf_reader.lines() {
        let line = line_result?;
        total_lines += 1;

        // Skip comment lines (GTF header)
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        // Parse the line using noodles-gtf
        let record: gtf::Record = match line.parse() {
            Ok(r) => r,
            Err(e) => {
                skipped_parse += 1;
                debug!("GTF parse error at line {}: {}", total_lines, e);
                continue;
            }
        };

        // Only parse exon records
        if record.ty() != "exon" {
            skipped_non_exon += 1;
            continue;
        }

        // Canonicalize the chromosome so chr1/1 and chrM/MT reconcile across sources.
        let chrom = crate::shared::contig::normalize_contig(record.reference_sequence_name());

        // Variant-guided filter: skip chromosomes without variants
        if !variant_chroms.contains(&chrom) {
            skipped_chrom += 1;
            continue;
        }

        // Convert coordinates: noodles Position is 1-based → 0-based exclusive-end
        let start_1based: i32 = usize::from(record.start()) as i32;
        let end_1based: i32 = usize::from(record.end()) as i32;
        let start = start_1based - 1; // 0-based inclusive
        let end = end_1based; // 0-based exclusive (GTF end is 1-based inclusive)

        // Parse strand via noodles. Keep unstranded ('.') distinct from '+' so it
        // propagates downstream as "no strand" (gene_strand = None, no enforcement)
        // rather than a false plus-strand call that would mis-orient strandedness
        // and splice-motif checks.
        let strand = match record.strand() {
            Some(noodles_gtf::record::Strand::Forward) => '+',
            Some(noodles_gtf::record::Strand::Reverse) => '-',
            _ => '.', // unstranded / unknown
        };

        // Extract transcript_id and gene_id from noodles-parsed attributes
        let attrs = record.attributes();
        let transcript_id = match attrs.get("transcript_id") {
            Some(id) => strip_gencode_version(id).to_string(),
            None => {
                skipped_parse += 1;
                continue;
            }
        };
        let gene_id = attrs
            .get("gene_id")
            .map(|id| strip_gencode_version(id).to_string())
            .unwrap_or_default();

        // Assign chromosome numeric ID
        let chrom_id = *chrom_map.entry(chrom).or_insert_with(|| {
            let id = next_chrom_id;
            next_chrom_id += 1;
            id
        });

        // Store exon with chrom_id for tree construction
        exons.push(ExonRecord {
            transcript_id: transcript_id.clone(),
            gene_id,
            chrom_id,
            start,
            end,
            strand,
        });

        // Track for intron derivation
        transcript_exons
            .entry(transcript_id.clone())
            .or_default()
            .push((start, end));
        // Track chrom_id per transcript (all exons of a transcript share the same chromosome)
        transcript_chrom_ids.entry(transcript_id).or_insert(chrom_id);
    }

    if exons.is_empty() {
        // Annotation is inert either way (splice distance, per-transcript counts and
        // strand all become no-ops), so make the *reason* loud and actionable. An
        // exon record that reached the chromosome filter either loaded or bumped
        // `skipped_chrom`; so `skipped_chrom > 0` means exons existed but matched no
        // variant chromosome (a naming mismatch), whereas `== 0` means the file had
        // no `exon` feature records at all (likely the wrong file or feature column).
        if skipped_chrom > 0 {
            warn!(
                "GTF parser: exon records exist but none on the variant chromosomes {:?} in {} \
                 ({} exon rows skipped by the chromosome filter) — likely a contig-naming \
                 mismatch (e.g. chr1 vs 1, chrM vs MT). RNA annotation will be inert.",
                variant_chroms, gtf_path, skipped_chrom,
            );
        } else {
            warn!(
                "GTF parser: no 'exon' feature records found in {} ({} lines read, {} non-exon, \
                 {} parse errors) — wrong file or feature column? RNA annotation will be inert.",
                gtf_path, total_lines, skipped_non_exon, skipped_parse,
            );
        }
    }

    info!(
        "GTF parser: {} lines read, {} exons loaded, {} skipped (non-exon: {}, chrom-filter: {}, parse-error: {})",
        total_lines, exons.len(), skipped_non_exon + skipped_chrom + skipped_parse,
        skipped_non_exon, skipped_chrom, skipped_parse,
    );

    let n_unstranded = exons.iter().filter(|e| e.strand == '.').count();
    if n_unstranded > 0 {
        debug!(
            "GTF parser: {} loaded exon records are unstranded ('.'); variants over \
             them will not have strandedness enforced or splice motifs oriented",
            n_unstranded,
        );
    }

    // (COITrees are not built here — they are rebuilt from `exons` by
    // GtfIndexBundle::into_index, so the cache stores only the portable intermediate.)

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

    Ok(GtfIndexBundle {
        format_version: CACHE_FORMAT_VERSION,
        exons,
        splice_sites,
        transcript_introns,
        chrom_map,
    })
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
        assert_eq!(idx.nearest_splice_distance("1", 100), 0); // exon start (0-based: 101-1=100)
        assert_eq!(idx.nearest_splice_distance("1", 200), 0); // exon end
        assert_eq!(idx.nearest_splice_distance("2", 550), i32::MAX); // filtered out

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
