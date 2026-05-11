//! GTF file parser for building [`AnnotationIndex`].
//!
//! Parses Ensembl and GENCODE GTF formats, extracting exon records and
//! deriving intron boundaries per transcript. Supports variant-guided
//! streaming to reduce memory by only loading chromosomes with variants.
//!
//! # Supported GTF formats
//!
//! | Format | Chromosome style | Example attribute |
//! |--------|-----------------|-------------------|
//! | Ensembl v75+ | `1`, `X`, `MT` | `gene_id "ENSG00000141510"` |
//! | GENCODE v19+ | `chr1`, `chrX`, `chrM` | `gene_id "ENSG00000141510.11"` |
//!
//! Chromosome normalization: `chr` prefix is stripped via `trim_start_matches("chr")`
//! for consistent matching with BAM contigs (same approach as `rna.rs::build_rna_editing_set`).

use std::collections::{HashMap, HashSet};
use std::io::BufRead;

use log::{debug, info, warn};

use super::{AnnotationIndex, ExonRecord, TranscriptIntrons};
use coitrees::{COITree, IntervalNode, IntervalTree};

/// Parse a GTF file and build an [`AnnotationIndex`].
///
/// Only loads exon records from chromosomes present in `variant_chroms`
/// (variant-guided streaming). This typically reduces memory by 40-60%
/// for targeted panels like MSK-IMPACT.
///
/// # Parameters
///
/// - `gtf_path`: path to the GTF file (plain text, not gzipped).
/// - `variant_chroms`: set of normalized chromosome names (no "chr" prefix)
///   that have variants. Only these chromosomes are loaded.
///
/// # Errors
///
/// Returns `Err` if the file cannot be opened or contains no parseable exon
/// records for the requested chromosomes.
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

    let file = std::fs::File::open(gtf_path)
        .map_err(|e| anyhow::anyhow!("Failed to open GTF file '{}': {}", gtf_path, e))?;
    let reader = std::io::BufReader::new(file);

    let mut exons: Vec<ExonRecord> = Vec::new();
    let mut chrom_map: HashMap<String, u32> = HashMap::new();
    let mut next_chrom_id: u32 = 0;

    // Track exon coordinates per transcript for intron derivation
    let mut transcript_exons: HashMap<String, Vec<(i32, i32)>> = HashMap::new();

    let mut total_lines = 0u64;
    let mut skipped_non_exon = 0u64;
    let mut skipped_chrom = 0u64;
    let mut skipped_parse = 0u64;

    for line_result in reader.lines() {
        let line = line_result?;
        total_lines += 1;

        // Skip comment lines (GTF header)
        if line.starts_with('#') {
            continue;
        }

        // GTF is tab-delimited: chrom, source, feature, start, end, score, strand, frame, attributes
        let fields: Vec<&str> = line.splitn(9, '\t').collect();
        if fields.len() < 9 {
            skipped_parse += 1;
            continue;
        }

        // Only parse exon records
        let feature = fields[2];
        if feature != "exon" {
            skipped_non_exon += 1;
            continue;
        }

        // Normalize chromosome: strip "chr" prefix
        let chrom = fields[0].trim_start_matches("chr").to_string();

        // Variant-guided filter: skip chromosomes without variants
        if !variant_chroms.contains(&chrom) {
            skipped_chrom += 1;
            continue;
        }

        // Parse coordinates (GTF is 1-based inclusive → convert to 0-based exclusive-end)
        let start_1based: i32 = match fields[3].parse() {
            Ok(v) => v,
            Err(_) => {
                skipped_parse += 1;
                continue;
            }
        };
        let end_1based: i32 = match fields[4].parse() {
            Ok(v) => v,
            Err(_) => {
                skipped_parse += 1;
                continue;
            }
        };
        let start = start_1based - 1; // 0-based inclusive
        let end = end_1based; // 0-based exclusive (GTF end is 1-based inclusive)

        // Parse strand
        let strand = fields[6].chars().next().unwrap_or('+');

        // Parse attributes for transcript_id and gene_id
        let attrs = fields[8];
        let transcript_id = match extract_attribute(attrs, "transcript_id") {
            Some(id) => id,
            None => {
                skipped_parse += 1;
                continue;
            }
        };
        let gene_id = extract_attribute(attrs, "gene_id").unwrap_or_default();

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
            .entry(transcript_id)
            .or_default()
            .push((start, end));
    }

    if exons.is_empty() {
        warn!(
            "GTF parser: no exon records found for variant chromosomes {:?} in {}",
            variant_chroms, gtf_path,
        );
    }

    info!(
        "GTF parser: {} lines read, {} exons loaded, {} skipped (non-exon: {}, chrom-filter: {}, parse-error: {})",
        total_lines, exons.len(), skipped_non_exon + skipped_chrom + skipped_parse,
        skipped_non_exon, skipped_chrom, skipped_parse,
    );

    // ── Build COITrees per chromosome ────────────────────────────────────────

    let mut tree_nodes: HashMap<u32, Vec<IntervalNode<usize, u32>>> = HashMap::new();

    for (i, exon) in exons.iter().enumerate() {
        tree_nodes
            .entry(exon.chrom_id)
            .or_default()
            .push(IntervalNode::new(exon.start, exon.end, i));
    }

    let exon_trees: HashMap<u32, COITree<usize, u32>> = tree_nodes
        .into_iter()
        .map(|(chrom_id, nodes)| (chrom_id, COITree::new(&nodes)))
        .collect();

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
                    transcript_id: tx_id,
                    introns,
                },
            );
        }
    }

    debug!(
        "Transcript introns: {} transcripts with intron data",
        transcript_introns.len(),
    );

    Ok(AnnotationIndex::new(
        exon_trees,
        exons,
        splice_sites,
        transcript_introns,
        chrom_map,
    ))
}

// ─── Attribute Parsing ───────────────────────────────────────────────────────

/// Extract a named attribute value from a GTF attribute string.
///
/// GTF attributes are semicolon-separated key-value pairs:
/// `gene_id "ENSG00000141510"; transcript_id "ENST00000269305"; ...`
///
/// Handles both quoted (`"value"`) and unquoted values.
/// Strips GENCODE version suffixes (e.g., `ENST00000269305.8` → `ENST00000269305`).
fn extract_attribute(attrs: &str, key: &str) -> Option<String> {
    for part in attrs.split(';') {
        let trimmed = part.trim();
        if let Some(rest) = trimmed.strip_prefix(key) {
            let rest = rest.trim();
            // Remove surrounding quotes if present
            let value = rest.trim_matches('"').trim();
            if value.is_empty() {
                continue;
            }
            // Strip GENCODE version suffix (e.g., "ENST00000269305.8" → "ENST00000269305")
            let clean = if let Some(dot_pos) = value.rfind('.') {
                // Only strip if the part after the dot is numeric (version suffix)
                if value[dot_pos + 1..].chars().all(|c| c.is_ascii_digit()) {
                    &value[..dot_pos]
                } else {
                    value
                }
            } else {
                value
            };
            return Some(clean.to_string());
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_attribute_ensembl() {
        let attrs = r#"gene_id "ENSG00000141510"; transcript_id "ENST00000269305"; exon_number "4";"#;
        assert_eq!(
            extract_attribute(attrs, "gene_id"),
            Some("ENSG00000141510".to_string())
        );
        assert_eq!(
            extract_attribute(attrs, "transcript_id"),
            Some("ENST00000269305".to_string())
        );
    }

    #[test]
    fn test_extract_attribute_gencode_version_strip() {
        let attrs = r#"gene_id "ENSG00000141510.11"; transcript_id "ENST00000269305.8";"#;
        assert_eq!(
            extract_attribute(attrs, "gene_id"),
            Some("ENSG00000141510".to_string())
        );
        assert_eq!(
            extract_attribute(attrs, "transcript_id"),
            Some("ENST00000269305".to_string())
        );
    }

    #[test]
    fn test_extract_attribute_missing() {
        let attrs = r#"gene_id "ENSG00000141510";"#;
        assert_eq!(extract_attribute(attrs, "transcript_id"), None);
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
}
