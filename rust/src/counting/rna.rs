//! RNA-seq-specific alignment filters and utilities.
//!
//! Provides RNA-aware preprocessing for the allele counting engine:
//!
//! - [`is_valid_rna_alignment`] — STAR/HISAT2 MAPQ filter with NH:i:1 rescue
//! - [`is_sense_strand`] — dUTP strand-specific filter
//! - [`build_rna_editing_set`] — REDIportal A-to-I editing site loader
//! - [`apply_consensus_splicing`] — Consensus intron snipping from local CIGAR N ops
//! - [`has_splice_junction`] — Check if a read spans a splice junction (CIGAR N)

use std::collections::HashMap;
use std::collections::HashSet;
use std::io::BufRead;

use log::{debug, trace};
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;


/// Check if an RNA-seq alignment passes quality filters.
///
/// Standard MAPQ threshold applies first. If it fails, reads with
/// NH:i:1 (uniquely mapping) are rescued — STAR assigns MAPQ=255 for
/// unique mappers in most cases, but novel splice junctions and
/// multi-mapping loci get MAPQ=0-3 despite being unique.
///
/// ## NH Tag Rescue Logic
///
/// STAR and HISAT2 use the NH:i tag to report the number of reported
/// alignments. NH:i:1 means the read mapped to exactly one location.
/// These reads are biologically informative even with low MAPQ.
///
/// ## Parameters
///
/// - `record`: BAM record to evaluate
/// - `min_mapq`: minimum MAPQ threshold (typically 1 for RNA)
pub fn is_valid_rna_alignment(record: &Record, min_mapq: u8) -> bool {
    let mapq = record.mapq();

    // Fast path: passes MAPQ threshold
    if mapq >= min_mapq {
        return true;
    }

    // Rescue uniquely-mapped reads (NH == 1) that fall below the MAPQ floor:
    // STAR/HISAT2 assign low MAPQ to reads at novel splice junctions even when
    // they map to exactly one locus. The NH tag's integer width varies by writer
    // (U8/U16/U32/I8/I16/I32), so accept any integer encoding rather than only a
    // couple of widths — otherwise valid unique mappers are silently dropped.
    if nh_tag_value(record) == Some(1) {
        trace!("uniquely-mapped read (NH=1) rescued below MAPQ floor (MAPQ={})", mapq);
        return true;
    }

    false
}

/// Read the `NH` (number-of-reported-alignments) tag as an integer, regardless of
/// the width the aligner encoded it with. Returns `None` when the tag is absent or
/// is not an integer type.
fn nh_tag_value(record: &Record) -> Option<i64> {
    use rust_htslib::bam::record::Aux;
    match record.aux(b"NH") {
        Ok(Aux::U8(v)) => Some(i64::from(v)),
        Ok(Aux::U16(v)) => Some(i64::from(v)),
        Ok(Aux::U32(v)) => Some(i64::from(v)),
        Ok(Aux::I8(v)) => Some(i64::from(v)),
        Ok(Aux::I16(v)) => Some(i64::from(v)),
        Ok(Aux::I32(v)) => Some(i64::from(v)),
        _ => None,
    }
}


/// RNA library strand protocol — controls how a read's genomic orientation folds to
/// the transcript strand it originated from. Mirrors featureCounts `-s`:
///
/// | Variant      | featureCounts | Convention            | Sense mate |
/// |--------------|---------------|-----------------------|------------|
/// | `Reverse`    | `-s 2`        | dUTP / fr-firststrand  | R2        |
/// | `Forward`    | `-s 1`        | fr-secondstrand        | R1        |
/// | `Unstranded` | `-s 0`        | none                   | —         |
///
/// `Reverse` is the default — it is the FORTE RNA pipeline default and what every
/// prior gbcms RNA run assumed (the direction was previously hardcoded). For
/// `Unstranded` there is no transcript strand, so both strandedness enforcement and
/// ASJD strand-discordance are disabled.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Strandedness {
    /// dUTP / fr-firststrand (featureCounts `-s 2`): R2 is sense.
    Reverse,
    /// fr-secondstrand (featureCounts `-s 1`): R1 is sense.
    Forward,
    /// Unstranded (featureCounts `-s 0`): no transcript strand.
    Unstranded,
}

impl Strandedness {
    /// Parse the FFI protocol token. Returns `Err` for an unrecognized value rather
    /// than silently defaulting, so a typo (`--strandedness reverze`) fails loudly at
    /// the boundary instead of mis-orienting every downstream strand call.
    pub fn from_protocol(token: &str) -> Result<Self, String> {
        match token {
            "reverse" => Ok(Self::Reverse),
            "forward" => Ok(Self::Forward),
            "unstranded" => Ok(Self::Unstranded),
            other => Err(format!(
                "invalid strandedness '{}' (expected 'reverse', 'forward', or 'unstranded')",
                other
            )),
        }
    }
}

/// Check if a read is on the sense strand relative to the gene, under the library's
/// strand protocol.
///
/// The read's transcript strand (see [`read_transcript_strand`]) is compared against
/// the gene's annotated strand. Returns `true` (pass) whenever strand cannot be
/// enforced — the gene strand is unknown, or the library is `Unstranded` — so
/// enforcement degrades to a no-op rather than silently dropping every read.
///
/// ## Parameters
///
/// - `record`: BAM record to evaluate
/// - `gene_strand`: annotated gene strand ('+' or '-'); `None` ⇒ pass all
/// - `strandedness`: the library strand protocol
pub fn is_sense_strand(
    record: &Record,
    gene_strand: Option<char>,
    strandedness: Strandedness,
) -> bool {
    match (gene_strand, read_transcript_strand(record, strandedness)) {
        (None, _) => true,                 // gene strand unknown — cannot enforce
        (_, None) => true,                 // unstranded library — no transcript strand
        (Some(gs), Some(ts)) => ts == gs,
    }
}

/// Infer the transcript strand a read originated from, given the library protocol.
///
/// Returned in genomic terms (`'+'`/`'-'`) for direct comparison against a gene
/// strand, or `None` for an `Unstranded` library (no transcript strand exists). Both
/// mates of a normal FR pair resolve to the *same* transcript strand, so a genuine
/// junction is not spuriously split across strands. Shared by [`is_sense_strand`] and
/// the ASJD strand-discordance detector so both fold orientation identically.
///
/// Folding (mirrors featureCounts `-s`):
/// - `Reverse` (dUTP, `-s 2`): R2 = genomic strand (sense); R1/single-end = flipped.
/// - `Forward` (`-s 1`): R1/single-end = genomic strand (sense); R2 = flipped.
/// - `Unstranded` (`-s 0`): `None`.
pub fn read_transcript_strand(record: &Record, strandedness: Strandedness) -> Option<char> {
    let genomic = if record.is_reverse() { '-' } else { '+' };
    let flipped = if genomic == '+' { '-' } else { '+' };
    // is_last_in_template() is R2; first-in-pair and single-end both read as R1.
    let is_r2 = record.is_last_in_template();
    match strandedness {
        Strandedness::Unstranded => None,
        Strandedness::Reverse => Some(if is_r2 { genomic } else { flipped }),
        Strandedness::Forward => Some(if is_r2 { flipped } else { genomic }),
    }
}


/// Check if a read spans a splice junction (contains CIGAR N op).
///
/// Reads with RefSkip (N) operations in their CIGAR string span
/// introns — these are informative for validating splice-site variants
/// and are counted as `splice_spanning_count` in BaseCounts.
#[inline]
pub fn has_splice_junction(record: &Record) -> bool {
    record.cigar().iter().any(|op| matches!(op, Cigar::RefSkip(_)))
}

/// Extract splice junction coordinates from a read's CIGAR string.
///
/// Each CIGAR `N` (RefSkip) operation represents a spliced-out intron.
/// Returns a Vec of `(intron_start, intron_end)` pairs in 0-based genomic
/// coordinates (exclusive end), matching the convention used by
/// [`AnnotationIndex::is_read_compatible`].
///
/// Returns an empty Vec for reads without splice junctions (e.g., DNA reads,
/// or RNA reads that fall entirely within a single exon).
///
/// # Example
///
/// A read aligned at pos=100 with CIGAR `50M200N50M` has:
/// - Match: [100, 150)
/// - Intron: [150, 350)  ← returned as (150, 350)
/// - Match: [350, 400)
///
/// Used by:
/// - **Per-transcript counting**: splice-junction compatibility check per transcript.
/// - **ASJD**: comparing junction distributions between REF/ALT reads.
pub fn extract_splice_junctions(record: &Record) -> Vec<(i64, i64)> {
    let mut junctions = Vec::new();
    let mut ref_pos = record.pos(); // 0-based start position

    for op in record.cigar().iter() {
        match op {
            // Operations that consume reference bases
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) | Cigar::Del(len) => {
                ref_pos += *len as i64;
            }
            // RefSkip (N) = splice junction
            Cigar::RefSkip(len) => {
                let intron_start = ref_pos;
                let intron_end = ref_pos + *len as i64;
                junctions.push((intron_start, intron_end));
                ref_pos = intron_end;
            }
            // Operations that do NOT consume reference bases
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    junctions
}
/// Build an O(1) lookup set of known RNA editing sites from REDIportal TABLE1.
///
/// Parses the REDIportal TABLE1 format (tab-delimited, with header row):
///
/// | Col | Header    | Example      | Notes                |
/// |-----|-----------|--------------|----------------------|
/// | 0   | Accession | EDHSAAAA0000 | Unique ID (skip)     |
/// | 1   | Region    | chr1         | Chromosome           |
/// | 2   | Position  | 87158        | 1-based genomic pos  |
/// | 3   | Ref       | T            | Reference base       |
/// | 4   | Ed        | C            | Edited base          |
/// | 5   | Strand    | -            | Strand               |
///
/// Strips "chr" prefix for consistent matching with BAM contigs.
/// Converts 1-based REDIportal positions to 0-based BAM coordinates.
///
/// ## Supported formats
///
/// - Plain text: `TABLE1_hg38_v3.txt`
/// - Gzipped:    `TABLE1_hg38_v3.txt.gz`
///
/// Auto-detects gzip by `.gz` file extension.
///
/// ## Performance
///
/// Loads ~15.7M sites in ~3-5 seconds (plain) or ~5-8 seconds (gzipped).
/// Memory: ~80 bytes per site ≈ 1.2 GB.
pub fn build_rna_editing_set(db_path: &str) -> anyhow::Result<HashSet<(String, i64, u8, u8)>> {
    let file = std::fs::File::open(db_path)
        .map_err(|e| anyhow::anyhow!("Failed to open RNA editing DB '{}': {}", db_path, e))?;

    // Auto-detect gzip by file extension
    let reader: Box<dyn BufRead> = if db_path.ends_with(".gz") {
        debug!("RNA editing DB: detected gzip format for {}", db_path);
        Box::new(std::io::BufReader::new(
            flate2::read::GzDecoder::new(file),
        ))
    } else {
        Box::new(std::io::BufReader::new(file))
    };

    let mut sites = HashSet::new();
    let mut is_header = true;
    let mut skipped = 0u64;

    for line_result in reader.lines() {
        let line = line_result?;

        // Skip header row (first non-empty line)
        if is_header {
            is_header = false;
            trace!("RNA editing DB: skipping header: {}", &line[..line.len().min(80)]);
            continue;
        }
        if line.is_empty() {
            continue;
        }

        // Split enough columns to reach Ed (col 4): Accession, Region, Position, Ref, Ed.
        let fields: Vec<&str> = line.splitn(6, '\t').collect();
        if fields.len() < 5 {
            skipped += 1;
            continue;
        }

        // Col 1 = Region (chrom), 2 = Position (1-based), 3 = Ref, 4 = Ed.
        let chrom = crate::shared::contig::normalize_contig(fields[1]);
        let pos_1based = match fields[2].parse::<i64>() {
            Ok(p) => p,
            Err(_) => {
                skipped += 1;
                continue;
            }
        };
        // Ref/Ed must be single bases to verify a SNV's substitution; skip otherwise.
        let (ref_b, ed_b) = (fields[3].as_bytes(), fields[4].as_bytes());
        if ref_b.len() != 1 || ed_b.len() != 1 {
            skipped += 1;
            continue;
        }
        // Convert 1-based REDIportal → 0-based BAM coordinates; store the
        // genomic-forward edit so the flag can require the variant's substitution
        // (Ref>Ed) to match, not merely its position.
        sites.insert((
            chrom,
            pos_1based - 1,
            ref_b[0].to_ascii_uppercase(),
            ed_b[0].to_ascii_uppercase(),
        ));
    }

    if skipped > 0 {
        debug!("RNA editing DB: skipped {} malformed lines", skipped);
    }
    debug!("Loaded {} RNA editing sites from {}", sites.len(), db_path);
    Ok(sites)
}


/// Apply consensus intron snipping to a reference context.
///
/// Examines the CIGAR strings of local reads to find the most common
/// RefSkip (N) operations — these represent splice junctions. The
/// intronic bases are physically removed from the reference context
/// to produce a mature mRNA-like reference for alignment.
///
/// This is critical for RNA-seq allele classification: without intron
/// snipping, the reference haplotype contains intronic sequence that
/// doesn't exist in the mRNA, causing systematic misalignment.
///
/// ## Algorithm
///
/// 1. Collect all N (RefSkip) operations from local reads
/// 2. Map N ops to reference context coordinates
/// 3. Find consensus introns (present in >50% of reads with N ops)
/// 4. Remove consensus intronic bases from ref_ctx
///
/// ## Parameters
///
/// - `ref_ctx`: reference context bytes
/// - `local_reads`: reads overlapping the variant region
/// - `ref_start`: genomic start position of ref_ctx
///
/// ## Returns
///
/// A new `Vec<u8>` with intronic bases removed. If no consensus introns
/// are found, returns a copy of the original ref_ctx.
pub fn apply_consensus_splicing(
    ref_ctx: &[u8],
    local_reads: &[&Record],
    ref_start: i64,
) -> Vec<u8> {
    let ref_end = ref_start + ref_ctx.len() as i64;

    // Collect intron positions from reads (genomic coords)
    let mut intron_counts: HashMap<(i64, i64), usize> = HashMap::new();
    let mut reads_with_n = 0usize;

    for record in local_reads {
        let mut has_n = false;
        let mut rpos = record.pos();
        for op in record.cigar().iter() {
            match op {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) | Cigar::Del(len) => {
                    rpos += *len as i64;
                }
                Cigar::RefSkip(len) => {
                    let intron_start = rpos;
                    let intron_end = rpos + *len as i64;
                    has_n = true;

                    // Only count introns that overlap the ref_ctx window
                    if intron_start < ref_end && intron_end > ref_start {
                        *intron_counts.entry((intron_start, intron_end)).or_insert(0) += 1;
                    }
                    rpos += *len as i64;
                }
                Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }
        if has_n {
            reads_with_n += 1;
        }
    }

    if reads_with_n == 0 || intron_counts.is_empty() {
        return ref_ctx.to_vec();
    }

    // Find consensus introns (present in >50% of reads with N ops)
    let threshold = reads_with_n / 2;
    let mut consensus_introns: Vec<(i64, i64)> = intron_counts
        .into_iter()
        .filter(|&(_, count)| count > threshold)
        .map(|(coords, _)| coords)
        .collect();

    if consensus_introns.is_empty() {
        return ref_ctx.to_vec();
    }

    // Sort by start position (descending) so we can remove from right to left
    // without invalidating earlier indices
    consensus_introns.sort_by(|a, b| b.0.cmp(&a.0));

    let mut result = ref_ctx.to_vec();
    for (intron_start, intron_end) in &consensus_introns {
        // Map genomic coordinates to ref_ctx indices
        let ctx_start = (*intron_start - ref_start).max(0) as usize;
        let ctx_end = ((*intron_end - ref_start) as usize).min(result.len());

        if ctx_start < ctx_end && ctx_start < result.len() {
            trace!(
                "apply_consensus_splicing: removing intron [{}, {}) from context (ctx indices [{}, {}))",
                intron_start, intron_end, ctx_start, ctx_end,
            );
            result.drain(ctx_start..ctx_end);
        }
    }

    debug!(
        "apply_consensus_splicing: snipped {} introns, context {} → {} bases",
        consensus_introns.len(),
        ref_ctx.len(),
        result.len(),
    );

    result
}


#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;
    use std::ffi::CString;

    /// Build a synthetic BAM record for testing.
    fn build_record(seq: &[u8], quals: &[u8], cigar: CigarString, pos: i64) -> Record {
        let mut record = Record::new();
        let name = CString::new("test_read").unwrap();
        record.set(name.as_bytes(), Some(&cigar), seq, quals);
        record.set_pos(pos);
        record.set_mapq(255);
        record
    }

    // ── is_valid_rna_alignment tests ──

    #[test]
    fn test_rna_alignment_passes_mapq() {
        let seq = b"ACGTACGT";
        let quals = vec![30u8; 8];
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let mut record = build_record(seq, &quals, cigar, 100);
        record.set_mapq(20);
        assert!(is_valid_rna_alignment(&record, 1));
    }

    #[test]
    fn test_rna_alignment_fails_mapq_no_nh() {
        let seq = b"ACGTACGT";
        let quals = vec![30u8; 8];
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let mut record = build_record(seq, &quals, cigar, 100);
        record.set_mapq(0);
        // No NH tag → fails
        assert!(!is_valid_rna_alignment(&record, 1));
    }

    #[test]
    fn test_rna_alignment_nh1_rescue() {
        let seq = b"ACGTACGT";
        let quals = vec![30u8; 8];
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let mut record = build_record(seq, &quals, cigar, 100);
        record.set_mapq(0);
        // Add NH:i:1 tag
        record.push_aux(b"NH", rust_htslib::bam::record::Aux::U8(1)).unwrap();
        assert!(is_valid_rna_alignment(&record, 1));
    }

    #[test]
    fn test_rna_alignment_nh2_not_rescued() {
        let seq = b"ACGTACGT";
        let quals = vec![30u8; 8];
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let mut record = build_record(seq, &quals, cigar, 100);
        record.set_mapq(0);
        // NH:i:2 → multi-mapped, NOT rescued
        record.push_aux(b"NH", rust_htslib::bam::record::Aux::U8(2)).unwrap();
        assert!(!is_valid_rna_alignment(&record, 1));
    }

    #[test]
    fn test_rna_alignment_nh1_rescue_across_aux_widths() {
        use rust_htslib::bam::record::Aux;
        // Aligners encode NH with varying integer widths; NH==1 must rescue in all.
        for nh in [Aux::U16(1), Aux::U32(1), Aux::I8(1), Aux::I16(1), Aux::I32(1)] {
            let cigar = CigarString(vec![Cigar::Match(8)]);
            let mut record = build_record(b"ACGTACGT", &[30u8; 8], cigar, 100);
            record.set_mapq(0);
            record.push_aux(b"NH", nh).unwrap();
            assert!(is_valid_rna_alignment(&record, 1));
        }
        // A wider-typed multi-mapper (NH=2) still must not rescue.
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let mut record = build_record(b"ACGTACGT", &[30u8; 8], cigar, 100);
        record.set_mapq(0);
        record.push_aux(b"NH", Aux::U16(2)).unwrap();
        assert!(!is_valid_rna_alignment(&record, 1));
    }

    #[test]
    fn test_editing_set_carries_edit_bases() {
        use std::io::Write;
        let path = std::env::temp_dir().join("gbcms_lo11_redi.txt");
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "Accession\tRegion\tPosition\tRef\tEd\tStrand").unwrap();
            writeln!(f, "EDH0001\t7\t100\tA\tG\t+").unwrap(); // + strand A>G
            writeln!(f, "EDH0002\t12\t200\tT\tC\t-").unwrap(); // - strand T>C (forward genome)
        }
        let sites = build_rna_editing_set(path.to_str().unwrap()).unwrap();
        let _ = std::fs::remove_file(&path);

        // 1-based REDIportal → 0-based internal coordinates (pos - 1).
        assert!(sites.contains(&("7".to_string(), 99, b'A', b'G')), "A>G + site present");
        assert!(sites.contains(&("12".to_string(), 199, b'T', b'C')), "T>C - site present");
        // A different substitution at the same coordinate is NOT catalogued.
        assert!(!sites.contains(&("7".to_string(), 99, b'A', b'T')), "A>T must not match");
    }

    // ── is_sense_strand tests ──

    #[test]
    fn test_sense_strand_no_gene_strand_passes() {
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let record = build_record(b"ACGTACGT", &[30u8; 8], cigar, 100);
        assert!(is_sense_strand(&record, None, Strandedness::Reverse));
    }

    #[test]
    fn test_sense_strand_r2_forward_plus_gene() {
        // R2 (sense) + forward read → transcript strand '+' → matches gene '+'
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let mut record = build_record(b"ACGTACGT", &[30u8; 8], cigar, 100);
        // Set as R2
        record.set_last_in_template();
        // Forward read (not reversed) → '+'
        assert!(is_sense_strand(&record, Some('+'), Strandedness::Reverse));
    }

    #[test]
    fn test_sense_strand_r1_reverse_plus_gene() {
        // R1 (antisense) + reverse read → genomic '-' → transcript '+' → matches gene '+'
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let mut record = build_record(b"ACGTACGT", &[30u8; 8], cigar, 100);
        record.set_first_in_template();
        record.set_reverse();
        assert!(is_sense_strand(&record, Some('+'), Strandedness::Reverse));
    }

    #[test]
    fn test_sense_strand_r1_forward_fails_plus_gene() {
        // R1 (antisense) + forward read → genomic '+' → transcript '-' → does NOT match gene '+'
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let mut record = build_record(b"ACGTACGT", &[30u8; 8], cigar, 100);
        record.set_first_in_template();
        // Forward (not reverse) → transcript strand = '-' (antisense flip)
        assert!(!is_sense_strand(&record, Some('+'), Strandedness::Reverse));
    }

    // ── Strandedness: protocol folds ──

    #[test]
    fn test_strandedness_from_protocol() {
        assert_eq!(Strandedness::from_protocol("reverse"), Ok(Strandedness::Reverse));
        assert_eq!(Strandedness::from_protocol("forward"), Ok(Strandedness::Forward));
        assert_eq!(
            Strandedness::from_protocol("unstranded"),
            Ok(Strandedness::Unstranded)
        );
        // Unknown token is a loud error, never a silent default.
        assert!(Strandedness::from_protocol("reverze").is_err());
    }

    #[test]
    fn test_transcript_strand_reverse_protocol() {
        // dUTP / reverse: R2 = genomic strand (sense); R1 = flipped.
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let mut r2_fwd = build_record(b"ACGTACGT", &[30u8; 8], cigar.clone(), 100);
        r2_fwd.set_last_in_template();
        assert_eq!(read_transcript_strand(&r2_fwd, Strandedness::Reverse), Some('+'));

        let mut r2_rev = build_record(b"ACGTACGT", &[30u8; 8], cigar.clone(), 100);
        r2_rev.set_last_in_template();
        r2_rev.set_reverse();
        assert_eq!(read_transcript_strand(&r2_rev, Strandedness::Reverse), Some('-'));

        let mut r1_fwd = build_record(b"ACGTACGT", &[30u8; 8], cigar, 100);
        r1_fwd.set_first_in_template();
        assert_eq!(read_transcript_strand(&r1_fwd, Strandedness::Reverse), Some('-'));
    }

    #[test]
    fn test_transcript_strand_forward_inverts_reverse() {
        // Forward (fr-secondstrand) is the mirror of reverse: R1 = sense, R2 = flipped.
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let mut r1_fwd = build_record(b"ACGTACGT", &[30u8; 8], cigar.clone(), 100);
        r1_fwd.set_first_in_template();
        assert_eq!(read_transcript_strand(&r1_fwd, Strandedness::Forward), Some('+'));

        let mut r2_fwd = build_record(b"ACGTACGT", &[30u8; 8], cigar, 100);
        r2_fwd.set_last_in_template();
        assert_eq!(read_transcript_strand(&r2_fwd, Strandedness::Forward), Some('-'));

        // The same read folds to opposite transcript strands under the two protocols.
        assert_ne!(
            read_transcript_strand(&r1_fwd, Strandedness::Forward),
            read_transcript_strand(&r1_fwd, Strandedness::Reverse),
        );
    }

    #[test]
    fn test_transcript_strand_unstranded_is_none() {
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let mut r2 = build_record(b"ACGTACGT", &[30u8; 8], cigar, 100);
        r2.set_last_in_template();
        assert_eq!(read_transcript_strand(&r2, Strandedness::Unstranded), None);
    }

    #[test]
    fn test_is_sense_strand_forward_protocol() {
        // Under forward, an R1 forward read is sense ('+') → matches gene '+', not '-'.
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let mut r1_fwd = build_record(b"ACGTACGT", &[30u8; 8], cigar, 100);
        r1_fwd.set_first_in_template();
        assert!(is_sense_strand(&r1_fwd, Some('+'), Strandedness::Forward));
        assert!(!is_sense_strand(&r1_fwd, Some('-'), Strandedness::Forward));
    }

    #[test]
    fn test_is_sense_strand_unstranded_passes_all() {
        // Unstranded: no transcript strand to enforce against → always pass, even
        // against a definite gene strand (a no-op, not a silent drop of every read).
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let mut r1_fwd = build_record(b"ACGTACGT", &[30u8; 8], cigar, 100);
        r1_fwd.set_first_in_template();
        assert!(is_sense_strand(&r1_fwd, Some('+'), Strandedness::Unstranded));
        assert!(is_sense_strand(&r1_fwd, Some('-'), Strandedness::Unstranded));
    }

    // ── has_splice_junction tests ──

    #[test]
    fn test_has_splice_junction_with_n() {
        let cigar = CigarString(vec![
            Cigar::Match(50),
            Cigar::RefSkip(1000),
            Cigar::Match(50),
        ]);
        let record = build_record(&[b'A'; 100], &[30u8; 100], cigar, 100);
        assert!(has_splice_junction(&record));
    }

    #[test]
    fn test_has_splice_junction_without_n() {
        let cigar = CigarString(vec![Cigar::Match(100)]);
        let record = build_record(&[b'A'; 100], &[30u8; 100], cigar, 100);
        assert!(!has_splice_junction(&record));
    }

    // ── apply_consensus_splicing tests ──

    #[test]
    fn test_splicing_no_introns() {
        let cigar = CigarString(vec![Cigar::Match(10)]);
        let record = build_record(&[b'A'; 10], &[30u8; 10], cigar, 0);
        let reads: Vec<&Record> = vec![&record];
        let result = apply_consensus_splicing(b"AAAAAGGGGG", &reads, 0);
        assert_eq!(result, b"AAAAAGGGGG");
    }

    #[test]
    fn test_splicing_consensus_intron() {
        // 3 reads all have the same N op: skip bases 5-7 (genomic)
        let cigar = CigarString(vec![
            Cigar::Match(5),
            Cigar::RefSkip(3),
            Cigar::Match(2),
        ]);
        let r1 = build_record(&[b'A'; 7], &[30u8; 7], cigar.clone(), 0);
        let r2 = build_record(&[b'A'; 7], &[30u8; 7], cigar.clone(), 0);
        let r3 = build_record(&[b'A'; 7], &[30u8; 7], cigar, 0);
        let reads: Vec<&Record> = vec![&r1, &r2, &r3];

        // ref_ctx: 0123456789 (10 bases), intron at [5,8)
        let ref_ctx = b"AAAAAXYZGG";
        let result = apply_consensus_splicing(ref_ctx, &reads, 0);
        // Should remove bases at positions 5,6,7 (XYZ)
        assert_eq!(result, b"AAAAAGG");
    }
}
