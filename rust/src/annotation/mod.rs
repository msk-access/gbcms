//! GTF-informed annotation index for RNA-seq splice masking and per-transcript counting.
//!
//! Provides [`AnnotationIndex`] — a thread-safe, read-only index of exon boundaries
//! built from a GTF file. Used by:
//!
//! - **P4a** (`counting/engine.rs`): Splice mask — suppress BAQ penalties near known
//!   exon boundaries, even at low coverage where consensus splicing fails.
//! - **P4b** (`counting/engine.rs`): Per-transcript counting — filter reads by
//!   splice-junction compatibility per transcript.
//! - **P4c** (`counting/engine.rs`): ASJD detection — compare junction usage between
//!   REF- and ALT-classified reads.
//!
//! # Architecture
//!
//! ```text
//! GTF file
//!   └─ gtf::parse_gtf()
//!       └─ AnnotationIndex
//!           ├── exon_trees: HashMap<u32, COITree>  (interval queries)
//!           ├── splice_sites: HashMap<u32, Vec<i32>>  (sorted boundary positions)
//!           ├── transcript_introns: HashMap<String, Vec<(i32, i32)>>
//!           └── chrom_map: HashMap<String, u32>
//! ```
//!
//! The index is built once in `count_bam_binned()` and shared via `Arc` across
//! rayon workers — the same pattern used for `editing_sites`.


mod gtf;

use std::collections::HashMap;

#[allow(unused_imports)] // IntervalTree needed for COITree::query trait
use coitrees::{COITree, IntervalTree};
use log::{debug, trace};

// Re-export the GTF parser for use by engine.rs (wired in P4a integration step)
#[allow(unused_imports)]
pub(crate) use gtf::parse_gtf;

// ─── Data Structures ─────────────────────────────────────────────────────────

/// Metadata for a single exon, stored in a flat Vec and referenced by COITree
/// node metadata (index into this Vec).
#[derive(Clone, Debug)]
pub struct ExonRecord {
    /// Ensembl/GENCODE transcript ID (e.g., "ENST00000269305").
    pub transcript_id: String,
    /// Ensembl/GENCODE gene ID (e.g., "ENSG00000141510").
    pub gene_id: String,
    /// Numeric chromosome ID (key into `AnnotationIndex::chrom_map`).
    pub chrom_id: u32,
    /// 0-based start position (inclusive).
    pub start: i32,
    /// 0-based end position (exclusive).
    pub end: i32,
    /// Strand: '+' or '-'.
    pub strand: char,
}

/// Intron structure for a single transcript, derived from sorted exons.
/// Used by P4b (per-transcript counting) and P4c (ASJD detection).
#[derive(Clone, Debug)]
pub struct TranscriptIntrons {
    /// Transcript ID.
    pub transcript_id: String,
    /// Numeric chromosome ID (key into chrom_map). Used by `is_junction_known()`
    /// to filter transcripts by chromosome before intron matching.
    pub chrom_id: u32,
    /// Introns as (start, end) pairs in genomic coordinates (0-based, exclusive end).
    /// Derived from consecutive exon boundaries: intron_start = exon_n.end,
    /// intron_end = exon_{n+1}.start.
    pub introns: Vec<(i32, i32)>,
}

/// Thread-safe, read-only index of exon boundaries from a GTF file.
///
/// Built once per `count_bam_binned()` call, shared via `Arc` across rayon workers.
/// Only loads chromosomes that have variants (variant-guided streaming) to reduce
/// memory footprint by 40-60% for typical panels.
pub struct AnnotationIndex {
    /// Chromosome → COITree for O(n+m) exon overlap queries.
    /// COITree metadata is an index into `self.exons`.
    exon_trees: HashMap<u32, COITree<usize, u32>>,

    /// Flat exon metadata, indexed by COITree node metadata.
    exons: Vec<ExonRecord>,

    /// Chromosome → sorted Vec of all exon start/end positions (deduplicated).
    /// Used for binary search in `nearest_splice_distance()`.
    splice_sites: HashMap<u32, Vec<i32>>,

    /// Transcript ID → intron boundaries.
    /// Used by P4b (per-transcript compatibility) and P4c (ASJD).
    transcript_introns: HashMap<String, TranscriptIntrons>,

    /// Chromosome name → numeric ID mapping (e.g., "1" → 0, "X" → 22).
    /// Normalized: no "chr" prefix.
    chrom_map: HashMap<String, u32>,
}

impl AnnotationIndex {
    /// Create a new AnnotationIndex from pre-parsed components.
    ///
    /// This is called by `gtf::parse_gtf()` — not directly by engine code.
    pub(crate) fn new(
        exon_trees: HashMap<u32, COITree<usize, u32>>,
        exons: Vec<ExonRecord>,
        splice_sites: HashMap<u32, Vec<i32>>,
        transcript_introns: HashMap<String, TranscriptIntrons>,
        chrom_map: HashMap<String, u32>,
    ) -> Self {
        let n_chroms = exon_trees.len();
        let n_exons = exons.len();
        let n_transcripts = transcript_introns.len();
        debug!(
            "AnnotationIndex built: {} chromosomes, {} exons, {} transcripts",
            n_chroms, n_exons, n_transcripts,
        );
        Self {
            exon_trees,
            exons,
            splice_sites,
            transcript_introns,
            chrom_map,
        }
    }

    // ─── P4a: Splice Mask ────────────────────────────────────────────────────

    /// Distance (bp) from `pos` to the nearest known exon boundary on `chrom`.
    ///
    /// Returns `i32::MAX` if:
    /// - The chromosome has no annotation (not in GTF, or filtered out by
    ///   variant-guided streaming).
    /// - The chromosome has no exon boundaries after dedup.
    ///
    /// Uses binary search on the pre-sorted `splice_sites` vec — O(log n).
    ///
    /// # Parameters
    ///
    /// - `chrom`: normalized chromosome name (no "chr" prefix).
    /// - `pos`: 0-based variant position.
    pub fn nearest_splice_distance(&self, chrom: &str, pos: i64) -> i32 {
        let chrom_id = match self.chrom_map.get(chrom) {
            Some(id) => *id,
            None => {
                trace!(
                    "nearest_splice_distance: chrom '{}' not in annotation index",
                    chrom
                );
                return i32::MAX;
            }
        };

        let sites = match self.splice_sites.get(&chrom_id) {
            Some(s) if !s.is_empty() => s,
            _ => return i32::MAX,
        };

        let pos_i32 = pos as i32;

        // Binary search for the insertion point
        match sites.binary_search(&pos_i32) {
            Ok(_) => 0, // Exact match — variant is AT an exon boundary
            Err(idx) => {
                // Check distance to neighbors on both sides
                let dist_left = if idx > 0 {
                    (pos_i32 - sites[idx - 1]).abs()
                } else {
                    i32::MAX
                };
                let dist_right = if idx < sites.len() {
                    (pos_i32 - sites[idx]).abs()
                } else {
                    i32::MAX
                };
                dist_left.min(dist_right)
            }
        }
    }

    // ─── P4b: Per-Transcript Counting ────────────────────────────────────────

    /// Get transcript IDs whose exons overlap the given position.
    ///
    /// Returns an empty Vec if the chromosome is not annotated or the position
    /// falls in an intergenic/intronic region.
    pub fn overlapping_transcripts(&self, chrom: &str, pos: i64) -> Vec<String> {
        let chrom_id = match self.chrom_map.get(chrom) {
            Some(id) => *id,
            None => return Vec::new(),
        };

        let tree = match self.exon_trees.get(&chrom_id) {
            Some(t) => t,
            None => return Vec::new(),
        };

        let pos_i32 = pos as i32;
        let mut transcript_ids = Vec::new();

        tree.query(pos_i32, pos_i32 + 1, |node| {
            // COITree metadata type varies by SIMD backend:
            //   nosimd: IntervalNode<usize, _> → field is usize
            //   NEON/AVX: Interval<&usize>     → field is &usize
            // Use Borrow<usize> to handle both uniformly.
            use std::borrow::Borrow;
            #[allow(noop_method_call)] // redundant on NEON/AVX, essential on nosimd
            let idx: &usize = node.metadata.borrow();
            let exon = &self.exons[*idx];
            // Deduplicate: a position may overlap multiple exons of the same transcript
            if !transcript_ids.contains(&exon.transcript_id) {
                transcript_ids.push(exon.transcript_id.clone());
            }
        });

        transcript_ids
    }

    /// Resolve the gene strand at a position from the exons overlapping it.
    ///
    /// Returns `Some('+')`/`Some('-')` when every stranded exon overlapping the
    /// position agrees. Returns `None` when the position is unannotated, overlaps
    /// only unstranded (`.`) exons, or overlaps exons on *both* strands (ambiguous,
    /// e.g. opposite-strand genes) — callers treat `None` as "do not enforce
    /// strandedness here" rather than guessing a direction.
    pub fn strand_at(&self, chrom: &str, pos: i64) -> Option<char> {
        let chrom_id = *self.chrom_map.get(chrom)?;
        let tree = self.exon_trees.get(&chrom_id)?;

        let pos_i32 = pos as i32;
        let mut seen_plus = false;
        let mut seen_minus = false;

        tree.query(pos_i32, pos_i32 + 1, |node| {
            // Metadata is an index into `self.exons` (see `overlapping_transcripts`).
            use std::borrow::Borrow;
            #[allow(noop_method_call)]
            let idx: &usize = node.metadata.borrow();
            match self.exons[*idx].strand {
                '+' => seen_plus = true,
                '-' => seen_minus = true,
                _ => {} // unstranded ('.') contributes no vote
            }
        });

        match (seen_plus, seen_minus) {
            (true, false) => Some('+'),
            (false, true) => Some('-'),
            _ => None, // unannotated, unstranded-only, or conflicting → no enforcement
        }
    }

    /// Get intron boundaries for a specific transcript.
    ///
    /// Returns `None` if the transcript ID is not in the index (should not
    /// happen if the ID came from `overlapping_transcripts()`).
    pub fn get_transcript_introns(&self, transcript_id: &str) -> Option<&TranscriptIntrons> {
        self.transcript_introns.get(transcript_id)
    }

    /// Check if a read's observed splice junctions are compatible with a
    /// transcript's annotated introns.
    ///
    /// A read is compatible if ALL its CIGAR N (RefSkip) operations match
    /// an annotated intron within ±`tolerance` bp. Reads without N ops are
    /// always compatible (ambiguous — no junction information).
    ///
    /// # Parameters
    ///
    /// - `observed_junctions`: splice junctions extracted from CIGAR N ops
    ///   via `rna::extract_splice_junctions()`.
    /// - `transcript`: the transcript to check compatibility against.
    /// - `tolerance`: maximum positional difference (bp) for a junction to
    ///   match an annotated intron. Default: 5 (matches STAR's junction
    ///   detection tolerance).
    pub fn is_read_compatible(
        &self,
        observed_junctions: &[(i64, i64)],
        transcript: &TranscriptIntrons,
        tolerance: i32,
    ) -> bool {
        if observed_junctions.is_empty() {
            return true; // No junctions = ambiguous = compatible
        }

        observed_junctions.iter().all(|(obs_start, obs_end)| {
            let os = *obs_start as i32;
            let oe = *obs_end as i32;
            transcript.introns.iter().any(|(t_start, t_end)| {
                (os - t_start).abs() <= tolerance && (oe - t_end).abs() <= tolerance
            })
        })
    }

    // ─── P4c: ASJD Helpers ───────────────────────────────────────────────────

    /// Check if an observed junction matches any annotated intron on the chromosome.
    ///
    /// Uses the pre-sorted `splice_sites` and `transcript_introns` to match within
    /// ±`tolerance` bp. Returns `true` if any annotated intron start/end pair
    /// matches the observed junction.
    pub fn is_junction_known(
        &self,
        chrom: &str,
        junction_start: i64,
        junction_end: i64,
        tolerance: i32,
    ) -> bool {
        let chrom_id = match self.chrom_map.get(chrom) {
            Some(id) => *id,
            None => return false,
        };

        // Only check transcripts on the target chromosome (chrom_id filter)
        for ti in self.transcript_introns.values().filter(|ti| ti.chrom_id == chrom_id) {
            for &(i_start, i_end) in &ti.introns {
                if (junction_start as i32 - i_start).abs() <= tolerance
                    && (junction_end as i32 - i_end).abs() <= tolerance
                {
                    return true;
                }
            }
        }

        false
    }

    // ─── Diagnostics ─────────────────────────────────────────────────────────

    /// Number of chromosomes loaded into the index.
    pub fn n_chromosomes(&self) -> usize {
        self.exon_trees.len()
    }

    /// Total number of exon records loaded.
    pub fn n_exons(&self) -> usize {
        self.exons.len()
    }

    /// Total number of transcripts with intron data.
    pub fn n_transcripts(&self) -> usize {
        self.transcript_introns.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use coitrees::{IntervalNode, IntervalTree};

    /// Build a minimal AnnotationIndex for testing.
    fn build_test_index() -> AnnotationIndex {
        // One chromosome "1" with two exons of one transcript:
        // Exon 1: [100, 200)
        // Exon 2: [300, 400)
        // → Intron: [200, 300)
        let exons = vec![
            ExonRecord {
                transcript_id: "ENST00000001".to_string(),
                gene_id: "ENSG00000001".to_string(),
                chrom_id: 0,
                start: 100,
                end: 200,
                strand: '+',
            },
            ExonRecord {
                transcript_id: "ENST00000001".to_string(),
                gene_id: "ENSG00000001".to_string(),
                chrom_id: 0,
                start: 300,
                end: 400,
                strand: '+',
            },
        ];

        let nodes: Vec<IntervalNode<usize, u32>> = exons
            .iter()
            .enumerate()
            .map(|(i, e)| IntervalNode::new(e.start, e.end, i))
            .collect();

        let tree = COITree::new(&nodes);
        let mut exon_trees = HashMap::new();
        exon_trees.insert(0u32, tree);

        let mut splice_sites = HashMap::new();
        splice_sites.insert(0u32, vec![100i32, 200, 300, 400]);

        let mut transcript_introns = HashMap::new();
        transcript_introns.insert(
            "ENST00000001".to_string(),
            TranscriptIntrons {
                transcript_id: "ENST00000001".to_string(),
                chrom_id: 0,
                introns: vec![(200, 300)],
            },
        );

        let mut chrom_map = HashMap::new();
        chrom_map.insert("1".to_string(), 0u32);

        AnnotationIndex::new(exon_trees, exons, splice_sites, transcript_introns, chrom_map)
    }

    // ── strand_at tests ──

    #[test]
    fn test_strand_at_resolves_and_disambiguates() {
        // chrom "1": a + exon [100,200), a - exon [500,600), an unstranded exon
        // [2000,2100), and an overlapping +/- pair [1000,1100)/[1050,1150) for the
        // ambiguous case.
        let exons = vec![
            ExonRecord { transcript_id: "tp".into(), gene_id: "gp".into(), chrom_id: 0, start: 100, end: 200, strand: '+' },
            ExonRecord { transcript_id: "tm".into(), gene_id: "gm".into(), chrom_id: 0, start: 500, end: 600, strand: '-' },
            ExonRecord { transcript_id: "ta".into(), gene_id: "ga".into(), chrom_id: 0, start: 1000, end: 1100, strand: '+' },
            ExonRecord { transcript_id: "tb".into(), gene_id: "gb".into(), chrom_id: 0, start: 1050, end: 1150, strand: '-' },
            ExonRecord { transcript_id: "tu".into(), gene_id: "gu".into(), chrom_id: 0, start: 2000, end: 2100, strand: '.' },
        ];
        let nodes: Vec<IntervalNode<usize, u32>> = exons
            .iter()
            .enumerate()
            .map(|(i, e)| IntervalNode::new(e.start, e.end, i))
            .collect();
        let mut exon_trees = HashMap::new();
        exon_trees.insert(0u32, COITree::new(&nodes));
        let mut chrom_map = HashMap::new();
        chrom_map.insert("1".to_string(), 0u32);
        let idx = AnnotationIndex::new(exon_trees, exons, HashMap::new(), HashMap::new(), chrom_map);

        assert_eq!(idx.strand_at("1", 150), Some('+'), "inside + exon");
        assert_eq!(idx.strand_at("1", 550), Some('-'), "inside - exon");
        assert_eq!(idx.strand_at("1", 1075), None, "overlapping opposite strands → ambiguous");
        assert_eq!(idx.strand_at("1", 2050), None, "unstranded exon → no enforcement");
        assert_eq!(idx.strand_at("1", 300), None, "intergenic gap → no annotation");
        assert_eq!(idx.strand_at("9", 150), None, "unknown chromosome");
    }

    // ── nearest_splice_distance tests ──

    #[test]
    fn test_at_exon_boundary() {
        let idx = build_test_index();
        assert_eq!(idx.nearest_splice_distance("1", 100), 0);
        assert_eq!(idx.nearest_splice_distance("1", 200), 0);
        assert_eq!(idx.nearest_splice_distance("1", 300), 0);
        assert_eq!(idx.nearest_splice_distance("1", 400), 0);
    }

    #[test]
    fn test_near_boundary() {
        let idx = build_test_index();
        assert_eq!(idx.nearest_splice_distance("1", 197), 3);
        assert_eq!(idx.nearest_splice_distance("1", 203), 3);
        assert_eq!(idx.nearest_splice_distance("1", 298), 2);
    }

    #[test]
    fn test_mid_exon() {
        let idx = build_test_index();
        assert_eq!(idx.nearest_splice_distance("1", 150), 50);
        assert_eq!(idx.nearest_splice_distance("1", 350), 50);
    }

    #[test]
    fn test_unknown_chrom() {
        let idx = build_test_index();
        assert_eq!(idx.nearest_splice_distance("X", 100), i32::MAX);
    }

    // ── overlapping_transcripts tests ──

    #[test]
    fn test_exonic_overlap() {
        let idx = build_test_index();
        let txs = idx.overlapping_transcripts("1", 150);
        assert_eq!(txs, vec!["ENST00000001"]);
    }

    #[test]
    fn test_intronic_no_overlap() {
        let idx = build_test_index();
        let txs = idx.overlapping_transcripts("1", 250);
        assert!(txs.is_empty());
    }

    #[test]
    fn test_intergenic_no_overlap() {
        let idx = build_test_index();
        let txs = idx.overlapping_transcripts("1", 50);
        assert!(txs.is_empty());
    }

    // ── is_read_compatible tests ──

    #[test]
    fn test_compatible_matching_junction() {
        let idx = build_test_index();
        let tx = idx.get_transcript_introns("ENST00000001").unwrap();
        // Junction [200, 300) matches the annotated intron exactly
        assert!(idx.is_read_compatible(&[(200, 300)], tx, 5));
    }

    #[test]
    fn test_compatible_within_tolerance() {
        let idx = build_test_index();
        let tx = idx.get_transcript_introns("ENST00000001").unwrap();
        // Junction [198, 302) is within ±5bp
        assert!(idx.is_read_compatible(&[(198, 302)], tx, 5));
    }

    #[test]
    fn test_incompatible_outside_tolerance() {
        let idx = build_test_index();
        let tx = idx.get_transcript_introns("ENST00000001").unwrap();
        // Junction [190, 310) is outside ±5bp
        assert!(!idx.is_read_compatible(&[(190, 310)], tx, 5));
    }

    #[test]
    fn test_no_junctions_always_compatible() {
        let idx = build_test_index();
        let tx = idx.get_transcript_introns("ENST00000001").unwrap();
        assert!(idx.is_read_compatible(&[], tx, 5));
    }
}
