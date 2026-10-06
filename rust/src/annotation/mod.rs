//! GTF-informed annotation index for RNA-seq splice masking and per-transcript counting.
//!
//! Provides [`AnnotationIndex`] — a thread-safe, read-only index of exon boundaries
//! built from a GTF file. Used by:
//!
//! - **Splice mask** (`counting/engine.rs`): suppress BAQ penalties near known
//!   exon boundaries — the GTF gives authoritative junction positions
//!   independent of read coverage.
//! - **Per-transcript counting** (`counting/engine.rs`): filter reads by
//!   splice-junction compatibility per transcript.
//! - **ASJD detection** (`counting/engine.rs`): compare junction usage between
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
//!           ├── intron_boundaries: HashMap<u32, Vec<(i32, char)>>  (derived: donor/acceptor sites + strand)
//!           ├── transcript_trees: HashMap<u32, COITree>  (derived: transcript spans, for strand)
//!           └── chrom_map: HashMap<String, u32>
//! ```
//!
//! The index is built once in `count_bam_binned()` and shared via `Arc` across
//! rayon workers — the same pattern used for `editing_sites`.


mod gtf;
mod gtf_line;

use std::collections::HashMap;

#[allow(unused_imports)] // IntervalTree needed for COITree::query trait
use coitrees::{COITree, IntervalNode, IntervalTree};
use log::{debug, trace};

pub(crate) use gtf::parse_gtf;

// ─── Data Structures ─────────────────────────────────────────────────────────

/// Metadata for a single exon, stored in a flat Vec and referenced by COITree
/// node metadata (index into this Vec).
#[derive(Clone, Debug)]
pub struct ExonRecord {
    /// Ensembl/GENCODE transcript ID (e.g., "ENST00000269305").
    pub transcript_id: String,
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
/// Used by per-transcript counting and ASJD detection.
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
    /// Used by per-transcript compatibility and ASJD.
    transcript_introns: HashMap<String, TranscriptIntrons>,

    /// Chromosome → sorted (position, strand) of every annotated intron
    /// boundary: true donor/acceptor sites only (no transcript termini).
    /// Derived in `new()` from `transcript_introns` + exon strands.
    intron_boundaries: HashMap<u32, Vec<(i32, char)>>,

    /// Chromosome → COITree of every transcript's span (first exon start to last
    /// exon end, introns included); metadata is the transcript's strand
    /// (`transcript_strands`). Derived in `new()` from `exons`.
    transcript_trees: HashMap<u32, COITree<usize, u32>>,
    transcript_strands: Vec<char>,

    /// Chromosome name → numeric ID mapping (e.g., "1" → 0, "X" → 22).
    /// Normalized: no "chr" prefix.
    chrom_map: HashMap<String, u32>,
}

/// Build the per-chromosome exon interval trees from the flat exon list; the
/// COITree metadata is the index into `exons`.
pub(crate) fn build_exon_trees(exons: &[ExonRecord]) -> HashMap<u32, COITree<usize, u32>> {
    let mut tree_nodes: HashMap<u32, Vec<IntervalNode<usize, u32>>> = HashMap::new();
    for (i, exon) in exons.iter().enumerate() {
        // COITree intervals are end-inclusive; exons are half-open.
        tree_nodes
            .entry(exon.chrom_id)
            .or_default()
            .push(IntervalNode::new(exon.start, exon.end - 1, i));
    }
    tree_nodes
        .into_iter()
        .map(|(chrom_id, nodes)| (chrom_id, COITree::new(&nodes)))
        .collect()
}

/// Per-chromosome sorted, deduplicated (position, strand) list of annotated
/// intron boundaries (each intron's start and exclusive end), with each
/// transcript's strand taken from its exons ('.' when unstranded or unknown).
fn derive_intron_boundaries(
    exons: &[ExonRecord],
    transcript_introns: &HashMap<String, TranscriptIntrons>,
) -> HashMap<u32, Vec<(i32, char)>> {
    let strand_of: HashMap<&str, char> =
        exons.iter().map(|e| (e.transcript_id.as_str(), e.strand)).collect();
    let mut out: HashMap<u32, Vec<(i32, char)>> = HashMap::new();
    for ti in transcript_introns.values() {
        let strand = strand_of.get(ti.transcript_id.as_str()).copied().unwrap_or('.');
        let v = out.entry(ti.chrom_id).or_default();
        for &(a, b) in &ti.introns {
            v.push((a, strand));
            v.push((b, strand));
        }
    }
    for v in out.values_mut() {
        v.sort_unstable();
        v.dedup();
    }
    out
}

/// Per-chromosome trees of transcript spans (first exon start to last exon end,
/// end-inclusive for COITree), with each span's strand in the returned list.
/// A span is one transcript's exons on one chromosome and strand: an ID the GTF
/// reuses on another chromosome or strand (PAR copies on X and Y, alternate
/// loci, version-stripped RefSeq IDs) gives a span per place, not one across
/// them. Copies on the same chromosome and strand still share one span.
fn derive_transcript_trees(exons: &[ExonRecord]) -> (HashMap<u32, COITree<usize, u32>>, Vec<char>) {
    let mut spans: HashMap<(&str, u32, char), (i32, i32)> = HashMap::new();
    for e in exons {
        let s = spans.entry((e.transcript_id.as_str(), e.chrom_id, e.strand)).or_insert((e.start, e.end));
        s.0 = s.0.min(e.start);
        s.1 = s.1.max(e.end);
    }
    let mut strands = Vec::with_capacity(spans.len());
    let mut nodes: HashMap<u32, Vec<IntervalNode<usize, u32>>> = HashMap::new();
    let mut sorted: Vec<_> = spans.into_iter().collect();
    sorted.sort_unstable(); // stable indices across runs
    for ((_, chrom, strand), (start, end)) in sorted {
        nodes.entry(chrom).or_default().push(IntervalNode::new(start, end - 1, strands.len()));
        strands.push(strand);
    }
    (nodes.into_iter().map(|(c, n)| (c, COITree::new(&n))).collect(), strands)
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
        let intron_boundaries = derive_intron_boundaries(&exons, &transcript_introns);
        let (transcript_trees, transcript_strands) = derive_transcript_trees(&exons);
        debug!(
            "AnnotationIndex built: {} chromosomes, {} exons, {} transcripts, {} intron boundaries",
            n_chroms,
            n_exons,
            n_transcripts,
            intron_boundaries.values().map(|v| v.len()).sum::<usize>(),
        );
        Self {
            exon_trees,
            exons,
            splice_sites,
            transcript_introns,
            intron_boundaries,
            transcript_trees,
            transcript_strands,
            chrom_map,
        }
    }

    // ─── Splice Mask ─────────────────────────────────────────────────────────

    /// Distance (bp, unsigned) from the span `[first, last]` to the nearest
    /// known exon boundary on `chrom`, exonic and intronic alike: the least
    /// distance from any base of the span, 0 when a boundary lies inside it
    /// (Ensembl VEP's overlap view). A one-base span is the distance from that
    /// base. A boundary is an exon's first base or its exclusive end (the first
    /// intron base), so an exon's last base is 1bp from its right edge.
    ///
    /// None when the contig has no annotation (not in the GTF, filtered out by
    /// variant-guided streaming, or no exon boundaries after dedup) — never a
    /// sentinel distance.
    ///
    /// Uses binary search on the pre-sorted `splice_sites` vec — O(log n).
    ///
    /// # Parameters
    ///
    /// - `chrom`: contig name in any naming; normalized here like every other
    ///   annotation lookup (`chr1` ~ `1`, `chrM` ~ `M` ~ `MT`).
    /// - `first`, `last`: 0-based, inclusive (a variant's REF bases).
    pub fn nearest_splice_distance(&self, chrom: &str, first: i64, last: i64) -> Option<i32> {
        let key = crate::shared::contig::normalize_contig(chrom);
        let chrom_id = match self.chrom_map.get(&key) {
            Some(id) => *id,
            None => {
                trace!(
                    "nearest_splice_distance: chrom '{}' not in annotation index",
                    chrom
                );
                return None;
            }
        };

        let sites = match self.splice_sites.get(&chrom_id) {
            Some(s) if !s.is_empty() => s,
            _ => return None,
        };

        let first = first as i32;
        let last = (last as i32).max(first);
        let idx = sites.partition_point(|&s| s < first);
        // The first boundary at or after `first` (0 when inside the span), and
        // the last one before it. `sites` is non-empty, so one of them exists.
        let right = sites.get(idx).map(|&s| (s - last).max(0));
        let left = idx.checked_sub(1).map(|j| first - sites[j]);
        right.into_iter().chain(left).min()
    }

    /// Every annotated intron boundary (an intron's first base or its exclusive
    /// end, either strand) in `[lo, hi]`, ascending.
    pub fn intron_boundaries_in(&self, chrom: &str, lo: i64, hi: i64) -> Vec<i64> {
        let Some(sites) = self.chrom_map.get(chrom).and_then(|id| self.intron_boundaries.get(id)) else {
            return Vec::new();
        };
        let lo = lo.clamp(i32::MIN as i64, i32::MAX as i64) as i32;
        let hi = hi.clamp(i32::MIN as i64, i32::MAX as i64) as i32;
        let idx = sites.partition_point(|&(p, _)| p < lo);
        let mut out: Vec<i64> = sites[idx..].iter().take_while(|&&(p, _)| p <= hi).map(|&(p, _)| p as i64).collect();
        out.dedup();
        out
    }

    /// Whether an annotated intron boundary (a true donor/acceptor site —
    /// transcript termini excluded) lies in `[lo, hi]` (inclusive, 0-based:
    /// an intron's first base or its exclusive end). With `strand` given,
    /// only boundaries of transcripts on that strand (or unstranded ones)
    /// count, so an antisense gene's splice sites cannot stand in for the
    /// variant's own. Binary search to the range, then a short scan.
    pub fn intron_boundary_in_range(
        &self,
        chrom: &str,
        lo: i64,
        hi: i64,
        strand: Option<char>,
    ) -> bool {
        let sites = match self.chrom_map.get(chrom).and_then(|id| self.intron_boundaries.get(id)) {
            Some(s) => s,
            None => return false,
        };
        let lo = lo.clamp(i32::MIN as i64, i32::MAX as i64) as i32;
        let hi = hi.clamp(i32::MIN as i64, i32::MAX as i64) as i32;
        let idx = sites.partition_point(|&(p, _)| p < lo);
        sites[idx..]
            .iter()
            .take_while(|&&(p, _)| p <= hi)
            .any(|&(_, s)| strand.is_none_or(|want| s == want || s == '.'))
    }

    // ─── Per-Transcript Counting ─────────────────────────────────────────────

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

        tree.query(pos_i32, pos_i32, |node| {
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

    /// Resolve the gene strand at a position: from the exons overlapping it, or, at
    /// a position no stranded exon covers (intronic, including splice sites), from
    /// the transcripts spanning it (first exon start to last exon end).
    ///
    /// Returns `Some('+')`/`Some('-')` when every stranded exon (else every
    /// stranded spanning transcript) agrees. Returns `None` when the position is
    /// unannotated, covered only by unstranded (`.`) features, or covered on *both*
    /// strands (overlapping opposite-strand genes, whose transcripts both carry the
    /// allele) — callers treat `None` as "do not enforce strandedness here" rather
    /// than guessing a direction.
    pub fn strand_at(&self, chrom: &str, pos: i64) -> Option<char> {
        let chrom_id = *self.chrom_map.get(chrom)?;
        let tree = self.exon_trees.get(&chrom_id)?;

        let pos_i32 = pos as i32;
        let mut seen_plus = false;
        let mut seen_minus = false;

        tree.query(pos_i32, pos_i32, |node| {
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

        if !seen_plus && !seen_minus {
            if let Some(spans) = self.transcript_trees.get(&chrom_id) {
                spans.query(pos_i32, pos_i32, |node| {
                    use std::borrow::Borrow;
                    #[allow(noop_method_call)]
                    let idx: &usize = node.metadata.borrow();
                    match self.transcript_strands[*idx] {
                        '+' => seen_plus = true,
                        '-' => seen_minus = true,
                        _ => {}
                    }
                });
            }
        }
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

    // ─── ASJD Helpers ────────────────────────────────────────────────────────

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

    /// Build a minimal AnnotationIndex for testing.
    fn build_test_index() -> AnnotationIndex {
        // One chromosome "1" with two exons of one transcript:
        // Exon 1: [100, 200)
        // Exon 2: [300, 400)
        // → Intron: [200, 300)
        let exons = vec![
            ExonRecord {
                transcript_id: "ENST00000001".to_string(),
                chrom_id: 0,
                start: 100,
                end: 200,
                strand: '+',
            },
            ExonRecord {
                transcript_id: "ENST00000001".to_string(),
                chrom_id: 0,
                start: 300,
                end: 400,
                strand: '+',
            },
        ];

        let exon_trees = build_exon_trees(&exons);

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
            ExonRecord { transcript_id: "tp".into(), chrom_id: 0, start: 100, end: 200, strand: '+' },
            ExonRecord { transcript_id: "tm".into(), chrom_id: 0, start: 500, end: 600, strand: '-' },
            ExonRecord { transcript_id: "ta".into(), chrom_id: 0, start: 1000, end: 1100, strand: '+' },
            ExonRecord { transcript_id: "tb".into(), chrom_id: 0, start: 1050, end: 1150, strand: '-' },
            ExonRecord { transcript_id: "tu".into(), chrom_id: 0, start: 2000, end: 2100, strand: '.' },
        ];
        let exon_trees = build_exon_trees(&exons);
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

    #[test]
    fn strand_at_an_intronic_position_comes_from_the_spanning_transcripts() {
        // tp '+' exons [100,200) [300,400); tm '-' exons [600,700) [900,1000); tb '-'
        // exons [250,260) [340,350) spans part of tp's intron and exon.
        let ex = |t: &str, s: i32, e: i32, st: char| ExonRecord {
            transcript_id: t.into(), chrom_id: 0, start: s, end: e, strand: st,
        };
        let exons = vec![
            ex("tp", 100, 200, '+'), ex("tp", 300, 400, '+'),
            ex("tm", 600, 700, '-'), ex("tm", 900, 1000, '-'),
            ex("tb", 250, 260, '-'), ex("tb", 340, 350, '-'),
        ];
        let exon_trees = build_exon_trees(&exons);
        let chrom_map = HashMap::from([("1".to_string(), 0u32)]);
        let idx = AnnotationIndex::new(exon_trees, exons, HashMap::new(), HashMap::new(), chrom_map);
        assert_eq!(idx.strand_at("1", 200), Some('+'), "donor +1: tp spans it");
        assert_eq!(idx.strand_at("1", 201), Some('+'), "donor +2");
        assert_eq!(idx.strand_at("1", 230), Some('+'), "deep intronic, before tb");
        assert_eq!(idx.strand_at("1", 270), None, "tp and tb both span it");
        assert_eq!(idx.strand_at("1", 320), Some('+'), "tp's exon wins over tb's span");
        assert_eq!(idx.strand_at("1", 345), None, "exons of both strands");
        assert_eq!(idx.strand_at("1", 800), Some('-'), "tm's intron");
        assert_eq!(idx.strand_at("1", 500), None, "between genes");
        assert_eq!(idx.strand_at("1", 1000), None, "one past tm's last exon");
    }

    #[test]
    fn a_transcript_id_reused_on_another_chromosome_spans_neither_gap() {
        // TD '-' has an exon on chrom 0 at [1100,1150) and, reused (a PAR copy or a
        // version-stripped ID), one on chrom 1 at [10,30): no span on chrom 0 may
        // reach back to 10.
        let ex = |c: u32, s: i32, e: i32| ExonRecord {
            transcript_id: "TD".into(), chrom_id: c, start: s, end: e, strand: '-',
        };
        let exons = vec![ex(0, 1100, 1150), ex(1, 10, 30)];
        let exon_trees = build_exon_trees(&exons);
        let chrom_map = HashMap::from([("1".to_string(), 0u32), ("2".to_string(), 1u32)]);
        let idx = AnnotationIndex::new(exon_trees, exons, HashMap::new(), HashMap::new(), chrom_map);
        assert_eq!(idx.strand_at("1", 51), None, "intergenic on chrom 1");
        assert_eq!(idx.strand_at("1", 1120), Some('-'));
        assert_eq!(idx.strand_at("2", 20), Some('-'));
    }

    #[test]
    fn an_exon_holds_its_edge_bases_and_not_the_bases_beside_it() {
        // build_test_index: one '+' transcript, exons [100,200) [300,400).
        let idx = build_test_index();
        for pos in [100, 199, 300, 399] {
            assert_eq!(idx.overlapping_transcripts("1", pos), vec!["ENST00000001"], "{pos}");
        }
        for pos in [99, 200, 299, 400] {
            assert!(idx.overlapping_transcripts("1", pos).is_empty(), "{pos}");
        }
    }

    // ── nearest_splice_distance tests ──

    #[test]
    fn test_intron_boundary_in_range() {
        // build_test_index: exons [100,200) [300,400) on '+' → intron [200,300).
        // splice_sites also holds termini 100 and 400; intron boundaries must not.
        let idx = build_test_index();
        assert!(idx.intron_boundary_in_range("1", 200, 200, None), "donor boundary, degenerate range");
        assert!(idx.intron_boundary_in_range("1", 298, 301, None), "acceptor boundary inside range");
        assert!(idx.intron_boundary_in_range("1", 150, 200, Some('+')), "same strand, inclusive upper end");
        assert!(!idx.intron_boundary_in_range("1", 150, 200, Some('-')), "antisense query ignores '+' introns");
        assert!(!idx.intron_boundary_in_range("1", 95, 105, None), "transcript start is not a splice site");
        assert!(!idx.intron_boundary_in_range("1", 395, 405, None), "transcript end is not a splice site");
        assert!(!idx.intron_boundary_in_range("1", 201, 299, None), "inside the intron, no boundary");
        assert!(!idx.intron_boundary_in_range("2", 0, 1000, None), "unannotated chromosome");
    }

    #[test]
    fn test_at_exon_boundary() {
        let idx = build_test_index();
        assert_eq!(idx.nearest_splice_distance("1", 100, 100), Some(0));
        assert_eq!(idx.nearest_splice_distance("1", 200, 200), Some(0));
        assert_eq!(idx.nearest_splice_distance("1", 300, 300), Some(0));
        assert_eq!(idx.nearest_splice_distance("1", 400, 400), Some(0));
    }

    #[test]
    fn test_near_boundary() {
        let idx = build_test_index();
        assert_eq!(idx.nearest_splice_distance("1", 197, 197), Some(3));
        assert_eq!(idx.nearest_splice_distance("1", 203, 203), Some(3));
        assert_eq!(idx.nearest_splice_distance("1", 298, 298), Some(2));
    }

    #[test]
    fn test_mid_exon() {
        let idx = build_test_index();
        assert_eq!(idx.nearest_splice_distance("1", 150, 150), Some(50));
        assert_eq!(idx.nearest_splice_distance("1", 350, 350), Some(50));
    }

    #[test]
    fn test_span_distance() {
        // Sites 100, 200, 300, 400: a span is as far as its nearest base, and 0
        // with a boundary inside it.
        let idx = build_test_index();
        assert_eq!(idx.nearest_splice_distance("1", 193, 196), Some(4), "ends 4bp short of a right edge");
        assert_eq!(idx.nearest_splice_distance("1", 190, 205), Some(0), "crosses a right edge");
        assert_eq!(idx.nearest_splice_distance("1", 290, 300), Some(0), "ends on a left edge");
        assert_eq!(idx.nearest_splice_distance("1", 304, 305), Some(4), "starts 4bp into an exon");
        assert_eq!(idx.nearest_splice_distance("1", 230, 260), Some(30), "mid-intron: the nearer end");
        assert_eq!(idx.nearest_splice_distance("1", 20, 30), Some(70), "before the first site");
        assert_eq!(idx.nearest_splice_distance("1", 450, 460), Some(50), "after the last site");
    }

    #[test]
    fn test_unknown_chrom() {
        let idx = build_test_index();
        assert_eq!(idx.nearest_splice_distance("X", 100, 100), None, "no sentinel distance");
    }

    #[test]
    fn test_distance_lookup_normalizes_contig_naming() {
        // The index is keyed by normalized names (the GTF parser applies
        // normalize_contig); callers pass the input's contig, so a chr-named or
        // chrM-named variant must still find its annotation.
        let idx = build_test_index();
        assert_eq!(idx.nearest_splice_distance("chr1", 298, 298), Some(2));
        let mut chrom_map = HashMap::new();
        chrom_map.insert("MT".to_string(), 0u32);
        let mut splice_sites = HashMap::new();
        splice_sites.insert(0u32, vec![100i32, 200]);
        let mt = AnnotationIndex::new(HashMap::new(), vec![], splice_sites, HashMap::new(), chrom_map);
        for name in ["chrM", "M", "MT", "chrMT"] {
            assert_eq!(mt.nearest_splice_distance(name, 198, 198), Some(2), "{name}");
        }
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
