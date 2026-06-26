---
name: rna-annotation
description: Reference for gbcms RNA mode — GTF-based transcript annotation (--gtf), the COITree AnnotationIndex, exon-boundary distance, per-transcript counts, and ASJD (allele-specific junction detection). Use when touching annotation/, rna.rs, or RNA output columns. Strandedness/sense-antisense gotchas live here.
---

# RNA Annotation & ASJD

## GTF-Based Annotation (`--gtf`)

When `--gtf` is provided in RNA mode, 17 additional columns are appended:
- 1 exon boundary distance (`exon_boundary_dist`)
- 2 per-transcript counts (`transcript_read_counts`, `transcript_variant_counts`)
- 14 ASJD fields (7 metrics × 2 alleles)

### Architecture

```
GTF file → gtf::parse_gtf() → AnnotationIndex
                                  ├── exon_trees: HashMap<u32, COITree>
                                  ├── splice_sites: HashMap<u32, Vec<i32>>
                                  ├── transcript_introns: HashMap<String, TranscriptIntrons>
                                  └── chrom_map: HashMap<String, u32>
```

Built once in `count_bam_binned()`, shared via `Arc` across rayon workers.

### ASJD (Allele-Specific Junction Detection)

6 diagnostic flags per variant comparing splice junction usage:
- `asjd_n_ref_known` / `asjd_n_ref_novel` / `asjd_n_alt_known` / `asjd_n_alt_novel`
- `asjd_n_ref_total` / `asjd_n_alt_total`
- `asjd_diagnostic`: semicolon-separated diagnostic string

### Gotchas (verify before trusting RNA outputs)
- `gene_strand` must be populated from the GTF or `is_sense_strand` returns true
  for all reads → `antisense_depth` stays 0 and `enforce_strandedness` is a no-op.
- The splice-motif classifier must account for minus-strand genes (revcomp), else
  canonical junctions read as OTHER / NON_CANONICAL.
- Contig naming: `chr`-strip alone does not reconcile `MT` vs `chrM`.

### Key Files
- `rust/src/annotation/mod.rs`: AnnotationIndex, COITree queries
- `rust/src/annotation/gtf.rs`: GTF parser
- `rust/src/counting/engine.rs`: P4a/P4b/P4c integration
- `src/gbcms/cli.py`: `--gtf` CLI option (DNA + RNA commands)
- `docs/reference/rna-annotation.md`: user documentation
