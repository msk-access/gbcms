---
name: rna-annotation
description: "Reference for gbcms RNA mode — GTF-based transcript annotation (--gtf), the COITree AnnotationIndex, exon-boundary distance, per-transcript counts, and ASJD (allele-specific junction detection). Use when touching annotation/, rna.rs, or RNA output columns. Strandedness/sense-antisense gotchas live here."
---

# RNA Annotation & ASJD

## GTF-Based Annotation (`--gtf`)

When `--gtf` is provided in RNA mode, 17 additional columns are appended:
- 1 exon boundary distance (`exon_boundary_dist`)
- 2 per-transcript counts (`transcript_read_counts`, `transcript_fragment_counts`)
- 14 ASJD fields (below)

### Architecture

```
GTF file → gtf::parse_gtf() → AnnotationIndex
                                  ├── exon_trees: HashMap<u32, COITree>
                                  ├── exons: Vec<ExonRecord>
                                  ├── splice_sites: HashMap<u32, Vec<i32>>
                                  ├── transcript_introns: HashMap<String, TranscriptIntrons>
                                  └── chrom_map: HashMap<String, u32>
```

Built once in `count_bam_binned()`, shared via `Arc` across rayon workers.

### ASJD (Allele-Specific Junction Detection)

14 fields per variant comparing the splice junction the REF and the ALT reads use
(`AsjdResult` in engine.rs, fields on `BaseCounts` in types.rs):
- test: `asjd_flag` (Fisher exact p < 0.05), `asjd_pval`, `asjd_qval` (BH, set after counting)
- dominant junction per allele: `asjd_ref_junction` / `asjd_alt_junction` ("start-end"),
  `asjd_ref_motif` / `asjd_alt_motif` (GT-AG, GC-AG, AT-AC, OTHER),
  `asjd_ref_known` / `asjd_alt_known` (matches an annotated intron)
- support: `asjd_n_ref_junc` / `asjd_n_alt_junc` (reads on the dominant junction),
  `asjd_n_ref_total` / `asjd_n_alt_total` (that allele's junction-spanning fragments)
- `asjd_diagnostic`: semicolon-separated flags (`LOW_REF_JUNC`, `LOW_ALT_JUNC`,
  `NOVEL_ALT_JUNC`, `STRAND_DISCORDANT`, ...)

Junction counts are **per fragment** (`JunctionTally` in `detect_asjd` dedups a
molecule's R1+R2 by QNAME hash to one vote) — required so the strand-discordance and
Fisher tests see independent observations, not mates. `LOW_REF_JUNC` and `LOW_ALT_JUNC`
gate on the `asjd_n_ref_total` and `asjd_n_alt_total` counts; `STRAND_DISCORDANT` and the
Fisher test gate on the `asjd_n_ref_junc` and `asjd_n_alt_junc` counts.

### Gotchas (verify before trusting RNA outputs)
- `gene_strand` must be populated from the GTF or `is_sense_strand` returns true
  for all reads → `antisense_depth` stays 0 and `enforce_strandedness` is a no-op.
- Under enforcement an antisense read is classified, tallied in `antisense_depth`
  (first-class REF/ALT over the anchor, as the sense tally), then dropped before any
  count (`antisense_excluded` in the per-read trace). `STRAND_DISCORDANT` only fires
  with `--no-strandedness` where the gene strand is resolved. `strand_at` uses exons
  only, so intronic loci (splice sites) and opposite-strand overlaps have no strand:
  enforcement is a no-op there (#185).
- Library strand protocol is selectable via `--strandedness` (`reverse` default /
  `forward` / `unstranded`); the read→transcript fold lives in
  `rna.rs::read_transcript_strand(record, Strandedness)`. `reverse` = dUTP/`-s 2`
  (FORTE default). Wrong protocol → every junction looks strand-discordant. Confirm
  with RSeQC `infer_experiment` or the upstream `featureCounts -s` flag.
- The splice-motif classifier must account for minus-strand genes (revcomp), else
  canonical junctions read as OTHER / NON_CANONICAL.
- Contig naming: `chr`-strip alone does not reconcile `MT` vs `chrM`.

### Key Files
- `rust/src/annotation/mod.rs`: AnnotationIndex, COITree queries
- `rust/src/annotation/gtf.rs`: GTF parser
- `rust/src/counting/engine.rs`: resolves `gene_strand` from the index, `detect_asjd`,
  `JunctionTally`, per-transcript counts
- `rust/src/counting/rna.rs`: `is_sense_strand`, `read_transcript_strand`, splice-junction helpers
- `src/gbcms/cli.py`: `--gtf` CLI option (DNA + RNA commands)
- `docs/reference/rna-annotation.md`: user documentation
