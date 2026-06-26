---
name: nextflow-pipeline
description: Reference for the gbcms Nextflow wrapper — module layout (dna/rna/filter_maf/merge/normalize/pipeline_summary), profiles and SLURM/iris resource routing, the samplesheet schema, multi-BAM merge options, and the places where Nextflow defaults DIVERGE from the CLI defaults. Use when editing anything under nextflow/ or debugging a pipeline run (OOM, queue routing, SLURM retries, output-name collisions).
---

# Nextflow Pipeline

## Layout

```
nextflow/
├── main.nf                       # entry
├── nextflow.config               # params + profiles + resource scopes
├── assets/samplesheet.csv
├── workflows/{dna,rna}.nf
└── modules/local/gbcms/
    ├── dna/main.nf               # GBCMS_DNA
    ├── rna/main.nf               # GBCMS_RNA
    ├── filter_maf/main.nf        # FILTER_MAF (include suffix in output name!)
    ├── merge/main.nf             # MERGE_COUNTS (multi-BAM)
    ├── normalize/main.nf         # standalone normalization
    └── pipeline_summary/main.nf  # PIPELINE_SUMMARY
```

## Resource / HPC notes (hard-won — see recent commits)

- Per-process resource blocks for `FILTER_MAF` and `PIPELINE_SUMMARY` exist to
  prevent OOM. Default memory bumped to 16GB.
- `iris` profile uses `process.queue` (not `params.partition`).
- Dynamic queue routing by resource request; **retry on SLURM null exit codes**.
- `FILTER_MAF` output names must include the `suffix` to avoid collision with
  multi-BAM samples.

## Multi-BAM merge

- `merge_counts` enables `MERGE_COUNTS` (requires `bam_type` in samplesheet).
- `merge_add_combined` computes `simplex_duplex_*` columns; `merge_legacy_naming`
  uses `t_{metric}_{type}` (genotype_variants compat).

## ⚠️ Nextflow vs CLI default divergences (trap)

`nextflow.config` is **not** a 1:1 mirror of CLI defaults. Known differences:
- `alignment_backend = 'pairhmm'` (NF) vs `sw` (CLI default).
- `min_mapq = 20`, `min_baseq = 20` set explicitly in NF.
- `enforce_strandedness = true` (NF) — note this is currently a no-op unless
  `gene_strand` is populated from the GTF.
- `apply_baq = false` in config, but RNA relies on the CLI default (on) unless
  explicitly set — modules pass `--apply-baq`/`--no-baq` conditionally.

When changing a CLI default or filter, update `nextflow.config` in the same change.

## Key Files
- `nextflow/nextflow.config`, `nextflow/modules/local/gbcms/*/main.nf`
- `nextflow/README.md`, `docs/nextflow/parameters.md`
