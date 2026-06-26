---
name: nextflow-cli-default-divergence
description: The Nextflow config defaults are NOT a 1:1 mirror of the gbcms CLI defaults; keep them in sync when changing either.
metadata: 
  node_type: memory
  type: project
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
---

`nextflow/nextflow.config` sets defaults that diverge from the CLI:
- `alignment_backend = 'pairhmm'` (NF) vs `sw` (CLI default).
- `min_mapq = 20`, `min_baseq = 20` set explicitly in NF.
- `enforce_strandedness = true` (NF) — currently a no-op unless `gene_strand` is
  populated (see [[engine-output-aware]]).
- `apply_baq = false` in config, but RNA relies on the CLI default (on) unless set;
  modules pass `--apply-baq` / `--no-baq` conditionally (DNA off, RNA on).

**Why it matters:** A user running via Nextflow gets different behavior than the bare
CLI, and a CLI default change can silently desync the pipeline.

**How to apply:** When changing a CLI default, a filter, or a quality knob, update
`nextflow.config` (and the relevant `modules/local/gbcms/*/main.nf`) in the same change.
Recent HPC lessons (OOM resource blocks, iris `process.queue`, SLURM null-exit retries,
FILTER_MAF suffix-in-output-name) live in the `nextflow-pipeline` skill.
