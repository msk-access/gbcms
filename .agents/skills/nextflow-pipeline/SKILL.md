---
name: nextflow-pipeline
description: "The Nextflow wrapper: module layout, profiles and SLURM/iris routing, samplesheet schema, multi-BAM merge, and where Nextflow defaults differ from the CLI. Use when editing anything under nextflow/ or debugging a pipeline run (OOM, queue routing, retries, output-name collisions)."
---

# Nextflow Pipeline

## Layout

```
nextflow/
├── main.nf                       # entry
├── nextflow.config               # params + profiles + resource scopes
├── assets/samplesheet.csv
├── workflows/{dna,rna}.nf
└── modules/local/
    ├── utils/main.nf             # asBool() and other helpers
    └── gbcms/
        ├── dna/main.nf           # GBCMS_DNA
        ├── rna/main.nf           # GBCMS_RNA
        ├── filter_maf/main.nf    # FILTER_MAF (include suffix in output name!)
        ├── merge/main.nf         # MERGE_COUNTS (multi-BAM)
        ├── normalize/main.nf     # standalone normalization
        └── pipeline_summary/main.nf  # PIPELINE_SUMMARY
```

## Resource / HPC notes (`process {}` in nextflow.config)

- `GBCMS_DNA|GBCMS_RNA`: 4 CPUs, 16 GB and 2 h, each times the attempt; the first attempt
  runs on `cmobic_short`, retries on `cmobic_cpu` (7-day limit). The process default is 6 GB.
- Own blocks exist so large inputs do not OOM: `MERGE_COUNTS` (4 GB x attempt), `FILTER_MAF`
  (8 GB x attempt, large multi-sample MAFs), `PIPELINE_SUMMARY` (2 GB).
- Retry (up to 3) on exit status null (a SLURM kill), 137, 143, 104, 134, 139 and 255.
- The `iris` profile sets `process.queue` (a bare `cpu` queue is rejected there), not
  `params.partition`; per-process `queue` directives still route first.
- `FILTER_MAF` output names must include the `suffix` (multi-BAM samples would collide).

## Samplesheet (`--input`; `docs/nextflow/samplesheet.md`)

Columns `sample,bam,bai,suffix,bam_type` (and an optional `tsb`); `bai`, `suffix` and `bam_type` may be empty. Output
suffix: the row's `suffix`, else `-<bam_type>`, else `params.suffix`. Rows sharing a `sample`
need a `suffix` or `bam_type`, or the workflow stops (their outputs would overwrite each other).

## Multi-BAM merge

- `merge_counts` enables `MERGE_COUNTS` (requires `bam_type` in samplesheet).
- `merge_add_combined` computes `simplex_duplex_*` columns; `merge_legacy_naming`
  uses `t_{metric}_{type}` (genotype_variants compat).

## ⚠️ Nextflow vs CLI default divergences (trap)

`nextflow.config` is **not** a 1:1 mirror of the CLI defaults, and the modules pass most of
them explicitly, so a CLI default change does not reach a pipeline run on its own.

Matches the CLI today (do not "fix"): `alignment_backend = 'pairhmm'` (the CLI default for
both DNA and RNA; `sw` is the opt-in), `min_baseq = 20`, the filter flags (always passed in
their explicit on/off form), `enforce_strandedness = true`, `strandedness = 'reverse'`,
`library_type = 'capture'`.

Where the DNA and RNA CLI defaults differ, RNA has its own param: `rna_min_mapq` (1),
`rna_apply_baq` (on), `rna_gap_*` (5e-3/0.25); `min_mapq` (20), `apply_baq` (off) and
`gap_*` are DNA's. `tests/test_nextflow_cli_defaults_contract.py` pins every module's
param defaults to its CLI's.

Caveat, not a divergence: `enforce_strandedness` does nothing for a variant with no
resolvable gene strand (no `--gtf`, an intergenic locus, genes of both strands over it); the
engine warns once, naming them.

When changing a CLI default or filter, update `nextflow.config` in the same change.

## ⚠️ Strict syntax (Nextflow ≥26.04 — the default parser)

The code must parse under both the 25.x legacy parser and the ≥26.04 strict parser.
The strict parser is stricter than Groovy — do **not** reintroduce any of these
(each is a hard parse error, not a warning):

- **No `while` loops** and **no `for` loops** — use `.readLines()`, `.each`, `.findAll`,
  etc. (`hasData` in `main.nf` reads this way).
- **No assignment inside a condition** (`while ((x = next()) != null)`), **no postfix
  `++`/`--`** — use `x += 1`.
- **No statements at a script's top level** — input validation, `if`-guards, and
  variable assignments must live inside `workflow {}` / a process / a function.
  Use `error 'msg'`, not the removed `exit`.
- **config files:** no function definitions (`def check_max(){}` → use
  `process.resourceLimits = [cpus:…, memory:…, time:…]`), and **no `try/catch` or
  `if` at the config top level** (this is why the remote nf-core `includeConfig` was
  removed — a fetch failure can no longer be caught; local `iris`/`slurm` profiles
  are self-sufficient, users add institutional configs via `-c`).
- **Deprecations (warn now, error later):** `Channel.x()` → `channel.x()`; declare
  explicit closure params instead of implicit `it`; prefix unused closure params `_`.
- **CLI params are Strings on ≥26.04** — `--flag false` arrives as the truthy String
  `"false"` (25.x coerced it to `Boolean false`), so `params.x ?` silently inverts.
  Route every boolean param through `asBool()` (`include { asBool } from
  '<...>/modules/local/utils/main'`) and cast numerics at use sites (see
  `resourceLimits` in nextflow.config). Also: never emit a bool CLI flag by
  omission-when-false if the gbcms CLI default is true — pass the explicit
  `--no-<flag>` form (see the filter block in the dna/rna modules).

**Always lint before merging** (bootstrap 26.04 once, it caches):
```bash
NXF_VER=26.04.0 nextflow lint nextflow/     # must report 0 errors
```

## Key Files
- `nextflow/nextflow.config`, `nextflow/modules/local/gbcms/*/main.nf`
- `nextflow/README.md`, `docs/nextflow/parameters.md`
