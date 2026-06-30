# Nextflow Parameters

Complete reference for all pipeline parameters.

## Required Parameters

| Parameter | Description |
|:----------|:------------|
| `--input` | Path to [samplesheet CSV](samplesheet.md) |
| `--variants` | Path to [VCF/MAF](../reference/input-formats.md) variants file |
| `--fasta` | Reference FASTA (with .fai index). Also used for CRAM decoding. |

## Mode

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `--mode` | `dna` | Analysis mode: `dna` (cfDNA/somatic) or `rna` (RNA-seq with transcriptome-aware filtering) |

## Output Options

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `--outdir` | `results` | Output directory |
| `--format` | `vcf` | Output format (`vcf` or `maf`) |
| `--suffix` | `''` | Suffix for output filenames |
| `--column_prefix` | `''` | Prefix for gbcms count columns in MAF output |
| `--preserve_barcode` | `false` | Keep original Tumor_Sample_Barcode from input MAF |

## mFSD Options (DNA only)

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `--mfsd` | `false` | Enable mFSD analysis — adds 34 mFSD columns to MAF and 7 `MFSD_*` INFO fields to VCF. See [mFSD metrics](../reference/counting-metrics.md#mfsd). |
| `--mfsd_parquet` | `false` | Write a companion `<sample>.fsd.parquet` with raw per-variant fragment size arrays. Requires `--mfsd`. |
| `--mfsd_report` | `false` | Generate an interactive HTML report with per-variant fragment size distributions. Implies `--mfsd` and `--mfsd_parquet`. See [mFSD Report](../reference/mfsd-report.md). |
| `--mfsd_report_min_alt` | `3` | Minimum ALT fragment count to include a variant in the HTML report. |
| `--mfsd_report_max_variants` | `20` | Maximum variants in the HTML report. Use `-1` for no limit. |

## Quality & Filtering Options

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `--min_mapq` | `20` | Minimum mapping quality |
| `--min_baseq` | `20` | Minimum base quality |
| `--fragment_qual_threshold` | `10` | Quality margin for [fragment consensus](../reference/counting-metrics.md#fragment-counting) — when R1/R2 disagree on a non-INDEL variant, the higher-quality allele wins only if the difference exceeds this. INDEL conflicts with structural CIGAR evidence bypass this threshold. |
| `--context_padding` | `5` | Minimum flanking bases for [Phase 3 alignment](../reference/allele-classification.md#phase-3-alignment-fallback) (auto-increased in repeats) |
| `--adaptive_context` | `true` | Dynamically increase context padding in [tandem repeat regions](../reference/variant-normalization.md#adaptive-context-padding) |
| `--filter_duplicates` | `true` | Filter duplicate reads |
| `--filter_secondary` | `true` | Filter secondary alignments |
| `--filter_supplementary` | `true` | Filter supplementary alignments |
| `--filter_qc_failed` | `true` | Filter QC failed reads |
| `--filter_improper_pair` | `false` | Filter improperly paired reads |
| `--filter_indel` | `false` | Filter reads with indels |
| `--filter_by_sample` | `false` | Filter multi-sample MAF by `Tumor_Sample_Barcode` ([details](samplesheet.md#multi-sample-maf-filtering)) |
| `--show_normalization` | `false` | Add `norm_*` columns showing left-aligned coordinates in output |
| `--rescue_mnp` | `false` | Enable [MNP rescue pass](../reference/architecture.md#mnp-rescue-pass-rescue-mnp-v430) — decomposes MNPs into individual SNPs for re-counting when `ad=0` |
| `--rescue_mnp_threshold` | `1.0` | Maximum disc/len ratio for MNP rescue eligibility (0.0–1.0). `1.0` = all MNPs eligible (C++ compatible). `0.5` = conservative sparse-only mode. `0.0` = disable rescue eligibility (diagnostics still emitted). |

## UMI & BAQ Options

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `--umi_tag` | `''` | UMI BAM tag for deduplication (e.g., `XM`, `RX`). When set, reads sharing the same UMI are grouped as a single observation. |
| `--apply_baq` | `false` (DNA) / `true` (RNA) | Apply Base Alignment Quality downgrade. Reduces false positives near indels and splice junctions (CIGAR `N`). See [BAQ Options](../cli/dna.md#baq-options). |

## RNA-Specific Options

These parameters are only used when `--mode rna` is specified.

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `--rna_editing_db` | `''` | Path to [REDIportal](http://srv00.recas.ba.infn.it/atlas/index.html) editing database file (e.g., `TABLE1_hg38_v3.txt`). Flags ALT sites that overlap known A→I RNA editing positions. |
| `--enforce_strandedness` | `true` | Enforce strand-specific library prep. Disable with `false` for unstranded RNA-seq libraries (`--no-strandedness` equivalent). |
| `--strandedness` | `'reverse'` | RNA library strand protocol: `'reverse'` (dUTP/fr-firststrand, featureCounts `-s 2` — the FORTE default), `'forward'` (fr-secondstrand, `-s 1`), or `'unstranded'` (`-s 0`). Sets the read→transcript-strand fold used by `--enforce_strandedness` and ASJD strand-discordance; `'unstranded'` disables both. |

!!! tip "RNA mode defaults"
    RNA mode uses different PairHMM gap penalties by default (`gap_open=5e-3`, `gap_extend=0.25`) to tolerate RT-induced stutter at homopolymers. These can be overridden via the alignment backend parameters below.

## GTF Index Caching (RNA)

Parsing a full Ensembl GTF takes ~9s. Across a cohort that cost is otherwise paid once **per sample**. The pipeline caches the parsed index so the GTF is parsed **once for the whole cohort** — the per-sample annotation load drops from ~9s to ~0.05s.

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `--gtf_cache` | `true` | When `--gtf` is set, pre-build the GTF index once per cohort so per-sample tasks skip the parse. Set `false` to disable (each sample re-parses). |

**This is automatic — no extra wiring needed.** When `--mode rna` is run with a `--gtf` (and `--gtf_cache` is left `true`), the pipeline runs a single `GBCMS_BUILD_GTF_CACHE` process up front, using the cohort's `--gtf` and `--variants`. Every per-sample `GBCMS_RNA` task then receives the prebuilt cache directory (via the Nextflow DAG, so it is guaranteed ready first) and loads the index in ~0.05s instead of re-parsing.

**Why the up-front build matters.** Nextflow launches up to `queueSize` tasks at once. If the cache did not already exist, every concurrently-launched sample would cold-miss and re-parse the GTF — so a plain cache saves nothing until a later wave, and nothing at all for a cohort smaller than `queueSize`. Building it once before the fan-out lets every sample start warm.

The pre-build uses the cohort `--variants` file. With `--filter_by_sample`, each sample's variants are a per-sample subset; if a subset's chromosome set differs from the cohort's, that sample simply falls back to a normal parse (the cache is keyed on the variant chromosome set). Caching is best-effort throughout — a missing, corrupt, stale-version, or unwritable cache always falls back to a fresh parse and never affects counts.

## Alignment Backend (Advanced)

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `--alignment_backend` | `pairhmm` | Phase 3 alignment backend: `pairhmm` (WFA2 + PairHMM, default) or `sw` (Smith-Waterman). See [CLI Reference](../cli/dna.md#alignment-backend). |
| `--llr_threshold` | `2.3` | PairHMM log-likelihood ratio threshold for confident calls |
| `--gap_open_prob` | `1e-4` | PairHMM gap-open probability for non-repeat regions |
| `--gap_extend_prob` | `0.1` | PairHMM gap-extend probability for non-repeat regions |
| `--gap_open_prob_repeat` | `1e-2` | PairHMM gap-open probability for tandem repeat regions |
| `--gap_extend_prob_repeat` | `0.5` | PairHMM gap-extend probability for tandem repeat regions |

## Merge Options

These parameters control the multi-BAM merge step. Requires a `bam_type` column in the [samplesheet](samplesheet.md#multi-bam-type-merging).

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `--merge_counts` | `false` | Enable multi-BAM merge — combine per-type MAFs into a single output with type-prefixed count columns. Requires `bam_type` in samplesheet. |
| `--merge_add_combined` | `true` | When both `duplex` and `simplex` inputs are present, compute additive `simplex_duplex_*` combined columns (20 columns including strand bias). |
| `--merge_legacy_naming` | `false` | Use `t_{metric}_{type}` naming (genotype_variants compatible). |

## Resource Limits

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `--max_cpus` | `16` | Maximum CPUs per job |
| `--max_memory` | `128.GB` | Maximum memory per job |
| `--max_time` | `240.h` | Maximum runtime per job |

## Performance Benchmarks

Representative metrics from cfDNA duplex BAM samples:

| Sample Type | BAM Size | Variants | Runtime | CPUs |
|:------------|:---------|:---------|:--------|:-----|
| ctDNA (plasma) | 1.3 GB | 608 | ~25s | 4 |
| Plasma control | 776 MB | 608 | ~20s | 4 |

## Execution Profiles

| Profile | Description |
|:--------|:------------|
| `docker` | Local with Docker containers |
| `singularity` | HPC with Singularity |
| `slurm` | SLURM cluster with Singularity |
| `local` | No container (requires local install) |

## Advanced

!!! tip "`task.ext.args` — Arbitrary CLI Arguments"
    Any CLI option not exposed as a Nextflow parameter can be passed via `task.ext.args` in your config:

    ```groovy
    process {
        withName: GBCMS_DNA {
            ext.args = '--verbose'
        }
    }
    ```

    See [DNA CLI Reference](../cli/dna.md) or [RNA CLI Reference](../cli/rna.md) for all available options.

## Related

- [Samplesheet](samplesheet.md) — Input format
- [Examples](examples.md) — Usage patterns
- [DNA CLI Reference](../cli/dna.md) — Underlying DNA command
- [RNA CLI Reference](../cli/rna.md) — Underlying RNA command
