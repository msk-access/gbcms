# gbcms Nextflow Workflow

This directory contains a Nextflow workflow for running `gbcms` on multiple samples in parallel.

## Prerequisites

1. [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (>=22.10.1)
2. [Docker](https://docs.docker.com/engine/installation/) or [Singularity](https://www.sylabs.io/guides/3.0/user-guide/)

## Quick Start

### 1. Prepare a samplesheet (CSV format)

**Basic samplesheet:**
```csv
sample,bam,bai
sample1,/path/to/sample1.bam,/path/to/sample1.bam.bai
sample2,/path/to/sample2.bam,
```

**With per-sample suffix (for multiple BAM types per sample):**
```csv
sample,bam,bai,suffix
sample1,/path/to/sample1.duplex.bam,,-duplex
sample1,/path/to/sample1.simplex.bam,,-simplex
sample1,/path/to/sample1.unfiltered.bam,,-unfiltered
sample2,/path/to/sample2.bam,,
```

**Output:**
- `sample1-duplex.vcf`
- `sample1-simplex.vcf`  
- `sample1-unfiltered.vcf`
- `sample2.vcf` (or `sample2{--suffix}.vcf` if global suffix set)

**Notes:**
- `bai` column is optional - will auto-discover `<bam>.bai` if not provided
- `suffix` column is optional - per-row suffix overrides global `--suffix` parameter

### 2. Run the workflow

**DNA mode (default) — local with Docker:**
```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    --outdir results \
    -profile docker
```

**RNA mode:**
```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.maf \
    --fasta reference.fa \
    --mode rna \
    --format maf \
    -profile docker
```

**SLURM cluster with Singularity:**
```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    --outdir results \
    -profile slurm
```

## Parameters

### Required
| Parameter | Description |
|-----------|-------------|
| `--input` | Path to samplesheet CSV |
| `--variants` | Path to VCF/MAF variants file |
| `--fasta` | Path to reference FASTA (with .fai index) |

### Mode
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--mode` | Analysis mode: `dna` or `rna` | `dna` |

### Optional
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--outdir` | Output directory | `results` |
| `--format` | Output format (`vcf` or `maf`) | `vcf` |
| `--suffix` | Suffix to append to output filenames | `''` (empty) |
| `--column_prefix` | Prefix for gbcms count columns in MAF output | `''` (empty) |
| `--preserve_barcode` | Keep original Tumor_Sample_Barcode from input MAF | `false` |
| `--show_normalization` | Add `norm_*` columns showing left-aligned coordinates | `false` |
| `--min_mapq` | Minimum mapping quality (DNA) | `20` |
| `--rna_min_mapq` | Minimum mapping quality (RNA) — `1` matches the CLI and keeps STAR multi-mapper primaries (MAPQ 3/1); raise to 20 for unique-only | `1` |
| `--min_baseq` | Minimum base quality | `20` |
| `--fragment_qual_threshold` | Quality margin for fragment consensus | `10` |
| `--context_padding` | Minimum flanking bases for alignment | `5` |
| `--adaptive_context` | Auto-increase context in tandem repeats | `true` |
| `--umi_tag` | UMI BAM tag for deduplication (e.g., `XM`, `RX`) | `''` (disabled) |
| `--apply_baq` | Apply BAQ recalibration (config default off; the RNA CLI enables it by default) | `false` |
| `--filter_duplicates` | Filter duplicate reads | `true` |
| `--filter_secondary` | Filter secondary alignments | `true` |
| `--filter_supplementary` | Filter supplementary alignments | `true` |
| `--filter_qc_failed` | Filter QC failed reads | `true` |
| `--filter_improper_pair` | Filter improperly paired reads | `false` |
| `--filter_indel` | Filter reads with indels | `false` |
| `--filter_by_sample` | Filter multi-sample MAF by Tumor_Sample_Barcode | `false` |
| `--alignment_backend` | Phase 3 alignment: `pairhmm` (default), `sw` (Smith-Waterman), or `hmm` (alias for `pairhmm`) | `pairhmm` |

### RNA-Specific
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--rna_editing_db` | Path to REDIportal editing database (TABLE1, tab-delimited) | `''` (disabled) |
| `--strandedness` | Library protocol: `reverse` (dUTP/`-s2`), `forward` (`-s1`), or `unstranded` (`-s0`) | `reverse` |
| `--enforce_strandedness` | Filter reads to the transcript's sense strand (requires `--gtf` for `gene_strand`). `--strandedness unstranded` disables it. | `true` |
| `--library_type` | `capture` (default) or `amplicon`. Amplicon treats R1/R2 as independent observations (no fragment consensus) and disables strandedness. | `capture` |
| `--gtf` | GTF for exon-boundary / per-transcript / ASJD annotation and `gene_strand` back-fill | `''` (disabled) |
| `--gtf_cache` | Pre-build the GTF index once per cohort (when `--gtf` is set) so per-sample tasks skip the ~9s parse | `true` |

### Feature Columns & Merge
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--mfsd` | Emit the 41 mFSD MAF columns / 13 VCF `MFSD_*` INFO fields (fragment-size distribution) | `false` |
| `--mfsd_parquet` | Also write the mFSD size arrays to Parquet (requires `--mfsd`) | `false` |
| `--mfsd_report` | Generate the HTML mFSD report | `false` |
| `--rescue_mnp` | Enable the MNP-rescue second pass | `false` |
| `--rescue_mnp_threshold` | Discordance threshold for MNP rescue | `1.0` |
| `--merge_counts` | Enable the `MERGE_COUNTS` process (multi-BAM merge; needs `bam_type` in the samplesheet) | `false` |
| `--merge_add_combined` | Compute `simplex_duplex_*` combined columns during merge | `true` |
| `--merge_legacy_naming` | Use `t_{metric}_{type}` naming (genotype_variants compat) | `false` |

See [Full Parameter Reference](../docs/nextflow/parameters.md) for all options including PairHMM tuning parameters.

### Resource Limits
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--max_cpus` | Maximum CPUs per process | `16` |
| `--max_memory` | Maximum memory per process | `128.GB` |
| `--max_time` | Maximum time per process | `240.h` |

## Profiles

- `-profile docker`: Use Docker containers (recommended for local)
- `-profile singularity`: Use Singularity images (recommended for HPC)
- `-profile slurm`: Run on SLURM cluster with Singularity (queue: `cmobic_cpu`)
- `-profile local`: No container (requires local `gbcms` install)
- `-profile debug`: Print hostname for debugging

## Customizing for Your Cluster

Edit `nextflow/nextflow.config` to customize the SLURM profile:

```groovy
slurm {
    process.executor       = 'slurm'
    process.queue          = 'your_queue_name'  // Change this
    singularity.enabled    = true
    singularity.autoMounts = true
    docker.enabled         = false
}
```

## Output

Results are published to `${params.outdir}/gbcms/`:
- VCF files: `<sample>.vcf`
- MAF files: `<sample>.maf`

(Output filenames use the sample id as the prefix; set `--suffix` to add an infix.)

Pipeline info and logs are in `${params.outdir}/pipeline_info/`.

## Pipeline Modules

| Module | Description |
|--------|-------------|
| `GBCMS_DNA` | DNA allele counting via `gbcms dna` |
| `GBCMS_RNA` | RNA allele counting via `gbcms rna` |
| `GBCMS_BUILD_GTF_CACHE` | Pre-build the shared GTF index once per cohort via `gbcms build-gtf-cache` (RNA + `--gtf`) |
| `GBCMS_NORMALIZE` | Variant normalization via `gbcms normalize` |
| `MERGE_COUNTS` | Merge per-BAM counts into simplex/duplex combined columns (`--merge_counts`) |
| `FILTER_MAF` | Pre-filter multi-sample MAF by sample |
| `PIPELINE_SUMMARY` | Aggregate filtering statistics |
