# Nextflow Examples

Common usage patterns for the gbcms Nextflow pipeline.

## DNA Mode (Default)

### Local with Docker

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    --outdir results \
    -profile docker
```

### SLURM Cluster

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    --outdir results \
    -profile slurm
```

## CRAM Files

CRAM files are transparently supported in both DNA and RNA modes. The `--fasta` parameter is required for CRAM decoding. The pipeline automatically discovers `.crai` index files.

```bash
nextflow run nextflow/main.nf \
    --input samplesheet_cram.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    --outdir results \
    -profile docker
```

## RNA Mode

### Basic RNA-seq

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.maf \
    --fasta reference.fa \
    --mode rna \
    --format maf \
    -profile docker
```

### RNA with Editing Database

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.maf \
    --fasta reference.fa \
    --mode rna \
    --format maf \
    --rna_editing_db /path/to/TABLE1_hg38_v3.txt \
    --enforce_strandedness \
    -profile docker
```

## MAF Output

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants input.maf \
    --fasta reference.fa \
    --format maf \
    -profile docker
```

## Strict Filtering

Enable additional read filters beyond the defaults (duplicates, secondary, supplementary, and QC-failed are already on by default):

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    --filter_improper_pair true \
    --filter_indel true \
    -profile docker
```

## UMI-Aware Counting

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    --umi_tag XM \
    -profile docker
```

## Multi-Sample MAF Filtering

Filter a multi-sample MAF so each BAM/CRAM only processes its own variants:

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants multi_sample.maf \
    --fasta reference.fa \
    --format maf \
    --filter_by_sample \
    -profile docker
```

With patient-level filtering via `tsb` column in samplesheet:

```csv
sample,bam,tsb
patient_A,/path/to/A.bam,P-0012345
patient_B,/path/to/B.cram,"P-0067890-T01,P-0067890-T02"
```

See [Samplesheet — Multi-Sample MAF Filtering](samplesheet.md#multi-sample-maf-filtering) for details.

## Resume Failed Run

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    -profile docker \
    -resume
```

## Custom Resource Limits

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    --max_cpus 8 \
    --max_memory 32.GB \
    -profile docker
```

## Related

- [Samplesheet](samplesheet.md) — Input format
- [Parameters](parameters.md) — All options
- [CLI Quick Start](../getting-started/quickstart.md) — For few samples
