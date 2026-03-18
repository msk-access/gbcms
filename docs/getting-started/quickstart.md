# CLI Quick Start

Process variants with the standalone CLI.

> **Many samples on HPC?** Use [Nextflow](../nextflow/index.md) instead.

---

## Basic Usage

### DNA / cfDNA

```bash
gbcms dna \
    --variants variants.vcf \
    --bam sample.bam \
    --fasta reference.fa \
    --output-dir results/
```

**Output:** `results/sample.vcf`

### RNA-seq

```bash
gbcms rna \
    --variants variants.vcf \
    --bam rna_sample:star_aligned.bam \
    --fasta reference.fa \
    --output-dir results/
```

**Output:** `results/rna_sample.vcf` (with RNA-specific columns)

---

## Common Options

### Output Format

```bash
# VCF output (default)
gbcms dna -v variants.vcf -b sample.bam -f ref.fa -o out/ --format vcf

# MAF output
gbcms dna -v variants.maf -b sample.bam -f ref.fa -o out/ --format maf
```

### Multiple Samples

```bash
# Using BAM list file
echo "sample1 /path/to/sample1.bam" > bam_list.txt
echo "sample2 /path/to/sample2.bam" >> bam_list.txt

gbcms dna \
    --variants variants.vcf \
    --bam-list bam_list.txt \
    --fasta reference.fa \
    --output-dir results/
```

### Custom Sample ID

```bash
gbcms dna \
    --variants variants.vcf \
    --bam MySample:sample.bam \
    --fasta reference.fa \
    --output-dir results/
```
**Output:** `results/MySample.vcf`

### Quality Filters

```bash
gbcms dna \
    --variants variants.vcf \
    --bam sample.bam \
    --fasta reference.fa \
    --output-dir results/ \
    --min-mapq 30 \
    --min-baseq 20 \
    --filter-duplicates \
    --filter-secondary
```

### Threading

```bash
gbcms dna ... --threads 8
```

---

## Complete Example

### DNA

```bash
gbcms dna \
    --variants variants.vcf \
    --bam TumorSample:tumor.bam \
    --fasta hg38.fa \
    --output-dir genotyped/ \
    --format vcf \
    --suffix .genotyped \
    --threads 8 \
    --min-mapq 30 \
    --min-baseq 20 \
    --filter-duplicates \
    --filter-secondary \
    --filter-supplementary
```

**Output:** `genotyped/TumorSample.genotyped.vcf`

### RNA with Editing Database

```bash
gbcms rna \
    --variants mutations.maf \
    --bam tumor_rna:aligned.bam \
    --fasta hg38.fa \
    --rna-editing-db TABLE1_hg38.txt.gz \
    --format maf \
    --output-dir results/
```

**Output:** `results/tumor_rna.maf` (with 5 RNA-specific columns)

---

## Docker

=== "DNA"

    ```bash
    docker run --rm -v $(pwd):/data ghcr.io/msk-access/gbcms:X.Y.Z \
        gbcms dna \
        --variants /data/variants.vcf \
        --bam /data/sample.bam \
        --fasta /data/reference.fa \
        --output-dir /data/results/
    ```

=== "RNA"

    ```bash
    docker run --rm -v $(pwd):/data ghcr.io/msk-access/gbcms:X.Y.Z \
        gbcms rna \
        --variants /data/variants.vcf \
        --bam /data/aligned.bam \
        --fasta /data/reference.fa \
        --output-dir /data/results/
    ```

---

## CLI Reference

```bash
gbcms dna --help
gbcms rna --help
```

| Option | Default | Description |
|:-------|:--------|:------------|
| `--variants` | Required | VCF or MAF file |
| `--bam` | Required | BAM file(s) |
| `--fasta` | Required | Reference FASTA |
| `--output-dir` | Required | Output directory |
| `--format` | vcf | Output format (vcf/maf) |
| `--min-mapq` | 20 (DNA) / 1 (RNA) | Minimum mapping quality |
| `--min-baseq` | 20 | Minimum base quality |
| `--threads` | 1 | Number of threads |

> 📖 See [DNA Reference](../cli/dna.md) and [RNA Reference](../cli/rna.md) for the complete list of options.

---

## Next Steps

- **[RNA Mode](../cli/rna.md)** — Transcriptome-aware counting
- **[Nextflow](../nextflow/index.md)** — Process many samples in parallel
- **[Architecture](../reference/architecture.md)** — How it works
- **[Allele Classification](../reference/allele-classification.md)** — How each variant type is counted
