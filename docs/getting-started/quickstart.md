# Quick Start

Get up and running in minutes. All examples below assume gbcms is [installed](installation.md).

!!! tip "Many samples on HPC?"
    Use the [Nextflow pipeline](../nextflow/index.md) instead of the CLI for parallel processing on a cluster.

---

## Basic Usage

=== "DNA / cfDNA"

    ```bash
    gbcms dna \
        --variants variants.vcf \
        --bam sample.bam \
        --fasta reference.fa \
        --output-dir results/
    ```

    **Output:** `results/sample.vcf`

=== "RNA-seq"

    ```bash
    gbcms rna \
        --variants variants.vcf \
        --bam rna_sample:star_aligned.bam \
        --fasta reference.fa \
        --output-dir results/
    ```

    **Output:** `results/rna_sample.vcf` (with 5 RNA-specific columns)

---

## Output Format

=== "DNA / cfDNA"

    ```bash
    # VCF output (default)
    gbcms dna -v variants.vcf -b sample.bam -f ref.fa -o out/

    # MAF output — preserves all input columns, appends gbcms counts
    gbcms dna -v variants.maf -b sample.bam -f ref.fa -o out/ --format maf
    ```

=== "RNA-seq"

    ```bash
    # VCF output (default)
    gbcms rna -v variants.vcf -b rna:aligned.bam -f ref.fa -o out/

    # MAF output — includes 5 RNA-specific columns
    gbcms rna -v variants.maf -b rna:aligned.bam -f ref.fa -o out/ --format maf
    ```

---

## Multiple Samples

=== "DNA / cfDNA"

    ```bash
    # BAM list file: each line is "sample_name /path/to/sample.bam"
    echo "tumor   /path/to/tumor.bam"  > bam_list.txt
    echo "normal  /path/to/normal.bam" >> bam_list.txt

    gbcms dna \
        --variants variants.vcf \
        --bam-list bam_list.txt \
        --fasta reference.fa \
        --output-dir results/
    ```

=== "RNA-seq"

    ```bash
    echo "rna_tumor /path/to/rna.bam" > bam_list.txt

    gbcms rna \
        --variants variants.vcf \
        --bam-list bam_list.txt \
        --fasta reference.fa \
        --output-dir results/
    ```

---

## Quality Filters

=== "DNA / cfDNA"

    ```bash
    gbcms dna \
        --variants variants.vcf \
        --bam sample.bam \
        --fasta reference.fa \
        --output-dir results/ \
        --min-mapq 30 \
        --min-baseq 20 \
        --filter-duplicates \
        --filter-secondary \
        --filter-supplementary
    ```

=== "RNA-seq"

    ```bash
    # Both DNA and RNA filter secondary, supplementary, and QC-failed by default
    gbcms rna \
        --variants variants.vcf \
        --bam rna:aligned.bam \
        --fasta reference.fa \
        --output-dir results/ \
        --min-baseq 20   # MAPQ default is 1 with NH:i:1 rescue
    ```

---

## Complete Example

=== "DNA with mFSD"

    ```bash
    gbcms dna \
        --variants variants.vcf \
        --bam TumorSample:tumor.bam \
        --fasta hg38.fa \
        --output-dir genotyped/ \
        --format maf \
        --suffix .genotyped \
        --threads 8 \
        --min-mapq 30 \
        --min-baseq 20 \
        --filter-duplicates \
        --filter-secondary \
        --filter-supplementary \
        --mfsd \
        --mfsd-parquet
    ```

    **Output:**
    - `genotyped/TumorSample.genotyped.maf` — allele counts + 41 mFSD columns
    - `genotyped/TumorSample.genotyped.fsd.parquet` — raw fragment size arrays

=== "RNA with Editing DB"

    ```bash
    gbcms rna \
        --variants mutations.maf \
        --bam tumor_rna:aligned.bam \
        --fasta hg38.fa \
        --rna-editing-db TABLE1_hg38.txt.gz \
        --format maf \
        --threads 8 \
        --output-dir results/
    ```

    **Output:** `results/tumor_rna.maf` — standard counts + 5 RNA columns:
    `rna_sense_depth`, `rna_antisense_depth`, `rna_alt_sense_count`,
    `rna_editing_site`, `rna_splice_spanning`

=== "Unstranded RNA-seq"

    ```bash
    gbcms rna \
        --variants variants.vcf \
        --bam unstranded:aligned.bam \
        --fasta reference.fa \
        --no-strandedness \
        --output-dir results/
    ```

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
        --bam /data/rna:aligned.bam \
        --fasta /data/reference.fa \
        --output-dir /data/results/
    ```

---

## Common CLI Options

| Option | Default | Description |
|:-------|:--------|:------------|
| `--variants` | Required | VCF or MAF file |
| `--bam` | Required | BAM file(s). Prefix with `name:` to set sample ID |
| `--bam-list` | — | File with BAM paths (one per line) |
| `--fasta` | Required | Reference FASTA |
| `--output-dir` | Required | Output directory |
| `--format` | `vcf` | Output format (`vcf` or `maf`) |
| `--min-mapq` | 20 (DNA) / 1 (RNA) | Minimum mapping quality |
| `--min-baseq` | 20 | Minimum base quality |
| `--threads` | 1 | Number of threads |

> 📖 Full option reference: [DNA CLI](../cli/dna.md) · [RNA CLI](../cli/rna.md)

---

## Related

- [DNA CLI Reference](../cli/dna.md) — mFSD, UMI, alignment backend options
- [RNA CLI Reference](../cli/rna.md) — Strandedness, editing DB, splice junction options
- [Nextflow Pipeline](../nextflow/index.md) — Process many samples in parallel on HPC
- [Allele Classification](../reference/allele-classification.md) — How the counting engine works
- [Troubleshooting](../resources/troubleshooting.md) — Common issues and solutions
