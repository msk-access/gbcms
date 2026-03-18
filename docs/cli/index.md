# CLI Reference

The `gbcms` command-line interface provides two primary commands for variant counting — one for DNA/cfDNA and one for RNA-seq — plus a normalization utility.

## Commands

| Command | Description |
|:--------|:------------|
| [**dna**](dna.md) | Count alleles in DNA/cfDNA BAM files |
| [**rna**](rna.md) | Count alleles in RNA-seq BAMs with transcriptome-aware filtering |
| [**normalize**](normalize.md) | Standalone variant normalization (no counting) |

## Quick Example

```bash
# DNA/cfDNA counting
gbcms dna \
    --variants mutations.maf \
    --bam sample1:sample1.bam \
    --bam sample2:sample2.bam \
    --fasta reference.fa \
    --output-dir results/

# RNA-seq counting
gbcms rna \
    --variants mutations.vcf \
    --bam rna_sample:star_aligned.bam \
    --fasta reference.fa \
    --output-dir results/
```

## Getting Help

```bash
gbcms --help
gbcms dna --help
gbcms rna --help
gbcms normalize --help
```

## Related

- [Quick Start](../getting-started/quickstart.md) — Common usage patterns
- [Nextflow Pipeline](../nextflow/index.md) — For processing many samples
- [Input Formats](../reference/input-formats.md) — VCF/MAF specifications
