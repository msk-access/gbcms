# gbcms

**Complete orientation-aware counting system for genomic variants**

[![Tests](https://github.com/msk-access/gbcms/workflows/Tests/badge.svg)](https://github.com/msk-access/gbcms/actions)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/msk-access/gbcms)

## Features

- 🚀 **High Performance**: Rust-powered core engine with multi-threading
- 🧬 **Complete Variant Support**: SNP, MNP, insertion, deletion, and complex variants (DelIns, SNP+Indel)
- 🧪 **WFA + PairHMM Phase 3**: Pangenomic fast-path WFA alignment with PairHMM fallback for complex multi-allelic classification
- 📊 **Orientation-Aware**: Forward and reverse strand analysis with fragment counting
- 📏 **mFSD (Mutant Fragment Size Distribution)**: Per-allele cfDNA fragment size profiling with KS test and log-likelihood ratio
- 🔬 **Statistical Analysis**: Fisher's exact test for strand bias (read-level and fragment-level)
- 📁 **Flexible I/O**: VCF and MAF input/output formats
- 🎯 **Quality Filters**: 8 configurable read and quality filtering options with heuristic BAQ
- 🧬 **RNA Mode**: Transcriptome-aware counting with strandedness, splice detection, and A-to-I editing
- 🔗 **UMI Support**: Molecule-level deduplication with UMI-aware fragment grouping
- 🔧 **Normalize Command**: Standalone variant normalization (left-align + REF validation) without counting

## Installation

**Quick install:**
```bash
pip install gbcms
```

**From source (requires Rust):**
```bash
git clone https://github.com/msk-access/gbcms.git
cd gbcms
pip install .
```

**Docker:**
```bash
docker pull ghcr.io/msk-access/gbcms:X.Y.Z  # Replace X.Y.Z with latest from PyPI
```

> 💡 Find the latest version on [PyPI](https://pypi.org/project/gbcms/) or [GHCR](https://github.com/msk-access/gbcms/pkgs/container/gbcms).

📖 **Full documentation:** https://msk-access.github.io/gbcms/

---

## Usage

`gbcms` can be used in two ways:

### 🔧 Option 1: Standalone CLI (1-10 samples)

**Best for:** Quick analysis, local processing, direct control

```bash
gbcms dna \
    --variants variants.vcf \
    --bam sample1.bam \
    --fasta reference.fa \
    --output-dir results/
```

**Output:** `results/sample1.vcf`

**Learn more:**
- 📘 [CLI Quick Start](https://msk-access.github.io/gbcms/getting-started/quickstart/)
- 📖 [CLI Reference — DNA](https://msk-access.github.io/gbcms/cli/dna/)
- 📖 [CLI Reference — RNA](https://msk-access.github.io/gbcms/cli/rna/)
- 📖 [CLI Reference — Normalize](https://msk-access.github.io/gbcms/cli/normalize/)

---

### 🔄 Option 2: Nextflow Workflow (10+ samples, HPC)

**Best for:** Many samples, HPC clusters (SLURM), reproducible pipelines

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    --mode dna \
    -profile slurm
```

**Features:**
- ✅ Automatic parallelization across samples
- ✅ SLURM/HPC integration
- ✅ Container support (Docker/Singularity)
- ✅ Resume failed runs

**Learn more:**
- 🔄 [Nextflow Workflow Guide](https://msk-access.github.io/gbcms/nextflow/)
- 📋 [Usage Patterns Comparison](https://msk-access.github.io/gbcms/getting-started/)

---

## Which Should I Use?

| Scenario | Recommendation |
|----------|----------------|
| 1-10 samples, local machine | **CLI** |
| 10+ samples, HPC cluster | **Nextflow** |
| Quick ad-hoc analysis | **CLI** |
| Production pipeline | **Nextflow** |
| Need auto-parallelization | **Nextflow** |
| Full manual control | **CLI** |

---

## Quick Examples

### CLI: DNA Single Sample
```bash
gbcms dna \
    --variants variants.vcf \
    --bam tumor.bam \
    --fasta hg19.fa \
    --output-dir results/ \
    --threads 4
```

### CLI: RNA-seq
```bash
gbcms rna \
    --variants variants.vcf \
    --bam rna_sample:aligned.bam \
    --fasta hg19.fa \
    --rna-editing-db TABLE1_hg38.txt.gz \
    --output-dir results/
```

### CLI: Normalize Variants
```bash
gbcms normalize \
    --variants variants.vcf \
    --fasta hg19.fa \
    --output-dir results/
```

### CLI: Multiple Samples (Sequential)
```bash
gbcms dna \
    --variants variants.vcf \
    --bam-list samples.txt \
    --fasta hg19.fa \
    --output-dir results/
```

### Nextflow: Many Samples (Parallel)
```bash
# samplesheet.csv:
# sample,bam,bai
# tumor1,/path/to/tumor1.bam,
# tumor2,/path/to/tumor2.bam,

nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta hg19.fa \
    --mode dna \
    --outdir results \
    -profile slurm
```

---

## Documentation

📚 **Full Documentation:** https://msk-access.github.io/gbcms/

**Quick Links:**
- [Installation](https://msk-access.github.io/gbcms/getting-started/installation/)
- [CLI Quick Start](https://msk-access.github.io/gbcms/getting-started/quickstart/)
- [Nextflow Workflow](https://msk-access.github.io/gbcms/nextflow/)
- [CLI Reference — DNA](https://msk-access.github.io/gbcms/cli/dna/)
- [CLI Reference — RNA](https://msk-access.github.io/gbcms/cli/rna/)
- [CLI Reference — Normalize](https://msk-access.github.io/gbcms/cli/normalize/)
- [Input Formats](https://msk-access.github.io/gbcms/reference/input-formats/)
- [Output Formats](https://msk-access.github.io/gbcms/reference/output-formats/)
- [Architecture](https://msk-access.github.io/gbcms/reference/architecture/)

---

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development guidelines.

To contribute to documentation, see the [`gh-pages` branch](https://github.com/msk-access/gbcms/tree/gh-pages).

---

## Citation

If you use `gbcms` in your research, please cite:

> Shah, R. et al. (2026). *gbcms: A high-performance orientation-aware genotype counting system for genomic variants.* Available at: https://github.com/msk-access/gbcms

**BibTeX:**
```bibtex
@software{pygbcms,
  author       = {Shah, Ronak and contributors},
  title        = {gbcms: A high-performance orientation-aware genotype counting system for genomic variants},
  year         = {2026},
  url          = {https://github.com/msk-access/gbcms},
  note         = {GitHub repository}
}
```

---

## License

AGPL-3.0 - see [LICENSE](LICENSE) for details.

---

## Support

- 🐛 **Issues:** https://github.com/msk-access/gbcms/issues
- 💬 **Discussions:** https://github.com/msk-access/gbcms/discussions
