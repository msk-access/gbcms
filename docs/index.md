# gbcms

> **Get Base Counts Multi-Sample** — High-performance variant counting from BAM files

[![Version](https://img.shields.io/pypi/v/gbcms)](https://pypi.org/project/gbcms/)
[![Python](https://img.shields.io/pypi/pyversions/gbcms)](https://pypi.org/project/gbcms/)
[![License](https://img.shields.io/github/license/msk-access/gbcms)](https://github.com/msk-access/gbcms/blob/main/LICENSE)

## What It Does

GBCMS extracts **allele counts** and **variant metrics** at specified positions in BAM files:

```mermaid
block-beta
    columns 3
    VCF["📄 VCF/MAF\nVariant positions"]:1
    Engine["⚡ gbcms\nPython + Rust"]:1
    Counts["📊 Allele Counts\nDP · RD · AD · VAF"]:1
    BAM["🗂️ BAM Files\n(1 to N samples)"]:1
    space:1
    Metrics["🧬 Fragment Counts\nStrand bias · mFSD"]:1

    VCF --> Engine
    BAM --> Engine
    Engine --> Counts
    Engine --> Metrics
```

## Visual Overview

<figure markdown="span">
  ![gbcms overview poster](assets/posters/gbcms-overview-poster.jpg){ loading=lazy width="100%" }
  <figcaption>gbcms end-to-end pipeline — click to enlarge</figcaption>
</figure>

### Detailed Overview (PDF)

![Detailed overview of gbcms](assets/posters/High_Performance_cfDNA_Variant_Counting_cmp.pdf){ type=application/pdf style="min-height:75vh;width:100%" }

### Key Metrics

| Metric | Formula | Description |
|:-------|:--------|:------------|
| **VAF** | `AD / (RD + AD)` | Variant Allele Frequency |
| **Strand Bias** | Fisher's exact test | Detect sequencing artifacts |
| **Fragment Counts** | Deduplicated pairs | PCR-aware counting |

---

## Quick Start

```bash
# Install
pip install gbcms

# DNA/cfDNA counting
gbcms dna --variants variants.vcf --bam sample.bam --fasta ref.fa --output-dir results/

# RNA-seq counting
gbcms rna --variants variants.vcf --bam rna:aligned.bam --fasta ref.fa --output-dir results/
```

**→ [Full Installation Guide](getting-started/installation.md)** | **→ [CLI Examples](getting-started/quickstart.md)**

---

## Choose Your Workflow

```mermaid
flowchart TD
    Start(["What data?"]):::start
    Start -->|"DNA / cfDNA\nWGS / WES / Panel"| DNA(["gbcms dna"]):::dna
    Start -->|"RNA-seq\n(STAR-aligned, dUTP)"| RNA(["gbcms rna"]):::rna

    DNA --> NsamD{"Many samples?
≥10 BAMs"}
    RNA --> NsamR{"Many samples?
≥10 BAMs"}

    NsamD -->|"No"| CLI(["🖥️ CLI"]):::cli
    NsamD -->|"Yes"| HPC{"HPC / SLURM?"}
    NsamR -->|"No"| CLI
    NsamR -->|"Yes"| HPC

    HPC -->|"Yes"| NF(["🔷 Nextflow"]):::nf
    HPC -->|"No"| CLI

    classDef start fill:#9b59b6,color:#fff,stroke:#7d3c98,stroke-width:2px;
    classDef dna fill:#27ae60,color:#fff,stroke:#1e8449,stroke-width:2px;
    classDef rna fill:#3498db,color:#fff,stroke:#2471a3,stroke-width:2px;
    classDef cli fill:#e67e22,color:#fff,stroke:#bf6516,stroke-width:2px;
    classDef nf fill:#2c3e50,color:#fff,stroke:#1a2530,stroke-width:2px;
```

| Workflow | Best For | Guide |
|:---------|:---------|:------|
| **CLI** | 1-10 samples, local/single server | [Quick Start](getting-started/quickstart.md) |
| **Nextflow** | 10+ samples, HPC/SLURM | [Nextflow Guide](nextflow/index.md) |

---

## Architecture

Python/Rust hybrid for maximum performance:

See **[Architecture Reference →](reference/architecture.md)** for full diagrams covering system layers, data flow, genomic binning, coordinate system, config hierarchy, and end-to-end sequence.

**→ [Technical Details](reference/architecture.md)** | **→ [How It Works](reference/allele-classification.md)**

---

## Documentation

| Section | Description |
|:--------|:------------|
| [Getting Started](getting-started/index.md) | Installation and first run |
| [CLI Reference](cli/index.md) | Command-line usage |
| [Nextflow Pipeline](nextflow/index.md) | HPC workflow |
| [How It Works](reference/architecture.md) | Architecture, algorithms, and formats |
| [Development](development/developer-guide.md) | Contributing guide |
