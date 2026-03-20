# Installation

## Quick Install

=== "PyPI (Recommended)"
    ```bash
    pip install gbcms
    gbcms --version
    ```
    
    !!! info "System Requirements"
        PyPI wheels require **glibc 2.34+** (Ubuntu 22.04+, RHEL 9+, Debian 12+).
        For older systems, see [Legacy Linux](#legacy-linux-rhel-8-hpc).

=== "Docker"
    ```bash
    docker pull ghcr.io/msk-access/gbcms:X.Y.Z
    docker run --rm ghcr.io/msk-access/gbcms:X.Y.Z gbcms --help
    ```

=== "From Source"
    ```bash
    git clone https://github.com/msk-access/gbcms.git
    cd gbcms
    pip install -e ".[dev]"
    ```

---

## Requirements

| Component | Requirement |
|:----------|:------------|
| Python | 3.10+ |
| OS | Linux (glibc 2.34+), macOS, Windows (WSL2) |
| Memory | 4GB+ (8GB for large BAMs) |

### For Nextflow Workflow

- Nextflow 21.10.3+
- Docker or Singularity

---

## Legacy Linux (RHEL 8 / HPC)

For RHEL 8, CentOS 8, or HPC systems with glibc < 2.34:

=== "Micromamba (Recommended)"
    ```bash
    # Create environment with build dependencies
    # Note: clangdev (not clang) provides headers needed by bindgen
    micromamba create -n gbcms_env python=3.13 clangdev rust -c conda-forge
    micromamba activate gbcms_env
    
    # Set libclang path for the Rust build
    export LIBCLANG_PATH=$CONDA_PREFIX/lib
    
    # Install from source
    git clone https://github.com/msk-access/gbcms.git
    cd gbcms
    pip install .
    ```

=== "Conda"
    ```bash
    conda create -n gbcms_env python=3.13 clangdev rust -c conda-forge
    conda activate gbcms_env
    export LIBCLANG_PATH=$CONDA_PREFIX/lib
    
    git clone https://github.com/msk-access/gbcms.git
    cd gbcms
    pip install .
    ```

=== "Singularity (HPC)"
    ```bash
    singularity pull docker://ghcr.io/msk-access/gbcms:X.Y.Z
    singularity exec gbcms_X.Y.Z.sif gbcms --help
    
    # With data binding
    singularity exec -B /path/to/data:/data gbcms_X.Y.Z.sif gbcms dna \
      --variants /data/variants.vcf --bam /data/sample.bam \
      --fasta /data/ref.fa --output-dir /data/results/
    ```

=== "Docker"
    ```bash
    docker pull ghcr.io/msk-access/gbcms:X.Y.Z
    docker run --rm -v $(pwd):/data ghcr.io/msk-access/gbcms:X.Y.Z gbcms --help
    ```

!!! note "Why not pip install?"
    PyPI wheels require glibc 2.34+. On RHEL 8 (glibc 2.28), pip falls back to 
    source compilation which requires Rust and clang headers. The conda environment
    provides these dependencies.

---

## Verification

```bash
# Check installation
gbcms --version
# Expected: X.Y.Z (your installed version)

# Test help
gbcms --help
```

---

## Docker Usage

```bash
docker run --rm \
  -v $(pwd):/data \
  ghcr.io/msk-access/gbcms:X.Y.Z \
  gbcms dna \
    --variants /data/variants.vcf \
    --bam /data/sample.bam \
    --fasta /data/reference.fa \
    --output-dir /data/results/
```

!!! tip "Docker Volume"
    Use `-v` to mount your data directory.

---

## Troubleshooting

See the **[full Troubleshooting Guide](../resources/troubleshooting.md)** for detailed solutions.

Quick checks:

```bash
# BAM index missing
samtools index sample.bam

# FASTA index missing
samtools faidx reference.fa

# Chromosome mismatch — compare names between FASTA, BAM, VCF/MAF
grep "^>" reference.fa | head -5
samtools view -H sample.bam | grep "^@SQ" | head -5
```

For glibc errors, installation failures, or Docker permission issues, see [Troubleshooting → Installation Issues](../resources/troubleshooting.md#installation-issues).

---

## Upgrade

```bash
pip install --upgrade gbcms
```

---

## Next Steps

- **[CLI Quick Start](quickstart.md)** — Command examples
- **[Nextflow Guide](../nextflow/index.md)** — HPC pipeline
