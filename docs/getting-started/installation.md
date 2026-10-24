# Installation

## Quick Install

=== "PyPI"
    ```bash
    pip install gbcms
    gbcms --version
    ```

    !!! info "One prebuilt wheel"
        A release publishes one wheel, for **CPython 3.11 on Linux x86_64 with glibc
        2.34+** (Ubuntu 22.04+, RHEL 9+, Debian 12+), and the source distribution. On
        any other platform or Python, `pip install gbcms` compiles the source
        distribution on your machine, which needs a Rust toolchain and a C compiler.
        See [Platform support](#platform-support).

=== "Docker"
    ```bash
    docker pull ghcr.io/msk-access/gbcms:X.Y.Z
    docker run --rm ghcr.io/msk-access/gbcms:X.Y.Z gbcms --help
    ```

=== "From Source"
    ```bash
    git clone https://github.com/msk-access/gbcms.git
    cd gbcms
    pip install .   # builds the Rust extension: needs a Rust toolchain and a C compiler

    # Contributors: an editable install with the test, lint and docs tools
    # (a PEP 735 dependency group, which needs pip>=25.1)
    pip install --upgrade "pip>=25.1"
    pip install -e . --group dev
    ```

---

## Requirements

| Component | Requirement |
|:----------|:------------|
| Python | 3.10+ (CI runs 3.10, 3.11 and 3.14 on every change, and 3.10 to 3.14 monthly) |
| OS | Linux (x86_64 wheel needs glibc 2.34+), macOS (built from source), Windows through WSL2 |
| Memory | 4GB+ (8GB for large BAMs) |

### For Nextflow Workflow

- Nextflow 24.04.0 or later (the pipeline manifest enforces it: `nextflowVersion = '!>=24.04.0'`)
- Docker or Singularity

---

## Platform support

What a release builds is exactly what `.github/workflows/release.yml` builds: the
`linux` job (`ubuntu-latest`, `--release -i python3.11`, `manylinux_2_34`) and the
`sdist` job. So `pip install gbcms` gives:

| Your platform | What `pip install gbcms` does |
|:--------------|:------------------------------|
| Linux x86_64, glibc 2.34+, Python 3.11 | Installs the prebuilt wheel (`cp311`, `manylinux_2_34_x86_64`) |
| Linux x86_64, glibc 2.34+, another Python (3.10, 3.12 to 3.14) | Builds the source distribution |
| Linux with glibc below 2.34 (RHEL 8, CentOS 8), or Linux aarch64 | Builds the source distribution; see [Legacy Linux](#legacy-linux-rhel-8-hpc) |
| macOS, Intel or Apple Silicon | Builds the source distribution (CI builds and tests it on macOS with Python 3.12) |
| Windows (native) | Not built or tested; use WSL2, which is the Linux row |

The Docker image (`ghcr.io/msk-access/gbcms`) is built by the same workflow for
`linux/amd64` and carries Python 3.11.

### Building from source

Every row above that builds from source needs, on your machine:

- a **Rust toolchain** (`rustup`, stable);
- a **C compiler** (the bundled htslib is compiled from C);
- **libclang**, which `bindgen` uses for the htslib headers (set `LIBCLANG_PATH` if it
  is not found);
- the development headers for **zlib, bzip2, xz, libcurl and OpenSSL**.

The release build installs `zlib-devel bzip2-devel xz-devel libcurl-devel
openssl-devel clang llvm-devel` (the Debian names are `zlib1g-dev libbz2-dev
liblzma-dev libcurl4-openssl-dev libssl-dev clang libclang-dev`). CI's macOS jobs
install `htslib`, `autoconf` and `automake` with Homebrew before building.

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
    The one prebuilt wheel requires glibc 2.34+. On RHEL 8 (glibc 2.28), pip falls back
    to source compilation, which requires Rust and clang headers. The conda environment
    provides these dependencies.

---

## Verification

```bash
# Check installation
gbcms --version
# Expected: gbcms vX.Y.Z (commit)   (the version, and the commit of the build when known)

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

## Related

- [Quick Start](quickstart.md) — Common usage patterns with DNA and RNA examples
- [CLI Reference — DNA](../cli/dna.md) — Full option reference for `gbcms dna`
- [CLI Reference — RNA](../cli/rna.md) — Full option reference for `gbcms rna`
- [Nextflow Pipeline](../nextflow/index.md) — Running many samples in parallel on HPC
- [Troubleshooting](../resources/troubleshooting.md) — Installation issues and common errors
