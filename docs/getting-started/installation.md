# Installation

## Quick Install

=== "PyPI"
    ```bash
    pip install gbcms
    gbcms --version
    ```

    !!! info "Two prebuilt wheels"
        A release publishes a wheel for **every CPython from 3.10 on Linux x86_64 with
        glibc 2.28+** (RHEL 8 and 9, Ubuntu 20.04+, Debian 10+), a wheel for **every CPython
        from 3.10 on Apple silicon Macs with macOS 11+**, and the source distribution.
        On any other platform, `pip install gbcms` compiles the source
        distribution on your machine, which needs a Rust toolchain and a C compiler.
        The wheel needs pip 20.3 or newer; on RHEL 8, install it into a virtual
        environment of a newer Python (see [RHEL 8](#legacy-linux-rhel-8-hpc)).
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
| OS | Linux (x86_64 wheel needs glibc 2.28+: RHEL 8 and newer), macOS (the wheel is for Apple silicon, macOS 11+; Intel builds from source), Windows through WSL2 |
| Memory | 4GB+ (8GB for large BAMs) |

### For Nextflow Workflow

- Nextflow 24.04.0 or later (the pipeline manifest enforces it: `nextflowVersion = '!>=24.04.0'`)
- Docker or Singularity

---

## Platform support

What a release builds is exactly what `.github/workflows/release.yml` builds: the
`linux` job (`ubuntu-latest`, `--release`, in the `manylinux_2_28` image, which is
AlmaLinux 8), the `macos` job (`macos-latest`, Apple silicon) and the `sdist` job. The
extension uses Python's stable ABI (`abi3-py310`), so each wheel job builds one wheel,
`cp310-abi3`, that every CPython from 3.10 installs. The `Wheels` workflow installs the
Linux wheel in Rocky Linux 8 and 9 containers and the macOS wheel on an Apple silicon
runner, with Python 3.10 to 3.14, and runs the test suite there. So `pip install gbcms`
gives:

| Your platform | What `pip install gbcms` does |
|:--------------|:------------------------------|
| Linux x86_64, glibc 2.28+ (RHEL 8 and newer), Python 3.10 to 3.14 | Installs the prebuilt wheel (`cp310-abi3`, `manylinux_2_28_x86_64`) |
| Linux with glibc below 2.28 (RHEL 7, CentOS 7), or Linux aarch64 | Builds the source distribution; see [Legacy Linux](#legacy-linux-rhel-8-hpc) |
| macOS on Apple silicon (arm64), macOS 11+, Python 3.10 to 3.14 | Installs the prebuilt wheel (`cp310-abi3`, `macosx_11_0_arm64`) |
| macOS on an Intel Mac | Builds the source distribution (not built in CI; Intel Macs are the last that macOS supports) |
| Windows (native) | Not built or tested; use WSL2, which is the Linux row |

The Docker image (`ghcr.io/msk-access/gbcms`) is built by the same workflow for
`linux/amd64` and `linux/arm64`, each on a runner of its own architecture, and joined
into one tag: a pull gets the image of your host's architecture, so on Apple silicon (and
on arm64 Linux) it runs natively, not under emulation. It carries Python 3.11. The
Nextflow `docker` profile still forces `--platform linux/amd64`, so a pipeline run
uses the same image everywhere; CI requires both architectures to write the same output
on a synthetic data set, which is what has to hold before that option is dropped.

### Building from source

Every row above that builds from source needs, on your machine:

- a **Rust toolchain** (`rustup`, stable);
- a **C compiler** (the bundled htslib is compiled from C);
- **libclang**, which `bindgen` uses for the htslib headers (set `LIBCLANG_PATH` if it
  is not found);
- the development headers for **zlib, bzip2 and xz**. The build needs no libcurl and
  no OpenSSL: gbcms reads local files only, so htslib is built without its remote-file
  support.

The release build installs `zlib-devel bzip2-devel xz-devel clang llvm-devel` (the
Debian names are `zlib1g-dev libbz2-dev liblzma-dev clang libclang-dev`). CI's macOS
jobs install `htslib`, `autoconf` and `automake` with Homebrew before building.

---

## Legacy Linux (RHEL 8 / HPC)

RHEL 8 (glibc 2.28) installs the prebuilt wheel with pip, on any Python from 3.10 to
3.14. RHEL 8's own Python is 3.6, so use a newer one, in a virtual environment (the
wheel needs pip 20.3 or newer):

```bash
sudo dnf install python3.11            # or: module load python/3.11, or conda
python3.11 -m venv gbcms_env
source gbcms_env/bin/activate
pip install --upgrade pip
pip install gbcms
gbcms --version
```

The other tabs are for glibc older than 2.28, for platforms without a wheel, and for
clusters that run containers.

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

!!! note "Singularity and Apptainer bind your HOME"
    Both runtimes bring your home folder and your environment into the container. The
    image sets `PYTHONNOUSERSITE=1`, so packages that you installed with `pip install
    --user` do not replace the ones in the image. A `PYTHONPATH` in your environment still
    reaches the container: unset it, or add `--cleanenv`, if it holds another gbcms or
    other copies of the dependencies. CI runs the image under Apptainer on every change
    that builds it, with a package in the caller's home as a trap, and requires output
    equal to the Docker run's.

!!! note "When pip builds from source"
    The one prebuilt wheel is for CPython 3.10 and newer and glibc 2.28+. On glibc
    older than 2.28, or on another platform, pip falls back to source compilation,
    which requires Rust and clang headers. The conda environments above provide these
    dependencies.

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
