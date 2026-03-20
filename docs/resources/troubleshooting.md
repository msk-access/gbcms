# Troubleshooting

Common issues and solutions for gbcms. Issues are grouped by phase — work top-to-bottom for systematic diagnosis.

---

## Installation Issues

??? question "Rust compilation fails: `error: linker 'cc' not found`"

    Install system build tools:

    === "macOS"
        ```bash
        xcode-select --install
        ```

    === "Ubuntu / Debian"
        ```bash
        apt-get install build-essential
        ```

    === "CentOS / RHEL"
        ```bash
        yum groupinstall "Development Tools"
        ```

??? question "htslib linking error: `cannot find -lhts`"

    The Rust backend compiles htslib statically. Ensure compression library headers are present:

    ```bash
    apt-get install zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev
    ```

??? question "`GLIBC_2.34 not found` when installing via pip"

    Your system has an older glibc (RHEL 8 / CentOS 8 / HPC). PyPI wheels require glibc 2.34+.
    Use the [Legacy Linux](../getting-started/installation.md#legacy-linux-rhel-8-hpc) install path
    (micromamba + clangdev + source install).

??? question "`LIBCLANG_PATH` not set error during build"

    Required for `bindgen` (htslib headers). Set:
    ```bash
    export LIBCLANG_PATH=$CONDA_PREFIX/lib
    ```
    Make sure you have `clangdev` (not `clang`) from conda-forge in your environment.

---

## BAM / Input Issues

??? question "BAI index not found for `sample.bam`"

    Create or verify the index:
    ```bash
    samtools index sample.bam
    # Expected: creates sample.bam.bai (or sample.bai)
    samtools quickcheck sample.bam  # Sanity check the BAM
    ```

??? question "Chromosome `X` not found in reference / chromosome mismatch"

    gbcms looks up reference sequence by the chromosome name in the variant file. Ensure
    chromosome naming is consistent across VCF/MAF, BAM, and FASTA:

    ```bash
    # Check names in FASTA
    grep "^>" reference.fa | head -5
    # Check names in BAM
    samtools view -H sample.bam | grep "^@SQ" | head -5
    # Check names in MAF/VCF
    head -20 variants.maf | cut -f5   # Chromosome column
    ```

    Common mismatch: `chr1` vs `1`. The FASTA and variant file must match.

??? question "`FETCH_FAILED` in `validation_status` for all variants on a chromosome"

    See chromosome mismatch above. Also check:
    - FASTA `.fai` index is present (`samtools faidx reference.fa`)
    - No truncated / corrupted FASTA

---

## Allele Count Issues

??? question "`alt_count = 0` for a variant I can see reads for in IGV"

    Check in this order:

    1. **`validation_status`** column — if `REF_MISMATCH` or `FETCH_FAILED`, the variant
       was excluded from counting. See [Normalization Issues](#normalization-issues).

    2. **Complex Del+SNV routing** — if your variant has deletion format (`REF` longer than `ALT`,
       `ALT` is a single base) but the anchor base also changes (e.g., `GC→T`, `AG→T`), it is a
       **complex Del+SNV** routed to `check_complex`, not `check_deletion`.
       Use `--trace` logging to confirm:
       ```bash
       RUST_LOG=trace gbcms dna --variants variants.maf --bam sample.bam \
           --fasta ref.fa --output-dir /tmp/debug/ 2>&1 | grep "check_complex\|check_deletion"
       ```
       See [Complex Indels](../reference/complex-indels.md) for detailed case studies.

    3. **Left-alignment shift** — for large deletions (≥5bp), BWA left-aligns the anchor;
       the CIGAR `D` position may be several bases right of the anchor.
       Run normalization to see the left-aligned coordinates:
       ```bash
       gbcms normalize --variants variants.maf --fasta ref.fa --output-dir /tmp/norm/
       # Check norm_pos vs original Start_Position
       ```

    4. **BAM coverage** — verify the BAM actually has reads at this position:
       ```bash
       samtools view -c sample.bam "chr17:7579309-7579322"
       ```

    5. **Verbose per-read trace** — most diagnostic tool:
       ```bash
       RUST_LOG=trace gbcms dna ... 2>&1 | grep "7579309"  # Filter to specific position
       ```

??? question "`ref_count = 0` for a large deletion"

    For large deletions (~100bp+), REF reads have clean M-only CIGARs. gbcms's
    `is_worth_realignment()` skips Phase 3 for clean CIGARs, and the **M-block REF fallback**
    in `check_complex` classifies these as REF.

    If `ref = 0`, check:
    - The variant was normalized correctly (M-block REF fallback only applies in `check_complex`)
    - The BAM slice has coverage at the anchor position (`samtools depth -a sample.bam -r chr22:30038094-30038095`)
    - The deletion is ≥50bp (interior REF guard applies; smaller deletions use the normal path)

    See [NF2 Case Study](../reference/complex-indels.md#case-2-nf2-large-deletion-ref-reads-invisible).

??? question "`alt_count` is lower than expected for a deletion ≥5bp"

    BWA left-alignment can shift the anchor further left than where the CIGAR `D` appears.
    The `has_nearby_length_match` Phase 3 fallback handles this for deletions ≥5bp.

    For deletions <5bp failing S3 sequence validation, CIGAR-definitive REF is used — this is
    intentional. A 1-4bp deletion in the wrong reference context is almost certainly spurious noise,
    not a left-alignment artifact.

    To diagnose deletion ≥5bp still returning low alt:
    ```bash
    RUST_LOG=trace gbcms dna ... 2>&1 | grep "has_nearby_length_match\|S3.*fail"
    ```

??? question "Large Δalt vs sign-out `t_alt` for MNP / ONP / DNP variants"

    Several expected causes — not gbcms bugs:

    | Cause | Pattern | Diagnosis |
    |:------|:--------|:----------|
    | **TRANS configuration** | Two MNP alleles on opposite chromosomes; sign-out counts 2× reads | VAF ~2× expected |
    | **MAF annotation artifact** | ONP defined at wrong position or missing adjacent SNV | `gbcms normalize` shows position shift |
    | **Sign-out duplication** | Variant captured by two MAF entries (DEL+INS pair both counted | Two nearby rows in MAF with same gene |

    Use `--verbose` to log per-variant read counts, or `--trace` for per-read detail.

??? question "`dp = 0` for a variant"

    - **BAM slice doesn't cover this sample**: `samtools view -c sample.bam "chrom:pos-pos"` returns 0
    - **BAM index missing**: see BAI issue above
    - **Chromosome mismatch**: see above
    - **Variant outside BAM region**: the BAM file only covers certain genomic regions

??? question "`DPF` is much lower than `DP` (large gap)"

    The gap `DPF − (RDF + ADF)` does **not** mean data loss — it counts fragments discarded
    due to R1/R2 disagreement within the quality threshold.

    ```
    DPF - (RDF + ADF) = ambiguous fragments (both R1 and R2 covered but base quality within --fragment-qual-threshold)
    ```

    A large gap at a specific locus is a **quality signal** indicating the site has many
    discordant fragments — typical of noisy or error-prone positions. Adjust with
    `--fragment-qual-threshold` (default 10).

??? question "`ref_count + alt_count > total_count` (INVARIANT violation)"

    This should **never** happen. If observed:

    1. Record the exact MAF row (`chrom`, `pos`, `ref`, `alt`, `Tumor_Sample_Barcode`)
    2. Run with `--trace` to capture per-read classification
    3. File a bug report at [github.com/msk-access/gbcms/issues](https://github.com/msk-access/gbcms/issues)
       with the MAF row, BAM details, and trace log

    The invariant `ref_count + alt_count ≤ total_count` is enforced throughout the counting engine.

---

## Normalization Issues

??? question "`REF_MISMATCH` in `validation_status` — what does it mean?"

    The REF allele in the variant file does not match the reference FASTA at the stated position.
    gbcms requires ≥90% base match to auto-correct; below that, the variant is excluded.

    Common causes:
    - Wrong reference genome build (GRCh37 data with GRCh38 reference or vice versa)
    - Chromosome naming mismatch (`chr1` vs `1`)
    - Trailing-base error in MAF annotation (fixed by tolerant validation if ≥90% match)

    Diagnose with:
    ```bash
    gbcms normalize --variants variants.maf --fasta ref.fa --output-dir /tmp/norm/
    grep "REF_MISMATCH" /tmp/norm/normalized.maf | head -5
    ```

??? question "`PASS_WARN_REF_CORRECTED` — is this a problem?"

    The REF allele was ≥90% match — corrected to FASTA sequence and counting proceeded normally.
    The original MAF REF is preserved in the `original_ref` column (when `--show-normalization` is set).

    This is **not** a problem in most cases. Common benign cause: 27bp EGFR exon 19 deletion
    with 26/27 bases correct and 1 trailing base that's a MAF annotation artifact.

    If you want to audit all corrected variants:
    ```bash
    grep "PASS_WARN_REF_CORRECTED" output.maf | wc -l
    ```

??? question "`PASS_WARN_HOMOPOLYMER_DECOMP` — what changed?"

    The variant overlapped a homopolymer, and the corrected allele (e.g., `CCCCCC→CCCCT` instead
    of `CCCCCC→T`) got more ALT support. The corrected allele counts were used.
    `used_decomposed = True` in the output.

    This is **correct behavior** — it means the variant caller collapsed a D(1)+SNV into a single
    complex variant, and gbcms detected and corrected for this.

---

## RNA Mode Issues

??? question "All counts are 0 in RNA mode"

    1. **Strandedness filter with missing gene annotation** — if your MAF lacks a `gene_strand`
       column, all reads fail the strandedness filter (defaulted to on).
       Disable with `--no-strandedness` for unstranded libraries or MAFs without strand annotation.

    2. **MAPQ filter too strict** — default MAPQ=1 with NH:i:1 rescue is correct for STAR.
       If your aligner assigns MAPQ differently, check `--min-mapq`.

    3. **BAM not STAR-aligned** — RNA mode assumes STAR BAMs with `NH` tag. Other aligners
       may lack the `NH:i:1` tag needed for MAPQ rescue. Check:
       ```bash
       samtools view sample.bam | head -1 | tr '\t' '\n' | grep "NH"
       ```

??? question "`rna_sense_depth` and `rna_antisense_depth` are both 0"

    The `gene_strand` column is missing from the MAF or all variants have `gene_strand = NA`.
    gbcms cannot assign sense/antisense without strand annotation.

    Solutions:
    - Add strand annotation to the MAF (from gene model GFF/GTF)
    - Use `--no-strandedness` to count all reads regardless of strand

??? question "RNA editing flag (`rna_editing_site_overlap`) is always False"

    Ensure `--rna-editing-db` points to a valid REDIportal TABLE1 file:
    ```bash
    # Verify file exists and is readable
    head -3 TABLE1_hg38.txt.gz | zcat | head -3
    
    # Download if missing
    wget http://srv00.recas.ba.infn.it/atlas/download/TABLE1_hg38.txt.gz
    ```

    Also verify chromosome naming matches (REDIportal uses `chr`-prefixed names).

---

## Nextflow Issues

??? question "Container not found: `Unable to find image 'ghcr.io/msk-access/gbcms:X.Y.Z'`"

    Pull manually:
    ```bash
    docker pull ghcr.io/msk-access/gbcms:X.Y.Z
    # or for Singularity
    singularity pull docker://ghcr.io/msk-access/gbcms:X.Y.Z
    ```

??? question "Task metrics error: `Command 'ps' required by nextflow`"

    Use a newer container image with `procps` installed. This occurs with minimal base images
    that lack the `ps` command needed for Nextflow resource tracking.

??? question "Docker permission denied"

    ```bash
    sudo usermod -aG docker $USER && newgrp docker
    ```

---

## Debugging Tools

gbcms has two levels of diagnostic output:

| Flag | Output | Use When |
|:-----|:-------|:---------|
| `--verbose` | Per-variant summary: read counts, classification decisions | First-pass diagnosis |
| `--trace` | Per-read detail: CIGAR walk, S1/S2/S3 checks, Phase 2/3 scores | Deep debugging a specific variant |

```bash
# Verbose: see per-variant decisions
RUST_LOG=info gbcms dna --verbose --variants variants.maf --bam sample.bam \
    --fasta ref.fa --output-dir /tmp/debug/ 2>&1 | tee debug.log

# Trace: per-read detail (slow — run on small BAMs or specific regions)
RUST_LOG=trace gbcms dna --trace --variants variants.maf --bam sample.bam \
    --fasta ref.fa --output-dir /tmp/debug/ 2>&1 | grep "7579309" > tp53_trace.log
```

### Inspecting Genomic Bin Construction

gbcms groups variants into **~10kb genomic bins** before BAM traversal to reduce `bam.fetch()` calls (see [Architecture → Genomic Binning](../reference/architecture.md#genomic-binning)). Enable `RUST_LOG=info` to observe bin construction:

```bash
RUST_LOG=info gbcms dna --variants variants.maf --bam sample.bam \
    --fasta ref.fa --output-dir /tmp/out/ 2>&1 | grep -i "bin"
# Example:
# [INFO  gbcms] Built 42 genomic bins from 1247 variants (window=10000bp, max_per_bin=200)
# [INFO  gbcms] Bin chr17:7571720-7579309 → 83 variants (1 fetch)
```

!!! tip "Expected bin counts"
    - **Targeted panels** (MSK-ACCESS, IMPACT): many variants per gene → very few bins per chromosome — optimal I/O
    - **Sparse WGS** variants: ~1 variant per bin ≈ per-variant I/O (expected, not a bug)
    - Unexpectedly many bins for a clustered MAF? Check for chromosome name mismatches (`chr1` vs `1`) — the binner treats different chromosome names as separate bins regardless of coordinate proximity.

---

## Related

- [Installation](../getting-started/installation.md) — Setup guide
- [Variant Normalization](../reference/variant-normalization.md) — Validation status reference
- [Allele Classification](../reference/allele-classification.md) — How variant types are classified
- [Complex Indels](../reference/complex-indels.md) — Complex Del+SNV, large DEL, shifted deletion case studies
- [Nextflow Parameters](../nextflow/parameters.md) — Resource configuration
