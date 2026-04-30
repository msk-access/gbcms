# gbcms dna

Count alleles at variant positions across one or more DNA/cfDNA BAM files.

!!! warning "Migrating from `gbcms run`"
    `gbcms run` is deprecated and hidden. Replace with `gbcms dna` — all arguments are identical. `gbcms run` will be removed in v4.1.0.

## Synopsis

```bash
gbcms dna [OPTIONS] --variants <FILE> --bam <NAME:PATH>... --fasta <FILE>
```

## Required Arguments

| Option | Description |
|:-------|:------------|
| `--variants`, `-v` | [VCF or MAF](../reference/input-formats.md) file with variant positions (`.vcf`, `.vcf.gz`, `.vcf.bgz`, or `.maf`). Unsupported extensions are rejected immediately. |
| `--bam`, `-b` | BAM file path (can repeat). Optionally prefix with `name:` for sample naming, e.g. `--bam tumor:tumor.bam`. If no name given, the filename stem is used. |
| `--bam-list`, `-L` | File containing BAM paths (one per line, optionally `sample_name path`). Alternative to repeated `--bam`. |
| `--fasta`, `-f` | Reference FASTA file (with .fai index) |
| `--lenient-bam` | Skip missing `--bam` paths and continue with remaining samples (default: exit immediately on first missing BAM). Note: a missing `--bam-list` file always fails regardless. |

## Output Options

| Option | Default | Description |
|:-------|:--------|:------------|
| `--output-dir`, `-o` | *required* | Output directory |
| `--format` | `vcf` | Output format (`vcf` or `maf`) |
| `--suffix` | `''` | Suffix for output filenames |
| `--column-prefix` | `''` | Prefix for gbcms count columns in MAF output. Only letters, digits, and underscores (`[A-Za-z0-9_]`) are allowed; invalid characters exit immediately. |
| `--preserve-barcode` | `false` | Keep original Tumor_Sample_Barcode from input MAF. No-op (with warning) when input is not MAF. |
| `--show-normalization` | `false` | Append `norm_*` columns showing left-aligned coordinates |
| `--context-padding` | `5` | Minimum flanking bases for haplotype construction. Range **1–50**, enforced at parse time. Auto-increased in repeat regions when `--adaptive-context` is enabled. |
| `--adaptive-context` | `true` | Dynamically increase context padding in [tandem repeat regions](../reference/variant-normalization.md#adaptive-context-padding) |
| `--threads` | `1` | Number of threads |

## mFSD Options

Mutant Fragment Size Distribution (mFSD) analysis compares insert-size distributions
for REF- vs ALT-classified fragments at each variant position, enabling detection of
short-fragment enrichment associated with tumor-derived cfDNA
(see [mFSD Metrics](../reference/counting-metrics.md#mfsd)).

| Option | Default | Description |
|:-------|:--------|:------------|
| `--mfsd` | `false` | Enable mFSD analysis. Adds 34 mFSD columns (KS test, LLR, mean sizes, pairwise comparisons, derived metrics) to MAF output and 7 `MFSD_*` INFO fields to VCF. |
| `--mfsd-parquet` | `false` | Write a companion `<sample>.fsd.parquet` with per-variant raw fragment size arrays (`ref_sizes`, `alt_sizes`). Enables downstream visualizations. **Requires `--mfsd`**. |
| `--mfsd-report` | `false` | Generate an interactive HTML report (`<sample>.mfsd_report.html`) with per-variant fragment size distributions, dual-axis histograms, and Fragment Origin Signal classification. **Implies `--mfsd` and `--mfsd-parquet`** (both are auto-enabled). See [mFSD Report](../reference/mfsd-report.md). |
| `--mfsd-report-min-alt` | `3` | Minimum ALT fragment count to include a variant in the HTML report. |
| `--mfsd-report-max-variants` | `20` | Maximum variants in the HTML report (selected by highest ALT count). Use `-1` for no limit. |

!!! tip
    To generate summary statistics, raw Parquet data, and an interactive HTML report:
    ```bash
    gbcms dna --mfsd-report --format maf \
        --variants variants.maf --bam tumor:tumor.bam --fasta hg19.fa -o ./results
    ```
    This produces `<sample>.maf` (with 34 mFSD columns), `<sample>.fsd.parquet`
    (raw fragment sizes), and `<sample>.mfsd_report.html` (interactive visualization).
    Using `--mfsd-report` auto-enables `--mfsd` and `--mfsd-parquet`.

## Filtering Options

| Option | Default | Description |
|:-------|:--------|:------------|
| `--min-mapq` | `20` | Minimum MAPQ |
| `--min-baseq` | `20` | Minimum BASEQ |
| `--filter-duplicates` | `true` | Filter duplicate reads |
| `--filter-secondary` | `true` | Filter secondary alignments |
| `--filter-supplementary` | `true` | Filter supplementary alignments |
| `--filter-qc-failed` | `true` | Filter QC failed reads |
| `--filter-improper-pair` | `false` | Filter improperly paired reads |
| `--filter-indel` | `false` | Filter reads with indels |
| `--fragment-qual-threshold` | `10` | Quality difference threshold for fragment consensus (see [Fragment Counting](../reference/counting-metrics.md#fragment-counting)) |

## BAQ Options

Base Alignment Quality (BAQ) heuristically downgrades base qualities near indels to prevent systematic errors from realignment artifacts.

| Option | Default | Description |
|:-------|:--------|:------------|
| `--apply-baq/--no-baq` | `off` | Enable BAQ quality downgrade near indels |

!!! info "When to Enable BAQ"
    Most modern pipelines (BQSR, fgbio consensus) already recalibrate base qualities. Enable BAQ only for legacy BAMs lacking quality recalibration, where bases near indels may have inflated quality scores that lead to false-positive allele calls.

## UMI Options

Unique Molecular Identifier (UMI) support for molecule-level deduplication.

| Option | Default | Description |
|:-------|:--------|:------------|
| `--umi-tag` | _(none)_ | BAM tag for UMI barcode (e.g. `RX`). Enables UMI-aware fragment grouping. |

!!! tip "UMI-Aware Fragment Counting"
    When `--umi-tag` is set, two reads are considered the same fragment only if they share both **QNAME** and **UMI barcode**. This prevents UMI-collapsed reads from different original molecules being incorrectly merged into a single fragment, which would deflate fragment-level allele counts.

    ```bash
    gbcms dna --umi-tag RX \
        --variants variants.vcf --bam sample.bam --fasta ref.fa -o results/
    ```

## Debugging Options

| Option | Default | Description |
|:-------|:--------|:------------|
| `--verbose`, `-V` | `false` | Enable verbose debug logging |
| `--trace`, `-T` | `false` | Enable per-read Rust trace logging (slow). Implies `--verbose`. Shows detailed per-read classification diagnostics. |

## Alignment Backend

Phase 3 (haplotype-based) classification uses a **two-stage pipeline**:

1. **WFA fast-path** (Wavefront Alignment, `wfa2lib-rs`) — edit-distance triage against the pangenomic haplotype matrix. Resolves ~70-80% of reads instantly at O(s²) cost where *s* is the edit distance. If REF and ALT scores differ clearly, the read is classified immediately.

2. **Marginalized PairHMM** (escalated only when WFA is ambiguous) — integrates per-base quality probabilities into alignment scoring, producing a log-likelihood ratio (LLR) confidence score. More sensitive in noisy, low-quality, or repeat-dense regions.

| Option | Default | Description |
|:-------|:--------|:------------|
| `--alignment-backend` | `pairhmm` | Phase 3 backend: `pairhmm` (WFA + PairHMM, default) or `sw` (Smith-Waterman only, no WFA triage). Invalid values are rejected at parse time. |
| `--llr-threshold` | `2.3` | PairHMM LLR threshold for confident calls (≈ ln(10)) |
| `--gap-open-prob` | `1e-4` | PairHMM gap-open probability for non-repeat regions |
| `--gap-extend-prob` | `0.1` | PairHMM gap-extend probability for non-repeat regions |
| `--repeat-gap-open-prob` | `1e-2` | PairHMM gap-open probability for tandem repeat regions |
| `--repeat-gap-extend-prob` | `0.5` | PairHMM gap-extend probability for tandem repeat regions |

!!! tip "pairhmm vs sw"
    `pairhmm` (default) uses WFA edit-distance triage first, then escalates to PairHMM only for ambiguous reads. This is faster than running PairHMM on every read and more accurate in low-quality or repeat-dense regions.

    `sw` runs Smith-Waterman on every Phase 3 read (no WFA pre-filter). Use only if you need exact reproducibility with versions <3.0.0.

## Examples

=== "Single BAM"

    ```bash
    gbcms dna \
        --variants mutations.vcf \
        --bam sample:sample.bam \ # (1)!
        --fasta reference.fa \
        --output-dir results/
    ```

    1.  The `sample:` prefix sets the output filename. Without it, the BAM filename stem is used.

=== "Multiple BAMs"

    ```bash
    gbcms dna \
        --variants mutations.maf \
        --bam tumor:tumor.bam \
        --bam normal:normal.bam \
        --fasta reference.fa \
        --format maf # (1)!
    ```

    1. MAF output preserves all input MAF columns and appends gbcms count columns.

=== "With Filtering"

    ```bash
    gbcms dna \
        --variants mutations.vcf \
        --bam sample:sample.bam \
        --fasta reference.fa \
        --filter-duplicates \
        --min-mapq 30
    ```

=== "With Normalization"

    ```bash
    gbcms dna \
        --variants mutations.maf \
        --bam sample:sample.bam \
        --fasta reference.fa \
        --format maf \
        --show-normalization # (1)!
    ```

    1. Appends `norm_chrom`, `norm_pos`, `norm_ref`, `norm_alt` columns showing left-aligned coordinates.

=== "With UMI Tags"

    ```bash
    gbcms dna \
        --variants mutations.vcf \
        --bam sample:umi_deduped.bam \
        --fasta reference.fa \
        --umi-tag RX \ # (1)!
        --output-dir results/
    ```

    1. The `RX` tag is the standard SAM tag for UMI barcodes (fgbio, gencore).

## Output

See [Output Formats](../reference/output-formats.md) for a complete column-level schema reference covering:

- **VCF** output: `##INFO` fields, `##FORMAT` fields, and annotated examples
- **MAF** output: VCF→MAF vs MAF→MAF column sets, `Tumor_Sample_Barcode` behaviour, and column prefix options
- **mFSD columns** (with `--mfsd`): all 34 mFSD fields
- **Normalization columns** (with `--show-normalization`)

## Related

- [Quick Start](../getting-started/quickstart.md) — Common patterns
- [gbcms rna](rna.md) — RNA-seq counting with transcriptome-aware filtering
- [gbcms normalize](normalize.md) — Standalone normalization (no counting)
- [Nextflow Pipeline](../nextflow/index.md) — For many samples
- [Input Formats](../reference/input-formats.md) — VCF/MAF specs
- [Output Formats](../reference/output-formats.md) — Complete column-level output reference
- [Variant Counting](../reference/allele-classification.md) — How each variant type is counted
