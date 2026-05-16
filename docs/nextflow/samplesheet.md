# Samplesheet Format

The samplesheet is a CSV file defining input samples for the Nextflow pipeline.

## Basic Format

```csv
sample,bam,bai
sample1,/path/to/sample1.bam,/path/to/sample1.bam.bai
sample2,/path/to/sample2.bam,
sample3,/path/to/sample3.bam,/path/to/sample3.bam.bai
```

## Columns

| Column | Required | Description |
|:-------|:---------|:------------|
| `sample` | Yes | Sample identifier (used in output filenames) |
| `bam` | Yes | Path to BAM or CRAM file |
| `bai` | No | Path to index file. Auto-discovers `.bam.bai`/`.bai` for BAM or `.cram.crai`/`.crai` for CRAM. |
| `suffix` | No | Per-sample output suffix |
| `bam_type` | No | BAM type label for [multi-BAM merging](#multi-bam-type-merging) (e.g., `duplex`, `simplex`) |
| `tsb` | No | Tumor_Sample_Barcode pattern(s) for [MAF filtering](#multi-sample-maf-filtering) |

## Index Auto-Discovery

If the `bai` column is empty, the pipeline auto-discovers the index based on the alignment file extension:

**BAM files** (`.bam`):

1. `<bam>.bai` (e.g., `sample.bam.bai`)
2. `<bam_without_extension>.bai` (e.g., `sample.bai`)

**CRAM files** (`.cram`, v5.3.0+):

1. `<cram>.crai` (e.g., `sample.cram.crai`)
2. `<cram_without_extension>.crai` (e.g., `sample.crai`)

!!! tip
    Leave the `bai` column empty to use auto-discovery. You can also provide an
    explicit index path of any naming convention.

!!! info "CRAM Support (v5.3.0+)"
    The `bam` column accepts both BAM and CRAM file paths.  The `--fasta` reference
    is automatically threaded to the Rust engine and Python layers for CRAM decoding.
    No additional configuration is needed.

## Per-Sample Suffix

For samples with multiple BAM types:

```csv
sample,bam,bai,suffix
sample1,/path/to/sample1.duplex.bam,,-duplex
sample1,/path/to/sample1.simplex.bam,,-simplex
sample1,/path/to/sample1.unfiltered.bam,,-unfiltered
```

Output files: `sample1-duplex.maf`, `sample1-simplex.maf`, `sample1-unfiltered.maf`

## Multi-Sample MAF Filtering

When using `--filter_by_sample` with a multi-sample MAF as `--variants`, each sample's variants are filtered by `Tumor_Sample_Barcode`.

### Exact Match (default)

When no `tsb` column is provided, the `sample` column is matched **exactly** against `Tumor_Sample_Barcode`:

```csv
sample,bam
P-0012345-T01-IM7,/path/to/tumor.bam
P-0067890-T01-IM7,/path/to/tumor2.bam
```

### Regex Match

Use the `tsb` column for pattern matching (e.g., patient-level filtering):

```csv
sample,bam,tsb
patient_A,/path/to/A.bam,P-0012345
patient_B,/path/to/B.bam,P-0067890
```

`P-0012345` matches any `Tumor_Sample_Barcode` containing that substring (e.g., `P-0012345-T01-IM7`, `P-0012345-T02-IM7`).

### Comma-Separated Multi-Select

Select specific samples with comma-separated patterns:

```csv
sample,bam,tsb
tumor_pair,/path/to/A.bam,"P-0012345-T01-IM7,P-0012345-T02-IM7"
```

Duplicate rows (matched by multiple patterns) are automatically deduplicated.

!!! note
    Samples with 0 matching variants are skipped and reported in `pipeline_summary.tsv`.

## Multi-BAM Type Merging

When `--merge_counts` is enabled, samples with a `bam_type` column are grouped by sample ID
and merged into a single output MAF with type-prefixed count columns.

```csv
sample,bam,bai,suffix,bam_type
sample1,/path/to/sample1.duplex.bam,,,duplex
sample1,/path/to/sample1.simplex.bam,,,simplex
```

When `bam_type` is set:

1. The `--column-prefix` is **auto-derived** from the type label (e.g., `duplex_`)
2. The output **suffix is auto-derived** as `-{bam_type}` (e.g., `sample1-duplex.maf`)
   unless an explicit `suffix` column value is provided
3. After per-BAM genotyping, the pipeline groups MAFs by sample and runs `gbcms merge`
4. Combined `simplex_duplex_*` columns are computed when both types are present

!!! tip "Suffix auto-derivation"
    When `bam_type` is set, you do **not** need to set `suffix` — it is automatically
    derived as `-{bam_type}`. Only set `suffix` explicitly if you need a custom value
    that differs from the `bam_type` label.

!!! tip "Minimum 2 BAM types per sample"
    Merge requires at least 2 inputs per sample. Samples with only 1 BAM type
    skip the merge step and output the single MAF directly.

See [Merge CLI](../cli/merge.md) for details on the merge algorithm.

## Related

- [Parameters](parameters.md) — All pipeline options
- [Examples](examples.md) — Common usage patterns
- [Troubleshooting](../resources/troubleshooting.md) — Common issues
