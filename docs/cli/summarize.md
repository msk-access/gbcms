# gbcms summarize

Summarize the QC flags of finished outputs: which flags fired, on how many rows, which rows need a look, and what to do about each.

## Synopsis

```bash
gbcms summarize [OPTIONS] <OUTPUT>...
```

## Description

`gbcms dna` and `gbcms rna` end their log with a QC summary ([QC Flags → End-of-run summary](../reference/qc-flags.md#end-of-run-summary)). `gbcms summarize` makes the same summary from outputs that are already written: an older run, many samples together, or a merged MAF. The text format is the block the run logs; the other formats are the same summary as JSON, a cohort table, or MultiQC content.

It reads each output once and changes nothing in it. The exit code is `0` whether or not rows need a look.

An input that is refused (not a `.maf` or `.vcf`), cannot be read, is not a gbcms output, or shares a sample name with an earlier input is **not checked**: the others are summarized, the block names it under "Not checked" with the reason ("1 of 3 samples not checked"), the closing line is a `WARNING`, and the exit code is `1`. A MultiQC report carries the same warning. With the text format that is all on standard output; the other formats print `not checked: <sample>: <reason>` on standard error.

Every flag, and what to do about it, is defined on the [QC Flags](../reference/qc-flags.md) page.

## Arguments

| Argument | Description |
|:---------|:------------|
| `OUTPUT...` | One or more finished outputs: `.maf` or `.vcf` from `gbcms dna`, `gbcms rna` or `gbcms merge`. Compressed files are not read. |

A sample's name is its file name without `.maf` or `.vcf` (a run writes `<sample><suffix>.maf`). Two outputs may not share a name.

## Options

| Option | Default | Description |
|:-------|:--------|:------------|
| `--format` | `text` | `text`, `json`, `tsv` or `multiqc` (see below) |
| `--output-dir`, `-o` | see below | Where `json`, `tsv` and `multiqc` files are written; created if missing |
| `--verbose`, `-V` | `false` | Enable debug logging |

## Formats

### text

The block the run logs, for one sample or the cohort block for several, then its closing line. It goes to standard output.

```bash
gbcms summarize results/*.maf
```

### json

One `<sample>.qc_summary.json` per output, beside the output or in `--output-dir`. Each file holds the sample's whole summary:

| Field | Content |
|:------|:--------|
| `schema` | `gbcms.qc_summary/1` |
| `sample`, `format`, `merged` | The sample, `maf` or `vcf`, and whether it is a merged MAF |
| `written_by`, `summarized_by` | The build that wrote the output, and the build that summarized it |
| `rows`, `need_look`, `need_look_share` | Rows, rows that need a look, and their share |
| `rows_by_tier` | Rows with a flag of each tier (`fix_input`, `check_call`, `informational`) |
| `flags` | Every flag, those that did not fire too (a zero is evidence): `name`, `family`, `tier`, `rows`, `share`, `action` and a `docs` link |
| `unlisted_flags` | Flags this build does not know, with their rows |
| `rows_to_look_at` | Every row that needs a look: `variant`, `gene` and its `flags` (input problems first) |

```bash
gbcms summarize --format json results/*.maf
```

### tsv

One cohort table, a row per sample: the rows, the rows that need a look and their share, the rows by tier, then a column per flag (zeros included) and the flags this build does not know. It goes to standard output, or to `gbcms_qc_summary.tsv` in `--output-dir`.

```bash
gbcms summarize --format tsv results/*.maf > cohort_qc.tsv
```

### multiqc

[MultiQC custom content](https://docs.seqera.io/multiqc/custom_content), written to `--output-dir` (the current directory by default):

- `gbcms_qc_summary_mqc.json`: a table with a row per sample. The columns are the rows, the rows that need a look, and the rows by tier, coloured by tier (a first-tier row red, a second-tier row amber, none green). A column per flag that fired is there too, hidden until switched on.
- `gbcms_qc_warning_mqc.json`: a "WARNING: gbcms rows need a look" section listing the samples that have such rows. It is written only when some do; a file left by an earlier run is removed when none do now.

MultiQC puts a custom-content section that its config's `custom_content: order` does not name at the top of the report. To pin the warning first, list `gbcms_qc_warning` before the other ids.

```bash
gbcms summarize --format multiqc --output-dir qc/ results/*.maf
multiqc qc/ -o report/
```

## Merged outputs

A merged MAF carries each input's flags under that input's prefix (`duplex_`, `simplex_`), and the summary reads them all: a row's flags are the union of its inputs'. It also derives [`COMBINED_NA`](../reference/qc-flags.md#merged-output) from the combined `simplex_duplex_*` cells, a flag that exists only in the summary.

## Nextflow

The pipeline runs this command over every output with `--qc_summary` ([parameters](../nextflow/parameters.md#qc-summary-options)).

## Related

- [QC Flags](../reference/qc-flags.md) — every flag, and the end-of-run summary
- [gbcms merge](merge.md) — merged MAFs
- [gbcms dna](dna.md) / [gbcms rna](rna.md) — the runs that log the summary
