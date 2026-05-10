# gbcms normalize

Left-align and validate variants without counting reads.

## Synopsis

```bash
gbcms normalize --variants <FILE> --fasta <FILE> --output <FILE>
```

## Description

The `normalize` subcommand applies the same variant preparation pipeline used
by `gbcms dna` — MAF anchor resolution, REF validation, and bcftools-style
left-alignment — but **without** performing any BAM counting.  The output is a
TSV file showing both the original and normalized coordinates for every variant.

This is useful for:

- **Debugging** — see exactly how each variant was transformed before counting
- **QC** — verify which variants fail REF validation
- **Preprocessing** — normalize a variant list before passing it to other tools

## Required Arguments

| Option | Description |
|:-------|:------------|
| `--variants` | VCF or MAF file with variant positions |
| `--fasta` | Reference FASTA file (with .fai index) |
| `--output` | Output TSV file path |

## Optional Arguments

| Option | Default | Description |
|:-------|:--------|:------------|
| `--threads`, `-t` | `1` | Number of threads |
| `--verbose`, `-V` | `false` | Enable debug logging |
| `--trace`, `-T` | `false` | Enable per-read Rust trace logging (slow). Implies `--verbose`. |

## Output Columns

| Column | Description |
|:-------|:------------|
| `chrom` | Chromosome |
| `original_pos` | Original 1-based position |
| `original_ref` | Original REF allele |
| `original_alt` | Original ALT allele |
| `norm_pos` | Left-aligned 1-based position |
| `norm_ref` | Left-aligned REF allele |
| `norm_alt` | Left-aligned ALT allele |
| `variant_type` | SNP, INSERTION, DELETION, or COMPLEX |
| `gbcms_status` | `PASS`, `PASS;WARN_HOMOPOLYMER_DECOMP`, `REF_MISMATCH`, or `FETCH_FAILED` |
| `was_anchor_resolved` | Whether MAF anchor resolution changed pos/ref/alt |
| `was_left_aligned` | Whether left-alignment shifted the variant |
| `was_normalized` | Combined: `True` if either anchor resolution or left-alignment changed the variant |

## Example

```bash
gbcms normalize \
    --variants mutations.maf \
    --fasta reference.fa \
    --output normalized.tsv \
    --threads 4
```

## Related

- [gbcms dna](dna.md) — Full counting pipeline (with `--show-normalization` flag)
- [Input Formats](../reference/input-formats.md) — VCF/MAF coordinate conventions
