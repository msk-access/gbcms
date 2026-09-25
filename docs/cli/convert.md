# gbcms convert

Convert variants between VCF and MAF without counting reads.

## Synopsis

```bash
gbcms convert --variants <VCF> --output <MAF>
gbcms convert --variants <MAF> --fasta <FILE> --output <VCF>
```

## Description

The `convert` subcommand writes the representation `gbcms dna` / `gbcms rna`
use for their own output, without a BAM:

- **VCF → MAF** — each ALT allele as [vcf2maf](https://github.com/mskcc/vcf2maf)
  writes it: leading bases REF and ALT share are trimmed (never trailing ones),
  `Start_Position` / `End_Position` and `Variant_Type` follow the trimmed alleles,
  and the record itself is kept in `vcf_pos` / `vcf_ref` / `vcf_alt`. The columns
  are those of a VCF-input `--format maf` row before its count columns. ALT
  alleles that cannot be counted (`*`, symbolic, breakends, `.`) are skipped the
  same way, with a WARNING per reason.
- **MAF → VCF** — each row as maf2vcf writes it: the reference base is prepended
  when an allele is `-`, or when the alleles differ in length and first base.
  Records are written in input order; contigs keep the MAF's naming.

It converts representation only: REF is not checked against the reference and
nothing is left-aligned ([gbcms normalize](normalize.md) does both). The rules
and worked examples are in [Output Formats](../reference/output-formats.md).

## Required Arguments

| Option | Description |
|:-------|:------------|
| `--variants`, `-v` | VCF (`.vcf`, `.vcf.gz`, `.vcf.bgz`) or MAF (`.maf`) |
| `--output`, `-o` | Output path: `.maf` for VCF input, `.vcf` for MAF input |

## Optional Arguments

| Option | Default | Description |
|:-------|:--------|:------------|
| `--fasta`, `-f` | — | Reference FASTA (indexed). Required for MAF input: maf2vcf's anchor base comes from it. A base it cannot supply is written as `N` and counted in a WARNING. |
| `--verbose`, `-V` | `false` | Enable debug logging |
| `--trace`, `-T` | `false` | Enable per-read Rust trace logging (slow). Implies `--verbose`. |

## Examples

```bash
# VCF -> MAF (vcf2maf coordinates)
gbcms convert --variants calls.vcf.gz --output calls.maf

# MAF -> VCF (maf2vcf records)
gbcms convert --variants mutations.maf --fasta reference.fa --output mutations.vcf
```

## Related

- [gbcms normalize](normalize.md) — REF validation and left-alignment, no counting
- [Output Formats](../reference/output-formats.md) — the VCF and MAF output columns
- [Input Formats](../reference/input-formats.md) — VCF/MAF input conventions
