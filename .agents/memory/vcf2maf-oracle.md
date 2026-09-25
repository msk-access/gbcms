---
name: vcf2maf-oracle
description: vcf2maf.pl / maf2vcf.pl are the VCF<->MAF representation oracle; how to run them locally and their quirks
metadata:
  type: reference
---

gbcms's VCF <-> MAF conversion (`CoordinateKernel.vcf_to_maf` / `maf_to_vcf` /
`maf_alleles`, `gbcms convert`) is defined as "what mskcc/vcf2maf writes" (#110).
The scripts are checked out at `~/Documents/Github/vcf2maf` (commit 589406f) and
need `samtools` on PATH (homebrew has it).

- VCF -> MAF: `perl vcf2maf.pl --input-vcf X --output-maf Y --ref-fasta F
  --inhibit-vep --tumor-id TUMOR --normal-id NORMAL --vcf-tumor-id TUMOR
  --vcf-normal-id NORMAL` (the VCF needs FORMAT/GT and those two sample columns).
- MAF -> VCF: `perl maf2vcf.pl --input-maf X --output-dir D --output-vcf D/x.vcf
  --ref-fasta F`.

Quirks (don't mistake them for gbcms bugs):
- maf2vcf dies on any allele outside `[ACGT-]` (N, lowercase), and puts rows
  whose REF mismatches the FASTA into `D/<name>.skipped.tsv`.
- maf2vcf writes an empty ALT (invalid VCF) for a REF == ALT row; gbcms FAILs
  it as `ALT_EQUALS_REF`.
- vcf2maf trims case-sensitively; gbcms case-insensitively (documented).
- maf2vcf writes a differing `Tumor_Seq_Allele1` as a second ALT; gbcms doesn't.

Acceptance record (2026-09-24): 102k unique real sign-out shapes matched both
oracles exactly except the two REF == ALT rows. The harnesses are local
(`~/test/gbcms/harness/t110/`) because the inputs are PHI; see
[[user-local-validation-data]]. Related oracle: [[pysam-validation-oracle]].
