---
name: rna-mapq-default
description: The RNA --min-mapq default of 1 is deliberate (keeps STAR's 2-4-locus reads); don't switch to unique-only (MAPQ 255) or to 0 without new evidence.
metadata:
  type: project
---

RNA mode keeps reads STAR placed at two to four loci (MAPQ 3 or 1), counted once at
the primary alignment (operator, 2026-10-03). Junction reads tie between a gene
and its processed pseudogene (spliced to the parent, unspliced to the retrocopy),
so they come out MAPQ 3.

**Why:** measured on FORTE (truth 33 samples, probes 3): unique-only cost 12 real
ALT reads in the truth set, 10 of them at PIK3CA E545K (exon 9 has a chr22
pseudogene copy), and 1.6% of junction fragments at the probes; MAPQ 0 added 6 ALT
and 33 REF reads, consistent with the samples, but admits 5+-locus reads. ASE
pipelines (GTEx, phASER) use unique-only plus WASP; that suits SNP allelic
imbalance, not a genotyper. STAR's MAPQ depends only on the number of loci.

**How to apply:** treat the default as decided; `--min-mapq 0` stays a deliberate
per-run choice for pseudogene-family genes (as PMS2 in DNA, [[mapq0-loci-pms2]]).
Measurement in `docs/reference/read-filters.md`; harness `~/test/gbcms/harness/g3/`.
