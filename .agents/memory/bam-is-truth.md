---
name: bam-is-truth
description: "Operator principle — the BAM (reads) is the truth; sign-out/t_alt is one more result to compare against, never the target. Score validations read-by-read."
metadata:
  node_type: memory
  type: feedback
  originSessionId: 890ffab2-2fbd-474b-915c-62b9807e66d3
  modified: 2026-09-24T01:57:54.544Z
---

When gbcms disagrees with sign-out, the question is "what do the reads show",
not "how do we match sign-out". Matching sign-out is evidence only when the
reads agree (e.g. TERT C250T: 93 reads carry the change, sign-out 93); where
they don't (a signed-out MNP no read carries), gbcms reporting the annotation as
absent is correct and sign-out's annotation is the thing that is wrong.

Sign-out alleles themselves can be wrong: curators sometimes move or edit
variants by hand and make mistakes (seen locally as garbage alleles like `LU`,
`SV`, `?N[]{}`, REF mismatches against the reference, and shifted
representations). The truth standard is what IGV shows for the reads. So a
read census must discover the haplotypes the reads carry, not only match the
annotated ones, and a recurrent unannotated haplotype is a finding (a missed
or mis-described allele), not noise — a finding for validation and
diagnostics, never ALT evidence for the row ([[count-the-given-allele]]).

**Why:** stated by the operator 2026-09-23 while reviewing the MNP rescue work
(T7), after rescue matched sign-out by adopting a germline SNP on ACCESS data.

The operator added (2026-09-25) that sign-out alleles are sometimes wrong from
manual edits, and that people use IGV as the truth.

**How to apply:** validate count-changing work against pysam/trace read counts
and matched normals, across every assay the change reaches (IMPACT tumour+normal,
ACCESS duplex+simplex via `gbcms merge`); treat any feature that reports a
different allele under a row's label (e.g. `--rescue-mnp`) as opt-in, flagged
and audited. See [[pysam-validation-oracle]].
