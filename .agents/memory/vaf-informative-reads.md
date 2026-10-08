---
name: vaf-informative-reads
description: gbcms VAF is AD/(RD+AD) over informative reads; per-sample total_count (DP) includes reads that cannot be judged — never advise AD/DP, least of all for indels in repeats.
metadata:
  type: feedback
---

gbcms's `vaf` column is `AD / (RD + AD)`: REF and ALT count only reads whose own bases
show that allele. The per-sample `total_count` is DP, every read overlapping the
variant, including the neither/uninformative ones (in a repeat, the reads that end
inside the tract and can't tell REF from the indel). So `alt_count / total_count` mixes
an informative numerator with an all-reads denominator and is biased low. Merge's
`simplex_duplex_total_count` is ref + alt (a different meaning of "total").

**Why:** on 2026-10-08 I told the operator to compute VAF as ALT/depth for the AR 9 bp
GGC-repeat deletion (15/191 = 7.9%). The operator questioned it. The informative VAF
is 15/31 ≈ 48% in the tumor and 0% in the matched normal (0/24 crossing reads): a
high-fraction somatic event the clinical sign-out (169 REF, 6.6%) understated about
sevenfold. The docs already said `VAF = AD/(RD+AD)` (docs/index.md,
docs/reference/architecture.md, allele-classification.md "25% instead of 50%").

**How to apply:** for any VAF statement, read gbcms's `vaf` (or AD/(RD+AD)), never
AD/DP or the passthrough `t_var_freq`/`t_depth` (input sign-out values). In repeats,
also say how many reads were informative (RD+AD vs DP) and that spanning bias
inflates a deletion's fraction. Check the documented definition before advising.
Related: [[pysam-validation-oracle]], [[genotyper-not-caller]].
