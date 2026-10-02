---
name: survey-several-tools
description: A community-practice survey covers several callers, genotypers and BAM tools, not just GATK.
metadata:
  type: feedback
---

When a decision needs community practice (how other tools treat a read shape, a
quality, a flag), survey several tools with sources, not one. The set to check:
GATK (HaplotypeCaller/Mutect2), samtools/bcftools mpileup (htslib), fgbio, VarDict,
Strelka2, freebayes, bam-readcount, LoFreq, and the original C++ GetBaseCounts.
Say plainly where a tool's handling could not be found.

**Why:** the operator (2026-10-02) noted, after a group-1 survey that cited only
GATK: "there are more tools than GATK to look at". One tool's choice is a data
point, not the field's practice; the validation standard asks for callers and
genotypers.

**How to apply:** in every decision table's evidence column, list what each
surveyed tool does (or "not documented"), with links. Related: [[genotyper-not-caller]].
