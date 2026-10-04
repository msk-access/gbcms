---
name: survey-several-tools
description: A community-practice survey covers several callers, genotypers and BAM tools (not just GATK) and the published literature, with sources.
metadata:
  type: feedback
---

When a decision needs community practice (how other tools treat a read shape, a
quality, a flag), survey several tools with sources, not one. The set to check:
GATK (HaplotypeCaller/Mutect2), samtools/bcftools mpileup (htslib), fgbio, VarDict,
Strelka2, freebayes, bam-readcount, LoFreq, and the original C++ GetBaseCounts.
Say plainly where a tool's handling could not be found. Also survey the published
literature (method and benchmarking papers, with DOI/PMID) for the topic, and add
domain tools where the topic is domain-specific (e.g. for RNA: GATK RNA best
practices and SplitNCigarReads, STAR, ASE counters such as ASEReadCounter, phASER
and WASP, RNA variant callers such as RNA-MuTect and SNPiR, editing callers such
as REDItools and JACUSA).

**Why:** the operator (2026-10-02) noted, after a group-1 survey that cited only
GATK: "there are more tools than GATK to look at". One tool's choice is a data
point, not the field's practice; the validation standard asks for callers and
genotypers.

The operator (2026-10-03) added, before group 3 (RNA): "we should also look at
literature and other tools of what people have done".

**How to apply:** in every decision table's evidence column, list what each
surveyed tool does (or "not documented"), with links. Related: [[genotyper-not-caller]].
