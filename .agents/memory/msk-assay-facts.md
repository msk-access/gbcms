---
name: msk-assay-facts
description: "How MSK BAMs are made and what that means for counting: IMPACT and ACCESS are ABRA2-realigned (ACCESS twice), WES is BWA-only, ACCESS XS1 (Marianas) and XS2 (fgbio) have different base-quality scales, production aligner settings (nucleo, FORTE), and no pipeline stage sets the QC-fail flag. Read before a placement-dependent or quality-dependent rule change."
metadata:
  type: project
---

## ABRA2 realignment (IMPACT and ACCESS)

IMPACT and ACCESS BAMs are ABRA2-realigned (94% of the junction-ALT reads in the RC set carry its
`YO` tag, the original BWA alignment; the RC set is the 6.5.0 release-candidate DNA runs, 28 of them,
kept locally under `rc650/`). ABRA2 moves a read onto an assembled contig only with strictly
fewer mismatches, and a singleton haplotype gets no contig, so a read it writes at the variant's
junction with a flank substitution is the given ALT plus a single-molecule error (C37 #245: 28 of 28
RC reads, spectrum-checked). The WES loci (`harness/c20/wes`, BWA plus
MarkDuplicates) are the no-realigner control: use them for any "where did the aligner put the gap"
question (operator, 2026-10-07).

Events of 50 bp to about 500 bp in IMPACT's small-variant files are called by VarDict or MuTect,
made easy by ABRA2's CIGAR signature (explicit long I/D ops on its contigs), not by independent
evidence. IMPACT's SV calls (the cBioPortal SV table) start at about 500 bp, so 50 to 500 bp events
have no second caller (operator, 2026-10-23). Judge such events by bases, read by read, against the
`YO` placement. C37 almost changed the strict path on the premise that a realigner forced reads onto
the known indel; the data showed ABRA2's placement was the better-informed one.
**Before a placement-dependent rule change,** measure the RC set (ABRA2), WES (BWA) and FORTE (STAR)
separately and check recurrence with a local haplotype spectrum (`harness/c37/spectrum.py`).

**ACCESS runs ABRA2 twice:** once on the standard (uncollapsed) BAM and once after UMI collapsing
(duplex and simplex). The post-collapse pass sometimes introduces an event that is not in the
molecules; the operator checks the standard BAM to see whether an event was really there
(2026-10-13). ABRA2 rewrites alignments, not bases: a rule that trusts the CIGAR counts an introduced
event; a rule that judges the read's own bases does not, unless the bases carry it. At ACCESS indel
rows treat collapsed-BAM ALT written only by the second pass as suspect, and when the standard BAM is
available compare its counts with the collapsed ones before accepting a count-changing rule. Such a
check belongs in validation tooling or `gbcms_diagnostic`, not a new column.

## ACCESS XS1 and XS2 base qualities

The operator confirmed (2026-10-26): "XS1 is Marianas and XS2 is fgbio". Measured on the plan670
ACCESS slices (2 XS1 and 3 XS2 samples; 3,000 reads per BAM):

| | XS1 (Marianas) | XS2 (fgbio) |
|---|---|---|
| headers | ABRA2 2.17, no bwa @PG | bwa 0.7.17 + ABRA2 2.22 |
| read tags | none from the collapser | MI, RX, aD/bD/cD, aE/bE/cE, ac/bc, ... |
| simplex quality | Q80 to Q90 (Q90 about 74%) | capped at Q45 (99.5%) |
| duplex quality | Q50 to Q90 (Q90 about 92% tumour, 48% normal) | Q68 / Q79 / Q90 |
| masked positions | N at Q5 (about 1.3% of simplex bases) | N at Q2 (about 0.4 to 0.5%) |

How the scales arise (source survey, 2026-10-26, #304): fgbio gives a single-strand consensus at most
Q45 (its pre-UMI error rate), a duplex base where the strands agree the sum of the two strand
qualities (Q90 = 45 + 45), and N at Q2 on a tie or a no-call; then nucleo's FilterConsensusReads
`--min-base-quality 30` masks every base below Q30 to N, which is why no XS2 base sits between Q3 and
Q29. Marianas sums the strand qualities, clamped to 5 to 90, and scales each strand's quality by its
read count, so a simplex read built from 3 or more reads reaches Q90. In XS1 only thresholds of 80 or
more tell bases apart.
**Why it matters:** a quality threshold above 45 masks an entire XS2 simplex BAM but none of an XS1
simplex BAM, so no single "high" `--min-baseq` means "both strands agree" across generations; the
duplex/simplex split is carried by the BAM files themselves. Thresholds 10 to 30 mask the same thing in
both (only the masked N positions), so ACCESS barely moves across 0 to 30, and a quality-weighted
mismatch cost is close to a no-op there (every base is Q45 or more).
**How to apply:** report ACCESS results split XS1 against XS2; sweep `--min-baseq` as on other arms
but add no ACCESS-only high threshold without asking what it should select; brief agents that ACCESS
qualities depend on the collapser.

## Production aligner settings

Read from the public repos on 2026-10-19 (full flag tables in the #275 comment of that date):
- **ACCESS** (`msk-access/nucleo`, CWL; submodules `uncollapsed_bam_generation`, `bam_collapsing`,
  `cwl-commandlinetools`, `cwl_subworkflows`; master 815a5794). Standard BAM: fastp 0.20.1, then
  `bwa mem -Y` (0.7.17; -K and -T unset in the CWL, the docs claim K=1000000 T=30), MergeBamAlignment
  (Picard defaults: trimmed bases come back as soft clips), MarkDuplicates, ABRA2 2.22 (`--mad 1000
  --sc 100,30,80,15 --mmr 0.1 --sga 8,32,48,1 --ca 10,1 --ws 800,700 --cons --no-edge-ci --nosort`),
  FixMate, BQSR. ABRA2 targets: dynamic (`bedtools genomecov -bg | bedtools merge -d 10` of the BAM)
  for collapsed BAMs, and for the standard BAM until uncollapsed_bam_generation@8e270032 (2024-11-07)
  switched it to a static BED input (`UBG_abra2_targets`, required; release 1.3.2 / nucleo 3.1.2
  include it). nucleo `feature/release_3.1.3` makes it optional (uncollapsed_bam_generation@735a35af):
  no BED, no `--targets` (the whole genome is eligible), not bedtools targets. Collapsed BAMs: fgbio
  consensus, `bwa mem -Y`, ABRA2 again (plus `--index`, targets from the BAM's own coverage), then
  FilterConsensusReads (min base quality 30).
- **IMPACT** has its own pipeline; the operator decided (2026-10-19) to use ACCESS's standard-BAM
  path as its stand-in for D12 ("should be similar for standard BAMs"). Production BWA sets
  `-K 1000000 -Y` (operator, 2026-10-19); `-T` is unset (the docs' T=30 is bwa's default).
- **FORTE** (`mskcc/forte`, Nextflow; main 0b3d9275): fastp 0.23.2, STAR 2.7.10a `--twopassMode
  Basic --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --outSAMattributes NH HI AS NM MD
  --outSAMstrandField intronMotif`, UMI dedup only (no MarkDuplicates). The config default genome is
  GRCh37 (Ensembl 75), but production runs GRCh38 (operator, 2026-10-19), as our FORTE arm does
  (Ensembl 111 GTF). FORTE's fill-out step calls the legacy C++ GetBaseCountsMultiSample
  (`gbcms:1.2.5`), not this package.
- For D12 (#275), re-align spike-ins with these tools and flags, and say which values are CWL
  defaults and which come from the docs where they disagree.

## No QC-fail flag in MSK data

Flagstat sweeps (2026-08-12) across ACCESS/CMO pipeline stages (uncollapsed_FM, collapsed_grouped,
standard `_cl_aln_srt_MD_IR_FX_BR`, collapsed, duplex, simplex; about 1.2B reads, XS1 and XS2
samples) and IMPACT standard BAMs found **zero** reads with FLAG 0x200. So `--filter-qc-failed` and
`--no-filter-qc-failed` change nothing on real MSK data today. The filter itself is correct end to end
(CLI to FFI to the gated check in `rust/src/shared/filters.rs`; no unconditional skip like bug #82),
and `test_qc_failed_filter_contract` in `tests/test_filters.py` (PR #86) pins it: flagged reads are
dropped by default, and `--no-filter-qc-failed` restores byte-identical counts. If an upstream
(Illumina chastity passthrough, fgbio consensus marking) starts setting 0x200, the filter becomes
load-bearing; the contract test already covers it.
