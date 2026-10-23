---
name: production-aligner-settings
description: "Where the production alignment settings live (nucleo for ACCESS/IMPACT standard BAMs, forte for STAR) and the decided stand-ins for the D12 spike-in arm; full flags are on #275."
metadata:
  node_type: memory
  type: reference
  originSessionId: a8529157-ff34-4400-8321-c28f0357971d
  modified: 2026-10-19T17:20:11.522Z
---

Production alignment settings, read from the public repos on 2026-10-19 (full flag
tables in the #275 comment of that date):

- **ACCESS** — `msk-access/nucleo` (CWL; tools in submodules `uncollapsed_bam_generation`,
  `bam_collapsing`, `cwl-commandlinetools`, `cwl_subworkflows`; read at master
  815a5794). Standard BAM: fastp 0.20.1 → `bwa mem -Y` (0.7.17; -K/-T unset in the CWL,
  the docs claim K=1000000 T=30) → MergeBamAlignment (Picard defaults: trimmed bases
  come back as soft clips) → MarkDuplicates → ABRA2 2.22 (`--mad 1000 --sc 100,30,80,15
  --mmr 0.1 --sga 8,32,48,1 --ca 10,1 --ws 800,700 --cons --no-edge-ci --nosort`) →
  FixMate → BQSR. ABRA2 targets: dynamic (`bedtools genomecov -bg | bedtools merge -d 10`
  of the BAM) for collapsed BAMs, and for the standard BAM until
  uncollapsed_bam_generation@8e270032 (2024-11-07) switched it to a static BED input
  (`UBG_abra2_targets`, required; release 1.3.2 / nucleo 3.1.2 include it). nucleo
  `feature/release_3.1.3` makes it optional (uncollapsed_bam_generation@735a35af): no BED
  → no `--targets` (ABRA2: whole genome eligible), not bedtools targets. Collapsed BAMs: fgbio consensus → `bwa mem -Y` → ABRA2
  again (+`--index`, targets from the BAM's own coverage) → FilterConsensusReads
  (min base quality 30). See [[abra2-collapse-artifacts]], [[abra2-junction-placement]].
- **IMPACT** has its own pipeline; the operator decided (2026-10-19) to use ACCESS's
  standard-BAM path as its stand-in for D12 ("should be similar for standard BAMs").
  Production BWA sets `-K 1000000 -Y` (operator, 2026-10-19); `-T` is unset (the docs'
  T=30 is bwa's default).
- **FORTE** — `mskcc/forte` (Nextflow; main 0b3d9275): fastp 0.23.2, STAR 2.7.10a
  `--twopassMode Basic --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1
  --outSAMattributes NH HI AS NM MD --outSAMstrandField intronMotif`, UMI dedup only (no
  MarkDuplicates). Config default genome is GRCh37 (Ensembl 75), but production runs
  GRCh38 (operator, 2026-10-19), as our FORTE arm does (Ensembl 111 GTF). FORTE's fill-out step
  calls the legacy C++ GetBaseCountsMultiSample (`gbcms:1.2.5`), not this package.

**How to apply:** for D12 (#275) re-align spike-ins with these tools and flags; say
which values are CWL defaults vs docs where they disagree.
