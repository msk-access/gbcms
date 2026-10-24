---
name: universal-haplotype-judge
description: Package 2 direction (2026-10-23): one read-to-haplotype judge for every event size, SNV to SV, graph-based candidates (colored de Bruijn graph, path-aware graph alignment) against a PairHMM baseline; replaces size-specific rules.
metadata:
  type: project
---

Operator, 2026-10-23: stop patching haplotype resolution shape by shape; find one
universal judge that works across all sizes, SVs included ("All sizes"). A read
counts for a row's ALT only where the best haplotype carrying that allele (alone, in
cis with other listed rows, with the locus's background changes) beats the best
lacking it by a margin; else depth only. Discovered alleles explain reads, never
get credit (invariant 7).

The scorer is chosen by research, not assumed: R1 survey (Fable agents) of colored
de Bruijn graphs (Rust: `debruijn`, GGCAT API, cuttlefish-rs), path-aware graph
alignment (Rust: POASTA, RecGraph, rspoa, handlegraph) and probabilistic linear
models (GATK/Mutect2, Octopus, STR models, WhatsHap/HapCUT2 realignment); R2 offline
bake-off judged against HG002 phased truth (never a BWA oracle the engine would
share); R3 decision with worked examples, then integration behind a flag.

**Why:** each size-specific shortcut (#294 windows, #300 sibling by bases, #246 listed
indels, #304/#305 compound hets, #290 imperfect repeats) broke at its edge; the
operator also wants SV-scale events ABRA2 writes into the small-variant files judged
by their bases (graph judges ignore the CIGAR; ABRA2 changes placement, not bases).

SV scope (operator, 2026-10-23): SVs come in a SEPARATE input file and get a separate
output file (the small-variant MAF unchanged); formats: cBioPortal SV text (`data_sv.txt`
style, local dump `~/Downloads/msk_impact/`), a standard SV VCF (e.g. SVtorm's; symbolic
alleles, BND, CIPOS/CIEND), and iCallSV text. Plan: normalize all to breakpoint pairs and
genotype each breakend as a junction (REF adjacency vs novel junction), translocations
included; output counts named to map onto data_sv's Tumor_Split/Paired_End/Read/Variant
counts. Phasing: (1) judge for sequence-resolved small-variant rows, (2) fragment evidence,
(3) SV file input. Small-variant, SV and CNV calls are separate at MSK but overlap. Formats found: iCallSV's
final columns are in rhshah/iCallSV `mergeFinalFiles.py` (Chr1/Pos1/Chr2/Pos2, SV_Type,
Connection_Type, BrkptType, ConsensusSequence = Delly's junction contig, Tumor/Normal
Reference/SplitReference/Variant/SplitVariant/ReadCount); SVtorm (jblancoheredia/SVtorm)
merges Delly, GRIDSS, Manta, SvABA and a recall via CrosSV into TSV (`*_CrosSV.tsv`,
`*_ANNOTE_SV_INN.tsv`; docs unwritten, example requested); a standard-VCF reader must take
Manta (symbolic + BND, SVINSSEQ, HOMSEQ), Delly (CT, CONSENSUS), GRIDSS/SvABA (BND only).

DECIDED (operator, 2026-10-23): package 2 adopts the judge for small-variant rows — S2
quality-masked edit scorer (read whole, haplotype ends free, masked base 0, mismatch 4,
gap 6+k), haplotypes = contexts (other listed rows within read reach) x site states (REF, X,
other alleles at X's site), ALT/REF only by a one-confident-base margin, ties depth-only,
no phase credit. Built in 6.7.0 behind a flag; default before the freeze only if it passes
every arm (PHI included), cost <= 2x develop, and the HPC gate; else opt-in, flip in 6.8.0.
R2 evidence (HG002 truth, develop -> judge, decidable carriers credited): #246 pairs SNV
91.9 -> 96.7%, INS 1-9 60 -> 97%; TR DEL 10-49 74 -> 98%, INS 50-299 44 -> 96%, INS 300-999
42 -> 96%; SV 50 bp-10 kb 96-99%; false ALT equal or lower except INS 10-49 (20 vs 9). PairHMM
scorer rejected for this (repeat gap model sinks perfect short repeat indels below the
margin). Prototype: harness plan670/r2/hapjudge (Rust), truth.py/evaluate.py beside it.
Build needs: speed (~300x develop now: k-mer triage, SIMD banded alignment e.g.
block-aligner, trimming), mRNA haplotypes for RNA, junction-based admission for long
deletions (RJ-24), discovered site alleles on the lacking side (REJ-20261007-002), guard
background changes as context rows, binning invariance via prep-level row lists.

**How to apply:** package 2 (#304 widened) is this judge. #246's prototype (c45c,
worktree c38, uncommitted; patch + README in harness plan670/parked_246) was parked for
this (operator, 2026-10-23); its red test and the HG002 #246 pair arm (giab/c246) are
acceptance cases. R1 survey verdict: per-haplotype alignment is the scorer (graph DP
can't beat it at <=16 haplotypes), k-mer triage in front, margin = one confident base at
--min-baseq, discovered alleles on the "lacking" side (REJ-20261007-002). R2 bake-off
approved (S0 rules, S1 per-haplotype PairHMM, S2 quality-masked edit cost, S3 k-mer
triage + S1), standalone Rust prototype in the harness; HG002 SV and GIAB TR arms approved. Phasing tools are oracles' sources, not the method:
population priors may never credit reads. Related: [[genotyper-not-caller]],
[[abra2-collapse-artifacts]], [[best-not-least-bad]], [[model-per-task]].
