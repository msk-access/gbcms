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
Reference/SplitReference/Variant/SplitVariant/ReadCount); SVtorm (jblancoheredia/SVtorm,
read 2026-10-24) calls Delly, GRIDSS, Manta and SvABA on tumour/normal pairs; SURVIVOR merges
(1 kb, 2 callers, no type or strand agreement); the "recall" is a second GRIDSS 2.13.2 run
(operator: SVtorm relies on GRIDSS for recalling) on fragments overlapping the panel target
BED (`params.intervals`; the merged candidates only steer its metrics), then a PON somatic
filter (SVSOMF); SURVIVOR filters the five VCFs; CrosSV reconciles them against the tumour
BAM (slop 100, per-patient prior BED) and appends recovered SVs. Output `*_CrosSV.tsv` /
`*_ANNOTE_SV_INN.tsv`, columns (SVRVOR2TSV.py): SV_BC, CHROM, POS, ID, REF, ALT, QUAL,
FILTER, CIEND, CIPOS, CHR2, END, SVLEN, SVMETHOD, SVTYPE, SUPP_VEC, SUPP, STRANDS, PE, SR,
CT, DP (site 2 = CHR2/END, orientation = STRANDS or CT); annotation is iAnnotateSV
(rhshah/iAnnotateSV, operator 2026-10-24): input `chr1 pos1 str1 chr2 pos2 str2` (1-based,
str 0 = plus, 1 = minus; chromosomes without chr), output adds gene/transcript/site 1 and 2 and
fusion; SVtorm's `*_ANNOTE_SV_INN.tsv` is its input, iCallSV uses it, and IMPACT data_sv
site descriptions use its vocabulary, so its six columns are the proposed internal
breakend-pair form (#323); a standard-VCF reader must take
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

**Ticketed (2026-10-24; both milestones moved to 6.9.0 on 2026-10-26, see the release plan below):** the SV work is parent #322 with
sub-issues #323 (input formats and the separate output), #324 (breakend-junction judge),
#325 (fragment evidence), #326 (truth, validation, small-variant overlap, community
standard) and #327 (RNA fusions on mRNA junctions). Decided (operator, 2026-10-24): when
the judge becomes the default, `--read-judge rules` stays one release as the fallback;
removing the previous classifier is #328 (moves with the flip). R5–R9 (#313–#317)
stay in 6.7.0: an R-ticket goes later only if every item is separate from 6.7.0 work.

**Scorer architecture under evaluation (2026-10-24; plan and pre-registered acceptance on
#304):** most judge faults came from enumerating contexts x states, so a graph scorer is being
prototyped on local branch `proto/c47-graph-scorer` (worktree c47 = 971f8f01 + the c46e2
fix): E1 per-site-state DP over a local sequence graph (context rows as bubbles), E3 coloured
k-mer triage (phase 1); E2 Lancet2's two-tier gap model and E4 k-mer bubble allele discovery
(phase 2, after the found-allele floor investigation). Graph = segments + edges so splices and
SV junctions (#322) are edges later; read-length independent; pluggable error model for long
reads. Coloured de Bruijn graphs are triage and discovery only (k-mers collapse repeat copies
and count k-mers, not reads). Pangenome stacks (vg, PanGenie) are comparators, not the
architecture: they remap reads, use population priors and report diploid genotypes
([[references]]).

**Decided (operator, 2026-10-25):** the graph is the judge's scorer (E1 failed pre-registered
criterion 1 only through the undecided "fits nothing" rule, applied exactly, and through
every-combination; it equals enumeration where enumeration is exact, is ~3x faster, removes
539 lines and 14 caps). No per-site switching: the graph is exact where enumeration is, and
switching to the rules per site would bring back a seam; speed comes from a cheap exit for
reads that fit nothing. k-mer triage (E3) dropped (exact variant saves 2-18%). found2 floor
confirmed on patient arms (no change at depth). Patient-data faults to fix on the graph:
rows inside a long deletion's bases collapse REF; MNV ALT misses; homopolymer and 20-49 bp
REF losses; gate cost 11.8x develop. Phases A-E with bars on #304 (A integrate, B fits-nothing
decision, C fix round, D MNV/SNV policy decision, E speed).

**Decided (operator, 2026-10-25), fits nothing:** xs + diagnostic: a read that fits no state
is depth only where the aligner also doubts it (XS >= AS); every row counts its fits-nothing
reads in gbcms_diagnostic. The absolute budget mostly caught structural reads (20+ bp clips,
10+ bp indels) the given rows lack, not noise, and cost panel SNV/MNV ~1 point; residual risk
is uniquely misplaced foreign reads at homozygous deletions in repeats (chr17 3.6 kb: 0.59).
MNV/SNV (O5 + MNV_SHARED diagnostic) decided the same day.

**Framing and name (operator, 2026-10-25):** the judge resolves which local sequence a read's
own bases spell at the row (not phasing, no linkage credit); contexts only explain the
read's other differences. User-facing name `--read-judge sequence`, "candidate sequences" in
docs (rename after phase C). Whether each context source (listed rows, guard background,
found alleles) earns its place is being measured, pre-registered on #304.

**Release plan (operator, 2026-10-26):**
- **6.7.0** ships the judge opt-in (#304).
- **6.8.0** only makes it the default, "meticulously and systematically": parent #331, with sub-issues
  - #332 tests and every output path,
  - #333 RNA,
  - #334 acceptance and the `--min-baseq` sweep,
  - #335 speed,
  - #336 the HPC gate,
  - #337 the switch-over.
- **6.9.0** takes the SV/fusion work (#322–#327) and retiring the rules judge (#328).
- The integration round (`src_c53d`) is pushed as branch `feature/c46-haplotype-judge`, with no PR yet. Its
  re-score on the corrected checks is recorded on #304.

**Decided (operator, 2026-10-26):**
- **The length-only insertion stays the ALT** (RJ-20 as amended on 2026-10-15, #278). The judge already follows it; the "port" proposed that morning came from stale docs.
- **Ties stay depth only.** This is GATK's uninformative read; develop's placement tie-break leans to REF.
- **xs stays.** Scaled's 162 extra REF reads were reads the engine never counts.
- **HG002 bars** (truth v2d): score only inside the benchmark regions; score long deletions on the reads the engine counts; compare homozygous rows with REF.
- **GATK's read evaluation** (the PairHMM, shared with its `--alleles` mode) is the fair comparator for read rules. Its discovery step, and its decision whether a variant is present, are not.

**How to apply:** package 2 (#304 widened) is this judge. #246's prototype (c45c,
worktree c38, uncommitted; patch + README in harness plan670/parked_246) was parked for
this (operator, 2026-10-23); its red test and the HG002 #246 pair arm (giab/c246) are
acceptance cases. R1 survey verdict: per-haplotype alignment is the scorer (graph DP
can't beat it at <=16 haplotypes), k-mer triage in front, margin = one confident base at
--min-baseq, discovered alleles on the "lacking" side (REJ-20261007-002). R2 bake-off
approved (S0 rules, S1 per-haplotype PairHMM, S2 quality-masked edit cost, S3 k-mer
triage + S1), standalone Rust prototype in the harness; HG002 SV and GIAB TR arms approved. Phasing tools are oracles' sources, not the method:
population priors may never credit reads. Related: [[genotyper-not-caller]],
[[msk-assay-facts]], [[working-agreements]]; the decision method is `.agents/rules/decision-protocol.md`.
