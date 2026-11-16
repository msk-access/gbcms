---
name: references
description: "Outside work to remember: Lancet2 (what to borrow for the read judge, what not); outside Rust tools for the judge kernel, SV and long-read work (POASTA, A*PA, spoars, kanpig, HiPhase, rustybam, bwa-mem3: ideas and comparators, not dependencies); the downstream mskcc-omics-workflows repos (not a release step; never propose the switch, the operator starts it)."
metadata:
  type: reference
---

## Downstream org repos (checked 2026-10-05)

gbcms has two downstream repos in the mskcc-omics-workflows org:
- **containers** (`github.com/mskcc-omics-workflows/containers`): one Dockerfile per version,
  `containers/gbcms/<X.Y.Z>/Dockerfile`. It clones the `msk-access/gbcms` tag, runs `maturin build
  --release --no-default-features` (gbcms has no crate features, so that flag changes nothing), and
  installs the wheel without our lock. It publishes `ghcr.io/mskcc-omics-workflows/gbcms:<X.Y.Z>`,
  which is not our `ghcr.io/msk-access/gbcms` image.
- **modules** (`github.com/mskcc-omics-workflows/modules`, `develop`):
  `modules/msk/gbcmsrs/{dna,rna,normalize,merge,buildgtfcache}` are nf-core-style wrappers. Files are
  staged as `path` inputs, and they pin the org image (6.3.1). `modules/msk/gbcms` is the old C++
  GetBaseCountsMultiSample 1.2.5.

**Rule:** do not touch these repos as a release step, and never propose bumping them. The operator
decided (2026-10-07, at the 6.6.0 cut) that bumping is premature. The reason: our own pipeline
(`nextflow/`) does not use them. It includes only `modules/local/gbcms/{dna,rna,merge,normalize,
filter_maf,pipeline_summary}` (no modules.json), on our image
`ghcr.io/msk-access/gbcms:${workflow.manifest.version}`. Moving our pipeline onto the org modules
comes first, and it is its own project. The operator starts that whole process (pipeline switch,
image, modules) when they are satisfied the org route is fit for production use. When that happens:
first the containers Dockerfile for the version, then the gbcmsrs modules adapted to that version's
CLI (since 6.3.1: drop `buildgtfcache` and the `rna` module's `gtf_cache` input, deprecated in 6.6.0
and removed in 6.7.0; widen the `gtf` pattern to `*.{gtf,gtf.gz}`; check the mFSD, observations and
merge changes). The modules are buehlere's (modules#259, 2026-08-27): ask whether to open a PR or
coordinate first.

## Lancet2 (read 2026-10-24)

NYGC, NAR Genom Bioinform 2026, lqag036: https://github.com/nygenome/Lancet2, docs
nygenome.github.io/Lancet2/latest, paper https://academic.oup.com/nargab/article/8/2/lqag036/8625920.
The operator asked whether the judge's graph direction is "like Lancet2, but for SV too".

**What it is:** a somatic **caller**. About 1 kb windows (20% overlap); a colored bidirected de Bruijn
graph of tumour and normal reads, k from 13 to 127 in steps of 6, rebuilt at a larger k on cycles or
a tangle (the "complexity guard": no haplotypes if no k is clean); paths enumerated (1M BFS visits
cap); SPOA merges haplotypes into a POA graph whose bubbles are the variants. **Genotyping is
linear:** every read is aligned to every haplotype independently with minimap2 (end bonus, large
z-drop and bandwidth; insertions to about 2 kb, deletions to about 500 bp); the best match wins (by
gap-compressed identity and score; one read supports a variant once, by exact sequence match). No
margin, no tie rule. The latest docs add a Dirichlet-multinomial PL and CMLOD; the paper describes
counts only. An Explainable Boosting Machine filter (trained on HCC1395) gives its precision. No
force-genotyping of known alleles. Small variants only: "broader applicability to larger structural
events will require future work"; no results above 50 bp in the paper.

**For gbcms:** our judge shares its architecture (haplotypes, then linear per-haplotype alignment) but
enumerates the given rows instead of assembling, hence our context-enumeration faults where Lancet2
hits tangled windows. The direction discussed: an edge-based graph (reference segments; bubbles =
small variants, edges = splices and SV breakend junctions), judged per adjacency so SNV-to-SV and RNA
share one model (`universal-haplotype-judge`, #322 to #327). **Worth measuring from Lancet2:** its
dual-affine gap model (a second, cheaper extension past about 20 bp: long events), alignment-free
allele discovery from a local cDBG (against CIGAR-spelled found alleles), a complexity guard.
**Not for us:** best-match without a margin, and a learned somatic filter (a genotyper reports every
given row; `genotyper-not-caller`).

## Outside tools for the judge's kernel, the SV and the long-read work (read 2026-11-15)

The operator shared these while S4 waited on review ("I think we can learn some from this"). None is a
dependency; each is an idea or a comparator. Why no aligner drops in as the judge's scorer: they give the
one best alignment to a whole graph (MSA, consensus, mapping), and the judge needs a best score for each
site state, with ties; none masks a low-quality base (invariant 2); a new scorer is a count decision, not
an exact speed-up (`speed-floor`).

**For the kernel (the plan's "later" items, `universal-haplotype-judge`):**
- **POASTA** (Broad, Rust, BSD-3): optimal affine partial-order alignment by A* with exact-match
  extension. Idea: exact pruning; reads with best score 0 take 85 to 87% of the judge's scoring time (M3),
  and a cheapest-first search finds such a path in about one pass while still returning the optimum.
- **A\*PA / A\*PA2** (Groot Koerkamp, Rust, MPL-2.0): the same pruning, proved for pairs; edit distance
  only, global only. A* pays off on long sequences, less on a 150 bp read against a few hundred bases.
- **spoars** (Rust copy of spoa 4.1.5, MIT): a SIMD engine bit-identical to its scalar one; read it for an
  exact integer SIMD kernel over a graph (`run_from`).
- **abPOA** (C, MIT): an adaptive band; not exact unless a bound proves the optimum lies inside.

**For the SV and long-read work (#322 to #327):**
- **kanpig** (Rust, MIT; the Truvari author): long-read SV genotyping on a local graph of every SV in a
  "neighbourhood" (like our context bubbles), matched by size and sequence similarity (0.90), so a long
  read's noisy copy of an insertion still matches. Not for us: it reads SV sequence from the CIGAR pileup
  (no split or soft-clipped reads, events over about 10 kb skipped: our breakend judge, #324, must cover
  them) and calls diploid genotypes. Its `AD` is a comparator on HG002 long reads (#326).
- **HiPhase** (PacBio; research use only; license not read): assigns reads to alleles by a graph WFA
  over small variants, SVs and repeats together, the same move as our 2026-10-25 graph decision; WFA
  suits HiFi's low error. Phasing is out (no phase credit). Ideas from the paper only.
- **rustybam** (Rust, MIT): CIGAR-aware PAF tools (`liftover`, `break-paf`, `trim-paf`). Harness use: an
  assembly-based SV truth for HG002 (minimap2, split at large indels, lift over the breakends) that no
  read aligner shares (#326).
- **trgt** (PacBio): tandem-repeat genotyping from HiFi; a comparator for repeat rows.
- **breakinator**: flags foldback and chimeric long-read artifacts; such reads must not count as
  junction support in #324.
- **CycSim**: a long-read simulator, for testing a long-read error model.
- **perbase** (Rust): fast per-base counts; an SNV comparator beside bam-readcount in `tools/`.
- **fgumi** (Fulcrum, Rust): UMI tools; if ACCESS consensus ever moves from fgbio to it, recheck the
  consensus base qualities and `--min-baseq` (`msk-assay-facts`).

**Upstream aligners:** **bwa-mem3** (Fulcrum; v0.17.0; Bioconda): `--compat=bwa-mem` writes records
byte-identical to bwa 0.7.19 and is several times faster; production IMPACT and ACCESS used 0.7.17, so
check 0.7.17 against 0.7.19 before using it for the harness's spike-in arm. Its Rust bindings
(bwa-mem3-rs) need a prebuilt genome index: no engine use. An aligner change in a pipeline is the
operator's call; measure gbcms on its BAMs before any placement-dependent rule change.

**Not useful:** scirs2-graph (network analysis, no sequence graphs; the judge has its own `Graph` in
`rust/src/counting/judge/graph.rs`; `petgraph` if a general container is ever needed); rs-vgaligner
(remaps reads to a pangenome graph: a comparator, like vg and PanGenie).
