---
name: references
description: "Outside work to remember: Lancet2 (what to borrow for the read judge, what not) and the downstream mskcc-omics-workflows repos (the containers image and the gbcmsrs modules: not a release step; never propose the switch, the operator starts it)."
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
