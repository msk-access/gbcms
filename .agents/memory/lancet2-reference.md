---
name: lancet2-reference
description: "Lancet2 (NYGC, NAR Genom Bioinform 2026, lqag036): a somatic CALLER — local colored de Bruijn graph assembly per ~1 kb window, SPOA to extract variants, reads realigned linearly to each haplotype with minimap2; small variants only (SVs named as future work). What to borrow and what not for gbcms's judge."
metadata:
  type: reference
---

Lancet2 (https://github.com/nygenome/Lancet2, docs nygenome.github.io/Lancet2/latest,
paper https://academic.oup.com/nargab/article/8/2/lqag036/8625920), read 2026-10-24 when the
operator asked whether the judge's graph direction is "like Lancet2, but for SV too".

**What it is:** a caller. ~1 kb windows (20% overlap); a colored bidirected de Bruijn graph
of tumour/normal reads, k 13 -> 127 in steps of 6, rebuilt at a larger k on cycles or a
tangle ("complexity guard": no haplotypes if no k is clean); paths enumerated (1M BFS
visits cap); SPOA merges haplotypes into a POA graph whose bubbles are the variants.
**Genotyping is linear:** every read aligned to every haplotype independently with minimap2
(end bonus, large z-drop and bandwidth; insertions to ~2 kb, deletions to ~500 bp), best
match wins (by gap-compressed identity and score; one read supports a variant once, by
exact sequence match). No margin, no tie rule. Latest docs add a Dirichlet-multinomial PL
and CMLOD; the paper describes counts only. An Explainable Boosting Machine filter (trained
on HCC1395) gives its precision. No force-genotyping of known alleles.
**Size:** small variants only; "broader applicability to larger structural events will
require future work"; no results above 50 bp in the paper.

**For gbcms:** our judge shares its architecture (haplotypes, then linear per-haplotype
alignment) but enumerates given rows instead of assembling, hence our context-enumeration
faults where Lancet2 hits tangled windows. The direction discussed: an edge-based graph
(reference segments; bubbles = small variants, edges = splices and SV breakend junctions),
judged per adjacency so SNV-to-SV and RNA share one model ([[universal-haplotype-judge]],
#322-#327). Worth measuring from Lancet2: its dual-affine gap model (a second, cheaper
extension past ~20 bp: long events), alignment-free allele discovery from a local cDBG (vs
CIGAR-spelled found alleles), a complexity guard. Not for us: best-match without a margin
and a learned somatic filter (a genotyper reports every given row; [[genotyper-not-caller]]).
