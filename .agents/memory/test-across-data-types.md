---
name: test-across-data-types
description: Validate a counting rule on every data type and input type it reaches, with more than one truth source; GIAB alone is not acceptance, and DNA and RNA may need different rules.
metadata:
  type: feedback
---

A count-changing rule is accepted on every data type it reaches (DNA tissue panel
tumour + matched normal, cfDNA duplex/simplex, WES/WGS without realignment, RNA) and
input type (aligner and realigner, paired/unpaired, clipped/split/hard-clipped
records), against more than one truth source. GIAB HG002 is one genome, one aligner
(BWA-MEM), germline hets only, and its 50+ bp deletions are mostly tandem repeats:
it cannot show false ALT (no hom-REF rows) or how ABRA2, STAR or cfDNA reads behave.

**Why:** operator, 2026-10-12, during C40: "we are relying too much on giab, we should
do testing across the board" and "DNA and RNA can have different rules where
necessary and we should take into account both and different input type". C40's
k, read-on and placement choices had been driven by HG002 plus panel sign-out (not
truth), with the matched-normal arm skipped; a placement candidate was found on HG002
alone.

**How to apply:** each decision table names its data types and truth sources and
says which are missing; normals (hom-REF) are run every time; where no truth exists,
build it (spike-ins). A rule may differ between DNA and RNA when their reads differ
(splices, clips holding the next exon); say so in the rule. Related:
[[survey-several-tools]], [[better-than-standard]], [[genotyper-not-caller]].
