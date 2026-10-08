---
name: abra2-junction-placement
description: "IMPACT/ACCESS BAMs are ABRA2-realigned; a junction insert there carries assembly evidence, so WES (BWA only) is the control for placement questions"
metadata:
  node_type: memory
  type: project
  originSessionId: 547cc6f6-ff25-4d74-a30c-2551fc0fe47a
  modified: 2026-10-07T14:39:53.897Z
---

MSK IMPACT and ACCESS BAMs are ABRA2-realigned (94% of RC junction-ALT reads carry its
YO tag = the original BWA alignment). ABRA2 moves a read onto an assembled contig only
with strictly fewer mismatches, and a singleton haplotype gets no contig, so a read it
writes at the variant's junction with a flank substitution is the given ALT plus a
single-molecule error (C37 #245: 28/28 RC reads, spectrum-checked). The WES loci
(`harness/c20/wes`, BWA + MarkDuplicates) are the no-realigner control: use them for any
"where did the aligner put the gap" question, as the operator suggested 2026-10-07.

**Why:** C37 almost changed the strict path on the premise that a realigner forced reads
onto the known indel; the data showed ABRA2's placement was the better-informed one.
**How to apply:** before a placement-dependent rule change, measure RC (ABRA2), WES (BWA)
and FORTE (STAR) separately and check recurrence with a local haplotype spectrum
(`harness/c37/spectrum.py`). Related: [[genotyper-not-caller]], [[count-the-given-allele]].
