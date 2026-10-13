---
name: abra2-collapse-artifacts
description: "ACCESS runs ABRA2 twice (standard BAM, and again after UMI collapsing); the second pass can introduce indel events, checked against the standard BAM"
metadata:
  node_type: memory
  type: project
  originSessionId: a8529157-ff34-4400-8321-c28f0357971d
  modified: 2026-10-13T16:46:38.377Z
---

ACCESS runs ABRA2 twice: once on the standard (uncollapsed) BAM and once after UMI
collapsing (duplex/simplex). The post-collapse pass sometimes introduces an event that is
not in the molecules; the operator checks the standard BAM to see whether an event was
really there (operator, 2026-10-13).

**Why:** ABRA2 rewrites alignments, not bases. A rule that trusts the CIGAR (an I/D op at
the junction) counts an introduced event; a rule that judges the read's own bases does not,
unless the bases carry it.
**How to apply:** at ACCESS indel rows, treat collapsed-BAM ALT written only by ABRA2's
second pass as suspect; when the standard BAM is available, compare its counts with the
collapsed ones before accepting a count-changing rule. Such a check belongs in validation
tooling or `gbcms_diagnostic`, not a new column. Related: [[abra2-junction-placement]],
[[genotyper-not-caller]], [[test-across-data-types]].
