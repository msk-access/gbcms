---
name: tolerant-deletion-deliberate
description: "The tolerant large-deletion match (skip seq check) was a deliberate sensitivity fix; don't naively revert it — add a partial sequence check instead."
metadata: 
  node_type: memory
  type: project
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
---

In `rust/src/counting/variant_checks.rs`, a large deletion whose observed length
differs from expected but passes a ≥50% reciprocal-overlap rule is accepted as ALT
with the deleted-sequence (S3) check skipped (`del_ok = true // tolerant match`).

**Why it matters:** This is a false-positive risk (review **CR-4**) — an unrelated SV
in a repeat region can be counted as ALT. BUT git history (commit `de7da8e4`, and the
v4.x "Interior REF guard removed" CHANGELOG entry) shows the tolerant path was added
*deliberately* to fix large-deletion ALT reads being misclassified/dropped (false
REF / lost sensitivity). Naively reinstating the strict check would regress that fix.

**How to apply:** Keep tolerant length acceptance; add a *partial* sequence-concordance
check over the overlapping span instead of skipping verification entirely. Preserve the
v4.x sensitivity fixtures as guardrails when changing this.
