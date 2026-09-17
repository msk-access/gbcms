---
name: tolerant-deletion-deliberate
description: "The ≥50bp tolerant length-acceptance for deletions is superseded by evidence (issue #91): real large pure dels show zero length wobble; wrong-length single-op pure indels are distinct alleles → partial_alt, never full ALT. Delins/split representations stay Phase-3."
metadata: 
  node_type: memory
  type: project
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
  modified: 2026-09-17T19:30:39.399Z
---

**Updated 2026-09-17 (supersedes the 2026-07 guidance "keep tolerant length
acceptance").** Read-level validation for issue #91 measured nine real pure
deletions (12–539bp) across five full clinical BAMs: **every supporting read is
an exact-length single `D` op — zero ±1–3bp wobble, zero splits**. The CR-4
reciprocal-overlap (≥0.5) premise has no observed support at read level, and its
sequence check is vacuous for same-anchor candidates (overlapping deleted bases
are identical reference bases by construction). Below the threshold, Phase-3
grants free ALT anyway (large-del ALT haplotype = short junction window).

**Three-regime taxonomy (all data-validated):**
1. Small pure indels in repeats — coexisting distinct-length populations are
   distinct alleles (e.g. 1bp germline slippage vs 2bp somatic). Exact length
   only; wrong length → `partial_alt`.
2. Delins/complex — split/mismatch-absorbed representations of one event;
   Phase-3 realignment is CORRECT (curated 74-indel set concordant with
   sign-out). Do not touch.
3. Large pure indels — dels are exact even at 539bp (band ≤3bp is insurance
   only); insertions show a truncation smear that is the same event → sequence
   containment (identity ≥~90%, non-low-complexity insert) keeps it.

**How to apply:** the fix (issue #91) replaces reciprocal-overlap with an
exact-length rule (±3bp band ≥50bp) for single-op pure dels, routes wrong-length
single-op pure indels to `partial_alt` (skip Phase-3), keeps Phase-3 for
multi-op/split reads and all delins. Preserve Fix 2/3/4 + PAX5 gains — see the
must-not-regress list in #91. Related: [[siblings-break-binned-legacy-parity]],
[[pysam-validation-oracle]].
