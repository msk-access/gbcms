---
name: identity-band-annotation-tolerance
description: "The ≥90% insert-identity tolerance exists for ANNOTATION-side representation errors (imperfect ALT strings), not just sequencing errors — BQ masking cannot substitute for it. Test both policies; never assume BQ-first."
metadata: 
  node_type: memory
  type: feedback
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
  modified: 2026-09-18T00:45:40.954Z
---

Operator correction (2026-09-17, during #91/#92 review): the ≥90% identity
tolerance for insert matching was designed to keep counting ALT reads **when
the ALT representation itself is imperfect** — an insert string off by a base,
an adjacent variant folded into the annotation, an incomplete definition. That
error source produces *position-consistent, full-confidence* mismatches shared
across reads; BQ-aware exact matching would demote every read at such a locus
to partial (`ad=0`). Read-side sequencing errors (scattered, low-BQ, per-read)
are the case BQ masking handles.

**Why:** the two error sources have opposite remedies, and the discriminating
signature (mismatch position-consistency across reads × BQ) is measurable —
so the policy choice is empirical, not aesthetic.

**How to apply:** when touching insert matching (the #92 error-tolerant band,
or `insert_truncation_match`), test BOTH policies (raw identity band vs
BQ-masked) against the divergence classes, run the offline mismatch-anatomy
survey on the local panels first, and consider the hybrid (BQ-mask first, band
on residual confident mismatches). Panels are regression validation for one
assay, never universal calibration — stutter/error rates vary by sample type,
assay, and duplex collapsing ([[tolerant-deletion-deliberate]]). No new CLI
flags or output columns for this (standing rule: new signal goes to logs,
existing diagnostics, or validation tooling).
