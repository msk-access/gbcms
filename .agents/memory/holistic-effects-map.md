---
name: holistic-effects-map
description: Every decision or fix gets an effects map — each place its change lands (all counting paths and outputs) — before implementation
metadata:
  type: feedback
---

Before implementing a decision or fix, map every place its effect lands and
state it in the PR (and the plan ticket):
- both counting paths (binned and the legacy parity oracle);
- per-transcript counts, ASJD, mFSD, and the observations export;
- MNP rescue, tract-cluster claiming, merge;
- the MAF and VCF writers and the diagnostics/flags;
- docs, tests and Nextflow.

For each, say "changes" or "unchanged, because …".

**Why:** the operator asked (2026-09-25) that decisions be made holistically,
"thinking about all the places our decisions have effects". A rule changed in
one path but not its twins, or in the counts but not a flag that reads them, is
the recurring failure mode. See [[legacy-parity-oracle]] and [[engine-output-aware]].

**How to apply:** trace the fields the change touches (e.g. `is_ref`/`is_alt`,
`has_nearby_evidence`, `qual`) to every consumer (grep them), then verify the
"unchanged" claims on real data as well as the "changes" ones.
