---
name: holistic-effects-map
description: Every decision or fix gets an effects map — each place its change lands (all counting paths and outputs) — before implementation
metadata:
  type: feedback
---

Before implementing a decision or fix, map every place its effect lands and
state it in the PR (and the plan ticket):
- both read loops (the main loop and per-transcript counting);
- per-transcript counts, ASJD, mFSD, and the observations export;
- MNP rescue, tract-cluster claiming, merge;
- the MAF and VCF writers and the diagnostics/flags;
- docs, tests and Nextflow.

For each, say "changes" or "unchanged, because …".

**Why:** the operator asked (2026-09-25) that decisions be made holistically,
"thinking about all the places our decisions have effects". A rule changed in
one path but not its twins, or in the counts but not a flag that reads them, is
the recurring failure mode. See [[engine-output-aware]].

**How to apply:** trace the fields the change touches (e.g. `is_ref`/`is_alt`,
`has_nearby_evidence`, `qual`) to every consumer (grep them), then verify the
"unchanged" claims on real data as well as the "changes" ones. A derived output
(`vaf`, a flag) is checked at every place it is computed before saying how a
change moves it: the per-sample writers compute `vaf` as alt / (REF + ALT),
`gbcms merge` as alt / total — a group 3 claim read from one site was wrong.
