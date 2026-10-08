---
name: bin-anchor-coverage
description: "A genomic bin's fetch-end must cover the anchor variant's full ref span, not just bin_start + window."
metadata: 
  node_type: memory
  type: project
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
---

In `build_genomic_bins` (`rust/src/counting/engine.rs`), `bin_end` is seeded at
`max(bin_start + window, the anchor's span end)`; the CR-1 bug seeded it at
`bin_start + window` and applied the span extension only to later variants.

**Why it matters:** A bin anchored by a large deletion/DelIns whose
`ref_allele.len()` exceeds `BIN_WINDOW` (10kb) under-fetches its right tail. Reads
aligned past the anchor are never cached; counts change only where such reads
count (clip-admitted carriers of a long event: every other counted read overlaps the
event's first base, which every fetch holds). This was review finding **CR-1**.

**How to apply:** Any binning change must enforce
`bin.end >= max over all bin variants of (v.pos + v.ref_allele.len() + window_pad_v)`,
including the anchor. The `build_genomic_bins` property test and
`tests/test_binning_invariance.py` (clip carriers past a long anchor) pin it; with the
bug reintroduced both fail (mutation-checked 2026-10-01). See [[engine-output-aware]] for the broader
"don't silently drop work" theme.
