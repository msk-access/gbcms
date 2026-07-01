---
name: bin-anchor-coverage
description: "A genomic bin's fetch-end must cover the anchor variant's full ref span, not just bin_start + window."
metadata: 
  node_type: memory
  type: project
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
---

In `build_genomic_bins` (`rust/src/counting/engine.rs`), `bin_end` is seeded as
`bin_start + window` from the anchor (leftmost) variant, and the `var_end + window/2`
extension is applied only to *subsequent* variants in the inner loop — never the anchor.

**Why it matters:** A bin anchored by a large deletion/DelIns whose
`ref_allele.len()` exceeds `BIN_WINDOW` (10kb) under-fetches its right tail. Reads
spanning the right breakpoint are never cached, so AD/ADF are undercounted and the
binned path diverges from legacy `count_bam` (which fetches per-variant). This is
review finding **CR-1** (CRITICAL).

**How to apply:** Any binning change must enforce
`bin.end >= max over all bin variants of (v.pos + v.ref_allele.len() + window_pad_v)`,
including the anchor, and must be covered by a binned↔legacy parity test that includes
a >10kb deletion and a complex DelIns. See [[engine-output-aware]] for the broader
"don't silently drop work" theme.
