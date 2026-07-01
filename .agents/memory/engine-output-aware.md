---
name: engine-output-aware
description: "Prefer plumbing intent (mfsd, output format, gene_strand) across the FFI boundary over computing-then-discarding or running silent no-ops in Rust."
metadata: 
  node_type: memory
  type: project
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
---

Several gbcms bugs share one root: the Rust engine doesn't know what the caller asked
for. `count_bam_binned` has no `mfsd` parameter, so mFSD stats + per-fragment size
arrays are built for *every* variant on *every* run (PF-1). `gene_strand` is never
populated, so RNA `enforce_strandedness` is a silent no-op (HI-7). Sub/mono-nucleosomal
fields are computed then dropped from VCF (ME-1).

**Why it matters:** Compute-then-discard wastes CPU/memory and silent no-ops produce
wrong-but-quiet results. Making the engine output-aware closes a whole class of these
at once.

**How to apply:** When adding engine behavior, thread the *intent* (a flag, the output
format, the strand) across the PyO3 boundary and gate the work on it, rather than
always computing and letting Python ignore it. See [[bin-anchor-coverage]] for the
sibling "don't silently drop reads" case.
