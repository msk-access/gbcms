---
name: siblings-break-binned-legacy-parity
description: sibling_variants intentionally diverge binned from legacy; parity tests must omit them
metadata: 
  node_type: memory
  type: project
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
---

Passing `sibling_variants` to `count_bam_binned` makes it diverge from legacy
`count_bam` **by design** — pangenomic sibling disambiguation is a binned-only
feature that legacy lacks. A read matching a "REF-here + ALT-sibling" haplotype is
classified "neither" by the binned path (with siblings) but REF by legacy.

**Why:** the load-bearing binned↔legacy parity invariant only holds *without*
siblings. The `count_both` test helper calls legacy with no siblings, so passing
siblings to it asserts an impossible parity and fails (e.g. a large-deletion REF
read: legacy rd=1, binned rd=0).

**How to apply:** parity tests (`count_both`) must NOT pass `sibling_variants`. To
exercise multi-variant binning in a parity test, rely on proximity-based bin
grouping (variants near each other share a bin regardless of siblings) — see
`tests/test_parity_large_deletion.py::test_parity_large_deletion_bin_anchor`. Test
sibling disambiguation separately, against expected values, not against legacy.

Related: [[pysam-validation-oracle]], [[bin-anchor-coverage]].
