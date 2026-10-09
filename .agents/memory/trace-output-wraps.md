---
name: trace-output-wraps
description: "gbcms --trace read-call lines: set COLUMNS=3000 (Rich wraps them); a read call is not a count (uncounted reads are logged too); the key is the prepared variant, so match rows by haplotype."
metadata:
  node_type: memory
  type: reference
  originSessionId: 547cc6f6-ff25-4d74-a30c-2551fc0fe47a
  modified: 2026-09-29T18:30:57.992Z
---

`gbcms dna|rna --trace` logs through a Rich console handler that wraps each line at
the terminal width, so a subprocess capture splits every `read call CHR:POS REF>ALT
read=… mate=… ref=… alt=… phase=… partial=N …` line and a regex over it silently
matches nothing (every read looks "not traced"). Run it with `COLUMNS=3000` in the
environment. `partial=` is a count (positions matching ALT), not a bool; `mate` is
1 for first-in-template, else 2.

**How to apply:** any per-read adjudication harness that diffs two builds' trace
calls ([[pysam-validation-oracle]], [[genotyper-not-caller]]) — set the env var and
assert a non-zero number of parsed calls before trusting a diff.

**Two more traps (2026-10-09, 6.7.0 plan measurements):**
- **A read call is not a count.** Every classified read is logged, including reads the
  engine then leaves out (a read that misses the variant position and no clip admits).
  At long delins the trace showed REF 2,403 where the MAF wrote 1,356, and a claimed
  "REF counted at both junctions" bias was only that. Take totals from the MAF; use
  the trace for which reads, not how many.
- **The key is the prepared variant** (VCF-anchored, left-aligned), not the input row:
  matching rows by position missed left-aligned deletions (841 REF reads at 100+ bp
  deletions looked untraced). Match a trace key to a row by the haplotype both give.
