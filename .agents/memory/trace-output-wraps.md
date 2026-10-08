---
name: trace-output-wraps
description: "gbcms --trace lines are wrapped by the Rich console handler; set COLUMNS=3000 before parsing \"read call\" lines."
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
