---
name: trace-output-wraps
description: "gbcms --trace read-call lines: set COLUMNS just above the longest allele (Rich wraps lines; 3000 wraps alleles over ~2.9 kb) but never a blanket 60000 (Rich pads every line to it: >450x slower); a read call is not a count; match rows by haplotype."
metadata:
  node_type: memory
  type: reference
  originSessionId: 547cc6f6-ff25-4d74-a30c-2551fc0fe47a
  modified: 2026-09-29T18:30:57.992Z
---

`gbcms dna|rna --trace` logs through a Rich console handler that wraps each line at
the terminal width, so a subprocess capture splits every `read call CHR:POS REF>ALT
read=… mate=… ref=… alt=… phase=… partial=N …` line and a regex over it silently
matches nothing (every read looks "not traced"). Run it with `COLUMNS=60000` in the
environment: `COLUMNS=3000` still wraps lines whose alleles pass ~2.9 kb (it dropped 3,226
calls on the HG002 pri_ins and pri_del50 arms, #265 measurement, 2026-10-16). `partial=` is a count (positions matching ALT), not a bool; `mate` is
1 for first-in-template, else 2.

**The cost of a wide COLUMNS (2026-10-24, the 6.7.0 speed review):** Rich pads every log
line to the full width, so `COLUMNS=60000` makes a trace more than 450x slower than the
untraced run (HG002 pri arm: 4.7 GB in 11 minutes against 1.3 s; a 15k-row SV arm ran over
5 hours and was abandoned). Set COLUMNS to the longest allele in the input plus ~400 (e.g.
12000 for alleles up to 10 kb) and trace only arms of a few thousand rows; the lasting fix
is a Rust-side read-call file instead of log lines (6.7.0 review ticket).

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
