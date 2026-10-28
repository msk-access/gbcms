---
name: trace-output-wraps
description: "gbcms --trace read-call lines: run the CLI through plan670/r2/plain_trace.py (plain log handler: no Rich padding or wrapping, identical calls); with plain gbcms, Rich wraps at COLUMNS and pads every line to it (wide COLUMNS is hundreds of times slower); a read call is not a count; match rows by haplotype."
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
line to the full width and renders on one Python thread, so a wide COLUMNS makes a trace
hundreds of times slower than the untraced run (HG002 pri arm at 60000: 4.7 GB in 11 minutes
against 1.3 s; the 15k-row SV arm ran over 5 hours at 60000 and again at 23200, on one core,
and was abandoned both times).

**Use the plain-trace wrapper (2026-10-24):** `~/test/gbcms/harness/plan670/r2/plain_trace.py`
runs the gbcms CLI with Rich's handler swapped for a plain `%(message)s` StreamHandler
(`<build venv>/bin/python plain_trace.py dna ... --trace`; `r2/run_plain.sh ARM LABEL PYTHON`
wraps it for the R2 arms). No COLUMNS needed, nothing wraps, nothing pads. Checked on the
GIAB TR arm (build c46c): the same MAF and the same 97,328 read calls as the Rich trace (the
Rich lines differ only by padding and the source-path column), in 1:59 against 6:19 at
COLUMNS=12000, at 380% CPU instead of ~225%. The lasting fix is still a Rust-side read-call
file instead of log lines (6.7.0 review ticket).

**How to apply:** any per-read adjudication harness that diffs two builds' trace
calls ([[pysam-validation-oracle]], [[genotyper-not-caller]]) — use the wrapper (or set
COLUMNS) and assert a non-zero number of parsed calls before trusting a diff.

**Two more traps (2026-10-09, 6.7.0 plan measurements):**
- **A read call is not a count.** Every classified read is logged, including reads the
  engine then leaves out (a read that misses the variant position and no clip admits).
  At long delins the trace showed REF 2,403 where the MAF wrote 1,356, and a claimed
  "REF counted at both junctions" bias was only that. Take totals from the MAF; use
  the trace for which reads, not how many.
- **The key is the prepared variant** (VCF-anchored, left-aligned), not the input row:
  matching rows by position missed left-aligned deletions (841 REF reads at 100+ bp
  deletions looked untraced). Match a trace key to a row by the haplotype both give.
