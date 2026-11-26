---
name: run-length-reads
description: "#334/#333 (reads that end inside a run beside or inside the change): A1's form is SETTLED as A1 narrowed (2026-11-26; do not re-measure A1 as built); an SNV row that moves only in RD or antisense depth beside a changed cis delins row is the multi-allelic REF guard, a known effect (do not re-trace)."
metadata:
  type: project
---

The plan is `docs/harness/plans/334-run-length-reads.md` (§3 A1 scope item 6; §7 from 2026-11-25 22:43). The
register's Open list (`docs/reference/read-judgment.md`) holds the class until it is decided.

## Settled: A1's form is A1 narrowed (operator, 2026-11-26: "Ok lets go with A1 narrowed")

The run-next rule fires unless the read's confident base at Q (the first base where "one more unit of the run"
and "the guard's change at the masked base" differ) shows the change, so RJ-23's guard credit stays as decided.
A1 as built (build o334a1, `.so` 1e30c02b) cut more false ALT only by overriding RJ-23 next to runs, a class of its
own. Evidence (HG002 RNA: change against develop; HG002 DNA: counts):

| Data | Measure | A1 as built | A1 narrowed |
|:--|:--|--:|--:|
| HG002 RNA | truth-ALT rows with ALT 0 (FN) | +46 | +18 |
| HG002 RNA | truth-REF rows with ALT 1+ (FP) | -31 | -9 |
| HG002 DNA | base partners, false ALT reads (rules; develop 726) | 399 | 671 |

On patients (pat334) the two forms differ by 4 reads in all.

**Why:** each form was built and measured on every data type (about a day of builds, an HPC batch and a trace).
**How to apply:** build later options (the SNV arm `snvn`, the MNV arm `mnvn`) on A1 narrowed (worktree
`agent-a3f063e0f69b90fd5`; its patch is kept in the harness, `hpc/stage/inputs/builds_patches/o334a1n.patch`). Re-open the form only if
the operator asks. Still open: the choice between the faithful options (A1n, B, the two arms). The two arms on A1n,
HG002 RNA against develop (plan §7 15:08): SNVn FN +86, FP1 -321, 4.67x CPU; MNVn FN +30, FP1 -9, 1.09x CPU. RNA VAF
is read with `d334/bin/rna_vaf334.py` (AD/(RD+AD)), never the scorer's ALT/DP `vaf_med`.

## Known effect: SNV rows move through the multi-allelic REF guard

When an option changes how a co-annotated complex row judges a read, an SNV row whose REF window holds that row's
change gains or loses that read as REF (`sibling_claims_ref`, engine.rs). Only RD and antisense depth move; AD
never does, so FN and FP do not change. Traced on HG002 RNA 2026-11-26 (job 14261698): all 92 changed SNV rows,
every moved read (B 742, A1 and A1n 24), tied to 41 `cis` delins rows in run category "grown".

**Why:** the no-run bar flagged these rows as outside the class, which the stop rule sends to the operator; the
trace cost an HPC job and an afternoon. **How to apply:** before tracing an SNV row that moves, check its fields
(AD unchanged) and whether an in-class sibling's change lies in its REF window: `d334/bin/norun334.py` sorts the
changed rows (92 of 92 traced rows match). Trace only rows it leaves unexplained. Related: [[decided-counting-rules]], [[genotyper-not-caller]].
