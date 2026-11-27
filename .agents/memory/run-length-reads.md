---
name: run-length-reads
description: "#334/#333 (reads that end inside a run beside or inside the change): DECIDED 2026-11-26 as RJ-36, A1 narrowed at every site (do not re-measure B, A2, the SNV or MNV arm, or A1 as built); an SNV row that moves only in RD or antisense depth beside a changed cis delins row is the multi-allelic REF guard, a known effect (do not re-trace)."
metadata:
  type: project
---

The rule is RJ-36 in `docs/reference/read-judgment.md` (with its evidence paragraph). The plan, with every run,
build and trace, is `docs/harness/plans/334-run-length-reads.md` (§3; §7 from 2026-11-25).

## Decided: A1 narrowed, at every site (operator, 2026-11-26: "A1n", then "Every site")

- Complex windows grow through every run (1 to 6 base unit, 2+ copies, doublets count) on either allele that
  overlaps or touches the change, to one base past its far end; a read short of it is depth only.
- The run-next rule spares a read whose confident base at Q shows the change, so RJ-23's guard credit stays as
  decided (A1 as built overrode it; that class stays open in the register, "RJ-23's guard credit").
- Every site: the judge holds the rules path's windows also where no run lies. The one existing test this changes,
  approved with the choice: `test_reads_ending_inside_a_delins_count_alike_on_either_allele`, (20, 20) to (10, 10)
  (the rules path already gave (10, 10)). The 4 Rust geometry tests were approved on 2026-11-25.
- Not taken, all measured on every data type: B, A2 (fails the speed floor), the SNV arm (3.6x to 4.7x CPU, FN +68
  on HG002 RNA past A1n), the MNV arm (loses oracle-ALT reads at MNV calls), A1 as built.

**Why:** the options took two days of builds, HPC batches and traces; the numbers are in RJ-36's evidence.
**How to apply:** land A1n's code from worktree `agent-a3f063e0f69b90fd5` (patch
`hpc/stage/inputs/builds_patches/o334a1n.patch` in the harness) on a feature branch from develop, with the test
change above and an outside review measured before any PR. Re-open an option only if the operator asks. RNA VAF
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
