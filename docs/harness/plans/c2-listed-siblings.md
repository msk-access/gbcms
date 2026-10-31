# C2: is a listed allele of the ALT's length always a candidate?

Status: draft, for the operator's approval
Rounds: 0 of 2 · New conditions: at most 1 (option B) · HPC batches: 0 of 4 (two builds, each with HG002 and patient arrays)

## 1. The class

- **Reads:** reads that fit the ALT (X) and a listed allele of X's length (a co-annotated sibling, S)
  equally: they do not read a confident base where X and S differ.
- **Evidence:** today two decided rules answer this class differently.
  - **RJ-25 (2026-10-14):** S is another allele only where at least 1 fragment and 5% of those
    nearest it spell S in full. Below that, the tie reads count ALT. The recorded reason: "a component
    of the event the reads never show apart from it (an SNV sign-out reports beside its MNP) is part
    of the event, not another allele".
  - **Rule X (#297, 2026-10-31):** a listed allele is always a candidate (input taken as correct), and
    a read that cannot tell X from it is partial evidence.
- **Invariants:** 7 (input taken as correct; count only reads that carry the given ALT), and both paths
  judge alike (rules `carrier.rs:2675–2680`; judge `haplotypes.rs:2120–2190`).
- **Principle:** one answer for "a listed allele the reads do not show apart from X", for every allele
  type.
- **Instance:** the audit's layer 3 (cluster C2) found the two rules side by side.

## 2. What already touches it

- RJ-5 (a read fitting two listed alleles equally carries neither for certain), RJ-25, RJ-28 (an
  MNV's carriers count at its component SNV's row, 2026-10-25), rule X in RJ-24 and RJ-25 on the #297
  branch. REJ-20261031-001 (option (b): a listed count decided absent from the reads).
- Input practice: MSK sign-outs list an SNV beside its MNP (the component), so S is often a component
  of X, and its molecules are X's molecules.

## 3. Options

Worked example (allele shapes only): X = CCG>TTA; S = C>T at its first base. 100 reads: 60 REF, 30
whole X, 8 with T at the first base and the other two bases below `--min-baseq`, 2 whole S (TCG).

| Option | The rule | Example (RD, AD, partial) | Invariant 7 | New conditions |
|:--|:--|:--|:--|:--|
| A. Rule X everywhere | Every listed allele is a candidate; RJ-25's sibling floor goes. | (60, 30, 10) | passes | −1 (removes the floor) |
| B. Component rule | A listed allele of X's length whose every change is one of X's changes, with X's base (a component), is part of X at X's row. Every other listed allele is a candidate. No share. | (60, 38, 2) | passes for non-components; for a component, the tie reads may be S's own molecules (2 of 62 here) | 0 (replaces the floor with a structural test) |
| C. Keep RJ-25's floor | As today; rule X's text names the exception and its mechanism. | (60, 38, 2); with S on 5% or more: (60, 30, 10) | a presence test on a listed allele, which rule X rejected | 0 |

The 2 whole-S reads are not ALT under any option (they read X's other bases as REF).

## 4. Measurement (fixed before the runs)

- **Data:** the panel's sign-out MNV rows with a listed component (real MSK input practice), the
  compound-het arm (a1, a2, both), HG002, the gate set, the normals (false ALT), both paths, BWA and
  ABRA2 where available. The spike-in has no listed-component MNV class: say so; add one only if A
  and B differ on the real rows.
- **Metric and bar:** per row and per class, AD, partial and VAF under A, B and C (develop); against
  oracle v3 (reads whose own bases carry X) and, on the compound-het arm, its truth. The rule that
  moves no read whose bases do not carry X into AD, and loses the fewest true carriers, wins. Report
  the rows where S is truly present (whole-S reads on 1 or more fragments).
- **Cannot see:** a listed component whose molecules exist only as tie reads (no whole-S read). No
  data set gives that truth; name it as the residual of option B.
- **Contract tests to keep:** RJ-25's and RJ-28's tests, the S1 rule X tests on the #297 branch.

## 5. Budget and stop rule

Two review rounds. Builds: A and B (C is develop: no build). HPC: one batch per build (HG002 and
patient arrays), plus one rerun. Stop and bring a class question if a third rule appears.

## 6. Work and models

- Sonnet builds A and B behind one measurement switch (never for merge), with red-first tests of the
  worked example on both paths, and the whole gate.
- The main session runs the HPC batches (after approval of this page).
- A fresh Fable or Opus reviewer reads this page, the measurement and the diff.

## 7. State

- 2026-10-31: drafted from the audit (layer 3, C2). Nothing is built.

## 8. Cost log

| Step | Agent and model | Tokens | Notes |
|:--|:--|--:|:--|
| | | | |
