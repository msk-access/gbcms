# #297: close the unit-run rule on one principle

Status: draft, for the operator's approval
Rounds: 1 of 2 (the final review of rule X) · New conditions: 0 of 0 (the options remove or align conditions) · HPC batches: 0 of 2

## 1. The class

- **Reads:** reads at a long deletion (50 bp or more) inside a tandem repeat, where another copy
  number of the same repeat (another "unit count") is a candidate.
- **Evidence:** a read is ALT only where its own confident bases carry the given ALT **and** tell it
  apart from every other candidate. The candidates are REF, every listed allele (always; rule X), and
  each found allele that passes the presence test (RJ-25). A read that fits the ALT and another
  candidate equally is partial evidence. A read whose bases show another candidate is partial
  evidence (RJ-6: another allele, never ALT).
- **Invariants:** 7 (input taken as correct; judge bases, not placement), 2 (one quality contract),
  1 (binning invariance), and both paths judge alike.
- **Principle:** one test for every candidate source and both paths: "do the read's own confident
  bases tell the given ALT from every candidate?"
- **Instance:** the final review of rule X (2026-10-31) found four places where the branch departs
  from this, three of them on the rules path only.

## 2. What already touches it

- The branch adds 13 conditions (the reviewer's inventory). None rests only on named rows. Two
  started from single rows: the listed-neighbour exception (HG002 4:3610329) and the nearest-count
  flag (two spike-in rows).
- Rules RJ-6, RJ-16, RJ-24, RJ-25, RJ-26, RJ-28 (a row's count must not depend on an unrelated input
  row), RJ-29. Rejected: option (b), REJ-20261031-001.
- The judge already follows the principle in F1, F2 and F3 below. The rules path is the 6.7.0
  default, and the judge becomes the default in 6.8.0 (#331).

## 3. Options

The review's findings, each with the principle's answer (option A) and the alternative (option B):

| Finding | What happens now | A: the principle on both paths | B: keep, and name the limit |
|:--|:--|:--|:--|
| F1 (high): the listed-neighbour exception (`junction.rs:403-411`) | One unrelated listed row (a 1-bp deletion 270 bp away, or an SNV no read carries) turns 6 tie reads into ALT on the rules path: (RD, AD, partial) (6, 0, 6) → (6, 6, 0). The judge stays (6, 0, 6). | A read counts as "the deletion with its listed neighbour" only where its confident bases show the neighbour's change. Keeps the HG002 4:3610329 test. | Not offered: it breaks invariant 7 and RJ-28. |
| F2: a listed count past the longest read (`:987`) | Rules (6, 12, 0): the listed count's 6 molecules count ALT. Judge (6, 0, 12). | Hold on both paths (partial), named `UNIT_COUNT_UNSPANNED`. Cost: such a row has AD 0 until a read can tell. | Keep counts on both paths. This credits reads that cannot tell, against rule X's reason. |
| F3: a read that spells an unlisted count below the presence floor | Rules count it ALT (AD 16 against the judge's 12 in the test). | Partial evidence on both paths (RJ-6). The floor still decides only the reads that cannot tell. | Keep: RJ-24 names it as a limit. |
| F4: deletion carriers with an unlisted 1-bp insertion read like one more unit deleted | Rules hold all 6 reads; judge (6, 3). | Open: measure first. Trace the test's reads on both paths; decide in round 2. | Name the limit on the rules path. |
| F6: the presence test reads raw qualities over the whole read (`:1024`); the hold reads BAQ-adjusted ones in the readable range (`:633`) | One condition, two quality rules (invariant 2). Matters for RNA clip edges and `--apply-baq`. | One rule: the hold's. | Not offered: it breaks invariant 2. |
| F7: an unlisted count that first differs 4 to 23 bases before the longest read's end | In reach (`:987`) but no read can spell it in full (`:1036`): never held, never named. | Reach uses the same 24 bases as the spelling test. | Name the gap. |

Routine fixes in all cases: every counting test asserts the four invariants (F5: four tests miss
some); the judge check at `:212` becomes exact; the spike-in spec's "or one unit inserted" line goes
(F8).

## 4. Measurement (fixed before the runs)

- **Contract tests:** the reviewer's probe (`scratchpad/rev297x/probe2.py`) becomes a red-first test
  of F1; one red-first test each for F2, F3, F6 and F7, on both paths; the S1 rule X tests and HG002
  4:3610329 stay green.
- **Bar:** on every arm, the rules path and the judge agree on every read the principle decides; no
  read whose bases do not carry the given ALT is ALT (the patient oracles v3 and v2 at the run's
  `--min-baseq`); HG002 truth classes do not lose true ALT reads beyond the F2 holds, which are named.
- **Runs:** HG002 c246, tr, sub and nocl locally, both paths, against develop (hfS, hfR) and rule X
  (n100). One HPC batch on every arm; one spike-in phase 1b rerun (S1 to S6, both aligners).
- **Cannot see:** real rows where a listed count differs only past every read (F2's cost); the spike-in
  has few. The HPC batch counts them by `UNIT_COUNT_UNSPANNED`.

## 5. Budget and stop rule

Two review rounds in all (this was round 1). No new condition: each fix tightens an existing
condition to the principle or removes one. If round 2 finds another departure, stop and bring it as
a class question. HPC: one build and one batch, plus one spike-in rerun.

## 6. Work and models

- Sonnet builds the red-first tests and the fixes, one commit per finding, with the whole gate.
- A fresh Opus or Fable reviewer (round 2) gets this page and the diff.
- The main session runs the HPC batch after approval.

## 7. State

- 2026-10-31: rule X is rebuilt (969ea25d) and merged with develop; the gate passes; the final
  review found F1 to F8. Nothing is built for them yet.

## 8. Cost log

| Step | Agent and model | Tokens | Notes |
|:--|:--|--:|:--|
| Rule X rebuild | Sonnet | 261,600 | |
| Final review (round 1) | Opus | 270,039 | found F1 to F8 |
