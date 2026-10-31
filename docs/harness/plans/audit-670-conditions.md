# Audit: the counting conditions added in 6.7.0, by class or by row

Status: decided 2026-10-31 (approved the same day: scope A; tools local for the trial; first run without a token budget)
Rounds: 0 of 2 · New conditions: 0 of 0 (an audit adds none) · HPC batches: 0 of 0

## 1. The class

- **What is audited:** every condition that changes how a read is counted, added or amended in 6.7.0.
  A condition is a threshold, an exception, a special case for one path, a judge state or window, a
  diagnostic name, or an environment switch (`GBCMS_*`). Both paths are in scope: the default rules
  path and the sequence judge (`--read-judge sequence`). Merged and in-flight work are both in scope.
- **The question for each condition:** is its reason a class of reads (a mechanism, stated without a
  row), or named rows? Does it keep the invariants, above all invariant 7 ("input taken as correct")?
- **Principle:** a rule follows from a mechanism, and a row is one instance of it (register: "fix the
  class, not the row").
- **Instance:** #297 gained about ten conditions in four review rounds; rule X (2026-10-31) removed
  eight of them.

## 2. What already touches it

- The register `docs/reference/read-judgment.md` (RJ-1 to RJ-30, 815 lines) and the judge design
  `docs/proposals/sequence-judge.md` (886 lines).
- The counting code: about 38,000 lines in `rust/src` (`counting/`, `judge/`, `normalize/`), 14
  environment switches, and 86 commits in `rust/src/counting` on develop since 2026-10-23.
- In flight: `fix/264-equal-windows` (batches z64 and zn64 on HPC) and `fix/297-narrow-unit-run` (the
  rule X rebuild, under way).
- `.agents/learnings/REJECTED.md` and the decision protocol, section 3.

## 3. Options (the scope)

| Option | Scope | Rows to sort (rough) | Recommended |
|:--|:--|:--|:--|
| A | 6.7.0 only: since develop `ca64099a` (2026-10-23), both paths, merged and in flight | about 100 to 300 | yes, first |
| B | Every condition in the counting code, whatever its age | about 3 times A | only if A finds drift beyond #297 |
| C | The judge only, since its first build | about half of A | no: the rules path is the 6.7.0 default |

None of these options changes a count, so the invariant check applies to the audit's findings, not to
the options.

## 4. Measurement (fixed before the runs)

- **Completeness:** every register amendment dated 2026-10-23 or later, and every added condition in
  the diff `ca64099a..origin/develop` and in the in-flight branches, is a row of the inventory. A count
  of the diff's changed condition lines against the inventory's rows checks it.
- **Sorting accuracy:** a fresh Haiku agent tags 20% of the rows blind. Agreement with the Sonnet tags
  is 90% or more, or the rubric is fixed and the sort runs again. The operator spot-checks 10 rows.
- **The output:** clusters of conditions. Each cluster has its class question, the invariant check, a
  class-level form, the conditions that form would remove, and the measurement that would settle it.
- **What this cannot see:**
  - a reason recorded only in conversation (not in a commit, the register or a test);
  - older code that interacts with a new condition;
  - count effects. The audit runs nothing; each change it proposes gets its own plan and measurement.

## 5. Budget and stop rule

- **First run: no token budget** (operator, 2026-10-31: "as we are running it for the first time run
  without budget"). The cost log records the tokens each layer uses; those numbers set the budget
  for later runs. The estimate was Sonnet about 300k and Fable about 400k.
- Layers 0 and 1 use no model tokens to run; a Sonnet agent writes the inventory script.
- Layer 3 (Fable) runs once, on the narrowed rows only.
- If the inventory has more than 400 rows, report the counts by module and by kind to the operator,
  then go on.
- The audit changes no code. A proposed change waits for its own plan.

## 6. Work and models

| Layer | Who | Input | Output (local, `~/test/gbcms/harness/plan670/audit/`) |
|:--|:--|:--|:--|
| 0. Inventory | a script that Sonnet writes (stdlib and git only) | the register, the diffs, the tests | `conditions.tsv`: id, kind, `file:line`, commit, date, ticket, rule, reason text, "cites rows" mark, path (rules, judge, both, shared) |
| 1. Code graph | `graphify rust/src --no-viz` (syntax parser only; no docs, so no model) | the Rust code | `graphify-out/` (kept out of git) |
| 2. Sort | Sonnet, with the rubric below; Haiku for the blind 20% | **one decision unit (commit) at a time** (operator, 2026-10-31): its subject, body and its rows' reason text; `group_by_commit.py` (no model) makes the units and drops the QC-summary rows | a tag per unit |
| 3. Reason | Fable, once | the rows tagged "rows" or "unclear", and every rule with four or more conditions; graph queries capped at 1,500 tokens | the clusters (section 4) |
| 4. Decide | the operator | the clusters | a decision per cluster; each change gets a plan |

**Rubric for layer 2** (run 2; run 1 agreed 12 of 17 = 71%, below the bar, so the rubric was fixed):
- **class:** the record states a mechanism that holds without any named row ("a read that ends inside
  the repeat", "insertions of 50 bp or more"). A cited row as evidence does not change this.
- **rows:** take the named rows out of the record; if no reason is left, it is "rows". This includes
  "to pass a bar" and "to move row X".
- **A revert or a narrowing:** tag the reason for it. A mechanism is "class"; "it failed bar X" or "it
  moved row Y" is "rows".
- **bug:** the change makes the code match a counting rule that already exists (a register row RJ-n,
  a decided rule or an invariant), and the record names that rule.
- **hygiene:** a lint, code-quality or refactor change with no count effect (new in run 2).
- **speed**, **naming** or **measure:** no count effect on the default behaviour by construction.
- **unclear:** the record says what changed but not why.
- **Input, the same for both sorters:** the unit's subject and its whole commit body (up to 3,000
  characters), and its rows' reason text. The blind sample of run 2 is drawn from the 66 units run 1
  did not use.

## 7. State

- 2026-10-31: all layers ran. Layer 3 (`~/test/gbcms/harness/plan670/audit/layer3_clusters.md`) put
  the 367 conditions of 26 units in 10 clusters by the reads they act on. Six are sound.
- **Decided (operator, 2026-10-31):**
  - **C2 (high):** reads that cannot tell X from a listed allele have two decided rules: RJ-25's
    sibling floor (1 fragment and 5%) and rule X (a listed allele is always a candidate). Plan and
    measure first: `plans/c2-listed-siblings.md`.
  - **C1 (medium):** three presence tests for one question, and a second copy of RJ-27's tolerance.
    Taken into #297's close plan.
  - **C3 (medium, judge):** the REF tie-break for every allele type. Trace and plan before the judge
    becomes the default (#331).
- **Recorded for later (low):** C5, write U049's mechanism into RJ-28, and name the silent
  `.take(32)` cap; C6, `FOUND_MIN_READS = 2` is still provisional; C8, check `GBCMS_JUDGE_BAND=full`
  is byte-identical on every arm; C9, name the two contig-edge tests and their mechanism in RJ-3.
- **Lessons for the next run:** fix the sorting bar on the split that decides layer 3, before the run;
  graphify's `explain` and `path` helped, its broad `query` output did not.

## 8. Cost log

| Layer | Agent and model | Tokens | Rows | Notes |
|:--|:--|--:|--:|:--|
| **Total** | | **1,258,000** | | below the baseline median (2.38 M); no review rounds, no HPC |
| 3. Reason | Fable, one pass | 229,341 | 26 units, 367 conditions | 10 clusters; 1 high (C2), 2 medium (C1, C3), the rest low or sound |
| 1. Code graph | `graphify update` and `cluster-only --no-label` (no model) | 0 | – | develop a065ab64: 1,582 nodes, 5,114 edges, 45 groups; `~/test/gbcms/harness/plan670/audit/code_develop/` |
| 2. Sort, run 2 | Sonnet (83 units); Haiku (a fresh blind 17) | 167,330; 119,051 | 83 | **bar not met on the full tag:** 15 of 17 (88%); both differences are between two "no named rows" tags (hygiene/class, class/bug). **Met on the split that decides layer 3** (rows or unclear, against the rest): 17 of 17. Operator (2026-10-31): go on, conservatively; layer 3 takes run 2's 24 units plus every unit either run tagged rows or unclear: 26 units, 367 conditions. 19 of 83 tags changed between runs, one across the class/rows line (U083). Lesson for later runs: fix the bar on the decision-relevant split before the run |
| 2. Sort, run 1 | Sonnet (83 units); Haiku (blind 17) | 237,119; 132,221 | 83 | agreement 12 of 17 (71%), below the 90% bar: the inputs differed (Sonnet read the full messages) and the class/rows line was unclear; rubric fixed, run 2 |
| 0b. Units | `group_by_commit.py` (no model) | 0 | 83 | 841 rows kept (110 QC-summary rows dropped) in 83 commits: merged 71, in flight 12; 11 name positions in their message; blind sample of 17 drawn with seed 20261031 |
| 0. Inventory | Sonnet writes `inventory.py` | 373,868 | 951 | merged 855, in flight 96; judge path 687 (the judge's first build landed in 6.7.0); 110 rows in QC-summary code that touches no read; the "cites rows" mark is noisy (376 of 419 match only the word "row") |
