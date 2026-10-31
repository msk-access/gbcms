# Audit: the counting conditions added in 6.7.0, by class or by row

Status: approved 2026-10-31 (scope A; tools local for the trial; first run without a token budget)
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
| 2. Sort | Sonnet, with the rubric below; Haiku for the blind 20% | one row at a time: its reason text only | a tag per row |
| 3. Reason | Fable, once | the rows tagged "rows" or "unclear", and every rule with four or more conditions; graph queries capped at 1,500 tokens | the clusters (section 4) |
| 4. Decide | the operator | the clusters | a decision per cluster; each change gets a plan |

**Rubric for layer 2:**
- **class:** the reason names a read shape or a mechanism without a position ("a read that ends inside
  the repeat", "insertions of 50 bp or more").
- **rows:** the reason names positions, rows or test-set shapes, or "to pass a bar".
- **bug:** the condition makes the code match a rule text that already exists.
- **speed** or **naming:** no count effect by construction.
- **unclear:** no reason is recorded.

## 7. State

- 2026-10-31: the plan is drafted. The baseline of #297, #264 and #369 is measured
  (`~/test/gbcms/harness/plan670/plan-first/baseline.md`; 5.82 M agent tokens for the three, 41 HPC
  batches). Nothing has run.

## 8. Cost log

| Layer | Agent and model | Tokens | Rows | Notes |
|:--|:--|--:|--:|:--|
| 1. Code graph | `graphify update` and `cluster-only --no-label` (no model) | 0 | – | develop a065ab64: 1,582 nodes, 5,114 edges, 45 groups; `~/test/gbcms/harness/plan670/audit/code_develop/` |
| 0. Inventory | Sonnet writes `inventory.py` | | | running |
