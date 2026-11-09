# Plan first: one page before a count-changing build

> **Trial (operator, 2026-10-31).** This page is on trial on the next three count-changing items. It
> becomes part of `.agents/rules/decision-protocol.md` only if the trial meets its bar (section
> "How the trial is judged"). If it does not, it goes to `.agents/learnings/REJECTED.md` with its
> numbers.
>
> **Result (2026-11-05): the bar was not met** (section "Result"). Operator, 2026-11-06: try it
> again on the 6.8.0 work, with a new bar fixed before the first count-changing build. Until that
> trial passes, the page is not part of the decision protocol.

The decision protocol and the register ("How a change to read judgment is made": spec first, decide,
then code) already say what a count-changing decision needs. Until now those steps were done when an
option went to the operator, after the code existed. Then each review round added a condition, and
each condition moved the boundary to the next row (#297: four rounds, about ten conditions, eight of
them removed by the final rule). A plan does the same steps **before the first build**, when they
cost the least.

## When a plan is needed

- **Needed:** any change that can move a count. This includes a new rule, a threshold, a path
  special case, a new state or window in the judge, and an option the operator must choose.
- **Not needed:** hygiene, platform, docs, and speed work that is count-neutral by construction.
  Such work still proves that it is neutral (the binning-invariance and identity checks).
- **Where:** `docs/harness/plans/<ticket>-<short-name>.md`. Link it from the ticket and the PR. HG002
  positions are public and may appear. Patient rows never appear.

## The page

Keep it to one page. Write it in plain English (ASD-STE100), with each label in brackets after its
plain meaning. Each section is short; a section that does not apply says "none" and why.

1. **The class.** State the mechanism without the row that showed it: which reads, which evidence a
   read must carry, and which invariants (`AGENTS.md`) the change touches. Give the principle in one
   sentence. Put the row that showed it at the end, as "instance".
2. **What already touches it.** List the rules (RJ-n, read to the end of the row), the code
   conditions on this class (`file:line`, from the condition inventory and the code graph), the
   entries in `REJECTED.md`, and the open and closed tickets with measurements.
3. **Options.** Give two or three class-level options. For each:
   - the rule in one sentence;
   - the invariant check, invariant 7 ("input taken as correct") first: pass or fail, and why;
   - binning invariance, the one quality contract, and whether both paths judge it the same way;
   - the conditions it adds (a number), and what it will not handle (the named limit).
4. **Measurement, fixed before the runs.** Give:
   - the class-level metric and its bar (rows are checks, not reasons);
   - the data types and truth sources (protocol section 7), and which of them are missing;
   - what the measurement cannot see (protocol section 6);
   - the contract tests the change must keep, and the HPC batches planned (how many, how large).
5. **Budget and stop rule.** At most two review rounds. At most the number of new conditions that
   the chosen option declares. Past either limit, stop: the work goes back to the operator as a class
   question with one class-level option. HPC batches: the planned number, plus one rerun.
6. **Work and models.** Use the cheapest layer first:
   - scripts and the code graph (`graphify`, code only), which cost no model tokens;
   - Sonnet to sort, build and test; Haiku for scans;
   - Fable or Opus for the review and for the class reasoning, on the narrowed input only.
   Each brief points to this page and does not retype it. An agent stops at the end of its phase.
   The next agent starts fresh from section 7.
7. **State.** At most ten lines, updated at the end of each phase: what is done, what is measured,
   and what is open. This is the hand-off to a fresh agent.
8. **Cost log.** Review rounds, conditions added and removed, HPC batches and tasks, agent tokens
   (from each agent's usage report), and days from the plan to the decision.

The operator approves the page before the first build. The approval is the date and "approved" in
the page header. Each reviewer gets the page and the diff, and answers two questions: is the defect
in this instance or in the mechanism, and does an option break an invariant?

## How the trial is judged

- **Baseline** (measured 2026-10-31 by `~/test/gbcms/harness/plan670/plan-first/baseline.py`, local;
  agent tokens are summed over context resets; the main session's own tokens are not visible):

  | Measure | #297 (unit-run rule) | #264 (long-insertion windows) | #369 (long-insertion band) |
  |:--|--:|--:|--:|
  | Agent tokens | 3.02 M | 2.38 M | 1.03 M |
  | Review rounds | 4 | 3 | 2 |
  | HPC builds, batches, array tasks | 6, 17, 496 | 5, 10, 174 | 2, 14, 320 |
  | Commits undone | 5 (of 11) | 0 | 0 |
  | Decision reversed later | yes (option (b)) | no | no |

- **Trial items:** the layered audit of rule conditions (the first plan), then the next two
  count-changing items.
- **Bar, fixed now** (per decision):
  - at most two review rounds;
  - no condition added and later removed, and no decision reversed;
  - HPC builds at most the options the plan measures, plus one final build; batches fewer than the
    baseline's median (14);
  - agent tokens below the baseline's median (2.38 M);
  - the operator can read the page in five minutes.
- **If the trial passes:** the page goes into the decision protocol and replaces text there. A
  count-changing PR must then link its plan, and the qa-check step looks for the link.

## Result (measured 2026-11-05; operator 2026-11-06)

Measured with `baseline.py`'s method (`~/test/gbcms/harness/plan670/plan-first/trial_results.md`
and `trial_config.json`): agent tokens summed over context resets; the main session's own tokens
are not visible.

| Bar | 1. Audit | 2. #297 close | 3. C2, then class plan B |
|:--|:--|:--|:--|
| Agent tokens below 2.38 M | 1.26 M: pass | 3.20 M: fail (2.64 M without the attributions marked unsure) | 6.01 M: fail (C2 1.56 M, B 4.45 M) |
| At most two review rounds | 0: pass | 2: pass | 4: fail (two for each plan page) |
| No decision reversed; no condition added and later removed | pass | pass | fail: C2's choice changed twice before B (2026-10-31, 2026-11-02); the refinement is a stopgap that #418 removes |
| HPC builds at most the options measured plus one final | 0: pass | 3 for 2: fail | 7 (C2 2, B 5): fail |
| HPC batches fewer than 14 | 0: pass | 8: pass | 13: pass |
| The operator reads the page in five minutes | not measured | not measured | not measured |

- **Better than the baseline:** no commit was undone (the baseline's #297 undid 5 of 11), and no
  merged decision was reversed.
- **Where the cost went:** the scorers used 2.75 M of item 3's 6.01 M tokens. The measurements
  found two real bugs (rq1); the plan did not prevent them.
- **Extra items (not in the bar):** #292, 3.08 M tokens (2.36 M without the prototypes made before
  its plan), 2 rounds, 4 builds, 8 batches; RNA package 3's lean path, 0.39 M tokens, no review
  round, no HPC, merged the day it was planned.
- **The code graph's own trial was not run:** no later plan answered "what already touches it" both
  ways. graphify stays outside the repository.
- **Verdict:** the bar was not met. **Operator (2026-11-06):** try it again on the 6.8.0 work. The
  new bar is fixed with the operator before the first 6.8.0 count-changing build.

## The code graph: a trial of its own (operator, 2026-10-31)

graphify (0.9.12) builds a graph of the Rust code with its syntax parser, at no model cost. In the
audit, `explain` and `path` helped; the broad `query` did not. Before it is wired into the project
(an MCP server in `.mcp.json`, a rebuild after each merge to develop), measure it:

- **Test:** in the next two plans, answer each "what already touches it" question twice: with the
  graph (`explain`, `path`, `affected`) and with grep or a scan agent. Record the tokens of each way,
  and any call site either way misses.
- **Bar, fixed now:** the graph saves 30% or more of the tokens and misses no call site.
- **If it passes:** the MCP server at project scope, pinned to 0.9.12; the free code-only rebuild
  after each merge to develop; `graphify-out/` kept out of git; one pointer line in the decision
  protocol, not a section in `CLAUDE.md`; a usage note (prefer `explain`, `path` and `affected`;
  `query` only with a context filter).
- **Until then:** the graph stays outside the repository (a `graphify-out/` in the repo would make
  the user-level skill treat every codebase question as a graph query first).

## Skeleton

```markdown
# <ticket>: <short name>
Status: draft | approved <date> | building | measured | decided <date>
Rounds: 0 of 2 · New conditions: 0 of <n> · HPC batches: 0 of <n>

## 1. The class
## 2. What already touches it
## 3. Options
## 4. Measurement (fixed before the runs)
## 5. Budget and stop rule
## 6. Work and models
## 7. State
## 8. Cost log
```
