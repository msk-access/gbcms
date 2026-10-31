# Plan first: one page before a count-changing build

> **Trial (operator, 2026-10-31).** This page is on trial on the next three count-changing items. It
> becomes part of `.agents/rules/decision-protocol.md` only if the trial meets its bar (section
> "How the trial is judged"). If it does not, it goes to `.agents/learnings/REJECTED.md` with its
> numbers.

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

- **Baseline:** the same measures on #297, #264 and #369 in this cycle, by the same script
  (`~/test/gbcms/harness/plan670/plan-first/baseline.py`, local).
- **Trial items:** the layered audit of rule conditions (the first plan), then the next two
  count-changing items.
- **Bar, fixed now:**
  - at most two review rounds per decision;
  - no condition added and later removed;
  - at most two HPC batches per decision;
  - fewer agent tokens per decision than the baseline's median;
  - the operator can read the page in five minutes.
- **If the trial passes:** the page goes into the decision protocol and replaces text there. A
  count-changing PR must then link its plan, and the qa-check step looks for the link.

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
