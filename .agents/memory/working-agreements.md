---
name: working-agreements
description: "How the operator wants me to work: keep the harness small, explain in Simplified Technical English with internal labels in brackets, pick each subagent's model by task."
metadata:
  type: feedback
---

## Keep the harness small

The harness is the durable advantage and its third property is that it stays small ("The Harness
Is the Moat", Justin Johnson; the claudelicious pattern, pinned in `docs/harness/00-plan.md`). Every
rule the model reads, every skill in its reach, is held in attention on every turn, and attention is
finite. Skills are curated, not accumulated: git-commit, code-review and release were folded into
`AGENTS.md` and `git-conventions.md`, and debug-error into the rust-python-ffi skill.

**How to apply:** before you add a rule, skill or memory, check that it does not duplicate existing
context; edit an existing artifact before you create one; keep one home per fact; push depth into
on-demand docs and rules. Targets (#317): `MEMORY.md` at most 5 KB, `AGENTS.md` at most 75 lines,
skill descriptions short, `CONTINUITY.md` holding only the current state ("Now" and "Next"; finished
sections go to `.agents/archive/`).

## Testing is never capped

Builds for testing, review rounds and HPC test batches are not limited by a number (operator,
2026-11-17: "why do we have limits, we should sumbit to HPC as required ... testing is core part so we
should not limit it"). A plan counts them in its cost log; it stops by content, not by count: an output
difference it did not accept, more new conditions than declared, or a change of design or decision.

**How to apply:** never write a "N at most" cap on builds, reviews or HPC batches in a plan, and never
stop a test because a count was reached. The PHI and HPC rules still hold: copy back counts, times and
profiles only, and log every HPC step. Home of the rule: `docs/harness/plan-first.md` §5.

## Explain so the operator can follow

The operator follows from the app, often after time away; a label helps only once its meaning is in
front of them.
- **Describe first, label in brackets** (2026-10-25: "when explaining dont just say RJ, Phase A etc,
  describe it and mention those in bracket"). Say what the thing is or does, then the label:
  "the fix round on the graph scorer (phase C)", "clipped bases inside the fragment are the read's own
  (RJ-13)", "the build after the fix round (`src_c50d`)". Tables name rows in words too. A bare label
  is fine only in code, commits and GitHub records, which keep their own conventions.
- **Build, job and run labels are labels too** (2026-10-30: "again just giving me bunch of
  abbrevations does not help, need to describe each and put abbreviations in brackets"). Status
  tables and relayed agent reports carry the most labels (`c61hs`, `hyR`, item A). Rewrite each one
  in words before it reaches the operator, and ask agents to report in this form. When work comes
  back after the operator thought it was done, say why it is back.
- **Write in ASD-STE100 Simplified Technical English** (2026-10-28):
  - Use 20 words at most per instruction sentence and 25 per description. Give one instruction per
    sentence, in the imperative.
  - Use active voice and simple tenses (present, past, future). Use articles ("the", "a").
  - Keep one word for one meaning: a "row" is always a row, a "read" always a read.
  - Do not use "-ing" forms as nouns or adjectives, except in technical names. Keep noun clusters to
    three words at most.
  - Write paragraphs of six sentences at most. Use vertical lists for steps and for complex text.
  - Technical names, tickets and labels stay as they are (code style), with the plain meaning first.
    Decisions still come with worked examples (`.agents/rules/decision-protocol.md`).

## Model per task

Set each delegated agent's `model` explicitly instead of inheriting the session's, and say which
(operator, 2026-10-12: "choose a model based on the task ... if it requires more reasoning and
thinking we can use fable and opus depending on which is suitable"):
- **Deep reasoning** (adversarial reviews, rule and spec design, decision analysis, hard debugging,
  holistic effect maps): Fable or Opus, whichever suits the task. Reviews of counting rules get the
  strongest reasoning model.
- **Research surveys and synthesis** (literature, other tools' methods, design options): Fable
  (operator, 2026-10-23: "We can use fable for research").
- **Routine implementation and test writing:** Sonnet.
- **Mechanical work** (file scans, listings, formatting, simple greps): Haiku.

When the session's work changes character (for example from measurement plumbing to rule design),
suggest switching the session model. `CLAUDE.md` carries the short form.

The rule holds for follow-up work too (2026-10-30: "Are we still delegating properly to different
models depending on the task complexity"). Do not hand routine aggregation, tables or scans to a
long-lived Opus agent because it is already there: start a Sonnet (or Haiku) agent with a full brief.
Before a count-changing option goes to the operator, a fresh Fable or Opus agent reviews it
adversarially, from the code and the numbers, not from the author's conclusion.

## Also standing

- Comments and log strings never cite ticket labels (`.agents/rules/code-quality.md`).
- A failing test's expectation is never edited without telling the operator first (`AGENTS.md`).
- **Decide by class, not by row** (2026-10-31: "i feel we regressed and made certain row level
  decisions vs category level"). Review findings do not go back to a builder as a fix list after the
  first round until the rule's conditions are counted (`.agents/rules/decision-protocol.md` §3).
- **Check each signed-out row that loses its ALT, then decide by class** (2026-11-09: "I dont
  understand how can Signout event got to AD 0 ... we should check and understand each"). A class
  count does not explain a signed-out call that falls to AD 0: look at its listed neighbours, the
  changes its ALT reads share and each read judgment's count (`docs/development/regression-panel.md`).
- **Measure what is built before deciding** (2026-11-19: "we should measure if we already built"). A decision
  never goes to the operator while a built option lacks a run on the row sets the recommended option ran on; the
  coverage check is in `.agents/rules/decision-protocol.md` §8.
- **Fix at the core, for accuracy; no patches unless asked** (2026-11-20: "we should solve problem at it core/source
  keeping accuracy in mind, and not patches unless i explicitly ask for patches"). Find the cause and change it where
  it lives (the judge, the read model, the shared contract), so the whole class is right and the counts stay accurate.
  A patch is a narrow change that hides a symptom: a special case, a row-level exception, a skip, a threshold moved to
  pass, a switch, a workaround around the cause. Do not propose or build one unless the operator asks for a patch by
  name. When only a patch is possible now, say so, say what the core fix is, and let the operator decide. Prototype
  switches in a measurement build are measurement tools, not the fix: the decided rule goes into the core code.
  **The test is accuracy against truth, not unchanged counts** (2026-11-20: "it should be about regression, we
  want to be closer to the truth with less FP and FN increasing accuracy"). A rule ships when it lowers false ALT
  (FP) and missed ALT (FN) on rows of known truth (HG002's benchmark, FORTE rows with matched DNA, catalogue-absent
  rows). A regression is any class that moves away from the truth; such a rule is fixed at the core until no class
  is worse. Shipped code has one behaviour, the most accurate; a kept option is a documented CLI option with tests,
  never a hidden environment switch.
- **Read an unknown word as a typo first** (2026-11-19: "*enough What about the 50 forte bams dont they have
  enough data?"). A word that names nothing in the notes is first read as a typo of a common word. If I still ask,
  the plain reading is the first option, and the data already in hand comes before new data.
- **Every failed task is reported** (2026-10-30: "when i see failed task should i probe you"). When a
  background task, HPC job or agent run fails, the next message says what failed and what follows:
  fixed and rerun, stopped on purpose, a red-first test failing as intended, or still open. The
  operator should never need to ask.

## Local first while CI minutes are short (operator, 2026-11-17)
"we are running currently out of github actioins minutes assigned to us and we should do most of it local untill we
are sure". Build, test and gate on this machine (the qa-check gate in full, the local gate arms); push a branch or
open a PR only when the work is sure, so CI runs once per finished change. HPC is not GitHub's and stays available
under its own scope.
