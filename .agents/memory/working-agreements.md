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

## Explain so the operator can follow

The operator follows from the app, often after time away; a label helps only once its meaning is in
front of them.
- **Describe first, label in brackets** (2026-10-25: "when explaining dont just say RJ, Phase A etc,
  describe it and mention those in bracket"). Say what the thing is or does, then the label:
  "the fix round on the graph scorer (phase C)", "clipped bases inside the fragment are the read's own
  (RJ-13)", "the build after the fix round (`src_c50d`)". Tables name rows in words too. A bare label
  is fine only in code, commits and GitHub records, which keep their own conventions.
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

## Also standing

- Comments and log strings never cite ticket labels (`.agents/rules/code-quality.md`).
- A failing test's expectation is never edited without telling the operator first (`AGENTS.md`).
