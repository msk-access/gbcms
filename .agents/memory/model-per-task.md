---
name: model-per-task
description: Choose each delegated agent's model by the task (Fable for research surveys, Fable or Opus for deep reasoning, Sonnet for routine work, Haiku for mechanical scans) and say which; suggest switching the session model when a phase changes.
metadata:
  type: feedback
---

Pick the model for each subagent explicitly rather than inheriting the session's:
- **Deep reasoning** (adversarial reviews, rule/spec design, decision analysis, hard
  debugging, holistic effect maps): Fable or Opus, whichever suits the task better.
- **Research surveys and synthesis** (literature, other tools' methods, design
  options): Fable (operator, 2026-10-23: "We can use fable for research").
- **Routine implementation and test writing:** Sonnet.
- **Mechanical work** (file scans, listings, formatting, simple greps): Haiku.
Name the model choice in the delegation. When the main session's work changes
character (e.g. from measurement plumbing to rule design), suggest switching the
session model.

**Why:** operator, 2026-10-12: "choose a model based on the task for example if it
requires more reasoning and thinking we can use fable and opus depending on which is
suitable and change the model based on requirements".

**How to apply:** every Agent call sets `model`; reviews of counting rules get the
strongest reasoning model. Related: [[harness-stay-small]].
