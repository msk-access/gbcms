---
name: no-ticket-labels-in-code
description: Code comments and log messages must explain what/why/how, never cite ticket labels (CR-/HI-/ME-/P4c/milestone IDs).
metadata:
  type: feedback
---

Comments and log strings must be self-contained for a reader with **no plan/PR
context**: state WHAT the code does, WHY, and HOW. Never cite ticket or milestone
labels (`CR-1`, `HI-11`, `ME-8`, `P4c`, `M2`, …) in code or logs — they are
ephemeral and mean nothing to a future contributor. Ticket context belongs in the
commit message, recoverable via `git blame`. (Plan/design docs are the opposite —
labels belong there.)

**Why:** A new person reading `// ME-8: padding fix` learns nothing about the actual
behavior; the label is noise. Ronak flagged this on 2026-06-27 after the mFSD PRs,
where both my new comments and the pre-existing `P4c` annotations carried labels.

**How to apply:** Convert any label into its rationale — `// ME-8: padding fix` →
`// exclude no-test variants; they would inflate the FDR family`. Going forward
(M3+) no labels in code/logs. Existing labels (mine + pre-existing `P4c`) are
tracked for cleanup as plan ticket DX-1. Promoted to a rule in
`.agents/rules/code-quality.md` (Comment & Log Hygiene). See [[harness-stay-small]].
