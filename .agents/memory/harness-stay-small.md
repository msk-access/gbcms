---
name: harness-stay-small
description: "The user wants the Claude harness (skills, rules, standing context) kept deliberately small and curated."
metadata: 
  node_type: memory
  type: feedback
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
---

The user adopted the "Claudelicious / The Harness Is the Moat" philosophy for this
repo: the harness is the durable advantage, and its third property is that **it stays
small**. Skills were curated, not accumulated — redundant ones (git-commit,
code-review, release) were dropped and folded into `CLAUDE.md` / `git-conventions.md`.

**Why:** "Every rule the model reads, every skill in its reach, is something it holds
in attention on every turn, and attention is finite." Compounding means keeping what
works and cutting what doesn't, not piling on scaffolding.

**How to apply:** When adding a rule/skill/memory, first check whether it duplicates
existing context. Prefer editing an existing artifact over creating a new one. Keep
`CLAUDE.md` tight; push depth into on-demand reference docs ([[claudelicious-reference]]).
