---
name: claudelicious-reference
description: The gbcms harness follows the Claudelicious pattern; where to find it and what it covers.
metadata: 
  node_type: memory
  type: reference
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
---

This project's Claude harness is modeled on **Claudelicious** by Justin Johnson
(Run Data Run), per the article "The Harness Is the Moat" (2026-06-25).

- Upstream: https://github.com/BioInfo/claudelicious (branch `main`, referenced 2026-06-26).
- In-repo plan & mapping: `docs/harness/00-plan.md`.
- Six harness pieces: rules/context, skills, hooks, memory, continuity, learning loop.
- We adopt the **shape and discipline**, not a clone. Layout is Claude Code-native
  (`CLAUDE.md`, `.claude/skills/`, `.claude/settings.json`). Following = pinned reference,
  not a submodule.
- License: claudelicious code MIT, prose CC BY-NC 4.0 — attribute when adapting prose.

Deferred (claudelicious docs 07–19): session search, second brain (Qdrant),
multi-provider routing, crons, always-on agents, the fleet. See [[harness-stay-small]].
