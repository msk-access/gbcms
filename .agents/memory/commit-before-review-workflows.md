---
name: commit-before-review-workflows
description: Adversarial-review subagents share the checkout and may git stash/checkout — commit work before launching them.
metadata:
  node_type: memory
  type: feedback
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
  modified: 2026-09-23T04:04:56.438Z
---

During the T1 cycle (2026-09-23) a review-workflow subagent ran `git stash` on
the shared checkout and silently reverted my uncommitted engine edits; I only
caught it because a later scripted edit's `assert` failed.

**Why:** workflow agents run in the same working tree unless
`isolation: 'worktree'` is set; any of them may mutate git state while
"verifying HEAD".

**How to apply:** commit (or at least stash-tag yourself) before launching
same-checkout review workflows, or pass `isolation: 'worktree'` for agents
that run git commands. After any workflow completes, check `git status` and
`git stash list` before further edits. Related: [[user-local-validation-data]].
