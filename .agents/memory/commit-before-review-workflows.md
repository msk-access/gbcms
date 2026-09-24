---
name: commit-before-review-workflows
description: Shared checkout hazards — review subagents may git stash; concurrent sessions write .agents/memory — commit before workflows, stage explicit paths only.
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

**Also (2026-09-23):** the memory dir is shared with concurrently running
sessions (it is symlinked into `.agents/memory/`), so `git add -A` sweeps
another session's in-flight memory edits into the current branch's commits
(happened in the T2 PR). Stage explicit paths only; check `git status` for
`.agents/memory/` changes you did not make and leave them unstaged.
