# CONTINUITY — where we left off

> Tactical state that must survive a closed laptop or a context summary.
> Update the **Now** and **Next** sections as work progresses.

_Last updated: 2026-06-26_

## Now
Cross-tool **harness** for gbcms is built and committed on **`feat/claude-harness`**.
**PR #21 → `develop` is open** (https://github.com/msk-access/gbcms/pull/21),
awaiting review. Following the Claudelicious pattern ("The Harness Is the Moat").

## Goal
Stand up the 6-piece harness spine (rules, skills, hooks, memory, continuity,
learning loop) in the Claude Code-native layout. Plan: `docs/harness/00-plan.md`.

## Done this session
- Branch `feat/claude-harness` cut.
- Built the 6-piece spine, then **refactored to cross-tool layout**:
  - **Standing context:** `AGENTS.md` canonical (cross-tool); `CLAUDE.md` = `@AGENTS.md` import.
  - **Skills (11):** canonical in `.agents/skills/` (Antigravity-native); `.claude/skills/`
    is a symlink. Added `fragment-counting`, `read-filters-qc`, `nextflow-pipeline`;
    dropped git-commit/code-review/release; workflows folded into skills.
  - **Hooks:** `.claude/settings.json` (time injection + `guard-destructive.sh`, tested);
    plus tool-agnostic `.githooks/pre-commit` (data-file + ruff/black guard).
  - **Memory:** canonical in-repo `.agents/memory/` (MEMORY.md + 6 facts); Claude recall
    dir symlinked to it (one source, still auto-recalled).
  - **Continuity:** this file. **Learning loop:** convention in `AGENTS.md`.

## Next
1. Get PR #21 reviewed and merged into `develop`.
2. Enable the git hook on other clones: `git config core.hooksPath .githooks`.
3. Verify the harness in a fresh session (success criteria in `docs/harness/00-plan.md`),
   ideally in both Claude Code and Antigravity (symlink traversal + memory recall).
4. Then resume the engineering backlog: `CODE_REVIEW_IMPLEMENTATION_PLAN.md`
   (45 tickets, milestones M1–M6; start M1 — count-correctness: CR-1, CR-3, CR-2, CR-4).

## Key decisions made
- Full 6-piece spine · **cross-tool layout** (AGENTS.md canonical, per-tool adapters) ·
  follow claudelicious upstream as a pinned reference.
- Skill library curated to 11 (no overlap) per the "stay small" discipline.
- `.agents/` is the neutral source of truth (rules, skills, memory); `.claude/` mirrors it.
- Deferred (claudelicious 07–19): session search, second brain, multi-provider,
  crons, always-on agents, fleet.
