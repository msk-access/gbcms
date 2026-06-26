# CONTINUITY — where we left off

> Tactical state that must survive a closed laptop or a context summary.
> Update the **Now** and **Next** sections as work progresses.

_Last updated: 2026-06-26_

## Now
Cross-tool Claude/Antigravity **harness** for gbcms, following the Claudelicious
pattern. The spine plus the security and learning layers are built. Two harness
PRs are **open against `develop`**; the earlier two are **merged**.

## Goal
A vendor-neutral harness (AGENTS.md canonical, in-repo memory) that stays small,
learns from corrections, and enforces guardrails deterministically. Plan & the
adopt-vs-defer reasoning: `docs/harness/00-plan.md`.

## PR status
- **#21 harness scaffold** — merged.
- **#22 remove uv.lock (pyproject = source of truth)** — merged.
- **#23 determinism/security hooks** — OPEN. `.claude/hooks/` (block-rm-rf,
  block-push-main, injection-guard, block-env-edits) + `.githooks/pre-push`
  (gitleaks + main block) + `.agents/rules/security.md`.
- **#24 structured learning loop + user-memory tier** — OPEN.
  `.agents/learnings/{ERRORS,LEARNINGS,REJECTED}.md` + private `user-*.md` tier.

## Next
1. Review/merge #23 and #24 into `develop` (independent; either order).
2. Per-clone setup: `git config core.hooksPath .githooks`; `brew install gitleaks`.
3. Verify the harness in a fresh session in **both** Claude Code and Antigravity
   (symlink traversal, memory recall, hooks firing).
4. Resume the engineering backlog: `CODE_REVIEW_IMPLEMENTATION_PLAN.md` (45 tickets,
   M1–M6; start M1 — count-correctness: CR-1, CR-3, CR-2, CR-4).

## Key decisions
- **Adopt-vs-defer = content-bound vs infra-bound.** Context/risk layers (memory,
  learning loop, rules, security hooks) adopted *because the mature codebase earned
  them*. Infra layers (agents, fleet, second brain, multi-provider, crons) held
  until a concrete trigger — logged `REJ-20260626-001`.
- **Skill-frontmatter discipline deferred** — needs per-skill analysis, not a blind
  `allowed-tools` sweep (would wedge knowledge skills).
- `.agents/` = neutral source of truth; `.claude/` mirrors it. `pyproject.toml` =
  dependency source of truth (no lockfile).
