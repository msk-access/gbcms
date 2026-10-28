---
description: Git conventions for gbcms development
alwaysApply: true
---

# Git Conventions

## Conventional Commits (MANDATORY)

All commits MUST follow conventional commits format:

```
<type>(<scope>): <description>

[optional body]
```

| Type | Use For |
|---|---|
| `feat` | New feature (counting, annotation, CLI option) |
| `fix` | Bug fix |
| `refactor` | Code restructuring without behavior change |
| `docs` | Documentation only (mkdocs, docstrings) |
| `test` | Adding/fixing tests |
| `chore` | Build, CI, dependency updates |
| `style` | Formatting (black, ruff) |
| `release` | Version bumps and release prep |

**Scopes** (optional but encouraged):
`rust`, `cli`, `pipeline`, `io`, `models`, `normalize`, `rna`, `gtf`, `asjd`,
`mfsd`, `rescue`, `docs`, `ci`, `nf`, `stubs`

**Examples:**
- `feat(rna): add GTF-based transcript annotation`
- `fix(rust): portable COITree metadata access via Borrow trait`
- `docs(cli): add --rescue-mnp-threshold boundary semantics`
- `release: bump version to 5.0.0`

## Branch Naming

```
feature/<scope>-<description>    # feature/rna-v5
fix/<scope>-<description>        # fix/coitree-metadata
docs/<description>               # docs/rna-annotation
release/<version>                # release/5.0.0
hotfix/<version>                 # hotfix/5.0.1
```

## Git Flow

```
main          ← stable releases (tagged X.Y.Z)
  └── develop      ← integration branch (target for PRs)
       └── feature/*   ← individual features/fixes
```

- Branch from `develop`: `git checkout -b feature/my-thing develop`
- Merge with `--no-ff`; never push directly to `main`
- Releases go through PR: `release/X.Y.Z` → `main`

## Issues

An issue's comments are its working record; the repo is what agents read.

- **Close only after the outcome is in the repo.** The decision goes in the read-rule
  register (`docs/reference/read-judgment.md`) or a proposal doc, with the numbers that
  decided it. Rejected ideas go in `.agents/learnings/REJECTED.md`, and shipped changes in
  `CHANGELOG.md`. The closing comment links each place.
- **Search closed issues before a fix, proposal or decision**, alongside `REJECTED.md`:
  `gh search issues --repo msk-access/gbcms --state closed "<terms>"`. Their comments hold
  the measurements behind past decisions.
- **Group related work as sub-issues** (operator, 2026-09-25): a parent issue with one sub-issue
  per ticket, nested as needed, never one combined issue or a checklist. Each cycle gets a
  tracking issue with every ticket under it. `gh issue create --parent N` for a new issue;
  `gh issue edit N --parent P` or `--add-sub-issue` for an existing one (gh 2.9x or newer).

## Merging

The repository disallows GitHub auto-merge. Merge only when the operator asks (2026-10-26: "Use gh
to merge", "when I ask"): check the PR's checks once with `gh pr checks <n> --repo
msk-access/gbcms`; if they are green, run `gh pr merge <n> --repo msk-access/gbcms --merge` (a
merge commit, as for earlier PRs); if they are pending, say so and wait to be asked again. Never
poll CI from the session. Keep Auto-fix on, so failures, conflicts and review comments wake the
session. After a merge, run the retirement checklist (`docs/harness/measurement-harness.md`).

## Shared checkout hazards

- **Review and workflow subagents may change git state.** They run in the session's current
  directory, in the same working tree unless started with `isolation: 'worktree'`. One ran `git
  stash` and silently reverted uncommitted engine edits (found only because a later scripted
  edit's `assert` failed); another checked out the reviewed SHA, so the next commit landed on a
  detached HEAD and `git push` printed "Everything up-to-date" while the remote stayed behind.
  Commit (or tag your stash) before you launch a same-checkout workflow, or give agents that run
  git commands `isolation: 'worktree'`; launch review workflows from the main checkout, never from
  a worktree. After any workflow check `git status`, `git stash list` and `git branch
  --show-current` (empty means detached), and verify a push by comparing `git ls-remote` with
  `HEAD`, not by push's exit code.
- **Concurrent sessions write `.agents/memory/`** (the Claude recall directory is a symlink to
  it), so `git add -A` sweeps another session's in-flight memory edits into your commit. Stage
  explicit paths only; leave `.agents/memory/` changes you did not make unstaged.
- **Gate exit codes.** `cmd | tail && git commit` commits when `cmd` fails (qa-check skill).

## Commit Discipline

- **Atomic commits**: One logical change per commit
- **Never commit**: `.venv/`, `__pycache__/`, `*.parquet` data files, `htmlcov/`
- **Run before every commit**: the qa-check gate (all of it, docs-only changes too)
- **Squash** trial-and-error commits before merging to `main`
