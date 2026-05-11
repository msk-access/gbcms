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

## Commit Discipline

- **Atomic commits**: One logical change per commit
- **Never commit**: `.venv/`, `__pycache__/`, `*.parquet` data files, `htmlcov/`
- **Run before every commit**: full lint suite (see code-quality.md)
- **Squash** trial-and-error commits before merging to `main`
