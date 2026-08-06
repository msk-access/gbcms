# Release Guide

This guide documents the complete release process for gbcms using git-flow workflow.

## Pre-Release Checklist

Before starting a release, ensure:

- [ ] All CI checks pass on `develop`
- [ ] All features for the release are merged to `develop`
- [ ] No blocking issues in milestone

---

## Version Locations

All these files must be updated with the new version (10 references total):

| File | Line | Format |
|:-----|:-----|:-------|
| `pyproject.toml` | 3 | `version = "X.Y.Z"` |
| `src/gbcms/__init__.py` | 11 | `__version__ = "X.Y.Z"` |
| `rust/Cargo.toml` | 3 | `version = "X.Y.Z"` |
| `nextflow/modules/local/gbcms/dna/main.nf` | 7 | `container "ghcr.io/msk-access/gbcms:X.Y.Z"` |
| `nextflow/modules/local/gbcms/build_gtf_cache/main.nf` | 4 | `container "ghcr.io/msk-access/gbcms:X.Y.Z"` |
| `nextflow/modules/local/gbcms/rna/main.nf` | 7 | `container "ghcr.io/msk-access/gbcms:X.Y.Z"` |
| `nextflow/modules/local/gbcms/normalize/main.nf` | 18 | `container "ghcr.io/msk-access/gbcms:X.Y.Z"` |
| `nextflow/modules/local/gbcms/merge/main.nf` | 7 | `container "ghcr.io/msk-access/gbcms:X.Y.Z"` |
| `nextflow/main.nf` | 53 | `gbcms vX.Y.Z — Nextflow Pipeline` |
| `nextflow/nextflow.config` | manifest | `version = 'X.Y.Z'` |
| `CHANGELOG.md` | Top section | `## [X.Y.Z] - YYYY-MM-DD` (new entry) |

!!! tip "Doc versions are now templated"
    Installation, quickstart, troubleshooting, and developer-guide docs use generic `X.Y.Z` notation. **No doc version bumps needed during release.**

!!! tip "Verify all references"
    After updating, run this to ensure no stale versions remain:
    ```bash
    grep -rn "OLD_VERSION" --include="*.py" --include="*.toml" --include="*.nf" --include="*.md" . \
      | grep -v ".git/" | grep -v "site/" | grep -v "CHANGELOG"
    ```

---

## Release Workflow

```mermaid
gitGraph LR:
   commit id: "ongoing develop work"
   branch release/X.Y.Z
   commit id: "bump versions (5 files)"
   commit id: "update CHANGELOG.md"
   checkout main
   merge release/X.Y.Z id: "PR merged" tag: "X.Y.Z"
   checkout develop
   merge release/X.Y.Z id: "back-merge"
```

!!! danger "Tags are bare `X.Y.Z` — NO `v` prefix"
    The `Release` workflow (`.github/workflows/release.yml`) triggers on the tag pattern
    `[0-9]+.[0-9]+.[0-9]+`. A `v`-prefixed tag (`v6.0.0`) **does not match** and will
    **silently fail to publish** — no PyPI, no Docker/GHCR, no docs deploy. Every existing
    release tag is bare (`5.3.0`, `5.2.0`, …); keep it that way. The `v` you see in
    `nextflow/main.nf`'s banner (`gbcms v6.0.0 — …`) is display text only, not the tag.

!!! info "Tag triggers CI"
    Pushing the bare tag `X.Y.Z` automatically triggers the CI pipeline which publishes to
    **PyPI**, **Docker/GHCR**, and deploys **gh-pages** docs.

---

## Step-by-Step Instructions

### 1. Create Release Branch

```bash
# From develop
git checkout develop
git pull origin develop

# Create release branch
git checkout -b release/X.Y.Z
```

### 2. Update Version Numbers

Update all version locations listed above. Use this command to verify:

```bash
# Check current versions
grep -E "^version|^__version__" pyproject.toml src/gbcms/__init__.py rust/Cargo.toml
grep "container\|gbcms v" nextflow/modules/local/gbcms/*/main.nf nextflow/main.nf
```

### 3. Update CHANGELOG.md

Add new section at top:

```markdown
## [X.Y.Z] - YYYY-MM-DD

### ✨ Added
- New feature description

### 🔧 Fixed
- Bug fix description

### 🔄 Changed
- Changes description
```

### 4. Run Pre-Release Checks

```bash
# Python linting + type checking
ruff check src/ tests/
black --check src/ tests/
mypy src/

# Rust linting + unit tests
cd rust && cargo clippy --all-targets -- -D warnings && cargo test && cd ..

# Integration tests
pytest -v
```

### 5. Commit and Push

```bash
git add -A
git commit -m "chore: bump version to X.Y.Z"
git push origin release/X.Y.Z
```

### 6. Create PR: release/X.Y.Z → main

- Title: `Release X.Y.Z`
- Describe changes from CHANGELOG
- Wait for CI to pass

### 7. Merge to main (creates tag)

After PR approval:
- **Merge commit** (do NOT squash) to `main`
- **Create tag**: `git tag X.Y.Z && git push origin X.Y.Z`

!!! warning "Do NOT squash-merge release PRs"
    Always use a **regular merge commit** for release PRs. Squash merging
    rewrites all commits into a single new SHA, which breaks shared ancestry
    between `main` and `develop`. This causes merge conflicts on every
    changed file during the Step 10 back-merge. Regular merge preserves
    commit history and makes the back-merge conflict-free.

### 8. CI Release Pipeline

The tag triggers `.github/workflows/release.yml`:

1. **Build wheels** (Linux x86_64, aarch64; macOS x86_64, arm64; Windows)
2. **Publish to PyPI** (via maturin)
3. **Build Docker image** → push to `ghcr.io/msk-access/gbcms:X.Y.Z`
4. **Deploy docs** → GitHub Pages (versioned via `mike` as `X.Y.Z` / `stable`)

### 9. Create the GitHub Release

!!! danger "The workflow does NOT create the GitHub Release"
    `release.yml` only publishes to PyPI / Docker / docs. The **Releases page** entry
    (with notes and the **Latest** badge) is a *separate* object you must create by hand,
    or the Releases page will keep showing the *previous* version even though the new tag
    exists and the packages published.

Create it from the CHANGELOG section on the (already-pushed) bare tag and mark it latest:

```bash
# Extract the [X.Y.Z] section from CHANGELOG.md into notes.md, then:
gh release create X.Y.Z \
  --title "X.Y.Z — <short summary>" \
  --notes-file notes.md \
  --latest --verify-tag
```

Verify with `gh release list` — the new version should show **Latest**.

### 10. Merge main back to develop

```bash
git checkout develop
git pull origin develop
git merge main
git push origin develop
```

### 11. Cleanup

```bash
# Delete local release branch
git branch -d release/X.Y.Z

# Delete remote release branch (optional)
git push origin --delete release/X.Y.Z
```

---

## Hotfix Workflow

For critical production fixes:

```bash
# Create hotfix from main
git checkout main
git checkout -b hotfix/X.Y.Z

# Fix, commit, push
git add -A
git commit -m "fix: critical issue description"
git push origin hotfix/X.Y.Z

# PR to main, then merge back to develop
```

---

## Automation Scripts

### git-flow-helper.sh

Interactive helper for git-flow operations:

```bash
./git-flow-helper.sh
# Options:
# 1) Create feature branch
# 2) Create release branch
# 3) Show git status
# 4) Cleanup merged branches
```

### Makefile Targets

| Target | Description |
|:-------|:------------|
| `make lint` | Run `ruff check` and `mypy` (Python only) |
| `make format` | Run `black` and `ruff --fix` |
| `make test` | Run `pytest` |
| `make test-cov` | Run tests with coverage report |
| `make docker-build` | Build Docker image locally |

!!! note
    `make lint` covers Python only. Always run `cargo clippy --all-targets -- -D warnings` separately to catch Rust linting errors before releasing.

---

## CI Workflows

| Workflow | Trigger | Purpose |
|:---------|:--------|:--------|
| `test.yml` | Push to develop/main, PR | Run tests |
| `release.yml` | Tag push `X.Y.Z` | Build wheels, publish PyPI, Docker |
| `deploy-docs.yml` | Push to main or develop (docs/) | Deploy versioned docs via `mike` (`stable` from main, `dev` from develop) |

---

## Troubleshooting

### PyPI Upload Fails

- Check if version already exists on PyPI (versions cannot be overwritten)
- Verify `PYPI_TOKEN` secret is set in GitHub repository

### Docker Build Fails

- Check `Dockerfile` paths match the new folder structure
- Verify rust/Cargo.toml version matches

### Docs Build Fails

- Verify `mkdocs-mermaid2-plugin` is installed in workflow
- Check snippet paths are correct (relative to root)

---

## Related

- [Developer Guide](developer-guide.md) — Setup, build commands, and project layout
- [Contributing](contributing.md) — Contribution workflow and code standards
- [Testing Guide](testing-guide.md) — Running and writing tests before a release
- [Changelog](changelog.md) — Version history
