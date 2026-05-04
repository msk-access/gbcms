# Release Guide

This guide documents the complete release process for gbcms using git-flow workflow.

## Pre-Release Checklist

Before starting a release, ensure:

- [ ] All CI checks pass on `develop`
- [ ] All features for the release are merged to `develop`
- [ ] No blocking issues in milestone

---

## Version Locations

All these files must be updated with the new version (7 references total):

| File | Line | Format |
|:-----|:-----|:-------|
| `pyproject.toml` | 3 | `version = "X.Y.Z"` |
| `src/gbcms/__init__.py` | 11 | `__version__ = "X.Y.Z"` |
| `rust/Cargo.toml` | 3 | `version = "X.Y.Z"` |
| `nextflow/modules/local/gbcms/dna/main.nf` | 7 | `container "ghcr.io/msk-access/gbcms:X.Y.Z"` |
| `nextflow/modules/local/gbcms/rna/main.nf` | 7 | `container "ghcr.io/msk-access/gbcms:X.Y.Z"` |
| `nextflow/modules/local/gbcms/normalize/main.nf` | 18 | `container "ghcr.io/msk-access/gbcms:X.Y.Z"` |
| `nextflow/main.nf` | 53 | `gbcms vX.Y.Z — Nextflow Pipeline` |
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
   merge release/X.Y.Z id: "PR merged" tag: "vX.Y.Z"
   checkout develop
   merge release/X.Y.Z id: "back-merge"
```

!!! info "Tag triggers CI"
    Pushing the tag `vX.Y.Z` automatically triggers the CI pipeline which publishes to **PyPI**, **Docker/GHCR**, and deploys **gh-pages** docs.

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
- Squash and merge to `main`
- **Create tag**: `git tag X.Y.Z && git push origin X.Y.Z`

### 8. CI Release Pipeline

The tag triggers `.github/workflows/release.yml`:

1. **Build wheels** (Linux x86_64, aarch64; macOS x86_64, arm64; Windows)
2. **Publish to PyPI** (via maturin)
3. **Build Docker image** → push to `ghcr.io/msk-access/gbcms:X.Y.Z`
4. **Deploy docs** → GitHub Pages (versioned via `mike` as `X.Y.Z` / `stable`)

### 9. Merge main back to develop

```bash
git checkout develop
git pull origin develop
git merge main
git push origin develop
```

### 10. Cleanup

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
