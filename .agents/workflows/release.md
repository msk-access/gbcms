---
description: PR-driven version bump, test, and publish workflow
---

# /release

This project strictly utilizes a PR-driven **Git Flow**. Direct merging to `main` is prohibited.

## Steps

1. **Initialize Release Branch**:
   ```bash
   git checkout develop
   git pull origin develop
   git checkout -b release/X.Y.Z
   ```

2. **Version Bump (9 locations)**:
   Update all these files from OLD → NEW version:

   | File | Format |
   |:-----|:-------|
   | `pyproject.toml` | `version = "X.Y.Z"` |
   | `src/gbcms/__init__.py` | `__version__ = "X.Y.Z"` |
   | `rust/Cargo.toml` | `version = "X.Y.Z"` |
   | `nextflow/nextflow.config` | `version = 'X.Y.Z'` |
   | `nextflow/main.nf` | `gbcms vX.Y.Z — Nextflow Pipeline` |
   | `nextflow/modules/local/gbcms/dna/main.nf` | `container "ghcr.io/msk-access/gbcms:X.Y.Z"` |
   | `nextflow/modules/local/gbcms/rna/main.nf` | `container "ghcr.io/msk-access/gbcms:X.Y.Z"` |
   | `nextflow/modules/local/gbcms/normalize/main.nf` | `container "ghcr.io/msk-access/gbcms:X.Y.Z"` |
   | `nextflow/modules/local/gbcms/merge/main.nf` | `container "ghcr.io/msk-access/gbcms:X.Y.Z"` |

3. **Verify No Stale Versions**:
   ```bash
   grep -rn "OLD_VERSION" --include="*.py" --include="*.toml" --include="*.nf" --include="*.lock" . \
     | grep -v ".git/" | grep -v "site/" | grep -v "CHANGELOG"
   ```

4. **Finalize CHANGELOG**:
   ```bash
   # Update date: ## [X.Y.Z] - YYYY-MM-DD
   sed -i '' "s/## \[X.Y.Z\] - Unreleased/## [X.Y.Z] - $(date +%Y-%m-%d)/" CHANGELOG.md
   ```

5. **Full Verification Matrix**:
   ```bash
   black --check src/ tests/
   ruff check src/ tests/
   mypy src/
   python -m pytest tests/ -q
   cd rust && cargo clippy -- -D warnings && cargo test && cd ..
   ```

6. **Commit and Push**:
   ```bash
   git add -A
   git commit -m "release: bump version to X.Y.Z"
   git push origin release/X.Y.Z
   ```

7. **Create PR**: `release/X.Y.Z` → `main`
   - Title: `Release X.Y.Z`
   - Body: copy from CHANGELOG
   - Wait for CI to pass

8. **After PR Merge — Tag**:
   ```bash
   git checkout main
   git pull origin main
   git tag X.Y.Z          # NO 'v' prefix — matches release.yml trigger
   git push origin X.Y.Z
   ```

9. **Back-merge**:
   ```bash
   git checkout develop
   git pull origin develop
   git merge main
   git push origin develop
   ```

10. **Cleanup**:
    ```bash
    git branch -d release/X.Y.Z
    git push origin --delete release/X.Y.Z
    ```

## CI Trigger

Tag `X.Y.Z` (no `v` prefix) triggers `.github/workflows/release.yml`:
- Build manylinux wheels via maturin-action
- Publish to PyPI via Trusted Publisher OIDC
- Build + push Docker image to `ghcr.io/msk-access/gbcms:X.Y.Z` + `:latest`
- Deploy docs to GitHub Pages
