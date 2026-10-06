# Release Guide

This guide documents the complete release process for gbcms using git-flow workflow.

## Pre-Release Checklist

Before starting a release, ensure:

- [ ] All CI checks pass on `develop`
- [ ] All features for the release are merged to `develop`
- [ ] No blocking issues in milestone

---

## Version Locations

One version source per ecosystem: **`rust/Cargo.toml`** for the Python package and the
extension (maturin converts `X.Y.Z-dev.N` to PEP 440 `X.Y.Z.devN`), and the **Nextflow
manifest** for the pipeline's container image and banner. A release edits three
files; the lock follows, and the rest is derived:

| Location | Role |
|:---------|:-----|
| `rust/Cargo.toml` | **edit**: `version` under `[package]`, the Python package and the extension |
| `rust/Cargo.lock` | **follows**: run `cargo check` after editing Cargo.toml, commit the lock |
| `nextflow/nextflow.config` | **edit**: `manifest.version`, the pipeline's image tag and banner |
| `CHANGELOG.md` | **edit**: the dated section `## [X.Y.Z] - YYYY-MM-DD — summary` (the summary becomes the GitHub Release title) |
| `pyproject.toml` | derived: `dynamic = ["version"]`, maturin reads Cargo.toml |
| `src/gbcms/__init__.py` | derived: `__version__` from the installed package metadata |
| `nextflow/modules/local/gbcms/*/main.nf` | derived: `container "ghcr.io/msk-access/gbcms:${workflow.manifest.version}"` |
| `nextflow/main.nf` | derived: the banner reads `workflow.manifest.version` |

On `develop` the manifest names the last **released** image (images publish only on
tags), so it trails the package version; at a release both equal `X.Y.Z`.

!!! tip "Check the sources"
    `python scripts/release.py check` runs on every PR (the derived locations still read
    their source, the lock matches Cargo.toml, the manifest is a released version no
    newer than the package). With `--tag X.Y.Z` it is the release workflow's first job:
    the tag must be a bare `X.Y.Z` and Cargo.toml, the lock and the manifest must all
    equal it, with a dated CHANGELOG section, or nothing builds or publishes.

    ```bash
    python scripts/release.py check --tag X.Y.Z
    ```

---

## Release Workflow

```mermaid
gitGraph LR:
   commit id: "ongoing develop work"
   branch release/X.Y.Z
   commit id: "bump the version (Cargo.toml, manifest)"
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
    `nextflow/main.nf`'s banner (`gbcms vX.Y.Z — …`, from the manifest) is display text only, not the tag.

!!! info "Tag triggers CI"
    Pushing the bare tag `X.Y.Z` triggers `release.yml`, which verifies the version
    sources and publishes to **PyPI** and **Docker/GHCR** and creates the **GitHub Release**.
    The docs deploy separately (`deploy-docs.yml`, on pushes to `main`/`develop` that touch
    `docs/`). Only a tag push publishes: a manual run builds and checks but publishes nothing.

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

### 2. Update the Version

Edit `version` in `rust/Cargo.toml` and `manifest.version` in
`nextflow/nextflow.config` to `X.Y.Z`, refresh the lock, and check:

```bash
cd rust && cargo check && cd ..        # updates the gbcms_rs entry in Cargo.lock
python scripts/release.py check --tag X.Y.Z
```

### 3. Update CHANGELOG.md

Add new section at top. The text after the em dash becomes the GitHub Release title
(`X.Y.Z — summary`), and the section body its notes:

```markdown
## [X.Y.Z] - YYYY-MM-DD — short summary

### ✨ Added
- New feature description

### 🔧 Fixed
- Bug fix description

### 🔄 Changed
- Changes description
```

### 3b. Remove cycle plans

Cycle plans (`CYCLE_X.Y.Z_PLAN.md`) live on `develop` and never ship. Delete
every one on the release branch, and point `CONTINUITY.md` elsewhere if it
references them:

```bash
git rm CYCLE_*_PLAN.md
```

The back-merge (step 10) then removes the finished cycle's plan from `develop`
too. The next cycle's plan is added to `develop` after the back-merge.

### 4. Run Pre-Release Checks

Refresh the Docker image's dependency lock first: the newest releases that satisfy
`pyproject.toml`. Run `latest-deps.yml` by hand first (Actions → Latest dependencies →
Run workflow): it tests those releases on every Python, and its summary shows what the
refresh changes. Commit the lock with the release.

```bash
uv pip compile pyproject.toml --upgrade --generate-hashes --python-version 3.11 --python-platform x86_64-manylinux_2_28 -o docker/requirements.lock
```

The image installs exactly these versions (`--require-hashes`) and runs `pip check`,
so a lock that no longer satisfies `pyproject.toml` fails the Docker build in CI.

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

!!! warning "Confirm the checks actually appeared — absent is not the same as passing"
    Opening the PR immediately after `git push` can race: GitHub occasionally fails to
    dispatch any workflow for the `pull_request` event, and the PR then shows only the
    GitBook statuses. Nothing is marked failed or pending — the checks are simply **not
    there**, which reads like "nothing to run" rather than "nothing ran". This happened on
    6.2.0; 6.1.0 got the full suite from the same steps, so it is intermittent, not config.

    ```bash
    gh pr checks <PR>            # expect Tests jobs + Nextflow Lint, not just GitBook
    gh run list --branch release/X.Y.Z
    ```

    If they are missing, re-fire rather than assuming:

    ```bash
    # verify the exact release SHA (no PR churn, no notifications)
    gh workflow run test.yml --ref release/X.Y.Z
    gh workflow run nextflow-lint.yml --ref release/X.Y.Z

    # or re-fire the pull_request event so checks attach to the PR itself
    gh pr close <PR> && gh pr reopen <PR>
    ```

    Pausing a beat between `git push` and PR creation makes the race far less likely.

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

0. **Verify** — `scripts/release.py check --tag X.Y.Z`; every other job waits for it
1. **Build one wheel** — `cp311`, `manylinux_2_34_x86_64` — plus an **sdist** (jobs `linux`
   and `sdist`)
2. **Publish to PyPI** (via maturin)
3. **Build Docker image** → push to `ghcr.io/msk-access/gbcms:X.Y.Z`
4. **Create the GitHub Release** — title and notes from the CHANGELOG section, the
   wheel and sdist attached with `SHA256SUMS`, and a build-provenance attestation for
   each artifact (verify with `gh attestation verify <file> --repo msk-access/gbcms`)

The docs are not part of this workflow: `deploy-docs.yml` publishes them (via `mike`, as
`X.Y.Z` / `stable`) when the release merge reaches `main`.

!!! warning "Re-run only the failed jobs"
    "Re-run all jobs" rebuilds the wheel, which is not bit-reproducible: PyPI keeps the
    file it already has (`skip-existing`), while the release page would get the new
    build, so its checksums and attestations would no longer match what PyPI serves.

!!! warning "One wheel, not a matrix"
    This list previously claimed Linux x86_64 + aarch64, macOS x86_64 + arm64, and Windows.
    `release.yml` has only `linux` and `sdist` jobs, and 6.2.0 published exactly
    `gbcms-X.Y.Z-cp311-cp311-manylinux_2_34_x86_64.whl` and `gbcms-X.Y.Z.tar.gz`.

    The practical consequence, worth knowing before telling a user to `pip install gbcms`:
    **everyone not on cp311 manylinux x86_64 builds from the sdist**, which needs a Rust
    toolchain. That includes macOS (Intel and Apple Silicon) and every Python other than
    3.11. Broadening the matrix is a real change to `release.yml`, not a docs fix — until
    then, this list should describe what actually ships.

### 9. Check the GitHub Release

The `github-release` job creates the Releases page entry once PyPI and the image are
published (or updates it, re-uploading the assets, if one already exists). Verify with
`gh release list`: the new version should show **Latest**. If the job failed, re-run it
from the Actions tab (failed jobs only). By hand, the page alone (without the assets,
checksums and attestations, which only the job adds) is:

```bash
python scripts/release.py notes --tag X.Y.Z --title-out title.txt --notes-out notes.md
gh release create X.Y.Z --title "$(cat title.txt)" --notes-file notes.md --latest --verify-tag
```

### 10. Merge main back to develop

```bash
git checkout develop
git pull origin develop
git merge main
git push origin develop
```

### 10b. Move develop to the next dev version

Right after the back-merge, bump `rust/Cargo.toml` on `develop` to the next dev version
(`X.Y+1.0-dev.0`), run `cargo check`, and add an empty `## [Unreleased]` section. Leave
the Nextflow manifest at `X.Y.Z`: it names the newest published image. Without the bump
every dev build stamps the released version in its provenance line; the version check
on `develop` fails until it is done.

```bash
# rust/Cargo.toml: version = "X.Y+1.0-dev.0"
cd rust && cargo check && cd ..
python scripts/release.py check --branch develop
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

# Fix, commit; then set the version to X.Y.Z in rust/Cargo.toml (cargo check) and the
# Nextflow manifest, add a dated CHANGELOG section, and check
python scripts/release.py check --tag X.Y.Z
git add -A
git commit -m "fix: critical issue description"
git push origin hotfix/X.Y.Z

# PR to main, tag X.Y.Z, then merge back to develop (keeping develop's dev version)
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
| `test.yml` | Push to develop/main, PR | Run tests: Ubuntu 3.10 (dependencies at their floors), 3.11 and 3.14, macOS 3.12; lint; Rust unit tests; build and run the image |
| `latest-deps.yml` | Monthly (the 1st), manual before a release | The suite on the newest dependency releases, Python 3.10–3.14, with compatible Rust updates; a failure opens or updates one `latest-deps` issue |
| `release.yml` | Tag push `X.Y.Z` | Verify the version sources, build the wheel and sdist, publish PyPI and Docker, create the GitHub Release |
| `nextflow-lint.yml` | Push/PR touching `nextflow/` | Strict-syntax lint; every process `nextflow inspect` resolves runs the manifest's image |
| `deploy-docs.yml` | Push to main or develop (docs/) | Deploy versioned docs via `mike` (`stable` from main, `dev` from develop) |

---

## Troubleshooting

### PyPI Upload Fails

- Check if version already exists on PyPI (versions cannot be overwritten)
- Verify `PYPI_TOKEN` secret is set in GitHub repository

### Docker Build Fails

- Check `Dockerfile` paths match the new folder structure
- Run `python scripts/release.py check --tag X.Y.Z` (the release workflow's first job)

### Docs Build Fails

- Verify `mkdocs-mermaid2-plugin` is installed in workflow
- Check snippet paths are correct (relative to root)

---

## Related

- [Developer Guide](developer-guide.md) — Setup, build commands, and project layout
- [Contributing](contributing.md) — Contribution workflow and code standards
- [Testing Guide](testing-guide.md) — Running and writing tests before a release
- [Changelog](changelog.md) — Version history
