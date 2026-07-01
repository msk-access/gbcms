---
name: add-feature
description: Procedure for adding a feature to gbcms following its Git Flow and Rust/Python architecture (branch from develop, plan the layer, implement, sync type stubs, add tests with the 4 counting invariants, update docs/CHANGELOG, QA, PR to develop). Use when starting new feature work.
---

# add-feature

1. **Create feature branch**:
   ```bash
   git checkout develop && git pull origin develop
   git checkout -b feature/<scope>-<description>
   ```

2. **Plan the feature**:
   - Identify the layer: Rust counting, Python orchestration, CLI, or docs
   - Check `.agents/rules/architecture.md` for the Rust/Python boundary
   - Check existing code for reusable patterns (don't re-implement)

3. **Implement**:
   - Rust: `rust/src/` → rebuild with `maturin develop --release`
   - Python: `src/gbcms/`
   - CLI option: add to `cli.py`, validate in `models/core.py`
   - Type stubs: update `src/gbcms/_rs.pyi` (the single stub for `gbcms._rs`)
   - If it changes a CLI default/filter: update `nextflow/nextflow.config` too

4. **Add tests**: assert all 4 counting invariants; cover edge cases.

5. **Update docs**: `docs/cli/*`, `docs/nextflow/parameters.md` (if NF param),
   `docs/reference/output-formats.md` (if schema changed), CHANGELOG `## [Unreleased]`.

6. **Run QA** (see the qa-check skill).

7. **PR to develop**:
   ```bash
   git push origin feature/<scope>-<description>
   gh pr create --base develop --title "feat(<scope>): <description>"
   ```
