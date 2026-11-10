---
name: add-feature
description: "Steps to add a feature: branch from develop, pick the layer (Rust or Python), implement, sync the type stub, test with the counting invariants, update docs and CHANGELOG, run QA, PR to develop. Use when starting new feature work."
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
   - Rust: `rust/src/` → rebuild with `maturin develop --release` from the repo root, with
     `VIRTUAL_ENV` set (rust-python-ffi, "Build & environment")
   - Python: `src/gbcms/`
   - CLI option: add to `cli.py`, validate in `models/core.py`
   - Type stubs: update `src/gbcms/_rs.pyi` (the single stub for `gbcms._rs`)
   - If it changes a CLI default/filter: update `nextflow/nextflow.config` too

4. **Add tests**: assert all 4 counting invariants; cover edge cases.

5. **Update docs**: `docs/cli/*`, `docs/nextflow/parameters.md` (if NF param),
   `docs/reference/output-formats.md` (if schema changed), CHANGELOG `## [Unreleased]`.

6. **Run QA** (see the qa-check skill).

   **Real-data acceptance** for count-changing work follows
   `.agents/rules/decision-protocol.md` (every data type, more than one truth source, DNA and RNA
   judged separately, scored against the reads' own bases). Build base and candidate in isolation
   (`git archive <sha>` and its own venv, the commit recorded beside each build; rust-python-ffi
   "Build & environment"; layout in `docs/harness/measurement-harness.md`). After the PR merges, run
   the retirement checklist there (delete the builds, gzip raw `--trace` logs).

7. **PR to develop**:
   ```bash
   git push origin feature/<scope>-<description>
   gh pr create --base develop --title "feat(<scope>): <description>"
   ```
