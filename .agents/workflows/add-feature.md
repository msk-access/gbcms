---
description: Add a new feature following gbcms architecture patterns
---

# /add-feature

1. **Create feature branch**:
   ```bash
   git checkout develop
   git pull origin develop
   git checkout -b feature/<scope>-<description>
   ```

2. **Plan the feature**:
   - Identify which layer: Rust counting, Python orchestration, CLI, or docs
   - Check architecture.md for Rust/Python boundary rules
   - Check existing code for reusable patterns

3. **Implement**:
   - Rust changes: `rust/src/` → rebuild with `maturin develop --release`
   - Python changes: `src/gbcms/`
   - CLI option: add to `cli.py`, validate in `models/core.py`
   - Type stubs: update `_rs.pyi` first, then sync `gbcms_rs.pyi`

4. **Add tests**:
   - New test file or extend existing
   - Assert all 4 counting invariants
   - Test edge cases

5. **Update docs**:
   - CLI docs: `docs/cli/dna.md` and/or `docs/cli/rna.md`
   - If Nextflow param: `docs/nextflow/parameters.md`
   - If output schema changed: `docs/reference/output-formats.md`
   - CHANGELOG: add entry under `## [Unreleased]`

6. **Run QA** (see qa-check.md)

7. **PR to develop**:
   ```bash
   git push origin feature/<scope>-<description>
   gh pr create --base develop --title "feat(<scope>): <description>"
   ```
