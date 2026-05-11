---
description: Automated code review checklist for a PR or diff
---

# /code-review

1. Identify files changed:
   ```bash
   git diff --name-only origin/develop..HEAD
   ```

2. For each changed file, review against:

   ### Architecture
   - [ ] Does this respect the Rust/Python boundary? (counting in Rust, orchestration in Python)
   - [ ] Is genomic binning used for BAM traversal? (not per-variant fetch)
   - [ ] Are type stubs synced? (`_rs.pyi` ↔ `gbcms_rs.pyi`)

   ### Code Quality
   - [ ] All functions have type hints and docstrings?
   - [ ] No `print()` — using `logging` module?
   - [ ] No hardcoded paths — using CLI args or config?
   - [ ] No `unwrap()` in Rust production paths?

   ### CLI Validation (4 Layers)
   - [ ] Parse-time: Typer `Enum`/`min`/`max` constraints?
   - [ ] Pre-model: file extension checks, cross-option deps?
   - [ ] Model-time: Pydantic `Field(ge=..., le=...)`, `@model_validator`?
   - [ ] No silent skips: fail-fast or explicit opt-out?

   ### Performance
   - [ ] No unnecessary copies?
   - [ ] Rayon parallelism per-bin, not per-variant?
   - [ ] COITree metadata access portable? (Borrow trait pattern)

   ### Testing
   - [ ] New features have tests?
   - [ ] All 4 counting invariants asserted? (DP >= RD+AD, strand consistency)
   - [ ] Edge cases covered?

   ### Security / Data
   - [ ] No PHI/PII logged or committed?
   - [ ] Sample IDs only — no patient identifiers?

3. Output findings as structured report with:
   - ✅ PASS / ⚠️ WARN / ❌ FAIL per category
   - Specific file:line references
   - Suggested fixes
