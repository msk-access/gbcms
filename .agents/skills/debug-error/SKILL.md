---
name: debug-error
description: "Triage procedure for a failing gbcms check — which log to read and which gbcms-specific cause to test first for a cargo compile error, a failing pytest or cargo test, or a red CI job (black/ruff version skew, stub mismatch, nosimd vs NEON). Use when something fails and you need the triage path. For how to build, rebuild or set up the venv and worktree environment see rust-python-ffi."
---

# debug-error

Triage only. Build, venv and stale-`.so` problems: see rust-python-ffi ("Build & environment");
when tests pass locally but a fix "has no effect", check that first.

1. **Identify the error type**:
   - Rust compilation: `cd rust && cargo check 2>&1`
   - Python runtime: the full traceback
   - CI failure: the GitHub Actions log of the failing job
   - Type error: `mypy src/`

2. **Rust compilation**:
   ```bash
   cd rust && cargo check 2>&1
   cargo check --target x86_64-unknown-linux-gnu 2>&1  # simulate CI
   ```
   - COITree metadata: use the `Borrow` trait for portable access (see rust-python-ffi)
   - PyO3 errors: check `_rs.pyi` matches Rust `#[pyo3(get)]`
     (`python -m pytest tests/test_rs_stub_parity.py`)

3. **CI-specific**:
   - Black/ruff pass locally but fail in CI: version skew. Run CI's versions (latest, from the
     `dev` group in `pyproject.toml`), not the venv's older ones, before trusting a clean result.
   - Black `target-version` in `pyproject.toml` must be `py310`
   - Ruff UP037: no quoted types in `.pyi` stubs
   - Rust nosimd vs NEON/AVX: test field-access patterns on both

4. **Test failures**:
   ```bash
   python -m pytest tests/<file>.py -v --tb=long
   cd rust && cargo test <name> -- --nocapture
   ```
   A failing counting test is a decision point, not something to edit away (AGENTS.md).
   Per-read evidence: run the CLI with `--trace` and read the `read call` lines.
