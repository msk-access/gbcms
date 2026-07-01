---
name: debug-error
description: Procedure for diagnosing gbcms build/test/CI errors — Rust compile failures, PyO3 stub mismatches, COITree SIMD (nosimd vs NEON/AVX) issues, CI black/ruff version traps, and test failures. Use when hitting an error and you need the gbcms-specific triage path.
---

# debug-error

1. **Identify error type**:
   - Rust compilation: `cd rust && cargo check 2>&1`
   - Python runtime: full traceback
   - CI failure: GitHub Actions log
   - Type error: `mypy src/`

2. **Rust compilation**:
   ```bash
   cd rust && cargo check 2>&1
   cargo check --target x86_64-unknown-linux-gnu 2>&1  # simulate CI
   ```
   - COITree metadata: use the `Borrow` trait for portable access (see rust-python-ffi)
   - PyO3 errors: check `_rs.pyi` matches Rust `#[pyo3(get)]`
   - OpenSSL: set `OPENSSL_NO_VENDOR=1` on Linux

3. **CI-specific**:
   - Black version mismatch: `target-version` in `pyproject.toml` must be `py310`
   - Ruff UP037: no quoted types in `.pyi` stubs
   - Rust nosimd vs NEON/AVX: test field-access patterns on both

4. **Test failures**:
   ```bash
   python -m pytest tests/<file>.py -v --tb=long
   cd rust && cargo test <name> -- --nocapture
   ```

5. **Environment**:
   ```bash
   mamba activate pyhbcms
   maturin develop --release
   gbcms --version
   ```
