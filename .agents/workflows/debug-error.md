---
description: Debug workflow for gbcms errors
---

# /debug-error

1. **Identify error type**:
   - Rust compilation: check `cargo check 2>&1`
   - Python runtime: check full traceback
   - CI failure: check GitHub Actions log
   - Type error: run `mypy src/`

2. **Rust compilation errors**:
   ```bash
   cd rust && cargo check 2>&1
   # For platform-specific issues:
   cargo check --target x86_64-unknown-linux-gnu 2>&1  # simulate CI
   ```
   - COITree metadata: use `Borrow` trait for portable access
   - PyO3 errors: check `_rs.pyi` stub matches Rust `#[pyo3(get)]`
   - OpenSSL: set `OPENSSL_NO_VENDOR=1` on Linux

3. **CI-specific failures**:
   - Black version mismatch: check `target-version` in `pyproject.toml` (must be `py310`)
   - Ruff UP037: no quoted types in `.pyi` stubs
   - Rust nosimd vs NEON/AVX: test field access patterns

4. **Test failures**:
   ```bash
   python -m pytest tests/<failing_file>.py -v --tb=long
   cd rust && cargo test <test_name> -- --nocapture
   ```

5. **Environment issues**:
   ```bash
   # Use pyhbcms mamba env for development
   mamba activate pyhbcms
   maturin develop --release
   gbcms --version
   ```
