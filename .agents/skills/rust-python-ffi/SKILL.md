---
name: rust-python-ffi
description: "Reference for the gbcms Rust/Python (PyO3) boundary and its build — maturin develop from the repo root, the stale .so trap, worktree venvs and build env vars, #[pyo3(get)] ↔ .pyi stub synchronization, BaseCounts being frozen from Python, and platform-portable COITree (nosimd vs NEON/AVX). Use when building or rebuilding the extension, editing PyO3 bindings or type stubs, or fixing platform/CI build issues. For triaging a failing test or CI job see debug-error."
---

# Rust/Python FFI Patterns

## Build & environment

`pyproject.toml` `[tool.maturin]` sets `manifest-path`, `python-source = "src"` and
`module-name = "gbcms._rs"`, so the build only works as configured from the repo root.

```bash
VIRTUAL_ENV=<venv> PATH=<venv>/bin:$PATH maturin develop --release   # repo root only
maturin develop   # fast, unoptimized: for iterating; use --release for tests and measurements
```

- **Repo root only.** `maturin develop -m rust/Cargo.toml` ignores `[tool.maturin]`, installs a
  useless top-level `_rs` package, and never refreshes `src/gbcms/_rs.cpython-*.so`, which the
  editable install imports first: the build "succeeds" and tests keep running the old engine.
- **Stale `.so` trap.** After every rebuild check freshness:
  `python -c "import gbcms._rs as m, os, time; print(m.__file__, time.ctime(os.path.getmtime(m.__file__)))"`
  must print `src/gbcms/_rs...so` with a just-now mtime. If an old `.so` shadows it, `trash` it
  and rebuild from the root. Never filter build output (`| grep error`) without this check:
  a missing `maturin` binary also looks like silence.
- **Always set `VIRTUAL_ENV`.** Without it maturin installs into the nearest `.venv` above the
  cwd; from a worktree under `.claude/worktrees/` that is the MAIN checkout's `.venv`, which
  then points at the worktree for every other session. Review-agent prompts say so too.
- **A worktree needs its own venv.** The main `.venv` imports gbcms from the main checkout, so
  tests run from a worktree would exercise the wrong code. Per worktree:
  `uv venv --python 3.10 $SP/venv`, then
  `VIRTUAL_ENV=$SP/venv uv pip install -r pyproject.toml --group dev maturin`, then the build
  line above from the worktree root, then
  `$SP/venv/bin/python -c "import gbcms; print(gbcms.__file__)"` must print the worktree path.
  Scratchpad venvs rot after a few days (macOS purges `/private/tmp`): rebuild, don't repair.
- **Linux wheels** (release.yml, manylinux) build with `OPENSSL_NO_VENDOR=1`.
- **Logging** is controlled by flags: `--verbose` (DEBUG, Rust records forwarded through
  `pyo3-log`) and `--trace` (per-read). `RUST_LOG` has no effect.

## PyO3 Bindings

- Rust structs with `#[pyclass]` expose fields via `#[pyo3(get)]` or `#[pyo3(get, set)]`
- Type stubs (`_rs.pyi`) must exactly match `#[pyo3(get)]` fields
  (`tests/test_rs_stub_parity.py` checks it)
- `BaseCounts` is frozen from Python (get-only fields, no mutators) — Python replaces a whole
  record rather than editing one (e.g. MNP rescue adopts the component SNV's `BaseCounts`)
- `src/gbcms/_rs.pyi` is the single stub for `gbcms._rs` — update it when fields/signatures change

## Platform-Portable Code

The `coitrees` crate uses different node types per SIMD backend:
- **nosimd** (CI Linux): `IntervalNode<usize, _>` → `node.metadata` is `usize`
- **NEON** (macOS ARM): `Interval<&usize>` → `node.metadata` is `&usize`

Pattern for portable metadata access:
```rust
use std::borrow::Borrow;
#[allow(noop_method_call)]
let idx: &usize = node.metadata.borrow();
let exon = &self.exons[*idx];
```
