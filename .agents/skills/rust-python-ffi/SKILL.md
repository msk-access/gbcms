---
name: rust-python-ffi
description: "Build and boundary reference for gbcms (Rust + Python via PyO3/maturin): the build recipe (repo root only, VIRTUAL_ENV set, CONDA_PREFIX unset), the stale .so trap, per-worktree venvs, #[pyo3(get)] and .pyi stub synchronization, BaseCounts being get-only from Python, COITree nosimd vs NEON, and triage when a build, cargo test, pytest or CI job fails. Use when building or rebuilding the extension, editing PyO3 bindings or type stubs, or when a Rust compile, test or CI check fails or a fix has no effect."
---

# Rust/Python FFI: build, boundary and triage

## Build & environment

`pyproject.toml` `[tool.maturin]` sets `manifest-path`, `python-source = "src"` and
`module-name = "gbcms._rs"`, so the build only works as configured from the repo root.

```bash
env -u CONDA_PREFIX VIRTUAL_ENV=<venv> PATH=<venv>/bin:$PATH maturin develop --release   # repo root only
maturin develop   # fast, unoptimized: for iterating; use --release for tests and measurements
```

- **Repo root only.** `maturin develop -m rust/Cargo.toml` ignores `[tool.maturin]`, installs a
  useless top-level `_rs` package, and never refreshes `src/gbcms/_rs.cpython-*.so`, which the
  editable install imports first: the build "succeeds" and tests keep running the old engine.
- **Stale `.so` trap.** After every rebuild check freshness:
  `python -c "import gbcms._rs as m, os, time; print(m.__file__, time.ctime(os.path.getmtime(m.__file__)))"`
  must print `src/gbcms/_rs...so` with a just-now mtime. If an old `.so` shadows it, `trash` it
  and rebuild from the root. Never filter build output (`| grep error`) without this check:
  a missing `maturin` binary also looks like silence (the main checkout's `.venv` has none).
- **Set `VIRTUAL_ENV`, unset `CONDA_PREFIX`.** Without `VIRTUAL_ENV` maturin installs into the
  nearest `.venv` above the cwd; from a worktree under `.claude/worktrees/` that is the MAIN
  checkout's `.venv`, which then points at the worktree for every other session (put it in
  every review-agent prompt; if it happens, rerun the build from the main root and check
  `gbcms.pth`). With both set, maturin stops: "Both VIRTUAL_ENV and CONDA_PREFIX are set".
- **A worktree needs its own venv**: the main `.venv` imports gbcms from the main checkout, so
  tests run from a worktree would exercise the wrong code. Per worktree, with `$SP` the session
  scratchpad: `uv venv --python 3.10 $SP/venv`, then
  `VIRTUAL_ENV=$SP/venv uv pip install -r pyproject.toml --group dev maturin` (runtime deps and
  the PEP 735 `dev` group: `test`, `docs` and the lint tools, unpinned, so CI's versions), then
  the build line above from the worktree root, then
  `$SP/venv/bin/python -c "import gbcms; print(gbcms.__file__)"` must print the worktree path.
- **Scratchpad venvs rot after a few days** (macOS purges `/private/tmp`; `maturin develop` then
  fails with "externally-managed-environment" or "Failed to read metadata"). Rebuild, don't
  repair: a fresh `uv venv`, the same install and build.
- **The main checkout's `.venv` has no maturin**: from its root use
  `VIRTUAL_ENV=$PWD/.venv uvx maturin develop --release` (with `CONDA_PREFIX` unset). Frozen
  measurement builds under `~/test` are not affected by a rotted scratchpad venv.
- **A baseline build of another commit**: `git archive <sha> | tar -x` into the scratchpad with
  its own venv (no change to any git state).
- `mkdocs build --strict` in a worktree: the git-revision-date plugin warns (`.git` is a file),
  so build with a temporary config without that plugin.
- **Linux wheels** (release.yml, manylinux) build with `OPENSSL_NO_VENDOR=1`.
- **Logging** is controlled by flags: `--verbose` (DEBUG, Rust records forwarded through
  `pyo3-log`) and `--trace` (per-read). `RUST_LOG` has no effect.

## When a check fails

1. **Tests pass but the fix "has no effect"**: the stale `.so` or the wrong venv (above).
2. **Rust compile error**: `cd rust && cargo check 2>&1`. CI's Linux job builds with the nosimd
   COITree backend (below); a pattern that compiles on this Mac's NEON backend can fail there.
3. **PyO3 or stub error**: `python -m pytest tests/test_rs_stub_parity.py`.
4. **Black or ruff pass locally, fail in CI**: lint-tool version skew (qa-check skill, "Lint
   tools"). `ruff` rule UP037: no quoted types in `.pyi` stubs. Black `target-version` in
   `pyproject.toml` is `py310`.
5. **A test fails**: `python -m pytest tests/<file>.py -v --tb=long`;
   `cd rust && cargo test <name> -- --nocapture`. A failing counting test is a decision point,
   never edited away (AGENTS.md). For per-read evidence run the CLI with `--trace` and read the
   `read call` lines.

## PyO3 Bindings

- Rust structs with `#[pyclass]` expose fields via `#[pyo3(get)]` or `#[pyo3(get, set)]`
- `src/gbcms/_rs.pyi` is the single stub for `gbcms._rs`; it must exactly match the
  `#[pyo3(get)]` fields and `#[pyo3(signature)]` params
  (`tests/test_rs_stub_parity.py` checks it)
- `BaseCounts` is get-only from Python (no setters, no mutators) — Python replaces a whole
  record rather than editing one (e.g. MNP rescue adopts the component SNV's `BaseCounts`)

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
