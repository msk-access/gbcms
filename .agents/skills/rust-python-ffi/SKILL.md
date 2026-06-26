---
name: rust-python-ffi
description: Reference for the gbcms Rust/Python (PyO3) boundary — maturin build commands, #[pyo3(get)] ↔ .pyi stub synchronization, BaseCounts copy-on-write, platform-portable COITree (nosimd vs NEON/AVX), and build env vars. Use when editing PyO3 bindings, type stubs, or fixing platform/CI build issues.
---

# Rust/Python FFI Patterns

## Building

```bash
maturin develop           # dev build (fast, unoptimized)
maturin develop --release # release build (use before integration tests)
```

## PyO3 Bindings

- Rust structs with `#[pyclass]` expose fields via `#[pyo3(get)]` or `#[pyo3(get, set)]`
- Type stubs (`_rs.pyi`) must exactly match `#[pyo3(get)]` fields
- `BaseCounts` is frozen from Python — mutations use `with_ad()` copy-on-write
- `src/gbcms/_rs.pyi` is authoritative; `src/gbcms_rs.pyi` mirrors it — keep both synced

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

## Environment Variables

```bash
GBCMS_LOG_LEVEL=DEBUG RUST_LOG=debug gbcms run ...
OPENSSL_NO_VENDOR=1  # required for Linux CI builds
```
