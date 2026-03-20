---
description: gbcms development rules — build, lint, testing, git workflow, CLI validation, mFSD rules, QC checklist
alwaysApply: true
---

# gbcms Development Guide

## Quick Start

```bash
git clone https://github.com/msk-access/gbcms.git && cd gbcms
python -m venv .venv && source .venv/bin/activate
maturin develop --release
gbcms --version
```

## Build Commands

```bash
maturin develop           # dev build (unoptimized, fast)
maturin develop --release # release build (use before integration tests)
maturin build --release --out dist  # wheel for distribution
```

## Full Lint Suite (All 4 Must Pass Before Every Commit)

```bash
source .venv/bin/activate
ruff check src/ tests/
black --check src/ tests/
mypy src/gbcms/ tests/
cargo clippy --manifest-path rust/Cargo.toml -- -D warnings
```

Fix in this order: ruff → black → mypy → clippy.

**pysam type patterns in test files:**
- Use `.cigartuples = [(op, len), ...]` (list of tuples) NOT `.cigar = (...)`
- Add `# type: ignore[assignment]` to `.query_qualities = [...]` lines
  (pysam stubs type it as `array[Any]` but lists work at runtime)

## Testing

```bash
pytest --no-cov -q          # fast (107+ tests)
pytest -v                   # verbose
pytest tests/test_mfsd_flag.py -v  # specific file
cargo test --manifest-path rust/Cargo.toml  # Rust inline tests
```

### Key Invariants (Assert in Every Counting Test)

```python
assert counts.dp >= counts.rd + counts.ad      # DP includes 'neither' reads
assert counts.dpf >= counts.rdf + counts.adf   # DPF includes discarded fragments
assert counts.rd == counts.rd_fwd + counts.rd_rev  # strand consistency
assert counts.ad == counts.ad_fwd + counts.ad_rev
```

### mFSD Column Count Assertions

```python
mfsd_cols = [c for c in maf.columns if c.startswith("mfsd_")]
assert len(mfsd_cols) == 34      # exact number of mFSD columns
assert len(no_mfsd_maf.columns) == 145   # without --mfsd
assert len(with_mfsd_maf.columns) == 179  # with --mfsd
```

## CLI Validation Rules (4 Layers — Enforce This Order)

1. **Parse-time (Typer)**: `Enum` for constrained choices, `min=`/`max=` for ranges. Typer rejects before Python code runs.
2. **Pre-model (cli.py command body)**: File-extension checks, cross-option dependencies. Log at `ERROR`, raise `typer.Exit(1)`.
3. **Model-time (Pydantic)**: `Field(ge=..., le=...)` and `@model_validator` in `models/core.py`.
4. **No silent skips**: Missing inputs → fail-fast or explicit opt-out (`--lenient-bam`).

The `cli.py` module docstring documents this order. Keep it in sync when adding options.

## mFSD Implementation Rules

- `--mfsd` gates output columns — writers check `self.mfsd`; absent when off, **not** NA-filled
- `--mfsd-parquet` requires `--mfsd` — validated at both CLI (pre-model) and Pydantic (`@model_validator`)
- `ref_sizes` and `alt_sizes` in `BaseCounts` are **Rust-internal** — no `#[pyo3(get)]`
- `write_fsd_parquet()` called from `pipeline.py` after `count_bam()`, not inside the counting engine
- Parquet uses ZSTD(1): `parquet = { default-features = false, features = ["arrow", "zstd"] }` in Cargo.toml — SNAPPY requires a separate `snap` feature

## Code Standards

### Python

- All public functions: type hints + Google-style docstrings
- `__all__` in every module
- Use `logging` — never `print()` or `console.print()` for status
- Pydantic models for all config — no raw dicts between layers
- `SimpleNamespace` only in test fixtures (`_zero_counts()`), never production

### Rust

- `///` doc comments on all public items
- `anyhow::Result` for internal errors; `pyo3::exceptions::PyIOError` for FFI boundary errors
- `log` crate for logging (`debug!`, `info!`, `warn!`)
- No `unwrap()` in production paths — use `?` or explicit error handling

## Git Workflow

```
main          ← stable releases (tagged)
  └── develop      ← integration branch (target for PRs)
       └── feature/*   ← individual features/fixes
```

- Branch from `develop`: `git checkout -b feature/my-thing develop`
- Merge with `--no-ff`; never push directly to `main`

## QC Checklist (Review Before Every Phase)

### Logging
- [ ] All status messages use proper log levels (debug/info/warning/error)
- [ ] No `print()` or `console.print()` for status
- [ ] Critical operations have timing/performance logging
- [ ] Errors include context (file path, sample name, variant position)

### Commenting
- [ ] All modules have module-level docstrings stating purpose and key design decisions
- [ ] All public functions have docstrings with Args/Returns
- [ ] Complex logic has inline comments explaining *why* (not what)
- [ ] No commented-out code

### No Code Duplication / Unused Code
- [ ] Common patterns extracted to helpers; no copy-paste between MafWriter and VcfWriter
- [ ] All imports are used (ruff)
- [ ] No dead functions or unreachable branches
- [ ] `__all__` exports match public API

### Tests
- [ ] New features have tests (new file or added to existing)
- [ ] All 4 counting invariants asserted
- [ ] `pytest --no-cov -q` passes (107+ tests)
- [ ] mFSD column counts verified if output.py changed

## Phase Implementation Workflow

```
1. Review affected files against QC checklist
2. Document findings in brain/ implementation plan
3. Fix QC issues FIRST before adding new features
4. Implement phase changes
5. Run full lint suite (ruff, black, mypy, clippy)
6. Run: pytest --no-cov -q
7. Run real-data integration test (for counting/mFSD/output changes)
8. Re-check QC standards
9. Commit with detailed message
```

## Real-Data Integration Test (22-BAM Suite)

For changes to counting engine, mFSD, shared/, or output writers:

```bash
maturin develop --release

# Build BAM list
for bam in ~/Downloads/bams_slice/P-*.bam; do
    name="$(basename "$bam" .bam)"; name="${name%_subset}"
    echo "$name $bam"
done > /tmp/gbcms_bam_list.txt

# Pass 1: no mFSD (regression baseline)
gbcms run \
  --variants ~/Downloads/bams_slice/filtered_indels.maf \
  --bam-list /tmp/gbcms_bam_list.txt \
  --fasta /path/to/Homo_sapiens_assembly19.fasta \
  --output-dir /tmp/gbcms_no_mfsd \
  --format maf --preserve-barcode --lenient-bam --threads 4

# Pass 2: mFSD + Parquet
gbcms run [same args] --mfsd --mfsd-parquet --output-dir /tmp/gbcms_with_mfsd

# Expected: 22 samples, 0 failed, exit 0
# No --mfsd: 145 columns, 0 .fsd.parquet files
# --mfsd:    179 columns (+34 mFSD), 22 .fsd.parquet files (62–208 KB each)
# Regression: alt_count/ref_count bit-for-bit identical to previous run
```

## Environment Variables

```bash
GBCMS_LOG_LEVEL=DEBUG RUST_LOG=debug gbcms run ...
```

## Release Process

1. `git checkout -b release/X.Y.Z develop`
2. Bump `version` in `pyproject.toml` and `rust/Cargo.toml`
3. Update `CHANGELOG.md`
4. PR to `main`, merge, tag: `git tag X.Y.Z && git push origin X.Y.Z`
5. CI: build → PyPI publish → Docker push to GHCR
6. Merge tag back to `develop`
