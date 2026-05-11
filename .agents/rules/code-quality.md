---
description: Code quality standards, review discipline, and QC checklist for gbcms
alwaysApply: true
---

# Code Quality Standards

## ⚠️ Step Review Discipline (MANDATORY)

**BEFORE implementing each step:**
- [ ] Review existing code for reuse — do NOT re-implement what already exists
- [ ] Identify shared patterns that should be extracted to helpers
- [ ] Check for functions that can be generalized instead of duplicated

**AFTER implementing each step:**
- [ ] Audit for code duplication — extract shared logic
- [ ] Audit for unused code — delete dead functions, imports, variables
- [ ] Audit for silent failures — every `except` must log or re-raise
- [ ] Update commenting — all new functions have docstrings, complex logic has "why" comments
- [ ] Update logging — all significant operations have structured log entries
- [ ] Update monitoring — timing/counts for data operations, progress for long loops

**AFTER implementing, compare code against implementation plan:**
- [ ] Re-read the relevant plan section(s)
- [ ] Diff what was implemented vs. what the plan specified
- [ ] If there are **gaps**: classify as intentional or unintentional, notify user
- [ ] If the implementation **diverged**: update plan to match reality

## Lint Suite — All Must Pass Before Every Commit

```bash
# Python (in order)
ruff check src/ tests/
black --check src/ tests/
mypy src/

# Rust
cd rust && cargo clippy --all-targets -- -D warnings && cargo test && cd ..
```

Fix in this order: ruff → black → mypy → clippy.

## Python Standards

- All public functions: type hints + Google-style docstrings
- `__all__` in every module
- Use `logging` — never `print()` or `console.print()` for status
- Pydantic models for all config — no raw dicts between layers
- `SimpleNamespace` only in test fixtures, never production

### pysam Type Patterns (tests)
- Use `.cigartuples = [(op, len), ...]` NOT `.cigar = (...)`
- Add `# type: ignore[assignment]` to `.query_qualities = [...]`

## Rust Standards

- `///` doc comments on all public items
- `anyhow::Result` for internal errors; `pyo3::exceptions::PyIOError` for FFI boundary
- `log` crate for logging (`debug!`, `info!`, `warn!`)
- No `unwrap()` in production paths — use `?` or explicit error handling
- Platform-portable code: test with both nosimd and SIMD backends

## No Duplication

- [ ] Common patterns extracted to helpers
- [ ] No copy-paste between MafWriter and VcfWriter
- [ ] All imports used (ruff enforces)
- [ ] No dead functions or unreachable branches
- [ ] `__all__` exports match public API

## No Silent Failures

- [ ] Every `try/except` logs or re-raises
- [ ] Every `if not data: return` logs what was skipped
- [ ] Empty results always logged
- [ ] Missing files logged with path context
- [ ] Errors include context (file path, sample name, variant position)

## Test Invariants (Assert in Every Counting Test)

```python
assert counts.dp >= counts.rd + counts.ad      # DP includes 'neither'
assert counts.dpf >= counts.rdf + counts.adf   # DPF includes discarded fragments
assert counts.rd == counts.rd_fwd + counts.rd_rev  # strand consistency
assert counts.ad == counts.ad_fwd + counts.ad_rev
```

## Test Suite: 255 Python + 143 Rust Tests

| File | What It Tests |
|:-----|:-------------|
| `test_accuracy.py` | SNP, insertion, deletion, complex, MNP accuracy; DP invariant |
| `test_shifted_indels.py` | Windowed indel detection (15 cases) |
| `test_fuzzy_complex.py` | Quality-aware masked complex matching (14+ cases) |
| `test_filters.py` | Read filter flags |
| `test_strand_counts.py` | Strand-specific counting |
| `test_alignment_backend.py` | SW vs PairHMM concordance |
| `test_fragment_consensus.py` | Fragment quality consensus, DPF invariant |
| `test_multi_allelic.py` | Sibling ALT exclusion |
| `test_dp_neither.py` | DP includes neither/third-allele reads |
| `test_normalization.py` | Left-alignment, REF validation, window expansion |
| `test_cli_sample_id.py` | CLI argument parsing, validation, error paths |
| `test_cli_dna_rna.py` | DNA/RNA CLI integration, GTF parameter wiring |
| `test_maf_reader.py` | MAF input parsing |
| `test_maf_preservation.py` | MAF column preservation |
| `test_pipeline_v2.py` | End-to-end integration |
| `test_phase2_output.py` | Phase 2 output column schema |
| `test_mfsd_flag.py` | mFSD flag gating, column counts, Parquet output |
| `test_rescue_mnp.py` | MNP rescue pass validation |
| `test_config_isolation.py` | Config validation, threshold ranges |
| `test_diagnostic_flags.py` | gbcms_diagnostic flag computation |
| `test_rna_output.py` | RNA output column schema |
