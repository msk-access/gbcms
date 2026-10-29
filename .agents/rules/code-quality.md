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
- [ ] Reproduce every claim before stating it (a review finding, an "unreachable" branch,
      behavior across inputs) with a run — synthetic BAM or CLI; otherwise say it is unverified
- [ ] Count-changing step: real-data acceptance on every data type it reaches (DNA tissue
      tumour + matched normal; cfDNA duplex/simplex via `gbcms merge`; WES/WGS without
      realignment; RNA) and more than one truth source — one dataset (e.g. GIAB) is a data
      point, not acceptance; DNA and RNA may need different rules (see add-feature step 6)
- [ ] A decision whose options all cost something is not settled by the least-bad one:
      name what the measurements cannot see (the risk the rule guards, a missing truth
      source), measure it, and look for a rule better on every axis before recommending
- [ ] Analysis aims past the standard: stress each candidate with synthetic shapes for
      what real data cannot show (the guarded risk, the realistic mix, the phantom where
      the allele is absent, natural quality), score what matters (VAF vs truth, phantom
      share), check the oracle trusts no alignment, and size a suspected defect with a
      count-neutral diagnostic on every arm before deciding
- [ ] Implementation carries the analysis: every stress shape becomes a red-first test,
      the rule lands in every path that judges the same question, a rule extended to
      the other allele carries each exclusion and threshold with a mirrored test, each
      new component gets its own unit test, and the measurements are rerun on the
      final build

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
- `SimpleNamespace` only in test fixtures and the one production placeholder, `_zero_counts()`
  (the zero-count row for a variant with no coverage); add no other

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

## Comment & Log Hygiene

Comments and log strings (`debug!`/`warn!`/`info!`, Python `logging`) are read by
someone with **no plan or PR context**. State what the code does, why, and how —
never cite ticket/milestone labels (`CR-1`, `HI-11`, `ME-8`, `P4c`, `M2`, …) in
code or logs. Those are ephemeral and meaningless to a future reader; the ticket
context belongs in the commit message and is recoverable via `git blame`. (This
applies to source only — plan/design docs are exactly where labels belong.)

- [ ] No ticket/milestone labels in comments or log strings
- [ ] A "why" comment states the reason, not the ticket that prompted it
  - ✗ `// ME-8: padding fix` → ✓ `// exclude no-test variants; they inflate the FDR family`
- [ ] Log messages are self-describing — a grep for the message explains the event

## Test Invariants (Assert in Every Counting Test)

```python
assert counts.dp >= counts.rd + counts.ad      # DP includes 'neither'
assert counts.dpf >= counts.rdf + counts.adf   # DPF includes discarded fragments
assert counts.rd == counts.rd_fwd + counts.rd_rev  # strand consistency
assert counts.ad == counts.ad_fwd + counts.ad_rev
```

A new test file starts with a name check: `git ls-files tests/<name>.py` prints nothing,
or create it with the Write tool (it refuses to overwrite a file it has not read). A
shell write (`open(p, "w")`, `>`) once replaced an existing contract's tests unseen.

## Test Suite

Counts drift; get them from the tools (`python -m pytest --collect-only -q | tail -1`,
`cd rust && cargo test`) and never write them down. Find a file with `ls tests`; the
families:

| Family | Files | What it holds |
|:-------|:------|:--------------|
| Contract tests | `tests/test_*_contract.py` | One decided rule each (wrong-length indels, splice evidence, carrier placement, RNA strandedness, ...): the shapes it covers, as red-first tests |
| Read-judgment spec | `tests/test_read_judgment_spec.py`, `tests/read_judgment_cases.py` | Executes the case table of `docs/reference/read-judgment.md`; a decided case fails when the engine breaks its rule |
| Classification oracle | `tests/census.py`, `tests/test_read_census.py` | Each read judged by its own bases; open decisions are strict xfails |
| Binning | `tests/test_binning_invariance.py`, `tests/helpers.py` (`count_checked`, `count_bam_checked`) | Counts do not depend on bin geometry; the Rust `build_genomic_bins` property test is the other half |
| Counting basics | `test_accuracy.py`, `test_shifted_indels.py`, `test_fuzzy_complex.py`, `test_filters.py`, `test_strand_counts.py`, `test_fragment_consensus.py`, `test_multi_allelic.py`, `test_dp_neither.py`, `test_alignment_backend.py` | Per-variant-type accuracy, filters, strands, fragments, siblings, backend concordance |
| Output and I/O | `test_phase2_output.py`, `test_rna_output.py`, `test_maf_reader.py`, `test_maf_preservation.py`, `test_mfsd_flag.py`, `test_observations.py`, `test_merge.py` | Column schemas and gating, readers, companion files, merge |
| CLI, config and pipeline | `test_cli_*.py`, `test_config_isolation.py`, `test_pipeline_*.py`, `test_e2e_*.py` | Argument parsing, validation, end to end |
| Stub parity | `tests/test_rs_stub_parity.py` | `_rs.pyi` matches the `#[pyo3(get)]` fields and signatures |
| Rust | `#[cfg(test)]` modules under `rust/src/` | Unit tests; `cargo test` |
