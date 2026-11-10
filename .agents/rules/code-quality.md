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
- [ ] Count-changing step: follow `.agents/rules/decision-protocol.md` (real-data acceptance on
      every data type and more than one truth source, the best rule rather than the least bad,
      red-first tests for every stress shape, the rule landed in every path that judges the
      same question, the measurements rerun on the final build)

**AFTER implementing, compare code against implementation plan:**
- [ ] Re-read the relevant plan section(s)
- [ ] Diff what was implemented vs. what the plan specified
- [ ] If there are **gaps**: classify as intentional or unintentional, notify user
- [ ] If the implementation **diverged**: update plan to match reality

## Lint Suite

All of it must pass before every commit: ruff → black (CI's version) → mypy → clippy → cargo
test → pytest, in that fix order. The commands, the exit-code rule and the tool-version notes
are in the qa-check skill (`.agents/skills/qa-check/SKILL.md`).

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

## Dependencies

`pyproject.toml` is the single source of truth: runtime dependencies in `[project.dependencies]`
(floors measured by bisection; CI's floors leg installs them with `uv --resolution lowest-direct`
and `scripts/check_floors.py` checks it), contributors' tools in the PEP 735
`[dependency-groups]` (`test`, `docs`, `dev`).
- If code, tests included, does `import X`, declare `X` where it is imported (runtime, or the
  `test` group when only tests import it, as `click` is); never rely on it arriving transitively.
  A transitive `click` once vanished from a fresh CI resolve and broke collection on every platform.
- There is no committed `uv.lock` (it is gitignored): CI, the Dockerfile and Nextflow build the
  wheel and resolve fresh, so a lock only drifts. Reintroduce one only with `uv sync --frozen` in
  CI and Docker and a `uv lock --check` drift gate: all or nothing.
- Reproducibility belongs in the Docker image: `docker/requirements.lock` (uv pip compile, hashes,
  linux/amd64, Python 3.11) is installed with `--require-hashes --no-deps`, then `pip check` gates
  drift. It is refreshed at release (release guide, step 4). The org's containers image does not
  use it (`.agents/memory/references.md`).

## Comment & Log Hygiene

Comments and log strings (`debug!`/`warn!`/`info!`, Python `logging`) are read by
someone with **no plan or PR context**. State what the code does, why, and how —
never cite ticket/milestone labels (`CR-1`, `HI-11`, `ME-8`, `P4c`, `M2`, …) in
code or logs (operator, 2026-06-27). Those are ephemeral and meaningless to a future reader; the ticket
context belongs in the commit message and is recoverable via `git blame`. (This
applies to source only — plan/design docs are exactly where labels belong.)

- [ ] No ticket/milestone labels in comments or log strings
- [ ] A "why" comment states the reason, not the ticket that prompted it
  - ✗ `// ME-8: padding fix` → ✓ `// exclude no-test variants; they inflate the FDR family`
- [ ] Log messages are self-describing — a grep for the message explains the event

## Test Invariants (Assert in Every Counting Test)

The four invariants (DP, DPF, strand consistency) are in AGENTS.md ("Counting test invariants").

A new test file starts with a name check: `git ls-files tests/<name>.py` prints nothing,
or create it with the Write tool (it refuses to overwrite a file it has not read). A
shell write (`open(p, "w")`, `>`) once replaced an existing contract's tests unseen.

## Test Suite

Keep tests minimal and high-signal: each is a maintenance contract, and a fix that removes duplication
is better than added code. Counts drift; get them from the tools (`python -m pytest --collect-only -q | tail -1`,
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
