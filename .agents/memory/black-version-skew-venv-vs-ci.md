---
name: black-version-skew-venv-vs-ci
description: Lint tools are stale locally vs CI (black AND clippy) — run CI's versions or you get false "clean" / phantom "drift".
metadata:
  type: feedback
---

Local lint tools here lag CI, so a locally-clean run can still fail CI (and vice-versa).
Two confirmed cases:

**black** — the repo `.venv` ships a **stale compiled `black` 25.9.0** (`.venv/bin/black`)
with **no pip** (uv-style venv), so `pip install` silently targets mambaforge. CI installs
**black 26.5.1** (`pyproject` dev pin `black>=25.9.0` → newest). They disagree on the
multiline `plotly_calls.append(f"""...""")` in `mfsd_report.py`: 25.9.0 wants a rewrap,
26.5.1 (CI) considers the original clean → trusting the venv black once produced a phantom
"drift" and a near-miss "fix" that would have *broken* CI.

**clippy** — local `rustup` stable lagged (rustc 1.91 / clippy 0.1.91) while CI's
`dtolnay/rust-toolchain@stable` was 1.96. The newer clippy flags lints the old one misses
(e.g. `clippy::unnecessary_sort_by` at `rna.rs:402`), so a locally-clean `cargo clippy`
failed CI once `-D warnings` was enabled.

**Why:** CI is the source of truth; a stale local toolchain gives false confidence.

**How to apply:** before trusting a lint gate, run the **CI-matching** version:
- black → `/Users/shahr2/mambaforge/bin/black --check src/ tests/` (26.5.x), not `.venv/bin/black`.
- clippy/cargo test → `rustup update stable` first, then `rustup run stable cargo clippy
  --all-targets -- -D warnings` (both feature configs) and `cargo test`.
`ruff`/`mypy` have not shown skew. See [[harness-stay-small]]. CI now runs `cargo test` +
clippy itself (both feature configs), so regressions surface there too.
