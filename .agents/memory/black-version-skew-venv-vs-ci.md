---
name: black-version-skew-venv-vs-ci
description: Lint-gate with mambaforge black 26.5.1 (matches CI), not the venv's stale 25.9.0 — they disagree and cause false "drift".
metadata:
  type: feedback
---

The repo `.venv` ships a **stale compiled `black` 25.9.0** (`.venv/bin/black`) and has
**no pip** (uv-style minimal venv), so `pip install` silently targets mambaforge instead.
CI installs **black 26.5.1** (`pyproject` dev pin `black>=25.9.0` → newest). 25.9.0 and
26.5.1 format differently — notably the multiline `plotly_calls.append(f"""...""")`
calls in `src/gbcms/report/mfsd_report.py`: 25.9.0 wants them rewrapped, 26.5.1 (and CI)
consider the original clean.

**Why:** trusting the venv's `black` once produced a phantom "drift" in `mfsd_report.py`
and a near-miss "fix" that would have *broken* CI lint (CI's 26.5.1 wants the original).

**How to apply:** for the black step of the lint gate, run the CI-matching binary
explicitly — `/Users/shahr2/mambaforge/bin/black --check src/ tests/` (26.5.1) — not the
activated-venv `black`. Confirm a file is *actually* drifted under 26.5.x before
reformatting it. `ruff`/`mypy`/the Rust gate are unaffected. See [[harness-stay-small]].
