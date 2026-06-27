# ERRORS

The terminal case for `environment` failures: the rule/skill/hook and its trigger
were correct and the failure was **external** (an upstream release, a missing
dependency, stale state, an API outage). These get logged HERE and the
environment gets fixed — they do **not** become a harness edit. Logging an
external failure as a rule learning makes a later pass "fix" something that was
never broken. Stop here for those.

Newest at the top.

---

## [ERR-20260626-001] CI `click` import failure from upstream dependency drift
- **What happened:** `tests/test_cli_dna_rna.py` failed at collection with
  `ModuleNotFoundError: No module named 'click'` across all CI platforms.
- **Rule/hook involved:** none at fault — the test and its trigger were correct.
- **Environment fix:** `click` was imported directly but only present
  transitively via `typer`; an upstream resolution change dropped it in CI's
  lock-bypassing wheel install. Declared `click>=8.0` in `pyproject.toml`.
  (Durable rule captured in memory: [[deps-pyproject-source-of-truth]].)
- **Date:** 2026-06-26
