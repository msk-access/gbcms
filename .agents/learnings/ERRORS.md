# ERRORS

The terminal case for `environment` failures: the rule/skill/hook and its trigger
were correct and the failure was **external** (an upstream release, a missing
dependency, stale state, an API outage). These get logged HERE and the
environment gets fixed — they do **not** become a harness edit. Logging an
external failure as a rule learning makes a later pass "fix" something that was
never broken. Stop here for those.

Newest at the top.

---

## [ERR-20261005-001] CI Codecov upload crashed on a TLS handshake failure
- **What happened:** on PR #217, `test (ubuntu-latest, 3.11)` failed at "Upload
  coverage to Codecov" (`write EPROTO ... ssl/tls alert handshake failure`, SSL
  alert 40) after "Run tests with coverage" passed. develop's previous run
  uploaded fine: an external, transient Codecov fault.
- **Rule/hook involved:** none at fault — the tests passed.
- **Environment fix:** the step already set `fail_ci_if_error: false`, but the
  action crashed before applying it; the step is now `continue-on-error: true`
  in `.github/workflows/test.yml`, so a Codecov outage cannot fail the job.
- **Date:** 2026-10-05

## [ERR-20260626-001] CI `click` import failure from upstream dependency drift
- **What happened:** `tests/test_cli_dna_rna.py` failed at collection with
  `ModuleNotFoundError: No module named 'click'` across all CI platforms.
- **Rule/hook involved:** none at fault — the test and its trigger were correct.
- **Environment fix:** `click` was imported directly but only present
  transitively via `typer`; an upstream resolution change dropped it in CI's
  lock-bypassing wheel install. Declared `click>=8.0` in `pyproject.toml`.
  (Durable rule captured in memory: [[deps-pyproject-source-of-truth]].)
- **Date:** 2026-06-26
