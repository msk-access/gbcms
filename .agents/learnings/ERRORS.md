# ERRORS

The terminal case for `environment` failures: the rule/skill/hook and its trigger
were correct and the failure was **external** (an upstream release, a missing
dependency, stale state, an API outage). These get logged HERE and the
environment gets fixed — they do **not** become a harness edit. Logging an
external failure as a rule learning makes a later pass "fix" something that was
never broken. Stop here for those.

Newest at the top.

---

## [ERR-20261116-001] reads of the two reference genomes hang from the session
- **What happened:** on 2026-11-16 the profile agent's runs of develop 473e8920 stopped in
  `prepare_variants`, opening `genome.fa.fai`: no error, no CPU, until the process was killed. Every
  `open()` of `~/Downloads/genome.fa` and of the OneDrive copy of `Homo_sapiens_assembly19.fasta`
  (`~/Library/CloudStorage/OneDrive-.../test_data/`) hung, `head -c` included; `stat` worked. The gate
  (`batch680.py`, `FA38`) and the profile arms read the same files; on 2026-11-15 they read them.
- **Rule/hook involved:** none at fault. The likely cause is outside the harness: a macOS privacy check for
  Downloads waiting for an answer on the Mac's screen, or a OneDrive file kept online only.
- **Environment fix:** the operator's choice: copy both files, with their `.fai`, to
  `~/test/gbcms/refs/` (outside Downloads and OneDrive) and point the manifests there, or allow the access
  on the Mac. Not worked round from the session (the Terminal panel runs outside the sandbox).
- **Cleared:** at 10:14 the same day both files read again (a 4 KB read of each `.fai`, a 64 KB read 2 GB
  into each FASTA, each under a 15 to 20 second limit), with nothing changed from the session. The cause
  is not known; the likely one is an answer to the Mac's privacy prompt. If it comes back, ask the
  operator first, and test with a time limit; never list the folder.
- **Date:** 2026-11-16

## [ERR-20261115-002] CI's Rocky 8 test job could not pull its container image, a second time
- **What happened:** on PR #440, `Rocky 8, Python 3.12, manylinux_2_28` (`wheels.yml`, step "Install and
  test the wheel") failed after 7 s: `docker run` pulled `quay.io/rockylinux/rockylinux:8`, and the
  download from quay's file server broke (`httpReadSeeker: failed open ... EOF`, exit 125). No test ran;
  the other Rocky 8 jobs of the run passed. The same failure as ERR-20261114-001, one day later.
- **Rule/hook involved:** none at fault. Red Hat's status page showed Quay operational, with no incident
  from 2026-11-14 to 2026-11-16: brief download errors, not an outage.
- **Environment fix:** the failed job was rerun (`gh run rerun <run> --failed`). Proposed, not done: since
  it repeats, pull the image in its own step with up to 3 tries before `docker run` (a workflow change,
  for the operator to approve).
- **Date:** 2026-11-15

## [ERR-20261115-001] a sign-in prompt ended the session and stopped its tracked background runs
- **What happened:** on 2026-11-15 the desktop app stopped the session for a sign-in (not for credits).
  Every tracked background task ended with it: the local #335 gate chain stopped after one arm of the
  s1c gate (905 default-path runs, all ok; no run left half-written, as each finished run has its DONE
  mark). The HPC jobs ran on, as they do not depend on the session. The operator: "this is not
  acceptable to abrutly stop for sign in, if i am out of credits it makes sense".
- **Rule/hook involved:** none at fault; the long-run rule (`docs/harness/measurement-harness.md`,
  "Make a long run visible") already says to start a run that must outlive the session with `nohup`,
  paired with a tracked waiter on its process.
- **Environment fix:** start a local run longer than about an hour detached (`nohup`, a new session)
  with a tracked waiter on its pid, so a stop of the app does not stop the run; the gate skips the
  runs with a DONE mark, so a restart loses nothing. The sign-in stop itself is the app's, outside the
  harness.
- **Date:** 2026-11-15

## [ERR-20261114-001] CI's Rocky 8 test job could not pull its container image
- **What happened:** on PR #438, `Rocky 8, Python 3.10, manylinux_2_28` (`wheels.yml`, step "Install and
  test the wheel") failed after 10 s: `docker run` pulled `quay.io/rockylinux/rockylinux:8`, and the
  pull broke ("short read: expected 76611732 bytes but got 35651585: unexpected EOF", exit 125). No
  test ran. The same image passed in the Python 3.12 and 3.13 jobs of the same run; the PR's code was
  not at fault.
- **Rule/hook involved:** none at fault.
- **Environment fix:** rerun the failed job once the run has ended (`gh run rerun <run> --failed`; a
  job of a run in progress "cannot be rerun"). If it repeats, check quay.io's status before changing
  any workflow.
- **Date:** 2026-11-14

## [ERR-20261103-001] CI's macOS build got no hosted arm64 runner
- **What happened:** on PR #407, `Build (macOS arm64)` (`wheels.yml`, `runs-on: macos-latest`) ended as
  cancelled after 15 minutes with no runner and no steps: "The job was not acquired by Runner of type
  hosted even after multiple attempts", with the notice "Due to capacity constraints, jobs targeting
  macOS arm64 runners may experience longer queue times". `gh pr checks` lists it as a failure. The
  macOS test job in the Tests workflow passed on the same commit; the PR's code was not at fault.
- **Rule/hook involved:** none at fault.
- **Environment fix:** rerun the failed jobs (`gh run rerun <run> --failed`). If it repeats, read the job's
  annotations (`gh api repos/<owner>/<repo>/check-runs/<job>/annotations`) before changing any workflow.
- **Date:** 2026-11-03

## [ERR-20261031-001] CI's macOS build cancelled by a GitHub runner-image brownout
- **What happened:** on PR #396, `Build (macOS arm64)` (`wheels.yml`, `runs-on: macos-14`) was cancelled
  with no steps run: "This job failed because of a scheduled macOS 14 runner image brownout. macOS 14
  runner images will be retired on November 2, 2026" (actions/runner-images#13518). Every other check
  passed; the PR's code was not at fault.
- **Rule/hook involved:** none at fault. The release and wheel workflows and their contract test
  (`tests/test_macos_wheel_contract.py`) pinned `macos-14`, which the provider retires.
- **Environment fix:** the macOS jobs run on `macos-latest` (arm64; operator, 2026-10-31); the test pins
  that label; the wheel tag stays `macosx_11_0_arm64` (set by `MACOSX_DEPLOYMENT_TARGET`, not the runner).
- **Date:** 2026-10-31

## [ERR-20261007-001] Local `mkdocs build --strict` aborted on "has no git logs"
- **What happened:** at the 6.6.0 cut, the strict docs build stopped on one
  warning, `[git-revision-date-localized-plugin] '.../docs' has no git logs`.
  The plugin runs `git log --diff-filter=r`, and Apple Git 2.50.1 (Apple
  Git-155) returns nothing for that filter on any path; `git log -- docs` works.
- **Rule/hook involved:** none at fault — the docs are clean, and CI's git
  (`deploy-docs.yml`, ubuntu) is unaffected.
- **Environment fix:** for a local check, point GitPython at a shim that drops
  `--diff-filter=r` (`GIT_PYTHON_GIT_EXECUTABLE=<shim> mkdocs build --strict`);
  the build then passes. Recheck after a git update.
- **Date:** 2026-10-07

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
