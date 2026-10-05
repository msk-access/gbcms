---
name: pipe-masks-gate-exit
description: Never chain `pytest ... | tail && git commit` — the pipe returns tail's exit code, so a failing test commits.
metadata:
  type: feedback
---

`cmd | tail -1 && next` runs `next` even when `cmd` failed: the pipeline's status is
`tail`'s. On 2026-10-05 (group 6) a fixture assertion failed and the commit went in
anyway (a0a90b4a), fixed in 679632e5.

**Why:** a gate that cannot fail is not a gate; the red-first discipline depends on it.

**How to apply:** redirect to a log and test `$?` (`pytest ... > log 2>&1; rc=$?;
tail -1 log; [ $rc -eq 0 ] && git commit ...`), or `set -o pipefail` in the shell
line. Same for `cargo test | grep`. Related: [[commit-before-review-workflows]].
