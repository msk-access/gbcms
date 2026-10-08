---
name: long-runs-visible
description: Long harness runs are started detached (nohup) to outlive the 30-min task limit, so always pair each with a tracked background waiter so the operator sees it running.
metadata:
  type: feedback
---

A detached `nohup` harness run (acceptance, probes, builds) does not appear in the
app's background-task list, so the operator saw "nothing running" while an
acceptance was mid-way (2026-10-02).

**Why:** the operator follows progress from the app; an invisible job reads as a
stalled session.

**How to apply:** start long runs with `nohup ... &` (they outlive the 30-minute
tracked-task limit) and, in the same step, start a tracked background waiter
(`until [ -f run.done ]; do sleep 20; done; <print summary>` with
run_in_background) that ends when the run does. Re-arm the waiter if it times
out. When a run is slow (e.g. full `--trace` on expressed genes), say why and
switch to a faster method rather than letting it grind. Runs that read a network
mount go one at a time, grouped by file, never in parallel (the local mounts
memory has why). Related:
[[harness-cleanup-after-merge]].
