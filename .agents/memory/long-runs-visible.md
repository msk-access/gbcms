---
name: long-runs-visible
description: Launch long harness runs as the tracked background task itself so the operator sees them; nohup only for runs that must outlive the session, then with a waiter on the runner's PID (never a log marker).
metadata:
  type: feedback
---

A detached `nohup` harness run (acceptance, probes, builds) does not appear in the
app's background-task list, so the operator saw "nothing running" while an
acceptance was mid-way (2026-10-02).

**Why:** the operator follows progress from the app; an invisible job reads as a
stalled session.

**How to apply:** start a long run as the tracked background task itself (Bash with
run_in_background; `| tee <log>` keeps a log), so the operator sees the run, not a
stand-in, and is told when it ends: tracked tasks ran 6+ hours on 2026-10-19, so the
old "30-minute limit" that led to `nohup` launches does not hold. Use `nohup ... &`
only for a run that must outlive the session, and then pair it with a tracked waiter
on the runner's **process**: `while kill -0 <pid> 2>/dev/null; do sleep 60; done;
<print the log's summary>`. Never wait on a marker line or file the run writes at its
end: a run that dies never writes it, and a restart reusing the log name overwrites
it, so the waiter hangs (twice on #294, 2026-10-19/20). One runner per log name; a run
and its waiter stop together; clear a stopped run's partial outputs before a restart
(the arm runners skip any output directory that exists). When a run is slow (e.g. full `--trace` on expressed genes), say why and
switch to a faster method rather than letting it grind. Runs that read a network
mount go one at a time, grouped by file, never in parallel (the local mounts
memory has why). Related:
[[harness-cleanup-after-merge]].
