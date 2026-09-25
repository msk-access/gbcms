---
name: macos-no-timeout-command
description: "This Mac has no `timeout`/`gtimeout`: `timeout N cmd && ok || fail` always prints fail — check SFTP mounts with a plain read, never through timeout."
metadata:
  node_type: memory
  type: feedback
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
  modified: 2026-09-24T20:03:55.007Z
---

macOS ships no GNU `timeout` (and coreutils' `gtimeout` isn't installed here), so
`timeout 20 head -c 100 <mount file> && echo readable || echo NOT readable`
reports "NOT readable" every time — command not found, not an I/O error. On
2026-09-24 this made two mount-health checks report the DMP/FORTE SFTP mounts as
down when the only real failure was one earlier `OSError: [Errno 5]`.

**Why:** a health check that cannot succeed is a silent failure of the check
itself; it cost a round-trip asking the operator to reconnect.

**How to apply:** check a mount with a plain read (`head -c 60 <file> >/dev/null
&& echo ok`) or `ls`, and treat `Errno 5` from Python/pysam as the real signal.
If a hang is the worry, run the read in the background instead of wrapping it in
`timeout`. See [[user-local-validation-data]] for the mount paths.
