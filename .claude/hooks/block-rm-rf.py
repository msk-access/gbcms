#!/usr/bin/env python3
"""PreToolUse / Bash hard-deny for `rm -rf`.

Matches `rm` invoked AS A COMMAND (at a command boundary) with recursive+force
flags. Anchoring on command position is the point: a naive substring check
hard-denies any command that merely *contains* "rm -rf" (grep, echo, a heredoc,
or a learnings note discussing it) — none of which delete anything. This still
blocks real invocations (rm -rf, rm -fr, rm -Rf, sudo rm -rf, `&& rm -rf`,
rm -r -f, rm --recursive --force, find | xargs rm -rf) while letting mentions pass.

FAIL-OPEN on parse error (never wedge Bash). Adapted from BioInfo/claudelicious
(MIT). Wire on PreToolUse / Bash.
"""
import json
import re
import sys

try:
    data = json.load(sys.stdin)
    command = data.get("tool_input", {}).get("command", "") or ""
except Exception:
    sys.exit(0)  # fail-open

RM_RF = re.compile(
    r"(?:^|[\n;&|(])\s*"
    r"(?:sudo\s+)?(?:xargs\s+(?:-\S+\s+)*)?(?:sudo\s+)?"
    r"\brm\s+"
    r"(?:"
    r"-[a-z]*r[a-z]*f[a-z]*"
    r"|-[a-z]*f[a-z]*r[a-z]*"
    r"|-r[a-z]*\s+-[a-z]*f"
    r"|-f[a-z]*\s+-[a-z]*r"
    r"|--recursive\b[^\n;&|]*?--force\b"
    r"|--force\b[^\n;&|]*?--recursive\b"
    r")",
    re.I,
)

if RM_RF.search(command):
    print(
        json.dumps(
            {
                "hookSpecificOutput": {
                    "hookEventName": "PreToolUse",
                    "permissionDecision": "deny",
                    "permissionDecisionReason": (
                        "rm -rf is blocked. Narrow the target, delete a specific "
                        "path, or use the session scratchpad for throwaway files."
                    ),
                }
            }
        )
    )
    sys.exit(2)

sys.exit(0)
