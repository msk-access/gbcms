#!/usr/bin/env python3
"""PreToolUse / Bash deny for direct pushes to main (gbcms Git Flow).

Git Flow here is feature/* -> develop -> main, and main is only updated via a
release PR. This denies a direct `git push ... main` at the model level (the
.githooks/pre-push hook enforces the same for human/other-tool pushes). Anchored
at a command boundary so a branch merely *named* in a string doesn't trip it.

FAIL-OPEN on parse error. Wire on PreToolUse / Bash.
"""
import json
import re
import sys

try:
    data = json.load(sys.stdin)
    cmd = data.get("tool_input", {}).get("command", "") or ""
except Exception:
    sys.exit(0)  # fail-open

# Isolate the `git push ...` segment (command-word anchored, up to the next
# separator), then require `main` as a standalone ref token WITHIN that segment.
# Scoping to the segment avoids matching "main" elsewhere in a compound command
# (e.g. a commit message); the token guard avoids "main-fix" / "origin/maint".
seg = re.search(r"(?:^|[\n;&|(])\s*git\s+push\b[^\n;&|]*", cmd)
target_main = seg and re.search(r"(?<![\w-])main(?![\w-])", seg.group(0))

if target_main:
    print(
        json.dumps(
            {
                "hookSpecificOutput": {
                    "hookEventName": "PreToolUse",
                    "permissionDecision": "deny",
                    "permissionDecisionReason": (
                        "Direct push to main is blocked (Git Flow). Push to a "
                        "feature/* branch and open a PR to develop; main is "
                        "updated only via a release PR."
                    ),
                }
            }
        )
    )
    sys.exit(2)

sys.exit(0)
