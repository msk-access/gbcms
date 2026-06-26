#!/usr/bin/env python3
"""PreToolUse / Edit|MultiEdit|Write hook: block edits to .env files.

Secrets belong in a secret manager, not a file that might get committed. Hard
deny on writing any `.env` or `.env.<suffix>`, exempting placeholder/template
env files (.env.example, .env.template, .env.sample) which carry no secrets and
are meant to be committed.

Adapted from BioInfo/claudelicious (MIT). Wire on PreToolUse / Edit|MultiEdit|Write.
"""
import json
import re
import sys

EXEMPT = (".example", ".template", ".sample")

try:
    data = json.load(sys.stdin)
except Exception:
    print("{}")
    sys.exit(0)

path = data.get("tool_input", {}).get("file_path", "") or ""

if re.search(r"(^|/)\.env(\.[^/]+)?$", path) and not path.lower().endswith(EXEMPT):
    print(
        json.dumps(
            {
                "hookSpecificOutput": {
                    "hookEventName": "PreToolUse",
                    "permissionDecision": "deny",
                    "permissionDecisionReason": (
                        "Editing .env files is blocked: keep secrets in a secret "
                        f"manager, not a committable file. Blocked path: {path}"
                    ),
                }
            }
        )
    )
else:
    print("{}")
