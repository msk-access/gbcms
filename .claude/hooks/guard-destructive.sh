#!/usr/bin/env bash
# PreToolUse(Bash) guardrail for the gbcms harness.
# Reads the tool-call JSON on stdin. Exit 2 blocks the call; exit 0 allows it.
# Deliberately narrow: parses the .command field and only blocks when the
# dangerous verb is the actual command word (so an `echo 'rm -rf /'` mention,
# or a path named elsewhere, does NOT trip it).
set -euo pipefail
input="$(cat)"

# Extract just the command string; fail-open (allow) if it can't be parsed.
cmd="$(printf '%s' "$input" | python3 -c 'import sys,json
try:
    print(json.load(sys.stdin).get("tool_input", {}).get("command", ""))
except Exception:
    pass' 2>/dev/null || true)"
[ -z "$cmd" ] && exit 0

# Recursive force-delete of a top-level/sensitive path, only when `rm` is the
# command word (start of a segment or after ; && || |).
if printf '%s' "$cmd" \
  | grep -Eq '(^|[;&|][[:space:]]*)rm[[:space:]]+-[a-zA-Z]*[rf][a-zA-Z]*[[:space:]]+([^[:space:]]+[[:space:]]+)*(/|~|/\*|\$HOME)([[:space:]]|$)'; then
  echo "BLOCKED: recursive delete of a top-level/sensitive path. Narrow the target or use the scratchpad." >&2
  exit 2
fi

# Direct push to main is prohibited (Git Flow: release/* -> main via PR).
if printf '%s' "$cmd" \
  | grep -Eq '(^|[;&|][[:space:]]*)git[[:space:]]+push([[:space:]]+[^;&|]*)?[[:space:]]main([[:space:]]|$)'; then
  echo "BLOCKED: direct push to main is prohibited. Use Git Flow (PR release/* -> main)." >&2
  exit 2
fi

exit 0
