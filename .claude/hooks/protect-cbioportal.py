#!/usr/bin/env python3
"""PreToolUse hook: keep protected study data (cBioPortal files) out of tool output.

Whatever a tool prints is sent to the model provider as part of the conversation, so
protected files are never read, listed, searched or edited by a tool. They are used only
through local scripts:
- the fixed summary script (`cbio_summary.py`: schema and counts, small counts hidden)
  may print to the console;
- any other command that names protected data must send all of its output, stdout and
  stderr, to a file under a protected directory, with no pipes, chains or substitutions.

Protected: every absolute path listed in ~/.claude/protected-paths.txt (kept outside the
repository), and cBioPortal file names (data_*.txt, meta_*.txt, case_lists) wherever they
sit. Wire on PreToolUse for Bash and for the file tools.
"""
import json
import os
import re
import sys

HOME = os.path.expanduser("~")
LIST = os.path.join(HOME, ".claude", "protected-paths.txt")
SUMMARY = os.path.join(HOME, "test", "gbcms", "harness", "cbio", "cbio_summary.py")
NAMES = re.compile(r"(?<![\w.-])(data_[A-Za-z0-9_]+\.txt|meta_[A-Za-z0-9_]+\.txt|case_lists)(?![\w-])")


def roots():
    try:
        with open(LIST) as fh:
            return [ln.strip().rstrip("/") for ln in fh if ln.strip() and not ln.startswith("#")]
    except OSError:
        return []


def expand(s):
    s = s.replace("${HOME}", HOME).replace("$HOME", HOME)
    return re.sub(r"(^|[\s=:'\"])~(?=/)", lambda m: m.group(1) + HOME, s)


def protected_in(text, rs):
    text = expand(text)
    return any(r in text for r in rs) or bool(NAMES.search(text.replace(SUMMARY, "")))


def deny(reason):
    print(json.dumps({"hookSpecificOutput": {
        "hookEventName": "PreToolUse",
        "permissionDecision": "deny",
        "permissionDecisionReason": reason,
    }}))
    sys.exit(0)


def bash_allowed(cmd, rs):
    c = expand(cmd).strip()
    # The summary script, run on its own: its output is schema and counts only.
    summary = re.escape(SUMMARY)
    if re.fullmatch(rf"(\S*/)?python3?(\s+-I)?\s+{summary}(\s+[^|;&<>`$()\n]+)?", c):
        return True
    # Any other command: one plain command, every output stream into a protected file.
    if re.search(r"[|;`\n]|&&|\$\(|<\(|\btee\b", c):
        return False
    m = re.search(r"(?:&>|>>?)\s*(\S+)", c)
    if not m or not any(m.group(1).startswith(r) for r in rs):
        return False
    if "&>" in c:
        return True
    err = re.search(r"2>>?\s*(\S+)", c)
    return "2>&1" in c[m.end():] or bool(err and any(err.group(1).startswith(r) for r in rs))


try:
    data = json.load(sys.stdin)
except Exception:
    print("{}")
    sys.exit(0)

tool = data.get("tool_name", "")
inp = data.get("tool_input", {}) or {}
rs = roots()
RULE = (
    "Protected study data (cBioPortal) stays on this machine: a tool may not read, list, "
    "search or edit it, because tool output is sent to the model provider. Run "
    f"`python3 -I {SUMMARY} --all` for schema and counts, or write a script and run it with "
    "all output (stdout and stderr) redirected to a file under a protected directory "
    "(no pipes or chains). For GitHub text that names these files, use --body-file."
)

if tool == "Bash":
    cmd = inp.get("command", "") or ""
    if protected_in(cmd, rs) and not bash_allowed(cmd, rs):
        deny(RULE)
else:
    fields = " ".join(str(inp.get(k, "")) for k in ("file_path", "path", "pattern", "glob", "notebook_path"))
    if protected_in(fields, rs):
        deny(RULE)

print("{}")
