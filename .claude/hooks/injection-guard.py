#!/usr/bin/env python3
"""PreToolUse / Bash injection-and-exfil guard (ask, never hard-deny).

Defense-in-depth for an agent with broad Bash + web access. Targets only the
genuinely anomalous shapes an injected instruction (web page, repo file, MCP
response) would use to exfiltrate secrets or run remote code. Everything else
passes. Returns 'ask' (human confirm), never a deny, so rare-but-legit cases
(a known `curl ... | sh` installer) stay possible behind one confirm.

FAIL-OPEN on any parse/internal error — this is one layer, not the only one
(block-rm-rf.py is the hard-deny sibling). Adapted from BioInfo/claudelicious
(MIT); host-specific gateway heuristics removed. Wire on PreToolUse / Bash.
"""
import json
import re
import sys


def ask(reason):
    print(
        json.dumps(
            {
                "hookSpecificOutput": {
                    "hookEventName": "PreToolUse",
                    "permissionDecision": "ask",
                    "permissionDecisionReason": reason,
                }
            }
        )
    )
    sys.exit(0)


try:
    data = json.load(sys.stdin)
    cmd = data.get("tool_input", {}).get("command", "") or ""
except Exception:
    sys.exit(0)  # fail-open

low = cmd.lower()

# 1a. Remote content piped straight into a shell (classic installer / injection).
if re.search(r"\b(curl|wget|fetch)\b[^|]*\|\s*(sudo\s+)?(sh|bash|zsh|fish)\b", low):
    ask(
        "Injection guard: remote content is piped into a shell (curl/wget | sh). "
        "Confirm this is a trusted installer and not an injected instruction."
    )

# 1b. Remote content piped into a BARE interpreter that executes it as code.
#     A fixed program (-c/-e/-m flag or a script path) treats piped bytes as
#     DATA (json parse, etc.) and is allowed; bare `-` is caught.
if re.search(
    r"\b(curl|wget|fetch)\b[^|]*\|\s*(sudo\s+)?"
    r"(python3?|node|ruby|perl)\b(?!\s+(-[A-Za-z]|[^\s|;&-]))",
    low,
):
    ask(
        "Injection guard: remote content is piped into a bare interpreter that "
        "would execute it as code. Confirm this is intentional."
    )

# 2. ANTHROPIC_BASE_URL set inline to an EXTERNAL literal URL (key-exfil vector).
#    Local/loopback is allowed.
m = re.search(r"\banthropic_base_url\s*=\s*['\"]?(https?://[^\s'\"]+)", low)
if m:
    host = re.sub(r"^https?://", "", m.group(1)).split("/")[0].split(":")[0]
    if host not in ("localhost", "127.0.0.1", "0.0.0.0"):
        ask(
            "Injection guard: ANTHROPIC_BASE_URL is set inline to an external "
            "URL (API-key exfiltration vector). Confirm this is intentional."
        )

# 3. A high-value SECRET source AND an OUTBOUND network sink in the SAME command
#    (the exfiltration shape). Reading a secret alone is fine; piping it out is
#    not. ssh/scp/rsync are excluded (authenticated, point-to-point, daily ops).
secret = re.compile(
    r"(~/\.ssh/|/\.ssh/id_|id_(rsa|ed25519|ecdsa)\b|\.password-store|"
    r"(^|\s)pass\s+show\b|\.aws/credentials|\.gcp|\.pem\b|\.p12\b|private[_-]?key)",
    re.I,
)
sink = re.compile(
    r"\b(curl|wget|nc|ncat|netcat|telnet|requests\.(post|put|get)|urllib|http\.client)\b",
    re.I,
)
if secret.search(cmd) and sink.search(cmd):
    ask(
        "Injection guard: a private key / credential source and an outbound "
        "network command appear together (exfiltration shape). Confirm intent."
    )

sys.exit(0)
