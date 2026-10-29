#!/usr/bin/env python3
"""PreToolUse / Bash guard for the model's use of ssh, scp and rsync (the HPC rules).

The operator opens the HPC connection (a Kerberos ticket, then `ssh -fN <host>`); the model
only reuses it, under the rules in .agents/rules/security.md ("HPC"). The code is generic;
the operator's own settings (host, working folders, partitions and their time limit, job-name
prefix, array cap, protected paths) stay on their machine, never in the repo. The hook reads
the first of: $HPC_GUARD_CONFIG, .claude/hpc-guard.local.json in the repo (gitignored), and
~/.claude/hpc-guard.json. Without one it does nothing; `hpc-guard.example.json` beside this
file shows the format.

- deny: any other host; agent or X11 forwarding left on (every ssh needs -a -x, every rsync
  `-e "ssh -a -x"`); port forwarding, jump hosts, proxy commands, another config file or a
  control command other than `-O check`; an interactive shell; command substitution or a
  redirect; a remote command not on the lists below; a path outside the working folders
  (remote) or the local roots; deleting, `find -delete/-exec`; protected paths (cBioPortal);
  `sbatch` without one absolute script, a prefixed job name, a `--time` and the partition
  that time calls for, or with `--wrap` or an array over the task cap; `scancel` other than
  by a prefixed job name; rsync `--delete`; sftp.
- ask the operator: every copy (rsync, scp), every read of a file or folder, and every
  `sbatch`, `scancel` and `mkdir`.
- allow: connection and job status only (`-O check`, squeue, sacct, sinfo, hostname, whoami).

FAIL-CLOSED for a command that names the HPC (a parse error denies it); fail-open otherwise.
Wire on PreToolUse / Bash.
"""

import json
import os
import re
import shlex
import sys


def _load_config():
    """The operator's settings: the first readable of the three places, or {} (guard off)."""
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    places = [
        os.environ.get("HPC_GUARD_CONFIG", ""),
        os.path.join(here, "hpc-guard.local.json"),
        os.path.expanduser("~/.claude/hpc-guard.json"),
    ]
    for place in places:
        if place and os.path.isfile(place):
            try:
                with open(place) as f:
                    return json.load(f) or {}
            except (OSError, ValueError):
                return {"host": "\0unreadable"}  # a broken config denies every HPC command
    return {}


CFG = _load_config()
HOST = CFG.get("host", "")
HPC_NAMES = tuple(CFG.get("hpc_names", [HOST] if HOST else []))
REMOTE_ROOTS = tuple(CFG.get("remote_roots", []))  # read and write
# Read-only: reads and copies from HPC may also reach these; writes stay in REMOTE_ROOTS.
READ_ROOTS = REMOTE_ROOTS + tuple(CFG.get("read_roots", []))
LOCAL_ROOTS = tuple(os.path.normpath(os.path.expanduser(p)) for p in CFG.get("local_roots", []))
SHORT_PARTITION = CFG.get("short_partition", "")
LONG_PARTITION = CFG.get("long_partition", "")
PARTITIONS = {p for p in (SHORT_PARTITION, LONG_PARTITION) if p}
SHORT_LIMIT_MIN = int(CFG.get("short_limit_minutes", 180))  # longer jobs go to the long partition
JOB_PREFIX = CFG.get("job_prefix", "")
MAX_ARRAY_TASKS = int(CFG.get("max_array_tasks", 20))
DENY_PATTERNS = [re.compile(p, re.I) for p in CFG.get("deny_patterns", [])]

STATUS = {"squeue", "sacct", "sinfo", "hostname", "whoami", "date", "id", "true", "cd", "pwd"}
READ = {
    "ls",
    "cat",
    "head",
    "tail",
    "wc",
    "du",
    "df",
    "find",
    "stat",
    "sha256sum",
    "md5sum",
    "grep",
    "sort",
    "uniq",
    "cut",
    "file",
    "test",
    "readlink",
    "realpath",
    "zcat",
}
CHANGE = {"sbatch", "scancel", "mkdir"}
SEPARATORS = {";", "&&", "||", "|", "&", "(", ")", "\n"}
SSH_ARG_OPTS = set("bcDEeFIiJLlmOoPpQRSWw")  # short options that take a value
SSH_DENY_FLAGS = set("AXYLRDJWwFM")
SSH_DENY_O = re.compile(
    r"^(ProxyJump|ProxyCommand|LocalForward|RemoteForward|DynamicForward|"
    r"ForwardAgent\s*=\s*yes|ForwardX11\s*=\s*yes|ForwardX11Trusted\s*=\s*yes|"
    r"PermitLocalCommand|LocalCommand|ControlPath|ControlMaster)",
    re.I,
)
RSYNC_VALUE_OPTS = {
    "-e",
    "--rsh",
    "--files-from",
    "--exclude-from",
    "--include-from",
    "--exclude",
    "--include",
    "--filter",
    "-f",
    "--log-file",
    "--chmod",
}
RSYNC_PATH_VALUES = {"--files-from", "--exclude-from", "--include-from", "--log-file"}


class Decision(Exception):
    def __init__(self, kind, reason):
        super().__init__(reason)
        self.kind, self.reason = kind, reason


def deny(reason):
    raise Decision("deny", reason)


def under(path, roots):
    p = os.path.normpath(os.path.expanduser(path))
    return any(p == r or p.startswith(r + "/") for r in roots)


def segments(cmd):
    """Split a shell command into simple commands, respecting quotes."""
    lex = shlex.shlex(cmd, posix=True, punctuation_chars=";&|()")
    lex.whitespace_split = True
    lex.commenters = ""
    out, cur = [], []
    for tok in lex:
        if tok in SEPARATORS or set(tok) <= set(";&|()"):
            if cur:
                out.append(cur)
            cur = []
        else:
            cur.append(tok)
    if cur:
        out.append(cur)
    return out


def options(words, short):
    """Parse --long=value, --long value and -s value / -svalue options up to the first argument.

    Returns ({long name: value}, [arguments]). `short` maps a short letter to its long name;
    a value-less long option (a flag) gets "".
    """
    opts, args, i = {}, [], 0
    while i < len(words):
        w = words[i]
        if args:
            args.append(w)
        elif w.startswith("--"):
            name, eq, val = w[2:].partition("=")
            if (
                not eq
                and i + 1 < len(words)
                and not words[i + 1].startswith("-")
                and name in short.values()
            ):
                i += 1
                val = words[i]
            opts[name] = val
        elif w.startswith("-") and len(w) > 1:
            name = short.get(w[1], w[1])
            val = w[2:]
            if not val and i + 1 < len(words):
                i += 1
                val = words[i]
            opts[name] = val
        else:
            args.append(w)
        i += 1
    return opts, args


def slurm_minutes(value):
    """A SLURM --time value in minutes (MM, MM:SS, HH:MM:SS, D-HH, D-HH:MM, D-HH:MM:SS), or None."""
    m = re.fullmatch(r"(?:(\d+)-)?(\d+)(?::(\d+))?(?::(\d+))?", value.strip())
    if not m:
        return None
    days, a, b, c = (int(x) if x else None for x in m.groups())
    if days is not None:  # D-HH[:MM[:SS]]
        return days * 1440 + a * 60 + (b or 0) + (1 if c else 0)
    if b is None:  # MM
        return a
    if c is None:  # MM:SS
        return a + (1 if b else 0)
    return a * 60 + b + (1 if c else 0)  # HH:MM:SS


def strip_benign(cmd):
    """Drop the redirects that write nowhere: 2>&1 and the ones to /dev/null."""
    return re.sub(r"\s*(2>&1|2>\s*/dev/null|>\s*/dev/null)", " ", cmd)


def check_remote_paths(words, roots=None):
    roots = roots or REMOTE_ROOTS
    for w in words:
        for part in re.split(r"[=,:]", w):
            if part.startswith(("/", "~")) and part != "/dev/null" and not under(part, roots):
                deny(f"remote path outside the allowed folders: {part}")


def check_remote(command):
    """Return 'allow' or 'ask' for a remote command; raise for a denial."""
    if not command.strip():
        deny("no interactive shell; give a command")
    if any(pat.search(command) for pat in DENY_PATTERNS):
        deny("this remote command names a protected path (cBioPortal)")
    if "$(" in command or "`" in command:
        deny("no command substitution in a remote command")
    stripped = strip_benign(command)
    if ">" in stripped or "<" in stripped:
        deny("no file redirect in a remote command; jobs write through sbatch --output")
    level = "allow"
    for words in segments(stripped):
        verb = os.path.basename(words[0])
        if verb in STATUS:
            pass
        elif verb in READ or verb in CHANGE:
            level = "ask"
        else:
            deny(f"remote command '{verb}' is not on the HPC list")
        check_remote_paths(
            words[1:], READ_ROOTS if verb in READ or verb in STATUS else REMOTE_ROOTS
        )
        joined = " ".join(words)
        if verb == "find" and re.search(r"\s-(delete|exec|execdir|ok|okdir|fprint)", " " + joined):
            deny("find may not delete or run commands")
        if verb == "sbatch":
            opts, args = options(
                words[1:], {"J": "job-name", "p": "partition", "a": "array", "t": "time"}
            )
            if "wrap" in opts:
                deny("sbatch --wrap is not allowed; submit a script")
            if not opts.get("job-name", "").startswith(JOB_PREFIX):
                deny(f"sbatch needs --job-name={JOB_PREFIX}<name>")
            if opts.get("partition") not in PARTITIONS:
                deny(f"sbatch needs --partition={SHORT_PARTITION} or {LONG_PARTITION}")
            minutes = slurm_minutes(opts.get("time", ""))
            if minutes is None:
                deny("sbatch needs --time (D-HH:MM:SS, HH:MM:SS or minutes)")
            if minutes <= SHORT_LIMIT_MIN and opts["partition"] != SHORT_PARTITION:
                deny(f"a job of {SHORT_LIMIT_MIN} minutes or less runs on {SHORT_PARTITION}")
            if minutes > SHORT_LIMIT_MIN and opts["partition"] != LONG_PARTITION:
                deny(f"a job over {SHORT_LIMIT_MIN} minutes runs on {LONG_PARTITION}")
            if "array" in opts:
                cap = re.search(r"%(\d+)$", opts["array"])
                if not cap or int(cap.group(1)) > MAX_ARRAY_TASKS:
                    deny(f"an array needs a %N limit of {MAX_ARRAY_TASKS} or less")
            if (
                len(args) < 1
                or not args[0].startswith("/")
                or not args[0].endswith((".sbatch", ".sh"))
            ):
                deny("sbatch needs one absolute script path in the working folders")
        if verb == "scancel":
            opts, args = options(words[1:], {"n": "name"})
            if not opts.get("name", "").startswith(JOB_PREFIX) or args:
                deny(f"scancel only by --name={JOB_PREFIX}<name>")
    return level


def check_ssh(words):
    i, host, flags, control = 1, None, set(), None
    while i < len(words) and host is None:
        w = words[i]
        if w.startswith("-") and len(w) > 1:
            opt, j = w[1:], 0
            while j < len(opt):
                c = opt[j]
                if c in SSH_DENY_FLAGS:
                    deny(f"ssh option -{c} is not allowed (forwarding, jumps or another config)")
                flags.add(c)
                if c in SSH_ARG_OPTS:
                    val = opt[j + 1 :]
                    if not val:
                        i += 1
                        val = words[i] if i < len(words) else ""
                    if c == "o":
                        if SSH_DENY_O.search(val):
                            deny(f"ssh -o {val} is not allowed")
                        if re.match(r"ForwardAgent\s*=\s*no", val, re.I):
                            flags.add("a")
                        if re.match(r"ForwardX11\s*=\s*no", val, re.I):
                            flags.add("x")
                    if c == "O":
                        control = val
                    break
                j += 1
        else:
            host = w
        i += 1
    if host != HOST:
        deny(f"ssh only to {HOST}, not {host}")
    if control is not None:
        if control == "check" and i >= len(words):
            return "allow"
        deny(f"ssh -O {control} is the operator's; only -O check is allowed")
    if not {"a", "x"} <= flags:
        deny("every ssh to HPC needs -a -x (no agent or X11 forwarding)")
    return check_remote(" ".join(words[i:]))


def check_local(path):
    if not path.startswith(("/", "~")):
        deny(f"local copy paths must be absolute and under ~/test/gbcms/harness: {path}")
    if not under(path, LOCAL_ROOTS):
        deny(f"local path outside ~/test/gbcms/harness: {path}")


def check_copy(words):
    tool = words[0]
    joined = " ".join(words)
    if any(pat.search(joined) for pat in DENY_PATTERNS):
        deny("this copy names a protected path (cBioPortal)")
    if tool == "rsync" and re.search(r"--delete|--remove-source-files", joined):
        deny("rsync may not delete files")
    rsh, remote_seen, i, positionals = None, False, 1, []
    while i < len(words):
        w = words[i]
        if w.startswith("-"):
            name, _, val = w.partition("=")
            if tool == "rsync" and name in RSYNC_VALUE_OPTS and not val:
                i += 1
                val = words[i] if i < len(words) else ""
            if tool == "rsync" and name in ("-e", "--rsh"):
                rsh = val
            elif tool == "rsync" and w.startswith("-e") and len(w) > 2 and "=" not in w:
                rsh = w[2:]
            if name in RSYNC_PATH_VALUES and val:
                check_local(val)
            if tool == "scp" and w in ("-F", "-J", "-o", "-S", "-i", "-P", "-l", "-c"):
                if w in ("-F", "-J", "-S"):
                    deny(f"scp option {w} is not allowed")
                i += 1
            i += 1
            continue
        positionals.append(w)
        i += 1
    for n, w in enumerate(positionals):
        is_dest = n == len(positionals) - 1
        m = re.match(r"^([^/:]+):(.*)$", w)
        if m:
            remote_seen = True
            if m.group(1) != HOST:
                deny(f"copies only to or from {HOST}, not {m.group(1)}")
            roots = REMOTE_ROOTS if is_dest else READ_ROOTS
            if not under(m.group(2), roots):
                where = "the working folders" if is_dest else "the readable folders"
                deny(
                    f"remote {'destination' if is_dest else 'source'} outside {where}: {m.group(2) or '(home)'}"
                )
        else:
            check_local(w)
    if not remote_seen:
        return None
    if tool == "rsync":
        r = shlex.split(rsh) if rsh else []
        if not r or os.path.basename(r[0]) != "ssh" or "-a" not in r or "-x" not in r:
            deny('rsync to HPC needs -e "ssh -a -x"')
    return "ask"


#: Commands that run another command from their arguments: a hidden ssh call could hide there.
EXECUTORS = {
    "bash",
    "sh",
    "zsh",
    "dash",
    "ksh",
    "fish",
    "eval",
    "exec",
    "command",
    "builtin",
    "xargs",
    "env",
    "nohup",
    "timeout",
    "gtimeout",
    "watch",
    "parallel",
    "nice",
    "time",
    "sudo",
    "su",
    "script",
    "expect",
    "screen",
    "tmux",
    "osascript",
    "python",
    "python3",
    "perl",
    "ruby",
    "node",
    "uv",
    "uvx",
}


def evaluate(cmd):
    names_hpc = any(n in cmd for n in HPC_NAMES)
    try:
        segs = segments(strip_benign(cmd) if names_hpc else cmd)
    except ValueError:
        if names_hpc:
            deny("could not parse this HPC command")
        return None
    level, saw = None, False
    for words in segs:
        k = 0
        while k < len(words) and re.match(r"^[A-Za-z_][A-Za-z0-9_]*=", words[k]):
            k += 1  # leading VAR=value assignments
        words = words[k:]
        if not words:
            continue
        tool = os.path.basename(words[0])
        if tool == "sftp":
            deny("use rsync, not sftp")
        if names_hpc and (tool in EXECUTORS or re.match(r"python3?(\.\d+)?$", tool)):
            deny("a command that names the HPC may not run through a shell, interpreter or wrapper")
        if tool == "ssh":
            lv = check_ssh(words)
        elif tool in ("scp", "rsync"):
            lv = check_copy(words)
            if lv is None:
                continue
        else:
            continue
        saw = True
        level = "ask" if lv == "ask" or level == "ask" else "allow"
    if saw and names_hpc:
        rest = strip_benign(cmd)
        if ">" in rest or "<" in rest:
            deny("no redirect in an HPC command; copy with rsync or scp so each copy is approved")
    return level


def main():
    broken = HOST == "\0unreadable"
    if not broken and (not HOST or not JOB_PREFIX or not REMOTE_ROOTS):
        sys.exit(0)  # no HPC configured on this machine: nothing to guard
    try:
        cmd = json.load(sys.stdin).get("tool_input", {}).get("command", "") or ""
    except Exception:
        sys.exit(0)
    try:
        if broken and re.search(r"\b(ssh|scp|rsync|sftp)\b", cmd):
            deny(
                "the HPC guard's settings file cannot be read; fix it before any ssh, scp or rsync"
            )
        level = evaluate(cmd) if not broken else None
    except Decision as d:
        level, reason = d.kind, d.reason
    else:
        reason = "this reads, copies or submits on HPC; the operator approves each one"
    if level in ("deny", "ask"):
        print(
            json.dumps(
                {
                    "hookSpecificOutput": {
                        "hookEventName": "PreToolUse",
                        "permissionDecision": level,
                        "permissionDecisionReason": "HPC rules: " + reason,
                    }
                }
            )
        )
        sys.exit(2 if level == "deny" else 0)
    sys.exit(0)


if __name__ == "__main__":
    main()
