#!/usr/bin/env python3
"""PreToolUse / Bash guard for the model's use of ssh, scp and rsync (the HPC rules).

The operator opens the HPC connection (a Kerberos ticket, then `ssh -fN <host>`); the model
only reuses it, under the rules in .agents/rules/security.md ("HPC"). The code is generic;
the operator's own settings (host, working folders, partitions and their time limit, job-name
prefix, array cap, protected paths) stay on their machine, never in the repo. The hook reads
the first of: $HPC_GUARD_CONFIG, .claude/hpc-guard.local.json in the repo (gitignored), and
~/.claude/hpc-guard.json. Without one it does nothing; a file that is present but unreadable
or incomplete denies every ssh, scp, rsync and sftp. `hpc-guard.example.json` beside this
file shows the format.

Every list here is an allow-list: what is not on it is denied.
- ssh: the configured host only, with options before it from a short list (-a -x -T -n -q -k
  -4 -6 -C, `-O check`, and -o settings that turn forwarding and ticket delegation off, batch
  mode and timeouts). Every ssh needs -a -x; every rsync `-e "ssh -a -x"`, same list.
- a remote command: one line with no variables, substitution, escapes, braces, redirects or
  control characters; commands joined by &&, ; or |, each called by name from the lists
  below. Paths stay in the working folders (writes) or the readable folders (reads), with no
  `..` or `~`. A relative name needs a `cd` into those folders first, joined by `&&`: the
  remote shell starts in the home folder.
- copies: rsync and scp with short option lists; remote paths absolute, in the folders and
  free of spaces and shell characters; local paths in the local roots.
- sbatch: one absolute script, a prefixed job name, a `--time` and the partition that time
  calls for, an array cap, paths in the folders, no `--wrap`. scancel: by a prefixed name.
- a command that names the HPC runs only through ssh, scp or rsync, or a text tool such as
  grep; no shell, interpreter or wrapper runs ssh, scp, rsync or sftp, and command
  substitutions are checked as commands too.

Decisions: allow connection and job status (`-O check`, squeue, sacct, sinfo, hostname,
whoami, pwd); ask the operator for every read, copy, mkdir, sbatch and scancel; deny the
rest. The hook reads the command line only: it is a guard, not a sandbox (see the rules).
FAIL-CLOSED for a command that names the HPC or runs ssh, scp, rsync or sftp (a parse error
or a fault in the hook denies it); fail-open otherwise. Wire on PreToolUse / Bash.
"""

import json
import os
import posixpath
import re
import shlex
import sys

KEYS = {
    "host",
    "hpc_names",
    "remote_roots",
    "read_roots",
    "local_roots",
    "short_partition",
    "long_partition",
    "short_limit_minutes",
    "job_prefix",
    "max_array_tasks",
    "deny_patterns",
}


def _problem(cfg):
    """Why a settings file cannot guard anything, or None. A file that is present is complete."""
    if not isinstance(cfg, dict):
        return "the settings are not a JSON object"
    unknown = sorted(set(cfg) - KEYS)
    if unknown:
        return f"unknown keys {unknown}"
    for key in ("host", "job_prefix", "short_partition", "long_partition"):
        if not isinstance(cfg.get(key), str) or not cfg[key].strip():
            return f"'{key}' must be a non-empty string"
    for key in ("remote_roots", "local_roots", "read_roots", "hpc_names", "deny_patterns"):
        value = cfg.get(key, [])
        if not isinstance(value, list) or not all(isinstance(v, str) and v for v in value):
            return f"'{key}' must be a list of non-empty strings"
    for key in ("remote_roots", "local_roots"):
        if not cfg.get(key):
            return f"'{key}' must name at least one folder"
    roots = cfg["remote_roots"] + cfg.get("read_roots", [])
    roots += [os.path.expanduser(p) for p in cfg["local_roots"]]
    for root in roots:
        if not root.startswith("/") or posixpath.normpath(root) == "/" or ".." in root.split("/"):
            return f"folder '{root}' must be absolute, below / and without '..'"
    for key in ("short_limit_minutes", "max_array_tasks"):
        value = cfg.get(key, 1)
        if type(value) is not int or value < 1:
            return f"'{key}' must be a whole number of 1 or more"
    for pattern in cfg.get("deny_patterns", []):
        try:
            re.compile(pattern)
        except re.error:
            return f"deny pattern '{pattern}' is not a valid regular expression"
    return None


def _load_config():
    """The operator's settings and why they are unusable; ({}, None) means there are none."""
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    named = os.environ.get("HPC_GUARD_CONFIG", "")
    if named and not os.path.isfile(named):
        return {}, f"$HPC_GUARD_CONFIG names no file: {named}"
    places = [
        named,
        os.path.join(here, "hpc-guard.local.json"),
        os.path.expanduser("~/.claude/hpc-guard.json"),
    ]
    for place in places:
        if place and os.path.isfile(place):
            try:
                with open(place) as f:
                    cfg = json.load(f)
            except (OSError, ValueError) as e:
                return {}, f"{place} cannot be read ({type(e).__name__})"
            problem = _problem(cfg)
            return (cfg, None) if problem is None else ({}, f"{place}: {problem}")
    return {}, None


CFG, PROBLEM = _load_config()
HOST = CFG.get("host", "")
HPC_NAMES = tuple({n.lower() for n in CFG.get("hpc_names", []) + [HOST] if n})
REMOTE_ROOTS = tuple(posixpath.normpath(p) for p in CFG.get("remote_roots", []))  # read, write
# Read-only: reads and copies from HPC may also reach these; writes stay in REMOTE_ROOTS.
READ_ROOTS = REMOTE_ROOTS + tuple(posixpath.normpath(p) for p in CFG.get("read_roots", []))
LOCAL_ROOTS = tuple(os.path.normpath(os.path.expanduser(p)) for p in CFG.get("local_roots", []))
SHORT_PARTITION = CFG.get("short_partition", "")
LONG_PARTITION = CFG.get("long_partition", "")
PARTITIONS = {p for p in (SHORT_PARTITION, LONG_PARTITION) if p}
SHORT_LIMIT_MIN = CFG.get("short_limit_minutes", 180)  # longer jobs go to the long partition
JOB_PREFIX = CFG.get("job_prefix", "")
MAX_ARRAY_TASKS = CFG.get("max_array_tasks", 20)
DENY_PATTERNS = [re.compile(p, re.I) for p in CFG.get("deny_patterns", [])]

# Remote commands, called by name. Status needs no approval; the others ask the operator.
STATUS = {"squeue", "sacct", "sinfo", "hostname", "whoami", "pwd", "true"}
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
    "readlink",
    "realpath",
    "zcat",
}
REMOTE_JOINS = {"&&", ";", "|"}
# What the remote shell would expand or read differently: variables and substitution,
# escapes, brace expansion, redirects, and control characters (a newline starts a command).
REMOTE_FORBIDDEN = re.compile(r"[\x00-\x1f\x7f$`\\{}<>]")
# Short options that take a value, per read command, so the value is not taken for a file
# name; PATH_LETTERS are the ones whose value is a file.
VALUE_LETTERS = {
    "grep": "efmABCdD",
    "head": "nc",
    "tail": "ncs",
    "cut": "bcdf",
    "sort": "ktS",
    "uniq": "fsw",
    "stat": "c",
    "du": "dBt",
    "df": "Bt",
    "ls": "Iw",
    "file": "efFmP",
    "mkdir": "m",
}
PATH_LETTERS = {"grep": "f", "file": "fm"}
FIND_VALUE_TESTS = {
    "-name",
    "-iname",
    "-path",
    "-ipath",
    "-wholename",
    "-iwholename",
    "-regex",
    "-iregex",
    "-regextype",
    "-type",
    "-xtype",
    "-maxdepth",
    "-mindepth",
    "-mtime",
    "-mmin",
    "-atime",
    "-amin",
    "-ctime",
    "-cmin",
    "-size",
    "-user",
    "-group",
    "-uid",
    "-gid",
    "-perm",
    "-printf",
    "-links",
    "-inum",
}
FIND_PATH_TESTS = {"-newer", "-anewer", "-cnewer", "-samefile", "-files0-from"}
FIND_ACTIONS = {"-delete", "-exec", "-execdir", "-ok", "-okdir", "-fprint", "-fprint0"}
FIND_ACTIONS |= {"-fprintf", "-fls"}


def _letter_map(pairs):
    """{'A': 'account', ...} from 'A=account ...': short SLURM options and their long names."""
    return dict(p.split("=") for p in pairs.split())


SBATCH_VALUES = _letter_map(
    "A=account a=array b=begin C=constraint c=cpus-per-task D=chdir d=dependency e=error"
    " F=nodefile G=gpus i=input J=job-name L=licenses M=clusters m=distribution N=nodes"
    " n=ntasks o=output p=partition q=qos S=core-spec t=time w=nodelist x=exclude"
)
SBATCH_FLAGS = _letter_map("H=hold k=no-kill O=overcommit Q=quiet s=oversubscribe v=verbose W=wait")
SBATCH_LONG_FLAGS = set(SBATCH_FLAGS.values()) | {
    "parsable",
    "exclusive",
    "requeue",
    "no-requeue",
    "test-only",
    "contiguous",
    "spread-job",
    "use-min-nodes",
}
SBATCH_READ_OPTIONS = {"export", "input"}  # their paths are inputs; every other path is written
SCANCEL_VALUES = _letter_map(
    "A=account M=clusters n=name p=partition q=qos R=reservation s=signal t=state u=user"
    " w=nodelist"
)
SCANCEL_FLAGS = _letter_map("b=batch f=full H=hurry Q=quiet v=verbose")

# ssh: forwarding off, no tty, quiet, no ticket delegation, address family, compression.
SSH_FLAGS = set("axTnqk46C")
SSH_SETTINGS = {
    "forwardagent": "no",
    "forwardx11": "no",
    "forwardx11trusted": "no",
    "gssapidelegatecredentials": "no",
    "clearallforwardings": "yes",
    "batchmode": "yes|no",
    "connecttimeout": r"\d+",
    "serveraliveinterval": r"\d+",
    "serveralivecountmax": r"\d+",
    "loglevel": "quiet|fatal|error|info|verbose",
}
SCP_FLAGS = set("46BCpqrv")
RSYNC_FLAGS = set("acDghilmnopqrstuvxzPSW")
RSYNC_LONG_FLAGS = {
    "archive",
    "recursive",
    "verbose",
    "quiet",
    "compress",
    "human-readable",
    "progress",
    "partial",
    "dry-run",
    "checksum",
    "itemize-changes",
    "stats",
    "update",
    "times",
    "perms",
    "owner",
    "group",
    "links",
    "devices",
    "specials",
    "protect-args",
    "mkpath",
    "size-only",
    "ignore-existing",
    "existing",
    "prune-empty-dirs",
    "whole-file",
    "sparse",
    "one-file-system",
    "no-motd",
}
RSYNC_LONG_VALUES = {"exclude", "include", "files-from", "exclude-from", "include-from"}
RSYNC_LONG_VALUES |= {"log-file", "rsh", "info", "timeout", "contimeout", "bwlimit", "chmod"}
RSYNC_LIST_FILES = {"files-from", "exclude-from", "include-from"}
REMOTE_SPEC = re.compile(r"^([^/:]+):(.*)$")
REMOTE_PATH_BAD = re.compile(r"[\s$`\\;&|<>(){}'\"!#\x00-\x1f\x7f]")

TRANSPORT = {"ssh", "scp", "rsync"}
# Local text tools that may name the HPC: they search or print text and never connect.
TEXT_TOOLS = {"grep", "egrep", "fgrep", "sed", "cat", "echo", "printf", "head", "tail", "wc"}
TEXT_TOOLS |= {"cut", "tr", "sort", "uniq", "diff", "ls", "jq", "column", "true", "false"}
TEXT_TOOLS |= {"test", "[", "cd", "pwd"}
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
    "sed",  # GNU sed's e command
    "sort",  # --compress-program
}
TRANSPORT_CALL = re.compile(r"(?:^|[\s;&|'\"(/=`])(?:ssh|scp|rsync|sftp)(?=[\s'\";&|)`]|$)")
ASSIGNMENT = re.compile(r"[A-Za-z_][A-Za-z0-9_]*=")
OPERATOR_CHARS = set(";&|()\n")


class Decision(Exception):
    def __init__(self, kind, reason):
        super().__init__(reason)
        self.kind, self.reason = kind, reason


def deny(reason):
    raise Decision("deny", reason)


def inside(path, roots):
    return any(path == r or path.startswith(r + "/") for r in roots)


def tokens(cmd):
    """Shell words and operators as (text, is_operator), read as a POSIX shell reads them.

    Quotes and escapes are removed from words; a quoted `;` stays a word. A run of `;&|` is
    one operator, and `(`, `)` and a newline are operators of their own.
    """
    out, word, quoted, i = [], [], False, 0

    def flush():
        nonlocal word, quoted
        if word or quoted:
            out.append(("".join(word), False))
        word, quoted = [], False

    while i < len(cmd):
        c = cmd[i]
        if c in " \t\r":
            flush()
        elif c in OPERATOR_CHARS:
            flush()
            j = i + 1
            if c in ";&|":
                while j < len(cmd) and cmd[j] in ";&|":
                    j += 1
            out.append((cmd[i:j], True))
            i = j - 1
        elif c == "'":
            j = cmd.find("'", i + 1)
            if j < 0:
                raise ValueError("no closing quotation")
            word.append(cmd[i + 1 : j])
            quoted, i = True, j
        elif c == '"':
            j = i + 1
            while j < len(cmd) and cmd[j] != '"':
                if cmd[j] == "\\" and j + 1 < len(cmd) and cmd[j + 1] in '"\\$`\n':
                    j += 1
                    if cmd[j] != "\n":
                        word.append(cmd[j])
                else:
                    word.append(cmd[j])
                j += 1
            if j >= len(cmd):
                raise ValueError("no closing quotation")
            quoted, i = True, j
        elif c == "\\":
            if i + 1 < len(cmd) and cmd[i + 1] != "\n":  # backslash-newline joins the lines
                word.append(cmd[i + 1])
            i += 1
        else:
            word.append(c)
        i += 1
    flush()
    return out


def segments(cmd, joins=False):
    """Split a shell command into simple commands; with joins, as (words, next operator)."""
    out, cur = [], []
    for text, is_op in tokens(cmd):
        if is_op:
            if cur:
                out.append((cur, text))
            cur = []
        else:
            cur.append(text)
    if cur:
        out.append((cur, None))
    return out if joins else [words for words, _ in out]


def substitutions(cmd):
    """The bodies of the `$(...)` and backtick substitutions the shell would run.

    Text in single quotes and escaped characters are literal and skipped; in double quotes a
    substitution still runs. A body's own nested quotes are not followed (it is read roughly,
    and a nested substitution is found when the body is checked in its turn).
    """
    bodies, i, single, double = [], 0, False, False
    while i < len(cmd):
        c = cmd[i]
        if single:
            single = c != "'"
        elif c == "\\":
            i += 1  # the next character is literal
        elif c == "'" and not double:
            single = True
        elif c == '"':
            double = not double
        elif c == "`":
            j = cmd.find("`", i + 1)
            j = len(cmd) if j < 0 else j
            bodies.append(cmd[i + 1 : j])
            i = j
        elif cmd.startswith("$(", i):
            depth, j = 1, i + 2
            while j < len(cmd) and depth:
                depth += {"(": 1, ")": -1}.get(cmd[j], 0)
                j += 1
            bodies.append(cmd[i + 2 : j - 1] if depth == 0 else cmd[i + 2 :])
            i = j - 1
        i += 1
    return bodies


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


def path_like(part):
    return "/" in part or part.startswith(("~", "."))


def remote_path(part, cwd, roots):
    """The absolute remote path a name stands for; it must lie in `roots`."""
    if part.startswith("~"):
        deny("the remote home folder is outside the allowed folders; give an absolute path")
    if ".." in part.split("/"):
        deny(f"no '..' in a remote path: {part}")
    if part.startswith("/"):
        path = posixpath.normpath(part)
    elif cwd is None:
        deny(
            f"'{part}' is relative, and the remote shell starts in the home folder: give an"
            " absolute path, or cd into an allowed folder first (cd <folder> && ...)"
        )
    else:
        path = posixpath.normpath(posixpath.join(cwd, part))
    if not inside(path, roots):
        where = "working" if roots == REMOTE_ROOTS else "readable"
        deny(f"remote path outside the {where} folders: {path}")
    return path


def bundle_value(verb, word):
    """The value-taking letter of a short-option bundle and its attached value, or (None, '')."""
    letters = VALUE_LETTERS.get(verb, "")
    for j, c in enumerate(word[1:], 1):
        if c in letters:
            return c, word[j + 1 :]
        if not c.isalnum():  # the rest may be a path attached to an option not on our list
            if re.search(r"[/~]|\.\.|\.[A-Za-z_]", word[j:]):
                deny(f"give an option's path as a separate word or as --name=path: {word}")
            break
    return None, ""


def check_words(verb, words, cwd, roots):
    """Check the paths among a remote command's words; return its operands: the file names,
    resolved, and '-' for standard input. A bare number is a name too: option values are
    consumed through VALUE_LETTERS, so a number left over is an operand."""
    names, i = [], 0
    pattern_first = verb == "grep" and not has_option(verb, words, "ef", ("regexp", "file"))
    while i < len(words):
        w = words[i]
        if w.startswith("--") and len(w) > 2:
            _, eq, value = w[2:].partition("=")
            for part in re.split(r"[,:]", value) if eq else ():
                if path_like(part):
                    remote_path(part, cwd, roots)
        elif w.startswith("-") and len(w) > 1:
            letter, value = bundle_value(verb, w)
            if verb == "find" and w in FIND_VALUE_TESTS | FIND_PATH_TESTS:
                letter = "f" if w in FIND_PATH_TESTS else "v"
            if letter is not None:
                if not value:
                    i += 1
                    value = words[i] if i < len(words) else ""
                if letter in PATH_LETTERS.get(verb, "") or (verb == "find" and letter == "f"):
                    remote_path(value, cwd, roots)
        elif pattern_first:
            pattern_first = False  # grep's first word is its pattern
        elif w == "-":
            names.append(w)  # standard input: an operand all the same
        elif verb in STATUS:
            if path_like(w):  # status commands take no files; a path stays in the folders
                remote_path(w, cwd, READ_ROOTS)
        else:
            names.append(remote_path(w, cwd, roots))
        i += 1
    return names


def has_option(verb, words, letters, longs):
    """Whether words hold a short option in `letters` or a long one in `longs` (or a prefix
    of it, as GNU tools accept)."""
    for w in words:
        if w.startswith("--") and len(w) > 2:
            name = w[2:].split("=")[0]
            if any(full.startswith(name) for full in longs):
                return True
        elif w.startswith("-") and len(w) > 1:
            for c in w[1:]:
                if c in letters:
                    return True
                if c in VALUE_LETTERS.get(verb, ""):
                    break
    return False


def check_read(verb, words, cwd):
    if verb == "find" and FIND_ACTIONS & set(words):
        deny(f"find may not use {sorted(FIND_ACTIONS & set(words))}: they delete, write or run")
    if verb == "sort" and has_option(
        verb, words, "oT", ("output", "temporary-directory", "compress-program")
    ):
        deny("sort may not write a file or run a program (-o, -T, --compress-program)")
    if verb == "file" and has_option(verb, words, "C", ("compile",)):
        deny("file -C writes a file")
    names = check_words(verb, words, cwd, READ_ROOTS)
    if verb == "uniq" and len(names) > 1:
        deny("uniq writes to its second operand ('-' counts); give it one input at most")
    recursive = verb == "grep" and has_option(
        verb, words, "rR", ("recursive", "dereference-recursive")
    )
    files = [n for n in names if n != "-"]
    if cwd is None and not files and (verb in ("ls", "du", "find") or recursive):
        deny(f"{verb} with no path reads the home folder; give a path in the allowed folders")


def check_cd(words, cwd):
    """The folder a remote cd moves to (None for the home folder)."""
    if any(w.startswith("-") for w in words) or len(words) > 1:
        deny("cd takes one folder and no options")
    return remote_path(words[0], cwd, READ_ROOTS) if words else None


def parse_options(words, values, flags, long_flags, stop):
    """SLURM-style options: ({long name: value}, [arguments]). Deny an option not in the maps
    and a long option given without `=` that is not a flag. With stop, parsing ends at the
    first argument (sbatch: the rest belong to the script)."""
    opts, args, i = {}, [], 0
    while i < len(words):
        w = words[i]
        if args and stop:
            args.append(w)
        elif w.startswith("--") and len(w) > 2:
            name, eq, value = w[2:].partition("=")
            if not eq and name not in long_flags:
                deny(f"give --{name}=value as one word")
            opts[name] = value
        elif w.startswith("-") and len(w) > 1:
            for j, c in enumerate(w[1:], 1):
                if c in flags:
                    opts[flags[c]] = ""
                    continue
                if c not in values:
                    deny(f"option -{c} is not on the HPC list for this command")
                value = w[j + 1 :]
                if not value:
                    i += 1
                    if i >= len(words):
                        deny(f"option -{c} needs a value")
                    value = words[i]
                opts[values[c]] = value
                break
        else:
            args.append(w)
        i += 1
    return opts, args


def check_sbatch(words):
    opts, args = parse_options(
        words, SBATCH_VALUES, SBATCH_FLAGS, SBATCH_LONG_FLAGS | {"wrap"}, stop=True
    )
    if any(len(name) >= 2 and "wrap".startswith(name) for name in opts):
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
    for name, value in opts.items():
        roots = READ_ROOTS if name in SBATCH_READ_OPTIONS else REMOTE_ROOTS
        for part in re.split(r"[,:=]", value):
            if path_like(part):
                remote_path(part, None, roots)
    if not args or not args[0].startswith("/") or not args[0].endswith((".sbatch", ".sh")):
        deny("sbatch needs one absolute script path in the working folders")
    remote_path(args[0], None, REMOTE_ROOTS)
    for word in args[1:]:  # the script's own arguments: inputs it may read
        for part in re.split(r"[,:=]", word):
            if path_like(part):
                remote_path(part, None, READ_ROOTS)


def check_scancel(words):
    opts, args = parse_options(
        words, SCANCEL_VALUES, SCANCEL_FLAGS, set(SCANCEL_FLAGS.values()), stop=False
    )
    if args or not opts.get("name", "").startswith(JOB_PREFIX):
        deny(f"scancel only by --name={JOB_PREFIX}<name>")


def check_remote(command):
    """Return 'allow' or 'ask' for a remote command; raise for a denial."""
    if not command.strip():
        deny("no interactive shell; give a command")
    if any(pat.search(command) for pat in DENY_PATTERNS):
        deny("this remote command names a protected path (cBioPortal)")
    stripped = strip_benign(command)
    bad = REMOTE_FORBIDDEN.search(stripped)
    if bad:
        deny(
            f"no {bad.group()!r} in a remote command: variables, substitution, escapes, braces,"
            " redirects and control characters (a newline starts a second command) change"
            " what the remote shell runs"
        )
    level, cwd, before = "allow", None, None
    for words, after in segments(stripped, joins=True):
        if after is not None and after not in REMOTE_JOINS:
            deny(f"join remote commands with &&, ; or | only, not {after!r}")
        if before == ";":
            cwd = None  # it runs even when a command before it failed, a cd among them
        verb, args = words[0], words[1:]
        if "/" in verb:
            deny(f"call a remote command by its name, not by a path: {verb}")
        if verb == "cd":
            target = check_cd(args, cwd)
            if before != "|" and after == "&&":
                cwd = target  # what follows runs only once this cd has succeeded
        elif verb in STATUS:
            check_words(verb, args, cwd, READ_ROOTS)
        elif verb in READ:
            level = "ask"
            check_read(verb, args, cwd)
        elif verb == "mkdir":
            level = "ask"
            check_words(verb, args, cwd, REMOTE_ROOTS)
        elif verb == "sbatch":
            level = "ask"
            check_sbatch(args)
        elif verb == "scancel":
            level = "ask"
            check_scancel(args)
        else:
            deny(f"remote command '{verb}' is not on the HPC list")
        before = after
    return level


def ssh_setting(value):
    """Check one `-o` setting against the short list; return its key."""
    m = re.fullmatch(r"\s*([A-Za-z0-9]+)\s*(?:=\s*|\s+)(\S+)\s*", value)
    key = m.group(1).lower() if m else ""
    if key not in SSH_SETTINGS or not re.fullmatch(SSH_SETTINGS[key], m.group(2), re.I):
        deny(f"ssh -o {value} is not on the HPC list")
    return key


def ssh_options(words, i):
    """Read ssh options from words[i]: (index after them, flags set, the -O command or None)."""
    flags, control = set(), None
    while i < len(words) and words[i].startswith("-") and len(words[i]) > 1:
        w = words[i]
        for j, c in enumerate(w[1:], 1):
            if c in SSH_FLAGS:
                flags.add(c)
                continue
            if c not in "oO":
                deny(
                    f"ssh option -{c} is not on the HPC list (forwarding, ticket delegation,"
                    " jumps, another config, user or port)"
                )
            value = w[j + 1 :]
            if not value:
                i += 1
                if i >= len(words):
                    deny(f"ssh -{c} needs a value")
                value = words[i]
            if c == "O":
                control = value
            else:
                key = ssh_setting(value)
                flags |= {"forwardagent": {"a"}, "forwardx11": {"x"}}.get(key, set())
            break
        i += 1
    return i, flags, control


def check_ssh(words):
    i, flags, control = ssh_options(words, 1)
    host = words[i] if i < len(words) else None
    if host != HOST:
        deny(f"ssh only to {HOST}, not {host}")
    rest = words[i + 1 :]
    if control is not None:
        if control == "check" and not rest:
            return "allow"
        deny(f"ssh -O {control} is the operator's; only -O check is allowed")
    if not {"a", "x"} <= flags:
        deny("every ssh to HPC needs -a -x (no agent or X11 forwarding)")
    if rest and rest[0].startswith("-"):
        deny("ssh options go before the host")
    return check_remote(" ".join(rest))


def check_local(path):
    full = os.path.normpath(os.path.expanduser(path))
    if not full.startswith("/") or not inside(full, LOCAL_ROOTS):
        deny(f"a local copy path must be absolute and in the local folders: {path}")


def check_remote_spec(spec, roots):
    host, path = spec.groups()
    if host != HOST:
        deny(f"copies only to or from {HOST}, not {host}")
    if REMOTE_PATH_BAD.search(path):
        deny(f"no spaces or shell characters in a remote path: {path}")
    if not path.startswith("/"):
        deny(f"a remote copy path must be absolute: {path or '(home)'}")
    remote_path(path, None, roots)


def check_rsh(rsh):
    words = shlex.split(rsh) if rsh else []
    if not words or os.path.basename(words[0]) != "ssh":
        deny('rsync to HPC needs -e "ssh -a -x"')
    end, flags, control = ssh_options(words, 1)
    if end != len(words) or control is not None:
        deny("the rsync -e value is ssh and its options only (no host, command or -O)")
    if not {"a", "x"} <= flags:
        deny('rsync to HPC needs -e "ssh -a -x"')


def check_copy(words):
    """Check an rsync or scp call; None when it names no remote host (a local copy)."""
    tool = os.path.basename(words[0])
    if any(pat.search(" ".join(words)) for pat in DENY_PATTERNS):
        deny("this copy names a protected path (cBioPortal)")
    problems, settings, lists, logs, positionals = [], [], [], [], []
    flags, values = (RSYNC_FLAGS, "e") if tool == "rsync" else (SCP_FLAGS, "lo")
    rsh, i = None, 1
    while i < len(words):
        w = words[i]
        if tool == "rsync" and w.startswith("--") and len(w) > 2:
            name, eq, value = w[2:].partition("=")
            if name in RSYNC_LONG_VALUES:
                if not eq:
                    i += 1
                    value = words[i] if i < len(words) else ""
                if name == "rsh":
                    rsh = value
                elif name in RSYNC_LIST_FILES:
                    lists.append(value)
                elif name == "log-file":
                    logs.append(value)
            elif eq or name not in RSYNC_LONG_FLAGS:
                problems.append(f"rsync option --{name} is not on the HPC list")
        elif w.startswith("-") and len(w) > 1:
            for j, c in enumerate(w[1:], 1):
                if c in flags:
                    continue
                if c not in values:
                    problems.append(f"{tool} option -{c} is not on the HPC list")
                    break
                value = w[j + 1 :]
                if not value:
                    i += 1
                    value = words[i] if i < len(words) else ""
                if c == "e":
                    rsh = value
                elif c == "o":
                    settings.append(value)
                elif not value.isdigit():
                    problems.append("scp -l takes a number (Kbit/s)")
                break
        else:
            positionals.append(w)
        i += 1
    if not any(REMOTE_SPEC.match(w) for w in positionals):
        return None
    if problems:
        deny(problems[0])
    if len(positionals) < 2:
        deny("a copy names its sources and one destination")
    for value in settings:
        ssh_setting(value)
    for n, w in enumerate(positionals):
        is_dest = n == len(positionals) - 1
        spec = REMOTE_SPEC.match(w)
        if spec:
            check_remote_spec(spec, REMOTE_ROOTS if is_dest else READ_ROOTS)
        else:
            check_local(w)
    for value in lists:
        spec = REMOTE_SPEC.match(value)
        if spec:
            check_remote_spec(spec, READ_ROOTS)
        else:
            check_local(value)
    for value in logs:
        check_local(value)
    if tool == "rsync":
        check_rsh(rsh)
    return "ask"


def names_hpc(text):
    low = text.lower()
    return any(n in low for n in HPC_NAMES)


def evaluate(cmd, depth=0):
    hpc = names_hpc(cmd)
    try:
        segs = segments(strip_benign(cmd) if hpc else cmd)
    except ValueError:
        if hpc or TRANSPORT_CALL.search(cmd):
            deny("could not parse this command, and it names the HPC or runs ssh, scp or rsync")
        return None
    hpc = hpc or any(names_hpc(w) for words in segs for w in words)  # names split by quotes
    level, saw = None, False
    for words in segs:
        k = 0
        while k < len(words) and ASSIGNMENT.match(words[k]):
            k += 1  # leading VAR=value assignments
        if k == len(words):
            continue
        tool = os.path.basename(words[k])
        if tool == "sftp":
            deny("use rsync, not sftp")
        if tool in EXECUTORS and TRANSPORT_CALL.search(" ".join(words[k + 1 :])):
            deny("ssh, scp and rsync run directly, never through a shell, interpreter or wrapper")
        if hpc and tool not in TRANSPORT | TEXT_TOOLS:
            deny(
                "a command that names the HPC runs only through ssh, scp or rsync, or a text"
                " tool such as grep; not through a shell, interpreter, wrapper or other transport"
            )
        if tool not in TRANSPORT:
            continue
        if k:
            deny("no variable assignments before ssh, scp or rsync")
        lv = check_ssh(words[k:]) if tool == "ssh" else check_copy(words[k:])
        if lv is None:
            continue
        saw = True
        level = "ask" if "ask" in (lv, level) else "allow"
    for body in substitutions(cmd) if depth < 4 else ():
        lv = evaluate(body, depth + 1)  # a substitution runs as a command of its own
        if lv is not None:
            saw = True
            level = "ask" if "ask" in (lv, level) else "allow"
    if saw and hpc:
        rest = strip_benign(cmd)
        if ">" in rest or "<" in rest:
            deny("no redirect in an HPC command; copy with rsync or scp so each copy is approved")
    return level


def main():
    try:
        cmd = json.load(sys.stdin).get("tool_input", {}).get("command", "") or ""
    except Exception:
        sys.exit(0)
    if not HOST and not PROBLEM:
        sys.exit(0)  # no HPC configured on this machine: nothing to guard
    reason = "this reads, copies or submits on HPC; the operator approves each one"
    try:
        if PROBLEM and TRANSPORT_CALL.search(cmd):
            deny(f"the HPC guard's settings are unusable ({PROBLEM}); fix them first")
        level = None if PROBLEM else evaluate(cmd)
    except Decision as d:
        level, reason = d.kind, d.reason
    except Exception as e:  # a fault in the guard fails closed for what may reach the HPC
        level = None
        if names_hpc(cmd) or TRANSPORT_CALL.search(cmd):
            level, reason = "deny", f"the guard could not check this command ({type(e).__name__})"
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
