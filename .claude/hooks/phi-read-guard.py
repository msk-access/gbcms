#!/usr/bin/env python3
"""PreToolUse hook: keep patient data (PHI) out of tool output.

Whatever a tool prints is sent to the model provider, so a tool may not read, list, search or
print patient data, nor the file names and paths in a patient folder (a BAM stem is a sample
ID). The rule is in .agents/rules/security.md ("Patient BAM slices"); this hook enforces it.

What is PHI comes from ~/.claude/phi-paths.txt (kept outside the repository, owned by the
operator; $PHI_PATHS_FILE names another file, for tests). One entry per line, `#` starts a
comment. An entry is an absolute root (every path under it is PHI) or a pattern with `*` for
one folder level or inside a folder name, such as `<root>/*/panel` (a patient arm folder under
any run label) or `<root>/out_*`. `sanitizer <absolute path>` names one more accepted sanitizer
script (the default is ~/test/gbcms/harness/hpc/kit/phi_sanitize.py). `phi-paths.example.txt`
beside this file shows the format. Without the file the hook does nothing; a file that is
present but unusable denies every Bash, Read, Grep, Glob and NotebookRead call until fixed.

A Bash command that names a PHI path (or runs in a PHI folder after `cd`) is allowed only in
these forms. Everything else is denied, with a reason that says what to do instead.
1. File operations that print no content: rsync, cp, mv, ln, mkdir, trash, without an option
   that prints names (-v, -P, -i, --progress, --list-only, ...). Patient data may go only into
   another PHI folder (rsync: or to another host); a copy to any other place would hide it from
   this list. For rsync, `--stats` plus `2>&1 | tail -N` (N up to 25) or
   `2>&1 | grep -E 'Number of files|Total transferred'` shows its summary lines.
2. A script run whose whole output (stdout and stderr) passes through the sanitizer, as the
   last stage and with nothing after it:
   `python3 -I <script>.py ... 2>&1 | python3 -I <sanitizer> [args]`.
3. One command whose stdout and stderr both go to a file under a PHI root
   (`cmd > <root>/out.txt 2>&1`), with no pipe, and no substitution.
Read, Grep, Glob and NotebookRead on a PHI path are denied. Write and Edit are allowed, also for a
script under a PHI root: they show only what the model wrote, and the tool refuses to overwrite
or edit a file that was not read first, which Read cannot do here.

A folder that holds a PHI folder (a run folder, a harness root) is guarded against the tools that
search below it: find, du, tree, rg, ag, fd, grep -r, ls -R, diff -r, and the Grep and Glob tools
(a Glob pattern says where it looks). A depth limit that stops at or above the PHI folder
(find -maxdepth, du -s or -d, tree -L) lists no content and is allowed.

Paths are read as the shell reads them: `~`, $HOME, quotes, `..`, braces, `*` patterns, letter
case, relative names (with the hook input's `cwd`, and `cd` inside the command) and symlinks.

A guard, not a sandbox: it reads the command line only. A script that builds a PHI path inside
itself, a variable the hook cannot see through, or a command built to hide its paths can still
pass; so can copying a folder that merely holds PHI folders into another such folder. Fail-closed
for a command that names a PHI path (a parse error or a fault in the hook denies it);
fail-open otherwise. Wire on PreToolUse for Bash and for Read|Grep|Glob|NotebookRead.
"""

import collections
import fnmatch
import json
import os
import posixpath
import re
import shlex
import sys

HOME = os.path.expanduser("~")
LIST_ENV = "PHI_PATHS_FILE"
LIST = os.path.join(HOME, ".claude", "phi-paths.txt")
DEFAULT_SANITIZER = os.path.join(HOME, "test", "gbcms", "harness", "hpc", "kit", "phi_sanitize.py")
GLOB = re.compile(r"[*?\[]")
SUB = "\x01"  # stands for a command substitution in the text the tokenizer reads
SYSTEM_BIN = {"/bin", "/usr/bin", "/usr/local/bin", "/opt/homebrew/bin", "/usr/sbin", "/sbin"}


class Deny(Exception):
    pass


def deny(reason):
    raise Deny(reason)


def real(path):
    try:
        return os.path.realpath(path)
    except (OSError, ValueError):
        return path


def expand_home(text):
    text = text.replace("${HOME}", HOME).replace("$HOME", HOME)
    return re.sub(r"""(^|[\s=:,'"(])~(?=/|$|[\s'"):;,])""", lambda m: m.group(1) + HOME, text)


def tilde(word):
    return re.sub(r"(^|[=:,])~(?=/|$)", lambda m: m.group(1) + HOME, word)


# ---------------------------------------------------------------- the PHI list


def _seg_regex(seg):
    if seg == "*":
        return r"""[^/\s'"`;|&<>()]+"""
    out = []
    for ch in seg:
        if ch == "*":
            out.append(r"""[^/\s'"`;|&<>()]*""")
        elif ch == "?":
            out.append(r"[^/\s]")
        else:
            out.append(re.escape(ch))
    return "".join(out)


def parse_list(text):
    """(entries as lower-case segment lists, sanitizer paths). Raises ValueError on a bad line."""
    entries, sanitizers = [], [DEFAULT_SANITIZER]
    for n, raw in enumerate(text.splitlines(), 1):
        line = re.sub(r"(^|\s)#.*$", "", raw).strip()
        if not line:
            continue
        m = re.match(r"sanitizer\s*[:=]?\s*(\S.*)$", line, re.I)
        if m:
            path = expand_home(m.group(1).strip())
            if not path.startswith("/"):
                raise ValueError(f"line {n}: the sanitizer must be an absolute path")
            sanitizers.append(posixpath.normpath(path))
            continue
        path = expand_home(line)
        if not path.startswith("/"):
            raise ValueError(f"line {n}: an entry must be an absolute path")
        segs = [s for s in posixpath.normpath(path).split("/") if s]
        if not segs:
            raise ValueError(f"line {n}: '/' cannot be an entry")
        if "**" in segs:
            raise ValueError(f"line {n}: use * for one folder level; ** is not supported")
        if GLOB.search(segs[0]):
            raise ValueError(f"line {n}: the first folder cannot be a pattern")
        entries.append(segs)
        k = next((i for i, s in enumerate(segs) if GLOB.search(s)), len(segs))
        prefix = "/" + "/".join(segs[:k])
        resolved = real(prefix)
        if resolved != prefix:  # the root itself is reached through a symlink
            entries.append([s for s in resolved.split("/") if s] + segs[k:])
    return [[s.lower() for s in e] for e in entries], sanitizers


def load_list():
    """(entries, sanitizers, text regex or None, problem or None)."""
    named = os.environ.get(LIST_ENV, "")
    path = named or LIST
    if not named and not os.path.lexists(path):
        return [], [DEFAULT_SANITIZER], None, None  # no list on this machine: nothing to guard
    try:
        with open(path) as fh:
            entries, sanitizers = parse_list(fh.read())
    except ValueError as e:
        return [], [DEFAULT_SANITIZER], None, f"{path}: {e}"
    except OSError as e:
        return [], [DEFAULT_SANITIZER], None, f"{path} cannot be read ({type(e).__name__})"
    if not entries:
        return [], sanitizers, None, None
    alts = ["/" + "/+".join(_seg_regex(s) for s in e) for e in entries]
    text_re = re.compile("(?:" + "|".join(alts) + r")(?![\w.@+\-])", re.I)
    return entries, sanitizers, text_re, None


ENTRIES, SANITIZERS, TEXT_RE, PROBLEM = load_list()


# ------------------------------------------------------------ path matching


def seg_match(p, e):
    """Whether a path folder name p (which may hold shell patterns) can be the entry's folder e."""
    if e == "*":
        return True
    if GLOB.search(p):
        return True if GLOB.search(e) else fnmatch.fnmatchcase(e, p)
    if GLOB.search(e):
        return fnmatch.fnmatchcase(p, e)
    return p == e


def rel(segs, entry):
    """'under' if the path is the entry or inside it, 'above' if it is a folder that holds the
    entry (a search below it reaches PHI), else None. A `**` in the path may reach anything below."""
    for i in range(min(len(segs), len(entry))):
        if segs[i] == "**":
            return "under"
        if not seg_match(segs[i], entry[i]):
            return None
    return "under" if len(segs) >= len(entry) else "above"


def rank(a, b):
    order = {None: 0, "above": 1, "under": 2}
    return a if order[a] >= order[b] else b


def segments_level(path):
    segs = [s for s in path.lower().split("/") if s]
    level = None
    for entry in ENTRIES:
        level = rank(level, rel(segs, entry))
        if level == "under":
            break
    return level


def path_level(joined):
    """Level of an absolute path as typed (lexically) and as the disk resolves it (symlinks)."""
    level = segments_level(posixpath.normpath(joined))
    if level != "under":
        level = rank(level, segments_level(real(joined)))
    return level


def brace_expand(s, limit=64):
    m = re.search(r"\{([^{}]*,[^{}]*)\}", s)
    if not m or limit <= 1:
        return [s]
    out = []
    for alt in m.group(1).split(","):
        out += brace_expand(s[: m.start()] + alt + s[m.end() :], limit)
        if len(out) >= limit:
            break
    return out[:limit]


def candidates(text):
    """The strings in one shell word that may be paths: the word, an option's value (--x=PATH,
    -xPATH), the parts of a list (a:b, a,b) and the brace expansions of each."""
    if SUB in text or len(text) > 1024:
        return []
    t = tilde(text)
    found = [t]
    if "=" in t:
        found.append(t.split("=", 1)[1])
    m = re.match(r"-[A-Za-z](.+)$", t)
    if m:
        found.append(m.group(1))
    for c in list(found):
        found += [p for p in re.split(r"[:,]", c) if p]
    out = []
    for c in found:
        out += brace_expand(c)
    return [c for c in dict.fromkeys(out) if c and not c.startswith("-")]


def resolve(cand, cwds):
    if cand.startswith("/"):
        return [cand]
    return [posixpath.join(cwd, cand) for cwd in cwds]


def word_level(text, cwds):
    best = None
    for cand in candidates(text):
        for joined in resolve(cand, cwds):
            level = path_level(joined)
            if level == "under":
                return level
            best = rank(best, level)
    return best


def all_under(text, cwds):
    """Whether a file name stands under a PHI root however the shell reads it."""
    cands = [c for c in candidates(text) if "=" not in c]
    return bool(cands) and all(path_level(j) == "under" for j in resolve(cands[0], cwds))


# ------------------------------------------------------------ shell reading

Word = collections.namedtuple("Word", "text dyn quoted")
Redir = collections.namedtuple("Redir", "fd op target dyn")
ASSIGN = re.compile(r"[A-Za-z_][A-Za-z0-9_]*=")


class Stage:
    def __init__(self):
        self.assigns, self.words, self.redirs, self.pipe_err = [], [], [], False

    def empty(self):
        return not (self.assigns or self.words or self.redirs)


def split_subst(cmd):
    """The command with each `$(...)`, backtick and `<(...)` replaced by SUB, and their bodies."""
    out, bodies = [], []
    i, n, single, double = 0, len(cmd), False, False
    while i < n:
        c = cmd[i]
        if single:
            single = c != "'"
            out.append(c)
        elif c == "\\" and i + 1 < n:
            out.append(c + cmd[i + 1])
            i += 1
        elif c == "'" and not double:
            single = True
            out.append(c)
        elif c == '"':
            double = not double
            out.append(c)
        elif c == "`":
            j = cmd.find("`", i + 1)
            if j < 0:
                raise ValueError("no closing backtick")
            bodies.append(cmd[i + 1 : j])
            out.append(SUB)
            i = j
        elif cmd.startswith("$(", i) or (not double and c in "<>" and cmd.startswith("(", i + 1)):
            depth, j = 1, i + 2
            while j < n and depth:
                depth += {"(": 1, ")": -1}.get(cmd[j], 0)
                j += 1
            if depth:
                raise ValueError("no closing parenthesis")
            bodies.append(cmd[i + 2 : j - 1])
            out.append(SUB)
            i = j - 1
        else:
            out.append(c)
        i += 1
    return "".join(out), bodies


def tokenize(s):
    """Tokens ('w', text, dyn, quoted), ('op', text), ('r', fd, op, target, dyn) and ('rp', fd, op)
    (a redirect whose target is the next word), and the text the shell ignores (comments,
    here-document bodies)."""
    toks, dropped, heredocs = [], [], []
    buf, quoted, dyn = [], False, False
    i, n = 0, len(s)

    def flush():
        nonlocal buf, quoted, dyn
        if buf or quoted:
            toks.append(("w", "".join(buf), dyn, quoted))
        buf, quoted, dyn = [], False, False

    while i < n:
        c = s[i]
        if c in " \t\r":
            flush()
        elif c == "\n":
            flush()
            toks.append(("op", "\n", None))
            i += 1
            for delim, strip in heredocs:
                while i < n:
                    j = s.find("\n", i)
                    line = s[i:] if j < 0 else s[i:j]
                    i = n if j < 0 else j + 1
                    if (line.lstrip("\t") if strip else line) == delim:
                        break
                    dropped.append(line)
            heredocs = []
            continue
        elif c == "#" and not buf and not quoted:
            j = s.find("\n", i)
            j = n if j < 0 else j
            dropped.append(s[i:j])
            i = j
            continue
        elif c == "'":
            j = s.find("'", i + 1)
            if j < 0:
                raise ValueError("no closing quotation")
            buf.append(s[i + 1 : j])
            quoted = True
            i = j + 1
            continue
        elif c == '"':
            j = i + 1
            while j < n and s[j] != '"':
                if s[j] == "\\" and j + 1 < n and s[j + 1] in '"\\$`\n':
                    if s[j + 1] != "\n":
                        buf.append(s[j + 1])
                    j += 2
                    continue
                if s[j] in "$" + SUB:
                    dyn = True
                buf.append(s[j])
                j += 1
            if j >= n:
                raise ValueError("no closing quotation")
            quoted = True
            i = j + 1
            continue
        elif c == "\\":
            if i + 1 < n and s[i + 1] != "\n":
                buf.append(s[i + 1])
                quoted = True
            i += 2
            continue
        elif c in ";|&()":
            flush()
            if c == ";":
                op = ";;" if s.startswith(";;", i) else ";"
            elif c == "|":
                op = "||" if s.startswith("||", i) else "|&" if s.startswith("|&", i) else "|"
            elif c == "&":
                if s.startswith("&>>", i):
                    toks.append(("rp", "both", ">>"))
                    i += 3
                    continue
                if s.startswith("&>", i):
                    toks.append(("rp", "both", ">"))
                    i += 2
                    continue
                op = "&&" if s.startswith("&&", i) else "&"
            else:
                op = c
            toks.append(("op", op, None))
            i += len(op)
            continue
        elif c in "<>":
            fd = None
            if buf and not quoted and "".join(buf).isdigit():
                fd = int("".join(buf))
                buf, quoted, dyn = [], False, False
            else:
                flush()
            for op in ("<<<", "<<-", "<<", "<&", "<>", "<") if c == "<" else (">>", ">&", ">|", ">"):
                if s.startswith(op, i):
                    break
            i += len(op)
            m = re.match(r"\d+-?|-", s[i:]) if op in (">&", "<&") else None
            if m:
                toks.append(("r", fd, op, m.group(0), False))
                i += m.end()
            elif op in ("<<", "<<-"):
                while i < n and s[i] in " \t":
                    i += 1
                m = re.match(r"""'([^']*)'|"([^"]*)"|\\?([^\s;|&<>()]+)""", s[i:])
                if not m:
                    raise ValueError("here-document without a delimiter")
                delim = next(g for g in m.groups() if g is not None)
                heredocs.append((delim, op == "<<-"))
                toks.append(("r", fd, op, delim, False))
                i += m.end()
            else:
                toks.append(("rp", fd, op))
            continue
        else:
            if c == "$" or c == SUB:
                dyn = True
            buf.append(c)
        i += 1
    flush()
    return toks, dropped


def parse(s):
    """(items, dropped text, groups): items are (stages, the operator after the pipeline)."""
    toks, dropped = tokenize(s)
    items, stages, cur, groups, k = [], [], Stage(), False, 0
    while k < len(toks):
        t = toks[k]
        if t[0] == "w":
            word = Word(t[1], t[2], t[3])
            if not cur.words and not t[3] and ASSIGN.match(t[1]):
                cur.assigns.append(word)
            else:
                cur.words.append(word)
        elif t[0] == "r":
            cur.redirs.append(Redir(t[1], t[2], t[3], t[4]))
        elif t[0] == "rp":
            k += 1
            if k >= len(toks) or toks[k][0] != "w":
                raise ValueError("a redirect without a target")
            cur.redirs.append(Redir(t[1], t[2], toks[k][1], toks[k][2]))
        else:
            op = t[1]
            if op in ("|", "|&"):
                if cur.empty():
                    raise ValueError("a pipe without a command")
                cur.pipe_err = op == "|&"
                stages.append(cur)
                cur = Stage()
            elif stages and cur.empty() and op == "\n":
                pass  # a pipe at the end of a line goes on in the next line
            else:
                groups = groups or op in ("(", ")")
                if not cur.empty():
                    stages.append(cur)
                if stages:
                    items.append((stages, op))
                stages, cur = [], Stage()
        k += 1
    if not cur.empty():
        stages.append(cur)
    if stages:
        items.append((stages, None))
    return items, dropped, groups


# ---------------------------------------------------------------- the rules

QUIET = {"cd", "pwd", "true", "false", ":", "test", "[", "[[", "sleep", "date"}  # print no file data
ECHO = {"echo", "printf"}  # print only what was typed, unless the shell expands a pattern
F1_VERBS = {"rsync", "cp", "mv", "ln", "mkdir", "trash"}
SEARCHERS = {"find", "du", "tree", "rg", "ag", "ack", "fd", "fdfind", "dust", "ncdu", "gdu"}
GREPS = {"grep", "egrep", "fgrep", "zgrep", "zegrep"}
PATTERN_FIRST = GREPS | {"rg", "ag", "ack", "fd", "fdfind"}
SUMMARY = {
    "number of files",
    "number of created files",
    "number of deleted files",
    "number of regular files transferred",
    "total file size",
    "total transferred",
    "total transferred file size",
    "total bytes sent",
    "total bytes received",
    "total size is",
    "literal data",
    "matched data",
    "file list size",
}
RSYNC_LETTERS = set("arltpogDHLkKuczhnqxXASWmRdEyIOUCJNF46s")
RSYNC_VALUE_LETTERS = set("efBT")
RSYNC_FLAGS = set(
    "archive recursive quiet compress human-readable dry-run checksum stats update times perms"
    " owner group links devices specials protect-args mkpath size-only ignore-existing existing"
    " prune-empty-dirs whole-file no-whole-file sparse one-file-system no-motd partial inplace"
    " copy-links copy-dirlinks keep-dirlinks hard-links numeric-ids append append-verify"
    " omit-dir-times no-perms no-owner no-group no-times no-links dirs relative no-relative"
    " executability acls xattrs ignore-times fuzzy delay-updates safe-links copy-unsafe-links"
    " preallocate cvs-exclude from0 secluded-args super fake-super ipv4 ipv6 mkpath".split()
)
RSYNC_VALUES = set(
    "exclude include exclude-from include-from files-from filter rsh partial-dir temp-dir bwlimit"
    " timeout contimeout max-size min-size block-size chmod chown usermap groupmap link-dest"
    " copy-dest compare-dest backup-dir suffix modify-window port address info skip-compress"
    " checksum-choice compress-level".split()
)
PY_EXE = re.compile(r"^python(\d+(\.\d+)*)?$")


def verb_of(stage):
    if not stage.words:
        return ""
    word = stage.words[0].text
    d, b = posixpath.split(word)
    return word if d and d not in SYSTEM_BIN else b


def recursive(verb, args):
    """Whether the command searches below the folders it is given."""
    if verb in SEARCHERS:
        return True
    letters = "".join(a[1:] for a in args if a.startswith("-") and not a.startswith("--"))
    longs = [a[2:].split("=")[0] for a in args if a.startswith("--")]
    if verb in GREPS:
        return bool(set(letters) & set("rR")) or any(x in longs for x in ("recursive", "dereference-recursive"))
    if verb == "ls":
        return "R" in letters or "recursive" in longs
    if verb == "diff":
        return "r" in letters or "recursive" in longs
    return False


def path_operands(verb, args):
    ops = [a for a in args if not a.startswith("-")]
    if verb in PATTERN_FIRST:
        return ops[1:]
    if verb == "find":
        first = []
        for a in args:
            if a.startswith("-") or a in ("(", "!", ")"):
                break
            first.append(a)
        return first
    return ops


def stage_texts(stage, typed_text_ok=False):
    """The words of a stage that may be paths. With typed_text_ok (echo, printf) a word counts only
    when it is written as an absolute path: the others are text, not names to find in the folder."""
    texts = [a.text for a in stage.assigns]
    for k, w in enumerate(stage.words):
        if typed_text_ok and k and not w.text.startswith("/"):
            continue
        if k or "/" in w.text:
            texts.append(w.text)
    for r in stage.redirs:
        if r.target and r.op not in ("<<", "<<-") and not (r.op in (">&", "<&") and re.match(r"[\d-]", r.target)):
            texts.append(r.target)
    return texts


def above_gap(joined):
    """How many folder levels lie between this path and the nearest PHI folder below it, or None."""
    best = None
    for p in {posixpath.normpath(joined), real(joined)}:
        segs = [s for s in p.lower().split("/") if s]
        for entry in ENTRIES:
            if rel(segs, entry) == "above":
                gap = len(entry) - len(segs)
                best = gap if best is None else min(best, gap)
    return best


def depth_cap(verb, args):
    """How many levels below its folders a searching command goes, or None when it has no limit."""
    num = lambda x: int(x) if x.isdigit() else None  # noqa: E731
    cap = None
    for i, a in enumerate(args):
        nxt = args[i + 1] if i + 1 < len(args) else ""
        if verb == "find" and a == "-maxdepth":
            cap = num(nxt)
        elif verb == "du":
            letters = a[1:] if a.startswith("-") and not a.startswith("--") else ""
            if a == "--summarize" or "s" in letters:
                cap = 0
            elif a in ("-d", "--max-depth"):
                cap = num(nxt)
            elif a.startswith("--max-depth="):
                cap = num(a.split("=", 1)[1])
        elif verb == "tree" and a == "-L":
            cap = num(nxt)
    return cap


def stage_level(stage, cwds):
    """'under' if the stage names PHI, 'cwd' if it only runs in a PHI folder, 'above' if it searches
    below a folder that holds PHI."""
    verb = verb_of(stage)
    args = [w.text for w in stage.words[1:]]
    if verb in QUIET:
        return None
    plain_echo = verb in ECHO and not any(GLOB.search(w.text) or SUB in w.text for w in stage.words[1:])
    level = None
    if any(TEXT_RE.search(t) for t in stage_texts(stage)):
        return "under"
    for t in stage_texts(stage, typed_text_ok=plain_echo):
        level = rank(level, word_level(t, cwds))
    if level == "under":
        return level
    if plain_echo:
        return None  # prints only what was typed
    cwd_levels = [path_level(c) for c in cwds]
    if verb and "under" in cwd_levels:
        return "cwd"
    if verb in F1_VERBS and level == "above":
        return "under"  # a folder that holds PHI folders is copied as a whole: check where it goes
    if verb and recursive(verb, args):
        operands = path_operands(verb, args)
        if not operands:
            for lv in cwd_levels:
                level = rank(level, lv)
        if level != "above":
            return None
        cap = depth_cap(verb, args)
        if cap is not None:  # a limit that stops at or above the PHI folder lists no content
            folders = [j for w in operands for c in candidates(w) for j in resolve(c, cwds)] or list(cwds)
            gaps = [g for g in (above_gap(j) for j in folders) if g is not None]
            if all(g >= cap for g in gaps):
                return None
        return "above"
    return None


def fd_state(stage, out1, out2, cwds):
    """Where stdout and stderr go once the redirects are applied in order."""
    state = {1: out1, 2: out2}

    def kind(target, dyn):
        if dyn:
            return "other"
        if target == "/dev/null":
            return "null"
        return "phi" if all_under(target, cwds) else "other"

    for r in stage.redirs:
        if r.op in ("<", "<<", "<<<", "<&", "<>"):
            continue
        if r.fd == "both" or (r.fd is None and r.op == ">&" and not re.match(r"[\d-]", r.target)):
            state[1] = state[2] = kind(r.target, r.dyn)
            continue
        fd = 1 if r.fd is None else r.fd
        if fd not in state:
            continue
        if r.op == ">&" and re.match(r"[\d-]", r.target):
            src = r.target.rstrip("-")
            state[fd] = state.get(int(src), "other") if src.isdigit() else "null"
        else:
            state[fd] = kind(r.target, r.dyn)
    return state


def is_sanitizer(word, cwds):
    wanted = {s.lower() for s in SANITIZERS}
    return any(posixpath.normpath(j).lower() in wanted for j in resolve(tilde(word), cwds))


def python_script(stage):
    """(script, its arguments) for `python3 -I script.py ...`, else None."""
    ws = stage.words
    if not ws or not PY_EXE.match(posixpath.basename(ws[0].text)):
        return None
    i, isolated = 1, False
    while i < len(ws) and ws[i].text.startswith("-") and len(ws[i].text) > 1:
        if ws[i].text == "-I":
            isolated = True
        elif ws[i].text not in ("-u", "-B", "-O", "-OO", "-s", "-S", "-E", "-q"):
            return None
        i += 1
    if not isolated or i >= len(ws) or ws[i].dyn:
        return None
    return ws[i].text, ws[i + 1 :]


def hint(verb):
    if verb in ("cat", "head", "tail", "less", "more", "sed", "awk", "zcat", "xxd", "od", "strings"):
        return f"`{verb}` prints patient file content."
    if verb in ("ls", "find", "du", "tree", "stat", "file", "wc", "fd"):
        return f"`{verb}` prints names, sizes or counts from a patient folder (a BAM name is a sample ID)."
    if verb in ("grep", "egrep", "fgrep", "rg", "ag", "ack", "zgrep"):
        return f"`{verb}` prints lines of patient files."
    if verb in ("rm", "rmdir", "unlink"):
        return f"`{verb}` is not on the PHI list; move the files to the Trash with `trash`."
    if verb in ("git", "gh"):
        return f"`{verb}` may not name a PHI path. Commit and issue text goes to GitHub: keep PHI paths out of it (git commit -F, gh --body-file)."
    return "This command names a PHI path and is not one of the allowed forms."


def no_subst(stage, what):
    if any(SUB in w.text for w in stage.words) or any(SUB in a.text for a in stage.assigns):
        deny(f"{what} may not use a command substitution")


def check_f2(first, second, last, cwds):
    run = python_script(first)
    if run is None or not run[0].endswith(".py"):
        deny(hint(verb_of(first)) + " Only a script run (python3 -I <script>.py ...) may feed the sanitizer.")
    no_subst(first, "the script call")
    run = python_script(second)
    if run is None or not is_sanitizer(run[0], cwds):
        deny("The last stage must be the sanitizer: python3 -I <sanitizer path>. Nothing else may follow the script.")
    if second.redirs or any(w.dyn for w in second.words) or any(a.dyn for a in second.assigns):
        deny("The sanitizer stage may not have a redirect, a variable or a substitution.")
    if not last:
        deny("Nothing may follow the sanitizer: no chain, pipe, redirect or second command after it.")
    state = fd_state(first, "pipe", "pipe" if first.pipe_err else "tty", cwds)
    for fd, name in ((1, "output"), (2, "error output")):
        if state[fd] not in ("pipe", "null", "phi"):
            deny(f"the script's {name} must reach the sanitizer: end the script call with 2>&1")


def check_f3(stage, cwds):
    no_subst(stage, "a command that names patient data")
    state = fd_state(stage, "tty", "tty", cwds)
    if any(state[fd] not in ("null", "phi") for fd in (1, 2)):
        deny(hint(verb_of(stage)) + " Its output (stdout and stderr) would reach the screen or a file outside the PHI folders.")


def option_fault(verb, a):
    return f"{verb} option {a} prints file names or is not on the PHI list"


def split_cp(verb, args):
    """(positional words, the -t folder or None) for cp, mv and ln; deny the options that print."""
    pos, target, ended, i = [], None, False, 0
    while i < len(args):
        a = args[i]
        if ended or not a.startswith("-") or a == "-":
            pos.append(a)
        elif a == "--":
            ended = True
        elif a.startswith("--"):
            name, eq, val = a[2:].partition("=")
            if len(name) >= 2 and any(x.startswith(name) for x in ("verbose", "interactive", "debug")):
                deny(option_fault(verb, a))
            if len(name) >= 2 and ("target-directory".startswith(name) or "suffix".startswith(name)):
                if not eq:
                    i += 1
                    val = args[i] if i < len(args) else ""
                if "target-directory".startswith(name):
                    target = val
        else:
            for j, ch in enumerate(a[1:], 1):
                if ch in "vi":
                    deny(option_fault(verb, a))
                if ch in "tS":
                    val = a[j + 1 :]
                    if not val:
                        i += 1
                        val = args[i] if i < len(args) else ""
                    if ch == "t":
                        target = val
                    break
        i += 1
    return pos, target


def check_rsh(value):
    words = shlex.split(value) if value else []
    if not words or posixpath.basename(words[0]) != "ssh":
        deny("rsync -e may run ssh only")


def check_rsync(args):
    """The positional words of an rsync call; deny options that print names, run a command or delete."""
    pos, ended, i = [], False, 0
    while i < len(args):
        a = args[i]
        if ended or not a.startswith("-") or a == "-":
            pos.append(a)
        elif a == "--":
            ended = True
        elif a.startswith("--"):
            name, eq, val = a[2:].partition("=")
            if name in RSYNC_FLAGS and not eq:
                pass
            elif name in RSYNC_VALUES:
                if not eq:
                    i += 1
                    val = args[i] if i < len(args) else ""
                if name == "rsh":
                    check_rsh(val)
                elif name == "info" and not re.fullmatch(r"stats\d?(,stats\d?)*", val):
                    deny("rsync --info may ask for stats only (other values print file names)")
            else:
                deny(option_fault("rsync", a))
        else:
            for j, ch in enumerate(a[1:], 1):
                if ch in RSYNC_LETTERS:
                    continue
                if ch not in RSYNC_VALUE_LETTERS:
                    deny(option_fault("rsync", a))
                val = a[j + 1 :]
                if not val:
                    i += 1
                    val = args[i] if i < len(args) else ""
                if ch == "e":
                    check_rsh(val)
                break
        i += 1
    return pos


def is_remote(text):
    return bool(re.match(r"^(rsync://|[\w.@-]+:)", text)) and not text.startswith(("/", "."))


def check_f1(stage, cwds, piped=False):
    """A file operation that prints no content: safe options, patient data goes only to
    another PHI folder or another host, and any redirect goes to /dev/null or a PHI folder."""
    verb = verb_of(stage)
    if any(w.dyn for w in stage.words) or any(a.dyn for a in stage.assigns):
        deny(f"{verb} on patient data may not use a variable or a substitution")
    if any(r.op in ("<<", "<<-", "<<<") for r in stage.redirs):
        deny(f"{verb} on patient data takes no here-document")
    state = fd_state(stage, "pipe" if piped else "tty", "pipe" if stage.pipe_err else "tty", cwds)
    if any(v == "other" for v in state.values()):
        deny(f"{verb} on patient data may not redirect into a file outside the PHI folders")
    args = [w.text for w in stage.words[1:]]
    if verb == "trash":
        if any(a.startswith("-") for a in args if a != "--"):
            deny("trash takes paths only (-l lists the Trash, -e empties it)")
        return
    if verb == "mkdir":
        if any(a == "--verbose" or (a.startswith("-") and not a.startswith("--") and "v" in a) for a in args):
            deny(option_fault(verb, "-v"))
        return
    if verb == "rsync":
        pos, target = check_rsync(args), None
        if len(pos) < 2:
            deny("rsync with one path lists it; name a source and a destination")
    else:
        pos, target = split_cp(verb, args)
        if not pos:
            deny(f"{verb} needs a source and a destination")
        if verb == "ln" and len(pos) == 1 and target is None:
            target = "."
    dest = target if target is not None else pos[-1]
    sources = pos if target is not None else pos[:-1]
    levels = [word_level(s, cwds) for s in sources]
    if "under" in levels or "above" in levels:
        # patient data stays in a PHI folder; a folder that merely holds PHI folders may go to one like it
        ok = ("under",) if "under" in levels else ("under", "above")
        if not (word_level(dest, cwds) in ok or (verb == "rsync" and is_remote(dest))):
            deny(
                f"{verb} would put patient data outside the PHI folders; copy it into another"
                " PHI folder, or add the destination to the PHI list first"
            )


def check_f1_filter(first, second, cwds):
    """`rsync ... | tail -N` and `rsync ... | grep -E 'summary line|summary line'`."""
    check_f1(first, cwds, piped=True)
    if second.redirs or any(w.dyn for w in second.words) or any(a.dyn for a in second.assigns):
        deny("the filter after rsync takes no redirect, variable or substitution")
    verb, args = verb_of(second), [w.text for w in second.words[1:]]
    if verb == "tail":
        one = re.fullmatch(r"-n?(\d+)", args[0]) if len(args) == 1 else None
        two = len(args) == 2 and args[0] == "-n" and args[1].isdigit()
        if not args or (one and int(one.group(1)) <= 25) or (two and int(args[1]) <= 25):
            return
        deny("after rsync, tail may take -n N (N up to 25) and no file")
    if verb == "grep" and len(args) >= 2 and all(re.fullmatch(r"-[Eihe]+", a) for a in args[:-1]):
        alts = [x.strip().lower() for x in re.split(r"\\?\|", args[-1])]
        if alts and all(a in SUMMARY for a in alts):
            return
    deny("after rsync, only tail or grep -E on its summary lines (Number of files, Total transferred) may follow")


def check_pipeline(stages, last, cwds):
    n = len(stages)
    verbs = [verb_of(s) for s in stages]
    if n == 1:
        if verbs[0] in F1_VERBS:
            check_f1(stages[0], cwds)
        else:
            check_f3(stages[0], cwds)
    elif n == 2 and verbs[0] == "rsync" and verbs[1] in ("tail", "grep"):
        check_f1_filter(stages[0], stages[1], cwds)
    elif n == 2:
        check_f2(stages[0], stages[1], last, cwds)
    else:
        deny("a pipe of more than two commands that names patient data is not allowed")


def cd_targets(stage, cwds):
    args = [w for w in stage.words[1:] if w.text not in ("-L", "-P", "--")]
    if not args:
        return [HOME]
    if len(args) != 1 or args[0].dyn or args[0].text == "-":
        return None
    out = []
    for joined in resolve(tilde(args[0].text), cwds):
        for p in (posixpath.normpath(joined), real(joined)):
            if p not in out:
                out.append(p)
    return out


def analyze(text, cwds, depth=0, strict=False):
    """Raise Deny unless the command is allowed; with strict, any PHI mention is denied."""
    if depth > 3:
        if TEXT_RE.search(text):
            deny("substitutions nested too deep to check")
        return
    flat, bodies = split_subst(text)
    items, dropped, groups = parse(flat)
    cwds0 = list(cwds)
    if dropped and TEXT_RE.search("\n".join(dropped)):
        deny("a comment or here-document names a PHI path; write a script that holds PHI paths with the Write tool, not inside a command")
    for idx, (stages, op) in enumerate(items):
        levels = [stage_level(s, cwds) for s in stages]
        if any(levels):
            if strict:
                deny("a command substitution may not name patient data")
            if "above" in levels and "under" not in levels:
                names = ", ".join(sorted({verb_of(s) for s, lv in zip(stages, levels) if lv == "above"}))
                deny(
                    f"`{names}` searches below a folder that holds a PHI folder, so it would list or read patient data."
                    " Search a folder that holds none, or stop above the PHI folder (find -maxdepth, du -s or -d, tree -L)."
                )
            if groups:
                deny("a group or subshell that names patient data is not supported; run simple commands")
            if op == "&":
                deny("patient data may not be used in a background command")
            try:
                check_pipeline(stages, idx == len(items) - 1, cwds)
            except Deny as d:
                if "cwd" in levels:
                    raise Deny(f"{d} This call runs in a PHI folder (an earlier `cd`): `cd` to a folder outside the PHI list first.")
                raise
        if len(stages) == 1 and verb_of(stages[0]) == "cd":
            targets = cd_targets(stages[0], cwds)
            if targets is not None and op in ("&&", ";", "\n", ";;"):
                cwds = targets if op == "&&" else cwds + [t for t in targets if t not in cwds]
    for body in bodies:
        analyze(body, cwds0 + [c for c in cwds if c not in cwds0], depth + 1, strict=True)


# ------------------------------------------------------------------- tools


def bash_reason(cmd, cwd):
    text = expand_home(cmd)
    cwds = [cwd] + ([real(cwd)] if real(cwd) != cwd else [])
    try:
        analyze(text, cwds)
    except Deny as d:
        return str(d)
    except Exception as e:  # unreadable shell syntax or a fault: fail closed when PHI may be involved
        if TEXT_RE.search(text) or any(path_level(c) == "under" for c in cwds):
            kind = "reads shell syntax it cannot follow" if isinstance(e, ValueError) else f"faulted ({type(e).__name__})"
            return f"the command may name patient data and the guard {kind}; split it into simple commands."
    return None


def text_of(value):
    return value if isinstance(value, str) else ""


def tool_reason(tool, inp, cwd):
    """Why a Read, Grep, Glob or NotebookRead call is refused, or None."""
    searching = tool in ("Grep", "Glob")
    base = next((text_of(inp.get(k)) for k in ("path", "file_path", "notebook_path") if inp.get(k)), "")
    pattern = text_of(inp.get("pattern"))
    paths = []
    if tool == "Glob" and pattern:
        # the pattern says where Glob looks: a folder above PHI is fine while the pattern stays out
        paths.append(posixpath.join(expand_home(base) if base else cwd, expand_home(pattern)))
    elif base or searching:
        paths.append(base or cwd)
    if tool == "Grep" and text_of(inp.get("glob")).startswith(("/", "~")):
        paths.append(text_of(inp.get("glob")))
    for p in paths:
        text = expand_home(p)
        if TEXT_RE.search(text):
            return f"`{tool}` on a PHI path shows patient content or names."
        level = None
        for cand in candidates(text) or [text]:
            for joined in resolve(cand, [cwd]):
                level = rank(level, path_level(joined))
        if level == "under" or (level == "above" and tool == "Grep"):  # a Glob pattern does not go below its last name
            where = "searches below this folder, which holds a PHI path, and shows" if level == "above" else "on a PHI path shows"
            return f"`{tool}` {where} patient content or names."
    return None


def check(tool, inp, cwd):
    """The reason to refuse this tool call, or None to allow it."""
    if tool not in ("Bash", "Read", "Grep", "Glob", "NotebookRead"):
        return None  # Write and Edit show only what the model wrote; the tool refuses unread files
    if PROBLEM:
        return f"the PHI list is unusable ({PROBLEM}); the operator must fix it first."
    if TEXT_RE is None:
        return None  # no list on this machine: nothing to guard
    try:
        if tool == "Bash":
            return bash_reason(text_of(inp.get("command")), cwd)
        return tool_reason(tool, inp, cwd)
    except Deny as d:
        return str(d)
    except Exception as e:  # a fault in the hook fails closed for what names a PHI path
        seen = expand_home(" ".join(str(v) for v in inp.values()))
        return f"the guard could not check this call ({type(e).__name__})." if TEXT_RE.search(seen) else None


HOW_TO = (
    "Tool output goes to the model provider, so a tool may not print patient data, nor file names"
    " or paths from a PHI folder. Do one of these instead:\n"
    "1. Run an aggregate script and send ALL its output through the sanitizer, as the last stage:"
    " python3 -I <script>.py ... 2>&1 | python3 -I {san} [--names-from FILE]\n"
    "2. Send ALL output (stdout and stderr) to a file under a PHI folder, with no pipe or chain:"
    " <command> > <PHI folder>/out.txt 2>&1 (the operator reviews that file; never read it with a tool).\n"
    "3. Move data without printing: rsync, cp, mv, ln, mkdir, trash (no -v, -P or -i). For rsync totals add"
    " --stats 2>&1 | grep -E 'Number of files|Total transferred'.\n"
    "Rule: .agents/rules/security.md (Patient BAM slices). The PHI list is ~/.claude/phi-paths.txt."
)


def message(reason):
    san = SANITIZERS[-1].replace(HOME, "~", 1) if SANITIZERS[-1].startswith(HOME + "/") else SANITIZERS[-1]
    return "PHI guard: " + reason + " " + HOW_TO.replace("{san}", san)


def emit(reason):
    if reason is None:
        print("{}")
        return
    print(json.dumps({"hookSpecificOutput": {
        "hookEventName": "PreToolUse",
        "permissionDecision": "deny",
        "permissionDecisionReason": message(reason),
    }}))


def main():
    try:
        data = json.load(sys.stdin)
    except Exception:
        data = None
    if not isinstance(data, dict):
        print("{}")
        return
    inp = data.get("tool_input")
    cwd = data.get("cwd")
    cwd = posixpath.abspath(cwd) if isinstance(cwd, str) and cwd else os.getcwd()
    emit(check(data.get("tool_name"), inp if isinstance(inp, dict) else {}, cwd))


if __name__ == "__main__":
    main()
