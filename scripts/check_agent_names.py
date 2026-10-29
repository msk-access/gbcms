#!/usr/bin/env python3
"""Fail when an agent skill or rule names code, a flag or a file that does not exist.

The skills (`.agents/skills/*/SKILL.md`) and the rules (`.agents/rules/*.md`) are read
on trust, so a renamed function or a deleted file stays in them until a reader acts on
it. This checks every backticked name:

* a path (`rust/src/shared/fragment.rs`, `engine.rs`, `tests/test_x.py::TestY`) must
  exist in the tree; a bare file name matches a tracked file by its name, a partial path
  by its tail;
* a command-line flag (`--umi-tag`, `--no-baq`) must appear in `src/` or the Nextflow
  files, as a complete name (`--mfsd-report-min` is not found in `--mfsd-report-min-alt`);
* an identifier (`check_allele_with_qual`, `window::allele_kind`, `BIN_WINDOW`, a CamelCase
  type, a `snake_case` field) must appear as a whole word in `src/` or `rust/src/`
  (the Nextflow and read-filters skills also search `nextflow/`);
* `NAME = value` for an upper-case constant must also carry that value in the source;
* the module map in `architecture.md` must list every Rust and Python module file by its
  full path (`counting/engine.rs` and `normalize/engine.rs` are two files), and only files
  that exist.

Plain words (`reverse`, `sw`), shell commands and phrases are not identifiers and are not
checked. `EXTERNAL` lists names that belong to another tool, each with its reason.
Fenced code blocks are skipped, except the module map.

Usage: python3 scripts/check_agent_names.py [--staged] [FILE ...]   (default: every skill and rule)
With --staged (the pre-commit hook) the files, the code they are compared with and the
module map are read from the git index, not the working tree, so the check judges what the
commit will hold. Exits 1 and prints one line per unresolved name.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: Names that belong to another tool or language, not to this code base.
EXTERNAL = {
    # Nextflow strict-syntax advice names Groovy/Nextflow constructs, not our code.
    "Channel.x",
    "channel.x",
    "asBool",  # nextflow/modules/local/utils/main.nf; the Nextflow roots hold it
    "resourceLimits",
    "includeConfig",
    "hasData",
    # Build and environment advice names other tools and their settings.
    "OPENSSL_NO_VENDOR",
    "VIRTUAL_ENV",
    "CONDA_PREFIX",
    "RUST_LOG",
    "NXF_VER",
    "pyo3-log",
    "manifest-path",
    "python-source",
    "module-name",
    "coitrees",
    "IntervalNode",
    "Interval",
    "Borrow",
    # RNA-protocol advice names external tools.
    "infer_experiment",
    # Names of the test-suite helpers appear in tests/, which is not searched.
    "count_checked",
    "count_bam_checked",
    "find_read_pos",
    # Git's own flags.
    "--no-ff",
    # Exists only on checkpoint-wheels/* branches, where the pre-commit hook admits it.
    "scripts/regression_panel/checkpoints.tsv",
}

#: Extra directories searched for a skill whose names live outside src/ and rust/src/.
EXTRA_ROOTS = {"nextflow-pipeline": ["nextflow"], "read-filters-qc": ["nextflow"]}
#: Files that exist but are not module-map entries (every package has them).
MAP_IGNORED = {"__init__.py", "__main__.py", "py.typed", "mod.rs"}
#: A path-like span starts with one of these top-level names, or ends in a file suffix.
TOP_LEVEL = (
    "src/",
    "rust/",
    "docs/",
    "tests/",
    "nextflow/",
    "scripts/",
    ".agents/",
    ".claude/",
    ".githooks/",
    ".github/",
)
SUFFIXES = (".rs", ".py", ".pyi", ".md", ".nf", ".toml", ".config", ".yml", ".yaml")

WORD = r"[A-Za-z_][A-Za-z0-9_]*"
SINGLE = re.compile(rf"^(?P<path>(?:{WORD}(?:::|\.))*{WORD})(?:\(.*\))?$")
CONST_VALUE = re.compile(r"^([A-Z][A-Z0-9_]+) = (\S+)$")
FENCE = re.compile(r"```.*?```", re.S)
SPAN = re.compile(r"`([^`\n]+)`")


CODE_SUFFIXES = (".rs", ".py", ".pyi", ".nf", ".config", ".toml")


class Source:
    """Where file names and contents come from: the working tree, or the git index.

    The index is what a commit will hold, so the pre-commit hook (`--staged`) reads it: a
    function added only in the working tree must not satisfy a skill that is committed
    without it, nor a skill fixed only in the working tree pass for the staged one.
    """

    def __init__(self, staged: bool = False) -> None:
        self.staged = staged
        self._blobs: dict[str, str] = {}
        self._names = self._list()

    def _list(self) -> list[str]:
        if self.staged:  # git is required: the index is what is read
            return self._git(["ls-files", "--cached"]).decode().splitlines()
        try:  # tracked files plus untracked files that are not ignored
            out = self._git(["ls-files", "--cached", "--others", "--exclude-standard"])
            return out.decode().splitlines()
        except (OSError, subprocess.CalledProcessError):
            # Without git (a source archive) the tree is walked, skipping build output.
            skip = {".git", "target", ".venv", "venv", "__pycache__", "node_modules"}
            return [
                str(p.relative_to(ROOT))
                for p in ROOT.rglob("*")
                if p.is_file() and not skip.intersection(p.relative_to(ROOT).parts)
            ]

    @staticmethod
    def _git(args: list[str], stdin: bytes | None = None) -> bytes:
        return subprocess.run(
            ["git", *args], cwd=ROOT, input=stdin, capture_output=True, check=True
        ).stdout

    @property
    def names(self) -> list[str]:
        return self._names

    def exists(self, rel: str) -> bool:
        """Whether a path is a file here (in the index when staged)."""
        return rel in set(self._names) if self.staged else (ROOT / rel).exists()

    def preload(self, files: list[str]) -> None:
        """Read the index blobs of `files` in one `git cat-file --batch` call."""
        if not self.staged:
            return
        todo = [f for f in files if f not in self._blobs]
        if not todo:
            return
        wanted = set(todo)
        listing = self._git(["ls-files", "-s", "-z"]).split(b"\0")
        oids = {}
        for entry in listing:
            meta, _, name = entry.partition(b"\t")
            if name.decode(errors="replace") in wanted:
                oids[name.decode(errors="replace")] = meta.split()[1].decode()
        raw = self._git(["cat-file", "--batch"], "".join(f"{o}\n" for o in oids.values()).encode())
        pos = 0
        for name in oids:
            end = raw.index(b"\n", pos)
            size = int(raw[pos:end].split()[2])
            self._blobs[name] = raw[end + 1 : end + 1 + size].decode(errors="replace")
            pos = end + 1 + size + 1

    def read(self, rel: str) -> str:
        """The text of a file here; empty when it is not."""
        if self.staged:
            self.preload([rel])
            return self._blobs.get(rel, "")
        try:
            return (ROOT / rel).read_text(errors="replace")
        except OSError:
            return ""

    def read_tree(self, roots: list[str]) -> str:
        """The text of every .rs, .py, .pyi, .nf, .config and .toml file under the roots."""
        chosen = [
            n
            for n in self._names
            if n.endswith(CODE_SUFFIXES)
            and any(n == r or n.startswith(r.rstrip("/") + "/") for r in roots)
        ]
        self.preload(chosen)
        return "\n".join(self.read(n) for n in chosen)


def is_identifier_like(word: str) -> bool:
    """A name worth looking up: it has an underscore, a CamelCase hump or is a constant."""
    if word.startswith("_") and len(word) > 1:
        return True
    if "_" in word.strip("_"):
        return True
    if re.search(r"[a-z][A-Z]", word):
        return True
    return word.isupper() and len(word) >= 4


def path_exists(span: str, source: Source) -> bool:
    """Whether a path-like span names a tracked file or directory."""
    names = source.names
    path = span.split("::")[0]
    if "*" in path or "<" in path or "$" in path or "{" in path:
        pattern = re.escape(path).replace(r"\*", ".*")
        pattern = re.sub(r"\\<.*?\\>", ".*", pattern)
        pattern = re.sub(r"\\\$\w+", ".*", pattern)
        return any(re.search(pattern, n) for n in names)
    if path.endswith("/"):
        return any(n.startswith(path) or ("/" + path) in n for n in names)
    return any(n == path or n.endswith("/" + path) for n in names) or source.exists(path)


def has_flag(flag: str, text: str) -> bool:
    """Whether `text` holds the complete flag, not a longer flag that starts with it."""
    return re.search(rf"(?<![\w-]){re.escape(flag)}(?![\w-])", text) is not None


def check_span(span: str, source: Source, text: str, flag_text: str) -> str | None:
    """Return why the span's name does not resolve, or None when it does or is not a name."""
    span = span.strip()
    if span.startswith("--"):
        flag = re.match(r"--[A-Za-z0-9][A-Za-z0-9-]*", span)
        if flag is None or flag.group(0) in EXTERNAL or span.endswith("-"):
            return None
        if "<" in span or " false" in span:
            return None
        return (
            None
            if has_flag(flag.group(0), flag_text)
            else f"flag {flag.group(0)} is not in the code"
        )
    # A file or directory: a known top-level directory, a file suffix, or a bare file name.
    head = span.split("::")[0].split(" ")[0]
    # Build products (the compiled extension) are not tracked; shared-library names are not paths.
    if head in EXTERNAL or head.startswith("~") or head.endswith(".so") or "..." in head:
        return None
    if (
        (head.endswith(SUFFIXES) or head.startswith(TOP_LEVEL))
        and re.fullmatch(r"[\w./*<>${}-]+", head)
        and not head.startswith(("-", "."))
    ):
        if not path_exists(span, source):
            return f"path {head} does not exist"
        if "::" in span and head.endswith(".py"):
            tail = span.split("::", 1)[1].split("::")[0].split("[")[0]
            if not re.search(rf"\b{re.escape(tail)}\b", source.read(head)):
                return f"{tail} is not in {head}"
        return None
    const = CONST_VALUE.match(span)
    if const and is_identifier_like(const.group(1)):
        name, value = const.groups()
        if not re.search(rf"\b{name}\b[^=\n]*=\s*{re.escape(value)}\b", text):
            return f"constant {name} is not {value} in the code"
        return None
    # `rna.rs::read_transcript_strand(record, X)` and `hash_molecule(qname, umi)`:
    # look up the callee. A phrase with a space before any bracket is not a name.
    m = SINGLE.match(span.split(" ")[0]) if "(" in span.split(" ")[0] else SINGLE.match(span)
    if m is None:
        return None
    full = m.group("path")
    if full in EXTERNAL or any(part in EXTERNAL for part in re.split(r"::|\.", full)):
        return None
    missing = []
    for part in re.split(r"::|\.", full):
        if part in EXTERNAL or f".{part}" in SUFFIXES or part in {"rs", "py", "nf", "toml"}:
            continue
        if not (is_identifier_like(part) or "(" in span):
            continue
        if not re.search(rf"\b{re.escape(part)}\b", text):
            missing.append(part)
    # A bare method call such as `.readLines()` belongs to another language.
    if span.startswith("."):
        return None
    return f"{', '.join(missing)} not found in the code" if missing else None


def map_paths(block: str) -> set[str]:
    """The full path of every file in the module-map tree, rebuilt from its indentation.

    Each level is four columns (`│   ` or four spaces) before `├── ` or `└── `; a name ending
    in `/` opens a directory (it may span segments, as `src/gbcms/`), and a line may list
    several files (`mod.rs, engine.rs`). Text after `#` is a comment.
    """
    stack: list[str] = []
    found: set[str] = set()
    for line in block.splitlines():
        entry = re.match(r"((?:│   |    )*)[├└]── (.*)$", line)
        if entry is None:
            continue
        depth = len(entry.group(1)) // 4
        body = entry.group(2).split("#")[0].strip()
        del stack[depth:]
        if body.endswith("/"):
            stack.append(body)
            continue
        for token in re.split(r",\s*", body):
            if re.fullmatch(r"[\w.]+\.(?:rs|py|pyi)", token):
                found.add("".join(stack) + token)
    return found


def module_map_problems(text: str, names: list[str]) -> list[str]:
    """Compare the module-map tree in architecture.md with the Rust and Python files.

    Full paths are compared: `counting/engine.rs` and `normalize/engine.rs` are two files,
    and neither may stand in for the other.
    """
    block = re.search(r"## Project Layout\s+```.*?\n(.*?)```", text, re.S)
    if block is None:
        return ["no module map under '## Project Layout'"]

    def ignored(path: str) -> bool:
        return path.rsplit("/", 1)[-1] in MAP_IGNORED

    listed = {f for f in map_paths(block.group(1)) if not ignored(f)}
    actual = {
        n
        for n in names
        if (n.startswith("rust/src/") and n.endswith(".rs"))
        or (n.startswith("src/gbcms/") and n.endswith((".py", ".pyi")))
    }
    actual = {f for f in actual if not ignored(f)}
    problems = [f"module map lists {f}, which is not a file" for f in sorted(listed - actual)]
    problems += [f"module map omits {f}" for f in sorted(actual - listed)]
    return problems


def read_doc(path: Path, source: Source) -> str:
    """A skill or rule's text: from the index when staged (a path outside the repo is read as is)."""
    if source.staged and path.is_relative_to(ROOT):
        return source.read(path.relative_to(ROOT).as_posix())
    return path.read_text()


def main(argv: list[str]) -> int:
    staged = "--staged" in argv
    files = [a for a in argv if a != "--staged"]
    source = Source(staged)
    if files:
        paths = [Path(a).resolve() for a in files]
    elif staged:
        paths = [
            ROOT / n
            for n in source.names
            if re.fullmatch(r"\.agents/(skills/[^/]+/SKILL\.md|rules/[^/]+\.md)", n)
        ]
    else:
        paths = sorted(ROOT.glob(".agents/skills/*/SKILL.md")) + sorted(
            ROOT.glob(".agents/rules/*.md")
        )
    code = source.read_tree(["src", "rust/src"])
    flag_text = code + "\n" + source.read_tree(["nextflow"])
    problems: list[str] = []
    for path in paths:
        label = f"{path.parent.name}/{path.name}" if path.name == "SKILL.md" else path.name
        extra = EXTRA_ROOTS.get(path.parent.name, [])
        text = code + "\n" + source.read_tree(extra) if extra else code
        body = read_doc(path, source)
        if body.startswith("---\n"):
            body = body.split("\n---\n", 1)[-1]
        for span in dict.fromkeys(SPAN.findall(FENCE.sub("", body))):
            why = check_span(span, source, text, flag_text)
            if why:
                problems.append(f"{label}: `{span}`: {why}")
        if path.name == "architecture.md":
            problems += [f"{label}: {p}" for p in module_map_problems(body, source.names)]
    for line in problems:
        print(line)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
