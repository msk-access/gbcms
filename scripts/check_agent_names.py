#!/usr/bin/env python3
"""Fail when an agent skill or rule names code, a flag or a file that does not exist.

The skills (`.agents/skills/*/SKILL.md`) and the rules (`.agents/rules/*.md`) are read
on trust, so a renamed function or a deleted file stays in them until a reader acts on
it. This checks every backticked name:

* a path (`rust/src/shared/fragment.rs`, `engine.rs`, `tests/test_x.py::TestY`) must
  exist in the tree; a bare file name matches a tracked file by its name, a partial path
  by its tail;
* a command-line flag (`--umi-tag`, `--no-baq`) must appear in `src/` or the Nextflow
  files;
* an identifier (`check_allele_with_qual`, `window::allele_kind`, `BIN_WINDOW`, a CamelCase
  type, a `snake_case` field) must appear as a whole word in `src/` or `rust/src/`
  (the Nextflow and read-filters skills also search `nextflow/`);
* `NAME = value` for an upper-case constant must also carry that value in the source;
* the module map in `architecture.md` must list every Rust and Python module file, and
  only files that exist.

Plain words (`reverse`, `sw`), shell commands and phrases are not identifiers and are not
checked. `EXTERNAL` lists names that belong to another tool, each with its reason.
Fenced code blocks are skipped, except the module map.

Usage: python3 scripts/check_agent_names.py [FILE ...]   (default: every skill and rule)
Exits 1 and prints one line per unresolved name.
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
    # Git's, gh's and pip's own names.
    "HEAD",
    "--no-ff",
    "--add-sub-issue",
    "--body-file",
    "--require-hashes",
    "--no-deps",
    # A gitignored directory that `git worktree add` creates.
    ".claude/worktrees/",
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


def tracked() -> list[str]:
    """Every file git tracks, plus untracked files that are not ignored.

    Without git (a source archive) the tree is walked instead, skipping build output.
    """
    try:
        out = subprocess.run(
            ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        return out.splitlines()
    except (OSError, subprocess.CalledProcessError):
        skip = {".git", "target", ".venv", "venv", "__pycache__", "node_modules"}
        return [
            str(p.relative_to(ROOT))
            for p in ROOT.rglob("*")
            if p.is_file() and not skip.intersection(p.relative_to(ROOT).parts)
        ]


def read_tree(roots: list[str], files: list[str]) -> str:
    """The text of every .rs, .py, .pyi, .nf, .config and .toml file under the roots."""
    chunks = []
    for name in files:
        if not name.endswith((".rs", ".py", ".pyi", ".nf", ".config", ".toml")):
            continue
        if any(name == r or name.startswith(r.rstrip("/") + "/") for r in roots):
            try:
                chunks.append((ROOT / name).read_text(errors="replace"))
            except OSError:
                continue
    return "\n".join(chunks)


def is_identifier_like(word: str) -> bool:
    """A name worth looking up: it has an underscore, a CamelCase hump or is a constant."""
    if word.startswith("_") and len(word) > 1:
        return True
    if "_" in word.strip("_"):
        return True
    if re.search(r"[a-z][A-Z]", word):
        return True
    return word.isupper() and len(word) >= 4


def path_exists(span: str, names: list[str]) -> bool:
    """Whether a path-like span names a tracked file or directory."""
    path = span.split("::")[0]
    if "*" in path or "<" in path or "$" in path or "{" in path:
        pattern = re.escape(path).replace(r"\*", ".*")
        pattern = re.sub(r"\\<.*?\\>", ".*", pattern)
        pattern = re.sub(r"\\\$\w+", ".*", pattern)
        return any(re.search(pattern, n) for n in names)
    if path.endswith("/"):
        return any(n.startswith(path) or ("/" + path) in n for n in names)
    return any(n == path or n.endswith("/" + path) for n in names) or (ROOT / path).exists()


def check_span(span: str, skill: str, names: list[str], text: str, flag_text: str) -> str | None:
    """Return why the span's name does not resolve, or None when it does or is not a name."""
    span = span.strip()
    if span.startswith("--"):
        flag = re.match(r"--[A-Za-z0-9][A-Za-z0-9-]*", span)
        if flag is None or flag.group(0) in EXTERNAL or span.endswith("-"):
            return None
        if "<" in span or " false" in span:
            return None
        return None if flag.group(0) in flag_text else f"flag {flag.group(0)} is not in the code"
    # A file or directory: a known top-level directory, a file suffix, or a bare file name.
    head = span.split("::")[0].split(" ")[0]
    # Build products (the compiled extension) are not tracked; shared-library names are not paths.
    if head in EXTERNAL or head.startswith("~") or head.endswith(".so") or "..." in head:
        return None
    if (
        (head.endswith(SUFFIXES) or head.startswith(TOP_LEVEL))
        and re.fullmatch(r"[\w./*<>${}-]+", head)
        and not head.startswith("-")
        and (not head.startswith(".") or head.startswith(TOP_LEVEL))
    ):
        if not path_exists(span, names):
            return f"path {head} does not exist"
        if "::" in span and head.endswith(".py"):
            tail = span.split("::", 1)[1].split("::")[0].split("[")[0]
            if not re.search(rf"\b{re.escape(tail)}\b", (ROOT / head).read_text()):
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


def module_map_problems(path: Path, names: list[str]) -> list[str]:
    """Compare the module-map tree in architecture.md with the Rust and Python files."""
    text = path.read_text()
    block = re.search(r"## Project Layout\s+```.*?\n(.*?)```", text, re.S)
    if block is None:
        return [f"{path.name}: no module map under '## Project Layout'"]
    listed = set()
    for line in block.group(1).splitlines():
        if not re.match(r"[│ ├└─]*[\w.]+(?:\.rs|\.py|\.pyi)\b", line):
            continue  # a directory line, or a comment line
        listed.update(re.findall(r"[\w]+\.(?:rs|py|pyi)\b", line.split("#")[0]))
    listed -= MAP_IGNORED
    actual = {
        n.rsplit("/", 1)[-1]
        for n in names
        if (n.startswith("rust/src/") and n.endswith(".rs"))
        or (n.startswith("src/gbcms/") and n.endswith((".py", ".pyi")))
    } - MAP_IGNORED
    problems = [f"module map lists {f}, which is not a file" for f in sorted(listed - actual)]
    problems += [f"module map omits {f}" for f in sorted(actual - listed)]
    return problems


def main(argv: list[str]) -> int:
    names = tracked()
    paths = (
        [Path(a).resolve() for a in argv]
        if argv
        else sorted(ROOT.glob(".agents/skills/*/SKILL.md"))
        + sorted(ROOT.glob(".agents/rules/*.md"))
    )
    code = read_tree(["src", "rust/src"], names)
    flag_text = code + "\n" + read_tree(["nextflow"], names)
    problems: list[str] = []
    for path in paths:
        label = f"{path.parent.name}/{path.name}" if path.name == "SKILL.md" else path.name
        extra = EXTRA_ROOTS.get(path.parent.name, [])
        text = code + "\n" + read_tree(extra, names) if extra else code
        body = path.read_text()
        if body.startswith("---\n"):
            body = body.split("\n---\n", 1)[-1]
        for span in dict.fromkeys(SPAN.findall(FENCE.sub("", body))):
            why = check_span(span, path.parent.name, names, text, flag_text)
            if why:
                problems.append(f"{label}: `{span}`: {why}")
        if path.name == "architecture.md":
            problems += [f"{label}: {p}" for p in module_map_problems(path, names)]
    for line in problems:
        print(line)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
