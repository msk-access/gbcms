#!/usr/bin/env python3
"""Release tooling: the version check and the release notes.

One version source per ecosystem:

* ``rust/Cargo.toml`` (``[package] version``): the Python package (maturin reads it,
  ``dynamic = ["version"]``) and the extension; ``rust/Cargo.lock`` follows
  ``cargo check``.
* ``nextflow/nextflow.config`` (``manifest.version``): the pipeline's container tag
  and banner. It names a released image, so on develop it trails the package.
* ``CHANGELOG.md``: the dated section ``## [X.Y.Z] - YYYY-MM-DD — summary``.

``check`` runs on every PR: the derived locations still read the sources, the lock
matches Cargo.toml, and the manifest names a released image (a dated CHANGELOG
section) older than a dev package; on develop, a released package version means the
post-release bump was missed. With ``--tag`` it gates the release workflow before anything publishes:
the tag is a bare ``X.Y.Z`` and Cargo.toml, the lock and the manifest all equal it,
with a dated CHANGELOG section. ``notes`` writes a release's title and body from its
CHANGELOG section (the GitHub Release page). ``version`` prints the package version,
which labels the deployed docs.

usage: release.py check [--tag X.Y.Z] [--branch NAME]
       release.py notes --tag X.Y.Z --title-out FILE --notes-out FILE
       release.py version
       release.py manifest-version
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Every place a version lives: the three a release edits, the lock that follows, and
# the derived ones the check keeps derived. The release guide's table lists exactly
# these paths (a test compares them).
LOCATIONS: list[tuple[str, str]] = [
    ("rust/Cargo.toml", "edit: `version` under `[package]`, the Python package and the extension"),
    ("rust/Cargo.lock", "follows: run `cargo check` after editing Cargo.toml, commit the lock"),
    ("nextflow/nextflow.config", "edit: `manifest.version`, the pipeline's image tag and banner"),
    ("CHANGELOG.md", "edit: the dated section `## [X.Y.Z] - YYYY-MM-DD — summary`"),
    ("pyproject.toml", 'derived: `dynamic = ["version"]`, maturin reads Cargo.toml'),
    ("src/gbcms/__init__.py", "derived: `__version__` from the installed package metadata"),
    ("nextflow/modules/local/gbcms/*/main.nf", "derived: the container tag reads the manifest"),
    ("nextflow/main.nf", "derived: the banner reads the manifest"),
]

PLAIN = re.compile(r"^\d+\.\d+\.\d+$")
SEMVER = re.compile(r"^(\d+)\.(\d+)\.(\d+)(-[0-9A-Za-z.]+)?$")
HEADING = re.compile(r"^## \[([^\]]+)\](?: - (\d{4}-\d{2}-\d{2}))?(?: — (.+))?\s*$", re.M)
ANY_HEADING = re.compile(r"^## ", re.M)
LINK_REF = re.compile(r"^\[[^\]]+\]:\s*\S+\s*$")
IMAGE = "ghcr.io/msk-access/gbcms"
DERIVED_TAG = IMAGE + ":${workflow.manifest.version}"
DERIVED_BANNER = "gbcms v${workflow.manifest.version}"
CONTAINER = re.compile(r"""^\s*container\s*[("']""")  # the directive, not containerOptions
# GitHub caps a release body at 125,000 characters; longer notes are cut at a section.
NOTES_LIMIT = 120_000

try:  # Python 3.11+
    import tomllib as _toml
except ModuleNotFoundError:  # pragma: no cover - 3.10 dev venvs
    try:
        import tomli as _toml  # type: ignore[no-redef]
    except ModuleNotFoundError:
        _toml = None  # type: ignore[assignment]


def _read(root: Path, rel: str) -> str:
    return (root / rel).read_text(encoding="utf-8")


def _toml_table(root: Path, rel: str, table: str) -> dict | None:
    """A TOML file's table, parsed (None when no TOML parser is available)."""
    if _toml is None:
        return None
    value = _toml.loads(_read(root, rel)).get(table, {})
    return value if isinstance(value, dict) else {}


def _section(text: str, header: str) -> str:
    """A TOML table's raw body (``[header]`` up to the next table): the fallback reader."""
    m = re.search(rf"^\[{re.escape(header)}\]\s*$(.*?)(?=^\[|\Z)", text, re.M | re.S)
    return m.group(1) if m else ""


def _code(text: str) -> list[str]:
    """A Nextflow file's lines that are not comments (``//`` and block comments)."""
    out, block = [], False
    for line in text.splitlines():
        stripped = line.strip()
        if block:
            block = "*/" not in stripped
            continue
        if stripped.startswith("/*"):
            block = "*/" not in stripped
            continue
        if stripped.startswith(("//", "*")):
            continue
        out.append(line)
    return out


def cargo_version(root: Path) -> str | None:
    package = _toml_table(root, "rust/Cargo.toml", "package")
    if package is not None:
        v = package.get("version")
        return v if isinstance(v, str) else None
    m = re.search(
        r"""^version\s*=\s*["']([^"']+)["']""",
        _section(_read(root, "rust/Cargo.toml"), "package"),
        re.M,
    )
    return m.group(1) if m else None


def lock_version(root: Path) -> str | None:
    """The ``gbcms_rs`` entry's version in Cargo.lock (any field order)."""
    for entry in _read(root, "rust/Cargo.lock").split("[[package]]"):
        if re.search(r'^name = "gbcms_rs"$', entry, re.M):
            m = re.search(r'^version = "([^"]+)"$', entry, re.M)
            return m.group(1) if m else None
    return None


def manifest_version(root: Path) -> str | None:
    code = "\n".join(_code(_read(root, "nextflow/nextflow.config")))
    block = re.search(r"^\s*manifest\s*\{(.*?)^\s*\}", code, re.S | re.M)
    if not block:
        return None
    m = re.search(r"""^\s*version\s*=\s*['"]([^'"]+)['"]""", block.group(1), re.M)
    return m.group(1) if m else None


def changelog_sections(root: Path) -> list[tuple[str, str | None, str | None, int, int]]:
    """``(version, date, summary, body_start, body_end)`` for each ``## [..]`` heading
    that parses; a body ends at the next ``## `` heading of any form."""
    text = _read(root, "CHANGELOG.md")
    starts = [m.start() for m in ANY_HEADING.finditer(text)] + [len(text)]
    out = []
    for h in HEADING.finditer(text):
        end = next(s for s in starts if s > h.start())
        out.append((h.group(1), h.group(2), h.group(3), h.end(), end))
    return out


def _base(version: str) -> tuple[int, int, int] | None:
    m = SEMVER.match(version)
    return (int(m.group(1)), int(m.group(2)), int(m.group(3))) if m else None


def _prerelease(version: str) -> bool:
    m = SEMVER.match(version)
    return bool(m and m.group(4))


def check(root: Path = ROOT, tag: str | None = None, branch: str | None = None) -> list[str]:
    """Problems with the version sources (empty when they agree). ``tag``: check
    against a release tag. ``branch``: the branch being checked (develop must move to
    the next dev version once a release is out)."""
    problems: list[str] = []

    for loc, _ in LOCATIONS:
        if not list(root.glob(loc)):
            problems.append(f"{loc}: listed as a version location but missing")

    project = _toml_table(root, "pyproject.toml", "project")
    if project is None:  # no TOML parser: read the table's text
        raw = _section(_read(root, "pyproject.toml"), "project")
        dynamic_version = bool(re.search(r'^dynamic\s*=\s*\[[^\]]*"version"', raw, re.M))
        static_version = bool(re.search(r"^version\s*=", raw, re.M))
    else:
        dynamic_version = "version" in project.get("dynamic", [])
        static_version = "version" in project
    if not dynamic_version:
        problems.append(
            'pyproject.toml: [project] must declare dynamic = ["version"] (from Cargo.toml)'
        )
    if static_version:
        problems.append(
            "pyproject.toml: a static [project] version duplicates Cargo.toml; remove it"
        )
    if re.search(r'__version__\s*=\s*["\']\d', _read(root, "src/gbcms/__init__.py")):
        problems.append(
            "src/gbcms/__init__.py: __version__ is a literal; read the package metadata"
        )

    for path in sorted((root / "nextflow").rglob("*.nf")):
        rel = path.relative_to(root).as_posix()
        code = _code(path.read_text(encoding="utf-8"))
        joined = "\n".join(code)
        if re.search(re.escape(IMAGE) + r":\d", joined):
            problems.append(f"{rel}: literal image tag; use {DERIVED_TAG}")
        for line in code:
            if CONTAINER.match(line) and DERIVED_TAG not in line:
                problems.append(f"{rel}: container must read the manifest ({DERIVED_TAG})")
        if re.search(r"gbcms v\d", joined):
            problems.append(f"{rel}: literal banner version; use {DERIVED_BANNER}")
    if DERIVED_BANNER not in _read(root, "nextflow/main.nf"):
        problems.append(f"nextflow/main.nf: the banner must read the manifest ({DERIVED_BANNER})")

    cargo, lock, manifest = cargo_version(root), lock_version(root), manifest_version(root)
    if cargo is None or _base(cargo) is None:
        problems.append(f"rust/Cargo.toml: no semver [package] version (found {cargo!r})")
    if lock != cargo:
        problems.append(f"rust/Cargo.lock: gbcms_rs is {lock}, Cargo.toml {cargo}; run cargo check")
    if manifest is None or not PLAIN.match(manifest):
        problems.append(
            f"nextflow/nextflow.config: manifest.version {manifest!r} is not a released X.Y.Z"
        )
    sections = changelog_sections(root)
    released = {v for v, date, *_ in sections if date}

    if tag is None:
        if not sections:
            problems.append("CHANGELOG.md: no ## [..] section")
        if manifest and PLAIN.match(manifest):
            # The manifest names an image, which exists only once its release is out.
            if manifest not in released:
                problems.append(
                    f"nextflow/nextflow.config: manifest {manifest} has no dated CHANGELOG "
                    "section, so its image is not published"
                )
            base = _base(cargo) if cargo else None
            mbase = _base(manifest)
            if (
                base
                and mbase
                and cargo
                and (mbase > base or (_prerelease(cargo) and mbase >= base))
            ):
                problems.append(
                    f"nextflow/nextflow.config: manifest {manifest} is not older than the "
                    f"package {cargo}; until {manifest} is released, name the last released image"
                )
        if branch == "develop" and cargo and PLAIN.match(cargo) and cargo in released:
            problems.append(
                f"rust/Cargo.toml: develop is still at the released {cargo}; move it to the next "
                "dev version (e.g. X.Y+1.0-dev.0) so dev builds don't stamp a released version"
            )
        return problems

    if not PLAIN.match(tag):
        problems.append(f"tag {tag!r} must be a bare X.Y.Z (the release workflow matches no other)")
    if cargo != tag:
        problems.append(f"rust/Cargo.toml: version {cargo} is not the tag {tag}")
    if manifest != tag:
        problems.append(
            f"nextflow/nextflow.config: manifest.version {manifest} is not the tag {tag}"
        )
    if tag not in released:
        problems.append(f"CHANGELOG.md: no dated section ## [{tag}] - YYYY-MM-DD")
    return problems


def notes(root: Path, tag: str) -> tuple[str, str]:
    """The Release page's title and body from the tag's CHANGELOG section: the body
    ends at the next ``## `` heading, without trailing link references, and is cut at
    a ``### `` section under GitHub's size limit with a link to the full CHANGELOG."""
    for version, _date, summary, start, end in changelog_sections(root):
        if version != tag:
            continue
        lines = _read(root, "CHANGELOG.md")[start:end].strip("\n").splitlines()
        while lines and (not lines[-1].strip() or LINK_REF.match(lines[-1])):
            lines.pop()
        body = "\n".join(lines).strip() + "\n"
        if len(body) > NOTES_LIMIT:
            more = (
                f"\n_The notes continue in [CHANGELOG.md](https://github.com/msk-access/gbcms/"
                f"blob/{tag}/CHANGELOG.md)._\n"
            )
            cut = body.rfind("\n### ", 0, NOTES_LIMIT - len(more))
            cut = cut if cut > 0 else body.rfind("\n", 0, NOTES_LIMIT - len(more))
            body = body[:cut].rstrip() + "\n" + more
        return (f"{tag} — {summary.strip()}" if summary else tag), body
    raise SystemExit(f"CHANGELOG.md has no section for {tag}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser(
        "check", help="check the version sources (with --tag: against a release tag)"
    )
    c.add_argument("--tag")
    c.add_argument("--branch", help="the branch being checked (develop: no released version)")
    sub.add_parser("version", help="print the package version (rust/Cargo.toml)")
    sub.add_parser("manifest-version", help="print the Nextflow manifest version")
    n = sub.add_parser("notes", help="write a release's title and body from the CHANGELOG")
    n.add_argument("--tag", required=True)
    n.add_argument("--title-out", required=True)
    n.add_argument("--notes-out", required=True)
    args = ap.parse_args(argv)
    if args.cmd == "version":
        print(cargo_version(ROOT) or "")
        return 0
    if args.cmd == "manifest-version":
        print(manifest_version(ROOT) or "")
        return 0
    if args.cmd == "notes":
        title, body = notes(ROOT, args.tag)
        Path(args.title_out).write_text(title + "\n", encoding="utf-8")
        Path(args.notes_out).write_text(body, encoding="utf-8")
        print(f"{title}\n({len(body.splitlines())} lines of notes)")
        return 0
    print(
        f"Cargo.toml {cargo_version(ROOT)} | Cargo.lock {lock_version(ROOT)} | "
        f"Nextflow manifest {manifest_version(ROOT)}" + (f" | tag {args.tag}" if args.tag else "")
    )
    problems = check(ROOT, args.tag, args.branch)
    for p in problems:
        print(f"  ✗ {p}")
    print("version sources agree" if not problems else f"{len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
