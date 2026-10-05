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
matches Cargo.toml, and the manifest is a released version no newer than the
package. With ``--tag`` it gates the release workflow before anything publishes:
the tag is a bare ``X.Y.Z`` and Cargo.toml, the lock and the manifest all equal it,
with a dated CHANGELOG section. ``notes`` writes a release's title and body from its
CHANGELOG section (the GitHub Release page).

usage: release.py check [--tag X.Y.Z]
       release.py notes --tag X.Y.Z --title-out FILE --notes-out FILE
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
SEMVER = re.compile(r"^(\d+)\.(\d+)\.(\d+)(?:-[0-9A-Za-z.]+)?$")
HEADING = re.compile(r"^## \[([^\]]+)\](?: - (\d{4}-\d{2}-\d{2}))?(?: — (.+))?\s*$", re.M)
IMAGE = "ghcr.io/msk-access/gbcms"
DERIVED_TAG = IMAGE + ":${workflow.manifest.version}"
DERIVED_BANNER = "gbcms v${workflow.manifest.version}"


def _read(root: Path, rel: str) -> str:
    return (root / rel).read_text(encoding="utf-8")


def _section(text: str, header: str) -> str:
    """A TOML table's body (``[header]`` up to the next table)."""
    m = re.search(rf"^\[{re.escape(header)}\]\s*$(.*?)(?=^\[|\Z)", text, re.M | re.S)
    return m.group(1) if m else ""


def cargo_version(root: Path) -> str | None:
    m = re.search(
        r'^version\s*=\s*"([^"]+)"', _section(_read(root, "rust/Cargo.toml"), "package"), re.M
    )
    return m.group(1) if m else None


def lock_version(root: Path) -> str | None:
    m = re.search(r'name = "gbcms_rs"\nversion = "([^"]+)"', _read(root, "rust/Cargo.lock"))
    return m.group(1) if m else None


def manifest_version(root: Path) -> str | None:
    block = re.search(r"manifest\s*\{(.*?)\n\}", _read(root, "nextflow/nextflow.config"), re.S)
    if not block:
        return None
    m = re.search(r"""^\s*version\s*=\s*['"]([^'"]+)['"]""", block.group(1), re.M)
    return m.group(1) if m else None


def changelog_sections(root: Path) -> list[tuple[str, str | None, str | None, int, int]]:
    """``(version, date, summary, body_start, body_end)`` for each ``## [..]`` heading."""
    text = _read(root, "CHANGELOG.md")
    heads = list(HEADING.finditer(text))
    out = []
    for k, h in enumerate(heads):
        end = heads[k + 1].start() if k + 1 < len(heads) else len(text)
        out.append((h.group(1), h.group(2), h.group(3), h.end(), end))
    return out


def _base(version: str) -> tuple[int, int, int] | None:
    m = SEMVER.match(version)
    return (int(m.group(1)), int(m.group(2)), int(m.group(3))) if m else None


def check(root: Path = ROOT, tag: str | None = None) -> list[str]:
    """Problems with the version sources (empty when they agree)."""
    problems: list[str] = []

    project = _section(_read(root, "pyproject.toml"), "project")
    if not re.search(r'^dynamic\s*=\s*\[[^\]]*"version"', project, re.M):
        problems.append(
            'pyproject.toml: [project] must declare dynamic = ["version"] (from Cargo.toml)'
        )
    if re.search(r"^version\s*=", project, re.M):
        problems.append(
            "pyproject.toml: a static [project] version duplicates Cargo.toml; remove it"
        )
    if re.search(r'__version__\s*=\s*["\']\d', _read(root, "src/gbcms/__init__.py")):
        problems.append(
            "src/gbcms/__init__.py: __version__ is a literal; read the package metadata"
        )

    for path in sorted((root / "nextflow").rglob("*.nf")):
        rel = path.relative_to(root).as_posix()
        text = path.read_text(encoding="utf-8")
        if re.search(re.escape(IMAGE) + r":\d", text):
            problems.append(f"{rel}: literal image tag; use {DERIVED_TAG}")
        for line in text.splitlines():
            if line.strip().startswith("container") and DERIVED_TAG not in line:
                problems.append(f"{rel}: container must read the manifest ({DERIVED_TAG})")
        if re.search(r"gbcms v\d", text):
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

    if tag is None:
        base = _base(cargo) if cargo else None
        if base and manifest and PLAIN.match(manifest) and _base(manifest) > base:  # type: ignore[operator]
            problems.append(
                f"nextflow/nextflow.config: manifest {manifest} is newer than the package {cargo}"
            )
        if not sections:
            problems.append("CHANGELOG.md: no ## [..] section")
        return problems

    if not PLAIN.match(tag):
        problems.append(f"tag {tag!r} must be a bare X.Y.Z (the release workflow matches no other)")
    if cargo != tag:
        problems.append(f"rust/Cargo.toml: version {cargo} is not the tag {tag}")
    if manifest != tag:
        problems.append(
            f"nextflow/nextflow.config: manifest.version {manifest} is not the tag {tag}"
        )
    dated = [s for s in sections if s[0] == tag]
    if not dated or not dated[0][1]:
        problems.append(f"CHANGELOG.md: no dated section ## [{tag}] - YYYY-MM-DD")
    return problems


def notes(root: Path, tag: str) -> tuple[str, str]:
    """The Release page's title and body from the tag's CHANGELOG section."""
    for version, _date, summary, start, end in changelog_sections(root):
        if version == tag:
            body = _read(root, "CHANGELOG.md")[start:end].strip("\n")
            return (f"{tag} — {summary.strip()}" if summary else tag), body.strip() + "\n"
    raise SystemExit(f"CHANGELOG.md has no section for {tag}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser(
        "check", help="check the version sources (with --tag: against a release tag)"
    )
    c.add_argument("--tag")
    n = sub.add_parser("notes", help="write a release's title and body from the CHANGELOG")
    n.add_argument("--tag", required=True)
    n.add_argument("--title-out", required=True)
    n.add_argument("--notes-out", required=True)
    args = ap.parse_args(argv)
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
    problems = check(ROOT, args.tag)
    for p in problems:
        print(f"  ✗ {p}")
    print("version sources agree" if not problems else f"{len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
