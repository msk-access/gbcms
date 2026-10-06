#!/usr/bin/env python3
"""Fail unless every runtime dependency is installed at its declared floor.

CI's floors leg installs `pyproject.toml`'s dependencies with uv's
`--resolution lowest-direct`, which keeps an installed version that already
satisfies a floor and moves to the next release when a floor is yanked; either way
the leg would test something other than the floors and stay green. This compares
each installed version with the `>=` version in `pyproject.toml`.
"""

import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from packaging.requirements import Requirement
from packaging.version import Version

try:  # Python 3.11+
    import tomllib
except ModuleNotFoundError:  # Python 3.10: the `test` group installs tomli there
    import tomli as tomllib  # type: ignore[no-redef,import-not-found]

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    deps = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["dependencies"]
    bad = []
    for req in map(Requirement, deps):
        floors = [s.version for s in req.specifier if s.operator == ">="]
        try:
            got = version(req.name)
        except PackageNotFoundError:
            bad.append(f"{req.name}: not installed")
            continue
        if len(floors) != 1:
            bad.append(f"{req.name}: needs exactly one >= floor, has {req.specifier}")
        elif Version(got) != Version(floors[0]):
            bad.append(f"{req.name}: installed {got}, floor {floors[0]}")
        else:
            print(f"{req.name} {got} (floor)")
    for line in bad:
        print(f"FLOOR NOT INSTALLED: {line}", file=sys.stderr)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
