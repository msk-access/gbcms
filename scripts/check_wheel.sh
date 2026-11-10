#!/usr/bin/env bash
# Install a built gbcms wheel the way a user does, and show that it works.
#
# Usage: scripts/check_wheel.sh WHEEL_DIR [PYTHON] [DIGEST_FILE]
#
#   WHEEL_DIR    a folder holding exactly one gbcms wheel
#   PYTHON       the interpreter to install into (default: python3.11)
#   DIGEST_FILE  where to write the fixture run's output digests (default: none)
#
# Run it from a checkout of the repository: the fixture script and the test suite are
# read from there. The checkout's own gbcms source is not imported; the suite runs
# against the installed wheel. The steps:
#   1. a fresh virtual environment, then `pip install --only-binary=:all:` of the wheel,
#      so a wheel that does not fit this machine fails here instead of pip quietly
#      building the source distribution;
#   2. `pip check`, `gbcms --version`, and a check that gbcms is the installed one;
#   3. a small counting run (scripts/wheel_fixture_run.py), whose output digests are the
#      proof that two installs count the same bytes;
#   4. the test suite (skipped when RUN_SUITE=0).
#
# The suite's hook tests commit in throwaway repositories, so they need git, bash, xargs
# (findutils), awk and gzip on the machine; a minimal container image has to install them.
#
# It creates a temporary folder (printed at the end) and leaves it for you to remove.
# It installs from PyPI, and the suite takes a few minutes: on a cluster, run it inside
# a job, not on the login node.
set -euo pipefail

wheel_dir=${1:?usage: check_wheel.sh WHEEL_DIR [PYTHON] [DIGEST_FILE]}
python=${2:-python3.11}
digest_file=${3:-}
run_suite=${RUN_SUITE:-1}
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

shopt -s nullglob
wheels=("$wheel_dir"/gbcms-*.whl)
if [ "${#wheels[@]}" -ne 1 ]; then
  echo "error: expected one gbcms wheel in $wheel_dir, found ${#wheels[@]}" >&2
  exit 1
fi
wheel=${wheels[0]}

work=$(mktemp -d)
venv=$work/venv

echo "== system"
cat /etc/os-release 2>/dev/null | grep -E '^(PRETTY_NAME|VERSION_ID)=' || uname -sr
getconf GNU_LIBC_VERSION 2>/dev/null || true
echo "== interpreter: $("$python" --version) ($(command -v "$python"))"
echo "== wheel: $(basename "$wheel")"

"$python" -m venv "$venv"
"$venv/bin/python" -m pip install --quiet --upgrade pip
"$venv/bin/python" -m pip install --only-binary=:all: "$wheel"
"$venv/bin/python" -m pip check

echo "== the installed gbcms"
"$venv/bin/python" - <<'PY'
import gbcms
import gbcms._rs as rs

print("gbcms", gbcms.__version__, "from", gbcms.__file__)
print("_rs from", rs.__file__)
for module in (gbcms, rs):
    assert "site-packages" in module.__file__, f"{module.__name__} is not the installed one"
PY
"$venv/bin/gbcms" --version

echo "== a counting run on the synthetic data set"
digests=$("$venv/bin/python" "$repo/scripts/wheel_fixture_run.py" "$work/fixture")
echo "$digests"
if [ -n "$digest_file" ]; then
  mkdir -p "$(dirname "$digest_file")"
  echo "$digests" > "$digest_file"
fi

if [ "$run_suite" = "1" ]; then
  echo "== the test suite, against the installed wheel"
  "$venv/bin/python" -m pip install --quiet --group "$repo/pyproject.toml:test"
  cd "$repo"
  # `pytest`, not `python -m pytest`, and no bytecode: nothing is imported from or
  # written into the checkout's src/. --no-cov: the wheel, not the source, is under test.
  # The venv's bin comes first on PATH: the hook tests make commits whose hooks call a bare
  # `python3`, which must be the venv's (with PyYAML), not the host's.
  PATH="$venv/bin:$PATH" PYTHONDONTWRITEBYTECODE=1 "$venv/bin/pytest" -q -p no:cacheprovider --no-cov tests/
fi

echo "== done; the working folder is $work"
