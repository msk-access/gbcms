"""A macOS wheel for Apple silicon is built, published and installed in CI.

Before, `pip install gbcms` on a Mac built the source distribution (Rust, clang, htslib).
Contracts:

- the release builds the arm64 wheel after the version check, with the one action the
  wheel checks also build with, and publishes it with the Linux wheel;
- the wheel checks install that wheel on an Apple silicon runner at both ends of the
  Python range, run the suite there, and compare the counting run's output digests with the
  Linux ones.
"""

import os
import subprocess
import zipfile
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
WF = ROOT / ".github" / "workflows"
ACTION = ROOT / ".github" / "actions" / "macos-wheel" / "action.yml"


def _jobs(name):
    return yaml.safe_load((WF / name).read_text())["jobs"]


def _needs(jobs, job):
    n = jobs[job].get("needs", [])
    return {n} if isinstance(n, str) else set(n)


def _ancestors(jobs, job):
    seen, todo = set(), list(_needs(jobs, job))
    while todo:
        j = todo.pop()
        if j not in seen:
            seen.add(j)
            todo.extend(_needs(jobs, j))
    return seen


def test_the_release_builds_and_publishes_the_apple_silicon_wheel():
    jobs = _jobs("release.yml")
    macos = jobs["macos"]
    # macos-latest is arm64 (macos-15 and later); macos-13 and any -intel runner are not.
    # macos-14 images retire on 2026-11-02.
    assert macos["runs-on"] == "macos-latest", macos["runs-on"]
    assert "verify" in _ancestors(jobs, "macos"), "nothing builds before the version check"
    assert "macos" in _needs(jobs, "release")
    uses = [s.get("uses") for s in macos["steps"]]
    assert "./.github/actions/macos-wheel" in uses, uses
    upload = next(s for s in macos["steps"] if "upload-artifact" in str(s.get("uses", "")))
    # the publish and release-page jobs download every `wheels-*` artifact
    assert upload["with"]["name"].startswith("wheels-"), upload["with"]


def test_the_action_targets_macos_11_and_records_the_commit():
    action = yaml.safe_load(ACTION.read_text())
    (step,) = [s for s in action["runs"]["steps"] if "maturin-action" in s.get("uses", "")]
    assert step["env"]["MACOSX_DEPLOYMENT_TARGET"] == "11.0"
    assert step["env"]["GBCMS_BUILD_COMMIT"] == "${{ inputs.commit }}"
    assert action["inputs"]["interpreter"]["default"] == "python3.10"


def test_the_action_checks_the_wheel_needs_no_newer_macos_than_its_target():
    # the tests run on macos-latest and would install a wheel that needs a newer macOS
    action = yaml.safe_load(ACTION.read_text())
    steps = action["runs"]["steps"]
    build = next(i for i, s in enumerate(steps) if "maturin-action" in s.get("uses", ""))
    (check,) = [s for s in steps[build + 1 :] if s.get("name") == "Check the wheel's macOS floor"]
    assert check["env"]["FLOOR"] == steps[build]["env"]["MACOSX_DEPLOYMENT_TARGET"]
    assert check["env"]["OUT"] == "${{ inputs.out }}"
    # the tag and every binary's load commands
    assert "tag <= floor" in check["run"] and "otool" in check["run"], check["run"]


def _load_commands(minos):
    # what `otool -l` prints for a Mach-O binary's minimum macOS (None: no such command)
    if minos is None:
        return "Load command 0\n      cmd LC_SEGMENT_64\n  cmdsize 72\n"
    return f"Load command 9\n      cmd LC_BUILD_VERSION\n  cmdsize 32\n platform 1\n    minos {minos}\n      sdk 15.0\n"


@pytest.mark.parametrize(
    "platform, minos, passes, why",
    [
        ("macosx_11_0_arm64", "11.0", True, None),
        ("macosx_14_0_arm64", "11.0", False, "macosx_14_0_arm64"),  # the tag needs a newer macOS
        ("macosx_11_0_arm64", "14.0", False, "(14, 0)"),  # the binary needs a newer macOS
        ("macosx_11_0_arm64", None, False, "no minimum macOS"),  # the binary names no minimum
    ],
)
def test_the_floor_check_rejects_a_wheel_that_needs_a_newer_macos(
    tmp_path, platform, minos, passes, why
):
    # the action's own script, on a made-up wheel, with a stand-in for otool (so it runs off macOS too)
    action = yaml.safe_load(ACTION.read_text())
    (check,) = [
        s for s in action["runs"]["steps"] if s.get("name") == "Check the wheel's macOS floor"
    ]
    dist = tmp_path / "dist"
    dist.mkdir()
    with zipfile.ZipFile(dist / f"gbcms-0.0.0-cp310-abi3-{platform}.whl", "w") as whl:
        whl.writestr("gbcms/_rs.abi3.so", b"made up")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    otool = bin_dir / "otool"
    otool.write_text("#!/bin/sh\ncat <<'EOF'\n" + _load_commands(minos) + "EOF\n")
    otool.chmod(0o755)
    env = {
        **os.environ,
        "OUT": str(dist),
        "FLOOR": check["env"]["FLOOR"],
        "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
    }
    run = subprocess.run(["bash", "-c", check["run"]], env=env, capture_output=True, text=True)
    assert (run.returncode == 0) == passes, (run.returncode, run.stdout, run.stderr)
    if why:
        assert why in run.stderr, run.stderr


def test_the_wheel_checks_install_the_macos_wheel_at_both_ends_of_the_range():
    jobs = _jobs("wheels.yml")
    assert jobs["macos-build"]["runs-on"] == "macos-latest"
    assert "./.github/actions/macos-wheel" in [s.get("uses") for s in jobs["macos-build"]["steps"]]
    test = jobs["macos-test"]
    assert test["runs-on"] == "macos-latest" and "macos-build" in _needs(jobs, "macos-test")
    assert {"3.10", "3.14"} <= set(test["strategy"]["matrix"]["python"])
    assert "scripts/check_wheel.sh" in "\n".join(str(s.get("run", "")) for s in test["steps"])
    assert "macos-test" in _needs(jobs, "compare")
