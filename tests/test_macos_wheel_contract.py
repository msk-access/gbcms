"""A macOS wheel for Apple silicon is built, published and installed in CI.

Before, `pip install gbcms` on a Mac built the source distribution (Rust, clang, htslib).
Contracts:

- the release builds the arm64 wheel after the version check, with the one action the
  wheel checks also build with, and publishes it with the Linux wheel;
- the wheel checks install that wheel on an Apple silicon runner at both ends of the
  Python range, run the suite there, and compare the counting run's output digests with the
  Linux ones.
"""

from pathlib import Path

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


def test_the_wheel_checks_install_the_macos_wheel_at_both_ends_of_the_range():
    jobs = _jobs("wheels.yml")
    assert jobs["macos-build"]["runs-on"] == "macos-latest"
    assert "./.github/actions/macos-wheel" in [s.get("uses") for s in jobs["macos-build"]["steps"]]
    test = jobs["macos-test"]
    assert test["runs-on"] == "macos-latest" and "macos-build" in _needs(jobs, "macos-test")
    assert {"3.10", "3.14"} <= set(test["strategy"]["matrix"]["python"])
    assert "scripts/check_wheel.sh" in "\n".join(str(s.get("run", "")) for s in test["steps"])
    assert "macos-test" in _needs(jobs, "compare")
