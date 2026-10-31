"""The Linux wheel installs on RHEL 8 (glibc 2.28), and CI proves it there.

Measured: the bundled curl (8.x) refuses OpenSSL older than 3.0, and RHEL 8 ships 1.1.1, so
the wheel was built on manylinux_2_34 (glibc 2.34, RHEL 9) and RHEL 8 built the source
distribution. gbcms opens only paths that exist, so htslib's remote-file support (the only
user of libcurl and OpenSSL) is not needed. Contracts:

- the extension links no libcurl and no OpenSSL (the lock holds neither crate);
- the release builds its Linux wheel in the manylinux_2_28 image, through the one action
  the wheel checks also build with;
- the wheel checks install that wheel in Rocky Linux 8 and 9 with pip allowed binary
  wheels only, run the suite there, and compare a counting run's output digests.
"""

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
WF = ROOT / ".github" / "workflows"
ACTION = ROOT / ".github" / "actions" / "linux-wheel" / "action.yml"


def _toml(path):
    try:
        import tomllib
    except ModuleNotFoundError:  # Python 3.10
        import tomli as tomllib
    return tomllib.loads(path.read_text())


def _jobs(name):
    return yaml.safe_load((WF / name).read_text())["jobs"]


def test_the_extension_links_no_libcurl_and_no_openssl():
    lock = {p["name"] for p in _toml(ROOT / "rust" / "Cargo.lock")["package"]}
    assert not lock & {"curl-sys", "openssl-sys", "openssl-src", "openssl"}, lock
    htslib = _toml(ROOT / "rust" / "Cargo.toml")["dependencies"]["rust-htslib"]
    assert htslib["default-features"] is False, htslib
    assert "curl" not in htslib["features"], htslib


def test_the_release_builds_the_linux_wheel_for_rhel_8():
    steps = _jobs("release.yml")["linux"]["steps"]
    build = next(s for s in steps if s.get("uses") == "./.github/actions/linux-wheel")
    assert build["with"]["manylinux"] == "manylinux_2_28", build
    assert ACTION.is_file()


def test_the_action_builds_in_the_requested_manylinux_image_without_openssl():
    action = yaml.safe_load(ACTION.read_text())
    (step,) = [s for s in action["runs"]["steps"] if "maturin-action" in s.get("uses", "")]
    assert step["with"]["manylinux"] == "${{ inputs.manylinux }}"
    assert "OPENSSL" not in step["with"]["docker-options"]
    assert "openssl" not in step["with"]["before-script-linux"].lower()
    assert "curl" not in step["with"]["before-script-linux"].lower()


def test_the_wheel_checks_install_the_wheel_in_rocky_linux_8_and_9():
    wf = yaml.safe_load((WF / "wheels.yml").read_text())
    jobs = wf["jobs"]
    legs = jobs["rocky"]["strategy"]["matrix"]["include"]
    built = {(leg["rocky"], leg["manylinux"]) for leg in legs}
    assert {("8", "manylinux_2_28"), ("9", "manylinux_2_28")} <= built, built
    assert "build" in jobs["rocky"]["needs"] and "rocky" in jobs["compare"]["needs"]
    script = (ROOT / "scripts" / "check_wheel.sh").read_text()
    # a wheel that does not fit must fail, not fall back to the source distribution
    assert "--only-binary=:all:" in script
    assert "wheel_fixture_run.py" in script and "pytest" in script
    # the checks run when the build recipe changes
    paths = wf["on"] if "on" in wf else wf[True]
    for event in ("pull_request", "push"):
        assert ".github/actions/linux-wheel/**" in paths[event]["paths"], event
        assert "rust/**" in paths[event]["paths"], event


def test_the_suite_in_a_minimal_image_has_what_the_hook_tests_run():
    """The suite's hook tests commit in throwaway repositories: the pre-commit hook calls
    xargs, awk, gzip and a bare python3. A minimal Rocky image lacks findutils, so the dnf
    line installs them, and the venv's bin leads PATH when pytest runs."""
    runs = "\n".join(str(step.get("run", "")) for step in _jobs("wheels.yml")["rocky"]["steps"])
    (dnf,) = [ln for ln in runs.splitlines() if ln.strip().startswith("dnf install")]
    for package in ("findutils", "gawk", "gzip", "git-core"):
        assert package in dnf.split(), (package, dnf)
    script = (ROOT / "scripts" / "check_wheel.sh").read_text()
    pytest_line = next(ln for ln in script.splitlines() if "bin/pytest" in ln and "-q" in ln)
    assert pytest_line.strip().startswith('PATH="$venv/bin:$PATH"'), pytest_line


def test_the_container_script_is_one_single_quoted_word():
    """The Rocky steps run `bash -c '...'`: one apostrophe in a comment inside it ends the
    script early and the rest runs on the runner."""
    steps = _jobs("wheels.yml")["rocky"]["steps"]
    run = next(s["run"] for s in steps if "bash -euo pipefail -c" in str(s.get("run", "")))
    assert run.count("'") == 2, run.count("'")
