"""Group 7 (D4 #139, PR A): Python dependencies, CI coverage and the Docker lock.

Measured (CYCLE_6.6.0_PLAN.md, "D4 PR A"): every Python from 3.10 to 3.14 passes on
the latest releases, but CI tested only 3.11 and 3.12; the declared floors crashed
(typer 0.9.0) or could not install (pysam 0.21.0 on arm64), the measured minimums
being typer 0.15.4 and pysam 0.22.0; two dev-dependency lists had drifted apart (and
the PEP 735 one makes `maturin develop` need pip >= 25.1); the image resolved its
dependencies afresh at every build and CI never ran it.

Contracts:
- the floors are the measured minimums, and a PR CI leg installs them;
- PR CI covers the ends (3.10 at the floors, 3.14), the image's 3.11 and macOS;
- one dev list: PEP 735 groups, `test` inside `dev`;
- a weekly job runs the latest releases on 3.10-3.14 and opens one issue on failure;
- the image installs a hash-pinned lock, checks it, and CI runs the image;
- the release guide refreshes the lock; the developer docs say pip >= 25.1.
"""

import re
from pathlib import Path

import pytest
import yaml
from packaging.requirements import Requirement
from packaging.version import Version

pytestmark = pytest.mark.xfail(strict=True, reason="D4 PR A is not implemented yet")

ROOT = Path(__file__).resolve().parents[1]
WF = ROOT / ".github" / "workflows"
LOCK = ROOT / "docker" / "requirements.lock"

FLOORS = {
    "pysam": "0.22.0",
    "typer": "0.15.4",
    "click": "8.0",
    "rich": "13.0.0",
    "pydantic": "2.0.0",
    "polars": "1.0.0",
}


def _toml():
    try:
        import tomllib
    except ModuleNotFoundError:  # Python 3.10
        import tomli as tomllib
    return tomllib.loads((ROOT / "pyproject.toml").read_text())


def _workflow(name):
    data = yaml.safe_load((WF / name).read_text())
    # PyYAML reads the bare key `on` as True.
    data["on"] = data.pop(True, data.get("on"))
    return data


def _runs(job):
    return "\n".join(step.get("run", "") for step in job["steps"])


# ── pyproject ────────────────────────────────────────────────────────────────


def test_the_floors_are_the_measured_minimums():
    reqs = {Requirement(r).name: Requirement(r) for r in _toml()["project"]["dependencies"]}
    assert set(reqs) == set(FLOORS)
    for name, floor in FLOORS.items():
        lows = [s.version for s in reqs[name].specifier if s.operator == ">="]
        assert lows == [floor], (name, str(reqs[name].specifier))


def test_the_classifiers_list_python_310_to_314():
    toml = _toml()["project"]
    assert toml["requires-python"] == ">=3.10"
    listed = {
        c.rsplit(":: ", 1)[1]
        for c in toml["classifiers"]
        if c.startswith("Programming Language :: Python :: 3.")
    }
    assert listed == {"3.10", "3.11", "3.12", "3.13", "3.14"}


def test_one_dev_dependency_list():
    toml = _toml()
    assert "dev" not in toml["project"].get("optional-dependencies", {})
    groups = toml["dependency-groups"]
    names = {
        g: {Requirement(r).name for r in reqs if isinstance(r, str)} for g, reqs in groups.items()
    }
    assert {"pytest", "pytest-cov", "pytest-mock", "pyarrow", "pyyaml"} <= names["test"]
    assert {"include-group": "test"} in groups["dev"]
    assert {"black", "ruff", "mypy", "types-pyyaml", "mkdocs", "mkdocs-material"} <= names["dev"]
    assert "pytest-benchmark" not in names["test"] | names["dev"]


# ── CI ───────────────────────────────────────────────────────────────────────


def test_pr_ci_covers_the_ends_the_image_and_the_floors():
    job = _workflow("test.yml")["jobs"]["test"]
    legs = {
        (leg["os"], str(leg["python-version"]), bool(leg.get("floors")))
        for leg in job["strategy"]["matrix"]["include"]
    }
    assert legs == {
        ("ubuntu-latest", "3.10", True),
        ("ubuntu-latest", "3.11", False),
        ("ubuntu-latest", "3.14", False),
        ("macos-latest", "3.12", False),
    }
    runs = _runs(job)
    assert "--resolution lowest-direct" in runs and "matrix.floors" in yaml.safe_dump(job)
    assert "--group test" in runs
    assert "matrix.os == 'Linux'" not in yaml.safe_dump(job), "a step that can never run"


def test_no_workflow_installs_the_old_dev_extra():
    for wf in WF.glob("*.yml"):
        assert ".[dev]" not in wf.read_text(), wf.name
    lint = _runs(_workflow("test.yml")["jobs"]["lint"])
    assert "--group dev" in lint
    assert "scipy-stubs" not in lint


def test_a_weekly_job_runs_the_latest_releases_and_reports_failure():
    wf = _workflow("latest-deps.yml")
    assert "schedule" in wf["on"] and "workflow_dispatch" in wf["on"]
    assert wf.get("permissions") == {"contents": "read"}
    jobs = wf["jobs"]
    pythons = set()
    for job in jobs.values():
        matrix = job.get("strategy", {}).get("matrix", {})
        pythons |= {str(p) for p in matrix.get("python-version", [])}
    assert pythons >= {"3.10", "3.11", "3.12", "3.13", "3.14"}
    text = (WF / "latest-deps.yml").read_text()
    assert "cargo update" in text and "--upgrade" in text
    reporters = [j for j in jobs.values() if j.get("permissions", {}).get("issues") == "write"]
    assert len(reporters) == 1
    rep = reporters[0]
    assert "failure()" in str(rep.get("if", "")) and "gh issue" in _runs(rep)


# ── Docker ───────────────────────────────────────────────────────────────────


def test_the_image_installs_a_hash_pinned_lock():
    text = LOCK.read_text()
    entries = re.findall(r"^([A-Za-z0-9._-]+)==([^\s\\]+)", text, re.M)
    assert entries, "the lock pins packages"
    blocks = re.split(r"\n(?=[A-Za-z0-9._-]+==)", text.split("\n", 1)[1] if text else "")
    for block in blocks:
        if re.match(r"[A-Za-z0-9._-]+==", block):
            assert "--hash=sha256:" in block, block.splitlines()[0]
    pinned = {n.lower().replace("_", "-"): Version(v) for n, v in entries}
    for r in _toml()["project"]["dependencies"]:
        req = Requirement(r)
        assert req.name in pinned, f"{req.name} is not in the lock"
        assert pinned[req.name] in req.specifier, (req.name, pinned[req.name], str(req.specifier))
    docker = (ROOT / "Dockerfile").read_text()
    assert "docker/requirements.lock" in docker
    assert re.search(r"pip install[^\n]*--require-hashes[^\n]*--no-deps[^\n]*-r", docker)
    assert re.search(r"pip install[^\n]*--no-deps[^\n]*\.whl", docker)
    assert "pip check" in docker


def test_ci_runs_the_image():
    job = _workflow("test.yml")["jobs"]["docker-test"]
    runs = _runs(job)
    assert re.search(r"docker run[^\n]*gbcms:test[^\n]*--version", runs)
    assert "pip check" in runs
    load = [s for s in job["steps"] if "docker/build-push-action" in s.get("uses", "")]
    assert load and load[0]["with"].get("load") is True


# ── Docs ─────────────────────────────────────────────────────────────────────


def test_the_release_guide_refreshes_the_lock_and_lists_the_weekly_job():
    guide = (ROOT / "docs" / "development" / "release-guide.md").read_text()
    assert re.search(r"uv pip compile[^\n]*--generate-hashes", guide)
    assert "docker/requirements.lock" in guide
    assert "latest-deps.yml" in guide


def test_the_developer_docs_say_pip_251_before_maturin_develop():
    for doc in ("CONTRIBUTING.md", "docs/development/developer-guide.md"):
        text = (ROOT / doc).read_text()
        first = text.index("maturin develop")
        assert re.search(r"pip>=25\.1", text[:first]), doc
