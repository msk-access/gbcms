"""Group 7 release-infrastructure contracts (D1 #136, D2 #137, P3 #152).

Operator decisions 2026-10-05: one version source per ecosystem (Cargo for the
Python package and the extension, the Nextflow manifest for the pipeline) with a
check on every PR that also gates the release workflow before anything publishes;
the GitHub Release page created by the workflow from the CHANGELOG section, with
the built artifacts, their checksums and build-provenance attestations; the
architecture page states that the bin window is a floor, not a maximum.
"""

import importlib.metadata
import importlib.util
import re
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "release.py"


def _release():
    spec = importlib.util.spec_from_file_location("release_tools", SCRIPT)
    assert spec and spec.loader, f"missing {SCRIPT}"
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── D1: one source per ecosystem ────────────────────────────────────────────


def test_the_python_package_takes_its_version_from_cargo():
    text = (ROOT / "pyproject.toml").read_text()
    project = text[text.index("[project]") : text.index("\n[", text.index("[project]") + 1)]
    assert re.search(r'^dynamic\s*=\s*\[[^\]]*"version"', project, re.M), "version must be dynamic"
    assert not re.search(r"^version\s*=", project, re.M), "no static version in [project]"


def test_the_module_version_is_the_installed_package_version():
    import gbcms

    assert gbcms.__version__ == importlib.metadata.version("gbcms")
    source = (ROOT / "src" / "gbcms" / "__init__.py").read_text()
    assert not re.search(r'__version__\s*=\s*["\']\d', source), "no version literal in __init__"


def test_nextflow_reads_its_image_tag_and_banner_from_the_manifest():
    modules = sorted((ROOT / "nextflow" / "modules" / "local" / "gbcms").glob("*/main.nf"))
    containers = [m for m in modules if "container" in m.read_text()]
    assert containers
    for m in containers:
        line = next(ln for ln in m.read_text().splitlines() if ln.strip().startswith("container"))
        assert "ghcr.io/msk-access/gbcms:${workflow.manifest.version}" in line, f"{m}: {line}"
    assert "gbcms v${workflow.manifest.version}" in (ROOT / "nextflow" / "main.nf").read_text()
    literal = re.compile(r"ghcr\.io/msk-access/gbcms:\d|gbcms v\d")
    stale = [p for p in (ROOT / "nextflow").rglob("*.nf") if literal.search(p.read_text())]
    assert not stale, stale


def test_the_version_check_passes_on_this_tree():
    assert _release().check(ROOT) == []


def _tree(tmp_path, cargo="6.6.0", lock=None, manifest="6.6.0", changelog=None, module=None):
    """A minimal release tree with the files the check reads."""
    (tmp_path / "rust").mkdir(parents=True)
    (tmp_path / "rust" / "Cargo.toml").write_text(
        f'[package]\nname = "gbcms_rs"\nversion = "{cargo}"\n'
    )
    (tmp_path / "rust" / "Cargo.lock").write_text(
        f'[[package]]\nname = "gbcms_rs"\nversion = "{lock or cargo}"\n'
    )
    (tmp_path / "pyproject.toml").write_text(
        '[project]\nname = "gbcms"\ndynamic = ["version"]\n\n[tool.maturin]\n'
    )
    (tmp_path / "src" / "gbcms").mkdir(parents=True)
    (tmp_path / "src" / "gbcms" / "__init__.py").write_text(
        'from importlib.metadata import version\n__version__ = version("gbcms")\n'
    )
    mod = tmp_path / "nextflow" / "modules" / "local" / "gbcms" / "dna"
    mod.mkdir(parents=True)
    (mod / "main.nf").write_text(
        module
        or 'process GBCMS_DNA {\n    container "ghcr.io/msk-access/gbcms:${workflow.manifest.version}"\n}\n'
    )
    (tmp_path / "nextflow" / "main.nf").write_text(
        'workflow {\n    log.info "gbcms v${workflow.manifest.version} — Nextflow Pipeline"\n}\n'
    )
    (tmp_path / "nextflow" / "nextflow.config").write_text(
        f"manifest {{\n    name = 'gbcms'\n    version = '{manifest}'\n}}\n"
    )
    (tmp_path / "CHANGELOG.md").write_text(
        changelog
        or "# Changelog\n\n## [Unreleased]\n\n## [6.6.0] - 2026-10-20 — graded evidence\n\n- x\n\n"
        "## [6.5.0] - 2026-09-25\n\n- y\n"
    )
    return tmp_path


def test_the_release_check_passes_when_every_source_agrees_with_the_tag(tmp_path):
    assert _release().check(_tree(tmp_path), tag="6.6.0") == []


@pytest.mark.parametrize(
    "kwargs,needle",
    [
        ({"cargo": "6.6.0-dev.0"}, "Cargo.toml"),  # a dev version under a release tag
        ({"lock": "6.5.0"}, "Cargo.lock"),  # lock not refreshed
        ({"manifest": "6.5.0"}, "nextflow.config"),  # pipeline still on the old image
        ({"changelog": "# Changelog\n\n## [Unreleased]\n\n- x\n"}, "CHANGELOG"),  # no section
        ({"module": 'process P {\n    container "ghcr.io/msk-access/gbcms:6.5.0"\n}\n'}, "main.nf"),
    ],
)
def test_the_release_check_names_what_disagrees(tmp_path, kwargs, needle):
    problems = _release().check(_tree(tmp_path, **kwargs), tag="6.6.0")
    assert problems and any(needle in p for p in problems), problems


def test_the_release_check_refuses_a_tag_that_is_not_a_plain_version(tmp_path):
    assert _release().check(_tree(tmp_path), tag="v6.6.0")


def test_the_release_guide_lists_exactly_the_checked_locations():
    guide = (ROOT / "docs" / "development" / "release-guide.md").read_text()
    section = guide[
        guide.index("## Version Locations") : guide.index(
            "\n## ", guide.index("## Version Locations") + 5
        )
    ]
    rows = {m.group(1) for m in re.finditer(r"^\| `([^`]+)` \|", section, re.M)}
    assert rows == {loc for loc, _ in _release().LOCATIONS}, rows


# ── D2: the Release page comes from the CHANGELOG ───────────────────────────


def test_release_notes_and_title_come_from_the_changelog_section(tmp_path):
    rel = _release()
    root = _tree(tmp_path)
    title, body = rel.notes(root, "6.6.0")
    assert title == "6.6.0 — graded evidence"
    assert body.strip() == "- x"
    (root / "CHANGELOG.md").write_text("## [6.5.0] - 2026-09-25\n\n### Fixed\n- y\n")
    title, body = rel.notes(root, "6.5.0")
    assert title == "6.5.0" and body.strip() == "### Fixed\n- y"


def _jobs(name):
    return yaml.safe_load((ROOT / ".github" / "workflows" / name).read_text())["jobs"]


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


def _runs(job):
    return "\n".join(str(s.get("run", "")) + str(s.get("uses", "")) for s in job.get("steps", []))


def test_the_release_workflow_publishes_nothing_before_the_version_check():
    jobs = _jobs("release.yml")
    assert "scripts/release.py check --tag" in _runs(jobs["verify"])
    for job in ("linux", "sdist", "release", "docker"):
        assert "verify" in _ancestors(jobs, job), job


def test_every_pr_runs_the_version_check():
    assert any("scripts/release.py check" in _runs(j) for j in _jobs("test.yml").values())


def test_the_nextflow_lint_job_inspects_the_resolved_containers():
    runs = "\n".join(_runs(j) for j in _jobs("nextflow-lint.yml").values())
    assert "nextflow inspect" in runs and "manifest" in runs


def test_the_release_workflow_creates_the_release_page_with_attested_artifacts():
    jobs = _jobs("release.yml")
    job = jobs["github-release"]
    assert {"release", "docker"} <= _ancestors(jobs, "github-release")
    steps = _runs(job)
    assert "scripts/release.py notes" in steps
    assert "actions/attest-build-provenance" in steps
    assert "sha256sum" in steps and "gh release" in steps
    perms = job.get("permissions", {})
    assert perms.get("contents") == "write" and perms.get("attestations") == "write"
    assert perms.get("id-token") == "write"


# ── P3: the bin window is a floor ───────────────────────────────────────────


def test_the_architecture_page_says_the_bin_window_is_a_floor():
    page = (ROOT / "docs" / "reference" / "architecture.md").read_text()
    row = next(ln for ln in page.splitlines() if ln.startswith("| `BIN_WINDOW`"))
    assert "Maximum span" not in row and "floor" in row.lower(), row
    assert "Parity testing" not in page


# ── review round ────────────────────────────────────────────────────────────


def test_develop_cannot_point_the_pipeline_at_an_unpublished_image(tmp_path):
    # 6.6.0 is newer than 6.6.0-dev.0: that image does not exist until the release.
    problems = _release().check(_tree(tmp_path, cargo="6.6.0-dev.0", manifest="6.6.0"))
    assert any("nextflow.config" in p for p in problems), problems


def test_a_develop_tree_on_the_last_released_image_passes(tmp_path):
    assert _release().check(_tree(tmp_path, cargo="6.6.0-dev.0", manifest="6.5.0")) == []


def test_develop_must_move_to_the_next_dev_version_after_a_release(tmp_path):
    # After the back-merge, develop at a released 6.6.0 would stamp every dev build
    # "gbcms v6.6.0". A release branch at 6.6.0 is fine.
    rel = _release()
    tree = _tree(tmp_path, cargo="6.6.0", manifest="6.6.0")
    assert any("develop" in p for p in rel.check(tree, branch="develop"))
    assert rel.check(tree, branch="release/6.6.0") == []


def test_the_release_check_refuses_a_prerelease_tag_even_when_cargo_matches(tmp_path):
    problems = _release().check(_tree(tmp_path, cargo="6.6.0-rc.1"), tag="6.6.0-rc.1")
    assert any("bare X.Y.Z" in p for p in problems), problems


def test_container_options_and_comments_are_not_version_locations(tmp_path):
    module = (
        "// gbcms v6.5.0 changed the merge columns\n"
        "process GBCMS_DNA {\n"
        '    container "ghcr.io/msk-access/gbcms:${workflow.manifest.version}"\n'
        '    containerOptions "--user root"\n'
        "}\n"
    )
    assert _release().check(_tree(tmp_path, module=module), tag="6.6.0") == []


def test_notes_stop_at_any_heading_and_drop_trailing_link_references(tmp_path):
    rel = _release()
    root = _tree(
        tmp_path,
        changelog="## [6.6.0] - 2026-10-20 — graded evidence\n\n- x\n\n"
        "## [2.1.1] - 2025-11-25 [YANKED]\n\n- old\n\n[6.6.0]: https://example.org/6.6.0\n",
    )
    _, body = rel.notes(root, "6.6.0")
    assert body.strip() == "- x", body
    root2 = _tree(
        tmp_path / "b", changelog="## [6.6.0] - 2026-10-20\n\n- x\n\n[6.6.0]: https://example.org\n"
    )
    assert rel.notes(root2, "6.6.0")[1].strip() == "- x"


def test_long_notes_are_cut_at_a_section_with_a_link_to_the_changelog(tmp_path):
    rel = _release()
    long = "\n".join(f"### Part {k}\n" + "- line\n" * 400 for k in range(60))
    root = _tree(tmp_path, changelog=f"## [6.6.0] - 2026-10-20\n\n{long}\n")
    _, body = rel.notes(root, "6.6.0")
    assert len(body) <= rel.NOTES_LIMIT and "CHANGELOG.md" in body.splitlines()[-1]


def test_every_listed_location_exists():
    for loc, _ in _release().LOCATIONS:
        assert list(ROOT.glob(loc)), loc


def test_no_workflow_pastes_a_ref_name_into_a_shell_script():
    for wf in ("release.yml", "test.yml", "nextflow-lint.yml"):
        for job in _jobs(wf).values():
            for step in job.get("steps", []):
                assert "github.ref_name" not in str(step.get("run", "")), (wf, step.get("name"))


def test_only_a_tag_push_publishes_and_the_default_token_is_read_only():
    wf = yaml.safe_load((ROOT / ".github" / "workflows" / "release.yml").read_text())
    assert wf.get("permissions") == {"contents": "read"}
    jobs = wf["jobs"]
    for job in ("release", "github-release"):
        assert "github.event_name == 'push'" in str(jobs[job].get("if", "")), job
    push = next(s for s in jobs["docker"]["steps"] if "build-push-action" in str(s.get("uses", "")))
    assert "github.event_name == 'push'" in str(push["with"]["push"])


def test_ci_installs_what_the_tests_import():
    # CI installs pyproject's `test` dependency group, which lists what the suite
    # imports beyond gbcms's own dependencies (yaml for these workflow tests).
    install = _runs(_jobs("test.yml")["test"])
    assert "--group test" in install
    toml = (ROOT / "pyproject.toml").read_text()
    group = re.search(r"^test = \[\n(.*?)^\]", toml, re.M | re.S).group(1)
    names = set(re.findall(r'^\s*"([A-Za-z0-9_-]+)', group, re.M))
    assert {"pyyaml", "pyarrow", "pytest-mock"} <= names, names
