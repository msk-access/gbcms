"""The image is published for amd64 and arm64, built natively, and CI runs both.

On Apple silicon, Docker runs an amd64-only image under emulation (the Nextflow `docker`
profile forces `--platform linux/amd64`). Contracts:

- the release builds each architecture on a runner of that architecture (no QEMU), pushes
  it by digest, and one job joins the digests into the release tag and `latest`;
- nothing is published except by a tag push, and the release page waits for the join;
- the test workflow builds and runs the image on both architectures, checks that the image
  is the runner's architecture, and requires the same counting-run output from both.
"""

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
WF = ROOT / ".github" / "workflows"


def _jobs(name):
    return yaml.safe_load((WF / name).read_text())["jobs"]


def _legs(job):
    return {leg["arch"]: leg["runner"] for leg in job["strategy"]["matrix"]["include"]}


def _runs(job):
    return "\n".join(str(step.get("run", "")) for step in job["steps"])


def test_each_architecture_is_built_on_its_own_runner():
    for name, job in (("release.yml", "docker"), ("test.yml", "docker-test")):
        legs = _legs(_jobs(name)[job])
        assert set(legs) == {"amd64", "arm64"}, (name, legs)
        assert legs["arm64"].endswith("-arm"), (name, legs)
        assert not legs["amd64"].endswith("-arm"), (name, legs)
    for name in ("release.yml", "test.yml"):
        text = (WF / name).read_text()
        assert "setup-qemu-action" not in text, f"{name} emulates an architecture"
        assert "platforms:" not in text and "--platform" not in text, name


def test_the_release_joins_the_digests_and_only_a_tag_push_publishes():
    jobs = _jobs("release.yml")
    build = jobs["docker"]
    push = next(s for s in build["steps"] if "build-push-action" in str(s.get("uses", "")))
    assert "push-by-digest=true" in push["with"]["outputs"]
    assert "github.event_name == 'push'" in str(push["with"]["push"])
    merge = jobs["docker-merge"]
    assert merge["needs"] == ["docker"] or merge["needs"] == "docker", merge["needs"]
    assert "github.event_name == 'push'" in merge["if"] and "refs/tags/" in merge["if"]
    runs = _runs(merge)
    assert "docker buildx imagetools create" in runs
    # a tag with one architecture is not multi-arch
    assert "-eq 2" in runs
    tags = next(s for s in merge["steps"] if "metadata-action" in str(s.get("uses", "")))
    assert (
        "type=ref,event=tag" in tags["with"]["tags"]
        and "type=raw,value=latest" in tags["with"]["tags"]
    )


def test_the_release_page_and_the_summary_wait_for_the_joined_image():
    jobs = _jobs("release.yml")
    assert "docker-merge" in jobs["github-release"]["needs"]
    assert "docker-merge" in jobs["summary"]["needs"]
    summary = _runs(jobs["summary"]) + str(jobs["summary"]["steps"])
    assert "needs.docker-merge.result" in summary


def test_ci_runs_the_image_on_both_architectures_and_compares_the_outputs():
    jobs = _jobs("test.yml")
    runs = _runs(jobs["docker-test"])
    assert "docker image inspect gbcms:test --format '{{.Architecture}}'" in runs
    assert "wheel_fixture_run.py" in runs
    compare = jobs["docker-compare"]
    assert compare["needs"] == "docker-test" or compare["needs"] == ["docker-test"]
    assert "diff -u amd64.txt arm64.txt" in _runs(compare)
    # Apptainer: the clusters run x86_64
    apptainer = [s for s in jobs["docker-test"]["steps"] if "pptainer" in s.get("name", "")]
    assert apptainer and all(s.get("if") == "matrix.arch == 'amd64'" for s in apptainer)
