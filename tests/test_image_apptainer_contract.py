"""CI runs the image the way a RHEL cluster does: through Apptainer.

RHEL clusters run the image with Singularity or Apptainer (the Nextflow `singularity`,
`slurm` and `iris` profiles), never with Docker. The two differ where it matters: the
container filesystem is read only, the user is the caller, and the caller's HOME and
environment come into the container. A Python in the image imports packages from the
caller's `~/.local` first, which broke the image for anyone with `pip install --user`
packages until the image turned the user site off. Contracts:

- the image sets PYTHONNOUSERSITE;
- the job that builds the image installs Apptainer, builds a SIF from the Docker image,
  runs the CLI and a counting run in it with a hostile package in the caller's HOME, and
  fails unless the outputs equal the Docker run's.
"""

import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def _image_job():
    wf = yaml.safe_load((ROOT / ".github" / "workflows" / "test.yml").read_text())
    return wf["jobs"]["docker-test"]


def _runs(job):
    return "\n".join(str(step.get("run", "")) for step in job["steps"])


def test_the_image_turns_the_users_python_packages_off():
    dockerfile = (ROOT / "Dockerfile").read_text()
    runtime = dockerfile[dockerfile.rindex("\nFROM ") :]
    assert re.search(r"^ENV PYTHONNOUSERSITE=1$", runtime, re.M), runtime


def test_ci_runs_the_image_through_apptainer_and_compares_with_docker():
    runs = _runs(_image_job())
    assert "ppa:apptainer/ppa" in runs and "apt-get install -y apptainer" in runs
    # the SIF is built from the image the job just built, as a cluster pulls one
    assert "apptainer build" in runs and "docker-daemon://gbcms:test" in runs
    # as Nextflow runs a task: exec in the SIF, through bash, with the data folders bound
    assert "apptainer exec" in runs and "-B " in runs and "/bin/bash -c" in runs
    # the counting run and the CLI start, in both runtimes
    assert "gbcms --version" in runs and "wheel_fixture_run.py" in runs
    assert "docker run" in runs
    # the outputs must be equal, and the check must be able to fail
    assert re.search(r"^\s*diff -u work/docker\S* work/apptainer\S*$", runs, re.M), runs


def test_ci_puts_a_hostile_package_in_the_callers_home_before_the_apptainer_run():
    runs = _runs(_image_job())
    decoy = runs.index(".local/lib/python3.11/site-packages/pydantic")
    assert decoy < runs.index("apptainer exec"), "the decoy must exist before the run"
    # Apptainer binds HOME by default; the Docker reference run does not see it
    assert "--no-home" not in runs and "--containall" not in runs
