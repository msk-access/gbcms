---
name: deps-pyproject-source-of-truth
description: pyproject.toml is the single source of truth for deps (runtime floors measured, dev tools in PEP 735 groups); only the Docker image installs a lock, so every directly-imported package must be declared.
metadata:
  type: project
---

`pyproject.toml` is the **single source of truth** for gbcms dependencies: runtime in
`[project.dependencies]` (floors measured by bisection; CI's Ubuntu 3.10 leg installs
them with `uv --resolution lowest-direct`), contributors' tools in the PEP 735
`[dependency-groups]` `test` and `dev` (since 2026-10-06; the `[dev]` extra is gone). `uv.lock` was removed (2026-06-26) and is
gitignored: nothing consumed it — CI (`test.yml`), the Dockerfile, and Nextflow all
build the wheel with `maturin` and run `uv pip install dist/*.whl` / `pip install`,
which **resolve fresh from wheel metadata and bypass any lockfile**. The lock had also
drifted to gbcms v2.8.0 while the project was 5.3.0 — false reproducibility.

**Why it matters:** A dependency that is only present *transitively* can vanish on an
upstream release and break PR CI even though releases and local dev are fine. This bit
us: `tests/test_cli_dna_rna.py` does `import click` (for the `click.Group` type from
`typer.main.get_command`), but `click` was undeclared — only pulled in via `typer`. A
fresh CI resolve stopped providing it → `ModuleNotFoundError: No module named 'click'`
at pytest collection, across all platforms. Fixed by adding `click>=8.0` to
`[project.dependencies]` (PR #21 / commit on develop).

**How to apply:**
- If code (incl. tests) does `import X`, declare `X` in `pyproject.toml` — never rely
  on it arriving transitively.
- Don't reintroduce a committed lock unless you also wire CI/Docker to install *from
  it* (`uv sync --frozen`) and add a `uv lock --check` drift gate — all or nothing.
- Reproducibility for clinical/validated builds belongs in the **Docker image**, and
  it has it: `docker/requirements.lock` (uv pip compile, hashes, linux/amd64, py3.11)
  is installed with `--require-hashes --no-deps`, then `pip check` gates drift. It is
  refreshed at release (release guide step 4). The org's containers image does not
  use it ([[downstream-org-containers-modules]]).
  See [[nextflow-cli-default-divergence]] for the related "keep deploy config in sync".
