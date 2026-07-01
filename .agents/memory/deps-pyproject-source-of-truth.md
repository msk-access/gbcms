---
name: deps-pyproject-source-of-truth
description: pyproject.toml is the single source of truth for deps; CI/Docker bypass any lockfile, so every directly-imported package must be declared.
metadata:
  type: project
---

`pyproject.toml` `[project.dependencies]` (+ the `[dev]` extra) is the **single
source of truth** for gbcms dependencies. `uv.lock` was removed (2026-06-26) and is
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
- Reproducibility for clinical/validated builds belongs in the **Docker image**
  (pinned final-stage install), not a dev lockfile every path bypasses.
  See [[nextflow-cli-default-divergence]] for the related "keep deploy config in sync".
