---
name: worktree-tests-need-isolated-venv
description: "In a .claude/worktrees/* checkout, the shared .venv (and mambaforge python) import gbcms from the MAIN checkout's src/ — tests silently exercise the wrong code. Build an isolated venv per worktree."
metadata:
  node_type: memory
  type: project
  originSessionId: 890ffab2-2fbd-474b-915c-62b9807e66d3
  modified: 2026-09-23T14:09:13.842Z
---

The main checkout's `.venv` has an editable `gbcms` install pointing at the
main checkout's `src/` (and its `_rs.so`). Running pytest
or the `gbcms` CLI with it from a worktree tests the main checkout's branch, not
the worktree's edits — green results mean nothing. Running `maturin develop`
against that shared venv from a worktree would instead repoint it for every
other session.

**How to apply:** per worktree, in the session scratchpad:
`uv venv --python 3.10 $SP/venv`, `VIRTUAL_ENV=$SP/venv uv pip install` the
runtime deps + pytest/pytest-cov/pytest-mock/pyarrow/black/ruff/mypy/
types-pyyaml/scipy-stubs/maturin (CI installs these unpinned → latest = CI's
versions), then from the worktree root
`VIRTUAL_ENV=$SP/venv PATH=$SP/venv/bin:$PATH maturin develop --release`
(repo root only — [[maturin-develop-repo-root-only]]). Verify with
`$SP/venv/bin/python -c "import gbcms; print(gbcms.__file__)"` → worktree path.
For `mkdocs build --strict`, the git-revision-date plugin warns in a worktree
(.git is a file) — build with a temp config minus that plugin. See
[[black-version-skew-venv-vs-ci]].
