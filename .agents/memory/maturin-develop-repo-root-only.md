---
name: maturin-develop-repo-root-only
description: "Run `maturin develop --release` from the repo root only. Passing `-m rust/Cargo.toml` bypasses [tool.maturin] in pyproject and installs a useless top-level `_rs` package while the stale src/gbcms/_rs.*.so keeps being imported."
metadata: 
  node_type: memory
  type: project
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
  modified: 2026-09-24T04:38:14.603Z
---

`pyproject.toml` `[tool.maturin]` sets `python-source = "src"` and
`module-name = "gbcms._rs"`. Invoking `maturin develop -m rust/Cargo.toml`
ignores that config (it reads only the crate), installs `site-packages/_rs/`,
and **never refreshes `src/gbcms/_rs.cpython-*.so`** — which the `gbcms.pth`
editable resolves first. Result: builds "succeed" while imports keep running a
stale extension (this masked a fix for a month in 2026-08/09; discovered while
fixing issue #89).

**How to apply:** always `VIRTUAL_ENV=$PWD/.venv uvx maturin develop --release`
from the repo root. After any rebuild, verify freshness:
`python -c "import gbcms._rs as m, os, datetime; print(m.__file__, datetime.datetime.fromtimestamp(os.path.getmtime(m.__file__)))"`
— the path must be `src/gbcms/_rs...so` with a just-now mtime. If an old .so
shadows, delete it and rebuild from the root.

The main checkout's `.venv` has **no maturin of its own**: `.venv/bin/maturin`
fails "no such file", and a `-q` build piped through `grep error` hides that, so
tests silently run the old engine (happened 2026-09-24; caught by the mtime).
Worktree venvs built with `pip install maturin` do have it. Never filter build
output without also checking the .so mtime.
