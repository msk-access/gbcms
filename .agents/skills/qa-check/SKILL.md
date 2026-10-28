---
name: qa-check
description: "Procedure for the gbcms pre-commit and pre-merge DEV QA gate — ruff, black (CI's version), mypy, clippy, cargo test, pytest (which includes the agent-skill checks) and a staged-data-file guard — with every exit code read. Use before every commit or merge, and when asked to QA a change. NOTE: this is developer QA (lint/test). For genomic read-level QC (filters, BQ gating, BAQ) see the read-filters-qc skill."
---

# qa-check

The one home of the lint and test gate. Fix order: ruff → black → mypy → clippy.

Run all of it for a docs-, skills- or Nextflow-only change too: contract tests read
`docs/` (every `docs/**/*.md`), `nextflow/`, the README and `.agents/` (skill frontmatter,
and every code name a skill or rule mentions), so a moved or reworded file can fail them.

```bash
ruff check src/ tests/;                          echo "ruff rc=$?"
uvx --from 'black>=26.5' black --check src/ tests/; echo "black rc=$?"   # CI's black, not the venv's
mypy src/;                                       echo "mypy rc=$?"
(cd rust && cargo clippy --all-targets -- -D warnings; echo "clippy rc=$?"; cargo test; echo "cargo test rc=$?")
python -m pytest tests/ --no-header -q > /tmp/pytest.log 2>&1; echo "pytest rc=$?"; tail -3 /tmp/pytest.log
git status && git diff --stat
git diff --cached --name-only | grep -Ei '\.(parquet|bam|bai|cram|crai|tsv|csv|maf|dat|vcf|bcf|sam|bed|fastq|fq)(\.(gz|bgz))?$' | grep -v '^tests/testdata/'   # must be empty (the pre-commit hook enforces it)
```

**Read every exit code; never gate on a pipe.** `pytest | tail && git commit` commits when
pytest fails, because the pipeline's status is `tail`'s (it happened: a failing fixture was
committed). Redirect to a log and test `$?`, or `set -o pipefail`. The same holds for
`cargo test | grep`.

**Lint tools** run at CI's versions. The shared `.venv`'s black (25.9) lags CI's latest and
disagrees on some formatting; CI installs the newest from the `dev` group. A fresh scratch venv
from that group (rust-python-ffi) has the right ones; `uvx --from 'black>=26.5' black` works
anywhere. For clippy, `rustup update stable` first (CI uses stable; a stale local clippy missed
a lint that failed CI under `-D warnings`). `ruff` and `mypy` have shown no skew.

Report on pass:
```
✅ QA PASS — ruff clean · black formatted · mypy 0 · pytest pass · cargo test pass · clippy clean · tree clean
```

(Exact test counts drift; assert "all pass," not a hardcoded number.)
