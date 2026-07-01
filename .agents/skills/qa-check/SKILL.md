---
name: qa-check
description: Procedure for the gbcms pre-merge DEV QA sweep — ruff, black, mypy, pytest, cargo test, clippy, and a staged-data-file guard. Use before merging or when asked to QA a change. NOTE: this is developer QA (lint/test). For genomic read-level QC (filters, BQ gating, BAQ) see the read-filters-qc skill.
---

# qa-check

Run the full lint/test gate. Fix order: ruff → black → mypy → clippy.

```bash
ruff check src/ tests/
black --check src/ tests/
mypy src/
python -m pytest tests/ --no-header -q
cd rust && cargo test && cargo clippy --all-targets -- -D warnings && cd ..
git status && git diff --stat
git diff --cached --name-only | grep -E '\.(parquet|tsv|csv|bam|bai)$'   # must be empty
```

Report on pass:
```
✅ QA PASS — ruff clean · black formatted · mypy 0 · pytest pass · cargo test pass · clippy clean · tree clean
```

(Exact test counts drift; assert "all pass," not a hardcoded number.)
