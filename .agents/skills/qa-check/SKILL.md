---
name: qa-check
description: "Procedure for the gbcms pre-merge DEV QA sweep — ruff, black, mypy, pytest, cargo test, clippy, and a staged-data-file guard. Use before merging or when asked to QA a change. NOTE: this is developer QA (lint/test). For genomic read-level QC (filters, BQ gating, BAQ) see the read-filters-qc skill."
---

# qa-check

Run the full lint/test gate. Fix order: ruff → black → mypy → clippy.

Run all of it for a docs-, skills- or Nextflow-only change too: contract tests read
`docs/` (every `docs/**/*.md`), `nextflow/` and the README, so a moved or reworded
file can fail them.

```bash
ruff check src/ tests/
black --check src/ tests/
mypy src/
python -m pytest tests/ --no-header -q
cd rust && cargo test && cargo clippy --all-targets -- -D warnings && cd ..
git status && git diff --stat
git diff --cached --name-only | grep -Ei '\.(parquet|bam|bai|cram|crai|tsv|csv|maf|dat|vcf|bed|fastq|fq)(\.(gz|bgz))?$' | grep -v '^tests/testdata/'   # must be empty (the pre-commit hook enforces it)
```

Report on pass:
```
✅ QA PASS — ruff clean · black formatted · mypy 0 · pytest pass · cargo test pass · clippy clean · tree clean
```

(Exact test counts drift; assert "all pass," not a hardcoded number.)
