---
description: Git commit workflow with pre-commit checks
---

# /git-commit

1. Stage changes:
   ```bash
   git add -A
   git status --short
   ```

2. Verify no data files staged:
   ```bash
   git diff --cached --name-only | grep -E '\.(parquet|tsv|csv|bam|bai|fasta|fa)$'
   ```

3. Run lint suite:
   ```bash
   ruff check src/ tests/
   black --check src/ tests/
   mypy src/
   ```

4. Run tests:
   ```bash
   python -m pytest tests/ --no-header -q
   ```

5. Commit with conventional format:
   ```bash
   git commit -m "<type>(<scope>): <description>"
   ```

6. Push:
   ```bash
   git push origin <branch>
   ```
