---
description: Full quality assurance sweep before merging to main
---

# /qa-check

1. Run ruff:
   ```bash
   ruff check src/ tests/
   ```

2. Check formatting:
   ```bash
   black --check src/ tests/
   ```

3. Type checking:
   ```bash
   mypy src/
   ```

4. Run Python tests (255 expected):
   ```bash
   python -m pytest tests/ --no-header -q
   ```

5. Run Rust tests (143 expected):
   ```bash
   cd rust && cargo test && cd ..
   ```

6. Rust linting:
   ```bash
   cd rust && cargo clippy -- -D warnings && cd ..
   ```

7. Check for uncommitted changes:
   ```bash
   git status
   git diff --stat
   ```

8. Verify no data files accidentally staged:
   ```bash
   git diff --cached --name-only | grep -E '\.(parquet|tsv|csv|bam|bai)$'
   ```

9. If all pass, report:
   ```
   ✅ QA PASS
   - ruff: clean
   - black: formatted
   - mypy: 0 errors
   - pytest: 255/255 passed
   - cargo test: 143/143 passed
   - clippy: clean
   - git: clean working tree
   ```
