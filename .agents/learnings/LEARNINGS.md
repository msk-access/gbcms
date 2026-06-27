# LEARNINGS

The work-order log. Memory captures the durable note ("don't repeat this"); this
log names the **promotion target** — the exact file to edit so the cause cannot
recur. The two load-bearing fields are **Cause** (roots the fix) and **Promotion
target** (file + change).

Attribution — classify every failure into exactly ONE cause before touching a file:
- `rule-body` — right rule/skill/hook fired, its logic was wrong/stale → edit that file's body
- `rule-trigger` — wrong one fired, or the right one stayed silent → edit its description/matcher
- `rule-permission` — wrong tool access (too broad or missing) → edit allowed-tools / matcher scope
- `environment` — rule+trigger were fine, failure was external → NOWHERE in the harness; log to ERRORS.md

Newest at the top. Before promoting an edit, grep `REJECTED.md` for the target +
topic — don't re-propose a vetoed dead end. Promote by delta edit, never by
regenerating a section.

---

## [LRN-20260626-001] rule-body | a guard hook must scope its target to the command segment, not the whole command string
- **Status:** resolved
- **Cause:** rule-body
- **Summary:** `block-push-main.py` matched the literal token `main` anywhere in
  the Bash command and denied a legitimate `git push` to a feature branch because
  the word "main" appeared in the **commit message** of the same compound command.
  A guard that keys on a token anywhere in the string false-positives under
  realistic compound commands.
- **Promotion target:** `.claude/hooks/block-push-main.py` — `DONE:` isolate the
  `git push …` segment (command-boundary anchored, up to the next separator) and
  require `main` as a standalone ref token *within that segment*
  (`(?<![\w-])main(?![\w-])`), so commit-message text and `main-fix`/`maint`
  branches don't trip it.
- **Related:** the claudelicious hook principle "anchor at command boundaries"
  (`.agents/rules/security.md`); [[deps-pyproject-source-of-truth]] is a sibling
  "the obvious check missed the real path" lesson.
