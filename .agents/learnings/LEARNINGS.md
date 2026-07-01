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

## [LRN-20260627-001] rule-body | code comments & logs must explain behavior, not cite ticket labels
- **Status:** resolved (rule promoted); cleanup of existing labels tracked as plan DX-1
- **Cause:** rule-body
- **Summary:** I annotated code comments and `debug!`/`warn!` strings with ticket
  labels (`CR-1`, `HI-11`, `ME-8`) — and the pre-existing code already did this with
  `P4c`. A reader with no plan/PR context can't decode those labels, so the comment
  explains nothing about the actual behavior. The code-quality rule said "add a why
  comment" but never said the why must be the *reason*, not the ticket that prompted it.
- **Promotion target:** `.agents/rules/code-quality.md` — `DONE:` added a "Comment &
  Log Hygiene" section: no ticket/milestone labels in code or logs; a why-comment
  states the reason (`✗ // ME-8: padding fix` → `✓ // exclude no-test variants; they
  inflate the FDR family`). Existing labels swept under plan ticket DX-1.
- **Related:** [[no-ticket-labels-in-code]] memory; sibling to the "update commenting"
  QC checklist item already in that rule file.

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
