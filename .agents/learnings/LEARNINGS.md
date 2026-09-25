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

## [LRN-20260923-002] rule-body | reproduce a claim before stating it
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** In the T7 rescue work I wrote "no_improvement is unreachable" into code,
  docs and the skill without trying a counterexample (a reviewer found one: indel-disrupted
  reads), and reported a review finding ("rescue label contig differs from the row") as
  general when running it showed it only happens for chr-prefixed MAF input. The step
  review checklist asked for comments/logs/monitoring but never for evidence behind claims.
- **Promotion target:** `.agents/rules/code-quality.md` — `DONE:` "AFTER implementing each
  step" gains "Reproduce every claim before stating it … otherwise say it is unverified".

## [LRN-20260923-001] rule-body | real-data acceptance must cover every assay the change reaches
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** The T7 rescue gate was validated on two IMPACT cohorts (68 samples) and
  presented as done; the ACCESS duplex/simplex study (the fragment-scored assay, combined
  by `gbcms merge`) then exposed a germline-SNP adoption (BRCA2, combined fragment ALT
  25 → 724) and duplex/simplex rescue conflicts. Each fix layered another heuristic
  (loosened gate → confirmed guard → error allowance) without first checking it against
  the operator's principle (the BAM is truth), which read as a rabbit hole. The add-feature
  procedure had lint/test QA but no real-data acceptance step at all.
- **Promotion target:** `.agents/skills/add-feature/SKILL.md` step 6 — `DONE:` "Real-data
  acceptance … every assay the change reaches (IMPACT tumour + matched normal; ACCESS
  duplex + simplex through `gbcms merge`) … score against the reads".
- **Related:** [[bam-is-truth]] memory.

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
