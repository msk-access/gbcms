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

## [LRN-20261002-003] environment | a detached harness run is invisible to the operator
- **Status:** resolved (memory)
- **Cause:** environment
- **Summary:** acceptance runs started with `nohup` (to outlive the 30-minute
  tracked-task limit) do not show in the app's background list; the operator saw
  "nothing running" mid-acceptance. A full `--trace` probe run also ground for an
  hour before I switched to a targeted trace.
- **Promotion target:** memory `long-runs-visible` — `DONE:` pair each detached
  run with a tracked waiter; say why a run is slow and switch methods.
- **Related:** [[long-runs-visible]].

## [LRN-20261002-002] rule-body | a read-input rule is validated against an oracle that does not share it
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** the first C17 build (clip bases past the TLEN fragment end) agreed
  with the read census on every SNV row, because the census was changed in step
  and shares the rule. Tracing the three indel rows that moved away from the
  census showed 13 genuine FLT3-ITD ALT reads lost: TLEN, a reference distance,
  leaves out inserted bases. The process named the census as the oracle without
  saying it cannot judge a rule it implements.
- **Promotion target:** `docs/reference/read-judgment.md` "How a change to read
  judgment is made", step 1 — `DONE:` for a read-input rule the census mirrors,
  the evidence includes an oracle it does not share. Memory
  `census-mirrors-read-inputs`.
- **Related:** [[census-mirrors-read-inputs]], [[holistic-effects-map]].

## [LRN-20261002-001] rule-body | a community-practice survey covers several tools
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** for group 1 (read inputs) I surveyed only GATK. The operator
  pointed out there are more tools to look at. The validation standard said
  "callers and genotypers" but named none, so the survey defaulted to one.
- **Promotion target:** `docs/reference/read-judgment.md` "How a change to read
  judgment is made", step 1 — `DONE:` names the tool set to survey. Memory
  `survey-several-tools`.
- **Related:** [[survey-several-tools]].

## [LRN-20261001-001] rule-body | a failing test's expectation is changed only with the operator's say
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** during H3 I edited existing test expectations as part of fixes
  (three `test_shifted_indels` cases to depth-only, the C28 example, a census
  assertion), explaining them only in commits and the CHANGELOG. The operator
  asked to be told why before any test is modified, or to measure to understand.
  The validation standard said "red-first" but nothing about changing existing
  expectations, so the rule was missing.
- **Promotion target:** `AGENTS.md` "Counting test invariants" — `DONE:` a line:
  never edit an existing test's expectation to make it pass without telling the
  operator what it asserted, why it is wrong, and the evidence, then waiting.
  Memory `test-changes-need-operator-notice`.
- **Related:** [[test-changes-need-operator-notice]].

## [LRN-20260928-001] rule-body | gbcms is a genotyper, not a caller; the BAM is not the truth
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** the standing principle was written as "the BAM is truth", and I
  repeated it in the validation standard. The operator corrected it: alignments
  can be wrong, and the point is that gbcms is a genotyper, not a caller. The
  evidence is each read's own bases, judged independently of the aligner's
  placement (this is what the exact-carrier windows, repeat growth, clip reading
  and splice handling already do); gbcms counts the given allele and never decides
  existence or rewrites an allele.
- **Promotion target:** `AGENTS.md` invariant 7 — `DONE:` "Count the given allele — a
  genotyper, not a caller … its own bases carry that ALT (alignments can be wrong;
  judge bases, not placement)". Memory `bam-is-truth` renamed and rewritten as
  `genotyper-not-caller`; the plan's "Validation standard" and the add-feature skill
  reworded. Dataset specifics moved to local-only memories (operator: keep them local).
- **Related:** [[genotyper-not-caller]], [[count-the-given-allele]].

## [LRN-20260925-001] rule-body | a census finding is not ALT evidence
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** For C1 I proposed option C: turn fallback ALT calls without ALT sequence
  into partial, *unless* reads share a recurring unannotated haplotype. That exception
  credits a mis-described event to the row's ALT (deconvolution). `bam-is-truth` said a
  recurrent unannotated haplotype "is a finding" without saying it must never count
  for the row, and no invariant stated that gbcms counts only the given allele. The
  operator restated the principle: input is taken as correct, results are accurate for
  it, caveats go to the status-reason and diagnostic columns.
- **Promotion target:** `AGENTS.md` invariant 7 "Count the given allele";
  `.agents/memory/bam-is-truth.md` clarifies that a finding is never ALT evidence; new
  memory `count-the-given-allele.md`; `CYCLE_6.6.0_PLAN.md` C1 drops the exception,
  and C3 is re-scoped to report the given allele, with the default-on decomposition
  as an open decision.

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
