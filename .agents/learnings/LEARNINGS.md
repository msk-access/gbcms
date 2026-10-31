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

**Memory names in older entries.** The memory files were merged in #317, so an older entry's
`[[name]]` or `.agents/memory/<name>.md` points here:
- `working-agreements`: harness-stay-small, describe-labels-in-brackets, model-per-task
- `genotyper-not-caller`: genotyper-not-caller, count-the-given-allele, identity-band-annotation-tolerance
- `decided-counting-rules`: tolerant-deletion-deliberate, rna-mapq-default, mapq0-loci-pms2, vaf-informative-reads
- `mfsd-graded-evidence`: mfsd-graded-evidence, fragmentomics-reference
- `msk-assay-facts`: abra2-junction-placement, abra2-collapse-artifacts, access-xs1-xs2-consensus, production-aligner-settings, qcfail-flag-absent-msk-data
- `references`: lancet2-reference, downstream-org-containers-modules
- `.agents/rules/decision-protocol.md`: decisions-worked-examples, best-not-least-bad, better-than-standard, survey-several-tools, test-across-data-types, fix-the-class-not-the-row, holistic-effects-map, external-review-measure-first, parameters-are-user-settable, read-rule-register-to-the-end
- `.agents/rules/security.md`: cbioportal-stays-local
- `.agents/rules/git-conventions.md`: merge-with-gh-when-asked, github-sub-issues, commit-before-review-workflows
- `.agents/rules/code-quality.md`: no-ticket-labels-in-code, deps-pyproject-source-of-truth; test-changes-need-operator-notice is in `AGENTS.md`
- `docs/harness/measurement-harness.md`: harness-cleanup-after-merge, long-runs-visible, hpc-partitions, macos-no-timeout-command, trace-output-wraps, harness-references, other-tools-harness, vcf2maf-oracle, giab-truth-arm, pysam-validation-oracle, cluster-aware-oracle, census-mirrors-read-inputs
- `.agents/skills/qa-check` and `rust-python-ffi`: pipe-masks-gate-exit, black-version-skew-venv-vs-ci, worktree-tests-need-isolated-venv, maturin-develop-repo-root-only
- Deleted as fixed or covered: wfa-bq-gate-contract, engine-output-aware, nextflow-cli-default-divergence (AGENTS.md invariants 2 and 3, the nextflow-pipeline skill), bin-anchor-coverage (`architecture.md`), claudelicious-reference (`docs/harness/00-plan.md`)

---

## [LRN-20261030-003] rule-body | the PHI sanitizer filtered stdout but not library warnings on stderr
- **Status:** resolved (rule promoted).
- **Cause:** rule-body.
- **Summary:**
  - An agent comparing gbcms with the C++ tool wrote its sanitizer first, as the rule says, but
    passed only stdout through it. htslib's "index older than data" warning on stderr printed WES
    patient BAM file names (sample IDs) into two tool outputs: one name, then twelve. No position or
    read name was shown. The agent reported it, masked its log and silenced the warnings.
  - The rule said "every printed line" but named no stderr, and library warnings are easy to miss.
- **Promotion target:** `.agents/rules/security.md` (patient BAM slices, the sanitizer bullet) —
  `DONE:` the sanitizer covers stderr; library warnings print file names; merge stderr into the
  filtered stream or keep it in a file that never prints.
- **Second slip, same day:** another agent's masking command did not mask and printed five DMP BAM
  stems from a region list. `DONE:` the sanitizer is tested on a made-up line holding each kind of
  identifier before it sees patient data.
- **Related:** the 2026-10-29 sanitizer-first rule (#372).

## [LRN-20261030-002] rule-trigger | follow-up work went to the agent at hand, not the model the task needs
- **Status:** resolved (rule promoted).
- **Cause:** rule-trigger.
- **Summary:**
  - The operator asked (2026-10-30) whether delegation still follows the model-per-task rule.
  - New agents got explicit models. But follow-up work went to the two long-lived agents because
    they were already running: an Opus agent was given the C++ timing tables (Sonnet work).
  - No fresh adversarial review was planned before count decisions; Greptile was the only outside read.
  - The rule fired when an agent was spawned, and stayed silent when work was handed to an agent
    already running.
- **Promotion target:** `.agents/memory/working-agreements.md` ("Model per task") — `DONE:` the rule
  holds for follow-up work; routine aggregation goes to a new Sonnet or Haiku agent; a fresh Fable or
  Opus agent reviews each count-changing option before it goes to the operator.
- **Related:** [[working-agreements]].

## [LRN-20261030-001] rule-trigger | status tables and relayed agent reports went out as bare labels
- **Status:** resolved (rule promoted).
- **Cause:** rule-trigger.
- **Summary:**
  - The operator asked again (2026-10-30) for plain descriptions with labels in brackets. My
    status table listed builds, HPC jobs and runs by label only (`c61hs`, `hyR`, `htM`, item A).
  - The rule existed in `working-agreements.md`, but it stayed silent in two places. The index
    line read each turn said only "labels in brackets", and the rule's examples covered rules and
    phases, not build, job or run labels.
  - Agent reports arrive full of labels, and I relayed them without a rewrite.
  - The table also showed two memory changes the operator thought were done, with no reason given.
- **Promotion target:**
  - `.agents/memory/MEMORY.md` (the working-agreements line) — `DONE:` names every label kind and
    says the rule holds in tables and relayed reports.
  - `.agents/memory/working-agreements.md` ("Explain so the operator can follow") — `DONE:` adds
    the bullet "Build, job and run labels are labels too": rewrite relays, ask agents for the same
    form, and say why work is back.
- **Related:** [[working-agreements]].

## [LRN-20261026-001] rule-body | an amended read rule was read from its superseded text
- **Status:** open: promotion waits on the operator re-confirming RJ-20.
- **Cause:** rule-body.
- **Summary:**
  - On #304 I said the sequence judge "breaks RJ-20" (an unreadable insert of the ALT's length should be neither + partial), queued a port, recorded it on #304 and in CONTINUITY, and wrote it into #329's Read Filters text.
  - RJ-20 was amended on 2026-10-15 (#278): such an insert inside the discrimination window *is* the ALT. Both paths follow the amended rule, and `tests/test_unreadable_insert_contract.py` pins it.
  - The register row states the old rule first and the amendment last. `allele-classification.md` (§Insertion, the wrong-length table, the flowchart) and `glossary.md` (`has_unreadable_insert`) were never updated. The fix round caught it before any code changed.
- **Promotion target:**
  - `docs/reference/read-judgment.md`: an amended row leads with the current rule, history after.
  - `docs/reference/allele-classification.md`, `glossary.md` and `read-filters.md` (#329): state the current rule.
  - Memory `read-rule-register-to-the-end`: read the whole row and the contract test's docstring before stating a rule (done).

## [LRN-20261024-001] rule-body | a docs-only change still runs the whole test suite
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** PR #321 (harness accuracy) moved a retired plan into `docs/archive/` and
  was pushed after the hook and path checks but without pytest, because it changed no
  code. Its CI failed on all four test legs: a docs contract walks every `docs/**/*.md`
  for stale mFSD wording, and the old plan carried it. Contract tests read `docs/`,
  `nextflow/` and the README, so no change is "docs only" to the suite.
  The plan moved to `.agents/archive/` (no test edited).
- **Promotion target:** `.agents/skills/qa-check/SKILL.md`: the full gate runs for a
  docs-, skills- or Nextflow-only change too, and a file moved into `docs/` is checked
  against the tests that walk it. (Done.)

## [LRN-20261021-001] rule-body | every option offered for a decision is measured, not estimated
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** On #294 three guard options were offered with only the current rule
  measured; the others' failure rates were reasoned from thresholds. The operator asked
  whether each was measured. Prototyping and running all of them on the same sets and
  HG002 changed the picture: the simulators wrote every fragment in one pair
  orientation (so the orientation option failed for an artifact of the data), option
  2's low-depth risk had no test set, and the tolerance amendment, measured, reversed my
  recommendation (then four existing contracts caught its false ALT on another allele's
  reads, fixed by forgiving one-off errors only).
- **Promotion target:** `.agents/memory/decisions-worked-examples.md`: every option
  offered is prototyped and measured on the same sets and truth arm, the sets mirroring
  a real library (both orientations, errors, low depth and depth). (Done.)

## [LRN-20261020-001] rule-body | a run's waiter watches its process, not a marker in its log
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** On #294 two background waiters polled a run's log for its closing
  marker. A restart reused the log name and overwrote it, and a killed run never wrote
  its marker, so both waiters hung until killed by hand; the operator asked whether the
  tasks were still running and said the waiting was launched wrongly.
  The operator then asked why only the waiter showed: the runs were launched detached
  on a stale premise (a 30-minute tracked-task limit; tracked tasks ran 6+ hours).
- **Promotion target:** `.agents/memory/long-runs-visible.md` "How to apply": launch the
  run itself as the tracked background task (tee a log); `nohup` only to outlive the
  session, then a waiter looping on `kill -0 <runner pid>`; one runner per log name; a
  run and its waiter stop together; clear partial outputs before a restart. (Done.)

## [LRN-20261014-001] rule-body | the best route holds in the analysis and in the implementation
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** On #277 the locus-aware rule looked finished after the first synthetic
  competitor arm. Five more synthetic shapes and a count-neutral diagnostic each found
  something the real arms and HG002 could not: a germline SNP in the grown repeat taken
  for a rival; ALT carriers aligned around the anchor counted as REF molecules (an
  alignment-trusting tally, which also made me call a real rival false); a rival's
  reference base read as a change REF carries; the MNP check judging the same question
  without the rule; and a pre-existing silent loss of REF reads to a germline SNP inside
  the grown event (RD -33..-50% on 4 real rows). The operator asked that both the
  analysis before implementation and the implementation itself take the best and most
  accurate route, past the current standard, testing and implementing wherever we can.
- **Promotion target:** `.agents/rules/code-quality.md` — `DONE:` two Step Review items:
  analysis stresses each candidate with synthetic shapes (guarded risk, realistic mix,
  phantom, natural quality), scores what matters, checks the oracle trusts no alignment,
  and sizes a suspected defect with a count-neutral diagnostic; implementation turns
  every stress shape into a red-first test, lands the rule in every path judging the
  same question, unit-tests each component, and reruns the measurements.
- **Related:** [[best-not-least-bad]], [[better-than-standard]], [[genotyper-not-caller]].

## [LRN-20261013-001] rule-body | a decision is not settled by picking the least-bad measured option
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** On #277 every candidate rule had a cost, and I recommended the one that
  cost least on HG002 (the first cut) with its repeat-position defect pinned as a known
  limit. The operator agreed but corrected the approach: the aim is the best and most
  accurate rule, not the least-bad. The comparison was incomplete: HG002 germline hets
  cannot show the risk the rule guards (another allele of the same length at the locus),
  so no candidate was measured where it matters.
- **Promotion target:** `.agents/rules/code-quality.md` — `DONE:` a Step Review
  Discipline item: when every option costs something, name what the measurements cannot
  see, measure it, and look for a rule better on every axis before recommending; only
  then recommend with the residual named.
- **Related:** [[best-not-least-bad]], [[better-than-standard]].

## [LRN-20261012-002] rule-body | a new test file's name is checked before it is written
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** I wrote new C39 placement tests to `tests/test_carrier_placement_contract.py`
  through a Python heredoc (`open(p, "w")`), which bypasses the Write tool's
  read-before-overwrite guard. The name was already the C37 (#245) contract; its five
  tests vanished and the suite still passed (1445). Caught only because the count fell
  short of the reviewer's 1453 plus the four new tests.
- **Promotion target:** `.agents/rules/code-quality.md` — `DONE:` under Test Invariants:
  a new test file starts with `git ls-files tests/<name>.py` empty, or is created with
  the Write tool; a shell write once replaced a contract unseen.
- **Related:** check the test count against the last known count after adding tests.

## [LRN-20261012-001] rule-trigger | count-changing ticket work skipped the every-data-type acceptance rule
- **Status:** resolved (rule promoted)
- **Cause:** rule-trigger
- **Summary:** C40's choices (k, read-on, the placement candidate) were measured on GIAB
  HG002 (one genome, BWA-MEM, germline hets, mostly tandem-repeat deletions) and on panel
  sign-out (not truth), with the matched-normal arm skipped and no ACCESS, WES or RNA
  run. The rule that acceptance covers every assay (LRN-20260923-001) lives in
  add-feature step 6, which plan-driven ticket work never loads, so it stayed silent.
  The operator also widened the tool survey to tools of other applications that judge
  reads alike, and noted DNA and RNA may need different rules.
- **Promotion target:** `.agents/rules/code-quality.md` (AFTER implementing each step) —
  `DONE:` a count-changing step needs real-data acceptance on every data type it reaches
  and more than one truth source; `.agents/skills/add-feature/SKILL.md` step 6 — `DONE:`
  adds WES/WGS, RNA, normals as hom-REF controls, spike-ins, DNA/RNA judged separately;
  memory `survey-several-tools` — `DONE:` SV/graph genotypers, realigners, split-read and
  RNA indel tools.
- **Related:** [[test-across-data-types]], [[survey-several-tools]], LRN-20260923-001.

## [LRN-20261009-001] rule-trigger | an external review is checked against tickets, REJECTED.md and measurements before tickets
- **Status:** resolved (rule promoted)
- **Cause:** rule-trigger
- **Summary:** asked to "verify what makes sense" in an outside review, I compared it
  with current behaviour only and offered to file tickets. The existing rule (grep
  `REJECTED.md` before re-proposing) never fired, because the proposals were not mine,
  and I neither mapped them to the open tickets nor measured their value. The operator
  asked to measure first; measuring showed most of the value was already covered
  (existing flags caught 39 of 43 low-explainability rows) and found an open ticket
  (C33 #214) already fixed.
- **Promotion target:** `AGENTS.md` (Learning loop) — `DONE:` one line: an external
  review or proposal is mapped to behaviour, open tickets, `REJECTED.md` and prior
  probes, and measured on our data, before any ticket is filed.
- **Related:** [[external-review-measure-first]] memory.

## [LRN-20261008-001] rule-body | a VAF claim follows the documented definition
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** asked how to read the AR 9 bp GGC-repeat deletion (6.6.0: RD 17, AD 15,
  DP 191), I said to compute VAF as ALT/depth (7.9%). The operator asked whether not
  using depth was right. It is: gbcms's `vaf` is AD/(RD+AD) (documented in
  docs/index.md and architecture.md); DP includes the reads that end inside the
  tract and cannot be judged, so AD/DP is biased low. The informative VAF is 48% in
  the tumor and 0% in the matched normal: a high-fraction somatic event the clinical
  6.6% understated. Same failure family as LRN-20260923-002: a claim stated without
  checking the documented definition.
- **Promotion target:** `DONE:` new memory `.agents/memory/vaf-informative-reads.md`
  and its MEMORY.md index line (the rule in the line itself, so it is in view each
  session). The docs were already right; nothing in the code to change. 6.7.0
  candidates noted to the operator: the two meanings of `total_count` (per-sample DP
  vs merged ref+alt) in the output docs, spanning bias in repeats (#253), and a
  diagnostic when RD+AD is a small share of DP.

---

## [LRN-20261005-002] rule-trigger | a review-agent prompt names the venv maturin installs into
- **Status:** resolved (rule promoted)
- **Cause:** rule-trigger
- **Summary:** the group 5 review agent ran `.rv/bin/maturin develop` from its
  worktree, which sits under `.claude/worktrees/` inside the main checkout; with
  no `VIRTUAL_ENV`, maturin installed into the main checkout's `.venv`, whose
  `gbcms.pth` then pointed at the worktree. My prompt gave that command. The same
  failure was recorded on 2026-10-03 in the memory file, but its index line (the
  part loaded every session) did not carry the rule, so it was not in view when I
  wrote the prompt. Repaired by `maturin develop` from the main root.
- **Promotion target:** `DONE:` `.agents/memory/MEMORY.md` — the
  worktree-venv index line now says every review-agent prompt sets
  `VIRTUAL_ENV=<scratch venv>` for `maturin develop`.
- **Related:** [[worktree-tests-need-isolated-venv]], [[maturin-develop-repo-root-only]].

## [LRN-20261005-001] rule-body | a derived output's behaviour is checked where each writer computes it
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** for group 3 I wrote in the docs, the changelog, the PR and an
  issue comment that R5 leaves `vaf` unchanged "because vaf is alt / total",
  read from `gbcms merge`'s formula without its inputs (merge's total is
  REF + ALT). Every writer computes vaf as alt / (REF + ALT), so R5 raises it at
  a splice-crossing deletion with carriers (no measured row moved: the affected
  rows had none). Found while mapping group 4; corrected in `3d44d240`. A
  follow-on claim that merge's VAF differs from the writers' was wrong for the
  same reason and was withdrawn before anything was filed.
- **Promotion target:** `DONE:` `.agents/memory/holistic-effects-map.md` — a
  derived output is checked at every place it is computed before stating how a
  change moves it.
- **Related:** [[holistic-effects-map]].

## [LRN-20261003-002] rule-body | runs over a network mount go one at a time
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** to speed up group 3's prototype acceptance I ran five builds'
  comparisons and three base runs in parallel over the same SFTP mount. Each
  probe sample crawled (about 25 minutes) and one read hit an NFS timeout. The
  operator: access one file on one mount at a time; serial is more responsive
  than several connections. The run guidance covered visibility and wedged
  mounts but not concurrency.
- **Promotion target:** `DONE:` `.agents/memory/long-runs-visible.md` (runs over a
  mount go one at a time, grouped by file) and the local mounts memory (the
  concurrency rule and the measurement); the group 3 harness now runs serially
  grouped by BAM (`serial_rest.py`).
- **Related:** [[long-runs-visible]].

## [LRN-20261003-001] rule-body | a practice survey includes the literature and the domain's own tools
- **Status:** resolved (rule promoted)
- **Cause:** rule-body
- **Summary:** the survey rule named general callers and genotypers only; for an
  RNA group the relevant practice lives in RNA tools (SplitNCigarReads, ASE
  counters, RNA variant and editing callers) and in method papers. The operator
  asked to "also look at literature and other tools of what people have done".
- **Promotion target:** `docs/reference/read-judgment.md` "How a change to read
  judgment is made", step 1 — `DONE:` adds the domain's tools and the published
  literature. Memory `survey-several-tools` updated.
- **Related:** [[survey-several-tools]].

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
  memory `count-the-given-allele.md`; `.agents/archive/CYCLE_6.6.0_PLAN.md` C1 drops the exception,
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
