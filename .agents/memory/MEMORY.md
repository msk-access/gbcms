# gbcms — durable memory index

One line per memory. Full content lives in the linked file. Keep this index tight.

## How we work (feedback)
- [Harness stays small](harness-stay-small.md) — curate skills/rules ruthlessly; every line is paid for each turn.
- [No ticket labels in code](no-ticket-labels-in-code.md) — comments/logs explain what/why/how, never CR-/HI-/ME-/P4c labels.
- [Identity band = annotation tolerance](identity-band-annotation-tolerance.md) — the ≥90% insert band exists for imperfect ALT representations; BQ masking can't replace it — test both policies. No new columns/flags: new signal goes to logs/diagnostics/validation tooling.
- [Commit before review workflows](commit-before-review-workflows.md) — same-checkout review agents may `git stash`; concurrent sessions write .agents/memory — commit first, stage explicit paths only.
- [Genotyper, not caller](genotyper-not-caller.md) — gbcms counts reads whose own bases carry the given allele; alignments can be wrong, so judge bases not placement; validate per read across data types.
- [Count the given allele](count-the-given-allele.md) — input is taken as correct; never deconvolute a wrong input into counts; explain via status-reason/diagnostic columns; only opt-in rescues differ.
- [Use GitHub sub-issues](github-sub-issues.md) — group related work as parent + sub-issues (cycle tracker, umbrella items), not combined issues or checklists.
- [Holistic effects map](holistic-effects-map.md) — before any fix/decision, map every place it lands (both counting paths, per-transcript/ASJD/mFSD/observations, rescue/clusters/merge, writers/flags, docs/tests/Nextflow).
- [Test changes need operator notice](test-changes-need-operator-notice.md) — never edit a failing test's expectation without first telling the operator what it asserts, why it's wrong, and the evidence; wait.
- [Survey several tools](survey-several-tools.md) — community practice means GATK, samtools/bcftools, fgbio, VarDict, Strelka2, freebayes, bam-readcount, LoFreq, GetBaseCounts, plus the literature and domain tools (RNA: SplitNCigarReads, STAR, ASE counters, RNA callers); say where handling isn't documented.

## References
- [Claudelicious harness](claudelicious-reference.md) — the upstream pattern this project's harness follows.
- [Fragmentomics (Tsui et al., MSK)](fragmentomics-reference.md) — the prior for mFSD and CH-LIKE; cite for S1 #153 / S2 #154 thresholds.

## Project facts (from the 2026-06-26 code review)
- [Bin fetch-end must cover the anchor variant](bin-anchor-coverage.md) — CR-1; pinned by the bin property test and binning-invariance tests.
- [WFA fast-path shares the base-quality gate](wfa-bq-gate-contract.md) — CR-2 cross-backend quality contract.
- [Engine should be output-aware](engine-output-aware.md) — plumb intent across FFI vs compute-then-discard.
- [Tolerant large-deletion match — superseded](tolerant-deletion-deliberate.md) — issue #91: real large dels are exact-length; wrong-length pure indels → partial_alt; delins stay Phase-3. The 50bp gate is an artifact-SIZE prior (artifacts are small; a ≥50bp op is real), not an event-rarity claim.
- [Nextflow defaults diverge from CLI](nextflow-cli-default-divergence.md) — keep nextflow.config in sync with CLI defaults.
- [MAPQ-0 loci (PMS2)](mapq0-loci-pms2.md) — pseudogene genes run at --min-mapq 0; keep MAPQ-0 alignments countable; validate read-admission changes at --min-mapq 0 too.

## Tooling / build
- [pyproject is the dep source of truth](deps-pyproject-source-of-truth.md) — CI/Docker bypass lockfiles; declare every directly-imported package.
- [Lint-tool version skew (local vs CI)](black-version-skew-venv-vs-ci.md) — venv black 25.9 lags CI 26.5 (rustc caught up to 1.96 on 2026-09-23); run CI's versions before trusting a clean/drift result.
- [Worktree tests need an isolated venv](worktree-tests-need-isolated-venv.md) — the shared .venv imports the MAIN checkout's gbcms; in a worktree build a scratch venv + maturin develop or tests exercise the wrong code.
- [maturin develop: repo root only](maturin-develop-repo-root-only.md) — `-m rust/Cargo.toml` bypasses [tool.maturin] and leaves a stale src/gbcms/_rs.so shadowing every rebuild.
- [Harness cleanup after merge](harness-cleanup-after-merge.md) — delete the ticket's `src_*` builds (rebuildable from BUILT_FROM; patch prototypes first), gzip traces, `cargo clean` when rust/target bloats.
- [Long runs visible](long-runs-visible.md) — pair every detached nohup harness run with a tracked background waiter so the operator sees it.
- [No `timeout` on this Mac](macos-no-timeout-command.md) — `timeout N cmd && ok || fail` always says fail; check SFTP mounts with a plain read.

## Validation / testing
- [QC-fail flag absent from MSK data](qcfail-flag-absent-msk-data.md) — no pipeline stage sets 0x200; filter verified correct but inert in practice (contract test pins it).
- [--trace output wraps](trace-output-wraps.md) — set COLUMNS=3000 before parsing `read call` trace lines, or every read looks untraced.
- [Census mirrors read inputs](census-mirrors-read-inputs.md) — the census shares read-input rules (fragment end, filters), so validate those against an oracle it doesn't share (the bases themselves, the mate's alignment, ITD rows first).
- [pysam validation oracle](pysam-validation-oracle.md) — use fetch()+get_reference_positions (not pileup) to cross-check gbcms counts; RD/AD match exact, DP includes neither.
- [vcf2maf oracle](vcf2maf-oracle.md) — VCF↔MAF representation is defined as vcf2maf/maf2vcf output; how to run them locally and their quirks (dies on N, empty ALT for REF==ALT).

## User (private, local-only — not committed)
- `user-*.md` memories (e.g. who the operator is, personal defaults) live in this
  recall dir but are gitignored, never published. See `.gitignore`.
