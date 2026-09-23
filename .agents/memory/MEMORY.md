# gbcms — durable memory index

One line per memory. Full content lives in the linked file. Keep this index tight.

## How we work (feedback)
- [Harness stays small](harness-stay-small.md) — curate skills/rules ruthlessly; every line is paid for each turn.
- [No ticket labels in code](no-ticket-labels-in-code.md) — comments/logs explain what/why/how, never CR-/HI-/ME-/P4c labels.
- [Identity band = annotation tolerance](identity-band-annotation-tolerance.md) — the ≥90% insert band exists for imperfect ALT representations; BQ masking can't replace it — test both policies. No new columns/flags: new signal goes to logs/diagnostics/validation tooling.
- [Commit before review workflows](commit-before-review-workflows.md) — same-checkout review agents may `git stash`; commit first, check git state after.

## References
- [Claudelicious harness](claudelicious-reference.md) — the upstream pattern this project's harness follows.

## Project facts (from the 2026-06-26 code review)
- [Bin fetch-end must cover the anchor variant](bin-anchor-coverage.md) — CR-1 invariant for binned↔legacy parity.
- [WFA fast-path shares the base-quality gate](wfa-bq-gate-contract.md) — CR-2 cross-backend quality contract.
- [Engine should be output-aware](engine-output-aware.md) — plumb intent across FFI vs compute-then-discard.
- [Tolerant large-deletion match — superseded](tolerant-deletion-deliberate.md) — issue #91: real large dels are exact-length; wrong-length pure indels → partial_alt; delins stay Phase-3. The 50bp gate is an artifact-SIZE prior (artifacts are small; a ≥50bp op is real), not an event-rarity claim.
- [Nextflow defaults diverge from CLI](nextflow-cli-default-divergence.md) — keep nextflow.config in sync with CLI defaults.

## Tooling / build
- [pyproject is the dep source of truth](deps-pyproject-source-of-truth.md) — CI/Docker bypass lockfiles; declare every directly-imported package.
- [Lint-tool version skew (local vs CI)](black-version-skew-venv-vs-ci.md) — venv black 25.9 lags CI 26.5 (rustc caught up to 1.96 on 2026-09-23); run CI's versions before trusting a clean/drift result.
- [Worktree tests need an isolated venv](worktree-tests-need-isolated-venv.md) — the shared .venv imports the MAIN checkout's gbcms; in a worktree build a scratch venv + maturin develop or tests exercise the wrong code.
- [maturin develop: repo root only](maturin-develop-repo-root-only.md) — `-m rust/Cargo.toml` bypasses [tool.maturin] and leaves a stale src/gbcms/_rs.so shadowing every rebuild.
- [Legacy count_bam parity oracle](legacy-parity-oracle.md) — feature-gated (`legacy-parity`, default on); shipped wheel omits it; mirror binned-path changes in both or parity breaks.

## Validation / testing
- [QC-fail flag absent from MSK data](qcfail-flag-absent-msk-data.md) — no pipeline stage sets 0x200; filter verified correct but inert in practice (contract test pins it).
- [pysam validation oracle](pysam-validation-oracle.md) — use fetch()+get_reference_positions (not pileup) to cross-check gbcms counts; RD/AD match exact, DP includes neither.
- [Siblings break binned↔legacy parity](siblings-break-binned-legacy-parity.md) — never pass sibling_variants to count_both; parity only holds without siblings.

## User (private, local-only — not committed)
- `user-*.md` memories (e.g. who the operator is, personal defaults) live in this
  recall dir but are gitignored, never published. See `.gitignore`.
