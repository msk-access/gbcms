# gbcms — durable memory index

One line per memory. Full content lives in the linked file. Keep this index tight.

## How we work (feedback)
- [Harness stays small](harness-stay-small.md) — curate skills/rules ruthlessly; every line is paid for each turn.
- [No ticket labels in code](no-ticket-labels-in-code.md) — comments/logs explain what/why/how, never CR-/HI-/ME-/P4c labels.

## References
- [Claudelicious harness](claudelicious-reference.md) — the upstream pattern this project's harness follows.

## Project facts (from the 2026-06-26 code review)
- [Bin fetch-end must cover the anchor variant](bin-anchor-coverage.md) — CR-1 invariant for binned↔legacy parity.
- [WFA fast-path shares the base-quality gate](wfa-bq-gate-contract.md) — CR-2 cross-backend quality contract.
- [Engine should be output-aware](engine-output-aware.md) — plumb intent across FFI vs compute-then-discard.
- [Tolerant large-deletion match is deliberate](tolerant-deletion-deliberate.md) — CR-4; don't regress the sensitivity fix.
- [Nextflow defaults diverge from CLI](nextflow-cli-default-divergence.md) — keep nextflow.config in sync with CLI defaults.

## Tooling / build
- [pyproject is the dep source of truth](deps-pyproject-source-of-truth.md) — CI/Docker bypass lockfiles; declare every directly-imported package.
- [Black version skew (venv vs CI)](black-version-skew-venv-vs-ci.md) — lint with mambaforge black 26.5.1 (=CI), not the venv's stale 25.9.0; they disagree and fake "drift".
- [Legacy count_bam parity oracle](legacy-parity-oracle.md) — feature-gated (`legacy-parity`, default on); shipped wheel omits it; mirror binned-path changes in both or parity breaks.

## Validation / testing
- [pysam validation oracle](pysam-validation-oracle.md) — use fetch()+get_reference_positions (not pileup) to cross-check gbcms counts; RD/AD match exact, DP includes neither.
- [Siblings break binned↔legacy parity](siblings-break-binned-legacy-parity.md) — never pass sibling_variants to count_both; parity only holds without siblings.

## User (private, local-only — not committed)
- `user-*.md` memories (e.g. who the operator is, personal defaults) live in this
  recall dir but are gitignored, never published. See `.gitignore`.
