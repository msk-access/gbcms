# gbcms — standing context (canonical, cross-tool)

> Canonical agent rules for this repo, read by Antigravity, Cursor, Codex, Copilot and (via the
> `CLAUDE.md` import) Claude Code. Every line is read on every turn: keep it to 75 lines. Depth
> lives in the linked docs, read on demand.

**What this is:** GBCMS (Get Base Counts Multi-Sample) — a high-performance, orientation-aware
genotype counting engine. Extracts allele counts and molecular metrics from BAM/CRAM at given
positions: SNP, MNP, insertion, deletion, complex DelIns; DNA (cfDNA) and RNA modes; mFSD.
**Architecture:** hybrid. **Python** = CLI/orchestration/IO/config (Typer + Pydantic + Rich).
**Rust** (`rust/`, PyO3 module `gbcms._rs`) = BAM traversal, read classification, fragment
tracking, alignment (SW / PairHMM / WFA), mFSD stats, GTF annotation, Parquet, Rayon per-bin
parallelism. Module map: `.agents/rules/architecture.md`.

## Load-bearing invariants (don't break these)

1. **Binning invariance.** Bin geometry is performance only: `count_bam_binned` gives identical
   counts (every field; BH q-values excepted when the set of rows changes) under any bin window or
   cap, one variant per bin, or one variant per call. Every bin's fetch holds each member's full
   read window, anchor included, and its guard window. Checks: `architecture.md` §"Binning invariance".
2. **One quality contract across alignment backends.** SW, PairHMM and the WFA fast-path apply the
   *same* base-quality gate. The fast path makes no definitive REF/ALT call on bases the fallback
   would reject.
3. **Engine should be output-aware.** Plumb intent (`mfsd`, output format, `gene_strand`) across the
   FFI boundary; do not compute-then-discard or run silent no-ops in Rust.
4. **0-based internal coords**, 1-based at the VCF/MAF boundary. Convert at the edge.
5. **Type stub tracks the bindings.** `src/gbcms/_rs.pyi` is the single stub for `gbcms._rs`; it
   matches the `#[pyo3(get)]` fields and `#[pyo3(signature)]` params exactly. There is no second stub.
6. **mFSD/RNA columns are gated** — absent when off, never NA-filled.
7. **Count the given allele — a genotyper, not a caller.** Input alleles are taken as correct: never
   credit a read to a row's ALT unless its own bases carry that ALT (alignments can be wrong; judge
   bases, not placement), and never infer what a wrong input "should" be. Explain mismatches via
   `gbcms_status_reason` / `gbcms_diagnostic`. Opt-in rescues (`--rescue-mnp`) are the only exception.

## Counting tests (assert in every counting test)
```python
assert counts.dp  >= counts.rd  + counts.ad      # DP includes 'neither'
assert counts.dpf >= counts.rdf + counts.adf     # DPF includes discarded fragments
assert counts.rd == counts.rd_fwd + counts.rd_rev    # strand consistency
assert counts.ad == counts.ad_fwd + counts.ad_rev
```
**Never edit an existing test's expectation to make it pass** without first telling the operator
what it asserted, why it is wrong, and the evidence (old against new) — then wait. A changed test can
hide a regression or quietly re-decide a decision. A failing test is a decision point (new tests, and
removing the strict-xfail marker a fix targets, are fine).

## Decisions, commits, learning

- **A count-affecting decision** follows `.agents/rules/decision-protocol.md` (read the register to
  its end, effects map, fix the class, several tools, every data type, worked example, measure first;
  an outside review is measured before any ticket).
- **Before every commit** run the whole `qa-check` gate, docs-only changes included (CI's black; read
  every exit code). Conventional commits + Git Flow (`feature/* → develop → main`, never push `main`):
  `.agents/rules/git-conventions.md`, `.agents/rules/code-quality.md`.
- **When I correct you,** do both in the same turn: (1) write it down — a durable note in
  `.agents/memory/` (index `MEMORY.md`), and the work-order in `.agents/learnings/LEARNINGS.md` with
  ONE root cause (rule-body / rule-trigger / rule-permission / environment) and a promotion target
  (exact file + change); external-only failures go to `ERRORS.md`, vetoed ideas to `REJECTED.md` (grep
  it before re-proposing); (2) fix at the source — name the `file:line` and make the edit now.

## Safety (enforced by hooks, not prose)

Determinism floor: `.claude/hooks/` (block-rm-rf, block-push-main, injection-guard, block-env-edits,
protect-cbioportal) + `.githooks/` (pre-commit: data-file and PHI guard, skill checks; pre-push:
gitleaks, main block). Once per clone: `git config core.hooksPath .githooks`; install gitleaks.
Intent, PHI and cBioPortal rules: `.agents/rules/security.md`.

## Pointers

- Deep rules: `.agents/rules/{architecture,code-quality,decision-protocol,git-conventions,security}.md`;
  skills: `.agents/skills/<name>/SKILL.md` (mirrored to `.claude/skills/`), by outcome.
- Where we left off: `CONTINUITY.md`. Durable facts: `.agents/memory/MEMORY.md`. Measurement harness:
  `docs/harness/measurement-harness.md`; harness plan and philosophy: `docs/harness/00-plan.md`.
- Open work: `CONTINUITY.md` "Next" (the 6.7.0 work packages), tracker #196, `PRODUCTION_READINESS.md`.
  Retired plans are in `.agents/archive/`.
