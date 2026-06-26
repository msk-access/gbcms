# gbcms — standing context (canonical, cross-tool)

> Canonical agent rules for this repo. Read by Antigravity, Cursor, Codex, Copilot,
> and (via `CLAUDE.md` import) Claude Code. Every line is read on every turn — keep it
> tight. Depth lives in the linked reference docs and is read on demand, not here.

**What this is:** GBCMS (Get Base Counts Multi-Sample) — a high-performance,
orientation-aware genotype counting engine. Extracts allele counts and molecular
metrics from BAM/CRAM at given genomic positions. Supports SNP, MNP, insertion,
deletion, complex DelIns; DNA (cfDNA) and RNA modes; mFSD fragment-size analysis.

**Architecture:** Hybrid. **Python** = CLI/orchestration/IO/config (Typer + Pydantic + Rich).
**Rust** (`rust/`, PyO3 module `gbcms._rs`) = BAM traversal, read classification,
fragment tracking, alignment (SW / PairHMM / WFA), mFSD stats, GTF annotation,
native Parquet, Rayon per-bin parallelism. Full module map: `.agents/rules/architecture.md`.

## Load-bearing invariants (don't break these)

1. **Binned ↔ legacy parity.** `count_bam_binned` must produce identical counts to
   legacy `count_bam`. Any binning change needs a parity test incl. large deletions
   and complex DelIns. (Bin fetch-end must cover the *anchor* variant's full ref span.)
2. **One quality contract across alignment backends.** SW, PairHMM, and the WFA
   fast-path must apply the *same* base-quality gate. The fast path must not make a
   definitive REF/ALT call on bases the fallback would reject.
3. **Engine should be output-aware.** Prefer plumbing intent (`mfsd`, output format,
   `gene_strand`) across the FFI boundary over computing-then-discarding or running
   silent no-ops in Rust.
4. **0-based internal coords**, 1-based at the VCF/MAF boundary. Convert at the edge.
5. **Type stubs stay synced.** `src/gbcms/_rs.pyi` is authoritative; `src/gbcms_rs.pyi`
   mirrors it. Both must match the `#[pyo3(get)]` fields exactly.
6. **mFSD/RNA columns are gated** — absent when off, never NA-filled.

## Counting test invariants (assert in every counting test)
```python
assert counts.dp  >= counts.rd  + counts.ad     # DP includes 'neither'
assert counts.dpf >= counts.rdf + counts.adf    # DPF includes discarded fragments
assert counts.rd == counts.rd_fwd + counts.rd_rev   # strand consistency
assert counts.ad == counts.ad_fwd + counts.ad_rev
```

## Before every commit (lint gate)
```bash
ruff check src/ tests/ && black --check src/ tests/ && mypy src/
cd rust && cargo clippy --all-targets -- -D warnings && cargo test && cd ..
```
Fix order: ruff → black → mypy → clippy. Conventional commits + Git Flow
(`feature/* → develop → main`, never push `main` directly). Details:
`.agents/rules/git-conventions.md`, `.agents/rules/code-quality.md`.

## Learning loop (when I correct you)
Do both, in the same turn:
1. **Write it down** — add/refresh a fact in `.agents/memory/` (index: `.agents/memory/MEMORY.md`)
   so the correction survives to next session.
2. **Fix at the source** — name the exact `file:line` and make the edit, so the
   mistake can't recur from stale guidance.

## Skills
Project skills live in `.agents/skills/` (cross-tool / Antigravity-native; mirrored to
`.claude/skills/` for Claude Code). Each is `<name>/SKILL.md` with `name` + `description`
frontmatter. Reach for them by outcome.

## Pointers
- Deep rules: `.agents/rules/{architecture,code-quality,git-conventions}.md`
- Skills: `.agents/skills/` · Harness plan & philosophy: `docs/harness/00-plan.md`
- Where we left off: `CONTINUITY.md` · Durable facts: `.agents/memory/MEMORY.md`
- Open work: `CODE_REVIEW_IMPLEMENTATION_PLAN.md` (45 tickets, milestones M1–M6)
