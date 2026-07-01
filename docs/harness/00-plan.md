# gbcms Harness Setup Plan

**Branch:** `feat/claude-harness` · **Date:** 2026-06-26 · **Following:** [BioInfo/claudelicious](https://github.com/BioInfo/claudelicious)

This plan adapts the **claudelicious** harness pattern ("run the model as a system, not a chat box") to the gbcms project. Per the source's own guidance — *"You cannot buy a harness… look at the shape of someone else's, see which parts map to your work"* — we adopt the **shape and discipline**, not a clone.

## Decisions (locked)

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Scope | **Full 6-piece spine** | Rules, skills, hooks, memory, continuity, learning loop — the claudelicious "week-1 spine." Defer only fleet/agents/multi-provider. |
| Layout | **Claude Code-native** | `CLAUDE.md` + `.claude/skills/` + `.claude/settings.json` so the harness actually auto-loads. Keep claudelicious's docs-first philosophy. |
| Following upstream | **Pin as reference** | Record upstream URL + commit in this doc; re-check by hand. No submodule coupling. |

## Guiding discipline (from the article)

1. **Shaped to you** — every rule/skill/memory encodes a gbcms-specific judgment call. Generic content is cut.
2. **It learns** — a correction becomes a durable memory + a work-order at the source, so a mistake doesn't recur.
3. **It stays small** — *"every line is paid for on every turn."* Standing context (`CLAUDE.md`) stays tight; depth lives in on-demand reference docs and curated skills. (Justin cut 134 skills → 48; we curate from day one.)

## Mapping: claudelicious piece → gbcms artifact

| # | Piece (claudelicious doc) | gbcms artifact (native layout) | Source material | Status |
|---|---------------------------|-------------------------------|-----------------|--------|
| 1 | Rules & context (01) | `CLAUDE.md` (small) + `.agents/rules/*` (deep, on-demand) | `.agents/rules/{architecture,code-quality,git-conventions}.md` | exists, not loaded |
| 2 | Skills (02) | `.claude/skills/<name>/SKILL.md` | 5 `.agents/skills/*` + 6 `.agents/workflows/*` | exists, not loaded |
| 3 | Hooks (03) | `.claude/settings.json` hooks | — | new |
| 4 | Memory (04) | `~/.claude/projects/<proj>/memory/` + `MEMORY.md` | code-review findings | dir empty |
| 5 | Continuity (06) | `CONTINUITY.md` (root) | current session state | new |
| 6 | Learning loop (05) | convention in `CLAUDE.md` + memory + work-order | — | new |

## Phase-by-phase manifest

### Phase 0 — Branch & audit ✅
- `git checkout -b feat/claude-harness` (done).
- Audited `.agents/`: 3 rules, 6 workflows, 5 skills. All gbcms-specific and worth keeping.

### Phase 1 — Rules / standing context
- **Create `CLAUDE.md`** (root, loaded every turn): project identity, hybrid Rust/Python boundary in one table, the handful of **load-bearing invariants** surfaced by the code review (binned↔legacy parity, cross-backend BQ-gating contract, "engine must be output-aware," type-stub sync, DP/DPF test invariants), the lint gate, and pointers to deep reference + memory + continuity. **Target: tight.**
- **Keep** `.agents/rules/{architecture,code-quality,git-conventions}.md` as on-demand depth that `CLAUDE.md` links to (read when needed, not every turn).

### Phase 2 — Skills (migrate + curate into `.claude/skills/`)
Convert each to a loadable skill (add `name:` + `description:` frontmatter; make action-oriented). Candidate set:
- From `.agents/skills/`: `counting-engine`, `mfsd-analysis`, `mnp-rescue`, `rna-annotation`, `rust-python-ffi` → **reference skills**.
- From `.agents/workflows/`: `add-feature`, `debug-error`, `qa-check`, `release`, `git-commit`, `code-review` → **procedure skills**.
- **Curate:** `git-commit`/`code-review` overlap with built-in `/code-review` and commit conventions — fold or drop to avoid the "stay small" trap. Final set decided during migration, not pre-committed.

### Phase 3 — Hooks (`.claude/settings.json`)
- **Time injection** (the article's first hook — the model has no clock).
- **One guardrail**: a destructive-op / protected-path guard (e.g., refuse recursive deletes outside scratch; protect `main`).
- Preserve existing `.claude/settings.local.json` permissions.

### Phase 4 — Memory (seed from the code review)
- Stand up `MEMORY.md` index in the project memory dir.
- Seed durable, gbcms-specific facts the review produced, e.g.:
  - Bin fetch-end must cover the anchor variant's full ref span (CR-1 invariant).
  - WFA fast-path must share the PairHMM base-quality gate (CR-2 contract).
  - The Rust engine should be output-aware (`mfsd`/format/`gene_strand` plumbed) — closes a class of compute-then-discard/no-op bugs.
  - Tolerant large-deletion match was a deliberate sensitivity fix (don't regress it; add partial seq-check instead).
- This is where the review stops being a one-off doc and becomes a **compounding asset**.

### Phase 5 — Continuity (`CONTINUITY.md`)
- Capture tactical state: harness setup in progress; `CODE_REVIEW_IMPLEMENTATION_PLAN.md` authored with 45 tickets, milestones M1–M6 pending; current branch.

### Phase 6 — Learning loop
- Document the convention in `CLAUDE.md`: *when corrected → (a) write/refresh a memory file, (b) emit a work-order naming the exact file:line to fix.* Wire it to the memory dir from Phase 4.

### Deferred (claudelicious "skip this week" — docs 07–19)
Session search (mneme/LanceDB), second brain (Qdrant), multi-provider routing, crons, always-on agents (Pulsar), the fleet. Valid later; out of scope for the spine.

## Success criteria (from claudelicious QUICKSTART)
Start a **fresh session** and verify the agent:
1. Reads `CLAUDE.md` and respects its constraints without being told.
2. Recalls a seeded memory fact.
3. Can invoke a migrated skill by outcome.
4. Reads `CONTINUITY.md` and knows where we left off.
5. On correction, writes it down (memory + work-order).

## Cross-tool portability (added after initial build)

The harness is deliberately **not** welded to one tool. `AGENTS.md` is the canonical
standing-context format (Linux Foundation Agentic AI Foundation standard; read by
Antigravity v1.20.3+, Cursor, Codex, Copilot, Windsurf, Zed, Junie, …). Principle:
**one neutral source of truth + thin per-tool adapters.**

| Piece | Canonical (neutral) | Claude adapter | Antigravity / others |
|-------|--------------------|----------------|----------------------|
| Standing context | `AGENTS.md` | `CLAUDE.md` = `@AGENTS.md` import | read `AGENTS.md` natively (`GEMINI.md` = Google override only) |
| Skills | `.agents/skills/<name>/SKILL.md` | `.claude/skills/` → symlink | Antigravity reads `.agents/skills/` |
| Enforcement | `.githooks/pre-commit` + `guard-destructive.sh` | `.claude/settings.json` wires the hook | git hook is tool-agnostic |
| Memory | `.agents/memory/*.md` (in-repo) | Claude recall dir → symlink | any tool ingests the same files |
| Continuity | `CONTINUITY.md` | referenced from `AGENTS.md` | referenced from `AGENTS.md` |
| Learning loop | convention in `AGENTS.md` | inherited | inherited |

Setup on a new clone: `git config core.hooksPath .githooks`. Verify in a fresh
session in each tool (success criteria above).

## Reference pin
- **Upstream:** https://github.com/BioInfo/claudelicious — branch `main`, as referenced 2026-06-26. (Resolve exact commit SHA on first sync.)
- **Article:** "The Harness Is the Moat," Justin Johnson / Run Data Run, 2026-06-25.
- **License note:** claudelicious code is MIT, prose is CC BY-NC 4.0 — we adapt patterns and attribute; we do not vendor prose without attribution.
- **AGENTS.md standard:** stewarded by the Linux Foundation Agentic AI Foundation; 28+ tools.
