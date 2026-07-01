---
name: wfa-bq-gate-contract
description: The WFA fast-path must apply the same base-quality gate as SW/PairHMM; it must not make a definitive call on bases the fallback would reject.
metadata: 
  node_type: memory
  type: project
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
---

`rust/src/counting/wfa_router.rs` aligns the raw read sequence with no base-quality
input — `read_quals`/`min_baseq` are never passed in. SW masks sub-`min_baseq` bases
to N and PairHMM weights by BQ, but the WFA fast-path can return a *confident* REF/ALT
call (edit distance 0) on a low-BQ discriminating base and short-circuit before the
quality-aware fallback runs.

**Why it matters:** Enabling the fast path then *changes counts* for low-BQ reads —
admitting low-quality ALT support exactly in the low-VAF cfDNA regime the BQ model
exists to protect. This is review finding **CR-2** (CRITICAL). Related: the off-target
branch uses global edit distance while PairHMM uses semiglobal, dropping deletion-ALT
reads (HI-3); PairHMM `p_correct` isn't clamped, so Q0 bases give NaN LLR (HI-4).

**How to apply:** One quality contract across SW / PairHMM / WFA. Thread
`read_quals` + `min_baseq` into `wfa_fast_path`; require the same
`usable_count >= MIN_USABLE_BASES` gate before any definitive call, else return None.
