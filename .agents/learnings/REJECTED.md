# REJECTED

The rejected-edit buffer (negative feedback). When a proposed change is vetoed
("no, leave that", "we decided against that"), log it here so no future session
re-proposes the dead end. Before promoting any learning, grep this file for the
target + topic; if a matching veto exists, surface "vetoed YYYY-MM-DD because X"
and ask whether something changed, rather than re-proposing.

Newest at the top.

---

## [REJ-20260627-001] HI-5: strip N bases from the read before the PairHMM (to make N "LLR-neutral near indels")
- **Target:** `rust/src/counting/pairhmm.rs` (`classify_by_marginalized_pairhmm` / `BQEmission`)
- **Proposed:** the M4 scope suggested removing N positions from the read (or otherwise
  tweaking N emission) so an N at an indel discriminating locus contributes a strictly
  zero LLR.
- **Reason vetoed:** an empirical probe disproved the premise. With the existing 0.25
  N emission **and** the marginalized (forward-sum) PairHMM, N at a deletion locus
  already yields LLR ≈ **0.41** — far below the 2.3 decision threshold, so it stays
  ambiguous (and N at a SNP locus is exactly 0). **Stripping the N makes the read
  mimic the deletion haplotype → LLR ≈ 8.3 → misclassifies as ALT.** So the proposed
  fix regresses, not improves. HI-5 is already handled correctly; the action taken was
  a regression test (`test_n_base_neutral_at_snp_and_indel`) that guards the behavior,
  not a code change. (HI-4 clamp and ME-5 `prob_emit_y=1.0` were the real fixes in the
  bundle.)
- **Date:** 2026-06-27

## [REJ-20260626-001] Adopt the infra-bound harness layers now (agents, fleet, second brain, multi-provider, crons)
- **Target:** the harness as a whole (`docs/harness/00-plan.md`)
- **Proposed:** build out claudelicious docs 07–19 alongside the spine.
- **Reason vetoed:** these are **infra-bound** — their value scales with running
  services / defined long-horizon jobs, not codebase maturity. Building them with
  nothing to hold is the "half-used scripts" anti-pattern. Promote one only when a
  concrete trigger appears (e.g. an always-on agent *iff* we decide to grind the
  M1–M6 backlog autonomously). Context-bound layers (memory, learning loop, rules)
  were adopted instead.
- **Date:** 2026-06-26
