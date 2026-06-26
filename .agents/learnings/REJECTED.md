# REJECTED

The rejected-edit buffer (negative feedback). When a proposed change is vetoed
("no, leave that", "we decided against that"), log it here so no future session
re-proposes the dead end. Before promoting any learning, grep this file for the
target + topic; if a matching veto exists, surface "vetoed YYYY-MM-DD because X"
and ask whether something changed, rather than re-proposing.

Newest at the top.

---

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
