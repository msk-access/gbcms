---
name: external-review-measure-first
description: An external review or proposal is mapped to current behaviour, open tickets, REJECTED.md and prior probes, and its value measured on our data, before any ticket is filed.
metadata:
  type: feedback
---

When the operator brings an outside review (another model's framework, a paper's
proposal), compare it against all four: current behaviour (verified in code/docs, not
from memory), the open tickets and parents of the cycle, `REJECTED.md`, and the
earlier review probes (rerun on the release build and now). Then measure each
proposal's value on our data before offering a ticket: only what measurement shows
relevant becomes one.

**Why:** 2026-10-09, an external "evidence accounting" review: my first pass checked
behaviour only and offered to file tickets; the operator asked to measure value first.
Measuring changed the picture: explainability flagged 43 rows but 39 already carried
a flag; R1/R2 disagreement was 23 of ~195k fragments; no panel insertion was 100+ bp;
hard clips mattered on BWA-only WGS (3.8% of breakpoint records) and not on ABRA2
panels (0%). And rerunning the old probes found C33 #214 already fixed by C39.

**How to apply:** build a map (proposal → status → ticket/REJ id), rerun prior probes
against a release-tag build, measure each open proposal's incremental value (what it
adds beyond existing flags/tickets) on panel, truth and BWA-only arms, then ask which
to file. Related: [[survey-several-tools]], [[better-than-standard]], [[genotyper-not-caller]].
