---
name: better-than-standard
description: Community practice is a floor, not a ceiling — where no tool sets a standard, or measurement shows a rule better than the field's, adopt the better rule and record why.
metadata:
  type: feedback
---

The survey of other tools informs a decision; it does not bind it. Where the
tools disagree, document nothing, or none handles the case, and where measuring
shows a rule that beats what they do, propose the better rule (with the
measurement) rather than matching the field.

**Why:** operator, 2026-10-04, during group 3: "if there are no standards or we
can do better than standards after measuring we should do that". Group 3's C32
was such a case: no surveyed tool confirms a junction-adjacent base against the
next exon (GATK's overhang fixer reads only the intron, JACUSA and SNPiR mask
fixed distances), and judging spliced reads on the spliced reference cut
spurious ALT 19 → 4 on the probes.

**How to apply:** in each decision table, show the field's practice and, beside
it, the rule the measurements support; when they differ, recommend the measured
rule and say plainly that it departs from practice and why. Still measure the
cost (real carriers lost, truth rows moved), and still survey first, so the
departure is informed. Related: [[survey-several-tools]], [[genotyper-not-caller]].
