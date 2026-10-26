---
name: describe-labels-in-brackets
description: "When explaining to the operator, describe each item in plain words first and put internal labels (phase A, RJ-13, c46e2, E1, O5, found2) in brackets after it; never lead with a bare label."
metadata:
  type: feedback
---

The operator asked (2026-10-25) not to explain with bare labels: "when explaining dont just
say RJ, Phase A etc, describe it and mention those in bracket". Over a long session the
labels (phases A-E, experiments E1-E4, options O1-O5, rules RJ-n, builds c46c/c48a/c50d,
found2/found3, xs) had piled up until status messages were hard to follow.

**Why:** the operator follows from the app, often after time away; a label is only useful
once its meaning is in front of them.

**How to apply:** in every message to the operator, say what the thing is or does, then the
label in brackets: "the fix round on the graph scorer (phase C)", "clipped bases inside the
fragment are the read's own (RJ-13)", "the build after the fix round (`src_c50d`)". Tables
name rows in words too. Labels without a description are fine only in code, commits and
GitHub records, which keep their own conventions. Related: [[no-ticket-labels-in-code]].
