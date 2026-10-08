---
name: test-changes-need-operator-notice
description: Never edit a failing test to make it pass without first telling the operator what it asserts, why it is wrong, and the evidence — then wait.
metadata:
  type: feedback
---

When a test fails (or a strict xfail turns green) because of a change, do not edit
the test's expectation on my own judgement. First tell the operator: what the test
asserted, why I think that expectation is wrong (the principle it conflicts with),
and the evidence (measured counts, adjudicated reads). Then wait for the decision.
Measuring first to understand the failure is always fine.

**Why:** the operator (2026-10-01) saw tests being changed as part of fixes without
being told (H3: three `test_shifted_indels` expectations moved to depth-only, the
C28 example swapped, a census assertion narrowed). A changed test can hide a
regression or quietly re-decide an earlier decision; the expectation is the
record of what was decided.

**How to apply:** a failing test is a decision point, not a chore. Present it with
the decision items (AskUserQuestion), quoting old vs new expectation. Adding new
tests and flipping a strict xfail that the change was meant to fix (its marker
removal is the point of the red-first commit) do not need asking; editing any
existing assertion does. Related: [[count-the-given-allele]], [[holistic-effects-map]].
