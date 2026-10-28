---
name: best-not-least-bad
description: When every measured option costs something, don't settle for the least-bad one; find what the measurements can't see and look for a rule better on every axis
metadata:
  type: feedback
---

In a hard read-judgment case the goal is the best and most accurate rule, not the
option that costs least on the data at hand. When every candidate has a cost, the
recommendation is not "keep the cheapest and pin its defect": first name what the
measurements cannot see (the risk the rule exists to guard, a truth source that would
show it), measure that, and search for a rule that is better on every axis. Only when
none exists, recommend with the residual named and the measurement that would settle it.

**Why:** operator, 2026-10-13, on #277 (an unreadable ALT core): I recommended keeping
the first cut and pinning its repeat-position defect because the stricter rules cost a
few HG002 rows. The operator agreed but said the point is "to find the best and the
most accurate solution and not given to mediocrity". HG002 germline hets cannot show
what that rule guards against (another allele of the same length at the locus), so the
comparison was incomplete.

**How to apply:** in every decision table add a row for "what this measurement cannot
see" and how to measure it; prefer building that measurement (spike-in, a hotspot
audit, a synthetic competitor) over choosing among under-measured options. Related:
[[better-than-standard]], [[test-across-data-types]], [[decisions-worked-examples]],
[[fix-the-class-not-the-row]].

**Both phases (operator, 2026-10-14):** the best and most accurate route, past the
current standard, holds in the analysis before implementation *and* in the
implementation itself, "testing and implementing as and where we can". On #277 each
new synthetic shape (germline SNP in the grown repeat, realistic low-VAF somatic on a
germline het, phantom with the allele absent, natural quality) found a flaw the real
arms could not, an alignment-trusting tally made me misjudge a real rival, and a
count-neutral diagnostic sized a pre-existing silent RD loss on every arm. So: stress
every candidate with such shapes and score VAF/phantom share, keep oracles free of
alignment trust, size suspected defects with diagnostics; then turn every shape into a
red-first test, land the rule in every path judging the same question, unit-test each
component, and rerun the measurements on the final build.
