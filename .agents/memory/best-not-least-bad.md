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
[[better-than-standard]], [[test-across-data-types]], [[decisions-worked-examples]].
