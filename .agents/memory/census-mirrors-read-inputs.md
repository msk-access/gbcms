---
name: census-mirrors-read-inputs
description: The read census mirrors decided read-input rules (fragment end, filters), so it cannot validate them; check such a rule against an oracle it does not share.
metadata:
  type: feedback
---

The read census (`tests/census.py`) is the oracle for read judgment, but for read
inputs (what a read contributes: its fragment end, filters) it is changed in step
with the engine, so it shares the rule and agrees with a wrong one. The first C17
build clipped bases past the TLEN boundary; the census agreed, yet 13 genuine
FLT3-ITD ALT reads were lost (TLEN, a reference distance, leaves out inserted
bases). Found only by tracing the rows that moved away from the census and
checking the clipped bases themselves.

**Why:** an oracle that implements the rule under test cannot falsify it.

**How to apply:** for a read-input rule, validate with evidence independent of the
rule: the bases in question (adapter motif at the boundary, mismatch rate to the
reference), the mate's actual alignment, and per-read traces at rows where the
engine and census move differently, ITD and indel rows first. Related:
[[genotyper-not-caller]], [[holistic-effects-map]].
