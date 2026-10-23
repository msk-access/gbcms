---
name: decisions-worked-examples
description: Present every count-rule decision with a real worked example drawn as a picture, the community standard, the RNA effect, and a truth check of whether the change moves counts toward truth.
metadata:
  type: feedback
---

When a decision changes how reads are counted, the operator decides from a worked
example, not from aggregate tables alone. Bring, in the first presentation:
- one real row drawn out: the read windows (or read classes) each option makes a read
  hold, with that row's counts and VAF under each option (an inline picture);
- the community standard (several tools, [[survey-several-tools]]);
- the RNA effect, even when RNA has few such rows (say how many);
- a truth check: does the AD/VAF change move toward a known truth (GIAB or germline
  hets at 0.5, [[genotyper-not-caller]]), not only toward the sign-out.
- every option offered measured, not estimated: prototype each one and run it on the
  same sets and truth arm as the recommended one, so the comparison table holds only
  runs; the synthetic sets must mirror a real library (both pair orientations, errors,
  low depth and depth) or an option keyed on that feature is not measured at all.

**Why:** 2026-10-08, 6.7.0 plan decisions 1 and 3 (#253, #254). The operator answered
the first round with: "I need to understand this better with real example and are
there any community standard around it and how does this effect RNA", and "visually if
possible; if the AD change is getting us closer to truth then it's OK". They decided
once the worked example (a 35 bp delins, rule by rule) and the HG002 truth were shown.

2026-10-21 (#294): offered three guard options with only the current one measured,
the others' costs reasoned; the operator asked "is each of the other option measured i
feel these are just estimated". Measuring them reversed one recommendation (the
tolerance amendment) and exposed a simulator flaw (every fragment one pair
orientation) that had made the orientation option look broken.

**How to apply:** build the truth arm (D9 GIAB, D11 germline normals) before the
decision round, and validate the measurement model on truth too: the first HG002 pass
exposed an oracle placement bug that the panel data had hidden.
