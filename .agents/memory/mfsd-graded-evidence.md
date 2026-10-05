---
name: mfsd-graded-evidence
description: mFSD fragment-size output is graded evidence toward CH or toward somatic, never a hard origin call; don't over-interpret it.
metadata:
  type: feedback
---

Report mFSD (ALT vs REF fragment sizes) as *increased confidence toward CH* or
*increased confidence toward somatic (tumor)*, with "no evidence either way" as the
honest answer when there are too few ALT fragments. Never a hard TUMOR-LIKE / CH-LIKE
origin call. Operator, 2026-10-05, during group 6 (S2 #154).

**Why:** per-variant size evidence is weak at the ALT counts ACCESS has (shift
detected in ~8% of draws at 5 ALT fragments, ~20% at 30), and the field agrees
(Munzur 2025: tumor and CH "not individually discriminable"; Marass/Tsui 2020: AUC
0.81 only at ≥20 reads). The matched buffy coat, not fragment size, decides CH.

**How to apply:** prefer a signed evidence measure with a neutral zone (e.g. a
likelihood ratio of tumor-shifted vs REF-like sizes, which stays near 1 at low n)
over thresholds on a non-significant test; word labels as "leans CH / leans
somatic / inconclusive"; keep the evidence subordinate to matched-normal data.
The fragment test is plasma-only (ALT vs REF fragments at the same locus, same
sample); the buffy coat is the answer key for validation, never an input, and its
fragment sizes are never used (sheared gDNA, not cfDNA). A matched-normal BAM input
is a possible later feature: depth-aware counts (0 ALT counts against CH only when
the normal could have seen it; ACCESS buffy coats run at ~half plasma depth), shown
beside the fragment evidence, not merged into it (operator, 2026-10-05).
Related: [[fragmentomics-reference]], [[genotyper-not-caller]], [[better-than-standard]].
