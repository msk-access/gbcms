---
name: mfsd-graded-evidence
description: "mFSD fragment-size output is graded evidence, never a hard origin call: the report can lean somatic or say it has no size evidence (nothing leans clonal hematopoiesis), the matched buffy coat decides CH, and the fragment test is plasma-only. Includes the Tsui et al. fragmentomics prior."
metadata:
  type: feedback
---

Report mFSD (ALT against REF fragment sizes) as graded evidence, with "no evidence either way" as the
honest answer when there are too few ALT fragments. Never a hard TUMOR-LIKE or CH-LIKE origin call
(operator, 2026-10-05, group 6, S2 #154). The former CH-LIKE class was removed in 6.6.0.

**As built** (`src/gbcms/report/mfsd_report.py`, `_classify_origin`; `docs/reference/mfsd-report.md`):
LEANS-SOMATIC, NO-SIZE-EVIDENCE or INSUFFICIENT. Nothing leans CH: on labeled ACCESS plasma the
somatic lean never fired on a white-cell (CH) variant, and a non-significant test is also the usual
result for tumour variants at these fragment counts (about 30% look REF-like below 50 fragments), so
it is reported as no size evidence, never as evidence for CH. Gene membership in the CH gene set is a
note, never a gate. The operator's 2026-10-05 wording allowed "increased confidence toward CH"; the code's
docstring gives the measurement above as the reason the report does not use it. Ask the operator
before you add a CH lean.

**Why:** per-variant size evidence is weak at the ALT counts ACCESS has (a shift is detected in about
8% of draws at 5 ALT fragments, about 20% at 30), and the field agrees (Munzur 2025: tumour and CH
"not individually discriminable"; Marass/Tsui 2020: AUC 0.81 only at 20 or more reads). The matched
buffy coat, not fragment size, decides CH.

**How to apply:**
- Prefer a signed evidence measure with a neutral zone (a likelihood ratio of tumour-shifted against
  REF-like sizes, which stays near 1 at low n) over thresholds on a non-significant test. Word labels
  as "leans somatic" or "inconclusive".
- Keep the evidence subordinate to matched-normal data. The fragment test is plasma-only (ALT against
  REF fragments at the same locus, same sample). The buffy coat is the answer key for validation,
  never an input, and its fragment sizes are never used (sheared gDNA, not cfDNA).
- A matched-normal BAM input is a possible later feature (operator, 2026-10-05): depth-aware counts
  (0 ALT counts against CH only when the normal could have seen it; ACCESS buffy coats run at about
  half plasma depth), shown beside the fragment evidence, not merged into it.

## Prior work

"Cell-Free DNA Fragmentomics: The New 'Omics' on the Block", Chiu, Heitzer, Lo, Mouliere, Tsui;
Clin Chem 2020;66(12):1480 (https://academic.oup.com/clinchem/article/66/12/1480/5998370). A Q&A
perspective, not a validation study: fragment size, ends and nucleosome footprints carry
tissue-of-origin signal; tumour, placental and transplant cfDNA runs about 25 bp shorter than
haematopoietic cfDNA; a route to telling tumour from clonal haematopoiesis. The operator noted
(2026-10-01) that this fragment-size testing was already being done at MSK with Dana Tsui. It is the
prior for mFSD (ALT against REF fragment-size comparison, S1 #153 per-fragment mean LLR, S2 #154
CH-LIKE must not rest on an underpowered KS; thresholds from a real-data power curve).
