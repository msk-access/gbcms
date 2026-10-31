---
name: speed-floor
description: "Speed floor: gbcms runs faster than the original C++ GetBaseCountsMultiSample, as a minimum; check it before a path becomes a default."
metadata:
  type: project
---

gbcms must run faster than the original C++ GetBaseCountsMultiSample in general, as a minimum
(operator, 2026-10-30: "we should make sure that we are running faster then original C++
implementation in general as a minimum").

**Why:** gbcms replaces the C++ tool in MSK pipelines. A replacement that is slower costs every
pipeline run, whatever it gains in counts.

**How to apply:**
- Measure against the C++ tool on the same BAMs, variant list, reference, filters, counts and
  threads, on every data type (panel, ACCESS, WES, HG002). The source is
  `~/Documents/Github/GetBaseCountsMultiSample`; give the flag mapping with the numbers.
- The default rules path meets the floor for every release. `--read-judge sequence` meets it before
  it becomes the default (6.8.0, #331; speed work #335).
- Report times with the rows (variants) and reads per data set, so a ratio such as "2× slower" has
  concrete numbers behind it.

Related: [[universal-haplotype-judge]], [[working-agreements]].
