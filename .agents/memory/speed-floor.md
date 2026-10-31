---
name: speed-floor
description: "Speed floor: gbcms runs faster than the original C++ GetBaseCountsMultiSample, as a minimum; production-size benchmark before each release, for the default path."
metadata:
  type: project
---

gbcms must run faster than the original C++ GetBaseCountsMultiSample in general, as a minimum
(operator, 2026-10-30: "we should make sure that we are running faster then original C++
implementation in general as a minimum").

**Why:** gbcms replaces the C++ tool in MSK pipelines. A replacement that is slower costs every
pipeline run, whatever it gains in counts.

**How to apply:**
- Measure against the production C++ binary 1.2.5 on HPC (run in place; the operator allowed it)
  on the same BAMs, variant list, reference, filters, counts and threads. Use production sizes, not
  the harness's small lists: the C++ tool spends about 8 s loading the reference, which hides the
  per-row cost on small lists. The kit is `kit/run_prod.sbatch` with `kit/manifests/pf71.*`
  (IMPACT and ACCESS, 1 BAM × 1k/5k/20k rows and a 9–10-BAM fillout, whole nodes, 3 repeats).
- Run that benchmark before each release, for whichever path is the default. Speed work goes to
  the path that is or will be the default; a speed-up only the outgoing path uses is out of scope
  (operator, 2026-10-31).
- Measured 2026-10-30 (develop fe34f443): the rules path meets the floor in all 8 cases (1.4–2.4×
  faster; ACCESS 1 BAM × 20,000 rows about equal); the sequence judge does not from 5,000 rows up
  (1.3–5.2× slower). The judge meets the floor before it becomes the default (6.8.0, #331; #335).
- Report times with the rows (variants) and reads per data set, so a ratio such as "2× slower" has
  concrete numbers behind it.

Related: [[universal-haplotype-judge]], [[working-agreements]].
