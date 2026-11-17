---
name: speed-floor
description: "Speed floor: gbcms runs faster than the original C++ GetBaseCountsMultiSample, as a minimum; production-size benchmark before each release, for the default path. Speed is never bought with accuracy: only exact (count-neutral) speed-ups."
metadata:
  type: project
---

gbcms must run faster than the original C++ GetBaseCountsMultiSample in general, as a minimum
(operator, 2026-10-30: "we should make sure that we are running faster then original C++
implementation in general as a minimum").

**Why:** gbcms replaces the C++ tool in MSK pipelines. A replacement that is slower costs every
pipeline run, whatever it gains in counts.

**Speed is never bought with accuracy** (operator, 2026-11-14: "we dont prioritze speed over
accuracy, we need to find solutions where speed is gained without loosing any accuracy"). A speed
gain must give the same output, every field, on every set. A cap, a sample, a skipped read or a
looser test that can change one read's call is not a speed option, not even as a fallback behind a
time bar. Where a decided rule costs time, the rule stays, and the work looks for exact speed-ups:
stop a count once its answer cannot change, compute a shared value once, skip work whose result
cannot reach the output. Show each one count-neutral (identical outputs on every set) before it
merges.

**Memory may rise to buy time, within a limit** (operator, 2026-11-16: "I am ok with memory
increase as long as we can reduce time and maintain accuracy", "Not unlimited increase", "I am ok
with 1.5 to 2x increase"). Report peak memory beside each timed cell, against the baseline before the
work began. Up to 2 times the baseline's peak is accepted when time falls and every output is
identical; name any cell above 1.5 times; above 2 times, bring it to the operator.

**How to apply:**
- Measure against the production C++ binary 1.2.5 on HPC (run in place; the operator allowed it)
  on the same BAMs, variant list, reference, filters, counts and threads. Use production sizes, not
  the harness's small lists: the C++ tool spends about 8 s loading the reference, which hides the
  per-row cost on small lists. The kit is `kit/run_prod.sbatch` with `kit/manifests/pf71.*`
  (IMPACT and ACCESS, 1 BAM × 1k/5k/20k rows and a 9–10-BAM fillout, whole nodes, 3 repeats).
- Time a speed option on every data type its change reaches before it is kept or dropped: RNA with
  a GTF for an RNA or GTF change, mFSD for an mFSD change, a fill-out for a per-BAM change. A part
  that was not timed is never a reason to drop the option (operator, 2026-11-15: "No we should test
  on RNA as well so keep s4").
- Run that benchmark before each release, for whichever path is the default. Speed work goes to
  the path that is or will be the default; a speed-up only the outgoing path uses is out of scope
  (operator, 2026-10-31).
- Measured 2026-10-30 (develop fe34f443): the rules path meets the floor in all 8 cases (1.4–2.4×
  faster; ACCESS 1 BAM × 20,000 rows about equal); the sequence judge does not from 5,000 rows up
  (1.3–5.2× slower). The judge meets the floor before it becomes the default (6.8.0, #331; #335).
- Measured 2026-11-06 for the 6.7.0 release (pf72, develop f7293618, rules path only): met in all
  8 cases, rules wall over C++ 0.39–0.94; ACCESS 1 BAM × 20,000 rows still the closest (0.94).
- Report times with the rows (variants) and reads per data set, so a ratio such as "2× slower" has
  concrete numbers behind it.

Related: [[universal-haplotype-judge]], [[working-agreements]].
