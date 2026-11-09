# CONTINUITY — where we left off

> Tactical state that must survive a closed laptop or a context summary.
> Update the **Now** and **Next** sections as work progresses. When a "Now" section is finished,
> move it to `.agents/archive/CONTINUITY_ARCHIVE.md` (older sections live there).

_Last updated: 2026-11-06_

## Now — 6.7.0, an interim release (operator, 2026-11-06)

**Decision (operator, 2026-11-06):** "it does not make sense to do speed if we are still making
changes to counts". 6.7.0 ships with the count changes merged. Every open 6.7.0 ticket moved to
6.8.0 (tracker #420); only the cycle tracker #196 stays, to close with the release.

**In 6.7.0 (merged to develop):**
- #297, one unit-count rule for long deletions in repeats (#404).
- #264, equal windows for long insertions (#396), and the judge-path follow-up (#413).
- RNA #292 option B: a splice lowers no base's quality (RJ-31), with the partial-span rule RJ-17 (#415).
- Class plan B, "tell what the read reaches" (#419; closed #394 and #398):
  - every listed allele that meets the site is a candidate, with no floor;
  - a read counts REF only where it tells REF from every candidate whose difference it reaches;
  - a long delins has a core (50 or more REF bases replaced);
  - an SNV row beside a listed indel or delins is counted by the judge, also under the rules.
  - HPC rq3 and rq4: every bar passes. rq4 equals rq3 on the rules path in all 1,078 runs; the
    judge path moves only #413's shape (partial to depth only at insertions of 50 bp or more).
- RNA package 3's lean path (#423; closed #343 to #346 and #359): ASJD makes the same sibling
  calls as the main counts; three pins; RJ-32.
- `--rescue-homopolymer` is deprecated (#421) and removed in 6.8.0 (#261): hp1 found that the
  twin wins at 1 of 20 candidate rows.
- The sequence judge, opt-in (`--read-judge sequence`), with its follow-ups.
- Platforms: the RHEL 8 wheel, one abi3 wheel, the macOS arm64 wheel, the multi-arch image and
  Apptainer in CI (#233 to #237); macOS CI on `macos-latest` (#399, #407).

**Release steps (in progress):**
1. The speed floor (the operator's standing rule, before each release): **met in all 8 cases.**
   pf72 on HPC, the release candidate rc670 (develop f7293618), C++ GetBaseCountsMultiSample 1.2.5
   against the rules path on pf71's production tasks, 3 repeats each (7650328; 7650331, 7652235).
   Rules wall time over C++'s: 0.39 to 0.57 for 1 BAM at 1,000 and 5,000 rows; IMPACT 1 × 20,000
   0.51; ACCESS 1 × 20,000 0.94 (as in October); fill-outs IMPACT 10 × 5,000 0.55 and ACCESS
   9 × 5,000 0.74. Against pf71 (develop fe34f443), five ratios are lower and three are slightly
   higher: IMPACT 1 × 5,000 0.462 (0.452), ACCESS 9 × 5,000 0.740 (0.725), ACCESS 1 × 20,000
   0.944 (0.942).
2. This record, the plan-first result and the plan pages (PR from `chore/plan-first`).
3. `release/6.7.0` from develop: the version, the CHANGELOG section, the gate; a PR to main. The
   operator merges; the tag publishes. Then main into develop, develop to 6.8.0-dev, and #196 and
   the milestone close.

**Plan-first:** the trial did not meet its bar (`docs/harness/plan-first.md`, "Result"). Operator
2026-11-06: try it again on the 6.8.0 work, with a new bar fixed before the first count-changing
build.

**Open:** the guard's masks depend on `--min-baseq` (#374); the SNV tie rule with a confident
third base (#367); Nextflow resources wait for the sweep (branch `feature/nextflow-resources`, local).

## Next — 6.8.0 (tracker #420)

- The judge becomes the default (#331), with #405, #406, #410, #414, #417 and #418; the judge
  meets the speed floor first (#335).
- Packages 4 to 7: evidence outside the aligned bases (#259 with #177, #179, #212); guard
  statistics and artifacts (#306, #180); retiring the homopolymer twin (#261 with #111, #112,
  #128); complex-event read-on (#274, #361). Also #290, #293, #178, #408 and #409.
- Moved from 6.7.0 on 2026-11-06: speed (#412, #312, #150); clean-up (#313 with #195, #314, #315,
  #316); the truth arms and the release gate (#260 with #266, #267, #268, #275); output and
  hygiene (#273, #132, #127, #138); the pre-release review (#318).
- The order is decided with the operator when 6.8.0 starts: count changes first, then speed and
  structure once the counts settle.
