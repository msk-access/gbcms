# #331: the fits-nothing rule's "aligner doubts it" test on realigned reads

Status: option decided 2026-11-10 (C with the placeholder pair, shared by both paths; operator: "I agree with this now"); approved 2026-11-10 (operator: "Approve the page"; the trial bar with scope A: "I am ok with A") · the first 6.8.0 count-changing item under the plan-first trial
Rounds: 0 of 2 · New conditions: 0 of 2 · Builds: 0 of 1 (the options were emulated; one final build) · HPC batches: 0 of 1

## 1. The class

- **Reads:** reads the sequence judge finds to fit no candidate sequence (best score below ten
  confident mismatches) whose AS and XS tags do not describe their current alignment, because a
  realigner moved or placed them and kept the aligner's tags. On MSK data: reads BWA left unmapped
  that ABRA2 placed on its contig (`YO:Z:N/A...`, BWA's unmapped default `AS:i:0 XS:i:0`, and a
  MAPQ that ABRA2 sets to the highest it has seen).
- **Evidence:** the aligner's doubt about the read's present placement. RJ-29: "only the aligner's
  own doubt about its placement marks a read that may belong to another locus".
- **Invariants:** 7 (a read is credited only where its own bases carry the ALT: unchanged, the rule
  only stops a read being made depth only); 1 (per-read tags: no bin dependency); 2 (the quality
  gate: unchanged); both paths: RJ-29's placement test has two copies, the judge's fit rule
  (`fits_elsewhere`) and the rules path's junction judge for deletions of 50 bp or more
  (`placed_here`); one shared test serves both.
- **Principle:** a tag that does not describe the read's present alignment is no evidence of doubt.
- **DNA and RNA:** DNA realigned by ABRA2 (IMPACT, ACCESS) holds the class. DNA aligned by BWA only
  (WES, HG002) holds no such read, so it is unchanged by construction. RNA (STAR, FORTE) is not
  realigned and STAR writes no XS score; its AS-0 reads carry no XS, so the placeholder pair (AS 0
  *and* XS 0) never matches them. The judge refuses RNA today. RNA's own doubt signal is `NH > 1`
  (MAPQ below 255), already read for admission in `rna.rs:62`; whether the fit rule weighs it under
  the judge is decided with #333 (J2), where the shared test takes it as one more input.
- **Instance:** 26 signed-out SNV rows below the sign-out on the 6.7.0 release panel (judge680/rbr).

## 2. What already touches it

- **Rule:** RJ-29 (read to its end; decided 2026-10-25, #304 phase B, option `xs`).
- **Code:** `rust/src/counting/judge/haplotypes.rs:108-113` (`fits_elsewhere`: XS >= AS, else
  MAPQ 0 where a tag is missing); `judge/decide.rs:1129` (the tags read: AS, XS, NM, MAPQ);
  `judge/decide.rs:265-292` (`fits_nothing`, `FITS_NOTHING` = ten mismatches); the measurement
  switch `GBCMS_JUDGE_FITS` (`decide.rs:238-256`, TEMPORARY). The rules path:
  `rust/src/counting/junction.rs:516-533` (`placed_here`, the same test), used at `:774` (another
  unit count's reader: partial where placed here, else depth only) and `:1204` (only reads placed
  here count toward the found-allele presence floor, RJ-25). No other code reads AS or XS; `rna.rs`
  reads `NH`; MAPQ is read by the admission filters.
- **Size on the tier (rules path):** at the 68 deletion rows of 50 bp or more, 247 reads placed from
  unmapped sit in the rows' windows (33 rows); 22 of the 13,287 reads with a 50+ bp deletion in
  their CIGAR were placed from unmapped. Its effect (partial against depth, the presence floor)
  cannot be emulated with a switch: the final build measures it against the 6.7.0 default.
- **Tickets:** #429 (this item, a sub-issue of #331); #304 (closed; the fit-rule options measured there), #331, #332.
  **REJECTED:** none.
- **Survey (2026-11-10, Fable):** no surveyed genotyper or caller trusts AS or XS after a
  realigner; fgbio and GATK3 IndelRealigner repair or delete the tags they invalidate; ABRA2 never
  writes AS or XS (mozack/abra2 `ReAligner.java:617-661`); `AS:i:0 XS:i:0` is BWA's unmapped default
  (lh3/bwa `bwamem.c`); ABRA2 takes an unmapped read only through a mapped mate within 1 kb and does
  not test the mate's MAPQ.

## 3. Options

| | Rule in one sentence | Invariant 7 | Conditions | Will not handle |
|:--|:--|:--|--:|:--|
| A | Keep: XS >= AS by the tags | pass | 0 | drops real carriers that ABRA2 placed |
| B0 | An AS of 0 on a mapped read is no score: the tags count as absent | pass | 1 | MAPQ is ABRA2's, so RJ-29's MAPQ-0 fallback never fires: such reads always count |
| C | ABRA2's `YO` beginning `N/A` (placed from unmapped) means the tags are no evidence | pass | 1 | tool-specific; same MAPQ issue |
| D | As C, and the read's doubt is its mate's: depth only if the mate's MAPQ is 0 | pass | 2 | needs the mate's MAPQ (an `MQ` tag or the mate's record) |
| C+P | The tags are no evidence where a realigner marks the read placed from unmapped (ABRA2 `YO` beginning `N/A`) or where they are BWA's unmapped placeholder pair (`AS:i:0` with `XS:i:0`) | pass | 2 | as C; the pair also covers a stripped `YO` and other realigners that place unmapped reads |

Where C and the placeholder pair differ (pg_tags.py, 2026-11-10, up to 25 BAMs and 20,000 reads per data
type): on every DNA arm `YO` beginning `N/A` and `AS:i:0 XS:i:0` select the same reads, except one
ACCESS-standard read (`YO:N/A`, AS rewritten later: C only). No read has the pair without `YO:N/A`.
WES and HG002 have neither. The header is not a reliable aligner record: IMPACT headers name only
MarkDuplicates, GATK PrintReads and samtools (no BWA, no ABRA2), although the reads carry ABRA2's
`YO`; WES, HG002, RNA (STAR) and most ACCESS headers name their aligner.
RNA: STAR writes AS (the pair's score; 0 at about 0.14% of mapped reads), NH and HI, and no XS here
(STAR writes XS only as a strand character, which the tag reader ignores). The judge refuses RNA
mode today; for #333 the aligner's doubt in RNA is NH > 1 (MAPQ below 255), a separate decision.

## 4. Measurement (fixed before the runs)

- **Class metric:** per option, the reads brought back from depth to ALT or REF, and the share of
  those that carry every covered marker of their row's ALT haplotype (differences shared by the
  reads the shipped judge calls ALT). Bar: 95% or more carry every covered marker.
- **Sign-out:** rows within the margin, per option and per stratum, beside 6.6.0; no stratum worse
  than A by more than 0.5 points and 2 rows.
- **False ALT:** the matched normals: ALT reads at rows whose sign-out normal ALT is 0; bar: no
  increase over A.
- **Data types:** IMPACT tumour and matched normals, ACCESS duplex and simplex (ABRA2 twice), WES
  and HG002 (BWA only: no tag of this kind, so unchanged by construction; checked by a scan), RNA
  (the judge refuses RNA; the rules path's RNA runs must stay byte-identical).
- **Both paths, on the final build:** the judge path must equal the emulation on the tier; the rules
  path is compared with the 6.7.0 default on the tier (changes only at deletions of 50 bp or more,
  each attributed read by read), on WES, HG002 and the RNA runs (identical).
- **Effects map (protocol section 2):** changes: the judge's calls of reads placed from unmapped
  that fit nothing (AD, RD, depth-only counts), the rules path's partial-against-depth calls and
  presence floor at long deletions, and what follows the calls (`vaf`, the fragment counts, mFSD's
  read sets at those rows, the observations export, `gbcms merge` sums). Unchanged, because the
  test is unchanged there: `UNEXPLAINED_READS` (it counts every read that fits nothing, called or
  not), per-transcript counts and ASJD (RNA), the writers' columns, Nextflow (no parameter). To
  check at the build: MNP rescue (opt-in) and tract-cluster claiming.
- **The default path, two routes (review, 2026-11-10):** the rules path's junction judge for
  deletions of 50 bp or more, and the SNV rows it routes to the judge (an SNV beside a listed indel
  or delins, `engine.rs:473-479`), where the fit rule and this fix apply under the default too. On
  the tier both routes are unchanged (13,711 of 13,711 rows equal to the 6.7.0 default).
- **What the fix equals in practice:** in DNA a read at MAPQ 0 is never admitted at the default
  `--min-mapq` (20), so RJ-29's MAPQ-0 fallback never fires on admitted reads; for this class the fix
  equals the fit rule off, and `UNEXPLAINED_READS` stays the visible signal.
- **Cannot see:** a foreign read (a paralog's) that ABRA2 anchored by its mate. D's mate test is the
  only guard; measure the mate MAPQ of every brought-back read; more matched normals (`across/`)
  before the final build. The rules-path half has no real-data firing on the tier (only the
  contract test pins it): a spike-in at a repeat deletion (#275) would show it. A read ABRA2 moved
  (a `YO` with a position) keeps its old AS and XS, which still describe a placement within 1 kb:
  a named gap, with no firing on the tier.
- **Contract tests to keep:** the RJ-29 tests in `tests/test_haplotype_judge_contract.py`; add one
  for the junction judge's placement test (a placed-from-unmapped read at a long deletion that
  spells another unit count is partial, not depth only).
- **HPC:** none needed for the decision (the tier is local); one gate batch for the final build.

## 4b. Related 6.8.0 tickets (what this class touches)
- **#331 J and its parts:** J1 #332 (the new contract tests), J2 #333 (RNA: STAR's `NH` joins the
  shared test there), J3 #334 (the fit rule depends on `--min-baseq` through its confident
  mismatches: in the sweep), J5 #336 (the HPC gate's baseline is the fixed build), J6 #337 (the
  migration notes name ABRA2's tags); #313 and #337 retire `GBCMS_JUDGE_FITS` with the other
  measurement switches (this fix adds none).
- **#259 B and #179 (mapping-bias diagnostic):** #179 compares how ALT and REF reads are mapped; a
  read placed from unmapped would look mapped worse (AS 0). It must use the same shared test.
- **#212 (the fragment end from the MC tag):** ABRA2 rewrites the mate fields of pairs it touches,
  and its README warns they may be inaccurate: measure MC against the mate's own CIGAR on reads
  with `YO` before #212's rule. **#177 (split reads, SA):** ABRA2 does not write SA; a moved read's
  SA describes its old alignment.
- **#180 (orientation evidence):** ABRA2 sets a placed-from-unmapped read's strand opposite its
  mate's, so its orientation is assigned, not observed.
- **#275 (spike-in truth arm):** the residual (a foreign read anchored by its mate) can be tested by
  spiking a paralog's reads into a realigned BAM.
- **#361 (unprepared 50+ bp deletions):** if routed to the junction judge, the shared test applies.
- Not touched: #405, #406, #410, #414, #417, #418, #425, #276, #290, #293, #306, #360, #274, #408,
  #409, #178, speed, structure, truth tooling, output.

## 5. Budget and stop rule
At most two review rounds and the chosen option's two conditions (a realigner's "placed from
unmapped" mark, BWA's placeholder pair), shared by both paths. Past either, stop and take a class
question back to the operator.

## 6. Work and models
Emulation from two traced arms of the shipped build (scripts, no model tokens); Sonnet for the
build and red-first test; Fable for the adversarial review of the option and the diff.

## 7. State
- 2026-11-10, measured by exact emulation (rerun/ ctl and fitsoff arms of src_rc670d on all 308
  DNA and normal runs of the 6.6.0 tier; the control equals the gate on 13,711 of 13,711 rows; the
  per-row identity holds at 50 of 50 changed rows). The rule fires at 50 rows in 3 runs; all 433
  reads it fires on were placed by ABRA2 from unmapped (AS 0).
- Within the margin (13,005 IMPACT-type rows): A 12,410; OFF, B0, BR, C 12,432; D 12,425; 6.6.0
  12,764. No stratum worse than A. ACCESS (636 rows) and normals (70 rows): no change.
- Markers: B0 and C bring 395 reads back to ALT; 387 (98%) carry every covered marker. D keeps 111
  more at depth; 109 of them carry every covered marker. Normals: clip admission adds 0 ALT reads.
- Decided C with the placeholder pair, one shared test for both paths (operator, 2026-11-10).
- Built (one final build, `src_fit680a`: 5e8a0527 plus the worktree's patch). Red-first: three
  contract tests failed as predicted, then passed; seven unit tests; the whole qa-check gate 0 on
  every tool (pytest 2,742 passed, cargo 520). Judge path: equal to the emulation on 13,711 of
  13,711 tier rows. Rules path: equal to the 6.7.0 default on 13,711 of 13,711 rows. RNA: 33 runs,
  94 of 94 rows equal.
- The adversarial review (Fable): no defect blocks the merge; the mechanism is sound. Fixed from it,
  count-neutral: the effects map (the routed SNV rows), the CHANGELOG entry, the RJ-24 row's
  summary of the test, three unit tests, and the fit rule's trace (the tags as written and a
  `placeholder` field).
- The final code equals the measured build on the 13 runs holding the changed rows (1,367 of
  1,367 rows, both paths, diagnostics included). The `across/` arm (155 runs: 132 IMPACT normal,
  10 ACCESS normal, ACCESS tumour and cohort; untraced, `run_across_nt.py`): no row changes in any
  kind, so the fix adds no ALT in a normal there, but the rule never fires on that arm either: the
  named residual stays unmeasured by real normals (a spike-in, #275, would measure it).
- Open: the PR.

## 8. Cost log
- Review rounds 1 of 2; conditions 2 declared, 0 added, 0 removed; builds 1 (the final); HPC
  batches 0; no decision reversed.
- Agents: the survey (Fable) 0.13 M tokens; the build (Sonnet) 0.13 M; the review (Fable) 0.13 M.
  The measurements and emulation: scripts.
