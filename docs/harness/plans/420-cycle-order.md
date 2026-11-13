# #420: the order of the 6.8.0 work

Status: decided 2026-11-12 (operator: "Ok lets do meticulously"): D1-A, D2-B, D3-B, D4-A, with
D3-B applied to D4's lists (section 2, "Decided"). Stage 1 started the same day. Sources: the tracker (#420)
and the 45 open tickets on its path, read with every comment on 2026-11-12; `CONTINUITY.md`; the
harness notes. Measured on 2026-11-12: the test suite under the judge (section 1.3) and the RJ-34
bar's rows after #433 (section 1.2). Measured on 2026-11-13: the judge's speed at production
sizes (the HPC batch pf73, section 1.4). Patient rows appear as counts and shapes only.

## 1. Where 6.8.0 stands

### 1.1 Done since 6.8.0 started

#430 (RJ-29: a realigner's placeholder tags are no evidence of the aligner's doubt; the first item
of the plan-first trial), #432 (RJ-33 and RJ-34 in the register), #433 (the held-back reads named:
`TIED_WITH`, and the found alleles a listed row explains), #434 (the docs after the merge).

### 1.2 The bar for the default (RJ-34), after #433

RJ-34: every signed-out row that 6.6.0 counts within the sign-out's margin and the judge counts
below it is either named (A) or has a cause on the fix list or a decision behind it (B). 37 rows
were in that group on 2026-11-11 (the 6.6.0 tier, 13,323 signed-out DNA rows).

| Group | Rows | With a name after #433 | Left |
|:--|--:|--:|:--|
| A: the reads carry another allele | 17 | 15 | 2: one row with no allele named (`ZERO_ALT`, `PARTIAL_DOMINANT` only); one unlisted 70>1 delins under the naming floor (#273) |
| B1: ties with a listed neighbour | 3 | 3 (`TIED_WITH`, n 28 to 68) | 0 |
| B2: judge faults (equal REF and ALT scores where the bases differ; 2 places) | 4 | 0 (one row names a found allele; all four carry `UNEXPLAINED_READS`) | 4: a fix |
| B3: 5 reads or fewer | 2 | 2 (`TIED_WITH`, n 22 and 33) | 0 |
| C: not yet read | 11 | the 3 checked carry `UNEXPLAINED_READS` only | 11: read them |

20 of the 37 rows carry a name. The fault and the 11 unread rows are the work left on this bar.

**After Stage 1.2 (2026-11-12; the 11 unread rows and the 2 unnamed A rows read by hand by an
Opus agent, every print through the sanitizer; counts and shapes only):**

| Group | Rows | Change | Left |
|:--|--:|:--|:--|
| A: the reads carry another allele | 24 | +6 SNV rows under a listed delins over the site (#398, decided in 6.7.0; named by `OBSERVED_ALLELE` and `TIED_WITH`); +1 ACCESS INS3 whose moved reads tie a found DEL1 that 54% of the reads write (named in part: `COEXISTING_ALLELE` names the DEL1, not these reads; the flags: [QC Flags](../../reference/qc-flags.md)) | 2 unnamed: a DEL6 with substitutions in cis, scanned at 4.8% against the 5% floor (#306); the 70>1 delins the scan's window cuts short (#273, #306) |
| B1: ties with a listed neighbour | 4 | +1 MNV2 whose first base is under Q20 (decided 2026-11-02) | 0 |
| B2: judge faults | 5 | +1 SNV at the second place (3 of its 7 moved reads are the listed INS2 row's) | 5: the fix (Stage 2.1) |
| B3: 5 reads or fewer | 3 | +1 SNV (4 reads) | 0 |
| B4: an open ticket | 1 | new: an SNV where a third-base state ties ALT for 10 reads from a dense mutation cluster (#367, the tie rule with a confident third base) | 1: #367 joins the path to the default |

The ACCESS INS3 row is the plan's open question: a tie with a found allele is not named by
`TIED_WITH` (partial or depth).

**The judge fault's mechanism (Stage 1.1, 2026-11-12; an Opus agent on made-up sequence, then
a trace build on the real run):** at the first place (two SNVs, a found 1,086-base deletion 10
bases on), `needed_by` (`judge/haplotypes.rs:1565`) asks whether the reads need the found
deletion as context of only the first 400 reads that reach it, in coordinate order
(`ATTEST_MAX_READS`), and keeps it only where 2 reads and 5% of them need it. The reads that
need it are its carriers with a short piece before the junction; they start last, so none is
weighed, the deletion is dropped, and every state ties for them. On the real run, with no cap
the judge gives 6.6.0's counts (201 and 250 ALT, against 148 and 206 with the cap); no other of
the run's 434 rows changes. On made-up sequence it fires for a found deletion of about 30 bp or
more with its junction within about half a read after X. The second place (a listed INS2 and
DEL2; 3 rows) has another trigger: on reads with many differences, REF's path skips a short
piece of the read and ties ALT. It is a separate scorer question.

### 1.3 The gate for the default (#331)

| Item | State | Number |
|:--|:--|:--|
| J1, tests under the judge (#332) | not started | 201 failures on develop 3be1765b (147 in the 6.7.0 round; that list is not kept, so the rise is not attributed). 103 are the judge refusing RNA. 98 are DNA tests in 22 files: the read-judgment spec cases 29, the masked-repeat contract 15, the long-deletion contract 10, the locus-alleles and complex-exact contracts 8 each. Neither binning-invariance failure breaks the invariant: one compares the rules path's routing, which a judge-only suite cannot do; the other is an RNA refusal. |
| J2, RNA under the judge (#333) | not started; the judge refuses RNA | the tier has 33 RNA runs of 411 |
| J3, acceptance (#334) | not started | the batch arms for `--min-baseq` 0, 10, 30 exist (batch680) |
| J4, speed (#335) | measured 2026-10-30 and 2026-11-13 (pf73), not met | met in 2 of 8 cells; table in 1.4 |
| J5, HPC gate (#336) | not started | |
| J6, switch-over (#337) | not started | 14 TEMPORARY judge switches to remove; `--rescue-homopolymer` to remove (#261) |
| Before #331 | #401 (a found allele's share), #405 (DP for off-anchor reads) | not started; batch680 holds their arms |
| Decide with the default | #410, #364, #414, #397 (the MNV and delins family), #418 (a stopgap to remove), #406 (memory) | measured in part |
| Acceptance rows (the judge already counts them right) | #417, #425, #290, #276, #293 | close at the flip |

### 1.3a J1's 98 DNA failures, sorted (Stage 1.3, 2026-11-12)

Sorted by a Sonnet agent from the test texts, the failure reasons and the register; no test was
edited. Checked: the 98 are exactly the 201 failures less the 103 RNA refusals.

| Group | Tests | What it is | What it needs |
|:--|--:|:--|:--|
| R, the judge breaks a decided rule | 3 | One gap: the judge has no candidate for a read's own indel beside the site (RJ-7, RJ-8, RJ-20) | a red test, then the fix |
| P, the expectation encodes the rules path | 19 | A homopolymer fixture where REF ties ALT one unit along (4); RJ-34's naming (2); decided by design on 2026-10-28/29 (3); VNTR and duplication fixtures (8); another allele is partial (2) | 3 operator notices, by group |
| H, harness | 23 | Spliced reads, for J2 (7); rules-path warnings and traces (6); strict xfails the judge now fixes (5); two FASTAs beside a BAM (2); 3 others | no expectation change |
| U, no decided rule settles it | 53 | U-W: a read that does not hold the whole window, counted by one confident telling base (23); U-S: a read that spells another complete allele, partial or neither (12); U-E: one confident mismatch inside the ALT event (7); U-F: another allele on 2 fragments inside the event (4); U-run: a run of changes from the window's last base (5); U-N: a listed neighbour beside a deletion (2, #297) | six decisions |

Check on U-F: RJ-26 says that a substitution inside the event on 2 or more fragments and 5% is
another allele; RJ-25's floor of 3 fragments is for alleles nobody listed. So U-F may be an R
case. It is read again when J1 is worked (Stage 4.1).

### 1.4 Speed (the floor: faster than the C++ GetBaseCountsMultiSample 1.2.5)

Measured on HPC with the pf71 kit's production tasks (whole nodes, untraced, 3 repeats, the
median wall time): on 2026-10-30 (pf71, develop fe34f443) and on 2026-11-13 (pf73, develop
3be1765b, build dv680). The rules path met the floor in all 8 cells both times (0.39 to 0.94 of
the C++ time). The judge met it only at 1,000 rows:

| Data | Rows × BAMs | Judge / C++, pf71 | Judge / C++, pf73 | Rules / C++, pf73 | C++ s | Rules s | Judge s | Judge less rules, s | C++ less rules, s | Cut needed in judge less rules |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| IMPACT | 1,000 × 1 | 0.75 | 0.76 | 0.42 | 14.6 | 6.1 | 11.1 | 5.0 | 8.5 | met |
| IMPACT | 5,000 × 1 | 1.31 | 1.36 | 0.46 | 29.9 | 13.8 | 40.7 | 26.9 | 16.0 | 40% |
| IMPACT | 20,000 × 1 | 2.84 | 2.98 | 0.55 | 46.1 | 25.3 | 137.5 | 112.2 | 20.8 | 81% |
| IMPACT | 5,000 × 10 | 1.63 | 1.64 | 0.58 | 235.7 | 135.8 | 386.3 | 250.4 | 99.8 | 60% |
| ACCESS | 1,000 × 1 | 0.82 | 0.80 | 0.39 | 21.2 | 8.3 | 17.0 | 8.7 | 12.9 | met |
| ACCESS | 5,000 × 1 | 2.06 | 2.33 | 0.62 | 30.3 | 18.9 | 70.6 | 51.6 | 11.4 | 78% |
| ACCESS | 20,000 × 1 | 5.15 | 5.27 | 0.94 | 46.2 | 43.2 | 243.2 | 200.0 | 3.0 | 98% |
| ACCESS | 5,000 × 9 | 2.85 | 2.96 | 0.75 | 117.2 | 88.2 | 346.9 | 258.6 | 29.0 | 89% |

"Judge less rules" is the time the judge takes beyond the rules path on the same task. "C++ less
rules" is the most it may take and still meet the floor, if the work both paths share (the
fetch, the pile-up, the output) stays at the rules path's time. The judge's peak memory is 0.12
to 0.32 of C++'s.

Worked example, ACCESS 1 BAM × 20,000 rows (pf73): C++ 46.2 s, the rules path 43.2 s, the judge
243.2 s. The judge must become 5.3 times faster to meet the floor there: its extra 200.0 s must
fall to 3.0 s, a cut of 98%. The extra grows faster than the rows. From 1,000 to 20,000 rows it
goes from 8.7 s to 200.0 s on ACCESS (about 10 ms for each added row) and from 5.0 s to 112.2 s
on IMPACT (about 5.6 ms), while C++ grows by 25.0 s and 31.5 s.

The levers measured so far give tens of percent, not 5 times: one cache slice per row (#150)
saves about 10% of the judge's CPU; per-variant read slicing (#312) saves 14% on an HG002 arm; the
SNV rows routed to the judge (#412) cost 2.25 and 2.81 times the rules path on two HG002 arms.
The judge changed after 2026-10-30 (the fit rule, the background, `TIED_WITH`); from pf71 to pf73
its ratio moved by -0.02 to +0.27. C++, the same program in both batches, moved by -4% to +18%
between pf72 and pf73 on the same tasks (node to node), so moves of this size are not attributed.

**Speed sizing, first results (Stage 1.4, 2026-11-12):**
- The survey (Fable; `speed_survey.md`, local): the options that cannot change a count by
  construction (one cache slice per row, one backward walk, a SIMD kernel, one read view per bin,
  shared site work, reused bands, build flags) give together about 2 to 2.9 times (an estimate:
  the gains overlap). 5 times needs a fast path for reads that spell a state exactly, whose ties,
  paths and ends must first be proved the same. Its basis: the HPC profile of 2026-10-31 on
  #150 (the read call 71.3% of the judge's CPU, the site build 10.9%, the shared work 14.8%).
- The judge's scoring time by read class, from the traced batch680 judge arm (44b02197; the
  `judge fits` line's `ns=`; traced and 10 jobs at a time, so only the shares are read; the clock
  covers the scoring call, not the whole read call):

  | Set | Reads under 2 µs (cache) | Best score 0 (exact) | Other scored reads |
  |:--|--:|--:|--:|
  | ACCESS (across, 155 runs) | 36.5% of reads, 0.4% of the time | 57.9%, 88.9% | 5.6%, 10.6% |
  | IMPACT (tier, 80 runs) | 42.3%, 0.7% | 49.0%, 73.7% | 8.7%, 25.7% |
  | HG002 (7 sets) | 0.5%, 0.0% | 73.6%, 63.5% | 25.9%, 36.5% |

  An exact read costs about 160 µs on ACCESS (89 µs on IMPACT): the whole programme over every
  state for a read whose best score is 0. On ACCESS the reads that are not exact hold 10.6% of
  the scoring time, under the survey's 15% line, so an exact-read fast path could carry the gap
  there if its results are proved identical. Not yet measured: the scoring call's share of the
  whole judge on the 20,000-row production cell. pf73 ran untraced, so it gives the whole time
  (the table above) and not that share; a profile of that cell gives it.
- The HPC batch pf73 (2026-11-13): the table above. 24 tasks, each with C++, the rules path and
  the judge; every run exited 0.
- The local HG002 timing did not run: the load average stayed between 6 and 16 for 2 h 40 min,
  above the agent's limit of 3.

### 1.5 The milestone

58 open tickets in 6.8.0. #364, #397 and #401 are under #331 but have no milestone. 36 of the 58
are on the path to the default (D4; #429 among them, to close now); 22 are not.

## 2. Decisions for the operator

### D1. What comes first in the count work

- **A (recommended): the judge fault first,** on a new ticket. At two SNVs 9 bases apart, 6.6.0
  counts 201 and 250 ALT reads (the sign-out 202 and 250); the judge counts 148 and 206. A found
  1,086-base deletion 10 bases on (184 reads) is the only extra item: its context ties the reads
  that hold its junction with at most 40 bases before it, with equal scores on every state (-43 to
  -55), so the site is never weighed. With the background switch off, the judge gives 6.6.0's
  counts there. The mechanism inside the scorer is not found. A scorer fault can move reads in any
  class, so class measurements made before its fix may have to be done again.
- **B: the "before #331" decisions first** (#405, #401), the fault after.

### D2. Speed: when to measure and design

- **A: as decided on 2026-10-31,** speed after the count work. Risk: the gap is 5.2 times at the
  worst case and the known levers give tens of percent. A design that needs a change of
  structure, found late, moves the release.
- **B (recommended): measure and design now, build after the counts settle.** One HPC batch (the
  pf71 kit, untraced, whole nodes) on develop 3be1765b; a local profile on HG002 (public data); a
  survey of options by a Fable agent. No speed code until the counts settle, so no count moves.
- **C: decide after B's numbers** whether the flip waits for the floor.

### D3. RNA in 6.8.0

- **A: the judge counts RNA in 6.8.0.** Build #333 (the judge on a spliced reference) with #409,
  #293 and #178. Not started. No measurement of the judge on RNA exists; 103 of J1's 201 failures
  are its refusal.
- **B (recommended): DNA flips in 6.8.0; RNA stays on the rules path** until #333 lands (6.9.0,
  before #328 removes the rules path). RNA users see no change in 6.8.0 (33 of the tier's 411
  runs). This changes #331's "done when" (all sub-issues closed): #333 moves to 6.9.0.
- **C: no flip until RNA is ready.**

### D4. The milestone

- **A (recommended): 6.8.0 keeps the path to the default; the rest moves to 6.9.0.**
  - Keep (36 of the 58, plus #364, #397, #401 and the fault ticket): #420; #429 (to close now);
    #331 with J1 to J6 (#332 to #337); #364, #397, #401,
    #405, #406, #410, #414, #418; #417, #425, #290, #276, #293; #412, #150, #312; #313 with #195;
    #261 with #111, #112, #128, #145, #146; #306, #408, #360; #266 and #275 (truth sources for J3);
    #318.
  - Move (22): #127, #132, #138, #177, #178 and #409 (with #333, under D3-B), #179, #180,
    #212, #258, #259, #260, #262, #263, #267, #268, #273, #274, #314, #315, #316, #361.
- **B: keep all 58.**

### Decided (2026-11-12)

The operator took the four recommended options: D1-A, D2-B, D3-B, D4-A.

D3-B changes D4-A's lists. The draft kept #333 (J2, RNA under the judge) with #331's sub-issues
and #293 (a per-read reading flank in a short exon; RNA, validated on FORTE) with the acceptance
rows. Under D3-B, RNA stays on the rules path in 6.8.0, so both move to 6.9.0 with #409 and #178.
The lists are then:
- Keep in 6.8.0: 34 of the 58, plus #364, #397, #401 and the fault ticket.
- Move to 6.9.0: 24 (the 22 above, #333 and #293).

#331's "done when" changes with D3-B: its item 3 ("RNA mode runs under the judge, validated on
FORTE") moves to 6.9.0 with #333. The change to #331's text was posted on 2026-11-13, after the
operator read its draft ("send it").

## 3. The plan in order (with the recommended options)

Each count-changing step has a plan page with its bar fixed first (`docs/harness/plan-first.md`,
the 6.8.0 trial). Its options are emulated from the batch680 arms where they can be; there is one
final build per step; every step tests every data type. The models are chosen per task: Fable for
research and adversarial reviews, Opus for deep reasoning, Sonnet for routine work, Haiku for
mechanical scans.

### Stage 0: housekeeping (no count; small)

| Step | What | Who |
|:--|:--|:--|
| 0.1 | Close #429, fixed by #430 (the closing comment shown first) | me, on "send it" |
| 0.2 | Add #425 to #420; give #364, #397 and #401 the milestone; the moves of D4 | me, on "send it" |
| 0.3 | Retire #430's builds and worktree; the housekeeping script for the old worktrees | the operator approves the script plainly |
| 0.4 | Greptile's review credits (it reviewed none of #432 to #434) | the operator |
| 0.5 | Record the trial's item 1: #430's page met every bar (1 round, no condition added, 1 build, no HPC, 0.39 M agent tokens) | me |

**Posted 2026-11-13 (operator: "send it"):** 0.1 and 0.2; #367 (the SNV tie rule with a confident
third base) got the 6.8.0 milestone and stays under #332, which is under #331; #333 moved under
#328; #331's item 3 and #420's body changed, each with a comment; the fault ticket is #435, a
sub-issue of #331. 0.3 and 0.4 are open. #259 (package B) is in 6.9.0 while #276, #290, #306 and
#360 under it stay in 6.8.0, as D4-A's lists give.

### Stage 1: the facts that set the rest (measurement only; no count code)

| Step | What | Done when | Model |
|:--|:--|:--|:--|
| 1.1 | The judge fault: a ticket, a plan page (the trial's item 2), a trace-only build at the two places | the mechanism is named, with a worked example | Opus; a Fable review |
| 1.2 | Read the 11 unread rows and the 2 unnamed A rows; then the 152 rows above the sign-out that are not SNV or DNP; then the 48 rows the judge brings into the margin | each row is in a class that is decided, on a ticket, or new | Opus for the classes; Haiku for the scans |
| 1.3 | Sort J1's 98 DNA failures: a decided rule the judge breaks, or the rules path's behaviour | a table for the operator; no test edited | Sonnet |
| 1.4 | Speed sizing (D2-B): one HPC batch, a local profile, a survey of options | each option has a measured or bounded gain | Sonnet for the runs; Fable for the survey |

### Stage 2: the count decisions before the default

| Step | What | Note |
|:--|:--|:--|
| 2.1 | Fix the judge fault (from 1.1) | before any class below (D1-A) |
| 2.2 | #405: which off-anchor reads count in DP (the trial's item 3) | two options, both emulated from the rules and judge arms |
| 2.3 | #401 with #306: a found allele's 5% share (with the tied plan's part 3) | one class: the share |
| 2.4 | The MNV and delins family: #410, #364, #414, #397 | the operator: "Decide it with this ticket's MNV family" |
| 2.5 | #418: the judge weighs every allele a REF read reaches | the stopgap goes at the flip |
| 2.6 | The new classes from 1.2 | each on its own page |
| 2.7 | A new baseline batch on the settled code (all four sets, local) | the arms for stages 3 and 4 |

The trial's verdict comes after 2.2: every bar on two of the three items, and bar 2 on all three.

### Stage 3: speed (count-neutral by construction)

#335 with #150, #312, #412 and #406, built in the order of the gains measured in 1.4. Each change
proves that it moves no count (binning invariance; every arm identical). Done when the floor holds
on the pf71 kit (HPC) and the judge stays within 2 times the rules path's wall time per arm.

### Stage 4: the default

| Step | What |
|:--|:--|
| 4.1 | J1 (#332): a decided rule that the judge breaks gets a red test, then the fix; an expectation that encodes the rules path changes only after an operator notice (one notice per group from 1.3) |
| 4.2 | J2 (#333), by D3 |
| 4.3 | J3 (#334): every data type and input; two truth sources (HG002, #266; the spike-in arm, #275, or the oracle: the operator chooses on the page); `--min-baseq` 0, 10, 20, 30 with #408 and #360 |
| 4.4 | J5 (#336): the HPC gate (each batch asked for first) |
| 4.5 | J6 (#337): the defaults flip (CLI, Nextflow); the rules path is deprecated; the 14 TEMPORARY switches go, each keeping its decided default (#313); `--rescue-homopolymer` goes (#261, which closes #111, #112, #128's half, #145, #146); the migration notes name `TIED_WITH`; the acceptance rows close |

### Stage 5: the release

#318 (the pre-release review); the speed floor on the pf71 kit for the new default; the release
gate; the tag (the operator pushes it).

## 4. Risks

1. **Speed.** 5.3 times at the worst case (pf73: ACCESS, 20,000 rows × 1 BAM), against levers of
   tens of percent. D2-B finds out early.
2. **The judge fault.** Its mechanism is not known; it may be a class, not two places.
3. **J1's size.** 98 DNA failures; each expectation change needs an operator notice. One notice
   per group keeps the operator's time bounded.
4. **RNA.** Without D3-B, #333 is the largest single build in the cycle, and it is not started.
5. **No bot review.** Greptile reviewed none of #432 to #434; each item's Fable review carries the
   load until the credits are back.
