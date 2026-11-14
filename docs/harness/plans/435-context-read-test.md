# #435: which reads the judge's context read-need test weighs

Status: approved 2026-11-13 (operator: "Approve the page"); measured 2026-11-13 and 2026-11-14; decided 2026-11-14: B (the B-even fallback withdrawn the same day: speed is never bought with accuracy) · the plan-first trial's item 2 (`plan-first.md`)
Rounds: 2 of 2 (the option review; the diff) · New conditions: 0 of the chosen option's · Builds: 3 of 3
(two measurement builds, `src_att680` and `src_att681`; the final, `src_fin682`) · HPC batches: 0 of 1

## 1. The class

- **Reads:** at a row X, reads whose call depends on a context item beside X: a found indel
  (no input lists it; 2 reads and 5% of the site's reads write it) or a guard background item
  that respells a candidate sequence over a read. Without the item, such a read's cheapest path
  skips the bases the item explains, every state ties, and the read is depth only.
- **Evidence:** the read-need test (`needed_by`). An item stays only where 2 reads and 5% of the
  reads weighed are REF or ALT with it and are not so without it. The test guards against an
  item that X's own reads write (a long insertion's copy of its flank): as a context it would
  tie them.
- **The defect:** the test weighs the first 400 reaching reads in coordinate order
  (`ATTEST_MAX_READS`). The decided text says "of those reaching it", with no cap. The reads that
  need an item after X start last, so at a deep site none of them is weighed (order bias). With
  no cap, a long deletion's 5% is taken over nearly every read over X, while only its carriers
  with a short piece before the junction need it (dilution).
- **Invariants:** 7, unchanged: a read is credited only by its own bases at X against every
  state; the test only decides which other differences a context explains. 1: every option reads
  the row's own reads; B, C and D also remove the dependence on read order. 2: the scoring and
  its quality gate are unchanged. Both paths: the rules path has no such test, but the SNV rows
  it routes to the judge (beside a listed indel or delins, `engine.rs:473-487`) use it under the
  default.
- **Principle:** a context stays where the reads that could need it do; which reads are weighed
  must not depend on where they start.
- **DNA and RNA:** DNA holds the class wherever a long found deletion lies near a deep row
  (IMPACT, ACCESS, WES; HG002 at long insertions). RNA is unchanged by construction: the judge
  refuses RNA (J2, #333), and routing is off in RNA mode (`engine.rs:478`).
- **Instance:** two signed-out SNVs 9 bases apart in one IMPACT tumour run, with a found
  1,086-base deletion 10 bases on (184 reads). ALT 148 and 206 under the cap; 201 and 250 with no
  cap; 6.6.0 201 and 250 (RJ-34's judge-fault group, place 1).

## 2. What already touches it

- **Rules:** no register row. The decided text is in `docs/proposals/sequence-judge.md:191-198`
  (2026-10-27: a respelling over a read; 2026-10-28: a found indel beside the site). RJ-8, read
  to its end: another indel outside the window is a separate event. RJ-26: 2 fragments and 5%
  "of those fitting the allele", a share over the allele's own reads. RJ-33: the CIGAR never
  decides a call. RJ-34: these reads are held back (B2, judge faults).
- **Code:** `judge/haplotypes.rs:1558` (the cap), `:1565-1586` (`needed_by`), `:2700-2724`
  (`settle`, one item at a time), `:2750-2753` (the items weighed: respellings over a read, and
  found indels beside the site), `:2455-2504` (`beside_found`, from `Candidates::found` with the
  reads that write each indel, `:940-943`); `judge/mod.rs:487-497` (the found scan: 2 reads and
  5% of the site's reads); `judge/decide.rs:368` (`decide`).
- **History:** 16a46a70 (2026-10-26, #304) made the test for respellings; its message names no
  measurement for the cap. d55c5041 (2026-10-28, #304) applied it to found indels beside the
  site, because taking every such indel as a context cost HG002 14 true ALT reads at insertions
  of 1 kb or more, 6 at 50-299 bp, and one wrong call (truth v2dbq20). Its case: chr16:34590143,
  a 1,795 bp insertion whose reads write a 1 bp insertion four bases before it on 639 of 2,393
  reads; as a context it tied 13 carriers.
- **Tests:** `judge/haplotypes.rs:4285` (a found indel beside the site that no read needs is no
  context); `judge/mod.rs:1142`, `:1200`, `:1288` (respellings); `tests/test_beside_site_indel_contract.py`.
- **Tickets:** #304 (closed; the test's origin); #401 (which reads a found allele's 5% share is
  taken over: the same family); #306 (a 5% cliff: one item's fate flips all its carriers); #335
  (speed: the reads the test scores). **REJECTED:** none.
- **Survey (2026-11-13, Fable; code and docs, sources in the scratchpad):** no haplotype tool
  weighs the first N reads by position. Caps are random or per alignment start (GATK
  HaplotypeCaller and Mutect2, DeepVariant, Octopus, ABRA2, freebayes); only htslib's pileup cap
  (bcftools, samtools) drops late-starting reads. A candidate's share is taken over the depth at
  its own position (freebayes 0.05, DeepVariant 0.06, VarDict 0.01, bcftools, Strelka2's
  binomial); none takes it over the candidate's own carriers. Every tool keeps a floor of 2
  reads. Only input alleles are kept without a test (GATK `--alleles`, Strelka2's candidate VCF).
  None has a "need" test or documents a rescue of reads after a candidate is pruned.

## 3. Options

| | Rule in one sentence | Invariant 7 | Conditions | Will not handle |
|:--|:--|:--|--:|:--|
| A | Keep: the first 400 reaching reads by position; 2 reads and 5% of them | pass | 1 (the cap) | late-starting readers at deep sites (the instance) |
| B | Weigh every reaching read (the decided text) | pass | 0 (the cap goes) | dilution: made-up, about 2,170 reads, 57 of 2,051 (2.8%) need it; every reaching read is scored twice per item |
| C | As B, and a found indel beside X is weighed over the reads that write it | pass | 1 (the found indel's read set) | a carrier soft-clipped at the junction is not weighed; respellings keep B's dilution; no surveyed tool takes a share over a candidate's carriers |
| D | A found indel beside X is always a context (d55c5041 undone); respellings as A | pass | 0 (the found-indel test goes) | X's own reads that write an indel: d55c5041's measured cost |

The field's share (the depth at the candidate's own position) is B's here: among the site's
reads, a read that reaches an item beside X covers the item's nearest base. So it is not a
separate option.

**C′, found by the measurement (2026-11-13):** as C, but the reads that need the found indel are
counted over every reaching read, and the 5% is of its writers: kept where 2 reads and 5% of its
writers' number need it. It keeps every item B or C keeps (a writer reaches the item, so 5% of
the writers is at most 5% of the reaching reads), and more: 10 of its 21 keeps on the truth sets
only D keeps, 4 of them with no needing writer (review round 1, section 3c). On every traced
found indel (5,632 unique in the A arm, 1,619 in the `--min-baseq` arms; the first count, 5,702
and 1,703, held duplicate lines) it leaves out none that A keeps. As a rule it can: only where an
indel has more than 400 writers and A keeps it on its first 400 reads; B, the decided text,
leaves such an item out too. 38 traced items have more than 400 writers; A keeps none of them.
Conditions: 3, not 1 (the writers' set; a found indel that also respells takes B's bar, worktree
`haplotypes.rs:2880-2883`; need and share from two read sets). Not on the approved page;
emulated from the arms where its decisions equal one arm's at every item of a row (below).

## 3b. Measured (2026-11-13; against A; build `src_att680`)

| | B | C | C′ (emulated) | D |
|:--|:--|:--|:--|:--|
| The instance's two rows (sign-out margin) | into it | into it | into it | into it |
| Made-up shapes A loses (of 8) | 3 fixed | 8 | 8 (decides as C) | 8 |
| Tier reads to ALT; with X's ALT bases | 194; 100% | 208; 100% | 214; 100% | 226; 99.1% |
| Tier reads lost (off ALT or REF) | 0 | 0 | 0 | 32 |
| Sign-out strata worse than A; normals' new ALT | 0; 0 | 0; 0 | 0; 0 | 0; 0 |
| HG002 truth v2dbq20 (c246, sv, tr, sub61), decidable truth-ALT reads gained / lost | 0 / 0 | 28 / **11** (read by read; first written 17, see 4d) | 0 / 0 on 17 of 21 rows; 4 rows at the final build | class 10-49 bp insertions −1.65 points, −8 reads |
| HG002, batch680's 7 sets, rows changed | 0 | 20 | 0 | 131 |
| Contract tests (Rust; suite under the judge, 201 failures) | pass; 201 | pass; 201 | at the final build | **fails** d55c5041's test; 201 |
| Default path (routed SNV rows); RNA | unchanged | unchanged | unchanged | 3 rows' diagnostics |
| Reads the test scores (against A) | 1.00-2.04; at most 45,040 at one item | 0.25-0.57; at most 8,042 | as B | 0.04-0.26 |
| `--min-baseq` 0, 10, 30 (items kept more, tier) | 4, 4, 0 | 7, 7, 2 | 8, 7, 8 | 636, 631, 390 |

C's 11 lost and 28 gained truth-ALT reads all lie at rows where C leaves out an item A keeps
(read by read); its keeping half moves no truth read. (Corrected in review round 1: C′ is not
"C without that half"; and C′'s "0 / 0" is empty, because none of the rows where C′ differs from
A is a scored truth row. See section 3c.) d55c5041's row
(chr16:34590143) is equal under A, C and D (AD 0, `ZERO_ALT` under the shipped judge too: the
long-insertion work, not this item).

## 3c. Review round 1: the options (2026-11-13; two Fable reviewers, independent)

R1 derived the options from the code and the history before it read sections 3 onward, then
attacked C′. R2 attacked the measurement. Reports and scripts: the session scratchpad
`review435/`. Each finding below was checked against the code or the traces (the last column).

| Finding | From | Severity | Checked |
|:--|:--|:--|:--|
| C′'s HG002 truth result is empty: of its 21 differing truth-set rows (15 distinct), 20 are not truth rows and 1 is outside the benchmark | R2 | blocks the claim | yes: each row against `truth_rows_v2d.tsv` |
| C′ keeps items only D keeps: 10 of 21 truth-set keeps; at 4 no writer needs the item; no measurement shows the needing non-writers are clipped carriers | R2 | changes it | yes: 6 + 4 kept by D only, 11 by C |
| On patient data C′ over B is 13 distinct reads at 2 rows; half of its keeps have 400 or fewer reaching reads (the new share base, not the cap) | R2 | changes it | aggregates only (R2) |
| The 5% of the writers decides 28 tier items and no HG002 item (a 2-read floor at WGS depth); at the tier rows it decides, D's keeping moves 22 REF reads to partial or uninformative | R1, R2 | the 5% works on IMPACT | yes: the A trace and D's read changes |
| C′'s bar rises with the writers and moves with other reads' CIGARs (made-up: 8 needing reads, kept at 129 writers, left out at 312, kept again when the 200 extra reads clip) | R1 | invariant 7 note | R1's runs (made-up) |
| C′ is 3 conditions, not 1 | R1 | changes the text | yes: worktree `haplotypes.rs:2880-2883` |
| A flip counts as need: a read ALT without the item and REF with it (`haplotypes.rs:1579-1582`) | R1 | mechanism note | yes |
| A's kept items cost truth reads too: at the 22 items C leaves out, leaving out gains 28 truth-ALT reads and loses 11 (read by read; first written 17) (the harm d55c5041 meant to stop, still under A, B and C′) | R1 | mechanism | from section 3b |
| Speed: B and C′ score 1.34 times A's reads (tier), 2.05 (truth sets), at worst 22,520 at one item; `sv` 2,412 s traced against 1,347 s under D; stopping early keeps every decision but saves little (most deep items are left out); the HPC trigger in section 4 is met | R2 | changes it | the code: `needed_by` returns a bool only |
| ACCESS and WES hold 5 items each, none differing; no option arm ran at `--min-baseq` 0, 10, 30 | R2 | coverage gap | aggregates only (R2) |
| Holds: `settle` tests each item alone; every changed read lies at a row with a differing decision; C′ never leaves out an A keep on the traces | R2 | passes | yes (property) |

Options the reviews add:
- **B-even:** the test weighs an even, fixed sample of up to 400 reaching reads spread over the
  item's span, not the first 400. It fixes the order bias at A's cost and keeps binning
  invariance; near the 5% the sample can decide otherwise than B. 1 condition (as A). Not measured.
- **E, the net effect (R1's step-1 pick):** keep the item where 2 reads need it and more reads
  need it than it harms (REF or ALT without it, not so with it). One rule for both kinds, no CIGAR,
  the same scoring as B; the only rule that sees the 28 truth reads above. Not measured; not the
  decided text.

Reading: the defect is in the mechanism (the test counts the reads an item helps, never the reads
it harms; every share is a stand-in for the harm). The measurement supports B over A (the order
bias, the decided text, 0 conditions, no read lost). It does not support C′ over B. Dilution
(5 of 8 made-up shapes; 13 patient reads at 2 rows) is the share family's question (#306, #401),
where C′ and E can be measured with a harm count in the trace and a negative shape for
d55c5041's harm.

## 4. Measurement (fixed before the runs)

- **The emulation and its basis:** `settle` weighs each item alone, against the site built with
  only the whole respellings left out. So an option that decides every item as A does gives A's
  site, read for read. One measurement build traces, per weighed item, the decision under A, B,
  C and D (count-neutral), and has a TEMPORARY switch that applies one option. Arms: A traced on
  every DNA run, judge path and default path; B, C and D only on the runs where a decision
  differs from A. The basis is checked first: on 10 runs with no differing decision, B, C and D
  must equal A on every row.
- **Control:** A's counts equal batch680's judge and rules arms on every row (develop 44b02197;
  #433 changed names only).
- **Class metric, per option against A:** the reads whose call changes. Bars: (1) of the reads
  brought to ALT, 95% or more carry X's ALT bases at the run's `--min-baseq` and write the kept
  item or are clipped at its junction; (2) every read moved off ALT is listed and attributed; on
  HG002 none is a truth ALT read.
- **Truth:** HG002 truth v2 (per row): no new wrong call against A; the three long-insertion rows
  of d55c5041 keep A's counts; germline hets move toward 0.5, not away.
- **Sign-out (tier):** rows within the margin per stratum, beside 6.6.0: no stratum worse than A
  by more than 0.5 points and 2 rows; the instance's two rows reach 6.6.0's ALT.
- **False ALT:** the matched normals (tier, across): no more ALT reads than A at rows whose
  sign-out normal ALT is 0.
- **Made-up sequence (#435's class table; scratchpad `fault/`):** the reads lost per shape and
  option, at about 330, 970 and 2,170 reads.
- **Data types:** IMPACT tumour and matched normals (tier, 378 DNA runs); ACCESS duplex, simplex
  and the cohort (across, 155 runs); WES (328, BWA); HG002 (batch680's 7 sets, and the 4 sets
  that hold truth v2d: c246, sv, tr and sub61; sv holds d55c5041's rows; BWA); RNA (the tier's
  33 runs on the default path: identical). Missing: a spike-in of a found deletion after X with
  known carriers (#275). cfDNA through `gbcms merge` follows the counts.
- **Parameters:** the decisions traced at `--min-baseq` 0, 10 and 30 on the tier (count-neutral);
  option arms where a decision differs.
- **Speed:** the reads the test scores, per option and data type, from the trace (a count, not a
  time, so the machine's load does not change it), and the worst site.
- **Contract tests:** the five in section 2; the suite under the judge (201 failures on develop
  3be1765b) gains none.
- **Effects map (protocol section 2):** changes: the judge's calls where an item's fate changes
  (AD, RD, depth only) and what follows them (`vaf`, the fragment counts, mFSD's read sets, the
  observations export, `gbcms merge` sums); the held-back reads RJ-34 names (`TIED_WITH`,
  `OBSERVED_ALLELE`, `UNEXPLAINED_READS`); the default path's routed SNV rows; a register row for
  the test and `sequence-judge.md`. Unchanged, because the test is the judge's alone: the rules
  path's own rows, RNA, per-transcript counts and ASJD, the writers' columns, Nextflow. To check
  at the build: MNP rescue (opt-in).
- **Cannot see:** whether a kept item is right where no truth exists (patient rows); the
  brought-back reads' own bases are the check, and a spike-in (#275) would measure it.
- **HPC:** none for the decision. If the trace shows B or C scoring many more reads than A, one
  speed batch on the production cells at the final build (asked first).

## 4b. Related 6.8.0 tickets

Touches: #401 and #306 (the same share family, Stage 2.3); #335 (the reads the test scores);
#332 (J1: the contract tests); #336 (the HPC gate's baseline is the fixed build); #313 and #337
(the measurement switch never merges). RJ-34's second judge-fault place (shape 2: REF's path
skips a short piece on reads with many differences) is a separate scorer ticket, not filed yet.
Not touched: #405, #410, #414, #417, #418, #425, #364, #397, RNA (#333).

## 4c. Measurement 2: the harm count, C′, E and B-even (planned 2026-11-13)

Operator, 2026-11-13, after review round 1: "Lets do 2 first then we decide but again plan it and
then do it" (2: a second measurement build before the decision). C′, E and B-even cannot be
emulated (E needs a count no trace holds; B-even a new read set), so plan-first bar 3 allows the
build. The final build follows the decision.

- **Build `src_att681`** (Sonnet; a second commit on `measure/435-attest`, never merged; frozen
  with the `c63/freeze.sh` recipe and `src_att680`'s 15 pins). Count-neutral at the default.
  - Trace, per weighed item, beside the present columns: `harm` (reaching reads REF or ALT without
    the item, neither with it), `flip` (REF or ALT both ways, but not the same), `gain` (REF or ALT
    with it, neither without; `need_b` = `gain` + `flip`), `need_e` and `weighed_e` (the even
    sample), and `respells` (a found indel taken as a respelling, worktree `:2880-2883`). With
    `GBCMS_JUDGE_ATTEST_TRACE=2`, one line per read that gains, flips or is harmed: its name, its
    class, and whether it writes the item or is clipped within 10 bases of it.
  - Switch values beside a, b, c, d: `cp` (C′), `e` (E: need over every reaching read, 2 reads and
    more needing than harmed; both item kinds), `be` (B-even: the reaching reads at positions
    floor(i·n/400), i < 400, in start order; 2 reads and 5% of them).
  - Tests first, red before the code: each new value parses; harm, flip and gain on a made-up
    site; the even sample's indices; the default equals the shipped judge.
- **Control:** at the default, `src_att681` equals `src_att680`'s A arm on every column of every
  row of every run; its new trace columns repeat the old ones where they overlap.
- **Arms (local):** A with the new trace on every DNA run (tier, across, WES, HG002's 7 sets, the 4
  truth sets; the judge and default paths), and on the tier at `--min-baseq` 0, 10, 30. Then `cp`,
  `e` and `be` on the runs where their decision differs from A, at each `--min-baseq` as well. The
  basis check again: 10 runs with no differing decision equal A on every row. Estimate: about
  2 hours of local runs.
- **Made-up shapes:** the 8 shapes in `fault/` and R1's four (`g40`, `g40x`, `g40xc`, `u48`), under
  a, b, cp, e and be. Two new ones (a Sonnet agent, made-up sequence only): (1) d55c5041's harm: a
  long insertion at X whose carriers' alignments write a 1 bp insertion 4 bases before it; (2) a
  found deletion beside X whose carriers are soft-clipped at the junction (C′'s stated case).
- **Bars (fixed now), per option against A and against B:**
  1. HG002 truth v2dbq20 on c246, sv, tr and sub61: decidable truth-ALT reads gained and lost, and
     no new wrong call. Report how many scored truth rows each option changes; where it is 0, the
     bar is empty and the page says so. The 22 items where C leaves out what A keeps are the only
     truth rows the class touches so far: report each option there.
  2. Reads brought to ALT carry X's ALT bases: for SNV and MNV rows the bases at `--min-baseq`;
     for indel and delins rows a check by bases (the read's own bases aligned against the REF and
     ALT sequences around X), not by CIGAR. Every read moved off ALT or REF is listed (aggregates
     for patient sets) and attributed to an item.
  3. Sign-out strata and the normals' new ALT, for every option as a whole.
  4. Made-up: reads lost per shape; on the d55c5041 shape no option loses an X carrier that A keeps.
  5. Contract tests: the Rust suite under each switch value; the five tests in section 2 pass and
     the judge suite (201 failures) gains none.
  6. Speed: the reads the test scores per option and data type, the worst item, and wall time on
     `sv` and the 3 deepest tier runs, one run at a time on this machine.
  7. Coverage: items each option decides otherwise than A, per data type; ACCESS and WES are named
     as thin.
- **What follows (fixed now):** if E gains truth-ALT reads at the 22 items (net above 0) with no new
  wrong call, and passes bars 2 to 5, E is the recommendation (it removes the stand-in for harm).
  Else, if C′ or E gains reads over B that carry X's ALT bases, loses none on any bar, and passes
  the d55c5041 shape, it is the recommendation, with its conditions stated. Else B. B-even replaces
  B only if B fails the speed batch at the final build. The operator decides.
- **Cannot see:** whether a kept item is right where no truth row exists (patient rows; the bases
  are the check); ACCESS at production depth (the final build's HPC speed batch can trace it,
  count-neutral, asked first).
- **Data rules:** as section 6. Any file that holds patient run labels, row positions or read
  names is a patient file, whatever its type (`cprime_emulation.json`, `att_reads_*.tsv`, the
  `TRACE=2` lines): every brief names such files, and every print from them goes through `./san`.
- **Budget:** builds 2 of 3 after this one (the measurement builds and the final). Review rounds:
  1 of 2 used; round 2 is the final diff. Agents: Sonnet for the build and the made-up shapes; Opus
  (this session) for the diff check, the scripts and the scoring. HPC: none.

## 4d. Measurement 2: results (2026-11-14; build `src_att681`, against A)

Control: the traced A arms of `src_att681` equal `src_att680`'s on every column of every row, and
the old trace lines are identical. Basis: on 10 runs with no differing decision, C′, E and B-even
equal A on every row and every read call. B's column is from section 3b. B-even keeps the same
4 tier items as B at every `--min-baseq` (trace). On its 2 traced tier runs, B-even's arm equals
B's arm on every row (868) and every read call (455,164).

| | B | B-even | C′ | E |
|:--|:--|:--|:--|:--|
| Items decided otherwise: tier; HG002's 7 sets; truth sets; ACCESS; WES | 4 kept; 0; 0; 0; 0 | 4 kept; 0; 0; 0; 0 | 8 kept; 0; 21 kept; 0; 0 | 40 kept; 10 left out; 33 kept and 8 left out; 0; 0 |
| Bar 1, HG002 truth: scored truth rows changed | 0 (empty) | 0 (empty) | 0 (empty) | 5, all on c246 |
| Bar 1: decidable truth-ALT reads gained / lost; new wrong calls | — | — | — | 31 / 2; truth-REF turned ALT 0, truth-ALT turned REF 2 (and 4 wrong REF and 6 wrong ALT fixed) |
| Bar 1 at the 22 items (C there: 28 / 11) | as A | as A | as A | keeps 18, leaves out 4: 23 / 0, no new wrong call |
| Bar 2, tier: reads to ALT; with X's ALT bases | 194; 194 | 194; 194 | 214; 200 (14 ties: C's delins row) | 223; 207 (16 ties) |
| Bar 2, HG002's 7 sets: ALT carriers moved off ALT; reads to ALT that fit REF | 0; 0 | 0; 0 | 0; 0 | **43** (32 at one 69 bp deletion, 16 in each of 2 sets); 10 |
| Bar 2, truth sets: reads to ALT, with X's ALT bases; ALT carriers moved off ALT | 0; 0 | 0; 0 | 25, 18; **2** | 67, 38; **15** |
| Tier REF reads lost (to partial or uninformative) | 0 | 0 | 0 | **24** (by their bases 22 carry REF, 2 tie) |
| Bar 3: sign-out strata worse than A; normals' new ALT | 0; 0 | 0; 0 | 0; 0 | 0; 0 |
| Bar 4: made-up shapes A loses, fixed (of 8); real negative check (D's 21 lost truth reads) | 3; pass | 3; pass | 8; pass | 8; pass |
| Bar 5: Rust suite; judge suite (201 failures) | pass; 201 | pass; 201 | pass; 201 | pass; 201 |
| Bar 6: reads scored against A (tier, truth sets); the most at one item | 1.34, 2.05; 45,040 | 1.00, 1.00; 800 | as B | as B |
| Bar 6: wall time, untraced, one run at a time: all 4; `sv`; the 3 deepest tier runs (A: 2,002 s; 1,720 s; 13, 188, 46 s) | 2,576 s (+29%); 2,290 s (+33%); 18, 185, 48 s | 2,003 s (+0%); 1,725 s (+0.3%); 13, 185, 46 s | 2,575 s; 2,289 s; 18, 185, 48 s | 2,573 s; 2,288 s; 19, 184, 48 s |
| `--min-baseq` 0, 10, 30: items kept otherwise (tier); ALT reads lost | 4, 4, 0; 0 | 4, 4, 0; 0, 0, — | 8, 7, 8; 0, 0, **11** | 37, 36, 37; 0, 0, 1 |

How to read it:
- **E** does what it was built for at the 22 items: it keeps 18, and leaves out the 4 where the
  item harms 26 reads and 13 need it; there it gains 23 truth-ALT reads and loses none. But at
  rows with no truth it loses ALT carriers. Worked example (HG002, public): the 69 bp deletion
  at 13:73108476 in two HG002 sets. A keeps a respelling, an 11 bp deletion at 73108605. The trace
  counts 2 reads that need it and 9 that it harms, so E leaves it out. In the output, 16 reads
  move off ALT (15 to uninformative, 1 to REF), and by their own bases all 16 carry the deletion.
  The per-item counts do not predict the output here (2 needing reads, 16 reads lost); the cause
  is not traced. E also moves 24 tier REF reads to partial or uninformative. They lie at 6 rows,
  each seen in both traced tier runs with equal counts. At every one of them E keeps an item that A
  leaves out: at 5 rows a found 1 bp deletion 9 to 27 bp from the row, which 4 to 20 reads need
  and 0 to 3 harm; at 1 row (a 3 bp insertion) two respellings 2 and 5 bp away. By their own
  bases, 22 of the 24 reads carry REF and 2 tie. At 2 of the 6 rows the trace counts harm 0, and
  the row still loses a REF read. So the harm count misses effects here too; the cause is not
  traced.
- **C′** keeps more than B, but it brings 14 tier reads to ALT whose bases tie (C's delins row),
  and on the truth sets it moves 2 ALT carriers off ALT; at `--min-baseq` 30 it loses 11 ALT reads
  on the tier (not checked by bases). Its truth result is empty.
- **B-even** equals B on every run measured, row by row and read by read. It scores no more reads
  than A (at most 800 at one item); B scores up to 45,040.
- **Time (bar 6):** B, C′ and E weigh every reaching read, so all three cost the same: +29% wall
  time over the 4 runs, +33% on `sv` (1,720 s to 2,290 s). B-even costs nothing (2,003 s against
  A's 2,002 s). One run per option on this machine (load average 1.7 to 2.9); A ran first and
  B-even last, and they agree to 0.3%, so drift is small. Whether +33% breaks the speed floor
  (faster than the original C++ GetBaseCountsMultiSample) is for the HPC speed batch at the final
  build; this machine has no C++ time for `sv`.
- The two base checks disagree at one row: c246 19:16560880 C>CA (a homopolymer, beside a
  respelling C>CAA). HG002 truth v2d calls 8 of the 10 reads E brings to ALT truth-ALT; this
  page's check fits them best as REF with the C>CAA respelling. At c246 20:18016289 A>C, truth
  calls 5 reads REF that the SNV check (by CIGAR placement) calls ALT; truth is the better judge
  there, because a delins in the TA repeat beside it explains the base.
- Correction to section 3b: C's truth-ALT reads lost are 11, read by read (3 to REF, 8 to
  partial), not 17; C also moves 5 truth-ALT reads off a wrong REF and 1 onto it.

By the rule fixed in section 4c: E gains truth reads at the 22 items with no new wrong call there,
but fails bar 2 (43 ALT carriers moved off ALT on HG002, 15 on the truth sets; 16 tier reads tie;
22 tier REF carriers lost). C′ loses 2 ALT carriers on the truth sets, so it fails "loses none".
So the rule gives **B**. B-even replaces B only if B fails the speed batch at the final build; on
this machine B costs +33% on `sv` and B-even costs nothing, with equal counts on every run
measured, so that batch is the deciding check. The operator decides.

## 4e. Two questions: a 1% share, and a rule per kind of item (2026-11-14)

Operator: "why not 1% then and why not mix solution bepending on the bin that it falls in".
Emulated from the traced A arms (`ATTEST2` lines) and the E and C′ arms; no new build. Each read
is judged as in section 4d (its own bases; HG002 truth where it scores the row).

- **1% instead of 5%, over every reaching read (B at 1%), against B:**
  - HG002's 7 sets: no item changes. At about 35 reads, 5% is about 2 reads, so the 2-read
    floor decides already.
  - Truth sets: 31 items change. E keeps all 31 the same way (C′ keeps 21 of them). None is at
    a scored truth row. By bases, 18 reads are won and 11 lost: 2 ALT carriers move off ALT,
    2 reads move to ALT without ALT bases, and 7 REF reads move to partial or uninformative
    without ALT bases.
  - Tier: 22 items change (16 found 1-9 bp, 2 found 10-49 bp, 1 found 50 bp or more, 3
    respellings). At 19 of the 22 rows, every item decides as in E's arm: 6 reads won, 16 REF
    carriers lost, 14 ties. 3 rows are not measured.
  - At `--min-baseq` 0, 10 and 30, the item counts are about the same (23, 22, 24).
  - So 1% gives C′'s effects plus E's keeping half. It loses reads, so it fails as C′ does.
    The new count also finds that C′ loses the same 11 truth-set reads; section 4d listed
    only its 2 ALT carriers.
- **One rule per kind of item:** E's arm against A, read by read, sorted by the kind of the item
  whose decision differs. Read won / lost:

  | Kind of item; the option | Won | Lost | Note |
  |:--|--:|--:|:--|
  | Found 50 bp or more; kept | 194 | 0 | B keeps these already |
  | Found 10-49 bp; kept | 7 | 0 | 14 ties (C's delins row) |
  | Found 1-9 bp; kept | 24 | 33 | |
  | Respelling, SNV or MNV; kept | 0 | 2 | 6 ties |
  | Found 1-9 bp; left out | 61 | 27 | |
  | Respelling 1-9 bp; left out | 8 | 12 | the bases check and truth disagree (section 4d) |
  | Respelling 10-49 bp; left out | 6 | 40 | 38 at the 69 bp deletion |

  Only found indels of 50 bp or more are all wins, and B keeps them already. Found 10-49 bp
  indels lose none, but bring 14 reads to ALT whose bases tie. Every other kind holds both
  wins and losses, so no rule by kind or size can keep the wins and drop the
  losses. The made-up dilution shape does not occur in the measured data: the one extra
  50 bp item that 1% keeps moves no read to ALT.
- **Why E's counts miss (the 69 bp deletion, HG002 `pri`, public):**
  - The row has 5 tested items: the respellings 73108605 (an 11 bp deletion) and 73108616
    (a 13 bp insertion), and 3 found indels.
  - The per-read trace marks 9 of the 16 lost carriers as harmed by 73108605 and 15 as harmed
    by 73108616. It marks none of the 16 as needing either item. Yet under A, which keeps
    73108605, all 16 are ALT.
  - The code explains this (`haplotypes.rs:2700-2724`). `settle` weighs each item against a
    site that holds every other tested item. The final site keeps only the items that pass. At
    this row, the 3 found indels and 73108616 are in the site during the test, but A's final
    site leaves them out (no read needs them).
  - So need and harm are counted in a site that no longer exists when the reads are judged.
    This holds for every option, A included; it hurts E most because E acts on the harm count.
    It is read from the trace and the code; I have not rerun the site.

## 5. Budget and stop rule

At most two review rounds, and the chosen option's conditions (B 0, C 1, D 0). Builds: one
measurement build (C and D cannot be emulated from existing arms; one build serves all four
options) and one final. HPC: at most one batch, plus one rerun. Past any limit, stop and take a
class question to the operator.

## 6. Work and models

Opus: the class and this page. Fable: the survey (done) and the adversarial review of the option
and the diff. Sonnet: the measurement build (the trace and the switch, count-neutral at A), then
the final build with red-first tests. Scripts: the emulation, the scoring and the read
attribution. Patient-derived prints go through `./san`.

## 7. State

- 2026-11-13: drafted. The mechanism is #435's; the code, the history and the tests are read; the
  survey is done. Approved the same day.
- Measurement build `src_att680` (9bbd2a9e: bd603dfe plus the switch and trace; branch
  `measure/435-attest`, never merged). Control: the traced A arms equal batch680's judge and
  rules arms on every column of every row of 868 and 901 runs, but the diagnostic; every
  diagnostic difference is a name #433 added after batch680 (TIED_WITH, COEXISTING_ALLELE,
  OBSERVED_ALLELE).
- The trace (judge path): the test ran 630 times on the tier (144 runs), 1,182 on HG002's 7 sets,
  5 on across, 6 on WES. Against A, B decides 4 items otherwise (2 runs), C 29 (7 runs; 22 on
  HG002 left out, 7 on the tier kept), D 1,478 (148 runs, all kept). Default path: D only, 16
  items in 6 tier runs. Reads the test scores against A: B 1.00 to 1.37, C 0.25 to 0.54, D 0.04
  to 0.26.
- Made-up sequence (#435's shapes): A loses reads in 8 shapes; B fixes 3 (the others fail on
  dilution, even at about 970 reads for 30 and 60 bp deletions); C and D fix all 8.
- Correction to section 4 (2026-11-13): of batch680's 7 HG002 sets only c246 holds truth v2d; the
  truth bar runs on c246, sv, tr and sub61 (a `truth` set added to the runner).
- The basis checked: on 10 runs with no differing decision, B, C and D equal A on every row and
  all 327,457 read calls. The option arms ran where a decision differs; results in section 3b.
  Not run: D's default path on the truth sets (D already fails the contract bar).
- Open for the decision: the operator's choice; C′'s 4 rows and contract tests at the final build.
- Review round 1 (2026-11-13; operator: "we should do an adversarial review of the recommendations
  and understand the options"): two Fable reviewers; findings in section 3c. The recommendation
  changes from C′ to B; dilution goes to the share family. The 4 open rows are 2 distinct SNV rows,
  not truth rows and not d55c5041's.
- Bar 2's check by bases (`batch680/att2_bases.py`: the read's whole sequence against REF and
  ALT around X, with the row's differing items in place or not, at BWA-MEM's costs: mismatch 4,
  gap open 6, extend 1), first on arm C (2026-11-13): HG002's indel reads all resolve; the tier's
  14 reads that C brings to ALT at its one indel row (a 6-base delins, 7 distinct reads, ending 5
  to 13 bases past X) cost the same under REF and ALT, with the item or without it. By their own
  bases they do not carry X's ALT. So C's tier bar 1 is 194 of 208 (93%), under the 95% bar; C′
  holds the same 14 reads (its own count comes with measurement 2's arms). B: 194 of 194.
- Made-up shapes for measurement 2 (Sonnet, scratchpad `fault2/`, options a to d): (1)
  d55c5041's harm was not reproduced: the 1 bp insertion is weighed as found, A, B and C leave it
  out (need 0), D keeps it and loses no ALT read, because the carriers' bases still fit ALT best;
  a tie needs reads that REF plus the 1 bp insertion explains as well as ALT, which the agent could
  not make without tuning past the real case. So bar 4's negative check uses HG002's real cost
  instead: the items D keeps where it loses truth reads (trace: need_b under 2) under each option.
  Aside: with an exact flank copy, 24 of 165 and 137 of 1,100 carriers are lost under every
  option (the long-insertion work). (2) Carriers clipped at the junction (150 bp deletion 25
  bases after an SNV): at 2,000 reads A loses 129 ALT reads and B, C and D keep them; the needing
  reads are all writers (need over every read = need over writers = 129), so the clipped carriers
  do not need the item, and C′'s stated case is not shown here either.
- Slips in round 1: R2 printed the head of `cprime_emulation.json` (tier run labels and row
  positions) to its own tool output; the brief listed the file without saying it holds patient
  labels. Nothing from it is in the reports. One of my commands moved the session's folder into
  the agent worktree for one step (a `cd`); nothing was written there.
- Measurement 2 (2026-11-13 and 2026-11-14; operator: "Lets do 2 first then we decide but again
  plan it and then do it"): build `src_att681` (93270f78 on `measure/435-attest`, never merged);
  control and basis pass; results in section 4d. The rule fixed in section 4c gives B. In the
  build, one existing test's expectation changed: `parse_attest("e")` was asserted unknown, and
  `e` is now an option, so the unknown case moved to `"f"` (the branch is never merged).
- Slips in measurement 2: I printed the first column of three run lists (the tier's run labels)
  without `./san`; `./san` passes those labels unchanged, so nothing it hides was shown. My
  script printed HG002 rows with alleles of more than 1,000 bases in full (public data; noise).
- Closed 2026-11-14 (section 4d): B-even equals B's arm row by row on its 2 traced tier runs (the
  basis arms share no runs, so that pair is not compared; each equals A); E's 24 lost tier REF
  reads are traced to their items; bar 6's wall times are in.
- Slip: my first B-even-against-B command read its exit code with bash's `PIPESTATUS` in zsh, so
  the code printed empty. I ran it again, writing to a file, and read exit code 0.
- Questions answered 2026-11-14 (section 4e): a 1% share and a rule per kind of item, both
  measured from the existing arms; both lose reads.
- **Decided 2026-11-14** (operator: "ok lets go with recommendation stays the same: B, with
  B-even taking its place if B fails the speed batch."): the test weighs every reaching read
  (register RJ-35). The final build is B, with red-first tests (Sonnet). The HPC speed batch at
  the final build (asked first) decides whether B-even replaces B; a B-even swap would be a
  fourth build, one over this page's budget of 3.
- Final build (2026-11-14; Sonnet): the cap and its const removed from `needed_by`; red-first
  tests: a Rust test (480 reaching reads in coordinate order, the 80 that need the deletion
  last) and `tests/test_read_need_every_reaching_read_contract.py` (the 100-base shape: AD 285
  and 285 under the cap, 369 and 378 after, equal to the rules path and to the reads' own
  bases; the 30-base dilution shape a strict xfail that an invariant failure cannot pass). Gates:
  clippy, cargo test, pytest (2,752 passed, 11 xfailed), the suite under the judge (201 failures,
  the same ids as the baseline).
- Review round 2 (Fable, the diff): nothing blocks. Fixed: the cost stated as averages with no
  bound per item (the worst: 22,520 reaching reads, each scored twice); RJ-35 says it qualifies
  RJ-8, and RJ-8 points to it; ACCESS, WES and `--min-baseq` 30 named as thin; the strict xfail;
  the Rust test in coordinate order. Checked and fine: binning invariance (the reads weighed are
  the row's own window from one uncapped fetch, `engine.rs:1327-1336`; the cap hid no
  dependence), no other user of the cap, every number in RJ-35, the CHANGELOG and the REJECTED
  entry. Open for the PR: speed has no bound per item; `call(kept, r)` is the same for every item
  in `settle` and could be cached (count-neutral; the speed ticket); the routed SNV rows of the
  default path pay the cost too, so the speed batch covers both paths.
- Control (plan-first bar 3): the frozen final build `src_fin682` (a `git stash create` of the
  working tree) equals the measured B arm on both of its tier runs (868 rows, every column) and
  A on the 10 basis runs (454 rows).
- Slips in the final build: the builder wrote two scratch logs to `/tmp`, and in its second pass
  lost `cargo test`'s exit code in a pipe; I ran the whole gate again.
- The speed bar (2026-11-14; operator: "Batch first, 10% bar"). The judge already misses the C++
  floor in 6 of pf73's 8 cells, so the floor cannot choose between B and B-even. One HPC batch
  (pf74) runs A (`builds/dv680`, develop 3be1765b: the same judge code as bd603dfe) and B
  (`builds/fin682`) on pf73's 8 cells, both paths, 3 repeats, in one batch. B-even replaces B if
  B's median wall time is more than 10% above A's in any cell, on either path. Then the PR, with
  the result and RJ-35's final text.
- pf74's result (2026-11-14; HPC jobs 10937759, 10937773, 10937778; 96 of 96 runs exit 0; timing
  lines only copied back; `kit/local/agg_pf74.py`): B/A median wall on the judge path 0.97 to 1.08
  (CPU 0.97 to 1.05), on the rules path 0.88 to 1.07 but for one cell. That cell, ACCESS 1,000 rows
  x 1 BAM on the rules path, is 1.145 (A 8.0 s, B 9.1 s; per repeat 1.13, 1.05, 1.28) with CPU
  0.993, so by the letter B failed the bar there. Peak memory 0.89 to 1.07. The one-machine +33% on
  `sv` came from the measurement build; the final build took 1,952 s there against A's 1,720 s
  (+13%, one run on a loaded machine).
- The bar withdrawn (operator, 2026-11-14: "we dont prioritze speed over accuracy, we need to find
  solutions where speed is gained without loosing any accuracy"). B-even is a sample, so it can
  miss the reads that need an item at a deep site, the defect the cap had; it is no speed option,
  and no time bar switches to it. B stays whatever its cost. The speed work keeps every count:
  exact speed-ups only, each shown count-neutral on every set (decision protocol §6; LRN-20261114-001).

## 8. Cost log

- Review rounds 2; conditions 0 added, 0 removed; HPC batches 1 (pf74); one recommendation reversed
  (C′ to B, by review round 1; the operator decides).
- Builds: 3 (the measurement builds `src_att680` and `src_att681`; the final `src_fin682`).
- Agents: the survey (Fable) 0.15 M tokens; the measurement build (Sonnet) 0.13 M; review round 1
  (Fable) R1 0.19 M, R2 0.21 M; measurement 2's build (Sonnet) 0.20 M and its made-up shapes
  (Sonnet) 0.18 M; the final build (Sonnet, two passes) 0.13 M; review round 2 (Fable) 0.15 M.
  Sum so far 1.34 M of the 2.4 M bar. The fault's
  mechanism (Stage 1.1, an Opus agent): not yet summed from its usage report.
