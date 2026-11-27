# Reads that end inside a run beside the change: both views, DNA and RNA (#334, #333)

> Approved 2026-11-24 (operator: "Lets do both views"; "Approve it").
> One decision for the DNA check (`334-germline-check.md`, mechanisms 1 and 3) and the RNA judge
> (`333-rna-judge.md` §7, cause 1): one class, measured on every data type both plans reach.

## 1. The class

- **Which reads.** A read whose bases beat REF and fit the ALT, but end before the far end of a run (a homopolymer
  or a tandem repeat) inside or beside the change. Its bases also fit the same change with another run length, an
  allele nobody listed. Also the reads that do reach past the run and show another length.
- **Which evidence.** Only the base past the run's far end tells the ALT from its other run lengths.
- **Invariants.** 7 (count the given allele; judge bases, not placement), 2 (one quality contract: masked bases in
  the run), 1 (binning invariance: the windows and fetches grow).
- **Instance** (HG002, public, b37):

  ```
  REF     C TTTTTT TTTTTTTTTTTTT AGATG   C, then 19 T
  truth   C TTTC   TTTTTTTTTTTTT AGATG   8:120818709 CTTTTTT>CTTTC: CTTTC, then 13 T
  del+1   C TTTC   TTTTTTTTTTTT  AGATG   CTTTTTTT>CTTTC: CTTTC, then 12 T
  read 1  C TTTC   TTTTTTTTT|            ends after 9 T: fits the truth and del+1 alike, not REF
  read 2  C TTTC   TTTTTTTTTTTTT AGATG   reads past the run: 13 T, the truth only
  ```

  Each row run alone: the judge credits read 1 to both rows (14 of its 15 ALT reads at the del+1 row end inside the
  run); the rules path calls it depth only. RNA: the judge credits such reads at short pure indels and at SNVs at a
  run's end (`333-rna-judge.md` §7, cause 1: 110 of the judge-only credits).

## 2. What already touches it

- **Rules.** RJ-1 (a read ending inside a pure indel's tract is depth only), RJ-3, RJ-4 (a complex event counts
  exact carriers), RJ-6 (an indel's length is part of its allele), RJ-9 (later bases overturn a fit), RJ-23 (complex
  events read over equal windows, the event grown through repeats; under the judge at a 50 bp+ insertion, a read
  that ties the ALT with another allele must hold both windows, or it is depth only), RJ-27 (a pure indel's windows
  grow through the repeat its change lengthens or shortens: "other lengths of that repeat are other alleles a read
  must span it to tell apart"), RJ-33 (bases, not placement).
- **Code** (develop 6ffde565). Rules: `carrier.rs:1763` (`grown_event`: repeats the change touches, not a run the
  ALT and REF share). Judge: `judge/decide.rs:369` (`decide`, the margin over REF and the listed or found alleles),
  `judge/decide.rs:904` (`conclude_margin`), `judge/decide.rs:988` (the tie-window check, only where
  `site.long_insertion()`).
- **Rejected.** REJ-20261019-001 (keep AD exact: another unit count is partial evidence, never AD);
  REJ-20261016-001 (grow a pure indel's windows on both alleles, or not at all; chosen instead: only on the allele
  that holds its bases); REJ-20261006-001 (one deciding base is enough for ALT where no run lies). Corrected
  2026-11-24: the first draft gave REJ-20261016-001 as the chosen rule.
- **Plans and tickets.** #334 (`334-germline-check.md`), #333 (`333-rna-judge.md`), #331 (the judge as default).
- **Out of this class, measured beside it unchanged.** The judge's margin credit for a read one confident base from
  the ALT and two or more from REF (#334 mechanism 2), and the base partner's false ALT (not yet traced). Each gets
  its own page.

- **The field's practice** (survey 2026-11-24, sources at the default branches of that date):
  - Likelihood tools take View A only against their candidates: a read whose two best alleles score alike is left
    out of AD (GATK `AlleleLikelihoods.java:1276-1281`, 0.2 log10; Octopus `assigned_depth.cpp:71`; Strelka2's
    TOR, "insufficient indel breakpoint overlap"). None adds a run length nobody listed, so with one given ALT
    GATK's rule gives View B. GATK's reference model does test every indel length on the REF side
    (`ReferenceConfidenceModel.java:583-595`).
  - bcftools takes a soft View A: a short tandem repeat that reaches a read's end lowers its indel quality
    (`bam2bcf_indel.c:562-585`).
  - freebayes keeps such reads apart as partial observations (PAO, PRO; `Parameters.cpp:214-216`).
  - CIGAR counters take View B: the original C++ GetBaseCountsMultiSample counts a read whose aligned bases equal
    the ALT over the variant's span and never looks at the run (`GetBaseCountsMultiSample.cpp:1996-2014`); VarDict
    likewise.
  - RNA practice adds nothing: ASEReadCounter is SNP-only, and the RNA path reuses HaplotypeCaller's rule.

## 3. Options

Every option fixes mechanism 3: a read that reaches past the run and shows another length is never ALT (RJ-9).

- **A1, the window (both paths).** Every event's windows (complex, short pure indel, SNV) grow through each run
  inside or touching the change, to one base past its far end, on the allele that holds it, as RJ-27 grows a pure
  indel's own repeat. A read short of that base is depth only; one that reaches it with another length is another
  allele. Invariant 7: pass. Both paths read the same window. Conditions: 2 (what a run is: a unit of 1 to 6 bases,
  2 or more copies; the window's growth). Named limits: imperfect repeats; a run longer than any read leaves every
  read depth only.
  - **Scope chosen 2026-11-25** (operator: "Narrow, SNVs as an arm", after the review in §7). It amends no
    decided rule:
    1. Pure indels: RJ-27 as decided. The windows grow through the repeat the change lengthens or shortens, on the
       allele that holds its bases. A run on the other allele only is beside the change.
    2. Complex events (delins, MNV): through every run on either allele that overlaps the change, starts at its
       first base past it, or ends at its last base before it (RJ-23's "its event, grown through repeats").
    3. SNVs: no window, as decided.
    4. The SNV arm (A1+SNV), measured as its own option: an SNV's windows grow through the run its base lengthens
       or shortens, on the allele that holds it. It amends the contract "SNVs have no window" (C10 #157) only if
       chosen.
    4a. The MNV arm (A1+MNV), measured as its own option (operator, 2026-11-25: "can we test both to see the
       numbers?"): a plain MNP read (no indel at the block) is judged on the block alone, in both paths (RJ-4;
       `check_mnp`). Under the arm, an MNP whose block a run overlaps or touches reads the block and each such run,
       to one base past its far end, in both paths. It amends RJ-4 only if chosen.
    5. Only the windows grow, with the fetch and the rules path's guard reach that must hold them (invariant 1).
       The judge's site, guard and masks stay as on develop. In the rules path the whole-read placement, the RNA
       exon-edge test and the BAQ span follow the grown windows; the guard zone and the long-layout switch stay on
       develop's event (corrected 2026-11-25, after the review and the HG002 trace in §7). A complex event's window
       keeps its two flank bases past a grown run (RJ-4, RJ-23).
    6. **Form chosen 2026-11-26: A1 narrowed** (operator: "Ok lets go with A1 narrowed", after the §7 entries of
       23:41, 11:25, 11:39 and 12:36). The run-next rule (§7, 18:11) is depth only for a read that shows, at the
       first base past the event's run where the guard masks it, a confident base that goes on with the run, unless
       the read's confident base at Q (the first base where "one more unit of the run" and "the guard's change at
       the masked base" differ) shows the change. RJ-23's guard credit stays as decided. Why: A1 as built's extra
       false-ALT cuts (base partners 399 against 671, rules path, HG002 DNA) come from overriding RJ-23 next to runs,
       which is its own class (§7, 18:11); A1 narrowed misses fewer true variants on HG002 RNA (FN +18 against
       +46); on patients the two forms differ by 4 reads. The SNV and MNV arms are built on this form (labels
       snvn, mnvn). A1 as built (o334a1, `.so` 1e30c02b) is not measured further. This settles A1's form only;
       the choice between the faithful options followed on 2026-11-26: A1n, at every site (RJ-36; §7, 16:47 and
       18:42).
- **A2, the run's other lengths as candidates (judge; A1's window in the rules path).** The judge adds the run's
  lengths one unit shorter and longer as candidates; a read that ties the ALT with one is depth only (RJ-23's tie
  rule, at every site). A read past the run with another length falls to that candidate. Above the field's
  practice where our data can test it: masked bases inside the run are judged as a stutter model would. Invariant 7:
  pass. The paths can differ where bases are masked (named). Conditions: 3 (the run, the candidate lengths, the tie
  rule at every site). Cost: more candidates per site (the speed floor applies).
- **B, ALT (both paths).** A read whose bases beat REF and fit the ALT counts ALT, whether or not it reaches past
  the run; run lengths nobody listed or found are not candidates. Invariant 7: pass on "input taken as correct";
  open on "its own bases carry that ALT" (they carry the change, not the run's length). Amends RJ-27's base rule and
  RJ-3 for reads that beat REF inside a pure indel's tract. Conditions: 2 in the rules path (credit those reads; read
  the exact-carrier window on through the run). Named limit: a wrong-length input gets every such read of the true
  allele. The field's CIGAR counters and GATK with one candidate behave this way.
- **Considered, not proposed:** A with these reads kept as partial evidence (`partial_alt`), as freebayes and
  Strelka2 keep them apart. RJ-23 decided such a tied read is "depth only, not partial evidence"; built only if the
  operator asks.

## 4. Measurement, fixed before the runs

- **Per data type and path, by the patient oracle's and the truth sets' per-read labels:**
  1. false ALT at wrong-length rows (the del+1 partners; RNA's REF-by-truth rows at tracts), by label;
  2. ALT at true calls: oracle-ALT reads credited and missed, "either" reads credited;
  3. REF reads that end inside the run (the REF side);
  4. VAF at truth hets beside runs (HG002 DNA and RNA: mean |VAF − 0.5|).
- **Bars** (a build is a faithful option; the operator chooses between faithful options):
  - A1, A2: no "either" read credited ALT or REF in the class; no read past the run with another length ALT.
  - B: every "either" read that beats REF credited ALT; no read past the run with another length ALT.
  - Every option: rows with no run inside or touching the change unchanged in every field; binning invariance
    identical; the decided-rule contract tests kept, except those the option amends, each named to the operator
    with what it asserted before its expectation changes.
- **Data types and truths.** DNA WES: HG002 exome (GIAB; the del+1 slice on this Mac). DNA tissue: IMPACT tumour
  and normal (pat334, patient oracle). cfDNA: ACCESS duplex, simplex and standard, tumour and buffy coat (pat334).
  RNA: HG002 RNA (public; #333's truth runs and internal truth) and FORTE (#333's sets). Missing: WGS, long reads,
  a patient oracle for RNA, a spike-in of complex events beside runs (built if the options come out close).
- **What it cannot see.** Reads the oracle cannot label (unk); mechanism 2 and the base partner, unchanged.
- **HPC.** Per option: the pat334 arms and oracle (471 tasks each) and #333's FORTE runs. HG002 DNA and RNA on
  this Mac.

## 5. Budget and stop rule

Testing is not capped. Stop and go back to the operator on an output difference outside the class, on more
conditions than an option declares, or on a change of design or decision.

## 6. Work and models

Scripts first (the oracle, the reports). Sonnet builds each option on its own branch from develop 6ffde565, with no
measurement switch in shipped code. Opus reviews each build (is the defect in the instance or the mechanism; does
the option break an invariant). Fable for the tool survey.

Build order. Every option grows the rules path's windows through the run (A1 and A2 make a read short of it depth
only; B credits one that beats REF), so that part is built once, first, as A1's rules half. A1's judge half, A2's
judge half and B's rules change are then built from it, in parallel, each in its own worktree.

Read with the options (written down at the approval, 2026-11-24, before any build):
- Mechanism 3 is the rules path's. In the judge, a read past the run with another length falls to mechanism 2
  (unchanged), except under A2, where that length is a candidate. The bar "no read past the run with another length
  ALT" holds for the rules path in every option, and for the judge under A2.
- "Both paths read the same window" (A1): the judge holds the rules path's windows for the event (the event with its
  flank bases, RJ-4, RJ-23, RJ-27), grown through runs. So a read ending at the change where only a single base
  follows (run band 1) is covered by the flank bases. Under A2 a single base counts as a run of one unit: its
  candidates are that base deleted or doubled.

## 7. State

- 2026-11-24: draft; the tool survey is done (§2). Approved the same day.
- 2026-11-24 21:50: the shared part (A1's rules half) is being built (Sonnet, its own worktree, develop 6ffde565;
  every build shares the main checkout's `rust/target`: 29 GB of disk free).
- 2026-11-24 22:41: stopped under §5 and back to the operator. The shared build fixes 6:108243120 (rules AD 42 to
  15) and fails 62 existing tests (57 Python, 5 Rust). The Opus review sorts them: 35 intended by A1 as written,
  5 judge counts moved, 18 guard faults, 4 outside the class. Three faults, all mine or the build's:
  1. The growth sits in the shared event (`carrier.rs` `event()`), so the judge's guard, masks, BAQ span, fetch
     and the RNA exon test move with it. It belongs in a window-only extent.
  2. The brief said "on each allele"; §3 says "on the allele that holds it". The build grows on both.
  3. A1 as written amends decided rules that §3 did not name: RJ-27's growth (2026-10-16: only the repeat the
     change lengthens or shortens; a run on the other allele only is beside the change), RJ-1's window
     (2026-10-15), and the contract "SNVs have no window" (C10 #157). Its run (2 or more copies of a 1 to 6 base
     unit) takes doublets: on random sequence about 97% of SNVs and 69% of 1 bp indels grow.
- 2026-11-25 10:48: the operator chose the narrow scope, SNVs as their own arm (§3, A1). The shared part is rebuilt
  to it in the same worktree.
- 2026-11-25 11:22: the Opus review of the narrow rebuild. 10 existing tests fail (3 Python, 7 Rust):
  1. 4 Rust geometry tests are intended by scope item 2 (no reads move). They go to the operator before any change.
  2. 2 build defects. The whole-read placement (`whole_haps`) still reads `event()`, so test_carrier_window_balance
     [48>46] gives AD 96 against RD 95. The rules-side guard zone grows through runs, which amends RJ-23
     (2026-10-12) and turns AD 997 into 0 in locus_alleles (5 tests).
  3. Both builds started from origin/main (0f5c1f0b, the 6.7.0 release merge), 42 commits behind develop: an
     agent worktree's base (LRN-20261125-001, `git-conventions.md`). The judge's counts match develop's on 19
     cases.
  The builder now moves the change onto develop 6ffde565 and fixes both defects.
- 2026-11-25 11:35: the shared part on develop 6ffde565 (branch feature/334-a1-rules, uncommitted). Both defects
  are fixed. pytest passes (2,776, 11 strict xfails still xfail); `cargo test` fails only the 4 geometry tests
  above. On the HG002 slice (public; 1,402 delins rows and 3,466 wrong-length rows):
  - the sequence judge is identical on every row and field;
  - the rules path changes 113 and 177 rows, all complex: AD falls on 90 and rises on 8, and DP rises on some.
  A stricter window should not raise AD: an Opus check traces those 8 rows before the options are built.
- 2026-11-25 12:16: the trace of the 8 rows (Opus, `--trace` on both builds, each row alone).
  1. Build defect, 2 rows (5:161263177, 3:4872813). The long layout is keyed on the grown window (`long = short >
     LONG_EVENT`, carrier.rs:2060). A 70-base (GAA)n run pushes it past 50, so the hold falls to 16 bases from the
     right, and 11 reads that start inside the run count ALT. Fix: key the layout on develop's event.
  2. One correct rise (8:144652597): a read reads past the GGG run.
  3. 4 false rows: each new ALT read's one telling base past the run is below the quality gate (Q13-16), so it fits
     the false row and the truth alike and counts on length alone. §1 names masked bases in the run as part of the
     class: a decision for the operator.
  4. 5:172625915: the run lies on the longer allele, past the shorter allele's window. RJ-23 reads the longer
     allele on only as far as a read reaches, so a read inside that run counts ALT. A named limit under the chosen
     scope, counted in the measurement.
  5. Of 55 sampled reads that left ALT, 36 end inside the run and 18 reach past it but stop before the two flank
     bases a complex event's window keeps (RJ-4, RJ-23), as on develop for every complex event.
  - The RNA exon-edge test and the BAQ span must follow the grown windows in the rules path. Otherwise a grown
    window can read into an intron, and a read's own gap in the grown part loses its BAQ sparing. Scope item 5 is
    corrected: they follow the windows in the rules path, measured on RNA.
  - The 4 geometry tests: changes found that keep each point (575 lib tests pass on a scratch copy). One changes an
    asserted literal with its fixture. They go to the operator.
- 2026-11-25 12:25, the operator's two answers:
  1. "Approve all 4" test changes: the contexts GGACTAGTACGA and GGACTTTACAGACGTCAG (literal ACGCC to ACGCA), seed
     38, and the flank-choice test's own extent read with `event_for_windows`.
  2. "Depth only": under A1 a read holds a grown window only when its base past the run is confident (at or above
     the minimum quality, not N). Only grown runs change. A1's judge half holds the same rule.
  The builder now fixes the layout key, the exon edge and the BAQ span, adds this rule, and applies the 4 changes.
- 2026-11-25 13:40: the fixes are in (pytest 2,778 pass, `cargo test` 576 pass, clippy clean). HG002 rerun:
  - the judge is identical on all 4,868 rows;
  - the long-layout rows no longer rise;
  - 3 rows still rise (7:100360042, 9:125617391, 19:45452269). 5 reads go from depth only to ALT, 3 of them behind
    junk left clips, starting inside the run. Their only bases that tell REF from ALT are below the quality gate.
  The builder now traces the cause and fixes it. The rise at 5:172625915 is the named limit (§7, 12:16, item 4).
- 2026-11-25 13:48: one cause fixed. The placement from the other flank weighed grown bases that lie in a junk
  clip; grown bases now give no placement evidence (`Windows.grown`). pytest 2,779 pass, `cargo test` 576 pass.
  Two reads still go from neither to ALT (7:100360042, 19:45452269). Each fits REF and ALT alike once its masked
  bases are left out. 8:144652597 falls to AD 111 (develop 114). An Opus audit now sorts every read whose call
  moved on HG002, and traces each unexplained move to its line.
- 2026-11-25 14:20: the audit (Opus). 2,580 reads moved on 482 rows; at least 2,446 (about 96%) as the scope
  intends. The rest come from four mechanisms, each a place where the grown span reaches a decision other than
  "does the read hold the window":
  1. the hold grows past the run plus two flank bases (`shared_past` from the grown end): 73 true calls lost;
  2. the masked-base check reads the base through the CIGAR, not the reading (misses clip bases and own indels);
  3. left growth moves the reading anchor across the read's own indel or clip;
  4. a sibling that is depth only for its grown run no longer blocks REF (RJ-25's tie rule).
  The geometry itself and the masked-base rule are sound. The builder now fixes all four. The audit runs again
  after.
- 2026-11-25 15:00: the four mechanisms are fixed. The shared part (A1's rules half) on HG002, against develop:
  - the judge is identical on all 4,868 rows;
  - in the rules path only complex rows change: the wrong-allele set 228 rows (216 delins, 12 MNV; AD down on 18,
    up on 1), the delins set 134 rows (AD down on 114, up on 1);
  - the two AD rises are 5:172625915 (the named limit) and 19:45452271 (one read; two checks found it correct);
  - 8:144652597 falls to AD 111 and RD 142 (develop 114 and 148): reads that end in the run leave both alleles.
  Not yet traced: 41 reads that move from neither to partial (no AD effect) and the reach of the guard. The review of
  A1's judge half traces them. The shared part is saved as a patch. Four builds start on it, in parallel (Sonnet):
  A1's judge half, A2's judge half, B's rules change, and the SNV arm (rules path, merged with A1's judge half
  later). Each gets an Opus review before any measurement.
  Open conflict for the operator: clarification 2 makes the judge hold the flank bases at every site, also at rows
  with no run (band 1). The bar "rows with no run unchanged" (§4) then fails for A1's judge half. The review counts
  those rows first.
- 2026-11-25 15:01: the operator: "Measure first". A1's judge half keeps the flank bases at every site, as
  clarification 2 says. The build reports the rows with no run that change, apart from the rest, with the reads
  that move. The operator then chooses between every site and only rows with a run.
- 2026-11-25 15:40: two builds are in (not yet reviewed).
  - A1's judge half (Sonnet). The judge now holds the rules windows at every site except long insertions, 50+ bp
    deletions and delins cores, which keep RJ-23 and RJ-24. On HG002 against develop's judge: the wrong-allele
    set 509 rows (AD down on 108), the delins set 397 rows (AD down on 369); AD never rises. Wall time unchanged.
    Rows with no run: 0 of 9 change. Rows where a run lies only inside the change: 19 change, only through the
    two flank bases; almost none of the moved reads reach the base past the change. 4 existing tests fail, all
    from "a read short of the window is depth only" (reads ending inside a delins; reads ending on an MNP
    block's last base). The rules path equals the shared part except on 32 MNV rows that the rules path routes
    to the judge (checked: rules equals judge on those rows in both builds).
  - The SNV arm (Sonnet), rules path only. The run definition (doublets included) gives windows to 371 of 423
    HG002 SNVs (88%; GRCh38 v5.0q truth over a chr1-3 slice). On those rows AD falls 9.5% and RD 9.8% against
    the shared part; the 52 SNVs with no run are unchanged; the judge is unchanged. 7 existing tests fail, all
    from "SNVs have no window", one of them a strict xfail (#409) that now passes.
  - A hazard found: the shared Rust target can hand one worktree another worktree's build. The builds now touch
    their sources before every build, and each `.so` hash is checked against the others.
  An Opus review of A1's judge half runs now; it also traces the 41 partial moves left from the shared part.
- 2026-11-25 15:46: B is in (Sonnet; not yet reviewed). In the rules path a read short of the grown window now
  counts by its bases: the build holds the event's own window plus two flank bases, cut to where the read reaches
  inside the growth, and drops develop's shared-bases extension for these rows. The judge is identical to develop
  on all 4,868 rows. Rules path on HG002 against the shared part: the wrong-allele set 498 rows (AD up on 97, down
  on 0), the delins set 339 rows (AD up on 248, down on 2); REF rises too (1:117122288 RD 92 to 163). Every changed
  row has a run touching the change. 19 Python and 4 Rust existing tests fail, among them RJ-4 and RJ-23 cases
  (reads ending where the alleles still agree: (0,0) to (6,6)). Open, by the builder: a read that starts inside a
  left-hand run follows the aligner's placement, and one with a single mismatch counts REF. An Opus review runs
  now. The SNV arm's review also runs.
- 2026-11-25 15:51: A2's judge half is in (Sonnet; not yet reviewed). It fails the gates: 15 Rust and 22 Python
  existing tests fail. The judge adds one-unit-shorter and one-unit-longer candidates for each run the change
  touches (a single base counts as a run of one unit): 3.6 per site on the wrong-allele set, at most 17. On HG002
  against develop's judge: the wrong-allele set 596 rows (AD down on 103, never up; partial_alt down on 441 rows,
  3,761 reads), the delins set 393 rows (AD down on 354). The truth row 8:120818713 goes AD 40 to 23, equal to the
  rules path; its del+1 row goes AD 15 to 1. The rules path equals the shared part in every field. The judge takes
  about twice the shared part's time (whole-second timing on a loaded machine). Ties the builder could not explain:
  reads reaching 17 to 120 bases past the change at 12:123343571, 7:99507353 and 21:34050940. An Opus review runs
  now, with a clean timing.
- 2026-11-25 16:10: the review of A1's judge half (Opus). Per read on 31 HG002 rows, 0 of 1,306 credited calls miss
  the window or carry a masked telling base. Of 677 moved reads, 473 move as the rule says; 134 end or start
  inside the event (clarification 2: a read short of the window is depth only, at every site); 20 sit behind junk
  clips the judge leaves unplaced; 2 miss by one low-quality end base the alignment leaves out. The builder's split
  of the "run only inside the change" rows was wrong: there, 22 of 25 moved reads end inside the event, and only 3
  fail on the two flank bases alone. Defects, now in a fix round (Sonnet):
  1. padded SNVs (a same-length row with padding, one base changed) got a grown window, in the judge and, through
     the rules path's judge routing (RJ-25), in the rules path (32 rows);
  2. the judge held the block, flank bases and grown run at plain MNPs, where the rules path reads the block only
     (`check_mnp`; RJ-4 counts only MNP reads with an indel at the block as complex): this failed 3 kept-base tests;
  3. the judge's window dropped the RNA exon-edge layout;
  4. minor: a telling column with no placed base was read as confident; "holds the window" differed from the rules
     path by one end base.
  The fix keeps plain MNPs as decided (block only) in both paths, since the scope amends no decided rule. Open for
  the operator: the class reaches plain MNVs beside runs too, and the narrow scope does not cover them. The 41
  partial moves left from the shared part could not be reproduced: the audit's logs give 189 reads on 26 rows, and
  partial_alt rises by 30 on 13 rows, with AD unchanged. They set only the "nearby" flag: event bases in a soft clip
  under a right-grown window (plausibly intended), and, at 1:25780892, a left-grown window crossing each read's own
  2-base insertion (another allele of the ALT's length).
- 2026-11-25 16:20: the review of the SNV arm (Opus). Fetch and guard reach hold (binning invariance identical);
  the 52 rows with no run and the judge are unchanged. Per read on 370 rows, 1,422 calls move; 37 the rule does not
  explain. The AD loss (9.5%) splits: doublet-only rows 11.8%, doublets plus 2-copy runs 50.9%, homopolymers of 3
  or more 22.8%, longer repeats 26.3%; a masked telling base gives 35% of it. Defects, now in a fix round (Sonnet):
  1. a read carrying the SNV and a listed neighbour in the grown run or flank was dropped (chr1:162767066 AD 40
     to 2 with every truth row listed; the judge gives 42);
  2. the sibling contest by bases ran for grown SNVs, so a listed MNV took the SNV's readers (RJ-28);
  3. soft-clip bases counted as ALT for SNVs (9 reads);
  4. judge-routed SNVs kept the arm's flag, so other rows' sibling checks moved;
  5. reads not holding the SNV became partial (10 reads);
  6. the near end of a run ending at the SNV was not a telling point;
  7. the brief's geometry (two flank bases on both sides) lost 79 reads the class does not need: the fix uses the
     plan's "one base past the far end" and no flank on a side with no run.
  The HG002 SNV check now lists every truth row in the slice, SNVs and indels.
- 2026-11-25 16:25: A1's judge half after its fix round (Sonnet). Padded SNVs grow no window (an MNP event grows
  only when more than one of its bases changes; `$SP/a1_shared_v5_delta.patch`, for every option); plain MNPs get
  no judge window, as decided; the judge's window takes the rules path's exon-edge layout; a telling column with
  no placed base reads the next aligned base, as the rules path does. Checked on HG002: the rules path equals the
  shared part in every field; the judge changes only delins rows (the wrong-allele set 453, AD down on 85; the
  delins set 396, AD down on 367); MNV and padded-SNV rows are identical to develop. pytest: 1 existing test fails
  (reads ending inside a delins, clarification 2), unedited. Wall time unchanged. Limits found: the sequence judge
  refuses RNA mode, so the RNA measurement of A1 is its rules half only; a leftover end base the alignment leaves
  out, if confident, still misses the window.
- 2026-11-25 16:30: the review of B (Opus). A credited read still beats the other allele by a confident base over
  the windows it is judged on; no read past the run with another length is ALT; 1:117122288's REF rise is all REF
  reads that pass. B's named limit shows at row level: 7:99507353 truth TT>G AD 17 to 72, wrong TTT>G AD 7 to 62.
  Defects, now in a fix round (Sonnet), which narrows B to what the class needs (only the run part is relaxed):
  1. the cut to a read's reach anchored a read starting inside a left-hand run at its aligned start, so the
     aligner placed it: 10 reads credited REF though their bases favour ALT, 16 though they fit both, 12 ALT lost;
  2. the relaxed hold dropped RJ-23's shared bases past the event that are not the run: 46 truth-row and 20
     wrong-row credits fit another (substitution-only) allele outside B's named limit;
  3. cut reads lost the RJ-23 masks.
  Tests B broke that are not its named amendments (RJ-23's shared bases, the long allele holding the short window
  whole, masked repeats, mFSD classes, reads ending inside the event) should pass again after the fix.
- 2026-11-25 16:31: the operator, on plain MNVs beside a run: "can we test both to see the numbers?". A1 stays
  block-only at plain MNPs; the MNV arm (§3, A1 scope item 4a) is built on A1 (Sonnet) and measured beside it.
- 2026-11-25 16:35: the review of A2's judge half (Opus): it fails, and the stop rule (§5) applies. In the class
  it does what it should: of 1,030 moved reads on 14 sites, 888 fit the ALT and a candidate exactly, 76 go to a
  winning candidate, 0 are unexplained; 8:120818713 AD 40 to 23, its del+1 row 15 to 1. Outside the class:
  1. the tie test (decide.rs:924-926) asks only whether the ALT or REF scores within the margin of a candidate, so
     a read one confident base from the ALT (score -4) ties a candidate one gap away (-7): confident full-span
     reads go to depth only (45 of 1,105 ties), and 758 of 792 partial reads that tell a third allele lose it
     (partial_alt 187 to 1 at 8:145151665);
  2. a candidate takes the ALT's margin before RJ-28 runs: MNP carriers with REF's base at one changed base go
     from partial to depth only;
  3. pure-indel candidates come from runs inside the indel's own bases, not RJ-27's repeat (a 39 bp insertion
     gets 24 candidates, the cap);
  4. pure indels get no single-base candidates;
  5. rules changed without being named: RJ-24 (50+ bp deletions get candidates), RJ-23's long-insertion tie, the
     contig-edge warning.
  Candidates per site: mean 6.5, at most 24 (34 sites at the cap), not 3.6. Speed, 3 runs each on a loaded
  machine: 1.8x (wrong-allele set) and 2.4x (delins set) develop's wall time, and 2.0x its memory, against the
  speed floor. 22 Python and 15 Rust existing tests fail. Back to the operator: fix A2 once and time it again, or
  drop it.
- 2026-11-25 16:56: the SNV arm after its fix round (Sonnet). A read whose own indel is a listed row's, sliding in
  the run, is judged on the SNV base alone; grown SNVs keep develop's sibling handling; the rule applies only when
  the SNV base is aligned, confident and not N; judge-routed SNVs take no window; both ends of a run are telling
  points; one base past each end of the runs' union (a side no run reaches keeps one base, so a longer run there
  is not ALT). HG002 with every truth row in the chr1-3 slice listed (571 SNVs, 136 others): the rules path changes
  459 SNV rows (AD down on 407, up on 0) and no other row; the judge is unchanged. Grown rows: AD -8.5%, RD -9.0%.
  Share of the AD loss: doublet only 14.9%, doublets plus 2-copy runs 42.2%, homopolymers of 3 or more 26.6%,
  longer repeats 16.2%. Per read: 1,645 moved, 0 the rule does not explain; masked telling base 840, short of the
  window 776. The review's examples recover (chr1:162767066 AD 36, against 40 on develop and 2 before the fix).
  2 existing tests fail: "SNVs have no window" (the arm's amendment) and the #409 strict xfail, which passes by
  accident (an SNV in a run). Limit: a listed indel is matched by length and sliding position, not by bases.
- 2026-11-25 16:59: the operator, on A2: "Fix once, re-time". One fix round (Sonnet): a candidate changes a call
  only by a tie the read's confident bases cannot break (depth only) or by a win over every class by the margin;
  otherwise the decision runs as on develop. Pure-indel candidates come from RJ-27's repeat only; RJ-24, RJ-23's
  long insertions and the contig edge stay as decided. Then a clean timing, interleaved with develop and the
  shared part. If A2 fails the speed floor after that, it is dropped.
- 2026-11-25 17:06: the MNV arm is in (Sonnet; not yet reviewed). Of 1,048 HG002 MNV rows, 932 have a run
  touching the block (58 are padded SNVs, 11 have no run), yet only 28 change in the rules path (AD down on 5,
  17 reads; RD down on 19, 109 reads; none up) and 21 in the judge (AD down on 7, 21 reads). The 11 no-run rows and
  the padded SNVs are unchanged; wall time unchanged. Departure from the brief: a left side with no run keeps one
  anchor base. It also grows the judge's guard scan for these rows, which scope item 5 keeps as on develop. 24 new
  existing-test failures, among them whole-carrier losses beside a germline deletion near the block ((0,40) and
  (40,0) where (40,40) is wanted), `--rescue-mnp` and RNA cases: an Opus review classifies them now.
- 2026-11-25 17:40: the review of the MNV arm (Opus). Binning invariance holds; no MNV row outside the class
  changes. Of the 24 test failures, 3 follow from the arm's rule and 21 from defects. The worst:
  1. with no right flank where only the left side grew, reads of another allele one base past the block count ALT
     (7:100371474 AD 2 to 240 with its neighbours left out; 359 such reads marked ALT in the full run);
  2. the telling check ignores the guard's masks, so a guard-masked telling base lets other-length reads count ALT;
  3. the exact-carrier rule's one forgiven shared-base substitution can fall on the telling base, so a read one
     unit longer ending there counts ALT (invariant 7);
  4. a background germline deletion at the anchor or telling base loses a whole haplotype ((0,40));
  5. the judge-guard change makes recurring indels next to the block stop being background;
  6. routing every read to the exact-carrier rule drops component carriers' partial evidence and breaks
     `--rescue-mnp` (RJ-25, RJ-28).
  The cause of most is my brief: it routed every read of such a row to the window rule. The fix round (Sonnet)
  keeps develop's MNP path for every read and adds the window as a gate on its REF and ALT calls, as A1's judge
  does; the judge's guard stays as on develop, and in the judge a read past the run with another length falls to
  the margin rule (A1's named limit). Defects 2 and 3 may sit in A1's rules half too: a check (Sonnet) runs them
  on delins and pure indels now.
- 2026-11-25 17:45: the check (Sonnet, synthetic, develop against A1, both paths). The guard-masked telling base
  is A1's defect too: when 12 reads carry an SNV at the telling base, the guard masks it, and 6 reads one T longer
  that end exactly there count ALT (AD 33 where 27 is right), for a delins and a pure insertion, in both paths; A1
  equals develop there. Without the mask they are not credited. The forgiven substitution on the telling base does
  not show for delins (MNP-specific). A fix round (Sonnet) makes a guard-masked telling base masked in both paths
  (the operator's "depth only" rule), as a delta for every option.
- 2026-11-25 17:52: B after its fix round (Sonnet). A read starting inside a left-hand run is placed by its bases,
  not its aligned start (REF credits that favour ALT 10 to 0; ALT lost through the cut 21 to 1); RJ-23's shared
  bases past the event are kept where develop reads them; only a read that stops short of the grown event's far
  end gets the relaxed hold; the masks stay. Every non-amendment test passes again; 9 fail, all B's named
  amendments (the run-window contract ×4, the "C>TA run10" spec cases ×5); `cargo test` 579 pass; binning
  invariance holds. HG002 rules path against the shared part: the wrong-allele set 417 rows (AD up on 82, down on
  0), the delins set 318 (AD up on 279, down on 0); the judge identical. Open: 46 new ALT credits (34 truth-row,
  12 wrong-row, of 1,644) still fit a substitution-only allele outside B's named limit, down from 66
  (16:1021001, 8 reads).
- 2026-11-25 17:59: A2 after its fix round (Sonnet). The counts are now faithful: no output difference outside the
  class; candidates enter only by a tie the read's confident bases cannot break or by a win by the margin; RJ-24,
  RJ-23's long insertions and the contig edge as decided; 1 existing Python test fails (the delins test shared
  with A1) and 11 Rust tests (6 count state graphs, 3 have carriers ending inside a run beside the change, 2 test
  helpers that compare against unscored candidates). HG002 judge against develop: the wrong-allele set 311 rows
  (AD down on 96, partial up on 84), the delins set 380 (AD down on 353); the rules path equals the shared part.
  Candidates: mean 5.4 and 4.7 per row, at most 47 (cap raised to 64, never reached). Speed, 5 interleaved runs
  under the build lock: 1.9x (7.31 s against 4.01 s) and 2.1x (3.16 against 1.50) develop's wall time, 2.2x and
  2.4x its memory. The time sits in the graph programme itself, so more states cost more time. Memory rises above
  2x and buys no time, so A2 fails the speed floor; as set at the operator's "Fix once, re-time", it is dropped
  from the measurement. Its build ($SP/a2j_v2_full.patch) and HG002 runs (d334/runs/a2j*) are kept.
- 2026-11-25 18:11: A1's telling-base fix (Sonnet) changes the design, so it stops for the operator (§5). The
  five wrong check cases now give (12, 27). The brief's premise did not hold: there the event already holds its
  run, so no telling base exists for the guard to mask. The fix adds a rule instead: at the first base past the
  event's run, where the guard masks it, a confident read that shows a base that goes on with the run is depth
  only, in both paths. Gates: 1 Python test fails (the delins test shared with A1), `cargo test` 577 pass,
  clippy, ruff and black clean; `.so` 1e30c02b, its own. HG002 against A1 before the fix (cmp_rows): 5, 2, 1
  and 1 rows move; the report's 10:37690930 did not move. Outside the class: on the wrong-allele set, 3 rows
  lose almost all ALT (8:145151666 187 to 2, 16:18865325 67 to 0, 16:28903765 19 to 0). Their reads span the
  run and carry the truth allele, which differs from the row's ALT at one base beside the deletion; develop
  counts them through RJ-23's guard (a change that recurs beside the event, tolerated on the ALT fragments).
  The rule overrides that only where the changed base goes on with a run, so it fixes part of another class.
  On the delins set, 8:121518866 (a truth row) loses 1 ALT and 3 REF reads: they show T or no base at the G
  past a 14-T run, where about 8 of 50 reads show T (noise at a long run's end, not a germline change); two of
  them span past it and their later bases are not the longer run. Open, a class of its own: RJ-23's guard
  credits a read whose truth allele differs from the row's ALT at a base beside the event (page later).
- 2026-11-25 18:31: the MNV arm after its fix round (Sonnet). Plain MNP reads take develop's `check_mnp` again;
  a gate over the block, each run beside it and the telling base decides held REF and ALT reads only (holds,
  depth only, or another allele); no substitution is forgiven on the run or the telling base; a germline
  deletion at the anchor or a telling base holds, as on develop. Gates: 4 Python tests fail (the 3 kept-base
  read-end tests the arm amends, and the delins test shared with A1), `cargo test` 578 pass, clippy, ruff and
  black clean; `.so` e6f4d995, its own. HG002 against A1 with the telling-base fix (cmp_rows): the rules path
  changes MNV rows only (wrong-allele set 20, AD down on 5; delins set 0). The judge changes 21 MNV rows (AD
  down on 7). By design, an MNV row's read that ends 1 to 3 bases past the block, short of the run's telling
  base, is depth only. Defect found in the check: the arm's judge files lack the telling-base fix's judge
  half, so its judge path equals A1 before the fix on 3 delins rows; the fix's final form goes in after the
  operator's decision on it.
- 2026-11-25 22:43, the operator's two answers:
  1. A1's run-next rule: "i would like to know number both ways before and after narrow". A1 as built stays
     (worktree and `.so` 1e30c02b unchanged); the narrowed rule is built beside it (Sonnet, label a1n): the rule
     fires unless the read's confident base at the first base that tells "one more unit of the run" from "the
     guard's change at the masked base" shows the change. Both are measured everywhere; the operator then
     chooses. The SNV and MNV arms take the chosen form after that.
  2. B: "we sould measure now", with the 46 unexplained credits as an open gap.
  Measurement now: develop, A1 as built and B on HG002 DNA on this Mac (the del+1 tree with the patient oracle,
  and the full row sets; Sonnet); the HPC batch for A1 and B is prepared on this Mac (Sonnet) and shown to the
  operator before it is submitted. HG002 RNA runs on HPC, where its BAM is, not on this Mac as §4 says.
- 2026-11-25 23:05: HG002 DNA, develop, A1 as built and B (Sonnet; `d334/m334/report.md`; a build is added with
  `m334/bin/m334_add.sh <label> <gbcms> <sha>`). Checks: develop's rerun equals the earlier oracle table in all
  1,957 cells and every earlier MAF in every field; each `.so` the same at start and end; 0 count-relation
  failures on 12,244 rows per build. On the 418 del+1 rows, each set alone, with the oracle (rules / judge):

  | Measure | develop | A1 as built | B |
  |:--|--:|--:|--:|
  | del+1 partners, false ALT reads | 1,032 / 1,189 | 811 / 706 | 1,222 / 1,189 |
  | of them "either" | 270 / 619 | 249 / 186 | 608 / 619 |
  | del+1 rows with ALT 2+ | 115 / 108 | 105 / 75 | 136 / 108 |
  | base partners, false ALT reads | 726 / 1,483 | 399 / 932 | 737 / 1,483 |
  | true calls, oracle-ALT reads missed (of 17,408) | 1,402 / 27 | 1,835 / 1,788 | 822 / 27 |
  | true calls, "either" reads credited | 797 / 2,368 | 711 / 733 | 1,373 / 2,368 |
  | true calls, REF reads ending inside the run | 562 / 2,055 | 403 / 424 | 1,024 / 2,055 |
  | gdelins hets (270), mean abs(VAF − 0.5) | 0.136 / 0.116 | 0.140 / 0.138 | 0.117 / 0.116 |

  B's judge equals develop's on every row. The bar "rows with no run unchanged" holds (0 of 11 such rows change)
  but tests little: with doublets counted, every delins row and tree row has a run. Not yet explained: A1's judge
  misses 1,761 more oracle-ALT reads than develop's (27 to 1,788).
- 2026-11-25 23:06: the HPC batch for B (o334b) and A1 as built (o334a1) is staged on this Mac and dry-run (Sonnet;
  `hpc/stage/batch_334opt.md`). Per option: a Linux build, the extension checked against dv6ff and the other
  options, then the pat334 arms with traces and the oracle (471 + 471 tasks and a sum), HG002 RNA rules (30
  tasks) and FORTE rules (#333's g333 set, 250 tasks and its aggregate). Each option writes its own run folders
  and reads the prep, chunks and pile-ups in place (no copy inside HPC): `pat334_arms.sbatch` and
  `pat334_oracle.sbatch` take an optional prep label, `g333_agg.py` a pile label; the old calls give the same
  output in the dry run. Only aggregates come back from patient data. An Opus trace of A1's extra missed ALT
  reads runs now. The batch waits for the operator.
- 2026-11-25 23:35: A1 narrowed (Sonnet; label a1n, `.so` dc87c364, its own; `$SP/a1n_delta.patch` against A1 as
  built). The run-next rule now fires unless the read's confident base at Q, the first base where "one more unit
  of the run" and "the guard's change at the masked base" differ, shows the change. Gates: 1 Python test fails (the
  delins test shared with A1), `cargo test`, clippy, ruff and black clean; 20 new tests in both paths (the five
  check cases still (12, 27); reads that show the change at Q (12, 33)). HG002 against A1 as built: 5 rows move,
  all among the 5 the rule moved, each back toward A1 before the rule (8:145151666 AD 2 to 185, develop 187;
  16:18865325 0 to 65, 67; 16:28903765 0 to 18, 19; 17:10555061 and 10:17204247 to A1 before the rule); the delins
  set is identical. The truth row 8:121518866 keeps its loss (RD 20 to 17): its reads' bases at Q are below Q20.
  Of 210 guarded sites, 155 have no Q inside the window, so the rule stays broad there. Staged for HPC as o334a1n
  (tarball fa2f268d).
- 2026-11-25 23:36: why A1's judge misses 1,761 more oracle-ALT reads (Opus). The oracle scores only the call and
  the two partners, so it is blind to a run one unit longer or shorter and to a change in a flank base. 1,330 of
  the reads are depth only in develop's own rules path too: the judge now holds the rules windows. By cause:
  1,124 end inside the grown run (all fit the longer run as well: A1 is right); 141 start inside the event (A1 is
  right, RJ-4); 249 read only the first flank base past the run (220 tell both run lengths; scope item 5's two
  flank bases drop them); 204 have the telling base below Q20 or N (200 tell both run lengths by later confident
  bases; the "Depth only" rule as written drops them); 20 short of RJ-23's shared bases; 13 flank bases at a site
  with no run (clarification 2, still open); 8 ALT bases behind a soft clip (a limit, `decide.rs:890-905`;
  develop's rules path drops them too); 2 the run-next rule (A1 is right: 14 T, not the ALT's 13). No group is
  a defect except the soft clip. The rules path's extra 433 have the same causes.
- 2026-11-25 23:41: HG002 DNA, A1 narrowed beside A1 as built (`m334/report.md`; `.so` dc87c364 at start and end;
  0 count-relation failures). Rules / judge:

  | Measure | develop | A1 as built | A1 narrowed | B |
  |:--|--:|--:|--:|--:|
  | del+1 partners, false ALT reads | 1,032 / 1,189 | 811 / 706 | 826 / 706 | 1,222 / 1,189 |
  | base partners, false ALT reads | 726 / 1,483 | 399 / 932 | 671 / 1,204 | 737 / 1,483 |
  | of them SIB (reads carrying the truth) | 521 / 1,144 | 239 / 747 | 511 / 1,019 | 520 / 1,144 |
  | true calls, oracle-ALT reads missed | 1,402 / 27 | 1,835 / 1,788 | 1,835 / 1,788 | 822 / 27 |
  | true calls, "either" reads credited | 797 / 2,368 | 711 / 733 | 711 / 733 | 1,373 / 2,368 |

  The narrowing changes nothing on the true calls. It gives back RJ-23's guard credits: reads that carry the truth
  and read past the masked base count ALT for a false allele one base from the truth (272 reads in each path on
  the base partners; the three wrong-allele rows of 23:35). So A1 as built has fewer false ALT reads, by
  overriding RJ-23 next to runs (the open class of 18:11); A1 narrowed keeps RJ-23 as decided.
- 2026-11-26 11:02: the operator: "Yes, all three". Batch 334opt on HPC (ledger 11:02 to 11:25): builds o334b (B),
  o334a1 (A1 as built), o334a1n (A1 narrowed), extensions f37c2e7f, e761c591, 95039bc4, all apart from dv6ff's
  f763393a. Rules path only for RNA (the judge refuses RNA in these builds).
- 2026-11-26 11:25: RNA results, against develop's dvR (same commit 6ffde565; 0 count-relation failures).
  HG002 RNA (`score_hg333.py`, 523,342 rows in every arm; develop reproduces FN 1,510 and FP2 5,077):

  | Change against develop | B | A1 as built | A1 narrowed |
  |:--|--:|--:|--:|
  | truth-ALT rows with ALT 0 at depth 10+ (FN) | −5 | +46 (INS +37, DEL +6, DELINS +3) | +18 (INS +10, DEL +6, DELINS +2) |
  | truth-REF rows ALT 1+ / 2+ (FP) | +4 / +1 | −31 / −23 | −9 / −4 |
  | ALT reads on truth-REF rows | +7 | −412 | −212 |
  | ALT reads on truth-ALT rows | +1,089 | −2,780 | −1,723 |

  FORTE (g333agg, 43 patients, aggregates): het carriers FN 1,762 in every option (develop 1,762), and the ref
  and sd sets unchanged. The wrong set (expect ALT 0): ALT reads at depth 10+ develop 25,733, B +31, A1 −221,
  A1n −187; rows ALT 1+ 7,312, B +10, A1 −45, A1n −24. The edge set (no truth) ALT reads B +738, A1 −200, A1n −43.
- 2026-11-26 11:25: the no-run bar on HG002 RNA (`$SP/rna_norun.py`, m334_report.py's run test on GRCh38): 0 of
  450 no-run rows change in every option. **Outside the class (stop rule):** SNV rows change, though no build has
  an SNV arm. A1 and A1n: 6 SNV rows (2 counts, 4 antisense depth only), all within 10 bases of a changed indel or
  complex row. B: 90 SNV rows (64 counts, 26 antisense depth only), 78 within 10 bases and 8 within 50; 4 are
  207 to 371 kb from any changed row (10:93702531 G>A and G>C, 11:117215716 T>A and T>C): antisense depth 1 and 3
  in develop and in A1, 0 in B, at total depth 0. BH q-values (asjd_qval) also move at 113 (B), 52 (A1), 9 (A1n)
  SNV rows, as invariant 1 allows when the row set's p-values change. Not explained yet; back to the operator.
- 2026-11-26 11:39: patients (pat334, 471 tasks per option, tumour and matched normal; the patient oracle; sums over
  roles and statuses, from `out/pat334_<label>/oracle/combined/`; the same 268 keys as develop's). Rules / judge:

  | Measure | develop | B | A1 as built | A1 narrowed |
  |:--|--:|--:|--:|--:|
  | IMPACT: del+1 partners, false ALT reads | 182 / 758 | 154 / 758 | 148 / 132 | 148 / 132 |
  | IMPACT: base partners, false ALT reads | 3,077 / 3,106 | 3,075 / 3,106 | 3,075 / 3,104 | 3,075 / 3,104 |
  | IMPACT: true calls, oracle-ALT reads missed (of 148,101) | 10,125 / 1,459 | 8,844 / 1,459 | 11,516 / 3,558 | 11,513 / 3,558 |
  | IMPACT: true calls, "either" reads credited | 1,156 / 2,598 | 1,174 / 2,598 | 1,134 / 1,214 | 1,134 / 1,214 |
  | ACCESS: del+1 partners, false ALT reads | 64 / 583 | 104 / 583 | 49 / 135 | 49 / 135 |
  | ACCESS: base partners, false ALT reads | 3,169 / 3,501 | 3,164 / 3,501 | 3,163 / 3,455 | 3,163 / 3,456 |
  | ACCESS: true calls, oracle-ALT reads missed (of 86,017) | 4,281 / 1,045 | 3,466 / 1,045 | 4,905 / 2,831 | 4,905 / 2,831 |
  | ACCESS: true calls, "either" reads credited | 1,065 / 1,843 | 1,212 / 1,843 | 1,018 / 1,141 | 1,018 / 1,141 |

  B's judge equals develop's. On patients the narrowing changes 4 reads in all (12 of 268 table rows): the base
  partners that HG002 DNA gave back (272 reads) are not in these rows. The oracle labels a read ALT by the call, the
  base partner and the del+1 partner only, so a read that ends inside the run can be oracle-ALT and still fit
  another run length; the HG002 trace (2026-11-25 23:36) found most of A1's extra misses of that kind.
- 2026-11-26 12:36: the SNV trace (the operator: "Yes, trace them"; job 14261698, `kit/trace334.sbatch`; public HG002
  RNA; 525 rows, every row within 500 bases of the 92 changed SNV rows; `$SP/trace334/*.log`). **Explained, inside the
  class by cause:**
  1. The subset runs give the full-chunk counts at all 92 SNV rows in every build, so the trace shows the real moves.
  2. The SNV's own call never changes. Every moved read (B 742; A1 and A1n 24 each) is one the multi-allelic REF
     guard (`sibling_claims_ref`) gives to a co-annotated delins row whose change lies where the SNV's REF is read.
     All 41 such rows are set `cis` (HG002's own SNV plus a run-length indel on one haplotype, written as one
     delins), in run category "grown". Only RD and antisense depth move at the SNV rows; AD never moves.
  3. B: the delins row's call goes from depth only (exact-carrier, uninformative) to ALT for 670 of these reads (B's
     "either" credit), so the guard takes them from the SNV's REF. A1 and A1n: 24 reads go from ALT to depth only at
     3 delins rows, so the SNV gets them back as REF (one row with counts, 2:105893996 C>T, RD 139 to 142).
  4. Truth (score_hg333's rows.tsv): the 92 changed SNV rows are 46 `wrong` and 46 `var`; 66 change RD (33 of each),
     the other 26 antisense depth only. AD never moves, so FN and FP do not change. At the 33 truth-ALT rows whose
     RD changes (all with a known zygosity), mean |VAF − expected|: het develop
     0.3962, B 0.3867 (13 rows closer, 1 farther), A1 and A1n 0.3963 (1 farther); hom develop 0.6958, B 0.6738
     (5 closer), A1 and A1n unchanged.
  5. The 4 "far" rows of 11:25 were an error in `rna_snvnear.py`: it measured distance to rows whose MAF counts
     change. Their siblings are 1 and 6 bases away (11:117215713 GTTT>GTA, 10:93702522 ACCGCCGCCG>ACCA); these
     siblings' own counts do not change, because the read is antisense or another sibling's ALT guard takes it.
  The trace does not tell A1 from A1n: both move the same 24 reads. Seen in passing, not an option effect: the
  guard checks the sibling alone, so 42 reads in B are ALT for the guard and demoted at the delins row itself
  (`sibling_claimed`); develop has the same asymmetry. Back to the operator for A1's form.
- 2026-11-26 14:25: the operator: "Ok lets go with A1 narrowed, but i think we should record this somehow so we dont
  need to do all these work again". A1's form is A1 narrowed (§3, A1 scope item 6). Recorded so it is not measured
  again: the register's Open list (`docs/reference/read-judgment.md`, this class), the memory note
  `.agents/memory/run-length-reads.md`, and the REF-guard effect as a known effect in
  `docs/harness/measurement-harness.md` ("Known effects"), and a no-run check that sorts each changed no-run row into
  qval-only, ref-guard or unexplained (`d334/bin/norun334.py`), so an SNV row that moves only through the REF guard
  is not a stop-rule item again. Checked: on HG002 RNA its ref-guard rows are the 92 traced rows exactly (B 90, A1 6,
  A1n 6; unexplained 0); on HG002 DNA's two sets, A1n against develop, 0 unexplained in both paths. The #334 tools
  and the earlier arm patches moved from the session scratchpad to `d334/bin/` and `d334/patches/`. The SNV and MNV arms are rebuilt on A1n
  (Sonnet, labels snvn and mnvn, new worktrees `334-snvn` and `334-mnvn`), then measured on HG002 DNA here and on
  HPC (pat334, HG002 RNA, FORTE) as the other options were.
- 2026-11-26 14:45: the MNV arm on A1 narrowed (Sonnet; label mnvn, worktree `334-mnvn`, `.so` 0769d6af, its own).
  The old arm was built on A1's 16:20 snapshot (its judge files equal that snapshot byte for byte); its own changes
  were merged onto A1n per file (`git merge-file`). The six conflicts in `carrier.rs` were all the run-next rule: A1n's
  narrowed rule kept, the arm's older copy dropped. One edit by hand: A1n's rule reads the allele the arm's gate names
  (`on_alt`, which is `result.is_alt` for every read the gate does not judge). A1n's judge files come in whole, so the
  18:31 defect is gone. Gates: `cargo test` 578 pass, clippy, ruff, black clean; pytest 4 expected failures (the
  delins test shared with A1, the 3 kept-base read-end tests the arm amends). A fifth failure was the arm's own new
  test (`test_mnv_run_window_contract.py`, the guard-masked telling base): it expected the judge to count 6 longer-run
  reads ALT, the 18:31 defect written down; corrected to (12, 27) in both paths, as A1n's own judge test of the same
  case (`test_run_window_judge.py`); not a develop test. HG002 against a1n: only MNV rows a run overlaps or touches
  change (wrong-allele set: rules 20, AD down on 5; judge 21, AD down on 7; delins set 0), as the old arm. The m334
  report (rules / judge): del+1 partners false ALT reads a1n 826 / 706, mnvn 811 / 664; del+1 rows ALT 2+ 105 / 75 to
  103 / 70; base partners, true-call misses and gdelins VAF unchanged; 0 count-relation failures. Staged as o334mnvn
  (tarball 65666369; its copy of the test file holds the old expectation, which no count reads). HPC build 14290659.
- 2026-11-26 14:46: the SNV arm on A1 narrowed (Sonnet; label snvn, worktree `334-snvn`, `.so` 121662dd, its own).
  The old arm was built on A1's shared part of 14:56 (`a1_shared_v4`, rules path only; blob ids equal in all four
  files); merged onto A1n per file. One conflict, in `event_for_windows`: A1n's `substituted`/`grows` test first, then
  the arm's SNV branch, then A1n's condition, so SNVs follow the arm, MNV and complex rows A1n, and padded
  one-substitution rows get no window as in A1n. A1n's run-next rule now also applies at the grown SNV windows. Gates:
  `cargo test` 578 pass, clippy, ruff, black clean; pytest 3 failures, all known: the delins test shared with A1, "SNVs
  have no window" (the arm's amendment) and the #409 strict xfail that passes by accident (an SNV in a run, as noted
  for the old arm on 2026-11-25). HG002 DNA: the wrong-allele and delins sets equal a1n in every field (they hold no
  SNV), so the m334 report rows equal a1n's. On the old arm's HG002 SNV set (423 rows), rules path: 347 rows change
  against a1n, all with a run holding the SNV's base (AD down on 309, up on 0; in-class AD 10,268 to 9,321, -9.2%; RD
  4,151 to 3,782, -8.9%); 0 SNVs with no run change; the judge equals a1n's. Against the old arm 199 rows differ (AD up
  on 88, down on 48): A1n's later fixes and its run-next rule at the SNV windows. Staged as o334snvn (tarball
  caa86c3b). HPC build 14291156.
- 2026-11-26 15:08: HG002 RNA, the two arms on A1 narrowed (HPC 14291416 mnvn and 14291441 snvn, 60 tasks COMPLETED
  0:0; rules path; scored against develop, 523,342 rows, `d334/score334.log`). FN, FP1, FP2, ALT reads at truth-REF
  rows and at truth-ALT rows, each against develop: B -5, +4, +1, +7, +1,089; A1 +46, -31, -23, -412, -2,780; A1n
  +18, -9, -4, -212, -1,723; SNVn +86, -321, -199, -4,858, -328,461; MNVn +30, -9, -4, -212, -7,677. SNVn moves only
  SNV rows past A1n: FN 810 to 878 (+68), FP1 -312, FP2 -195. MNVn moves only MNP rows past A1n: FN 11 to 23 (+12),
  ALT reads -5,951; the set has no truth-REF MNP row, so its FP cannot fall. No-run check
  (`d334/bin/norun334_rna_check2.log`, `norun334_rna_snvarm.log`): 0 unexplained rows in every arm; with `--snv-arm`
  SNVn changes 103,922 of 415,962 snv-run rows, and its 63 changed SNV rows with no run are q-value only. Speed (CPU,
  `d334/rna_times334.log`): MNVn 1.09x develop (0.98x A1n); SNVn 6,489 s, 4.67x develop (4.21x A1n), worst chunk c10
  8.56x A1n, a speed-floor item (not timed against the C++ tool); peak RSS 2.56 to 2.62 GB in every arm. The scorer's
  `vaf_med` is ALT/DP, so it falls when a read becomes depth only (SNV rows 0.610 to 0.583; MNP `dvaf_med` -0.027 to
  -0.071) whether or not the REF-to-ALT ratio moves; measure 4 uses the decided VAF, AD/(RD+AD)
  (`d334/bin/rna_vaf334.py`, log beside it; truth rows with DP 10+ in develop). Het rows beside a run (26,315), mean
  |VAF - 0.5|: develop 0.1613, B 0.1613, A1 0.1616, A1n 0.1614, SNVn 0.1630, MNVn 0.1615. Hom rows beside a run
  (15,914), mean |VAF - 1|: 0.0324, 0.0324, 0.0320, 0.0321, 0.0308, 0.0321. At het SNVs a run holds (20,920) SNVn
  takes RD down 7.8% and AD down 5.7% (pooled VAF 0.458 to 0.463), and mean |VAF - 0.5| goes from 0.1485 to 0.1505;
  at hom ones (13,156) from 0.0104 to 0.0088. No option moves a truth row that is not beside a run.
- 2026-11-26 16:05: the patient and FORTE batch for both arms on A1 narrowed (the operator: "Send it and run for
  both"; 16 jobs, all COMPLETED 0:0, no task log with an error; `d334/pat334_arms.log`, `d334/forte334_arms.log`,
  aggregates only). The record so far is committed on `docs/335-source-fixes` (f5b7c260, d590bd89; gate passed).
  Patients (pat334): the sets hold DELINS and MNV rows only, so the SNV arm equals A1n in every cell, and the patients
  cannot measure it. The MNV arm against A1n changes MNV calls only. Oracle-ALT reads credited ALT (rules / judge):
  IMPACT tumour somatic -3,335 / -8,778 (of 113,083 / 120,819), ACCESS tumour somatic -1,707 / -3,674 (of 58,535 /
  60,442), germline -66 to -178. "Either" reads credited ALT: -220 / -291 and -39 / -123. Oracle-REF reads credited
  ALT: 4 to 0 (IMPACT tumour somatic). Partner false ALT: at most 14 fewer (ACCESS tumour base, rules -16 against
  A1n's -2). FORTE: MNVn equals A1n in every stratum (both aggregate tables byte-identical). SNVn against develop
  (148 strata differ from A1n): truth-REF rows with ALT 1+ at DP 10+ 20,314 to 18,143 (-2,171), ALT 2+ 9,012 to
  7,897; het rows with ALT 1+ 98,989 to 97,943 (-1,046); ALT reads at hets -4.1%; wrong-allele rows with ALT 1+
  7,312 to 6,750 (A1n 7,288); no count-relation failure. FORTE CPU: develop 7,336 s, A1n 7,888 (1.08x), MNVn 8,160
  (1.11x), SNVn 26,598 (3.63x). All four options (B, A1n, SNVn, MNVn) are now measured on every data type the plan
  names; the operator chooses.
- 2026-11-26 16:47: the operator: "A1n". A1 narrowed is the rule for this class; B and the two arms are not taken.
  Before the register entry, the one existing test A1n changes (the 4 Rust geometry tests were approved on
  2026-11-25, "Approve all 4"): `tests/test_haplotype_judge_contract.py::test_reads_ending_inside_a_delins_count_
  alike_on_either_allele`. It asserts that the judge counts the 10 REF and 10 ALT reads ending 3 bases into a 24>6
  delins by the bases they hold: (rd, ad) = (20, 20). Under A1n: (10, 10), those 20 reads depth only. The test's
  "unique sequence" holds a GG doublet at the reads' end, a run by the decided definition. But the same case with a
  run-free stretch (`d334/delins_probe/probe.py`, built from a square-free word) gives develop rules (10, 10),
  develop judge (20, 20), A1n rules (10, 10), A1n judge (10, 10). So the change is not the class: it is A1n's judge
  holding the rules path's window at every site (clarification 2), the open choice of 2026-11-25 15:01 ("every
  site" or "only rows with a run"). On real data the two differ little: HG002, 0 of 9 rows with no run change in
  either path (13 reads at no-run sites in the review of 2026-11-25 16:10). Asked of the operator.
- 2026-11-26 18:42: the operator: "Every site" (the question named the test change to (10, 10)). The judge holds
  the rules path's windows at every complex site, also where no run lies, and the test's expectation changes with
  the code. Recorded as RJ-36 in `docs/reference/read-judgment.md` (rule, the options not taken, and an evidence
  paragraph with the four data types); the register's Open bullet for this class is replaced by the class that
  stays open, RJ-23's guard credit for a read that carries another allele beside the event (18:11). Next: land
  A1n's code on a feature branch from develop, with the test change and an outside review measured before any PR.

## 8. Cost log

| Date | Reviews | Conditions (added, removed) | HPC batches (tasks) | Agent tokens |
|:--|--:|:--|:--|--:|
| 2026-11-24 | 0 | 0, 0 | 0 | 156,990 (Fable, the survey) |
| 2026-11-24 | 1 | 0, 0 | 0 | 458,970 (Sonnet, the shared build: 346,426; Opus, its review: 112,544) |
| 2026-11-25 | 1 | 0, 0 | 0 | 140,760 (Sonnet, the narrow rebuild: 25,415; Opus, its review: 115,345) |
| 2026-11-25 | 1 | 0, 0 | 0 | 186,682 (Sonnet, the port: 24,177; Opus, the HG002 trace: 162,505) |
| 2026-11-25 | 0 | 1, 0 (the masked base past a grown run, operator) | 0 | 134,147 (Sonnet, two rounds of fixes: 45,265 and 88,882) |
| 2026-11-25 | 0 | 0, 0 | 0 | 1,497,913 (Sonnet, the four option builds: A1 judge 320,951; A2 judge 434,858; B 416,812; SNV arm 325,292) |
| 2026-11-25 | 1 | 0, 0 | 0 | 116,667 (Opus, the review of A1's judge half) |
| 2026-11-25 | 1 | 0, 0 | 0 | 112,981 (Opus, the review of the SNV arm) |
| 2026-11-25 | 0 | 0, 0 | 0 | 231,870 (Sonnet, the fix round of A1's judge half) |
| 2026-11-25 | 1 | 0, 0 | 0 | 173,336 (Opus, the review of B) |
| 2026-11-25 | 1 | 0, 0 | 0 | 180,110 (Opus, the review of A2's judge half) |
| 2026-11-25 | 0 | 0, 0 | 0 | 381,817 (Sonnet, the fix round of the SNV arm) |
| 2026-11-25 | 0 | 0, 0 | 0 | 419,423 (Sonnet, the MNV arm build) |
| 2026-11-25 | 1 | 0, 0 | 0 | 140,273 (Opus, the review of the MNV arm) |
| 2026-11-25 | 0 | 0, 0 | 0 | 145,078 (Sonnet, the telling-base check of A1) |
| 2026-11-25 | 0 | 0, 0 | 0 | 428,039 (Sonnet, the fix round of B) |
| 2026-11-25 | 0 | 0, 0 | 0 | 345,386 (Sonnet, the fix round and timing of A2) |
| 2026-11-25 | 0 | 0, 0 | 0 | 226,074 (Sonnet, the fix round of A1's telling base) |
| 2026-11-25 | 0 | 0, 0 | 0 | 388,178 (Sonnet, the fix round of the MNV arm) |
| 2026-11-25 | 0 | 0, 0 | 0 | 207,733 (Sonnet, the HG002 DNA measurement of develop, A1 and B) |
| 2026-11-25 | 0 | 0, 0 | 0 | 346,587 (Sonnet, the HPC batch for the options, staged and dry-run) |
| 2026-11-25 | 0 | 0, 0 | 0 | 304,598 (Sonnet, A1 narrowed) |
| 2026-11-25 | 1 | 0, 0 | 0 | 111,955 (Opus, the trace of A1's missed ALT reads) |
| 2026-11-26 | 0 | 0, 0 | 1 (batch 334opt: 3 builds; per option 471 arms, 471 oracle, 1 combine, 30 HG002 RNA, 250 FORTE, 1 aggregate: 3,672 tasks) | 0 |
| 2026-11-26 | 0 | 0, 0 | 1 (trace334: 4 tasks) | 0 |
| 2026-11-26 | 0 | 0, 0 | 0 | 197,252 (Sonnet, the MNV arm on A1 narrowed) |
| 2026-11-26 | 0 | 0, 0 | 0 | 176,819 (Sonnet, the SNV arm on A1 narrowed) |
