# #335 phase 2: shared scans, the judge's work and a narrower fetch (S2, S3, S5)

Status: approved 2026-11-16 (operator: "Approve as drafted": the options, the order and the budget; D2 decided (a), accept S5's log change; D3 open) · plan-first (`plan-first.md`); Stage 3 of the 6.8.0 order (`420-cycle-order.md`) ·
base develop 473e8920 (S1 #439 and S4 #440, #441 merged) · phase 1: `335-exact-speedups.md`
Rounds: 2 of 5 (one review per option; S2's, S3a's) · New count conditions: 0 (speed changes add none) · Builds: 3 of 9 (`fin473`, `s2`, `s3a`) ·
HPC batches: 0 of 3 · Reruns: 0 of 1 (phase 1's reserved rerun for M2's two IMPACT cells moves here)

## 1. The class

- **The work:** the time a run spends finding, decoding and scoring the same reads again, or reading
  records that no row needs. Three sources, each measured (phase 1, HPC M2 and M3, develop 1b11ffea):
  - **The scans (the rules path):** each row, and each part working for it, walks the whole read cache of
    its bin to find its few reads (the row loop 1 to 23% of CPU, the guard 1 to 18%).
  - **The fetch (both paths):** BAM fetch and decode is 26 to 82% of the rules path's CPU (78 to 82% on
    IMPACT) and 4 to 39% of the judge path's (most at 1,000 rows). S1 (`libdeflate`) cut part of it.
  - **The judge's programme:** one banded programme per state, per read, per row, with exact ties. At
    ACCESS 20,000 rows x 1 BAM it is 67% of all CPU: the fill 30%, the walks back through the cells 21%,
    the set-up per state 15% (M2 `perf`). Reads whose best score is 0 hold 87% of the scoring clock (M3).
- **The principle:** remove repeated or wasted work at its source. Never change what a read is called, an
  output field or a `--verbose` log line, unless the operator decides a log change (as for S4).
- **Invariants:** 1, binning invariance (S2 and S5 change which records a rule walks; each rule keeps its
  own test); 2, one quality contract (the judge's fast path gives the same scores, ties, bands and paths);
  3, output-aware (compute only what is read). 6 and 7 are not touched: no call changes.
- **Instance:** none (a class of cost, not a row).

## 2. What already touches it

- Phase 1 (`335-exact-speedups.md`): S1 and S4 merged, S1b dropped, M1 to M3 measured, the gate (18 arms,
  every output field, four whole-log arms).
- Research of 2026-11-16 (local, `scratchpad/speed2/plan2/`), on develop 473e8920:
  - the scans and the fetch (Opus; `s2_s5.md`; the cache and fetch code read in full);
  - the judge (Opus; `s3_judge.md`; `judge/` read in full; M2 `perf` and the M2 account);
  - the outside survey (Fable; `external2.md`; htslib, samtools, the C++ GBCMS, GATK, Octopus, vg, BWA,
    Edlib, Parasail, cgranges, COITrees, rust-lapper, rust-bio, rust-htslib, noodles; sources opened);
  - the local profiles (Sonnet): blocked, no result (section 4, M0).
- Outside, in short:
  - Octopus and BWA skip the programme when a read equals the haplotype; only the zero-mismatch case is
    exact (vg's 4-mismatch shortcut is not).
  - GATK's haplotype-prefix reuse is switched off in its main loop, and GATK tests it to 1e-9, not bit
    for bit. Parasail tests scores and whole trace tables against its scalar code; bwa-mem2 checks that
    its whole output is identical.
  - Only sorted starts with a running maximum end, and the cgranges-style tree (rust-bio's
    `ArrayBackedIntervalTree`), give a row's reads in cache order. COITrees, superintervals and
    rust-lapper do not: the first-mate rule would then pick another read.
  - htslib's multi-region iterator returns each record once, in file order, but it is not in
    rust-htslib's safe API, and for gaps under about 16 kb it reads the same compressed blocks as one
    wide fetch.
- Earlier claims, checked against today's code: `decided_states` gives each state `decided`'s score and
  band (proved); `placed_alone` holds everything the programmes read (proved); J2's cut 4 is **not**
  exact (the skipped nested build writes DEBUG lines); the 2026-11-14 survey's R4 (fetch the one-row
  regions) is **not** exact (finding F1, section 3.8).
- Tickets: #335 (J4, the judge's speed), #150 (deep-bin fetch), #312 (per-variant read slicing), #412 (SNV
  rows routed to the judge), #406 (memory). Rules: decision protocol §6, `speed-floor`, REJ-20261114-002.

## 3. Options (each its own build and PR, so each gain is its own)

| Option | Source it removes | Path | Gain (E, unless M) | Effort | Exact |
|:--|:--|:--|:--|:--|:--|
| **S2** shared scans, with **V** | the whole-cache walks; a CIGAR copy at each look; bases decoded before the position test | both | rules 25-40% dense, 2-6% sparse; judge 5-15% (V 5-9%); prototype -40 to -50% (M, older code) | medium | proved |
| **S3a** the judge's small items (P, N, F, M) | paths never read; the read-need test's repeats; the found-allele check's repeats | judge | 8-10% of the judge | low-medium | proved |
| **S3b** the certificate and the call on demand (C, D) | full programmes for the best-score-0 class | judge | 27-44% of the judge (central 37%) | high | if Theorem T holds in code (oracle) |
| **S3c** shared band, copy check and first columns (B) | per-state set-up repeated across a site's states | judge | 6-7% | medium | if the bands are equal (oracle) |
| **S5** a narrower fetch | records no row reads; bins that overlap | both | rules 35-60% at 1,000 rows, 10-25% at 5,000, 0-13% at 20,000 | medium | rows: proved; **log: changes** |
| **K** a SIMD fill (later) | the scalar fill that stays after C | judge | about 2% after C (9-15% without) | high | 32-bit lanes; 16-bit if n + C < 1,170 |

### 3.1 S2: shared scans (with V)
- **The source and the class:** 13 walks over the whole bin cache per row (the row loop, the
  `OBSERVED_ALLELE` scan, the tolerance scan, the judge's site filter, four guard passes, two junction
  walks, three RNA walks); 45 places that copy a CIGAR (`cigar()` copies; `cigar_cached()` borrows); 4
  places that decode a read's bases before its position test; `clip_to_fragment` decodes the CIGAR up
  to four times per record. The judge decodes a read's view and CIGAR again in each of its parts (V).
- **The change:** per bin, after the cache is built, two arrays: the largest read end so far and the
  smallest read start still to come. A rule with the interval [lo, hi) reads only the run from the first
  record whose running end passes lo to the last record whose remaining start is below hi, and keeps
  its own test inside the run. Each cached record decodes its CIGAR once; every helper borrows it. The
  guard tests the position first. V: one read view and one CIGAR per record per row for the judge.
- **The proof:** a record before the run ends at or before lo; a record after it starts at or after hi.
  So its span cannot meet [lo, hi), and each rule's own test needs the span to meet its interval (one
  lemma per rule, listed in `s2_s5.md` §2.4). The run is one slice of the cache, so the order is kept.
  The proof holds for any cache order. The CIGAR sums are the same integers in the same order.
- **Worked example (made-up numbers, to show the mechanism; M5 measures the real ones):** a bin of 4,000
  cached reads with 20 rows. Today each row looks at all 4,000 reads (80,000 looks), and each look decodes
  a CIGAR list. With S2 a row reads only the run that can meet its window, about 300 reads (about 6,000
  looks), and no look builds a list.
- **The oracle:** for each rule, the run's accepted records equal the whole-cache filter, index for index,
  on random clusters with clipped reads; `ref_end` from the raw CIGAR equals `cigar().end_pos()` on random
  CIGARs (op codes 0 to 8).
- **The log:** no `--verbose` change. In RNA, `--trace` shows fewer filter lines (`--trace` is not in the
  gate's log arms).

### 3.2 S3a: the judge's small proved items (P, N, F, M)
- **P:** a path is traced only where a reader reads it; one walk gives both the path and the steps (2-3%).
- **N:** the read-need test's exact cuts 1, 2, 3, 5 and 6: one `kept` call per read per `settle`; stop when
  the need is met or cannot be met (the test is monotone in the count); score the states together
  (`decided_states`, proved equal); score alike reads once, weighed by their number (4.5-6% before S3b).
- **F:** the found-allele check scores the given graphs once per read, not once per found allele (0.4% at
  ACCESS 20,000 x 1; more where found alleles are many, as on HG002).
- **M:** small items: the per-read clock and level check only when tracing; keys without per-read lists
  (about 1%).
- **Out:** J2's cut 4 (DEBUG lines go missing); returning `kept` as the final site (DEBUG lines go
  missing); reusing the site build's scores in the read call (another graph set, another band).

### 3.3 S3b: the certificate and the call on demand (C, D)
- **The source:** for every scored read, every state runs the full banded programme, the edge search and
  two tracebacks. Most state scores are read only through comparisons: a margin of 4, a floor, the first
  best (`decide.rs:336-384`, `476-481`, `1099-1103`).
- **The change:** C, a truncated programme (the certificate): per column, a bit set of the cells whose
  value is at or above a level L. At L = -6 no gap fits, so two bit sets per column are enough. It gives
  the exact score, end, edge flag and path of every state that scores 0 or -4: every state of a
  best-score-0 read at an SNV row. D: the call asks for a state's exact value only where a proved bound
  does not settle its decision; such a state runs the full programme, as today.
- **The proof:** Theorem T (`s3_judge.md` §5.2): the truncated recurrence gives every cell at or above L
  its exact value, so every state that scores at or above L gets the full programme's score, end, edge
  flag and path. Open obligations, as tests: O1, the certificate equals `run_reference` and `run_from` on
  random graphs and reads (exact, one mismatch, indels, soft ends, masked bases), every returned field;
  O2 and O3, the list of every reader of a state's score is complete, and each is a comparison that the
  bound settles.
- **Worked examples** (`s3_judge.md` §5.3, §5.5):
  - A read of 150 bases at an SNV row that spells ALT exactly. The row has four states: REF, ALT and the
    two other bases. Today: four full programmes, four edge searches, eight tracebacks. With C: four
    certificates (two bit sets per column). ALT scores 0; REF and the two other bases score -4 (one
    mismatch at the row's base). All four are exact, so no full programme runs.
  - A read at a one-base deletion row that spells the deletion. The deletion state scores 0 (certified).
    REF needs a gap, so it scores -7 or less: the certificate gives only that bound. `lost` needs the
    losing state's exact score and path, so REF still runs the full programme (or the deeper certificate
    at L = -8, §5.4).
  - A read whose best score is below 0 (a sequencing error): each state runs exactly where a bound does
    not settle its decision, so the worst case is today's work.
- **Measure first:** M4 (section 4) says how many state runs C covers and what it costs.

### 3.4 S3c: one band, one copy check and shared first columns (B)
- A site's states share the band, the check for copies off the band (Lemma K), and REF's first columns.
  Exact if the bands are equal array for array (oracle: shared against per-state, every field). 6-7%.

### 3.5 S5: a narrower fetch per bin
- **The source:** a bin fetches its whole span, so a sparse bin decodes many records no row reads (8.5
  times what rows read at 1,000 rows, older Mac measure). Neighbouring bins overlap by about 450 bases,
  and by up to a window after a cap split (a record decoded up to 5 times in dense hotspots).
- **The change:** per bin, the merged union of the members' rule intervals (S2's functions, cut to
  today's bin), fetched one region at a time with rust-htslib's `fetch`. A record the region before
  already returned is dropped. If records arrive out of order, the bin falls back to the wide fetch.
- **The proof:** with S2's lemma each rule reads only records that meet its own interval. Those records are
  all in the union and arrive in file order, so each rule reads the same records in the same order. The
  regions are the rules' own intervals, not the one-row read windows (F1 shows those are not enough).
- **The log changes (D2, decided 2026-11-16: (a); the operator: "accept the S5 log change (a)"):**
  1. the DEBUG line per bin, "Bin tid=... : N reads cached (M filtered ...)": the counts of the regions;
  2. the WARN lines per BAM that count records with no bases or no qualities: records outside every row's
     reach are no longer counted. Today these counts depend on the bin geometry; with S5 they do not;
  3. `--trace`: fewer filter and clip lines;
  4. a corrupt record outside every row's reach no longer stops the run (develop stops with an error).
- **Worked example:** two rows 6,000 bases apart share one bin (`BIN_WINDOW` is 10,000, `engine.rs:108`),
  so today the bin decodes every record of those 6,000 bases and 5,000 more past the last row. With S5 it
  decodes only the records that overlap each row's rule intervals: the read window of an SNV row (its base
  and 5 more each side), up to 300 bases each side for a row the exact-carrier rule judges.
- **Not used:** htslib's multi-region iterator (the repo's first `unsafe` code; CRAM not proved).

### 3.6 K: a SIMD fill (decided after S3b is measured)
- One read per lane, the same graph for all lanes, scalar traceback. 32-bit lanes are exact by
  construction (integer add and max do not round); 16-bit lanes are exact while n + C < 1,170, checked
  per read. After C only about 13% of the scoring is left, so K may not be worth its risk.

### 3.7 The floor arithmetic (judge path, CPU per read-row pair; `s3_judge.md` §10)

| Cell | Floor (the C++ tool) | Today (M2, before S1) | After S2, S3a-c (central) | Against the floor |
|:--|--:|--:|--:|:--|
| ACCESS 20,000 x 1 | 12.2 us | 65.0 us | 24.6 us (19.8 to 35) | about 2.0 times over |
| ACCESS 5,000 x 9 | 21.4 us | 68.8 us | 26.5 us | about 1.24 times over |
| ACCESS 5,000 x 1 | 32.5 us | 70.2 us | 27.6 us | meets (0.85) |

- Exact options give about 2.6 times on CPU at ACCESS 20,000 x 1 (range 1.9 to 3.3); that cell needs 5.3.
  Five of the eight production cells would meet the floor; ACCESS 5,000 x 9 and IMPACT 20,000 x 1 would
  be about 1.2 times over; ACCESS 20,000 x 1 about 2 times over.
- To close the rest, every one of these is needed, each exact but none yet measured: the certificate at a
  tenth of the scalar cost on 95% or more of state runs; no exact work for the other-base states unless a
  decision needs it; the read call's glue at about 2 us; the site build at about 1.5 us per pair; S2 and S5
  on the shared part. Even then the cell is near the line. D4 (after M4) is yours.

### 3.8 Two findings outside #335 (each for its own ticket, by the decision protocol)
- **F1, a gap in invariant 1:** the long-deletion junction judge's scan for other alleles
  (`junction.rs:1140-1141`) reads reads that start up to `mh + 3 - pad` bases past the row's read window.
  Production bins reach 5,000 bases past every row (`engine.rs:203-205`) and no CLI option sets the
  window, so CLI runs are not affected. The one-row geometry (window 1, cap 1) misses those reads. Found
  from the code; not yet reproduced with reads.
- **F2, a cache that can differ between runs:** the verdict cache's key (`decide.rs:1309-1314`) leaves out
  qualities, MAPQ, tags and mate fields. A record duplicated byte for byte except those gets the first
  copy's call only while the thread's budget has room (`judge/mod.rs:292-301`), and the room depends on
  which bins ran on that thread before. Rare (a normal BAM has no such duplicate). Found from the code;
  not reproduced.

**Out (not exact; on record):** caps or samples of reads; the 4-mismatch gapless shortcut; adaptive bands,
X-drop, Z-drop; a band cap without proof; 8- or 16-bit cells without a bound; tolerance testing; a pruned
score that a reader reads later; interval indexes that lose the cache order; J2's cut 4; reusing
site-build scores; WFA or bit-parallel distance as the scorer; any floating-point reordering.

## 4. Measurement, fixed before the runs

- **M0, the reference files (blocked 2026-11-16):** from this session, reads of `~/Downloads/genome.fa`
  and of the OneDrive copy of `Homo_sapiens_assembly19.fasta` hang. The gate reads them too. Fix: copy
  both, with their `.fai`, to `~/test/gbcms/refs/`, check the sizes, and point the manifests there.
- **M5, no build:** local profiles on develop 473e8920 (HG002 `sv` and `pri`, the panel; both paths;
  shares only); the counts from the `--verbose` bin lines: records fetched per bin against records in the
  union of its rows' reaches (S5), looks against reads used (S2).
- **M4, a counter build (never merged):** per programme run: the state's class, the band step, the band's
  rows, the score (the share at 0 and -4, at -7 and -8, and below); per row the site type and states;
  `build_core`'s parts. Locally on HG002 and the panel; on HPC at ACCESS 20,000 x 1 and 5,000 x 9, with
  the reserved rerun for M2's two IMPACT cells on Gold 6348 nodes (your yes first).
- **The gate for every option (bar: 0 differences):** phase 1's 18 arms against a new baseline (develop
  473e8920), every output field of every row, the four whole-log arms; plus the option's oracle tests,
  `cargo test` and `pytest`. S5's log arms are compared under D2's rule: every line identical, except
  (1) the per-bin DEBUG line "Bin tid=T S-E: N reads cached (M filtered: ...)", which keeps its bin (tid,
  start, end) and may change its counts, and (2) the two per-BAM WARN lines "skipped N record(s) stored
  without bases (SEQ '*')" and "... without base qualities (QUAL '*')", which may change N or be absent.
  The gate gets this mask with the S5 build. `--trace` is not in the log arms. The fourth change (a
  corrupt record outside every row's reach no longer stops the run) gets its own test with a made-up
  BAM. The CHANGELOG names all four.
- **The gain (HPC, your yes first):** pf77 times S2, S3a and S5 against develop; pf78 times S3b and S3c.
  The rules of pf75: the 8 production cells, both paths, 3 repeats, order rotated; CPU is the main number.
  An option is kept if the gate holds and it cuts CPU on a path; a cell worse by more than 3% of CPU
  comes back to you. Peak memory is reported beside each cell. It may rise when time falls and every
  output is identical, up to 2 times develop's peak before phase 2 (`fin473`) on every timed cell, for
  the phase as a whole; a cell above 1.5 times is named in the report, and one above 2 times comes back
  to you (operator, 2026-11-16: "I am ok with memory increase as long as we can reduce time and
  maintain accuracy", "Not unlimited increase", "I am ok with 1.5 to 2x increase"; this replaces the
  earlier 10% limit). S2 and S5 reach RNA too, so pf77 holds the FORTE RNA cells (rules path).
- **What this cannot see:** whether ACCESS 20,000 x 1 reaches the floor before S3b is built (M4 says how
  much S3b can cover); timings on this Mac are not speed results.

## 5. Budget and stop rule

- **Builds: 9 at most:** the develop baseline (`fin473`), M4's counter build, S2, S3a, S3b, S3c, S5, one
  combined, and one spare for a fix (phase 1 needed two).
- **Reviews: 5 at most:** one adversarial review (Fable) of each option's diff against its proof.
- **HPC: 3 batches** (M4, pf77, pf78) plus one rerun.
- An option with any output difference is dropped, not mended by accepting the difference. Past any
  limit, stop and bring it back to you. K, and any option not listed here, needs its own page.

## 6. Work and models

- Opus: each option's code design and its proof (S3b first). Sonnet: the builds, the gate arms, the
  timings, the profiles. Fable: the adversarial review of each diff. Haiku: mechanical scans (every
  `cigar()` site, every whole-cache walk).
- Every brief carries the PHI rules in full, also for a code-only agent (2026-11-16: a code-only agent
  listed the HPC copy-back folders to find M2's files; no identifier was printed).
- **Order:** M0 and M5 (no build) → S2 with V → S3a → M4 → S3b → S3c → S5 (after D2) → pf77 and pf78 →
  the K decision and D4. Each kept option lands as its own PR to develop with its gate evidence.

## 7. State

- 2026-11-16: research done (four agents; the profiles blocked by M0); the page drafted. M0 cleared at
  10:14 (ERR-20261116-001); M5 restarted at 10:15. The S3b worked example corrected (four states at an SNV
  row, all certified; the deletion-row case added). Approved at 11:14 as drafted; D2 (S5's log) and D3
  (F1, F2 tickets) open. S2 with V started (branch `feature/335-s2-shared-scans`, from 473e8920).
  M5 done (11:20; local shares, Mac, 4 threads): on the judge path the site build is 45-58% of CPU on
  HG002 and a 35-row panel (about 14% at ACCESS on HPC), and the read-need test 5% (HG002 primary), 31%
  (HG002 SV) and 57% (the panel), so N's gain depends on the data type; the site build outside the
  read-need test (40% at HG002 primary) has no option yet, and M4's `build_core` profile is the first
  step. Rules path: S2's scans 14-19%, fetch and decode 3-21%. The counts half of M5 (records fetched
  against records used) waits for S5. M6 (pf75's table, no new run): one sample, as Nextflow runs it
  (one BAM, 4 CPUs), keeps its CPUs 86-98% busy on the judge path and 53-94% on the rules path, so a CPU
  cut is a wall-time cut on the judge path.
  11:50: D2 decided, (a): S5 may change the four log outputs of §3.5; the gate's rule for them is in §4.
- 2026-11-16, S2 with V (branch `feature/335-s2-shared-scans`; Opus wrote it, Fable reviewed it):
  - The code: a bin's `ReadCache` keeps each record's decoded CIGAR and the two bounds; 13 walks read
    only their run; 45 CIGAR copies borrow the kept list; the judge's call takes the view its site was
    built from. A record whose CIGAR holds an op code 9-15 is not decoded at the cache build and lies in
    every run, so a corrupt file stops at the same record, line and message as develop (found in our own
    check: the first version stopped earlier, before develop's DEBUG lines, and at another panic site).
  - The review (round 1, Fable): no input that breaks exactness; test gaps (narrowed intervals that
    passed, five walks without an old-against-new test, a stale copy of each rule's test) and one memory
    point (views kept for rows with no site), all fixed. Of 35 mutations, 29 fail a test and 6 cannot
    change what a rule reads. The production changes after the review (two extracted functions, shared
    predicates, the two view drops) were read line by line: each keeps develop's order of tests.
  - The baseline `fin473` (develop 473e8920, build 1): 0 rows differ from `fin682` on all 18 arms; the
    only log difference is S4's approved one on `logs_rnafill` (350 and 70 lines, as on 2026-11-15).
  - The gate on `s2` (build 2: 473e8920 plus the patch, sha256 7b29efad8b5bfac1, test lines added after
    it only): **0 differences on every arm**, every row and the four whole logs (`logs_rnafill` 6,900
    lines). Local tests: clippy and cargo test (563 passed), pytest (2,769 passed), ruff, black, mypy, all
    rc=0. The gain is pf77's (HPC; your yes first).
- 2026-11-16, S3a (branch `feature/335-s3a-judge-cuts`, from develop dcba8a94; Opus wrote it, Fable
  reviewed it):
  - The code: P (one walk back gives a state's path and its steps); N's cuts 1, 2, 3, 5 and 6 in the
    read-need test (each group's call on `kept` once per `settle`; the count stops once its answer is
    fixed; a read's states scored together; reads alike in all the scorer reads scored once, weighed by
    their number); F (the found-allele check scores the given graphs once per read); M (the read's key
    built once, the clock read only when tracing). Left out: tracing only the paths a reader reads (the
    budget charges each path's capacity), two M items that save nothing, and cut 4.
  - The review (round 2, Fable): no input that changes an output, a `--verbose` or `--trace` line, an
    error or a budget choice. Two test gaps (nothing pinned the read's start in N's grouping key; nothing
    kept the read-need test off the thread's pool) closed with tests only. Mutations: 24 of the author's
    25 fail a test (the 25th is the control the proof says changes nothing), and 7 of the review's 8
    (the 8th moves only the trace clock).
  - The gate on `s3a` (build 3: dcba8a94 plus the patch, sha256 283fb6c559e5b8c4, test lines added after
    it only) against `s2`: **0 differences on every arm**, the four whole logs included.
  - Memory (HG002 SV judge arm, 15,054 rows, 4 threads, one run each on this Mac): CPU 1,852 to 1,631 s
    (0.88), wall 1,001 to 908 s (0.91), peak memory 1,406 to 1,924 MB (1.37); the 13,473 data rows
    identical. The rise is N's one score table per state during `settle` (develop reused one); the budget
    counts neither. Within the operator's limit above (at most 2 times; 2026-11-16).

## 8. Cost log

- Research agents (2026-11-16): the scans and the fetch (Opus) 0.40 M tokens; the judge (Opus) 0.60 M; the
  outside survey (Fable) 0.22 M; the profiles (Sonnet) 0.16 M, blocked.
- Slips: the profile agent ran `ls -lO`, `xattr` and `grep` on Downloads before it understood the hang
  (each hung and was killed). The judge agent listed `hpc/out/` and `hpc/out/m1profJ/` with `ls` and
  `find` (no sample ID or BAM name in its transcript or report, checked by pattern); my brief did not
  carry the copy-back rule (LRN-20261116-002).
