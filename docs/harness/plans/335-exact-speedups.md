# #335: exact speed-ups across the whole run, both paths

Status: approved 2026-11-14 (operator: "Approve the page") · plan-first (`plan-first.md`); Stage 3 of the
6.8.0 order (`420-cycle-order.md`) · base develop 1b11ffea (#438 merged)
Rounds: 0 of 2 · New count conditions: 0 (speed changes add none) · Builds: 8 of 8 (`src_m1prof`, `src_s1`, `src_s1b`, `src_s4`, `src_s1c`, `src_s1d`, `src_s4b`, `src_s4c`) · HPC batches: 3 of 3 (M2, pf75 and pf76 done; pf76, RNA, added by the operator) · Reruns: 0 of 1

## 1. The class

- **The work:** the time a run spends, on both paths, on work that is repeated or wasted. Examples:
  each row looks through its bin's whole read cache to find its own few reads, a read's CIGAR is
  decoded into a new list at each look, the judge scores the same read against the same graph more
  than once, a file is read again for each BAM.
- **The principle:** remove repeated or wasted work; never change what any read is called. Speed is
  never bought with accuracy (operator, 2026-11-14; decision protocol §6; REJ-20261114-002). A change
  is kept only if every output field of every row is identical on every set.
- **Invariants it touches:** 1, binning invariance (a cache slice or a narrower fetch must hold every
  read a one-row bin holds); 2, one quality contract (a faster judge kernel keeps the same gate and
  the same ties); 3, the output-aware engine (the insert size is computed when mFSD is off, then
  dropped); 6, gated columns. Invariant 7 is not touched: no call changes.
- **The size of the gap (pf73, HPC, develop 3be1765b):** the default path takes 0.39 to 0.94 of the
  C++ tool's time; the judge 0.76 to 5.27. At ACCESS 20,000 rows × 1 BAM the judge's extra time over
  the default path is 200 s and must fall to 3 s (a 5.3 times cut of the whole judge). pf74 (today):
  the read-need change costs the judge 0.97 to 1.08 of its earlier wall time.
- **Where the time goes (deep ACCESS, 20,000 rows; HPC `perf` and timers, 2026-10-31, older code):**
  the default path: BAM fetch and decode 27%, the row loop over the cache 23%, the `OBSERVED_ALLELE`
  scan 21%, the guard 18%, the read call 6% (the CIGAR lists 30% and malloc and free 13% of all
  CPU). The judge: the read call 71%, the site build 11%, the shared work 15% (`graph::run_from`
  31%). **Today (this Mac, develop's code; sampling profiles, shares only):** on HG002 `sv`,
  `run_from` is 81% of the judge's CPU and the site build 45%; on the IMPACT panel, the site build is
  58% of the judge's CPU (`run_from` 51% inside it: the read-need test and the found-allele check
  score reads there) and the read call 40%. On the default path at `sv`, 67% of the CPU is judge
  work, for the SNV rows routed to the judge (beside a listed indel or delins). So judge speed-ups
  also speed the default path, and the site build is the judge's largest part on panel data.
- **Instance:** none (a class of cost, not a row).

## 2. What already touches it

- Tickets: #335 (judge speed; the four levers of 2026-10-31), #150 (deep-bin fetch; the profile of
  2026-10-31), #312 (per-variant read slicing; protoA), #412 (SNV rows routed to the judge), #406
  (memory; mask copies).
- Prototypes, never merged (local `plan670/perf/`): one cache slice per row plus no CIGAR list per
  look (protoA: HG002 `sv`, default path, -40% at 1 thread and -50% at 8, MAF identical; Mac, older
  code); a faster allocator (protoC, mimalloc: 6 to 17% more, peak memory +12 to 21%); work split
  inside a bin (protoB: 0 to 6%); a lazy `polars` import (about 0.2 s a run).
- Merged and exact: the judge's alike-read caches, the column-wise graph programme, shared first
  columns (#304, 2026-10-28; no timing of each).
- Rules and records: decision protocol §6 (exact speed-ups only); `speed-floor` memory;
  REJ-20261114-002 (no cap or sample for speed); LRN-20261114-001.
- Research of 2026-11-14 (local, `scratchpad/speed2/`): the inventory (Sonnet), the code survey
  (Opus; 30 candidates with `file:line` and proofs), the outside survey (Fable; 8 exact techniques
  with sources), the profiles (Sonnet; HG002 and the IMPACT panel, both paths). Checked against the
  code: the Rust log level is opened to TRACE (`lib.rs:50`); the insert size is computed for every
  read (`engine.rs:1999`); the GTF and editing sites load inside each counting call
  (`engine.rs:504-540`); `rust-htslib` 1.0.1 has a `libdeflate` feature, off here; `cigar()` clones
  even a cached CIGAR, `cigar_cached()` borrows.
- Corrected: the 2026-11-12 survey's option 6 (reuse the site build's read scores in the read call)
  is not exact: a read's band can differ between the two (`haplotypes.rs:3622-3624`, `3897-3905`).

## 3. Options (each a separate measurement build, so each gain is its own)

| Option | What it removes | Path | Expected CPU cut | Effort | Proof that it is exact |
|:--|:--|:--|:--|:--|:--|
| **S1. Build settings and libraries** | LTO and one code unit (`[profile.release]`); `libdeflate` for BAM blocks; the log level set from Python's | both | 5 to 15%; `libdeflate` 8 to 12% of the default path at 20,000 rows, more at 1,000 (E) | low | no source change in counting; Rust does not reorder or contract floating-point work; inflate is lossless and checked by CRC |
| **S1b. A faster allocator** (mimalloc) | heap time (the Mac's `free` is slow) | both | 6 to 17% (M, protoC, default path, older code) | low | no source change; costs peak memory (+12 to 21%, M): your call |
| **S2. Shared scans** | each row's look through the whole bin cache (one slice per row for every consumer); a read's end read from the raw CIGAR, start tested first; guard loops test the position before they decode bases; one CIGAR decode per record, borrowed | both, mostly the default path | the default path 40 to 65% on deep cells (E; protoA -40 to -50%, M); the judge about 10% (E) | medium | the slice drops only records each consumer's own test rejects, in the same order; a binning-invariance arm; an old-against-new unit oracle on random reads |
| **S3. Judge work done once** | the read-need test's six cuts (the call with every item in place once per read; stop when the answer is fixed; and others); the found-allele check scored once per read; alike reads keyed without their names; one backward walk in place of three; bands reused | the judge, and the default path's routed rows | 10 to 30% of the judge (E), up to the site build's share (45 to 58% on HG002 and the panel, M); measured first (M1, M2) | medium | pure functions called once; tests monotone in the count; the reference programme (`run_reference`) as an oracle on random graphs |
| **S4. Python, RNA and small items** | the BAM properties step (a Python loop over 20,000 records per BAM); the GTF and editing sites read once per run, not per BAM and again in the MNP rescue; no insert size when mFSD is off; the lazy `polars` import | both; RNA | 0.3 to 0.6 s per BAM; seconds per BAM in RNA (E) | low | the same tallies and the same index; mFSD's columns are absent when off, so the insert size never reaches an output |
| **S5. A narrower fetch per bin** | reads fetched for gaps between a bin's rows | both | the default path 30 to 55% at 1,000 rows, 10 to 17% at 5,000, about 0 at 20,000 (E) | medium | binning invariance: every row's reads are still fetched once, in file order |

(E) is an estimate from older profiles; (M) is measured. No option has been built or timed on
today's code.

**Later, in its own plan (high effort; the judge's floor needs them):** an integer SIMD kernel for
`run_from` (32-bit lanes are exact; 16-bit lanes only with a proved bound), and a fast path for
reads that spell one state exactly (its paths, ties and ends must first be proved the same). The
exact options above sum to about 2 to 3 times on the judge (E); the floor needs 5.3 times at the
worst cell. M3 tells whether the fast path can carry the rest.

**Out (not exact; on record):** any cap or sample of reads; early exits that skip a read whose result
is read later; any floating-point reordering (a SIMD PairHMM in `bio`, other sums in mFSD or Fisher,
`f32`, fast-math); WFA or bit-parallel distance as the scorer; a CPU target above the wheel's
baseline; `panic = "abort"`; the htslib decode threads (they compete with the bin threads).

## 4. Measurement, fixed before the runs

- **M1, before any option (measurement only, never merged):** fix the profile kit (`tp_summary.py`
  reads no stack from this Xcode; the profile agent's `tp2.py` does); a phase-timer build on develop
  (the #150 timers, plus named frames for the read-need test, the found-allele check and the routed
  rows), so S3's shares are measured, not estimated.
- **M2, one HPC profile batch (your yes first):** the 8 production cells, both paths, `perf` on the
  node, with the phase-timer build: the shares on today's code at production size.
- **M3, the read classes:** on the ACCESS and IMPACT production cells, the share of the judge's
  scoring time on reads answered by a cache, reads with best score 0, and the rest (the survey's
  deciding number for the fast path).
- **The gate for every option (bar: 0 differences):** every output field of every row equal to
  develop 1b11ffea on HG002's 7 sets (both paths), the IMPACT tier, ACCESS, WES and FORTE RNA (GTF on
  and off), a 9- or 10-BAM fill-out, mFSD on (with Parquet), the MNP rescue on, `gbcms merge`, and
  the binning-invariance arms (one row per bin against the default bins). Plus `cargo test`,
  `pytest`, and an old-against-new oracle for each changed function. From 2026-11-15 the gate also compares whole
  logs (three log arms, `--verbose`: HG002 and the fill-out on both paths, the FORTE RNA runs with the
  GTF; times, durations, paths and source lines masked), as a log line is output too.
- **The gain, on HPC (one batch, pf75; your yes first):** each kept option against develop in the
  same task, on the 8 production cells, both paths, 3 repeats, order rotated. CPU time is the main
  number (pf74: wall time moves up to 15% node to node; CPU 5%); wall time and peak memory beside it.
  An option is kept if the gate holds and it cuts CPU on a path; any cell worse by more than 3% CPU,
  or peak memory up by more than 10%, comes back to you.
- **Every data type the option reaches (operator, 2026-11-15):** an option is timed on each data type
  its change reaches before it is kept or dropped. pf75's cells are DNA only, so S4's RNA part (the GTF
  and editing sites read once per run) and S1 on RNA get their own batch, pf76: RNA with the GTF, one BAM
  and a fill-out, develop against S1 and S4, 3 repeats, the same rules. Its steps (each needs the
  operator's yes for the HPC part):
  1. R1, read only: the labels free; the GTF (`inputs/refs/grch38/Homo_sapiens.GRCh38.111.gtf`, 1.46 GB)
     present; the count and sizes (by `stat`) of the FORTE STAR BAMs; the queue; space.
  2. A copy job: 10 whole FORTE STAR BAMs, read one at a time, to `inputs/prod/rna/b01.bam` to `b10.bam`
     (`select_copy_prod.sbatch` with an `rna` type: largest first after the top 5% are set aside, as for
     IMPACT; the map to the sources stays in `inputs/prod/.private`).
  3. A rows job: positions in GTF exons that `b01` covers, PHI-free lists of 1,000, 5,000 and 20,000 rows
     (`build_prod_rows.sbatch` and `make_prod_rows.py` with the GRCh38 reference).
  4. The S1 build again, as fixed for the log lines (S1 with the log cap held to each counting call).
  5. Two arrays: RNA 1 BAM x 1,000, 5,000 and 20,000 rows, and RNA 10 BAMs x 5,000 rows, with the GTF, on the
     rules path only (`gbcms rna` refuses the sequence judge: it does not judge spliced reads yet); plus DNA ACCESS 5,000 x 9 (S4's memory flag, and the fixed S1 on DNA). Builds develop, the
     fixed S1 and S4 in each task, rotated by repeat. `run_prod.sbatch` gains the tool `rnarules`
     (`gbcms rna ... --gtf`).
  6. The copy-back of the timing lines only, and `agg_pf76.py`.
- **What this cannot see:** whether exact options reach the floor (likely not on the judge's worst
  cells: M3 says); timings on this Mac are not speed results (load 5 to 19 today).

## 5. Budget and stop rule

Builds: M1's phase-timer build, one per option (S1, S1b, S2, S3, S4, S5) and one combined: 8 at most.
Review rounds: 2 at most (Fable, adversarial, on exactness: each option's diff against its proof).
HPC: 3 batches (M2, pf75 and pf76; the operator added pf76 on 2026-11-15) plus one rerun. An option with any output difference is dropped, not
mended by accepting the difference. Past any limit, stop and bring it back to you.

## 6. Work and models

Opus: each option's design and its proof. Sonnet: the builds, the gate arms and the timings. Fable:
the adversarial review of each diff. Haiku: mechanical scans (every `cigar()` call site). Scripts:
the batch680 compare, the pf aggregators. Each kept option lands as its own PR to develop, in the
order of measured gain over effort (S1 and S4 first, then S2, S3, S5), each with its gate evidence.

## 7. State

- 2026-11-14: research done; the page approved. Four builds frozen, each develop 1b11ffea plus its patch.
- M1 (`src_m1prof`): outputs equal to develop's with the timers off and on (HG002, both paths, 0 cells
  differ). Shares (this Mac): the read-need test is 33% of the judge's account on `sv` and 56% on the
  IMPACT panel; on `sv`'s default path, the 8.4% of rows routed to the judge hold 70% of it. M3: the
  cache class is near 0 on HG002 and the panel; best score 0 is 25 to 76% of the scoring clock.
- S1, S1b, S4: 0 cells differ from develop on HG002 (both paths; S4 also with mFSD and Parquet; S1's
  trace identical). S4's properties step saves only about 0.03 s per BAM.
- The gate (`batch680/gate335.py`, 14 arms): develop's baseline passed every arm; a repeat equals it;
  one row per bin equals the default bins on both paths.
- 2026-11-15, the gate on the options: **S1 and S4 passed** (16 comparisons, 0 rows differ, no run failed).
  **S1b did not pass:** 0 rows differ, but one judge run (an IMPACT tier run, 11 rows) ended with a
  segmentation fault (exit -11) after it had written its output and its QC summary. Twenty reruns of that run
  passed on S1b, develop and S1: the crash is rare (1 in about 3,500 S1b runs; 0 in about 3,500 runs each of
  develop and S1). S1b is S1 plus mimalloc, so the likely cause is the allocator at process exit. S1b is out of
  pf75; a crash is a failed run, and the approval covers builds that pass the gate.
- M2 (HPC, done; 6 of 8 cells usable): on the Gold 6348 nodes the timers cost 1 to 12 s of kernel time. The two
  IMPACT cells that ran on Xeon 6747P nodes (20,000 × 1 and 5,000 × 10, both paths) spent 517 to 2,247 s in the
  kernel on top of normal user time. pf73 and pf74's builds, with no timers, show no such cost on 6747P nodes,
  so it is the timers' clock reads there; those two cells' shares are not used. Shares of the judge's CPU
  account on the clean cells: the read call 44 to 69%; the site build 5 to 17% (the read-need test 3 to
  12%); BAM fetch and decode 4 to 39% (the largest at 1,000 rows). The rules path: fetch and decode 26 to
  82% (78 to 82% on IMPACT), the guard 1 to 18%, the row loop 1 to 23%. M3: the judged cache answers 0.1% of
  the judge's scoring clock; reads with best score 0 hold 85 to 87% of it in every clean cell.
- pf75 (HPC, done 2026-11-15; tools 11032837, builds 11032851 and 11032870, arrays 11032890 and 11032906):
  develop, S1 and S4, 6 lines a task, builds rotated by repeat; 144 lines, all exit 0, 3 good repeats in every
  cell. Option CPU over develop's (median of 3):
  - **S1 cuts CPU on both paths, in all 16 cells: kept.** Rules path 0.52 to 0.78 (IMPACT 0.52 to 0.60;
    ACCESS 0.53 to 0.78); judge path 0.76 to 0.96 (smallest gain on the 20,000-row cells). Wall time follows
    (0.53 to 0.96); peak memory 0.95 to 1.09. Worked example: IMPACT 5,000 rows x 10 BAMs, rules path,
    518 CPU s on develop, 270 on S1; the judge path 1,501 and 1,244.
  - **S4 does not cut CPU in every cell: not kept by the rule.** CPU 0.96 to 1.02 (below 1 in 6 of 8 rules
    cells and 5 of 8 judge cells; the most, 0.96, at IMPACT 1,000 x 1 rules, 0.8 s). Its planned 0.3 to 0.6 s
    per BAM did not show (the properties step is 0.03 s per BAM); its RNA part is not in these cells. It
    comes back to the operator for peak memory: ACCESS 5,000 x 9, judge, 1.14 (S4 3,713, 2,796 and 3,348 MB;
    develop 2,627, 2,956 and 2,938). One repeat is far from the other two, so it may be noise; nothing here
    shows which.
  - Against the C++ floor (pf73's ratios times pf75's, so about +/-10%): the rules path with S1 is 0.25 to
    0.74 of C++'s wall; the judge path with S1 is under 1 only at 1,000 x 1 (0.58 and 0.63), 1.14 at IMPACT
    5,000 x 1, and 2.2 to 5.0 at the deep ACCESS and 9- and 10-BAM cells. S1 alone does not bring the judge
    to the floor; the read call (M2: 44 to 69%) is next.
- 2026-11-15, a flaw in S1 (found before any PR): S1 set the Rust log cap when the module loaded, and the CLI
  loads the module before it sets the log level. So the cap stayed at WARNING until the first counting call,
  and the two INFO lines of the variant preparation were missing from every log (shown on the repo's test
  data; the other lines matched). The gate compared output rows only, so it passed. Fixed as `s1c`: the
  cap is set at the start of each counting call and opened again at its end; outside counting the logging is
  develop's. The counting loop is the one pf75 timed. The gate now has log arms (section 4); `s1c` goes
  through the whole gate, and pf76 times it again.
- 2026-11-15, the operator, on S4: "No we should test on RNA as well so keep s4". S4 stays; pf76 times it
  (and S1) on RNA. The memory flag at ACCESS 5,000 x 9 is read again in pf76's fill-out.
- 2026-11-15, the operator: the plan's one rerun waits until pf75 is done. If pf75 does not need it, it goes to
  M2's two IMPACT cells, pinned to Gold 6348 nodes.
- pf76 (HPC, done 2026-11-15; rows 11075567, arrays 11075595, 11075618 and 11075401): develop, the fixed S1
  (`s1c`) and S4. The FORTE BAMs were read in place (operator: "Can we not just read from original sources"):
  neutral-name soft links to 10 FORTE STAR BAMs near the median size (30.9 to 40.5 GiB; median of 58: 33.5),
  one RNA task at a time, no warm read; the RNA cells on the rules path only (`gbcms rna` refuses the judge).
  54 lines, all exit 0, 3 good repeats in every cell. CPU over develop's (median of 3):

  | Cell | s1c | S4 |
  |:--|--:|--:|
  | RNA 1,000 x 1 | 0.71 | 0.90 |
  | RNA 5,000 x 1 | 0.72 | 0.91 |
  | RNA 20,000 x 1 | 0.67 | 0.91 |
  | RNA 5,000 x 10 | 0.73 | 0.87 |
  | DNA ACCESS 5,000 x 9, rules | 0.69 | 1.00 |
  | DNA ACCESS 5,000 x 9, judge | 0.93 | 0.99 |

  - **s1c** matches pf75's S1 on DNA (0.68 and 0.93): the log fix did not change the speed. Peak memory 0.96
    to 1.07.
  - **S4 cuts CPU on RNA in every cell: kept.** Worked example: RNA, 5,000 rows x 10 BAMs, 368 CPU s on
    develop, 318 on S4 (the GTF and editing sites read once per run, not once per BAM). On DNA S4 is within
    noise (0.96 to 1.02 in pf75). Its pf75 memory flag was noise: the same cell is 0.96 here (1.14 in pf75).
  - Failures on the way: the first rows job could not open the link from inside its container (the
    container saw only the work folder), so the first RNA arrays were cancelled before any task ran; the
    rows job now binds the link target's folder read-only. My staging brief asked for the judge on RNA,
    which `gbcms rna` refuses; the RNA cells were cut to the rules path before any run.
- 2026-11-15, the gate on `s1c` (local, against develop 1b11ffea): **0 differences on every arm**, the three
  log arms included (HG002 and the 9-BAM fill-out on both paths, 8 logs each; the FORTE RNA runs with the GTF,
  33 logs, 2,287 lines). The app's sign-in stop ended the first chain; the restart ran detached and skipped
  the finished runs. One judge run (a truth set) finished just before the stop but has no time line in the
  arm's manifest (the line is written in job order, and the stop came first); its output has its DONE mark
  and was compared.
- 2026-11-15, Greptile on the S1 PR (#439, 4 of 5): two valid findings, both fixed (build `s1d`). (1) Counting
  releases the GIL, so another thread can lower a level and call `reset_log_caching()` during a count; the cap
  held the old level until the call ended. The reset now opens the cap. A new test lowers the engine logger
  to TRACE from inside the call and expects TRACE lines in the same call; it failed before the fix. (2) A test
  set the root logger outside `caplog` and left it at WARNING for later tests; it now sits in its own
  `caplog.at_level`. The log arms on `s1d`: 0 differences (the reset runs once, before any count, on the CLI
  path). My slip: I reported the review as not posted without reading it (LRN-20261115-002).
- 2026-11-15, S4's log arms: 0 differences too. They do not yet test S4's change: each RNA log run reads one
  BAM, so the GTF is read once in both builds. Before the S4 PR the gate gets an RNA fill-out log run (several
  BAMs, the GTF on), where S4 should show one GTF load, not one per BAM.
- 2026-11-15, S1 merged (#439, develop 1992fcb4); Greptile's second review 5 of 5, no new findings.
- 2026-11-15, the RNA fill-out log arm (`logs_rnafill`, now in the gate): every FORTE run of the tier in one
  `gbcms rna` run, 3 of them again under another tag (36 BAMs), on the union of their rows (57), with the GTF
  and a made-up editing database (39 sites at the union's SNVs). Develop and a repeat: equal, line for line.
  `s1d`: 0 differences. `s4`: 0 of 2,052 rows differ; the log differs only where BAMs 2 to 36 reuse the parse
  (6 line kinds only in develop, the parse's timing and DEBUG lines and the loader's own line; 2 only in S4,
  the reuse said at DEBUG; 35 times each). One local run, not a speed measure: 56 CPU s on develop, 4 on S4.
- 2026-11-15, the operator, on S4's log (option B: on each BAM after the first, the reuse said at INFO and the
  parse's warnings said again; a one-BAM log unchanged): "I am ok with B if there is no issue in accuracy".
  The check, on the 35 BAMs that read the index from memory (1,995 rows): every field equals develop's, and
  the GTF columns hold values (the exon distance in 1,995 rows, 44 values; the transcript counts in 1,855;
  the known-junction flags True in 19; `asjd_flag` is False in every row, so these rows do not test it). The
  same BAM read again from memory equals its fresh parse (3 BAMs, 57 rows each). Built as `s4b` (develop
  1992fcb4 plus S4 with option B).
- 2026-11-15, the gate on `s4b` (local, against develop 1b11ffea; 6,585 s): **every row identical on all 18
  arms** (0 rows differ, no run failed), and 0 log differences on `logs`, `logs_rna` and `logs_judge`. On
  `logs_rnafill` the log differs only as option B says: on each of the 35 BAMs after the first, develop's 10
  lines of the parse and the load (the two "Loaded N RNA editing sites" lines, "Loading GTF annotation from",
  the parser's summary, its timing, its DEBUG lines and "built index from") are replaced by the 2 reuse lines
  (7,180 lines compared; 350 only in develop, 70 only in `s4b`). The gate script reports that arm as not zero,
  as it must for a log that changes by design.
- 2026-11-15, the S4 PR (#440): CI's Rocky 8 job could not pull its image a second time (quay's file
  server; ERR-20261115-002); the rerun passed, 24 of 24 checks. Greptile (4 of 5): two valid findings, both
  fixed (build `s4c`). (1) The cache key held the path as given: a Python caller that moves to another
  working folder and passes the same relative name, a file of the same size and modification time, got the
  first folder's index. The key is now the resolved path and, on Unix, the device and inode (a file renamed
  over the path is read again too). A test reproduces the case; it failed on `s4b`'s key. (2) The quality
  tally kept every sampled quality (about 800 MB for 20,000 records of 20 kb, with its copy); it now tallies
  1 MiB at a time. The arms these reach, on `s4c` against develop 1b11ffea: every row identical (`rules`,
  905 runs, and the four log arms); logs as on `s4b` (0 differences on `logs`, `logs_rna`, `logs_judge`; the
  fill-out's planned reuse lines only). Builds are now 8 of 8: S2, S3, S5 and a combined build need more,
  the operator's call.
- 2026-11-15, a slip: reading the fill-out log, I cut lines to 200 characters before the sanitizer. A cut ended
  inside a FORTE slice's path, before `.bam`, and one folder name of the path was printed (no BAM name, read
  name or sample ID). The sanitizer now masks every path under a folder of patient slices, whole or cut
  short, and the PHI rule says to cut only after it (LRN-20261115-003).
- 2026-11-16, #440 merged (8a3966d8; CI on d442c720 all passed). Greptile's third review (5 of 5) came 89
  seconds before the merge, unread by me then; its one finding is valid and my fix caused it. The key held
  only the resolved file, but the loaders read the name: the editing database picks gzip from a `.gz`
  ending, and the GTF's warnings print the name. In one process, a second name for a database loaded
  under the first reused its sites, so a name that fails in a fresh process passed (plain bytes named
  `.gz`: "invalid gzip header"; gzip bytes without `.gz`: "stream did not contain valid UTF-8"), and a GTF
  warning said again printed the first name. A CLI run uses one name, so its reuse was right. The fix
  (option B, the operator's "Ok"): the key holds the name as given and the resolved file, so every reuse
  gives what a fresh call under that name gives, its log included. The operator approved the change of two
  lines of an existing Rust test (one file under two spellings, and a link and its target, each gave one key;
  now two). Three new tests (both database cases, and the GTF warning under a second name) failed on the
  merged build and pass now. No gate build: within one CLI run every call passes the same name, so the
  reuse is the one `s4c` measured (the two-BAM test still reuses each parse once). Builds stay 8 of 8.
  Learning: LRN-20261116-001.
- Slips: the agents' worktrees started on `main` (0f5c1f0b), not develop; each moved to 1b11ffea. On
  2026-11-15 I listed M2's copy-back lanes with `find` (they hold only cell folders and timing files, and no
  identifier was printed, but the rule counts every lane as a patient folder); later scripts build lane paths
  from the known cell layout.
- 2026-11-16, S4's fix merged (#441, develop 473e8920). Phase 2 (S2, S3, S5, the judge's kernel) has its own
  page: `335-phase2.md`.

## 8. Cost log

- Research agents: inventory (Sonnet) 0.20 M tokens; code survey (Opus) 0.44 M; outside survey
  (Fable) 0.18 M; profiles (Sonnet) 0.16 M.
- Build agents (Sonnet): M1 0.34 M; S1 and S1b 0.25 M; S4 0.26 M; the gate 0.23 M; pf76 staging 0.19 M. Review rounds 0;
  builds 8 (with `s1c`, `s1d`, `s4b` and `s4c`); HPC batches 3 (M2, pf75, pf76); reruns 0.
