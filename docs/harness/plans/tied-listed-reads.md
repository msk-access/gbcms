# Naming the reads the judge holds back

Status: building 2026-11-12 · option B, widened, parts 1 and 2 with the own-deletion condition
(operator, 2026-11-11: "Yes, widen it"; 2026-11-12: "Yes, with it") · ticket: none
yet (a sub-issue of E #262, beside #273) · serves RJ-34 (`docs/reference/read-judgment.md`, decided
2026-11-11) · count-neutral by construction (the plan-first page is used for its form; this item
is not one of the trial's count-changing items, and it keeps a cost log)
Rounds: 0 of 2 · New conditions: 0 of 3 · Builds: 2 measurement builds used (trace-only, never merged;
operator, 2026-11-12: "Ok lets do it", then "Yes, one more trace"), the final build not yet · HPC
batches: 0 of 1 · Part 3 and the open question moved to #306 (operator, 2026-11-12)

## 1. The class

- **Reads:** reads the sequence judge holds back from AD (depth only, or partial_alt) whose bases
  do carry an allele, while no field of the row says which. Three kinds:
  - **Tied with a listed row.** The row's ALT and another listed row's allele read the same over the
    read, so the read fits both exactly. RJ-25 credits a tie to neither. The listed row is either a
    candidate state at the site (a listed allele the row's site meets) or a change the best REF path
    carries beside the site.
  - **A listed allele at the site that the reads fit better.** The reads carry another listed row's
    allele at the same site (for example, a listed delins over a listed deletion). They tell it from
    the row's allele, so they are that row's reads, partial here (RJ-25: a listed allele at the site
    is always another allele). The allele-naming scan leaves listed alleles out ("already an input
    row": `observed.rs:76` and `:196`, `known_alleles` at `:206-212`), so the row says only
    `PARTIAL_DOMINANT`.
  - **An unlisted allele under the naming floor.** The read carries an allele nobody listed (for
    example, a deletion that ends at an SNV's site, so the molecule is a delins). The judge does not
    keep it, because it is under the presence test (3 fragments and 5%, RJ-25), and the scan names
    an allele only on 5% of the scanned reads (`observed.rs:32-34`, `named` at `:200-204`).
- **Evidence the diagnostic must report:** how many reads, and which allele they carry. It never
  reports which row "should" have them.
- **Invariants:** 7 (no count changes; the flag names what the reads carry, which is what
  `gbcms_diagnostic` is for); 1 (a per-read decision, so the flag is the same under any bin
  geometry); 2 (no change to the quality gate); 3 (computed only where the judge runs, from scores it
  has already computed); 5 (a new `#[pyo3(get)]` field changes `_rs.pyi`).
- **Principle:** RJ-34. Where the judge holds back reads that carry a change at the row, the row
  names what they carry and how many.
- **DNA and RNA:** DNA under the judge. On the rules path, only the SNV rows it sends to the judge
  (#398). RNA: the judge refuses RNA mode today (#333 decides it), so RNA output is unchanged.
- **Worked example (made-up bases).** The reference reads `GAT TGCAG TTA`, with the site at the C. A
  listed SNV makes it `GAT TGAAG TTA`. A listed insertion after `GAT` adds `TGAAG`, a copy of the
  site's region that carries the A, so its molecule reads `GAT TGAAG TGCAG TTA`. The read
  `GATTGAA` is in both molecules, so it fits both exactly and the judge counts it for neither (depth
  only). A read that goes on to `TTA` is the SNV's. One that goes on to `TGC` is the insertion's.
  Today the SNV row says only `ZERO_ALT`.
- **Instance:** the 6.6.0 tier, 2026-11-11 (develop 44b02197; the batch's judge arm).
  - Tied with a listed row (from the whole-tier measure: the judge cuts AD to half or less, or
    to 0, against the 6.7.0 default):
    - an SNV beside a listed 53-base insertion whose copy holds its ALT base: 42 to 0;
    - an MNV beside a listed 24-base insertion whose copy holds its two bases: 13 to 0;
    - a 3-base deletion with a listed MNV 6 bases on: 82 to 18 (and 42 to 15 in a second run);
    - an MNV with a listed 3-base insertion right before it: 45 to 7.
  - A listed allele at the site that the reads fit better (signed-out rows the judge counts below
    the sign-out; in the 6.7.0 default too, so the counts are decided, and only the name is
    missing):
    - a 1-base deletion under a listed 3>2 delins: sign-out 173, the judge 0, partial 177; 175 reads
      the rules path calls "other" are partial under the judge;
    - a 3-base deletion under a listed 5>2 delins: sign-out 21, the judge 2, partial 24; 24 reads.
    Each deletion and its delins are different sequences (not one allele written two ways).
  - An unlisted allele under the floor:
    - an SNV with a 69-base deletion ending at the site, on about 1% of the reads (8 reads): 9 to 1
      (and 7 to 1 in a second run);
    - a 6-base deletion row with no listed row within 150 bases: 11 of 247 reads (4%) partial, and
      the judge weighed no found allele there. What those 10 reads carry is not yet read.

  An edit-distance check that does not use the judge's scorer agrees at the tie rows. Of the reads
  held back, 42 of 42, 60 of 64, 35 of the 35 it read, and 11 of 13 fit the row's allele and the
  listed allele equally, most of them exactly (distance 0). The others fit the listed allele better,
  except one read at the first MNV row that fits only the MNV (the judge counts 0 there; the check
  counts 1).

## 2. What already touches it

- **Rules:** RJ-34 (name every held-back read; the bar for #331) and RJ-33 (the CIGAR never breaks
  a tie); RJ-25 (read to its end: a tie is never credited; a listed allele at the site is always
  another allele; the presence test for found alleles); RJ-23 (a read that ties two alleles must
  hold both alleles' windows, or it is depth only); RJ-28 (`MNV_SHARED` and `INDEL_SHARED`, the
  model: one molecule shown at the rows that share it, no count changed).
- **Code:**
  - the call: `rust/src/counting/judge/decide.rs:303-319` (`decide`, the classes' best scores and
    `DELTA` = 4, `score.rs:26`) and `:322` (`tied_with`: another allele tied with ALT or REF,
    the third class beaten, is partial; every other tie is depth);
  - which state won for each class: `decide.rs:509-517` (`score_once`, `best: [(score, state)]`);
    the listed changes a path takes: `site.taken_changes` (used at `decide.rs:496-501`);
  - the model to copy: `rust/src/counting/engine.rs:2086-2117` (`mnv_shared`, `indel_shared`) and
    `rust/src/types.rs:546-561`;
  - the tokens: `src/gbcms/pipeline.py:1186-1212` (`_diagnostic_flags`, documented at
    `:997-1087`); the VCF writer puts the field in INFO (`io/output.py:1239`); the "check" list:
    `src/gbcms/qc_summary.py:126-151`;
  - the naming scan: `rust/src/counting/observed.rs` (`MIN_CARRIERS` 3 and `MIN_FRACTION` 0.05 at
    `:32-34`; `named` at `:200-204`: n >= 3, n more than the given allele's carriers, and 5% of
    the scanned reads; listed alleles left out by `known_alleles`, `:206-212`); the judge's own
    presence test: `carrier::is_found` (`carrier.rs:292`, `GUARD_MIN_READS` 3, `GUARD_REF_SHARE`
    0.05).
- **The tokens' meanings** (`docs/reference/qc-flags.md`): `OBSERVED_ALLELE(…:n/0)`, no read
  carries the given allele and n reads carry this one; `COEXISTING_ALLELE(…:n/m)`, the given allele
  is present on m reads and more reads (n) carry another allele, "a VAF caveat". The widened option
  keeps both meanings.
- **Tickets:** #273 (explainability share: it counts OTHER apart from uninformative reads, so tied
  reads, which are uninformative by design, are outside what it measures); #262 (diagnostics
  parent); #306 (the 5% share is a cliff); #410 (a listed MNV and a listed insertion inside it); #417
  (a listed change beside a deletion moves its carriers); #331 J and J6 #337 (the migration notes
  must name the new flag and the naming changes).
- **REJECTED.md:** two entries near this class; neither vetoes naming.
  - REJ-20261102-001 vetoed breaking a tie by elimination (the own-gap rule for every state). This
    plan breaks no tie.
  - REJ-20261031-001 vetoed deciding a listed allele absent (with a `UNIT_COUNT_ABSENT` name): a
    listed allele is always a candidate. This plan decides nothing absent and changes no count.

## 3. Options

| | Rule in one sentence | Invariant 7 | Conditions | Will not handle |
|:--|:--|:--|--:|:--|
| A | Name each listed row whose allele ties the ALT on held-back reads: `TIED_WITH(chrom:pos:REF>ALT:n/m)`, n such reads, m the row's AD | pass (no count) | 1 | a listed allele the reads fit better; the unlisted allele |
| B, as drafted | A, and name a found allele the judge did not keep where at least 3 reads carry it and they outnumber the row's AD, with `OBSERVED_ALLELE` / `COEXISTING_ALLELE` (no 5% floor for naming) | pass (no count) | 2 | a listed allele the reads fit better (the scan leaves listed alleles out) |
| **B, widened** | B, and let `OBSERVED_ALLELE` / `COEXISTING_ALLELE` also name a listed allele at the row's site that held-back reads carry better than the row's allele | pass (no count) | 3 | an allele fewer than 3 reads carry; an unlisted allele on fewer reads than the row's AD (measured, below) |
| C | #273's share: `LOW_EXPLAINED(x)` where (RD + AD) / informative reads is under a threshold | pass (no count) | 1 and a threshold | tied reads (they are uninformative, so the share does not see them); it does not name the allele |

- **The widened option's three parts:**
  1. **Ties:** `TIED_WITH` (A).
  2. **A listed allele the reads fit better:** `OBSERVED_ALLELE` when no read carries the row's
     allele, else `COEXISTING_ALLELE` when more reads carry the listed one: the scan stops leaving
     out a listed allele at the row's site when the judge holds its reads back at this row.
  3. **An unlisted allele:** named from 3 reads, the presence test's own mechanism (three
     molecules, not an error), without the 5% share. The rule keeps "more reads than the row's
     allele", because `COEXISTING_ALLELE` means that. **Moved to #306 (operator, 2026-11-12, on
     the emulation in section 7):** the naming scan reads only the visible end of a deletion that
     runs past its window, so it names an allele the reads do not carry, and without the 5% share
     it flags noise-sized counts on matched normals. This item keeps parts 1 and 2.
- **Measured before it is decided (an open question for the operator; moved to #306 with part
  3):** rows where 3 or more
  held-back reads carry an unlisted allele but fewer than carry the row's allele. The tokens'
  meanings do not cover them, and RJ-34 asks that they be named. The emulation counts them per data
  type, and how many are signed-out rows the judge counts below the sign-out. Then the operator
  chooses: leave them (few, and the row's AD is the larger), or a form for "fewer than the row's
  allele".
- **Both paths:** the flags run only where the judge runs. Under the 6.7.0 default that is the
  routed SNV rows; with the judge as the default, every DNA row.
- **Binning invariance:** the flags come from the read's own states and the reads at the row, so
  they are the same under any bin window or cap. The binning-invariance check covers
  `gbcms_diagnostic` too.
- **Quality contract:** unchanged. The flags read the scores the judge already computed.
- **Decided (operator, 2026-11-12, "Yes, with it"), the rule the final build implements:**
  - **Part 1, `TIED_WITH(chrom:pos:REF>ALT:n/m)`:** a held-back read (call depth only or another
    allele; not "fits nothing") whose best ALT score is within `DELTA` of the best of REF and
    OTHER names listed row L when, for every class tying ALT (within `DELTA` of it), L is that
    class's own listed state, or L is a listed neighbour on that class's winning path such that,
    with L left out of every state's graph, ALT alone still ties the read's best score and that
    class falls `DELTA` or more behind ALT. A read whose own alignment writes a deletion of 20
    bases or more names nothing, unless that deletion is L itself (operator, 2026-11-12: "Yes,
    exempt them"; then it names L) or the row's own allele (2026-11-12, the same reasoning: a
    listed row explains the deletion; to be confirmed by the operator). n: the reads that name L, in the row's depth; m: the row's
    AD. A tie whose partner is a found allele is not named.
  - **Part 2:** a held-back read whose OTHER winner is a listed allele at the row's site, beating
    ALT and REF by `DELTA` or more, counts toward that allele; named `OBSERVED_ALLELE` where n >= 3
    and AD is 0, `COEXISTING_ALLELE` where n >= 3 and n > AD.
  - Part 3 and the open question: #306.
  - n1's count is accepted as measured (8 of its 68 reads fit the named MNV better than ALT; the
    two scorers' difference there is RJ-31 work, not this flag's).
- **Recommendation (2026-11-11): B, widened.** It names all three kinds the tier shows, adds one token
  (`TIED_WITH`) and no column, keeps the existing tokens' meanings, and follows the form of
  `MNV_SHARED` and `INDEL_SHARED`. C stays #273's question. It is a different measure, and it cannot
  see this class.

## 4. Measurement (fixed before the runs)

- **Emulate first, with no build:** the batch's judge arm (`batch680/judge`, the whole stderr of
  every run, the judge's own lines included), and a rerun with `--keep-judge` where a line is
  missing (local).
  - Part 1's rows: the held-back reads whose best ALT is within `DELTA` of the best non-ALT class,
    at rows with a listed row in reach.
  - Part 2's rows: held-back reads whose best state is a listed allele at the row's site, by
    `DELTA` or more.
  - Part 3's rows: from the judge's found-allele lines ("written by N of M reads, told on F
    fragments") where N >= 3 and N is more than the row's AD; and, where the judge weighed no found
    allele (the 6-base deletion row), what the held-back reads carry by the independent check.
  - The open question's rows: as part 3, with N >= 3 and N at most the row's AD.
- **Class metric and bar:**
  - (1) every field except `gbcms_diagnostic` is identical on every row of every arm (count-neutral
    proof);
  - (2) the rows in section 1 carry the flag: the tie rows with n within 10% of the independent
    check's tie count (`read_oracle.py`, harness); the two listed-delins rows name the delins with
    n within 10% of the reads the check puts on it; the two unlisted-allele rows name the allele
    (the 6-base deletion row once its reads are read);
  - (3) a sample of 10 newly flagged rows per data type, read by the independent check: at least 9
    of 10 flags name an allele the reads carry;
  - (4) the number of newly flagged rows per data type is reported, and the open question's count.
    On matched normals, a flag is expected only where a listed tumour row ties a neighbour in the
    normal's reads too;
  - (5) no flag from a scorer fault: at the place where a listed insertion and a listed deletion
    lie 3 bases apart, the judge scores REF and ALT equally on 7 reads per row while the
    independent check gives the row's change the same gain as on the reads it counts. `TIED_WITH`
    must not fire there; the fault is fixed first, or the flag waits for it.
- **Data types:** IMPACT tumour and matched normals, ACCESS duplex and simplex, WES and HG002 (BWA
  only; HG002 positions are public and may be shown), RNA (the judge refuses it: byte-identical).
- **Effects map:** changes only `gbcms_diagnostic` (MAF and VCF writers), the QC summary's check
  list, and anything that copies the field (`gbcms merge`, the observations export: check that
  they pass the token through). `OBSERVED_ALLELE` and `COEXISTING_ALLELE` appear on more rows (parts
  2 and 3): anything that filters on them sees more rows, with the same meaning. Unchanged by
  construction: every count, `vaf`, the fragment counts, mFSD read sets, Nextflow (no parameter).
- **Cannot see:** a tie with an allele that is neither listed nor found (two reads or fewer). Those
  reads stay unnamed depth, as `UNEXPLAINED_READS` reads stay today.
- **Contract tests:**
  - red-first: a made-up SNV beside an insertion whose copy holds its base gives `TIED_WITH` naming
    the insertion, with counts unchanged;
  - red-first: a made-up deletion under a listed delins at the same site, whose reads carry the
    delins, gives `OBSERVED_ALLELE` naming the delins, with counts unchanged;
  - red-first: a made-up SNV with a deletion ending at it, on 2% of the reads, gives
    `OBSERVED_ALLELE` naming the delins;
  - keep `tests/test_diagnostic_flags.py`, `test_indel_shared_contract.py`,
    `test_observed_allele_contract.py` and `test_haplotype_judge_contract.py`.
- **HPC:** none for the decision (the tier is local); at most one gate batch for the final build.

## 4b. Related 6.8.0 tickets

- **#273:** this plan takes the tie class out of it. #273's share still needs its own decision.
- **#410, #417:** two of the tie rows are near these shapes: an MNV with a listed insertion right
  after it, and a deletion with a listed MNV past its end. The flag shows these rows' reads while
  those tickets decide their counts.
- **#306:** part 3 and the open question moved there (operator, 2026-11-12). They name alleles under
  the 5% cliff, so they belong with the cliff's own decision. #306 also needs the naming scan to
  widen over a deletion that runs past its window (the 69-base deletion row, section 7).
- **#337 (J6):** the migration notes name `TIED_WITH` and the two naming changes.
- **A separate question, not this item (it changes `partial_alt`):** the same evidence is filed in
  two places. A read that ties ALT and another allele, with plain REF beaten, is partial when the
  other allele is a found allele (the judge enters it as another state) and depth only when it is a
  listed neighbour (the judge puts it on the REF path). The emulation counts the reads it would
  move; the operator decides it on its own.
- **Also seen, not this item:** at two of the tie rows, the judge counts the neighbour's carriers
  REF (46 and 33 reads, which the 6.7.0 default puts in partial_alt). By their best alignment those
  reads hold the reference at the row's bases. Record it with #417, not here.

## 5. Budget and stop rule

- At most two review rounds, and the widened option's three conditions.
- Past either limit, stop and take one class question back to the operator.
- Builds: one final build; the row sets are emulated first. Measurement builds (trace-only, never
  merged): two used, each on the operator's word.

## 6. Work and models

- The emulation: scripts (no model tokens).
- The build and the red-first tests: Sonnet.
- The review of the diff and of the flags' wording: Fable.

## 7. State

- 2026-11-11, measured: the whole-tier measure (13 rows new in 6.8.0; 5 with no record of the reads
  the judge holds back) and the read-by-read look at those 5; the signed-out rows the judge counts
  below the sign-out, by cause (37 rows; RJ-34's table), which showed the listed-allele kind and
  the 6-base deletion row (tier harness: `absent_rows.py`, `missed_signout.py`, `row_picture.py`,
  `read_oracle.py`).
- 2026-11-11, operator: option B, widened. While writing it, the existing rule's "more reads than
  the given allele" turned out to be part of `COEXISTING_ALLELE`'s meaning (`named`), so dropping it
  is now the measured open question above, not a rule.
- 2026-11-12, measured: the emulation. The batch's logs give each read's best score per class, not
  which state won it, so a trace-only build (develop 44b02197 plus two trace lines: the state that
  won each class and the rows its path takes; every allele the naming scan counted before its
  floor) reran the judge arm on all four local sets (868 runs, 3,520 CPU-seconds).
  - Control: all 904 output files (MAF, VCF, Parquet) identical to the batch's judge arm; all 868
    logs identical line for line once the two new lines and the time and version lines are left
    out.
  - Bar 2, ties (the check run on the same reads, every held-back read at the row):
    - the MNV with a listed 3-base insertion before it: 42 against 38 to 39 (pass);
    - the 3-base deletion with a listed MNV 6 bases on: 68 against 60 (fail, 13%);
    - the MNV beside a listed 24-base insertion: 14 against 11 to 12 (fail);
    - the SNV beside a listed 53-base insertion: 0 against 42 (fail).
    At the last one, on 75 reads every class's best path takes the insertion and scores the same
    (the read ends in the insertion's copy). The trace cannot tell these reads from reads that
    carry a listed neighbour and never reach the site: that needs the scores of ALT and REF
    without the neighbour, which the trace does not log. At the deletion, the judge ties 8 reads
    the check tells apart (cause not read; a likely cause, not verified: the lowered qualities
    beside a read's own indel, RJ-31).
  - Bar 2, a listed allele the reads fit better: the 1-base deletion names the listed 3>2 delins on
    174 reads (the check 175); the 3-base deletion names the listed 5>2 delins on 22 (24). Pass.
  - Bar 2, an unlisted allele: fail at the SNV with a 69-base deletion ending at it (both runs). The
    aligner writes the 69-base deletion on 8 reads, but the deletion runs past the scan's reference
    window, so the scan does not widen over it and reads its visible end: it would name a 3>1
    delins (8 reads; 5>1 on 10 reads in the second run), an allele the reads do not carry. The
    judge's own found-allele search sees long events, but it logs only those over its floor (the
    deletion is on about 1% of the reads). At the 6-base deletion row the scan names a 7>1 delins
    11 bases on (8 reads), not yet read.
  - Bar 5: pass. No `TIED_WITH` at the two rows of the scorer-fault place.
  - Bar 4, newly flagged rows (main data types; the reruns with other arguments left out):
    - IMPACT tumour: 140 of 12,977 rows. `TIED_WITH` 78 rows (46 with 3 or more reads), 976 reads
      (423 depth only today, 553 partial); part 2, 59 rows; part 3, 20 rows.
    - IMPACT normals: 21 of 631 rows: part 3, 19 rows, 18 of them at under 2% of the scanned reads
      with no read carrying the given allele (3 or more reads of several hundred). Bar 4 expects a
      flag on normals only where a listed tumour row ties a neighbour, so without the 5% share
      part 3 puts `OBSERVED_ALLELE` on matched normals from noise-sized counts.
    - ACCESS: duplex 14 and simplex 16 of 326 rows each; normals 2 of 10 each (part 2).
    - WES: 1 of 342 rows (part 3). HG002: 705 of 4,982 rows, nearly all `TIED_WITH` in the set of
      adjacent insertion and SNV pairs built for this shape.
  - The open question (an unlisted allele on 3 or more reads, fewer than the given allele's):
    IMPACT tumour 605 rows (385 with 3 to 5 reads; 309 under 1% of the scanned reads; the
    commonest is another base at an SNV's site, 200 rows); ACCESS 23; WES 3; normals 2; HG002 210.
  - The separate partial-or-depth question: of the 976 tied tumour reads, 423 are depth only
    today (a listed neighbour on the REF path) and 553 are partial (a listed state at the site).
  - Bar 3 (a sample of new flags read by reads) is not run: it waits for the rule.
  - Tools (harness): `name_heldback.py` (the emulation; local rows table) and `name_bars.py` (bars
    2 and 5); the build is `src_trace680` (its patch beside it).
- 2026-11-12, operator: part 3 and the open question to #306; one more trace for part 1 ("Yes, one
  more trace").
- 2026-11-12, measured: the second trace. A second trace-only build (the first plus one line: on
  each held-back read, each class's best score with each listed or found neighbour its winning
  paths take left out of every graph) reran the judge arm (868 runs, 4,049 CPU-seconds, none
  failed).
  - Control: every output file identical to the first trace's; every log identical once the new
    line and the judge's time per read are left out; the scores laid again with nothing left out
    equal the judge's own on every read scored once.
  - The part 1 rule measured: a held-back read names neighbour L when, for every class tying ALT,
    L is that class's own listed state, or L is on its winning path and, with L left out, ALT
    alone still ties the read's best score and that class falls 4 or more behind ALT.
  - Bar 2, read by read. The check is the same as before (edit distance), with one change: its
    haplotype for the neighbour now holds L alone, as the flag names it (the first count applied
    every listed row within 150 bp together):
    - the SNV beside the listed 53-base insertion: 42 named, 42 ties by the check (pass; 0 before);
    - the MNV with the listed 3-base insertion before it: 44 named, 44 ties (pass);
    - the MNV beside the listed 24-base insertion: 14 named, 14 ties, 13 the same reads (pass);
    - the 3-base deletion with the listed MNV past its end: 68 named, 60 ties, all 60 named; the
      other 8 named reads fit the listed MNV better than ALT. The flag names the right allele,
      but its count is 13% over the check's ties (fail by the letter of the bar). On those 8 reads
      ALT's winning path takes no neighbour and scores a perfect 0, where the check counts 2 edits
      for ALT: the two scorers differ on these reads (cause not read; the lowered qualities beside
      an indel, RJ-31, are one candidate).
    - Part 2 as before (174 of 175; 22 of 24). At the first of those rows, 5 reads also tie ALT
      with the same listed delins (the check: 4 ties, 1 fits the delins better).
    - Bar 5: pass.
  - Bar 3 (a seeded sample of 10 flagged rows per data type; every flag at those rows; a flag
    holds where at least half its reads fit the named allele by the check):
    - IMPACT tumour 10 of 10; HG002 10 of 11; IMPACT normals 3 of 3; the reruns with other
      arguments (cohort 13 of 13, Smith-Waterman 11 of 11). Every part 2 flag held.
    - ACCESS: duplex 13 of 16, simplex 14 of 24 (fail). Every failure is a `TIED_WITH` flag of 1 to
      6 reads at two loci of one sample, a dense cluster of listed indels. The failing reads write
      deletions of their own (50 to 85 bases) that no row holds; by the check, REF fits them as
      well as anything.
    - A condition measured for that: a tied read whose own alignment writes a deletion of 20
      bases or more names nothing. With it: ACCESS duplex 13 of 14, simplex 14 of 15 (pass); IMPACT
      tumour loses 29 tied reads (1 flag); HG002 4 reads; at the bar 2 rows, 1 of the 44 ties at
      the 3-base insertion row.
    - Every flag, not a sample: IMPACT tumour 150 of 152; HG002 665 of 740 (90%). HG002's failures
      are adjacent insertion and SNV pairs: on 198 of their 318 reads the check finds ALT fits
      better than the neighbour, where the judge scores them within 4.
  - Bar 4, newly flagged rows (parts 1 and 2): IMPACT tumour 123 of 12,978 (`TIED_WITH` 81 rows, 47
    with 3 or more reads, 1,017 reads: 465 depth only today, 552 partial; part 2, 59 rows); IMPACT
    normals 2 of 562 and 0 of 70; ACCESS duplex 11 and simplex 8 of 316; ACCESS normals 2 of 10
    each; WES 0 of 342; HG002 737 of 4,986.
  - A limit of this trace, not of the rule: it keeps the first 4 neighbours per read. At 18 dense
    loci of the IMPACT tumour runs (101 rows), a tie's partner lay past them, so the tie counts
    there may be low. The final build tests every neighbour.
  - Tools (harness): `name_heldback.py` (the second rule; the first is `name_heldback_v1.py`),
    `bar3_check.py`, `name_bars.py`, `batch680/control_trace2.py`; the build is `src_trace680b`.
- 2026-11-12, operator: "Yes, with it": the final build, parts 1 and 2, the own-deletion condition
  in part 1, n1's count accepted as measured (section 3, "Decided").
- 2026-11-12, built and measured: the final build (branch `feature/tied-with-flag`; the red-first
  contract tests by a Sonnet agent, 8 tests in `tests/test_tied_with_contract.py`; the diff
  reviewed by a Fable agent). Measured twice on the judge arm (868 runs each): before the review's
  fixes and after them (row indices in place of a 64-row limit; at a long insertion, the class
  must tie ALT on the first graphs too; a neighbour placed at both ends of a run left out whole).
  The two builds' outputs are identical on every run.
  - Bar 1: pass. Every MAF and VCF field but `gbcms_diagnostic` equal to develop's judge arm on
    all 30,703 MAF and 1,462 VCF rows; the other 36 files identical; no token removed. 926 rows
    gain a token (`TIED_WITH` 835, `COEXISTING_ALLELE` 100, `OBSERVED_ALLELE` 23; see [QC Flags](../../reference/qc-flags.md)).
  - The build counts the held-back reads in the row's depth: a read off the anchor that the
    judge does not admit is in no count, so it names nothing (the emulation counted every judged
    read). Read by read, this gives the build's n at every bar row.
  - Bar 2: the 53-base insertion row 42 (the check 42); the 3-base insertion row 38 and the
    24-base insertion row 12 (the check's ties among the reads in depth: 38 and 12 to 13); the
    3-base deletion row 68 (accepted); the delins rows 171 (175) and 22 (24). Pass.
  - Bar 3 (the emulation's sample, with the own-deletion condition): IMPACT tumour 10 of 10, HG002
    10 of 11, ACCESS 27 of 29, normals 3 of 3. Pass.
  - Bar 4: IMPACT tumour 113 of 12,978 rows (`TIED_WITH` 72 rows, 841 reads; `OBSERVED_ALLELE`
    13; `COEXISTING_ALLELE` 42); IMPACT normals 2 of 562 and 0 of 70; ACCESS duplex 11 and simplex
    8 of 316, normals 2 of 10 each; WES 0 of 342; HG002 713 of 4,986.
  - Bar 5: pass (no flag at the two scorer-fault rows).
  - Time: the judge arm took 3,223 s against develop's 3,169 s (+1.7%, one run each, traces on).
  - The review found no count-changing bug. Open for the operator: in IMPACT tumour, 15 of the
    30 tied reads the own-deletion condition drops carry the named listed deletion itself (the
    condition silences reads that carry the allele it would name); in ACCESS every dropped read
    carries another deletion, the case the condition was decided for.
- 2026-11-12, operator: "Yes, exempt them": a read whose own long deletion is the named row
  names it (section 3). The red-first test for it showed the same reads at the other row of the
  pair, where their deletion is that row's own allele; the build exempts that case too (a listed
  row explains the deletion; to be confirmed by the operator).
- 2026-11-12, measured: the build with the exemption (868 runs). Bar 1 passes as before. Four
  output files change against the build before it: one IMPACT tumour run (and its
  Smith-Waterman rerun) gains two `TIED_WITH` flags carrying the 15 reads that carry the named
  deletion; one ACCESS simplex row (and its mFSD rerun) gains one read through the row's own
  deletion. IMPACT tumour: 113 rows flagged, `TIED_WITH` on 74 rows and 856 reads. This run's
  time is not a measurement (the machine's load average was over 50 while it ran); the clean
  comparison is the build before it against develop: median per-run time ratio 0.98, total
  +1.7%.
- Open: the operator's word on a test accessor widened to count both graph caches (section 8);
  the order of the two PRs (RJ-34's register text is on `docs/aligner-evidence`); the ticket and
  the PR (drafts to be shown before they are posted).

## 8. Cost log

Rounds 1 (the Fable review of the diff); conditions 1 added (the own-deletion condition,
operator, 2026-11-12); builds 2 trace-only measurement builds (never merged) and the final build,
measured twice (before and after the review's fixes); HPC batches 0; agents 2 (Sonnet: the
red-first contract tests; Fable: the review). One existing Rust test's accessor, `kept_graphs()`,
now counts the bytes of both of the site's graph caches (the test's expectations unchanged; for
the operator, section 7).
