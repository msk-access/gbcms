# Site candidates and REF: which alleles compete at the ALT's site, and when a read counts REF

Status: approved 2026-11-02 for its measurement step (one class by the operator's decision of 2026-11-02: C2 round 2, #394, #398, the delins REF geometry)
Rounds: 0 of 2 · New conditions: 0 of 1 (the delins core) · HPC batches: 0 of 4 (lanes on the builds in hand, one build per option, one final)

## 1. The class

- **Reads:** reads at a row X that fit REF or X and also fit a listed allele (S) over the bases they
  hold: S's span meets X's site, and the read reads no confident base where the two differ, or ends
  before that base.
- **Evidence:** which listed alleles are candidates at X's site; what a REF read must read against
  them; and whether a long event's core (the bases X deletes or replaces, past the junction windows)
  holds candidates or contexts.
- **Invariants:** 7 (every listed allele is present; AD only on the bases); 2 (one `--min-baseq`
  gate for the telling base); 1 (per-read tests, no bin dependency); both paths alike (today they
  differ in three places).
- **Principle:** one site test says which alleles compete, one reach rule says what a REF read must
  tell, one core rule for long events; both paths.
- **Instance:** C2 round 2 (the rules path's REF loss), #394, #398, the review's delins row.

## 2. What already touches it

- **Rules (read to the end):** RJ-5 (a read spelling another listed allele carries it; fitting
  two equally, neither for certain); RJ-9 (a read ending inside a long complex event is judged by
  the bases it has); RJ-24 (at a pure deletion of 50 bp or more the core is a context); RJ-25 (the
  1-fragment-and-5% sibling floor; the REF clause of 2026-10-14: REF only where told from every
  REF-length allele, else neither); RJ-28 (a read carrying X's change inside a longer listed allele
  is X's ALT).
- **Rules path:** `rust/src/counting/carrier.rs:2671-2688` (`spells`: a sibling is the allele's
  length and differs nowhere past `FLANK` 2 beside the event; `seen()`: the floor);
  `carrier.rs:1177-1182` (`tells_others`, `ref_untold`); `junction.rs:58` (`LONG_DELETION` 50, pure
  deletions only).
- **Judge:** `rust/src/counting/judge/haplotypes.rs:1811-1826` (`place`, `meets`);
  `haplotypes.rs:1834-1848` (`deleted_core`: pure deletions of `LONG_EVENT` 50 only, past
  `JUNCTION_WINDOW` 4); `haplotypes.rs:1747-1759` (`reaches_state`); `decide.rs:644-663`
  (`ref_by_reach`, and the own-gap tie-break: `GBCMS_JUDGE_OWNGAP=all` measured, kept narrow, no
  class reason).
- **Tickets:** #394, #398 (open); #304 (the reach rule lives in `docs/proposals/sequence-judge.md`,
  not the register); #331; #277, #300, #285, #297. **REJECTED:** REJ-20261031-001 (a listed allele is
  never decided absent); REJ-20261013-001, REJ-20261009-003 (read-on through a long deletion's
  bases). Nothing vetoes a site test or a reach rule.
- **Measured:** C2 rounds 1 and 2 (`~/test/gbcms/harness/plan670/plan-first/c2_results*.md`):
  rules `site` lane c2RS (build aa095d07), judge c2SA (cb62a63f). Aggregates only.

## 3. Options

Shared by A and B: **candidates** are every listed allele whose span meets X's site (the judge's
`meets`); an ALT-side tie is never credited (RJ-25). The site test's two edges are named, not
changed: a length change's left anchor base counts and the base after the span does not; a
substitution against an indel is tested on the indel's left-aligned span, not its run. The lane
counts the alleles each edge turns into a context; a count above 0 on the normals is a round-2
class question.

| | A. Tell every candidate | B. Tell what the read reaches |
|:--|:--|:--|
| REF rule | REF only where a confident base tells REF from every REF-length candidate; short of it, neither (RJ-25's clause on both paths; the judge's reach and own-gap tie-breaks go) | REF where the read's confident bases tell REF from X and from every REF-length candidate whose difference its aligned span reaches; a candidate past the read's end does not contest it (RJ-9); reached and not told, neither. The judge's rule today, written into RJ-25 and into the rules path |
| Long cores | the core of a pure deletion or delins of 50 bp or more, past the 4-base junction windows, is a context (RJ-24 extended to delins). Measured (rtb1): it gives back none of the ties, because all 10,352 lie in a junction window. To keep them REF, A needs a rule for each junction: B's reach rule for long events, and a third REF rule, which triggers the stop rule (section 5) | RJ-24 as is, extended to delins of 50 bp or more (a read that reaches and spells S inside the core is REF with a context); the ties need no core rule |
| Invariants 7, 2, 1; both paths | pass; alike | pass (AD equals A's); alike |
| Conditions | −3 (the floor on listed alleles, `spells`' scope, `ref_by_reach`), +1 (the delins core); own-gap stays narrow (decided 2026-10-27; rtb1 confirms it) | −3 (the floor, `spells`' scope, `ref_untold` for unreached candidates), +1 (the delins core); own-gap alike |
| Expected effect, rules | c2RS, plus #398's reads out of AD (69 at one row, 1,228 at 14 SNV rows: partial). Not as expected: the 3 delins-50+ rows' REF reads do not come back, because S is within 4 bases of the delins's end, where the core rule cannot reach it, so the three v2 failures stay | c2RS, plus #398's reads out of AD, plus the 1,313 reach reads back (the normals' delins RD, judge 4,303 against rules 3,558, closes); the three v2 failures go |
| Expected effect, judge | 79 rows as c2SA (measured on rtbT), and **12,780 more REF reads to neither** (measured, rtb1; normals 3,927). That is about 21 times the estimate of 600: 10,352 at 50+ events' junction windows, 2,428 under 50 bp | unchanged: c2SA (PASS 20 of 20 views; HG002 classes equal to C; 0 wrong-ALT added, 2 left) |
| Will not handle | the same geometry under 50 bp (16 rows, 673 rules-path REF reads); a component with only tie reads (46 of 57 rows) | S-spellers in a core under 50 bp stay partial (RJ-5); the component residual; the 7 carriers untellable at the run's `--min-baseq` |

**C. Reach on both sides** (B, and a read fitting X that does not reach where X and an X-length
candidate differ counts X): fails invariant 7 (credits reads whose bases do not carry the whole
ALT) and reverses RJ-25's MNP clause and RJ-28. Not measured further.

**Worked example** (shapes only): a 60>2 delins X with a listed SNV S at the 20th replaced base;
100 reads over the left junction: 20 carriers, 30 REF reads that reach S's base and read REF there,
50 REF reads that end before it.

```
REF  ...GATTACA c1 c2 .. c19 [c20] c21 .. c60 TGGCAT...
X    ...GATTACA AG TGGCAT...                      the junction tells X from REF and from S
S    ...GATTACA c1 c2 .. c19 [ T ] c21 .. c60 TGGCAT...

r1 carrier     GATTACA AG TGGCAT >                X: AD under every option
r2 REF, short  GATTACA c1 c2 c3 c4 c5 >           REF or S: no base tells them
r3 REF, long   GATTACA c1 .. c19 c20 c21 .. >     REF: reads c20 as REF's
```

| (RD, AD, neither), VAF | develop (floor) | rules `site` | judge today | A | B |
|:--|:--|:--|:--|:--|:--|
| 60>2 delins | (80, 20, 0) 0.20 | (30, 20, 50) 0.40 | (80, 20, 0) 0.20 | (80, 20, 0) 0.20 | (80, 20, 0) 0.20 |
| the same at 30>2 | (80, 20, 0) 0.20 | (30, 20, 50) 0.40 | (80, 20, 0) 0.20 | (30, 20, 50) 0.40 | (80, 20, 0) 0.20 |

Every r2 is a REF-side molecule (X's molecules hold no c bases), so A's "neither" at 30>2 moves VAF
by a geometry X's reads never meet: the bias the review named.

**Recommendation: B.** It is the judge's measured behaviour (every view passes, HG002 truth classes
unchanged, no wrong-ALT added), it settles #398 and the rules path's three v2 failures with no
length cut, and it removes rules rather than adding one. Its cost is named: RJ-25's REF clause
narrows to the bases a read reaches, a decided clause the operator reverses explicitly, with the
register updated. A is the register as written, at the price of a threshold-bound rule and the bias
under 50 bp.

## 4. Measurement (fixed before the runs)

- **Metric and bar, by class (rows are checks):** (1) AD: no read whose bases do not carry X enters
  AD, the fewest left (#398's reads leave the rules path's AD). (2) RD truth: at HG002 het delins
  rows with a phased SNV inside the span listed as a row, VAF within 0.05 of 0.5 on both paths (A is
  expected to miss under 50 bp). (3) The 20 patient view lines (v3 and v2) PASS on both paths; REF
  reads lost per X shape reported with v2 (REF) beside v3 (undecided), as neither oracle decides the
  reach question; the normals' AD unchanged. (4) (AD, RD) alike on both paths on every changed row
  but those a named window difference explains (today 14 of 78).
- **Data and truth:** panel, ACCESS (XS1, XS2, IMPACT), gate set, their normals; HG002 (8 arms,
  v2dbq20 and v2db) plus the delins-with-phased-SNV rows; the compound-het arm; oracles v3 and v2
  per view; `--min-baseq` 10, 20, 30, 0 on the final build. **Missing:** WES and germ (no listed
  neighbour: 0 rows); RNA (the FORTE slices: count the rows with a listed allele at the site, say
  how many); a spike-in with a listed SNV inside a delins, built only if HG002 holds under 10 such
  rows.
- **Cannot see:** whose molecules the reach reads are (only the HG002 phased rows decide it); the
  component-only residual.
- **Contract tests to keep:** `test_listed_sibling_floor_contract.py` (its floor expectation is
  part of the decision, told first), `test_haplotype_judge_contract.py`,
  `test_own_gap_option_contract.py`, `test_long_deletion_contract.py`,
  `test_unreadable_insert_contract.py`, `test_mnv_shared_confidence_contract.py`,
  `test_indel_shared_contract.py`, the rule-X tests on #297's branch.
- **Batches:** (0) on the builds in hand, no build: a trace lane on c2SA's build counting
  `ref_by_reach` calls by the other allele's kind and X's shape; `GBCMS_JUDGE_OWNGAP=all` against
  narrow; the HG002 delins-with-SNV rows on c2RC, c2RS and c2SA; the site-edge counts. (1) One
  build per option, both paths behind one switch (never for merge), patient and HG002 arrays.
  (2) One final build.

## 5. Budget and stop rule

Two review rounds. New conditions: 1. Builds: 2 plus 1 final; batches 4 plus 1 rerun. Past any
limit, or when a lane asks for a third REF rule, stop and bring the class question with one option.

## 6. Work and models

Scripts first (the round-2 scorers under `scratchpad/c2score/r2/`, extended for the shape table
and the HG002 rows; no model tokens). Sonnet builds A and B behind one switch, red-first from the
worked example on both paths, then the whole gate. The main session runs HPC after approval. A
fresh Fable or Opus reviewer per round: this page, the lanes, the diff; instance or mechanism, and
does an option break an invariant.

## 7. State

- 2026-11-02: drafted from C2 round 2, #394, #398 and the review's delins geometry. Not approved;
  nothing built; no lane run.
- 2026-11-02: batch 0 approved and run (rtb1: build ad64f930, never for merge; 1,641 runs; the
  aggregation's first job failed safely, the second ran). Results:
  `~/test/gbcms/harness/plan670/plan-first/rtb1_results.md`. Review 1 (Fable): ready with
  conditions. They are met here: A's cost corrected to 12,780 REF reads; A's long cores give back
  none of the ties; A's fix for long events is a third REF rule (the stop rule); own-gap general
  struck (decided 2026-10-27; in REJECTED.md). The site edges (16 normal rows) are round 2, shared by
  A and B. Batch 1 is not needed before the choice.
- 2026-11-02: **the operator chose B** ("tell what the read reaches"). This reverses RJ-25's REF
  clause explicitly; the register is rewritten in the build's PR. Next: one build for the rules path
  behind a switch (patient sets and HG002), then the site-edge question (round 2), the final build,
  and the delins spike-in for metric (2).

## 8. Cost log

| Step | Agent and model | Tokens | Notes |
|:--|:--|--:|:--|
| Plan page | Fable | from the usage report | read-only; one new file, not committed |
| Batch 0 (rtb1) | main session (HPC), scripts | 0 agent tokens | 1 build, 1,641 runs, 2 aggregation jobs (1 failed safely) |
| Review 1 | Fable | 191,902 | ready with conditions; one high finding (A's cost and core claim) |
