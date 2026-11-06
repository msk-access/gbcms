# RNA package 3, the rest: #293, #178, #359, #343 (#344, #345, #346)

Status: decided 2026-11-05 (operator: the lean path); built and merged in #423 (f7293618)
Rounds: 0 of 2 (Greptile only) · New conditions: 0 of 1 (class A moved) · HPC batches: 0 of 2 (none needed: measured on the local synthetic set)

Develop holds #415 (RJ-31). Three classes and one set of pins; the cheapest safe path comes first.

## 1. The classes

- **A. Which flank places a read when exon edges cut both flanks** (#293). Reads: RNA reads at a
  pure indel or delins that fills an exon of about 5 to 16 bases; each holds one flank, the one on
  its own spliced reference. Today the flank is chosen once per variant (`carrier.rs:172-184`,
  `:1870-1878`); a read that does not hold it is depth only. Two per-read cases exist: spliced at the
  left edge, from the left (`:1502-1510`); unspliced at the exon start, from the right (`:836-848`).
  Invariants: 7, equal windows at one flank (RJ-23), binning, one gate. Principle: a read is placed by
  the flank base it reads (RJ-3), whichever flank that is. Instance: the #291 review fuzz.
- **B. Every tally reads one verdict** (#359; pins #344, #345). Reads: RNA reads with a GTF at a row
  with a co-annotated sibling. The main loop (`engine.rs:1809-1821`) and the per-transcript counts
  (`:3323-3331`) apply both sibling guards (`sibling_claims_alt`, `sibling_claims_ref`); the
  junction-usage test (ASJD, `detect_asjd`, `:3823`) takes the raw verdict. Only `asjd_*` can move.
  Rules path only: the sequence judge has no sibling contest and refuses RNA. Principle: one verdict
  per read in every tally (invariant 7, RJ-5).
- **C. A catalogued editing site inside a read's judged bases** (#178; pin #346). The catalogue
  (`--rna-editing-db`, REDIportal) flags SNV rows only (`engine.rs:2104-2130`). Inside the windows
  the reads already decide (RJ-26): one mismatch is an error; a change recurring on 2 fragments and
  5% is another allele; one recurring at a fixed distance is masked by the guard (`carrier.rs:57-60`).
  An edit inside an insertion's bases has no coordinate, so no catalogue reaches it (#346).

## 2. What already touches it

- **Rules (read to their ends):** RJ-3 (the flank base; the short-piece limit, #409), RJ-17, RJ-19
  (read on through a short exon), RJ-23, RJ-27, RJ-31, RJ-5 and RJ-25 (siblings), RJ-26.
- **Rejected, not proposed again:** a splice-aware Phase 3 and hard-filtering junction-adjacent
  variants (REJ-20260924-001); a junction-keyed BAQ rule (REJ-20260925-002); `rna_antisense_depth`
  as "reads that reached classification" (REJ-20260925-003).
- **Tickets:** #289 and #291 (the flank tests; one reverted), #292 (its limits #409, #408), #385,
  #309 (#359's origin), #94 and #114 (the pins' origin). Pins in place: `test_rna_*_contract.py`
  (five files), `test_antisense_depth_contract.py`, `test_asjd2_contract.py`.

## 3. Options

**Class A, worked example.** A 2-base deletion at bases 3 to 4 of an 8-base exon (`E`; deleted
bases `d`; `L` the previous exon on the mRNA, `R` the next).

```
  mRNA             L L L L L L [E E d d E E E E] R R R R
  (i)   spliced in L L L L L L [E E d d E E              holds L only
  (ii)  exon start             [E E d d E E E E] R R R R  holds R only
  (iii) through    L L L L L L [E E d d E E E E] R R R R  holds both
  (iv)  pre-mRNA   intron . . .[E E d d E E E E]. . .     reads the genome
```

| Option | Rule | Inv. 7; equal windows | Conditions | Will not handle |
|:--|:--|:--|:--|:--|
| **A0. One flank per variant**, documented | as today: (i) or (ii) is depth only, by the edges | passes | 0 | one read group at every such row; the loss unmeasured |
| **A1. Per-read flank** (#293) | a read reads from the flank it holds on its spliced reference; (iii) from the left; (iv) as the genome | passes: bases decide at either flank; a REF and an ALT molecule of one shape read from the same flank | +1 (replaces the variant-level `cut_right` test) | RJ-3's short pieces (#409); #385's shape |
| **A2. The judge** (#331, 6.8.0) | the per-haplotype judge reads the whole read on its spliced reference: no flank to choose | passes | 0 now | 6.7.0's rules path keeps A0 |

Recommendation: measure first. Build A1 only if the loss is visible; else A0 now, A2 in 6.8.0 and
#293 moved. A1 is larger than it looks: it touches every spliced read's layout, and the #289
right-flank attempt was reverted once.

**Class B, worked example.** Row X (a deletion), sibling S (an MNP inside X's window). 10 spliced
fragments are ALT for X by bases, 4 claimed by S. Main AD 6; ASJD `n_alt_total` 10 today.

| Option | Rule | Inv. 7 | Conditions |
|:--|:--|:--|:--|
| **B1. ASJD applies both guards** | the two calls the per-transcript loop makes, made in `detect_asjd`: `n_alt_total` 6 | passes | 0 |
| **B2. Leave it, document** | the raw split | fails in effect: a sibling's molecule votes for X's junction | 0 |

Recommendation: B1; the measurement sizes its note. Pins #344 and #345 land with it.

**Class C, worked example.** A 3-base insertion `TAC`; a catalogued A>G site two bases past it,
read `G` by 6 of 40 carriers. RJ-26: a shared change in the flank is a phased neighbour; the carriers
keep ALT. Inside the event (#346's table): one edited carrier stays ALT; two or more at 5% are
another allele; an edit on the one telling base leaves the read without evidence.

| Option | Rule | Inv. 7 | Conditions | Will not handle |
|:--|:--|:--|:--|:--|
| **C0. No mask; document, pin** (#346) | the reads decide: RJ-26 and the guard | passes | 0 | an edit under the guard's floor on the telling base (depth only: right by bases) |
| **C1. Catalogue mask in windows** (#178) | a catalogued position is masked in RNA mode, the row flag kept | passes as a mask, never a credit | +1 (the first rule keyed on an external list) | inserted bases; uncatalogued sites; a catalogue on the wrong build |

Recommendation: C0 now; C1 to 6.8.0 unless the measurement shows a residual the guard misses.

## 4. Measurement (fixed before the runs)

One build of develop with three **count-neutral diagnostics** (trace lines, no column); one HPC
array over the 58 FORTE STAR slices (aggregates only), bars 20 and 30, GTF on:

- **A.** Rows: the #292 indel arm's generated rows kept to exons of 5 to 16 bases, plus the truth
  cohort's rows there; the #292 oracle classes them clean, low, het. Metric: reads that hold the other
  flank and are depth only for want of the chosen one, by class (i) to (iv), as a share of the row's
  informative reads. Local first: the synthetic set (`~/test/gbcms/harness/plan670/m292/syn`, a
  short-exon stratum added: hets at 0.5 and hom-REF rows, both pair orientations) gives het VAF error
  and false ALT. Bar to build A1: the synthetic het VAF error at short-exon rows exceeds the plain
  stratum's by more than the seed-to-seed spread, and the FORTE depth-only share at such rows exceeds
  the control band's spread (#292: 0.158%). Below it: A0 and A2.
- **B.** Rows: the truth cohort's rows with a sibling (`count_sibling_rows.py`, the R1 scratch
  folder). Metric: fragments ASJD counts that the guards exclude; rows where `asjd_flag` or
  `asjd_pval` moves. Local first on the local FORTE sample; HPC only if it has no sibling rows.
- **C.** Needs REDIportal for GRCh38 staged on HPC (public; the operator stages it). Metric: counted
  rows with a catalogued site inside the windows; of them, rows the guard already masks; reads a mask
  would move. Bar to build C1 in 6.7.0: such reads at any truth row, or at over 1% of RNA indel rows.
- **Truth and data:** RNA only (FORTE STAR, the synthetic set). DNA is unchanged by construction (no
  GTF, no edges); the final batch proves it by identity. HG002: none, it is DNA. Missing: a second
  aligner (HISAT2), a second cohort. **Cannot see:** whether a class (i) read at a DNA-absent row is
  a novel-junction piece (#409); edits below the guard's floor.
- **Contract tests kept:** the files above, the splice-skip exclusion, the four counting invariants on
  every RNA test. **Final batch:** the chosen build, the same arrays; identity on every row and column
  the change does not name (`asjd_*` only for B1).

## 5. Budget and stop rule

Two review rounds for the package. New conditions: at most 1 (A1). A fix keyed on a row, a fixed
exon length or a catalogue entry stops the work and brings the class question back with one option.
HPC: 2 batches plus one rerun. If the diagnostic batch is not aggregated within one week, A and C
move to 6.8.0 and B1 lands on the local measurement.

## 6. Work and models

1. **Pins first (tests only): #344, #345, #346.** Sonnet, red-first, the RNA runner's invariants.
   Re-run #346's table on develop before pinning: its numbers predate RJ-31; a number that differs is
   a decision point, not an edit. Pins that might show behaviour that is not the decided one: #344's
   exon-skipping transcript (a read with no junction fits both transcripts; which DP holds it is not
   in the register); #345's no-GTF case (a no-op with a warning: the counts are the contract). Report
   both to the operator before the tests go green.
2. **Diagnostics and the batch.** Haiku for the scans (sibling rows, short exons in the GTF); Sonnet
   for the trace lines and the aggregation script; the main session submits the array.
3. **B1.** Sonnet: the two guard calls in `detect_asjd`, a red-first test (ASJD and main AD agree at
   a row with a sibling), a note under RJ-5, a CHANGELOG line sized by the measurement.
4. **A1, if the bar is met.** Opus for the layout change (class reasoning in `carrier.rs`); Sonnet
   for tests; Fable or Opus for the one review. Else #293 to 6.8.0 with the numbers.
5. **C0.** Sonnet: one doc sentence (the catalogue flags SNV rows; the reads decide inside windows;
   inserted bases are out of reach) and #346's pin. #178 closes or moves with the numbers.

Each brief points to this page; an agent stops at the end of its phase.

## 7. State

Draft 2026-11-05; decided and built the same day.

**Decided 2026-11-05 (operator): the lean path.** 6.7.0 takes the three pins (#344, #345, #346; tests
only) and class B for #359 (ASJD makes the same two sibling calls as the main counts; 0 conditions;
only `asjd_*` may move), measured on the local synthetic set with no HPC batch. Classes A (#293) and
C (#178) move to 6.8.0 (#420): under the judge default, A's flank question goes away, and C's mask
would be the first rule keyed on an outside list. The build is on `fix/rna-pins-asjd-guards`.

**Done (#423, merged 2026-11-05):** the three pins; B1 (`detect_asjd` makes the two sibling calls; a note under RJ-5);
RJ-32 decided for #344's case (a junction-free read counts in the DP and RD of each overlapping spliced
transcript; a single-exon transcript gets no entry). Synthetic set, measured against develop before and after
class plan B (#419): 0 of 1,072 rows move without siblings; with one sibling per row 12 of 2,144 rows move, in
`asjd_*` only. #343, #344, #345, #346 and #359 are closed. #293 and #178 are on 6.8.0.

## 8. Cost log

Rounds 0 (Greptile's one doc finding, fixed); conditions added 0, removed 0; HPC batches 0; agents: this page
(Fable) and the build (tokens in the trial record, `docs/harness/plan-first.md`); days from the plan to the merge 0.
