# #333: RNA under the sequence judge
Status: draft, revision 2 after the review (2026-11-17; operator: "What did the rules branch do and what does the
research tell then plan it and do an adverarial review of the plan and then let us decide"; "we dont just hold
to community standards we try to go above it as we can test with the data at hand"; "We can route depending on
the cli call") · plan-first (`plan-first.md`); decision protocol (`.agents/rules/decision-protocol.md`) · base
develop 6ffde565 · in 6.8.0 as an opt-in (R1 (a), 2026-11-17; `335-phase2.md` §3.9-§3.10)
Reviews: 1 (Fable, `scratchpad/speed2/rna333/review/review.md`: revise, 8 findings, all taken) · New count
conditions: 0 of 11 declared · Builds: 0 · HPC batches: 0 (testing is counted, not capped)

This page is longer than one page: #333 is the cycle's largest build and carries eleven options.

## 1. The class

The judge refuses RNA (`engine.rs:444-450`, `cli.py:854-859`), and a spliced read gives it no view
(`decide.rs:1270-1277`: depth only). RNA reads differ from DNA reads in six ways the judge must handle:
they are spliced (CIGAR `N`); they carry clips and short overhangs at exon edges; they are admitted at
`--min-mapq 1` with STAR's multi-mappers; antisense reads are judged for a tally but not counted; a row's
reads are called again by the per-transcript counts and ASJD; and many rows are SNVs at sites with
catalogued RNA edits.

**The principle:** a spliced read is judged by its own bases on the reference spliced along its own
junctions, with the rules path's answer kept wherever no decision changes it. The field's practice is a
floor, not a ceiling (decision protocol §5): wherever the data at hand can test a better rule, the plan
measures it beside the field's.

**What changes where.** Under `gbcms rna --read-judge sequence` (the opt-in) the judge counts RNA by the
routes in E10. Under `gbcms rna` with its default (`--read-judge rules`) and under `gbcms dna`, no output
changes in 6.8.0 until the freeze decides RNA's default class by class.

**Invariants:** 7, pass with a limit: the read's junction set is where its bases sit, taken as given, never a
splice the judge chooses; where STAR placed the junction wrongly (#385, #409) the judge inherits it. 2, one
quality contract: the same BAQ (`effective_quals`, `baq_applies`) in every view and on both paths. 1, binning
invariance: a row's candidate set, its frames and its caches depend on the row and its reads only (E0). 6:
RNA columns gated.

## 2. What already touches it

Sources: the rules-path map (`scratchpad/speed2/rna333/rules_path.md`, Opus), the outside-tools survey
(`rna333/survey.md`, Fable), M12 (`335-phase2.md` §3.10), the cache audit (`335-source-fixes.md` §4.2), the
review (`rna333/review/review.md`, Fable).

- **The splice frame exists on the rules path.** Windows on the reference spliced along the read's junctions
  (`carrier.rs:1531-1611` `spliced_windows`, `:1618-1667` `splice_segments`; RJ-19, 2026-10-04, amended
  RJ-27). M12: the band scorer on such a linear frame equals brute force (1,656 of 1,656 reads) at DNA's
  cost, and splice edges in one graph break the band (7 of 36 reads at intron 60; 6x-114x cells). M12 tested
  the scorer only: one SNV at an acceptor, one intron, 220-base exons; not frame construction, short exons,
  two-junction reads, events at a junction or a row with several frames (E0).
- **REJ-20260924-001** (`REJECTED.md:440`): scoring a read across a junction against a spliced haplotype
  (B2b of #94), vetoed 2026-09-24 as not demanded: every covered junction-adjacent indel classified
  correctly through the rules path's structural paths. **What changed:** that was a fallback for the rules
  path; the judge has no structural path, so without a frame every spliced read is depth only, and 41.3%
  of the reads at FORTE's truth rows are spliced near their row (M12's sizing, 46,787 reads). The demand
  is now the judge's whole RNA input.
- **Outside tools** (survey): GATK, Strelka2 and DeepVariant split reads at `N` and lose the junction's
  flank. ABRA2's RNA mode builds one reference per junction combination and found 27 of 27 EGFR exon 19
  deletions in RNA, but its gain came from choosing among junctions and moving reads, which this plan
  forbids; the precedent here is the spliced reference, not the choice.
- **The rules path's CIGAR-shape rules:** E2's table.
- **The three calls.** The main counts, per-transcript counts and ASJD each classify a record
  (`engine.rs:1859`, `:3392`, `:4006`), each with its own admission: the anchor or read-from-right rule and
  the antisense tally; junction compatibility within 5 bp (RJ-32); no anchor gate for ASJD (ASJD-2, RJ-5).
- **Admission and doubt.** `--min-mapq 1` with the `NH:i:1` rescue (`rna.rs:22-75`); RJ-29's doubt test
  needs XS, which STAR does not write; the 2026-10-03 decision kept multi-mappers (unique-only lost 12 ALT
  reads, 10 of them at PIK3CA E545K).
- **The gap model.** RNA's PairHMM gap defaults (RT stutter) reach only the rules path's fallback; the judge
  scores with a fixed gap model (`score.rs:21-24`); FORTE's homopolymer slippage matched IMPACT DNA's rate
  (2026-10-27, `docs/reference/rna-splice-handling.md:115-118`).
- **Decided before this plan:** BAQ in RNA (RJ-31: an `N` lowers nothing); editing sites flag only (C0,
  2026-11-05); strandedness (#114, #309, #345); the decomposed twin deprecated in 6.7.0 and removed in
  6.8.0 (2026-11-05).
- **Open tickets in the class:** #385 (a deletion of an exon's first bases; STAR writes 21+ bp as `N`),
  #409 (RJ-3: a short anchor at an annotated acceptor next to a novel one), #293 (short exons), #178 (editing
  sites), D9 (the spliced cache's key), D6 (F2-a: the judge's verdict map removed; `335-source-fixes.md`).
- **REJECTED.md, also:** REJ-20260925-002 (no N-penalty BAQ rule at annotated junctions); REJ-20260925-003.
- **Validation so far:** FORTE truth (33 samples, 94 rows: 68 SNVs, 26 others); T9 (978 exon-edge rows); C1
  (7,224 delins rows); #292's 58 STAR slices, 40 with matched IMPACT DNA; a synthetic set with modelled
  STAR-like CIGARs at overhang 3 and 5. No judge run on RNA exists.
- **Tests that pin the refusal** (their expectations change, with an operator notice first):
  `tests/test_alignment_backend.py:288` `test_cli_rna_refuses_the_sequence_judge` and `:298`
  `test_engine_refuses_the_sequence_judge_in_rna_mode` (the opt-in lifts the refusal);
  `tests/test_ref_reach_contract.py:622` `test_a_routed_row_in_rna_mode_is_counted_by_the_rules_with_one_warning`
  (changes only if RNA's default routing changes at the freeze).

## 3. Options

Each option gives the field's practice where one exists and, beside it, an option above it where our data
can test it.

**E0, a row's candidate set when its reads have several frames (new).** The site is one per row and a
function of its reads (`judge/mod.rs:168-181`); the found-allele floor is a share of the site's reads
(`GUARD_REF_SHARE * reads.len()`, `mod.rs:588-602`). (a, recommended) decide the candidate set once per row
on the genome, from all the row's counted reads, spliced and unspliced (contexts, found alleles, site
states, as the rules path does); each frame only respells that set on its spliced reference; the floor's
denominator is all the row's counted reads; caches per (row, frame). Binning: pass (per row). Conditions: 1
(the respelling). Limit: an allele that only spliced reads spell across a junction is not found (E4). (b) a
site per frame: an allele on 3 spliced and 3 unspliced reads can fall under the floor in both. (c) one
pooled site with genome positions mapped into each frame.

**E1, the frame (D10).** (b, recommended) one spliced linear reference per junction set, built with
`splice_segments`; the read's own junctions fix its frame; unspliced reads use the genome frame.
Invariant 7: pass with the limit in §1. Binning: pass. Conditions: 1. (a) splice edges in one graph: a new
stand test and an edge-following reach, up to 114x the cells per spliced read.

**E2, the rules path's CIGAR-shape rules.** Each is new code in the judge; a ported rule is re-derived on
the judge's inputs (protocol §2). Every rule kept withholds a call and never credits ALT (the safe side of
invariant 7). One line per rule:

| Rule (rules path) | Form in the judge | Register | Invariant 7 | Keep | Pinning test | The oracle scores it |
|:--|:--|:--|:--|:--|:--|:--|
| Splice-skip triage: a read whose `N` covers every differing position counts nowhere, DP included (`variant_checks.rs:266-342`) | front test, already in front of the judge (`engine.rs:2924`) | 6.4.0 (#94), RJ-17 | the read holds no base at the row | permanent | `test_rna_splice_contract.py` | no base, no call: agrees |
| A splice through the differing bases: depth only (`carrier.rs:774-807`) | front test (new) | RJ-17 | withholds | temporary: the frame alone should give depth only (both candidates spell the same) | `test_rna_partial_span_contract.py` | agrees if both spell alike |
| A deletion written beside a splice: depth only (`carrier.rs:821-829`) | front test (new) | RJ-17 | placement doubt | temporary, measured against the frame | `test_rna_splice_contract.py` | disagrees by design |
| A partial exon piece of a deletion: depth only (`carrier.rs:3416-3471`, `partial_span_piece`) | front test (new) | RJ-17 | placement doubt | temporary, measured | `test_rna_partial_span_contract.py` | disagrees by design |
| A junction entering the event: ALT only over "REF spliced at the read's own junction" by one bar-quality base (`carrier.rs:1250-1262`, `:1449-1497` `own_alignment_weight`) | a candidate sequence in the site (new) | RJ-19 | by bases | permanent | 4 alternative-donor and NAGNAG tests in `test_rna_group3_contract.py` | agrees |
| The exon-start cut (`carrier.rs:842-869`) | a read-view rule (new) | #291 | placement (annotation) | temporary, measured | `test_rna_exon_start_run_contract.py` | disagrees |
| Clip edges and `EDGE_SLACK` 5 (`carrier.rs:186-226`) | already in the judge's read view (`decide.rs:1278`, `mod.rs:482-485`) | RJ-18 | — | permanent | `carrier.rs:4449`; `test_rna_group3_contract.py` | agrees |
| No admission by clips (`engine.rs:1960-1966`) | admission (new) | RJ-18 | — | permanent | `test_clip_carrier_contract.py::test_rna_reads_are_not_admitted_by_clips` | — |
| One reading flank per row (`carrier.rs:1950-1975`) | **dropped**: a whole-read judge has no flank | #293 | — | an accepted difference | measured on the #291 fuzz shapes and a short-exon stratum | — |
| A 50+ bp deletion's spliced read: depth only (`junction.rs:715-720`) | front test (new) | RJ-24, RJ-17 | withholds | until #385 | `test_rna_splice_contract.py::test_large_deletion_near_junction` | — |

(a, recommended) keep the table's "keep" column for the opt-in; each "temporary" rule is measured
against the frame alone and dropped only by its own decision. Conditions: 7. (b) drop the temporary rules
now: each a count change with its own register entry.

**E3, unspliced short overhangs (3.3% of reads at FORTE rows).** The field: no tool counts them on the
junction; GATK clips them away. (b, the target, above the field) a second candidate placement on a
junction's spliced reference, with three limits: the junction carried by another read's `N` or the
annotation, never inferred; the spliced placement must win by the one-confident-base margin, else the
overhang's bases give no vote; bases matching the intron stay unspliced. Invariant 7: the junction is
another read's or the annotation's, so this is a placement the judge weighs by bases: a limit to state.
Binning: pass (junctions from the row's reads). Conditions: 1. (a) judged on the genome, as today, until (b)
passes its bars.

**E4, contexts across an intron.** Today the guard stops at the first junction (`carrier.rs:3131-3146`),
siblings and rows in reach come from prep on the genome within 300 bp (`normalize/engine.rs:282-318`), and
masks move into a frame only where they sit in place. (a) the same on the judge (E0 (a)). (b) a listed row
on the next exon becomes a context for the spliced reads that reach it: its evidence is the row's spliced
reads only (the bin's fetch is genomic; a wider fetch would collide with S5). Binning: pass if found from
the row's reads only. Decided after M13.

**E5, SNV rows under the opt-in.** The rules read an SNV from its one aligned base (RJ-31,
`variant_checks.rs:424-458`). (a) judged whole, DNA's method, measured as its own class (68 of the 94 truth
rows). (b) the one-base read kept unless RJ-25 routes the row (an SNV beside a listed indel or delins); also
faster on most RNA rows. Both are measured; E10 sets the route. Under the default `--read-judge rules`, RJ-25
stays off in RNA (`engine.rs:480-481`) until the freeze.

**E6, one verdict per record per row (a speed item, not a count option).** The judge is pure, and its
`judged` map already gives the three calls one verdict where it has room. Once D6 (F2-a) lands, the reuse
keys on the record's identity, its effective qualities, its clip edges, its frame and the rules. Each view
keeps today's admission. (b), the survey's partition (per-transcript and ASJD counts as a split of the main
count's fragments), changes ASJD's counts: its own decision, measured.

**E7, RNA's doubt signal for reads that fit nothing (RJ-29).** (a, recommended) `NH` above 1 where the tag
is present; no signal where it is absent, said once in the log; `NH` joins the read view's tags
(`haplotypes.rs:97-103`) and so the cache key. MAPQ is not used: its scale is the aligner's (STAR's 255;
rewritten BAMs carry 60). Bar: the PIK3CA E545K reads keep their calls; a paralog stratum measures the
gain. (b) none in 6.8.0, with the residual named.

**E8, base quality near junctions.** The field: htslib skips BAQ on any read with an `N`; GATK, REDItools and
others use position rules because they score on the unspliced genome; the rules path lowers nothing at an
`N` (RJ-31). (a, above the field) measure the judge's error rate by distance to a junction after spliced
scoring, on matched DNA's germline hets and the STAR-aligned synthetic set, and add a rule only where the
data shows a bias. (b) no rule, unmeasured.

**E9, editing sites (#178).** The field: callers remove or tag catalogued sites; given-site counters do
nothing. The rules path flags an SNV row whose catalogued edit matches (C0). (a, above the field) keep the
flag, and measure with matched DNA what no tool reports: how often a catalogued edit sits inside a row's
read reach, and whether hyper-edited reads fall to "neither" through the fit rule; an A>G in RNA that the
matched DNA lacks is an edit. (b) the flag only, unmeasured.

**E10, routing by the CLI call (operator, 2026-11-17: "We can route depending on the cli call").** The engine
routes per row: `judge_of` (`engine.rs:1300`) sends a row to the judge where prep marked it (`judge_routed`,
RJ-25), else to the call's own `--read-judge`; `gbcms dna` and `gbcms rna` each have that option
(`cli.py:391`, `:791`). (a, recommended) a route per CLI call and row class (SNV, MNP, indel under 50,
complex, 50+ deletion), set in prep as RJ-25's mark is. Until the freeze it acts only under the opt-in; at the
freeze each RNA class's default route is set by its bars, and a class that fails stays on the rules, said
once in the log. RNA's fronts (E1-E3, E7, RNA admission) run only under `gbcms rna`; the scorer is shared.
Binning: pass (fixed in prep). (b) one switch per call.

**E11, the gap model on RNA.** (a, recommended) DNA's gap model, by the evidence that FORTE's slippage
matches IMPACT DNA's; RNA's RT-stutter defaults stop reaching the counts under the judge. (b) an RNA gap
model, measured.

**E12, a repeat-aware gap cost in the judge's scorer (planned 2026-11-23; shared with DNA, so not in the
opt-in).** The research is `333-gap-survey.md` (Fable; public sources and code only).
- **The class.** The scorer has one cost order everywhere (`score.rs:21-26`): a gap of k bases costs 6+k, a
  confident mismatch 4, the margin 4. In a homopolymer of 8 bases or more, a 1-base slip is as likely as a
  confident base error or more likely (survey §3.2). So in a long run the scorer takes "a confident base is
  wrong" before "the run is 1-2 bases longer". Read check J (§7, 2026-11-22) measured two effects on HG002 RNA:
  1. Missed ALT. An ALT molecule with a 1-base slip in the run the substitution makes. Row 13 (C>A; A9 C A2
     becomes A12): an A13 read scores ALT −7, REF −11, and `indel_unexplained` (`decide.rs:939-951`) sends it
     to Other. The judge counts ALT for 106 of 1,729 ALT+1 reads and 54 of 458 ALT−1 reads.
  2. False ALT. A REF molecule with a slip. Row 45 (G>A; A3 G A10 becomes A14): an A3 G A12 read scores
     ALT −4 (with a listed AA insertion), REF −8 (a 2-base gap), so ALT by the margin. 6 reads likely, 54
     unproven.
- **What the other tools do** (survey §1, §4):
  - The likelihood tools (GATK, DRAGstr, Octopus, Platypus, Strelka2, HipSTR) make a slip cheaper in long
    runs and weigh each read; none has a class like `partial_alt`.
  - Lancet2's genotyper keeps every gap strict (12+3k against a mismatch of 4) to stop "allele bleeding". Its
    two-tier gap is only in its assembly step.
  - GATK and DRAGstr take the cost from the read's own bases, so a slip costs the same on every haplotype.
    Octopus and Platypus take it from the haplotype, which gives the allele with the longer run an advantage.
  - On the judge's scale (cost = 4 / Q_mm × Phred, Q_mm 25), a 1-base slip on a PCR library costs about 5, 4,
    3, 2, 2, 1 at L = 6, 8, 10, 12, 15, 20 (PCR-free 6, 5, 4, 3, 3, 2); today it costs 7 at every L. Each
    further unit costs about 1.6. Contractions are 5-14 times as frequent as expansions with PCR. No public
    table by tract length exists for human RNA-seq.
- **Constraints from our rules:**
  1. At an indel row the length is the allele (RJ-6, RJ-27, the decided stutter rule of #92). Inside the
     row's discrimination window a gap keeps its full cost. Worked example: a 1-base deletion row in an A12
     run. A REF read scores REF 0 and ALT −7 today, REF by 7. With a slip cost of 2 inside the window, REF
     wins by 2, under the margin of 4: every REF and ALT read of a long run becomes depth only.
  2. The cost comes from the read's own tract at the gap (the GATK and DRAGstr form), never from the
     haplotype. In a tract that both haplotypes share, a slip then costs the same on both, and the bases decide
     (invariant 7).
  3. The mismatch cost, the margin, the quality gate (invariant 2) and the masked-base rule do not change.
     Only a gap inside a read tract of period 1-3, at or above the measured threshold, changes.
  4. One fixed table in code: no environment switch, and not learned per sample. A per-sample model learns
     an MSI-high tumour's real slips as error, and a panel or one RNA sample has too few long tracts (DRAGstr
     needs 50 loci per cell). A second table for `gbcms rna` is plumbed intent (invariant 3), chosen only if
     M1 shows that RNA differs from DNA.
  5. The cost depends on the read and the row's window only, so bin geometry does not touch it (invariant 1).
- **Options:**
  - (a, measure first) The read-tract cost of constraint 2, outside the row's window; `indel_unexplained`
    unchanged. Row 45: ALT −4, REF about −4 (a 2-base slip in A12), so no class wins by 4 and the read is depth
    only, not ALT. Row 13's ALT+1 reads stay Other. In a tract the row does not change, today's costs already
    cancel; there (a) changes only whether a read with a slip passes the fit floor (`decide.rs:361`), so fewer
    reads are depth only.
  - (b) (a), and at a substitution row a 1-unit slip in the run the substitution makes is not "another
    allele" where the read's confident base at the row decides. This amends the scope of `indel_unexplained`
    and RJ-7 at substitution rows; the indel-row rule of #92 does not change. Row 13's A13 read: ALT about −2,
    REF about −6, so ALT. A REF+1 read keeps its confident REF base and stays REF.
  - (c) No change (Lancet2's position). Row 45's false ALT stays.
  - Open design point: is a listed neighbour row's window (row 45's AA insertion) strict too? Measured on the
    context rows of read check J.
- **Measurement, before any table is written:**
  - M1. The slip rate per spanning read by tract length (6-20), period (1-3) and direction, at truth-REF
    tracts (no truth variant in the tract or its flanks), per data type: HG002 RNA (GIAB, public); HG002 DNA if
    the operator approves the download; FORTE RNA and IMPACT tumour DNA at the internal truth's `d_ref` rows
    (aggregates only; a new kind of job on clinical BAMs, so the operator approves it first). The table is the
    measured rates by the survey's conversion; the survey's values are the check.
  - M2. A prototype build per option on `r333g3`: HG002 RNA (all hg333 rows, and the 242 rows of read check J),
    the FORTE internal truth (`truth333_agg`), the 18 DNA gate arms, IMPACT, ACCESS (UMI consensus: fewer
    slips) and HG002 DNA. The four count relations on every count.
  - M3. Rust tests: a deletion row in A12 (REF and ALT reads keep their calls); a substitution beside a shared
    A12 tract, read with a slip (the cost cancels); row 45's shape (depth only); row 13's shape (Other under
    (a), ALT under (b)); a masked base at the row (depth only).
  - M4. Speed: the read's tract lengths are found once per read, in time linear in its length; timed on every
    data type (the speed floor).
  - Bars: on truth rows, false ALT and missed ALT do not rise on any data type, and one of them falls; at
    indel rows inside tracts, RD and AD equal `r333g3`'s (constraint 1); the depth-only share (DP − RD − AD)
    reported per data type (cheap gaps make more ties); VAF at hom rows nearer 1.0.
- **Risks:** MSI-high tumours (at indel rows the strict window keeps their real slips apart; at substitution
  rows (b) credits more ALT molecules with a slip, which is right for the given allele); PCR against PCR-free
  (2.5-9 times the slip rate; IMPACT and ACCESS are PCR capture, HG002 GIAB DNA is PCR-free); UMI consensus
  (ACCESS has fewer slips, so a PCR table is permissive there; that moves only ties); more ties in disordered
  repeats ("allele bleeding").
- **Protocol:** a new Open row in the register; the effects map is the classes above (substitution rows in or
  beside a tract of 8 or more, indel rows in tracts (unchanged by design), reads under the fit floor by a
  slip); several tools (the survey); every data type (M2); worked examples (rows 13 and 45, the A12 deletion);
  measure first (M1). The change is shared with DNA, so #333's DNA identity gate cannot hold for it: it needs
  its own ticket and branch, and the operator's decision after M1-M3.

**Declared conditions (the ceiling): 11.** E0 the respelling (1); E1 the frame selection (1); E2's new forms
(7: four front tests, one candidate, one view rule, one admission); E7 the doubt signal (1); E10 the route
table (1), which replaces `mode != "rna"` in RJ-25's routing at the freeze. E3 (b), if chosen, adds 1.

## 4. Measurement (fixed before the runs)

**Accepted differences.** Under the opt-in, differences between the judge and the rules on RNA are the
measurement, not a stop. On DNA (both judges) and on RNA's default path, no difference is accepted.

**Truth, in order of independence from the judge:**
1. **Molecule truth from STAR itself:** reads simulated from transcripts with known alleles, aligned by STAR at
   FORTE's settings (annotated junctions from 1 base), with strata for short exons, a novel acceptor next to an
   annotated one, a deletion of an exon's first bases, a paralog, hyper-edited reads and 3-8 base overhangs.
   This decides E1-E3 and E7.
2. **Matched IMPACT DNA** for 40 of #292's 58 FORTE slices: germline hets (an allele at about 0.5 in DNA
   should show in RNA, allele-specific expression aside), somatic rows (DNA ALT reads against RNA ALT credit),
   edits (an A>G in RNA that DNA lacks).
3. **FORTE truth** (94 rows; only 26 non-SNV rows, so MNP, complex and 50+ deletion classes have almost no RNA
   truth: named as a gap).
4. **The bases-on-the-spliced-reference oracle** (pysam; the run's `--min-baseq`): a second check on real data;
   it shares E1 (b) and E3 (a) with the judge, so it finds scorer bugs, not frame errors, and disagrees with
   E2's front tests by design.

**Per option:**

| Option | Metric | Bar | Truth | Set |
|:--|:--|:--|:--|:--|
| E0 | found alleles recalled; false contexts; memory | (a) loses no found allele the rules path finds; peak within 2x `fin473` | 1, 3 | synthetic, tier FORTE |
| E1 | ALT recall and false ALT per class | ALT recall at or above the rules', false ALT at or below | 1, 2, 3 | all |
| E2 | for each temporary rule: rows where the frame alone differs | each difference explained, with a worked example | 1, 4 | synthetic, T9, C1 |
| E3 | intron-against-exon ties; ALT recall at overhangs | (b) adopted only if recall rises and false ALT does not | 1, 2 | synthetic overhang stratum, #292 slices |
| E4 | rows with a listed row across an intron within reach (M13); state meetings across a junction | decided by M13's count | 1, 3 | tier FORTE, 10 FORTE BAMs |
| E5 | SNV ALT recall and false ALT, (a) against (b) | the better of the two on truth 1 and 2 | 1, 2, 3 | all |
| E6 | calls per record per row; CPU | one call; CPU not above the rules' | — | gate RNA arms |
| E7 | fits-nothing reads by `NH`; E545K | E545K reads keep their calls | 1 (paralog), 3 | synthetic, FORTE |
| E8 | error rate by distance to a junction | a rule only where a bias shows | 1, 2 | synthetic, #292 slices |
| E9 | catalogued edits within reach; hyper-edited reads to "neither" | measured and reported | 1, 2 | #292 slices |
| E10 | each class's bars above | a class routes to the judge at the freeze only if it passes | 1-3 | all |
| E11 | indel ALT recall at homopolymers | DNA's model at or above the RNA defaults | 1, 3 | synthetic, FORTE |
| E12 | slip rate by tract (M1); false ALT and missed ALT on truth rows; indel-row counts in tracts | false ALT and missed ALT not higher on any data type, one of them lower; indel rows in tracts unchanged | 2, GIAB (HG002) | every data type, DNA included |

**Arms:** the gate's RNA arms (the tier's 33 FORTE runs, the RNA fill-out, without the GTF, the whole RNA
log); FORTE truth; T9 (978); C1 (7,224); #292's 58 STAR slices; the synthetic set; a `--min-baseq` sweep (0, 10,
20, 30; #331 item 4); BAQ on and off under the judge; a **DNA identity gate** (shared judge code changes: all
18 DNA arms, 0 differences; a DNA read with `N` goes from depth only to judged, so DNA arms holding such
reads are counted apart); binning invariance on RNA (the gate's `bins1` kind: every row its own bin); the
four count relations on every count.

**Two bar lists.** *Opt-in (to ship `--read-judge sequence` for RNA in 6.8.0):* the two refusal tests changed
after the operator's notice and J1's other RNA tests resolved; the DNA identity gate; binning invariance on
RNA; every output path follows the verdict (per-transcript, ASJD, antisense depth, splice-spanning,
observations); no failure on CRAM, unpaired, hard-clipped or amplicon input. *Default, per class at the
freeze:* the class's rows of the table above; the speed floor (faster than the C++ GetBaseCountsMultiSample
on the 10 FORTE BAMs of the HPC speed batch pf76: 1,000 rows by 1 BAM and 5,000 rows by 10 BAMs; peak memory
within 2x); the HPC gate.

**What this cannot see:** a second aligner (no HISAT2 data) or cohort; strand protocols other than FORTE's;
STAR's wrong splices where no simulated stratum holds them; the non-SNV classes' real-data truth (26 rows).

**HPC:** M13 (10 FORTE BAMs, counts only, one job); the STAR alignment of the synthetic set (one array, its
size set by the strata); the #292 slices under the judge (58 tasks); the speed cells (pf76's 10 BAMs); the
HPC gate.

## 5. Budget and stop rule

Testing is counted, not capped. Stop and go back to the operator on a difference §4 does not accept (DNA, or
RNA's default path), on more conditions than the 11 declared, or on a change of design or decision.

## 6. Work and models

1. Before any build: E0's choice; the extended M12 probe (frame construction, a short exon, a two-junction
   read, an event at a junction, a row with several frames: Opus); the STAR-aligned synthetic set and the
   oracle (Opus designs, Sonnet builds); M13 (HPC); D6 (F2-a) and D9 (the spliced cache) landed.
2. E0, E1 and E2 (a) behind the opt-in (Opus writes, Fable reviews the diff, Sonnet gates), with the two
   refusal tests changed after the operator's notice.
3. The arms; E3, E5, E7, E8, E9 and E11 measured; the decisions to the operator with worked examples.
4. The speed items on spliced frames (`335-phase2.md`: S3b, S3c, K) and the speed cells on HPC; the freeze's
   routes (E10).

## 7. State

- 2026-11-17: research done (rules-path map, survey, M12, FORTE sizing); revision 1; review (Fable, revise, 8
  findings); revision 2.
- **Decided (operator, 2026-11-17: "I believe we should do it all then measure it and then do another review to fix
  thing to the source and class if anything needs fixing"):** E0 (a), E1 (b) (D10), E2 (a), E10 (a), E11 (a), and the
  two refusal tests change when the opt-in lands. Every measured option (E3 (b), E5 (a) and (b), E7 (a), E8, E9) is
  built as a switch so all are measured on the same build; then a second review; then fixes at the source and the
  class. **Local first** ("we are running currently out of github actioins minutes ... we should do most of it
  local untill we are sure"): no push and no PR until the work is sure; the qa-check gate and the gate runs are local.

- 2026-11-17: **built** (scratch, nothing committed): the truth tools (`scratchpad/speed2/truth333`: small reference of 17
  panel genes plus a made-up paralog, 65,492 simulated pairs over 81 rows in 8 strata, the oracle, the scorer; 41 tests);
  FORTE's STAR container on HPC (STAR 2.7.9a by its conda record, though the module declares 2.7.10a) aligned the set
  (job 11727992: 98.3% unique, 90,069 splices); the core build `r333` (develop + F2-a + D9 + E0, E1, E2, E6, the opt-in;
  579 Rust and 2,771 Python tests pass, lint clean). **First measurement** on the STAR set: develop's rules and `r333`'s
  rules give identical rows (the default path is untouched); the judge differs on 36 of 81 rows (+30 ALT, +451 REF reads);
  no false ALT at hom-REF rows on any arm; ALT recall against molecule truth: oracle 0.821, rules 0.681, judge 0.682, the
  gap in indels (INS 0.68 against 0.96, DEL 0.39 against 0.53, exon-first-bases deletions 0.013 against 0.21). Running:
  the read-level gap analysis, the switch build (E3 (b), E5, E7, E10), `r333` on the gate's 18 arms against develop.

- 2026-11-17: **the gap, read by read** (Opus; `scratchpad/speed2/truth333/runs/star_obs/analysis.md`; per-read calls
  from the trace, equal to the MAF on 162 of 162 row-arms): the oracle calls 20,974 reads ALT; the rules lose 3,580, the
  judge 3,552; neither path counts a false ALT. Causes (rules / judge): RNA admits no read by its clip (`engine.rs:1967`)
  2,409 / 2,420; RNA clip edges leave the whole clip unread (`carrier.rs:239-252`) 530 / 609 in depth; left-alignment of
  an insertion into the intron in prep (`normalize/engine.rs:1076`, one row) 469 / 472; deletions abutting a junction
  (#385) about 330; the exon-start cut 140 / 37. 2,669 of the rules' losses hold the event in a STAR soft clip. Judge-only:
  a `read_from_right` parity bug (11 reads; `decide.rs:1545-1551`); 424 truth-ALT reads off the anchor called REF, hidden
  today by the depth gate (R013, R044); R044's 157-158 false REF on both paths (STAR's junction 2 bases upstream: the
  invariant-7 limit). The set holds its alleles near junctions by design, so these are causes, not FORTE's rates.

- 2026-11-17: **the switch build** `r333s` (`scratchpad/speed2/r333s`, patch sha256 3d37097e; frozen `plan670/src_r333s`; 591
  Rust and 2,775 Python tests, lint clean; with no switch set, equal to `r333` field for field in 288 of 288 made-up cases;
  every switch on leaves DNA and RNA's default unchanged). **Measured on the STAR set** (ALT recall against molecule truth,
  all rows; REF recall): judge default 0.682, 0.946; E3 (b) overhangs 0.682, 0.946 (INS 0.690, DEL 0.390, exon-start
  deletions 0.021 against 0.013; SNV 0.863 against 0.865, novel acceptor 0.906 against 0.909); E5 one-base SNVs 0.682, 0.946;
  every class to the rules 0.681, 0.931 (the rules path); E7 `NH` doubt 0.682 (no read fits nothing here); oracle 0.821,
  0.984. No arm has a false ALT. The switches move ALT recall by 0.005 at most; the clip rules (RJ-18) hold the gap. The
  judge's REF recall beats the rules' (0.946 against 0.931).

- 2026-11-18: **second review** (Fable; text in `scratchpad/speed2/review333/review.md`; probes in `review333/probe/`):
  revise before any commit. Corrections to this section: "no arm has a false ALT" holds only on the STAR set, which has no
  row of the class the review found: **under the judge's default, a REF read whose 1-4 base overhang the aligner wrote into
  the intron is counted ALT at an insertion row whose inserted bases equal the next exon's bases** (probe: judge AD 12 where
  the rules give 0; E3 (b) gives 0 and recovers the REF reads). (The review's claim that the judge adds 15 counted false REF
  reads was checked by a record-level join and is wrong: secondary alignments inherited their primary's call; the
  judge adds 2, with 164 counted false REF in all against the rules' 168; `truth333/runs/star_obs/rejoin.tsv`.) Also found:
  **an indel whose repeat run crosses an exon edge counts 0 REF and 0 ALT on both paths, status PASS, no reason**
  (rules path included: a defect in today's RNA counts; R041 is one instance, the class is wider); ASJD under the judge
  takes off-anchor REF calls the main counts drop (R013: 667 against 481); the `read_from_right` gap is a class; the
  count-level scorer cannot see false ALT or REF; the set cannot decide E0's respelling, E4, E5 (b), E7 or E11; the tests
  named brute force share the judge's scorer and skip junctions that meet the event; `GBCMS_RNA_SHAPE_TESTS` takes any
  value but `off` as on and no output records the switches; the builds' comments and logs carry ticket labels; the docs
  are not updated. Default paths untouched and binning invariance held under every switch (the review's probes).

- 2026-11-18: **the DNA identity gate passed**: `r333` (develop + F2-a + D9 + the core) on all 18 gate arms against develop,
  0 differences in every field and every whole log (24 comparisons): DNA under both judges and RNA's default path, D9
  included, are unchanged. The truth tools' v2 set (162 rows, 135,389 pairs; the review's missing strata; read-level
  scoring of false ALT and false REF for every arm; a record-level join) is staged for FORTE's STAR on HPC.

- 2026-11-18: **fix-now round** (the review's "(a)" items; no decided rule changed). Code `r333f` (patch sha256 31cb64c2;
  frozen `plan670/src_r333f`; 593 Rust and 2,783 Python tests, 19 strict xfails, lint clean; `gbcms dna` and RNA's default
  equal develop field for field): the judge's RNA admission now admits a read from the right flank when its own aligned
  bases hold the anchor, decided on its frame (11 truth-ALT reads recovered on the STAR set, no other of 52,780 calls
  changes); the shape switch parsed strictly and logged; labels removed from code, logs and help; docs and CHANGELOG
  updated; an independent-scorer fuzz (3,332 reads: entering junctions, clips, qualities around `--min-baseq`, both
  mates); findings 1-3 pinned as strict xfails; a once-per-row log note under the opt-in for rows whose spliced reads all
  fall to depth only across an exon edge. Open from the round: with the shape tests off, finding 2's class also gives
  false ALT (REF reads spliced at the canonical junction, 1-6 reads per seed); finding 1 reaches the rules path too (AD
  12 on one donor-side case); in DNA a clipped off-anchor read whose aligned bases hold the anchor is still held to its
  fragment (a DNA count change, the operator's); the R015 class (reads starting in a deletion's bases). Truth tools v2:
  162 rows, 135,389 pairs, read-level false ALT and false REF for every arm, a record-level join (it confirmed the gap
  analysis's false-REF numbers and refuted the review's). The v2 set aligns on HPC (job 11911044).

- 2026-11-18: **v2 measured, read by read** (`truth333/runs_v2/star_r333f`; build `r333f`; FORTE's STAR on 135,389 pairs;
  103,085 calls resolved, every arm's counted calls equal its MAF on 162 of 162 rows). All rows: ALT recall rules 0.592,
  judge 0.607, E3 (b) 0.608, E5 one-base 0.607, E7 `NH` 0.607, every class to the rules 0.592, all switches 0.608, oracle
  0.796; REF recall 0.820, 0.838, 0.837, 0.838, 0.838, 0.820, 0.837, 0.977; counted false ALT 203, 202 on every judge arm
  (0 from REF reads whose own bases tell the alleles apart; rules 1); counted false REF 514 against 496-500. By stratum
  (rules, judge, oracle): edge insertions 0.063, 0.145, 0.634 with 265 and 263 false REF (truth-ALT reads counted REF:
  STAR aligns the inserted bases as the next exon's); indel runs across exon edges 0.082 on every arm against 0.482
  (finding 2); exon-start deletions 0.013, 0.013, 0.176; homopolymer slippage 0.875, 0.888, 0.967 with 106 false ALT on
  every arm, all slips; the made-up processed paralog 0.399 with 88 foreign false ALT on every arm, rules included;
  low VAF 0.845, 0.847, 0.975; pairs in reach 0.909, 0.910, 0.958. Finding 1's false ALT did not appear with STAR at
  FORTE's settings (annotated junctions spliced from 1 base); it stands for other aligners and novel junctions.

- **Decided (operator, 2026-11-18: "Yes, both", to "Take the ready decisions and approve prototyping the six rules?"):**
  the RNA judge ships in 6.8.0 as an opt-in only, RNA's default stays the rules path; E5 (a) SNV rows read whole; E7 (b) no
  doubt signal for now; E3 (b) kept as a switch, off, decided with FORTE data; E10 every class to the judge under the
  opt-in, the default's routes set per class at the freeze. **Approved to prototype as switches and measure:** admitting
  RNA reads by their clips (RJ-18), after the judge's off-anchor REF calls on truth-ALT reads are fixed; reading clip bases
  up to the exon edge; indel runs across an exon edge (three candidates side by side: an ALT state spliced at the read's
  own junction, the event's placements stopped at the exon edge, no left-alignment into the intron in RNA prep);
  RJ-19's mirror at a junction the aligner placed off; ASJD's admission under the judge; a doubt signal for foreign
  reads (a processed paralog the reference lacks).

- 2026-11-18: the edge-run defect (both paths) posted as #445 after the operator's "Send it"; body read back equal to the draft.

- 2026-11-18: **local measurements** (Opus; `scratchpad/speed2/views333`, `cram333`, `forte_truth`). Per-transcript counts
  take exactly the main counts' reads on both paths (rebuilt 1,403 of 1,403 entries); ASJD takes reads off the anchor on
  both paths (under the judge +221 REF fragments, 146 of them false REF; R013 645 against 482). CRAM gives every count
  field equal to BAM on both paths and on develop. Amplicon mode runs on both paths; one loss under the judge (R160: 2
  truth-REF records to neither). FORTE truth (33 runs, 94 rows, 68 with depth; counts only, through the sanitizer): both
  paths find ALT on the same 27 rows; the judge adds 56 ALT reads on deletion rows, all fitting ALT better by an
  independent base check. **New class, both paths (invariant 7):** a read whose own bases carry a different, unlisted
  allele is credited to the row's ALT when REF fits worse and no candidate holds the true allele (below the found-allele
  floor): FORTE, judge 1 row (2 reads) at a constructed other-allele DELINS; simulated v2 with the same construction,
  reads credited ALT whose bases favour the true allele: rules 80, judge 58. Also: on FORTE SNV rows the judge drops REF
  calls on reads with the REF base at Q20 or more (-18 at positives; 33 of them within 2 bases of a read end or junction).

- 2026-11-18: **prototype build 1** `r333g1` (patch sha256 35de7e4f; frozen `plan670/src_r333g1`; 600 Rust and 2,784
  Python tests, 19 strict xfails; with no switch set equal to `r333f`, DNA and RNA's default equal to develop). On v2, read
  by read: `GBCMS_RNA_CLIP_ADMIT=on` ALT recall 0.607 to 0.625 (+656 reads), false ALT 202 to 205 (2 paralog, 1 slip; 0
  from REF reads whose bases tell the alleles apart), false REF unchanged; it carries the fix for the judge's off-anchor
  REF calls on truth-ALT reads (393 to 8; costs 41 true REF votes at deletion rows an exon edge cuts).
  `GBCMS_RNA_CLIP_READ=edge` 0.614 (+253; on the rules path the part up to the edge applies too: +294). `GBCMS_RNA_ASJD_ADMIT=main`
  ASJD REF total at R013 645 to 482 (the rules path's), main counts untouched. All three: 0.673 (+2,491), false ALT 205 (0
  telling), false REF 496. Exon-start deletions and indel runs across edges barely move (0.013; 0.082 to 0.099).

- 2026-11-18: **prototype build 2** `r333g2` (patch sha256 fb51e618; frozen `plan670/src_r333g2`; 607 Rust and 2,783 Python
  tests, 19 strict xfails; with no switch set equal to `r333f`, DNA and RNA's default equal to develop; `prepare_variants`
  gains a keyword-only `gtf_path`, the stub updated). On v2, read by read (ALT recall / false ALT (telling) / false REF /
  REF recall): judge default .607 / 202 (0) / 496 / .838; `GBCMS_RNA_MOVED_JUNCTION=mirror` .607 / 202 / 68 / .838 (428
  truth-ALT reads from REF to depth, no true call lost); `GBCMS_RNA_FOREIGN_DOUBT=on` .607 / 138 / 496 / .838 (paralog 88
  to 24; costs 16 of 37,653 truth-ALT and 15 of 71,843 truth-REF reads); `GBCMS_RNA_EDGE_RUN=stop_at_edge` .664 / 202 / 636
  / .906 (edge runs .082 to .402, oracle .482; `alt_state` the same with 54 more false REF; `no_left_align` changes no call);
  **stop_at_edge + mirror + foreign: .663 / 138 (0) / 77 / .906**; on the rules path the same three: .647 / 139 (1) / 104 /
  .888 against .592 / 203 / 514 / .820. All apply on both paths. Open: alt_state or stop_at_edge (they differ only at
  novel junctions); edge-run deletions written inside their junction stay depth only (crediting them would choose a
  splice); the mirror's 30-base reach and the foreign test's threshold need real data.

- 2026-11-18: **the six rules together**, build `r333g12` (`r333g1` and `r333g2` merged over `r333f` by `git merge-file`,
  six conflicts resolved by keeping both sides; patch sha256 cfe75364; frozen `plan670/src_r333g12`; clippy clean, 614 Rust
  tests, ruff, black, mypy, 2,784 Python tests and 19 strict xfails pass). Merge check: every arm that sets one build's
  switches gives that build's own per-read calls and MAF (the calls' only difference is the trace's `ns=` timing field);
  DNA (both judges) and RNA's default equal develop on 192 of 192 made-up checks, DNA also with all seven switches set;
  the opt-in with no switch set equals `r333f`; a misspelt switch stops the call. On v2, read by read (ALT recall / false
  ALT (telling) / false REF / REF recall): rules .592 / 203 (1) / 514 / .820; judge .607 / 202 (0) / 496 / .838; **all six
  under the judge (clip admit, clip read, ASJD main, stop_at_edge, mirror, foreign) .730 / 141 (0) / 77 / .906**; with
  alt_state .730 / 141 / 80 / .906; with E3's overhang placement too .727 / 142 (1) / 86 / .904 (E3 adds nothing on
  STAR); all six on the rules path .655 / 140 (1) / 104 / .888; oracle .796 / 333 / 1,287 / .977. The two builds' gains
  add (.607 + .066 + .056 = .729). By stratum (judge default, six, oracle): edge insertions .145, .530, .634 (false REF
  263, 8, 490); edge runs .082, .419, .482; exon-start deletions .013, .013, .176 (untouched: the #385 class); junction
  .675, .760, .833; low VAF .847, .952, .975; paralog false ALT 88, 26, 90. Tables: `scratchpad/speed2/r333g12/tables.md`.
- 2026-11-18: **third adversarial review** (Fable, build `r333g12`; `scratchpad/speed2/review333b/review.md`; every
  number rerun by the parent, `verify/`). Sound on every probe: clip read, ASJD main, the mirror. Defects: (1) clip admit
  admits reads whose clip fits nothing (400 junk-clipped pairs over an SNV row with 6 REF reads: RD 6 / AD 0 becomes
  12 / 4); develop's DNA judge does the same, which is the "about 1% wrong reads" the operator accepted on 2026-10-24 to
  26 when DNA kept `admit` over `clean:0` (8% of true ALT); for RNA the rule is new, so its own cost is measured; (2)
  the foreign-read doubt drops the reads of the haplotype that carries an unlisted germline change in cis (AD 18 to 9
  with two SNPs, to 4 with a DNP); (3) it reads the call's other rows, so one row per call differs (invariant 1); (4)
  alt_state turns ALT reads written with the deleted base inside the N into REF (19 of 19); stop_at_edge counts their
  REF reads and none of their ALT (v2's edge-run deletions VAF 0.06 against 0.41); (5) stop_at_edge and the mirror do
  nothing without `--gtf`; (6) the judge core drops a call at an SNV row one base from a wrong last base (a one-base gap
  −7 beats two mismatches −8; DNA too); (7) the unlisted-allele class is unchanged; (8) the overhang false ALT stands
  (only `GBCMS_RNA_OVERHANG=spliced` clears it); (9) all six under the judge take 11.4 s against the judge's 6.1 s and
  develop's 2.1 s (the mirror +2.7 s, the foreign test +0.8 s); (10) 7 of 24 mutants pass every test; (11) no_left_align
  changes nothing. Pre-existing on develop's RNA rules path: clips holding another exon's bases give 2 false ALT at a
  4-base deletion. Per rule on v2 (judge; ALT recall / false ALT / false REF): default .607 / 202 / 496; stop_at_edge
  .664 / 202 / 636; mirror .607 / 202 / 68; foreign .607 / 138 / 496.
- 2026-11-18: **fix round on `r333g12`**, three Opus builds side by side, each its own files: `r333g3a` (alt_state and
  no_left_align removed; edge-run candidates `stop_balanced` and `own_junction` beside stop_at_edge; RNA clip admit
  `clean:K` on every row class, DNA unchanged; a GTF-needing switch without `--gtf` stops the call; tests for the
  surviving mutants and the CLI); `r333g3b` (the foreign test without the call's other rows; candidates `min:N`,
  `share:F`; the mirror's speed by exact speed-ups only); `r333g3c` (the scorer's read-end gap, `GBCMS_JUDGE_END_GAP`,
  DNA and RNA, measured on v2 and HG002).
- 2026-11-18: **M14 and M15 on HPC** (Opus agent, stopped at handback; the main session takes over). Slip: the first HPC
  build of `r333f` reused develop's Rust extension (both sha256 f763393a): `tar` restored file times older than the shared
  cargo target's outputs, so cargo compiled nothing, and every run on it was develop or failed; all its jobs were cancelled
  by name. Rebuilt as `r333fb` (extension 669e6237; its source unpacked locally equals `frozen_r333f` file for file).
  Source fix: `hpc/kit/build.sbatch` touches every unpacked source file (HPC copy sha256 9e8b2643). Second slip: the
  agent's ledger times for 12:10, 12:25 and 12:40 do not match the clock (by sacct the runs went in at 15:20-15:24);
  corrected in the ledger. Running at 15:52: M14 on the 58 FORTE slices (39 of 58 done), M14 on the 40 FORTE BAMs with
  #292's matched IMPACT DNA counts (job reads the BAMs in place, counts only), M15 timing (pf76's cells; the small cell
  done). First M15 repeat, 1,000 rows x 1 BAM: C++ 38.4 CPU s / 11.7 s wall / 3,816 MB; develop rules 9.5 / 5.1 / 642;
  `r333fb` rules 8.0 / 4.4 / 646; judge 20.9 / 7.6 / 728 (one repeat, not a median).
- 2026-11-18: **M15 medians** (3 repeats, Gold 6348; CPU s / wall s / peak MB): 1,000 rows x 1 FORTE BAM: C++ 37.6 /
  10.7 / 3,816; develop rules 8.3 / 4.4 / 642; `r333fb` rules 8.0 / 4.4 / 630; judge 21.1 / 7.7 / 725 (0.56x the C++
  tool's CPU). **5,000 rows x 10 BAMs: C++ 671.9 / 200.7 / 23,981; develop rules 218.2 / 66.8 / 3,286; `r333fb` rules
  219.6 / 67.5 / 3,235; judge 808.9 / 223.2 / 3,472: the RNA judge with no switch is 1.20x the C++ tool's CPU and 1.11x
  its wall, below the speed floor** (3.71x develop's CPU; memory 1.06x). `hpc/out/m15b_time.tsv`.
- 2026-11-18: **fix build B** `r333g3b` (Opus; patch sha256 01238a81; frozen `plan670/src_r333g3b`; 618 Rust and 2,812
  Python tests, 19 strict xfails, clippy, ruff, black, mypy clean; defaults equal `r333g12` on 20 of 20 set comparisons
  and 6,140 made-up ones). The foreign test reads the row and the read alone (the call's other rows dropped); values
  `on` (= `min:2`), `min:N`, `share:F`, `dense:K`. Haplotype probe (2 SNPs, a DNP, 3 SNPs in cis): every `share` value
  keeps AD 18 / RD 17; `on` 9 / 8 / 4 / 3. Per call (parent rerun): `on` and `share:0.1,dense:4` give one call = one row
  per call; `on` now counts listed rows in cis as differences (AD 17 to 10, 10 to 0). v2 (judge; ALT recall / false ALT
  / false REF; paralog false ALT): default .607 / 202 / 496, 88; `on` .607 / 133 / 496, 19 (cost 24 truth-ALT, 81
  truth-REF reads); `min:3` 68; `share:0.2` 36 (cost 13 / 15); `share:0.1,dense:4` 19 (cost 13 / 15); five plus
  `share:0.1,dense:4` .730 / 136 / 77, 21. Limits: a sparse copy (1 difference in 60 bases) at a share of 0.33 or more is
  missed by every value; a haplotype below the share loses its ALT reads (allelic imbalance, clustered somatic changes);
  4+ real SNPs in 100 bases (HLA-like) lose the exemption (not measured). **The mirror's memo** (its verdict kept per row
  in the row's exon-edge object, keyed on every input it reads): MAFs byte-equal to `r333g12` (parent rerun, the five
  switches under both judges); judge plus mirror 8.8 to 6.3-6.9 s, rules plus mirror 5.0 to 2.4-2.6 s.
- 2026-11-18: **fix build A** `r333g3a` (Opus; patch sha256 a361e776; frozen `plan670/src_r333g3a`; 621 Rust and 2,789
  Python tests, 19 strict xfails, gate clean; defaults equal `r333g12` on 2,068 of 2,068 MAFs; `normalize/engine.rs`,
  `_rs.pyi` and `pipeline.py` back to develop's text, checked). alt_state and no_left_align removed (they now stop the
  call). A GTF-reading switch without `--gtf` stops the call. Edge-run deletions, three values beside stop_at_edge:
  `stop_balanced` (spliced reads at the edge depth only), `own_junction` (a read spliced at the annotated junction made
  the deletion's length longer, or a slide of it, is ALT where no splice motif on the GTF's strand lets REF splice there
  and the read fits ALT by a confident base on the annotated frame), `motif` (`own_junction`, balanced to depth only
  where the shorter intron has a motif too: NAGNAG). Edge-run probe (parent rerun): `own_junction` and `motif` REF-only
  19/0, ALT-only 0/19 in every shape. v2 (ALT recall / false ALT / false REF / REF recall): six with stop_at_edge .730 /
  141 / 77 / .906; **six with `own_junction` or `motif` .766 / 141 / 77 / .906** (edge-run deletions .092 to .472,
  exon-start deletions .013 to .450, the #385 class); `stop_balanced` .728 / 141 / 75 / .839 (exon-start REF recall .885
  to .095). Cost: +3.6 s on v2 (two full-programme scores per read that carries the longer intron). RNA clip admit
  `clean:K` (every row class, RNA only): junk-clip probe RD 6 / AD 0 for every K; six with `clean:0` .728, `clean:1` .722
  (low VAF .952 to .929), `clean:2` .716. Open: edge-run insertions; R013 (STAR clips its ALT reads).
- 2026-11-18: **fix build C** `r333g3c` (Opus; patch sha256 291c9f7e; frozen `plan670/src_r333g3c`; 624 Rust and 2,788
  Python tests, 19 strict xfails, gate clean; with the switch unset equal to `r333g12` on v2, the dense set and 100 made-up
  cases; the rules paths equal with it set). The scorer's end bar `GBCMS_JUDGE_END_GAP=barK`: a gap needs more than K
  confident read bases on each side to the read's ends (masked bases do not count), in every programme of the scorer,
  under `--read-judge sequence` only, DNA and RNA. Read-end drop probe (parent rerun): bar4 brings back 8 of 8 REF and 8 of
  8 ALT in every case on the judge and the six; the rules path unchanged. v2 (read level): judge off .6072 / 202 / 496,
  bar4 .6075 / 202 / 496 (14 right calls gained, 2 lost); bar2 9 lost, 3 gained. HG002 (public; tandem repeats, adjacent
  pairs, the SV subset; `--truth v2dbq20`): no class fails; bar4 gains 20 right calls and loses 13, bar2 gains 13 and loses
  4. Cost with the bar on: +5-8% wall. Open: a junk tail longer than K past the wrong base keeps the gap; FORTE not
  measured (the 33 REF drops within 2 bases of a read end or junction).
- 2026-11-18: **merged fix build** `r333g3` (`r333g3a` and `r333g3b` merged over `r333g12` by `git merge-file`: two
  conflicts in engine.rs resolved by keeping A's GTF-edge plumbing and B's share tally, A's listed-rows attachment dropped
  as B removed it; one prototype test line changed from `foreign_doubt: true` to `ForeignDoubt::ON` (a type change);
  `r333g3c` merged on top with no conflict). Patch sha256 f5c0577a, 23 files against develop; frozen
  `plan670/src_r333g3` (0 build warnings; its source equals the gated tree file for file). Gate: clippy clean, 635 Rust
  tests (625 + C's 10), ruff, black, mypy clean, 2,821 Python tests (2,784 + 5 + 28 + 4), 19 strict xfails. M16 (HPC
  agent, Opus): `r333g3` on FORTE (the 58 slices, the 40 BAMs with matched DNA) with the candidate set and timing.
- 2026-11-18: **`r333g3` on v2** (defaults: DNA and RNA under both judges equal `r333g12` on v2 and the dense set, 8 of 8).
  The candidate set "best" = clip admit `clean:0`, clip read `edge`, ASJD `main`, edge run `motif`, mirror, foreign
  `share:0.1,dense:4`. Read level (ALT recall / false ALT / false REF / REF recall): rules .592 / 203 / 514 / .820; judge
  .607 / 202 / 496 / .838; `r333g12`'s six .730 / 141 / 77 / .906; **best .765 / 136 / 77 / .906**; best with the end bar
  .766 / 136 / 77 / .906; with the foreign test off .766 / 205; with `clean:1` .759 / 135; with clip admit off .707 / 133;
  with clip admit `on` .767 / 136; with stop_at_edge .728; with own_junction .765 (equal to motif here); with no edge rule
  .672 / 136 / 68 / .837; the best switches on the rules path .692 / 135 / 104 / .888. By stratum (judge, best): edge runs
  .082, .590; exon-start deletions .013, .450; edge insertions .145, .528 (false REF 263, 8); junction .675, .789; low VAF
  .847, .952 (`clean:1` .929); pairs .910, .949; paralog false ALT 88, 21 (foreign off 90). Time on v2 (1 thread, median
  of 3, wall / peak MB): develop rules 1.86 / 65; rules 1.97 / 66; judge 5.66 / 98; `r333g12` six 10.74 / 104; best 11.49
  / 102; best with the bar 11.86 / 104 (the mirror's memo saves about 2.5 s, the motif rule costs about 3.6 s: two
  full-programme scores per read that carries the longer intron, scored again by the per-transcript and ASJD callers).
  Tables: `scratchpad/speed2/r333g3/tables.md`.
- 2026-11-18 21:55: **machine restart at 22:00; where to resume.** The session scratchpad (`/private/tmp/...`, cleared by
  a restart) is copied, without build caches, venvs, source copies or made-up case BAMs, to
  `~/test/gbcms/harness/scratch_kept/session_e9e65c2a_20261118/` (7.1 GB): the reviews (`speed2/review333/review.md`,
  `speed2/review333b/review.md` and probes), every patch (`speed2/r333g3/r333g3.patch` f5c0577a and the fix builds'),
  the truth set and scripts (`speed2/truth333/`, incl. `star_v2_after_r333g3.sh`), tables, timing. The frozen builds
  stay in `plan670/src_*` and `frozen_src/`. **M14 done** (copied back 19:56: `hpc/out/m14bagg/`, summary
  `hpc/out/m14_summary.md`, read it through `plan670/judge680/rbr/san`). **M16 running on HPC** (the HPC agent dies with
  the restart; the jobs do not): `auto-gbcms-m16s` 12137764, 12137841; `auto-gbcms-m16d` 12137868, 12137913; timing
  `auto-gbcms-m16t.s` 12137928, `.m` 12137948; aggregate `auto-gbcms-m16agg` 12144643 (afterany; script
  `kit/rna333g3_agg.py`). Next: sacct those ids (a ledger line each step), copy back the m16agg tables and the m16t
  times (times and node only), print through the sanitizer, then bring the decisions (the six rules with the fixes, the
  end bar, the foreign test, RNA judge speed: 1.20x the C++ tool's CPU on 5,000 rows x 10 BAMs) with a visual.
- 2026-11-19: **M16 on FORTE** (all jobs exit 0; aggregate verify changed=0, rows withheld=0; `hpc/out/m16agg/`, read
  through the sanitizer). Arms: develop's rules (dvR), the judge (fJ, M14), `r333g3`'s judge (gJ), the candidate set gB
  (clip admit `clean:0`, clip read, ASJD main, motif, mirror, foreign `share:0.1,dense:4`), gBb (gB + end bar), gBf (gB,
  foreign off), gB1 (gB, `clean:1`). **Identity:** gJ equals fJ on every line (DNA-check SNV rows on 40 BAMs, the 58
  slices' SNV and indel rows: 0 differ). **False ALT on real data, which v2 never showed:** DNA-check rows where the
  matched normal is hom-REF and the oracle sees no ALT (36,613 rows, about 15.7 M REF reads): ALT reads / rows with ALT
  / rows with 2+: dvR 0 / 0 / 0; fJ 1 / 1 / 0; **gB 531 / 205 / 108; gBb 670 / 244 / 127; gBf 533 / 207 / 108; gB1 120
  / 49 / 22**. The slices' oracle-clean SNV rows (173,508): dvR 0; fJ 22 on 10 rows; gB 2,203 on 905; gBb 2,935 on
  1,080; gB1 368 on 202. So clip admit `clean:0` makes most of it (gB to gB1: −411 reads), the end bar adds 140, the
  foreign test none. **Gains:** germline hets (559 rows, DNA VAF 0.2-0.8, RNA depth 10+): fJ to gB +4,095 ALT reads on
  300 rows (gB1 +1,909); mean RNA-DNA VAF −0.0247 (fJ) to −0.0155 (gB), −0.0201 (gB1); rows more than 0.2 off 73, 67,
  69. Somatic (3 rows): ALT 842 to 850. Indels on the slices (fJ to gB): edge-run deletions +172 ALT (clean rows +19 on
  6), deletions flush with an exon edge +1,068 ALT (het rows +512 on 15; clean +18 on 8) with 17,614 REF reads to depth
  where motif balances an ambiguous splice; other deletions +1,972 (clean +42 on 20). **REF:** the judge alone drops
  6,439 REF reads on 1,837 clean hom-REF rows (the read-end class); the foreign test withholds 14,924 more REF reads on
  4,445 of them; the end bar brings back 2,984 and drops 1,300. **CPU** (summed; DNA-check / slices SNV / slices indel):
  dvR 3,812 / 2,992 / 5,344 s; fJ 32,757 / 56,459 / 15,805; gB 75,406 / 91,387 / 22,224; gBf 41,792 / 67,954 / 18,748;
  gBb 77,046 / 94,298 / 23,061 (the foreign test's share tally is the largest cost on real BAMs). Open before any
  choice: which switch makes the false ALT left under `clean:1` (clip read, the edge rule, clip admit itself), clip
  admit `clean:2` and `clean:3`, and the end bar alone, on real data: M17.
- 2026-11-19: **M16 timing** (median of 3, Gold 6348, beside M15; CPU s / wall s / CPU against the C++ tool): 1,000 rows
  x 1 BAM: judge 21.9 / 8.0 / 0.58; candidate set 35.1 / 11.3 / 0.93; with the end bar 35.7 / 11.4 / 0.95. **5,000 rows x
  10 BAMs: judge 834.1 / 226.0 / 1.24; candidate set 1,494.0 / 403.6 / 2.22; with the end bar 1,510.7 / 408.0 / 2.25**
  (C++ 671.9 / 200.7; develop rules 218.2 / 66.8). The candidate set is twice the C++ tool's time on the large cell.
- 2026-11-19: **M17 submitted** (`kit/rna333g4_*`, copies of M16's scripts with six arms more, each from the candidate
  set with the foreign test off, one switch changed: the end bar alone on the judge; clip admit off; clip read unread;
  `clean:2`; `clean:3`; the edge rule off). Jobs 12220952, 12220999 (`auto-gbcms-m17d`, the 40 BAMs), 12221011, 12221013,
  12221104, 12221170 (`auto-gbcms-m17s`, the 58 slices), aggregate 12221195 (`auto-gbcms-m17agg`, afterany, to
  `runs/m17agg`). Next: sacct, copy back `runs/m17agg`, print through the sanitizer; then the decisions with a visual.
- 2026-11-19: **fourth adversarial review started** (Fable; the measurement plan and its real-data reading: truth labels,
  row set, power, M17's design, the aggregator, unchecked invariants on real data, speed, what each rule's evidence
  supports; probes in `scratchpad/review333c/`). **M18 approved** (operator, 08:57: "Yes, counts only", to reading the 40
  matched IMPACT tumour BAMs in place, counts only, at panel sites): the slices' SNV and indel rows inside the panel,
  counted in each matched tumour (`gbcms dna`, develop), to label RNA rows tumour-ALT (sensitivity, by DNA VAF band) or
  tumour hom-REF (false ALT, the strict truth: the RNA is the tumour's), and to recheck M16's DNA-check "absent" rows
  against the tumours (HPC agent, Opus). Normals are not reread.
- 2026-11-19: **fourth adversarial review** (Fable; `harness/scratch_kept/session_e9e65c2a_20261118/review333c/review.md`;
  the parent checked the two code claims and reran the two probes). The M16 numbers are right; their reading is weaker:
  (1) the clean rows' false-ALT baseline is zero by selection (no aligned Q20+ mismatch, ALT a fixed transversion); on
  all DNA hom-REF rows (245,492, 510 M REF reads) the switches add +4.35% ALT (dvR 906,357 to gBf 935,567), almost the
  share they add at germline hets (+3.89%): clip admission is a flat "about 4% more ALT wherever ALT is seen", strongest
  on A>G (antisense +16%); (2) the #292 oracle reads aligned bases only (`edge_sites.py:107`), so "oracle clean" cannot
  judge a clip rule; (3) `clean:K` only cuts the read's last K bases (a local aligner clips a last-base mismatch), and a
  1-2 base overhang from another exon placed on the row is credited ALT under `clean:0` and `clean:1`, not under
  `clean:2` (invariant 7: placement, not bases, for K below 3); (4) the edge rule's ALT has no DNA truth, and it raises
  REF reads at edge-run deletion rows by 42%; (5) rows recur across samples and sit 0-9 bases from junctions; 3 somatic
  rows give no power; (6) the foreign test removes ALT about 10x more where DNA lacks it than at hets, but its paralog
  class label is circular and it is 70-82% of the added CPU; (7) the end bar adds false ALT on RNA; (8) the judge alone
  costs 19x the rules on SNV slice rows with 1.5% fewer ALT and no better het VAF: the real SNV gain is clip admission;
  (9) the aggregator drops a task for every arm when one arm is missing and reads no ASJD or per-transcript column.
  Decidable today: the end bar off for RNA; clip admit not at `clean:0`. Proposed: M-A (the patients' own IMPACT somatic
  calls as rows: needs the operator), M-B (a fair hom-REF row set), M-C (clip tallies in the oracle), M-D (DNA at the
  edge-rule rows), M-E (aggregator columns), M-F (mirror off, a binning arm, a foreign-off timing arm).
- 2026-11-19: **the operator chose "All, with M-A"**, then for M-A's rows **"Derive from tumour DNA"** (no mutation
  lists; cBioPortal files never go to HPC): a scan of each matched tumour BAM inside the panel (ALT 3+, VAF 0.02+, depth
  30+), confirmed by `gbcms dna`, lifted to GRCh38 and counted in the patient's FORTE BAM; the union of all patients'
  sites, hom-REF in another patient's tumour, as negative controls; by tumour VAF band (0.35+, 0.1-0.35, 0.02-0.1).
  **Slip found and stopped:** the repo's HPC rules say an array runs at most 20 tasks at once and only the main session
  runs HPC commands; M14, M16 and M17 ran 40-60 at once (M17's by the main session), and the M14, M16 and M18 agents ran
  HPC commands. From now on `%20` and the main session runs every step; the M18 agent hands back its scripts; a script
  agent writes M-A, M-B, M-C, M-E and M-F locally. Work-order `.agents/learnings/LEARNINGS.md` LRN-20261119-001 (the
  guard hook to enforce both; the operator decides).
- 2026-11-19: **M18 submitted by the main session** (scripts by a subagent, reviewed: reads only the pairs list's
  tumour BAMs, in place; prints through the sanitizer). `kit/dna333_tumour.sbatch` job 12243296 (40 indices, `%20`);
  aggregate `kit/rna333_m18_agg.sbatch` job 12243450 (afterany on M18 and M17's aggregate 12221195; to `runs/m18agg`).
  It covers M-D: the slices' indel rows inside the panel (coding exons plus 20 bases; not UTR edges) carry M16's edge
  classes. Built-in check: M16's clean absent rows under gB must give back 531 / 205 / 108. Next: copy back
  `runs/m18agg`, print with `kit/local/sum_m18.py` through the sanitizer.
- 2026-11-19: **M-A, M-B, M-C, M-E, M-F submitted by the main session** (scripts by a subagent, written locally with
  dry runs; reviewed: every job's output goes through phi_sanitize, the tumour jobs read only the pairs list's tumour
  lines, none deletes or copies; four scripts and the aggregate job changed by the main session to take absolute
  paths under the work folder, which the HPC guard asks for). M18's 40 tasks finished exit 0 (its aggregate waits on
  M17's). Jobs (every array `%20` at most): `mEagg16` 12250628 (the new aggregator on M16 and M17); `mFs` 12250634,
  12250655 (slices: gR identity, rS the rules path with clip read, motif and mirror, gB1f, gBm0 mirror off); `mFb`
  12250675 (one row per bin and one thread, 5 slices); M-A `mAs` 12250914 (tumour scan inside the panel), `mAu`
  12250946 (union, lift), `mAt` 12250963 (tumour counts at the union), `mAr` 12251108, 12251117 (RNA arms there);
  M-B (12,000 random panel rows, 3,000 per distance stratum, made locally from public references) `mBt` 12250951
  (tumour counts), `mBr` 12251004, 12251007 (RNA arms); `mFd` 12251050, 12251055 (M-F arms at the DNA-check rows);
  M-C `mCf` 12251144, 12251150 (clip tallies, full BAMs), `mCs` 12251180 (slices); timing `mFt.s` 12251192, `mFt.m`
  12251208 (foreign-off arm beside gJ and gB); final aggregate `mEagg` 12251371 (to `runs/mEagg`). The full-BAM jobs
  run one array after another. Next: copy back `runs/mEagg16`, `runs/m17agg`, `runs/m18agg`, `runs/mEagg` and the
  `mFt*` times; print through the sanitizer.
- 2026-11-19: **results of M17, M18, M-A, M-B, M-C, M-E and M-F** (every job exit 0; aggregates only, printed through
  the sanitizer: `hpc/out/m17agg`, `m18agg`, `mEagg16`, `mEagg`, `mFtG`, `mFtB`, `mFtBf`; local summariser
  `kit/local/sum_mE.py`). **Checks:** under the candidate set, one row per bin, a bin cap and one thread give 0 different
  MAF lines (5 slices, SNV and indel); develop's rules equal `r333g3`'s (1 slice indel line, as in M14); fJ equals gJ;
  the four count relations hold in every arm and part; per-transcript counts never pass the main counts; 3,282 arm runs,
  0 missing. **Row sets:** M-A 6,076 sites from the 40 tumour scans: 31,555 row-samples where the patient's own tumour
  carries the allele (5,284 distinct positions), 163,316 negative controls hom-REF at depth 50+ in every tumour of the
  patient; M-B 12,000 random panel positions, 438,899 hom-REF row-samples; M18 rechecks M16's "absent": 225,066 of
  245,492 rows have no DNA ALT read in any tumour or normal. **Clip admission** (arms from the candidate set with the
  foreign test off; `clean:3` stands for off, as on M18 where off and `clean:3` differ by 40 of 19,722 reads):
  `clean:0` / `clean:1` / `clean:2` against `clean:3`: own-allele ALT +104,842 (+3.4%) / +47,660 (+1.5%) / +870; RNA
  minus DNA VAF −0.0285 to −0.0240 / −0.0265 / −0.0284; controls ALT +5,157 (+11.5%) / +2,057 (+4.6%) / +421; M-B
  false ALT 188 ppm to 243 / 200 / 189, rows ALT 2+ 5,521 to 7,411 / 5,920 / 5,530; M16's clean strict hom-REF rows
  492 / 110 / 18 (`clean:3` 16, off 1). M-C (an oracle that reads soft clips): clipped ALT bases that pass `clean:K`, own
  alleles K0 119,882, K1 55,158, K2 339; M-B 18,712, 2,870, 34; a junction other reads carry at the clip start holds
  1.5% of the own-allele clip ALT and 24% of M-B's. **Clip read** (`edge` against `unread`, under `clean:0`): hets +1,505
  ALT, clean hom-REF rows +54, CPU +12-14%; on the rules path (rS) the SNV rows do not change (M-A own alleles −46 of
  3.09 M ALT). **ASJD:** rows where ASJD's ALT is above the main ALT under today's judge: DNA-check 3,710, slices 3,183,
  M-A 812, M-B 1,508; 0 under `main`; develop's rules keep 213 slice indel rows (rS 231). **Edge rule** (`motif` against
  off): SNV rows unchanged; slice indel rows tumour hom-REF +3,506 REF and 0 ALT; DNA-absent indel rows +1 ALT; no
  tumour-ALT edge-run row in any set. **Mirror** (gB against gBm0): DNA-absent rows −1,544 ALT (0.17%), −35,545 REF to
  depth; hets −5 ALT; CPU within 1.5%. **Foreign test** as M16. **End bar** on the judge alone: clean-row false ALT 1 to
  45, tumour hom-REF +67 ALT. **The judge against the rules on SNV rows** (fJ against dvR): own alleles ALT −0.14%, REF
  −0.56%, 131 fewer rows with ALT; at DNA VAF 0.02-0.1 ALT −19% and 86 fewer rows with ALT; controls ALT −2.8%; M-B
  −2.1%; CPU 6-19x (sums over shared nodes). **Timing** (median of 3, 5,000 rows x 10 BAMs, CPU against the C++ tool):
  judge 1.25x, candidate set without the foreign test 1.36x, with it 2.23x. Slip: the first summary counted M-C's clip
  reads twice (the table holds a per-distance row and a total for each read); ratios were right, counts doubled; fixed in
  `sum_mE.py` (total rows only). Decisions page: `~/.agent/diagrams/gbcms-333-rna-decisions.html`.
- 2026-11-19: **M-G submitted** (the operator: "Why not also measure D2 first and then we decide what to keep i feel we
  should measure if we already built"). One rule out of the recommended set R (`clean:1`, clip read `edge`, ASJD
  `main`, `motif`, the mirror, the foreign test off; arm gB1f) at a time: gB1fr (clip read unread: D2), f0 (clip admit
  off), f0r (both off: clip read alone), gB1fm0 (the mirror off), gB1fe (the edge rule off, slices only). Jobs
  (every array `%20` at most): `mGs` 12328727, 12328766 (slices); `mGb` 12328772, 12328824 (M-B); `mGa` 12328861,
  12328872 (M-A); `mGd` 12328917, 12328942 (DNA-check rows); D2's timing `mGt.s` 12328980, `mGt.m` 12328999 (the judge,
  R, R with clip read unread; both cells, 3 repeats); aggregate `mGagg` 12329014. The full-BAM arrays run one after
  another. Rule promoted: `.agents/rules/decision-protocol.md` §8 coverage check (LRN-20261119-002).
- 2026-11-19: **the build carries five more switches never run on real data** (`GBCMS_RNA_ROUTE`, `_SNV`, `_DOUBT`,
  `_OVERHANG`, `_SHAPE_TESTS`; `ROUTE=snv=rules` is D9's alternative, built in an earlier round), and other built values
  (edge rule `own_junction`, `stop_at_edge`, `stop_balanced`; foreign test `on`, `share:0.2`; end bar `bar2`; the end bar
  in DNA) lack real-data runs. **The operator: "We should do all"** (M-H, the local input-type checks, and (a) the DNA
  judge with and without the end bar on the 40 matched tumour BAMs, counts only, at panel sites; (b) a tumour scan for
  indels near exon edges; (c) UCSC's segmental-duplication track for paralog truth rows, the file confirmed first:
  `genomicSuperDups.txt.gz`, hg38, 4.2 MB, not yet downloaded). **M-H submitted** (one built option added to R at a
  time; 13 arms run first on the simulated v2 set, every one exit 0; one-base SNVs and the SNV route give the same v2
  MAF): slices `mHs` 12332132, 12332149, 12340195, 12340216, 12340461, 12340513; M-B `mHb` 12340549, 12340573,
  12340602, 12340616; M-A `mHa` 12340654, 12340675, 12340718, 12340737; DNA-check `mHd` 12340780, 12340859, 12340919,
  12340973, 12341011, 12341098; timing (the judge, R, R with the SNV route, R with NH doubt) `mHt.s` 12341134, `mHt.m`
  12341173; aggregate `mHagg` 12341533 over every arm so far. Scripts for (a) and (b) are being written by two Opus
  agents (local only, dry-run on simulated trees; the main session submits); the input-type checks run locally
  (Sonnet, simulated data).
- 2026-11-19: **input types under R, on the simulated v2 set** (Sonnet; `scratchpad/inputs333/`; the parent rechecked
  the CRAM and thread identity and the relations): rules, judge and R each give MAFs equal line for line on CRAM and
  BAM and on 1 and 4 threads; the four count relations hold on every row of all 21 runs (3,402 rows); amplicon mode,
  unpaired reads and hard-clipped reads run with exit 0; R never counts a hard-clipped base (with 30% of soft clips made
  hard, R's ALT falls on 68 rows and rises on none). The kept `cram333/syn.v2.cram` holds other reads than the v2 STAR
  BAM; a fresh CRAM was made for this check.
- 2026-11-19: **D7 on real DNA submitted** (approval (a)): `gbcms dna --read-judge sequence` on every matched tumour BAM
  at the M-A and M-B rows (b37), with no switch, with `GBCMS_JUDGE_END_GAP=bar4` and `bar2`, and develop's judge as an
  identity check (scripts by an Opus agent, dry-run on a simulated tree, reviewed: tumour lines only, in place, prints
  through the sanitizer). Jobs `d7b` 12345792, `d7b7` 12345859 (index 7), aggregate `d7bagg` 12345909.
- 2026-11-19: **paralog truth rows submitted** (approval (c); the operator's "Sure" to the download): UCSC hg38
  `genomicSuperDups.txt.gz` (4.2 MB, public) gives 286 panel exon positions where a duplication's other copy differs at
  that base with 15 equal bases each side; ALT = the other copy's base (the base a mis-placed paralog read carries);
  bands by the pair's identity: 0.99+ 24, 0.98-0.99 9, below 253; 279 lift to b37 (`kit/local/rna333h_sd_rows.py`, public
  files only). Tumour counts `mSt` 12353827; RNA arms (rules, judge, R, R with the foreign test at share:0.1/dense:4, on,
  share:0.2, and NH doubt) `mSr` 12353949, 12353993 after M-H's timing; aggregates `mSall`, `mShi`, `mSmid`, `mSlo`
  12354046, 12354061, 12354103, 12354120. Limit: few rows at identity 0.99+, where paralog reads mis-map most.
- 2026-11-19: **M-I submitted** (approval (b); scripts by an Opus agent, dry-run on a simulated tree, reviewed): exon
  edges lifted to b37 from public files (`mIp` 12354413); a scan of every tumour BAM for 1-50 bp indels within 30 bases
  of an exon edge, depth 30+, ALT 3+, VAF 0.02+ (`mIs` 12354436, 12354443); union and lift to GRCh38 (`mIu` 12354458);
  tumour counts (`mIt` 12354475, 12354486); RNA arms rules, judge, R, R with the edge rule off, `own_junction`,
  `stop_at_edge`, `stop_balanced` (`mIr` 12354515, 12354542, after M-H's timing); aggregate by edge class (`mIagg`
  12354656).
- 2026-11-19: **wider row sets** (the operator: "We should try to make our test comprehensive we have ebooks data to do
  it", then "*enough What about the 50 forte bams dont they have enough data? They might not be with impact but on their own").
  Slip: I read the typo "ebooks" as a data-set name and offered only new data sets (LRN-20261119-003). The RNA is whole
  transcriptome; only the DNA labels stop at the panel. Sizes from public files (`kit/local/hg002_cover.py`, Ensembl
  111 exons): exonic paralog rows inside HG002's v5.0q benchmark 124,274 (0.99+ 4,506; 0.98-0.99 12,384; below
  107,384) against the panel's 286 (24); exonic bases inside the benchmark 148.4 M (0.99+ segdups 1.25 M) against the
  panel's 1.0 M (2,129). HG002 v5.0q exonic truth (one person; genome / panel): het SNV 96,127 / 277; indels 29,204 /
  19; indels within 30 bases of an exon edge 1,642 / 9, of which 146 have a repeat run that reaches the edge (0 on the
  panel); 50+ bp events 844 / 0; cis MNP candidates 1,476 / 1; cis delins candidates 3,729 / 3. Public HG002 RNA
  (Fable survey): GIAB Google Illumina mRNA (poly-A, stranded, 2x150, about 107 M pairs; FASTQ 14.7 GB, a STAR 2.6.1d
  GRCh38 BAM 9.8 GB); UNC replicates (stranded mRNA and total RNA, 56-74 GB each); HG004 has the same design.
  Proposal: the 50 FORTE BAMs genome-wide with labels that need no DNA (absent from a population catalogue =
  hom-REF; common alleles with carriers judged by their own reads; wrong-allele rows), and HG002 for what those labels
  cannot decide (MNP, delins, 50+ bp events, 0.99+ paralogs, a second aligner). The operator: "We can start both".
- 2026-11-19: **both started, locally.** The HPC guard refused a listing of the shared pipeline-resource folder (outside
  the readable folders; the operator's setting, not changed). Public sizes read: GATK's AF-only gnomAD GRCh38
  (`af-only-gnomad.hg38.vcf.gz` 3,184,275,189 bytes, `.tbi` 2,443,190); GIAB HG002 GM24385 mRNA FASTQ R1 7,321,569,102
  and R2 7,375,938,930 bytes, `md5.in` 390. HG002 rows made (`kit/local/h333_rows.py`, public files only; 523,342
  distinct rows): every exonic truth allele (SNV 158,988; INS 14,722; DEL 14,482; 50+ bp 844), cis MNP 1,382 and
  delins 1,215, random hom-REF 24,000 (6,000 per distance stratum), paralog 124,274 (hom-REF 120,997), wrong-allele
  186,948. Slip: I first called 233 paralog rows mislabelled; I had left out the panel rows when comparing counts
  (the rerun is equal). FORTE genome-wide maker written (`kit/local/g333_rows.py`); it waits for the catalogue.
- 2026-11-19: **download job submitted** (the operator: "But we can subit a job to download"): `auto-gbcms-dl333`
  12393133 (`kit/dl333_fetch.sbatch`) fetches the HG002 FASTQs (MD5 against GIAB's `md5.in`) and the AF-only gnomAD
  GRCh38 file (size check) into `inputs/public/` on HPC, after a 390-byte probe of internet access from the node; it
  reports tabix, pysam and free space, so the FORTE row maker can run on HPC. The first submission was refused by the
  HPC guard (the job name must be on the sbatch line); nothing was submitted by it.
- 2026-11-20: **downloads checked** (job 12393133, exit 0, 9 min 45 s): both HG002 FASTQs match GIAB's MD5s, the
  gnomAD file and index match the server's sizes. tabix, bgzip and samtools are in the operator's samtools
  environment (1.23 on a node; the kit's `site.sh` names them; slip LRN-20261120-001: I first reported "no tabix").
  The FORTE row maker takes `--tabix` and keeps the alleles of the last catalogue position (a missing flush).
- 2026-11-20: **the operator's rule** (working agreements): "we should solve problem at it core/source keeping
  accuracy in mind, and not patches unless i explicitly ask for patches". The rules decided here go into the core
  judge code; the prototype's switches are measurement tools only.
- 2026-11-20: **HPC cleanup** (the operator deletes; about 48 GB freed; every #333 run kept). Slip
  (LRN-20261120-003): a state-filtered `sacct` listed no running job, and the operator emptied `gbcms_tmp` while
  M-H's slice task 12340216_26 ran (11 h in; it timed out at its 12-hour limit an hour later).
- 2026-11-20: **M-G, D7, the paralog rows and M-I finished; M-H not yet.** Every job of M-G, D7, the paralog rows and
  M-I ended exit 0; their aggregates (tables and READMEs) and M-G's and M-H's timing are copied back. M-H: the full-BAM
  arms and the timing ended exit 0; the slices' last group (Rrt, Rej, Res, Reb; 12340513) runs its 15 slow-slice
  tasks since 12:10, and the aggregate `mHagg` waits. The overhang arm (Rov, `RNA_OVERHANG=spliced`) timed out on
  slow slice 7 (over 12 h against about 2.5 h for the other arms) and took 1.4-3x the time of the other arms on the
  other four slow slices. **Timing** (median of 3, CPU against the C++ tool, 5,000 rows x 10 BAMs; 1,000 x 1 in
  brackets): the judge 1.21-1.25 (0.60); R 1.35-1.36 (0.62); R with clip read off 1.25 (0.60), so clip read costs
  about 9% CPU; R with the SNV route 0.35 (0.23); R with NH doubt 1.35 (0.63). Only the SNV route is under the speed
  floor. **Identity finding, not reported before:** the build's rules path with no switch (gR) differs from develop's
  rules (dvR) on 1 slice indel line (the counts, ASJD and per-transcript columns), in M-E and M-G alike.
- 2026-11-20: **results of M-G, D7, the paralog rows and M-I** (aggregates through the sanitizer; tables extracted by
  a Sonnet agent and spot-checked; no missing task and no broken count relation in any set). **D2, clip read** (R
  against R with clip read off): own alleles (M-A) +2,675 ALT (+0.09%); DNA-check hets +822 (+0.74%); false ALT M-B
  +413 (+1.1%), M-A hom-REF +102, DNA clean-absent rows 94 to 121; CPU +9%. **D1** (R against clip admission off):
  own alleles +49,737 (+1.6%); M-B false ALT +2,383 (+7.0%); DNA clean-absent rows 1 to 121 (22 rows with ALT 2+).
  **D5, the mirror** under R: own alleles -50 ALT, M-B -16. **D7, the end bar on real DNA** (the DNA judge on the 40
  tumour BAMs): bar4 raises false ALT at M-B hom-REF rows from 1,356 to 5,714 and at M-A hom-REF from 11,062 to
  12,881, and own-allele ALT by 0.02%; bar2 about half of that. The new build's DNA judge equals develop's on every
  line (58 MAFs, both sets). **D6, the foreign test** at paralog rows: hom-REF false ALT 2,994 under R, -372 (-12.4%)
  with share:0.1/dense:4, -480 with on, -421 with share:0.2, all from bands below 0.99 identity (the 0.99+ band holds
  73 false ALT under every arm); true ALT at paralog rows that carry the variant -43 (-0.24%); NH doubt changes
  nothing. **D4, the edge rule** (M-I: 16 tumour indel rows whose repeat run reaches an exon edge, 259 row-samples):
  ALT unchanged (1,655); REF 57,749 to 75,099. At deletions REF barely moves (55,478 to 55,156); at insertions it rises
  8.8x (2,271 to 19,943) and the RNA VAF falls from 0.171 to 0.023. own_junction equals motif on every row;
  stop_at_edge and stop_balanced drop 234 to 263 true ALT at flush-edge rows. The rule was built for deletions; why it
  credits REF at insertions needs a read-level look before D4 is decided.
- 2026-11-20: **DNA truth deleted, then restored.** My cleanup list held `runs/p292b` and `runs/p292c` (#292's DNA
  check and slice indel rows), which #333's arm, aggregate, bins and DNA tumour jobs read (LRN-20261120-004). M-H's
  last slice array (12340513, 20 tasks) ran without its indel set and ended exit 0 (`rows=missing`); the M-H
  aggregate (12341533) ran before the restore, ended exit 0 with `fatal=0`, and has no DNA labels: it is not valid.
  Restored by job 12550026 from the volume's snapshot of midnight that day: 2,080 and 2,467 files, every file equal
  by path, size and sha256. Both kit steps fail open on a missing input (the arm script skips, the aggregator returns
  no labels): a class to fix at the source.
- 2026-11-20: **read-level checks** (Opus agent scripts, dry-run on a simulated tree, reviewed; jobs 12550699 and
  12550707, exit 0). **The edge rule at insertions:** 17,641 of the 17,672 changed reads sit at one locus,
  chr9:8,331,574 (a 21-bp insertion and two variants of it, 33 row-samples, 7 bases into the intron at an exon edge
  of a minus-strand gene). Each is a spliced read whose N spans the insertion site and whose exon bases match the
  reference (for example `86M11615N10M`). Edge rule off: the judge leaves it out of DP. Motif: it judges the read
  at another placement of the insertion, inside the exon, and credits REF (`other_placement_judged`, `fits_ref`).
  Both alleles can splice to the same mRNA, so these bases cannot tell them apart; the DNA disagrees (for example
  339 REF / 149 ALT in DNA against 40 REF / 0 ALT under R). The other four positions: 31 reads at a 34-bp insertion
  (REF under R, neither with the rule off), none at a 1-bp insertion. **The rules-path line:** slice 30,
  chr2:178,678,414 CTT>C (a TT deletion at the run of four Ts at an exon start); one read (`28M207N68M`, its junction
  ending inside the run; after the junction `A` then `TT` then `GCAGG`, where REF has `C TTTT GCAGG` and ALT
  `C TT GCAGG`): develop's rules call it neither, the prototype's rules (no switch) ALT, both through the exact-carrier
  rule. The prototype changes shared code that the rules path reads; the core implementation must keep the rules
  path identical to develop.
- 2026-11-20: **M-H results** (aggregate `mHagg2`, valid: the fail-closed aggregator found every DNA truth file;
  6,736 arm runs; only the overhang arm misses slice 7). Each arm is R plus one built option; deltas against R.
  **Foreign test** (share:0.1/dense:4, on, share:0.2): own alleles -4,337 to -10,063 ALT, false ALT at M-B -597 to
  -652 and at M-A hom-REF -1,520 to -2,367: it removes about twice as many true ALT reads as false ones on panel rows,
  at +64% CPU. **NH doubt:** no change (own alleles -108, M-B -31). **End bar** bar4 / bar2: false ALT M-B +1,732 /
  +914, DNA clean-absent rows 121 to 261 / 158, own alleles +239 / +43. **Overhang spliced:** own alleles -2,655,
  REF -238,689 at M-B rows (depth lost broadly), 1.4-4.6x slower. **Shape tests off:** within a few counts of R
  everywhere, at equal CPU. **SNVs by one base** and **the SNV route** give identical counts, equal to develop's rules
  on SNV rows: against R, own alleles -45,242 ALT (-1.44%), but false ALT M-B -1,621 (-4.4%), M-A hom-REF -1,261
  (-2.7%), DNA clean-absent rows 121 to 0; at own alleles of DNA VAF 0.02-0.1 the route finds ALT in 922 rows against
  R's 864 (+58, +6.7%) and has 103 more rows at depth 10+; 0.35x the C++ tool's CPU against R's 1.36x. **Edge
  variants:** own_junction equals motif on every row; stop_at_edge equals motif on these sets; stop_balanced loses
  8,801 REF at slice indel hom-REF rows (M-I: it and stop_at_edge also drop true ALT at flush-edge rows).
- 2026-11-22: **truth runs, HG002 and FORTE, complete** (the operator's accuracy rule of 2026-11-20: a rule ships when
  it lowers false ALT and missed ALT on rows of known truth). HG002 RNA against GIAB v5.0q: 30 chunks, 270 runs, 523,342
  rows in every arm; missed ALT (truth-ALT row at depth 10+ with ALT 0) / false ALT (truth-REF row with ALT 2+): rules
  1,510 / 5,077; the judge alone 1,203 / 5,073; R 1,097 / 5,217; R without clip admission 1,193 / 5,084; R with the SNV
  route 1,175 / 5,265; R with the foreign test 1,217 / 4,972. FORTE (g333agg, 43 patients, exit 0, nothing withheld):
  missed carriers (pile-up rule P3) / false ALT row-samples (depth 10+, ALT 2+, at hom-REF, paralog-absent, noncarrier
  and wrong-allele rows): rules 1,762 / 25,725; the judge alone 1,833 / 25,416; R 1,819 / 26,687; R without clip
  admission 1,835 / 25,388; R with the SNV route 1,764 / 25,770; R with the foreign test 1,988 / 26,313. What each rule
  adds (HG002 missed / false; FORTE missed / false): clip admission -96 / +133, -16 / +1,299; clip read -13 / +20,
  -3 / +312; the edge rule 0 / 0, 0 / -18; the mirror 0 / -1, +1 / -8; the foreign test +120 / -245, +169 / -374 (it
  removes 1.22 M of 7.74 M truth-ALT reads at HG002, MNP rows -57%, rows with variants in cis -52%); the SNV route
  +78 / +48 at HG002. No arm beats the rules on both sets; the judge alone comes closest, with one class to read: 71
  FORTE carriers it misses (none within 1 base of a splice site). CPU summed: HG002 the judge alone 83.5x the rules, R
  89.8x; chunks c11, c12 and c19 233-262x, where one bin holds one thread (c11: 33,103 s CPU, 32,990 s wall, on 4 CPUs;
  ordinary chunks keep 3.2-3.3 busy); FORTE 5.7x and 6.4x. Peak memory on HG002: R 4.21 GB against 2.57 GB (c12). The
  decisions page is updated (D1, D2, D6, D7: out of the code; D4: fix the spanning-junction class first; D5 on with no
  switch; D8 opt-in with exact speed-ups at deep loci; D9 one judge). Scores: `out/score_hg333_full`, `out/g333agg`.
- 2026-11-22: **read checks A, B, C** (the operator: "i approve all three"; traced one-row runs, each row-sample's
  totals checked against its own MAF; FORTE read lists stay on HPC, only sums come back). **A, the 71 FORTE carriers
  the judge misses** (jobs rlC1, rlCw): all SNVs, 32 patients, 1,781 reads ALT under the rules and neither under the
  judge, 1,750 of them the judge's "another allele". At 45 row-samples the reads' own alignments write an indel beside
  the SNV (a 1-base deletion 3-5 bases before it, or a 1-3-base insertion 1-2 bases before it, by 5% to over 50% of
  the site's reads), and the judge keeps it as an allele: 1,512 moved reads fit that indel allele exactly (score 0),
  and 1,498 of them lose 7-20 points on the SNV allele. Their bases spell the indel, not the SNV; the rules see ALT
  because the aligner did not write the indel as a gap, so the read's shifted bases sit on the row as mismatches or in
  a soft clip (1,715 of the 1,750 have NM 1-4; the row is within 10 bases of the read end in only 39% of the 1,436
  listed moved reads, so "near the read end" does not hold). 110 fit the indel allele with
  one mismatch (they could carry both; no state holds the SNV and a found indel together). At the other 26
  row-samples the moved reads fit no allele well (-7 or worse on all three). The judge never credits a carrier the
  rules miss (reverse 0). So these are not judge misses; P3, built from the same aligned bases, labels them carriers.
  **B, clipped reads at FORTE expect-0 rows** (rlC2, no cap: 1,351 row-samples, 2,923 runs, 42 patients): clip
  admission credits 3,074 ALT reads at 895 row-samples, 90% of them through an end clip of 1-5 bases whose one base at
  the row is ALT, with one clean base past it, 75% with no annotated exon edge in reach: read-end mismatches the
  aligner clipped. Clip read credits 1,218 ALT reads at 677 row-samples, 93% with the clip starting 1-5 bases from an
  annotated exon edge (exon overhangs). **C, HG002 wrong-allele rows R credits and the judge does not** (rlC3, public,
  126 rows): every read that changes is a clipped read; of 198 listed R-only ALT reads, 151 spell the asked (wrong)
  allele with their own bases, 130 of them reads starting inside a poly-A or poly-T run with the flank before the run
  clipped (85 of their 93 indel rows hold a run of 10 bases or more: a run-length error seen through the clip), and 36
  spell no allele (mostly reads ending inside a long run, where clip read hands the clip to the judge and the longest
  allele fits best). The judge alone calls 151 of the 198 ALT too and then does not count them, as clip admission is
  off (f0's traces show the same: the clip is substituted and judged, then dropped; computed then discarded).
  Results: `out/rlC1`, `out/rlC2`, `out/rlC3`.
- 2026-11-22: **the DNA side of check A** (job rlCd; the operator: "i approve both, read low-VAF rows and the DNA";
  the matched IMPACT tumour BAMs in place, tumour lines only, counts only). Check A's SNV rows and the 46 found
  alleles the judge keeps at 45 row-samples (the trace's own count, equal to rlCw's), lifted to b37 by sequence (103
  of 103, no REF mismatch), counted with develop's `gbcms dna` in 39 tumour BAMs of 29 patients (8 row-samples are
  speed-batch slots with no DNA). The DNA cannot decide the class: depth is 0 at all 39 tumour row-samples with a kept
  found allele, at the SNV and at the indel, and at 59 of the 63 SNV rows; the sites are off the IMPACT panel. The 4
  SNV rows with depth (10-29 to 100+) carry the SNV at VAF 0.65+, and none of them has a found allele. Results:
  `out/rlCd`.
- 2026-11-22: **check L, own alleles at tumour VAF 0.02-0.1** (jobs rlC4, rlCw4; the same approval; check A's method
  at M-A's rows, each one-row VCF also holding the union rows within 150 bases, so the judge sees the M-A run's
  neighbours; the why aggregate now reads each trace at the row's own site, unchanged on one-row traces). The same
  10,217 own-allele row-samples as M-G; ALT reads at any depth: rules 11.9k, the judge 9,604 (-19.0%; M-G at depth
  10+: 11,747 against 9,524). 252 row-samples in 38 patients have more ALT under the rules (112 with judge ALT 0:
  85 of these have 1-2 rules ALT, 16 have 6+); 30 have more under the judge (+198 reads). The traced runs give the M-A
  counts at 252 of 252. Read moves: ALT to neither 2,451 and ALT to REF 110; REF to ALT 98, neither to ALT 15; REF
  to neither 7,110 (the judge drops REF reads too). This is mostly not check A's class: 15 row-samples keep a found
  deletion (16 alleles, 1-3 bases, 1-5 bases before the row, 5-50% of the reads), with about 500 moved reads; 344 of
  them fit the ALT and the deletion equally (both score 0, REF -4 to -6), and the judge calls them another allele.
  The main part: the judge calls 1,426 of the 2,489 moved reads Depth, most with all three scores -7 or worse (944
  with NM 5+): reads that fit no allele, which the rules count as ALT by the one base at the row. Of the 2,451 ALT to
  neither reads, the bases of 2,263 (92%) spell neither allele over the window (1,354 with 3+ other differences, 294
  at MAPQ 1-9, 531 through a clip); 188 spell the ALT, and the judge drops them. Rows with ALT at depth 10+ fall from
  922 to 836 at rows the tumour DNA confirms. Open: read the examples (bases and positions, no IDs) of the 16
  row-samples with rules ALT 6+ and judge ALT 0 and of the 188 reads that spell the ALT. Results: `out/rlC4`.
- 2026-11-22: **check L's reads, one by one** (job rlCx, four runs; the operator: "Read examples in detail"; bases,
  offsets, CIGARs and calls through the sanitizer, read on HPC; the job writes no file). All 252 row-samples and all
  2,561 lost reads (rules ALT, judge not ALT) are classified; the counts equal rlC4's. **Where:** 226 of the 252
  row-samples (90%) sit in a tandem repeat or run of 8+ bases (12+: 141, 8-11: 85); 15 hold a found indel at the
  row's site. **How the judge scores them:** all three tie 686 (295 by the unexplained-indel rule); ALT ties another
  allele 512; ALT ties REF 432; ALT best by 1-3, under the margin of 4, 466; REF fits better 245 (110 to REF); another
  allele fits better 180; ALT best by 4+ 40, all by the unexplained-indel rule. Four mechanisms explain the examples:
  1. **An indel allele cut by the read end.** The aligner writes a read that ends inside a repeat as the SNV. Example:
     an X-chromosome CT repeat, T>C, rules ALT 37 and 54, judge 0, at about 1,000 reads, with a found 2-base deletion
     (CTT>C) in 10-20% of the reads. Every rules-ALT read starts 7-9 bases before the row (96M, NM 1, the one
     difference +0:T>C), and its first bases fit the deletion allele equally (alt 0, other 0, ref -4). No read that
     spans the repeat carries T>C: judge ALT 0 against REF 673 and 983. The same at a CAG repeat with a one-unit
     contraction (G>A, rules 199 and 93, judge 14 and 2) and at a 15-base deletion (G>T, rules 70 and 34, judge 2 and
     0). Of the 190 reads whose own bases spell the ALT over the window, 95 are this case: an aligned end within 10
     bases of the row and a found indel at the site.
  2. **An indel behind a clip.** At an (AT)8 repeat next to a T13 run (T>A, about 20 row-samples), reads clip 3-8
     bases, for example GGTA at -18..-15 where the genome reads TAAT and has GGTA two bases to the left. The read is
     REF with two T fewer; the aligner kept its right end, so the A of the last AT unit lands on the row's T. Every
     haplotype scores -8: one 2-base gap that covers the row's base, so the read has no base at the row. Many of these
     reads have MAPQ 3.
  3. **Reads from another sequence.** 1,287 lost reads fit no allele (Depth, all three scores -7 or worse): reads
     with 6-10 differences, often the same shifted haplotype in many reads (one GC-rich locus at 3 row-samples, where
     the judge's REF falls too; an X CAG/CAA region, rules ALT 29, judge 0, about one mismatch every 3 bases).
  4. **ALT with a run of the wrong length** (the 40 reads where ALT fits best by 4+). 32 sit at SNV rows on the
     junction of two runs, where the ALT base joins or extends a run (A>T at A11|T10 gives A10 T11; A>T at
     ...TATA|T13 gives T15; T>C at CCC|T|CCCCCCC gives C11). These reads carry the ALT base and a run 1-2 bases
     longer than ALT's: 25 are written as the SNV plus a 1-2-base insertion at the row (26M1I69M, NM 2), 2 as a
     second substitution next to it, 5 hide the change behind a clip. Example: alt -7, ref -11, other -11; the best haplotype still needs a 1-base gap in the
     row's run, so the unexplained-indel rule calls another allele. The other 8 reads have 5-10 mismatches (poly-G
     read tails at one GC-rich locus, a TTC repeat, a G-rich run).
  **At a run junction the aligner decides what the rules count.** A read with ALT's run one base longer is written as
  the SNV plus an insertion, and the rules count it ALT by the base at the row. A read with ALT's run one base
  shorter is written as a deletion of the row's base, and the rules do not count it. So the rules count ALT+1 and
  never ALT-1; the judge counts neither. That is invariant 7 (judge bases, not placement) and the decided rule for
  wrong-length indels (`partial_alt`, not ALT). A read that starts inside the run cannot show the run's length, and
  the judge counts it ALT (alt -4, ref -8).
  **What this means for the measure.** By its own bases, almost every lost read does not carry the given SNV alone:
  it is an indel allele cut by the read end, an indel behind a clip, a read from another sequence, or a run of
  another length. The "rows the tumour DNA confirms" label (922 to 836 rows with ALT at depth 10+) is a DNA call at
  tumour VAF 0.02-0.1 at these same repeats, so it is not independent truth for what the RNA reads carry. On known
  truth (HG002 RNA, public, `out/score_hg333_full`), the judge lowers both errors at SNV rows: FN 810 to 797, FP 4,236
  to 4,119 (ALT 1+) and 2,975 to 2,898 (ALT 2+). Open: score HG002's SNV rows in repeats of 8+ bases and at run
  junctions by arm, to test mechanisms 1, 2 and 4 on truth. Results: the rlCx job logs, on HPC only.
- 2026-11-22: **the four mechanisms on truth** (HG002 RNA, public, this Mac; the operator: "Sure do it";
  `kit/local/score_hg333_rep.py`, which reuses score_hg333.py's rows, labels and metrics unchanged and adds each SNV
  row's reference context from GRCh38). 462,948 SNV rows scored (54 of 466,196 differ from the FASTA's REF and are
  left out; the all-SNV line then equals `out/score_hg333_full`: rules FN 810, FP 4,236 and 2,975). Rules against
  the judge, truth-ALT FN and truth-REF FP (ALT 1+):

  | SNV rows | truth-ALT rows | FN rules → judge | FP rules → judge |
  |---|---|---|---|
  | tract 12+ | 3,460 | 39 → 38 | 51 → 47 |
  | tract 8-11 | 59,812 | 320 → 313 | 1,585 → 1,538 |
  | tract under 8 | 95,689 | 451 → 445 | 2,600 → 2,534 |
  | a truth indel 0-10 bases away | 2,551 | 129 → 118 | 31 → 17 |
  | a truth indel 11-30 bases away | 2,609 | 22 → 20 | 35 → 25 |
  | the REF base in a run of 8+ | 1,429 | 18 → 13 | 10 → 7 |
  | the ALT base makes a run of 8+ (joins or extends one) | 1,252 | 3 → 7 | 22 → 19 |

  Mechanisms 1-3 hold on truth: in every repeat band and beside a truth indel the judge has fewer FN and fewer FP
  (beside an indel, ALT reads on truth-REF rows fall 410 to 186). Mechanism 4 is the one class that moves away from
  the truth: where the ALT base joins or extends a run of 8+, FN rise 3 to 7 (of 1,252) and the median VAF deficit
  (ALT/depth, not the decided VAF; see the correction in the next entry) grows (ALT run 12+ joining two runs: -0.18 to -0.31; 12+ extending: -0.22 to -0.34; 8-11: -0.10 to -0.13 and
  -0.09 to -0.10), with FP 22 to 19. The rules already fall short there (they lose ALT-1 reads as deletions); the
  judge also drops ALT+1. That is the decided stutter rule (`bam-evidence-caveats.md`: "Stutter in repeats: counted
  exactly; named by a diagnostic", measured on DNA MSI-high indel loci, #92), applied to RNA SNVs whose ALT makes a
  run, where RT stutter is more frequent. The class is small (1,252 of 158,961 truth-ALT SNV rows, 0.8%). Changing it
  is a count-affecting decision (decision protocol); not proposed here. Results: `out/score_hg333_rep`.
- 2026-11-22: **read check J, the junction rows read by read** (HG002 RNA, public; the operator: "Read HG002 junction
  reads"; jobs 13243112, 13244677, 13245744, 13246469, 13246801; `kit/hg333_rlcj.py`, `kit/hg333_rlcj.sbatch`,
  `kit/local/hg333_rlcj_rows.py`, `kit/local/hg333_rlcj_sum.py`). The 242 rows of the class where the rules (dvR) and
  the judge (fJ) count different ALT (230 truth-ALT, 12 truth-REF), each run as its own VCF with its neighbour rows
  within 320 bases: the judge takes the input's rows within 300 bases as context, so a one-row VCF changed fJ's counts
  at 110 of 242 rows; with the neighbours both arms equal hg333 at 242 of 242. Each read is classed by its own bases
  between the genome's flanks around the ALT run and the runs beside it (clips included), by run lengths against the
  REF and the ALT allele.
  - **What the judge drops** (truth-ALT reads ALT under dvR and not under fJ: 2,891): bases that match neither
    allele's base order 1,311 (45%; 73% hold a base under Q20); no flank, or the read starts or ends in the run, 744
    (26%; the missing flank is the side sequenced after the run: right flank 176 forward against 35 reverse, left flank
    261 reverse against 25 forward); the ALT allele with its run 1 base longer or shorter 351 (12%; 275 with no other
    run changed); 2-3 bases or more 267 (9%); ambiguous 93 (3%); a REF molecule 64 (2%); the exact ALT run 61 (2%,
    most with another run changed, which `indel_unexplained` sends to Other). fJ calls 98% of them neither and 2% REF.
    The other way: 117 reads.
  - **How each arm calls stutter reads** (all classed, counted reads): ALT+1 dvR ALT 395 of 1,729 (23%), fJ 106 (6%);
    ALT-1 dvR 115 of 458 (25%), fJ 54 (12%); exact REF reads ALT 0 in both. Row 13 (C>A at chr1:85,653,613; A9 C A2
    becomes A12): 44 reads ALT+1 against 22 exact ALT; an ALT+1 read (144M6S) has its extra A written as four shifted
    bases after the row, dvR counts ALT by the base at the row, fJ scores alt -7, ref -11, `indel_unexplained`,
    neither. Row 86 (G>T at chr6:31,197,235; T6 (GT)n becomes T8 (GT)n-1): 30 reads are T8 (GT)n, the ALT allele
    with one GT unit more or REF with two more T (each a 2-base insertion away); dvR counts one by placement (87M63S), fJ
    scores -52 on all three, neither.
  - **On truth**: at the 92 hom-ALT rows with depth (every molecule ALT, no allelic imbalance) REF reads fall 303 to 38
    and VAF (AD/(RD+AD)) is nearer 1.0 under fJ at 28 rows, under dvR at 1 (63 equal; mean |VAF - 1| 0.027 to 0.004);
    ALT reads 11,854 to 10,967 (the judge abstains, it does not call REF). At the 12 truth-REF rows ALT reads fall 604
    to 521 and rows with ALT 1+ 11 to 8 (reads ALT under dvR only 84, under fJ only 1). At the 124 het rows RNA allelic
    imbalance hides the truth (mean |VAF - 0.5| 0.21 and 0.26). FN 0 to 4: row 125 (dvR's ALT reads are 3 ALT+1),
    row 128 (1 ALT+2, 1 other), row 166 (a 30-base run, 8 reads other), row 86 (dvR's one ALT read is the read above).
    What-if, fJ plus the ALT+1 (or ALT+1 and ALT-1) reads it calls neither: FN 4 to 3, FP rows 8 to 8; FP at the 414
    truth-REF rows where both arms agree is not measured.
  - **Correction** to the entry above: its median VAF deficit used ALT/depth (`score_hg333.py`), which counts abstained
    (neither) reads against the ALT; the decided VAF is AD/(RD+AD). By that VAF the judge is nearer the truth at this
    class's hom rows.
  - **Open**: on truth-ALT rows fJ counts ALT for 141 of 538 counted reads the classifier calls REF+1 (the REF base
    present), and for none on truth-REF rows. The classifier ignores the neighbour rows the judge uses, which may explain
    it; check it per read with the judge's trace before any decision (invariant 7).
  - **Result**: here the judge follows invariant 7 and the decided stutter rule; most of the rules' extra ALT is credit by
    placement from reads whose bases cannot separate the alleles (ALT stutter is 12% of it). No count change proposed.
    Results: `out/rlCj2` (report.tsv, summary.txt; the earlier classifiers' reports kept as report.flank.tsv and
    report.rle1.tsv).
- 2026-11-22: **read check J, part 2: the REF-type reads the judge counts ALT** (the operator: "Yes, check them";
  job 13249085, `kit/hg333_rlcj_refalt.py`, public; every REF-type read of the 242 rows with its window, the input rows
  in it and the judge's `judge fits` fields; class equal to the report's at 4,026 of 4,026 reads). The judge counts ALT
  for 216 such reads, all on truth-ALT rows:
  - 138 are exactly the ALT allele in cis with one or two listed input rows (for example row 234, hom A>C: the read is
    G6 A C10, the listed GA>AC; my classifier compared only with the row's own REF and ALT). The judge is right; the
    classifier was wrong.
  - 14 more are on hom-ALT rows (every molecule ALT): right in outcome.
  - 64 on het rows are not proven either way. 6 (rows 45, 221, 118) hold a confident REF base at its own place: at row
    45 (G>A, A3 G A10 becoming A14) the read is A3 G A12, REF with 2 more A; no listed row puts A after the G, so REF
    pays the 2-base gap (-8), while ALT with a listed AA insertion before it pays only the G as a mismatch (-4), and
    the judge counts ALT. 54 hold the REF base 1-3 bases off the row with the runs beside it changed (28 confident, 26
    below the bar; 39 of them at rows 150 and 151, two het G>A 4 bases apart), a REF molecule with two run slips or ALT
    with substitution errors; the judge's margin is the one base at the row. 4 others.
  - **Cause, one cost order:** a gap costs 6+k anywhere, a confident mismatch 4, so in a long run the scorer ranks "a
    confident base is wrong" above "a run is 1-2 bases longer". With `indel_unexplained` the same order loses ALT
    molecules with a slip (part 1: ALT+1 alt -7, ref -11, Other) and here credits REF molecules with a slip. A
    repeat-aware gap cost (E12, §3) moves the false ALT toward the truth (option a); the missed ALT also needs a change
    to `indel_unexplained` (option b). Count-affecting and shared with DNA: the decision protocol, every data type,
    measured first. (Corrected 2026-11-23: this line first said "the planned E2, Lancet2's two-tier gap model" and
    that one change moves both. E2 is another item; Lancet2's two-tier gap is in its assembly only, and its genotyper
    keeps gaps strict.)
  - The open item of part 1 is closed: mostly the classifier; a small real class remains (6 reads likely false ALT, 54
    unproven, of 216).
- 2026-11-23: **the internal truth: labels that do not come from the RNA reads** (the operator: "Build thw internal
  truth and plan and research the gap model").
  - **v0, FORTE, from aggregates already here** (the mHagg2 aggregate; local scorer `kit/local/truth333_forte.py`;
    labels from the matched IMPACT DNA; RNA depth 10+ by dvR; like for like). Missed ALT = ALT 0 at a truth-ALT
    row-sample; false ALT = ALT 2+ at a truth-REF row-sample.

    | Truth-ALT set | Row-samples | dvR | fJ | gB1f | Rrt |
    |:--|--:|--:|--:|--:|--:|
    | germline het and hom | 572 | 0 | 0 | 0 | 0 |
    | tumour variant, DNA VAF 0.35+ | 12,113 | 239 | 242 | 238 | 239 |
    | tumour variant, DNA VAF 0.1-0.35 | 702 | 230 | 248 | 243 | 230 |
    | tumour variant, DNA VAF 0.02-0.1 | 6,071 | 5,121 | 5,193 | 5,164 | 5,121 |
    | missed ALT, VAF 0.1+ and germline | 13,387 | 469 | 490 | 481 | 469 |

    | Truth-REF set | Row-samples | dvR | fJ | gB1f | Rrt |
    |:--|--:|--:|--:|--:|--:|
    | DNA-check absent (a stress set: chosen where RNA shows ALT) | 245,485 | 171,510 | 169,623 | 169,581 | 171,233 |
    | tumour hom-REF, part ma | 105,273 | 4,573 | 4,431 | 4,653 | 4,571 |
    | tumour hom-REF, part mb | 340,627 | 5,629 | 5,479 | 5,899 | 5,626 |
    | slices, absent, SNV | 4,231 | 4,045 | 3,987 | 3,961 | 4,032 |
    | slices, absent, indel | 173 | 170 | 170 | 170 | 170 |
    | slices, hom-REF, SNV | 3,616 | 1 | 1 | 2 | 1 |
    | slices, hom-REF, indel | 1,752 | 68 | 69 | 71 | 71 |
    | false ALT without the stress set | 455,672 | 14,486 | 14,137 | 14,756 | 14,471 |

    - Reading: against DNA truth the judge (fJ) has the fewest false ALT (349 fewer than the rules, −2.4%) and 21 more
      missed ALT of 13,387 (18 of them at DNA VAF 0.1-0.35). The RNA switches (gB1f) add 270 false ALT; Rrt is the
      rules within 15. The slices' absent SNV rows have ALT 2+ in 96% of row-samples in every arm, like the stress set;
      how they were chosen is not checked here.
    - Mean RNA minus DNA VAF at germline hets: dvR −0.0232, fJ −0.0236, gB1f −0.0188.
    - The gap: germline truth-ALT is thin (571 het, 1 hom). The kit below adds DNA labels at the genome-wide run's
      panel rows.
    - The DNA side: d7bagg's labels are the rules' own counts, so they are not truth for DNA. An independent DNA truth
      needs the matched normals (not covered by the approvals) or HG002 DNA (a download).
  - **The kit for panel-row DNA labels** (not run: HPC was unreachable):
    - `truth333_rows.py` (public): g333's SNV rows of the sets het, wrong, ref and sd inside the panel territory, with
      their repeat bands (tract, ALT run, REF run; the definitions of `score_hg333_rep.py`, copied unchanged).
    - `truth333_dna.sbatch OUT BUILD KIND`: the rows lifted to b37, then a plain pile-up (`g333_pile.py`, base and
      mapping quality 20) of each index's matched tumours. KIND normal is marked not covered.
    - `truth333_agg.py` and `truth333_agg.sbatch`: aggregates only. The DNA label per index and row: `d_hom` (every
      deep tumour, depth 50+, at VAF 0.90+ with ALT 10+), `d_het` (VAF 0.15-0.85 with ALT 5+), `d_ref` (ALT 0),
      `d_amb`, `d_low`, `nodata`, `nolift`. Truth classes: het rows `d_hom` alt_hom, `d_het` alt_het, `d_ref` ref;
      wrong rows (a base no catalogue allele holds) `d_ref` ref, split by the partner het row's carrier state (the
      given-allele test); ref and sd rows `d_ref` ref, `d_het` or `d_hom` alt_found. Missed ALT, false ALT (ALT 1+ and
      2+), VAF and the four count relations, per set, truth class and repeat band, like for like.
  - **Dry run on simulated data** (seven made-up b37 DNA BAMs: a low-depth contig, an LOH tumour, a 5% contaminated
    tumour; four simulated RNA arms): rows 368 read, 237 kept (het 73, ref 80, sd 11, wrong 73), 235 lifted. An
    independent checker (separate code) recomputes the labels and the metrics from the piles and MAFs: 28 of 28 metric
    keys and 15 of 15 label cells equal; every het and ref row has the label its design predicts (LOH and
    contamination `d_amb`, low depth `d_low`). A rerun keeps every pile (cached); a stale pile is redone and the
    tables stay equal; the aggregate job gives the same tables and keeps the earlier ones as `.prev.<time>`.
  - **A fix found by the dry run:** a pile folder with no pile, or with a pile of other rows, was left out without a
    word, so the label came from the other tumours alone (an LOH tumour left out turns `d_amb` into `d_het`). The
    slot now fails and the aggregate exits 4.
  - **Approvals** (the operator, 2026-11-23): the tumour job, "Yes, approve it"; the normals and an HG002 DNA download,
    "you can read mathc normal and downloading HG002 DNA but no PHI"; the download's names and sizes, "Yes, download".
  - **The panel-row run** (jobs 13284559 tumour, 13284561 normal, 13284578 aggregate; all exit 0; aggregates and
    sanitized logs in `out/tr333agg`, summary `kit/local/truth333_sum.py`). Rows inside the panel: 1,001 (het 132,
    ref 678, sd 63, wrong 131), 999 lifted to b37; 58 tumour and 40 normal piles, 0 failed; 40 patients, every arm
    present at every row-sample; the four count relations 0 on every count. At RNA depth 10+ by dvR:

    | Truth class (DNA label) | Row-samples | dvR | fJ | gB1f | Rrt |
    |:--|--:|--:|--:|--:|--:|
    | missed ALT, germline het (`d_het`) | 1,058 | 38 | 38 | 38 | 38 |
    | missed ALT, germline hom (`d_hom`) | 445 | 0 | 0 | 0 | 0 |
    | mean RNA VAF, het / hom | | 0.447 / 0.993 | 0.448 / 0.994 | 0.451 / 0.994 | 0.447 / 0.993 |
    | false ALT 2+, common SNP rows, DNA hom-REF | 1,454 | 83 | 83 | 83 | 83 |
    | false ALT 2+, random exonic bases, DNA REF | 18,555 | 252 | 245 | 270 | 251 |
    | false ALT 2+, paralog rows, DNA REF | 1,105 | 38 | 38 | 38 | 38 |
    | false ALT 2+, a base no catalogue allele holds; patient carries the real allele | 1,366 | 23 | 22 | 24 | 23 |
    | the same; patient does not carry it | 1,307 | 14 | 14 | 14 | 14 |
    | the same; carrier state unknown | 515 | 6 | 6 | 6 | 6 |
    | false ALT 2+, all truth-REF | 24,302 | 416 | 408 | 435 | 415 |

    - Reading: on these rows (SNVs in coding exons, few in long runs) the four arms are almost equal. The judge has 8
      fewer false ALT than the rules (−1.9%), the RNA switches 19 more (+4.6%), Rrt 1 fewer; missed ALT is equal. This
      agrees in direction with v0. The given-allele test passes in every arm: where the patient carries the real
      allele, the wrong base gets ALT 2+ in 1.7% of row-samples, against 1.1% where the patient does not.
    - Missed ALT at germline hets is the same in every arm, so it is not a counting difference (allele-specific
      expression is the likely cause): 24 of 391 in runs of 8-11 (6.1%), 14 of 667 in shorter runs (2.1%).
    - Not truth: the paralog rows' DNA ALT (`alt_found`, 63 row-samples, 47 of 54 in runs of 8-11, RNA VAF 0.13) looks
      like paralog reads in the DNA; left out of the counts above.
    - Labels: tumour and normal agree at 4,023 of the 5,280 het-row row-samples (`d_het`/`n_het` 1,572, `d_hom`/`n_hom`
      633, `d_ref`/`n_ref` 1,818); 630 are `d_amb` where the normal is REF (the tumour has a low ALT fraction at a
      common SNP), not yet explained.
    - What this cannot see: E12's class. Runs of 12 or more hold 66 truth-REF row-samples and no truth-ALT. E12's M1
      and M2 need the HG002 DNA exome (job 13284844, downloading) and the repeat rows.
  - **HG002 DNA in hand** (job 13284844, `kit/dl333d_fetch.sbatch`; public): GIAB's Oslo University Hospital exome
    (Agilent SureSelect v5, HiSeq 2500, bwa-mem to GRCh37, duplicates marked; 10,339,495,604 bytes) and the GIAB
    v4.2.1 GRCh37 benchmark (VCF, tbi, BED), in `inputs/public/giab_hg002_dna_exome/` and
    `inputs/public/giab_hg002_v421_grch37/`. Sizes equal to the server's, MD5s equal to GIAB's index, b37 reference
    names and lengths. A capture library with PCR, like IMPACT: the DNA truth for the judge against the rules, and
    E12's M1 on DNA.
  - **The HG002 DNA truth score** (2026-11-23; the operator: "Yes, run it"; public data only). Kit
    `kit/hgd333_rows.py`, `hgd333_slip.py`, `hgd333_prep.sbatch`, `hgd333_arms.sbatch`, scorer
    `kit/local/score_hgd333.py`; jobs 13288607 (prep) and 13289325 (282 arm tasks), all exit 0; copy-back in
    `pub/hgd333/`. A dry run on a made-up tree with planted tracts and slips passed first (tracts 792 of 792, truth rows
    73 of 73, slip rates 5.2% and 1.5% against 5% and 1.5% planted).
    - **Rows.** Territory: the v4.2.1 benchmark BED inside the covered exome (depth 10+ at mapping quality 20,
      chromosomes 1-22), 110.05 Mb. Layers in separate VCFs: truth alleles (var, 126,124), phased merged MNP and
      DELINS (cis, 362, all hom: most genotypes are not phased), a wrong allele beside each truth allele (wrong,
      125,676), random REF bases far from truth (ref, 30,000), a one-unit deletion and insertion at REF repeat tracts
      of period 1-3 and length 6+ (tdel and tins, 32,612 each), and the base after a homopolymer changed to the run's
      base (tsnv, 9,468). Arms: dvR (develop, `--read-judge rules`), dvJ (develop, `--read-judge sequence`), fJ
      (the #333 branch, `--read-judge sequence`); `gbcms dna` defaults.
    - **dvJ and fJ are equal** on all 356,854 rows, in every count field: the #333 branch changes no DNA count.
      The four count relations hold on every row of every arm.
    - **Result at depth 10+ (dvR → dvJ):**

      | Truth class | Rows | dvR | dvJ | Change |
      |:--|--:|--:|--:|--:|
      | missed ALT, het | 76,526 | 457 | 334 | −123 |
      | missed ALT, hom | 48,388 | 62 | 26 | −36 |
      | missed ALT, multi-allelic | 1,205 | 35 | 22 | −13 |
      | mean VAF, het / hom | | 0.485 / 0.995 | 0.488 / 0.996 | |
      | false ALT 2+, wrong allele, partner het | 76,526 | 1,029 | 754 | −275 |
      | false ALT 2+, wrong allele, partner hom | 48,388 | 1,143 | 669 | −474 |
      | false ALT 2+, wrong allele, partner multi | 756 | 176 | 137 | −39 |
      | false ALT 2+, random REF bases | 30,000 | 55 | 56 | +1 |
      | false ALT 2+, one-unit deletion at a REF tract | 32,612 | 4,345 | 4,498 | +153 |
      | false ALT 2+, one-unit insertion at a REF tract | 32,612 | 1,731 | 2,005 | +274 |
      | false ALT 2+, SNV after a homopolymer (E12's class) | 9,468 | 313 | 559 | +246 |
      | false ALT 2+, all truth-REF | 230,362 | 8,792 | 8,678 | −114 |

    - Reading: the judge finds 172 more truth alleles (−31% missed ALT) and gives 788 fewer false ALT at wrong
      alleles. It gives 673 more ALT at REF tracts. At tdel and tins rows most ALT reads are real slip reads, which
      the stutter rule counts by design (#92), so those rows are not false calls by themselves; the check there is
      AD against the reads that carry the slip (below). The tsnv rise is E12's class: an SNV after a long run gets
      ALT from reads that slip, mostly in runs of 12+ (ALT run 12+ that extends the tract 88 → 236; that joins two
      tracts 198 → 282).
    - **E12's M1 on DNA** (`hgd333_slip.py`: reads of mapping quality 20+, primary, not duplicates, that span the REF
      tract with 5 anchor bases; one-unit deletions and insertions per 1,000 spanning reads; tracts with events in
      25%+ of 8+ reads left out in brackets):

      | Tract | Deletions | Insertions |
      |:--|--:|--:|
      | A/T run of 6 | 0.28 | 0.01 |
      | A/T run of 8 | 1.65 | 0.71 |
      | A/T run of 10 | 16.9 | 4.2 |
      | A/T run of 12 | 41.7 | 12.6 |
      | A/T run of 14 | 80 (75) | 26 |
      | A/T run of 15 | 103 (92) | 37 |
      | A/T run of 16 | 158 (113) | 46 |
      | C/G run of 10 | 55 | 6.5 |
      | period 2, 12 bases | 3.7 | 0.2 |
      | period 2, 16 bases | 20 | 0.3 |
      | period 2, 20 bases | 41 | – |
      | period 3, 15 bases | 1.5 | – |
      | period 3, 21+ bases | 9 | – |

      On the judge's scale (cost = 0.16 × Phred of the rate) a one-unit slip in an A/T run costs about 5.7, 4.2, 2.7,
      2.0 and 1.4 at lengths 6, 8, 10, 12 and 15. The survey's PCR values are 5, 4, 3, 2 and 2; today's cost is 7 at
      every length. Deletions are 3-4 times as frequent as insertions. C/G runs slip more than A/T runs of the same
      length. So M1 confirms the gap model's premise on DNA: the flat cost 7 is too high for runs of 10+.
    - **AD against M1's clean slip reads** at REF tract rows (period 1, high tracts left out): AD is 1.2-2.0 times
      M1's count for dvR and 1.24-2.5 times for dvJ, and the ratio rises with run length (deletions, runs 12-14: 1.31
      against 1.42; insertions, runs 15-19: 2.02 against 2.51). RD is 1.02-1.09 times the spanning reads in both. Part
      of the excess is M1's 5-base anchor (reads that end in the tract are not in M1); the rest needs a read check.
    - **Time** (summed over 94 chunks; CPU / wall / peak memory): dvR 12.9 min / 4.4 min / 0.24 GB; dvJ 95.9 / 25.9
      / 0.29; fJ 95.7 / 25.9 / 0.30. The judge takes 7.4 times the rules' CPU on these exome rows.
    - **Limits.** One sample, one library (PCR capture, bwa-mem, no ABRA2 realignment, unlike IMPACT and ACCESS). cis
      rows are hom only. Truth outside the v4.2.1 BED is not scored. Tract rows are REF by truth; most of their ALT
      is slip reads, so false ALT there mostly measures the stutter floor (the read check below measures the rest).
  - **The read check at REF tract rows** (2026-11-23; the operator: "Yes, run it", then "Yes, trace them"; public
    data only). Kit `kit/hgd333_reads.py` (select, report, selftest: 18 bases cases, 6 masked cases, 5 clip cases),
    `hgd333_reads.sbatch`; jobs 13316156 (traced runs, 6 tasks) and the reports (13322742 the first to pass its checks;
    13340383, 13342413, 13345884, 13346360 and 13346986 add the trace dump and the masked check; reports 1-5 failed
    their checks, below); copy-back in `pub/rc333/`.
    - **Rows.** Period-1 REF tract rows of the truth score, depth 10+, high tracts left out, where AD in either arm
      exceeds M1's clean slip reads (tsnv: AD 2+), up to 30 per length band: 161 tdel, 160 tins, 125 tsnv. Each set
      reran in both arms with `--trace`, with the rows' neighbours within 600 bases (the judge reads every listed row
      within 300 bases as context, so a row's counts depend on its listed neighbours; this is by design, and binning
      invariance holds for a fixed set of rows).
    - **Checks** (all pass, in every report from 13322742 on): the traced runs equal the arm runs on 446 of 446 rows;
      the per-read sums (REF, ALT, DP) equal the engine's at 446 of 446 rows in both arms; M1 recounted from the reads
      equals the truth score's at 446 of 446; every call line finds its read.
    - **Correction (my error).** The first reading of this check was wrong. Its bases view placed each read by its
      CIGAR and counted every mismatch, also on bases below the quality floor (`--min-baseq` 20), which the engine
      reads as unknown. The trace of its 60 "firm" contradicting credits showed the flaw: 33 of the 60 fit the ALT once
      their low-quality bases are masked. That view gave 17.6-36.1% of ALT reads as not carried, the judge's share
      always the higher, and it said that most of the judge's added rows with ALT 2+ at tins and tsnv (18 of 24, 23 of
      40) were wrong credits. The masked check below replaces that table, its reading and its count of the judge's
      additions; the true counts are 6 of 24 and 14 of 40.
    - **The masked check.** Each read counted ALT in either arm is aligned by its bases alone (its CIGAR plays no part)
      to REF, the ALT and the other lengths of the run (one and two units each way; RJ-27 names these as alleles that a
      read must tell apart). It uses the judge's scores: a mismatch costs 4, a gap of k bases costs 6+k, and a base
      below quality 20 (or N) costs nothing. Soft-clipped bases count (RJ-13). The window is 30 bases each side of the
      tract. The ALT "holds" when it beats every other haplotype by 4, the judge's margin. The check gives the judge's
      own REF-to-ALT score gap at 57 of the 60 traced reads; the other 3 have bases past the check's window (for
      example 100S26M).

      | Set | Arm | ALT reads | ALT holds | Ties REF / another length | Fits worse than REF / another length | Not held | Rows AD 2+ (under 2 reads that hold) |
      |:--|:--|--:|--:|:--|:--|--:|:--|
      | tdel | dvR | 563 | 509 | 19 / 17 | 5 / 13 | 9.6% | 104 (4) |
      | tdel | dvJ | 652 | 582 | 33 / 27 | 5 / 5 | 10.7% | 108 (5) |
      | tins | dvR | 310 | 218 | 22 / 42 | 2 / 26 | 29.7% | 62 (14) |
      | tins | dvJ | 370 | 290 | 4 / 73 | 0 / 3 | 21.6% | 86 (20) |
      | tsnv | dvR | 266 | 190 | 0 / 64 | 0 / 12 | 28.6% | 78 (25) |
      | tsnv | dvJ | 396 | 274 | 0 / 118 | 0 / 4 | 30.8% | 118 (39) |

      (At tsnv "another length" is mostly the one-unit insertion: it reads as the SNV up to the second base past the
      run.)
    - **Where the arms differ** (reads counted ALT, by the arm that counts them):

      | Set | Both arms: hold / not | Rules only: hold / not | Judge only: hold / not |
      |:--|:--|:--|:--|
      | tdel | 509 / 22 | 0 / 32 | 73 / 48 |
      | tins | 217 / 40 | 1 / 52 | 73 / 40 |
      | tsnv | 189 / 58 | 1 / 18 | 85 / 64 |

      The judge removes 104 of the rules' ALT credits, and 102 of them do not hold. It adds 383, and 231 of them hold:
      slip reads that the aligner wrote as a mismatch at a read end or as a soft clip, and at tsnv reads that the rules
      do not count and the judge admits by their soft clip (82 of 140). Of the 152 added credits that do not hold, 110
      tie another length of the run. The judge adds 4, 24 and 40 rows with ALT 2+, and 1, 6 and 14 of them have fewer
      than 2 reads that hold the ALT. So most of the judge's added ALT is real slip reads (the stutter floor), not
      wrong credits. This reverses the first reading at tins and tsnv.
    - **The runs beside the tract.** A second verdict adds the one-unit slips of the run before the tract and of the
      run after it as candidates. Many reads that hold the ALT then tie such a slip: the one base that tells them apart
      is below the floor, or the read has no base there. At tdel this is 223 of the rules' and 295 of the judge's ALT
      reads. Most are beside a single base (154 and 201), whose loss is not a repeat's slip; beside a run of 2+ bases
      there are 69 and 94 (tins 8 and 20, tsnv 4 and 4). Both arms credit these reads. No decided rule names this case:
      RJ-27 names the run's own lengths, and RJ-3 a read's first base.
    - **The engine causes** (traced on the 60 reads with their trace lines and the code; counted by the check above):
      1. *Judge: a read that cannot tell the ALT from another length of the run is credited ALT* (110 of the
         judge-only credits). At a short pure indel the judge decides by the margin over REF, the ALT and the listed or
         found alleles (`decide()`, `judge/decide.rs:369`, called from `conclude_margin`, `judge/decide.rs:904`). The
         run's other lengths are not candidates. The window check that makes a read inside the run depth only runs
         only at insertions of 50 bp or more (`judge/decide.rs:988`, RJ-23). RJ-27 (decided) asks every pure indel's
         read to hold the tract and read past it; the rules meet it through RJ-3's mask (`carrier.rs:748-765`).
         Example, tins k=173 (126M, every base q30+): the read's first base sits on the flank T and reads A, then come
         7 A's and the right flank. It fits INS1 and INS2 alike. The rules mask that base and make the read depth only;
         the judge gives ALT (0 against REF −4). At tsnv both arms credit the SNV to reads that end before the second
         base past the run (58 reads in both arms, 60 more in the judge).
      2. *Rules: a read that carries another indel next to the tract is credited ALT.* `own_event`
         (`carrier.rs:3534`) places the read's own length change at each offset and counts only confident mismatches,
         up to one (`carrier.rs:3616`, `carrier.rs:3633`). The change is inside when any placement with the fewest
         errors touches the window, the tract plus one base each side (`carrier.rs:3646-3648`, `carrier.rs:3677`;
         `window.rs:185`), so a deletion of the flank base is inside. `spells_alt` then forgives one confident
         mismatch against the ALT (`carrier.rs:3681-3685`), and the read is `Given`, so it counts ALT. Example, tdel
         k=34 (107M1D19M, every base q30+): the read deletes one T of the TTT before the A tract. It fits that
         deletion exactly and the ALT with one mismatch. The rules give ALT; the judge gives Depth (−8 against −11,
         under the margin). In tins k=165 (47M1I78M) the aligner put the insertion 17 bases before the tract, past
         bases of quality 10-29; masked, the read fits an 11th C in the run beside the tract best, and the rules give
         ALT. Of the 102 rules-only credits that do not hold, 49 have an indel in the CIGAR and 50 more a soft clip; I
         have traced the indel path only. The rules also credit the one-unit ALT to reads that carry two units (12 tdel,
         26 tins; tdel k=94, 62M2D64M at the tract, which the judge calls another allele). I have not traced that path.
      3. *Rules: the bases are read where the CIGAR put them* (the judge's 231 added credits that hold). Example, tins
         k=277 (97M1I28M, every base q30+): the read begins two bases before the tract and shows 14 A's for REF's 13.
         The aligner wrote the extra A as a mismatch at the flank base. The rules call REF; the judge calls ALT (−7
         against −11). This is RJ-33's known limit of the rules path, and the judge removes it.
      4. *Both arms: a masked base where the tract meets the run beside it* (the runs beside, above; the worked
         example below).

      Two traced reads follow decided rules and are not engine faults. tdel k=153: the aligner scores another place
      higher (AS 29, XS 32), but the read fits the candidates within 10 mismatches (best −18, above −40), so the
      decided `xs` fit rule leaves it to the margin (`judge/decide.rs:352`); it fits a two-unit deletion best
      (cause 1). tdel k=85 (100S26M): the judge reads the soft clip (RJ-13), past the check's window.
    - **Worked example, corrected** (tdel k=6: a 6-T run after a run of 10 C's; both arms credit both reads as the T
      deletion). The first version of this example came from the flawed bases view; with the qualities it reads:

      ```
      ref   CTCCCCCCCCCC|TTTTTT|GTTTTGATCTT
      read  CT-CCCCCCCCC|TTTTTT|GTTTTTATCTT   CIGAR 19M1D104M3S, mapq 60
      qual  33 333333333|123333|12333111112   (one digit: quality // 10)
      read  CT-CCCCCCCCC|TTTTTT|GTTTTGATCTT   CIGAR 19M1D107M,   mapq 60
      qual  33 313233332|133333|11231112233
      ```

      In both reads the first T after the C run, and the G after the T run, have quality 10-19, so the engine reads
      them as unknown. Masked, each read shows `C9 ? T5 ? T4`. The ALT (`C10 T5 G`), a one-C deletion (`C9 T6 G`) and
      a deletion of the G fit it equally (cost 0); REF fits only with a gap (cost 7). So the reads hold the ALT over
      REF and over the T run's other lengths, and tie the slips of the runs beside it. The rules move the read's
      deletion from the C run into the tract over the masked base, with no confident error, so the read is `Given`;
      the judge gives ALT 0 against REF −7.
    - **Harness slips found by the checks** (none changes an engine count): the report first fetched reads only from
      the tract start, so it missed reads that the engine admits by a soft clip over the anchor (their aligned bases
      end before it); and htslib writes a warning ("the index file is older than the data file") to the same stderr
      with no line end, so 11-21 trace lines per run began with it, and the parser, which matched at the line start,
      dropped them. `rna333h_rl_reads.py` parses the same way; the earlier read checks that import it need the same
      count. The bases view's flaw is the correction above; my oracle script also had a sign error (the judge's gap
      taken as REF − ALT), found and fixed before its numbers were used.
    - **Limits.** A sample chosen for AD's excess, not every row. One library. REF-by-truth rows only, so this check
      measures false ALT, not missed ALT. The masked check uses the judge's own scores, not an independent model; its
      candidates are REF, the run's lengths and (second verdict) the runs beside, and its window is 30 bases each side.
    - **Measurement switches in shipped code** (seen in the trace; on the "drop the dead switches" list):
      `GBCMS_JUDGE_CLIPSUB` (`judge/decide.rs:1423`), `GBCMS_JUDGE_NEITHER` (`judge/decide.rs:1010`),
      `GBCMS_JUDGE_REFBASE` (`judge/decide.rs:1041`, `judge/decide.rs:1073`), `GBCMS_JUDGE_FITS` and
      `GBCMS_JUDGE_FITS_MEASURE` (`judge/decide.rs:269-327`, marked temporary). The runs here clear every `GBCMS_*`
      variable, so these switches did not act.
  - **Next:** a decision on cause 1 (the judge applies RJ-27's window at short pure indels and at an SNV at a run's
    end, or takes the run's other lengths as candidates; since 2026-11-24 one decision with the DNA check's "either"
    reads, both views measured on DNA and RNA: `334-run-length-reads.md`), measured first on these rows (false ALT)
    and on the truth
    indels and SNVs at tracts (missed ALT), under the decision protocol; cause 2 likewise, unless the judge becomes the
    default first; the operator's view on the runs beside the tract (cause 4); trace the rules' two-unit credits and
    their soft-clip credits; the htslib-warning count in the earlier read checks' traces; the high tracts; E12's M1 on
    HG002 RNA; E12's M2 prototypes (a ticket and branch of their own, the operator's decision); the 630
    `d_amb`/`n_ref` row-samples.
  - **The germline het and hom check of MNVs and complex events** (2026-11-23, the operator: "the point to switch to
    judge was the hope to be better at mnvs and complex events"): its own record, `334-germline-check.md`. In short,
    the judge lowers missed ALT and false REF in every MNV and DELINS stratum. It raises false ALT on complex alleles
    one deleted base from the truth, the same near-miss credit class as cause 1.

## 8. Cost log

- Research: M12 (Opus) 0.37 M tokens; rules-path map (Opus) 0.48 M; survey (Fable) 0.19 M; cache audit (Opus)
  0.40 M. Review of revision 1 (Fable) 0.24 M. Gap-model survey for E12 (Fable, 2026-11-23): tokens not recorded.
  HPC (2026-11-23): the panel-row truth run (jobs 13284559, 13284561, 13284578), the HG002 DNA download (13284844),
  the HG002 DNA truth score (13288607 prep, 13289325 arms: CPU 12.9 min rules, 95.9 and 95.7 min judge), the read
  check at REF tract rows (13316156 traced runs; reports 13316186, 13318536, 13319409, 13320130, 13322149 failed their
  checks, 13322742 passed; the trace dump and the masked check: reports 13340383, 13342413, 13345884, 13346360,
  13346986, 2 CPUs each; the last two ran 41 s and 39 s).
