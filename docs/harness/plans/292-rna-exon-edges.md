# 292: reads at exon edges — does a splice lower a base's quality?

Status: measured; reviewed (Fable, "sound with conditions"); conditions met; decided 2026-11-02 (option B); merged 2026-11-03 in #415
Rounds: 1 of 2 · New conditions: 0 of 0 (option B removes one) · HPC batches: 5 of 6 (one final build left)

## 1. The class

- **Reads:** RNA reads spliced at a junction, judged at a base within five read bases of it: pre-mRNA
  (unspliced), spliced with a long anchor (5 or more aligned bases past the junction), or a short one.
- **Evidence:** the read's own base and its quality. Today heuristic BAQ lowers the five read bases each
  side of every junction by 20 (`rust/src/shared/baq.rs:21-23`, `:137-150`), so a base passes the run's
  `--min-baseq` only where its quality was 20 above the bar: Q40 data passes at bar 20, nothing at bar
  30. With a GTF, BAQ is off within 5 bp of an annotated edge (the exon-edge exception, R1 #106,
  `engine.rs:3530-3562`); 87% of the measured rows sit at annotated edges.
- **Invariants:** one quality gate on every path (2); ALT only by the read's own bases (7); binning
  (BAQ is per read); DNA reads have no junction.
- **Principle:** a base is judged by its own quality; only a mechanism that misplaces it (a piece too
  short for its bases to place it, RJ-3) justifies a mask.
- **Instance:** the #291 review: the penalty masked the bases that place a read's own indel (fuzz: 63
  wrong-allele calls, 53 under BAQ).

## 2. What already touches it

- **Rules (read to their ends):** RJ-17 (a splice is not coverage; a partial deletion piece is depth
  only, amended 2026-10-30), RJ-18 (a clip at an edge is the next exon), RJ-19 (the spliced reference),
  RJ-23 (a cut flank reads from the other), RJ-27 ("with BAQ masking the bases at its junction": B
  rewrites that clause), RJ-3 (a read is placed by a flank base it reads).
- **Rejected:** a read-level "no penalty at annotated junctions" (REJ-20260925-002); splice-aware Phase 3
  and hard-filtering junction-adjacent variants (REJ-20260924-001).
- **Tickets:** #106, #289/#291, #385 (queued), #293, #178; the ALT-side insertion gap (queued). Dropped
  2026-10-30: (c) GATK's overhang clip, (d) the short-anchor mask (section 3; the note records the
  decision, not the reason). Code conditions and contract tests: in the PR text.

## 3. Options

**Survey.**

| Tool | At a junction | Source |
|:--|:--|:--|
| GATK SplitNCigarReads | splits at N; clips an overhang only with more than 1 mismatch in at most 40 bases; no quality change, HaplotypeCaller reads full quality (B is nearer this than A) | tool doc: `--max-mismatches-in-overhang 1`, `--max-bases-in-overhang 40` |
| htslib BAQ (`sam_prob_realn`) | skips any read with a reference skip whole: `samtools` is the field's practice | htslib `sam.c`; Li 2011, Bioinformatics 27:1157 |
| ASEReadCounter | pileup with quality floors (defaults 0); no junction rule named; GATK runs SplitNCigarReads first | tool doc |
| phASER, WASP | quality floors; WASP re-maps for allele bias; no junction rule documented | READMEs; `read-filters.md:265-275` |
| RNA-MuTect | variant-level filters (realignment, editing, panel of normals); a splice-distance filter not confirmed (code page unreachable) | Yizhak 2019, Science 364:eaaw0726 |
| SNPiR | drops intronic candidates within 4 bp of a junction: the variant-level filter REJ-20260924-001 rejected | README, `filter_intron_near_splicejuncts.pl` |
| STAR | splices a novel junction from 5 bases past it (`alignSJoverhangMin`), an annotated one from 3 (`alignSJDBoverhangMin`); FORTE runs `--alignSJDBoverhangMin 1`; `outSJfilterOverhangMin 30 12 12 12` gates the SJ file only | `parametersDefault`; `msk-assay-facts` |
| HISAT2 | no anchor minimum documented | manual |

**The runs:** 58 FORTE slices (STAR), SNV rows at the first ten bases of exons beside sampled junctions,
163,882 rows in reach (bases 0-4) and 168,243 as the control (5-9); every rule, GTF and BAQ on and off,
bar 20 and 30; an oracle with no gbcms rule classes rows (clean, low, het by where the ALT sits); the
DNA check labels rows by the matched IMPACT DNA (40 pairs); the synthetic set, 127,104 reads with truth.
Numbers: no GTF, BAQ on, bar 20 unless said; with a GTF every option already counts as B at annotated edges.

| Option | The rule | Inv. 7 | One gate, both paths | Conditions | Will not handle |
|:--|:--|:--|:--|:--|:--|
| **A. Keep the penalty** (`penalty`) | A splice lowers the five read bases each side by 20; the exception lifts it at annotated edges with a GTF. | passes | fails in effect: strength is the quality bin and the bar | 0 (keeps 3) | its target class: the ALT sits in long anchors |
| **B. No junction penalty on any path** (`nopenalty`; `samtools` differs by 7 ALT, 1,476 REF reads) | A splice never lowers a base's quality; the read's indels keep the indel penalty. | passes | passes | −1 | a short piece STAR writes at an annotated acceptor when the molecule spliced at a novel one 2-6 bases on (RJ-3's class) |
| **C. Spare the exact-carrier rule only** (`spare`) | That rule judges without the penalty; every other path keeps it. | passes | fails: two gates for one base; 0 moves on every SNV row | +1 | SNVs and every read that rule does not judge |
| **Bar 10 and 0** (synthetic sweep, 2 seeds; no GTF) | Bar 10: A and B are near (het VAF error 0.0474 and 0.0467; novel-acceptor ALT 2.1% and 2.2%). Bar 0: B equals BAQ off, row for row; A differs at 8 of 1,072 rows (12 of 14,667 ALT reads). | – | – | – | Clean-row ALT at bar 0: 104 (A), 105 (B), from the evidence floor of 0 (`carrier.rs:1170`), not the penalty. B's extra read is one Q2 or Q11 base |
| **Indel and delins rows at exon edges** (58 FORTE slices; 41,112 rows in reach; no GTF) | B against A in reach, ALT and REF reads: bar 20, −3 and +31,014; bar 30, +1,432 (+3.3%) and +1.89M (+4.5%); bar 10, +179 and +38,751; bar 0, +259 and −220. | – | – | – | Clean rows that gain ALT under B: 1, 5, 45, 96 (bars 20, 30, 10, 0); with a GTF, at most 1. No germline indel in the DNA check |

**B against A, in reach** (the control band's reads move by under 0.1%):
- 136,257 of 163,882 rows move: ALT +14,269 reads, REF +5.70M (+1.9%). Clean rows (87,441): ALT 0
  under every rule. Low rows: the ALT share in reach is 0.155% under B, the control band's 0.158%: no
  excess of mismatches in an exon's first five bases. Hets (717): ALT +1,235, REF +1,552; the VAF gap
  to the unspliced reads 0.0718 (A), 0.0721 (B).
- **DNA truth.** At the 219 germline hets, A hides 984 ALT reads (1.9%) at bar 20 and 38,541 (76%) at
  bar 30. At the 79 DNA-absent rows whose ALT is in spliced reads only, A hides 52 of 2,248 (2.3%); in
  the oracle's 137-row class of that shape (58 slices), 4,309 of 4,873 ALT reads sit in long anchors,
  out of any junction rule's reach. At 100,278 DNA-absent low rows B adds 16,243 ALT reads (3.9%, the
  background rate) and 5.03M REF.
- **Synthetic.** Het VAF error 0.0610 (A), 0.0475 (B), 0.0591 (C); A leaves 7,325 SNV reads depth only.
  Wrong calls 102 (A), 184 (B); the exact-carrier rule's fall 52 to 33 (indels), 50 to 18 (delins).
  False ALT on hom-REF rows: 0 under every rule at the plain and overhang strata; at the novel-acceptor
  stratum (hom-REF rows too) 133 of 6,097 reads (2.2%) under B, 0 under A, which hides half the depth
  there (2,967 of 6,097).
- **(c) and (d):** (c) loses 8,515 ALT reads at hets against 899 gained and misses 128 of the 133;
  (d) at N = 5 closes them (1 of 5,965) by masking 1.98M REF reads A keeps (`m292/syn/score_p292a.txt`,
  `p292a/agg/moves_by_class.tsv`).

**Indel and delins rows, B against A, in reach** (85,135 rows, 88% at annotated edges; in reach: clean
37,297, low 3,748, het 67):
- **Bar 20.** 3,699 of 41,112 rows move; ALT up 56, down 59. Clean rows hold ALT under A too (80 reads, 46
  rows); under B 76 reads, 43 rows: one row gains 2 reads, four lose. Low rows: ALT share 0.124% under
  both rules, the control 0.108%. Hets: ALT 3,167 to 3,158; VAF gap to the unspliced reads 0.0730 (A),
  0.0728 (B).
- **Bar 30.** 23,065 rows move. The control moves too (REF +0.9%; the SNV control, under 0.1%): cause not
  tested. Delins go the other way: ALT 280 to 173.
- **Bars 10 and 0.** The clean-row gains are insertions (42 of 45, 96 of 96); the control gains 5 rows. At
  bar 0 the clean rows hold 146,917 ALT reads under A, 99.9% at donors: the synthetic bar-0 shape, not the
  penalty. The oracle counts Q20 bases only.
- **DNA check.** 196 rows (34 slices) with ALT in spliced reads only: 173 DNA-absent, 23 unknown, 0
  germline, 0 somatic. At the 50 DNA-absent rows in reach, A and B show 235 ALT reads at bar 20 (A hides
  0); at bar 30 A hides 7 of 234 (3.0%).
- **Tasks.** None missing: 58 of 58 complete at 4 bars, 0 failed (tasks 14 and 30 timed out; one rerun).
  18 slices have no matched DNA; 6 of the 40 have no row to ask about.

**Worked example** (per-row means from the DNA check; no position): a germline het SNV at an acceptor's
third base (S).

```
  intron ....................][exon 2:  . . S . . . . .      S = the row's base
  pre-mRNA        ==========================S=========      unspliced: S at its own quality
  spliced, long   =====N . . . . . . . . N[= = S = = = = =   12 bases past the junction; S is 2 bases in
  spliced, short  =========N . . . . . . N[= = S              3 bases past the junction: a piece of 3
```

| Row class (rows) | A, bar 20: REF, ALT (VAF) | B, bar 20 | A, bar 30 | B, bar 30 |
|:--|:--|:--|:--|:--|
| germline het, ALT in unspliced reads too (219) | 304, 232 (0.433) | 311, 237 (0.432) | 68, 55 (0.448) | 302, 231 (0.433) |
| DNA-absent, ALT in spliced reads only (79) | 42, 28 (0.396) | 43, 28 (0.401) | 37, 0.7 (0.019) | 41, 28 (0.401) |

At bar 20 A and B agree within 2% and A guards nothing at the false class; at bar 30 A hides three
quarters of a true het's reads: a splice filter whose strength is the quality bin.

**Recommendation: B.** The penalty removes 2% of the reads in reach at MSK's bar and no false ALT the
DNA shows; its strength comes from the instrument, not a mechanism; B is nearer GATK's practice than A.
**The limit, sized:** the synthetic set models STAR at overhang 3 and 5 (`syn292.py:155`) while FORTE
runs `--alignSJDBoverhangMin 1`, so 2.2% is a floor; but the limit exists today at bar 20 under A too
(52 of 2,248 hidden), so B adds no new limit at MSK's bar. **The residual, named:** RJ-3's class, a
piece too short for its own bases to place it is depth only (`shortanchor` with its cutoff from STAR's
novel-junction minimum, not a fixed N), to be ticketed. Why (d) was dropped is not on record; the
measured reason: its cutoff is the aligner's setting (FORTE 1 at annotated junctions, 5 at novel) and
at N = 5 it masks 1.98M REF reads A keeps, so it needs its own sweep.

## 4. Measurement (done, with what is missing)

- **Metric and bar, per class:** clean rows keep ALT 0; the low rows' ALT share in reach within the
  control band's spread; DNA-germline hets lose no ALT read the base's own quality admits; false ALT
  added at DNA-absent rows, by class, against the background rate; synthetic false ALT by stratum.
- **Data and truth:** RNA only: 58 FORTE STAR slices at bar 20 and 30 (40 with matched IMPACT DNA); the
  synthetic set (overhang 3 and 5; true and shifted GTF; the fuzz is stratum f). DNA arms are unchanged
  by construction: one DNA identity batch is owed before merge.
- **Done 2026-11-03:** the bar at 10 and 0 (local synthetic sweep); indel and delins rows at exon edges
  on the FORTE slices, A and B with the DNA check (HPC).
- **Still missing:** other RNA aligners (HISAT2; STAR at FORTE's overhang 1 in the synthetic set; long
  reads); the DNA check for (c) and (d); a second RNA cohort; strandedness at one setting; a germline
  indel for the DNA check.
- **Cannot see:** whether a short-anchor ALT read at a DNA-absent row is a novel-acceptor piece or a
  pseudogene read; whether the ALT reads B adds at clean insertion rows (bars 10, 0) are base errors.
  To build: #385's shifted-junction census; REF alternative acceptors at 5-30% usage in the synthetic
  set.

## 5. Budget and stop rule

Two review rounds (1 used). New conditions: 0; a fix that adds a junction-keyed mask or a path special
case stops the work and brings the class question back with one option. Builds: 1. HPC: the indel
batch (done; its one rerun is used), 1 final-build batch.

## 6. Work and models

- **Sonnet:** remove BAQ's splice arm, the switch and its lane; the 30 measurement tests become contract
  tests of B, red-first; a register row ("a splice never lowers a base's quality") amending RJ-27; the
  effects map in the PR. The code, test and doc list (14 files from `grep SplitNCigarReads`, the CLI
  help, the Nextflow lines): `292-rna-exon-edges.builder.md`. Nextflow flags unchanged.
- **Main session:** the HPC batches after approval. **Fable or Opus:** the second review, if any.

## 7. State and cost (tail)

2026-10-17 measure decided; 10-30 four rules built (e45a3bb3), the partial-span rule committed
(b247db31), (c) and (d) dropped; 10-31 arm aggregated (58 of 58); 11-02 DNA check aggregated (40
pairs), page drafted, reviewed (round 1), revised; 11-03 bar-10 and 0 sweep recomputed, indel arm and
DNA check aggregated (58 of 58 tasks; 40 DNA pairs), both rows filled. HPC before 11-03: 2 builds, 4
batches, 2 aggregations, 2 task reruns (jobs 5430389/90, 5438438/39, 5441284/99, 5513308/11, 5561674,
6012647). The indel batch added 1 batch, 2 aggregations, 1 rerun (jobs 6074312/406/413/447; tasks 14 and
30: 6158990 with 6159043/61/91). Agents: Opus (prototypes, rule, scorers), Fable (page, review), Sonnet
(11-03 rows); tokens in their usage reports.

**Decided 2026-11-03 (operator): B, no junction penalty on any path.** Next: the build (this page's
builder list; red-first; register and docs); before the PR, explain why delins ALT falls under B at
bar 30 (280 to 173). Ticketed: the bar-0 evidence floor (#408), RJ-3's short anchors at a novel
acceptor (#409). Still open: #385 and the insertion gap behind it.
