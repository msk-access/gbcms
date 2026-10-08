# 6.7.0 cycle — plan

> The 6.7.0 work: tracker **#196** and the **6.7.0 milestone**. Sources: the 6.6.0
> release gate's follow-ups (#253, #254, #255, adjudicated on #252), the items deferred
> at the 6.6.0 cut (#244, #246's sibling shapes, the panel's GIAB and TEMPO arms, the
> `build-gtf-cache` removal), and the issues the 6.6.0 triage moved here. Branch-from:
> `develop` at `6.7.0-dev.0`.
>
> **Discipline (unchanged from 6.6.0).**
> - Measure first; red-first tests (strict xfail before the fix); an adversarial review
>   before merge; real-data acceptance per read on the data-type matrix (below).
> - Read judgment is decided spec-first: `docs/reference/read-judgment.md` (RJ register)
>   changes first, the operator decides, code follows.
> - Counts must not depend on bin geometry; classification is checked against the read
>   census, and read-input rules against an oracle the census does not share.
> - No new output columns: new signal goes to `gbcms_diagnostic`, logs or validation
>   tooling. Aggregates and allele shapes only in issues, PRs and this plan.
> - Check `.agents/learnings/REJECTED.md` before proposing an alternative.

**Status (2026-10-08).** Draft for operator review; nothing is implemented. The
decision tables below need the operator before any ticket starts.

Priority: **H** high, **M** medium, **L** low. **[counts]** can change counts.
**[decide]** needs an operator decision first. **[new]** has no issue yet.

## Summary

| ID | Ticket | Pri | Flags | Issue | Parent (proposed) |
|:--|:--|:-:|:--|:--|:--|
| C39 | Long delins: read every complex event by one-sided, equal-length junction windows (no padding loss, no REF double count) | H | [counts] [decide] | #253 | A |
| C41 | REF and VAF of indels in repeats against germline truth (allele-level vs by-base, spanning bias) | H | [decide] [new] | from #253 | A |
| C40 | Long deletions: one window at the POS breakpoint for both alleles, clipped bases read (RJ-13) | H | [counts] [decide] | #254 | B |
| C42 | Per-locus gap parameters (stutter model) | L | [decide] | #244 | A |
| C38 | Sibling rows: SNV sibling double-counts a compensating-mismatch carrier | L | [counts] [decide] | #246 | A |
| C33 | Exact-carrier REF and ALT windows read from different anchors | L | [counts] | #214 | A (folds into C39) |
| C6 | Error-tolerant exact-length insertion matching (identity band) | L | [counts] | #143 | A |
| C24 | Local-alignment fallback tail reads stale scores | L | [counts] | #195 | A |
| R3 | RNA: catalogued editing positions inside carrier windows | L | [counts] | #178 | A |
| C15 | C12 follow-ups: clipped carriers of pure deletions, anchors in the clip | M | [counts] | #173 | B |
| C18 | Split-read evidence for long events (supplementary alignments) | M | [counts] | #177 | B |
| C7 | Clip-borne ITD carriers | L | [counts] | #144 | B |
| C31 | Fragment end from the mate's unclipped 5' end (MC tag) | L | [counts] | #212 | B |
| O5 | Mapping-bias diagnostic (ALT reads mapped or clipped worse than REF) | M | | #179 | B |
| D8 | `attribute.py prepare/check/attribute` skip runs whose latest attempt failed | M | | #255 | C |
| D9 | Panel arm: GIAB HG002 truth (public, PHI-free) | H | [decide] [new] | — | C |
| D10 | Panel arm: WES (BWA only) at the panel's paired long-event loci | H | [decide] [new] | — | C |
| D11 | Panel arm: IMPACT germline hets in matched normals | M | [new] | — | C |
| C3 | Homopolymer decomposition: fix or retire `--rescue-homopolymer` | M | [counts] [decide] | #111 | D |
| — | Decomposed-allele wins (umbrella; #145, #146) | M | | #112 | D |
| M1 | Merge rows whose flavors report different alleles | M | | #128 | D |
| O10 | End-of-run QC summary | M | | #227 | E |
| O3 | MNP rescue in fillouts: flag a germline component | L | | #132 | E |
| O6 | Read-orientation evidence for oxoG/FFPE (decide) | L | [decide] | #180 | E |
| I5 | Nextflow `GBCMS_CONVERT` module | L | | #127 | E |
| H4 | Remove `build-gtf-cache` and `--gtf-cache-dir` (deprecated in 6.6.0) | M | [decide] [new] | — | F |
| D3 | mkdocs-material 2.0 | L | | #138 | F |
| P1 | Deep-bin fetch reduction | L | | #150 | F |
| C34 | GTF transcripts keyed by ID alone | L | | #216 | F |
| I7 | MAF dash insertion at `Start_Position` 0 | L | | #219 | F |
| D7 | Platforms: RHEL first, macOS second (D7a–D7e) | H | | #232 (#233–#237) | G |

## What was measured for this plan (2026-10-08)

Local harness `~/test/gbcms/harness/plan670/` (PHI stays there; aggregates here).
- **Data.** The 6.6.0 gate panel (`tier_660`): every delins row (137) and every
  pure deletion of 50 bp or more (60) in its IMPACT tumour runs, 133 BAMs; and 195
  germline-tagged insertions and deletions from the sign-out dump, genotyped in the
  patients' IMPACT normals (one per patient).
- **Engine.** develop (the 6.6.0 counting code) with per-read traces, on local slices
  (one serial pass per BAM). The slices reproduce the HPC gate's full-BAM output: all
  130 rows without a neighbouring variant are identical; rows are run with their
  neighbours (every panel variant within 300 bp), as the gate ran them.
- **Oracle.** Each read's own bases, soft clips included, placed ungapped on each
  allele's haplotype at the left and the right anchor; a base below Q20 or an N
  matches anything. It reports the flank a read holds on each side of the event and
  where its confident mismatches fall, independent of how the aligner wrote it.
- **Surveys** (sources in each table; full reports kept with the harness): complex
  anchors and spanning bias (18 tools plus Varlociraptor and biastools), clipped
  junction reads (SV callers and genotypers), stutter and gap models (21 tools), GIAB
  data sources.

## Decisions

### 1. C39 #253 — how a read must hold a long delins

**Today (RJ-4, RJ-9).** A complex event is read as whole windows: the event (grown
through repeats) with two flank bases, the shorter allele's window padded with
reference flank until both are equal, the padding split evenly between the two sides.
So a 28>6 delins needs about 13 clean flank bases on *each* side of the ALT. Past 50
bases, windows become junction windows read from either flank.

**Measured** (137 delins rows; "clean" = no recurrent neighbouring change, below):

| Rows | Rule | AD | RD | VAF ÷ 6.6.0 (median) | AD within 10% of sign-out |
|:--|:--|--:|--:|--:|--:|
| 14 clean, span 20+ (whole windows today) | 6.6.0 | 1,203 | 4,117 | 1.00 | 0/14 |
| | whole allele, 2 flank bases, no padding (what the sign-out counts track) | 1,712 | 4,569 | **1.22** | 9/14 |
| | one-sided junction windows, read from the left flank | 1,742 | 5,804 | **0.98** | 9/14 |
| | one-sided junction windows, read from the right flank | 1,712 | 5,770 | **1.00** | 9/14 |
| | junction windows from either flank (today's rule past 50 bases) | 1,774 | 6,984 | **0.87** | 9/14 |
| 4 clean, window over 50 (junctions today) | 6.6.0 | 458 | 2,403 | 1.00 | 2/4 |
| | one-sided, left / right | 459 / 460 | 1,351 / 1,350 | **1.40 / 1.56** | 3/4, 2/4 |
| 70 clean, span under 20 | 6.6.0 | 9,557 | 24,053 | 1.00 | 52/70 |
| | one-sided, right | 10,227 | 26,076 | 1.00 | 68/70 |
| 49 with a recurrent neighbour | 6.6.0 | 2,401 | 14,964 | 1.00 | 19/49 |
| | whole allele, 2 flank bases | 6,853 | 15,125 | 1.04 | 23/49 |
| | one-sided + observed-flank guard | 1,030–1,100 | 13,500–13,640 | 0.97–0.99 | 8–12/49 |

- **Padding costs reads, not bias.** On long delins the padded windows drop 31% of the
  ALT carriers and 29% of the REF reads that one-sided windows judge, and the VAF is
  the same (left 0.98, right 1.00). The carriers lost are the ones the gate saw:
  sign-out agreement 0/14 → 9/14. (On short delins the oracle's left reading under-reads
  9 of 70 rows, a limit of its ungapped placement where the aligner shifted the anchor;
  the build's own left and right readings are compared there.)
- **The two biased alternatives behave as predicted.** Whole alleles without padding
  favour the shorter allele (VAF ×1.22). Reading junctions from either flank favours
  the longer one (×0.87), because a REF read can prove REF at either end while the
  shorter ALT has one window.
- **6.6.0's long-event rule has that second bias today.** On the 4 clean delins with a
  window over 50, it counts REF at both junctions: VAF reads 0.64–0.71 of the
  one-sided value.
- **Recurrent neighbours.** In 49 rows the well-anchored ALT carriers share a confident
  change at a fixed distance (up to 25 bp) from the event. In 47 it is on ALT reads
  only, so the reads carry a larger allele than the one given. In 2 it is also on REF
  reads (likely a germline SNP in cis). 6.6.0's padding withholds part of these; a
  minimal-flank rule would credit the larger allele (6,853 vs 2,401 AD).

**Community practice** (s1 survey):

| Tool | Complex/delins rule | Source |
|:--|:--|:--|
| GATK HC / Mutect2 | read scored against whole haplotypes; counts in AD when the best allele beats the next by log10 0.2; no anchor rule, no padding | `AlleleLikelihoods.java`, `DepthPerAlleleBySample.java` |
| Strelka2 | ≥5 bases on each side of *one* breakpoint; posterior ≥0.51 | `starling_read_align_score_indels.cpp` |
| Octopus | reads assigned to the MAP haplotype; ties dropped | `read_assigner.cpp` |
| freebayes | full observations (RO/AO) must span the whole haplotype window; partial support reported apart (PRO/PAO) | `AlleleParser.cpp`, README |
| Paragraph | each breakpoint genotyped independently, ≥16 bases anchoring each side | `doc/graph-counting.md` |
| GangSTR, ExpansionHunter, Varlociraptor | model the chance a read spans each allele (∝ read length − allele window) | Mousavi 2019; `FragLogliks.cpp`; Köster 2020 |
| bam-readcount, LoFreq, VarDict, DeepVariant, Platypus | by base; no complex allele, or CIGAR-based | per-tool source |
| GetBaseCounts (C++) | delins classed as INS/DEL by length; deleted sequence never checked | `GetBaseCountsMultiSample.cpp` |

No surveyed tool pads windows to equal length. Likelihood callers need only the
discriminating bases at one breakpoint; three tools correct spanning bias with a model.

**Options.**

| | Rule | AD (long delins) | VAF | Cost / risk |
|:--|:--|:--|:--|:--|
| A | Keep 6.6.0 | baseline | unbiased under 50 bases; REF double-counted past 50 | loses ~30% of carriers; long-event VAF ×0.64–0.71 |
| B | Whole allele, minimal flank, no padding | +42% | ×1.22 (shorter allele favoured) | biased; credits larger alleles at recurrent loci |
| C | **One-sided, equal-length junction windows for every complex event**, plus an observed-flank guard | +45% (RD +41%) | ×0.98–1.00 | needs a per-variant pass for the guard |
| D | B plus a spanning-corrected VAF (GangSTR-style weights) | +42% | corrected | a new column; counts and VAF disagree |

- **C in detail.** Read every complex event from one flank, inward. Take the shorter
  allele's whole window (2 flank bases each side) and the same number of the longer
  allele's bases. Then every further base the read holds must fit that allele (RJ-9's
  read-on, unchanged). Both alleles are judged at the same anchor over windows of the
  same length, so neither is favoured by where reads start. The same rule then serves
  short and long events (no 50-base switch).
- **The guard.** Where the locus's well-anchored carriers show a recurrent change on ALT
  reads only, both windows grow (equally) past it, so a read must show it is the given
  allele and not the larger one (invariant 7). A recurrent change also present on REF
  reads (a germline SNP) is masked instead, not grown past. The diagnostic names the
  larger allele (`OBSERVED_ALLELE`, unchanged).
- **The anchor flank.** The POS (left) flank, the breakpoint pure deletions count REF at
  today. A read whose own clip holds that flank is read there (RJ-13). Reads holding
  only the right junction become depth for both alleles alike.

**Recommendation: C.** It departs from practice (no tool reads equal one-sided
windows), but it measures better than each alternative. It keeps the unbiased VAF the
padding was there for, reaches the carriers likelihood callers and the sign-out count,
and removes today's long-event REF double count. D needs a column and lets counts and
VAF disagree. The guard is the one sub-decision:
- **(a) guard on:** strict, invariant 7. AD on the 47 larger-allele rows falls from
  2,401 to about 1,065.
- **(b) guard off:** like 6.6.0's padding, which guards partly by accident.

I recommend (a).
**Spec:** a new RJ-23 amending RJ-4 and RJ-9. C33 #214 (REF and ALT windows read from
different anchors) is closed by construction: one anchor for both.
**Effects map.**
- `carrier.rs`: `windows()`, `Event::long` and `room_needed` go; the guard pass
  is new.
- Both read loops use `carrier::classify` (the main loop and per-transcript
  counting), so they change together. Spliced windows (`cut_at`, RJ-19) are rebuilt
  one-sided. ASJD, mFSD and observations follow the per-read calls.
- MNP rescue and clusters: unchanged (MNPs keep their path).
- The census (`tests/census.py`), the spec test and the clip-carrier contracts change
  in step. Existing expectations that move go to the operator first.
- Validation: the RC set, WES, FORTE and the GIAB arm. The two one-sided readings must
  agree (an unbiasedness check), and germline delins hets must read 0.5.

### 2. C41 (from #253) — REF and VAF of indels in repeats against truth

**Today.** REF counts by allele: a read is REF only when it reads past the tract's
first difference from its flank (RJ-1), and `vaf` = AD/(RD+AD). Pileup tools and the
sign-out count REF by base at the indel's first position.

**Measured (pilot):** germline hets in IMPACT normals (truth 0.5), engine `vaf` vs a
by-base VAF (REF = reference base at the base after the anchor, as GetBaseCounts and
pileup tools count it) vs a crossing census (reads crossing the repeat on their
allele).

| Germline hets (193) | n | `vaf` median (IQR) | by-base VAF | crossing census | informative share |
|:--|--:|:--|:--|:--|:--|
| Deletions in repeats | 100 | 0.499 (0.465–0.518) | 0.444–0.473 | 0.504–0.513 | 0.88–0.93 |
| Deletions, unique | 28 | 0.477 (0.456–0.516) | 0.479 | 0.483 | 0.97 |
| Insertions in repeats | 35 | 0.468 (0.449–0.500) | 0.440–0.454 | not measured | 0.85–0.92 |
| Insertions, unique | 30 | 0.481 (0.455–0.505) | 0.462 | not measured | 0.95 |

- **Repeats.** Allele-level VAF holds 0.5 for deletions in repeats, whatever the unit
  (homopolymer 0.502, 2-bp 0.501, 3-bp 0.490, 4-bp 0.498) or tract (under 10 bp 0.500,
  10–14 bp 0.496).
- **By base.** VAF reads 3–6 points low, and falls further with unit size (4-bp STR
  0.444) and tract length (10–14 bp 0.451). That's the reads ending inside the tract,
  counted REF.
- **The small shortfall in unique sequence and for insertions** (0.47–0.48) shows in both
  counts. It is mapping bias against ALT reads (Lunter 2011), not repeat handling; O5
  #179 is where that becomes a diagnostic.
- **Pilot limits.** Tracts are at most 14 bp (the reported germline findings are short
  frameshifts), and the data is ABRA2-realigned panel DNA only. Long tracts,
  trinucleotide repeats, longer indels and WES/WGS need GIAB (D9). The worked example
  is the origin event in #253: a 9 bp deletion in a 51 bp repeat reads 15/31 by allele
  against 6.6% by base.

**Community practice** (s1): GATK, Strelka2, Octopus, bcftools (`--ambig-reads drop`,
default) and freebayes (RO/AO) leave reads that end inside the tract out of AD.
bam-readcount, LoFreq, VarDict, DeepVariant's AlleleCounter, Platypus and the C++
GetBaseCounts count them REF by base. Only GangSTR, ExpansionHunter and Varlociraptor
model spanning bias. Literature: biastools (Lin 2024) shows the shorter allele
over-represented, growing with length; lobSTR alleles to ~45 bp show minimal bias
(Willems 2014).

**Options.** (a) keep allele-level REF and VAF, and document how they differ from IGV
and by-base counts; (b) a spanning-length correction, which needs a column or
non-integer counts (out); (c) a diagnostic when few reads can be judged (e.g. under
half the depth informative), through `gbcms_diagnostic`.

**Recommendation: (a).** The pilot supports allele-level REF and VAF. Document
how they differ from IGV and by-base counts, using the origin event as the worked example
and naming per-sample `total_count` (DP) against merged `simplex_duplex_total_count`
(ref + alt). Decide (c) on GIAB: add a diagnostic only if long tracts there show VAF
moving with the informative share.
The full measurement is the GIAB arm (D9): HG002 v5.0q hets with the repeat motif tag,
by unit (homopolymer, 2–6 bp, trinucleotide) and tract length. It needs your OK to
fetch region slices (decision 7).

### 3. C40 #254 (+ C15 #173 item 2, C18 #177) — long deletions

**Today.** A pure deletion's ALT is credited from the CIGAR's D op (RJ-2, RJ-6, the
≥50 bp band). Clipped carriers are not read. REF counts at the POS breakpoint with one
margin base. RJ-13 (adopted 2026-10-02, built in 6.7.0) says clipped bases inside the
fragment are judged by the same rules, and split reads join their molecule once.

**Measured** (60 deletions of 50+ bp, IMPACT, ABRA2-realigned):
- Junction carriers (the oracle, 4+ clean bases each side past the microhomology):
  6,137 written as a D op, 217 as a soft clip (3.4%; 87 with 10+ junction bases).
  Clipped carriers are in 28 of 60 rows; 6.6.0 credits none.
- Soft clips that start at the POS breakpoint (169 that can be told): 142 hold the
  junction (117 with 10+ clean bases past the microhomology, 12 with 8–9, 13 with
  4–7), 23 are neither allele (adapter, chimera, another allele), 4 are REF.
- REF: no 6.6.0 REF read at 100+ bp deletions holds only the far breakpoint (8,766
  hold the POS breakpoint cleanly; the other 1,901 show a confident mismatch there), so
  REF is not pooled across both breakpoints.
- VAF: windows of equal length at the POS breakpoint for both alleles, whatever the
  representation, agree with each other for k = 4–16 clean bases (VAF ÷ 6.6.0
  1.04–1.08). 6.6.0 reads 4–7% low: ALT reads with a short far side get clipped
  and go uncredited, while REF reads need one or two deleted bases.
- Reads judged under equal windows, against 6.6.0's AD:

  | Length | k = 4 | k = 8 | k = 10 |
  |:--|--:|--:|--:|
  | 50–99 bp | −3% | −30% | −33% |
  | 100+ bp | −15% | −26% | −36% |

  Crediting only the clipped carriers at k = 10 adds 1.2–1.4% to AD.
- WES (BWA only) is not measured here. C18's record shows 0 ALT on WES at 68–106 bp
  deletions the realigned panel counts, and BWA-MEM writes a D only when the far side
  outscores the gap (s2). The paired WES loci are decision 7's ask.

**Community practice** (s2):

| Tool | Clip use | Anchor past the junction | REF at the breakpoints | Source |
|:--|:--|:--|:--|:--|
| Manta | realigns the clip to the junction | 16 each side outside microhomology (8 in tier 2) | pooled from both breakpoints | Manta docs, source |
| Delly | split-read realignment | 13 each side | pooled (`RR`) | source |
| Paragraph | graph alignment | ⌊L/10⌋+1 (16 at 150 bp) | per breakpoint, genotyped apart | `doc/graph-counting.md` |
| SVTyper | where the clip falls | — | pooled (`RS`) | source |
| SvABA | realigns; score must beat the original | 5–8 | — | 2018 paper, v2.5 source |
| Pindel / MAVIS / STAR | pattern growth / junction match | 8 / 6 exact / 5 (3 annotated) | — | docs |
| VarDict | soft-clip consensus realigned | ≥7 for large deletions | by base | `VariationRealigner.java` |
| GRIDSS, CREST, LUMPY, Socrates | place the clip genome-wide | 20–25 | — | docs |
| GetBaseCounts (C++) | none | exact D only | one base at anchor+1 | source |

Tools that test a known junction use 5–16 bases. The 17–25 of the genome-wide placers
is a uniqueness limit a genotyper does not need. Pooling REF from both breakpoints puts a
het near 1/3 (the survey's reading of Manta, Delly and SVTyper).

**Options.**

| | Rule | IMPACT AD | VAF | WES |
|:--|:--|:--|:--|:--|
| A | Keep | baseline | 4–7% low | ≈0 ALT for deletions longer than about half a read (C18: 0 at 68–106 bp) |
| B | Credit clipped junction carriers at k = 10, REF unchanged | +1.2–1.4% | still low | most carriers |
| C | **One window at the POS breakpoint for both alleles, any representation, k = 4 clean bases each side past the microhomology; a low-complexity far flank raises k to 10 for both** | −3% (50–99), −15% (100+) | unbiased | all carriers that hold the window |
| D | C at k = 8 | −26 to −30% | unbiased | as C |

**Recommendation: C**, built with C18's split reads (one molecule, once) and C15's
anchors in the clip.
- It is RJ-13 as written: the same rule for the bases however the aligner wrote them,
  and the same window for REF and ALT.
- k = 4 is enough on this data: 23 junk clips at the breakpoint predict under 0.1 false
  carriers at 4^−4.
- The complexity guard covers flanks where chance matches are not 4^−k.
- B is the field's floor, but it keeps the REF-favouring bias.
- Scope: deletions of 50+ bp (RJ-6's regime). The 20–49 bp range is measured on WES
  first.

**Spec:** a new RJ-24 amending RJ-2 and RJ-6 for 50+ bp deletions.

**Effects map.**
- `check_deletion` and the windowed scan: the fetch must hold END ± read length for
  right-anchored clips, which the binning contract adds per variant.
- Fragment consensus: one molecule, once (SA).
- The census changes in step; validate against the bases, not the census (it shares
  read inputs).
- Validation on RC, WES (paired), GIAB large deletions, and at `--min-mapq 0`.

### 4. C42 #244 — per-locus gap parameters

**Measured.** Phase 3 (the PairHMM) decided 0 of ~421,000 read calls at the 197 delins
and long-deletion rows and their neighbours. From the C36 traces at insertion rows: 50
of 94,530 on RC (0.05%), 27 of 9,349 on FORTE DNA and 103 of 31,992 on FORTE RNA
(0.3%). At the 193 germline indels (mostly in repeats): 0 of 100,598. Where Phase 3
still decides, its calls disagree with the reads' bases: on FORTE RNA, 49 of the 103
reads it decided were called REF while their bases carry another allele.

**Community practice** (s3):
- **The field's floor:** repeat context indexed by unit and copy number, not by bp
  (GATK's PCR indel model, DRAGstr, Octopus, Strelka2). gbcms's bp-only curve gives
  A×10 and (CA)×5 the same gap open, where DRAGstr's table separates them by about 13
  Phred.
- **Fitted per sample:** DRAGstr (grid-search ML by period, thin bins pooled) and
  Strelka2 (a mixture fitted at two anchors, log-linear between).
- **Per-locus EM:** only in STR genotypers (HipSTR, at least 100 reads).
- **Stutter is deletion-dominated** (3–14×), and duplex consensus suppresses it 10–100×.

**Options.**
- **(a)** Close #244 as measured: no gap model can move more than the ≤0.3% of calls
  Phase 3 makes.
- **(b)** A static unit-aware table (DRAGstr's defaults) as a cheap floor.
- **(c)** A per-sample fitted table.
- **(d)** Per-locus EM with shrinkage, as filed.

**Recommendation: (a)**, with one replacement sub-issue: measure the reads Phase 3
still decides on every arm (RC, WES, FORTE, GIAB) and route the "another allele by its
bases" ones by their bases (neither, with partial evidence), as RJ-15, RJ-20 and RJ-21
did for their shapes ("never Phase 3's closer haplotype"). Revisit (b) only if an arm
shows Phase 3 deciding more than 1% of a row's reads.

### 5. C38 #246 — sibling rows

**Measured in 6.6.0 (2026-10-07):** 0 reads of the filed shape on RC. FORTE has no SNV
rows and WES no sibling rows. A complex row and its component deletion double-count
348 reads in 4 rows (kept and documented). The panel holds 616 SNV-within-5-bp-of-indel
rows (44 selected) and 35 MNP-near-indel rows.

**Recommendation.** Re-measure the double credit on the panel's SNV-near-indel and
complex-overlaps-component rows with the C39 harness, and close C38 as documented
unless a read is credited at both rows. If one is, the candidate rule stands: contest
an SNV-row ALT whose read is exactly an indel sibling's allele (C2's sibling rules,
RJ-5, extended to ALT).

### 6. H4 — remove `build-gtf-cache` and `--gtf-cache-dir`

Deprecated in 6.6.0 (a no-op with a warning, announced for removal in 6.7.0). Remove
the command, the `rna` option, the Nextflow `--gtf_cache` parameter, their docs and
tests, with a CHANGELOG breaking-change note.

**Effects map.**
- `src/gbcms/cli.py` and `models/core.py`.
- `nextflow/main.nf`, `nextflow.config` and `nextflow/README.md`.
- Docs: `docs/cli/{index,rna}.md`, `docs/nextflow/parameters.md` and
  `docs/reference/rna-annotation.md`.
- Tests: `tests/test_build_gtf_cache.py` and the group-7 GTF contract test.

The org's `gbcmsrs/buildgtfcache` module is not a release step (the operator moves it
when switching). **Recommendation:** remove as announced.

### 7. Truth arms (D9, D10, D11) and the asks they need

- **D9 GIAB (public, PHI-free; numbers can go on GitHub).**
  - Truth: HG002 v5.0q (T2T Q100-based; GRCh37 and GRCh38; repeat motif and period per
    variant), with v4.2.1 as the cross-check, the tandem-repeat benchmark v1.0.1
    (GRCh38, ≥5 bp TR indels) and stvar for 50+ bp deletions.
  - BAMs (bwa-mem, indexed, region fetch over HTTPS):
    - WGS: the DeepVariant case-study HG002 NovaSeq PCR-free 35x, GRCh38 (46 GB).
    - WES: the IDT exome at 100x (3.5 GB).
  - A few thousand sites at ±2 kb is about 0.3–0.5 GB of transfer, into the
    harness. **Ask:** OK to fetch these region slices.
- **D10 WES at paired loci.** The panel's delins and 50+ bp deletion loci in the
  TEMPO WES recaptures of the same libraries (BWA, no realignment). **Ask:** HPC
  slices (a region list and a one-line `samtools` loop, from the harness). TEMPO's
  cohort folder holds no germline calls, so WES truth comes from GIAB.
- **D11 IMPACT germline normals.** The pilot above (195 hets picked, 193 genotyped),
  kept as an arm.
- **D8 #255.** `prepare`, `check` and `attribute` read the same time records as
  `compare_panel.py` and skip runs whose latest attempt failed. Red-first test in
  `tests/test_regression_panel_tools.py`; the RUNBOOK's manual set-aside goes.

### 8. Issue structure

Tracker #196 → seven parents, each with its sub-issues (moved from #196's flat list):
- **A. Allele windows and spanning bias:** #253 (C39), C41 [new], #244, #246, #214,
  #143, #195, #178.
- **B. Evidence outside the aligned bases (RJ-13):** #254, #173, #177, #144, #212, #179.
- **C. Regression panel: truth arms and tooling:** #255, D9, D10, D11 [new].
- **D. Homopolymer twin: fix or retire `--rescue-homopolymer`:** #111, #112 (#145,
  #146), #128.
- **E. Diagnostics and reporting:** #227, #132, #180, #127.
- **F. Hygiene and infrastructure:** H4 [new], #138, #150, #216, #219.
- **G. Platforms:** #232 (#233–#237), already a parent.

**Ask:** OK to file the new parents and tickets and re-parent the existing issues.

## Tickets by parent

### A. Allele windows and spanning bias
- **C39 #253** — decision 1. Red-first battery: synthetic delins at every span and
  read start (VAF at uniform starts must be 0.5 for a het), long events, recurrent
  neighbours, spliced reads.
- **C41** — decision 2; measured on D9 and D11.
- **C42 #244** — decision 4.
- **C38 #246** — decision 5.
- **C33 #214** — closed by C39 (one anchor for both windows); its failing probe joins
  C39's battery.
- **C6 #143** — identity band for exact-length insertions. Measure on the two
  long-insertion loci and the GIAB insertions; test both policies (raw band and
  BQ-masked).
- **C24 #195** — the stale semiglobal scores in the local-alignment tail. A Rust test,
  then the fix (no reach on prepared variants).
- **R3 #178** — catalogued editing sites masked inside RNA carrier windows. Measure on
  FORTE first.

### B. Evidence outside the aligned bases (RJ-13)
- **C40 #254** — decision 3, with C15 item 2 and C18.
- **C15 #173** — items 2 (C40) and 3 (reads clipped before the flank, measured on WES).
- **C18 #177** — split reads join their molecule (primary + SA, counted once).
- **C7 #144** — clip-borne ITD carriers: the clip validated against the duplication.
- **C31 #212** — fragment end from the mate's unclipped 5' end (MC tag); decided with
  C40, since clips there are now read.
- **O5 #179** — per-variant MAPQ, clip and outside-depth rates of ALT-like vs REF reads;
  a diagnostic only.

### C. Regression panel
- **D8 #255**, **D9**, **D10**, **D11** — decision 7. The panel's next gate adds the
  arms; the GIAB results are the shareable numbers.

### D. Homopolymer twin
- **C3 #111**, with **#112** (#145, #146) and **M1 #128**: fix or retire
  `--rescue-homopolymer` (opt-in since O4). Measure how often the twin wins on the
  panel's twin-candidate stratum before deciding; retiring it closes #145, #146 and
  #128's decomposition half.

### E. Diagnostics and reporting
- **O10 #227** end-of-run QC summary; **O3 #132**; **O6 #180** (decide: orientation
  counts as a diagnostic only); **I5 #127** (only if a pipeline needs conversion
  without counting).

### F. Hygiene and infrastructure
- **H4** (decision 6), **D3 #138**, **P1 #150** (gated on the binning-invariance
  tests), **C34 #216**, **I7 #219**.

### G. Platforms
- **D7 #232**: D7a Apptainer in CI, D7b RHEL 8 wheel, D7c Python 3.10–3.14 wheels,
  D7d arm64 image, D7e macOS wheels.

## Validation standard (every count-affecting change, before merge)

Carried from 6.6.0, with this cycle's additions.

| Axis | What it exposes |
|:--|:--|
| Synthetic contracts: red-first, fuzz, uniform-start VAF | rule logic, symmetry, quality encodings |
| Panel DNA, ABRA2-realigned (RC set, the gate panel) | the production default |
| DNA without realignment: WES at the paired loci (D10) | clips, split reads |
| RNA (FORTE) | splicing, BAQ, aligner clips |
| GIAB HG002 WGS/WES (D9) | germline hets at 0.5; shareable numbers |
| IMPACT germline normals (D11) | hets at 0.5 on the production assay |

**New this cycle (C39, C40): an unbiasedness check.** Rows are counted by the rule
read from the left flank and from the right flank. A rule is accepted when the two
readings' VAFs agree, and when germline hets read 0.5 across span and tract length.

## Suggested order

1. Decisions 1–8 (operator).
2. Measurement infrastructure: D8 #255; H4; D9 and D10 slices (with the asks); D11.
3. C39 (spec RJ-23 → red-first → build → review → acceptance on every arm).
4. C40 with C18 and C15 (spec RJ-24, same path), then C31.
5. C41's findings into the docs (and the diagnostic, if decided); C42's
   Phase-3 measurement and routing.
6. The carried items by parent: D (homopolymer twin), A's small counting items (C6,
   C24, R3, C38's re-measure), E, F, G (D7 RHEL first).
7. Release: the panel gate (6.6.0's method, from 6.6.0, with the new arms).
