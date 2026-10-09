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

**Status (2026-10-08).** Draft for operator review; nothing is implemented.
- **Decided (operator, 2026-10-08):** 2 (C41: keep allele-level REF and VAF, add docs),
  4 (close #244 and replace it with the Phase-3 sub-issue), 5 (re-measure #246 and close
  it unless seen) and 6 (remove `build-gtf-cache`). Go-aheads: the GIAB region fetch,
  the WES slices on HPC, filing the issue structure, and this PR.
- **Decided (operator, 2026-10-08, after the worked examples and the GIAB truth):**
  - 1 (C39): one-sided, equal-length junction windows for every complex event (POS
    flank in DNA; in RNA the flank not cut by an exon edge).
  - 3 (C40): one window at the POS breakpoint for both alleles at 50+ bp deletions,
    whatever the representation, k = 4 clean bases each side past the microhomology
    (10 in low complexity), REF read at the POS breakpoint.
  - The shared sub-decision: a change recurring next to the event on ALT reads only is
    a larger allele (withheld; `OBSERVED_ALLELE` names it), for delins and pure
    deletions alike.
- **Issues (2026-10-08):** parents A–F are #258–#263 under #196; new tickets #264–#269;
  #244 closed as measured.

Priority: **H** high, **M** medium, **L** low. **[counts]** can change counts.
**[decide]** needs an operator decision first. **[new]** has no issue yet.

## Summary

| ID | Ticket | Pri | Flags | Issue | Parent (proposed) |
|:--|:--|:-:|:--|:--|:--|
| C39 | Long delins: read every complex event by one-sided, equal-length junction windows (no padding loss) | H | [counts] [decided] | #253 | A #258 |
| C41 | REF and VAF of indels in repeats against germline truth (allele-level vs by-base, spanning bias) | H | [decided] | #264 | A #258 |
| C40 | Long deletions: one window at the POS breakpoint for both alleles, clipped bases read (RJ-13) | H | [counts] [decided] | #254 | B #259 |
| C42 | Per-locus gap parameters: closed as measured; replaced by Phase-3 routing | L | [decided] | #244 → #265 | A #258 |
| C38 | Sibling rows: SNV sibling double-counts a compensating-mismatch carrier | L | [counts] [decided: re-measure] | #246 | A #258 |
| C33 | Exact-carrier REF and ALT windows read from different anchors | L | [counts] | #214 | A #258 (folds into C39) |
| C6 | Error-tolerant exact-length insertion matching (identity band) | L | [counts] | #143 | A #258 |
| C24 | Local-alignment fallback tail reads stale scores | L | [counts] | #195 | A #258 |
| R3 | RNA: catalogued editing positions inside carrier windows | L | [counts] | #178 | A #258 |
| C15 | C12 follow-ups: clipped carriers of pure deletions, anchors in the clip | M | [counts] | #173 | B #259 |
| C18 | Split-read evidence for long events (supplementary alignments) | M | [counts] | #177 | B #259 |
| C7 | Clip-borne ITD carriers | L | [counts] | #144 | B #259 |
| C31 | Fragment end from the mate's unclipped 5' end (MC tag) | L | [counts] | #212 | B #259 |
| O5 | Mapping-bias diagnostic (ALT reads mapped or clipped worse than REF) | M | | #179 | B #259 |
| D8 | `attribute.py prepare/check/attribute` skip runs whose latest attempt failed | M | | #255 | C #260 |
| D9 | Panel arm: GIAB HG002 truth (public, PHI-free) | H | first run 2026-10-08 | #266 | C #260 |
| D10 | Panel arm: WES (BWA only) at the panel's paired long-event loci | H | slices on HPC | #267 | C #260 |
| D11 | Panel arm: IMPACT germline hets in matched normals | M | pilot run | #268 | C #260 |
| C3 | Homopolymer decomposition: fix or retire `--rescue-homopolymer` | M | [counts] [decide] | #111 | D #261 |
| — | Decomposed-allele wins (umbrella; #145, #146) | M | | #112 | D #261 |
| M1 | Merge rows whose flavors report different alleles | M | | #128 | D #261 |
| O10 | End-of-run QC summary | M | | #227 | E #262 |
| O3 | MNP rescue in fillouts: flag a germline component | L | | #132 | E #262 |
| O6 | Read-orientation evidence for oxoG/FFPE (decide) | L | [decide] | #180 | E #262 |
| I5 | Nextflow `GBCMS_CONVERT` module | L | | #127 | E #262 |
| H4 | Remove `build-gtf-cache` and `--gtf-cache-dir` (deprecated in 6.6.0) | M | [decided] | #269 | F #263 |
| D3 | mkdocs-material 2.0 | L | | #138 | F #263 |
| P1 | Deep-bin fetch reduction | L | | #150 | F #263 |
| C34 | GTF transcripts keyed by ID alone | L | | #216 | F #263 |
| I7 | MAF dash insertion at `Start_Position` 0 | L | | #219 | F #263 |
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
| | one-sided junction windows, read from the left flank | 1,761 | 5,790 | **0.99** | 9/14 |
| | one-sided junction windows, read from the right flank | 1,820 | 5,763 | **1.00** | 9/14 |
| | junction windows from either flank (today's rule past 50 bases) | 1,869 | 6,984 | **0.87** | 9/14 |
| 4 clean, window over 50 (junctions today) | 6.6.0 | 458 | 1,356 | 1.00 | 2/4 |
| | one-sided, left / right | 466 / 465 | 1,351 / 1,350 | 1.03 / 1.12 | 2/4 |
| 70 clean, span under 20 | 6.6.0 | 9,557 | 24,053 | 1.00 | 52/70 |
| | one-sided, left / right | 10,342 / 10,304 | 26,086 / 26,017 | 1.00 | 66/70 |
| 49 with a recurrent neighbour | 6.6.0 | 2,401 | 14,964 | 1.00 | 19/49 |
| | one-sided, no guard | 6,909–6,944 | 16,188–16,370 | 1.02–1.04 | 23/49 |
| | one-sided + observed-flank guard | 1,154–1,199 | 13,485–13,633 | 0.99–1.01 | 12/49 |

**Against truth (D9, HG002 germline hets, VAF 0.5; WGS 35x, bwa-mem, no realignment).**
The delins are rebuilt from the v5.0q records dipcall splits (same-phase records within
10 bp), and isolated from other variants.

| HG002 het delins | Rule | Rows | VAF median (IQR) | Rows within 0.35–0.65 | AD | RD |
|:--|:--|--:|:--|--:|--:|--:|
| span 20+ (250) | 6.6.0 | 237 | 0.464 (0.357–0.544) | 166 | 2,731 | 3,258 |
| | no padding | 245 | 0.469 (0.363–0.551) | 176 | 3,147 | 3,663 |
| | one-sided, left / right | 246 | **0.470 / 0.481** (0.38–0.55) | **186 / 188** | 3,425 / 3,473 | 3,885 / 3,865 |
| | either flank | 246 | 0.471 (0.381–0.561) | 183 | 3,751 | 4,087 |
| span 5–19 (150) | every rule | 149 | 0.500 | 124–130 | 2,087–2,306 | 2,143–2,446 |

- **Padding costs reads, not bias.** On the panel's long delins the padded windows
  drop a third of the reads that one-sided windows judge, ALT and REF alike, with the
  same VAF (0.99–1.00). The carriers lost are the ones the gate saw: sign-out agreement
  0/14 → 9/14. On HG002, one-sided windows judge 25% more ALT reads, sit closest to 0.5,
  and put the most rows in the 0.35–0.65 band.
- **The two biased alternatives behave as predicted** on the panel, where the length
  change is large and reads are 101 bp. Whole alleles without padding favour the
  shorter allele (VAF ×1.22). Reading junctions from either flank favours the longer
  one (×0.87), because a REF read can prove REF at either end while the shorter ALT
  has one window. On HG002 the biases are small: its delins mostly change length by a
  few bases, and its reads are 151 bp.
- **Read-start arithmetic predicts the row-level effects.** A worked example is a
  35 bp → 1 bp delins with 101-bp reads. The padded and one-sided rules each give REF
  and ALT the same number of read starts that can judge them (63 and 63; 97 and 97).
  No padding gives 63 and 97, predicting VAF 0.53; 0.54 was measured. Either flank
  gives 132 and 97, predicting 0.35; 0.36 was measured.
- **Past 50 bases, 6.6.0 shows no such bias in its counts.** Its rule reads junctions
  from either flank, but only reads over the variant's first base (or admitted by a
  clip) are counted, so a REF read holding only the far junction is classified and
  never counted. On the 4 clean long delins its written REF (1,356) matches the
  one-sided reading (1,351). (Corrected 2026-10-09: an earlier draft read the engine's
  per-read trace, which also logs those uncounted reads, and reported VAF ×0.64–0.71.)
- **Recurrent neighbours.** In 49 rows the well-anchored ALT carriers share a confident
  change at a fixed distance (up to 25 bp) from the event. In 47 it is on ALT reads
  only: the reads carry a larger haplotype than the given allele. In 2 it is also on
  REF reads (likely a germline SNP in cis). 6.6.0's padding withholds part of these.
  Pure deletions meet the same shape (decision 3), and 6.6.0 credits them there under
  RJ-8, so the two paths disagree today; see sub-decision (a)/(b) below.

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
| A | Keep 6.6.0 | baseline | unbiased | loses ~30% of carriers on long delins |
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
- **The anchor flank.** One flank per variant, the same for both alleles. In DNA it's the
  POS (left) flank, the breakpoint pure deletions count REF at today. A read whose own
  clip holds that flank is read there (RJ-13). Reads holding only the other junction
  become depth for both alleles alike. Any allele-independent choice keeps the windows
  unbiased; the panel's left and right readings agree (0.99 and 1.00).
- **RNA.** The rule is the same for RNA's complex events: its windows are built over the
  reference spliced at the read's junctions (RJ-19), and exon-edge clips are not read
  (RJ-18). One addition: where the POS flank is cut by an exon edge, reads come from the
  exon side, so the variant reads from the other flank. The choice is structural (from
  the exon edges), so it still treats both alleles alike.
  - Reach today: the FORTE truth arm has 94 rows, one of them a delins, so its counts
    barely move.
  - Community: RNA callers add nothing specific for complex events. GATK's RNA best
    practice is SplitNCigarReads, then the DNA likelihood rule; the ASE counters
    (ASEReadCounter, phASER) count SNVs only.
  - Validated with the RJ-19 spliced contracts, plus synthetic delins at exon edges.

**Recommendation: C.** It departs from practice (no tool reads equal one-sided
windows), but it measures better than each alternative. It keeps the unbiased VAF the
padding was there for, reaches the carriers likelihood callers and the sign-out count,
and replaces two rules (padded windows, either-flank junctions) with one. D needs a column and lets counts and
VAF disagree. One sub-decision remains: a change that recurs on every ALT read next to
the event, on ALT reads only. It applies to delins and pure deletions alike (3 of the
panel's 31 long-deletion rows have one, 3–7 bases from the junction):
- **(a) a larger allele (invariant 7):** grow both windows past it, so a read must show
  it is the given allele. AD on those rows falls (delins 2,401 → about 1,180). The
  diagnostic names what the reads carry (`OBSERVED_ALLELE`).
- **(b) a separate event (RJ-8):** credit the given allele, as 6.6.0 does for pure
  deletions today. Delins AD on those rows rises to about 6,900, pure deletions keep
  theirs.

Either way the two paths should agree. I lean (a), because these neighbours are on ALT
reads only (somatic haplotypes larger than annotated). But it's a policy call: (b)
matches what the sign-out credits.
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
- **What the AD drop is** (6.6.0 ALT reads not kept at k = 4):

  | | 50–99 bp | 100+ bp |
  |:--|--:|--:|
  | Kept | 92% | 82% |
  | Under 4 bases past the junction on one side, past the microhomology | 5.8% | 4.5% |
  | A confident change next to the junction | 2.3% | 13.8% |

  - The short-junction reads (the ~5%) show too little to tell the alleles apart. Equal
    windows drop the matching REF reads too, so VAF holds.
  - The 13.8% sits in 3 of 31 rows (lengths 113, 279, 327 bp). In each, nearly every ALT
    read has the exact-length D at the given position plus one recurring confident
    change 3 or 7 bases from the junction. That's the recurring-neighbour question of
    decision 1, sub-decision (a)/(b). Under (b) those reads stay ALT, and the drop at
    100+ bp is about 5%, like the 50–99 bp rows.
- **Which breakpoint REF is read at matters on capture data, row by row.** On 48 panel
  deletions with REF at both breakpoints, VAF with REF at the left against the right
  has median 1.02, so no bias on average. But it spans 0.69–1.86 (10–90%), and 19 of 48
  rows differ by more than 25%: capture depth varies along the event. On WGS (HG002)
  it doesn't.
- **Against truth (D9, HG002 germline het deletions of 50+ bp, VAF 0.5; WGS 35x,
  bwa-mem, no realignment), 280 rows:**

  | Rule | Rows judged | VAF median (IQR) | Within 0.35–0.65 | AD | RD |
  |:--|--:|:--|--:|--:|--:|
  | 6.6.0 | 170 | **0.000** (0.000–0.000) | 0 | 31 | 3,409 |
  | Equal windows at POS, k = 4 | 221 | **0.500** (0.368–0.590) | 140 | 2,826 | 2,936 |
  | Equal windows at POS, k = 8 | 213 | 0.500 (0.376–0.592) | 134 | 2,608 | 2,718 |
  | Equal windows at POS, k = 10 | 213 | 0.500 (0.369–0.592) | 128 | 2,497 | 2,627 |

  - 6.6.0 credits no ALT at 287 of 300 such hets: without realignment the carriers are
    clips (1,719 of the 2,497 junction carriers at k = 10).
  - C18's WES record says the same at 68–106 bp, and BWA-MEM writes a D only when the
    far side outscores the gap (s2). This is the WES/WGS failure #254 predicted.
  - Equal windows read 0.500 at every k; k only trades reads (k = 4 judges 13% more
    than k = 10).

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
| A | Keep | baseline | 4–7% low | ≈0 ALT: HG002 WGS VAF 0.000 at 287/300 het deletions |
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
- Truth backs it: on HG002 (BWA only) equal windows read 0.500 where 6.6.0 reads 0.000.
- Scope: deletions of 50+ bp (RJ-6's regime). The 20–49 bp range is measured on WES
  first.
- Sub-option for events longer than a read, on capture data: REF read at the POS
  breakpoint (as today; per-row VAF can swing with depth along the event), or the mean of
  both breakpoints' REF rounded to a read (steadier; a derived count). I lean POS (a read
  count stays a read count) and putting the swing in the docs.

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

Filed 2026-10-08: A #258, B #259, C #260, D #261, E #262, F #263; C41 #264, Phase-3
routing #265, D9 #266, D10 #267, D11 #268, H4 #269.

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

1. Decisions 1–8 (operator): all decided 2026-10-08.
2. Measurement infrastructure: D8 #255; H4; D9 and D10 slices (with the asks); D11.
3. C39 (spec RJ-23 → red-first → build → review → acceptance on every arm).
4. C40 with C18 and C15 (spec RJ-24, same path), then C31.
5. C41's findings into the docs (and the diagnostic, if decided); C42's
   Phase-3 measurement and routing.
6. The carried items by parent: D (homopolymer twin), A's small counting items (C6,
   C24, R3, C38's re-measure), E, F, G (D7 RHEL first).
7. Release: the panel gate (6.6.0's method, from 6.6.0, with the new arms).
