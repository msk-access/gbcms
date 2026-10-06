# 6.6.0 cycle — plan

> Every open finding after the 6.5.0 cut, whatever its priority. Each ticket has
> a GitHub issue in the **6.6.0 milestone**, under the tracking issue **#140**
> (grouped work as sub-issues). Sources: the
> open issues (#92, #106, #111, #112, #114); the leftovers in `CONTINUITY.md` and
> `CYCLE_6.5.0_PLAN.md`; the #116 adversarial reviews and the plan-vs-code
> verification of 2026-09-25; and the deferred items of
> `CODE_REVIEW_IMPLEMENTATION_PLAN.md`. Branch-from: `develop` after the 6.5.0
> cut.
>
> **Discipline (unchanged from 6.5.0).**
> - Every count-affecting ticket is measured first and red-first
>   (xfail-strict battery committed before the fix).
> - Counts must not depend on bin geometry (binning invariance), and
>   classification is checked against the read census (T1/T2, 2026-10-01; until
>   then, changes were mirrored in the legacy parity oracle).
> - Each ticket gets an adversarial review before merge, then real-data
>   acceptance on the local truth sets. Patient data stays local; issues and
>   PRs stay PHI-free.
> - No new output columns: new signal goes to `gbcms_diagnostic`, logs, or
>   validation tooling.
> - Check `.agents/learnings/REJECTED.md` before proposing an alternative.

Priority: **H** high, **M** medium, **L** low. **[counts]** marks a ticket that
can change counts. **[decide]** marks a ticket that needs an operator decision (**[decided]**: the decision is recorded under the ticket)
before implementation. **[6.7.0]** marks a ticket moved to the 6.7.0 milestone and
**[closed]** one closed by the triage of 2026-09-30 (see "Triage" below). **[in review]**
marks a ticket with an open PR.

## Summary

| ID | Ticket | Pri | Flags | Issue |
|:--|:--|:-:|:--|:--|
| C1 | Complex variants count exact carriers (was: partial-ALT in the SW local fallback) | L | [counts] [done] | #141 (#92) |
| C2 | REF fragments at grouped rows (main vs per-transcript) | M | [counts] [done] | #119 |
| C3 | Homopolymer decomposition arbitration redesign | M | [counts] [6.7.0] | #111, #145 (#112) |
| C4 | Reference windows near contig ends | M | [counts] [done] | #142 (#92) |
| C5 | Long insertions exceed the pangenomic matrix cap | L | [counts] [closed] | #120 |
| C6 | Error-tolerant exact-length insertion matching | L | [counts] [6.7.0] | #143 (#92) |
| C7 | Rescue for clip-borne ITD carriers | L | [counts] [6.7.0] | #144 (#92) |
| C8 | One-base-REF delins without a shared anchor | L | [counts] [done] | #121 |
| C9 | Count a MAF deletion at Start 1 | L | [counts] group 4 [done] | #122 |
| C10 | Reads ending inside an indel's repeat tract counted REF | H | [counts] [done] | #157 |
| C11 | Phase-3 context misses tandem duplications longer than the repeat finder's motifs | M | [counts] [closed] | #159 |
| C12 | Count carriers whose allele lies in soft-clipped bases (complex variants) | H | [counts] [done] | #167 |
| C13 | BAQ spares the variant's own indel evidence | H | [counts] [done] | #166 |
| C14 | Records without bases (SEQ `*`) crash the SNP path | M | [counts] [done] | #172 |
| C15 | C12 follow-ups: RNA, clipped pure deletions, anchors in the clip | M | [counts] RNA part group 3 (RJ-18, built); the rest [6.7.0] | #173 |
| C16 | Stray ALT calls at RNA exon-edge probes | L | [counts] [decided] [done] | #174 |
| T1 | Test architecture: retire the legacy parity path for binning-invariance tests | M | [done] | #170 |
| T2 | Read census as the classification oracle in tests | M | [done] | #171 |
| H3 | Code-quality sweep of the cycle's code: duplication, unused code, silent failures, comments, logging, monitoring | M | [done] | #204 |
| C17 | A read ends at its fragment end: adapter-like read-through bases are clipped in every read | M | [counts] [decided] [done] | #176 |
| C18 | Split-read evidence for long events (supplementary alignments) | M | [counts] [6.7.0] | #177 |
| C19 | Absent base qualities (QUAL `*`, read as 0xFF) overflow fragment consensus | M | [counts] [decided] [done] | #182 |
| C20 | ALT carriers ending inside an indel's repeat tract credited from the CIGAR gap (C10's ALT side) | M | [counts] [done] | #188 |
| C21 | Homopolymer insertion carriers placed elsewhere in the run counted REF (S3 anchor-base test) | H | [counts] [done] | #189 |
| C22 | Same-length non-equivalent deletions ≥5bp near a deletion row reach Phase 3, which calls them ALT | M | [counts] [done] | #191 |
| C23 | Distinct alleles in long-period repeats (motif > 6bp) keep REF: "in a repeat" is decided by `repeat_span` | L | [counts] [done] | #192 |
| C24 | Local-alignment fallback tail reads stale semiglobal scores (rare no-reference path; partial_alt only) | L | [counts] [6.7.0] | #195 |
| C25 | Exact-carrier long-event junction windows: a read holding one junction decides the call (REF for anchor-keeping reads; ALT for reads carrying only the substitution) | M | [counts] [decided] [done] | #199 |
| C26 | Which of a read's other indels decide its REF call: a short one inside the window counts REF, a ≥5bp one outside it withdraws REF | M | [counts] [decided] [done] | #200 |
| C27 | A read spelling the ALT across several indel ops is judged by its ops, not its bases | L | [counts] [decided] [done] | #201 |
| C28 | A read deleting a pure deletion's anchor falls back to Phase 3, which credits the closer haplotype | M | [counts] [decided] [done] | #202 |
| C29 | `check_complex`'s inline query walk counts hard clips (a hard-clipped read's anchor quality is read from the wrong base) | L | [counts] [done] | #207 |
| C30 | Two pure-indel tests disagree with the allele-kind classification on lowercase or unprepared alleles | L | [counts] group 4 [done] | #208 |
| C31 | Measure the fragment end from the mate's unclipped 5' end (MC tag), so a mate's clipped 5' end does not clip molecule bases as adapter | L | [counts] [6.7.0] | #212 |
| C32 | Junction-placement guard for spliced reads at exon-edge events (needs reference access in the engine) | M | [counts] group 3 (RJ-19) [done] | #213 |
| C34 | GTF transcripts keyed by ID alone: an ID reused on another chromosome merges (introns, per-transcript counts, ASJD); 0 such IDs in Ensembl 111 | S | [6.7.0] | #216 |
| C33 | Exact-carrier REF and ALT windows read from different anchors can miss a contradicting base | L | [counts] [6.7.0] | #214 |
| R3 | RNA: catalogued editing positions inside carrier windows | L | [counts] [6.7.0] | #178 |
| R4 | Gene strand unresolved at intronic loci (splice sites) and opposite-strand overlaps | M | [counts] group 3 (decided) [done] | #185 |
| R5 | C10's informative rule counts a read's splices as reference coverage (RNA reads spliced inside a repeat tract) | L | [counts] group 3 (RJ-17) [done] | #198 |
| O5 | Mapping-bias diagnostic (ALT reads mapped or clipped worse than REF) | M | [6.7.0] | #179 |
| O6 | Read-orientation evidence for oxoG/FFPE artifacts | L | [decide] [6.7.0] | #180 |
| O7 | Unmapped mates (flag 0x4) placed at a variant count in `mq0_count` | L | [decided] [done] | #183 |
| O8 | `OBSERVED_ALLELE`/`COEXISTING_ALLELE` read antisense reads under enforcement (no NH rescue) | L | group 3 [done] | #186 |
| O9 | VCF output of a whole-contig MAF deletion at Start 1 writes REF past the contig end | L | [decided] group 5 (built; from 6.7.0) | #220 |
| R1 | Span-aware exon-edge BAQ rule | L | [counts] [decided] [done] | #106 |
| R2 | RNA strandedness gating observability | M | [decided] [done] | #114 |
| I1 | MAF allele base check | M | [decided] group 4 [done] | #123 |
| I2 | `End_Position` optional | L | group 4 [done] | #124 |
| I3 | VCF→MAF `Tumor_Seq_Allele1` | L | [decided] group 4 [done] | #125 |
| I4 | maf2vcf's second ALT from `Tumor_Seq_Allele1` | L | [decided] group 4 [done] | #126 |
| I5 | Nextflow `convert` module | L | [6.7.0] | #127 |
| I6 | `REF_MISMATCH` rows say where the given REF sits | L | [decided] group 4 [done] | #218 |
| I7 | A MAF dash insertion at Start 0 is skipped, not a FAIL row | L | [6.7.0] | #219 |
| M1 | Merge rows whose flavors report different alleles | M | [6.7.0] | #128 |
| M2 | Merge inputs from different gbcms versions | M | [decided] group 5 (built) | #129 |
| M3 | Decomposed-allele hardening (observations, list length) | M | | #147; #146 [6.7.0] (#112) |
| M4 | Merge sums NA/nan count cells as 0 silently (the documented warning was never implemented) | M | [decided] group 5 (built) | #194 |
| M5 | Merge leaves the annotation columns of a row only a later input has empty | L | [decided] group 5 (built) | #221 |
| M6 | Merge takes mFSD and RNA columns from the first input only | M | [decided] group 5 (built) | #223 |
| M1a | Merge: combined columns NA where the flavors' MNP rescue outcomes differ (rescue half of M1, from 6.7.0) | L | [decided] group 5 (built) | #224 |
| O1 | UMI and no-bases warnings repeated by the rescue recount | L | [decided] group 5 (built) | #130 |
| O2 | Run-start summary of enabled options | L | [decided] group 5 (built) | #131 |
| O3 | Rescue in fillouts without the MNP | L | [6.7.0] | #132 |
| H1 | Writers closed when a write fails | L | [decided] group 5 (built) | #148 |
| H2 | `is_indel` in preparation | L | group 4 [done] | #149 |
| P1 | Deep-bin fetch reduction (M5b) | L | [6.7.0] | #150 |
| P2 | Bin cost-sort (PF-2) | L | [closed] | #151 |
| P3 | Document the bin-span soft floor (LO-3) | L | | #152 |
| S1 | Mean LLR per fragment (CR-5) | L | [decided] | #153 |
| S2 | `MIN_FOR_KS` floor (ME-9) | L | [decided] | #154 |
| D1 | CI version-consistency check | M | | #136 |
| D2 | Release workflow creates the GitHub Release | M | | #137 |
| D3 | mkdocs-material 2.0 | L | [6.7.0] (pin: D3a #225, group 5) | #138 |
| D3a | Pin the docs toolchain below MkDocs 2.0 | L | [decided] group 5 (built) | #225 |
| D4 | Dependency upgrade audit: does anything break on current releases? | M | | #139 |
| D5 | Coverage-driven regression panel (replaces the 56-sample matrix) | M | | #155 |
| D6 | One QC-flags reference page | M | [decided] group 5 (built) | #156 |

## Counting correctness

### C1 — Partial-ALT evidence in the Smith-Waterman local fallback (#141, under #92) · L [counts] [done]
**What the code does.** Under the Smith-Waterman backend (`--alignment-backend
sw`, or its fallback when the default method fails), a read at a *complex*
variant (both alleles longer than one base, unequal lengths) is classified in
two tries:
1. semiglobal alignment of the whole read to the REF and ALT haplotypes;
2. only when the first try is weak or too close to call, local alignment
   (which ignores badly matching read ends). This exists for sign-out alleles
   that are slightly wrong or incomplete.

The second try's scores decide REF / ALT. But the "this read shows some ALT
evidence" flag, which feeds `partial_alt`, is still computed from the first
try's discarded scores. That inconsistency is what #92 reported
(`rust/src/counting/alignment.rs`).

**Measured (2026-09-25; local data, aggregates).**
- 129 real reads took the second try: 40 signed-out complex variants in 40
  IMPACT samples, SW backend. Each read was compared with what it carries
  (ALT- or REF-specific 8-mers).
- **The originally planned fix (flag from the second try's scores) does not
  help.** The flag is right for 36 of 67 REF/tie reads today, and for 32 with
  that fix. Both fail for the same reason: the flag fires when the two scores
  are *close*, not when the read *contains ALT sequence*. For example, 28 tie
  reads showing neither allele get flagged.
- **A second finding.** 22 of the 62 reads the second try calls ALT (35%)
  carry no ALT-specific sequence, and they count in `alt_count`. Some may carry
  a real event whose sign-out allele is wrong; others may be noise.
- **Scope.** The default `pairhmm` backend never reaches this code, and the SW
  fallback under it fired 0 times on the traced real runs. The affected users
  are those who choose `--alignment-backend sw`: about 3 reads per complex
  variant.

**Options.**

| | Change | Measured | Risk |
|---|---|---|---|
| A | Flag from the second try's scores (the original plan) | Slightly worse (32 vs 36 correct) | Adds error; ruled out |
| B | Flag by read content: only when the read carries ALT-specific sequence | Removes the ~31 wrong flags | `partial_alt` drops for complex variants under SW |
| C | B, plus: second-try ALT calls without ALT sequence become partial | Also corrects up to 35% of those ALT calls | `alt_count` changes |
| D | Leave it; document that `partial_alt` is score-based under SW | None | Known inaccuracy stays; ruled out |

**Principle (operator, 2026-09-25): count the given allele.** A read counts ALT only
if it carries the given ALT. A recurring unannotated haplotype is a diagnostic
finding (`PARTIAL_DOMINANT` already surfaces heavy partial evidence), never ALT
for the row. This removes the exception C originally had ("unless reads share a
recurring unannotated haplotype"), and makes C the principled option: the local
fallback exists to credit reads when the sign-out allele is slightly wrong,
which is the deconvolution the principle rules out.

**Decision pending: B vs C.** It rests on a verification round with two
questions:
1. **Real carrier or noise?** Are the ALT calls without ALT sequence real
   carriers (they share a recurring unannotated haplotype) or scattered noise?
2. **Independent check of B.** B was scored with the same 8-mer rule it would
   use. So each fallback read is linked back to its BAM read and re-judged
   C2-style: rebuilt across the variant window and matched to REF / ALT /
   OTHER by edit distance.

**Round 1 (2026-09-25): inconclusive.**
- Only 22 of the 129 fallback reads could be linked to a single BAM read. The
  traced sequence also matches the overlapping mate, and the link demanded
  exactly one hit.
- 21 of those 22 did not cover the judging window (the MAF span ±10), so the
  independent judge could not rule on them.
- B's apparent 15/15 on REF/tie reads is therefore vacuous: B raised no flag,
  and the judge had no verdict to compare with.
- It does suggest something to check: the fallback may fire mostly on reads
  that do not span the event.

**Round 2 design.**
- Link by read name, so the two mates of a fragment count as one.
- Judge over the discrimination window built for C10 (the event ±1 base),
  not ±10.
- Report how many fallback reads are uninformative under C10's rule. If most
  are, C10 already settles them (depth only), and the B-vs-C question shrinks
  to the reads that span the event.

**Round 2 (2026-09-25): points to a different fix.**
- **Linking still fails for most reads:** 24 of 129 linked. The trace prints
  only the read's sequence over the variant window, which many reads share, so
  exact linking needs the read name in the trace line.
- **Of the 24 linked, 16 do not span the event ±1.** The local fallback fires
  mostly on reads that end in or at the event.
- **A pure k-mer rule would demote real carriers:** 1 ALT call without an
  ALT-specific 8-mer was judged ALT by the census.

**Recommended direction (replaces B vs C).** Apply C10's rule to complex
variants: a read is REF or ALT only if it spans the whole event (span ±1).
Under "count the given allele", only such a read can carry the given ALT, and
only such a read can rule it out. Reads that end inside the event count
toward depth only.
- This settles most fallback reads, whose local-alignment rescue is then moot.
- The `partial_alt` flag for spanning reads follows what the read contains.
- To confirm: add the read name to the fallback trace, then rerun the
  verification with exact linking on the 40 variants.

**Round 3 (2026-09-25; read names in the trace, #162): 110 of 129 fallback reads linked exactly.**

| Fallback reads | Count | Calls | Against what the read shows |
|---|---|---|---|
| Don't span the event (±1) | 77 (70%) | 30 ALT, 10 REF, 37 tie | Unjudgeable: the read doesn't show the whole allele |
| Span the event | 33 | 23 ALT, 4 REF, 6 tie | ALT 23/23 right, REF 3/4; ties lost 2 ALT and 2 REF reads that had a clear answer |

- 30 of the fallback's 53 ALT calls come from reads that don't span the event.
- The 8-mer content rule (B) would demote 6 of the 23 true ALT carriers, so it
  is ruled out as a detector. Round 1's "35% without ALT sequence" came largely
  from non-spanning reads and detector misses.

**How the standards do it** (VarDictJava and GATK Mutect2 source; see #141):
- **Neither counts exactly the given allele.**
  - VarDict matches its own CIGAR-derived allele text, then rescue moves
    near-matching reads in (soft-clip consensus with up to 3 mismatches).
    Its REF is a one-base count, and it has no given-allele mode.
  - Mutect2 credits each read's best haplotype by likelihood, even when the read
    carries extra events. `--alleles` injects the given allele, so reads carrying
    an unassembled, slightly different allele can be credited to it.
- **Mutect2 counts a read if its bases tell the alleles apart**, even when it
  ends partway through the event. That is option (b) below. It drops
  uninformative reads from both AD and DP.

**Recommended: (a) REF and ALT only from reads that span the event.**
- **Option (a)** is symmetric and follows "count the given allele": only a
  spanning read shows a delins's whole inserted sequence. On these reads it
  removes 30 ALT and 10 REF calls from non-spanning reads. 3 of those ALT
  reads keep their fragment's ALT through a spanning mate.
- **Option (b)**, Mutect2-style, would keep reads that show part of the ALT,
  which the principle rules out. Counting REF from partial reads while
  requiring full ALT would bias VAF down.
- **Spanning ties:** decide them by a full haplotype comparison of the read
  over the event (edit distance, as the census does), not by 8-mers.
- **Expected differences:** gbcms REF will be below VarDict's one-base REF at
  long events and repeats; document this for D5 comparisons against
  VarDict-called sign-out.

**Rounds 4 and 5 (2026-09-25): soft clips and exact carriers.**
- **Soft clips count.** 11 of the 12 ALT-called fallback reads clipped at the
  event carry the whole given ALT in their clipped bases. So "spans" is judged
  on the read's own sequence, clipped bases included, not on aligned extent.
- **Given alleles are almost always right, but the fallback picks the
  imperfect reads.**
  - Across every read at the 40 variants, 96% of ALT-like reads carry the given
    ALT exactly.
  - Only 2 of the fallback's 23 aligned-spanning ALT calls do. It fires where
    end-to-end alignment was poor, then credits near-matches.

**Default backend (pairhmm), read by read** (a named trace per classified read,
#162; 38 variants; anchor-overlapping reads judged over the event ±2):

| Default call | Exact given allele | Partial (ends within 2 of the event) | Carries neither exactly |
|---|---|---|---|
| ALT (4,185) | 3,767 (90%) | 384 (9%) | 27 |
| REF (13,390) | 11,017 (82%) | 2,278 (17%) | 69 |

- Near-match crediting is rare on the default backend.
- Its issue is C10's, on complex variants: it counts reads that end at or
  inside the event, and credits them to REF almost twice as often as to ALT.

**Decision (2026-09-25): an exact-carrier rule for complex variants, on both
backends.**
- ALT: the read's own bases (aligned or soft-clipped; bases below min BQ
  masked) contain the given ALT with 2 reference bases of flank on each side.
- REF: the same test against REF. Anything else is neither.
  - A read closer to ALT counts as `partial_alt`.
  - A recurring off allele is named in `gbcms_diagnostic`.
- Tolerance comes from base quality only: one quality rule across backends.
- On the 38 variants it keeps 90% of ALT calls and 82% of REF calls, and VAF
  rises by a median of 3.8%.
  - 12 variants move by more than 10% (0.84× to 2.84×).
  - Acceptance adjudicates each mover read by read before release, as C10's
    outliers were.
- **Design notes:**
  - Extend the flank past any repeat the event's ends can slide into, reusing
    C10's shift-region machinery.
  - The local-alignment fallback no longer decides complex calls, so #92's
    stale flag goes with it.
  - It builds on #160 (`counting/window.rs`).
  - **Compare placement-independently, in canonical form.** Widen each read
    over any indel whose shift region touches the event, and compare left-aligned
    minimal alleles, as O4's scan does (#164). O4 showed that a fixed-window
    comparison misreads carriers whose indel the aligner placed elsewhere in a
    repeat: 4 of 5 "mis-described" complex variants were such artifacts. So the
    90%-of-ALT-kept figure above is probably conservative.

**Implemented (2026-09-26, `feature/c1-exact-carrier`).** The exact-carrier rule
(`rust/src/counting/carrier.rs`) for delins, Del+SNV and MNP reads with an indel
in or right beside the block, on both backends:
- **Windows.** The event (union of the left- and right-first trims, grown through
  tandem repeats touching it on either allele, unit 1–6 or the length change)
  plus 2 flank bases, read at the event's own position from both anchors. The
  shorter window is padded to equal length. Over 50 bases, equal junction
  windows at each end read the shorter allele where it fits; a mismatch at any
  junction rules the read out.
- **Outcomes.** A read decides only if it holds both windows; one that cannot is
  depth only (no `partial_alt`, no mFSD class). The previous classifier remains
  only for unprepared variants and events at a contig end; prep widens
  `event_ref` as needed (≤16,384 per side, never past a contig end).
- **Multi-allelic guard.** A read matching the row's ALT and a sibling's exactly
  (they differ only at bases the read has masked) is neither row's AD.
- **Review.** Two independent adversarial rounds on synthetic data found 11
  defects: position-free matching, run-edge placement, long-event gaps and bias,
  partial from truncated reads, the MNP gate's reach, prep at contig ends, and
  records without bases. Each was fixed red-first (`test_complex_exact_contract.py`,
  35 tests). Fuzz: 700 random delins under left and right gap placement and
  ~400 MNPs, all clean. VAF 0.495–0.505 on every shape with uniform read starts.

**Effects map.**
- **Changes:** `ref_count`, `alt_count`, `partial_alt`, `any_alt`, the fragment
  counts, mFSD classes, the observations export, ASJD for RNA complex rows,
  `PARTIAL_DOMINANT`, VCF `PAD`/`AAD`, `mnp_confirmed_alt` (MNP rescue keeps a
  haplotype its carriers show whole), and `event_ref` (the observed-allele scan
  can see a wider reference). The guard's exact-tie rule covers all grouped rows.
- **Unchanged:** SNVs, pure indels outside such ties, and MNP reads without an
  indel at the block.

**Acceptance (local data; aggregates only).**
- **Backends:** `sw` equals `pairhmm` on 25 of 25 traced complex variants.
- **Read by read** (25 signed-out IMPACT complex variants; truth = equal-length
  padded windows, reads holding both): exact REF reads are called REF 7,457 of
  7,468 times (11 neither); exact ALT reads ALT 2,018 of 2,020 (2 partial). No
  REF↔ALT flip. ALT calls are 98.8% exact-ALT reads (develop 89.8%); REF calls
  98.5% exact-REF (develop 85.8%). Reads not holding the windows: REF 2.9%, ALT
  2.8%.
- **RC set:** 8 of 1,154 rows change, all complex or in a complex row's group.
  Movers against a BAM census: 30>14, 24>3, 2>17 and 4>3 match within 0.007 VAF;
  3>5 in a repeat matches an edit-distance census (0.546 vs 0.523); 2>3 beside a
  sibling insertion found the multi-allelic tie (fixed: 1 ALT, census 0–3).
  Its two grouped deletion rows lose 1 and 3 `partial_alt` reads to the sibling.
- **RNA (FORTE; ~7,150 synthetic delins probes in expressed genes, 2 samples;
  REF share of depth, develop → C1):** genomic windows reached into the intron,
  so spliced reads were neither: 0.987 → 0.017 at 0–1 bp from an exon edge,
  0.982 → 0.484 at 2–4 bp. Windows now end at a read's own junction: 0.923 and
  0.906, the same as mid-exon (0.909). The ~8% below develop everywhere is
  reads ending inside the window (depth only by design). Pure indels (C10) show
  no exon-edge effect (1–3%, flat). Default RNA BAQ removes indel and complex
  ALT reads when a base's quality minus 20 falls below min BQ (#166); FORTE's
  Q40 bins escape it.
- **Open (operator):** read inclusion, not the rule. Aligners clip ALT reads
  that end or start within ~10 bases of an event, and clipped reads are outside
  depth, so a 50% sample reads 0.45–0.46 VAF at small events (0.37–0.44 at 20
  bases). Counting soft-clip carriers in depth would remove that bias.

**Priority: done in 6.6.0** (was L; the decision widened C1 to both backends).

### C2 — REF fragments at grouped rows (#119) · M [counts] [decided]
**Finding.** At rows in a multi-allelic or tract-cluster group, the main counts
record a fragment as REF *before* the REF-side sibling guard runs. So a read
excluded from `ref_count` (it carries a sibling's ALT) still counts in
`ref_count_fragment`. The per-transcript counts exclude it from both. This
predates 6.5.0: the guard is the original multi-allelic one, and T1 widened
the groups.
**Decision.** Align the main counts with per-transcript (exclude the fragment
too), or document the difference. Recommended: align, since REF testimony from
a read that carries a sibling's ALT is vacuous for this row's fragments too.
**Direction.** Run the REF-side guard before fragment evidence is recorded,
in both the binned and legacy paths (`count_bam` parity).
**Measure first.** RDF deltas on the T1 acceptance clusters (ACCESS BRCA2
duplex and simplex, the complex-cluster IMPACT samples, MSI-high).
**Acceptance.** Only grouped rows' RDF moves, and downward; `dpf ≥ rdf + adf`
holds; the parity suite stays green.

**Decision (2026-09-25, operator).** Window-aware rule: a read is excluded from REF at A only when a sibling's event lies inside A's discrimination window (its repeat tract, or its span in unique sequence, ±1 base). The same rule applies to read and fragment counts. This is GATK's overlap rule with the tract standing in for representation shifts. It uses the same window as C10. Evidence on #119: the IGV-lens census shows 1,142 reads wrongly dropped from `ref_count` and 439 fragments wrongly kept in `ref_count_fragment`.

### C3 — Homopolymer decomposition arbitration redesign (#111; #112 item 1 is #145) · M [counts]
**Finding.** When a delins looks like a miscollapsed homopolymer event, two
permissive classifiers compete: the called allele and a corrected allele
(`REF[..len-1] + X`). The margins are thin, and at 4 of the 11 real twin loci
most reads carry a third allele that both claim. When the corrected allele
wins, per-transcript counts, ASJD and `NON_DISCRIMINATING_LOCUS` still
describe the original (#112 item 1).
**Decision (2026-09-25, operator): count the given allele; the twin becomes opt-in;
reads' own allele is named by a diagnostic (O4).**
- **Why it was added** (commit 94e06f70, 2.6.0): callers sometimes merge a
  1bp deletion plus an SNV in a homopolymer into one inflated delins. At SOX2,
  `CCCCCC→T` was signed out with ALT 3, while 79 reads carried
  `CCCCCC→CCCCT` (the 1bp deletion plus C→T).
- **But the code builds a different allele than intended.** The docs'
  arithmetic hid it: `C×(len−1)+T` is `CCCCCT`, a same-length SNV at the run's
  end, not `CCCCT`. So even SOX2 is won by tolerance, not by an exact match.
- **Its "self-validating" claim fails.** At 4 of the 11 real twin loci the
  twin claims 93–97% of the called allele's exact carriers.
- **What changes:**
  - By default, the dual count is off: the row counts the given allele.
  - `--rescue-homopolymer` keeps today's twin dual count, flagged
    `WARN_HOMOPOLYMER_DECOMP`, like `--rescue-mnp`. Nextflow follows the CLI
    default.
  - O4's `OBSERVED_ALLELE` names what the reads carry. At SOX2 that would
    say `CCCCT`, 79 reads exact, against 3 for the given allele.

Earlier analysis:
**Conflict with "count the given allele" (2026-09-25).** The
decomposition is on by default. Every eligible deletion is counted twice, as
given and as a corrected allele. When the corrected allele gets more ALT reads,
its counts are reported under the row's label (`WARN_HOMOPOLYMER_DECOMP`).
That is deconvolution of a possibly wrong input, on by default, unlike the
opt-in `--rescue-mnp`.
**Recommended direction** (replaces the 6.5.0 plan § T11 direction):
- Count the given allele, always.
- Measure the corrected shapes' exact haplotype support (the census method) as
  a diagnostic, and flag when the reads carry a corrected shape or neither.
- If corrected counts are still wanted, put them behind an opt-in flag that is
  audited per row, like `--rescue-mnp`.

**Truth set.** The census harness at the 11 real twin loci
(`~/test/gbcms/harness/t9t10/`, local).
**Acceptance.** The census outcome is reproduced at all 11 loci, and no change
at non-twin loci. RNA: the transcript and ASJD columns agree with the reported
allele. Merge's mixed-winner check (M1) lands in the same PR.

### C4 — Reference windows near contig ends (#142, under #92) · M [counts]
**Finding.** The left-align wide window is not clamped to the contig length,
so an indel within about 100bp of a contig end is not left-aligned (a WARN
only). Prep's `ref_context` fetch fails the same way near the end, which is
one cause of `SW_FALLBACK`.
**Direction.** Clamp every reference fetch window to the contig length from
the `.fai` (left-align and `ref_context`).
**Tests.** Indels 10bp and 50bp from a contig end: they left-align, get a
context, and are scored by the pangenomic matrix, not SW.
**Acceptance.** No change away from contig ends (the RC set is
byte-identical).

**As built (2026-09-30).**
- `fetch_window` returns the part of `[start, end)` on the contig. Callers index
  from `start`, so only the end moves, and a short window means the contig end.
- It serves the five window fetches: left-align, `ref_context`, the adaptive
  repeat scan, the shift region and the event reference.
- Exact fetches keep `fetch_region` and still fail past the end: REF validation (a
  clamped REF could be "corrected" to a shorter one), the MAF anchor, the
  homopolymer next base and the splice motifs.
- **Found while testing: position dependence.** A complex variant near a contig end
  had no event reference, so the tolerant classifier that C1 replaced judged it.
  The same bases and reads gave 0/0/0 mid-contig and REF 5 / ALT 7 at the end, two
  reads carrying another allele among the ALT.
- The first count test had used reads ending inside the tract. Those reads showed
  C10's ALT side: a carrier ending inside the tract is credited from its CIGAR gap.
  Filed as C20 #188.
- **The adversarial review found two defects in the first clamp:**
  - Cost: prep ran 63× slower on GRCh38 with alt/decoy contigs, because every
    window looked the contig length up and the index lookup clones every record
    name. Now only a window that fails is clamped: 12.2 µs per indel on 3,366
    contigs, develop 12.1.
  - Aliases: a FASTA holding a contig under two names at different lengths
    clamped to the other record's length and moved an indel 59bp. The contig is
    now resolved by the name the fetch reads.
- The review also led to more test cases: an exact REF-validation guard, an alias
  FASTA, clipped carriers dropping out of depth (DP 5 vs 10 on develop), and
  informative reads in the position-independence test.
- Checked and left as is: the exact-carrier windows' padding at a contig end
  always fits the other margin (an event exhausting it is past LONG_EVENT).
- The review also found C21 #189: homopolymer insertion carriers placed elsewhere
  in the run are counted REF. It predates C4 and is reproduced mid-contig.

**Measured (2026-09-30).** All runs were local, develop vs the final branch head.
- **Prevalence:** among about 516,000 unique signed-out events, none lies within
  16kb of a contig end. The three that appear to are malformed rows whose
  `End_Position` exceeds the contig length. So for panel data the fix is
  defensive; it matters for WGS, chrM's control region and GRCh38 alt/decoy
  contigs.
- **RC set:** byte-identical, 33 FORTE truth samples, the T9 probes (3 samples)
  and 28 DNA runs (2,132 rows). Develop logged no window-fetch warning on any of
  them.
- **Cost:** prep on a 3,366-contig FASTA takes 12.2 µs per indel (develop 12.1).

### C5 — Long insertions exceed the pangenomic matrix cap (#120) · L [counts]
**Finding.** `MAX_HAP_LEN = 400` (`pangenome.rs`). ALT haplotypes longer than
that (insertions of roughly 300–390bp and up, depending on padding) cannot
use the matrix, so their reads are scored by SW (`SW_FALLBACK`).
**Measure first.** How many signed-out insertions are 250bp or longer, and
their `SW_FALLBACK` counts.
**Direction.** If demand exists, size the cap from the variant (insert length
+ padding), bounded by a documented memory limit.

### C6 — Error-tolerant exact-length insertion matching (#143, under #92) · L [counts]
**Finding.** Exact-length insertions whose bases confidently mismatch count as
`partial_alt`. Some are true ALT molecules with a sequencing error inside the
insert (observed at 1–6 reads on two long-insertion loci in the local data).
**Direction.** A backend-consistent identity band, like the truncation rule's
≥90% / non-low-complexity gates. The band exists for imperfect ALT
representations and BQ masking cannot replace it, so test both policies
(`.agents/memory/identity-band-annotation-tolerance.md`).
**Count-the-given-allele check (2026-09-25).** A tolerant band credits reads
whose insert confidently differs from the given ALT. That is accurate only
when the differences are read errors. Measure both the existing truncation band
and any new band at the loci where they admit reads: scattered mismatches mean
errors (count them); one recurring alternative insert means a different allele
or a wrong input (do not count it; flag what the reads carry).
**Acceptance.** The two loci recover their carriers under both backends, and
ladder and tract rows don't gain AD.

### C7 — Rescue for clip-borne ITD carriers (#144, under #92) · L [counts]
**Finding.** When every carrier of an insertion is represented as a soft clip,
the insertion counts no ALT. `CLIP_CANDIDATES(n)` (6.5.0) makes it visible.
**Direction.** The design filed on #92: validate clip sequence against the
insert's duplication. Rare: it matters only when all carriers are clips.

### C8 — One-base-REF delins without a shared anchor (#121) · L [counts]
**Finding.** Counting sends every one-base-REF variant to `check_insertion`,
including a delins such as `C>TA`, whose anchor base changes. This is
deliberate (the wrong-length rule covers anchor-substituting Ins+SNV), and
the 17-shape battery counts correctly.
**Direction.** Verify first: widen the battery (repeat flanks, reads with
errors, both backends). If any case miscounts, route these variants to
`check_complex`, as N×1 deletions with a substituted anchor already are.
**Found by C21's review (2026-09-30): a miscount.** `A>CCC` (in `A CC T`):
reads that keep the anchor A and insert CC at the junction count ALT (5 of 5 on
develop and the C21 branch), because the strict path compares only the inserted
bases. C21 already stops such an ALT matching a placement elsewhere; the strict
path needs the routing above.

### C9 — Count a MAF deletion at Start 1 (#122) · L [counts]
**Finding.** Preparation cannot anchor a MAF deletion at Start 1 (there is no
base before it), so the row is `FETCH_FAILED`. The VCF writer writes the VCF-spec
form with the base after the event instead.
**Direction.** Resolve such rows to the same base-after form in preparation
and count them. Telomere-only.

**Decision (2026-10-05, operator).** Base-after form by default (no flag): the
MAF output keeps the input row, so a result is looked up by it, and the row
counts as the same event given as VCF already does (REF reads; ALT cannot show at
a contig start). With it: VCF output of MAF input carries the MAF row
(`MAF_START`, `MAF_REF`, `MAF_ALT`) on every record. 0 such rows in the
sign-out dump (1.13M rows).

### C10 — Reads ending inside an indel's repeat tract are counted REF (#157) · H [counts]
**Finding** (while measuring C2). A read that starts or ends inside an indel's
repeat tract cannot show whether the tract carries the indel: its bases match
both haplotypes, so the aligner places no gap and IGV cannot tell. gbcms counts
it as REF.
- **Synthetic:** a 10-A homopolymer with a 1bp deletion. With 10 REF and 10 ALT
  reads spanning the tract, plus 10 of each ending mid-tract, gbcms reports VAF
  25% where the informative reads give 50%.
- **Real data:** 295 signed-out indel rows in the 6.5.0 RC runs. gbcms VAF ÷
  informative VAF has a median of 0.91 (STR) / 0.92 (homopolymer) / 0.97
  (unique), and 33 rows fall below 0.8.

**Standard.** GATK's AD counts only informative reads.
**Direction.** A read that doesn't span the indel's discrimination window
(its tract, or its span in unique sequence, ±1 base) counts toward depth only:
neither REF nor ALT. This is C2's window. Mirror it in the legacy parity oracle.
**Measure first.** Deltas on the RC set and the D5 panel; adjudicate a sample
of rows read by read in IGV. SNVs are unchanged.

**Implementation refinement (2026-09-25).** Two sharper definitions came out of
writing the tests; both keep the operator's rule.
- **The window uses the indel's shift-equivalence region**, not the tract from
  the repeat finder: the stretch the event slides over without changing the
  haplotype. In a homopolymer the two are the same. It also counts partial STR
  copies, and a one-base deletion inside a dinucleotide STR (which cannot
  slide) stays local instead of taking the whole STR.
- **A read is informative when it spans one side of the event**: a flank base
  through one base past the first base where REF and ALT differ, reading
  inward from that flank. The extra base is the aligner margin: one terminal
  mismatch costs less than a clip, so an ALT read ending on the differing base
  would otherwise look REF. For a homopolymer this is the tract ±1, as decided.
  A deletion longer than a read still gets REF reads from either junction; a
  whole-span rule would give it none.
- C2 keeps the region ±1: its question is whether another event sits where
  this row's alleles differ.
- The rule applies to pure indels only. Substitution-bearing events have no
  shift region, and their first base already discriminates.
- **Prep measures the region** (`Variant.shift_region`) over its own fetch sized
  to the event. The slide over `ref_context` was cut short for tandem
  duplications: the context is padded for 1–6bp motifs, and 12 of 82 RC
  insertion rows have longer regions (up to 58bp). A synthetic 30bp
  duplication read rd 30 instead of 10.

**Status (2026-09-25): implemented with C2 on one branch; RC acceptance done.**
Local data, aggregates only: the 6.5.0 RC runs, 1,060 DNA rows and 94 RNA rows.
- **Scope held.**
  - SNV and MNP rows: 0 of 825 changed.
  - `alt_count` and `total_count`: never moved.
  - `alt_count_fragment`: +1 to +2 on 3 rows, where a mate's vacuous REF call no longer contests the other mate's ALT.
- **Convergence** (gbcms VAF ÷ informative-read VAF, median):

  | Context | Before | After |
  |---|---|---|
  | STR | 0.910 | 0.983 |
  | Homopolymer | 0.919 | 0.994 |
  | Unique | 0.973 | 1.000 |

  Rows below 0.8 went from 33 to 20.
- **Outliers are census limits, not gate errors.**
  - Of the 13 rows above 1.2, 8 are tandem duplications whose true region (16–58bp) is far longer than the census's tract window. Re-censused over the true region, gbcms REF matches read for read (for example 295/295, 119/119, 1,275/1,275; all within 4%).
  - The other 5 are low-ALT rows that were already above 1.2: an ALT-side difference outside C10.
  - The lowest row is a 113bp deletion: the census demands a whole-span read, while gbcms takes REF from either junction, as decided.

### C11 — Phase-3 context misses tandem duplications (#159) · M [counts]
**Finding** (while landing C10). Prep pads `ref_context` from the 1–6bp motif
repeat span, capped at 50. A tandem duplication slides over its whole
duplicated segment, so the context often ends inside it:
- 12 of 82 insertion rows on the RC set, with true regions up to 58bp;
- a synthetic 30bp duplication got an 11-base context.

This breaks #91's rule that the haplotype window holds the whole tract. C10's
REF rule no longer depends on the context. ALT-side and Phase-3 classification
(the pangenomic matrix, WFA, S3, the AD-claiming windows) still do.
**Measure first.** On the RC ITD rows and the D5 FLT3-ITD stratum, compare
ALT/partial calls with the context padded to cover `shift_region` against
today's. Watch `MAX_HAP_LEN` (C5).
**Direction** (if calls change): pad `ref_context` to the shift region plus
flank, bounded by the matrix cap.

### C14 — Records without bases (SEQ `*`) (#172) · M [counts]
**Finding.** A BAM record with no sequence reached the classifiers.
- The SNV check indexed the empty sequence and panicked.
- Heuristic BAQ sliced its empty qualities and panicked (RNA, or DNA with
  `--apply-baq`).
- The insertion and deletion checks counted it REF and depth, or a fragment when
  it was a kept secondary, from its CIGAR alone.
- The MNP and complex checks already skipped it.

**As built (2026-09-29).** Both counting paths send every fetched record through
the shared `ReadFilter::passes`, which now drops a record without bases after the
flag filters. The record counts in neither depth, fragments nor `mq0_count`. The
per-bin tallies are summed across bins and warned once per counting pass. The
`--rescue-mnp` re-count is a second pass and warns again; O1 (#130) tracks the
same quirk for the UMI warning.

**Measured (2026-09-29).** All runs were local, develop vs the branch.
- At defaults, every RC input is byte-identical: 33 FORTE truth samples, the T9
  probes (3 samples) and 28 DNA runs.
- DNA with `--no-filter-secondary`: 28/28 byte-identical, and develop did not
  panic.
- No input held a record without bases. The fix is defensive in MSK data, like
  the QC-fail filter.

**Found by the review, not in scope.** Both were reproduced and filed:
- C19 (#182): absent QUAL (0xFF) overflows the u8 quality margin in fragment
  consensus. It panics in debug builds and wraps in release.
- O7 (#183): unmapped-flag (0x4) mates placed at a variant count in
  `mq0_count`.
- O1 (#130) now also covers the no-bases warning repeating in the rescue
  re-count.

### C21 — Indel carriers written elsewhere in their repeat (#189) · H [counts]
**Finding.** The windowed checks accepted a shifted indel placement by a proxy,
not by its haplotype:
- an insertion when the reference base before it equalled the anchor base: never
  true inside the repeat, so carriers written elsewhere counted REF (synthetic
  `G AAAAA T`: REF 10 / ALT 0 instead of 5 / 5), and true by chance for some
  placements of another haplotype, which counted ALT;
- a deletion when its removed bases equalled the given ones: a rotated STR
  placement under 5bp counted REF, the same bases deleted outside the repeat ALT.

**Measured first (2026-09-30).** Placement of repeat indels, near the RC set's
indel rows:
- DNA panels (realigned): 37 of 15,047 insertions and 1 of 12,244 deletions not
  left-aligned; 2 of 82 insertion rows have shifted equivalent carriers (87 reads,
  against 13,355 left-placed).
- FORTE RNA (STAR): 162 of 170 insertions not left-aligned, but the truth set has
  3 insertion rows and none has carriers; the RNA probes (under Measured) test
  STAR's placements directly.
- WES without realignment (BWA-MEM): 0 of 620 insertions, 0 of 1,367 deletions.

**As built (2026-09-30).**
- S3 is now the haplotype: `same_insertion_haplotype` (`X + S = S + Y`, or
  `S + X = Y + S` left of the junction; never for an anchor-substituting ALT) and
  `same_deletion_haplotype` (the stretch between the placements repeats with
  period `len`), read from the event reference (`event_ref`), else `ref_context`.
- The placement must be the read's only change across the discrimination window
  (`only_change_in_window`, both sides): another gap, insertion or splice there
  makes a distinct allele. Window bases past the read's end are not required,
  as on the strict path (C20 owns that rule).
- `window::scan_window` reaches every placement: every junction of an insertion,
  a deletion's starts up to `hi - len`, besides `max(5, repeat_span + 2)`
  (`repeat_span` counts motifs of up to 6 bases, so a longer duplication slid past
  it). The splice triage and the tract-cluster grouping pad use the same reach.
- The variant's inserted bases at a non-equivalent junction: a distinct allele
  inside the discrimination window, REF (a separate event) outside it, as before.
- The backward-boundary check (an insertion right before the anchor base) now
  requires the same haplotype; it credited the variant's bases there as ALT.
- The deletion flag `has_wrong_length_nearby` is now `has_distinct_allele_nearby`.

**Two adversarial reviews.** The first (on the first fix) found no check that the
read aligns the stretch between the placements, anchor-substituting ALTs matched
as pure insertions, REF withdrawn from reads whose same-bases insertion lies
outside the window, the backward-boundary and window gaps, and lost coverage of
the AD-claiming guard (a guard now covers its Test 1, an unreadable insert at a
multi-allelic site). The second (on the follow-up) found the deletion scan
reaching into its own deleted span (a different same-length deletion there went to
Phase 3 and came back ALT), and a read-side check too narrow on insertions and
absent on deletions (a +AA read, a cancelled deletion, a split −4 read and a
deletion after a splice over the anchor counted ALT). All fixed red-first. Two
RNA splice contract tests pinned ALT for reads spliced over the anchor with a D
written after the junction (their bases are those of a reference read spliced two
bases later); they now count toward depth only.

**Measured (2026-09-30).** All runs were local, develop vs the final branch
head, every changed read adjudicated by its own bases (read where they sit between
the aligned flank bases, bases below Q20 matching anything).
- DNA (28 runs, 1,060 rows): 13 rows change; ALT +46, partial −62, REF +13,
  depth unchanged.
  - Two insertion rows gain 44 and 10 ALT: carriers written 1–6 junctions right of
    the left-aligned position, each holding the ALT haplotype. Their last inserted
    base is an N (Q5), which failed develop's base check; 38 of them lay past the
    old window. One read with a 6bp insertion in the 7bp row's region becomes
    partial.
  - A BRCA2 cluster (2bp, 14bp and 33bp deletion rows, two samples): the 2bp rows
    move 6 and 18 reads partial → REF. They delete the same two bases 6bp before
    the anchor, another haplotype outside the window, which develop matched by
    its bases and the sibling guard then demoted. The 33bp rows regain the same
    reads as REF, which develop's false match had claimed for the 2bp sibling.
  - Seven 1bp indel rows lose 8 ALT reads (to REF or neither) that delete or insert
    the same base outside the run; one ends before the anchor.
  - Every read that gained ALT holds the ALT haplotype; none that lost ALT does.
- FORTE RNA: the 33 truth samples and the T9 probes are byte-identical.
- RNA probes of STAR's repeat insertions near the truth loci (23 probes): ALT 11 →
  183, REF 22,102 → 21,929. An independent census counts 176 shifted equivalent
  sense carriers; all are ALT on the branch (4 were ALT on develop too, where the
  anchor-base proxy happened to hold).
- The first follow-up build, which scanned the whole deleted span, moved 52 and 123
  reads at the 33bp rows from REF to partial (other deletions inside the span);
  the final reach removes that.

**Found by the reviews, not in scope:**
- C22 (#191): a same-length deletion ≥5bp that gives another haplotype still goes
  to Phase 3, which calls it ALT (70/160 at 5bp, 150/160 at 8bp, synthetic).
- C23 (#192): "in a repeat" for a distinct allele is `repeat_span >= 2`, so events
  in repeats with motifs over 6bp keep REF + partial.
- C20 (#188) gains the windowed paths (a shifted carrier ending inside the run is
  credited ALT from its gap) and the strict path's lack of an only-change check (a
  read with the indel at the junction plus another change in the window counts
  ALT, develop too).
- C8 (#121): the strict path credits an anchor-substituting insertion (`A>CCC`)
  from its inserted bases alone.

### Cluster 1 — pure-indel read judgment: C20 (#188), C22 (#191), C23 (#192), C8 (#121); C11 (#159) closed · [counts] [done]
One rule, one PR: a carrier holds the ALT haplotype across the window, by its own
bases.

**Measured first (2026-09-30).** Develop, every read call traced and judged by its
own bases, on the RC panels (295 pure-indel rows, 42,594 ALT reads), WES without
realignment (80 loci from paired IMPACT/TEMPO libraries, sliced on the cluster)
and FORTE RNA (18 pure-indel rows):
- C20: 25 RC ALT reads credited without discriminating bases (13 fit both
  alleles, 12 align neither flank), 0 WES, 0 RNA. C10's reference-coordinate
  windows must not be mirrored onto ALT reads: 2,738 RC carriers (mostly long
  insertions) span neither window yet their bases discriminate.
- C20 (strict path, another indel in the window): 2 RC reads. C22: 6 RC reads,
  all another allele by their bases. C23: 0 reads (21 RC, 17 WES, 4 RNA rows
  qualify). C11: 0 ALT calls change (17 RC, 14 WES, 1 RNA rows) — closed with the
  measurement.
- C8: 769 anchor-changing one-base-REF rows in the 1.13M signed-out set; at 10
  WES loci 156 of 463 ALT reads (34%) keep the REF anchor (another allele).

**Decisions (operator).**
- Scope (2026-09-30): C8 plus the small fixes in one PR; close C11.
- No margin base on the ALT side (2026-10-01): a non-informative ALT read keeps
  ALT when it reads, unmasked, a base where the alleles differ. C10's margin
  guards CIGAR-only REF calls against a hidden terminal mismatch; this check reads
  the deciding base itself, on a read the CIGAR already calls ALT.
- C8 (2026-10-01): accept the exact-carrier rule's semantics for anchor-changing
  variants. Investigated first: the reads that lose REF either do not carry the
  REF anchor (another allele — 52 at one WES locus carry the anchor change without
  the insertion) or do not hold the rule's window, which grows through repeats and
  is padded to the ALT's length so neither allele is favoured by where reads end
  (a sharp line: withdrawn reads hold ≤6 bases past an 8-base event, kept ones ≥7).

**As built.**
- `alt_needs_the_window` + `window::alt_bases_discriminate`: the ALT side of C10,
  by bases (reference two bases past the window, soft clips not read, masked bases
  skipped). `ClassifyResult::ref_uninformative` → `uninformative`.
- Strict path: an indel at the junction is ALT only with no other I/D across the
  window (below 50bp; at ≥50bp the large-deletion band's tolerance applies).
- C22: a same-length deletion that fails the placement test counts ALT when its
  bases spell the ALT across the window (`window::read_spells_alt`: anchored on the
  read's nearest aligned bases outside, and readable wherever the haplotype its own
  CIGAR proposes differs from the ALT); otherwise, at ≥5bp, a distinct allele.
- A distinct allele is neither + partial where the event slides (`window::slides`)
  or the read has an indel inside the discrimination window; REF + partial only in
  unique sequence with the window clear (C23).
- C8: one-base-REF ALTs that change the anchor go to the exact-carrier rule.

**Reviews and the read-level check.** The first review found the C22 decision
ignoring bases, the ALT rule skipping events ≥58bp, masked-base and clip handling,
and the band asymmetry; all fixed red-first. Rendering real reads then showed the
spell-the-ALT check accepting an N where two deletion alleles differ (a CG>C row
whose reads deleted the C before the run) — fixed red-first. The second review
(on the final three fixes) found, all fixed red-first:
- deletion carriers ending inside the region counted ALT: a carrier's extent
  counts its own gap, so it spans C10's reference windows with up to L−1 fewer
  bases than it needs. The ALT side now uses the windows read on the ALT
  haplotype (a deletion's are an insertion's), then the bases;
- the walk saw an I/D only right after an aligned block, so `M I D M` for a
  deletion row and `M D I M` for an insertion row counted REF while the same pair
  written the other way was judged;
- the ≥50bp band stopped short of the discrimination window for a deletion
  sliding through a long repeat (a D(60) with an insertion 80 bases on counted
  ALT);
- a base past the REF stretch fetched for the ALT check counted as a difference;
- two test fixtures built homopolymers instead of random sequence.
It also found, pre-existing, that a read with a short indel inside the window
counts REF, and that cluster 1's in-window rule flips such reads to neither when
an unrelated op appears: a decision, filed as C26 #200. Not defects under the
rule (noted): a truncated +AA read written `I(A)`+clip vs `I(AA)` is decided by
placement between +A and +AA; a shifted exact-length ≥50bp deletion does not get
the band's tolerance.

A third review (on the second-review fixes) found three defects in them, all
fixed red-first:
- a deletion followed by a same-length insertion of other bases went to Phase 3,
  which called it ALT; its length change is never the ALT's;
- from 58bp an insertion's REF stretch stopped short, because prep's event
  reference kept a 60bp margin, so carriers masked past the first difference
  lost ALT. The margin now reaches the insert's length plus three;
- a deleted anchor followed by an insertion became partial whatever it held; its
  bases now decide.
Its fuzz of 480 random I/D/X pairs near the anchor found no other new false ALT,
and binned and legacy agreed on every probe. Pre-existing, filed:
- reads spelling the ALT across several indel ops (`D1 D1` for −AA, shifted
  pairs), as C27 #201 (6.7.0);
- the ALT-side window counting a splice N, noted on R5 #198.

The reviewer's fuzz, rerun on the final head over 2,400 reads (five seeds), moved
only reads whose bases are another allele (to partial), with no new false ALT.
Its 130 ALT calls on reads whose bases are not exactly the ALT are all ALT on
develop too:
- 42 carry the ALT plus mismatches, which is correct;
- 88 delete the anchor base and fall back to Phase 3, whose PairHMM credits the
  closer of two haplotypes. Filed as C28 #202 (measure first).

**Acceptance (final head, every changed read adjudicated by its own bases;
harness `~/test/gbcms/harness/c20/`, local).**
- RC DNA: 15 of 1,060 rows change (ALT −39, REF −107, partial −15, depth 0).
  - ALT withdrawn: carriers ending inside the repeat or before the deciding
    base; reads carrying other deletions (24–54bp, or a 33bp one 15 bases off)
    that Phase 3 had called ALT for a 33bp deletion; 6 reads at a 6bp deletion
    whose bases fit both alleles.
  - REF −118 at the two anchor-changing rows (C8).
  - Two consequences that are not improvements, both from rules that predate
    this PR:
    - At a co-annotated +AAG insertion, 16 reads move from partial to REF. They
      delete across the anchor (10) or carry a D1 in the window (6), and only the
      33bp row's false ALT had kept them out of REF. Recorded on C28 #202 and
      C26 #200.
    - At a 6bp deletion, 5 reads whose bases are REF across the window lose REF
      to partial. A same-length deletion of other bases lies outside the window,
      and C22 flags it wherever it lies in the scan window. C26 #200 now asks the
      question in both directions.
- WES: 10 of 80 rows change, all C8 (ALT −12, REF −399, partial +152, depth +62).
  At one locus 146 reads keep the REF anchor and carry only the insertion, so
  they count partial.
- FORTE RNA: the truth set and the STAR repeat-insertion probes are unchanged. One
  T9 probe (a 4bp deletion) loses 2 ALT reads, which carry a 3bp deletion ending
  just after it.

## Test architecture

### T1 (#170) + T2 (#171) — design note (2026-10-01) · [approved] [in review]
**Facts.**
- **The legacy path.** `count_bam` (200 lines) and `count_single_variant` (501)
  sit behind the `legacy-parity` feature: on by default, off in the shipped
  wheel. 283 of 850 tests call it. 234 compare legacy with binned; 49 run legacy
  only, so they never test production.
- **What parity caught.**
  - Its one bin/fetch bug, CR-1, was found by code review. The synthetic parity
    suite could not have caught it: every synthetic contig (≤2kb) fits in one
    10kb bin.
  - The two divergences between the loops (mq0 tallied before vs after the RNA
    strand filter; the mFSD N heuristic) were also found in review, in fields
    parity does not compare.
  - Every classification bug this cycle sat in the shared classifier: CR-2, CR-4,
    #91, C12, C21, #166, cluster 1. Parity cannot see those.
  - It is also blind to BAQ, UMI, RNA, mFSD, ASJD and siblings.
- **Cost.** 35 commits edited `count_single_variant` twice, 24 of them since the
  June review. Its tests take 0.07s of the suite's 13.4s, and the extra CI steps
  about 4s, so retiring it saves the double edits, not time.
- **Binning.**
  - There is no fixed grid. A bin starts at the leftmost unbinned variant and
    ends at max(start + window, the anchor's span end). Variants join while they
    start before the end and the bin holds fewer than the cap.
  - The fetch pads by max(5, repeat_span + 2), and each variant filters the
    cached reads to its own window.
  - `build_genomic_bins` already takes the window. Window 1 with cap 1 reproduces
    the legacy fetch exactly, through the production loop.
  - One variant per call works today; MNP rescue relies on it. Only the BH
    q-values depend on which rows are in a call.

**Design.**
1. **A test hook.**
   - Two optional keyword arguments on `count_bam_binned` and
     `count_bam_binned_observations`: `bin_window` and `bin_max_variants`.
   - `None` means the production constants; a value below 1 raises `ValueError`.
   - The pipeline never passes them. The stub carries them.
2. **One helper replaces the parity helpers.** `count_both` / `count_one_both`
   become `count_checked`.
   - It runs production geometry and window 1 / cap 1 (the legacy fetch, now
     through the production loop), and compares every field: integers exact,
     floats to 1e-9 relative, NaN equal to NaN.
   - Siblings, RNA, BAQ, UMI, mFSD and ASJD are now covered, because there is
     one loop.
   - The 49 legacy-only tests move to production through the same helper, and
     `test_accuracy`'s four legacy/binned twins merge.
3. **Binning invariance**, in `tests/test_binning_invariance.py` (replacing
   `test_parity_large_deletion.py`). Counts must agree under these geometries:
   - production;
   - window 1 / cap 1;
   - a tiny window (7–20bp), where anchors longer than the window reproduce CR-1;
   - cap 2–3;
   - one call per row (q-values excepted);
   - shuffled input with 4 threads;
   - decoy variants that move bin starts (q-values excepted).

   Fixtures:
   - a dense SNV/MNP/indel/delins cluster with reads starting at every offset;
   - a 50–60bp deletion anchor with an SNV inside, and reads past its
     breakpoint;
   - a long repeat tract on a non-anchor member;
   - contig start and end;
   - siblings and a decomposed twin;
   - RNA with a GTF, strandedness, ASJD, BAQ, mFSD and a UMI tag;
   - the repo's test BAM with a ~1kb window, so it splits into bins.

   Sorted observation rows must match too. A Rust property test on
   `build_genomic_bins` checks that every variant lands in exactly one bin and
   every bin's fetch holds each member's window. As a mutation check,
   re-introduce the CR-1 bug locally.

   The mutation check (2026-10-01) found that BAM-level counts can show CR-1 in
   only one read shape. Every counted read overlaps the event's first base, which
   every fetch holds, except a DNA read admitted by its soft-clipped bases. So
   the matrix includes left-clipped carriers aligned only after a 60bp delins:
   with CR-1 back, AD drops from 20 to 0 under window 1. The Rust property test
   fails too, on the geometry itself.
4. **The read census (T2)**, in `tests/census.py`. It ports the harness's
   per-read judge:
   - bases below `min_baseq` are masked;
   - it anchors on the read's aligned bases just outside the event (pure indels:
     the shift region plus a three-base margin; anchor-changing one-base REF:
     two flank bases each side);
   - soft-clipped bases are read where they sit;
   - reading stops at a splice N (SPLICED when the N covers every deciding
     base).

   Verdicts: ALT, REF, FITS_BOTH, CONTRADICTS_BOTH, NOT_ANCHORED, SPLICED.
   `census(bam, ref, variant, ...)` gives per-read verdicts and REF/ALT/neither
   counts with fragment bounds. `assert_matches(counts, census)` checks RD and
   AD exactly, DP, the fragment bounds and the four counting invariants.

   It is used in three ways:
   - on the binning-invariance fixtures;
   - in a property test over generated pure-indel reads, where engine RD/AD must
     equal the census for the settled read shapes (one indel plus mismatches).
     The open decisions become strict xfails: C26 #200, C28 #202 and C27 #201;
   - by new contract tests, with existing ones moving over when touched.

   As a mutation check, revert one cluster-1 fix locally; the property test must
   fail.
5. **Delete:**
   - `count_bam` and `count_single_variant`, and the feature in `Cargo.toml`,
     `mod.rs` and `lib.rs`;
   - their stub entries;
   - the helpers `count_one`, `count_one_both` and `count_both`
     (`PARITY_FIELDS` → `COUNT_FIELDS`);
   - the two `TestBinnedParity` classes;
   - the duplicate CI clippy and `cargo test` steps, and `--no-default-features`
     in `release.yml` and the Dockerfile;
   - the parity memories, marked superseded;
   - doc lines in AGENTS.md, the rules, skills, CONTRIBUTING and the
     developer/testing guides.
6. **AGENTS.md invariant 1 becomes binning invariance.** Bin geometry is
   performance only. `count_bam_binned` must give identical counts under any bin
   window or cap, one variant per bin, or one variant per call, so every bin's
   fetch must hold each member's full window. Classification is checked against
   the read census, never a second engine.

**Order.**
1. Add the hook, the invariance tests and the census, green on the current code
   and mutation-checked.
2. Migrate the helpers.
3. Delete the legacy path.
4. Update the docs.

Production output must stay byte-identical: the engine does not change, and the
RC/WES/FORTE acceptance must show 0 changed rows.

**As built (2026-10-01; approved by the operator as written).**
- Step 1: the hook, the Rust bin property test, and
  `tests/test_binning_invariance.py`. Mutation-checked: re-introducing CR-1 fails
  the property test and the clip-carrier fixture.
- Step 2: the census (`tests/census.py`, `tests/test_read_census.py`).
  Mutation-checked: disabling the ALT-side window, or equivalent insertion
  placements, fails it. Writing it settled the census's conventions to match the
  decided rules: the anchors are the tract's flank bases; a one-sided reading runs
  on along each haplotype past the window; wrong-length reads run past where they
  differ from the ALT, because a two-haplotype census cannot see a third allele in
  a truncated read.
- Step 3: `count_checked` everywhere. The 49 legacy-only tests now test
  production, and all hold. Duplicate twins were removed.
- Step 4: the legacy path deleted (841 lines out, 59 in), along with the feature,
  the stub entry, the CI steps and `--no-default-features`.
- Step 5: docs, rules, skills and memory. The two parity memories were deleted,
  since the repo now records the rule.

**Found on the way, for H3 #204:** `mfsd_ref_llr` varies in its last ulps
between identical runs (float summation in HashMap order); output rounds it to
4 decimals.

## RNA


### C12 — Carriers whose allele lies in soft-clipped bases (complex variants) (#167) · H [counts] [decided]
**Finding.** A read enters DP, RD and AD only if its *aligned* span covers VCF POS
(the anchor-overlap gate, both counting loops in `engine.rs`). Aligners soft-clip
ALT reads whose allele sits within a few bases of a read end; REF reads at the same
positions align fully. So clipped ALT carriers drop out while the REF reads beside
them count, and VAF reads low. The exact-carrier rule already reads clipped bases
from an aligned anchor: at one repeat locus it called 9 clipped reads ALT, and the
gate then dropped them.

**Measured (local data; aggregates).** Signed-out IMPACT complex variants, and the
same libraries recaptured on WES by TEMPO (bwa, no indel realignment):
- Synthetic, 50% sample, ALT reads clipped when fewer than K bases lie past the
  event: VAF 0.45–0.46 at K = 10, 0.37–0.44 at K = 20 (0.50 at K = 0).
- IMPACT (ABRA2-realigned): 19 exact-ALT carriers outside depth over 40 variants
  in one set; 270 over the 40 paired loci, 202 of them at one 18>8 delins.
- WES, same 40 loci: 420 exact-ALT carriers outside depth against 2 REF; pooled VAF
  0.174 → 0.204 if counted. WES − IMPACT VAF, median −0.021 → −0.007 with them.
- Long deletions (25–106 bp): 0 ALT reads on WES at all five, against 5–291 on
  IMPACT. Without realignment their carriers are clipped or split, and the
  measurement above (whole windows in the clip) could not see them either.
- C1 itself agrees with a position-aware read census on both platforms
  (median |VAF difference| 0.002 IMPACT, 0.000 WES).

**Proposal (for review).**
1. *Admission.* For variants the exact-carrier rule judges (delins, Del+SNV, MNP
   reads with an indel at the block), a read whose aligned span does not cover POS
   is admitted when the rule **decides** it, REF or ALT, from its own bases with at
   least one window anchor aligned. It then counts in DP, RD/AD and fragments like
   any read. A clipped read the rule cannot decide stays out, as today, so DP gains
   only informative reads. REF and ALT are admitted by the same test.
2. *Long events.* Junction windows read on into the clip, so a split read's primary
   alignment, whose clip holds the far side of the deletion, is judged at the
   junction it shows. Supplementary alignments stay filtered: they are the same
   molecule.
3. *Guards.* Clipped bases past the fragment end (|TLEN| shorter than the read) are
   adapter, not allele: masked. Reads without a defined fragment (unpaired, mate
   unmapped, TLEN 0) are not admitted by their clips (as GATK does). Masked clipped
   bases match both windows, so an all-low-quality clip cannot decide a read.
4. *Anchor in the clip.* If an aligner clips before the flank, no anchor is aligned
   and the read stays out. Extrapolating the alignment into the clip is a later
   step, taken only if the acceptance shows those reads matter.
5. *Scope.* DNA only at first: STAR clips short junction overhangs, whose bases come
   from the next exon (measure on FORTE before enabling RNA). Pure insertions stay
   with C7 (#144); pure deletions and SNVs are unchanged.
6. *Both counting paths* (binned and the legacy parity oracle) and the
   per-transcript path apply the same admission. No new column: a trace line per
   admitted read and a debug count per variant.

**Acceptance.** The truth is per BAM and per read, not cross-platform agreement.
- *Per read:* every admitted clipped read carries exactly the allele it is counted
  for, by a position-aware census that reads clipped bases and uses the rule's
  repeat growth and equal-length windows; no REF↔ALT flip.
- *Symmetry:* synthetic uniform-start samples with ALT reads clipped at K = 10 and
  20 read 0.50 VAF; REF reads clipped the same way are admitted alike.
- *Same library, secondary:* IMPACT and WES are expected to move toward each other,
  not to match. They legitimately differ by bait design (a long deletion removes
  target sequence, so capture efficiency differs by allele), depth and duplicate
  rates (sampling noise), sequencing run and instrument (base-quality profiles),
  aligner versions, ABRA2 and BQSR, and MAPQ distributions. Report each locus's
  WES − IMPACT gap against its binomial interval from the two depths, and
  adjudicate read by read only the loci outside it after the change.
- *Regression:* RC set, only complex rows (and their grouped rows) move; FORTE
  unchanged (RNA off); run time within noise.

**Implemented and accepted (2026-09-28, `feature/c12-clip-carriers`; aggregates).**
- 40 paired loci (signed-out IMPACT complex variants; the same libraries on TEMPO
  WES): WES − IMPACT VAF median −0.021 → −0.009, mean −0.043 → −0.018. Long
  deletions return on WES (79>4: 0 → 139 ALT reads, VAF 0 → 0.174; 106>7: 0 → 12;
  68>3: 0 → 3); IMPACT unchanged there.
- Every admitted read is ALT (none REF on these loci). An independent anchored
  check (the read's own bases from the aligned edge of the tract ±2 into its clip)
  confirms 233 on IMPACT and 365 on WES; the rest are reads that cannot hold the
  check's longer stretch, and one WES read differs one base beyond the 2-base flank.
- A defect the first run exposed, fixed red-first: unclipped REF reads starting
  inside a long deletion were admitted (IMPACT VAF fell 0.324 → 0.213 at 79>4).
  Admission now requires a window that reads the read's clipped bases.
- RC set (28 DNA runs): 11 of ~1,060 rows move, 10 MNP (clipped carriers at the
  block, +1 to +10 ALT) and 1 complex; no SNV or indel row. An independent count of
  exact clipped carriers is never below C12's change (C12 admits fewer where
  windows run past the fragment end: adapter).

**Decisions (2026-09-28, operator): as proposed.** (a) DP gains only decided
clipped reads; (b) RNA off until measured on FORTE; (c) admitted reads are named in
the trace log with a per-variant debug count, no new column or flag.
### R1 — Span-aware exon-edge BAQ rule (#106, 6.5.0 § T9) · L [counts] [decided]
**Finding.** The exon-edge BAQ exception keys on the variant's first base, so
a multi-base variant reaching an exon's right edge keeps BAQ. MNP counts move
at most 0.17%, and deletions are unaffected.
**Decision.** Should `exon_boundary_dist` become the span distance too, so the
rule and the column share one definition? Recommended: yes. It changes the
column only for multi-base variants.
**Direction.** Option A of the 6.5.0 plan: distance to the REF span, 0 when a
boundary lies inside it. Limit the "left edges are safe" guard to spans that
start at or after the exon start.
**Acceptance.** SNV rows byte-identical; multi-base edge probes change as
designed; affected sign-out rows adjudicated read by read.

**Decision (2026-09-25, operator).** Span-based distance (option A), VEP-style: 0 when a boundary lies inside the REF span. `exon_boundary_dist` takes the same definition, and the change is documented as a column change. Evidence on #106: 41% of signed-out indels change distance, and 1 of the 94 RNA truth rows crosses the 5bp window.

**As built (2026-09-29).**
- `nearest_splice_distance(chrom, first, last)` returns the least distance from
  any base of the span, 0 when a boundary lies inside it. The engine measures the
  prepared variant's REF, which is left-aligned, and a pure deletion's REF starts
  at its anchor base. BAQ and `exon_boundary_dist` both use this value.
- The adversarial review found that a `--rescue-mnp` component re-count resolved
  the rule on its own base. A far component therefore kept BAQ while the row
  skipped it, and a rescued row reported the component's distance. Now each
  component is counted over its MNP's span (`Variant.boundary_span`).
- Left edges need no guard: the span distance covers a span that starts in the
  intron and reaches an exon's first bases.

**Measured (2026-09-29).** All runs were local, one BAM at a time, develop vs the
branch. Only aggregates are recorded here.
- **FORTE truth cohort** (33 samples, 94 rows): all 68 SNV rows are
  byte-identical. The column changes on 8 multi-base rows, and 1 row flips the
  window with unchanged counts. No count changes.
- **T9 exon-edge probes** (326 probes × 3 samples):
  - 139 MNP rows move, all where the window flips (none elsewhere). No deletion
    moves.
  - Summed changes: RD +489, partial +6, AD −2; depth unchanged.
  - Per read: in 135 of 139 rows, the change equals what the reads' own bases
    predict, read at sequencer quality vs under emulated BAQ. In the other 4 the
    census tally is off by one read. Every read gbcms re-called there was
    inspected, and each new call is what its own bases say:
    - mostly low quality → REF: under BAQ, the only readable base (a Q24 REF base
      4bp from the junction) fell to Q4;
    - one REF → neither: a Q24 third-allele base that BAQ had hidden.
- **Column:** 978 of 978 probe rows equal an independent span distance of the
  engine-normalized variant. 61 probe deletions are left-aligned by prep.
- **`--rescue-mnp` on the probes:** 4 of 243 MNP rows change outcome.
  - 1 goes from rescued to `haplotype_confirmed`: once the edge base is read, a
    read shows the whole MNP.
  - 3 are newly rescued on one-read components, which is rescue's existing
    behaviour once the edge base is read.
- **Final code:** byte-identical to the first fix on all 36 RNA inputs at
  defaults.
- **DNA** (28 runs, 1,060 rows): byte-identical.

### R2 — RNA strandedness gating observability (#114) · M [decided]
**Finding.** At RNA defaults (strandedness enforced), antisense reads are
filtered before the sense/antisense tally. So `rna_antisense_depth` is always
0, and ASJD's `STRAND_DISCORDANT` is effectively unreachable.
**Decision.**
1. Should `rna_antisense_depth` count antisense reads before the filter, as
   `mq0_count` already does? Or be documented as "reads that reached
   classification"?
2. Should `STRAND_DISCORDANT` measure all junction reads, or be documented as
   a `--no-strandedness` diagnostic?

**Acceptance.** A red test at defaults for the chosen behaviour, and the doc
note in output-formats; counts unchanged.

**Decision (2026-09-25, operator).**
1. `rna_antisense_depth` counts antisense reads before the strand filter, like `mq0_count`. REF/ALT are unchanged.
2. `STRAND_DISCORDANT` is documented as a `--no-strandedness` diagnostic.

Evidence on #114: 143 antisense reads (0.5%) at 18% of truth loci are reported as 0 today; antisense junction reads are 0.1%.

**As built (2026-09-29).**
- Under enforcement the strand filter marks an antisense read instead of skipping
  it. A read that is not first-class, or not over the anchor, stops there. Any other
  is classified as a sense read would be (MAPQ rule, BAQ, classification, sibling
  guards). When it is REF or ALT it is tallied in `rna_antisense_depth`, the sense
  tally's definition. It is then dropped before any count, fragment, mFSD, distance
  or observation.
- The column therefore means the same with `--no-strandedness`. It counts REF/ALT
  reads, not every read at the locus, so it reports fewer than #114's 143 raw
  antisense reads.
- Binned path only: strandedness is parity-exempt, and the legacy comment is
  updated.
- `STRAND_DISCORDANT` is documented as a `--no-strandedness` diagnostic where the
  gene strand is resolved.
- The adversarial review found no counting defect; its oracle probe matched every
  column. It led to:
  - an oracle test per variant type;
  - doc precision: the column's exceptions, ANT descriptions, the strand table
    inverted for `reverse`, and "gene strand from the MAF" corrected to the GTF;
  - two filed follow-ups: R4 #185 (no gene strand at intronic loci and
    opposite-strand overlaps, so enforcement is a no-op there; reproduced) and
    O8 #186 (the observed-allele diagnostic reads antisense reads and skips the
    NH rescue).

**Measured (2026-09-29).** All runs were local, one BAM at a time, develop vs the
branch.
- **FORTE truth cohort** (33 samples, 94 rows):
  - no column changes except `rna_antisense_depth`;
  - the column goes from 0 to 203 reads, on 15 rows;
  - it equals its `--no-strandedness` value on 94/94 rows.
- **Independent per-read census** of the 68 truth SNV rows (the reads' own bases,
  strand and MAPQ rules, BAQ emulated): 68/68 rows exact, 181 = 181 reads.
- **T9 probes** (978 rows): no other column changes. The column goes from 0 to
  3,678 reads (0.12% of sense depth, on 393 rows) and equals `--no-strandedness`
  on 978/978.
- **DNA** (28 runs): byte-identical.
- **Final head** (with the review's early skip) vs the tested build: byte-identical
  on every input, both modes and DNA.

## Input and representation

### I1 — MAF allele base check (#123) · M [decided]
**Finding.** `MafReader` does not check allele bases. A MAF ALT with an IUPAC
code (e.g. `R`) passes preparation and counts 0 ALT silently. The VCF reader
skips such alleles since 6.5.0. The same check closes a cosmetic
Python/Rust label mismatch on non-ASCII alleles.
**Decision.** A `FAIL` row with a new reason (e.g. `NON_SEQUENCE_ALLELE`) —
recommended, because MAF→MAF output must keep every input row — or a skip
with a WARN like VCF.
**Tests.** MAF rows with `R`, `.`, and lowercase (valid) alleles.

**Decision (2026-09-25, operator).** A visible `FAIL` row with a reason, so MAF→MAF output keeps every row. Lowercase bases are valid. Evidence on #123: the automated pipeline's 19.8M calls have 0 invalid alleles; the 9 in the 1.13M curated sign-out rows come from hand edits.

**Decision (2026-10-05, operator).** VCF output writes such a row as the
symbolic record `<NON_SEQUENCE>` (REF the reference base at POS, declared in a
`##ALT` line): it keeps every row and stays valid VCF. Re-measured: 7 invalid
alleles in the current dump.

### I2 — `End_Position` optional (#124) · L
**Finding.** `MafReader` requires an integer `End_Position` but nothing uses
it; rows without one are skipped with a WARN. maf2vcf converts them.
**Direction.** Parse it when present; don't require it. Update the
required-columns table.

**Decision (2026-10-05, operator).** Accept rows without it (as maf2vcf does).
The sign-out dump has it, consistent with Start and REF, on all but one of 1.13M
rows.

**Merge (2026-10-05, operator).** End_Position is not part of a variant's
identity (genotyping never reads it), so `gbcms merge` joins on contig, Start and
alleles, fills End_Position from the inputs that have each row (a later-only row
would otherwise lose it), and logs an input that writes it differently. No effect
on binning or counts (merge runs on outputs); lookup by the input is kept (each
row keeps its End_Position). Measured: no two sign-out rows (per sample) differ
only in End_Position.

### I6 — `REF_MISMATCH` rows say where the given REF sits (#218) · L [decided]
**Finding.** REF validation compares only at the stated position (tolerant
validation, 2.6.1: ≥90% → `WARN_REF_CORRECTED`, FASTA REF used; below →
`REF_MISMATCH`, FAIL, zero counts). Nothing has looked elsewhere. Sign-out dump
(1.13M rows): 19 `WARN_REF_CORRECTED` (none hides a shifted REF), 154
`REF_MISMATCH`, of which 109 come from the 134 legacy ANNOVAR-annotated rows;
87 with a REF of 3+ bases match exactly 1–3 bp away (47 at Start−1 with the VCF
anchor still in the alleles); ±10 adds 9 at a higher chance rate; 2-base REFs
match nearby by chance about one time in three.
**Survey.** `bcftools norm --check-ref` exits, warns, excludes or sets/fixes REF
(swapping REF/ALT); it never relocates. maf2vcf skips the row. None explains it.

**Decision (2026-10-05, operator): B.** The row stays FAIL and uncounted (count
the given allele); `gbcms_diagnostic` gives `REF_AT_OFFSET(k)`, every exact
offset within ±3 bp for a REF of 3+ bases, nearest first, from the row's own
position and alleles. An opt-in rescue that counts at the offset would first
need a per-read check on the BAMs.

### I3 — VCF→MAF `Tumor_Seq_Allele1` (#125) · L [decided]
**Finding.** For VCF input, `Tumor_Seq_Allele1`, `Strand` and
`Variant_Classification` are empty. vcf2maf fills `Tumor_Seq_Allele1` from the
genotype.
**Decision.** Fill `Tumor_Seq_Allele1` with REF (the heterozygous
convention), or keep it empty (gbcms genotypes no sample GT). The other two
stay empty: gbcms does not annotate.

**Decision (2026-09-25, operator).** Fill `Tumor_Seq_Allele1` with REF. That is MSK's own convention (100% of 1.13M sign-out rows) and maf2vcf's round trip. No hom-alt inference from VAF. Evidence on #125.

### I4 — maf2vcf's second ALT from `Tumor_Seq_Allele1` (#126) · L [decided]
**Finding.** maf2vcf writes a `Tumor_Seq_Allele1` that differs from both REF
and Allele2 as a second ALT. gbcms reads one variant allele per row. This is
documented.
**Decision.** Keep (recommended: one row, one allele) or genotype the second
allele as its own row.

**Decision (2026-09-25, operator).** Keep one row, one allele (vcf2maf's reading), and document that cBioPortal picks Allele1 for such rows. There are 0 such rows in MSK data. Evidence on #126.

### I5 — Nextflow `convert` module (#127) · L
**Direction.** A `GBCMS_CONVERT` module, only if pipelines need conversion
without counting; the `dna` / `rna` output already converts.

## Merge and outputs

### M1 — Merge rows whose flavors report different alleles (#128) · M
**Finding.** Two cases, both with warnings, can sum different alleles in the
`simplex_duplex_*` columns:
- decomposition: duplex's corrected allele won, simplex's original did (#112
  item 2), and there is no check;
- MNP rescue: the flavors disagree (warned per row since 6.5.0).

**Direction.** Leave the combined columns empty (NA) for such rows, and write a
conflicts file beside the merged MAF (the 6.5.0 T7 follow-up). The
decomposition check lands with C3.

### M2 — Merge inputs from different gbcms versions (#129) · M
**Finding.** 6.5.0 changed VCF-input MAF coordinates. Merging a pre-6.5.0 MAF
with a 6.5.0 one joins badly; the only sign is an INFO line about unmatched
rows. The same happens when some inputs carry `vcf_pos` / `vcf_ref` /
`vcf_alt` and others don't.
**Direction.** WARN when the inputs' `#gbcms vX` provenance lines disagree
across the 6.5.0 representation boundary, or when the VCF-record columns are
present in only some inputs.
**Measured (2026-10-05).** The same 102k-record VCF genotyped by 6.4.0 and by
develop: 6.3% of rows change their variant key; merge exits 0 with 12,771
half-empty rows; its "n/m variants have no <type> counts" INFO counts rows with
REF count 0, not rows an input lacks. 6.5.0 vs develop join fully, but every dev
build since 6.5.0 says `#gbcms v6.5.0` while counts differ. The merged MAF has no
provenance line. No surveyed merger compares producer versions (bcftools, Picard
and htsjdk keep the first `##source` silently); structural mismatches (samples,
contigs) are errors.
**Decision (2026-10-05, operator): B.** (1) The merged MAF carries its own
version and command lines and one line per input with that input's version line.
(2) WARN when the inputs' version lines differ. (3) Refuse the merge when an input
has the pre-6.5.0 VCF-input shape (`vcf_pos` without `vcf_ref`/`vcf_alt`) and
another does not. (4) The INFO line counts the rows each input lacks. (5) Develop
builds say their own version: `6.6.0.dev0` plus the commit.

### M3 — Decomposed-allele hardening (#112 items 3–4: #146, #147) · M
**Finding.** The observations export cannot tell that a variant's per-molecule
rows describe the corrected allele. The `decomposed` list is not
length-validated or padded as `sibling_variants` is: a short list panics the
binned engine and silently drops rows in the legacy path.
**Direction.** Validate and pad `decomposed` like `sibling_variants`, and mark
the winning allele in the observations output (no new MAF columns).

## Observability

### O1 — UMI and no-bases warnings repeated by the rescue recount (#130) · L
With `--rescue-mnp`, the component recount calls the counting pass again with
the same `--umi-tag`, so the "tag never seen" WARN can repeat for one BAM.
C14 (#172) added a second per-pass warning (records stored without bases), which
repeats the same way with a smaller number. Suppress both in the recount; the
main pass's warnings stand for the BAM.
**Measured.** Reproduced (synthetic): the no-bases and `--umi-tag` warnings print
twice per BAM under `--rescue-mnp`; the no-bases count is an upper bound
(overlapping bins). Field: warn once per type (GATK OneShotLogger, htslib,
bcftools); end-of-run counts (GATK, LoFreq); none scopes per input file.
**Decision (2026-10-05, operator): B.** The recount runs quiet (intent across the
FFI), so each BAM warns once, from the main pass; records without bases or
qualities are counted once each (exact, not an upper bound).

### O2 — Run-start summary of enabled options (#131) · L
One INFO block at run start naming the enabled options and what they imply
(e.g. `--rescue-mnp` replaces counts on rescued rows). A 6.5.0 T7 follow-up.
**Measured.** The run start logs the command, version, mode, BAM count and four
settings. Ten real BAMs (IMPACT 4, ACCESS 4, FORTE 2; first 200k records):
duplicates flagged 12–25% in IMPACT, none in ACCESS consensus or FORTE (both
legitimate); binned BQ in FORTE (2/9/24/40) and one ACCESS simplex; bases below
Q20 0.8–2.1%, below a threshold between bins up to 12% (IMPACT Q30) or ~100%
(IMPACT Q40); no ALT contigs, hard-clipped primaries or records without bases;
the proper-pair flag absent on every IMPACT read, but the engine never reads it.
Field: GATK and fgbio log every resolved argument; no tool warns on unmarked
duplicates, BQ binning or hard clips.
**Decision (2026-10-05, operator): B.** One INFO block with every resolved option
(grouped) and a one-line implication for options that change counts or columns;
one INFO facts line per BAM from its first records; WARN only when `--min-baseq`
removes more than 10% of sampled bases, the header has ALT contigs, or more than
1% of primaries are hard-clipped. No duplicate or proper-pair warning.

### O3 — Rescue in fillouts without the MNP (#132) · L
In a fillout of other timepoints or normals, where the MNP itself is absent,
rescue can adopt a germline component. This is documented in `--rescue-mnp`'s
help. Consider a diagnostic when the adopted component is present in the
matched normal. An upstream-annotation question, from 6.5.0 T7.

### O4 — Name the allele the reads carry (`OBSERVED_ALLELE`) · H [decided]
**Why.** Under "count the given allele", a mis-described input gets honest,
low counts. The caveat must say what the reads carry, or the low VAF misleads.
- SOX2: a homopolymer twin case.
- C1's rounds: one mis-described complex variant in 35.

**Rule.**
- For each indel, MNP or complex row, rebuild every spanning read across the
  event (C10's shift region for pure indels, the span otherwise) plus 5 flank
  bases. Aligned bases only in the first cut. Bases below min BQ disqualify
  the read.
- Count identical sequences. The most frequent sequence that is neither REF
  nor the given ALT is the observed allele: n reads carry it exactly, m the
  given ALT.
- Emit in `gbcms_diagnostic` when n ≥ 3, n > m, n ≥ 5% of the scanned reads,
  and the allele is not already an input row. The allele is in canonical VCF
  form (1-based POS).
  - `OBSERVED_ALLELE(chrom:pos:REF>ALT:n/0)`: no read carries the given allele
    exactly, so the input is likely mis-described.
  - `COEXISTING_ALLELE(chrom:pos:REF>ALT:n/m)`: the given allele is present
    (m > 0) beside a more frequent one. Operator, 2026-09-25: distinguish the
    two rather than choose one.
- As built (#164), the comparison is placement-independent: it widens over any
  read indel whose shift region touches the core, and compares canonical
  alleles.
- No new columns, and no count changes.

**Tests.**
- SOX2-shaped synthetic: given `CCCCCC>T`, reads carry `CCCCT`. The diagnostic
  names it; the counts stay the given allele's.
- A correct allele gets no flag.
- Scattered sequencing errors (< 3 identical) get no flag.

**Acceptance.** The RC set: flags only where the reads carry a different
allele, each checked read by read. C1's mis-described variant is flagged.

## Hygiene

### H1 — Writers closed when a write fails (#148, under #133) · L
`_write_output` does not close the writer (or its reference handle) when a
write raises. Use context managers.
**Measured.** A write failure exits 1 but leaves the output at its final path
(truncated) and an unclosed handle. Nine output kinds; the GTF cache already
writes temp-then-rename. Field: htslib tools check close and exit non-zero but
leave the partial file; GATK/Picard write in place; only workflow managers clean
up (Snakemake, Nextflow).
**Decision (2026-10-05, operator): B.** Every output is written to a temp file
beside it, flushed, fsynced, closed with errors checked, and renamed over the
final name; on failure the temp is removed and the final path untouched.

### H2 — `is_indel` in preparation (#149, under #133) · L
`is_indel` reduces to `ref_len != alt_len`; its second clause is exactly
`is_mnp`. Simplify it.

### H3 (#204) — code-quality sweep of the cycle's code · [after T1]
Operator request (2026-10-01). Audited the same day by three read-only reviews;
the ranked work list is on #204. In brief:
- **A. Silent failures that change counts** (test first, measured):
  - S1: an ALT call stands on placement when the reference near a pure indel is
    unavailable (contig ends, unprepared variants);
  - S2: anchor-changing insertions fall back silently from the exact-carrier
    rule in long runs, because prep never widens their reference;
  - S4, QUAL `*` read as Q255, joins C19 #182.
- **B. Degraded modes with no warning:** unprepared variants (S3).
- **C. Monitoring:** carry the deciding rule on `ClassifyResult`; tally per
  variant on the Phase stats line; warn once per variant on unjudged ALT and
  carrier fallback; INFO totals per pass; prep summary counts. No new columns.
- **D. Logging:** levels, misleading texts, 1-based loci and qnames, named traces
  for every decision path.
- **E. Duplication:** CIGAR walks, end-of-walk resolution, the read loops'
  repeated rules, helpers.
- **F. Unused or dead code:** a dead parameter; work computed and then discarded
  (the decomposed twin, `observed_allele` in observations mode); dead branches.
- **G. Comments:** stale ones, and ticket labels.
- **H. Need a decision:** the allele-kind enum vs RNA triage; one fragment class
  shared by observations and mFSD; uppercasing `ref_context` (changes SW on
  soft-masked FASTAs); `mfsd_ref_llr` ulp nondeterminism; for a read that starts
  on a pure indel's left flank base, the windows start at the flank while the
  read-by-bases ALT rule needs an aligned base outside the window (found writing
  the census, T2).

Pure refactors must leave the acceptance output byte-identical.

**Decisions (operator, 2026-10-01).**
- Allele kinds: one enum. A>CCC triages like a delins, `[pos, pos+1)`. A read
  spliced over the anchor that resumes at pos+1 leaves DP (DP-only, RNA-only,
  measured).
- Fragment class: share the classifier and keep both behaviours (no output
  change).
- `ref_context`: uppercase it at prep (fixes SW on soft-masked FASTAs).
- mFSD: sum in a deterministic order.
- Reads starting on a flank: measured (one row, a 66bp duplication regaining 11
  carriers that start on the anchor and hold the whole insert), then decided: the
  read-by-bases rule reads from the flank.
- S2's widening stays; C25 #199 moves into 6.6.0 (its junction windows can make a
  false ALT, now reached by anchor-changing insertions in long runs too).
- A deciding base past a contig edge: depth only (circular contigs), counted and
  warned like other unjudged reads.
- `observe_molecules` requires a reference FASTA (argument or config); without one
  it is a `ValueError`, as for the CLI. This replaces the S3 warning.

**Plan: two PRs.**
- A: the count-affecting fixes (S1, S2, S3, the decisions above), test-first and
  measured.
- B: refactors, dead code, logging, monitoring and comments, byte-identical to A.

**PR A as built (2026-10-01).**
- S1: the ALT read-by-bases rule reads clamped reference and reads from a flank
  the read starts on. Unjudged reads (no reference, or a deciding base past a
  contig edge) are depth only, counted and warned once per variant.
- S2: `carrier::judges` is shared by prep and the dispatcher, so anchor-changing
  insertions and the decomposed twin get a widened reference, and a fallback
  warns.
- `observe_molecules` requires a reference.
- A>CCC uses the delins triage span.
- Reference windows are upper case.
- mFSD sums run in a fixed order.

Reviews:
- The first review found a twin with no reference, warning texts, and C25's ALT
  side (filed and xfailed). All fixed.
- The flank review found that a carrier starting (or ending) on a masked flank
  was ALT, in both the fallbacks and the ALT windows: the flank base must now be
  read (unmasked, the reference's). Also fixed:
  - the census anchors on a flank the same way;
  - C28's example spelled the ALT one base along, so it now holds neither
    allele, with xfails in both directions;
  - a reference that is not a file is refused, in the CLI config and in
    `observe_molecules`.
- Real data after the flank review: develop vs PR A, 140 of 144 files
  byte-identical. The 66bp duplication still gains 11 ALT reads. Three
  deletions lose 5 ALT reads (fragments unchanged); each read's flank base was at
  BQ 9–15, and its bases fit both alleles.
- (Before the flank review, the first measurement: 143 of 144 files
  byte-identical; only the 66bp duplication changed.)

**PR B as built (2026-10-01).** Byte-identical to PR A on prepared input. Measured: PR A vs PR B, 144 of 144 files byte-identical (RC DNA, FORTE RNA, WES), every MAF cell compared. An adversarial review plus a develop-vs-branch fuzz (about 2.2M values) found no count change beyond the declared ones.
- Logging:
  - per-read WARN/DEBUG lines move to trace;
  - rows judged degraded are warned once per pass;
  - truthful texts;
  - read names on per-read traces;
  - named traces for every decision path.
- Monitoring:
  - the deciding rule on each read (`rule=` on `read call`);
  - per-variant tallies on the Phase stats line;
  - one INFO totals line per pass;
  - prep counts its short fetches and warns once.
- Dead code:
  - prep builds the twin only with `rescue_homopolymer`;
  - per-transcript `clip_admission`, `anchor_kept`, the unreachable
    empty-allele guards (the guard moves to the dispatcher), never-written
    `BaseCounts` fields and unused setters are removed.
- Duplication: one helper each for
  - the reference end and soft clips;
  - the scan pad and read window;
  - the read loops' qualities, scoring and molecule key;
  - the insertion and deletion end-of-walk resolution;
  - the molecule classifier (decision 2);
  - the allele kind (decision 1), wherever the existing rules already agreed;
  - plus named constants for the band.
- Not consolidated, because merging would change behaviour:
  - the indel-in-window predicates;
  - the haplotype builders' reference sources;
  - `observed::canonical` vs left-alignment.
- Left for later: the FASTA growth loops (their stop rules differ), and the
  per-variant carrier windows (a performance change).
- Comments: stale ones fixed; ticket labels removed from code and log text.
- Follow-ups filed:
  - C29 #207 (6.6.0): `check_complex`'s inline query walk counts hard clips;
  - C30 #208 (6.7.0): two pure-indel tests disagree with the allele kind on
    lowercase or unprepared input.

## Performance (M5 leftovers)

### P1 — Deep-bin fetch reduction (M5b) (#150, under #134) · L
Deep cfDNA bins read 150k+ reads to count a few variants. Narrowing the fetch
is the only remaining cfDNA lever. It must keep binning invariance, so scope it
behind the binning-invariance tests (T1). Investigation first.

### P2 — Bin cost-sort (PF-2) (#151, under #134) · L
Niche: cfDNA has no long-pole bin. Cheap if a skewed workload appears.

### P3 — Document the bin-span soft floor (LO-3) (#152, under #134) · L
Doc only. Never cap the span: a cap risks re-breaking the bin-anchor
invariant (`.agents/memory/bin-anchor-coverage.md`).

## Statistics (accepted deviations)

### S1 — Mean LLR per fragment (CR-5) (#153, under #135) · L [decided]
Report the LLR per fragment rather than the sum. It changes displayed values,
so coordinate with report consumers. Only if a consumer needs cross-variant
comparability.

**Decision (2026-09-25, operator).** Report the per-fragment mean LLR (with n alongside) in the existing LLR columns; display only. Evidence on #153: on real cfDNA the sum grows about 22× with n, while the mean is flat.

### S2 — `MIN_FOR_KS` floor (ME-9) (#154, under #135) · L [decided]
Raise the floor above 5 only if a power analysis justifies it. With exact
small-N KS and the `ks_valid` gate, 5 is defensible.

**Decision (2026-09-25, operator).** Keep `MIN_FOR_KS` at 5: the exact test is calibrated at every n measured (false-positive rate about 0.05). Change the report so **CH-LIKE no longer rests on a non-significant KS at low power**. It requires an ALT-fragment count at which the test has useful power, derived from the real-data power curve (about 10% at 5–10 fragments). Evidence on #154.

### Group 6 — measured (2026-10-05; decisions pending)
Operator direction (2026-10-05): mFSD is **graded evidence** (increased confidence
toward CH or toward somatic, or none), never a hard origin call; the fragment test
is **plasma-only** (ALT vs REF at the same locus, same sample); the buffy coat is
the answer key for validation, never an input; the regime that matters is **2–5
duplex ALT fragments** (where CH lives); production runs mFSD on the duplex BAM.
Harness `~/test/gbcms/harness/g6/` (local; PHI stays there).

Data and labels:
- ACCESS plasma calls from the signed-out cBioPortal dump, labeled by the
  patient's buffy coat: CH when the buffy coat has ≥ 3 ALT reads (VAF ≥ 0.2× the
  plasma VAF); not CH when it has 0 where it could have seen the variant (≥ 5
  expected ALT reads at the plasma VAF; ACCESS buffy coats run at about half the
  plasma depth, median 1,013 vs 2,106 fragments). Germline excluded everywhere
  (`Mutation_Status` not SOMATIC, buffy VAF ≥ 0.25, gnomAD AF ≥ 1%, or tagged
  GERMLINE in the patient's IMPACT samples). The dump's plasma counts are duplex +
  simplex fragments, its normal counts the buffy coat's `unfilter` BAM (both
  matched our own genotyping on 211 and 195 of 214 pilot rows). The dump's
  `COMMENTS` and `SO_COMMENTS` are empty, so no sign-out CH annotation is used.
- Pilot: 39 patients, a CH call (DNMT3A/TET2/ATM/PPM1D) plus tumor calls in the
  same plasma sample; 136 concordant labels (97 tumor, 39 CH, 3 of them TI-CH).

Measured:
- KS p-value: conservative, never over-calling (null false-positive rate
  3.2–4.6% at 5–20 ALT fragments vs a permutation test's 4.1–5.2%); most real
  rows use the asymptotic branch (348 of 512), which overstates p 1.7–5× at small
  n (survey simulation). Exact at any depth is cheap without overflow (shares of
  in-band lattice paths, computed by anti-diagonal).
- Today's classes (pilot, pooled flavors): CH rows CH-LIKE 29/39; tumor rows
  TUMOR-LIKE 22/97, CH-LIKE 3/97 — 3 of 17 tumor variants in CH-associated genes
  (12 TP53) are CH-LIKE, and a tumor TP53 variant can never be TUMOR-LIKE (gene
  membership gates the size evidence). At 2–5 fragments the report says
  INSUFFICIENT.
- Significant shortening (KS or short-fragment excess, BH q < 0.05) never fired on
  a CH row (0/39) and found 44% of tumor rows (16% at 5–19 ALT fragments, 75% at
  ≥ 50). "REF-like sizes" fired on about 30% of tumor rows below 50 fragments: no
  basis for a CH lean from a non-significant test.
- Where ALT fragments fall within their locus's REF density (mid-rank percentile;
  1.0 = same as REF), duplex: tumor 2.0× in the shortest decile, 0.6–0.8× through
  the mono-nucleosome core, about 1.1× in the longest; CH flat (about 0.8–1.3).
  Tumor fragments avoid the REF core in both directions.
- Exact percentile tests (Irwin–Hall null, valid from one fragment, no model):
  "tails" (sum |2u − 1|) flags tumor 12/15/16/20% at 2/3/4/5 duplex fragments vs
  CH 4/4/4/3% (chance); "shorter" (sum u) 13–16% vs 4%. A likelihood ratio from a
  tumor percentile curve learned leaving the patient out (an e-value under "ALT
  like REF") leans at 3:1 in 5–20% of draws at 2–5 fragments, right about 3–7×
  as often as wrong; AUC 0.65–0.72 vs 0.60–0.63 for today's fixed-Gaussian LLR.
  Size carries about 0.07–0.10 nats per fragment, so 2–5 fragments give weak
  evidence (expected likelihood ratio about 1.5 at 5): "no evidence" is the usual
  honest answer.
- Within the CH-associated genes (pilot, duplex): buffy evidence strong → REF-like
  fragments (tails 3%); weak (1–2 buffy reads) 41%; none 48% — a gene confound
  remains (strong mostly DNMT3A/TET2/PPM1D/ATM, none mostly TP53).
- Pooling duplex + simplex fragments detects more (44% vs 25–29% per flavor), but
  per-fragment information is the same in duplex alone.
- C31 #212 (fragment end from the mate's unclipped 5' end): 0.01–0.39% of
  fragments have a clipped 5' end and ≤ 0.026% cross 150 bp; 3 of 7 local ACCESS
  BAMs carry no MC tags. Negligible for mFSD genome-wide; indel loci not yet checked.
- Empty fragment classes write mean size 0 and LLR 0 where the docs say NA (the N
  class is empty on every row measured).
- Survey (2026-10-05): no published method makes a validated per-variant call from
  1–5 fragments on fragment features alone (plasmaCHORD reaches AUC 0.84 at 3–5
  reads mostly from gene, VAF and age); end motifs and nucleosome positions show no
  documented per-fragment gain; a patent (Illumina US12054712B2) claims a
  fragment-size likelihood-ratio test for CH vs cancer — for the operator's
  tech-transfer office.

Held-out validation pairs (2026-10-05; feed the separate CH model, not the gbcms
decisions): calls in the CH-associated genes at native 2–5 duplex ALT fragments, 199
patients not in the pilot, 128 CH vs 87 not CH, the tumor curve frozen from the pilot.
CH fragments sit like REF (deciles 0.68–1.28; no CH shortening); not-CH fragments 2.75x
in the shortest decile. Tails test p < 0.05: CH 5%, not CH 29%; curve LR leans somatic
5% vs 33%, leans CH 16% vs 3%; AUC 0.76. Within TP53 (36 vs 48): leans somatic 3% vs
31%, leans CH 17% vs 0%. The 6.5 report called 18 of 87 not-CH calls CH-LIKE.

**Decided (2026-10-05, operator): keep mFSD in gbcms lean; the CH-vs-tumor
prediction is a separate model.** gbcms group 6 builds:
1. S1 as decided: `mfsd_alt_llr` / `mfsd_ref_llr` are the mean per fragment (n is
   `mfsd_*_count`); an empty class writes NA for its mean size and LLR (the docs
   already say so).
2. S2: the report no longer calls CH-LIKE from a non-significant test, and gene
   membership no longer gates the somatic direction (a tumor TP53 variant can lean
   somatic). Classes: **leans somatic** (significant ALT shortening) / **no size
   evidence** / **insufficient**; the CH-gene flag shown as a note beside them.
3. The KS p-value exact at any depth (shares of in-band lattice paths, by
   anti-diagonal; no overflow), replacing the asymptotic branch above n·m 10,000.
4. `mfsd_alt_confidence` tiers renamed so a 5-fragment class is not "HIGH".
5. `<sample>.fsd.parquet` (each variant's REF and ALT sizes) stays the interface for
   models; no new evidence columns in gbcms.

The CH model lives in its own private repository (operator, 2026-10-05): inputs the
gbcms MAF + `fsd.parquet`; features such as each fragment's REF percentile and size,
the exact tails/shorter tests, a tumor-curve likelihood ratio (cohort or per-sample
adapted; +0.02–0.03 AUC on the pilot), VAF, gene and fragment counts; trained on
buffy-labeled ACCESS (about 760 strong vs 4,400 none in the CH genes), cross-fitted by
patient; evaluated forensic-style (Tippett plots, Cllr, calibration); output a graded
likelihood ratio or posterior, never a hard call. The validation pairs are its first
held-out set. The repository is local only for now (no remote; operator, 2026-10-05).
CMO-CH (white-blood-cell genomic DNA, UMI, ACCESS methodology; 596 targets, 0.14 Mb,
full TP53/TET2/ASXL1/DNMT3A/PPM1D/CHEK2/ATM/SF3B1/SRSF2/U2AF1/U2AF2; standard
~1,000x or high sensitivity >10,000x raw; calls to 0.1% VAF; many false positives from
the small, deep panel) is sheared gDNA, so its fragment sizes are not cfDNA and cannot
train the size model. It is the deep blood truth for patients who also have ACCESS
plasma: it resolves the calls the ACCESS buffy coat (~1,000x) cannot (1–2 reads, or too
shallow). Its own false-positive problem (real CH vs artifact, read-level features) is a
possible later task for the same repository.

Later (6.7.0, idea, not filed — operator 2026-10-05): an optional matched-normal BAM
for mFSD, read for depth-aware counts only (0 ALT counts against CH only when the
normal could have seen it), shown beside the plasma fragment evidence, not merged
into it; never the normal's fragment sizes (sheared genomic DNA, not cfDNA).
Pooling duplex + simplex fragments at merge time is a second candidate.

## Release and docs infrastructure

### D1 — CI version-consistency check (#136) · M
The 6.4.0 cut missed the Nextflow banner (`nextflow/main.nf`, still `v6.3.1`
until 6.5.0). A CI check that the release guide's 11 version references agree
would catch it.

### D2 — Release workflow creates the GitHub Release (#137) · M
6.4.0 has a tag and published artifacts but no GitHub Release page. Create the
page from the CHANGELOG section in the release workflow.

### Group 7 — measured and decided (2026-10-05, operator)
- **D1 #136: single source + check.** Measured: 11 hand-edited references in two
  formats (`6.6.0.dev0` / `6.6.0-dev.0`); Nextflow trails on purpose at the last
  released image (images publish only on tags); nothing checks the tag against the
  version before PyPI/Docker publish (a `6.6.0` tag on `6.6.0.dev0` would publish a
  dev wheel; `skip-existing` hides a re-tag). Decided: Python takes the Cargo
  version through maturin (`dynamic = ["version"]`, `__version__` from the
  package metadata); the Nextflow modules and banner read
  `workflow.manifest.version`; hand edits drop to `rust/Cargo.toml`, the Nextflow
  manifest and the CHANGELOG (`Cargo.lock` follows `cargo check`). A check script
  runs on every PR (single source intact, lock = manifest) and gates the release
  workflow before anything publishes (tag = Cargo = Nextflow manifest, a dated
  CHANGELOG section); a test keeps the release guide's table equal to the checked
  locations. `gbcms --version` and the Nextflow lint must stay identical.
- **D2 #137: automatic Release page.** Measured: the 6.4.0 / 6.5.0 pages were made by
  hand (body = the CHANGELOG section, title `X.Y.Z — summary`, no assets). Decided:
  the release workflow creates the page after PyPI and Docker succeed, body from the
  CHANGELOG section, title from its heading (`## [X.Y.Z] - YYYY-MM-DD — summary`),
  wheels + sdist + SHA256 sums attached, with GitHub build-provenance attestations.
- **D4 #139: floors + weekly latest + Docker lock.** Re-surveyed 2026-10-05: Python
  majors rich 15, pytest 9, mypy 2.4, black 26.10 (CI pins 26.5); Rust majors pyo3
  0.29, rust-htslib 1.0.1, bio 4.2, arrow/parquet 60, noodles-gtf 0.58, statrs 0.19,
  bincode 3 (GTF cache format + `CACHE_FORMAT_VERSION`). The dev venv lacks
  `pytest-mock` and `types-pyyaml`. Decided: keep floors; a weekly CI job on the
  latest releases; Docker built from a lock. PR A: Python, CI, Docker lock. PR B: the
  Rust majors one at a time, the RC set byte-identical after each (isolated builds).
- **D4 GTF loading (noodles-gtf 0.58, bincode 3): no cache, one fast parser.**
  bincode was a default in M5a, never compared. Surveyed tools: those with their own
  databases ship them prebuilt (SnpEff, the VEP cache, CTAT); most parse the GTF
  each run (bcftools csq, featureCounts, Arriba, regtools); bgzip+tabix is the
  standard index (VEP GTF mode, igv.js, JBrowse), but Ensembl, GENCODE and iGenomes
  ship none. Measured on GRCh38.111 (`harness/gtfidx/`):
  - **noodles was the cost.** It parses every line before the feature check:
    5.9 s (6.3 s with the index build; 9.4 s in `build-gtf-cache`). A byte-level loop
    checks the feature column first: 0.7 s (1.1 s with the build) for the same
    16 chromosomes, 1.3 s for the whole genome, 1.7 s from `.gtf.gz`. Peak memory
    is unchanged (357 MB for 16 chromosomes).
  - **Identical exons.** Checksums match on both chromosome sets. A per-line diff
    against noodles on Ensembl 111 and GENCODE v50 comprehensive, v50 basic and
    v47lift37 (26.1M lines, 12.0M exons) finds no line that differs.
  - **Edge corpus (34 lines).** noodles reads one space between key and value
    (`transcript_id  "T1"` gives the ID ` "T1"`, quotes included). It keeps
    start > end, keeps an empty `transcript_id` (all such exons become one
    transcript), and keeps coordinates past i32, which wrap negative.
  - **The cache saves 1–1.4 s a sample on Ensembl.** A hit costs 0.9 s in the CLI
    (~0.3 s of it is startup), and the cache is 97 MB per chromosome set. It brings a
    build step, a Nextflow process and the bincode format, and bincode 3.0 is a
    `compile_error!` tombstone.
  - **Shipped loader (review round), `gbcms rna` load time:** Ensembl 111 1.6 s for
    16 chromosomes (noodles 8.9 s), 2.0 s whole genome (9.3 s), 2.4 s from
    `.gtf.gz`; GENCODE v50 basic 4.0 s and comprehensive 6.5 s (`.gtf.gz`, whole
    genome, 1.3 GB peak). Slower than the prototype above, which skipped UTF-8 and
    full attribute checks; zlib-rs replaces miniz_oxide for gzip.
  - **Acceptance (develop vs the final build, plain GTF and a BGZF copy in the same
    line order): 0 changed rows in any column.** FORTE over the mount, one BAM at a
    time: the RC truth cohort (33 samples, 94 rows) and the T9 probes on the 3 T6
    samples (978 rows). A full local FORTE BAM: the C1 splice probes (2,408), T9
    (326), T10 (40) and the C16 masked-window probes (56); every row has a boundary
    distance and 2,778 of 2,830 have transcripts. The review rebuilt the index both
    ways and found it byte-identical (Ensembl 111, GENCODE v50 basic).
  - **Why the cache was added (M5a, #55) and what is left of it.** RNA runs were
    ~90% GTF parse (~8 s vs ~40 ms of counting), repeated per Nextflow task, and a
    pre-build step was needed because concurrent tasks all missed the cache. Now:
    a warm-cache run took 0.86 s / 1.02 s (16 / all chromosomes) against 2.05 s /
    2.5 s plain and 2.4 s / 2.9 s from `.gtf.gz`. The cache would save 1.2–1.9 s a
    task, while the pre-build was a serial job every task waited on. Per-task reads
    were 96–125 MB with the cache; the `.gtf.gz` is 61 MB (the plain GTF 1.46 GB),
    so the Nextflow page recommends it. Memory is the same either way.
  - **`.gtf.gz` is advertised, but it crashes.** The config validator accepts
    `.gtf.gz`, and the parser fails with "stream did not contain valid UTF-8".
  - **Tabix gives exactly the same answers.** All seven queries
    (strand, nearest splice distance, intron boundaries ×2, overlapping transcripts,
    transcript introns, junction known) agree on 5,757 positions (the probes plus
    boundary, random and gene-desert positions) and 111,911 junction probes:
    0 mismatches. The fetches are a variant window, each found transcript's span
    from its transcript record, a nearest-boundary window that grows until the best
    boundary lies inside it, and exons ending within tolerance of a junction's
    start. Building the index takes 6.6 s (61 MB + 318 KB). The Python cost is
    1.9 ms a variant; the Rust cost is not measured.
  Decided:
  - Replace noodles-gtf with a byte-level parser that reads `.gtf` and
    `.gtf.gz`. It keeps noodles' grammar field by field, except that it allows
    whitespace runs between key and value, trims unquoted values, and rejects
    start > end, coordinates past i32 and an empty `transcript_id` (each counted).
  - Deprecate the cache in 6.6.0: `--gtf-cache-dir`, `build-gtf-cache` and the
    Nextflow `--gtf_cache` are accepted and ignored with a warning. The Nextflow
    step, the cache module and bincode go now; the flags go in 6.7.0.
  - Tabix stays out. The design above is the answer if memory ever becomes the
    limit.
  - A PR of its own, before PR B, which then has no noodles or bincode migration.
  - A public GTF API (a narrow wrapper over `_rs`, like `observe_molecules`) is a
    candidate, not a ticket: no consumer reads GTFs today (mulligan,
    ch-fragmentomics), and the annotation answers ship as output columns. The
    parser (`gtf_line.rs`) and the index's seven queries are self-contained if one
    appears.
- **D4 PR A: Python, CI, Docker lock (2026-10-06, operator).** Measured on macOS
  arm64 against `develop` (`harness/d4/py/`):
  - **Every Python works.** On 3.10, 3.11, 3.12, 3.13 and 3.14 with the latest
    releases (pysam 0.24.1, typer 0.27.2, click 8.5.0, rich 15.0.0, pydantic 2.13.5,
    polars 1.44.2), the suite passes (1,294 each). CI tests 3.11 and 3.12; the
    declared minimum (3.10) and the newest (3.13, 3.14) were untested.
  - **The declared floors were false.** At `lowest-direct` the CLI crashes: typer
    0.9.0 cannot read its `list[...]` options. pysam 0.21.0 has no macOS arm64 wheel,
    and its sdist fails on current setuptools (`pkg_resources`); its Linux wheels
    pass the suite (review). click is no runtime dependency at all: typer 0.2x
    bundles its own, and only a test imports it.
  - **Measured minimums.** Bisected one package at a time, the others at latest:
    typer 0.15.4 (older typer breaks with click >= 8.2, which `click>=8.0` allows)
    and pysam 0.22.0. click 8.0, rich 13.0, pydantic 2.0 and polars 1.0 hold. All
    minimums together pass on 3.11, and on 3.10 once a test helper stops calling
    `issubclass()` on a generic alias.
  - **Two dev lists, drifted apart.** `maturin develop` needs pip >= 25.1 because
    `[dependency-groups] dev` makes it run `pip install --group`. That table and the
    `dev` extra (what CI installs) diverged: the group lacked pytest-mock,
    types-pyyaml, pyyaml and mkdocs; the extra lacked pyarrow, which the tests
    import. The developer guide's `python -m venv` + `maturin develop` fails on
    Python's bundled pip.
  - **The image has no lock.** Its dependencies resolve fresh at each build, and CI
    never runs anything inside the image. A hash-pinned lock for linux/amd64 and
    Python 3.11 is 16 packages. One Ubuntu step in `test.yml` can never run.
  Decided:
  - Floors become `typer>=0.15.4` and `pysam>=0.22.0`, and click moves to the
    `test` group. A PR CI leg (Ubuntu, Python 3.10) installs the floors first
    (`uv --resolution lowest-direct`) and `scripts/check_floors.py` checks them.
  - On every PR, CI tests Ubuntu 3.10 (at the floors), 3.11 (the image's) and 3.14,
    plus macOS 3.12. The classifiers list 3.10–3.14.
  - One dev list: PEP 735 groups `test`, `docs` (also what `deploy-docs.yml`
    installs) and `dev` (which includes both). CI installs them with uv, and the docs
    say pip >= 25.1.
  - The image installs a hash-pinned lock (`--require-hashes --no-deps`, then
    `pip check`), refreshed as a release-guide step, and CI runs the image.
  - A monthly job (and a manual run before each release) runs the latest Python
    releases on 3.10–3.14, with semver-compatible Rust updates. On failure it opens or
    updates one issue. It is monthly, not weekly or quarterly: PR CI already installs
    the latest releases, and GitHub disables schedules after 60 days without activity.
- **D4 PR B: the Rust majors (2026-10-06).** One dependency per commit. Each
  step passes clippy, `cargo test` and the Python suite. The local FORTE acceptance
  (`harness/d4/rust/`) compares every column and every Parquet table. It runs RNA
  mode with the GTF (the STAR and UMI-dedup BAMs × the C1, T9, T10, C16 and 806
  read-observed indel probes) and DNA mode (mFSD and observations Parquet × both
  backends): 9,536 rows. It was run after step 1, after steps 1–5, and on the final
  tree. The cluster acceptance (28 RC DNA, 33 RNA truth, 3 probe samples, 80 WES
  loci; 2,212 rows) was run on the final tree.

  | Dependency | From → to | Outcome |
  |:--|:--|:--|
  | semver-compatible (127 lock entries) | — | upgraded |
  | statrs | 0.18 → 0.19.1 | upgraded; Fisher's exact test bit-identical under both versions on every table with cells 0–40 |
  | rust-htslib | 0.51 → 1.0.1 | upgraded; no code change, same hts-sys 2.2.1 |
  | arrow / parquet | 53 → 60 | upgraded; the Parquet companions are identical by content (string column-chunk statistics truncate at 64 bytes now) |
  | pyo3 (+ pyo3-log) | 0.27 → 0.29.3 (0.13.4) | upgraded: `allow_threads` → `detach`; `Variant`/`BaseCounts` opt in to `from_py_object`, the output-only types opt out; `gil_used = true` keeps the module's GIL behaviour |
  | wfa2lib-rs | rev 7d6ec921 → cb8aa3d6 | upgraded, `default-features = false`: clap, tracing and mimalloc (the benchmark binary's) leave the build |
  | bio | 3.0 (4.2.1 available) | **held** (operator, 2026-10-06); see below |
  | noodles-gtf, bincode, serde | — | removed in #230 |
  | flate2 | 1.1.10 | zlib-rs backend since #230 |
  | coitrees 0.4, bio-types 1.0.4 | — | already current |

  **bio 4, measured and held.**
  - **Smith-Waterman.** The gap rule changed (rust-bio#660): a k-base gap scores
    `open + extend·(k−1)`. With `SW_GAP_OPEN` −6, every score is identical on 199,949
    random pairs, for semiglobal, local, global and clipped alignments.
  - **PairHMM.** The forward-algorithm fix (rust-bio#701) costs each extended
    inserted base its emission. Our gap probabilities are symmetric, and
    haplotype-side gaps emit with probability 1.
    - For uniform-quality reads, the LLR moves at most 0.004 at Q37 and 0.27 at
      Q20, with no call crossing ±2.3, on insertion carriers and on REF reads at
      deletions alike.
    - Reads whose inserted bases are N or Q2 move across the threshold (synthetic:
      1.8 → 5.0). The review's side-by-side model crossed on 1,349 of 22,320 such
      cases.
  - **Real data.** No row changes: 9,536 local rows and 2,212 cluster rows with the
    bio 4 tree.
  - **Why held.** The engine already credits ALT to reads whose inserted bases are
    all N or below min BQ at non-repeat insertions, under both versions. That breaks
    invariant 7 and invariant 2. So bio 4 lands with a gate that requires the
    inserted bases to carry the ALT, re-measured.
- **Fisher's exact test overflows (found in the D4 PR B review, operator: fix in
  6.6.0, own PR).** `fisher_exact_2x2` returns p = 0 for any table above ~1,030
  reads: statrs's binomial overflows, the pmf turns NaN, and the sum stays 0.
  Examples: (258,258,258,258) gives 0 against an exact p of 1.0; (30,25,500,520)
  gives 0 against 0.49.
  - **Affected:** `strand_bias_p_value` and `fragment_strand_bias_p_value`
    (`SB_PVAL`/`FSB_PVAL`) at ACCESS depth, and the ASJD junction test (p < 0.05).
  - **Fix:** log-space hypergeometric, validated against an exact reference on a
    grid that includes deep tables; measure the real rows whose p-values and ASJD
    flags change.
- **N or low-BQ inserted bases credited ALT (found in the review, operator: measure
  in 6.6.0, then decide).** Under the PairHMM backend, reads whose inserted bases
  are all N or all Q2 are called ALT at non-repeat insertions, under both bio
  versions. Trace the RC DNA set for ALT calls that rest only on such bases. Fix in
  6.6.0 if they occur, otherwise in 6.7.0 with bio 4.
- **Platforms (operator, 2026-10-06):** Red Hat first, macOS second, filed for 6.7.0
  as D7 #232. The sub-issues: the image under Apptainer in CI (#233), a wheel that
  installs on RHEL 8 (#234), wheels for Python 3.10–3.14 (#235), a native arm64
  image (#236) and macOS wheels (#237). In 6.6.0, Docker/Singularity remains the
  RHEL path.
- **S3 #239: Fisher's exact test at depth (2026-10-06, operator).** Measured on the
  cached D4 acceptance outputs (`develop`), each row's strand-bias table recomputed
  exactly with scipy:
  - **Overflow.** `fisher_exact_2x2` returns p = 0 once C(n, a + c) passes f64's
    range (statrs's binomial overflows): both strands deep, from ~1,030 reads on a
    strand-balanced table; a table with one thin strand stays correct. Rows printing p = 0 where the exact p
    is above 0, out of all rows:
    - RC cfDNA: read test 316 of 1,060 (298 of them have an exact p ≥ 0.05),
      fragment test 235;
    - FORTE RNA mode: 233 of 7,272 (read), 134 (fragment);
    - FORTE DNA mode: 164 of 2,264 (read), 70 (fragment);
    - WES 3 of 80, RNA truth 1 of 94, RNA probes 11 of 978.
  - **ASJD.** 5 of the 7 local junction flags rest on p = 0.
  - **The tie tolerance was absolute (1e-10).** Strongly biased tables are floored
    near 1e-10 (e.g. (250, 50, 20, 280) prints 2.1e-10; the exact p is 7.0e-90). On
    real rows at n ≤ 1,030 only 3 differ beyond the printed precision, all at
    p < 1e-9.
  - **Survey.** R's `fisher.test` uses a relative tolerance of 1 + 1e-7 and scipy
    1 + 1e-14; both work in log space. GATK's FisherStrand scales tables above 400
    reads down to 200 before testing.
  - **Merge** recomputes combined strand bias through the same binding.

  Decided:
  - The exact two-sided p on the raw counts (no GATK-style normalization; the odds
    ratio column gives the effect size), with R's tie rule: the sum over tables with
    probability ≤ the observed × (1 + 1e-7). Computed in log space.
  - Validated against an exact integer oracle; then measure the real rows and ASJD
    flags that change.

  Review: no blockers. A recurrence started at the low end of the range drifts as
  its sum grows (~n ln 2 on balanced tables). Against the exact oracle on 3-SD
  tables (p ≈ 0.003): 5.8e-11 at n = 10⁵, 9.3e-10 at 10⁶. It now starts at the mode
  and stops where terms underflow (e^-750): 2.0e-15 and 5.0e-14 on those tables,
  and u32-max cells take 40 ms instead of a 34 GB allocation.
  Tests now pin R's tie rule against scipy's on three near-tie tables.
- **C35 #240: inserted bases nobody can read (2026-10-06, measured).** Develop, the
  reviewer's synthetic case first, then real data (`harness/c35/`):
  - **Both backends, four paths.** A same-length insert whose bases are all N or all
    below min BQ is "unverifiable" on the strict path and goes to Phase 3, which calls
    ALT on length alone (REF pays for the gap) under PairHMM and SW alike. Shifted
    placements reach Phase 3 the same way. `read_bases_fit_alt` (split ops, a deleted
    anchor) needs one readable base anywhere between the flanks, so flank bases
    suffice. `insert_truncation_match` ignores BQ. The post-splice twin already calls
    the same read partial (pinned by the RNA splice contract). At a homopolymer only
    the repeat gap penalties keep the LLR under 2.3 with bio 3; bio 4 crosses it.
  - **RC DNA (28 samples, 82 insertion rows).** Same-length carriers: 11,487 clean
    and 660 partly masked at the junction, 172 readable shifted; 157 with no readable
    inserted base (22 junction, 135 shifted), of which develop calls 127 ALT, all in
    Phase 3. 105 of those are one shape: a 1bp insertion into a homopolymer run of
    3–8 with a duplex-masked N inside the run, the aligner writing the N as the
    insertion. The rest (22) are low-BQ letters. Three rows of one sample lose 62, 15
    and 11 ALT reads (AD 430, 299, 111); the other 23 rows 0–8. 127 fragments have
    no readable insert on either read. Residuals: 0 truncation candidates without a
    readable base (3,412 readable), 10 reads with the insert in a soft clip at the
    junction, all bases unreadable (not traced).
  - **FORTE (local).** DNA mode 6 such reads (3 ALT today), RNA mode 14 (8 ALT), all
    low-BQ letters, no N. BAQ spares the variant's own insertion, so RNA BAQ does not
    add to them except within 5bp of a splice.
  - **Survey.** Likelihood callers (GATK HC/Mutect2, bcftools, Strelka2) credit ALT on
    length: N matches anything, low BQ is floored (Q6/Q7). Sequence-keyed counters
    never credit an N insert: VarDict and freebayes drop it, bam-readcount and LoFreq
    report it as its own allele, and GetBaseCountsMultiSample (`baseCountIndelDMP`,
    exact string match) counts it as neither (DP only); it checks BQ on the anchor
    only, so low-BQ letters that match count ALT there (`--generic_counting` drops
    them). No literature on counting indels from N-masked consensus reads.

  Decided (operator, 2026-10-06), 6.6.0, own PR:
  - An insertion's ALT needs one of the read's own inserted bases read: not N and at
    or above min BQ (one quality gate, invariant 2; its own bases, invariant 7). A
    read whose inserted bases are all unreadable is partial (depth and partial_alt),
    on every path: the strict junction, shifted placements, split ops, a deleted
    anchor, truncations. Partly masked inserts whose readable bases match stay ALT.
  - No new output signal: partial_alt carries them, a trace line names the cause.
  - The read census gets the same rule as policy.
  - Inserts in a soft clip at the junction (10 RC reads): a sub-issue.
  - Then bio 4 is measured on top.

  Review (one round, adversarial): a blocker (a deletion outside the flanks offset
  the deleted anchor in `read_bases_fit`), an unreadable insert flagged anywhere in
  the scan window (outside the discrimination window it is a separate event, RJ-8),
  and the census netting deletions only between its anchors. All fixed red-first.
  Acceptance on the final tree (develop's outputs as the base): RC DNA 20 of 1,060
  rows, ALT −127 (exactly the measured reads), partial +126, REF +1, fragment ALT
  −103, depth unchanged; WES 8 of 80 rows, ALT −21; FORTE DNA −8 ALT (both backends),
  RNA −9; RNA truth and probes unchanged.
  Residuals measured after the review:
  - The soft-clip reads: develop calls all 10 RC reads (and 12 FORTE) neither, never
    ALT. Nothing to fix; recommended: no sub-issue (operator to confirm).
  - The truncation band reads letters at any quality. Both policies on the RC
    truncation candidates: 3,387 pass both, 0 pass only by low-BQ letters, 1 passes
    only with masked identity (an N counted as a mismatch), 24 fail both; FORTE 17
    fail both. Recommended (operator to confirm): the band stays as it is (the
    identity-band rule), with a readable base required.
- **P3 #152: doc the floor.** Measured by replaying the bin rule: `BIN_WINDOW` (10 kb)
  is a floor, not a maximum — each member extends the end by its span plus half a
  window, and the 200-variant cap stops dense inputs. Per-sample signed-out lists
  (141,845 samples, 1.04M bins): median 10 kb, p99 11.6 kb, max 29 kb (1.6% past the
  window, 0 cap splits); one input holding the ACCESS cohort union (14,296 variants):
  p99 36 kb, max 47 kb (4.7x; 19 cap splits); the IMPACT union (424,572): p99 23 kb,
  max 44 kb (1,458 cap splits). The architecture page also still cites the retired
  parity suite.
- **D5 #155:** built last, against the release candidate (after D4).
- **After the 6.6.0 tag (operator, 2026-10-05): the downstream org repos.**
  - `mskcc-omics-workflows/containers` gets `containers/gbcms/6.6.0/Dockerfile`,
    which builds `ghcr.io/mskcc-omics-workflows/gbcms:6.6.0`.
  - `mskcc-omics-workflows/modules` bumps the `gbcmsrs` modules from 6.3.1:
    - drop `buildgtfcache` and the `rna` module's `gtf_cache` input (they still
      work under 6.6.0, with a warning; 6.7.0 removes them);
    - widen the `gtf` pattern to `*.{gtf,gtf.gz}`;
    - check the 6.4–6.6 CLI and output changes, then update the nf-test snapshots.
- 6.7.0 overlap: D3 #138 stays (D3a pins below MkDocs 2); P1 #150 stays (P3's spans
  inform it).

### D3 — mkdocs-material 2.0 (#138) · L
The docs build prints mkdocs-material's MkDocs 2.0 incompatibility notice. Pin
`mkdocs-material<2`, or plan the migration.

### D4 — Dependency upgrade audit: does anything break on current releases? (#139) · M
**Why.** The runtime Python dependencies are floors only (`pysam>=0.21`,
`typer>=0.9`, `click>=8.0`, `rich>=13`, `pydantic>=2`, `polars>=1`), and CI and
Docker install without the lockfile. So a fresh install already gets the latest
releases, whether or not we have tested them.

**Survey (2026-09-25, read-only).** Python, dev venv vs latest:
- **Runtime:** pysam 0.23.3 → 0.24.1, typer 0.19.2 → 0.27.2, click 8.3.0 → 8.5.0,
  **rich 14.1 → 15.0 (major)**, pydantic 2.11 → 2.13, polars 1.42 → 1.44.
- **Dev:** black 25.9 → 26.5 (CI already runs 26.5), ruff 0.13 → 0.16,
  **mypy 1.18 → 2.3 (major)**, **pytest 8.4 → 9.1 (major)**, mkdocs-material
  9.7.4 → 9.7.7 (see D3).
- **Build:** maturin is capped `<2.0`. Docker base is `python:3.11-bookworm` /
  `-slim-bookworm`. Nextflow is `>=24.04` (the config notes the strict parser
  ≥26.04).
- **Venv hygiene:** a stale editable `py-gbcms` 2.8.0 (the old package name)
  sits in the dev venv; remove it.

Rust, `Cargo.toml` vs crates.io — mostly semver-major moves:
- **pyo3 0.27 → 0.29** (+ pyo3-log 0.13.2 → 0.13.4).
- **rust-htslib 0.51 → 1.0**, **bio 3 → 4** (bio-types 1.0.4).
- **arrow / parquet 53 → 60**, **noodles-gtf 0.37 → 0.57**, **statrs 0.18 →
  0.19**.
- **bincode 1.3 → 3.0.** The GTF index cache (M5a) is bincode-serialized, so a
  bump changes the on-disk format. `CACHE_FORMAT_VERSION` must be bumped with it
  so old caches are rebuilt, not misread.
- rayon 1.10 → 1.12 and thiserror / flate2 are compatible; coitrees is
  current; wfa2lib-rs is pinned to a git rev.

**Method.**
1. Work in a scratch worktree with an isolated venv (`maturin develop` from the
   worktree root; see the worktree and maturin memory notes).
2. **Python:** install the latest of everything from `pyproject.toml` (no
   lock). Run the full gate (ruff, black, mypy, pytest, mkdocs). Record each
   break against the package that caused it.
3. **Rust:** `cargo update` (semver-compatible) first, then the gates. Then the
   majors one at a time: pyo3 → rust-htslib → bio → arrow/parquet →
   noodles-gtf → statrs → bincode. After each: clippy, `cargo test`, a maturin
   build, pytest (with the binning-invariance tests). The API migrations are the real work.
4. **Output identity:** on the fully upgraded build, the RC set (28 DNA + 33
   RNA runs) must be byte-identical and the binning-invariance tests green. Any count
   change is explained (e.g. an htslib or bio behaviour change) before it is
   adopted.
5. **Policy:** for each dependency decide upgrade, cap (why), or blocked
   (issue). Add a scheduled CI job that installs the latest dependencies weekly
   and runs the suite, so a breaking release is caught before users meet it.
   Check the Docker base image and the Nextflow floor the same way.

**Acceptance.** A per-dependency table (upgraded / capped / blocked), gates
green on the upgraded set, and the RC set byte-identical.

### D5 — Coverage-driven regression panel (#155) · M
**Why.** The HPC 56-sample IMPACT matrix is hand-picked, so it covers only what
it happened to include. The operator runs the release comparison after 6.6.0
(deferred from 6.5.0), so it should be a panel chosen for coverage.

**Design.** Tag every signed-out variant in the local cBioPortal dump with the
strata that exercise distinct code paths:
- allele shape and size: SNV / DNP / TNP / ONP; insertions and deletions of 1,
  2–5, 6–19, 20–49 and 50+ bases; delins with and without a shared first base;
- repeat context of every non-SNV: homopolymer (6+), STR (3+ copies of a 2–6bp
  unit), unique;
- co-annotated structure: tract-cluster candidates (length-changing variants
  within 50bp), same-position multi-allelic rows, SNVs between two indels,
  homopolymer-twin candidates, FLT3-ITDs;
- VAF and depth bins; X, Y and MT; REF == ALT and placeholder rows;
- sample strata: assay (the IMPACT panels, IMPACT-HEME, ACCESS duplex +
  simplex), MSI-high, TMB ≥ 20.

Start from the existing acceptance samples, then let a greedy set cover pick
the fewest BAMs that bring every stratum to its target (default 25 variants,
or 4 samples for sample strata) within a run budget. Genotype each chosen
sample on all its signed-out variants, and fill a subset's tumor variants into
their matched normals. A local harness generates the selection (it names
samples, so it stays local); the operator runs it on HPC.

**Comparison.** Two per run:
- version against version, every changed cell attributed to a ticket (as the
  6.5.0 RC check did);
- gbcms against the sign-out counts, concordance by stratum. The BAM is the
  truth, so disagreements are adjudicated read by read.

**Acceptance.** A coverage table (every stratum at its target or at the most
the data has), the comparison scripts, and the panel run on HPC as the 6.6.0
release gate.

**Arms beyond the IMPACT/ACCESS panels (added 2026-09-28).** The panel BAMs are
ABRA2-realigned, deep and targeted, and the C1 work showed that each of those
properties hides behaviour other data exposes. Add:
- **RNA (FORTE):** matched IMPACT DNA and FORTE RNA for the same patients, plus
  truth-free probes at exon-boundary distances in expressed genes (REF-side
  read retention; the harness from the C1 splice check).
- **WES/WGS (TEMPO):** normals at germline hets, where VAF should read 0.5.
  TEMPO runs bwa mem without indel realignment and unbinned BQSR qualities, and
  already fills out with the C++ GetBaseCountsMultiSample (`--maq 0`).
- **Public, PHI-free:** GIAB HG002/3/4 WGS and WES BAMs (bwa, NovaSeq 4-bin
  qualities; one Oslo WES set unbinned), genotyped at truth-set indels. The
  truth sets decompose complex variants and MNPs, so complex events need
  re-merged nearby records. Results can go straight into issues and PRs.
- Report per arm: clip-covered reads outside depth, depth-only fraction in
  repeat strata, and run time at WGS scale.

### D6 — One QC-flags reference page (#156) · M
**Finding.** Every QC flag is documented in `output-formats.md`, but it is
spread out:
- the status reasons and ten `gbcms_diagnostic` flags each sit in one large
  table cell;
- the ASJD flags have their own subsection, the rescue outcomes live in the
  `gbcms_rescue` cell, and the VCF equivalents in the INFO section;
- related signals are on other pages: `rna_editing_site`, `mfsd_ch_flag`,
  `mfsd_ks_valid`, and the mFSD report's TUMOR-LIKE / CH-LIKE / INSUFFICIENT
  classes;
- the glossary has only the status table.

**Direction.** One page, `docs/reference/qc-flags.md`, with a table per family
(verdict and status reasons; `gbcms_diagnostic`; ASJD; rescue outcomes; QC
columns; mFSD classes). Each row says when the flag is set, the mode (DNA, RNA,
opt-in), its MAF column and VCF field, and what to do about it.
`output-formats.md` and the glossary link to it. A test asserts that every flag
string the code emits appears on the page, so the page stays complete.
**Measured.** All 50 flag strings the code emits are documented, across 6+ pages,
the diagnostics in one large table cell; no completeness test. Field: GATK's
header registry throws on an undocumented key; Picard and nf-core lint test
completeness; no variant tool tests flag values against its user docs.
**Decision (2026-10-05, operator): A** as directed, with the test scanning the
flag strings Rust and Python emit.

## Reviewed, no action

- vcf2maf trims case-sensitively; gbcms compares bases case-insensitively (the
  VCF spec's view). Documented.
- `vcf_alt` is each row's own allele, where vcf2maf writes the whole ALT
  column. Documented; gbcms writes one row per ALT.
- Docs examples show old version strings (`gbcms v5.3.0`); the release guide
  does not bump docs at release.

## Validation standard (every count-affecting change, before merge)

Adopted 2026-09-28. gbcms is a genotyper, not a caller: it counts reads whose own
bases carry the given allele. Alignments can be wrong (indels placed elsewhere in a
repeat, alleles clipped at read ends, long events split), so the evidence is each
read's bases, judged independently of the aligner's placement. The operator: "always
try to do things as generalized as possible; use different datasets to make the
architecture general." Each gap this cycle hid a defect only another data type
showed: realignment hid clipped carriers, Q40 bins hid the RNA BAQ loss, DNA panels
hid the RNA exon-edge collapse.

| Axis | What it exposes |
|:--|:--|
| Synthetic contracts: red-first tests; fuzz (random delins/MNP, left/right gap placement); uniform-start VAF; Q37/Q40/binned/unbinned qualities | rule logic, symmetry, quality encodings |
| Panel DNA, indel-realigned: tumour, normal, duplex, simplex | the production default |
| DNA without indel realignment: WES (the same libraries recaptured, paired with the panel), WGS | clipping, split reads, lower depth |
| RNA: spliced alignments, BAQ on; truth-free exon-edge probes | splicing, BAQ, aligner clips |
| Public reference data (GIAB WGS/WES) | shareable numbers; germline hets at 0.5 |
| Scale: WGS-sized variant lists | run time, memory |

- Truth is per read: an independent, position-aware census of each read's own
  bases across the event's repeat tract (clipped bases included, equal-length
  windows, masked low-quality bases). Cross-platform agreement is secondary:
  platforms differ by capture, depth, run, aligner and realignment.
- When the right behaviour is not known, measure it on this matrix and survey
  community practice (callers and genotypers, with sources) before deciding;
  record the decision in the tracking issue. The caveats this addresses are listed
  in `docs/reference/bam-evidence-caveats.md`.
- PRs, issues and plans carry aggregates and allele shapes only.
- Not yet covered: WGS, public reference data, other aligners and quality bins
  (D5 arms).

**Read judgment is decided spec-first (operator, 2026-10-01).** Read-judgment
tickets kept moving each other's boundaries (the same BRCA2 cluster came back in
#92, C2, C21 and C26), because each was decided against develop alone. Now:
- `docs/reference/read-judgment.md` is the spec: read shape × variant → call,
  with a decision register (RJ-n). A ticket that touches read judgment changes
  the spec first; the operator decides; code follows. A new decision names the
  ones it amends.
- `tests/test_read_judgment_spec.py` executes it (decided cells hold the rule's
  call; open cells pin today's), so a change shows in review the day it happens.
- Related tickets are decided together from that table (C25, C26, C27, C28 were).
- An existing test expectation changes only with the operator's say (AGENTS.md).

**C25 + C26 as built (2026-10-01).** Decided together with C27 and C28 from the
spec: RJ-7 (another indel inside the window: not REF, partial), RJ-8 (outside:
a separate event, REF, no partial), RJ-9 (long complex events read through the
read). C27 and C28 adopted in principle (judge their bases), to land after a
prototype is measured. Evidence: the census agrees on every synthetic shape; 105
changed DNA/WES rows: REF develop 60,619, rule 57,218, census 57,199; FORTE rows:
every row moves no more reads than its census "contradicts both" count. C25
changed nothing on real data.

**Group 1 as built (2026-10-02): what a read contributes.** RJ-10 to RJ-13
(branch `feature/g1-read-inputs`).
- C19: records without base qualities are dropped and warned; the consensus
  margin saturates. O7: unmapped records are dropped; mates of unmapped reads and
  mapped MAPQ-0 alignments still count. C29: `find_read_pos` for the anchor.
- C17 took three builds, each measured:
  - soft-clip plus masked qualities: an adversarial review found masked adapter
    filling complex windows as matches (false ALT), and RF pairs losing bases;
  - hard clip at the TLEN boundary, inward pairs only: the FLT3 ITD lost 13 ALT
    reads, because TLEN, a reference distance, leaves out inserted bases (the
    census agreed, since it shares the rule);
  - final (operator): clip only adapter-like bases (soft-clipped, or at most two
    aligned past the boundary, none inserted). One aligned base past the boundary
    mismatches the reference 75%, three 9%, ten or more 0.7%.
- Measured (develop vs final, every MAF cell; default MAPQ and `--min-mapq 0`):
  182 rows (ALT −57, REF −81), 184 at MAPQ 0. Every row is attributed to
  read-through reads. SNV rows equal the census; indel rows move toward it. All
  57 lost ALT reads had their ALT on a removed overhang; the first base past the
  boundary is the adapter's A in 97% of clipped reads. PMS2 keeps its MAPQ-0
  reads. Watchlist: base, the hard-clip build and the final build recorded.
- Known limit (second review): soft-clipped bases past a TLEN that understates
  the molecule (mate's 5' end clipped) are clipped as adapter. No count moved on
  real data (only complex/MNP windows read clips); 203 of 1,782 such reads had a
  mate 5' clip. Filed as C31 #212 (6.7.0, with C15/C7): measure from the mate's
  unclipped 5' end via the MC tag, as fgbio does.

## Triage (2026-09-30, operator)

The reviews of each ticket kept adding genuine follow-ups (17 of 50 open items
by C21). The operator went through every open issue by cluster (which items share
code or a principle, and which must come first) and cut the scope:
- **Kept in 6.6.0**, gating the cut, in the order below.
- **Moved to 6.7.0**, each with its reason on the issue:
  - evidence outside the aligned bases (C15 #173, C7 #144, C18 #177, O5 #179),
    whose prerequisite C17 stays;
  - the homopolymer twin, re-scoped as "fix or retire `--rescue-homopolymer`" now
    that O4 made it opt-in (C3 #111, #145, #146, M1 #128, umbrella #112);
  - C6 #143 (its policy follows cluster 1's quality rule), R3 #178, O3 #132, O6
    #180, I5 #127, P1 #150 (gated on T1's binning-invariance tests), and D3
    #138's migration (the pin stays).
- **Closed:** C5 #120 (measured: 0 of 50,836 signed-out insertions are 250bp or
  longer, and short reads cannot span the cap anyway; reopen for long reads); P2
  #151 (no long-pole bins in cfDNA); the umbrellas #133, #134, #135 (their items
  sit directly under #140); #92 (below).
- **Decided:** T1 retires the legacy parity path for binning-invariance tests,
  confirmed by a design note first; T2 brings the read census into the test
  suite with it (#170, #171).
- **#92's checklist** (13 silent-failure findings from an earlier audit, never
  ticked) was verified against develop: 10 were fixed by later PRs (mostly PR
  #95), so #92 is closed with each item's commit. Still live: M4 #194 (merge sums
  NA count cells as 0; decided: a missing value makes the combined cell NA, with a
  per-column warning), kept in 6.6.0; C24 #195 (stale scores in the
  local-alignment fallback tail), moved to 6.7.0.
- The deferred issues sit under a 6.7.0 tracker, #196.

## Suggested order

Refreshed 2026-10-02 (operator): the remaining work is grouped by the rule each
item touches, across 6.6.0 and 6.7.0, so related decisions are made once. Done:
cluster 1 (#203), T1 + T2 (#205), H3 (#206, #209), C25 + C26 (#210), group 1 (#211).

Each group starts from one table (affected cases, today's call, the proposal,
the evidence) and an effects map (the earlier decisions, spec cells, tests and
watchlist loci it touches). Its decisions cover the whole rule family: 6.7.0
items in the group are decided as policy now and built later under it.
Read-judgment groups extend `docs/reference/read-judgment.md`; the read census
changes in step whenever the read inputs do.

1. **What a read contributes** (bases and qualities before any rule reads them):
   C19 #182 (absent base qualities), O7 #183 (unmapped mates in mq0_count), C17
   #176 (read-through bases past the fragment end), C29 #207 (hard clips in
   `check_complex`'s query walk). Policy for 6.7.0: C15 #173 (clipped carriers),
   C7 #144 (clip-borne ITDs), C18 #177 (supplementary and split reads).
2. **Read judgment, from the spec:** C28 #202 and C27 #201 (adopted in principle:
   judge their bases), C16 #174 (stray ALT at RNA exon edges; traced first).
   Policy for 6.7.0: C24 #195 (Phase-3 fallback tail), C6 #143 (identity band).
3. **RNA strand and splices:** R4 #185 with O8 #186, R5 #198 (splices counted as
   REF coverage), C15 #173's RNA part (clips are not evidence) and C32 #213
   (spliced reads read across their junction), brought in by the operator
   2026-10-03. Policy for 6.7.0: R3 #178 (editing sites in windows).
4. **Input and representation:** C9 #122, I1 #123, I2 #124, I3 #125, I4 #126,
   #147, H2 #149, C30 #208 (lowercase or unprepared alleles).
5. **Merge, outputs, observability:** M4 #194, M5 #221, M6 #223, M1a #224, O9 #220, D3a #225, M2 #129, H1 #148, O1 #130, O2
   #131, D6 #156. Policy for 6.7.0: M1 #128, #146.
6. **Statistics:** S1 #153 (decided), S2 #154.
7. **Release:** D1 #136, D2 #137, D4 #139 before the cut, P3 #152, D5 #155 last.

Stay in 6.7.0 (no shared rule with the groups): the decomposition line (#111,
#112, #145), performance and infrastructure (#150, #127, #138), and the
diagnostics O5 #179 and O6 #180.

The 6.6.0 cut is gated on groups 1–7. The release comparison is the D5 panel run
on HPC, with the same release-candidate check as 6.5.0: every changed cell
attributed to a ticket.
