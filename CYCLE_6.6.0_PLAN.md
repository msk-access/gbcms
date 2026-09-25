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
> - Changes to classification, filtering or fragment consensus are mirrored
>   in the legacy parity oracle.
> - Each ticket gets an adversarial review before merge, then real-data
>   acceptance on the local truth sets. Patient data stays local; issues and
>   PRs stay PHI-free.
> - No new output columns: new signal goes to `gbcms_diagnostic`, logs, or
>   validation tooling.
> - Check `.agents/learnings/REJECTED.md` before proposing an alternative.

Priority: **H** high, **M** medium, **L** low. **[counts]** marks a ticket that
can change counts. **[decide]** marks a ticket that needs an operator decision (**[decided]**: the decision is recorded under the ticket)
before implementation.

## Summary

| ID | Ticket | Pri | Flags | Issue |
|:--|:--|:-:|:--|:--|
| C1 | Partial-ALT evidence in the SW local fallback | L | [counts] [decide] | #141 (#92) |
| C2 | REF fragments at grouped rows (main vs per-transcript) | M | [counts] [decided] | #119 |
| C3 | Homopolymer decomposition arbitration redesign | M | [counts] | #111, #145 (#112) |
| C4 | Reference windows near contig ends | M | [counts] | #142 (#92) |
| C5 | Long insertions exceed the pangenomic matrix cap | L | [counts] | #120 |
| C6 | Error-tolerant exact-length insertion matching | L | [counts] | #143 (#92) |
| C7 | Rescue for clip-borne ITD carriers | L | [counts] | #144 (#92) |
| C8 | One-base-REF delins without a shared anchor | L | [counts] | #121 |
| C9 | Count a MAF deletion at Start 1 | L | [counts] | #122 |
| C10 | Reads ending inside an indel's repeat tract counted REF | H | [counts] | #157 |
| R1 | Span-aware exon-edge BAQ rule | L | [counts] [decided] | #106 |
| R2 | RNA strandedness gating observability | M | [decided] | #114 |
| I1 | MAF allele base check | M | [decided] | #123 |
| I2 | `End_Position` optional | L | | #124 |
| I3 | VCF→MAF `Tumor_Seq_Allele1` | L | [decided] | #125 |
| I4 | maf2vcf's second ALT from `Tumor_Seq_Allele1` | L | [decided] | #126 |
| I5 | Nextflow `convert` module | L | | #127 |
| M1 | Merge rows whose flavors report different alleles | M | | #128 |
| M2 | Merge inputs from different gbcms versions | M | | #129 |
| M3 | Decomposed-allele hardening (observations, list length) | M | | #146, #147 (#112) |
| O1 | UMI warning repeated by the rescue recount | L | | #130 |
| O2 | Run-start summary of enabled options | L | | #131 |
| O3 | Rescue in fillouts without the MNP | L | | #132 |
| H1 | Writers closed when a write fails | L | | #148 (#133) |
| H2 | `is_indel` in preparation | L | | #149 (#133) |
| P1 | Deep-bin fetch reduction (M5b) | L | | #150 (#134) |
| P2 | Bin cost-sort (PF-2) | L | | #151 (#134) |
| P3 | Document the bin-span soft floor (LO-3) | L | | #152 (#134) |
| S1 | Mean LLR per fragment (CR-5) | L | [decided] | #153 (#135) |
| S2 | `MIN_FOR_KS` floor (ME-9) | L | [decided] | #154 (#135) |
| D1 | CI version-consistency check | M | | #136 |
| D2 | Release workflow creates the GitHub Release | M | | #137 |
| D3 | mkdocs-material 2.0 | L | | #138 |
| D4 | Dependency upgrade audit: does anything break on current releases? | M | | #139 |
| D5 | Coverage-driven regression panel (replaces the 56-sample matrix) | M | | #155 |
| D6 | One QC-flags reference page | M | | #156 |

## Counting correctness

### C1 — Partial-ALT evidence in the Smith-Waterman local fallback (#141, under #92) · L [counts] [decide]
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
| C | B, plus: second-try ALT calls without ALT sequence become partial, unless reads share a recurring unannotated haplotype there | Also corrects up to 35% of those ALT calls | `alt_count` changes |
| D | Leave it; document that `partial_alt` is score-based under SW | None | Known inaccuracy stays; ruled out |

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

**Effects map (B/C).**
- **Changes:** `partial_alt`, `any_alt`, `PARTIAL_DOMINANT`, VCF `PAD`/`AAD`.
  C also changes `alt_count`/`ref_count`, the fragment counts, and the
  REF/ALT-based consumers (observations export, ASJD partitions, mFSD fragment
  classes). All of it is SW backend only, and identical in the production and
  legacy paths (shared `classify_by_alignment`).
- **Unchanged:** the default `pairhmm` path; SNVs and MNPs; MNP rescue (MNPs
  only); tract-cluster claiming (complex rows exempt).

**Acceptance.** The default backend is byte-identical on the RC set. Under
SW, version-against-version on the D5 panel's complex variants, with a sample
of changed reads adjudicated in IGV.

**Priority: L** (was H). The originally planned fix doesn't improve accuracy,
doing it right needs the verification round, and it affects only a
non-default backend. C10 and C2 affect every repeat indel and every clustered
row on the default path, so they go first.

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
**Direction** (6.5.0 plan § T11).
- Arbitrate on exact haplotype support, the census method: count the reads
  that carry each candidate exactly (the called allele, both corrected shapes,
  "other").
- Report the allele the reads carry, and flag when it is neither.
- Make every consumer follow the reported allele.
- Revisit `repeat_span` for any corrected allele kept.

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

### C9 — Count a MAF deletion at Start 1 (#122) · L [counts]
**Finding.** Preparation cannot anchor a MAF deletion at Start 1 (there is no
base before it), so the row is `FETCH_FAILED`. The VCF writer writes the VCF-spec
form with the base after the event instead.
**Direction.** Resolve such rows to the same base-after form in preparation
and count them. Telomere-only.

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

## RNA

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

### I2 — `End_Position` optional (#124) · L
**Finding.** `MafReader` requires an integer `End_Position` but nothing uses
it; rows without one are skipped with a WARN. maf2vcf converts them.
**Direction.** Parse it when present; don't require it. Update the
required-columns table.

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

### M3 — Decomposed-allele hardening (#112 items 3–4: #146, #147) · M
**Finding.** The observations export cannot tell that a variant's per-molecule
rows describe the corrected allele. The `decomposed` list is not
length-validated or padded as `sibling_variants` is: a short list panics the
binned engine and silently drops rows in the legacy path.
**Direction.** Validate and pad `decomposed` like `sibling_variants`, and mark
the winning allele in the observations output (no new MAF columns).

## Observability

### O1 — UMI warning repeated by the rescue recount (#130) · L
With `--rescue-mnp`, the component recount calls the counting pass again with
the same `--umi-tag`, so the "tag never seen" WARN can repeat for one BAM.
Suppress it in the recount.

### O2 — Run-start summary of enabled options (#131) · L
One INFO block at run start naming the enabled options and what they imply
(e.g. `--rescue-mnp` replaces counts on rescued rows). A 6.5.0 T7 follow-up.

### O3 — Rescue in fillouts without the MNP (#132) · L
In a fillout of other timepoints or normals, where the MNP itself is absent,
rescue can adopt a germline component. This is documented in `--rescue-mnp`'s
help. Consider a diagnostic when the adopted component is present in the
matched normal. An upstream-annotation question, from 6.5.0 T7.

## Hygiene

### H1 — Writers closed when a write fails (#148, under #133) · L
`_write_output` does not close the writer (or its reference handle) when a
write raises. Use context managers.

### H2 — `is_indel` in preparation (#149, under #133) · L
`is_indel` reduces to `ref_len != alt_len`; its second clause is exactly
`is_mnp`. Simplify it.

## Performance (M5 leftovers)

### P1 — Deep-bin fetch reduction (M5b) (#150, under #134) · L
Deep cfDNA bins read 150k+ reads to count a few variants. Narrowing the fetch
is the only remaining cfDNA lever. It is parity-sensitive, so scope it behind
the binned↔legacy parity gate. Investigation first.

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

## Release and docs infrastructure

### D1 — CI version-consistency check (#136) · M
The 6.4.0 cut missed the Nextflow banner (`nextflow/main.nf`, still `v6.3.1`
until 6.5.0). A CI check that the release guide's 11 version references agree
would catch it.

### D2 — Release workflow creates the GitHub Release (#137) · M
6.4.0 has a tag and published artifacts but no GitHub Release page. Create the
page from the CHANGELOG section in the release workflow.

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
   build, pytest, and the parity suite. The API migrations are the real work.
4. **Output identity:** on the fully upgraded build, the RC set (28 DNA + 33
   RNA runs) must be byte-identical and the parity oracle green. Any count
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

## Reviewed, no action

- vcf2maf trims case-sensitively; gbcms compares bases case-insensitively (the
  VCF spec's view). Documented.
- `vcf_alt` is each row's own allele, where vcf2maf writes the whole ALT
  column. Documented; gbcms writes one row per ALT.
- Docs examples show old version strings (`gbcms v5.3.0`); the release guide
  does not bump docs at release.

## Suggested order

1. **Decisions: done** (2026-09-25; recorded under each ticket and on its
   issue). All the recommendations were accepted.
2. **Count-affecting, measured first:** C10 and C2 first (they share the
   discrimination window: one branch), then R2, R1, and C1 once its B-vs-C
   check is in — one branch each otherwise. C1 moved down on 2026-09-25: its
   originally planned fix did not improve accuracy on real reads, and it
   affects only the non-default SW backend.
3. **The decomposition redesign:** C3 with M1's decomposition check and M3.
4. **Hardening:** I1, I2, I3, C9, M1 (rescue conflicts), M2, O1, H1, H2.
5. **Investigations and enhancements:** C5, C6, C7, C8, P1, then P2, O2, O3,
   I4, I5, S1, S2.
6. **Infrastructure and docs:** D1, D2, D4, D5 and D6 (before the 6.6.0 cut),
   D3, P3. D4 goes early if a dependency release breaks users first.

The 6.6.0 cut is gated on steps 2–4 plus D1, D2, D4 and D6, and the release
comparison is the D5 panel run on HPC, with the same release-
candidate check as 6.5.0: every changed cell attributed to a ticket.
