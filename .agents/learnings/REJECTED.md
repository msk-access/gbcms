# REJECTED

The rejected-edit buffer (negative feedback). When a proposed change is vetoed
("no, leave that", "we decided against that"), log it here so no future session
re-proposes the dead end. Before promoting any learning, grep this file for the
target + topic; if a matching veto exists, surface "vetoed YYYY-MM-DD because X"
and ask whether something changed, rather than re-proposing.

Newest at the top.

---

## [REJ-20261023-001] Choose the guard's rule by whether the ALT is present; a cis rescue for SNP-only reads (#294)
- **Target:** the exact-carrier rule's repair of the neighbouring-change guard at a one-unit
  insertion in a VNTR whose units differ at single bases (RJ-27 amendment of 2026-10-20)
- **Proposed:** (1) choose the masking rule by whether the ALT is present at the locus; (2)
  a cis rescue: count SNP-only reads as ALT where other fragments link the SNP to the
  insertion.
- **Reason vetoed:** (1) the ALT's presence is not observable at the VNTR, and the rule lost
  the ALT in 6 contract tests. (2) It credits a read on other fragments' evidence, so under
  invariant 7 it could only ever be an opt-in rescue; it would gain at most 6 ALT at HG002
  chr22:39653046 and nothing at chr21:41865095 or at the #294 row, where no fragment links
  the SNV to a base that tells the insertion apart. Not built as a rule. (2) was built once, as
  a temporary measurement option on #264 (a masked change read as the reference's on REF's
  haplotypes where REF's molecules lack it): at HG002's two duplications it counted 18 carriers,
  10 of them on other fragments' evidence; the option was removed on 2026-10-31. Each option
  offered was measured before it was decided; the record is
  https://github.com/msk-access/gbcms/issues/294#issuecomment-6364311473 ("Considered and
  not adopted"), and the numbers are in docs/reference/read-judgment.md (evidence behind the
  RJ-27 amendment of 2026-10-20 to 2026-10-22).
- **Date:** 2026-10-23 (the decision record; the options were measured 2026-10-17 to 2026-10-22)

## [REJ-20261020-001] Leave a guard mask unapplied where it would leave no base telling the alleles apart (#294 option 1)
- **Target:** the exact-carrier rule's neighbouring-change guard (RJ-23, RJ-27) at a
  one-unit insertion in a VNTR whose units differ at single bases
- **Proposed:** a mask that would leave the windows no base telling the alleles apart is
  not applied: the change is part of what tells them apart (decided 2026-10-17 on #294,
  extending "the given allele's own bases are never masked").
- **Reason vetoed:** it counts a germline SNP's haplotype as the insertion (operator,
  2026-10-20, replaced by G3: the windows read on, the placed-apart check, reads read
  whole). Where the SNP turns the next unit into the inserted one, both haplotypes read
  alike over a short read, so un-masking credits the SNP's reads as ALT. HG002's 28-bp
  VNTR insertion (SNV in trans): RD/AD 10/0 → 1/8, both phased ALT reads from the SNV's
  haplotype (G3: 0/0, no read reaches the telling base). BWA-placed synthetic VNTR (28-bp
  units): the SNP alone, no insertion, 68/0 → 34/31 (all 31 ALT false; G3 34/0); the
  insertion and the SNP in trans 73/0 → 0/67 (G3 0/0); the insertion alone 68/0 → 27/34
  (G3 the same). Frozen build `plan670/src_c44u`; shapes pinned by
  `tests/test_masked_repeat_reads_contract.py`. Breaks invariant 7.
- **Date:** 2026-10-20

## [REJ-20261019-001] Count a read carrying the given repeat-unit change plus more units as ALT (and partial)
- **Target:** pure indels in repeats: a read of another unit count than the given one
  (today partial_alt only, never AD; option A of #285)
- **Proposed:** where the read's bases carry the given unit change and further units of
  the same repeat, count it in AD and also in partial_alt (conservative: no ALT read lost).
- **Reason vetoed:** keep AD exact (operator, 2026-10-19). On HG002 truth hets it gains
  almost nothing (upper bound, every partial read added: mean |VAF−0.5| STR deletions
  0.090 → 0.092, STR insertions 0.093 → 0.094, 20+ bp insertions 0.159 → 0.151); it
  credits ALT in every sample carrying another length of the repeat (germline STR lengths
  in normals at fill-out; #294's SNV-only VNTR shape); "that exact unit" is undefined in a
  perfect repeat; STR genotypers (HipSTR, GangSTR, DRAGstr) weigh ±1-unit reads in a
  stutter model, not in the allele's count. Breaks invariant 7.
- **Date:** 2026-10-19

## [REJ-20261018-001] Accept another breakpoint that gives the reads' bases to within one base (#285 first posting)
- **Target:** RJ-24's run clause: a run of changes past a long deletion's junction that
  rejoins the flank
- **Proposed:** a run is a phased neighbour unless another breakpoint gives the reads'
  consensus bases "one base aside" (posted in the refinement of 2026-10-18).
- **Reason vetoed:** operator correction, 2026-10-18
  (https://github.com/msk-access/gbcms/issues/285#issuecomment-6252631315). An adversarial
  review measured the cost: across 40 random contexts a phased DNP on the first two bases
  past the junction counted ALT in only 21; in the other 19 the deletion shifted right
  matched the reads at all but one base. With an exact match it counts ALT in 36 of 40 (the
  four that do not: the deletion shifted right exactly, 2, and DNPs covering every base
  where the alleles' windows differ, 2). The exact rule's cost: another breakpoint that also
  carries an SNV inside the run counts as the given deletion with a neighbour. A second
  review (F1, #296 folded in,
  https://github.com/msk-access/gbcms/issues/285#issuecomment-6258175767) found a
  same-length deletion 3 bp right counted ALT in 31 of 150 random contexts when each change
  was masked alone; an allele's recurring changes are now checked together against other
  breakpoints. Option A of the thread is REJ-20261019-001.
- **Date:** 2026-10-18

## [REJ-20261016-002] The CIGAR path's identity band for truncated inserts under the exact-carrier rule (#278 D2)
- **Target:** pure insertions judged by their bases (RJ-20, RJ-27): a read whose insert is a
  shorter slice of the ALT's
- **Proposed:** keep the 90% identity band (4 bp or more, at least 90% identity to the best
  window of the expected insert, both sequences non-low-complexity;
  `insert_truncation_match`) so a truncated copy counts ALT.
- **Reason vetoed:** operator, 2026-10-16 (D2;
  https://github.com/msk-access/gbcms/issues/278#issuecomment-6223388424; the same figures
  closed #143, https://github.com/msk-access/gbcms/issues/143#issuecomment-6225682238). It
  added 158 ALT reads on HG002 long insertions with truth concordance unchanged, and on the
  panel 13 ALT calls the bases oracle agrees with against 111 it does not. A truncated copy
  is another allele (partial evidence). The band stays only where the rule cannot judge (no
  prepared reference, the CIGAR path). #143 (error-tolerant insertion matching) was
  resolved with RJ-26 in #289.
- **Note:** this supersedes the 2026-09-22 recommendation to keep the band as primary
  (REJ-20260922-001).
- **Date:** 2026-10-16

## [REJ-20261016-001] Grow the pure-indel windows through repeats on both alleles, or not at all (#278)
- **Target:** RJ-27's window growth
- **Proposed:** growth on both alleles, or no growth.
- **Reason vetoed:** the three variants were measured on the same inputs (operator,
  2026-10-16; https://github.com/msk-access/gbcms/issues/278#issuecomment-6223388424).
  Chosen: growth only on the allele holding the indel's bases (ALT for insertions, REF for
  deletions). HG002 long insertions, mean |VAF − 0.5|: .1498 chosen, .1510 both alleles,
  .1540 none. HG002 deletions in band: 563, 563, 571. Oracle-correct panel reads, REF:
  489,324, 488,127, 492,962; ALT: 151,958, 151,533, 154,045. Reads of another allele counted
  ALT: 897, 851, 1,166. No growth reads more oracle-correct panel reads (REF 492,962, ALT
  154,045) and 571 deletions in band, but counts 1,166 reads of another allele as ALT
  against 897 and is farther from 0.5 at long insertions (.1540 against .1498). Growth on
  both alleles counts fewer reads of another allele as ALT (851) but is farther from 0.5 at
  long insertions (.1510) and reads fewer oracle-correct panel reads (REF 488,127, ALT
  151,533).
- **Date:** 2026-10-16

## [REJ-20261014-001] An unreadable ALT core always carries only the length (#277 first cut, and its own-base variants)
- **Target:** the exact-carrier rule's reading of an ALT core whose bases are below min BQ
  (RJ-20), and the MNP check's partly read reads
- **Proposed:** a carrier that reads none of the ALT's own bases turns partial at every
  locus (c504a454); variants: own bases count only where they differ from REF's at that
  read position (71edff1a), and the same for REF (b77873c9).
- **Reason vetoed:** the locus-aware rule removes the same false ALT without the cost
  (operator, 2026-10-14). The unconditional form demoted carriers where no other allele
  of the ALT's length exists: synthetic no-competitor shapes lose 286 carriers vs 202,
  normals lose 13 reads an independent check calls ALT, HG002 delins shifts; it left the
  MNP check's phantom ALT untouched (natural quality: 15 phantom reads where the allele
  is absent). Locus-aware: 0 false ALT in every competitor and phantom shape, no-
  competitor loss as develop, normals/WES/HG002 unchanged. Frozen builds
  `plan670/src_c43a` (first cut), `src_c43b`, `src_c43c`, `src_c43m` (locus-aware).
- **Date:** 2026-10-14

## [REJ-20261013-002] Keep the CIGAR checks for pure indels; let bases overrule disagreeing or clipped reads (C43 hybrid)
- **Target:** pure-indel read judgment (`check_insertion`, `check_deletion`)
- **Proposed:** keep the CIGAR verdict, and hand a read to the exact-carrier rule only
  when the CIGAR path does not call it ALT and it shows 2+ confident changes around the
  event or a soft clip there. It kept every decided spec case.
- **Reason vetoed:** a partial fix (operator, 2026-10-13, chose every indel by its bases,
  C43 #276). HG002: insertions 10+ bp 287/603 in band (bases rule 332), insertions 1-9
  bp 518/669 (534), REF still 1.37x and 1.22x the REF reads: the REF credited to reads
  ending 1-4 bases past a tract is untouched, and a duplicated copy aligned onto the
  next copy often shows one mismatch, under the trigger. Two rules per read.
- **Date:** 2026-10-13

## [REJ-20261013-001] Grow a long deletion's continuation until four bases tell the alleles apart; let the guard read on to tell a germline run from another breakpoint
- **Target:** the C40 junction judge and its guard (RJ-24)
- **Proposed:** for a germline change on one of the few bases telling the alleles
  apart: the continuation grows from 4 until 4 of its bases differ between REF and ALT
  (at most 16), and the guard assigns a read to REF when it fits REF's deleted bases 20
  bases past the windows, masking a run there as the REF haplotype's own.
- **Reason vetoed:** no gain on any arm (a re-test of REJ-20261009-002's diff rule with
  the guard added): HG002 131 -> 129/195 in 0.35-0.65, panel 160 -> 159/196 within 10%
  of the reads, WES and normals unchanged; one panel row gained 58 REF reads. The panel
  shortfall it aimed at was mostly the sibling guard reading the whole deleted span
  (fixed in C40, 5aaf76a4). Prototype builds `plan670/src_pj4`, `src_pj5`.
- **Date:** 2026-10-13
- **Revisited in part (2026-10-19, #285):** the guard now tells a fragment's allele by
  the 12 bases past the windows (with a shift test, reading on through 24), for either
  allele, without growing the continuation, and masks a run past the junction on
  either allele's reads only where its bases rejoin the allele's continuation and no
  other breakpoint gives them. Measured gain where this entry found none: a panel
  167-bp deletion whose REF molecules carry a germline SNP on the first deleted base
  went RD 117 -> 319 (sign-out 352). Decided by the operator as part of #285 (runs
  inside the window); the continuation growth stays rejected.

## [REJ-20261009-003] Read-on for the longer allele at 50+ bp deletions (C40 option A, RJ-9 style)
- **Target:** the C40 junction judge for 50+ bp pure deletions (RJ-24)
- **Proposed:** as decided in the 6.7.0 plan: REF at the POS breakpoint must also be
  clean through every deleted base the read covers past the equal window.
- **Reason vetoed:** it checks REF over 50-300 bases and ALT over about 4, so ordinary
  errors cost REF only. Measured (bases-only oracle, k = 4): HG002 50+ bp het
  deletions in 0.35-0.65 100/166 with read-on vs 103/166 without (mean |VAF-0.5| 0.168
  vs 0.159, median 0.538 vs 0.515); IMPACT within 10% of sign-out 29/60 vs 40/60 (RD
  +14% without). Adopted instead: both alleles judged on the equal window only
  (option B). Equal exposure (both alleles checked over everything they cover, R2)
  was middling: HG002 104/166, panel 34/60, and it hurt delins 5-19 (99 vs 105/147).
  Delins (RJ-9) show the same direction; that is a measured follow-up, not vetoed.
- **Date:** 2026-10-09

## [REJ-20261009-002] A dynamic k at 50+ bp deletions from a reference signature
- **Target:** the C40 junction judge (k, the clean bases past the microhomology)
- **Proposed:** raise k where the sequence signals it: until the REF and ALT
  continuations differ at m positions (diff2/3/4), until the junction has no other
  placement near the breakpoint (uniq), or both; floored at 4 so it only rises.
- **Reason vetoed:** no gain against truth. HG002 50+ bp het deletions, the same 191
  rows: fixed k = 4 in 0.35-0.65 131 (mean |VAF-0.5| 0.132); diff3 127 (0.138); uniq
  129 (0.137); diff3+uniq 126 (0.141). Panel: 30-31/60 for every rule. A larger k only
  costs reads. The error a dynamic k aimed at was read-on's one-sided exposure
  (REJ-20261009-003), not k; tolerant read-on under a low-complexity signature (R3)
  did not help either (HG002 low-complexity rows 12/28, as without it).
- **Date:** 2026-10-09

## [REJ-20261009-001] k = 10 for both alleles where the discriminating bases are low complexity (C40 guard)
- **Target:** the C40 junction judge
- **Proposed:** in the 6.7.0 plan's option C: a low-complexity far flank (the engine's
  test: the longest tandem repeat covers >= 80% of the stretch) raises k from 4 to 10
  for both alleles.
- **Reason vetoed:** worse against truth. On the 39 HG002 deletions it applies to, k =
  10 put 12/31 rows in 0.35-0.65 vs 15/32 at k = 4 (overall 144/230 vs 147/231); it
  never applies on the IMPACT panel. Those deletions read high (median 0.61-0.64) at
  any k: 27% of their REF reads show a confident mismatch in the repetitive deleted
  bases (9% elsewhere) and truth has no variant there, so the cause is error
  exposure in repeats (C41 #264), not chance matches.
- **Date:** 2026-10-09

## [REJ-20261008-001] Padded whole windows, whole alleles without padding, or junctions from either flank for complex events (#253)
- **Target:** how a complex event's windows are read (RJ-23, C39)
- **Proposed:** the alternatives measured beside the decided one-sided, equal-length
  windows: (a) padded whole windows (6.6.0: the shorter allele's window padded to equal
  length); (b) whole alleles without padding; (c) junctions read from either flank (6.6.0's
  rule past 50 bases).
- **Reason vetoed:** operator, 2026-10-08
  (https://github.com/msk-access/gbcms/issues/253#issuecomment-6073704478; correction
  https://github.com/msk-access/gbcms/issues/253#issuecomment-6074829582). Padding costs
  reads, not bias: on the panel's 137 delins rows one-sided windows judge about 46% more ALT
  and 41% more REF at the same VAF (×0.99-1.00), and sign-out agreement goes 0/14 → 9/14. On
  HG002 truth (400 germline het delins, span 20+) one-sided reads VAF 0.470-0.481 with 25%
  more ALT reads, against 6.6.0's 0.464. The alternatives are biased as read-start counts
  predict: whole alleles without padding favour the shorter allele (panel VAF ×1.22),
  junctions from either flank favour the longer one (×0.87). Community practice: no tool
  pads windows; the likelihood callers (GATK, Strelka2, Octopus) count a read once its
  discriminating bases at one breakpoint are read; only GangSTR, ExpansionHunter and
  Varlociraptor correct spanning bias, and with a model. The claim that 6.6.0's past-50-base
  rule reads VAF about 0.65× of one-sided was withdrawn on 2026-10-09: the trace logs reads
  the engine then never counts; on the four long delins 6.6.0's written REF is 1,356 against
  the one-sided reading's 1,351 (VAF ratio 1.03).
- **Date:** 2026-10-08

---

## [REJ-20261007-002] Credit ALT from placements scored by readable mismatches (C37 option C′)
- **Target:** insertion read judgment — `insert_placements` / the windowed scan
  (`rust/src/counting/variant_checks.rs`), RJ-21.
- **Proposed:** score every junction a read's same-length insert could sit at by its
  readable mismatches (not only those a fits-gated slide reaches) and call ALT when a
  minimum-cost placement shows the ALT.
- **Reason vetoed:** operator, after the C37 measurement (#245): at a BRCA2 insertion row (G>GA)
  the reads it would turn from REF into ALT carry a recurring other allele (G>A plus an
  A insertion, 11+ reads in two samples); AD would rise 23→28 and 28→40. A per-read
  rule cannot tell a singleton error from a recurring allele, so placement scoring may
  decide "another allele versus separate event", never ALT.
- **Date:** 2026-10-07

## [REJ-20261007-001] Judge the strict path's junction insert over its best-aligned placements (C37 option B)
- **Target:** `resolve_anchor_insertion_candidate` (`rust/src/counting/variant_checks.rs`)
  — the junction ALT with a flank substitution a slide absorbs into another insert.
- **Proposed:** make the strict path judge like C36's windowed scan (only the
  fewest-mismatch placements count), so the junction ALT gives way to the slid insert.
- **Reason vetoed:** operator, after the C37 measurement (#245): all 28 RC reads are
  the given ALT plus a single-molecule substitution (24 in no other read, 5 an SNV on
  the haplotype, no other allele recurring); ABRA2 put them at the junction because a
  singleton gets no contig of its own; GATK, ABRA2, bcftools and the literal counters
  all credit them ALT. The change would drop real carriers (BRCA2 AAG AD 13→10, 43→39).
- **Date:** 2026-10-07

## [REJ-20261006-001] Give ALT the REF margin: no ALT call on a read whose last base is the deciding base
- **Target:** pure-indel read judgment (RJ-1/RJ-3; `tests/census.py` policy rules) — the
  one-anchor ALT verdict.
- **Proposed:** require one base past the deciding base for ALT as well as REF, so a
  terminal sequencing error cannot turn a REF read into an ALT read.
- **Reason vetoed:** operator, after the C35 BAM check (#242): at the 18 changed RC rows,
  reads ending on the deciding base show a third base (always an error) 0 times in 348
  ACCESS reads and once in 71 IMPACT reads, predicting ~0 and ~0.5 false ALT against 32
  and 13 terminal ALT reads. The REF margin guards a systematic aligner bias (an ALT read
  near its end written REF with a cheaper mismatch); a false terminal ALT needs a random
  Q20+ error, the same one-base evidence an SNV ALT call rests on. A symmetric margin
  would drop ~45 real carriers there. Kept; the third-base control becomes a validation
  check, not a rule.
- **Date:** 2026-10-06

## [REJ-20261002-001] A strict (C22-style) rule and the pure-indel ALT-by-bases rule for a read deleting the anchor (#202)
- **Target:** RJ-15: a read whose own deletion covers a pure indel's anchor
- **Proposed:** two prototypes. Strict, as C22: an unmasked base wherever the alignment's
  haplotype differs from the ALT. The pure-indel ALT-by-bases rule, which stops at the first
  deciding base.
- **Reason vetoed:** operator, 2026-10-02, who adopted a third prototype
  (https://github.com/msk-access/gbcms/issues/202#issuecomment-5974777463). Measured against
  develop and the read census on the same inputs: strict, ALT −26, 22 of them reads with a
  masked base at the deleted anchor that the census and develop count ALT; ALT-by-bases,
  false ALT at insertion rows (a BRCA2 row 9 → 35, census 17). Adopted: ALT net 0 on RC DNA
  (±1 read at 6 rows, all masked-base edge cases) and −2 on WES; REF +57 at the BRCA2
  cluster, every row toward the census; summed distance to the census on the changed indel
  rows REF 344 → 251, ALT 40 → 44. The old route, Phase 3's closer-haplotype credit, is gone
  from this case (a fuzz of 2,400 random compositions found 88 ALT calls on reads whose
  bases are neither allele). Evidence: docs/reference/read-judgment.md (RJ-14 and RJ-15).
- **Date:** 2026-10-02

## [REJ-20261001-001] Keep the 1-4bp "noise" policy for another indel in a pure indel's window (#200 option b)
- **Target:** which of a read's other indels decide its REF call (RJ-7, RJ-8)
- **Proposed:** keep the policy that a 1-4bp windowed deletion failing S3 is alignment noise
  ("CIGAR definitive"), so such a read stays REF, and make the post-walk rule consistent
  with it (no flip to neither + partial when an unrelated op appears).
- **Reason vetoed:** operator, 2026-10-01, decided in #210
  (https://github.com/msk-access/gbcms/issues/200#issuecomment-5950921993). A read that
  deletes or inserts inside the discrimination window does not hold REF there, so under
  "count the given allele" it is another allele. The read census, which judges bases only,
  agrees: on the 105 changed DNA and WES rows it counts 57,199 REF reads, develop 60,619,
  the decided rule 57,218; 155 rows change (REF −4,327, partial +4,324; ALT and depth
  unchanged). Under the kept policy each slippage allele's carriers add to the other rows'
  RD at the same tract: synthetic (REF, ALT, partial) of (10, 0, 0) for −A carriers on a +A
  row, +A carriers on a −A row, a D(2) inside the span of a unique −8 row, and a D(1) right
  after the anchor of a unique +10bp insertion. The rule is stricter than GATK's AD by
  decision.
- **Date:** 2026-10-01

## [REJ-20260925-003] A pre-filter STRAND_DISCORDANT; documenting rna_antisense_depth as "reads that reached classification" (#114)
- **Target:** RNA strandedness observability: `rna_antisense_depth` and ASJD
  `STRAND_DISCORDANT`
- **Proposed:** (1) leave `rna_antisense_depth` as the reads that reached classification and
  only document it; (2) measure `STRAND_DISCORDANT` on all junction reads, before the strand
  filter.
- **Reason vetoed:** operator, 2026-09-25
  (https://github.com/msk-access/gbcms/issues/114#issuecomment-5834318189,
  https://github.com/msk-access/gbcms/issues/114#issuecomment-5834502322). (1) At defaults
  the column was always 0 while the census shows 143 antisense reads (0.5%) at 13 of the 71
  FORTE truth loci (18%), 24 carrying the ALT; counting before the filter, like `mq0_count`,
  keeps the column meaningful and changes no REF or ALT count. (2) Before any filter,
  antisense junction evidence is about 0.1% (125,541 sense and 138 antisense junction reads
  near the loci), so a pre-filter flag would be close to noise; it stays a
  `--no-strandedness` diagnostic. Numbers and standards: docs/reference/read-filters.md
  (strandedness filter).
- **Date:** 2026-09-25

## [REJ-20260925-002] A read-level BAQ rule: no N-penalty at annotated junctions (#106 option B)
- **Target:** the RNA exon-edge BAQ exception (R1)
- **Proposed:** skip the N-penalty in the heuristic BAQ for Ns whose ends are annotated
  junctions (within the junction tolerance) and keep the indel penalties, instead of a
  variant-level rule on the REF span's distance to the nearest boundary.
- **Reason vetoed:** it re-penalizes novel junctions near a variant, which is ASJD's
  splice-creation signal, restores indel penalties at edges, and is a larger change to the
  main counts; the 6.5.0 plan rejected it unless the span-aware rule (option A) left residual
  edge masking, and at acceptance A left none: 135 of the 139 moved rows equal what the
  reads' own bases predict, and the other 4 were inspected read by read (decided 2026-09-25,
  https://github.com/msk-access/gbcms/issues/106#issuecomment-5834502040; accepted
  https://github.com/msk-access/gbcms/issues/106#issuecomment-5896307927). Archived plan
  `.agents/archive/CYCLE_6.5.0_PLAN.md` § T9, option B. Numbers:
  docs/reference/rna-splice-handling.md (why the exception reads the whole REF span).
- **Date:** 2026-09-25

## [REJ-20260925-001] Flag from the local scores; an 8-mer content rule; that rule with demotion; leave it documented (#141 options A to D)
- **Target:** the Smith-Waterman local fallback's partial-ALT flag and ALT calls for complex
  variants
- **Proposed:** A: thread the second try's (local) scores into the tail (the original plan).
  B: flag only when the read carries ALT-specific sequence (an 8-mer content rule). C: B
  plus turning the second try's ALT calls without ALT sequence into partial. D: leave it and
  document that `partial_alt` is score-based under SW.
- **Reason vetoed:** operator, 2026-09-25, who replaced all four with the exact-carrier rule
  for complex variants
  (https://github.com/msk-access/gbcms/issues/141#issuecomment-5839231028). A: flagging from
  the local scores was right for 32 of 67 REF/tie reads against 36 as it stood; both flag on
  score closeness, not read content. B: the 8-mer rule would demote 6 of the 23 true ALT
  carriers (one ALT call without an ALT-specific 8-mer was judged ALT by the census), and
  the first round's 22 of 62 second-try ALT calls (35%) without ALT sequence came largely
  from non-spanning reads and detector misses
  (https://github.com/msk-access/gbcms/issues/141#issuecomment-5837829090). C: rested on B's
  detector. D: keeps a known inaccuracy. Of the 110 fallback reads linked exactly, 77 (70%)
  did not span the event; the rule keeps 90% of ALT calls and 82% of REF calls on the
  default backend, VAF +3.8% median. Archived plan `.agents/archive/CYCLE_6.6.0_PLAN.md` §
  C1 (options table); numbers in docs/reference/read-judgment.md (evidence behind RJ-4).
- **Date:** 2026-09-25

## [REJ-20260924-001] Splice-aware Phase 3 across junctions (B2b of #94); hard-filtering junction-adjacent variants
- **Target:** RNA: scoring a read across a splice junction against a spliced haplotype (a
  coordinate-mapped spliced-haplotype Phase 3)
- **Proposed:** B2b: a genomic→spliced coordinate map with junction-compatible extraction,
  so Phase 3 scores junction reads against a spliced haplotype. Weighed beside it: hard-
  filtering junction-adjacent variants, as several pipelines do.
- **Reason vetoed:** closed as not demanded (2026-09-24,
  https://github.com/msk-access/gbcms/issues/94#issuecomment-5815889505). Measured demand
  across the RNA truth cohort was zero: every covered junction-adjacent indel classified
  correctly through the structural paths (post-N inspection, windowed scan, span-aligned REF
  testimony) (https://github.com/msk-access/gbcms/issues/94#issuecomment-5779396531). No
  mainstream tool scores reads across junctions (GATK's SplitNCigarReads plus
  HaplotypeCaller cuts reads at every N; long-read RNA pipelines document junction-adjacent
  indels as a weak spot and several hard-filter splice-adjacent variants; RNA allele
  counters do not support indels), so B2b would be novel territory
  (https://github.com/msk-access/gbcms/issues/94#issuecomment-5779885342). Hard-filtering
  such variants, as those pipelines do, is not an option for a genotyper: what cannot be
  judged is exposed through diagnostics, never by filtering. The consensus intron snipping
  that preceded it had failed on consensus instability and the missing coordinate map and
  was removed. Design note and reopen condition (a delins carrier at a junction in truth
  data): docs/reference/rna-splice-handling.md. RJ-19 later built read-derived spliced
  windows for the exact-carrier rule; that is not B2b.
- **Date:** 2026-09-24

## [REJ-20260923-006] First designs of the cluster contest: a phase-rank gate, a per-sibling span cost with tie-demotion, ALT-versus-REF alone (#92, #99)
- **Target:** exclusive assignment at co-annotated tract clusters (6.5.0; the first form of
  RJ-5)
- **Proposed:** (1) a phase-rank gate: demote an ALT read by the phase that called it; (2) a
  per-sibling span cost, a sibling with the lower cost claiming the read and ties demoted;
  (3) ALT-versus-REF alone as the test.
- **Reason vetoed:** each was falsified by real-data acceptance before the shipped
  three-test contest (https://github.com/msk-access/gbcms/issues/92#issuecomment-5802243054;
  archived plan `.agents/archive/CYCLE_6.5.0_PLAN.md` § T1): (1) cross-attribution on
  ABRA-realigned consensus reads flows through Phase-3 PairHMM, not the windowed scan, and
  the rank inverts on delins; (2) anchor-exact siblings stole delins carriers, which zeroed
  complex rows; (3) same-tract same-length deletions and ladders both beat REF. No counts
  were recorded for these three. The canonical-identity design of the adjudication
  (https://github.com/msk-access/gbcms/issues/92#issuecomment-5788067318) was the starting
  point, replaced by an engine-level sibling contest because the contest covers delins. The
  shipped contest's acceptance is in docs/reference/read-judgment.md (evidence behind RJ-5's
  cluster contest).
- **Date:** 2026-09-23

## [REJ-20260923-005] Make dynamic_sw_gap_extend real, a locus-adaptive Smith-Waterman gap-extend (#92)
- **Target:** the Smith-Waterman gap-extend penalty (`dynamic_sw_gap_extend` in
  `rust/src/counting/alignment.rs`)
- **Proposed:** make the intended tight-to-free relaxation for deep repeat tracts real,
  mapping the curve so large spans reach 0.
- **Reason vetoed:** the SW fallback never fires on real data. Traced runs on an ACCESS
  duplex sample and the heaviest MSI-high sample (42 tract indels; 31 Phase-3
  classifications) made zero Smith-Waterman fallback invocations; per the measurement
  comment, under the default PairHMM backend the pangenomic matrix builds whenever prep
  supplies `ref_context`, and on the explicit SW backend the function was a documented
  constant −1 anyway (the logistic is capped at 0.5, so it rounded to −1 for every
  `repeat_span`). Making it real would change SW alignment scores on a path that never
  fires, at the cost of a validation pass on the sw backend and the matrix-build-failure
  fallback. Replaced by a documented constant, no behaviour change
  (https://github.com/msk-access/gbcms/issues/92#issuecomment-5788174498; landed in #102).
- **Date:** 2026-09-23

## [REJ-20260923-004] Drop the contig name from MNP rescue labels
- **Target:** `src/gbcms/pipeline.py` — `RESCUED_COMPONENT(...)` / `gbcms_rescue` labels
- **Proposed:** write `pos(REF>ALT)` only, since the row names the contig.
- **Reason vetoed:** operator — labels follow the output file's own contig naming.
- **Date:** 2026-09-23

## [REJ-20260923-003] Return the MNP rescue gate to the original strict `ad == 0` (T7 "path B")
- **Target:** `src/gbcms/pipeline.py` — `_rescue_mnp_pass` candidate gate
- **Proposed:** keep only the rescue bug fixes, drop the partial-dominance gate.
- **Reason vetoed:** operator chose path A (simplify the gate to "rescue only when no read
  shows the whole MNP") — the TERT C250T-shaped rows must remain rescuable.
- **Date:** 2026-09-23

## [REJ-20260923-002] Error allowance `ceil(partial_alt × 10^(−min_baseq/10))` for confirmed whole-MNP reads
- **Target:** `src/gbcms/pipeline.py` — `_confirmed_error_allowance`
- **Proposed:** tolerate that many fully-read full-MNP reads as sequencing error.
- **Reason vetoed:** real ACCESS data: BRCA2 AAG>TAC duplex rescued to a germline SNP
  (6 real whole-MNP reads ≤ allowance 8); KRAS ACC>CCA simplex rescued with 1 real
  whole-MNP read (allowance rounds up to 1). Replaced by "no read shows the whole MNP".
- **Date:** 2026-09-23

## [REJ-20260923-001] Report-only MNP rescue (component counts only in `gbcms_rescue`, row keeps MNP counts)
- **Target:** `src/gbcms/pipeline.py` — `_rescue_mnp_pass`
- **Proposed:** never replace counts; put the component split in the audit or flags.
- **Reason vetoed:** operator — it would get lost; rescue is opt-in and exists to put the
  component count in the count columns, flagged (`RESCUED_COMPONENT`) and warned.
- **Date:** 2026-09-23

## [REJ-20260922-002] A repeat-aware gate for windowed deletions under 5bp (#92 class 2)
- **Target:** the windowed deletion scan's `<5bp` noise rule (`variant_checks.rs`)
- **Proposed:** also flag a sub-5bp windowed deletion that falls inside the variant's tract
  interval, so the existing repeat resolution applies (tract: neither + `partial_alt`;
  unique flank: plain REF), as for insertions of any size (proposed 2026-09-18, #92).
- **Reason vetoed:** not needed, measured on real MSI-high truth data: 6 MSI-high samples
  (MSI scores 45-53) carrying 246 signed-out small indels at classic MSI target tracts
  (https://github.com/msk-access/gbcms/issues/92#issuecomment-5780293921). Median
  ad-versus-sign-out relative error +1.4% (IQR +0.6% to +2.9%) across all 247 covered loci;
  35/247 loci show 3 or more same-tract wrong-length ops, and across them 481 stutter reads
  are flagged `partial_alt` against 35 not flagged (93%); residual unflagged stutter about
  0.14 reads per locus; off-length sub-5bp ops 11.7 per 1000 dp at MSI-high loci against 0
  at the MSS control. Left-alignment collapses same-tract ops onto the shared anchor, where
  the anchor-position wrong-length rule already fires, so real stutter is not displaced from
  it. Documented in docs/reference/allele-classification.md (note after the windowed-scan
  safeguards).
- **Date:** 2026-09-22

## [REJ-20260922-001] A BQ-first policy for exact-length insertions with mismatches, in place of the identity band (#92 class 2)
- **Target:** insert matching in the 6.5.0 design (the >=90% identity band,
  `insert_truncation_match`)
- **Proposed:** match inserts only on quality-masked bases (BQ-first), with no identity
  band.
- **Reason vetoed:** measured on IMPACT raw reads and ACCESS consensus reads
  (https://github.com/msk-access/gbcms/issues/92#issuecomment-5779396851,
  https://github.com/msk-access/gbcms/issues/92#issuecomment-5787996247): the dominant
  non-exact class is zero-mismatch truncations (about 24% of carriers), where no quality
  policy acts; mismatching candidates are about 2% of carriers, mostly Q20-35 singletons at
  unique positions that BQ-first would have demoted; the class where both policies agree
  (all-low-BQ) is 6 reads. On ACCESS consensus reads the classes were exact, length
  truncations (2-14 per locus), N-carrying (0-10 duplex, up to 83 simplex) and zero
  same-length mismatches, so BQ-first had nothing to act on. The band was kept. Numbers:
  docs/reference/bam-evidence-caveats.md (sequencing errors inside or beside the event).
- **Revisited in part (2026-10-16, #278 D2):** for prepared insertions the exact-carrier
  rule judges quality-masked bases with one substitution tolerated (RJ-26), and the identity
  band itself was then tested and not adopted (REJ-20261016-002); the band survives only in
  the CIGAR fallback. The both-policies principle (annotation-side against read-side errors)
  stands as history; the input allele is now taken as correct (invariant 7).
- **Date:** 2026-09-22

## [REJ-20260917-001] The reciprocal-overlap tier and the split-suspect escape hatch to Phase 3 for 50+ bp deletions (#91)
- **Target:** wrong-length deletions at the anchor (`check_deletion`,
  `rust/src/counting/variant_checks.rs`)
- **Proposed:** keep (a) the CR-4 tier, a 50% reciprocal overlap with a sequence check
  counting as ALT, and (b) an escape hatch sending a split-suspect read to Phase 3 for
  net-matching.
- **Reason vetoed:** (a) at a shared anchor the tier admitted any deletion sharing half the
  length (a D(60) passed for a 100 bp variant) while its sequence check compared reference
  against reference, identically true by construction; on the large-deletion panel every
  length 40-120bp counted as full ALT, and read-level validation found zero length wobble on
  nine real pure deletions of 12-539bp. (b) removed after the phase-2 review: the
  first-draft band summed net deleted bases anywhere in the region, so a D(2) at the anchor
  + M(53) + D(49) read passed for an expected D(50), and Phase 3 with narrow padding
  promoted such not-the-event reads to ALT
  (https://github.com/msk-access/gbcms/issues/91#issuecomment-5721200517). The
  placement-aware band that replaced them was itself retired by RJ-24 (deletions of 50 bp or
  more are judged by their bases at the POS breakpoint). Treating delins as in scope of the
  wrong-length rule was also rejected then (regime 2); complex events now go through the
  exact-carrier rule (RJ-4).
- **Date:** 2026-09-17

---

## [REJ-20260701-001] Gate structural-ALT INDEL win on the REF mate's base quality (ME-12b)
- **Target:** `rust/src/shared/fragment.rs` — `FragmentEvidence::resolve` (~:206)
- **Proposed:** before a structural ALT (matching CIGAR I/D op) wins a fragment,
  require it not be contradicted by a high-BQ REF on the other mate — i.e. compare
  the structural ALT against the REF mate's base quality.
- **Reason vetoed:** for INDELs, the REF mate reports **anchor** base quality (the
  base before the ins/del), which measures "how confident is the anchor base call,"
  **not** "how confident is the INDEL detection." The two are orthogonal, so gating
  the structural ALT on that BQ compares a semantically meaningless quantity. The
  CIGAR I/D op is the only signal that discriminates the alleles. The unconditional
  structural-priority rule was validated on a DNMT3A duplex dataset (all 7 conflict
  fragments were genuine D-op evidence, MAPQ=60, correct length). Adding the gate
  would discard true-positive INDEL fragments. Phase-3 (complex / wrong-length)
  classifications are `is_structural=false` and already fall through to
  quality-weighted arbitration, so the concern doesn't apply there either.
- **Note:** ME-12 part (a) — relabeling the stale "Majority Rule" comment — WAS done.
- **Date:** 2026-07-01

---

## [REJ-20260627-001] HI-5: strip N bases from the read before the PairHMM (to make N "LLR-neutral near indels")
- **Target:** `rust/src/counting/pairhmm.rs` (`classify_by_marginalized_pairhmm` / `BQEmission`)
- **Proposed:** the M4 scope suggested removing N positions from the read (or otherwise
  tweaking N emission) so an N at an indel discriminating locus contributes a strictly
  zero LLR.
- **Reason vetoed:** an empirical probe disproved the premise. With the existing 0.25
  N emission **and** the marginalized (forward-sum) PairHMM, N at a deletion locus
  already yields LLR ≈ **0.41** — far below the 2.3 decision threshold, so it stays
  ambiguous (and N at a SNP locus is exactly 0). **Stripping the N makes the read
  mimic the deletion haplotype → LLR ≈ 8.3 → misclassifies as ALT.** So the proposed
  fix regresses, not improves. HI-5 is already handled correctly; the action taken was
  a regression test (`test_n_base_neutral_at_snp_and_indel`) that guards the behavior,
  not a code change. (HI-4 clamp and ME-5 `prob_emit_y=1.0` were the real fixes in the
  bundle.)
- **Date:** 2026-06-27
- **Amended 2026-10-22 (operator):** the premise was measured on a flankless G run, through
  a helper scoring with the read's end gaps free (the PairHMM's orientation then). With
  flank bases, as prep's contexts have, an N at the deleted base of a REF-length read leans
  REF by its length (LLR about -7.1, before and after the read was aligned whole): a base
  the deletion has none of. That stands; the veto stands (stripping the N gives +8.6, ALT).
  The test now asserts both on the production path.

## [REJ-20260626-001] Adopt the infra-bound harness layers now (agents, fleet, second brain, multi-provider, crons)
- **Target:** the harness as a whole (`docs/harness/00-plan.md`)
- **Proposed:** build out claudelicious docs 07–19 alongside the spine.
- **Reason vetoed:** these are **infra-bound** — their value scales with running
  services / defined long-horizon jobs, not codebase maturity. Building them with
  nothing to hold is the "half-used scripts" anti-pattern. Promote one only when a
  concrete trigger appears (e.g. an always-on agent *iff* we decide to grind the
  M1–M6 backlog autonomously). Context-bound layers (memory, learning loop, rules)
  were adopted instead.
- **Date:** 2026-06-26
