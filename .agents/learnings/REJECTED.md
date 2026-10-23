# REJECTED

The rejected-edit buffer (negative feedback). When a proposed change is vetoed
("no, leave that", "we decided against that"), log it here so no future session
re-proposes the dead end. Before promoting any learning, grep this file for the
target + topic; if a matching veto exists, surface "vetoed YYYY-MM-DD because X"
and ask whether something changed, rather than re-proposing.

Newest at the top.

---

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

---

## [REJ-20261007-002] Credit ALT from placements scored by readable mismatches (C37 option C′)
- **Target:** insertion read judgment — `insert_placements` / the windowed scan
  (`rust/src/counting/variant_checks.rs`), RJ-21.
- **Proposed:** score every junction a read's same-length insert could sit at by its
  readable mismatches (not only those a fits-gated slide reaches) and call ALT when a
  minimum-cost placement shows the ALT.
- **Reason vetoed:** operator, after the C37 measurement (#245): at BRCA2 13:32906888
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
