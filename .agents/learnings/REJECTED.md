# REJECTED

The rejected-edit buffer (negative feedback). When a proposed change is vetoed
("no, leave that", "we decided against that"), log it here so no future session
re-proposes the dead end. Before promoting any learning, grep this file for the
target + topic; if a matching veto exists, surface "vetoed YYYY-MM-DD because X"
and ask whether something changed, rather than re-proposing.

Newest at the top.

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
