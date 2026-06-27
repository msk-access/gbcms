# CONTINUITY — where we left off

> Tactical state that must survive a closed laptop or a context summary.
> Update the **Now** and **Next** sections as work progresses.

_Last updated: 2026-06-27_

## Now
Working the code-review remediation plan (`CODE_REVIEW_IMPLEMENTATION_PLAN.md`)
ticket by ticket, one PR each, with review-before/after discipline and real-data
validation on MSK slice + ctDNA-duplex BAMs.

**M1 (count-correctness) — COMPLETE & merged:**
CR-1 (#25 bin anchor coverage), CR-3 (#26 empty-allele guards), CR-2 (#27 WFA BQ
gate), CR-4 (#28 tolerant-deletion partial seq-check). Plus #29 Docker CI
resilience. All independently cross-checked with pysam (ALT 571/571 exact).

**M2 (statistical integrity) — in progress:**
CR-5 (closed-form LLR, removes ±∞; exact small-N KS, Rust-only, no scipy) on
branch `feature/cr-5-mfsd-llr-ks`. Real-data: 0 non-finite LLR, all KS p in [0,1].

## Next
1. Open/merge the CR-5 PR.
2. Continue M2: HI-10 (report mislabels the fragment-size LLR as the PairHMM LLR),
   HI-11 (no FDR on the 6 KS p-values), ME-8 (ASJD BH padding), ME-9 (MIN_FOR_KS
   vs report min_alt), ME-10 (sub_nuc NaN), ME-11 (report isinf guard).
3. Then M3 (RNA): HI-7 gene_strand, HI-8 splice fragment size, HI-9 motif strand, …

## Open follow-ups (tracked, not lost)
- BAM-level binned↔legacy parity CI gate, incl. large deletions (the gap that hid CR-1).
- Prep-time empty-allele validation (loud, once-per-variant) — complements CR-3.
- HI-3: WFA off-target global-edit-distance + fixed threshold.
- LO-2: feature-gate the legacy per-variant traversal (keep as parity oracle, exclude
  from the shipped wheel) — **revisit after this plan**, per maintainer.

## Key decisions
- Stats stay in **Rust** (KS/LLR/Fisher in `mfsd.rs`/`shared/stats.rs`); **no scipy**
  dependency — exact KS is a self-contained Rust DP, validated against baked
  SciPy reference constants.
- Legacy `count_bam` (per-variant) is kept for now as the binned↔legacy **parity
  oracle**; production uses `count_bam_binned` only.
- Tests kept minimal/high-signal (each is a maintenance contract); fixes that reduce
  duplication are preferred over adding code.
