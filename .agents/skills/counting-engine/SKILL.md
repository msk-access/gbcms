---
name: counting-engine
description: "Map of the gbcms counting core: how check_allele_with_qual dispatches a read (splice triage, SNV/MNP checks, the exact-carrier rule, the 50 bp deletion junction judge, CIGAR checks as the residual fallback), the key files, and genomic binning. Use when modifying engine.rs, carrier.rs, junction.rs or variant_checks.rs, debugging a count discrepancy, or reasoning about binning invariance. For fragment/consensus see fragment-counting; for filters/quality see read-filters-qc."
---

# Counting Engine Map

What call each read shape gets, and the decision behind it (register RJ-n), is specified in
`docs/reference/read-judgment.md`; how the engine computes it is
`docs/reference/allele-classification.md`. Read those first; this page maps the code.
A fix to a call starts from the mechanism, not the failing row: state it without the row,
test the other shapes of its class, derive any cutoff from it, and justify any narrowing
(read-judgment.md, "How a change to read judgment is made").

## Dispatch: `check_allele_with_qual` (engine.rs), in order
1. **Splice triage** (`splice_skip_triage`, variant_checks.rs): a CIGAR `N` spanning every
   discriminating position gives `no_coverage`; the read is left out of DP/DPF.
2. **`window::allele_kind(ref, alt)`** picks the route from the alleles, never from the
   `variant_type` label: Snv, Mnp, Insertion, Deletion or Complex.
3. **SNV** → `check_snp`. **MNP** → `check_mnp` (masked per position, no alignment); an MNP
   read with an indel or clip at the block takes the complex route (6).
4. **Pure insertion, and a pure deletion step 5 does not take** → the exact-carrier rule (`by_bases` →
   `carrier::classify`): a read counts when its own bases carry the allele.
5. **Prepared deletion of 50 bp or more** (`junction::judges`, `LONG_DELETION`) →
   `junction::classify`: the read's bases at the POS breakpoint, however the aligner wrote
   the gap.
6. **Complex** (delins, Del+SNV, an insertion that changes the anchor) → `classify_complex`,
   which is the exact-carrier rule too.
7. **Residual fallback.** Only where a rule cannot judge (no prepared reference, or one
   without the event's flank) do the CIGAR checks run: `check_insertion`, `check_deletion`,
   `check_complex` (windowed scan, edit distance, then "Phase 3": WFA triage, then PairHMM or
   Smith-Waterman). Such reads set `carrier_fallback` and are counted and warned once per
   variant. The fallback decides essentially no reads: Phase 3 decided 0 of ~4M read calls
   on every runnable arm (#265).

The `read call` trace line names the deciding rule (`deciding_rule`: `exact_carrier`,
`junction`, `carrier_fallback`, `alt_by_bases`, ...).

## Key files (`rust/src/counting/`)
- `engine.rs`: binning (`build_genomic_bins`), the shared read loop, dispatch, fragment
  resolution, ASJD, per-transcript counts; `ref_needs_the_window` / `alt_needs_the_window`.
- `carrier.rs`: the exact-carrier rule, the neighbouring-change guard (`Guard`,
  `tells_others_aligned`), `sibling_by_bases`. `junction.rs`: the long-deletion judge.
- `window.rs`: `allele_kind`, read/guard/discrimination windows, informative-read tests.
- `observed.rs`: `observed_allele` (the OBSERVED_ALLELE diagnostic; changes no count).
- `variant_checks.rs`: the CIGAR checks and windowed scans. `utils.rs`: `ClassifyResult`,
  `ClassifyPhase`, haplotype and masked-compare helpers.
- `alignment.rs` (SW), `pairhmm.rs`, `pangenome.rs`, `wfa_router.rs`. The CLI default
  backend is `pairhmm`; all backends and the WFA fast path share one base-quality gate.

## Binning (load-bearing: AGENTS.md invariant 1)
Variants group into bins of at least `BIN_WINDOW` (10 kb) and at most `BIN_MAX_VARIANTS`
(200); one `bam.fetch()` per bin, rayon over bins. Bin geometry is performance only: each
bin's fetch must hold every member's full read window, anchor included, and its guard
window (`window::guard_window`). Keep `tests/test_binning_invariance.py` and the
`build_genomic_bins` property test green. A co-annotated sibling's ALT read is not REF for
this row, and a tie between rows is neither row's AD (RJ-5).
