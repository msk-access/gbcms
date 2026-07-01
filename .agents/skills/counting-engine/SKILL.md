---
name: counting-engine
description: Reference for the gbcms counting core — genomic binning (10kb bins, one fetch per bin), the 4-phase variant-check pipeline, and alignment-backend dispatch. Use when modifying engine.rs or variant_checks.rs, debugging count discrepancies, or reasoning about binning and binned↔legacy parity. For fragment/consensus see fragment-counting; for filters/quality see read-filters-qc.
---

# Counting Engine Patterns

## Genomic Binning

Variants are grouped into ~10kb bins for efficient BAM traversal:
- `BIN_WINDOW = 10_000` — one `bam.fetch()` per bin (not per variant)
- `BIN_MAX_VARIANTS = 200` — forced split when exceeded
- Padding: `max(repeat_span + 2, 5)` to capture nearby evidence
- Rayon `par_iter()` parallelizes across bins

**Invariant (load-bearing):** a bin's fetch-end must cover the *anchor* (leftmost)
variant's full ref span, not just `bin_start + window`. A bin anchored by a large
deletion/DelIns whose `ref_allele.len()` exceeds the window will otherwise
under-fetch its right tail and undercount AD/ADF. Any binning change needs a
binned↔legacy parity test including large deletions and complex DelIns.

## Variant Check Pipeline

```
Phase 1: Simple match (exact REF/ALT at position)
Phase 2: Windowed scan (±5bp, expanding for repeats)
Phase 2.5: Edit distance fallback (Levenshtein, >1 edit margin)
Phase 3: Full alignment (SW or PairHMM, via wfa_router fast-path)
```

## Multi-Allelic Handling

At overlapping loci, sibling ALT alleles are excluded from each other's counting.

## Alignment dispatch

- Default backend SW (`--alignment-backend sw`); PairHMM via `hmm`.
- A WFA2 fast-path (`wfa_router.rs`) triages before the fallback. The fast path
  must apply the **same base-quality gate** as SW/PairHMM — it must not make a
  definitive REF/ALT call on bases the fallback would reject.

## Key Files
- `rust/src/counting/engine.rs`: binning, main loop, rayon dispatch
- `rust/src/counting/variant_checks.rs`: 4-phase check pipeline
- `rust/src/counting/alignment.rs`: Smith-Waterman
- `rust/src/counting/pairhmm.rs`: PairHMM backend
- `rust/src/counting/wfa_router.rs`: WFA2 fast-path + routing
