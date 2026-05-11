# Counting Engine Patterns

## Genomic Binning

Variants are grouped into ~10kb bins for efficient BAM traversal:
- `BIN_WINDOW = 10_000` — one `bam.fetch()` per bin (not per variant)
- `BIN_MAX_VARIANTS = 200` — forced split when exceeded
- Padding: `max(repeat_span + 2, 5)` to capture nearby evidence
- Rayon `par_iter()` parallelizes across bins

## Variant Check Pipeline

```
Phase 1: Simple match (exact REF/ALT at position)
Phase 2: Windowed scan (±5bp, expanding for repeats)
Phase 2.5: Edit distance fallback (Levenshtein, >1 edit margin)
Phase 3: Full alignment (SW or PairHMM)
```

## Fragment Consensus

- Always on — quality-weighted R1/R2 consensus
- QNAME hashing (u64) for fragment tracking
- Discarded ambiguous fragments counted in DPF, not RDF/ADF
- `--library-type amplicon` bypasses fragment consensus

## Multi-Allelic Handling

At overlapping loci, sibling ALT alleles are excluded from each other's counting.

## Key Files
- `rust/src/counting/engine.rs`: binning, main loop, rayon dispatch
- `rust/src/counting/variant_checks.rs`: 4-phase check pipeline
- `rust/src/shared/fragment.rs`: FragmentEvidence, consensus
- `rust/src/counting/alignment.rs`: Smith-Waterman
- `rust/src/counting/pairhmm.rs`: PairHMM backend
