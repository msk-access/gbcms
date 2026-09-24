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

## Wrong-length pure indels (issue #91 — load-bearing)

A read whose CIGAR proves a pure indel of a DIFFERENT length at the anchor is a
**distinct allele** → `neither + has_nearby_evidence` (counts `partial_alt`),
never REF, never ALT, and **never Phase 3** — the haplotype window is
length-blind in repeat tracts, and alignment scoring promotes wrong-sequence
inserts to ALT. Same-event escapes that stay ALT:
- exact-length sequence-verified (strict/windowed),
- ≥50bp deletion **placement-aware band** (≤3 expected-span bases retained,
  ≤3 changed outside — covers wobble + split `D+M+D`; the 50bp gate is an
  artifact-SIZE prior: artifacts are small, a ≥50bp op is a real deletion),
- insertion **truncation containment** (≥4bp, ≥90% identity, both sequences
  non-low-complexity).
Same-length S3-fail / unverifiable-bases candidates keep Phase-3 arbitration
via `has_shifted_same_length` (partial propagated on non-ALT). Windowed
wrong-length (dels ≥5bp only — 1–4bp windowed Ds are noise → plain REF; ins any
size): repeat tract → neither+partial; unique context → REF+partial.
Delins stay `check_complex` (Phase-3 realignment is CORRECT for them).
Contract battery: `tests/test_wrong_length_contract.py` (parity-oracle).

## Splice-aware evidence (issue #94 cluster A — load-bearing)

A read testifies only through **aligned bases (or a D op) at the
discriminating positions**. CIGAR `N` = asserted splicing = no observation:
`splice_skip_triage` (runs before every checker) returns
`ClassifyResult::no_coverage` (`covers_locus=false`) when the N spans ALL of
them — deletion span `[pos+1, pos+ref_len)`, insertion flanks `[pos, pos+2)`,
else `[pos, pos+ref_len)` — and the engine **excludes the read from DP/DPF
entirely** (samtools-pileup semantics; mirrored binned + legacy + per-transcript;
`splice_skip_excluded=` in the Phase-stats debug line). `D` is deletion
evidence, `N` never is (no D-vs-N representation flip). Indel ops directly
after an N get the same anchor/windowed inspection as after an M (shared
helpers `resolve_anchor_{deletion,insertion}_candidate` /
`scan_windowed_*_candidate` — M-arm and N-arm must never diverge). Phase 3
never scores across a splice: `extract_raw_read_window` → `None` on N-crossing
windows; `check_complex` reconstruction → neither.
**Span-aligned REF testimony** (pure deletions only): anchor spliced out
(inside the read's N) + FULL deleted span covered by M ops + no competing
indel flags → REF, qual from the first span base. Partial span → neither.
Contract battery: `tests/test_rna_splice_contract.py`.

## Multi-Allelic Handling

At overlapping loci, sibling ALT alleles are excluded from each other's counting.

## Alignment dispatch

- Default backend PairHMM (`--alignment-backend hmm`, pangenomic WFA→PairHMM);
  SW via `sw`. Under PairHMM, SW runs only when the haplotype matrix cannot be
  built — counted (`sw_fallback_reads`), WARNed per variant, and flagged
  `SW_FALLBACK(n)`. SW penalties are constants (`SW_GAP_OPEN`/`SW_GAP_EXTEND`
  in `alignment.rs`) that CLI gap flags never reach.
- A WFA2 fast-path (`wfa_router.rs`) triages before the fallback. The fast path
  must apply the **same base-quality gate** as SW/PairHMM — it must not make a
  definitive REF/ALT call on bases the fallback would reject.

## Key Files
- `rust/src/counting/engine.rs`: binning, main loop, rayon dispatch
- `rust/src/counting/variant_checks.rs`: 4-phase check pipeline
- `rust/src/counting/alignment.rs`: Smith-Waterman
- `rust/src/counting/pairhmm.rs`: PairHMM backend
- `rust/src/counting/wfa_router.rs`: WFA2 fast-path + routing
