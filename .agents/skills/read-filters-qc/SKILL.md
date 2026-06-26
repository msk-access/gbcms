---
name: read-filters-qc
description: Reference for gbcms read-level (genomic) QC — BAM flag filters (duplicate/secondary/supplementary/qc-failed/improper-pair/indel), MAPQ and base-quality gating, heuristic BAQ near indels/splice junctions, and QC diagnostics (mq0_count, n_count, alt/ref distance-to-end). Use when changing filtering or quality gating. NOTE: this is genomic QC — for the dev lint/test sweep see the qa-check skill.
---

# Read Filters & Quality Gating (genomic QC)

## Phase 0 universal filters (`shared/filters.rs`)

`ReadFilter::passes()` checks, in order: duplicate → secondary → supplementary →
qc-failed → improper-pair → indel(CIGAR I/D). Rejections are tallied in
`FilterCounts` for debug logging.

**Defaults** (CLI): `filter_duplicates`, `filter_secondary`, `filter_supplementary`,
`filter_qc_failed` = **true**; `filter_improper_pair`, `filter_indel` = **false**.
(Nextflow `nextflow.config` mirrors these — keep them in sync.)

**MAPQ is NOT in `filters.rs`.** It is mode-specific and stays in the engine:
RNA does NH:i:1 rescue of low-MAPQ unique reads; `mq0_count` tracks MAPQ==0 reads.

## Base-quality gating

- `min_baseq` (default 20) gates discriminating bases.
- SW masks sub-`min_baseq` bases to `N`; PairHMM weights by BQ. The WFA fast-path
  must apply the **same** gate (don't make a definitive call on bases the fallback
  would reject).
- `usable_count >= 3` high-quality bases required for a confident alignment call.

## Heuristic BAQ (`shared/baq.rs`, Li 2011)

- Downgrades BQ within `BAQ_RADIUS = 5` bp of any Ins/Del/RefSkip by
  `BAQ_PENALTY = 20` (clamped at 0). Allocation-free when the read has no indel/N.
- For RNA, the RefSkip (N) penalty substitutes for GATK SplitNCigarReads overhang
  clipping.
- Default: **DNA off, RNA on** (`--apply-baq` / `--no-baq`). The Nextflow modules
  encode this asymmetry — DNA passes `--apply-baq` only if true; RNA passes
  `--no-baq` only if false.

## QC diagnostic outputs

- `mq0_count`: reads with MAPQ==0 at the locus.
- `n_count`: reads with an N at ≥1 discriminating position (NAD) — duplex masking / failure.
- `alt_dist_end_median` / `ref_dist_end_median`: median base distance to read end.

## Key Files
- `rust/src/shared/filters.rs`: `ReadFilter`, `FilterCounts`
- `rust/src/shared/baq.rs`: `apply_heuristic_baq`
- `rust/src/counting/engine.rs`: Phase 0 wiring, MAPQ/mq0, BQ gate
- `tests/test_filters.py`
