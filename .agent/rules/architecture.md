---
description: gbcms architecture reference — Python/Rust boundary, module table, mFSD/Parquet design, algorithmic features, invariants
alwaysApply: true
---

# gbcms Architecture Reference

## System Overview

gbcms is a **Python/Rust hybrid** tool for counting alleles at variant positions across BAM files. Python handles CLI, orchestration, and I/O; Rust handles all BAM traversal, read classification, fragment tracking, Fisher tests, mFSD analysis, and Parquet writing.

## Layer Diagram

```
┌─────────────────────────────────────────────────────────────┐
│                     CLI Layer (Typer)                       │
│  cli.py — argument parsing, 4-layer validation, logging     │
└──────────────────────────┬──────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────┐
│                   Orchestration Layer                        │
│  pipeline.py — sample iteration, progress, Parquet dispatch  │
└──────────────────────────┬──────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────┐
│                       I/O Layer                              │
│  io/input.py  — VcfReader, MafReader, ReferenceChecker      │
│  io/output.py — VcfWriter, MafWriter  (mFSD-gated columns)  │
└──────────────────────────┬──────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────┐
│                     Core / Models                            │
│  core/kernel.py — CoordinateKernel (VCF/MAF → 0-based)      │
│  models/core.py — GbcmsConfig, OutputConfig, AlignmentConfig │
└──────────────────────────┬──────────────────────────────────┘
                           │
┌─────────────────────────────────────────────────────────────┐
│                   Rust Engine (gbcms._rs)                    │
│                                                             │
│  counting/                                                   │
│    engine.rs         — Main loop, DP gating, read iteration  │
│    variant_checks.rs — check_snp/mnp/ins/del/complex         │
│    alignment.rs      — Smith-Waterman backend                │
│    pairhmm.rs        — PairHMM backend (--alignment-backend) │
│    mfsd.rs           — Mutant Fragment Size Distribution      │
│    utils.rs          — ClassifyResult, haplotype helpers     │
│    parquet_writer.rs — write_fsd_parquet() via ZSTD Parquet  │
│                                                             │
│  shared/                                                      │
│    fragment.rs       — FragmentEvidence, QNAME/UMI hashing   │
│    stats.rs          — Fisher's exact test (strand bias)     │
│    bam_utils.rs      — median_qual, find_read_pos            │
│    filters.rs        — ReadFilter struct (BAM flag checks)   │
│    baq.rs            — BAQ heuristic (Li 2011)               │
│                                                             │
│  normalize/                                                  │
│    engine.rs         — Normalization pipeline                │
│    left_align.rs     — bcftools-style left-alignment         │
│    decomp.rs         — Homopolymer decomposition detection   │
│    fasta.rs          — Reference sequence fetcher            │
│    repeat.rs         — Tandem repeat detection               │
│                                                             │
│  types.rs            — Variant, BaseCounts (PyO3 bindings)   │
└─────────────────────────────────────────────────────────────┘
```

## Python/Rust Boundary Rules

| Rust ✓ | Python ✓ |
|--------|----------|
| BAM traversal (rust-htslib) | CLI / Typer commands |
| Read classification (all variant types) | Config validation (Pydantic) |
| Fragment tracking (QNAME hashing) | Orchestration / progress (Rich) |
| Fisher's exact test | VCF/MAF I/O |
| mFSD analysis (KS test, LLR) | Workflow coordination |
| Native Parquet writing (Arrow + ZSTD) | Logging setup |
| Normalization (left-align, decomp) | |
| Rayon parallelism per-variant | |

## Module Reference

| Module | Purpose |
|:-------|:--------|
| `cli.py` | CLI entry, 4-layer validation, logging setup |
| `pipeline.py` | Orchestration, progress, Parquet dispatch |
| `normalize.py` | Standalone normalization workflow (no counting) |
| `io/input.py` | VcfReader, MafReader, ReferenceChecker |
| `io/output.py` | VcfWriter, MafWriter with mFSD column gating |
| `core/kernel.py` | CoordinateKernel — VCF/MAF → 0-based conversion |
| `models/core.py` | GbcmsConfig, OutputConfig, AlignmentConfig (Pydantic) |
| `utils/logging.py` | Structured logging setup |
| `rust/counting/engine.rs` | Main counting loop, DP gating, Rayon parallelism |
| `rust/counting/variant_checks.rs` | check_snp/mnp/ins/del/complex, windowed scan |
| `rust/counting/alignment.rs` | Smith-Waterman implementation |
| `rust/counting/pairhmm.rs` | PairHMM backend, LLR scoring |
| `rust/counting/mfsd.rs` | Fragment size distribution analysis (KS test, LLR) |
| `rust/counting/parquet_writer.rs` | write_fsd_parquet(), Arrow/ZSTD native Parquet |
| `rust/shared/fragment.rs` | FragmentEvidence, quality-weighted consensus, QNAME hashing |
| `rust/shared/stats.rs` | Fisher's exact test (strand bias) |
| `rust/shared/bam_utils.rs` | median_qual, find_read_pos |
| `rust/shared/filters.rs` | ReadFilter struct, FilterCounts (universal BAM flag checks) |
| `rust/shared/baq.rs` | Heuristic BAQ (Base Alignment Quality, Li 2011) |
| `rust/normalize/` | Left-alignment, decomp, fasta, repeat (7 modules) |
| `rust/types.rs` | Variant, BaseCounts PyO3 bindings |

## Key Design Decisions

1. **Rust for counting**: rust-htslib for BAM; Rayon for per-variant parallelism.
2. **0-based internal coordinates**: 1-based in VCF/MAF externally; converted at boundary.
3. **mFSD is opt-in** (`--mfsd`): Writers gate 34 MAF cols and 7 VCF INFO fields behind `self.mfsd`. When off, columns are **absent** (not NA-filled).
4. **Rust-native Parquet** (`--mfsd-parquet`): `write_fsd_parquet()` in Rust via `arrow`/`parquet` crates with ZSTD(1). No `pyarrow`. `ref_sizes`/`alt_sizes` are **Rust-internal** — no `#[pyo3(get)]`.
5. **ZSTD requires explicit feature**: `parquet = { default-features = false, features = ["arrow", "zstd"] }`. SNAPPY and ZSTD are both disabled by `default-features = false`.
6. **4-layer CLI validation**: Parse-time → Pre-model → Model-time → No silent skips.
7. **Fragment counting always on**: Quality-weighted consensus; discards on ambiguity counted in DPF not RDF/ADF.
8. **Windowed indel detection**: ±5bp scan (expanding to `max(5, repeat_span + 2)`) with 3-layer safeguards.
9. **Dual alignment backends**: SW (default) or PairHMM (`--alignment-backend hmm`).
10. **Parquet dispatched from pipeline.py**: `write_fsd_parquet()` is called after `count_bam()`, not inside the counting engine.

## Critical Algorithmic Features

| Feature | Location | Description |
|:--------|:---------|:------------|
| Windowed indel scan | `variant_checks.rs` | ±5bp + repeat expansion, 3-layer safeguards |
| Masked complex comparison | `variant_checks.rs` | Quality-aware masking, 3-case comparison |
| Phase 2.5 edit distance | `variant_checks.rs` | Levenshtein fallback, >1 edit margin guard |
| Phase 3 alignment | `alignment.rs`, `pairhmm.rs` | Dual-trigger SW or PairHMM |
| Fragment consensus | `fragment.rs` | u64 QNAME hash, quality-weighted, discard ambiguous |
| Multi-allelic exclusion | `engine.rs` | Sibling ALT exclusion at overlapping loci |
| Adaptive context | `normalize/repeat.rs` | Dynamic padding in tandem repeat regions |
| Strand bias | `stats.rs` | Fisher's exact test at read + fragment level |
| mFSD analysis | `counting/mfsd.rs` | KS test, LLR, pairwise comparisons across 4 classes |
| Native Parquet | `parquet_writer.rs` | Arrow ListArray for variable-length size arrays |

## Configuration Model

`GbcmsConfig` (Pydantic) contains:
- `OutputConfig`: `output_dir`, `format`, `suffix`, `column_prefix`, `preserve_barcode`, `show_normalization`, **`mfsd`**, **`mfsd_parquet`**
- `ReadFilters`: duplicate/secondary/supplementary/qc/improper/indel flags
- Thresholds: `min_mapq` (20), `min_baseq` (20), `fragment_qual_threshold` (10), `context_padding` (5)
- `AlignmentConfig`: `backend` (sw|hmm), `llr_threshold`, gap probabilities

## Test Suite (15 files, 107+ tests)

| File | What It Tests |
|:-----|:-------------|
| `test_accuracy.py` | SNP, insertion, deletion, complex, MNP accuracy; DP invariant |
| `test_shifted_indels.py` | Windowed indel detection (15 cases) |
| `test_fuzzy_complex.py` | Quality-aware masked complex matching (14+ cases) |
| `test_filters.py` | Read filter flags |
| `test_strand_counts.py` | Strand-specific counting |
| `test_alignment_backend.py` | SW vs PairHMM concordance |
| `test_fragment_consensus.py` | Fragment quality consensus, DPF invariant |
| `test_multi_allelic.py` | Sibling ALT exclusion |
| `test_dp_neither.py` | DP includes neither/third-allele reads |
| `test_normalization.py` | Left-alignment, REF validation, window expansion |
| `test_cli_sample_id.py` | CLI argument parsing, validation, error paths |
| `test_maf_reader.py` | MAF input parsing |
| `test_maf_preservation.py` | MAF column preservation |
| `test_pipeline_v2.py` | End-to-end integration |
| `test_mfsd_flag.py` | mFSD flag gating, column counts, Parquet output (8 tests) |

## Invariants (Always Assert in Tests)

- `dp >= rd + ad` — DP includes 'neither' reads
- `dpf >= rdf + adf` — DPF includes discarded ambiguous fragments
- `rd == rd_fwd + rd_rev` — strand consistency
- `ad == ad_fwd + ad_rev` — strand consistency
- MAF column count: **145** without `--mfsd`, **179** with `--mfsd` (+34)
- VCF INFO fields: exactly 7 `MFSD_*` fields added only when `--mfsd` is set
