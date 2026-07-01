---
description: gbcms architecture — project layout, Rust/Python boundary, module map
alwaysApply: true
---

# Architecture

## Project Layout

```
gbcms/
├── src/gbcms/           # Python package
│   ├── cli.py           # Typer CLI (DNA + RNA commands)
│   ├── pipeline.py      # Orchestration, progress, Parquet dispatch
│   ├── normalize.py     # Standalone normalization workflow
│   ├── _rs.pyi          # Primary Rust type stubs (authoritative)
│   ├── io/
│   │   ├── input.py     # VcfReader, MafReader, ReferenceChecker
│   │   └── output.py    # VcfWriter, MafWriter (mFSD/RNA column gating)
│   ├── models/
│   │   └── core.py      # GbcmsConfig, OutputConfig, AlignmentConfig (Pydantic)
│   ├── report/
│   │   └── mfsd_report.py  # HTML mFSD report generator
│   └── utils/
│       └── logging.py   # Structured logging setup
├── src/gbcms_rs.pyi     # Legacy stub (must stay synced with _rs.pyi)
├── rust/                # Rust crate (gbcms_rs)
│   └── src/
│       ├── lib.rs       # PyO3 module entry
│       ├── types.rs     # Variant, BaseCounts PyO3 bindings
│       ├── counting/
│       │   ├── engine.rs        # Main counting loop, genomic binning, Rayon par_iter()
│       │   ├── variant_checks.rs # check_snp/mnp/ins/del/complex, windowed scan
│       │   ├── alignment.rs     # Smith-Waterman
│       │   ├── pairhmm.rs       # PairHMM backend, LLR scoring
│       │   ├── pangenome.rs     # Haplotype matrix for complex phase
│       │   ├── wfa_router.rs    # WFA2 fast-path alignment
│       │   ├── rna.rs           # RNA validation, strandedness, splice junctions, editing sites
│       │   ├── mfsd.rs          # Fragment size distribution (KS test, LLR)
│       │   ├── parquet_writer.rs # Arrow/ZSTD native Parquet
│       │   └── fragment.rs      # Re-export shim for shared::fragment
│       ├── annotation/          # v5.0.0: GTF-informed annotation
│       │   ├── mod.rs           # AnnotationIndex (COITree, splice sites, transcript introns)
│       │   └── gtf.rs           # GTF parser (variant-guided streaming)
│       ├── shared/
│       │   ├── fragment.rs      # FragmentEvidence, QNAME/UMI hashing
│       │   ├── stats.rs         # Fisher's exact test, strand bias
│       │   ├── bam_utils.rs     # median_qual, find_read_pos
│       │   ├── filters.rs       # ReadFilter, FilterCounts
│       │   └── baq.rs           # Heuristic BAQ (Li 2011)
│       └── normalize/           # Left-alignment, decomp, fasta, repeat
├── nextflow/            # Nextflow pipeline wrapper
├── tests/               # 255 Python tests (15+ files)
├── docs/                # MkDocs documentation
└── .agents/             # Agent rules, workflows, skills
```

## Rust/Python Boundary

| Rust ✓ | Python ✓ |
|--------|----------|
| BAM traversal (rust-htslib) | CLI / Typer commands |
| Read classification (all variant types) | Config validation (Pydantic) |
| Fragment tracking (QNAME hashing) | Orchestration / progress (Rich) |
| Fisher's exact test | VCF/MAF I/O |
| mFSD analysis (KS test, LLR) | Workflow coordination |
| Native Parquet writing (Arrow + ZSTD) | Logging setup |
| Normalization (left-align, decomp) | HTML report generation |
| Rayon parallelism per-bin (10kb windows) | |
| GTF annotation index (COITree) | |
| ASJD detection | |

## Key Design Decisions

1. **Rust for counting**: rust-htslib for BAM; Rayon for per-**bin** parallelism (`par_iter()` over ~10kb genomic bins).
   - **`--threads` is the TOTAL thread budget per process.** Multi-sample parallelism is Nextflow's job (gbcms runs as N concurrent processes, each pinned to `task.cpus`), so every parallel section must stay within `--threads` — all rayon pools are sized from `shared::resolve_thread_budget(threads)` (which also guards `num_threads(0)`=all-cores), and any future htslib decode threads must **subdivide** this budget, never add to it. No `par_iter` may run on rayon's global pool.
2. **0-based internal coordinates**: 1-based in VCF/MAF externally; converted at boundary.
3. **mFSD is opt-in** (`--mfsd`) — gated at *both* layers (output-aware engine, invariant #3). Writers gate 34 MAF cols / 7 VCF INFO fields behind `self.mfsd` (absent when off, not NA-filled); the **binned engine also gates the compute** — `mfsd` is plumbed `OutputConfig.mfsd → count_bam_binned → count_variant_from_cache`, and when off the per-fragment size arrays, the `compute_mfsd_stats` stats, and the post-counting mFSD BH-FDR pass are all skipped (no compute-then-discard, no held `ref_sizes`/`alt_sizes`). The legacy parity oracle always computes mFSD (mFSD ∉ `PARITY_FIELDS`, so this never breaks parity); both paths share `compute_mfsd_stats`.
4. **Rust-native Parquet** (`--mfsd-parquet`): `write_fsd_parquet()` via `arrow`/`parquet` crates with ZSTD(1). No `pyarrow`.
5. **4-layer CLI validation**: Parse-time (Typer) → Pre-model (cli.py) → Model-time (Pydantic) → No silent skips.
6. **Fragment counting always on**: Quality-weighted consensus; discards counted in DPF not RDF/ADF.
7. **Windowed indel detection**: ±5bp scan expanding to `max(5, repeat_span + 2)`.
8. **Dual alignment backends**: SW (default) or PairHMM (`--alignment-backend hmm`).
9. **Genomic binning**: ~10kb bins, one `bam.fetch()` per bin, max 200 variants/bin.
10. **COITree for annotation**: Platform-portable metadata access via `Borrow` trait (nosimd vs NEON/AVX backends). The tree *layout* is arch-specific, so the **GTF disk cache (`--gtf-cache-dir`, M5a) never serializes the trees** — it persists only the parsed intermediate (`GtfIndexBundle`: exon records, splice sites, introns, chrom map; `bincode`, version-tagged) and rebuilds the trees via `build_exon_trees` on load. Caching is best-effort (missing/corrupt/stale/unwritable → log + plain parse); keyed on GTF identity + variant chroms. Concurrent cohorts must pre-warm via `gbcms build-gtf-cache` (else the first wave all cold-miss); the per-sample `count_bam_binned` runs then load in ~0.05s instead of re-parsing (~9s).
11. **Diagnostic flags**: `gbcms_diagnostic` and `gbcms_rescue` are strongly-typed Rust fields, not dynamic attributes.

## Legacy `count_bam` parity oracle (`legacy-parity` feature)

**What it is.** `count_bam` (in `engine.rs`, with its helper `count_single_variant`)
is a *second, independent* implementation of the counting logic that fetches reads
**per variant**. Production never calls it — the pipeline only uses `count_bam_binned`
(one `bam.fetch()` per ~10kb bin). `count_bam` exists **solely as the parity oracle**:
the tests cross-check that the optimized binned path produces identical counts to the
straightforward per-variant path.

**Build matrix.** Gated behind the default `legacy-parity` Cargo feature
(`rust/Cargo.toml`) — present in dev/test builds, absent from the shipped wheel:

| Build | Command | `count_bam` present? |
|-------|---------|----------------------|
| dev / local | `maturin develop` (default feature) | ✅ yes |
| `cargo test` | default feature | ✅ yes |
| CI test (`test.yml`) | `maturin build --release` (default) | ✅ yes |
| **shipped wheel** (Dockerfile, `release.yml`) | `maturin build --release --no-default-features` | ❌ no |

So the parity tests run against a build that has it; production ships without it. If
you build `--no-default-features` locally, `gbcms._rs.count_bam` is missing and the
parity tests fail — that's expected.

**Maintenance contract — read before changing the counting core.** The binned↔legacy
parity tests (`count_both` in `tests/helpers.py`; `test_filters.py`,
`test_parity_large_deletion.py`, `test_multi_allelic.py`, …) assert the two paths
return identical `PARITY_FIELDS`. Therefore:

- **If you change read classification, filtering, fragment consensus, or fetch-window
  logic in the binned path** (`count_bin_shared` / `count_variant_from_cache`), you
  **must mirror the same change in `count_single_variant`**, or the parity tests fail.
  The duplication is deliberate — two independent implementations are what make the
  cross-check meaningful.
- **Exempt:** RNA-only / mFSD / ASJD / strandedness features live only in the binned
  path and are *not* in `PARITY_FIELDS`; do **not** add them to `count_single_variant`.
- **Parity holds only without `sibling_variants`** — pangenomic sibling disambiguation
  is binned-only and intentionally diverges from legacy. Never pass siblings to
  `count_both`. (See the `siblings-break-binned-legacy-parity` memory.)

**When to remove it.** `count_bam` + `count_single_variant` can be deleted once the
binned path is trusted enough that the cross-check is no longer needed ("remove after
parity sign-off"). Until then, keep it gated and in sync. If a change is genuinely
impractical to mirror, remove `count_bam` and its parity tests *with maintainer
sign-off* — never let the two paths silently diverge.

## Type Stub Synchronization

- `src/gbcms/_rs.pyi` is **authoritative** — edit here first
- `src/gbcms_rs.pyi` is a **mirror** — must stay synced
- Both include: `BaseCounts`, `PreparedVariant`, `with_ad()`, diagnostic fields, ASJD fields, nucleosomal fraction
