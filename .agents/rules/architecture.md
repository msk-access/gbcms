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
│   ├── convert.py       # Standalone VCF <-> MAF conversion
│   ├── _rs.pyi          # Primary Rust type stubs (authoritative)
│   ├── core/
│   │   └── kernel.py    # Coordinates; the one VCF<->MAF rule (vcf2maf/maf2vcf) + type labels
│   ├── io/
│   │   ├── input.py     # VcfReader (skips uncountable ALTs), MafReader
│   │   ├── output.py    # VcfWriter, MafWriter (mFSD/RNA column gating)
│   │   └── reference.py # Reference bases for MAF -> VCF anchors
│   ├── models/
│   │   └── core.py      # GbcmsConfig, OutputConfig, AlignmentConfig (Pydantic)
│   ├── report/
│   │   └── mfsd_report.py  # HTML mFSD report generator
│   └── utils/
│       └── logging.py   # Structured logging setup
├── rust/                # Rust crate (gbcms_rs)
│   └── src/
│       ├── lib.rs       # PyO3 module entry
│       ├── types.rs     # Variant, BaseCounts PyO3 bindings
│       ├── counting/
│       │   ├── engine.rs        # Main counting loop, genomic binning, Rayon par_iter()
│       │   ├── variant_checks.rs # check_snp/mnp/ins/del/complex, windowed scan
│       │   ├── carrier.rs       # Exact-carrier rule for complex variants (delins, Del+SNV, structural MNP)
│       │   ├── observed.rs      # OBSERVED_ALLELE / COEXISTING_ALLELE: the allele the reads carry
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
│       │   ├── gtf.rs           # GTF parser (variant-guided streaming)
│       │   └── cache.rs         # GTF disk cache (GtfIndexBundle, bincode) — M5a
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
3. **mFSD is opt-in** (`--mfsd`) — gated at *both* layers (output-aware engine, invariant #3). Writers gate 41 MAF cols / 13 VCF INFO fields behind `self.mfsd` (absent when off, not NA-filled); the **binned engine also gates the compute** — `mfsd` is plumbed `OutputConfig.mfsd → count_bam_binned → count_variant_from_cache`, and when off the per-fragment size arrays, the `compute_mfsd_stats` stats, and the post-counting mFSD BH-FDR pass are all skipped (no compute-then-discard, no held `ref_sizes`/`alt_sizes`).
4. **Rust-native Parquet** (`--mfsd-parquet`): `write_fsd_parquet()` via `arrow`/`parquet` crates with ZSTD(1). No `pyarrow`.
5. **4-layer CLI validation**: Parse-time (Typer) → Pre-model (cli.py) → Model-time (Pydantic) → No silent skips.
6. **Fragment counting always on**: Quality-weighted consensus; discards counted in DPF not RDF/ADF.
7. **Windowed indel detection**: ±5bp scan expanding to `max(5, repeat_span + 2)`.
8. **Dual alignment backends**: SW (default) or PairHMM (`--alignment-backend hmm`).
9. **Genomic binning**: ~10kb bins, one `bam.fetch()` per bin, max 200 variants/bin.
10. **COITree for annotation**: Platform-portable metadata access via `Borrow` trait (nosimd vs NEON/AVX backends). The tree *layout* is arch-specific, so the **GTF disk cache (`--gtf-cache-dir`, M5a) never serializes the trees** — it persists only the parsed intermediate (`GtfIndexBundle`: exon records, splice sites, introns, chrom map; `bincode`, version-tagged) and rebuilds the trees via `build_exon_trees` on load. Caching is best-effort (missing/corrupt/stale/unwritable → log + plain parse); keyed on GTF identity + variant chroms. Concurrent cohorts must pre-warm via `gbcms build-gtf-cache` (else the first wave all cold-miss); the per-sample `count_bam_binned` runs then load in ~0.05s instead of re-parsing (~9s).
11. **Diagnostic flags**: `gbcms_diagnostic` and `gbcms_rescue` are strongly-typed Rust fields, not dynamic attributes.

## Binning invariance

**The rule.** The engine groups variants into bins, fetches each bin's reads once, and
counts every member from that cache. Which variants share a bin, and how far it
reaches, must never change a count: every field of every row is the same under any
bin window or per-bin cap, one variant per bin (window 1, cap 1, the per-variant
fetch through the production loop), one variant per call, shuffled input on several
threads, or extra rows that move where bins start. The Benjamini-Hochberg q-values
(`mfsd_qval_alt_ref`, `asjd_qval`) are computed over the rows of a call, so they are
excepted when the set of rows changes. It follows that every bin's fetch must hold
each member's read window `[pos − pad, pos + ref_len + pad)`, `pad = max(5,
repeat_span + 2)`, the anchor's included.

**How it is checked.**
- `build_genomic_bins` property test (Rust): every variant in exactly one bin whose
  fetch holds its window, over random clusters × window {1, 7, 100, 10k} × cap
  {1, 2, 3, 200}.
- `tests/test_binning_invariance.py`: the geometries above on synthetic DNA pairs
  (plain; BAQ + UMI + mFSD), clip carriers aligned past a long anchor (the one read
  shape in which an under-fetched anchor changes a count: every other counted read
  overlaps the event's first base), an RNA locus with a GTF, and the real test BAM;
  observation rows too.
- `tests/helpers.py` `count_checked` / `count_bam_checked`: every counting test runs
  production bins and one variant per bin and compares every field.
- `bin_window` / `bin_max_variants` are test arguments of `count_bam_binned` and
  `count_bam_binned_observations`; the pipeline passes neither.

**Classification** is checked against the read census (`tests/census.py`,
`tests/test_read_census.py`): each read judged by its own bases across the event's
tract, with the decided rules (REF needs one base past the first difference, ALT only
that base). Open decisions are strict xfails there. A second engine sharing the
classifier cannot see classification bugs, which is why the per-variant `count_bam`
oracle was retired (#170).

## Type Stub Synchronization

- `src/gbcms/_rs.pyi` is the **single** stub for the `gbcms._rs` extension module —
  edit it whenever the `#[pyo3(get)]` fields or `#[pyo3(signature)]` params change
- It includes: `BaseCounts`, `PreparedVariant`, diagnostic fields, ASJD fields, nucleosomal fractions
- (Historical: a top-level `src/gbcms_rs.pyi` used to mirror this; it stubbed a module
  nothing imports and was removed in LO-1. Don't reintroduce it.)
