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
2. **0-based internal coordinates**: 1-based in VCF/MAF externally; converted at boundary.
3. **mFSD is opt-in** (`--mfsd`): Writers gate 34 MAF cols and 7 VCF INFO fields behind `self.mfsd`. When off, columns are **absent** (not NA-filled).
4. **Rust-native Parquet** (`--mfsd-parquet`): `write_fsd_parquet()` via `arrow`/`parquet` crates with ZSTD(1). No `pyarrow`.
5. **4-layer CLI validation**: Parse-time (Typer) → Pre-model (cli.py) → Model-time (Pydantic) → No silent skips.
6. **Fragment counting always on**: Quality-weighted consensus; discards counted in DPF not RDF/ADF.
7. **Windowed indel detection**: ±5bp scan expanding to `max(5, repeat_span + 2)`.
8. **Dual alignment backends**: SW (default) or PairHMM (`--alignment-backend hmm`).
9. **Genomic binning**: ~10kb bins, one `bam.fetch()` per bin, max 200 variants/bin.
10. **COITree for annotation**: Platform-portable metadata access via `Borrow` trait (nosimd vs NEON/AVX backends).
11. **Diagnostic flags**: `gbcms_diagnostic` and `gbcms_rescue` are strongly-typed Rust fields, not dynamic attributes.

## Type Stub Synchronization

- `src/gbcms/_rs.pyi` is **authoritative** — edit here first
- `src/gbcms_rs.pyi` is a **mirror** — must stay synced
- Both include: `BaseCounts`, `PreparedVariant`, `with_ad()`, diagnostic fields, ASJD fields, nucleosomal fraction
