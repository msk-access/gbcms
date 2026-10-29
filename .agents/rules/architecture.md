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
│   ├── merge.py         # `gbcms merge`: per-BAM-type MAFs (duplex, simplex) -> one, simplex_duplex_* columns
│   ├── observations.py  # Per-molecule observation export (public API; companion Parquet)
│   ├── rescue_audit.py  # The gbcms_rescue audit contract, shared by the MNP rescue pass and merge
│   ├── qc_summary.py    # End-of-run QC summary: the flag tier table, read from the written output, one log block
│   ├── qc_export.py     # `gbcms summarize`: the QC summary as text, JSON, a cohort TSV or MultiQC content
│   ├── _rs.pyi          # Primary Rust type stubs (authoritative)
│   ├── core/
│   │   └── kernel.py    # Coordinates; the one VCF<->MAF rule (vcf2maf/maf2vcf) + type labels
│   ├── io/
│   │   ├── input.py     # VcfReader (skips uncountable ALTs), MafReader
│   │   ├── output.py    # VcfWriter, MafWriter (mFSD/RNA column gating)
│   │   ├── reference.py # Reference bases for MAF -> VCF anchors
│   │   ├── atomic.py    # Atomic output files (temp file + fsync + rename)
│   │   └── batch.py     # Polars batch I/O for post-processing (merge, mFSD report)
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
│       │   ├── variant_checks.rs # splice triage; check_snp/mnp/ins/del/complex CIGAR checks, windowed scan
│       │   ├── carrier.rs       # Exact-carrier rule: pure indels, delins, Del+SNV, structural MNP; neighbouring-change guard
│       │   ├── junction.rs      # Junction judge for pure deletions of 50 bp or more
│       │   ├── window.rs        # allele_kind, read/guard/discrimination windows, informative-read tests
│       │   ├── observed.rs      # OBSERVED_ALLELE / COEXISTING_ALLELE: the allele the reads carry
│       │   ├── utils.rs         # ClassifyResult/ClassifyPhase, haplotype and masked-compare helpers
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
│       │   ├── gtf.rs           # GTF loader (plain or gzip/BGZF, variant-guided streaming)
│       │   └── gtf_line.rs      # One GTF line's exon fields (noodles-gtf's grammar)
│       ├── shared/
│       │   ├── fragment.rs      # FragmentEvidence, QNAME/UMI hashing
│       │   ├── stats.rs         # Fisher's exact test, strand bias
│       │   ├── bam_utils.rs     # median_qual, find_read_pos
│       │   ├── contig.rs        # Contig-name normalization across BAM/GTF/VCF (chrM vs MT)
│       │   ├── filters.rs       # ReadFilter, FilterCounts
│       │   └── baq.rs           # Heuristic BAQ (Li 2011)
│       └── normalize/           # Variant preparation (PreparedVariant)
│           ├── mod.rs, engine.rs, types.rs  # entry, preparation engine, PreparedVariant
│           ├── left_align.rs    # Left-alignment
│           ├── decomp.rs        # Homopolymer-twin detection (a delins miscollapsed in a run)
│           ├── fasta.rs         # Reference windows
│           └── repeat.rs        # Tandem-repeat detection and adaptive padding
├── nextflow/            # Nextflow pipeline wrapper
├── tests/               # pytest suite: contract tests, census.py, helpers.py (`ls tests`)
├── docs/                # MkDocs documentation
└── .agents/             # Agent rules, skills, memory, learnings; archive/ holds retired plans
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
8. **Dual alignment backends**: PairHMM (the CLI default, `--alignment-backend pairhmm`; `hmm` is an alias; WFA2 triage runs in front of it) or Smith-Waterman (`sw`, opt-in). They are the residual fallback behind the exact-carrier and junction rules (counting-engine skill).
9. **Genomic binning**: ~10kb bins, one `bam.fetch()` per bin, max 200 variants/bin.
10. **COITree for annotation**: Platform-portable metadata access via `Borrow` trait (nosimd vs NEON/AVX backends). The index is built per run from the GTF (plain or gzip/BGZF; ~2 s for a whole Ensembl GTF, ~6.5 s for GENCODE comprehensive; each line's feature and chromosome read before anything else is parsed); there is no disk cache (the 6.6.0 build-gtf-cache command and its cache-dir option are gone).
11. **Diagnostic flags**: `gbcms_diagnostic` and `gbcms_rescue` are strongly-typed Rust fields, not dynamic attributes.
12. **Statistics stay in Rust** (KS, LLR, Fisher in `counting/mfsd.rs` and `shared/stats.rs`); there is no scipy dependency. The exact KS p-value is a self-contained Rust DP, checked against SciPy reference values baked into the tests as literals.

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
repeat_span + 2)`, the anchor's included, and its guard window
(`window::guard_window`): for a variant the exact-carrier rule judges, a whole read's
reach (300 bases) each side of its REF span, the fragments the guard's placed-apart
check reads.

**How it is checked.**
- `build_genomic_bins` property test (Rust): every variant in exactly one bin whose
  fetch holds its read window and its guard window, over random clusters × window
  {1, 7, 100, 10k} × cap {1, 2, 3, 200}.
- `tests/test_binning_invariance.py`: the geometries above on synthetic DNA pairs
  (plain; BAQ + UMI + mFSD; siblings, a decomposed twin, overlapping mates), clip
  carriers aligned past a long anchor (the one read shape in which an under-fetched
  anchor changes a count: every other counted read overlaps the event's first base),
  an RNA locus with a GTF and antisense reads, and the real test BAM; observation rows
  too. Each geometry is asserted to split the bins (the engine's bin log), and the
  decoys to change them, so no check is vacuous. Not exercised: ASJD (the RNA fixture
  has no junction-borne ALT) and the `--mfsd-parquet` size arrays (no getters).
- `tests/helpers.py` `count_checked` / `count_bam_checked`: every counting test also
  counts with one variant per bin and compares every field. Most such calls hold one
  variant, so they check the per-variant fetch window only; multi-variant binning,
  threads and bin starts are covered by `test_binning_invariance.py`.
- `bin_window` / `bin_max_variants` are test arguments of `count_bam_binned` and
  `count_bam_binned_observations`; the pipeline passes neither.

**Classification** is checked against the read census (`tests/census.py`,
`tests/test_read_census.py`): each read judged by its own bases across the event's
tract, with the decided rules (REF needs one base past the first difference, ALT only
that base). Open decisions are strict xfails there. A second engine sharing the
classifier cannot see classification bugs, which is why the per-variant counting
oracle was retired (#170).

## Type Stub Synchronization

- `src/gbcms/_rs.pyi` is the **single** stub for the `gbcms._rs` extension module —
  edit it whenever the `#[pyo3(get)]` fields or `#[pyo3(signature)]` params change
- It includes: `BaseCounts`, `PreparedVariant`, diagnostic fields, ASJD fields, nucleosomal fractions
- There is no second stub (AGENTS.md invariant 5)
