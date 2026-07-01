# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### 🐛 Fixed

- **Contig-naming mismatches no longer silently produce zero counts.** When the BAM's
  contigs were named differently from the variants (UCSC `chr1` vs Ensembl/b37 `1`), the
  binning step found no matching contig, built 0 bins, and returned zero counts for every
  variant while the run still exited 0 — a systematic failure masked as success. `resolve_tid`
  now reconciles the naming via `normalize_contig` (so `chr1`↔`1`, `chrM`↔`MT` match), and a
  genuinely-absent chromosome now logs a loud `WARN` (once per chromosome) instead of being
  skipped silently. (MSK ACCESS runs were unaffected — BAMs and MAFs both use b37 naming.)
  Surfaced by the new end-to-end MAF test.

- **A run now exits non-zero when a sample fails (HI-1).** `Pipeline.run()` catches
  per-sample errors and returns normally, so `gbcms dna`/`rna` previously exited `0` even
  when every BAM failed (e.g. a Rust panic surfaced as `PyErr`) — masking systematic
  failure as success under Nextflow. The CLI now exits **code 1** when any sample actually
  fails (all-failed or partial). An **empty variant set is not a failure**: a sample that
  legitimately has no variants called still exits `0`, so per-sample workflows don't fail on
  it. Successful runs are unchanged (exit 0).

- **Rejected variants keep their reason in the output.** When *every* variant is rejected
  during preparation (e.g. a contig mismatch → `FAIL_FETCH_FAILED`), the run no longer
  short-circuits with no output — it now writes each variant with its `FAIL_*` reason in the
  `gbcms_status` column (and zero counts), so the reasons are in the output file, not only
  the log. (Partial rejections already did this; this closes the all-rejected gap.)

### 🔄 Changed

- **Per-transcript count columns now use `|` between transcripts, not `;` (ME-2).**
  `transcript_read_counts` / `transcript_fragment_counts` (MAF) and `TXRC` / `TXFC` (VCF)
  separate transcripts with `|` — e.g. `ENST…:AD,RD,DP|ENST…:AD,RD,DP`. `;` is the VCF INFO
  field separator, so the VCF path already converted to `|` while MAF emitted raw `;`; the
  engine now emits `|` directly so all three (MAF, VCF, the documented header) agree.
  **Consumers that split the MAF column on `;` must switch to `|`.**

- **`--threads` is now a hard, validated thread budget.** It is the *total* worker
  budget for one gbcms process (multi-sample parallelism is Nextflow's job — N
  concurrent processes, each pinned to `task.cpus`). `--threads 0` is now rejected
  loudly instead of silently becoming rayon's all-cores default (which would
  oversubscribe a small SLURM allocation); all rayon pools are sized through
  `resolve_thread_budget`, and the resolved budget is logged. Counts are unaffected.

### 📦 Packaging

- **The shipped wheel no longer exports the legacy `count_bam` parity oracle.** The
  per-variant `count_bam` (the binned↔legacy parity oracle) is now behind a default
  `legacy-parity` Cargo feature; release builds use `--no-default-features` to omit it.
  Production always used `count_bam_binned`, so this only trims a test-only symbol from
  the wheel. Dev/test builds keep it (default on) so the parity suite still runs.

### ⚡ Performance

- **mFSD is now computed only when requested (`--mfsd`/`--mfsd-parquet`/`--mfsd-report`).**
  Previously the engine always built the per-fragment size arrays and ran the full
  KS/LLR/delta statistics for every variant, then discarded them unless an mFSD output
  was selected. The binned production path now gates this work on the `mfsd` flag
  (plumbed from `OutputConfig.mfsd`), so the default `--dna` run no longer holds the
  per-variant `ref_sizes`/`alt_sizes` arrays across the whole sample — the dominant mFSD
  memory cost, and the one that multiplies under Nextflow fan-out (N concurrent
  processes). The size-array statistics are also factored into a single shared
  `compute_mfsd_stats` helper used by both the binned and legacy paths. **Counts are
  unaffected:** validated on 3,040 real cfDNA variants — all 246 non-mFSD columns are
  byte-identical with mFSD on vs off, and the 41 mFSD columns appear only when enabled.

- **GTF annotation can be cached on disk (`--gtf-cache-dir`).** Parsing a full Ensembl
  GTF takes ~9s, and under a Nextflow cohort that cost was paid once *per sample*. When
  `--gtf-cache-dir` is set, the parsed intermediate (exon records, splice sites, introns,
  chrom map) is serialized once and reused on later runs over the same GTF + variant set,
  dropping the annotation load from ~9s to **~0.05s** (validated on the full GRCh38.111
  GTF; cold-vs-warm output byte-identical across 47 variants and all annotation columns).
  The architecture-specific interval trees are *not* cached — they are rebuilt from the
  cached records on load, so a cache file is portable across x86/ARM. Caching is
  best-effort: a missing/corrupt/stale-version/unwritable cache logs and falls back to a
  normal parse, never affecting correctness.

### ✨ Added

- **VCF now emits the mFSD sub-/mono-nucleosomal fields (ME-1).** `MFSD_SUB_NUC_REF_FRAC`,
  `MFSD_SUB_NUC_ALT_FRAC`, `MFSD_SUB_NUC_ENRICHMENT`, `MFSD_MONO_NUC_REF_FRAC`,
  `MFSD_MONO_NUC_ALT_FRAC` were computed under `--mfsd` and written to MAF but silently
  dropped from VCF — they're now in the VCF INFO too (VCF↔MAF parity), taking the VCF mFSD
  surface to 13 INFO fields. Also corrected long-standing count drift in the docs/CLI help:
  `--mfsd` adds **41** MAF columns (not 34) and **13** VCF INFO fields (not 7).

- **`gbcms build-gtf-cache` — pre-warm the GTF index cache for a cohort.** A dedicated
  command (`--gtf --variants --gtf-cache-dir`, no BAM) that parses the GTF once and writes
  the cache, so a fan-out of per-sample `gbcms rna --gtf-cache-dir <same-dir>` jobs all
  start warm. This is what makes the cache pay off under concurrency: without a pre-warm
  step, samples launched together all cold-miss and each re-parses the GTF. Run it as a
  single Nextflow process upstream of the per-sample fan-out so the whole cohort parses
  the GTF exactly once. See `docs/nextflow/parameters.md`.

- **`--strandedness` for RNA mode (`reverse` | `forward` | `unstranded`).** The
  read→transcript-strand fold was previously hardcoded to dUTP/reverse
  (fr-firststrand, featureCounts `-s 2`); it is now selectable to support forward
  (fr-secondstrand, `-s 1`) and unstranded (`-s 0`) libraries. The fold drives both
  `--enforce-strandedness` filtering and ASJD strand-discordance detection;
  `unstranded` disables both (as do amplicon libraries). Default is `reverse`,
  preserving prior behavior and matching the FORTE pipeline default. Unknown values
  are rejected loudly at both the model and FFI layers. Validated on a real
  reverse-stranded RNA sample: `reverse + forward = unstranded` read counts exactly
  (100% strand partition) at ACTB and GAPDH.

### 🔄 Changed

- **ASJD junction evidence is now counted per fragment, not per mate.** A molecule
  whose R1 and R2 overlap on a short cDNA insert can have both mates carry the same
  splice junction (the genomic insert looks large only because it spans the intron) —
  on real reverse-stranded RNA this is 35.6% of fragment×junction incidences and is
  *not* a UMI/PCR duplicate. Counting both mates inflated the per-junction strand
  tally ~1.38× and could fire spurious `STRAND_DISCORDANT` at low depth. `detect_asjd`
  now dedups each fragment (by QNAME) to one vote per allele-total and per junction,
  matching the per-fragment treatment the rest of the engine already uses. Mates always
  fold to the same transcript strand (verified 0/319k disagreements), so the dedup is
  unambiguous. This shifts `asjd_n_*_junc` / `asjd_n_*_total` and a small number of
  low-count `STRAND_DISCORDANT` flags (18/58,240 junctions genome-wide on the test
  sample), always in the corrective direction.
- **Empty-allele variants are now rejected at prep time.** A structurally empty REF
  or ALT (`""`) is malformed input — the internal representation is VCF-style
  (anchor-based), and MAF dash alleles arrive as the literal `-`, never empty. Such
  variants previously fell through to counting and silently produced zero counts;
  they are now rejected during `prepare_variants` with a `FAIL_EMPTY_ALLELE` status
  and a warning, so they are surfaced in the output rather than quietly dropped to
  all-zero. Legitimate MAF `-` alleles are unaffected.
- **Supplementary/secondary alignments no longer count toward read-level depth.**
  They share a QNAME with their primary, so counting them at read level
  double-counts DP/RD/AD at the anchor. They are now skipped unconditionally in
  both the binned and legacy counting paths, independent of `--no-filter-secondary`
  / `--no-filter-supplementary` (those flags now only govern whether such records
  enter the read cache, not whether they are counted). Affects only the
  non-default flag combination; default behavior is unchanged.

### 🐛 Fixed

- **Read-level supplementary/secondary double-count** under `--no-filter-secondary`
  / `--no-filter-supplementary` (see Changed above).
- **`check_complex` trailing-insertion handling clarified and guarded.** The Phase-1
  reconstruction deliberately includes an insertion at the exclusive REF end (a
  trailing insertion that belongs to the ALT, e.g. REF=`AB`, ALT=`ABC`); documented
  the rationale and added regression tests so the boundary condition isn't "fixed"
  into a misclassification.

## [5.3.0] - 2026-05-16

### ✨ Added

- **CRAM Support**: Full CRAM file support across all commands (`dna` and `rna`). The `--bam` argument now transparently accepts CRAM files.
- **CRAM Index Auto-Discovery**: The Nextflow pipeline automatically discovers `.crai` or `.cram.crai` indexes in the same directory as the CRAM files.
- **Reference Binding**: Rust engine's `IndexedReader` correctly initializes CRAM references using the provided `--fasta`.

### 🔄 Changed

- **Provenance Comment Lines**: Output MAF files now include `#gbcms vX.Y.Z` and `#command ...` provenance lines before the header row in both DNA and RNA modes.
- **VCF Header Provenance**: VCF headers now include `##source`, `##gbcms_command`, `##reference`, and `##contig` metadata in both DNA and RNA modes.
- **Strand Bias Sanitization**: Handled edge cases in Fisher's exact test where structural strand arrays contain all zeros, returning `.`, `1.0`, or `NA` instead of `NaN`/`Inf`.
- **Merge Engine Robustness**: `gbcms merge` uses Polars `comment_prefix` natively to skip provenance comment lines without failure.

### 📚 Documentation

- Updated `cli/rna.md` to note CRAM and provenance support for RNA mode.
- Updated `cli/merge.md` to explain MAF comment line compatibility.
- Updated `cli/index.md` and diagrams to explicitly mention `BAM/CRAM Files`.
- Expanded `reference/input-formats.md` with explicit CRAM and `.crai` index requirements.
- Added DNA vs RNA comparison snippet in `reference/output-formats.md` demonstrating provenance headers.

### 🧪 Tests

- Expanded testing suite (now 329 tests), providing integration coverage for RNA provenance headers and CRAM pipeline data flow.
- Added test validation for merge engine skipping MAF provenance comment lines safely.

## [5.2.0] - 2026-05-15

### ✨ Added

- **`gbcms merge` command**: New CLI command for merging per-BAM-type genotyped
  MAFs (e.g., duplex + simplex) into a single output with type-prefixed count
  columns and optional additive combined metrics.
- **Merge engine** (`src/gbcms/merge.py`): Polars-based lazy join engine with
  3-phase combined column computation:
  - Phase 1: Additive sums (12 read + fragment + strand count columns)
  - Phase 2: Derived totals and VAFs
  - Phase 3: Fisher's exact strand bias (read + fragment level) via Rust
- **Batch I/O module** (`src/gbcms/io/batch.py`): Centralized Polars-based
  `read_maf`, `scan_maf`, `read_parquet`, `write_maf` for batch operations.
- **`MergeConfig` Pydantic model**: Validated configuration with ≥2 input
  enforcement, file existence checks, and `add_combined`/`legacy_naming` options.
- **`fisher_exact_2x2` PyO3 wrapper**: Exposed the existing Rust Fisher's exact
  test as `gbcms._rs.fisher_exact_2x2()` for Python-side strand bias computation
  in the merge engine, ensuring numerical parity with the primary counting engine.
- **Nextflow `MERGE_COUNTS` process**: New module at
  `nextflow/modules/local/gbcms/merge/main.nf` that calls `gbcms merge` with
  `--input type:path` arguments built from grouped BAM type channels.
- **Nextflow `bam_type` samplesheet column**: Optional column enabling automatic
  `--column-prefix` derivation and `groupTuple`-based merge orchestration.
- **Nextflow merge parameters**: `merge_counts`, `merge_add_combined`,
  `merge_legacy_naming` in `nextflow.config`.
- **Polars dependency**: `polars>=1.0.0` added as a core dependency.

### 🔧 Fixed

- **Fragment consensus: INDEL structural evidence priority** — When R1 and R2
  of a duplex fragment disagree on an insertion or deletion, the read with direct
  CIGAR evidence (I/D op) now wins the fragment consensus unconditionally, instead
  of comparing anchor base qualities (which are identical for both reads and
  caused systematic discard). This recovers ~2-5% of INDEL fragment-level evidence
  that was previously lost, critical for low-VAF cfDNA detection.
  - Insertions: `adf` increases by ~5% for conflict fragments (validated on RB1 INS)
  - Deletions: `adf` increases by ~2% for conflict fragments (validated on DNMT3A DEL)
  - Single-molecule events: `adf=0→1` recovery (validated on RUNX1 INS, ad=1)
  - SNP behavior: unchanged (no structural flags on SNP classifications)
  - Phase 3 alignment returns: unchanged (non-structural, quality comparison retained)

### 🧪 Tests

- **[NEW]** `tests/test_merge.py` — 24 tests covering variant key join,
  multi-type merge, combined columns (read + fragment + strand), Fisher strand
  bias (biased + balanced), column order validation, asymmetric row counts,
  null-fill, legacy naming, and CLI input format validation.
- **[NEW]** `tests/test_batch_io.py` — 8 tests covering `read_maf`, `scan_maf`,
  `write_maf`, `read_parquet`, and error handling for missing files.
- **[NEW]** `tests/test_indel_fragment_consensus.py` — 11 tests covering INS/DEL
  conflict recovery (structural ALT priority), wrong-length INDEL Phase 3 dispatch,
  agreement paths, singleton reads, SNP regression (tie behavior unchanged), and
  REF agreement at INDEL sites. All 4 counting invariants asserted per test.
- **[NEW]** `rust/src/shared/fragment.rs` `#[cfg(test)]` — 11 Rust unit tests for
  `FragmentEvidence::resolve()` structural priority and `observe()` sticky flag logic.
- **298 Python tests** (up from 255): 32 merge, 8 batch, 11 INDEL consensus tests.
- **161 Rust tests** (up from 150): 11 new fragment consensus unit tests, 0 Clippy warnings.

### 📚 Documentation

- **[NEW]** `docs/cli/merge.md` — CLI reference for `gbcms merge` covering usage,
  combined columns schema, Nextflow integration, and Fisher strand bias.
- **Updated** `docs/cli/index.md` — Added `merge` command to commands table.
- **Updated** `docs/reference/output-formats.md` — Added Merged MAF Output section
  with type-prefixed columns and 20 combined metrics.
- **Updated** `docs/reference/architecture.md` — Merge engine in system overview
  diagram, `batch.py` and `merge.py` in module tree, `MergeConfig` in config diagram.
- **Updated** `src/gbcms/_rs.pyi` — Added `fisher_exact_2x2` type stub.
- **Updated** `src/gbcms/__init__.py` — Exported `merge_mafs` and `MergeConfig`.

## [5.1.0] - 2026-05-11

### ⚠️ Breaking Changes

- **MAF column order reordered** (v5.1 schema):
  `any_alt`, `partial_alt`, `n_count` moved from end to immediately after
  `alt_count` for discoverability. Read strand counts (`ref_count_forward`, etc.)
  now precede derived strand bias statistics. Read and fragment metric layers
  fully separated (no interleaving). Downstream MAF parsers using positional
  indexing must be updated.
- **VCF FORMAT fields restructured** (VCF 4.2 spec compliance):
  - `DP` is now a single integer (total depth), was `ref,alt` pair.
  - `AD` is now `Number=R` with `ref,alt` totals (VCF spec), was `fwd,rev`.
  - `RD` and `RDF` removed. Replaced by `ADF` (forward strand `ref_fwd,alt_fwd`)
    and `ADR` (reverse strand `ref_rev,alt_rev`) following bcftools convention.
  - New `FAD`, `FADF`, `FADR` for fragment-level strand-by-allele counts.
  - `FAF` renamed from position after `VAF` to after fragment group.
  - Downstream VCF parsers expecting old `GT:DP:RD:AD:RDF:ADF:VAF:FAF:...` must
    be updated to `GT:DP:AD:ADF:ADR:VAF:FAD:FADF:FADR:FAF:...`.
- **VCF INFO field order changed**: `AAD`, `PAD`, `NAD` now appear immediately
  after `GR` (before strand bias fields), matching the diagnostic proximity
  principle.

### 🔧 Fixed

- **Wrong-length insertion Phase 3 fallback** (PAX5-class discordance fix):
  `check_insertion` now routes wrong-length insertions at the strict anchor
  position to Phase 3 (SW/PairHMM) for haplotype arbitration, mirroring
  `check_deletion`'s existing behavior. Previously, a read with `I(1)` at the
  anchor for an expected `I(2)` was silently classified as REF with no
  diagnostic signal. Now classified as `partial_alt` via `has_nearby_evidence`,
  triggering the `PARTIAL_DOMINANT` diagnostic flag.
- **Wrong-length insertion windowed scan tracking**: Wrong-length insertions
  within the ±window range now set `has_nearby_length_match`, routing to
  Phase 3 fallback. Previously silently ignored.
- **Insertion `!found_ref_coverage` haplotype fallback**: Added Phase 3
  fallback for insertion reads where the CIGAR walk found no evidence but
  the read spans the anchor (e.g., unusual CIGAR geometry, soft-clip at
  anchor). Mirrors `check_deletion`'s existing `!found_ref_coverage` path.
- **Wrong-length deletion windowed scan tracking**: Wrong-length deletions
  in the windowed scan now set `has_nearby_length_match` for Phase 3
  fallback. Covers two previously silent-drop cases:
    - Small deletions (≥5bp, <50bp): always flagged (1-4bp excluded as
      homopolymer noise)
    - Large deletions (≥50bp) with low reciprocal overlap (<50%):
      previously silently dropped, now flagged for Phase 3 arbitration

## [5.0.0] - 2026-05-11

### ⚠️ Breaking Changes

- **New RNA output columns**: When `--gtf` is provided, 17 additional columns are
  appended to RNA output (1 exon boundary, 2 per-transcript, 14 ASJD). Downstream
  parsers expecting a fixed column count must be updated.
- **Fragment counts in amplicon mode**: When `--library-type amplicon` is used,
  fragment counts (`dpf`, `rdf`, `adf`) will approximate read counts (`dp`, `rd`,
  `ad`) because R1/R2 fragment consensus is bypassed. This is expected behavior,
  not a bug.

### ✨ Added

- **GTF-based transcript annotation** (`--gtf`): RNA mode can now load an
  Ensembl/GENCODE GTF file to enable:
    - **Exon boundary distance** (`exon_boundary_dist`): Signed distance to the
      nearest exon boundary for splice-proximal variant filtering.
    - **Per-transcript counting** (`transcript_read_counts`,
      `transcript_fragment_counts`): ALT counts stratified by overlapping
      transcript, resolving ambiguity at multi-transcript loci.
    - **Aberrant Splice Junction Detection (ASJD)**: 14 columns comparing
      observed read splice junctions against annotated transcript splice sites.
      Flags novel vs annotated junctions with REF/ALT stratification.
- **`--library-type` CLI flag**: RNA library preparation method selector.
  `capture` (default, IDT xGen-style) or `amplicon`. Available on `gbcms rna` only.
- **Amplicon mode — fragment consensus bypass**: When `--library-type amplicon`
  is set, the Rust counting engine XORs `mol_hash` with a read-specific tag
  (`0x1` for R1, `0x2` for R2), treating each read as an independent molecule.
  This prevents incorrect R1/R2 merging in amplicon libraries where both reads
  are PCR duplicates of the same template.
- **Amplicon auto-strandedness override**: Both the CLI (`cli.py`) and model
  validator (`GbcmsRnaConfig.model_validator`) auto-disable `enforce_strandedness`
  when `library_type="amplicon"`, with a warning. This ensures safe defaults for
  both CLI and API users.
- **Per-read strand tracking in junction accumulator**: RNA reads now carry
  per-read strand orientation through the splice junction accumulator, enabling
  accurate sense/antisense counting at splice-spanning positions.
- **COITree annotation index** (`rust/src/annotation/`): New Rust module
  providing O(log n + k) exon overlap queries via cache-oblivious interval
  trees. Built once at startup, shared immutably across Rayon threads via
  `Arc<AnnotationIndex>`.
- **Splice mask construction**: Per-transcript `HashSet<(chrom, start, end)>`
  for O(1) splice site lookup during ASJD computation.
- **BH multiple-testing correction** (`shared/stats.rs`): Benjamini-Hochberg
  FDR correction added for strand bias p-values across variants.
- **Nextflow `library_type` parameter**: Registered in `nextflow.config` and
  threaded through `rna/main.nf` for amplicon mode support in pipeline runs.

### 🔧 Fixed

- **`lib.rs` dead-code comment**: Updated annotation type comment to reflect
  that annotation types are now fully wired (previously noted as unused).

### 🧪 Tests

- **255 Python tests** (up from 238): 17 new tests covering:
    - `test_config_isolation.py`: 11 tests for `library_type` field/validator,
      GTF field, DNA isolation, amplicon auto-strandedness, and `rescue_mnp_threshold`
      range validation (default, shared, >1.0 rejection, <0.0 rejection).
    - `test_cli_dna_rna.py`: 4 tests for `--library-type` and `--gtf` option
      isolation between DNA and RNA commands.
    - `test_diagnostic_flags.py`: 2 tests for `MNP_RESCUE_ELIGIBLE` threshold-based
      eligibility (conservative 0.50 mode).
- **143 Rust tests** (up from 119): 24 new tests covering:
    - GTF parsing and chromosome normalization
    - COITree overlap queries
    - Splice distance computation
    - BH FDR correction
    - Fisher's exact test edge cases
- **0 Clippy warnings** (strict `-D warnings` mode).

### 📚 Documentation

- **[NEW]** `docs/reference/rna-annotation.md` — GTF requirements, annotation
  index architecture, exon boundary distance, per-transcript counting, and
  ASJD detection reference.
- **Updated** `docs/cli/rna.md` — `--gtf` and `--library-type` options, updated
  pipeline diagram with annotation layer, amplicon example tab, amplicon override
  warning, updated DNA vs RNA comparison table.
- **Updated** `docs/reference/output-formats.md` — GTF-aware MAF columns (17),
  amplicon mode behavioral note.
- **Updated** `docs/reference/architecture.md` — AnnotationIndex in system
  overview diagram, `annotation/` module in tree, config diagram with
  `library_type` and `gtf` fields, `stats.rs` BH correction note.
- **Updated** `docs/development/developer-guide.md` — `annotation/` module in
  project structure diagram.
- **Updated** `mkdocs.yml` — "RNA Annotation" added to navigation.

## [4.2.0] - 2026-05-10

### ⚠️ Breaking Changes

- **`validation_status` → `gbcms_status`**: MAF column and Python API parameter
  renamed. `gbcms_status` now uses semicolon-separated multi-value format
  (e.g., `PASS;WARN_REF_CORRECTED`, `PASS;MULTI_ALLELIC`). The first token
  is always `PASS` or `FAIL_*`.
- **VCF `VS` → `GS`/`GD`/`GR`**: VCF INFO key `VS` replaced by `GS`
  (status), `GD` (diagnostic), and `GR` (rescue). Downstream VCF parsers
  must update field references.
- **Status value format**: Old underscore-joined statuses
  (`PASS_WARN_HOMOPOLYMER_DECOMP`, `PASS_MULTI_ALLELIC`) replaced with
  semicolon-separated (`PASS;WARN_HOMOPOLYMER_DECOMP`, `PASS;MULTI_ALLELIC`).

### ✨ Added

- **`gbcms_diagnostic` column (MAF) + `GD` INFO key (VCF)**: Post-counting
  diagnostic flags computed automatically. Flags include:
    - `ZERO_ALT`: No confirmed ALT reads despite successful counting.
    - `PARTIAL_DOMINANT`: More structural/partial evidence than confirmed ALT.
    - `MNP_DISC_RATIO(n/m)`: MNP discriminating position ratio (always emitted for MNPs).
    - `MNP_RESCUE_ELIGIBLE`: MNP qualifies for rescue (disc/len ≤ `--rescue-mnp-threshold`).
    - `HIGH_N_FRACTION(f)`: N-base fraction exceeding 5% at discriminating positions.
- **`--rescue-mnp` CLI flag**: Enables MNP rescue pass for multi-base
  substitutions. When `ad=0` and `MNP_RESCUE_ELIGIBLE` is flagged, decomposes the
  MNP into individual SNP positions and re-counts via `count_bam_binned`.
  Available in both `gbcms dna` and `gbcms rna` modes.
- **`--rescue-mnp-threshold` CLI flag**: Maximum disc/len ratio for MNP rescue
  eligibility (0.0–1.0, default: 1.0). At 1.0, all MNPs are eligible (C++ gbcms
  compatible). Set to 0.5 for conservative sparse-only mode.
- **`gbcms_rescue` column (MAF) + `GR` INFO key (VCF)**: Conditional —
  only present when `--rescue-mnp` is enabled. Contains structured audit trail:
  `method=decomposed;original_alt=0;positions=chr:pos(R>A):count,...`.
  Failed rescues include `outcome=no_signal`.
- **`has_nearby_evidence` (Rust)**: New `ClassifyResult` field propagating
  structural evidence from variant checkers (INS/DEL/Complex) and alignment
  backends (SW/PairHMM). Enables `partial_alt` counting for INDELs.
- **`partial_alt` now populated for INDELs**: Previously always 0 for
  SNP/INDEL. Now fired when the counting engine detects nearby structural
  evidence (right-length INDEL, non-zero ALT alignment score).
- **Diagnostic flag summary logging**: `info`-level log of diagnostic flag
  distribution per sample (e.g., `ZERO_ALT=12, PARTIAL_DOMINANT=3`).

### 🔧 Fixed

- **`partial_alt` description**: Documentation corrected to reflect that
  `partial_alt` is now populated for all variant types, not just MNP/Complex.

### 📝 Notes

- **Invariant 1 breakage (rescue only)**: After rescue, `any_alt = ad + partial_alt`
  no longer holds for rescued variants. `ad` is updated with the best decomposed
  SNP count while `any_alt` and `partial_alt` retain original MNP-level values
  as forensic evidence.
- **Rescue strategy**: Python-side post-processing using decomposed SNP
  re-counting. Coordinate shift strategy is reserved for a future release.

### 🧪 Tests

- **[NEW]** `tests/test_diagnostic_flags.py` — 17 tests covering all 4
  diagnostic flags, multi-flag combinations, FAIL exclusion, boundary
  conditions, parametric formatting, and `rescue_mnp_threshold`-based
  eligibility gating (permissive 1.0 and conservative 0.50 modes).
- **[NEW]** `tests/test_rescue_mnp.py` — 13 tests covering config defaults,
  conditional column/INFO presence, candidate identification, guard rails
  (skip non-MNP, skip ad>0, skip FAIL), audit trail format, no-signal cases,
  column count with rescue, and invariant breakage verification.
- Updated `test_normalization.py` (12 refs), `test_pipeline_v2.py` (4 refs),
  `test_phase2_output.py` (column count 24→26), `test_multi_allelic.py` for
  `gbcms_status` rename.
- All 238 Python tests pass; all 119 Rust tests pass.

### 📚 Documentation

- **Updated** `docs/reference/output-formats.md` — `gbcms_status`,
  `gbcms_diagnostic`, `gbcms_rescue` columns; `partial_alt` description
  corrected; prefix behavior updated.
- **Updated** `docs/reference/architecture.md` — MNP rescue architecture,
  data flow diagram, Python-vs-Rust design rationale, invariant impact.
- **Updated** `docs/development/developer-guide.md` — Rescue debugging,
  extension guidelines, test fixture requirements.
- **Updated** `docs/reference/variant-normalization.md` — status field name
  and multi-value format.
- **Updated** `docs/reference/allele-classification.md` — multi-allelic
  status format.
- **Updated** `docs/cli/normalize.md` — `gbcms_status` column name.
- **Updated** `docs/resources/troubleshooting.md` — status field references.


## [4.1.0] - 2026-05-04

### ⚠️ Breaking Changes

- **`gbcms run` command removed**: The deprecated `gbcms run` alias (introduced
  in v4.0.0 as a transitional shim) has been removed. Use `gbcms dna` instead —
  all arguments are identical.
- **Nextflow config defaults aligned with CLI**: `filter_secondary`,
  `filter_supplementary`, `filter_qc_failed` changed from `false` to `true`;
  `enforce_strandedness` changed from `false` to `true`;
  `alignment_backend` changed from `'hmm'` to `'pairhmm'`. Pipelines relying
  on the old Nextflow defaults may see behavior changes.

### ✨ Added

- **Physical fragment sizing**: Rust-native fragment size calculation using
  aligned read positions instead of TLEN, improving accuracy for supplementary
  alignments and soft-clipped reads.
- **`--mfsd-report` flag**: Generates an interactive HTML report with per-variant
  fragment size distribution analysis, dual-axis histograms, and Fragment Origin
  Signal classification (TUMOR-LIKE / CH-LIKE / AMBIGUOUS / INSUFFICIENT).
  Implies `--mfsd` and `--mfsd-parquet`.
- **Variant navigator**: STRiDE-inspired sticky navigation bar for multi-variant
  mFSD reports — dropdown selector, prev/next buttons, Focus/Show All toggle,
  keyboard shortcuts (←/→). Automatically hidden for single-variant reports.
- **`--mfsd-report-min-alt`**: Minimum ALT fragment count to include a variant
  in the report (default: 3).
- **`--mfsd-report-max-variants`**: Maximum variants per report (default: 20).
- **Theme toggle**: Light/dark mode switching in HTML reports.
- **Print compliance**: Reports render audit-ready when printed (navigator hidden,
  all variants at full opacity, branded footer included).
- **RNA BAQ default**: `--apply-baq` now defaults to `True` for `gbcms rna`.
  RNA pipelines typically lack upstream BQSR, so BAQ penalizes bases near
  splice junctions and indels to reduce false-positive variant calls.
- **BAQ trace logging**: Per-read BAQ adjustments logged at `trace` level,
  showing indel and splice junction counts per read.
- **Nextflow `cache = 'lenient'`**: Ensures `-resume` works correctly on
  GPFS/Spectrum Scale filesystems where inode metadata changes during file
  pool migration.
- **Nextflow `manifest` block**: Pipeline metadata (name, version, author,
  homepage) for Nextflow Tower and nf-core registry compatibility.
- **Nextflow SLURM job naming**: `clusterOptions` adds descriptive job names
  (`nf-GBCMS_DNA_sampleid`) for `squeue` readability.
- **Nextflow `executor.queueSize`**: Caps concurrent SLURM submissions at 100.
- **nf-core institutional configs**: Auto-loads site-specific profiles (iris,
  jax, sanger, etc.) via `nf-core/configs`.
- **Extended trace fields**: IO diagnostics (`rchar`, `wchar`, `syscr`, `syscw`,
  `read_bytes`, `write_bytes`) added to execution trace.
### 🔧 Fixed

- **Dual-axis gridline artifact**: Fixed Plotly `yaxis2` overlay creating a
  duplicate x-axis line in mFSD histograms by standardizing `mirror: false`,
  `rangemode: 'tozero'`, and `showline` controls.
- **BAQ RefSkip early-exit bug**: `apply_heuristic_baq()` silently skipped
  reads containing only splice junctions (CIGAR `N`) but no indels. The
  early-exit gate now includes `Cigar::RefSkip`, ensuring splice-spanning
  reads receive the BAQ quality penalty.
- **Nextflow config defaults**: Aligned `nextflow.config` with CLI defaults —
  `filter_secondary`, `filter_supplementary`, `filter_qc_failed` corrected
  to `true`; `enforce_strandedness` corrected to `true`; `alignment_backend`
  corrected to `pairhmm`.
- **Nextflow RNA BAQ wiring**: Added `--no-baq` / `--apply-baq` argument
  to RNA module (`rna/main.nf`), which was previously missing entirely.
- **Nextflow RNA strandedness**: Fixed `strandedness_arg` logic — now passes
  `--no-strandedness` only when disabled (was incorrectly passing
  `--enforce-strandedness` as an additive flag).

### 📚 Documentation

- **[NEW]** `docs/reference/mfsd-report.md` — mFSD interactive report reference
  covering Fragment Origin Signal classification, interactive features, output
  columns, and print compliance.
- **[NEW]** `docs/reference/rna-splice-handling.md` — RNA splice-junction handling
  guide with dual-mechanism comparison (consensus intron snipping vs BAQ), GATK
  SplitNCigarReads comparison, defense-in-depth analysis (5 layers), and visual
  splice bleed examples.
- **Updated** BAQ documentation across 7 files: `read-filters.md`, `cli/dna.md`,
  `cli/rna.md`, `glossary.md`, `abbreviations.md`, `nextflow/parameters.md`,
  `architecture.md` — all now include BAQ_RADIUS/BAQ_PENALTY constants,
  mode-specific defaults, and guidance on when to enable BAQ for DNA.
- **Updated** `mkdocs.yml` — added "RNA Splice Handling" to navigation.
- **Updated** `docs/cli/dna.md` — added `--mfsd-report`, `--mfsd-report-min-alt`,
  `--mfsd-report-max-variants` to CLI reference.
- **Updated** `docs/nextflow/parameters.md` — added Nextflow params for mFSD report;
  BAQ default now shows mode-specific values.
- **Updated** `docs/development/release-guide.md` — version locations table corrected
  from 5 to 7 references (reflecting v4.0.0 Nextflow module split).
- **Updated** `nextflow/nextflow.config` — added `mfsd_report`, `mfsd_report_min_alt`,
  `mfsd_report_max_variants` pipeline parameters.
- **Updated** `nextflow/modules/local/gbcms/dna/main.nf` — wired `--mfsd-report` flags
  through the DNA module with HTML report output channel (`emit: mfsd_report`).

### 🧪 Tests

- **[NEW]** `tests/test_mfsd_report.py` — 13 unit tests covering report creation,
  navigator presence/absence, Plotly integration, theme toggle, branding, summary
  cards, min_alt/max_variants filtering, and error handling. Uses synthetic test
  fixtures with no patient identifiers.
- `mfsd_report.py` coverage: 0% → 91%.
- **[NEW]** `tests/test_config_isolation.py` — BAQ default assertions for RNA
  (`apply_baq=True`) and DNA (`apply_baq=False`) modes.

### 🧹 Chores

- Deleted ad-hoc analysis scripts from `scripts/` (compare_tlen_vs_physical,
  concordance, plot_fsd_distributions, plot_fsd_histogram).
- Deleted `scripts/*_test/` directories containing test artifacts.
- Added `.gitignore` patterns for `*.parquet`, `*.mfsd_report.html`, and
  `scripts/*_test/` directories.
- Ruff B904 fix: `raise ... from None` in `mfsd_report.py`.
- Ruff UP015 fix: removed unnecessary `"r"` mode argument from `open()`.

### 🧬 N-Base Diagnostic Integration (Phases 0–3)

#### ✨ Added

- **Diagnostic output columns**: `any_alt`, `partial_alt`, `n_count` appended to
  all MAF output (24 gbcms DNA columns total, up from 21).
- **VCF diagnostic tags**: `AAD` (Any ALT Depth), `PAD` (Partial ALT Depth),
  `NAD` (N-base Depth) emitted in both INFO and FORMAT sections.
- **N-base defense-in-depth**: Explicit N-base guards in `check_snp`, `check_mnp`,
  and `check_complex` — N bases are classified as uninformative regardless of
  reported base quality, preventing silent evidence inflation from duplex-masked
  positions (fgbio) or sequencer failure.
- **Structural invariants**: `any_alt = AD + partial_alt`, `any_alt >= AD`,
  `DP >= RD + AD + partial_alt + n_count` enforced and documented.
- **`trace!`-level diagnostics**: N-base detection, n_count accumulation, and
  partial_alt counting logged at trace level for production debugging.
- **ALT-contains-N rejection**: Variants where the ALT allele contains N are
  rejected with `FAIL_ALT_CONTAINS_N` validation status and `warn!`-level log.

#### 🔄 Changed

- **MNP quality strategy**: Replaced all-or-nothing `min(BQ across block)` gate
  with **masked per-position evaluation** — each discriminating position (REF ≠ ALT)
  is independently assessed; low-BQ and N bases are masked but unmasked positions
  still vote. This recovers reads in GC-rich regions (e.g., TERT promoter) where
  a single low-quality position previously dropped the entire read.
- **MNP ThirdAllele handling**: Mixed-vote reads now track `positions_matching_alt`
  for diagnostic partial_alt counting instead of being silently discarded.

#### 🧪 Tests

- **[NEW]** `tests/test_mnp_concordance.py` — 6 tests for MNP concordance with
  C++ gbcms on production duplex BAMs.
- **[NEW]** `tests/test_phase2_output.py` — 5 tests for diagnostic column presence,
  invariant validation, and N-count sanity on fixture data.
- **26+ Rust unit tests** for N-base masking, MNP per-position evaluation, partial
  match tracking, invariant enforcement, and edge cases (all-N reads, mixed BQ/N).
- Updated `tests/test_column_count_delta_is_three` to assert 24 gbcms DNA columns.

#### 📚 Documentation

- **Updated** `docs/reference/counting-metrics.md` — diagnostic columns, invariant
  tables, VCF tag definitions, N-base handling section.
- **Updated** `docs/reference/output-formats.md` — AAD/PAD/NAD in VCF header,
  INFO/FORMAT tables, annotated example; any_alt/partial_alt/n_count in MAF table.
- **Updated** `docs/reference/allele-classification.md` — SNP flowchart with N guard,
  MNP section rewritten for masked per-position algorithm, Complex N-base note.
- **Updated** `docs/reference/architecture.md` — structural invariants and diagnostic
  output fields in Formulas section.
- **Updated** `docs/development/developer-guide.md` — regression invariant checklist.
- **Updated** `docs/development/testing-guide.md` — Phase 2 test suite, silent
  failures matrix, invariant table, updated test counts.

## [4.0.1] - 2026-03-24

### 🔧 Fixed

- **`was_normalized` flag accuracy**: Split into granular `was_anchor_resolved`
  and `was_left_aligned` flags with backward-compatible `was_normalized` getter.
  Fixes 1150 false negatives (anchor resolution not tracked) and 58 false
  positives (unnecessary anchor+trim round-trip for non-dash complex variants).
  No impact on BAM counting — display/logging only.
- **Left-alignment false positives**: Fixed case-sensitive `modified` check in
  `left_align_variant` to use `eq_ignore_ascii_case`, preventing soft-masked
  FASTA bases from triggering spurious normalization flags.
- **Non-dash anchor resolution**: Narrowed MAF anchor resolution guard to
  dash-allele-only variants. Non-dash complex/deletion MAF variants
  (e.g., `GG>A`) no longer enter the unnecessary anchor+trim cycle.
- **PairHMM pangenome panic**: Fixed unsigned integer underflow
  (`range end index 18446744073709551615`) in pangenomic haplotype
  construction caused by left-to-right delta-adjusted coordinate math.
  Rewrote `build_haplotype_matrix` with right-to-left variant application
  algorithm that eliminates coordinate drift by construction, plus
  power-set sibling combinatorics for true multi-haplotype evaluation.
  Only affects `--alignment-backend hmm`.

### 🔄 Changed

- Normalization logging now shows granular breakdown:
  `"X normalized (Y anchor-resolved, Z left-aligned)"` in both Rust engine
  and Python pipeline/normalize logs.
- `gbcms normalize` TSV output now includes `was_anchor_resolved` and
  `was_left_aligned` columns before the existing `was_normalized` column.
- CLI and reference docs updated with new column descriptions.

## [4.0.0] - 2026-03-20

### ⚠️ Breaking Changes

- **Nextflow module split**: Single `run/main.nf` replaced by three dedicated
  modules — `dna/main.nf`, `rna/main.nf`, `normalize/main.nf`. Consumer
  pipelines must update `include` paths and use `GBCMS_DNA`, `GBCMS_RNA`, or
  `GBCMS_NORMALIZE` process names.
- **Rust `shared/` module**: Common BAM utilities, BAQ, filters, fragment
  logic, and statistics extracted from `counting/` into `shared/`. Any Rust
  consumer crate linking against gbcms internals must update import paths.

### ✨ Added

- **Phase 3 WFA+PairHMM unification** (`feat: unify check_complex Phase 3`):
  Complex indel classification now routes through a unified pangenomic pipeline
  — fast-path WFA alignment with PairHMM fallback. Haplotype matrix
  construction via `pangenome.rs`, WFA routing via `wfa_router.rs`.
  Significantly improves classification accuracy on complex multi-allelic
  variants.
- **RNA mode output columns** (`fix(rna): pass mode= to VcfWriter/MafWriter`):
  `gbcms rna` now correctly emits RNA-specific columns in both VCF and MAF
  output. Previously, all RNA columns (`SEN`, `ANT`, `ASEN`, `RED`, `SPL` in
  VCF; `rna_sense_depth`, `rna_antisense_depth`, `rna_alt_sense_count`,
  `rna_editing_site`, `rna_splice_spanning` in MAF) were silently absent
  regardless of mode. Regression tests added.
- **`gbcms normalize` Nextflow module**: New `normalize/main.nf` for standalone
  variant normalization without counting in Nextflow pipelines.
- **Output Formats reference doc**: `docs/reference/output-formats.md` — complete
  column-level schema reference for VCF and MAF output under all mode/flag
  combinations (DNA vs RNA, with mFSD, with normalization columns).

### 🔧 Fixed

- **Complex indel classification** (`fix: correctly classify complex indels`):
  Fixes for Phase 3 dispatch cases 2, 3, and 4 — previously misclassified
  complex variants where `ref_len ≠ alt_len` and the CIGAR structure doesn't
  map cleanly to pure insertion or deletion.
- **`rna_editing_db` log leakage** (`fix: gate rna_editing_db from DNA mode`):
  DNA mode no longer emits a log line referencing `rna_editing_db`.
- **CliRunner terminal width** (`fix: widen CliRunner terminal`): CI help
  output truncation resolved — `CliRunner(mix_stderr=False, terminal_width=120)`
  prevents Typer/Rich from hiding options in the middle of the params list.
- **Filter defaults documentation**: Secondary, supplementary, and QC-failed
  filters are on by default for DNA mode — corrected in `cli/dna.md`,
  `cli/rna.md`, and `read-filters.md`.
- **RNA mode output pipe wiring** (critical silent bug): `Pipeline._write_output()`
  now passes `mode=self.config.mode` to both `VcfWriter` and `MafWriter`.

### 🏗️ Refactored

- **Rust `shared/` module extraction**: `bam_utils`, `baq`, `filters`, `fragment`,
  `stats` extracted from the `counting/` directory into a new top-level
  `shared/` module, enabling reuse across `counting/` and `normalize/`.
- **`parquet_writer.rs` relocated**: Moved from `counting/` into `shared/` during
  module extraction.

### 📚 Documentation

Major documentation overhaul — 28+ files updated:
- Complete MkDocs plugin utilization pass (tabbed, details, admonitions,
  mermaid, code annotations, glightbox) across all reference pages
- WFA fast-path Phase 3 documented in `allele-classification.md`
- Complex indels guide: RNA compatibility and exon-boundary limitation (D6)
  documented with cross-links
- Architecture module tree corrected for `shared/` and new Nextflow modules
- Filter defaults corrected across `cli/dna.md`, `cli/rna.md`, `read-filters.md`
- Mermaid diagrams fixed: raw unicode escapes removed, `\\n` → `<br/>`,
  backslash-escaped quotes removed
- NEW: `docs/reference/output-formats.md` — authoritative output schema reference
- **Versioned docs assets**: old opaque-named binary files replaced with
  `{name}_{version}.{ext}` convention (`overview_4.0.0.pdf`,
  `allele_classification_4.0.0.pdf`, `read_filter_4.0.0.jpg`); 5 stale files
  deleted; poster references corrected to match page content

### 🧹 Chores

- All Clippy `-D warnings` resolved across `engine.rs`, `rna.rs`,
  `variant_checks.rs`, `pairhmm.rs`
- `ruff`, `black`, `mypy` all pass with 0 errors (38 source files checked)
- Auto-generated mermaid SVGs removed from git (added to `.gitignore`)
- `fallback_to_build_date=true` added to `git-revision-date-localized` plugin
- `test_cli_dna_rna.py`: mypy `attr-defined` fixed by annotating `_click_app`
  as `click.Group`; `Set[str | None]` → `Set[str]` via `if p.name is not None`

### 🧪 Tests

- `tests/test_pipeline_rna.py` (NEW): 7 pipeline-level integration tests for
  RNA mode output — VCF INFO headers, INFO values, FORMAT field, MAF column
  headers, MAF values, and negative DNA assertions
- `tests/test_rna_output.py`: 4 write round-trip tests added (VCF INFO values,
  RED flag on/off, MAF RNA column values)
- `tests/test_maf_preservation.py`: `test_vcf_to_maf_always_uses_sample_name` added
- `tests/test_cli_dna_rna.py`: mypy fixes (click.Group cast, None guard)
- All test `MockCounts` helpers updated with RNA fields

## [3.0.0] - 2026-03-05

### ⚠️ Breaking Changes
- **Package renamed**: `py-gbcms` → `gbcms`. Update your dependencies:
  ```bash
  pip uninstall py-gbcms
  pip install gbcms
  ```
  A final `py-gbcms==3.0.0` deprecation stub on PyPI re-exports `gbcms` and
  issues a `DeprecationWarning` for smooth migration.

### ✨ Added
- **mFSD native integration**: `--mfsd` and `--mfsd-parquet` flags output a
  Parquet file with 31 MAF columns + 7 VCF INFO fields via the Rust native
  Parquet writer. Compression: ZSTD level 1.
- **CLI validation hardening**: 12 validation gaps resolved — fail-fast BAM
  accessibility check, `--lenient-bam` flag for permissive mode, `.vcf.bgz`
  accepted as a valid variant input extension.
- **`py-gbcms` deprecation stub**: `compat/py-gbcms/` shim published to PyPI
  as `py-gbcms==3.0.0` for backwards compatibility during migration.

### 🔧 Fixed
- Parquet output compression switched from SNAPPY to ZSTD(1).
- Resolved `mypy` `AlignedSegment.cigar` attribute errors in test suite.
- Lint and stale-comment cleanup post Phase 6 audit.

### 📚 Documentation
- Full documentation sync with codebase after Phase 6 (mFSD) merge.
- PDF generation guide added to developer documentation.
- `.agent/rules/` directory added; stale `.antigravity` docs removed.

### 🏗️ CI
- Added `mkdocs-print-site-plugin` to CI docs install step and `pyproject.toml`.

## [2.8.0] - 2026-02-23

### ✨ Added
- **PairHMM alignment backend**: Alternative Phase 3 alignment via `--alignment-backend hmm` with probabilistic scoring using base quality probabilities. Configurable LLR threshold (default 2.3 ≈ ln(10)) and gap probabilities for repeat/non-repeat regions. 6 new CLI options. Exposed as first-class Nextflow params in `nextflow.config`
- **MNP min-BQ-across-block quality strategy**: MNP quality now assessed using min(BQ) across the entire block, matching C++ GBCMS `baseCountDNP`. Low-quality MNP reads now fall through to `check_complex` for masked comparison instead of being silently skipped
- **Per-phase ClassifyResult counters**: Diagnostic counters track how many reads are resolved in each classification phase (Phase 1/2/2.5/3)
- **`--trace` flag**: Two-tier Rust logging — `--verbose` for debug, `--trace` for per-read classification diagnostics via pyo3-log

### 🔧 Fixed
- **Phase 2/2.5 overcounting**: Complex variants with `REF >> ALT` now skip Phase 2 Case B and Phase 2.5 when `ref_len > 2 × alt_len` — short ALT trivially matches, edit distance is biased toward shorter allele
- **Phase 3 bypass for pure DEL/INS**: Removed `is_worth_realignment` prefilter — CIGAR structure is ground truth for pure deletions/insertions, prefilter was overcounting
- **DP anchor overlap**: Now uses single-position check (`read_start ≤ variant.pos`), matching Mutect2, VarDictJava, and samtools mpileup standard
- **MNP LowQuality routing**: LowQuality reads now fall through to `check_complex` for masked comparison instead of being skipped entirely
- **S3 underflow guard**: Guards against negative `ctx_offset` in deletion S3 validation

### 🏗️ Refactored
- **Rust module structure**: Split `counting.rs` (1904 LOC) and `normalize.rs` (1221 LOC) into idiomatic module directories: `counting/` (7 modules) and `normalize/` (7 modules)
- **AlignmentBackend threading**: `AlignmentBackend` enum threaded through all Phase 3 call sites

### 📚 Documentation
- **Comprehensive audit**: 28 fixes across 23 files — all GitBook URLs → MkDocs, version templating (X.Y.Z), undocumented CLI options, 5 mermaid diagrams updated for post-2.7.0 logic, architecture module tree refreshed, PairHMM documented end-to-end
- Deleted stale `nextflow/CHANGES.md`

### 🧹 Chores
- Fix all cargo clippy warnings
- Fix ruff lint issues, black formatting
- Update test expectations for new behavior
- CI: skip test workflow for docs-only changes

### 🧪 Tests
- 8 Python alignment backend integration tests
- 5 SW-vs-PairHMM concordance tests
- 10 MNP unit tests
- Multi-allelic isolation, DP/neither, fragment consensus, normalization tests

## [2.7.0] - 2026-02-19

### ✨ Added
- **Phase 2.5 edit distance fallback**: When read reconstruction length matches neither REF nor ALT (e.g., incomplete MAF definition), Levenshtein distance discriminates the closest allele with >1 edit margin safety guard
- **Phase 3 local SW fallback**: Complex variants where semiglobal alignment produces confident-but-wrong calls (e.g., EPHA7 `TCC→CT`) are now rescued via local Smith-Waterman that soft-clips mismatched flanks. Dual-trigger requires both score reversal and ≥2-point margin

### 🔧 Fixed
- **Allele-based dispatch** (`check_allele_with_qual`): Routes by `ref_len × alt_len` instead of unreliable `variant_type` string labels. SNP (1×1), insertion (1×N), deletion (N×1), MNP (N×N equal), complex (N×M unequal) — eliminates misrouting when callers emit inconsistent type annotations
- **SW semiglobal argument order**: Fixed `ref_hap`/`alt_hap` argument swap in semiglobal alignment that was scoring reads against the wrong haplotype
- **Haplotype trimming removed**: Eliminated shared symmetric trim that caused `slice index starts at 7 but ends at 6` panics on asymmetric indels; replaced with validated per-haplotype bounds
- **MNP fallback**: MNP reads now correctly fall through to SW alignment instead of silently returning "neither" on partial mismatches
- **Dual-count guard**: Prevents a single read from being counted as both REF and ALT when SW scores are exactly equal
- **Soft-clip restriction**: Soft-clipped bases no longer incorrectly contribute to variant region reconstruction
- **Strand bias orientation**: Strand bias (Fisher's exact) now couples to the winning allele, not the raw alignment orientation
- **Interior REF quality proxy**: Reads falling entirely within a large deletion (>50bp) now use median base quality instead of 0
- **Interior REF guard removed**: Eliminated the `has_large_cigar_del` guard that massively overcounted REF for large deletions by misclassifying ALT-supporting reads

### 🧹 Chores
- **Clippy**: Removed unused `has_large_cigar_del` variable
- **Tests**: Updated `test_fuzzy_complex::TestLengthMismatch` expectation to reflect Phase 3 local SW fallback behavior

## [2.6.1] - 2026-02-19

### 🔧 Fixed
- **Per-haplotype trimming**: Fixed `slice index starts at 7 but ends at 6` panic in `counting.rs` on asymmetric indels. Replaced shared symmetric trim with independent per-haplotype `trim_haplotype()` function that calculates bounds safely for each allele

### ✨ Added
- **Tolerant REF validation**: Variants with ≥90% REF match against the FASTA are now counted (status `PASS_WARN_REF_CORRECTED`) instead of being silently rejected. The FASTA REF is used for haplotype construction. Variants with <90% match are still rejected as `REF_MISMATCH`

### 📚 Documentation
- **Visual posters**: Added overview, normalization, and read-filter/counting-metrics posters (JPG) to reference documentation pages with lightbox support
- **Embedded PDFs**: Added inline PDF viewer for allele classification guide and detailed overview presentation via `mkdocs-pdf` plugin
- **Variant normalization**: Updated REF validation docs with 3-tier flowchart, `PASS_WARN_REF_CORRECTED` status, and EGFR exon 19 real-world example

### 🔧 CI
- **`deploy-docs.yml`**: Added `mkdocs-pdf` to docs CI pip install dependencies

## [2.6.0] - 2026-02-18

### ✨ Added
- **Adaptive context padding**: Dynamically increases `ref_context` flanking in tandem repeat regions (homopolymer through hexanucleotide). Formula: `max(default, repeat_span/2 + 3)`, capped at 50bp. Enabled by default (`--adaptive-context/--no-adaptive-context`)
- **`gbcms normalize` command**: Standalone variant normalization (left-align + REF validate) without counting, outputs TSV with original and normalized coordinates
- **Nextflow parameters**: `fragment_qual_threshold`, `context_padding`, `show_normalization`, `adaptive_context` now configurable in `nextflow.config`
- **Docs restructure**: Split monolithic `variant-counting.md` into 4 focused pages: Variant Normalization, Allele Classification, Counting Metrics, Read Filters
- **HPC install docs**: Micromamba-based source install with Python 3.13

### 🔧 Fixed
- **Interior REF guard** for large deletions (>50bp): Reads falling entirely within a deleted region are now correctly classified as REF instead of ALT by Smith-Waterman
- **Windowed reciprocal overlap**: Improved shifted indel detection using bidirectional overlap scoring
- **Complex variant counting** (EPHA7 `TCC→CT`): Fixed base quality extraction for all variant-type handlers (`check_insertion`, `check_deletion`, `check_mnp`, `check_complex`)
- **MAF VCF-style conversion**: Corrected complex variant handling in MAF→internal coordinate conversion
- **Lint**: Fixed ruff I001/E402/B905 and black formatting in `pipeline.py`

### 🔄 Changed
- **Dead code removed**: `GenomicInterval` class, `Variant.interval` property, `fragment_counting` config field
- **`fetch_single_base()` refactored**: Delegates to `fetch_region()`, removing 33 lines of duplicated chr-prefix retry logic
- **Release guide**: Updated version locations table with exact line numbers and verification command
- **Nextflow pipeline diagram** added to docs index

## [2.5.0] - 2026-02-12

### ✨ Added
- **`--preserve-barcode` flag**: Keeps original `Tumor_Sample_Barcode` from input MAF instead of overriding with BAM sample name (MAF→MAF workflows)
- **`--column-prefix` parameter**: Controls prefix for gbcms count columns in MAF output (default: none; use `--column-prefix t_` for legacy compatibility) ⚠️
- **`CoordinateKernel`**: Centralized MAF↔internal 0-based coordinate conversion with variant-type-aware logic for SNP, insertion, deletion, and complex variants
- **Nextflow `FILTER_MAF` module**: Per-sample MAF variant filtering by `Tumor_Sample_Barcode` supporting exact match, regex, and multi-select (comma-separated) modes
- **Nextflow `PIPELINE_SUMMARY` module**: Aggregated per-sample filtering statistics with formatted console output
- **Nextflow `--filter_by_sample` parameter** and samplesheet `tsb` column for multi-sample MAF workflows
- **Nextflow documentation**: Samplesheet `tsb` column guide, `--filter_by_sample` parameter reference, multi-sample MAF filtering examples

### 🔧 Fixed
- **Fragment quality extraction** (critical): All variant-type handlers (`check_insertion`, `check_deletion`, `check_mnp`, `check_complex`) now return actual `base_qual` from CIGAR walk instead of 0 — fixes systematic ALT undercount for indels in fragment-level consensus
- **FILTER_MAF heredoc conflict**: Restructured script from Python shebang to `python3 << 'PYEOF'` pattern, resolving `SyntaxError` from bash syntax in Python context
- **FILTER_MAF string quoting**: Changed to single-quoted Python strings for Nextflow variable interpolation to prevent CSV-parsed double-quote conflicts
- **`splitCsv` quote handling**: Added `quote:'"'` parameter for correct RFC 4180 parsing of comma-separated TSB values within quoted CSV fields
- **mypy `no-redef` error**: Removed redundant type annotation in `output.py` else branch

### 🔄 Changed
- **MAF output column prefix default**: Changed from `t_` to empty string (no prefix). Use `--column-prefix t_` for legacy `t_ref_count` / `t_alt_count` style columns ⚠️
- **`MafWriter` refactored**: MAF→MAF path preserves all original columns verbatim; VCF→MAF path builds row from GDC fieldnames with `CoordinateKernel` coordinate conversion
- **Nextflow `GBCMS_RUN` input**: Variants bundled into sample tuple `(meta, bam, bai, variants)` instead of separate channel
- **Nextflow `GBCMS` workflow**: Simplified to 2-channel interface (`ch_samples`, `ch_fasta`) from 3 channels
- **Nextflow config**: Added `column_prefix`, `preserve_barcode`, `filter_by_sample` parameters

## [2.4.0] - 2026-02-10

### ✨ Added
- **Fragment Consensus Engine**: `FragmentEvidence` struct with u64 QNAME hashing and quality-weighted R1/R2 consensus; ambiguous fragments are discarded (not assigned to REF)
- **`--fragment-qual-threshold`**: New CLI option (default 10) controlling consensus quality difference for fragment conflict resolution
- **Windowed Indel Detection**: ±5bp positional scan with 3-layer safeguards (sequence identity, closest match, reference context validation)
- **Quality-Aware Complex Matching**: Masked comparison that ignores bases below `--min-baseq`; 3-case comparison (equal-length, ALT-only, REF-only) with ambiguity detection
- **Variant Counting Guide**: New `docs/reference/variant-counting.md` with algorithm diagrams for all variant types (~700 lines)
- **MAF Normalization Docs**: Added indel normalization and coordinate handling to `docs/reference/input-formats.md`
- **47 Tests**: Up from 16 — added `test_shifted_indels.py` (15), `test_fuzzy_complex.py` (14), `test_fragment_consensus.py` (4)

### 🔄 Changed
- **`--min-baseq` default**: `0` → `20` (Phred Q20) — activates quality masking by default for improved accuracy on low-coverage samples ⚠️
- **`--version` flag**: Added to CLI (`gbcms --version`)
- **Deploy-Docs Workflow**: Replaced `mkdocs gh-deploy` with `mike` for multi-version documentation; deploys `stable` (tagged version) from main and `dev` from develop branch; added `extra.version.provider: mike` to `mkdocs.yml` with version switcher widget

### 🔧 Fixed
- **Fragment double-counting bug**: R1+R2 pairs previously counted as two independent observations; now collapsed via quality-weighted consensus
- **MAF Input Hardening**: Graceful handling of missing/malformed fields with warnings instead of crashes
- **CI Release Pipeline**: Stabilized manylinux builds — migrated from manylinux_2_28 to manylinux_2_34, resolved OpenSSL/CURL vendor conflicts via `docker-options` pattern
- **Type stubs**: `_rs.pyi` and `gbcms_rs.pyi` synced with Rust bindings (added `ref_context`, `ref_context_start`, `fragment_qual_threshold`)
- **Linting**: All files pass black, ruff, and mypy

### 📚 Documentation
- Architecture comparison table updated with windowed indels and masked comparison
- Nextflow config and docs updated with new `min_baseq` default
- Testing guide expanded with Phase 2a/2b test files
- HPC/RHEL 8 installation instructions updated with `clangdev` header management
- Release guide updated with docs version locations
- `.antigravity` project files updated with current Rust LOC (~1270) and test counts (47)

## [2.3.0] - 2026-02-06

### ✨ Added
- **Nextflow BAI Auto-Discovery**: Checks `.bam.bai` and `.bai` extensions automatically
- **Documentation Modernization**: Hierarchical navigation, glightbox, panzoom, abbreviations
- **Performance Benchmarks**: cfDNA duplex sample metrics in documentation
- **RHEL 8 Installation Guide**: Conda-based source installation for legacy Linux

### 🔄 Changed
- **Dockerfile**: Added `procps`, `bash`, OCI labels, `maturin[patchelf]`, selective COPY
- **Nextflow Config**: `--platform linux/amd64`, shell config, local profile, observability (trace/report/timeline/dag)
- **MkDocs**: Switched to `navigation.sections`, 20+ abbreviations with hover tooltips
- **GitHub Actions**: Consolidated deploy-docs workflows, added caching and PR validation
- **CI Wheels**: Migrated from `manylinux_2_28` to `manylinux_2_34` (AlmaLinux 9 with OpenSSL 3.0+)

### 🔧 Fixed
- **Nextflow**: Empty `--suffix` argument no longer causes failures
- **Admonitions**: Converted GitHub-style alerts to MkDocs syntax
- **CI Build**: Resolved `curl-sys` OpenSSL version conflict by switching to manylinux_2_34

## [2.2.0] - 2026-02-04

### ✨ Added
- **Multi-platform Wheel Publishing**: Maturin-based CI builds for Linux (x86_64, aarch64), macOS (Intel, Apple Silicon), and Windows
- **Structured Logging**: New `utils/logging.py` module with Rich console output, timing utilities, and log file support
- **Mermaid Diagrams**: Architecture documentation with interactive flowcharts
- **Release Guide**: Comprehensive `docs/RELEASE.md` with git-flow workflow

### 🔄 Changed
- **Folder Restructure**: Moved Rust code to `rust/` (bundled as `gbcms._rs`)
- **Config Hierarchy**: Nested Pydantic models (`ReadFilters`, `QualityThresholds`, `OutputConfig`) for better organization
- **Code Quality**: Added `__all__` exports, docstrings, and type hints across all modules
- **StrEnum**: Modern enum pattern with Python 3.10 backport

### 📚 Documentation
- New `docs/ARCHITECTURE.md` with system diagrams
- New `docs/DEVELOPMENT.md` (developer guide)
- New `docs/TESTING.md` (testing guide)
- Updated MkDocs with mermaid2 plugin and snippet includes

## [2.1.2] - 2025-11-25

### 🔧 Fixed
- **PyPI Distribution**: Fixed source distribution size issue by correctly excluding large files (tests, docs, etc.) via `pyproject.toml` configuration.

## [2.1.1] - 2025-11-25 [YANKED]

!!! warning "Yanked Release"
    This release was yanked from PyPI due to a source distribution size limit error. Use 2.1.2 instead.

### 🔧 Fixed
- **PyPI Distribution**: Added MANIFEST.in (failed to work with Hatchling) to reduce source distribution size
- **Documentation**: Added comprehensive Installation guide
- **Documentation**: Unified Contributing guide (merged code + docs contributions)
- **Documentation**: Added Changelog to documentation navigation

## [2.1.0] - 2025-11-25

### ✨ Added

#### Nextflow Workflow
- **Production-ready Nextflow workflow** for processing multiple samples in parallel
- **SLURM cluster support** with customizable queue configuration
- **Per-sample suffix support** via optional `suffix` column in samplesheet
- **Docker and Singularity profiles** for containerized execution
- **Automatic BAI index discovery** with validation
- **Resume capability** for failed workflow runs
- **Resource management** with automatic retry and scaling
- **Comprehensive documentation** in `docs/NEXTFLOW.md` and `nextflow/README.md`

#### Documentation
- **Usage pattern comparison** guide (`docs/WORKFLOWS.md`) for choosing between CLI and Nextflow
- **MkDocs integration** for beautiful GitHub Pages documentation
- **Local documentation preview** with live reload (`mkdocs serve`)
- **Staging deployment** from `develop` branch for testing docs
- **Production deployment** from `main` branch
- **Reorganized documentation structure** with clear CLI vs Nextflow separation
- **CLI Quick Start guide** (`docs/quick-start.md`)

### 🔧 Changed
- **Documentation workflow**: docs now live on `main` branch with automated deployment
- **GitBook integration**: configured to read from `main` branch
- **Nextflow module**: improved parameter passing with meta.suffix support

### 📝 Documentation
- Complete Nextflow workflow guide with SLURM examples
- Per-sample suffix usage examples
- Git-flow documentation workflow guide
- Local preview instructions
- Updated README with clear usage pattern separation

## [2.0.0] - 2025-11-21

### 🚀 Major Rewrite

Version 2.0.0 represents a complete rewrite of py-gbcms with a focus on performance, correctness, and modern architecture.

### ✨ Added

#### Core Features
- **Rust-based Counting Engine**: Hybrid Python/Rust architecture for 20x+ performance improvement
- **Strand Bias Statistics**: Fisher's exact test p-values and odds ratios for both reads (`SB_PVAL`, `SB_OR`) and fragments (`FSB_PVAL`, `FSB_OR`)
- **Fragment-Level Counting**: Majority-rule fragment counting with strand-specific counts (`RDF`, `ADF`)
- **Variant Allele Fractions**: Read-level (`VAF`) and fragment-level (`FAF`) allele fraction calculations
- **Thread Control**: Explicit control over parallelism via `--threads` argument (default: 1)

#### Input/Output
- **VCF Output Format**: Standard VCF with comprehensive INFO and FORMAT fields
- **MAF Output Format**: Extended MAF with custom columns for strand counts and statistics
- **Column Preservation**: Input MAF columns are preserved in output
- **Multiple BAM Support**: Process multiple samples via `--bam-list` or repeated `--bam` arguments
- **Sample ID Override**: Explicit sample naming via `--bam sample_id:path` syntax

#### Filters
- `--filter-duplicates`: Filter duplicate reads (default: enabled)
- `--filter-secondary`: Filter secondary alignments
- `--filter-supplementary`: Filter supplementary alignments
- `--filter-qc-failed`: Filter reads that failed QC
- `--filter-improper-pair`: Filter improperly paired reads
- `--filter-indel`: Filter reads with indels in CIGAR

#### CLI & Usability
- **Modern CLI**: Built with Typer and Rich for beautiful terminal output
- **Progress Tracking**: Real-time progress bars and status indicators
- **Direct Invocation**: Use `gbcms run` instead of `python -m gbcms.cli`
- **Output Customization**: `--suffix` flag for output filename customization
- **Flexible Input**: Support for both VCF and MAF input formats

#### Infrastructure
- **Docker Support**: Production-ready multi-stage Dockerfile with optimized layers
- **Type Safety**: Full type annotations with mypy support
- **Type Stubs**: Provided `.pyi` stub file for Rust extension
- **Comprehensive Tests**: Extended test suite with accuracy and filter validation
- **CI/CD**: GitHub Actions workflows for testing, linting, and releases

### 🔄 Changed

#### Architecture
- Migrated from pure Python to hybrid Python/Rust architecture
- Core counting logic implemented in Rust using `rust-htslib`
- Data parallelism over variants with per-thread BAM readers

#### Output Formats
- **VCF FORMAT fields**: Strand-specific counts now use comma-separated values (e.g., `RD=5,3` for forward,reverse)
- **MAF columns**: Standardized column names (`t_ref_count_forward`, `t_alt_count_reverse`, etc.)
- **Coordinate System**: Internal 0-based indexing with correct conversion for VCF (1-based) and MAF output

#### Performance
- **Speed**: 20x+ faster than v1.x on typical datasets
- **Memory**: Efficient per-thread BAM readers with minimal overhead
- **Scalability**: Configurable thread pool for optimal resource usage

#### Dependencies
- **Python**: Updated to require Python ≥3.10
- **Rust**: pyo3 0.27.1, rust-htslib 0.51.0, statrs 0.18.0
- **Python Packages**: pysam ≥0.21.0, typer ≥0.9.0, rich ≥13.0.0, pydantic ≥2.0.0

### 🗑️ Removed

- **Legacy Python Counting**: Pure Python implementation removed in favor of Rust
- **Old CLI**: Deprecated `python -m gbcms.cli` entry point
- **Unused Dependencies**: Removed `cyvcf2` and `numba` (no longer needed)
- **Pre-commit Hooks**: Removed in favor of explicit linting in CI

### 🐛 Fixed

- Correct handling of complex variants (MNPs, DelIns)
- Proper strand assignment for fragment counting
- Reference validation against FASTA for all variant types
- Thread-safe BAM access with per-thread readers

### 📚 Documentation

- Complete rewrite of all documentation
- New guides: `INSTALLATION.md`, `CLI_FEATURES.md`, `INPUT_OUTPUT.md`
- Comprehensive API documentation
- Docker usage examples
- Contributing guidelines updated

### 🔧 Technical Details

#### Rust Components
- `gbcms._rs`: PyO3-based Rust extension (bundled in wheel)
- Fisher's exact test via `statrs` crate
- Rayon-based parallelism with configurable thread pools
- Safe memory management with Rust's ownership model

#### Testing
- 16 comprehensive test cases
- Accuracy validation with synthetic BAM files
- Filter validation for all read flag combinations
- Integration tests with real-world data

### ⚠️ Breaking Changes

Version 2.0.0 is **not backward compatible** with 1.x. Key breaking changes:

1. **CLI syntax**: Use `gbcms run` instead of `python -m gbcms.cli`
2. **Output format**: VCF/MAF column structures have changed
3. **Default behavior**: Only duplicate filtering enabled by default (was: all filters)
4. **Dependencies**: Requires Rust toolchain for installation from source
5. **Python version**: Minimum Python 3.10 (was: 3.8)

### 📦 Installation

```bash
# From PyPI (includes pre-built wheels)
pip install gbcms

# From source (requires Rust)
pip install git+https://github.com/msk-access/gbcms.git

# Docker
docker pull ghcr.io/msk-access/gbcms:2.0.0
```

### 🙏 Acknowledgments

This rewrite was designed and implemented with a focus on correctness, performance, and modern best practices in bioinformatics software development.

---

## [1.x] - Legacy

Previous versions (1.x) used a pure Python implementation. See git history for details.
