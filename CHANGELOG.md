# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed — release infrastructure (#136, #137, #152)

- **One version source per ecosystem (#136).** The Python package takes its version
  from `rust/Cargo.toml` (maturin, `dynamic = ["version"]`; `gbcms.__version__` reads the
  installed metadata), and the Nextflow modules and banner read `manifest.version`. A
  release edits `rust/Cargo.toml`, the Nextflow manifest and the CHANGELOG instead of
  11 places in two formats. `scripts/release.py check` runs on every PR, and with
  `--tag` it is the release workflow's first job: nothing builds or publishes unless
  the tag is a bare `X.Y.Z` equal to Cargo.toml, its lock and the manifest, with a
  dated CHANGELOG section (a tag on a dev version would have published a dev wheel
  under the release). On `develop` the manifest must name a published image older
  than the dev package, and the package must move to the next dev version after a
  release. Only a tag push publishes (a manual run builds and checks), the tag reaches
  shell steps through `env:`, and the workflow's default token reads only. The
  Nextflow lint job checks every process `nextflow inspect` resolves runs the
  manifest's image; `gbcms --version` and the resolved images are unchanged.
- **The release workflow creates the GitHub Release (#137).** After PyPI and the image
  are published, it creates (or updates) the Releases entry: the title and notes from
  the tag's CHANGELOG section (`## [X.Y.Z] - YYYY-MM-DD — summary`; cut at a section,
  with a link to the full file, past GitHub's body limit), the wheel and sdist attached
  with `SHA256SUMS`, and a build-provenance attestation for each artifact.
- **The bin window is documented as the floor it is (#152).** `BIN_WINDOW` (10 kb) is a
  minimum span that dense variants chain past, bounded by the 200-variant cap; measured
  on real inputs (per-sample lists max 29 kb; one input holding a cohort's variants max
  47 kb) and never capped by span. The page no longer cites the retired parity suite.

### Fixed — an insertion's ALT needs one of the read's own inserted bases read (#240)

- **Reads whose inserted bases nobody can read counted ALT.** A read carrying an
  insertion of the right length whose bases are all N (fgbio masks a duplex
  disagreement to N at Q2) or all below `--min-baseq` was credited ALT on length
  alone: the strict path found the bases unverifiable and handed the read to Phase 3,
  where REF pays for the gap, under both backends. Shifted placements went the same
  way, and the deleted-anchor and split-op reading accepted a readable flank base.
  On the RC DNA set that was 127 reads, mostly 1bp insertions into homopolymer runs
  with a duplex-masked N inside the run (the aligner writes the N as the insertion);
  one row went from 430 ALT reads to 368, most rows lost 0–8. On the FORTE sample, 3
  reads in DNA mode and 8 in RNA mode, all low-quality letters.
- **The rule (RJ-20).** A read whose inserted bases are all unreadable is neither,
  with partial evidence, on every path: at the junction, at another placement inside
  the discrimination window, across several ops, after a deleted anchor, and as a
  truncation. Outside the window such an insert is a separate event, as a readable one
  is. "Readable" means not N and at or above `--min-baseq`, the gate SNV bases pass.
  A read whose readable inserted bases match the ALT stays ALT, however many are
  masked. The read census holds the same rule. The post-splice path already counted
  such reads this way.
- **Community practice.** Likelihood callers (GATK, bcftools, Strelka2) credit ALT on
  length; sequence-keyed counters do not: VarDict and freebayes drop such a read,
  bam-readcount and LoFreq report it as its own allele, and GetBaseCountsMultiSample
  counts it as neither. gbcms counts the given allele from the read's own bases.

### Fixed — Fisher's exact test at any depth (#239)

- **Strand-bias p-values were 0 at depth.** `strand_bias_p_value`,
  `fragment_strand_bias_p_value` (VCF `SB_PVAL`, `FSB_PVAL`) and merge's combined
  columns were 0, as if strongly biased, whenever the binomial C(n, forward total)
  passed f64's range. That needs both strands well covered: from about 1,030 reads
  on a strand-balanced table. Tables with one thin strand stayed correct. On the
  RC cfDNA rows that was 316 of 1,060 read-level p-values, 298 of them truly
  ≥ 0.05. The ASJD junction test (p < 0.05) read the same zeros.
- **The fix.** The test is now the exact two-sided p of R's `fisher.test`, computed
  in log space from the ratio of neighbouring tables, walking out from the most
  likely table: no factorial, nothing overflows, relative error near 10⁻¹³ at any
  depth, and tables too unlikely to count are never visited.
  - Ties follow R's rule (probability ≤ observed × (1 + 10⁻⁷)). That replaces an
    absolute 10⁻¹⁰ tolerance which floored strongly biased tables near 10⁻¹⁰; one
    printed 2.1 × 10⁻¹⁰ where the exact p is 7.0 × 10⁻⁹⁰.
  - The p = 1 guards (≤ 1 ALT observation, empty or full margins) and the odds ratio
    are unchanged.
  - The p-value stays the exact p of the raw counts. At deep coverage, read the odds
    ratio for the size of a bias.
- **Measured on the local FORTE sample (9,536 rows).**
  - What moved: 397 read-level and 205 fragment-level strand-bias p-values, and 5
    ASJD p-values, of which one flag clears. All counts and Parquet tables are
    unchanged.
  - Every reported strand-bias p now equals scipy's exact value.
- **statrs** leaves the direct dependencies.

### Changed — Rust dependencies on their current releases (#139)

The Rust dependencies move to their current releases, one per commit
(`CYCLE_6.6.0_PLAN.md`, "D4 PR B"). No output changes. bio stays at 3.0 for now.
- **pyo3 0.29** (from 0.27). The API migration: `allow_threads` is now `detach`,
  and `Variant` and `BaseCounts` opt in to the by-value conversion that Python
  passes back. The module keeps the GIL on a free-threaded interpreter, as before;
  gbcms is not built or tested free-threaded.
- **rust-htslib 1.0, arrow/parquet 60**, and the semver-compatible updates.
  (statrs went to 0.19 here, then left the dependencies with #239.)
  - parquet 60 truncates column-chunk min/max statistics of strings longer than 64
    bytes. The data are unchanged.
- **wfa2lib-rs moves to upstream's current commit without its default features.**
  Those only built its benchmark binary, so clap, tracing and mimalloc leave the
  build.
- **bio is held at 3.0.** bio 4.1 fixes how the PairHMM scores the bases of an
  extended gap.
  - For reads whose inserted bases are N or below the base-quality floor, that moves
    gbcms's LLR across the ±2.3 threshold (synthetic reads: 1.8 → 5.0).
  - gbcms already credits ALT to such reads at non-repeat insertions, so bio 4 waits
    for a gate that requires the inserted bases to carry the ALT.
  - bio 4's other change, the Smith-Waterman gap rule, needs the open penalty moved
    from −5 to −6 to keep every score; that ships with the bio upgrade.
- **Unchanged output:**
  - Local FORTE sample (RNA and DNA modes, both alignment backends, 9,536 rows
    including 806 indels observed in the reads): compared after step 1, after steps
    1–5, and for the final tree. No column and no Parquet table changes.
  - RC sets on the cluster (28 cfDNA DNA, 33 RNA truth, 3 probe samples and 80 WES
    loci; 2,212 rows): byte-identical.

### Changed — dependencies, CI coverage and the image (#139)

Measured first (`CYCLE_6.6.0_PLAN.md`, "D4 PR A"). Nothing changes at runtime.
- **The dependency floors are true.** At the declared floors the CLI crashed: typer
  0.9.0 cannot read its `list[...]` options. pysam 0.21.0 has no macOS arm64 wheel
  and an sdist that fails on current setuptools. Bisected one package at a time, the
  minimums are `typer>=0.15.4` (older typer does not cap click and breaks with
  click >= 8.2) and `pysam>=0.22.0`; rich 13.0, pydantic 2.0 and polars 1.0 held.
  All of them together pass the suite on Python 3.10 and 3.11. A new CI leg installs
  the floors on every PR and checks that they are what got installed.
- **click is no longer a runtime dependency.** typer bundles its own copy, and gbcms
  never imports click; the CLI tests do, so it moved to the `test` group, and the
  image no longer installs it.
- **CI tests the supported range.** Python 3.10–3.14 all pass on the latest
  releases, but CI tested only 3.11 and 3.12. PRs now run Ubuntu 3.10 (at the
  floors), 3.11 (the image's) and 3.14, plus macOS 3.12, and the classifiers list
  3.10–3.14. A monthly workflow (`latest-deps.yml`, also run by hand before a release)
  runs every version on the newest releases with the semver-compatible Rust updates. On failure it opens or updates
  one `latest-deps` issue.
- **One dev-dependency list.** The `dev` extra and a PEP 735 `dev` group had drifted
  apart. The group lacked pytest-mock, types-pyyaml, pyyaml and mkdocs; the extra
  lacked pyarrow, which the tests import. They are now the dependency groups `test`
  and `dev` (`dev` includes `test` and `docs`, the site's tools, which
  `deploy-docs.yml` now installs instead of its own list), kept out of the published
  metadata.
  `pip install gbcms[dev]` no longer exists. `maturin develop` installs `dev`, which
  needs pip >= 25.1 (or uv); the developer docs now say so, since the old
  `python -m venv` steps failed on Python's bundled pip. `make setup` installs it too.
- **The image installs a hash-pinned lock.** It resolved its dependencies afresh at
  each build. `docker/requirements.lock` (linux/amd64, Python 3.11) is installed
  with `--require-hashes --no-deps`, then `pip check`, and is refreshed at each
  release (release guide step 4). CI now runs the image (`--version`, `pip check`).

### Changed — GTF loading (#139)

Measured on GRCh38.111 (`CYCLE_6.6.0_PLAN.md`, "D4 GTF loading"). Only `--gtf` runs
are affected, and their output is unchanged for well-formed GTFs.
- **The GTF loads in 2 s instead of 9 s, plain or `.gtf.gz`.** noodles-gtf parsed
  every line before reading its feature column. A byte-level parser reads each
  line's feature and chromosome first and checks only the exon lines it keeps.
  Measured in `gbcms rna` on Ensembl 111: 1.6 s for 16 chromosomes (8.9 s before),
  2.0 s for the whole genome (9.3 s), 2.4 s from `.gtf.gz`. GENCODE v50 loads in
  4.0 s (basic) and 6.5 s (comprehensive) from `.gtf.gz`. Peak memory falls (407 MB
  vs 470 MB for 16 chromosomes): `gene_id`, stored for every exon and never read,
  is gone.
  - It loads exactly the exons noodles did: no line differs across Ensembl 111 and
    GENCODE v50 comprehensive, v50 basic and v47lift37 (26.1M lines).
  - gzip and BGZF are detected from the file's content. The config already accepted
    `.gtf.gz`, but the parser failed on it with "stream did not contain valid UTF-8".
- **Malformed exon lines are rejected and warned about.** The grammar is noodles-gtf's,
  column by column, except in four places where noodles loaded a line wrongly:
  - A whitespace run may now separate a key from its value. noodles read
    `transcript_id  "T1"` as the ID ` "T1"`, quotes included.
  - Unquoted values are trimmed.
  - A line with start after end is rejected; noodles kept it as a negative-length exon.
  - A coordinate past 2,147,483,647 is rejected; noodles wrapped it negative.
  - An exon with an empty `transcript_id` is rejected, as a missing one already was.
    noodles merged every such exon into one transcript.
  - A `transcript_id` that strips to nothing (`"."`, `".1"`) is rejected the same way.
  - An exon line that is not UTF-8 text is rejected; it used to stop the run.
  - Each of these lines used to move `exon_boundary_dist` with a bogus boundary.
    Rejected lines now get one warning with their count and the first one's line
    number and reason.
- **The GTF index cache is deprecated.** On Ensembl it saved 1–1.4 s a sample: a hit
  cost 0.9 s in the CLI (~0.3 s of it start-up) against a 1.6–2.0 s load now. It
  brought a build step, a Nextflow process and a 97 MB file per chromosome set. It is
  no longer used:
  - `--gtf-cache-dir`, `gbcms build-gtf-cache` and the Nextflow `--gtf_cache` are
    accepted and ignored with a warning, and are removed in 6.7.0.
  - The `GBCMS_BUILD_GTF_CACHE` process, the cache module and its binding are gone.
  - The noodles-gtf, bincode and serde dependencies are gone too. bincode 3.0 is a
    `compile_error!` release, so there was no upgrade path.

### Changed — mFSD statistics (#153, #154)

Measured first on ACCESS plasma labeled by the patient's buffy coat (operator
decisions 2026-09-25, refined 2026-10-05): mFSD is graded, plasma-only evidence,
never an origin call; the CH-vs-tumor prediction belongs to a separate model that
reads `<sample>.fsd.parquet`. Only `--mfsd` output changes.
- **The fragment-size LLR is the mean per fragment (#153).** `mfsd_alt_llr` and
  `mfsd_ref_llr` (VCF `MFSD_ALT_LLR`, `MFSD_REF_LLR`) were sums, which grew about
  22x with depth on real cfDNA; the mean stays comparable, with n in
  `mfsd_*_count`. An empty fragment class writes `NA` for its mean size and LLR (it
  wrote 0, where the docs already said `NA`; the N class is empty on most rows).
- **The mFSD report grades the evidence (#154).** Its classes are `LEANS-SOMATIC`
  (ALT fragments significantly shorter than REF: KS q < 0.05 with the ALT ECDF
  above the REF ECDF where they differ most), `NO-SIZE-EVIDENCE` and
  `INSUFFICIENT`, replacing `TUMOR-LIKE` / `CH-LIKE` / `AMBIGUOUS`. `CH-LIKE` rested
  on a non-significant KS, which is also the usual result for tumor variants at
  these fragment counts: on 39 labeled samples it called 3 of 17 tumor variants in
  CH-associated genes CH-LIKE, and about 30% of tumor variants below 50 ALT
  fragments looked REF-like. Nothing leans CH now. Gene membership is a note, not
  a gate (a tumor TP53 variant could never be `TUMOR-LIKE`), and the enrichment
  > 1.3 gate is gone: on duplex fragments the new rule leans somatic on 22 of 97
  tumor variants and 0 of 39 white-cell variants (the old rule found 16). The
  direction is read from the KS gap, not from the sub-nucleosomal share, which
  missed an ALT 31 bp shorter with no fragment under 150 bp. A significantly
  longer ALT is not somatic evidence. Variants below `--mfsd-report-min-alt` are
  left out of the report, as before.
- **The KS p-value is exact** up to 10⁷ lattice cells, with Stephens' corrected
  asymptotic series above. The switch at `n·m` 10,000 sent most real ALT-vs-REF
  pairs (348 of 512) to the uncorrected series, which overstated p 1.7–2.3x at 5
  ALT fragments near p = 0.05 and up to 45x for a strong shift. The exact value
  compares lattice points with the observed deviation in integers (a float band
  could count the observed point as not reaching D near the cap), computes p
  directly as the share of paths leaving the band (a tiny p keeps its digits), and
  visits only the cells inside the band.
- **A variant on a contig absent from the BAM stays out of the mFSD q-values.**
  With `--mfsd` its fields held 0.0, which read as a KS test with p = 0 and
  deflated every real variant's BH q-value; it now reads as no test (`NA`).
- **`mfsd_alt_confidence` is `TESTABLE` / `SPARSE` / `NONE`** (was `HIGH` / `LOW` /
  `NONE`): it names how much ALT data there is (≥ 5, 1–4, 0 fragments); at 5
  fragments a real size shift is detected only about 8% of the time.

### Changed — merge, outputs and observability (#194, #221, #223, #224, #129, #148, #130, #131, #156, #220, #225)

Measured first (operator decisions 2026-09-30 for #194, 2026-10-05 for the rest);
community practice surveyed (bcftools, Picard/htsjdk, GATK, samtools/htslib,
maftools, genotype_variants, Snakemake, Nextflow, bam-readcount, LoFreq, fgbio,
DeepVariant, Strelka2). Where no tool sets a standard, gbcms now does more.
- **A missing count is not a zero in `gbcms merge` (#194).** When either flavor's
  count cell is missing or not a finite number (empty, `NA`, `nan`, `inf`, text)
  in a row it has, the combined `simplex_duplex_*` cell is `NA`, as are the
  totals, VAFs and strand bias built from it, and merge warns once per column
  with the number of rows (it summed them as 0, silently). A row an input lacks
  still counts 0 for it. Combined VAFs are written with four decimals, as the
  writers write theirs (`f"{v:.4f}"`), and the combined strand-bias p-values and
  odds ratios as theirs (`1.7045e-01`, `3.0000`). No real output carries such cells (0 of 180,348 count
  cells measured); they come from edited files or other tools.
- **Merge keeps every gbcms column per input (#223).** The mFSD and RNA columns
  (63 of the 89 columns the writers can emit) were taken from the first input
  only, unprefixed, and every later input's were dropped: with `--mfsd`, which
  the pipeline allows alongside merge, a merged ACCESS output showed the duplex
  BAM's fragment sizes as the sample's. They are prefixed per input
  (`duplex_mfsd_ref_mean`), the set taken from the writer, with a test that it
  matches in every mode. Inputs written with `--column-prefix` (`duplex_` as the
  pipeline runs it, `t_`) carry their counts under that prefix and their status,
  strand-bias, mFSD and RNA columns unprefixed: merge renamed nothing for them, so
  the second input's status and strand bias were dropped too, and `t_` counts were
  never combined. Every column is now found under the writer's name.
- **A row only a later input has keeps its annotations (#221):** its first-input
  columns (gene, sample barcode, classification...) come from the earliest
  later input that has the row (matched by that input's own row, so one variant
  listed for two samples keeps each sample's); they were empty. No work when
  every row is in the first input, as in the pipeline.
- **Merge reads its inputs' provenance (#129).** The merged MAF starts with its
  own `#gbcms`/`#command` lines and one `#input` line per input with that
  input's version line (it had none). Merge warns when the inputs come from
  different versions or builds (builds differ only when both name a commit), or
  when only some say which version wrote them, and stops when one is a VCF-input
  MAF from before 6.5.0 (`vcf_pos` without `vcf_ref`/`vcf_alt`, its version line
  missing or older than 6.5.0; a later output can carry vcf2maf's `vcf_pos`) and
  another is not: the
  same 102k VCF records genotyped by 6.4.0 and by this version differ in 6.3% of
  their rows, which merged into 12,771 half-empty rows with exit 0. The INFO line
  that claimed "n/m variants have no <type> counts" counted rows whose REF count
  was 0; it now counts the rows each input lacks. No surveyed merger compares
  producer versions.
- **Merge no longer sums two alleles where the flavors' MNP rescue differs
  (#224, the rescue half of #128, pulled in from 6.7.0).** When one flavor
  reports a rescued component and the other the MNP (or rescue ran on one flavor
  only), the row's combined `simplex_duplex_*` cells are `NA`; they added the two
  alleles' counts, with only a warning. A row one flavor lacks is not mixed. 20
  real ACCESS pairs had none; the warning names each such row.
- **Builds name their commit.** Provenance lines (`#gbcms`, VCF `##source`, run
  logs, `gbcms --version`) read `gbcms v6.6.0.dev0 (9c371263)`; develop now
  carries a `.devN` version, since every development build since 6.5.0 reported
  6.5.0. The commit comes from git, or from `GBCMS_BUILD_COMMIT` (set by the
  Dockerfile and the release workflow).
- **VCF output of a whole-contig deletion is valid (#220, pulled in from
  6.7.0).** A MAF deletion at Start 1 spanning its whole contig has no reference
  base before or after it; it was written with a REF padded past the contig end
  with `N` (`ACGTACN > N`). It is now the symbolic `<NON_SEQUENCE>` record (REF
  the base at POS), in counting runs and `gbcms convert`; the row was already FAIL
  (`FETCH_FAILED`). None in real data.
- **No partial output files (#148).** Every output (MAF, VCF, the merged MAF,
  `convert` and `normalize` files, both Parquet files, the mFSD report) is
  written to `.<name>.partial`, fsynced, and renamed into place; a failed run,
  including a failure while closing (a full disk), leaves nothing at the output
  path and no temp file, and the observations Parquet goes into place only after
  the MAF/VCF. A symlinked output keeps its link (its target is replaced), a
  replaced file keeps its permission bits, and a device or FIFO (`/dev/stdout`)
  is written in place. It left a truncated file
  under the final name (htslib tools, GATK and Picard do too; only workflow
  managers clean up).
- **Per-BAM warnings once (#130).** The `--rescue-mnp` re-count no longer repeats
  the records-without-bases, records-without-qualities and absent `--umi-tag`
  warnings, and those records are counted once each (the count was an upper
  bound: overlapping bins fetched a record more than once).
- **The run start says what the run does (#131).** One INFO block lists every
  resolved option (generated from the configuration, so none is left out) and
  what each count- or column-changing option implies; one INFO line per BAM
  gives facts from its first 20,000 records (duplicates flagged, base-quality
  values, the share below `--min-baseq`, ALT contigs, hard-clipped primaries). It
  warns only when `--min-baseq` removes more than 10% of sampled bases, the
  header has ALT contigs, or more than 1% of primaries are hard-clipped: on 10
  MSK BAMs the default removes 0.8–2.1% and none has ALT contigs or hard-clipped
  primaries. Unmarked duplicates are reported, not warned: ACCESS consensus and
  FORTE RNA BAMs carry none by design. The CLI's four-setting `Config:` line is
  replaced by the block. Implications are stated against the mode's defaults (an
  RNA run's are RNA's), RNA amplicon and strandedness settings included; turning
  off the secondary or supplementary filter lets those alignments join fragment
  evidence, never read counts. A BAM the facts cannot read is left to counting.
- **The docs toolchain is pinned below MkDocs 2.0 (#225, the pin half of #138,
  pulled in from 6.7.0):** `mkdocs>=1.6,<2` and `mkdocs-material>=9.5,<10` in the
  dev extras and the docs workflow (both installed unpinned), with a test; the
  MkDocs 2.0 migration stays in 6.7.0.
- **One page for every QC flag, each defined once (#156):**
  `docs/reference/qc-flags.md`, a table per family (status reasons, diagnostics,
  rescue outcomes, ASJD, QC columns, mFSD classes, VCF record shapes) with mode,
  MAF column, VCF field and what to do. Each family is a snippet section that the
  pages needing the table include (normalization, architecture, RNA annotation,
  output formats, mFSD report); the rest link to it. Tests fail when a flag the
  code emits is missing from the page, when a flag is defined anywhere else (a
  table row or bullet), when a line lists three or more flags without linking it,
  or when an include names a missing section; the docs build now fails on a
  missing snippet or anchor. Consolidating found copies that had drifted: the mFSD
  report's class table and flowchart (thresholds 1.0 and the raw KS p-value; the
  code uses 1.3/1.2, the FDR q-value and the CH-gene check), its `mfsd_ch_flag`
  ("CH-like profile"; it marks a CH gene), output formats' `STRAND_DISCORDANT`
  (from before intronic loci got a strand), the VCF `GS` row (the pre-6.0
  combined status; `GSR` was missing), and the report's own tooltips (raw p).

### Changed — input and representation (#122, #123, #124, #125, #126, #147, #208, #149, #218)

Measured on the MSK sign-out dump (1,133,044 rows; operator decisions 2026-09-25
and 2026-10-05).
- **Non-sequence MAF alleles are FAIL rows (#123).** An allele that is not a base
  sequence (an IUPAC code such as `R`, `.`, a stray character from a hand edit)
  made the row count 0 silently; it is now a FAIL row, `NON_SEQUENCE_ALLELE`,
  kept in MAF output. Lowercase bases are bases; `-` is a MAF dash allele, so a
  `-` given as non-MAF input (the observations API) is one too. VCF output,
  which cannot carry such an allele (nor an empty one), writes the symbolic
  record `<NON_SEQUENCE>` (REF the reference base at POS, declared in the
  header), in counting runs and in `gbcms convert`. 7 such rows in the sign-out
  data.
- **A `REF_MISMATCH` row says where the given REF sits (#218).** The row stays
  FAIL and uncounted; `gbcms_diagnostic` (VCF `GD`, and `gbcms normalize`'s
  new `gbcms_diagnostic` column) gives `REF_AT_OFFSET(k)` when the REF (3+
  bases) matches the reference exactly within 3 bases of its position, every
  such offset listed nearest first. In the sign-out data 87 of 154
  `REF_MISMATCH` rows sit 1–3 bases off, mostly legacy ANNOVAR-annotated
  indels at Start−1 whose alleles still carry the VCF anchor base. No tool
  surveyed moves or explains such a row (`bcftools norm --check-ref` exits,
  warns, excludes or fixes REF in place; maf2vcf skips it).
- **VCF output of MAF input names its MAF row.** Every record carries
  `MAF_START`, `MAF_REF` and `MAF_ALT`: the row's Start and alleles as
  written (a placeholder such as `0` too; percent-encoded, `.` when empty), so
  a result can be looked up by its input, as VCF input's MAF output carries
  `vcf_pos`, `vcf_ref` and `vcf_alt`. Also in `gbcms convert`.
- **A MAF deletion at Start 1 is counted (#122).** It has no base before it; it
  is resolved to the VCF spec's position-1 form (the base after it), as gbcms's
  VCF output already wrote it, and counts as the same event given as VCF does
  (it was `FETCH_FAILED`). None in the sign-out data.
- **`End_Position` is optional (#124).** gbcms places a variant by
  `Start_Position` and its alleles; rows without an integer `End_Position` were
  skipped and are now read (maf2vcf converts them). `gbcms merge` no longer
  joins on it: inputs that write it differently for one variant join into one
  row (the first input's `End_Position`, the difference logged), and a row only
  a later input has keeps that input's. Where an input lists one variant twice
  with different `End_Position`, it joins on `End_Position` too, so each row
  pairs with its own counterpart. The sign-out data has it, consistent
  with Start and REF, on every row but one; no two of its rows differ only in
  `End_Position`.
- **VCF input's MAF output fills `Tumor_Seq_Allele1` (#125)** with the reference
  allele: MSK's sign-out convention on every row, and maf2vcf's reading of an
  empty one. Previously empty.
- **One row, one allele (#126):** a `Tumor_Seq_Allele1` that differs from both
  REF and Allele2 is not a second allele (vcf2maf's reading; cBioPortal picks
  Allele1). Documented and tested; no such rows in the sign-out data.
- **Engine API (#147):** a `decomposed` or `sibling_variants` list shorter than
  the variants is padded; a longer one raises `ValueError` (a short decomposed
  list panicked; a long list was cut silently).
- **One allele-kind rule (#208):** the AD-claiming guard's pure-indel test and
  the pure-indel windows use `allele_kind` (case-insensitive anchor; a
  multi-base shared prefix is complex, as the dispatcher counts it). No prepared
  sign-out row changes; it reaches only lowercase or unprepared engine input.
- `is_indel` in preparation is `ref_len != alt_len` (#149; no change).


### Changed — RNA: strand at intronic loci, splices are not coverage, clips and junctions (#185, #198, #186, #173, #213)

Five RNA rules (operator, 2026-10-04), each measured on FORTE against the
previous build: the truth cohort (94 signed-out rows), the T9 exon-edge indel
probes (978) and the C1 splice probes (7,224 delins at exon edges and mid-exon,
where no read carries the ALT).
- **Gene strand at intronic positions (#185).** A position no stranded exon
  covers (intronic, splice sites included) takes the strand of the transcripts
  spanning it, all agreeing; where both strands' genes cover a position there is
  no strand and every read counts (both genes' transcripts carry the allele).
  The exon index read the first and last intron bases as exonic (an end-inclusive
  interval tree built from half-open exons): fixed, for the strand and the
  per-transcript counts. The unresolved-strand warning names the variants.
  Truth: no row changes; T9: 95 antisense REF reads move to
  `rna_antisense_depth` on 62 rows.
- **A splice is not reference coverage (#198, RJ-17).** A pure-indel read is
  informative only when one aligned block between splices spans the window, REF
  and ALT alike. Truth: no row changes; T9: 24,703 REF reads at splice-crossing
  deletions become depth only (their bases fit both alleles). `vaf` (alt over
  REF plus ALT) at such a deletion with carriers rises to the VAF among the
  reads that show the event; no measured row's `vaf` moved (the T9 rows have no
  carriers). A deletion written right after a read's splice (no aligned flank)
  is depth only.
- **Diagnostics read the counted reads (#186).** `OBSERVED_ALLELE` and
  `COEXISTING_ALLELE` take n/m over the reads the counts read: no antisense read
  under enforcement, and the RNA mapping rule's unique mappers. Counts unchanged.
- **An RNA read's clip at an exon edge is not allele evidence (#173, RJ-18).**
  STAR clips a junction overhang it cannot splice, so a clip reaching an exon
  edge or junction end (annotated, or one the reads splice at), or ending within
  five bases of one, holds the next exon's bases and is not read. Any other clip
  is the read's own bases and is read as in DNA (a mid-exon MNP carrier with its
  second base clipped still counts ALT). Splice probes: spurious ALT 19 → 16.
- **A spliced read is read across its junctions (#213, RJ-19).** The
  exact-carrier windows of a read spliced near the event are built over the
  reference spliced at its own junctions (the far exons read from the FASTA,
  every junction the windows reach followed), not cut at the exon edge; a
  junction starting inside the event splices the haplotypes at its edge, so a
  delins carrier counts however the gap is written; such a read counts ALT only
  when its bases also beat REF spliced at its own junction (an alternative
  donor or acceptor), and is otherwise depth only. Splice probes: spurious ALT
  19 → 4; REF −0.95% at probes 0–1 bp from the exon edge (reads reaching one or
  two bases past the junction hold no spliced flank), −0.11% at 2–4 bp.
- **All five together**, against the branch point: splice probes spurious ALT
  19 → 2, REF −0.39%, partial −453; T9 106 rows (REF −24,813, mostly R5); truth: 4 rows'
  per-transcript columns (the exon-index fix), no count; RC DNA and WES
  byte-identical (before the review follow-ups; RNA-only changes since).
- **Survey:** GATK splits RNA reads at N and counts a piece only when its bases
  favour an allele, and runs HaplotypeCaller with `-dont-use-soft-clipped-bases`;
  bcftools never uses spliced reads for indels; phASER, WASP and ASEReadCounter
  never read clips; REDItools resolves a site's strand from the annotation
  spanning it and leaves mixed strands undetermined; allele counters take every
  tally over one filtered read set.

### Documented — why the RNA mapping-quality default is `--min-mapq 1`

No behaviour change. The RNA default keeps reads STAR placed at two to four loci
(MAPQ 3 or 1), counted once at their primary alignment, because junction reads tie
between a gene and its processed pseudogene. Measured on FORTE: counting unique
alignments only would cost real ALT reads at genes with pseudogenes (PIK3CA E545K:
10 of 159; 12 across the truth set) and 1.6% of junction fragments at the probes;
`--min-mapq 0` would add 6 ALT and 33 REF reads across the truth set. The STAR
MAPQ scale and the measurement are in `docs/reference/read-filters.md`, which also
corrects a note that STAR gives novel junctions low MAPQ (its MAPQ depends only on
the number of loci).

### Changed — an exact-carrier ALT call needs quality-weighted evidence (#174)

A read's bases across a complex variant's windows are now weighed by their
quality (a base matches with 1 − e and mismatches with e/3, e its error
probability), and an ALT call needs them to favour ALT over REF by at least what
one base read at `--min-baseq` gives (about 2.5 log10 at the default 20). Bases
matching both alleles cancel: the evidence is the weight of the bases that
mismatch REF less those that mismatch ALT, including a long event's bases past
its junction windows.
Clearly read bases must still match the ALT, as before. A low-quality base now
counts for little instead of fitting either allele.
- **Why:** reads whose window was mostly low quality could be called ALT on one
  clear sequencing error. Synthetic probes where no read carries the ALT found
  2.0 spurious ALT reads per million in FORTE RNA, 1.8 in IMPACT, one in WES and
  none in ACCESS duplex or simplex.
- **Measured:** on the probes, RNA spurious ALT 70 → 15–19, IMPACT 7 → 1, WES
  1 → 0; on every complex DNA/WES row with ALT reads, 1–2 of 1,185 real ALT
  reads lost. Two simpler rules were measured and set aside: requiring every
  event base read (46 real reads lost) or at most one masked (15 lost).
- **Not in this change:** spliced reads whose junction an aligner placed a few
  bases late at an exon edge (about 9 of the remaining RNA probe calls), and a
  REF molecule with one clear error just outside the window its ALT reading is
  anchored away from; 6.7.0.

### Changed — reads deleting the anchor, and the ALT written across several ops, are judged by their bases (#202, #201)

- **A read whose own deletion covers a pure indel's anchor** (#202) is judged by
  its bases between its nearest aligned flanks: ALT when they equal the ALT
  (masked bases fit, at least one base read), with the aligner's placement of
  its gap only a tie-break; REF when they equal the REF (the reference written
  as the anchor deleted and re-inserted); otherwise neither, with partial
  evidence. Before, Phase 3 credited such a read to whichever of REF and ALT was
  closer, so reads holding another allele counted ALT or REF; at insertion rows
  they counted REF.
- **The ALT written across several insertion or deletion ops** (#201), such as a
  deletion written as two or a 1bp deletion in a run written as D2 + I1, counts
  ALT when the read's bases spell it. Before, it counted partial.
- **Measured:** develop vs this branch on RC DNA, FORTE RNA and WES, every MAF
  cell compared.
  - The split-op rule changed no row (144 of 144 files byte-identical).
  - The anchor rule changes 33 rows: ALT net 0 on RC DNA (±1 read at 6 rows),
    −2 on WES; REF +57 and partial +210 at the BRCA2 cluster, where reads that
    were a co-annotated deletion's false ALT now count REF where their bases
    across the row's window are REF.
  - Against the read census, the changed indel rows move toward it: summed
    distance REF 344 → 251, ALT 40 → 44 (masked-base edge cases).

### Changed — what a read contributes: adapter read-through, absent qualities, unmapped records (#176, #182, #183, #207)

The read-judgment spec gains RJ-10 to RJ-12, which decide what a read brings
before any rule reads it.
- **A read ends at its fragment end** (#176). When the insert is shorter than the
  read, the bases past the mate's 5' end are adapter. They are now hard-clipped as
  the read enters counting, as if the read had been trimmed, so they are neither
  bases nor reach in any rule.
  - Only adapter-like bases are clipped: soft-clipped, or at most two aligned past
    the boundary (an aligner's chance extension into adapter), none inserted.
  - A read whose bases go on aligning past the boundary, or hold an insertion
    there, keeps them. TLEN is a reference distance, so it leaves out an ITD's
    inserted bases and a mate's clipped 5' bases.
  - Only an inward-facing pair defines a fragment; an outward one (TLEN negative
    on the forward read) is left alone.
  - In MSK data the first base past the boundary is A (the adapter's first base)
    in 97% of clipped reads.
- **Records without base qualities** (QUAL `*`) are dropped by the read filter
  and warned once per BAM (#182). Before, they voted as Q255, and the fragment
  consensus margin could overflow; it now saturates.
- **Unmapped records** (flag 0x4) are dropped by the read filter (#183). An
  unmapped mate placed at a variant counted in `mq0_count`, and at
  `--min-mapq 0` as a read (an ALT read when it carried the ALT). A read whose
  mate is unmapped still counts, and mapped MAPQ-0 alignments stay countable at
  `--min-mapq 0` (pseudogene loci such as PMS2).
- **Fixed:** the complex classifier read a hard-clipped read's anchor quality
  from the wrong base (#207).
- **Measured:** develop vs this branch on RC DNA, FORTE RNA and WES, every MAF
  cell compared, at the default MAPQ and at `--min-mapq 0`.
  - 182 rows change at the default MAPQ. RC DNA: 160 rows, ALT −57, REF −66. WES:
    13 rows, REF −12. FORTE: 9 probe rows, REF −3; the truth set is unchanged.
  - Every changed row is explained by reads that read through. All 57 lost ALT
    reads had their ALT on a removed adapter overhang; 52 of the 54 at SNVs
    showed the adapter's A.
  - On SNV rows the change equals the read census's (ALT −53, REF −31). On the
    57 changed pure-indel rows the engine moves toward the census: summed
    distance REF 167 to 125, ALT 69 to 62.
  - At `--min-mapq 0`: 184 rows. The two extra rows lose one REF read each to
    the same clip on reads below MAPQ 20. The PMS2 row keeps its MAPQ-0 reads
    (REF 5,439 at MAPQ 0 vs 5,414 at the default). No row moved by an unmapped
    record.
  - The FLT3 ITD a first build lost 13 ALT reads at is unchanged.

### Changed — which reads count REF for a pure indel; long complex events read through the read (#200, #199)

Read judgment now has one spec, `docs/reference/read-judgment.md`:
a table of read shapes and their calls, with the decision behind each, executed
by `tests/test_read_judgment_spec.py`. Changes to how reads are judged are
decided against it.
- **Another indel inside the discrimination window means not REF** (#200). A
  read carrying another insertion or deletion inside the window (not the ALT
  at another placement) is neither, with partial evidence: its bases are not
  REF there. This extends the sibling REF guard (#119) from co-annotated rows
  to any indel. Most visible at deep slippage loci: a BRCA2 cluster where an
  unannotated 1bp deletion in an A run lies inside two rows' windows, and
  homopolymer runs.
- **Another indel outside the window is a separate event** (#200). Of any
  length, it leaves the read REF where its bases across the window are REF,
  with no partial evidence. Before, a 5bp-or-longer one withdrew REF in a
  repeat, and kept REF with partial evidence in unique sequence.
- **Long complex events read through the read** (#199). The junction windows
  of an event too long for one read are read on inward as far as the read
  reaches. A read whose later bases contradict an allele is no longer called
  that allele: for example, a substitution-only read reaching the end of the
  run after a C>TA or CA>T, or a read that keeps the anchor and changes the
  run length.
- **Measured:** develop vs this branch, on RC DNA, FORTE RNA and WES, with every
  MAF cell compared.
  - 112 of 144 files are byte-identical. 155 rows change, all from the
    inside-window rule: REF −4,327, partial +4,324. ALT and depth are
    unchanged everywhere. The long-event change moved no row.
  - Checked against the read census (bases only): on the 105 changed DNA/WES
    rows it counts 57,199 REF reads, develop 60,619, this branch 57,218. On
    every FORTE row, the reads moved out of REF are no more than those whose
    bases contradict both alleles.
  - The largest moves:
    - a BRCA2 cluster: the +AAG row goes 1,413 → 408 REF in one sample
      (census 368), the 12bp-deletion row 1,213 → 402;
    - a 1bp-deletion row in a homopolymer: 1,009 → 698;
    - a FORTE T-run probe: about −1.5% REF.

### Changed — code-quality sweep: logging, monitoring, dead code, duplication (#204)

Counts are unchanged on prepared input: 144 of 144 acceptance files (RC DNA,
FORTE RNA, WES) are byte-identical to the previous build, every MAF cell compared.
- **Logging.**
  - Per-read WARN and DEBUG lines move to trace.
  - Rows the engine can only judge degraded are warned once per counting pass,
    with their count: indels and complex variants without a prepared reference
    context, MNPs whose REF equals ALT, and rows with an empty allele.
  - The exact-carrier fallback warning names its reason.
  - Misleading texts are corrected; for example, an exact-length deletion is
    no longer called "wrong-length".
  - Per-read traces name their read; the Phase stats line is 1-based.
  - Every decision path has a named trace: why an ALT call is kept, and each
    route to the exact-carrier rule.
- **Monitoring** (no new columns).
  - Each read's `read call` trace names the rule that decided it (`rule=`).
  - The Phase stats line counts these rules per variant: reads withdrawn as
    uninformative, ALT kept by its bases, exact-carrier judged or fell back,
    sibling-guard exclusions and clip admissions.
  - One INFO line per counting pass gives the totals.
  - Prep counts what it could not fetch in full (capped or missing shift
    regions, missing or short event references, missing reference contexts,
    capped left-alignments) and warns once when any occurred.
- **API.**
  - `prepare_variants` takes `rescue_homopolymer` (default off). It builds the
    homopolymer twin only when the twin will be counted; the pipeline and
    `observe_molecules` pass the flag.
  - `BaseCounts.singleton_alt_count` and `duplex_alt_count` are removed. Nothing
    ever wrote them; they were always 0.
  - `Variant.shift_region` and `event_ref` are read-only from Python.
  - A row with an empty allele counts neither and is warned once per pass; no
    read can carry an empty allele. Prep fails such rows, so pipeline counts are
    unchanged. Two places did see them:
    - `observe_molecules`, which passes FAIL rows through to keep rows
      positional: such a row's molecules are now OTHER (with an empty ALT they
      counted REF);
    - unprepared rows passed straight to the engine, which could reach the MNP
      check, where `len - 1` underflowed.
- **Dead code and duplication.**
  - Unused parameters, branches and setters are removed.
  - Observations and mFSD share one molecule classifier.
  - Six reference-end CIGAR walks, the soft-clip rule, the scan pad and read
    window, the read loops' quality, scoring and molecule-key code, and the
    insertion and deletion checks' end-of-walk resolution each have one helper.
  - The large-deletion band's literals are named constants.
  - One allele-kind classification (SNV, MNP, insertion, deletion, complex)
    drives the dispatcher, splice triage's span, the exact-carrier rule's scope
    and prep's variant-type label, which each spelled out the same rule.
- **Comments.** Stale comments are corrected, and ticket labels are removed
  from code comments and log text.

### Changed — `observe_molecules` requires a reference FASTA (#204) — breaking

- `observe_molecules()` takes `reference_fasta` from the argument or from `config`.
  Neither (or an empty string) raises `ValueError`, as the CLI does; so does a path
  that is not a file. Without a reference, the variants were neither normalized nor given a reference context.
  The indel and complex-variant rules then ran degraded, with nothing said:
  - carriers written elsewhere in a repeat were not recognised;
  - reads were judged against the bare event;
  - complex variants went to the previous classifier.
- `ObservationResult.variant_status` is always set (a list; it was `None` without a
  reference).
- The CLI's config refuses a reference that is not a file. An empty `--fasta ""`
  is `.`, a directory, which passed the existence check and failed in the engine
  with a FASTA-index error.

### Fixed — missing data never changes a count silently (#204)

- **ALT calls without a reference to read against.** An ALT call on a pure indel
  from a read that spans neither window stood on placement alone when the
  reference around the event was unavailable. At a contig end the ALT side's
  reading stretch failed to fetch, so carriers whose bases fit both alleles
  counted ALT while the same REF reads were withdrawn. Now:
  - the rule reads the reference that exists, stopping at the contig end;
  - a read with a side to read from but no prepared reference (an unprepared
    variant, or one that failed prep) counts toward depth only, and is counted
    and warned once per variant;
  - a read lying inside the tract has nothing to read from either way, so it is
    depth only;
  - a read whose deciding base lies past a contig edge stays depth only (a
    circular contig continues at its start), but is counted and warned with the
    others.
- **Reads starting on a flank.** The ALT read-by-bases rule reads from a flank base
  the read starts on (rightwards), or ends on (leftwards). The ALT side's windows
  already started there, and the flank is a base both alleles share. Carriers that
  start on the anchor and whose bases discriminate count ALT; before, they were
  withdrawn.
  - The read must read that flank base: unmasked (BQ at the threshold) and the
    reference's. For a deletion sliding through a repeat, the flank is all that
    tells a carrier from REF placed one base along the run. This now holds for the
    ALT windows too: a window that starts (or ends) on the read's own end base
    needs that base read, where before the read's extent alone made it
    informative. A carrier with a masked flank is no longer ALT.
- **The exact-carrier rule's reference.** Prep widens the reference for exactly the
  variants the rule judges: one predicate, `carrier::judges`, shared with the
  dispatcher. Anchor-changing insertions (C>TA) and the homopolymer twin now get
  it, so a long run no longer sends them silently to the previous classifier. A
  read the rule still cannot judge is counted and warned once per variant.
  - In a long run, the rule's junction windows have a known gap (#199): a read
    holding only the left junction can decide the call. Both sides are strict
    xfails.
- **Splice triage for anchor-changing insertions.** An insertion that also changes
  its anchor (A>CCC) uses the delins span. An RNA read spliced over the anchor
  and resuming at the next base leaves depth; it counted toward DP as neither.
- **Soft-masked references.** Every reference window prep reads is upper case.
  Over a soft-masked (lower-case) FASTA region, left-alignment could not shift.
  The SW and PairHMM scorers saw reference bases as mismatches: on a synthetic
  deletion locus REF went 6 → 12 under SW, and partial 6 → 0 under both.
  Prepared alleles there are now upper case.
- **mFSD determinism.** The fragment-size statistics are bit-identical run to run.
  The size vectors arrived in hash order, and the float sums depended on it: 17
  distinct `mfsd_ref_llr` values over 20 identical runs. The `--mfsd-parquet` size
  arrays are now in a fixed order. The MAF and VCF columns round far coarser and
  are unchanged.
- Three older tests expected ALT from 10bp all-A reads inside an all-A context.
  Their bases fit both alleles, so they are depth only.
- **Measured** on RC DNA, FORTE RNA and WES (develop vs branch, every MAF cell
  compared):
  - 140 of 144 files are byte-identical; the other fixes target cases this data
    does not hold;
  - four rows change, each adjudicated read by read:
    - a 66bp tandem duplication gains 11 ALT reads (8 fragments). Each starts on
      the anchor and holds the whole 66bp insert past the first base where the
      alleles differ, so they are carriers by their bases;
    - three deletions (14bp, and 1bp in two G runs) lose 5 ALT reads in all; no
      fragment counts change. Each read starts or ends on a flank base read at
      BQ 9–15, so its bases fit both alleles.

### Changed — binning invariance and a read census replace the legacy parity path (#170, #171)

- The per-variant `count_bam` (with `count_single_variant`, about 700 lines that
  duplicated the binned loop) is removed, along with the `legacy-parity` Cargo
  feature. Production never called it, so output is unchanged: the RC, FORTE and
  WES acceptance runs are byte-identical.
  - It shared the classifier, so it could not see classification bugs. Every one
    found this cycle came from judging reads by their bases.
  - Its one bin bug, the anchor's fetch end, was found in review. Synthetic
    fixtures fit in one bin, so parity could not have caught it.
  - It cost every counting change a second edit.
- `count_bam_binned` and `count_bam_binned_observations` take two optional test
  arguments, `bin_window` and `bin_max_variants` (default: the production
  constants; below 1: `ValueError`). Counts must not depend on them:
  - A Rust property test checks that every variant lands in exactly one bin whose
    fetch holds its read window.
  - `tests/test_binning_invariance.py` compares every field under the per-variant
    fetch (window 1, cap 1), tiny windows, small caps, one call per row, shuffled
    input on 4 threads, and decoy rows. Each geometry is asserted to split the bins.
    Its fixtures are synthetic DNA (plain, and BAQ+UMI+mFSD; siblings, a decomposed
    twin, overlapping mates), clip carriers past a long anchor, an RNA locus with a
    GTF and antisense reads, and the real test BAM.
  - Every counting test now runs through `count_checked`, which also counts with
    one variant per bin, so the 49 tests that ran only the legacy path test
    production. Most of those calls hold one variant and so check the per-variant
    fetch window; the binning suite covers multi-variant bins.
- `tests/census.py` judges each read by its own bases across a pure indel's tract,
  with the decided rules: REF needs one base past the first difference, ALT only
  that base. The census finds each tract by its own slide, and a test checks that
  prep agrees. `tests/test_read_census.py` checks the engine against it on reads
  generated around ten pure indels (ending anywhere past the anchor, some
  soft-clipped). The open decisions are strict xfails: #200, #201, #202.
- Both new suites were mutation-checked. Re-introducing the bin-end bug fails
  the property test and the clip-carrier fixture. Disabling the ALT-side window,
  or the equivalent-placement rule, fails the census test.

### Fixed — pure-indel reads count ALT only where their own bases hold the ALT (#188, #191, #192, #121)

- C10 (#160) takes REF only from reads whose bases settle the allele. ALT reads
  had no such check, so a read the CIGAR calls ALT that ends inside a repeat,
  where its bases fit both alleles, counted ALT. An ALT read now stands on its
  CIGAR when it spans C10's windows as seen from the ALT haplotype. For a
  deletion these are an insertion's windows, because a carrier's reference extent
  counts its own gap: a −AA carrier in a run of six A's could span the REF
  windows while holding only `G AAAA`. Otherwise the read keeps ALT only when it
  reads, unmasked, a base where the alleles differ, with every unmasked base it
  reads fitting the ALT. Failing both, it counts as depth only, as on the REF
  side. The check reads bases, not reference coordinates: 2,738 RC carriers,
  mostly of long insertions, span neither reference window, but their bases
  discriminate, and they keep ALT. There is no margin base. C10's margin guards a
  REF call made from the CIGAR alone, whereas this check reads the deciding base
  itself.
- An indel written at the junction counts ALT only when the read has no other
  insertion or deletion across the variant's discrimination window. A read with
  +A and −T for a +A row, or a split +AA, is another allele (neither + partial in
  a repeat). At 50bp or more the large-deletion band's length tolerance applies,
  and the band now counts changes across the whole window. Before, it stopped
  short of the window for a deletion sliding through a long repeat. The order an
  aligner writes an I/D pair in no longer matters. Before, a deletion written
  right after an insertion (`M I D M`), or an insertion right after a deletion
  (`M D I M`), was not inspected, and the read counted REF. A same-length
  insertion of other bases beside another indel in the window is another allele.
  Phase 3 called such reads ALT, although their length change is never the ALT's.
  A deleted anchor followed by an insertion is judged by the read's bases.
- A same-length deletion placed elsewhere that gives another haplotype now
  counts ALT when its bases spell the ALT across the window. Compensating
  mismatches can make the read the ALT allele even though its gap sits elsewhere.
  The read is anchored on its nearest aligned bases outside the window. Wherever
  the haplotype its own CIGAR proposes differs from the ALT, the read base must be
  unmasked and match the ALT: an N where two deletion alleles differ does not
  count. Otherwise, at 5bp or more, it is a distinct allele. These reads went to
  Phase 3, which called them ALT; on the RC set every read Phase 3 reached was
  another allele by its bases. Under 5bp the read stays REF as before. A deletion
  of another length within the ≥50bp band keeps Phase-3 arbitration (#191).
- A distinct allele counts neither + partial where the event slides (its shift
  region is wider than the event) or where the read has an indel inside the
  discrimination window. It counts REF + partial only in unique sequence with the
  window clear. Before, only a repeat of a ≤6bp motif gave neither + partial, so
  carriers of a long-period duplication with another allele counted REF (#192).
- A one-base-REF variant whose ALT changes the anchor base (`A>CCC`) went to the
  insertion check, which matched the inserted bases and ignored the substituted
  anchor. A read keeping the REF anchor and carrying only the insertion counted
  ALT. Such variants now go to the exact-carrier rule (#141): a carrier holds the
  whole ALT, anchor included, across a window that grows through repeats and is
  padded to the ALT's length. Reads that do not hold that window are depth only,
  so a read's call does not depend on where it ends (#121).
- The Phase-3 context for a slid indel is not changed (#159). Sizing it by the
  shift region changed 0 ALT calls on the RC, WES and RNA arms, so it was closed
  with the measurement.
- Measured on the RC panels (realigned), WES without realignment (80 loci from
  paired IMPACT/TEMPO libraries) and FORTE RNA; develop vs branch, every changed
  read adjudicated by its own bases:
  - DNA: 15 of 1,060 rows change (ALT −39, REF −107, partial −15, depth
    unchanged).
    - ALT: carriers that end inside the repeat or before the deciding base (7
      at a 1bp deletion, 11 at a 66bp duplication ending with the insert, 5
      at two 14bp deletions). Also reads carrying other deletions (24–54bp, or
      a 33bp one 15 bases off) that Phase 3 had called ALT for a 33bp
      deletion, and 6 reads at a 6bp deletion whose bases fit both alleles.
    - REF: −118 at the two anchor-changing rows (#121), from reads that do not
      carry the REF anchor or do not hold the exact-carrier window.
    - At a +AAG insertion co-annotated with that 33bp deletion (two samples),
      16 reads move from partial to REF. They carry a large deletion across the
      anchor, or a D1 inside the window, and only the 33bp row's false ALT kept
      them out of REF. The insertion row's own rules count them (#200, #202).
    - At a 6bp deletion, 5 reads whose bases are REF across the window lose REF
      to partial. They carry a same-length deletion of other bases outside the
      window (#200).
  - WES: 10 of 80 rows change, all anchor-changing (#121): ALT −12, REF −399,
    partial +152, depth +62. At one locus 146 reads keep the REF anchor and
    carry only the insertion, so they count partial, not ALT. At others, reads
    that carry the whole ALT gain ALT.
  - FORTE RNA: the truth set and the STAR repeat-insertion probes are unchanged.
    One 4bp deletion probe loses 2 ALT reads (a 3bp deletion ending just after).
- Three adversarial reviews; every changed read was adjudicated. Follow-ups:
  C26 #200 (which of a read's other indels decide its REF call), C28 #202 (a
  read deleting the anchor falls back to Phase 3's closer haplotype), C27 #201,
  C25 #199 and R5 #198.

### Fixed — carriers of an indel written elsewhere in its repeat count ALT, not REF (#189)

- An indel in a repeat can be written at any junction of its shift region with the
  same haplotype. The input is left-aligned, but aligners may not be: near the
  FORTE truth set's indel rows STAR placed 162 of 170 repeat insertions away from
  the left-aligned position (the DNA panels, realigned: 37 of 15,047; WES without
  realignment: 0 of 620). The windowed checks accepted
  a shifted placement only by a proxy:
  - an insertion when the reference base before it equalled the anchor base. That
    is never true inside the repeat, so carriers written elsewhere counted **REF**
    (in a synthetic `G AAAAA T` locus, REF 10 / ALT 0 instead of 5 / 5), and it was
    true by chance for some placements that give another haplotype, which counted
    ALT;
  - a deletion when the bases it removes equalled the given ones, so a rotated STR
    placement (removing `AC` for `CA`) under 5bp counted REF and the same bases
    deleted outside the repeat counted ALT;
  - the scan reached `max(5, repeat_span + 2)` bases from the anchor, and
    `repeat_span` counts motifs of up to 6 bases only, so a longer duplication's
    carriers written past that reach counted REF.
- A placement is now accepted when it gives the variant's haplotype (the same
  bases, or a rotation, elsewhere in the shift region) and is the read's only
  change across the variant's discrimination window. The scan reaches every
  placement: every junction of an insertion, a deletion's starts up to its last
  placement. Other placements:
  - the variant written elsewhere with another gap, insertion or splice across the
    window (a +AA read for a +A row, a split −4 read for a −2 row, a deletion
    cancelled by an insertion): a distinct allele (neither + partial in a
    repeat). A read spliced over the anchor with a deletion written just after
    the junction shows only where the aligner ended the splice: depth only (it
    counted ALT);
  - the variant's inserted bases where they give another haplotype, inside the
    variant's discrimination window: a distinct allele, as a wrong-length
    insertion is (neither + partial in a repeat, REF + partial in unique
    sequence). Outside the window the read shows the window as reference: REF, as
    before;
  - the variant's bases inserted before its anchor base were counted ALT by the
    backward-boundary check; they are now another haplotype (REF). An
    anchor-substituting ALT (`A>CCC`) is never matched by a placement elsewhere
    (C8, #121, covers the strict path);
  - a same-length deletion that gives another haplotype: REF under 5bp, as before;
    at 5bp or more it still goes to Phase 3 (C22, #191, revisits that route).
- Measured on the RC set, develop vs branch, every changed read adjudicated by its
  own bases:
  - DNA: 13 of 1,060 rows change (ALT +46, partial −62, REF +13, depth
    unchanged). Two insertion rows gain 44 and 10 ALT from carriers written up to
    6 junctions right, each holding the ALT haplotype (their last inserted base is
    an N). Two 2bp deletion rows move 6 and 18 reads from partial to REF: the
    reads delete the same two bases 6bp before the anchor, another haplotype
    outside the window. Their co-annotated 33bp rows regain the same reads as REF,
    which develop had claimed for the 2bp row. Seven 1bp indel rows lose 8 ALT
    reads that delete or insert the same base outside the run (one ends before
    the anchor).
  - FORTE RNA: the 33 truth samples and the T9 probes are unchanged (the truth
    set's 3 insertion rows have no carriers). At 23 probes built from STAR's
    repeat insertions near the truth loci, ALT goes from 11 to 183; an independent
    census counts 176 shifted equivalent carriers, all ALT on the branch.
- This also changes grouped rows (#99): a tract-mate deleting the same bases at
  another position is another haplotype, not this row's allele, so where its
  change lies outside this row's window its carriers are REF here instead of being
  matched and then demoted to partial. Tract-cluster grouping now reaches as far
  as the scan.

### Fixed — indels near a contig end are left-aligned and get their reference windows (#142)

- Prep pads its reference windows on both sides but clamped them only at the
  contig start. The FASTA reader rejects a window that passes the contig end, so
  near a contig end each window failed within its own reach:
  - left-alignment within about 100bp (a WARN; the variant kept its input
    position, so reads aligned at the left-aligned position could fall out of
    depth);
  - `ref_context` within its padding (5–50bp), so there were no Phase-3
    haplotypes;
  - the shift region within 256bp (up to 16kb in long repeats), and the event
    reference within 60bp, for every variant type;
  - a complex variant fell back to the tolerant classifier that the exact-carrier
    rule (#141) replaced, so the same reads were judged differently near a contig
    end. In a synthetic case the same bases and reads gave REF 0 / ALT 0
    mid-contig but REF 5 / ALT 7 at the contig end, two reads carrying another
    allele among the ALT.
- Every window (left-align, `ref_context`, the adaptive repeat scan, the shift
  region, the event reference) now holds the part of the contig it covers. SNVs and
  MNPs near a contig end therefore also get their event reference, used by the
  observed-allele diagnostic. The contig is resolved by the same name the fetch
  reads, so a FASTA holding a contig under two names cannot mix their lengths.
  Exact fetches (REF validation, a MAF anchor) stay exact, so a REF running past
  the end still fails validation instead of being "corrected" to a shorter one.

### Changed — `rna_antisense_depth` counts antisense reads at RNA defaults (#114)

- With strandedness enforced (the RNA default), antisense reads were dropped
  before the sense/antisense tally, so `rna_antisense_depth` (VCF `ANT`) was
  always 0. Such a read is now classified as a sense read would be, and tallied
  when it is a first-class REF or ALT read over the anchor. It is then dropped, so
  REF, ALT, depth and every other count are unchanged. The column means the same
  with `--no-strandedness`.
- `STRAND_DISCORDANT` is documented as a `--no-strandedness` diagnostic: under
  enforcement, antisense reads never reach the junction tally wherever the gene
  strand is resolved.
- Docs: the read-filters strand table is corrected for the default `reverse`
  protocol, and the gene strand is documented as coming from the `--gtf` exons (not
  the MAF).
- The per-read trace line names each excluded antisense read
  (`antisense_excluded=true`).

### Fixed — records stored without bases (SEQ `*`) are not counted (#172)

- A BAM record with no sequence, such as a secondary alignment kept with
  `--no-filter-secondary` or a primary in a stripped BAM, is now dropped where
  every read enters counting. It counts in neither depth, fragments nor
  `mq0_count`. Each counting pass logs one WARNING when any were skipped.
  - Before: the SNV check panicked on such a record, and so did heuristic BAQ (RNA,
    or DNA with `--apply-baq`). The insertion and deletion checks counted it as REF
    and depth (or as a fragment) from its CIGAR alone.
  - Counts change only where such records reach counting. With the default
    filters, that means primary records without bases.

### Changed — the exon-edge BAQ rule and `exon_boundary_dist` measure from the REF span (#106)

- `exon_boundary_dist` (VCF `EBD`) is now the distance from the nearest annotated
  exon boundary to any base of the variant's REF span, and `0` when a boundary lies
  inside the span (as in Ensembl VEP's overlap test). It was measured from the first
  base only. The span is the normalized (left-aligned, VCF-style) REF: SNVs and
  insertions (a one-base REF) are unchanged, and a pure deletion's span starts at
  its anchor base (VEP drops the anchor, so its distance to a boundary on the left
  is one more). This is a column change for multi-base variants: the distance is
  smaller for about 4 in 10 signed-out multi-base variants and indels.
- RNA mode with `--gtf`: the BAQ exception at exon edges uses the same distance, so
  a multi-base variant reaching into an exon's last five bases skips BAQ, as an SNV
  there does. The main counts, per-transcript counts, ASJD and each `--rescue-mnp`
  component re-count resolve the rule the same way: a component is counted over
  its MNP's span, and a rescued row reports the MNP's distance (it reported the
  adopted component's).
  - Why: heuristic BAQ lowers the five bases before every splice junction, below
    `--min-baseq` at Q30 or Q37. An MNP starting more than 5bp from a right exon
    edge but reaching into those bases lost its edge-side bases in every spliced
    read. So no spliced carrier showed the whole MNP (`--rescue-mnp` then adopted a
    component where the reads show the annotated haplotype), and a spliced read
    carrying only the edge-side change voted REF on its other bases.
  - Counts change only for multi-base variants that reach within 5bp of an edge
    while their first base does not. DNA is unaffected (no GTF).

### Changed — a row counts the given allele; the homopolymer twin is opt-in (#163)

- The homopolymer twin is no longer dual-counted by default. It was the corrected
  allele gbcms counted for a delins like `CCCCCC>T` (`CCCCCT`), reporting whichever
  form had more ALT reads. The row now counts the allele it is given.
  `--rescue-homopolymer` (Nextflow `--rescue_homopolymer`) restores the dual count,
  flagged `WARN_HOMOPOLYMER_DECOMP` as before.
  - Why: both forms accept near-matches, so the winner could report another
    allele's reads under the row's label. At the case the twin was built for
    (SOX2), the reads carry `CCCCT`, not the twin. The twin won by tolerance.
  - It won at 2 of 11 real twin loci, and on 0 of the 6.5.0 RC rows, so counts
    change only where it used to win.
- `observe_molecules()` takes `rescue_homopolymer` and threads the twin the same
  way.

### Changed — complex variants count exact carriers only (#141)

- A delins, a deletion whose anchor also changes (Del+SNV), or an MNP read with
  an indel now counts a read as REF or ALT only when the read's own bases carry
  that allele across the whole event, with two reference flank bases each side.
  Soft clips count; bases below `--min-baseq` match anything, the only tolerance.
  - Why: local alignment, likelihoods and the edit-distance margin credited
    near-matches, so reads carrying another allele, or ending inside the event,
    were counted for the given one.
  - The windows are read at the event's own position, grown through repeats
    touching the event on either allele, and padded to equal length so neither
    allele is favoured by where reads start. Events over 50 bases are judged at
    both junctions, reading a short ALT whole; a mismatch at either rules a read
    out.
  - A read that cannot hold both windows is depth only (no allele, no
    `partial_alt`, no mFSD class). A read that holds them, matches neither and is
    closer to ALT counts in `partial_alt`.
  - An MNP read with an indel in or right beside the block is judged the same way
    (an aligner may write a shifted block as an insertion and a deletion); an
    indel further off leaves the base-by-base comparison in place.
  - Prep fetches enough reference to hold an event grown through a long repeat.
  - RNA: a spliced read's windows end at its own junction, so a delins at an
    exon's first or last bases counts spliced reads; a read spliced through
    the event counts toward depth only.
- Complex variants (DNA) count reads whose allele lies in soft-clipped bases. An
  aligner clips an ALT read near its end while the REF reads beside it align in
  full, so those carriers used to fall outside depth and VAF read low. A read the
  exact-carrier rule decides from its own bases, inside a well-defined fragment
  (past the fragment end a clip is adapter), now counts in DP and REF/ALT;
  undecided clipped reads stay out. Not in RNA mode (#167).
- BAQ (the RNA default, and DNA with `--apply-baq`) no longer penalizes a read's
  own insertion or deletion when that indel is the variant being counted. At Q37
  the penalty put every ALT read of a small indel or delins below min BQ while
  REF reads kept full quality; FORTE's Q40 bins hid it. Splice junctions and
  indels elsewhere in the read are penalized as before (#166).
- In a multi-allelic group, a read that matches this row's ALT and a sibling's
  ALT exactly (the two differ only at bases the read has masked) is neither
  row's AD: it counts in `partial_alt`.
- `--alignment-backend pairhmm` and `sw` now give identical counts for these
  variants: the rule uses no alignment scoring.

### Added

- Two `gbcms_diagnostic` flags name the allele the reads carry when it is not the
  given one (canonical VCF form; n reads carry it exactly, m the given allele):
  - `OBSERVED_ALLELE(chrom:pos:REF>ALT:n/0)`: no read carries the given allele
    exactly, so the input is likely mis-described.
  - `COEXISTING_ALLELE(chrom:pos:REF>ALT:n/m)`: the given allele is present, but a
    different allele in the same stretch is more frequent (e.g. a germline indel
    or stutter in a repeat). A caveat for reading the VAF.
  - Both need n ≥ 3, n > m, at least 5% of the scanned reads, and an allele that
    is not already an input row.
  - The scan covers every variant type: each spanning read is compared over the
    event plus one base each side, so a germline SNP beside the event does not
    count.
  - Counts are unchanged; the flag says what the reads show so the input can be
    checked.
  - Prep stores the event's reference bases on `Variant.event_ref`.

### Fixed

- **Indel REF counts use only reads that can tell the alleles apart** (#157).
  A read that starts or ends inside an indel's repeat tract matches both
  alleles: the aligner places no gap either way. It was counted REF from anchor
  coverage alone, which biased repeat-indel VAF down (a 50% homopolymer deletion
  read 25% in a synthetic test; on the 6.5.0 RC data the median gbcms-to-
  informative VAF ratio was 0.91 for STRs and 0.92 for homopolymers). Such reads
  now count toward depth only, as in GATK's AD. A deletion longer than a read
  still gets REF reads from either junction. A tandem duplication (an ITD)
  slides over its whole duplicated segment, so REF needs a read across all of
  it. Prep measures that region over its own fetch (`Variant.shift_region`);
  the repeat context kept for alignment is often too short (a 30bp duplication
  had an 11-base context).
  See "Informative Reads for Indels" in the allele-classification reference.
- **Grouped rows: REF reads and REF fragments exclude the same molecules** (#119).
  The sibling-ALT guard dropped a read from `ref_count` when it carried any
  co-annotated sibling's ALT, even one far from this row, and only after its
  fragment had been recorded as REF. The guard now applies only to siblings
  whose change lies inside the row's discrimination window. It runs before
  fragment evidence in every counting path, so `ref_count` and
  `ref_count_fragment` agree.
  - Carriers of a sibling outside the window count REF again: IGV shows them as
    REF at this row.
  - Carriers of a sibling inside the window now also leave `ref_count_fragment`.
- What moves with these fixes:
  - `ref_count`, `ref_count_fragment` and their strand forms, VAF, strand
    bias, and VCF `AD`/`ADF`/`ADR`/`FAD`;
  - values built from REF + ALT: the merged `simplex_duplex_total_count`, and
    RNA `rna_sense_depth`/`rna_antisense_depth` (these count REF and ALT reads;
    the docs now say so);
  - per-transcript REF counts and the ASJD REF partition;
  - `partial_alt`: sibling carriers outside the window count REF, not partial.
    A carrier of an in-window sibling still counts partial, even when it ends
    inside the tract.
  - mFSD: a fragment whose reads all end inside the tract lands in no class. It
    used to land in REF, and must not look like a third allele (NonREF). In the
    observation export it is `OTHER`.
- `alt_count_fragment` can rise by a molecule or two: a mate whose REF call was
  vacuous (it ended in the tract, or carried a sibling) no longer contests the
  other mate's ALT. This happened on 3 of 1,060 rows in the 6.5.0 RC DNA runs,
  by +1 to +2.
- Unchanged: `alt_count`, `total_count` (DP) and `total_count_fragment`. SNV and
  MNP rows change only when they share a site with a co-annotated indel. None
  did in the RC runs.

## [6.5.0] - 2026-09-25

### ⚠️ Breaking Changes — VCF ↔ MAF representation follows vcf2maf / maf2vcf (#110)

> VCF-input MAF output changes coordinates, alleles and `Variant_Type` for many
> indel and MNP shapes and gains two columns; MAF-input VCF output changes
> `POS`/`REF`/`ALT` for indels. Counts are unchanged for every countable
> allele shape; only the edge rows below (Allele1 fallback, `ALT_EQUALS_REF`,
> skipped ALTs) change. Re-genotype every flavor
> (e.g. duplex and simplex) with the same version before `gbcms merge`.

- VCF input → MAF writes each record exactly as vcf2maf does. It trims the
  leading bases REF and ALT share, never trailing ones; `Start_Position`,
  `End_Position` and `Variant_Type` follow the trimmed alleles. Before, the
  conversion followed a length-only type label, so 9 of 17 allele shapes
  differed from vcf2maf. For example:
  - `TTAC>A` became `702–704 TAC>-` (vcf2maf: `701–704 TTAC>A`, DEL);
  - `C>TA` became `->A` (vcf2maf: `821 C>TA`, INS);
  - `TCT>TCG` stayed a 3bp TNP (vcf2maf: SNP `T>G` at the changed base);
  - `TC>TCGG` stayed untrimmed (vcf2maf: `- > GG`).

  Bases are compared case-insensitively, so only mixed-case alleles can differ
  from vcf2maf.
- New MAF columns `vcf_ref` / `vcf_alt` keep the VCF record itself, alongside
  `vcf_pos` (vcf2maf's names). Each row's `vcf_alt` is its own allele;
  multi-allelic records give one row per ALT. These columns appear only for VCF
  input.
- MAF input → VCF output writes each row as maf2vcf does. When an allele is `-`,
  or the alleles differ in length and first base, the reference base before
  them is prepended from `--fasta`. At position 1, the base after the event is
  appended instead (VCF spec). Before, `-` alleles were written into the VCF
  (`462 AA>-`), which is not valid VCF.
- MAF alleles are read as maf2vcf reads them:
  - `Tumor_Seq_Allele1` is the variant allele when `Tumor_Seq_Allele2` is empty
    or the reference (older MAFs). Such rows were genotyped as REF against
    itself; the reader now logs how many rows it read this way.
  - Alleles made only of `-`, `?` or `0` are read as `-`.
- A variant whose ALT equals its REF (any case; `-` for both in a MAF) is a
  `FAIL` row with the new reason `ALT_EQUALS_REF`. Before, it passed and was
  counted, meaninglessly.
- `--show-normalization` `norm_*` MAF columns follow the left-aligned alleles,
  not a type label (a delins was written with its first base stripped).
- Type labels come from the alleles in both readers and in variant preparation.
  INSERTION / DELETION now requires a shared anchor base; a delins such as
  `TTAC>A` or `C>TA` is COMPLEX. `gbcms normalize` reports these labels.
  Counting never read them.
- The VCF reader skips alleles it cannot count, with a WARNING and per-reason
  totals: ALT `*`, symbolic `<...>`, breakends, a missing `.`, other
  non-sequence alleles, and a REF that is not a base sequence. Before, `*` and
  `<DEL>` were genotyped as nonsense rows, `.` was dropped silently, and an
  empty REF (read as `.`) was genotyped. ALT `N` still reaches preparation,
  which reports it as a `FAIL` row.
- `gbcms merge` adds the VCF record (`vcf_pos` / `vcf_ref` / `vcf_alt`) to its
  join key when every input is VCF-derived, because two VCF records can trim
  to one MAF record. It also warns when an input repeats a join key; each such
  row joined every matching row of the other inputs, silently.
- The mFSD HTML report finds each variant's MAF row by the key its Parquet uses
  (the VCF record for VCF input). VCF-input indels, and some other shapes, were
  reported without their MAF row (no gene, no statistics).
- Counting is otherwise unchanged: every countable allele shape counts
  identically from VCF and MAF input.
- Docs: for VCF input, `Strand`, `Variant_Classification` and
  `Tumor_Seq_Allele1` are empty, and `vcf_id` is empty for a `.` ID. They had
  been documented otherwise. `End_Position` is listed as a required MAF column.

### Added — `gbcms convert`

- `gbcms convert` converts VCF → MAF (vcf2maf's coordinates, with `vcf_pos` /
  `vcf_ref` / `vcf_alt`) or MAF → VCF (maf2vcf's records; needs an indexed
  `--fasta`) without counting. It uses the same conversion as the `dna` / `rna`
  output.

### Fixed — homopolymer-decomposed twin: strand and per-sample flag (#107)

- `WARN_HOMOPOLYMER_DECOMP` is per sample again. The prepared variants are
  shared by every sample of a run, so once one sample's corrected allele won,
  every *later* sample's row carried the flag, whatever its reads showed.
  Output depended on sample order.
- In RNA with a GTF, the corrected allele takes its original's gene strand, so
  `--enforce-strandedness` applies to it. Before, where it won, its counts
  included antisense reads (measured: 0.065% of its depth at synthetic
  probes).
- The docs describe the corrected allele as built: the run with its last base
  replaced, the same length as REF (`CCCCCC→CCCCCT`). They had described a 1bp
  deletion plus the change (`CCCCCT` was written `CCCCT`). The SOX2 example's
  chromosome is corrected to 3.

### Fixed — output keeps the input's contig naming (#103)

- `chr`-named input (e.g. hg38 with UCSC names) no longer loses its naming in the
  output. VCF input → MAF writes `Chromosome`/`vcf_region` as the input names them
  (was the stripped `1`), and VCF output writes `CHROM` in the input's naming with
  `##contig` lines declared under the same names (reference lengths kept) — before,
  records read `1` under a `##contig=<ID=chr1>` header, a malformed VCF that htslib
  rejects (`Contig '1' is not defined in the header`). Rescue labels and the mFSD
  Parquet follow the same naming. Counting is unchanged: contigs are still reconciled
  internally between the variant file, FASTA and BAM. One INFO line per distinct
  naming pair names the two conventions when they differ. Unprefixed (b37,
  Ensembl) input against a reference and BAM that use the same naming is
  unaffected.
- The reference FASTA now reconciles the mitochondrion's spellings as the BAM side
  already did. Before, a `chrM` (or `M`, `MT`, `chrMT`) variant against a FASTA
  naming it differently was rejected `FAIL` / `FETCH_FAILED` with zero counts.
- `gbcms merge` joins inputs whose `Chromosome` naming differs (`chr1` ~ `1`,
  `chrM` ~ `MT`) into one row instead of two half-empty ones. It writes each
  contig one way, the first input's name where it has the contig, and logs each
  later input's difference.
  - An input that names one contig two ways now gets a WARN, since a variant
    listed under both names is repeated.
  - Input columns that collide with merge's join helpers are rejected with the
    column named.

### Fixed — `gbcms merge` row order is deterministic

- Merged rows came out in a different order on every run: the full outer join
  guarantees no order (5 runs on the same real inputs gave 5 orders), so merged
  MAFs could not be diffed across runs. Rows now follow the inputs: the first
  input's rows as it lists them, then rows only a later input has, in that
  input's order. Row contents are unchanged.

### Fixed — one BAQ rule across RNA views; deterministic ASJD junctions

- Per-transcript counts and ASJD now apply the main counts' exon-boundary BAQ
  exception: BAQ is skipped at variants within 5bp of an annotated exon boundary,
  where its CIGAR-N penalty lands on exactly the reads that splice there.
  Before, both applied BAQ regardless. Measured on three junction-rich RNA
  samples at 2010 exon-edge probes:
  - per-transcript REF counts were ~3% low overall, and some transcripts lost
    all their spliced REF reads. The median gap to the main counts drops from
    2.5% to 0.05%, and the per-transcript mismatch rate now equals the main
    counts';
  - ASJD saw no junction at some edge variants, including known junctions
    carried by thousands of fragments.

  Main counts are unchanged. Away from exon edges, per-transcript counts are
  unchanged. ASJD rows change there only where the tie rule below applies, or
  through `asjd_qval`, which is corrected across the sample.
- ASJD's dominant junction no longer depends on hash order. On a tie, the same
  input could report a significant divergence on one run and none on the next.
  A tie is now not a divergence: when a REF top junction and an ALT top junction
  are the same splice event (both ends within 5bp, ASJD's tolerance
  throughout), each allele reports its own of the pair and no test is run.
  Otherwise each allele reports the tied junction the other allele's fragments
  use most, and the leftmost only on a further tie. Coordinates alone could
  otherwise decide the Fisher cells and `asjd_flag`.
- `exon_boundary_dist` is found in any contig naming. A `chrM` variant against
  an `MT`-named GTF read `2147483647`, and so did any contig the GTF does not
  annotate. The first now reads the real distance, so the exon-edge rule
  applies there. The second is now empty. The column is documented as the
  unsigned distance it has always been; it was previously documented as
  signed.

### Added — observability for silent fallbacks (no count changes)

- `SW_FALLBACK(n)` in `gbcms_diagnostic`, plus one WARN per affected variant
  naming the reason and outcome: under the default `pairhmm` backend, n
  depth-contributing reads could not be evaluated by the pangenomic haplotype
  matrix because the variant's reference context is missing (prep's fetch
  failed, e.g. an indel near a contig end) or does not contain it, or the ALT
  haplotype exceeds 400bp (very long insertions). They are scored by the
  Smith-Waterman fallback
  where it can run and otherwise end NEITHER — previously with no trace on the
  row. It never fired on the traced real runs. Never set under
  `--alignment-backend sw`.
- `CLIP_CANDIDATES(n)` in `gbcms_diagnostic`: an insertion with no confirmed
  ALT where n (≥ 2) reads carry a ≥ 8bp soft clip within the insert's
  duplication reach — clip-represented carriers (typically a tandem-duplication
  ITD) the engine cannot claim.
- `--umi-tag TAG` that no processed read carries now logs one WARN per counting
  pass over a BAM (a `--rescue-mnp` recount can repeat it): fragment grouping
  silently fell back to read names.

### Changed — Smith-Waterman gap-extend is a documented constant

- `dynamic_sw_gap_extend` rounded to −1 for every repeat span, so its repeat
  relaxation never engaged. It is replaced by `SW_GAP_OPEN = -5` /
  `SW_GAP_EXTEND = -1` constants; alignment scores are unchanged.

### Added — ASJD-2 splice-disruption markers (RNA + GTF)

- `asjd_diagnostic` gains two markers that read the population the
  splice-aware evidence rule excludes — fragments whose CIGAR `N` spans the
  variant — at variants within two bases of an annotated intron boundary
  on the gene's strand, at or above ASJD's own junction-evidence floors:
  `RETENTION_DOMINANT(n)` (spliced-over fragments dominate a junction-free
  classified population, so `vaf` is the retention-population VAF) and
  `NOVEL_JUNC_AT_SPLICE_LOSS(n@start-end)` (an anchored, unannotated junction
  on those fragments outnumbers confirmed ALT — the mutant allele's exon-skip
  or alternative-site outcome). Both are population comparisons with no tuned
  rates; they previously surfaced only as `LOW_REF_JUNC;LOW_ALT_JUNC`
  (issue #97). No new columns.
- Docs now state the `asjd_*_junction` coordinate convention (0-based,
  half-open intron).

### ⚠️ Changed — exclusive assignment at co-annotated tract clusters (user-visible)

- **Tract-cluster grouping.** `assign_multi_allelic_groups` gains a second
  membership criterion: length-changing variants (indels/delins — never
  SNV/MNP) whose scan windows (`max(5, repeat_span+2)` each side) overlap
  join a group transitively, even when their REF spans never touch. Members
  whose spans truly intersect keep the `MULTI_ALLELIC` reason tag;
  window-only members get the new `TRACT_CLUSTER` tag.
- **AD-claiming contest.** At a grouped locus, anchor-exact (Phase 0)
  evidence is never contested; every other ALT read is demoted to
  `partial_alt`/`any_alt` (excluded from AD *and* ADF) when any of three
  tests fires: the read's window does not favor the row's ALT haplotype
  strictly over its REF haplotype (foreign flank events); the call is
  alignment-phase on an anchor-preserved pure indel (no matching
  structural op anywhere — unannotated ladder absorption in tracts;
  complex/MNP rows exempt); or a co-annotated sibling's ALT haplotype
  explains the read strictly better over a window covering both spans
  (equal cost = equivalent representations, both rows keep the read).
  True carriers — pure indels, delins/complex, shifted self
  representations, noisy but real reads — keep AD, except a pure-indel
  carrier resolved only by alignment (e.g. soft-clipped), which surfaces as
  `partial_alt`.
- **REF-side symmetry.** Sibling-claimed reads are excluded from RD on every
  path (now including per-transcript); in the main counts they surface as
  `partial_alt` instead of vanishing silently.
- Validated against a signed-out ACCESS hypermutation cluster: per-row
  fragment ALT counts land exactly on (or within counting-basis of)
  sign-out where the previous windowed counting over-attributed 2–5×;
  panels without co-annotated clusters are byte-identical.

### ⚠️ Changed — MNP rescue reports a coherent component genotype (opt-in `--rescue-mnp`)

- **Candidate gate is partial dominance.** Rescue now considers PASS
  `MNP_RESCUE_ELIGIBLE` MNPs with `partial_alt > alt_count` instead of only
  `alt_count == 0`. Masked per-position evaluation counts a component carrier
  whose other discriminating base is low-BQ as full ALT, so a single such
  read blocked rescue of MNPs whose carriers hold only one component (a
  TERT promoter GAGGG>AAGGA row sat at alt 1 / partial 88).
- **Rescued rows adopt the component's full counts.** Previously only
  `alt_count` (and read VAF) changed, leaving strand, fragment, strand-bias,
  mFSD and diagnostic columns at MNP values and breaking
  `alt_count = alt_count_forward + alt_count_reverse` in the output. A rescued
  row now carries the winning component SNV's counts in every column (all
  counting invariants hold, including `any_alt = ad + partial_alt`), and
  `gbcms_diagnostic` is recomputed from them. Fragment-level consumers
  (ACCESS, `gbcms merge --add-combined`) now see rescued counts.
- **`gbcms_rescue` format** (downstream parsers): every entry has
  `outcome=` (`rescued`, `skipped_grouped`, `haplotype_confirmed`, `no_improvement`,
  `ref_validation_failed`) and the MNP's own `original_ref`/`original_alt`/
  `original_partial` (`original_alt` was hard-coded `0`); rescued rows add
  `adopted=`. A position whose synthetic SNV failed preparation reads
  `ref_fail` instead of a silent `0`. `outcome=no_signal` is gone.
- **MNPs whose haplotype the BAM shows are kept** (`outcome=haplotype_confirmed`,
  no re-count): the engine now counts MNP ALT reads in which every
  discriminating base was read (internal `BaseCounts.mnp_confirmed_alt`, not an
  output column), and rescue fires only when that count is **zero** — no read
  in the BAM shows the annotated MNP. Validation with matched normals found
  rescue adopting germline het SNPs where reads carried the real MNP (IMPACT:
  36 such reads; ACCESS duplex: 6); such rows now keep the MNP's counts. The
  audit gains `original_confirmed`.
- **Rescue labels use the output's contig naming** (an input MAF's own
  `Chromosome`, e.g. `chr1`, rather than the stripped internal name).
- **`gbcms merge` warns on mixed rescue:** when duplex and simplex rescue
  outcomes differ (one rescued, or different components), each row is named in
  a WARNING with a per-run count — with or without `--add-combined`, since the
  two flavors' columns then describe different alleles in one row (with it, the
  `simplex_duplex_*` columns add them); counts are unchanged. When only one
  flavor was genotyped with `--rescue-mnp`, merge says so once and names every
  row rescued in that flavor (previously that case was silently skipped).
- **Grouped MNPs are skipped** (`outcome=skipped_grouped`) so rescue cannot
  hand back reads that exclusive assignment gave a co-annotated sibling.
- **Fixed: rescue audit leaked across samples.** In a multi-BAM run a later
  sample's row could show an earlier sample's `gbcms_rescue`.
- **Rescued rows say so:** `gbcms_diagnostic` (MAF column, VCF `GD`) gains
  `RESCUED_COMPONENT(chrom:pos:REF>ALT)` on every rescued row, each rescue
  logs a WARNING (the component and the MNP's own counts), and enabling
  `--rescue-mnp` logs a WARNING that rescued rows report a component.
- **Fixed: VCF `GR` was split by parsers.** The audit's positions list was
  comma-separated and VCF parsers (htslib/pysam) split a comma-bearing
  Number=1 INFO value, silently truncating the audit. Positions are now
  joined with `+`; the MAF column and VCF `GR` carry the same content (VCF
  writes `;` as `|`).
- Per-sample INFO outcome summary; warnings for anomalies (a component SNV
  failing preparation). With `--mfsd-parquet`, a rescued row's record keeps
  the MNP's coordinates but carries the adopted component's fragment sizes
  (logged per sample).
- `BaseCounts.with_ad()` removed from the Python bindings (no remaining
  callers).
- Default runs (`--rescue-mnp` off) are unchanged.

## [6.4.0] - 2026-09-22

> The changes below alter reported counts at wrong-length indel loci and at
> spliced RNA positions — the minor version bump this release delivers.

### ⚠️ Changed — span-aligned REF testimony at spliced deletion loci (user-visible)

- **Structural REF evidence survives fragment consensus at quality zero.**
  Span-aligned REF's carried quality is the first exon base after the
  junction — exactly where BAQ can zero it — and `FragmentEvidence` gated
  REF existence on quality > 0, silently dropping such fragments from
  `rdf` (the mirror of the structural-ALT divergence fixed in the
  splice-aware change below). Span-REF observations now carry a structural
  flag that keeps them alive in consensus; a structural ALT in the same
  fragment still wins. On the FORTE locus this recovers one real fragment
  (`rdf` 726→727 = `dpf`).

- **Junction reads that observe every deleted-span base now count REF** at a
  pure-deletion locus whose anchor base is spliced out (typical for deletion
  annotations left-aligned onto the last intronic base at an acceptor). The
  anchor base is not the discriminating fact for a deletion — the span is —
  and these reads demonstrate the span is present. Previously they counted
  DP as neither (conservative cluster-A behavior): at a validated FORTE
  acceptor locus, 823 of the 824 anchor-spliced junction reads convert
  (`rd` 24→847 of `dp` 848, fragment `rdf` 23→727 of `dpf` 727, DP
  unchanged); the one
  residual read's exon2 alignment covers only one of the two deleted-span
  bases (`95M92N1M`) and honestly stays neither. Guards: the FULL span
  must be aligned (partial coverage stays neither), and reads carrying a
  competing shifted/wrong-length indel candidate keep their existing
  arbitration paths. Read- and fragment-level, mirrored in both engine
  paths via the shared checker.

### ⚠️ Changed — consensus intron snipping of ref_context removed (user-visible)

- **`ref_context` is always genomic now.** The RNA-mode step that drained
  consensus introns from the variant's reference context before Phase 3 had
  no coordinate map: the context shrank while `ref_context_start` stayed
  genomic, so shifted-indel sequence verification, the large-deletion band's
  context guard, and haplotype offsets all read garbage right of a snipped
  intron — so exon-contained reads were scored against haplotypes whose
  offsets no longer matched their genomic coordinates, and pre-mRNA /
  intron-retention reads against a haplotype missing bases their sequence
  genuinely contains. Removing it
  fixes windowed shifted-deletion matching near acceptors (contract test
  `test_windowed_deletion_after_junction_in_repeat_rna`, committed red) and
  erases the binned-vs-legacy RNA divergence around it; junction reads
  needing alignment-based scoring stay conservatively `neither` under the
  splice-aware evidence rule. Re-validated on local (non-repo) clinical RNA
  data — the b37 dedup and FORTE hg38 smoke loci report unchanged counts. A
  coordinate-mapped spliced-haplotype rework remains tracked in issue #94,
  gated on real-data measurement; it must keep pre-mRNA / intron-retention
  reads genomically scored.
- **`mq0_count` now tallies before the RNA strandedness filter in the binned
  path**, matching the legacy path: an antisense MAPQ-0 read is still a
  physical read at the locus, and the two paths previously diverged on this
  diagnostic in stranded RNA mode.

### ⚠️ Changed — splice-aware evidence in RNA counting (user-visible)

- **Reads spliced over a variant no longer count depth or REF.** A CIGAR `N`
  spanning every discriminating position (a deletion's deleted span, an
  insertion's junction flanks, an SNV/MNP's REF bases) is asserted splicing —
  no observation — so the read classifies neither AND is excluded from
  DP/fragment depth, matching samtools pileup's zero coverage inside an N gap.
  Previously such reads counted definitive REF (deletions) or inflated DP with
  zero aligned bases (intronic positions): on a junction-dense real RNA locus,
  DP at an intronic SNV drops from ~826 span-overlapping reads to the 2 with
  aligned bases there. Excluded totals are logged per variant at debug level
  (`Phase stats … splice_skip_excluded=`).
- **The REF/ALT call no longer flips on the aligner's D-vs-N representation
  choice.** A `D(100)` at the expected span is deletion evidence (ALT); the
  same gap as `N(100)` is splicing (excluded).
- **Indels directly after a splice junction are now examined.** An indel op
  whose only neighbor is a splice `N` (`M-N-D-M` / `M-N-I-M` — an event at an
  exon boundary reached through the junction) was structurally invisible to the
  strict and windowed CIGAR scans; carriers now classify through the same
  anchor/windowed inspection as `M`-adjacent ops (evidence attributed to the
  nearest aligned/inserted base when the anchor itself is spliced out).
- **Phase 3 no longer scores across a splice.** Raw-window extraction refuses
  windows that an `N` overlaps, and `check_complex`'s reconstruction classifies
  such reads neither instead of stitching exon arms into a junction-chimeric
  sequence in which the missing intron reads as deletion evidence.
- **The exclusion is not silent.** Deletion-type loci where splice-skip
  exclusions exceed confirmed ALT are flagged `SPLICE_SKIP_DOMINANT(n)` in
  `gbcms_diagnostic` — RNA aligners write large deletions as splices (STAR:
  ≥ `alignIntronMin`, default 21bp), so `alt_count`=0 there may mean the
  carriers exist as junction reads. Per-variant totals also appear in the
  debug-level `Phase stats` line (`splice_skip_excluded=`).
- **Shifted representations across a junction stay in classification.** A
  read whose N covers the annotated span but which carries an I/D op inside
  the scan window (repeat-tract shift) is deferred to the windowed scans
  instead of being excluded.
- **Fragment consensus recognizes structural ALT independent of base
  quality.** BAQ can stack splice and indel penalties to quality 0 on the
  base carried with a junction-adjacent I/D op; `FragmentEvidence::resolve`
  previously treated qual-0 ALT evidence as no evidence, reporting `ad > 0`
  with `adf = 0` at the same locus.
- **`check_complex` refuses string comparison for structurally anomalous
  junction reads** (splice N inside the context window): with the indel
  shifted outside the variant span, the span reconstruction is clean REF
  sequence and Phase 2 would absorb an ALT carrier into `rd`.
- ASJD junction tallies at intronic splice-region variants now reflect only
  allele-informative reads (spliced reads with no observation at the locus
  no longer contribute fabricated REF junctions); expect `LOW_REF_JUNC`
  where the REF tally was previously fed by such reads.
- Contract battery: `tests/test_rna_splice_contract.py` (13 cases incl. a
  legacy-parity case with N-CIGAR reads and one strict xfail pinning the
  known D6 ref_context coordinate corruption for the cluster-B fix).
  DNA-mode classification is untouched (reads without `N` ops never enter
  the triage), and the legacy parity oracle mirrors every engine-loop
  change.

### ⚠️ Changed — orchestration fail-fast (user-visible)

- **`--bam-list` entries that do not exist now fail the run** (exit 1) unless
  `--lenient-bam` is given — matching the long-documented fail-fast promise. An
  unreadable list file is equally fatal: previously both cases silently ran a
  partial sample set and exited 0.
- **Duplicate sample names are now a hard error** (from `--bam`, `--bam-list`, or
  a mix): the later BAM silently replaced the earlier one, which was then never
  processed. Use `sample_id:path` (or two-column list entries) to disambiguate
  deliberate same-stem inputs.
- **Ragged MAF rows now raise in the batch readers** (`gbcms merge` inputs)
  instead of having their overflow fields silently truncated — the gbcms count
  columns are the trailing columns, exactly what truncation dropped.

### 🔧 Fixed

- **A rejected (FAIL) variant no longer crashes `--mfsd`, RNA `--gtf`, or
  `--mfsd-parquet` runs**: the zero-count stub now carries every column the
  writers read (mFSD q-value and nucleosomal fractions, RNA/GTF annotation and
  ASJD fields), and the mFSD Parquet excludes rejected variants (they have no
  fragment data) with the exclusion logged.
- **`gbcms merge` no longer coerces float-formatted count strings to 0**: a
  pandas/R round-trip renders integers as `12.0`, which the combined-column sums
  silently nulled to 0; genuinely non-numeric values (e.g. `NA`) are counted and
  warned about per column.
- **Re-genotyping warns about replaced columns**: an input MAF already carrying
  gbcms output column names (e.g. `ref_count` from a previous run) has those
  values refreshed — now with a warning naming every colliding column, and the
  writer's docstring no longer claims originals are never overwritten.
- **Failed-sample reports name the exception type and log the traceback**
  (`KeyError: some_key` instead of a bare `'some_key'`).

### ⚠️ Changed — wrong-length pure-indel evidence is a distinct allele (issue #91)

- **A read whose CIGAR proves a pure indel of a DIFFERENT length at the variant anchor now
  counts as `partial_alt`, never as REF or the queried ALT.** In repeat tracts, coexisting
  distinct-length indel populations are distinct slippage alleles; the previous behavior
  routed them to the length-blind Phase-3 haplotype window, which promoted them to ALT and
  inflated VAF several-fold at tract loci (a homopolymer locus with 7 true-ALT reads
  reported 259). `alt_count` now reflects exact support; the distinct-allele evidence
  surfaces in `partial_alt`/`any_alt`, and `PARTIAL_DOMINANT` flags loci where it dominates.
- **Placement-aware ≥50bp deletion band** replaces the reciprocal-overlap tier: a
  wrong-length D at the anchor is the annotated event only when the read deletes essentially
  the whole expected span (≤3 retained bases, ≤3 changed outside) — accepting breakpoint
  wobble and split `D+M+D` representations while rejecting displaced net-matches. The old
  ≥50% overlap rule accepted any deletion sharing half the annotated length (its sequence
  check compared reference against reference, identically true by construction).
- **Insertion truncation containment**: a shorter insert that is a ≥90%-identity slice of a
  non-low-complexity expected insert still counts as ALT (sequencer truncation smear of long
  insertions); both sequences must be non-low-complexity so repeat-tract slippage is never
  mistaken for truncation.
- **Same-length insertions with confidently mismatching bases** (≥`--min-baseq`) are a
  third allele → `partial_alt`; they were previously absorbed into `ref_count` with no
  signal. Unverifiable cases (all inserted bases low-quality) go to Phase-3 arbitration.
- **Windowed wrong-length ops** (deletions ≥5bp — 1–4bp windowed Ds remain
  CIGAR-definitive alignment noise; insertions at any size): repeat tract → `partial_alt`;
  unique context → REF with the stray op surfaced as `partial_alt`. Same-length S3-fail candidates keep the Phase-3
  left-alignment rescue (TP53-class), now via the split `has_shifted_same_length` flag.
- **Repeat scan anchors at the first changed base** (not the shared VCF anchor, which sits
  one base left of a left-aligned tract), and adaptive context padding covers the whole
  tract plus flank — Phase-3 haplotype windows in repeat regions are no longer too narrow
  to distinguish tract lengths.
- Delins/complex variants are unaffected: they route to `check_complex`, whose Phase-3
  realignment correctly resolves split and mismatch-absorbed representations of one event.

### 🔧 Fixed

- **`--trace` never emitted a single per-read Rust trace line**: pyo3-log's default filter
  capped forwarding at DEBUG, and the Python side configured a logger name (`gbcms_rs`)
  that pyo3-log never uses (`_rs.…`). Both fixed; `_rs.reset_log_caching()` exposed so
  enabling trace after Rust code has logged is not silently ignored.
- **Left-alignment failures are loud**: wide-window FASTA fetch failure warns with the
  error (reachable near contig ends) instead of silently skipping normalization; the
  2500bp expansion cap binding without convergence warns; non-UTF-8 alleles from a corrupt
  reference keep the variant fully unnormalized instead of emitting a shifted position with
  reverted alleles; a stray `println!` on ref_context fetch failure removed.
- **`--enforce-strandedness` warns when variants have no resolved gene strand** (no GTF,
  uncovered locus, contig mismatch) instead of silently not enforcing.
- **GTF annotation warns per variant chromosome with zero loaded exons**, on both the
  text-parse and the `--gtf-cache-dir` cache-hit paths (a warm cohort cache previously
  silenced the warning after the first sample).
- Per-read `debug!` diagnostics demoted to `trace!` (WFA router, marginalized PairHMM).

### ✨ Added

- Wrong-length contract battery (`tests/test_wrong_length_contract.py`, 20 tests through
  the binned↔legacy parity oracle) and an end-to-end `PARTIAL_DOMINANT` reporting-chain
  test driving the full CLI.
- Documentation: "Wrong-Length Pure Indels" rule reference in
  `docs/reference/allele-classification.md` and a worked case study (Case 4 in
  `docs/reference/complex-indels.md`);
  `partial_alt` semantics updated in output docs; fictional `RUST_LOG`/`GBCMS_LOG_LEVEL`
  controls removed from docs (logging is `--verbose`/`--trace`).

## [6.3.1] - 2026-08-25

### 🔧 Fixed

- **Nextflow: turning off a default-on read filter now actually reaches the engine.**
  `--filter_duplicates false`, `--filter_secondary false`, `--filter_supplementary false`,
  and `--filter_qc_failed false` were silent no-ops: the modules omitted the flag when
  false, and the CLI default (on) won. The DNA/RNA modules now always pass the explicit
  on/off form (`--filter-x` / `--no-filter-x`) for all six read filters.
- **Nextflow ≥26.04: boolean and numeric CLI params are coerced explicitly.** The strict
  parser hands CLI overrides to the script as Strings (`--flag false` arrives as the
  truthy String `"false"`; 25.x coerced it to `Boolean false`), which silently inverted
  every boolean param check. A shared `asBool()` helper
  (`nextflow/modules/local/utils/main.nf`) now guards every boolean param decision, and
  `process.resourceLimits` casts `max_cpus`/`max_memory`/`max_time`. Verified end-to-end
  on 25.10 and 26.04: identical generated commands both ways for all six filters.

### ✨ Added

- **Contract tests pinning filter and phasing behavior** (no production code changes):
  - `test_qc_failed_filter_contract` — three-way qc-failed filter proof at read and
    fragment level: flagged reads dropped by default, `--no-filter-qc-failed` restores
    byte-identical counts, no effect on unflagged BAMs.
  - `test_partial_haplotype_contract` — the cis-phasing contract for DNP/ONP variants:
    a multi-nt ALT counts toward `alt_count` only when the full haplotype co-occurs on
    one read; partial-haplotype reads are surfaced via `partial_alt`/`any_alt`
    (`any_alt == ad + partial_alt` asserted end-to-end), and a spurious multimer with
    no read support reports all three as 0.

## [6.3.0] - 2026-08-07

### ✨ Added

- **`observe_molecules()` takes the settings it actually uses.** New keyword arguments
  `filters`, `quality`, `alignment`, `umi_tag`, `threads`, `apply_baq` and `library_type`
  sit alongside the existing `config=`, which keeps working; an individual argument
  overrides the matching `config` field.

  Previously the only way to change a setting was to build a whole `GbcmsDnaConfig`, whose
  four required fields — `variant_file`, `bam_files`, `reference_fasta`, `output` — have
  **zero overlap** with the seven this entry point reads. Two of them must name files that
  exist on disk. So adjusting one filter meant fabricating an output directory and a variant
  path that are never touched: paths that look load-bearing and are not.

  Found by the first real consumer (mulligan) needing a non-default read filter.

  `umi_tag` is the one argument where `None` is an explicit choice rather than an omission —
  it means "group by read pair, not UMI family" and overrides `config` accordingly. Every
  other argument treats `None` as "not supplied". Without that distinction, a caller passing
  `umi_tag=None` alongside a config that sets one would silently inherit it and change what
  counts as a molecule.

### 🐛 Fixed

- **`--filter-secondary` and `--filter-supplementary` now change the result when turned
  off.** They previously could not change **any** output: an unconditional skip in the
  counting loop discarded those records *before* fragment evidence, regardless of the
  flags. Measured across all six read filters on real MSK-ACCESS data, these were the only
  two that did nothing (`duplicates`, `improper_pair` and `indel` all worked).
  ([#82](https://github.com/msk-access/gbcms/issues/82))

  **Default runs are byte-identical** — both filters still default to on, so such records
  never reach the cache. Verified on real data: total depth unchanged at 1202 across 40
  loci, and binned↔legacy parity holds.

  Turning a filter **off** now admits those records to **fragment**-level evidence
  (`dpf`/`rdf`/`adf`) and to the per-molecule observation export. They still never
  contribute to read-level `dp`/`rd`/`ad`: a supplementary shares a QNAME with its primary
  and is the same physical read, so counting both would report depth 2 where one read
  exists. Fragments are immune either way — `hash_molecule` keys on QNAME, so a primary and
  its supplementary collapse into one.

  What this fixes concretely: a locus reached **only** by a supplementary segment used to
  report `dpf=0` — not a filtered read but a wrong answer, since a molecule demonstrably
  covers it. That is what made a molecule spanning a large deletion invisible to
  cross-locus phasing.

  New consequence to be aware of: with the filter off, an admitted supplementary raises
  `dpf` while leaving `dp` unchanged, so the two stop moving together. The observation
  export reconciles with `dpf`, and that invariant is verified under both settings.

## [6.2.0] - 2026-08-06

### ✨ Added

- **`min_mapq` on `Observation`** and in the observations Parquet schema — the worst MAPQ
  among the reads backing each molecule. Minimum rather than best, deliberately: a molecule
  is only as trustworthy as its least confidently placed read, so a paired fragment reports
  its worse mate and one cleanly-placed read cannot launder an ambiguous one.
  `255` is the SAM spec's *unavailable*, distinct from `0`, which is a real value meaning
  the read mapped ambiguously.

  Exported rather than left to inference because the `--min-mapq` threshold *bounds* mapping
  confidence without measuring it: assuming every surviving molecule sat at the gate would
  charge them all the same error, swamping the base-quality term and flattening exactly the
  per-molecule differences a weighted linkage statistic exists to express. Consumers doing
  read-backed phasing need the measured value.

  `FragmentEvidence` gained a matching `min_mapq` field, tracked for **every** read reaching
  fragment evidence — including reads that are neither REF nor ALT, which still count toward
  `dpf` and still describe the fragment's placement.

### 🔧 Changed

- **`alignment_backend` now defaults to `pairhmm` in `gbcms._rs`**, matching the CLI,
  `Pipeline`, and `observe_molecules`, which all already defaulted to it. A direct `_rs`
  caller that omitted the argument previously got Smith-Waterman — a different classifier
  than the identical call made through any supported entry point (measured: 5 of 104
  molecule calls differ on real ACCESS deletions). **No production path changes**: the CLI,
  Nextflow, and the Python API have always passed the backend explicitly, so counts are
  unchanged.

### 🐛 Fixed

- **An unrecognized `alignment_backend` token now raises `ValueError`** instead of silently
  falling back to Smith-Waterman. `"PairHMM"`, `"smith-waterman"`, or any typo previously
  produced plausible counts from a classifier the caller never asked for, with nothing
  logged. Validation now happens before any I/O, so a bad token fails immediately rather
  than after a full read pass. Unreachable from the CLI (typer validates the enum first);
  reachable from `_rs`. ([#79](https://github.com/msk-access/gbcms/issues/79))
- **`_rs.pyi` corrected to match the extension, and the match is now enforced.**
  `count_bam` and `count_bam_binned` each declared nine defaults (`min_mapq` … `threads`)
  that the pyo3 signature does not grant, so mypy accepted calls that fail at runtime with
  *missing 9 required positional arguments*; `prepare_variants` and `write_fsd_parquet`
  named parameters that do not exist at runtime (`reference_fasta` → `fasta_path`,
  `output_path` → `path`), so any keyword call using the documented name raised `TypeError`.
  New `tests/test_rs_stub_parity.py` compares the stub against pyo3's `__text_signature__`,
  turning AGENTS.md invariant 5 from a review item into a CI gate — it catches all four of
  the above.
- Removed the `Default` impl on `AlignmentBackend`, which still named `SmithWaterman` long
  after `pairhmm` became the real default everywhere else. Nothing called it; it served only
  to mislead.

## [6.1.0] - 2026-08-06

### ✨ Added

- **`--observations-parquet` CLI flag** on `gbcms dna` and `gbcms rna`, writing
  `<sample>.observations.parquet` alongside the normal counts output. Follows the
  `--mfsd-parquet` pattern (flag → `OutputConfig` → pipeline calls the native writer). This is
  what makes the capability reachable for Nextflow / ACCESS-Pipeline / CWL consumers, which
  invoke gbcms through the CLI rather than the Python API. Counts output is unchanged, and
  nothing is written unless the flag is passed.
- **Parquet sink for observations.** `observations_path=` on
  `count_bam_binned_observations()` / `gbcms.observe_molecules()` writes the rows to
  Parquet **from Rust**, so they never cross the FFI boundary. A panel- or genome-wide run
  produces 10^6–10^7 observations, where materializing that many Python objects — not the
  counting — is the bottleneck; the same reasoning behind `write_fsd_parquet`. The file is
  self-describing (`variant_index, chrom, pos, ref, alt, molecule_hash, allele, best_qual`),
  since a bare index means nothing once the data outlives the call that produced it.
  Verified on real data: identical row content to the in-memory path.
- **Per-molecule observation export.** `count_bam_binned_observations()` returns the same
  counts as `count_bam_binned()` plus the per-molecule allele calls the engine normally
  aggregates away — one `Observation` per fragment per variant, carrying the molecule
  identity, the resolved allele (REF / ALT / N / OTHER), and its best supporting base
  quality. Because a molecule observed at two variants in one call shares a
  `molecule_hash`, callers can link alleles **across loci** (read-backed phasing, allelic
  imbalance, duplex-masking QC); gbcms itself does no such linking.
  - **Public API:** `gbcms.observe_molecules()` → `ObservationResult`. This is the
    supported surface; `gbcms._rs` remains internal and is not covered by SemVer.
  - **Counting is untouched** — both entry points share one core, `count_bam_binned`'s
    signature and return are unchanged, and observations are off unless requested
    (no allocation when off). Binned↔legacy parity is unaffected.
  - Rows are sorted by `(variant_index, molecule_hash)`, so output is deterministic
    despite parallel bin processing and hash-map iteration order.
  - Decomposed variants emit **one** set of rows: observations follow the same
    higher-`ad` arbitration as the counts, so the losing allele form is never exported.

### 🔄 Changed

- **Nextflow pipeline is now compatible with the strict syntax parser** that is the
  default from **Nextflow 26.04**, while still running on 25.x. Verified with
  `nextflow lint` on 26.04.0 (0 errors) and config resolution on both 25.10 and 26.04.
  - `nextflow.config`: replaced the legacy `check_max()` helper (function definitions
    are illegal in strict config) with the native `process.resourceLimits` directive;
    `--max_cpus/--max_memory/--max_time` still apply. Bumped the `nextflowVersion`
    floor to `>=24.04.0` (when `resourceLimits` landed).
  - `main.nf`: moved input validation into the `workflow` body (top-level statements
    are disallowed) using `error` instead of the removed `exit`; rewrote the
    `hasData` helper without a `while` loop (removed in strict syntax); `Channel` →
    `channel`.
  - **Removed the remote nf-core institutional-config include.** The strict config
    syntax forbids the `try/catch` (and `if`) guards that made a fetch failure
    non-fatal. The local `iris`/`slurm` profiles are self-sufficient; users who want
    an institutional config can still supply one with `-c`. The `custom_config_base`
    / `custom_config_version` params were removed.

## [6.0.0] - 2026-07-01

> [!WARNING]
> **Breaking output changes — downstream parsers of gbcms output must review before upgrading.**
> This major release changes several output behaviours:
> - **`gbcms_status` is now two fields** — a verdict (`gbcms_status` = `PASS`/`FAIL`) and a
>   `|`-separated `gbcms_status_reason`. VCF gains a `GSR` INFO key. Filters on
>   `FAIL_FETCH_FAILED` etc. must switch to `gbcms_status == "FAIL"` + a reason check.
> - **Per-transcript count columns use `|` between transcripts, not `;`** (RNA + `--gtf`).
> - **Supplementary/secondary alignments no longer count toward read-level depth** — `DP`
>   may drop for BAMs with those alignments and the corresponding filters disabled.
> - **Empty-allele variants are rejected at prep time** (verdict `FAIL`, reason `EMPTY_ALLELE`).
> - **Nextflow RNA now defaults `min_mapq=1`** (via `rna_min_mapq`) instead of inheriting the
>   DNA `20`, matching the `gbcms rna` CLI — RNA depths increase vs. the 5.x pipeline.
> - **`--mfsd` MAF/VCF schema grew** to 41 MAF columns / 13 VCF `MFSD_*` INFO fields
>   (gated behind `--mfsd`, absent otherwise).

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

- **Nextflow RNA runs no longer silently drop STAR multi-mapper reads.** The pipeline used a
  single global `min_mapq = 20` for both modes, so RNA counting ran at MAPQ 20 instead of the
  `gbcms rna` CLI default of **1** — silently dropping STAR's 2–4-locus multi-mapper primaries
  (which STAR encodes as MAPQ 3/1, vs 255 for uniquely-mapped reads). On a real FORTE RNA
  sample that was ~16% of reads genome-wide (chr7) and up to ~6% of depth at an individual
  locus. A new `rna_min_mapq` param (default **1**) now drives the RNA module; DNA keeps
  `min_mapq = 20`. Set `--rna_min_mapq 20` to restore unique-only RNA counting.

### 🔄 Changed

- **`gbcms_status` is split into a verdict + a reason field (breaking output change).**
  The status is now two fields with a consistent grammar in both formats: a verdict
  (`gbcms_status` = exactly `PASS` or `FAIL`) and a `|`-separated reason list
  (`gbcms_status_reason`, empty for a clean PASS). MAF gains the `gbcms_status_reason`
  column; VCF gains a `GSR` INFO key. Reason tags dropped their verdict prefix
  (`FAIL_REF_MISMATCH` → verdict `FAIL` + reason `REF_MISMATCH`; `PASS;WARN_REF_CORRECTED`
  → verdict `PASS` + reason `WARN_REF_CORRECTED`). Reasons now **stack** with `|`
  (`WARN_REF_CORRECTED|WARN_HOMOPOLYMER_DECOMP`), which also fixes a bug where a
  REF-correction warning was silently dropped on the success path and where
  `WARN_HOMOPOLYMER_DECOMP` overwrote any existing reason. `|` is used (never `;`/`,`,
  both VCF-INFO-unsafe), so the reason string is byte-identical in the MAF column and the
  VCF `GSR` — no format-specific conversion. **Consumers filtering on `FAIL_FETCH_FAILED`
  etc. must switch to `gbcms_status == "FAIL"` + a `gbcms_status_reason` check.**

- **Fragment-count consensus is labeled accurately (ME-12).** The stale "Majority Rule"
  comment on the `BaseCounts` fragment fields (`DPF`/`RDF`/`ADF`) is replaced with
  "quality-weighted consensus with an INDEL structural-priority override and a discard
  band," pointing at `FragmentEvidence::resolve` — which is what the code has always done.
  No behavior change. The alternative of gating a structural-ALT INDEL win on the REF mate's
  base quality was considered and rejected (anchor BQ is orthogonal to INDEL-detection
  confidence; see `REJECTED.md`).

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

### 📖 Documentation

- **Repo-wide documentation accuracy sweep.** A four-way audit against the current code
  corrected the docs where they had drifted: the mFSD schema (now **41 MAF columns / 13 VCF
  `MFSD_*` INFO fields**, up from a stale 34/7, with the q-value + nucleosomal-fraction +
  CH-flag columns added); the mono-nucleosomal range (150–200 bp); `gbcms_status` values (the
  real `FAIL_*` prefixes + the two missing statuses); crashing examples (`normalize` uses
  `--output <file.tsv>`, not `--output-dir`; corrected RNA MAF column names in the quickstart);
  the Nextflow README defaults (on-by-default filters, `pairhmm` backend, Nextflow ≥22.10.1) and
  its missing param surface (`--strandedness`, `--gtf`, `--library_type`, `--mfsd*`, `--merge_*`);
  and previously-undocumented features (`build-gtf-cache` / `--gtf-cache-dir`, exit-code
  semantics, the `--threads` hard budget, and contig auto-reconciliation). Also fixed the
  `--rna-editing-db` CLI help (it loads a REDIportal **TABLE1**, not a VCF) and several
  developer-doc references to files/paths that had moved.

### 🧹 Internal

- **Single type stub for the Rust extension (LO-1).** Removed the orphan
  `src/gbcms_rs.pyi`, which stubbed a top-level `gbcms_rs` module that isn't importable
  (the extension is `gbcms._rs`) and had already drifted from the real one (missing the
  nucleosomal fields and `fisher_exact_2x2`). The co-located `src/gbcms/_rs.pyi` (shipped
  via `py.typed`) is now the single source of truth, and its `count_bam` stub gained the
  missing `reference_fasta` parameter so it matches the real `#[pyo3(signature)]`.

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
