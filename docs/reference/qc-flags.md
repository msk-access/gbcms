# QC Flags

Every flag, reason and class gbcms writes, in one place, with where it appears
and what to do about it. The tables are grouped by family; each column's full
definition stays in [Output Formats](output-formats.md). A test checks that every
flag string the code can emit is on this page.

**Where flags appear.** MAF columns are listed first, then the VCF field. In VCF,
`;` inside a value is written `|`.

| Family | MAF column | VCF field | Mode |
|:-------|:-----------|:----------|:-----|
| [Verdict and status reasons](#verdict-and-status-reasons) | `gbcms_status`, `gbcms_status_reason` | INFO `GS`, `GSR` | all |
| [Diagnostics](#diagnostics) | `gbcms_diagnostic` | INFO `GD` | all (some RNA only) |
| [MNP rescue](#mnp-rescue) | `gbcms_rescue` | INFO `GR` | `--rescue-mnp` |
| [ASJD splice flags](#asjd-splice-flags) | `asjd_flag`, `asjd_diagnostic`, `asjd_*_motif` | not in VCF | RNA with `--gtf` |
| [QC columns](#qc-columns) | `rna_editing_site`, `mfsd_ks_valid`, `mfsd_ch_flag` | INFO `RED`; mFSD INFO | RNA; `--mfsd` |
| [mFSD report classes](#mfsd-report-classes) | — (HTML report) | — | `--mfsd-report` |
| [Record shapes in VCF output](#record-shapes-in-vcf-output) | — | ALT `<NON_SEQUENCE>`, INFO `MAF_*`, `NORM_*` | MAF input; `--show-normalization` |

## Verdict and status reasons

`gbcms_status` is the verdict: `PASS` (counted) or `FAIL` (not counted; zero
counts). `gbcms_status_reason` lists `|`-separated reasons; empty for a clean
PASS. Reasons stack (`WARN_REF_CORRECTED|MULTI_ALLELIC`). Details:
[Variant Normalization](variant-normalization.md).

| Reason | Verdict | Set when | What to do |
|:-------|:--------|:---------|:-----------|
| `WARN_REF_CORRECTED` | PASS | The REF matches the reference at ≥90% of its bases; the reference bases are used | Check the input's REF; counts are for the corrected REF |
| `WARN_HOMOPOLYMER_DECOMP` | PASS | `--rescue-homopolymer` only: the corrected homopolymer allele won | The row reports that allele's counts, not the given one |
| `MULTI_ALLELIC` | PASS | Another input row is a different ALT at the same locus (sibling-ALT exclusion) | None; each row counts its own ALT |
| `TRACT_CLUSTER` | PASS | Another length-changing row shares the repeat tract (exclusive ALT assignment) | None; read the rows together |
| `REF_MISMATCH` | FAIL | The REF matches the reference at under 90% of its bases | Check the build and coordinates; see `REF_AT_OFFSET` below |
| `FETCH_FAILED` | FAIL | The reference region could not be read (contig missing, past its end) | Check the contig naming and the FASTA |
| `EMPTY_ALLELE` | FAIL | An allele is empty | Fix the input row (a MAF dash allele is `-`) |
| `NON_SEQUENCE_ALLELE` | FAIL | An allele is not a base sequence (an IUPAC code, a stray character, `-` outside MAF input) | Fix the input row |
| `ALT_EQUALS_REF` | FAIL | ALT equals REF | Fix the input row |
| `ALT_CONTAINS_N` | FAIL | The ALT contains `N` | Fix the input row |

## Diagnostics

`gbcms_diagnostic` (VCF `GD`): what the reads, or the reference, show about the
row. Semicolon-separated in MAF; `|` in VCF. Set after counting on a PASS row; on
a FAIL row only `REF_AT_OFFSET`. Parametric flags count under their base name in
the run log.

| Flag | Mode | Set when | What to do |
|:-----|:-----|:---------|:-----------|
| `ZERO_ALT` | all | No read carries the ALT | None, or check the input allele (see `OBSERVED_ALLELE`) |
| `PARTIAL_DOMINANT` | all | `partial_alt` exceeds `alt_count`: more reads carry part of the event or another length | Inspect: often a coexisting or mis-described allele |
| `MNP_DISC_RATIO(n/m)` | all | Every MNP row: `n` positions where REF and ALT differ, of its `m` bases | Describes the MNP's shape |
| `MNP_RESCUE_ELIGIBLE` | all | `n/m` is at most `--rescue-mnp-threshold`: the MNP is a rescue candidate | With `--rescue-mnp`, rescue runs on it |
| `HIGH_N_FRACTION(f)` | all | More than 5% of the depth (`f`) has `N` at the locus | A masking hotspot (duplex); counts are thinner |
| `CLIP_CANDIDATES(n)` | DNA | An insertion with no confirmed ALT where `n` (≥ 2) reads carry a long soft clip in the insert's reach | Inspect in IGV: carriers may be clipped |
| `OBSERVED_ALLELE(chrom:pos:REF>ALT:n/0)` | all | No read carries the given allele; `n` reads carry this one | The input allele is likely mis-described |
| `COEXISTING_ALLELE(chrom:pos:REF>ALT:n/m)` | all | The given allele is carried (`m`), but more reads (`n`) carry another in the same stretch | A caveat for reading the VAF (germline indel, stutter) |
| `SW_FALLBACK(n)` | all | `n` depth reads could not be scored by the haplotype matrix and fell back to Smith-Waterman (or none) | Counts came partly from another scorer |
| `SPLICE_SKIP_DOMINANT(n)` | RNA | More reads splice over a deletion-type locus than confirm ALT | Carriers may exist as junction reads |
| `NON_DISCRIMINATING_LOCUS` | all | A nearby germline combination makes REF and ALT the same sequence | The locus cannot be genotyped |
| `RESCUED_COMPONENT(chrom:pos:REF>ALT)` | `--rescue-mnp` | The row's counts are that component SNV's | Read `gbcms_rescue` |
| `REF_AT_OFFSET(k)` | all (FAIL rows) | A `REF_MISMATCH` row's REF (3+ bases) matches the reference exactly `k` bases from its position (every such offset within 3, nearest first) | Correct the input's coordinates; the row is not moved or counted |

## MNP rescue

`gbcms_rescue` (VCF `GR`), with `--rescue-mnp` only: an audit trail
`method=decomposed;outcome=<outcome>;original_ref=…`. See
[Architecture → MNP Rescue Pass](architecture.md#mnp-rescue-pass-rescue-mnp-v430).

| Outcome | Set when | Row's counts |
|:--------|:---------|:-------------|
| `rescued` | A component SNV has more ALT than the MNP | The adopted component's (`RESCUED_COMPONENT`) |
| `skipped_grouped` | The MNP is grouped with co-annotated rows, whose reads it is assigned against | The MNP's |
| `haplotype_confirmed` | At least one read shows every changed base | The MNP's |
| `no_improvement` | No component beats the MNP | The MNP's |
| `ref_validation_failed` | A component failed REF validation | The MNP's |

## ASJD splice flags

RNA with `--gtf`. `asjd_flag` is `True` when REF and ALT junction usage differ
(Fisher p < 0.05; `asjd_qval` is BH-corrected). `asjd_diagnostic` lists QC flags,
semicolon-separated; counts are per fragment. Motifs (`asjd_ref_motif`,
`asjd_alt_motif`) are `GT-AG`, `GC-AG`, `AT-AC`, `OTHER` or `UNKNOWN`. Details:
[Output Formats → ASJD](output-formats.md#asjd-diagnostic-flags).

| Flag | Set when | What to do |
|:-----|:---------|:-----------|
| `LOW_ALT_JUNC` | Fewer than 5 ALT fragments with a junction | Too little ALT evidence for ASJD |
| `LOW_REF_JUNC` | Fewer than 10 REF fragments with a junction | Too little REF baseline |
| `NOVEL_ALT_JUNC` | The ALT's dominant junction differs from REF's and is unannotated | Candidate aberrant splicing |
| `NON_CANONICAL_MOTIF` | That junction's motif is not GT-AG, GC-AG or AT-AC | Likely a mapping artifact |
| `STRAND_DISCORDANT` | Mixed transcript-strand support (minority ≥ 0.30) on the ALT junction | Likely an alignment artifact |
| `MULTI_JUNCTION` | ALT fragments use more than 2 junctions | Complex splicing event |
| `RETENTION_DOMINANT(n)` | Near an intron boundary, `n` ≥ 10 fragments splice over the locus and outnumber the classified ones | `vaf` is the VAF among intron-retaining reads |
| `NOVEL_JUNC_AT_SPLICE_LOSS(n@start-end)` | An unannotated junction anchored at the boundary is carried by `n` ≥ 5 spliced-over fragments, more than confirm ALT | The mutant's splicing outcome is visible while `alt_count` is not |

## QC columns

| Column | Mode | Set when | What to do |
|:-------|:-----|:---------|:-----------|
| `rna_editing_site` (VCF `RED`) | RNA, `--rna-editing-db` | The locus is a known A-to-I editing site | An A>G / T>C call there is likely editing |
| `mfsd_ks_valid` | `--mfsd` | Both ALT and REF have ≥ 5 fragments in the size window | When `False`, the KS statistics are NA |
| `mfsd_ch_flag` | `--mfsd` | The variant is in a clonal-hematopoiesis gene | Read with the CH-LIKE class |
| `mfsd_alt_confidence` | `--mfsd` | `HIGH` with 5 or more ALT fragments in the size window, `LOW` with 1–4, `NONE` with none | Weigh the mFSD statistics by it |

## mFSD report classes

The `--mfsd-report` HTML classifies each variant's fragment sizes (see
[mFSD Report](mfsd-report.md)):

| Class | Meaning |
|:------|:--------|
| `TUMOR-LIKE` | ALT enriched for sub-nucleosomal fragments (enrichment > 1.3, or ALT short fragments where REF has none), KS q < 0.05, gene not in the CH set: a ctDNA signal |
| `CH-LIKE` | A CH gene, enrichment < 1.2, KS q > 0.05: ALT sizes mirror REF, a clonal-hematopoiesis signal |
| `AMBIGUOUS` | Neither: mixed signals |
| `INSUFFICIENT` | ALT fragments below `--mfsd-report-min-alt`, or no defined statistics |

In VCF, the mFSD summary is in INFO: `MFSD_REF_COUNT`, `MFSD_ALT_COUNT`,
`MFSD_REF_LLR`, `MFSD_ALT_LLR`, `MFSD_DELTA_ALT_REF`, `MFSD_KS_ALT_REF`,
`MFSD_PVAL_ALT_REF`, `MFSD_QVAL_ALT_REF`, `MFSD_SUB_NUC_REF_FRAC`,
`MFSD_SUB_NUC_ALT_FRAC`, `MFSD_SUB_NUC_ENRICHMENT`, `MFSD_MONO_NUC_REF_FRAC` and
`MFSD_MONO_NUC_ALT_FRAC` ([Output Formats → INFO](output-formats.md#info-fields)).

## Record shapes in VCF output

| Field | Set when | Meaning |
|:------|:---------|:--------|
| ALT `<NON_SEQUENCE>` | MAF input, a row whose allele is not a base sequence (or empty) | The row is FAIL (reason in `GSR`); REF is the reference base at POS |
| INFO `MAF_START`, `MAF_REF`, `MAF_ALT` | MAF input | The MAF row the record came from, as written: look a result up by its input |
| INFO `NORM_POS`, `NORM_REF`, `NORM_ALT` | `--show-normalization` | The left-aligned record gbcms counted |
