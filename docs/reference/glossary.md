# Glossary

Technical terms used throughout the documentation.

## Metrics

| Term | Definition |
|:-----|:-----------|
| **DP** | Total depth — all mapped, quality-filtered reads overlapping the variant anchor position, including REF, ALT, and 'neither'. `DP ≥ RD + AD`. |
| **VAF** | Variant Allele Frequency — `AD / (RD + AD)` |
| **Strand Bias** | Fisher's exact test for read direction imbalance |
| **AD** | Alternate allele depth (supporting reads) |
| **RD** | Reference allele depth |
| **Informative read** | For an indel, a read that covers a flank of the event's shift-equivalence region and runs one base past the first base where REF and ALT differ. Only informative reads count REF; an ALT read that is not informative counts only when its own bases discriminate (a truncated long insertion's carrier does). The others count toward DP only ([details](allele-classification.md#informative-reads-for-indels)). |

## File Formats

| Term | Definition |
|:-----|:-----------|
| **VCF** | Variant Call Format — standard variant file |
| **MAF** | Mutation Annotation Format — annotation-rich variant file |
| **BAM** | Binary Alignment Map — compressed alignment file |
| **BAI** | BAM Index — enables random access to BAM |
| **FASTA** | Reference genome sequence file |
| **FAI** | FASTA Index — enables random access to FASTA |

## Quality Scores

| Term | Definition |
|:-----|:-----------|
| **MAPQ** | Mapping Quality — confidence in read alignment |
| **BASEQ** | Base Quality — confidence in base call |

## cfDNA Terms

| Term | Definition |
|:-----|:-----------|
| **cfDNA** | Cell-free DNA — circulating DNA in plasma |
| **ctDNA** | Circulating tumor DNA — tumor-derived cfDNA |
| **Duplex** | Reads from both strands of original molecule |

## Alignment Backends

| Term | Definition |
|:-----|:-----------|
| **PairHMM** | Default Phase 3 alignment backend for all modes. Uses pair hidden Markov model with base quality probabilities for probabilistic alignment scoring. Selected via `--alignment-backend pairhmm`. |
| **Smith-Waterman (SW)** | Alternative Phase 3 alignment backend. Uses edit-distance scoring to align reads against REF and ALT haplotypes. Confident call requires ≥2 score margin. Selected via `--alignment-backend sw`. |
| **LLR** | Log-Likelihood Ratio — PairHMM confidence metric. `LLR = ln(P(read|ALT)) - ln(P(read|REF))`. Default threshold: 2.3 (≈ ln(10), i.e., 10:1 odds). |

## Fragment Metrics

| Term | Definition |
|:-----|:-----------|
| **DPF** | Fragment depth — total unique fragments (including discarded) |
| **RDF** | Reference fragment count — fragments resolved to REF |
| **ADF** | Alternate fragment count — fragments resolved to ALT |
| **Fragment consensus** | Quality-weighted method to resolve R1/R2 disagreements within a fragment |

## Validation Status

Status is two fields: `gbcms_status` (verdict `PASS`/`FAIL`) and `gbcms_status_reason`
(`|`-separated reason tags, empty for a clean PASS).

The reasons and what each means: [QC Flags → Verdict and status reasons](qc-flags.md#verdict-and-status-reasons).

## RNA Terms

| Term | Definition |
|:-----|:-----------|
| **Sense strand** | mRNA coding strand; matches gene transcription direction |
| **Antisense strand** | Non-coding strand; template for transcription |
| **A-to-I editing** | Post-transcriptional RNA modification by ADAR enzymes; adenosine → inosine (read as guanosine) |
| **NH tag** | BAM tag indicating number of reported alignments for a read |
| **Splice junction** | CIGAR `N` operation representing an intron skip in RNA-seq reads |
| **REDIportal** | Database of known human A-to-I RNA editing sites |

## Engine Terms

| Term | Definition |
|:-----|:-----------|
| **BAQ** | Base Alignment Quality — heuristic quality downgrade (−20 BQ within 5 bp) near indels and splice junctions. Off by default for DNA (upstream BQSR), on by default for RNA (no upstream BQ calibration). See [Read Filters: BAQ](read-filters.md#baq-quality-downgrade). |
| **UMI** | Unique Molecular Identifier — barcode for molecule-level deduplication |
| **Genomic binning** | Variant grouping by chromosome for cache-efficient counting (`count_bam_binned`) |
| **Binning invariance** | The rule that counts never depend on how variants are grouped into bins (window, cap, one variant per bin or per call); checked by `tests/test_binning_invariance.py` and every `count_checked` call |
| **Read census** | The tests' classification oracle (`tests/census.py`): each read judged by its own bases across the variant's tract |
| **Complex Del+SNV** | A deletion-format variant (`len(REF) > len(ALT) == 1`) where the anchor base also substitutes (`alt[0] ≠ ref[0]`). Examples: `GC→T`, `AG→T`. Judged by the [exact-carrier rule](allele-classification.md#the-exact-carrier-rule) instead of `check_deletion`. |
| **Complex Ins+SNV** | An insertion-format variant (`len(REF) == 1 < len(ALT)`) where the anchor base also substitutes (`alt[0] ≠ ref[0]`). Examples: `C→TA`, `A→CCC`. Judged by the exact-carrier rule instead of `check_insertion` (since 6.6.0, #121). |
| **has_in_band_mismatch** | Flag for an in-band large deletion (≥50bp, another length within 3bp) whose deleted bases differ or cannot be verified: Phase-3 arbitration, with partial evidence on non-ALT. The only windowed candidate left for Phase 3 since 6.6.0. |
| **has_shifted_same_length** (removed in 6.6.0) | Routed a windowed same-length insertion of other bases to Phase-3 arbitration, whose closer haplotype let length win ALT and absorbed other alleles into REF. Such a read is now judged by its bases (RJ-21); the deletion side's in-band case is `has_in_band_mismatch`. |
| **has_unreadable_insert** | Flag for an insertion of the variant's length none of whose bases can be read (each N or below `min_baseq`, or the insert runs past the read end), at the anchor or inside the discrimination window: the read carries the length, not the sequence. Neither + `partial_alt` unless the variant is also written readably in the read; never Phase 3 (RJ-20, since 6.6.0). Outside the window it is a separate event, as a readable insert there is. |
| **has_distinct_allele_nearby** | Flag for a distinct-allele candidate in the scan window: a WRONG-length op (deletions only when ≥5bp — 1–4bp windowed Ds are alignment noise → plain REF; insertions at any size), the variant's inserted bases where they give another haplotype inside its discrimination window, or the variant written elsewhere with another change (gap, insertion, splice) across that window. Resolved without Phase 3: an indel of the read inside the discrimination window → neither + `partial_alt` (RJ-7); outside it a separate event → REF, no partial evidence (RJ-8). |
| **another_allele_in_window** | Flag for a windowed insertion of the variant's length, with a readable base of its own, that is not the ALT by its bases at any placement it can take and can sit inside the discrimination window, wherever the aligner wrote it (it slides a junction when the base it places onto the reference fits it): another allele there → neither + `partial_alt` (RJ-21, since 6.6.0). |
| **left-alignment shift** | BWA repositions an indel's anchor to the leftmost equivalent position. For partial-repeat regions, this places the anchor several bases left of where the CIGAR `D` operation appears in reads. |
| **Exact carrier** | A read whose own bases (soft clips included, bases below `--min-baseq` masked) carry the whole given allele of a complex variant with two reference flank bases each side, read at the event's own position. Only exact carriers count as REF or ALT for delins, Del+SNV and structural MNP reads. See [the exact-carrier rule](allele-classification.md#the-exact-carrier-rule). |
| **is_worth_realignment** | Predicate in the previous complex classifier (`check_complex`, used only for a variant with no reference holding the event) returning `true` when a read's CIGAR contains indel evidence in the variant window. If `false`, M-block anchor coverage check classifies clean REF reads for deletion-direction variants. |

## Related

- [Architecture](architecture.md) — System design
- [Input Formats](input-formats.md) — File specifications
- [Complex Indels](complex-indels.md) — Real-world case studies (Del+SNV, large DEL, shifted deletion)
- [RNA CLI Reference](../cli/rna.md) — RNA-specific options
