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

| `gbcms_status` | `gbcms_status_reason` | Meaning |
|:-------|:--------|:--------|
| **PASS** | *(empty)* | Variant validated against reference, ready for counting |
| **PASS** | **WARN_REF_CORRECTED** | REF allele ≥90% match; corrected to FASTA REF |
| **PASS** | **WARN_HOMOPOLYMER_DECOMP** | Variant spans a homopolymer; dual-counted with corrected allele (corrected won) |
| **PASS** | **MULTI_ALLELIC** | Variant overlaps another variant at the same locus; sibling ALT exclusion active |
| **PASS** | **TRACT_CLUSTER** | Variant shares a repeat-tract scan window with a co-annotated length-changing variant (spans need not touch); sibling exclusion and exclusive AD assignment active |
| **FAIL** | **REF_MISMATCH** | REF allele does not match the reference genome at the stated position |
| **FAIL** | **FETCH_FAILED** | Could not fetch the reference region (chromosome not found, etc.) |
| **FAIL** | **EMPTY_ALLELE** | Empty REF or ALT (malformed indel) |
| **FAIL** | **ALT_EQUALS_REF** | ALT equals REF (any case; `-` for both in a MAF): no change to count |
| **FAIL** | **ALT_CONTAINS_N** | ALT allele contains an `N` base |

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
| **Parity test** | Test verifying `count_bam` and `count_bam_binned` produce identical results |
| **Complex Del+SNV** | A deletion-format variant (`len(REF) > len(ALT) == 1`) where the anchor base also substitutes (`alt[0] ≠ ref[0]`). Examples: `GC→T`, `AG→T`. Routed to `check_complex` instead of `check_deletion`. |
| **has_shifted_same_length** | Flag routing a same-length indel candidate to Phase-3 arbitration (partial evidence propagated on non-ALT). Deletions: a windowed same-length D (≥5bp) whose deleted bases fail the S3 check (BWA left-alignment shifts). Insertions: a windowed same-length I whose inserted bases fail the sequence comparison, or a strict-anchor same-length I whose bases are unverifiable (all below `min_baseq`, or the insert runs past the read end). |
| **has_wrong_length_nearby** | Flag for a windowed WRONG-length op — a distinct-allele candidate. Set for deletions only when the op is ≥5bp (1–4bp windowed Ds are alignment noise → plain REF); insertions at any size. Resolved without Phase 3: repeat tract → neither + `partial_alt`; unique context → REF + `partial_alt`. |
| **left-alignment shift** | BWA repositions an indel's anchor to the leftmost equivalent position. For partial-repeat regions, this places the anchor several bases left of where the CIGAR `D` operation appears in reads. |
| **is_worth_realignment** | Predicate in `check_complex` returning `true` when a read's CIGAR contains indel evidence in the variant window. If `false`, M-block anchor coverage check classifies clean REF reads for deletion-direction variants. |

## Related

- [Architecture](architecture.md) — System design
- [Input Formats](input-formats.md) — File specifications
- [Complex Indels](complex-indels.md) — Real-world case studies (Del+SNV, large DEL, shifted deletion)
- [RNA CLI Reference](../cli/rna.md) — RNA-specific options
