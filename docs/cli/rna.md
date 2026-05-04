# gbcms rna

Count alleles in RNA-seq BAMs with transcriptome-aware filtering.

RNA mode extends the core DNA counting engine with filters and metrics specific to RNA-seq data: STAR-aware MAPQ rescue, dUTP strandedness filtering, splice junction tracking, and A-to-I RNA editing site flagging.

## When to Use


| Scenario | Command |
|:---------|:--------|
| cfDNA (ACCESS, IMPACT) | `gbcms dna` |
| WGS / WES / Panel | `gbcms dna` |
| STAR + dUTP RNA-seq | `gbcms rna` |
| Unstranded RNA-seq | `gbcms rna --no-strandedness` |

## Synopsis

```bash
gbcms rna [OPTIONS] --variants <FILE> --bam <NAME:PATH>... --fasta <FILE>
```

## RNA Pipeline Overview

```mermaid
flowchart LR
    subgraph Inputs ["Inputs"]
        direction TB
        BAM["BAM Files"]
        VCF["VCF/MAF"]
        FASTA["Reference"]
        DB["REDIportal DB"]:::optional
    end

    subgraph Filters ["RNA Filters"]
        direction TB
        NH["NH:i:1 MAPQ Rescue"]:::filter
        Strand["Strandedness Filter"]:::filter
    end

    subgraph Counting ["Allele Classification"]
        direction TB
        Classify["Phase 1-3 Classification"]
        Sense["Sense/Antisense Tracking"]
        Splice["Splice Junction Detection"]
        Edit["RNA Editing Flagging"]
    end

    BAM --> NH --> Strand --> Classify
    VCF --> Classify
    FASTA --> Classify
    Classify --> Sense
    Classify --> Splice
    Classify --> Edit
    DB -.->|"optional"| Edit
    Sense --> Out(["RNA VCF/MAF"]):::pass
    Splice --> Out
    Edit --> Out

    classDef filter fill:#e67e22,color:#fff,stroke:#bf6516,stroke-width:2px;
    classDef optional fill:#95a5a6,color:#fff,stroke:#7f8c8d,stroke-width:1px,stroke-dasharray:4;
    classDef pass fill:#27ae60,color:#fff,stroke:#1e8449,stroke-width:2px;
    class NH,Strand filter;
```

---

## Required Arguments

!!! info "Shared Arguments"
    RNA mode shares all [required arguments](dna.md#required-arguments), [output options](dna.md#output-options), [filtering options](dna.md#filtering-options), [BAQ options](dna.md#baq-options), [UMI options](dna.md#umi-options), and [debugging options](dna.md#debugging-options) with DNA mode. See the [`gbcms dna` reference](dna.md) for full descriptions.

!!! tip "BAQ is On by Default in RNA Mode"
    Unlike DNA mode, `--apply-baq` defaults to **on** for RNA. BAQ penalizes bases within 5 bp of both indels and splice junctions (CIGAR `N`), serving as a lightweight alternative to GATK SplitNCigarReads. Disable with `--no-baq` if your upstream workflow already runs SplitNCigarReads or BQSR. See [RNA Splice-Junction Handling](../reference/rna-splice-handling.md).

| Option | Description |
|:-------|:------------|
| `--variants`, `-v` | VCF or MAF file |
| `--bam`, `-b` | BAM file(s) |
| `--fasta`, `-f` | Reference FASTA (with .fai index) |
| `--output-dir`, `-o` | Output directory |

---

## RNA-Specific Options

These options are **only available** on `gbcms rna`, not on `gbcms dna`.

### Strandedness Filtering

| Option | Default | Description |
|:-------|:--------|:------------|
| `--enforce-strandedness/--no-strandedness` | `true` | Filter reads by dUTP strand orientation relative to gene strand |

!!! info "Biological Context: dUTP Stranded Libraries"
    In dUTP-stranded RNA-seq, the second strand (synthesized with dUTP) is degraded, so sequenced reads reflect the **antisense** strand of the original mRNA. The strandedness filter uses the variant's `gene_strand` annotation (from the input MAF) to determine whether each read's orientation is consistent with the expected transcript direction.

    **Disable** with `--no-strandedness` for unstranded RNA-seq libraries where read orientation is random.

```mermaid
flowchart TD
    Read(["📖 RNA Read"]) --> HasGS{"gene_strand annotated?"}
    HasGS -->|No| Pass(["✅ Pass — no filter"]):::pass
    HasGS -->|Yes| Check{"Read orientation\nvs gene strand?"}
    Check -->|"Consistent"| Sense(["✅ Sense — counted in DP/AD/RD"]):::pass
    Check -->|"Inconsistent"| Anti(["📊 Antisense — counted in rna_antisense_depth"]):::anti

    classDef pass fill:#27ae60,color:#fff,stroke:#1e8449,stroke-width:2px;
    classDef anti fill:#e67e22,color:#fff,stroke:#bf6516,stroke-width:2px;
```

### RNA Editing Database

| Option | Default | Description |
|:-------|:--------|:------------|
| `--rna-editing-db` | _(none)_ | Path to REDIportal TABLE1 file of known A-to-I RNA editing sites (`.gz` supported) |

!!! info "Biological Context: A-to-I RNA Editing"
    ADAR (Adenosine Deaminase Acting on RNA) enzymes convert adenosine to inosine in double-stranded RNA regions. Inosine is read as guanosine by the sequencing machinery, creating apparent **A→G variants** that are biological RNA modifications — not somatic mutations. Flagging these sites prevents false-positive variant calls.

    The REDIportal database catalogs >16 million known human A-to-I editing sites from multiple tissues and conditions.

!!! tip "Obtaining REDIportal"
    Download from [REDIportal](http://srv00.recas.ba.infn.it/atlas/):
    ```bash
    # Download TABLE1 format (tab-separated, 1-based coordinates)
    wget http://srv00.recas.ba.infn.it/atlas/download/TABLE1_hg38.txt.gz
    ```
    The file is loaded once at startup and stored as a hash set for O(1) lookup per variant.

---

## Quality Thresholds

RNA mode uses **different MAPQ default** to accommodate STAR aligner behavior:

| Option | DNA Default | RNA Default | Rationale |
|:-------|:------------|:------------|:----------|
| `--min-mapq` | 20 | **1** | STAR assigns MAPQ 255 to unique, 0–3 to multi-mappers. Low MAPQ threshold with NH:i:1 rescue captures novel splice junctions. |
| `--min-baseq` | 20 | 20 | Same |
| `--fragment-qual-threshold` | 10 | 10 | Same |

---

## RNA Read Filters

RNA mode extends the standard [read filter cascade](../reference/read-filters.md) with two additional checks:

### NH:i:1 MAPQ Rescue

STAR assigns MAPQ=0 to reads that map to multiple loci. However, reads with `NH:i:1` (Number of Hits = 1) are uniquely mapped despite low MAPQ — they were assigned low scores because STAR hadn't observed the splice junction before.

```mermaid
flowchart TD
    Read(["📖 RNA Read"])
    Read -->|"read.mapq"| MAPQ{"MAPQ ≥ min_mapq?\n(default: 1)"}
    MAPQ -->|"Yes — MAPQ=255\nunique alignment"| Pass(["✅ Pass"]):::pass
    MAPQ -->|"No — MAPQ=0\n(check NH tag)"| NH{"NH:i:1?"}
    NH -->|"Yes — novel junction\nunique but low MAPQ"| Rescue(["✅ Rescued"]):::rescue
    NH -->|"No — NH>1\nmulti-mapper"| Drop(["❌ Discard"]):::drop

    classDef pass fill:#27ae60,color:#fff,stroke:#1e8449,stroke-width:2px;
    classDef rescue fill:#f39c12,color:#fff,stroke:#d68910,stroke-width:2px;
    classDef drop fill:#e74c3c,color:#fff,stroke:#c0392b,stroke-width:2px;
```

!!! info "Why MAPQ=1 Default?"
    Setting `--min-mapq 1` together with NH rescue ensures:

    - **Uniquely mapped** reads (MAPQ=255): ✅ always pass
    - **Novel junction** reads (MAPQ=0, NH=1): ✅ rescued
    - **True multi-mappers** (MAPQ=0, NH>1): ❌ filtered

### Filter Defaults

RNA mode uses the **same filter defaults** as DNA mode — duplicates, secondary, supplementary, and QC-failed reads are all filtered by default in both modes:

| Filter | DNA Default | RNA Default |
|:-------|:------------|:------------|
| `--filter-secondary` | on | on (same) |
| `--filter-supplementary` | on | on (same) |
| `--filter-qc-failed` | on | on (same) |

!!! caution "Known Limitation: Variants at Exon-Intron Boundaries"
    For variants positioned very close to a splice junction, the D6 consensus intron-snipping step (which removes introns from the `ref_context` before Phase 3 alignment) requires **>50% of local reads to carry a CIGAR `N` op** at the same position. If the variant sits at a boundary where only ~40% of reads span the junction with an `N` op — e.g., because many soft-clipped reads end just before the junction — the intron is **not snipped**.

    In this case Phase 3 aligns the read against an unspliced (intron-containing) reference context. The alignment score is lower than it would be against the mature mRNA sequence, and the read may fall back to `neither` rather than being classified as REF or ALT. This affects `dp` (depth) but not correctness for reads that *do* span the junction cleanly.

    **Workaround:** Use `--trace` logging and grep for `D6 splice` to confirm whether intron snipping fired for the variant of interest:
    ```bash
    RUST_LOG=trace gbcms rna ... 2>&1 | grep "D6 splice\|<chrom>:<pos>"
    ```

---

## Alignment Backend

RNA mode uses the same **two-stage Phase 3 pipeline** as DNA:

1. **WFA fast-path** (`wfa2lib-rs`) — edit-distance triage; resolves ~70-80% of reads instantly.
2. **Marginalized PairHMM** (escalated when WFA is ambiguous) — uses per-base quality probabilities for confident LLR scoring.

RNA uses **relaxed gap penalties** to tolerate reverse transcriptase (RT) stutter artifacts:

| Option | DNA Default | RNA Default | Rationale |
|:-------|:------------|:------------|:----------|
| `--alignment-backend` | `pairhmm` | `pairhmm` | Same (WFA + PairHMM) |
| `--gap-open-prob` | `1e-4` | **`5e-3`** | RT introduces more gap artifacts than DNA polymerase |
| `--gap-extend-prob` | `0.1` | **`0.25`** | RT stutter extends gaps more frequently |

!!! info "Biological Context: RT Stutter"
    Reverse transcriptase (used in cDNA synthesis) has lower fidelity than DNA polymerases and frequently introduces insertions/deletions at homopolymer runs and microsatellite regions. The relaxed gap penalties prevent these artifacts from being classified as ALT-supporting reads.

!!! info "Complex Indel Support in RNA Mode"
    All complex indel classification improvements — Del+SNV routing (`GC→T` pattern), large deletion REF visibility, and left-alignment shifted deletion rescue — apply identically to RNA mode. The allele classification engine (`variant_checks.rs`) has no mode branches; RNA reads go through the same Phase 1/2/3 pipeline as DNA reads.

    See [Complex Indels](../reference/complex-indels.md) for case studies and read-level examples.

    **Note on left-alignment shift (Fix 4):** This fix was motivated by BWA's indel left-alignment behaviour, where the variant anchor can sit several bases left of where the CIGAR `D` op appears in reads. STAR does not left-align indels in the same way, so this scenario is less common in RNA BAMs — but it still applies when STAR places a CIGAR `D` slightly offset from the annotated position, or when the input MAF was derived from a BWA-aligned DNA call set being counted against an RNA BAM.

---

## RNA-Specific Output Columns

??? info "RNA-Specific Output Columns (click to expand)"

    When running in RNA mode, additional columns capture transcriptome-specific metrics. These columns are **absent** from DNA mode output.

    ### MAF Columns (5 additional)

    | Column | Type | Description |
    |:-------|:-----|:------------|
    | `rna_sense_depth` | u32 | Reads aligning to the gene **sense** strand |
    | `rna_antisense_depth` | u32 | Reads aligning to the gene **antisense** strand |
    | `rna_alt_sense_count` | u32 | ALT-classified reads on the sense strand |
    | `rna_editing_site` | bool | Variant overlaps a known A→I editing site from `--rna-editing-db` |
    | `rna_splice_spanning` | u32 | ALT-classified reads containing splice junctions (CIGAR `N` operations) spanning the variant |

    ### VCF INFO Fields (5 additional)

    | Field | Type | Description |
    |:------|:-----|:------------|
    | `SEN` | Integer | Sense strand depth |
    | `ANT` | Integer | Antisense strand depth |
    | `ASEN` | Integer | ALT sense strand count |
    | `RED` | Flag | Known A→I RNA editing site overlap (flag; present if true) |
    | `SPL` | Integer | Splice-spanning ALT read count |

    ### VCF FORMAT Fields (4 additional)

    Per-sample fields added to the FORMAT column:

    | Field | Type | Description |
    |:------|:-----|:------------|
    | `SEN` | Integer | Per-sample sense depth |
    | `ANT` | Integer | Per-sample antisense depth |
    | `ASEN` | Integer | Per-sample ALT sense count |
    | `SPL` | Integer | Per-sample splice-spanning count |

    !!! note "RED is INFO-only"
        The `RED` editing-site flag appears only in the VCF INFO column (as a flag field), not in the per-sample FORMAT column.

---

## Examples

=== "Basic RNA"

    ```bash
    gbcms rna \
        --variants mutations.vcf \
        --bam rna_sample:star_aligned.bam \ # (1)!
        --fasta reference.fa \
        --output-dir results/
    ```

    1. STAR-aligned BAMs should have the `NH` tag for multi-mapping rescue.

=== "With Editing DB"

    ```bash
    gbcms rna \
        --variants mutations.maf \
        --bam tumor_rna:aligned.bam \
        --fasta hg38.fa \
        --rna-editing-db TABLE1_hg38.txt.gz \ # (1)!
        --format maf \
        --output-dir results/
    ```

    1.  REDIportal TABLE1 format. Download from [REDIportal](http://srv00.recas.ba.infn.it/atlas/).

=== "Unstranded Library"

    ```bash
    gbcms rna \
        --variants mutations.vcf \
        --bam unstranded:aligned.bam \
        --fasta reference.fa \
        --no-strandedness \ # (1)!
        --output-dir results/
    ```

    1.  Disables dUTP strand filtering for unstranded RNA-seq protocols.

=== "With UMI Tags"

    ```bash
    gbcms rna \
        --variants mutations.vcf \
        --bam rna_sample:umi_tagged.bam \
        --fasta reference.fa \
        --umi-tag RX \
        --output-dir results/
    ```

---

## Differences from DNA Mode

| Feature | `gbcms dna` | `gbcms rna` |
|:--------|:------------|:------------|
| **MAPQ default** | 20 | 1 (with NH rescue) |
| **Gap-open probability** | 1e-4 | 5e-3 (RT tolerance) |
| **Gap-extend probability** | 0.1 | 0.25 (RT tolerance) |
| **Secondary filter** | on | on (same) |
| **Supplementary filter** | on | on (same) |
| **QC-failed filter** | on | on (same) |
| **BAQ default** | off | **on** (no upstream BQSR; splice penalty) |
| **Strandedness filter** | N/A | enabled (dUTP) |
| **RNA editing flagging** | N/A | optional (`--rna-editing-db`) |
| **Output columns** | Standard (DP, RD, AD, etc.) | Standard **+ 5 RNA columns** |
| **mFSD analysis** | ✅ available | ❌ not available |

!!! note "No mFSD in RNA Mode"
    Mutant Fragment Size Distribution analysis is specific to cfDNA and is not available in RNA mode. The `--mfsd` and `--mfsd-parquet` options are only present on `gbcms dna`.

---

## Related

- [gbcms dna](dna.md) — DNA/cfDNA counting
- [Quick Start](../getting-started/quickstart.md) — Common patterns
- [Output Formats](../reference/output-formats.md) — Complete VCF and MAF column-level schema (RNA columns, INFO/FORMAT fields)
- [Read Filters](../reference/read-filters.md) — Filter cascade details
- [Counting & Metrics](../reference/counting-metrics.md) — How reads become counts
- [Allele Classification](../reference/allele-classification.md) — How each variant type is counted
