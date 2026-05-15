# Output Formats

gbcms writes one output file per BAM sample. The output format, column
composition, and sample-naming strategy all depend on the CLI flags used.

!!! info "Quick Reference"
    Output file path: `{--output-dir}/{sample_name}{--suffix}.{vcf|maf}`

    `sample_name` is set by the `name:` prefix on `--bam` (e.g. `--bam tumor:tumor.bam`)
    or falls back to the BAM filename stem.

---

## How the Output Path Is Decided

The diagram below shows every decision point from CLI flags to the final
output column set. Follow your input type and desired output format to see
exactly what you get.

```mermaid
flowchart TD
    Input(["Input variants"]):::start

    Input --> InputType{"Input type?"}
    InputType -->|"VCF / VCF.GZ"| VCFIn["VCF-origin<br/>no metadata"]:::vcf
    InputType -->|MAF| MAFIn["MAF-origin<br/>full row metadata"]:::maf

    FmtChoice{"--format?"}
    VCFIn --> FmtChoice
    MAFIn --> FmtChoice

    FmtChoice -->|vcf| VCFWriter["VcfWriter"]:::writer
    FmtChoice -->|maf| MAFWriter["MafWriter"]:::writer

    VCFWriter --> ModeVCF{"Mode?"}
    ModeVCF -->|dna| DNAVCF["VCF: standard INFO + FORMAT"]:::dna
    ModeVCF -->|rna| RNAVCF["VCF: + SEN ANT ASEN RED SPL"]:::rna

    MAFWriter --> ModeMAF{"Mode?"}
    ModeMAF -->|dna| DNAMAFPath{"Input?"}
    ModeMAF -->|rna| RNAMAFPath{"Input?"}

    DNAMAFPath -->|"VCF-origin"| DNAVMAF["GDC MAF columns + gbcms counts"]:::dna
    DNAMAFPath -->|"MAF-origin"| DNAMMAF["All original columns + gbcms counts"]:::dna

    RNAMAFPath -->|"VCF-origin"| RNAVMAF["GDC MAF columns + gbcms counts<br/>+ 5 rna_* columns"]:::rna
    RNAMAFPath -->|"MAF-origin"| RNAMMAF["All original columns + gbcms counts<br/>+ 5 rna_* columns"]:::rna

    classDef start fill:#9b59b6,color:#fff,stroke:#7d3c98,stroke-width:2px
    classDef vcf fill:#2471a3,color:#fff,stroke:#1a5276,stroke-width:2px
    classDef maf fill:#117a65,color:#fff,stroke:#0e6655,stroke-width:2px
    classDef writer fill:#7d6608,color:#fff,stroke:#6d5f07,stroke-width:2px
    classDef dna fill:#1a5276,color:#fff,stroke:#154360,stroke-width:2px
    classDef rna fill:#1e8449,color:#fff,stroke:#196f3d,stroke-width:2px
```

---

## VCF Output (`--format vcf`)

A standards-compliant VCFv4.2 file with one row per variant per sample.

### File Header

The `##fileformat`, `##source`, and `##INFO`/`##FORMAT` meta-lines are
written once. RNA-specific meta-lines are only included when running
`gbcms rna` — the header is self-describing.

=== "DNA mode"

    ```
    ##fileformat=VCFv4.2
    ##source=gbcms
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
    ##INFO=<ID=GS,Number=1,Type=String,Description="gbcms normalization/counting status">
    ##INFO=<ID=GD,Number=1,Type=String,Description="gbcms post-counting diagnostic flags">
    ##INFO=<ID=GR,Number=1,Type=String,Description="gbcms rescue audit trail">
    ##INFO=<ID=AAD,Number=1,Type=Integer,Description="Any ALT Depth (any_alt = ad + partial_alt)">
    ##INFO=<ID=PAD,Number=1,Type=Integer,Description="Partial ALT Depth">
    ##INFO=<ID=NAD,Number=1,Type=Integer,Description="N-base Depth (duplex masking QC)">
    ##INFO=<ID=SB_PVAL,Number=1,Type=Float,Description="Fisher strand bias p-value">
    ##INFO=<ID=SB_OR,Number=1,Type=Float,Description="Fisher strand bias odds ratio">
    ##INFO=<ID=FSB_PVAL,Number=1,Type=Float,Description="Fisher fragment strand bias p-value">
    ##INFO=<ID=FSB_OR,Number=1,Type=Float,Description="Fisher fragment strand bias odds ratio">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (ref,alt)">
    ##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Allelic depths on forward strand (ref_fwd,alt_fwd)">
    ##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Allelic depths on reverse strand (ref_rev,alt_rev)">
    ##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele fraction (read level)">
    ##FORMAT=<ID=FAD,Number=R,Type=Integer,Description="Fragment allelic depths (ref_frag,alt_frag)">
    ##FORMAT=<ID=FADF,Number=R,Type=Integer,Description="Fragment depths on forward strand">
    ##FORMAT=<ID=FADR,Number=R,Type=Integer,Description="Fragment depths on reverse strand">
    ##FORMAT=<ID=FAF,Number=1,Type=Float,Description="Variant allele fraction (fragment level)">
    ##FORMAT=<ID=AAD,Number=1,Type=Integer,Description="Any ALT depth (alt + partial_alt)">
    ##FORMAT=<ID=PAD,Number=1,Type=Integer,Description="Partial ALT depth">
    ##FORMAT=<ID=NAD,Number=1,Type=Integer,Description="N-base depth">
    #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  <sample_name>
    ```

=== "RNA mode"

    ```
    ##fileformat=VCFv4.2
    ##source=gbcms
    ##INFO=<ID=DP,...>
    ##INFO=<ID=GS,...>
    ##INFO=<ID=GD,...>
    ##INFO=<ID=GR,...>
    ##INFO=<ID=AAD,...>
    ##INFO=<ID=PAD,...>
    ##INFO=<ID=NAD,...>
    ##INFO=<ID=SB_PVAL,...>
    ##INFO=<ID=SB_OR,...>
    ##INFO=<ID=FSB_PVAL,...>
    ##INFO=<ID=FSB_OR,...>
    ##INFO=<ID=SEN,Number=1,Type=Integer,Description="Reads on the transcript sense strand">
    ##INFO=<ID=ANT,Number=1,Type=Integer,Description="Reads on the antisense strand">
    ##INFO=<ID=ASEN,Number=1,Type=Integer,Description="ALT reads on the transcript sense strand">
    ##INFO=<ID=RED,Number=0,Type=Flag,Description="Locus is a candidate A-to-I RNA editing site">
    ##INFO=<ID=SPL,Number=1,Type=Integer,Description="ALT reads spanning a splice junction (CIGAR N)">
    ##FORMAT=<ID=GT,...>
    ##FORMAT=<ID=DP,...>
    ##FORMAT=<ID=AD,...>
    ##FORMAT=<ID=ADF,...>
    ##FORMAT=<ID=ADR,...>
    ##FORMAT=<ID=VAF,...>
    ##FORMAT=<ID=FAD,...>
    ##FORMAT=<ID=FADF,...>
    ##FORMAT=<ID=FADR,...>
    ##FORMAT=<ID=FAF,...>
    ##FORMAT=<ID=AAD,...>
    ##FORMAT=<ID=PAD,...>
    ##FORMAT=<ID=NAD,...>
    ##FORMAT=<ID=SEN,Number=1,Type=Integer,Description="Sense strand depth">
    ##FORMAT=<ID=ANT,Number=1,Type=Integer,Description="Antisense strand depth">
    ##FORMAT=<ID=ASEN,Number=1,Type=Integer,Description="ALT sense strand count">
    ##FORMAT=<ID=SPL,Number=1,Type=Integer,Description="Splice-spanning ALT count">
    #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  <sample_name>
    ```

---

### Fixed Fields

| Column | Source | Notes |
|:-------|:-------|:------|
| `CHROM` | Variant chromosome | Preserved from input |
| `POS` | Variant position | 1-based (VCF convention) |
| `ID` | Original VCF `ID` field | `.` when input is MAF (no `ID` column) |
| `REF` | Reference allele | From input; validated against FASTA |
| `ALT` | Alternate allele | From input |
| `QUAL` | `.` | Always missing — gbcms does not perform variant calling |
| `FILTER` | `.` | Not set |

---

### INFO Fields

The `INFO` column is a semicolon-separated list of `KEY=VALUE` pairs.

=== "Always present (DNA + RNA)"

    | Field | Type | Description |
    |:------|:-----|:------------|
    | `DP` | Integer | Total read depth at position |
    | `GS` | String | gbcms normalization/counting status. Pipe-separated multi-value in VCF (e.g., `PASS\|WARN_REF_CORRECTED`). Semicolons in MAF. |
    | `GD` | String | Post-counting diagnostic flags. Pipe-separated in VCF (e.g., `ZERO_ALT\|PARTIAL_DOMINANT`). Semicolons in MAF. `.` if none. |
    | `GR` | String | Rescue audit trail. Pipe-separated key=value pairs. `.` if no rescue attempted. |
    | `AAD` | Integer | Any ALT Depth — reads with ALT evidence at ≥1 discriminating position. Invariant: `AAD = AD + PAD` |
    | `PAD` | Integer | Partial ALT Depth — reads matching ALT at some but not all discriminating positions. Populated for all variant types including INDELs (via Phase 3 structural evidence propagation). |
    | `NAD` | Integer | N-base Depth — reads with N base at ≥1 discriminating position (duplex masking QC metric) |
    | `SB_PVAL` | Float | Fisher's exact test p-value for read-level strand bias |
    | `SB_OR` | Float | Fisher's exact test odds ratio for read-level strand bias |
    | `FSB_PVAL` | Float | Fragment-level strand bias p-value |
    | `FSB_OR` | Float | Fragment-level strand bias odds ratio |

=== "RNA mode only"

    | Field | Type | Description |
    |:------|:-----|:------------|
    | `SEN` | Integer | Total reads on the transcript **sense** strand |
    | `ANT` | Integer | Total reads on the **antisense** strand |
    | `ASEN` | Integer | ALT reads on the sense strand |
    | `SPL` | Integer | ALT reads spanning a splice junction (reads with `N` CIGAR op) |
    | `RED` | Flag | Present when the locus overlaps a known A-to-I RNA editing site (requires `--rna-editing-db`) |

=== "--mfsd only"

    | Field | Type | Description |
    |:------|:-----|:------------|
    | `MFSD_DELTA_ALT_REF` | Float | mean(ALT) − mean(REF) fragment size delta (bp) |
    | `MFSD_KS_ALT_REF` | Float | 2-sample KS D-statistic (ALT vs REF fragments) |
    | `MFSD_PVAL_ALT_REF` | Float | KS test p-value (ALT vs REF) |
    | `MFSD_ALT_LLR` | Float | Log-likelihood ratio for ALT fragments vs healthy/tumor Gaussian model |
    | `MFSD_REF_LLR` | Float | Log-likelihood ratio for REF fragments |
    | `MFSD_ALT_COUNT` | Integer | ALT-classified fragments in 50–1000 bp size window |
    | `MFSD_REF_COUNT` | Integer | REF-classified fragments in 50–1000 bp size window |

=== "--gtf only (RNA mode)"

    These fields are emitted **only** when `gbcms rna --gtf <file>` is provided.

    | Field | Type | Description |
    |:------|:-----|:------------|
    | `EBD` | Integer | Distance to nearest annotated exon boundary (`.` when no GTF) |
    | `TXRC` | String | Per-transcript read counts. Format: `ENST:AD,RD,DP\|ENST:AD,RD,DP` |
    | `TXFC` | String | Per-transcript fragment counts. Format: `ENST:ADF,RDF,DPF\|ENST:ADF,RDF,DPF` |
    | `ASJD` | Flag | Allele-Specific Junction Divergence detected |
    | `ASJDP` | Float | ASJD raw Fisher exact p-value |
    | `ASJDQ` | Float | ASJD BH-corrected q-value |
    | `ASJDRJ` | String | REF dominant junction (`start-end`) |
    | `ASJDAJ` | String | ALT dominant junction (`start-end`) |
    | `ASJDRM` | String | REF splice motif (GT-AG/GC-AG/AT-AC/OTHER/UNKNOWN) |
    | `ASJDAM` | String | ALT splice motif |
    | `ASJDRK` | Integer | REF junction in GTF (1/0) |
    | `ASJDAK` | Integer | ALT junction in GTF (1/0) |
    | `ASJDNR` | Integer | REF reads on dominant junction |
    | `ASJDNA` | Integer | ALT reads on dominant junction |
    | `ASJDD` | String | ASJD diagnostic flags (pipe-separated) |

=== "--show-normalization only"

    | Field | Type | Description |
    |:------|:-----|:------------|
    | `NORM_POS` | Integer | Left-aligned VCF position (1-based) after normalization |
    | `NORM_REF` | String | Left-aligned REF allele |
    | `NORM_ALT` | String | Left-aligned ALT allele |

---

### FORMAT Fields

=== "DNA mode"

    `FORMAT` column: `GT:DP:AD:ADF:ADR:VAF:FAD:FADF:FADR:FAF:AAD:PAD:NAD`

    | Tag | Values | Description |
    |:----|:-------|:------------|
    | `GT` | `0/0` or `0/1` | Diploid genotype — `0/1` when any ALT reads present |
    | `DP` | integer | Total read depth (single integer, VCF spec) |
    | `AD` | `ref,alt` | Allelic depths — ref_total,alt_total (Number=R) |
    | `ADF` | `ref_fwd,alt_fwd` | Forward strand per allele (bcftools convention) |
    | `ADR` | `ref_rev,alt_rev` | Reverse strand per allele |
    | `VAF` | float | Variant allele fraction at read level |
    | `FAD` | `ref_frag,alt_frag` | Fragment allelic depths (Number=R) |
    | `FADF` | `ref_frag_fwd,alt_frag_fwd` | Fragment forward strand per allele |
    | `FADR` | `ref_frag_rev,alt_frag_rev` | Fragment reverse strand per allele |
    | `FAF` | float | Variant allele fraction at fragment level |
    | `AAD` | integer | Any ALT Depth (reads with any ALT evidence) |
    | `PAD` | integer | Partial ALT Depth (partial ALT match only) |
    | `NAD` | integer | N-base Depth (reads with N at discriminating position) |

=== "RNA mode"

    `FORMAT` column: `GT:DP:AD:ADF:ADR:VAF:FAD:FADF:FADR:FAF:AAD:PAD:NAD:SEN:ANT:ASEN:SPL`

    All DNA fields above (including `AAD`, `PAD`, `NAD`), plus:

    | Tag | Values | Description |
    |:----|:-------|:------------|
    | `SEN` | integer | Sense-strand read depth |
    | `ANT` | integer | Antisense-strand read depth |
    | `ASEN` | integer | ALT count on sense strand |
    | `SPL` | integer | Splice-junction-spanning ALT count |

---

### Annotated Example

```vcf
#CHROM  POS     ID      REF  ALT  QUAL  FILTER  INFO                                              FORMAT           sample1
chr7    55174772  rs121913527  T    A    .     .     DP=312;GS=PASS;GD=.;AAD=22;PAD=0;NAD=3;SB_PVAL=2.4000e-01;SB_OR=1.3000;FSB_PVAL=3.1000e-01;FSB_OR=1.1000  GT:DP:AD:ADF:ADR:VAF:FAD:FADF:FADR:FAF:AAD:PAD:NAD  0/1:312:290,22:145,10:145,12:0.0705:145,5:72,5:73,6:0.0735:22:0:3 # (1)!
```

1. `DP=312` total reads; `GS=PASS` normalization status; `GD=.` no diagnostic flags; `AAD=22` reads with any ALT evidence; `PAD=0` no partial matches (SNP — always 0); `NAD=3` reads with N at variant position. FORMAT `DP=312` total depth (single int). `AD=290,22` → 290 REF + 22 ALT reads. `ADF=145,10` → forward strand. `ADR=145,12` → reverse strand. `VAF=0.0705` (read level). `FAD=145,5` → fragment counts. `FAF=0.0735` (fragment level).

---

## MAF Output (`--format maf`)

A tab-separated file following GDC MAF conventions. One row per variant per sample.

### Two Output Paths

The set of columns in the first row of the header depends on whether the
**input** was a VCF or a MAF.

=== "VCF → MAF"

    gbcms generates a GDC-compatible MAF from scratch, since VCF records
    have no MAF metadata. The following **fixed** headers are always present:

    | Column | Description |
    |:-------|:------------|
    | `Hugo_Symbol` | Empty — not populated from VCF input |
    | `Chromosome` | Chromosome name |
    | `Start_Position` | 1-based MAF start position |
    | `End_Position` | 1-based MAF end position |
    | `Strand` | `+` |
    | `Variant_Classification` | Derived from variant type |
    | `Variant_Type` | `SNP`, `INS`, `DEL`, or `ONP` |
    | `Reference_Allele` | MAF-style REF (`-` for pure insertions) |
    | `Tumor_Seq_Allele1` | Reference allele (same as `Reference_Allele`) |
    | `Tumor_Seq_Allele2` | MAF-style ALT (`-` for pure deletions) |
    | `Tumor_Sample_Barcode` | BAM sample name (from `--bam name:path`) |
    | `Matched_Norm_Sample_Barcode` | Empty |
    | `vcf_id` | Original VCF `ID` field (rsID or `.`) |
    | `vcf_pos` | Original VCF 1-based `POS` |
    | `vcf_region` | `chr:pos` tracking field |

    Then all [gbcms count columns](#gbcms-count-columns) are appended.

=== "MAF → MAF"

    All original input MAF columns are preserved **exactly** (values never
    overwritten, column order never changed). gbcms count columns are
    appended after the last original column.

    !!! important "Column Pass-Through Guarantee"
        Every column in your input MAF — including custom lab-specific columns
        like `patient_id`, `assay_version`, pipeline provenance fields, etc. —
        appears unchanged in the output. Only new gbcms columns are added.

---

### `Tumor_Sample_Barcode` Behaviour

!!! caution "rsIDs in Tumor_Sample_Barcode?"
    If you see rsIDs (e.g. `rs121913527`) in `Tumor_Sample_Barcode`, the
    likely cause is that your **input MAF** already has rsIDs in that column
    **and you ran with `--preserve-barcode`**. The fix is either to not use
    `--preserve-barcode`, or to pre-clean the input MAF.

| Input | `--preserve-barcode` | `Tumor_Sample_Barcode` value |
|:------|:--------------------:|:-----------------------------|
| VCF → MAF | any | BAM `sample_name` (always — VCF has no barcode) |
| MAF → MAF | `false` (default) | BAM `sample_name` overwrites original |
| MAF → MAF | `true` | **Original value** from input MAF row |

---

### gbcms Count Columns

These columns are **always** appended regardless of input format.

=== "Default (no prefix)"

    | Column | Type | Description |
    |:-------|:-----|:------------|
    | `gbcms_status` | String | Normalization/counting status. Semicolon-separated multi-value. First token is always `PASS` or `FAIL_*`. Examples: `PASS`, `PASS;WARN_REF_CORRECTED`, `FAIL_REF_MISMATCH`. |
    | `gbcms_diagnostic` | String | Post-counting diagnostic flags. Semicolon-separated. Empty string when no diagnostics. Examples: `ZERO_ALT`, `PARTIAL_DOMINANT;MNP_DISC_RATIO(2/5);MNP_RESCUE_ELIGIBLE`. |
    | `gbcms_rescue` | String | **Conditional** — only present when `--rescue-mnp` is enabled. Structured audit trail for MNP decomposition rescue. Format: `method=decomposed;original_alt=0;positions=chr:pos(R>A):count,...`. Empty when no rescue was attempted. Failed rescues include `outcome=no_signal`. |
    | `ref_count` | Integer | REF read depth |
    | `alt_count` | Integer | ALT read depth |
    | `any_alt` | Integer | Any ALT Depth — reads with ALT evidence at ≥1 discriminating position. Invariant: `any_alt = alt_count + partial_alt` |
    | `partial_alt` | Integer | Partial ALT Depth — reads matching ALT at some but not all discriminating positions. Populated for all variant types including INDELs (via Phase 3 structural evidence propagation). |
    | `n_count` | Integer | N-base Depth — reads with N base at ≥1 discriminating position (duplex masking QC metric) |
    | `total_count` | Integer | Total read depth (DP) |
    | `vaf` | Float | Read-level variant allele fraction |
    | `ref_count_forward` | Integer | REF reads on forward strand |
    | `ref_count_reverse` | Integer | REF reads on reverse strand |
    | `alt_count_forward` | Integer | ALT reads on forward strand |
    | `alt_count_reverse` | Integer | ALT reads on reverse strand |
    | `strand_bias_p_value` | Float | Fisher's exact test p-value (read-level) |
    | `strand_bias_odds_ratio` | Float | Fisher's exact test odds ratio (read-level) |
    | `ref_count_fragment` | Integer | REF fragment count |
    | `alt_count_fragment` | Integer | ALT fragment count |
    | `total_count_fragment` | Integer | Total fragment count |
    | `vaf_fragment` | Float | Fragment-level variant allele fraction |
    | `ref_count_fragment_forward` | Integer | REF fragments on forward strand |
    | `ref_count_fragment_reverse` | Integer | REF fragments on reverse strand |
    | `alt_count_fragment_forward` | Integer | ALT fragments on forward strand |
    | `alt_count_fragment_reverse` | Integer | ALT fragments on reverse strand |
    | `fragment_strand_bias_p_value` | Float | Fragment-level strand bias p-value |
    | `fragment_strand_bias_odds_ratio` | Float | Fragment-level strand bias odds ratio |

=== "With `--column-prefix t_`"

    All count columns above (except `gbcms_status`, `gbcms_diagnostic`, `gbcms_rescue`, and strand bias) are
    prefixed with `t_`:

    | Column | Example |
    |:-------|:--------|
    | `t_ref_count` | `80` |
    | `t_alt_count` | `20` |
    | `t_total_count` | `100` |
    | `t_vaf` | `0.2000` |
    | `t_ref_count_fragment` | `45` |
    | ... | ... |

    Use `--column-prefix t_` for downstream tools that expect the legacy
    `t_ref_count` / `t_alt_count` column naming.

=== "With custom prefix"

    Any prefix matching `[A-Za-z0-9_]` is accepted:

    ```bash
    gbcms dna --column-prefix plasma_ ...
    # → plasma_ref_count, plasma_alt_count, plasma_total_count, ...
    ```

    !!! note
        `gbcms_status`, `gbcms_diagnostic`, `gbcms_rescue`, and the four `strand_bias_*` columns are **never
        prefixed** — they are always unique even when count columns share a prefix.

---

### RNA-Specific MAF Columns

!!! info "RNA mode only"
    These 5 columns are appended **only** when using `gbcms rna`. They do not
    appear at all in DNA mode output.

| Column | Type | Description |
|:-------|:-----|:------------|
| `rna_sense_depth` | Integer | Total reads on the transcript sense strand at this position |
| `rna_antisense_depth` | Integer | Total reads on the antisense strand |
| `rna_alt_sense_count` | Integer | ALT reads on the sense strand |
| `rna_editing_site` | Boolean | `True` if the locus overlaps a known A-to-I editing site (requires `--rna-editing-db`) |
| `rna_splice_spanning` | Integer | ALT reads whose alignment spans a splice junction (`N` CIGAR operation) |

---

### GTF-Aware MAF Columns (v5.0.0)

!!! info "RNA mode + `--gtf` only"
    These columns are appended **only** when using `gbcms rna --gtf <file>`. They are
    completely absent without the `--gtf` flag — no empty/NA placeholders.

#### Exon Boundary Distance

| Column | Type | Description |
|:-------|:-----|:------------|
| `exon_boundary_dist` | Integer | Signed distance to the nearest exon boundary. Positive = exonic (distance from exon edge inward); negative = intronic (distance from nearest exon edge outward). `0` = exactly at an exon boundary. |

#### Per-Transcript Counts

| Column | Type | Description |
|:-------|:-----|:------------|
| `transcript_read_counts` | String | Semicolon-separated per-transcript read-level count triplets. Format: `ENST...:AD,RD,DP;ENST...:AD,RD,DP`. Example: `ENST00000269305:11,140,162;ENST00000445888:7,95,108`. Empty when no GTF or no overlapping transcripts. |
| `transcript_fragment_counts` | String | Same format as `transcript_read_counts` but with fragment-level counts: `ENST...:ADF,RDF,DPF`. Fragment counts ≤ read counts for each transcript. |

#### Aberrant Splice Junction Detection (ASJD)

| Column | Type | Description |
|:-------|:-----|:------------|
| `asjd_flag` | Boolean | `True` when allele-specific junction divergence is detected (Fisher p < 0.05) |
| `asjd_pval` | Float | Raw Fisher exact test p-value comparing REF vs ALT junction usage |
| `asjd_qval` | Float | Benjamini-Hochberg corrected q-value (FDR control across all variants) |
| `asjd_ref_junction` | String | Dominant REF junction coordinates (`start-end`), empty if no junction |
| `asjd_alt_junction` | String | Dominant ALT junction coordinates (`start-end`), empty if no junction |
| `asjd_ref_motif` | String | Splice motif at REF junction: `GT-AG`, `GC-AG`, `AT-AC`, `OTHER`, or `UNKNOWN` |
| `asjd_alt_motif` | String | Splice motif at ALT junction (same categories) |
| `asjd_ref_known` | Boolean | `True` if the REF dominant junction matches a GTF-annotated intron |
| `asjd_alt_known` | Boolean | `True` if the ALT dominant junction matches a GTF-annotated intron |
| `asjd_n_ref_junc` | Integer | REF reads on the dominant junction |
| `asjd_n_alt_junc` | Integer | ALT reads on the dominant junction |
| `asjd_n_ref_total` | Integer | Total REF reads with any splice junction |
| `asjd_n_alt_total` | Integer | Total ALT reads with any splice junction |
| `asjd_diagnostic` | String | Semicolon-separated QC flags (see [Diagnostic Flags](#asjd-diagnostic-flags)) |

##### ASJD Diagnostic Flags

| Flag | Condition | Meaning |
|:-----|:----------|:--------|
| `LOW_ALT_JUNC` | `asjd_n_alt_junc < 5` | Insufficient ALT junction evidence |
| `LOW_REF_JUNC` | `asjd_n_ref_junc < 10` | Insufficient REF baseline |
| `NOVEL_ALT_JUNC` | `asjd_alt_known == false` | ALT uses unannotated junction |
| `NON_CANONICAL_MOTIF` | ALT motif not GT-AG/GC-AG/AT-AC | Likely mapping artifact |
| `STRAND_DISCORDANT` | ALT junction minority strand ≥ 30% | dUTP artifact |
| `MULTI_JUNCTION` | ALT reads use > 2 junctions | Complex splicing event |

---

### Library Type Behavioral Note (v5.0.0)

!!! warning "Amplicon Mode"
    When `--library-type amplicon` is used, fragment counts (`dpf`, `rdf`, `adf`,
    `ref_count_fragment`, `alt_count_fragment`) will approximate read counts (`dp`, `rd`, `ad`,
    `ref_count`, `alt_count`). This is expected — amplicon mode bypasses R1/R2 fragment
    consensus merging, treating each read as an independent observation.

    This does **not** affect DNA mode output — `library_type` is an RNA-only parameter.

??? note "mFSD MAF Columns (`--mfsd` only)"

    34 columns are appended when `--mfsd` is set. They are completely absent
    without the flag (not NA-filled):

    | Column | Type | Description |
    |:-------|:-----|:------------|
    | `mfsd_ref_count` | Integer | REF-classified fragments in 50–1000 bp window |
    | `mfsd_alt_count` | Integer | ALT-classified fragments |
    | `mfsd_nonref_count` | Integer | Non-REF, non-ALT fragments |
    | `mfsd_n_count` | Integer | Fragments with no valid insert size |
    | `mfsd_alt_llr` | Float | Log-likelihood ratio (ALT fragments; positive = tumor-like) |
    | `mfsd_ref_llr` | Float | Log-likelihood ratio (REF fragments) |
    | `mfsd_ref_mean` | Float | Mean fragment size for REF class (bp) |
    | `mfsd_alt_mean` | Float | Mean fragment size for ALT class (bp) |
    | `mfsd_nonref_mean` | Float | Mean fragment size for non-REF class (bp) |
    | `mfsd_n_mean` | Float | Mean fragment size for N class (bp) |
    | `mfsd_delta_alt_ref` | Float | mean(ALT) − mean(REF) delta (bp) |
    | `mfsd_ks_alt_ref` | Float | KS D-stat (ALT vs REF) |
    | `mfsd_pval_alt_ref` | Float | KS p-value (ALT vs REF) |
    | `mfsd_delta_alt_nonref` | Float | mean(ALT) − mean(non-REF) delta |
    | `mfsd_ks_alt_nonref` | Float | KS D-stat (ALT vs non-REF) |
    | `mfsd_pval_alt_nonref` | Float | KS p-value |
    | `mfsd_delta_ref_nonref` | Float | mean(REF) − mean(non-REF) delta |
    | `mfsd_ks_ref_nonref` | Float | KS D-stat |
    | `mfsd_pval_ref_nonref` | Float | KS p-value |
    | `mfsd_delta_alt_n` | Float | mean(ALT) − mean(N) delta |
    | `mfsd_ks_alt_n` | Float | KS D-stat |
    | `mfsd_pval_alt_n` | Float | KS p-value |
    | `mfsd_delta_ref_n` | Float | mean(REF) − mean(N) delta |
    | `mfsd_ks_ref_n` | Float | KS D-stat |
    | `mfsd_pval_ref_n` | Float | KS p-value |
    | `mfsd_delta_nonref_n` | Float | mean(non-REF) − mean(N) delta |
    | `mfsd_ks_nonref_n` | Float | KS D-stat |
    | `mfsd_pval_nonref_n` | Float | KS p-value |
    | `mfsd_error_rate` | Float | non-REF fraction of valid mFSD fragments |
    | `mfsd_n_rate` | Float | N-class fraction |
    | `mfsd_size_ratio` | Float | mean(ALT) / mean(REF) |
    | `mfsd_quality_score` | Float | 1 − error_rate − n_rate |
    | `mfsd_alt_confidence` | String | `HIGH` (≥5 ALT fragments), `LOW` (1–4), or `NONE` |
    | `mfsd_ks_valid` | Boolean | `True` when both ALT and REF have ≥5 fragments for reliable KS test |

??? note "Normalization MAF Columns (`--show-normalization` only)"

    | Column | Type | Description |
    |:-------|:-----|:------------|
    | `{prefix}norm_Start_Position` | Integer | Left-aligned MAF start position |
    | `{prefix}norm_End_Position` | Integer | Left-aligned MAF end position |
    | `{prefix}norm_Reference_Allele` | String | Left-aligned REF allele |
    | `{prefix}norm_Tumor_Seq_Allele2` | String | Left-aligned ALT allele |

    The `{prefix}` matches `--column-prefix` (default: no prefix).

---

## Merged MAF Output (`gbcms merge`)

When multiple BAM types (e.g., duplex, simplex) are genotyped separately and merged
via `gbcms merge`, the output MAF contains type-prefixed columns plus optional
combined metrics.

### Type-Prefixed Columns

Each input MAF's gbcms count columns are prefixed with the BAM type label:

| Input Label | Example Columns |
|:------------|:---------------|
| `duplex` | `duplex_ref_count`, `duplex_alt_count`, `duplex_vaf`, ... |
| `simplex` | `simplex_ref_count`, `simplex_alt_count`, `simplex_vaf`, ... |

Annotation columns (e.g., `Hugo_Symbol`, `Chromosome`) are taken from the first
input and **not** duplicated.

### Combined `simplex_duplex_*` Columns

When both `simplex` and `duplex` inputs are present (and `--no-combined` is not
set), 20 combined columns are appended. Duplex and simplex consensus molecules
are distinct — counts are **additive** across BAM types with no double-counting.

| Phase | Columns | Count | Method |
|:------|:--------|:------|:-------|
| **Additive sums** | Read counts, strand counts, fragment counts, fragment strand counts | 12 | `simplex_{x} + duplex_{x}` |
| **Derived totals** | `total_count`, `total_count_fragment` | 2 | `ref + alt` |
| **Derived VAFs** | `vaf`, `vaf_fragment` | 2 | `alt / total` (0/0 → 0.0) |
| **Strand bias** | `strand_bias_p_value`, `strand_bias_odds_ratio`, `fragment_strand_bias_p_value`, `fragment_strand_bias_odds_ratio` | 4 | Rust Fisher exact 2×2 test |

!!! note "Schema-Aware"
    If strand-level columns are absent from the input MAFs (e.g., older gbcms versions),
    only the available metrics are computed. Missing columns are logged and skipped —
    the pipeline does not fail.

---

## Per-Sample File Naming

```
{--output-dir}/{sample_name}{--suffix}.{vcf|maf}
```

| Component | Source |
|:----------|:-------|
| `sample_name` | `name` from `--bam name:path`; falls back to BAM filename stem |
| `--suffix` | Literal string appended before the extension (e.g. `.genotyped`) |
| Extension | `vcf` or `maf` depending on `--format` |

**Examples:**

```bash
--bam tumor:tumor.bam --suffix .fillout --format maf
# → tumor.fillout.maf

--bam tumor.bam --format vcf
# → tumor.vcf  (stem = "tumor")
```

??? note "Companion Parquet file (`--mfsd-parquet`)"
    When `--mfsd-parquet` is also set (alongside `--mfsd`), a second file is
    written alongside the main output:

    ```
    {--output-dir}/{sample_name}{--suffix}.fsd.parquet
    ```

    It contains per-variant raw fragment size arrays (`ref_sizes`, `alt_sizes`)
    for downstream mFSD visualisations (density plots, empirical CDF comparisons).
    Written natively by Rust — no pyarrow dependency required.

---

## Missing Values

| Format | Missing value sentinel |
|:-------|:----------------------|
| MAF columns | `NA` |
| VCF INFO numeric fields | `.` (VCF spec) |

A value is `NA`/`.` when the count supporting it is zero (e.g. `mfsd_alt_mean`
when no ALT fragments were observed) or when the input variant was rejected
during preparation (all counts are zero-filled in that case).

---

## Related

- [Input Formats](input-formats.md) — VCF and MAF input requirements
- [Counting & Metrics](counting-metrics.md) — How counts are computed from reads
- [gbcms dna](../cli/dna.md) — DNA mode CLI reference
- [gbcms rna](../cli/rna.md) — RNA mode CLI reference
- [Variant Normalization](variant-normalization.md) — How variants are prepared before counting
- [Allele Classification](allele-classification.md) — How each read is classified as REF/ALT/neither
