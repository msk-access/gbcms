# RNA Splice-Junction Handling

## The Problem: Alignment Artifacts at Exon Boundaries

RNA-seq aligners (STAR, HISAT2) produce splice-aware alignments with
CIGAR `N` (RefSkip) operations spanning introns. Bases near these splice
junctions may be misaligned — the aligner sometimes extends a few bases
into the intron ("overhang"), creating false variant calls at exon
boundaries.

```
Reference:   ...EXON1=====][INTRON...intron...INTRON][=====EXON2...
Read:        ...EXON1=====XX                        YYYY=====EXON2...
                          ^^                        ^^^^
                     overhang (2bp)             overhang (4bp)
                     may be misaligned          may be misaligned
```

## The Evidence Rule: What a RefSkip Means

Before any artifact handling, gbcms applies one semantic rule to every
spliced read: **a read testifies about a variant only through aligned
bases (or a `D` op) at the discriminating positions.** A CIGAR `N`
asserts spliced-out reference — unlike `D`, it is not a claim that the
molecule lacks those bases; it is a claim that the read observes
nothing there (samtools pileup reports zero coverage inside an N gap).

Concretely (`splice_skip_triage` in `variant_checks.rs`, applied before
per-type classification):

- A read whose `N` spans **every** discriminating position — the
  deleted span of a deletion, both junction flanks of an insertion,
  every REF base of an SNV/MNP/complex — classifies **neither** and is
  **excluded from DP and fragment depth entirely**. It is not "neither
  at the locus"; it is not at the locus. An intronic position inside a
  spliced-out intron therefore gets depth only from pre-mRNA reads,
  matching pileup. (At an anchor-preserved deletion this is deliberately
  stricter than pileup depth at POS — the anchor base may be aligned,
  but the event is still unobserved, and an unobservant read in DP only
  deflates VAF.) One caveat this creates: RNA aligners represent large
  deletions as splices (STAR writes any deletion ≥ `alignIntronMin`,
  default 21bp, as `N`), so genuine large-deletion carriers can arrive
  N-represented and be excluded. Deletion loci where such exclusions
  exceed confirmed ALT are flagged `SPLICE_SKIP_DOMINANT(n)` in
  `gbcms_diagnostic` — inspect them in IGV before trusting `AD=0`.
- The excluded population is not discarded silently. With `--gtf`, ASJD
  reads it at splice-site variants and reports what it shows in
  `asjd_diagnostic`: `RETENTION_DOMINANT(n)` when spliced-over fragments
  dominate a junction-free classified population (the counts then describe
  intron-retaining reads only — `vaf` is the retention-population VAF), and
  `NOVEL_JUNC_AT_SPLICE_LOSS(n@start-end)` when an anchored, unannotated
  junction on those reads outnumbers confirmed ALT (the mutant allele's
  exon-skip or alternative-site outcome). See
  [RNA Annotation → Diagnostic Flags](rna-annotation.md#diagnostic-flags).
- A `D` op remains deletion evidence; an `N` op never is. The same
  100bp gap counts ALT when the aligner writes `D(100)` and counts
  nothing when it writes `N(100)` — the call must not flip on the
  aligner's representation choice. (Pre-rule, the N form fell through
  the anchor fast path as definitive REF.)
- An indel op **directly after** a splice `N` (`M-N-D-M`, `M-N-I-M` —
  an event at an exon boundary reached through the junction) gets the
  same anchor/windowed inspection as an op after an `M` block. The
  anchor base itself may be spliced out; the evidence is attributed to
  the nearest aligned or inserted base.
- **Span-aligned REF testimony**: at a pure-deletion locus whose anchor
  base is spliced out, a junction read whose aligned bases cover the
  **entire deleted span** demonstrates the deletion is absent — it counts
  REF (the span, not the anchor, is the discriminating fact; presence,
  not per-base identity, same coverage standard as the anchor fast
  path). Partial span coverage stays neither, and a read carrying a
  competing indel candidate in the scan window keeps its arbitration
  path. Typical population: deletion annotations left-aligned to the
  last intronic base at an acceptor — at a validated FORTE locus this
  reclassifies 823 of 824 anchor-spliced junction reads to REF (the one
  residual covers only part of the span and stays neither).
- Phase 3 never scores across a splice: `extract_raw_read_window`
  refuses windows that an `N` overlaps (a contiguous slice would stitch
  the exon arms into a junction-chimeric sequence in which the missing
  intron reads as deletion evidence), and `check_complex`'s
  reconstruction classifies such reads neither instead of
  string-comparing them.

Reads without `N` ops never enter this triage — DNA-mode classification
is untouched. Per-variant exclusion counts are logged at debug level in
the `Phase stats` line (`splice_skip_excluded=`).

## Community Approach: GATK SplitNCigarReads

GATK addresses this by **physically modifying the BAM** before variant
calling:

1. **Splits reads** at each CIGAR `N` boundary into separate alignment
   records (e.g., `45M 5000N 55M` → Fragment 1 `45M` + Fragment 2
   `55M`).
2. **Hard-clips overhangs** if they contain >1 mismatch within a ≤40bp
   window.
3. **Reassigns MAPQ** 255 → 60 for downstream caller compatibility.

This produces a clean BAM where no read contains `N` operators, at the
cost of losing the original read structure.

## gbcms Approach (No BAM Modification)

gbcms handles splice-junction artifacts **without modifying the BAM**:
the evidence rule above governs what a spliced read may testify to, and
heuristic BAQ handles overhang misalignment.

!!! note "Removed: consensus intron snipping of `ref_context`"
    An earlier step drained consensus introns from `ref_context` in place
    so Phase 3 could score junction reads against a mature-mRNA haplotype.
    It had no coordinate map — the context shrank while
    `ref_context_start` stayed genomic — so every position-indexed check
    right of a snipped intron (shifted-deletion sequence verification, the
    large-deletion band's context guard, haplotype offsets) read garbage,
    so exon-contained reads were scored against haplotypes whose offsets
    no longer matched their genomic coordinates, and pre-mRNA /
    intron-retention reads against a haplotype missing bases their
    sequence genuinely contains. It was removed; `ref_context` is always
    genomic.
    Splice-aware Phase-3 scoring, if real-data measurement shows it is
    needed, requires an explicit genomic→spliced coordinate map with
    junction-compatible extraction (tracked in issue #94).

### Heuristic BAQ Near Splice Junctions (`--apply-baq`)

`apply_heuristic_baq()` in `baq.rs` penalizes base qualities within
5bp of CIGAR `N`/RefSkip operations:

- Subtracts 20 from BQ at bases near splice boundaries
- Bases that drop below `--min-baseq` (default 20) are effectively
  filtered from variant evidence
- This mimics GATK SplitNCigarReads' overhang clipping without
  removing bases — a softer approach that preserves information

**Exon-boundary exception (with `--gtf`).** At a variant within 5bp of an
annotated exon boundary (`exon_boundary_dist`), the penalty would land on
exactly the reads that splice there, which are the evidence at an exon edge.
BAQ is therefore skipped for that variant. The main counts, the per-transcript
counts and ASJD all use this one rule, so at an exon edge a transcript's
counts match the variant's counts, as they do elsewhere. Away from a boundary,
BAQ applies in every view. A DEBUG line names each variant form where it is
skipped.

**Default behavior:**

| Mode | `--apply-baq` default | Rationale |
|:-----|:---------------------|:----------|
| `gbcms dna` | `False` (off) | DNA BAMs typically go through BQSR or consensus calling (fgbio/Marianas), which already recalibrates BQ |
| `gbcms rna` | `True` (on) | RNA BAMs typically do **not** go through BQSR or consensus calling; raw sequencer BQ values are preserved |

## When to Use `--no-baq` in RNA Mode

Disable BAQ if your upstream workflow already handles splice-boundary
quality adjustment:

- **SplitNCigarReads**: If you run GATK SplitNCigarReads before gbcms,
  the overhang problem is already resolved. Use `--no-baq` to avoid
  double-penalization.
- **BQSR**: If your pipeline runs GATK Base Quality Score Recalibration,
  base qualities are already machine-learning-calibrated. BAQ on top
  may over-penalize.
- **Consensus-based deduplication**: If your pipeline uses fgbio or
  similar consensus calling (not typical for RNA), BQ values reflect
  consensus agreement. BAQ is unnecessary.

## When BAQ Is Most Valuable

BAQ provides the most benefit when the upstream workflow does **not**
run:

- GATK SplitNCigarReads (overhang clipping)
- GATK BQSR (quality recalibration)
- Consensus-based deduplication (e.g., fgbio)

In these pipelines, raw sequencer BQ values are preserved all the way
to gbcms. A base near a splice boundary may have BQ 30+ even if it is
a misalignment artifact. The -20 BAQ penalty brings it below the
`min_baseq=20` threshold, effectively filtering it from variant
evidence.

## Technical Details

### BAQ Parameters

| Parameter | Value | Description |
|:----------|:------|:------------|
| `BAQ_RADIUS` | 5 bp | Bases within 5bp of indel/splice boundary are penalized |
| `BAQ_PENALTY` | 20 | Subtracted from BQ (clamped to 0, never negative) |

### Performance

BAQ is applied lazily — `apply_heuristic_baq()` returns `None` for
reads without indels or splice junctions (no allocation). For RNA BAMs
where most reads have splice junctions, the cost is one CIGAR walk per
read plus a quality vector copy, which is negligible compared to PairHMM
alignment.

### Trace Logging

When `--trace` is enabled, BAQ logs per-read adjustments:

```
TRACE BAQ: adjusting read A00000:123:AAAA:1:1:1234:5678 (0 indels, 2 splice junctions)
```

---

## Defense in Depth: Splice Bleed Protection

"Splice bleed" occurs when a read has only 2–3 bases on one side of a
splice junction. The aligner may force those bases to align to the
wrong exon, creating false variant evidence at exon boundaries.

### Worst-Case Example

A known SNP (G→A) sits at the first base of Exon 2:

```
                   Exon boundary
                       │
                       │  Known SNP: G→A
                       │  ↓
Genome:    ...INTRON────┤──[G]CGTACGTACGT...
                        │   Exon 2

Read A:             ════╗══[A]═══════════════     → ALT ✅ real variant
                        ╚                         (CIGAR: 10M 3000N 40M)

Read B:    ════════════─╗──[A]──                  → ALT 🟡 splice bleed!
                        ╚                         (CIGAR: 30M 3000N 3M)
                                 ↑
                                 Only 3 bases in Exon 2.
                                 The [A] may be from intronic
                                 sequence, not a real ALT.

Read C:    ═════════════╗══[G]═══════════════     → REF ✅
                        ╚                         (CIGAR: 20M 3000N 30M)
```

### How gbcms Catches Read B: 5 Defense Layers

For a splice-bled base to produce a false positive, **all 5 layers
must fail simultaneously** — which is exceedingly rare:

```
┌────────────────────────────────────────────────────────────────────┐
│                       Defense Layers                               │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  Layer 1: MAPQ filter (--min-mapq, default 1 for RNA)             │
│  └── STAR gives MAPQ=0 to reads it can't confidently place.       │
│      Reads with very short exonic overlap are often multi-mapped.  │
│                                                                    │
│  Layer 2: Base quality filter (--min-baseq, default 20)           │
│  └── Bled bases typically receive low BQ from the aligner because  │
│      the alignment evidence is ambiguous near the splice site.     │
│                                                                    │
│  Layer 3: PairHMM quality-weighted classification                 │
│  └── Even if a bled base passes the BQ filter, PairHMM evaluates  │
│      the FULL read against REF/ALT haplotypes. The 2–3 uncertain  │
│      bases are outweighed by 30+ properly aligned bases on the     │
│      other side. LLR typically falls below threshold → NEITHER.    │
│                                                                    │
│  Layer 4: splice_spanning_count (tracked in output)               │
│  └── Downstream QC can flag variants where ALL ALT evidence comes  │
│      from splice-spanning reads (rna_splice_spanning in MAF/VCF).  │
│                                                                    │
│  Layer 5: BAQ splice-proximity penalty (--apply-baq, on for RNA)  │
│  └── Subtracts 20 from BQ within 5bp of splice junctions. A base  │
│      with BQ 30 near a splice site becomes BQ 10, which is below  │
│      --min-baseq 20 → excluded from allele evidence entirely.      │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

### Why All 5 Layers Must Fail

For a splice-bled base to cause a false ALT call:

1. The variant must sit within **1–2 bp of the exon boundary**
2. The aligner must give the bled base **MAPQ ≥ 1** (Layer 1 pass)
3. The aligner must give the bled base **BQ ≥ 20** (Layer 2 pass)
4. BAQ must fail to bring BQ below 20 — but BAQ subtracts 20 from
   bases within 5bp of the splice junction, so a base needs **BQ ≥ 40**
   to survive (Layer 5 pass)
5. PairHMM must give the read **LLR ≥ 2.3** despite only 2–3 bases
   of evidence on the bled side (Layer 3 pass)
6. The wrong base must **happen to match** the known ALT allele
   (1/3 chance for SNPs)

In practice, Layers 2 and 5 together require BQ ≥ 40 at the bled
position, which is exceptionally unlikely for a 2–3 base overhang.
