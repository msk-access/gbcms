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
  the nearest aligned or inserted base. A windowed candidate (not at the
  anchor) counts ALT only when it is the read's one change across the
  variant's discrimination window: a `D` written just after an `N` that
  covers part of that window is where the aligner ended the `N` (the read's
  bases are those of a reference read spliced a few bases later), so such a
  read counts toward depth only (it counted ALT from its CIGAR before
  6.6.0).
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
  intron reads as deletion evidence).
- A complex variant's exact-carrier windows for a read spliced in their
  flank are built over the reference spliced at the read's own
  junctions (RJ-19, C32 #213): the read's bases past a splice are the next
  exon's, so both alleles' haplotypes continue there (the far exon read from
  the FASTA), and repeat growth and flank are measured on that spliced
  sequence. A length change that pushes an exon's last bases past the
  junction shows against the next exon instead of fitting a window cut at the
  exon edge (on FORTE splice probes, spurious ALT 19 → 4), and the read must
  hold the spliced flank as any read holds its flank (reads reaching one or
  two bases past the junction hold none: −0.95% of REF at probes 0–1 bp from
  the edge). A read spliced through the bases where the alleles differ shows
  neither and counts toward depth only.
- A splice is not reference coverage (RJ-17, R5 #198): a pure-indel read is
  informative only when one aligned block between splices spans the window,
  REF and ALT alike. A read spliced inside an indel's change interval is
  depth only (the event may sit in the skipped intron), and so is a deletion
  written right after a read's splice (its bases equal REF spliced at an
  acceptor further on).
- An RNA read's soft clip that reaches an exon edge or a junction end
  (annotated, or one the variant's reads splice at), or whose aligned bases end
  within five bases of one, is not allele evidence (RJ-18, C15 #173): STAR clips
  a junction overhang it cannot splice, so such a clip holds the next exon's
  bases. Any other clip is the read's own bases, read as in DNA.
- A junction that starts inside a complex event's differing bases and runs past
  them (or ends inside them) splices the haplotypes at the event's edge (RJ-19),
  so a delins carrier counts whether the aligner wrote `X D N` or `X N`.
- A pure indel's read of another length is another allele in RNA as in DNA, judged by its bases on
  its own spliced reference; homopolymer slippage on the FORTE sample matched IMPACT DNA's rate at
  each run length (about 1% at 6 bases, 5% at 8, 10% at 10), so `PARTIAL_DOMINANT` reads the same
  in both modes (measured 2026-10-27 on the FORTE RNA slices).

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
    #94 closed splice-aware Phase-3 scoring as not demanded (B2b;
    [closing comment](https://github.com/msk-access/gbcms/issues/94#issuecomment-5815889505);
    reopen if truth data produces a delins carrier at a junction; see the next
    note). `repeat_span` stays genomic: it sizes the read window, the bins and the
    CIGAR checks' scan, and is not used by the exact-carrier rule, which grows each
    spliced read's event on its own spliced reference (RJ-19, RJ-27). Widening the
    read window 20 bases changed no count on 1,132 FORTE rows, 112 of whose tracts
    meet an exon edge (measured 2026-10-27 on the FORTE RNA slices).

!!! info "Closed as not demanded: splice-aware Phase-3 scoring (B2b of #94)"
    B2b was the coordinate-mapped spliced-haplotype Phase 3: scoring a read across a
    junction against a haplotype spliced at that junction. It was closed with #94 on
    2026-09-24 ([closing comment](https://github.com/msk-access/gbcms/issues/94#issuecomment-5815889505)).

    - **Demand:** measured demand across the RNA truth cohort was zero. Every covered
      junction-adjacent indel classified correctly through the structural paths (post-N
      inspection, windowed scan, span-aligned REF testimony); no case needed alignment scoring
      across a splice ([comment of 2026-09-22](https://github.com/msk-access/gbcms/issues/94#issuecomment-5779396531)).
    - **Community check** ([comment](https://github.com/msk-access/gbcms/issues/94#issuecomment-5779885342)):
      no mainstream tool scores reads across splice junctions. GATK's RNA best practice
      (SplitNCigarReads plus HaplotypeCaller) cuts reads at every N and calls on the exon
      fragments, from aligned bases only: the destructive equivalent of the splice-aware evidence
      rule above. Long-read RNA pipelines apply the same split-style transformation and still
      document junction-adjacent indels as their weak spot (recall drops near junctions), and
      several simply hard-filter splice-adjacent variants. RNA allele counters do not support
      indels at all. So gbcms is at or beyond community practice for the junction population, and
      B2b would be novel territory, which raises the evidence bar for building it. (References: the
      GATK RNAseq short-variant workflow; Genome Biology 13059-023-02923-y on long-read RNA
      alignment and junction-indel recall; ScanNeo2 #192 on artifact inflation without splice
      handling.)
    - **Design note, if demand appears:** prefer annotation-driven transcript-space projection
      (GTF transcript coordinates, a deterministic genomic↔transcript map, per-read
      isoform-compatibility gating) over read-derived consensus splicing. The removed consensus
      snipping failed on consensus instability and the missing coordinate map; transcript-space
      construction avoids both, and keeps pre-mRNA and intron-retention reads genomically scored.
    - **Reopen if** truth data produces a delins carrier at a junction (the entry in
      `.agents/learnings/REJECTED.md` is REJ-20260924-001).
    - **Later, RJ-19** (C32 #213, 2026-10-04, [Read Judgment](read-judgment.md)) built read-derived
      spliced windows for the exact-carrier rule: a spliced read's windows are made over the
      reference spliced at that read's own junctions. It uses no alignment scoring, so it is not
      B2b, and it came from synthetic splice probes (spurious ALT at exon-edge events, 19 → 4), not
      from a truth-data case, so B2b stays closed. Its junctions come from the read, not the
      annotation; the design note's preference for annotation-driven projection concerns any
      future Phase-3 scoring across a splice.

### Heuristic BAQ Near Splice Junctions (`--apply-baq`)

`apply_heuristic_baq()` in `baq.rs` penalizes base qualities within
5bp of CIGAR `N`/RefSkip operations:

- Subtracts 20 from BQ at bases near splice boundaries
- Bases that drop below `--min-baseq` (default 20) are effectively
  filtered from variant evidence
- This mimics GATK SplitNCigarReads' overhang clipping without
  removing bases — a softer approach that preserves information

**Exon-boundary exception (with `--gtf`).** At a variant whose REF span comes
within 5bp of an annotated exon boundary (`exon_boundary_dist`), the penalty
would land on exactly the reads that splice there, which are the evidence at an
exon edge. BAQ is therefore skipped for that variant. The distance is measured
from the whole REF span, so an MNP or deletion reaching into an exon's last
five bases is covered even when its first base is farther out. The main counts,
the per-transcript counts, ASJD and each `--rescue-mnp` component re-count all
use this one rule, so at an exon edge a transcript's counts match the variant's
counts, as they do elsewhere. Away from a boundary, BAQ applies in every view.
A DEBUG line names each variant form where it is skipped.

!!! info "Why the exception reads the whole REF span (measured for #106, R1)"
    Before 6.6.0 the distance was measured from the variant's first base, so a multi-base
    variant reaching an exon's right edge kept BAQ and lost its edge-side bases in every spliced
    read. Local data, aggregates only; the numbers are from
    [the measurement](https://github.com/msk-access/gbcms/issues/106#issuecomment-5816857647),
    [the evidence for the decision](https://github.com/msk-access/gbcms/issues/106#issuecomment-5834304141)
    and [the acceptance](https://github.com/msk-access/gbcms/issues/106#issuecomment-5896307927);
    the rounds are in the archived plans, `.agents/archive/CYCLE_6.5.0_PLAN.md` § T9 and
    `.agents/archive/CYCLE_6.6.0_PLAN.md` § R1.

    - **Demand.** The engine's own preparation over 74,617 unique signed-out multi-base events,
      against GENCODE exon boundaries from all transcripts: the geometry (first base more than 5bp
      from every site, REF span within 5bp) is present in 6,130 of 61,123 deletions (10%; 4,887
      cross an exon boundary, mostly >20bp) and 99 of 13,359 MNPs (0.7%; nearly all end 1-5bp
      short of a right edge).
    - **Effect.** An A/B of develop against a span-aware prototype on 3 junction-rich RNA samples,
      326 probes (one per junction, at least 300bp apart so length-changing probes cannot form
      clusters). Deletions: no effect, 0 of 510 probe-sample rows moved apart from a single
      fragment count (deletion calls are structural). MNPs ending 4-5bp short of a right edge (DNP,
      TNP, 4bp): 139 of 183 rows changed, RD +0.03-0.17% (at most 1% per probe), per-transcript and
      ASJD moving alike, AD and DP unchanged.
    - **Why counts barely move, and what does matter.** Counted exactly with pysam, a median
      90-92% of the reads covering such an MNP have a splice N within BAQ's 5bp radius of an MNP
      base. REF and ALT calls vote on the unmasked bases, so counts hold. But an MNP read counts as
      confirmed (`mnp_confirmed_alt`) only when no base is masked, so about 9 in 10 real carriers
      there could never be confirmed, and `--rescue-mnp`'s `haplotype_confirmed` safeguard was
      blind at right exon edges in RNA.
    - **The column follows the rule.** The 6.5.0 measurement first recommended keeping
      `exon_boundary_dist` POS-based (a span-based column would change about 10% of edge deletions
      for no counting benefit). The evidence of 2026-09-25 reversed it: among 97,246 unique
      signed-out multi-base and indel variants, the span-based distance differs from the POS-based
      one for 41%, 6.1% flip in or out of the 5bp window (98 of them MNPs), and on the FORTE truth
      cohort 26 of the 94 rows are multi-base or indels, the distance differs for 8 and the window
      flips for 1. Operator decision, 2026-09-25: span distance for both.
    - **Standards.** No single spec. Ensembl VEP applies its splice windows by overlap (any changed
      base counts); ANNOVAR measures from the nearest end of the span and flags spans that cross
      the boundary; SpliceAI uses POS only. There is no standard for BAQ near junctions: htslib
      skips BAQ for spliced reads, and GATK4 does not use BAQ.
    - **Not adopted: a read-level rule** (no N-penalty at annotated junctions). It would
      re-penalize novel junctions near a variant, which is ASJD's splice-creation signal, restore
      indel penalties at edges, and change the main counts more (REJ-20260925-002 in
      `.agents/learnings/REJECTED.md`).
    - **Acceptance (2026-09-29).** FORTE truth cohort (33 samples, 94 rows): 0 count changes, all
      68 SNV rows byte-identical, the column changes on 8 multi-base rows and 1 row flips the 5bp
      window with unchanged counts. Exon-edge probes (326 × 3 samples): counts move on 139 MNP
      rows, all where the window flips, and no deletion moves; summed RD +489, partial +6, AD −2,
      depth unchanged. In 135 of the 139 moved rows the change equals what the reads' own bases
      predict (read at sequencer quality against emulated BAQ); in the other 4 every re-called read
      was inspected and each new call is what its bases say (mostly a Q24 REF base 4bp from the
      junction that BAQ had lowered to Q4; one Q24 third-allele base BAQ had hidden). Column: 978
      of 978 probe rows equal an independent span distance of the engine-normalized variant.
      `--rescue-mnp` on the probes: 4 of 243 MNP rows change outcome (1 from rescued to
      `haplotype_confirmed`, 3 newly rescued on one-read components). DNA (28 runs, 1,060 rows)
      byte-identical.

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
