# Reading allele evidence from a BAM: caveats

gbcms is a **genotyper, not a caller**. Given an allele, it counts the reads whose
own bases carry it (and REF, and what else is there, as diagnostics). It does not
decide whether a variant exists, and it never rewrites a given allele.

The BAM is not the truth either. Counting from it rests on two assumptions:

1. **A read is placed where it came from.**
2. **Its bases are the molecule's sequence.**

Both fail in known ways. This page lists each way: how it shows in a BAM, what
gbcms does today, and what it should do. Each "should do" is one of:

- **decided**, when the genotyper principle settles it;
- **proposed**, when there is a design;
- **to measure**, when we don't know yet.

Each caveat also records what the community's tools do, and how to see the
caveat in gbcms's output.

!!! info "How we decide what gbcms should do"
    1. **Principle first.** Judge a read by its own bases, not by the aligner's
       placement ([counting the given allele](allele-classification.md#the-exact-carrier-rule)).
       Expose what cannot be judged through diagnostics, never by filtering.
    2. **If the right behaviour isn't clear, measure it.** Use the validation matrix:
       synthetic data, realigned panel DNA, DNA without realignment, RNA, and
       public reference data.
    3. **Survey the community.** Record what callers and genotypers do (GATK,
       Strelka2, VarDict, bcftools, freebayes and others), with sources.
    4. **Record the decision** in the tracking issue before implementing.

## Summary

| Caveat | gbcms now | Should do | Tracking |
|:--|:--|:--|:--|
| Indel placed differently in a repeat, or one event written as several indels | Judged by bases (anchored windows, repeat growth) | Decided (done); test oracle | #171, C11 #159 |
| Allele soft-clipped at a read end | Counted in DNA by every rule that reads clips (complex, MNP, pure indels, deletions of 50 bp or more); an RNA clip at an exon edge or junction end is not evidence | Decided (done); clip-borne ITD insertions open | C15 #173 (closed), C7 #144 |
| Long event split across alignments | The primary's clip is read at the junction | Proposed | C18 #177 |
| Indel realignment (ABRA2) fixes or forces placements | Invariant by design (bases, not placement) | Decided; verify | D5 #155 |
| Aligner settings change clipping; hard clips drop bases | Clips read; hard-clipped bases lost; warned at run start | Decided (done: warn) | O2 #131, D5 #155 |
| Reads from paralogs or segmental duplications | MAPQ ≥ 20; `mq0_count` | Proposed (diagnostic, 6.7.0); validate on paralog-rich loci | O5 #179 |
| Reference bias: ALT reads map worse, clip, or go unmapped | Counts what maps | Proposed (diagnostic in `gbcms_diagnostic`; no filter, no new column; 6.7.0) | O5 #179 |
| Alt contigs and decoys split coverage | Contig naming handled; ALT contigs warned at run start | Decided (done: warn) | O2 #131 |
| Base-quality calibration: binned, BQSR, BAQ | One masking rule; BAQ spares the variant's own indel; run start reports the share below `--min-baseq` | Decided (done: warn above 10%) | O2 #131 |
| Adapter read-through | Clipped at the fragment end in every read | Decided (done) | C17 #176 |
| A mate's soft-clipped 5' end makes the fragment end (TLEN) short | The read's own soft clips past that end are clipped as adapter | Proposed (measure from the mate's unclipped end with MC) | C31 #212 |
| Exact-length insertion with a sequencing error inside; truncated copies | One error inside the insert is tolerated; a truncated copy counts `partial_alt` | Decided (done); the identity band for truncated inserts was not adopted | C6 #143 |
| A REF molecule with one error beside the event reads ALT | Both alleles' windows are read from one anchor | Decided (done, RJ-23) | C33 #214 |
| Stutter in repeats | Counted exactly; named by a diagnostic | Decided (done) | — |
| OxoG / FFPE orientation artifacts | Strand counts only | To decide: optional F1R2/F2R1 counts or an orientation-imbalance diagnostic; measure first on FFPE-rich data (6.7.0) | O6 #180 |
| Unmarked duplicates; UMI consensus N's | Duplicate flag honoured; N masked; run start reports the flagged share | Decided (done: report, no warning) | O2 #131 |
| Overlapping mates | Fragment counts | Decided (done) | — |
| RNA: splicing at exon edges | Windows built over the reference spliced at each read's junctions; a splice is not coverage; ASJD | Decided (done) | C32 #213, R5 #198 |
| RNA: STAR clips hold next-exon bases | A clip reaching an exon edge or junction end is not evidence; other clips are read | Decided (done) | C15 #173 |
| RNA: editing (A-to-I) inside a window | Flag only | Proposed | R3 #178 |
| RNA: allele-specific expression and NMD | RNA VAF reported as is | Decided (interpretation) | — |
| RNA: strandedness | Gating exists; observability open | Decided | R2 #114 |
| RNA: a GTF transcript ID reused on another chromosome or strand | Exons and introns merge by ID | Proposed | C34 #216 |
| Variant not prepared against the reference (direct API use, failed reference fetch) | The alignment fallback's tail reads stale scores (`partial_alt` only) | Open (6.7.0) | C24 #195 |
| Reads ending inside the event; low depth | Depth only until a read holds the shorter allele's window | Decided (done) | C39 #253 |
| Purity, clonality, contamination, CH | Counts, not calls | Decided (interpretation) | — |

## 1. Where the aligner put things, not what the molecule is

### Indel placement in repeats; one event written as several indels
- **In the BAM:** the same haplotype aligned with its indel anywhere along a
  repeat, or a delins written as an insertion plus a deletion.
- **gbcms now:**
    - windows are anchored on flank bases outside the event;
    - the event grows through repeats on either allele, and each window is read
      from both anchors;
    - pure indels use their shift region (C10).

  The placement does not change the call.
- **Should do:** decided and done. Make an independent, position-aware read census
  the test oracle (#171). Tandem duplications longer than the repeat finder's
  motifs remain open (#159).
- **Community:** GATK HaplotypeCaller and Mutect2 reassemble haplotypes and
  realign reads to them; Strelka2 realigns to candidate indels.
- **See it:** `--trace` read calls; `partial_alt`; `COEXISTING_ALLELE`.

### Alleles soft-clipped at read ends
- **In the BAM:** an ALT read whose allele sits near its end is clipped from there,
  while the REF reads beside it align in full. Excluding clipped reads biases VAF
  down.
- **gbcms now:** in DNA, a read whose clipped bases carry the allele is counted
  when the rule that judges the variant reads them, inside a well-defined
  fragment (RJ-13; [carriers in soft-clipped bases](counting-metrics.md)):
    - complex variants and MNPs, and pure indels (those under 50 bp are judged by
      their bases with soft-clipped bases read, clip-admitted where the rule decides
      them), and pure deletions of 50 bp or more;
    - a read clipped before the flank is placed from the other flank (RJ-23).
- **Should do:** decided; done for everything above (C15 #173, closed after an
  audit on 2026-10-23). Open: clipped duplication carriers, the clip-borne ITD
  insertions (C7 #144); split reads stay on #177. In RNA, clipped bases may come
  from the next exon (STAR clips a junction overhang it cannot splice), so an RNA
  read's clip that reaches an exon edge or junction end is not allele evidence;
  other clips are the read's own bases (RJ-18; measured on FORTE splice probes: 3
  of 19 spurious ALT reads removed).
- **Clip-borne ITD carriers, measured** (#92,
  [2026-09-22](https://github.com/msk-access/gbcms/issues/92#issuecomment-5780022052) and
  [2026-09-23](https://github.com/msk-access/gbcms/issues/92#issuecomment-5788174498)):
    - a 30bp ITD with sign-out t_alt=7 had zero I-op reads in the BAM. Its carriers exist only as
      soft clips whose clip start sits 9-47bp downstream of the anchor (ITD alignments diverge at
      the duplication boundary, not at the insertion point). Phase 3 then admitted soft clips only
      within ±1bp of the variant position, so they never reached haplotype scoring: ad=0,
      partial=0, `ZERO_ALT`. The candidate fix widens the clip-admission window for insertions to
      the anchor plus the insert length (a structural bound, not a tuned constant). It changes
      counts and risks clip noise, so it needs its own red-first, measured cycle;
    - prevalence: 12 further ITD loci (15-63bp, t_alt 11-409) are 12/12 I-op-dominant (I carriers
      16-597 against 1-10 clip candidates per locus). A different timepoint of the same patient as
      the clip-only outlier shows 587 I-op carriers: the outlier was a low-carrier sample (t_alt=7)
      where the few molecules happened to land clip-represented. Verdict: the rescue (C7 #144) stays
      filed at low priority, with `CLIP_CANDIDATES(n)` as the observability stopgap. The truth
      batch these come from is in `PRODUCTION_READINESS.md`.
- **Community** (verified in source):
    - GATK HaplotypeCaller and Mutect2 revert clips and realign reads, but only
      with a well-defined fragment;
    - Strelka2 unrolls edge clips;
    - VarDict rescues clip consensus;
    - samtools/bcftools mpileup, bam-readcount and freebayes ignore clipped bases;
    - bcftools annotates soft-clip bias (`SCBZ`).
- **See it:** trace `admitted … by its soft-clipped bases`; `CLIP_CANDIDATES` for
  insertions.

### Long events split across alignments
- **In the BAM:** without realignment, a long-deletion carrier becomes a primary
  alignment (often clipped) plus a supplementary on the far side (`SA` tag).
- **gbcms now:** the primary's clip is read at the junction. Supplementary
  alignments are filtered, since they are the same molecule. Carriers whose far
  side is hard-clipped or too short stay invisible.
- **Should do:** proposed. Count a molecule as ALT when its alignments jump exactly
  across the given breakpoints (primary + `SA`, counted once) (C18 #177).
- **Community:** structural-variant callers (Manta, Delly) use split reads and
  discordant pairs; small-variant callers rely on assembly or realignment. To
  survey in detail.
- **See it:** `ZERO_ALT` at a long deletion with depth.

### Indel realignment (ABRA2) and aligner settings
- **In the BAM:**
    - a realigner converts clips into indels, and can force a placement;
    - aligner settings change clipping (a lower mismatch penalty, `-Y` soft-clipped
      supplementaries);
    - hard clips drop bases entirely; the run start warns when more than 1% of a
      BAM's sampled primary alignments are hard-clipped (none in 10 MSK BAMs).
- **gbcms now:** judgements rest on bases, not placement, so they are invariant by
  design. Hard-clipped bases cannot be read.
- **Should do:** decided. Verify invariance on realigned and non-realigned pairs of
  the same libraries, and on other aligners (D5 #155). Warn when primaries are
  hard-clipped (O2 #131).
- **Community:** GATK4 replaced indel realignment with local assembly.
- **See it:** compare runs across pipelines; `ZERO_ALT` at loci the other
  pipeline counts.

## 2. Real bases, but from somewhere else

### Paralogs, pseudogenes and segmental duplications
- **In the BAM:** reads from a homologous region carry its differences as if they
  were variants.
- **gbcms now:** minimum MAPQ 20; `mq0_count` is reported. There is no paralog
  awareness: gbcms reports what is at the locus.
- **Should do:** proposed (O5 #179, 6.7.0): a diagnostic in `gbcms_diagnostic`,
  never a filter, no new column. The target is settled by measurement first: how
  often ALT reads at low-mappability and segmental-duplication loci are
  paralog-derived, validated on paralog-rich loci and the paired panel/WES loci.
- **Community:**
    - GIAB stratifies benchmarks by low mappability and segmental duplications;
    - callers annotate mapping-quality rank-sum statistics (GATK `MQRankSum`, to
      verify) or bias Z-scores (bcftools, to verify).
- **See it:** `mq0_count`; per-allele MAPQ once O5 lands.

### Reference bias
- **In the BAM:** reads carrying a long indel or a cluster of changes map with
  lower MAPQ, get clipped, or go unmapped. Baits designed on the reference can also
  capture REF molecules better.
- **gbcms now:** counts every read that maps; clipped-carrier admission recovers
  part of the loss.
- **Should do:** proposed. Per variant, compare the MAPQ, clip rate and
  outside-depth rate of ALT-like reads against REF reads (and the MQ0 fraction),
  and flag a clear imbalance in `gbcms_diagnostic`; no filtering, no new column
  (O5 #179). Validate on the paired panel/WES loci (realignment removes most of the
  bias) and on paralog-rich loci.
- **Community:** clip- and MAPQ-bias annotations, as above; graph and
  alt-aware references reduce the bias upstream. To survey.
- **See it:** `alt_dist_end_median` against `ref_dist_end_median`; per-allele
  MAPQ and clip rate once O5 lands.

### Alt contigs and decoys
- **In the BAM:** alt-aware GRCh38 alignment can move reads onto alt haplotype
  contigs, and decoys absorb reads, so coverage at the primary locus drops.
- **gbcms now:** contig naming is reconciled; nothing else.
- **Should do:** done (O2 #131): the run start warns when a BAM's header has ALT
  contigs (`*_alt`, `HLA-*`): without ALT-aware alignment, reads where ALT
  haplotypes exist get low MAPQ and `--min-mapq` drops them. None of 10 sampled
  MSK BAMs (IMPACT, ACCESS, FORTE) has any.
- **Community:** alt-aware aligners (bwa's `.alt` handling) and GATK's GRCh38
  guidance. To survey.

## 3. Bases that aren't the molecule's sequence

### Base-quality calibration
- **In the BAM:** Illumina binning (for example 2/12/23/37 or 2/12/24/40), BQSR
  recalibration, and heuristic BAQ all change which bases pass `--min-baseq`.
- **gbcms now:**
    - one masking rule: bases below min BQ, and N, match anything;
    - BAQ spares the variant's own indel (#166).
- **Should do:** done (O2 #131): the run start reports each BAM's base-quality
  values (binned when 8 or fewer) and the share of sampled bases below
  `--min-baseq`, and warns above 10%. Measured on 10 MSK BAMs: 0.8–2.1% at the
  default Q20; a threshold between bins removes more (FORTE's 2/9/24/40 at Q25:
  4%; IMPACT at Q30: 5–12%; IMPACT at Q40: nearly all).
- **Community:** GATK BQSR; bcftools and samtools BAQ (`-B`/`-E`).
- **See it:** `n_count`; the per-read trace.

### Sequencing errors inside or beside the event
- **In the BAM:** an exact-length insertion read with one sequencing error inside the
  inserted bases; a REF molecule with one clear error just outside the window of the
  allele it is anchored to; a truncated copy of an insertion (the read holds part of
  the insert).
- **gbcms now:**
    - one substitution in a read is an error on REF and ALT reads alike: an
      exact-length insertion with one error inside is judged by its bases, the error
      tolerated where the other allele mismatches the same base, unless it recurs
      inside the event (RJ-26; C6 #143);
    - both alleles' windows are read from one anchor, so an error any reading sees
      counts against the allele it contradicts and a REF molecule with one flank
      error is not credited ALT (RJ-23; C33 #214, a limit of the first exact-carrier
      rule, where the two windows were read from different anchors);
    - a truncated copy counts `partial_alt`, not ALT. An identity band for truncated
      inserts (the truncation rule's 90% identity and non-low-complexity gates) was
      measured and not adopted for the exact-carrier rule (operator, 2026-10-16): it
      added 158 ALT reads on HG002 long insertions with truth concordance unchanged;
      on the panel it added 13 ALT calls the bases oracle agrees with, against 111 it
      does not. The band stays only where the rule cannot judge (no prepared
      reference, the CIGAR path).
- **Measured before the decision** (#92, [2026-09-22](https://github.com/msk-access/gbcms/issues/92#issuecomment-5779396851)
  and [2026-09-23](https://github.com/msk-access/gbcms/issues/92#issuecomment-5787996247)), when the
  question was whether a base-quality-first policy could replace the 90% identity band:
    - IMPACT raw reads, a mismatch-anatomy survey over the insertion loci of 4bp or more (about 900
      carrier reads at 38bp and 19bp loci): the dominant non-exact class was zero-mismatch
      truncations (about 24% of carriers; 219 candidates), pure length containment that no quality
      policy touches. Mismatching candidates were about 2% of carriers, their bases mostly high
      quality (Q20-35) yet every pattern a singleton at a unique position with a unique substitution
      (isolated errors inside long inserts, not a systematic representation shift); the class where
      both policies agree (all-low-BQ mismatches) was 6 reads.
    - ACCESS consensus reads, at the two band-eligible insertion loci of the truth panel (6bp at
      about 261 duplex candidates; 7bp at about 266): exactly three classes, exact (the
      overwhelming majority), length truncations (2-14 per locus) and N-carrying candidates (0-10
      duplex, up to 83 simplex). Same-length sequence mismatches: zero, on both flavors, at both
      loci. That is what the consensus SOPs predict (v2 masks bases under Q30 to N; v1 marks a
      disagreement as N), so a base-quality-first policy has nothing to act on in consensus data.
    - A supporting truth point from the first FLT3-ITD run (36bp ITD): ad=129 against sign-out
      t_alt=128, with the truncation-smear class (13 reads) separated into `partial_alt`.
    - The band was kept then (REJ-20260922-001 in `.agents/learnings/REJECTED.md`); the exact-carrier
      rule later tested it and did not adopt it (REJ-20261016-002).
- **Should do:** decided and done (the band: decided against). See
  [Read Judgment](read-judgment.md).
- **Community:** not surveyed for this caveat.
- **See it:** `partial_alt`; the `rule=` of each `--trace` read call.

### Adapter read-through
- **In the BAM:** when the fragment is shorter than the read, bases past the mate's
  5' end are adapter.
- **gbcms now:** a read ends at its fragment end. Bases past it are hard-clipped
  as the read enters counting, as if trimmed, so they are neither bases nor reach,
  in every rule (RJ-10, C17 #176), when they look like adapter (soft-clipped, or
  at most two aligned). Bases that go on aligning, or an insertion past the
  boundary, are kept: TLEN understates a molecule carrying an ITD.
- **A limit (C31 #212):** the boundary is TLEN, measured between the pair's
  *aligned* 5' ends. When the mate's 5' end is soft-clipped, TLEN understates the
  molecule, and a read's own soft-clipped bases past the boundary can be molecule
  bases that are clipped as adapter. The adversarial review of the rule built a
  failing case: an anchor-changing insertion that BWA-MEM soft-clips near the read
  end, with a mate whose 5' clip holds the same bases; R1 loses its ALT (AD 8 to 4).
  Measured on the group 1 data, when only complex and MNP windows read soft clips:
    - 203 of 1,782 clipped reads with soft-clipped bases past the boundary have a
      mate with a 5' soft clip (about 2,100 bases; for insertions, 71 of 416 reads);
    - no count moved on RC DNA, FORTE or WES;
    - the `CLIP_CANDIDATES` diagnostic can drop for such reads.
- **Should do:** the boundary rule is decided and done. The mate-clip refinement is
  proposed (C31 #212): measure the fragment end from the mate's *unclipped* 5' end
  when the MC tag (mate CIGAR) is present, and keep today's rule without MC. MC is
  a per-record tag, so the count stays binning-invariant (using the mate record
  from the bin would not be). MC is on all WES reads and part of the RC DNA here;
  FORTE and the probes have none. Decide together with C15 #173 and C7 #144, and
  measure on ITD loci first (the FLT3 ITD row is on the watchlist).
- **Community:** GATK hard-clips adapter at the insert-size boundary; fgbio
  ClipBam clips bases past the mate (`--clip-bases-past-mate`), measuring the
  mate's end from its unclipped coordinates. Adapter trimming upstream (fastp,
  cutadapt) is common but not universal.

### Stutter and polymerase slippage in repeats
- **In the BAM:** reads with one repeat unit more or fewer than the molecule.
- **gbcms now:** counted exactly. They are not the given allele, so they become
  `partial_alt`, and `COEXISTING_ALLELE` names a frequent one.
- **Should do:** decided and done (report, don't correct).
- **Measured** (#92, 2026-09-22,
  [comment](https://github.com/msk-access/gbcms/issues/92#issuecomment-5780293921)): on 6 MSI-high
  samples (MSI scores 45-53) the wrong-length rule already separates real stutter. Median
  ad-versus-sign-out relative error is +1.4% (IQR +0.6% to +2.9%) across 247 covered loci; at the
  35 loci with 3 or more same-tract wrong-length ops, 481 stutter reads are flagged `partial_alt`
  against 35 not flagged (93%). A repeat-aware gate for deletions under 5bp was therefore not
  built; the numbers are in [Allele Classification](allele-classification.md) (the note after the
  windowed-scan safeguards).
- **Community:** callers model stutter (for example HipSTR, GATK's STR models);
  a genotyper reports it.

### Orientation artifacts (oxoG G>T, FFPE C>T)
- **In the BAM:** artifacts concentrated in one read orientation (F1R2 or F2R1),
  which strand counts do not separate.
- **gbcms now:** forward and reverse counts per allele.
- **Should do:** to decide (O6 #180, 6.7.0): optional per-allele F1R2/F2R1 counts, or
  an orientation-imbalance diagnostic only (the project prefers no new columns).
  Measure first on FFPE-rich data: orientation skew at oxoG/FFPE-prone contexts
  against true variants on the validation matrix.
- **Community:** Mutect2's read-orientation model uses F1R2 counts (to verify).

### Duplicates and UMI consensus
- **In the BAM:** PCR duplicates are marked or not; UMI consensus masks
  disagreeing bases as N.
- **gbcms now:** the duplicate flag is honoured; N matches anything; UMI families
  are supported.
- **Should do:** decided (O2 #131): report, don't warn. Each BAM's run-start line
  gives the share of its first records flagged duplicate. ACCESS consensus and
  FORTE RNA BAMs carry no duplicate flags by design (0 in 6 sampled; IMPACT
  12–25%), so a warning would fire on every such run.
- **Counting basis against sign-out** (measured, #92,
  [2026-09-22](https://github.com/msk-access/gbcms/issues/92#issuecomment-5780022052)): a 66bp ITD
  read 334 against sign-out 887. That is a counting-basis difference, not a miss: the deduplicated
  BAM holds only about 380 I-bearing reads in total (233 exact I66 plus a truncation ladder), and
  sign-out's 887/1550 exceeds that population, implying pre-dedup or pair-level counting upstream.
  Read-level ad=334 (exact reads plus the truncations the identity band then admitted; since
  RJ-27 a truncation is `partial_alt`) was faithful to the BAM.
- **Community:** Picard MarkDuplicates; fgbio consensus calling.

### Overlapping mates
- **In the BAM:** both mates of a short fragment cover the variant.
- **gbcms now:** read counts count both; fragment counts count the molecule once.
- **Should do:** decided and done.

## 4. RNA

- **Splicing and exon edges:**
    - windows end at the read's own junction;
    - a read spliced through the event counts toward depth only;
    - ASJD reports junction usage for events at splice sites.

  See [RNA splice handling](rna-splice-handling.md). Decided and done.
- **STAR clips:** a short overhang into the next exon is clipped, so its bases
  come from that exon, not the intron. A clip that reaches an exon edge or a
  junction end (annotated, or one the variant's reads splice at), or whose aligned
  bases end within five bases of one, is not allele evidence; any other clip is
  the read's own bases and is read as in DNA (RJ-18, C15 #173). Decided and done.
- **RNA editing (A-to-I):** a catalogued editing site inside a window makes an ALT
  carrier look like partial evidence. Proposed: treat the edited base as matching
  inside windows (R3 #178). Today `rna_editing_site` flags only variants that are
  themselves editing sites. An edit inside an insertion's inserted bases has no genomic
  coordinate, so neither `--rna-editing-db` nor a catalogue mask (proposed above) reaches it:
  an isolated one is one error (RJ-26); one recurring on 2 or more fragments and 5% of the
  carriers makes them another allele (partial evidence). One of 118 FORTE insertion rows
  showed one (6 reads) (measured 2026-10-27 on the FORTE RNA slices).
- **Allele-specific expression and NMD:** RNA VAF is expression-weighted and is
  not a DNA VAF. It is reported as is: an interpretation caveat.
- **Strandedness:** antisense reads are gated when enforced, and `rna_antisense_depth`
  still tallies them (R2 #114). An intronic position takes the strand of the
  transcripts spanning it (R4 #185). Loci with no gene strand (no `--gtf`, no gene,
  or genes of both strands over the locus) are not gated: every read there counts as
  sense. A library that is not stranded needs `--strandedness unstranded`; see
  [`gbcms rna`](../cli/rna.md#strandedness-filtering).
- **A GTF transcript ID used at more than one place (C34 #216):** the GTF parser
  groups exons by transcript ID alone, after stripping a numeric version suffix, so
  an ID used on more than one chromosome or strand merges into one transcript. Exons
  of different chromosomes are sorted together, which makes introns that do not exist
  and misses real ones, and the transcript's introns sit on one chromosome only. That
  feeds the intron boundaries (`RETENTION_DOMINANT`, the exon-edge BAQ rule, the RNA
  clip edges of RJ-18), per-transcript counts and junction compatibility, and ASJD.
  The strand spans are already keyed by transcript, chromosome and strand. It can
  bite on GTFs that reuse an ID (UCSC/RefSeq: the same NM_ on chrX and chrY in the
  PARs, alternate loci) and on versions that differ only by suffix (NM_x.1 and
  NM_x.2); GENCODE's `_PAR_Y` IDs are not stripped and stay distinct. Measured on the
  production GTF (Ensembl GRCh38.111, used for FORTE): 0 of 252,989 transcript IDs
  occur on more than one chromosome or strand, so no production row is affected
  today, and RJ-18's clip edges also take every junction the reads splice at, which
  covers a missing annotated edge in most cases. Proposed: key transcripts by (ID,
  chromosome, strand), report a per-transcript count under its ID with the locus where
  an ID occurs twice, and warn once when the GTF reuses an ID.

## 5. Not the BAM, but easy to misread from counts

- **Reads ending inside the event, and low depth:** a read that ends before
  holding the shorter allele's window (its event with two flank bases) counts
  toward depth only, since it cannot show either allele; one that holds it and
  ends inside the longer allele is judged by the bases it has (RJ-23). At low
  depth, informative reads are few, so compare `alt_count + ref_count` with depth.
- **Variants not prepared against the reference (direct API use, or a failed
  reference fetch) (C24 #195):** when the semiglobal alignment is judged unreliable,
  `classify_by_alignment` re-scores with a local alignment, but its tail (the REF
  branch's nearby-evidence check and the tie branch) still reads the semiglobal
  scores. REF and ALT calls use the local scores correctly; only
  `has_nearby_evidence`, so `partial_alt`, can be affected. It runs only where the
  exact-carrier rule returns nothing (no `ref_context` or `event_ref` holding the
  event and its flank), which prepared variants have. Not reproduced end to end (it
  needs a Rust test). Open for 6.7.0: use the fallback's own scores in the tail, or
  retire the fallback and flag the leftover path through `SW_FALLBACK`.
- **Purity, clonality, contamination, clonal haematopoiesis, germline vs
  somatic:** counts are evidence, not calls. A matched normal and, for cfDNA,
  [mFSD](mfsd-report.md) help interpret them.
