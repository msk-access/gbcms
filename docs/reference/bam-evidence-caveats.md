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
| Allele soft-clipped at a read end | Counted for complex/MNP in DNA | Decided (extend) | #173, C7 #144 |
| Long event split across alignments | The primary's clip is read at the junction | Proposed | C18 #177 |
| Indel realignment (ABRA2) fixes or forces placements | Invariant by design (bases, not placement) | Decided; verify | D5 #155 |
| Aligner settings change clipping; hard clips drop bases | Clips read; hard-clipped bases lost | Decided; warn on hard clips | O2 #131, D5 #155 |
| Reads from paralogs or segmental duplications | MAPQ ≥ 20; `mq0_count` | To measure | O5 #179 |
| Reference bias: ALT reads map worse, clip, or go unmapped | Counts what maps | Proposed (diagnostic) | O5 #179 |
| Alt contigs and decoys split coverage | Contig naming handled | Proposed (warn) | O2 #131 |
| Base-quality calibration: binned, BQSR, BAQ | One masking rule; BAQ spares the variant's own indel | Decided (done); warn on bins | O2 #131 |
| Adapter read-through | Masked only for clipped-carrier admission | Proposed (every read) | C17 #176 |
| Stutter in repeats | Counted exactly; named by a diagnostic | Decided (done) | — |
| OxoG / FFPE orientation artifacts | Strand counts only | To measure (decide) | O6 #180 |
| Unmarked duplicates; UMI consensus N's | Duplicate flag honoured; N masked | Proposed (warn) | O2 #131 |
| Overlapping mates | Fragment counts | Decided (done) | — |
| RNA: splicing at exon edges | Splice-aware windows; ASJD | Decided (done) | — |
| RNA: STAR clips hold next-exon bases | Clipped-carrier admission off in RNA | To measure | #173 |
| RNA: editing (A-to-I) inside a window | Flag only | Proposed | R3 #178 |
| RNA: allele-specific expression and NMD | RNA VAF reported as is | Decided (interpretation) | — |
| RNA: strandedness | Gating exists; observability open | Decided | R2 #114 |
| Reads ending inside the event; low depth | Depth only | Decided (done) | — |
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
- **gbcms now:** complex variants and MNPs in DNA count a read whose clipped bases
  carry the exact allele, inside a well-defined fragment
  ([carriers in soft-clipped bases](counting-metrics.md)).
- **Should do:** decided. Extend to pure deletions and to clip-borne ITD
  insertions (#173, #144). In RNA, clipped bases may come from the next exon, so
  measure before enabling there.
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
    - hard clips drop bases entirely.
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
- **Should do:** to measure. Find how often low-mappability loci carry
  paralog-derived ALT reads, then decide on a diagnostic (O5 #179).
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
- **Should do:** proposed. A diagnostic when ALT reads map or clip clearly worse
  than REF reads; no filtering (O5 #179).
- **Community:** clip- and MAPQ-bias annotations, as above; graph and
  alt-aware references reduce the bias upstream. To survey.
- **See it:** `alt_dist_end_median` against `ref_dist_end_median`; per-allele
  MAPQ and clip rate once O5 lands.

### Alt contigs and decoys
- **In the BAM:** alt-aware GRCh38 alignment can move reads onto alt haplotype
  contigs, and decoys absorb reads, so coverage at the primary locus drops.
- **gbcms now:** contig naming is reconciled; nothing else.
- **Should do:** proposed. Warn at run start when the reference has alt contigs
  (O2 #131).
- **Community:** alt-aware aligners (bwa's `.alt` handling) and GATK's GRCh38
  guidance. To survey.

## 3. Bases that aren't the molecule's sequence

### Base-quality calibration
- **In the BAM:** Illumina binning (for example 2/12/23/37 or 2/12/24/40), BQSR
  recalibration, and heuristic BAQ all change which bases pass `--min-baseq`.
- **gbcms now:**
    - one masking rule: bases below min BQ, and N, match anything;
    - BAQ spares the variant's own indel (#166).
- **Should do:** decided and done. Warn when qualities look binned and
  `--min-baseq` falls between bins (O2 #131).
- **Community:** GATK BQSR; bcftools and samtools BAQ (`-B`/`-E`).
- **See it:** `n_count`; the per-read trace.

### Adapter read-through
- **In the BAM:** when the fragment is shorter than the read, bases past the mate's
  5' end are adapter.
- **gbcms now:** those bases are masked only when admitting a clipped carrier.
- **Should do:** proposed. Mask them in every read (C17 #176).
- **Community:** GATK checks for a well-defined fragment before trusting clipped
  bases. Adapter trimming upstream (fastp, cutadapt) is common but not universal.

### Stutter and polymerase slippage in repeats
- **In the BAM:** reads with one repeat unit more or fewer than the molecule.
- **gbcms now:** counted exactly. They are not the given allele, so they become
  `partial_alt`, and `COEXISTING_ALLELE` names a frequent one.
- **Should do:** decided and done (report, don't correct).
- **Community:** callers model stutter (for example HipSTR, GATK's STR models);
  a genotyper reports it.

### Orientation artifacts (oxoG G>T, FFPE C>T)
- **In the BAM:** artifacts concentrated in one read orientation (F1R2 or F2R1),
  which strand counts do not separate.
- **gbcms now:** forward and reverse counts per allele.
- **Should do:** to measure. Decide between optional F1R2/F2R1 counts and an
  orientation-imbalance diagnostic (O6 #180).
- **Community:** Mutect2's read-orientation model uses F1R2 counts (to verify).

### Duplicates and UMI consensus
- **In the BAM:** PCR duplicates are marked or not; UMI consensus masks
  disagreeing bases as N.
- **gbcms now:** the duplicate flag is honoured; N matches anything; UMI families
  are supported.
- **Should do:** proposed. Warn when duplicates appear unmarked (O2 #131).
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
  come from that exon, not the intron. Clipped-carrier admission is off in RNA.
  To measure, including a junction-aware reading (#173).
- **RNA editing (A-to-I):** a catalogued editing site inside a window makes an ALT
  carrier look like partial evidence. Proposed: treat the edited base as matching
  inside windows (R3 #178). Today `rna_editing_site` flags only variants that are
  themselves editing sites.
- **Allele-specific expression and NMD:** RNA VAF is expression-weighted and is
  not a DNA VAF. It is reported as is: an interpretation caveat.
- **Strandedness:** antisense reads are gated when enforced, and `rna_antisense_depth`
  still tallies them (R2 #114). Loci without a resolved gene strand (intronic,
  opposite-strand overlaps) are not gated yet (R4 #185).

## 5. Not the BAM, but easy to misread from counts

- **Reads ending inside the event, and low depth:** such reads count toward depth
  only, since they cannot show either allele. At low depth, informative reads are
  few, so compare `alt_count + ref_count` with depth.
- **Purity, clonality, contamination, clonal haematopoiesis, germline vs
  somatic:** counts are evidence, not calls. A matched normal and, for cfDNA,
  [mFSD](mfsd-report.md) help interpret them.
