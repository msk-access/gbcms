---
name: decided-counting-rules
description: "Operator-decided counting behaviours to keep: wrong-length indels are another allele (partial_alt), complex events and 50+ bp deletions are judged by their bases, the 50 bp gate is an artifact-size prior, the RNA --min-mapq default of 1, pseudogene loci run at --min-mapq 0, and vaf = AD/(RD+AD). Read the register (docs/reference/read-judgment.md) to its end before you state any of them."
metadata:
  type: project
---

The register `docs/reference/read-judgment.md` is the owner of every rule below; this page keeps the
decisions that are easy to undo by accident.

## Wrong-length indels, complex events, large deletions

A row's route follows `window::allele_kind`, never the `variant_type` label; dispatch is
`check_allele_with_qual` in `rust/src/counting/engine.rs` (checked against the code, 2026-10-27;
the counting-engine skill maps it).
1. **Small pure indels in repeats:** a read carrying another length is another allele, `partial_alt`,
   never ALT or REF (RJ-6). Pure insertions and deletions under 50 bp are judged by their bases
   (RJ-27). No length tolerance is kept: read-level validation for #91 (2026-09-17) found every
   supporting read of nine real pure deletions (12 to 539 bp, five full clinical BAMs) an
   exact-length single `D` op, with no 1 to 3 bp wobble and no splits.
2. **Complex events** (delins, anchor-changing indels, MNP reads with an indel at the block):
   `classify_complex` counts a read only when its own bases carry the whole allele with two flank
   bases (RJ-4). `check_complex` and Phase 3 run only where the rule cannot judge the row (no
   prepared reference, or one without the event's flank); those reads set `carrier_fallback` and are
   counted and warned once per variant. Phase 3 decided 0 of about 4M read calls on every runnable
   arm (#265; the figure is in `.agents/archive/CONTINUITY_ARCHIVE.md`, not re-measured). A read whose own deletion covers
   the anchor is judged by its bases between its flanks, never Phase 3's closer haplotype (RJ-15).
3. **Large pure indels:** deletions of 50 bp or more go to `junction::classify`, which reads the POS
   breakpoint however the aligner wrote the gap (RJ-24); the 3-base placement band is gone. For
   insertions, a truncated copy is another allele under the exact-carrier rule (#278 D2); the CIGAR
   fallback's 90% band is the only place a truncation still counts ALT (`genotyper-not-caller`).

Superseded: the 2026-09-17 "delins stay Phase 3, do not touch" and the placement-aware 50+ bp band
(by RJ-4 #141, RJ-15 #202, RJ-24 #254 and #285, #278 D2).

**Why 50 bp (operator-confirmed):** `LONG_DELETION = 50` in `junction.rs` is an artifact-SIZE prior,
not an event-rarity claim. Slippage, stutter and alignment artifacts make small spurious indel ops
(overwhelmingly under 50 bp, mostly 1 to a few bp in repeat context); a read essentially never gets a
50 bp or longer `D` by artifact, so a big observed op is trusted as a real deletion in that molecule.
Below 50 bp no length tolerance is safe (the under-5 bp windowed noise gates of the CIGAR fallback
encode the same prior).

## RNA `--min-mapq` defaults to 1

RNA mode keeps reads STAR placed at two to four loci (MAPQ 3 or 1), counted once at the primary
alignment (operator, 2026-10-03). Junction reads tie between a gene and its processed pseudogene
(spliced to the parent, unspliced to the retrocopy), so they come out MAPQ 3. Do not switch to
unique-only (MAPQ 255) or to 0 without new evidence.
- Measured on FORTE (truth 33 samples): unique-only cost 12 real ALT reads in the truth set, 10 of
  them at one hotspot in an exon with a chr22 pseudogene copy, and 1.6% of junction fragments at the
  probes; MAPQ 0 added 6 ALT and 33 REF reads, consistent with the samples, but admits 5-plus-locus
  reads. ASE pipelines (GTEx, phASER) use unique-only plus WASP, which suits SNP allelic imbalance,
  not a genotyper. STAR's MAPQ depends only on the number of loci.
- Measurement: `docs/reference/read-filters.md`; harness `~/test/gbcms/harness/g3/`.

## Pseudogene loci run at `--min-mapq 0`

Some genes (PMS2, with its pseudogene PMS2CL) are genotyped with `--min-mapq 0` because their reads
multi-map and come out MAPQ 0 (operator, 2026-10-02). gbcms must keep mapped MAPQ-0 alignments
countable there (Mutect2 drops them by default: MappingQualityNotZeroReadFilter). Reads whose mate
is unmapped also count (the gbcms default; GATK the same; samtools and bcftools mpileup skip them as
anomalous pairs unless `-A`). At `--min-mapq 0` the MAPQ filter no longer screens MAPQ-0 junk: an
unmapped mate (flag 0x4) carrying a CIGAR once reached DP and ALT in a synthetic test (#183); the
filter now drops unmapped records first (`shared/filters.rs`).
**Any change to read admission or read inputs (filters, clips, fragment rules) is measured at the
default MAPQ and at `--min-mapq 0`, PMS2 rows included when the data has them.**

## `vaf` uses informative reads

`vaf` is `AD / (RD + AD)`: REF and ALT count only reads whose own bases show that allele. The
per-sample `total_count` is DP, every read overlapping the variant, including the neither and
uninformative ones (in a repeat, the reads that end inside the tract and cannot tell REF from the
indel). So `alt_count / total_count` mixes an informative numerator with an all-reads denominator and
is biased low. Merge's `simplex_duplex_total_count` is REF + ALT, a different meaning of "total".
Never advise AD/DP or the passthrough `t_var_freq` / `t_depth` (input sign-out values); read gbcms's
`vaf`. In repeats, also say how many reads were informative (RD + AD against DP), and that spanning
bias inflates a deletion's fraction. Why (2026-10-08): I advised ALT/depth for a 9 bp repeat-tract
deletion (7.9%); the informative VAF was about 48% in the tumour and 0% in the matched normal, a
high-fraction somatic event the clinical sign-out understated about sevenfold. The docs already said
`VAF = AD/(RD+AD)` (`docs/index.md`, `docs/reference/architecture.md`, `allele-classification.md`).
