# Read Judgment

How gbcms decides what one read says about one variant: REF, ALT, partial
evidence of another allele, or depth only. This page is the **spec**: a table of
read shapes and the call each gets, with the decision behind every call. The
algorithms that implement it are in [Allele Classification](allele-classification.md).

The table is executable. `tests/read_judgment_cases.py` builds every shape
(synthetic, four reads each) and `tests/test_read_judgment_spec.py` checks the
engine's call for each one. A **decided** case fails when the engine breaks its
rule. An **open** case pins today's call until its decision lands, so a change to
it shows in review instead of slipping in with an unrelated fix.

## How a change to read judgment is made

1. **Spec first.** A change to how reads are judged starts here: the shapes it
   affects, their call today, the proposed call, and the evidence (the read
   census's verdict on each shape, real-data counts adjudicated per read, and
   what other tools do: GATK, samtools/bcftools, fgbio, VarDict, Strelka2,
   freebayes, bam-readcount, LoFreq and the original GetBaseCounts, plus the
   domain's own tools (for RNA: GATK's RNA practice, STAR, ASE counters, RNA
   variant and editing callers) and the published literature, saying where a
   tool's handling is not documented). The census changes in step with a rule
   about read inputs (what a read contributes), so for such a rule it cannot be
   the evidence: check the bases in question themselves and the mates'
   alignments, ITD and indel rows first. Other tools' practice is a floor, not
   a ceiling: where none sets a standard, or the measurements support a better
   rule, propose that rule and say where and why it departs from practice.
2. **Decide.** The operator decides; the decision gets an entry in the register
   below, naming any earlier decision it amends.
3. **Then code.** The implementation turns the case table's open cells into
   decided ones. A test expectation that changes is part of the decision, never a
   side effect of the fix.

## Principles

- **Count the given allele; judge bases, not placement.** A read is ALT only when
  its own bases carry the ALT, and REF only when its bases show REF where the
  alleles differ. Alignments can be wrong, so the CIGAR's placement of an indel
  never decides a call on its own (AGENTS.md invariant 7).
- **Informative reads only.** A read that cannot tell the alleles apart (one that
  starts or ends inside an indel's repeat tract, or holds no window of a complex
  variant) counts toward depth only.
- **One oracle.** For pure indels, the read census (`tests/census.py`) judges each
  read by its bases alone; the engine's pure-indel counts are checked against it.

## Decision register

### Decided

| ID | Rule | Source |
|:--|:--|:--|
| RJ-1 | REF only from reads that span an informative window; a read starting or ending inside the tract is depth only. | C10 #157 |
| RJ-2 | ALT only where the read's own bases hold the ALT: carriers ending inside the tract are depth only; a read with another indel in the window on the strict path, or a same-length deletion spelling another allele, is another allele; reads that keep the anchor of an anchor-changing insertion are another allele. | #161, cluster 1 (#188, #191, #192, #121) |
| RJ-3 | Reads starting (or ending) on a pure indel's flank are read from the flank when they read that base (unmasked, the reference's). A deciding base past a contig edge leaves the read depth only, counted and warned. | H3 #204 |
| RJ-4 | Complex variants (delins, anchor-changing indels, MNP reads with an indel at the block) count exact carriers: the read's bases across the whole event with two flank bases; long events by junction windows. | C1 #141 |
| RJ-5 | A REF read that is ALT for a co-annotated sibling whose event lies inside this row's discrimination window is not REF here (partial evidence), for reads and fragments alike. | C2 #119 |
| RJ-6 | Wrong-length pure indels are another allele (partial evidence); deletions of 50bp or more match within a 3-base band. | #91 |
| RJ-7 | A read with another insertion or deletion inside the discrimination window (not the ALT at another placement) is not REF: neither, with partial evidence. Extends RJ-5 from annotated siblings to any indel. | C26 #200, operator 2026-10-01 |
| RJ-8 | Another insertion or deletion outside the window, of any length, is a separate event: the read is REF where its bases across the window are REF, with no partial evidence. | C26 #200, operator 2026-10-01 |
| RJ-9 | A long complex event's junction windows are read on inward as far as the read reaches; a read whose later bases contradict an allele is not that allele. A read ending inside the event is judged by the bases it has (one base-identical to REF counts REF). | C25 #199, operator 2026-10-01 |
| RJ-10 | A read ends at its fragment end: bases past it (read-through into adapter, insert shorter than the read) are neither bases nor reach, in every read and rule, the census included. They are hard-clipped as the read enters counting, as if trimmed (a masked base would still fill a window as a match). Only an inward-facing pair defines a fragment (TLEN positive on the forward read), and only adapter-like bases are clipped: soft-clipped, or at most two aligned past the boundary, none inserted. More are the molecule's: TLEN, a reference distance, leaves out an ITD's inserted bases and a mate's clipped 5' bases. | C17 #176, operator 2026-10-02 (adapter-like bases only, same day) |
| RJ-11 | A record with absent base qualities (QUAL `*`, stored as 0xFF) is dropped by the read filter and warned once per BAM, as a record without bases is. | C19 #182, operator 2026-10-02 |
| RJ-12 | An unmapped record (flag 0x4) is not an alignment: dropped by the read filter. A mapped read whose mate is unmapped still counts, and mapped MAPQ-0 alignments stay countable (`--min-mapq 0`, pseudogene loci such as PMS2). | O7 #183, operator 2026-10-02 |
| RJ-13 | Soft-clipped bases inside the fragment are the read's own bases, judged by the same rules (an RNA read's clips are not evidence, RJ-18); split reads (SA) join their molecule and count once across the given breakpoints. Policy adopted now, built in 6.7.0. | C15 #173, C7 #144, C18 #177, operator 2026-10-02 |
| RJ-14 | The ALT written across several insertion or deletion ops in the discrimination window counts ALT when the read's bases between its nearest aligned flanks are the ALT (masked bases fit, at least one base read; for an insertion, one of its own inserted bases, RJ-20), the ops' placement only a tie-break as in RJ-15: judged by its bases before the wrong-length and one-change rules read its ops. Amends RJ-2 and RJ-7. | C27 #201, operator 2026-10-02 |
| RJ-15 | A read whose own deletion covers the anchor is judged by its bases between its nearest aligned flanks: ALT, or REF, when they equal that allele (masked bases fit; at least one base read; for an insertion's ALT, one of its own inserted bases, RJ-20), the aligner's placement of its gap only a tie-break (a reference written as the anchor deleted and re-inserted is REF); otherwise neither, with partial evidence. Never Phase 3's closer haplotype. Deletion and insertion rows alike. | C28 #202, operator 2026-10-02 |
| RJ-16 | An exact-carrier ALT call needs evidence: its clearly read bases fit the ALT (as before), and its bases, each weighed by its quality (a base matches with 1 − e and mismatches with e/3), favour ALT over REF by at least what one base read at `--min-baseq` gives (about 2.5 log10 at 20). Bases matching both alleles cancel, so the evidence is the weight of the bases that mismatch REF less those that mismatch ALT, over the windows and, for a long event, the read's bases past them. A low-quality base counts for little instead of fitting either allele. | C16 #174, operator 2026-10-03 |
| RJ-17 | A splice is not reference coverage: a pure-indel read is informative (RJ-1, RJ-2) only when one of its aligned blocks between splices spans the window, REF and ALT alike. A read spliced inside the change interval shows none of it past the splice (the event may sit in the skipped intron) and is depth only, as a read ending there is; a deletion written right after the read's splice has no aligned flank (its bases equal REF spliced at an acceptor further on) and is depth only. Amends RJ-1 and RJ-2. | R5 #198, operator 2026-10-04 |
| RJ-18 | An RNA read's soft clip that reaches an exon edge or a junction end (annotated, or one the variant's reads splice at), or whose aligned bases end within five bases of one, is not allele evidence: STAR soft-clips a junction overhang it cannot splice, so such a clip holds the next exon's bases. Any other clip is the read's own bases (a local aligner clips a mismatching end) and is read as in DNA. A clip admits no RNA read. Amends RJ-13 for RNA. | C15 #173, operator 2026-10-04 (refined after review, same day) |
| RJ-19 | A spliced read is judged against the haplotypes spliced at its own junctions: the exact-carrier windows (event, growth, flank, padding) are built over the reference spliced where the read splices, the far exon's bases read from the reference, so bases a length change pushes past the junction show against the next exon, and the read must hold the spliced flank as any read holds its flank. Every junction the windows reach is followed (a short exon between two is read through). A junction starting inside the bases where the alleles differ and running past them (or ending inside them) splices the haplotypes at the event's edge, so the bases decide however the aligner wrote the gap; such a read is also REF spliced where it splices (an alternative donor or acceptor, or a splice site the event straddles), so it counts ALT only when its bases favour ALT over that reading too, by one `--min-baseq` base's evidence, and is otherwise depth only. Any other splice through those bases leaves the read depth only. | C32 #213, operator 2026-10-04 |
| RJ-20 | An insertion's ALT needs one of the read's own inserted bases read: not N and at or above `--min-baseq`, the one quality gate every base passes. A read whose inserted bases are all N (fgbio masks a duplex disagreement to N at Q2) or all below `--min-baseq` carries the insertion's length, not its sequence: neither, with partial evidence, on every path (at the junction, at another placement inside the discrimination window, across several ops, after a deleted anchor, a truncation), never Phase 3, where REF pays for the gap so length alone wins ALT. An unreadable insert outside the discrimination window cannot be the variant: a separate event (RJ-8), as a readable one is. In a run the aligner writes the masked base as the insertion, so the read's own inserted bases are its I-op bases, net of the reference bases its deletions beside them re-insert (counted, not placed: a masked re-inserted anchor still nets one readable base). A read whose readable inserted bases match the ALT is ALT however many are masked; a truncation needs one readable base, and its identity band reads its letters at any quality as before. The read census holds the same rule. | C35 #240, operator 2026-10-06 |
| RJ-21 | A same-length insertion of other readable bases near the variant is judged by the read's bases, never Phase 3's closer haplotype (where REF pays for the gap, so length won ALT against the bases, and another allele could be absorbed into REF), and where it can sit, not only where the aligner wrote it: the insert slides a junction when the read base it places onto the reference is read and matches it (an absorbed sequencing error or a compensating mismatch slides back; an unreadable inserted base never leaves the insert, RJ-20). At the placement nearest the variant's junction it can reach: at the junction it is judged as the strict path judges an insert there (ALT when its readable bases are the ALT's, another allele otherwise, partial when none is readable); a read whose bases across the window spell the ALT is ALT; otherwise it is another allele when it can sit inside the discrimination window (neither, with partial evidence), and a separate event when it cannot (RJ-8). Both backends give the same call. The census trusts the aligned flank, so it cannot judge a slid insert or a compensating mismatch on the flank. | C36 #243, operator 2026-10-06 |

### Open

None.

Evidence behind RJ-7 to RJ-9 (2026-10-01): on 105 changed DNA and WES rows the
read census counts 57,199 REF reads; develop counted 60,619 and the decided rules
57,218. The largest moves are deep slippage loci (a BRCA2 cluster where an
unannotated 1bp deletion in an A run sits inside two annotated rows' windows; a
T run). RJ-9 changed no row on RC, FORTE or WES.

Evidence behind RJ-10 to RJ-13 (2026-10-02). A survey of every harness input
(RC DNA, FORTE, WES) found 1,197 reads with event bases past their fragment end
at 384 rows (0.03% of reads), 61 unmapped records placed on events (27 rows), and
no records with absent qualities or hard-clipped admitted reads. Masking adapter
bases alone (measured: 244 rows, ALT −63, REF −128) left reach: a read whose
molecule ends inside a repeat still counted REF on its adapter, hence RJ-10 clips.
An adversarial review found soft-clipping with masked qualities still left reach
in complex and MNP windows (a forward molecule ending on the event's first base
counted ALT, unlike the same read trimmed), and an outward pair's TLEN read as a
fragment; the clip now removes the bases and needs an inward pair. That build
lost 13 ALT reads at an FLT3 ITD: TLEN, a reference distance, left out the
molecule's 66 inserted bases. On the 12,239 read-through reads of the changed
rows, one aligned base past the boundary mismatched the reference 75% of the time
(adapter; 9,174 reads), two 30%, ten or more 0.7% (the molecule; 97 reads), and
446 reads held inserted bases past it (three aligned bases mismatched 9%, so
they are mostly the molecule). Hence only adapter-like bases are clipped. The
final rule changes 182 rows (ALT −57, REF −81); the first base past the
boundary is A, the adapter's first base, in 97% of clipped reads, and 52 of the
54 lost SNV ALT reads showed that A.
At `--min-mapq 0` an unmapped mate carrying a CIGAR counted as an ALT read. Other
tools: GATK hard-clips adapter at the insert-size boundary and drops reads whose
bases and qualities differ in length (WellformedReadFilter) and unmapped reads
(MappedReadFilter); fgbio ClipBam clips bases past the mate (soft, soft with
mask, or hard); samtools and bcftools mpileup always discard unmapped reads, null
overlapping mates' duplicate bases by quality 0, and bcftools caps base quality at
max-BQ; bam-readcount filters nothing by default; VarDict, Strelka2, freebayes,
LoFreq and GetBaseCounts document none of these cases.

Evidence behind RJ-14 and RJ-15 (2026-10-02). RJ-14 changed no row on RC DNA,
FORTE or WES (144 of 144 files byte-identical); it decides synthetic shapes the
read census calls ALT. RJ-15 was prototyped three ways and measured against
develop and the census on the same inputs:
- strict, as C22 (an unmasked base wherever the alignment's haplotype differs
  from the ALT): ALT −26, of which 22 are reads with a masked base at the deleted
  anchor that the census and develop count ALT (the same bases placed one base
  along take the strict path and count ALT);
- the pure-indel ALT-by-bases rule, which stops at the first deciding base:
  false ALT at insertion rows (a BRCA2 row 9 → 35, census 17);
- the rule adopted: ALT net 0 on RC DNA (±1 read at 6 rows, all masked-base
  edge cases), −2 on WES; REF +57 at the BRCA2 cluster, every row toward the
  census (reads that were a co-annotated row's false ALT now count REF where
  their bases across this row's window are REF). Summed distance to the census
  on the changed indel rows: REF 344 → 251, ALT 40 → 44.

Evidence behind RJ-16 (2026-10-02/03). Synthetic delins probes in covered
sequence, where no read carries the ALT, found spurious ALT reads in data whose
reads are not consensus-collapsed: FORTE RNA 2.0 per million reads (70 at the
exon-edge and mid-exon probes), IMPACT 1.8 per million, one in WES; none in ACCESS
duplex or simplex. Each was a read whose clearly read bases fit the ALT while most
of its window was low quality (masked bases fit either allele). Three rules were
measured on the probes and on every complex DNA/WES row with ALT reads (1,185
real ALT reads): every event base read (RNA probes 70 → 16, DNA 8 → 0; 46 real
reads lost), at most one masked event base (RNA 31, DNA 0; 15 lost), and the
quality-weighted evidence adopted (thresholds 2–3 log10: RNA 15–19, IMPACT 1,
WES 0; 1–2 real reads lost). A survey: GATK HaplotypeCaller and Mutect2 weigh
every base by its quality in each read's likelihood (bases below Q18 down to Q6)
and count a read toward AD only if its best allele leads by 0.2 log10; Strelka2
(indel posterior ≥ 0.51) and bcftools (per-read indel quality) use per-read
evidence too; the counting tools use hard base-quality cutoffs; none treats a
low-quality base as fitting either allele.

An adversarial review of group 2 found three defects, fixed before merge: a read
whose bases are exactly REF, written as the anchor deleted and re-inserted, had
counted partial (RJ-15 now counts it REF, as its bases say); the evidence had
compared REF and ALT readings over different read bases, so a shared flank
base's quality could flip a call, and had skipped a long event's bases past its
junction windows (both fixed by the mismatch-weight form above). A pre-existing
limit remains: a REF molecule with one clear error just outside the window that
its ALT reading is anchored away from can still count ALT (6.7.0).

Evidence behind RJ-17 to RJ-19 (2026-10-03/04), on FORTE RNA: the truth cohort
(94 signed-out rows), the T9 exon-edge indel probes (978 rows) and the C1 splice
probes (7,224 delins rows at exon-edge distances 0–4 and mid-exon controls, where
no read carries the ALT; run on local slices of the probe regions, byte-identical
to the full BAMs). Each rule was a separate prototype against the same base.
- RJ-17 changed no truth row. At the T9 probes it moved 24,703 REF reads to depth
  on 30 rows, 24,695 of them normally spliced reads ending at the exon edge inside
  a 12–60 bp splice-crossing deletion: their bases fit both alleles (the deletion's
  first discriminating base or its margin lies past the splice). `vaf` is
  `alt_count / (ref_count + alt_count)`, so at such a deletion with carriers it
  rises to the VAF among the reads that show the event (unspliced ones), as
  `RETENTION_DOMINANT` describes; the T9 rows have no carriers and the truth set
  has no such row, so no measured `vaf` moved.
- RJ-18, first as "an RNA clip is never evidence", changed no truth row; at the
  splice probes spurious ALT 19 → 16, and 4,290 of 31.6 million REF reads
  (0.01–0.02% in every stratum) became depth only. Refined (below): a clip is
  read unless it reaches an exon edge or junction end; the three spurious ALT
  reads stay removed and 337 of those REF reads come back.
- RJ-19 changed no truth or T9 row; at the splice probes spurious ALT 19 → 4
  (exon edge 0–1 bp 14 → 2, 2–4 bp 5 → 2). REF −0.95% at edge 0–1 bp probes, from
  reads reaching only one or two bases past the junction, which hold no spliced
  flank; −0.11% at 2–4 bp, −0.01% at the controls. Re-judging the residual
  spliced ALT reads' loci against haplotypes spliced at each read's junction,
  1 of 7,461 reads fit the ALT. A first prototype that replaced the genomic
  window's intronic bases one for one with the next exon's (instead of building
  the windows over the spliced reference) required bases a repeat grew into the
  intron and doubled the REF loss.
An adversarial review of group 3 found: transcript spans keyed by ID alone (an
ID reused on another chromosome, as PAR copies or version-stripped RefSeq IDs,
gave intergenic sites a strand; spans are now per transcript, chromosome and
strand); spliced windows that spliced only the nearest junction (REF reads
through a 3 bp exon counted 0; every junction reached is now followed); a
worker silently falling back when the reference could not be opened (now an
error); a delins carrier whose aligner wrote its junction starting inside the
event counted depth only while the same bases written X D N counted ALT (now
spliced at the event's edge); and RJ-18's first form dropping mid-exon clips, which are a local
aligner's clipped mismatching ends, the read's own bases (an MNP carrier whose
second base lies in such a clip counted ALT before group 3 and not after; now
read again: splice probes 261 rows, REF +337, partial +26, no ALT). A second
review of those follow-ups found the gap form calling REF molecules ALT or
partial: a read spliced at an alternative donor or acceptor inside the event,
or at a splice site the event straddles, skips the event's bases and so fits
the shorter ALT spliced at the event's edge (a pure-REF probe showed 50% VAF
with its kept base masked). Such a read is REF as aligned at its own junction;
it now counts ALT only over that reading too, else depth only (the 22 splice-
probe rows the gap form had first changed, REF +7 and partial +19, were such
molecules and are back to depth only; the probes hold no carriers). With every
rule, against the branch point: splice probes spurious ALT 19 → 2, REF −0.39%
(−121,786, almost all at probes 0–1 bp from an exon edge), partial −453.

A survey: GATK's RNA workflow splits reads at N (SplitNCigarReads), so each piece
ends at its exon edge and HaplotypeCaller counts a piece toward AD only when its
bases favour an allele (a read partly overlapping a tandem repeat is
uninformative), and runs with `-dont-use-soft-clipped-bases`; its overhang fixer
compares a short overhang only with the intron's reference. bcftools never uses
a spliced read for indels. phASER, WASP and ASEReadCounter never read clips.
JACUSA masks bases within 6 nt of a read's own junction, SNPiR and REDItools drop
sites within 4 bp of annotated junctions; none confirms a junction-adjacent base
against the next exon, which RJ-19 does.

Other tools: GATK, Strelka2 and freebayes judge a read by its bases against
haplotypes, so a split ALT counts ALT; bam-readcount, LoFreq and the original
GetBaseCounts count CIGAR ops. A read deleting the anchor goes to GATK's
spanning-deletion allele or another candidate where one exists, is its own
haplotype in freebayes, a separate allele in VarDict and depth only in
GetBaseCounts; a closer-haplotype credit, as Phase 3 gave, appears only in the
likelihood tools without such a bucket.

## The cases

Each group below is a set of shapes in `read_judgment_cases.py`, run at pure
indels in homopolymers, dinucleotide repeats, duplications and unique sequence,
and at anchor-changing events before an A run.

| Group | Shapes | Status |
|:--|:--|:--|
| REF reads | spanning the tract; ending inside it | decided (RJ-1) |
| Carriers | exact carrier; the ALT at another placement in the tract | decided (RJ-2) |
| Complex, whole windows | REF, substitution only, anchor kept with a length change, exact carrier; ending inside the run, on the base after it, past it (10-A run) | decided (RJ-4) |
| Siblings | a co-annotated SNV inside the span; outside the window | decided (RJ-5) |
| Other indels inside the window | D1 or I1 near the anchor, D2 after it | decided (RJ-7) |
| Other indels outside the window | D1, D5 or I1 past the tract; a carrier with one | decided (RJ-8) |
| Complex, long events | the same read haplotypes before a 60-A run | decided (RJ-9) |
| The ALT across ops | a deletion written as two; a 1-base deletion in a run written D2 + I1 (at the anchor, inside the run); a 1-base insertion written I2 + D1 | decided (RJ-14) |
| Read inputs | read-through adapter base on an SNV (ALT, REF); absent qualities; a hard-clipped read at the previous complex classifier (with an unclipped control) | decided (RJ-10, RJ-11; C29 #207 is a bug fix) |
| Anchor deleted | the anchor deleted, with or without an insertion; an insertion row's anchor deleted | decided (RJ-15) |
| Unreadable inserts | +1 in a run with its inserted base N at the anchor or inside the run, or a letter below `--min-baseq`; a 10-base insert all N, all below `--min-baseq`, half masked (ALT); the anchor deleted and re-inserted with the insert unreadable | decided (RJ-20) |
| Inserts of other bases | a 10-base insert of other bases one junction left and four right (REF); an 8-base one inside its duplicated tract (another allele); the ALT written one junction off with a compensating mismatch (ALT); the ALT with an error at its first inserted base, written one junction right (another allele) | decided (RJ-21) |

Run `python tests/read_judgment_cases.py` for the full table: every case's call,
with the read census's verdict next to each pure-indel case.
