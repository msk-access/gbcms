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
   freebayes, bam-readcount, LoFreq and the original GetBaseCounts, saying where
   a tool's handling is not documented). The census changes in step with a rule
   about read inputs (what a read contributes), so for such a rule it cannot be
   the evidence: check the bases in question themselves and the mates'
   alignments, ITD and indel rows first.
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
| RJ-13 | Soft-clipped bases inside the fragment are the read's own bases, judged by the same rules; RNA exon-edge clips are excluded until measured; split reads (SA) join their molecule and count once across the given breakpoints. Policy adopted now, built in 6.7.0. | C15 #173, C7 #144, C18 #177, operator 2026-10-02 |
| RJ-14 | The ALT written across several insertion or deletion ops in the discrimination window counts ALT when the read's bases between its nearest aligned flanks are the ALT (masked bases fit, at least one base read), the ops' placement only a tie-break as in RJ-15: judged by its bases before the wrong-length and one-change rules read its ops. Amends RJ-2 and RJ-7. | C27 #201, operator 2026-10-02 |
| RJ-15 | A read whose own deletion covers the anchor is judged by its bases between its nearest aligned flanks: ALT, or REF, when they equal that allele (masked bases fit; at least one base read), the aligner's placement of its gap only a tie-break (a reference written as the anchor deleted and re-inserted is REF); otherwise neither, with partial evidence. Never Phase 3's closer haplotype. Deletion and insertion rows alike. | C28 #202, operator 2026-10-02 |
| RJ-16 | An exact-carrier ALT call needs evidence: its clearly read bases fit the ALT (as before), and its bases, each weighed by its quality (a base matches with 1 − e and mismatches with e/3), favour ALT over REF by at least what one base read at `--min-baseq` gives (about 2.5 log10 at 20). Bases matching both alleles cancel, so the evidence is the weight of the bases that mismatch REF less those that mismatch ALT, over the windows and, for a long event, the read's bases past them. A low-quality base counts for little instead of fitting either allele. | C16 #174, operator 2026-10-03 |

### Open

None. The C16 junction-placement guard (a spliced read whose aligner placed its
junction a few bases late can show the next exon's bases over an exon-edge
event) is deferred to 6.7.0 with C15 and RJ-13's RNA exon-edge clips: telling
it from a genuine carrier needs the reference at the splice's far end.

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

Run `python tests/read_judgment_cases.py` for the full table: every case's call,
with the read census's verdict next to each pure-indel case.
