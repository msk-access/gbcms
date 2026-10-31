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

   **Fix the class, not the row.** Start from the mechanism: what about the
   read, the alleles or the judge makes the call wrong, stated without the row
   that showed it, and the rule that follows from the principles below. The
   row is one instance, so the cases cover the other shapes of the same
   mechanism (SNV, MNV, indel and delins; repeat and unique sequence; read
   ends, strands and qualities). A cutoff comes from a mechanism (an error rate
   per telling base, a quality, a read length), never from where a measurement
   splits or from an oracle's own cut, and is checked across the values users
   set. A rule narrower than its mechanism says why, or is recorded as a known
   partial with the general rule ticketed. Judge it on data the bug was not
   found in, against a truth that does not share the rule. Before proposing,
   look for earlier work on the same mechanism: this register, `REJECTED.md`,
   and closed issues, whose comments hold the measurements behind past
   decisions.
2. **Decide.** The operator decides; the decision gets an entry in the register
   below, stating the mechanism it follows from and naming any earlier decision
   it amends.
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
- **The sequence judge resolves sequence, not haplotypes** (`--read-judge sequence`,
  #304; operator 2026-10-25). It asks which local sequence a read's own bases spell at
  the row: REF, the row's ALT, or another allele at its site (the row's *candidate
  sequences*). Contexts (other listed rows, the guard's recurring background, alleles the
  reads spell away from the site) only explain the read's other differences, so they
  never pass for the row's allele. No linkage, no phasing: a variant linked to the row
  elsewhere never credits a read.

## Decision register

### Decided

| ID | Mechanism | Rule | Source |
|:--|:--|:--|:--|
| RJ-1 | A read that ends inside a repeat tract cannot show how many units its molecule holds, so its bases fit REF and the indel alike. | REF only from reads that span an informative window; a read starting or ending inside the tract is depth only. Amended 2026-10-15 (RJ-27): for a pure indel the window is its tract with one base each side. | C10 #157; #278 |
| RJ-2 | Alignments can be wrong and other alleles share the locus, so only the read's own bases across the change show the allele its molecule carries. | ALT only where the read's own bases hold the ALT: carriers ending inside the tract are depth only; a read with another indel in the window on the strict path, or a same-length deletion spelling another allele, is another allele; reads that keep the anchor of an anchor-changing insertion are another allele. | #161, cluster 1 (#188, #191, #192, #121) |
| RJ-3 | A read starting or ending at an indel is placed by the flank base that anchors it; a flank base it does not read cannot place it. | Reads starting (or ending) on a pure indel's flank are read from the flank when they read that base (unmasked, the reference's). A deciding base past a contig edge leaves the read depth only, counted and warned. Amended 2026-10-15 (RJ-27): one flank base decides, and a read whose first base is that flank base but unread or not the reference's is placed from its other flank. | H3 #204; #278 |
| RJ-4 | A complex event changes several bases at once and a partial match fits other alleles too, so only a read showing the whole event carries it. | Complex variants (delins, anchor-changing indels, MNP reads with an indel at the block) count exact carriers: the read's bases across the whole event with two flank bases; long events by junction windows. | C1 #141 |
| RJ-5 | Listed alleles overlap: a read whose bases spell another listed allele inside this row's window carries that allele, and one fitting two listed alleles equally carries neither for certain. | A REF read that is ALT for a co-annotated sibling whose event lies inside this row's discrimination window is not REF here (partial evidence), for reads and fragments alike. Amended 2026-10-22 (#300, operator): an ALT read the exact-carrier rule judged is contested by the siblings on its bases, wherever it could lie, never re-tested against REF on its aligned slice (that slice of a carrier the aligner wrote a unit off reads as REF, so listing a sibling zeroed the row's AD). The read, every readable base of it, is placed ungapped anywhere on this row's ALT and on each sibling's ALT, over the whole-read reach, its confident mismatches counted up to three: a sibling placing it with fewer, or as few, claims it (neither row's AD where the two tie: it is either allele, whether or not any read shows the sibling; partial evidence); a sibling spelling this row's own allele is no competitor; a read no allele places within three is judged by the aligned sibling tests as before. Other ALT calls keep the aligned guard. | C2 #119; #300 |
| RJ-6 | An indel's length is part of its allele, so a read carrying another length carries another allele. | Wrong-length pure indels are another allele (partial evidence); deletions of 50bp or more match within a 3-base band (retired by RJ-24: such deletions are judged by their bases at the POS breakpoint). | #91 |
| RJ-7 | A read with its own length change where the alleles differ does not show REF there: its molecule carries something else. | A read with another insertion or deletion inside the discrimination window (not the ALT at another placement) is not REF: neither, with partial evidence. Extends RJ-5 from annotated siblings to any indel. On the exact-carrier rule (2026-10-15, #287): a read whose bases fit REF over the windows, but whose own alignment changes the length inside them with an indel that does not recur as background (on 3 or more fragments and 20% of those spanning the windows, the event's own representation aside), is judged by its bases read on without gaps to its end or next splice; two or more confident mismatches against the reference show the indel (an insertion whose first bases repeat the flank, which the windows do not reach), and the read is another allele. Bases decide: a read whose bases still fit REF (one ending inside a run) stays REF, and a germline indel beside the event keeps its reads' call. | C26 #200, operator 2026-10-01; C39a #287 |
| RJ-8 | Bases outside the discrimination window do not tell the alleles apart, so a change there says nothing about which allele the read carries. | Another insertion or deletion outside the window, of any length, is a separate event: the read is REF where its bases across the window are REF, with no partial evidence. | C26 #200, operator 2026-10-01 |
| RJ-9 | Every base a read holds is evidence: a fit at the junction is overturned by later bases that contradict it. | A long complex event's junction windows are read on inward as far as the read reaches; a read whose later bases contradict an allele is not that allele. A read ending inside the event is judged by the bases it has (one base-identical to REF counts REF). | C25 #199, operator 2026-10-01 |
| RJ-10 | Bases past the fragment end are adapter, not the molecule, and can fit an allele by chance. | A read ends at its fragment end: bases past it (read-through into adapter, insert shorter than the read) are neither bases nor reach, in every read and rule, the census included. They are hard-clipped as the read enters counting, as if trimmed (a masked base would still fill a window as a match). Only an inward-facing pair defines a fragment (TLEN positive on the forward read), and only adapter-like bases are clipped: soft-clipped, or at most two aligned past the boundary, none inserted. More are the molecule's: TLEN, a reference distance, leaves out an ITD's inserted bases and a mate's clipped 5' bases. | C17 #176, operator 2026-10-02 (adapter-like bases only, same day) |
| RJ-11 | A record without base qualities cannot pass the quality gate every base is judged by. | A record with absent base qualities (QUAL `*`, stored as 0xFF) is dropped by the read filter and warned once per BAM, as a record without bases is. | C19 #182, operator 2026-10-02 |
| RJ-12 | An unmapped record has no placement at the locus; mapping quality is the user's filter (`--min-mapq`), not evidence about the bases. | An unmapped record (flag 0x4) is not an alignment: dropped by the read filter. A mapped read whose mate is unmapped still counts, and mapped MAPQ-0 alignments stay countable (`--min-mapq 0`, pseudogene loci such as PMS2). | O7 #183, operator 2026-10-02 |
| RJ-13 | An aligner clips the bases it cannot place, often the variant itself; they are still the molecule's bases. | Soft-clipped bases inside the fragment are the read's own bases, judged by the same rules (an RNA read's clips are not evidence, RJ-18); split reads (SA) join their molecule and count once across the given breakpoints. Policy adopted now, built in 6.7.0. | C15 #173, C7 #144, C18 #177, operator 2026-10-02 |
| RJ-14 | An aligner can write one allele as several operations; the bases between the aligned flanks, not the operations, show the allele. | The ALT written across several insertion or deletion ops in the discrimination window counts ALT when the read's bases between its nearest aligned flanks are the ALT (masked bases fit, at least one base read; for an insertion, one of its own inserted bases, RJ-20), the ops' placement only a tie-break as in RJ-15: judged by its bases before the wrong-length and one-change rules read its ops. Amends RJ-2 and RJ-7. | C27 #201, operator 2026-10-02 |
| RJ-15 | Where an aligner writes a gap is a choice among equal spellings; the bases between the flanks are the molecule's. | A read whose own deletion covers the anchor is judged by its bases between its nearest aligned flanks: ALT, or REF, when they equal that allele (masked bases fit; at least one base read; for an insertion's ALT, one of its own inserted bases, RJ-20), the aligner's placement of its gap only a tie-break (a reference written as the anchor deleted and re-inserted is REF); otherwise neither, with partial evidence. Never Phase 3's closer haplotype. Deletion and insertion rows alike. Amended 2026-10-15 (RJ-27): a deletion whose every base between the flanks is masked is ALT where its confident flanks lie the ALT's length apart. | C28 #202, operator 2026-10-02; #278 |
| RJ-16 | A fit resting on low-quality bases is as likely an error as the allele, so a call needs at least one confident base's worth of evidence. | An exact-carrier ALT call needs evidence: its clearly read bases fit the ALT (as before), and its bases, each weighed by its quality (a base matches with 1 − e and mismatches with e/3), favour ALT over REF by at least what one base read at `--min-baseq` gives (about 2.5 log10 at 20). Bases matching both alleles cancel, so the evidence is the weight of the bases that mismatch REF less those that mismatch ALT, over the windows and, for a long event, the read's bases past them. A low-quality base counts for little instead of fitting either allele. | C16 #174, operator 2026-10-03 |
| RJ-17 | A spliced read skips the intron's bases, so it shows nothing of an event there; only an aligned block across the window reads it. | A splice is not reference coverage: a pure-indel read is informative (RJ-1, RJ-2) only when one of its aligned blocks between splices spans the window, REF and ALT alike. A read spliced inside the change interval shows none of it past the splice (the event may sit in the skipped intron) and is depth only, as a read ending there is; a deletion written right after the read's splice has no aligned flank (its bases equal REF spliced at an acceptor further on) and is depth only. Amends RJ-1 and RJ-2. | R5 #198, operator 2026-10-04 |
| RJ-18 | STAR clips a junction overhang it cannot splice, so a clip at an exon edge holds the next exon's bases, not the reference beside it. | An RNA read's soft clip that reaches an exon edge or a junction end (annotated, or one the variant's reads splice at), or whose aligned bases end within five bases of one, is not allele evidence: STAR soft-clips a junction overhang it cannot splice, so such a clip holds the next exon's bases. Any other clip is the read's own bases (a local aligner clips a mismatching end) and is read as in DNA. A clip admits no RNA read. Amends RJ-13 for RNA. | C15 #173, operator 2026-10-04 (refined after review, same day) |
| RJ-19 | A spliced read's molecule is the mRNA, so its bases are compared with the alleles as spliced, not with the genome across the intron. | A spliced read is judged against the haplotypes spliced at its own junctions: the exact-carrier windows (event, growth, flank) are built over the reference spliced where the read splices, the far exon's bases read from the reference, so bases a length change pushes past the junction show against the next exon, and the read must hold the spliced flank as any read holds its flank. Every junction the windows reach is followed (a short exon between two is read through). Amended 2026-10-16 (RJ-27): so is every junction their reading on reaches, so a read spliced just past the windows reads on along its next exon, not the intron. A junction starting inside the bases where the alleles differ and running past them (or ending inside them) splices the haplotypes at the event's edge, so the bases decide however the aligner wrote the gap; such a read is also REF spliced where it splices (an alternative donor or acceptor, or a splice site the event straddles), so it counts ALT only when its bases favour ALT over that reading too, by one `--min-baseq` base's evidence, and is otherwise depth only. Any other splice through those bases leaves the read depth only. | C32 #213, operator 2026-10-04 |
| RJ-20 | A masked insert shows its length but not its sequence; the length tells the ALT apart only where no other allele at the locus has it. | **Current rule** (amended 2026-10-15 with RJ-27, #278; re-confirmed 2026-10-26, #304). "Readable" is the one quality gate every base passes: not N and at or above `--min-baseq`. A read whose inserted bases are all unreadable, all N (fgbio masks a duplex disagreement to N at Q2) or all below `--min-baseq`, carries the insertion's length, not its sequence. Inside the discrimination window, an unreadable insert of the ALT's length is the ALT where no other allele of that length is at the locus (it matches both through masked bases otherwise: neither). One of another length is another allele: neither, with partial evidence. One written outside the window cannot be the variant: a separate event (RJ-8). It does not slide (RJ-21 judges a readable insert by its placements). A read whose readable inserted bases match the ALT is ALT however many are masked. In a run the aligner writes the masked base as the insertion, so the read's own inserted bases are its I-op bases, net of the reference bases its deletions beside them re-insert (counted, not placed: a masked re-inserted anchor still nets one readable base). A truncation needs one readable base, and its identity band reads its letters at any quality. The read census holds this rule. The CIGAR checks, which judge only an insertion the exact-carrier rule cannot (no prepared reference), still count an unreadable insert of the ALT's length neither, with partial evidence (`has_unreadable_insert`), never Phase 3, where REF pays for the gap so length alone would win ALT. **History:** C35 (2026-10-06) made such a read neither with partial evidence on every path; the 2026-10-15 amendment made it the ALT under the exact-carrier rule. | C35 #240, operator 2026-10-06; #278; #304 |
| RJ-21 | Inside a repeat an insertion's placement is the aligner's choice; the read's bases over every placement it can take show the allele. | A windowed insertion of the variant's length with a readable base of its own is judged by the read's bases, never Phase 3's closer haplotype (where REF pays for the gap, so length won ALT against the bases, and another allele could be absorbed into REF), and over the best placements it can take, not only the aligner's: the insert slides a junction when the base it places onto the reference fits it (the same base, or a masked one), and a slide that takes in a read mismatch resolves it, so an absorbed sequencing error or a compensating mismatch slides back; only the placements with the fewest mismatches count. It is ALT when a placement shows the ALT by read bases (at the variant's junction as the strict path judges an insert there, its one readable base possibly a flank base taken in, as RJ-20 allows; or the variant's haplotype elsewhere in a repeat, the read's only change across the window) or when the read's bases across the window spell the ALT; otherwise it is another allele (neither, with partial evidence) when a placement sits inside the discrimination window, and a separate event (RJ-8) when none does. A readable insert written past the scan window is scanned when a placement reaches the window. Whether it can be read is its own inserted bases' question (RJ-20): an insert with no readable base does not slide, so it is judged where the aligner wrote it while a readable one with a masked edge base slides (an asymmetry RJ-20's homopolymer decision requires). Both backends give the same call. Amends RJ-8 (an insert's placement is where it can sit) and RJ-20 (the unreadable outside-window case). The census trusts the aligned flank: it misreads a readable insert written just outside the window that can sit inside it, and a compensating mismatch on the flank. | C36 #243, operator 2026-10-06 |
| RJ-22 | An insert that can sit inside the window may be an allele there, and one read cannot tell its own error from a recurring other allele. | A readable insert of the variant's length is never a separate event while one of its fewest-mismatch placements sits inside the discrimination window, every junction it can reach scored by the read's readable mismatches there (through the read's contiguous aligned bases, up to 64 junctions each way), not only the junctions a slide reaches base by base: readable mismatches beside the insert (the anchor read as another base, say) block that slide, but are resolved once the insert passes the base that does not fit, so a placement inside the window can be as good as the written one. Such a read carries another allele (neither, with partial evidence), never REF; and never ALT, though that placement may show the ALT's bases, since a per-read rule cannot tell one molecule's error from a recurring other allele. The strict path keeps its junction ALT when a slide could absorb a flank substitution into another insert: the read is the given ALT plus one substitution, as a REF read with a substitution in the window stays REF. An insert with no readable base is unchanged (RJ-20), and the same haplotype written inside a run still slides to the junction as RJ-21 judges it (not seen in the measured data). Amends RJ-8 and RJ-21. The census trusts the aligned flank, so it reads such a carrier as REF; and it needs the window to be one allele exactly, so it reads a junction ALT with a substitution in the window as neither (its first limit). Amended 2026-10-15 (RJ-26): a substitution the carriers share inside the event is another allele where it recurs; one read's is an error. | C37 #245, operator 2026-10-07; #278 |
| RJ-23 | Equal-length windows at one flank give REF and ALT molecules as many reads able to judge them, so the VAF is unbiased; a change both haplotypes carry is a separate event, not evidence. | Every complex event is read by one pair of equal-length junction windows at one flank, the same for both alleles: the shorter allele's whole window (its event, grown through repeats, with two flank bases each side) against as many of the longer allele's bases, read inward from that flank, then on through the rest of the longer allele as far as the read reaches (RJ-9). Windows of equal length at one anchor give a REF and an ALT molecule the same number of read starts that can judge them, so the VAF is unbiased without padding the shorter window (which cost a third of the reads); reading junctions from either flank instead would let the longer allele count at both ends. Past 50 bases nothing changes in practice: only reads over the variant's first base (or admitted by a clip) are counted, so the far junction never decided a call. The flank is the POS flank; an RNA variant whose POS flank is cut by an exon edge reads from the other one (amended 2026-10-16, RJ-27: an exon start in the flank or at the event's first base cuts the left flank, an exon end only inside it, so an exon end where the event starts leaves the left flank whole), and when both alleles pass 50 bases the windows run from the flank nearer the first difference (an event ending a long run is out of reach from the other); a read judged from the right counts though it starts past POS, as the other allele's molecule starting there would. A change that recurs at a fixed distance beside the event is a separate event, as another indel outside the window is (RJ-8): it is masked in both windows, so a read carrying the given allele counts ALT; one on the ALT reads alone (a haplotype larger than the given allele: a phased neighbour, or an event the input under-describes) is named in `gbcms_diagnostic`, never withheld. Amended 2026-10-12: the windows run past the event through any bases the shorter allele's flank shares with the longer allele's next bases, and on through two bases from where they first differ (as a long deletion's run through its microhomology), so a read carrying only a substitution of the event's bases is neither allele; a read whose flank is clipped, under its own deletion, or aligned at the event's other copy (as mismatches over the event's far bases) is placed from the other flank as the allele its bases favour: with each allele's window placed from that flank at its own length, its confident bases (at or above the minimum quality, not N) must mismatch one allele by at least one base read at the minimum quality more than the other, so a REF molecule behind a poor clip shows no placement; it counts where that allele fits, else carries another allele (partial evidence when the ALT's), and one placed through aligned bases counts only where its allele fits, though it starts past POS. Amended 2026-10-16 (#279, operator: spanning carriers only): a read aligned straight across the windows (no indel of its own there) that fits the shorter allele read from its reading flank is confirmed from the other flank, read whole at the longer allele's length (the shorter allele extended on the reading side by the reference beyond it). This holds only where its own bases place that flank unambiguously: they continue the reference past it, one error aside, and fit it worse by two bases shifted by the length change either way. Where its confident bases there favour the longer allele by two bases, the aligner wrote that allele's other copy over the reading flank (an insert, or a deleted sequence, whose end resembles the flank, aligned onto it as mismatches, the rest clipped). A read holding the longer allele's window from there and reaching past it, as one read from the reading flank does, counts as that allele. One starting inside the event shows only the other junction, which the shorter allele's molecules have no counterpart of: depth only, never the shorter allele. Counting such reads would give the longer allele read starts the shorter one cannot have. A fit of the longer allele reads its whole event and stands. A recurring change the REF fragments carry too (a germline SNP or MNP, on 3 or more of them: one REF read's error is no background) is background, masked run or not; on the ALT fragments alone an isolated change is masked and named in the larger allele, and a run of changes (the base past the window included) is another allele, beside which no larger allele is named (the tallies cannot spell it: an MNP, or an indel shifting the bases past it; the observed allele names it from the reads' alignments). Amends RJ-4, RJ-9 and RJ-19 (no padding). Amended 2026-10-14: a change the REF molecules carry inside the repeat the event is grown through, off the given allele's own bases (a germline SNP in the run), is background too, found from the REF and ALT molecules told by their bases (RJ-25, from the variant's own reads) as beside the event (on 3 or more REF molecules and 20% of those reading it, or recurring on the ALT molecules and on 3 or more and 5% of the REF molecules); it is masked in both windows where it sits on each allele (left of the given allele at the same place, right of it shifted by the length change), taking the reference's base and the change's only (a third base there is another allele, 2026-10-14), so a REF read carrying it counts REF and a carrier carrying it counts ALT; the scan beside the event lets it fit as well, so a read carrying it is still read on past the event. A change there on the ALT molecules alone stays unmasked for an event changing the length: inside the repeat it may be the event's own sequence, another allele. An MNP has nothing to place, so the run it was grown through is beside it: an isolated change on its ALT molecules there is a phased neighbour, masked and named as beside any event (found from the molecules placed by their bases, so a carrier aligned with a clip or an indel at the block counts once its fellows show it; carriers all clipped short of the windows' reach cannot show it), and a run of them is another allele. The given allele's own bases are never masked (a germline change merged into the input's allele is named by `COEXISTING_ALLELE`). Amended again 2026-10-14 (operator): an MNP's bases it keeps inside its core sit on both alleles as the grown part's do, so background there (on 3 or more and 20% of the REF molecules, as above) is masked too, in cis or in trans with the MNP: the germline haplotype's REF reads count REF and its carriers ALT, on the MNP check's route and the exact-carrier rule's alike (VCF gives such a change its own record, its haplotype REF at the MNP's row); a change there on the ALT molecules alone is another allele. Inside a delins the bases kept depend on the alignment, so a change there stays as before (reads showing it are neither, on either side). | C39 #253, operator 2026-10-08, 2026-10-09, 2026-10-12; C43 #277, operator 2026-10-14 |
| RJ-24 | A long deletion's evidence is its junction; reading on through the deleted bases would expose REF to errors ALT never meets. | A pure deletion of 50 bp or more is judged by its bases at the POS breakpoint, however the aligner wrote it (a D op, a soft clip, a split read). Both alleles are read from the left flank over one window of equal length: 4 clean flank bases, the microhomology, then 4 clean bases of the allele's own continuation (REF: the deleted bases; ALT: the right flank). A read matching one window is that allele; matching both or neither, it is neither (partial evidence when it fits ALT better). Nothing past the window is checked, for either allele: read-on through the deleted bases would expose REF to errors ALT never meets. A read whose left flank is not aligned where the windows start (in its soft clip, under its own deletion, or written over the deleted bases' end as mismatches) is placed from the right flank, which both alleles share, as the ALT's window sits there (RJ-13; REF's window ends inside the deleted bases, which no ALT molecule holds), when its confident bases (at or above the minimum quality, not N) favour the left flank over the deleted bases REF has before the right flank by one base read at the minimum quality, and by more than nothing (the base next to the junction always differs, the deletion placed leftmost; a masked flank base weighs nothing for ALT, and the deleted base it is weighed against stays read, so a germline SNP at the anchor still shows the junction; an N flank shows nothing at any minimum quality). It is ALT where it fits, else another allele; written over the deleted bases' end, only where it fits. In DNA such a read counts toward depth when the rule decides it; a hard-clipped record is never evidence on its own. A read that fits REF but whose own alignment deletes or inserts 5 bases or more inside the deleted span carries another allele (C26): partial evidence. A co-annotated sibling claims a REF read (RJ-5) only from the anchor through REF's window's last deleted base: further into the deleted bases the windows read nothing and no ALT molecule holds them, and the flank before the anchor both windows read, so a sibling there would take REF reads alone. A change recurring at a fixed offset of the windows is a separate event, masked at that offset in both windows (as RJ-23): in the flank and microhomology, which both alleles share, one the REF reads carry too (a germline SNP or MNP, a run included, on 3 or more REF reads); an isolated one of the ALT reads' own; and past the junction an isolated one on the REF reads, where a run on reads fitting REF within a base is another allele (another breakpoint), not background; a run of changes on the ALT reads alone, the base either side of the windows included, is another allele; the isolated changes on the ALT reads alone are named in `gbcms_diagnostic` with the deletion, and beside a run none is; and masks never take every offset where the alleles differ. The 3-base length band of RJ-6 is retired. The same in RNA, after splice triage (a deletion written as a splice N stays excluded). Amended 2026-10-18 (#285): a run past the junction, on the ALT reads or on reads fitting REF, is a phased neighbour, not another breakpoint, where the bases after it rejoin that allele's continuation (the 12 after it read by 3 or more of its fragments, no change recurring on any) and no other breakpoint gives the reads' bases, the run as most show them: neither the deletion at another start or of another length (up to twice its length: prep holds as many bases past a deletion of 50 bp or more as it deletes), nor the allele with one more insertion or deletion from the junction through the run. The reads' bases are a consensus, so another breakpoint must give them exactly. A neighbour is masked in both windows as an isolated change is, and on the ALT reads alone named with the deletion through the run's last base, past the windows included. An allele's own recurring changes are first taken together (those in the windows, and on past them each within 12 bases of the last; #296, folded in): where another breakpoint gives the reads' bases exactly, the 12 after the last change read by 3 or more of its fragments, none of them is masked, isolated or run, and nothing is named, so a same-length deletion a few bases right, or another unit count of a VNTR whose units differ at one base, is another allele. A change on the carriers that another input row lists (a substitution giving that base there; prep records the input's substitution rows inside the reference it holds for the deletion) is that row's allele, not the reads' own: the input is taken as correct, so they carry the given deletion with a listed neighbour, masked and not named, and it is never evidence of another breakpoint (operator, 2026-10-19). A run at the junction whose bases are the deleted bases' own is the deletion shifted right, and a run before a repeat that a deletion of another length gives whole (another unit count of a VNTR) is that deletion: each stays another allele, as does a run whose bases keep mismatching after it. For the guard a fragment is read once, its mates merged base by base (the better-quality one's where both read a base), and taken as an allele by its bases past the windows where they place it, as RJ-27 reads on: the 12 there fit that allele within one base and the other by two or more, fit neither allele shifted by up to 8 bases nearly as well (a read of another length reads the same in a repeat), and the 24 bases it reads on along that allele fit it within one all told (another allele whose bases past the windows happen to read the same goes on to differ; a read ending sooner is not placed so, as it may end before it differs); else by its four continuation bases, which a run may make look like the other allele's. Past the windows it is read up to its next splice. Named limits: another breakpoint that also carries an SNV inside the run counts as the given deletion with a neighbour; where a deletion is one unit of a VNTR whose unit is about its length, the given deletion with a phased change and another unit count read the same over every base the guard reads, so such carriers are another allele unless the input lists the change (HG002's 90-bp deletion with its phased SNV: VAF 0.11 with the deletion alone in the input, 0.58 with the SNV row, truth 0.5); where the units are identical as far as a read reaches, a read of another unit count still counts ALT (#297); a carrier whose run lies inside the windows and that carries two phased changes 13 to 24 bases past them is not placed and counts partial; a copy farther past the deletion than its length is not looked for; in RNA a run within 12 bases of a splice stays another allele unless 3 fragments read past it unspliced. Amends RJ-2 and RJ-6; builds RJ-13 for these deletions; RJ-5 reads REF's window only. | C40 #254, operator 2026-10-09, 2026-10-13; #285, operator 2026-10-17, 2026-10-18 |
| RJ-25 | A read counts for an allele only where its bases tell it from every other allele that may be at the locus; one that cannot fits either. | A read that fits the ALT counts ALT only when, for every other allele of the ALT's length at the locus, it reads a base where the two differ, as the ALT's base and at or above `--min-baseq`; otherwise it is that allele or the ALT, neither with partial evidence. The other alleles are the co-annotated siblings that spell one (a haplotype of the ALT's length differing from it inside the event and nowhere past two bases beside it) and those read in full on 3 or more fragments and on 5% of the fragments spelling the allele they are nearest, REF or the ALT (the larger when as near: an allele recurs above the errors of the molecules it could be an error of; nearness and spelling are read past the guard's masks, and for an MNP past the bases beside its own), by the variant's own reads' bases placed from the left flank as the windows read them, not where the aligner wrote the gap (the event and two flank bases each side read; for an event changing the length, the ALT placed that way fitting the eight bases past it but for two, and no worse than REF placed the same way). An allele whose every difference from the ALT is a change from the reference the REF molecules carry too, at a base the event keeps (on 3 or more and 5% of those reading it: a separate event on both haplotypes), is not another allele, co-annotated or not; the reference's own base is no such change, so an allele lacking part of the event (the germline haplotype beside a somatic change written as one MNP) is one. Nor is one differing from the ALT only where the guard masks (background, or the ALT reads' own isolated change): that is the ALT with a separate event, which no read could tell from it; nor, for an MNP, one differing only beside the MNP's own bases: with nothing to place, that is the ALT with a phased neighbour, named as beside any event. REF and ALT molecules are told by their bases as well (read through eight bases past the event, at most two changes beside it, and for an event changing the length fewer mismatches over the event and past it than the other allele placed the same way), so a carrier aligned with its gap outside the event, or shifted into a run past it, is never a REF molecule, and where the bases past the event repeat themselves the event's own bases tell the two apart. Where the reads show no other allele nothing changes: a carrier reading only the event's length, or an MNP read with a differing base below `--min-baseq`, counts ALT. Every read the MNP check calls ALT follows the same rule (one leaving a differing base unread, and one spelling another allele where the ALT keeps the reference's base, which that check does not compare), and a decomposed form counted beside the row reads the same other alleles. Amends RJ-4, and RJ-5 for siblings of the ALT's length (judged by the bases a read holds, not only by a claim). A pure insertion keeps RJ-20 until C43 judges pure indels by their bases. Amended 2026-10-14 (operator): REF alike, a read that fits REF counts REF only when, for every other allele of REF's length at the locus, it reads a base telling the two apart (found as for the ALT, from the REF molecules' window; not another allele: one whose every difference from REF is a change from the ALT's base that the ALT molecules carry too at a base the event keeps, a separate event on both haplotypes, the ALT's own base being no such change; at a base the event changes the ALT molecules may be a third allele's); otherwise it is neither, without partial evidence, save a read that is ALT for a co-annotated sibling there, partial evidence as RJ-5 says. A co-annotated sibling is another allele only where the fragments spelling it in full (past the masks) are at least one and 5% of those spelling the allele it is nearest: a component of the event the reads never show apart from it (an SNV sign-out reports beside its MNP) is part of the event, not another allele, and at depth one error read spelling it shows nothing (the annotation lowers the 3-fragment floor, not the share). Amended 2026-10-29 (operator, #304): under the sequence judge each allele the reads spell is a site state written as the bases it changes, past the guard's masks (the allele's own base where a mask takes a carrier's), as its spelling is read above: it holds no masked base as its change, so the masks stay on every state, and one that so written is a candidate the site already holds (the ALT with a listed neighbour beside it) is that candidate; spellings apart only at a mask are one allele. Written over the window with its carriers' bases, an allele held a masked base as its change and took from X and REF every read carrying the masked change, by that base alone (a 17-base duplication whose carriers all carry a listed SNV past the copies, with recurring changes masked beside it: 293 carriers per run partial). | C43 #277, operator 2026-10-14; #304, operator 2026-10-29 |
| RJ-26 | One mismatch on one read is most likely a sequencing error and costs both alleles alike; the same change recurring across fragments is a real allele. | One substitution in a read is an error, on REF and ALT reads alike: a read whose bases fit an allele but for one confident mismatch, where the other allele's reading mismatches that read base too (past the shorter allele's window, the reference's continuation, which both alleles share), counts for the allele its other bases show; the error costs both alleles alike in its evidence (RJ-16). An allele tolerates it only where the read writes that allele's length (no length change for REF, the ALT's for the ALT) or a separate event judged on both. Inside the event, or at a base the input names (the anchor), the substitution is no error when it recurs: on 2 or more fragments and 5% of those fitting the allele it is another allele at the locus (partial evidence when nearer the ALT), found from the variant's own reads. In the flank a shared substitution is a phased neighbour (RJ-23) and the reads keep their call, so the ALT written far off with a substitution in its flank is ALT. Amends RJ-16, RJ-22 and RJ-23. | #278, operator 2026-10-15 (D3, D4) |
| RJ-27 | A pure indel in a repeat is defined by its length over the tract, so a read must hold the tract and read past it, and its own indel is placed where its bases fit, not where the aligner wrote it. | A pure insertion or deletion prepared against the reference is judged by its reads' bases (the exact-carrier rule), the CIGAR checks only where the rule cannot judge it. Its windows are its tract, where REF and ALT differ over every placement (the anchor no part of it), grown through the repeat its change lengthens or shortens (on the allele holding its bases: the ALT for an insertion, REF for a deletion; other lengths of that repeat are other alleles a read must span it to tell apart), with one base each side; a repeat on the other allele only is beside the change (operator, 2026-10-16). A read holds the shorter allele's window and reads on through the longer one's, so one ending inside the tract is depth only (RJ-1) and one ending a base past it decides; a read whose first base is its anchoring flank but unread or not the reference's is placed from its other flank (RJ-3). The read's own indel inside the windows, or the nearest within 64 bases, is placed by its bases, from the windows' first base to its end, next splice or next indel, where they fit the reference with it best (the fewest mismatches over the same reach, one at most): inside the discrimination window (its tract with a base each side) a change of the ALT's length spelling the ALT (one error aside) is the ALT, the windows anchored where its bases put them, never REF; any other is another allele (partial evidence); placed only outside, it is a separate event, applied to both windows (RJ-7, RJ-8, RJ-21, RJ-22). An insert with no readable base sits where written (RJ-20); inside the window, of the ALT's length, it is the ALT. A read whose readings fit both alleles, or miss both alike, is read on past the windows to its next indel; one placed from its reading flank that fits the ALT though its bases place no length change (it writes none, or one they place nowhere) is read on too, and stays ALT only where its further bases fit the ALT better than REF. Where masks leave no base telling the alleles apart and both windows read on to one that does, a read short of it is held as far as before and counts where its own bases fit one allele alone and exactly (a clean REF read showing the reference's base where the masks take either). A fit resting on fewer bases than tell the alleles from those beside them is confirmed past the windows, where another allele's bases, shifted by its other length, show (one base past the tract may be another allele's extra base, the flank's twin): for the ALT where the read's own bases place no change or another one, a pure indel's REF where the read writes none, and either allele's fit on a tolerated mismatch, the read's bases from its first mismatch (or the window's end) must not fit the allele shifted by up to 8 bases better, by two confident bases, than as they stand where two or more mismatch: a shifted read is another allele, while a phased change beside the event (a DNP on every carrier) mismatches at its own bases under any shift and keeps the call. A deletion whose every base between the flanks is masked is ALT where its confident flanks lie the ALT's length apart (RJ-15). A change the REF molecules carry inside a deletion's deleted bases, single or several, is background on REF (RJ-23); the given allele's own bases are never masked. A spliced read is judged on its spliced reference (RJ-19), spliced inside the windows or past them within their reading on, so it reads on along its mRNA (with BAQ masking the bases at its junction, the next exon's bases decide), and read from the left where the exon edge cutting the flank is its own junction; an exon end where the event starts leaves the left flank whole (RJ-23); an unspliced read starting at an exon start that cuts the left flank, or within the five bases of overhang an aligner writes into the intron as mismatches, reads the windows up to that edge, as a spliced read does (on its mRNA the bases before it are the previous exon's, so a run crossing the junction is read from inside), while a pre-mRNA read through the intron reads its bases as they are. Soft-clipped bases are read (RJ-13). A read whose insert is a shorter slice of the ALT's, the reference after it, shows another length: another allele, not a truncated ALT (the CIGAR path's identity band stays only where the rule cannot judge; operator, 2026-10-16). Mates that disagree resolve by quality, as for every bases-judged call (none is structural). A call whose held bases are all below `--min-baseq` takes the quality of the bases it read. Amended 2026-10-20 (#294, operator; for every event the rule judges, complex ones too): where the guard's masks hide every base telling the alleles apart, prep fetches on until the windows, every base the guard could mask masked, read on to one from either flank, so the rule judges the variant instead of the CIGAR checks; a read short of that base is held as above. A change the REF fragments show as their most common change there (background, or on more than 1% of them) that gives the reads the ALT's own base, shifted onto them where the aligner wrote no length change (the base an ALT read placed from the event's far flank shows past the shorter allele's tract; a germline SNP may give the same base), is told by the fragments placed apart from the event: those reading it whose bases, where the reference and that shifted ALT differ past it (aligned with it, no length change between), show the reference's more often than the shifted ALT's, overlapping 300 bases each side of the variant's REF span (or its read window, where wider). Where three or more of them carry it, and 3% of them (above the errors at one base: random ones, and most systematic ones), a haplotype holds it, read as the ALT over a short read: masked as background, at whatever share of the REF fragments; unless the carriers' pair orientations (read 1 forward and read 2 reverse, or the reverse), against the balance of all the reads placed apart, are too one-sided for a haplotype's (a binomial chance under 0.1%, where each orientation holds a fifth of the reads placed apart): a haplotype's molecules are read in both orientations, an oxoG or FFPE artifact sits in one (operator, 2026-10-21; a library read in one orientation, as an amplicon's, weighs none; eight carriers all one way are a haplotype's chance, twelve are not). Where the change is left unmasked though two or more of them carry it, over 1%, `gbcms_diagnostic` names it, POSSIBLE_MINORITY_HAPLOTYPE (the counts stand; one carrier is an error's). Otherwise, where 15 or more read it, a background change there is the ALT's own base on reads placed through the event: not masked, and spelling no other allele. With fewer, it is masked where it is background. An error both readings share decides nothing: a read whose bases place no length change and whose two readings miss exactly the same confident bases is read on with those bases unread, never counted for the allele that tolerates the error for the length the aligner wrote. Read on past the windows, where the alleles differ at most bases (the other shifted by its length), one confident error past them the other allele's reading does not share is forgiven where that reading misses three or more, the read writes that allele's length (as within the windows), and the error is a one-off: recurring as RJ-25 counts a tolerated substitution (on two fragments and 5% of those fitting that allele) it is another allele's base, never forgiven (operator, 2026-10-21). Within the windows only a shared error is tolerated (2026-10-15): there one base may be all that tells another allele of the same length. A read the aligner placed a unit or a copy off is read whole: one whose alignment is not clean (another event inside the discrimination window, or one no placement fits, RJ-7; or, placing no length change, two or more confident changes from the reference within 300 bases of the variant's REF span) and not spliced is placed ungapped on each allele, every readable base of it, over 300 bases past either end of where the change can be placed (at least the prepared reference; prep holds it apart from the reference the windows read), the guard's masks taking their bases as in the windows. An allele fits where its windows read exactly and at most one confident base past them mismatches; the read is that allele where it fits so and both the other allele (its event's bases open, weighed wherever the read holds as many of its bases as the shorter allele's windows, from either flank) and every other allele at the locus (RJ-25, each taking the masks where its base is one they take) miss by more. A read fitting the shorter allele's whole windows exactly and the longer allele only from inside its event is depth only (RJ-23); otherwise, and where its own bases place its indel outside the windows (RJ-8) or the ALT's change inside them, or the windows judge its change as written, the read keeps its call. A call decided past the windows (where the masks made both windows read on, or a read was read on past readings that tell nothing) rests on bases read where the aligner wrote the read, and in a repeat the other allele may hold them a copy over: it stands only where the other allele, the whole read placed ungapped anywhere on it over the same reach as the reference spells it, misses more confident bases than the called allele does (no mask is laid there: a base the guard masks beside the event may be the other allele's own shifted onto the reads written through it, as where a run meets another, and would let that allele fit them; HG002's 16-bp TC insertion at chr11:5794390, operator, 2026-10-22); where the other allele fits it with at most one miss and no more than the called one, the read is depth only (operator, 2026-10-22: HG002's 89-bp duplication on chr21, whose reads written through the copy fit REF where written and the insertion exactly a copy over, 2 of them counted REF). Every bin's fetch holds the reads the guard reads (its guard window). This replaces the decision of 2026-10-17 (#294) that a mask leaving no telling base is not applied: that counted a germline SNP's haplotype as the insertion (REJ-20261020-001). Named limits: the windows anchor on the exact repeat's edge (#290), so in an imperfect repeat (units differing at single bases) a read lying inside it, or read from a flank the repeat runs on past, is placed by the aligner, not by its bases, where the windows decide it (decided past them, it is placed anywhere as above): a carrier written a unit off fits REF as written, a germline SNP's read can read as one more unit, and a read whose bases, read whole, reach no unique sequence on either side keeps its aligned call; where 15 or more fragments placed apart read it, a haplotype carrying the change that shows on fewer than three of them, or under 3%, is taken for errors, and its reads, reading as the ALT over a short read, count ALT (with 30 to 45 of them, an SNP haplotype on 5% to 10% of the molecules in about a third of runs and on 15% in a tenth; at depth, any under 3% of the molecules, its reads adding to AD in proportion); these rows are named, POSSIBLE_MINORITY_HAPLOTYPE; a systematic error giving that base on 3% of the reads placed apart in both orientations, or in one orientation on fewer carriers than tell it from a haplotype (about ten), masks it, and the ALT's reads short of the next telling base are depth only; an ALT molecule carrying a compensating indel just past the windows, written as mismatches and with the windows' telling base unread, reads as REF with one error there and counts REF. Amends RJ-1, RJ-3, RJ-7, RJ-8, RJ-15, RJ-20, RJ-21, RJ-22 and RJ-23; builds RJ-13 for pure indels. | #278, operator 2026-10-15 (D1, D5, D6, D7, D8; groups B and C), 2026-10-16 (growth, D2); #294, operator 2026-10-20 |
| RJ-28 | A read's bases at a substitution's site show its change whatever else its molecule carries beside it, so a read carrying a longer allele that holds the change carries the change too; crediting it to one listed row alone would make a row's count depend on which rows were listed, not on the read's bases. | Under the sequence judge (`--read-judge sequence`), a read whose bases carry a substitution row's change at its site is that row's ALT, whatever else is listed: a candidate sequence holding the row's change there (an MNV listed beside its component SNV, or an allele the reads spell over a window) is the row's ALT, so a read carrying the MNV counts at both rows. At the MNV's row a read carrying only one of its components is another allele (partial evidence, RJ-25). The SNV's row names, per such MNV, how many of its ALT reads spell the MNV whole: `MNV_SHARED(chrom:pos:REF>ALT:n/m)` in `gbcms_diagnostic`, the overlap shown without changing AD or adding a column. Amended 2026-10-30 (operator, #304), its indel twin: a listed indel or delins that holds the row's whole change (the row with listed neighbours: a delins that is two listed deletions together) counts its carriers as the row's ALT, so one molecule counts at each row that holds it; the row names, per such listed row, how many of its ALT reads spell it whole: `INDEL_SHARED(chrom:pos:REF>ALT:n/m)` in `gbcms_diagnostic`, no count changed, no column added. Worked row (a patient panel row; allele shapes only): a 14-base deletion, a listed 7-base deletion four bases before it, and a listed 24>3 delins that is the two together; (RD, AD, partial) (747, 0, 129) → (747, 129, 0), the 7-base deletion counting AD 130 and the delins AD 125, so one molecule counts at all three; oracle v3 calls the 129 reads ALT, v2 a sibling or "either"; the rules path (746, 3, 126). | #304, operator 2026-10-25 (phase D: O5, with the O3 diagnostic); operator 2026-10-30 (the indel twin) |
| RJ-29 | A read can fit none of the candidate sequences because it carries structure the given rows do not describe (a long clip, a large indel), not because it comes from elsewhere, so a misfit alone does not make it foreign; only the aligner's own doubt about its placement marks a read that may belong to another locus. | Under the sequence judge, a read whose best score over every one of the row's candidate sequences (REF, the ALT, the other alleles at its site, each with every context) falls below ten confident mismatches fits none of them. It is depth only where its aligner also doubts its placement (XS at least AS; MAPQ 0 where no XS is written); a read the aligner places here alone is judged by its bases at the row, as any other. The row names how many of its depth reads fit no candidate sequence: `UNEXPLAINED_READS(n/m)` in `gbcms_diagnostic`, so a row that may count foreign reads (a paralog's, reads placed at a homozygous deletion in a repeat) is visible without a new column. | #304, operator 2026-10-25 (phase B: the fit rule, xs, with its diagnostic) |
| RJ-30 | A read's body can be a copy of a stretch inside a long listed insertion (an expanded repeat, a duplication) while its aligner placed it beside the insertion; all of an insertion's bases stand at its anchor's coordinate, so a band around the read's own placement that is not widened for the insertion cannot reach that path and stands on the placement's worse score. | Under the sequence judge, a listed insertion longer than the band's slack (16 bases) inside a read's reach (its span and its length each side) but past the window its band widens for (the slack and two bases each side) sends every state of the read to the full programme over its reach, and the linear aligner drops its band for the read. A deletion does not (a path to a far one pays the gap to it); insertions of 16 bases or fewer are left to the band. The site build uses the same scores. | #369, operator 2026-10-30 (the whole read); `docs/proposals/sequence-judge.md` §10 item 8 |

### Open

None.

Evidence behind RJ-4 (2026-09-25 to 09-28), #141. The decision is
[the #141 comment of 2026-09-25](https://github.com/msk-access/gbcms/issues/141#issuecomment-5839231028);
the rounds are in the archived plan, `.agents/archive/CYCLE_6.6.0_PLAN.md` § C1. The ticket began as a
stale-score defect in the Smith-Waterman local fallback and became the exact-carrier rule once the
reads were adjudicated one by one (local data, aggregates only):
- 129 real reads took the fallback (40 signed-out complex variants in 40 IMPACT samples, SW
  backend); 110 were linked exactly to their BAM read. 77 (70%) did not span the event (±1): 30 ALT,
  10 REF, 37 tie. Of the 33 that did, ALT was right 23/23 and REF 3/4, and ties lost 2 ALT and 2 REF
  reads that had a clear answer. 30 of the fallback's 53 ALT calls came from reads that do not span
  the event.
- The first plan, flagging from the second try's (local) scores, was right for 32 of 67 REF/tie
  reads against 36 as it stood: both flag on score closeness, not read content. A content rule on
  ALT-specific 8-mers would demote 6 of the 23 true ALT carriers. The first round's 22 of 62
  second-try ALT calls (35%) without ALT-specific sequence came largely from non-spanning reads
  and detector misses.
- Every read at the 40 variants: of the ALT-like reads (at most 1 edit of the ALT over the event
  ±1), 3,276 (96%) carry the given ALT exactly over the event ±6, 54 are 1 edit off and 82 more than
  1 edit off, and for 79 of those the same off sequence recurs (3 or more reads). Only 2 of the
  fallback's 23 aligned-spanning ALT calls carry the given ALT exactly: it fires where end-to-end
  alignment was poor, then credits near-matches. Soft clips: 11 of the 12 ALT-called fallback reads
  clipped at the event carry the whole given ALT in their clipped bases, and none of the 16 REF- or
  tie-called clipped reads carries either allele.
- Default backend, 38 variants, anchor-overlapping reads judged over the event ±2: ALT calls 4,185
  (3,767 exact, 90%; 384 partial, 9%; 27 carry neither exactly), REF calls 13,390 (11,017 exact, 82%;
  2,278 partial, 17%; 69 neither). The rule keeps 90% of ALT calls and 82% of REF calls, VAF rises
  by a median of 3.8%, and 12 variants move by more than 10% (0.84× to 2.84×).
- Acceptance (25 signed-out IMPACT complex variants, truth = equal-length padded windows, reads
  holding both): `sw` equals `pairhmm` on 25 of 25 traced variants; exact REF reads are called REF
  7,457 of 7,468 times (11 neither) and exact ALT reads ALT 2,018 of 2,020 (2 partial), with no
  REF↔ALT flip; ALT calls are 98.8% exact-ALT reads (develop 89.8%), REF calls 98.5% exact-REF
  (develop 85.8%). On the RC set 8 of 1,154 rows change, all complex or in a complex row's group.
  RNA (FORTE, about 7,150 synthetic delins probes in expressed genes, 2 samples; REF share of depth,
  develop → the rule): genomic windows reached into the intron, so spliced reads were neither,
  0.987 → 0.017 at 0-1 bp from an exon edge and 0.982 → 0.484 at 2-4 bp; with windows ending at a
  read's own junction, 0.923 and 0.906 against 0.909 mid-exon.
- Other tools (VarDictJava and GATK Mutect2 source): neither counts exactly the given allele.
  VarDict matches its own CIGAR-derived allele text and its rescue moves near-matching reads in
  (soft-clip consensus with up to 3 mismatches); its REF is a one-base count and it has no
  given-allele mode. Mutect2 credits each read's best haplotype by likelihood, even when the read
  carries extra events, so under `--alleles` a read carrying a slightly different, unassembled allele
  can be credited to the given one; it counts a read whose bases tell the alleles apart even when
  the read ends partway through the event.
- Not adopted: flagging from the second try's scores, the 8-mer content rule, and that rule with
  demotion of ALT calls lacking ALT sequence (REJ-20260925-001 in
  `.agents/learnings/REJECTED.md`).

Evidence behind RJ-5's first form, the cluster contest of 6.5.0 (2026-09-22/23), #92 and #99
([adjudication](https://github.com/msk-access/gbcms/issues/92#issuecomment-5788067318),
[acceptance](https://github.com/msk-access/gbcms/issues/92#issuecomment-5802243054); design in
`.agents/archive/CYCLE_6.5.0_PLAN.md` § T1). At a hypermutation cluster in a repeat tract of BRCA2
(five annotated deletions among four 1-3bp inserts in an AG-rich stretch of about 150bp; ACCESS
duplex and simplex fragments against sign-out), each annotated row counted the union of the carriers
compatible with it: per-row counts of 39-107 against sign-out 21-45, and a sum of about 397 against
about 150 distinct ALT molecules. Every deletion op in the cluster's molecules was canonically
left-aligned, one vote per molecule, and assigned only on exact canonical identity:
- no equivalence collisions: the five annotations canonicalize to distinct alleles, so exclusive
  assignment is well-defined;
- exclusive assignment against sign-out: 28/28, 32/32 and 21/21 exact on the three 2bp deletions;
  28 vs 22 (33bp) and 41 vs 45 (14bp), within a few molecules;
- the locus is messier than annotated: 1,299 molecules carry an unannotated 1bp deletion in the
  tract (plus D11/D17 ladders) and 373 molecules carry 2 or more deletion ops in the window.

The contest that shipped in #99 (anchor-exact evidence never contested; every other ALT read
demoted to `partial_alt` if the read's own window does not favour the row's ALT over REF, if it is an
alignment-phase call on a pure indel, or if a sibling's ALT explains it strictly better) was
accepted on: the three 2bp deletions at 28/28, 32/32 and 21/21; the 33bp and 14bp deletions from 76
and 94 down to 37 and 51 against sign-out 22 and 45 (the remainder is same-length shifted ops, kept
so the shifted-representation rescue still works; the canonical-identity table above and these
figures differ for the two longer deletions because they are two different models); the other
signed-out variants of the truth panel (252 across four samples) unchanged; five IMPACT samples with
signed-out complex (delins) clusters, where every pure-deletion row matches exactly, delins rows keep
their own carriers, and the MAF and VCF input paths agree on every row; group-free panels (FLT3-ITD,
MSI, the RNA cohort) byte-identical. The release-candidate check of 2026-09-24 (28 DNA runs, 1030 of
1060 rows byte-identical; the 30 changed rows were 29 cluster rows and one `CLIP_CANDIDATES` flag)
reproduced it: 28/28, 32/32, 21/21 exact; 33bp 37 v 22; 14bp 51 v 45. Three earlier designs failed
real data and are not adopted (REJ-20260923-006 in `.agents/learnings/REJECTED.md`). Test 3 was later
replaced for exact-carrier calls by #300 (the amendment above).

Evidence behind RJ-6 (2026-09-17/18), #91
([phase 2](https://github.com/msk-access/gbcms/issues/91#issuecomment-5721200517),
[phase 5](https://github.com/msk-access/gbcms/issues/91#issuecomment-5722764306)). Before the rule a
read proving a wrong-length pure indel got a definitive call: on the large-deletion panel (unique
context, same-anchor read deletions against annotated 100 and 60bp) every length 40-120bp counted
as full ALT, and a 2bp deletion at the edge of a 10bp homopolymer read about 51% VAF where
exact-length support was about 1.4%. The synthetic matrix used 8 REF + 6 wrong-allele reads per
scenario. Real-data targets, read-truth verified and all exact on the final build: a FLT3-class 2bp
deletion ad 259 → 7 plus 252 partial (VAF 51% → 2.7%); an AR-class 9bp GGC-tract deletion ad 24 →
15 plus 9; an RNA 2bp-deletion locus ad 32 → 15 plus 17; a SMARCB1-class 19bp insertion stays 310.
The 74-row curated indel panel: 72/74 byte-identical after the rule (the two changed rows were
delins, window effects), 69/74 on the final build (exactly 5 changed rows); the 5-sample full-BAM
large-deletion and delins panel 14/15 byte-identical (the one difference a diagnostic `partial_alt`
read, with rd/ad/dp identical). A same-length insertion whose bases confidently mismatch the
expected insert (a third allele) had been absorbed into `rd`: ad unchanged everywhere, three
insertion rows moved a total of 8 provably mismatched reads from `rd` to partial, the largest 6 reads
on one row. The first draft of the placement band summed net deleted bases anywhere in the region, so
a D(2) at the anchor + M(53) + D(49) read passed it for an expected D(50); the shipped band was
placement-aware (retired by RJ-24). Four adversarial passes found 69 findings (19, 21, 7 and 21,
plus 1 more code finding), each fixed on the branch or tracked on #92. Not adopted: the
reciprocal-overlap tier and the split-suspect escape hatch to Phase 3 (REJ-20260917-001 in
`.agents/learnings/REJECTED.md`).

Evidence behind RJ-7 to RJ-9 (2026-10-01): on 105 changed DNA and WES rows the
read census counts 57,199 REF reads; develop counted 60,619 and the decided rules
57,218. The largest moves are deep slippage loci (a BRCA2 cluster where an
unannotated 1bp deletion in an A run sits inside two annotated rows' windows; a
T run). RJ-9 changed no row on RC, FORTE or WES.

Behind RJ-7's choice (#200, [decision](https://github.com/msk-access/gbcms/issues/200#issuecomment-5950921993)):
the ticket offered (a) any insertion or deletion inside the discrimination window makes the read
not REF, and (b) keeping the 1-4bp "noise" policy and making the post-walk rule consistent with it
(no flip on an unrelated op); (a) was chosen and (b) is not adopted (REJ-20261001-001 in
`.agents/learnings/REJECTED.md`). The synthetic shapes, (REF, ALT, partial) on develop: −A carriers
on a +A row in a homopolymer (10, 0, 0); +A carriers on a −A row (10, 0, 0); a D(2) inside the 8bp
span of a unique −8 row (10, 0, 0); a D(1) right after the anchor of a unique +10bp insertion row
(10, 0, 0). Real-data cases from the cluster-1 acceptance: at one RC insertion locus (+AAG in an AAAG
repeat, two samples), six reads carrying a large deletion upstream and a D1 four bases into the
insertion's window were partial on develop, kept out of the row's REF only because a co-annotated
33bp deletion's Phase 3 wrongly called them ALT; at an RC 6bp deletion in a CAG/CAGCAA repeat, five
reads whose bases are REF across the discrimination window lost REF (REF + partial on develop,
partial on the branch), each carrying a same-length deletion of other bases 23bp away, outside the
window but inside the scan window. The ticket's remark that GATK counts such reads as neither
allele: as understood, GATK's AD credits a read to the allele whose haplotype it fits clearly
better, so RJ-7 is stricter than GATK, by decision, because it follows the reads' bases, which the
census checks. On the measured change, 112 of 144 files are byte-identical and 155 rows change (REF
−4,327, partial +4,324; ALT and depth unchanged); on every FORTE row the reads moved out of REF are no
more than those whose bases contradict both alleles.

Rarity behind RJ-9 (#199, [decision](https://github.com/msk-access/gbcms/issues/199#issuecomment-5950922254)):
769 anchor-changing one-base-REF rows among 1.13M signed-out; 133 with an ALT of 6 bases or more. The
two false calls it removed, three synthetic reads each: CA>T (Del+SNV) before a 60-A run, reads
carrying only C>T, were ALT (on develop at every run length); C>TA (Ins+SNV) before the same run was
ALT too. Reads that keep the anchor with the extra A and end inside an 80-A run are base-identical to
the REF reads, so they count REF like them; placement (the CIGAR's insertion) does not decide.

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

A fuzz of 2,400 random insertion, deletion and mismatch compositions near the anchor over eight
pure-indel rows (#202 body; rerun over 2,400 reads, five seeds, on the final head of the cluster-1
change) gave 130 ALT calls on reads whose bases are not exactly the ALT, all the same on develop: 42
carry the ALT plus mismatches, which is correct, and 88 delete the anchor base and fall to Phase 3,
whose PairHMM credits the closer of two haplotypes. Those 88 are the route RJ-15 closes. Real data
from the cluster-1 acceptance: at one RC insertion locus (+AAG in an AAAG repeat, two samples),
reads deleting 24-54bp across the anchor were partial on develop and counted REF + partial once a
co-annotated 33bp deletion's false ALT was withdrawn, 10 reads across the two samples. On a synthetic
−8 row in unique sequence, five reads that only delete the anchor and the next two bases (`D3`) hold
neither allele, yet the engine counted 5 REF. Not adopted: the strict and the ALT-by-bases
prototypes (REJ-20261002-001 in `.agents/learnings/REJECTED.md`).

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

Evidence behind the RJ-23 amendment of 2026-10-16 (#279). On HG002's long insertions,
reads were compared base by base. Some fit the ALT placed from the right flank better,
by two confident bases or more, than REF over their bases left of it: carriers of an
insert whose end resembles the left flank, aligned over it. Counting all of them ALT
moved rows in band from 357 to 411 at read level. It also overshot from 50 bp on:
50-149 bp to a median VAF of 0.611, and 150 bp or more to 0.857, because carriers that
start inside the insert give the ALT read starts REF cannot have. Counting only those
reaching past the insert's start, and making the rest depth rather than REF, kept
RJ-23's one reading flank: 385 in band at read level, median 0.462. Implemented with
the other flank anchored only where the read's bases place it unambiguously (a C run
or a CAG repeat past it let deletion carriers written without their gap through).
Against the merged build: HG002 10+ bp insertions in band 334 -> 343 (mean
|VAF - 0.5| 0.157 -> 0.153), short insertions 552 -> 555, deletions unchanged. On the
gate panel, 7 insertion rows change and no oracle-ALT read is lost; WES, normals and
ACCESS, and FORTE are unchanged.

Evidence behind RJ-29 (2026-10-25), #304 phase B, HG002 truth and aggregate patient
arms. The reads an absolute budget sends to depth are mostly structural, not noisy: 30-99%
carry a soft clip of 20 or more bases (1-12% of kept reads do), on the panel 65% a CIGAR
indel of 10 or more, only 0-6% are paralog-like (XS at least AS), and their qualities
predict about 0.06 errors per read at every misfit level. Oracle agreement, panel / gate
subset: the budget applied exactly .9262 / .9624, no budget .9323 / .9637, a budget scaled
by read length .9313 / .9636, xs .9323 / .9637 (xs equals no budget on every patient arm
within 0.01 points, and drops chr20's paralog-like REF reads). The exact budget loses 1,114
true ALT and 713 true REF calls on the panel to avoid 18 wrong ALT and 261 wrong REF. The
residual risk is foreign reads the aligner places uniquely at homozygous deletions in
repeats (HG002 chr17 3.6 kb: VAF 0.59 under xs, 0.97 under the exact budget);
UNEXPLAINED_READS flags such rows. Tools with an absolute test (GATK, Paragraph,
Graphtyper) score against haplotypes assembled from the reads, which explain each read's
structure; the given rows here do not.

Evidence behind RJ-28 (2026-10-25), #304 phase D, on the 6.6.0 gate set (aggregates).
45 component-SNV rows sat beside 42 MNV rows; in 31 of the 45 the MNV accounts for at least
90% of the SNV's ALT reads. Within sign-out: develop 37 of 45 (37 of 38 informative), the
judge counting the MNV's carriers as another allele 6 of 45; SNV median VAF 0.148 against
0.000; the multi-allelic-same-start stratum 0.839 against 0.613. On a synthetic DNP the
exclusive reading credited the same reads 20, 6 or 20 times depending only on which rows
were listed: it zeroed the SNV only when the MNV was listed alone, not when it was unlisted
and not when both components were listed. Tools that count per base give the SNV's row the
read: GetBaseCounts 1.2.5 (the sign-out and FORTE's fill-out), bam-readcount, mpileup,
LoFreq, Strelka2, DeepVariant, and GATK HaplotypeCaller with given alleles. Exclusive tools
(VarDict, FreeBayes, Mutect2, Lancet2) get there by emitting one record, which a genotyper
of given rows cannot do; gnomAD keeps component SNV counts and lists MNVs separately. A
molecule carrying the MNV enters both rows' mFSD size distributions (two KS tests on one
event in the BH family); MNV_SHARED marks those rows. FORTE lists no MNV rows, and no truth
set lists an MNV with its components.

Evidence behind the #278 decision of 2026-10-15
([comment](https://github.com/msk-access/gbcms/issues/278#issuecomment-6195173189)), a prototype on every
MNP, delins and pure-indel row of the 6.6.0 gate panel (sliced locally); these are the prototype's
numbers, before the final build below:
- Group A, what the read can see. RJ-1 (a read ending inside or just past the tract is depth only,
  for both alleles) moves 27.8k REF and 8.0k ALT reads on the panel and removes the REF excess the
  CIGAR rule shows on HG002: deletions' median VAF 0.455 → 0.500, short insertions 0.429 → 0.500,
  10+ bp insertions 0.150 → 0.43. RJ-20 (an insertion carrier whose inserted bases are unreadable
  counts ALT where no other allele of that length is at the locus, as RJ-25 decides for complex
  events) gives +560 ALT on the panel, and AD +1 in all across the matched normals.
- Group B, errors in the read. RJ-22 keeps about 1.8k ALT reads, and the matching REF reads, where
  the junction ALT has a substitution in the run.
- Group C, another event in the read. Outside the windows (RJ-8) a read whose indel's best placement
  lies outside keeps its call as a separate event; that check, ahead of the read-on of #287, touches
  0.7k REF reads and 4 tests.
- Found on the prototype and fixed on the branch: a germline change joining the run a pure indel
  slides through left REF molecules reading as the ALT; in one matched normal that gave RD 437 → 0
  and VAF 1.0. With background in the run masked around the indel's own bases, and the windows
  reading on where masks hide the alleles' difference, that normal reads RD 398 / AD 8, with
  sign-out at 8.

Evidence behind RJ-26 and RJ-27 (2026-10-15/16), #278. Other tools' read-level judges
(GATK HaplotypeCaller and Mutect2, Octopus) decide an indel read by one base past the
repeat and take a single mismatch as an error unless it recurs; the recurrence floor
(two fragments) is freebayes' and bcftools' default. Octopus, built and run on the same
synthetic shapes, agreed with the rule per read in 73 of 87 cases (its disagreements:
reads it re-aligns past the window). On the 6.6.0 gate panel, against a bases-only
oracle (each read placed ungapped on REF and on the ALT, an allele called where its
placement beats the other's by more than one confident base), the full local gate (411
runs) changes 1,062 rows of the flagged strata against develop; the rows within the
oracle's tolerance are sign-out 749, develop 838, this rule 867 (any discriminating
base; 794, 921 and 968 under a two-flank oracle, 14 rows farther from it than develop
and 58 closer). On HG002's het indels (truth VAF 0.5) against the CIGAR checks: 10+ bp
insertions in band (0.35-0.65) 110 -> 334 (median VAF 0.150 -> 0.424), short
insertions 465 -> 552, deletions 529 -> 563. The window growth (on the allele holding
the indel's bases) was chosen over growth on both alleles and over none (REJ-20261016-001 in
`.agents/learnings/REJECTED.md`): HG002 long
insertions mean |VAF - 0.5| .1498 against .1510 and .1540, deletions in band 563
against 563 and 571, oracle-correct panel reads 489,324 REF and 151,958 ALT against
488,127 / 151,533 and 492,962 / 154,045, reads of another allele counted ALT 897
against 851 and 1,166. The identity band for truncated inserts (D2) added 158 ALT reads
on HG002 long insertions with truth concordance unchanged, and on the panel 13 ALT
calls the oracle agrees with against 111 it does not: not adopted (REJ-20261016-002 in
`.agents/learnings/REJECTED.md`). An adversarial
review's third-allele shapes (21, pinned) led to the confirmation past the windows: a
count threshold lost 516 oracle-ALT panel reads (carriers of one row with a phased DNP
eight bases on); the shift test keeps them and removes 91 oracle-disagreeing calls for
36 correct ones. Normals and ACCESS: AD +3 against develop over 1,402 rows. FORTE:
11 of 25 pure-indel rows (all deletions) RD -230, AD +225, most of the AD a 12-bp
deletion STAR soft-clipped on 173 of its 175 carriers. RNA (#291): on insertions and deletions at both edges of an exon, spliced and exon
reads, BAQ on and off, without a GTF, every spliced read counts for its allele once
it reads on along its mRNA (before: spliced REF reads at indels two to four bases
before an exon's end neither, partial under BAQ; develop lost the spliced ALT reads
under BAQ); on the FORTE slices one row changes, with the GTF or without (RD +1, a
spliced REF read the intron had left undecided).

Evidence behind the RJ-27 amendment of 2026-10-20 to 2026-10-22 (#294,
[decision record](https://github.com/msk-access/gbcms/issues/294#issuecomment-6364311473)), what short
reads can tell at three HG002 rows. Every read (137) was fitted against the true haplotypes and
against the same pair with the insertion swapped: 53 fit every haplotype; 21 can only be sorted by
knowing the SNP's phase; at the #294 row no fragment reaches both the SNV and a base that tells the
insertion apart. Develop's extra ALT calls at two of the rows (8 of 13) were reads showing only a SNP
that happens to be in cis; in trans the same reading is the error the ticket describes. A whole-read
PairHMM with the nearby germline changes as haplotypes reached what these reads' bases can show: all
63 reads whose bases tell were right, with none wrong or missed. Two caveats: gbcms's own PairHMM had
its end gaps reversed (#299), and in gbcms the haplotypes would have to come from the guard. The
comment called it a candidate for later; #304 took it up for package 2. Considered and not
adopted: choosing the rule by whether the ALT is present, and a cis rescue (REJ-20261023-001 in
`.agents/learnings/REJECTED.md`).

Evidence behind RJ-25 and the RJ-23 grown-event amendment (2026-10-14), measured on a
synthetic arm of 13 shapes (each event on its own contig: REF, the given ALT and another
allele of its length; 2x150 pairs aligned by BWA-MEM; NovaSeq-binned base qualities, a
share of molecules with the telling bases or the whole core at Q2 in the stress shapes),
on the 6.6.0 gate panel (1,170 rows), its matched normals and ACCESS runs (1,402), TEMPO
WES (686) and 800 HG002 germline het delins (GIAB v5.0q), every changed read judged by
its own bases:
- Another allele of the ALT's length: develop credited 165 of its reads as the ALT
  under the stress masks, and at natural quality 15 where the given allele is absent
  (a phantom ALT, on a germline allele at 50%); for a somatic MNP at 10% on a germline
  het's haplotype, a third of its AD (VAF 0.207 against 0.151 by informative reads).
  The rule: none in any shape, VAF 0.137. With no other allele it loses what develop
  loses (202 carriers over three shapes); demoting every read that left the ALT's bases
  unread lost 286 there, 13 reads at the normals that their bases show to be the ALT,
  and moved HG002 (REJ-20261014-001).
- On real data, against develop, the final build changes 18 gate-panel rows (RD +197,
  AD -22), 9 matched-normal and ACCESS rows (RD +347), 3 WES rows (RD +28) and 3 HG002
  rows (RD +9); the FORTE RNA cohort's one delins is unchanged. By rule:
  - Another allele of the ALT's length: AD -43 at 10 panel rows, each beside another
    allele the reads show (MNP sub-alleles in a hypermutated gene, co-annotated MNP
    siblings, the germline haplotype under an MNP that merged it with a somatic change,
    and a 3-to-36 duplication whose carriers all carry two SNVs in the copy: the given
    ALT has no exact carrier and reads 0 with `PARTIAL_DOMINANT`). Of these, 20 carriers
    at 4 rows were demoted only for a co-annotated component the reads never show
    apart from its MNP, and count again; 2 more at one row (a GC>AG beside a
    co-annotated GC>AA on about 5% of its molecules: AD 123 -> 125, sign-out 126) count
    under the seen share. Of the panel's MNP ALT reads 3.3% leave a differing base
    unread, and 0.6% read nothing that tells them from an allele the reads show.
  - REF alike: 87 panel reads at 13 rows, 62 at 8 normal and ACCESS rows, 33 at 2 WES
    rows and 3 at HG002 rows leave REF. Of the panel's, 83 read the base telling REF
    from another allele below `--min-baseq` and 4 never reach it; of WES's, 31 sit at
    the two rows whose given REF holds a homozygous germline change, where no molecule
    carries the given REF (named by `OBSERVED_ALLELE`).
  - The grown event: +284 REF reads at 2 panel rows, +409 at a matched normal, +61 at a
    WES row and +12 at an HG002 row.
  - The REF side's separate event, RJ-5's partial evidence for a REF read carrying a
    co-annotated sibling, and the background at the bases an MNP keeps change no real
    row; the synthetic arm is unchanged by them.
- Every MNP and delins row of the 6.6.0 gate panel (sliced during the HPC downtime; the
  6.6.0 release reproduces all 73,649 of the gate's recorded non-VCF cells on the
  slices, 23 RNA per-transcript lists in another order): of 742 MNP and 180 delins rows
  over 235 runs, 39 MNP rows change (RD -203, AD -49) and 13 delins (RD +230, AD -7);
  none of the 655 neighbouring rows. By cause: another allele the reads show, or REF
  alike, at 24 rows (RD -112, AD -55; e.g. CC>AA beside its component CC>CA on about
  40 reads: AD 33 -> 10, each of the 23 reading the telling base below `--min-baseq`);
  the given REF on almost no molecule (a homozygous germline change merged into the
  input's allele) at 12 rows (RD to 0 or 1); rows whose reads carry another allele
  (`ZERO_ALT`, `OBSERVED_ALLELE`) at 9, losing develop's 32 partly read ALT reads; 43
  clipped carriers with a phased neighbour at 5 MNP rows, whose plain fellows the MNP
  check already counted; and the grown event at 2 delins (RD +284). Against sign-out
  VAF 14 rows move closer and 38 farther: sign-out's counter credits reads that cannot
  tell the given allele from another the reads show, counts a merged germline input's
  somatic part alone, and leaves clipped carriers out.
- Decided on the fifth review (operator, 2026-10-14), each measured on every arm:
  - The seen share: at depth one sequencing error spelled a co-annotated component in
    full and made it another allele. In a Monte Carlo (CC>GT co-annotated with C>G, no
    true component, 0.1% substitution errors, 5% of bases below min BQ), 3 of 60 loci
    flipped at 200x and 9 of 60 at 1000x, each losing about 5% of RD and AD (47 RD at
    1000x); with the share, none at either depth. A real sibling under 5% of its nearest
    allele is not looked for: the reads that fit both, its low-quality ones, count for
    the row, at most about 0.25% of depth.
  - The bases an MNP keeps: GCA>ACT with a germline T at the kept base gave RD 20 / AD
    10 in trans (the germline haplotype's REF reads neither: VAF 0.33 for 0.20), and in
    cis carriers written with an indel at the block were neither (AD 5) while plain ones
    counted. Masked as background, both give RD 40 / AD 10 on both routes. VCF gives such
    a change its own record, its haplotype REF at the MNP's row; the original GBCMS
    (strict span) calls it neither. 13 of the panel's 32 MNPs keep an interior base;
    none has such a change.
- Other tools: GATK HaplotypeCaller and Mutect2 score each read against haplotypes and
  count it in AD only when one allele beats the next by 0.2 log10, so a read that cannot
  tell the given allele from an assembled germline haplotype is left out, and counts
  where no such haplotype exists; Octopus credits reads to called haplotypes and drops
  exact ties; the original GBCMS drops a DNP read with any base below its threshold,
  other allele or not (as the first cut did). A two-allele tally (REF against ALT) cannot
  see the other allele: it credits all 69 carriers of the duplication above. The read
  census, given the true allele set, agrees with the rule.
- Inside the grown event: REF reads holding the given REF allele and a germline SNP in
  the repeat the event grows through were neither, silently: 2 panel rows (about 140 REF
  reads each, RD down 33-47%; VAF 0.43 read as 0.62 at one), one matched normal (410
  reads) and one WES row (64). HG002 has 3 such rows, none affected. Masked, they count
  REF.
- Limits: an allele on fewer than 3 fragments is not found (at low depth the call is
  develop's); in a long event few reads hold every event base from the left flank, so an
  allele the reads show is found less often there; and a read that cannot tell two
  alleles apart counts toward neither, so beside another allele as common as the ALT,
  with many partly read reads, the VAF reads low, as GATK's AD does. REF alike: a REF
  read that ends before an allele's telling base deep in a long deleted span is neither
  (a synthetic 28-to-6 delins with a het SNP in trans: RD 91 with the SNP at the 5th
  deleted base, 73 at the 25th), the read-start asymmetry equal windows remove
  elsewhere; on real data, 4 reads. Inside a delins the bases kept depend on the alignment, so a germline change
  there is not masked: reads showing it are neither on either side. The older sibling
  claim still takes reads that fit a co-annotated sibling no read shows (#282).

Evidence behind RJ-24 (2026-10-09), measured for the 6.7.0 plan and C40 on the 6.6.0 gate
panel's 60 deletions of 50+ bp (IMPACT, ABRA2-realigned) and 300 HG002 het deletions
(GIAB v5.0q, truth 0.5; WGS 35x, BWA-MEM only), every read judged by its own bases:
- 6.6.0 credits a carrier only from its CIGAR: on HG002 the carriers are clips (1,958 of
  the junction carriers), and 6.6.0 reads VAF 0.000 at 287 of 300 hets; equal windows
  at the POS breakpoint read 0.500.
- The 3-base length band credited no read on the panel (598 reads written with another
  length are neither by their bases) and one on HG002, whose bases show the junction.
- k = 4 is enough: a dynamic k (raised until the alleles' continuations differ at m
  positions, or until the junction has no other placement nearby) and k = 10 in
  low-complexity sequence did no better against truth (REJ-20261009-001, -002).
- Read-on cost accuracy: REF was checked over every deleted base it covered and ALT
  over four, so errors removed REF only. Without it, HG002 103/166 rows in 0.35-0.65
  (100 with it), the panel 40/60 within 10% of sign-out (29) (REJ-20261009-003).
- Recurring changes inside the windows touch 12 of 60 panel deletions (and 7 of 213
  HG002): in them, 38% of an allele's REF reads, or 72% of its ALT reads, carry the
  change; unmasked, those rows would lose most of one allele.
- Low-complexity deletions (40 of 300 on HG002, none on the panel) still read high
  (median 0.57): REF reads meet more errors in a repeat than ALT reads in the flank (C41
  #264).
- Surveyed tools that test a known junction anchor 5-16 bases past the microhomology
  (Manta, Delly, Paragraph, SvABA, Pindel); genome-wide placers need 20-25 for
  uniqueness, which a genotyper given the breakpoint does not.
- ABRA2 and clips on the panel (60 deletions of 50+ bp; #254,
  [decision](https://github.com/msk-access/gbcms/issues/254#issuecomment-6073704863)): ABRA2 writes
  97% of carriers as a D op; 3.4% are clips, in 28 of 60 rows, and 6.6.0 credits none of them (the
  6.7.0 plan counts 6,137 junction carriers written as a D op and 217 as a soft clip, 87 of them
  with 10+ junction bases). 6.6.0 reads 4-7% low there, because ALT reads with a short far side are
  clipped while REF reads need only a base or two. Under the rule about 5% of ALT reads become
  depth: they show fewer than 4 bases past the junction, and their REF counterparts drop alike. 3
  rows lose more: their ALT reads all carry a change 3-7 bases from the junction, a larger allele
  (the shared sub-decision with #253). The junction-carrier counts on HG002 differ between the
  sources: the decision comment gives 1,719 of 2,497 as clips (k = 10, "without realignment every
  carrier is a clip"), the first bullet of this list 1,958; the basis of the difference is not
  recorded.
- Junk risk: 23 junk clips at the breakpoint predict under 0.1 false carriers at 4^−4.
- Which breakpoint REF is read at: pooling REF from both breakpoints, as Manta, Delly and SVTyper do
  (Paragraph genotypes each breakpoint on its own), puts a het near 1/3. On capture data row-level VAF
  swings with depth along the event: REF at the left against the right breakpoint has median 1.02
  over 48 panel deletions, and spans 0.69-1.86 for the middle 80% of rows, 19 of 48 differing by
  more than 25%. On WGS (HG002) it does not. REF stays at the POS breakpoint (a far-breakpoint REF
  read gives HG002 0.43 against 0.50).
- Spec shapes added for the rule ([comment](https://github.com/msk-access/gbcms/issues/254#issuecomment-6085365571)):
  deletions of 500 bp and 1 kb (only one breakpoint is held and it is credited as junction
  evidence); a split carrier, a primary soft-clipped at the junction plus a hard-clipped
  supplementary, counted once; a hard-clipped record with no primary in the fetch, not evidence; REF
  at the POS breakpoint only. Not added, because measured negligible: R1/R2 disagreement (23 of
  about 195k panel fragments) and insertions longer than a read (none of 459 panel insertions
  reaches 100 bp).
- Correction of the ticket's premise (2026-10-15, interim 6.7.0 gate adjudication;
  [comment](https://github.com/msk-access/gbcms/issues/254#issuecomment-6191561532)): the "7 and 4
  exact-junction reads" at the two deletions over 1 kb were runs of N in clipped tails, with 0-1
  confident bases in the 30-mer. The panel's 1,023-bp rows are truncations: their reads carry
  1,046-1,118-bp deletions from the same POS (the #285 bullet below says 1,086-1,118-bp; both are
  as the sources give them), and no read holds the given junction, so 0 is right in 6.6.0 and
  develop alike. What the rule delivered stands: clipped carriers whose own bases hold the junction
  are credited (+184 reads across the gate's 50+ bp stratum). The limit the adjudication found is
  right-clipped carriers whose aligned bases end before the junction windows, which get no read call
  (#286).
- The sibling guard (2026-10-13): over the whole deleted span, co-annotated SNVs
  inside it withdrew 2,667 REF reads at 10 of the panel's 60 deletions, where no ALT
  molecule could carry them (one row's VAF 0.17 became 0.30); limited to the anchor
  through REF's window, 378 reads at 2 rows, and 41 -> 45/60 rows within 10% of
  sign-out. A growing continuation and a guard reading on past the windows gained
  nothing on any arm (REJ-20261013-001).
- Runs past the junction (2026-10-19, #285 with #296): of the 587 local deletions of
  50+ bp (panel 60, normals and ACCESS 77, WES 150, HG002 300), one carries a run past
  the junction whose bases rejoin the flank: a phased GG>AA 3-4 bp past a 61-bp
  deletion on all 48 carriers (0 of 444 other reads), which the run clause made
  partial (ALT 0; sign-out 49). Against develop the panel changes at three rows: that
  one (ALT 0 -> 48, the larger allele named), a 167-bp deletion whose REF molecules
  carry a germline SNP on the first deleted base where the ALT has the same base (RD
  117 -> 319, AD 293 -> 284, sign-out 352/291: they tie on the four continuation bases
  and are told by the bases past the windows, the guard reading on that
  REJ-20261013-001 tried with a growing continuation, here without it), and a 67-bp
  deletion beside its sibling two bases right (AD 47 -> 46); the gate panel's 4,606
  pure indels change at the same three, and normals and ACCESS, WES, HG002's 300 and
  FORTE not at all. The prepared reference now holds as many bases past the deletion
  as it deletes, so at three 1,023-bp panel rows the observed allele names the reads'
  1,086-1,118-bp deletions. Two adversarial reviews: a VNTR unit count, a same-length
  deletion shifted 2-3 bp right (31 of 150 random contexts counted ALT before the
  joint check; develop's 2-bp shift about 23%), mates read in the wrong order, and
  short reads of another allele; each is a test. A phased DNP at the junction counts
  ALT in 36 of 40 random contexts (the other four: the deletion shifted right exactly,
  or a DNP on every base the windows tell the alleles by), elsewhere 40 of 40. The
  one-base-aside tolerance first posted (2026-10-18) counted it ALT in only 21 of the 40: in the
  other 19 the deletion shifted right matched the reads at all but one base. The operator's
  correction was an exact match, at the cost that another breakpoint that also carries an SNV
  inside the run counts as the given deletion with a neighbour (REJ-20261018-001 in
  `.agents/learnings/REJECTED.md`). HG002's
  90-bp deletion, one unit of a VNTR of about 90-bp units whose carriers hold a phased
  SNV equal to the next unit, reads VAF 0.11 with the deletion alone in the input and
  0.58 with the SNV row (truth 0.5). The breakpoint search takes 0.1 s on a 1 Mb
  deletion.

Evidence behind RJ-23 (2026-10-08), measured for the 6.7.0 plan on the 6.6.0 gate
panel's 137 delins rows (IMPACT, every read judged by its own bases beside the engine's
call) and on 400 HG002 germline het delins (GIAB v5.0q, truth VAF 0.5; WGS 35x, bwa-mem,
no realignment):
- The padded windows cost reads, not bias: on long delins (span 20+) one-sided windows
  judge about 46% more ALT and 41% more REF reads at the same VAF, and the gate's
  sign-out agreement there goes 0/14 → 9/14. On HG002 (span 20+) one-sided windows read
  0.470–0.481 (6.6.0 0.464) with 25% more ALT reads.
- The alternatives are biased as read-start counts predict: whole alleles without
  padding favour the shorter allele (panel VAF ×1.22); junctions read from either
  flank favour the longer one (×0.87). On the four long delins (windows past 50
  bases) the written counts match the one-sided reading (REF 1,356 against 1,351):
  the reads that hold only the far junction are classified but never counted
  (REJ-20261008-001 in `.agents/learnings/REJECTED.md`).
- Recurrent neighbours: in 47 of the 137 rows the ALT reads share a confident change
  3–25 bases from the event that REF reads lack: the reads carry a larger allele.
  Withholding those reads (as first decided) was measured and reversed on
  2026-10-09: on HG002 it fired at 5 of 400 het delins, each larger allele ending on
  a phased truth SNV in cis, and took their VAF from 0.42 (5/5 within 0.35–0.65) to
  0.09 (0/5); on the panel it cut AD at the 49 rows with a recurrent neighbour by 54%
  (sign-out agreement 25 → 13 of 49). A germline neighbour in cis sits on every ALT read and no REF read, so
  reads cannot tell it from an under-described event; the counts stay the given
  allele's and the diagnostic names the larger one.
- Reading from the right needs the POS gate lifted for the reads it judges: a REF
  read starting inside the REF event holds the right windows as the ALT read
  starting at the same place does, but only the ALT read was admitted (by its clip),
  so long shapes read AD 94–97 against RD 29–51 from equal molecules.

The 2026-10-12 amendments, measured on every arm (the gate panel, matched normals of
every panel tumour, ACCESS duplex and simplex with their normals, TEMPO WES, HG002) and
on synthetic shapes:
- A one-base shorter allele whose flank repeats the longer allele's next bases left
  a window that ended before the alleles differ: in one matched normal a germline
  substitution at a delins's first base counted 235 ALT reads (6.6.0, padded: 0). Of
  the 2,399 reads the windows turned ALT on the panel against 6.6.0, 1,134 fit the ALT
  throughout and 385 were REF with one substitution, 329 of them at that one row.
- Placement from the other flank needs the flank shown: a REF molecule starting
  inside the event behind a low-quality or N clip fitted ALT through its masked
  flank. Carriers an aligner writes at the event's other copy (no clip) were left
  depth only while REF molecules from the same starts were judged; every surveyed
  genotyper re-scores a read's bases whatever its CIGAR (Delly, Manta, Paragraph,
  GraphTyper2, GATK, Strelka2, bcftools, ABRA2).
- The placement is judged on confident bases, each allele placed from the same
  anchor at its own length. A first gate weighed the read's bases against its own
  alignment, at every quality: synthetic REF molecules behind a Q12 clip (a NovaSeq
  quality bin), aligned plainly from inside the event, counted ALT (51 of 200), and
  carriers of a long event read from the right lost the ALT when clipped at the
  window's last base (the gate scored the bases past the window).
- A run of changes beside the event on the ALT fragments alone is another allele
  (count the given allele); an isolated recurring change stays a separate event,
  masked and named, as decided on 2026-10-09. A run the REF fragments carry too is
  background, as a germline SNP is: treated as another allele, a germline MNP in the
  flank left both alleles uncountable (synthetic: RD 6, AD 6 to 0, 0). Beside a run
  the guard names no larger allele: named from the isolated changes alone, it left
  out the run every ALT read carries.
- No surveyed tool pads windows. GATK, Strelka2 (five bases each side of one
  breakpoint) and Octopus count a read once its discriminating bases are read;
  Paragraph genotypes each breakpoint apart; GangSTR, ExpansionHunter and
  Varlociraptor model spanning bias instead of equalising windows.

Evidence behind RJ-22 (2026-10-07): every read call at insertion rows traced on
develop, each read judged by its local haplotype, on RC (ABRA2-realigned), FORTE
(STAR) and the WES loci (BWA only, the control with no realigner).
- The strict path's junction ALT with a substitution a slide absorbs: RC 28 of
  12,015 junction ALT reads (FLT3 ITDs, a BRCA2 AAG row), WES 1 of 1,647, FORTE 0.
  24 of the 29 substitutions appear in no other read and 5 also in non-carriers (an
  SNV on the haplotype); at every row the given ALT is the dominant non-REF
  haplotype and no other allele recurs. BWA had soft-clipped 19 of the 28 and
  written 9 with the insert elsewhere in the tract; ABRA2 moves a read only onto a
  contig with strictly fewer mismatches, and a singleton gets no contig of its own.
- The carrier counted REF, measured with the fixed build against develop's outputs:
  RC 19 reads at 2 rows, both at one BRCA2 insertion row, G>GA (an 8-A haplotype, the anchor
  read as A plus an inserted A, in 11 or more reads across two samples, which ABRA2
  wrote one and two junctions left); WES 2 (one fragment, a singleton anchor
  substitution); FORTE 2, one written past the scan window at the far end of a GCT
  tract. All REF to partial; no ALT or depth changed. Scoring the placements for ALT
  as well would credit the recurring allele (AD 23 → 28 and 28 → 40), so they decide
  only another allele against a separate event.

A survey: GATK prunes haplotypes seen in fewer than two reads and force-calls the
given alleles beside the discovered ones, so a read with one molecule's
substitution supports the given ALT; ABRA2 places it on the given allele's contig;
bcftools treats every same-length insertion as one allele;
GetBaseCountsMultiSample and bam-readcount count the CIGAR at the junction, and
count a carrier with the gap elsewhere as REF.

## The cases

Each group below is a set of shapes in `read_judgment_cases.py`, run at pure
indels in homopolymers, dinucleotide repeats, duplications and unique sequence,
and at anchor-changing events before an A run.

| Group | Shapes | Status |
|:--|:--|:--|
| REF reads | spanning the tract; ending inside it | decided (RJ-1) |
| Carriers | exact carrier; the ALT at another placement in the tract | decided (RJ-2) |
| Complex, whole windows | REF, substitution only, anchor kept with a length change (another allele, read past the run), exact carrier; ending inside the run, on the base after it, past it (10-A run) | decided (RJ-4, RJ-27) |
| Siblings | a co-annotated SNV inside the span; outside the window | decided (RJ-5) |
| Other indels inside the window | D1 or I1 near the anchor, D2 after it (placed by its bases, onto the anchor for a 10-base insertion); an MNP's REF reads with a 24-bp duplication of the flank written as an insertion right after the block, on a minority of reads (another allele) or as germline background (REF) | decided (RJ-7) |
| Other indels outside the window | D1, D5 or I1 past the tract; a carrier with one | decided (RJ-8) |
| Complex, long events | the same read haplotypes before a 60-A run | decided (RJ-9) |
| Complex, one-sided windows | a 28-to-6 delins: exact carriers and REF reads with a few flank bases on one side (eight on the right, where the event grows through a repeat; four on the left); a 60-to-2 delins: REF reads starting inside the event, exact carriers; carriers ending before a change that recurs on ALT reads only; carriers with a recurring change on the window's last base, and with a run of changes from it on; a REF read starting inside the event behind a low-quality clip; a carrier written at the event's other copy as plain matches; a 17-to-1 delins whose ALT flank repeats the deleted bases' next two: exact carriers, REF reads, and a REF read with the event's first base substituted | decided (RJ-23) |
| Other alleles at the locus | a 28-to-6 delins with another allele of its length on four fragments, one inserted base apart: carriers with that base below `--min-baseq` (another allele), carriers reading it (ALT), the other allele's reads (another allele); the same with the other allele on two fragments, or absent (the unread carriers count ALT), or co-annotated; an MNP joining a somatic change to a germline SNP beside it: reads showing only the germline base (another allele), the somatic one (ALT), and with no germline allele (ALT) | decided (RJ-25) |
| Other alleles of REF's length | GC>AT beside the germline GT: REF reads showing only the G (neither); a 28-to-6 delins with a germline change inside its deleted bases on a fifth of the REF molecules: REF reads with that base unread (neither); with no other allele, partly read REF reads (REF); GCA>ACT beside GTA on under a fifth of the REF molecules (neither), on 23% (REF: background at the base it keeps); a germline T at that base in cis with the MNP (GTA reads REF, ATT carriers ALT, those with an indel at the block too) or in trans (GTA reads REF); the T co-annotated (GTA reads partial evidence, RJ-5); a delins with a germline change at a base the ALT keeps, in cis: REF reads with that base unread (REF) | decided (RJ-25, 2026-10-14) |
| Co-annotated components | CC>GT co-annotated with C>G, no fragment spelling G then C in full: carriers reading only the G (ALT); one fragment spelling it beside 20 REF fragments (another allele); beside 100, one (ALT: an error's share) and five (another allele) | decided (RJ-25, 2026-10-14) |
| Background in the grown event | an insertion grown through a CA repeat with a germline SNP inside the repeat: REF reads (REF) and carriers (ALT) carrying it; the same SNP on the carriers alone (another allele); a third base where the background is masked, in the run or at a base an MNP keeps (another allele, on either side); a change at a base the event changes, carried by a third allele's molecules (no separate event) | decided (RJ-23, 2026-10-14) |
| Long deletions | 60 bp (unique), 60 bp with a 3-base microhomology, 1 kb: exact D carriers with four (or three) clean bases past the junction; right- and left-clipped carriers; a split carrier with a hard-clipped supplementary; a hard-clipped primary; REF reads holding four (or three) deleted bases past the microhomology; a REF read with an error inside the window, or 20 bases past it; a deletion of another length; one written as two with bases kept between; a read starting inside the deletion; carriers and REF reads with a recurring change in the window; a REF read starting inside the deletion behind a low-quality or N clip; carriers with a recurring change in the flank (written as a D or a clip) or on the windows' last base, and with a run of changes in the flank or from the last base on; REF reads, and carriers and REF reads, with a run of changes in the flank; left-clipped carriers with two, or four, low-quality flank bases; a carrier written over the deleted bases' end without a clip; REF reads carrying a co-annotated SNV past REF's window, and inside it; carriers and REF reads with a run from the window's last base on whose bases rejoin the allele's after it, or keep mismatching it; carriers with a run of three from the last base on; carriers and REF reads with a phased DNP inside the window; carriers of the deletion two bases right; carriers of a deletion one base longer, and of the given one with a phased DNP, before an A run; carriers with a DNP at the junction whose first base is the deleted one's; carriers of another unit count of a VNTR, alone and beside the given deletion's; overlapping mates whose first ends inside the bases a run must rejoin over; reads of another allele reading like ALT for 12 bases past the windows, then differing, and reads of it ending before they differ; carriers of the same deletion two or three bases right, alone and beside the given deletion's; carriers of another unit count of a VNTR whose units differ at one base; carriers whose change another input row lists, a row listing another base there, and none | decided (RJ-24; runs #285; #296) |
| The ALT across ops | a deletion written as two; a 1-base deletion in a run written D2 + I1 (at the anchor, inside the run); a 1-base insertion written I2 + D1 | decided (RJ-14) |
| Read inputs | read-through adapter base on an SNV (ALT, REF); absent qualities; a hard-clipped read at the previous complex classifier (with an unclipped control) | decided (RJ-10, RJ-11; C29 #207 is a bug fix) |
| Anchor deleted | the anchor deleted, with or without an insertion; an insertion row's anchor deleted | decided (RJ-15) |
| Unreadable inserts | +1 in a run with its inserted base N at the anchor or inside the run, or a letter below `--min-baseq`; a 10-base insert all N, all below `--min-baseq`, half masked; the anchor deleted and re-inserted with the insert unreadable (ALT, each) | decided (RJ-20, RJ-27) |
| Inserts of other bases | a 10-base insert of other bases one junction left and four right (REF); an 8-base one inside its duplicated tract (another allele); the ALT written one junction off with a compensating mismatch (ALT); the ALT with an error at its first inserted base, written one junction right (another allele) | decided (RJ-21) |
| Inserts behind a substitution | +1 in a run, the ALT at the junction with a substitution in the run that all four reads carry (another allele: it recurs in the event); the ALT with its anchor substituted, written one junction left, in a run and as a 10-base insert (another allele) | decided (RJ-22, RJ-26) |

Run `python tests/read_judgment_cases.py` for the full table: every case's call,
with the read census's verdict next to each pure-indel case.
