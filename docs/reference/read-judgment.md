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
   census's verdict on each shape, and real-data counts adjudicated per read).
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

### Open

| ID | Question | Today | Proposal | Evidence |
|:--|:--|:--|:--|:--|
| C27 | **The ALT written across several ops** (a deletion split in two). | Partial | ALT (its bases hold the ALT). Adopted in principle (operator, 2026-10-01); lands after a prototype is measured. | Census: ALT on every such shape. |
| C28 | **A read deleting the anchor.** | Phase 3's closer haplotype: ALT or REF. | Judged by its bases: neither unless they hold an allele. Adopted in principle (operator, 2026-10-01); lands after a prototype is measured. | Census: contradicts both on every such shape. |

Evidence behind RJ-7 to RJ-9 (2026-10-01): on 105 changed DNA and WES rows the
read census counts 57,199 REF reads; develop counted 60,619 and the decided rules
57,218. The largest moves are deep slippage loci (a BRCA2 cluster where an
unannotated 1bp deletion in an A run sits inside two annotated rows' windows; a
T run). RJ-9 changed no row on RC, FORTE or WES.

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
| The ALT across ops | a deletion written as two | open (C27) |
| Anchor deleted | the anchor deleted, with or without an insertion | open (C28) |

Run `python tests/read_judgment_cases.py` for the full table: every case's call,
with the read census's verdict next to each pure-indel case.
