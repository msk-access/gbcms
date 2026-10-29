---
name: genotyper-not-caller
description: "Operator principle: gbcms is a genotyper, not a caller. It takes the input allele as correct and counts the reads whose own bases carry it; the alignment can be wrong, so judge bases, not placement; never infer what a wrong input should have been; explain mismatches in the status-reason and diagnostic columns; rescues are opt-in. Includes where the identity band still applies."
metadata:
  type: feedback
---

gbcms is a **genotyper, not a caller** (AGENTS.md invariant 7). Given an allele, it reports how many
reads carry it, how many carry REF, and what else is there as diagnostics. It does not decide
whether a variant exists, discover alleles for a row, or rewrite a wrong input.

## Count the given allele

gbcms assumes the input allele is correct and reports accurate counts for exactly that allele. When
the reads suggest the input is wrong (a mis-described, shifted or incomplete allele), it does not
work out what the allele should have been and count that. It counts the given allele honestly (often
low ALT, high partial) and explains the situation in the existing columns: `gbcms_status` and
`gbcms_status_reason` for the verdict, `gbcms_diagnostic` for what the reads show (for example
`PARTIAL_DOMINANT`, `OBSERVED_ALLELE`). The only exceptions are opt-in rescues the user asks for,
flagged and audited per row: `--rescue-mnp` and `--rescue-homopolymer` (both off by default,
checked in `src/gbcms/cli.py` 2026-10-28). A recurrent unannotated haplotype is a finding for
validation and diagnostics, never ALT evidence for the row.

**The boundary** (sharpened with the operator, 2026-09-25):
1. **Normalize** anything the reference fully determines, and flag it: REF bases
   (`WARN_REF_CORRECTED`), anchoring, left-alignment, and equivalent representations that give the
   same ALT haplotype.
2. **Tolerate read-side noise**: sequencing errors, BQ masking, how the aligner writes the reads.
3. **Never infer a different ALT haplotype** for the counts.

**Risk:** an honest low VAF for a mis-described allele misleads anyone who ignores
`gbcms_diagnostic`, so the caveat must be specific: say what the reads carry, not only "partial".

**Test for any counting change:** does it credit a read to the row's ALT when the read does not
carry the given ALT haplotype? If yes, it is out, unless it is an opt-in rescue. Tolerating how
reads are written (aligner representation, sequencing errors, truncated inserts in reads that end
inside them) is reading the given allele accurately; tolerating a wrong input is not.

## The BAM is not the truth either

Alignments can be wrong: an indel placed elsewhere in a repeat, an allele soft-clipped at a read end,
a long event split across alignments, a misaligned read. The evidence is each read's **own bases**,
judged independently of how the aligner represented them. That is why the rules anchor windows on
flank bases, grow events through repeats, read clipped bases, and treat junctions and splices
explicitly, rather than trusting the CIGAR. Sign-out alleles are one more result to compare against,
and they can be wrong (manual edits, shifted or garbage alleles; people use IGV).

Operator quotes: 2026-09-23 (reads, not sign-out, are what we score against); 2026-09-25 ("go for the
right and accurate results ... explain with all caveats using the failure reason and diagnostic
columns"); 2026-09-28 ("BAM is not necessarily the truth, it can have alignment errors, but the point
is we are not a caller but a genotyper"; "always try to do things as generalized as possible; use
different datasets to make the architecture general").

**How to validate:** per read, against an independent, position-aware census of the reads' own bases
(clipped bases included, equal-length windows, masked low-quality bases), across data types before
merge: synthetic contracts and fuzz with varied base-quality encodings; realigned panel DNA (tumour,
normal, duplex, simplex); DNA without indel realignment (WES, WGS); RNA (spliced, BAQ); public
reference data where available. Cross-platform agreement is secondary. When the right behaviour is
not known, measure it on that matrix and survey community practice before deciding
(`.agents/rules/decision-protocol.md`). The caveats are in `docs/reference/bam-evidence-caveats.md`.

## The identity band (CIGAR fallback only)

The 90% insert-identity band, once the primary tolerance for an ALT string that was itself imperfect,
now applies only in `insert_truncation_match` (`rust/src/counting/variant_checks.rs`, called from
`resolve_anchor_insertion_candidate` inside `check_insertion`). `check_allele_with_qual` reaches
`check_insertion` for a pure insertion only when the exact-carrier rule (`by_bases`) returns nothing:
no prepared reference, or one without the event's flank. Such reads set `carrier_fallback` and are
counted and warned. Gates: observed insert at least 4 bp and shorter than expected, at least 90%
identity to the best window, both sequences non-low-complexity, at least one readable inserted base.
Nothing else uses an identity band (checked in the code, 2026-10-27). On every prepared pure insertion
a truncated copy is another allele (`partial_alt`), and one sequencing error inside an exact-length
insert is tolerated by RJ-26 (a substitution is an error unless it recurs on 2 or more fragments and
5% of those fitting the allele), not by a band.

Why the band was dropped as primary (#278 D2, operator, 2026-10-16): for truncated inserts under the
exact-carrier rule it added 158 ALT reads on HG002 long insertions with truth concordance unchanged,
and on the panel 13 ALT calls the bases oracle agrees with against 111 it does not (RJ-27 evidence in
`docs/reference/read-judgment.md`; #143). It sits on the boundary above: right only when the
mismatches are read errors. Measure it: scattered mismatches are errors; one recurring alternative
insert is a different allele or a wrong input, and must not count. Earlier measurements (2026-09-22/23,
#92): IMPACT raw reads at the three covered insertion loci of 4 bp or more (about 900 carriers) had
zero-mismatch truncations in about 24% of carriers and mismatching candidates in about 2% (mostly
Q20 to Q35 singletons at unique positions; all-low-BQ mismatches in 6 reads); ACCESS consensus reads
at the two band-eligible loci had exact truncations of 2 to 14 per locus, N-carrying reads (0 to 10
duplex, up to 83 simplex) and no same-length mismatches, so a BQ-first policy has nothing to act on.

**Standing preference:** no new CLI flags or output columns for a change here; new signal goes to
logs, existing diagnostics or validation tooling. When you touch insert matching, test raw identity
against BQ-masked policies on the divergence classes. Panels are regression validation for one assay,
never universal calibration.
