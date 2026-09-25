---
name: count-the-given-allele
description: "Operator principle — gbcms takes the input allele as correct and counts exactly it; never infer or deconvolute what a wrong input should have been; explain mismatches through gbcms_status_reason / gbcms_diagnostic. Opt-in rescues (--rescue-mnp) are the only exception."
metadata:
  type: feedback
---

gbcms assumes the input allele is correct and reports accurate counts for exactly
that allele. When the reads suggest the input is wrong (a mis-described,
shifted or incomplete allele), gbcms does not work out what the allele should
have been and count that instead. It counts the given allele honestly (often low
ALT, high partial) and explains the situation with caveats in the existing
columns: `gbcms_status` / `gbcms_status_reason` for the verdict and
`gbcms_diagnostic` for what the reads show (e.g. `PARTIAL_DOMINANT`).

The only exceptions are opt-in rescues the user asks for, like `--rescue-mnp`,
which are flagged and audited per row.

A recurrent unannotated haplotype is a finding for validation work (see
[[bam-is-truth]]), never ALT evidence for the row.

**Why:** stated by the operator 2026-09-25 ("go for the right and accurate
results … explain with all caveats using the failure reason and diagnostic
columns"). The prompt was a C1 option I proposed that kept ALT calls when reads
shared a recurring unannotated haplotype, which is deconvolution.

**The boundary** (sharpened with the operator, who invited critique of the
principle, 2026-09-25):
1. **Normalize** anything the reference fully determines, and flag it: REF
   bases (`WARN_REF_CORRECTED`), anchoring, left-alignment, and equivalent
   representations that give the same ALT haplotype.
2. **Tolerate read-side noise**: sequencing errors, BQ masking, and how the
   aligner writes the reads.
3. **Never infer a different ALT haplotype** for the counts.

The ≥90% insert identity band ([[identity-band-annotation-tolerance]]) sits on
this boundary. It is right only when the mismatches are read errors. Measure
it: scattered mismatches are errors; one recurring alternative insert is a
different allele or a wrong input, and must not count.

**Risk of the principle:** an honest low VAF for a mis-described allele misleads
anyone who ignores `gbcms_diagnostic`. So the caveat must be specific: say what
the reads carry, rather than only "partial".

**How to apply:** for any counting change, ask "does this credit a read to the
row's ALT when the read does not carry the given ALT haplotype?" If yes, it is
out, unless it is an opt-in rescue.
- Tolerating how *reads* are written (aligner representation, sequencing
  errors, truncated inserts in reads that end inside them) is fine; that is
  reading the given allele accurately.
- Tolerating a *wrong input* is not.
- Default-on features that report another allele's counts under the row's
  label (homopolymer decomposition, `WARN_HOMOPOLYMER_DECOMP`) conflict with
  this. They are open decisions in the 6.6.0 plan.
