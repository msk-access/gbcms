---
name: genotyper-not-caller
description: "Operator principle — gbcms is a genotyper, not a caller: it counts reads whose own bases carry the given allele; the alignment can be wrong, so judge bases, not placement; sign-out is a comparison. Validate per read across data types."
metadata:
  type: feedback
---

gbcms is a **genotyper, not a caller**. Given an allele, it reports how many reads
carry it (and REF, and what else is there as diagnostics); it does not decide
whether a variant exists, discover alleles for a row, or rewrite a wrong input
([[count-the-given-allele]]).

The BAM is **not** the truth either: alignments can be wrong (an indel placed
elsewhere in a repeat, an allele soft-clipped at a read end, a long event split
across alignments, a misaligned read). The evidence is each read's **own bases**,
judged independently of how the aligner represented them. That is why the counting
rules anchor windows on flank bases, grow events through repeats, read clipped
bases and treat junctions and splices explicitly, rather than trusting the CIGAR.

Sign-out alleles are one more result to compare against, and they can be wrong
(manual edits, shifted or garbage alleles). A recurrent unannotated haplotype is a
finding for validation and diagnostics, never ALT evidence for the row.

**Why:** operator, 2026-09-23 (reads, not sign-out, are what we score against),
2026-09-25 (sign-out alleles can be wrong; people use IGV), and 2026-09-28: "BAM is
not necessarily the truth, it can have alignment errors, but the point is we are
not a caller but a genotyper"; and "always try to do things as generalized as
possible; use different datasets to make the architecture general."

**How to apply:** validate count-changing work per read against an independent,
position-aware census of the reads' own bases (clipped bases included,
equal-length windows, masked low-quality bases), across data types before merge:
synthetic contracts and fuzz with varied base-quality encodings; realigned panel
DNA (tumour, normal, duplex, simplex); DNA without indel realignment (WES/WGS);
RNA (spliced, BAQ); public reference data when available. Cross-platform agreement
is secondary. When the right behaviour is not known, measure it on that matrix and
survey community practice (callers/genotypers, with sources) before deciding (operator,
2026-09-28); the caveats are in `docs/reference/bam-evidence-caveats.md`. Treat any feature that reports a different allele under a row's
label (e.g. `--rescue-mnp`) as opt-in, flagged and audited. The concrete local
datasets are in a local-only memory. See [[pysam-validation-oracle]] and
[[holistic-effects-map]].
