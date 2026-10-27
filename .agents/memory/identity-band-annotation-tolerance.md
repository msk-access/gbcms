---
name: identity-band-annotation-tolerance
description: "The >=90% insert-identity band now applies only in the CIGAR fallback (insert_truncation_match, when the exact-carrier rule cannot judge an insertion); #278 D2 measured it for truncated inserts and did not adopt it. Standing preference kept: no new columns or flags - new signal goes to logs, diagnostics or validation tooling."
metadata: 
  node_type: memory
  type: feedback
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
  modified: 2026-10-27T00:00:00.000Z
---

**Superseded:** the 2026-09-22 "keep the band as primary, no BQ-first switch", by #278 D2
(operator, 2026-10-16): the band, tested for truncated inserts under the exact-carrier
rule, added 158 ALT reads on HG002 long insertions with truth concordance unchanged and,
on the panel, 13 ALT calls the bases oracle agrees with against 111 it does not. Recorded
in `docs/reference/read-judgment.md` (RJ-27 evidence) and on #143.

**Where the band applies now (checked in the code, 2026-10-27).** Only
`insert_truncation_match` (`rust/src/counting/variant_checks.rs`), called from
`resolve_anchor_insertion_candidate` inside `check_insertion`. `check_allele_with_qual`
reaches `check_insertion` for a pure insertion only when `by_bases` (the exact-carrier
rule) returns nothing: no prepared reference, or one without the event's flank. Such
reads set `carrier_fallback` and are counted and warned. The gates: observed insert >=4bp
and shorter than expected, >=90% identity to the best window, both sequences
non-low-complexity, at least one readable inserted base. Nothing else uses an identity
band (the engine, carrier, window and junction modules hold none). On every prepared
pure insertion a truncated copy is another allele (`partial_alt`), and one sequencing
error inside an exact-length insert is tolerated by RJ-26 (one substitution is an
error unless it recurs on 2 or more fragments and 5% of those fitting the allele), not by
a band.

**Why the band existed (history).** The 90% tolerance was meant to keep counting ALT when
the ALT string itself was imperfect (an insert off by a base, an adjacent variant folded
in): position-consistent, full-confidence mismatches shared across reads, which BQ
masking cannot absorb, while read-side errors are scattered and low-BQ. Since then
the project takes the input allele as correct (AGENTS.md invariant 7) and explains a
mismatch through `gbcms_status_reason` / `gbcms_diagnostic`.

**Measured before the change (2026-09-22/23, #92).** IMPACT raw reads at the three
covered >=4bp insertion loci (~900 carriers): zero-mismatch truncations ~24% of carriers
(219 candidates); mismatching candidates ~2%, mostly Q20-35 singletons at unique
positions; all-low-BQ mismatches 6 reads. ACCESS consensus reads at the two band-eligible
loci: exact, length truncations (2-14 per locus), N-carrying (0-10 duplex, up to 83
simplex) and zero same-length mismatches, so a BQ-first policy has nothing to act on
there. Numbers and links: `docs/reference/bam-evidence-caveats.md` (sequencing errors
inside or beside the event).

**Still valid - how to apply.** No new CLI flags or output columns for a change here
(standing rule): new signal goes to logs, existing diagnostics or validation tooling.
When touching insert matching, test raw identity against BQ-masked policies on the
divergence classes. Panels are regression validation for one assay, never universal
calibration ([[tolerant-deletion-deliberate]]).
