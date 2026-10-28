---
name: tolerant-deletion-deliberate
description: "Wrong-length pure indels are another allele (partial_alt), never ALT. Complex events go through the exact-carrier rule (classify_complex), 50+ bp deletions through the junction judge; Phase 3 decides essentially no reads (#265). The older tolerant length band and 'delins stay Phase-3' are superseded."
metadata: 
  node_type: memory
  type: project
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
  modified: 2026-10-27T00:00:00.000Z
---

**Superseded:** the 2026-09-17 "delins stay Phase-3, do not touch" and the placement-aware
50+ bp band, by RJ-4 (#141), RJ-15 (#202), RJ-24 (#254, #285) and #278 D2.

**Current rule (checked against the code, 2026-10-27).** Dispatch is
`check_allele_with_qual` in `rust/src/counting/engine.rs`; a row's route follows
`window::allele_kind`, never the `variant_type` label.
1. **Small pure indels in repeats** (regime 1): a read carrying another length is another
   allele, `partial_alt`, never ALT or REF (RJ-6). Pure insertions and deletions under
   50 bp are judged by their bases (RJ-27).
2. **Complex events** (regime 2: delins, anchor-changing indels, MNP reads with an indel
   at the block): `classify_complex` counts a read only when its own bases carry the
   whole allele with two flank bases (RJ-4). `check_complex` and Phase 3 run only where
   the rule cannot judge the row (no prepared reference, or one without the event's
   flank); those reads set `carrier_fallback` and are counted and warned once per
   variant. Phase 3 decided 0 of ~4M read calls on every runnable arm (#265; the figure
   is in CONTINUITY.md, not re-measured here). A read whose own deletion covers the
   anchor is judged by its bases between its flanks, never Phase 3's closer haplotype
   (RJ-15).
3. **Large pure indels** (regime 3). Deletions of 50 bp or more: `junction::classify`
   reads the POS breakpoint, however the aligner wrote the gap (RJ-24); the 3-base
   placement band is gone. Insertions: a truncated copy is another allele under the
   exact-carrier rule (#278 D2); the CIGAR fallback's 90% containment band is the only
   place a truncation still counts ALT (see [[identity-band-annotation-tolerance]]).

**Why real data stays the guide:** read-level validation for #91 (2026-09-17) measured
nine real pure deletions (12-539bp) across five full clinical BAMs: every supporting read
is an exact-length single `D` op, zero +-1-3bp wobble, zero splits. That is why a wrong
length is another allele and no length tolerance is kept.

**Why the 50bp threshold (operator-confirmed):** an artifact-SIZE prior, not an
event-rarity claim. Slippage, stutter and alignment artifacts make SMALL spurious indel
ops (overwhelmingly <50bp, mostly 1-few bp in repeat context); a read essentially never
acquires a >=50bp D by artifact, so a big observed op is trusted as a real deletion in
that molecule. Below 50bp no length tolerance is safe (the <5bp windowed noise gates of
the CIGAR fallback encode the same prior). `LONG_DELETION = 50` in `junction.rs` now
selects the junction judge.

Related: [[pysam-validation-oracle]].
