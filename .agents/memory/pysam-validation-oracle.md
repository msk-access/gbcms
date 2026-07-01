---
name: pysam-validation-oracle
description: For cross-validating gbcms counts against a BAM, use pysam fetch()+get_reference_positions, not pileup() (which under-counted ~10% here).
metadata:
  type: reference
---

When independently checking gbcms counts against the raw BAM:

- **Use `bam.fetch(chrom, pos0, pos0+1)` + `read.get_reference_positions(full_length=True)`**
  to locate the read base aligned to a reference position. `bam.pileup()`
  under-counted depth by ~10% in a real ctDNA duplex BAM **even with `max_depth`
  raised to 1e7** (an artifact of pileup's stepper/overlap handling), which made
  gbcms look like it over-counted depth — it did not.
- **Match gbcms `dna` filters:** MAPQ≥20, BQ≥20, exclude flags
  `0x4|0x100|0x200|0x400|0x800` (unmapped/secondary/qcfail/dup/supplementary);
  improper-pair and indel are NOT filtered by default.
- **`DP = RD + AD + neither`** ([[bin-anchor-coverage]] sibling invariant): an
  independent `RD+AD` over high-BQ M-aligned reads matches gbcms exactly; gbcms's
  DP is higher by the *neither* reads (low-BQ / soft-clipped at the locus), which
  it correctly counts in depth. So compare RD and AD for exactness; expect DP ≥ a
  naive aligned-read count.

**Result (2026-06-27, ctDNA duplex, 571 SNPs):** ALT 571/571 exact, REF 571/571
within ±1%, DP 90% within ±1% (rest = neither reads). Reusable check:
`scratchpad/xcheck.py` (local, uses patient data — never committed).
