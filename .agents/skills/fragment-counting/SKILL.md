---
name: fragment-counting
description: Reference for fragment-level counting in gbcms — R1/R2 quality-weighted consensus, FAD/FAF (fragment alt depth/fraction), QNAME/UMI deduplication, duplex vs singleton UMI families, strand consistency, and amplicon split. Use when touching shared/fragment.rs or any DPF/RDF/ADF logic, or debugging fragment-level VAF.
---

# Fragment-Level Counting

## Consensus model (R1/R2)

Fragment counting is **always on**. When the two mates of a fragment overlap the
locus, evidence is collapsed to one molecule:
- **Quality-weighted consensus**, not a read-count majority vote. Per allele,
  evidence is the max base-quality observation across mates.
- **Structural ALT wins**: if either mate carries a CIGAR I/D op matching the ALT
  indel, ALT wins (correct for indels). Caution: a single mate with a spurious
  indel CIGAR can override a high-BQ REF mate.
- **Discard band**: if mates disagree and neither allele's quality clears the
  threshold, the fragment is **discarded** — counted in DPF, excluded from RDF/ADF.
- `fragment_qual_threshold` (default 10) governs the discard band.

## Deduplication / hashing

- `hash_qname` → u64 (SipHash) keys a fragment; `hash_molecule` adds the UMI tag.
- The key is **QNAME(+UMI) only** — not orientation-aware (orientation is resolved
  downstream in `FragmentEvidence`). Collision probability is negligible at realistic
  per-locus depth.
- Overlapping mate pairs map to one `mol_hash` → counted once at fragment level.

## FAD / FAF and the invariants

```
assert dpf >= rdf + adf          # DPF includes discarded fragments
assert adf == adf_fwd + adf_rev  # fragment strand consistency
assert rdf == rdf_fwd + rdf_rev
```

## UMI families

- `singleton_alt_count`: ALT from singleton UMI families (no mate confirmation).
- `duplex_alt_count`: ALT from duplex families (both strands confirmed).
- `--umi-tag` selects the tag used for molecule hashing.

## Amplicon mode

`--library-type amplicon` **bypasses** R1/R2 consensus: mates are XOR-split into
separate fragments (both are PCR duplicates of one template). Disagreeing mates then
count as two independent fragments rather than discarding — FAF for the same physical
disagreement differs between capture and amplicon libraries by design.

## Read-level vs fragment-level (caution)

Supplementary/secondary records share a QNAME, so they collapse correctly at the
fragment level, but they increment **read-level** DP/RD/AD per record. With their
filters disabled this can double-count read-level depth at the anchor — see
read-filters-qc.

## Key Files
- `rust/src/shared/fragment.rs`: FragmentEvidence, `observe`/`resolve`, hashing
- `rust/src/counting/engine.rs`: fragment resolution loop, amplicon split
- `tests/test_fragment_consensus.py`, `tests/test_indel_fragment_consensus.py`
