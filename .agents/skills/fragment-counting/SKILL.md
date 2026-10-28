---
name: fragment-counting
description: "Reference for fragment-level counting in gbcms — R1/R2 quality-weighted consensus, the DPF/RDF/ADF fragment counts (FAD/FAF in VCF), QNAME/UMI molecule hashing (--umi-tag), duplex/simplex merging, strand consistency, and amplicon split. Use when touching shared/fragment.rs or any DPF/RDF/ADF logic, or debugging fragment-level VAF."
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
- **Discard band**: if mates disagree and the better allele's quality does not beat the
  other's by more than the threshold, the fragment is **discarded** — counted in DPF,
  excluded from RDF/ADF.
- `fragment_qual_threshold` (`--fragment-qual-threshold`, default 10) governs the discard band.

## Deduplication / hashing

- `hash_qname` → u64 (SipHash) keys a fragment; `hash_molecule` adds the UMI tag.
- The key is **QNAME(+UMI) only** — not orientation-aware (orientation is resolved
  downstream in `FragmentEvidence`). Collision probability is negligible at realistic
  per-locus depth.
- Overlapping mate pairs map to one `mol_hash` → counted once at fragment level.

## Fragment counts (DPF/RDF/ADF; FAD/FAF in VCF) and invariants

```
assert dpf >= rdf + adf          # DPF includes discarded fragments
assert adf == adf_fwd + adf_rev  # fragment strand consistency
assert rdf == rdf_fwd + rdf_rev
```

## UMI molecules, duplex and simplex

- gbcms does not see UMI families, and cannot tell a duplex read from a singleton one.
  `--umi-tag` (e.g. `RX`) only changes the molecule key (`hash_molecule(qname, umi)` via
  `molecule_key` in engine.rs): reads that share a QNAME but carry different UMIs are
  different molecules.
- Duplex and simplex are separate BAMs, counted separately and combined afterwards:
  `gbcms merge` (merge.py) joins the per-BAM-type MAFs and adds the `simplex_duplex_*`
  columns (`bam_type` in the Nextflow samplesheet).

## Amplicon mode

`gbcms rna --library-type amplicon` (an RNA option; DNA has none) **bypasses** R1/R2
consensus: mates are XOR-split into separate fragments (both are PCR duplicates of one
template). Disagreeing mates then count as two independent fragments rather than
discarding — FAF for the same physical disagreement differs between capture and amplicon
libraries by design. It also forces `--strandedness unstranded`.

## Read-level vs fragment-level (caution)

Supplementary/secondary records share a QNAME, so they collapse correctly at the
fragment level, but they increment **read-level** DP/RD/AD per record. With their
filters disabled this can double-count read-level depth at the anchor — see
read-filters-qc.

## Key Files
- `rust/src/shared/fragment.rs`: FragmentEvidence, `observe`/`resolve`, hashing
- `rust/src/counting/engine.rs`: fragment resolution loop, amplicon split
- `tests/test_fragment_consensus.py`, `tests/test_indel_fragment_consensus.py`
