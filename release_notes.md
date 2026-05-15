## [5.1.0] - 2026-05-11

### ⚠️ Breaking Changes

- **MAF column order reordered** (v5.1 schema):
  `any_alt`, `partial_alt`, `n_count` moved from end to immediately after
  `alt_count` for discoverability. Read strand counts (`ref_count_forward`, etc.)
  now precede derived strand bias statistics. Read and fragment metric layers
  fully separated (no interleaving). Downstream MAF parsers using positional
  indexing must be updated.
- **VCF FORMAT fields restructured** (VCF 4.2 spec compliance):
  - `DP` is now a single integer (total depth), was `ref,alt` pair.
  - `AD` is now `Number=R` with `ref,alt` totals (VCF spec), was `fwd,rev`.
  - `RD` and `RDF` removed. Replaced by `ADF` (forward strand `ref_fwd,alt_fwd`)
    and `ADR` (reverse strand `ref_rev,alt_rev`) following bcftools convention.
  - New `FAD`, `FADF`, `FADR` for fragment-level strand-by-allele counts.
  - `FAF` renamed from position after `VAF` to after fragment group.
  - Downstream VCF parsers expecting old `GT:DP:RD:AD:RDF:ADF:VAF:FAF:...` must
    be updated to `GT:DP:AD:ADF:ADR:VAF:FAD:FADF:FADR:FAF:...`.
- **VCF INFO field order changed**: `AAD`, `PAD`, `NAD` now appear immediately
  after `GR` (before strand bias fields), matching the diagnostic proximity
  principle.

### 🔧 Fixed

- **Wrong-length insertion Phase 3 fallback** (PAX5-class discordance fix):
  `check_insertion` now routes wrong-length insertions at the strict anchor
  position to Phase 3 (SW/PairHMM) for haplotype arbitration, mirroring
  `check_deletion`'s existing behavior. Previously, a read with `I(1)` at the
  anchor for an expected `I(2)` was silently classified as REF with no
  diagnostic signal. Now classified as `partial_alt` via `has_nearby_evidence`,
  triggering the `PARTIAL_DOMINANT` diagnostic flag.
- **Wrong-length insertion windowed scan tracking**: Wrong-length insertions
  within the ±window range now set `has_nearby_length_match`, routing to
  Phase 3 fallback. Previously silently ignored.
- **Insertion `!found_ref_coverage` haplotype fallback**: Added Phase 3
  fallback for insertion reads where the CIGAR walk found no evidence but
  the read spans the anchor (e.g., unusual CIGAR geometry, soft-clip at
  anchor). Mirrors `check_deletion`'s existing `!found_ref_coverage` path.
- **Wrong-length deletion windowed scan tracking**: Wrong-length deletions
  in the windowed scan now set `has_nearby_length_match` for Phase 3
  fallback. Covers two previously silent-drop cases:
    - Small deletions (≥5bp, <50bp): always flagged (1-4bp excluded as
      homopolymer noise)
    - Large deletions (≥50bp) with low reciprocal overlap (<50%):
      previously silently dropped, now flagged for Phase 3 arbitration

