# Complex Indels

Real-world case studies of how gbcms handles the three most common complex indel classification scenarios. Each case shows the biological event, why a naive approach fails, and exactly which engine component resolves it.

> 📖 See [Allele Classification](allele-classification.md) for the complete algorithm reference.

---

## Overview

Complex indels fall into three categories that require special handling beyond the standard CIGAR walk:

| Case | Example | Root Cause | Fix |
|:-----|:--------|:-----------|:----|
| [Del+SNV routing](#case-1-complex-delsnv-sox9-gcт) | SOX9 `GC→T`, ABL1 `AG→T` | Deletion-format variant with anchor substitution sent to `check_deletion` | Route to `check_complex` when `alt[0] ≠ ref[0]` |
| [Clean-CIGAR REF reads](#case-2-nf2-large-deletion-ref-reads-invisible) | NF2 ~100bp DEL | REF reads skipped by `is_worth_realignment()` → misclassified as 'neither' | M-block anchor coverage REF fallback in `check_complex` |
| [Left-alignment shift](#case-3-tp53-12bp-left-alignment-shifted-deletion) | TP53 12bp DEL | BWA anchor 3bp left of CIGAR `D` position → S3 validates wrong bases | `has_nearby_length_match` (del_len ≥ 5) → Phase 3 SW |

---

## Case 1: Complex Del+SNV — SOX9 (`GC→T`) { #case-1-complex-delsnv-sox9 }

### The Variant

**Gene:** SOX9 · **Position:** chr17:70119766 · **REF:** `GC` · **ALT:** `T`

This is a deletion-format variant (`len(REF)=2 > len(ALT)=1`), but it is **not** a pure deletion:
the anchor base G also changes to T. This makes it a complex Del+SNV — one net base is deleted
AND the anchor substitutes simultaneously.

### Why the Old Code Returned alt=0

The dispatch logic routed all `N×1` variants to `check_deletion`. `check_deletion`'s CIGAR walk:
1. Finds a `D(1)` at the anchor position in reads carrying the deletion
2. Applies **Safeguard S3**: checks whether the anchor base in the read matches the expected anchor (`G`)
3. The ALT read has anchor base `T` (not `G`) — S3 rejects it
4. `found_ref_coverage = true` → classified as **REF**

Result: all DEL reads classified as REF → `alt = 0`.

### The Fix (Fix 2 — `engine.rs`)

Before routing `N×1` variants, `engine.rs` now checks `alt_allele[0] != ref_allele[0]`:

```
REF = GC  →  ref[0] = G
ALT = T   →  alt[0] = T
Since G ≠ T → route to check_complex (not check_deletion)
```

`check_complex` builds `REF_haplotype = ...GC...` and `ALT_haplotype = ...T...` and runs Phase 3
Smith-Waterman. ALT reads align better to `ALT_hap` → classified as **ALT**.

### Read-Level Example

```
Variant: chr17:70119766  GC → T
Reference:  ... A  G  C  T  A  T ...
                   ↑
               anchor pos (G)

ALT read CIGAR: 5M 1D 4M
    ... A [T] ── T  A  T ...
         ↑  ↑
         T  D(1): anchor substitutes G→T, C is deleted
    → check_deletion S3: base at anchor = T ≠ G → REJECTED (old)
    → check_complex Phase 3: ALT haplotype wins        → ALT ✅ (new)

REF read CIGAR: 9M
    ... A  G  C  T  A  T ...
    → check_complex Phase 3: REF haplotype wins        → REF ✅
```

### Affected Variants

The same fix handles any Del+SNV with this anatomy:

| Gene | REF | ALT | del_len | anchor change |
|:-----|:----|:----|:-------:|:-------------:|
| SOX9 | `GC` | `T` | 1 | G→T |
| ABL1 | `AG` | `T` | 1 | A→T |

---

## Case 2: NF2 Large Deletion — REF Reads Invisible

### The Variant

**Gene:** NF2 · **Position:** chr22:30038094 · **Type:** ~100bp deletion

### Why the Old Code Returned ref=0

The flow for large deletions reaching `check_complex`:

1. `check_deletion` routes to `check_complex` when no windowed CIGAR match is found
2. `check_complex` Phase 1 reconstructs the read's sequence for the variant region
3. **`is_worth_realignment()`** — returns `false` for reads with clean M-only CIGARs
   (no I/D/S/N operations in the variant window)
4. Phase 3 SW is skipped for clean CIGARs
5. Old code: return `Neither` when Phase 3 is skipped → REF reads become 'neither'

REF reads for large deletions naturally have clean M CIGARs — they don't have a deletion
because they support the reference. Skipping Phase 3 without a fallback silently discards all REF evidence.

### The Fix (Fix 3 — `variant_checks.rs` `check_complex`)

After `is_worth_realignment()` returns `false`, gbcms now checks:

```
if ref_len > alt_len           # deletion-direction variant
   AND any M-block covers anchor_pos:
   → classify as REF
```

The M-block coverage check confirms the read spans the anchor position, establishing REF support
without needing haplotype alignment.

### Read-Level Example

```
Variant: chr22:30038094  A[100bp]→A (100bp deletion)
Anchor position: 30038094

REF read CIGAR: 150M  (spans anchor at 30038094)
    M-block: [30038040 ── 30038190]  ← covers anchor ✅
    is_worth_realignment = false (clean CIGAR)
    Old: return Neither
    New: M-block covers anchor AND ref_len(101) > alt_len(1) → REF ✅

ALT read CIGAR: 50M 100D 50M  (carries the deletion)
    Deletion at anchor → classified via windowed scan or Phase 3 → ALT ✅
```

### Impact

| Metric | Before Fix | After Fix |
|:-------|:----------:|:---------:|
| `ref_count` | 0 | 56 |
| `alt_count` | 49 | 49 (unchanged) |
| `dp` | 49 | 105 (correct) |

---

## Case 3: TP53 12bp Left-Alignment Shifted Deletion

### The Variant

**Gene:** TP53 · **Position:** chr17:7579309 · **REF:** `GACCGTGCAAGT` (12bp) · **ALT:** `-` (deletion)

### Root Cause: BWA Left-Alignment

BWA left-aligns indels to the leftmost equivalent position in the genome. For this 12bp deletion in
a region with a partial repeat (`ACC` motif), BWA places the anchor 3bp **left** of where the
CIGAR `D(12)` operation actually appears in the reads.

```
Left-aligned anchor position (what MAF/gbcms sees):  chr17:7579309
Actual CIGAR D(12) position in reads:                chr17:7579312  (+3bp right)
```

### Why the Old Code Returned alt=0

`check_deletion` windowed scan:
1. Finds `D(12)` within the ±5bp window (at position 7579312)
2. Passes **S1**: length matches (12bp)
3. Applies **S3**: compares reference bases at `7579312` against expected deleted sequence `GACCGTGCAAGT`
4. Reference at `7579312` is `CGTGCAAGT...`, not `GACCGTGCAAGT` → **S3 fails**
5. Old code: S3 fail → Continue → read not counted → `found_ref_coverage = true` → REF

All 221 ALT reads fail S3 → `alt = 0`, `ref = 775`.

### The Fix (Fix 4 — `variant_checks.rs` `check_deletion`)

When a windowed `D(12)` matches in **length** but fails S3, and `del_len ≥ 5`:

```
has_nearby_length_match = true
```

After the CIGAR walk, if `has_nearby_length_match AND found_ref_coverage`:
- Route to `check_complex` for Phase 3 SW arbitration
- Phase 3 builds correct haplotypes using the left-aligned anchor
- SW correctly aligns ALT reads to the ALT haplotype → **ALT** ✅

### Why del_len ≥ 5?

Short (1–4bp) same-length Dels failing S3 are almost certainly **unrelated spurious deletions**
in the wrong reference context — not left-alignment artifacts. CIGAR-definitive REF is correct
for these cases.

Large deletions (≥5bp) failing S3 are more likely to be genuine left-alignment shifts: BWA's
left-alignment can move the anchor up to `repeat_span` positions left, and 5bp is the minimum
spanning a small repeat motif (e.g., dinucleotide repeat `ACAC`).

### Read-Level Example

```
Variant (left-aligned): chr17:7579309  GACCGTGCAAGT → -
Anchor (gbcms):          pos = 7579309

ALT read CIGAR: 48M 12D 52M
    D(12) at pos 7579312  (3bp right of anchor)
    S1: length 12 == 12 ✅
    S3: ref[7579312..7579324] = "CGTGCAAGT..." ≠ "GACCGTGCAAGT" ❌
    Old: Continue → REF (incorrect)
    New: del_len(12) ≥ 5 → has_nearby_length_match = true
    → Phase 3 SW: ALT haplotype aligns correctly → ALT ✅

REF read CIGAR: 150M  (no deletion)
    → M-block covers anchor → REF ✅
```

### Impact

| Metric | Before Fix | After Fix | Sign-out |
|:-------|:----------:|:---------:|:-------:|
| `ref_count` | 775 | 554 | — |
| `alt_count` | 0 | **221** | 222 (Δ=-1) |
| `dp` | 775 | 775 (unchanged) | — |

The Δalt=-1 is explained by one read landing in the ties (Phase 3 score margin < 2).

---

## Diagnostic Commands

### Identify which checker handles a variant

```bash
# --trace shows check_deletion vs check_complex routing per read
RUST_LOG=trace gbcms dna \
    --variants variants.maf \
    --bam sample.bam \
    --fasta ref.fa \
    --output-dir /tmp/debug/ \
    2>&1 | grep "7579309\|check_deletion\|check_complex" | head -30
```

### Verify left-alignment shift

```bash
# See the normalized position vs original
gbcms normalize \
    --variants variants.maf \
    --fasta ref.fa \
    --output-dir /tmp/norm/ \
    --show-normalization

# If norm_pos ≠ Start_Position, left-alignment shifted the anchor
grep "TP53" /tmp/norm/normalized.maf | cut -f5,6,11,12  # pos, ref, norm_pos, norm_ref
```

### Check CIGAR D position vs anchor

```python
import pysam

bam = pysam.AlignmentFile("sample.bam")
anchor = 7579308  # 0-based

for read in bam.fetch("chr17", anchor - 10, anchor + 15):
    for op, length in read.cigartuples:
        if op == 2:  # D (deletion)
            print(f"D({length}) at ref_pos ~ {read.reference_start}")
```

---

## Related

- [Allele Classification](allele-classification.md) — Full dispatch and algorithm reference
- [Variant Normalization](variant-normalization.md) — How left-alignment works
- [Troubleshooting](../resources/troubleshooting.md) — Common count issues and diagnosis steps
- [Read Filters](read-filters.md) — Which reads reach the allele checker
