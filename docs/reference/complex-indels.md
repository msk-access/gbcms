# Complex Indels

Real-world case studies of how gbcms handles the three most common complex indel classification scenarios. Each case shows the biological event, why a naive approach fails, and exactly which engine component resolves it.

> 📖 See [Allele Classification](allele-classification.md) for the complete algorithm reference.

---

## Overview

Complex indels fall into four categories that require special handling beyond the standard CIGAR walk:

| Case | Example | Root Cause | Fix |
|:-----|:--------|:-----------|:----|
| [Del+SNV routing](#case-1-complex-delsnv-sox9) | SOX9 `GC→T`, ABL1 `AG→T` | Deletion-format variant with anchor substitution sent to `check_deletion` | Route to `check_complex` when `alt[0] ≠ ref[0]` |
| [Clean-CIGAR REF reads](#case-2-nf2-large-deletion-ref-reads-invisible) | NF2 ~100bp delins (large deletion with a short replacement) | REF reads skipped by `is_worth_realignment()` → misclassified as 'neither' | M-block anchor coverage REF fallback in `check_complex` |
| [Left-alignment shift](#case-3-tp53-12bp-left-alignment-shifted-deletion) | TP53 12bp DEL | BWA anchor 3bp left of CIGAR `D` position → S3 validates wrong bases | `has_shifted_same_length` (del_len ≥ 5) → Phase 3 arbitration |
| [Wrong-length pure indels](#case-4-wrong-length-pure-indels-distinct-alleles) | Homopolymer 1bp-vs-2bp DEL; repeat-ladder DELs | Every tract-touching indel counted as the annotated event → VAF inflated several-fold | Wrong length → `partial_alt`, never REF/ALT; placement-aware ≥50bp band keeps split representations |

!!! note "Cases 1 and 2 are now judged by the exact-carrier rule"
    Del+SNV and delins variants no longer go through `check_complex`'s haplotype scoring or its
    M-block REF fallback. A read is REF or ALT only when its own bases carry that allele across
    the whole event, read at the event's own position; a large delins is judged at its two
    junctions. See [the exact-carrier rule](allele-classification.md#the-exact-carrier-rule).
    The case studies below record why the routing exists and what the earlier fixes found.

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

**Gene:** NF2 · **Position:** chr22:30038094 · **Type:** ~100bp deletion with a short (7bp) replacement — a large **delins** (`ref_len > alt_len > 1`), dispatched to `check_complex`

!!! note "Why this case is check_complex's problem, not check_deletion's"
    A PURE large deletion never reaches this fallback: its clean anchor-spanning REF reads
    are classified REF directly inside `check_deletion` (M coverage at the anchor with no
    D op is definitive). The invisible-REF failure below is specific to deletion-direction
    variants **dispatched to `check_complex`** — delins like this one, and Del+SNV.

### Why the Old Code Returned ref=0

The flow for deletion-direction delins in `check_complex`:

1. The variant dispatches to `check_complex` (`ref_len > alt_len > 1`)
2. `check_complex` Phase 1 reconstructs the read's sequence for the variant region
3. **`is_worth_realignment()`** — returns `false` for reads with clean M-only CIGARs
   (no I/D/S/N operations in the variant window)
4. Phase 3 is skipped for clean CIGARs
5. Old code: return `Neither` when Phase 3 is skipped → REF reads become 'neither'

REF reads for large deletion-direction variants naturally have clean M CIGARs — they don't
have a deletion because they support the reference. Skipping Phase 3 without a fallback
silently discards all REF evidence.

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
Variant: chr22:30038094  [~100bp REF] → [7bp ALT]  (large delins)
Anchor position: 30038094

REF read CIGAR: 150M  (spans anchor at 30038094)
    M-block: [30038040 ── 30038190]  ← covers anchor ✅, anchor BQ ≥ min_baseq ✅
    is_worth_realignment = false (clean CIGAR)
    Old: return Neither
    New: deletion-direction (ref_len > alt_len) + M coverage + anchor BQ → REF ✅

ALT read (deletion + replacement in CIGAR)
    is_worth_realignment = true → Phase 3 haplotype alignment → ALT ✅
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
has_shifted_same_length = true
```

After the CIGAR walk, if `has_shifted_same_length AND found_ref_coverage`:
- Route to Phase-3 arbitration (WFA+PairHMM under the default backend)
- Phase 3 builds correct haplotypes using the left-aligned anchor
- ALT reads align to the ALT haplotype → **ALT** ✅ (with `partial_alt`
  evidence propagated whenever Phase 3 does not confirm ALT)

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
    New: del_len(12) ≥ 5 → has_shifted_same_length = true
    → Phase 3: ALT haplotype aligns correctly → ALT ✅

REF read CIGAR: 150M  (no deletion)
    → M-block covers anchor → REF ✅
```

### Impact

| Metric | Before Fix | After Fix | Sign-out |
|:-------|:----------:|:---------:|:-------:|
| `ref_count` | 775 | 554 | — |
| `alt_count` | 0 | **221** | 222 (Δ=-1) |
| `dp` | 775 | 775 (unchanged) | — |

The Δalt=-1 is explained by one read landing in Phase 3's arbitration tie zone
(below the confidence threshold — the LLR ambiguity band under the default
`pairhmm` backend, or the score margin < 2 under `sw`) → counted neither.

---

## Case 4: Wrong-Length Pure Indels — Distinct Alleles { #case-4-wrong-length-pure-indels-distinct-alleles }

### The Scenario

A pure indel annotated in a repeat tract — a 2bp deletion in a homopolymer, or one rung of a
trinucleotide-repeat deletion ladder. The BAM carries **coexisting indel populations of
different lengths** at the same tract: for example, hundreds of reads with a 1bp germline
slippage deletion alongside a handful carrying the annotated 2bp somatic deletion.

### Why the Old Code Over-Counted

Wrong-length indels at the anchor were routed to Phase-3 haplotype arbitration. Phase 3's
haplotype window is **length-blind inside repeat tracts**: a read whose tract is one base
shorter than reference aligns almost equally well to a haplotype whose tract is two bases
shorter, and the log-likelihood ratio confidently picks ALT. Every tract-touching indel was
counted as the annotated event — inflating VAF several-fold at such loci (a locus with a
handful of true-ALT reads could report a VAF an order of magnitude too high). Separately, the
old ≥50bp "reciprocal overlap" tier accepted **any** deletion sharing half the annotated length
at the anchor (a `D(60)` counted for a 100bp deletion) because its sequence check compared
reference bases against reference bases — identically true by construction.

### The Fix (issue #91)

A read whose CIGAR proves a pure indel of a **different length** at the anchor is evidence of a
**different allele** → neither + `partial_alt`, never REF and never the queried ALT. Three
same-event escapes keep real support counted as ALT:

1. **Exact-length, sequence-verified** indels (unchanged).
2. **Placement-aware ≥50bp band**: the read deletes essentially the whole expected span
   (≤3 span bases retained, ≤3 changed outside) — covers breakpoint wobble and split
   `D+M+D` representations. Net-matching but *displaced* deletions (M ops across the span)
   are rejected.
3. **Insertion truncation containment**: an observed insert of ≥4bp, strictly shorter than
   expected, ≥90% identical to a slice of the expected insert, with **both** sequences
   non-low-complexity (sequencers truncate long insertions — the observed smear of shorter
   prefixes is the same event; 1–3bp fragments and repeat-only slices stay distinct-allele
   partial evidence).

Validated on read-level truth: the headline homopolymer locus went from counting 259 ALT reads
to the 7 that genuinely carry the annotated 2bp deletion (VAF 51% → 2.7%), with the 1bp-allele
population surfaced as `partial_alt`; a trinucleotide-ladder locus and an RNA locus landed
exactly on their read-truth counts; long-insertion loci kept their full ALT support via
containment.

### Reading the Output

```
alt_count    = reads carrying the annotated event (exact / in-band / truncation)
partial_alt  = distinct-allele evidence at the same tract
gbcms_diagnostic contains PARTIAL_DOMINANT when partial_alt > alt_count
```

A `PARTIAL_DOMINANT` pure-indel locus usually carries a coexisting allele the annotation does
not describe. Inspect per-read lengths with `--trace`.

---

## Diagnostic Commands

### Identify which checker handles a variant

```bash
# --trace shows per-read classification (check_deletion vs check_complex
# routing, wrong-length branches, Phase-3 fallbacks). RUST_LOG has no
# effect on gbcms — logging routes through Python's logging module.
gbcms dna --trace \
    --variants variants.maf \
    --bam sample.bam \
    --fasta ref.fa \
    --output-dir /tmp/debug/ \
    2>&1 | grep "check_deletion\|check_complex" | head -30
```

### Verify left-alignment shift

```bash
# See the normalized position vs original (writes a TSV with columns:
# chrom, original_pos, original_ref, original_alt, norm_pos, norm_ref,
# norm_alt, variant_type, gbcms_status, gbcms_status_reason,
# was_anchor_resolved, was_left_aligned, was_normalized)
gbcms normalize \
    --variants variants.maf \
    --fasta ref.fa \
    --output /tmp/norm/normalized.tsv

# Rows where left-alignment shifted the anchor (original_pos != norm_pos)
awk -F'\t' 'NR > 1 && $2 != $5' /tmp/norm/normalized.tsv
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
