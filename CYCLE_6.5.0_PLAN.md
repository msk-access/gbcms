# 6.5.0 cycle — implementation plan

> Branch-from: `develop` (post-#96). Every count-affecting ticket is red-first
> (xfail-strict battery committed before the fix), adversarially reviewed
> before commit, and validated on the local truth data named in its
> Evidence line. No new output columns anywhere — new signal goes to
> `gbcms_diagnostic`, logs, or validation tooling. Evidence links live on
> issues #92 / #94 / #97.

## T1 — Cluster exclusive assignment (count-affecting; the cycle's core)

**Problem (measured).** Co-annotated same-type indels sharing a repeat tract
are each counted independently; in repeat context a molecule carrying one
allele also passes the windowed S3 check of a tract-mate, so per-row `ad`
counts *compatible* molecules, not *assigned* ones (BRCA2 cluster: per-row
39–107 vs sign-out 21–45; Σad ≈ 397 vs ~150 distinct ALT molecules).
Adjudication proved exclusive canonical assignment reproduces sign-out
(28/28, 32/32, 21/21 exact; #92 comment of 2026-09-22).

**Design.**
1. *Grouping (prep, Python):* extend sibling detection with a second
   criterion — same-type variants whose scan windows overlap (window =
   `max(5, repeat_span+2)` each side, the engine's own formula). Groups are
   transitive closures. Emit group id into the existing sibling plumbing
   (`sibling_variants` already flows binned+legacy).
2. *Claiming (Rust, shared checker):* when a variant has same-type
   tract-mates, canonicalize each CIGAR indel candidate (left-align op
   against ref_context — same loop `normalize` uses) and compare canonical
   (pos, deleted-bases | inserted-bases) against the variant's own canonical
   form AND each tract-mate's:
   - matches self → structural ALT (unchanged);
   - matches a tract-mate → this read is that allele: for the current row it
     falls under the existing distinct-allele rule (`neither +
     has_nearby_evidence`, never REF, never ALT);
   - matches neither → current wrong-length / windowed logic unchanged.
   Windowed S3 shift-tolerance is therefore *scoped*: it may only rescue
   representations that canonicalize to SELF.
3. *No cross-variant runtime state needed* — each row decides from its own
   read walk plus the (already-passed) sibling list; binned/legacy parity
   holds by construction (shared checker).

**Files.** `src/gbcms/prep.py`/sibling grouping site; `rust/src/counting/
variant_checks.rs` (canonical helper + claiming in
`resolve_anchor_*_candidate` / `scan_windowed_*_candidate`);
`rust/src/normalize/` (reuse left-align); tests.

**Tests (red-first).** Synthetic BRCA2-like cluster: 3×D2 (AG/AC/AG) + D33 +
D14 interleaved in an AG tract, reads per allele → xfail-strict pins:
per-row ad = own molecules only; Σad ≤ distinct ALT molecules; tract-mate
carriers appear as `partial_alt`. Green guards: two *distant* same-length
dels (no window overlap) unaffected; a genuinely shifted representation of
SELF still rescued by S3; INS twin of the cluster. Parity case via
`count_both` with sibling groups (both paths).

**Acceptance.** Battery green; real-data: ACCESS harness rerun — BRCA2
cluster rows within a few molecules of sign-out (table from #92), all other
252 variants unchanged; FLT3/MSI/FORTE suites byte-identical (no cluster →
no change).

**Risks.** Canonicalization must use genomic ref_context (guaranteed since
B1); grouping too eager (window overlap on dense SNV+indel mix — group
same-type only); left-align loop bounds at ctx edges (reuse the guarded
normalize code, don't re-implement).

## T2 — ASJD-2: retention + novel-junction markers (#97)

**Problem (measured).** ASJD junction tallies go 0/0 exactly when a splice
variant's effect is total: (a) allele-specific intron retention (TP53 case:
440 ALT M-through, 1,696 spliced WT reads excluded as no-observation) and
(b) splice-destroying deletions expressed as novel junctions (B2M 360bp del
→ 41-read exon-2-skip junction among excluded reads).

**Design (diagnostic strings only, structural rules).**
1. *Retention-dominant marker* — at variants within splice regions
   (`exon_boundary_dist ≤ 2`, GTF-gated): if allele-classified reads with a
   junction = 0 while `splice_skip_excluded > dp` (excluded-spliced
   population dominates the classified one), append
   `RETENTION_DOMINANT(n_excluded)` to `asjd_diagnostic`. Population
   comparison, no tuned rate.
2. *Novel-junction marker* — deletion-type loci where
   `splice_skip_excluded > ad`: collect the junction spectrum of the
   excluded reads (the engine already walks their CIGARs in triage; add a
   bounded per-variant junction counter behind the existing GTF gate),
   test top junctions against `AnnotationIndex.splice_sites`; if a
   non-annotated junction with depth > ad exists whose donor or acceptor
   falls inside/adjacent to the deleted span, append
   `NOVEL_JUNC_AT_SPLICE_LOSS(n@donor-acceptor)`.
3. Both markers RNA+GTF only; absent otherwise (column gating rule).

**Files.** `rust/src/counting/engine.rs` (triage-adjacent junction counter;
ASJD assembly at the `detect_asjd` call site), `annotation/` lookup,
`types.rs` internal fields (not columns), docs `rna-annotation.md` +
`rna-splice-handling.md`.

**Tests (red-first).** Synthetic TP53-retention geometry (boundary SNV, ALT
M-through population + excluded spliced population) → RETENTION_DOMINANT;
synthetic B2M geometry (deletion destroying a donor + exon-skip N among
excluded reads) → NOVEL_JUNC_AT_SPLICE_LOSS; guards: ordinary exonic SNV
near junction (no marker), annotated-junction traffic (no NOVEL marker).
Acceptance on real data: the two FORTE cases produce their markers; the
33-sample cohort produces no spurious markers elsewhere.

## T3 — CLIP_CANDIDATES diagnostic (observability for clip-borne ITDs)

**Problem (measured).** Clip-only ITD carriers are invisible (±1bp Phase-3
clip window); prevalence survey: 12/12 loci I-op-dominant → rescue is a
sensitivity-tail case, deprioritized; analysts still need a pointer at the
affected zeros.

**Design.** At insertion-type loci with `ad == 0`: during the existing read
walk, count soft-clips `≥8bp` whose clip start lies within
`[anchor − (ins_len+10), anchor + (ins_len+10)]` (duplication reach —
structural bound from the insert's own length). If count ≥ 2, append
`CLIP_CANDIDATES(n)` to `gbcms_diagnostic`. Counter on BaseCounts
(internal), threaded like `splice_skip_excluded`; mirrored binned+legacy.

**Tests.** Red-first: the 30bp clip-only geometry (clips at duplication
boundary, zero I ops) → flag with n; guards: I-dominant locus (no flag even
with stray clips), SNV locus (never), `ad>0` locus (never). Real-data
acceptance: the P-0098981-T08 locus flags; the 12 survey loci do not.

## T4 — UMI-tag warn

**Design.** When `--umi-tag TAG` is set and, at the end of a sample's run,
zero processed reads carried TAG: one WARN per sample ("umi-tag TAG never
seen; fragment grouping fell back to QNAME"). Count during observe-side
extraction (a bool per worker, OR-reduced); no per-read logging. Warn, not
fail (mixed-BAM workflows legitimate).

**Tests.** BAM without RX + `--umi-tag RX` → warning in caplog, counts
unchanged; BAM with RX → no warning.

## T5 — dynamic_sw_gap_extend cleanup

**Evidence.** SW fallback fired zero times on traced real runs (ACCESS
duplex, MSI-high); function is arithmetically constant (−1).

**Design.** Delete `dynamic_sw_gap_extend`; use a named constant
`SW_GAP_EXTEND: i32 = -1` with a doc comment recording the measurement and
that CLI gap flags intentionally do not reach SW. Update
allele-classification.md SW-gap section + counting-engine skill. No
behavior change; clippy/tests must stay green untouched.

## T6 — Per-transcript BAQ decision

**Question.** Per-transcript counting applies BAQ unconditionally (no
exon-boundary suppression), making tx counts slightly conservative near
junctions vs the main counts (documented in code).

**Decision procedure (measure, then choose).** On 3 junction-rich FORTE
samples, diff per-transcript counts with/without suppression plumbed
(one-off instrumented build). If deltas are negligible → keep current
behavior, upgrade the code comment to a documented decision + doc note. If
material → small fix: thread `exon_boundary_dist` into
`count_per_transcript`'s BAQ call, mirroring the main loop, with a red
tx-count pin first.

## Order & discipline

T1 (own branch, own review) → T2 (own branch; needs the FORTE geometries) →
T3+T4+T5 (one observability/cleanup branch) → T6 (measurement decides).
Each branch: battery red → implement → suites + clippy + lint gate →
sonnet adversarial review → real-data acceptance named above → merge to
develop. 6.5.0 cut only after T1's ACCESS rerun and T2's cohort recheck.
