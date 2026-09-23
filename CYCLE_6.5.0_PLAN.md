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

**Design (as implemented — revised in-session; complex-variant support is
the key improvement, so claiming moved from a checker-level canonical
comparison to an engine-level sibling contest that covers delins).**
1. *Grouping (Rust prep, `assign_multi_allelic_groups`):* second membership
   criterion — **length-changing** variants (ref_len ≠ alt_len; never
   SNV/MNP, which would drain rd via the REF-side guard) whose scan windows
   overlap (window = `max(5, repeat_span+2)` each side, uncapped — the
   engine's own formula). Groups are transitive closures and may be
   non-contiguous in position order (visited-marker sweep; look-ahead bound
   derived from the input's max pad). Honest reason tags: true span
   intersection → `MULTI_ALLELIC`, window-only membership → `TRACT_CLUSTER`.
   Existing sibling plumbing (`sibling_variants`) carries the wider groups —
   zero new FFI.
2. *AD-claiming guard (Rust engine, `sibling_claims_alt`):* at a grouped
   locus, anchor-exact Phase 0 evidence is never contested; every other
   ALT faces three demotion tests, each decisive for a measured failure
   mode (haplotypes/reconstructions via `window_haplotypes` +
   `reconstruct_span`, the latter factored out of check_complex Phase 1):
   **(1) REF test** — over the read-covered context window the read must
   explain strictly better with the row's ALT haplotype than its REF
   haplotype (catches foreign flank events); **(2) probabilistic pure
   indel** — a Phase 3 ALT on an anchor-preserved pure indel row has no
   matching structural op anywhere and is ambiguity in a tract (catches
   unannotated ladder absorption: PairHMM prefers D14 over REF for a D11
   read); complex/MNP rows are exempt because Phase 3 is their carriers'
   normal path; **(3) strictly-better sibling** — over a window covering
   both spans, a sibling ALT haplotype with strictly lower cost claims
   the read; equal cost = equivalent representations (e.g. a delins
   double-annotated as an insertion) and both rows keep it. Earlier
   iterations falsified by real data: phase-rank gate (cross-attribution
   is Phase 3 on ABRA consensus reads; rank inverts on delins),
   per-sibling span-cost with tie-demote (anchor-exact siblings stole
   delins carriers — complex rows zeroed), ALT-vs-REF alone (same-tract
   same-length dels and ladders both beat REF). No sibling
   re-classification calls remain. Downgrade happens at classification
   time (before fragment evidence), so AD **and ADF** exclude demoted
   reads consistently; they surface as `partial_alt`/`any_alt`. Applied
   at binned, legacy, and per-transcript sites (the per-transcript site
   also gains the REF-side sibling guard). Reads not fully spanning the
   event, or splice-poisoned windows, keep their classification;
   isolated variants are untouched.
3. *REF-side guard symmetry:* sibling-claimed reads dropped from `rd` now
   surface as `partial_alt` too (distinct-allele evidence, same category as
   the wrong-length rule) instead of vanishing silently.

**Files.** `rust/src/normalize/engine.rs` (grouping rewrite + summary log);
`rust/src/counting/engine.rs` (claiming helper + 3 sites + REF-guard
partial); `src/gbcms/_rs.pyi` (reason-tag docs); tests.

**Tests (red-first, committed 81a761d/510ac71, flipped green).** Synthetic
four-deletion tract cluster (already-canonical, canonically distinct; D is
the window-only probe) + complex twin (pure del + delins CAG>T): pins
per-row ad = own molecules; Σad ≤ distinct ALT molecules; tract-mate
carriers surface as partial; delins sibling claims windowed carriers. Green
guards: distant same-sequence dels unaffected; shifted SELF representation
still rescued; MAF and VCF front-door agreement for both geometries (parity
oracle not used — siblings are parity-exempt).

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

**Design (as implemented — revised by measuring the excluded populations
on the three FORTE cases before coding).** Both markers live in
`detect_asjd` (it already re-classifies every read in the window), reading
a new excluded-population tally (reusing `JunctionTally`, fragment-deduped,
junctions over the locus only) — no new BaseCounts fields, no columns.
Measured populations (classified / excluded fragments): TP53 donor SNV
372 / 1,466 (dominant excluded junction annotated; 27 anchored novel);
RECQL acceptor SNV 63 / 303 (annotated 206; anchored novel exon-skip 50+36);
B2M exon-removing 360bp deletion 8,262 / 42 (dominant excluded junction is
the novel E1→E3 skip, 41).
1. *Gate (shared):* the variant's REF span reaches within two bases of an
   annotated intron boundary on the gene's strand
   (`intron_boundary_in_range(s-1, e+1, gene_strand)` — true donor/acceptor
   sites only; review found the first cut's `splice_sites` lookup admitted
   transcript termini and antisense genes' sites). Replaces the planned
   `exon_boundary_dist ≤ 2`: B2M's deletion starts 17bp from a boundary but
   contains both of exon 2's splice sites. Both markers also speak only above
   ASJD's own junction-evidence floors (10 REF-side / 5 ALT-side fragments).
2. *`RETENTION_DOMINANT(n)`:* excluded fragments > classified fragments AND
   classified fragments mostly junction-free. Describes observability (the
   genotyping reads are the retention population; `vaf` is retention-VAF),
   not causation — allele specificity is read from `vaf` (TP53 99%, RECQL
   53%).
3. *`NOVEL_JUNC_AT_SPLICE_LOSS(n@start-end)`:* applies to any variant type
   at the gate, not only deletions (the RECQL SNV shows the same exon-skip
   signal). Top unannotated junction on the excluded fragments that is
   anchored to an annotated site (±`JUNCTION_TOLERANCE`) and is not the
   deletion itself written as a splice (same length, at the locus) must be
   carried by more fragments than confirm ALT. The plan's "donor/acceptor
   inside the deleted span" test was dropped: an exon skip's endpoints lie
   outside the deletion by construction.
4. Both computed before ASJD's no-junction early return; RNA+GTF only.

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

## T5 — SW-under-PairHMM: keep-but-warn + diagnostic, and gap-extend cleanup

**Evidence.** SW has two roles: the explicit `--alignment-backend sw`
backend (kept — concordance tests and the cross-backend quality contract
depend on it), and a silent last-resort fallback under the PairHMM backend
when the pangenomic matrix cannot be built. The fallback fired zero times
on traced real runs (ACCESS duplex, MSI-high): prep always supplies
ref_context and post-B1 it is always genomic, so matrix construction cannot
fail on well-formed input. `dynamic_sw_gap_extend` is arithmetically
constant (−1).

**Design (operator decision 2026-09-22: keep-but-warn, log AND flag).**
1. The under-PairHMM fallback is retained but made loud: when it fires,
   (a) WARN once per variant naming the matrix-build failure reason
   (missing ref_context / offset out of bounds), and (b) count fallback
   reads in a new internal `BaseCounts::sw_fallback_reads` (threaded like
   `splice_skip_excluded`; NOT an output column; mirrored binned+legacy;
   `_rs.pyi` updated).
2. `_compute_diagnostics` appends `SW_FALLBACK(n)` to `gbcms_diagnostic`
   when `sw_fallback_reads > 0` — the row itself tells the analyst these
   counts came partly from a different scorer and upstream input was
   malformed.
3. Delete `dynamic_sw_gap_extend`; named constant `SW_GAP_EXTEND: i32 = -1`
   with a doc comment recording the measurement, the two SW roles, and that
   CLI gap flags intentionally do not reach SW. Update
   allele-classification.md SW-gap section + counting-engine skill.

**Tests.** Red-first: a variant engineered with broken context — as
implemented, a `ref_context` window that ends at the variant (an *absent*
context never reaches the fallback: both SW sites sit inside
`if let Some(ref_context)`); measured: the fallback's own SW also cannot
build haplotypes there, so those reads end as NEITHER — the flag is their
only trace under the PairHMM backend →
`SW_FALLBACK(n)` flag + caplog WARN, counts still produced; guards: normal
variants never flag (whole existing battery doubles as the guard); explicit
SW backend never flags (it is chosen, not fallen into). Parity: flag counter
equal via both paths on the synthetic case.

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
