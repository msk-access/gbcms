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
acceptance: the clip-only ZERO_ALT ITD locus (t_alt=7) flags; the 12 survey loci do not.

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

**Tests.** Red-first: a variant engineered with broken context (e.g.
ref_context absent via direct `_rs` call) under the PairHMM backend →
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

## T7 — MNP rescue: gate matches intent, rescued row is one coherent genotype (opt-in path only)

**Intent (v4.3.0, unchanged).** `--rescue-mnp` exists for sign-out MNPs
whose carriers hold only a component of the annotated haplotype: gbcms
correctly reports the annotated allele as absent (carriers land in
`partial_alt`), sign-out reports the component count. Rescue re-counts each
discriminating position as an SNV and reports the best component — the
architecture doc's worked example is the TERT promoter GAGGG>AAGGA itself.
Opt-in; the flag's contract *is* component counting.

**Problem (measured, T1 acceptance).** Isolated TERT GAGGG>AAGGA
(disc 2/5): ad 1, partial 88, rd 486 vs sign-out t_alt 93 / t_ref 484.
Read-verified: 90 reads carry only the first G>A (C250T), 1 only the second,
0 the full haplotype — engine exact (pinned: `TestONPCarrierShapes`). Rescue
did not fire because the gate is `ad == 0`, and masked per-position
evaluation makes that brittle: one C250T read with BQ<20 at the second
position votes on the first alone → full ALT → ad 1. Gate intent is
"haplotype effectively absent"; code implements "exactly zero".

Four further defects in the rescue pass (verified end-to-end on synthetic
BAMs):
- **Stale row.** `with_ad()` replaces `ad` only. A rescued row reported
  alt_count 9 beside alt_count_forward/reverse 0, alt_count_fragment 0,
  `ZERO_ALT` — violating `ad == ad_fwd + ad_rev` in the written output; VCF
  AD/VAF move while ADF/ADR/FAD/FAF, strand bias, mFSD, any_alt/partial_alt
  stay at MNP values. Fragment columns are never rescued, so ACCESS
  (scored on fragments) and `merge --add-combined` never see rescue.
- **Cross-sample leak.** `prepared` is shared across all BAMs of a run and
  `gbcms_rescue` is never reset: a later sample with 10 genuine full-haplotype
  reads printed the earlier sample's rescue string.
- **T1 bypass.** Synthetic SNVs are counted with no siblings; a grouped MNP
  that lost reads to a sibling (now surfaced as `partial_alt`) would get them
  back — breaks Σ per-row ad ≤ distinct ALT molecules. Reachable today when a
  grouped MNP loses every read (ad 0).
- **Audit hard-codes `original_alt=0`.**

**Status (2026-09-23).** Implemented on `feature/mnp-rescue-gate`: red battery
041c733 → fix 3844052 → docs 3a06dc0. Real-data acceptance (below) pending.

**Design (Python orchestration; the only Rust/stub change deletes
`BaseCounts.with_ad`, which has no callers left).**
1. *Gate:* candidate when PASS, MNP, `MNP_RESCUE_ELIGIBLE`, and
   `partial_alt > ad` (the existing `PARTIAL_DOMINANT` condition), replacing
   `ad == 0`. Population comparison, no tuned rate; ad 0 with partial > 0 is
   still covered.
2. *Replace only if the component beats the haplotype:* best SNV `ad > ad`,
   else `outcome=no_improvement` and keep the MNP counts. The masked-full-
   carrier case needs partial ≤ ad, so the gate excludes it — but review
   found (and a synthetic battery confirmed) another source: reads with an
   indel inside the block go to the complex path, count REF with
   nearby-indel evidence (`partial_alt`), and no single-base count calls
   them ALT. Legitimate decline, logged at DEBUG; decision logic lives in
   the pure `_resolve_mnp_rescue` (unit-tested).
3. *Adopt the winning SNV's whole `BaseCounts`,* not an ALT-side graft. Every
   count, fragment, strand, strand-bias, mFSD and RNA column then comes from
   one real counting pass, so all four counting invariants hold on rescued
   rows by construction and Invariant-1 breakage disappears. Rejected
   alternative — graft ALT fields onto the MNP record (needs a Rust copy
   method + stats recompute): mixes two classifications, and fragment
   consensus lets a neither-read abstain, so an R1-partial/R2-masked-REF
   fragment is MNP-REF *and* SNV-ALT → can break `dpf ≥ rdf + adf`.
   Consequence to document: rescued rd/dp are the SNV's (reads carrying only
   another component count REF at the winning position; reads covering the
   position but not the whole block count in dp).
4. *Audit carries the MNP forensics* (they no longer live in the count
   columns): `method=decomposed;outcome=O;original_ref=R;original_alt=A;original_partial=P[;adopted=chr:pos(R>A)][;positions=chr:pos(R>A):<ad|ref_fail>,...]`
   with outcomes rescued / skipped_grouped / no_improvement /
   ref_validation_failed — `original_alt` keeps its key (now the true value,
   not 0); a component whose synthetic SNV fails preparation reads
   `ref_fail` (was a silent 0).
5. *Diagnostics describe the row as written:* recompute `gbcms_diagnostic`
   for rescued rows from the adopted counts (keeps `MNP_DISC_RATIO`/
   `MNP_RESCUE_ELIGIBLE`); a non-empty `gbcms_rescue` is the record of why.
   Removes the ZERO_ALT-beside-alt_count>0 contradiction.
6. *Skip grouped rows* (`multi_allelic_group` set): `outcome=skipped_grouped`;
   exclusive assignment owns those reads.
7. *Reset `gbcms_rescue` per sample* before candidate selection (fixes the
   leak).
8. *mFSD Parquet (found in review):* a rescued row's `--mfsd-parquet`
   record keeps the MNP's coordinates with the adopted component's
   fragment sizes (consistent with the row's mFSD columns; no new column —
   documented and logged per sample, join to `gbcms_rescue`).
9. *Found in pre-implementation review:* both count calls now share
   `_engine_kwargs()` (the component re-count had a verbatim copy of the
   main call's ~25 engine settings); diagnostics split into per-row
   `_diagnostic_flags` so rescued rows reuse the same code; unused `snp_map`
   and the quadratic per-position index lookup removed.

**Files.** `src/gbcms/pipeline.py` (`_rescue_mnp_pass`, per-sample reset,
diagnostics recompute for rescued rows); `rust/src/types.rs` + `_rs.pyi`
(`with_ad` removed); `cli.py` help + `models/core.py`
field descriptions (both commands); tests; docs `cli/dna.md`,
`nextflow/parameters.md`, `reference/architecture.md` (candidate table,
invariant-impact table → "invariants hold"), `reference/output-formats.md`
(`gbcms_rescue` format), `development/developer-guide.md`; `mnp-rescue`
skill; CHANGELOG (behavior change on the opt-in path + `gbcms_rescue`
format change for downstream parsers).

**Tests (red-first).** Battery on the TERT geometry via the CLI:
(a) component carriers + one masked stray full read → rescued, row equals
the winning SNV's counts, four invariants asserted on the written row,
audit carries original ref/alt/partial; (b) cis carriers (ad > partial) →
untouched; (c) indel-disrupted carriers (1bp insertion inside the block)
→ `no_improvement`, counts untouched, plus outcome-resolution unit tests
(ref_validation_failed and leftmost tie-break are unreachable through a BAM);
(d) grouped MNP → `skipped_grouped`; (e) two BAMs in one run → second
sample's `gbcms_rescue` reflects only itself; (f) rescue off → output
byte-identical to today. Replace `test_rescue_skips_nonzero_ad` with the
new gate's pin. Carrier-shape engine pin already written
(`TestONPCarrierShapes`, test_mnp_concordance.py). Parity oracle unaffected
(Python post-pass).

**Acceptance.** TERT row with `--rescue-mnp`: alt ≈ 88–90 vs sign-out 93,
audit shows the 1295250/1295254 split; the sample's seven dinucleotide MNPs
unchanged (ad > partial); complex-cluster, ACCESS, FLT3/MSI/FORTE harnesses
byte-identical with the flag off, and with it on change only rows whose audit
says rescued.

**Risks.** Germline-component MNPs (somatic change merged with a het SNP)
report the germline component — inherent to the flag since v4.3.0, now
reachable in tumors where the haplotype is present but dominated; document
in the flag help. `PARTIAL_DOMINANT` for MNPs also counts third-allele
reads that match ALT at ≥1 position — noise-level at real depth, and item 2
still requires the component to beat the haplotype.

## Order & discipline

T1 (own branch, own review) → T2 (own branch; needs the FORTE geometries) →
T3+T4+T5 (one observability/cleanup branch) → T6 (measurement decides).
T7 is independent (Python-only, opt-in path) — own small branch
(`feature/mnp-rescue-gate`), can land anytime; its leak fix may go first.
Each branch: battery red → implement → suites + clippy + lint gate →
sonnet adversarial review → real-data acceptance named above → merge to
develop. 6.5.0 cut only after T1's ACCESS rerun and T2's cohort recheck.
