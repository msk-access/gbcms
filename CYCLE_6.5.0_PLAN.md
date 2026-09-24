# 6.5.0 cycle — implementation plan

> Branch-from: `develop` (post-#96). Every count-affecting ticket is red-first
> (xfail-strict battery committed before the fix), adversarially reviewed
> before commit, and validated on the local truth data named in its
> Evidence line. No new output columns anywhere — new signal goes to
> `gbcms_diagnostic`, logs, or validation tooling. Evidence links live on
> issues #92 / #94 / #97 (T9/T10: #106 / #107).

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

**Measured (material → fixed).** 3 junction-rich RNA samples, 2010 probes
0–4bp inside annotated exon edges plus 397 at 25bp, two builds of one commit:
- Main counts: identical.
- Per-transcript RD: +3.0% at edge probes (up to 3x for single pairs).
- Median gap between best-transcript RD and main RD: 2.5% → 0.05%
  (p90: 6.4% → 1.45%).
- The recovered reads' mismatch rate is ~Q30, and the best-transcript
  mismatch rate afterwards equals the main counts'.
- ASJD: 5800 rows changed. Some edge variants had whole known junctions,
  thousands of fragments, invisible.

**As implemented.**
- One rule, `baq_applies` (+ `BAQ_BOUNDARY_SUPPRESS_BP`). The main counts
  resolve it once per variant. The per-transcript and ASJD passes resolve it at
  their own variant's position and receive it as `use_baq`.
- A DEBUG line names each skipped variant.
- The measurement surfaced a second defect: ASJD took each partition's
  dominant junction in hash order. Ties flipped junction and p-value between
  runs of the same input, which also moved 4 "control" probes. Now
  `top_junctions` breaks ties: when REF's and ALT's tied-top sets overlap,
  both report the shared junction and no test is run; otherwise the leftmost.

Review round:
- The tie match is tolerance-aware. One `same_junction` predicate (±5bp) now
  serves the tie-break and Step 4: an exact-only match could still call a
  2bp-shifted tie a divergence.
- `nearest_splice_distance` normalizes its contig, as every other annotation
  lookup does, and returns `Option`. `chrM` against an `MT` GTF and
  unannotated contigs read 2147483647 before.
- The per-transcript and ASJD passes reuse the main counts' distance.
- Non-shared ties are broken by evidence (`most_supported`): each allele takes
  the tied junction the other allele uses most, then the leftmost. A leftmost-
  only rule let coordinates flip the flag, e.g. p 0.003 vs 0.16 on the same
  reads.
- Shared RNA fixtures moved to `tests/rna_fixtures.py`.
- Pre-existing and out of scope, flagged separately: the BAQ rule and
  `exon_boundary_dist` key on `pos`, not the variant span; the decomposed twin
  never receives `gene_strand`.

## T7 — MNP rescue (`--rescue-mnp`, opt-in): report only what the BAM shows

**Principle (operator, 2026-09-23).** The BAM is the truth; sign-out is one
more result to compare against. Default output (rescue off) reports the
annotated MNP exactly as the reads support it and is unchanged by T7.

**What rescue is for.** Some signed-out MNPs are really one SNV (or two SNVs
on different molecules) annotated as a single multi-base change. The reads
then carry only a *component*: gbcms correctly counts the annotated MNP as
~absent and puts the carriers in `partial_alt`. Rescue re-counts each changed
position as an SNV and, when the annotated MNP is absent, reports the
best-supported component under the MNP's row, flagged and audited. It
relabels; it never invents reads.

### Where we are (PR #101, open against develop)

| Area | State |
|---|---|
| Default output (rescue off) | Unchanged — 35/35 recorded real-data runs identical to develop |
| Bug fixes | Done (list below) |
| Rescue gate | Works on IMPACT; **too generous on ACCESS** (below) — item 12a |
| Merge of rescued duplex/simplex | Can add two different alleles — item 12b |
| Review nits | Contig label, test invariants, indel test — items 12c–12d |

**Done and staying, whatever the gate decision** (history: items 1–11 below):
- Rescued rows are one coherent genotype: all count, fragment, strand,
  strand-bias and mFSD columns come from the adopted component's own counting
  pass (previously only `alt_count` changed).
- Audit leak across samples fixed (reset per sample).
- VCF `GR` no longer split by parsers (positions joined with `+`).
- Grouped MNPs skipped, so rescue cannot undo T1 exclusive assignment.
- Failed positions read `ref_fail`, not a silent 0.
- Every rescued row flagged `RESCUED_COMPONENT(...)`, a WARNING per rescue and
  one when the flag is enabled; MAF and VCF carry identical content.
- Engine counts `mnp_confirmed_alt`: MNP ALT reads in which every changed base
  was read (none low-quality, none N) — reads that *show* the whole MNP.

### What real data showed (2026-09-23)

| Data | Result |
|---|---|
| IMPACT, 2 cohorts, 68 samples / 140 MNP rows | 13 rescues, all exact vs sign-out, all somatic in the matched normal, all with **0** reads showing the whole MNP |
| Same, partial-read make-up | 97–100% genuine single-change carriers |
| `--apply-baq` on all 7 candidate rows | No outcome, count or confirmed-read change |
| ACCESS, 20 samples / 22 MNP rows (duplex + simplex) | **2 wrong rescues**, each also a duplex/simplex conflict: BRCA2 AAG>TAC duplex adopted a **germline** SNP (normal VAF 0.43; 6 reads showed the whole MNP but the allowance was 8) → combined fragment ALT 25 → 724; KRAS ACC>CCA simplex rescued with 1 whole-MNP read (allowance rounded up to 1) although duplex shows 4 |

Cause: the error allowance `ceil(partial × 10^(−q/10))` is loose twice over —
it counts any error rather than an error to the one ALT base (~1/3), and
rounding up always excuses one read.

### Remaining work — item 12 "simplify and close"

**Decision (operator, 2026-09-23): path A** — simplify the gate, steps
12a–12g. (Path B, returning to `ad == 0`, was declined; logged in
`REJECTED.md`.)

**12a. One rule instead of an allowance.** Rescue only when *no* read shows
the whole MNP (`mnp_confirmed_alt == 0`) and partial reads outnumber full-ALT
reads. Delete `_confirmed_error_allowance` and its test. Why: every correct
rescue seen has exactly 0; every wrong one had ≥1. Sequencing error making a
fake whole-MNP read needs a specific substitution at high quality at another
position — expected well under one read at real depths. Expected on the data
above: the 13 IMPACT rescues unchanged; GRIN2A, TP53, BRCA2 and KRAS all kept
as `haplotype_confirmed`; both ACCESS conflicts disappear.

**12b. Merge warns on mixed rescue.** `gbcms merge` (with or without
`--add-combined` — review finding: the two flavors' columns describe different
alleles in one row either way; the combined-column note is added only when
those columns are written): when the
duplex and simplex rows disagree (one rescued and one not, or rescued to
different components), log a WARNING per row and a per-run count naming the
rows. Counts unchanged in this PR. Follow-up ticket: leave the combined
columns empty (NA) for such rows and write a conflicts file beside the merged
MAF.

**12c. Labels use the output file's contig naming** (operator decision;
dropping the contig was declined). `RESCUED_COMPONENT(...)` and audit labels
take the contig exactly as the row writes it: the input MAF's `Chromosome`
for MAF input (today they show the stripped internal name — row `chr1`,
label `1`); for VCF input, whatever the output writes, which T8 makes the
input's own naming. Never surfaced on b37 data, which has no `chr`.

**12d. Tests.** Four counting invariants in the two engine counter tests; an
end-to-end test that an indel row with `partial_alt > ad` is untouched by
`--rescue-mnp`; a merge test for the 12b warning; update pins for 12a/12c.

**12e. Docs.** Gate rule, merge warning, labels in architecture /
output-formats / dna.md / skill / CHANGELOG; BAQ caveat ("near indels BAQ can
lower qualities so reads stop counting as showing the whole MNP").

**12f. Durable harnesses.** The validation harnesses lived in a session
scratchpad and were lost on a restart. Recreated in a local-only directory
outside the repo (they carry patient-data paths), with their own branch and
develop venvs: flag-off parity on a seeded real-data panel (IMPACT tumours +
ACCESS duplex/simplex, all signed-out variants — the earlier recorded runs were
lost too), the IMPACT cohort with matched normals (seeds 1 and 2 reproduce the
two T7 cohorts), per-read traces, and the ACCESS merge study (also checks the
merge warning fires exactly on disagreeing rows). The partial-read make-up and
BAQ studies were one-off investigations; their results are recorded above.

**12g. Validation (one BAM at a time).** Flag-off parity (35 runs, must stay
identical); both IMPACT cohorts (13 rescues must be unchanged); ACCESS merge
study (expect 0 conflicts and no germline adoption).

**Item 12 result (2026-09-23, met).** Code: red battery 69e3722a → fix
b422a9ef → docs be566216. Real data (harness `~/test/gbcms/harness/t7`, one BAM
at a time): flag-off parity 16/16 runs identical (8 IMPACT tumours + 4 ACCESS
duplex/simplex pairs, all signed-out variants); IMPACT cohort seed 1 — 5
rescued (all exact vs sign-out, somatic in normal, 0 whole-MNP reads), GRIN2A
and TP53 kept (33 / 253 whole-MNP reads); seed 2 — 8 rescued, all exact,
somatic, 0 whole-MNP reads; ACCESS 20 samples — 0 duplex/simplex conflicts, 0
merge warnings, BRCA2 and KRAS kept in both flavors (combined fragment ALT
unchanged: 25 and 4); traces agree.

**Out of scope (tickets):** merge NA + conflicts file; a run-start summary of
enabled options and their implications; fillout samples where the MNP is
absent but a germline component is present (upstream annotation). The
always-on `mnp_confirmed_alt` counter stays (one check per MNP ALT read) — an
accepted deviation from "engine should be output-aware".

### History (items 1–11, all merged on the branch)
1. Gate `partial_alt > ad` instead of `ad == 0` (a single masked read blocked TERT).
2. Replace only if the component beats the MNP; else `no_improvement` (reads with an indel inside the block).
3. Adopt the component's whole `BaseCounts`.
4. Audit keeps the MNP's own counts; `outcome`; `ref_fail`.
5. Diagnostics recomputed for rescued rows.
6. Grouped rows skipped.
7. Audit reset per sample.
8. mFSD Parquet provenance logged/documented.
9. Shared `_engine_kwargs`, per-row `_diagnostic_flags`, dead code removed.
10. `mnp_confirmed_alt` and the confirmed-haplotype gate (fixed GRIN2A germline adoption on IMPACT).
11. `RESCUED_COMPONENT` flag, warnings, `+`-joined positions (VCF parse fix).

**Tests.** `tests/test_rescue_mnp.py` (end-to-end MAF and VCF battery,
invariants on written rows, audit/outcome units),
`tests/test_mnp_concordance.py` (`TestONPCarrierShapes`, confirmed-read engine
pins). Harnesses (local, patient data): flag-off parity, IMPACT cohort with
matched normals, per-read traces, partial-read make-up, BAQ effect, ACCESS
duplex/simplex merge.

## T8 — Output keeps the input's contig naming (#103; own branch — changes default output)

**Problem (reproduced 2026-09-23, synthetic).** Readers strip `chr`
(`CoordinateKernel.normalize_chromosome`) and `Variant` keeps no original
name. VCF input → MAF writes `Chromosome=1` for `chr1`; VCF output writes
records as `1` under a `##contig=<ID=chr1>` header taken from the reference
`.fai` — a malformed VCF (htslib: "Contig '1' is not defined in the
header"). MAF input keeps its `Chromosome` column (row passthrough), but labels
built from the internal name (the rescue audit) show `1`. Counting is correct:
contig reconciliation between variant file, FASTA and BAM works. Unseen so far
because MSK b37 (`1`…`22`) and Ensembl-named hg38 have no `chr`.

**Design.** Keep the input's contig name on `Variant` at read time (VCF and
MAF readers) and write it everywhere a contig is written: MAF `Chromosome` and
`vcf_region` for VCF input, VCF `CHROM`, and derived labels (rescue flag and
audit — T7 12c). VCF `##contig` lines must declare the names the records use.
Log once per run (INFO) when the input's naming differs from the reference or
BAM naming, naming both, so the reconciliation is visible.

**As implemented.** `Variant.original_chrom` (set by both kernel readers) with an
`output_chrom` accessor; `chrom` stays normalized, so counting is untouched. Writers,
rescue labels (the old `_output_contig` helper is gone) and the mFSD Parquet use
`output_chrom`. `##contig` lines: each `.fai` contig under the input's name(s) for it
(reference length kept); input contigs absent from the `.fai` are declared without a
length. The `.fai` is read once per run (was once per sample). One INFO line per
distinct naming pair (reference, then any BAM). Observations Parquet keeps the
internal name — it echoes loci from inside the engine (documented). Names are
paired by `CoordinateKernel.contig_key`, the Python mirror of the engine's
`normalize_contig` (`chrM` ~ `MT`, not just a `chr` strip); a `.fai` listing one
contig under two aliases declares the input's name once. `gbcms merge` joins on
the same key (review finding: inputs counted from differently named variant files
otherwise split into two rows), keeps the first input's name, and logs differences.

Second review round:
- The FASTA fetch folded no mitochondrial alias, so `chrM` vs an `MT` FASTA was
  `FETCH_FAILED`. This predates T8, but T8's log claimed reconciliation, and the
  first chrM test asserted names, not counts. `fetch_region` now also tries
  `MITO_SPELLINGS`, and the test matrix asserts counts.
- Merge writes each contig one way. Rows only a later input has had kept its
  spelling.
- Merge warns when one input names a contig two ways, and rejects helper-named
  input columns.
- Merge's row order is deterministic (found by the merge parity harness). Rows
  follow the inputs via per-input row numbers, which works on every polars 1.x.

**Tests.** `chr`-named input × {VCF, MAF} × {VCF, MAF output}: output names
equal the input's; VCF output parses with no undefined-contig warning; b37
naming unchanged. **Acceptance:** flag-off parity on the real-data sets stays
byte-identical (all b37, no `chr`).

## T9 — Span-aware exon-edge BAQ rule (#106; count-affecting, RNA + GTF only)

**Finding** (adversarial review of T6; pre-existing, unchanged by T6). The rule
keys on `variant.pos`:
- `baq_applies` takes `exon_boundary_dist`, the distance from the variant's
  *first* base to the nearest annotated boundary.
- Heuristic BAQ penalizes read bases within 5bp of a CIGAR N. For an exon
  ending at `E`, bases `E-5..E-1` of every read spliced at `E` lose 20 BQ.
- So a multi-base variant whose `pos` is more than 5bp from `E`, but whose REF
  span reaches within 5bp (an MNP at `E-7..E-4`), keeps BAQ. Its edge-side
  bases are masked in every spliced read, and MNP confirmation collapses
  (`MnpResult` confirmed=false). All three views agree on the rule, and it is
  the rule that is wrong.
- With `--rescue-mnp`, the MNP row (BAQ on) and its component recount at `E-4`
  (BAQ off) are counted under two BAQ regimes.
- Left edges are safe: `pos` is the leftmost base.

**Options.**
- **A (recommended). Span-aware variant rule.** The distance is from the
  nearest annotated boundary to the REF span `[pos, pos+len(REF))`, and 0 when
  a boundary falls inside it. `baq_applies` uses it, so there is still one
  resolved value for main, per-transcript, ASJD and each rescue recount. SNVs
  are unchanged.
- **B. Read-level.** Skip the N-penalty in `apply_heuristic_baq` for Ns whose
  ends are annotated junctions (±`JUNCTION_TOLERANCE`), and keep indel
  penalties. Rejected unless A's measurement leaves residual edge masking:
  - It re-penalizes *novel* junctions near a variant, which is ASJD's
    splice-creation signal.
  - It restores indel penalties at edges.
  - It is a larger change to the main counts.

**Decision for the operator.** Does `exon_boundary_dist` become the span
distance too?
- **Recommended after measuring: no.** The column stays `pos`-based and the
  rule uses the span internally. The column change would touch about 10% of
  signed-out deletions near exon edges, and deletions gain nothing from the
  rule. It would be a user-visible change with no counting benefit.
- The first recommendation was "yes", made before the measurement below.

**Measure first (as T6).**
1. **Demand:** signed-out multi-base variants (MNP, DNP, deletion, delins) in
   the local RNA truth cohort whose span reaches within 5bp of an annotated edge
   while `pos` does not.
2. **Effect:** on the T6 FORTE harness, add MNP and deletion probes with the
   span end 0–4bp inside right exon edges and `pos` 6–10bp away, plus controls.
   Report main, per-transcript and ASJD deltas and `mnp_confirmed_alt`, with
   and without the change.

**Red battery.**
- An MNP whose span reaches E1's donor edge, with `pos` more than 5bp away. Its
  spliced carriers count as confirmed ALT in the main counts, per-transcript
  equals main, and ASJD sees them.
- With `--rescue-mnp`, the MNP row and its component recount resolve the rule
  identically.
- Guards: SNV rows unchanged; left-edge geometries unchanged; BAQ still applies
  away from edges.

**Fix.** Add a span distance to `AnnotationIndex` (the binary search to the
range, as `intron_boundary_in_range` already does). `count_variant_from_cache`
resolves it once, and the per-transcript and ASJD passes keep reusing it.

**Acceptance.**
- The FORTE harness: SNV rows byte-identical to develop, multi-base edge probes
  change as designed, and repeat runs byte-identical.
- The RNA truth cohort: any affected sign-out rows are adjudicated read-by-read
  (the BAM is truth).
- RNA-only, so exempt from the legacy parity oracle. Confirm `PARITY_FIELDS` is
  untouched.

**Measured (2026-09-24, local data; no data in repo).**
- **Demand.** Across 74,617 unique signed-out multi-base events:
  - 6,130 deletion loci (10% of 61,123) have the geometry, 4,887 of them
    crossing an exon boundary;
  - 99 MNP loci (0.7%) have it.
- **Effect.** An A/B of develop vs a span-aware prototype on 3 FORTE RNA
  samples, 326 probes (231 in the T9 class). The harness was clean: controls
  byte-identical, the column unchanged, invariants held, repeat runs
  byte-identical.
  - **Deletions: no effect.** 0 rows moved across 7 categories and 510
    probe-sample rows (crossing at 12, 30 and 60bp, and ending short), apart
    from one fragment count. Deletion calls are structural.
  - **MNPs, ending 4–5bp short of a right edge:** RD +0.03% to +0.17%, with
    per-transcript and ASJD moving alike. AD and DP are unchanged.
  - **MNP confirmation.** Counted exactly with pysam: a median 90–92% of the
    reads covering such an MNP have a splice N within BAQ's radius of an MNP
    base. REF and ALT calls vote on the unmasked bases, which is why the
    counts barely move. But `MnpResult::Alt` confirmation needs no masked
    base. So under the current rule about 9 in 10 real carriers there can never
    be MNP-confirmed, and `--rescue-mnp`'s `haplotype_confirmed` safeguard is
    largely blind at right exon edges in RNA.
- **Conclusion.** The fix matters for MNP confirmation, which only the opt-in
  RNA rescue path consumes. It is not needed for counts, and not needed for
  deletions at all.
  - Priority: low.
  - The red battery should pin `mnp_confirmed_alt` for an edge MNP's spliced
    carriers.

## T10 — Decomposed twin inherits what the original resolved (#107; count-affecting)

**Finding** (adversarial review of T6; pre-existing).
- The homopolymer-decomposed twin (normalize/engine.rs, Step 5) is built with
  `gene_strand: None`.
- `count_bam_binned_core` resolves gene strand from the GTF only for
  `variants`, and the `--enforce-strandedness` WARN counts only `variants`.
- So with RNA + GTF + `--enforce-strandedness` at a decomposition locus, the
  twin admits antisense reads (`is_sense_strand(.., None, ..)` is true). This
  inflates its AD and biases the dual-count arbitration
  (`counts_decomp.ad > counts_orig.ad`) toward it.
- The reported row then carries antisense evidence that per-transcript and
  ASJD, which classify the stranded original, exclude.
- No warning covers the twin: a silent no-op of an RNA default.

**Audit in the same ticket.** The twin also gets `repeat_span: 0` ("repeat
info is not critical"), yet by construction it sits in a homopolymer.
`repeat_span` sets the fetch-window pad (`max(5, repeat_span + 2)`) and the
wrong-length rule's repeat branch. List every field the original carries after
prep and resolution. The twin inherits each unless a documented reason says
otherwise.

**Red battery.**
- RNA + GTF + `--enforce-strandedness`, at a homopolymer-decomposable deletion
  (the `check_homopolymer_decomp` geometry).
  - Antisense reads support the decomposed form, and sense reads support the
    original.
  - The antisense reads are not counted.
  - The arbitration is decided by sense reads.
  - The per-transcript counts are consistent with the row.
- A `repeat_span` case, if the audit finds a behavioral difference.

**Fix.**
- `repeat_span`: set it where the twin is built, in `prepare_variants`, to the
  original's computed value. Today it is computed just *after* the twin
  (Step 5 builds the twin; the `repeat_span` computation follows), so compute
  it first. The legacy oracle then dual-counts the same twin, so binned↔legacy
  parity holds (invariant 1).
- Strand: after strand resolution in `count_bam_binned_core`, the twin takes
  its original's resolved strand (same contig and position). The legacy path
  has no strand.

**Acceptance.**
- A strand-only change leaves DNA byte-identical: re-run the T8 b37 parity
  harness.
- If `repeat_span` inheritance changes DNA counts, it is count-affecting for
  DNA:
  - measure it on the T1/ACCESS truth harnesses;
  - keep the legacy oracle's dual-count on the *same* twin, so binned↔legacy
    parity holds (invariant 1).
- RNA: rows at decomposition loci change only by the excluded antisense reads.
  Measure how often decomposition fires and wins on the FORTE cohort.

**Measured (2026-09-24, local data; no data in repo).** Develop vs a prototype
twin that inherits strand and `repeat_span`, with log-only instrumentation
reporting both alleles' counts at every dual-count.
- **Frequency.** A decomposed twin exists at only 13 of 61,123 signed-out
  deletion loci (0.02%), all in annotated exons.
- **DNA,** 11 of those loci on their own samples' BAMs. The other 2 samples
  are absent from the local key files.
  - **The twin wins at 2 of 11.** At one locus the original's ALT count is 5
    and the twin's is 78, so decomposition does real work.
  - **`repeat_span: 0` undercounts the twin's ALT** at 7 of 11 loci, by up to
    15 reads (about 7%). The truncated scan window misses carriers, so the
    arbitration compares a handicapped twin.
  - One output row changed; there were no winner flips at these loci.
  - The original's counts are identical between builds at 11 of 11, so the
    difference is `repeat_span` alone (DNA has no strand).
- **RNA,** 40 homopolymer-decomposable probes × 3 FORTE samples, strandedness
  on and off.
  - The missing strand leaks antisense reads into the twin at 14 of 120 loci:
    52 reads, 0.065% of the twin's RD.
  - No output change and no winner flip.
  - `repeat_span` had no effect on these REF-only loci, as expected: it acts on
    windowed ALT evidence.
- **Conclusion.** `repeat_span` inheritance is the count-affecting part
  (DNA, at real twin loci). The strand part is a consistency fix with
  negligible measured effect. Priority: T10 before T9.

**As implemented (narrowed after review and a read census).**
- The review showed the measured `repeat_span` gain was confounded with the
  twin's shape. The read census then showed the `repeat_span` change has no
  benefit, so it was dropped. The twin stays at 0 (unique sequence), with the
  measured reason in the code. Census at the 11 real twin loci:
  - the called delins is exact at 4 loci, yet the twin claims 93–97% of those
    same reads;
  - the twin as built (`C^(L-1)X`) is exact at 1 locus;
  - the documented D(1)+SNV form is exact at 1 locus (the 5 → 78 win);
  - at 4 loci most reads carry other alleles of the run (same-length MNPs and
    rewrites) that both classifiers claim.
  - The `repeat_span` change narrowed the original's margin where the original
    is right (by 2–15 reads), and flipped no outcome.
- Shipped:
  - strand inheritance, outside the GTF block, so a strand supplied upstream
    reaches the twin;
  - `WARN_HOMOPOLYMER_DECOMP` re-derived per sample (the prepared variants are
    shared by all samples of a run; B's row carried A's flag);
  - docs and comments that describe the corrected allele as built.
- Red battery (`tests/test_decomposed_twin_contract.py`): strand, the
  per-sample flag, and a parity guard with invariants.

## T11 — Homopolymer decomposition arbitration (#111; redesign; measure first)

**Finding** (T10 review and read census; pre-existing). The dual-count pits
two permissive classifiers against each other. Reads at real loci often carry
a third allele that neither describes, and "more ALT wins" hangs on thin
margins (the twin claims 93–97% of the called allele's own carriers). The
twin's allele (`REF[..len-1] + X`) does not match its documentation
(`C^(L-2)X`).

**Direction.** Arbitrate on exact haplotype support, the census method.
- Reduce each read to its sequence across the run plus flanks (from the
  alignment) and count the reads that carry each candidate exactly: the called
  allele, the corrected alleles (both shapes), and "other".
- Report the allele the reads carry, and say so when it is neither form.
- Revisit `repeat_span` for any corrected allele kept.

**Measure first.**
- The census harness at the 11 real twin loci is the truth set.
- Twins are rare: 13 of 61,123 signed-out deletion loci. Priority: after the
  6.5.0 cut unless the operator decides otherwise.
- Related, filed separately:
  - #112: when the corrected allele wins, per-transcript/ASJD, merge,
    observations and input validation don't follow the winner;
  - #110: a VCF delins with a one-base ALT is treated as a pure deletion.

## Order & discipline

T1 (own branch, own review) → T2 (own branch; needs the FORTE geometries) →
T3+T4+T5 (one observability/cleanup branch) → T6 (measurement decides).
T8 (contig naming, #103) follows T7 on its own branch.
T7 is independent (Python-only, opt-in path) — own small branch
(`feature/mnp-rescue-gate`), can land anytime; its leak fix may go first.
Each branch: battery red → implement → suites + clippy + lint gate →
sonnet adversarial review → real-data acceptance named above → merge to
develop. 6.5.0 cut only after T1's ACCESS rerun and T2's cohort recheck.
T10 then T9, each on its own branch; both were measured on 2026-09-24.
- **T10 before the 6.5.0 cut.** It is small, and it is count-affecting for DNA
  at real twin loci: `repeat_span` biases the dual-count arbitration.
- **T9 is low priority.** Deletions are unaffected, and MNP counts move at most
  0.17%. It restores MNP confirmation at right exon edges, which only the
  opt-in RNA rescue safeguard uses. It can follow the cut. The measurement
  recommends keeping `exon_boundary_dist` `pos`-based.
- **T11 is a redesign.** It measures first on the T10 census harness and
  follows the cut. Twins are rare (0.02% of deletion loci).
