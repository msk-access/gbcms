# gbcms — Code Review Remediation Plan

**Version reviewed:** 5.3.0 (`develop`) · **Date:** 2026-06-26 · **Reviewer:** Principal engineering deep-dive (Phases 1–5)

This plan enumerates **every** finding from the architecture/algorithm/statistics/concurrency review — Critical, High, Medium, and Low — as actionable tickets. Each ticket carries: the detailed issue, exact code locations, **git archaeology** (why the code is the way it is, where recoverable), root cause, a step-by-step fix, a test plan, and regression risk.

## How to read this

Each ticket has a stable ID (`CR-*`, `HI-*`, `ME-*`, `LO-*`). IDs are referenced in the milestone plan at the end. "Git archaeology" deconvolutes the originating commit and the author's apparent intent, so we don't re-break a deliberate fix while patching its side effect.

### Severity legend
- **CRITICAL** — silently produces wrong allele counts or clinical labels, or panics on valid-ish input. Fix before next release.
- **HIGH** — wrong results in a definable regime (RNA mode, low-BQ reads, deep sites), or a feature that is silently inert.
- **MEDIUM** — correctness-adjacent (output parity, multiplicity, over-conservative stats) or meaningful performance tax.
- **LOW** — hygiene, docs, micro-bias, dead code, defensive hardening.

### Master index

| ID | Sev | Area | One-line |
|----|-----|------|----------|
| CR-1 | CRIT | Binning | Bin anchored by >10 kb deletion under-fetches its tail → AD/ADF undercount |
| CR-2 | CRIT | Alignment | WFA fast-path classifies with no base-quality gate; counts depend on backend |
| CR-3 | CRIT | variant_checks | `len()-1` underflow + slice panic on empty ALT/REF |
| CR-4 | CRIT | variant_checks | Tolerant large-deletion match accepts ALT with no sequence check |
| CR-5 | CRIT | mFSD stats | LLR ±∞ from one tail fragment; KS p-value invalid at n=5–20 (drives clinical labels) |
| HI-1 | HIGH | Pipeline | Sample failures swallowed; exit code 0 even when all samples fail |
| HI-2 | HIGH | Fragment | Read-level DP/AD/RD double-count supplementary/secondary when their filters are off |
| HI-3 | HIGH | Alignment | `OFF_TARGET_THRESHOLD=20` not read-length-aware; drops deletion-ALT reads |
| HI-4 | HIGH | PairHMM | `ln(0)` at Q0 → NaN LLR; `p_correct` not clamped |
| HI-5 | HIGH | PairHMM | N-base emission not truly LLR-neutral near indels |
| HI-6 | HIGH | variant_checks | `check_complex` insertion off-by-one at exclusive REF end |
| HI-7 | HIGH | RNA | `gene_strand` never populated → sense/antisense meaningless, `enforce_strandedness` a no-op |
| HI-8 | HIGH | RNA / mFSD | Spliced-read fragment size inflated by full intron length (N not discounted) |
| HI-9 | HIGH | RNA | Splice-motif classifier ignores strand → all minus-strand canonical junctions mislabeled |
| HI-10 | HIGH | mFSD stats | LLR model hardcoded/uncited; report mislabels it as the PairHMM LLR |
| HI-11 | HIGH | mFSD stats | 6 KS p-values per variant get no FDR correction at any scope |
| PF-1 | HIGH | Perf | mFSD stats + size arrays computed for every variant on every run (no `mfsd` flag in Rust) |
| PF-3 | HIGH | Perf/IO | htslib decode threads never enabled (`set_threads` absent) |
| ME-1 | MED | Output | `sub_nuc`/`mono_nuc` fields computed but dropped from VCF |
| ME-2 | MED | Output | Transcript-count delimiter `;` vs documented `|`; MAF writes raw `;` |
| ME-3 | MED | variant_checks | Soft-clip boundary differs between `check_complex` Phase 1 and `extract_raw_read_window` |
| ME-4 | MED | Alignment | Haplotype parity collisions degrade locus to universal ties |
| ME-5 | MED | PairHMM | `prob_emit_y = 0.999` unjustified/asymmetric vs quality-derived emit_x |
| ME-6 | MED | RNA/GTF | Contig normalization is `chr`-strip only → `MT`/`chrM` never reconcile |
| ME-7 | MED | RNA | `STRAND_DISCORDANT` uses raw genomic strand without R1/R2 folding → false flags |
| ME-8 | MED | RNA | ASJD BH includes `p=1.0` padding in `n` → over-conservative FDR |
| ME-9 | MED | mFSD | `MIN_FOR_KS=5` vs report `min_alt=3` → NaN-KS variants get classified |
| ME-10 | MED | mFSD | `sub_nuc_enrichment` NaN is ambiguous (empty ALT vs empty REF) |
| ME-11 | MED | Report | Formatters guard `isnan` but not `isinf` → `+Inf` LLR renders literal "inf" |
| ME-12 | MED | Fragment | "Majority rule" is structural-ALT-wins consensus; single spurious indel CIGAR overrides high-BQ REF |
| PF-2 | MED | Perf | One rayon task per bin, no depth cap → ultra-deep bin serializes the tail |
| PF-4 | MED | Perf | rayon pool rebuilt per `count_bam_binned` call (×2 with MNP rescue) |
| ME-13 | MED | Perf/IO | Adjacent bins' padded fetch ranges overlap → reads fetched/decoded twice |
| LO-1 | LOW | Build | Duplicated `.pyi` stubs; `count_bam` stub missing params |
| LO-2 | LOW | API | Legacy `count_bam` (one fetch per variant) still exported |
| LO-3 | LOW | Binning | Bin span can chain well past 10 kb |
| LO-4 | LOW | Fragment | u64 QNAME hash stores no key (collision prob negligible) |
| LO-5 | LOW | Docs | "orientation-aware hashing" is a misnomer |
| LO-6 | LOW | Fragment | FSB strand taken from best-quality read, not deterministic |
| LO-7 | LOW | Util | `median_qual` upper-median bias on even counts |
| LO-8 | LOW | variant_checks | Deletion left-anchor invariant undocumented/undefended |
| LO-9 | LOW | GTF | Empty/exon-less GTF degrades silently; weak diagnostics |
| LO-10 | LOW | GTF | `.`/unstranded exon defaults to `+` silently |
| LO-11 | LOW | RNA | Editing-site flag is positional-only; doesn't verify A>I base change |
| LO-12 | LOW | RNA | NH-tag rescue only matches `U8`/`I32` Aux variants |
| LO-13 | LOW | Report | KDE bandwidth clamp (5.0 bp) arbitrary/undocumented |
| LO-14 | LOW | Perf | `threads` not clamped to `available_parallelism()` |
| LO-15 | LOW | Alignment | `usable_count >= 3` magic number duplicated in 3 places |
| DX-1 | LOW | Maintainability | Code comments & log strings cite ephemeral ticket labels (`CR-1`/`HI-11`/`ME-8`/`P4c`); strip them and state the reason instead |

---

# CRITICAL

## CR-1 — Bin anchored by a large deletion/DelIns under-fetches its right tail

**Locations:** `rust/src/counting/engine.rs:156` (seed), `:182-183` (extend), consumer window `:1060-1062`.

**Issue.** `build_genomic_bins` seeds the fetch end from the *anchor* (leftmost) variant as `bin_end = bin_start + window` (window = 10 000). The extension `bin_end = bin_end.max(var_end + window/2)` that accounts for a variant's true reference span runs **only inside the inner `j` loop, i.e. only for subsequent variants** (`j > i`). It is never applied to the anchor variant itself. Meanwhile `count_variant_from_cache` computes each variant's overlap window as `v_end = variant.pos + ref_allele.len() + window_pad` (`:1062`). So for a bin whose anchor is a deletion/DelIns with `ref_allele.len() > BIN_WINDOW` (or even `> window − padding`), reads aligning in `[bin_start + 10 000, anchor.pos + ref_len]` are **never fetched into `read_cache`** and therefore cannot be classified. For a large deletion these are exactly the reads spanning the right breakpoint and carrying the D-op ALT evidence → **AD/ADF undercount**, and divergence from the legacy `count_single_variant` path (which fetches per-variant `[pos − pad, pos + ref_len + pad]` at `:1570` and has no gap).

**Git archaeology.** Introduced in `761b98b` (2026-03-18, *"feat: nextflow DNA/RNA module split, normalize module"*) — a large multi-purpose refactor that also introduced WFA and `gene_strand`. The `var_end + window/2` line was clearly written to cover *trailing* variants added to a bin; the author did not revisit the anchor's simple `bin_start + window` seed. This is an oversight inside a big commit, **not** a deliberate trade-off — there is no comment justifying a 10 kb cap on the anchor's span.

**Root cause.** The fetch interval does not enforce the invariant it must: `bin.end ≥ maxᵥ(v.pos + v.ref_allele.len() + window_pad_v)` over *all* variants in the bin, including the anchor.

**Fix.**
1. Seed the end with the anchor's true span:
   ```rust
   let anchor = &variants[first_idx];
   let mut bin_end = (bin_start + window)
       .max(anchor.pos + anchor.ref_allele.len() as i64 + window / 2);
   ```
2. Keep the per-`j` extension as-is (it already covers trailing variants).
3. Add a debug assert after bin construction: for every `idx` in the bin, `bin.end >= variants[idx].pos + variants[idx].ref_allele.len() + max(5, repeat_span+2)` and `bin.start <= variants[idx].pos - max(5, repeat_span+2)`.

**Tests.** New parity test: a synthetic BAM with reads spanning a 12 kb deletion's right breakpoint; assert `count_bam_binned` AD == `count_bam` (legacy) AD. Add a unit test on `build_genomic_bins` asserting the end-coverage invariant for a single >10 kb variant and for a 200-variant cluster.

**Regression risk.** Low. Strictly widens fetch for affected bins; no change for bins whose anchor span < window. Slightly larger `read_cache` for those bins.

**Effort.** S (≈ ½ day incl. tests).

---

## CR-2 — WFA fast-path classifies with no base-quality gate (backend-dependent counts)

**Locations:** `rust/src/counting/wfa_router.rs:54` (`wfa_fast_path` signature — no quals param), decision at `:93-109`; entry `variant_checks.rs:61,89`; contrast PairHMM gate `pairhmm.rs:454-458` and SW masking `alignment.rs:233-245`.

**Issue.** The PairHMM and SW paths are base-quality-aware: SW masks sub-`min_baseq` bases to `N` (`alignment.rs:233`), and both skip reads with `usable_count < 3` high-quality bases. The WFA fast-path aligns the **raw `read_seq`** with no quality input at all — `read_quals`/`min_baseq` are never passed in — and `med_qual` is only *attached* to the result, never used to gate. Consequently a read whose discriminating base is Q2 but happens to match ALT at edit distance 0 is returned as **confident ALT** at `wfa_router.rs:98` and short-circuits out of `pangenomic_classify` (`variant_checks.rs:89`) *before* PairHMM runs. The same read, routed through PairHMM, would be quality-gated to "neither." **Enabling the fast path therefore changes counts** — and it admits low-quality ALT support precisely in the low-VAF cfDNA regime the BQ model exists to protect.

A closely related defect in the same function (tracked together): the off-target branch at `:99` uses **global** edit distance (`align_end2end`) while PairHMM uses semiglobal with free end gaps. For a deletion ALT (haplotype much shorter than the read), the global distance to *both* classes can exceed `OFF_TARGET_THRESHOLD` purely from length, so true deletion-ALT reads are dropped as NEITHER. (Threshold value itself is HI-3.)

**Git archaeology.** WFA introduced in `761b98b` (2026-03-18); CHANGELOG (≈ line 544) markets it as *"Phase 3 WFA+PairHMM unification … significantly improves classification accuracy on complex multi-allelic variants."* The unification focused on *alignment topology*, not on porting the quality model into the fast path — the BQ gate lives in the older PairHMM/SW code and was never threaded into the new fast path. Oversight, not intent.

**Root cause.** The fast path and the fallback do not share a quality contract; the fast path can make *definitive* calls without the gate that the fallback enforces.

**Fix.**
1. Add `read_quals: &[u8]` and `min_baseq: u8` to `wfa_fast_path`.
2. Before any definitive REF/ALT return (`:95,:98`), require the discriminating position(s) to pass `min_baseq` and the read to satisfy the same `usable_count >= MIN_USABLE_BASES` gate (see LO-15) used by PairHMM. If it fails, return `None` so the fallback (or "neither") decides.
3. For the off-target branch, only short-circuit when `|read_len − hap_len|` is small; otherwise return `None`. Or normalize the distance by the length delta before comparing to the threshold.
4. Thread `read_quals`/`min_baseq` from `pangenomic_classify` (already has the read).

**Tests.** (a) A read with a Q2 ALT-discriminating base must classify identically (`neither`) with WFA on and off. (b) A deletion-ALT read whose ALT haplotype is 30 bp shorter than the read must be classified ALT (not off-target). (c) Backend-parity test: random reads classified by WFA-on vs WFA-off must agree for all reads passing the BQ gate.

**Regression risk.** Medium — will *reduce* some ALT counts that were previously (incorrectly) admitted. Communicate as a correctness fix; expect small VAF changes at low-BQ sites. Gate behind a release note.

**Effort.** M (1–2 days incl. parity harness).

---

## CR-3 — `len()-1` underflow + slice panic on empty ALT/REF

**Locations:** `rust/src/counting/variant_checks.rs:1052-1053` (insertion), `:1428-1430` (deletion).

**Issue.** `check_insertion` computes `expected_ins_len = variant.alt_allele.len() - 1` and slices `&variant.alt_allele.as_bytes()[1..]`. `check_deletion` mirrors this with `ref_allele`. If a variant reaches these with an empty ALT (insertion) or empty REF (deletion) — a malformed VCF record, or a non-left-anchored representation — `len() - 1` underflows `usize` to `usize::MAX` and the `[1..]` slice panics ("slice start out of range"). A panic here becomes an opaque `PyRuntimeError` with no locus context (caught at `engine.rs:728`), so a single bad record kills the sample with an uninformative message.

**Git archaeology.** This is the standard VCF left-anchor assumption (REF/ALT share a leading anchor base). The code predates the binning refactor and assumes upstream normalization guarantees the anchor. `normalize.py`/`prepare_variants` usually does, but there is no defensive check at the Rust boundary, and CLI users can pass pre-normalized or hand-built variant lists.

**Root cause.** The VCF left-anchor invariant is assumed, never asserted, at the function entry.

**Fix.** Guard at entry of each function:
```rust
if variant.alt_allele.is_empty() { return ClassifyResult::neither(...); }   // check_insertion
if variant.ref_allele.is_empty() { return ClassifyResult::neither(...); }   // check_deletion
```
Prefer surfacing a one-time `warn!` with `chrom:pos` so malformed input is visible rather than silently "neither." Optionally validate in `prepare_variants` and reject with a clear Python-side error listing the offending records.

**Tests.** Unit tests passing a `Variant` with empty ALT / empty REF to each function; assert `neither` and no panic. Integration test feeding a malformed MAF row; assert graceful skip + warning, sample still completes.

**Regression risk.** None (pure hardening).

**Effort.** XS (≈ 1–2 h).

---

## CR-4 — Tolerant large-deletion match accepts ALT with zero sequence verification

**Locations:** `rust/src/counting/variant_checks.rs:1598-1599` (windowed/tolerant path), `:1488-1507` (strict path analogue).

**Issue.** When an observed deletion's length differs from the expected (`del_len_usize != expected_del_len`) but a ≥50 bp / ≥50 %-reciprocal-overlap rule passes, the code sets `del_ok = true // tolerant match — length differs, skip seq check`, **bypassing the deleted-sequence (S3) verification entirely**. So a large deletion of a *different* sequence at a nearby position, as long as its length overlaps the expected by ≥50 %, is accepted as ALT evidence. In segmental-duplication / repeat-dense regions this can absorb an unrelated structural variant as ALT for the target deletion → **false-positive AD**.

**Git archaeology.** This is a **deliberate** trade-off, and the history explains why. The v4.x line *"Interior REF guard removed: Eliminated the `has_large_cigar_del` guard that massively overcounted REF for large deletions by misclassifying ALT-supporting reads"* (CHANGELOG ≈ 710) plus *"Large deletions (≥50bp) with low reciprocal overlap (<50%): previously silently dropped, now flagged for Phase 3 arbitration"* (≈ 151) show the team was fixing a **false-REF / dropped-ALT** problem on large deletions. The tolerant branch (commit `de7da8e4`, 2026-02-22) was the pragmatic way to capture imperfect-length ALT reads. The cost — skipping the sequence check — was accepted to recover sensitivity. The gap is that it traded a false-negative for a potential false-positive without a partial sequence check.

**Root cause.** Tolerant length acceptance was bolted on without a corresponding *partial* sequence concordance check over the overlapping span.

**Fix.** Keep the tolerant length acceptance (don't regress the sensitivity win), but require the **overlapping portion** of the deleted reference to match:
1. Compute the overlap span between observed deletion `[obs_start, obs_end)` and expected `[exp_start, exp_end)`.
2. Verify the read's flanking sequence is concordant with `expected_del_seq` over the shared coordinates (or, minimally, require the breakpoints to be within a tolerance proportional to `repeat_span`).
3. Set `del_ok = (overlap_fraction >= 0.5) && partial_seq_match`.
Apply identically to the strict path at `:1488`.

**Tests.** (a) A genuine large deletion with ±N bp ragged breakpoints in a repeat → still counted (sensitivity preserved). (b) An unrelated large deletion of different sequence with ≥50 % length overlap → **not** counted (false-positive closed). Add both to the large-deletion fixture set.

**Regression risk.** Medium — must verify the sensitivity cases from the v4.x fix still pass. Pull those fixtures forward as guardrails.

**Effort.** M (1–2 days; the partial-overlap concordance is the substantive part).

---

## CR-5 — mFSD statistics: ±∞ LLR and small-N invalid KS p-value (drive clinical labels)

**Locations:** `rust/src/counting/mfsd.rs:145-153` (LLR `log(0)`), `:280-294` (KS p-value), consumed by `mfsd_report.py:_classify_origin` (TUMOR-LIKE / CH-LIKE at p<0.05).

**Issue — two coupled statistical defects in the values that drive the fragment-origin call:**

**(a) LLR can be ±∞ from a single fragment.** The per-fragment term guards only the denominator: `ratio = if p_healthy < EPSILON { INFINITY } else { p_tumor/p_healthy }; ratio.ln()`. For a tail fragment where `p_healthy` underflows below `f64::EPSILON`, `ratio = +∞` and `ln(∞) = +∞`, pinning the whole `mfsd_alt_llr` sum to `+Infinity`. There is no symmetric guard for `p_tumor` underflow (→ `ln(0) = −∞`). One outlier fragment dominates the metric.

**(b) KS p-value uses the asymptotic Kolmogorov series**, `Q(λ)=2Σ(−1)^{k−1}e^{−2k²λ²}` with `λ=D·√(nm/(n+m))` (`:280`). This is the large-sample limit and is materially inaccurate for `n_alt = 5–20`, the exact low-input cfDNA regime this module targets. The reported p — which thresholds the clinical TUMOR-LIKE/CH-LIKE label — is not trustworthy at the N where decisions are made.

**Git archaeology.** Both from `a9d8a00` (2026-03-04, *"native mFSD integration — 31 MAF cols + 7 VCF INFO fields"*). The `log(0)` comment literally reads *"Guard: clamp denominator to avoid log(0)"* — the author **knew** about `log(0)` but fixed only one side. The LLR Gaussian params are documented *"calibrated to the MSK-ACCESS cohort"* but **no fitting code is committed** (see HI-10). `MIN_FOR_KS=5`'s comment *"insufficient data rather than a spurious result"* shows small-N awareness, but 5 is a floor, not a power-justified threshold (see ME-9).

**Root cause.** (a) The metric is computed as a ratio of densities (which underflow) instead of a difference of log-densities (closed-form, finite). (b) The asymptotic approximation is applied outside its valid range.

**Fix.**
1. **Closed-form log-space LLR.** The Gaussian log-ratio is `0.5·(z_h² − z_t²) + ln(σ_h/σ_t)` with `z = (x−μ)/σ` — finite for all finite `x`. Replace the pdf division entirely:
   ```rust
   let zt = (x - tumor_mu)/tumor_sigma;
   let zh = (x - healthy_mu)/healthy_sigma;
   let log_ratio = 0.5*(zh*zh - zt*zt) + (healthy_sigma/tumor_sigma).ln();
   sum += log_ratio;   // optionally clamp to ±MAX_PER_FRAG as defense in depth
   ```
2. **Report mean LLR per fragment**, not the raw sum, so magnitude is comparable across variants of different ALT depth.
3. **Small-N KS.** For small `n·m`, compute the exact two-sample KS null (Hodges lattice-path, as SciPy `method="exact"`). If exact is too costly, apply Stephens' finite-sample correction before the series: `λ_eff = (√Nₑ + 0.12 + 0.11/√Nₑ)·D`, `Nₑ = nm/(n+m)`. Document the chosen method and its valid range.

**Tests.** (a) A single 500 bp ALT fragment must not pin LLR to ∞; assert finiteness. (b) Golden-value test of the closed-form log-ratio vs a high-precision reference at several sizes. (c) KS p-value vs SciPy `ks_2samp(method="exact")` for n=m∈{5,8,12,20}; assert within tolerance. The D-statistic merge-walk (`:255-268`) is already correct — keep its tests.

**Regression risk.** Medium — LLR magnitudes and KS p-values will shift (correctly). Re-baseline any golden mFSD outputs; coordinate with whoever consumes the report labels.

**Effort.** M (1–2 days incl. exact-KS implementation + validation).

---

# HIGH

## HI-1 — Sample failures swallowed; exit code 0 even when all samples fail

**Status: Done (2026-07-01).** Added `_exit_on_sample_failure(result)` (cli.py), called by
both `dna` and `rna` *outside* the command's `try/except` (so its `typer.Exit` isn't
swallowed). It exits **code 1** only when a sample actually **failed** (is in
`failed_samples` — a Rust panic surfaced as `PyErr`, an unreadable BAM, …), logging
all-failed vs partial. An **empty variant set is deliberately NOT a failure**: a sample can
legitimately have no variants called, and per-sample Nextflow tasks must not fail on that —
those runs process zero samples with no `failed_samples` and exit `0`. The two duplicated
`Pipeline().run()` blocks collapsed to `result = Pipeline(config).run()` + the shared helper
(no duplication). Tests: `tests/test_hi1_exit_codes.py` (helper unit branches + end-to-end
empty-variant-file→0 + mocked success/failure through the CLI).

**Locations:** `src/gbcms/pipeline.py:491-493` (catch/append/return), `src/gbcms/cli.py:482` (result ignored).

**Issue.** `_process_sample` catches `Exception`, logs `logger.error`, appends to `_failed_samples`, and returns. `run()` logs a summary but `cli.py` never inspects `failed_samples`, so a run where every BAM failed inside Rust (e.g. a panic surfaced as `PyErr`) exits `0`. Under Nextflow this masks systematic failure as success.

**Root cause.** No propagation of partial/total failure to the process exit code.

**Fix.** In `cli.py` after `pipeline.run()`, inspect returned stats: `raise typer.Exit(code=1)` when `failed_samples` is non-empty (configurable: `--allow-partial-failure` to downgrade to a warning), and always exit non-zero when `samples_processed == 0`. Print a clear per-sample failure summary.

**Tests.** CLI test with one deliberately broken BAM → exit 1 and the other samples still written. Test `--allow-partial-failure` downgrades to 0 with warning.

**Regression risk.** Low (behavioral; document in release notes — CI pipelines that ignored failures will now see non-zero).

**Effort.** S.

---

## HI-2 — Read-level DP/AD/RD double-count supplementary/secondary when their filters are disabled

**Locations:** `rust/src/counting/engine.rs:1205-1213` (DP), `:1308-1349` (RD/AD); filter gate `shared/filters.rs:67-76`; defaults `cli.py:237-238` (both `True`).

**Issue.** Fragment-level counts collapse correctly (supplementary/secondary share QNAME → same `mol_hash`). But **read-level** DP/RD/AD increment once per record. With `--no-filter-supplementary` (or secondary), a supplementary/chimeric segment overlapping the same anchor as its primary inflates DP and AD/RD, skewing read-level VAF = AD/DP and breaking the `DP ≥ RD+AD+…` intent. Safe by default (both filters `True`), but the failure is silent when toggled.

**Git archaeology.** CHANGELOG (≈ 340, 391) shows filter defaults were deliberately set to `True` for secondary/supplementary/qc-failed — the team chose safe defaults rather than fixing the read-level path, leaving the hazard latent when a user overrides.

**Root cause.** Read-level counting is per-record; supplementary/secondary are not first-class fragment observations and should never increment read-level depth at the anchor.

**Fix.** In the Phase-1 read loop, unconditionally skip read-level DP/RD/AD increments for `record.is_supplementary() || record.is_secondary()`, independent of the user flag (retain them only for haplotype evidence if desired). Keep the flags governing whether such records contribute at all.

**Tests.** Synthetic read with a supplementary segment overlapping the anchor, filters off → DP increments once. Add to the counting fixture matrix.

**Regression risk.** Low; only affects the non-default flag combination.

**Effort.** S.

---

## HI-3 — `OFF_TARGET_THRESHOLD = 20` is not read-length-aware

**Locations:** `rust/src/counting/wfa_router.rs:32`, used at `:99`.

**Issue.** A fixed edit-distance cutoff of 20 means 40 % divergence for a 50 bp read but <10 % for a 250 bp read (ref_context can reach ~400 bp per `pangenome.rs:44`). Legitimate on-target reads on long windows exceed it; the doc comment ("too divergent to classify reliably") gives no derivation.

**Git archaeology.** Same `761b98b` WFA-introduction commit; a placeholder magic number that shipped.

**Fix.** Make it length-relative: `let off_target = max(10, read_len/10) as i32;` (or derive from expected per-base error rate). Document the basis. Coordinate with CR-2's length-aware off-target gate (they touch the same branch).

**Tests.** Reads of 50/150/250 bp at fixed % divergence classify consistently relative to the threshold.

**Regression risk.** Low-medium; re-validate off-target rate on a panel.

**Effort.** S.

---

## HI-4 — PairHMM `ln(0)` at Q0 → NaN LLR; clamp `p_correct`

**Locations:** `rust/src/counting/pairhmm.rs:80-83, 100-102`; downstream comparison `:382-384, :492-494`.

**Issue.** `p_correct = 1.0 − 10^(−q/10)` is `0.0` at Q0 → `LogProb(ln 0) = −∞`. A single Q0 *matching* base zeroes that path; for both REF and ALT this drives `ll_alt, ll_ref → −∞` and `llr = −∞ −(−∞) = NaN`. `NaN > thr` and `NaN < −thr` are both false, so it falls through to the ambiguous branch and returns "neither." Not a crash, but Q0 matching bases silently collapse reads, and the intermediate NaN is fragile. The `usable_count < 3` gate does not protect (it counts bases `≥ min_baseq`; Q0 bases can still be emitted).

**Fix.** Clamp quality-derived error to a known range, GATK-style: `let p_err = (10f64.powf(-(q as f64)/10.0)).clamp(1e-6, 1.0 - 1e-6); let p_correct = 1.0 - p_err;`. Apply at both emit sites.

**Tests.** A read with a Q0 matching base must not produce NaN LLR; assert finite and that classification matches the same read at Q2.

**Regression risk.** Low (numeric hardening).

**Effort.** XS.

---

## HI-5 — N-base emission is not truly LLR-neutral near indels

**Locations:** `rust/src/counting/pairhmm.rs:74-76, 95-99`.

**Issue.** N bases return `Mismatch(LogProb::from(Prob(0.25)))`, intended to be "net-zero LLR." This cancels only when the N aligns to the same column in both haplotypes; under semiglobal alignment with free end gaps, an N near the indel can align to different positions in the REF-optimal vs ALT-optimal paths, leaking a nonzero LLR. `prob_emit_y` returns a fixed `0.999`, making the gap bookkeeping for N asymmetric.

**Fix.** Make N position-independent and symmetric across both haplotypes (or exclude N positions from the aligned read prior to DP). Add a unit test asserting `|LLR| < 1e-6` for an N at the discriminating locus of a nontrivial indel. The current `test_classify_pairhmm_n_base_neutral` only asserts "not confident ALT" — strengthen it.

**Regression risk.** Low; mostly affects duplex-masked N-heavy reads.

**Effort.** S-M.

---

## HI-6 — `check_complex` insertion off-by-one at the exclusive REF end

**Locations:** `rust/src/counting/variant_checks.rs:586-598`.

**Issue.** Phase-1 reconstruction includes inserted read bases whenever `ref_pos >= start_pos && ref_pos <= end_pos` — inclusive of `end_pos`, which is the *exclusive* REF end. An insertion sitting exactly at `end_pos` is appended into the reconstructed haplotype though it lies outside `[start_pos, end_pos)`, inflating `recon_len` and potentially flipping a length-based REF/ALT decision for complex variants.

**Fix.** Use `ref_pos >= start_pos && ref_pos < end_pos` for the insertion-inclusion condition (or explicitly justify including the trailing insertion and adjust the ALT/REF length comparison accordingly). See ME-3 — unify with the window helper.

**Tests.** Complex variant with an insertion immediately after the REF span; assert correct REF vs ALT classification.

**Regression risk.** Medium for complex variants — re-run the complex/DelIns fixture suite.

**Effort.** S.

---

## HI-7 — `gene_strand` never populated → sense/antisense meaningless, `enforce_strandedness` a silent no-op

**Locations:** `rust/src/types.rs:33` (defaults `None`), `src/gbcms/pipeline.py:225-233` & `normalize.py:69` (never set), `rust/src/counting/rna.rs:95-99` (`is_sense_strand` returns `true` when `None`), consumers `engine.rs:1126, 1353-1360, 2292, 2591`.

**Issue.** `Variant.gene_strand` is never populated from the GTF anywhere in the pipeline. With `gene_strand == None`, `is_sense_strand` returns `true` for **all** reads, so `sense_depth` accumulates everything, `antisense_depth` stays 0, and `enforce_strandedness` filters nothing. The `AnnotationIndex` *has* per-exon strand (`mod.rs:58`) but exposes no "strand at position" API, and the engine never derives it. The entire sense/antisense subsystem and the dUTP strandedness filter are inert.

**Git archaeology.** The field was added in `761b98b` (2026-03-18) with the comment *"None in DNA mode — zero cost."* The RNA strandedness machinery (CHANGELOG 5.0.0 ≈ 185–191: amplicon auto-disable, *"per-read strand tracking in junction accumulator"*) was built around it, and `enforce_strandedness` was deliberately flipped to `true` (≈ 342, 392). But the population step was never wired — the per-read strand tracking added for ASJD is a *different* mechanism and does not feed `gene_strand`. This is an incomplete feature, and because the default was flipped to `true`, it's a silent no-op in production RNA runs.

**Root cause.** Missing GTF→variant strand derivation; no lookup API on `AnnotationIndex`.

**Fix.**
1. Add `AnnotationIndex::gene_strand_at(chrom, pos) -> Option<char>` (resolve from overlapping exons; on conflicting-strand overlap return `None` and `warn!`).
2. In RNA mode, populate `variant.gene_strand` before counting — either in Rust right after building the index, or in `pipeline.py` from the parsed GTF.
3. Emit a `warn!` when `enforce_strandedness == true` but `gene_strand` is `None` for a variant, so the no-op is never silent.

**Tests.** RNA fixture on a known minus-strand gene: assert `antisense_depth > 0` and that `enforce_strandedness` drops antisense reads. Test conflicting-strand overlap → `None` + warning.

**Regression risk.** Medium-high for RNA mode — this *activates* a filter that was inert. RNA outputs will change. Stage behind a clearly-noted minor release; provide `--no-enforce-strandedness` escape hatch (already exists).

**Effort.** M.

---

## HI-8 — Spliced-read fragment size inflated by full intron length

**Locations:** `rust/src/counting/mfsd.rs:181-217` (`calc_physical_insert_size`).

**Issue.** Physical size = `|TLEN| − Σ D + Σ I`, iterating CIGAR for `Del`/`Ins` only; **`Cigar::RefSkip` (N) falls into `_ => {}`** and is ignored. BAM TLEN is the genomic reference span, which for a splice-spanning pair includes the entire intron. So spliced reads report fragment sizes inflated by the summed intron length (often hundreds–thousands of bp), corrupting every RNA mFSD statistic (means, KS, LLR, nucleosomal fractions). The ref-context splicing path (`rna.rs:292`) *does* discount introns — the two paths are inconsistent.

**Git archaeology.** `a9d8a00`/`d358537` mFSD work was validated on cfDNA indel BAMs (per doc comments), which have no N ops — the splice case was never exercised.

**Fix.** Add `Cigar::RefSkip(n) => skip_bp += *n as i32,` and compute `physical = abs_tlen − del_bp − skip_bp + ins_bp`. Update the pathological-clamp guard accordingly. Add a regression test with a spliced read pair.

**Regression risk.** Low for DNA (no N ops); changes RNA mFSD numbers (correctly).

**Effort.** S.

---

## HI-9 — Splice-motif classifier ignores transcript strand

**Locations:** `rust/src/counting/engine.rs:2483-2521` (`classify_splice_motif`), diagnostic `:2721-2723`.

**Issue.** Donor/acceptor dinucleotides are read as forward-genomic 2 bp and matched against `GT-AG`/`GC-AG`/`AT-AC`. For a minus-strand gene a canonical intron appears as `CT…AC` on the forward genome → classified `OTHER` → spurious `NON_CANONICAL_MOTIF` on every real minus-strand junction. Strand is never passed in (and per HI-7 isn't even available).

**Git archaeology.** `996a1ad` (*"real splice motif classification via FASTA + NON_CANONICAL_MOTIF flag"*) — forward-strand-only implementation.

**Fix.** Resolve junction/gene strand (via HI-7's `gene_strand_at`) and, for minus strand, reverse-complement: donor = revcomp(2 bp at `junction_end−2`), acceptor = revcomp(2 bp at `junction_start`). Equivalently accept a junction as canonical if either the forward read or its reverse complement matches a canonical motif.

**Tests.** Minus-strand canonical junction → `GT-AG`, no `NON_CANONICAL_MOTIF`. Depends on HI-7.

**Regression risk.** Medium for RNA ASJD outputs.

**Effort.** S-M (after HI-7).

---

## HI-10 — LLR model hardcoded/uncited; report mislabels it as the PairHMM LLR

**Locations:** `rust/src/counting/mfsd.rs:45-66` (`LlrModelParams::human`), `src/gbcms/report/mfsd_report.py:45-48, 733`.

**Issue.** `P_healthy=N(167,30)`, `P_tumor=N(145,35)` are hardcoded, not sample-fit, not CLI-configurable, and the two densities are near-degenerate (means 22 bp apart, overlapping σ). The report describes this value as *"Log-likelihood ratio from PairHMM alignment"* and *"positive = tumor-like"* — but it is the **fragment-size Gaussian** LLR, unrelated to the PairHMM classifier LLR (`hmm_llr_threshold=2.3` at `engine.rs:228`). Two distinct quantities are conflated in user-facing text.

**Git archaeology.** `a9d8a00` comment claims *"calibrated to the MSK-ACCESS cohort"* but no fitting code or citation is committed.

**Fix.** (a) Correct the report glossary/tooltip to describe it as a fragment-size Gaussian LLR, distinct from the HMM LLR. (b) Expose `healthy_mu/sigma`, `tumor_mu/sigma` as CLI/config options. (c) Ideally fit `P_healthy` from the sample's own REF fragment-size distribution so the ratio is sample-relative. (d) Add a literature/cohort citation for the defaults.

**Regression risk.** Low (config + docs); fitting from sample REF changes values — gate behind a flag.

**Effort.** S (text+config) / M (with sample-fit).

---

## HI-11 — Six KS p-values per variant get no multiplicity correction

**Locations:** `rust/src/counting/engine.rs:1483-1488, 1988-1993`; existing `benjamini_hochberg` at `shared/stats.rs:189` (applied only to ASJD, RNA mode, `engine.rs:711-723`).

**Issue.** Each variant runs 6 KS tests; the report classifies on raw `ks_pval < 0.05` across many variants with no within-variant or across-variant FDR. Expected false TUMOR-LIKE calls scale with variant count.

**Fix.** Collect the classification-driving comparison's p-values (`mfsd_pval_alt_ref`) across all variants → one BH pass → store `mfsd_qval_alt_ref`; have `_classify_origin` threshold on the q-value. Reuse the verified `benjamini_hochberg`. Optionally correct the other 5 comparisons too. Document the correction scope.

**Tests.** Synthetic set of N null variants; assert corrected false-positive rate near α. Reuse the BH unit test.

**Regression risk.** Medium — fewer TUMOR-LIKE calls; re-baseline report expectations.

**Effort.** S-M.

---

## PF-1 — mFSD stats + size arrays computed for every variant on every run

**Status: Done (PR #54, 2026-06-30).** Added an `mfsd` flag plumbed
`OutputConfig.mfsd → count_bam_binned → count_variant_from_cache`; the size-array
reservation, the extracted `compute_mfsd_stats` helper, and the post-counting mFSD
BH-FDR pass are all gated on it. The legacy oracle still computes mFSD (∉ `PARITY_FIELDS`).
Validated count-neutral on 3,040 real cfDNA variants (all 246 non-mFSD columns
byte-identical with mFSD on vs off). Kept for **per-process memory ×N** under fan-out, not
CPU (see §"M5 — empirical scoping").

**Locations:** signature `rust/src/counting/engine.rs:466` (no `mfsd` param), compute block `:1459-1505`, arrays `:1504-1505`; `types.rs:255` ("Populated in all runs").

**Issue.** Because no `mfsd` flag reaches Rust, the engine unconditionally runs 6 KS tests, 2 LLR passes, and 4 `calc_fraction_in_range` scans per variant, and fills `ref_sizes`/`alt_sizes: Vec<u32>` (≈ DP entries each) on every `BaseCounts` — even on a plain `gbcms dna` run that never uses them. At thousands of variants × deep cfDNA, that's wasted CPU plus tens-to-hundreds of MB of FFI-resident memory held through output writing.

**Fix.** Add `mfsd: bool` to `count_bam_binned`/`count_bam`; gate the entire block `:1459-1505` behind `if mfsd`. Plumb `self.config.output.mfsd` from `pipeline.py:427`. This is also the enabling change for ME-1 and HI-7-style "engine is output-aware" hygiene.

**Tests.** Assert size arrays empty and mFSD fields default when `mfsd=false`; unchanged when `true`. Add a memory/timing micro-benchmark on a 5 000-variant panel.

**Regression risk.** Low.

**Effort.** S.

---

## PF-3 — htslib decode threads never enabled

**Status: Dropped (M5 empirical scoping, 2026-06-30).** Fetch is ~89% of cfDNA counting
but already parallelized across ~201 bins at ~84% efficiency — there are no idle cores to
feed. Under N concurrent 4-core Nextflow processes, per-reader decode threads would
oversubscribe the node. Only helps single runs with fewer deep bins than cores. Not worth
the oversubscription/parity risk; the thread-budget contract (LO-14) explicitly forbids
adding decode threads on top of the rayon budget. See §"M5 — empirical scoping".

**Locations:** `rust/src/counting/engine.rs:346, 577, 625` (reader opens); `threads` consumed only by rayon at `:320-323, :607-612`.

**Issue.** No reader calls `reader.set_threads(...)`, so BGZF/CRAM block decompression is single-threaded per reader. Fine when bins ≥ workers, but on the ultra-deep long-pole bin (PF-2) decompression runs serially on one core.

**Fix.** Call `reader.set_threads(k)` on each thread-local reader, budgeting `rayon_workers × k ≤ threads` to avoid oversubscription (couple with LO-14). Document the split. Particularly benefits the long-pole bin.

**Regression risk.** Low; verify total thread count stays within `task.cpus`.

**Effort.** S.

---

# MEDIUM

## ME-1 — `sub_nuc`/`mono_nuc` fields computed but dropped from VCF
**Status: Done (2026-07-01).** Chose VCF↔MAF parity (they're computed under `--mfsd` and
already in MAF — output them, don't discard). Added 5 INFO fields to the VCF mFSD block —
`MFSD_SUB_NUC_REF_FRAC`, `MFSD_SUB_NUC_ALT_FRAC`, `MFSD_SUB_NUC_ENRICHMENT`,
`MFSD_MONO_NUC_REF_FRAC`, `MFSD_MONO_NUC_ALT_FRAC` (header + data row), taking VCF from 8 to
**13** mFSD INFO fields. Also corrected pre-existing count drift in docs/CLI help (MAF is
**41** mFSD columns, not 34; VCF is **13**, not 7) and added the new rows to the VCF field
tables. Tests: `test_vcf_writer_mfsd_info_when_enabled` (13 header lines) +
`test_vcf_writer_mfsd_subnuc_values_in_data_row`.

**Locations:** computed `engine.rs:1494-1502`, MAF `output.py:583-587`, VCF block `output.py:916-925` (omits them). **Fix:** either add `MFSD_SUB_NUC_*`/`MFSD_MONO_NUC_*` INFO fields to `_write_header`/`write`, or skip computing them when output is VCF (couple with PF-1's gate). Decide intended VCF surface. **Effort:** S.

## ME-2 — Transcript-count delimiter `;` vs documented `|`
**Status: Done (2026-07-01).** Standardized on `|` (the VCF-INFO-safe choice — `;` is the
VCF INFO separator, which is why the writer had to repair `;`→`|` for VCF only, leaving MAF
inconsistent). The engine now joins transcripts with `|` (`engine.rs`), so MAF, VCF, and the
documented `ENST:AD,RD,DP|…` header all agree; the redundant Python `;`→`|` repair for
TXRC/TXFC was dropped. Other `;`-joined fields (gbcms_status, diagnostics) are unchanged.
Docs updated (rna-annotation.md, output-formats.md MAF section, the Rust doc-comment). Test:
`tests/test_rna_output.py::test_transcript_counts_use_pipe_delimiter_in_vcf`.

**Locations:** join `engine.rs:2385`; header says `|` `output.py:819`; VCF repairs `:955`; MAF writes raw `;` `:485`. **Fix:** join with `|` in Rust and drop the Python `replace`, or standardize on `;` and fix the header. One delimiter. **Effort:** XS.

## ME-3 — Soft-clip boundary differs between `check_complex` Phase 1 and `extract_raw_read_window`
**Locations:** `variant_checks.rs:609` (`< end_pos`) vs `alignment.rs:105` (±1 bp window). A clip carrying ALT just upstream of the anchor is dropped in Phase 1 but picked up in Phase 3 → path-dependent classification. **Fix:** extract one shared boundary helper used by both Phase 1 reconstruction and Phase 3 raw extraction (also resolves HI-6 consistently). **Effort:** S.

## ME-4 — Haplotype parity collisions degrade locus to universal ties
**Locations:** `pangenome.rs:309-336`; consumers assume even-index=REF/odd=ALT (`wfa_router.rs:70-85`, `pairhmm.rs:464-480`). When a REF-class haplotype byte-string equals an ALT-class one, every read scores an exact tie → all "neither." **Fix:** detect when a haplotype appears in both classes (locus genuinely non-discriminating) and short-circuit explicitly rather than emitting universal ties. **Effort:** S-M.

## ME-5 — `prob_emit_y = 0.999` unjustified/asymmetric
**Locations:** `pairhmm.rs` `prob_emit_y`. A fixed 0.999 emission for every haplotype base in an insertion context is arbitrary and asymmetric vs the quality-derived emit_x, biasing insertion-vs-deletion gap arms. **Fix:** document/derive it; if intent is "base certainly present," use `LogProb::ln_one()`; if an emission over {A,C,G,T}, use 0.25. **Effort:** S.

## ME-6 — Contig normalization is `chr`-strip only (`MT`/`chrM`)
**Locations:** `gtf.rs:99-102`, `rna.rs:247`, lookups `engine.rs:1375, 2225, 2557`. `trim_start_matches("chr")` maps `chrM→M` but a GRCh37 `MT` BAM never matches → all chrM annotation/editing/ASJD silently empty. **Fix:** central `normalize_contig()` mapping `chrM/M ↔ MT` (+ common aliases, case-fold), applied at parse and every lookup, including caller-normalized `variant_chroms`. **Effort:** S.

## ME-7 — `STRAND_DISCORDANT` uses raw genomic strand without R1/R2 folding
**Locations:** `engine.rs:2618-2634, 2725-2739`. Partitioning junction reads by `record.is_reverse()` without dUTP R1/R2 folding makes a normal FR pair contribute one forward + one reverse → minority fraction ≈ 0.5 → false `STRAND_DISCORDANT`. **Fix:** convert each read to inferred transcript strand (same logic as `is_sense_strand`) before tallying, and dedupe by QNAME. Depends on HI-7. **Effort:** S.

## ME-8 — ASJD BH includes `p=1.0` padding in `n`
**Locations:** `engine.rs:712-715`; BH uses `n = pvalues.len()` (`stats.rs:190,209`). Untested variants carry `asjd_pval = 1.0`, inflating `n` → over-conservative q-values, real junctions lost. **Fix:** run BH only over variants with `asjd_n_alt_junc + asjd_n_ref_junc > 0`, scatter q-values back, leave others at NaN/1.0 sentinel. **Effort:** S.

## ME-9 — `MIN_FOR_KS=5` vs report `min_alt=3`
**Locations:** `mfsd.rs:33`, `mfsd_report.py:197, 106`. Variants with 3–4 ALT fragments are shown/classified but always have NaN KS; `_classify_origin` treats NaN as non-significant → can label CH-LIKE on a missing test. Also, at n=m=5 the test has essentially no resolution. **Fix:** raise `MIN_FOR_KS` to a power-justified value (≥8 typical; document min detectable D); make report `min_alt ≥ MIN_FOR_KS`; have `_classify_origin` return `INSUFFICIENT` on NaN rather than "not significant." **Effort:** S.

## ME-10 — `sub_nuc_enrichment` NaN is ambiguous
**Locations:** `engine.rs:1496-1500, 2001-2005`. Enrichment is NaN whether ALT is empty or REF is empty; the report renders both as "N/A," and classification can't distinguish "no ALT" from "no REF." **Fix:** gate on both counts and emit `INSUFFICIENT` (distinct from generic NaN) when ALT is empty; document the asymmetry. **Effort:** XS.

## ME-11 — Report formatters don't guard `isinf`
**Locations:** `mfsd_report.py:325-338` (`_fmt_val`/`_fmt_pval`). With CR-5 unresolved, a `+Inf` LLR renders literal "inf" into clinical HTML. **Fix:** `if math.isnan(v) or math.isinf(v): return "N/A"`. Defense-in-depth even after CR-5. **Effort:** XS.

## ME-12 — "Majority rule" is structural-ALT-wins consensus
**Locations:** `shared/fragment.rs:194-222`. When R1/R2 disagree: structural ALT (indel CIGAR) wins unconditionally (`:206`), else higher-BQ wins past threshold, else discard (in DPF, not RDF/ADF). A single mate with a spurious indel CIGAR op overrides a high-BQ REF mate. **Fix:** (a) relabel docs/`BaseCounts` comments from "Majority Rule" to "quality-weighted consensus with discard band"; (b) optionally require the structural ALT not be contradicted by a high-BQ REF on the other mate before it wins. **Effort:** S (docs) / M (logic change + fixtures).

## PF-2 — One rayon task per bin, no depth cap → long-pole tail
**Status: Deferred — niche (M5 empirical scoping, 2026-06-30).** cfDNA has no long-pole (busiest of ~201 bins = 3% of bin-work, top-5 = 13%); RNA skew is extreme but trivial in absolute time (~40ms of bin-work total). The cost-sort is a cheap future win, not a milestone driver.

**Locations:** `engine.rs:621` (`bins.par_iter()`), `:142-198`. `BIN_MAX_VARIANTS` caps variants, not depth; a 1-variant ultra-deep bin is one indivisible mega-task that serializes the tail. **Fix:** cheap win — `bins.sort_by_key(|b| Reverse(estimated_cost(b)))` (variants × span, or a depth proxy) before `par_iter` (longest-processing-time-first). Stretch — split very deep bins by read sub-range. **Effort:** S (sort) / M (split).

## PF-4 — rayon pool rebuilt per `count_bam_binned` call
**Status: Dropped (M5 empirical scoping, 2026-06-30).** Core gbcms processes one sample per invocation, so the pool is built once per run; only the MNP-rescue 2nd `count_bam_binned` call rebuilds it. Moot at cohort scale, where parallelism is Nextflow's N concurrent processes, not repeated in-process calls.

**Locations:** `engine.rs:607-612` (and `:320-323`). New OS thread pool per sample, ×2 with MNP rescue (`pipeline.py:697`). **Fix:** build once (`OnceCell`/`build_global`) and `install`, keyed by thread count. Pool scoping is otherwise correct (no oversubscription bug). **Effort:** S.

## ME-13 — Adjacent bins' padded fetch ranges overlap → double I/O
**Status: Dropped (M5 empirical scoping, 2026-06-30).** Measured **0.00% re-fetch** on real cfDNA — adjacent bins never overlapped in practice — and this is the single highest parity-risk change in M5. Not worth it.

**Locations:** `engine.rs:156, 182-183, 189-194`. Bins partition variants but inflate fetch on both ends; adjacent same-chrom bins overlap → overlap reads fetched/decoded twice. Bounded (padding vs 10 kb) but real for clustered MAF input. **Fix:** clamp each bin's fetch `start` to `max(start, prev_bin_end_on_same_tid)` so overlaps are fetched once. Correctness unaffected. **Effort:** S.

---

# LOW

## LO-1 — Duplicated `.pyi` stubs; `count_bam` stub missing params
`src/gbcms/_rs.pyi` and `src/gbcms_rs.pyi` are byte-duplicates kept in sync by hand; both `count_bam` stubs omit `reference_fasta`/sibling params present in the real signature (`engine.rs:272`). **Fix:** generate one stub, symlink/re-export the other; regenerate from the actual `#[pyo3(signature)]`. **Effort:** XS.

## LO-2 — Legacy `count_bam` (one fetch per variant) still exported
`lib.rs:15` exports `count_bam`; pipeline only calls `count_bam_binned`. It fetches per variant (`engine.rs:1574`) — the O(seeks) pattern binning replaced — and doubles maintenance for parity. **Fix:** gate behind a `test-only` feature or remove after parity sign-off; keep one as the parity oracle in tests. **Effort:** S.

**Status: Resolved.** Gated behind the default `legacy-parity` Cargo feature
(`rust/Cargo.toml`): `count_bam` + `count_single_variant` + their re-export/registration
are `#[cfg(feature = "legacy-parity")]`. Dev/test builds keep it (default on) so the
binned↔legacy parity tests still run; release wheels build `--no-default-features`
(Dockerfile, release.yml) to exclude it from the shipped surface. Safe now that the
parity oracle is codified in CI (the parity gate, formerly an open follow-up).

## LO-3 — Bin span can chain well past 10 kb
**Status: Document, don't cap (M5 empirical scoping, 2026-06-30).** Chaining is bounded in practice (≤ ~2.5× the 10 kb window); a hard cap risks re-breaking CR-1's anchor coverage. Document the soft floor rather than capping.

`engine.rs:167,183` — repeated `var_end + window/2` extension can grow a bin far beyond 10 kb until the 200-variant cap trips. Over-inclusion (correctness-safe) inflating per-bin fetch/compute. **Fix:** cap `bin_end` growth or document 10 kb as a soft floor. (Interacts with CR-1 and ME-13.) **Effort:** S.

## LO-4 — u64 QNAME hash stores no key
`shared/fragment.rs:244-249`; map keyed on SipHash u64 with no stored key (`engine.rs:1036`). Birthday collision ≈ 5e-10 at realistic depth — negligible. **Fix:** none required; revisit only if per-locus fragment counts grow orders of magnitude (then store the QNAME for tiebreak). **Effort:** n/a (note only).

## LO-5 — "orientation-aware hashing" is a misnomer
`shared/fragment.rs:20-22` comments imply orientation is in the hash key; it is QNAME(+UMI) only, with orientation resolved downstream. **Fix:** correct the doc comments. **Effort:** XS.

## LO-6 — FSB strand from best-quality read
`shared/fragment.rs:226-238` — the fragment's FSB strand is taken from whichever mate had the higher BQ at the locus (quality-driven, not deterministic) for symmetric overlaps. Generally fine; minor FSB noise. **Fix:** document, or pick a deterministic tiebreak (e.g. R1) on equal quality. **Effort:** XS.

## LO-7 — `median_qual` upper-median bias
`bam_utils.rs:66` returns `filtered[len/2]` (upper median on even counts), biasing reported quality slightly high. Cosmetic. **Fix:** average the two middle elements for even counts if exactness matters. **Effort:** XS.

## LO-8 — Deletion left-anchor invariant undocumented
`variant_checks.rs:1427-1430` — anchor logic assumes `variant.pos` is the left-anchor base (never inside the deletion). Load-bearing and undocumented; a non-left-anchored deletion breaks it. **Fix:** document the invariant and assert it (ties to CR-3's guard). **Effort:** XS.

## LO-9 — Empty/exon-less GTF degrades silently
`gtf.rs:162-167` — parser ingests only `exon` records; a GTF with genes/transcripts but no/malformed exons yields an empty `AnnotationIndex` with only a warn. Can't distinguish "no exons" from "wrong build." **Fix:** track and report counts of `gene`/`transcript`/`exon` lines seen; optionally hard-error on RNA-mode-with-GTF-but-zero-exons behind a flag. **Effort:** S.

## LO-10 — `.`/unstranded exon defaults to `+`
`gtf.rs:117-120` — silent default hides malformed strand. Minor until HI-7 makes strand consumed. **Fix:** warn on unexpected strand tokens; carry `Option<char>` rather than defaulting. **Effort:** XS.

## LO-11 — Editing-site flag is positional-only
`engine.rs:1374-1384`, loader discards Ref/Ed bases (`rna.rs:247-250`). Any variant at an editing coordinate gets flagged regardless of being A>G/T>C. Flag-only (no count impact). **Fix:** store the edited base pair and require A>I consistency before flagging, or document positional-only semantics. **Effort:** S.

## LO-12 — NH-tag rescue only matches `U8`/`I32`
`rna.rs:48-67` — NH:i:1 rescue ignores `U16/I8/I16/U32` Aux encodings → some uniquely-mapped low-MAPQ reads dropped in RNA mode. **Fix:** match all integer Aux variants and compare `== 1`, or use an integer-coercing accessor. **Effort:** XS.

## LO-13 — KDE bandwidth clamp (5.0 bp) arbitrary
`mfsd_report.py:160-173` — zero-variance size class clamps bandwidth to 5.0 bp, undocumented in the methodology. **Fix:** document the clamp and its rationale (or derive from data range). **Effort:** XS.

## LO-14 — `threads` not clamped to cores

**Status: Done (PR #53, 2026-06-30).** Recast as a hard, validated **thread-budget
contract**: `--threads` is the total worker budget for one process (`min=1` on the CLI),
all rayon pools are sized through `shared::resolve_thread_budget` (which also guards the
`num_threads(0)`=all-cores foot-gun), and any future htslib decode threads must
*subdivide* this budget, never add to it — so gbcms stays a good citizen under Nextflow
fan-out (N concurrent 4-core processes). Documented in `architecture.md` design decision #1.

`pipeline.py:439` passes `config.threads` verbatim to Rust; direct-CLI users or multiple manual processes can oversubscribe (compounds with PF-3 if htslib threads added). **Fix:** clamp to `available_parallelism()`; budget `rayon × htslib ≤ cpus`. Nextflow path already passes `task.cpus` (safe). **Effort:** XS.

## LO-15 — `usable_count >= 3` magic number duplicated
`pairhmm.rs:341, 455`, `alignment.rs:242` — three independent `< 3` literals, no shared constant or rationale. **Fix:** single named `MIN_USABLE_BASES` with a documented basis, shared across SW/PairHMM and applied to WFA per CR-2. **Effort:** XS.

---

# Cross-cutting recommendations

1. **Make the Rust engine output-aware.** Several findings (PF-1, ME-1, HI-7's silent no-op) share one root: the engine doesn't know what the user requested. Thread `mfsd: bool`, output format, and `gene_strand` across the FFI boundary. One structural change closes a whole class of compute-then-discard and silent-no-op bugs.
2. **Add a binned-vs-legacy parity gate to CI.** A parity test over a panel that includes large deletions, complex DelIns, and clustered variants would have caught CR-1 and CR-4 directly. Keep legacy `count_bam` (LO-2) solely as the parity oracle until this exists.
3. **One quality contract across alignment backends.** CR-2, HI-4, HI-5, LO-15 all stem from the fast path and fallback not sharing base-quality semantics. Centralize the BQ gate and apply it uniformly to WFA/SW/PairHMM.
4. **Statistical re-baseline.** CR-5, HI-10, HI-11, ME-8, ME-9 change reported p-values/LLRs/labels. Plan a single coordinated re-baseline of golden mFSD/ASJD outputs and notify report consumers.

# M5 — empirical scoping (measured on real data, 2026-06-30)

Before committing M5 effort, the tickets were measured on representative MSK data — a
deep duplex cfDNA panel (608 variants, ~16M reads, b37) and a real RNA-seq sample
(GRCh38 + Ensembl GTF) — at 1–8 threads, with temporary per-bin `fetch_ms`/`classify_ms`
instrumentation in `count_bin_shared`. **The measurements substantially re-prioritized the
milestone, mostly downward.**

### Single-sample findings
- **cfDNA counting is fetch-bound (89% fetch / 11% classify)** but already **~84%
  parallel-efficient** with **no long-pole** (busiest of 201 bins = 3% of bin-work; top-5
  = 13%). Deep bins read 150k+ reads each. t=1 → t=8 scales 3.65×.
- **Per-variant classification is ~0.5ms** (5× more variants over the same regions added
  only +1.3s); 94% SNVs, so alignment-backend cost is negligible.
- **RNA counting bin-work is trivial (~40ms total)**; the RNA run is dominated by the **~8s
  GTF parse**. Its one deep bin is *classify/ASJD*-bound, not fetch-bound.

### The multi-sample reframing (load-bearing)
**Core gbcms processes ONE sample per invocation; Nextflow provides multi-sample
parallelism** (`GBCMS_DNA|GBCMS_RNA` = 4 cores/sample, up to 100 concurrent jobs, GTF
passed per task). So the cohort runs as N concurrent 4-core processes:
- The GTF is **re-parsed every sample** → N × ~8s of identical work at cohort scale.
- Per-process **memory** (mFSD arrays) and **thread count** multiply by the concurrent
  process count — the node, not the single run, is the constraint.

### Revised per-ticket verdict → outcome
| Ticket | Measured verdict | Outcome |
|--------|------------------|---------|
| **LO-14** thread clamp | **Top priority** — recast as a hard `--threads` budget (rayon + any decode threads ≤ `--threads`) so gbcms never oversubscribes under fan-out. | ✅ **Done (#53)** |
| **PF-1** mFSD gate | **Keep — for memory, not CPU.** Per-process RSS ×N concurrent processes is the real win; also the "output-aware" enabler. | ✅ **Done (#54)** |
| **M5a** GTF index cache | **New — the biggest cohort lever.** Parse the GTF once, cache the parsed intermediate, reuse across the cohort (~9s → ~0.05s/sample). | ✅ **Done (#55)** |
| **PF-3** htslib decode threads | **Drop / liability.** Fetch already parallelized (84% eff.); decode threads oversubscribe the node under fan-out. | ❌ Dropped |
| **ME-13** overlap re-fetch | **Drop.** 0.00% re-fetch measured; highest parity risk in M5. | ❌ Dropped |
| **PF-4** pool reuse | **Moot.** 1 sample/process; only the MNP-rescue 2nd call rebuilds. | ❌ Dropped |
| **PF-2** sort long-pole | **Marginal.** cfDNA has no long-pole; RNA skew is trivial in absolute time. | ⏸ Niche flag (deferred) |
| **LO-3** bin-span cap | **Document, don't cap.** Chaining bounded (≤2.5× window); a cap risks re-breaking CR-1. | 📝 Document only |

### New investigation the data surfaced
- **M5b — Deep-bin fetch reduction (not started).** Deep cfDNA bins read 150k+ reads to
  count a handful of variants; the one-fetch-per-bin design is the cost. Narrowing it is
  the only remaining cfDNA lever, but it is parity-sensitive — scope carefully behind the
  binned↔legacy parity gate.

**Bottom line.** M5 as originally scoped optimizes within-sample counting that is already
efficient and is *not* the cohort bottleneck. The three items with real cohort ROI —
**LO-14 (budget) + PF-1 (memory) + M5a (GTF cache) — are shipped (#53 / #54 / #55).** What
remains is optional: **PF-2** (niche cost-sort), **LO-3** (a doc note), and the **M5b**
investigation. Dropped: **ME-13, PF-3, PF-4**.

# Suggested milestones

| Milestone | Tickets | Theme |
|-----------|---------|-------|
| **M1 — Count correctness (release-blocking)** | CR-1, CR-3, CR-2, CR-4, HI-2, HI-6 | Stop silent miscounts & panics |
| **M2 — Statistical integrity** | CR-5, HI-10, HI-11, ME-8, ME-9, ME-10, ME-11 | Trustworthy clinical labels |
| **M3 — RNA mode correctness** | HI-7, HI-8, HI-9, ME-6, ME-7, LO-9, LO-10, LO-11, LO-12 | Make RNA features actually work |
| **M4 — Alignment robustness** | HI-3, HI-4, HI-5, ME-4, ME-5, LO-15 | Backend parity & numerics |
| **M5 — Performance/IO** | ✅ LO-14 (#53) + PF-1 (#54) + M5a GTF cache (#55); ❌ dropped ME-13/PF-3/PF-4; ⏸ PF-2 niche; 📝 LO-3 doc; 🔍 M5b open | **Re-scoped on real-data measurement** — see §"M5 — empirical scoping" |
| **M6 — Hygiene & contracts** | HI-1, ME-1, ME-2, ME-12, LO-1, LO-2, LO-4, LO-5, LO-6, LO-7, LO-8, LO-13, DX-1 | Exit codes, output parity, docs |

Recommended order: **M1 → M2 → M3** are correctness; do them first and in that order. **M4** can parallelize with M2/M3. **M5** is independent (pure perf). **M6** is continuous cleanup.

---

# Deferred ledger — one-stop index (as of 2026-06-30)

Everything consciously **not done yet**, in one place so nothing is silently lost. Details
live in the linked ticket sections above; this is the index.

### Dropped after real-data measurement (revive only with new evidence)
- **ME-13** overlap re-fetch — 0.00% re-fetch measured; highest parity risk in M5.
- **PF-3** htslib decode threads — oversubscribes the node under Nextflow fan-out.
- **PF-4** rayon pool reuse — moot at 1 sample/process.

### Deferred / optional (M5)
- **PF-2** bin cost-sort — niche (no cfDNA long-pole); cheap if a skewed workload appears.
- **LO-3** bin-span cap — **document only**, never cap (a cap risks re-breaking CR-1).
- **M5b** deep-bin fetch reduction — *not started*; the only remaining cfDNA lever, parity-sensitive.

### Accepted deviations — correctness shipped, refinement deferred (see §"Accepted deviations")
- **CR-5** — report **mean** LLR per fragment vs raw sum (display nicety; sign/labels unaffected).
- **ME-9** — raise `MIN_FOR_KS` 5→≥8 (largely moot after CR-5's exact small-N KS).

### M6 backlog — deferred to the next milestone
- **HIGH:** **HI-1** exit-code propagation (all-samples-fail still returns 0).
- **MED:** ME-1 (sub/mono-nuc computed but dropped from VCF), ME-2 (`;` vs documented `|` delimiter), ME-12 ("majority rule" is structural-ALT-wins — relabel).
- **LOW:** LO-1 (dedupe `.pyi` stubs), LO-4, LO-5, LO-6, LO-7, LO-8, LO-13 (docs / minor numerics).
- *(Already done in M6: LO-2 feature-gate, DX-1 label sweep.)*

### Test / CI gaps (surfaced by the 2026-06-30 coverage review; fold into M6)
- **⚠ CI ran no `cargo test` and no `cargo clippy`.** The ~210 Rust `#[cfg(test)]` unit tests
  (mFSD math, GTF-cache roundtrip, thread-budget guard, alignment/CIGAR, WFA, contig
  normalization) and clippy never ran in CI, so a Rust-logic regression could ship green. Note
  the binned↔legacy **parity** assertions are *Python-level* (`helpers.count_both`) and **do**
  already run under pytest against the default-features wheel — the gap was the Rust side.
  **Fixed by #58:** a `rust-test` job (`cargo test`, default *and* `--no-default-features`, on
  macOS since the pyo3 `extension-module` crate can't link a test binary on Linux) + `cargo
  clippy -D warnings` (both feature configs) added to the `lint` job.
- **CRAM** counting has zero tests (only BAM fixtures).
- **HI-1** needs a test that an all-samples-fail run returns non-zero (couples with the fix).
- **`--threads N`** has a reject-`0` test but no positive cap-enforcement test.
- Real-cohort integration (mFSD / strandedness / ASJD / GTF cache) is validated only
  out-of-repo on MSK data; no in-repo cohort fixture.

### Known issue — not yet a ticket
- **DNA-mode Rust INFO logs not always captured to stderr** (pyo3-log routing). Cosmetic; M6 hygiene.

---

# Accepted deviations (as-built through M1–M4)

The shipped M1–M4 work departs from the plan's letter in three places. Each is a
deliberate, documented choice: the safety-critical intent of every ticket is met,
and what remains deferred is a refinement that changes neither the counts nor the
clinical labels. Recorded here so the plan reflects what actually shipped.

## CR-5 — mean LLR per fragment (deferred refinement)
**Done.** The two critical defects are fixed. The LLR is computed in closed-form
log-space (`mfsd.rs:130-146`: `0.5·(z_h²−z_t²) + ln(σ_h/σ_t)`), finite for all finite
inputs — a single tail fragment can no longer pin it to ±∞. Small-N KS uses the exact
two-sample null (`mfsd.rs:278-325`, Hodges lattice-path) below a size threshold, with
the asymptotic series only where it is accurate.
**Deferred.** Plan item 2 — *reporting the mean LLR per fragment* rather than the raw
sum (`calc_llr` still sums, `mfsd.rs:111-119`). This is a cross-variant comparability
nicety, not a correctness issue: the sign of the statistic and the per-variant label
are unaffected.
**Promotion target.** If adopted, divide by `lengths.len()` in `calc_llr_with_params`
and re-baseline golden LLR magnitudes (it shifts displayed values — coordinate with
report consumers).

## ME-7 — QNAME dedup of junction reads (implemented)
**Done.** Both halves shipped.
1. *Transcript-strand fold* — FR mates splitting a genuine junction across both genomic
   strands (the false `STRAND_DISCORDANT` bug) is fixed by folding to inferred
   transcript strand before tallying (`read_transcript_strand`), now selectable via
   `--strandedness`.
2. *Per-fragment dedup* — ASJD junction evidence is now counted once per fragment, not
   per mate (`JunctionTally` in `detect_asjd`). This was driven by measuring the real
   reverse-stranded FORTE BAM (3.6M reads, 117,862 junctions):
   - Mate-overlap at junctions is **common, not rare**: 35.6% of fragment×junction
     incidences (the earlier "uncommon" guess was wrong). Mechanism: a molecule's R1
     and R2 overlap on a short cDNA insert and both carry the same `N` op for a shared
     (often large) intron — not a UMI/PCR duplicate (UMI dedup leaves 0 flagged dups,
     2 records/QNAME). It inflates junction totals ~1.38× genome-wide.
   - Mates always fold to the same transcript strand (0/319k disagreements), so the
     dedup is unambiguous (first-seen wins).
   - Per-variant impact is small (each variant's ASJD window sees few reads): 18 of
     58,240 junctions flip `STRAND_DISCORDANT` (0.031%), all low-count and all
     corrective. The fix is principled (per-fragment matches the rest of the engine)
     rather than high-magnitude.
   Probes: `scratchpad/me7_junction_dedup_validation.py`, `…dedup_design_probe.py`.

## ME-9 — raise MIN_FOR_KS / report min_alt (now largely moot)
**Done.** The safety goal — never classifying a variant on a KS test that did not run —
is met. `_classify_origin` returns `INSUFFICIENT` when `ks_valid` is false
(`mfsd_report.py:111-112`), and the BH family excludes NaN-D variants (HI-11 / ME-8).
**Deferred (and largely moot).** Raising `MIN_FOR_KS` from 5 to ≥8 and forcing
`min_alt ≥ MIN_FOR_KS`. The original concern was that the *asymptotic* KS has no
resolution at n=5; CR-5's exact small-N KS removes it — at n=5 the p-value is now the
exact null, not an out-of-range approximation. With exact KS plus the `ks_valid`
INSUFFICIENT gate, the floor of 5 is defensible and report `min_alt=3` only governs
what is *displayed* (always gated to INSUFFICIENT when the test could not run).
**Promotion target.** If a power analysis later justifies a higher floor, bump
`MIN_FOR_KS` (`mfsd.rs:33`) and match the report `min_alt` default; document the
minimum detectable D at the chosen n.

---

# DX — Developer experience / maintainability

## DX-1 — Code comments and logs cite ephemeral ticket labels

**Problem.** Comments and log strings across the codebase carry ticket/milestone
labels — pre-existing `P4c` markers in the RNA/ASJD engine paths, and `CR-*`/`HI-*`/
`ME-*` labels added during the M1–M2 remediation. A label like `// ME-8: padding fix`
or `debug!("P4c BH-FDR: …")` is meaningless to a contributor with no plan/PR context;
it documents *which ticket touched the line*, not *what the code does or why*.

**Fix.** Sweep code comments and log strings; replace each label with the underlying
reason (e.g. `// ME-8: padding fix` → `// exclude no-test variants; they inflate the
FDR family`). Ticket context stays in the commit message (`git blame`). Scope: source
only — this plan and other design docs are exactly where labels belong, so leave them.
Start points: `rust/src/counting/engine.rs` (`P4a`/`P4b`/`P4c`, `CR-1`, `HI-11`),
`rust/src/counting/{mfsd,variant_checks,wfa_router}.rs`, `src/gbcms/io/output.py`,
`src/gbcms/report/mfsd_report.py`, the `.pyi` stubs. `grep -rnE '(CR|HI|ME|LO|PF)-[0-9]|P4[abc]'`
over `src/` and `rust/src/` enumerates the sites.

**Prevention (done).** Promoted to a standing rule — `.agents/rules/code-quality.md`
§"Comment & Log Hygiene" — so new code doesn't reintroduce labels. Logged as
LRN-20260627-001.

**Severity:** LOW (no behavior change). **Milestone:** M6.

**Status: Resolved** (PR #43). Source comments and logs swept across the 10 flagged
files; `grep -rnE '(CR|HI|ME|LO|PF)-[0-9]|P4[abc]'` over `src/` and `rust/src/` returns
zero. Design docs (this plan included) intentionally retain labels.
