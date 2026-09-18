# CONTINUITY — where we left off

> Tactical state that must survive a closed laptop or a context summary.
> Update the **Now** and **Next** sections as work progresses.

_Last updated: 2026-09-18_

## Now
**Issue #94 clusters A + B1 — DONE on `fix/rna-splice-aware-counting`
(battery 5bc8f47 red → 56e703b green → B1 commit; not yet pushed/PR'd).** The rule: a read testifies only through aligned bases (or a D
op) at the discriminating positions; N over all of them → neither + excluded
from DP/DPF (splice_skip_triage + covers_locus, mirrored binned/legacy/
per-transcript); D-vs-N never flips the call; post-N anchor/windowed
inspection via helpers shared with the M-arm; Phase 3 never scores across a
splice; FragmentEvidence::resolve counts structural ALT at qual 0 (BAQ
stacking made ad>0/adf=0 diverge); SPLICE_SKIP_DOMINANT(n) diagnostic for
STAR's N-represented large deletions. 19-agent adversarial review: 15
confirmed findings all fixed or deliberately accepted+documented. Validated:
b37 dedup RNA (DP drop pysam-exact) and FORTE hg38 GAPDH junction smoke
(intronic DP 826→2; use `~/Downloads/pipeline_resources/gunzip_gtf/
Homo_sapiens.GRCh38.111.gtf` for FORTE — Ensembl contigs, found 2026-09-18).

**#94 cluster B1 — DONE (same branch, uncommitted→committed today):**
consensus splicing of ref_context REMOVED (no coordinate map = corruption;
its Phase-3 consumer unreachable; exon-contained reads mis-scored).
ref_context is always genomic; the xfail flipped green; mq0_count now
tallied before the strandedness filter in binned (matches legacy); two new
regression pins (≥50bp band guard near junction; Phase-3 pangenome
variant+sibling matrix DNA↔RNA mode-equivalence). Sonnet-model adversarial
review (Fable subagents hit the monthly spend cap — model override in the
workflow script is the mechanism): all confirmed findings fixed. Real data:
b37 + FORTE counts unchanged.

**Next: #94 B2 (evidence-gated) + remaining items** — measurement on FORTE
GAPDH acceptor deletion: 824 junction reads end `neither` via the splice
guards, and they observe the deleted span WITH ALIGNED BASES → the dominant
recoverable population needs only a STRUCTURAL rule (span-aligned REF
testimony at deletion loci whose anchor is spliced out), NOT the
coordinate-mapped spliced-haplotype machinery (B2b, delins-carrier tail
only; must keep pre-mRNA/intron-retention reads genomically scored). Then:
strandedness gating decision remains open (order now consistent for mq0);
RNA-BAQ measurement (junction-adjacent inserted bases are BAQ-unverifiable
without --gtf suppression — pinned in the insertion contract test);
per-transcript partial_alt is locus-level-only by design (no new columns).
PR after B2a decision.

**Issue #91 — wrong-length pure-indel fix: MERGED to develop (#93, 75062d6);
issue closed.** Ships with the next minor release — the release branch cuts
**6.4.0** per the CHANGELOG [Unreleased] callout. Three-regime taxonomy
validated on internal data (local pointers: gitignored
`user-local-validation-data` memory; baselines:
`~/Downloads/gbcms_test/wrong_length_baselines/` — keep for #92 validation).

**Next up:** class-1 xfail battery (#92) -> **#94 RNA splice-aware correctness**
(20 confirmed findings, 6 high: spliced-around reads count definitive REF for
junction-abutting indels; D6 consensus-splicing corrupts ref_context offsets;
insert_truncation_match lacks the reliable-base gate [our Phase-2, invariant-2
violation]; RNA-BAQ vs unverifiable-bases branch). Holistic requirements now
standing for every item: RNA lens on BOTH builds (b37 dedup RNA + v47lift37;
hg38 FORTE STAR + v50 — all local, verified), real-data validation (ACCESS
duplex + consensus ladders local; MSI-H + FLT3-ITD being pulled), community
head-to-head harness (C++ GBCMS, bam-readcount, GATK) scored vs pysam truth.
Identity-band design intent: annotation-representation tolerance — test both
policies (see memory + #92 comments). No new output columns, ever.
**Issue #92**: class-1 DONE — PR #95 open (battery a4aa891 red → e23c4c3
green; 11 review findings fixed incl. ASJD-sentinel high; real ACCESS/RNA
smoke). Next after #95 merges: #94 RNA splice-aware correctness, then
remaining class-2 (insert band both-policies + in-tract sub-5bp gate with
offline mismatch-anatomy survey first), class-3 observability.
(Original class-1 note: the xfail-strict battery (encode every deterministic repro as a
failing test first: _zero_counts crash, --bam-list fail-fast, OSError swallow,
duplicate samples, MafWriter collisions, merge coercion, ragged rows, traceback
logging), then fix red→green on one branch. Class 2 (measured, count-affecting:
error-tolerant exact-length insert band; inert dynamic_sw_gap_extend decision)
and class 3 (UMI warn, per-transcript BAQ suppression question) each get their
own branch + validation round.

- **Phase 0 DONE:** target-contract battery `tests/test_wrong_length_contract.py`;
  baselines captured; MAF/VCF input equivalence confirmed on the bug case.
- **Phase 1 DONE** (760a1ea): repeat scan anchors at the first changed base;
  padding `span + default_pad`.
- **Phase 2 DONE** (586218a + 6a52414): wrong-length rule — placement-aware
  ≥50bp deletion band (≤3 retained in span, ≤3 changed outside; covers split
  reps, rejects displaced net-matches), insertion truncation containment
  (dual low-complexity gates), windowed flags split (same-length S3-fail keeps
  Phase-3; wrong-length → repeat: neither+partial / unique: REF+partial), NO
  Phase-3 for wrong-length (length-blind + tiny-context promotion risk, proven
  by adversarial review). Plus: `--trace` was entirely dead (pyo3-log Debug cap
  + wrong logger name "gbcms_rs" vs "_rs") — fixed, with `_rs.reset_log_caching()`.
  Contract 19/19; pytest 448; all acceptance targets exact (FLT3 7/252,
  AR 15/9, RNA 15/17, SMARCB1 310); 74-set + fullbam5 panels explained-identical.
  Phase-2 detail comment on #91.
- **Phase 3 DONE** (a96b8a6): 36-agent silent-failure audit → 21 confirmed
  findings. Fixed here: same-length wrong-seq insertion no longer absorbed
  into rd (third allele → neither+partial; Phase-3 only for unverifiable
  bases — SW provably promotes confident wrong-seq to ALT); phase3 "SW
  fallback" logs corrected; enforce-strandedness + GTF-coverage warns (incl.
  cache-hit path); left-align failure loudness (fetch fail, cap bind, UTF-8);
  per-read debug→trace; dynamic_sw_gap_extend divergence doc; PARTIAL_DOMINANT
  e2e reporting-chain test. Measured: ad unchanged everywhere; 8 reads
  rd→partial on 3 insertion rows. Separable findings (2 high: _zero_counts
  stub crash, --bam-list fail-fast violation) → issue #92.
- **Phase 4 DONE** (f261d9b): docs rewritten + 24-agent docs-vs-code
  verification (21 discrepancies fixed, much pre-existing rot: MNP no-fallback,
  Phase-2.5 routing, interior guard, fictional RUST_LOG/GBCMS_LOG_LEVEL, false
  DP invariant, broken normalize command, invalid windowed examples). New code
  finding → #92: dynamic_sw_gap_extend is a constant −1 (logistic caps at 0.5;
  relaxation never engages). CHANGELOG [Unreleased] with minor-bump callout;
  counting-engine skill updated; MAF↔VCF equivalence test committed
  (tests/test_e2e_maf_vcf_equivalence.py).
- **Phase 5 DONE:** final-build validation matrix all-exact (FLT3 7/252,
  AR 15/9, RNA 15/17; 74-set 69/74 identical with exactly the 5 predicted
  rows; fullbam5 14/15 with the 1 predicted row; CLI matrix + large panel
  conform, 2 script rows attributed to un-pinned-anchor construction).
  Review reconciliation: 69 findings across 4 adversarial passes, all fixed
  or in #92. Branch PHI scan clean. Version: release branch cuts 6.4.0
  (CHANGELOG callout shipped). **PR #93 open to develop** (PHI/accuracy
  verified body). Remaining: CI green -> merge; then #92 follow-ups
  (class-1 xfail battery first).

Previous state (code-review remediation) below for history.

### Previous: code-review remediation
Working the code-review remediation plan (`CODE_REVIEW_IMPLEMENTATION_PLAN.md`)
ticket by ticket, one PR each, with review-before/after discipline and real-data
validation on MSK cfDNA-duplex (b37) and STAR/FORTE RNA (GRCh38) samples.

**M1 (count correctness) — COMPLETE & merged:** CR-1, CR-3, CR-2, CR-4 (+#29 Docker
CI). pysam-cross-checked (ALT 571/571 exact).

**M2 (statistical integrity) — COMPLETE & merged:** CR-5 (closed-form LLR, exact
small-N KS, no scipy), HI-10, HI-11 (BH-FDR `mfsd_qval_alt_ref`), ME-8 (ASJD BH
padding), ME-9 (classify only on a valid KS test; `MIN_FOR_KS` stays 5), ME-10 (NaN
enrichment disambiguation), ME-11 (isinf guard).

**M3 (RNA correctness) — COMPLETE & merged (#34, #35, #36):** LO-12 (NH integer
widths), HI-8 (intron-discounted fragment size), LO-9 (GTF diagnostics), ME-6
(`normalize_contig` M/MT), HI-7 (populate `gene_strand` — the keystone), LO-10
(unstranded→None), HI-9 (strand-aware splice motifs), ME-7 (dUTP-folded
STRAND_DISCORDANT), LO-11 (base-aware editing). Validated end-to-end on a real
reverse-stranded RNA sample; dUTP sense/antisense split matched samtools 2/2.

**M4 (alignment robustness) — COMPLETE & merged:** LO-15 (shared `MIN_USABLE_BASES`),
HI-3 (WFA off-target global edit distance + length-aware threshold), HI-4/HI-5/ME-5
(PairHMM numerics), ME-4 (haplotype-parity collision: `warn!` + `NON_DISCRIMINATING_LOCUS`
flag, #40). Defensive — byte-identical counts on real cfDNA.

**M1–M4 audit + follow-ups — COMPLETE & merged:** cross-checked the plan against
the code (4-agent sweep). Closed the real gaps — HI-2 + HI-6 (#41), MIN_FOR_KS coupling
+ dead tooling (#42), DX-1 label sweep (#43), test/doc gaps + BH-family DRY (#44),
prep-time empty-allele rejection (#45), large-deletion binned↔legacy parity gate (#46),
as-built plan-deviation notes (#47).

**ME-7 fully resolved (RNA strandedness) — merged:** dynamic `--strandedness`
(reverse/forward/unstranded, #48) and per-fragment ASJD junction dedup (#50) — both
data-driven from the real FORTE BAM. ME-7 is now *implemented*, not a deviation.

**LO-2:** the legacy `count_bam` parity oracle is feature-gated (`legacy-parity`, default
on) out of the shipped wheel.

**M5 — mostly shipped, re-scoped on real-data measurement.** The tickets were measured on
a real deep cfDNA panel + RNA sample (counting is already ~84% parallel-efficient; the
cohort runs as N concurrent 4-core Nextflow processes, so the *node* — memory, GTF
re-parse ×N — is the constraint, not the single run). The three items with real cohort ROI
are merged: **LO-14** (hard `--threads` budget, #53), **PF-1** (output-aware mFSD gating
for per-process memory, #54), and **M5a** (GTF index cache + `build-gtf-cache` pre-warm +
automatic Nextflow wiring; ~9s → ~0.05s per sample, #55). The rest were **dropped as
low-ROI/liabilities**: ME-13 (0% re-fetch measured), PF-3 (decode threads oversubscribe
under fan-out), PF-4 (moot at 1 sample/process). See the plan's §"M5 — empirical scoping".

## Next
**M6 (Hygiene & contracts) — COMPLETE.** All M6 tickets landed: HI-1 (exit non-zero on
sample failure; empty variant set still exits 0), ME-1 (sub/mono-nuc reach VCF), ME-2
(`|` transcript delimiter), ME-12 (fragment-consensus relabel; #64), LO-1 (single `_rs`
stub + `count_bam` signature; #65), LO-2 + DX-1 (earlier), and the LO-4/5/6/7/8/13
documentation sweep (#65). CI now runs `cargo test` + clippy (both feature configs).
Out of band this session: contig-mismatch reconcile + warn (silent zero-count bug, #63)
and a synthetic DNA→MAF end-to-end test (#62).

**All six milestones (M1–M6) are now complete and merged**, except the explicitly
optional/deferred M5 leftovers below.

Remaining (all optional):
1. **M5 leftovers:** **PF-2** (bin cost-sort — niche; cfDNA has no long-pole),
   **LO-3** (document the bin-span soft floor, don't cap — doc only), and **M5b** (deep-bin
   fetch reduction — the only remaining cfDNA lever, parity-sensitive; not started).
2. **RNA LO polish not in M6 scope:** LO-9/10/11/12 landed in M3; no open RNA LO items.

## Open follow-ups (tracked, not lost)
- ✅ BAM-level binned↔legacy parity gate, incl. large deletions — **done (#46)**.
- ✅ Prep-time empty-allele validation (loud, once-per-variant) — **done (#45)**.
- ✅ DX-1 source label sweep — **done (#43)**; standing rule prevents reintroduction.
- ✅ LO-2: feature-gate the legacy `count_bam` parity oracle out of the shipped wheel
  — **done** (default `legacy-parity` feature; release builds `--no-default-features`).
- Accepted plan deviations (CR-5 mean-LLR, ME-9 thresholds) are documented in
  `CODE_REVIEW_IMPLEMENTATION_PLAN.md` § "Accepted deviations" with promotion targets —
  pick up only if a consumer needs them. (ME-7 is now implemented, not a deviation.)

## Key decisions
- Stats stay in **Rust** (KS/LLR/Fisher in `mfsd.rs`/`shared/stats.rs`); **no scipy**
  dependency — exact KS is a self-contained Rust DP, validated against baked
  SciPy reference constants.
- Legacy `count_bam` (per-variant) is kept as the binned↔legacy **parity oracle** but
  feature-gated (`legacy-parity`, default on) so the **shipped wheel excludes it**
  (release builds `--no-default-features`); production uses `count_bam_binned` only.
- Tests kept minimal/high-signal (each is a maintenance contract); fixes that reduce
  duplication are preferred over adding code.
- Source comments/logs explain what/why/how — **never** ticket labels (`CR-`/`HI-`/
  `ME-`/`P4c`/…); that context lives in the commit message. Rule in
  `.agents/rules/code-quality.md`; existing labels swept under DX-1.
