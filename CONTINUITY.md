# CONTINUITY — where we left off

> Tactical state that must survive a closed laptop or a context summary.
> Update the **Now** and **Next** sections as work progresses.

_Last updated: 2026-09-24_

## Now
**6.5.0 cycle complete on develop; release-candidate validation PASSED
(2026-09-24).** Merged to develop: T1 cluster exclusive assignment (#99), T2
ASJD-2 markers (#100), T3–T5 observability (#102), T7 MNP rescue (#101), T8
contig naming + merge (#104), T6 one BAQ rule across RNA views (#105), T10
decomposed twin strand + per-sample flag (#113). 6.4.0 is released (#98, tag
`6.4.0`) and main is merged back into develop.

RC check (develop 3ed1a4d3 vs the 6.4.0 release, every changed cell attributed;
harness `~/test/gbcms/harness/rc650/`, local only, README there):
- DNA, 28 runs (FLT3 IMPACT, ACCESS duplex+simplex incl. BRCA2, complex-cluster,
  MSI-high): headers identical; 1030/1060 rows byte-identical; 30 changed = 29 T1
  cluster rows + 1 T3 flag; 0 unexplained. BRCA2 cluster reproduces the T1
  record: 28/28, 32/32, 21/21 exact; 33bp 37 v 22; 14bp 51 v 45.
- RNA, rebuilt FORTE truth cohort (33 samples / 94 rows, sequence-anchored hg38
  lift): zero main-count changes; the six T2 marker rows of the T2 record and one
  T6 exon-edge row; 0 unexplained.

The T1/T2 acceptance harnesses lost in the 2026-09-23 scratchpad wipe are rebuilt
there; per-ticket harnesses: `~/test/gbcms/harness/{t6,t7,t8,t9t10}/`.

**Next → release 6.5.0** (same procedure as 6.4.0): cut `release/6.5.0` from
develop — version bump in the 11 files the 6.4.0 cut touched (pyproject,
`src/gbcms/__init__.py`, `rust/Cargo.toml` + lock, `nextflow/nextflow.config`, the
five `nextflow/modules/local/gbcms/*/main.nf` container tags) + CHANGELOG cut →
PR → main, gated on: container build+push (operator) → the HPC 56-sample IMPACT
matrix vs the 6.4.0 baseline (script at `~/test/gbcms/dev_regression/` on HPC,
repointed at the rc container) + RNA smokes → tag `v6.5.0` → merge → back-merge
develop. Optional: head-to-head vs C++ GBCMS 1.2.4/1.2.5 (on HPC).

**After the cut:** T9 span-aware exon-edge BAQ rule (#106, low priority), T11
decomposition arbitration redesign (#111), #112 (consumers of a decomposed
winner), #110 (VCF delins with a one-base ALT treated as a pure deletion — decide
whether it rides 6.5.0), #114 (RNA strandedness gating: `rna_antisense_depth`
always 0 and `STRAND_DISCORDANT` unreachable at defaults — an open decision
carried from #94), and #92's three remaining items.

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
