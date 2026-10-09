# CONTINUITY — where we left off

> Tactical state that must survive a closed laptop or a context summary.
> Update the **Now** and **Next** sections as work progresses.

_Last updated: 2026-10-08_

## Now
**6.6.0 is released** (2026-10-08): #252 merged to main (merge commit 67c45989, tree
identical to the gated candidate 09612675), annotated tag `6.6.0`, main back-merged
into develop, develop at `6.7.0-dev.0`. The release workflow publishes PyPI, GHCR, the
docs and the GitHub Release page (now automatic). The docs deploy on main failed at
6.6.0 (the stable step read pyproject's version, dynamic since #229); #256 reads it
from Cargo.toml via `release.py version` and deployed the 6.6.0 docs.

The gate (#155) measured from **6.3.1**, the last release compared on HPC; 6.4.0 and
6.5.0 were untested midpoints, so their merges were checkpoints inside one interval
(operator decision, 2026-10-07). Results are on #252 (aggregates). Local record:
`~/test/gbcms/harness/regression_panel/gate_660/ADJUDICATION.md` (PHI); the HPC kit is
`kit_660/` (local copy beside it; on HPC under the operator's test folder).
- 6.6.0 411/411 runs; 6.3.1 407 (4 mfsd crash on a rejected variant, fixed in 6.4.0 #95).
- 73,649 changed cells attributed over 40 checkpoints, none unattributed, reduction exact.
- Concordance with sign-out 98.4% -> 98.1%; the 19 strata that fell were adjudicated
  read by read: no 6.6.0 defect. Peak memory per task 4.74 -> 0.66 GB (median).
- The origin events (AR 9 bp GGC-repeat deletion, FLT3 2 bp A-run deletion, IMPACT BAM
  SW*211-T): 6.6.0 counts exactly the read census (AR 17/15, FLT3 183/7 REF/ALT); 6.3.1
  gave 167/24 and 249/259 (it credited FLT3's 1 bp deletion carriers). Local record:
  `~/test/gbcms/harness/origin_events/README.md`.
- Lessons kept in the kit's RUNBOOK: run scripts from a Python 3.11 venv (system python3
  is older); the manylinux_2_34 wheels run inside `python:3.11-slim-bookworm` on RHEL 8
  (glibc 2.28); index FORTE's genome.fa in the kit; set aside failed runs' partial MAFs
  before `attribute.py prepare` (#255); peak memory from `sacct` MaxRSS (no GNU time).

Working rules (carried into 6.7.0): measure first, red-first tests, an adversarial
review, real-data acceptance per read on a matrix of data types (gbcms is a
genotyper, not a caller), mount runs one at a time, and community practice as a
floor, not a ceiling. The mskcc-omics-workflows containers and gbcmsrs modules are not
a release step: the operator starts that switch when satisfied it is production-ready.

## Next → 6.7.0
Tracker #196 and the 6.7.0 milestone hold the work. From the 6.6.0 gate:
- #253 long delins under the discrimination window (spanning 20+ bp: 59% -> 6%
  agreement with sign-out; short-anchored exact carriers withdrawn by design);
- #254 clipped exact-junction reads of deletions longer than a read;
- #255 `attribute.py prepare` skips runs whose latest attempt failed.
Also deferred: per-locus gap parameters (#244), the sibling shapes measured for #246,
the panel's GIAB and TEMPO arms, and the `build-gtf-cache` removal.

**Plan drafted (2026-10-08), awaiting the operator's decisions 1–8:** `CYCLE_6.7.0_PLAN.md`
on `feature/6.7.0-plan` (not pushed). Measured on the gate panel's delins and 50+ bp
deletion rows, a 193-row germline-het pilot in IMPACT normals, and four surveys; local
harness `~/test/gbcms/harness/plan670/` (PHI). Headlines: one-sided equal-length junction
windows keep VAF unbiased and recover the long-delins carriers (AD +45%, sign-out 0/14 →
9/14); past 50 bases 6.6.0's written counts already match one-sided reading; long
deletions read 4–7% low (clipped ALT, REF needs 1–2 bases); allele-level VAF holds 0.5
at germline hets in repeats, by-base 3–6 points low; Phase 3 decides ≤0.3% of reads
(#244: close). Plan approved 2026-10-08; implementing (D8 #270, H4 #271, C39 next).

### Previous: 6.6.0 release
`release/6.6.0` was cut from develop c01749ab on 2026-10-07: version bump, the cycle
plan removed (private copy in `~/test/gbcms/plans/`), the Docker lock refreshed
(polars 2.0.0). The commit guard learned to admit the checkpoint list (09612675).
Contents: the CHANGELOG's 6.6.0 section (tracker #140).

### Previous: 6.5.0 release
**6.5.0 is released** (2026-09-25): #118 merged to main, bare tag `6.5.0`, and
main back-merged into develop. The release workflow publishes to PyPI, GHCR and
the docs; the GitHub Release page is created by hand (release guide step 9).
Cycle plans live on develop (PHI-free) and are removed on the release branch at
each cut, so they never ship. `CYCLE_6.5.0_PLAN.md` was removed at this cut; a
private copy is kept locally in `~/test/gbcms/plans/`.

6.5.0 contents:
- T1 cluster exclusive assignment (#99).
- T2 ASJD-2 markers (#100).
- T3–T5 observability (#102).
- T7 MNP rescue (#101).
- T8 contig naming + merge (#104).
- T6 one BAQ rule across RNA views (#105).
- T10 decomposed twin strand + per-sample flag (#113).
- T12 VCF ↔ MAF representation per vcf2maf / maf2vcf, plus `gbcms convert`
  (#116, closes #110). A breaking output change: see the CHANGELOG.
  - #116 also fixed what verifying the plan against the code found: T1's
    grouping let an SNV between two window-joined deletions join their cluster;
    merge by VCF record dropped later-only rows' coordinates; a MAF deletion at
    Start 1 wrapped u64 in debug builds.

RC check (develop vs the 6.4.0 release, every changed cell attributed; harness
`~/test/gbcms/harness/rc650/`, local only):
- DNA, 28 runs (FLT3 IMPACT, ACCESS duplex+simplex incl. BRCA2, complex-cluster,
  MSI-high): 1030/1060 rows byte-identical; 30 changed = 29 T1 cluster rows + 1
  T3 flag; 0 unexplained. BRCA2 reproduces the T1 record.
- RNA, FORTE truth cohort (33 samples / 94 rows): zero main-count changes; the six
  T2 marker rows and one T6 exon-edge row; 0 unexplained.
- Re-checked with #116's code, and finally on the exact 6.5.0 code (an isolated
  build of develop 851119f3, `~/test/gbcms/harness/t110/real/final_rc.py`):
  byte-identical, DNA 28/28 runs and RNA 33/33.

**Then → 6.6.0.** The plan was `CYCLE_6.6.0_PLAN.md` (#117; removed at the 6.6.0
cut, a private copy in `~/test/gbcms/plans/`). It had 35 tickets, each with a GitHub
issue in the **6.6.0 milestone**, all under the tracking issue **#140**. Grouped
work is filed as sub-issues: #92 → #141–#144, #112 → #145–#147, #133 → #148–#149,
#134 → #150–#152, #135 → #153–#154. The other top-level issues are #106, #111, #114,
#119–#139, #155 and #156. Start with the operator decisions (C2, R1, R2, I1, I3, I4,
S1/S2), then the count-affecting tickets, each measured first.

Release comparison: the HPC run was not done before the 6.5.0 tag. It became the
6.6.0 release gate, on the coverage-driven regression panel (D5, #155) instead
of the hand-picked 56 samples. The selection harness is local
(`~/test/gbcms/harness/regression_panel/`, PHI) and the operator carries its
output to HPC. If a comparison feeds VCF input and compares MAF output by
`Start_Position`, key it on `vcf_pos`/`vcf_ref`/`vcf_alt` instead. Optional: a
head-to-head against C++ GBCMS 1.2.4/1.2.5 (on HPC).

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

### Previous: M6 (hygiene & contracts)
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
- Binning invariance replaces the legacy `count_bam` parity oracle (T1 #170,
  2026-10-01): a second engine sharing the classifier never caught a classification
  bug; counts must not depend on bin geometry, and the read census checks
  classification.
- Tests kept minimal/high-signal (each is a maintenance contract); fixes that reduce
  duplication are preferred over adding code.
- Source comments/logs explain what/why/how — **never** ticket labels (`CR-`/`HI-`/
  `ME-`/`P4c`/…); that context lives in the commit message. Rule in
  `.agents/rules/code-quality.md`; existing labels swept under DX-1.
