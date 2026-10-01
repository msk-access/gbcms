# CONTINUITY — where we left off

> Tactical state that must survive a closed laptop or a context summary.
> Update the **Now** and **Next** sections as work progresses.

_Last updated: 2026-10-01_

## Now
**6.6.0 cycle in progress** (plan: `CYCLE_6.6.0_PLAN.md` on develop, tracker #140,
every item a sub-issue). Merged to develop so far:
- #160 C10 + C2: REF only from informative reads; window-aware sibling REF guard.
- #161 count-the-given-allele principle (AGENTS.md invariant 7) and the plan.
- #162 named per-read trace lines (read-level validation).
- #164 O4: `OBSERVED_ALLELE` / `COEXISTING_ALLELE`; homopolymer twin opt-in.
- #165 C1: complex variants count exact carriers (anchored windows, repeat
  growth, equal-length windows, junctions for long events, splice-aware windows,
  multi-allelic exact ties).
- #168 C12: carriers whose allele lies in soft-clipped bases (complex/MNP, DNA).
- #169 C13 (#166): BAQ spares the variant's own indel evidence.
- #175 pending list, validation standard, the genotyper principle, and the
  BAM-caveats reference (`docs/reference/bam-evidence-caveats.md`).
- #181 R1 (#106): the exon-edge BAQ rule and `exon_boundary_dist` measure from
  the REF span; `--rescue-mnp` components are counted over their MNP's span.
  FORTE truth cohort 0 count changes; probes adjudicated per read (plan R1 entry;
  harness `~/test/gbcms/harness/r1/`, local).
- #184 C14 (#172): records without bases (SEQ `*`) are dropped in the shared read
  filter (they panicked, or counted REF from the CIGAR); none in MSK data.
  Follow-ups C19 #182 and O7 #183.
- #187 R2 (#114): rna_antisense_depth counts antisense reads under strandedness
  enforcement (classified, tallied, dropped; REF/ALT unchanged); STRAND_DISCORDANT
  documented as a --no-strandedness diagnostic. Follow-ups R4 #185, O8 #186.
- #190 C4 (#142): reference windows clamp only when they pass a contig end (by
  the record read); RC set byte-identical. Follow-up C20 #188.

- #193 C21 (#189): the windowed indel checks accept a placement by its haplotype
  when it is the read's only change across the discrimination window; the scan
  and tract-cluster grouping reach every placement in the shift region. Two
  adversarial reviews; acceptance adjudicated per read (harness
  `~/test/gbcms/harness/c21/`, local).

**Merged — cluster 1, #203 (C20 #188, C22 #191, C23 #192, C8 #121; C11 #159
closed with its measurement):** pure-indel reads count ALT only where their own bases hold
the ALT; strict-path reads with another indel in the window, same-length deletions
whose bases spell another allele, and anchor-keeping reads of `A>CCC` variants are
another allele. Operator decisions (2026-10-01): no margin base on the ALT side;
accept the exact-carrier rule's semantics for anchor-changing variants. Branch
`feature/cluster1-pure-indel-read-judgment`; harness `~/test/gbcms/harness/c20/`
(local). Follow-ups: R5 #198 and C27 #201 (6.7.0); C25 #199 (moved into 6.6.0); C26 #200 (6.6.0,
decide: which of a read's other indels decide its REF call, inside vs outside
the window); C28 #202 (6.6.0,
measure: anchor-deleting reads fall back to Phase 3's closer haplotype).

**In review — T1 #170 + T2 #171 (test architecture):** the legacy per-variant
`count_bam` path is retired. Binning invariance replaces parity: a Rust bin property
test, `tests/test_binning_invariance.py`, and `count_checked` (production vs one
variant per bin) in every counting test. A read census (`tests/census.py`) checks
pure-indel classification; the open decisions C26/C27/C28 are strict xfails there.
Branch `feature/t1-t2-test-architecture`; byte-identity harness
`~/test/gbcms/harness/t1/` (local). Next: H3 #204, the code-quality sweep (audited;
work list on the issue).

**Triage (2026-09-30):** every open issue was gone through by cluster with the
operator. 6.6.0 keeps 32 work items (plus the tracker #140 and umbrella #92); 16
moved to a new 6.7.0 milestone (each with its reason on the issue); C5 and P2
closed; umbrellas #133–#135 closed. See the plan's "Triage" and "Suggested
order". #92's checklist was verified: 10 items fixed by later PRs, #92 closed;
M4 #194 (merge NA cells, 6.6.0) and C24 #195 (6.7.0) filed. Deferred issues sit
under the 6.7.0 tracker #196.

**Validation standard** (adopted 2026-09-28): gbcms is a genotyper, not a caller;
every count-affecting change is accepted per read, against the reads' own bases, on
a matrix of data types (synthetic, realigned panel DNA, DNA without realignment,
RNA, public reference data), not one assay. See the plan's "Validation standard".
Where the right behaviour is unknown, measure it on that matrix and survey
community practice before deciding; record the decision in the issue.

## Next (in order; the plan's "Suggested order" is canonical)
1. T1 #170 (retire the legacy path; binning invariance) with T2 #171 (census
   oracle): design note first. Then a code-quality sweep of the cycle's code
   (duplication, unused code, silent failures, comments, logging, monitoring) as
   its own PR, H3 #204 (operator request, 2026-10-01).
2. C25 #199 (moved into 6.6.0), C26 #200 (measure first, decide) and C28 #202 (measure first); C17 #176
   (measure first); R4 #185 with O8 #186; C16 #174 (trace first).
3. Small batches (any order): C19+O7; I1+I2+C9; I3+I4; hygiene (H1, H2, #147,
   mkdocs pin, P3); O1+O2; S1+S2; M2+M4.
4. Release: D1, D2; D4 before the cut; D6 late; D5 last.


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

**Next → 6.6.0.** The plan is `CYCLE_6.6.0_PLAN.md` on develop (#117). It has
35 tickets, each with a GitHub issue in the **6.6.0 milestone**, all under the
tracking issue **#140**. Grouped work is filed as sub-issues: #92 → #141–#144,
#112 → #145–#147, #133 → #148–#149, #134 → #150–#152, #135 → #153–#154. The
other top-level issues are #106, #111, #114, #119–#139, #155 and #156. Start with the operator
decisions (C2, R1, R2, I1, I3, I4, S1/S2), then the count-affecting tickets,
each measured first.

Release comparison: the HPC run was not done before the 6.5.0 tag. The operator
runs it after 6.6.0, on the coverage-driven regression panel (D5, #155) instead
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
