# Proposal: the sequence judge (package 2, #304)

**Status:** v1 — architecture decided (operator, 2026-10-23; record on #304); this document
fixes how it enters the engine. **Type:** read judgment, behind `--read-judge sequence`
(default `rules`; `haplotype`, the judge's name until 2026-10-26, is taken as `sequence`
with a warning). **Compat:** default output byte-identical; the flag changes counts.
**Version:** 6.7.0 (default flip only if every arm, the cost target and the HPC gate pass).
All file:line references are to develop 336c07cd. Internal names (`judge::haplotypes`,
`ReadJudge::Haplotype`) keep the old word until the module is restructured (R6, #314).

## 0. The principle (operator, 2026-10-25)
The judge resolves the read's sequence at the row: which local sequence its own bases spell
there, REF, the row's ALT or another allele at its site. Contexts (the other listed rows, the
guard's recurring background, the alleles the reads spell away from the site) only explain
the read's *other* differences, so they never pass for the row's allele (#246: an insertion
beside an SNV whose first inserted base reads as the SNV's ALT). The judge does no linkage
and no phasing: a variant linked to the row elsewhere on the molecule never credits a read.
The sequences a read is weighed against are the row's *candidate sequences*.

## 1. The rule
For a row X and each admitted read:
- **Candidate sequences:** every *context* of the other listed rows within read reach of X
  (and the guard's recurring background changes, which explain reads but are never
  credited), crossed with every *state* of X's site: REF, X, each other listed allele whose
  REF span overlaps X's (an insertion's span is its anchor), and each recurring discovered
  allele there (RJ-25; REJ-20261007-002).
- **Scorer:** the read, every readable base, against each candidate sequence with its ends
  free: affine edit cost, confident mismatch 4, gap 6 + k, a base below `--min-baseq` or N
  costing 0 (the one quality gate).
- **Decision:** ALT where the best candidate sequence carrying X beats the best without it
  by one confident base (4); REF where the best with REF at the site beats every other site state by
  as much; partial evidence where another site allele wins so; otherwise depth only. Ties are
  never credited. Because every context appears in every state, only the read's bases at X
  decide; a linked variant elsewhere never credits it.
- **Phase C refinements (2026-10-25):**
  - A state holding X's change at X's site (for a substitution X, an allele of its length
    over X's bases that reads X's bases there, as an MNV listed beside its component SNV)
    is X's ALT (RJ-28); `MNV_SHARED` names the reads the two rows share.
  - A same-length state's changes beside a substitution X are also contexts of every state,
    so a read carrying X with such a change is X's.
  - A substitution inside a long (50+ bp) deletion's bases, past a 4-base junction window at
    each end, is a context on REF's side, never another allele at X's site (RJ-24); inside
    a shorter deletion it stays another allele (RJ-5).
  - A tie of REF with another allele, X beaten, is partial evidence only where the read
    reaches past where the two candidate sequences differ; a read short of it shows REF.
  - An MNP's ALT call needs every changed base covered and none read as another base; one
    below `--min-baseq` weighs nothing (RJ-25), and `mnp_confirmed` needs all read.
  - A context row is one whose REF span meets X's span widened by the longest read, so a
    long deletion starting within reach (and running past it) is a context too.
- **The fit rule (phase B, decided 2026-10-25; RJ-29):** a read whose best score over every
  candidate sequence is below ten confident mismatches fits none of them; it is depth only
  where its aligner also doubts it (XS ≥ AS, or MAPQ 0 with no XS), otherwise its margin
  decides it. `UNEXPLAINED_READS(n/m)` names how many of a row's depth reads fit none. The
  other options measured stay behind `GBCMS_JUDGE_FITS` (measurement only).

Measured (HG002 phased truth, R2): see #304. The prototype is `plan670/r2/hapjudge` (local).

## 2. Where it enters
| Seam | Where | Change |
|---|---|---|
| Flag | `cli.py:387/755` → `models/core.py:159` (`AlignmentConfig`) → `pipeline.py:887 _engine_kwargs` → `#[pyo3(signature)]` `engine.rs:924/985` → `count_bam_binned_core` (385) → `count_bin_shared` (1076) → `count_variant_from_cache` (1477); stub `_rs.pyi:225–296` (parity test) | `read_judge: "rules" \| "haplotype"`, plumbed like `alignment_backend`; also the MNP-rescue recount (`pipeline.py:1340`) and Nextflow (`dna/main.nf:60`) |
| Rows in reach | prep, a whole-input post-pass after `assign_listed_neighbours` (`normalize/engine.rs:131`) | new `Variant::reach_rows`: every other PASS row whose REF span lies within X's ± `READ_REACH` (300, as `carrier::WHOLE_READ_REACH`); computed also for the decomposed twin (1066) and the MNP-rescue rows (prepared apart, `pipeline.py:1306`) |
| Reference | prep, beside `whole_ref` (1104–1105) | `Variant::judge_ref`: X ± (`READ_REACH` + the widest reach row), for every kind (SNVs and deletions ≥ 50 bp have no `whole_ref` today) |
| Classifier | `engine.rs:1662`, after `splice_skip_triage` | under the flag, `judge::classify` returns a `ClassifyResult`; the sibling guards (1691, 1699) are bypassed (the judge weighs the site's alleles itself) |
| Background | the guard, forced before the loop (`engine.rs:1540`; `carrier::guard_for` 2655) | `Guard.background` → context rows; `Guard.others` → site states; for SNVs (no guard) the site holds REF, X and the listed alleles |

Nothing else changes: the read cache, fragment evidence (`fragment.rs:147/258`), RD/AD/DP,
strand, mFSD and diagnostics consume the `ClassifyResult` as today.

## 3. Mapping a decision to `ClassifyResult` (`utils.rs:70`)
- **REF / ALT:** `is_ref` / `is_alt`; `qual` = the lowest quality among the read's bases where
  the winning and the best losing candidate sequence differ (the bases that decided it; at least
  `--min-baseq`), so mate resolution (`resolve`, `fragment.rs:258`) weighs a decision by the
  bases it rests on. `is_structural`: to be measured both ways (§10): `false` for every
  call, or `true` where the read's CIGAR writes X's exact op.
- **Partial:** neither; `partial_match_count = 1` (another site allele won).
- **Depth only:** neither, `uninformative = true`.
- **Admission:** `covers_locus` from splice triage as today; `clip_admissible` (DNA) where a
  soft clip carries bases over X's span (RJ-13, clipped bases are the read's);
  `read_from_right` where the read's first readable base lies past X's left flank. A read is
  admitted to DP when its aligned span, clips included in DNA, reaches X's REF span; for a
  deletion of 50 bp or more, its anchor (RJ-24): reads wholly inside the deleted bases are not
  admitted, as today.

## 4. The scorer and cost
- Semi-global affine DP (read global, candidate sequence local) with a wildcard for masked
  bases, banded around the read's aligned diagonal mapped into each candidate sequence, the
  band widened by its length change. Per read, candidate sequences are trimmed to its reach.
- **Triage first:** a 31-mer index of the variant's candidate sequences; a read whose fully
  confident k-mers fall in one site state's sequences only, with no unknown k-mer over X's span ± k,
  is decided without alignment (R2 round 1: ~40% of reads).
- **Target:** within 2× develop's runtime on the panel, ACCESS (deep) and WES arms. The
  prototype is ~300× develop; SIMD banding (e.g. the `block-aligner` crate, which takes
  position-specific scores for the mask) is the next lever if the plain banded DP is not
  enough. `wfa2lib-rs` (in-tree) aligns ends-free and affine but has no wildcard.

## 5. Contexts and states
- Rows in reach ≤ 4 (outside the site): every non-overlapping subset is a context.
- More: singletons, pairs a read shows together (per-row evidence on ≥ 2 reads), then
  read-supported sets, to 32 candidate sequences; the empty context always.
- States: REF, X, the site's listed alleles, and for an SNV every other base; for an MNP or
  an indel the guard's discovered alleles there. A state holding X's change is of X's class
  (RJ-28); a substitution in a long deletion's bases past its junction windows is a context
  on REF's side (§1).
- Deterministic order everywhere (rows by position; ties to fewer rows); per variant only.

## 6. RNA
Spliced reads are judged on candidate sequences built along their mRNA: the read's segments
(`carrier::splice_segments`, 1471) give the spliced reference, rows mapped through
`SpliceMap` (1364); a splice through X's site leaves the read depth only (as `carrier.rs:726`).
Strandedness and antisense handling are unchanged (engine.rs:1637/1721/1924). Measured with
package 3 on FORTE.

## 7. What it replaces, and what stays
Under the flag, `check_snp`, `check_mnp`, the exact-carrier rule, `junction::classify` and
`classify_complex` are not called for the row's own counts; the guard still runs (it supplies
background and discovered alleles). Per-transcript counts (3096) and ASJD (3704) call the same
`judge::classify`, so every output agrees. `Tolerance` and the error census are not used by
the judge (a recurring change becomes a context row instead).

## 8. Tests
- `tests/test_haplotype_judge_contract.py` (red-first): #246 in trans and in cis, a compound
  het, a repeat unit count, a low-quality telling base, no phase credit, a soft-clipped long
  insertion, admission at a long deletion's anchor.
- Every counting contract runs with the flag too where its shape is the judge's
  (parametrized); `tests/census.py` stays the default rules' check — the judge's independent
  checks are the R2 truth arms (HG002 phased truth) and the BWA read checks on local arms.
- Binning invariance (`tests/test_binning_invariance.py`) with the flag: every geometry, one
  call per row, decoys in the call.
- Rust unit tests: the scorer (wildcard, ends-free, affine), the construction, the decision.

## 9. Acceptance (before the default flips)
Every arm against develop and truth: the HG002 strata and R2 arms (#246 pairs, tandem
repeats, SVs), the compound hets; panel, ACCESS and normals, WES, the gate set against
sign-out, FORTE (with package 3); the IMPACT 50–500 bp events read by read (ABRA2 op, `YO`
placement, bases). Cost within 2× develop. The HPC gate.

## 10. Decided and open points (operator, 2026-10-23)
1. **Decided:** a judged call's `qual` is the lowest quality among its deciding bases.
2. **Decided:** a third base at an SNV is partial evidence: SNV sites hold all four bases as
   states; MNP and indel sites take the guard's discovered alleles.
3. **To measure, then decide:** `is_structural` under the flag — retired (every call rests on
   its bases; mates resolve by quality) or kept for reads whose CIGAR writes X's exact op.
   The build carries a measurement-only switch for both; the arms decide, and the switch is
   removed before the merge.
