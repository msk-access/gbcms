# Proposal: per-molecule observation export for gbcms

**Status:** v4 — **design complete, all decisions ruled.** Code-grounded analysis of `develop`,
plus maintainer rulings (one-invocation/one-BAM, Parquet path co-designed, emit N/Other, join key,
`_rs` internal + narrow public wrapper) and measurements against real ACCESS v1/v2 BAMs.
Ready to become a gbcms issue or a PR-1 branch.
**Type:** additive engine capability · **Compat:** existing API byte-identical · **Version:** MINOR (6.1.0)

> v1 of this draft proposed a flag on `count_bam_binned` that widened its return to a tuple,
> and emitted inside `count_variant_from_cache`. **Both were wrong** and are corrected below:
> a tuple return breaks every caller unconditionally (a pyo3 fn has one static return type),
> and `count_variant_from_cache` runs **twice** for decomposed variants, so emitting there
> double-emits contradictory calls with zero parity signal. All claims below are verified
> against `py-gbcms@develop` with file:line.

## 1. Do we need this at all?

**Yes — but the reason is single-source-of-truth, not performance.** For a targeted consumer
(a handful of loci) a consumer-side htslib pass is tractable. The real argument: the per-molecule
allele call is the output of a classifier stack gbcms **already runs and discards** at the
finalization loop (`engine.rs:1706`) — `hash_molecule` identity (`fragment.rs:277`),
`FragmentEvidence::observe` consensus + BAQ/N/structural handling, `resolve()` with structural
priority + qual-diff discard band (`fragment.rs:204`), fed by Smith-Waterman realignment. A
consumer that **re-derives** allele calls produces phasing calls that need not agree with the
`counts.ad/adf` it is phasing on — reporting cis/trans from numbers that silently disagree with
the counts in the same table.

**The guarantee is only real if the export is bound to counts by test (§5).** Ship it without
the fragment-level `adf/rdf/dpf` invariant test and you get the maintenance cost and none of the
consistency — at which point the consumer's own read pass was the honest choice.

**Seam (hard):** this is a *generic* capability — emit the resolved per-molecule call. **Zero
phasing/join/reversion logic in gbcms** (no "molecules covering ≥2 loci" pre-filter — that is
phasing, and belongs to the consumer). A generic variant-subset/min-coverage param is fine.

## 2. The record

```rust
#[pyclass]
pub struct Observation {
    #[pyo3(get)] pub variant_index: u32,  // join key; semantics-free; stable WITHIN one call
    #[pyo3(get)] pub molecule_hash: u64,  // hash_molecule(qname, umi); see caveats
    #[pyo3(get)] pub allele: u8,          // 0=Ref 1=Alt 2=N 3=Other (+ Python IntEnum)
    #[pyo3(get)] pub best_qual: u8,       // winning allele's best base qual; 0 for N/Other
}
```

Corrections vs v1:
- **Drop `is_duplex`.** `duplex_alt_count`/`singleton_alt_count` (`types.rs:284/287`) are
  **declared but never assigned** anywhere in `rust/src`; `FragmentEvidence` has no strand-family
  field. Populating it means new fgbio-consensus-tag logic in the engine → classification-adjacent,
  seam-breaking, parity-pressuring. **Epsilon = `10^(-best_qual/10)`** instead. (Honest caveat:
  `best_qual` is max-across-mates, an optimistic per-fragment proxy — document it.)
- **Tier is by-BAM, not a field.** The cfDNA workflow emits **four** BAMs — `duplex`, `simplex`,
  `all-unique` (= simplex ∪ duplex ∪ unique molecules), and `standard` (uncollapsed) — and gbcms is
  already run **per-BAM with counts merged** (`merge.py`; the `D_*`/`SD_*` fragment columns). Calling
  is done on `duplex`; evidence uses `duplex + simplex`. So the consumer invokes the export **once per
  BAM** and stamps the tier / ε floor (duplex `~1e-6`, simplex `~3.4e-5`, raw/standard `~1e-3`) from
  the BAM it chose — the engine needs **no duplex-family tracking** and the record needs **no tier
  field**. `best_qual` refines ε *within* a tier.
  **Set relationships matter (consumer-side, out of scope for gbcms but worth stating):** `duplex` and
  `simplex` are **disjoint** molecule sets, so their union is safe — that is why the workflow combines
  exactly those two. `all-unique` and `standard` are **supersets**, so mixing a BAM with its superset
  re-observes the same molecule and inflates linkage evidence. The consumer must never combine
  overlapping BAMs.
  **The standard (non-cfDNA) BAM is a first-class case, not a fallback:** tumor/WES/WGS is one
  uncollapsed BAM, `raw` ε — the same export, one invocation, no tier bookkeeping at all.
  (Verified: the read-end XOR at `engine.rs:1527` is guarded by `if amplicon_mode`, so a
  non-amplicon BAM correctly merges R1/R2 into one molecule by qname; `hash_molecule(qname,
  None)` works with no UMI.)
- **✅ RESOLVED — ONE invocation, one BAM, for every context.** cfDNA → `all-unique`;
  standard tumor/WES/WGS → the standard BAM. No union, no per-BAM bookkeeping, no special-casing:
  the cfDNA path and the standard path are the *same* call.
  **Why duplicates are not a concern:** the cfDNA BAMs are produced **after UMI collapsing**, so
  one read-pair *is* one molecule. Grouping by qname (gbcms's default when no UMI tag is given) is
  therefore exactly right, and PCR duplicates were already removed upstream by the collapse.
  **ε policy — and this is version-dependent (verified against real ACCESS v1 and v2 BAMs,
  tag keys + tier composition only, 3000 reads/BAM):**

  | | collapsing tool | tags on `all-unique` | tier recoverable from one BAM? |
  |---|---|---|---|
  | **ACCESS v2** | fgbio (`_collapsed_grouped`, bwa) | full set: `MI, RX, aD, bD, cD, cE, aE/bE, …` | **yes** |
  | **ACCESS v1** | Marianas | `AS, MQ, NM, RG, XS` only — **no fgbio tags** | no |

  In v2, `aD>0 && bD>0` **exactly reproduces the file partition**: the `duplex` BAM is 100%
  duplex by that rule and the `simplex` BAM is 100% single-strand, so the tags discriminate
  perfectly *within* `all-unique`. So for v2, **one invocation on `all-unique` yields both the
  complete molecule set AND full tier resolution** — the operational simplicity and the ε
  hierarchy, not a trade between them.
  **ACCESS v1 (Marianas) — verified in detail.** No per-read tier marker exists: a single `@RG`,
  aligner tags only (`AS/XS/NM/MQ`), and `YA/YM/YO/YX` turn out to be **ABRA2 realignment
  annotations** (values like `1:2361920:850M2I473M`), present in *both* duplex and simplex at
  0.2–0.4% — not tier. Tier is therefore recoverable **only from which file** a read is in. Also
  verified in a 20 kb window: `duplex ∩ simplex = 0` (**disjoint**, so the union is safe),
  `duplex ⊆ all-unique` and `simplex ⊆ all-unique` both **100%** (qnames are preserved identically
  across the collapsed BAMs), and composition is ~**2.6% duplex / 13.5% simplex / 84% singletons**.
  Three viable v1 policies, **none of which need a gbcms change**:
  **(a)** observe `duplex` + `simplex`, tier by file — matches ACCESS's own evidence convention
  ("call on duplex, evidence from duplex+simplex"); excludes singletons. *Recommended.*
  **(b)** one invocation on `all-unique` with a **uniform conservative ε** — simplest and uniform
  with v2/standard; underweights the 2.6% duplex molecules (the highest-value ones). Errs **toward
  indeterminate, never toward false-cis** (ADR-0016 §5), and costs little because "confident CIS
  needs duplex OR ≥2 concordant molecules" already forces ≥2 for non-duplex evidence.
  **(c)** one pass on `all-unique` + a cheap **qname-set lookup** at the two loci against the
  duplex/simplex BAMs to stamp tier — viable precisely because qnames are 100% preserved; full tier
  *and* all molecules, at the cost of consumer-side complexity.
  v0.3 calibration then **measures** what each policy's evidence is actually worth.

  **⚠ Corrects an earlier claim in this document.** "gbcms never needs a tier feature" holds only
  if we accept the v1-style uniform ε. To get tier from a *single* v2 invocation, **gbcms must read
  `aD`/`bD`** — which it currently does not (verified: no fgbio tag reads in `rust/src`). Crucially,
  that is a far smaller and safer change than feared: stamping a **tier on the exported
  `Observation`** is export-only and **cannot affect counts → parity-safe by construction**, unlike
  populating `duplex_alt_count` in `BaseCounts` (`types.rs:287`), which would change counts and
  pressure the parity suite. Scope: read two integer tags where `FragmentEvidence` is already
  built, and add one field to the record.
  **There are exactly four BAMs — `simplex`, `duplex`, `all-unique`, `standard` — and there is
  NO simplex+duplex BAM.** "SimplexDuplex" (ACCESS-Pipeline `python_tools/constants.py:331`,
  and the `SD_*` MAF columns) is a **counts-level reporting merge** of the two disjoint BAMs, not
  a file: gbcms is run on each BAM and the counts merged (`merge.py`). That merge existing at all
  is itself the evidence that `Simplex` and `Duplex` are disjoint sets — you would not merge
  counts from a BAM with its own superset.
- **⚠ Superseded background — why `all-unique` was rejected for phasing.** `all-unique` mixes
  duplex + simplex + singleton molecules in **one** file, so pointing at it yields observations
  with **no tier information** — and **gbcms does not read fgbio duplex tags anywhere**
  (verified: no `aD`/`bD`/`cD`/`cE` reads in `rust/src`; `duplex_alt_count`/`singleton_alt_count`
  at `types.rs:284/287` are declarations only). Every molecule would have to be weighted at the
  **worst** tier (~1e-3), discarding the duplex sensitivity that makes low-AF reversion calling
  work. Three ways out:
  **(a) observe `duplex` + `simplex` separately** — two invocations, disjoint sets, tier from BAM
  identity, ε hierarchy preserved; loses only *singletons*, which are the weakest evidence anyway.
  **No gbcms change.** ← recommended for PR-1.
  **(b) observe `all-unique` + add per-molecule tier to gbcms** — one file, singletons included at
  an honest high ε; requires a **new engine feature** (read fgbio tags, populate a tier) that is
  classification-adjacent and needs its own parity story. Possible follow-on.
  **(c) `all-unique` with no tier** — rejected: silently throws away the duplex advantage.
- **`allele` is 4-state, REF first-class.** `resolve()` returns `(is_ref, is_alt)` with REF a real
  state — **back-mutation/reversion IS the REF signal**, so REF must be explicit, not "absence of
  ALT." Split `NEITHER` into `N` vs `Other` via the already-tracked `has_n_base` (`fragment.rs:53`,
  routed at `engine.rs:1738`) — free, and QC-valuable (duplex-masking vs third-allele/contamination).
- **`best_qual = 0` for N/Other** (no called allele) — specified so consumers are deterministic.
- **`molecule_hash` semantics — documented on the API:** per-read-end (not per-molecule) in
  **amplicon** mode (`mol_hash ^= is_read1?0x1:0x2`, `engine.rs:1527`); comparable **only** within a
  single call with fixed `umi_tag`/`library_type`. Surface `library_type`/umi context with the result.
- **Optional QC-only fields** (`is_forward`, `insert_size: Option<i32>` — keep the `None`-when-TLEN=0
  distinction, `fragment.rs:49`) only if a QC consumer needs them; not in the hot phasing core.

Serves phasing (join + ε), LOH/allelic-imbalance, and duplex-masking/contamination QC — no phasing
semantics in gbcms.

## 3. API — new pyfunction + shared core (reject the flag+tuple)

```rust
// EXISTING symbol — signature + return BYTE-IDENTICAL; calls core, discards obs
#[pyfunction] pub fn count_bam_binned(...) -> PyResult<Vec<BaseCounts>>

// NEW sibling — the only entry point that emits. `observations_path`:
//   None       → return observations in-memory (targeted phasing; ~128 KB)
//   Some(path) → stream them to Parquet FROM RUST (genome-wide/ACCESS; avoids the PyO3 cliff),
//                returning an empty obs Vec + a row count
#[pyfunction] pub fn count_bam_binned_observations(...same params..., observations_path: Option<&str>)
    -> PyResult<(Vec<BaseCounts>, Vec<Observation>)>

// PRIVATE core refactored out of the current body
fn count_bam_binned_core(..., emit_obs: bool, obs_sink: ObsSink) -> Result<(Vec<BaseCounts>, Vec<Observation>)>
```

**Co-designed for both scales (per the maintainer's "also genome-wide/ACCESS-MAF").** The
`observations_path` param unifies the in-memory and at-scale paths in **one** function/core: with no
path, observations return over FFI (fine for a few loci); with a path, the core writes Parquet
directly (mirroring `parquet_writer.rs`/`write_fsd_parquet`) so 10⁷ rows **never cross the FFI
boundary** — the exact cliff `write_fsd_parquet` avoids. This is not a separate `write_obs_parquet`
(that would take the obs *back* through FFI, defeating the purpose); the write happens inside the
Rust core. Per-bin streaming into the Parquet writer (vs. collect-then-write) is a later memory
optimization, not required for v1.

Why: a pyo3 fn has **one** static return type, so widening `count_bam_binned` to a tuple breaks
**all** callers *unconditionally* (the `emit=false` default can't protect the return) — verified
sites: `pipeline.py:465` (`zip(valid_indices, counts_list)`), `pipeline.py:757`
(`snp_counts[offset]`), parity harness `helpers.py:165`, `test_pipeline_v2.py:123`, ~10 parity
fixtures, `_rs.pyi:193/225`. A **new pyfunction** keeps the existing return byte-identical and
matches house style (`count_bam`, `count_bam_binned`, `write_fsd_parquet` are already separate
pyfunctions, `lib.rs:16-20`). `Observation` is a **pyclass** (not a bare tuple/`u8`) so fields can
be added without an arity break. **`_rs`-layer only** — do **not** thread into `Pipeline`/CLI
(`write_fsd_parquet` precedent: an `_rs` export called post-count, not a CLI concept).
Add a **signature-parity test** so the two ~32-param wrappers can't drift.

### 3a. Public wrapper — RULED: keep `_rs` internal, expose one narrow supported API

`_rs` stays **internal** (88-field `BaseCounts`, 35-param `count_bam_binned` — a surface nobody
should promise to freeze). Consumers get one small, supported entry point instead:

```python
# src/gbcms/observations.py — re-exported in gbcms.__all__
def observe_molecules(
    bam: str | Path,
    variants: Sequence[Variant],           # already-public gbcms.Variant
    *,
    config: GbcmsDnaConfig | None = None,  # already-public config; defaults when None
    observations_path: str | Path | None = None,   # None → in-memory; path → Parquet
) -> ObservationResult: ...

@dataclass(frozen=True)
class ObservationResult:
    observations: list[Observation]   # empty when streamed to Parquet
    path: Path | None                 # where written, if any
    n_rows: int                       # always the true row count
```

**Observations only — deliberately no counts.** Counts already have a public path (`Pipeline`), and
`BaseCounts` is precisely the 88-field surface we want to keep internal. Consistency does not require
one *invocation*, only one *engine*: gbcms is deterministic, so `Pipeline` counts and these
observations agree by construction. It also matches the consumer's model — mulligan consumes VAF/AD
from its input MAF (ADR-0005) and needs only the per-molecule calls.

**Total public commitment:** one function (4 params) + `Observation` (4 fields) + `ObservationResult`
(3 fields), reusing the already-public `Variant` / `GbcmsDnaConfig`. Small enough to actually keep
stable — which is the point of the ruling. Adding a field to `Observation` later (e.g. a `tier` from
fgbio `aD`/`bD`) stays non-breaking under attribute access.
`_rs.count_bam_binned_observations` remains the internal implementation the wrapper delegates to;
state in CONTRIBUTING that **`gbcms._rs` is internal and not covered by SemVer**.

## 4. Consequences to existing workflows

| Surface | Option C (new fn) | v1 draft (tuple) |
|---|---|---|
| `pipeline.py:465/:757` | none | **breaks** |
| Parity harness + ~10 fixtures | none | **breaks all** |
| `_rs.pyi` stub | +1 additive | forced edit + mypy break |
| External MSK (ACCESS-Pipeline, cwl, nextflow) | **none** — all shell to the `gbcms` CLI in container `ghcr.io/msk-access/gbcms:6.0.0`; none touch `_rs` | **breaks** direct-`_rs` callers |
| Legacy `count_single_variant` oracle | untouched, non-emitting | untouched |

**Blast radius under Option C ≈ zero.** MINOR bump (6.1.0), `[Unreleased] → Added`.

## 5. Parity & determinism — the load-bearing part

**"Parity green" is necessary-not-sufficient:** the suite compares only `BaseCounts` fields
(`helpers.py`), so observations are invisible to it. Three holes it cannot see:

1. **Decomposed double-emit (most dangerous).** `count_variant_from_cache` runs **twice** per
   decomposed variant (`engine.rs:1084` original, `~1094` twin), arbitrated at `~1105`. Both run
   the `1706` loop → emitting there pushes two contradictory allele sets per molecule (incl. the
   *losing* form); counts keep only the winner, so **parity stays green while the export is
   corrupt.** → Emit into a returned per-call `Vec`; `count_bin_shared` keeps only the **winning**
   invocation's observations at arbitration.
2. **Threading nondeterminism.** Order varies from the cross-bin `try_reduce` merge (`engine.rs:816`,
   work-stealing) and `HashMap::iter` (`RandomState`) at `1706`. The hash *value* is deterministic
   (`DefaultHasher`, `fragment.rs:277`); *order* is not; counts are order-invariant so parity never
   flags it. → **Stable `sort_by_key(|o| (variant_index, molecule_hash))`** after merge.
3. **RNA transcript loop is the wrong site — but RNA-mode runs are still covered.**
   `count_per_transcript` (`engine.rs:2411`) resolves per-transcript `tx_fragments` at `~2567` → K
   duplicate-key rows per molecule; it is also a *reporting* path (returns formatted
   `"ENST…:AD,RD,DP"` strings, not data). The `1706` path **already runs in RNA mode**.
   → **Emit only from `1706`.** Consequence: **RNA-mode observations work in PR-1** — same record,
   keyed by molecule, which is what RNA-read phasing needs (spliced reads co-locate loci that are far
   apart genomically). Only **isoform-resolved** export is deferred (key `(variant, transcript,
   molecule)`, which the flat record cannot represent) — that answers *which isoform* a molecule
   supports, a different question from phasing.

**Concurrency shape:** the v1 `Option<&mut Vec>` sink through the rayon `par_iter` (`engine.rs:748`)
is **not implementable** (shared `&mut` can't cross the closure; a `Mutex` serializes hot bins). Use
a **per-bin owned `Vec<Observation>`** carried alongside the existing `Vec<(usize, BaseCounts)>`
(`count_bin_shared`, `engine.rs:970`) and merged in the same `try_reduce`.

**No independent oracle → the invariant test is the only guard, and must be fragment-level:**
> per variant: `#obs(Alt)==counts.adf` ∧ `#obs(Ref)==counts.rdf` ∧ `len(obs)==counts.dpf`.

Read-level (`ad/rd/dp`) would spuriously fail — a read dropped by the multi-allelic sibling guard
still carries REF into the fragment (`engine.rs:1542`, before the `continue`), so it counts in `rdf`
and emits `Ref`. **Tests:** (1) fragment-level `adf/rdf/dpf`↔obs invariant; (2) decomposed case;
(3) multi-allelic sibling case; (4) determinism (two runs → identical sorted Vec); (5) flag-off
regression (`BaseCounts` frozen); (6) all obs comparisons order-insensitive; (7) signature-parity.

## 5a. Validation on real data (post-implementation)

Run against real MSK-ACCESS v2 cfDNA (GRCh37) and a STAR RNA BAM (GRCh38), on sites
discovered from the BAMs themselves. Statistics only — no identifiers or sequence.

| check | DNA (ACCESS v2) | RNA (STAR) |
|---|---|---|
| fragment-level invariant (`#ALT==adf`, `#REF==rdf`, `len==dpf`) | **0 failures / 39,561 rows** over 40 variants | **0 failures / 3,274 rows** |
| determinism (6 runs × threads 1/4/8) | identical, sorted, unique keys | identical |
| counts vs `count_bam_binned` | identical | identical |
| cross-locus molecule linking | 7,607 / 18,050 molecules at >1 variant | 463 / 487 |

**Answers to the two narrowing questions raised in review** — measured, not argued:

- **`mode` (and RNA annotation) are irrelevant to the export.** On the STAR BAM,
  `mode="dna"`, `mode="rna"`, and `mode="rna"` **with the real 1.5 GB GRCh38.111 GTF loaded**
  (1.65 M exons, annotation active) all produced **byte-identical** observations — same rows,
  same molecule hashes, same alleles, same qualities. RNA-specific work (ASJD, per-transcript
  counts) lands in RNA *counts*, not in the per-fragment resolution the export reads. So the
  wrapper not forwarding `mode`/`gtf_path` costs the export nothing.
- **`sibling_variants` and `alignment_backend` did not change the export** on real SNVs
  (15,686 rows, identical allele mix for sw vs pairhmm, and with vs without siblings) — nor on
  real deletions. The backend is forwarded anyway (the config default is `pairhmm` while the
  Rust default is `sw`, so not forwarding it would have made the *same config object* count
  differently here than through `Pipeline`).
- **Normalization measurably matters for indels.** On real deletion sites the public wrapper
  (normalize + decompose) recovered **13 ALT molecules vs 11** without it, and the
  un-normalized path logged `check_deletion: ref_context is None … cannot validate deleted
  bases`. This is the defect the first implementation shipped with.
- **`N` rows appear only below the base-quality gate.** At the default `min_baseq=20`: 0 N /
  300 OTHER. At `min_baseq=0`: **202 N / 97 OTHER** — the same molecules, reclassified. N bases
  carry low quality, so they are filtered before reaching fragment evidence unless the gate is
  relaxed. Worth documenting: the N/OTHER split is real but conditional on `min_baseq`.

**End-to-end cross-locus use.** Linking molecules by `molecule_hash` across two nearby real
variants and weighting each by its base quality yields a usable phase statistic — e.g. a pair
38 bp apart with 256/402 molecules spanning both loci. Note that interpreting *ALT+REF* is the
caller's job, not gbcms's: against a homozygous background every molecule reads ALT+REF, which
is not "trans" in any meaningful sense. The engine reports what each molecule showed; what it
means is the consumer's problem — exactly the seam this proposal preserves.

## 6. Performance & memory

- **Flag-off: zero** (mirror the mFSD gate `Vec::with_capacity(0)` at `engine.rs:1700`).
- Memory ≈ Σ per-variant deduped fragment count (the `1706` map), **not** `#frags × #variants`.
  - Targeted phasing (~2 loci): ~**128 KB** → in-memory return is ideal.
  - ACCESS panel MAF: ~10⁶–10⁷ obs ≈ **30–320 MB**; genome-wide ≈ **~160 MB**, ×Nextflow fan-out.
- **The real cliff is PyO3 materialization** (10⁷ PyObjects, doubled peak RSS) — exactly what
  `write_fsd_parquet` was built to avoid. Because the maintainer wants the **genome-wide/ACCESS-MAF**
  path too, the `observations_path` param (§3) is **co-designed into v1**: no path → in-memory
  (targeted phasing); path → Parquet written from Rust (no FFI crossing). One function serves both.

## 7. Alternatives

| | Option | Compat | Parity | Rayon-safe | Seam | Verdict |
|---|---|---|---|---|---|---|
| A | flag + tuple return (v1) | ❌ breaks all | fixtures break | sink not impl. | ok | **reject** |
| B | `CountResult` pyclass, one entry | ❌ same return break | fixtures break | ok | ok | reject now |
| **C** | **new fn + shared core** | ✅ byte-identical | green by construction | ✅ per-bin Vec | ✅ `_rs`-only | **PICK** |
| D | Parquet sink (`observations_path`) | ✅ additive | green | ✅ serial | ok | **folded into C** (at-scale mode) |
| E | consumer keeps own htslib pass | ✅ no change | n/a | n/a | ✅ | reject — calls **drift** from counts (§1) |

## 8. Plan & open decisions

**PR-1 (in-memory return; covers DNA *and* RNA mode):** refactor body into
`count_bam_binned_core(emit_obs, obs_sink)`; add `Observation` pyclass + `.pyi` + Python `IntEnum`;
emit at `1706`, keep winning decomposed set at `~1105`, per-bin Vec in `try_reduce`, stable sort; add
`count_bam_binned_observations` (with an inert `observations_path=None`) in `lib.rs`; **add the
public wrapper `gbcms.observe_molecules` + `ObservationResult` (§3a), re-export in `__all__`, and
document in CONTRIBUTING that `gbcms._rs` is internal / not SemVer-covered**; document
`molecule_hash` + by-BAM tier + the ACCESS `filter_improper_pair` note; tests 1–7 **plus a wrapper
smoke test**; CHANGELOG; MINOR bump.
**PR-2 (DONE — Parquet sink active):** wire `observations_path=Some(...)` to the Rust-side Parquet
writer (mirroring `parquet_writer.rs`) for the genome-wide/ACCESS cliff — additive, the signature is
already in place from PR-1.
**PR-3 (separate):** **isoform-resolved** RNA export (`(variant, transcript, molecule)` key) — for
"which isoform does this molecule support", *not* needed for RNA phasing (PR-1 covers that).

**Resolved by the maintainer:** *tier* = by-BAM, **conditional on a tier-homogeneous BAM** (§2 ⚠);
*workload* = both scales → Parquet path co-designed into the signature now (§3); RNA-mode covered by
PR-1, isoform-resolved deferred (§5.3).

**Also resolved:** *which BAM(s)* — **one invocation per tier-homogeneous BAM** (standard = 1;
cfDNA = duplex + simplex, disjoint, unioned). **No gbcms tier feature is needed, now or later** (§2).
With that settled, **nothing further changes gbcms scope** — the remaining items are shape choices.

**All open decisions are now RULED — the design is complete.**

1. **Is `gbcms._rs` part of the SemVer public contract?** → **RULED: no — keep `_rs` internal AND add
   a narrow public wrapper** (`gbcms.observe_molecules`, §3a). `_rs` is documented as internal in
   CONTRIBUTING and stays free to evolve; consumers depend on a 4-param function + two small records.
   Best of both: a supported surface for the consumer, no frozen 88-field/35-param commitment.
2. **Emit N/Other rows?** → **Emit.** **Measured on real ACCESS v1 *and* v2 BAMs** (exact
   classification vs the b37 reference, `ignore_orphans=False` — see the ⚠ below):

   | BAM | REF | N | non-REF |
   |---|---|---|---|
   | v2 all-unique (phasing input) | 98.72% | **0.33%** | 0.96% |
   | v2 duplex | 99.60% | 0.40% | 0.001% |
   | v1 all-unique | 98.71% | **1.13%** | 0.16% |
   | v1 duplex | 99.97% | 0.03% | 0.007% |
   | v2 / v1 standard (uncollapsed) | 99.67 / 99.81% | **0.003 / 0.002%** | 0.33 / 0.19% |

   Most "non-REF" is *true ALT* (germline/somatic), which the export labels **ALT**, not Other. So
   **N/Other rows are ≈0.03–1.1% of the total** across every platform — dropping them saves ~1% at
   most, nowhere near enough to trade away an invariant. Two corroborating findings: N is ~**0.002%**
   in uncollapsed BAMs vs **0.03–1.1%** in consensus BAMs, i.e. N marks **strand-discordant**
   molecules and those rows carry real QC information; and non-REF runs **0.001% (v2 duplex) vs
   0.33% (standard)**, a ~300× gap that independently confirms the ε hierarchy.

   > **Note (measured): ACCESS v1 BAMs are 100% paired but 0.0% `PROPER_PAIR`** (all-unique, duplex
   > *and* standard; v2 is 93.7–100%). Any proper-pair filter therefore discards **100% of v1 reads**
   > silently. **gbcms already handles this correctly** — `filter_improper_pair` defaults to
   > **`False`** in both the CLI (`cli.py:280`, `:673`) and the config model (`models/core.py:94`),
   > so out-of-the-box runs are safe; it is caller-controlled and simply should not be enabled for
   > ACCESS. No gbcms change needed. (Worth stating in the export's docs anyway: the same trap in
   > pysam — whose `pileup(ignore_orphans=True)` default is the *opposite* choice — produced an empty
   > first measurement here, and the consumer's own read pass must likewise never filter on it.) (ii) The `len==dpf` invariant is **load-bearing** (no second implementation cross-checks
   the export). (iii) **Span breaks without them**: a molecule that is Ref at locus 1 and N at locus 2
   would emit one row and be indistinguishable from one that never spanned, destroying the
   "looked-and-saw-nothing" vs "couldn't-look" distinction (ADR-0006). If memory bites on the Parquet
   path, add an opt-in `min_informative` filter *then*.
3. **Join key?** → **Both, split by path.** In-memory: **`variant_index`** — compact, unambiguous
   within the call, and already a first-class internal concept (the `all_counts[vi]` scatter tracks
   original index through binning). Parquet: **additionally echo `(chrom,pos,ref,alt)`** so a
   persisted file is self-describing once it outlives the call; columnar dictionary-encoding makes
   those columns nearly free on disk.

**Criticality note:** allele-specific CN is unavailable for cfDNA (FACETS is IMPACT-only; the cfDNA
CNV caller is coverage/loess with no `lcn` — verified), so **LOH-based phasing cannot run there**.
For the primary cfDNA workload this export is therefore the *only* phasing evidence, not a
consistency bonus — which is the argument for doing it properly (the fragment-level invariant test)
rather than cheaply.

**Verify before coding:** (a) `count_variant_from_cache` returning `(BaseCounts, Vec<Observation>)`
threads cleanly through both decomposed call sites + arbitration; (b) `best_*_qual`/orientation in
scope at `1706`; (c) no existing determinism test asserts on returned-data order beyond `BaseCounts`.
