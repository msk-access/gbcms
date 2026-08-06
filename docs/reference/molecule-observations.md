# Per-Molecule Observations

The allele each individual molecule carried at each variant — the layer beneath `BaseCounts`.

## Overview

Standard output tells you **how many** molecules supported REF and ALT at a position. Some
analyses need to know **which molecule** carried **which allele** — because the same molecule
observed at two positions links them.

gbcms already resolves that per-fragment call while counting, then aggregates it away. The
observation export surfaces it instead of discarding it. Counting is unaffected: the same
classifier produces both, so the rows can never disagree with the counts beside them.

!!! tip "Quick Start"
    ```bash
    gbcms dna --observations-parquet \
      -v variants.maf -b sample.bam -f ref.fa -o output/
    ```
    Produces the usual VCF/MAF **plus** `<sample>.observations.parquet` with one row per
    molecule per variant.

**What it enables.** A molecule seen at two variants in the same run carries the same
`molecule_hash`, so joining on it answers questions counts cannot:

- **Read-backed phasing** — are two variants on the *same* molecule (cis) or different ones (trans)?
- **Allelic imbalance** — per-molecule allele ratios rather than pooled totals.
- **Duplex-masking QC** — how many molecules were strand-discordant (`N`) at a locus.

!!! warning "gbcms does not interpret"
    The export reports what each molecule showed. It performs **no linking, phasing, or
    calling** — deciding what the pattern *means* belongs to the consumer. Against a
    homozygous background, for instance, every molecule reads ALT+REF, which is not "trans"
    in any meaningful sense.

---

## CLI Options

| Flag | Default | Description |
|:-----|:--------|:------------|
| `--observations-parquet` | `false` | Write a companion `<sample>.observations.parquet` with the per-molecule allele call at each variant. Available on `dna` and `rna`. Counts output is unchanged. |

Nothing is written unless the flag is passed, and no per-molecule rows are built when it is
absent — the cost is zero when off.

---

## Output Schema

One row per **fragment** per variant. Every row echoes its locus, so the file stands alone
once written.

| Column | Type | Description |
|:-------|:-----|:------------|
| `variant_index` | `uint32` | Index into the input variant list |
| `chrom` | `string` | Chromosome |
| `pos` | `int64` | 0-based position |
| `ref` | `string` | Reference allele |
| `alt` | `string` | Alternate allele |
| `molecule_hash` | `uint64` | Molecule identity — **the join key across variants** |
| `allele` | `uint8` | `0`=REF, `1`=ALT, `2`=N, `3`=OTHER |
| `best_qual` | `uint8` | Best base quality supporting the called allele; `0` for N/OTHER |
| `min_mapq` | `uint8` | Worst MAPQ among the reads backing this molecule; `255` = unavailable |

Rows are sorted by `(variant_index, molecule_hash)`, so output is byte-reproducible run to run
despite parallel processing. Compression is Zstandard level 1.

### Allele values

| Value | Meaning |
|:------|:--------|
| `0` REF | Molecule carried the reference allele. **First-class** — reference observations are themselves signal for some analyses (e.g. back-mutation). |
| `1` ALT | Molecule carried the alternate allele. |
| `2` N | Ambiguous base. In consensus BAMs this marks a **strand-discordant** molecule — diagnostic, not noise. |
| `3` OTHER | A third allele, or no fragment consensus. |

!!! info "N rows depend on `--min-baseq`"
    N bases carry low quality, so they are usually filtered before reaching fragment
    evidence. Measured on real cfDNA: at the default `--min-baseq 20` a locus set yielded
    0 N / 300 OTHER rows; at `--min-baseq 0` the same molecules gave 202 N / 97 OTHER.

### `min_mapq` semantics

The **minimum** MAPQ across the reads that contributed evidence to the molecule — not the
best, and not an average. A molecule is only as trustworthy as its least confidently placed
read, so a consumer weighting evidence by error probability needs the pessimistic bound.
For a paired fragment this means the worse mate's value: one cleanly-placed mate cannot
launder an ambiguously-placed one.

Two values are easy to confuse and mean opposite things:

| Value | Meaning |
|:------|:--------|
| `0` | A **real** MAPQ — the read mapped ambiguously (multi-mapping). Trust it least. |
| `255` | **Unavailable** (SAM spec). No mapping quality was recorded; not a confidence claim. |

!!! note "Why this is exported rather than inferred from `--min-mapq`"
    The filter threshold *bounds* mapping confidence without measuring it. Assuming every
    surviving molecule sat at the threshold would charge them all the same error — swamping
    the base-quality term and flattening precisely the per-molecule differences that make a
    weighted statistic worth computing. Reads that clear a `--min-mapq 20` gate are mostly
    MAPQ 60; treating them as 20 discards that.

Reads below `--min-mapq` never reach fragment evidence, so this value is always at or above
the threshold you ran with.

#### Why `min_mapq` is a minimum but `best_qual` is a maximum

The two aggregate in opposite directions on purpose — they describe different failure modes.

**MAPQ asks whether the read is in the right place.** If *any* read backing the molecule is
badly placed, the molecule's locus assignment is suspect, so the weakest link governs.

**Base quality asks how strong the evidence for the called allele is.** Reads that agree
reinforce each other: if two reads both call ALT at Q40 and Q20, the chance both are wrong is
roughly the product of their error rates — *better* than either alone. Taking the maximum
therefore **understates** the combined evidence rather than flattering it.

Reads that *disagree* are not handled by either number: fragment resolution arbitrates them,
and a fragment with no consensus is reported as `OTHER` with `best_qual = 0`. So the case a
"worst base quality" column would be for is already visible in the `allele` field.

There is deliberately no `min_baseq` column for the same reason there is a `min_mapq` one:
bases below `--min-baseq` are filtered before reaching fragment evidence, so such a column
would mostly re-report the threshold you passed — whereas MAPQ genuinely varies above its
threshold and is worth measuring.

### `molecule_hash` semantics

`molecule_hash` is `hash_molecule(qname, umi)` and is comparable **only within a single run**
using one `--umi-tag` and library type. Two caveats:

- In **amplicon** mode the read number is folded into the hash, so the key is per-read-end
  rather than per-fragment.
- It is a 64-bit hash with no stored qname. Collisions are negligible per locus but aggregate
  over very large cross-locus joins.

---

## Python API

For targeted work — a handful of loci — the rows come back in memory:

```python
from gbcms import Variant, VariantType, observe_molecules

variants = [Variant(chrom="17", pos=7676153, ref="T", alt="A",
                    variant_type=VariantType.SNP)]

result = observe_molecules("sample.bam", variants, reference_fasta="ref.fa")

for obs in result.observations:
    print(obs.variant_index, obs.molecule_hash, obs.allele, obs.best_qual, obs.min_mapq)
```

| Parameter | Description |
|:----------|:------------|
| `bam` | Indexed BAM/CRAM. |
| `variants` | Variants to observe. `variant_index` indexes into this list **positionally** — it is never filtered or reordered, so the join key stays meaningful even when a variant fails validation. |
| `reference_fasta` | Reference for normalization. **Strongly recommended for indels** (see below). |
| `is_maf` | Set when variants came from a MAF, whose `-` alleles need anchor resolution. |
| `config` | `GbcmsDnaConfig` for filters, quality gates, alignment backend, threads. |
| `observations_path` | Write to Parquet instead of returning rows (see [Scale](#scale)). |

`ObservationResult` carries `observations`, `path`, `n_rows`, and `variant_status` — one status
per input variant. Anything other than `PASS` means normalization rejected that variant (e.g.
`REF_MISMATCH`) and its rows should be discarded.

!!! warning "Always pass a reference for indels"
    Normalization supplies left-alignment, `ref_context`, and the **decomposed** form of
    complex indels — often the form the reads actually carry. Without it those molecules are
    scored `OTHER` and the variant appears to have **zero ALT support**. Measured on real
    ACCESS deletions: 13 ALT molecules with a reference, 11 without.

!!! note "`gbcms._rs` is internal"
    `observe_molecules` is the supported surface. `gbcms._rs` is an implementation detail and
    is **not covered by SemVer** — see [Contributing](../development/contributing.md).

### Alignment backend

Observations are emitted **after** classification, so they inherit whatever the counting
backend decided — the export has no opinion of its own. Rows always reconcile with the counts
beside them regardless of backend, but the underlying calls can differ between `sw` and
`pairhmm` on ambiguous indels (measured: 5 of 104 molecule calls on real ACCESS deletions).
Compare observation sets only when they were produced under the same backend.

`pairhmm` is the default at **every** layer — the CLI, `Pipeline`, `observe_molecules`, and
the `gbcms._rs` entry points alike. An unrecognized backend token raises `ValueError` rather
than quietly falling back, so a typo cannot silently score your data with the other
classifier.

WFA is **not** a third choice: it is an edit-distance triage step in front of PairHMM that
resolves clear-cut reads directly and defers only ambiguous ones to the full model.

---

## Scale

Two paths, one capability:

| Workload | Approach |
|:---------|:---------|
| **Targeted** (a few loci) | Omit `observations_path` — rows return in memory. A two-locus query is ~128 KB. |
| **Panel / genome-wide** | Pass `observations_path=` (or use `--observations-parquet`). Rows are written **from Rust** and never become Python objects. |

At panel scale a run yields 10⁶–10⁷ observations. Materializing that many Python objects costs
GC pressure and roughly doubles peak RSS — under Nextflow per-sample fan-out that is where a
worker dies. Counting was never the bottleneck; the FFI boundary is. Writing from Rust avoids
it entirely, and the columnar format compresses well (39,561 real rows → 227 KB).

---

## Read Filtering Caveat

Molecules are grouped by `hash_molecule(qname, umi)`, so reads must not be pre-filtered in a
way that removes whole molecules.

!!! danger "Do not enable `--filter-improper-pair` on consensus BAMs"
    Some UMI-collapsed pipelines emit reads with **no `PROPER_PAIR` flag at all** — measured
    at 0.0% across every BAM type in one ACCESS version. Filtering on that flag discards
    **100% of reads**, silently. gbcms defaults this filter to `false`, so out-of-the-box runs
    are safe; simply do not turn it on for such data.

---

## Guarantees

- **Counts are unchanged.** Both entry points share one core, so enabling observations cannot
  alter a single count. The existing `count_bam_binned` signature and return are untouched.
- **Rows reconcile with fragment-level counts.** For every variant, `#ALT == adf`,
  `#REF == rdf`, and `total rows == dpf`. This is asserted in the test suite — note it is
  *fragment*-level, not read-level, because a read excluded by the multi-allelic sibling guard
  still contributes REF evidence to its fragment.
- **Deterministic output.** Stable sort by `(variant_index, molecule_hash)`; verified identical
  across repeated runs at 1, 4 and 8 threads on real BAMs.
- **Decomposed variants emit one set of rows.** When a decomposed form wins arbitration, the
  observations follow the counts — the losing allele form is never exported.

---

## Related

- [DNA Command](../cli/dna.md) — CLI flag reference
- [RNA Command](../cli/rna.md) — same flag in RNA mode
- [Counting & Metrics](counting-metrics.md) — the aggregate counts these rows underlie
- [Output Formats](output-formats.md) — VCF and MAF output structure
- [Variant Normalization](variant-normalization.md) — why a reference matters for indels
- [Read Filters](read-filters.md) — filter semantics, including `--filter-improper-pair`
