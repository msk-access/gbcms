# Regression Panel (the release gate)

Before a release is tagged, the release candidate and the previous release genotype
the same panel of real samples on HPC, and the two outputs are compared:

- **version against version:** every changed cell is attributed to the merge (or the
  interval of merges) where it changed;
- **each version against the sign-out counts,** by stratum. The BAM is the truth, so
  disagreements are adjudicated read by read.

The panel is chosen for code-path coverage, not by hand. Every signed-out variant is
tagged with the strata that exercise distinct code paths, and a greedy set cover
picks the fewest samples that bring every stratum to its target (25 variants; 4
samples for a sample stratum), starting from the existing acceptance samples. The
selection reads the clinical sign-out dump and names samples and BAM paths, so the
selector and its output (the *tier folder*) stay on the selection host and on HPC.
The tools that run and compare a panel (`scripts/regression_panel/`) hold no patient
data.

## Strata

| Group | Strata |
|:--|:--|
| Allele shape and size | SNV, DNP, TNP, ONP; insertions and deletions of 1, 2–5, 6–19, 20–49 and 50+ bases (deletions also 100+); delins with and without a shared first base, growing or shrinking, spanning 20+ bases |
| Repeat context (every non-SNV) | homopolymer (6+), STR (3+ copies of a 2–6 bp unit), unique; an indel in a run beside another run; an insertion duplicating the reference after it (6+ bases) |
| Co-annotated structure | length-changing variants within 50 bp (tract clusters); same-position rows; an SNV between two indels, or inside an indel's window; an MNP within 10 bp of an indel; a delins overlapping a deletion or insertion row; homopolymer-twin candidates; FLT3-ITDs |
| Signal | sign-out VAF (< 2%, 2–10%, 10–40%, ≥ 40%) and depth (< 100, 100–500, 500–2,000, ≥ 2,000) |
| Loci | X, Y, MT; PMS2 (run at `--min-mapq 0`); REF == ALT and placeholder rows |
| Samples | assay (each IMPACT panel, IMPACT-HEME, ACCESS duplex + simplex), MSI-high, TMB ≥ 20, ACCESS indels in homopolymers (duplex-masked bases) |

## Arms and configurations

| Arm | Runs |
|:--|:--|
| `dna` | every panel BAM at production defaults (MAF out, PairHMM, MAPQ 20); ACCESS as duplex and simplex |
| `mfsd` | the ACCESS BAMs with `--mfsd --mfsd-parquet --observations-parquet` |
| `mq0` | the PMS2 rows of samples that have them, at `--min-mapq 0` |
| `sw` | ~25 samples chosen for indel strata, `--alignment-backend sw` |
| `vcf` | a few samples, `--format vcf` |
| `normal` | matched normals genotyped on their tumor's variants (fillout) |
| `cohort` | a few BAMs genotyped on one MAF holding several samples' variants, a variant shared by two samples listed twice (fillout-style input) |
| `rna` | the FORTE RNA truth cohort, RNA mode with the Ensembl GTF (hg38) |

## The tier folder

What a selector hands to HPC (any selector that writes this works):

- `runs.tsv`, one gbcms call per row: `run_id`, `tag`, `arm`, `mode` (`dna`/`rna`),
  `root` (`dmp`: b37 BAMs; `forte`: hg38), `bam_relpath` (under the root),
  `variants` (relative to the tier folder), `extra_args`.
- `mafs/<tag>.maf`: the sample's signed-out variants, with the sign-out counts as
  `signout_t_ref_count`, `signout_t_alt_count`, `signout_n_ref_count`,
  `signout_n_alt_count` (never a gbcms output name, so no `--column-prefix` can
  overwrite them) and the variant's strata as `panel_strata` (`;`-separated). gbcms
  passes both through to its MAF output.
- `coverage.tsv` (every stratum: available, target, covered) and `panel.tsv`.

The strata approximate the code paths they name (a tract cluster is length-changing
variants within 50 bp, not the engine's own grouping); they choose samples, and the
comparison then judges every row whatever its strata.

## Running a build

`run_panel.sh` runs one build over a run list. A run that exits 0 leaves
`<run>.ok`; re-running skips those and retries every other run, so a run that failed
after writing its MAF is not taken as finished. Every attempt is a row of
`times.<shard>.tsv` (exit code, wall time, peak memory, end time), and a run's latest
attempt is its outcome. `submit_panel.sbatch` runs it as a SLURM array (0-based; edit
its CONFIG block):

```bash
sbatch scripts/regression_panel/submit_panel.sbatch   # BUILD=6.5.0, the release image
sbatch scripts/regression_panel/submit_panel.sbatch   # BUILD=6.6.0, the candidate
```

The previous release runs from its image
(`singularity exec -B <every path the run reads> docker://ghcr.io/msk-access/gbcms:X.Y.Z gbcms`).
The candidate runs from its checkpoint wheel, the last row of `checkpoints.tsv`
(below): the commit the release ships.

## Comparing

```bash
python scripts/regression_panel/compare_panel.py TIER OUT 6.5.0 6.6.0
```

matches rows by locus (and occurrence, for a locus listed twice), never by position,
leaves failed runs out (their latest attempt did not exit 0) and writes `gate.txt`
and the reports it summarises: `version_summary.tsv` (with rows or records in one
version only),
`version_cells.tsv`, `version_rows.tsv` (every changed cell; VCF output field by
field: `FILTER`, `INFO:<id>`, `FORMAT:<id>`), `header_diff.tsv` (columns, or INFO and
FORMAT IDs, only in one version), `concordance_by_stratum.tsv` (each version against
the sign-out ALT count, per stratum: median delta and the fraction within
max(2 reads, 10%), at read and fragment level and at the **matched** level, the one
the sign-out counted: reads for IMPACT, duplex + simplex fragments for ACCESS),
`discordant.tsv` (the largest disagreements at the matched level, and the top 20 per
stratum), `normals.tsv`, `times.tsv` and `parquet.tsv` (needs `pyarrow`). PMS2's
concordance at `--min-mapq 0` comes from the mq0 arm (strata prefixed `mq0:`); an
ACCESS sample missing one flavour is left out of concordance and counted. `gate.txt`
ends with `GATE: clean` or `GATE: ATTENTION` and the reasons (failed or missing runs, a
stratum that fell more than 5 points, an arm 1.5× slower or larger), and the script
then exits 1. The rna arm's FORTE MAFs carry no sign-out counts, so it is compared
version against version only; the vcf arm's counts are the dna arm's.

## Attributing every changed cell

A release changes many cells at once (6.6.0's read-judgment work changed every
indel row of the acceptance set), so attribution is mechanical: one build per
count-affecting merge since the previous release, run on the changed rows only.

1. **The checkpoints.** `checkpoints.py REPO 6.5.0 develop > checkpoints.tsv` lists
   one build per commit on develop (in practice each merged PR) that touches counting
   code, skipping those whose own acceptance was byte-identical (`NEUTRAL` in the
   script). Each row names every commit in its interval, so a change landing in a
   skipped one is still named. Row 0 is the previous release built as a wheel, so a
   difference between its image and a wheel build is named as such; the last row is
   the commit the release ships.
2. **The wheels.** One manylinux wheel per commit, built as the release builds it
   (the *Checkpoint wheels* workflow; each wheel is an artifact named `wheel-<sha>`).
   Before a release the workflow is not yet on the default branch, where GitHub
   requires a manually run workflow to be, so start it with a push: a branch
   `checkpoint-wheels/X.Y.Z` that commits the list as
   `scripts/regression_panel/checkpoints.tsv`. Once it is on the default branch,
   *Run workflow* with the `sha` column works too. Each artifact carries its commit's
   `requirements.lock` when the commit has one; install with it, so every checkpoint
   gets the dependencies its commit pinned:
   `python3.11 -m venv cp_<sha> && cp_<sha>/bin/pip install --require-hashes -r requirements.lock && cp_<sha>/bin/pip install --no-deps gbcms-*.whl`.
3. **The reduced runs.** `attribute.py prepare TIER OUT 6.5.0 6.6.0 ATTRIB` writes, for
   every run with changes (a changed cell, or a row in one version only), a variant file
   holding the changed rows and every row within 100 bp of a kept row's span, closed
   transitively (siblings classify together), and `ATTRIB/runs.tsv`. The mfsd and rna
   arms keep their full variant files: their Benjamini–Hochberg q-values are computed
   across the run's rows. Like `compare_panel.py`, every step reads the builds' time
   records. It leaves out, and names, a run whose latest attempt failed in either
   version: a crash partway through writing leaves a partial MAF that would otherwise read
   as changed rows. No manual setting aside is needed.
4. **Each build on them**, into `ATTRIB/out`: the two versions first, then every
   checkpoint as `cp_<sha>` (`run_panel.sh TIER ATTRIB/runs.tsv cp_<sha> cp_<sha>/bin/gbcms … ATTRIB/out`).
5. **`attribute.py check`** confirms the reduction changed nothing: every reduced row
   equals the full run's, in both versions. `gbcms_status_reason` is exempt on context
   rows only (it names the co-annotation group, which the reduction may cut at its
   edge), never on a changed row. It exits 1 when a row differs, or when a reduced run's
   latest attempt failed (it names the run to re-submit).
6. **`attribute.py attribute TIER OUT 6.5.0 6.6.0 ATTRIB checkpoints.tsv`** traces each
   changed cell through the checkpoints: `attribution_cells.tsv` (each cell and the
   intervals where it changed; a row in one version only is traced as its presence,
   column `(row)`), `attribution_summary.tsv` (per checkpoint, arm and column),
   `unattributed.tsv` (cells whose trail misses a checkpoint's output or does not end at
   the candidate's value). A reduced run whose latest attempt failed counts as missing
   output, whatever its MAF holds, and is named to re-submit. It exits 1 when any cell
   is unattributed.

## The gate

The candidate is tagged when:

- **no run fails** in either build, the reduction check is exact, and **every changed
  cell is attributed** (`unattributed.tsv` empty);
- **concordance by stratum** is reported for both versions; any stratum whose fraction
  within max(2 reads, 10%) of the sign-out (matched level) falls by more than 5 points,
  and the top discordant rows, are adjudicated read by read (the read census, or the
  reads themselves) before the tag;
- **run time and peak memory** per arm stay within 1.5× of the previous release, or
  the excess is explained (wall time summed over the runs paired across the builds that
  took 10 s or more, at least 3 per arm; peak memory over paired runs).
