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

## Running a build

`run_panel.sh` runs one build over a run list; re-running skips finished runs, and
each run's wall time and peak memory go to `times.<shard>.tsv`.
`submit_panel.sbatch` runs it as a SLURM array (edit its CONFIG block):

```bash
sbatch scripts/regression_panel/submit_panel.sbatch   # BUILD=6.5.0, the release image
sbatch scripts/regression_panel/submit_panel.sbatch   # BUILD=6.6.0, the candidate
```

The previous release runs from its image
(`singularity exec docker://ghcr.io/msk-access/gbcms:X.Y.Z gbcms`). The candidate
runs from its checkpoint wheel, the last row of `checkpoints.tsv` (below).

## Comparing

```bash
python scripts/regression_panel/compare_panel.py TIER OUT 6.5.0 6.6.0
```

writes `gate.txt` and the reports it summarises: `version_summary.tsv`,
`version_cells.tsv`, `version_rows.tsv` (every changed cell), `header_diff.tsv`,
`concordance_by_stratum.tsv` (each version against the sign-out ALT count, per
stratum, at read and fragment level: median delta and the fraction within
max(2 reads, 10%)), `discordant.tsv` (the largest disagreements, and the top 20 per
stratum), `normals.tsv`, `times.tsv` and `parquet.tsv` (needs `pyarrow`). ACCESS rows
are compared with the sign-out after summing duplex and simplex.

## Attributing every changed cell

A release changes many cells at once (6.6.0's read-judgment work changed every
indel row of the acceptance set), so attribution is mechanical: one build per
count-affecting merge since the previous release, run on the changed rows only.

1. **The checkpoints.** `checkpoints.py REPO 6.5.0 develop > checkpoints.tsv` lists
   one build per merge into develop that touches counting code, skipping merges whose
   own acceptance was byte-identical (`NEUTRAL` in the script). Each row names every
   merge in its interval, so a change landing in a skipped merge is still named.
2. **The wheels.** Actions → *Checkpoint wheels* → *Run workflow*, with the `sha`
   column of `checkpoints.tsv`: one manylinux wheel per commit, built as the release
   builds it. Install each into its own venv on HPC
   (`python3.11 -m venv cp_<sha> && cp_<sha>/bin/pip install gbcms-*.whl`).
3. **The reduced runs.** `attribute.py prepare TIER OUT 6.5.0 6.6.0 ATTRIB` writes, for
   every run with changes, a variant file holding the changed rows and every row within
   100 bp (siblings classify together), and `ATTRIB/runs.tsv`.
4. **Each build on them**, into `ATTRIB/out`: the two versions first, then every
   checkpoint as `cp_<sha>` (`run_panel.sh TIER ATTRIB/runs.tsv cp_<sha> cp_<sha>/bin/gbcms … ATTRIB/out`).
5. **`attribute.py check`** confirms the reduction changed nothing: every reduced row
   equals the full run's, in both versions.
6. **`attribute.py attribute TIER OUT 6.5.0 6.6.0 ATTRIB checkpoints.tsv`** traces each
   changed cell through the checkpoints: `attribution_cells.tsv` (each cell and the
   intervals where it changed), `attribution_summary.tsv` (per checkpoint, arm and
   column), `unattributed.tsv` (cells whose last checkpoint is not the candidate's
   value, or that never step).

## The gate

The candidate is tagged when:

- **no run fails** in either build, the reduction check is exact, and **every changed
  cell is attributed** (`unattributed.tsv` empty);
- **concordance by stratum** is reported for both versions; any stratum whose fraction
  within max(2 reads, 10%) of the sign-out falls by more than 5 points, and the top
  discordant rows, are adjudicated read by read (the read census, or the reads
  themselves) before the tag;
- **run time and peak memory** per arm stay within 1.5× of the previous release, or
  the excess is explained.
