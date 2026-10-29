# Measurement harness: layout and conventions

The measurement harness is where count-changing work is measured before it merges (the
method is `.agents/rules/decision-protocol.md`). It lives on the operator's Mac at
`~/test/gbcms/harness/`, outside the repository, because most of its inputs are patient BAM
slices. This page records its **structure and rules only**: no sample IDs, coordinates or
file contents. Describe a file by its pattern, never by a patient file name.

- Patient slices (panel, ACCESS, WES, normals, the gate) may be worked on and tested with. Nothing
  from them is used downstream. Everything a tool prints goes to the model provider, so report
  aggregates and allele shapes, never sample IDs or read names (`.agents/rules/security.md`).
- `cbio/` holds cBioPortal-derived work and is protected: no tool reads, lists or searches it
  (`.agents/rules/security.md`, hook `protect-cbioportal.py`).
- HG002 (`giab/`) is public: its numbers may go to GitHub.

## Top level

| Path | What it is |
|:--|:--|
| `plan670/` | The 6.7.0 plan's harness: nearly all current work (below). |
| `regression_panel/` | The D5 panel kit (`kit_660/`: `RUNBOOK.md`, `tier_660/` run list, scripts) and its run, compare and attribution scripts; the repo's `scripts/regression_panel/` is the committed part. |
| `cbio/` | Protected cBioPortal work; the one console script is `cbio_summary.py` (schema and counts). |
| `obs_ab.py`, `REMOVED_BUILDS.md` | Shared constants (`FA`, the b37 FASTA path) and the log of retired builds (each one's record is kept in `plan670/patches/<tree>/`). |
| `c<NN>/`, `g<N>/`, `t<NN>/`, `h<N>/`, `s<N>/`, `r<N>/`, `d<N>/`, `rc650/`, `spikein/`, `tempo/` | Older work, named for the plan item (`c37`), group (`g3`, RNA), ticket (`t110`) or release candidate it served. |

## Arms in `plan670/`

An arm is one input set a build is run on. Each arm directory holds one `out_<label>/` per build run.

| Directory | Arm | Reference FASTA |
|:--|:--|:--|
| `slices/` | Panel slices (IMPACT and ACCESS) for the plan rows, one `<run_id>.bam` per run of the 6.6.0 gate tier (`regression_panel/kit_660/tier_660/`). | b37 |
| `hpc_kit_670/` | The HPC slicing kit (README, region lists, `slice.sh`) and what it brought back: `slices_indels` (the gate panel's indel, MNP and delins rows), `slices_all`, `slices_std` (ACCESS standard BAMs), `slices_forte` (RNA). | b37; `slices_forte` GRCh38 |
| `gate_pi/`, `gate_cx/` | The gate panel's pure-indel rows and complex rows (MNPs, delins): one run per sliced BAM, with every panel variant within 300 bp (co-annotated siblings change counts). | b37 |
| `gate_local/` | The full local gate: run scripts per build, and `attrib/` (every changed cell attributed to the merge where it changed). | b37 |
| `across/` | Other data types: matched normals (hom-REF controls for somatic rows, hets for germline ones), ACCESS duplex and simplex, IMPACT cohort runs. | b37 |
| `germ/` | Germline rows genotyped in their normal slices. | b37 |
| `wes/` | WES (the TEMPO set; BWA only, no realigner) at the panel's long-event loci, with its BAM slices in `slices/`: the control for any "where did the aligner put the gap" question. | b37 |
| `giab/` | HG002 truth arm (public), by event class: `pri*` (priority sites), `c246` (adjacent pairs), `tr` (tandem repeats), `sv`. | GRCh38 |
| `c43syn/`, `c43_rna/` | Synthetic sets (they must mirror a real library: both pair orientations, read errors, low and normal depth) and the RNA (FORTE) slices. | synthetic; RNA GRCh38 |
| `tools/` | Ten other tools on the same slices (below). | per arm |
| `r2/` | The read-judge bake-off: HG002 truth builders, `evaluate.py`, `plain_trace.py`, and the `hapjudge` prototype. | GRCh38 |
| `perf/`, `rjspec/`, `survey/` | Speed measurements; the read-judgment spec cases (synthetic, PHI-free) run through other tools, with `reading.md` and `compare.md`; the tool and literature surveys S1 to S5 (markdown). | |
| `m<issue>/` | Measurements for one issue, named by its number (`m294` holds #294's oracles). | |
| `sp_tools/` | A mirror of the session scratch scripts (scripts only; paths inside still point at the old scratchpad, and outputs were not copied): `c<NN>/freeze.sh` (the build scripts), `retire_builds.sh`, `lane_slot.sh`, the patient oracles (`xdata/`, `oraclev2/`, `oraclev3/`), the truth scripts (`truthv2/`, `v2d/`) and small probes (`bq/`, `mnv/`, `phi/`). | |
| `parked_<issue>/`, `set_aside_<date>/` | A parked prototype (patch and README), and outputs set aside. | |

A measurement run against the wrong FASTA fails REF checks or silently skips rows.

- **b37:** `Homo_sapiens_assembly19.fasta` (contigs `1`, `2`, ...) in the OneDrive
  `nucleo_qc_generation/test_data` folder; the arm scripts carry the path as `obs_ab.FA`.
- **GRCh38:** `~/Downloads/genome.fa` (contigs without `chr`; HG002 rows say `chr10`, and gbcms
  reconciles `chr1` with `1`); RNA uses the Ensembl 111 GTF.

Name the b37 FASTA in every brief for a panel arm: an agent once defaulted to `genome.fa`.

## Builds, outputs and labels

A **label** names one build: the plan item plus an iteration (`c53d` is the fourth build of item
C53), or a release or commit (`v660`, `dev<hash prefix>`). Everything is keyed by it.

| Pattern | Meaning |
|:--|:--|
| `src_<label>/` | A frozen build, made by the package's freeze script (`sp_tools/c<NN>/freeze.sh <label> <commit> "<note>"`): a `git archive` extraction of the commit in `src/`, a release wheel in `dist/`, `build.log` and its own `.venv/`. It is not a registered git worktree, so a build needs no `git worktree remove`. About 200 MB each (the cargo target is kept in the scratchpad, not in the tree). Older builds (`plan670/build_frozen.sh`, which built from a temporary worktree) have no `src/`: only the wheel, the venv and the record. |
| `src_<label>/BUILT_FROM` | The commit it was built from (the first word is the SHA) and a note. A prototype with local edits records its base commit and keeps a `<label>.patch` in the tree. |
| `out_<label>/` | The outputs of running that build: per run, a directory with the MAF, and `<run_id>.trace.gz` (the `read call` lines only). A label ending `nt` is an untraced run. |
| `accept_<label>.sh` | One ticket's acceptance runner. `accept_c39.sh` is the template (panel, germline, HG002, the oracles, the summaries); `accept_one.sh <label>` runs a single build. |
| `run*.py <label> <gbcms> [parallel]`, `*_summary.py <labels>`, `oracle*.py <label>`, `... --diff <A> <B>` | The runner, summary and comparison scripts. A runner skips a run whose done marker exists (`<run_id>.trace.gz`, or a `DONE` file beside the MAF) and reruns one without it; move a finished run's outputs away to redo it. |

Build and test each candidate in isolation (`git archive <sha>` and its own venv: the
rust-python-ffi skill), and record the commit beside the build. Never measure with the shared
`.venv`.

**Lane slots.** At most two harness lanes run at once on the machine. `sp_tools/lane_slot.sh
<command>` waits for a free slot (`plan670/.slots/slot1`, `slot2`, held with `lockf`) and runs
the command; `--exclusive` holds both, for wall-clock timing runs.

## Oracles and truths

A count is judged by the reads' own bases, never by the aligner's placement
(`genotyper-not-caller` memory). Each oracle below trusts no alignment, but check that it does
not share the rule under test.

| Oracle | Use |
|:--|:--|
| Bases-only oracle (`oracle.py`, `oracle_c40.py`) | Judges each read by its own bases against REF and ALT, whatever the aligner did. Alignment scores are rescored over the whole read from the CIGAR (`r2/bwa_rescore.py`, shared with the patient oracle): BWA's local `AS` leaves a read-end mismatch uncharged. |
| Patient oracle (v3: `sp_tools/oraclev3/xoracle_v3.py truth <arm> [--min-baseq N]`; it imports v2, `sp_tools/xdata/xoracle_v2.py`) | Per-haplotype indices and an exact DP; a label needs both. Must take the run's `--min-baseq`, never the default. |
| Cluster-aware oracle (`m294/oracle_sibling.py`) | For any change that can move a read between co-annotated rows: scores against REF, the row's ALT and every sibling's ALT (rows of the same sample within 30 bp). A REF-versus-ALT check (`oracle_changed.py`, right for isolated rows) calls a removed sibling molecule "lost": on #300 it reported 157 panel and 204 gate reads lost, and the sibling-aware check found 556 removals were the sibling's molecules or ties, and only 17 gained calls wrong. `oracle_sibling.py` imports `oracle_changed.py`, which runs its own check at import, so a run takes twice as long. A read whose sibling allele contains this row's change is this row's ALT (a sibling MNV that contains the row's SNV). |
| HG002 truth (`r2/truth_v2.py --tag v2d`, `truth_v2d_bases.py`) | See the scoring rules below. |
| pysam oracle | Cross-check counts with `bam.fetch()` and `read.get_reference_positions(full_length=True)`, not `pileup()`, which under-counted depth by about 10% on a duplex BAM even with `max_depth` raised. Match the run's filters (MAPQ, BQ, flags `0x4 0x100 0x200 0x400 0x800`; improper pair and indel are not filtered). `RD` and `AD` match exactly; `DP = RD + AD + neither`, so DP is higher (571 SNPs on a ctDNA duplex BAM: ALT 571 of 571 exact, REF within 1%, DP within 1% for 90%, the rest being neither reads). |
| Read census (`tests/census.py`, in the repo) | Mirrors the decided read-input rules (fragment end, filters), so it cannot validate them. Check such a rule against evidence it does not share: the bases in question (an adapter motif at the boundary, the mismatch rate to the reference), the mate's own alignment, and per-read traces at rows where engine and census move differently, ITD and indel rows first. |
| vcf2maf oracle | VCF to MAF representation is "what mskcc/vcf2maf writes" (#110). `vcf2maf.pl` and `maf2vcf.pl` live at `~/Documents/Github/vcf2maf` and need `samtools` on `PATH`. Its quirks are not gbcms bugs: `maf2vcf` dies on an allele outside `[ACGT-]`, writes an empty ALT for REF equal to ALT (gbcms FAILs it as `ALT_EQUALS_REF`), and sends rows whose REF mismatches the FASTA to a `.skipped.tsv`; `vcf2maf` trims case-sensitively, gbcms case-insensitively; `maf2vcf` writes a differing `Tumor_Seq_Allele1` as a second ALT, gbcms does not. Acceptance (2026-09-24, harness `t110/`): about 102k unique real sign-out shapes matched both oracles exactly except the two REF-equals-ALT rows. |

**HG002 truth arm** (public; `giab/`, `r2/`). Truth is HG002 v5.0q (GRCh38; NIST
`release/AshkenazimTrio/HG002_NA24385_son/v5.0q/`): single-ALT phased hets inside the benchmark
BED, no other truth variant within 30 bp; dipcall splits a delins into nearby records, so
same-phase records within 10 bp (one an indel) are merged into one delins. Reads are the
DeepVariant case-study BAM (bwa-mem, no realignment, 151 bp), fetched as region slices with
`samtools view -M -L`. Fetching is latency-bound (about 1.5 s to the first byte, 20 KB/s per
connection): fetch in parallel chunks (48 chunks did 700 sites in minutes); too many connections
at once can fail with "Destination address required", so retry the failed set. Run without
`--trace` for truth counts; the bases-only oracle does the per-read work.

Scoring rules (operator, 2026-10-26; use `--truth v2dbq20` for every HG002 bar):
- Bars count reads decidable by their bases at the run's `--min-baseq` (the threshold comes from
  the run, never hard-coded). A read that fits its truth sequence and another candidate equally is
  not decidable by bases, and reads only phasing can decide leave the bars.
- Homozygous rows are compared with REF.
- Bars count rows inside the GIAB v5.0q benchmark, picked by variant size: the smvar BED below
  50 bp, the stvar BED at 50 bp or more (the smvar BED excludes every 50 bp+ variant).
- Long deletions: count the reads holding the anchor, plus reads placed inside the deleted bases
  whose own bases hold the junction.
- A class fails only past 0.5 points and 2 reads.
- Truth v2 aligned each row on its own (rows 1 bp apart share contigs in one combined index, so a
  missing alignment could decide a label: about 9% of the adjacent-pair het labels); v2c rescored
  from the CIGAR; v2d added the homozygous and benchmark rules. Old-truth numbers are build-to-build comparisons, not
  absolutes. The truth scripts read the default `--min-baseq` from `src/gbcms/cli.py`
  (`cli_default_min_baseq`), so their pattern must allow the option to span lines.
- HG002 is one genome, one aligner, germline hets: it cannot show false ALT, or how ABRA2, STAR or
  cfDNA reads behave (`.agents/rules/decision-protocol.md`).

## Traces

`--trace` logs one `read call CHR:POS REF>ALT read=... mate=... ref=... alt=... phase=...
partial=N` line per classified read, through a Rich console handler that wraps and pads every
line at the terminal width. To trace, run the CLI through `r2/plain_trace.py` (a plain
`%(message)s` handler: nothing wraps or pads, the calls are identical, about 3x faster), as
`<build venv>/bin/python plain_trace.py dna ... --trace` or `r2/run_plain.sh ARM LABEL PYTHON`.
Plain `gbcms --trace` needs `COLUMNS=60000` in the environment (`COLUMNS=3000` still wraps
alleles past about 2.9 kb) and is hundreds of times slower than an untraced run. Assert a
non-zero count of parsed calls before trusting a diff of two traces.

- `partial=` is a count of positions matching the ALT, not a flag; `mate` is 1 for first in
  template, else 2.
- A read call is not a count: reads the engine then leaves out are logged too (at long delins the
  trace showed 2,403 REF reads where the MAF wrote 1,356). Take totals from the MAF; use the trace
  for which reads.
- The key is the prepared variant (VCF-anchored, left-aligned), not the input row: match a trace
  key to a row by the haplotype both give.

## Other tools

`tools/` (PHI for panel, WES and normals) runs ten tools in amd64 containers (`images.txt`: GATK with
HaplotypeCaller and Mutect2, bcftools, freebayes, VarDict, Octopus, bam-readcount, Paragraph,
GraphTyper, Delly, SVTyper) and
compares them with gbcms and with an independent per-read check (`readtruth.py`: BWA-MEM of every
read to the REF and ALT local haplotypes). Entry points: `make_sites.py`, `run_tools.py`,
`compare.py`, `summarize_tools.py`, `per_read.py`, `build_report.py`.

- Run **one container at a time**: five emulated containers at once hung Docker's file share
  (every container at 0% CPU; `docker desktop restart` clears it). Mount only the harness folder;
  batch per-BAM tools in one container (a start costs about 15 s).
- The IMPACT and ACCESS slices carry no proper-pair flag: bcftools skips such reads by default
  (use `-A`); gbcms does not filter on it either.
- Mutect2 (GATK 4.6.2.0) crashes on forced alleles longer than its region padding
  (`--assembly-region-padding` must cover the longest deletion; slow under emulation). Delly 2.x
  genotyping reads sites only in its own format (`SVMETHOD=EMBL.DELLY*`, `CONSENSUS`/`CONSBP`, and
  `SR`/`PE` above 0, or it writes nothing) and needs every BAM of a run to share contigs; Delly
  and SVTyper estimate insert sizes from the BAM and fail on small slices. Paragraph 2.3's VCF
  writer aborts (an htslib assertion): read `genotypes.json.gz`. GraphTyper copies output with
  permissions the share refuses: write to `/tmp` in the container. Octopus 0.7.4 hangs under
  emulation (AVX).
- GATK, freebayes and VarDict decompose a given delins: credit them with the records that
  together spell the allele (`compare.from_components`).
- The per-read check reads deletions' VAF a little low: a REF read starting just before POS holds
  deleted bases and is clearly REF; an ALT read starting just before the junction is mostly
  shared right flank and stays unclear.

## Rules for running

- **Make a long run visible.** Start it as the tracked background task itself (Bash with
  `run_in_background`; `| tee <log>` keeps a log), so the operator sees the run and is told when it
  ends. Tracked tasks ran 6+ hours on 2026-10-19. Use `nohup ... &` only for a run that must
  outlive the session, and pair it with a tracked waiter on the runner's **process**:
  `while kill -0 <pid> 2>/dev/null; do sleep 60; done; <print the log's summary>`. Never wait on a
  marker line or file the run writes at its end (a run that dies never writes it, and a restart
  that reuses the log name overwrites it), and never on `pgrep -f "<pattern>"` (the waiter's own
  command line holds the pattern, so it finds itself). One runner per log name; a run and its
  waiter stop together. When a run is slow (full `--trace` on expressed genes), say why and switch
  to a faster method.
- **Network mounts** (the DMP and FORTE SFTP mounts): runs that read one go one at a time, grouped
  by file, never in parallel. This Mac has no `timeout` or `gtimeout`: `timeout N cmd && ok || fail`
  always prints "fail" (command not found). Check a mount with a plain read
  (`head -c 60 <file> >/dev/null && echo ok`) or `ls`; `Errno 5` from Python or pysam is the real
  signal.
- **HPC (SLURM).** Never run work on the head node. Setup and analysis scripts run in an
  interactive allocation (`srun --pty bash`, or `salloc` and then check `hostname`, since a plain
  `salloc` can leave the shell on the head node); submitting with `sbatch` from the head node is
  fine. `cmobic_short`: up to 3 h, with `--qos=priority`, for interactive work
  (`--mem=40G --time=02:59:00 --partition=cmobic_short --qos=priority`). `cmobic_cpu`: 7-day limit,
  for anything longer (the regression panel's `submit_panel.sbatch` and the Nextflow retry tier use
  it). Whether compute nodes reach GitHub or PyPI is not known (2026-10-07): an offline wheelhouse
  is the fallback for installs.
- **Disk.** Frozen builds, raw traces and `rust/target` are what fill it. Check `df -h` on
  `/System/Volumes/Data` (the Data volume, not `/`) before a large run.
- **Review agents** run in the session's current directory and may change git state; see
  `.agents/rules/git-conventions.md` (shared checkout hazards).

## Retirement checklist (after each merge)

Run it once the ticket's PR is merged and its issue is closed. Acceptance builds pile up one or two
per ticket and are needed only until the merge; everything in a build is rebuildable from the
commit in its `BUILT_FROM`.

1. **Builds.** Run `sp_tools/retire_builds.sh [--trash] <tree> ...` (for example `src_c63g`). For
   each tree it keeps the record (every file except `.venv`, `dist`, `target` and `src`) in
   `plan670/patches/<tree>/`. Where the tree holds its source (`src/`), it saves that source's edits
   against the `BUILT_FROM` commit as `edits.patch` and checks that the patch re-applies; a tree
   whose patch fails is not trashed. Only with `--trash` does it log the tree in `REMOVED_BUILDS.md`
   and move it to the Trash; without `--trash` it keeps the record and checks the patch. Builds are
   not registered worktrees, so no `git worktree remove` applies to them. The disk is freed only
   once the operator empties the Trash.
2. **Traces.** `gzip` raw `--trace` logs (about 400x smaller). Keep outputs, scripts and logs: they
   are the acceptance evidence.
3. **Cargo.** `cargo clean` in `rust/` when `rust/target` exceeds a few GB (the debug profile
   bloats past 6 GB).
4. **Worktrees and branches.** Remove the merged branch's worktree (`git worktree remove`) and
   delete the merged local branch, and the detached `agent-*` worktrees that review agents leave,
   once nothing in them is unpushed (`git status`, `git log @{u}..`). Never remove one that
   another session is using.
5. **Notes.** Move the finished section of `CONTINUITY.md` to `.agents/archive/`, update its
   "Now", and check `.agents/memory/MEMORY.md` against `.agents/memory/` for a decision that
   changed.

Leave global caches (pip, uv, the cargo registry) and the app's worktree pool alone.
