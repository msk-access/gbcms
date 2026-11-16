# CONTINUITY — where we left off

> Tactical state that must survive a closed laptop or a context summary.
> Update the **Now** and **Next** sections as work progresses. When a "Now" section is finished,
> move it to `.agents/archive/CONTINUITY_ARCHIVE.md` (older sections live there).

_Last updated: 2026-11-16_

## Now — #335, exact speed-ups (Stage 3; plan `docs/harness/plans/335-exact-speedups.md`)

Started 2026-11-14 on develop 1b11ffea (#438 merged: the judge fault fix RJ-35; #435 closed). The rule
(operator, 2026-11-14): speed is never bought with accuracy; a change is kept only if every output field
of every row is identical on every set, and the whole `--verbose` log too (from 2026-11-15).
- **Merged:**
  - S1, build settings (#439, develop 1992fcb4): LTO, one code unit, `libdeflate`, the Rust log cap held to
    each counting call. CPU over develop's on HPC (pf75, pf76): the rules path 0.52 to 0.78, the judge path
    0.76 to 0.96, RNA 0.67 to 0.73.
  - S4, RNA and small items (#440, 8a3966d8; its fix #441, 473e8920): the GTF and the editing database read
    once per run, keyed by the name as given and the resolved file; no insert size when mFSD is off; the
    sampled base qualities tallied in pieces; `import gbcms` without polars. RNA CPU 0.87 to 0.91; DNA
    within noise.
- **Dropped:** S1b (mimalloc): 0 rows differ, but 1 judge run in about 3,500 ended with a segmentation
  fault at exit.
- **Left:** S2 (shared scans), S3 (the judge's work done once), S5 (a narrower fetch per bin), and the
  judge's kernel (its own plan). The budget is spent: builds 8 of 8, HPC batches 3 of 3; the one rerun
  (0 of 1) is held for M2's two IMPACT cells on Gold 6348 nodes. More needs the operator's approval.
- **The floor, with S1** (pf73 times pf75, about +/-10%): the rules path 0.25 to 0.74 of the C++ tool's
  wall time; the judge path under 1 only at 1,000 rows x 1 BAM, 1.14 at IMPACT 5,000 x 1, 2.2 to 5.0 at
  the deep ACCESS and the 9- and 10-BAM cells. M2: the read call is 44 to 69% of the judge's CPU; M3: reads
  with best score 0 hold 85 to 87% of the scoring clock.
- **Harness:** the gate `batch680/gate335.py` (18 arms with three log arms and the RNA fill-out; baseline
  `fin682`, develop 1b11ffea); the frozen builds `src_s1c`, `src_s1d`, `src_s4b`, `src_s4c` and others.
  After #335 closes, the retirement checklist (`docs/harness/measurement-harness.md`).
- **Next (operator, 2026-11-16):** plan S2, S3 and S5 "meticulously including research and solutions for
  the source and class". The draft: `docs/harness/plans/335-phase2.md` (research in the session's
  `scratchpad/speed2/plan2/`). Approved 2026-11-16 11:14 as drafted (PR #442); D2 decided 11:50, (a):
  accept S5's log change; D3 (F1, F2 tickets) open. Work order: S2 with V (gate 0 differences on all
  18 arms; its PR), S3a, M4, S3b, S3c, S5. M0 cleared at 10:14 (ERR-20261116-001). M5 done (local shares, `plan2/measure/MEASURE.md`): the
  site build is 45-58% of the judge here (14% at ACCESS on HPC). M6 (pf75): one sample keeps its 4 CPUs
  86-98% busy on the judge path, 53-94% on the rules path. When S2 is done: move the scratchpad's
  `dev473/` (an unused develop build, 1.3 GB) to the Trash (operator 12:22: "Clean up when its done").
- **Two findings outside #335** (from the code, not reproduced; own tickets): F1, the junction judge's
  scan reads past the one-row geometry's fetch (invariant 1, test geometry only); F2, the verdict cache's
  key leaves out qualities, so a rare duplicated record can get another call between runs.

## Now — 6.8.0 starts (tracker #420)

**6.7.0 is released (2026-11-10):** tag `6.7.0` on the merge 0f5c1f0b (#426); PyPI `gbcms 6.7.0`
(abi3 wheels for Linux x86_64 and Apple silicon, and the sdist); `ghcr.io/msk-access/gbcms:6.7.0`
(amd64 and arm64); the GitHub Release. #196 and the 6.7.0 milestone are closed. develop is at
6.8.0-dev.0. The 6.7.0 record (the gate, the decisions, the release) is in the archive.

**The theme (#420):** make the sequence read judge the default (#331), meticulously and
systematically. The operator decides the order at the start: count changes first, then speed and
structure once the counts settle. The judge meets the speed floor before it is the default (#335).

**First task — the judge's read-by-read check against 6.6.0** (the 6.7.0 release panel; CHANGELOG
6.7.0, "Release gate"). With `--read-judge sequence`, 456 rows flip against 6.6.0: 408 leave the
sign-out margin, 348 of them above the sign-out (138 DNP rows).
- Done: the DNP rows above the sign-out are a sensitivity gain, not a fault. The judge credits 1,544
  reads that 6.6.0 did not; 1,543 of them carry the DNP only in the soft clip, and 97% carry both ALT
  bases at Q20 or above (RJ-13: soft-clipped bases are the read's own).
- Done (2026-11-10 to 2026-11-12; harness notes `judge680/rbr/STATE.md`): the 58 SNV rows above
  the sign-out (57 agree with the judge, 1 ties); the 60 rows below it (37 SNV, 23 other). Merged:
  #430 (RJ-29: a realigner's placeholder tags are no evidence of the aligner's doubt), #432 (RJ-33
  and RJ-34 in the register), #433 (RJ-34: `TIED_WITH`, and the found alleles a listed row
  explains; ticket #431; plan `docs/harness/plans/tied-listed-reads.md`).
- Left, in this order: the judge faults at 2 places (equal REF and ALT scores where the bases
  differ; RJ-34's measure); the 152 rows above the sign-out that are not SNV or DNP (DEL 64,
  INS 53, MNP3-11 22, DELINS 13); the 48 rows the judge brings into the margin. Attribute each
  class read by read before any rule change (`.agents/rules/decision-protocol.md`).
- Acceptance rows for #331: the three shapes the rules count as 0 and the judge counts (6.7.0
  "Known issues": #417, the tandem duplication through #417, and #425).

**Harness (local, `~/test/gbcms/harness/plan670/`):**
- `gate_local/`: the release panel and its scripts (`judge_check.py`, `judge_verdicts.py`,
  `mnp_bases.py`, `strata3.py`, `class_map.py`, `class_split.py`, `zero_alt_*.py`,
  `run_judge670.sh`, which builds the judge arm).
- `m294/oracle_sibling_fast.py`: the oracle (BWA-MEM against REF, ALT and sibling haplotypes), with
  the arms `flip` and `flipj`.
- `r2/plain_trace.py` for traces (`docs/harness/measurement-harness.md`, "Traces").
- `gates/gate_wt.sh WORKTREE VENV LOGDIR`: the whole qa-check gate, one exit code per tool.
- `judge680/`: the tool survey for #417 and #425.
- `batch680/`: develop 44b02197's 18 arms over the tier, ACCESS, WES and HG002 (whole logs kept),
  and each later ticket's arms beside them (`README.md`; first reads in `READS.tsv`).
- `judge680/rbr/`: the read-by-read notes (`STATE.md`), labels and scripts; every patient-derived
  print goes through `./san`.
- Frozen builds: `src_<label>/.venv` holds an installed wheel (for example `src_v660` for 6.6.0 and
  `src_rc670d` for the candidate f7293618). `build_frozen.sh` makes them from checkouts in
  `frozen_src/`; the panel scripts at f7293618 are in `frozen_src/rc670d/`. `run_v660.sh` and
  `build_proto_*.sh` are finished one-off runs that still name the 6.7.0 session's scratchpad: do
  not rerun them as they are.
- The PHI rules hold: the patient files (for example `zero_alt_map.tsv`, `flips_judge.tsv`) stay
  local; never list a gate folder; take paths from the manifests inside a sanitized script.

**Plan-first:** the 6.7.0 trial did not meet its bar (`docs/harness/plan-first.md`, "Result").
Operator 2026-11-06: try it again on the 6.8.0 work, with a new bar fixed before the first
count-changing build.

**At the 6.8.0 release:** the app's auto-mode check stops an agent's push of a release tag (the push
publishes to PyPI and GHCR). The agent tags the merge commit locally and runs
`scripts/release.py check --tag X.Y.Z`; the operator pushes the tag; the agent then checks the
workflow, PyPI, GHCR and the Release page (release guide, steps 8 and 9).

**Open:** the guard's masks depend on `--min-baseq` (#374); the SNV tie rule with a confident third
base (#367); Nextflow resources wait for the sweep (branch `feature/nextflow-resources`, local).
#420's sub-issue list does not name #425 yet.

**Housekeeping (operator approval first):** 37 worktrees under `.claude/worktrees/` and
`../py-gbcms-d4`, most on merged branches: list them, and retire the merged ones to the Trash
(never `rm -rf`), then `git worktree prune`. The `release/6.7.0` branch on GitHub (deletion is
optional, release guide step 11).

## Next — 6.8.0 (tracker #420; the order decided 2026-11-12)

The plan page: `docs/harness/plans/420-cycle-order.md` (the measurements, the options, the stages).
- **Decided (operator, 2026-11-12: "Ok lets do meticulously"):** D1-A, the judge fault first;
  D2-B, measure and design speed now and build it after the counts settle; D3-B, DNA flips in
  6.8.0 and RNA stays on the rules path until #333 (6.9.0); D4-A, 6.8.0 keeps the path to the
  default (34 tickets, plus #364, #397, #401 and the fault ticket #435) and 24 tickets move to 6.9.0.
- **Stage 0, housekeeping (posted 2026-11-13):** #429 closed; #425 added to #420; the milestone
  moves; #331's "done when" item 3 (RNA) moved with #333 under #328; the fault ticket #435. Open:
  Greptile's review credits (the operator) and the old-worktree script.
- **Stage 1, measurement only (started 2026-11-12):** the fault's mechanism (1.1); the 13 open
  rows of RJ-34's bar, then the 152 and the 48 (1.2); J1's 98 DNA failures sorted (1.3); speed
  sizing: a local HG002 profile, a survey, one HPC batch (1.4). Done: the survey (2026-11-12); the
  batch pf73 (2026-11-13; build dv680, develop 3be1765b): the judge at 0.76 to 5.27 of the C++
  time, the floor met only at 1,000 rows; at ACCESS 20,000 × 1 its extra time over the rules path
  must fall from 200 s to 3 s (plan page, 1.4). Open: the local HG002 profile; the scoring call's
  share of the judge on the 20,000-row cell (a profile of that cell).
- **Stage 2, the count decisions:** the fault fix (the trial's item 2); #405 (item 3); #401 with
  #306; the MNV and delins family (#410, #364, #414, #397); #418; the new classes; a baseline batch.
- **Stage 3, speed** (#335 with #150, #312, #412, #406). **Stage 4, the default** (J1 and J3 to
  J6: #332, #334 to #337; J2, RNA, is #333 in 6.9.0). **Stage 5, the release** (#318, the floor, the gate, the tag).
