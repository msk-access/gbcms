# CONTINUITY — where we left off

> Tactical state that must survive a closed laptop or a context summary.
> Update the **Now** and **Next** sections as work progresses. When a "Now" section is finished,
> move it to `.agents/archive/CONTINUITY_ARCHIVE.md` (older sections live there).

_Last updated: 2026-11-12_

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

## Next — 6.8.0 (tracker #420)

- The judge becomes the default (#331), with #405, #406, #410, #414, #417, #418 and #425; the judge
  meets the speed floor first (#335).
- Packages 4 to 7: evidence outside the aligned bases (#259 with #177, #179, #212); guard
  statistics and artifacts (#306, #180); retiring the homopolymer twin (#261 with #111, #112,
  #128); complex-event read-on (#274, #361). Also #290, #293, #178, #408 and #409.
- Moved from 6.7.0 on 2026-11-06: speed (#412, #312, #150); clean-up (#313 with #195, #314, #315,
  #316); the truth arms and the release gate (#260 with #266, #267, #268, #275); output and
  hygiene (#273, #132, #127, #138); the pre-release review (#318).
- The order is decided with the operator when 6.8.0 starts: count changes first, then speed and
  structure once the counts settle.
