# CONTINUITY — where we left off

> Tactical state that must survive a closed laptop or a context summary.
> Update the **Now** and **Next** sections as work progresses. When a "Now" section is finished,
> move it to `.agents/archive/CONTINUITY_ARCHIVE.md` (older sections live there).

_Last updated: 2026-11-03_

## Now — 6.7.0: the count-changing class plans, one at a time (2026-11-03)

**Merged to develop (2026-11-02/03):**
- #297, the unit-run rule (rule X: a listed allele is always a candidate) with Greptile's three fixes (#404).
- #264, long-insertion windows (#396), and its follow-up (#413, judge path only):
  - one mask copy per read;
  - read-on by path columns;
  - the tie window check;
  - the full class rule: a change a read cannot place is a plain alternative, and copies are made only for the changes its best path spans;
  - the kept graph sets are charged to the cache budget.
  - z68 measured it: every HG002 class equals develop, the patient views pass 24 of 24, and the rules path is byte-identical.
- RNA #292 option B (operator 2026-11-02/03): no junction penalty on any path (RJ-31; it reverses the penalty of 2026-05-04), with the decided partial-span rule RJ-17 (#415).
  - Measured on FORTE: delins at bar 30 go from 280 to 173 (p292d). The reads losing ALT are spliced, with the event 0-4 bases from the junction, and their own bases read REF.
  - DNA mode is byte-identical.
- CI: macOS on `macos-latest` (#399); the macOS wheel's minimum version is checked on the tag and on each binary (#407).
- Agent rules: read the bot review before ready or merge (#403); the gate covers only the tree it ran on, and cherry-pick and merge commit too; never list a patient folder; own-gap for every state is rejected (#411).
- #402 (the PHI read guard) is closed unmerged. It is being rebuilt as a separate, tool-agnostic project (agent-guards; DESIGN.md and KICKOFF.md in `~/Documents/Github/agent-guards`, the operator's own session).

**In progress: the class plan for site candidates and REF** (`docs/harness/plans/site-candidates-ref.md` on `chore/plan-first`):
- Operator 2026-11-02: option B, "tell what the read reaches". It reverses RJ-25's REF clause.
- Branch `fix/ref-what-the-read-reaches` (local, f3c6565e), on both paths:
  - every listed allele at the site is a candidate (no floor);
  - REF where the read tells REF from every candidate it reaches;
  - the delins core by REF bases replaced (50 or more).
- Rules path only: SNV rows with a listed indel or delins at the site are judged by the sequence judge (#398; operator 2026-11-03).
  - It removes 1,230 patient wrong-ALT reads.
  - Speed is accepted for now and ticketed (#412). Patient lanes run 1.12x slower.
- rq1 failed three rules-path bars. Two causes were found and fixed: the candidate test now reads the event's own flank, and reach is measured on the allele's unmasked difference. The third, the delins gap between the paths, is older than B and is ticketed for #331 (#414).
- **Running on HPC:** rq2 (6731328; 6731367, 6731384). Next steps:
  1. copy back rq2R and rq2S, and score them against n105 and rq1. The bars: the panel INS-1 class and HG002 INS 1-9 on the rules path; 24 of 24 views on the judge path; no wrong-ALT;
  2. a fresh review;
  3. the PR. Close #394 (superseded by B) and #398 with it.

**Plan-first trial:** branch `chore/plan-first` (local) holds the template and the plans (audit, #297 close, C2, #292, site candidates). The measurement records are in `~/test/gbcms/harness/plan670/plan-first/`. Next: write the trial's result against its baseline bar, then open its PR.

**Tickets opened 2026-11-02/03:**
- #405: DP for off-anchor reads.
- #406: the mask-copy cap per state graph, with the measured memory.
- #408: at `--min-baseq` 0, a base with no evidence counts.
- #409: a piece too short to place its own bases, at a read's end or between two junctions (RNA B's limit).
- #410: a listed MNV and an insertion carried together.
- #412: the routing's speed.
- #414: the delins gap between the paths.

**Cleanup waiting for the operator:** the worktrees and builds of the merged PRs (macfloor, rules1103, rna292b, f264b) and the RNA measurement builds, to the Trash, by the retirement checklist.

**Open:**
- the guard's masks depend on `--min-baseq` (#374);
- the SNV tie rule with a confident third base (#367);
- Nextflow resources wait for the sweep (branch `feature/nextflow-resources`, local).

**Next:**
1. Class plan B through rq2 to its PR.
2. Packages 4 → 7 (table below; the clipped-bases package #212/#177 is the next plan-first item), the truth sets (#275) and the HPC gate.
3. Speed: #335 for the judge, #412 for the routing, #406's memory, and the C++ floor on the production profile before 6.7.0 ships.
4. Clean-up.
5. #331 for 6.8.0.

## Next — 6.7.0 work packages (operator, 2026-10-23: the recommended order; #304 in 6.7.0)
Tickets that answer the same question, or share code paths and measurements, land together.
- **Count-changing packages** run one at a time. Each starts from develop on its own branch,
  goes red-first, gets an adversarial review, and is measured on every arm with the checks
  above.
- **Output-only and hygiene items** run alongside in their own worktrees.
- **Then** D9–D12 on HPC, and one 6.7.0 freeze and gate.

| # | Package | Tickets | Why together | First step |
|:--|:--|:--|:--|:--|
| 1 | Co-annotated alleles (**done**: #305 closed, #246 parked into 2) | #305, #246 | the AD-claiming guard and RJ-5/RJ-25; one measurement set: the cluster-aware check, the compound-het arm, the panel and gate clusters | #305 read by read on the data in hand; then #246, whose SNV sibling still claims the read through its own uncontested call |
| 2 | Repeats: unit counts and placement (**merged**: the sequence judge, #304, opt-in in 6.7.0 and the default in 6.8.0, #331; #297 in #404; #264 in #396 and #413) | #290, #301, #297, #144, #304; #264 (docs); #246 | one question, in four code paths: which length of a repeat a read carries when the aligner placed it elsewhere. One acceptance metric: VAF at germline hets in repeats (long insertions 0.441 on HG002) | decide the architecture first. Measure #304 (a whole-read haplotype judge) on every arm as a prototype against targeted fixes, then build; #264's documentation with it |
| 3 | RNA read judgment (**#292 merged** in #415, B: no junction penalty; #293 and #178 next) | #292, #293, #178 | the spliced-read path; FORTE | #292's three prototypes, measured (decided 2026-10-17) |
| 4 | Evidence outside the aligned bases | #177, #212, #179 | read admission and clipped bases (RJ-10, RJ-13); WES and HG002 without realignment | #212 first: the fragment end bounds what a clip may hold |
| 5 | Guard statistics and artifacts | #306, #180 | telling a recurring change from an error or an artifact (fragment shares, pair orientation; #302 added the binomial orientation veto) | count the rows near the 5% line; decide #180's form |
| 6 | Homopolymer twin and decomposition | #111, #112 (#145, #146), #128 | the opt-in `--rescue-homopolymer` path | decide fix or retire (#111); retiring closes #112, #145 and #146 |
| 7 | Complex-event read-on | #274 | RJ-9 for delins | measure in the engine, then decide |
| P | Output and reporting (parallel) | #273, #227, #132, #127 | no count change | #273, then #227 (it feeds the QC summary); close #127 unless a pipeline needs it |
| H | Hygiene (parallel) | #284, #216, #219, #138, #150, #195 | no count change | #195's fallback path is one Phase 3 never reaches: fix it or retire it |
| D7 | Platforms (parallel) | #232 (#233–#237) | infrastructure | as planned |
| G | Truth arms and the gate | #266, #267, #268, #275 | the HPC gate | the D12 settings on #275; HPC when available |

**Order (operator, 2026-10-23):** 1 → 2 (decision, then build) → 3 → 4 → 5 → 6 → 7, with P,
H and D7 alongside; then G, the freeze and the gate. The packages are recorded on #196;
the A–F parents stay as they are.
