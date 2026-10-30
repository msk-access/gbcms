# CONTINUITY — where we left off

> Tactical state that must survive a closed laptop or a context summary.
> Update the **Now** and **Next** sections as work progresses. When a "Now" section is finished,
> move it to `.agents/archive/CONTINUITY_ARCHIVE.md` (older sections live there).

_Last updated: 2026-10-30_

## Now — 6.7.0 package 2: the sequence judge and its follow-ups merged (2026-10-30)

**Merged to develop** (all opt-in under `--read-judge sequence`; the default path is unchanged):
- the sequence judge itself (#363);
- the call examined on the traced path, path-align (#370), and the clean-up that tells preparation the judge (#373);
- the repeat-aware band floor (#375);
- found alleles written past the guard's masks, the indel twin of RJ-28 and the shared-indel note `INDEL_SHARED` (#377);
- path steps kept only for the states a call can examine (#378);
- one 128 MB memory budget per thread for the row caches and the work space (#379);
- harness: the HPC guard hook (#362), the PHI sanitizer-first rule (#372), labels after their plain meaning (#376).

**Final develop check** (HPC, develop fe34f443, every arm): the default path is byte-identical to develop
a2ca57c1; the judge path is byte-identical to the measured lanes (the floor's and the found-allele
change's moves, nothing else). Two HG002 reads (three calls across arms) move only with the floor and
the found-allele change together; both were traced, and both final calls are right: each read fits a
long listed insertion nearby equally on every haplotype, so it is depth only.

**HPC works** (kit, ledger and copies in `~/test/gbcms/harness/hpc/`): Linux builds equal the Mac
builds byte for byte; about 2 minutes per task; arrays at most 20 at a time (the operator's guard
setting); clinical BAMs one at a time, the work-folder copies in parallel.

**Open:**
- **The band can miss a path through a long listed insertion** (found 2026-10-30): a state's banded
  score stands above the floor while a much better path through a long listed insertion lies outside
  its band. On two small HG002 row sets, 30 reads at 3 rows differ from no-band calls. It predates
  #375 and #377. Not yet sized or ticketed.
- The guard's masks depend on `--min-baseq` (#374); the SNV tie rule with a confident third base (#367).
- `MNV_SHARED` keeps a 32-row cap (needs 33 listed MNV rows over one SNV; not realistic). Noted only.
- Nextflow resources (branch `feature/nextflow-resources`, local): waits for the sweep; per-assay
  profiles held (operator, 2026-10-30).

**Next (package 2):** size and decide the band miss; then re-check and close #290, #301, #297, #144,
#246 and #264's docs with the judge; then packages 3→7 (table below), and #331 for 6.8.0.

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
| 2 | Repeats: unit counts and placement (**building**: the sequence judge, #304; opt-in in 6.7.0, the default in 6.8.0, #331) | #290, #301, #297, #144, #304; #264 (docs); #246 | one question, in four code paths: which length of a repeat a read carries when the aligner placed it elsewhere. One acceptance metric: VAF at germline hets in repeats (long insertions 0.441 on HG002) | decide the architecture first. Measure #304 (a whole-read haplotype judge) on every arm as a prototype against targeted fixes, then build; #264's documentation with it |
| 3 | RNA read judgment | #292, #293, #178 | the spliced-read path; FORTE | #292's three prototypes, measured (decided 2026-10-17) |
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
