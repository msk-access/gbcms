# CONTINUITY — where we left off

> Tactical state that must survive a closed laptop or a context summary.
> Update the **Now** and **Next** sections as work progresses. When a "Now" section is finished,
> move it to `.agents/archive/CONTINUITY_ARCHIVE.md` (older sections live there).

_Last updated: 2026-10-31_

## Now — 6.7.0 package 2: the sequence judge and its follow-ups merged (2026-10-31)

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

**Merged 2026-10-30/31:** the band rule for long listed insertions, whole-read variant (#388, #369);
the stale fallback scores (#382, #195); the docs build notice (#381); the RHEL 8 wheel without
OpenSSL (#383, #234) and one abi3 wheel for Python 3.10–3.14 (#387, #235); the measured differences
from GetBaseCountsMultiSample 1.2.5 (#386); record updates (#380, #384). Closed with the judge's
fix: #301, #144 (status note: its count rests partly on reads RJ-23 calls depth only) and #246.

**Order (operator, 2026-10-31):** all count-changing work first, then speed, then code clean-up.
The platform PRs (#389 Apptainer, #236 arm64 image, #237 macOS wheel) and the spike-in truth set
(#275) continue alongside, as they change no count.

**Speed floor** (`.agents/memory/speed-floor.md`): production-size benchmark on HPC against the C++
binary 1.2.5 (`kit/run_prod.sbatch`, `pf71.*`). The default path meets it in all 8 cases (1.4–2.4×
faster; ACCESS 1 BAM × 20,000 rows about equal); the judge does not from 5,000 rows up (1.3–5.2×
slower). The plan for the judge is on #335, after the counting work; #150's profile is the start.

**HPC works** (kit, ledger and copies in `~/test/gbcms/harness/hpc/`): Linux builds equal the Mac
builds byte for byte; up to 300 jobs at once across our arrays (operator, 2026-10-30); clinical BAMs
one at a time, the work-folder copies in parallel; containers through Singularity
(`module load singularity`), cached in the work folder.

**In flight (counts):**
- #297, the narrow unit-run rule (fix/297-narrow-unit-run, three review fixes in 8be9f3c7, HPC n99):
  ties stay partial evidence; merges after the spike-in truth (phase 1b, #275) shows its benefit.
- #264, RJ-23's equal windows in the judge with RJ-27's symmetric read-on (fix/264-equal-windows,
  HPC y64); the per-read symmetry bar replaces the 0.5 VAF bar; the cis option stays opt-in; REF
  credited on large unscored clips is its own fix before #331.
- RNA package 3 (#292): (c) and (d) dropped; the partial-span carrier rule committed; the DNA check
  of spliced-only ALT rows running on HPC; #385 (RJ-17 exon-start deletions) and the ALT-side
  insertion gap queued.
- Recorded, not fixed: the rules path's masked-comparison false ALT (on #331); the band miss for
  listed insertions inside a long deletion's span (draft, not filed).

**Open:** the guard's masks depend on `--min-baseq` (#374); the SNV tie rule with a confident third
base (#367); Nextflow resources wait for the sweep (branch `feature/nextflow-resources`, local).

**Next (package 2):** #297 and #264 to their decisions; then packages 3→7 (table below), the truth
sets and the HPC gate; then speed (#335) and clean-up; #331 for 6.8.0.

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
