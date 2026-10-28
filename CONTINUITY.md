# CONTINUITY — where we left off

> Tactical state that must survive a closed laptop or a context summary.
> Update the **Now** and **Next** sections as work progresses. When a "Now" section is finished,
> move it to `.agents/archive/CONTINUITY_ARCHIVE.md` (older sections live there).

_Last updated: 2026-10-26_

## Now — 6.7.0 package 2: the sequence judge integrated, fix round running (2026-10-26)

**Release plan (operator, 2026-10-26).**
- **6.7.0** ships the judge opt-in (`--read-judge sequence`, #304).
- **6.8.0** does one thing: make it the default, meticulously and systematically. The parent issue is #331, with sub-issues:
  - #332 tests and contracts, and every output path;
  - #333 RNA;
  - #334 acceptance and the `--min-baseq` sweep;
  - #335 speed;
  - #336 the HPC gate;
  - #337 the switch-over.
- **6.9.0 (new)** takes SV/fusion genotyping (#322–#327) and retiring the rules judge (#328).
- The 6.6.0 milestone is closed.

**Judge branch** `feature/c46-haplotype-judge`: pushed as a backup, no PR. The integration build is frozen as `plan670/src_c53d`. It has:
- the graph scorer;
- the fix round for long deletions and MNVs;
- the MNV/SNV rule with `MNV_SHARED`;
- the fits-nothing rule (xs) with `UNEXPLAINED_READS` (RJ-29);
- clipped read ends at SNV/MNV rows admitted (RJ-13);
- the name `--read-judge sequence`, with `haplotype` as a deprecated alias;
- a Nextflow `read_judge` parameter (DNA).

Leaving out context items that respell a candidate sequence (`GBCMS_JUDGE_BACKGROUND`) was measured and not adopted: the background stays as is. The rules path is byte-identical to develop on 616 runs.

**Truths, corrected 2026-10-26:**
- **HG002 truth v2c:** BWA alignments rescored end to end from the CIGAR; score with `--truth v2cq20`.
- **Patient oracle v2** (`xdata/xoracle_v2.py`): per-haplotype indices plus an exact DP, and a label needs both.
- **Bars:** count reads decidable at the run's `--min-baseq`; long deletions by breakpoint reads; a class fails only past 0.5 points AND 2 reads; wrong-call bars past 2 extra reads.
- **Scratch scripts** are mirrored in `plan670/sp_tools/`.

**Re-score of `src_c53d`** (v2 and v2c; details on #304):
- **No agreement bar fails on any patient arm:** panel +2.68 pt, ACCESS XS1 +5.61, XS2 +7.88, normals +2.50, WES +4.46, gate +0.66.
- **Gate groups below clinical sign-out:** in all nine, the oracle sides with the judge.
- **Open:**
  - wrong-call bars: gate MNVs +743, delins, gate unique 1 bp insertions other→ALT +114, unique deletions of 50+ bp other→REF;
  - chrX: 56 true-ALT reads lost; SNVs near an indel: 21;
  - true-REF reads called partial at tract clusters and 36 bp duplications;
  - HG002 SNVs in tandem repeats: −40 true-REF reads;
  - deletions of 1 kb+: scaled calls 162 more true-REF reads REF on HG002, while xs keeps 307 (panel) and 78 (gate) more correct calls.

**Fix round 1 (`src_c54a`): none of its four items was a judge bug.**
- **Insertion readability:** the judge already follows RJ-20 as amended on 2026-10-15 (#278): an unreadable insert of the ALT's length inside the window is the ALT, on both paths. The "port" rested on stale reference text, and the operator is re-confirming the rule.
- **HG002 SNVs in tandem repeats:** the −40 reads are ties that only phasing separates. Develop's placement tie-break would make 443 wrong calls against the judge's 2.
- **chr18:12145098:** a truth error, not a judge error. The row is outside the benchmark, and v2c's homozygous shortcut never compares a read with REF.
- **Deletions of 1 kb+:** scaled's 162 "gains" are reads the engine never counts (they start inside the deleted bases). xs stays.

**In flight:**
- **Fix round 2:** the six diagnosed judge bugs, and a measured switch for "a confident third base inside X's span is another allele".
- **Oracle v3:** the patient check, with six fixes.

**Decided 2026-10-24 to 2026-10-26:**
- **Clip admission:** admit. `clean:0` costs 8% of true ALT on patient data to avoid about 1% wrong reads.
- **Fits-nothing:** xs; a refinement for deletions of 1 kb+ is pending.
- **MNV/SNV:** O5 plus `MNV_SHARED`; oracles score containment as the row's ALT.
- **Contexts:** listed rows, the guard background (as is) and found alleles all stay.
- **`--min-baseq`:**
  - it masks bases, never reads (#329); MSK runs at 20;
  - the sweep runs at 0, 10, 20 and 30, chosen per arm by its quality bins;
  - ACCESS is reported apart by collapser: XS1 Marianas, XS2 fgbio;
  - a quality-scaled cost (cap 30 or 40) is adopted only if it wins at every threshold;
  - the below-threshold "floor" cost is dropped.
- **PHI:**
  - tool output goes to the model provider;
  - patient slices are fine for testing, never downstream;
  - cBioPortal files never enter tool output (#330). The only console channel is `~/test/gbcms/harness/cbio/cbio_summary.py`; protected paths live in `~/.claude/protected-paths.txt`.
- **PRs:** the repo disallows auto-merge, so merge with `gh` when the operator asks.
- **Length-only insertion stays the ALT** (RJ-20 as amended on 2026-10-15, #278; re-confirmed 2026-10-26): an insert of the ALT's length inside the window whose bases are all masked is the ALT where no other allele of that length is at the locus.
- **Ties stay depth only:** a read that fits REF and another candidate equally counts in DP, not in RD or AD, as GATK's uninformative read does.
- **HG002 bars, truth v2d (done):** compare homozygous rows with REF; the benchmark BED by variant size; long deletions on anchor or junction-holding reads; count only reads decidable by their own bases. Ties that only phasing can decide go to #340.
- **Both "keep" decisions measured:**
  - length-only inserts: 0 truth-REF reads, 0 contradicting mates;
  - ties as depth: VAF error 0.079, against 0.096 for develop.

**Next (package 2):**
1. Decisions from the two reports.
2. The `--min-baseq` sweep and the scaled-cost option.
3. Speed: at most 2× develop. Candidates: Lancet2's two-tier gap model and k-mer allele discovery.
4. The judge PR to develop, opt-in. Then re-check and close #290, #301, #297, #144, #246 and #264's docs.
5. Then packages 3→7, with P, H and D7 alongside; then G, the freeze and the gate.

**Main checkout:** on develop once this PR merges. It was on 336c07cd, with these memory notes uncommitted.

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
