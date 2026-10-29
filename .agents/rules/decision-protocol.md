---
description: How a count-affecting decision is made, measured and presented to the operator (worked example, several tools, every data type, fix the class, effects map)
alwaysApply: false
---

# Decision protocol (count-affecting changes)

Read this before you propose, decide or implement a change to how reads are counted. The spec-first
steps and "fix the class, not the row" are in `docs/reference/read-judgment.md`, "How a change to
read judgment is made"; this rule holds the rest. The operator decides; a decision gets a register
entry (`docs/reference/read-judgment.md`) before code. A failing test is a decision point, never
edited away (AGENTS.md).

## 1. Read the rule before you state it

Before you say what a decided rule says, read its whole row in `docs/reference/read-judgment.md`
(the register appends amendments at the end of a row), then the docstring of its contract test.
`allele-classification.md` and `glossary.md` can lag an amendment. When a rule is amended, put the
current rule at the start of its row and update the reference pages in the same change.

Why (2026-10-26): I told the operator that the read judge broke the insertion readability rule
(RJ-20). The row ends with an amendment of 2026-10-15 (RJ-27, #278): an unreadable insert of the
ALT's length inside the window is the ALT. The judge and the rules path both follow it, and
`tests/test_unreadable_insert_contract.py` pins it. I had read the stale reference page and the
start of the row, queued a wrong "port", and wrote it into #304 and #329.

## 2. Map the effects

Before you implement a decision or fix, state in the PR and the plan ticket every place its effect
lands, each as "changes" or "unchanged, because ...":
- both read loops (the main loop and per-transcript counting);
- per-transcript counts, ASJD, mFSD and the observations export;
- MNP rescue, tract-cluster claiming and merge;
- the MAF and VCF writers, the diagnostics and the flags;
- docs, tests and Nextflow.

Trace the fields the change touches (`is_ref`/`is_alt`, `has_nearby_evidence`, `qual`) to every
consumer (grep them), and verify the "unchanged" claims on real data as well as the "changes". Check
a derived output (`vaf`, a flag) at every place it is computed, inputs included, before you say how
a change moves it: `vaf` is alt / (REF + ALT) in the writers, and merge's alt / total uses a merged
total that is REF + ALT, not the per-sample depth (a claim read from one formula's surface was wrong).
A rule ported from one judge to another is re-derived on the target's own inputs: the exact-carrier
guard's germline-run masking, ported as-is into the junction guard, masked another breakpoint's run
as background (47 WES reads turned REF).

## 3. Fix the class, not the row

The full text is in `read-judgment.md`. Additions:
- Treat a revert, or a narrowing made to pass a bar, as a signal that the mechanism was wrong. Ask
  that before you narrow further. (The judge's second fix round: 3 of 6 fixes reverted, a tie rule
  narrowed to substitution rows and reverted anyway, a 20 bp cut placed where the measurements
  split.)
- Sweep any constant users can change (section 7).
- Report each fix as mechanism, general rule, this instance, what else it touches (section 2).

## 4. Survey several tools, with sources

A decision that needs community practice surveys several tools, not one, and says plainly where a
tool's handling could not be found (links; "not documented"). The set:
- Callers and genotypers: GATK (HaplotypeCaller, Mutect2), samtools/bcftools mpileup (htslib),
  fgbio, VarDict, Strelka2, freebayes, bam-readcount, LoFreq, the original C++ GetBaseCounts, and
  Octopus.
- Tools of other applications that judge reads alike (operator, 2026-10-12): SV genotypers
  (SVTyper, Paragraph, GraphTyper2, Delly genotyping, Manta, smoove, SVJedi, GATK-SV),
  graph and realignment genotypers (vg, BayesTyper, ABRA2, GATK IndelRealigner), split-read callers
  (Pindel, ScanIndel), RNA indel and junction tools (transIndel, RNAIndel, STAR chimeric output,
  Arriba), RNA fusion callers (STAR-Fusion, FusionCatcher, JAFFA, CICERO), and read simulators and
  spike-in tools for truth (BAMSurgeon, wgsim, ART).
- The published literature (method and benchmarking papers, with DOI or PMID).
- For RNA, the domain's tools: GATK RNA practice and SplitNCigarReads, STAR, ASE counters
  (ASEReadCounter, phASER, WASP), RNA variant callers (RNA-MuTect, SNPiR) and editing callers
  (REDItools, JACUSA). Survey DNA and RNA practice separately where they differ.
- Most SV callers work above a size floor and small-variant callers cap indel length: survey each
  tool's working size range and its size-stratified accuracy, so a 50 to 300 bp gap is visible, not
  assumed covered.

The evidence column of every decision table lists what each surveyed tool does.

## 5. The standard is a floor, not a ceiling

Where the tools disagree, document nothing, or none handles the case, or measurement shows a rule
better than the field's, propose the better rule with the measurement (operator, 2026-10-04). Show
the field's practice beside the rule the measurements support, say plainly that it departs from
practice and why, and still measure the cost (real carriers lost, truth rows moved) and survey
first. (Example: no surveyed tool confirms a junction-adjacent base against the next exon; judging
spliced reads on the spliced reference cut spurious ALT from 19 to 4 on the probes.)

## 6. The best rule, not the least bad

When every measured option costs something, do not keep the cheapest and pin its defect (operator,
2026-10-13, on #277). First name what the measurements cannot see (the risk the rule guards, a truth
source that would show it), measure that, and look for a rule better on every axis. Only when none
exists, recommend with the residual named and the measurement that would settle it. Every decision
table gets a row for "what this measurement cannot see" and how to measure it; prefer building that
measurement (a spike-in, a hotspot audit, a synthetic competitor) over choosing among
under-measured options.

The same holds in the implementation (operator, 2026-10-14): stress each candidate with synthetic
shapes real data cannot show (a germline SNP in the grown repeat, a low-VAF somatic on a germline
het, a phantom with the allele absent, natural quality), score VAF and phantom share, keep oracles
free of alignment trust, and size a suspected defect with a count-neutral diagnostic on every arm.
Then make every shape a red-first test, land the rule in every path that judges the same question,
give each new component its own unit test, carry each exclusion and threshold of a rule extended to
the other allele with a mirrored test, and rerun the measurements on the final build.

## 7. Every data type, more than one truth source

A count-changing rule is accepted on every data type and input type it reaches, against more than
one truth source (operator, 2026-10-12):
- Data types: DNA tissue panel tumour with its matched normal, cfDNA duplex and simplex (through
  `gbcms merge`), WES and WGS without realignment, RNA.
- Input types: the aligner and realigner, paired or unpaired, clipped, split and hard-clipped
  records.
- GIAB HG002 is one genome, one aligner (BWA-MEM), germline hets only, and its 50+ bp deletions are
  mostly tandem repeats: it cannot show false ALT (no hom-REF rows) or how ABRA2, STAR or cfDNA
  reads behave, so it is never acceptance alone.
- Each decision table names its data types and truth sources and says which are missing. Run the
  normals (hom-REF) every time. Where no truth exists, build it (spike-ins). DNA and RNA may need
  different rules where their reads differ (splices, clips holding the next exon); say so in the
  rule.

**Parameters are user-settable** (operator, 2026-10-26). `--min-baseq`, `--min-mapq`, the gap
penalties and every margin exposed as a flag are user settings; their defaults are only defaults.
Truth and oracle scripts take the run's own value, never a hard-coded default (the patient oracle
and truth v2's `v2q20` both hard-coded 20 once). Write bars as "decidable at the run's
`--min-baseq`". Sweep a count-changing rule whose behaviour depends on a parameter over plausible
values (for `--min-baseq`: 10, 20, 30 and 0) in its acceptance. MSK runs at 20, so the default
stays; the sweep checks that the rule holds at the values other runs need, and is not a search for
a new default. ACCESS qualities depend on the collapser, so report ACCESS split by generation
(`msk-assay-facts` memory).

## 8. Present the decision with a worked example

The operator decides from a worked example, not from aggregate tables alone. In the first
presentation bring:
- one real row drawn out as a picture: the read windows or read classes each option makes a read
  hold, with that row's counts and VAF under each option;
- the community standard (section 4);
- the RNA effect, even when RNA has few such rows (say how many);
- a truth check: does the AD or VAF change move toward a known truth (GIAB, or germline hets at
  0.5), not only toward the sign-out;
- every option measured, not estimated: prototype each one and run it on the same sets and truth
  arm as the recommended one, so the comparison table holds only runs. The synthetic sets must
  mirror a real library (both pair orientations, errors, low depth and depth), or an option keyed on
  that feature is not measured at all. (On #294 only the current guard option was measured; measuring
  the others reversed one recommendation and exposed a simulator that used one pair orientation.)
- the decided-rule contract tests the change must keep, beside the arms. A pre-registered bar names
  both: an arm-only bar would have removed the judge's found alleles away from the site, because no
  arm moved, while six contract tests of a decided rule failed without them.

Build the truth arm before the decision round and validate the measurement model on truth too (the
first HG002 pass exposed an oracle placement bug the panel data had hidden). The measurement
harness is `docs/harness/measurement-harness.md`.

## 9. An outside review: measure first

When the operator brings an outside review (another model's framework, a paper's proposal), compare
it against all four: current behaviour (verified in code and docs, not from memory), the open and
closed tickets and the cycle's parents (closed ones' comments hold past measurements),
`.agents/learnings/REJECTED.md`, and the earlier review probes (rerun on the release build and
now). Then measure each proposal's value on our data before you offer a ticket: only what
measurement shows relevant becomes one. Build a map (proposal, status, ticket or REJ id), rerun
prior probes against a release-tag build, and measure the incremental value (what it adds beyond
existing flags and tickets) on the panel, truth and BWA-only arms.

Why (2026-10-09): an external "evidence accounting" review. Measuring changed the picture:
explainability flagged 43 rows but 39 already carried a flag; R1/R2 disagreement was 23 of about
195k fragments; no panel insertion was 100+ bp; hard clips mattered on BWA-only WGS (3.8% of
breakpoint records) and not on ABRA2 panels (0%); and rerunning the old probes found C33 (#214)
already fixed by C39.
