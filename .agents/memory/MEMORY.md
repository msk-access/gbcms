# gbcms — durable memory index

One line per memory. Full content lives in the linked file. Keep this index tight.

## How we work (feedback)
- [Harness stays small](harness-stay-small.md) — curate skills/rules ruthlessly; every line is paid for each turn.
- [No ticket labels in code](no-ticket-labels-in-code.md) — comments/logs explain what/why/how, never CR-/HI-/ME-/P4c labels.
- [Describe labels in brackets](describe-labels-in-brackets.md) — when explaining, say what a thing is first and put phase/RJ/build/option labels in brackets; never lead with a bare label.
- [Identity band = annotation tolerance](identity-band-annotation-tolerance.md) — the ≥90% insert band exists for imperfect ALT representations; BQ masking can't replace it — test both policies. No new columns/flags: new signal goes to logs/diagnostics/validation tooling.
- [Commit before review workflows](commit-before-review-workflows.md) — same-checkout review agents may `git stash`; concurrent sessions write .agents/memory — commit first, stage explicit paths only.
- [Genotyper, not caller](genotyper-not-caller.md) — gbcms counts reads whose own bases carry the given allele; alignments can be wrong, so judge bases not placement; validate per read across data types.
- [Count the given allele](count-the-given-allele.md) — input is taken as correct; never deconvolute a wrong input into counts; explain via status-reason/diagnostic columns; only opt-in rescues differ.
- [Use GitHub sub-issues](github-sub-issues.md) — group related work as parent + sub-issues (cycle tracker, umbrella items), not combined issues or checklists.
- [Holistic effects map](holistic-effects-map.md) — before any fix/decision, map every place it lands (both counting paths, per-transcript/ASJD/mFSD/observations, rescue/clusters/merge, writers/flags, docs/tests/Nextflow).
- [Fix the class, not the row](fix-the-class-not-the-row.md) — state the mechanism first (no row in it), enumerate the class as tests, derive cutoffs from a mechanism, justify narrowing, judge on data the bug wasn't found in.
- [Read the rule register to the end](read-rule-register-to-the-end.md) — before stating a decided rule, read its whole read-judgment.md row (amendments are appended) and the contract test docstring; reference pages can lag.
- [Test changes need operator notice](test-changes-need-operator-notice.md) — never edit a failing test's expectation without first telling the operator what it asserts, why it's wrong, and the evidence; wait.
- [Survey several tools](survey-several-tools.md) — community practice means GATK, samtools/bcftools, fgbio, VarDict, Strelka2, freebayes, bam-readcount, LoFreq, GetBaseCounts, the literature, domain tools (RNA), and tools of other applications that judge reads alike (SV/graph genotypers, realigners, split-read and RNA indel tools); say where handling isn't documented.
- [Test across data types](test-across-data-types.md) — accept a count-changing rule on every data type (tissue tumour+normal, cfDNA, WES/WGS, RNA) and input type, with more than one truth source; GIAB alone is not acceptance; DNA and RNA may need different rules.
- [Parameters are user-settable](parameters-are-user-settable.md) — --min-baseq etc. are user flags: truth/oracles take the run's value, never the default; sweep quality-sensitive rules across values.
- [cBioPortal stays local](cbioportal-stays-local.md) — tool output goes to the model provider: never read cBioPortal files with a tool; only the fixed schema-and-counts script prints; finer work writes to protected local files (hook-enforced). Patient slices are fine for testing.
- [Better than the standard](better-than-standard.md) — the survey is a floor, not a ceiling: where no tool sets a standard or a measured rule beats the field's, adopt it and say why.
- [Best, not least-bad](best-not-least-bad.md) — in analysis and implementation alike: name what the data can't see, stress it with synthetic shapes and diagnostics, then make every shape a red-first test and land the rule in every path.
- [Model per task](model-per-task.md) — set each subagent's model by the task: Fable for research surveys, Fable/Opus for deep reasoning (reviews, rule design, decisions), Sonnet for routine work, Haiku for mechanical scans; suggest a session switch when the phase changes.
- [External review: measure first](external-review-measure-first.md) — map an outside review to behaviour, tickets, REJECTED.md and old probes, and measure its value, before filing tickets.
- [mFSD is graded evidence](mfsd-graded-evidence.md) — fragment size gives increased confidence toward CH or toward somatic (or none), never a hard origin call; don't over-interpret.
- [Decisions with worked examples](decisions-worked-examples.md) — count-rule decisions come with a real row drawn out, the community standard, the RNA effect, a truth check, and every option offered measured (not estimated).

## References
- [Claudelicious harness](claudelicious-reference.md) — the upstream pattern this project's harness follows.
- [Fragmentomics (Tsui et al., MSK)](fragmentomics-reference.md) — the prior for mFSD and the CH-vs-tumor question; cite for S1 #153 / S2 #154.
- [Downstream org repos](downstream-org-containers-modules.md) — mskcc-omics-workflows containers image (6.3.1) and gbcmsrs modules; not a release step; the operator starts the switch when satisfied it's production-ready (2026-10-07).
- [Production aligner settings](production-aligner-settings.md) — nucleo (ACCESS; IMPACT's stand-in) BWA+ABRA2, forte STAR 2.7.10a; FORTE's fill-out still calls the C++ gbcms 1.2.5; full flags on #275.
- [Lancet2](lancet2-reference.md) — NYGC somatic caller: cDBG assembly per ~1 kb window, linear minimap2 genotyping per haplotype, small variants only; borrow its two-tier gap model, cDBG allele discovery, complexity guard; not its margin-free best match.

## Project facts (from the 2026-06-26 code review)
- [Universal haplotype judge](universal-haplotype-judge.md) — DECIDED 2026-10-23: package 2 adopts a per-haplotype read judge (edit scorer, context x site-state closure, one-base margin) for small-variant rows, flag-first in 6.7.0; `rules` kept one release after the flip (#328); SV/fusion work is 6.8.0 (#322–#327).
- [Bin fetch-end must cover the anchor variant](bin-anchor-coverage.md) — CR-1; pinned by the bin property test and binning-invariance tests.
- [WFA fast-path shares the base-quality gate](wfa-bq-gate-contract.md) — CR-2 cross-backend quality contract.
- [Engine should be output-aware](engine-output-aware.md) — plumb intent across FFI vs compute-then-discard.
- [Tolerant large-deletion match — superseded](tolerant-deletion-deliberate.md) — issue #91: real large dels are exact-length; wrong-length pure indels → partial_alt; delins stay Phase-3. The 50bp gate is an artifact-SIZE prior (artifacts are small; a ≥50bp op is real), not an event-rarity claim.
- [Nextflow defaults diverge from CLI](nextflow-cli-default-divergence.md) — keep nextflow.config in sync with CLI defaults.
- [RNA MAPQ default is 1](rna-mapq-default.md) — keeps STAR's 2–4-locus reads (pseudogene-tied junction reads); unique-only cost PIK3CA E545K ALT; decided 2026-10-03.
- [MAPQ-0 loci (PMS2)](mapq0-loci-pms2.md) — pseudogene genes run at --min-mapq 0; keep MAPQ-0 alignments countable; validate read-admission changes at --min-mapq 0 too.
- [ABRA2 junction placement](abra2-junction-placement.md) — IMPACT/ACCESS are ABRA2-realigned (YO = original BWA alignment; junction inserts carry assembly evidence); WES (BWA only) is the placement control; 50-500 bp events are called off ABRA2's CIGARs with no SV-caller second opinion.
- [ABRA2 after collapsing](abra2-collapse-artifacts.md) — ACCESS runs ABRA2 twice (standard, then post-collapse); the second pass can introduce events: check the standard BAM; CIGAR-trusting rules count them, bases-judging rules don't.

## Tooling / build
- [pyproject is the dep source of truth](deps-pyproject-source-of-truth.md) — floors measured (CI floors leg), dev tools in PEP 735 groups, only the image installs a lock (hash-pinned + pip check); declare every directly-imported package.
- [Lint-tool version skew (local vs CI)](black-version-skew-venv-vs-ci.md) — venv black 25.9 lags CI 26.5 (rustc caught up to 1.96 on 2026-09-23); run CI's versions before trusting a clean/drift result.
- [Worktree tests need an isolated venv](worktree-tests-need-isolated-venv.md) — the shared .venv imports the MAIN checkout's gbcms; in a worktree build a scratch venv + maturin develop. Every review-agent prompt says `VIRTUAL_ENV=<scratch venv> maturin develop`, or maturin installs into the main .venv. Scratchpad venvs rot after ~3 days (macOS purges /private/tmp): rebuild, don't repair.
- [maturin develop: repo root only](maturin-develop-repo-root-only.md) — `-m rust/Cargo.toml` bypasses [tool.maturin] and leaves a stale src/gbcms/_rs.so shadowing every rebuild.
- [Harness cleanup after merge](harness-cleanup-after-merge.md) — delete the ticket's `src_*` builds (rebuildable from BUILT_FROM; patch prototypes first), gzip traces, `cargo clean` when rust/target bloats.
- [Merge with gh when asked](merge-with-gh-when-asked.md) — repo disallows GitHub auto-merge; when the operator asks, check CI once and `gh pr merge --merge` if green; never poll.
- [Long runs visible](long-runs-visible.md) — launch long runs as the tracked background task itself; nohup only to outlive the session, then a waiter on the runner's PID (never a log marker); network-mount runs one at a time, grouped by file.
- [HPC partitions](hpc-partitions.md) — never run on the head node (srun/salloc shell, check hostname); cmobic_short ≤3 h + priority QOS, cmobic_cpu (7-day) for longer.
- [No `timeout` on this Mac](macos-no-timeout-command.md) — `timeout N cmd && ok || fail` always says fail; check SFTP mounts with a plain read.
- [Pipes mask gate exit codes](pipe-masks-gate-exit.md) — `pytest | tail && git commit` commits on failure; capture `$?` or use pipefail.

## Validation / testing
- [QC-fail flag absent from MSK data](qcfail-flag-absent-msk-data.md) — no pipeline stage sets 0x200; filter verified correct but inert in practice (contract test pins it).
- [--trace output wraps](trace-output-wraps.md) — trace through `plan670/r2/plain_trace.py` (no Rich: no COLUMNS, no wrap, no padding; identical calls, 3x faster); plain gbcms wraps at COLUMNS and pads to it; a read call is not a count (take totals from the MAF); match keys by haplotype.
- [Census mirrors read inputs](census-mirrors-read-inputs.md) — the census shares read-input rules (fragment end, filters), so validate those against an oracle it doesn't share (the bases themselves, the mate's alignment, ITD rows first).
- [VAF uses informative reads](vaf-informative-reads.md) — `vaf` = AD/(RD+AD); per-sample total_count (DP) includes unjudgeable reads; never advise AD/DP (repeat indels: reads ending inside the tract).
- [pysam validation oracle](pysam-validation-oracle.md) — use fetch()+get_reference_positions (not pileup) to cross-check gbcms counts; RD/AD match exact, DP includes neither.
- [Cluster-aware oracle](cluster-aware-oracle.md) — score changes at co-annotated rows against every allele of the cluster (`m294/oracle_sibling.py`); REF-vs-ALT checks call removed sibling molecules "lost"; a sibling MNV that contains the row's change is the row's ALT (decided 2026-10-25).
- [ACCESS XS1 vs XS2 qualities](access-xs1-xs2-consensus.md) — XS1 = Marianas (simplex Q80–90), XS2 = fgbio (simplex capped Q45, duplex 68/79/90); split ACCESS results by generation; no single high threshold means "duplex-agreed" on both.
- [Harness references](harness-references.md) — panel/ACCESS/WES slices are b37 (Homo_sapiens_assembly19 on OneDrive), HG002 and FORTE GRCh38 (~/Downloads/genome.fa); brief agents with the right one.
- [Other tools on the slices](other-tools-harness.md) — the 10-tool comparison harness (plan670/tools): one container at a time, IMPACT/ACCESS slices have no proper-pair flag (bcftools -A), Delly/Mutect2/Paragraph quirks, decomposed delins.
- [vcf2maf oracle](vcf2maf-oracle.md) — VCF↔MAF representation is defined as vcf2maf/maf2vcf output; how to run them locally and their quirks (dies on N, empty ALT for REF==ALT).
- [GIAB truth arm](giab-truth-arm.md) — HG002 v5.0q hets on the DeepVariant 35x WGS BAM; score with truth v2d (`--truth v2d`/`v2dq20`: hom rows vs REF, benchmark by variant size, decidable by bases); fetch slices in parallel; run untraced.

## User (private, local-only — not committed)
- `user-*.md` memories (e.g. who the operator is, personal defaults) live in this
  recall dir but are gitignored, never published. See `.gitignore`.
