# gbcms - durable memory index

Read each turn. Open a file when its trigger fits. Standing rules that live elsewhere are listed
here with their home, so none is lost. Keep this index under 5 KB.

## Standing operator rules (home in brackets)
- PHI: cBioPortal files never enter tool output; patient slices only for testing, never downstream (rules/security.md)
- A failing test's expectation is never edited without telling the operator first (AGENTS.md)
- A count-rule decision: read the register to its end, map the effects, fix the class, survey several tools, test every data type, show a worked example, measure first (rules/decision-protocol.md)
- Merge a PR with `gh` only when asked; use sub-issues; commit before review workflows; stage explicit paths (rules/git-conventions.md)
- Read every gate exit code; CI's black; never pipe a gate (skills/qa-check). Build, venv, stale `.so` (skills/rust-python-ffi)
- Long runs as tracked tasks; never on the HPC head node; oracles, traces, FASTA per arm, build cleanup after a merge (docs/harness/measurement-harness.md)
- No ticket labels in code or logs (rules/code-quality.md)
- [Working agreements](working-agreements.md) - keep the harness small; explain in STE; say what each thing is, then its label in brackets (tickets, rules, branches, builds, jobs, runs), in tables and relayed agent reports too; model per task.
- [Genotyper, not caller](genotyper-not-caller.md) - count the given allele; judge bases, not placement; rescues are opt-in; the insert band survives only in the CIGAR fallback; no new columns or flags.

## Decided behaviour and facts
- [Decided counting rules](decided-counting-rules.md) - wrong-length indels are `partial_alt`; complex events and 50+ bp deletions are judged by bases; RNA `--min-mapq` 1; pseudogene loci at 0; `vaf` = AD/(RD+AD).
- [mFSD is graded evidence](mfsd-graded-evidence.md) - leans somatic or says nothing; never an origin call; the buffy coat decides CH.
- [Universal haplotype judge](universal-haplotype-judge.md) - DECIDED 2026-10-23: one per-haplotype read judge; opt-in in 6.7.0, default in 6.8.0 (#331); SV and fusions 6.9.0 (#322-#327).
- [Speed floor](speed-floor.md) - faster than the original C++ GetBaseCountsMultiSample as a minimum; the judge meets it before it is the default.
- [MSK assay facts](msk-assay-facts.md) - ABRA2 realigns IMPACT and ACCESS (ACCESS twice); XS1 Marianas against XS2 fgbio qualities; production aligner flags; no QC-fail flags.
- [References](references.md) - Lancet2 (borrow its gap model, not its best match); downstream org repos (never propose the switch; the operator starts it).

## Private
- `user-*.md` files are local-only and gitignored (who the operator is, personal defaults).
