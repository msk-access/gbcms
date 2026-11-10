# 292: builder's list (companion to `292-rna-exon-edges.md`)

The PR text draws on this page. Every place that describes BAQ's junction penalty, or says it mimics
GATK SplitNCigarReads, changes with option B (a splice never lowers a base's quality; the read's indels
keep the indel penalty). Line numbers are from the measurement branch (`feature/292-rna-judgment-measure`,
b247db31) and may drift; re-run the grep before editing.

## Code

- `rust/src/shared/baq.rs:9-12` (module doc: "similar purpose to SplitNCigarReads' overhang clipping"),
  `:137-150` (the splice arm of the penalty), `:142` (the comment), `:166-171` (the overhang constants
  and `overhang_bases`, measurement only: remove with the switch), `short_anchor_bases` (same).
- `rust/src/counting/junction_rule.rs` (the whole switch, `GBCMS_RNA_JUNCTION`: remove) and its lane in
  `rust/src/counting/carrier.rs:718-721`; `rust/src/counting/engine.rs:2674` (the junction view for the
  overhang rule) and the `junctions` field threaded through `effective_quals` (`engine.rs:3089`).
- `src/gbcms/models/core.py:523` (the comment on overhang clipping).
- `src/gbcms/cli.py:619-621` (the `--apply-baq/--no-baq` help for `gbcms rna`: "near indels and splice
  junctions" becomes "near indels").

## Tests

- Expectations that change (part of the decision): `baq.rs::the_decided_penalty_lowers_five_bases_each_side_of_a_junction`;
  the junction cases of `tests/test_baq.py`; the "masked junction" case of
  `tests/test_rna_spliced_read_on_contract.py`.
- The 30 measurement tests in `tests/test_rna_junction_rule_measure.py` become contract tests of B
  (red-first on today's engine), then the switch goes.
- Keep green: `test_boundary_baq_contract.py`, `test_span_exon_edge_baq_contract.py` (the exon-edge
  exception stays), `test_baq_own_evidence_contract.py`, `test_rna_partial_span_contract.py` (68),
  `test_rna_splice_contract.py`.

## Docs and reference pages

`grep -rn SplitNCigarReads docs src nextflow rust/src README.md CHANGELOG.md .agents` finds 14 files.
Four stay as they are (history or a survey list): `.agents/learnings/REJECTED.md:373`,
`.agents/learnings/LEARNINGS.md:261`, `.agents/rules/decision-protocol.md:70`, `CHANGELOG.md:3334`.
Change the rest, plus the lines that describe the junction penalty without naming GATK:

- `docs/cli/rna.md:97` (BAQ on by default for RNA; "indels and splice junctions").
- `docs/cli/dna.md:100-120` (the BAQ section and its table; the RNA note at :120).
- `docs/nextflow/parameters.md:70` (`--rna_apply_baq`: "the penalty also applies near splice junctions").
- `nextflow/README.md:106` (`--rna_apply_baq`; no flag change, the description only).
- `docs/reference/read-filters.md:154-156` (the `--min-baseq` note), `:332-347` (the BAQ section: the
  radius table names splice junctions), `:354` (the mode table), `:371` (the RNA paragraph).
- `docs/reference/rna-splice-handling.md:124-135` (the SplitNCigarReads section), `:177`, `:210` ("this
  mimics SplitNCigarReads' overhang clipping"), `:291`, `:306`; the gbcms-approach section's "heuristic
  BAQ handles overhang misalignment" sentence.
- `docs/reference/read-judgment.md`: a new register row for B (amending RJ-27's "with BAQ masking the
  bases at its junction" clause and the RJ-17 amendment note of 2026-10-30), and the survey line at
  `:376`.
- `docs/reference/architecture.md:509` (the page summary), `docs/reference/bam-evidence-caveats.md`
  (its BAQ-and-junction caveat).
- `.agents/skills/read-filters-qc/SKILL.md:36` ("the RefSkip penalty substitutes for SplitNCigarReads").
- `CHANGELOG.md`: a new entry; the exon-edge exception (R1) is unchanged and stays documented.

## Effects map to state in the PR (changes / unchanged, because)

Both read loops (main `engine.rs:1298`, per-transcript `:1319`): change alike, one `effective_quals`.
The rules path and the sequence judge (`engine.rs:2810`): change alike. ASJD allele tallies and mFSD:
change with the calls, no formula change. MNP rescue: reads the same qualities. Writers, diagnostics,
flags, Nextflow flags: unchanged. DNA: unchanged by construction (no junction); one DNA identity batch
(counts byte-equal) before merge.
