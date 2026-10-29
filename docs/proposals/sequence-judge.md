# Proposal: the sequence judge (package 2, #304)

**Status:** v1 — architecture decided (operator, 2026-10-23; record on #304); this document
fixes how it enters the engine. **Type:** read judgment, behind `--read-judge sequence`
(default `rules`; `haplotype`, the judge's name until 2026-10-26, is taken as `sequence`
with a warning). **Compat:** default output byte-identical; the flag changes counts.
**Version:** 6.7.0 (default flip only if every arm, the cost target and the HPC gate pass).
All file:line references are to develop 336c07cd. Internal names (`judge::haplotypes`,
`ReadJudge::Haplotype`) keep the old word until the module is restructured (R6, #314).

## 0. The principle (operator, 2026-10-25)
The judge resolves the read's sequence at the row: which local sequence its own bases spell
there, REF, the row's ALT or another allele at its site. Contexts (the other listed rows, the
guard's recurring background, the alleles the reads spell away from the site) only explain
the read's *other* differences, so they never pass for the row's allele (#246: an insertion
beside an SNV whose first inserted base reads as the SNV's ALT). The judge does no linkage
and no phasing: a variant linked to the row elsewhere on the molecule never credits a read.
The sequences a read is weighed against are the row's *candidate sequences*.

## 1. The rule
For a row X and each admitted read:
- **Candidate sequences:** every *context* of the other listed rows within read reach of X
  (and the guard's recurring background changes, which explain reads but are never
  credited), crossed with every *state* of X's site: REF, X, each other listed allele whose
  REF span overlaps X's (an insertion's span is its anchor), and each recurring discovered
  allele there (RJ-25; REJ-20261007-002).
- **Scorer:** the read, every readable base, against each candidate sequence with its ends
  free: affine edit cost, confident mismatch 4, gap 6 + k, a base below `--min-baseq` or N
  costing 0 (the one quality gate).
- **Decision:** ALT where the best candidate sequence carrying X beats the best without it
  by one confident base (4); REF where the best with REF at the site beats every other site state by
  as much; partial evidence where another site allele wins so; otherwise depth only. Ties are
  never credited. Because every context appears in every state, only the read's bases at X
  decide; a linked variant elsewhere never credits it.
- **Phase C refinements (2026-10-25):**
  - A state holding X's change at X's site (for a substitution X, an allele of its length
    over X's bases that reads X's bases there, as an MNV listed beside its component SNV)
    is X's ALT (RJ-28); `MNV_SHARED` names the reads the two rows share.
  - A same-length state's changes beside a substitution X are also contexts of every state,
    so a read carrying X with such a change is X's.
  - A substitution inside a long (50+ bp) deletion's bases, past a 4-base junction window at
    each end, is a context on REF's side, never another allele at X's site (RJ-24); inside
    a shorter deletion it stays another allele (RJ-5).
  - A tie of REF with another allele, X beaten, is partial evidence only where the read
    reaches past where the two candidate sequences differ; a read short of it shows REF.
  - An MNP's ALT call needs every changed base covered and none read as another base; one
    below `--min-baseq` weighs nothing (RJ-25), and `mnp_confirmed` needs all read.
  - A context row is one whose REF span meets X's span widened by the longest read, so a
    long deletion starting within reach (and running past it) is a context too.
- **The fit rule (phase B, decided 2026-10-25; RJ-29):** a read whose best score over every
  candidate sequence is below ten confident mismatches fits none of them; it is depth only
  where its aligner also doubts it (XS ≥ AS, or MAPQ 0 with no XS), otherwise its margin
  decides it. `UNEXPLAINED_READS(n/m)` names how many of a row's depth reads fit none. The
  other options measured stay behind `GBCMS_JUDGE_FITS` (measurement only).
- **Second fix round (2026-10-26), each red-first:**
  - A substitution's read whose confident bases at X's telling bases show some of X's
    changes and REF's bases at others is another allele (RJ-28's components), whichever of
    REF and X the margin favours, or where they tie.
  - At a substitution X, a recurring guard indel with which REF reads as X over a read
    reaching X's site (a repeat X lengthens, written as more units beside X) is left out
    where no reads need it: the single-item, REF-to-X part of the background refinement,
    for indels only.
  - A tie of REF with another allele, X beaten, is REF where the read writes an indel of its
    own across that allele's change (as where it never reaches it).
  - Measured and dropped: an indel on an MNV's last base taken as at X's site (its REF
    losses were an oracle v3 artefact); containment for states that change the length (a
    carrier of a listed delins ending in X's base is that delins's, oracle v3); and X with
    a listed change on a tandem duplication's inserted copy as a candidate (23 true-REF
    reads lost at HG002's 1 kb+ insertions, truth v2db; open for the operator); and a tie of
    X with another allele the read never reads apart from X credited to X (it credited
    reads the oracles call either).
  - `GBCMS_JUDGE_NEITHER=other` (measurement only; default `margin`): a covered telling base
    read confidently as neither allele makes a REF or ALT call on an MNV or a delins another
    allele.
- **Third measurement round (2026-10-27; measurement options, defaults unchanged):**
  decided 2026-10-27: margin; a confident third base below the found floor is an error (RJ-26);
  the per-read forms' gains come from the oracles' strict core. Beside-X: the default stays
  `site`; `beside` and `apart` are not adopted.
  - The confident third base (`GBCMS_JUDGE_NEITHER`). Mechanism: a read is REF only where its
    bases show REF's where the alleles differ, so a confident base neither allele has there
    makes a REF or ALT call neither. Errors put such bases on telling bases at the reads' own
    mismatch rate (on the patient arms, single-read neither bases per telling base at or below
    the mismatch rate elsewhere in the same reads), so the switch moves reads in proportion to the bases read:
    about 0.1% of calls at MNVs, 1% at delins of 20-49 bases, 2% at 50 or more; and where
    REF's molecules hold bases X's never do, the moves fall on REF and raise VAF (0.2-0.4
    points at 20 bases or more), RJ-24's mechanism. Forms: `other` (every telling base);
    `other:N` (alleles under N bases: 20, oracle v3's own core cut, is no mechanism);
    `recur:K` (where K other reads show the base: reads, not fragments, and two reads are
    mostly one fragment's mates; recurring bases are damage-like, C>A and C>T, as single
    ones are); `window` (only where both alleles' reads are read alike: a substitution's
    changed bases, not the bases an MNV keeps; a length change's junction windows, 4 telling
    bases from either end, as RJ-24 reads). Against `margin` (patient arms on oracles v2 and
    v3, HG002 truth v2dbq20 and the delins arm's het band): `other` fails delins of 20 bases
    or more on the panel, cohort, WES and gate arms; `recur:1` and `recur:2` fail 50+; `other:20`
    and `window` pass every arm on both oracles (`window` synthesised from the `other` lane,
    exact for the reads); HG002's truth arms change no class under any. No form breaks a
    decided contract. A base the read's haplotype holds, or the guard's background explains,
    is no neither base.
  - Beside an insertion (`GBCMS_JUDGE_ENDS`). `beside`: an allele removing bases whose
    placements end at an insertion X's point is a context. It moves 120 reads at 5 gate
    insertion rows from partial evidence to REF, the listed deletion's carriers (oracle v3
    calls them sibling carriers: 108 more wrong calls), and on HG002 loses ALT reads at
    insertions of 50-299 bases. Its general form, `apart` (a length change with a placement
    clear of every placement of a length-changing X is beside it), breaks 24 decided
    contracts: an allele's carriers are another allele wherever it is written (RJ-5, RJ-6,
    RJ-7). The default stays `site`.
- **Fourth measurement round (2026-10-27; the general forms of the second round's two
  narrowed fixes, measurement options, defaults unchanged, open for the operator):**
  - The background respelled (`GBCMS_JUDGE_BACKGROUND=respell`; `GBCMS_JUDGE_RESPELL_TRACE=1`
    names each item weighed). Mechanism: a context with which one candidate's sequence spells
    another's over the reads leaves those reads unable to tell the two apart. General rule:
    any background item or found allele (the masks excepted, decided), alone or joined with
    others into one change at X's site, in either direction, for any state and any X, is left
    out where no read needs it. Against the default: every patient arm passes on both oracles;
    HG002's truth arms gain (c246's 1-9 bp insertions 34 ALT and 15 REF reads, tr's 10-49 bp
    10 ALT, the 1 kb insertions 5 ALT and one wrong call fewer), no class failing. Nearly all
    of it is the step to indel and delins X (HG002: 289 rows, 417 reads changing call, 76 of
    them right by truth after and 2 before); the guard-indel rule's own items change nothing,
    and substitutions, other states and the state-to-REF direction move 5 reads at most. No
    decided contract changes.
  - The own gap (`GBCMS_JUDGE_OWNGAP`). Mechanism: a read that reaches a state only through an
    indel of its own at that state's change does not carry it. `all` (every state, wins and
    ties) breaks four decided contracts and fails classes on every patient arm: dropping REF
    or X in their own tie hands the read to the other by elimination, though its bases show
    neither where they differ. `others` (another allele's states only) breaks none, but fails
    delins under 20 bases and gate insertions on oracle v3: where another allele wins by the
    margin, a read writing an indel of its own beside that allele's change still carries the
    allele's bases (25 of 26 such panel reads wrong). Neither replaces the decided tie rule.
  - decided 2026-10-27 (operator): respell is the default, replacing the narrow guard-indel
    rule; the tie fix stays narrow, `all` and `others` off.
- **Decided-rule fix round (operator, 2026-10-28; each red-first, class tests in
  `tests/test_kept_base_change_contract.py`, `test_listed_sibling_floor_contract.py`,
  `test_contig_edge_judge_contract.py`, `test_beside_site_indel_contract.py`):**
  - *A change at a base X keeps* (REF's and X's alleles alike there, inside the bases X
    changes; for a delins, aligned from either end) that the REF molecules carry too is a
    separate event, not another allele (RJ-25): it is a context lying with REF, and X with it
    is a state of X's class, so carriers holding it stay ALT and REF molecules REF. It
    qualifies where it recurs on both haplotypes (3 fragments and 5% of the REF molecules,
    and of X's, or a state spelling X with it), or where the guard's background holds it and
    the input does not list REF with it (a co-annotated sibling stays one, RJ-5). Until now
    the judge took X's span from its first to its last change as X's own bases: such a
    change was dropped as meeting X, a mask there taken for X's base, and an ALT call needed
    X's base there while a REF call needed none. A base X keeps now tells nothing on either
    side, and a mask there holds on X's haplotypes as on REF's.
  - *A co-annotated sibling of X's length at X's site* is a state only where the reads'
    own alignments spell it in full on at least one fragment and on 5% of those spelling
    the allele it is nearest (RJ-25, 2026-10-14): a component never shown apart is part of
    the event, and at depth one fragment is an error's share.
  - *A pure indel whose placements touch a contig edge* (its repeat running to the edge, or,
    operator 2026-10-29, the indel itself at the contig's first or last base: a MAF insertion
    at Start 0, the base-after form at position 1, and an insertion after the last base): a
    read reaching the edge can only be called for the longer allele (the shorter fits it as a
    read ending early), a call resting on the shorter allele holding no base past the edge,
    where a circular contig's runs on (RJ-3). An ALT call there (an insertion's) is depth
    only, counted and warned, as under the rules; reads alike share the judge's verdict and
    the edge is weighed on each. A REF call there (a deletion's) stands: two decided tests
    count it REF (`test_contig_end_n_read_contract`, the MAF deletion at Start 1), where the
    rules' exact-carrier route makes it depth (a run of six A's ending the contig, five REF
    reads reading to the edge: (RD, DP) (0, 5) under the rules, (5, 5) under the judge);
    carried to #332 (J1) with these numbers, no change here (operator, 2026-10-29). At a MAF insertion at
    Start 0 the REF reads tell nothing either (they fit the ALT placed past the insertion), so
    the row is (RD, AD) (0, 0) under the judge, against the rules' (6, 6) that develop's
    #219 test asserts. Decided by design (operator, 2026-10-29), as the two delins
    germline-change tests: the test stays as it is, its expectation changing at the default
    flip (#331).
  - *An indel the reads write beside X's site* (two reads and 5%, every placement clear of
    the site, the repeat X changes included) is a context, a separate event (RJ-8), however
    a read's aligner wrote it, where the reads need it for their calls (as a respelling
    over a read is, 2026-10-27: with it, two reads and 5% of those reaching it are REF or
    ALT that are not so without it). Taken without that test (the first measurement of this
    round), it cost HG002 14 true ALT reads at insertions of 1 kb or more and 6 at 50-299
    bp: the event's own reads write such indels (a 1,795-bp insertion at chr16:34590143,
    whose reads write a 1-bp insertion four bases before it on 639 of 2,393 reads, absent
    from GIAB), and as a context it tied 13 carriers (ALT -15 against REF -19 without it,
    -12 against -12 with it). Known partial: where no read writes it beside the site, the
    judge has no candidate holding it (placing a read's own indel by its bases, RJ-27, is not
    yet the judge's).
  - *The other allele's base against a third base (`GBCMS_JUDGE_REFBASE`, measurement
    only, default `margin`; operator approved measuring it 2026-10-28; class tests in
    `tests/test_ref_base_option_contract.py`):* at a base an MNP or a delins changes, a
    read's confident base that is the other allele's (at a delins, the other allele's base
    at the same telling index counted from either end, where the read's own gap may put
    it) is that allele's evidence, not an error, so `partial` makes such a REF or ALT call
    another allele (RJ-25: an allele lacking part of the event; RJ-28 at an MNP, where the
    default already does so); a third base, neither allele's, is an error the margin
    forgives (RJ-26), so `partial` keeps an MNP carrier with one ALT instead of depth only.
    Bases X keeps tell nothing under both. The evidence it answers: a mate test of 72
    patient delins ALT reads the margin forgave (19 showing REF's base; under 20 bases, 10
    of 11 decided carry the odd base on both mates). Measured against the default (oracle
    v3): delins under 20 bp gain (panel +23 reads, gate subset +22, wrong calls −23 and
    −22, mostly REF reads showing the ALT's base, which v3 calls another allele); delins of
    20 bp or more lose a few (panel −3, WES −3); MNVs gain wrong calls (panel +6, gate subset
    +28: carriers with a third base at a changed base, which v3 calls another allele); HG002
    delins: 21 reads on 19 rows move (of the informative ones, 1 right and 4 wrong by phase),
    the het band at 20-49 bp 9 → 7 rows, no truth-arm read; a patient normal's 3>1 delins
    (RD, AD, partial) (828, 1, 5) → (827, 0, 7). Decided (operator, 2026-10-29): the default
    stays `margin`; the switch stays, off by default; a delins-only form (`partial` at a
    delins, `margin` at an MNP) is measured at the default flip (#331).
  - *Decided (operator, 2026-10-28): a germline change at a base a delins keeps.* Carried on
    both haplotypes, it is a separate event (RJ-25), and each read counts by the bases where
    REF and ALT differ: under the judge a delins carrier with the change counts ALT and a REF
    read with it counts REF, as at an MNP (`test_locus_alleles_contract.py`: the delins test
    `test_a_germline_change_inside_a_delins_the_carriers_carry_too_leaves_ref_reads_ref`
    gives (RD, AD) (46, 16) under the judge, the rules' exact reading (26, 6); its MNP twin
    `test_a_germline_change_the_carriers_carry_too_leaves_ref_reads_ref` expects (40, 10) on
    both). The delins expectation updates at the default flip (#331); the rules path keeps it.
  - *Carried to #332 (J1), no change here (operator, 2026-10-29):* at an MNP an ALT call with
    a confident third base at a base X changes is depth only (the phase C rule above) while a
    REF call with one stays REF, and a delins forgives one either way (RJ-26, round 3). On
    CC>GT and GCA>ACT, 20 REF reads and 10 carriers, one more read each: a carrier with a
    confident third base at a changed base is partial under the rules and depth under the
    judge; a REF read with one is depth under the rules and REF under the judge; at CAGTC>GAGC
    both judges count such a carrier ALT and such a REF read REF.
    `GBCMS_JUDGE_REFBASE=partial` measured one rule for both (above).
- **The background refinement (measured 2026-10-26, not adopted; `GBCMS_JUDGE_BACKGROUND`,
  default `asis`):** `respell` leaves out a context item of the guard's background or the
  found alleles that respells the site's own candidate sequences: one that turns one
  candidate into another whole (chr5:44868216 A>AT, X with a recurring 1 bp T deletion
  spelling REF's SNV haplotype), items that together spell one change meeting X
  (chr16:23024630 TA>AC, a 7 bp deletion and an inserted C spelling a carrier of X with a
  6 bp deletion beside it), and one that respells only over a read where no reads need it
  for their call. The guard's masks stay: removing them (`off`) breaks 26 decided-rule
  contract tests (a masked germline SNP whose haplotype reads as X never counts ALT), and
  a duplication's carriers cannot be told from such a SNP's at its one telling base. Against
  its pre-registered bar it kept every class and every contract test but recovered under
  half of what removing the background gains: none of the gate set's 36 bp insertion rows
  (INS 20-49 +7.7 points with `off`), 28% of the HG002 insertion reads, 21 of 442 gate
  wrong calls. It stays behind the switch, off. (The fourth round measured it as the general
  form of the guard-indel rule; it was adopted as the default on 2026-10-27, `asis` and `off`
  kept for measurement.)

Measured (HG002 phased truth, R2): see #304. The prototype is `plan670/r2/hapjudge` (local).

## 2. Where it enters
| Seam | Where | Change |
|---|---|---|
| Flag | `cli.py:387/755` → `models/core.py:159` (`AlignmentConfig`) → `pipeline.py:887 _engine_kwargs` → `#[pyo3(signature)]` `engine.rs:924/985` → `count_bam_binned_core` (385) → `count_bin_shared` (1076) → `count_variant_from_cache` (1477); stub `_rs.pyi:225–296` (parity test) | `read_judge: "rules" \| "sequence"` (`"haplotype"`, the former name, taken as `"sequence"`), plumbed like `alignment_backend`; also the MNP-rescue recount (`pipeline.py:1340`) and Nextflow (`dna/main.nf:60`) |
| Rows in reach | prep, a whole-input post-pass after `assign_listed_neighbours` (`normalize/engine.rs:131`) | new `Variant::reach_rows`: every other PASS row whose REF span lies within X's ± `READ_REACH` (300, as `carrier::WHOLE_READ_REACH`); computed also for the decomposed twin (1066) and the MNP-rescue rows (prepared apart, `pipeline.py:1306`) |
| Reference | prep, beside `whole_ref` (1104–1105) | `Variant::judge_ref`: X ± (`READ_REACH` + the widest reach row), for every kind (SNVs and deletions ≥ 50 bp have no `whole_ref` today) |
| Classifier | `engine.rs:1662`, after `splice_skip_triage` | under the flag, `judge::classify` returns a `ClassifyResult`; the sibling guards (1691, 1699) are bypassed (the judge weighs the site's alleles itself) |
| Background | the guard, forced before the loop (`engine.rs:1540`; `carrier::guard_for` 2655) | `Guard.background` → context rows; `Guard.others` → site states; for SNVs (no guard) the site holds REF, X and the listed alleles |

Nothing else changes: the read cache, fragment evidence (`fragment.rs:147/258`), RD/AD/DP,
strand, mFSD and diagnostics consume the `ClassifyResult` as today.

## 3. Mapping a decision to `ClassifyResult` (`utils.rs:70`)
- **REF / ALT:** `is_ref` / `is_alt`; `qual` = the lowest quality among the read's bases where
  the winning and the best losing candidate sequence differ (the bases that decided it; at least
  `--min-baseq`), so mate resolution (`resolve`, `fragment.rs:258`) weighs a decision by the
  bases it rests on. `is_structural`: `false` for every call (decided, §10.3).
- **Partial:** neither; `partial_match_count = 1` (another site allele won).
- **Depth only:** neither, `uninformative = true`.
- **Admission:** `covers_locus` from splice triage as today; `clip_admissible` (DNA) where a
  soft clip carries bases over X's span (RJ-13, clipped bases are the read's);
  `read_from_right` where the read's first readable base lies past X's left flank. A read is
  admitted to DP when its aligned span, clips included in DNA, reaches X's REF span; for a
  deletion of 50 bp or more, its anchor (RJ-24): reads wholly inside the deleted bases are not
  admitted, as today.

## 4. The scorer and cost
- Semi-global affine DP (read global, candidate sequence local) with a wildcard for masked
  bases, banded around the read's aligned diagonal mapped into each candidate sequence, the
  band widened by its length change. Per read, candidate sequences are trimmed to its reach.
- **Triage first:** a 31-mer index of the variant's candidate sequences; a read whose fully
  confident k-mers fall in one site state's sequences only, with no unknown k-mer over X's span ± k,
  is decided without alignment (R2 round 1: ~40% of reads).
- **Target:** within 2× develop's runtime on the panel, ACCESS (deep) and WES arms. The
  prototype is ~300× develop; SIMD banding (e.g. the `block-aligner` crate, which takes
  position-specific scores for the mask) is the next lever if the plain banded DP is not
  enough. `wfa2lib-rs` (in-tree) aligns ends-free and affine but has no wildcard.

## 5. Contexts and states
- Rows in reach ≤ 4 (outside the site): every non-overlapping subset is a context.
- More: singletons, pairs a read shows together (per-row evidence on ≥ 2 reads), then
  read-supported sets, to 32 candidate sequences; the empty context always.
- States: REF, X, the site's listed alleles, and for an SNV every other base; for an MNP or
  an indel the guard's discovered alleles there. A state holding X's change is of X's class
  (RJ-28); a substitution in a long deletion's bases past its junction windows is a context
  on REF's side (§1).
- Deterministic order everywhere (rows by position; ties to fewer rows); per variant only.

## 6. RNA
Spliced reads are judged on candidate sequences built along their mRNA: the read's segments
(`carrier::splice_segments`, 1471) give the spliced reference, rows mapped through
`SpliceMap` (1364); a splice through X's site leaves the read depth only (as `carrier.rs:726`).
Strandedness and antisense handling are unchanged (engine.rs:1637/1721/1924). Measured with
package 3 on FORTE.

## 7. What it replaces, and what stays
Under the flag, `check_snp`, `check_mnp`, the exact-carrier rule, `junction::classify` and
`classify_complex` are not called for the row's own counts; the guard still runs (it supplies
background and discovered alleles). Per-transcript counts (3096) and ASJD (3704) call the same
`judge::classify`, so every output agrees. `Tolerance` and the error census are not used by
the judge (a recurring change becomes a context row instead).

## 8. Tests
- `tests/test_haplotype_judge_contract.py` (red-first): #246 in trans and in cis, a compound
  het, a repeat unit count, a low-quality telling base, no phase credit, a soft-clipped long
  insertion, admission at a long deletion's anchor.
- Every counting contract runs with the flag too where its shape is the judge's
  (parametrized); `tests/census.py` stays the default rules' check — the judge's independent
  checks are the R2 truth arms (HG002 phased truth) and the BWA read checks on local arms.
- Binning invariance (`tests/test_binning_invariance.py`) with the flag: every geometry, one
  call per row, decoys in the call.
- Rust unit tests: the scorer (wildcard, ends-free, affine), the construction, the decision.

## 9. Acceptance (before the default flips)
Every arm against develop and truth: the HG002 strata and R2 arms (#246 pairs, tandem
repeats, SVs), the compound hets; panel, ACCESS and normals, WES, the gate set against
sign-out, FORTE (with package 3); the IMPACT 50–500 bp events read by read (ABRA2 op, `YO`
placement, bases). Cost within 2× develop. The HPC gate.

## 10. Decided and open points (operator, 2026-10-23)
1. **Decided:** a judged call's `qual` is the lowest quality among its deciding bases.
2. **Decided:** a third base at an SNV is partial evidence: SNV sites hold all four bases as
   states; MNP and indel sites take the guard's discovered alleles.
3. **Decided (operator, 2026-10-29):** `is_structural` under the flag is retired: every call
   rests on its bases (invariant 7: judge bases, not placement), and mates resolve by
   quality. Keeping it for reads whose CIGAR writes X's exact op (`GBCMS_JUDGE_STRUCTURAL=keep`,
   removed in #363) was measured on every arm and not adopted: read counts identical
   everywhere; fragment counts changed on 40 gate-subset rows, all 1-bp homopolymer indels
   (ADF +76, RDF −8; 82 fragments there have one mate ALT and the other REF, which `keep`
   resolves to the ALT mate's structure), and on 7 HG002 rows (ADF +9); no truth set scores
   fragments.

## 11. Open tickets the judge resolves (measured 2026-10-28)
On the merged branch (develop 6efeb5bb merged in), against its own default path, which gives
develop's counts. HG002 (public) and synthetic shapes; HG002 truth is v2dbq20 (reads whose own
bases at the run's `--min-baseq` decide them). Each is resolved for users only with the default
flip (#331): the default rules path keeps today's behaviour.
- **#246** (an SNV beside an insertion credited with the insertion's carriers): resolved. The
  listed insertion is a context in every candidate sequence, so a carrier whose first inserted
  base an aligner wrote on the SNV's base fits REF at the SNV. #246 arm (4,578 rows), SNV class:
  REF reads called ALT 243 → 0, ALT reads called REF 322 → 0, true ALT credited 92.9% → 99.3%,
  true REF 94.9% → 99.5%, Σ|AD − truth| 1,702 → 734; insertions 1-9 bp true ALT 59.6% → 98.6%.
- **#144** (clip-borne duplication carriers): resolved on HG002. Clipped bases inside the
  fragment are the read's own (RJ-13) and admit it, and carriers written through count by their
  bases. The 89-bp duplication: AD 2 → 7 (bases-only ceiling 10); the 98-bp: 6 → 16 (ceiling 16).
  The panel's ITD rows were not read read by read.
- **#290** (another length of the repeat written upstream counts REF): resolved by its decided
  rule (another allele wherever written). The found alleles make the expansion a state of the
  site: its synthetic counts 6/6/6 (RD/AD/partial) written 0 to 60 bases upstream, where the rules
  give 6/6/6 up to 30 bases and 12/6/0 at 60. The five strict xfails pinning RJ-27's named limit
  pass under the judge; in the 40-bp VNTR compound het the rules call 22 of the insertion
  haplotype's reads REF, the judge none.
- **#301** (a 33-bp deletion row counting a co-located 22-bp deletion's reads): resolved. The
  other deletion is a state whether listed or not. Against a per-read bases oracle at the site
  (12 reads the 33-bp allele's, 14 the 22-bp's, 21 ties, 3 fitting none): with the 33-bp row
  listed alone the rules count 11 of the 22-bp allele's reads ALT, the judge none (AD 25 → 14,
  12 of them its own); the counts no longer depend on which rows are listed.
- **#297** (a long deletion's unit count in a repeat): partly. Its synthetic (three 60-bp
  units; 4 REF, 4 row and 4 two-unit carriers of 100 bp): units differing at one or two bases
  give 4/4/4 under both judges (since #285); identical units give 4/8/0 under both, the
  two-unit carriers counted ALT where no read shows the other unit count; with 4 more two-unit
  carriers of 190 bp that read past the repeat, the rules give 4/12/0 and the judge 4/0/12.
  Open for the operator: the identical-units case with no read past the repeat (count the given
  allele, or depth only as #297 proposes), and the phased-change shape of
  `test_a_change_the_input_lists_beside_the_deletion_is_its_neighbour`, which the judge counts
  ALT (no read writes the two-unit allele).
- **#264** (REF and VAF of indels in repeats): the long-insertion shortfall it handed to #304.
  HG002 het median VAF (rows in 0.35-0.65): tandem-repeat insertions 10-49 bp 0.471 (18/41) →
  0.500 (23/44), 50-299 bp 0.238 (6/34) → 0.429 (17/41), 300-999 bp 0.200 → 0.368, 1 kb+
  0.241 → 0.348; SV-subset insertions 50-299 bp 0.364 → 0.474, 1 kb+ 0.125 → 0.444. Its
  decided documentation (allele-level REF and VAF) is still open on #264.
