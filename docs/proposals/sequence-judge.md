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
  - Beside an insertion (`GBCMS_JUDGE_ENDS`, since removed). `beside`: an allele removing bases whose
    placements end at an insertion X's point is a context. It moves 120 reads at 5 gate
    insertion rows from partial evidence to REF, the listed deletion's carriers (oracle v3
    calls them sibling carriers: 108 more wrong calls), and on HG002 loses ALT reads at
    insertions of 50-299 bases. Its general form, `apart` (a length change with a placement
    clear of every placement of a length-changing X is beside it), breaks 24 decided
    contracts: an allele's carriers are another allele wherever it is written (RJ-5, RJ-6,
    RJ-7). The default stayed `site`, and the switch is removed.
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
  - *A pure insertion of 50 bp or more* (RJ-23, operator 2026-10-30, #264) is read over
    equal windows from the anchor, as the rules read it: a REF or ALT call holds, on its
    allele's sequence as the read aligns to it, the anchor and as many bases on as REF's
    window holds (the run the insert slides over and the base past it); a read inside the
    insert, or holding only its far junction, holds none and is depth only. Where the guard's
    masks hide the base telling the windows apart (a change past the run lands, on X's
    haplotype, past the insert), both read on to the next base telling them apart and a read
    short of it is depth only, REF and ALT alike (RJ-27's read-on); a carrier ending at a
    change in cis with the insertion stays depth only (RJ-27's named limit). Reading such a
    change on REF's haplotypes as the reference's where REF's molecules show they lack it was
    measured with a temporary option and not adopted (REJ-20261023-001); the option is removed.
    `GBCMS_JUDGE_LONG_INS=window` keeps the windows without the read-on and `off` restores the
    judge before them (measurement only).
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
- **Sixth measurement round (2026-10-28/29; the `--min-baseq` sweep and two measurement
  options, measured on `feature/c46-found-floor`, not merged; one fix adopted):**
  - The sweep (the judge and the rules path at `--min-baseq` 0, 10, 20 and 30 on every arm,
    each scored on oracles v2 and v3 and HG002 truth v2db built at that threshold; ACCESS XS1
    and XS2 apart). A base below the threshold is unknown and costs nothing (decided). The
    judge's wrong definite calls fall as the threshold rises (oracle v3: panel 696, 653, 365,
    157; gate 2001, 1689, 1160, 519): at 0 and 10 a Q2-Q19 base counts as a confident one, so
    one low-quality error can decide a read (gate MNVs: +532 wrong calls against the rules
    path at 0, +135 at 20). ACCESS XS1 and XS2 read the same at 0 to 30 (their bases are Q30
    or more, their unresolved ones N); HG002 the same at 0 and 10 (no base between Q2 and Q11).
    At 30 far fewer reads are decidable on IMPACT (31% fewer on the panel, 45% on the gate).
    Against the rules path the judge agrees with the oracles on far more reads at every
    threshold (panel about +9,600, gate +14,000), with fewer wrong calls overall; the classes
    over the bar's wrong-call limit are as at 20 (delins, MNVs, some 2-49 base indels).
  - The quality-scaled cost (`GBCMS_JUDGE_COST`; not adopted, operator 2026-10-29): a base
    at or above `--min-baseq` costing 4 × min(Q, C) / C, so a read whose telling base is below
    the cap C is depth only. At 20 it fails the bar on every IMPACT, WES and normal arm and on
    HG002 (cap 30: 39,500 agreeing reads lost on the panel, 343,000 on the gate, for 2 and 89
    fewer wrong calls; HG002 loses true ALT and REF reads in every class); cap 40 breaks 177
    decided contract tests. The flat cost stays.
  - The quality-aware found floor (`GBCMS_JUDGE_FLOOR`; the tested forms `tail`, `llr` and
    `learned` not adopted, operator 2026-10-29): a spelling is another allele when it recurs
    on more fragments than errors at its bases' qualities allow (LoFreq's exact tail, Mutect2's
    log odds, or change rates learned beside the event). No form passes the bar at 20. They
    find more alleles (gate IMPACT 248 spellings at 122 rows against 159 at 90) and gain where
    the residual lives (ACCESS XS1 delins −14 wrong calls, gate MNVs −20), but they act on
    every recurring spelling near ALT or REF: true ALT reads become partial at IMPACT indels
    read at Q29-Q40, and at a GC-repeat insertion multi-change spellings become alleles (gate
    INS 2-19 repeat +18 wrong calls on v3). `llr:30` is the best of them and fails the same two
    gate classes; the learned rates change one allele only and give back that XS1 gain. The
    share floor (RJ-25) stays.
  - **Decided (operator, 2026-10-29): the judge writes the alleles the reads spell past the
    guard's masks (RJ-25, amended; class tests in
    `tests/test_found_allele_spelling_contract.py`).** Mechanism: the guard finds the other
    alleles the reads spell at the locus, reading each spelling past its masks (a masked
    base, recurring background, is no difference), and wrote each over its window with its
    carriers' bases. The judge took that row as a site state, so a masked base became the
    allele's change, read literally for every state: a read carrying the masked change fit the
    found allele better than X or REF by that base alone. Now the judge takes each allele as
    the bases it changes, past the masks (the allele's own base where a mask takes a
    carrier's), so the masks stay on every state; one that so written is a candidate the site
    already holds (X with a listed neighbour) is that candidate, and spellings apart only at a
    mask are one allele (`Guard::judge_others`; the rules' windows keep the window spelling,
    which they read past the masks already). Written as its changes, such an allele no longer
    spans the guard's window, so the check that drops an other state spelling X (or REF) with
    listed neighbours reads every neighbour that, or one of whose equivalent placements, meets
    the state's place, not only those wholly inside it. Two HG002 sites showed why (each a
    regression of a first build, found in the measurement): at chr10:12910609 A>ATTT (an
    insertion where an A run meets a T run, listed A>T and AT>TT on the next two bases) the
    allele, X with both, is AA>TTTTT, which the AT>TT runs one base past (read inside only, 15
    carriers tied X with it, 13 ALT by the truth set); at chr10:83775882 A>AC (at
    `--min-baseq` 0; a C inserted in an A run after CAAAA units, listed T>TA before the units
    and A>C on the next base) the allele of REF's length, REF with both, starts at the units'
    first C, while the T>TA is written before it and slides into it (15 REF reads, REF by the
    truth set, tied REF with it). With the wider check both rows are develop's.
    - Worked example (a patient panel tumour, the same run in the panel and twice in the gate
      set; read by read, operator 2026-10-29): X a 17-base tandem duplication; a listed SNV
      seven bases past the copies on every carrier (111 of 111 ALT-spelling fragments), and
      three recurring unlisted changes just past it, masked. The guard's other allele was X
      with that SNV, spelled with a carrier's base at a mask: under the old writing every
      carrier fit it at least as well as X, and 293 carriers per run were partial (ALT by both
      oracles); written past the mask it is X with the listed SNV, a candidate already, and
      they count ALT (panel INS 2-19 repeat agreement 58.1% → 99.65%). The bar's 2 more wrong
      calls on v3 (6 on v2) at that event are the oracles' limits: two reads carry the
      duplication with one confident substitution inside the first copy, each on one fragment
      (an error, RJ-26; v3's clean test allows no mismatch in the 18-base event window), and
      four reads v3 calls ALT, as the judge does, while v2 ties them with a sibling (it holds
      no ALT haplotype with the listed neighbours). Those four (and 93 more per run both
      oracles call ALT) do not read through the duplicated tract, so they cannot tell two
      copies from three: the spanning question, #371 (under #332), not this rule's.
    - Class tests (red-first; each fails on develop): every read carries a germline change one
      base before a 7>4 delins (masked). X's carriers all carry a listed SNV beside the event:
      (RD, AD, partial) (60, 0, 30) → (60, 30, 0). REF's other allele beside the mask:
      (0, 100, 300) → (280, 100, 20). X's other allele beside the mask: (300, 0, 100) →
      (300, 80, 20). Unit tests: the allele trimmed to its change, past the mask; two spellings
      apart only at the mask one allele; the two HG002 shapes above, an other state that is X
      with neighbours running past its changes, and one that is REF with a neighbour written
      before it, dropped (3 state graphs under the narrower checks, 2 now).
    - Contracts: the judge-switch suite (`GBCMS_TEST_READ_JUDGE=sequence`) keeps develop's 139
      baseline failures with identical assertion values; the default suite passes.
    - Measured against develop on every arm (panel, cohort IMPACT, ACCESS XS1 and XS2 apart,
      WES, the normals, the gate set; the branch head after merging develop 7e7ea3c7), both
      oracles at the run's threshold and HG002 truth v2db at that threshold, at `--min-baseq`
      0, 10, 20 and 30:
      - The writing moves one event, at 20 only: the worked example's carriers, partial → ALT
        (panel +240 agreeing reads on v3, INS 2-19 repeat 58.1% → 99.65%; the gate set, which
        holds the run twice, +480); the row's (RD, AD, partial, VAF) (375, 0, 250, 0) → (375,
        247, 3, 0.40), fragments (RDF, ADF) (306, 0) → (306, 199). Elsewhere it moves under 20
        reads per arm (gate MNVs: 9 reads at 5 rows, 6 of them REF reads made partial). Cohort
        IMPACT, ACCESS XS1 and XS2, WES and the normals are unchanged at every threshold.
      - The wider spelled-state check moves 4 rows in the panel's runs (5 in the gate set) at
        every threshold, all where the input lists one molecule several ways: 247 panel reads
        partial → ALT at 0, 10 and 20, 153 at 30; ALT by v3, a listed sibling or "either" by
        v2 (the item below).
      - The bar as registered (>2 more wrong definite calls, or >0.5 points and >2 reads fewer
        agreeing): at 20 v3 fails only the worked example's event (gate INS 2-19 repeat +4
        wrong calls, 2 per run: the two one-substitution reads), and v2 fails it (+6 per run)
        and DEL 2-19 unique (the multi-listed rows, −1.24 points on the panel, −0.97 on the
        gate's IMPACT runs); each read above as the oracle's limit. At 0, 10 and 30 v3 passes
        every class on every arm, and v2 fails DEL 2-19 unique alone (the same rows). Every
        other class on every arm passes on both oracles.
      - HG002 passes at every threshold: at 20, c246 INS 1-9 +2 ALT reads and INS 1k+ +1 ALT
        and +1 wrong call (nocl, sub); at 0 and 10, c246 INS 1-9 +29 ALT and +3 REF reads, tr
        INS 10-49 +5 ALT; at 30, INS 1k+ +3 ALT. Every read with a truth label moves toward it
        but one (a REF read at c246 INS 10-49 at 0, partial → depth). The two first builds'
        regressions (chr10:12910609 at 20, chr10:83775882 at 0, above) are gone.
      - Combined check with develop's path-alignment change (#370, develop a2ca57c1) merged
        in, at 20 (batch c61hs): on the judge path the merged head moves the same reads
        against develop a2ca57c1 as the branch did against develop before it, read for read
        (panel 572, the gate set 944 of its 947: of the other 3, path-alignment already makes
        two partial evidence, as the change did, and the third stays depth only; HG002 c246 20, tr 8, ch_a1 and ch_both 8 each, sub 22 and nocl 28), and
        one more: an HG002 read at chr11:11247039 G>GTCTC (no truth label there), partial
        evidence → depth only. The scores and the HG002 classes are those above; the
        shared-indel note (INDEL_SHARED) is on the same 8 panel and 10 gate-set rows. On the
        default path the merged head is byte-identical to develop a2ca57c1 on every patient
        arm, the gate set and the 8 HG002 arms.
      - Identity: on the default path (`--read-judge rules`) the branch head is byte-identical
        to develop 7e7ea3c7 on every patient arm, the gate set and the 8 HG002 arms (1,070
        runs); on the judge path the merge of develop changed nothing (the merged head before
        the placement fix, against the build measured before the merge: 0 runs differ), and
        the placement fix changes nothing at 20 (0 runs differ), only the HG002 rows above at 0
        and 10.
      - Known limit, outside this rule's class (#374, J1e, under #332): at 0, 10 and 30 the
        worked example's carriers stay partial under develop and under this rule alike, (RD,
        AD, partial) (376, 0, 250) at 0 and 10, (354, 0, 246) at 30. The guard's masks depend
        on the threshold: at 20 it masks six recurring changes beside the event, among them a
        G>A two bases past the listed SNV that every carrier carries (the ALT reads' own
        isolated change); at 0 and 10 it masks two, at 30 five, without that G>A. The allele
        the reads spell there is X with the listed SNV and that unmasked G>A, a spelling with a
        change of its own, not a masked base. The rules path, which this change leaves alone,
        counts that row's carriers partial at 20 too: (308, 0, 237).
    - Residuals at HG002 (each a read or two, the bar passing): at chr8:80072317, a compound
      het of a 12- and a 14-base TG insertion, the 14-base row's AD 12 → 10: the reads spell
      the 14-base allele with a G>A at a masked base on 3 fragments (RJ-25's floor: 3 and 5% of
      the nearest allele's), which written past the mask is that change alone, and the two
      reads carrying it are partial evidence, as RJ-25 decides; the window spelling had held
      the carriers' base at the mask, so it missed them. At chr13:112014990 and chr11:11246937
      (nocl, no truth label at those rows) a few reads move between partial, depth and REF.
    - **Decided (operator, 2026-10-30): the wider spelled-state check reaches listed siblings,
      RJ-28's indel twin (RJ-28 amended), named by `INDEL_SHARED`.** The
      check above drops any other state that spells X (or REF) with listed neighbours, not only
      the guard's alleles: a listed delins that is X with a listed neighbour is X's too. On the
      panel it moves 4 rows in 4 runs (247 reads, 5 rows in the gate set), all inputs listing
      one molecule several ways. The worked row (allele shapes and counts only; offsets from
      X's first base):
      ```
      REF       …T TCTGAT TCTG ACGAAGATGGTGAA …   (offsets −11 … +13)
      molecule  …. ...... TCTG .............. …   (. deleted: both deletions together)
      X         ACGAAGATGGTGAA>-          at  0   (a 14-base deletion)
      listed    TTCTGAT>-                 at −11  (a 7-base deletion)
      listed    TCTGATTCTGACGAAGATGGTGAA>CTG at −10 (a 24>3 delins: the two together)
      ```
      At X's site the delins is a state (it overlaps X). Its sequence is X with the listed
      7-base deletion beside it, a context, so it is X's; the narrower check read only the
      neighbours wholly inside the delins's place, and the 7-base deletion starts one base
      before it. Develop: X (RD, AD, DP, partial) (747, 0, 898, 129), ZERO_ALT and
      PARTIAL_DOMINANT, while the 7-base deletion counts AD 130 and the delins AD 125; this
      change: X (747, 129, 898, 0), the other two unchanged, so one molecule counts at all
      three rows (VAF 0.15, 0.14 and 0.14). Oracle v3 calls the 129 reads ALT; v2 a sibling or
      "either" (it holds no ALT haplotype with the listed neighbours). The rules path counts
      X (746, 3, 898, 126) there; at the three 3-base deletion rows of the same kind it counts
      48, 35 and 27 ALT reads where develop's judge counts none and this change 51, 35 and
      32. It follows the decided rule of 2026-10-25 for
      substitutions (RJ-28: a candidate sequence holding the row's change there, an MNV listed
      beside its component SNV, is the row's ALT, so a read carrying the MNV counts at both
      rows, `MNV_SHARED` naming the overlap) and the judge's own rule for a state spelling X
      with listed neighbours. Decided as its indel twin: a listed indel or delins that holds
      the row's whole change counts its carriers as the row's ALT, so one molecule counts at
      each row that holds it, and the row names, per such listed row, how many of its ALT
      reads spell it whole: `INDEL_SHARED(chrom:pos:REF>ALT:n/m)` in `gbcms_diagnostic` (no
      count, no column; red-first tests in `tests/test_indel_shared_contract.py`: the worked
      shape named at both deletions and not at the delins, and a listed deletion holding more
      than the row's change not named). An ALT read spells such a row whole when its
      haplotype holds X with the neighbours that make it and it reads every base from the
      base before the row's span to the base after it, confidently. Measured at 20 against
      the build without it: counts and read calls byte-identical on every arm and HG002
      arm; the token on 8 panel rows (the four rows and their listed deletions; at the
      worked row 118 of the 129 ALT reads spell the delins whole, at its 7-base deletion 117
      of 130), 10 in the gate set (five such rows and their listed deletions), none in the
      cohort, WES, normals or HG002.
    - Notes, not changed here: the rules path counts the worked example's carriers partial
      at 20 too ((RD, AD, partial) (308, 0, 237); RJ-25 there takes X with the listed SNV as
      another allele); and at an SNV listed right beside a delins its carriers all carry, the
      judge counts them partial on develop and with this change alike (seen while building
      the class tests).
    - Community standard: a variant is written by the bases it changes, left-aligned and
      trimmed (the parsimonious representation of [vt normalize](https://pmc.ncbi.nlm.nih.gov/articles/PMC4481842)
      (Tan, Abecasis and Kang 2015, DOI 10.1093/bioinformatics/btv112); GATK
      [LeftAlignAndTrimVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360036350152-LeftAlignAndTrimVariants)
      trims the bases common to all alleles; `bcftools norm`); freebayes writes haplotype
      alleles over a window and its authors point to decomposing them into allelic
      primitives ([vcfallelicprimitives](https://manpages.debian.org/bookworm/libvcflib-tools/vcfallelicprimitives.1))
      where pointwise alleles are wanted. A neighbouring change is its own record in every
      such tool; none documents a masked background base, which is gbcms's guard's own.
    - Effects: changes only under `--read-judge sequence` (the rules' windows read the window
      spelling as before, `Guard::others`; the rules path's outputs are byte-identical to
      develop's on every patient arm, the gate set and the HG002 arms, 1,070 runs and 8 arms).
      Under the judge every path builds a row's site through `judge::site_for` (the main and
      per-transcript counts and ASJD alike; the decomposed form and the MNP-rescue recount
      build theirs the same way), so all follow: the per-transcript counts and ASJD reuse the
      read's verdict (the judge refuses RNA, so RNA is unchanged), mFSD and the observations
      export (`count_bam_binned_observations`, the same core) take the calls, merge and the
      writers carry the counts (no column changes), and `UNEXPLAINED_READS` follows the depth
      reads. Nextflow and the CLI are unchanged (no new flag).
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
4. **Decided (operator, 2026-10-29):** a read's call is examined on the alignment its score
   came from. The call reads which bases sit at X's telling bases, whether the read holds the
   anchor, and the deciding quality off the read's alignment to the haplotype its best path
   spells. That was a second, linear alignment, banded its own way; it is now the path the
   graph programme traced (`graph::kept_steps`, `Site::alignment_on_path`). The linear aligner
   stays only where no state's path spells the haplotype or keeps its steps: `Site::reaches_state`
   (REF's haplotype with another state's contexts, for REF tied with another allele), and a
   path keeping no steps (a graph not laid out in order, which the builder never lays out, or
   an empty read).
   - **Mechanism.** On tied placements in a run both choose the same alignment (32,064
     synthetic alignments agreed). They part where the two bands differ (the linear band can
     miss the path the graph finds: a read beside a listed 14-base deletion in an 18-base
     repeat aligned at -20 against its state's -17) and where a deletion ties across repeat
     units (each aligner breaks the tie its own way, and the deciding quality followed the
     second one: 21 where X's base reads 25). Red-first tests pin both shapes and the tie.
   - **Effects.** Both read loops call the same `judge::classify`, so every output that
     consumes its `ClassifyResult` follows: RD/AD and fragment counts (mates resolve by
     quality), mFSD, the observations export and the MNP-rescue recount. RNA, per-transcript
     counts and ASJD are unchanged: RNA mode refuses the sequence judge. The writers,
     diagnostics, flags and Nextflow are unchanged in code.
   - **Measured on every arm** (develop dc16364c against this change; oracles v2 and v3;
     HG002 truth v2dbq20):

     | Arm | Reads changed | Oracle v2 (agree, wrong) | Oracle v3 (agree, wrong) |
     |:--|--:|:--|:--|
     | panel (IMPACT tumour) | 37 (2 SNV rows) | +36, +0 | +0, +0 |
     | across (IMPACT normal, ACCESS XS1, XS2) | 0 | 0 | 0 |
     | WES | 0 | 0 | 0 |
     | germline normals | 0 | 0 | 0 |
     | gate subset | 69 (18 rows) | +36 SNV, +2 MNV, -2 DEL 20-49 unique; +0 | +2 MNV, -2 DEL 20-49 unique; +0 |
     | HG002, 8 arms | 830 (317 rows; most in SV rows the truth leaves undecided) | | |

     HG002 truth v2dbq20: tr SNVs +1 ALT, no wrong calls. The decided-contract suite (the
     judge switch) is identical to develop's. Cost (Mac CPU, develop's judge path against this
     change, timed together): 2 to 10% lower on the patient arms (panel 10%, across 4%, WES 5%,
     normals 2%), 4% on the gate subset (3.19× to 3.07× the rules path), 11 to 28% lower on
     HG002 (no second alignment per examined state).
   - **Worked rows (HG002, public).** tr 20:19346218 A>G: one read scores ALT 0, REF -4,
     OTHER -4 either way; the linear alignment does not place its base at X (depth only), the
     traced path does (ALT); truth ALT. c246 17:80743837 A>T: one read scores ALT -4, REF 0,
     OTHER -4; depth only becomes REF; truth REF. On the panel, two SNV rows' 36 reads go from
     depth only to ALT (oracle v2 ALT, v3 ambiguous).
   - **Other tools.** GATK HaplotypeCaller genotypes from PairHMM likelihoods summed over every
     alignment of a read to each haplotype (no single alignment decides), its realigned output
     a separate Smith-Waterman alignment; freebayes calls from reads' literal sequences over a
     window; Strelka2's handling is not documented.
5. **Decided (operator, 2026-10-29): the five speed options.** Each was built as an
   off-by-default switch (`feature/c46-speed-opts`, not merged) and measured on every arm
   against the judge head: oracles v2 and v3, HG002 truth v2dbq20, the decided-contract suite,
   Mac CPU time.

   | Option | Speed | Counts | Decision |
   |:--|:--|:--|:--|
   | The call examined on the traced path (`GBCMS_JUDGE_PATHALIGN`) | −6% gate subset, −3 to −10% patient arms, −9 to −33% HG002 | small; agrees with oracle v2 and truth | adopted as the default (#370) |
   | No other-base states at an SNV for a read whose base at X is REF's or X's (`GBCMS_JUDGE_SKIPOTHER`) | 0 to −3% | where a third base ties X or REF the tie rule (RJ-25) no longer sees it: oracle v2 fails on wrong calls (panel +12, gate subset +38), v3 and truth gain | not a speed lever; a rule question (#367, under #332); the switch stays off |
   | The narrow band widened only by the rows the read crosses (`GBCMS_JUDGE_BANDROWS=crossed`) | none | HG002 tr SNVs +3 wrong calls (fails the bar) | rejected |
   | The floor from the band's own half-width (`GBCMS_JUDGE_BANDFLOOR=half`) | up to −5% | gate subset −22 agreeing, +4 wrong (v2); widens the band gap | rejected |
   | The full programme for every read (`GBCMS_JUDGE_BAND=full`) | ×2.75 gate subset, ×1.1 to ×2.9 elsewhere | exact; small gains | measurement only: the yardstick of the band's rule |

   The full programme showed the band gap: a band standing on a worse score than the read's
   best, on 2.3 reads in 10,000 on the gate subset, 32 on the panel and up to 40 on HG002.

6. **Chosen (operator, 2026-10-29): a band's floor counts the repeat a read may slide along
   (#369).**
   - **Mechanism.** A band stood on its edge alone for any score down to what an indel leaving
     it costs (16 bases each side: −23). A read its aligner placed a repeat unit off does not
     pay that: it slides a unit along the repeat at no cost and leaves the band on a short gap,
     so the band stood on a worse score than the read's best. The floor now takes the longest
     tandem-repeat period p in the read's reach (its span ± its length) off the slack: a band
     stands on its edge alone only down to the cost of an indel of slack − p + 1 bases, and
     only on a perfect score once p reaches the slack (`score::stand_floor`). A repeat is a
     stretch repeating at period p (to 64) for two periods and at least 20 bases
     (`score::tandem_tracts`); a long run counts its multi-unit shifts. It holds for the narrow
     and wide graph bands and for the linear aligner's band. Red-first tests: a period-13 read
     placed a unit off scored −19 in its band against its best −16 (now −16); 20 of 912
     synthetic shapes of the class fell short, none now.
   - **Effects.** As for the path-align decision: both read loops call `judge::classify`, so
     read and fragment counts, mFSD, the observations export and the MNP-rescue recount
     follow; RNA mode refuses the sequence judge; the writers, flags and Nextflow are unchanged
     in code.
   - **Measured on every arm** (develop dc16364c against this change on it, and against the
     full programme on this change; the judge cleanup since is identity-neutral):

     | Arm | Reads changed | Oracle v2 (agree, wrong) | Oracle v3 (agree, wrong) | Reads the full programme changes |
     |:--|--:|:--|:--|--:|
     | panel (IMPACT tumour) | 0 | 0 | 0 | 6 |
     | across (IMPACT normal, ACCESS XS1, XS2) | 0 | 0 | 0 | 0 |
     | WES | 0 | 0 | 0 | 1 |
     | germline normals | 0 | 0 | 0 | 0 |
     | gate subset | 2 (1 DEL row) | +2 DEL 1 repeat, +0 | +2 DEL 1 repeat, +0 | 177 |
     | HG002, 8 arms | 666 (401 rows) | | | 1,508 |

     HG002 truth v2dbq20: tr INS 1-9 +2 ALT, INS 10-49 +1, SNV +2; no wrong calls (the full
     programme: the same, and DEL 50-299 +1 on nocl and sub). The looser truths agree: v2db
     one wrong call fewer on INS 50-299 (nocl, sub); v2cq20 tr +17 ALT and 6 wrong calls fewer,
     c246 +1 REF. Of the 666 HG002 reads, 600 move to the full programme's call and 38 away
     from it; on truth-decided reads the two differ by 5. Most moves are REF to depth only (319,
     truth undecided): in a long repeat the read fits the other states as well one unit along,
     where it does not reach X. The decided-contract suite is identical to develop's.

     On develop with path-align (a2ca57c1) the floor moves the same 2 reads on the patient arms
     and 666 on HG002, 618 of them as on dc16364c. The others differ through path-align: 27 it
     already moves to the floor's call, 12 the floor moved only through the linear aligner's
     band (which no longer decides a call), 25 start from a call path-align changed, 14 move
     now because the call is read off the traced path the floor changes, and 9 move from or to
     another call. HG002 truth v2dbq20: tr INS 1-9 +2 ALT, INS 10-49 +1, SNV +1 (path-align
     makes the other); no wrong calls under any truth set.
   - **The band gap left** (reads the band scores below the full programme, per 10,000 judged;
     the band rule traced beside the full programme on one build, `GBCMS_JUDGE_FREESHIFT=off`
     for the floor before):

     | Arm | Reads judged | The floor before: short (per 10k), margin's call moved | This floor |
     |:--|--:|--:|--:|
     | panel | 421,745 | 1,352 (32.1), 43 | 1,343 (31.8), 43 |
     | across | 411,166 | 329 (8.0), 0 | 329 (8.0), 0 |
     | WES | 125,144 | 133 (10.6), 2 | 125 (10.0), 2 |
     | germline normals | 100,598 | 0 | 0 |
     | gate subset | 12,151,735 | 2,797 (2.3), 208 | 2,743 (2.3), 206 |
     | HG002 delins | 16,662 | 7 (4.2), 1 | 4 (2.4), 0 |
     | HG002 compound hets (3 arms) | 46,462 | 8 (1.7), 0 | 0 |
     | HG002 c246 | 130,386 | 307 (23.5), 197 | 133 (10.2), 98 |
     | HG002 tr | 97,328 | 267 (27.4), 118 | 109 (11.2), 45 |
     | HG002 sub | 151,157 | 584 (38.6), 286 | 371 (24.5), 195 |
     | HG002 nocl | 473,505 | 1,908 (40.3), 923 | 1,228 (25.9), 621 |

     On the patient arms what is left is reads that fit no state well (best below −23: long
     indels' sites, and 2 panel SNV rows whose 43 such reads the full programme leaves depth
     only, 6 of them ALT under the band) and reads short on a losing state only: among reads
     that fit well (best −12 or better) the margin's call moves for none on the panel, across
     and WES, and for 6 on the gate subset. On HG002 it is repeats the exact-tract rule does
     not see (units varying in length or base, as a CCTTC/CCCTTC run at 13:113419038, where a
     REF state scores −58 in its band and −4 in full) and copies of the read away from X.
   - **Worked rows (HG002, public).** tr 1:193564034 CTATA>C, a 4-base deletion inside a 69-base
     tract of a 34-base unit (the floor there 0): one read scored ALT −20 and REF −10 in its
     band and was called REF; it now scores ALT 0, as in the full programme, and is called ALT;
     truth ALT. c246 13:60977860 G>A, in a 171-base compound TG/TA repeat: one read scored REF
     0, ALT and OTHER −4 (REF); one unit along it fits ALT and OTHER perfectly (0/0/0), as in
     the full programme: depth only; truth undecided.
   - **Cost** (Mac CPU, develop's judge path against this change, timed together): panel +12%,
     across +4%, WES −1%, normals 0, gate subset +2%; HG002 +16 to +38%, the compound-het arms
     (where no read changes) +29 to +38%. The cost is reads in long repeats going on to the
     wider band or the full programme, not the reads whose calls move. The full programme for
     every read costs ×1.35 on the panel, ×2.75 on the gate subset and up to ×2.9 on HG002.
   - **Other tools.** GATK's PairHMM scores every read against every haplotype over the full
     matrix, so it has no band to stand on; htslib's BAQ uses a fixed band of 10 (realn.c). The
     repeat-aware floor is the judge's own rule for when its band may stand in for the full
     programme.
7. **Decided (operator, 2026-10-29): one memory budget per thread, 128 MB.**
   - **Mechanism.** A row's caches (the verdicts of reads alike, the programmes of reads alike
     in bases and placement with their path steps, and the call each read got) grow with its
     reads, and deep rows hold many. The judge's work space (the graph programmes' score
     tables and the aligner's values) grows only to the band a read needs and is kept on the
     thread. The two now share one budget per thread (`THREAD_BUDGET`, 128 MB): the work space
     takes what the reads need, and the caches of the rows alive on the thread the rest. Once
     they reach it, a row adds no more entries and its later reads are judged afresh (the same
     calls; the time a cache saves forgone). A read that alone needs more work space than the
     budget still gets it. A row's caches are freed at its end, and so is work space past a
     quarter of the budget. The bytes counted are those held: the work space's arrays and
     the band each workspace keeps for its path (a band states share, once); each cache's hash
     table whole (its buckets are allocated as it grows); and a kept programme as later reads
     alike grow it. `GBCMS_JUDGE_BUDGET_MB` sets the budget for measurement only
     (`none`: no limit); `GBCMS_JUDGE_CACHE_TRACE=1` traces each row's cache at its end.
   - **Count-neutral.** A test judges a deletion's reads with no limit, with nothing cached and
     with budgets that close the caches part way: every call is the same. On every arm the
     build at its budget and at 0 (nothing cached) gives byte-identical MAFs and read calls to
     the judge without the budget: on HPC (8 HG002 arms, 1,070 patient runs) and on the Mac
     (every arm), first on the judge branch and again on develop a2ca57c1 (HPC runs hzS at the
     budget and hz0 with nothing cached, against develop's judge path c62tS; the default path,
     run hzR, against develop's c62uR).
   - **The knee** (HPC, Linux, whole nodes, Xeon Gold 6348 and Xeon 6747P; CPU against no limit,
     runs compared on one node model; the Mac is the check):

     | Data type | CPU at 32 / 64 / 128 / 256 / 512 MB (no limit 1.000) | Peak RSS, no limit → 128 MB | Knee (CPU within 2%) |
     |:--|:--|:--|:--|
     | IMPACT tumour, IMPACT normal, germline normals, WES, ACCESS XS1 and XS2, gate subset | within 2% at every budget (the shortest ACCESS runs, under 2 s of CPU, vary by up to 2.4% at any budget) | 85 to 250 MB, unchanged | 32 MB |
     | HG002 c246 | 1.020 / 1.013 / 1.008 / 1.005 / 1.002 | 586 → 382 MB | 64 MB |
     | HG002 tr | 1.026 / 1.016 / 1.003 / 0.992 / 0.998 | 476 → 345 MB | 64 MB |
     | HG002 sub | 1.036 / 1.036 / 1.015 / 1.003 / 0.995 | 1,219 → 739 MB | 128 MB |
     | HG002 nocl | 1.041 / 1.027 / 1.014 / 1.005 / 1.000 | 2,808 → 1,645 MB | 128 MB |
     | HG002 delins and compound hets | within 3% at every budget | 123 to 134 MB, unchanged | 32 to 128 MB |

     128 MB is the smallest budget within 2% of no limit's CPU on every data type. The Mac
     agrees: patient types flat, HG002 nocl's knee 128 MB (CPU 1.011; peak 2,653 → 2,067 MB),
     sub 64 MB, the others 32 MB. Measured again on develop a2ca57c1, with path-align's path
     steps in the caches (HPC; 64, 128 and 256 MB and none): 128 MB stays within 2% of no
     limit's CPU on every data type (HG002 nocl 1.009, sub 1.013, tr 1.012, c246 1.002; the
     patient types within 1%, the shortest ACCESS runs, under 1 s of CPU, within 2.1%), and
     cuts HG002 peaks from 2,805 to 1,788 MB (nocl), 1,121 to 748 (sub), 545 to 397 (c246) and
     460 to 329 MB (tr). With the bytes counted as held (the bands, the tables whole, a kept
     programme's growth; Greptile's findings on #379) and on develop dc215a60 (HPC, HG002, runs
     haK128 against haKnone): nocl 1.018 (2,833 to 1,657 MB), sub 1.009 (1,111 to 733), c246
     1.008 (513 to 390), tr 1.005 (458 to 340), delins 1.009; the compound-het arms, about 4 s
     of CPU, vary by up to 5% as at every budget before. Calls are byte-identical at 128 MB and
     with nothing cached on every arm (runs haS and ha0).
   - **Out of memory on HPC, from the harness.** Two gate-subset runs at one thread went over the
     4 GB their task asked for. gbcms itself peaks at about 250 MB there; the harness's runner
     held the traced output in its Python parent. The HPC kit's runner now streams it.
   - **Cache hit rates** (Mac, traced; verdict cache / programme cache): ACCESS duplex 42 to 53%
     / 41 to 66%; IMPACT 2 to 9% / 32 to 41%; WES 2% / 26%; HG002 about 0% / 0.6 to 2.7%. The
     budget closed rows' caches only on HG002 (nocl 33 of 13,367 rows, sub 24 of 5,630).
     On develop a2ca57c1 the caches also hold the path steps: 500 to 670 KB a row on HG002 and
     1,000 to 8,500 KB on the gate subset, the hit rates unchanged; the budget closes rows on
     HG002 nocl (29) and sub (21) only.
   - **Effects.** Speed and memory only: no output changes. RNA mode refuses the sequence judge,
     and the default rules path does not run it. Nextflow's per-task memory is unchanged in code.

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
