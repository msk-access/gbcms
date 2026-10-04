"""Group 3 (RNA) contracts.

C15 #173 (RNA part): an RNA read's soft clip is not allele evidence. STAR cannot
splice an overhang shorter than its minimum (3 nt annotated, 5 nt novel) or one
with mismatches; it soft-clips it, and the clip holds the next exon's bases
(GATK's RNA workflow runs HaplotypeCaller with ``-dont-use-soft-clipped-bases``;
phASER, WASP and ASEReadCounter never read clips).

C32 #213: a spliced read is judged against the haplotypes spliced at its own
junction. Cutting its windows at the exon edge drops the bases past the junction,
which tell a length-changing allele from REF there.

R4 #185: the gene strand at an intronic position comes from the transcripts
spanning it (first exon start to last exon end), with the exons' agreement rule:
every stranded transcript agrees, else no strand. Splice-site variants (donor
+1/+2, acceptor -1/-2) are intronic. Where both strands' exons or transcripts
cover a position, both genes' transcripts carry the allele: no strand, every read
counts (REDItools' annotation mode leaves such sites undetermined).

O8: ``OBSERVED_ALLELE`` / ``COEXISTING_ALLELE`` name what the reads carry, with
n/m read counts set beside ``alt_count`` and ``ref_count``. They must read the
reads those counts read: under strandedness enforcement (the RNA default) no
antisense read, and in RNA the mapping rule the counts use (a read below
``--min-mapq`` with ``NH:i:1`` is a unique mapper and counts). Allele tallies in
other tools are taken over the same filtered read set as their primary counts
(bam-readcount, ASEReadCounter ``otherBases``, VarDict).
"""

import pytest
from helpers import make_read
from rna_fixtures import (
    E1,
    E2,
    READ_LEN,
    SENSE,
    mk_ref,
    run_rna,
    write_bam,
    write_fasta,
    write_gtf,
    write_vcf,
)

ANTISENSE = 0  # a forward single-end read is antisense to a '+' gene under dUTP
POS = E1[0] + 100  # 0-based SNV position, mid-exon


def _carriers(ref, base, n, prefix, flag=SENSE, mapq=60, tags=()):
    out = []
    for i in range(n):
        s = POS - 50 + (i % 7)
        seq = ref[s:POS] + base + ref[POS + 1 : s + READ_LEN]
        r = make_read(f"{prefix}{i}", seq, s, ((0, READ_LEN),), flag=flag, mapq=mapq)
        for k, v in tags:
            r.set_tag(k, v)
        out.append(r)
    return out


def _alleles(ref):
    others = [b for b in "ACGT" if b != ref[POS]]
    return others[0], others[1]  # the given ALT, another allele


def _run(tmp_path, ref, reads, mode_extra=()):
    fa, gtf = write_fasta(tmp_path, ref), write_gtf(tmp_path)
    bam = write_bam(tmp_path, ref, reads)
    alt, _ = _alleles(ref)
    vcf = write_vcf(tmp_path, [(POS + 1, ref[POS], alt)])
    (row,) = run_rna(tmp_path, vcf, bam, fa, gtf, extra=mode_extra)
    return row


def _flags(row, name):
    return [f for f in row["gbcms_diagnostic"].split(";") if f.startswith(name)]


def test_observed_alleles_read_no_antisense_read(tmp_path):
    """15 antisense reads carry another allele, 5 sense reads the given ALT and
    10 sense reads REF. Under enforcement the counts exclude the antisense reads,
    so nothing the counted reads carry outnumbers the given ALT: no flag, the same
    as the BAM without the antisense reads."""
    ref = mk_ref()
    alt, other = _alleles(ref)
    sense = _carriers(ref, alt, 5, "a") + _carriers(ref, ref[POS], 10, "r")
    row = _run(tmp_path, ref, sense + _carriers(ref, other, 15, "x", flag=ANTISENSE))
    assert (int(row["alt_count"]), int(row["ref_count"])) == (5, 10)
    (tmp_path / "sense_only").mkdir()
    alone = _run(tmp_path / "sense_only", ref, sense)
    assert _flags(row, "COEXISTING_ALLELE") == _flags(alone, "COEXISTING_ALLELE") == []
    assert _flags(row, "OBSERVED_ALLELE") == _flags(alone, "OBSERVED_ALLELE") == []


def test_observed_alleles_read_the_unique_mappers_the_counts_read(tmp_path):
    """15 sense reads at MAPQ 0 with NH:i:1 (unique mappers the RNA counts keep)
    carry another allele; 5 sense reads carry the given ALT. The counts read all
    20, so the diagnostic names the other allele over the same 20: 15/5."""
    ref = mk_ref()
    alt, other = _alleles(ref)
    reads = _carriers(ref, alt, 5, "a") + _carriers(ref, other, 15, "u", mapq=0, tags=(("NH", 1),))
    row = _run(tmp_path, ref, reads)
    assert int(row["total_count"]) == 20
    (flag,) = _flags(row, "COEXISTING_ALLELE")
    assert flag.endswith(f":{ref[POS]}>{other}:15/5)")


def test_dna_observed_alleles_read_every_counted_read(tmp_path):
    """Guard: DNA has no strand filter, so the same 15 forward reads carrying
    another allele are counted and named, 15/5."""
    from glob import glob

    from helpers import read_maf_output
    from rna_fixtures import runner

    from gbcms.cli import app

    ref = mk_ref()
    alt, other = _alleles(ref)
    reads = _carriers(ref, alt, 5, "a") + _carriers(ref, other, 15, "x", flag=ANTISENSE)
    fa, bam = write_fasta(tmp_path, ref), write_bam(tmp_path, ref, reads)
    vcf = write_vcf(tmp_path, [(POS + 1, ref[POS], alt)])
    out = tmp_path / "dna"
    res = runner.invoke(
        app,
        ["dna", "-v", str(vcf), "-b", f"S:{bam}", "-f", str(fa), "-o", str(out), "--format", "maf"],
    )
    assert res.exit_code == 0, res.output
    (row,) = list(read_maf_output(glob(str(out / "*.maf"))[0]))
    (flag,) = _flags(row, "COEXISTING_ALLELE")
    assert flag.endswith(f":{ref[POS]}>{other}:15/5)")


# ── C15 (RNA): an exon-edge clip holding the next exon's bases ──────────────
EDGE = E1[1]  # 300: E1's exclusive end, the donor


def _edge_case():
    """An MNP over E1's bases 297-298 (E1 ends at 300). Reads align over 297
    with the ALT's first base there and soft-clip the rest; the clip holds the
    next exon's bases, set here to read as the ALT's second base and its flank."""
    ref = list(mk_ref())
    p0 = EDGE - 3
    alt = next(b for b in "TACG" if b != ref[p0]) + next(b for b in "GTCA" if b != ref[p0 + 1])
    ref[E2[0]] = alt[1]
    ref[E2[0] + 1 : E2[0] + 3] = ref[p0 + 2 : p0 + 4]
    return "".join(ref), p0, alt


def test_an_rna_soft_clip_is_not_allele_evidence(tmp_path):
    """Reads aligned over the MNP's first base (the ALT's) whose next 6+ bases,
    the next exon's, are soft-clipped. Read into the clip they fit the ALT; the
    clip is not evidence in RNA, so they are depth only, not ALT."""
    ref, p0, alt = _edge_case()
    reads = []
    for i in range(6):
        clip = 6 + i
        s = p0 + 1 - (READ_LEN - clip)
        seq = ref[s:p0] + alt[0] + ref[E2[0] : E2[0] + clip]
        reads.append(make_read(f"c{i}", seq, s, ((0, p0 + 1 - s), (4, clip)), flag=SENSE))
    fa, gtf, bam = write_fasta(tmp_path, ref), write_gtf(tmp_path), write_bam(tmp_path, ref, reads)
    vcf = write_vcf(tmp_path, [(p0 + 1, ref[p0 : p0 + 2], alt)])
    (row,) = run_rna(tmp_path, vcf, bam, fa, gtf)
    assert int(row["total_count"]) == 6
    assert int(row["alt_count"]) == 0


# ── C32: a spliced read is read across its own junction ─────────────────────
def _junction_ref():
    """E1 ends at 300 and E2 starts at 500. The next exon's first base equals
    E1's last base (as at the FORTE case), and its second differs from its first."""
    ref = list(mk_ref())
    ref[E2[0]] = ref[EDGE - 1]
    if ref[E2[0] + 1] == ref[E2[0]]:
        ref[E2[0] + 1] = "A" if ref[E2[0]] != "A" else "C"
    if ref[EDGE - 1] == ref[EDGE - 2]:
        raise AssertionError("fixture: E1's last two bases must differ")
    return "".join(ref)


def _junction_case():
    ref = _junction_ref()
    p0 = EDGE - 3
    err = next(b for b in "TACG" if b != ref[p0])
    alt = err + ref[p0 + 1] + ref[EDGE - 1]  # 2 bases -> 3: the exon one base longer
    return ref, p0, err, alt


def _spliced_run(tmp_path, ref, reads, row):
    fa, gtf, bam = write_fasta(tmp_path, ref), write_gtf(tmp_path), write_bam(tmp_path, ref, reads)
    (out,) = run_rna(tmp_path, write_vcf(tmp_path, [row]), bam, fa, gtf)
    return out


def test_a_spliced_ref_read_with_one_error_is_not_alt(tmp_path):
    """REF reads spliced at 300 -> 500 with one error (the ALT's first base) at
    297. Cut at the exon edge, their bases fit the ALT; read across the junction
    the ALT's extra base would push E1's last base before the next exon, and the
    reads show the next exon instead: one mismatch against either allele, not ALT."""
    ref, p0, err, alt = _junction_case()
    reads = []
    for i in range(6):
        s = p0 - 40 - i
        left = EDGE - s
        seq = ref[s:p0] + err + ref[p0 + 1 : EDGE] + ref[E2[0] : E2[0] + READ_LEN - left]
        reads.append(
            make_read(
                f"e{i}", seq, s, ((0, left), (3, E2[0] - EDGE), (0, READ_LEN - left)), flag=SENSE
            )
        )
    out = _spliced_run(tmp_path, ref, reads, (p0 + 1, ref[p0 : p0 + 2], alt))
    assert int(out["total_count"]) == 6
    assert int(out["alt_count"]) == 0


def test_a_spliced_alt_carrier_is_alt(tmp_path):
    """Guard: carriers of the same ALT (E1 one base longer: err, ref[298],
    ref[299], ref[299]) spliced at the donor count ALT, read across the junction."""
    ref, p0, err, alt = _junction_case()
    reads = []
    for i in range(6):
        s = p0 - 40 - i
        m1 = EDGE - 1 - s  # ref[s:299] aligned (the error a mismatch)
        seq = ref[s:p0] + err + ref[p0 + 1 : EDGE - 1] + ref[EDGE - 1] + ref[EDGE - 1]
        right = READ_LEN - len(seq)
        seq += ref[E2[0] : E2[0] + right]
        cigar = ((0, m1), (1, 1), (0, 1), (3, E2[0] - EDGE), (0, right))
        reads.append(make_read(f"a{i}", seq, s, cigar, flag=SENSE))
    out = _spliced_run(tmp_path, ref, reads, (p0 + 1, ref[p0 : p0 + 2], alt))
    assert int(out["alt_count"]) == 6


# ── R4: the gene strand at intronic positions ───────────────────────────────
INTRON1 = (E1[1], E2[0])  # [300, 500)


def _strand_reads(ref, pos, n_sense, n_anti):
    out = []
    for i in range(n_sense + n_anti):
        s = pos - 50 + (i % 9)
        flag = SENSE if i < n_sense else ANTISENSE
        out.append(make_read(f"s{i}", ref[s : s + READ_LEN], s, ((0, READ_LEN),), flag=flag))
    return out


def _strand_run(tmp_path, ref, pos, reads, extra_gtf=""):
    fa, bam = write_fasta(tmp_path, ref), write_bam(tmp_path, ref, reads)
    gtf = write_gtf(tmp_path)
    if extra_gtf:
        gtf.write_text(gtf.read_text() + extra_gtf)
    alt = next(b for b in "ACGT" if b != ref[pos])
    (row,) = run_rna(tmp_path, write_vcf(tmp_path, [(pos + 1, ref[pos], alt)]), bam, fa, gtf)
    return row


def _depths(row):
    return int(row["ref_count"]), int(row["rna_sense_depth"]), int(row["rna_antisense_depth"])


@pytest.mark.parametrize(
    "pos",
    [
        pytest.param(INTRON1[0], id="donor+1"),
        pytest.param(INTRON1[0] + 1, id="donor+2"),
        pytest.param(INTRON1[1] - 2, id="acceptor-2"),
        pytest.param(INTRON1[1] - 1, id="acceptor-1"),
        pytest.param(INTRON1[0] + 100, id="deep-intronic"),
    ],
)
def test_intronic_positions_take_the_spanning_transcripts_strand(tmp_path, pos):
    """10 sense and 6 antisense unspliced reads over an intronic position of a
    '+' gene: under enforcement the antisense reads count nowhere but in
    rna_antisense_depth, as at an exonic position."""
    ref = mk_ref()
    row = _strand_run(tmp_path, ref, pos, _strand_reads(ref, pos, 10, 6))
    assert _depths(row) == (10, 10, 6)


def test_exonic_positions_are_unchanged(tmp_path):
    """Guard: the exonic case the rule extends."""
    ref = mk_ref()
    pos = E1[0] + 100
    assert _depths(_strand_run(tmp_path, ref, pos, _strand_reads(ref, pos, 10, 6))) == (10, 10, 6)


def test_opposite_strand_exons_leave_no_strand(tmp_path):
    """Guard: a '-' gene's exon over the '+' gene's exonic position: no strand,
    every read counts (both genes' transcripts carry the allele)."""
    ref = mk_ref()
    pos = E1[0] + 100
    anti = f'chr1\tTEST\texon\t{pos - 20}\t{pos + 30}\t.\t-\t.\tgene_id "G2"; transcript_id "T2";\n'
    row = _strand_run(tmp_path, ref, pos, _strand_reads(ref, pos, 10, 6), anti)
    assert _depths(row) == (16, 16, 0)


def test_both_strands_spanning_an_intronic_position_leave_no_strand(tmp_path):
    """Guard: a '-' transcript spanning the same intronic position (its exons
    on either side): no strand, every read counts."""
    ref = mk_ref()
    pos = INTRON1[0] + 100
    anti = "".join(
        f'chr1\tTEST\texon\t{a}\t{b}\t.\t-\t.\tgene_id "G2"; transcript_id "T2";\n'
        for a, b in ((pos - 40, pos - 30), (pos + 30, pos + 40))
    )
    row = _strand_run(tmp_path, ref, pos, _strand_reads(ref, pos, 10, 6), anti)
    assert _depths(row) == (16, 16, 0)


@pytest.mark.parametrize("pos", [INTRON1[0], INTRON1[1] - 1], ids=["donor+1", "acceptor-1"])
def test_an_intronic_position_lies_in_no_transcripts_exon(tmp_path, pos):
    """The exon index is end-inclusive (coitrees) but was built from half-open
    exons and queried over two bases, so the first and last intron bases read as
    exonic: per-transcript counts listed T1 there. No exon holds them."""
    ref = mk_ref()
    row = _strand_run(tmp_path, ref, pos, _strand_reads(ref, pos, 10, 6))
    assert row["transcript_read_counts"] == ""


@pytest.mark.parametrize("pos", [INTRON1[0] - 1, INTRON1[1]], ids=["exon1-last", "exon2-first"])
def test_an_exons_edge_bases_lie_in_its_transcript(tmp_path, pos):
    """Guard: E1's last base and E2's first base are T1's."""
    ref = mk_ref()
    row = _strand_run(tmp_path, ref, pos, _strand_reads(ref, pos, 10, 6))
    assert row["transcript_read_counts"].startswith("T1:")


# ── R5: a read spliced inside an indel's change interval ────────────────────
RUN = (EDGE - 4, EDGE + 4)  # an A-run over E1's last 4 bases and the intron's first 4


def _run_ref():
    ref = list(mk_ref())
    ref[RUN[0] : RUN[1]] = "A" * (RUN[1] - RUN[0])
    for p in (RUN[0] - 1, RUN[1]):
        if ref[p] == "A":
            ref[p] = "C"
    return "".join(ref)


def _run_reads(ref, n_spliced, n_through, carrier=False):
    reads = []
    for i in range(n_spliced):
        s = EDGE - 40 - i
        left = EDGE - s
        ins = "A" if carrier else ""
        right = READ_LEN - left - len(ins)
        seq = ref[s:EDGE] + ins + ref[E2[0] : E2[0] + right]
        cigar = ((0, left),) + (((1, 1),) if carrier else ()) + ((3, E2[0] - EDGE), (0, right))
        reads.append(make_read(f"sp{i}", seq, s, cigar, flag=SENSE))
    for i in range(n_through):
        s = RUN[0] - 40 - i
        reads.append(make_read(f"th{i}", ref[s : s + READ_LEN], s, ((0, READ_LEN),), flag=SENSE))
    return reads


def _run_row(tmp_path, ref, reads):
    fa, gtf, bam = write_fasta(tmp_path, ref), write_gtf(tmp_path), write_bam(tmp_path, ref, reads)
    anchor = RUN[0] - 1
    vcf = write_vcf(tmp_path, [(anchor + 1, ref[anchor], ref[anchor] + "A")])
    (row,) = run_rna(tmp_path, vcf, bam, fa, gtf)
    return int(row["ref_count"]), int(row["alt_count"]), int(row["total_count"])


def test_a_read_spliced_inside_the_run_is_not_ref(tmp_path):
    """A +A insertion in an A-run that crosses E1's end: reads spliced at 300
    show 4 of the run's 8 A's and none past the splice (the insertion may sit in
    the intron's part). They are depth only; reads through the run are REF."""
    ref = _run_ref()
    assert _run_row(tmp_path, ref, _run_reads(ref, 5, 5)) == (5, 0, 10)


def test_a_carrier_spliced_inside_the_run_is_not_alt(tmp_path):
    """Guard: the same reads with I(A) at the junction. Their bases fit both
    alleles (an A more in the exon's part, or the intron's part spliced out), so
    they are depth only, as the unspliced read ending at the same base is."""
    ref = _run_ref()
    assert _run_row(tmp_path, ref, _run_reads(ref, 5, 5, carrier=True)) == (5, 0, 10)


def test_a_splice_past_the_window_leaves_the_call(tmp_path):
    """Guard: an insertion in a run ending 10 bases before E1's end; reads
    spliced at 300 hold the whole window in their block: REF."""
    ref = list(mk_ref())
    lo, hi = EDGE - 18, EDGE - 10
    ref[lo:hi] = "A" * (hi - lo)
    for p in (lo - 1, hi):
        if ref[p] == "A":
            ref[p] = "C"
    ref = "".join(ref)
    reads = []
    for i in range(5):
        s = EDGE - 60 - i
        left = EDGE - s
        seq = ref[s:EDGE] + ref[E2[0] : E2[0] + READ_LEN - left]
        reads.append(
            make_read(
                f"g{i}", seq, s, ((0, left), (3, E2[0] - EDGE), (0, READ_LEN - left)), flag=SENSE
            )
        )
    fa, gtf, bam = write_fasta(tmp_path, ref), write_gtf(tmp_path), write_bam(tmp_path, ref, reads)
    vcf = write_vcf(tmp_path, [(lo, ref[lo - 1], ref[lo - 1] + "A")])
    (row,) = run_rna(tmp_path, vcf, bam, fa, gtf)
    assert (int(row["ref_count"]), int(row["total_count"])) == (5, 5)


def test_a_spliced_ref_read_needs_only_the_spliced_flank(tmp_path):
    """Guard: E1 ends ...GG with the intron's G (the planted GT) continuing the
    run, so genomic windows grow into the intron; the next exon starts with C and
    ends the run. REF reads spliced at 300 holding 3 of the next exon's bases show
    the spliced windows whole: REF."""
    ref = list(mk_ref())
    p0 = EDGE - 3
    ref[p0 + 1] = ref[p0 + 2] = "G"  # REF = ref[297] G, E1's last base G, then the intron's G
    if ref[p0] == "G":
        ref[p0] = "T"
    ref[E2[0]] = "C"
    ref = "".join(ref)
    alt = next(b for b in "CAT" if b != ref[p0]) + "AT"
    reads = []
    for i in range(6):
        left = READ_LEN - 3
        s = EDGE - left
        seq = ref[s:EDGE] + ref[E2[0] : E2[0] + 3]
        reads.append(make_read(f"r{i}", seq, s, ((0, left), (3, E2[0] - EDGE), (0, 3)), flag=SENSE))
    out = _spliced_run(tmp_path, ref, reads, (p0 + 1, ref[p0 : p0 + 2], alt))
    assert (int(out["ref_count"]), int(out["total_count"])) == (6, 6)


# ── Review cases (2026-10-04): junction chains, left splices, gap forms, clips ──
N_LEN = E2[0] - EDGE  # intron 1


def _other(ref, avoid):
    return next(b for b in "ACGT" if b not in avoid)


def _spliced_ref_reads(ref, n, start_back=40):
    reads = []
    for i in range(n):
        s = EDGE - start_back - i
        left = EDGE - s
        seq = ref[s:EDGE] + ref[E2[0] : E2[0] + READ_LEN - left]
        cig = ((0, left), (3, N_LEN), (0, READ_LEN - left))
        reads.append(make_read(f"r{i}", seq, s, cig, flag=SENSE))
    return reads


def _counts(tmp_path, ref, reads, row):
    out = _spliced_run(tmp_path, ref, reads, row)
    return int(out["ref_count"]), int(out["alt_count"])


@pytest.mark.parametrize("exon_len", [3, 20])
def test_a_spliced_read_is_read_through_a_short_exon(tmp_path, exon_len):
    """REF reads spliced 300 -> 500, through a short exon [500, 500 + len), then
    -> 800, at a delins whose windows reach past the short exon: their bases are
    the transcript, REF (every junction the windows reach is followed)."""
    ref = mk_ref()
    p0 = EDGE - 3
    alt = "".join(_other(ref, (ref[p0], ref[p0 + 1], ref[p0 - 1], ref[p0 + 2])) for _ in range(9))
    reads = []
    for i in range(6):
        s = EDGE - 40 - i
        left = EDGE - s
        rest = READ_LEN - left - exon_len
        seq = ref[s:EDGE] + ref[E2[0] : E2[0] + exon_len] + ref[800 : 800 + rest]
        cig = ((0, left), (3, N_LEN), (0, exon_len), (3, 800 - E2[0] - exon_len), (0, rest))
        reads.append(make_read(f"r{i}", seq, s, cig, flag=SENSE))
    assert _counts(tmp_path, ref, reads, (p0 + 1, ref[p0 : p0 + 2], alt))[0] == 6


@pytest.mark.parametrize("alt_len", [1, 3, 5])
def test_spliced_windows_before_the_event(tmp_path, alt_len):
    """A delins at E2's second and third bases; reads spliced 300 -> 500 start
    in E1: REF reads are REF and carriers ALT, read across the junction on the
    left."""
    ref = mk_ref()
    p0 = E2[0] + 1
    alt = "".join(
        _other(ref, (ref[p0], ref[p0 + 1], ref[p0 - 1], ref[p0 + 2])) for _ in range(alt_len)
    )
    reads = _spliced_ref_reads(ref, 6, start_back=30)
    d = alt_len - 2
    for i in range(4):
        s = EDGE - 30 - i
        left = EDGE - s
        right = READ_LEN - left
        body = ref[E2[0] : p0] + alt + ref[p0 + 2 : p0 + 2 + READ_LEN]
        seq = ref[s:EDGE] + body[:right]
        if d > 0:
            cig = ((0, left), (3, N_LEN), (0, 3), (1, d), (0, right - 3 - d))
        else:
            cig = ((0, left), (3, N_LEN), (0, 2), (2, -d), (0, right - 2))
        reads.append(make_read(f"a{i}", seq, s, cig, flag=SENSE))
    assert _counts(tmp_path, ref, reads, (p0 + 1, ref[p0 : p0 + 2], alt)) == (6, 4)


@pytest.mark.parametrize(
    "form",
    [
        "D+N",
        pytest.param(
            "shifted-N",
            marks=pytest.mark.xfail(
                strict=True, reason="gap form: a junction starting inside the event"
            ),
        ),
    ],
)
def test_a_delins_carrier_counts_however_the_gap_is_written(tmp_path, form):
    """E1's last four bases replaced by one base X. The carriers' bases are
    ...X then the next exon, written as X 3D N or as X with the junction starting
    inside the event (a donor three bases up); REF reads are spliced at the
    annotated donor. The call follows the bases, not the gap's placement: ALT."""
    ref = mk_ref()
    p0 = EDGE - 4
    x = _other(ref, (ref[p0], ref[E2[0]], ref[EDGE - 1]))
    reads = _spliced_ref_reads(ref, 6)
    for i in range(6):
        s = EDGE - 40 - i
        m1 = p0 + 1 - s
        seq = ref[s:p0] + x + ref[E2[0] : E2[0] + READ_LEN - m1]
        if form == "D+N":
            cig = ((0, m1), (2, 3), (3, N_LEN), (0, READ_LEN - m1))
        else:
            cig = ((0, m1), (3, N_LEN + 3), (0, READ_LEN - m1))
        reads.append(make_read(f"a{i}", seq, s, cig, flag=SENSE))
    assert _counts(tmp_path, ref, reads, (p0 + 1, ref[p0:EDGE], x)) == (6, 6)


@pytest.mark.xfail(strict=True, reason="C15: a mid-exon clip is the read's own bases")
def test_a_mid_exon_clip_is_allele_evidence(tmp_path):
    """An MNP 150 bases from either exon edge: REF reads align in full; ALT reads
    have the MNP's second base and the rest soft-clipped (a local aligner clips a
    mismatching end). No exon edge or junction is within the clip's reach, so it
    holds the same exon's bases, the read's own: ALT."""
    ref = mk_ref()
    p0 = E1[0] + 50
    alt = _other(ref, (ref[p0],)) + _other(ref, (ref[p0 + 1],))
    end = p0 + 5
    s = end - READ_LEN
    reads = []
    for i in range(6):
        reads.append(make_read(f"r{i}", ref[s:end], s, ((0, READ_LEN),), flag=SENSE))
        seq = ref[s:p0] + alt + ref[p0 + 2 : end]
        reads.append(make_read(f"a{i}", seq, s, ((0, p0 + 1 - s), (4, end - p0 - 1)), flag=SENSE))
    assert _counts(tmp_path, ref, reads, (p0 + 1, ref[p0 : p0 + 2], alt)) == (6, 6)
