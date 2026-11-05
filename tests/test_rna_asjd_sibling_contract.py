"""ASJD reads the same sibling verdict as the main counts (RJ-5).

At an RNA row with a co-annotated sibling, the main counts and the per-transcript counts
give a read that is the sibling's molecule to the sibling: a REF read that is ALT for a
sibling whose change lies where this row's REF is read is partial evidence, not REF; an ALT
match the sibling explains as well or better is partial evidence, not ALT. The junction-usage
test (ASJD) took the classifier's REF and ALT as they were, so it counted those reads in
its REF and ALT partitions: a sibling's molecule voted for this row's junction. It asks the
same two questions now, so every tally holds one verdict per read. Only the `asjd_*` columns
move; the main counts and the per-transcript counts are unchanged (and pinned here).

Committed red before the change; flipped green with it. Geometry as `test_rna_pins_contract`:
the three-exon gene of `rna_fixtures`, reads spliced exon 1 to exon 2, sense under dUTP.
"""

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

DONOR, ACCEPTOR = E1[1], E2[0]


def _spliced(ref, s, seq_left, cigar_left, name, quals=None):
    """A read from `s`: `seq_left` and `cigar_left` cover the bases up to the donor (their
    reference length is DONOR - s), then N over intron 1 and the rest in exon 2."""
    n_query = sum(ln for op, ln in cigar_left if op in (0, 1))
    right = READ_LEN - n_query
    seq = seq_left + ref[ACCEPTOR : ACCEPTOR + right]
    cig = (*cigar_left, (3, ACCEPTOR - DONOR), (0, right))
    q = None
    if quals is not None:
        q = list(quals) + [30] * (len(seq) - len(quals))
    return make_read(name, seq, s, cig, flag=SENSE, quals=q)


def _run(tmp_path, ref, rows, reads):
    return run_rna(
        tmp_path,
        write_vcf(tmp_path, rows),
        write_bam(tmp_path, ref, reads),
        write_fasta(tmp_path, ref),
        write_gtf(tmp_path),
    )


def _tx(row):
    """(AD, RD, DP) of the one transcript."""
    name, counts = row["transcript_read_counts"].split(":")
    assert name == "T1"
    return tuple(int(x) for x in counts.split(","))


# ── the REF guard: a sibling's SNV on the deleted base ──────────────────────────
RUN = 260  # an 8-A run at [260, 268) in exon 1


def _run_ref():
    ref = list(mk_ref())
    ref[RUN - 1] = "C"
    ref[RUN : RUN + 8] = list("A" * 8)
    ref[RUN + 8] = "G"
    return "".join(ref)


def test_a_sibling_snv_on_the_deleted_base_is_not_ref_in_asjd(tmp_path):
    """Row X deletes the run's first A (anchor C); sibling S is A>G on that base. Six reads
    are REF, 4 carry X's deletion, 8 carry S's G (no gap: REF for X by its alignment, but S's
    molecules). Main RD 6 and the transcript's RD 6, 8 reads partial evidence; ASJD's REF
    partition holds the 6, not the 14."""
    ref = _run_ref()
    s_of = lambda i: 230 + (i % 5)  # noqa: E731
    reads = []
    for i in range(6):
        s = s_of(i)
        reads.append(_spliced(ref, s, ref[s:DONOR], ((0, DONOR - s),), f"ref{i}"))
    for i in range(4):  # X's deletion: the run's first A gone
        s = s_of(i)
        seq = ref[s:RUN] + ref[RUN + 1 : DONOR]
        cig = ((0, RUN - s), (2, 1), (0, DONOR - RUN - 1))
        reads.append(_spliced(ref, s, seq, cig, f"del{i}"))
    for i in range(8):  # S's molecules
        s = s_of(i)
        seq = ref[s:RUN] + "G" + ref[RUN + 1 : DONOR]
        reads.append(_spliced(ref, s, seq, ((0, DONOR - s),), f"sib{i}"))
    rows = [(RUN, ref[RUN - 1 : RUN + 1], ref[RUN - 1]), (RUN + 1, "A", "G")]
    x, s_row = _run(tmp_path, ref, rows, reads)
    assert (int(x["ref_count"]), int(x["alt_count"]), int(x["partial_alt"])) == (6, 4, 8)
    assert _tx(x) == (4, 6, 18)
    assert int(x["asjd_n_alt_total"]) == 4
    assert int(x["asjd_n_ref_total"]) == 6, "a sibling's molecule votes in ASJD's REF partition"
    assert int(s_row["asjd_n_alt_total"]) == 8


# ── the ALT guard: a read the sibling explains as well as this row's ALT ─────────
P = 262  # TT at [262, 264), flanked as in the delins contract


def _delins_ref():
    ref = list(mk_ref())
    ref[P - 10 : P + 12] = list("GCAGTCAGGA" + "TT" + "AGCGTCAGGT")
    return "".join(ref)


def test_a_read_the_sibling_explains_as_well_is_not_the_rows_alt_in_asjd(tmp_path):
    """TT>GTC beside a co-annotated C inserted after the second T: the two ALTs differ only
    at the first base. 20 reads carrying TTC with that base below min BQ match both exactly:
    not the delins's ALT. 10 reads carrying GTC stay its ALT. Main AD 10, the transcript's
    AD 10; ASJD's ALT partition holds the 10, not the 30."""
    ref = _delins_ref()
    reads = []
    for i in range(20):  # TTC, the first T's base below min BQ
        s = P - 40 + (i % 20)
        k = P - s
        seq = ref[s : P + 2] + "C" + ref[P + 2 : DONOR]
        cig = ((0, k + 2), (1, 1), (0, DONOR - P - 2))
        q = [30] * len(seq)
        q[k] = 5
        reads.append(_spliced(ref, s, seq, cig, f"ins{i}", quals=q))
    for i in range(10):  # GTC
        s = P - 40 + (i % 10)
        k = P - s
        seq = ref[s:P] + "GTC" + ref[P + 2 : DONOR]
        cig = ((0, k + 2), (1, 1), (0, DONOR - P - 2))
        reads.append(_spliced(ref, s, seq, cig, f"cx{i}"))
    rows = [(P + 1, "TT", "GTC"), (P + 2, "T", "TC")]
    cx, _ = _run(tmp_path, ref, rows, reads)
    assert (int(cx["alt_count"]), int(cx["ref_count"])) == (10, 0)
    assert _tx(cx)[0] == 10
    assert int(cx["asjd_n_alt_total"]) == 10, "a sibling's molecule votes in ASJD's ALT partition"
