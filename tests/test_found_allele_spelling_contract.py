"""The sequence judge's found alleles are written past the guard's masks (#304, RJ-25).

The guard finds the other alleles the reads spell at the locus (RJ-25: on 3 or more
fragments and 5% of those spelling the allele they are nearest), reading each spelling past
its masks: a masked base (a recurring change beside the event, background) is no difference.
The sequence judge takes each found allele as a site state. Written over the guard's window
with its carriers' bases, a found allele would hold a masked base as its own change, and the
judge would read that base literally for every state: a read carrying the masked change then
fits the found allele better than X or REF by that base alone. So the judge writes each
found allele as the bases it changes, past the masks (the allele's own base where a mask
takes a carrier's): the masks stay on every state, and a found allele that so written is a
candidate the site already holds (X with a listed neighbour beside it) is that candidate.

The shapes, each with a germline change one base before the event on every read (masked as
background, inside the guard's window):
- X's carriers all carry a listed neighbour beside the event: the found allele is X with
  that neighbour, which every candidate holds as a context. The carriers count ALT.
- REF's reads and a real other allele of REF's length beside the mask: REF's reads count
  REF, the other allele's partial evidence.
- X's carriers and a real other allele of X's length beside the mask: X's carriers count
  ALT, the other allele's partial evidence.

The patient case behind them (a 17-base tandem duplication, a listed SNV past the copies on
every carrier, recurring changes masked inside the window): 293 carriers per run made
partial by the window spelling; the record is in docs/proposals/sequence-judge.md.
Counts are (ref_count, alt_count, total_count, partial_alt) per row in input order.
"""

import glob
import random

import pysam
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()
READ_LEN = 100
S = 299  # the delins' first REF base, 0-based (VCF position 300)
DEL = 7  # REF bases it replaces
INS = "TTCT"  # the bases it writes
GERM = S - 2  # a germline change one base before the event (the base before S is the anchor's left)
LISTED = S + DEL  # a listed SNV on the first base after the event
INSIDE_REF = S + 3  # a change inside the deleted bases (REF's other allele)
INSIDE_ALT = 1  # a change at the second written base (X's other allele)
COUNT_COLUMNS = (
    "total_count",
    "ref_count",
    "alt_count",
    "total_count_fragment",
    "ref_count_fragment",
    "alt_count_fragment",
    "ref_count_forward",
    "ref_count_reverse",
    "alt_count_forward",
    "alt_count_reverse",
)


def _contig(seed=371):
    rng = random.Random(seed)
    c = [rng.choice("ACGT") for _ in range(1200)]
    # the delins' REF bases share none with what it writes at either end
    c[S : S + DEL] = "GAGAGAG"
    return "".join(c)


REF = _contig()


def _other(base):
    return next(b for b in "ACGT" if b != base)


def _change(seq, at, base=None):
    s = list(seq)
    s[at] = base if base is not None else _other(seq[at])
    return "".join(s)


def _ref_read(name, hap, start, flag=0):
    return make_read(name, hap[start : start + READ_LEN], start, ((0, READ_LEN),), flag=flag)


def _carrier(name, hap, start, ins=INS, flag=0):
    """A carrier of the delins, its gap written at the event: `hap` is REF's coordinates
    with the read's other changes, `ins` the bases it writes."""
    left = S - start
    right = READ_LEN - left - len(ins)
    seq = hap[start:S] + ins + hap[S + DEL : S + DEL + right]
    return make_read(name, seq, start, ((0, left), (1, len(ins)), (2, DEL), (0, right)), flag=flag)


def _starts(n, k):
    # starts spread so every read covers the event and its flank, both strands
    return [(S - 70 + (i * 7 + k) % 41, 16 if i % 2 else 0) for i in range(n)]


def _run(tmp_path, rows, reads, ref=REF):
    fa = tmp_path / "ref.fasta"
    fa.write_text(">1\n" + ref + "\n")
    pysam.faidx(str(fa))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": len(ref)}]}
    raw, bam = tmp_path / "raw.bam", tmp_path / "reads.bam"
    with pysam.AlignmentFile(raw, "wb", header=header) as fh:
        for a in reads:
            fh.write(a)
    pysam.sort("-o", str(bam), str(raw))
    pysam.index(str(bam))
    vcf = tmp_path / "rows.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(f"1\t{p}\t.\t{r}\t{a}\t.\t.\t.\n" for p, r, a in rows)
    )
    out = tmp_path / "out"
    args = ["dna", "-v", str(vcf), "-b", f"S:{bam}", "-f", str(fa), "-o", str(out)]
    result = runner.invoke(
        app, args + ["--format", "maf", "--read-judge", "sequence"], env={"COLUMNS": "400"}
    )
    assert result.exit_code == 0, result.output
    got = list(read_maf_output(glob.glob(str(out / "*.maf"))[0]))
    for r in got:
        n = {k: int(r[k]) for k in COUNT_COLUMNS}
        assert n["total_count"] >= n["ref_count"] + n["alt_count"]
        assert n["total_count_fragment"] >= n["ref_count_fragment"] + n["alt_count_fragment"]
        assert n["ref_count"] == n["ref_count_forward"] + n["ref_count_reverse"]
        assert n["alt_count"] == n["alt_count_forward"] + n["alt_count_reverse"]
    return [
        (int(r["ref_count"]), int(r["alt_count"]), int(r["total_count"]), int(r["partial_alt"]))
        for r in got
    ]


X = (S, REF[S - 1 : S + DEL], REF[S - 1] + INS)  # VCF (1-based): the anchor, then the delins
GERMLINE = _change(REF, GERM)


def test_carriers_with_a_listed_neighbour_count_alt_past_a_masked_change(tmp_path):
    """Every read carries the germline change one base before the event (masked); X's 30
    carriers all carry the listed SNV on the base after it. The reads spell one other allele
    of X's length, X with that SNV (on all 30), with the germline base in its window: written
    past the mask, that allele is X with a listed neighbour, a candidate already, and the
    carriers count ALT. Written with their germline base, it would take every carrier from X
    by that base alone. (The SNV's own row is not this rule's: its carriers carry the delins
    beside it, the same on either writing.)"""
    snv = (LISTED + 1, REF[LISTED], _other(REF[LISTED]))
    with_snv = _change(GERMLINE, LISTED, snv[2])
    reads = [_ref_read(f"r{i}", GERMLINE, st, fl) for i, (st, fl) in enumerate(_starts(60, 0))]
    reads += [_carrier(f"a{i}", with_snv, st, flag=fl) for i, (st, fl) in enumerate(_starts(30, 3))]
    (x_rd, x_ad, x_dp, x_part), _ = _run(tmp_path, [X, snv], reads)
    assert (x_rd, x_ad, x_part) == (60, 30, 0)
    assert x_dp == 90


def test_ref_reads_count_ref_beside_another_allele_past_a_masked_change(tmp_path):
    """Every read carries the germline change (masked); 20 of the 300 REF-length fragments
    also carry a change inside the deleted bases: another allele of REF's length (20 of 280,
    above 5%). Written past the mask, it is that change alone: its carriers are partial
    evidence and the other 280 count REF. Written with their germline base, every REF read
    would tie with it by that base."""
    odd = _change(GERMLINE, INSIDE_REF)
    reads = [_ref_read(f"r{i}", GERMLINE, st, fl) for i, (st, fl) in enumerate(_starts(280, 0))]
    reads += [_ref_read(f"o{i}", odd, st, fl) for i, (st, fl) in enumerate(_starts(20, 5))]
    reads += [
        _carrier(f"a{i}", GERMLINE, st, flag=fl) for i, (st, fl) in enumerate(_starts(100, 3))
    ]
    ((rd, ad, dp, part),) = _run(tmp_path, [X], reads)
    assert (rd, ad, part) == (280, 100, 20)
    assert dp == 400


def test_carriers_count_alt_beside_another_allele_of_their_length_past_a_masked_change(tmp_path):
    """Every read carries the germline change (masked); 20 of X's 100 carriers write a
    different second base: another allele of X's length (20 of 80, above 5%). Written past
    the mask it is that base alone: its carriers are partial evidence and the other 80
    count ALT."""
    written = INS[:INSIDE_ALT] + _other(INS[INSIDE_ALT]) + INS[INSIDE_ALT + 1 :]
    reads = [_ref_read(f"r{i}", GERMLINE, st, fl) for i, (st, fl) in enumerate(_starts(300, 0))]
    reads += [_carrier(f"a{i}", GERMLINE, st, flag=fl) for i, (st, fl) in enumerate(_starts(80, 3))]
    reads += [
        _carrier(f"o{i}", GERMLINE, st, ins=written, flag=fl)
        for i, (st, fl) in enumerate(_starts(20, 6))
    ]
    ((rd, ad, dp, part),) = _run(tmp_path, [X], reads)
    assert (rd, ad, part) == (300, 80, 20)
    assert dp == 400
