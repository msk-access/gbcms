"""`INDEL_SHARED`: RJ-28's indel twin under the sequence judge (#304, operator 2026-10-30).

A listed indel or delins that holds a row's whole change (the row with listed neighbours
beside it: a delins that is two listed deletions together) counts its carriers as the
row's ALT, so one molecule counts at each row that holds it. As `MNV_SHARED` does for an
SNV beside an MNV holding its change, the row names, per such listed allele, how many of
its ALT reads spell it whole: `INDEL_SHARED(chrom:pos:REF>ALT:n/m)` in
`gbcms_diagnostic`, the overlap shown without changing a count or adding a column.

Shapes:
- the worked row (a patient panel row, allele shapes only): a 14-base deletion X, a listed
  7-base deletion Y four bases before it, and a listed 24>3 delins written from Y's
  second base that is the two together. Every carrier carries both deletions: X and Y
  count them ALT and each names the delins; the delins row names nothing.
- a listed deletion that holds X's bases and two more (no listed row makes up the rest):
  another allele at X's site, its carriers partial at X, and X names nothing.
Counts per row: (ref_count, alt_count, partial_alt).
"""

import glob
import random

import pysam
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()
READ_LEN = 100


def _contig(seed, length=1500):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _deleting(name, ref, start, cuts):
    """A read of `ref` from `start` whose molecule lacks each `[a, b)` of `cuts` (0-based,
    ascending), its gaps written there, READ_LEN bases long."""
    seq, cigar, at, left = [], [], start, READ_LEN
    for a, b in cuts:
        n = a - at
        seq.append(ref[at:a])
        cigar += [(0, n), (2, b - a)]
        left -= n
        at = b
    seq.append(ref[at : at + left])
    cigar.append((0, left))
    return make_read(name, "".join(seq), start, tuple(cigar))


def _run(tmp_path, ref, rows, reads):
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
    args = [
        "dna",
        "-v",
        str(vcf),
        "-b",
        f"S:{bam}",
        "-f",
        str(fa),
        "-o",
        str(out),
        "--format",
        "maf",
    ]
    result = runner.invoke(app, args + ["--read-judge", "sequence"], env={"COLUMNS": "400"})
    assert result.exit_code == 0, result.output
    got = list(read_maf_output(glob.glob(str(out / "*.maf"))[0]))
    for r in got:
        dp, rd, ad = int(r["total_count"]), int(r["ref_count"]), int(r["alt_count"])
        assert dp >= rd + ad
        assert rd == int(r["ref_count_forward"]) + int(r["ref_count_reverse"])
        assert ad == int(r["alt_count_forward"]) + int(r["alt_count_reverse"])
    return got


def _counts(r):
    return int(r["ref_count"]), int(r["alt_count"]), int(r["partial_alt"])


def _tokens(r, name):
    return [t for t in r["gbcms_diagnostic"].split(";") if t.startswith(name)]


def _worked(seed=373):
    """The worked shape on a random contig: Y deletes [y0, y0+7), X deletes [x0, x0+14),
    four bases kept between; L is the delins from y0 through X's last base writing the
    kept bases. The delins shares no first or last base with what it writes."""
    rng = random.Random(seed)
    while True:
        ref = _contig(rng.randrange(1 << 30))
        y0 = 700
        y1, x0 = y0 + 7, y0 + 11
        x1 = x0 + 14
        if ref[y0] != ref[y1] and ref[x1 - 1] != ref[x0 - 1]:
            break
    x = (x0, ref[x0 - 1 : x1], ref[x0 - 1])  # VCF: the anchor at 0-based x0-1
    y = (y0, ref[y0 - 1 : y1], ref[y0 - 1])
    delins = (y0 + 1, ref[y0:x1], ref[y1:x0])
    return ref, (y0, y1, x0, x1), x, y, delins


def test_a_delins_that_is_two_listed_deletions_is_named_at_each(tmp_path):
    """20 carriers of both deletions (spanning the delins whole) and 30 REF reads: X and
    Y each count the 20 ALT and name the delins, 20 of 20; the delins counts them too
    and names nothing."""
    ref, (y0, y1, x0, x1), x, y, delins = _worked()
    reads = [_deleting(f"c{i}", ref, y0 - 60 + i, [(y0, y1), (x0, x1)]) for i in range(20)]
    reads += [
        make_read(f"r{i}", ref[y0 - 50 + i : y0 + 50 + i], y0 - 50 + i, ((0, READ_LEN),))
        for i in range(30)
    ]
    got = {
        (r["Reference_Allele"], r["Tumor_Seq_Allele2"]): r
        for r in _run(tmp_path, ref, [x, y, delins], reads)
    }
    rx = next(
        r
        for k, r in got.items()
        if k[0].endswith(x[1][1:]) and len(k[0]) - len(k[1].strip("-")) == 14
    )
    ry = next(
        r
        for k, r in got.items()
        if k[0].endswith(y[1][1:]) and len(k[0]) - len(k[1].strip("-")) == 7
    )
    rl = next(r for k, r in got.items() if k[0] == delins[1])
    want = f"INDEL_SHARED(1:{delins[0]}:{delins[1]}>{delins[2]}:20/20)"
    assert _counts(rx)[1:] == (20, 0)
    assert _counts(ry)[1:] == (20, 0)
    assert _counts(rl)[1] == 20
    assert _tokens(rx, "INDEL_SHARED") == [want]
    assert _tokens(ry, "INDEL_SHARED") == [want]
    assert _tokens(rl, "INDEL_SHARED") == []


def test_a_listed_deletion_holding_more_than_the_rows_change_is_not_named(tmp_path):
    """X deletes three bases; a listed deletion takes them and the two after (no listed row
    makes up those two): another allele at X's site. 15 X carriers, 10 carriers of the
    longer deletion, 30 REF reads: X counts its 15 ALT and the 10 partial, and names
    nothing."""
    ref = _contig(374)
    x0 = 700
    x = (x0, ref[x0 - 1 : x0 + 3], ref[x0 - 1])
    longer = (x0, ref[x0 - 1 : x0 + 5], ref[x0 - 1])
    reads = [_deleting(f"x{i}", ref, x0 - 50 + i, [(x0, x0 + 3)]) for i in range(15)]
    reads += [_deleting(f"l{i}", ref, x0 - 45 + i, [(x0, x0 + 5)]) for i in range(10)]
    reads += [
        make_read(f"r{i}", ref[x0 - 50 + i : x0 + 50 + i], x0 - 50 + i, ((0, READ_LEN),))
        for i in range(30)
    ]
    got = {
        (r["Reference_Allele"], r["Tumor_Seq_Allele2"]): r
        for r in _run(tmp_path, ref, [x, longer], reads)
    }
    rx = next(r for k, r in got.items() if len(k[0]) - len(k[1].strip("-")) == 3)
    assert _counts(rx)[1:] == (15, 10)
    assert _tokens(rx, "INDEL_SHARED") == []
