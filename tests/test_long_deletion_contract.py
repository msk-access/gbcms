"""Pure deletions of 50 bp or more are judged by their bases at the POS breakpoint
(RJ-24, C40 #254; docs/reference/read-judgment.md).

The engine's calls on every long-deletion read shape of the read-judgment spec
must equal the read census's, which judges the same windows from the bases alone
(clips and the hypothesis placement included), and must not depend on how the
variants are binned (`count_checked`). Two kinds of shape the census does not model are
left to the spec: changes recurring across reads (masked in the engine) and a
read's own 5+ bp indel inside the deleted span (C26, from its alignment).
A recurring change on the carriers alone names the larger allele they carry.
Synthetic and PHI-free.
"""

import glob

import pytest
import read_judgment_cases as rj
from census import assert_matches, census
from helpers import count_checked, make_read, read_maf_output, write_contig
from typer.testing import CliRunner

from gbcms import _rs
from gbcms.cli import app

A = rj.A
CENSUS_LIMITS = {
    "carriers with a recurring change two bases past the microhomology",
    "carriers with a recurring change in the flank",
    "right-clipped carriers with a recurring change in the flank",
    "carriers with a recurring change on the window's last base",
    "REF reads with a recurring change inside the window",
    "REF reads with a run of changes in the flank",
    "carriers and REF reads with a run of changes in the flank",
    "REF read with its own 20-bp deletion inside the deleted span",
}
SHAPES = [
    (v, shape)
    for v, shapes in rj.LONG_DEL_SHAPES.items()
    for shape in shapes
    if shape not in CENSUS_LIMITS
]


def _prepared(tmp_path, v, shape):
    n, mh = rj.LONG_DEL[v]
    contig = rj._long_contig(n, mh)
    fa, bam = write_contig(tmp_path, contig, rj._long_del_reads(n, mh, shape), "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", A - 1, contig[A - 1 : A + n], contig[A - 1], "DEL")],
        fa,
        5,
        False,
        1,
        True,
    )
    return contig, bam, pv.variant


@pytest.mark.parametrize("v, shape", SHAPES, ids=[f"{v} | {s}" for v, s in SHAPES])
def test_the_engine_matches_the_census(tmp_path, v, shape):
    contig, bam, variant = _prepared(tmp_path, v, shape)
    (c,) = count_checked(bam, [variant])
    assert_matches(c, census(bam, contig, variant), f"{v} | {shape}")


def _canonical(pos: int, ref: str, alt: str, contig: str) -> tuple[int, str, str]:
    """observed::canonical: trim shared last bases (extending left when an allele
    empties), then shared leading bases keeping one; 0-based position."""
    while True:
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        elif not ref or not alt:
            pos, ref, alt = pos - 1, contig[pos - 1] + ref, contig[pos - 1] + alt
        else:
            break
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        pos, ref, alt = pos + 1, ref[1:], alt[1:]
    return pos, ref, alt


@pytest.mark.parametrize(
    "shape, at",
    [
        ("carriers with a recurring change two bases past the microhomology", "past"),
        ("right-clipped carriers with a recurring change in the flank", "flank"),
    ],
)
def test_a_recurring_change_on_the_carriers_names_the_larger_allele(tmp_path, shape, at):
    """Every carrier shows a change no REF read can (two bases past the junction,
    where REF reads read deleted bases; or in the flank, on carriers alone). The
    change is masked, so the carriers count ALT, and the diagnostic names the
    deletion with the change: from the windows' first base through the change, and
    at least through the junction."""
    n, mh = rj.LONG_DEL["del 60"]
    contig = rj._long_contig(n, mh)
    reads = rj._long_del_reads(n, mh, shape) + rj._long_del_reads(n, mh, "REF read, long flanks")
    for i, r in enumerate(reads):
        r.query_name = f"m{i}"
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        f"##fileformat=VCFv4.2\n##contig=<ID=1,length={len(contig)}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"1\t{A}\t.\t{contig[A - 1 : A + n]}\t{contig[A - 1]}\t.\t.\t.\n"
    )
    out = tmp_path / "out"
    res = CliRunner().invoke(
        app, ["dna", "-v", str(vcf), "-b", f"S:{bam}", "-f", fa, "-o", str(out), "--format", "maf"]
    )
    assert res.exit_code == 0, res.output
    (row,) = read_maf_output(glob.glob(str(out / "*.maf"))[0])
    assert (int(row["ref_count"]), int(row["alt_count"])) == (4, 4)
    first = A - 4  # the windows' first base, 0-based
    if at == "past":
        changed = A + n + mh + rj.RECUR_PAST
        ref_al = contig[first : changed + 1]
        alt_al = contig[first:A] + contig[A + n : changed] + rj._other(contig[changed])
    else:
        changed = A - 2  # through the junction: the first base past it
        ref_al = contig[first : A + n + 1]
        alt_al = (
            contig[first:changed]
            + rj._other(contig[changed])
            + contig[changed + 1 : A]
            + contig[A + n]
        )
    pos, ref_al, alt_al = _canonical(first, ref_al, alt_al, contig)
    assert len(ref_al) - len(alt_al) == n
    larger = f"OBSERVED_ALLELE(1:{pos + 1}:{ref_al}>{alt_al}:4/0)"
    assert larger in row["gbcms_diagnostic"].split(";"), row["gbcms_diagnostic"]


def test_a_germline_snp_at_the_anchor_keeps_clipped_carriers(tmp_path):
    """Every read, REF and ALT, shows a germline SNP at the base before the deletion,
    where the flank's other three bases equal the deleted bases before the right flank:
    the guard masks the SNP, and a carrier clipped at the junction (placed from the
    right flank) still shows the junction by that base against the deleted tail.
    D-written and clipped carriers count alike."""
    n, mh = rj.LONG_DEL["del 60"]
    c = list(rj._long_contig(n, mh))
    c[A - 4 : A - 1] = c[A + n - 4 : A + n - 1]
    contig = "".join(c)
    assert contig[A - 1] != contig[A + n - 1]  # still placed leftmost

    def snp(seq, s):
        return seq[: A - 1 - s] + "A" + seq[A - s :]  # the contig has no A

    hap = contig[:A] + contig[A + n :]
    reads = []
    for i in range(4):
        s = A - 50 - i
        reads.append(
            make_read(
                f"d{i}", snp(hap[s : s + 120], s), s, ((0, A - s), (2, n), (0, 120 - (A - s)))
            )
        )
        reads.append(make_read(f"r{i}", snp(contig[s : s + 120], s), s, ((0, 120),)))
        k = 20 + i  # bases before the junction, in the clip
        seq = snp(hap[A - k : A - k + 120], A - k)
        reads.append(make_read(f"c{i}", seq, A + n, ((4, k), (0, 120 - k))))
    for r in reads:
        r.flag = 0x1 | 0x2 | 0x20 | 0x40
        r.next_reference_id, r.next_reference_start, r.template_length = (
            0,
            r.reference_start + 250,
            400,
        )
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", A - 1, contig[A - 1 : A + n], contig[A - 1], "DEL")],
        fa,
        5,
        False,
        1,
        True,
    )
    (got,) = count_checked(bam, [pv.variant])
    assert (got.rd, got.ad) == (4, 8)


@pytest.mark.parametrize("min_baseq", [0, 1])
def test_an_n_clip_shows_no_flank_at_any_minimum_quality(tmp_path, min_baseq):
    """REF reads starting inside the deleted span behind an N clip show nothing of the
    left flank: at --min-baseq 0 or 1 one base's evidence weighs nothing, and a tie
    places no read."""
    n, mh = rj.LONG_DEL["del 60"]
    contig = rj._long_contig(n, mh)
    reads = rj._long_del_reads(n, mh, "REF read starting inside the deletion behind an N clip")
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", A - 1, contig[A - 1 : A + n], contig[A - 1], "DEL")],
        fa,
        5,
        False,
        1,
        True,
    )
    (got,) = count_checked(bam, [pv.variant], min_baseq=min_baseq)
    assert got.ad == 0
