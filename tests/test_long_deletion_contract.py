"""Pure deletions of 50 bp or more are judged by their bases at the POS breakpoint
(RJ-24, C40 #254; docs/reference/read-judgment.md).

The engine's calls on every long-deletion read shape of the read-judgment spec
must equal the read census's, which judges the same windows from the bases alone
(clips and the hypothesis placement included), and must not depend on how the
variants are binned (`count_checked`). Two shapes the census does not model are
left to the spec: a change recurring across reads (masked in the engine) and a
read's own 5+ bp indel inside the deleted span (C26, from its alignment).
A recurring change on the carriers alone names the larger allele they carry.
Synthetic and PHI-free.
"""

import glob

import pytest
import read_judgment_cases as rj
from census import assert_matches, census
from helpers import count_checked, read_maf_output, write_contig
from typer.testing import CliRunner

from gbcms import _rs
from gbcms.cli import app

A = rj.A
CENSUS_LIMITS = {
    "carriers with a recurring change two bases past the microhomology",
    "REF reads with a recurring change inside the window",
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


def test_a_recurring_change_on_the_carriers_names_the_larger_allele(tmp_path):
    """Every carrier shows a change two bases past the junction; no REF read can
    (they read the deleted bases there). The change is masked, so the carriers
    count ALT, and the diagnostic names the deletion with the change."""
    n, mh = rj.LONG_DEL["del 60"]
    contig = rj._long_contig(n, mh)
    shape = "carriers with a recurring change two bases past the microhomology"
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
    changed = A + n + mh + rj.RECUR_PAST  # the carriers' change, 0-based
    alt_base = rj._other(contig[changed])
    larger = f"1:{A}:{contig[A - 1 : changed + 1]}>{contig[A - 1]}{contig[A + n : changed]}{alt_base}:4/0"
    # named in canonical form: the shared anchor trimmed
    pos, ref_al, alt_al = (
        A,
        contig[A - 1 : changed + 1],
        contig[A - 1] + contig[A + n : changed] + alt_base,
    )
    while len(ref_al) > 1 and len(alt_al) > 1 and ref_al[0] == alt_al[0]:
        pos, ref_al, alt_al = pos + 1, ref_al[1:], alt_al[1:]
    larger = f"OBSERVED_ALLELE(1:{pos}:{ref_al}>{alt_al}:4/0)"
    assert larger in row["gbcms_diagnostic"].split(";"), row["gbcms_diagnostic"]
