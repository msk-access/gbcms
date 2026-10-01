"""Binning invariance: bin geometry is performance only (#170).

The engine groups variants into bins, fetches each bin's reads once, and counts
every member from that cache. Which variants share a bin, and how far a bin
reaches, must never change a count: under any bin window or per-bin cap, one
variant per bin, one variant per call, shuffled input on several threads, or
extra variants that move where bins start, every field of every row is the
same. (The Benjamini-Hochberg q-values are computed over the rows in one call,
so they are excepted when the set of rows changes.)

Window 1 with a cap of 1 is the old per-variant fetch, reproduced through the
production loop. Every read the engine counts overlaps the event's first base,
which every bin's fetch holds, except a DNA read admitted by its soft-clipped
bases: a carrier aligned only after a long event. Those make a bin that stopped
short of its anchor's own span (fixed in 6.4) visible here under a small window;
the Rust property test on `build_genomic_bins` pins the geometry itself.

Fixtures: dense synthetic DNA clusters (SNVs, an MNP, insertions, deletions, a
delins, a homopolymer and an STR, a 55bp deletion with variants inside and past
its breakpoint, multi-allelic siblings, contig start and end) as read pairs,
counted plain and with BAQ, a UMI tag and mFSD on; an RNA locus with a GTF,
spliced reads, strandedness and BAQ; and the repository's real test BAM with
its 155 VCF variants.
"""

import random
from pathlib import Path

import pysam
import pytest
from helpers import ROW_SET_FIELDS, assert_same_counts, make_read, write_contig
from rna_fixtures import E1, E2, SENSE, mk_ref, spliced, through, write_bam, write_fasta, write_gtf

from gbcms import _rs as gbcms_rs

ARGS = {
    "min_mapq": 20,
    "min_baseq": 20,
    "filter_duplicates": True,
    "filter_secondary": True,
    "filter_supplementary": True,
    "filter_qc_failed": False,
    "filter_improper_pair": False,
    "filter_indel": False,
}

# Each geometry against production. Window 1 / cap 1 is the per-variant fetch.
GEOMETRIES = [
    {"bin_window": 1, "bin_max_variants": 1},
    {"bin_window": 7},
    {"bin_window": 13},
    {"bin_window": 20},
    {"bin_max_variants": 2},
    {"bin_max_variants": 3},
    {"bin_window": 50, "bin_max_variants": 2},
]

READ, FRAG = 100, 260
TESTDATA = Path(__file__).parent / "testdata"


# ── Synthetic DNA ─────────────────────────────────────────────────────────────
def _dna_contig():
    rng = random.Random(170)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    parts = [
        rand(900),
        "G" + "A" * 9 + "T",  # homopolymer at 900..909
        rand(90),
        "G" + "CA" * 8 + "T",  # STR at 1000..1017
        rand(1382),
    ]
    return "".join(parts)


DNA = _dna_contig()


def _v(pos, ref_len, alt):
    return (pos, DNA[pos : pos + ref_len], alt)


def _alt_base(pos):
    return "T" if DNA[pos] != "T" else "G"


# (pos, ref, alt), 0-based. Two ALTs at 700 are multi-allelic siblings.
DNA_VARIANTS = [
    _v(3, 1, _alt_base(3)),  # contig start
    _v(620, 1, _alt_base(620)),
    _v(640, 1, _alt_base(640)),
    _v(641, 1, _alt_base(641)),
    _v(660, 2, _alt_base(660) + _alt_base(661)),  # MNP
    (700, DNA[700], "T" if DNA[700] != "T" else "C"),
    (700, DNA[700], "G" if DNA[700] != "G" else "A"),
    _v(720, 1, DNA[720] + "TTG"),  # unique insertion
    _v(760, 5, DNA[760]),  # unique deletion
    _v(900, 1, "GA"),  # homopolymer insertion
    _v(1000, 3, "G"),  # STR deletion (-CA)
    _v(1300, 56, DNA[1300]),  # 55bp deletion
    _v(1320, 1, _alt_base(1320)),  # inside the 55bp deletion
    _v(1360, 1, _alt_base(1360)),  # just past its right breakpoint
    _v(1400, 2, "TTT" if DNA[1400:1402] != "TT" else "GGG"),  # delins
    _v(len(DNA) - 5, 1, _alt_base(len(DNA) - 5)),  # contig end
]


def _alt_read(contig, start, variant, length):
    """A read of `length` bases from `start` carrying `variant`, written as the
    aligner would: M for substitutions, M D/I M after the shared prefix."""
    pos, ref, alt = variant
    hap = contig[:pos] + alt + contig[pos + len(ref) :]
    seq = hap[start : start + length]
    k = 0
    while k < min(len(ref), len(alt)) and ref[k] == alt[k]:
        k += 1
    if len(ref) == len(alt):
        return seq, ((0, len(seq)),)
    left = pos + k - start
    dele, ins = len(ref) - k, len(alt) - k
    ops = [(0, left)]
    if dele:
        ops.append((2, dele))
    if ins:
        ops.append((1, ins))
    rest = len(seq) - left - ins
    if rest > 0:
        ops.append((0, rest))
    return seq, tuple(ops)


def _pairs(contig, variants, seed, step=3, umi=False):
    """Read pairs from every `step`-th start: each fragment is REF or carries one
    variant fully inside both mates' reach; a few low-quality bases."""
    rng = random.Random(seed)
    reads = []
    for n, s in enumerate(range(0, len(contig) - FRAG, step)):
        inside = [v for v in variants if s < v[0] and v[0] + len(v[1]) + 6 < s + READ]
        variant = rng.choice(inside) if inside and rng.random() < 0.45 else None
        frag = []
        for mate, start in ((1, s), (2, s + FRAG - READ)):
            if variant is not None and mate == 1:
                seq, cigar = _alt_read(contig, start, variant, READ)
            else:
                seq, cigar = contig[start : start + READ], ((0, READ),)
            quals = [30 if rng.random() > 0.03 else 10 for _ in seq]
            flag = 0x1 | 0x2 | (0x40 | 0x20 if mate == 1 else 0x80 | 0x10)
            r = make_read(f"f{n}", seq, start, cigar, flag=flag, quals=quals)
            r.next_reference_id = 0
            r.next_reference_start = s + FRAG - READ if mate == 1 else s
            r.template_length = FRAG if mate == 1 else -FRAG
            if umi:
                r.set_tag("RX", f"U{n % 7}")
            frag.append(r)
        reads.extend(frag)
    return reads


def _prepared(fa, rows):
    pvs = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", p, r, a, "X") for p, r, a in rows], fa, 5, False, 1, True
    )
    assert all(pv.gbcms_status == "PASS" for pv in pvs), [pv.gbcms_status_reason for pv in pvs]
    variants = [pv.variant for pv in pvs]
    decomposed = [pv.decomposed_variant for pv in pvs]
    # Co-annotated rows at one position are each other's siblings, as the pipeline groups them.
    siblings = [
        [w for j, w in enumerate(variants) if j != i and w.pos == v.pos]
        for i, v in enumerate(variants)
    ]
    return variants, decomposed, siblings


# ── The invariance checks ─────────────────────────────────────────────────────
def _check_invariance(bam, variants, decomposed, siblings, per_row=True, **kw):
    def run(vs=variants, ds=decomposed, ss=siblings, threads=1, **geometry):
        return gbcms_rs.count_bam_binned(
            bam, vs, ds, **ARGS, threads=threads, sibling_variants=ss, **kw, **geometry
        )

    base = run()
    for c in base:
        assert c.dp >= c.rd + c.ad and c.dpf >= c.rdf + c.adf
        assert c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
    assert sum(c.ad for c in base) > 0, "the fixture must exercise ALT reads"

    for geometry in GEOMETRIES:
        assert_same_counts(run(**geometry), base, f"geometry {geometry}")

    # Shuffled input on four threads: results come back in input order.
    order = list(range(len(variants)))
    random.Random(7).shuffle(order)
    shuffled = run(
        [variants[i] for i in order],
        [decomposed[i] for i in order],
        [siblings[i] for i in order],
        threads=4,
    )
    back = [None] * len(order)
    for k, i in enumerate(order):
        back[i] = shuffled[k]
    assert_same_counts(back, base, "shuffled, 4 threads")

    if per_row:
        one_by_one = [
            run([v], [d], [s])[0] for v, d, s in zip(variants, decomposed, siblings, strict=True)
        ]
        assert_same_counts(one_by_one, base, "one call per row", skip=ROW_SET_FIELDS)
    return base, run


@pytest.mark.parametrize(
    "features",
    [{}, {"apply_baq": True, "umi_tag": "RX", "mfsd": True}],
    ids=["plain", "baq-umi-mfsd"],
)
def test_dna_counts_do_not_depend_on_bin_geometry(tmp_path, features):
    reads = _pairs(DNA, DNA_VARIANTS, seed=1, umi="umi_tag" in features)
    fa, bam = write_contig(tmp_path, DNA, reads, "dna")
    variants, decomposed, siblings = _prepared(fa, DNA_VARIANTS)
    base, run = _check_invariance(bam, variants, decomposed, siblings, **features)

    # Decoy rows interleaved before, among and after the targets move where bins
    # start; the targets' counts stay.
    decoy_rows = [
        _v(p, 1, _alt_base(p))
        for p in range(480, 1500, 37)
        if p not in {r[0] for r in DNA_VARIANTS}
    ]
    decoys, _, _ = _prepared(fa, decoy_rows)
    with_decoys = run(
        decoys[: len(decoys) // 2] + variants + decoys[len(decoys) // 2 :],
        [None] * (len(decoys) // 2) + decomposed + [None] * (len(decoys) - len(decoys) // 2),
        [[]] * (len(decoys) // 2) + siblings + [[]] * (len(decoys) - len(decoys) // 2),
    )
    targets = with_decoys[len(decoys) // 2 : len(decoys) // 2 + len(variants)]
    assert_same_counts(targets, base, "with decoys", skip=ROW_SET_FIELDS)


def test_observation_rows_do_not_depend_on_bin_geometry(tmp_path):
    """The per-molecule rows come from the same loop: same rows, same calls."""
    reads = _pairs(DNA, DNA_VARIANTS, seed=2)
    fa, bam = write_contig(tmp_path, DNA, reads, "obs")
    variants, decomposed, siblings = _prepared(fa, DNA_VARIANTS)

    def rows(**geometry):
        counts, obs = gbcms_rs.count_bam_binned_observations(
            bam, variants, decomposed, **ARGS, threads=1, sibling_variants=siblings, **geometry
        )
        return counts, sorted(
            (o.variant_index, o.molecule_hash, o.allele, o.best_qual, o.min_mapq) for o in obs
        )

    base_counts, base_rows = rows()
    assert base_rows
    for geometry in GEOMETRIES[:3]:
        counts, got = rows(**geometry)
        assert_same_counts(counts, base_counts, f"observations, geometry {geometry}")
        assert got == base_rows, f"observation rows differ under {geometry}"


def test_clip_carriers_past_a_long_anchor_do_not_depend_on_bin_geometry(tmp_path):
    """A 60bp REF to a 10bp ALT anchoring its bin, carriers aligned only after the
    event with the event clipped at their start: counted reads that never overlap
    the anchor, so a bin whose fetch stopped short of the anchor's own span (the
    6.4 bug the Rust property test pins) would lose them under a small window."""
    ref, alt, hap = _clip_case()
    reads = [
        _paired(make_read(f"r{i}", ref[s : s + READ], s, ((0, READ),)), FRAG)
        for i, s in enumerate(range(CLIP_POS - 60, CLIP_POS - 50))
    ]
    after = CLIP_POS + 60  # first REF base after the event; on the haplotype, CLIP_POS + 10
    for i, hs in enumerate(range(CLIP_POS - 40, CLIP_POS - 20)):
        clip = CLIP_POS + 10 - hs
        read = make_read(f"l{i}", hap[hs : hs + READ], after, ((4, clip), (0, READ - clip)))
        reads.append(_paired(read, FRAG, reverse=True))
    fa, bam = write_contig(tmp_path, ref, reads, "clip")
    variants, decomposed, siblings = _prepared(fa, [(CLIP_POS, ref[CLIP_POS : CLIP_POS + 60], alt)])
    base, _ = _check_invariance(bam, variants, decomposed, siblings)
    assert (base[0].rd, base[0].ad) == (10, 20), "the clipped carriers must count"


CLIP_POS = 400


def _clip_case():
    rng = random.Random(60)
    ref = "".join(rng.choice("ACGT") for _ in range(900))
    alt = ""
    while not alt or alt[0] == ref[CLIP_POS] or alt[-1] == ref[CLIP_POS + 59]:
        alt = "".join(rng.choice("ACGT") for _ in range(10))
    return ref, alt, ref[:CLIP_POS] + alt + ref[CLIP_POS + 60 :]


def _paired(read, tlen, reverse=False):
    """One mate of a proper pair spanning `tlen` bases."""
    read.flag = 0x1 | 0x2 | 0x40 | (0x10 if reverse else 0x20)
    read.next_reference_id = 0
    read.next_reference_start = read.reference_start + (0 if reverse else max(0, tlen - READ))
    read.template_length = -tlen if reverse else tlen
    return read


# ── RNA ───────────────────────────────────────────────────────────────────────
def test_rna_counts_do_not_depend_on_bin_geometry(tmp_path):
    ref = mk_ref()
    donor, acceptor = E1[1], E2[0]
    snvs = [E1[1] - 30, E1[1] - 12, E1[1] - 3, E2[0] + 4, E2[0] + 25]
    reads = spliced(ref, donor, acceptor, 30, "s")
    for i, p in enumerate(snvs):
        alt = "T" if ref[p] != "T" else "G"
        reads += through(ref, p, alt, 6, f"a{i}_") + through(ref, p, ref[p], 6, f"r{i}_")
    ins_at = E1[1] - 20
    for i in range(6):
        s = ins_at - 50 + i
        seq = ref[s : ins_at + 1] + "AC" + ref[ins_at + 1 : s + 98]
        reads.append(
            make_read(
                f"i{i}",
                seq,
                s,
                ((0, ins_at + 1 - s), (1, 2), (0, len(seq) - (ins_at + 1 - s) - 2)),
                flag=SENSE,
            )
        )
    fa = write_fasta(tmp_path, ref, contig="1")
    bam = write_bam(tmp_path, ref, reads)
    gtf = write_gtf(tmp_path, contig="1")
    rows = [(p, ref[p], "T" if ref[p] != "T" else "G") for p in snvs] + [
        (ins_at, ref[ins_at], ref[ins_at] + "AC")
    ]
    variants, decomposed, siblings = _prepared(str(fa), rows)
    _check_invariance(
        str(bam),
        variants,
        decomposed,
        siblings,
        mode="rna",
        apply_baq=True,
        enforce_strandedness=True,
        gtf_path=str(gtf),
    )


# ── The repository's real test BAM ────────────────────────────────────────────
EXTRACT_START = 11_180_000  # the bundled FASTA holds chr1 from here (0-based), 20kb


def test_real_bam_counts_do_not_depend_on_bin_geometry():
    """41,556 real alignments: the VCF's three variants plus an SNV every 97bp of
    the 20kb region (REF from the bundled extract; unprepared, since the extract's
    coordinates are not the BAM's), split into many bins by a 1kb window."""
    variants = []
    for line in (TESTDATA / "integration_test_variants.vcf").read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        variants.append(gbcms_rs.Variant(f[0], int(f[1]) - 1, f[3], f[4], "X"))
    extract = pysam.FastaFile(str(TESTDATA / "integration_test_reference.fa")).fetch("1").upper()
    for i in range(300, len(extract) - 300, 97):
        alt = "T" if extract[i] != "T" else "G"
        variants.append(gbcms_rs.Variant("1", EXTRACT_START + i, extract[i], alt, "SNP"))
    assert len(variants) > 150
    bam = str(TESTDATA / "sample1_integration_test.bam")
    none, empty = [None] * len(variants), [[] for _ in variants]

    def run(**geometry):
        return gbcms_rs.count_bam_binned(
            bam, variants, none, **ARGS, threads=2, sibling_variants=empty, **geometry
        )

    base = run()
    # Real data at these loci is reference (the VCF's ALTs are below 0.2% and
    # filtered here): depth, REF and strand counts carry the check.
    assert sum(c.dp for c in base) > 10_000 and sum(c.rd for c in base) > 10_000
    for geometry in (
        {"bin_window": 1000},
        {"bin_window": 1, "bin_max_variants": 1},
        {"bin_max_variants": 3},
    ):
        assert_same_counts(run(**geometry), base, f"real BAM, geometry {geometry}")


def test_a_bin_geometry_below_one_is_rejected(tmp_path):
    fa, bam = write_contig(tmp_path, DNA[:400], [], "empty")
    (v,) = [
        pv.variant
        for pv in gbcms_rs.prepare_variants(
            [gbcms_rs.Variant("1", 100, DNA[100], "T" if DNA[100] != "T" else "G", "X")],
            fa,
            5,
            False,
            1,
            True,
        )
    ]
    for bad in ({"bin_window": 0}, {"bin_max_variants": 0}):
        with pytest.raises(ValueError, match="at least 1"):
            gbcms_rs.count_bam_binned(bam, [v], [None], **ARGS, threads=1, **bad)
