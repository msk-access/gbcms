"""Each entry point is told the read judgment it serves, and the sequence judge never counts a
row prepared without the rows a read is weighed against (#304).

Mechanism. Preparation records each row's rows within a read's reach (`Variant.reach_rows`,
from the whole input) for the sequence judge alone. Told nothing, it prepares for the default
rules, which never read them (invariant 3: the intent is plumbed, nothing is computed to be
discarded). `observe_molecules` takes a read judgment (its alignment config's, "rules" unless
set) and tells both preparation and counting. Counting under the sequence judge refuses a row
prepared for the rules: judged without the rows in its reach, its reads would be weighed
against too few candidate sequences, silently. A row built by hand, never prepared, lists no
other row and is counted as before.
"""

import random

import pytest
from helpers import SUITE_READ_JUDGE, count_bam_checked, make_read, write_contig

import gbcms
import gbcms.observations
from gbcms import _rs
from gbcms.models.core import AlignmentConfig, VariantType
from gbcms.models.core import Variant as PyVariant

READ = 120


def _contig():
    rng = random.Random(304)
    return "".join(rng.choice("ACGT") for _ in range(800))


def _alt(contig, p):
    return "A" if contig[p] != "A" else "C"


def _rows(contig):
    """An SNV and a deletion 30 bases on, whose REF reaches it: each is in the other's reach."""
    return [
        _rs.Variant("1", 300, contig[300], _alt(contig, 300), "SNP"),
        _rs.Variant("1", 329, contig[329:340], contig[329], "DEL"),
    ]


def _fixture(tmp_path):
    """Six reads over both rows; every other one carries the SNV's ALT."""
    contig = _contig()
    reads = []
    for i in range(6):
        s = 240 + 10 * i
        seq = list(contig[s : s + READ])
        if i % 2:
            seq[300 - s] = _alt(contig, 300)
        reads.append(make_read(f"r{i}", "".join(seq), s, ((0, READ),), flag=16 if i % 4 else 0))
    fa, bam = write_contig(tmp_path, contig, reads)
    return contig, fa, bam


def _prepare(contig, fa, **kw):
    return _rs.prepare_variants(_rows(contig), fa, 5, False, 1, True, **kw)


def _count(entry, bam, fa, variants, decomposed=None, **kw):
    decomposed = decomposed if decomposed is not None else [None] * len(variants)
    args = (bam, variants, decomposed, 20, 20, True, True, True, False, False, False, 1)
    return entry(*args, reference_fasta=fa, **kw)


def _checked(counts):
    for c in counts:
        assert c.dp >= c.rd + c.ad
        assert c.dpf >= c.rdf + c.adf
        assert c.rd == c.rd_fwd + c.rd_rev
        assert c.ad == c.ad_fwd + c.ad_rev
    return counts


def test_preparation_told_nothing_prepares_for_the_rules(tmp_path):
    """The engine's own default, not the judge-switch suite's: told no read judgment,
    preparation records no row's rows in reach, as told "rules"."""
    prepare = getattr(_rs.prepare_variants, "__wrapped__", _rs.prepare_variants)
    contig, fa, _ = _fixture(tmp_path)
    untold = prepare(_rows(contig), fa, 5, False, 1, True)
    assert [p.gbcms_status for p in untold] == ["PASS", "PASS"]
    assert [p.variant.reach_rows for p in untold] == [[], []]


@pytest.mark.parametrize(
    "entry", [_rs.count_bam_binned, _rs.count_bam_binned_observations], ids=lambda e: e.__name__
)
def test_the_sequence_judge_refuses_rows_prepared_for_the_rules(tmp_path, entry):
    contig, fa, bam = _fixture(tmp_path)
    prepared = [p.variant for p in _prepare(contig, fa, read_judge="rules")]
    with pytest.raises(ValueError, match=r"prepared for the rules.*read_judge \"sequence\""):
        _count(entry, bam, fa, prepared, read_judge="sequence")


def test_the_sequence_judge_refuses_a_twin_prepared_for_the_rules(tmp_path):
    """The decomposed twin is weighed against its parent's rows: a twin prepared for the rules
    beside rows prepared for the judge is refused too."""
    contig, fa, bam = _fixture(tmp_path)
    judged = [p.variant for p in _prepare(contig, fa, read_judge="sequence")]
    ruled = [p.variant for p in _prepare(contig, fa, read_judge="rules")]
    with pytest.raises(ValueError, match="prepared for the rules"):
        _count(_rs.count_bam_binned, bam, fa, judged, [None, ruled[1]], read_judge="sequence")


def test_rows_count_under_the_judge_they_were_prepared_for(tmp_path):
    """Prepared for the judge, rows count under it; prepared for either, they count the same
    under the rules (which read no row's rows in reach), in any bin geometry."""
    contig, fa, bam = _fixture(tmp_path)
    judged = [p.variant for p in _prepare(contig, fa, read_judge="sequence")]
    ruled = [p.variant for p in _prepare(contig, fa, read_judge="rules")]
    args = (bam, judged, [None, None], 20, 20, True, True, True, False, False, False, 1)
    under_judge = _checked(count_bam_checked(*args, reference_fasta=fa, read_judge="sequence"))
    assert under_judge[0].ad == 3
    rules = [
        _checked(count_bam_checked(bam, v, *args[2:], reference_fasta=fa, read_judge="rules"))
        for v in (judged, ruled)
    ]
    assert [(c.dp, c.rd, c.ad) for c in rules[0]] == [(c.dp, c.rd, c.ad) for c in rules[1]]
    assert rules[0][0].ad == 3


def test_a_row_built_by_hand_is_counted_under_the_judge(tmp_path):
    """Never prepared, a row lists no other row; it is counted, as before, not refused."""
    contig, fa, bam = _fixture(tmp_path)
    hand = _rows(contig)[:1]
    (snv,) = _checked(_count(_rs.count_bam_binned, bam, fa, hand, read_judge="sequence"))
    assert snv.ad == 3


def _py_rows(contig):
    return [
        PyVariant(
            chrom="1", pos=300, ref=contig[300], alt=_alt(contig, 300), variant_type=VariantType.SNP
        ),
        PyVariant(
            chrom="1",
            pos=329,
            ref=contig[329:340],
            alt=contig[329],
            variant_type=VariantType.DELETION,
        ),
    ]


@pytest.mark.parametrize(
    "given, told",
    [
        ({}, SUITE_READ_JUDGE),  # the alignment config's default: "rules" unless the suite's
        ({"read_judge": "rules"}, "rules"),
        ({"read_judge": "sequence"}, "sequence"),
        ({"alignment": AlignmentConfig(read_judge="sequence")}, "sequence"),
        ({"alignment": AlignmentConfig(read_judge="sequence"), "read_judge": "rules"}, "rules"),
    ],
    ids=["default", "rules", "sequence", "alignment", "argument-over-alignment"],
)
def test_observe_molecules_tells_preparation_and_counting_its_read_judge(
    tmp_path, monkeypatch, given, told
):
    contig, fa, bam = _fixture(tmp_path)
    seen = {}

    def spy(name):
        real = getattr(gbcms.observations, name)

        def call(*args, **kwargs):
            seen[name] = kwargs.get("read_judge")
            return real(*args, **kwargs)

        monkeypatch.setattr(gbcms.observations, name, call)

    spy("prepare_variants")
    spy("count_bam_binned_observations")
    result = gbcms.observe_molecules(bam, _py_rows(contig), reference_fasta=fa, **given)
    assert seen == {"prepare_variants": told, "count_bam_binned_observations": told}
    assert result.variant_status == ["PASS", "PASS"]
    alt = [o for o in result.observations if o.variant_index == 0 and o.allele == 1]
    assert len(alt) == 3


def test_observe_molecules_refuses_an_unknown_read_judge(tmp_path):
    contig, fa, bam = _fixture(tmp_path)
    with pytest.raises(ValueError, match="read_judge"):
        gbcms.observe_molecules(bam, _py_rows(contig), reference_fasta=fa, read_judge="sequnce")
