"""Per-transcript count columns list their transcripts in ID order.

The two columns (``transcript_read_counts``, ``transcript_fragment_counts``) read the
overlapping transcripts off an interval tree whose visiting order depends on the CPU's
SIMD backend and on how equal intervals fall, so the same run on two machines wrote the
same entries in another order. The order is now a property of the transcript IDs alone:
ascending, in both columns, whatever order the GTF lists the transcripts in.
"""

import pytest
from rna_fixtures import (
    E1,
    E2,
    mk_ref,
    run_rna,
    through,
    write_bam,
    write_fasta,
    write_vcf,
)

POS = E1[0] + 100  # 0-based SNV position, mid-exon of every transcript below

# Each transcript has two exons, so it has an intron and gets its own entry. Their first
# exons start and end in no alphabetical order, and ENST_B's first exon has the very same
# interval as ENST_A's (a tie the tree can break either way).
EXONS = {
    "ENST_C": ((150, 260), E2),
    "ENST_A": ((100, 300), E2),
    "ENST_E": ((120, 1000), (1100, 1150)),
    "ENST_B": ((100, 300), E2),
    "ENST_D": ((190, 210), E2),
}
SORTED_IDS = sorted(EXONS)


def _gtf(tmp_path, transcripts, name):
    """A GTF of the given transcript IDs, in the order given."""
    path = tmp_path / name
    path.write_text(
        "".join(
            f'chr1\tTEST\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "G1"; transcript_id "{t}";\n'
            for t in transcripts
            for s, e in EXONS[t]
        )
    )
    return path


def _ids(cell):
    return [entry.rsplit(":", 1)[0] for entry in cell.split("|") if entry]


def _run(tmp_path, transcripts, name):
    ref = mk_ref()
    alt = next(b for b in "ACGT" if b != ref[POS])
    reads = through(ref, POS, alt, 4, "a") + through(ref, POS, ref[POS], 9, "r")
    (row,) = run_rna(
        tmp_path,
        write_vcf(tmp_path, [(POS + 1, ref[POS], alt)]),
        write_bam(tmp_path, ref, reads),
        write_fasta(tmp_path, ref),
        _gtf(tmp_path, transcripts, f"{name}.gtf"),
        outname=name,
    )
    return row


@pytest.mark.parametrize(
    "order",
    [
        pytest.param(list(EXONS), id="gtf-order"),
        pytest.param(SORTED_IDS, id="id-order"),
        pytest.param(SORTED_IDS[::-1], id="reverse-id-order"),
    ],
)
def test_per_transcript_columns_list_transcripts_in_id_order(tmp_path, order):
    """Five transcripts overlap the variant; whatever order the GTF lists them in,
    both columns hold one entry per transcript, in ascending ID order."""
    row = _run(tmp_path, order, "run")
    assert _ids(row["transcript_read_counts"]) == SORTED_IDS
    assert _ids(row["transcript_fragment_counts"]) == SORTED_IDS


def test_the_cells_are_the_same_bytes_for_any_gtf_order(tmp_path):
    """The two cells are identical whichever order the GTF lists the transcripts in,
    so a byte comparison between two machines needs no set comparison."""
    cells = []
    for i, order in enumerate([list(EXONS), SORTED_IDS[::-1], SORTED_IDS]):
        row = _run(tmp_path, order, f"run{i}")
        cells.append((row["transcript_read_counts"], row["transcript_fragment_counts"]))
    assert cells[0] == cells[1] == cells[2]
    assert all(cells[0]), "the cells are not empty"
