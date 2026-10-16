"""The read-judgment cases: the executable table behind
`docs/reference/read-judgment.md`.

Each case is one read shape at one variant: four reads of that shape, counted
alone, and the call they get as (ref_count, alt_count, partial_alt). Every case
carries the decision it follows. A decided case states the decided call; an open
case pins today's call until its decision lands, so a change to either shows in
review (`tests/test_read_judgment_spec.py`). Changing an expectation is a
decision: it needs the operator's say first (AGENTS.md).

Synthetic and PHI-free. `python tests/read_judgment_cases.py` prints the table,
with each read's own verdict where the read census can judge it (pure indels).
"""

from __future__ import annotations

import random
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from census import census, tract  # noqa: E402
from helpers import make_read, write_contig  # noqa: E402

from gbcms import _rs  # noqa: E402

A, L, N_READS = 400, 900, 4
ARGS = {
    "min_mapq": 20,
    "min_baseq": 20,
    "filter_duplicates": True,
    "filter_secondary": True,
    "filter_supplementary": True,
    "filter_qc_failed": False,
    "filter_improper_pair": False,
    "filter_indel": False,
    "threads": 1,
}

# Pure indels at A (motif written from the anchor), as in the read census.
PURE = {
    "hp-A": ("GAAAAAT", "GA", "G"),
    "hp-AA": ("GAAAAAAT", "GAA", "G"),
    "hp+A": ("GAAAAAT", "G", "GA"),
    "u-8": ("GATCGGATAC", "GATCGGATA", "G"),
    "u+10": ("GT", "G", "GACGTTGCATC"),
    "ca-CA": ("G" + "CA" * 5 + "T", "GCA", "G"),
    "dup-8": ("G" + "ACGTTGCA" * 2 + "T", "GACGTTGCA", "G"),
    "dup+8": ("G" + "ACGTTGCA" * 2 + "T", "G", "GACGTTGCA"),
}
# Anchor-changing events before an A run: a whole-window event (10 A's) and a
# long one judged by junction windows (60 A's).
COMPLEX = {
    "C>TA run10": (10, "C", "TA"),
    "C>TA run60": (60, "C", "TA"),
    "CA>T run60": (60, "CA", "T"),
}

# Delins in unique sequence, read by one-sided junction windows (RJ-23): one whose
# windows are whole (28 reference bases replaced by 6) and a long one (60 by 2).
DELINS = {
    "delins 28>6": (28, "ATATGA"),
    "delins 60>2": (60, "AT"),
    "delins 17>1": (17, "A"),
}
# Delins whose contig is shaped: "delins 17>1" deletes 17 bases for one, and the ALT's
# right flank repeats the deleted bases' next two (the alleles agree for two bases
# past the one-base ALT before they differ).
SHARED_PAST_EVENT = 2
# How far past the event a recurring change sits on the ALT reads of the guard case.
RECUR_AT = 7


# What a read contributes, at an SNV (C>A at A + 1) and an unprepared Del+SNV
# (GCT>A at A, which the previous complex classifier judges): read-through bases
# past the fragment end, absent base qualities, hard clips.
# The Illumina TruSeq adapter a read runs into past a short insert.
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" * 2

READ_INPUTS = {
    "read-through: the base on the SNV is adapter, showing the ALT",
    "read-through: the base on the SNV is adapter, showing REF",
    "read-through: the read aligns on past its TLEN end (R2's 5' end clipped), showing the ALT",
    "absent base qualities (QUAL '*'), showing the ALT",
    "previous complex classifier: REF read, no clip",
    "previous complex classifier: REF read with a leading hard clip",
}


@dataclass(frozen=True)
class Case:
    id: str
    group: str
    variant: str
    shape: str


def _contig(motif: str, seed: int = 189) -> str:
    rng = random.Random(seed)
    seq = [rng.choice("CGT") for _ in range(L)]
    seq[A : A + len(motif)] = list(motif)
    return "".join(seq)


def _complex_contig(run: int) -> str:
    rng = random.Random(61)
    left = "".join(rng.choice("CGT") for _ in range(A))
    right = "".join(rng.choice("CGT") for _ in range(L - A - run - 2))
    return left + "C" + "A" * run + "G" + right


def _build(contig: str, start: int, length: int, events) -> tuple[str, tuple]:
    """A read of `length` bases from `start`, with events (pos, "I", bases) or
    (pos, "D", n) placed before their reference positions."""
    seq, cigar, p = "", [], start
    for pos, op, x in sorted(events, key=lambda e: e[0]):
        if pos > p:
            seq += contig[p:pos]
            cigar.append([0, pos - p])
            p = pos
        if op == "I":
            seq += x
            cigar.append([1, len(x)])
        elif op == "X":  # a substitution: one base
            seq += x
            cigar.append([0, 1])
            p += 1
        else:
            cigar.append([2, x])
            p += x
    seq += contig[p : p + length]
    cigar.append([0, length])
    out, q = [], 0
    for op, n in cigar:
        if q >= length:
            break
        if op == 2:
            out.append([2, n])
            continue
        take = min(n, length - q)
        out.append([op, take])
        q += take
    merged: list[list[int]] = []
    for op, n in out:
        if merged and merged[-1][0] == op == 0:
            merged[-1][1] += n
        else:
            merged.append([op, n])
    return seq[:length], tuple(map(tuple, merged))


def _own(ref: str, alt: str):
    return (A + 1, "D", len(ref) - 1) if len(ref) > len(alt) else (A + 1, "I", alt[1:])


def _pure_shapes(contig: str, ref: str, alt: str):
    """Read shapes at a pure indel: (events, read start offset, read length)."""
    lo, hi = tract(A, ref, alt, contig)
    own = _own(ref, alt)
    other = "C" if contig[hi + 4] != "C" else "G"
    in_d1 = (lo + 1, "D", 1) if hi - lo > 1 else (A + 2, "D", 1)
    # The discrimination window: the tract plus one base each side.
    w_lo, w_hi = lo - 1, hi + 1

    def placements(event) -> list[int]:
        """The positions of the event's equivalent placements: slid each way while
        the haplotype stays the same (a read's own indel is placed by its bases)."""
        pos, op, x = event
        out = [pos]
        if op == "D":
            p = pos
            while contig[p - 1] == contig[p + x - 1]:
                p -= 1
                out.append(p)
            p = pos
            while contig[p] == contig[p + x]:
                p += 1
                out.append(p)
        else:
            p, s = pos, x
            while contig[p - 1] == s[-1]:
                p, s = p - 1, s[-1] + s[:-1]
                out.append(p)
            p, s = pos, x
            while contig[p] == s[0]:
                p, s = p + 1, s[1:] + s[0]
                out.append(p)
        return out

    def where(event) -> str:
        """Inside the window when one of its placements is (RJ-22)."""
        _, op, x = event

        def inside(p: int) -> bool:
            return (w_lo < p < w_hi) if op == "I" else (p < w_hi and p + x > w_lo)

        return "inside" if any(inside(p) for p in placements(event)) else "outside"

    def is_the_alt(event) -> bool:
        """The variant itself at another placement in its tract."""
        pos, op, x = event
        if op == "D":
            return len(ref) > len(alt) and x == len(ref) - 1 and lo <= pos and pos + x <= hi
        ins = alt[1:]
        return len(alt) > len(ref) and len(x) == len(ins) and x in ins + ins and lo <= pos <= hi

    shapes = {
        "REF read spanning the tract": ([], 0, 100),
        "REF read ending inside the tract": ([], 0, (lo + 1) - (A - 50)),
        "exact carrier": ([own], 0, 100),
        "carrier with another indel outside the window": ([own, (hi + 6, "D", 1)], 0, 100),
    }
    for label, ev in (
        ("D1 near the anchor", in_d1),
        ("I1 near the anchor", (lo + 1, "I", contig[lo + 1])),
        ("D2 after the anchor", (A + 2, "D", 2)),
        ("D1 past the tract", (hi + 4, "D", 1)),
        ("D5 past the tract", (hi + 4, "D", 5)),
        ("I1 past the tract", (hi + 4, "I", other)),
    ):
        if is_the_alt(ev):
            shapes[f"the ALT at another placement: {label}"] = ([ev], 0, 100)
        else:
            shapes[f"another indel {where(ev)} the window: {label}"] = ([ev], 0, 100)
    if len(ref) > 2:  # a deletion of two or more bases written as two ops
        k = len(ref) - 1
        shapes["the ALT written as two deletions"] = (
            [(A + 1, "D", 1), (A + 2, "D", k - 1)],
            0,
            100,
        )
    run_base = contig[A + 1]
    in_run = contig[A + 1 : A + 6] == run_base * 5
    if in_run and len(ref) == 2 and len(alt) == 1:  # -1 in a run as a shifted pair
        shapes["the ALT written as D2 + I1 at the anchor"] = (
            [(A + 1, "D", 2), (A + 3, "I", run_base)],
            0,
            100,
        )
        shapes["the ALT written as D2 + I1 inside the run"] = (
            [(A + 2, "D", 2), (A + 4, "I", run_base)],
            0,
            100,
        )
    if in_run and len(ref) == 1 and len(alt) == 2:  # +1 in a run as a shifted pair
        shapes["the ALT written as I2 + D1 at the anchor"] = (
            [(A + 1, "I", run_base * 2), (A + 1, "D", 1)],
            0,
            100,
        )
    if len(ref) > 3:  # the anchor deleted
        shapes["anchor deleted, then an insertion"] = ([(A, "D", 3), (A + 5, "I", "T")], 0, 100)
        shapes["anchor deleted"] = ([(A, "D", 3)], 0, 100)
    if len(alt) > 2:  # an insertion row whose read deletes the anchor
        shapes["anchor deleted"] = ([(A, "D", 2)], 0, 100)
    return shapes


def _unreadable_shapes(contig: str, ref: str, alt: str):
    """Carriers whose inserted bases are masked: (events, (kind, first, count)), the
    read's inserted bases [first, first + count) masked as N at Q2 ("N") or kept at
    Q5 ("low")."""
    ins = alt[1:]
    if len(ins) == 1:  # +1 in a run: the aligner writes the masked base as the insertion
        return {
            "an unreadable insert at the anchor (N)": ([(A + 1, "I", "N")], ("N", 0, 1)),
            "an unreadable insert inside the run (N)": ([(A + 3, "I", "N")], ("N", 0, 1)),
            "an insert below min BQ inside the run": ([(A + 3, "I", ins)], ("low", 0, 1)),
        }
    k = len(ins)
    return {
        "an unreadable insert (N)": ([(A + 1, "I", ins)], ("N", 0, k)),
        "an insert below min BQ": ([(A + 1, "I", ins)], ("low", 0, k)),
        "a partly masked insert": ([(A + 1, "I", ins)], ("N", 0, k // 2)),
        "anchor deleted, the insert unreadable": (
            [(A, "D", 1), (A + 1, "I", contig[A] + ins)],
            ("N", 1, k),
        ),
    }


def _other_bases_shapes(contig: str, ref: str, alt: str):
    """Same-length inserts of other bases near a pure insertion, and the ALT written
    one junction off with a compensating mismatch: read-judgment events."""
    ins = alt[1:]
    other = "CTTAGCCTAG"[: len(ins)]
    if len(ins) == 8:  # dup+8: inside its duplicated tract
        return {"a same-length insert of other bases inside the tract": [(A + 3, "I", other)]}
    return {
        "a same-length insert of other bases one junction left": [(A, "I", other)],
        "a same-length insert of other bases four junctions right": [(A + 5, "I", other)],
        "the ALT written one junction off with a compensating mismatch": [
            (A + 1, "X", ins[0]),
            (A + 2, "I", ins[1:] + contig[A + 1]),
        ],
        # Its first inserted base misread as the reference base after the anchor:
        # the aligner writes it one junction right with no mismatch; it slides back.
        "the ALT with an error at its first inserted base, written one junction right": [
            (A + 2, "I", ins[1:] + contig[A + 1]),
        ],
    }


def _behind_substitution_shapes(contig: str, ref: str, alt: str):
    """The ALT with a substitution next to its insert: at the junction with a
    substitution in the run (the strict path keeps it), and with the anchor
    substituted, written one junction left (no slide carries it back): events."""
    ins = alt[1:]
    out = {
        # The aligner's tie: one mismatch either way, at the anchor; the first
        # inserted base (the read's anchor, then the ALT's first bases) goes left.
        "the ALT with its anchor substituted, written one junction left": [
            (A, "I", "T" + ins[:-1]),
            (A, "X", ins[-1]),
        ],
    }
    if len(set(ins)) == 1:  # in a run: a substitution the slide would absorb
        out["the ALT at the junction with a substitution in the run"] = [
            (A + 1, "I", ins),
            (A + 3, "X", "C"),
        ]
    return out


def _mask_inserted(read, mask) -> None:
    """Mask the read's inserted bases [first, first + count), in CIGAR order."""
    kind, first, count = mask
    seq, quals = list(read.query_sequence), list(read.query_qualities)
    q, k = 0, 0
    for op, n in read.cigartuples:
        if op == 1:
            for i in range(q, q + n):
                if first <= k < first + count:
                    seq[i], quals[i] = ("N", 2) if kind == "N" else (seq[i], 5)
                k += 1
        if op in (0, 1, 4):
            q += n
    read.query_sequence = "".join(seq)
    read.query_qualities = quals


def _complex_reads(run: int, ref_allele: str, alt_allele: str):
    """Read haplotypes at an anchor-changing event before an A run, ending inside
    the run, on the base after it, or past it: (events, last reference base)."""
    end_run = A + 1 + run  # the G after the run
    haps = {
        "REF": [],
        "substitution only": [(A, "X", "T")],
        "anchor kept, one A more": [(A + 1, "I", "A")],
        "anchor kept, one A fewer": [(A + 1, "D", 1)],
        "exact carrier": [(A, "X", "T")]
        + ([(A + 1, "I", "A")] if alt_allele == "TA" else [(A + 1, "D", 1)]),
    }
    ends = {
        "ends inside the run": A + run // 2,
        "ends on the base after the run": end_run,
        "ends past the run": end_run + 20,
    }
    return haps, ends


def _grows(r: str, b: int, unit: int, left: bool) -> bool:
    """Whether a tandem repeat (unit 1..`unit`, at least one whole unit) runs from
    `b` leftward (`left`) or rightward: carrier.rs's repeat growth, so a shaped case
    can avoid it."""
    for k in range(1, unit + 1):
        i, n = b, 0
        if left:
            while i >= 1 and i - 1 + k < len(r) and r[i - 1] == r[i - 1 + k]:
                i, n = i - 1, n + 1
        else:
            while i < len(r) and i >= k and r[i] == r[i - k]:
                i, n = i + 1, n + 1
        if n >= k:
            return True
    return False


def _delins_contig(variant: str) -> str:
    """The delins case's contig: random C/G/T; for "delins 17>1", the first seed
    where the ALT's right flank repeats the deleted bases after the first for
    SHARED_PAST_EVENT bases, then differs, and no repeat grows the event on either
    allele (the windows stay the short allele's own)."""
    if variant != "delins 17>1":
        return _contig("")
    n, alt = DELINS[variant]
    for seed in range(189, 1189):
        seq = list(_contig("", seed))
        for k in range(1, SHARED_PAST_EVENT + 1):
            seq[A + k] = seq[A + n + k - 1]
        k = SHARED_PAST_EVENT + 1
        if seq[A + k] == seq[A + n + k - 1]:
            continue
        ref = "".join(seq)
        hap = ref[:A] + alt + ref[A + n :]
        unit = max(6, n - len(alt))
        if not (
            _grows(ref, A, unit, True)
            or _grows(hap, A, unit, True)
            or _grows(ref, A + n, unit, False)
            or _grows(hap, A + len(alt), unit, False)
        ):
            return ref
    raise RuntimeError("no seed shapes the delins 17>1 contig")


def _delins_reads(n: int, alt: str, shape: str, contig: str = ""):
    """Reads of a delins shape: (start, query length, events), the delins replacing
    reference [A, A + n) by `alt`, starts and ends varying per read so no window
    length is special; or ("raw", position, bases, cigar, qualities) for a read the
    builder cannot write (a low-quality clip, bases aligned where the aligner put a
    carrier at the event's other copy)."""
    carrier = [(A, "I", alt), (A, "D", n)]
    m = len(alt)
    hap = contig[:A] + alt + contig[A + n :]
    out = []
    for i in range(N_READS):
        far = 60 + 3 * i
        if shape == "exact carrier, eight flank bases right":
            # eight: the 28-to-6 event grows three bases right through a TCG repeat
            out.append((A - far, far + m + 8, carrier))
        elif shape == "exact carrier, four flank bases left":
            out.append((A - 4, 4 + m + far, carrier))
        elif shape == "exact carrier, long flanks":
            out.append((A - far, far + m + 60, carrier))
        elif shape == "REF read, four flank bases left":
            out.append((A - 4, 4 + n + far, []))
        elif shape == "REF read, four flank bases right":
            out.append((A - far, far + n + 4, []))
        elif shape == "REF read starting inside the event":
            out.append((A + n // 2, n - n // 2 + far, []))
        elif shape == "carriers ending before a change that recurs on ALT reads only":
            # four carriers end four bases past the event (inside its right flank, as
            # it grows three bases through a TCG repeat: depth only); four more read
            # on past a change seven bases out, which every one of them shows
            out.append((A - far, far + m + 4, carrier))
            sub = (A + n + RECUR_AT - 1, "X", "A")
            out.append((A - far, far + m + 40, carrier + [sub]))
        elif shape == "REF read with the event's first base substituted":
            # the ALT's own base where REF's first is, then REF's bases on: a
            # substitution, not the delins (its next bases are not the ALT's flank)
            out.append((A - far, far + n + 40, [(A, "X", alt)]))
        elif shape == "carriers with a recurring change on the window's last base":
            # the 28-to-6 event grows three bases right, so its window's last base is
            # the right flank's fifth: masked (an isolated change), the carriers count
            out.append((A - far, far + m + 40, carrier + [(A + n + 4, "X", "A")]))
        elif shape == "carriers with a run of changes from the window's last base on":
            # the window's last base and the next: a run (another allele), not masked
            out.append(
                (A - far, far + m + 40, carrier + [(A + n + 4, "X", "A"), (A + n + 5, "X", "A")])
            )
        elif shape == "REF read starting inside the event behind a low-quality clip":
            # twelve REF bases before the event's end, soft-clipped at Q2, the rest
            # aligned from its first right-flank base: no flank is read
            seq = contig[A + n - 12 : A + n - 12 + 100]
            quals = [2] * 12 + [30] * 88
            out.append(("raw", A + n, seq, ((4, 12), (0, 88)), quals))
        elif shape == "carrier written at the event's other copy as plain matches":
            # an aligner placing a carrier by its right part: its flank and ALT bases
            # aligned over the REF event's last bases as mismatches, no gap
            k = 8 + i  # flank bases before the event, in the read
            seq = hap[A - k : A - k + 100]
            out.append(("raw", A - k + (n - m), seq, ((0, 100),), [30] * 100))
        else:
            raise KeyError(shape)
    return out


# Pure deletions of 50 bp or more, judged at the POS breakpoint (RJ-24): (deleted
# bases, microhomology). The contig is longer, so a 1 kb deletion fits; reads are 100
# bases, each the forward read of a proper pair (its fragment admits clipped bases).
LONG_DEL = {"del 60": (60, 0), "del 60 mh3": (60, 3), "del 1000": (1000, 0)}
LONG_L, LONG_READ = 2000, 100
# Where a recurring change sits past the junction, and inside REF's window.
RECUR_PAST = 1


def _long_contig(n: int, mh: int) -> str:
    """Random C/G/T with the deletion [A, A + n) sliding exactly `mh` bases right and
    none left."""
    rng = random.Random(254)
    seq = [rng.choice("CGT") for _ in range(LONG_L)]
    for i in range(mh):
        seq[A + i] = seq[A + n + i]
    if seq[A + mh] == seq[A + n + mh]:
        seq[A + mh] = next(b for b in "CGT" if b != seq[A + n + mh])
    if seq[A - 1] == seq[A + n - 1]:
        seq[A - 1] = next(b for b in "CGT" if b != seq[A + n - 1])
    return "".join(seq)


def _other(base: str) -> str:
    return "A" if base != "A" else "C"


def _long_del_reads(n: int, mh: int, shape: str) -> list:
    """The four reads of a long-deletion shape, as pysam records (pairs' forward reads;
    a supplementary record where the shape has one)."""
    contig = _long_contig(n, mh)
    hap = contig[:A] + contig[A + n :]  # ALT: genome g >= A + n sits at g - n
    win = mh + 4  # past the junction: the microhomology and four clean bases
    reads = []

    def add(name, seq, pos, cigar, flag=0, quals=None):
        r = make_read(
            name, seq, pos, tuple(cigar), flag=0x1 | 0x2 | 0x20 | 0x40 | flag, quals=quals
        )
        r.next_reference_id = 0
        r.next_reference_start = pos + 250
        r.template_length = 400
        reads.append(r)

    def carrier(i, s, right, subs=()):
        """An exact D carrier from `s` holding `right` bases past the junction, with a
        third base at each genome position in `subs`."""
        left = A - s
        seq = list(hap[s : A + right])
        for g in subs:
            k = g - s if g < A else g - n - s
            seq[k] = _other(seq[k])
        add(f"c{i}", "".join(seq), s, [(0, left), (2, n), (0, right)])

    def ref_read(i, s, length, sub=None):
        seq = list(contig[s : s + length])
        if sub is not None:
            seq[sub - s] = _other(seq[sub - s])
        add(f"r{i}", "".join(seq), s, [(0, length)])

    for i in range(N_READS):
        if shape == "exact D carrier, long flanks":
            carrier(i, A - 50 - i, LONG_READ - 50 - i)
        elif shape == "exact D carrier, four bases past the microhomology":
            carrier(i, A - (LONG_READ - win), win)
        elif shape == "exact D carrier, three bases past the microhomology":
            carrier(i, A - (LONG_READ - win + 1), win - 1)
        elif shape.startswith("right-clipped carrier"):
            s = A - 70 - i
            seq = hap[s : s + LONG_READ]
            if shape.endswith("a recurring change in the flank"):
                seq = seq[: A - 2 - s] + _other(seq[A - 2 - s]) + seq[A - 1 - s :]
            add(f"c{i}", seq, s, [(0, A - s), (4, LONG_READ - (A - s))])
        elif shape.startswith("left-clipped carrier"):
            c = 20 + i  # bases before the junction, in the clip
            quals = [30] * LONG_READ
            if shape.endswith("two low-quality flank bases"):
                quals[c - 4 : c - 2] = [2, 2]  # the base next to the junction stays read
            elif shape.endswith("a low-quality flank"):
                quals[c - 4 : c] = [2] * 4
            seq = hap[A - c : A - c + LONG_READ]
            add(f"c{i}", seq, A + n, [(4, c), (0, LONG_READ - c)], quals=quals)
        elif shape == "split carrier: a soft-clipped primary and a hard-clipped supplementary":
            s = A - 70 - i
            seq = hap[s : s + LONG_READ]
            add(f"c{i}", seq, s, [(0, A - s), (4, LONG_READ - (A - s))])
            add(f"c{i}", seq[A - s :], A + n, [(5, A - s), (0, LONG_READ - (A - s))], flag=0x800)
        elif shape == "hard-clipped primary aligned from the right flank":
            c = 20 + i
            add(f"c{i}", hap[A : A - c + LONG_READ], A + n, [(5, c), (0, LONG_READ - c)])
        elif shape == "REF read, long flanks":
            ref_read(i, A - 50 - i, LONG_READ)
        elif shape == "REF read, four deleted bases past the microhomology":
            ref_read(i, A + win - LONG_READ, LONG_READ)
        elif shape == "REF read, three deleted bases past the microhomology":
            ref_read(i, A + win - 1 - LONG_READ, LONG_READ)
        elif shape == "REF read with an error 20 bases past the window":
            ref_read(i, A - 50 - i, LONG_READ, sub=A + win + 20)
        elif shape == "REF read with an error inside the window":
            # an error is one read's: the other three are clean (one shared by every
            # read is a recurring change, masked)
            ref_read(i, A - 50 - i, LONG_READ, sub=A + mh + 1 if i == 0 else None)
        elif shape == "a deletion two bases shorter at the same start":
            s = A - 50 - i
            other = contig[:A] + contig[A + n - 2 :]
            add(
                f"c{i}",
                other[s : s + LONG_READ],
                s,
                [(0, A - s), (2, n - 2), (0, LONG_READ - (A - s))],
            )
        elif shape == "the deletion written as two with two bases kept between":
            s = A - 50 - i
            kept = contig[:A] + contig[A + 30 : A + 32] + contig[A + n :]
            rest = LONG_READ - (A - s) - 2
            add(
                f"c{i}",
                kept[s : s + LONG_READ],
                s,
                [(0, A - s), (2, 30), (0, 2), (2, n - 32), (0, rest)],
            )
        elif shape == "REF read with its own 20-bp deletion inside the deleted span":
            s = A - 50 - i
            other = contig[: A + 20] + contig[A + 40 :]
            add(
                f"r{i}",
                other[s : s + LONG_READ],
                s,
                [(0, A + 20 - s), (2, 20), (0, LONG_READ - (A + 20 - s))],
            )
        elif shape == "REF read starting inside the deletion":
            ref_read(i, A + 10 + i, LONG_READ)
        elif shape.startswith("REF read starting inside the deletion behind"):
            # six bases before the right flank, clipped: they would be the flank of a
            # carrier placed from the right, but are read too poorly to show it
            s = A + n - 6
            seq, quals = contig[s : s + LONG_READ], [30] * LONG_READ
            if "N clip" in shape:
                seq = "N" * 6 + seq[6:]
            else:
                quals[:6] = [2] * 6
            add(f"r{i}", seq, A + n, [(4, 6), (0, LONG_READ - 6)], quals=quals)
        elif shape == "carriers with a recurring change two bases past the microhomology":
            carrier(i, A - 50 - i, LONG_READ - 50 - i, subs=(A + n + mh + RECUR_PAST,))
        elif shape == "carriers with a recurring change in the flank":
            carrier(i, A - 50 - i, LONG_READ - 50 - i, subs=(A - 2,))
        elif shape == "carriers with a run of changes in the flank":
            carrier(i, A - 50 - i, LONG_READ - 50 - i, subs=(A - 3, A - 2))
        elif shape == "REF reads with a run of changes in the flank":
            seq = list(contig[A - 50 - i : A - 50 - i + LONG_READ])
            for g in (A - 3, A - 2):
                seq[g - (A - 50 - i)] = _other(seq[g - (A - 50 - i)])
            add(f"r{i}", "".join(seq), A - 50 - i, [(0, LONG_READ)])
        elif shape == "carriers and REF reads with a run of changes in the flank":
            carrier(i, A - 50 - i, LONG_READ - 50 - i, subs=(A - 3, A - 2))
            seq = list(contig[A - 50 - i : A - 50 - i + LONG_READ])
            for g in (A - 3, A - 2):
                seq[g - (A - 50 - i)] = _other(seq[g - (A - 50 - i)])
            add(f"r{i}", "".join(seq), A - 50 - i, [(0, LONG_READ)])
        elif shape == "carrier aligned over the deleted bases' end, no clip":
            # its bases before the junction written over the last deleted bases, as
            # mismatches: an aligner placing it by its right part, without a clip
            c = 6 + i
            add(f"c{i}", hap[A - c : A - c + LONG_READ], A + n - c, [(0, LONG_READ)])
        elif shape == "carriers with a recurring change on the window's last base":
            carrier(i, A - 50 - i, LONG_READ - 50 - i, subs=(A + n + mh + 3,))
        elif shape == "carriers with a run of changes from the window's last base on":
            carrier(i, A - 50 - i, LONG_READ - 50 - i, subs=(A + n + mh + 3, A + n + mh + 4))
        elif shape == "REF reads with a recurring change inside the window":
            ref_read(i, A - 50 - i, LONG_READ, sub=A + mh + RECUR_PAST)
        elif shape == "REF reads with a run of changes from the window's last base on":
            seq = list(contig[A - 50 - i : A - 50 - i + LONG_READ])
            for g in (A + mh + 3, A + mh + 4):
                seq[g - (A - 50 - i)] = _other(seq[g - (A - 50 - i)])
            add(f"r{i}", "".join(seq), A - 50 - i, [(0, LONG_READ)])
        else:
            raise KeyError(shape)
    return reads


_LONG_DEL_ALL = (
    "exact D carrier, long flanks",
    "exact D carrier, four bases past the microhomology",
    "exact D carrier, three bases past the microhomology",
    "right-clipped carrier, the junction in the clip",
    "left-clipped carrier aligned from the right flank",
    "left-clipped carrier with two low-quality flank bases",
    "left-clipped carrier with a low-quality flank",
    "split carrier: a soft-clipped primary and a hard-clipped supplementary",
    "hard-clipped primary aligned from the right flank",
    "REF read, long flanks",
    "REF read, four deleted bases past the microhomology",
    "REF read, three deleted bases past the microhomology",
    "REF read with an error 20 bases past the window",
    "REF read with an error inside the window",
    "a deletion two bases shorter at the same start",
    "the deletion written as two with two bases kept between",
    "REF read with its own 20-bp deletion inside the deleted span",
    "REF read starting inside the deletion",
    "REF read starting inside the deletion behind a low-quality clip",
    "REF read starting inside the deletion behind an N clip",
    "carriers with a recurring change two bases past the microhomology",
    "REF reads with a recurring change inside the window",
    "carriers with a recurring change in the flank",
    "right-clipped carriers with a recurring change in the flank",
    "carriers with a run of changes in the flank",
    "REF reads with a run of changes in the flank",
    "carriers and REF reads with a run of changes in the flank",
    "REF reads with a run of changes from the window's last base on",
    "carriers with a recurring change on the window's last base",
    "carriers with a run of changes from the window's last base on",
    "carrier aligned over the deleted bases' end, no clip",
)
LONG_DEL_SHAPES = {
    "del 60": _LONG_DEL_ALL,
    "del 60 mh3": (
        "exact D carrier, long flanks",
        "exact D carrier, four bases past the microhomology",
        "exact D carrier, three bases past the microhomology",
        "left-clipped carrier aligned from the right flank",
        "REF read, long flanks",
        "REF read, four deleted bases past the microhomology",
        "REF read, three deleted bases past the microhomology",
    ),
    "del 1000": (
        "exact D carrier, long flanks",
        "right-clipped carrier, the junction in the clip",
        "left-clipped carrier aligned from the right flank",
        "REF read, long flanks",
        "REF read starting inside the deletion",
    ),
}


DELINS_SHAPES = {
    "delins 28>6": (
        "exact carrier, eight flank bases right",
        "exact carrier, four flank bases left",
        "exact carrier, long flanks",
        "REF read, four flank bases left",
        "REF read, four flank bases right",
        "REF read starting inside the event",
        "carriers ending before a change that recurs on ALT reads only",
    ),
    "delins 60>2": (
        "exact carrier, eight flank bases right",
        "exact carrier, long flanks",
        "REF read, four flank bases left",
        "REF read starting inside the event",
    ),
    "delins 17>1": (
        "exact carrier, long flanks",
        "REF read, four flank bases right",
        "REF read with the event's first base substituted",
    ),
}
# RJ-23 amendments (2026-10-12), on the 28-to-6 delins.
DELINS_SHAPES["delins 28>6"] += (
    "carriers with a recurring change on the window's last base",
    "carriers with a run of changes from the window's last base on",
    "REF read starting inside the event behind a low-quality clip",
    "carrier written at the event's other copy as plain matches",
)

# Other alleles of the ALT's length at the locus (RJ-25), and a change inside the repeat
# an event is grown through (RJ-23, 2026-10-14). Each case is its shape's four reads
# with the reads of the locus they are judged beside.
# The 28-to-6 delins, with another allele one inserted base apart.
OTHER_INSERT = "ATCTGA"
# An MNP joining a somatic change (G>A) to a germline SNP beside it (C>T): the germline
# allele alone is GT.
MNP = ("GC", "AT", "GT")
# An insertion grown through a CA repeat: G then CACACA, a germline A>G at the repeat's
# fourth base.
CA_RUN = ("GCACACAT", "G", "TCA", 4)
LOCUS = {
    "delins 28>6, another allele": (
        "carriers with the telling base below min BQ, the other allele on four fragments",
        "carriers reading the telling base, the other allele on four fragments",
        "carriers with the telling base below min BQ, the other allele on two fragments",
        "carriers with the telling base below min BQ, no other allele",
        "carriers with the telling base below min BQ, the other allele co-annotated",
    ),
    "MNP GC>AT, germline GT": (
        "reads showing only the germline base, the germline allele on four fragments",
        "reads showing only the somatic base, the germline allele on four fragments",
        "reads showing only the germline base, no germline allele",
    ),
    "G>TCA in a CA run, germline SNP": (
        "REF reads with the SNP, beside REF reads without it",
        "carriers with the SNP, beside REF reads with it",
        "carriers with the SNP on the carriers alone",
    ),
}


def _quals_at(read, at: list[int], q: int = 5) -> None:
    """Set the bases at query positions `at` to quality `q`, keeping their letters."""
    quals = list(read.query_qualities)
    for i in at:
        quals[i] = q
    read.query_qualities = quals


def _locus_reads(variant: str, shape: str):
    """A locus case: (contig, reads, rows (pos, ref, alt), co-annotated rows). Reads
    of the shape are r0..r3; the locus's other reads follow."""
    reads = []

    def add(tag, contig, s, length, events, low=()):
        seq, cig = _build(contig, s, length, events)
        r = make_read(f"{tag}{len(reads)}", seq, s, cig)
        if low:
            _quals_at(r, [i - s for i in low])
        reads.append(r)

    if variant == "delins 28>6, another allele":
        n, alt = DELINS["delins 28>6"]
        contig = _delins_contig("delins 28>6")
        telling = next(k for k in range(len(alt)) if alt[k] != OTHER_INSERT[k])
        others = 4 if "four fragments" in shape else 2 if "two fragments" in shape else 0
        masked = "below min BQ" in shape
        for i in range(N_READS):
            far = 60 + 3 * i
            seq, cig = _build(contig, A - far, far + len(alt) + 40, [(A, "I", alt), (A, "D", n)])
            r = make_read(f"r{i}", seq, A - far, cig)
            if masked:
                _quals_at(r, [far + telling])  # the inserted base after `far` flank bases
            reads.append(r)
        for i in range(others):
            far = 61 + 3 * i
            add("o", contig, A - far, far + len(alt) + 40, [(A, "I", OTHER_INSERT), (A, "D", n)])
        rows = [(A, contig[A : A + n], alt)]
        co = [(A, contig[A : A + n], OTHER_INSERT)] if "co-annotated" in shape else []
        return contig, reads, rows, co
    if variant == "MNP GC>AT, germline GT":
        ref, alt, germ = MNP
        contig = _contig(ref)
        only = "germline" if "only the germline base" in shape else "somatic"
        for i in range(N_READS):
            s = A - 50 + i
            seq = contig[s:A] + alt + contig[A + 2 : s + 100]
            r = make_read(f"r{i}", seq, s, ((0, 100),))
            # showing only the germline base: the somatic base (A) below min BQ
            _quals_at(r, [A - s if only == "germline" else A + 1 - s])
            reads.append(r)
        if "on four fragments" in shape:
            for i in range(N_READS):
                s = A - 49 + i
                reads.append(
                    make_read(f"g{i}", contig[s:A] + germ + contig[A + 2 : s + 100], s, ((0, 100),))
                )
        return contig, reads, [(A, ref, alt)], []
    motif, ref, alt, snp = CA_RUN
    contig = _contig(motif)
    snp_base = "G"
    carrier = [(A, "X", alt[0]), (A + 1, "I", alt[1:])]
    with_snp = [(A + snp, "X", snp_base)]
    if shape == "REF reads with the SNP, beside REF reads without it":
        shape_events, beside = with_snp, []
    elif shape == "carriers with the SNP, beside REF reads with it":
        shape_events, beside = carrier + with_snp, with_snp
    else:
        shape_events, beside = carrier + with_snp, []
    for i in range(N_READS):
        s = A - 50 + i
        seq, cig = _build(contig, s, 100, shape_events)
        reads.append(make_read(f"r{i}", seq, s, cig))
    for i in range(N_READS):
        add("b", contig, A - 49 + i, 100, beside)
    return contig, reads, [(A, ref, alt)], []


def cases() -> list[Case]:
    out: list[Case] = []
    for v, (motif, ref, alt) in PURE.items():
        for shape in _pure_shapes(_contig(motif), ref, alt):
            group = _pure_group(shape)
            out.append(Case(f"{v} | {shape}", group, v, shape))
    for v, (run, r, a) in COMPLEX.items():
        haps, ends = _complex_reads(run, r, a)
        for h in haps:
            for e in ends:
                group = "C25 long events" if run > 40 else "C1 exact carriers"
                out.append(Case(f"{v} | {h}, {e}", group, v, f"{h}, {e}"))
    for v, shapes in DELINS_SHAPES.items():
        for shape in shapes:
            out.append(Case(f"{v} | {shape}", "C39 one-sided windows", v, shape))
    for v, shapes in LONG_DEL_SHAPES.items():
        for shape in shapes:
            out.append(Case(f"{v} | {shape}", "C40 long deletions", v, shape))
    for v, shapes in LOCUS.items():
        group = (
            "C43 background in the grown event"
            if "CA run" in v
            else "C43 other alleles at the locus"
        )
        for shape in shapes:
            out.append(Case(f"{v} | {shape}", group, v, shape))
    out.append(
        Case(
            "u-8 | sibling SNV inside the span", "C2 siblings", "u-8", "sibling SNV inside the span"
        )
    )
    out.append(
        Case(
            "u-8 | sibling SNV outside the window",
            "C2 siblings",
            "u-8",
            "sibling SNV outside the window",
        )
    )
    for shape in sorted(READ_INPUTS):
        out.append(Case(f"read inputs | {shape}", "read inputs", "read inputs", shape))
    for v in ("hp+A", "u+10"):
        motif, ref, alt = PURE[v]
        for shape in _unreadable_shapes(_contig(motif), ref, alt):
            out.append(Case(f"{v} | {shape}", "C35 unreadable inserts", v, shape))
    for v in ("u+10", "dup+8"):
        motif, ref, alt = PURE[v]
        for shape in _other_bases_shapes(_contig(motif), ref, alt):
            out.append(Case(f"{v} | {shape}", "C36 inserts of other bases", v, shape))
    for v in ("hp+A", "u+10"):
        motif, ref, alt = PURE[v]
        for shape in _behind_substitution_shapes(_contig(motif), ref, alt):
            out.append(Case(f"{v} | {shape}", "C37 inserts behind a substitution", v, shape))
    return out


def _read_input_bam(d: Path, shape: str):
    """The reads of a read-input case, the contig, and the variant they cover."""
    contig = _contig("GCT")
    reads = []
    if shape.startswith("read-through"):
        # R1's molecule ends just before the SNV and it reads on into adapter: the
        # adapter's first base is aligned on the SNV by chance (showing the ALT or
        # REF), the rest soft-clipped. With "aligns on", its bases past TLEN's end
        # are the genome's instead (R2's 5' end clipped): the molecule goes on.
        base = "A" if "ALT" in shape else contig[A + 1]
        on = "aligns on" in shape
        for i in range(N_READS):
            s, end = A - 60 + i, (A - 9 if on else A + 1)
            if on:
                r1 = list(contig[s : s + 100])
                r1[A + 1 - s] = base
                seq, cig = "".join(r1), ((0, 100),)
                r2, r2_cig = contig[end - 80 : end + 20], ((0, 80), (4, 20))
            else:
                tail = ADAPTER[1 : 100 - (A + 2 - s) + 1]
                seq, cig = contig[s : A + 1] + base + tail, ((0, A + 2 - s), (4, len(tail)))
                r2, r2_cig = contig[end - 100 : end], ((0, 100),)
            a = make_read(f"f{i}", seq, s, cig, flag=99)
            b = make_read(f"f{i}", r2, end - len(r2) + (20 if on else 0), r2_cig, flag=147)
            a.template_length, b.template_length = end - s, s - end
            a.next_reference_id = b.next_reference_id = 0
            a.next_reference_start, b.next_reference_start = b.reference_start, s
            reads += [a, b]
        return contig, reads, _rs.Variant("1", A + 1, contig[A + 1], "A", "SNP")
    if shape.startswith("absent"):
        for i in range(N_READS):
            s = A - 50 + i
            seq = list(contig[s : s + 100])
            seq[A + 1 - s] = "A"
            a = make_read(f"q{i}", "".join(seq), s, ((0, 100),))
            a.query_qualities = None  # QUAL '*': the BAM stores 0xFF
            reads.append(a)
        return contig, reads, _rs.Variant("1", A + 1, contig[A + 1], "A", "SNP")
    clip = "hard clip" in shape
    for i in range(N_READS):
        s = A - 60 + i
        cig = ((5, 20), (0, 70)) if clip else ((0, 70),)
        reads.append(make_read(f"h{i}", contig[s : s + 70], s, cig))
    # Unprepared (no reference context): the exact-carrier rule cannot judge it.
    return contig, reads, _rs.Variant("1", A, contig[A : A + 3], "A", "COMPLEX")


def _pure_group(shape: str) -> str:
    if shape.startswith("REF read"):
        return "C10 REF reads"
    if shape == "exact carrier" or shape.startswith("the ALT at another placement"):
        return "invariant 7 carriers"
    if "inside the window" in shape:
        return "C26 inside the window"
    if "outside the window" in shape:
        return "C26 outside the window"
    if shape.startswith("the ALT written as"):
        return "C27 the ALT across ops"
    return "C28 anchor deleted"


def run(case: Case, workdir: Path) -> tuple[tuple[int, int, int], str]:
    """Count the case's reads; returns ((rd, ad, partial), census verdicts)."""
    d = Path(tempfile.mkdtemp(dir=workdir))
    if case.variant in PURE:
        motif, ref, alt = PURE[case.variant]
        contig = _contig(motif)
        siblings: list = []
        mask = None
        if case.group == "C2 siblings":
            snv_at = A + 3 if "inside" in case.shape else tract(A, ref, alt, contig)[1] + 6
            snv_alt = "C" if contig[snv_at] != "C" else "G"
            events, s0, length = [(snv_at, "X", snv_alt)], 0, 100
        elif case.group == "C35 unreadable inserts":
            (events, mask), s0, length = _unreadable_shapes(contig, ref, alt)[case.shape], 0, 100
        elif case.group == "C36 inserts of other bases":
            events, s0, length = _other_bases_shapes(contig, ref, alt)[case.shape], 0, 100
        elif case.group == "C37 inserts behind a substitution":
            events, s0, length = _behind_substitution_shapes(contig, ref, alt)[case.shape], 0, 100
        else:
            events, s0, length = _pure_shapes(contig, ref, alt)[case.shape]
        reads = []
        for i in range(N_READS):
            s = A - 50 + i + s0
            seq, cig = _build(contig, s, length - i, events)
            reads.append(make_read(f"r{i}", seq, s, cig))
            if mask:
                _mask_inserted(reads[-1], mask)
        fa, bam = write_contig(d, contig, reads, "c")
        rows = [_rs.Variant("1", A, ref, alt, "X")]
        if case.group == "C2 siblings":
            rows.append(_rs.Variant("1", snv_at, contig[snv_at], snv_alt, "SNP"))
        pvs = _rs.prepare_variants(rows, fa, 5, False, 1, True)
        vs = [p.variant for p in pvs]
        if case.group == "C2 siblings":
            siblings = [[vs[1]], [vs[0]]]
        counts = _rs.count_bam_binned(
            bam, vs, [None] * len(vs), sibling_variants=siblings or [[] for _ in vs], **ARGS
        )
        c = counts[0]
        verdicts = Counter(x.value for x in census(bam, contig, vs[0]).counts.elements())
        return (c.rd, c.ad, c.partial_alt), ", ".join(
            f"{k} {n}" for k, n in sorted(verdicts.items())
        )
    if case.variant == "read inputs":
        contig, reads, v = _read_input_bam(d, case.shape)
        _, bam = write_contig(d, contig, reads, "c")
        (c,) = _rs.count_bam_binned(bam, [v], [None], **ARGS)
        return (c.rd, c.ad, c.partial_alt), ""
    if case.variant in DELINS:
        n, alt = DELINS[case.variant]
        contig = _delins_contig(case.variant)
        reads = []
        for i, read in enumerate(_delins_reads(n, alt, case.shape, contig)):
            if read[0] == "raw":
                _, s, seq, cig, quals = read
                r = make_read(f"r{i}", seq, s, cig, flag=0x1 | 0x2 | 0x20 | 0x40, quals=quals)
                r.next_reference_id, r.next_reference_start, r.template_length = 0, s + 150, 300
                reads.append(r)
                continue
            s, length, events = read
            seq, cig = _build(contig, s, length, events)
            reads.append(make_read(f"r{i}", seq, s, cig))
        fa, bam = write_contig(d, contig, reads, "c")
        (pv,) = _rs.prepare_variants(
            [_rs.Variant("1", A, contig[A : A + n], alt, "X")], fa, 5, False, 1, True
        )
        (c,) = _rs.count_bam_binned(bam, [pv.variant], [None], **ARGS)
        return (c.rd, c.ad, c.partial_alt), ""
    if case.variant in LOCUS:
        contig, reads, rows, co = _locus_reads(case.variant, case.shape)
        fa, bam = write_contig(d, contig, reads, "c")
        given = [_rs.Variant("1", p, r, a, "X") for p, r, a in rows + co]
        vs = [pv.variant for pv in _rs.prepare_variants(given, fa, 5, False, 1, True)]
        siblings = [[w for w in vs if w is not v] for v in vs]
        counts = _rs.count_bam_binned(bam, vs, [None] * len(vs), sibling_variants=siblings, **ARGS)
        c = counts[0]
        return (c.rd, c.ad, c.partial_alt), ""
    if case.variant in LONG_DEL:
        n, mh = LONG_DEL[case.variant]
        contig = _long_contig(n, mh)
        _, bam = write_contig(d, contig, _long_del_reads(n, mh, case.shape), "c")
        fa = str(d / "c.fa")
        (pv,) = _rs.prepare_variants(
            [_rs.Variant("1", A - 1, contig[A - 1 : A + n], contig[A - 1], "DEL")],
            fa,
            5,
            False,
            1,
            True,
        )
        (c,) = _rs.count_bam_binned(bam, [pv.variant], [None], **ARGS)
        verdicts = Counter(x.value for x in census(bam, contig, pv.variant).counts.elements())
        return (c.rd, c.ad, c.partial_alt), ", ".join(
            f"{k} {m}" for k, m in sorted(verdicts.items())
        )
    run_len, r, a = COMPLEX[case.variant]
    contig = _complex_contig(run_len)
    haps, ends = _complex_reads(run_len, r, a)
    h, e = case.shape.rsplit(", ", 1)
    events, last = haps[h], ends[e]
    reads = []
    for i in range(N_READS):
        s = A - 40 - 3 * i  # same last base, different starts
        gap = sum((len(x) if op == "I" else -x) for _, op, x in events if op in "ID")
        seq, cig = _build(contig, s, last + 1 - s + gap, events)
        reads.append(make_read(f"r{i}", seq, s, cig))
    fa, bam = write_contig(d, contig, reads, "c")
    (pv,) = _rs.prepare_variants([_rs.Variant("1", A, r, a, "X")], fa, 5, False, 1, True)
    (c,) = _rs.count_bam_binned(bam, [pv.variant], [None], **ARGS)
    return (c.rd, c.ad, c.partial_alt), ""


# The decision each group follows: decided (the call is the rule) or open (the
# call is today's, pinned until the decision lands).
DECISIONS = {
    "C10 REF reads": "decided: REF only from informative reads (C10 #157)",
    "invariant 7 carriers": "decided: ALT only where the read's bases hold it (#161, #188, #191, #192, #204)",
    "C1 exact carriers": "decided: complex variants count exact carriers (C1 #141)",
    "C2 siblings": "decided: a sibling's event inside the window withdraws REF (C2 #119)",
    "C26 inside the window": "decided: another indel inside the window, not REF (C26 #200)",
    "C26 outside the window": "decided: another indel outside the window is a separate event (C26 #200)",
    "C27 the ALT across ops": "decided: the ALT across several ops counts ALT, by its bases (C27 #201)",
    "C28 anchor deleted": "decided: a read deleting the anchor is judged by its bases (C28 #202)",
    "C25 long events": "decided: junction windows read on through the read (C25 #199)",
    "C39 one-sided windows": "decided: one pair of equal-length junction windows at one flank (C39 #253)",
    "C40 long deletions": "decided: 50+ bp deletions judged by their bases at the POS breakpoint (C40 #254)",
    "read inputs": "decided: a read contributes its molecule's bases, with qualities (C17 #176, C19 #182; C29 #207 a fix)",
    "C35 unreadable inserts": "decided: an unreadable insert of the ALT's length at its junction is the ALT (RJ-20 #278, superseding C35 #240)",
    "C36 inserts of other bases": "decided: a same-length insert of other bases is judged by its bases, never Phase 3 (C36 #243)",
    "C37 inserts behind a substitution": "decided: a carrier behind a substitution is another allele, not REF (C37 #245); the junction ALT with a substitution stays ALT unless it recurs inside the event: another allele (D4 #278)",
    "C43 other alleles at the locus": "decided: ALT only where the read tells it from every other allele of its length at the locus (C43 #277)",
    "C43 background in the grown event": "decided: a change the REF molecules carry inside the grown event is background (C43 #277)",
}

# Decided calls the engine does not give yet: the red-first cases of a rule being
# built, strict xfails in test_read_judgment_spec.py until it lands.
PENDING: set[str] = set()

# (ref_count, alt_count, partial_alt) for the case's four reads.
EXPECT = {
    "hp-A | REF read spanning the tract": (4, 0, 0),
    "hp-A | REF read ending inside the tract": (0, 0, 0),
    "hp-A | exact carrier": (0, 4, 0),
    "hp-A | carrier with another indel outside the window": (0, 4, 0),
    "hp-A | the ALT at another placement: D1 near the anchor": (0, 4, 0),
    "hp-A | another indel inside the window: I1 near the anchor": (0, 0, 4),
    "hp-A | another indel inside the window: D2 after the anchor": (0, 0, 4),
    "hp-A | another indel outside the window: D1 past the tract": (4, 0, 0),
    "hp-A | another indel outside the window: D5 past the tract": (4, 0, 0),
    "hp-A | another indel outside the window: I1 past the tract": (4, 0, 0),
    "hp-AA | REF read spanning the tract": (4, 0, 0),
    "hp-AA | REF read ending inside the tract": (0, 0, 0),
    "hp-AA | exact carrier": (0, 4, 0),
    "hp-AA | carrier with another indel outside the window": (0, 4, 0),
    "hp-AA | another indel inside the window: D1 near the anchor": (0, 0, 4),
    "hp-AA | another indel inside the window: I1 near the anchor": (0, 0, 4),
    "hp-AA | the ALT at another placement: D2 after the anchor": (0, 4, 0),
    "hp-AA | another indel outside the window: D1 past the tract": (4, 0, 0),
    "hp-AA | another indel outside the window: D5 past the tract": (4, 0, 0),
    "hp-AA | another indel outside the window: I1 past the tract": (4, 0, 0),
    "hp-AA | the ALT written as two deletions": (0, 4, 0),
    "hp+A | REF read spanning the tract": (4, 0, 0),
    "hp+A | REF read ending inside the tract": (0, 0, 0),
    "hp+A | exact carrier": (0, 4, 0),
    "hp+A | carrier with another indel outside the window": (0, 4, 0),
    "hp+A | another indel inside the window: D1 near the anchor": (0, 0, 4),
    "hp+A | the ALT at another placement: I1 near the anchor": (0, 4, 0),
    "hp+A | another indel inside the window: D2 after the anchor": (0, 0, 4),
    "hp+A | another indel outside the window: D1 past the tract": (4, 0, 0),
    "hp+A | another indel outside the window: D5 past the tract": (4, 0, 0),
    "hp+A | another indel outside the window: I1 past the tract": (4, 0, 0),
    "u-8 | REF read spanning the tract": (4, 0, 0),
    "u-8 | REF read ending inside the tract": (4, 0, 0),
    "u-8 | exact carrier": (0, 4, 0),
    "u-8 | carrier with another indel outside the window": (0, 4, 0),
    "u-8 | another indel inside the window: D1 near the anchor": (0, 0, 4),
    "u-8 | another indel inside the window: I1 near the anchor": (0, 0, 4),
    "u-8 | another indel inside the window: D2 after the anchor": (0, 0, 4),
    "u-8 | another indel outside the window: D1 past the tract": (4, 0, 0),
    "u-8 | another indel outside the window: D5 past the tract": (4, 0, 0),
    "u-8 | another indel outside the window: I1 past the tract": (4, 0, 0),
    "u-8 | the ALT written as two deletions": (0, 4, 0),
    "u-8 | anchor deleted, then an insertion": (0, 0, 4),
    "u-8 | anchor deleted": (0, 0, 4),
    "u+10 | anchor deleted": (0, 0, 4),
    "dup+8 | anchor deleted": (0, 0, 4),
    "hp-A | the ALT written as D2 + I1 at the anchor": (0, 4, 0),
    "hp-A | the ALT written as D2 + I1 inside the run": (0, 4, 0),
    "hp+A | the ALT written as I2 + D1 at the anchor": (0, 4, 0),
    "u+10 | REF read spanning the tract": (4, 0, 0),
    "u+10 | REF read ending inside the tract": (4, 0, 0),
    "u+10 | exact carrier": (0, 4, 0),
    "u+10 | carrier with another indel outside the window": (0, 4, 0),
    "u+10 | another indel outside the window: D1 near the anchor": (4, 0, 0),
    "u+10 | another indel outside the window: I1 near the anchor": (4, 0, 0),
    "u+10 | another indel inside the window: D2 after the anchor": (0, 0, 4),
    "u+10 | another indel outside the window: D1 past the tract": (4, 0, 0),
    "u+10 | another indel outside the window: D5 past the tract": (4, 0, 0),
    "u+10 | another indel outside the window: I1 past the tract": (4, 0, 0),
    "ca-CA | REF read spanning the tract": (4, 0, 0),
    "ca-CA | REF read ending inside the tract": (0, 0, 0),
    "ca-CA | exact carrier": (0, 4, 0),
    "ca-CA | carrier with another indel outside the window": (0, 4, 0),
    "ca-CA | another indel inside the window: D1 near the anchor": (0, 0, 4),
    "ca-CA | another indel inside the window: I1 near the anchor": (0, 0, 4),
    "ca-CA | the ALT at another placement: D2 after the anchor": (0, 4, 0),
    "ca-CA | another indel outside the window: D1 past the tract": (4, 0, 0),
    "ca-CA | another indel outside the window: D5 past the tract": (4, 0, 0),
    "ca-CA | another indel outside the window: I1 past the tract": (4, 0, 0),
    "ca-CA | the ALT written as two deletions": (0, 4, 0),
    "dup-8 | REF read spanning the tract": (4, 0, 0),
    "dup-8 | REF read ending inside the tract": (0, 0, 0),
    "dup-8 | exact carrier": (0, 4, 0),
    "dup-8 | carrier with another indel outside the window": (0, 4, 0),
    "dup-8 | another indel inside the window: D1 near the anchor": (0, 0, 4),
    "dup-8 | another indel inside the window: I1 near the anchor": (0, 0, 4),
    "dup-8 | another indel inside the window: D2 after the anchor": (0, 0, 4),
    "dup-8 | another indel outside the window: D1 past the tract": (4, 0, 0),
    "dup-8 | another indel outside the window: D5 past the tract": (4, 0, 0),
    "dup-8 | another indel outside the window: I1 past the tract": (4, 0, 0),
    "dup-8 | the ALT written as two deletions": (0, 4, 0),
    "dup-8 | anchor deleted, then an insertion": (0, 0, 4),
    "dup-8 | anchor deleted": (0, 0, 4),
    "dup+8 | REF read spanning the tract": (4, 0, 0),
    "dup+8 | REF read ending inside the tract": (0, 0, 0),
    "dup+8 | exact carrier": (0, 4, 0),
    "dup+8 | carrier with another indel outside the window": (0, 4, 0),
    "dup+8 | another indel inside the window: D1 near the anchor": (0, 0, 4),
    "dup+8 | another indel inside the window: I1 near the anchor": (0, 0, 4),
    "dup+8 | another indel inside the window: D2 after the anchor": (0, 0, 4),
    "dup+8 | another indel outside the window: D1 past the tract": (4, 0, 0),
    "dup+8 | another indel outside the window: D5 past the tract": (4, 0, 0),
    "dup+8 | another indel outside the window: I1 past the tract": (4, 0, 0),
    "C>TA run10 | REF, ends inside the run": (0, 0, 0),
    "C>TA run10 | REF, ends on the base after the run": (0, 0, 0),
    "C>TA run10 | REF, ends past the run": (4, 0, 0),
    "C>TA run10 | substitution only, ends inside the run": (0, 0, 0),
    "C>TA run10 | substitution only, ends on the base after the run": (0, 0, 0),
    "C>TA run10 | substitution only, ends past the run": (0, 0, 0),
    "C>TA run10 | anchor kept, one A more, ends inside the run": (0, 0, 0),
    "C>TA run10 | anchor kept, one A more, ends on the base after the run": (0, 0, 0),
    # one base from each allele (an A more than REF, C for ALT's T): neither (RJ-23)
    "C>TA run10 | anchor kept, one A more, ends past the run": (0, 0, 4),
    "C>TA run10 | anchor kept, one A fewer, ends inside the run": (0, 0, 0),
    "C>TA run10 | anchor kept, one A fewer, ends on the base after the run": (0, 0, 0),
    "C>TA run10 | anchor kept, one A fewer, ends past the run": (0, 0, 4),
    "C>TA run10 | exact carrier, ends inside the run": (0, 0, 0),
    # holds the REF window's 15 bases: the T, all eleven A's and the G after (RJ-23)
    "C>TA run10 | exact carrier, ends on the base after the run": (0, 4, 0),
    "C>TA run10 | exact carrier, ends past the run": (0, 4, 0),
    "C>TA run60 | REF, ends inside the run": (4, 0, 0),
    "C>TA run60 | REF, ends on the base after the run": (4, 0, 0),
    "C>TA run60 | REF, ends past the run": (4, 0, 0),
    "C>TA run60 | substitution only, ends inside the run": (0, 4, 0),
    "C>TA run60 | substitution only, ends on the base after the run": (0, 0, 0),
    "C>TA run60 | substitution only, ends past the run": (0, 0, 0),
    "C>TA run60 | anchor kept, one A more, ends inside the run": (4, 0, 0),
    "C>TA run60 | anchor kept, one A more, ends on the base after the run": (0, 0, 0),
    "C>TA run60 | anchor kept, one A more, ends past the run": (0, 0, 4),
    "C>TA run60 | anchor kept, one A fewer, ends inside the run": (4, 0, 0),
    "C>TA run60 | anchor kept, one A fewer, ends on the base after the run": (0, 0, 0),
    "C>TA run60 | anchor kept, one A fewer, ends past the run": (0, 0, 4),
    "C>TA run60 | exact carrier, ends inside the run": (0, 4, 0),
    "C>TA run60 | exact carrier, ends on the base after the run": (0, 4, 0),
    "C>TA run60 | exact carrier, ends past the run": (0, 4, 0),
    "CA>T run60 | REF, ends inside the run": (4, 0, 0),
    "CA>T run60 | REF, ends on the base after the run": (4, 0, 0),
    "CA>T run60 | REF, ends past the run": (4, 0, 0),
    "CA>T run60 | substitution only, ends inside the run": (0, 4, 0),
    "CA>T run60 | substitution only, ends on the base after the run": (0, 0, 0),
    "CA>T run60 | substitution only, ends past the run": (0, 0, 0),
    "CA>T run60 | anchor kept, one A more, ends inside the run": (4, 0, 0),
    "CA>T run60 | anchor kept, one A more, ends on the base after the run": (0, 0, 0),
    "CA>T run60 | anchor kept, one A more, ends past the run": (0, 0, 4),
    "CA>T run60 | anchor kept, one A fewer, ends inside the run": (4, 0, 0),
    "CA>T run60 | anchor kept, one A fewer, ends on the base after the run": (0, 0, 0),
    "CA>T run60 | anchor kept, one A fewer, ends past the run": (0, 0, 4),
    "CA>T run60 | exact carrier, ends inside the run": (0, 4, 0),
    "CA>T run60 | exact carrier, ends on the base after the run": (0, 4, 0),
    "CA>T run60 | exact carrier, ends past the run": (0, 4, 0),
    # C39 #253: one-sided, equal-length junction windows (RJ-23).
    "delins 28>6 | exact carrier, eight flank bases right": (0, 4, 0),
    "delins 28>6 | exact carrier, four flank bases left": (0, 4, 0),
    "delins 28>6 | exact carrier, long flanks": (0, 4, 0),
    "delins 28>6 | REF read, four flank bases left": (4, 0, 0),
    "delins 28>6 | REF read, four flank bases right": (4, 0, 0),
    "delins 28>6 | REF read starting inside the event": (0, 0, 0),
    # the change is a separate event, masked: the carriers past it count ALT (RJ-23)
    "delins 28>6 | carriers ending before a change that recurs on ALT reads only": (0, 4, 0),
    "delins 60>2 | exact carrier, eight flank bases right": (0, 4, 0),
    "delins 60>2 | exact carrier, long flanks": (0, 4, 0),
    "delins 60>2 | REF read, four flank bases left": (4, 0, 0),
    "delins 60>2 | REF read starting inside the event": (0, 0, 0),
    # RJ-23 amendments (2026-10-12)
    "delins 17>1 | exact carrier, long flanks": (0, 4, 0),
    "delins 17>1 | REF read, four flank bases right": (4, 0, 0),
    # the windows run past the bases the alleles share after the event: a read with
    # the event's first base substituted, then REF's bases, is neither allele
    "delins 17>1 | REF read with the event's first base substituted": (0, 0, 0),
    # an isolated recurring change is a separate event, masked; a run of changes
    # (the base past the window included) is another allele, closer to ALT
    "delins 28>6 | carriers with a recurring change on the window's last base": (0, 4, 0),
    "delins 28>6 | carriers with a run of changes from the window's last base on": (0, 0, 4),
    # placed from the other flank, a read's flank must favour that placement over its
    # own alignment by one read base; a carrier the aligner wrote at the event's
    # other copy is placed and counts though it starts past POS
    "delins 28>6 | REF read starting inside the event behind a low-quality clip": (0, 0, 0),
    "delins 28>6 | carrier written at the event's other copy as plain matches": (0, 4, 0),
    "u-8 | sibling SNV inside the span": (0, 0, 4),
    "u-8 | sibling SNV outside the window": (4, 0, 0),
    "read inputs | absent base qualities (QUAL '*'), showing the ALT": (0, 0, 0),
    "read inputs | previous complex classifier: REF read with a leading hard clip": (4, 0, 0),
    "read inputs | previous complex classifier: REF read, no clip": (4, 0, 0),
    "read inputs | read-through: the base on the SNV is adapter, showing REF": (0, 0, 0),
    "read inputs | read-through: the base on the SNV is adapter, showing the ALT": (0, 0, 0),
    "read inputs | read-through: the read aligns on past its TLEN end (R2's 5' end clipped), showing the ALT": (
        0,
        4,
        0,
    ),
    # The read's own inserted bases unreadable (N at Q2, or letters at Q5): the
    # insert's length at the ALT's junction is the ALT (RJ-20, superseding C35 #240).
    "hp+A | an unreadable insert at the anchor (N)": (0, 4, 0),
    "hp+A | an unreadable insert inside the run (N)": (0, 4, 0),
    "hp+A | an insert below min BQ inside the run": (0, 4, 0),
    "u+10 | an unreadable insert (N)": (0, 4, 0),
    "u+10 | an insert below min BQ": (0, 4, 0),
    "u+10 | a partly masked insert": (0, 4, 0),
    "u+10 | anchor deleted, the insert unreadable": (0, 4, 0),
    # C36 #243: same-length inserts of other bases, judged by their bases.
    "u+10 | a same-length insert of other bases one junction left": (4, 0, 0),
    "u+10 | a same-length insert of other bases four junctions right": (4, 0, 0),
    "u+10 | the ALT written one junction off with a compensating mismatch": (0, 4, 0),
    "u+10 | the ALT with an error at its first inserted base, written one junction right": (
        0,
        0,
        4,
    ),
    "dup+8 | a same-length insert of other bases inside the tract": (0, 0, 4),
    # C37 #245: a readable insert whose fewest-mismatch placement sits inside the
    # window is another allele (RJ-22); the junction ALT with a substitution in the
    # run that all four reads carry is another allele too (D4: it recurs in the event).
    "hp+A | the ALT at the junction with a substitution in the run": (0, 0, 4),
    "hp+A | the ALT with its anchor substituted, written one junction left": (0, 0, 4),
    "u+10 | the ALT with its anchor substituted, written one junction left": (0, 0, 4),
    "del 60 | exact D carrier, long flanks": (0, 4, 0),
    "del 60 | exact D carrier, four bases past the microhomology": (0, 4, 0),
    # three junction bases past the microhomology: too few to judge (RJ-24)
    "del 60 | exact D carrier, three bases past the microhomology": (0, 0, 0),
    "del 60 | right-clipped carrier, the junction in the clip": (0, 4, 0),
    "del 60 | left-clipped carrier aligned from the right flank": (0, 4, 0),
    # placed from the right, a read's flank must favour the junction over its own
    # alignment by one read base: the base next to the junction suffices; a flank
    # read too poorly to show it is no junction
    "del 60 | left-clipped carrier with two low-quality flank bases": (0, 4, 0),
    "del 60 | left-clipped carrier with a low-quality flank": (0, 0, 0),
    "del 60 | split carrier: a soft-clipped primary and a hard-clipped supplementary": (0, 4, 0),
    "del 60 | hard-clipped primary aligned from the right flank": (0, 0, 0),
    "del 60 | REF read, long flanks": (4, 0, 0),
    "del 60 | REF read, four deleted bases past the microhomology": (4, 0, 0),
    # REF needs as many deleted bases as ALT needs right-flank bases (RJ-24)
    "del 60 | REF read, three deleted bases past the microhomology": (0, 0, 0),
    # nothing past the window is checked, for either allele (RJ-24)
    "del 60 | REF read with an error 20 bases past the window": (4, 0, 0),
    "del 60 | REF read with an error inside the window": (3, 0, 0),
    # another allele, one mismatch from each window more than the other: neither
    "del 60 | a deletion two bases shorter at the same start": (0, 0, 0),
    # another allele, closer to ALT than REF over the windows: partial
    "del 60 | the deletion written as two with two bases kept between": (0, 0, 4),
    # its own deletion inside the deleted span: another allele (C26), partial
    "del 60 | REF read with its own 20-bp deletion inside the deleted span": (0, 0, 4),
    "del 60 | REF read starting inside the deletion": (0, 0, 0),
    # placed from the right, a read must show the flank: a clip read too poorly to
    # show it is no junction (REF molecules there show the same right flank)
    "del 60 | REF read starting inside the deletion behind a low-quality clip": (0, 0, 0),
    "del 60 | REF read starting inside the deletion behind an N clip": (0, 0, 0),
    # a recurring change is a separate event: masked, the carriers count
    "del 60 | carriers with a recurring change two bases past the microhomology": (0, 4, 0),
    "del 60 | REF reads with a recurring change inside the window": (4, 0, 0),
    "del 60 | carriers with a recurring change in the flank": (0, 4, 0),
    "del 60 | right-clipped carriers with a recurring change in the flank": (0, 4, 0),
    "del 60 | carriers with a recurring change on the window's last base": (0, 4, 0),
    # a run of changes is another allele (here nearer ALT): partial, whether it lies
    # inside the window or runs on past its edge
    "del 60 | carriers with a run of changes in the flank": (0, 0, 4),
    "del 60 | carriers with a run of changes from the window's last base on": (0, 0, 4),
    # a run the REF reads carry too is background (a germline MNP), masked as a
    # germline SNP is: both alleles count
    "del 60 | REF reads with a run of changes in the flank": (4, 0, 0),
    "del 60 | carriers and REF reads with a run of changes in the flank": (4, 4, 0),
    # a run past the junction, where the alleles differ, is another allele even on
    # reads that fit REF elsewhere (a germline run is background only in the bases
    # both alleles share): neither
    "del 60 | REF reads with a run of changes from the window's last base on": (0, 0, 0),
    # placed from the right flank, its bases before the junction show the left flank,
    # not the deleted bases its alignment puts there: ALT, as a REF read from the
    # same start is REF
    "del 60 | carrier aligned over the deleted bases' end, no clip": (0, 4, 0),
    "del 60 mh3 | exact D carrier, long flanks": (0, 4, 0),
    "del 60 mh3 | exact D carrier, four bases past the microhomology": (0, 4, 0),
    "del 60 mh3 | exact D carrier, three bases past the microhomology": (0, 0, 0),
    "del 60 mh3 | left-clipped carrier aligned from the right flank": (0, 4, 0),
    "del 60 mh3 | REF read, long flanks": (4, 0, 0),
    "del 60 mh3 | REF read, four deleted bases past the microhomology": (4, 0, 0),
    "del 60 mh3 | REF read, three deleted bases past the microhomology": (0, 0, 0),
    "del 1000 | exact D carrier, long flanks": (0, 4, 0),
    "del 1000 | right-clipped carrier, the junction in the clip": (0, 4, 0),
    "del 1000 | left-clipped carrier aligned from the right flank": (0, 4, 0),
    "del 1000 | REF read, long flanks": (4, 0, 0),
    "del 1000 | REF read starting inside the deletion": (0, 0, 0),
    # Other alleles at the locus (RJ-25): each count holds the shape's four reads and
    # the locus's others (another allele's reads are another allele).
    "delins 28>6, another allele | carriers with the telling base below min BQ, the other allele on four fragments": (
        0,
        0,
        8,
    ),
    "delins 28>6, another allele | carriers reading the telling base, the other allele on four fragments": (
        0,
        4,
        4,
    ),
    "delins 28>6, another allele | carriers with the telling base below min BQ, the other allele on two fragments": (
        0,
        4,
        2,
    ),
    "delins 28>6, another allele | carriers with the telling base below min BQ, no other allele": (
        0,
        4,
        0,
    ),
    "delins 28>6, another allele | carriers with the telling base below min BQ, the other allele co-annotated": (
        0,
        0,
        4,
    ),
    "MNP GC>AT, germline GT | reads showing only the germline base, the germline allele on four fragments": (
        0,
        0,
        8,
    ),
    "MNP GC>AT, germline GT | reads showing only the somatic base, the germline allele on four fragments": (
        0,
        4,
        4,
    ),
    "MNP GC>AT, germline GT | reads showing only the germline base, no germline allele": (0, 4, 0),
    # Background in the grown event (RJ-23, 2026-10-14): REF reads beside the shape's.
    "G>TCA in a CA run, germline SNP | REF reads with the SNP, beside REF reads without it": (
        8,
        0,
        0,
    ),
    "G>TCA in a CA run, germline SNP | carriers with the SNP, beside REF reads with it": (4, 4, 0),
    "G>TCA in a CA run, germline SNP | carriers with the SNP on the carriers alone": (4, 0, 4),
}


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmp:
        for case in cases():
            got, verdicts = run(case, Path(tmp))
            print(f"{case.group}\t{case.id}\t{got}\t{verdicts}", flush=True)
