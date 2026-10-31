# Spike-in truth for #297: a third unit count at long tandem-repeat deletions (D12 phase 1b)

Spec for the spike-in arm (D12, #275). It defines the loci, the alleles, the variant allele
fractions (VAF) and what to measure, so that the operator can decide #297's narrow rule
(RJ-24, amended 2026-10-30) on truth. HG002 cannot show the benefit: it is diploid, so no
HG002 site holds REF, the given deletion and another unit count together.

## What the rule does (short)

A deletion of 50 bp or more is judged at its junction (RJ-24). Where the bases after the
deletion are a tandem run of its unit and the deleted copy differs from the run's copies, the
junction windows stop at the copy's first differing base. Another unit count (one more unit
deleted, or one unit inserted) then reads like the given deletion until the run's copies end.
Call the distance from the windows' first base to that first difference **D**.

- Where the other count is present at the site (listed in the input, or spelled in full on 3
  or more fragments and 5% of the fragments that spell the given deletion there), an ALT read
  must reach D and show the given count there. A read that ends sooner is "either allele":
  partial evidence, not ALT.
- Where no other count is present, nothing changes.
- REF reads keep their windows (REF differs from both deletions at the junction).
- Where D is longer than the longest read at the site, nothing is held and the row is named
  in `gbcms_diagnostic`: `UNIT_COUNT_UNSPANNED(change:bases/longest)`.

## Loci

Select from the reference only, with the arm's own code (not gbcms code), in the territory
that phase 1 uses (inside the GIAB v5.0q benchmark, depth 20+, no HG002 truth variant within
250 bp of the run).

- A tandem run whose unit gives a deletion of 50 bp or more: one unit of 50 to 150 bp, or
  several units of a shorter unit (7 to 49 bp) that together make 50 bp or more.
- For each locus, compute from the reference: **mh**, how far the given deletion slides
  right (its microhomology); and **D**, the first base where the haplotype "given deletion"
  and the haplotype "one more unit deleted" differ, counted from 4 bases before the junction.
- Classes (40 loci each, 20 at least):

| Class | Condition | What it tests |
|:--|:--|:--|
| T1 | D − mh from 12 to 60 | the hold, short |
| T2 | D − mh from 61 to 120 | the hold, long |
| T3 | D from 121 to 145 (2×150 reads) | the hold near the read length |
| U | D above 150 | no hold; the row is named (`UNIT_COUNT_UNSPANNED`) |
| C | identical copies (the deletion slides through the run: D − mh below 12) | control: the windows already hold the run; nothing changes |

The 772 identical-unit loci in phase 1's notes are class C. Classes T1 to T3 need copies that
differ (most real VNTRs): the deleted copy and the next copy differ within the first 10 bases
of the junction (small mh), while the run's later copies match until D.

## Alleles and VAFs

Three alleles per locus: REF, the **given** deletion (the row in the input), and the
**other** count. Two kinds, half the loci each: the given deletion removes one unit (the copy
that differs) and the other removes one more (**more**); or the given removes two units and
the other one (**fewer**). An insertion of a unit keeps the differing copy, so the junction
tells it from the given deletion: use it only as a control (it must change nothing). Molecule
VAFs (each locus gets one setting; split the loci of each class over the settings):

| Setting | Given | Other | REF | What it tests |
|:--|--:|--:|--:|:--|
| S1 | 0.50 | 0 | 0.50 | control: no other count; counts equal develop's |
| S2 | 0.25 | 0.25 | 0.50 | the core case: both present beside REF |
| S3 | 0.40 | 0.10 | 0.50 | the other count near the 3-fragment floor |
| S4 | 0 | 0.40 | 0.60 | the given allele absent: any AD is a false ALT (develop counts the other count's short reads ALT) |
| S5 | 0.50 | 0.50 | 0 | two alleles, no REF (a compound het) |
| S6 | 0.20 | 0.05 | 0.75 | the other count below the floor: the rule must not fire |

Plus the VAF-0 arm (nothing spiked), as in phase 1.

Two input forms for every row: the other count **not listed** (it must be found from the
reads) and **listed** (a second input row). At HG002's depth (about 35 molecules) few reads
reach D, so the found path rarely fires; the listed form tests the hold there. The found path
needs depth: run S2 to S4 and S6 again on a deeper background (a deep public HG002 set, or a
patient background on HPC only, phase B).

## Runs

- Alignment: BWA-MEM only, and ABRA2 (phase 1's re-alignment).
- Read length: 2×150. If the arm can trim, also 2×100 (it moves loci between T3 and U, and so
  checks the read-length limit).
- gbcms: develop and the #297 branch, the rules and the sequence judge (`--read-judge
  sequence`), `--trace`.

## Truth and scores

Truth per read, by its bases against three haplotypes (REF, given, other), with the arm's DP
(method A for development, B for the gate): **given** where the read fits the given haplotype
better than both others by 0.9 of an edit; **other** likewise; **either** where it fits the
given and the other alike and REF worse (it ends before D); **REF**; **none**.

Per class, setting, input form, alignment and path:

1. False ALT: reads whose truth is "other" or "either" that count ALT; AD in S4 (must be 0
   where the other count is listed or found).
2. Lost ALT: reads whose truth is "given" that do not count ALT.
3. Partial: reads counted as partial evidence, by truth.
4. VAF: AD / (RD + AD) against the exact molecule fraction given / (given + REF). REF keeps its
   windows while an ALT read must reach D, so beside a present other count the VAF reads low by
   the share of ALT reads that end before D (RJ-24's named limit; the synthetic contract test
   `test_ad_reads_low_beside_a_present_count` pins 0.379 for a true 0.5 at D = 64 with 150-bp
   reads). **Bar:** by class and setting (S2, S3, S5 listed and found), the measured VAF lies
   within 0.05 of the limit's prediction (the exact molecule VAF times the share of the given
   allele's reads that reach D, over that share plus the REF reads holding the windows); a
   larger gap is a defect, not the limit. The operator decides whether the predicted bias
   itself is acceptable, by class.
5. RD against develop: must be equal (the rule changes no REF call).
6. The diagnostics: `UNIT_COUNT_UNSPANNED` on every class-U row and on no row of T1 to T3
   (2×150); the other count named (`OBSERVED_ALLELE` or `COEXISTING_ALLELE` at the deletion's
   anchor) on every row where it is present (listed, or found) and on no S1 row, nor on an S6
   row of the found form.
7. S1 and S6 rows, and class C: counts equal to develop's.

## What this cannot see

UMI consensus (ACCESS), cfDNA fragment sizes, ABRA2 at panel depth and RNA. These need patient
backgrounds on HPC (phase B). Report them as not measured.
