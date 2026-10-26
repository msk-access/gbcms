---
name: giab-truth-arm
description: The HG002 truth arm (D9) — public data, how sites are chosen, how to fetch region slices fast, and the pitfalls met building it.
metadata:
  type: reference
---

Harness: `~/test/gbcms/harness/plan670/giab/` (public data; numbers can go on GitHub).
- **Truth:** HG002 v5.0q smvar and stvar (GRCh38, NIST FTP `release/AshkenazimTrio/
  HG002_NA24385_son/v5.0q/`), single-ALT phased hets inside the benchmark BED, no
  other truth variant within 30 bp. dipcall splits a delins into nearby records, so
  same-phase records within 10 bp (one an indel) are merged into one delins.
- **Reads:** the DeepVariant case-study `HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam`
  (bwa-mem, no realignment, 151 bp), fetched as region slices with `samtools view -M -L`.
  The local Ensembl-named `~/Downloads/genome.fa` works: gbcms reconciles `chr1` ↔ `1`.
- **Fetching is latency-bound** from this Mac (about 1.5 s to the first byte, ~20 KB/s
  per connection): one connection did ~15 regions a minute; 48 parallel chunks fetched
  700 sites in minutes. Too many connections at once can fail ("Destination address
  required"); retry the failed set.
- **Run without `--trace`** for truth counts (a traced run over 270k reads took 20+
  minutes); the bases-only oracle does the per-read work.

Related: [[decisions-worked-examples]], [[pysam-validation-oracle]].

**Flaw found and fixed 2026-10-26:** `r2/truth.py` maps each read against all rows'
haplotype contigs in one combined BWA index. Rows 1 bp apart (#246 pairs) have identical
contigs, so a read's best hit can land on the neighbour's contig, leaving no or only a partial
alignment to its own row; a missing alignment scored 0 and could decide the label (~9% of
het-row labels on #246, ~42% on GIAB TR rested on one haplotype only). Truth v2: per-row
alignment, end to end, two methods independent of the judge, missing = unphased. Until v2
lands, treat old-truth numbers as build-to-build comparisons, not absolutes.
**Truth v2 (done 2026-10-26):** `r2/truth_v2.py` (+ `truthv2_lib.py`, `truthv2_dp.c`), files
`truth_{reads,rows}_v2.tsv` beside the old; `evaluate.py --truth v2`. Scoring rules (operator):
bars on reads decidable at the run's `--min-baseq` (default 20: `v2q20`; full v2 reported as the gate's cost; the threshold comes from the run, never hard-coded); long deletions
breakpoint reads only (strict v2; `v2dp` beside); a class fails only past 0.5 points AND 2
reads. All earlier decisions held on v2.

**Truth v2c supersedes v2 (2026-10-26):** v2 scored BWA with its `AS` tag, which is the local
best and leaves a read-end mismatch uncharged (24% of v2 scores sat above the true end-to-end
optimum). That under-charged exactly the read-end bases that clip admission rests on. v2c
rescores every alignment over the whole read from its CIGAR (`r2/bwa_rescore.py`, shared with
the patient oracle `xoracle_v2.py`). Files: `truth_{reads,rows}_v2c.tsv` plus
`truth_settings_v2c.tsv` (the runs' `--min-baseq`). Score with `--truth v2c` / `v2cq20`; `v2cq20`
refuses a run whose threshold differs. Result: about 4,200 more decided reads (half with X
within 4 bases of the read end), no read moved between REF and ALT, no verdict changed. The
scratchpad copies of these tools are mirrored in `plan670/sp_tools/` (scripts only).
**Truth v2d and the scoring rules (2026-10-26; operator):** `r2/truth_v2.py --tag v2d` writes
`truth_*_v2d.tsv` (v2c is untouched). Score with `--truth v2d` / `v2dq20`.
- **Homozygous rows** are compared with REF. Under v2c's shortcut chr18:12145098's reads were
  labelled ALT; they are REF.
- **Bars count rows inside the GIAB v5.0q benchmark,** picked by variant size: the smvar BED
  below 50 bp, the stvar BED at 50 bp or more. The smvar BED excludes every 50 bp+ variant by
  construction.
- **Long deletions:** count reads holding the anchor, plus reads placed inside the deleted bases
  whose own bases hold the junction.
- **Reads only phasing can decide leave the bars.** A read whose bases fit its truth sequence and
  another candidate equally is not "decidable by bases".

Use v2d for every HG002 bar from now on.
