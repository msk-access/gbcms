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
