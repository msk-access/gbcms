---
name: mapq0-loci-pms2
description: Pseudogene loci like PMS2 are run at --min-mapq 0; MAPQ-0 alignments must stay countable, and read-admission changes are validated at --min-mapq 0 too.
metadata:
  type: project
---

Some genes (PMS2, with its pseudogene PMS2CL) are genotyped with `--min-mapq 0`
because their reads multi-map and come out MAPQ 0 (operator, 2026-10-02). gbcms
must keep mapped MAPQ-0 alignments countable there (Mutect2, by contrast, drops
them by default: MappingQualityNotZeroReadFilter). Reads whose mate is unmapped
also count (gbcms default; GATK the same; samtools/bcftools mpileup skip them as
anomalous pairs unless -A).

**Why:** at `--min-mapq 0` the MAPQ filter no longer screens MAPQ-0 junk: an
unmapped mate (flag 0x4) carrying a CIGAR reached DP and ALT in a synthetic test
(O7 #183), where at the default MAPQ it only inflated mq0_count.

**How to apply:** any change to read admission or read inputs (filters, clips,
fragment rules) is measured at the default MAPQ and at `--min-mapq 0`, including
PMS2 rows when the data has them. Related: [[holistic-effects-map]].
