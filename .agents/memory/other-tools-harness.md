---
name: other-tools-harness
description: Running GATK, bcftools, freebayes, VarDict, Paragraph, GraphTyper, Delly, SVTyper on the harness slices — the setup and the traps
metadata:
  type: reference
---

`~/test/gbcms/harness/plan670/tools/` (local; PHI for panel, WES, normals) runs ten
other tools in amd64 containers on this Mac and compares them with gbcms and with an
independent per-read check (`readtruth.py`: BWA-MEM of every read to the REF and ALT
local haplotypes). Entry points: `make_sites.py`, `run_tools.py`, `compare.py`,
`summarize_tools.py`, `per_read.py`, `build_report.py` (2026-10-12 audit report).

Traps found the first time:
- Run **one container at a time**: five emulated containers at once hung Docker's file
  share (every container at 0% CPU; `docker desktop restart` clears it). Mount only the
  harness folder; batch per-BAM tools in one container (a start costs ~15 s).
- The **IMPACT and ACCESS slices carry no proper-pair flag**: bcftools skips such reads
  by default and sees nothing (use `-A`; gbcms does not filter on it either).
- GATK 4.6.2.0 **Mutect2 crashes** on forced alleles longer than its region padding
  (`--assembly-region-padding` must cover the longest deletion; slow under emulation).
- **Delly 2.x genotyping** reads sites only in its own format (`SVMETHOD=EMBL.DELLY*`,
  `CONSENSUS`/`CONSBP`, and `SR`/`PE` > 0 or it writes nothing) and needs every BAM of a
  run to share contigs. Delly and SVTyper estimate insert sizes from the BAM and fail on
  small slices. Paragraph 2.3's VCF writer aborts (htslib assertion): read
  `genotypes.json.gz`. GraphTyper copies output with permissions the share refuses:
  write to `/tmp` in the container. Octopus 0.7.4 hangs under emulation (AVX).
- GATK, freebayes and VarDict **decompose a given delins**: credit them with the
  records that together spell the allele (`compare.from_components`).
- The per-read check reads deletions' VAF a little low: a REF read starting just
  before POS holds deleted bases and is clearly REF; an ALT read starting just before
  the junction is mostly shared right flank and stays unclear.

Related: [[survey-several-tools]], [[test-across-data-types]], [[pysam-validation-oracle]].
