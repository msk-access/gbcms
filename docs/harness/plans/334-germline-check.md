# #334 — the germline het and hom check of MNVs and complex events (HG002 exome)

Status: measured on public data (2026-11-23); the patient part waits on the operator's sample-to-BAM mapping file.
Issue: #334 (J3, acceptance at the default: every data type and input, two truth sources). Parent: #331 (the
sequence read judge as the default).

## Why

The operator, 2026-11-23: "the point to switch to judge was the hope to be better at mnvs and complex events", and
"I think we should also have germline events check / Het and Homzygous". The bar is the operator's policy: "no
count-chaning it should be about regression, we want to be closer to the truth with less FP and FN increasing
accuracy". A path passes a class only when it lowers both missed ALT and false ALT on rows of known truth.

## What was run

- **Truth:** GIAB HG002 v5.0q small variants, GRCh37, phased by assembly (MD5s equal GIAB's `checksum.md5`).
- **Reads:** the GIAB HG002 exome (Agilent SureSelect v5, bwa-mem), in the exome territory (`runs/hgd333/cover`).
- **Paths:** develop 6ffde565 (build dv6ff), `gbcms dna` defaults (`--min-baseq 20`, `--min-mapq 20`), with
  `--read-judge rules` (dvR) and `--read-judge sequence` (dvJ).
- **Rows** (`kit/hgd334g_rows.py`; label hgd334i, 34,931 rows):
  - one event per haplotype: touching SNVs join into an MNV; an indel and another event within 3 bases join into a
    DELINS;
  - zygosity: hom (both haplotypes give the same sequence over the event's ambiguity span), het (the other
    haplotype is REF there), multi (the other haplotype changes a base there);
  - context: alone (no other event within 30 bp), cis (another event on the same haplotype), trans (only on the
    other haplotype);
  - events within 30 bp of a record with an unknown haplotype make no rows;
  - truth-ALT layers: gsnv and gind (capped 3,000 per stratum), every gmnv (915) and gdelins (1,402) row;
  - false-ALT layers: gtrans (91: two het events on opposite haplotypes written as one allele) and gwrong (3,466:
    beside each MNV and DELINS row, one ALT base set to a third base (base), and, for a DELINS, one more REF base
    deleted (del+1)). A false allele that either haplotype carries is dropped (128).
- **Scorer:** `kit/local/score_hgd334g.py` (selftest 0 failures). Depth 10+ by the rules path. Truth-ALT: missed
  (AD 0), false REF at hom and multi rows (RD 2+), and, row by row, closer or farther from the expected AD/DP.
  Truth-REF: false ALT (AD 1+ and 2+). The four count relations hold on every row of both paths.

## Results (hgd334i, judge against rules)

Truth-ALT rows (31,374 at depth 10+):

| Layer | rows | missed, rules | missed, judge | false REF 2+, rules | false REF 2+, judge | judge closer / farther |
|:--|--:|--:|--:|--:|--:|--:|
| DELINS het | 440 | 34 | 18 | . | . | 191 / 74 |
| DELINS hom | 204 | 10 | 5 | 6 | 13 | 143 / 5 |
| DELINS multi | 758 | 87 | 50 | 229 | 198 | 254 / 144 |
| MNV het | 483 | 32 | 29 | . | . | 104 / 32 |
| MNV hom | 340 | 2 | 1 | 17 | 8 | 24 / 4 |
| MNV multi | 92 | 1 | 1 | 1 | 1 | 16 / 10 |
| indel het | 7,466 | 255 | 167 | . | . | 1,953 / 462 |
| indel hom | 6,349 | 97 | 62 | 1,093 | 822 | 2,413 / 146 |
| SNV het | 7,335 | 99 | 102 | . | . | 125 / 107 |
| all | 31,374 | 774 | 561 | 2,271 | 1,838 | 5,897 / 1,340 |

False-ALT rows (3,557 at depth 10+); the clean subset is the partners with no other event within 30 bp:

| False allele | rows | ALT 1+, rules | ALT 1+, judge | ALT 2+, rules | ALT 2+, judge |
|:--|--:|--:|--:|--:|--:|
| DELINS, third base, clean | 420 | 54 | 48 | 34 | 30 |
| DELINS, one more base deleted, clean | 418 | 101 | 120 | 69 | 90 |
| MNV, third base, clean | 521 | 13 | 13 | 6 | 6 |
| DELINS, beside other events | 1,713 | 418 | 440 | 243 | 272 |
| two trans events as one DELINS | 51 | 17 | 20 | 6 | 17 |
| two trans events as one MNV | 40 | 1 | 1 | 0 | 0 |

What it says:

1. **Missed ALT and false REF: the judge is better in every MNV and DELINS stratum.** It misses 213 fewer rows in
   all, and 433 fewer rows carry 2+ false REF reads. The one stratum that is worse on one measure is DELINS hom
   (false REF 2+ on 13 rows against 6, with 5 fewer missed).
2. **False ALT: the judge is worse on complex alleles one deleted base from the truth.** On the clean del+1 rows it
   credits 2+ ALT reads on 90 of 418 against the rules' 69. Example, 8:120818709: the truth is `CTTTTTT>CTTTC`; the
   false allele deletes one more T, so it is the truth with 12 T where the truth has 13 after the change. The judge
   credits 15 reads against 1. (Corrected by the trace, §"The del+1 trace": 14 of the 15 end inside the T run, so
   they fit the truth and the false allele alike; the reads that carry the truth and span the run are another allele
   under the judge.)
3. **Both paths credit that false allele often** (24% of the clean del+1 rows at 1+). The causes are in §"The del+1
   trace".
4. **Rows that both paths miss are mapping, not counting.** Most of the 26 het MNVs in cis with no ALT read lie in
   paralogous regions (LPA at 6:161.05 Mb; the chr19 clusters at 43.9, 50.47 and 54.7 Mb), with RD close to DP.

By the policy, the judge does not yet pass for complex events: it lowers missed ALT but raises false ALT on
near-miss complex alleles.

## The del+1 trace (hgd334d, 2026-11-24)

- **Method.**
  - Job hgd334d_slice (13573795) cut the HG002 exome's reads at the 418 clean del+1 rows (each ±1,000 bases,
    385,308 reads). A frozen build of develop 6ffde565 on the operator's Mac (`~/test/gbcms/harness/d334/`) gives the
    hgd334i arms' counts on all 418 rows, both paths.
  - A pat334-shaped tree holds three sets, each run alone with `--trace`: the truth (call), the base partner and the
    del+1 partner. `--trace` changes no count.
  - The patient oracle (`kit/pat334_oracle.py`: BWA per haplotype and the exact DP; a label only where both agree)
    judged each traced read by its own bases. The truth is a sibling haplotype there (SIB).
- **Run alone, as the patient arms run each set, the false ALT is larger.** In hgd334i the base partner was listed
  beside the del+1 row, and a listed sibling takes some reads.

  | Set (418 rows, depth 10+) | ALT 2+, rules | ALT 2+, judge |
  |:--|--:|--:|
  | del+1 partner, alone | 115 | 108 |
  | base partner, alone | 39 | 84 |
  | truth (call), alone | 399 | 411 |
  | del+1 with the base partner listed (hgd334i) | 69 | 90 |

- **The reads each path credits as ALT at the del+1 rows (alone), by the oracle's label:**

  | Path | ALT calls | carry the truth (SIB) | fit the truth and the ALT alike (either) | carry the ALT | other or unk |
  |:--|--:|--:|--:|--:|--:|
  | rules | 1,032 | 495 | 270 | 140 | 127 |
  | judge | 1,189 | 291 | 619 | 154 | 125 |

  The reads that carry the ALT lie almost all where a run of 10 bases or more follows the change: slips of the
  truth's run. Their bases carry the ALT, so they count (§"Oracles and truths").
- **The truth and the del+1 allele differ only by the length of the run that follows the change** (one base fewer),
  so the difference shows at the run's far end: at 8:120818709, 13 T against 12, 13 bases past the change.
- **Three mechanisms:**
  1. The judge credits reads that end inside that run (the 619 "either" reads; at 8:120818709, 14 of its 15 ALT
     reads). The rules path calls them uninformative there.
  2. The judge credits reads that carry the truth, one edit from the ALT and two or more from REF (the 291 SIB reads;
     at 1:154306330 `CAAACT>CAAAAA`, 15 reads where no run follows). The margin decides (`judge/decide.rs`, the module
     header: ALT where the best ALT haplotype beats the rest by one confident base), and the truth is not listed.
  3. The rules path's exact-carrier rule credits reads whose bases past its window contradict the ALT (453 SIB reads;
     at 6:108243120 `AAAAACA>GGGA`, 107 reads with 13 A after `GGG` where the ALT has 12). Its window is the event
     and one base each side: the false allele does not change the run's length against REF, so the window is not
     grown through it (`carrier.rs`, `grown_event`).
- **What the decided rules say.** RJ-4: a complex event counts exact carriers, the read's bases across the whole
  event. RJ-6: an indel's length is part of its allele. RJ-9: a fit at the junction is overturned by later bases that
  contradict it. RJ-27 grows a pure indel's windows through the repeat its change lengthens or shortens, because
  "other lengths of that repeat are other alleles a read must span it to tell apart". No rule grows a complex event's
  window through a run that its ALT and REF share, and the judge's margin (decided 2026-10-24) has no length test.

## The same check on IMPACT and ACCESS (pat334, 2026-11-24)

- **Method.** The pat334 arms (471 tasks; tumour and, for 458 tasks, the matched normal) ran each set alone with
  `--trace`: the calls, the base partners and the del+1 partners. The patient oracle (jobs 13597285 and 13597320, the
  same `pat334_oracle.py`) judged each traced read by its own bases; the combine (13597324) summed the counts and
  passed the table checks. Only that table (counts, closed labels) left HPC. Both partners are made-up false alleles,
  expected REF in every BAM, so every ALT read on them is a false ALT.
- **Reads on the true calls (oracle ALT reads credited, and missed).** The judge finds almost every read whose bases
  carry the ALT; the rules path misses 5 to 9% of them in the tumours.

  | Family, role | Oracle-ALT reads | Missed, rules | Missed, judge | "either" reads credited, rules | judge |
  |:--|--:|--:|--:|--:|--:|
  | IMPACT, tumour | 138,328 | 9,278 | 1,432 | 1,041 | 1,953 |
  | ACCESS, tumour | 85,847 | 4,156 | 928 | 1,063 | 1,832 |
  | IMPACT, normal | 9,773 | 847 | 27 | 115 | 645 |
  | ACCESS, normal | 170 | 125 | 117 | 2 | 11 |

- **False ALT reads on the partners (each run alone).**

  | Family, role | del+1 rows | del+1, rules (SIB / either) | del+1, judge (SIB / either) | base partner, rules | judge |
  |:--|--:|--:|--:|--:|--:|
  | IMPACT, tumour | 70 | 157 (130 / 19) | 510 (178 / 320) | 2,086 | 2,105 |
  | ACCESS, tumour | 210 | 59 (24 / 20) | 568 (132 / 413) | 3,120 | 3,453 |
  | IMPACT, normal | 68 | 25 (14 / 8) | 248 (54 / 192) | 991 | 1,001 |
  | ACCESS, normal | 210 | 5 (1 / 0) | 15 (2 / 9) | 49 | 48 |

- **The three mechanisms of the HG002 trace hold on both assays.**
  1. The judge's "either" reads are most of its del+1 false ALT (320 of 510 and 413 of 568 in the tumours). By run
     band, IMPACT's lie mostly where no run follows (259 in band 1: the read ends at the change itself).
  2. The judge's SIB reads on del+1 (178 and 132) are the margin's credits; the truth is not listed.
  3. The rules path's false ALT on del+1 is almost all its exact-carrier rule (149 of 149 and 43 of 44 SIB or
     either reads in the tumours).
- **The base partner is the largest false-ALT block, the same size in both paths** (about 0.3% of its reads in the
  tumours). Every such read of the rules path comes from its exact-carrier rule (2,061 of 2,061 and 3,085 of 3,085);
  the judge's fit mechanism 2 (a read with the truth is one base from the partner's ALT and two or more from REF).
  Not yet traced read by read.
- **Open: the ACCESS normals.** Both paths miss most of their 170 oracle-ALT reads on the calls (125 and 117). Not yet
  traced.

## Slips and limits

- **Two builder faults, found by their own signs, fixed at the source, each with a selftest case that the faulty
  builder fails.**
  1. The first builder (hgd334g) took a false allele that the other haplotype carries in another representation as
     false, and judged het, hom and multi on the event's written bases only.
  2. Its first fix (hgd334h) missed insertions that the other haplotype anchors on the base before a repeat. The
     sign was multi/trans insertions falling from 746 to 144.
  hgd334g and hgd334h stay on HPC; only hgd334i is scored here.
- **The HG002 exome's index.** Every run prints htslib's "index file is older than the data file". The BAM and its
  index are links to GIAB's own files, and GIAB's index differs byte for byte from a fresh samtools index. The
  functional test (reads fetched through each index, per contig and over every row window) is in §"Index check".
- **Partners beside other events.** A false allele beside other events can still be carried by reads of a haplotype
  with another change nearby; read the clean subset first.
- **One exome, one library.** IMPACT and ACCESS (UMI consensus, ABRA2 realigned) are in §"The same check on IMPACT
  and ACCESS"; RNA and FORTE are not yet checked.

## Index check

Job idx334 (13366951, samtools 1.23) built a fresh index of the same BAM in `/scratch` and fetched reads through
each index:

| Check | Index in use (GIAB's) | Fresh index |
|:--|--:|--:|
| `samtools quickcheck` | ok | . |
| bytes | differ from byte 9 (the first contig's bin count) | . |
| contigs with a different read count | 0 | . |
| reads in the 34,931 row windows (row ± 300) | 12,370,750 | 12,370,750 |

The two files differ because another indexer wrote GIAB's, but both fetch the same reads. The htslib warning comes
from the file times only; it does not change a count.

## Next

1. The trace is done on HG002 (§"The del+1 trace") and on IMPACT and ACCESS (§"The same check on IMPACT and
   ACCESS"). Next: the operator decides whether a read whose bases fit the ALT and another run length alike (the
   "either" reads) counts as ALT or as depth only. It is one decision with the RNA judge's open cause 1
   (`333-rna-judge.md` §7, "Next": the judge credits a read that cannot tell the ALT from another length of the run,
   at short pure indels and at an SNV at a run's end). The operator chose to measure both views (2026-11-24: "Lets
   do both views"); the page for that one decision is `334-run-length-reads.md`. Then, under `decision-protocol.md`,
   the class-level options for
   the three mechanisms and the base partner's false ALT: a register entry, the effects map, prototypes in both paths,
   and measurement on DNA (HG002, IMPACT, ACCESS, WES) and RNA (HG002 RNA, FORTE) together. The fix belongs to
   invariant 7's class (a read is ALT only when its bases carry the ALT).
2. The patient part (IMPACT and ACCESS). The selection is done. `kit/local/select334_cbio.py` (selftest 0 failures)
   read the cBio files on the operator's machine. It wrote these files to `cbio/out/select334/` for the operator to
   review: the sample list, the MNV and DELINS rows of the chosen samples, a mapping template, and counts in
   `summary.txt`. It takes one sample per patient in each assay family:
   - up to 150 samples per family with somatic MNV or DELINS calls (5+ ALT reads, depth 30+);
   - up to 200 samples per family with germline-tagged MNV or DELINS rows (a sample already chosen for its
     somatic calls first).

   A slip in the first selection: its germline part checked neither the family nor the patient. A germline row
   from a study that is neither IMPACT nor ACCESS gave a sample, and a patient could give several. The operator
   saw C-* (CMO) IDs in the sample list (2026-11-23: "They should not have any cmo ids, are these from WES repo
   or Forte?"). The germline part now applies both checks, with a selftest case that the first code fails.
   `kit/local/idforms334_cbio.py` (selftest 0 failures) writes, per study folder, the samples by ID form and
   family and the selected samples by study, to `cbio/out/select334_idforms.txt` for the operator to read. A
   second fix (2026-11-24): only a DMP ID can be IMPACT or ACCESS, so a research sample with an IMPACT or ACCESS
   panel is OTHER (a selftest case that the first code fails). The selection ran again into
   `cbio/out/select334_dmp/`. It writes `ids.tsv` (one column, `sample_id`) as the input of the operator's
   mapping script, and `rows.tsv` now carries each sample's family and normal.

   The operator's `make_bam_mapping.py` (ACCESS and IMPACT keys only) makes the mapping: one CSV row per input
   tumour ID, 13 columns. They hold the study, the tumour's ACCESS duplex, simplex, standard and unfiltered BAMs
   and IMPACT standard BAM, and the matched normal's ID and BAMs of the same kinds (ACCESS and IMPACT apart). A
   missing or ambiguous match is left empty, never chosen at random. The kit reads that CSV as it is
   (`runs/pat334/in/bam_mapping.csv`) and does not use the unfiltered BAMs. The cBio files do not pair a tumour
   with its normal: `Matched_Norm_Sample_Barcode` is empty in every study's mutation files (0% of rows in the
   schema summary), and no clinical sample file has a normal column. So the normals come from the operator's
   mapping (2026-11-24), and the kit no longer reads `rows.tsv`'s `normal_sample_id`. A normal is not used when
   its ID is empty or the tumour's own, or its BAM lacks a file or index or is the tumour's BAM; each case is
   counted.

   The operator chose to copy `rows.tsv` up to a private folder in the work folder (`runs/pat334/in/`, mode 700).
   The HPC guard refuses a copy whose source is under a cBioPortal-protected root, so the operator ran that copy
   (2026-11-23: 3,048 rows on HPC, in folders that only the owner can open).
   It waits on the operator's `bam_mapping.csv`. Then the plan is: rows on the calls and a
   false-allele layer beside them, counted in the tumour and in the matched normal. The normal gives each
   germline call its het or hom stratum. It is not a REF control: an ALT read in the normal at a somatic site
   can be real (clonal haematopoiesis, tumour in the normal, an error that gives the ALT bases), and gbcms does
   not model these (operator, 2026-11-24: "the point is to measure what we see in BAM"). A read's bases decide
   whether it is ALT. The operator chose both measures (2026-11-24: "lets do both patient oracle and counts").
   `pat334_score.py` now counts only, with no expected allele: per arm, rows with AD 0, 1+ and 2+, the mean
   VAF, the AD, RD and DP sums, and the dvJ − dvR deltas, the same for tumour and normal. The patient oracle
   (v3, `scratch_kept/oraclev3/xoracle_v3.py`) judges each read by its bases. It is to be ported into the kit,
   with bwa from the work folder's `containers/bwa_0.7.17.sif` and the arms' per-read traces.
   The sanitizer runs first, and reports carry aggregates only.

   The kit is written and on HPC (2026-11-23; sha256 equal on both sides). It has three jobs:
   - `pat334_prep.sbatch` runs `pat334_rows.py`. It writes one numbered folder per sample, a task list with one
     task per sample and BAM kind, and the names file for the sanitizer. IMPACT has the standard BAM. ACCESS has
     the standard, duplex and simplex BAMs, each with the buffy coat's BAM of the same kind.
   - `pat334_arms.sbatch` runs `gbcms dna` with the rules path (dvR) and the sequence judge (dvJ) on the tumour BAM
     and the normal's BAM. Each set of rows runs on its own: the calls, the base partners and the del+1 partners.
     Rows in one VCF are counted as alleles of one site (MULTI_ALLELIC) and change each other's counts. In
     production, a wrong input allele comes alone, without the true call beside it.
   - `pat334_score.sbatch` runs `pat334_score.py`, which writes aggregates per family, BAM kind, chemistry, role,
     status, class and layer. A score is marked checked only when the sanitizer's check passes: exit 0 and
     `changed=0`.

   A local dry run on made-up data passed (the test BAMs, the test FASTA padded to genome positions, made-up IDs):
   prep, 4 tasks with 12 runs each, and the score. Every exit code was 0, every sanitizer check was `changed=0`, and
   no row broke a count relation. The dry run found two faults, both fixed before the copy. First, the sanitizer
   check always exited 0, so a job could not stop on it; it now exits 1 when masking would change a line (its
   selftest covers this). Second, calls and partners were listed in one VCF.

   A second dry run (2026-11-24) used the CSV mapping: prep made 4 tasks, all 42 runs exited 0, the score was
   checked, and no row broke a count relation. The two changed kit files are on HPC (sha256 equal on both sides).
   A third dry run (2026-11-24) used the 13-column CSV: 4 tasks, 42 runs exited 0, 3 tasks took their normal
   from the CSV, the score was checked, and no row broke a count relation. `pat334_rows.py` is on HPC (sha256
   equal on both sides).

   The operator's mapping (`bam_mapping_wnormal.csv`, 14 columns, 2026-11-24) passed `kit/check334_mapping.py`
   (job 13561672; counts only): 458 rows for the 458 `rows.tsv` samples, and every normal is of the same patient
   and assay, with its BAMs present and none equal to the tumour's. Prep (job 13561924, `pat334_prep.sbatch pat334
   bam_mapping_wnormal.csv`) made 471 tasks:

   | Family | samples in rows | with a tumour BAM | tasks | tasks with a normal | somatic calls (MNV / DELINS) | germline calls (MNV / DELINS) |
   |:--|--:|--:|--:|--:|--:|--:|
   | ACCESS | 150 | 125 | 375 (3 kinds) | 366 | 94 / 70 | 0 / 0 |
   | IMPACT | 308 | 96 | 96 | 92 | 832 / 24 | 5 / 46 |

   212 IMPACT samples and 25 ACCESS samples are not in the operator's keys. The germline part is small: 51
   germline calls, all IMPACT.

   The arms run once, with gbcms's per-read traces, after the oracle port is ready (the operator, 2026-11-24:
   "Once, with traces"), so each clinical BAM is read once for both the counts and the oracle. The arms and
   oracle jobs are shown to the operator before they are submitted.
