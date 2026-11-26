# Repeat-aware gap models: a survey for the gbcms read judge

Date: 2026-11-23. Sources: public web pages and public source code, plus `rust/src/counting/judge/score.rs`
and `decide.rs` in this repository (read only). No patient data, no HPC.

How to read this report:
- "Q" is a Phred value: Q = -10 log10(p). Q20 is p = 0.01. Q30 is p = 0.001.
- "L" is the tract length in bases for a homopolymer. "n" is the number of repeat units for period 2 or more.
- "Read from source" means I read the code or the table itself. "From the paper text" means a fetched page.
- "Not verified" means I could not read the number at its source. These items are listed again in section 7.

The judge today (from `score.rs`): match 0, confident mismatch -4 (`MISMATCH`), gap of k bases -(6+k)
(`GAP_OPEN` 6, `GAP_EXTEND` 1), decision margin 4 (`DELTA`). A base under `--min-baseq` or N costs 0.

---

## 1. Repeat-aware indel and gap models in other tools

### 1.0 Summary table

| Tool | Form of the model | Default numbers | Fixed or learned | Costs taken from | Same cost on every haplotype? |
|---|---|---|---|---|---|
| GATK HaplotypeCaller / Mutect2, PCR indel model [`--pcr-indel-model`] | Indel open quality lowered by the number of repeat units at each read base | Flat open Q45, extend Q10; CONSERVATIVE default; floor Q10; units capped at 20; period up to 8 | Fixed formula | The read's bases | Yes (the penalty belongs to the read) |
| GATK DRAGstr [`--dragstr-params-path`] | Tables GOP, GCP, API by period (1-8) and repeat length (1-20) | Default table below; GOP capped at 40 inside the read | Learned per sample by `CalibrateDragstrModel`; default table when data are too few | The read's bases | Yes |
| DRAGEN | Same idea as DRAGstr (GOP, GCP, and an indel prior by STR period and length) | Not published as numbers in the pages I read | Learned per sample, tumour and normal apart | Not verified | Not verified |
| Octopus | Open and extend penalty arrays by motif class (A/T homopolymer, C/G homopolymer, dinucleotide, trinucleotide) and number of units; plus an SNV cap near repeats | Tables per library and sequencer (below) | Fixed built-in tables; a custom table can be profiled from data | The haplotype's sequence | No. A haplotype with a longer run has cheaper gaps |
| Strelka2 | Indel error rate by STR context (period s, count r) | Adaptive: homopolymer 4.9e-3 (r=2) to 4.5e-2 (r>=16); dinucleotide 1.0e-2 (r=2) to 1.8e-2 (r>=9). Older ramp: 5e-5 to 3e-4 | Germline: learned per sample. Somatic: fixed | The allele context (REF context for REF-to-indel, indel context for the reverse) | No, by design (two directions, two contexts) |
| bcftools / samtools mpileup | Indel quality = min(open + ext*(k-1), h*k/L) for a k-base indel in an L-base homopolymer, L >= 3 | `-o` 40, `-e` 20, `-h` 500 (100 up to release 1.12) | Fixed; per-platform profiles | The reference run | Yes (one quality per read and indel type) |
| Platypus | Gap open by homopolymer length; table of per-base indel error rates | 2.9e-5 (L 1-4) up to 1.4e-3 (L 10), then +4.3e-4 per base; extend 3; "nucprior" 2 | Fixed | The haplotype's sequence | No |
| Dindel | Indel error rate by homopolymer length in a read-to-haplotype HMM | Table S1 (not verified) | Fixed, from 1000 Genomes pilot data | The haplotype | No |
| HipSTR | Stutter model: P(up) u, P(down) d, geometric step size | Default model 0.05 / 0.05 / geometric 0.9 in frame; 0.01 / 0.01 / 0.9 out of frame | Learned per locus by EM; default when reads are few | The locus | Yes (one model per locus) |
| GangSTR | Stutter up, down, step | 0.05, 0.05, 0.90 | Fixed defaults | The locus | Yes |
| Lancet2 | Dual-affine gap in the assembly MSA step only; single affine and strict in the genotyping step | MSA: mismatch -6, gaps 6+2k or 26+1k. Genotyping: match +1, mismatch -4, gap 12+3k | Fixed | - | Yes. No repeat awareness in the gap cost |
| Mutect2 filter [`PolymeraseSlippageFilter`] | Not a read model. A site filter for a 1-unit change in an STR of 8 or more bases | Slippage rate 0.1; minimum STR length 8 | Fixed | - | - |

### 1.1 GATK HaplotypeCaller and Mutect2

Read from source (GATK master):
- Default indel qualities: `DEFAULT_INSERTION_DELETION_QUAL = 45` for both insertion and deletion open
  (`ReadUtils.java`). Flat gap continuation `gcpHMM = 10` (`LikelihoodEngineArgumentCollection.java`).
  https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/utils/read/ReadUtils.java
  https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/LikelihoodEngineArgumentCollection.java
- `--pcr-indel-model` default is `CONSERVATIVE`. The rate factors are HOSTILE 1.0, AGGRESSIVE 2.0,
  CONSERVATIVE 3.0. NONE turns the model off. The model is also turned off when DRAGstr parameters are given.
- The formula (`PairHMMLikelihoodCalculationEngine.java`):
  `Q(n) = max(10, round(40 - exp(n / (rateFactor * pi)) + 1))`, with `INITIAL_QSCORE = 40`,
  `MIN_ADJUSTED_QSCORE = 10`, n = number of repeat units found at the read base, capped at
  `MAX_REPEAT_LENGTH = 20`, repeat unit up to `MAX_STR_UNIT_LENGTH = 8` bases.
- It is applied to the read: for each read base, `readInsQuals[i-1]` and `readDelQuals[i-1]` become the
  minimum of the existing quality and `Q(n)`. Insertion and deletion get the same value.
  https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/PairHMMLikelihoodCalculationEngine.java
  https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/ReadLikelihoodCalculationEngine.java

Values of the formula (computed from the formula above; n = units, so n = L for a homopolymer):

| n (units) | 4 | 6 | 8 | 10 | 12 | 15 | 20 |
|---|---|---|---|---|---|---|---|
| CONSERVATIVE (default) | 39 | 39 | 39 | 38 | 37 | 36 | 33 |
| AGGRESSIVE | 39 | 38 | 37 | 36 | 34 | 30 | 17 |
| HOSTILE | 37 | 34 | 28 | 17 | 10 | 10 | 10 |

Note: the default model is mild. It moves Q45 only to Q38 at 10 units. It does not use the period: 10 units of
a dinucleotide get the same quality as a 10-base homopolymer.

Because the penalty is a property of the read base, each haplotype is scored with the same gap penalties.
The read is counted for an allele only if it is "informative": `LOG_10_INFORMATIVE_THRESHOLD = 0.2`
(the best allele must beat the next by 0.2 in log10 likelihood; `AlleleLikelihoods.java`).
https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/utils/genotyper/AlleleLikelihoods.java

Mutect2 uses the same engine and the same `--pcr-indel-model` argument collection. In addition,
`FilterMutectCalls` has a site-level slippage filter (`PolymeraseSlippageFilter.java`): when the reference STR
has at least `min-slippage-length` = 8 bases and the ALT differs by exactly 1 repeat unit, it tests the ALT
count against a slippage artifact with `pcr-slippage-rate` = 0.1.
https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/filtering/PolymeraseSlippageFilter.java
https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/filtering/M2FiltersArgumentCollection.java

Genotype-given-alleles: current GATK has `--alleles` ("The set of alleles to force-call regardless of
evidence", `AssemblyBasedCallerArgumentCollection.java`). The given alleles enter the haplotype set; the read
likelihood model is the same PairHMM with the same PCR indel model or DRAGstr. There is no separate rule for
known alleles.
https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/AssemblyBasedCallerArgumentCollection.java

### 1.2 DRAGstr in GATK, and `CalibrateDragstrModel`

Read from source (`DragstrParams.java`, `DragstrHyperParameters.java`, `DragstrPairHMMInputScoreImputator.java`):
- Three tables by period (rows 1-8) and repeat length in units (columns 1-20): GOP (gap open penalty), GCP
  (gap continuation penalty), API (a prior for an indel allele at such a locus; used in genotyping, not in the
  read score).
- GCP is not learned: `GCP = 10 / period` (10.0, 5.0, 3.33, 2.5, ...). So an extra base costs Q10 in a
  homopolymer, and one whole extra unit costs about Q10 for every period.
- In the read, GOP is capped: `MAX_GOP_IN_READ = 40`. The last read base gets GOP 45 and GCP 10.
- The penalty is found from the read's own bases (`DragstrReadSTRAnalyzer`: most repeated period and its
  number of repeats at each read base). Insertion open and deletion open use the same GOP.
- Estimation (`CalibrateDragstrModel`): per sample, on STR sites sampled from the reference
  (`ComposeSTRTableFile`). Grid search: GOP candidates 10 to 50 in steps of 0.25, API candidates 0 to 40 in
  steps of 1. Defaults: at least 50 loci per (period, length) cell (cells with fewer are joined with a
  neighbour length), minimum depth 10, minimum MAPQ 20, padding 5 bases, het-to-hom ratio 2.0.
- The default table is used "when there is not enough data for estimation". The source comment says the
  default values "were copied from DRAGstr documentation provided by Illumina".

Default GOP (Phred), the rows that matter here:

| Period | n=1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11-20 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 45 | 45 | 45 | 45 | 45 | 45 | 40.5 | 33.5 | 28 | 24 | 21.75 |
| 2 | 39.5 | 39.5 | 39.5 | 39.5 | 36 | 30 | 27.25 | 25 | 24.25 | 24.75 | 26.25 (26.75 at 20) |
| 3 | 38.5 | 41 | 41 | 41 | 41 | 37.5 | 35.25 | 34.75 | 34.75 | 33.25 | 33.25 down to 28.5 at 15, then 29 |
| 4 | 37.5 | 39 | 39 | 37.75 | 34 | 34 | 30.25 | 30.25 | 30.25 | 30.25 | 30.25, then 31.75 from 16 |

Default API for period 1 (Phred prior of a real indel allele): 39, 39, 37, 35, 32, 26, 20, 16, 12, 10, 8, 7, 7,
6, 6, 5, 5, 4, 4, 4. This is a prior on real variation, not an error rate; I give it because it shows that the
model treats a real length change in a long run as likely too.

https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/utils/dragstr/DragstrParams.java
https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/dragstr/DragstrHyperParameters.java
https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/DragstrPairHMMInputScoreImputator.java
https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/dragstr/CalibrateDragstrModel.java

### 1.3 DRAGEN's STR model

From the DRAGEN somatic methods preprint (page text): the "sample-specific indel error model" was first
described for the germline pipeline; it "estimates the rates of insertion and deletion errors due to PCR
amplification"; the parameters "depend on both the reference short tandem repeat (STR) period and the repeat
length"; three sets are optimised (gap open penalties, gap continuation penalties, and the probability of an
indel variant of any non-zero length) by maximising the probability of the observed read pileups; the
estimates feed the HMM that computes read likelihoods; tumour and normal can be estimated apart.
https://www.biorxiv.org/content/10.1101/2023.03.23.534011v2.full
The Nature Biotechnology paper says the HMM takes previously estimated per-sample noise levels as input
(from a search summary; page not fetched). https://www.nature.com/articles/s41587-024-02382-1
Not verified: the DRAGEN command-line option names, whether the penalty is taken from the read or from the
haplotype, and the numbers (the GATK default table above is the only published copy I found). The GATK blog
page on DRAGEN-GATK returned HTTP 403.

### 1.4 Octopus

Read from source (`error_model_factory.cpp`, `basic_repeat_based_indel_error_model.cpp/.hpp`,
`repeat_based_indel_error_model.cpp`, `repeat_based_snv_error_model.cpp`, branch develop):
- Indel model: four gap-open arrays indexed by the number of repeat units (index 0, 1, 2, ...):
  A/T homopolymer, C/G homopolymer, dinucleotide, trinucleotide (period 3 values are also used for longer
  periods). CG and GC dinucleotides get 2 less when the value is over 7. Past the end of an array the last
  value is used.
- Extend penalties (same for all libraries): homopolymer `{3,3,3,3,3,3,4,5,6,6,8,8,7,6,5,4,3}` by units;
  dinucleotide and trinucleotide `{3,3,5,4,3,2}`.
- Tables exist for library {PCR, PCR-free, 10X, MDA} x sequencer {HiSeq-2000, -2500, -4000, X10, NovaSeq,
  BGISEQ-500, PacBio, PacBio CCS}. The model is chosen by `--sequence-error-model` (for example `PCR.NOVASEQ`).
  A custom model can be made with `profiler.py` from `--data-profile` output (Octopus guide).
  https://luntergroup.github.io/octopus/docs/guides/errorModels
- The penalties are set per haplotype (`do_set_penalties(const Haplotype& haplotype, ...)`, exact tandem
  repeats of period 1 to 5 found in the haplotype sequence). So two haplotypes with different run lengths give
  different gap costs to the same read.

Gap-open penalty (Phred) by number of units, NovaSeq tables:

| Units | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 | 15 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| PCR, A/T homopolymer | 43 | 41 | 40 | 36 | 34 | 24 | 16 | 12 | 10 | 7 |
| PCR, C/G homopolymer | 42 | 40 | 37 | 33 | 27 | 17 | 12 | 9 | 6 | 4 |
| PCR, dinucleotide | 38 | 37 | 32 | 26 | 21 | 16 | 14 | 13 | 11 | 10 |
| PCR, trinucleotide | 37 | 32 | 26 | 22 | 20 | 18 | 17 | 15 | 13 | 12 |
| PCR-free, A/T homopolymer | 43 | 42 | 41 | 38 | 35 | 29 | 24 | 21 | 19 | 15 |
| PCR-free, C/G homopolymer | 43 | 40 | 37 | 33 | 28 | 19 | 15 | 12 | 7 | 5 |
| PCR-free, dinucleotide | 39 | 38 | 33 | 29 | 26 | 23 | 21 | 20 | 18 | 15 |
| PCR-free, trinucleotide | 39 | 34 | 30 | 27 | 25 | 22 | 21 | 21 | 19 | 13 |

The PCR against PCR-free difference for an A/T run of 10 is Q16 against Q24: about 6 times in probability.

- SNV model near repeats (`repeat_based_snv_error_model.cpp`): a cap on the mismatch penalty by the run
  length next to the base, per period 1, 2, 3. Homopolymer caps by run length (index = length):
  PCR-free `{125,125,60,55,50,30,20,15,12,12,10,10,10,10,8,7,6,...}`; PCR
  `{125,125,60,55,38,23,16,14,11,10,9,8,7,7,6,6,...}`. So after a run of 8 the cap is Q12 (PCR-free) or Q11
  (PCR); after a run of 10 it is Q10 or Q9. The model builds a "mask" that is the haplotype shifted by one
  base, forward and reverse. My reading: the cheap mismatch is for a read base that copies the neighbouring
  run base (a slip that looks like a substitution), not for any mismatch. Not verified: I did not read the
  HMM code that uses the mask. The SNV model is off for PacBio.

https://github.com/luntergroup/octopus/blob/develop/src/core/models/error/error_model_factory.cpp
https://github.com/luntergroup/octopus/blob/develop/src/core/models/error/basic_repeat_based_indel_error_model.hpp
https://github.com/luntergroup/octopus/blob/develop/src/core/models/error/repeat_based_indel_error_model.cpp
https://github.com/luntergroup/octopus/blob/develop/src/core/models/error/repeat_based_snv_error_model.cpp

### 1.5 Strelka2

Read from source and from the methods text in the repository (branch v2.9.x):
- "Indel sequencing errors are modeled as a process that occurs independently in each read, with some fixed
  probably of an indel error occurring as a function of the local sequence STR context." The context is the
  pair (s, r): STR pattern size and repeat count.
- Germline: learned per sample from allele counts at a subset of sites, with a two-state mixture per locus
  (clean and noisy). Clean-state error is fixed at 1e-8. To keep the parameter count low, homopolymers
  (s=1, 2 <= r <= 16) and dinucleotides (s=2, 2 <= r <= 9) are "log-linear in r"; only the end points
  (1,2), (1,16), (2,2), (2,9) are estimated; values at (1,16) and (2,9) are used for longer tracts.
  Insertion and deletion rates are set to their geometric mean for calling.
- Somatic: "a simpler non-adaptive approach ... error rates are a function of the homopolymer context length
  r only". The ramp is log-linear from 5e-5 (r=1) to 3e-4 (r>=16). The text also says the somatic indel
  model multiplies the rates by 100 "to heuristically account for assembly and mapping errors".
- Fixed "adaptiveDefault" values in the code: non-STR 8e-3; homopolymer 4.9e-3 (r=2) to 4.5e-2 (r>=16);
  dinucleotide 1.0e-2 (r=2) to 1.8e-2 (r>=9).
- Two directions with two contexts (`getIndelErrorRate`): `refToIndelErrorProb` uses the REF repeat count;
  `indelToRefErrorProb` uses the indel allele's repeat count. During calling the reversion rate is 1.8 times
  the insertion or deletion rate.
- The demo documentation turns the estimation off for small inputs (`--disableSequenceErrorEstimation`)
  (from a search summary of the user guide).

Values of the ramps (computed from the constants above):

| r (homopolymer length) | 4 | 6 | 8 | 10 | 12 | 15 | 20 |
|---|---|---|---|---|---|---|---|
| adaptiveDefault rate | 6.7e-3 | 9.2e-3 | 1.3e-2 | 1.7e-2 | 2.4e-2 | 3.8e-2 | 4.5e-2 |
| as Phred | 21.7 | 20.3 | 19.0 | 17.6 | 16.2 | 14.2 | 13.5 |
| logLinear ramp rate | 7.2e-5 | 9.1e-5 | 1.2e-4 | 1.5e-4 | 1.9e-4 | 2.7e-4 | 3.0e-4 |

https://github.com/Illumina/strelka/blob/v2.9.x/src/c%2B%2B/lib/calibration/IndelErrorModel.cpp
https://github.com/Illumina/strelka/blob/v2.9.x/docs/methods/primary/methods.tex
Paper: https://doi.org/10.1038/s41592-018-0051-x (the bioRxiv page returned HTTP 429; the paper text itself
was not fetched).

### 1.6 bcftools and samtools mpileup

Read from source (`bam2bcf_indel.c`, `mpileup.c`, `doc/bcftools.txt`, branch develop):

```
q  = openQ + extQ * (abs(l) - 1)
qh = l_run >= 3 ? tandemQ * abs(l) / l_run : 1000
seqQ = min(q, qh)
```
`l` is the indel length, `l_run` is the homopolymer length on the reference. Defaults: `-o/--open-prob` 40,
`-e/--ext-prob` 20, `-h/--tandem-qual` 500.
Documentation of `-h`: "Coefficient for modeling homopolymer errors. Given an l-long homopolymer run, the
sequencing error of an indel of size s is modeled as INT*s/l [500] ... Bcftools 1.12 and earlier had a default
of 100, which was tuned around more error prone instruments."
Profiles (`-X`): `pacbio-ccs` 300 (open 25, ext 1), `ont-sup` 110 (open 15, ext 1), `ultima` 250 (open 20,
ext 10), `1.12` 100. `--del-bias` skews deletions against insertions (default 1.0). `--poly-mqual` uses the
lowest base quality in the homopolymer.

Sequence quality of a 1-base indel (computed from the formula):

| L | 4 | 6 | 8 | 10 | 12 | 15 | 20 |
|---|---|---|---|---|---|---|---|
| `-h 500` (now) | 40 | 40 | 40 | 40 | 40 | 33 | 25 |
| `-h 100` (to 1.12) | 25 | 17 | 12 | 10 | 8 | 7 | 5 |

The form is 1/L: the quality falls in proportion to the run length, and rises in proportion to the indel
size. A 2-base indel in a 10-run costs twice a 1-base indel. The source has a note that the model does not
check that the inserted bases agree with the run ("FIXME: check if the inserted sequence is consistent with
the homopolymer run").
https://github.com/samtools/bcftools/blob/develop/bam2bcf_indel.c
https://github.com/samtools/bcftools/blob/develop/mpileup.c
https://github.com/samtools/bcftools/blob/develop/doc/bcftools.txt

### 1.7 Platypus and Dindel

Platypus, read from source (`chaplotype.pyx`, `align.c`):
```
per_base_indel_errors = [2.9e-5, 2.9e-5, 2.9e-5, 2.9e-5, 4.3e-5, 1.1e-4, 2.4e-4, 5.7e-4, 1.0e-3, 1.4e-3]
                        + [1.4e-3 + 4.3e-4*(n-10) for n in range(11,50)]
homopolq = Phred of (idx+1) * q        # tract rate = length x per-base rate
```
The gap-open penalty is written on the haplotype (`annotateWithGapOpen`: "Annotate this haplotype with a
context-specific gap open penalty, using the homopolymer model"). Gap extend is 3 and "nucprior" is 2 (Phred).

Tract rate and Phred (computed from the list):

| L | 4 | 6 | 8 | 10 | 12 | 15 | 20 |
|---|---|---|---|---|---|---|---|
| rate per read | 1.2e-4 | 6.6e-4 | 4.6e-3 | 1.4e-2 | 2.7e-2 | 5.3e-2 | 1.1e-1 |
| Phred | 39.4 | 31.8 | 23.4 | 18.5 | 15.7 | 12.7 | 9.4 |

https://github.com/andyrimmer/Platypus/blob/master/src/cython/chaplotype.pyx
https://github.com/andyrimmer/Platypus/blob/master/src/c/align.c

Dindel (Albers et al. 2011), from the paper text: "We parameterize these rates by the length of the
homopolymer run at each position"; background indel error 3.2e-5 per site per read (MAPQ > 10); in 10-base
runs the data fit "a rate of 0.5 indels per nucleotide across the data set", about 30 times the background;
up to 7 bases the spectrum looks like real variation, longer runs are dominated by errors. Reads are
realigned to candidate haplotypes (default 8) in windows of at least 120 bases. Priors: 1/1000 for a SNP
site, 1/10,000 for an indel site.
https://pmc.ncbi.nlm.nih.gov/articles/PMC3106329
Not verified: Dindel's Supplemental Table S1 (the per-length rates). The Platypus list starts at 2.9e-5 and
has the shape Dindel describes; I could not confirm that it is a copy of Table S1.

### 1.8 HipSTR and GangSTR

HipSTR (Willems et al. 2017), from the paper text and the README:
- Per-locus stutter model with three parameters: u (stutter adds units), d (stutter removes units), rho
  (geometric step size). P(delta = 0) = 1 - u - d; P(delta > 0) = u * rho * (1 - rho)^(delta - 1);
  P(delta < 0) = d * rho * (1 - rho)^(-delta - 1).
- Learned per locus with EM over all samples. The default needs 100 reads at the locus ("a good minimum
  threshold for learning stutter models").
- `--def-stutter-model`: "PGEOM=0.9 and UP=DOWN=0.05 for in-frame artifacts and PGEOM=0.9 and UP=DOWN=0.01
  for out-of-frame artifacts". So with the default model a 1-unit slip is 0.05 * 0.9 = 0.045 in each
  direction (Q13.5), and a 2-unit slip is 0.0045 (Q23.5).
- In the read-to-haplotype HMM the STR block "assumes that a read's sequence differs from the underlying
  haplotype by at most one indel whose magnitude D is a multiple of the repeat unit length M". A stutter
  deletion "can arise anywhere within the STR region"; inserted bases are "periodic copies of the STR sequence
  directly preceding the insertion".
https://pmc.ncbi.nlm.nih.gov/articles/PMC5482724
https://github.com/tfwillems/HipSTR

GangSTR, from the README: `--stutterup` 0.05, `--stutterdown` 0.05, `--stutterprob` 0.90 ("The geometric
parameter for modeling the stutter step size distribution").
https://github.com/gymreklab/GangSTR

### 1.9 Lancet2

From the documentation and the source (branch main):
- The dual-affine gap is in Phase 1 only (the SPOA multiple alignment of assembled contigs to the reference):
  match 0, mismatch -6, short-gap model open -6 extend -2, long-gap model open -26 extend -1. "The two models
  intersect at exactly 20 bp: solving 6 + 2L = 26 + 1L yields L = 20." The short model is described as the
  one for "homopolymer stutters, PCR slippage"; the long model lets the aligner cross large real indels.
- Phase 2 (reads to haplotypes, the genotyping step) is single affine and strict (`scoring_constants.h`):
  `SCORING_MATCH = 1`, `SCORING_MISMATCH = 4`, `SCORING_GAP_OPEN = 12`, `SCORING_GAP_EXTEND = 3`, and the
  second gap tier is set equal to the first (`q2 = q`, `e2 = e`). A 1-base gap costs 15, that is 3.75
  mismatches. gbcms today: 7, that is 1.75 mismatches.
- The stated reason (`genotyper.h`): "A read should perfectly match its parent haplotype with zero biological
  gaps. We use heavy gap/mismatch penalties to punish sequencer noise ... and prevent 'allele bleeding' (noisy
  reads mapping to the wrong allele)."
- Repeat awareness: none in the gap costs. Two other parts know about low complexity. (1) The read score is
  `combined = (global_score - local_raw_score - sc_penalty) + (local_pbq_score * local_identity)`; the guide
  says "a high PBQ score from a fragmented, low-identity alignment (e.g., STR artifact) is discounted toward
  zero". (2) A sequence-complexity score from longdust (`longdust_scorer.h`) is computed around variants.
https://nygenome.github.io/Lancet2/latest/guides/variant_discovery_genotyping.html
https://github.com/nygenome/Lancet2/blob/main/src/lancet/caller/scoring_constants.h
https://github.com/nygenome/Lancet2/blob/main/src/lancet/caller/genotyper.h
Not verified: the paper text (doi 10.1093/nargab/lqag036) was not fetched; the statements are from the
project's own documentation and code.

### 1.10 Genotypers of known alleles

- GATK `--alleles`: see 1.1. Same read model as discovery.
- Graphtyper: realigns reads from a region to a pangenome graph and genotypes the graph's variants (from a
  search summary of the papers). Not verified: its gap costs and any homopolymer rule.
  https://www.nature.com/articles/ng.3964 and https://pmc.ncbi.nlm.nih.gov/articles/PMC6881350/
- PanGenie: "alignment-free, k-mer-based"; it uses read k-mer counts and known haplotype paths, so it has no
  gap model at all. It "cannot detect new variant alleles not represented in the pangenome graph".
  https://pangenie.readthedocs.io/en/latest/
- HipSTR and GangSTR are genotypers of a known locus with an explicit stutter model (1.8).
- I found no genotyper of given SNV and indel alleles that has a repeat-aware gap cost and also a hard
  "wrong length is not ALT" rule. The likelihood tools weigh the read; they do not classify it (section 4).

---

## 2. Empirical slip rates

### 2.1 What the literature gives, and what it does not

There is no single public table of "P(slip of +/-1 unit per read) by tract length" for Illumina DNA with PCR,
PCR-free, and RNA-seq. The numbers below are the nearest things. The learned tool tables (Octopus by library
and sequencer, DRAGstr default, Platypus) are the only sources that give a value for each length.

| Source | Data | What was measured | Numbers |
|---|---|---|---|
| Gymrek lab analysis, 2014 (blog, and bioRxiv 043448) | HiSeq 2000, 100 bp paired; one PCR+ genome at 25.9x, one PCR-free at 37.6x; lobSTR calls with 5 or more reads | Fraction of reads at an STR that are not the called allele ("stutter rate"), by motif period, all tract lengths together | Period 1: 0.17 (PCR) and 0.027 (PCR-free), 6.3 times. Period 2: 0.038 and 0.011. Period 3: 0.011 and 0.0029. Period 4: 0.0069 and 0.0016. Period 5: 0.0074 and 0.0015. Average 0.051 and 0.011, 4.6 times |
| Fungtammasan et al. 2015 (Genome Research) | HiSeq 100 bp; PCR+ 150x, PCR-free 245x cut to 150x; plasmids on MiSeq | STR error rate by repeat number (Figure 1; values are in the figure only) | PCR+ over PCR-free: about 2.5 times for mononucleotides once the tract reaches about 11 bp; about 3.6 times for dinucleotides from about 16 bp. In plasmids with 7 or more repeats: 4.2 times (mono) and 9.3 times (di). At 5-unit mononucleotides the rate is "10^-6 to 10^-5". No difference for very short repeats (5 units or fewer) |
| Dindel (Albers 2011) | Illumina GA II, 1000 Genomes pilot | Indel error by homopolymer length | Background 3.2e-5 per site per read; 10-base runs about 30 times the background |
| Platypus list (1.7) | Fixed table in code | Per-read tract rate | L 6: 6.6e-4; L 8: 4.6e-3; L 10: 1.4e-2; L 12: 2.7e-2; L 15: 5.3e-2; L 20: 1.1e-1 |
| Octopus tables (1.4) | Built-in, by library and sequencer | Gap-open Phred | A/T run of 10: Q16 (PCR) and Q24 (PCR-free); of 15: Q10 and Q19; of 20: Q7 and Q15 |
| Shinde et al. 2003 (NAR) | Taq polymerase, cloned tracts, single-molecule PCR | Slip rate per repeat per PCR cycle | (A)n: 1.5e-2 per repeat per cycle above a threshold of 8 units (rate = 1.5e-2 * (n - 8) - 2.3e-3 per template per cycle); no mutation seen at 8 or fewer after 60 cycles. (CA)n: 3.6e-3 per repeat per cycle above 4 units |
| NAR 2025 (gkaf131), plasmids with known tracts | Illumina PCR-free and PCR, iSeq 100 and NextSeq 500 | Fraction of reads with the original allele | PCR-free keeps "90% on average" of original alleles for tracts of 15 to 25 units; (AT)20 89%, (AT)25 52%. Stutter is seen from the first PCR cycle, an average drop of 5 points against PCR-free (1.2 for (AC)15 to 8.9 for (AT)25). A/T and AT repeats lose more than AC. UMI consensus reduces the length error by no more than 2 times for 20 or more units, and by almost 8 times for 15 units |
| Fang et al. 2014 (bioRxiv 006148) | NA12878 WGS, PCR and PCR-free | Indel call quality | PCR data gives about 6.3 times more low-quality indels; poly-A/T is 80% of the low-quality indels that only the PCR data has, and 62% of those that only the PCR-free data has |
| Stoler and Nekrutenko 2021 (NAR Genom Bioinform) | 1943 public datasets, 7 Illumina instruments | Substitution errors only (indels were not measured) | Median substitution error: HiSeq X Ten 0.087%, NovaSeq 6000 0.109%, HiSeq 2500 0.112%, NextSeq 500 0.429%. Substitution errors after a homopolymer of 3 or more are 0.7% to 5.3% of all errors and rise with the homopolymer length |

URLs:
- https://gymreklab.com/science/2014/04/04/stutternoise.html and https://www.biorxiv.org/content/10.1101/043448.full
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4417121
- https://pmc.ncbi.nlm.nih.gov/articles/PMC3106329
- https://pmc.ncbi.nlm.nih.gov/articles/PMC149199
- https://academic.oup.com/nar/article/53/5/gkaf131/8044956
- https://www.biorxiv.org/content/10.1101/006148v1.full
- https://pmc.ncbi.nlm.nih.gov/articles/PMC8002175

A consistency check between the sources: the Shinde rate for (A)12 with Taq is 1.5e-2 * (12 - 8) = 6% per
template per cycle (computed from the paper's formula). That is far above what modern library polymerases and
few cycles give; the Octopus PCR value for an A/T run of 12 is Q12 (6% per read in total). Treat Shinde as an
upper limit and as evidence for the shape: no slips under a threshold near 8 bases, then a rise in proportion
to the length.

### 2.2 The shape by length

All length-resolved sources agree on the shape:
- Up to about 5 or 6 bases a homopolymer behaves like any other sequence (Q40 to Q45).
- The cost falls quickly between 6 and 10 bases, about 4 to 5 Phred per added base (Octopus PCR A/T: 34, 30,
  24, 20, 16 for 6 to 10; DRAGstr default: 45, 40.5, 33.5, 28, 24 for 6 to 10).
- After 10 to 12 bases the fall is slow (Octopus PCR: 12 at 12, 10 at 15, 7 at 20; DRAGstr default is flat at
  21.75 from 11).
- For periods 2 and 3 the same number of units costs less than a homopolymer of that many bases at small
  counts, and more at large counts (Octopus PCR at 10 units: homopolymer 16, dinucleotide 14, trinucleotide
  17; at 6 units: 34, 21, 20). For the same tract length in bases, a longer period slips less (Gymrek: 0.17,
  0.038, 0.011, 0.0069 for periods 1 to 4).
- C/G homopolymers are worse than A/T in the Octopus tables (run of 10: Q12 against Q16 with PCR). Stoler
  and Nekrutenko report that the base preference differs by instrument (NovaSeq and X Ten: A/T worse).

### 2.3 RNA-seq

- Reverse transcription as a whole: an RT indel error rate of 8.8e-6 per site (+/- 1.4e-6) and a sequencing
  indel error rate of 1.7e-5 in the same experiment; RT substitution errors 1.14e-4 (Gout et al. 2013, PNAS,
  C. elegans, replicated cDNA). "Transcription indels tend to occur in homopolymeric nucleotide runs", and all
  were 1 base long; 18 insertions and 8 deletions.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC3832031
- A second fidelity study (RNA journal 2024) reports RT deletions of 6e-6 to 2.8e-5 per base, not by
  homopolymer length (from a search summary; the page returned HTTP 403).
  https://rnajournal.cshlp.org/content/30/9/1246.full
- Slippage by the RNA polymerase itself gives real molecules of the wrong length, before any library step.
  In the bacterial endosymbiont Buchnera, transcripts of genes with poly(A) tracts show slippage in 19% of
  clones at a 9A tract and 40% to 70% at 10A to 12A tracts; variants were both shorter (8-9 A) and longer
  (11-18 A). The M-MLV control on a synthetic A11 RNA had a significantly lower error rate (P < 0.0001).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC2567471/
  From a search summary: in E. coli, slippage is "confined to single-nucleotide tracts of roughly eight or
  more". https://academic.oup.com/nar/article/42/9/5823/1256349
- Not found: a table of slip rate by tract length for human RNA-seq (RT plus PCR plus Pol II). The bacterial
  numbers do not transfer to human Pol II. For RNA the per-length rate has to come from the project's own
  measurement (HG002 RNA against GIAB, already in hand).

### 2.4 Are +1 and -1 equal?

No. With PCR, contractions are more frequent than expansions.
- Taq: contractions over expansions 14 to 1 for (CA)n and about 5 to 1 for (A)n (expansion share 0.158)
  (Shinde 2003).
- Fungtammasan 2015: "The '1-unit-del' was the most common error category for all STRs", then 1-unit
  insertion, then deletions of more than one unit, then insertions of more than one unit.
- Gymrek: stutter "has a strong bias toward decreasing the number of repeats"; PCR-free has "a slightly lower
  downward bias".
- NAR 2025: "The stutter background was increasingly biased toward deletions, as the number of PCR cycles
  increased."
- Raz et al. 2019 (NAR): the best stutter model lets a molecule lose up to 3 units or gain 1 unit per step,
  with a rate linear in the tract length ("contraction-biased"). Coefficients are in figures only (not
  verified). https://pmc.ncbi.nlm.nih.gov/articles/PMC6412005
- RNA: both directions are seen (Buchnera); in C. elegans transcription indels, insertions were more frequent
  (18 against 8). The gbcms measurement (REF molecules read 1-2 bases longer) agrees with this.

Most tools still use one value for both directions: GATK and DRAGstr (same GOP), Strelka2 (geometric mean of
the two), HipSTR and GangSTR defaults (0.05 and 0.05). HipSTR's learned model and bcftools `--del-bias`
allow a difference.

---

## 3. A conversion to the judge's scale

### 3.1 The unit

The judge's costs are integers with mismatch = 4. Let Q_mm be the Phred value that a confident mismatch
stands for. Then

    cost(p) = round( k * Q(p) ),   k = 4 / Q_mm,   Q(p) = -10 log10(p)

| Q_mm | k (cost per Phred) | cost per factor of 10 |
|---|---|---|
| 20 | 0.200 | 2.0 |
| 25 | 0.160 | 1.6 |
| 30 | 0.133 | 1.33 |

A check that the present costs are on this scale: GATK's flat indel open is Q45. With Q_mm = 25, cost =
0.16 * 45 = 7.2. The judge's 1-base gap is 6 + 1 = 7. So the present gap cost is the GATK flat default, if a
confident mismatch is read as Q25. The present decision margin (`DELTA` = 4) is then 25 Phred, a likelihood
ratio of about 300. GATK's own margin for counting a read is 0.2 log10 (2 Phred).
With Q_mm = 20 the same gap of 7 is Q35; with Q_mm = 30 it is Q52.

Each added base of a gap costs 1 today. On the Q_mm = 25 scale that is Q6.25 per base. DRAGstr uses Q10 per
base in a homopolymer (cost 1.6) and Q10 per whole unit for any period; HipSTR and GangSTR use a factor of
10 per added unit (geometric 0.9), also Q10 (cost 1.6). So the present extend cost of 1 is near the low end
of these but in range.

### 3.2 Worked table: homopolymers

Slip quality for one unit (+/-1 base), by source. "PCR" column set = Octopus PCR.NovaSeq A/T and Platypus;
"PCR-free" = Octopus PCR-free.NovaSeq A/T; DRAGstr = the default table (capped at 40 in reads).

| L | Octopus PCR A/T | Platypus | DRAGstr default | Octopus PCR-free A/T | GATK CONSERVATIVE | Slip probability, PCR library (range) |
|---|---|---|---|---|---|---|
| 6 | 34 | 31.8 | 40 (45) | 35 | 39 | 4e-4 to 7e-4 |
| 8 | 24 | 23.4 | 33.5 | 29 | 39 | 4e-3 to 5e-3 |
| 10 | 16 | 18.5 | 24 | 24 | 38 | 1.4e-2 to 2.5e-2 |
| 12 | 12 | 15.7 | 21.75 | 21 | 37 | 2.7e-2 to 6e-2 |
| 15 | 10 | 12.7 | 21.75 | 19 | 36 | 5e-2 to 1e-1 |
| 20 | 7 | 9.4 | 21.75 | 15 | 33 | 1.1e-1 to 2e-1 |

Cost on the judge's scale, cost = k * Q, before rounding (present cost of a 1-base gap: 7):

| L | PCR library, Q (mean of Octopus PCR and Platypus) | cost at Q_mm 20 | cost at Q_mm 25 | cost at Q_mm 30 | PCR-free, Q (Octopus) | cost at Q_mm 25 |
|---|---|---|---|---|---|---|
| 6 | 33 | 6.6 | 5.3 | 4.4 | 35 | 5.6 |
| 8 | 24 | 4.8 | 3.8 | 3.2 | 29 | 4.6 |
| 10 | 17 | 3.4 | 2.7 | 2.3 | 24 | 3.8 |
| 12 | 14 | 2.8 | 2.2 | 1.9 | 21 | 3.4 |
| 15 | 11 | 2.2 | 1.8 | 1.5 | 19 | 3.0 |
| 20 | 8 | 1.6 | 1.3 | 1.1 | 15 | 2.4 |

Reading the table at Q_mm = 25 for a PCR library: a 1-base slip costs about 5 at L 6, about 4 at L 8 (the
same as one confident mismatch), about 3 at L 10, about 2 at L 12 to 15, and about 1 at L 20. The cross-over
"slip = one confident mismatch" is at L 8 for a PCR library and at L 10 for a PCR-free library. This agrees
with the measured problem class ("runs of 8 or more").

The DRAGstr default and the GATK CONSERVATIVE model are much less permissive (cost 3.5 to 3.8 at L 10 to 20
for DRAGstr; 5.3 to 6.2 for GATK CONSERVATIVE at Q_mm = 25). They are the conservative end of the range.

A second unit: DRAGstr and HipSTR add about Q10 per added unit (cost about 1.6 at Q_mm = 25); the Octopus
homopolymer extend penalty is Q3 to Q8 by length (Q6 to Q8 at 8 to 12 units, Q4 at 15); bcftools doubles the
quality for a 2-base indel. So a 2-unit slip in a run of 10 is about Q25 to Q27 (cost about 4) on a PCR
library.

### 3.3 Worked table: periods 2 and 3

n = number of units. Tract length in bases = n * period.

| Period | n | Tract (bases) | Octopus PCR | Octopus PCR-free | DRAGstr default | Strelka2 adaptiveDefault | Slip probability, PCR (Octopus) | Cost at Q_mm 25, PCR | Cost at Q_mm 25, PCR-free |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 3 | 6 | 37 | 38 | 39.5 | 19.6 | 2e-4 | 5.9 | 6.1 |
| 2 | 4 | 8 | 32 | 33 | 39.5 | 19.3 | 6e-4 | 5.1 | 5.3 |
| 2 | 5 | 10 | 26 | 29 | 36 | 18.9 | 2.5e-3 | 4.2 | 4.6 |
| 2 | 6 | 12 | 21 | 26 | 30 | 18.5 | 8e-3 | 3.4 | 4.2 |
| 2 | 8 | 16 | 16 | 23 | 25 | 17.8 | 2.5e-2 | 2.6 | 3.7 |
| 2 | 10 | 20 | 14 | 21 | 24.75 | 17.4 | 4e-2 | 2.2 | 3.4 |
| 3 | 3 | 9 | 32 | 34 | 41 | - | 6e-4 | 5.1 | 5.4 |
| 3 | 4 | 12 | 26 | 30 | 41 | - | 2.5e-3 | 4.2 | 4.8 |
| 3 | 5 | 15 | 22 | 27 | 41 | - | 6e-3 | 3.5 | 4.3 |
| 3 | 6 | 18 | 20 | 25 | 37.5 | - | 1e-2 | 3.2 | 4.0 |
| 3 | 8 | 24 | 18 | 22 | 34.75 | - | 1.6e-2 | 2.9 | 3.5 |
| 3 | 10 | 30 | 17 | 21 | 33.25 | - | 2e-2 | 2.7 | 3.4 |

Measured averages to compare (Gymrek, all tract lengths together): period 2, 0.038 with PCR (Q14) and 0.011
PCR-free (Q20); period 3, 0.011 (Q20) and 0.0029 (Q25); period 4, 0.0069 (Q22) and 0.0016 (Q28).

The Strelka2 adaptiveDefault values are the "noisy-state" rates and are high at small counts (Q20 at 2 to 4
units). They are outliers against the other three sources at small counts and should not set the short end.

A slip of one unit of period 2 or 3 is a gap of 2 or 3 bases. On the present scale that gap costs 8 or 9. The
table gives a cost for the whole unit, not per base.

### 3.4 Limits of the conversion

- The tool tables are HMM gap-open penalties. Octopus and Platypus were tuned inside a likelihood model with
  per-base qualities; the integer judge has one mismatch cost. The mapping is correct to about +/-1 cost unit,
  not better.
- The choice of Q_mm moves every cost by the same factor (0.133 to 0.2 per Phred). The ratio "slip cost /
  mismatch cost" does not depend on it: it is Q_slip / Q_mm. What does depend on Q_mm is where the ratio
  passes 1.
- The values are for DNA. Section 2.3 gives the reason that RNA needs its own measurement.

---

## 4. Symmetry and decision

### 4.1 The likelihood tools do not enforce allele-length exactness; they weigh the read

In GATK, Octopus, Platypus, Strelka2 and HipSTR a read is not classified first. Each read has a likelihood
under each haplotype, and a slip is one more term in that likelihood.

- A REF tract of L and an ALT of L-1 (a deletion row). A read with L-1 bases: P(read | ALT) is about 1,
  P(read | REF) is p_slip(L). The read supports ALT by 1 / p_slip. When the slip gets cheaper, the read says
  less; it never flips. In a long run a single read of the ALT length is weak evidence. This is the design of
  all of them, and it is why DRAGstr also has the API table (a prior for a real indel at that tract).
- A read with L-2 bases (the wrong length for both). P(read | ALT) = one slip; P(read | REF) = a 2-unit slip.
  The difference is the extend term: Q10 in GATK (`gcpHMM = 10`) and in DRAGstr. Q10 is 1.0 in log10, above
  GATK's informative threshold of 0.2. So GATK counts this read in the ALT's AD. (Derived from the constants
  read in the source; I did not run GATK.) The likelihood tools give a wrong-length read to the nearest
  allele, with a small weight. None of them has a class like gbcms `partial_alt`.
- HipSTR is the most explicit: the read's length difference from each candidate allele is a stutter event
  with probability u * rho * (1 - rho)^(delta - 1) (or d for the other direction). The observed length of
  each read is kept and reported apart from the genotype (per-read sizes in the VCF), so "this read has
  length x" and "this read supports allele y" are two different outputs.
- Strelka2 models both directions with their own context: REF-to-indel error with the REF repeat count,
  indel-to-REF error with the indel allele's repeat count, and the reversion rate is 1.8 times the forward
  rate when calling. Alleles other than REF and the two candidates get the error probability of their type.
- bcftools: `indelQ` is capped by `seqQ` ("if (indelQ > seqQ) indelQ = seqQ"), so in a long run each read's
  indel evidence is capped by the homopolymer term; again weight, not class.

### 4.2 Tools that keep the gap strict on purpose

Lancet2's genotyping step is the published design nearest to "length is the allele": "A read should perfectly
match its parent haplotype with zero biological gaps. We use heavy gap/mismatch penalties to punish sequencer
noise ... and prevent 'allele bleeding' (noisy reads mapping to the wrong allele)." Its 1-base gap is 15
against a mismatch of 4. It deals with repeats by discounting the local score with the local identity, not by
making gaps cheaper.

Mutect2 keeps the read model mild (CONSERVATIVE) and handles slippage at the site: a 1-unit change in an STR
of 8 or more reference bases must have more ALT reads than a slippage rate of 0.1 explains.

### 4.3 Where the cost comes from decides whether it cancels

This is the design point that matters most for a read with an SNV allele plus a slip.

- GATK and DRAGstr take the penalty from the read's own bases. The same read has the same penalties against
  every haplotype. A slip in a run that is the same on REF and ALT then costs the same on both, the cost
  cancels in the comparison, and the SNV base decides. The size of the slip cost does not matter for the
  call; it matters only for whether the read is explained at all.
- Octopus and Platypus take the penalty from the haplotype. If the two haplotypes have different runs, the
  same slip costs less on the haplotype with the longer run. That is correct as a likelihood (a longer run
  slips more), but it gives the longer-run haplotype an advantage for any read with a length error.
- Strelka2 uses the context of the allele on each side, by design.

In the worked example of the brief (REF `AAA G AAAAAAAAAA`, ALT G>A) the SNV changes the run: REF has A3 and
A10, ALT has one run of A14 (A16 with the listed AA insertion). A haplotype-derived cost is not the same on
the two sides here; a read-derived cost is.

Arithmetic for the example, on today's scale and with the table of 3.2 (PCR library, Q_mm = 25; a 2-unit slip
in a run of 10 costs about 4):
- Today: ALT side 1 mismatch = -4; REF side 2-base gap = -8; ALT wins by 4 = `DELTA`; the read is ALT.
- With the 2-unit slip at about -4 on the REF side: ALT -4, REF -4. No side wins by 4. The read is not ALT.
  It is also not REF by the margin rule.
- For REF to win by the margin of 4, the ALT side would have to cost 8 or more, or the REF side 0. A cheaper
  slip alone removes the false ALT; it does not by itself make the read REF.
- The asymmetry in the example is that the ALT side pays nothing for 2 extra bases, because a listed
  insertion explains them, while the REF side pays for the same 2 bases as a gap. In the likelihood tools a
  listed allele is not free: it has a prior (DRAGstr API, Strelka2 theta, Dindel 1/10,000 per indel site).
  That is the published treatment of "an allele that explains the extra bases": the alternative explanation
  competes with the slip explanation on a common scale.

### 4.4 What the designs say about a read with the SNV allele and a 1-unit slip nearby

- Independent events. Strelka2: "Indel sequencing errors are modeled as a process that occurs independently
  in each read". PairHMM tools multiply the mismatch term and the gap term. So a read with the SNV base at a
  confident quality and one slip in a nearby run is "SNV allele, and one slip": it supports the SNV allele
  with the weight of the base, as long as the slip costs the same on the two haplotypes compared.
- The SNV base is judged by its base quality; some tools lower that weight near a run. Octopus caps the
  mismatch penalty next to a run (Q10 after a run of 10, PCR-free), and Strelka2 has a filter for a
  "phasing noise artifact" across a homopolymer: two candidate haplotypes of the same length with one
  mismatch, all support on one strand, the mismatch next to a homopolymer. Both are about a mismatch that is
  the copy of the neighbouring run base. In the gbcms example the confident base is the G that interrupts the
  run, which is the opposite case: the read base is not a copy of its neighbours.
- RNA-seq variant work removes sites in or next to homopolymers instead of modelling them (common practice in
  RNA editing and RNA variant pipelines, for example a filter on homopolymer runs of 5 or more). Not verified
  in this survey: I did not fetch a source for the exact filter values.

### 4.5 What none of the tools gives

No tool in this survey has both (a) a cheaper slip in long runs and (b) a hard rule that a read of the wrong
length is not the indel allele. The two published positions are: weigh everything (GATK, DRAGEN, Octopus,
Strelka2, Platypus, HipSTR), or keep gaps strict everywhere (Lancet2). The gbcms rule (`partial_alt` for a
wrong length at an indel row) is a third position. The nearest published separation of the two questions is
HipSTR's: the read's observed length and the allele it supports are reported apart. For an SNV row the run
length is a nuisance that the tools let cancel; for an indel row the run length is the allele, and the tools
reduce the read's weight but still let it lean to the nearest allele.

---

## 5. Risks

### 5.1 MSI-high tumours

- In an MSI-high tumour, a length change in a microsatellite is real and somatic. A slip model calls it noise.
  Mutect2's slippage filter is the plain example: its fixed slippage rate of 0.1 at STRs of 8 or more bases
  removes real 1-unit indels with a low allele fraction.
- Scale of the effect: across 5,930 exomes and more than 200,000 microsatellite loci, MSI was found in 14 of
  18 cancer types, mostly at mononucleotide repeats; one locus was unstable in 119 of 171 MSI-H tumours and
  11 of 446 MSS tumours (Hause et al. 2016, from a search summary).
  https://krishna.gs.washington.edu/documents/hause_nature_medicine_2016.pdf
  Also Cortes-Ciriano et al. 2017 (7,919 tumour-normal pairs, 23 cancer types)
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5467167 and Maruvka et al. 2017 (6,747 tumours)
  https://www.nature.com/articles/nbt.3966.pdf . Not verified: per-length numbers from these papers.
- The DRAGstr answer to this is to keep two things apart: the error (GOP) and the prior of a real indel at
  such a tract (API). Both rise with the tract length (API for a homopolymer of 10 is Q10).
- For a tool with a `partial_alt` class, the consequence is in the class counts, not in ALT: in an MSI-high
  sample more reads at an indel row will have a length that is neither REF nor ALT. A per-sample learned
  model (DRAGstr, Strelka2 germline) would learn the tumour's real instability as "error" and make slips
  cheaper in exactly the sample where they are real. DRAGEN estimates tumour and normal apart.

### 5.2 Low-complexity regions

- Tract detection is not uniform across tools: exact repeats only (Octopus, period 1 to 5), period up to 8 and
  20 units (GATK), homopolymers only (Platypus, bcftools, Strelka2 somatic). Interrupted and compound repeats
  get no discount or a wrong one. bcftools notes that its model does not check that the inserted bases agree
  with the run.
- A slip can be told only in a read that spans the whole tract with anchors on both sides. Strelka2 extends
  its regions "so that the first and last locus are not within a homopolymer or short tandem repeat (STR)
  region ... because alignments that do not fully span such repeats are often erroneous"; anchors must be
  outside homopolymers of 3 or more and STRs with units of 2 to 50. HipSTR uses reads that span the repeat.
- In a run, a slip and a mismatch are two explanations of the same read (the gbcms example). Cheap gaps make
  more alignments tie. Lancet2's reason for strict gaps is this "allele bleeding"; its guide also warns that
  "a noisy read might accumulate a high raw DP score simply by traversing a chaotic repeat locus."
- Substitution errors also rise after homopolymers (Stoler and Nekrutenko: 0.7% to 5.3% of all substitution
  errors, rising with length; Octopus caps mismatch quality there). A model that makes only the gap cheaper
  moves the balance between "slip" and "base change" explanations.

### 5.3 A model from one sample type used on another

- Library: PCR against PCR-free differs by 2.5 to 9 times in error rate (Fungtammasan), 6.3 times for
  homopolymers and 4.6 times on average (Gymrek), and by 8 Phred in the Octopus tables at an A/T run of 10.
  A table tuned on PCR-free WGS is too strict for capture panels with PCR; the reverse is too permissive.
- Number of PCR cycles and polymerase: stutter grows with each cycle and A/T and AT repeats lose more than AC
  (NAR 2025); polymerases differ (Raz 2019, qualitative).
- UMI and duplex consensus: consensus removes late-cycle errors but not all. Single-UMI consensus reduced
  the length error by no more than 2 times for tracts of 20 or more units and by almost 8 times at 15 units
  (NAR 2025). A table from raw reads is too permissive for consensus reads; one from consensus reads is too
  strict for raw reads.
- Sequencer: Octopus has a table per instrument; bcftools has a profile per platform; base preference (A/T
  against C/G) differs by instrument.
- RNA against DNA: RNA adds reverse transcription and transcription slippage, and both directions occur
  (2.3, 2.4). No public per-length table for human RNA-seq was found.
- Data needs of a learned model: DRAGstr needs 50 loci per (period, length) cell at depth 10 or more and
  falls back to the default table; Strelka2 turns the estimation off for small inputs; HipSTR needs 100 reads
  per locus. A targeted panel or one RNA sample has few long tracts with enough depth.
- Direction: most models are symmetric, the data are not (2.4), and the direction differs between DNA with
  PCR (shorter) and the RNA cases above (both, or longer).

---

## 6. Key numbers in one place

- GATK: indel open Q45, extend Q10; `--pcr-indel-model` CONSERVATIVE = `max(10, round(41 - exp(n / (3*pi))))`,
  Q38 at 10 units, Q33 at 20; HOSTILE reaches the floor Q10 at 12 units. Informative threshold 0.2 log10.
  Mutect2 slippage filter: STR of 8 bases or more, rate 0.1.
- DRAGstr default GOP, period 1: Q45 up to 6 bases, 40.5, 33.5, 28, 24 at 7 to 10, 21.75 from 11. GCP = 10 /
  period. GOP capped at Q40 in the read. Learned per sample on a grid (GOP 10 to 50, step 0.25).
- Octopus A/T homopolymer, NovaSeq: PCR Q34 / 24 / 16 / 12 / 10 / 7 and PCR-free Q35 / 29 / 24 / 21 / 19 / 15
  at L = 6 / 8 / 10 / 12 / 15 / 20.
- Platypus: 6.6e-4, 4.6e-3, 1.4e-2, 2.7e-2, 5.3e-2, 1.1e-1 at the same lengths.
- Strelka2: homopolymer 4.9e-3 to 4.5e-2 (r 2 to 16), dinucleotide 1.0e-2 to 1.8e-2 (r 2 to 9); older ramp
  5e-5 to 3e-4.
- bcftools: indel quality = min(40 + 20(k-1), 500 k / L); `-h` was 100 up to release 1.12.
- HipSTR default stutter 0.05 up, 0.05 down, geometric 0.9; GangSTR the same.
- Lancet2 genotyping: match +1, mismatch -4, gap 12 + 3k; dual affine (6+2k, 26+1k, cross-over 20 bases) only
  in the assembly MSA.
- Measured stutter (fraction of reads, HiSeq 2000): period 1, 0.17 with PCR and 0.027 PCR-free; period 2,
  0.038 and 0.011; period 3, 0.011 and 0.0029; period 4, 0.0069 and 0.0016.
- Taq: (A)n 1.5e-2 per repeat per cycle above 8 units, contractions 5 times expansions; (CA)n 3.6e-3 above 4
  units, contractions 14 times expansions.
- RNA: RT indel 8.8e-6 per site; Buchnera transcripts 19% slipped at A9 and 40-70% at A10 to A12.
- Conversion: cost = (4 / Q_mm) * Q. At Q_mm = 25 the present 1-base gap (7) is Q44. A 1-base slip on a PCR
  library is about cost 5 at L 6, 4 at L 8, 3 at L 10, 2 at L 12 to 15, 1 at L 20; PCR-free 6, 5, 4, 3, 3, 2.

---

## 7. Not verified

- DRAGEN: option names, default numbers, and whether penalties come from the read or the haplotype. The GATK
  blog page on DRAGEN-GATK returned HTTP 403. The DRAGEN Nature Biotechnology paper was not fetched.
- Strelka2 paper text (bioRxiv returned HTTP 429). The statements come from the methods source in the
  repository and from the code.
- Dindel Supplemental Table S1; whether the Platypus list is a copy of it.
- Fungtammasan 2015 Figure 1 (error rate for each repeat number) and its supplement: only the fold changes
  and the text values were read.
- Raz 2019 model coefficients; HipSTR's learned stutter rates in PCR and PCR-free data; Gymrek's rates by
  tract length (the blog and the preprint give rates by period only).
- Octopus: how its HMM uses the SNV mask; how the built-in tables were trained (the guide names a profiler
  script only).
- Lancet2 paper text (doi 10.1093/nargab/lqag036).
- Graphtyper's gap costs and repeat handling; Paragraph was not surveyed.
- RNA: RT slip rate by homopolymer length for any enzyme; human Pol II slippage rate by tract length; the RNA
  journal 2024 paper and the PNAS page (HTTP 403; the PNAS numbers were read from the PMC copy).
- MSI papers: per-length instability numbers (search summaries only).
- The homopolymer filters used by RNA editing and RNA variant pipelines (4.4).
- GATK's treatment of a read 2 units off (4.1) is derived from the constants, not from a run.
