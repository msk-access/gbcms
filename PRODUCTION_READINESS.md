# gbcms production readiness

Operator, 2026-10-02. Production scope: **DNA fillouts (ACCESS duplex and simplex,
IMPACT), FORTE RNA, and WES**, at default settings, with no opt-in feature
(`--rescue-mnp`, `--rescue-homopolymer` and `--mfsd` stay off).

The gate is not "every issue closed". New issues come from deeper validation (37
were closed in September, many filed by our own adversarial reviews), so a zero
count would never arrive. gbcms is production-ready when all four hold:

1. **No known silent error at default settings on in-scope data.** Every issue
   that can silently change a count, drop a row, or mislabel an output for
   DNA/RNA/WES fillouts is fixed (section 1).
2. **Validation evidence on the release candidate.** The D5 regression panel runs
   on HPC across every in-scope assay, with every changed cell attributed and
   concordance reported by stratum (section 2).
3. **Known limits published with measured bounds.** Each deferred count issue is
   on `docs/reference/bam-evidence-caveats.md` with how often it bites (section 3).
4. **A reproducible release.** Version references agree, the release page is
   generated, dependencies are pinned and audited (section 4).

Opt-in features (`--rescue-mnp`, `--rescue-homopolymer`, `--mfsd`) are off by
default in the CLI and the Nextflow pipeline, and production keeps them off
(operator, 2026-10-02): their issues (section 5) are not in the gate. A pipeline
that turns one on must clear its section 5 issues first.

## 1. Blockers: silent errors at default settings

| Issue | What can go wrong in production | Data | Group | Status |
|:--|:--|:--|:--|:--|
| C28 #202 | A read deleting a pure indel's anchor is credited REF or ALT by Phase 3's closer haplotype | DNA, WES, RNA | 2 | Decided (RJ-15) and built; final acceptance with group 2 |
| C27 #201 | A read writing the ALT across several indel ops counts partial, not ALT | all | 2 | Decided (RJ-14) and built; no real-data change |
| C16 #174 | Spurious single ALT reads in heavily masked windows (2.0 per million reads in RNA, 1.8 in IMPACT, none in ACCESS) | RNA, IMPACT, WES | 2 | Decided (RJ-16): quality-weighted evidence, at least one `--min-baseq` base's worth; RNA probes 70 → 15–19, IMPACT 7 → 1, 1–2 of 1,185 real ALT reads lost. The residual RNA junction cases are C32 (group 3) |
| R4 #185 | Strandedness not enforced at intronic loci and opposite-strand overlaps | RNA | 3 | Decided and built: intronic positions take the spanning transcripts' strand; both strands, no strand (every read counts); exon index off-by-one fixed. Truth 0 rows; final acceptance with group 3 |
| R5 #198 | Reads spliced inside a repeat tract counted as REF coverage | RNA | 3 | Decided (RJ-17) and built: informative only within one aligned block. Truth 0 rows; final acceptance with group 3 |
| O8 #186 | `OBSERVED_ALLELE`/`COEXISTING_ALLELE` count a different read set from the counts beside them | RNA | 3 | Decided and built: the diagnostic reads the counted reads. Counts unchanged |
| C32 #213 | A spliced read's windows cut at the exon edge: spurious ALT at exon-edge events (19 per 7,224 splice probes after C16) | RNA | 3 | Decided (RJ-19) and built: windows over the reference spliced at the read's junctions (all reached; a junction entering the event spliced at its edge). With C15: splice probes 19 → 2 spurious ALT reads |
| C15 #173 (RNA) | An RNA read's soft clip (the next exon's bases) read as allele evidence | RNA | 3 | Decided (RJ-18) and built: a clip reaching an exon edge or junction end is not evidence; other clips are read |
| #123 | A non-sequence ALT (IUPAC `R`) counts 0 with no warning | MAF input | 4 | Built: FAIL `NON_SEQUENCE_ALLELE`, kept in MAF output; `<NON_SEQUENCE>` in VCF output (7 sign-out rows) |
| C30 #208 | Lowercase or unprepared alleles judged inconsistently between two pure-indel paths | all input | 4 | Built: both use `allele_kind`; 0 prepared sign-out rows affected |
| I6 #218 | A `REF_MISMATCH` row says nothing about a REF written a few bases off (87 of 154 sign-out mismatches) | MAF and VCF input | 4 | Built: FAIL kept, uncounted; `gbcms_diagnostic` `REF_AT_OFFSET(k)` |
| I3 #125 | VCF to MAF: `Tumor_Seq_Allele1` empty (decide REF or empty) | VCF input | 4 | Built: the reference allele |
| I4 #126 | A maf2vcf second ALT: one allele per row (documented) or genotype both | MAF from VCF | 4 | Built: one row, one allele (documented, tested; 0 sign-out rows) |
| M4 #194 | `gbcms merge` sums NA count cells as 0 silently | ACCESS merged fillouts | 5 | Built: combined cell NA, warned per column (0 such cells in 180,348 real ones) |
| M5 #221 | `gbcms merge`: a row only a later input has gets empty annotation columns (not in the pipeline: 0 of 40 ACCESS pairs) | hand-made merges | 5 | Built: filled from the earliest later input |
| M6 #223 | `gbcms merge` takes mFSD and RNA columns from the first input only | merged `--mfsd` fillouts | 5 | Built: every gbcms column per input, set taken from the writer |
| M2 #129 | Merging outputs of different gbcms versions or representations without a warning | merged fillouts | 5 | Built: provenance kept, versions warned, pre-6.5.0 VCF shape refused (6.3% of rows split before) |
| H1 #148 | Writers (and the reference handle) not closed when a write fails: partial files | all | 5 | Built: every output written atomically |
| D6 #156 | No single reference for status reasons, diagnostics, ASJD and rescue flags | all | 5 | Open (docs) |

Done in 6.6.0 so far, each with real-data acceptance: read inputs (adapter
read-through, absent qualities, unmapped records; #211), read judgment (C27,
C28, C16; #215), C25/C26 (#210), the pure-indel read-judgment cluster (#203),
C10/C2, C1, C12, C13, C14, C21, R1, R2.

## 2. Validation evidence: the release gate

The D5 panel (#155) is the gate run, chosen for coverage rather than by hand
(`CYCLE_6.6.0_PLAN.md` § D5). Arms:

- **IMPACT and IMPACT-HEME**, **ACCESS duplex and simplex**: signed-out variants
  stratified by allele shape and size, repeat context, co-annotated structure,
  VAF and depth, X/Y/MT, MSI-high and TMB-high samples.
- **FORTE RNA**: matched IMPACT DNA and FORTE RNA for the same patients, plus the
  truth-free exon-edge probes.
- **WES (TEMPO normals)**: VAF at germline heterozygous sites should read 0.5.
- **GIAB HG002/3/4 WES and WGS**: genotyped at truth-set indels; PHI-free.

Each run reports:

- version against version: every changed cell attributed to a ticket;
- gbcms against sign-out counts and truth sets, concordance by stratum,
  disagreements adjudicated read by read;
- the read census on the synthetic read-judgment spec: every decided cell holds;
- the spurious-ALT rate on negative-control probes (DNA and RNA), per million
  reads, with a threshold the operator sets;
- run time and memory at production scale.

## 3. Documented limits (6.7.0), with the bound to measure

| Issue | Limit | Bound to report from D5 |
|:--|:--|:--|
| C7 #144 | Clip-borne ITD carriers are not rescued | FLT3-ITD ALT reads recovered by an orthogonal count |
| C18 #177 | Long events split across alignments (supplementary reads) are not joined | ALT reads lost at 50+ bp events |
| C15 #173 | Clipped carriers of pure deletions (DNA; an RNA read's clips are not evidence, RJ-18) | Clip-only carriers outside counts, by stratum |
| C31 #212 | A mate's clipped 5' end makes TLEN short (soft clips past it are clipped) | Rows changed (none in the group 1 data) |
| C6 #143 | Exact-length insertions with a sequencing error count partial | Carriers at long-insertion loci |
| R3 #178 | Catalogued RNA editing positions inside carrier windows | Rows whose window holds an editing site |
| O6 #180 | No read-orientation evidence for oxoG/FFPE artifacts (IMPACT FFPE) | Strand-skewed ALT at C>T/G>T rows |
| O5 #179 | No mapping-bias diagnostic | ALT vs REF MAPQ and clipping skew |
| C24 #195 | Stale scores in the no-reference fallback (partial only) | Rows on that path (unprepared input only) |
| C33 #214 | A REF molecule with one clear error just outside the window its ALT reading is anchored away from can count ALT (pre-existing) | Spurious ALT at the synthetic probes after RJ-16 |
| C34 #216 | A GTF that reuses a transcript ID on another chromosome merges its introns (per-transcript counts, ASJD, exon-edge rules) | Transcript IDs reused across chromosomes or strands in the production GTF (Ensembl 111: 0) |

## 4. Release and reproducibility

| Issue | Need |
|:--|:--|
| D1 #136 | CI check that every release version reference agrees |
| D2 #137 | The release workflow creates the GitHub Release page |
| D4 #139 | Dependency upgrade audit on current releases; pins where needed |
| D5 #155 | The regression panel itself (section 2) |

## 5. Gated only if production enables an opt-in feature

| Feature (default off) | Issues |
|:--|:--|
| `--rescue-mnp` | #132 (germline component in fillouts), M1 #128 (merge of flavors reporting different alleles) |
| `--rescue-homopolymer` | #111 (arbitration margins, third alleles), #112, #145, #146, #147, M1 #128 |
| `--mfsd` | S1 #153, S2 #154 (both decided) |

## 6. Not production-gating (scheduled, no effect on outputs)

#122 (a MAF deletion at Start 1 fails loudly), #124, #130, #131, #149, #152, #127,
#138, #150.

## Decisions

1. Opt-in features: none in production (operator, 2026-10-02).
2. C16: an exact-carrier ALT call needs quality-weighted evidence, one
   `--min-baseq` base's worth (RJ-16, operator 2026-10-03); the junction cases
   moved to group 3 (C32 #213, decision 4).
3. RNA mapping quality: the default stays `--min-mapq 1` (operator, 2026-10-03).
   Unique-only (as allele-specific expression pipelines filter) would cost real
   ALT reads at genes with pseudogenes (PIK3CA E545K: 10 of 159) and 1.6% of
   junction totals; `--min-mapq 0` stays a deliberate choice for pseudogene-family
   genes. Measurement in `docs/reference/read-filters.md`.
4. Group 3, RNA (operator, 2026-10-04): intronic positions take the spanning
   transcripts' strand and both-strand positions none (R4); a splice is not
   reference coverage (R5, RJ-17); diagnostics read the counted reads (O8); an
   RNA read's clips are not evidence (C15, RJ-18); a spliced read is judged on
   windows over the reference spliced at its junctions (C32, RJ-19).

Open:

1. The spurious-ALT threshold per million reads for the negative-control probes.
2. Whether any section 3 limit (FLT3-ITD sensitivity in particular) must move to
   section 1 for the clinical use of fillouts.

## Order

The remaining 6.6.0 groups already hold every section 1 and section 4 item:
group 2 (C28, C27, C16; merged), group 3 (R4, R5, O8, C15 RNA, C32; merged), group 4 (#123, C30, #125, #126, I6; built),
group 5 (M4, M2, H1, D6), group 7 (D1, D2, D4, then the D5 panel on HPC as the
gate run). Group 6 (S1, S2) stays in 6.6.0 but does not gate production (mFSD is
off).
