# #335 source fixes: the BAM sort-order check (A) and the judge's verdict cache (B)

Status: draft, revision 3 after the second review (2026-11-17; operator: "I would solve it at the source ...
make a plan accordingly, review it and then fix things based on the review"; "this is also for s3"; on the
second review: "Accept all defaults") · plan-first (`plan-first.md`); decision protocol
(`.agents/rules/decision-protocol.md`) · base develop 6ffde565 · found in #335 phase 2 (`335-phase2.md`:
S5's review, finding 1; S3b's design stop)
Reviews: 2 (this page, Fable, twice) · New count conditions: 0 · Builds: 0 · HPC batches: 1 (M7)
(testing is counted, not capped)

## 1. The class

Two inputs make a count depend on something other than the reads and the options:
- **A, the file's order.** gbcms reads each bin through the BAM/CRAM index and never checks that the file
  is sorted by position. htslib's region reader stops at the first record that starts past the region's
  end (`hts.c:4141`). So on a file whose records are out of order, records later in the file are missed
  with no warning: today at bin ends, with S5 at every region end. Three routes make such an input:
  - `samtools index` (1.23) indexes an unsorted CRAM, and a region query then returns 0 of 9 reads;
  - on a shuffled CRAM, a failed `samtools index` (rc 1) still leaves a `.crai`, and region queries with it
    return 1 of 40, 0 of 21, 0 of 22 and 0 of 22 reads (review 2);
  - the Nextflow samplesheet's `bai` column (`docs/nextflow/samplesheet.md:20`, `nextflow/main.nf:98-99`)
    takes an index from any path, so an index from another file reaches the engine.
- **B, the judge's verdict cache (finding F2).** The sequence judge keeps each read's call per row, keyed
  on the read's name, flags, position, CIGAR and bases. The call also reads base qualities, MAPQ, the
  AS/XS/NM/YO tags, TLEN and the mate's contig. So a second copy of a read that differs only there gets
  the first copy's call while the thread's 128 MB budget has room, and its own call when it has not. The
  room depends on what else the thread holds, so the count can change between runs, and every speed-up
  that changes the scorer's memory (S3b, the SIMD fill K, faster walks W) moves it. This is why S3b
  stopped (`scratchpad/speed2/s3b/DESIGN.md` §0).

**The principle:** an output depends only on the reads and the options. Fix each at its source: check
the input's contract while records are read; keep only caches whose answer equals a fresh computation.

**Invariants:** 1 (binning invariance) holds on valid input, which no check changes; 3 (output-aware); 7
is not touched (no allele is inferred; a broken file is refused, not mended).

## 2. What already touches it

- **Outside tools** (Fable survey, `scratchpad/speed2/source/survey.md`, sources cited there): no tool
  checks order when it reads a region (htslib, samtools, bcftools, pysam, rust-htslib, GATK4, the C++
  GetBaseCountsMultiSample); only full-pass validators (Picard ValidateSamFile, an index rebuild) see
  every record. The SAM spec makes `@HD`/`SO` optional. No surveyed tool caches a per-read result on a key
  that holds less than the result reads.
- **The code** (Opus research, `scratchpad/speed2/source/code.md`): every count opens the file once in
  `count_bam_binned_core` (`engine.rs:657-668`) before any bin; the only record loop is
  `engine.rs:1188-1206` (S5: `read_regions`, `read_whole`, the fallback). Of the judge's caches only the
  verdict map is impure; the caches outside the judge are not yet read (§4.2).
- **The review of revision 1** (Fable, `scratchpad/speed2/review_source/review.md`). Taken in revision 2:
  the header check refused a valid file; the index check was not researched; A1 does not close class A;
  D7's premise needed evidence (M8); F2-a changes two tests and `--trace`; the effects map, data types and
  register entry were missing. Not taken in revision 2, and how revision 3 handles each:

  | Review 1 item | Revision 3 |
  |:--|:--|
  | A2's cost not measured | taken: timed in M10, beside A7 |
  | read on past the stop record to the end of the inflated BGZF block | taken into A7, which reads on past the stop through the public API |
  | `.crai` slice starts in order per contig (CRAM, index only) | taken as option A8; not measured (it sees a low record moved late, not a high record moved early) |
  | F2-a's own time and peak memory | taken: M11, before D6 |
  | the comment at `decide.rs:1313-1314` (per-transcript counts and ASJD reuse the main call) goes stale | taken: F2-a rewrites it; the RNA consequence goes to #333 (§4.2) |
  | no real CRAM in any measurement | not taken, because the 29 production files M7 counted are all BAM; the tests make CRAMs (§7) |
  | ACCESS not split by generation (XS1 Marianas, XS2 fgbio) | not taken, because order and duplicate records do not depend on the consensus caller's quality model, and M7 found 0 key-sharing records in all 9 ACCESS BAMs |
  | no `--min-baseq` sweep for the F2 tests | taken: the F2 test runs at `--min-baseq` 0, 10, 20 and 30 (§7) |
  | A2, A5, A6, F2-b, F2-c, D1, D2 described, not run | A2 and A7 measured in M10; A5 and A6 stay later prototypes; F2-b and F2-c not run, because F2-a removes the defect at lower cost; D1 and D2 not run, because D0 holds on M8 |
  | no real HG002 row with exact copies drawn out for B | taken: F2-a's worked example adds one HG002 row's counts before and after, from M11 (HG002 is public) |
- **The review of revision 2** (Fable, `scratchpad/speed2/review_source2/review.md`, 8 findings; the
  operator accepted every proposed fix on 2026-11-17). Each is in this revision: A7 and M10 (finding 1);
  M9 and D8 (2); M7's record counts and the raw-record key (4); the missing effects, docs and tests (5,
  8); (c)'s rule on read errors (6); M11, the caches outside the judge and the bytes figure (7); this
  table (3).
- **S5:** its gate passed (2026-11-17: 0 differences on all 18 arms and the four whole logs, S5's log
  rule); its review fixes made the two unsorted-file cases tests that wait for A.

## 3. Measurements (counts only)

| Measure | Result |
|:--|:--|
| Local arms' BAMs (16: patient panel, ACCESS, WES, the gate tier, FORTE x2; HG002 x10) | 16 of 16 `SO:coordinate`; first 10,000 records in order in all 16 |
| Index older than its BAM by file time | 1 of 16 (WES), a valid file (htslib warns) |
| Test BAMs | `tests/testdata` 2 of 2 `SO:coordinate`; 37 of 38 Python test writers set SO; S5's Rust test BAMs have no `@HD` |
| A valid file that says `SO:queryname` (review: records in position order) | `samtools index` rc 0; a fetch returns 20 of 20; develop counts it |
| Sorted made-up files under checks (c) and (d) | 0 breaks: 12 files and 2,400 fetches (research); an edge file with pos -1, ties, no-reference-op records, secondary, supplementary, a 5 kb splice, placed-unmapped, empty contigs, 500 fetches (review) |
| S5's layout with the sorted twin's index (60 made-up files) | whole-bin fetches (develop; S5's CRAM path): (d) stops 51 (in 10 nothing was missed yet), 9 miss records with no word. Per-row regions (S5's BAM path; M10's prototype): (d) stops 37, 22 short with no word |
| Check (d) per fetch against the wide read (A7), review 2's probe (`review_source2/probe/wide_bam.py`), counted again by M10 | the probe's own truth count was wrong (`fetch(until_eof=True)` reads on from where the last fetch stopped), so its 111 is low. M10, with the engine's htslib 1.19.1, a BAM with two starts swapped and the sorted twin's index, 216 fetches: 125 missed records; (d) fired on 28 fetches, on none of the 125; A7 fired on 166, all 125 among them; 0 false fires on the sorted twin for both. CRAMs `moved` and `high_early` (own index): A7 fired on 12 of 12 and 128 of 128 missed fetches, (d) on 0 |
| Records sharing the verdict map's key, local arms | patient arms 0 groups; HG002 arms 1,212-4,839 groups of exact copies; copies that differ in a field the call reads: 0 |
| F2 on develop (made-up BAM; one ALT read and its copy with quality 5 at the ALT base) | AD 13 at `GBCMS_JUDGE_BUDGET_MB=0`, AD 14 at 128 MB |
| The judge's scratch pools (review: garbage written into every reused array before each programme) | 571 of 571 library tests pass in two modes; a late-garbage control fails (the test is live) |
| M7 (HPC, job 11643184): the census on 29 production BAMs (33 lines: 10 IMPACT and 9 ACCESS from pf73, with their rows; pf76's 10 FORTE) | 29 of 29 `SO:coordinate`, an index no older than the file, **the first 10,000 records** in order (the rest of each file was not read). The 19 DNA BAMs (23 row sets, records starting within 600 bases of a row): every line read records (kept share 0.75-0.82 IMPACT, 0.97-1.00 ACCESS); 0 records share the verdict map's key. The census keys the raw record; the engine keys it after `clip_to_fragment` (`engine.rs:1246-1249`) |
| M8 (local): develop (`s3a`) at `GBCMS_JUDGE_BUDGET_MB` 0 and unlimited against 128 MB on the gate's 7 judge arms | 0 differences at 0 MB and at unlimited: all 7 arms, 1,693 runs each, every field, and the 8 whole logs. each run's wall seconds, summed over the 7 arms (local, 4 at a time; not a speed measure; the gate prints them as "CPU-seconds", but `batch680.py` times each run by the wall clock): 10,447 at 0 MB, 9,529 at 128 MB, 9,535 unlimited (the rescue arm 1,458, 1,092, 1,044) |
| The caches outside the judge (Opus audit, `scratchpad/speed2/caches/audit.md`; 24 caches) | 23 pure; **1 impure: the spliced cache** (`carrier.rs:192-198`, used at `1542-1609`), RNA only. Its key (position, REF, ALT, splice segments) leaves out the caller's layout (masks, background, other alleles, read direction), and one cache per row (`engine.rs:1287`) serves both the row's rules and the decomposed twin's rules (`engine.rs:1385-1390`, which bring their own guard). Probe (rerun 2026-11-17): the twin's read is ALT with a fresh cache and neither once a co-listed sibling with the twin's alleles filled it under the row's rules; the reverse order flips the sibling from neither to ALT. No budget decides it: the call order does. Run at the read-call level; the engine path is traced in code, not run end to end |
| M9 (HPC, job 11694480, 2026-11-17, 1 h 2 min): the M7 BAMs in whole | all 29 production BAMs pass the whole-file order pass (htslib's own rule through `samtools index -o /dev/null`; nothing written); over the 33 row sets, 124,813 fetches (IMPACT 54,450, ACCESS 26,392, FORTE 43,971): check (d) fired 0 times, A7 0 times, narrow and wide fetches equal. Writers by @PG: IMPACT Picard MarkDuplicates then GATK (htsjdk), ACCESS bwa and ABRA2, FORTE STAR: the whole-file pass covers the htsjdk-written IMPACT files that §7's proof does not |
| M10 (local): A1 and A7, prototype on S5 (`a7.patch`, sha256 7af7947d; frozen `src_a7`; `GBCMS_ORDER_CHECK` = `a7` (default), `a1` or `off`, measurement only) | 597 Rust tests pass, clippy clean. Per file (60 BAMs, 20 CRAMs and the probe files all have fewer than 10,000 records, so (c) stops every one at open). Per fetch, as if the disorder lay past record 10,000: whole-bin fetches, A1 stops 51 of 60 BAMs, A7 59 (fetches that missed: A1 0 of 223, A7 182); per-row regions, A1 37 of 60, A7 60 (A1 0 of 168, A7 153); 20 CRAMs: both stop 20 (A1 0 of 80, A7 80). A7 leaves fetches short in 14 of the 60: the moved record's new start sits in a later 16 kb index window than the record before it, so the stale index starts the fetch past it. The shuffled CRAM with the `.crai` its failed index left: only (c) stops it. False fires on sorted files: 0 for (c), A1 and A7 (12 random files, 1,949 fetches, 37,576 records; the twins; the whole S5 suite run under A7). Cost on a made-up file (one 248 Mb contig, 40,000 rows, 28,761 fetches, 1 thread): A7 reads 1.87 M records past region ends against 139 k counted, +0.5 to 0.8 s user time over A1 and off (3.4-3.9 s); (c) 2-4 ms per count call. The gate's 18 arms (2026-11-17, against S5, beside a control run of S5 under the same load): **A7 passes**, 0 differences on all 18 arms, every field and the whole logs; **A1 passes** once its RNA fill-out log arm is run alone (side by side, three gates rewrote one shared editing-site file, so A1's run loaded it twice: a harness race, fixed in `batch680.py` by writing the file only when its content changes). CPU time over all arms 1.002 of S5 for both (tier 1.014 A7, 1.011 A1; WES 1.011, 1.009; HG002 0.950, 0.943; the 1.3-1.5 on the `logs` and `multi` arms are runs of 5-7 CPU-seconds); peak memory per run median 1.001, at most 1.24. The speed measure is pf77 on HPC |
| M11 (local, 2026-11-17): the F2-a build (`f2a_v2.patch`, sha256 0600dd4e; the review's removal on clean develop) | **gate passed**: 0 differences on all 7 judge arms (1,693 runs, every field) and in the 8 whole logs, against `s3a`. Beside a control run of `s3a` under the same load (`GBCMS_GATE_RSS=1`): CPU time 10,725 s against 10,785 s (0.994); peak memory per run, median ratio 0.995-1.011 by arm, at most 1.13 (the log arm's small runs, 302 against 267 MB); the largest peaks fall (1,267 against 1,432 MB; 1,511 against 1,614 MB). The HG002 row drawn out: open |

## 4. Options

### 4.1 A: the sort-order contract

| Option | What it does | Valid input | Broken input | Cost |
|:--|:--|:--|:--|:--|
| **A1: two record checks, each stops the sample with a clear error** | (c) prefix: the first 10,000 records follow htslib's own index rule (one block per contig, starts never fall within a contig, unplaced last); (d) per fetch: starts never fall among the records one fetch returns (dropped records included). (d) replaces S5's fallback | no change (§7) | stops where (c) or (d) sees disorder; the message names the contig and two positions and says "sort it (samtools sort) and index it again". (d) never sees the records its own fetch missed (§3): it fires only when another fetch returns the moved record | one more open and 10,000 records per count call; one compare per record. (c) on CRAM untimed |
| **A7: A1 plus the wide read** | each fetch reads from the region's start to the contig's end with the same public calls (`fetch`, `records()`, `engine.rs:1188`), stops by itself at the first start at or past the region's end, and reads 64 records more; (d) also checks those 64 | no change: the counted records are the same (probe: equal records on the sorted twin; 0 false fires in 216 fetches) | fires on most fetches that missed records (M10: 182 of 223 whole-bin, 153 of 168 per-row, 80 of 80 CRAM); stops 59 of 60 BAMs by bins, 60 of 60 by regions | 65 records past each region's end (M10's made-up file: 1.87 M extra against 139 k counted, +0.5-0.8 s); the gate's arms in M10. No `unsafe` and no chunk list |
| A2: A1 plus an opt-in full order scan | reads the whole file once before counting | off by default: no change | proves order | one full read of each file (M10); a new option (the operator has preferred no new flags) |
| A5: read each fetch's chunks to their end (no early stop, as noodles does) | sees all disorder inside the chunks the index names | no change | closes the case of a file with its own index (an unsorted CRAM) | up to one 16 kb index bin of extra records per fetch; needs htslib's chunk list (low-level or `unsafe`); A7 is its public-API form; a later prototype |
| A6: an index offset check (BAM only) | per contig, the start offset of the index's first chunk lands on a record of that contig | passes only with start offsets read from the index htslib loaded (review: end offsets and an own parser refuse every valid file) | catches an index from another file | 0.35-1.9 s on 3,366 contigs (local); `unsafe` FFI; a later prototype |
| A8: `.crai` slice starts in order per contig (CRAM only) | reads the CRAM index only | no change | sees a low record moved late, not a high record moved early | index-only, cheap; not measured; a later option |
| A3, A4, the header check | refuse non-`coordinate` SO, an older index, or `SO:queryname` | each refuses valid files (§3) | — | rejected |
| A0: keep S5's fallback | re-reads a bin whole on disorder seen in one region | — | matches develop's wrong counts; repeats `--trace` lines | rejected |

**(c)'s rule on errors and short files.** (c) stops checking at the first read error, at the file's end,
or after 10,000 records. It raises only the order error. A read error is left to the read path that
meets it today, so every error stays identical (a truncated block that no bin fetches still gives no
error). A file with fewer than 10,000 records, or with unplaced records only, is checked as far as it
goes. The switch that turns the checks off in tests is a `cfg(test)` Rust item, not a PyO3 parameter, so
`src/gbcms/_rs.pyi` does not change. gbcms's own errors stay identical; on a truncated BAM, (c) adds
htslib's own stderr lines (two `[E::bgzf_read...]` lines and one more "EOF marker is absent" warning;
M10's smoke run: exit 0, the same counts).

**A7's rule on errors.** htslib keeps a BGZF read error on the reader (`bgzf.c:998`), so a look-ahead
that reads into a damaged block would fail every later fetch on that worker. The prototype reopens the
reader: after an unreadable record before the region's end, the region's own fetch reads on from the
records already handed on; an unreadable record past the region's end only ends the look-ahead. The
truncated-BAM test fails without the reopen and passes with it.

**A1 and A7 detect disorder; neither proves order.** Only a full pass (A2) proves it. After A7, a broken
file is refused, or is still counted wrongly with no word where a low-start record moved late by more
than 64 records, or where the record sits in chunks a stale index does not name (M10: 14 of 60 files
leave some fetches short; other fetches still stop 59 of 60 by bins, 60 of 60 by regions). Under A1 by
regions, S5's two waiting tests stay silent (DP 0 and DP 2): each region's own fetch stops at the moved
record.

**Worked example (S5's review construction).** A sorted BAM with two SNV rows 5,000 bases apart in one
bin; one read of row A swaps its start with the first read of row B; the index from the sorted copy is
kept. Develop counts row A: DP 7, RD 3, AD 4; S5 by regions: DP 0, silently. With A1, a fetch that
returns the moved record stops the sample, and a file this small (fewer than 10,000 records) stops at
open by (c). Past record 10,000, per fetch (M10): by S5's regions, A1 stops 37 of 60 and A7 60 of 60.

### 4.2 B: the verdict cache and the budget

| Option | What changes | Cost | Exact |
|:--|:--|:--|:--|
| **F2-a (recommended): remove the verdict map** | the copy in F2's corner gets its own call (develop's value at budget 0); exact copies get the same call through the judged map (keyed on the whole read view); `--trace` prints the per-read judge lines once per copy (HG002 arms: 1,212-4,839 more groups); the `GBCMS_JUDGE_CACHE_TRACE` line loses its verdict-map count; the other caches get the map's room | saves per judged read the key's build and hash; frees about 253-313 bytes per entry (64 + 48 + 1 B plus the key's heap of 140-200 B for 150-250 bp reads); time and peak memory measured in M11 | yes |
| F2-b: complete the key | add qualities, MAPQ, AS/XS/NM/YO, TLEN, mate contig | +150-300 bytes per entry; kept complete by hand | yes, while complete |
| F2-c: key by record identity | the record's address | never hits in DNA | yes |

**Worked example (F2).** Two copies of one read reach an SNV row's call: same name, flags, position,
CIGAR and bases; the copy has base quality 5 at the ALT base (below `--min-baseq` 20). Today: with room,
the copy takes the first copy's call (ALT), AD 14; without room, its own (depth only), AD 13. After F2-a:
AD 13 at every budget. No such pair is in the 16 local arms or the 19 production DNA BAMs (§3). M11
adds one real HG002 row with exact copies, before and after.

**RNA (#333, 6.9.0).** The judge refuses RNA today (`engine.rs:444`), and 6.8.0 keeps RNA on the rules
path (operator, 2026-11-12, D3-B in `420-cycle-order.md`). So F2-a is DNA only in 6.8.0. When #333 lets
the judge read RNA, the per-transcript counts and ASJD call each record again in the same row (the
comment at `decide.rs:1313-1314` says they reuse the main counts' call); with no verdict map that is up to
3 calls per record per row. #333 then needs a pure reuse: a key on the record's identity, its qualities
and the read rules (F2-c with the rules in its key), measured on FORTE.

**The caches outside the judge.** The class fix covers every cache. The audit (Opus, 2026-11-17,
`scratchpad/speed2/caches/audit.md`) read 24 caches, 11 outside the judge and 13 inside it: 23 are pure,
each with its reason (among them the guard's whole-read haplotypes, `carrier.rs:2332-2354`; `Tolerance`,
`carrier.rs:655-701`, which inherits the spliced cache's fault through the windows it reads;
`CachedFasta`, `normalize/fasta.rs:10-51`; the GTF and editing-site memos); every cache whose admission a
budget decides is pure; the spliced cache is not (below, D9).

**The spliced cache (RNA rules path, in 6.8.0).** The audit's one impure cache is on the path RNA
uses in 6.8.0, so the class fix must cover it in this release. Options: **S-a (recommended): key the cache
on everything its windows read** (the layout's masks, background, other alleles, `from_right`,
`cut_left`, and the raw junction ends `own_edge` reads, `carrier.rs:1589`), exact for every caller; S-b:
give the decomposed twin its own cache in `decomp_rules` (fixes the case found, not the class); S-c: no
cache (time measured on the RNA arms). Measure first: the S-a build on the gate's RNA arms (the RNA tier's
33 runs, the RNA fill-out, the FORTE arms) and on the 10 FORTE BAMs, every field; a count change is the
decision protocol's case (register entry, worked example).

**The budget after F2-a.** Every cache is then pure (the review found no state that survives a read and
changes a call), so the budget's choice of entries changes time and memory only. Options: **D0
(recommended), keep today's 128 MB rule** (M8: unlimited room takes the same summed wall time as 128 MB, no room 9.6%
more); D2, charge sizes computed from the inputs (the same choices under every scorer, but needs an
upper-bound proof per entry); D1, count caches only (peak memory can rise; measure HG002 SV). The 128 MB
knee (`mod.rs:95-98`) was measured with the verdict map charged; M11 measures it without.

## 5. Effects map

| Place | A1 or A7 on valid input | A1 or A7 on a broken file | F2-a |
|:--|:--|:--|:--|
| MAF, VCF, mFSD, observation Parquet, every field | unchanged | the sample has no output where a check fires; else wrong as today | unchanged except F2's corner (RD/AD and strands, fragment counts and VAF where the copy decides its fragment, `partial_alt`, the flags read from them, the observation rows and mFSD class of that fragment) |
| `--verbose` log | unchanged | one ERROR line through the sample's failure; S5's per-bin DEBUG line of the failing bin and the per-BAM SEQ/QUAL warnings are not printed | unchanged |
| `--trace`, `GBCMS_JUDGE_CACHE_TRACE` | unchanged | the fallback's repeated lines are gone | judge lines once per exact copy; the cache line loses a count |
| Errors | none new; a read error stays where it is today ((c)'s rule, §4.1) | new: the order error, in Phase 0 of a bin or at open; with several threads, which bin reports first can differ between runs | none |
| `--rescue-mnp` recount | unchanged: its own count call opens the file and fetches other rows' bins, so (c) runs again and (d) checks those fetches too | (d) or A7 can fire in the recount after the main count passed; the sample then fails as a whole | unchanged |
| Tract-cluster claiming | unchanged: it claims from the same records | the order error stops the sample before any claim | unchanged |
| `gbcms merge` | unchanged | a part whose sample failed has no output to merge; what merge does then is traced in A's build before code | unchanged |
| CLI and Nextflow | unchanged | the sample fails, other samples go on, exit 1; Nextflow's `errorStrategy` is `finish` for exit 1 (`nextflow/nextflow.config:128`), so one refused BAM ends the run after the running tasks finish (D8) | unchanged |
| Data types | DNA, RNA, BAM, CRAM, fill-out (many BAMs per call: the open-time cost per BAM): one check path for all | the same | DNA only in 6.8.0 (§4.2, RNA) |
| Time, memory | +1 open and 10,000 records per count call (twice with `--rescue-mnp`); one compare per record; A7: up to 64 records per fetch region (M10) | — | a little less time per judged read; more room for the pure caches; M11 |
| Tests | the new tests in §7 | the same | two existing tests change (D6); the new tests in §7 |
| Docs | "coordinate-sorted and indexed; gbcms detects disorder it reads but cannot prove order" in `docs/reference/input-formats.md:160-161` (D5's home: it already says "Coordinate-sorted"), `docs/getting-started/index.md`, `docs/cli/dna.md`, `docs/cli/rna.md`, `docs/reference/molecule-observations.md`; the order error in `docs/resources/troubleshooting.md` (beside `:59`); the `bai` column must index this file (`docs/nextflow/samplesheet.md:20`); S5's lines; the CHANGELOG | — | the register entry (`docs/reference/read-judgment.md`); `.agents/rules/architecture.md:59` (the per-row verdict caches); the cache trace line in `docs/proposals/sequence-judge.md`; the CHANGELOG |

## 6. What this cannot see

- A: with A7, a low-start record moved late by more than 64 records, and a record in chunks a stale index
  does not name; with A1 alone, every fetch that missed records (§3) unless another fetch returns the
  moved record; a CRAM's moved record past the first 10,000 records unless a fetch returns it. M7 read
  only the first 10,000 records of each production BAM; M9 reads them in whole.
- B: F2's corner in BAMs we have not counted (M7 found 0 in 19 production DNA BAMs, keyed on raw records;
  HG002's copies are exact); RNA, which the judge does not read in 6.8.0; time and memory until M11.

## 7. Proof and tests

- **A on valid input.** htslib reads a fetch's chunks forward in file order and returns one contig's
  records (`hts.c:3415-3428`, `4101-4152`, `4141-4143`; CRAM `cram_decode.c:3480-3497`). A sorted file
  holds each contig's records in one block with starts that never fall, so (c), (d) and A7 never fire.
  (c) is htslib's own index rule (`hts.c:2482-2497`), which every BAM that htslib indexed meets; an index
  made by another tool (htsjdk) carries no such proof, which M9's whole-file pass covers. The engine's
  htslib is 1.19.1 (rust-htslib); the probes used 1.21 (pysam) and 1.23 (samtools).
- **B on valid input.** `classify` is the verdict map's only reader; the judge's other caches are pure
  (§4.2), and the caches outside the judge are read before D6; so removing the map changes a call only in
  F2's corner.
- **Tests (A):** random sorted BAMs and CRAMs: no check fires, and the counts equal the counts without the
  checks; S5's layout with the sorted twin's index stops, by regions and whole; the review's probe layout
  (two starts swapped) stops under A7 on the fetch that missed; an unsorted CRAM with its own `.crai`
  stops where the moved record is within the first 10,000 or in a fetch; a shuffled CRAM with the `.crai`
  a failed `samtools index` left; the prefix rule (contig blocks, ties, unplaced last); (c) on a truncated
  BAM and on a BAM with fewer than 10,000 records: today's error and today's counts; `SO:queryname` with
  records in order counts as today; S5's two waiting tests and its fallback test rewritten to expect the
  error (S5's own new tests). The repo has no CRAM counting test today: these are the first.
- **S5's tests that A changes** (uncommitted, S5's own; the operator hears before they land): the fallback
  test becomes `a_bin_whose_regions_give_records_out_of_order_stops_with_the_order_error` (it asserted the
  bin is read whole; it asserts the order error under A1 and A7, by regions and whole); the two waiting
  tests are un-ignored and assert the error under A7, by whole bins under A1 and from (c), and silence
  under A1 by regions (DP 0 and DP 2, as S5 counts today); two helpers change signatures only.
- **Tests (CLI):** a refused sample on the model of `tests/test_hi1_exit_codes.py`: the failed sample, no
  MAF, the observation Parquet discarded, exit 1, the other samples' outputs intact; the four count
  relations on every count that is written.
- **Tests (B):** the F2 case gives AD 13 at 0 MB and 128 MB, at `--min-baseq` 0, 10, 20 and 30, as a Rust
  in-module test (the budget override `BUDGET_FOR_TESTS` is `cfg(test)`, and the env var is read once per
  process, `mod.rs:106-119`); the review's garbage-in-the-pools test kept as a repo test; the gate's judge
  arms at 0 MB, 128 MB and unlimited with identical outputs (M8). F2-a rewrites the comments at
  `decide.rs:1313-1314` and `decide.rs:1608-1609` and drops the import it leaves unused (`mod.rs:29`).
- **Contract tests to keep:** every counting test (the four count relations); `tests/test_binning_invariance.py`;
  S5's oracles and property test; the judge's budget tests (`mod.rs:2181-2234`).
- **The gate (each fix its own build):** the 18 arms against develop, every field and the four whole
  logs: 0 differences on valid input.

## 8. Decisions for the operator

- **D5 (A), after M9 and M10:** A1, or A7 (A1 plus the wide read), with the stated limit: each detects
  disorder it reads and does not prove order. The choice rests on A7's cost against S5's gain (M10) and on
  M9 (no check fires on a production BAM). A5, A6 and A8 stay later options; A2 only if you want an opt-in
  full scan. No header check.
- **D6 (B), after M11 and the read of the caches outside the judge:** F2-a, with a register entry. It
  removes the existing test `a_records_key_is_develops_in_its_bytes_and_its_capacity` (`decide.rs:1968`:
  the verdict map's key has develop's bytes and capacity; the key no longer exists) and changes
  `a_call_takes_the_view_its_site_was_read_from_and_only_that_one` (`mod.rs:2338`: it calls `drop_view`,
  which goes with the map; the dropped-view branch and its `dropped > 500` check go, every other assertion
  stays). Evidence: the real removal (review 2, `review_source2/tree_f2a`): only these two tests change,
  570 library tests pass, the F2 case gives AD 13 at every budget. DNA only in 6.8.0; #333 inherits the
  RNA reuse rule (§4.2).
- **D7 (B, the rule for S3b):** after F2-a, a speed change may change which entries the budget admits, as
  long as every output, `--verbose` line and error stays identical; time and peak memory are measured
  (peak within 2 times `fin473`). D0 keeps 128 MB. Evidence: M8 (0 MB, 128 MB and unlimited identical on
  1,693 runs each); the judge writes every DEBUG and WARN line in the site's build, never on the per-read
  path a cache hit skips (review 2).
- **D9 (B, the spliced cache; RNA, in 6.8.0), after its measurement:** S-a, the key holds every input
  the windows read; a register entry if any RNA count changes.
- **D8 (A, Nextflow):** one refused BAM ends a Nextflow run (exit 1, `errorStrategy` `finish`).
  Recommended: keep today's rule; the error names the file, and a broken input set should stop the run.
  The other option: finish the other samples, then fail (a new exit code and a config change).

- **Decided (operator, 2026-11-18: "Yes, both"):** D5, A7 (A1 plus the wide read); D6, F2-a with its register entry, and
  the operator's OK to remove `a_records_key_is_develops_in_its_bytes_and_its_capacity` and to drop the `drop_view` part of
  `a_call_takes_the_view_its_site_was_read_from_and_only_that_one`; D7, outputs fixed while the entries the budget holds may
  vary, D0 at 128 MB; D8, keep today's rule (a refused BAM ends a Nextflow run; the error names the file); D9, S-a (the
  spliced cache keyed on every input its windows read).

## 9. Order of work, and models

1. M9 (HPC, under the standing scope; counts only), M10 and M11 (local); the read of the caches outside
   the judge (Opus); then D5, D6 and D8 come to you with the numbers (D7 can be decided now).
2. A on its own branch from develop (Opus writes, Fable reviews the diff, Sonnet gates); S5 rebased on it,
   its fallback replaced by (d) (and A7 if D5 chooses it).
3. B (F2-a, D0) on its own branch with its register entry (the same roles); then D7; then S3b in the read
   call.
4. pf77 times A, B, S2, S3a and S5 on HPC.

## 10. State

- 2026-11-17: research (survey, code, local census); revision 1; review (Fable, 7 findings); revision 2.
  S5's gate passed. M7 submitted (job 11643184).
- 2026-11-17: M7 reported (29 production BAMs sorted and in order in their first 10,000 records; 0 F2
  pairs in the 19 DNA BAMs). M8 reported (0 differences at 0 MB and unlimited).
- 2026-11-17: review 2 (Fable, 8 findings; verdict revise). The operator accepted every proposed fix
  ("Accept all defaults"); revision 3. D5 and D6 now wait for M9, M10 and M11; D8 is new.
- 2026-11-17: M9, M10 and M11 reported (all clean: §3). D5, D6, D7, D8 and D9 have their measurements.
- 2026-11-17: M9 submitted (job 11694480); M10 (Opus) and M11 running. The cache audit (Opus) found the
  spliced cache impure (RNA rules path): D9 is new.

## 11. Cost log

- Research: survey (Fable) 0.21 M tokens; code (Opus) 0.39 M; S3b design (Opus) 0.36 M; review of
  revision 1 (Fable) 0.25 M; review of revision 2 (Fable) 0.28 M. HPC: M7 (1 job, 18 min). Local: M8
  (two gate runs of the 7 judge arms, about 95 min each).
