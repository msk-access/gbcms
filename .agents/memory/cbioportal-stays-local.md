---
name: cbioportal-stays-local
description: "cBioPortal files never enter a tool output: no cat/head/Read/grep of them; work only through scripts that run locally and write to local files; patient slices are fine for testing (operator, 2026-10-26)."
metadata:
  type: feedback
---

On 2026-10-26 the operator asked what reaches Anthropic's servers (every tool input and
output, for the main session and every subagent) and decided:
- **Patient BAM slices** (panel, ACCESS, WES, germline normals, gate): fine for testing,
  "as long as its not used for anything downstream". The model-improvement setting is off.
- **cBioPortal files** (data_sv.txt, data_mutations*.txt, data_clinical*.txt, …): "I need
  that to stay here nothing sent from there, I need you to only write scripts and run them
  here to understand and tease what can be used".

**Why:** a tool's output is sent to Anthropic as part of the conversation. cBioPortal
files hold cohort-wide patient-level rows, unlike a test slice.

**How to apply:**
- Never Read, cat, head, grep or awk a cBioPortal file, or print its paths, rows, IDs,
  genes or positions.
- The only console channel is the fixed summary script, approved by the operator as
  "schema and counts": `python3 -I ~/test/gbcms/harness/cbio/cbio_summary.py --all`. It
  prints column names, row counts, fill rates, distinct counts, counts per value for a
  fixed set of categorical columns, and binned read counts, with small counts hidden.
- Anything finer: write a script and send all of its output (stdout and stderr) to a file
  under `~/test/gbcms/harness/cbio/out/`, then let the operator review it. Never read that
  directory with a tool.
- Enforced by `.claude/hooks/protect-cbioportal.py`. The protected roots are listed in
  `~/.claude/protected-paths.txt`, outside the repo, built from a local search: 101
  files in 18 study directories.
- An audit on 2026-10-26, counts only, found no cBioPortal data row in this project's
  session transcripts.
- Brief every agent that touches SV inputs with this rule.
- Public format docs (cBioPortal's file-format pages) are fine to read.

Related: [[parameters-are-user-settable]], [[harness-references]].
