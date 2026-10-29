---
description: gbcms security posture — destructive ops, secrets, outward-facing actions. Enforced by hooks, not the model's judgment.
alwaysApply: true
---

# Security Posture

A rule states intent; a **hook enforces it** (a rule in a config competes with the
model's habits and loses under load). The deterministic floor lives in
`.claude/hooks/` (AI layer) and `.githooks/` (tool-agnostic). See
`docs/harness/00-plan.md` and `.agents/skills/` for wiring.

## Destructive ops
- No `rm -rf` — narrow the target or use the session scratchpad. (Enforced:
  `.claude/hooks/block-rm-rf.py`.)
- No direct push to `main`; Git Flow is `feature/* → develop → main` via PR.
  (Enforced: `.claude/hooks/block-push-main.py` + `.githooks/pre-push`.)
- Ask before bulk/irreversible operations.

## Secrets
- Never commit secrets, API keys, credentials, or PHI/patient identifiers — this
  is a clinical-genomics tool; treat sample IDs and paths as sensitive.
- Keep keys in a secret manager, not in `.env` files or hardcoded in scripts.
  (Enforced: `.claude/hooks/block-env-edits.py`.)
- Every push is secret-scanned before it leaves the machine. (Enforced:
  `.githooks/pre-push` runs gitleaks.)
- Scan tracked **code and fixtures**, not only data files — secrets hide in test
  fixtures and hardcoded example lists.
- The pre-commit hook refuses data files outside `tests/testdata/`: parquet, BAM/CRAM
  and their indexes, TSV/CSV/MAF, `.dat`, and VCF/BCF/SAM/BED/FASTQ (`.fq`), plain or `.gz`/`.bgz`
  (a VCF's header and `#CHROM` line name the sample). It also refuses a staged
  `.agents/skills/*/SKILL.md` whose frontmatter does not parse (step 3), and a staged skill or
  rule that names code, a flag or a file that does not exist (step 3b).
- Test fixtures (`tests/testdata/`) are **synthetic**. The pre-commit hook lets data
  files in there only after scanning each staged one (decoded where binary or
  compressed, so a `.vcf.gz` header is read) for MSK
  sample, CMO and anonymized-BAM identifiers, under 2 MB; a slice of a real BAM
  never goes there, renamed or not — its reads are the patient's. (Enforced:
  `.githooks/pre-commit` step 1b.) The one other data-extension file it admits is
  `scripts/regression_panel/checkpoints.tsv` (commit SHAs and PR labels), only in
  `checkpoints.py`'s shape and with no identifier in it (step 1a).

## Tool output leaves the machine
Everything a tool prints, for the session and every agent, is sent to the model provider as part
of the conversation. The operator decided (2026-10-26, after asking what reaches the provider):
- **Patient BAM slices** (panel, ACCESS, WES, germline normals, the gate) are fine for testing,
  "as long as it's not used for anything downstream". The model-improvement setting is off. Report
  aggregates and allele shapes, never sample IDs or read names.
- **cBioPortal study files** (`data_*.txt`, `meta_*.txt`, case lists: cohort-wide patient-level
  rows) stay on this machine: "I need you to only write scripts and run them here to understand
  and tease what can be used". Enforced by `.claude/hooks/protect-cbioportal.py`; the protected
  roots are in `~/.claude/protected-paths.txt`, outside the repo, built from a local search.
  - Never Read, cat, head, grep or awk a cBioPortal file, and never print its paths, rows, IDs,
    genes or positions.
  - The only console channel is the fixed summary script (schema and counts, small counts hidden):
    `python3 -I ~/test/gbcms/harness/cbio/cbio_summary.py --all`. It prints column names, row
    counts, fill rates, distinct counts, counts per value for a fixed set of categorical columns,
    and binned read counts.
  - Anything finer: write a script and send all of its output (stdout and stderr) to a file under
    `~/test/gbcms/harness/cbio/out/`, no pipes or chains, then let the operator review it. Never
    read that directory with a tool. For GitHub text that names these files, use `--body-file`.
  - Brief every agent that touches SV inputs with this rule. Public format documentation (the
    cBioPortal file-format pages) is fine to read.
  - An audit on 2026-10-26, counts only, found no cBioPortal data row in this project's session
    transcripts.

## Outward-facing actions
- **Ask before public.** The default destination for any artifact is private. A
  public host, a deploy, a published file — each needs explicit, per-artifact
  approval that names the destination. (Reflex: making a private repo public
  leaves pre-scrub secrets in git *history* — build a clean export, don't flip.)
- **No comms without approval** — no email/message/PR-to-public goes out without
  showing the draft and getting an explicit "send it."

## Injection is the real threat, not fat-fingering
With broad Bash + live web access, the realistic threat is an *injected*
instruction (web page, repo file, MCP response) chaining an allowed command into
exfiltration. Keep `curl`/`wget` off any blanket allowlist. (Defense-in-depth:
`.claude/hooks/injection-guard.py` asks on the anomalous remote-pipe / secret-exfil shapes.)
