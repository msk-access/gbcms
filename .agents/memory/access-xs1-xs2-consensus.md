---
name: access-xs1-xs2-consensus
description: "MSK-ACCESS has two collapsing generations with different base-quality scales: XS1 = Marianas (simplex Q80-90, duplex mostly Q90, masked = N/Q5), XS2 = fgbio (simplex capped Q45, duplex Q68/79/90, masked = N/Q2); split ACCESS results by generation and never read a quality threshold as meaning the same thing on both."
metadata:
  type: project
---

The operator confirmed (2026-10-26): "XS1 is Marianas and XS2 is fgbio". Measured on the
plan670 ACCESS slices (`across/slices/*_{dna,normal}_{duplex,simplex}.bam`; 2 XS1 and 3 XS2
samples; 3,000 reads per BAM):

| | XS1 (Marianas) | XS2 (fgbio) |
|---|---|---|
| headers | ABRA2 2.17, no bwa @PG | bwa 0.7.17 + ABRA2 2.22 |
| read tags | none from the collapser | MI, RX, aD/bD/cD, aE/bE/cE, ac/bc… |
| simplex quality | Q80–90 (Q90 ~74%) | capped at Q45 (99.5%) |
| duplex quality | Q50–90 (Q90 ~92% tumour, ~48% normal) | Q68 / Q79 / Q90 |
| masked positions | N at Q5 (~1.3% of simplex bases) | N at Q2 (~0.4–0.5%) |

How the scales arise (source survey, 2026-10-26, #304): fgbio gives a single-strand consensus
at most Q45 (its pre-UMI error rate), a duplex base where the strands agree the sum of the two
strand qualities (Q90 = 45+45), and N at Q2 on a tie or a no-call; then nucleo's
FilterConsensusReads `--min-base-quality 30` masks every base below Q30 to N. That is why no
XS2 base sits between Q3 and Q29. Marianas sums the strand qualities, clamped to 5–90, and
scales each strand's quality by its read count, so a simplex read built from 3 or more reads
reaches Q90. In XS1, only thresholds of 80 or more tell bases apart.

**Why:** a quality threshold above 45 masks an entire XS2 simplex BAM but none of an XS1
simplex BAM, so no single "high" `--min-baseq` means "both strands agree" across generations;
the duplex/simplex split is carried by the BAM files themselves. Thresholds 10–30 mask the
same thing in both (only the masked N positions), so ACCESS barely moves across 0–30, and a
quality-weighted mismatch cost is close to a no-op there (every base is ≥Q45).

**How to apply:** report ACCESS results split XS1 vs XS2; sweep `--min-baseq` as on other arms
but don't add an ACCESS-only high threshold without asking what it is meant to select; brief
agents that ACCESS qualities depend on the collapser. Related: [[parameters-are-user-settable]],
[[abra2-collapse-artifacts]], [[test-across-data-types]], [[harness-references]].
