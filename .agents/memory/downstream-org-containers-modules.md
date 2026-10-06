---
name: downstream-org-containers-modules
description: After each gbcms release, the mskcc-omics-workflows containers repo builds the org image and the modules repo's gbcmsrs modules pin it — both need a bump, done after the tag.
metadata:
  type: reference
---

gbcms has two downstream repos in the mskcc-omics-workflows org (checked 2026-10-05):

- **containers** (`github.com/mskcc-omics-workflows/containers`): one Dockerfile per
  version, `containers/gbcms/<X.Y.Z>/Dockerfile`. It clones the `msk-access/gbcms`
  tag, runs `maturin build --release --no-default-features` (gbcms has no crate
  features, so that flag changes nothing), and installs the wheel without our
  lock. It publishes `ghcr.io/mskcc-omics-workflows/gbcms:<X.Y.Z>`, which is not our
  `ghcr.io/msk-access/gbcms` image.
- **modules** (`github.com/mskcc-omics-workflows/modules`, `develop`):
  `modules/msk/gbcmsrs/{dna,rna,normalize,merge,buildgtfcache}` are nf-core-style
  wrappers. Files are staged as `path` inputs, and they pin the org image (6.3.1).
  `modules/msk/gbcms` is the old C++ GetBaseCountsMultiSample 1.2.5.

**How to apply:** after a release tag (the operator wants it then, not before),
first add the containers Dockerfile for the new version, then bump the gbcmsrs
modules and adapt them to that version's CLI. For 6.6.0: drop `buildgtfcache` and
the `rna` module's `gtf_cache` input (deprecated in 6.6.0, removed in 6.7.0), widen
the `gtf` pattern to `*.{gtf,gtf.gz}`, and check the 6.4–6.6 CLI and output changes
(mFSD, observations, merge). Related: [[nextflow-cli-default-divergence]].
