---
name: downstream-org-containers-modules
description: The mskcc-omics-workflows containers repo builds an org gbcms image (6.3.1) that the modules repo's gbcmsrs modules pin; bumping them is not a release step (premature, 2026-10-07), only when the operator asks.
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

**How to apply:** don't touch these repos as a release step. The operator decided
2026-10-07, at the 6.6.0 cut, that bumping them is premature; leave them until the
operator asks. The reason: our own pipeline (`nextflow/`) doesn't use them. It
includes only `modules/local/gbcms/{dna,rna,merge,normalize,filter_maf,pipeline_summary}`
(no modules.json), on our image `ghcr.io/msk-access/gbcms:${workflow.manifest.version}`.
Moving our pipeline onto the org modules comes first, and it is its own project.
The operator starts that whole process (pipeline switch, image, modules) when they
are satisfied the org route is fit for production use; never propose it as a step.
When that happens: first the containers Dockerfile for the version, then the gbcmsrs
modules adapted to that version's CLI (since 6.3.1: drop
`buildgtfcache` and the `rna` module's `gtf_cache` input, deprecated in 6.6.0 and
removed in 6.7.0; widen the `gtf` pattern to `*.{gtf,gtf.gz}`; check the mFSD,
observations and merge changes). The modules are buehlere's (modules#259,
2026-08-27); ask whether to open a PR or coordinate first.
Related: [[nextflow-cli-default-divergence]].
