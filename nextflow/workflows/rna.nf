/*
========================================================================================
    IMPORT MODULES / SUBWORKFLOWS
========================================================================================
*/

include { GBCMS_RNA } from '../modules/local/gbcms/rna/main'

/*
========================================================================================
    RUN RNA WORKFLOW
========================================================================================
*/

workflow GBCMS_RNA_WF {
    take:
    ch_samples    // channel: [ val(meta), path(bam), path(bai), path(variants) ]
    ch_fasta      // channel: [ path(fasta), path(fai) ]
    ch_gtf_cache  // value channel: prebuilt GTF index cache dir, or [] when caching is off

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Run gbcms rna (each task reuses the cohort's prebuilt GTF cache)
    //
    GBCMS_RNA (
        ch_samples,
        ch_fasta,
        ch_gtf_cache
    )
    ch_versions = ch_versions.mix(GBCMS_RNA.out.versions)

    emit:
    counts   = GBCMS_RNA.out.counts
    versions = ch_versions
}
