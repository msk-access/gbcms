/*
========================================================================================
    IMPORT MODULES / SUBWORKFLOWS
========================================================================================
*/

include { GBCMS_DNA } from '../modules/local/gbcms/dna/main'

/*
========================================================================================
    RUN DNA WORKFLOW
========================================================================================
*/

workflow GBCMS_DNA_WF {
    take:
    ch_samples  // channel: [ val(meta), path(bam), path(bai), path(variants) ]
    ch_fasta    // channel: [ path(fasta), path(fai) ]

    main:
    
    ch_versions = Channel.empty()

    //
    // MODULE: Run gbcms dna
    //
    GBCMS_DNA (
        ch_samples,
        ch_fasta
    )
    ch_versions = ch_versions.mix(GBCMS_DNA.out.versions)

    emit:
    counts   = GBCMS_DNA.out.counts
    versions = ch_versions
}
