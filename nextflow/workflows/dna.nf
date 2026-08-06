/*
========================================================================================
    IMPORT MODULES / SUBWORKFLOWS
========================================================================================
*/

include { GBCMS_DNA }     from '../modules/local/gbcms/dna/main'
include { MERGE_COUNTS }  from '../modules/local/gbcms/merge/main'

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
    
    ch_versions = channel.empty()

    //
    // MODULE: Run gbcms dna (per-BAM genotyping)
    //
    GBCMS_DNA (
        ch_samples,
        ch_fasta
    )
    ch_versions = ch_versions.mix(GBCMS_DNA.out.versions)

    //
    // MODULE: Merge per-BAM-type MAFs (optional, requires bam_type in samplesheet)
    //
    // When merge_counts is enabled and samples have bam_type labels:
    //   1. Filter to MAF outputs from typed BAMs (bam_type != null)
    //   2. Group by sample ID using groupTuple
    //   3. Run MERGE_COUNTS to produce a single merged MAF per sample
    //
    if (params.merge_counts) {

        // Extract typed MAF outputs: [ sample_id, bam_type, maf_path ]
        ch_typed_mafs = GBCMS_DNA.out.counts
            .filter { meta, _counts -> meta.bam_type != null }
            .map { meta, counts ->
                // counts may be a list — find the MAF file
                def maf = counts instanceof List
                    ? counts.find { f -> f.name.endsWith('.maf') }
                    : (counts.name.endsWith('.maf') ? counts : null)
                if (maf == null) {
                    log.warn "Sample ${meta.id} (${meta.bam_type}): no MAF output found for merge — skipping"
                    return null
                }
                return [ meta.id, meta.bam_type, maf ]
            }
            .filter { item -> item != null }

        // Group by sample ID: [ sample_id, [type1, type2, ...], [maf1, maf2, ...] ]
        ch_grouped = ch_typed_mafs
            .groupTuple(by: 0)
            .filter { sample_id, types, _mafs ->
                if (types.size() < 2) {
                    log.warn "Sample ${sample_id}: only ${types.size()} BAM type(s) — need ≥2 for merge, skipping"
                    return false
                }
                log.info "Sample ${sample_id}: merging ${types.size()} BAM types: ${types.join(', ')}"
                return true
            }

        MERGE_COUNTS ( ch_grouped )
        ch_versions = ch_versions.mix(MERGE_COUNTS.out.versions)
    }

    emit:
    counts   = GBCMS_DNA.out.counts
    versions = ch_versions
}
