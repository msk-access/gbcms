process MERGE_COUNTS {
    tag "$sample_id"
    label 'process_low'

    publishDir "${params.outdir}/gbcms/merged", mode: params.publish_dir_mode

    container "ghcr.io/msk-access/gbcms:6.1.0"

    input:
    tuple val(sample_id), val(bam_types), path(mafs)

    output:
    tuple val(sample_id), path("${sample_id}.merged.maf"), emit: merged_maf
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    // Build --input type:path arguments from parallel arrays
    // bam_types = ['duplex', 'simplex'], mafs = ['sample1-duplex.maf', 'sample1-simplex.maf']
    def input_args = [bam_types, mafs].transpose()
        .collect { type, maf -> "--input ${type}:${maf}" }
        .join(' ')

    // Optional: disable combined columns
    def combined_arg = params.merge_add_combined ? "" : "--no-combined"

    // Optional: use legacy naming (t_{metric}_{type})
    def legacy_arg = params.merge_legacy_naming ? "--legacy-naming" : ""

    """
    gbcms merge \\
        ${input_args} \\
        --output ${sample_id}.merged.maf \\
        ${combined_arg} \\
        ${legacy_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gbcms: \$(python -c "import gbcms; print(gbcms.__version__)" 2>/dev/null || echo "0.0.0")
    END_VERSIONS
    """
}
