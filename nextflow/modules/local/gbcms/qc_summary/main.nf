/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QC_SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Summarize the QC flags of every finished gbcms output once all samples are done:
    which flags fired, on how many rows, which rows need a look.

    Input:  every *.maf / *.vcf the run wrote (collected), merged MAFs included
    Output: gbcms_qc_summary_mqc.json (and gbcms_qc_warning_mqc.json when any row needs
            a look): MultiQC custom content; <name>.qc_summary.json per output

    Runs only with --qc_summary. It reads the outputs and changes none of them.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process QC_SUMMARY {
    label 'process_low'

    publishDir "${params.outdir}/gbcms/qc_summary", mode: params.publish_dir_mode

    container "ghcr.io/msk-access/gbcms:${workflow.manifest.version}"

    input:
    path outputs  // collected gbcms outputs: *.maf and *.vcf of every sample

    output:
    path "*_mqc.json"       , emit: multiqc
    path "*.qc_summary.json", emit: json
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gbcms summarize --format multiqc --output-dir . ${outputs}
    gbcms summarize --format json --output-dir . ${outputs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gbcms: \$(python -c "import gbcms; print(gbcms.__version__)" 2>/dev/null || echo "0.0.0")
    END_VERSIONS
    """
}
