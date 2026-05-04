process GBCMS_RNA {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/gbcms", mode: params.publish_dir_mode

    container "ghcr.io/msk-access/gbcms:4.1.0"

    input:
    tuple val(meta), path(bam), path(bai), path(variants)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path("*.{vcf,maf}"),  emit: counts
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = params.format ?: 'vcf'
    
    // Use per-sample suffix from meta, fallback to global params.suffix
    def suffix = meta.suffix ?: params.suffix
    def suffix_arg = suffix ? "--suffix ${suffix}" : ""

    // Column prefix for MAF output count columns
    def col_prefix = params.column_prefix ?: ''
    def col_prefix_arg = col_prefix ? "--column-prefix ${col_prefix}" : ""

    // Preserve original Tumor_Sample_Barcode from input MAF
    def preserve_barcode_arg = params.preserve_barcode ? "--preserve-barcode" : ""

    // Show normalization columns in output
    def show_norm_arg = params.show_normalization ? "--show-normalization" : ""

    // Adaptive context padding in repeat regions
    def adaptive_arg = params.adaptive_context ? "" : "--no-adaptive-context"

    // Alignment backend — always pass explicitly
    def backend_arg = "--alignment-backend ${params.alignment_backend}"
    def hmm_args = params.alignment_backend in ['hmm', 'pairhmm'] ? \
        "--llr-threshold ${params.llr_threshold} --gap-open-prob ${params.gap_open_prob} --gap-extend-prob ${params.gap_extend_prob} --repeat-gap-open-prob ${params.gap_open_prob_repeat} --repeat-gap-extend-prob ${params.gap_extend_prob_repeat}" : ""
    
    // Construct filter arguments
    def filters = ""
    if (params.filter_duplicates)    filters += " --filter-duplicates"
    if (params.filter_secondary)     filters += " --filter-secondary"
    if (params.filter_supplementary) filters += " --filter-supplementary"
    if (params.filter_qc_failed)     filters += " --filter-qc-failed"
    if (params.filter_improper_pair) filters += " --filter-improper-pair"
    if (params.filter_indel)         filters += " --filter-indel"

    // UMI tag (e.g., 'XM', 'RX')
    def umi_arg = params.umi_tag ? "--umi-tag ${params.umi_tag}" : ""

    // RNA-specific: REDIportal editing database
    def editing_db_arg = params.rna_editing_db ? "--rna-editing-db ${params.rna_editing_db}" : ""

    // RNA-specific: dUTP strandedness enforcement
    // CLI default for RNA is --enforce-strandedness (true), so pass --no-strandedness only if disabled
    def strandedness_arg = params.enforce_strandedness ? "" : "--no-strandedness"

    // BAQ: CLI default for RNA is --apply-baq (on). Pass --no-baq only if user explicitly disables.
    def baq_arg = params.apply_baq == false ? "--no-baq" : ""

    """
    gbcms rna \\
        --variants ${variants} \\
        --bam ${prefix}:${bam} \\
        --fasta ${fasta} \\
        --output-dir . \\
        --format ${format} \\
        ${suffix_arg} \\
        ${col_prefix_arg} \\
        ${preserve_barcode_arg} \\
        ${show_norm_arg} \\
        ${adaptive_arg} \\
        ${backend_arg} \\
        ${hmm_args} \\
        --threads ${task.cpus} \\
        --min-mapq ${params.min_mapq} \\
        --min-baseq ${params.min_baseq} \\
        --fragment-qual-threshold ${params.fragment_qual_threshold} \\
        --context-padding ${params.context_padding} \\
        ${filters} \\
        ${umi_arg} \\
        ${baq_arg} \\
        ${editing_db_arg} \\
        ${strandedness_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gbcms: \$(python -c "import gbcms; print(gbcms.__version__)" 2>/dev/null || echo "0.0.0")
    END_VERSIONS
    """
}
