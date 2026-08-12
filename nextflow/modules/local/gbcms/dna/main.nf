include { asBool } from '../../utils/main'

process GBCMS_DNA {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/gbcms", mode: params.publish_dir_mode

    container "ghcr.io/msk-access/gbcms:6.3.0"

    input:
    tuple val(meta), path(bam), path(bai), path(variants)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path("*.{vcf,maf}"),          emit: counts
    tuple val(meta), path("*.observations.parquet"),         emit: observations_parquet, optional: true
    tuple val(meta), path("*.fsd.parquet"),         emit: fsd_parquet,  optional: true
    tuple val(meta), path("*.mfsd_report.html"),    emit: mfsd_report,  optional: true
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = params.format ?: 'vcf'
    
    // Use per-sample suffix from meta, fallback to global params.suffix
    def suffix = meta.suffix ?: params.suffix
    def suffix_arg = suffix ? "--suffix ${suffix}" : ""

    // Column prefix for MAF output count columns.
    // When bam_type is set (multi-BAM merge), auto-derive prefix from it
    // (e.g., bam_type='duplex' → --column-prefix duplex_).
    // Otherwise, fall back to global params.column_prefix.
    def col_prefix = meta.bam_type ? "${meta.bam_type}_" : (params.column_prefix ?: '')
    def col_prefix_arg = col_prefix ? "--column-prefix ${col_prefix}" : ""

    // Boolean params go through asBool(): the ≥26.04 strict parser hands CLI
    // overrides to the script as Strings, and the String "false" is truthy.

    // Preserve original Tumor_Sample_Barcode from input MAF
    def preserve_barcode_arg = asBool(params.preserve_barcode) ? "--preserve-barcode" : ""

    // Show normalization columns in output
    def show_norm_arg = asBool(params.show_normalization) ? "--show-normalization" : ""

    // MNP rescue pass (v4.3.0 — decomposes MNPs into SNPs for re-counting)
    def rescue_mnp_arg = asBool(params.rescue_mnp) ? "--rescue-mnp --rescue-mnp-threshold ${params.rescue_mnp_threshold}" : ""

    // Adaptive context padding in repeat regions
    def adaptive_arg = asBool(params.adaptive_context) ? "" : "--no-adaptive-context"

    // Alignment backend — always pass explicitly
    def backend_arg = "--alignment-backend ${params.alignment_backend}"

    // Per-molecule observation export (companion Parquet; counts unchanged).
    def observations_arg = asBool(params.observations_parquet) ? "--observations-parquet" : ""
    def hmm_args = params.alignment_backend in ['hmm', 'pairhmm'] ? \
        "--llr-threshold ${params.llr_threshold} --gap-open-prob ${params.gap_open_prob} --gap-extend-prob ${params.gap_extend_prob} --repeat-gap-open-prob ${params.gap_open_prob_repeat} --repeat-gap-extend-prob ${params.gap_extend_prob_repeat}" : ""
    
    // Construct filter arguments. Always pass the explicit on/off form: omitting
    // a flag falls back to the CLI default (true for duplicates/secondary/
    // supplementary/qc-failed), which silently re-enabled filters a user had
    // turned off in Nextflow.
    def filters = ""
    filters += asBool(params.filter_duplicates)    ? " --filter-duplicates"    : " --no-filter-duplicates"
    filters += asBool(params.filter_secondary)     ? " --filter-secondary"     : " --no-filter-secondary"
    filters += asBool(params.filter_supplementary) ? " --filter-supplementary" : " --no-filter-supplementary"
    filters += asBool(params.filter_qc_failed)     ? " --filter-qc-failed"     : " --no-filter-qc-failed"
    filters += asBool(params.filter_improper_pair) ? " --filter-improper-pair" : " --no-filter-improper-pair"
    filters += asBool(params.filter_indel)         ? " --filter-indel"         : " --no-filter-indel"

    // UMI tag (e.g., 'XM', 'RX')
    def umi_arg = params.umi_tag ? "--umi-tag ${params.umi_tag}" : ""

    // BAQ: CLI default for DNA is --no-baq (off). Pass --apply-baq only if user explicitly enables.
    def baq_arg = asBool(params.apply_baq) ? "--apply-baq" : ""

    // mFSD analysis (off by default — must opt in)
    def mfsd_arg         = asBool(params.mfsd)         ? "--mfsd"          : ""
    def mfsd_parquet_arg = asBool(params.mfsd_parquet) ? "--mfsd-parquet"  : ""

    // mFSD interactive HTML report
    def mfsd_report_arg  = asBool(params.mfsd_report)  ? "--mfsd-report --mfsd-report-min-alt ${params.mfsd_report_min_alt} --mfsd-report-max-variants ${params.mfsd_report_max_variants}" : ""

    """
    gbcms dna \\
        --variants ${variants} \\
        --bam ${prefix}:${bam} \\
        --fasta ${fasta} \\
        --output-dir . \\
        --format ${format} \\
        ${suffix_arg} \\
        ${col_prefix_arg} \\
        ${preserve_barcode_arg} \\
        ${show_norm_arg} \\
        ${rescue_mnp_arg} \\
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
        ${mfsd_arg} \\
        ${mfsd_parquet_arg} \\
        ${mfsd_report_arg} \\
        ${observations_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gbcms: \$(python -c "import gbcms; print(gbcms.__version__)" 2>/dev/null || echo "0.0.0")
    END_VERSIONS
    """
}
