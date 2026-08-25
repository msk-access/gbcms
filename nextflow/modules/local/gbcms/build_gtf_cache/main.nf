process GBCMS_BUILD_GTF_CACHE {
    label 'process_medium'

    container "ghcr.io/msk-access/gbcms:6.3.1"

    input:
    tuple path(gtf), path(variants)

    output:
    path "gtf_cache"     , emit: cache_dir
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Pre-warm the GTF index cache ONCE for the whole cohort, before the per-sample
    // GBCMS_RNA tasks fan out. Without this, concurrently-launched samples all
    // cold-miss and each re-parses the GTF (~9s); building it up front lets every
    // sample load the prebuilt index in ~0.05s. --variants must be the cohort
    // variants file the per-sample runs use, so the cache key lines up.
    """
    mkdir -p gtf_cache
    gbcms build-gtf-cache \\
        --gtf ${gtf} \\
        --variants ${variants} \\
        --gtf-cache-dir gtf_cache \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gbcms: \$(python -c "import gbcms; print(gbcms.__version__)" 2>/dev/null || echo "0.0.0")
    END_VERSIONS
    """
}
