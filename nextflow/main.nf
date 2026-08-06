#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Input validation lives inside the workflow body below: the strict syntax of
// Nextflow ≥26.04 disallows statements at the script top level (only includes,
// params, processes, workflows, and functions are permitted here).

/*
========================================================================================
    IMPORT LOCAL MODULES/WORKFLOWS
========================================================================================
*/

include { GBCMS_DNA_WF }          from './workflows/dna'
include { GBCMS_RNA_WF }          from './workflows/rna'
include { GBCMS_BUILD_GTF_CACHE } from './modules/local/gbcms/build_gtf_cache/main'
include { FILTER_MAF }            from './modules/local/gbcms/filter_maf/main'
include { PIPELINE_SUMMARY }      from './modules/local/gbcms/pipeline_summary/main'

// Helper: Check if a MAF file has at least one data row (not just header/comments).
// A function (not a closure assigned to a var): the strict syntax of Nextflow ≥26.04
// disallows the assignment-in-while-condition and postfix ++ the closure form used.
def hasData(mafFile) {
    // > 1 non-comment, non-blank line = header + at least one data row.
    // No while-loop: strict syntax (Nextflow ≥26.04) removed while loops.
    def rows = mafFile.readLines().findAll { line -> !line.startsWith('#') && line.trim() }
    return rows.size() > 1
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    // Validate required inputs (moved here from the script top level for the
    // strict syntax of Nextflow ≥26.04; `error` replaces the removed `exit`).
    if (!params.input)    { error 'Input samplesheet not specified! Use --input' }
    if (!params.variants) { error 'Variants file not specified! Use --variants' }
    if (!params.fasta)    { error 'Reference FASTA not specified! Use --fasta' }
    if (!(params.mode in ['dna', 'rna'])) {
        error "Invalid mode '${params.mode}'. Must be 'dna' or 'rna'. Use --mode dna|rna"
    }

    log.info """
    ============================================================
      gbcms v6.1.0 — Nextflow Pipeline
      Mode:     ${params.mode.toUpperCase()}
      Variants: ${params.variants}
      Output:   ${params.outdir}
      Backend:  ${params.alignment_backend}
    ============================================================
    """.stripIndent()

    //
    // STEP 1: Parse samplesheet
    //
    channel
        .fromPath(params.input)
        .splitCsv(header:true, sep:',', quote:'"')
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            
            // Optional: BAM type label for multi-BAM merge (e.g., 'duplex', 'simplex').
            // When set, auto-derives --column-prefix in the DNA module so counts
            // are pre-prefixed, and enables groupTuple-based merging in the DNA workflow.
            meta.bam_type = row.containsKey('bam_type') && row.bam_type ? row.bam_type : null

            // Optional: explicit Tumor_Sample_Barcode pattern for MAF filtering
            meta.tsb = row.containsKey('tsb') && row.tsb ? row.tsb : null

            // Per-sample suffix: explicit row.suffix > auto-derived from bam_type > global params.suffix
            // When bam_type is set but suffix is not, auto-derive suffix as "-{bam_type}"
            // to disambiguate output filenames (e.g., sample1-duplex.maf vs sample1-simplex.maf).
            if (row.containsKey('suffix') && row.suffix) {
                meta.suffix = row.suffix
            } else if (meta.bam_type) {
                meta.suffix = "-${meta.bam_type}"
            } else {
                meta.suffix = params.suffix
            }
            
            def alignment = file(row.bam, checkIfExists: true)
            
            // Handle index: if provided use it, otherwise auto-discover.
            // Supports both BAM (.bai) and CRAM (.crai) index conventions.
            def idx
            if (row.bai) {
                idx = file(row.bai, checkIfExists: true)
            } else if (row.bam.endsWith('.cram')) {
                // CRAM index conventions: .cram.crai and .crai
                def crai_path1 = "${row.bam}.crai"
                def crai_path2 = row.bam.replaceAll(/\.cram$/, '.crai')
                def crai1 = file(crai_path1)
                def crai2 = file(crai_path2)
                
                if (crai1.exists()) {
                    idx = crai1
                } else if (crai2.exists()) {
                    idx = crai2
                } else {
                    error "CRAI index not found for ${row.bam}. Searched: ${crai_path1}, ${crai_path2}"
                }
            } else {
                // BAM index conventions: .bam.bai and .bai
                def bai_path1 = "${row.bam}.bai"
                def bai_path2 = row.bam.replaceAll(/\.bam$/, '.bai')
                def bai1 = file(bai_path1)
                def bai2 = file(bai_path2)
                
                if (bai1.exists()) {
                    idx = bai1
                } else if (bai2.exists()) {
                    idx = bai2
                } else {
                    error "BAI index not found for ${row.bam}. Searched: ${bai_path1}, ${bai_path2}"
                }
            }
            
            return [ meta, alignment, idx ]
        }
        .set { ch_samplesheet }

    // Validate: duplicate sample IDs without suffix/bam_type would produce
    // identical output filenames, silently overwriting each other.
    ch_samplesheet
        .map { meta, _bam, _bai -> "${meta.id}${meta.suffix ?: ''}" }
        .collect()
        .map { keys ->
            def dupes = keys.countBy { k -> k }.findAll { _k, v -> v > 1 }
            if (dupes) {
                error "Duplicate output keys detected: ${dupes.keySet().join(', ')}. " +
                      "Multiple rows share the same sample ID without a distinguishing " +
                      "'suffix' or 'bam_type' column. This would cause output files to " +
                      "overwrite each other. Add a 'bam_type' or 'suffix' column to disambiguate."
            }
        }
    
    // Prepare reference inputs
    ch_variants_file = file(params.variants)
    ch_fasta_file = file(params.fasta)
    ch_fai_file = file("${params.fasta}.fai")
    ch_fasta_tuple = [ ch_fasta_file, ch_fai_file ]

    //
    // STEP 2: Conditional MAF filtering (shared by DNA and RNA)
    //
    if (params.filter_by_sample && ch_variants_file.name.endsWith('.maf')) {

        // Pair each sample with the shared MAF for filtering
        ch_to_filter = ch_samplesheet.map { meta, _bam, _bai -> [ meta, ch_variants_file ] }

        FILTER_MAF( ch_to_filter )

        // Skip samples with 0 matching variants
        ch_filtered_valid = FILTER_MAF.out.maf
            .filter { meta, maf ->
                if (!hasData(maf)) {
                    log.warn "Sample ${meta.id}: 0 variants after MAF filtering — skipping"
                    return false
                }
                return true
            }

        // Join filtered MAF back with BAM info
        ch_ready = ch_samplesheet
            .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
            .join( ch_filtered_valid.map { meta, maf -> [ meta.id, maf ] } )
            .map { _id, meta, bam, bai, variants -> [ meta, bam, bai, variants ] }

        // Collect ALL stats (including skipped samples) for summary
        PIPELINE_SUMMARY(
            FILTER_MAF.out.stats
                .map { _meta, stats -> stats }
                .collect()
        )

    } else {
        // No filtering — all samples get the full variants file
        ch_ready = ch_samplesheet
            .map { meta, bam, bai -> [ meta, bam, bai, ch_variants_file ] }
    }

    //
    // STEP 3: Run the appropriate workflow based on mode
    //
    if (params.mode == 'rna') {
        // M5a: pre-warm the GTF index cache ONCE for the cohort (shared --gtf +
        // --variants), then broadcast the prebuilt cache dir to every per-sample
        // GBCMS_RNA task so none of them re-parse the GTF (~9s each). Without this
        // up-front build, concurrently-launched samples all cold-miss. Disabled (no
        // GTF, or --gtf_cache false) => [] => the per-sample runs parse as before.
        if (params.gtf && params.gtf_cache) {
            GBCMS_BUILD_GTF_CACHE( [ file(params.gtf), ch_variants_file ] )
            ch_gtf_cache = GBCMS_BUILD_GTF_CACHE.out.cache_dir.first()
        } else {
            ch_gtf_cache = channel.value([])
        }

        GBCMS_RNA_WF (
            ch_ready,
            ch_fasta_tuple,
            ch_gtf_cache
        )
    } else {
        GBCMS_DNA_WF (
            ch_ready,
            ch_fasta_tuple
        )
    }
}
