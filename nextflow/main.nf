#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

if (!params.input)    { exit 1, 'Input samplesheet not specified! Use --input' }
if (!params.variants) { exit 1, 'Variants file not specified! Use --variants' }
if (!params.fasta)    { exit 1, 'Reference FASTA not specified! Use --fasta' }
if (!(params.mode in ['dna', 'rna'])) {
    exit 1, "Invalid mode '${params.mode}'. Must be 'dna' or 'rna'. Use --mode dna|rna"
}

/*
========================================================================================
    IMPORT LOCAL MODULES/WORKFLOWS
========================================================================================
*/

include { GBCMS_DNA_WF }     from './workflows/dna'
include { GBCMS_RNA_WF }     from './workflows/rna'
include { FILTER_MAF }       from './modules/local/gbcms/filter_maf/main'
include { PIPELINE_SUMMARY } from './modules/local/gbcms/pipeline_summary/main'

// Helper: Check if a MAF file has at least one data row (not just header/comments)
def hasData = { file ->
    def dataLineCount = 0
    file.withReader { reader ->
        String line
        while ((line = reader.readLine()) != null && dataLineCount < 2) {
            if (!line.startsWith('#') && line.trim()) {
                dataLineCount++
            }
        }
    }
    return dataLineCount > 1
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    log.info """
    ============================================================
      gbcms v5.3.0 — Nextflow Pipeline
      Mode:     ${params.mode.toUpperCase()}
      Variants: ${params.variants}
      Output:   ${params.outdir}
      Backend:  ${params.alignment_backend}
    ============================================================
    """.stripIndent()

    //
    // STEP 1: Parse samplesheet
    //
    Channel
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
        .map { meta, bam, bai -> "${meta.id}${meta.suffix ?: ''}" }
        .collect()
        .map { keys ->
            def dupes = keys.countBy { it }.findAll { k, v -> v > 1 }
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
        ch_to_filter = ch_samplesheet.map { meta, bam, bai -> [ meta, ch_variants_file ] }

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
            .map { id, meta, bam, bai, variants -> [ meta, bam, bai, variants ] }

        // Collect ALL stats (including skipped samples) for summary
        PIPELINE_SUMMARY(
            FILTER_MAF.out.stats
                .map { meta, stats -> stats }
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
        GBCMS_RNA_WF (
            ch_ready,
            ch_fasta_tuple
        )
    } else {
        GBCMS_DNA_WF (
            ch_ready,
            ch_fasta_tuple
        )
    }
}
