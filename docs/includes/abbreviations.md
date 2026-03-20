<!-- Abbreviation definitions for gbcms documentation -->
<!-- These provide hover tooltips for technical acronyms -->
<span hidden>abbreviations</span>

*[VAF]: Variant Allele Frequency - Fraction of reads supporting the variant allele
*[MAF]: Mutation Annotation Format - Standard mutation annotation file format
*[VCF]: Variant Call Format - Standard variant call file format
*[BAM]: Binary Alignment Map - Compressed sequence alignment file format
*[BAI]: BAM Index - Index file for rapid BAM access
*[FASTA]: FASTA format - Reference genome sequence file
*[FAI]: FASTA Index - Index file for rapid FASTA access
*[MAPQ]: Mapping Quality - Confidence score for read alignment
*[BASEQ]: Base Quality - Phred-scaled base call confidence
*[cfDNA]: Cell-free DNA - Circulating DNA fragments in blood plasma
*[ctDNA]: Circulating tumor DNA - Tumor-derived cfDNA fragments
*[HPC]: High-Performance Computing - Cluster computing environment
*[SIF]: Singularity Image Format - Container image format for HPC
*[TSV]: Tab-Separated Values - Plain text tabular data format
*[GBCMS]: Get Base Count Multi Sample - Multi-BAM variant counting tool
*[CLI]: Command Line Interface - Terminal-based program interface
*[CI]: Continuous Integration - Automated build and test system
*[CD]: Continuous Deployment - Automated release system
*[BAQ]: Base Alignment Quality - Heuristic quality downgrade near indels
*[UMI]: Unique Molecular Identifier - Molecular barcode for deduplication
*[PairHMM]: Pair Hidden Markov Model - Probabilistic alignment using base quality scores
*[SW]: Smith-Waterman - Edit-distance based local sequence alignment algorithm
*[LLR]: Log-Likelihood Ratio - Confidence metric for PairHMM alignment
*[ADAR]: Adenosine Deaminase Acting on RNA - Enzyme catalyzing A-to-I RNA editing
*[RT]: Reverse Transcriptase - Enzyme converting RNA to cDNA
*[NH]: Number of Hits - BAM tag indicating number of reported alignments
*[dUTP]: Deoxyuridine Triphosphate - Used in strand-specific RNA-seq protocols
*[STAR]: Spliced Transcripts Alignment to a Reference - RNA-seq aligner
*[BQSR]: Base Quality Score Recalibration - GATK quality recalibration tool
*[DP]: Total Depth - All mapped reads overlapping the variant position
*[RD]: Reference Depth - Reads supporting the reference allele
*[AD]: Alternate Depth - Reads supporting the alternate allele
*[DPF]: Fragment Depth - Total unique fragments at a variant position
*[RDF]: Reference Fragment Depth - Fragments resolved to reference allele
*[ADF]: Alternate Fragment Depth - Fragments resolved to alternate allele
*[mFSD]: Mutant Fragment Size Distribution - Insert-size analysis for REF vs ALT fragments
*[MSI]: Microsatellite Instability - Instability in short tandem repeat regions
*[MNP]: Multi-Nucleotide Polymorphism - Multiple adjacent base substitutions
*[SNP]: Single Nucleotide Polymorphism - Single base substitution
*[CIGAR]: Compact Idiosyncratic Gapped Alignment Report - BAM alignment operation string
*[WFA]: Wavefront Alignment - Edit-distance alignment algorithm used as Phase 3 fast-path triage (wfa2lib-rs); resolves ~70-80% of reads before PairHMM escalation
*[BWA]: Burrows-Wheeler Aligner - Standard short-read DNA aligner (mem/aln)
*[PyO3]: Python-Rust interop library for building native Python modules
