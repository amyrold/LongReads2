#!/usr/bin/env nextflow

/*
========================================================================================
    LongReads2: 16S rRNA Variation Analysis Pipeline
========================================================================================
    Github : https://github.com/amyrold/LongReads2
    Author : Aaron Myrold
----------------------------------------------------------------------------------------
    Description:
    Nextflow pipeline to analyze within-strain and between-strain variation in
    bacterial 16S rRNA genes using long-read sequencing data.
========================----------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*
========================================================================================
    WORKFLOW PARAMETERS
========================================================================================
*/

// Input parameters
params.species = null
params.max_genomes = null
params.outdir = 'results'

// Process parameters
params.min_16s_length = 1400
params.max_16s_length = 1600
params.blast_evalue = 1e-10
params.use_barrnap = true  // Use Barrnap for 16S detection (faster, more accurate)

// Resource parameters
params.cpus = 4
params.memory = '8.GB'

// Help message
params.help = false

def helpMessage() {
    log.info"""
    =========================================
    LongReads2 Pipeline v1.0
    =========================================

    Usage:
    nextflow run main.nf --species 'Escherichia coli' --outdir results

    Required arguments:
      --species         Bacterial species name (quoted if contains spaces)
      --outdir          Output directory path

    Optional arguments:
      --max_genomes     Maximum number of genomes to analyze (default: all)
      --min_16s_length  Minimum 16S length in bp (default: 1400)
      --max_16s_length  Maximum 16S length in bp (default: 1600)
      --use_barrnap     Use Barrnap for 16S detection (default: true)
      --cpus            Number of CPUs per process (default: 4)
      --memory          Memory per process (default: 8.GB)

    Profiles:
      -profile docker       Use Docker containers (default)
      -profile singularity  Use Singularity containers (for HPC)
      -profile test         Run test dataset

    Example:
      nextflow run main.nf \\
        --species 'Escherichia coli' \\
        --max_genomes 10 \\
        --outdir results \\
        -profile docker
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.species) {
    log.error "ERROR: --species parameter is required"
    helpMessage()
    exit 1
}

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { DOWNLOAD_METADATA } from './modules/local/download_metadata/main.nf'
include { FILTER_LONGREADS } from './modules/local/filter_longreads/main.nf'
include { DOWNLOAD_GENOMES } from './modules/local/download_genomes/main.nf'
include { EXTRACT_16S } from './modules/local/extract_16s/main.nf'
include { ALIGN_SEQUENCES } from './modules/local/align_sequences/main.nf'
include { CALCULATE_DISTANCES } from './modules/local/calculate_distances/main.nf'
include { VISUALIZE_RESULTS } from './modules/local/visualize_results/main.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    // Print pipeline information
    log.info """
    =========================================
    LongReads2 Pipeline
    =========================================
    Species        : ${params.species}
    Max genomes    : ${params.max_genomes ?: 'all'}
    Output dir     : ${params.outdir}
    16S detection  : ${params.use_barrnap ? 'Barrnap' : 'BLAST'}
    Min 16S length : ${params.min_16s_length} bp
    Max 16S length : ${params.max_16s_length} bp
    =========================================
    """.stripIndent()

    //
    // STEP 1: Download metadata for species from NCBI
    //
    DOWNLOAD_METADATA(
        params.species
    )

    //
    // STEP 2: Filter metadata for long-read sequencing technologies
    //
    FILTER_LONGREADS(
        DOWNLOAD_METADATA.out.metadata
    )

    //
    // STEP 3: Download genome assemblies
    //
    DOWNLOAD_GENOMES(
        FILTER_LONGREADS.out.accessions,
        params.max_genomes
    )

    //
    // STEP 4: Extract 16S rRNA sequences
    //
    EXTRACT_16S(
        DOWNLOAD_GENOMES.out.genomes,
        params.use_barrnap,
        params.min_16s_length,
        params.max_16s_length
    )

    //
    // STEP 5: Align 16S sequences with MAFFT
    //
    ALIGN_SEQUENCES(
        EXTRACT_16S.out.sequences
    )

    //
    // STEP 6: Calculate distance matrices
    //
    CALCULATE_DISTANCES(
        ALIGN_SEQUENCES.out.alignment,
        EXTRACT_16S.out.metadata
    )

    //
    // STEP 7: Generate visualizations and statistics
    //
    VISUALIZE_RESULTS(
        CALCULATE_DISTANCES.out.within_distances,
        CALCULATE_DISTANCES.out.between_distances,
        CALCULATE_DISTANCES.out.matrix,
        EXTRACT_16S.out.metadata
    )

}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """
    =========================================
    Pipeline completed successfully!
    Results saved to: ${params.outdir}
    =========================================
    """.stripIndent()
}

workflow.onError {
    log.error """
    =========================================
    Pipeline execution failed!
    Error message: ${workflow.errorMessage}
    =========================================
    """.stripIndent()
}
