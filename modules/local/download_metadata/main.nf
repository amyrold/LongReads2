/*
========================================================================================
    DOWNLOAD_METADATA: Download genome metadata from NCBI
========================================================================================
    Downloads metadata for all complete genomes of a specified bacterial species
    from NCBI using the datasets CLI tool.
========================================================================================
*/

process DOWNLOAD_METADATA {
    tag "$species"
    label 'process_low'

    publishDir "${params.outdir}/1_raw_data", mode: 'copy'

    input:
    val species

    output:
    path "metadata.tsv", emit: metadata

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "Downloading metadata for ${species}..."

    # Download genome metadata using NCBI datasets
    # Filter for complete genomes only
    datasets summary genome taxon '${species}' \\
        --assembly-level complete \\
        --as-json-lines | \\
    dataformat tsv genome \\
        --fields accession,assminfo-sequencing-tech \\
        > metadata.tsv

    # Count how many genomes were found
    GENOME_COUNT=\$(tail -n +2 metadata.tsv | wc -l | tr -d ' ')
    echo "Found \${GENOME_COUNT} complete ${species} genome assemblies"

    # Check if we found any genomes
    if [ \$GENOME_COUNT -eq 0 ]; then
        echo "ERROR: No complete genomes found for ${species}"
        exit 1
    fi
    """

    stub:
    """
    echo -e "Assembly Accession\\tAssembly Sequencing Tech" > metadata.tsv
    echo -e "GCF_000005845.2\\tPacBio" >> metadata.tsv
    echo -e "GCF_000008865.2\\tOxford Nanopore" >> metadata.tsv
    """
}
