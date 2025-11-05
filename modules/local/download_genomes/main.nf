/*
========================================================================================
    DOWNLOAD_GENOMES: Download genome assemblies from NCBI
========================================================================================
    Downloads genome assemblies using NCBI datasets and combines them into a
    single multi-FASTA file. Optionally limits the number of genomes downloaded.
========================================================================================
*/

process DOWNLOAD_GENOMES {
    tag "Download genomes"
    label 'process_medium'

    publishDir "${params.outdir}/1_raw_data", mode: 'copy'

    input:
    path accessions
    val max_genomes

    output:
    path "genomes.fasta", emit: genomes
    path "download_log.txt", emit: log

    script:
    def limit = max_genomes ? "head -n ${max_genomes}" : "cat"
    """
    #!/bin/bash
    set -euo pipefail

    echo "Starting genome download..." | tee download_log.txt

    # Create empty output file
    touch genomes.fasta

    # Prepare accessions list (skip header, apply limit)
    # Using temp file to avoid SIGPIPE from broken pipes
    tail -n +2 ${accessions} > accessions_full.txt
    ${limit} accessions_full.txt > accessions_to_download.txt || true

    # Read accessions and download
    COUNTER=0
    while read -r line; do
        # Extract accession (first column)
        ACC=\$(echo "\$line" | cut -f1)

        echo "Downloading \$ACC..." | tee -a download_log.txt
        COUNTER=\$((COUNTER + 1))

        # Download genome using datasets
        datasets download genome accession \$ACC \\
            --include genome \\
            --filename \${ACC}.zip 2>> download_log.txt

        # Extract and append to multi-FASTA (overwrite files without prompting)
        unzip -o -q \${ACC}.zip
        cat ncbi_dataset/data/\${ACC}/*.fna >> genomes.fasta
        rm -rf \${ACC}.zip ncbi_dataset

        echo "  Downloaded genome \$COUNTER: \$ACC" | tee -a download_log.txt
    done < accessions_to_download.txt

    # Count total sequences
    SEQ_COUNT=\$(grep -c "^>" genomes.fasta || echo 0)
    echo "\\nTotal sequences downloaded: \$SEQ_COUNT" | tee -a download_log.txt

    if [ \$SEQ_COUNT -eq 0 ]; then
        echo "ERROR: No sequences downloaded" | tee -a download_log.txt
        exit 1
    fi
    """

    stub:
    """
    echo ">GCF_000005845.2_test_genome" > genomes.fasta
    echo "ATCGATCGATCGATCG" >> genomes.fasta
    echo "Stub: Downloaded 1 test genome" > download_log.txt
    """
}
