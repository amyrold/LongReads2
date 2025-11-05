/*
========================================================================================
    FILTER_LONGREADS: Filter metadata for long-read sequencing technologies
========================================================================================
    Filters the metadata to keep only genomes sequenced with long-read technologies
    (PacBio and Oxford Nanopore). Excludes hybrid assemblies that mix technologies.
========================================================================================
*/

process FILTER_LONGREADS {
    tag "Filter for long reads"
    label 'process_low'

    publishDir "${params.outdir}/1_raw_data", mode: 'copy'

    input:
    path metadata

    output:
    path "longreads.tsv", emit: accessions

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import sys

    # Read metadata
    df = pd.read_csv('${metadata}', sep='\\t')
    initial_count = len(df)
    print(f"Initial genome count: {initial_count}")

    # Filter for long-read technologies
    # Keep only PacBio or Nanopore, exclude hybrid assemblies
    exclude_terms = [';', 'and', ',', '+', 'Illumina', 'Ion Torrent']

    longread_genomes = []

    for idx, row in df.iterrows():
        acc = row['Assembly Accession']
        tech = str(row['Assembly Sequencing Tech'])

        # Check if it's pure PacBio or Nanopore
        if ('PacBio' in tech or 'Nanopore' in tech):
            # Make sure it doesn't contain excluded terms
            if not any(term in tech for term in exclude_terms):
                longread_genomes.append({
                    'Accession': acc,
                    'Sequencing_Technology': tech
                })

    # Create filtered dataframe
    filtered_df = pd.DataFrame(longread_genomes)

    if len(filtered_df) == 0:
        print("ERROR: No long-read genomes found after filtering")
        sys.exit(1)

    # Save to file
    filtered_df.to_csv('longreads.tsv', sep='\\t', index=False)

    # Print summary
    filtered_count = len(filtered_df)
    pacbio_count = sum('PacBio' in tech for tech in filtered_df['Sequencing_Technology'])
    nanopore_count = sum('Nanopore' in tech for tech in filtered_df['Sequencing_Technology'])

    print(f"\\nFiltering Results:")
    print(f"  Total genomes after filtering: {filtered_count}")
    print(f"  PacBio: {pacbio_count}")
    print(f"  Oxford Nanopore: {nanopore_count}")
    print(f"  Retention rate: {filtered_count/initial_count*100:.1f}%")
    """

    stub:
    """
    echo -e "Accession\\tSequencing_Technology" > longreads.tsv
    echo -e "GCF_000005845.2\\tPacBio" >> longreads.tsv
    echo -e "GCF_000008865.2\\tOxford Nanopore" >> longreads.tsv
    """
}
