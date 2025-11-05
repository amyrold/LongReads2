/*
========================================================================================
    CALCULATE_DISTANCES: Calculate pairwise distances from alignment
========================================================================================
    Calculates pairwise distances using alignment-based methods instead of
    raw edit distances. This is more biologically meaningful and accounts for
    insertions/deletions properly.

    Improvements over original:
    - Uses p-distance (proportion of differing sites) from aligned sequences
    - Separates within-genome and between-genome comparisons
    - Calculates multiple distance metrics
    - Provides distance matrix for downstream analysis
========================================================================================
*/

process CALCULATE_DISTANCES {
    tag "Calculate distances"
    label 'process_medium'

    publishDir "${params.outdir}/3_output", mode: 'copy'

    input:
    path alignment
    path metadata

    output:
    path "within_distances.csv", emit: within_distances
    path "between_distances.csv", emit: between_distances
    path "distance_matrix.csv", emit: matrix
    path "distance_summary.txt", emit: summary

    script:
    """
    #!/usr/bin/env python3
    from Bio import AlignIO
    from Bio.Phylo.TreeConstruction import DistanceCalculator
    import pandas as pd
    import numpy as np
    import csv
    from itertools import combinations

    print("Loading alignment...")
    alignment = AlignIO.read('${alignment}', 'fasta')
    num_seqs = len(alignment)

    print(f"Calculating pairwise distances for {num_seqs} sequences...")

    # Load metadata to determine genome of origin for each sequence
    metadata_df = pd.read_csv('${metadata}')
    genome_map = {}
    for _, row in metadata_df.iterrows():
        genome_map[row['sequence_id']] = row['genome_accession']

    # Calculate distance matrix using p-distance (proportion of differences)
    # This is more appropriate than raw edit distance for aligned sequences
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)

    # Convert to full matrix
    seq_names = [record.id for record in alignment]
    n = len(seq_names)
    dist_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i != j:
                dist_matrix[i][j] = dm[seq_names[i], seq_names[j]]

    # Save distance matrix
    matrix_df = pd.DataFrame(dist_matrix, index=seq_names, columns=seq_names)
    matrix_df.to_csv('distance_matrix.csv')

    print("Separating within-genome and between-genome distances...")

    # Separate within-genome and between-genome comparisons
    within_distances = []
    between_distances = []

    for i, j in combinations(range(n), 2):
        seq1 = seq_names[i]
        seq2 = seq_names[j]
        distance = dist_matrix[i][j]

        # Get genome accessions (handle case where sequence not in metadata)
        genome1 = genome_map.get(seq1, seq1.split('_16S_')[0])
        genome2 = genome_map.get(seq2, seq2.split('_16S_')[0])

        if genome1 == genome2:
            # Within-genome comparison
            within_distances.append({
                'seq1': seq1,
                'seq2': seq2,
                'genome': genome1,
                'distance': distance
            })
        else:
            # Between-genome comparison
            between_distances.append({
                'seq1': seq1,
                'seq2': seq2,
                'genome1': genome1,
                'genome2': genome2,
                'distance': distance
            })

    # Save within and between distances
    if within_distances:
        within_df = pd.DataFrame(within_distances)
        within_df.to_csv('within_distances.csv', index=False)
    else:
        # Create empty file if no within-genome comparisons
        pd.DataFrame(columns=['seq1', 'seq2', 'genome', 'distance']).to_csv('within_distances.csv', index=False)

    if between_distances:
        between_df = pd.DataFrame(between_distances)
        between_df.to_csv('between_distances.csv', index=False)
    else:
        # Create empty file if no between-genome comparisons
        pd.DataFrame(columns=['seq1', 'seq2', 'genome1', 'genome2', 'distance']).to_csv('between_distances.csv', index=False)

    # Calculate summary statistics
    with open('distance_summary.txt', 'w') as f:
        f.write("Distance Calculation Summary\\n")
        f.write("============================\\n\\n")
        f.write(f"Total sequences: {num_seqs}\\n")
        f.write(f"Total pairwise comparisons: {len(within_distances) + len(between_distances)}\\n")
        f.write(f"Within-genome comparisons: {len(within_distances)}\\n")
        f.write(f"Between-genome comparisons: {len(between_distances)}\\n\\n")

        if within_distances:
            within_vals = [d['distance'] for d in within_distances]
            f.write(f"Within-genome distance statistics:\\n")
            f.write(f"  Mean: {np.mean(within_vals):.6f}\\n")
            f.write(f"  Median: {np.median(within_vals):.6f}\\n")
            f.write(f"  Std Dev: {np.std(within_vals):.6f}\\n")
            f.write(f"  Min: {np.min(within_vals):.6f}\\n")
            f.write(f"  Max: {np.max(within_vals):.6f}\\n\\n")

        if between_distances:
            between_vals = [d['distance'] for d in between_distances]
            f.write(f"Between-genome distance statistics:\\n")
            f.write(f"  Mean: {np.mean(between_vals):.6f}\\n")
            f.write(f"  Median: {np.median(between_vals):.6f}\\n")
            f.write(f"  Std Dev: {np.std(between_vals):.6f}\\n")
            f.write(f"  Min: {np.min(between_vals):.6f}\\n")
            f.write(f"  Max: {np.max(between_vals):.6f}\\n\\n")

        if within_distances and between_distances:
            f.write(f"Comparison:\\n")
            f.write(f"  Mean within-genome distance: {np.mean(within_vals):.6f}\\n")
            f.write(f"  Mean between-genome distance: {np.mean(between_vals):.6f}\\n")
            f.write(f"  Ratio (between/within): {np.mean(between_vals)/np.mean(within_vals):.2f}x\\n")

    print(f"Distance calculations complete!")
    print(f"  Within-genome: {len(within_distances)} comparisons")
    print(f"  Between-genome: {len(between_distances)} comparisons")
    """

    stub:
    """
    echo "seq1,seq2,genome,distance" > within_distances.csv
    echo "seq1,seq2,0.001" >> within_distances.csv
    echo "seq1,seq2,genome1,genome2,distance" > between_distances.csv
    echo "seq1,seq3,genome1,genome2,0.005" >> between_distances.csv
    echo "seq1,seq2" > distance_matrix.csv
    echo "Stub: Distance calculation complete" > distance_summary.txt
    """
}
