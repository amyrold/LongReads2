/*
========================================================================================
    ALIGN_SEQUENCES: Multiple sequence alignment of 16S sequences
========================================================================================
    Aligns 16S rRNA sequences using MAFFT. This is a major improvement over the
    original pipeline which used raw edit distances.

    Benefits:
    - Proper handling of insertions/deletions
    - More accurate distance calculations
    - Enables phylogenetic analysis
    - Better handling of sequence variation
========================================================================================
*/

process ALIGN_SEQUENCES {
    tag "Align 16S sequences"
    label 'process_high'

    publishDir "${params.outdir}/2_filtered_data", mode: 'copy'

    input:
    path sequences

    output:
    path "16s_alignment.fasta", emit: alignment
    path "alignment_stats.txt", emit: stats

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "Aligning 16S sequences with MAFFT..."

    # Count sequences
    SEQ_COUNT=\$(grep -c "^>" ${sequences})
    echo "Number of sequences to align: \$SEQ_COUNT"

    if [ \$SEQ_COUNT -lt 2 ]; then
        echo "ERROR: Need at least 2 sequences for alignment"
        exit 1
    fi

    # Run MAFFT with appropriate algorithm based on sequence count
    # For <200 sequences: use L-INS-i (most accurate)
    # For 200-500 sequences: use E-INS-i
    # For >500 sequences: use FFT-NS-2 (faster)

    if [ \$SEQ_COUNT -lt 200 ]; then
        echo "Using MAFFT L-INS-i (accurate, for <200 sequences)..."
        mafft --localpair --maxiterate 1000 \\
            --thread ${task.cpus} \\
            ${sequences} > 16s_alignment.fasta
    elif [ \$SEQ_COUNT -lt 500 ]; then
        echo "Using MAFFT E-INS-i (for 200-500 sequences)..."
        mafft --ep 0 --genafpair --maxiterate 1000 \\
            --thread ${task.cpus} \\
            ${sequences} > 16s_alignment.fasta
    else
        echo "Using MAFFT FFT-NS-2 (fast, for >500 sequences)..."
        mafft --retree 2 --maxiterate 0 \\
            --thread ${task.cpus} \\
            ${sequences} > 16s_alignment.fasta
    fi

    echo "Alignment complete!"

    # Calculate alignment statistics
    python3 << 'EOF'
from Bio import AlignIO
from Bio.Align import AlignInfo
import sys

# Read alignment
alignment = AlignIO.read('16s_alignment.fasta', 'fasta')

# Calculate statistics
num_seqs = len(alignment)
align_length = alignment.get_alignment_length()

# Calculate conservation
summary = AlignInfo.SummaryInfo(alignment)

# Count gaps
total_gaps = 0
gap_positions = []
for i in range(align_length):
    column = alignment[:, i]
    gaps = column.count('-')
    total_gaps += gaps
    if gaps > 0:
        gap_positions.append((i, gaps))

gap_percentage = (total_gaps / (num_seqs * align_length)) * 100

# Calculate pairwise identity
identities = []
for i in range(num_seqs):
    for j in range(i+1, num_seqs):
        seq1 = str(alignment[i].seq)
        seq2 = str(alignment[j].seq)
        matches = sum(a == b and a != '-' for a, b in zip(seq1, seq2))
        valid_positions = sum(a != '-' and b != '-' for a, b in zip(seq1, seq2))
        if valid_positions > 0:
            identity = (matches / valid_positions) * 100
            identities.append(identity)

avg_identity = sum(identities) / len(identities) if identities else 0
min_identity = min(identities) if identities else 0
max_identity = max(identities) if identities else 0

# Write statistics
with open('alignment_stats.txt', 'w') as f:
    f.write(f"Multiple Sequence Alignment Statistics\\n")
    f.write(f"======================================\\n\\n")
    f.write(f"Number of sequences: {num_seqs}\\n")
    f.write(f"Alignment length: {align_length} bp\\n")
    f.write(f"Total gaps: {total_gaps} ({gap_percentage:.2f}%)\\n")
    f.write(f"Number of positions with gaps: {len(gap_positions)}\\n\\n")
    f.write(f"Pairwise Identity Statistics:\\n")
    f.write(f"  Average: {avg_identity:.2f}%\\n")
    f.write(f"  Minimum: {min_identity:.2f}%\\n")
    f.write(f"  Maximum: {max_identity:.2f}%\\n\\n")
    f.write(f"Gap-heavy positions (>50% gaps):\\n")
    for pos, gaps in gap_positions:
        if gaps > num_seqs * 0.5:
            f.write(f"  Position {pos}: {gaps}/{num_seqs} gaps ({gaps/num_seqs*100:.1f}%)\\n")

print(f"Alignment statistics written to alignment_stats.txt")
print(f"Average pairwise identity: {avg_identity:.2f}%")
EOF

    echo "Alignment statistics calculated!"
    """

    stub:
    """
    echo ">seq1" > 16s_alignment.fasta
    echo "ATCG-ATCG" >> 16s_alignment.fasta
    echo ">seq2" >> 16s_alignment.fasta
    echo "ATCGGATCG" >> 16s_alignment.fasta
    echo "Stub: Alignment complete" > alignment_stats.txt
    """
}
