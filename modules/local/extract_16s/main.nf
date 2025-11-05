/*
========================================================================================
    EXTRACT_16S: Extract 16S rRNA sequences from genomes
========================================================================================
    Extracts 16S rRNA sequences using either Barrnap (recommended) or BLAST.
    Barrnap is more sensitive and doesn't require a reference sequence.

    Improvements over original:
    - Uses Barrnap for better detection (species-agnostic)
    - Proper strand handling
    - Length filtering with flexible ranges
    - Quality control and metadata tracking
========================================================================================
*/

process EXTRACT_16S {
    tag "Extract 16S sequences"
    label 'process_medium'

    publishDir "${params.outdir}/2_filtered_data", mode: 'copy'

    input:
    path genomes
    val use_barrnap
    val min_length
    val max_length

    output:
    path "16s_sequences.fasta", emit: sequences
    path "16s_metadata.csv", emit: metadata
    path "extraction_log.txt", emit: log

    script:
    if (use_barrnap)
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    import subprocess
    import csv
    import sys

    print("Extracting 16S sequences using Barrnap...")

    # Run Barrnap to identify 16S rRNA genes
    # Barrnap is more accurate than BLAST for rRNA detection
    print("Running Barrnap...")
    barrnap_cmd = [
        'barrnap',
        '--kingdom', 'bac',  # bacterial
        '--threads', '${task.cpus}',
        '${genomes}'
    ]

    with open('barrnap_output.gff', 'w') as f:
        subprocess.run(barrnap_cmd, stdout=f, check=True)

    print("Barrnap complete. Parsing results...")

    # Parse Barrnap GFF output and extract sequences
    genome_dict = SeqIO.to_dict(SeqIO.parse('${genomes}', 'fasta'))

    extracted_sequences = []
    metadata_records = []
    copy_counter = {}  # Track copy number per genome

    with open('barrnap_output.gff', 'r') as gff:
        for line in gff:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\\t')
            if len(fields) < 9:
                continue

            # Parse GFF fields
            seq_id = fields[0]
            feature = fields[2]
            start = int(fields[3]) - 1  # Convert to 0-based
            end = int(fields[4])
            score = float(fields[5]) if fields[5] != '.' else 0
            strand = fields[6]
            attributes = fields[8]

            # Only process 16S_rRNA features (check attributes field)
            if '16S_rRNA' not in attributes:
                continue

            # Get sequence
            if seq_id not in genome_dict:
                # Try without version number
                seq_id_base = seq_id.split('.')[0]
                matching = [k for k in genome_dict.keys() if k.startswith(seq_id_base)]
                if matching:
                    seq_id = matching[0]
                else:
                    continue

            sequence = genome_dict[seq_id].seq[start:end]

            # Handle reverse strand
            if strand == '-':
                sequence = sequence.reverse_complement()

            seq_length = len(sequence)

            # Apply length filter
            if seq_length < ${min_length} or seq_length > ${max_length}:
                continue

            # Generate unique ID for this 16S copy
            genome_acc = seq_id.split()[0]  # Get accession without description
            if genome_acc not in copy_counter:
                copy_counter[genome_acc] = 0
            copy_counter[genome_acc] += 1

            seq_id_unique = f"{genome_acc}_16S_{copy_counter[genome_acc]}"

            # Create SeqRecord
            record = SeqRecord(
                sequence,
                id=seq_id_unique,
                description=f"length={seq_length} strand={strand} score={score:.2f}"
            )
            extracted_sequences.append(record)

            # Store metadata
            metadata_records.append({
                'sequence_id': seq_id_unique,
                'genome_accession': genome_acc,
                'copy_number': copy_counter[genome_acc],
                'start': start,
                'end': end,
                'length': seq_length,
                'strand': strand,
                'score': score
            })

    # Write sequences to FASTA
    with open('16s_sequences.fasta', 'w') as f:
        SeqIO.write(extracted_sequences, f, 'fasta')

    # Write metadata to CSV
    with open('16s_metadata.csv', 'w', newline='') as f:
        if metadata_records:
            writer = csv.DictWriter(f, fieldnames=metadata_records[0].keys())
            writer.writeheader()
            writer.writerows(metadata_records)

    # Write log
    with open('extraction_log.txt', 'w') as f:
        total_16s = len(extracted_sequences)
        num_genomes = len(copy_counter)
        avg_copies = total_16s / num_genomes if num_genomes > 0 else 0

        f.write(f"16S Extraction Summary\\n")
        f.write(f"=====================\\n")
        f.write(f"Method: Barrnap\\n")
        f.write(f"Total 16S sequences extracted: {total_16s}\\n")
        f.write(f"Number of genomes: {num_genomes}\\n")
        f.write(f"Average copies per genome: {avg_copies:.2f}\\n")
        f.write(f"Length filter: {${min_length}}-${max_length} bp\\n\\n")

        f.write(f"Copy number distribution:\\n")
        for genome, copies in sorted(copy_counter.items()):
            f.write(f"  {genome}: {copies} copies\\n")

    print(f"Extracted {total_16s} 16S sequences from {num_genomes} genomes")

    if total_16s == 0:
        print("ERROR: No 16S sequences found")
        sys.exit(1)
    """
    else
    """
    #!/bin/bash
    set -euo pipefail

    echo "Extracting 16S sequences using BLAST..."

    # Create BLAST database
    makeblastdb -in ${genomes} -dbtype nucl -out genome_db

    # Download 16S reference (E. coli)
    echo "Downloading 16S reference..."
    esearch -db nucleotide -query "J01859" | efetch -format fasta > 16s_ref.fasta

    # Run BLAST
    blastn -query 16s_ref.fasta \\
        -db genome_db \\
        -out blast_results.txt \\
        -outfmt "6 qacc sacc pident sstart send length sseq sstrand evalue" \\
        -evalue ${params.blast_evalue} \\
        -num_threads ${task.cpus}

    # Process BLAST results with Python
    python3 << 'EOF'
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv

# Parse BLAST results
blast_hits = []
with open('blast_results.txt', 'r') as f:
    for line in f:
        fields = line.strip().split('\\t')
        if len(fields) >= 9:
            blast_hits.append({
                'genome': fields[1],
                'start': int(fields[3]),
                'end': int(fields[4]),
                'length': int(fields[5]),
                'strand': fields[7],
                'evalue': float(fields[8])
            })

# Load genomes
genome_dict = SeqIO.to_dict(SeqIO.parse('${genomes}', 'fasta'))

# Extract sequences
extracted = []
metadata = []
copy_counter = {}

for hit in blast_hits:
    genome_id = hit['genome']
    start = min(hit['start'], hit['end']) - 1
    end = max(hit['start'], hit['end'])

    if genome_id not in genome_dict:
        continue

    seq = genome_dict[genome_id].seq[start:end]
    if hit['strand'] == 'minus':
        seq = seq.reverse_complement()

    # Length filter
    if len(seq) < ${min_length} or len(seq) > ${max_length}:
        continue

    # Unique ID
    if genome_id not in copy_counter:
        copy_counter[genome_id] = 0
    copy_counter[genome_id] += 1

    seq_id = f"{genome_id}_16S_{copy_counter[genome_id]}"

    record = SeqRecord(seq, id=seq_id, description="")
    extracted.append(record)

    metadata.append({
        'sequence_id': seq_id,
        'genome_accession': genome_id,
        'copy_number': copy_counter[genome_id],
        'length': len(seq),
        'strand': hit['strand']
    })

# Write outputs
SeqIO.write(extracted, '16s_sequences.fasta', 'fasta')

with open('16s_metadata.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['sequence_id', 'genome_accession', 'copy_number', 'length', 'strand'])
    writer.writeheader()
    writer.writerows(metadata)

with open('extraction_log.txt', 'w') as f:
    f.write(f"Method: BLAST\\n")
    f.write(f"Total 16S sequences: {len(extracted)}\\n")
    f.write(f"Number of genomes: {len(copy_counter)}\\n")
EOF
    """

    stub:
    """
    echo ">test_genome_1_16S_1" > 16s_sequences.fasta
    echo "ATCGATCGATCG" >> 16s_sequences.fasta
    echo "sequence_id,genome_accession,copy_number,length,strand" > 16s_metadata.csv
    echo "test_genome_1_16S_1,test_genome_1,1,1500,+" >> 16s_metadata.csv
    echo "Stub: Extracted 1 16S sequence" > extraction_log.txt
    """
}
