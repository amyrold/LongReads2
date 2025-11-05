# LongReads2 Workflow Documentation

## Pipeline Overview

LongReads2 analyzes 16S rRNA variation within and between bacterial strains using long-read sequencing data. The pipeline consists of seven main processes executed sequentially.

## Workflow Diagram

```
┌─────────────────────────┐
│  DOWNLOAD_METADATA      │  Download genome metadata from NCBI
│  Species → metadata.tsv │
└───────────┬─────────────┘
            │
            ▼
┌─────────────────────────┐
│  FILTER_LONGREADS       │  Filter for PacBio/Nanopore only
│  metadata → longreads   │
└───────────┬─────────────┘
            │
            ▼
┌─────────────────────────┐
│  DOWNLOAD_GENOMES       │  Download genome assemblies
│  accessions → genomes   │
└───────────┬─────────────┘
            │
            ▼
┌─────────────────────────┐
│  EXTRACT_16S            │  Extract 16S using Barrnap/BLAST
│  genomes → 16S seqs     │
└───────────┬─────────────┘
            │
            ▼
┌─────────────────────────┐
│  ALIGN_SEQUENCES        │  Multiple sequence alignment (MAFFT)
│  16S seqs → alignment   │
└───────────┬─────────────┘
            │
            ▼
┌─────────────────────────┐
│  CALCULATE_DISTANCES    │  Compute pairwise distances
│  alignment → distances  │
└───────────┬─────────────┘
            │
            ▼
┌─────────────────────────┐
│  VISUALIZE_RESULTS      │  Generate plots and statistics
│  distances → figures    │
└─────────────────────────┘
```

## Process Details

### 1. DOWNLOAD_METADATA
**Input**: Species name (string)
**Output**: metadata.tsv
**Container**: ncbi-tools

Downloads complete genome metadata from NCBI using the datasets API. Includes accession numbers and sequencing technology information.

**Key parameters**:
- Assembly level: complete
- Output fields: accession, sequencing technology

### 2. FILTER_LONGREADS
**Input**: metadata.tsv
**Output**: longreads.tsv
**Container**: ncbi-tools

Filters metadata to retain only genomes sequenced with long-read technologies (PacBio or Oxford Nanopore). Excludes hybrid assemblies that mix technologies.

**Filtering criteria**:
- Must contain "PacBio" OR "Nanopore"
- Must NOT contain: Illumina, semicolons, "and", "+"

### 3. DOWNLOAD_GENOMES
**Input**: longreads.tsv, max_genomes (optional)
**Output**: genomes.fasta
**Container**: ncbi-tools

Downloads genome assemblies from NCBI and combines them into a single multi-FASTA file.

**Key parameters**:
- `max_genomes`: Limit number of downloads (useful for testing)

### 4. EXTRACT_16S
**Input**: genomes.fasta, use_barrnap, min_length, max_length
**Output**: 16s_sequences.fasta, 16s_metadata.csv
**Container**: sequence-analysis

Extracts 16S rRNA sequences from genomes using either Barrnap (recommended) or BLAST.

**Barrnap mode (default)**:
- Species-agnostic detection using HMM models
- More sensitive than BLAST
- Proper strand handling
- Assigns unique IDs: `GENOME_16S_COPYNUMBER`

**BLAST mode (fallback)**:
- Uses E. coli 16S reference (J01859)
- Less sensitive but faster
- May miss divergent 16S sequences

**Length filtering**:
- Default: 1400-1600 bp
- Configurable via parameters

### 5. ALIGN_SEQUENCES
**Input**: 16s_sequences.fasta
**Output**: 16s_alignment.fasta, alignment_stats.txt
**Container**: sequence-analysis

Performs multiple sequence alignment using MAFFT. Algorithm selection is automatic based on sequence count:
- <200 sequences: L-INS-i (most accurate)
- 200-500 sequences: E-INS-i (balanced)
- >500 sequences: FFT-NS-2 (fastest)

**Why alignment is important**:
The original pipeline used raw edit distances, which don't account for insertions/deletions properly. Alignment-based distances are more biologically meaningful.

### 6. CALCULATE_DISTANCES
**Input**: 16s_alignment.fasta, 16s_metadata.csv
**Output**: within_distances.csv, between_distances.csv, distance_matrix.csv
**Container**: sequence-analysis

Calculates pairwise distances from the alignment using p-distance (proportion of differing sites).

**Distance categorization**:
- **Within-genome**: Comparisons between 16S copies from the same genome
- **Between-genome**: Comparisons between different genomes

**Output format**:
- `within_distances.csv`: seq1, seq2, genome, distance
- `between_distances.csv`: seq1, seq2, genome1, genome2, distance
- `distance_matrix.csv`: Full pairwise distance matrix

### 7. VISUALIZE_RESULTS
**Input**: within_distances.csv, between_distances.csv, distance_matrix.csv, metadata
**Output**: Multiple PNG plots, statistical_report.txt
**Container**: visualization

Generates comprehensive visualizations and statistical analyses:

**Plots generated**:
1. `distance_distribution.png` - Violin/box plot of distances
2. `distance_density.png` - Density plot showing distribution overlap
3. `copy_number_distribution.png` - Histogram of 16S copies per genome
4. `shared_unique_sequences.png` - Proportion of identical vs different sequences
5. `distance_heatmap.png` - Pairwise distance heatmap (if ≤100 sequences)
6. `mds_plot.png` - Multidimensional scaling visualization

**Statistical tests**:
- Wilcoxon rank-sum test (non-parametric)
- Welch's t-test (parametric)
- Effect size (Cohen's d)

## Data Flow

```
Species name
    ↓
Metadata TSV (all genomes)
    ↓
Filtered TSV (long-read only)
    ↓
Multi-FASTA (genomes)
    ↓
16S FASTA + metadata CSV
    ↓
Aligned 16S FASTA
    ↓
Distance matrices (within/between)
    ↓
Plots + statistical report
```

## Key Improvements Over Original

### Methodological Improvements
1. **Better 16S detection**: Barrnap uses HMM models vs single BLAST reference
2. **Proper alignment**: MAFFT alignment before distance calculation
3. **Better distances**: p-distance from alignment vs raw edit distance
4. **Flexible parameters**: Configurable length filters and detection methods

### Workflow Improvements
1. **Modular design**: Each step is a separate, reusable process
2. **Containerization**: Reproducible environments with Docker/Singularity
3. **Resource management**: Automatic scaling and retry logic
4. **Execution tracking**: Timeline, resource usage, and DAG reports

### Analysis Improvements
1. **Enhanced statistics**: Multiple tests, effect sizes
2. **Better visualizations**: 6 plots vs 3, publication-quality
3. **MDS analysis**: Explore sequence relationships in 2D space
4. **Comprehensive reporting**: Detailed summary statistics

## Resource Requirements

### Typical resource usage by process:

| Process | CPUs | Memory | Time |
|---------|------|--------|------|
| DOWNLOAD_METADATA | 1 | 2 GB | 10 min |
| FILTER_LONGREADS | 1 | 2 GB | 5 min |
| DOWNLOAD_GENOMES | 2 | 4 GB | 1-4 hours |
| EXTRACT_16S | 4 | 8 GB | 30 min - 2 hours |
| ALIGN_SEQUENCES | 8 | 16 GB | 30 min - 2 hours |
| CALCULATE_DISTANCES | 4 | 8 GB | 10 min - 1 hour |
| VISUALIZE_RESULTS | 2 | 4 GB | 5-15 min |

**Total pipeline time**: 2-8 hours (depending on genome count)

## Resuming Failed Runs

Nextflow caches completed steps. To resume after failure:

```bash
nextflow run main.nf -resume --species 'Escherichia coli' --outdir results
```

Only failed/incomplete processes will re-run.

## Debugging

View detailed logs:
```bash
# View work directory for a specific process
cat .nextflow/log

# Check process work directory
cd work/XX/XXXX...  # Hash from timeline report
cat .command.log    # Process stdout/stderr
```

## Common Issues

### No genomes found
- Check species name spelling
- Try broader taxonomy (e.g., "Escherichia" instead of "Escherichia coli")
- Some species may have no complete long-read assemblies

### Memory errors
- Increase memory: `--memory '32.GB'`
- Reduce genome count: `--max_genomes 20`

### Container pull failures
- Ensure Docker/Singularity is running
- Pre-pull containers: `docker pull longreads2/ncbi-tools:latest`

### No 16S found
- Try BLAST mode: `--use_barrnap false`
- Relax length filter: `--min_16s_length 1200`
