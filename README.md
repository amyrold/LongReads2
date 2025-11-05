# LongReads2

## Project Origin
This project is a complete architectural redesign of the original [LongReads](https://github.com/amyrold/LongReads) bioinformatics pipeline, originally developed in collaboration with Niru Shanbhag, Asad Shahzad, and Japani Doan. The original collaborative version remains available for reference. LongReads2 represents a seperate implementation using current best practices in bioinformatics workflow management.

## Description
LongReads2 is a Nextflow-based bioinformatics pipeline that compares the variation of the 16S rRNA region within and between strains of a specified bacterial species. Based on a user-specified genome, LongReads2 downloads metadata and filters to genomes sequenced using Long Read technology (PacBio and Oxford Nanopore). It then analyzes and visualizes the intra-genomic and inter-genomic variation of the different 16S regions seen in each genome.

## Major Improvements in LongReads2

### Architectural Improvements
- **Nextflow Pipeline Management**: Complete migration from Python scripts to Nextflow for robust workflow orchestration
- **Containerized Modules**: Each process runs in isolated Docker/Singularity containers for reproducibility
- **Simplified Installation**: Single installation method eliminates complexity of multiple deployment options
- **Modular Design**: Workflow broken into discrete, reusable processes
- **Enhanced Portability**: Runs consistently across different computing environments (local, HPC, cloud)
- **Improved Error Handling**: Automatic retry logic and comprehensive logging
- **Resource Management**: Dynamic resource allocation based on dataset size

### Methodological Improvements
- **Better 16S Detection**: Uses Barrnap (HMM-based) instead of single-reference BLAST for species-agnostic detection
- **Multiple Sequence Alignment**: Implements MAFFT alignment before distance calculation (original used raw edit distances)
- **Improved Distance Metrics**: Calculates alignment-based p-distances instead of raw edit distances
- **Flexible Parameters**: Configurable length filters, detection methods, and resource allocation
- **Enhanced Visualizations**: 6 publication-quality plots including dendrograms, heatmaps, and MDS plots
- **Comprehensive Statistics**: Multiple statistical tests (Wilcoxon, t-test, effect sizes)

## Prerequisites

- **Nextflow** (≥22.10.0): [Installation Guide](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- **Container Engine**: Docker or Singularity
- **Java** (≥11): Required for Nextflow

## Installation

LongReads2 uses a streamlined installation process with Nextflow handling all dependencies automatically.

1. **Install Nextflow**:
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

2. **Clone the repository**:
```bash
git clone https://github.com/[yourusername]/LongReads2
cd LongReads2
```

3. **Test the installation**:
```bash
nextflow run main.nf --help
```

That's it! Nextflow will automatically pull all required containers on first run.

## Quick Start

### First Run Test (5 genomes, ~30 minutes)

```bash
nextflow run main.nf \
  --species 'Escherichia coli' \
  --max_genomes 5 \
  --outdir test_results \
  -profile docker
```

### Standard Analysis (20 genomes)

```bash
nextflow run main.nf \
  --species 'Escherichia coli' \
  --max_genomes 20 \
  --outdir results
```

See [Quick Start Guide](docs/QUICKSTART.md) for detailed instructions.

## Container Options

LongReads2 supports multiple container engines:

**Docker (default)**:
```bash
nextflow run main.nf --species 'Escherichia coli' -profile docker
```

**Singularity (for HPC environments)**:
```bash
nextflow run main.nf --species 'Escherichia coli' -profile singularity
```

## Usage

### Basic Command Structure
```bash
nextflow run main.nf [options]
```

### Required Parameters
- `--species`: Bacterial species name (quoted if contains spaces)
- `--outdir`: Output directory path

### Optional Parameters
- `--max_genomes`: Maximum number of genomes to analyze (default: all available)
- `--min_length`: Minimum sequence length filter (default: 1000)
- `--blast_evalue`: BLAST e-value threshold (default: 1e-10)
- `--cpus`: Number of CPUs per process (default: 4)
- `--memory`: Memory allocation per process (default: '8.GB')

### Examples

**Basic analysis:**
```bash
nextflow run main.nf \
  --species 'Escherichia coli' \
  --outdir results
```

**Limited genome analysis with custom resources:**
```bash
nextflow run main.nf \
  --species 'Salmonella enterica' \
  --max_genomes 50 \
  --outdir salmonella_analysis \
  --cpus 8 \
  --memory '16.GB'
```

**Production run with all available genomes:**
```bash
nextflow run main.nf \
  --species 'Bacillus subtilis' \
  --outdir bacillus_complete \
  --cpus 16 \
  --memory '32.GB'
```

### Execution Profiles

**Local execution:**
```bash
nextflow run main.nf --species 'E. coli' -profile local
```

**HPC/SLURM execution:**
```bash
nextflow run main.nf --species 'E. coli' -profile slurm
```

**Cloud execution (AWS):**
```bash
nextflow run main.nf --species 'E. coli' -profile aws
```

## Pipeline Architecture

LongReads2 is built as a modular Nextflow pipeline with the following processes:

### Core Workflow Modules

1. **`DOWNLOAD_METADATA`** - Fetch species metadata from NCBI
2. **`FILTER_LONGREADS`** - Filter for long-read sequencing platforms
3. **`DOWNLOAD_GENOMES`** - Download filtered genome assemblies
4. **`EXTRACT_16S`** - BLAST-based 16S rRNA extraction
5. **`CALCULATE_DISTANCES`** - Compute edit distances within/between genomes
6. **`GENERATE_MATRIX`** - Create distance matrices
7. **`VISUALIZE_RESULTS`** - Generate plots and summary statistics

### Container Strategy
Each process runs in isolated containers:
- `longreads2/ncbi-tools`: NCBI datasets and BLAST
- `longreads2/biopython`: Sequence processing
- `longreads2/r-analysis`: Statistical analysis and visualization

## Output Structure

The pipeline generates four main directories:

### 1_raw_data
- `metadata.tsv` - Accession numbers and sequencing technology (pre-filtering)
- `longreads.tsv` - Accession numbers and sequencing technology (post-filtering)
- `wgs.fasta` - Multi-FASTA file of filtered genomes

### 2_filtered_data
- `trim.csv` - 16S rRNA BLAST results with copy numbers
- `trim.fasta` - Multi-FASTA file of 16S rRNA sequences

### 3_output
- `between.csv` - Edit distances between genomes
- `within.csv` - Edit distances within genomes
- `matrix.csv` - Complete edit distance matrix
- `output.txt` - Summary statistics
- `copies_per_genome.png` - Distribution of 16S copy numbers
- `inter_intra_shared_unique.png` - Shared vs. unique 16S sequences
- `boxplots.png` - Edit distance distributions

### 4_blast
- `16S.***` - BLAST database files
- `myresults.csv` - Raw BLAST results

## Dependencies

### Runtime Requirements
- **Nextflow** ≥22.10.0
- **Java** ≥11 (for Nextflow)
- **Container Engine**: Docker or Singularity

### Containerized Dependencies
All bioinformatics tools are automatically managed through containers:

**NCBI Tools Container:**
- ncbi-datasets CLI
- ncbi-dataformat CLI
- BLAST+ suite

**Sequence Analysis Container:**
- Python 3.9+
- BioPython
- pandas, numpy
- editdistance

**Visualization Container:**
- R 4.2+
- ggplot2, dplyr, readr
- Custom visualization scripts

## Configuration

### Resource Configuration
Customize resource allocation in `nextflow.config`:

```groovy
process {
    withName: 'DOWNLOAD_GENOMES' {
        cpus = 4
        memory = '8.GB'
        time = '2h'
    }
    withName: 'EXTRACT_16S' {
        cpus = 8
        memory = '16.GB'
        time = '4h'
    }
}
```

### Container Registry
By default, containers are pulled from Docker Hub. Configure custom registry:

```groovy
docker {
    registry = 'your-registry.com'
}
```

## Development

### Adding New Modules
1. Create module in `modules/local/`
2. Define process in `modules/local/newmodule/main.nf`
3. Add module to main workflow
4. Update container definitions

### Testing
Run the test dataset:
```bash
nextflow run main.nf -profile test
```

### Building Custom Containers
See [containers/README.md](containers/README.md) for instructions on building containers locally.

## Documentation

### Learning Resources
- **[Quick Start Guide](docs/QUICKSTART.md)** - Get started in 10 minutes
- **[Workflow Documentation](docs/WORKFLOW.md)** - Detailed pipeline architecture and methodology
- **[Container Guide](containers/README.md)** - Building and managing Docker/Singularity containers

### Understanding the Pipeline
The pipeline consists of 7 main processes:
1. **DOWNLOAD_METADATA** - Fetch genome metadata from NCBI
2. **FILTER_LONGREADS** - Filter for PacBio/Nanopore sequencing
3. **DOWNLOAD_GENOMES** - Download genome assemblies
4. **EXTRACT_16S** - Extract 16S rRNA using Barrnap or BLAST
5. **ALIGN_SEQUENCES** - Multiple sequence alignment with MAFFT
6. **CALCULATE_DISTANCES** - Compute pairwise distances
7. **VISUALIZE_RESULTS** - Generate plots and statistical reports

See [WORKFLOW.md](docs/WORKFLOW.md) for detailed process descriptions.

## Troubleshooting

### Common Issues

**Container pull failures:**
```bash
# Pre-pull containers
docker pull longreads2/ncbi-tools:latest
docker pull longreads2/sequence-analysis:latest
docker pull longreads2/visualization:latest
```

**Memory issues:**
```bash
# Increase memory allocation
nextflow run main.nf --memory '32.GB'
```

**Resume failed runs:**
```bash
nextflow run main.nf -resume --species 'Your Species' --outdir results
```

**No 16S sequences found:**
```bash
# Try BLAST mode instead of Barrnap
nextflow run main.nf --species 'Your Species' --use_barrnap false

# Or relax length filters
nextflow run main.nf --species 'Your Species' --min_16s_length 1200
```

For more troubleshooting tips, see the [Quick Start Guide](docs/QUICKSTART.md).

## Contributing

This is a personal portfolio project, but feedback and suggestions are welcome through issues.

## Citation

If you use this tool in your research, please cite the original collaborative work and acknowledge this extended version:

```
Original LongReads pipeline developed by Aaron Myrold, Niru Shanbhag, Asad Shahzad, and Japani Doan.
Extended version available at: https://github.com/[yourusername]/LongReads-Extended
```

## License

[Add your chosen license here]

## Contact

[Your contact information]

---

**Original Collaborative Team:** Aaron Myrold, Niru Shanbhag, Asad Shahzad, Japani Doan