# LongReads2

## Project Origin
This project is a complete architectural redesign of the original [LongReads](https://github.com/amyrold/LongReads) bioinformatics pipeline, originally developed in collaboration with Aaron Myrold, Niru Shanbhag, Asad Shahzad, and Japani Doan. The original collaborative version remains available for reference. LongReads2 represents a modern reimplementation using current best practices in bioinformatics workflow management.

## Description
LongReads2 is a Nextflow-based bioinformatics pipeline that compares the variation of the 16S rRNA region within and between strains of a specified bacterial species. Based on a user-specified genome, LongReads2 downloads metadata and filters to genomes sequenced using Long Read technology (PacBio and Oxford Nanopore). It then analyzes and visualizes the intra-genomic and inter-genomic variation of the different 16S regions seen in each genome.

## Major Architectural Improvements in LongReads2
- **Nextflow Pipeline Management**: Complete migration from Python scripts to Nextflow for robust workflow orchestration
- **Containerized Modules**: Each process runs in isolated Docker/Singularity containers for reproducibility
- **Simplified Installation**: Single installation method eliminates complexity of multiple deployment options
- **Modular Design**: Workflow broken into discrete, reusable processes
- **Enhanced Portability**: Runs consistently across different computing environments (local, HPC, cloud)
- **Improved Error Handling**: Automatic retry logic and comprehensive logging
- **Resource Management**: Dynamic resource allocation based on dataset size

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

Run a test analysis with 10 *E. coli* genomes:

```bash
nextflow run main.nf \
  --species 'Escherichia coli' \
  --max_genomes 10 \
  --outdir results
```

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

## Troubleshooting

### Common Issues

**Container pull failures:**
```bash
# Pre-pull containers
docker pull longreads2/ncbi-tools:latest
docker pull longreads2/biopython:latest
docker pull longreads2/r-analysis:latest
```

**Memory issues:**
```bash
# Increase memory allocation
nextflow run main.nf --max_memory '64.GB'
```

**Resume failed runs:**
```bash
nextflow run main.nf -resume
```

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