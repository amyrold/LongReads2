# LongReads2 Quick Start Guide

This guide will get you up and running with LongReads2 in under 10 minutes.

## Prerequisites

You need three things installed:
1. **Nextflow** (≥22.10.0)
2. **Java** (≥11) - required by Nextflow
3. **Docker** OR **Singularity** - for containers

## Installation

### Step 1: Install Nextflow

```bash
# Download and install Nextflow
curl -s https://get.nextflow.io | bash

# Move to a directory in your PATH
sudo mv nextflow /usr/local/bin/
# OR add to local bin
mkdir -p ~/bin && mv nextflow ~/bin/

# Verify installation
nextflow -version
```

### Step 2: Clone the Repository

```bash
git clone https://github.com/amyrold/LongReads2
cd LongReads2
```

### Step 3: Build Docker Containers (Optional)

If you want to use your own containers:

```bash
cd containers
./build_containers.sh
cd ..
```

Otherwise, the pipeline will pull pre-built containers automatically.

## Running Your First Analysis

### Test Run (5 genomes, ~30 minutes)

```bash
nextflow run main.nf \
  --species 'Escherichia coli' \
  --max_genomes 5 \
  --outdir test_results \
  -profile docker
```

This will:
1. Download 5 E. coli genomes
2. Extract 16S sequences
3. Align them with MAFFT
4. Calculate distances
5. Generate plots and statistics

### Small Analysis (10-20 genomes, 1-2 hours)

```bash
nextflow run main.nf \
  --species 'Salmonella enterica' \
  --max_genomes 20 \
  --outdir salmonella_results
```

### Full Analysis (all available genomes)

```bash
nextflow run main.nf \
  --species 'Bacillus subtilis' \
  --outdir bacillus_full \
  --cpus 8 \
  --memory '16.GB'
```

## Understanding the Output

After the pipeline completes, you'll find:

```
results/
├── 1_raw_data/
│   ├── metadata.tsv          # All genomes (before filtering)
│   ├── longreads.tsv          # Long-read genomes only
│   └── genomes.fasta          # Downloaded genomes
├── 2_filtered_data/
│   ├── 16s_sequences.fasta    # Extracted 16S sequences
│   ├── 16s_metadata.csv       # 16S copy information
│   ├── 16s_alignment.fasta    # Aligned sequences
│   └── alignment_stats.txt    # Alignment quality metrics
├── 3_output/
│   ├── within_distances.csv   # Intra-genomic distances
│   ├── between_distances.csv  # Inter-genomic distances
│   ├── distance_matrix.csv    # Full distance matrix
│   ├── distance_summary.txt   # Statistical summary
│   ├── *.png                  # Visualization plots
│   └── statistical_report.txt # Comprehensive stats
└── pipeline_info/
    ├── execution_timeline.html # When each process ran
    ├── execution_report.html   # Resource usage
    └── pipeline_dag.svg        # Workflow diagram
```

### Key Files to Check

**Quick summary**:
- `3_output/statistical_report.txt` - Main results
- `3_output/distance_distribution.png` - Main figure

**Detailed analysis**:
- `3_output/distance_matrix.csv` - For custom analyses
- `2_filtered_data/16s_alignment.fasta` - For phylogenetic trees
- `pipeline_info/execution_report.html` - Performance metrics

## Common Customizations

### Using Different Species

```bash
nextflow run main.nf \
  --species 'Mycobacterium tuberculosis' \
  --outdir myco_results
```

Note: Species names should match NCBI taxonomy. Use quotes if the name has spaces.

### Using BLAST Instead of Barrnap

```bash
nextflow run main.nf \
  --species 'Escherichia coli' \
  --use_barrnap false \
  --outdir results
```

### Relaxing Length Filters

If you're not finding many 16S sequences:

```bash
nextflow run main.nf \
  --species 'Your Species' \
  --min_16s_length 1200 \
  --max_16s_length 1700 \
  --outdir results
```

### Running on HPC with SLURM

```bash
nextflow run main.nf \
  --species 'Escherichia coli' \
  --outdir results \
  -profile slurm,singularity
```

Edit `nextflow.config` to set your SLURM account:
```groovy
slurm {
    process {
        queue = 'your_partition'
        clusterOptions = '--account=your_account'
    }
}
```

## Resuming Failed Runs

If the pipeline fails or is interrupted:

```bash
nextflow run main.nf -resume \
  --species 'Escherichia coli' \
  --outdir results
```

Nextflow will skip completed steps and only re-run failed processes.

## Getting Help

### View all parameters
```bash
nextflow run main.nf --help
```

### Check pipeline execution
```bash
# View last run log
cat .nextflow.log

# View detailed timeline
open results/pipeline_info/execution_timeline.html
```

### Troubleshooting

**Issue**: "No genomes found"
- Check species name spelling
- Try: `--species 'Escherichia'` instead of `--species 'Escherichia coli'`

**Issue**: "No 16S sequences found"
- Try BLAST mode: `--use_barrnap false`
- Relax filters: `--min_16s_length 1200`

**Issue**: Out of memory
- Increase: `--memory '32.GB'`
- Reduce genomes: `--max_genomes 10`

**Issue**: Docker permission denied
- Add yourself to docker group: `sudo usermod -aG docker $USER`
- Log out and back in

## Next Steps

1. **Explore the visualizations** in `3_output/`
2. **Read the workflow documentation**: `docs/WORKFLOW.md`
3. **Customize the analysis** for your specific needs
4. **Compare multiple species** by running the pipeline multiple times

## Example: Comparing E. coli Strains

```bash
# Run analysis
nextflow run main.nf \
  --species 'Escherichia coli' \
  --max_genomes 30 \
  --outdir ecoli_analysis

# Check results
cat ecoli_analysis/3_output/statistical_report.txt

# View main figure
open ecoli_analysis/3_output/distance_distribution.png
```

Look for:
- **Mean within-genome distance**: How variable are 16S copies within a single genome?
- **Mean between-genome distance**: How different are 16S sequences between strains?
- **Statistical significance**: Is there a significant difference?

## Tips for Success

1. **Start small**: Use `--max_genomes 5` for testing
2. **Monitor progress**: Check `tail -f .nextflow.log`
3. **Check resources**: View `results/pipeline_info/execution_report.html`
4. **Save your commands**: Keep a record of parameters used
5. **Use resume**: Don't restart from scratch after failures

## Advanced: Custom Profiles

Create a custom profile in `nextflow.config`:

```groovy
profiles {
    mylab {
        process {
            cpus = 16
            memory = '64.GB'
            queue = 'long'
        }
        docker.enabled = true
    }
}
```

Use it:
```bash
nextflow run main.nf --species 'E. coli' -profile mylab
```
