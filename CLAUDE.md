# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

LongReads2 is a Nextflow DSL2 pipeline for analyzing 16S rRNA variation within and between bacterial strains using long-read sequencing data. This is a complete rewrite from Python/R scripts to a production-ready workflow with containerized components.

**Key improvement**: Uses Barrnap (HMM-based) + MAFFT alignment + p-distance instead of the original BLAST + raw edit distance approach.

## Running the Pipeline

### Basic Test Run (5 genomes, ~30 min)
```bash
nextflow run main.nf \
  --species 'Escherichia coli' \
  --max_genomes 5 \
  --outdir test_results \
  -profile docker
```

### Resume Failed Runs
```bash
nextflow run main.nf -resume --species 'Your Species' --outdir results
```

### Using Different Profiles
- `docker` - Docker containers (default)
- `singularity` - Singularity for HPC
- `slurm` - SLURM scheduler + Singularity
- `test` - Pre-configured test with 5 E. coli genomes
- `local` - No containers (assumes tools installed locally)

### Key Parameters
- `--species` - Bacterial species (quoted if spaces)
- `--max_genomes` - Limit number of genomes (default: all)
- `--use_barrnap` - Use Barrnap vs BLAST (default: true)
- `--min_16s_length` / `--max_16s_length` - Length filters (default: 1400-1600 bp)
- `--cpus` / `--memory` - Resource allocation per process

### View Help
```bash
nextflow run main.nf --help
```

## Building Containers

### Build All Containers
```bash
cd containers
./build_containers.sh
```

### Build Individual Container
```bash
docker build -t longreads2/ncbi-tools:latest containers/ncbi-tools/
docker build -t longreads2/sequence-analysis:latest containers/sequence-analysis/
docker build -t longreads2/visualization:latest containers/visualization/
```

### Test Containers
```bash
docker run --rm longreads2/ncbi-tools:latest datasets --version
docker run --rm longreads2/sequence-analysis:latest barrnap --version
docker run --rm longreads2/visualization:latest R --version
```

## Architecture

### Pipeline Flow (7 Sequential Processes)
```
DOWNLOAD_METADATA → FILTER_LONGREADS → DOWNLOAD_GENOMES → EXTRACT_16S →
ALIGN_SEQUENCES → CALCULATE_DISTANCES → VISUALIZE_RESULTS
```

### Module Organization
Each process lives in `modules/local/<process_name>/main.nf` with:
- Input/output channel definitions
- Container specification (via nextflow.config)
- Script block (bash, Python, or R)
- Optional stub block for dry-runs

### Container-to-Process Mapping

**ncbi-tools container** (used by):
- DOWNLOAD_METADATA - NCBI datasets API
- FILTER_LONGREADS - Python filtering script
- DOWNLOAD_GENOMES - Batch download with datasets CLI

**sequence-analysis container** (used by):
- EXTRACT_16S - Barrnap or BLAST for 16S detection
- ALIGN_SEQUENCES - MAFFT multiple sequence alignment
- CALCULATE_DISTANCES - BioPython distance calculations

**visualization container** (used by):
- VISUALIZE_RESULTS - R with ggplot2, pheatmap, vegan

### Key Data Transformations

1. **EXTRACT_16S**: Assigns unique IDs as `{GENOME_ACCESSION}_16S_{COPY_NUMBER}`
2. **ALIGN_SEQUENCES**: MAFFT algorithm auto-selected based on sequence count:
   - <200 seqs: L-INS-i (accurate)
   - 200-500 seqs: E-INS-i (balanced)
   - >500 seqs: FFT-NS-2 (fast)
3. **CALCULATE_DISTANCES**: Separates within-genome vs between-genome comparisons based on genome accession prefix

### Output Structure
```
{outdir}/
├── 1_raw_data/          # Downloaded genomes & metadata
├── 2_filtered_data/     # 16S sequences & alignment
├── 3_output/            # Distance matrices & plots
└── pipeline_info/       # Execution reports (timeline, DAG, trace)
```

## Adding New Modules

1. Create `modules/local/new_process/main.nf`
2. Define process with inputs/outputs using `.emit` for named channels
3. Add container specification in `nextflow.config`:
   ```groovy
   withName: 'NEW_PROCESS' {
       container = 'longreads2/container-name:latest'
       cpus = 4
       memory = 8.GB
   }
   ```
4. Include in `main.nf`:
   ```groovy
   include { NEW_PROCESS } from './modules/local/new_process/main.nf'
   ```
5. Wire into workflow block with channel connections

## Configuration Details

### Resource Management
`nextflow.config` uses `check_max()` helper function to enforce limits:
- `params.max_cpus` = 16
- `params.max_memory` = 128.GB
- `params.max_time` = 240.h

Process-specific resources scale with `task.attempt` for retry logic.

### Error Handling
Exit codes 143, 137, 104, 134, 139 trigger retry (up to 2 attempts).
All other errors finish immediately.

### Execution Reports
Automatically generated in `{outdir}/pipeline_info/`:
- `execution_timeline.html` - When each process ran
- `execution_report.html` - Resource usage per process
- `pipeline_dag.svg` - Workflow diagram
- `execution_trace.txt` - Raw trace data

## Methodological Details

### Why Barrnap over BLAST?
- Species-agnostic (uses bacterial HMM models)
- More sensitive for divergent 16S sequences
- Provides strand information directly
- ~10x faster than BLAST

### Why MAFFT Alignment?
Original used raw edit distances which don't account for insertions/deletions properly.
MAFFT alignment enables biologically meaningful p-distance calculations.

### Distance Metrics
- **p-distance**: Proportion of differing sites in aligned sequences
- Used instead of raw Levenshtein distance for proper gap handling
- Standard metric in phylogenetics for closely related sequences

## Debugging

### View Process Work Directory
```bash
# Find work directory hash from error or .nextflow.log
cd work/a1/b2c3d4...
cat .command.sh      # Script that ran
cat .command.log     # stdout/stderr
cat .command.err     # stderr only
cat .exitcode        # Exit code
```

### Common Issues

**No 16S found**:
- Try BLAST mode: `--use_barrnap false`
- Relax filters: `--min_16s_length 1200`

**Memory errors**:
- Increase: `--memory '32.GB'`
- Reduce dataset: `--max_genomes 10`

**Container pull failures**:
- Pre-pull: `docker pull longreads2/sequence-analysis:latest`
- Check Docker daemon is running

### Nextflow Cache
Work directory (`work/`) caches all process results. With `-resume`, Nextflow reuses cached outputs if inputs unchanged. Cache is NOT persistent across different parameters (species, genome count, etc.).

## Documentation Structure

- `README.md` - Main user documentation
- `docs/QUICKSTART.md` - 10-minute getting started guide
- `docs/WORKFLOW.md` - Detailed process descriptions & data flow
- `docs/IMPLEMENTATION_NOTES.md` - Technical decisions & learning outcomes
- `containers/README.md` - Container build instructions
- `PROJECT_SUMMARY.md` - High-level project summary

## Important Notes

- **DSL2 syntax**: All processes use `nextflow.enable.dsl=2`
- **Channel emission**: Outputs use `.emit` for named channels in workflow
- **Container strategy**: Three specialized containers instead of one monolithic image
- **Cross-platform support**: Containers auto-detect architecture (ARM64/AMD64) using Docker's `TARGETARCH`
- **publishDir mode**: Uses 'copy' to preserve outputs even after work directory cleanup
- **Python in processes**: Uses heredoc syntax `python3 << 'EOF'` for embedded scripts
- **R in processes**: Uses `#!/usr/bin/env Rscript` shebang for embedded R

## Platform Compatibility

All Docker containers automatically detect the build architecture:

- **Apple Silicon (M1/M2/M3/M4)**: Builds with ARM64 binaries
- **Intel/AMD**: Builds with AMD64 binaries
- **AWS Graviton**: Uses ARM64 binaries

The `ncbi-tools` Dockerfile uses `ARG TARGETARCH` to conditionally download platform-specific NCBI datasets binaries. No manual configuration needed - the same Dockerfile works everywhere.
