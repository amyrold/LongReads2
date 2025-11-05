# LongReads2 Docker Containers

This directory contains Dockerfiles for the three containers used in the LongReads2 pipeline.

## Containers

### 1. ncbi-tools
**Purpose**: Genome metadata download, filtering, and genome downloading

**Tools included**:
- NCBI datasets CLI
- NCBI dataformat CLI
- NCBI BLAST+
- NCBI Entrez Direct
- Python 3 with pandas, numpy, BioPython

**Used by processes**:
- DOWNLOAD_METADATA
- FILTER_LONGREADS
- DOWNLOAD_GENOMES

### 2. sequence-analysis
**Purpose**: 16S sequence extraction, alignment, and distance calculation

**Tools included**:
- Barrnap (rRNA detection)
- MAFFT (multiple sequence alignment)
- HMMER
- Python 3 with BioPython, pandas, numpy, scipy

**Used by processes**:
- EXTRACT_16S
- ALIGN_SEQUENCES
- CALCULATE_DISTANCES

### 3. visualization
**Purpose**: Statistical analysis and visualization

**Tools included**:
- R 4.3
- tidyverse, ggplot2, dplyr
- pheatmap, vegan, ape
- ggtree (Bioconductor)

**Used by processes**:
- VISUALIZE_RESULTS

## Platform Support

All containers are **cross-platform compatible** and automatically detect your system architecture:

- ✅ **ARM64**: Apple Silicon (M1/M2/M3/M4), AWS Graviton
- ✅ **AMD64**: Intel/AMD processors (x86_64)

Docker automatically downloads the correct binaries for your platform during build. No configuration needed!

## Building Containers

### Build all containers locally
```bash
cd containers
./build_containers.sh
```

### Build and push to Docker Hub
```bash
# Login to Docker Hub
docker login

# Build and push
./build_containers.sh --push --registry your_dockerhub_username
```

### Build individual containers
```bash
docker build -t longreads2/ncbi-tools:latest ncbi-tools/
docker build -t longreads2/sequence-analysis:latest sequence-analysis/
docker build -t longreads2/visualization:latest visualization/
```

## Testing Containers

Test that containers work correctly:

```bash
# Test NCBI tools
docker run --rm longreads2/ncbi-tools:latest datasets --version
docker run --rm longreads2/ncbi-tools:latest blastn -version

# Test sequence analysis
docker run --rm longreads2/sequence-analysis:latest barrnap --version
docker run --rm longreads2/sequence-analysis:latest mafft --version

# Test visualization
docker run --rm longreads2/visualization:latest R --version
```

## Converting to Singularity

For HPC environments, convert Docker images to Singularity:

```bash
singularity build ncbi-tools.sif docker://longreads2/ncbi-tools:latest
singularity build sequence-analysis.sif docker://longreads2/sequence-analysis:latest
singularity build visualization.sif docker://longreads2/visualization:latest
```

## Container Sizes

Approximate sizes:
- ncbi-tools: ~800 MB
- sequence-analysis: ~1.2 GB
- visualization: ~1.5 GB

## Updating Containers

When updating tool versions or adding new dependencies:

1. Edit the appropriate Dockerfile
2. Rebuild the container: `docker build -t longreads2/CONTAINER:latest CONTAINER/`
3. Test the container
4. Push to registry: `docker push longreads2/CONTAINER:latest`
5. Update the version tag in `nextflow.config` if needed

### Cross-Platform Architecture Detection

The `ncbi-tools` container uses Docker's `TARGETARCH` variable to automatically download the correct binaries:

```dockerfile
ARG TARGETARCH
RUN if [ "$TARGETARCH" = "arm64" ]; then \
        ARCH="arm64"; \
    else \
        ARCH="amd64"; \
    fi && \
    curl -o /usr/local/bin/datasets "https://...linux-${ARCH}/datasets"
```

This ensures the same Dockerfile works on both Apple Silicon and Intel/AMD machines without modification.
