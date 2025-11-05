# Container for NCBI tools: datasets, dataformat, BLAST, and Entrez
FROM ubuntu:22.04

LABEL maintainer="Aaron Myrold"
LABEL description="NCBI tools for genome download and BLAST"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    unzip \
    python3 \
    python3-pip \
    ncbi-blast+ \
    ncbi-entrez-direct \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip3 install --no-cache-dir \
    pandas \
    numpy \
    biopython

# Download and install NCBI datasets and dataformat
# Docker automatically sets TARGETARCH (amd64, arm64, etc.)
ARG TARGETARCH
RUN if [ "$TARGETARCH" = "arm64" ]; then \
        ARCH="arm64"; \
    else \
        ARCH="amd64"; \
    fi && \
    curl -o /usr/local/bin/datasets "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-${ARCH}/datasets" && \
    curl -o /usr/local/bin/dataformat "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-${ARCH}/dataformat" && \
    chmod +x /usr/local/bin/datasets /usr/local/bin/dataformat

# Set working directory
WORKDIR /work

# Default command
CMD ["/bin/bash"]
