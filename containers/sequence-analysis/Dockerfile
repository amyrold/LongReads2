# Container for sequence analysis: Barrnap, MAFFT, BioPython
FROM ubuntu:22.04

LABEL maintainer="Aaron Myrold"
LABEL description="Sequence analysis tools for 16S extraction, alignment, and distance calculation"

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    git \
    build-essential \
    python3 \
    python3-pip \
    mafft \
    perl \
    bedtools \
    && rm -rf /var/lib/apt/lists/*

# Install HMMER from source (not available in ARM64 Ubuntu repos)
RUN wget http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz && \
    tar -xzf hmmer-3.4.tar.gz && \
    cd hmmer-3.4 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf hmmer-3.4 hmmer-3.4.tar.gz

# Install Barrnap for rRNA detection
RUN wget https://github.com/tseemann/barrnap/archive/refs/tags/0.9.tar.gz && \
    tar -xzf 0.9.tar.gz && \
    cd barrnap-0.9 && \
    cp bin/barrnap /usr/local/bin/ && \
    mkdir -p /usr/local/db && \
    cp -r db/* /usr/local/db/ && \
    cd .. && \
    rm -rf barrnap-0.9 0.9.tar.gz

# Set Barrnap database path
ENV BARRNAP_DB=/usr/local/db

# Install Python packages for sequence analysis
RUN pip3 install --no-cache-dir \
    pandas \
    numpy \
    biopython \
    scipy \
    editdistance

# Set working directory
WORKDIR /work

# Default command
CMD ["/bin/bash"]
