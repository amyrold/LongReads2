# Container for visualization and statistical analysis with R
FROM rocker/r-ver:4.3

LABEL maintainer="Aaron Myrold"
LABEL description="R environment for statistical analysis and visualization"

# Install system dependencies for R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libwebp-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c( \
    'tidyverse', \
    'ggplot2', \
    'dplyr', \
    'readr', \
    'pheatmap', \
    'vegan', \
    'ape', \
    'BiocManager' \
    ), repos='https://cloud.r-project.org/')"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('ggtree', 'Biostrings'))"

# Set working directory
WORKDIR /work

# Default command
CMD ["R"]
