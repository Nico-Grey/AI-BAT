# filepath: /ai-bats-multiarch/ai-bats-multiarch/docker/Dockerfile
# =============================================================================
# AI-BATS Multi-Architecture Docker Container Configuration
# =============================================================================
#
# This Dockerfile creates a containerized environment for the AI-BATS
# (Artificial Intelligence-driven Brown Adipose Tissue Scoring) application.
# The container includes R, Shiny, Python, and all required dependencies for
# protein expression analysis and machine learning workflows.
#
# This version supports building images for both x86 and ARM architectures.
#
# Author: AI-BATS Team
# Date: July 2025
# Version: 1.0
#
# Usage:
#   docker buildx build --platform linux/amd64,linux/arm64 -t my_r_shiny_app .
#   docker run -p 3838:3838 my_r_shiny_app
# =============================================================================

# Use the Shiny-verse base image for R and Shiny
ARG TARGETPLATFORM
FROM rocker/shiny-verse:4.3.3 AS base

# Install system dependencies required for some R packages
RUN apt-get update && apt-get install -y \
    libnetcdf-dev \
    libhdf5-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libglpk-dev \
    libgmp-dev \
    curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /app

# Install Bioconductor manager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

# Install specific versions of CRAN packages
RUN R -e "install.packages(c( \
    'magrittr'='2.0.3', 'dplyr'='1.1.4', 'data.table'='1.15.4', 'ggplot2'='3.5.1', 'tidyr'='1.3.1', 'purrr'='1.0.2', \
    'RColorBrewer'='1.1-3', 'gridExtra'='2.3', 'tibble'='3.2.1', 'shiny'='1.10.0', 'DT'='0.33', 'png'='0.1-8', 'grid'='4.3.2' ), \
    repos='https://cloud.r-project.org/')"

# Install Bioconductor packages for bioinformatics analysis
RUN R -e "BiocManager::install(c('Biostrings', 'ComplexHeatmap', 'MSnbase', 'sva', 'mixOmics', 'AnnotationDbi', 'IRanges', 'S4Vectors', 'GenomeInfoDb', 'SummarizedExperiment', 'GenomicRanges', 'DESeq2', 'edgeR', 'limma', 'org.Hs.eg.db', 'org.Mm.eg.db', 'clusterProfiler', 'HarmonizR', 'biomaRt', 'imputeLCMD'), ask = FALSE, update = FALSE, dependencies = TRUE)"

# Copy application source code and data files into the container
COPY ./ /app

# Install Python via Miniconda for machine learning components
#RUN wget \
#    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
#    && mkdir /root/.conda \
#    && bash Miniconda3-latest-Linux-x86_64.sh -b \
#    && rm -f Miniconda3-latest-Linux-x86_64.sh 

ARG TARGETARCH
ENV CONDA_DIR=/opt/conda
RUN set -e; \
    case "${TARGETARCH}" in \
      amd64) CONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" ;; \
      arm64) CONDA_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh" ;; \
      *) echo "Unsupported TARGETARCH: ${TARGETARCH}" >&2; exit 1 ;; \
    esac; \
    echo "Downloading ${CONDA_URL} for ${TARGETARCH}"; \
    curl -fsSL "${CONDA_URL}" -o /tmp/conda.sh; \
    bash /tmp/conda.sh -b -p "${CONDA_DIR}"; \
    rm /tmp/conda.sh; \
    ln -s "${CONDA_DIR}/etc/profile.d/conda.sh" /etc/profile.d/conda.sh; \
    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh" >> /etc/bash.bashrc
ENV PATH="${CONDA_DIR}/bin:${PATH}"


# Copy Python environment specification
COPY ./environment.yml /app/environment.yml

RUN conda tos accept --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --channel https://repo.anaconda.com/pkgs/r

# Create conda environment and initialize conda shell integration
RUN conda env create -f /app/environment.yml -y \
    && conda clean -afy \
    && conda init

# Activate the conda environment by default
RUN /opt/conda/bin/conda run -n aibats python python_scripts/download_all_models.py
#SHELL ["/opt/conda/bin/conda","run","-n","aibats","/bin/bash","-c"]
#RUN apt-get install -y git
#RUN pip install "tabpfn-extensions[all] @ git+https://github.com/PriorLabs/tabpfn-extensions.git"
#RUN conda activate aibats \ && python python_scripts/download_all_models.py

# Activate environment by default and update PATH
RUN echo "conda activate aibats" >> ~/.bashrc
ENV PATH=/opt/conda/envs/aibats/bin:$PATH

# Expose Shiny's default port
EXPOSE 3838  

# Set the command to run the Shiny app when the container starts
CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]