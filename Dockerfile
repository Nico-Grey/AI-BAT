# =============================================================================
# AI-BATS Docker Container Configuration
# =============================================================================
# 
# This Dockerfile creates a containerized environment for the AI-BATS
# (Artificial Intelligence-driven Brown Adipose Tissue Scoring) application.
# The container includes R, Shiny, Python, and all required dependencies for
# protein expression analysis and machine learning workflows.
#
# Author: AI-BATS Team
# Date: July 2025
# Version: 1.0
#
# Container Components:
# 1. Base R/Shiny environment (rocker/shiny-verse)
# 2. System dependencies for R packages
# 3. R packages (CRAN and Bioconductor)
# 4. Python environment via Miniconda
# 5. Application files and configuration
# 6. Network and runtime configuration
#
# Usage:
#   docker build -t my_r_shiny_app .
#   docker run -p 3838:3838 my_r_shiny_app
# =============================================================================

# =============================================================================
# BASE IMAGE SELECTION
# =============================================================================
# Use rocker/shiny-verse as base image - provides R 4.3.3 with Shiny and tidyverse
# This image includes pre-installed R packages commonly used in data science

# Use the Shiny-verse base image for R and Shiny
FROM rocker/shiny-verse:4.3.3

# =============================================================================
# SYSTEM DEPENDENCIES INSTALLATION
# =============================================================================
# Install system-level libraries required for R packages compilation
# These libraries support various R packages for data processing and analysis
# - libnetcdf-dev: NetCDF library for scientific data formats
# - libhdf5-dev: HDF5 library for hierarchical data storage
# - libxml2-dev: XML parsing library
# - libcurl4-openssl-dev: HTTP/HTTPS client library with SSL support
# - libssl-dev: OpenSSL development libraries
# - libglpk-dev: GNU Linear Programming Kit
# - libgmp-dev: GNU Multiple Precision Arithmetic Library

# Install system dependencies required for some R packages
RUN apt-get update && apt-get install -y \
    libnetcdf-dev \
    libhdf5-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libglpk-dev \
    libgmp-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# =============================================================================
# WORKING DIRECTORY SETUP
# =============================================================================
# Set the working directory inside the container where application files will reside

# Set the working directory
WORKDIR /app

# =============================================================================
# R PACKAGE MANAGER INSTALLATION
# =============================================================================
# Install BiocManager for managing Bioconductor packages

# Install Bioconductor manager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

# =============================================================================
# CRAN PACKAGES INSTALLATION
# =============================================================================
# Install specific versions of R packages from CRAN
# Version pinning ensures reproducibility and compatibility
# Packages include:
# - Data manipulation: magrittr, dplyr, data.table, tidyr, purrr, tibble
# - Visualization: ggplot2, RColorBrewer, gridExtra, png, grid
# - Web application: shiny, DT

# Install specific versions of CRAN packages
RUN R -e "install.packages(c( \
    'magrittr'='2.0.3', 'dplyr'='1.1.4', 'data.table'='1.15.4', 'ggplot2'='3.5.1', 'tidyr'='1.3.1', 'purrr'='1.0.2', \
    'RColorBrewer'='1.1-3', 'gridExtra'='2.3', 'tibble'='3.2.1', 'shiny'='1.10.0', 'DT'='0.33', 'png'='0.1-8', 'grid'='4.3.2' ), \
    repos='https://cloud.r-project.org/')"

# =============================================================================
# BIOCONDUCTOR PACKAGES INSTALLATION
# =============================================================================
# Install Bioconductor packages for bioinformatics analysis
# These packages provide specialized functionality for:
# - Biological sequence analysis (Biostrings)
# - Mass spectrometry data (MSnbase)
# - Batch effect removal (sva)
# - Multivariate omics analysis (mixOmics)
# - Gene annotation and pathway analysis
# - Data harmonization and imputation

# Install core Bioconductor packages for data structures and visualization
RUN R -e "BiocManager::install(c('Biostrings', 'ComplexHeatmap'), ask = FALSE, update = FALSE, dependencies = TRUE)"

# Install packages for mass spectrometry and batch correction
RUN R -e "BiocManager::install(c('MSnbase', 'sva', 'mixOmics'), ask = FALSE, update = FALSE, dependencies = TRUE)"

# Install core Bioconductor infrastructure packages
RUN R -e "BiocManager::install(c('AnnotationDbi', 'IRanges', 'S4Vectors', 'GenomeInfoDb'), ask = FALSE, update = FALSE, dependencies = TRUE)"

# Install genomics and differential expression packages
RUN R -e "BiocManager::install(c('SummarizedExperiment', 'GenomicRanges', 'DESeq2', 'edgeR', 'limma'), ask = FALSE, update = FALSE, dependencies = TRUE)"

# Install annotation databases and specialized analysis packages
RUN R -e "BiocManager::install(c('org.Hs.eg.db', 'org.Mm.eg.db', 'clusterProfiler', 'HarmonizR', 'biomaRt', 'imputeLCMD'), ask = FALSE, update = FALSE, dependencies = TRUE)"

# =============================================================================
# APPLICATION FILES DEPLOYMENT
# =============================================================================
# Copy application source code and data files into the container
# Files include:
# - BATpipeline.R: Main analysis pipeline script
# - dataset.rds: Example/test dataset for analysis
# - general_functions.R: Utility functions for analysis
# - app.R: Shiny web application interface

# Copy all necessary files into the container
COPY BATpipeline.R /app/BATpipeline.R
COPY dataset.rds /app/dataset.rds
COPY general_functions.R /app/general_functions.R
COPY app.R /app/app.R

# =============================================================================
# PYTHON ENVIRONMENT SETUP
# =============================================================================
# Install Python via Miniconda for machine learning components
# This provides a complete Python environment with conda package management

# Install Python and required packages
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 

# =============================================================================
# PYTHON PATH CONFIGURATION
# =============================================================================
# Add conda to system PATH for Python package management

# Add conda to PATH
ENV PATH="/root/miniconda3/bin:${PATH}"

# =============================================================================
# CONDA ENVIRONMENT SETUP
# =============================================================================
# Create Python environment from environment.yml file
# This installs all Python dependencies for machine learning analysis

# Copy Python environment specification
COPY environment.yml /app/environment.yml

# Create conda environment and initialize conda shell integration
RUN conda env create -f environment.yml -y \
    && conda clean -afy \
    && conda init

# Activate the conda environment by default
SHELL ["conda", "run", "-n", "aibats", "/bin/bash", "-c"]
RUN pip install "tabpfn-extensions[all] @ git+https://github.com/PriorLabs/tabpfn-extensions.git"
# =============================================================================
# PYTHON ENVIRONMENT ACTIVATION
# =============================================================================
# Configure shell to automatically activate the AI-BATS Python environment

# Activate environment by default and update PATH
RUN echo "source activate aibats" > ~/.bashrc
ENV PATH=/opt/conda/envs/aibats/bin:$PATH
# =============================================================================
# NETWORK CONFIGURATION
# =============================================================================
# Expose port 3838 for Shiny application access
# This is the default port used by Shiny Server

# Expose Shiny's default port
EXPOSE 3838  

# =============================================================================
# CONTAINER STARTUP COMMAND
# =============================================================================
# Define the default command to run when container starts
# Launches the Shiny application on all network interfaces (0.0.0.0)
# This allows external access to the containerized application

# Set the command to run the Shiny app when the container starts
CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]
