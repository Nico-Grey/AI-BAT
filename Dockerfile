# Use the Shiny-verse base image for R and Shiny
FROM rocker/shiny-verse:4.3.3

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

# Set the working directory
WORKDIR /app

# Install Bioconductor manager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

# Install CRAN packages in a separate step to use Docker caching
RUN R -e "install.packages(c( \
    'magrittr', 'dplyr', 'data.table', 'ggplot2', 'tidyr', 'purrr', \
    'RColorBrewer', 'gridExtra', 'tibble'), \
    repos='https://cloud.r-project.org/')"

# Install Bioconductor packages in **small groups** for better caching
RUN R -e "BiocManager::install(c('Biostrings', 'ComplexHeatmap'), ask = FALSE, update = TRUE, dependencies = TRUE)"
RUN R -e "BiocManager::install(c('MSnbase', 'sva', 'mixOmics'), ask = FALSE, update = TRUE, dependencies = TRUE)"
RUN R -e "BiocManager::install(c('AnnotationDbi', 'IRanges', 'S4Vectors', 'GenomeInfoDb'), ask = FALSE, update = TRUE, dependencies = TRUE)"
RUN R -e "BiocManager::install(c('SummarizedExperiment', 'GenomicRanges', 'DESeq2', 'edgeR', 'limma'), ask = FALSE, update = TRUE, dependencies = TRUE)"
RUN R -e "BiocManager::install(c('org.Hs.eg.db', 'org.Mm.eg.db', 'clusterProfiler', 'HarmonizR', 'biomaRt', 'imputeLCMD'), ask = FALSE, update = TRUE, dependencies = TRUE)"

# Copy project files into the container
COPY BATpipeline.R /app/BATpipeline.R
COPY dataset.rds /app/dataset.rds
COPY general_functions.R /app/general_functions.R

# Command to run your R script when the container starts
CMD ["Rscript", "/app/BATpipeline.R"]
