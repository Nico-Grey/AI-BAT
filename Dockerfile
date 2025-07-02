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
    'RColorBrewer', 'gridExtra', 'tibble', 'shiny', 'DT', 'png', 'grid'), \
    repos='https://cloud.r-project.org/')"
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

# Install specific versions of CRAN packages
RUN R -e "install.packages(c( \
    'magrittr'='2.0.3', 'dplyr'='1.1.4', 'data.table'='1.15.4', 'ggplot2'='3.5.1', 'tidyr'='1.3.1', 'purrr'='1.0.2', \
    'RColorBrewer'='1.1-3', 'gridExtra'='2.3', 'tibble'='3.2.1', 'shiny'='1.10.0', 'DT'='0.33', 'png'='0.1-8', 'grid'='4.3.2' ), \
    repos='https://cloud.r-project.org/')"

# Install specific versions of Bioconductor packages
RUN R -e "BiocManager::install(c('Biostrings', 'ComplexHeatmap'), ask = FALSE, update = FALSE, dependencies = TRUE)"
RUN R -e "BiocManager::install(c('MSnbase', 'sva', 'mixOmics'), ask = FALSE, update = FALSE, dependencies = TRUE)"
RUN R -e "BiocManager::install(c('AnnotationDbi', 'IRanges', 'S4Vectors', 'GenomeInfoDb'), ask = FALSE, update = FALSE, dependencies = TRUE)"
RUN R -e "BiocManager::install(c('SummarizedExperiment', 'GenomicRanges', 'DESeq2', 'edgeR', 'limma'), ask = FALSE, update = FALSE, dependencies = TRUE)"
RUN R -e "BiocManager::install(c('org.Hs.eg.db', 'org.Mm.eg.db', 'clusterProfiler', 'HarmonizR', 'biomaRt', 'imputeLCMD'), ask = FALSE, update = FALSE, dependencies = TRUE)"

# Copy all necessary files into the container
COPY BATpipeline.R /app/BATpipeline.R
COPY dataset.rds /app/dataset.rds
COPY general_functions.R /app/general_functions.R
COPY app.R /app/app.R

# Expose Shiny's default port
EXPOSE 3838  

# Set the command to run the Shiny app when the container starts
CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]
