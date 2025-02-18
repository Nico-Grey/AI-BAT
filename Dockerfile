# Use the Shiny-verse base image for R and Shiny
FROM rocker/shiny-verse:4.3.3

# Install system dependencies required for some R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /app

# Install Bioconductor manager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

# Install all CRAN packages
RUN R -e "install.packages(c( \
    'magrittr', 'dplyr', 'data.table', 'ggplot2', 'tidyr', 'HarmonizR', 'biomaRt', \
    'purrr', 'RColorBrewer', 'gridExtra', 'tibble', 'imputeLCMD'), \
    repos='https://cloud.r-project.org/')"

# Install all Bioconductor packages separately
RUN R -e "BiocManager::install(c( \
    'Biostrings', 'ComplexHeatmap', 'MSnbase', 'sva', 'mixOmics'), \
    ask = FALSE, update = TRUE)"

# Copy your R script into the container
COPY BATpipeline.R /app/BATpipeline.R

# Copy datasets into the container
COPY dataset.rds /app/dataset.rds


# Copy additional R scripts into the container
COPY general_functions.R /app/general_functions.R


# Command to run your R script when the container starts
CMD ["Rscript", "/app/BATpipeline.R"]
