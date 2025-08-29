# BATpipeline.R - Brown Adipose Tissue Analysis Pipeline

## Overview

BATpipeline.R is the core R script that performs comprehensive analysis of protein expression data for brown and white adipose tissue classification. This pipeline handles data preprocessing, normalization, batch correction, dimensionality reduction, and machine learning-based scoring for adipose tissue browning analysis.

## Features

### 🔬 Data Processing
- **Multi-format Support**: CSV, RDS, RData file compatibility
- **Automatic Normalization**: Log transformation and median centering
- **Batch Correction**: ComBat and HarmonizR integration
- **Missing Value Imputation**: LCMD-based imputation methods

### 🧮 Analysis Methods
- **Dimensionality Reduction**: PCA, t-SNE, UMAP
- **Machine Learning**: Random Forest, LDA, Logistic Regression
- **Feature Selection**: LASSO regularization
- **Cross-validation**: Stratified k-fold validation

### 📊 Visualization
- **Quality Control Plots**: Before/after normalization
- **Dimensionality Reduction Plots**: PCA, t-SNE, UMAP visualizations
- **Score Distributions**: Brown vs. white tissue scoring
- **Feature Importance**: Model-specific feature rankings

## Dependencies

### Required R Packages
```r
# Data manipulation
library(magrittr)
library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(purrr)

# Bioinformatics
library(Biostrings)
library(biomaRt)
library(MSnbase)
library(imputeLCMD)

# Statistical analysis
library(mixOmics)
library(sva)
library(HarmonizR)

# Visualization
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

# Python integration
library(reticulate)
```

## Installation

### Prerequisites
1. Install R (version 4.0 or higher)
2. Install required packages:
```r
# CRAN packages
install.packages(c("magrittr", "dplyr", "data.table", "tidyr", 
                   "tibble", "purrr", "ggplot2", "RColorBrewer", 
                   "gridExtra", "reticulate"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "biomaRt", "MSnbase", 
                       "imputeLCMD", "mixOmics", "sva"))
```

3. Install Python dependencies for integration:
```bash
pip install numpy pandas scikit-learn
```

## Usage

### Command Line Execution
```bash
Rscript BATpipeline.R <input_path> <output_data_path> <plot_dir> <param_path>
```

### Parameters
- `input_path`: Path to input data file (CSV/RDS/RData)
- `output_data_path`: Path for output processed data
- `plot_dir`: Directory for saving generated plots
- `param_path`: Path to parameters file (RDS format)

### Example
```bash
Rscript BATpipeline.R data/protein_data.csv results/processed_data.rds plots/ params.rds
```

## Data Format Requirements

### Input Data Structure

#### Expected Format
The input should be a named list containing datasets with the following structure:
```r
data_structure <- list(
  "dataset1_batch1" = list(
    intensity = matrix,  # rows = genes/proteins, cols = samples
    metadata = data.frame  # sample information
  ),
  "dataset2_batch2" = list(
    intensity = matrix,
    metadata = data.frame
  )
)
```

#### Intensity Matrix
- **Rows**: Gene/protein identifiers
- **Columns**: Sample identifiers
- **Values**: Expression intensities or protein abundances

#### Metadata Requirements
- Sample identifiers matching matrix columns
- Batch information for correction
- Tissue type labels (Brown/White/Intermediate)

### File Format Support

#### CSV Files
- First column: Gene/protein identifiers
- Subsequent columns: Sample expression values
- Header row required

#### RDS Files
- Serialized R objects
- Must contain properly structured list as described above

#### RData Files
- R workspace files
- First object in workspace used as input data

## Processing Pipeline

### 1. Data Loading and Validation
```r
# File format detection
ext <- tolower(tools::file_ext(input_path))

# Format-specific loading
if (ext == "rds") {
  data_stes <- readRDS(input_path)
} else if (ext == "csv") {
  data_stes <- data.table::fread(input_path)
}
```

### 2. Normalization Assessment
```r
# Check maximum values to determine normalization need
max_vals <- sapply(data_stes, function(x) max(x$intensity, na.rm = TRUE))
datasets_to_normalize <- names(data_stes)[which(max_vals > 10000)]
```

### 3. Log Transformation and Median Centering
```r
# For datasets requiring normalization
log_data <- log1p(data_stes[[i]]$intensity)
med_col <- apply(log_data, 2, function(x) median(x, na.rm = TRUE))
normalized_data <- log_data - matrix(rep(med_col, nrow(log_data)), 
                                   nrow = nrow(log_data), byrow = TRUE)
```

### 4. Gene Union and Matrix Preparation
```r
# Combine all genes across datasets
all_genes <- Reduce(union, lapply(data_stes, function(x) rownames(x$intensity)))

# Create unified matrices
all_matrices <- lapply(data_stes, function(x) {
  x$normalized[all_genes, ]
})
```

### 5. Batch Correction
Multiple batch correction methods are supported:
- **ComBat**: Empirical Bayes batch correction
- **HarmonizR**: Advanced harmonization methods

### 6. Quality Control Visualization
- Before/after normalization comparisons
- Batch effect assessment plots
- Distribution analysis

## Output Files

### Processed Data
- **Format**: RDS file containing processed expression matrix
- **Content**: Normalized, batch-corrected expression data
- **Structure**: Matrix with genes as rows, samples as columns

### Visualization Outputs
Generated plots include:
- `normalization_comparison.png`: Before/after normalization
- `batch_correction_pca.png`: PCA before/after batch correction
- `quality_control_boxplots.png`: Expression distribution summaries
- `correlation_heatmap.png`: Sample correlation matrix

### Log Files
- Processing steps and parameters
- Quality control metrics
- Error messages and warnings

## Configuration

### Normalization Parameters
```r
# Threshold for normalization decision
normalization_threshold <- 10000

# Log transformation method
log_method <- "log1p"  # log(1 + x)

# Centering method
centering_method <- "median"  # column-wise median centering
```

### Batch Correction Settings
```r
# ComBat parameters
combat_par_prior <- TRUE
combat_mean_only <- FALSE

# HarmonizR parameters
harmonizr_algorithm <- "ComBat"
harmonizr_adjust_bio <- TRUE
```

## Integration with Python

The pipeline supports Python integration for advanced machine learning:

```r
# Python environment setup
library(reticulate)
use_python("/path/to/python")

# Call Python analysis scripts
system("python script_modularized/browning_pipeline.py")
```

## Troubleshooting

### Common Issues

#### Memory Errors
```r
# Increase memory limit
options(java.parameters = "-Xmx8g")
memory.limit(size = 8000)
```

#### File Format Problems
- **Issue**: "Unsupported extension" error
- **Solution**: Ensure file has .csv, .rds, or .rdata extension

#### Missing Data Handling
```r
# Check for missing values
missing_summary <- sapply(data_list, function(x) {
  sum(is.na(x$intensity))
})
```

#### Batch Effect Issues
- **Issue**: Strong batch effects persist
- **Solution**: Try different batch correction methods or parameters

### Performance Optimization

#### Large Datasets
```r
# Use data.table for faster operations
library(data.table)
setDTthreads(threads = 4)

# Parallel processing
library(parallel)
cl <- makeCluster(detectCores() - 1)
```

#### Memory Management
```r
# Clear intermediate objects
rm(large_object)
gc()

# Use sparse matrices for large datasets
library(Matrix)
sparse_matrix <- Matrix(dense_matrix, sparse = TRUE)
```

## Advanced Usage

### Custom Normalization
```r
# Implement custom normalization function
custom_normalize <- function(intensity_matrix) {
  # Your normalization logic here
  return(normalized_matrix)
}
```

### Adding New Batch Correction Methods
```r
# Template for new batch correction
new_batch_correction <- function(data_matrix, batch_labels) {
  # Implementation
  return(corrected_matrix)
}
```

## Quality Control Metrics

### Normalization Assessment
- Coefficient of variation reduction
- Distribution symmetry improvement
- Outlier detection and handling

### Batch Correction Validation
- Principal component analysis
- Silhouette analysis
- Batch effect quantification

## Support and Maintenance

### Debugging
Enable verbose output:
```r
options(verbose = TRUE)
debug(function_name)
```

### Logging
```r
# Custom logging function
log_message <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg))
}
```

### Version Compatibility
- R version: 4.0+
- Bioconductor version: 3.12+
- Python version: 3.7+ (for integration)

## Performance Benchmarks

### Typical Processing Times
- Small dataset (< 1000 genes, < 100 samples): 1-2 minutes
- Medium dataset (1000-5000 genes, 100-500 samples): 5-10 minutes  
- Large dataset (> 5000 genes, > 500 samples): 15-30 minutes

### Memory Requirements
- Small dataset: < 1GB RAM
- Medium dataset: 2-4GB RAM
- Large dataset: 8GB+ RAM recommended
