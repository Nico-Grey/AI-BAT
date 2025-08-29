library(magrittr)
library(dplyr)
library(data.table)
library(Biostrings)
library(ggplot2)
library(tidyr)
library(mixOmics)
library(sva)
library(HarmonizR)
library(biomaRt)
library(sva)
library(purrr)
library(RColorBrewer)
library(gridExtra)
library(tibble)
library(purrr)
library(MSnbase)
library(imputeLCMD)
library(reticulate)
print("Library uploaded")

# Define working directory
working_dir <- getwd()


# Source external script directly (no subfolder)
source(file.path(working_dir, "general_functions.R"))



# ---- parse command‐line args ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript BATpipeline.R <input_path> <output_data_path> <plot_dir> <param_path>")
}

input_path  <- args[1]
output_path <- args[2]
plot_dir    <- args[3]


# Folder for plot saving
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# ---- load the dataset from Shiny mount ----
ext <- tolower(tools::file_ext(input_path))
if (ext == "rds") {
  data_stes <- readRDS(input_path)
} else if (ext == "csv") {
  data_stes <- data.table::fread(input_path, stringsAsFactors = FALSE)
} else if (ext %in% c("rdata", "rda")) {
  temp_env <- new.env()
  load(input_path, envir = temp_env)
  # assume the first object is the dataset
  obj_name  <- ls(temp_env)[1]
  data_stes <- temp_env[[obj_name]]
} else {
  stop(sprintf("Unsupported extension: %s", ext))
}

# Check normalization
max_vals <- c()
for(i in 1:length(data_stes)){
  max_vals <- c(max_vals, max(data_stes[[i]]$intensity, na.rm = TRUE))
}

datasets_to_normalize <- names(data_stes)[which(max_vals > 10000)]

# Normalize each data set individually
for(i in names(data_stes)){
  if(i %in% datasets_to_normalize){
    # Calculate the median per column
    log_data <- log1p(data_stes[[i]]$intensity)
    med_col <- apply(log_data, 2, function(x) median(x, na.rm = T))
    
    data_stes[[i]]$normalized <- (log_data - matrix(rep(med_col, nrow(log_data)), nrow = nrow(log_data), byrow = T))
    
    rownames(data_stes[[i]]$normalized) <- rownames(data_stes[[i]]$intensity)
    
  } else {
    # If no normalization is needed, just copy the original data
    data_stes[[i]]$normalized <- data_stes[[i]]$intensity
  }
}
print("Normalization Complete")


##Batch Correction and Data Plotting
all_genes <- lapply(data_stes, function(x){
  return(rownames(x$intensity))
})

all_genes <- Reduce(union, all_genes)

n <- names(data_stes)
mat_names <- unlist(lapply(n, function(x) strsplit(x, "_")[[1]][1]))
all_matrices <- lapply(data_stes, function(x){
  m <- x$normalized[all_genes,]
  rownames(m) <- all_genes
  return(m)
})

for(i in 1:length(mat_names)){
  colnames(all_matrices[[i]]) <- paste0(mat_names[i],"-",colnames(all_matrices[[i]]))
}

all_matrices_mat <- Reduce(cbind, all_matrices)
common_prot <- all_matrices_mat[complete.cases(all_matrices_mat), ]



df_descroption <- data.frame(ID = colnames(all_matrices_mat), sample = 1:ncol(all_matrices_mat))
batch_vector <- c()
for (i in 1:length(all_matrices)){
  batch_vector <- c(batch_vector, rep(i, ncol(all_matrices[[i]])))
}
df_descroption$batch <- batch_vector
write.csv(df_descroption, file = file.path(getwd(), "all_normalized_description.csv"), row.names = FALSE)
batch_data <- read.csv(file.path(getwd(), "all_normalized_description.csv"), sep = ",", header = TRUE)

batch_names <- c(
  ##Add function for batch names
)
mat_names <- c(
 ##Add function for matrix names
)


# Loop through each dataset in data_stes
for (dataset_name in names(data_stes)) {
  # Check if 'meta' exists and contains the column 'protein_type'
  if ("meta" %in% names(data_stes[[dataset_name]]) &&
      "protein_type" %in% colnames(data_stes[[dataset_name]]$meta)) {
    # Rename 'protein_type' to 'tissue'
    colnames(data_stes[[dataset_name]]$meta)[
      colnames(data_stes[[dataset_name]]$meta) == "protein_type"
    ] <- "tissue"
  }
}


all_tissue_types <- c()
for (l in names(data_stes)) {
  tissue_type <- data_stes[[l]][["meta"]][["tissue"]]
  all_tissue_types <- c(all_tissue_types, tissue_type)
}
all_tissue_types <- factor(all_tissue_types, levels = unique(all_tissue_types))


# Perform PCA - Before Batch Correction
pca_before <- prcomp(t(common_prot), scale.=TRUE)
pca_data_before <- data.frame(
  PC1 = pca_before$x[, 1], 
  PC2 = pca_before$x[, 2], 
  Batch = factor(batch_data$batch, levels = names(batch_names), labels = mat_names)
)
palette <- brewer.pal(n = length(unique(batch_data$batch)), name = "Set3")
# 1) Build the plot and store it in a variable
p_before <- ggplot(pca_data_before, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = palette) +
  ggtitle("PCA Plot Before Batch Correction") +
  theme_minimal() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  ) +
  labs(color = "Batch")

# 2) Save it into your mounted folder
ggsave(
  filename = file.path(plots_dir, "pca_before_batch_correction.png"),
  plot     = p_before,
  width    = 8,      # inches → 800 px at dpi=100
  height   = 6,      # inches → 600 px at dpi=100
  units    = "in",
  dpi      = 100
)



combat_exprs <- ComBat(dat=common_prot, batch=batch_data$batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)


# Perform PCA - After Batch Correction
pca_after <- prcomp(t(combat_exprs), scale.=TRUE)
pca_data_after <- data.frame(
  PC1 = pca_after$x[, 1], 
  PC2 = pca_after$x[, 2], 
  Batch = factor(batch_data$batch, levels = names(batch_names), labels = batch_names)
)
pca_after <- ggplot(pca_data_after, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.8) +  
  scale_color_manual(values = palette) +
  ggtitle("PCA Plot After Batch Correction") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  labs(color = "Batch")

ggsave(
  filename = file.path(plots_dir, "pca_after_batch_correction.png"),
  plot     = pca_after,
  width    = 8,      # inches → 800 px at dpi=100
  height   = 6,      # inches → 600 px at dpi=100
  units    = "in",
  dpi      = 100
)
print("Batch Correction Complete")

##Projections
# auto-generate a named character vector of “prefixes”
patterns_list <- setNames(
  vapply(names(data_stes),
         function(nm) strsplit(nm, "_")[[1]][1],
         FUN.VALUE = character(1)),
  names(data_stes)
)



#Function for projection
project_X_to_Y_proteins <- function(X, Y, k = 8, d = 4){
  
  # X: Uncorrected matrix with proteins in the rows and samples in the columns
  # Y: Corrected matrix with proteins in the rows and samples in the columns
  # number of proteins in X should be larger than number of proteins in Y
  # k: the number of nearest neighbors to be considered
  # d: Number of PCs to use
  #Select all the proteins that need to be projected
  
  all_missing_proteins <- setdiff(rownames(X), rownames(Y))
  name_lengths <- nchar(all_missing_proteins)
  all_missing_proteins <- all_missing_proteins[name_lengths > 0]
  all_common_proteins <- intersect(rownames(X), rownames(Y))
  Y_proj <- Y
  
  #Calculate meaningful distance metric
  pr_X <- princomp(X)
  
  #We have a problme with d. When d is larger than the column of the data set we get an out of bound
  X_dist <- as.matrix(dist(pr_X$scores[,01:d], method = 'euclidean'))
  
  
  for (j in all_missing_proteins){
    
    protein_to_project <- j
    
    #Find the k closest neighbors that are also available in Y
    closest_proteins_index <- sort(X_dist[protein_to_project, all_common_proteins])[1:k]
    closest_proteins <- names(sort(X_dist[protein_to_project, all_common_proteins])[1:k])
    
    #Calculate the batch vector
    batch_vec <- X[closest_proteins,] - Y[closest_proteins,]
    
    #Project protein into Y batch space for each remaining neighbor
    projected_vec <- matrix(0, nrow = nrow(batch_vec), ncol = ncol(batch_vec))
    for (i in 1:nrow(batch_vec)){
      projected_vec[i,] <- X[protein_to_project,] - batch_vec[i,]
    }
    
    #Use the average projected vector
    average_proj <- matrix(colMeans(projected_vec), nrow = 1)
    
    # average_proj is a 1×#samples matrix
    rownames(average_proj) <- protein_to_project
    
    #Put the new protein into the batch corrected matrix
    Y_proj <- rbind(Y_proj, average_proj)
    
  }
  
  return(Y_proj)
  
}

# Iterate over the names of data_stes
for (dataset_name in names(data_stes)) {
  
  X_count <- as.matrix(data_stes[[dataset_name]][["normalized"]])
  X_count[is.na(X_count)] <- 0
  
  # Extract the relevant pattern for the current dataset name
  pattern <- patterns_list[[dataset_name]]  # Match pattern with the dataset name
  
  # Use grep to extract columns based on the pattern
  Y_count <- combat_exprs[, grep(pattern, colnames(combat_exprs), ignore.case = TRUE, value = TRUE)]
  
  # Project the data using your custom function
  data_stes[[dataset_name]][["projection"]] <- as.data.frame(project_X_to_Y_proteins(X_count, Y_count))
  
}

##Projection Matrix and Imputation
all_genes <- lapply(data_stes, function(x){
  return(rownames(x$projection))
})

all_genes <- Reduce(union, all_genes)

n <- names(data_stes)
mat_names <- unlist(lapply(n, function(x) strsplit(x, "_")[[1]][1]))
all_matrices <- lapply(data_stes, function(x){
  m <- x$projection[all_genes,]
  rownames(m) <- all_genes
  return(m)
})

for(i in 1:length(mat_names)){
  colnames(all_matrices[[i]]) <- paste0(mat_names[i],"-",colnames(all_matrices[[i]]))
}

all_matrices_proj <- Reduce(cbind, all_matrices)

print("Projection Complete")

#HeatMap
# Replace NA with a specific value to identify them later
all_matrices_mat_na <- all_matrices_proj %>% 
  mutate(across(everything(), ~replace(., is.na(.), -Inf)))

# Gather the data into long format for ggplot
long_data <- all_matrices_mat_na %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Column", value = "Value", -Row)

# Create a new variable to identify the NA and zero values
long_data <- long_data %>%
  mutate(
    ValueType = case_when(
      Value == -Inf ~ "NA",
      Value == 0 ~ "Zero",
      TRUE ~ "Other"
    )
  )

# Define custom colors
custom_colors <- c("NA" = "grey", "Zero" = "red", "Other" = "blue")

# Create the heatmap
heatmap_one <- ggplot(long_data, aes(x = Column, y = Row, fill = ValueType)) +
  geom_tile() +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Heatmap of NA and Zero Values", fill = "Value Type") +
  theme_minimal() +
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_blank())

ggsave(
  filename = file.path(plots_dir, "Heatmap of NA.png"),
  plot     = heatmap_one,
  width    = 8,      # inches → 800 px at dpi=100
  height   = 6,      # inches → 600 px at dpi=100
  units    = "in",
  dpi      = 100
)



missing_percentage <- rowSums(is.na(all_matrices_proj)) / ncol(all_matrices_proj) * 100

rows_to_impute <- missing_percentage <= 22
matrix_to_impute <- as.matrix(all_matrices_proj[rows_to_impute, ])

msnset_object <- MSnSet(exprs = matrix_to_impute)

imputed_result <- imputeLCMD::impute.QRILC(exprs(msnset_object))

imputed_matrix <- imputed_result[[1]]


# Verify
sum(is.na(imputed_matrix))   # Should be 0 if imputation was successful

all_tissue_types <- c()
for (l in names(data_stes)) {
  tissue_type <- data_stes[[l]][["meta"]][["tissue"]]
  all_tissue_types <- c(all_tissue_types, tissue_type)
}
all_tissue_types <- factor(all_tissue_types, levels = unique(all_tissue_types))


# Perform PCA - After Batch Correction
pca_after <- prcomp(t(imputed_matrix), scale.=TRUE)
pca_data_after <- data.frame(
  PC1 = pca_after$x[, 1], 
  PC2 = pca_after$x[, 2], 
  Batch = factor(batch_data$batch, levels = names(batch_names), labels = batch_names)
)
imputed_pca <- ggplot(pca_data_after, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.8) +  
  scale_color_manual(values = palette) +
  ggtitle("PCA Plot After Imputation 18%") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  labs(color = "Batch")

ggsave(
  filename = file.path(plots_dir, "PCA Plot After Imputation 18%.png"),
  plot     = imputed_pca,
  width    = 8,      # inches → 800 px at dpi=100
  height   = 6,      # inches → 600 px at dpi=100
  units    = "in",
  dpi      = 100
)



pca_data_after <- data.frame(
  PC1 = pca_after$x[, 1], 
  PC2 = pca_after$x[, 2], 
  Batch = all_tissue_types
)
tissue_impute <- ggplot(pca_data_after, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.8) +  
  scale_color_manual(values = palette) +  
  ggtitle("PCA Tissue 18% Imputation") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  labs(color = "Tissues")
ggsave(
  filename = file.path(plots_dir, "Tissue Plot After Imputation 18%.png"),
  plot     = tissue_impute,
  width    = 8,      # inches → 800 px at dpi=100
  height   = 6,      # inches → 600 px at dpi=100
  units    = "in",
  dpi      = 100
)



all_tissue_types <- c()
for (l in names(data_stes)) {
  tissue_type <- data_stes[[l]][["meta"]][["diet"]]
  all_tissue_types <- c(all_tissue_types, tissue_type)
}
all_tissue_types <- factor(all_tissue_types, levels = unique(all_tissue_types))

pca_data_after <- data.frame(
  PC1 = pca_after$x[, 1], 
  PC2 = pca_after$x[, 2], 
  Batch = all_tissue_types
)
ggplot(pca_data_after, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3, alpha = 0.8) +  
  scale_color_manual(values = palette) +  
  ggtitle("PCA Diet 18% Imputation") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  labs(color = "Diet")

# AI model training
## function

# 5) Save the final processed object back to Shiny
saveRDS(imputed_matrix, file = output_data_path)

# 6) Export the “last” ggplot to the plot file
#    (If you want a specific plot, either wrap it in ggsave() directly
#     or use last_plot() as below.)

ggsave(filename = output_plot_path,
       plot     = last_plot(),
       device   = "png",
       width    = 10,    # inches, adjust as you like
       height   = 7,
       units    = "in")

