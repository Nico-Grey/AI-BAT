# general_functions.R

# This file includes utility functions used throughout the application.

# Function to load dataset
load_dataset <- function(file_path) {
  if (file.exists(file_path)) {
    return(readRDS(file_path))
  } else {
    stop("Dataset not found at the specified path.")
  }
}

# Function to preprocess data
preprocess_data <- function(data) {
  # Example preprocessing steps
  data <- na.omit(data)  # Remove missing values
  data <- data[data$expression > 0, ]  # Filter out non-positive expression values
  return(data)
}

# Function to generate summary statistics
generate_summary <- function(data) {
  summary_stats <- summary(data)
  return(summary_stats)
}

# Function to plot data
plot_data <- function(data, x_col, y_col) {
  library(ggplot2)
  ggplot(data, aes_string(x = x_col, y = y_col)) +
    geom_point() +
    theme_minimal() +
    labs(title = "Data Plot", x = x_col, y = y_col)
}