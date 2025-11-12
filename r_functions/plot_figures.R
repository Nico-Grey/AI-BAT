pca_plot <- function(input_data, meta_data, file_name, plot_dir) {
  # ensure plot_dir exists
  if (is.null(plot_dir) || !nzchar(plot_dir)) {
    plot_dir <- file.path(tempdir(), "plots")
  }
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Defensive checks
  if (is.null(input_data) || nrow(input_data) == 0) {
    message("pca_plot: input_data is NULL or has 0 rows; nothing to plot.")
    return(NULL)
  }
  
  # Remove rows with any NA
  input_data <- input_data[complete.cases(input_data), , drop = FALSE]
  if (nrow(input_data) == 0) {
    message("pca_plot: all rows removed by complete.cases(); nothing to plot.")
    return(NULL)
  }
  
  # PCA
  pr <- try(prcomp(t(input_data), center = TRUE, scale. = TRUE), silent = TRUE)
  if (inherits(pr, "try-error")) {
    message("pca_plot: prcomp() failed: ", pr)
    return(NULL)
  }
  pc_df <- as.data.frame(pr$x[, 1:2, drop = FALSE])
  pc_df$sample <- rownames(pc_df)
  
  # Merge meta if possible
  if (!is.null(meta_data) && "file_name" %in% colnames(meta_data)) {
    pc_df <- merge(pc_df, meta_data, by.x = "sample", by.y = "file_name", all.x = TRUE)
  } else {
    message("pca_plot: meta_data missing or has no 'file_name' column — continuing without merge.")
    pc_df$batch <- sapply(strsplit(as.character(pc_df$sample), "-"), `[`, 1)
  }
  
  # safe batch extraction if not present
  if (!"batch" %in% colnames(pc_df)) {
    pc_df$batch <- sapply(strsplit(as.character(pc_df$sample), "-"), `[`, 1)
  }
  
  # Plot
  p1 <- ggplot(pc_df, aes(x = PC1, y = PC2, color = batch)) +
    geom_point() +
    ggtitle(file_name) +
    theme_minimal()
  
  # Save (and verify)
  out_path <- file.path(plot_dir, paste0(file_name, ".png"))
  # message("pca_plot: saving to ", out_path)
  tryCatch({
    ggsave(filename = out_path, plot = p1, width = 8, height = 6, dpi = 150, device = "png")
  }, error = function(e) {
    message("pca_plot: ggsave failed: ", e$message)
  })
  
  # Verify
  if (file.exists(out_path)) {
    # message("pca_plot: saved OK: ", out_path)
    return(normalizePath(out_path))
  } else {
    message("pca_plot: save did not produce a file. plot_dir contents: ",
            paste(list.files(plot_dir), collapse = ", "))
    return(NULL)
  }
}
