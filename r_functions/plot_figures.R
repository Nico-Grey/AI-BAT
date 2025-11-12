# ---- PCA plot ----
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

# ---- PCA browning score ----
plot_browning_score <- function(sample, ai_results, plot_dir) {
  df <- ai_results
  df$pca_browning_score_PC1_scaled <- scales::rescale(df$pca_browning_score_PC1, to = c(0, 100))
  
  # Output path
  out_path <- file.path(plot_dir, paste0("Browning_score_", sample, ".png"))
  
  # Open PNG device
  png(filename = out_path, width = 1200, height = 600, res = 150)
  
  plot(1, ylim=c(1,4), xlim=c(-30,100), type="n", axes=F, xlab="", ylab="")
  
  cols <- colorRampPalette(c("lightblue", "cyan", "orange"))(200)
  for(i in 1:200){
    rect(xleft = i * 0.5,xright = i * 0.5 + 0.5,ybottom = 1.5,ytop = 2,col=cols[i],border = NA)
  }
  
  text(x=5,y=1,label = "1 (White)",cex = 1,col="lightblue")
  text(x=95,y=1,label = "100 (Brown)",cex = 1,col="orange")
  text(x=-15,y=2.3,label = "Browning score:")
  points(x=df$pca_browning_score_PC1_scaled,
         y=rep(2.3,nrow(df)),
         pch = 16,
         col = "grey")
  
  sel <- df[df$X %in% sample, ]
  if (nrow(sel) > 0) {
    arrows(sel$pca_browning_score_PC1_scaled, 2.3,
           sel$pca_browning_score_PC1_scaled, 2.7,
           length = 0.1, lwd = 1.5, col = "black")
    text(sel$pca_browning_score_PC1_scaled, 3.0,
         labels = paste0(sel$X, "\n", round(sel$pca_browning_score_PC1_scaled, 2)),
         cex = 0.8, font = 2)
  }
  
  text(10, 3.5, "The browning score for your sample is:", cex = 1.2)
  dev.off()
}

