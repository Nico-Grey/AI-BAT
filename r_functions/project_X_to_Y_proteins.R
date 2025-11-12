project_X_to_Y_proteins <- function(X, Y, k = 2, d = 4) {
  # X: Uncorrected matrix with proteins in rows and samples in columns
  # Y: Corrected matrix with proteins in rows and samples in columns
  # k: number of nearest neighbors
  # d: number of principal components to use
  
  # --- Input checks ---
  if (is.null(rownames(X)) || is.null(rownames(Y))) {
    stop("Both X and Y must have row names corresponding to proteins.")
  }
  
  all_missing_proteins <- setdiff(rownames(X), rownames(Y))
  all_missing_proteins <- all_missing_proteins[nchar(all_missing_proteins) > 0]
  all_common_proteins  <- intersect(rownames(X), rownames(Y))
  
  if (length(all_common_proteins) == 0) {
    stop("No common proteins between X and Y — cannot project.")
  }
  
  Y_proj <- Y
  
  # --- Handle single-sample case ---
  if (ncol(X) == 1 || ncol(Y) == 1) {
    message("Single-sample mode detected — skipping PCA and using direct differences.")
    
    for (j in all_missing_proteins) {
      if (!(j %in% rownames(X))) next
      
      # Find k nearest neighbors among common proteins (by absolute difference)
      diffs <- abs(X[j, 1] - X[all_common_proteins, 1])
      closest_proteins <- names(sort(diffs))[1:min(k, length(diffs))]
      
      # Batch difference (common proteins)
      batch_vec <- X[closest_proteins, 1, drop = FALSE] - Y[closest_proteins, 1, drop = FALSE]
      
      # Project missing protein
      projected_val <- X[j, 1] - mean(batch_vec, na.rm = TRUE)
      
      avg_proj <- matrix(projected_val, nrow = 1)
      rownames(avg_proj) <- j
      colnames(avg_proj) <- colnames(Y)
      
      Y_proj <- rbind(Y_proj, avg_proj)
    }
    
    return(Y_proj)
  }
  
  # --- Multi-sample case ---
  
  # Principal component analysis on X
  pr_X <- princomp(X)
  
  # Limit d to available components
  d <- min(d, ncol(pr_X$scores))
  
  # Compute Euclidean distances between proteins
  X_dist <- as.matrix(dist(pr_X$scores[, 1:d, drop = FALSE], method = "euclidean"))
  
  for (j in all_missing_proteins) {
    if (!(j %in% rownames(X_dist))) next
    
    # Find k closest neighbors among common proteins
    dist_vec <- X_dist[j, all_common_proteins]
    closest_proteins <- names(sort(dist_vec, decreasing = FALSE))[1:min(k, length(dist_vec))]
    
    # Calculate batch correction vector
    batch_vec <- X[closest_proteins, , drop = FALSE] - Y[closest_proteins, , drop = FALSE]
    
    # Project protein into Y batch space
    projected_vec <- matrix(0, nrow = nrow(batch_vec), ncol = ncol(batch_vec))
    for (i in 1:nrow(batch_vec)) {
      projected_vec[i, ] <- X[j, , drop = FALSE] - batch_vec[i, ]
    }
    
    # Average projected vector
    avg_proj <- matrix(colMeans(projected_vec, na.rm = TRUE), nrow = 1)
    rownames(avg_proj) <- j
    colnames(avg_proj) <- colnames(Y)
    
    # Append to Y_proj
    Y_proj <- rbind(Y_proj, avg_proj)
  }
  
  return(Y_proj)
}
