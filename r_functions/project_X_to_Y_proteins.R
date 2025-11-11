project_X_to_Y_proteins <- function(X, Y, k = 2, d = 4) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  if (is.null(dim(X))) X <- matrix(X, ncol = 1, dimnames = list(names(X), "Sample1"))
  if (is.null(dim(Y))) Y <- matrix(Y, ncol = 1, dimnames = list(names(Y), "Sample1"))
  
  rn <- rownames(X)
  valid_rows <- !is.na(rn) & nzchar(rn)
  X <- X[valid_rows, , drop = FALSE]
  rownames(X) <- rn[valid_rows]
  
  common_cols <- intersect(colnames(X), colnames(Y))
  if (length(common_cols) == 0) {
    warning("No common sample columns between X and Y → returning Y as-is")
    return(as.data.frame(Y))
  }
  Xs <- X[, common_cols, drop = FALSE]
  Ys <- Y[, common_cols, drop = FALSE]
  
  if (ncol(Xs) != ncol(Ys)) stop("Projection aborted: mismatch in sample counts")
  
  all_missing_proteins <- setdiff(rownames(Xs), rownames(Ys))
  all_missing_proteins <- all_missing_proteins[nchar(all_missing_proteins) > 0]
  all_common_proteins  <- intersect(rownames(Xs), rownames(Ys))
  Y_proj <- Ys
  
  pr_X <- try(princomp(Xs), silent = TRUE)
  if (inherits(pr_X, "try-error") || is.null(pr_X$scores)) {
    # fallback: use rows as-is
    for (j in all_missing_proteins) {
      avg <- matrix(rowMeans(Xs[j, , drop = FALSE], na.rm = TRUE), nrow = 1)
      rownames(avg) <- j
      colnames(avg) <- colnames(Ys)
      Y_proj <- rbind(Y_proj, avg)
    }
    return(as.data.frame(Y_proj))
  }
  if (is.null(rownames(pr_X$scores))) rownames(pr_X$scores) <- rownames(Xs)
  d_use <- min(d, ncol(pr_X$scores))
  X_dist <- as.matrix(dist(pr_X$scores[, 1:d_use, drop = FALSE], method = "euclidean"))
  
  for (protein_to_project in all_missing_proteins) {
    if (!protein_to_project %in% rownames(X_dist)) next
    neighbors <- sort(X_dist[protein_to_project, all_common_proteins, drop = FALSE])[1:min(k, length(all_common_proteins))]
    closest_proteins <- names(neighbors)
    if (length(closest_proteins) == 0) next
    
    batch_vec <- Xs[closest_proteins, , drop = FALSE] - Ys[closest_proteins, , drop = FALSE]
    projected_vec <- matrix(0, nrow = nrow(batch_vec), ncol = ncol(batch_vec))
    for (ii in seq_len(nrow(batch_vec))) {
      projected_vec[ii, ] <- Xs[protein_to_project, ] - batch_vec[ii, ]
    }
    average_proj <- matrix(colMeans(projected_vec, na.rm = TRUE), nrow = 1, ncol = ncol(Ys))
    rownames(average_proj) <- protein_to_project
    colnames(average_proj) <- colnames(Ys)
    Y_proj <- rbind(Y_proj, average_proj)
  }
  return(as.data.frame(Y_proj))
}