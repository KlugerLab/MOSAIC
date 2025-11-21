run_MOSAIC_two_modality <- function(rna_seurat, atac_seurat, assay1 = "SCT", assay2="ADT",sample_meta="sample_id",condition_meta="condition",block_method="cosine"){
  # Global scale
  rna_seurat <- ScaleData(rna_seurat,features=rownames(rna_seurat))
  atac_seurat <- ScaleData(atac_seurat,features=rownames(atac_seurat))
  
  rna.list <- SplitObject(rna_seurat, split.by = sample_meta)
  atac.list <- SplitObject(atac_seurat, split.by = sample_meta)
  
  dataset_list <- list()
  for (i in 1:length(rna.list)){
    gene_seurat <- rna.list[[i]]
    gene_matrix <- gene_seurat[[assay1]]$scale.data
    gene_matrix <- t(as.matrix(gene_matrix))
    
    atac_seurat <- atac.list[[i]]
    atac_matrix <- atac_seurat[[assay2]]$scale.data
    atac_matrix <- t(as.matrix(atac_matrix))
    
    dataset_list[[i]] <- list(gene_matrix, atac_matrix)
  }
  
  eigen_output_cosine_list <- list()
  r_cosine_list <- c()
  block_matrix_list <- list()
  for (i in 1:length(dataset_list)){
    print(paste0("Processing sample ",i))
    matrix_set <- dataset_list[[i]]
    
    if (block_method == "cosine"){
      eigen_output_cosine <- eigen_block_cosine_two(matrix_set[[1]],matrix_set[[2]])
    } else if (block_method == "gaussian"){
      print("gaussian")
      eigen_output_cosine <- eigen_block_gaussian_two(matrix_set[[1]],matrix_set[[2]])
    } else if (block_method == "emd"){
      print("emd")
      eigen_output_cosine <- eigen_block_emd_two(matrix_set[[1]],matrix_set[[2]])
    } else {
      print("block method not found")
    }
    
    r_cosine_list <- c(r_cosine_list,eigen_output_cosine[[2]])
    eigen_output_cosine_list[[i]] <- eigen_output_cosine[c(1,2)]
    block_matrix_list[[i]] <- eigen_output_cosine[[3]]
  }
  r <- round(median(r_cosine_list))
  print(paste0("r: ", r))
  names(eigen_output_cosine_list) <- names(rna.list)
  names(block_matrix_list) <- names(rna.list)
  
  sample_name <- names(eigen_output_cosine_list)
  
  
  sample_condition_map <- unique(data.frame(
    sample_id = rna_seurat[[sample_meta]],
    condition = rna_seurat[[condition_meta]]
  ))
  
  colnames(sample_condition_map) <- c("sample_id","condition")
  
  annotation <- data.frame(
    Condition = sample_condition_map$condition
  )
  rownames(annotation) <- sample_condition_map$sample_id
  
  annotation <- annotation[names(eigen_output_cosine_list), , drop = FALSE]
  
  
  normalized_r_eigenvectors <- list()
  for (i in 1:length(eigen_output_cosine_list)){
    r_to_use <- r_cosine_list[i]
    if (r_to_use == 1){
      r_to_use <- 2
    }
    if (block_method == "gaussian"){
      original_matrix <- eigen_output_cosine_list[[i]][[1]]$vectors[,2:r_to_use]
    } else {
      original_matrix <- eigen_output_cosine_list[[i]][[1]]$vectors[,1:r_to_use]
    }
    normalized_r_eigenvectors[[i]] <- original_matrix
  }
  
  
  names(normalized_r_eigenvectors) <- names(eigen_output_cosine_list)
  
  sample_projection_matrix <- list()
  for (i in 1:length(eigen_output_cosine_list)){
    r_to_use <- r_cosine_list[i]
    if (block_method == "gaussian"){
      U_U_transpose <- normalized_r_eigenvectors[[i]] %*% diag(eigen_output_cosine_list[[i]][[1]]$values[2:r_to_use]) %*% t(normalized_r_eigenvectors[[i]])
    } else {
      U_U_transpose <- normalized_r_eigenvectors[[i]] %*% diag(eigen_output_cosine_list[[i]][[1]]$values[1:r_to_use]) %*% t(normalized_r_eigenvectors[[i]])
    }
    sample_projection_matrix[[i]] <- U_U_transpose
  }
  
  sum_projection_matrix <- Reduce("+", sample_projection_matrix)
  V_result <- RSpectra::eigs_sym(sum_projection_matrix, 50, which = "LA")
  r2 <- find_elbow_kneedle(V_result$values[1:50])
  V_matrix <- V_result$vectors[,1:r2] %*% diag(V_result$values[1:r2])
  
  projected_combed_list <- list()
  for (i in 1:length(sample_projection_matrix)){
    projected_combed_list[[i]] <- block_matrix_list[[i]] %*% V_matrix
  }
  names(projected_combed_list) <- names(eigen_output_cosine_list)
  
  return(list(projected_combed_list=projected_combed_list,eigen_output_cosine_list=eigen_output_cosine_list,
              annotation=annotation,block_matrix_list=block_matrix_list,eigenvalues=V_result$values))
}

eigen_block_cosine_two <- function(matrix1,matrix2,whether_full=TRUE){
  A_B <- cosine_distance(matrix1,matrix2)
  dim1 <- dim(matrix1)[2]
  dim2 <- dim(matrix2)[2]
  block_matrix <- rbind(
    cbind(matrix(0, nrow=dim1, ncol=dim1), 
          A_B),
    cbind(t(A_B), 
          matrix(0, nrow=dim2, ncol=dim2))
  )
  if(whether_full){
    A_A <- cosine_distance(matrix1,matrix1)
    B_B <- cosine_distance(matrix2,matrix2)
    block_matrix <- rbind(
      cbind(A_A, 
            A_B),
      cbind(t(A_B), 
            B_B)
    )
  }
  print("Finish calculating coupling matrix")
  eigen_output <- eigs_sym(block_matrix, 50, which = "LA")
  print("Finish eigendecomposition")
  #r <- max(which(eigen_output$values[1:49]/eigen_output$values[2:50]>1.05))
  r <- find_elbow_kneedle(eigen_output$values)
  
  return(list(eigen_output,r,block_matrix))
}

find_elbow_kneedle <- function(eigenvalues) {
  # 1. Create the x and y coordinates for the curve.
  # We use the log of eigenvalues for better numerical stability and to emphasize the drop-off.
  for (i in 1:length(eigenvalues)){
    if (eigenvalues[i] <= 0){
      eigenvalues[i] <- eigenvalues[i-1]
    }
  }
  
  x <- 1:length(eigenvalues)
  
  y <- log(eigenvalues)
  
  # Handle cases with too few points to form a curve.
  if (length(x) < 3) {
    return(length(x))
  }
  
  # 2. Normalize the data to the [0, 1] range.
  # This prevents the result from being skewed by the scale of the axes.
  x_scaled <- (x - min(x)) / (max(x) - min(x))
  y_scaled <- (y - min(y)) / (max(y) - min(y))
  
  # 3. Calculate the perpendicular distance of each point from the line
  # connecting the first and last points of the normalized curve.
  distances <- perpendicular_distance(
    xp = x_scaled,
    yp = y_scaled,
    x1 = x_scaled[1],
    y1 = y_scaled[1],
    x2 = x_scaled[length(x_scaled)],
    y2 = y_scaled[length(y_scaled)]
  )
  
  # 4. The elbow is the index of the point with the maximum distance.
  elbow_index <- which.max(distances)
  
  return(elbow_index)
}


cosine_distance <- function(A, B) {
  A <- t(A)
  B <- t(B)
  A_norm <- sqrt(rowSums(A^2))
  B_norm <- sqrt(rowSums(B^2))
  
  A_norm[A_norm == 0] <- 1
  B_norm[B_norm == 0] <- 1
  
  A_normalized <- A / A_norm
  B_normalized <- B / B_norm
  
  similarity <- tcrossprod(A_normalized, B_normalized)
  return(similarity)
}

perpendicular_distance <- function(xp, yp, x1, y1, x2, y2) {
  # This function calculates the distance from each point (xp, yp) to the line
  # defined by the two points (x1, y1) and (x2, y2).
  
  # Line equation in the form Ax + By + C = 0
  A <- y1 - y2
  B <- x2 - x1
  C <- -A * x1 - B * y1
  
  # Calculate the perpendicular distance for each point using the formula:
  # |Ax_p + By_p + C| / sqrt(A^2 + B^2)
  distance <- abs(A * xp + B * yp + C) / sqrt(A^2 + B^2)
  return(distance)
}

plot_eigen <- function(eigenvalues){
  df <- data.frame(
    y = eigenvalues[1:50],
    x = 1:50
  )
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(size = 1) +
    theme_minimal() +
    labs(
      x = "n",
      y = "eigenvalues") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  return(p)
}

cal_peranova_per_feature_parallel_no_zero <- function(coembed_matrix_cosine_positive_list, n_sample, groups, n_cores = NULL,dist_method="cosine",num_permutation=999,whether_permanova=TRUE) {
  # Register parallel backend
  existing_connections <- showConnections()
  for(i in seq_len(nrow(existing_connections))) {
    try(close(getConnection(i)), silent = TRUE)
  }
  
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 1  # Leave 1 core free
  }
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  n_features <- dim(coembed_matrix_cosine_positive_list[[1]])[1]
  
  # Export variables and functions to workers
  clusterExport(cl, c("coembed_matrix_cosine_positive_list", "n_sample", "groups","cosine_distance","remove_zero_rows"),
                envir = environment())
  
  # Process each feature in parallel
  results <- foreach(idx = 1:n_features, .packages = c("vegan","stats","cluster","transport")) %dopar% {
    if (idx %% 100 == 0) {
      cat("Processing index:", idx, "\n")
    }
    
    combined_coembed_matrix <- matrix(0, nrow = n_sample, ncol = dim(coembed_matrix_cosine_positive_list[[1]])[2])
    for (i in 1:n_sample){
      combined_coembed_matrix[i,] <- coembed_matrix_cosine_positive_list[[i]][idx, ]
    }
    rownames(combined_coembed_matrix) <- names(coembed_matrix_cosine_positive_list)
    remove_result <- remove_zero_rows(combined_coembed_matrix)
    combined_coembed_matrix <- remove_result[[1]]
    remove_idx <- remove_result[[2]]
    
    # Compute similarity matrix
    dim_idx <- dim(combined_coembed_matrix)[1]
    similarity_matrix <- matrix(0, nrow = dim_idx, ncol = dim_idx)
    for (i in 1:dim_idx) {
      for (j in 1:dim_idx) {
        x <- combined_coembed_matrix[i, ]
        y <- combined_coembed_matrix[j, ]
        if (dist_method == "cosine"){
          similarity_matrix[i, j] <- cosine_distance(x, y)
        } else if (dist_method == "emd"){
          similarity_matrix[i, j] <- wasserstein1d(x, y, p = 1)
        } else {
          similarity_matrix[i, j] <- as.numeric(dist(rbind(x,y)))
        }
      }
    }
    rownames(similarity_matrix) <- rownames(combined_coembed_matrix)
    colnames(similarity_matrix) <- rownames(combined_coembed_matrix)
    
    if (dist_method == "cosine"){
      dissimilarity_matrix <- as.dist(1 - (similarity_matrix))
    } else {
      dissimilarity_matrix <- as.dist(similarity_matrix)
    }
    
    
    
    # Run PERMANOVA
    #dissimilarity_matrix <- as.dist(1 - similarity_matrix)
    if (whether_permanova){
      condition_numeric <- as.integer(as.factor(groups[!remove_idx]))
      silhouette_scores <- silhouette(condition_numeric, dissimilarity_matrix)
      average_silhouette_score <- mean(silhouette_scores[, "sil_width"])
      
      permanova_result <- adonis2(dissimilarity_matrix ~ groups[!remove_idx], permutations = num_permutation)
      
      list(
        pr = permanova_result$`Pr(>F)`[1],
        r2 = permanova_result$R2[1],
        similarity_matrix = similarity_matrix,
        permanova_result = permanova_result,
        group = groups[!remove_idx],
        F_statistics = permanova_result$`F`[1],
        average_silhouette_score = average_silhouette_score
      )
    } else {
      list(
        similarity_matrix = similarity_matrix,
        group = groups[!remove_idx]
      )
    }
  }
  
  stopCluster(cl)  # Clean up
  
  if (whether_permanova){
    # Extract results from the list
    pr_list <- sapply(results, function(x) x$pr)
    r2_list <- sapply(results, function(x) x$r2)
    similarity_matrix_list <- lapply(results, function(x) x$similarity_matrix)
    permanova_result_list <- lapply(results, function(x) x$permanova_result)
    group_list <- lapply(results, function(x) x$group)
    F_stats_list <- lapply(results, function(x) x$F_statistics)
    sil_score_list <- lapply(results, function(x) x$average_silhouette_score)
    
    
    return(list(
      pr_list = pr_list,
      r2_list = r2_list,
      similarity_matrix_list = similarity_matrix_list,
      permanova_result_list = permanova_result_list,
      group_list = group_list,
      F_stats_list = F_stats_list,
      sil_score_list = sil_score_list
    ))
  } else {
    similarity_matrix_list <- lapply(results, function(x) x$similarity_matrix)
    group_list <- lapply(results, function(x) x$group)
    
    return(list(
      similarity_matrix_list = similarity_matrix_list,
      group_list = group_list
    ))
  }
  
}

remove_zero_rows <- function(mat) {
  # Find which rows have all elements equal to zero
  zero_rows <- apply(mat, 1, function(row) all(row == 0))
  
  # Return the matrix without the all-zero rows
  return(list(mat[!zero_rows, , drop = FALSE],zero_rows))
}

calculate_matrix_residual <- function(X,Y){
  ssr <- sum((Y-X)^2)
  return(ssr)
}

cosine_similarity <- function(vector, matrix) {
  n_rows <- nrow(matrix)
  cos_sim <- matrix(0, nrow = n_rows, ncol = 1)
  if (all(vector == 0)) {
    return(cos_sim)
  }
  vector_norm <- sqrt(sum(vector^2))
  matrix_norms <- sqrt(rowSums(matrix^2))
  non_zero_rows <- !apply(matrix, 1, function(row) all(row == 0))
  if (any(non_zero_rows)) {
    matrix_subset <- matrix[non_zero_rows, , drop = FALSE]
    matrix_norms_subset <- matrix_norms[non_zero_rows]
    dot_products <- matrix_subset %*% vector
    cos_sim[non_zero_rows] <- dot_products / (matrix_norms_subset * vector_norm)
  }
  return(cos_sim)
}

cosine_vector_safe <- function(vector_a, vector_b) {
  if (all(vector_a == 0) || all(vector_b == 0)) {
    return(0)  # or NA, or any other value that makes sense for your application
  }
  
  # Normal calculation
  dot_product <- sum(vector_a * vector_b)
  norm_a <- sqrt(sum(vector_a^2))
  norm_b <- sqrt(sum(vector_b^2))
  
  return(dot_product / (norm_a * norm_b))
}


calculate_empirical_pvalue <- function(test_value, data_list, alternative = "greater") {
  if (alternative == "greater") {
    # P-value for right-tailed test (is test_value unusually large?)
    p_value <- mean(data_list >= test_value)
  } else if (alternative == "less") {
    # P-value for left-tailed test (is test_value unusually small?)
    p_value <- mean(data_list <= test_value)
  } else if (alternative == "two.sided") {
    # P-value for two-tailed test (is test_value unusually extreme in either direction?)
    # First find which tail the test value is in
    if (test_value > median(data_list)) {
      p_value <- 2 * mean(data_list >= test_value)
    } else {
      p_value <- 2 * mean(data_list <= test_value)
    }
    # Cap at 1 for extreme cases
    p_value <- min(p_value, 1)
  }
  
  return(p_value)
}


plot_mds_cluster_arrow <- function(similarity_matrix, title, cluster, size=1, 
                                   custom_colors = NULL, dis=FALSE, idx1=1, idx2=2,
                                   add_ellipse=TRUE, ellipse_level=0.66,
                                   add_arrow_from_to = NULL) { 
  
  # Convert similarity to dissimilarity
  if (!dis) {
    dissimilarity_matrix <- 1 - similarity_matrix
  } else {
    dissimilarity_matrix <- similarity_matrix
  }
  
  # Perform Classical MDS
  mds_result <- cmdscale(dissimilarity_matrix, k = 5) 
  
  umap_df <- data.frame(
    UMAP1 = mds_result[, idx1],
    UMAP2 = mds_result[, idx2],
    cluster = as.factor(cluster)
  )
  
  if (is.null(custom_colors)) {
    num_clusters <- length(unique(cluster))
    custom_colors <- RColorBrewer::brewer.pal(min(8, num_clusters), name = "Set2")
    names(custom_colors) <- levels(umap_df$cluster)
  } else {
    umap_df$cluster <- factor(umap_df$cluster, levels = names(custom_colors))
  }
  
  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster, fill = cluster)) +
    geom_point(size = size) +
    theme_minimal() +
    labs(title = title, x = paste0("MDS_", idx1), y = paste0("MDS_", idx2)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = custom_colors)
  
  # Add confidence ellipse if requested
  if (add_ellipse) {
    p <- p + 
      stat_ellipse(aes(fill = cluster), 
                   geom = "polygon", 
                   alpha = 0.2, 
                   level = ellipse_level, 
                   type = "t", 
                   linetype = 2, 
                   linewidth = 0.5) + 
      scale_fill_manual(values = custom_colors)
  }
  
  # --- Add Arrow Between Centroids ---
  if (!is.null(add_arrow_from_to) && length(add_arrow_from_to) == 2) {
    
    start_group <- add_arrow_from_to[1]
    end_group <- add_arrow_from_to[2]
    
    if (start_group %in% umap_df$cluster && end_group %in% umap_df$cluster) {
      
      # --- THIS IS THE FIX ---
      # Use base R 'aggregate' to calculate centroids.
      # This avoids all dplyr piping/masking conflicts.
      centroids_all <- aggregate(cbind(UMAP1, UMAP2) ~ cluster, data = umap_df, FUN = mean)
      
      # Subset the resulting data frame
      arrow_data_df <- centroids_all[centroids_all$cluster %in% c(start_group, end_group), ]
      
      # Check if both groups were found
      if (nrow(arrow_data_df) == 2) {
        
        # Create a data frame for the arrow segment
        arrow_data <- data.frame(
          x_start = arrow_data_df$UMAP1[arrow_data_df$cluster == start_group],
          y_start = arrow_data_df$UMAP2[arrow_data_df$cluster == start_group],
          x_end = arrow_data_df$UMAP1[arrow_data_df$cluster == end_group],
          y_end = arrow_data_df$UMAP2[arrow_data_df$cluster == end_group]
        )
        
        # Add the arrow layer to the plot
        p <- p + 
          geom_segment(
            data = arrow_data,
            aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
            arrow = arrow(type = "closed", length = unit(0.15, "inches")), # <-- Smaller arrowhead
            color = "darkgrey",                                            # <-- Darker color
            linewidth = 0.88,                                            # <-- Slightly thicker
            inherit.aes = FALSE  
          )
      } else {
        warning("Could not find centroids for both start and end groups.")
      }
    } else {
      warning("Arrow groups specified in 'add_arrow_from_to' not found in 'cluster' column.")
    }
  }
  # --- End of Arrow Code ---
  
  return(p)
}
