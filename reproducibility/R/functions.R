# Helper functions shared by the analysis scripts.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(RSpectra)
  library(ggplot2)
  library(RColorBrewer)
  library(cluster)
  library(vegan)
  library(foreach)
  library(doParallel)
  library(Rcpp)
})

`%ni%` <- Negate(`%in%`)

plot_mds_cluster <- function(similarity_matrix,title,cluster,size=1, custom_colors = NULL,dis=FALSE,idx1=1,idx2=2) {

  
  # similarity to dissimilarity
  if (!dis) {
    dissimilarity_matrix <- 1 - similarity_matrix} else {
    dissimilarity_matrix <- similarity_matrix
    }
  
  
  # classical MDS
  mds_result <- cmdscale(dissimilarity_matrix, k = 5)  # 2D projection

  
  umap_df <- data.frame(
    UMAP1 = mds_result[,idx1],
    UMAP2 = mds_result[,idx2],
    cluster = as.factor(cluster)
  )
  
  if (is.null(custom_colors)) {
    num_clusters <- length(unique(cluster))
    custom_colors <- RColorBrewer::brewer.pal(min(8, num_clusters), name = "Set2")
    names(custom_colors) <- levels(umap_df$cluster)
  } else {
    umap_df$cluster <- factor(umap_df$cluster, levels = names(custom_colors))
  }

  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
    geom_point(size = size) +
    theme_minimal() +
    labs(title = title, x = paste0("MDS_",idx1), y = paste0("MDS_",idx2)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = custom_colors)
  return(p)
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

eigen_block_cosine_three <- function(TF_matrix,peak_matrix,gene_matrix){
    TF_peak <- cosine_distance(TF_matrix,peak_matrix)
    TF_gene <- cosine_distance(TF_matrix,gene_matrix)
    peak_gene <- cosine_distance(peak_matrix,gene_matrix)
    TF_TF <- cosine_distance(TF_matrix,TF_matrix)
    peak_peak <- cosine_distance(peak_matrix,peak_matrix)
    gene_gene <- cosine_distance(gene_matrix,gene_matrix)
    
    
    block_matrix <- rbind(
      cbind(TF_TF,
            TF_peak,
            TF_gene),

      cbind(t(TF_peak),
            peak_peak,
            peak_gene),

      cbind(t(TF_gene),
            t(peak_gene),
            gene_gene)
    )
    
    print("finish block")
    eigen_output <- eigs_sym(block_matrix, 50, which = "LA")
    r <- find_elbow_kneedle(eigen_output$values)
    
    return(list(eigen_output,r,block_matrix))
}

eigen_block_cosine_one <- function(matrix){
  block_matrix <- cosine_distance(matrix,matrix)
  print("finish calculating cosine")
  eigen_output <- eigs_sym(block_matrix, 50, which = "LA")
  print("finish eigendecomposition")
  r <- find_elbow_kneedle(eigen_output$values)
  return(list(eigen_output,r,block_matrix))
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
  print("finish calculating cosine")
  eigen_output <- eigs_sym(block_matrix, 50, which = "LA")
  print("finish eigendecomposition")
  r <- find_elbow_kneedle(eigen_output$values)
  
  return(list(eigen_output,r,block_matrix))
}

gaussian_kernel_distance <- function(A, B, sigma) {
    A <- t(A)
    B <- t(B)
    A_norm <- rowSums(A^2)
    B_norm <- rowSums(B^2)
    squared_dist <- outer(A_norm, rep(1, nrow(B))) + 
      outer(rep(1, nrow(A)), B_norm) - 
      2 * tcrossprod(A, B)
    K <- exp(-squared_dist / (2 * sigma^2))
    return(K)
}

choose_sigma <- function(A, B) {
  A <- t(A)
  B <- t(B)
  A_norm <- rowSums(A^2)
  B_norm <- rowSums(B^2)
  squared_dist <- outer(A_norm, rep(1, nrow(B))) + 
    outer(rep(1, nrow(A)), B_norm) - 
    2 * tcrossprod(A, B)
  squared_dist[squared_dist < 0] <- 0
  dist <- sqrt(squared_dist)
  return(median(dist))
}

eigen_block_gaussian_three <- function(TF_matrix,peak_matrix,gene_matrix){
    sigma <- choose_sigma(TF_matrix,peak_matrix)
    TF_peak <- gaussian_kernel_distance(TF_matrix,peak_matrix,sigma)
    print(paste0("finish TF_peak, sigma=",sigma))
    sigma <- choose_sigma(TF_matrix,gene_matrix)
    TF_gene <- gaussian_kernel_distance(TF_matrix,gene_matrix,sigma)
    print(paste0("finish TF_gene, sigma=",sigma))
    sigma <- choose_sigma(peak_matrix,gene_matrix)
    peak_gene <- gaussian_kernel_distance(peak_matrix,gene_matrix,sigma)
    print(paste0("finish peak_gene, sigma=",sigma))
    
    sigma <- choose_sigma(TF_matrix,TF_matrix)
    TF_TF <- gaussian_kernel_distance(TF_matrix,TF_matrix,sigma)
    print(paste0("finish TF_TF, sigma=",sigma))
    sigma <- choose_sigma(gene_matrix,gene_matrix)
    gene_gene <- gaussian_kernel_distance(gene_matrix,gene_matrix,sigma)
    print(paste0("finish gene_gene, sigma=",sigma))
    sigma <- choose_sigma(peak_matrix,peak_matrix)
    peak_peak <- gaussian_kernel_distance(peak_matrix,peak_matrix,sigma)
    print(paste0("finish peak_peak, sigma=",sigma))
    
    
    block_matrix <- rbind(
      cbind(TF_TF,
            TF_peak,
            TF_gene),

      cbind(t(TF_peak),
            peak_peak,
            peak_gene),

      cbind(t(TF_gene),
            t(peak_gene),
            gene_gene)
    )
    
    print("finish block")
    eigen_output <- eigs_sym(block_matrix, 50, which = "LA")
    r <- find_elbow_kneedle(eigen_output$values)
    
    return(list(eigen_output,r,block_matrix))
}

eigen_block_gaussian_two <- function(TF_matrix,peak_matrix){
  sigma <- choose_sigma(TF_matrix,peak_matrix)
  TF_peak <- gaussian_kernel_distance(TF_matrix,peak_matrix,sigma)
  
  sigma <- choose_sigma(TF_matrix,TF_matrix)
  TF_TF <- gaussian_kernel_distance(TF_matrix,TF_matrix,sigma)
  
  sigma <- choose_sigma(peak_matrix,peak_matrix)
  peak_peak <- gaussian_kernel_distance(peak_matrix,peak_matrix,sigma)
  
  print(paste0("finish TF_peak, sigma=",sigma))
  
  block_matrix <- rbind(
    cbind(TF_TF,
          TF_peak),
    
    cbind(t(TF_peak),
          peak_peak)
  )
  
  print("finish block")
  eigen_output <- eigs_sym(block_matrix, 50, which = "LA")
  r <- find_elbow_kneedle(eigen_output$values)
  
  return(list(eigen_output,r,block_matrix))
}

plot_eigen <- function(eigenvalues){
  df <- data.frame(
    y = eigenvalues[1:100],
    x = 1:100
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
  zero_rows <- apply(mat, 1, function(row) all(row == 0))
  
  return(list(mat[!zero_rows, , drop = FALSE],zero_rows))
}

calculate_matrix_residual <- function(X,Y){
  ssr <- sum((Y-X)^2)
  return(ssr)
}

get_genes_expressed_across_samples <- function(seurat_obj, sample_id_col, min_expr = 0) {
  if (!sample_id_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", sample_id_col, "not found in Seurat object metadata"))
  }
  
  sample_ids <- unique(seurat_obj@meta.data[[sample_id_col]])
  
  genes_by_sample <- list()
  
  counts <- GetAssayData(seurat_obj, slot = "counts")
  
  # features with total count >= min_expr in every sample
  for (sample in sample_ids) {
    cells_in_sample <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[sample_id_col]] == sample, ])
    
    sample_expr <- rowSums(counts[, cells_in_sample, drop = FALSE])
    
    expressed_genes <- names(sample_expr[sample_expr >= min_expr])
    
    genes_by_sample[[sample]] <- expressed_genes
  }
  
  genes_expressed_across_all <- Reduce(intersect, genes_by_sample)
  
  return(genes_expressed_across_all)
}

eigen_block_emd_two <- function(TF_matrix,peak_matrix,cal_r=50){
  TF_TF <- pairwise_normalized_emd(TF_matrix,TF_matrix)
  TF_peak <- pairwise_normalized_emd(TF_matrix,peak_matrix)
  peak_peak <- pairwise_normalized_emd(peak_matrix,peak_matrix)

  print(paste0("finish TF_peak"))
  
  block_matrix <- rbind(
    cbind(TF_TF,
          TF_peak),
    
    cbind(t(TF_peak),
          peak_peak)
  )
  
  print("finish block")
  eigen_output <- eigs_sym(block_matrix, cal_r, which = "LA")
  
  r <- find_elbow_kneedle(eigen_output$values)
  return(list(eigen_output,r,block_matrix))
}

eigen_block_emd_three <- function(TF_matrix,peak_matrix,gene_matrix){

  TF_peak <- pairwise_normalized_emd(TF_matrix,peak_matrix)
  print(paste0("finish TF_peak"))
  
  TF_gene <- pairwise_normalized_emd(TF_matrix,gene_matrix)
  print(paste0("finish TF_gene"))
  
  peak_gene <- pairwise_normalized_emd(peak_matrix,gene_matrix)
  print(paste0("finish peak_gene"))
  
  TF_TF <- pairwise_normalized_emd(TF_matrix,TF_matrix)
  print(paste0("finish TF_TF"))
  
  gene_gene <- pairwise_normalized_emd(gene_matrix,gene_matrix)
  print(paste0("finish gene_gene"))
  
  peak_peak <- pairwise_normalized_emd(peak_matrix,peak_matrix)
  print(paste0("finish peak_peak"))
  
  
  block_matrix <- rbind(
    cbind(TF_TF,
          TF_peak,
          TF_gene),
    
    cbind(t(TF_peak),
          peak_peak,
          peak_gene),
    
    cbind(t(TF_gene),
          t(peak_gene),
          gene_gene)
  )
  
  
  print("finish block")
  eigen_output <- eigs_sym(block_matrix, 50, which = "LA")
  r <- find_elbow_kneedle(eigen_output$values)
  
  return(list(eigen_output,r,block_matrix))
}

pairwise_normalized_emd <- function(A, B, k = 100) {
  if (is.vector(A)) A <- as.matrix(A)
  if (is.vector(B)) B <- as.matrix(B)

  # Precompute quantiles (k points)
  probs <- seq(0, 1, length.out = k)
  quantile_fun <- function(X) {
    apply(X, 2, function(col) quantile(col, probs = probs, type = 7, names = FALSE))
  }
  A_quant <- quantile_fun(A)  # k x p1 matrix
  B_quant <- quantile_fun(B)  # k x p2 matrix

  sds_A <- pmax(apply(A, 2, sd), 1e-9)
  sds_B <- pmax(apply(B, 2, sd), 1e-9)

  aa <- normalized_emd_cpp(A_quant, B_quant, sds_A, sds_B)
  return(exp(-aa))
}

find_elbow_kneedle <- function(eigenvalues) {
  # log eigenvalue curve
  for (i in 1:length(eigenvalues)){
    if (eigenvalues[i] <= 0){
      eigenvalues[i] <- eigenvalues[i-1]
    }
  }
  
  x <- 1:length(eigenvalues)
  
  y <- log(eigenvalues)
  
  if (length(x) < 3) {
    return(length(x))
  }
  
  # normalize both axes to [0, 1]
  x_scaled <- (x - min(x)) / (max(x) - min(x))
  y_scaled <- (y - min(y)) / (max(y) - min(y))
  
  # distance of each point to the chord between the first and last points
  distances <- perpendicular_distance(
    xp = x_scaled,
    yp = y_scaled,
    x1 = x_scaled[1],
    y1 = y_scaled[1],
    x2 = x_scaled[length(x_scaled)],
    y2 = y_scaled[length(y_scaled)]
  )
  
  # elbow = point with the largest distance (at least 2)
  elbow_index <- which.max(distances)
  output <- max(elbow_index,2)
  
  return(output)
}

perpendicular_distance <- function(xp, yp, x1, y1, x2, y2) {
  # distance from (xp, yp) to the line through (x1, y1) and (x2, y2)
  
  A <- y1 - y2
  B <- x2 - x1
  C <- -A * x1 - B * y1
  
  distance <- abs(A * xp + B * yp + C) / sqrt(A^2 + B^2)
  return(distance)
}


cppFunction('
NumericMatrix normalized_emd_cpp(
    NumericMatrix A_quant, NumericMatrix B_quant, 
    NumericVector sds_A, NumericVector sds_B
) {
  int k = A_quant.nrow();
  int p1 = A_quant.ncol();
  int p2 = B_quant.ncol();
  NumericMatrix D(p1, p2);
  
  for (int i = 0; i < p1; i++) {
    for (int j = 0; j < p2; j++) {
      double total = 0.0;
      for (int q = 0; q < k; q++) {
        total += std::abs(A_quant(q, i) - B_quant(q, j));
      }
      double emd_raw = total / k;
      double pooled_sd = (sds_A[i] + sds_B[j]) / 2 + 1e-9;
      D(i, j) = emd_raw / pooled_sd;
    }
  }
  return D;
}
')

