# Supplemental Figure S2: null calibration of the DC test under sample-specific heterogeneity (simulation)
# Output:   results/S2/null_calib_final/null_calib_results.rds, typeI_summary.csv
# Usage:    Rscript fig3/figS2_null_calibration.R --n_rep 10 --sigma 0.1,0.5,1.0

source("config.R")
source("R/functions.R")

## Global null; condition labels permuted across samples, per-feature null

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(irlba)
  library(RSpectra)
  library(pracma)
  library(coop)
  library(foreach)
  library(doParallel)
  library(vegan)
  library(ggplot2)
  library(dplyr)
})

opt_list <- list(
  make_option("--out_dir", type = "character",
              default = file.path(RESULTS_DIR, "S2", "null_calib_final")),
  make_option("--n_rep", type = "integer", default = 1L),
  make_option("--sigma", type = "character", default = "0,0.1,0.25,0.5,1.0",
              help = "Comma-separated sample_noise_level values; 0.5 is the Fig 3 setting"),
  make_option("--B", type = "integer", default = 999L, help = "Sample-label permutations"),
  make_option("--alpha", type = "double", default = 0.05),
  make_option("--n_cores", type = "integer", default = N_CORES),
  make_option("--summary_only", action = "store_true", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = opt_list))
sigma_levels <- as.numeric(strsplit(opt$sigma, ",")[[1]])
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

## Simulator and MOSAIC pipeline (as in the Figure 3 benchmarks)

generate_svd_simulated_data <- function(
    simulation_type = "DC", # "DC" for Differential Connectivity, "DE" for Differential Expression
    add_batch = TRUE,
    de_prop = 0.2,
    de_fc = 2.5, # Fold change for DE simulation
    seed = 1,
    whether_rescale = TRUE,
    sample_noise_level = 0.5,
    dia_matrix = diag(seq(5, 1, length.out = r))
) {
  
  set.seed(seed)
  
  # Basic parameters
  n_samples <- 20
  n_conditions <- 2
  n_batch <- 4
  cells_per_sample <- 1000
  r <- 10
  px <- 1000
  py <- 800
  pz <- 600
  batch_shift_mean <- 10
  batch_shift_sd <- 0.2
  batch_scale_mean <- 1.2
  batch_scale_sd <- 0.1
  n_clusters <- 10
  
  sample_metadata <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    condition = rep(c("A", "B"), each = n_samples / n_conditions),
    batch = rep(c("batch_1", "batch_2", "batch_3", "batch_4"), each = n_samples / n_batch)
  )
  
  # Global U (cell embeddings)
  cluster_centers_U <- matrix(rnorm(n_clusters * r), nrow = n_clusters, ncol = r)
  
  U_list <- list()
  cell_cluster_assignments_list <- list()
  for (j in 1:n_samples) {
    n_cells <- cells_per_sample + sample(-200:200, 1)
    cluster_assignments <- sample(1:n_clusters, n_cells, replace = TRUE)
    initial_points <- matrix(0, nrow = n_cells, ncol = r)
    for (i in 1:n_cells) {
      cluster <- cluster_assignments[i]
      initial_points[i, ] <- cluster_centers_U[cluster, ] + rnorm(r, sd = 0.1)
    }
    
    U <- initial_points
    U_list[[j]] <- U
    cell_cluster_assignments_list[[j]] <- cluster_assignments
  }
  
  # Global V (feature loadings), same cluster centres
  set.seed(seed)
  cluster_centers_V <- cluster_centers_U
  
  create_global_V <- function(p, n_clusters, r, cluster_centers) {
    cluster_assignments <- sample(1:n_clusters, p, replace = TRUE)
    initial_points <- matrix(0, nrow = p, ncol = r)
    for (i in 1:p) {
      cluster <- cluster_assignments[i]
      initial_points[i, ] <- cluster_centers[cluster, ] + rnorm(r, sd = 0.1)
    }
    # Return the loading matrix and cluster assignments
    list(V_global = initial_points, feature_clusters = cluster_assignments)
  }
  
  V_x_data <- create_global_V(px, n_clusters, r, cluster_centers_V)
  Vx_global <- V_x_data$V_global
  x_feature_cluster_assignments <- V_x_data$feature_clusters
  
  V_y_data <- create_global_V(py, n_clusters, r, cluster_centers_V)
  Vy_global <- V_y_data$V_global
  y_feature_cluster_assignments <- V_y_data$feature_clusters
  
  V_z_data <- create_global_V(pz, n_clusters, r, cluster_centers_V)
  Vz_global <- V_z_data$V_global
  z_feature_cluster_assignments <- V_z_data$feature_clusters
  
  # Sample-specific V: global V plus noise
  create_sample_V_list <- function(V_global, n_samples) {
    V_list <- list()
    p <- nrow(V_global)
    r <- ncol(V_global)
    for (i in 1:n_samples) {
      V_sample <- V_global + matrix(rnorm(p * r, 0, sample_noise_level), nrow = p, ncol = r)
      V_list[[i]] <- V_sample
    }
    return(V_list)
  }
  
  Vx_list <- create_sample_V_list(Vx_global, n_samples)
  Vy_list <- create_sample_V_list(Vy_global, n_samples)
  Vz_list <- create_sample_V_list(Vz_global, n_samples)
  
  # DE / DC features
  set.seed(seed)
  de_features_x <- sample(1:px, round(px * de_prop))
  de_features_y <- sample(1:py, round(py * de_prop))
  de_features_z <- sample(1:pz, round(pz * de_prop))
  
  # Initialize final data lists
  X_list <- list()
  Y_list <- list()
  Z_list <- list()
  
  # Define diagonal sigma matrices for each modality
  identical_sigma <- dia_matrix
  Sx <- identical_sigma
  Sy <- identical_sigma
  Sz <- identical_sigma
  
  # DC or DE signal
  
  if (simulation_type == "DC") {
    # DC: rewire selected features in condition B
    
    # 1. Map each selected feature to a feature from another cluster
    constrained_shuffle <- function(feature_indices, cluster_assignments) {
      shuffled_indices <- feature_indices
      for (i in seq_along(feature_indices)) {
        original_idx <- feature_indices[i]
        original_cluster <- cluster_assignments[original_idx]
        
        # Create a pool of candidates from different clusters
        candidate_pool <- feature_indices[cluster_assignments[feature_indices] != original_cluster]
        
        # If no candidates, fall back to any other feature
        if (length(candidate_pool) == 0) {
          candidate_pool <- feature_indices[feature_indices != original_idx]
        }
        
        shuffled_indices[i] <- sample(candidate_pool, 1)
      }
      return(shuffled_indices)
    }
    
    shuffled_cols_x <- constrained_shuffle(de_features_x, x_feature_cluster_assignments)
    shuffled_cols_y <- constrained_shuffle(de_features_y, y_feature_cluster_assignments)
    shuffled_cols_z <- constrained_shuffle(de_features_z, z_feature_cluster_assignments)
    
    # Apply permutation to V matrices for Condition B samples
    Vx_list_de <- Vx_list
    Vy_list_de <- Vy_list
    Vz_list_de <- Vz_list
    
    for (i in 1:n_samples) {
      if (sample_metadata$condition[i] == "B") {
        Vx_list_de[[i]][de_features_x, ] <- Vx_list[[i]][shuffled_cols_x, ]
        Vy_list_de[[i]][de_features_y, ] <- Vy_list[[i]][shuffled_cols_y, ]
        Vz_list_de[[i]][de_features_z, ] <- Vz_list[[i]][shuffled_cols_z, ]
      }
    }
    
    # 2. Generate intermediate expression matrices
    for (i in 1:n_samples) {
      X_list[[i]] <- U_list[[i]] %*% Sx %*% t(Vx_list_de[[i]])
      Y_list[[i]] <- U_list[[i]] %*% Sy %*% t(Vy_list_de[[i]])
      Z_list[[i]] <- U_list[[i]] %*% Sz %*% t(Vz_list_de[[i]])
    }
    
    # 3. Rescale condition B means to those of a condition A sample
    
    if (whether_rescale){
      target_means_x <- colMeans(U_list[[1]] %*% Sx %*% t(Vx_list[[1]]))[de_features_x]
      target_means_y <- colMeans(U_list[[1]] %*% Sy %*% t(Vy_list[[1]]))[de_features_y]
      target_means_z <- colMeans(U_list[[1]] %*% Sz %*% t(Vz_list[[1]]))[de_features_z]
      
      for (i in 1:n_samples) {
        if (sample_metadata$condition[i] == "B") {
          # Rescale X
          current_means_x <- colMeans(X_list[[i]])[de_features_x]
          diff_x <- target_means_x - current_means_x
          X_list[[i]][, de_features_x] <- sweep(X_list[[i]][, de_features_x], 2, diff_x, "+")
          
          # Rescale Y
          current_means_y <- colMeans(Y_list[[i]])[de_features_y]
          diff_y <- target_means_y - current_means_y
          Y_list[[i]][, de_features_y] <- sweep(Y_list[[i]][, de_features_y], 2, diff_y, "+")
          
          # Rescale Z
          current_means_z <- colMeans(Z_list[[i]])[de_features_z]
          diff_z <- target_means_z - current_means_z
          Z_list[[i]][, de_features_z] <- sweep(Z_list[[i]][, de_features_z], 2, diff_z, "+")
        }
      }
    }
    
  } else if (simulation_type == "DE") {
    # DE: fold change in condition B
    
    # 1. Generate all data without any DE signal first
    for (i in 1:n_samples) {
      X_list[[i]] <- U_list[[i]] %*% Sx %*% t(Vx_list[[i]])
      Y_list[[i]] <- U_list[[i]] %*% Sy %*% t(Vy_list[[i]])
      Z_list[[i]] <- U_list[[i]] %*% Sz %*% t(Vz_list[[i]])
    }
    
    # 2. Apply a multiplicative fold change for Condition B samples
    for (i in 1:n_samples) {
      if (sample_metadata$condition[i] == "B") {
        X_list[[i]][, de_features_x] <- X_list[[i]][, de_features_x] * de_fc
        Y_list[[i]][, de_features_y] <- Y_list[[i]][, de_features_y] * de_fc
        Z_list[[i]][, de_features_z] <- Z_list[[i]][, de_features_z] * de_fc
      }
    }
    
  } else {
    stop("simulation_type must be either 'DC' or 'DE'")
  }
  
  # Batch effects
  if (add_batch) {
    for (i in 1:n_samples) {
      
      batch_num <- as.numeric(sub("batch_", "", sample_metadata$batch[i]))
      
      # Modality X: per-feature scale and shift
      batch_scale_x_vec <- rnorm(px, batch_scale_mean, batch_scale_sd)
      batch_shift_x_vec <- rnorm(px, batch_shift_mean, batch_shift_sd) * (batch_num / n_batch)
      
      X_list[[i]] <- sweep(X_list[[i]], 2, batch_scale_x_vec, `*`) 
      X_list[[i]] <- sweep(X_list[[i]], 2, batch_shift_x_vec, `+`) 
      
      # Modality Y
      batch_scale_y_vec <- rnorm(py, batch_scale_mean, batch_scale_sd)
      batch_shift_y_vec <- rnorm(py, batch_shift_mean, batch_shift_sd) * (batch_num / n_batch)
      
      Y_list[[i]] <- sweep(Y_list[[i]], 2, batch_scale_y_vec, `*`) 
      Y_list[[i]] <- sweep(Y_list[[i]], 2, batch_shift_y_vec, `+`) 
      
      # Modality Z
      batch_scale_z_vec <- rnorm(pz, batch_scale_mean, batch_scale_sd)
      batch_shift_z_vec <- rnorm(pz, batch_shift_mean, batch_shift_sd) * (batch_num / n_batch)
      
      Z_list[[i]] <- sweep(Z_list[[i]], 2, batch_scale_z_vec, `*`) 
      Z_list[[i]] <- sweep(Z_list[[i]], 2, batch_shift_z_vec, `+`) 
    }
  }
  
  
  # Cell x feature matrices and metadata
  create_modality_matrix <- function(modality_matrix_list, feature_prefix, cell_cluster_assignments_list,
                                     feature_cluster_assignment, de_features) {
    n_feature <- dim(modality_matrix_list[[1]])[2]
    n_samples <- length(modality_matrix_list)
    
    all_cells <- c()
    all_features <- paste0(feature_prefix, "-", 1:n_feature)
    all_data <- list()
    all_metadata <- list()
    
    cell_offset <- 0
    
    for (i in 1:n_samples) {
      n_cell_sample <- dim(modality_matrix_list[[i]])[1]
      mod_matrix <- modality_matrix_list[[i]]
      
      sample_cells <- paste0("sample", i, "_cell_", 1:n_cell_sample)
      all_cells <- c(all_cells, sample_cells)
      
      sample_meta <- data.frame(
        cell_id = sample_cells,
        sample_id = rep(sample_metadata$sample_id[i], n_cell_sample),
        condition = rep(sample_metadata$condition[i], n_cell_sample),
        batch = rep(sample_metadata$batch[i], n_cell_sample),
        true_cell_type = cell_cluster_assignments_list[[i]]
      )
      all_metadata[[i]] <- sample_meta
      
      sparse_matrix <- as(mod_matrix, "dgCMatrix")
      all_data[[i]] <- sparse_matrix
      
      cell_offset <- cell_offset + n_cell_sample
    }
    combined_sparse <- do.call(rbind, all_data)
    
    features_df <- data.frame(
      feature_id = all_features,
      feature_cluster = feature_cluster_assignment,
      is_de = ((1:length(all_features)) %in% de_features)
    )
    
    combined_metadata <- do.call(rbind, all_metadata)
    
    rownames(combined_sparse) <- all_cells
    colnames(combined_sparse) <- all_features
    rownames(features_df) <- all_features
    rownames(combined_metadata) <- all_cells
    
    return(list(
      matrix = combined_sparse,
      features = features_df,
      cells = combined_metadata
    ))
  }
  
  X_data <- create_modality_matrix(X_list, "feature-x", cell_cluster_assignments_list, x_feature_cluster_assignments, de_features_x)
  Y_data <- create_modality_matrix(Y_list, "feature-y", cell_cluster_assignments_list, y_feature_cluster_assignments, de_features_y)
  Z_data <- create_modality_matrix(Z_list, "feature-z", cell_cluster_assignments_list, z_feature_cluster_assignments, de_features_z)
  
  # Seurat objects
  seurat_X <- CreateSeuratObject(
    counts = t(X_data$matrix), 
    data = t(X_data$matrix), 
    meta.data = X_data$cells,
    assay = "originalexp"
  )
  seurat_X[["originalexp"]] <- AddMetaData(
    object = seurat_X[["originalexp"]],
    metadata = X_data$features
  )
  
  seurat_Y <- CreateSeuratObject(
    counts = t(Y_data$matrix),
    data = t(Y_data$matrix),
    meta.data = Y_data$cells,
    assay = "originalexp"
  )
  seurat_Y[["originalexp"]] <- AddMetaData(
    object = seurat_Y[["originalexp"]],
    metadata = Y_data$features
  )
  
  seurat_Z <- CreateSeuratObject(
    counts = t(Z_data$matrix),
    data = t(Z_data$matrix),
    meta.data = Z_data$cells,
    assay = "originalexp"
  )
  seurat_Z[["originalexp"]] <- AddMetaData(
    object = seurat_Z[["originalexp"]],
    metadata = Z_data$features
  )
  
  # Return a list of all relevant objects for analysis
  return(list(
    seurat_X = seurat_X,
    seurat_Y = seurat_Y,
    seurat_Z = seurat_Z,
    X_list = X_list,
    Y_list = Y_list,
    Z_list = Z_list,
    sample_metadata = sample_metadata
  ))
}

run_our_method_three <- function(rna_seurat, atac_seurat,adt_seurat, sample_meta="sample_id",condition_meta="condition",block_method="cosine"){
  # Global scale
  rna_seurat <- ScaleData(rna_seurat,features=rownames(rna_seurat))
  atac_seurat <- ScaleData(atac_seurat,features=rownames(atac_seurat))
  adt_seurat <- ScaleData(adt_seurat,features=rownames(adt_seurat))
  
  rna.list <- SplitObject(rna_seurat, split.by = sample_meta)
  atac.list <- SplitObject(atac_seurat, split.by = sample_meta)
  adt.list <- SplitObject(adt_seurat, split.by = sample_meta)
  
  dataset_list <- list()
  for (i in 1:length(rna.list)){
    print(i)
    gene_seurat <- rna.list[[i]]
    gene_matrix <- gene_seurat@assays$originalexp$scale.data
    gene_matrix <- t(as.matrix(gene_matrix))
    
    atac_seurat <- atac.list[[i]]
    atac_matrix <- atac_seurat@assays$originalexp$scale.data
    atac_matrix <- t(as.matrix(atac_matrix))
    
    adt_seurat <- adt.list[[i]]
    adt_matrix <- adt_seurat@assays$originalexp$scale.data
    adt_matrix <- t(as.matrix(adt_matrix))
    
    dataset_list[[i]] <- list(gene_matrix, atac_matrix, adt_matrix)
  }
  
  eigen_output_cosine_list <- list()
  r_cosine_list <- c()
  block_matrix_list <- list()
  for (i in 1:length(dataset_list)){
    print(i)
    matrix_set <- dataset_list[[i]]
    
    if (block_method == "cosine"){
      print("cosine")
      eigen_output_cosine <- eigen_block_cosine_three(matrix_set[[1]],matrix_set[[2]],matrix_set[[3]])
    } else if (block_method == "gaussian"){
      print("gaussian")
      eigen_output_cosine <- eigen_block_gaussian_three(matrix_set[[1]],matrix_set[[2]],matrix_set[[3]])
    } else if (block_method == "emd"){
      print("emd")
      eigen_output_cosine <- eigen_block_emd_three(matrix_set[[1]],matrix_set[[2]],matrix_set[[3]])
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
  
  sample_name <- names(eigen_output_cosine_list)
  
  
  sample_condition_map <- unique(data.frame(
    sample_id = rna_seurat[[sample_meta]],
    condition = rna_seurat[[condition_meta]]
  ))
  
  
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
  
  plot_eigen(V_result$values)
  
  r2 <- find_elbow_kneedle(V_result$values[1:50])
  
  print(paste0("r2: ",r2))
  
  V_matrix <- V_result$vectors[,1:r2] %*% diag(V_result$values[1:r2])
  
  projected_combed_list <- list()
  for (i in 1:length(sample_projection_matrix)){
    projected_combed_list[[i]] <- block_matrix_list[[i]] %*% V_matrix
  }
  names(projected_combed_list) <- names(eigen_output_cosine_list)
  
  return(list(projected_combed_list=projected_combed_list,eigen_output_cosine_list=eigen_output_cosine_list,annotation=annotation))
}

## Pseudo-F (adonis2-consistent), vectorised over the columns of the group indicator A
sq_dist <- function(X) {
  cp <- tcrossprod(X); sq <- diag(cp)
  D2 <- outer(sq, sq, "+") - 2 * cp
  D2[D2 < 0] <- 0; diag(D2) <- 0
  D2
}
pseudoF_perms <- function(D2, A, Ac, nA, nB, S, G = 2) {
  SS_T <- sum(D2) / (2 * S)
  tA  <- colSums(A  * (D2 %*% A))
  tAc <- colSums(Ac * (D2 %*% Ac))
  SS_W <- tA / (2 * nA) + tAc / (2 * nB)
  SS_B <- SS_T - SS_W
  (SS_B / (G - 1)) / (SS_W / (S - G))
}
## F for every feature and every column of A
F_matrix <- function(emb_list, A, n_cores) {
  S <- nrow(A); Ac <- 1 - A; nA <- colSums(A); nB <- S - nA
  feat <- rownames(emb_list[[1]]); Fn <- length(feat); d <- ncol(emb_list[[1]])
  E <- array(NA_real_, dim = c(S, d, Fn))
  for (s in seq_len(S)) E[s, , ] <- t(emb_list[[s]])
  chunks <- split(seq_len(Fn), cut(seq_len(Fn), n_cores, labels = FALSE))
  out <- foreach(idx = chunks, .combine = rbind,
                 .export = c("sq_dist", "pseudoF_perms")) %dopar% {
    m <- matrix(NA_real_, length(idx), ncol(A))
    for (k in seq_along(idx))
      m[k, ] <- pseudoF_perms(sq_dist(E[, , idx[k]]), A, Ac, nA, nB, S)
    m
  }
  rownames(out) <- feat
  out
}

## One replicate
run_one <- function(sigma, rep, opt) {
  seed <- 1000L + rep
  t0 <- Sys.time()
  sim <- generate_svd_simulated_data(
    simulation_type = "DE", add_batch = FALSE, de_prop = 0, de_fc = 1,
    seed = seed, whether_rescale = FALSE, sample_noise_level = sigma,
    dia_matrix = diag(seq(5, 1, length.out = 10)))
  seu <- list(rna = sim$seurat_X, atac = sim$seurat_Y, adt = sim$seurat_Z)

  ## Observed embedding
  obs <- run_our_method_three(seu$rna, seu$atac, seu$adt, "sample_id", "condition",
                              block_method = "cosine")
  emb <- obs$projected_combed_list
  ann <- obs$annotation
  stopifnot(identical(names(emb), rownames(ann)))
  S <- length(emb); grp <- as.character(ann$Condition)
  stopifnot(S == 20L, all(table(grp) == 10L))
  feat <- rownames(emb[[1]])
  for (s in 2:S) stopifnot(identical(rownames(emb[[s]]), feat))

  ## Observed assignment, then B random assignments excluding it and its complement
  base <- as.integer(grp == grp[1])
  set.seed(seed * 7L + round(sigma * 1000))
  perms <- list(); keys <- character(0)
  comp_key <- paste(1 - base, collapse = "")
  while (length(perms) < opt$B) {
    a <- integer(S); a[sample(S, sum(base))] <- 1L
    k <- paste(a, collapse = "")
    if (k %in% c(paste(base, collapse = ""), comp_key, keys)) next
    perms[[length(perms) + 1]] <- a; keys <- c(keys, k)
  }
  A <- cbind(base, do.call(cbind, perms)); colnames(A) <- NULL

  Fmat  <- F_matrix(emb, A, opt$n_cores)          # features x (1 + B)
  F_obs <- Fmat[, 1]; F_perm <- Fmat[, -1, drop = FALSE]

  ## Check against adonis2 on a few features
  chk <- vapply(seq_len(min(5, length(feat))), function(f) {
    D <- as.dist(sqrt(sq_dist(t(sapply(emb, function(m) m[f, ])))))
    as.numeric(adonis2(D ~ grp, permutations = 9)$F[1])
  }, numeric(1))
  dev <- max(abs(chk - F_obs[seq_along(chk)]))
  if (dev > 1e-8) stop("pseudo-F mismatch vs adonis2: ", dev)

  ## Per-feature null
  p_new <- (1 + rowSums(F_perm >= F_obs)) / (opt$B + 1)

  ## Per-feature null heterogeneity
  q95 <- apply(F_perm, 1, quantile, probs = 0.95)
  het <- data.frame(sigma = sigma, rep = rep,
                    q95_median = median(q95), q95_IQR = IQR(q95),
                    q95_relIQR = IQR(q95) / median(q95))

  arms <- list(new = p_new)
  metrics <- do.call(rbind, lapply(names(arms), function(a) {
    p <- arms[[a]]
    data.frame(sigma = sigma, rep = rep, arm = a,
               alpha_hat = mean(p < opt$alpha),
               n_BH = sum(p.adjust(p, "BH") < opt$alpha),
               n_feat = length(p))
  }))
  message(sprintf("[sigma=%.2f rep=%d] alpha_hat %s | BH %s | %.1f min", sigma, rep,
                  paste(sprintf("%s=%.3f", metrics$arm, metrics$alpha_hat), collapse = " "),
                  paste(sprintf("%s=%d", metrics$arm, metrics$n_BH), collapse = " "),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  set.seed(seed)
  list(metrics = metrics, het = het, pvals = arms, F_obs = F_obs,
       F_perm_sub = sample(as.vector(F_perm), min(1e5, length(F_perm))))
}

## Main loop
res_file <- file.path(opt$out_dir, "null_calib_results.rds")
if (!opt$summary_only) {
  cl <- makeCluster(opt$n_cores); registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
  store <- if (file.exists(res_file)) readRDS(res_file) else list()
  for (rep in seq_len(opt$n_rep)) for (sigma in sigma_levels) {
    key <- sprintf("s%.3f_r%d", sigma, rep)
    if (!is.null(store[[key]])) { message("skip ", key); next }
    store[[key]] <- run_one(sigma, rep, opt)
    saveRDS(store, res_file)                          # checkpoint
  }
} else {
  store <- readRDS(res_file)
}

## Summaries
metrics <- do.call(rbind, lapply(store, `[[`, "metrics"))
het     <- do.call(rbind, lapply(store, `[[`, "het"))
write.csv(metrics, file.path(opt$out_dir, "typeI_by_replicate.csv"), row.names = FALSE)
write.csv(het,     file.path(opt$out_dir, "null_heterogeneity_by_replicate.csv"), row.names = FALSE)

ci <- function(x) { n <- length(x); if (n < 2) return(c(NA, NA))
  m <- mean(x); se <- sd(x) / sqrt(n); m + c(-1, 1) * qt(0.975, n - 1) * se }
summ_all <- metrics %>% group_by(sigma, arm) %>%
  summarise(n_rep = n(), alpha_mean = mean(alpha_hat),
            alpha_lo = ci(alpha_hat)[1], alpha_hi = ci(alpha_hat)[2],
            alpha_min = min(alpha_hat), alpha_max = max(alpha_hat),
            BH_mean = mean(n_BH), BH_max = max(n_BH), BH_any_prop = mean(n_BH > 0), .groups = "drop")
print(as.data.frame(summ_all), digits = 3)
write.csv(summ_all, file.path(opt$out_dir, "typeI_summary.csv"), row.names = FALSE)
