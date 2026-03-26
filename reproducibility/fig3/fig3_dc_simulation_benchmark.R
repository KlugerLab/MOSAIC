# Figure 3: DC detection benchmark on simulated data (MOSAIC vs DESeq2, edgeR, MAST, Seurat-Wilcoxon)
# Output:   results/fig3/benchmark_results.csv, figures/fig3A.png, fig3B.png, fig3C.png, fig3D.png

source("config.R")
source("R/functions.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(muscat)     # aggregateData(), pbDS() for DESeq2 / edgeR pseudobulk
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(MAST)       # used by Seurat::FindMarkers(test.use = "MAST")
  library(transport)  # loaded on PERMANOVA workers
  library(PRROC)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
})

out_dir <- file.path(RESULTS_DIR, "fig3")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

generate_svd_simulated_data <- function(
    simulation_type = "DC", # "DC" or "DE"
    add_batch = TRUE,
    de_prop = 0.2,
    de_fc = 2, # Fold change for DE simulation
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
  
  # Batches 1-2 contain only condition A samples, batches 3-4 only condition B
  sample_metadata <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    condition = rep(c("A", "B"), each = n_samples / n_conditions),
    batch = rep(c("batch_1", "batch_2", "batch_3", "batch_4"), each = n_samples / n_batch)
  )
  
  # ---- Global U (cell embeddings) ----
  # Shared cluster centers for cells across all samples
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

    # Cell embedding matrix U (not orthonormalized)
    U <- initial_points
    U_list[[j]] <- U
    cell_cluster_assignments_list[[j]] <- cluster_assignments
  }
  
  # ---- Global V (feature loadings) ----
  # Shared cluster centers for features across all modalities
  set.seed(seed)
  cluster_centers_V <- cluster_centers_U
  
  create_global_V <- function(p, n_clusters, r, cluster_centers) {
    cluster_assignments <- sample(1:n_clusters, p, replace = TRUE)
    initial_points <- matrix(0, nrow = p, ncol = r)
    for (i in 1:p) {
      cluster <- cluster_assignments[i]
      initial_points[i, ] <- cluster_centers[cluster, ] + rnorm(r, sd = 0.1)
    }
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
  
  # sample-specific V: global V plus noise
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
  
  # ---- Select DE / DC features ----
  set.seed(seed)
  de_features_x <- sample(1:px, round(px * de_prop))
  de_features_y <- sample(1:py, round(py * de_prop))
  de_features_z <- sample(1:pz, round(pz * de_prop))
  
  # Initialize final data lists
  X_list <- list()
  Y_list <- list()
  Z_list <- list()
  
  # Diagonal sigma matrices, identical for all modalities
  identical_sigma <- dia_matrix
  Sx <- identical_sigma
  Sy <- identical_sigma
  Sz <- identical_sigma
  
  # ---- DC or DE signal ----
  if (simulation_type == "DC") {
    # DC: permute loadings of selected features in condition B
    
    # each DC feature takes the loading row of a DC feature from another module (with replacement)
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
    
    # expression matrices
    for (i in 1:n_samples) {
      X_list[[i]] <- U_list[[i]] %*% Sx %*% t(Vx_list_de[[i]])
      Y_list[[i]] <- U_list[[i]] %*% Sy %*% t(Vy_list_de[[i]])
      Z_list[[i]] <- U_list[[i]] %*% Sz %*% t(Vz_list_de[[i]])
    }
    
    # optionally match condition B means of DC features to condition A
    
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
    # DE: multiplicative fold change in condition B
    
    for (i in 1:n_samples) {
      X_list[[i]] <- U_list[[i]] %*% Sx %*% t(Vx_list[[i]])
      Y_list[[i]] <- U_list[[i]] %*% Sy %*% t(Vy_list[[i]])
      Z_list[[i]] <- U_list[[i]] %*% Sz %*% t(Vz_list[[i]])
    }
    
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
  
  # ---- Condition-confounded affine perturbation (optional) ----
  if (add_batch) {
    for (i in 1:n_samples) {
      
      batch_num <- as.numeric(sub("batch_", "", sample_metadata$batch[i]))
      
      # modality X
      # per-feature scale and shift
      batch_scale_x_vec <- rnorm(px, batch_scale_mean, batch_scale_sd)
      batch_shift_x_vec <- rnorm(px, batch_shift_mean, batch_shift_sd) * (batch_num / n_batch)
      
      X_list[[i]] <- sweep(X_list[[i]], 2, batch_scale_x_vec, `*`) 
      X_list[[i]] <- sweep(X_list[[i]], 2, batch_shift_x_vec, `+`) 
      
      # modality Y
      batch_scale_y_vec <- rnorm(py, batch_scale_mean, batch_scale_sd)
      batch_shift_y_vec <- rnorm(py, batch_shift_mean, batch_shift_sd) * (batch_num / n_batch)
      
      Y_list[[i]] <- sweep(Y_list[[i]], 2, batch_scale_y_vec, `*`) 
      Y_list[[i]] <- sweep(Y_list[[i]], 2, batch_shift_y_vec, `+`) 
      
      # modality Z
      batch_scale_z_vec <- rnorm(pz, batch_scale_mean, batch_scale_sd)
      batch_shift_z_vec <- rnorm(pz, batch_shift_mean, batch_shift_sd) * (batch_num / n_batch)
      
      Z_list[[i]] <- sweep(Z_list[[i]], 2, batch_scale_z_vec, `*`) 
      Z_list[[i]] <- sweep(Z_list[[i]], 2, batch_shift_z_vec, `+`) 
    }
  }
  
  
  # ---- Seurat objects ----
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
    matrix_set <- dataset_list[[i]]
    
    if (block_method == "cosine"){
      eigen_output_cosine <- eigen_block_cosine_three(matrix_set[[1]],matrix_set[[2]],matrix_set[[3]])
    } else if (block_method == "gaussian"){
      eigen_output_cosine <- eigen_block_gaussian_three(matrix_set[[1]],matrix_set[[2]],matrix_set[[3]])
    } else if (block_method == "emd"){
      eigen_output_cosine <- eigen_block_emd_three(matrix_set[[1]],matrix_set[[2]],matrix_set[[3]])
    } else {
      stop("block method not found")
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

run_deseq2 <- function(seurat_rna, seurat_atac, seurat_adt, labels, feature_names) {
  
  # DESeq2 on one modality
  run_deseq2_single <- function(seurat_obj, labels_subset) {
    seurat_obj$group_id <- as.factor(seurat_obj$condition)
    sce <- as.SingleCellExperiment(seurat_obj)
    sce$cluster_id <- "cluster1"
    pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
    metadata(pb)$experiment_info <- 1
    pb@assays@data$cluster1 <- pb@assays@data$cluster1 + abs(min(pb@assays@data$cluster1)) + 1
    ds_pb <- pbDS(pb, method = "DESeq2")
    
    scores <- -log10(ds_pb[["table"]]$B$cluster1$p_adj.glb)
    if (any(is.infinite(scores))) {
      max_score <- max(scores[is.finite(scores)], na.rm = TRUE)
      scores[is.infinite(scores)] <- max_score + 1
    }
    
    pr_obj <- pr.curve(scores.class0 = scores, weights.class0 = labels_subset, curve = TRUE)
    return(pr_obj)
  }
  
  # Split labels
  n_rna <- nrow(seurat_rna)
  n_atac <- nrow(seurat_atac)
  n_adt <- nrow(seurat_adt)
  
  labels_1 <- labels[1:n_rna]
  labels_2 <- labels[(n_rna + 1):(n_rna + n_atac)]
  labels_3 <- labels[(n_rna + n_atac + 1):(n_rna + n_atac + n_adt)]
  
  pr_obj_1 <- run_deseq2_single(seurat_rna, labels_1)
  pr_obj_2 <- run_deseq2_single(seurat_atac, labels_2)
  pr_obj_3 <- run_deseq2_single(seurat_adt, labels_3)
  
  return(list(
    pr_obj_1 = pr_obj_1,
    pr_obj_2 = pr_obj_2,
    pr_obj_3 = pr_obj_3,
    mean_auprc = mean(c(pr_obj_1$auc.integral, pr_obj_2$auc.integral, pr_obj_3$auc.integral))
  ))
}

run_edger <- function(seurat_rna, seurat_atac, seurat_adt, labels, feature_names) {
  
  run_edger_single <- function(seurat_obj, labels_subset) {
    seurat_obj$group_id <- as.factor(seurat_obj$condition)
    sce <- as.SingleCellExperiment(seurat_obj)
    sce$cluster_id <- "cluster1"
    pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
    metadata(pb)$experiment_info <- 1
    pb@assays@data$cluster1 <- pb@assays@data$cluster1 + abs(min(pb@assays@data$cluster1)) + 1
    ds_pb <- pbDS(pb, method = "edgeR")
    
    scores <- -log10(ds_pb[["table"]]$B$cluster1$p_adj.glb)
    if (any(is.infinite(scores))) {
      max_score <- max(scores[is.finite(scores)], na.rm = TRUE)
      scores[is.infinite(scores)] <- max_score + 1
    }
    
    pr_obj <- pr.curve(scores.class0 = scores, weights.class0 = labels_subset, curve = TRUE)
    return(pr_obj)
  }
  
  n_rna <- nrow(seurat_rna)
  n_atac <- nrow(seurat_atac)
  n_adt <- nrow(seurat_adt)
  
  labels_1 <- labels[1:n_rna]
  labels_2 <- labels[(n_rna + 1):(n_rna + n_atac)]
  labels_3 <- labels[(n_rna + n_atac + 1):(n_rna + n_atac + n_adt)]
  
  pr_obj_1 <- run_edger_single(seurat_rna, labels_1)
  pr_obj_2 <- run_edger_single(seurat_atac, labels_2)
  pr_obj_3 <- run_edger_single(seurat_adt, labels_3)
  
  return(list(
    pr_obj_1 = pr_obj_1,
    pr_obj_2 = pr_obj_2,
    pr_obj_3 = pr_obj_3,
    mean_auprc = mean(c(pr_obj_1$auc.integral, pr_obj_2$auc.integral, pr_obj_3$auc.integral))
  ))
}

run_mast <- function(seurat_rna, seurat_atac, seurat_adt, labels, feature_names) {
  
  run_mast_single <- function(seurat_obj, labels_subset, feature_names_subset) {
    de_result <- FindMarkers(seurat_obj, slot = "data", group.by = "group_id", 
                             ident.1 = "A", ident.2 = "B", test.use = "MAST", logfc.threshold = 0)
    de_result$gene <- rownames(de_result)
    
    # -log10 adjusted P per feature, 0 if not returned
    scores <- numeric(length(feature_names_subset))
    names(scores) <- feature_names_subset
    
    for (gene in feature_names_subset) {
      if (gene %in% de_result$gene) {
        scores[gene] <- -log10(de_result[gene, "p_val_adj"])
      } else {
        scores[gene] <- 0
      }
    }
    
    # Handle infinite values
    if (any(is.infinite(scores))) {
      max_score <- max(scores[is.finite(scores)], na.rm = TRUE)
      scores[is.infinite(scores)] <- max_score + 1
    }
    
    pr_obj <- pr.curve(scores.class0 = scores, weights.class0 = labels_subset, curve = TRUE)
    return(pr_obj)
  }
  
  n_rna <- nrow(seurat_rna)
  n_atac <- nrow(seurat_atac)
  n_adt <- nrow(seurat_adt)
  
  labels_1 <- labels[1:n_rna]
  labels_2 <- labels[(n_rna + 1):(n_rna + n_atac)]
  labels_3 <- labels[(n_rna + n_atac + 1):(n_rna + n_atac + n_adt)]
  
  seurat_rna$group_id <- as.factor(seurat_rna$condition)
  seurat_atac$group_id <- as.factor(seurat_atac$condition)
  seurat_adt$group_id <- as.factor(seurat_adt$condition)
  
  pr_obj_1 <- run_mast_single(seurat_rna, labels_1, feature_names[1:n_rna])
  pr_obj_2 <- run_mast_single(seurat_atac, labels_2, feature_names[(n_rna + 1):(n_rna + n_atac)])
  pr_obj_3 <- run_mast_single(seurat_adt, labels_3, feature_names[(n_rna + n_atac + 1):length(feature_names)])
  
  return(list(
    pr_obj_1 = pr_obj_1,
    pr_obj_2 = pr_obj_2,
    pr_obj_3 = pr_obj_3,
    mean_auprc = mean(c(pr_obj_1$auc.integral, pr_obj_2$auc.integral, pr_obj_3$auc.integral))
  ))
}

run_wilcox <- function(seurat_rna, seurat_atac, seurat_adt, labels, feature_names) {
  
  run_wilcox_single <- function(seurat_obj, labels_subset, feature_names_subset) {
    de_result <- FindMarkers(seurat_obj, slot = "data", group.by = "group_id",
                             ident.1 = "A", ident.2 = "B", test.use = "wilcox", logfc.threshold = 0)
    de_result$gene <- rownames(de_result)
    
    # -log10 adjusted P per feature, 0 if not returned
    scores <- numeric(length(feature_names_subset))
    names(scores) <- feature_names_subset
    
    for (gene in feature_names_subset) {
      if (gene %in% de_result$gene) {
        scores[gene] <- -log10(de_result[gene, "p_val_adj"])
      } else {
        scores[gene] <- 0
      }
    }
    
    # Handle infinite values
    if (any(is.infinite(scores))) {
      max_score <- max(scores[is.finite(scores)], na.rm = TRUE)
      scores[is.infinite(scores)] <- max_score + 1
    }
    
    pr_obj <- pr.curve(scores.class0 = scores, weights.class0 = labels_subset, curve = TRUE)
    return(pr_obj)
  }
  
  n_rna <- nrow(seurat_rna)
  n_atac <- nrow(seurat_atac)
  n_adt <- nrow(seurat_adt)
  
  labels_1 <- labels[1:n_rna]
  labels_2 <- labels[(n_rna + 1):(n_rna + n_atac)]
  labels_3 <- labels[(n_rna + n_atac + 1):(n_rna + n_atac + n_adt)]
  
  seurat_rna$group_id <- as.factor(seurat_rna$condition)
  seurat_atac$group_id <- as.factor(seurat_atac$condition)
  seurat_adt$group_id <- as.factor(seurat_adt$condition)
  
  pr_obj_1 <- run_wilcox_single(seurat_rna, labels_1, feature_names[1:n_rna])
  pr_obj_2 <- run_wilcox_single(seurat_atac, labels_2, feature_names[(n_rna + 1):(n_rna + n_atac)])
  pr_obj_3 <- run_wilcox_single(seurat_adt, labels_3, feature_names[(n_rna + n_atac + 1):length(feature_names)])
  
  return(list(
    pr_obj_1 = pr_obj_1,
    pr_obj_2 = pr_obj_2,
    pr_obj_3 = pr_obj_3,
    mean_auprc = mean(c(pr_obj_1$auc.integral, pr_obj_2$auc.integral, pr_obj_3$auc.integral))
  ))
}

run_benchmark <- function(
    simulation_type = "DC",
    add_batch = TRUE,
    de_prop = 0.2,
    de_fc = 2,
    seed = 1,
    whether_rescale = TRUE,
    methods = c("custom", "DESeq2", "edgeR", "MAST", "wilcox"),
    verbose = TRUE,
    dia_type = "medium"
) {
  
  if (dia_type == "medium"){
    dia_matrix = diag(seq(5, 1, length.out = 10))
  } else if (dia_type == "strong"){
    dia_matrix = diag(seq(10, 5, length.out = 10))
  } else if (dia_type == "weak"){
    dia_matrix = diag(seq(3, 0.5, length.out = 10))
  }
  
  if (verbose) {
    cat(sprintf("\n=== Running Simulation ===\n"))
    cat(sprintf("Type: %s |Method: %s | Batch: %s | DE prop: %.2f | Seed: %d\n |diag_type: %s", 
                simulation_type, methods, add_batch, de_prop, seed, dia_type))
  }
  
  # Generate simulated data
  sim_data <- generate_svd_simulated_data(
    simulation_type = simulation_type,
    add_batch = add_batch,
    de_prop = de_prop,
    de_fc = de_fc,
    seed = seed,
    whether_rescale = whether_rescale,
    dia_matrix = dia_matrix
  )
  
  seurat_rna <- sim_data$seurat_X
  seurat_atac <- sim_data$seurat_Y
  seurat_adt <- sim_data$seurat_Z
  
  # Extract ground truth labels
  feature_names <- c(rownames(seurat_rna), rownames(seurat_atac), rownames(seurat_adt))
  
  rna_de <- rownames(seurat_rna[["originalexp"]][[]])[which(seurat_rna[["originalexp"]][["is_de"]] == TRUE)]
  atac_de <- rownames(seurat_atac[["originalexp"]][[]])[which(seurat_atac[["originalexp"]][["is_de"]] == TRUE)]
  adt_de <- rownames(seurat_adt[["originalexp"]][[]])[which(seurat_adt[["originalexp"]][["is_de"]] == TRUE)]
  
  true_differential_idx <- which(feature_names %in% c(rna_de, atac_de, adt_de))
  labels <- as.integer(1:length(feature_names) %in% true_differential_idx)
  
  # Initialize results list
  results <- list(
    simulation_params = list(
      simulation_type = simulation_type,
      add_batch = add_batch,
      de_prop = de_prop,
      de_fc = de_fc,
      seed = seed,
      whether_rescale = whether_rescale
    ),
    auprc_results = list()
  )
  
  # MOSAIC (custom_cosine); custom_emd is the EMD-kernel variant, not plotted
  if ("custom" %in% methods) {
    block_methods <- c("cosine", "emd")
    
    for (block_method in block_methods) {
      method_name <- paste0("custom_", block_method)
      
      if (verbose) cat(sprintf("Running custom method (%s)...\n", block_method))
      
      tryCatch({
        custom_result <- run_our_method_three(seurat_rna, seurat_atac, seurat_adt, 
                                           "sample_id", "condition", block_method = block_method)
        
        n_sample <- length(custom_result$projected_combed_list)
        peranova_result <- cal_peranova_per_feature_parallel_no_zero(
          custom_result$projected_combed_list, n_sample,
          groups = custom_result$annotation[names(custom_result$eigen_output_cosine_list), "Condition"],
          dist_method = "eu", num_permutation = 999, n_cores = N_CORES
        )
        
        # rank features by observed PERMANOVA pseudo-F
        scores <- unlist(peranova_result$F_stats)
        if (any(is.infinite(scores))) {
          max_score <- max(scores[is.finite(scores)], na.rm = TRUE)
          scores[is.infinite(scores)] <- max_score + 1
        }
        
        # AUPRC over all features pooled across modalities
        pr_obj <- pr.curve(scores.class0 = scores, weights.class0 = labels, curve = TRUE)
        results$auprc_results[[method_name]] <- pr_obj$auc.integral
        
        if (verbose) cat(sprintf("  Custom (%s) AUPRC: %.4f\n", block_method, pr_obj$auc.integral))
        
      }, error = function(e) {
        if (verbose) cat(sprintf("  Custom method (%s) failed: %s\n", block_method, e$message))
        results$auprc_results[[method_name]] <- NA
      })
    }
  }
  

  # Run DESeq2
  if ("DESeq2" %in% methods) {
    if (verbose) cat("Running DESeq2...\n")
    tryCatch({
      deseq2_result <- run_deseq2(seurat_rna, seurat_atac, seurat_adt, labels, feature_names)
      results$auprc_results$DESeq2 <- deseq2_result$mean_auprc
      results$auprc_results$DESeq2_modality1 <- deseq2_result$pr_obj_1$auc.integral
      results$auprc_results$DESeq2_modality2 <- deseq2_result$pr_obj_2$auc.integral
      results$auprc_results$DESeq2_modality3 <- deseq2_result$pr_obj_3$auc.integral
      
      if (verbose) cat(sprintf("  DESeq2 mean AUPRC: %.4f\n", deseq2_result$mean_auprc))
    }, error = function(e) {
      if (verbose) cat(sprintf("  DESeq2 failed: %s\n", e$message))
      results$auprc_results$DESeq2 <- NA
    })
  }
  
  # Run edgeR
  if ("edgeR" %in% methods) {
    if (verbose) cat("Running edgeR...\n")
    tryCatch({
      edger_result <- run_edger(seurat_rna, seurat_atac, seurat_adt, labels, feature_names)
      results$auprc_results$edgeR <- edger_result$mean_auprc
      results$auprc_results$edgeR_modality1 <- edger_result$pr_obj_1$auc.integral
      results$auprc_results$edgeR_modality2 <- edger_result$pr_obj_2$auc.integral
      results$auprc_results$edgeR_modality3 <- edger_result$pr_obj_3$auc.integral
      
      if (verbose) cat(sprintf("  edgeR mean AUPRC: %.4f\n", edger_result$mean_auprc))
    }, error = function(e) {
      if (verbose) cat(sprintf("  edgeR failed: %s\n", e$message))
      results$auprc_results$edgeR <- NA
    })
  }
  
  # Run MAST
  if ("MAST" %in% methods) {
    if (verbose) cat("Running MAST...\n")
    tryCatch({
      mast_result <- run_mast(seurat_rna, seurat_atac, seurat_adt, labels, feature_names)
      results$auprc_results$MAST <- mast_result$mean_auprc
      results$auprc_results$MAST_modality1 <- mast_result$pr_obj_1$auc.integral
      results$auprc_results$MAST_modality2 <- mast_result$pr_obj_2$auc.integral
      results$auprc_results$MAST_modality3 <- mast_result$pr_obj_3$auc.integral
      
      if (verbose) cat(sprintf("  MAST mean AUPRC: %.4f\n", mast_result$mean_auprc))
    }, error = function(e) {
      if (verbose) cat(sprintf("  MAST failed: %s\n", e$message))
      results$auprc_results$MAST <- NA
    })
  }
  
  # Run Wilcoxon
  if ("wilcox" %in% methods) {
    if (verbose) cat("Running Wilcoxon...\n")
    tryCatch({
      wilcox_result <- run_wilcox(seurat_rna, seurat_atac, seurat_adt, labels, feature_names)
      results$auprc_results$wilcox <- wilcox_result$mean_auprc
      results$auprc_results$wilcox_modality1 <- wilcox_result$pr_obj_1$auc.integral
      results$auprc_results$wilcox_modality2 <- wilcox_result$pr_obj_2$auc.integral
      results$auprc_results$wilcox_modality3 <- wilcox_result$pr_obj_3$auc.integral
      
      if (verbose) cat(sprintf("  Wilcoxon mean AUPRC: %.4f\n", wilcox_result$mean_auprc))
    }, error = function(e) {
      if (verbose) cat(sprintf("  Wilcoxon failed: %s\n", e$message))
      results$auprc_results$wilcox <- NA
    })
  }
  
  return(results)
}

run_multiple_benchmarks <- function(param_grid, methods = c("custom", "DESeq2", "edgeR", "MAST", "wilcox"), 
                                    verbose = TRUE) {
  
  results_list <- list()
  
  for (i in 1:nrow(param_grid)) {
    if (verbose) {
      cat(sprintf("\n\n========================================\n"))
      cat(sprintf("Running benchmark %d of %d\n", i, nrow(param_grid)))
      cat(sprintf("========================================\n"))
    }
    
    # Extract parameters for this run
    params <- as.list(param_grid[i, ])
    
    # Run benchmark
    result <- do.call(run_benchmark, c(params, list(methods = methods, verbose = verbose)))
    
    # Convert to data frame row
    result_row <- data.frame(
      run_id = i,
      simulation_type = params$simulation_type,
      add_batch = params$add_batch,
      dia_type = params$dia_type,
      de_prop = params$de_prop,
      de_fc = ifelse(is.null(params$de_fc), 2, params$de_fc),
      seed = params$seed,
      whether_rescale = ifelse(is.null(params$whether_rescale), TRUE, params$whether_rescale)
    )
    
    # Add AUPRC results
    for (method_name in names(result$auprc_results)) {
      result_row[[paste0("auprc_", method_name)]] <- result$auprc_results[[method_name]]
    }
    
    results_list[[i]] <- result_row
  }
  
  # Combine all results
  results_df <- do.call(rbind, results_list)
  
  return(results_df)
}

# ---- Run grid ----
# 12 DC runs: perturbation (add_batch) x mean correction (whether_rescale) x signal strength
param_grid_dc <- expand.grid(
  simulation_type = c("DC"),
  add_batch = c(TRUE, FALSE),
  whether_rescale = c(TRUE,FALSE),
  de_prop = c(0.2),
  de_fc = 2,
  seed = 190,
  dia_type = c("strong","medium","weak"),
  stringsAsFactors = FALSE
)

# mean shift only (DE) run for Fig 3A
param_grid_de <- data.frame(
  simulation_type = "DE",
  add_batch = FALSE,
  whether_rescale = TRUE,
  de_prop = 0.2,
  de_fc = 2,
  seed = 190,
  dia_type = "medium",
  stringsAsFactors = FALSE
)

param_grid <- rbind(param_grid_dc, param_grid_de)   # 13 runs; DE is run 13

results_df <- run_multiple_benchmarks(
  param_grid = param_grid,
  methods = c("custom", "DESeq2", "edgeR", "MAST", "wilcox"),
  verbose = TRUE
)

# Save results
write.csv(results_df, file.path(out_dir, "benchmark_results.csv"), row.names = FALSE)

# ---- Figure 3A: AUPRC by signal type ----
# no perturbation, medium signal
filtered_results <- results_df %>%
  filter(add_batch == FALSE, dia_type == "medium")

# Reshape data for plotting
plot_data <- filtered_results %>%
  dplyr::select(simulation_type,
                whether_rescale,
                custom_cosine = auprc_custom_cosine,
                custom_emd = auprc_custom_emd,
                DESeq2 = auprc_DESeq2,
                edgeR = auprc_edgeR,
                MAST = auprc_MAST,
                wilcox = auprc_wilcox) %>%
  mutate(rescale_label = ifelse(simulation_type == "DE", "Mean Shift Only",
                                ifelse(whether_rescale, "Pure Connectivity Signal (Rewiring Only)", "Confounded Signal (Rewiring + Mean Shift)")))

# Convert to long format using reshape2::melt
plot_data_long <- melt(plot_data, 
                       id.vars = c("simulation_type", "whether_rescale", "rescale_label"),
                       variable.name = "method",
                       value.name = "auprc")

plot_data_long$rescale_label <- factor(plot_data_long$rescale_label, 
                                       levels = c("Confounded Signal (Rewiring + Mean Shift)", "Pure Connectivity Signal (Rewiring Only)",
                                                  "Mean Shift Only"))

plot_data_long <- plot_data_long[which(!plot_data_long$method == "custom_emd"),]

plot_data_long$method <- fct_recode(plot_data_long$method, 
                                    "MOSAIC" = "custom_cosine")
plot_data_long$method <- fct_recode(plot_data_long$method, 
                                    "Seurat-Wilcoxon" = "wilcox")

# Create grouped bar plot
p <- ggplot(plot_data_long, aes(x = method, y = auprc, fill = rescale_label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(auprc, 3)), 
            position = position_dodge(width = 0.8),
            vjust = -0.5,
            size = 3) +
  scale_fill_manual(values = c("Confounded Signal (Rewiring + Mean Shift)" = "#f1a226",
                               "Pure Connectivity Signal (Rewiring Only)" = "#298c8c",
                               "Mean Shift Only" = "#7f8c8d"),
                    labels = c("Confounded Signal\n(Rewiring + Mean Shift)",
                               "Pure Connectivity Signal\n(Rewiring Only)",
                               "Mean Shift Only")) +
  labs(
    title = "",
    subtitle = "",
    x = "Method",
    y = "AUPRC",
    fill = "Signal Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, NA))

png(file.path(FIG_DIR, "fig3A.png"),width = 8.5, height = 3.5, units = "in", res = 600)
print(p)
dev.off()

# ---- Figure 3B: with / without perturbation ----
# pure connectivity, medium signal
filtered_results <- results_df %>%
  filter(simulation_type == "DC", whether_rescale == TRUE, dia_type == "medium")

# Reshape data for plotting
plot_data <- filtered_results %>%
  dplyr::select(add_batch,
                custom_cosine = auprc_custom_cosine,
                custom_emd = auprc_custom_emd,
                DESeq2 = auprc_DESeq2,
                edgeR = auprc_edgeR,
                MAST = auprc_MAST,
                wilcox = auprc_wilcox) %>%
  mutate(batch_label = ifelse(add_batch, "With perturbation", "Without perturbation"))

# Convert to long format
plot_data_long <- melt(plot_data, 
                       id.vars = c("add_batch", "batch_label"),
                       variable.name = "method",
                       value.name = "auprc")

plot_data_long$batch_label <- factor(plot_data_long$batch_label, 
                                     levels = c("Without perturbation", "With perturbation"))

plot_data_long <- plot_data_long[which(!plot_data_long$method == "custom_emd"),]

plot_data_long$method <- fct_recode(plot_data_long$method, 
                                    "MOSAIC" = "custom_cosine")
plot_data_long$method <- fct_recode(plot_data_long$method, 
                                    "Seurat-Wilcoxon" = "wilcox")

# Create slope data
slope_data <- dcast(plot_data_long, method ~ batch_label, value.var = "auprc")
names(slope_data) <- c("method", "Without", "With")

# Create slope graph
p <- ggplot(slope_data, aes(x = 1, xend = 2, y = Without, yend = With)) +
  geom_segment(aes(color = method), size = 1.5, alpha = 0.7) +
  geom_point(aes(x = 1, y = Without, color = method), size = 4) +
  geom_point(aes(x = 2, y = With, color = method), size = 4) +
  geom_text_repel(aes(x = 1, y = Without, label = method, color = method), 
                  hjust = 1, nudge_x = -0.1, direction = "y", 
                  segment.size = 0.2, size = 4, fontface = "bold",
                  box.padding = 0.5, point.padding = 0.5) +
  scale_x_continuous(limits = c(0.2, 2.7), breaks = c(1, 2), 
                     labels = c("Without\nperturbation", "With\nperturbation")) +
  scale_y_continuous(limits = c(0, max(c(slope_data$Without, slope_data$With)) * 1.05)) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "",
    subtitle = "",
    y = "AUPRC",
    x = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "gray30", hjust = 0.5),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    axis.text.x = element_text(size = 12, face = "bold")
  )

png(file.path(FIG_DIR, "fig3B.png"),width = 4, height = 4, units = "in", res = 600)
print(p)
dev.off()

# ---- Figure 3C / 3D: signal strength sweep ----
# pure connectivity, without (3C) and with (3D) perturbation
signal_strength_files <- c(
  "FALSE" = "fig3C.png",    # Fig 3C
  "TRUE"  = "fig3D.png"   # Fig 3D
)

for (batch_setting in c(FALSE, TRUE)) {
  
  filtered_results <- results_df %>%
    filter(simulation_type == "DC", whether_rescale == TRUE, add_batch == batch_setting)
  
  # Reshape data for plotting
  plot_data <- filtered_results %>%
    dplyr::select(dia_type,
                  custom_cosine = auprc_custom_cosine,
                  custom_emd = auprc_custom_emd,
                  DESeq2 = auprc_DESeq2,
                  edgeR = auprc_edgeR,
                  MAST = auprc_MAST,
                  wilcox = auprc_wilcox)
  
  # Convert to long format
  plot_data_long <- melt(plot_data, 
                         id.vars = c("dia_type"),
                         variable.name = "method",
                         value.name = "auprc")
  
  # Remove custom_emd
  plot_data_long <- plot_data_long[which(!plot_data_long$method == "custom_emd"),]
  
  # Rename methods
  plot_data_long$method <- fct_recode(plot_data_long$method, 
                                      "MOSAIC" = "custom_cosine")
  plot_data_long$method <- fct_recode(plot_data_long$method, 
                                      "Seurat-Wilcoxon" = "wilcox")
  
  # Set factor levels for dia_type
  plot_data_long$dia_type <- factor(plot_data_long$dia_type, 
                                    levels = c("strong", "medium", "weak"))
  
  # Create wide format for slope graph
  slope_data <- dcast(plot_data_long, method ~ dia_type, value.var = "auprc")
  names(slope_data) <- c("method", "Strong", "Medium", "Weak")
  
  # Create slope graph
  p <- ggplot(slope_data) +
    # Lines from Strong to Medium
    geom_segment(aes(x = 1, xend = 2, y = Strong, yend = Medium, color = method), 
                 size = 1.5, alpha = 0.7) +
    # Lines from Medium to Weak
    geom_segment(aes(x = 2, xend = 3, y = Medium, yend = Weak, color = method), 
                 size = 1.5, alpha = 0.7) +
    # Points
    geom_point(aes(x = 1, y = Strong, color = method), size = 4) +
    geom_point(aes(x = 2, y = Medium, color = method), size = 4) +
    geom_point(aes(x = 3, y = Weak, color = method), size = 4) +
    # Labels on the left (Strong)
    geom_text_repel(aes(x = 1, y = Strong, label = method, color = method), 
                    hjust = 1, nudge_x = -0.1, direction = "y", 
                    segment.size = 0.2, size = 4, fontface = "bold",
                    box.padding = 0.5, point.padding = 0.5) +
    scale_x_continuous(limits = c(0.2, 3.8), breaks = c(1, 2, 3), 
                       labels = c("Strong", "Medium", "Weak")) +
    scale_y_continuous(limits = c(0, max(c(slope_data$Strong, slope_data$Medium, slope_data$Weak)) * 1.05)) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      title = "",
      subtitle = "",
      y = "AUPRC",
      x = "Signal Strength"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11, color = "gray30", hjust = 0.5),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray90"),
      axis.text.x = element_text(size = 12, face = "bold")
    )
  
  png(file.path(FIG_DIR, signal_strength_files[[as.character(batch_setting)]]),
      width = 4, height = 4, units = "in", res = 600)
  print(p)
  dev.off()
}
