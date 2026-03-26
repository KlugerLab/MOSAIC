#' Simulate Multi-Modal Single-Cell Data
#'
#' Creates synthetic multi-modal single-cell data following the generative model
#' X = U * S * t(V), where U is a sample-specific cell embedding matrix, S is a
#' diagonal singular value matrix, and V is a sample-specific feature loading
#' matrix. The number of modalities (1 to 3) is determined by the length of
#' \code{n_features}. The simulation generates 20 samples in 2 conditions
#' (10 per condition), with 10 latent cell types.
#'
#' Two simulation modes are supported:
#' \itemize{
#'   \item DC (Differential Connectivity): Feature loadings in condition B are
#'     permuted across clusters for a fraction of features, rewiring their
#'     inter-feature relationships while optionally preserving mean expression.
#'   \item DE (Differential Expression): A multiplicative fold change is applied
#'     to a fraction of features in condition B, altering abundance without
#'     changing connectivity.
#' }
#'
#' The returned Seurat objects include per-feature metadata (\code{feature_cluster}
#' and \code{is_de}) accessible via
#' \code{seurat_obj[["originalexp"]][[]]$feature_cluster}, which provides ground
#' truth for benchmarking.
#'
#' @param simulation_type Character. \code{"DC"} for differential connectivity
#'   or \code{"DE"} for differential expression.
#' @param n_features Integer vector specifying the number of features per
#'   modality. Its length determines the number of modalities (1 to 3).
#'   Default: \code{c(1000, 800, 600)} (3 modalities).
#' @param signal_prop Numeric between 0 and 1. Proportion of features affected
#'   by DC or DE signal in each modality. Default: \code{0.2}.
#' @param fold_change Numeric. Fold change applied to affected features in
#'   condition B. Only used when \code{simulation_type = "DE"}.
#'   Default: \code{2}.
#' @param seed Integer. Random seed for reproducibility.
#' @param rescale Logical. For DC simulation only: if \code{TRUE} (default),
#'   rescale mean expression of DC features in condition B to match condition A,
#'   ensuring a pure connectivity signal without abundance change.
#' @param sample_noise_level Numeric. Standard deviation of Gaussian noise added
#'   to each sample's feature loading matrix. Controls between-sample
#'   variability. Default: \code{0.5}.
#' @param feature_noise_level Numeric. Standard deviation of Gaussian noise
#'   around feature cluster centers in the loading space. Controls how tight
#'   feature clusters are. Default: \code{0.1}.
#' @param cell_noise_level Numeric. Standard deviation of Gaussian noise around
#'   cell type cluster centers in the cell embedding space. Controls cell type
#'   separation. Default: \code{0.1}.
#' @param nonlinear Logical. If \code{TRUE}, apply a sigmoid transformation
#'   (1 / (1 + exp(-x))) to the expression matrices after the linear generative
#'   step. Simulates nonlinear gene regulation. Default: \code{FALSE}.
#' @param singular_values Numeric vector of length \code{r} (latent rank,
#'   default 10) specifying the diagonal entries of the singular value matrix.
#'   If \code{NULL} (default), uses \code{seq(5, 1, length.out = 10)}.
#' @param add_batch Logical. If \code{TRUE} (default), add per-feature
#'   batch-specific scale and shift effects across 4 batches.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{seurat_list}}{A named list of Seurat objects, one per modality
#'     (named \code{"modality_1"}, \code{"modality_2"}, etc.). Each uses assay
#'     \code{"originalexp"}. Feature metadata includes \code{feature_cluster}
#'     (integer 1-10) and \code{is_de} (logical). Cell metadata includes
#'     \code{sample_id}, \code{condition} ("A" or "B"), \code{batch}, and
#'     \code{true_cell_type}.}
#'   \item{\code{raw_matrix_list}}{A named list (one per modality) of lists
#'     containing 20 raw expression matrices (cells x features), one per
#'     sample.}
#'   \item{\code{sample_metadata}}{Data frame with columns \code{sample_id},
#'     \code{condition}, and \code{batch}.}
#' }
#'
#' @seealso \code{\link{run_MOSAIC}} to run the MOSAIC pipeline on the
#'   simulated data.
#'
#' @examples
#' \dontrun{
#' # --- Simulate 2-modality DC data with nonlinearity ---
#' sim <- simulate_multimodal_data(
#'   simulation_type = "DC",
#'   n_features = c(1000, 800),
#'   signal_prop = 0.2,
#'   seed = 42,
#'   nonlinear = TRUE
#' )
#'
#' # Ground truth
#' is_dc <- sim$seurat_list[[1]][["originalexp"]][[]]$is_de
#'
#' # Run MOSAIC
#' result <- run_MOSAIC(
#'   sim$seurat_list,
#'   assays = rep("originalexp", 2),
#'   sample_meta = "sample_id",
#'   condition_meta = "condition"
#' )
#'
#' # --- Simulate 3-modality DE data ---
#' sim_de <- simulate_multimodal_data(
#'   simulation_type = "DE",
#'   n_features = c(1000, 800, 600),
#'   signal_prop = 0.3,
#'   fold_change = 3,
#'   seed = 123
#' )
#'
#' # --- Simulate 1-modality with custom singular values ---
#' sim_1mod <- simulate_multimodal_data(
#'   n_features = c(500),
#'   singular_values = seq(10, 1, length.out = 10),
#'   seed = 1
#' )
#' }
#'
#' @export
simulate_multimodal_data <- function(
    simulation_type = "DC",
    n_features = c(1000, 800, 600),
    signal_prop = 0.2,
    fold_change = 2,
    seed = 1,
    rescale = TRUE,
    sample_noise_level = 0.5,
    feature_noise_level = 0.1,
    cell_noise_level = 0.1,
    nonlinear = FALSE,
    singular_values = NULL,
    add_batch = TRUE
) {

  n_modalities <- length(n_features)
  if (n_modalities < 1 || n_modalities > 3) {
    stop("n_features must have length 1, 2, or 3.")
  }

  set.seed(seed)

  # Basic parameters
  n_samples <- 20
  n_conditions <- 2
  n_batch <- 4
  cells_per_sample <- 1000
  r <- 10
  batch_shift_mean <- 10
  batch_shift_sd <- 0.2
  batch_scale_mean <- 1.2
  batch_scale_sd <- 0.1
  n_clusters <- 10

  if (is.null(singular_values)) {
    singular_values <- seq(5, 1, length.out = r)
  }
  sigma_matrix <- diag(singular_values)

  sample_metadata <- data.frame(
    sample_id = paste0("sample_", seq_len(n_samples)),
    condition = rep(c("A", "B"), each = n_samples / n_conditions),
    batch = rep(paste0("batch_", seq_len(n_batch)), each = n_samples / n_batch),
    stringsAsFactors = FALSE
  )

  apply_non_linearity <- function(mat) {
    1 / (1 + exp(-mat))
  }

  # ---- Global U (Cell Embeddings) ----
  cluster_centers_U <- matrix(stats::rnorm(n_clusters * r), nrow = n_clusters, ncol = r)

  U_list <- list()
  cell_cluster_assignments_list <- list()
  for (j in seq_len(n_samples)) {
    n_cells <- cells_per_sample + sample(-200:200, 1)
    cluster_assignments <- sample(seq_len(n_clusters), n_cells, replace = TRUE)
    initial_points <- matrix(0, nrow = n_cells, ncol = r)
    for (i in seq_len(n_cells)) {
      cluster <- cluster_assignments[i]
      initial_points[i, ] <- cluster_centers_U[cluster, ] +
        stats::rnorm(r, mean = 0, sd = cell_noise_level)
    }
    U_list[[j]] <- initial_points
    cell_cluster_assignments_list[[j]] <- cluster_assignments
  }

  # ---- Global V (Feature Loadings) per modality ----
  set.seed(seed)
  cluster_centers_V <- cluster_centers_U

  create_global_V <- function(p, n_clusters, r, cluster_centers) {
    cluster_assignments <- sample(seq_len(n_clusters), p, replace = TRUE)
    initial_points <- matrix(0, nrow = p, ncol = r)
    for (i in seq_len(p)) {
      cluster <- cluster_assignments[i]
      initial_points[i, ] <- cluster_centers[cluster, ] +
        stats::rnorm(r, sd = feature_noise_level)
    }
    list(V_global = initial_points, feature_clusters = cluster_assignments)
  }

  create_sample_V_list <- function(V_global, n_samples) {
    V_list <- list()
    p <- nrow(V_global)
    r <- ncol(V_global)
    for (i in seq_len(n_samples)) {
      V_sample <- V_global + matrix(stats::rnorm(p * r, mean = 0, sd = sample_noise_level),
                                    nrow = p, ncol = r)
      V_list[[i]] <- V_sample
    }
    V_list
  }

  # Build per-modality V data
  V_data_list <- list()
  V_global_list <- list()
  feature_cluster_list <- list()
  V_sample_list <- list()

  for (k in seq_len(n_modalities)) {
    vd <- create_global_V(n_features[k], n_clusters, r, cluster_centers_V)
    V_data_list[[k]] <- vd
    V_global_list[[k]] <- vd$V_global
    feature_cluster_list[[k]] <- vd$feature_clusters
    V_sample_list[[k]] <- create_sample_V_list(vd$V_global, n_samples)
  }

  # ---- Select signal features per modality ----
  set.seed(seed)
  signal_features_list <- list()
  for (k in seq_len(n_modalities)) {
    signal_features_list[[k]] <- sample(seq_len(n_features[k]),
                                         round(n_features[k] * signal_prop))
  }

  # Initialize expression matrix lists per modality
  expr_list <- vector("list", n_modalities)
  for (k in seq_len(n_modalities)) {
    expr_list[[k]] <- list()
  }

  # ---- Simulation Logic ----
  if (simulation_type == "DC") {
    constrained_shuffle <- function(feature_indices, cluster_assignments) {
      shuffled_indices <- feature_indices
      for (i in seq_along(feature_indices)) {
        original_idx <- feature_indices[i]
        original_cluster <- cluster_assignments[original_idx]
        candidate_pool <- feature_indices[cluster_assignments[feature_indices] != original_cluster]
        if (length(candidate_pool) == 0) {
          candidate_pool <- feature_indices[feature_indices != original_idx]
        }
        shuffled_indices[i] <- sample(candidate_pool, 1)
      }
      shuffled_indices
    }

    # Compute shuffled column indices per modality
    shuffled_cols_list <- list()
    for (k in seq_len(n_modalities)) {
      shuffled_cols_list[[k]] <- constrained_shuffle(
        signal_features_list[[k]], feature_cluster_list[[k]]
      )
    }

    # Apply permutation to V matrices for condition B
    V_sample_de_list <- V_sample_list
    for (k in seq_len(n_modalities)) {
      for (i in seq_len(n_samples)) {
        if (sample_metadata$condition[i] == "B") {
          V_sample_de_list[[k]][[i]][signal_features_list[[k]], ] <-
            V_sample_list[[k]][[i]][shuffled_cols_list[[k]], ]
        }
      }
    }

    # Generate expression matrices
    for (k in seq_len(n_modalities)) {
      for (i in seq_len(n_samples)) {
        expr_list[[k]][[i]] <- U_list[[i]] %*% sigma_matrix %*% t(V_sample_de_list[[k]][[i]])
      }
    }

    if (nonlinear) {
      for (k in seq_len(n_modalities)) {
        for (i in seq_len(n_samples)) {
          expr_list[[k]][[i]] <- apply_non_linearity(expr_list[[k]][[i]])
        }
      }
    }

    if (rescale) {
      for (k in seq_len(n_modalities)) {
        sf <- signal_features_list[[k]]
        target_means <- colMeans(U_list[[1]] %*% sigma_matrix %*%
                                   t(V_sample_list[[k]][[1]]))[sf]
        for (i in seq_len(n_samples)) {
          if (sample_metadata$condition[i] == "B") {
            current_means <- colMeans(expr_list[[k]][[i]])[sf]
            expr_list[[k]][[i]][, sf] <- sweep(expr_list[[k]][[i]][, sf], 2,
                                                target_means - current_means, "+")
          }
        }
      }
    }

  } else if (simulation_type == "DE") {
    for (k in seq_len(n_modalities)) {
      for (i in seq_len(n_samples)) {
        expr_list[[k]][[i]] <- U_list[[i]] %*% sigma_matrix %*% t(V_sample_list[[k]][[i]])
      }
    }

    if (nonlinear) {
      for (k in seq_len(n_modalities)) {
        for (i in seq_len(n_samples)) {
          expr_list[[k]][[i]] <- apply_non_linearity(expr_list[[k]][[i]])
        }
      }
    }

    for (k in seq_len(n_modalities)) {
      sf <- signal_features_list[[k]]
      for (i in seq_len(n_samples)) {
        if (sample_metadata$condition[i] == "B") {
          expr_list[[k]][[i]][, sf] <- expr_list[[k]][[i]][, sf] * fold_change
        }
      }
    }

  } else {
    stop("simulation_type must be either 'DC' or 'DE'")
  }

  # ---- Batch Effects ----
  if (add_batch) {
    for (k in seq_len(n_modalities)) {
      pk <- n_features[k]
      for (i in seq_len(n_samples)) {
        batch_num <- as.numeric(sub("batch_", "", sample_metadata$batch[i]))
        batch_scale_vec <- stats::rnorm(pk, batch_scale_mean, batch_scale_sd)
        batch_shift_vec <- stats::rnorm(pk, batch_shift_mean, batch_shift_sd) * (batch_num / n_batch)
        expr_list[[k]][[i]] <- sweep(expr_list[[k]][[i]], 2, batch_scale_vec, `*`)
        expr_list[[k]][[i]] <- sweep(expr_list[[k]][[i]], 2, batch_shift_vec, `+`)
      }
    }
  }

  # ---- Create Seurat Objects ----
  create_modality_matrix <- function(modality_matrix_list, feature_prefix,
                                      cell_cluster_assignments_list,
                                      feature_cluster_assignment, de_features) {
    n_feature <- ncol(modality_matrix_list[[1]])
    n_samp <- length(modality_matrix_list)

    all_cells <- c()
    all_features <- paste0(feature_prefix, "-", seq_len(n_feature))
    all_data <- list()
    all_metadata <- list()

    for (i in seq_len(n_samp)) {
      n_cell_sample <- nrow(modality_matrix_list[[i]])
      sample_cells <- paste0("sample", i, "_cell_", seq_len(n_cell_sample))
      all_cells <- c(all_cells, sample_cells)

      sample_meta_df <- data.frame(
        cell_id = sample_cells,
        sample_id = rep(sample_metadata$sample_id[i], n_cell_sample),
        condition = rep(sample_metadata$condition[i], n_cell_sample),
        batch = rep(sample_metadata$batch[i], n_cell_sample),
        true_cell_type = cell_cluster_assignments_list[[i]],
        stringsAsFactors = FALSE
      )
      all_metadata[[i]] <- sample_meta_df
      all_data[[i]] <- methods::as(modality_matrix_list[[i]], "dgCMatrix")
    }

    combined_sparse <- do.call(rbind, all_data)
    features_df <- data.frame(
      feature_id = all_features,
      feature_cluster = feature_cluster_assignment,
      is_de = (seq_len(length(all_features)) %in% de_features),
      stringsAsFactors = FALSE
    )
    combined_metadata <- do.call(rbind, all_metadata)

    rownames(combined_sparse) <- all_cells
    colnames(combined_sparse) <- all_features
    rownames(features_df) <- all_features
    rownames(combined_metadata) <- all_cells

    list(matrix = combined_sparse, features = features_df, cells = combined_metadata)
  }

  .make_seurat <- function(mod_data) {
    mat_t <- Matrix::t(mod_data$matrix)
    obj <- Seurat::CreateSeuratObject(
      counts = mat_t,
      data = mat_t,
      meta.data = mod_data$cells,
      assay = "originalexp"
    )
    obj[["originalexp"]] <- Seurat::AddMetaData(
      object = obj[["originalexp"]],
      metadata = mod_data$features
    )
    obj
  }

  # Build output
  seurat_out <- list()
  raw_out <- list()
  modality_names <- paste0("modality_", seq_len(n_modalities))

  for (k in seq_len(n_modalities)) {
    mod_data <- create_modality_matrix(
      expr_list[[k]],
      paste0("feature-", letters[k]),
      cell_cluster_assignments_list,
      feature_cluster_list[[k]],
      signal_features_list[[k]]
    )
    seurat_out[[modality_names[k]]] <- .make_seurat(mod_data)
    raw_out[[modality_names[k]]] <- expr_list[[k]]
  }

  list(
    seurat_list = seurat_out,
    raw_matrix_list = raw_out,
    sample_metadata = sample_metadata
  )
}
