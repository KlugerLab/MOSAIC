#' Run MOSAIC Multi-Omics Co-Embedding
#'
#' The main entry point of the MOSAIC framework. For each sample, MOSAIC
#' constructs a sample-specific coupling matrix that captures intra- and
#' cross-modality feature interactions using cosine similarity. It then performs
#' spectral decomposition on each per-sample coupling matrix, aggregates the
#' resulting projection matrices across all samples, and applies a second-level
#' spectral decomposition to obtain a shared latent space. Finally, every
#' sample's coupling matrix is projected into this shared space, yielding a
#' per-sample feature embedding that can be used for downstream differential
#' connectivity analysis, subgroup detection, or clinical prediction.
#'
#' The function automatically handles 1, 2, or 3 modalities based on the length
#' of \code{seurat_list}. Each Seurat object should contain cells from multiple
#' samples (individuals), identified by a shared metadata column
#' (\code{sample_meta}). Condition labels (\code{condition_meta}) are extracted
#' and returned in the annotation table for downstream analyses.
#'
#' @param seurat_list A named list of 1 to 3 Seurat objects, one per modality.
#'   Each object must contain cells from all samples, with sample membership
#'   stored in the metadata column specified by \code{sample_meta}. The data
#'   should be normalized (e.g. via \code{Seurat::NormalizeData}) but does not
#'   need to be scaled; MOSAIC will call \code{Seurat::ScaleData} internally.
#'   Examples:
#'   \itemize{
#'     \item One modality: \code{list(RNA = seurat_rna)}
#'     \item Two modalities: \code{list(RNA = seurat_rna, ATAC = seurat_atac)}
#'     \item Three modalities: \code{list(RNA = seurat_rna, ATAC = seurat_atac,
#'       ADT = seurat_adt)}
#'   }
#' @param assays Character vector of assay names to use for each Seurat object,
#'   in the same order as \code{seurat_list}. For example,
#'   \code{c("RNA", "ATAC")} for a two-modality analysis. If \code{NULL}
#'   (default), the default assay of each Seurat object is used.
#' @param sample_meta Character string specifying the column name in the Seurat
#'   metadata that contains sample (individual) identifiers. All Seurat objects
#'   in \code{seurat_list} must share the same sample IDs in this column.
#'   Default: \code{"sample_id"}.
#' @param condition_meta Character string specifying the column name in the
#'   Seurat metadata that contains condition or group labels (e.g. "control" vs
#'   "disease"). Used to build the annotation table returned in the output.
#'   Default: \code{"condition"}.
#' @param n_eigen Integer specifying the number of eigenvalues/eigenvectors to
#'   compute in each spectral decomposition step. A larger value retains more
#'   spectral information but increases computation time. The kneedle algorithm
#'   (\code{\link{find_elbow_kneedle}}) is applied to automatically select the
#'   effective dimensionality from the top \code{n_eigen} eigenvalues.
#'   Default: \code{50}.
#' @param verbose Logical. If \code{TRUE} (default), print progress messages
#'   including the number of modalities and samples detected, per-sample
#'   processing status, and the selected ranks. Set to \code{FALSE} to
#'   suppress all messages.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{mosaic_embed_list}}{A named list (one entry per sample) of
#'     feature embedding matrices. Each matrix has dimensions
#'     \emph{n_features x r}, where \emph{n_features} is the total number of
#'     features across all modalities (stacked) and \emph{r} is the
#'     automatically selected latent dimensionality. Row \emph{i} of each matrix
#'     is the embedding of feature \emph{i} in that sample. Feature order is:
#'     all features from modality 1, then modality 2, then modality 3.}
#'   \item{\code{sample_eigen_list}}{A named list (one entry per sample)
#'     containing the per-sample eigendecomposition output (\code{$eigen}: an
#'     \code{eigs_sym} result with \code{$values} and \code{$vectors};
#'     \code{$r}: the automatically selected rank for that sample).}
#'   \item{\code{annotation}}{A data frame with one row per sample. Row names
#'     are sample IDs and column \code{Condition} contains the condition labels
#'     from \code{condition_meta}.}
#'   \item{\code{coupling_matrix_list}}{A named list (one entry per sample) of
#'     the raw coupling matrices (features x features). Useful for inspecting
#'     per-sample feature interaction structure or for downstream neighborhood
#'     analysis.}
#'   \item{\code{eigenvalues}}{Numeric vector of length \code{n_eigen}
#'     containing the eigenvalues from the aggregated (second-level) spectral
#'     decomposition. Can be visualized with \code{\link{plot_eigen}}.}
#' }
#'
#' @seealso
#' \code{\link{find_elbow_kneedle}} for the dimensionality selection algorithm,
#' \code{\link{run_DC_test}} for differential
#' connectivity testing on the output,
#' \code{\link{compute_module_similarity}} for subgroup detection,
#' \code{\link{plot_eigen}} for visualizing the eigenvalue spectrum.
#'
#' @examples
#' \dontrun{
#' # ---- One modality (RNA only) ----
#' result <- run_MOSAIC(
#'   list(RNA = seurat_rna),
#'   assays = c("RNA"),
#'   sample_meta = "sample_id",
#'   condition_meta = "condition"
#' )
#'
#' # Inspect eigenvalue spectrum
#' plot_eigen(result$eigenvalues)
#'
#' # Per-sample feature embedding for first sample
#' head(result$mosaic_embed_list[[1]])
#'
#' # ---- Two modalities (RNA + ATAC) ----
#' result <- run_MOSAIC(
#'   list(RNA = seurat_rna, ATAC = seurat_atac),
#'   assays = c("RNA", "ATAC"),
#'   sample_meta = "sample_id",
#'   condition_meta = "condition"
#' )
#'
#' # ---- Three modalities (RNA + ATAC + ADT) ----
#' result <- run_MOSAIC(
#'   list(RNA = seurat_rna, ATAC = seurat_atac, ADT = seurat_adt),
#'   assays = c("RNA", "ATAC", "ADT"),
#'   sample_meta = "sample_id",
#'   condition_meta = "condition"
#' )
#'
#' # ---- Downstream: Differential Connectivity ----
#' dc <- run_DC_test(
#'   result$mosaic_embed_list,
#'   n_sample = length(result$mosaic_embed_list),
#'   groups = result$annotation$Condition
#' )
#' }
#'
#' @export
run_MOSAIC <- function(seurat_list,
                       assays = NULL,
                       sample_meta = "sample_id",
                       condition_meta = "condition",
                       n_eigen = 50,
                       verbose = TRUE) {

  n_modalities <- length(seurat_list)
  if (n_modalities < 1 || n_modalities > 3) {
    stop("seurat_list must contain 1, 2, or 3 Seurat objects.")
  }

  # Default assays
  if (is.null(assays)) {
    assays <- sapply(seurat_list, Seurat::DefaultAssay)
  }
  if (length(assays) != n_modalities) {
    stop("Length of 'assays' must match length of 'seurat_list'.")
  }

  # Scale data for each modality
  for (i in seq_along(seurat_list)) {
    seurat_list[[i]] <- Seurat::ScaleData(
      seurat_list[[i]],
      features = rownames(seurat_list[[i]]),
      verbose = FALSE
    )
  }

  # Split by sample
  split_lists <- lapply(seurat_list, function(s) {
    Seurat::SplitObject(s, split.by = sample_meta)
  })
  n_samples <- length(split_lists[[1]])
  sample_names <- names(split_lists[[1]])

  if (verbose) {
    message(sprintf("MOSAIC: %d modality/modalities, %d samples detected",
                    n_modalities, n_samples))
  }

  # Extract matrices per sample: cells x features (transposed from features x cells)
  dataset_list <- vector("list", n_samples)
  for (i in seq_len(n_samples)) {
    matrices <- vector("list", n_modalities)
    for (m in seq_len(n_modalities)) {
      s <- split_lists[[m]][[i]]
      mat <- s[[assays[m]]]$scale.data
      mat <- t(as.matrix(mat))
      matrices[[m]] <- mat
    }
    dataset_list[[i]] <- matrices
  }

  # Eigendecomposition per sample
  eigen_output_list <- vector("list", n_samples)
  r_list <- numeric(n_samples)
  coupling_matrix_list <- vector("list", n_samples)

  for (i in seq_len(n_samples)) {
    if (verbose) message(sprintf("Processing sample %d / %d", i, n_samples))
    matrices <- dataset_list[[i]]

    eigen_out <- switch(
      as.character(n_modalities),
      "1" = .eigen_block_cosine_one(matrices[[1]], n_eigen = n_eigen),
      "2" = .eigen_block_cosine_two(matrices[[1]], matrices[[2]], n_eigen = n_eigen),
      "3" = .eigen_block_cosine_three(matrices[[1]], matrices[[2]], matrices[[3]], n_eigen = n_eigen)
    )

    r_list[i] <- eigen_out$r
    eigen_output_list[[i]] <- list(eigen = eigen_out$eigen, r = eigen_out$r)
    coupling_matrix_list[[i]] <- eigen_out$block_matrix
  }

  r <- round(stats::median(r_list))
  names(eigen_output_list) <- sample_names
  names(coupling_matrix_list) <- sample_names

  if (verbose) message(sprintf("Per-sample rank (median): %d", r))

  # Build annotation
  sample_condition_map <- unique(data.frame(
    sample_id = seurat_list[[1]][[sample_meta, drop = TRUE]],
    condition = seurat_list[[1]][[condition_meta, drop = TRUE]],
    stringsAsFactors = FALSE
  ))

  annotation <- data.frame(
    Condition = sample_condition_map$condition,
    stringsAsFactors = FALSE
  )
  rownames(annotation) <- sample_condition_map$sample_id
  annotation <- annotation[sample_names, , drop = FALSE]

  # Build per-sample projection matrices
  sample_projection_matrix <- vector("list", n_samples)
  for (i in seq_len(n_samples)) {
    r_to_use <- r_list[i]
    if (r_to_use == 1) r_to_use <- 2

    idx_range <- 1:r_to_use

    evecs <- eigen_output_list[[i]]$eigen$vectors[, idx_range, drop = FALSE]
    evals <- eigen_output_list[[i]]$eigen$values[idx_range]

    sample_projection_matrix[[i]] <- evecs %*% diag(evals, nrow = length(evals)) %*% t(evecs)
  }

  # Sum projection matrices and second eigendecomposition
  sum_projection_matrix <- Reduce("+", sample_projection_matrix)
  V_result <- RSpectra::eigs_sym(sum_projection_matrix, n_eigen, which = "LA")

  r2 <- find_elbow_kneedle(V_result$values[seq_len(n_eigen)])
  if (verbose) message(sprintf("Global rank: %d", r2))

  V_matrix <- V_result$vectors[, 1:r2, drop = FALSE] %*%
    diag(V_result$values[1:r2], nrow = r2)

  # Project coupling matrices into shared latent space
  mosaic_embed_list <- vector("list", n_samples)
  for (i in seq_len(n_samples)) {
    mosaic_embed_list[[i]] <- coupling_matrix_list[[i]] %*% V_matrix
  }
  names(mosaic_embed_list) <- sample_names

  if (verbose) message("Done.")

  return(list(
    mosaic_embed_list = mosaic_embed_list,
    sample_eigen_list = eigen_output_list,
    annotation = annotation,
    coupling_matrix_list = coupling_matrix_list,
    eigenvalues = V_result$values
  ))
}


# ---------------------------------------------------------------------------
# Internal: cosine similarity between columns of two cell x feature matrices
# Returns a features_A x features_B similarity matrix
# ---------------------------------------------------------------------------
.cosine_distance <- function(A, B) {
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


# ---------------------------------------------------------------------------
# Internal: eigendecomposition for 1 modality
# ---------------------------------------------------------------------------
.eigen_block_cosine_one <- function(matrix, n_eigen = 50) {
  block_matrix <- .cosine_distance(matrix, matrix)
  eigen_output <- RSpectra::eigs_sym(block_matrix, n_eigen, which = "LA")
  r <- find_elbow_kneedle(eigen_output$values)
  list(eigen = eigen_output, r = r, block_matrix = block_matrix)
}


# ---------------------------------------------------------------------------
# Internal: eigendecomposition for 2 modalities
# ---------------------------------------------------------------------------
.eigen_block_cosine_two <- function(matrix1, matrix2, whether_full = TRUE, n_eigen = 50) {
  A_B <- .cosine_distance(matrix1, matrix2)
  dim1 <- ncol(matrix1)
  dim2 <- ncol(matrix2)

  if (whether_full) {
    A_A <- .cosine_distance(matrix1, matrix1)
    B_B <- .cosine_distance(matrix2, matrix2)
    block_matrix <- rbind(
      cbind(A_A, A_B),
      cbind(t(A_B), B_B)
    )
  } else {
    block_matrix <- rbind(
      cbind(matrix(0, nrow = dim1, ncol = dim1), A_B),
      cbind(t(A_B), matrix(0, nrow = dim2, ncol = dim2))
    )
  }

  eigen_output <- RSpectra::eigs_sym(block_matrix, n_eigen, which = "LA")
  r <- find_elbow_kneedle(eigen_output$values)
  list(eigen = eigen_output, r = r, block_matrix = block_matrix)
}


# ---------------------------------------------------------------------------
# Internal: eigendecomposition for 3 modalities
# ---------------------------------------------------------------------------
.eigen_block_cosine_three <- function(matrix1, matrix2, matrix3, n_eigen = 50) {
  A_B <- .cosine_distance(matrix1, matrix2)
  A_C <- .cosine_distance(matrix1, matrix3)
  B_C <- .cosine_distance(matrix2, matrix3)
  A_A <- .cosine_distance(matrix1, matrix1)
  B_B <- .cosine_distance(matrix2, matrix2)
  C_C <- .cosine_distance(matrix3, matrix3)

  block_matrix <- rbind(
    cbind(A_A, A_B, A_C),
    cbind(t(A_B), B_B, B_C),
    cbind(t(A_C), t(B_C), C_C)
  )

  eigen_output <- RSpectra::eigs_sym(block_matrix, n_eigen, which = "LA")
  r <- find_elbow_kneedle(eigen_output$values)
  list(eigen = eigen_output, r = r, block_matrix = block_matrix)
}
