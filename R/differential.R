#' Run Differential Connectivity Test
#'
#' The core function for MOSAIC's differential connectivity (DC) analysis.
#' For each feature, this function collects the feature's embedding vector
#' across all samples, builds a sample-by-sample distance matrix from these
#' vectors, and then runs PERMANOVA (\code{vegan::adonis2}) to test whether
#' the embedding structure differs significantly between conditions. It also
#' computes silhouette scores to quantify the separation between condition
#' groups in the feature's embedding space.
#'
#' Samples where a feature has an all-zero embedding (e.g. due to dropout) are
#' automatically excluded on a per-feature basis. The function runs in parallel
#' across features using the \code{foreach} / \code{doParallel} framework.
#'
#' To obtain calibrated p-values, run this function on both the real data and a
#' label-shuffled null, then compare the F-statistics using
#' \code{\link{calculate_empirical_pvalue}}.
#'
#' @param mosaic_embed_list Named list of per-sample projected feature matrices
#'   (features x latent dims), as returned by
#'   \code{\link{run_MOSAIC}()$mosaic_embed_list}. Each element is a matrix
#'   where row \emph{i} is feature \emph{i}'s embedding in that sample.
#' @param n_sample Integer. Number of samples (length of the input list).
#' @param groups Character or factor vector of condition labels, one per sample,
#'   in the same order as the names of \code{mosaic_embed_list}. Typically
#'   obtained from \code{result$annotation$Condition}.
#' @param n_cores Integer. Number of parallel cores to use. If \code{NULL}
#'   (default), uses \code{parallel::detectCores() - 1}.
#' @param dist_method Character string specifying the distance metric for the
#'   per-feature sample distance matrix:
#'   \itemize{
#'     \item \code{"euclidean"} (default): Euclidean distance.
#'     \item \code{"cosine"}: Cosine dissimilarity (1 - cosine similarity).
#'   }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{pvalue_list}}{Numeric vector of PERMANOVA p-values, one per
#'     feature.}
#'   \item{\code{r2_list}}{Numeric vector of PERMANOVA R-squared values (effect
#'     sizes), one per feature.}
#'   \item{\code{F_stats_list}}{Numeric vector of PERMANOVA F-statistics, one
#'     per feature. Used as the test statistic for empirical p-value
#'     calibration.}
#'   \item{\code{silhouette_score_list}}{List of average silhouette scores per
#'     feature.}
#'   \item{\code{similarity_matrix_list}}{List of per-feature sample-by-sample
#'     similarity/distance matrices. Each matrix has dimension
#'     \emph{n_valid_samples x n_valid_samples} (after zero-row removal).}
#'   \item{\code{group_list}}{List of condition label vectors per feature (after
#'     zero-row removal).}
#' }
#'
#' @seealso
#' \code{\link{run_MOSAIC}} to generate the input embedding,
#' \code{\link{calculate_empirical_pvalue}} to compute calibrated p-values
#' by comparing observed F-statistics against a shuffled null distribution.
#'
#' @examples
#' \dontrun{
#' # Step 1: Run MOSAIC embedding
#' result <- run_MOSAIC(
#'   list(RNA = seurat_rna, ADT = seurat_adt),
#'   assays = c("SCT", "ADT"),
#'   sample_meta = "sample_id",
#'   condition_meta = "time"
#' )
#'
#' # Step 2: Run DC analysis on real data
#' n_sample <- length(result$mosaic_embed_list)
#' dc_result <- run_DC_test(
#'   result$mosaic_embed_list,
#'   n_sample = n_sample,
#'   groups = result$annotation$Condition
#' )
#'
#' # Step 3: Run on shuffled labels for null distribution
#' shuffle_dc <- run_DC_test(
#'   shuffle_mosaic$mosaic_embed_list,
#'   n_sample = n_sample,
#'   groups = shuffle_mosaic$annotation$Condition
#' )
#'
#' # Step 4: Compute empirical p-values
#' F_obs <- unlist(dc_result$F_stats_list)
#' F_null <- unlist(shuffle_dc$F_stats_list)
#' pvalues <- sapply(F_obs, function(x)
#'   calculate_empirical_pvalue(x, F_null))
#'
#' # Features with significant DC
#' dc_features <- which(pvalues < 0.05)
#' }
#'
#' @export
run_DC_test <- function(
    mosaic_embed_list,
    n_sample,
    groups,
    n_cores = NULL,
    dist_method = "euclidean"
) {
  # Clean up any existing connections
  existing_connections <- showConnections()
  for (i in seq_len(nrow(existing_connections))) {
    try(close(getConnection(i)), silent = TRUE)
  }

  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 1
  }
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  n_features <- nrow(mosaic_embed_list[[1]])

  # Helper functions to export to workers
  .cosine_distance_local <- .cosine_distance
  .remove_zero_rows_local <- .remove_zero_rows

  parallel::clusterExport(cl, c(
    "mosaic_embed_list", "n_sample", "groups",
    ".cosine_distance_local", ".remove_zero_rows_local"
  ), envir = environment())

  pkg_list <- c("vegan", "stats", "cluster")

  idx <- NULL  # avoid R CMD check NOTE
  results <- foreach::foreach(idx = seq_len(n_features), .packages = pkg_list) %dopar% {
    if (idx %% 100 == 0) {
      cat("Processing feature:", idx, "\n")
    }

    # Collect embedding for this feature across samples
    combined_coembed_matrix <- matrix(
      0, nrow = n_sample,
      ncol = ncol(mosaic_embed_list[[1]])
    )
    for (i in seq_len(n_sample)) {
      combined_coembed_matrix[i, ] <- mosaic_embed_list[[i]][idx, ]
    }
    rownames(combined_coembed_matrix) <- names(mosaic_embed_list)

    # Remove zero rows
    remove_result <- .remove_zero_rows_local(combined_coembed_matrix)
    combined_coembed_matrix <- remove_result[[1]]
    remove_idx <- remove_result[[2]]

    # Pairwise distance/similarity
    dim_idx <- nrow(combined_coembed_matrix)
    similarity_matrix <- matrix(0, nrow = dim_idx, ncol = dim_idx)

    for (i in seq_len(dim_idx)) {
      for (j in seq_len(dim_idx)) {
        x <- combined_coembed_matrix[i, ]
        y <- combined_coembed_matrix[j, ]
        if (dist_method == "cosine") {
          similarity_matrix[i, j] <- .cosine_distance_local(
            matrix(x, nrow = 1), matrix(y, nrow = 1)
          )
        } else {
          similarity_matrix[i, j] <- as.numeric(stats::dist(rbind(x, y)))
        }
      }
    }
    rownames(similarity_matrix) <- rownames(combined_coembed_matrix)
    colnames(similarity_matrix) <- rownames(combined_coembed_matrix)

    if (dist_method == "cosine") {
      dissimilarity_matrix <- stats::as.dist(1 - similarity_matrix)
    } else {
      dissimilarity_matrix <- stats::as.dist(similarity_matrix)
    }

    condition_numeric <- as.integer(as.factor(groups[!remove_idx]))
    silhouette_scores <- cluster::silhouette(condition_numeric, dissimilarity_matrix)
    average_silhouette_score <- mean(silhouette_scores[, "sil_width"])

    permanova_result <- vegan::adonis2(
      dissimilarity_matrix ~ groups[!remove_idx],
      permutations = 999
    )

    list(
      pvalue = permanova_result$`Pr(>F)`[1],
      r2 = permanova_result$R2[1],
      similarity_matrix = similarity_matrix,
      group = groups[!remove_idx],
      F_statistics = permanova_result$`F`[1],
      average_silhouette_score = average_silhouette_score
    )
  }

  parallel::stopCluster(cl)

  list(
    pvalue_list = sapply(results, function(x) x$pvalue),
    r2_list = sapply(results, function(x) x$r2),
    F_stats_list = sapply(results, function(x) x$F_statistics),
    silhouette_score_list = lapply(results, function(x) x$average_silhouette_score),
    similarity_matrix_list = lapply(results, function(x) x$similarity_matrix),
    group_list = lapply(results, function(x) x$group)
  )
}
