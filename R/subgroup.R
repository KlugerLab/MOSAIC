#' Compute Module-Level Sample Similarity Matrix
#'
#' Given a set of feature indices defining a module (e.g. a cluster of
#' co-regulated features), computes a sample-by-sample cosine similarity matrix
#' based on the concatenated feature embeddings within that module. This
#' captures how similar two samples are in terms of the connectivity structure
#' of the module's features.
#'
#' This is the first step in MOSAIC's unsupervised subgroup detection pipeline:
#' compute module similarity, then test for a balanced partition with
#' \code{\link{find_partition_hclust}}.
#'
#' @param projected_list Named list of per-sample projected feature matrices
#'   (features x latent dims), as returned by
#'   \code{\link{run_MOSAIC}()$mosaic_embed_list}. Each element is a matrix
#'   where row \emph{i} is feature \emph{i}'s embedding in that sample.
#' @param feature_idx Integer vector of feature indices (row indices) belonging
#'   to the module. For example, \code{which(feature_clusters == module_id)}.
#'
#' @return A symmetric numeric matrix of dimension \emph{n_samples x n_samples},
#'   with cosine similarity values (range -1 to 1). Row and column names
#'   correspond to sample names from \code{projected_list}.
#'
#' @seealso \code{\link{find_partition_hclust}} for testing whether the
#'   similarity matrix reveals distinct subgroups,
#'   \code{\link{plot_mds_cluster}} for visualizing the sample similarity.
#'
#' @examples
#' \dontrun{
#' # Run MOSAIC
#' result <- run_MOSAIC(
#'   list(RNA = seurat_rna, ATAC = seurat_atac),
#'   assays = c("RNA", "GeneACT"),
#'   sample_meta = "sample_name",
#'   condition_meta = "condition"
#' )
#'
#' # Suppose feature_clusters is a vector of module assignments
#' # Compute similarity for module 5
#' module_idx <- which(feature_clusters == 5)
#' sim_mat <- compute_module_similarity(
#'   result$mosaic_embed_list,
#'   feature_idx = module_idx
#' )
#'
#' # Visualize
#' plot_mds_cluster(sim_mat, "Module 5",
#'                  cluster = result$annotation$Condition)
#' }
#'
#' @export
compute_module_similarity <- function(projected_list, feature_idx) {
  n_sample <- length(projected_list)
  sim_mat <- matrix(0, nrow = n_sample, ncol = n_sample)

  for (i in seq_len(n_sample)) {
    for (j in i:n_sample) {
      x <- as.vector(projected_list[[i]][feature_idx, , drop = FALSE])
      y <- as.vector(projected_list[[j]][feature_idx, , drop = FALSE])
      s <- .cosine_sim_vectors(x, y)
      sim_mat[i, j] <- s
      sim_mat[j, i] <- s
    }
  }

  rownames(sim_mat) <- names(projected_list)
  colnames(sim_mat) <- names(projected_list)
  return(sim_mat)
}


#' Find a Balanced Two-Way Partition via Hierarchical Clustering
#'
#' Performs hierarchical clustering on a distance matrix and cuts the
#' dendrogram into two groups. A partition is considered "balanced" if both
#' groups contain at least \code{max(3, ceiling(n * min_group_frac))} samples.
#' For balanced partitions, the average silhouette width is computed to
#' quantify separation quality.
#'
#' This function is used in MOSAIC's subgroup detection pipeline to test
#' whether a feature module defines meaningful patient subtypes: compute a
#' module similarity matrix with \code{\link{compute_module_similarity}},
#' convert to distance, and test for a balanced partition here.
#'
#' @param dist_mat A \code{\link[stats]{dist}} object representing pairwise
#'   distances between samples. Typically obtained from a module similarity
#'   matrix via \code{as.dist(1 - sim_mat)}.
#' @param min_group_frac Numeric value between 0 and 0.5 specifying the minimum
#'   fraction of total samples required in each group. The actual minimum group
#'   size is \code{max(3, ceiling(n_samples * min_group_frac))}. Default:
#'   \code{0.25}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{groups}}{Named integer vector of group assignments (1 or 2),
#'     one per sample. \code{NULL} if no balanced partition was found.}
#'   \item{\code{silhouette}}{Average silhouette width (higher = better
#'     separation). \code{NA} if the partition is imbalanced.}
#'   \item{\code{hclust}}{The \code{\link[stats]{hclust}} object, which can be
#'     passed to \code{pheatmap} for consistent dendrogram ordering.}
#'   \item{\code{balanced}}{Logical indicating whether a balanced partition was
#'     found.}
#' }
#'
#' @seealso \code{\link{compute_module_similarity}} to create the input
#'   similarity matrix, \code{\link{plot_mds_cluster}} to visualize subgroups.
#'
#' @examples
#' \dontrun{
#' # Compute module similarity
#' sim_mat <- compute_module_similarity(
#'   result$mosaic_embed_list,
#'   feature_idx = which(feature_clusters == module_id)
#' )
#'
#' # Test for subgroups
#' partition <- find_partition_hclust(as.dist(1 - sim_mat))
#'
#' if (partition$balanced) {
#'   cat("Silhouette:", partition$silhouette, "\n")
#'   cat("Group sizes:", table(partition$groups), "\n")
#'
#'   # Visualize
#'   plot_mds_cluster(sim_mat, "Subgroups",
#'                    cluster = partition$groups)
#' }
#'
#' # Permutation test for significance
#' n_perm <- 1000
#' null_sils <- numeric(n_perm)
#' for (p in seq_len(n_perm)) {
#'   rand_idx <- sample(total_features, length(module_idx))
#'   rand_sim <- compute_module_similarity(projected_list, rand_idx)
#'   rand_part <- find_partition_hclust(as.dist(1 - rand_sim))
#'   null_sils[p] <- ifelse(rand_part$balanced, rand_part$silhouette, 0)
#' }
#' pval <- (sum(null_sils >= partition$silhouette) + 1) / (n_perm + 1)
#' }
#'
#' @export
find_partition_hclust <- function(dist_mat, min_group_frac = 0.25) {
  n_samples <- attr(dist_mat, "Size")
  min_group_size <- max(3, ceiling(n_samples * min_group_frac))

  hc <- stats::hclust(stats::dist(as.matrix(dist_mat), method = "euclidean"),
                       method = "complete")
  groups <- stats::cutree(hc, k = 2)
  group_sizes <- table(groups)

  if (min(group_sizes) < min_group_size) {
    return(list(groups = NULL, silhouette = NA, hclust = hc, balanced = FALSE))
  }

  sil_scores <- cluster::silhouette(groups, dist_mat)
  avg_sil <- mean(sil_scores[, "sil_width"])

  list(groups = groups, silhouette = avg_sil, hclust = hc, balanced = TRUE)
}
