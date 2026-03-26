#' Plot Eigenvalue Scree Plot
#'
#' Displays the top eigenvalues as a scatter plot to help visualize the spectral
#' gap and inform dimensionality selection. The elbow in the scree plot
#' corresponds to the rank automatically selected by
#' \code{\link{find_elbow_kneedle}}.
#'
#' @param eigenvalues Numeric vector of eigenvalues, typically from
#'   \code{run_MOSAIC()$eigenvalues} or from a per-sample eigendecomposition.
#' @param n Integer. Number of eigenvalues to plot. Default:
#'   \code{min(50, length(eigenvalues))}.
#'
#' @return A \code{ggplot} object that can be further customized or saved with
#'   \code{ggsave()}.
#'
#' @examples
#' # Basic usage
#' eigenvalues <- c(10, 5, 3, 2, 1.5, 1.2, 1.1, 1.05, 1.01, 1.0)
#' plot_eigen(eigenvalues)
#'
#' # Plot only the first 5
#' plot_eigen(eigenvalues, n = 5)
#'
#' \dontrun{
#' # After running MOSAIC
#' result <- run_MOSAIC(list(RNA = seurat_rna),
#'                      sample_meta = "sample_id",
#'                      condition_meta = "condition")
#' plot_eigen(result$eigenvalues)
#' }
#'
#' @seealso \code{\link{find_elbow_kneedle}} for automatic dimensionality
#'   selection from the same eigenvalues.
#'
#' @export
plot_eigen <- function(eigenvalues, n = NULL) {
  if (is.null(n)) n <- min(50, length(eigenvalues))
  eigenvalues <- eigenvalues[seq_len(n)]

  df <- data.frame(x = seq_len(n), y = eigenvalues)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "n", y = "eigenvalues") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "none"
    )
}


#' MDS Plot with Cluster Coloring
#'
#' Performs classical multidimensional scaling (MDS) on a similarity or
#' dissimilarity matrix and displays a 2D scatter plot colored by group labels,
#' with optional confidence ellipses. This is the primary visualization for
#' MOSAIC's per-feature sample similarity matrices and module-level subgroup
#' structure.
#'
#' @param similarity_matrix A square numeric matrix of pairwise similarities or
#'   distances. By default, this is treated as a similarity matrix and converted
#'   to dissimilarity via \code{1 - x}. Set \code{dis = TRUE} if the input is
#'   already a dissimilarity/distance matrix.
#' @param title Character string for the plot title.
#' @param cluster Vector of group labels (character, factor, or numeric), one
#'   per row/column of \code{similarity_matrix}. Used for point coloring and
#'   ellipses.
#' @param size Numeric. Point size. Default: \code{1}.
#' @param custom_colors Named character vector mapping group labels to colors.
#'   Names must match the levels of \code{cluster}. Example:
#'   \code{c("Control" = "#56B4E9", "Disease" = "#E69F00")}. If \code{NULL}
#'   (default), the \code{Set2} palette from RColorBrewer is used.
#' @param dis Logical. If \code{TRUE}, \code{similarity_matrix} is treated as a
#'   dissimilarity matrix directly (no \code{1 - x} conversion). Default:
#'   \code{FALSE}.
#' @param idx1 Integer. MDS dimension to use for the x-axis. Default: \code{1}.
#' @param idx2 Integer. MDS dimension to use for the y-axis. Default: \code{2}.
#' @param add_ellipse Logical. If \code{TRUE} (default), add t-distribution
#'   confidence ellipses around each group.
#' @param ellipse_level Numeric between 0 and 1. Confidence level for the
#'   ellipses. Default: \code{0.66}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' # Create a toy similarity matrix
#' set.seed(42)
#' sim <- matrix(runif(100, 0.5, 1), 10, 10)
#' sim <- (sim + t(sim)) / 2  # make symmetric
#' diag(sim) <- 1
#' rownames(sim) <- colnames(sim) <- paste0("S", 1:10)
#' groups <- rep(c("A", "B"), each = 5)
#'
#' plot_mds_cluster(sim, "Example MDS", cluster = groups)
#'
#' # With custom colors and no ellipses
#' plot_mds_cluster(sim, "Custom Colors", cluster = groups,
#'                  custom_colors = c("A" = "#56B4E9", "B" = "#E69F00"),
#'                  add_ellipse = FALSE, size = 3)
#'
#' \dontrun{
#' # Visualize a DC feature's sample similarity
#' plot_mds_cluster(
#'   dc_result$similarity_matrix_list[[feature_idx]],
#'   title = "STAT5B Connectivity",
#'   cluster = paste0("Day_", dc_result$group_list[[feature_idx]]),
#'   custom_colors = c("Day_0" = "#56B4E9", "Day_7" = "#E69F00")
#' )
#' }
#'
#' @seealso \code{\link{compute_module_similarity}} for generating similarity
#'   matrices for subgroup detection.
#'
#' @export
plot_mds_cluster <- function(similarity_matrix, title, cluster, size = 1,
                             custom_colors = NULL, dis = FALSE, idx1 = 1, idx2 = 2,
                             add_ellipse = TRUE, ellipse_level = 0.66) {
  if (!dis) {
    dissimilarity_matrix <- 1 - similarity_matrix
  } else {
    dissimilarity_matrix <- similarity_matrix
  }

  mds_result <- stats::cmdscale(dissimilarity_matrix, k = 5)

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

  p <- ggplot2::ggplot(umap_df, ggplot2::aes(
    x = .data$UMAP1, y = .data$UMAP2,
    color = .data$cluster, fill = .data$cluster
  )) +
    ggplot2::geom_point(size = size) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = paste0("MDS_", idx1),
      y = paste0("MDS_", idx2)
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::scale_color_manual(values = custom_colors)

  if (add_ellipse) {
    p <- p +
      ggplot2::stat_ellipse(
        ggplot2::aes(fill = .data$cluster),
        geom = "polygon",
        alpha = 0.2,
        level = ellipse_level,
        type = "t",
        linetype = 2,
        linewidth = 0.5
      ) +
      ggplot2::scale_fill_manual(values = custom_colors)
  }

  return(p)
}
