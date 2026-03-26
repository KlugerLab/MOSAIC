#' @keywords internal
"_PACKAGE"

#' @importFrom Seurat ScaleData SplitObject DefaultAssay CreateSeuratObject AddMetaData
#' @importFrom RSpectra eigs_sym
#' @importFrom Matrix Matrix
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal theme_classic labs theme
#'   element_text stat_ellipse scale_color_manual scale_fill_manual
#' @importFrom rlang .data
#' @importFrom methods as
#' @importFrom parallel detectCores makeCluster stopCluster clusterExport
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach `%dopar%`
#' @importFrom vegan adonis2
#' @importFrom cluster silhouette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats cmdscale cor dist hclust cutree kmeans median rnorm
#'   p.adjust aggregate as.dist
NULL
