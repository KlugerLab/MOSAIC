# Supplemental Figure S8: pseudobulk PCA on the HIV+ subgroups
# Input:    PFC_L23IT_RNA, PFC_L23IT_ATAC, results/fig5/fig5_session_objects.rds
# Output:   figures/figS8A.png, figures/figS8B.png
# Requires: fig5/fig5_figS6_S7_hiv_subgroup_discovery.R

source("config.R")
source("R/functions.R")

library(Seurat)
library(pheatmap)
library(cluster)

fig5_objects <- readRDS(file.path(RESULTS_DIR, "fig5", "fig5_session_objects.rds"))
MODULE_ID           <- fig5_objects$MODULE_ID
module_sim_matrices <- fig5_objects$module_sim_matrices
module_hclust_list  <- fig5_objects$module_hclust_list
module_assignments  <- fig5_objects$module_assignments

## L2/3 IT objects: samples with < 10 cells and CTR removed
rna_seurat_subset <- readRDS(PFC_L23IT_RNA)
atac_seurat_subset <- readRDS(PFC_L23IT_ATAC)

DefaultAssay(atac_seurat_subset) <- "GeneACT"
atac_seurat_subset[["ATAC"]] <- NULL

sample_counts <- table(rna_seurat_subset$sample_name)
samples_to_remove <- names(sample_counts[sample_counts < 10])

rna_seurat_subset <- subset(rna_seurat_subset, subset = sample_name %ni% samples_to_remove)
atac_seurat_subset <- subset(atac_seurat_subset, subset = sample_name %ni% samples_to_remove)

rna_seurat_subset <- subset(rna_seurat_subset, subset = condition_new %ni% c("CTR"))
atac_seurat_subset <- subset(atac_seurat_subset, subset = condition_new %ni% c("CTR"))

stopifnot(setequal(unique(as.character(rna_seurat_subset$sample_name)), fig5_objects$hiv_samples))

## Pseudobulk PCA

compute_pseudobulk_pca_similarity <- function(seurat_obj, assay = "RNA", n_pcs = 30) {
  # Aggregate to pseudobulk per sample
  DefaultAssay(seurat_obj) <- assay

  # Get sample IDs
  sample_name <- unique(seurat_obj$sample_name)

  # Build pseudobulk matrix: samples x features
  pseudobulk <- do.call(rbind, lapply(sample_name, function(s) {
    cells <- which(seurat_obj$sample_name == s)
    if (length(cells) > 1) {
      colMeans(GetAssayData(seurat_obj, slot = "data")[, cells])
    } else {
      GetAssayData(seurat_obj, slot = "data")[, cells]
    }
  }))
  rownames(pseudobulk) <- sample_name

  # PCA
  pca <- prcomp(pseudobulk, scale. = TRUE, center = TRUE)

  # Distance in PC space
  pca_coords <- pca$x[, 1:min(n_pcs, ncol(pca$x))]
  dist_mat <- dist(pca_coords, method = "euclidean")

  # Convert to similarity (for heatmap comparability with MOSAIC)
  dist_matrix <- as.matrix(dist_mat)
  max_dist <- max(dist_matrix)
  sim_mat <- 1 - dist_matrix / max_dist

  return(list(sim_mat = sim_mat, dist_mat = dist_mat, pca = pca))
}

# Run on RNA
rna_result <- compute_pseudobulk_pca_similarity(rna_seurat_subset, assay = "RNA")

# Run on ATAC (gene activity scores)
atac_result <- compute_pseudobulk_pca_similarity(atac_seurat_subset, assay = "GeneACT")

## Heatmaps with subgroup annotation
# Get Module 27 ordering and groups
m27_sim <- module_sim_matrices[[as.character(MODULE_ID)]]
m27_hc <- module_hclust_list[[as.character(MODULE_ID)]]
m27_groups <- module_assignments[[as.character(MODULE_ID)]]
m27_order <- m27_hc$order
ordered_samples <- rownames(m27_sim)[m27_order]

plot_with_mosaic_order <- function(sim_mat, ordered_samples, groups, title_text) {
  common <- intersect(ordered_samples, rownames(sim_mat))
  ordered_common <- ordered_samples[ordered_samples %in% common]
  sim_ordered <- sim_mat[ordered_common, ordered_common]

  annot <- data.frame(
    Subgroup = ifelse(groups[ordered_common] == 1, "Group1", "Group2"),
    row.names = ordered_common,
    Condition = rep("HIV", length(ordered_common))
  )

  annot_colors <- list(
    Subgroup = c("Group1" = "#4DBBD5", "Group2" = "#E64B35"),
    Condition = c("HIV" = "#FC8D62")
  )

  # Quantile breaks, simple sequence if quantiles collapse
  vals <- sim_ordered[upper.tri(sim_ordered)]
  breaks <- unique(quantile(vals, probs = seq(0, 1, length.out = 101)))
  if (length(breaks) < 3) {
    breaks <- seq(min(vals) - 0.01, max(vals) + 0.01, length.out = 100)
  }

  pheatmap(
    sim_ordered,
    annotation_row = annot,
    annotation_col = annot,
    annotation_colors = annot_colors,
    breaks = breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = title_text,
    fontsize = 10
  )
}

# Silhouette of the Module 27 partition in each pseudobulk PCA space
rna_common <- intersect(names(m27_groups), rownames(rna_result$sim_mat))
rna_dist <- as.dist(1 - rna_result$sim_mat[rna_common, rna_common])
rna_sil <- mean(silhouette(m27_groups[rna_common], rna_dist)[, "sil_width"])

atac_common <- intersect(names(m27_groups), rownames(atac_result$sim_mat))
atac_dist <- as.dist(1 - atac_result$sim_mat[atac_common, atac_common])
atac_sil <- mean(silhouette(m27_groups[atac_common], atac_dist)[, "sil_width"])

cat(sprintf("  MOSAIC Module %d: %.3f\n", MODULE_ID, fig5_objects$module_silhouette))
cat(sprintf("  Pseudobulk PCA RNA: %.3f\n", rna_sil))
cat(sprintf("  Pseudobulk PCA ATAC: %.3f\n", atac_sil))

pp1 <- plot_with_mosaic_order(rna_result$sim_mat, ordered_samples, m27_groups,
                              sprintf("Pseudobulk PCA on RNA (silhouette=%.3f)", rna_sil))
pp2 <- plot_with_mosaic_order(atac_result$sim_mat, ordered_samples, m27_groups,
                              sprintf("Pseudobulk PCA on ATAC (silhouette=%.3f)", atac_sil))

# Supplemental Figure S8A
png(file.path(FIG_DIR, "figS8A.png"),width = 6.2, height = 5.5, units = "in", res = 600)
print(pp1)
dev.off()

# Supplemental Figure S8B
png(file.path(FIG_DIR, "figS8B.png"),width = 6.2, height = 5.5, units = "in", res = 600)
print(pp2)
dev.off()
