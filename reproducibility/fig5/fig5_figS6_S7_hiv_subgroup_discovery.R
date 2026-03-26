# Figure 5A, 5B, 5E, Supplemental Figures S6, S7: HIV+ subgroup discovery
# Input:    PFC_L23IT_RNA, PFC_L23IT_ATAC
# Output:   results/fig5/fig5_session_objects.rds, figures/fig5A.png, fig5B.png, fig5E.png, figS6.png, figS7.png

source("config.R")
source("R/functions.R")

dir.create(file.path(RESULTS_DIR, "fig5"), showWarnings = FALSE, recursive = TRUE)

rna_seurat_subset <- readRDS(PFC_L23IT_RNA)
atac_seurat_subset <- readRDS(PFC_L23IT_ATAC)

DefaultAssay(atac_seurat_subset) <- "GeneACT"
atac_seurat_subset[["ATAC"]] <- NULL

sample_counts <- table(rna_seurat_subset$sample_name)
samples_to_remove <- names(sample_counts[sample_counts < 10])
print(samples_to_remove)

rna_seurat_subset <- subset(rna_seurat_subset, subset = sample_name %ni% samples_to_remove)
atac_seurat_subset <- subset(atac_seurat_subset, subset = sample_name %ni% samples_to_remove)

rna_seurat_subset <- subset(rna_seurat_subset, subset = condition_new %ni% c("CTR"))
atac_seurat_subset <- subset(atac_seurat_subset, subset = condition_new %ni% c("CTR"))

counts <- GetAssayData(rna_seurat_subset, assay = "RNA", slot = "counts")
sample_ids <- unique(as.character(rna_seurat_subset$sample_name))
list_of_expressed_genes <- lapply(sample_ids, function(id) {
  cells_in_sample <- which(rna_seurat_subset$sample_name == id)
  sample_counts <- counts[, cells_in_sample, drop = FALSE]
  expressed_genes <- rownames(sample_counts)[rowSums(sample_counts) > 3]
  return(expressed_genes)
})
genes_to_keep <- Reduce(intersect, list_of_expressed_genes)
rna_seurat_subset <- FindVariableFeatures(rna_seurat_subset,nfeatures=10000)
genes_to_keep <- genes_to_keep[which(genes_to_keep %in% VariableFeatures(rna_seurat_subset))]
length(genes_to_keep)

counts <- GetAssayData(atac_seurat_subset, assay = "GeneACT", slot = "counts")
sample_ids <- unique(as.character(atac_seurat_subset$sample_name))
list_of_expressed_genes <- lapply(sample_ids, function(id) {
  cells_in_sample <- which(atac_seurat_subset$sample_name == id)
  sample_counts <- counts[, cells_in_sample, drop = FALSE]
  expressed_genes <- rownames(sample_counts)[rowSums(sample_counts) > 5]
  return(expressed_genes)
})

atac_seurat_subset <- FindVariableFeatures(atac_seurat_subset,nfeatures=3000)
peak_to_keep <- Reduce(intersect, list_of_expressed_genes)
peak_to_keep <- peak_to_keep[which(peak_to_keep %in% VariableFeatures(atac_seurat_subset))]
length(peak_to_keep)

seurat_atac <- atac_seurat_subset
seurat_atac <- subset(seurat_atac, features = peak_to_keep)

seurat_rna <- rna_seurat_subset
seurat_rna <- subset(seurat_rna, features = genes_to_keep)

seurat_rna <- NormalizeData(seurat_rna)
seurat_atac <- NormalizeData(seurat_atac)

feature_names <- c(rownames(seurat_rna@assays$RNA$counts),paste0("ATAC-",rownames(seurat_atac@assays$GeneACT$counts)))

## MOSAIC embedding

run_our_method_two <- function(rna_seurat, atac_seurat, assay1 = "SCT", assay2="ADT",sample_meta="sample_id",condition_meta="condition",block_method="cosine"){
  # Global scale
  rna_seurat <- ScaleData(rna_seurat,features=rownames(rna_seurat))
  atac_seurat <- ScaleData(atac_seurat,features=rownames(atac_seurat))
  
  rna.list <- SplitObject(rna_seurat, split.by = sample_meta)
  atac.list <- SplitObject(atac_seurat, split.by = sample_meta)
  
  dataset_list <- list()
  for (i in 1:length(rna.list)){
    print(i)
    gene_seurat <- rna.list[[i]]
    gene_matrix <- gene_seurat[[assay1]]$scale.data
    gene_matrix <- t(as.matrix(gene_matrix))
    
    atac_seurat <- atac.list[[i]]
    atac_matrix <- atac_seurat[[assay2]]$scale.data
    atac_matrix <- t(as.matrix(atac_matrix))
    
    dataset_list[[i]] <- list(gene_matrix, atac_matrix)
  }
  
  eigen_output_cosine_list <- list()
  r_cosine_list <- c()
  block_matrix_list <- list()
  for (i in 1:length(dataset_list)){
    print(i)
    matrix_set <- dataset_list[[i]]
    
    if (block_method == "cosine"){
      print("cosine")
      eigen_output_cosine <- eigen_block_cosine_two(matrix_set[[1]],matrix_set[[2]],)
    } else if (block_method == "gaussian"){
      print("gaussian")
      eigen_output_cosine <- eigen_block_gaussian_two(matrix_set[[1]],matrix_set[[2]])
    } else if (block_method == "emd"){
      print("emd")
      eigen_output_cosine <- eigen_block_emd_two(matrix_set[[1]],matrix_set[[2]])
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
  names(block_matrix_list) <- names(rna.list)
  
  sample_name <- names(eigen_output_cosine_list)
  
  sample_condition_map <- unique(data.frame(
    sample_id = rna_seurat[[sample_meta]],
    condition = rna_seurat[[condition_meta]]
  ))
  
  colnames(sample_condition_map) <- c("sample_id","condition")
  
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
  
  return(list(projected_combed_list=projected_combed_list,eigen_output_cosine_list=eigen_output_cosine_list,annotation=annotation,block_matrix_list=block_matrix_list))
}

no_shuffle_result <- run_our_method_two(seurat_rna,seurat_atac,sample_meta="sample_name",assay1 = "RNA",assay2="GeneACT",condition_meta="condition_new",block_method = "cosine")

# Per-feature sample-by-sample similarity matrices
n_sample <- length(no_shuffle_result$projected_combed_list)
peranova_result <- cal_peranova_per_feature_parallel_no_zero(no_shuffle_result$projected_combed_list,n_sample,
                                                             groups = no_shuffle_result$annotation[names(no_shuffle_result$eigen_output_cosine_list),"Condition"],
                                                             n_cores = N_CORES,
                                                             dist_method = "cosine",num_permutation=9999,whether_permanova = FALSE)

library(dynamicTreeCut)
library(cluster)       # silhouette
library(lsa)           # cosine
library(ggplot2)
library(pheatmap)

## Parameters
MIN_MODULE_SIZE   <- 20    # minimum features per module
MIN_GROUP_FRAC    <- 0.25  # each group must have >= this fraction of samples
SIGCLUST_NSIM     <- 1000  # simulations for sigclust
PADJ_THRESHOLD    <- 0.05  # significance after BH correction
DEEP_SPLIT        <- 2     # 0-4, lower = fewer larger modules
N_PERM <- 1000             # random feature sets per module (permutation null)
MODULE_ID         <- 27    # module characterized in Figure 5
## Step 1: Feature stratification profiles
similarity_matrix_per_feature <- peranova_result$similarity_matrix_list

stopifnot(length(feature_names) == length(similarity_matrix_per_feature))
names(similarity_matrix_per_feature) <- feature_names

similarity_vector_per_feature <- lapply(similarity_matrix_per_feature, function(mat) {
  mat[upper.tri(mat, diag = FALSE)]
})

feature_names <- names(similarity_matrix_per_feature)
strat_matrix <- do.call(cbind, similarity_vector_per_feature)
colnames(strat_matrix) <- feature_names
feature_correlation_matrix <- cor(strat_matrix, method = "spearman")

## Step 2: Feature modules
dist_matrix <- as.dist(1 - feature_correlation_matrix)
hclust_features <- hclust(dist_matrix, method = "ward.D2")

feature_clusters <- cutreeDynamic(
  dendro         = hclust_features,
  distM          = as.matrix(dist_matrix),
  minClusterSize = MIN_MODULE_SIZE,
  method         = "hybrid",
  deepSplit      = DEEP_SPLIT
)
names(feature_clusters) <- feature_names

module_ids <- sort(unique(feature_clusters[feature_clusters > 0]))
n_modules <- length(module_ids)

cat(sprintf("Dynamic tree cut (deepSplit=%d) identified %d modules (+ %d unassigned)\n",
            DEEP_SPLIT, n_modules, sum(feature_clusters == 0)))
cat(sprintf("Module sizes: %s\n",
            paste(sapply(module_ids, function(m) sum(feature_clusters == m)), collapse = ", ")))

## UMAP of features (Figure 5A, Supplemental Figure S6)
library(uwot)
pca_res <- prcomp(feature_correlation_matrix, scale. = TRUE)

pcs <- pca_res$x[, 1:30]

umap_res <- umap(pcs, n_neighbors = 30, min_dist = 0.3)

# plot
clusters <- as.factor(feature_clusters)

n <- length(unique(clusters))

df <- data.frame(
  UMAP1 = umap_res[,1],
  UMAP2 = umap_res[,2],
  cluster = clusters
)

pal36 <- c(
  "#4E79A7","#A0CBE8","#F28E2B","#FFBE7D","#59A14F","#8CD17D",
  "#B6992D","#F1CE63","#499894","#86BCB6","#E15759","#FF9D9A",
  "#79706E","#BAB0AC","#D37295","#FABFD2","#B07AA1","#D4A6C8",
  "#9D7660","#D7B5A6","#2E91E5","#1CA71C","#FB0D0D","#DA16FF",
  "#222A2A","#B68100","#EB663B","#750D86","#511CFB","#00A08B",
  "#FB00D1","#FC0080","#B2828D","#6C7C32","#778AAE","#862A16"
)

cluster_ids <- sort(unique(feature_clusters))
cluster_colors <- setNames(pal36[1:length(cluster_ids)], cluster_ids)

p1 <- ggplot(df, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = 0.7) +
  scale_color_manual(values = cluster_colors, name = "Feature modules") +
  theme_classic(base_size = 16) +
  theme(
    legend.title = element_text(size=13),
    legend.key.size = unit(0.4, "cm")
  )

df$highlight <- ifelse(df$cluster == MODULE_ID, "Feature Module 27", "Others")

p2 <- ggplot(df, aes(UMAP1, UMAP2)) +
  
  geom_point(
    data = subset(df, highlight == "Others"),
    aes(color = highlight),
    size = 0.6
  ) +
  
  geom_point(
    data = subset(df, highlight == "Feature Module 27"),
    aes(color = highlight),
    size = 0.9
  ) +
  
  scale_color_manual(
    values = c(
      "Others" = "grey85",
      "Feature Module 27" = "#EB663B"
    ),
    name = "Feature Module"
  ) +
  
  theme_classic(base_size = 16)

# Supplemental Figure S6
png(file.path(FIG_DIR, "figS6.png"),width = 8, height = 5.5, units = "in", res = 600)
print(p1)
dev.off()

# Figure 5A
png(file.path(FIG_DIR, "fig5A.png"),width = 8, height = 5.5, units = "in", res = 600)
print(p2)
dev.off()

## Step 3: Helper functions

compute_module_similarity <- function(projected_list, feature_idx) {
  n_sample <- length(projected_list)
  sim_mat <- matrix(0, nrow = n_sample, ncol = n_sample)
  for (i in 1:n_sample) {
    for (j in i:n_sample) {
      x <- as.vector(projected_list[[i]][feature_idx, ])
      y <- as.vector(projected_list[[j]][feature_idx, ])
      s <- cosine(x, y)
      sim_mat[i, j] <- s
      sim_mat[j, i] <- s
    }
  }
  rownames(sim_mat) <- names(projected_list)
  colnames(sim_mat) <- names(projected_list)
  return(sim_mat)
}

# Two-group partition of the samples for one module.
find_partition_hclust <- function(dist_mat, min_group_frac) {
  n_samples <- attr(dist_mat, "Size")
  min_group_size <- max(3, ceiling(n_samples * min_group_frac))
  
  hc <- hclust(dist(dist_mat, method = "euclidean"), method = "complete")
  groups <- cutree(hc, k = 2)
  group_sizes <- table(groups)
  
  if (min(group_sizes) < min_group_size) {
    return(list(groups = NULL, silhouette = NA, hclust = hc, balanced = FALSE))
  }
  
  sil_scores <- silhouette(groups, dist_mat)
  avg_sil <- mean(sil_scores[, "sil_width"])
  
  return(list(groups = groups, silhouette = avg_sil, hclust = hc, balanced = TRUE))
}

## Steps 4-5: Score and test each module
projected_list <- no_shuffle_result$projected_combed_list
n_sample <- length(projected_list)

module_results <- data.frame(
  module        = integer(),
  n_features    = integer(),
  silhouette    = numeric(),
  perm_pval     = numeric(),
  group1_size   = integer(),
  group2_size   = integer(),
  stringsAsFactors = FALSE
)
module_assignments <- list()
module_sim_matrices <- list()
module_hclust_list <- list()   # hclust objects (row/column order of heatmaps)

cat(sprintf("\nScreening %d modules (n_samples = %d, min_group_size = %d)...\n",
            n_modules, n_sample, max(3, ceiling(n_sample * MIN_GROUP_FRAC))))

total_features <- dim(projected_list[[1]])[1]
all_feature_indices <- 1:total_features  # for permutation test

for (m in module_ids) {
  chosen_idx <- which(feature_clusters == m)
  n_feat <- length(chosen_idx)
  cat(sprintf("  Module %d: %d features ... ", m, n_feat))
  
  # Compute module-specific sample similarity
  sim_mat <- compute_module_similarity(projected_list, chosen_idx)
  dist_mat <- as.dist(1 - sim_mat)
  
  # Two-group hierarchical clustering partition
  partition <- find_partition_hclust(dist_mat, MIN_GROUP_FRAC)
  
  if (!partition$balanced) {
    cat("no balanced partition, skipping\n")
    module_results <- rbind(module_results, data.frame(
      module = m, n_features = n_feat,
      silhouette = NA, perm_pval = NA,
      group1_size = NA, group2_size = NA
    ))
    next
  }
  
  groups <- partition$groups
  sil_score <- partition$silhouette
  group_sizes <- table(groups)
  
# Permutation null: silhouette of random feature sets of the same size
  null_sils <- numeric(N_PERM)
  for (p in 1:N_PERM) {
    rand_idx <- sample(all_feature_indices, n_feat, replace = FALSE)
    rand_sim <- compute_module_similarity(projected_list, rand_idx)
    rand_dist <- as.dist(1 - rand_sim)
    rand_hc <- hclust(rand_dist, method = "ward.D2")
    rand_groups <- cutree(rand_hc, k = 2)
    rand_sil <- mean(silhouette(rand_groups, rand_dist)[, "sil_width"])
    null_sils[p] <- rand_sil
  }
  perm_pval <- (sum(null_sils >= sil_score) + 1) / (N_PERM + 1)
  
  cat(sprintf("sil=%.3f, perm_p=%.4f, groups=%d/%d\n",
              sil_score, perm_pval, group_sizes[1], group_sizes[2]))
  
  module_results <- rbind(module_results, data.frame(
    module = m, n_features = n_feat,
    silhouette = sil_score, perm_pval = perm_pval,
    group1_size = group_sizes[1], group2_size = group_sizes[2]
  ))
  
  module_assignments[[as.character(m)]] <- groups
  module_sim_matrices[[as.character(m)]] <- sim_mat
  module_hclust_list[[as.character(m)]] <- partition$hclust
}

## Step 6: Multiple testing correction and ranking
module_results$padj <- p.adjust(module_results$perm_pval, method = "BH")
module_results <- module_results[order(module_results$padj, -module_results$silhouette), ]

sig_mask <- !is.na(module_results$padj) & module_results$padj < PADJ_THRESHOLD
significant_modules <- module_results[sig_mask, ]

cat(sprintf("\n=== RESULTS ===\n"))
cat(sprintf("Total modules: %d\n", n_modules))
cat(sprintf("Balanced partitions found: %d\n", sum(!is.na(module_results$silhouette))))
cat(sprintf("Significant (padj < %.2f): %d\n", PADJ_THRESHOLD, nrow(significant_modules)))

cat("\nAll scored modules (ranked):\n")
print(module_results[!is.na(module_results$silhouette), ])

if (nrow(significant_modules) > 0) {
  cat("\nSignificant modules:\n")
  print(significant_modules)
}

## Supplemental Figure S7: heatmaps of the top modules other than MODULE_ID

library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(grid)

# Get modules with valid p-values
scored_modules <- module_results[!is.na(module_results$padj), ]
scored_modules <- scored_modules[order(scored_modules$padj, -scored_modules$silhouette), ]

cat(sprintf("Plotting %d modules with valid p-values...\n", nrow(scored_modules)))

heatmap_grobs <- list()

for (i in 1:nrow(scored_modules)) {
  m <- scored_modules$module[i]
  m_char <- as.character(m)
  
  if (!(m_char %in% names(module_sim_matrices))) next
  if (!(m_char %in% names(module_assignments))) next
  if (!(m_char %in% names(module_hclust_list))) next
  
  sim_mat <- module_sim_matrices[[m_char]]
  groups <- module_assignments[[m_char]]
  hc <- module_hclust_list[[m_char]]
  
  annot <- no_shuffle_result$annotation
  annot$Subgroup <- ifelse(groups == 1, "Group1", "Group2")
  
  annot_colors <- list(
    Subgroup = c("Group1" = "#4DBBD5", "Group2" = "#E64B35")
  )
  
  if ("Condition" %in% colnames(annot)) {
    conditions <- unique(annot$Condition)
    cond_palette <- c(
      "HIV"    = "#FC8D62",
      "HIVOUD" = "#E5C494",
      "OUD"    = "#8DA0CB",
      "CTR"    = "#66C2A5"
    )
    annot_colors$Condition <- cond_palette[conditions]
  }
  
  vals <- sim_mat[upper.tri(sim_mat)]
  breaks <- quantile(vals, probs = seq(0, 1, length.out = 101), na.rm = TRUE)
  breaks <- unique(breaks)
  
  m_row <- scored_modules[scored_modules$module == m, ]
  title_text <- sprintf(
    "Module %d | %d features | Silhouette Score =%.3f ",
    m, m_row$n_features, m_row$silhouette
  )
  
  p <- pheatmap::pheatmap(
    sim_mat,
    annotation_row = annot,
    annotation_col = annot,
    annotation_colors = annot_colors,
    breaks = breaks,
    cluster_rows = hc,
    cluster_cols = hc,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = title_text,
    fontsize = 8,
    fontsize_row = 6,
    annotation_legend = FALSE,
    silent = TRUE
  )
  
  heatmap_grobs[[m_char]] <- p$gtable
}

# Top nine modules excluding MODULE_ID (shown in Figure 5B)
other_module_ids <- setdiff(names(heatmap_grobs), as.character(MODULE_ID))
cat(sprintf("Module %d is ranked %d of %d scored modules\n", MODULE_ID,
            match(as.character(MODULE_ID), names(heatmap_grobs)), length(heatmap_grobs)))
heatmap_grobs_use <- heatmap_grobs[head(other_module_ids, 9)]

png(file.path(FIG_DIR, "figS7.png"),width =11.5, height = 12, units = "in", res = 600)

grid.arrange(
  grobs = heatmap_grobs_use,
  layout_matrix = cbind(
    matrix(1:9, nrow = 3, byrow = TRUE),
    rep(NA, 3)   # empty column on the right
  ),
  widths = c(1,1,1,0.3)
)
dev.off()
## Characterization of MODULE_ID

library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

## Module features and sample partition
chosen_idx <- which(feature_clusters == MODULE_ID)
module_features <- feature_names[chosen_idx]

# Separate RNA and ATAC features
atac_mask <- grepl("^ATAC-", module_features)
rna_features_in_module <- module_features[!atac_mask]
atac_features_in_module <- gsub("^ATAC-", "", module_features[atac_mask])

cat(sprintf("Module %d: %d total features (%d RNA, %d ATAC)\n",
            MODULE_ID, length(module_features),
            length(rna_features_in_module), length(atac_features_in_module)))

# Compute sample similarity matrix for this module
projected_list <- no_shuffle_result$projected_combed_list
n_sample <- length(projected_list)
sim_mat <- matrix(0, nrow = n_sample, ncol = n_sample)
for (i in 1:n_sample) {
  for (j in i:n_sample) {
    x <- as.vector(projected_list[[i]][chosen_idx, ])
    y <- as.vector(projected_list[[j]][chosen_idx, ])
    s <- cosine(x, y)
    sim_mat[i, j] <- s
    sim_mat[j, i] <- s
  }
}
rownames(sim_mat) <- names(projected_list)
colnames(sim_mat) <- names(projected_list)

# Hierarchical clustering into 2 groups
dist_mat <- as.dist(1 - sim_mat)
hc <- hclust(dist(dist_mat, method = "euclidean"), method = "complete")  
groups <- cutree(hc, k = 2)

# Silhouette
sil <- silhouette(groups, dist_mat)
avg_sil <- mean(sil[, "sil_width"])
cat(sprintf("Silhouette: %.3f, Group sizes: %d / %d\n",
            avg_sil, sum(groups == 1), sum(groups == 2)))

# Define group membership
group1_samples <- names(groups[groups == 1])
group2_samples <- names(groups[groups == 2])

## Heatmap (Figure 5B)
annotation <- no_shuffle_result$annotation
annotation$Subgroup <- ifelse(groups == 1, "Group1", "Group2")

annotation_colors <- list(
  Subgroup = c("Group1" = "#4DBBD5", "Group2" = "#E64B35")
)
if ("Condition" %in% colnames(annotation)) {
  conditions <- unique(annotation$Condition)
  cond_palette <- c("HIV" = "#FC8D62", "HIVOUD" = "#E5C494",
                    "OUD" = "#8DA0CB", "CTR" = "#66C2A5")
  annotation_colors$Condition <- cond_palette[intersect(conditions, names(cond_palette))]
}

p3 <- pheatmap(
  sim_mat,
  annotation_row = annotation[, "Subgroup", drop = FALSE],
  annotation_col = annotation[, "Subgroup", drop = FALSE],
  annotation_colors = annotation_colors["Subgroup"],
  cluster_rows = hc,
  cluster_cols = hc,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = sprintf("Module %d | sil=%.3f | %d/%d",
                 MODULE_ID, avg_sil, sum(groups == 1), sum(groups == 2))
)

# Cross-tabulation with condition
cat("\nSubgroup x Condition cross-tabulation:\n")
print(table(annotation$Subgroup, annotation$Condition))

# Figure 5B
png(file.path(FIG_DIR, "fig5B.png"),width = 6.2, height = 5.5, units = "in", res = 600)
print(p3)
dev.off()

## Differential analysis between subgroups

# Module 27 subgroup assignment
subgroup_map <- setNames(
  as.character(annotation$Subgroup),
  rownames(annotation)
)

# Check sample IDs
stopifnot(
  all(unique(rna_seurat_subset$sample_name) %in% names(subgroup_map)),
  all(unique(atac_seurat_subset$sample_name) %in% names(subgroup_map))
)

# Add subgroup metadata to RNA
rna_seurat_subset$subgroup <- unname(
  subgroup_map[as.character(rna_seurat_subset$sample_name)]
)

# Add subgroup metadata to ATAC
atac_seurat_subset$subgroup <- unname(
  subgroup_map[as.character(atac_seurat_subset$sample_name)]
)

# --- RNA DE ---
Idents(rna_seurat_subset) <- "subgroup"
de_rna <- FindMarkers(rna_seurat_subset,
                      ident.1 = "Group1", ident.2 = "Group2",
                      test.use = "MAST",
                      logfc.threshold = 0.25, min.pct = 0.1)
de_rna$gene <- rownames(de_rna)
de_rna$significance <- "Not Significant"
de_rna$significance[de_rna$p_val_adj < 0.01 & de_rna$avg_log2FC > 1] <- "Up in Group1"
de_rna$significance[de_rna$p_val_adj < 0.01 & de_rna$avg_log2FC < -1] <- "Down in Group1"

cat(sprintf("\nRNA DE: %d up, %d down (padj<0.01, |log2FC|>1)\n",
            sum(de_rna$significance == "Up in Group1"),
            sum(de_rna$significance == "Down in Group1")))

p_volcano_rna <- ggplot(de_rna, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up in Group1" = "#D73027",
                                "Down in Group1" = "#4575B4",
                                "Not Significant" = "grey")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_minimal() +
  labs(title = sprintf("RNA DE: Module %d Subgroups", MODULE_ID),
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-value")
print(p_volcano_rna)

# --- ATAC DA ---
Idents(atac_seurat_subset) <- "subgroup"
de_atac <- FindMarkers(atac_seurat_subset,
                       ident.1 = "Group1", ident.2 = "Group2",
                       test.use = "MAST",
                       logfc.threshold = 0.25, min.pct = 0.1)
de_atac$gene <- rownames(de_atac)
de_atac$significance <- "Not Significant"
de_atac$significance[de_atac$p_val_adj < 0.01 & de_atac$avg_log2FC > 1] <- "Up in Group1"
de_atac$significance[de_atac$p_val_adj < 0.01 & de_atac$avg_log2FC < -1] <- "Down in Group1"

cat(sprintf("ATAC DA: %d up, %d down (padj<0.01, |log2FC|>1)\n",
            sum(de_atac$significance == "Up in Group1"),
            sum(de_atac$significance == "Down in Group1")))

p_volcano_atac <- ggplot(de_atac, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up in Group1" = "#D73027",
                                "Down in Group1" = "#4575B4",
                                "Not Significant" = "grey")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_minimal() +
  labs(title = sprintf("ATAC DA: Module %d Subgroups", MODULE_ID),
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-value")
print(p_volcano_atac)

## Pathway enrichment; "Down in Group1" = up in HIV-Group2
up_genes <- de_rna$gene[de_rna$significance == "Up in Group1"]
down_genes <- de_rna$gene[de_rna$significance == "Down in Group1"]

if (length(up_genes) > 5) {
  entrez_up <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (nrow(entrez_up) > 0) {
    reactome_up <- enrichPathway(gene = entrez_up$ENTREZID, pvalueCutoff = 0.1, readable = TRUE)
    if (nrow(as.data.frame(reactome_up)) > 0) {
      p_pathway_up <- dotplot(reactome_up, showCategory = 15, font.size = 8) +
        ggtitle(sprintf("Pathways: Up in Group1 (Module %d)", MODULE_ID))
      print(p_pathway_up)
    } else {
      cat("No significant Reactome pathways for upregulated genes.\n")
    }
  }
}

if (length(down_genes) > 5) {
  entrez_down <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (nrow(entrez_down) > 0) {
    reactome_down <- enrichPathway(gene = entrez_down$ENTREZID, pvalueCutoff = 0.1, readable = TRUE)
    if (nrow(as.data.frame(reactome_down)) > 0) {
      p_pathway_down <- dotplot(reactome_down, showCategory = 15, font.size = 8) +
        ggtitle(sprintf("Pathways: Down in Group1 (Module %d)", MODULE_ID))
      print(p_pathway_down)
    } else {
      cat("No significant Reactome pathways for downregulated genes.\n")
    }
  }
}

# Figure 5E: Reactome enrichment of genes upregulated in HIV-Group2
p_E <- p_pathway_down

library(scales)
p_E <- p_E +
  labs(title = NULL) +
  scale_size(range = c(5, 11)) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  scale_color_gradient(low = "#E0605E", high = "#3B7DBF",
                       breaks = breaks_pretty(n = 3),
                       labels = label_scientific(digits = 2)) +
  guides(color = guide_colorbar(order = 1, barheight = unit(3.5, "cm")),
         size  = guide_legend(order = 2)) +
  theme(axis.text.y  = element_text(size = 17, colour = "black"),
        axis.text.x  = element_text(size = 15),
        axis.title.x = element_text(size = 17),
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14),
        legend.position = "inside",
        legend.position.inside = c(0.78, 0.30),
        legend.background = element_blank())
ggsave(file.path(FIG_DIR, "fig5E.png"), p_E, width = 5.6, height = 11, dpi = 600)

## Save objects used by figS8_pseudobulk_pca.R and fig5D_cross_region_validation.R
fig5_objects <- list(
  MODULE_ID               = MODULE_ID,
  hiv_samples             = names(projected_list),   # 18 HIV donors analysed
  samples_to_remove       = samples_to_remove,       # samples with < 10 cells
  feature_names           = feature_names,
  feature_clusters        = feature_clusters,
  module_results          = module_results,
  module_sim_matrices     = module_sim_matrices,
  module_assignments      = module_assignments,
  module_hclust_list      = module_hclust_list,
  module_silhouette       = avg_sil,                 # silhouette of the MODULE_ID partition
  groups                  = groups,                  # MODULE_ID partition (1 / 2)
  group1_samples          = group1_samples,
  group2_samples          = group2_samples,
  rna_features_in_module  = rna_features_in_module,
  atac_features_in_module = atac_features_in_module,
  de_rna                  = de_rna,
  de_atac                 = de_atac,
  up_genes                = up_genes,                # up in HIV-Group1
  down_genes              = down_genes               # up in HIV-Group2 (stress signature)
)
saveRDS(fig5_objects, file.path(RESULTS_DIR, "fig5", "fig5_session_objects.rds"))
