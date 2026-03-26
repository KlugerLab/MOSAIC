# Figure 5C, Supplemental Figure S9: RNA and ATAC volcano plots between the two HIV+ subgroups
# Input:    PFC_L23IT_RNA, PFC_L23IT_ATAC
# Output:   figures/fig5C.png, figures/figS9.png

source("config.R")
source("R/functions.R")
dir.create(file.path(FIG_DIR, "other"), showWarnings = FALSE, recursive = TRUE)

library(Seurat)
library(pheatmap)
library(cluster)
library(lsa)            # cosine
library(uwot)           # umap
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

rna_seurat_subset <- readRDS(PFC_L23IT_RNA)
atac_seurat_subset <- readRDS(PFC_L23IT_ATAC)

DefaultAssay(atac_seurat_subset) <- "GeneACT"
atac_seurat_subset[["ATAC"]] <- NULL

sample_counts <- table(rna_seurat_subset$sample_name)
samples_to_remove <- names(sample_counts[sample_counts < 10])
print(samples_to_remove)

rna_seurat_subset <- subset(rna_seurat_subset, subset = sample_name %ni% samples_to_remove)
atac_seurat_subset <- subset(atac_seurat_subset, subset = sample_name %ni% samples_to_remove)

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

## Sample similarity on the full embeddings

n_sample <- length(no_shuffle_result$projected_combed_list)
similarity_matrix <- matrix(0, nrow = n_sample, ncol = n_sample)
for (i in 1:n_sample){
  for (j in 1:n_sample){
    x <- no_shuffle_result$projected_combed_list[[i]]
    y <- no_shuffle_result$projected_combed_list[[j]]
    similarity_matrix[i,j] <- calculate_matrix_residual(x,y)
  }
}

rownames(similarity_matrix) <- names(no_shuffle_result$projected_combed_list)
colnames(similarity_matrix) <- names(no_shuffle_result$projected_combed_list)

p <- pheatmap(
  similarity_matrix,
  annotation_row = no_shuffle_result$annotation,
  annotation_col = no_shuffle_result$annotation,
  main = "Heatmap of Sample Similarities",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE
)
condition <- no_shuffle_result$annotation$Condition
condition_numeric <- as.integer(as.factor(condition))
silhouette_scores <- silhouette(condition_numeric, 1-similarity_matrix)
average_silhouette_score <- mean(silhouette_scores[, "sil_width"])
print(average_silhouette_score)
plot_mds_cluster(similarity_matrix,"",condition,size=2,dis = FALSE)+ggtitle(paste0("sil score: ",average_silhouette_score))

peranova_result <- cal_peranova_per_feature_parallel_no_zero(no_shuffle_result$projected_combed_list,n_sample,
                                                             groups = no_shuffle_result$annotation[names(no_shuffle_result$eigen_output_cosine_list),"Condition"],
                                                             dist_method = "cosine",num_permutation=9999,whether_permanova = FALSE)

similarity_matrix_per_feature <- peranova_result$similarity_matrix_list
similarity_vector_per_feature <- lapply(similarity_matrix_per_feature, function(mat) mat[upper.tri(mat, diag = FALSE)])

feature_correlation_matrix <- cor(do.call(cbind, similarity_vector_per_feature), method = "spearman")
dist_matrix <- as.dist(1 - feature_correlation_matrix)

library(ComplexHeatmap)
library(circlize)

# Feature clusters: ward.D2, cutree k = 30
hclust_result <- hclust(dist_matrix, method = "ward.D2")
feature_clusters <- cutree(hclust_result, k = 30)

p_list <- list()
for (ing in c(9)){
  print(ing)
  chosen_idx <- which(feature_clusters == ing)
  n_sample <- length(no_shuffle_result$projected_combed_list)
  similarity_matrix <- matrix(0, nrow = n_sample, ncol = n_sample)
  for (i in 1:n_sample){
    for (j in 1:n_sample){
      x <- no_shuffle_result$projected_combed_list[[i]][chosen_idx,]
      y <- no_shuffle_result$projected_combed_list[[j]][chosen_idx,]
      similarity_matrix[i,j] <- cosine(as.vector(x),as.vector(y))
    }
  }
  
  rownames(similarity_matrix) <- names(no_shuffle_result$projected_combed_list)
  colnames(similarity_matrix) <- names(no_shuffle_result$projected_combed_list)
  
  p <- pheatmap(
    similarity_matrix,
    annotation_row = no_shuffle_result$annotation,
    annotation_col = no_shuffle_result$annotation,
    main = paste0("Heatmap of Sample Similarities",ing),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE
  )
  
  condition <- no_shuffle_result$annotation$Condition
  condition_numeric <- as.integer(as.factor(condition))
  p <- plot_mds_cluster(similarity_matrix,"",condition,size=2,dis = FALSE)+ggtitle(paste0("sil score: ",average_silhouette_score))
  print(p)
  p_list[[ing]] <- p
}

# Group1 = HIV-Group2 of the Module 27 partition (10 donors)
group1 <- c("10288_HIV_PFC","10305_HIV_PFC","7200776574_HIV_PFC","7101847783_HIV_PFC","10280_HIV_PFC",
            "1059_HIV_PFC","1133_HIV_PFC","1068_HIV_PFC","1137_HIV_PFC","1040_HIV_PFC")

annotation <- no_shuffle_result$annotation
annotation$Group <- ifelse(rownames(annotation) %in% group1, "Group1","Group2")

annotation_colors <- list(
  Condition = c("HIV" = "#FC8D62"),
  Group = c("Group1"="#4dbbd5","Group2"="#c1e8ee")
)

p <- pheatmap(
  similarity_matrix,
  annotation_row = annotation,
  annotation_col = annotation,
  main = paste0("Heatmap of Sample Similarities",ing),
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE
)

png(file.path(FIG_DIR, "other", "Subgroup_discovery_cluster9_heatmap.png"),width = 5.5, height = 4.7, units = "in", res = 600)
p
dev.off()

# Module scores
gene_list <- feature_names[which(feature_clusters==9)]

rna_seurat_subset <- AddModuleScore(rna_seurat_subset, 
                                    features = list(gene_list), 
                                    name = "GeneModule")

atac_seurat_subset<- AddModuleScore(atac_seurat_subset, 
                                    features = list(gene_list), 
                                    name = "GeneModule")

rna_seurat_subset$subgroup <- ifelse(rna_seurat_subset$sample_name %in% group1, "Group1", "Group2")
rna_seurat_subset$subgroup <- factor(rna_seurat_subset$subgroup, levels = c("Group1", "Group2"))

atac_seurat_subset$subgroup <- ifelse(atac_seurat_subset$sample_name %in% group1, "Group1", "Group2")
atac_seurat_subset$subgroup <- factor(atac_seurat_subset$subgroup, levels = c("Group1", "Group2"))

p <- ggplot(atac_seurat_subset@meta.data, aes(x = subgroup, y = GeneModule1, fill = subgroup)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_minimal() +
  labs(title = "Gene Module Score by Condition",
       x = "Condition",
       y = "Module Score") +
  theme(legend.position = "none")

## Figure 5C: RNA DE (MAST), Group1 (HIV-Group2) vs Group2
Idents(rna_seurat_subset) <- "subgroup"

de_results <- FindMarkers(
  rna_seurat_subset,
  ident.1 = "Group1",
  ident.2 = "Group2",
  test.use = "MAST",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

de_results$gene <- rownames(de_results)

# Add significance label
de_results$significance <- "Not Significant"
de_results$significance[de_results$p_val_adj < 0.01 & de_results$avg_log2FC > 1] <- "Up-regulated in HIV-Group1"
de_results$significance[de_results$p_val_adj < 0.01 & de_results$avg_log2FC < -1] <- "Down-regulated in HIV-Group1"

# View top results
head(de_results[order(de_results$p_val_adj), ], 10)

x_range <- range(de_results$avg_log2FC, na.rm = TRUE)
y_range <- range(-log10(de_results$p_val_adj), na.rm = TRUE)

# Create volcano plot
p <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up-regulated in HIV-Group1" = "#D73027", 
                                "Down-regulated in HIV-Group1" = "#4575B4", 
                                "Not Significant" = "grey")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  theme_minimal() +
  xlim(x_range) +  # Set x limits
  labs(
    title = "Differential Expression: Block 1 vs Block 2",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significance"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )

print(p)

# Figure 5C
png(file.path(FIG_DIR, "fig5C.png"),width = 8, height = 4.7, units = "in", res = 600)
print(p)
dev.off()

## Supplemental Figure S9: ATAC gene activity DA (MAST)
Idents(atac_seurat_subset) <- "subgroup"

de_results2 <- FindMarkers(
  atac_seurat_subset,
  ident.1 = "Group1",
  ident.2 = "Group2",
  test.use = "MAST",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

de_results2$gene <- rownames(de_results2)

# Add significance label
de_results2$significance <- "Not Significant"
de_results2$significance[de_results2$p_val_adj < 0.01 & de_results2$avg_log2FC > 1] <- "Up-regulated in HIV-Group1"
de_results2$significance[de_results2$p_val_adj < 0.01 & de_results2$avg_log2FC < -1] <- "Down-regulated in HIV-Group1"

# View top results
head(de_results2[order(de_results2$p_val_adj), ], 10)

# Create volcano plot
p <- ggplot(de_results2, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up-regulated in HIV-Group1" = "#D73027", 
                                "Down-regulated in HIV-Group1" = "#4575B4", 
                                "Not Significant" = "grey")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  theme_minimal() +
  xlim(x_range) +  # Set x limits
  labs(
    title = "Differential Expression: Block 1 vs Block 2",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significance"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )

print(p)
# Supplemental Figure S9
png(file.path(FIG_DIR, "figS9.png"),width = 8, height = 4.7, units = "in", res = 600)
print(p)
dev.off()

# Reactome enrichment of genes up in Group1 (HIV-Group2)
gene_list <- rownames(de_results[which(de_results$significance %in% c("Up-regulated in HIV-Group1")),])

entrez_ids <- bitr(gene_list, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)
reactome_results <- enrichPathway(
  gene = entrez_ids$ENTREZID,
  pvalueCutoff = 0.1,
  readable = TRUE
)

# Plot results
p <- dotplot(reactome_results, showCategory = 15, font.size = 8) +
  theme(axis.text.y = element_text(size = 8, margin = margin(r = 5)))

png(file.path(FIG_DIR, "other", "Subgroup_discovery_rna_up_pathway.png"),width = 5, height = 6, units = "in", res = 600)
p
dev.off()

# Feature UMAP with cluster 9 highlighted
mds_df <- umap(feature_correlation_matrix)
mds_df <- as.data.frame(mds_df)
colnames(mds_df) <- c("Dim1", "Dim2")

chosen_idx <- which(feature_clusters==9)
mds_df <- mds_df %>%
  mutate(highlight = factor(
    ifelse(1:dim(mds_df)[1] %in% chosen_idx, "Selected Cluster", "Other"), 
    levels = c("Other", "Selected Cluster")
  ))

p <- ggplot(mds_df, aes(x = Dim1, y = Dim2, color = highlight)) +
  geom_point(aes(size = highlight)) + 
  theme_minimal() +
  labs(
    title = "UMAP Plot of Features in Selected Cluster",
    x = "UMAP_1",
    y = "UMAP_2",
    color = ""
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.05, 0.05),
    legend.justification = c("left", "bottom"),
    legend.background = element_blank(),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  ) +
  scale_color_manual(values = c("Selected Cluster" = "#e64b35", "Other" = "grey70")) +
  scale_size_manual(values = c("Selected Cluster" = 1, "Other" = 0.1), guide = "none")

png(file.path(FIG_DIR, "other", "Subgroup_discovery_feature_embedding_cluster9_plot.png"),width = 6, height = 4.8, units = "in", res = 600)
p
dev.off()
