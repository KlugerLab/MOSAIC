# Supplemental Figure S5B: MDS of PFC donors (HIV vs CTR) on one feature module
# Input:    PFC_L23IT_RNA, PFC_L23IT_ATAC
# Output:   figures/figS5B.png

source("config.R")
source("R/functions.R")
dir.create(file.path(FIG_DIR, "other"), showWarnings = FALSE, recursive = TRUE)

library(pheatmap)
library(foreach)
library(doParallel)
library(Seurat)
library(cluster)
library(circlize)
library(vegan)
library(Signac)
library(Seurat)
library(Matrix)
library(uwot)
library(Rtsne)
library(ggplot2)
library(ggrepel)
library(viridis) 
library(gridExtra)
library(patchwork)
library(dbscan)
library(RSpectra)
library(distances)
library(FNN)
library(lsa)
library(coop)
library(dplyr)

rna_seurat_subset <- readRDS(PFC_L23IT_RNA)
atac_seurat_subset <- readRDS(PFC_L23IT_ATAC)

DefaultAssay(atac_seurat_subset) <- "GeneACT"
atac_seurat_subset[["ATAC"]] <- NULL

sample_counts <- table(rna_seurat_subset$sample_name)
samples_to_remove <- names(sample_counts[sample_counts < 10])
print(samples_to_remove)

rna_seurat_subset <- subset(rna_seurat_subset, subset = sample_name %ni% samples_to_remove)
atac_seurat_subset <- subset(atac_seurat_subset, subset = sample_name %ni% samples_to_remove)

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

# Color function
col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))
library(pheatmap)
library(ComplexHeatmap)

pca_result <- prcomp(dist_matrix, center = TRUE, scale. = TRUE)
pca_coords <- pca_result$x[, 1:15] 

mds_df <- umap(pca_coords)
mds_df <- as.data.frame(mds_df)
colnames(mds_df) <- c("Dim1", "Dim2")

ggplot(mds_df, aes(x = Dim1, y = Dim2)) +
  geom_point(color = "blue", size = 0.5) +
  labs(
    title = "MDS Plot from Correlation Matrix",
    x = "MDS Dimension 1",
    y = "MDS Dimension 2"
  ) +
  theme_minimal() +
  coord_fixed()

calculate_feature_sil_from_similarity <- function(similarity_matrix_per_feature, k) {
  
  sil_scores <- sapply(similarity_matrix_per_feature, function(sim_matrix) {
    
    # distance = 1 - similarity
    dist_matrix <- as.dist(1 - sim_matrix)
    
    # Error handling
    if (nrow(as.matrix(dist_matrix)) < k || any(is.na(dist_matrix))) {
      return(NA)
    }
    
    # k-medoids clustering
    set.seed(123)
    pam_result <- pam(dist_matrix, k = k, diss = TRUE)
    
    # mean silhouette width
    sil_result <- silhouette(pam_result$clustering, dist_matrix)
    
    return(summary(sil_result)$avg.width)
  })
  
  return(sil_scores)
}

sil_score_per_feature <- calculate_feature_sil_from_similarity(similarity_matrix_per_feature,2)

umap_df <- data.frame(
  UMAP1 = mds_df[,1],
  UMAP2 = mds_df[,2],
  cluster = as.numeric(sil_score_per_feature)
)
umap_df <- umap_df %>% 
  arrange(cluster)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.1) +
  theme_minimal() +
  labs(title = "", x = paste0("MDS_",1), y = paste0("MDS_",2)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_gradient(low = "blue", high = "red") 

threshold <- quantile(sil_score_per_feature, 0.9, na.rm = TRUE)

umap_df <- data.frame(
  UMAP1 = mds_df[,1],
  UMAP2 = mds_df[,2],
  sil_score = as.numeric(sil_score_per_feature)
)

umap_df <- umap_df %>% 
  arrange(sil_score)

# Flag the top 10% of features
umap_df <- umap_df %>%
  mutate(highlight = factor(
    ifelse(sil_score >= threshold, "Top 10%", "Other"), 
    levels = c("Other", "Top 10%")
  ))

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = highlight)) +
  geom_point(aes(size = highlight)) + 
  theme_minimal() +
  labs(
    title = "MDS Plot of Features Highlighting Top 10%",
    x = "MDS_1",
    y = "MDS_2",
    color = "Feature Group"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("Top 10%" = "red", "Other" = "grey70")) +
  scale_size_manual(values = c("Top 10%" = 0.2, "Other" = 0.1), guide = "none")

chosen_idx <- which(sil_score_per_feature >= threshold)

# Feature clusters: ward.D2, cutree k = 20; cluster 10 used for S5B
hclust_result <- hclust(dist_matrix, method = "ward.D2")
feature_clusters <- cutree(hclust_result, k = 20)

p_list <- list()
for (ing in c(10)){
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
  silhouette_scores <- silhouette(condition_numeric, 1-similarity_matrix)
  average_silhouette_score <- mean(silhouette_scores[, "sil_width"])
  print(average_silhouette_score)
  p <- plot_mds_cluster(similarity_matrix,"",condition,size=2,dis = FALSE)+ggtitle(paste0("sil score: ",average_silhouette_score))
  print(p)
  p_list[[ing]] <- p
}

gene_list <- feature_names[which(feature_clusters==10)]

plot_mds_cluster_arrow <- function(similarity_matrix, title, cluster, size=1, 
                                   custom_colors = NULL, dis=FALSE, idx1=1, idx2=2,
                                   add_ellipse=TRUE, ellipse_level=0.66,
                                   add_arrow_from_to = NULL) { 
  
  # Convert similarity to dissimilarity
  if (!dis) {
    dissimilarity_matrix <- 1 - similarity_matrix
  } else {
    dissimilarity_matrix <- similarity_matrix
  }
  
  # Perform Classical MDS
  mds_result <- cmdscale(dissimilarity_matrix, k = 5) 
  
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
  
  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster, fill = cluster)) +
    geom_point(size = size) +
    theme_minimal() +
    labs(title = title, x = paste0("MDS_", idx1), y = paste0("MDS_", idx2)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = custom_colors)
  
  # Add confidence ellipse if requested
  if (add_ellipse) {
    p <- p + 
      stat_ellipse(aes(fill = cluster), 
                   geom = "polygon", 
                   alpha = 0.2, 
                   level = ellipse_level, 
                   type = "t", 
                   linetype = 2, 
                   linewidth = 0.5) + 
      scale_fill_manual(values = custom_colors)
  }
  
  # Arrow between group centroids
  if (!is.null(add_arrow_from_to) && length(add_arrow_from_to) == 2) {
    
    start_group <- add_arrow_from_to[1]
    end_group <- add_arrow_from_to[2]
    
    if (start_group %in% umap_df$cluster && end_group %in% umap_df$cluster) {
      
      # Group centroids
      centroids_all <- aggregate(cbind(UMAP1, UMAP2) ~ cluster, data = umap_df, FUN = mean)
      
      # Subset the resulting data frame
      arrow_data_df <- centroids_all[centroids_all$cluster %in% c(start_group, end_group), ]
      
      # Check if both groups were found
      if (nrow(arrow_data_df) == 2) {
        
        # Create a data frame for the arrow segment
        arrow_data <- data.frame(
          x_start = arrow_data_df$UMAP1[arrow_data_df$cluster == start_group],
          y_start = arrow_data_df$UMAP2[arrow_data_df$cluster == start_group],
          x_end = arrow_data_df$UMAP1[arrow_data_df$cluster == end_group],
          y_end = arrow_data_df$UMAP2[arrow_data_df$cluster == end_group]
        )
        
        # Add the arrow layer to the plot
        p <- p + 
          geom_segment(
            data = arrow_data,
            aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
            arrow = arrow(type = "closed", length = unit(0.15, "inches")),
            color = "darkgrey",
            linewidth = 0.88,
            inherit.aes = FALSE  
          )
      } else {
        warning("Could not find centroids for both start and end groups.")
      }
    } else {
      warning("Arrow groups specified in 'add_arrow_from_to' not found in 'cluster' column.")
    }
  }
  
  return(p)
}

p <- plot_mds_cluster_arrow(
  similarity_matrix,
  "",
  condition,
  size = 2,
  dis = FALSE,
  add_ellipse = TRUE,
  ellipse_level=0.68,
  custom_colors = c("CTR" = "#66C2A5", "HIV" = "#FC8D62"),
  add_arrow_from_to = c("CTR", "HIV")
)

# Supplemental Figure S5B
png(file.path(FIG_DIR, "figS5B.png"),width = 6, height = 4, units = "in", res = 600)
print(p)
dev.off()

annotation_colors <- list(
  Condition = c("CTR" = "#66C2A5", "HIV" = "#FC8D62")
)

p <- pheatmap(
  similarity_matrix,
  annotation_row = no_shuffle_result$annotation,
  annotation_col = no_shuffle_result$annotation,
  annotation_colors = annotation_colors, 
  main = paste0("Heatmap of Sample Similarities",ing),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE
)

png(file.path(FIG_DIR, "other", "Subgroup_detection_validation_heatmap_cluster10_plot.png"),width = 4.7, height = 4, units = "in", res = 600)
p
dev.off()

mds_df <- umap(feature_correlation_matrix)
mds_df <- as.data.frame(mds_df)
colnames(mds_df) <- c("Dim1", "Dim2")

ggplot(mds_df, aes(x = Dim1, y = Dim2)) +
  geom_point(color = "blue", size = 0.5) +
  labs(
    title = "MDS Plot from Correlation Matrix",
    x = "MDS Dimension 1",
    y = "MDS Dimension 2"
  ) +
  theme_minimal() +
  coord_fixed()

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
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.background = element_blank(),  # No background
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  ) +
  scale_color_manual(values = c("Selected Cluster" = "#e64b35", "Other" = "grey70")) +
  scale_size_manual(values = c("Selected Cluster" = 1, "Other" = 0.1), guide = "none")

png(file.path(FIG_DIR, "other", "Subgroup_detection_validation_feature_embedding_cluster10_plot.png"),width = 6, height = 4.8, units = "in", res = 600)
p
dev.off()

library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)  

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
dotplot(reactome_results, showCategory = 15, font.size = 8) +
  theme(axis.text.y = element_text(size = 8, margin = margin(r = 5)))

go_results <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# View results
head(go_results@result)

dotplot(go_results, showCategory = 15, font.size = 8) +
  theme(axis.text.y = element_text(size = 8, margin = margin(r = 5)))

rna_seurat_subset <- AddModuleScore(rna_seurat_subset, 
                             features = list(gene_list), 
                             name = "GeneModule")

VlnPlot(rna_seurat_subset, 
        features = "GeneModule1", 
        group.by = "condition_new",
        pt.size = 0) +
  labs(title = "Gene Module Score by Condition")

p <- ggplot(rna_seurat_subset@meta.data, aes(x = condition_new, y = GeneModule1, fill = condition_new)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = c("CTR" = "#66C2A5", "HIV" = "#FC8D62")) +
  theme_minimal() +
  labs(title = "Gene Module Score by Condition",
       x = "Condition",
       y = "Module Score") +
  theme(legend.position = "none")

png(file.path(FIG_DIR, "other", "Subgroup_detection_validation_module_score.png"),width = 3, height = 4, units = "in", res = 600)
p
dev.off()

chosen_idx <- which(feature_clusters ==10)
hub_score <- rowSums(feature_correlation_matrix[chosen_idx,chosen_idx])/length(chosen_idx)
feature_hub <- feature_names[chosen_idx][order(hub_score,decreasing = TRUE)]
hub_score[order(hub_score,decreasing = TRUE)][1:10]

top_scores <- hub_score[order(hub_score, decreasing = TRUE)][1:10]
top_features <- feature_hub[1:10]
plot_data <- data.frame(
  Feature = top_features,
  Score = top_scores
)
plot_data$Feature <- factor(plot_data$Feature, levels = plot_data$Feature)

p <- ggplot(plot_data, aes(x = Feature, y = Score)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Top 10 Hub Scores",
       x = "Feature",
       y = "Hub Score")

png(file.path(FIG_DIR, "other", "Subgroup_detection_validation_hubscore.png"),width = 4, height = 4, units = "in", res = 600)
p
dev.off()

## Sub-clusters within cluster 10
chosen_feature_correlation_matrix <- feature_correlation_matrix[chosen_idx,chosen_idx]
chosen_dist_matrix <- as.dist(1 - chosen_feature_correlation_matrix)

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

hc <- hclust(chosen_dist_matrix, method = "complete")
k <- 5
clusters <- cutree(hc, k = k)
mds_result <- cmdscale(chosen_dist_matrix, k = 2)

mds_coords <- as.data.frame(mds_result)
colnames(mds_coords) <- c("Dim1", "Dim2")
mds_coords$Cluster <- as.factor(clusters)

p <- ggplot(mds_coords, aes(x = Dim1, y = Dim2, color = Cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "Feature Clusters (MDS Plot)") +
  theme_minimal()

png(file.path(FIG_DIR, "other", "temporary.png"),width = 6, height = 4.7, units = "in", res = 600)
p
dev.off()

chosen_idx_cluster <- chosen_idx[which(clusters == 5)]

n_sample <- length(no_shuffle_result$projected_combed_list)
similarity_matrix <- matrix(0, nrow = n_sample, ncol = n_sample)
for (i in 1:n_sample){
  for (j in 1:n_sample){
    x <- no_shuffle_result$projected_combed_list[[i]][chosen_idx_cluster,]
    y <- no_shuffle_result$projected_combed_list[[j]][chosen_idx_cluster,]
    similarity_matrix[i,j] <- cosine(as.vector(x),as.vector(y))
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

hc <- hclust(chosen_dist_matrix, method = "ward.D2")

plot(hc, main = "Hierarchical Clustering Dendrogram", xlab = "")

library(cluster)

pam_result <- pam(chosen_dist_matrix, k = 5, diss = TRUE)

print("K-medoid (PAM) Clusters:")
print(pam_result$clustering)

pca_result <- prcomp(chosen_dist_matrix, center = TRUE, scale. = TRUE)
pca_coords <- pca_result$x[, 1:15] 

mds_df <- umap(pca_coords)
mds_df <- as.data.frame(mds_df)
colnames(mds_df) <- c("Dim1", "Dim2")

ggplot(mds_df, aes(x = Dim1, y = Dim2)) +
  geom_point(color = "blue", size = 0.5) +
  labs(
    title = "MDS Plot from Correlation Matrix",
    x = "MDS Dimension 1",
    y = "MDS Dimension 2"
  ) +
  theme_minimal() +
  coord_fixed()

cluster_coembed <- dbscan(as.data.frame(mds_df),eps=0.08,minPts =5)
table(cluster_coembed$cluster)
