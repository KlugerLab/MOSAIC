# Figure 6A: lineage screen, MOSAIC severity prediction (Lasso, repeated stratified 5-fold CV) in five immune lineages
# Input:    COVID_H5AD
# Output:   results/fig6/fig6A_result_list.RData, figures/fig6A.png

source("config.R")
source("R/functions.R")

library(schard)
library(MASS)
library(glmnet)
library(pROC)
library(caret)

dir.create(file.path(RESULTS_DIR, "fig6"), showWarnings = FALSE, recursive = TRUE)

run_our_method_one <- function(rna_seurat, assay1 = "RNA",sample_meta="sample_id",condition_meta="condition",block_method="cosine"){
  # Global scale
  rna_seurat <- ScaleData(rna_seurat,features=rownames(rna_seurat))
  
  rna.list <- SplitObject(rna_seurat, split.by = sample_meta)
  
  dataset_list <- list()
  for (i in 1:length(rna.list)){
    print(i)
    gene_seurat <- rna.list[[i]]
    gene_matrix <- gene_seurat[[assay1]]$scale.data
    gene_matrix <- t(as.matrix(gene_matrix))
    
    dataset_list[[i]] <- list(gene_matrix)
  }
  
  eigen_output_cosine_list <- list()
  r_cosine_list <- c()
  block_matrix_list <- list()
  for (i in 1:length(dataset_list)){
    print(i)
    matrix_set <- dataset_list[[i]]
    
    if (block_method == "cosine"){
      print("cosine")
      eigen_output_cosine <- eigen_block_cosine_one(matrix_set[[1]])
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

## Data preprocessing

seurat_obj <- schard::h5ad2seurat(COVID_H5AD)

celltypes <- c("Mono","B","CD4","CD8","NK")
result_list <- list()

for (cell_type in celltypes){
  print(cell_type)
  
  seurat_obj_subset <- seurat_obj[,which(seurat_obj$majorType == cell_type)]
  seurat_obj_subset <- seurat_obj_subset[,which(!seurat_obj_subset$`CoVID-19 severity`=="control")]
  
  sample_counts <- table(seurat_obj_subset$sampleID)
  samples_to_remove <- names(sample_counts[sample_counts < 10])
  print(paste0("Sample to remove: ", samples_to_remove))
  
  seurat_obj_subset <- subset(seurat_obj_subset, subset = sampleID %ni% samples_to_remove)
  # Feature selection
  counts <- GetAssayData(seurat_obj_subset, assay = "RNA", slot = "counts")
  sample_ids <- unique(as.character(seurat_obj_subset$sampleID))
  list_of_expressed_genes <- lapply(sample_ids, function(id) {
    cells_in_sample <- which(seurat_obj_subset$sampleID == id)
    sample_counts <- counts[, cells_in_sample, drop = FALSE]
    expressed_genes <- rownames(sample_counts)[rowSums(sample_counts) > 3]
    return(expressed_genes)
  })
  genes_to_keep <- Reduce(intersect, list_of_expressed_genes)
  
  print(paste0("Number of genes: ", length(genes_to_keep)))
  seurat_obj_subset[["RNA"]] <- subset(seurat_obj_subset[["RNA"]], features = genes_to_keep)
  seurat_rna <- NormalizeData(seurat_obj_subset)
  
  print("finish preprocessing data")
  
  ## MOSAIC embedding
  suppressMessages(
    no_shuffle_result <- run_our_method_one(
      seurat_rna,
      sample_meta = "sampleID",
      condition_meta = "CoVID-19 severity",
      block_method = "cosine"
    )
  )
  
  ## Prediction
  
  # Patient-level connectivity profiles (vectorized feature embeddings)
  no_shuffle_vectors <- lapply(no_shuffle_result$projected_combed_list, as.vector)
  X <- do.call(rbind, no_shuffle_vectors)
  label <- no_shuffle_result$annotation$Condition
  y <- as.numeric(factor(label)) - 1
  
  # Lasso, repeated stratified CV
  n_repeats <- 10  # Run 5-Fold CV 10 times
  n_folds <- 5
  
  result_list[[cell_type]] <- rep(0,n_repeats)
  
  for (seed in 1:(n_repeats)) {
    set.seed(seed)
    
    # Stratified folds
    folds <- createFolds(factor(y), k = n_folds, list = FALSE)
    
    # Temporary storage for this round
    preds_mos_round <- rep(NA, length(y))
    
    for (i in 1:n_folds) {
      train_idx <- which(folds != i)
      test_idx  <- which(folds == i)
      
      # Lasso
      set.seed(seed * 1000 + i)
      cv_mos <- cv.glmnet(X[train_idx, ], y[train_idx], 
                          family="binomial", alpha=1, standardize=TRUE)
      preds_mos_round[test_idx] <- predict(cv_mos, X[test_idx, ], s="lambda.min", type="response")
    }
    auc_mos <- auc(roc(y, as.numeric(preds_mos_round), quiet=TRUE))
      
    result_list[[cell_type]][seed] <- auc_mos
  }
      
    print(paste("Stratified Lasso 5-Fold CV AUC-", cell_type,":",mean(result_list[[cell_type]])))
  
}

save(result_list,file=file.path(RESULTS_DIR, "fig6", "fig6A_result_list.RData"))

## Plot

library(ggplot2)
library(dplyr)

# Mean and SD of AUC per lineage

cell_types <- names(result_list)
plot_data <- data.frame(
  CellType = cell_types,
  Mean_AUC = sapply(result_list, mean, na.rm=TRUE),
  SD_AUC   = sapply(result_list, sd, na.rm=TRUE)
)

# Sort by AUC
plot_data <- plot_data %>% arrange(desc(Mean_AUC))

# Lock the order for plotting
plot_data$CellType <- factor(plot_data$CellType, levels = plot_data$CellType)

# Create a "Highlight" column to color Monocytes differently
plot_data$Highlight <- ifelse(plot_data$CellType == "Mono", "Driver", "Other")

# Print to check numbers
print(plot_data)

# Display names
display_names <- c(Mono = "Monocyte", CD4 = "CD4+ T", CD8 = "CD8+ T", B = "B", NK = "NK")
levels(plot_data$CellType) <- unname(display_names[levels(plot_data$CellType)])

## Figure 6A
p_6A <- ggplot(plot_data, aes(x = CellType, y = Mean_AUC, fill = Highlight)) +
  # Bar Plot
  geom_bar(stat = "identity", width = 0.7, color = "black", alpha = 0.8) +
  
  # Error bars (mean +/- SD)
  geom_errorbar(aes(ymin = Mean_AUC - SD_AUC, ymax = Mean_AUC + SD_AUC), 
                width = 0.2, position = position_dodge(0.9)) +
  
  # Monocytes highlighted
  scale_fill_manual(values = c("Driver" = "salmon2", "Other" = "#B09C8599")) +
  
  geom_text(aes(label = sprintf("%.3f", Mean_AUC),  # Format to 3 decimals
                y = Mean_AUC + SD_AUC + 0.02),      # Place just above error bar
            vjust = 0, size = 4) +
  
  # Y-axis range
  coord_cartesian(ylim = c(0.5, 1.0)) +
  scale_y_continuous(breaks = seq(0.5, 1.0, 0.1)) +
  
  # Aesthetics
  theme_classic(base_size = 14) +
  labs(
    title = "     Unbiased Screen of Immune Drivers",
    y = "Prediction Accuracy (AUC)",
    x = ""
  ) +
  
  # No legend
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(size = 14,hjust = 0.5))

# Display
print(p_6A)

png(file.path(FIG_DIR, "fig6A.png"),width = 5, height = 5.6, units = "in", res = 600)
print(p_6A)
dev.off()
