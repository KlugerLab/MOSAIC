# Figure 6B, 6C: monocyte severity prediction with abundance (pseudobulk), connectivity (MOSAIC) and integrated models
# Input:    COVID_H5AD
# Output:   results/fig6/fig6BC_objects.rds, figures/fig6B.png, figures/fig6C.png

source("config.R")
source("R/functions.R")

library(schard)
library(MASS)
library(glmnet)
library(ggplot2)
library(pROC)
library(caret)
library(ggrepel)

dir.create(file.path(RESULTS_DIR, "fig6"), showWarnings = FALSE, recursive = TRUE)

## Data preprocessing

seurat_obj <- schard::h5ad2seurat(COVID_H5AD,use.raw = TRUE)

gc()

seurat_obj_subset <- seurat_obj[,which(seurat_obj$majorType == "Mono")]
seurat_obj_subset <- seurat_obj_subset[,which(!seurat_obj_subset$`CoVID-19 severity`=="control")]

sample_counts <- table(seurat_obj_subset$sampleID)
samples_to_remove <- names(sample_counts[sample_counts < 10])
print(samples_to_remove)

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
length(genes_to_keep)

seurat_rna <- subset(seurat_obj_subset, features = genes_to_keep)
seurat_rna <- NormalizeData(seurat_rna)

## MOSAIC embedding

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
no_shuffle_result <- run_our_method_one(seurat_rna,sample_meta="sampleID",condition_meta="CoVID-19 severity",block_method = "cosine")

## Prediction

# MOSAIC input
no_shuffle_vectors <- lapply(no_shuffle_result$projected_combed_list, as.vector)
X <- do.call(rbind, no_shuffle_vectors)
label <- no_shuffle_result$annotation$Condition
y <- as.numeric(factor(label)) - 1

# Pseudobulk expression input
pseudobulk_list <- AggregateExpression(seurat_rna, 
                                       group.by = "sampleID", 
                                       assays = "RNA", 
                                       slot = "counts",
                                       normalization.method = NULL) # Use RAW counts
pseudobulk_counts <- pseudobulk_list$RNA

# Normalize pseudobulk samples (CP10K + log)
pb_obj <- CreateSeuratObject(counts = pseudobulk_counts)
pb_obj <- NormalizeData(pb_obj, normalization.method = "LogNormalize", scale.factor = 10000)
pseudobulk_norm <- GetAssayData(pb_obj, slot = "data")
X_pseudo <- t(as.matrix(pseudobulk_norm))

sample_order <- rownames(X_pseudo)
sample_meta <- unique(seurat_rna@meta.data[, c("sampleID", "CoVID-19 severity")])
rownames(sample_meta) <- sample_meta$sampleID
label_pseudo <- sample_meta[sample_order, "CoVID-19 severity"]
y_pseudo <- ifelse(label_pseudo %in% c("severe/critical"), 1, 0)

# Align samples of the two inputs
common_samples <- intersect(rownames(X), rownames(X_pseudo))

print(paste("MOSAIC Samples:", nrow(X)))
print(paste("Expression Samples:", nrow(X_pseudo)))
print(paste("Common Samples:", length(common_samples)))

# Subset and sort both matrices
X <- X[common_samples, , drop=FALSE]
X_pseudo <- X_pseudo[common_samples, , drop=FALSE]

if(!is.null(names(label))) {
  y_aligned <- label[common_samples]
} else {
  y_map <- setNames(label, rownames(no_shuffle_result$annotation)) 
  y_aligned <- y_map[common_samples]
}

y <- ifelse(y_aligned %in% c("severe/critical"), 1, 0)

# Check alignment
if(!all(rownames(X) == rownames(X_pseudo))) stop("Critical Fail: IDs still don't match!")
print("Alignment Successful.")

# Combined matrix
X_combined <- cbind(X, X_pseudo)
print(paste("Combined Dimensions:", nrow(X_combined), "Samples x", ncol(X_combined), "Features"))

X_mosaic   <- X
X_expr     <- X_pseudo
y          <- y

## CV configuration
n_repeats <- 10  # Run 5-Fold CV 10 times
n_folds   <- 5

# Storage for Statistics (AUCs per repeat)
results_mosaic <- c()
results_expr   <- c()
results_comb   <- c()

# Storage for Plotting (Concatenated Predictions across all repeats)
all_labels      <- c() 
all_preds_mos   <- c()
all_preds_exp   <- c()
all_preds_comb  <- c()

round_exp_list <- list()
round_mos_list <- list()

print("Starting Repeated Stratified CV...")

## Repeated stratified CV
for (seed in 1:n_repeats) {
  set.seed(seed)
  
  # Create Stratified Folds
  folds <- createFolds(factor(y), k = n_folds, list = FALSE)
  
  # Temp vectors for this specific round (length = N samples)
  round_mos <- rep(NA, length(y))
  round_exp <- rep(NA, length(y))
  
  for (i in 1:n_folds) {
    train_idx <- which(folds != i)
    test_idx  <- which(folds == i)
    
    # Train MOSAIC
    set.seed(seed * 1000 + i)
    cv_mos <- cv.glmnet(X_mosaic[train_idx, ], y[train_idx], 
                        family="binomial", alpha=1, standardize=TRUE)
    round_mos[test_idx] <- predict(cv_mos, X_mosaic[test_idx, ], s="lambda.min", type="response")
    
    # Train expression model
    set.seed(seed * 1000 + i)
    cv_exp <- cv.glmnet(X_expr[train_idx, ], y[train_idx], 
                        family="binomial", alpha=1, standardize=TRUE)
    round_exp[test_idx] <- predict(cv_exp, X_expr[test_idx, ], s="lambda.min", type="response")
  }
  
  # Per-round AUC
  auc_mos <- auc(roc(y, as.numeric(round_mos), quiet=TRUE))
  auc_exp <- auc(roc(y, as.numeric(round_exp), quiet=TRUE))
  
  # Ensemble (Late Fusion)
  round_comb <- (round_mos + round_exp) / 2
  auc_comb   <- auc(roc(y, as.numeric(round_comb), quiet=TRUE))
  
  # Store Statistics
  results_mosaic <- c(results_mosaic, auc_mos)
  results_expr   <- c(results_expr, auc_exp)
  results_comb   <- c(results_comb, auc_comb)
  
  # Store Predictions for Plotting (Concatenate)
  all_labels     <- c(all_labels, y)
  all_preds_mos  <- c(all_preds_mos, round_mos)
  all_preds_exp  <- c(all_preds_exp, round_exp)
  all_preds_comb <- c(all_preds_comb, round_comb)
  
  round_exp_list <- append(round_exp_list, list(round_exp))
  round_mos_list <- append(round_mos_list, list(round_mos))
  
  # Progress
  print(paste("Round", seed, "| Mosaic:", round(auc_mos,3), "| Expr:", round(auc_exp,3), "| Comb:", round(auc_comb,3)))
}

## Summary
print("------------------------------------------------")
print(paste("Mean MOSAIC AUC:", round(mean(results_mosaic), 3), "+/-", round(sd(results_mosaic),3)))
print(paste("Mean Expression AUC:", round(mean(results_expr), 3), "+/-", round(sd(results_expr),3)))
print(paste("Mean Combined AUC:", round(mean(results_comb), 3), "+/-", round(sd(results_comb),3)))

# Paired t-tests
t_test_expr <- t.test(results_comb, results_expr, paired=TRUE, alternative="greater")
pval_expr <- format.pval(t_test_expr$p.value, digits=3)
print(paste("P-value (Combined > Expression):", pval_expr))
cat(sprintf("Figure 6B p-value (one-sided paired t-test, Integrated > Pseudobulk, 10 repeat AUCs): %s\n", pval_expr))

t_test_mos <- t.test(results_comb, results_mosaic, paired=TRUE, alternative="greater")
pval_mos <- format.pval(t_test_mos$p.value, digits=3)
print(paste("P-value (Combined > MOSAIC):", pval_mos))

## Figure 6B: aggregated ROC curves
# Calculate ROC objects from the concatenated data
roc_exp_full  <- roc(all_labels, all_preds_exp)
roc_mos_full  <- roc(all_labels, all_preds_mos)
roc_comb_full <- roc(all_labels, all_preds_comb)

# ROC data
roc_data <- rbind(
  data.frame(FPR = 1 - roc_exp_full$specificities, TPR = roc_exp_full$sensitivities, Model = "Pseudobulk (Abundance)"),
  data.frame(FPR = 1 - roc_mos_full$specificities, TPR = roc_mos_full$sensitivities, Model = "MOSAIC (Connectivity)"),
  data.frame(FPR = 1 - roc_comb_full$specificities, TPR = roc_comb_full$sensitivities, Model = "Integrated (Connectivity + Abundance)")
)
roc_data$Model <- factor(roc_data$Model, levels = c("Integrated (Connectivity + Abundance)", "MOSAIC (Connectivity)", "Pseudobulk (Abundance)"))

# AUC text
stats_text <- paste0(
  "AUC Scores\n",
  "Integrated:  ", sprintf("%.3f", mean(results_comb)), "\n",
  "MOSAIC:      ", sprintf("%.3f", mean(results_mosaic)), "\n",
  "Pseudobulk:  ", sprintf("%.3f", mean(results_expr)), "\n\n",
  "Integrated vs Pseudobulk\n",
  "P < 1.09 x 10^-6"
)

p <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Model, linetype = Model, linewidth = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray80") +
  geom_line() +
  
  # Colors and styles
  scale_color_manual(values = c("#56B4E9", "#E69F00", "grey60")) +
  scale_linetype_manual(values = c("solid", "solid", "solid")) +
  scale_linewidth_manual(values = c(1.5, 1.2, 0.8)) +
  
  # Annotation, bottom right
  annotate("text", x = 1, y = 0.35, label = stats_text, 
           hjust = 1, vjust = 0, size = 2.8, color = "black", lineheight = 1.0) +
  
  labs(title = "Predictive Performance for Monocytes", x = "False Positive Rate", y = "True Positive Rate") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    # Legend, bottom right
    legend.position = c(1, 0.0), 
    legend.justification = c("right", "bottom"), 
    
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.key.width = unit(1.0, "cm"),
    legend.text = element_text(size = 10)
  )

print(p)

png(file.path(FIG_DIR, "fig6B.png"),width = 5.4, height = 4.6, units = "in", res = 600)
print(p)
dev.off()

## Figure 6C: patient risk scores

# Risk score = out-of-fold predicted probability averaged over repeats; threshold 0.5

plot_data <- data.frame(
  Expr_Risk = as.numeric(rowMeans(do.call(cbind, round_exp_list))),
  Mos_Risk  = as.numeric(rowMeans(do.call(cbind, round_mos_list))),
  Outcome   = ifelse(y == 1, "Severe", "Moderate")
)

mos_exp_cor <- cor(plot_data$Expr_Risk,plot_data$Mos_Risk)

# Categories
plot_data$Category <- "Consensus" # Default status

# Connectivity-rescued
plot_data$Category[plot_data$Expr_Risk < 0.5 & plot_data$Mos_Risk > 0.5 & plot_data$Outcome == "Severe"] <- "Type_II_Rescue"

# Abundance false positive
plot_data$Category[plot_data$Expr_Risk > 0.5 & plot_data$Mos_Risk < 0.5 & plot_data$Outcome == "Moderate"] <- "Specificity_Win"

# Abundance-rescued
plot_data$Category[plot_data$Expr_Risk > 0.5 & plot_data$Mos_Risk < 0.5 & plot_data$Outcome == "Severe"] <- "Type_I_Rescue"

# Labels
plot_data$Category <- factor(plot_data$Category, 
                             levels = c("Consensus", "Type_II_Rescue", "Type_I_Rescue", "Specificity_Win"),
                             labels = c("Consensus Severe",             
                                        "Connectivity-Rescued",         
                                        "Abundance-Rescued",            
                                        "Abundance False Positive"))       

# Colors by clinical outcome
color_map <- c("Moderate" = "#56B4E9", "Severe" = "#D55E00")

p_final <- ggplot(plot_data, aes(x = Expr_Risk, y = Mos_Risk)) +
  
  # Background zones
  annotate("rect", xmin=0, xmax=0.5, ymin=0.5, ymax=1.05, fill="#FFFFE0", alpha=0.3) +
  annotate("rect", xmin=0.5, xmax=1.05, ymin=0, ymax=0.5, fill="#E6E6FA", alpha=0.3) +
  
  # Threshold lines
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey60") +
  
  # Consensus points
  geom_point(data = subset(plot_data, Category == "Consensus Severe"), 
             aes(color = Outcome, shape = Category), size = 3, alpha = 0.6) +
  
  # Discordant points
  geom_point(data = subset(plot_data, Category != "Consensus Severe"), 
             aes(color = Outcome, shape = Category), size = 3.5, stroke = 1.2) +
  
  # Region labels
  annotate("text", x = 0.05, y = 0.95, label = "Type II: Connectivity-High", 
           hjust = 0, fontface = "bold", color = "#D55E00", size = 3.5) +
  annotate("text", x = 0.95, y = 0.10, label = "Type I: Abundance-High", 
           hjust = 1, fontface = "bold", color = "darkblue", size = 3.5) +
  
  # Scales
  scale_color_manual(values = color_map) +
  scale_shape_manual(values = c("Connectivity-Rescued" = 17,
                                "Consensus Severe" = 19,
                                "Abundance-Rescued" = 15,
                                "Abundance False Positive" = 1)) +
  # Labels and theme
  labs(
    title = "Orthogonal Mechanisms of Severity",
    subtitle = paste0("Pearson R = ",round(mos_exp_cor,2)," (Distinct Patient Stratification)"),
    x = "Abundance Risk Score (Pseudobulk)",
    y = "Connectivity Risk Score (MOSAIC)",
    color = "Clinical Outcome",
    shape = "Patient Subtype"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 11),
    
    # Legend
    legend.position = "right",
    legend.title = element_text(face="bold", size=10),
    legend.text = element_text(size=9),
    legend.key.height = unit(0.8, "cm"), 
    
    # Clean up grid
    panel.grid.minor = element_blank()
  )

print(p_final)
png(file.path(FIG_DIR, "fig6C.png"),width = 8, height = 6, units = "in", res = 600)
print(p_final)
dev.off()

## Save objects used by fig6DE_predictive_features.R
# cv_mos, cv_exp: cv.glmnet fits of the last fold of the last repeat
fig6BC_objects <- list(
  X_mosaic            = X_mosaic,
  X_expr              = X_expr,
  y                   = y,
  cv_mos              = cv_mos,
  cv_exp              = cv_exp,
  mosaic_feature_name = rownames(no_shuffle_result$projected_combed_list[[1]]),
  dim_r               = dim(no_shuffle_result$projected_combed_list[[1]])[2],
  results_mosaic      = results_mosaic,
  results_expr        = results_expr,
  results_comb        = results_comb,
  risk_scores         = plot_data
)
saveRDS(fig6BC_objects, file.path(RESULTS_DIR, "fig6", "fig6BC_objects.rds"))
