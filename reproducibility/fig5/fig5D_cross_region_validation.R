# Figure 5D: cross-region validation of the HIV-Group2 stress signature (PFC and INS)
# Input:    PFC_COMBINED_QS2, INS_COMBINED_QS2, PFC_L23IT_RNA, results/fig5/fig5_session_objects.rds
# Output:   figures/fig5D.png
# Requires: fig5/fig5_figS6_S7_hiv_subgroup_discovery.R

source("config.R")
source("R/functions.R")

library(qs2)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)

fig5_objects <- readRDS(file.path(RESULTS_DIR, "fig5", "fig5_session_objects.rds"))
MODULE_ID              <- fig5_objects$MODULE_ID
group1_samples         <- fig5_objects$group1_samples
group2_samples         <- fig5_objects$group2_samples
down_genes             <- fig5_objects$down_genes
rna_features_in_module <- fig5_objects$rna_features_in_module

pfc <- qs_read(PFC_COMBINED_QS2, nthreads = N_CORES)
DefaultAssay(pfc) <- "RNA"
pfc[["GeneACT"]] <- NULL
pfc[["ATAC"]] <- NULL

pfc <- subset(pfc, subset = condition %in% c("HIVOUD", "OUD"))
pfc <- subset(pfc, subset = cell_types_predicted %in% c("L2/3 IT"))
pfc$sample_name <- pfc$sample_ids

sample_counts <- table(pfc$sample_name)
samples_to_remove <- names(sample_counts[sample_counts < 10])
print(samples_to_remove)
pfc <- subset(pfc, subset = sample_name %ni% samples_to_remove)

# CTR and HIV PFC L2/3 IT cells (same object as in the discovery analysis)
rna_ctr <- readRDS(PFC_L23IT_RNA)

sample_counts <- table(rna_ctr$sample_name)
samples_to_remove <- names(sample_counts[sample_counts < 10])
print(samples_to_remove)

rna_ctr <- subset(rna_ctr, subset = sample_name %ni% samples_to_remove)

shared_features <- intersect(rownames(pfc), rownames(rna_ctr))
length(shared_features)
pfc_shared <- subset(pfc, features = shared_features)
rna_ctr_shared <- subset(rna_ctr, features = shared_features)

pfc_counts <- GetAssayData(pfc, assay = "RNA", slot = "counts")[shared_features, , drop = FALSE]
ctr_counts <- GetAssayData(rna_ctr, assay = "RNA", slot = "counts")[shared_features, , drop = FALSE]

pfc_clean <- CreateSeuratObject(
  counts = pfc_counts,
  project = "PFC"
)
meta <- pfc@meta.data[colnames(pfc_clean), , drop = FALSE]
meta <- as.data.frame(meta)
pfc_clean$condition_new <- meta$condition
pfc_clean$sample_name <- meta$sample_name

ctr_clean <- CreateSeuratObject(
  counts = ctr_counts,
  meta.data = rna_ctr@meta.data,
  project = "CTR"
)

pfc_clean$dataset <- "pfc"
ctr_clean$dataset <- "rna_ctr"
combined <- merge(
  x = pfc_clean,
  y = ctr_clean,
  add.cell.ids = c("PFC", "CTR"),
  project = "combined_RNA"
)

combined[["RNA"]] <- JoinLayers(combined[["RNA"]])
combined <- NormalizeData(combined)

combined$subgroup <- combined$condition_new
idx_group1 <- which(combined$sample_name %in% group1_samples)
combined$subgroup[idx_group1] <- "HIV_Group1"
idx_group2 <- which(combined$sample_name %in% group2_samples)
combined$subgroup[idx_group2] <- "HIV_Group2"

if (length(rna_features_in_module) > 0) {
  rna_genes_valid <- intersect(rna_features_in_module, rownames(combined))
  if (length(rna_genes_valid) > 0) {
    combined <- AddModuleScore(combined,
                                        features = list(down_genes),
                                        name = "RNA_Module")
    combined$subgroup <- factor(
      combined$subgroup,
      levels = c("CTR", "OUD", "HIVOUD","HIV_Group2", "HIV_Group1")
    )
    p_rna_score <- ggplot(combined@meta.data,
                          aes(x = subgroup, y = RNA_Module1, fill = subgroup)) +
      geom_boxplot(outlier.size = 0.5) +
      theme_minimal() +
      labs(title = sprintf("RNA Module Score (Module %d)", MODULE_ID),
           x = "Subgroup", y = "Module Score") +
      theme(legend.position = "none")
    print(p_rna_score)
  }
}

# Pseudobulk module score: average per sample
meta <- combined@meta.data
pseudobulk_scores <- meta %>%
  group_by(sample_name, subgroup) %>%
  summarise(mean_module_score = mean(RNA_Module1, na.rm = TRUE), .groups = "drop")

p <- ggplot(pseudobulk_scores, aes(x = subgroup, y = mean_module_score, fill = subgroup)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("CTR" = "#66C2A5", "HIV_Group2" = "#E64B35",
                               "HIV_Group1" = "#4DBBD5", "HIVOUD" = "#FC8D62",
                               "OUD" = "#E5C494")) +
  theme_minimal() +
  labs(title = sprintf("Pseudobulk Module Score (Module %d)", MODULE_ID),
       x = "Subgroup", y = "Mean Module Score per Sample") +
  theme(legend.position = "none")
print(p)

png(file.path(FIG_DIR, "pfc_module_score.png"),width = 6.2, height = 5, units = "in", res = 600)
print(p)
dev.off()

ins <- qs_read(INS_COMBINED_QS2, nthreads = N_CORES)
DefaultAssay(ins) <- "RNA"
ins[["GeneACT"]] <- NULL
ins[["ATAC"]] <- NULL

ins <- subset(ins, subset = cell_types_predicted %in% c("L2/3 IT"))

sample_counts <- table(ins$sample_ids)
samples_to_remove <- names(sample_counts[sample_counts < 10])
print(samples_to_remove)
ins <- NormalizeData(ins)

library(stringr)
ins$sample_name <- sub("_.*", "", ins$sample_ids)

group1_ids <- sub("_.*", "", group1_samples)
group2_ids <- sub("_.*", "", group2_samples)

ins$subgroup <- ins$condition
idx_group1 <- which(ins$sample_name %in% group1_ids)
ins$subgroup[idx_group1] <- "HIV_Group1"
idx_group2 <- which(ins$sample_name %in% group2_ids)
ins$subgroup[idx_group2] <- "HIV_Group2"

ins <- subset(ins, subset = subgroup %ni% c("HIV"))

ins$subgroup <- factor(
  ins$subgroup,
  levels = c("CTR", "OUD", "HIVOUD","HIV_Group2", "HIV_Group1")
)

if (length(rna_features_in_module) > 0) {
  rna_genes_valid <- intersect(rna_features_in_module, rownames(ins))
  if (length(rna_genes_valid) > 0) {
    ins <- AddModuleScore(ins,
                               features = list(down_genes),
                               name = "RNA_Module")
    
    p_rna_score <- ggplot(ins@meta.data,
                          aes(x = subgroup, y = RNA_Module1, fill = subgroup)) +
      geom_boxplot(outlier.size = 0.5) +
      theme_minimal() +
      labs(title = sprintf("RNA Module Score (Module %d)", MODULE_ID),
           x = "Subgroup", y = "Module Score") +
      theme(legend.position = "none")
    print(p_rna_score)
  }
}

# Pseudobulk module score: average per sample
meta <- ins@meta.data
pseudobulk_scores <- meta %>%
  group_by(sample_name, subgroup) %>%
  summarise(mean_module_score = mean(RNA_Module1, na.rm = TRUE), .groups = "drop")

p <- ggplot(pseudobulk_scores, aes(x = subgroup, y = mean_module_score, fill = subgroup)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("CTR" = "#FC8D62", "HIV_Group2" = "#B8A500",
                               "HIV_Group1" = "#00BA38", "HIVOUD" = "#619CFF",
                               "OUD" = "#F564E3")) +
  theme_minimal() +
  labs(title = sprintf("Pseudobulk Module Score (Module %d)", MODULE_ID),
       x = "Subgroup", y = "Mean Module Score per Sample") +
  theme(legend.position = "none")
print(p)

pfc_hiv <- subset(combined, subset = condition_new %in% c("HIV"))
ins_hiv <- subset(ins, subset = subgroup %in% c("HIV_Group2","HIV_Group1"))

if (length(rna_features_in_module) > 0) {
  rna_genes_valid <- intersect(rna_features_in_module, rownames(pfc_hiv))
  if (length(rna_genes_valid) > 0) {
    pfc_hiv <- AddModuleScore(pfc_hiv,
                               features = list(down_genes),
                               name = "RNA_Module")
    pfc_hiv$subgroup <- factor(
      pfc_hiv$subgroup,
      levels = c("HIV_Group2", "HIV_Group1")
    )
    p_rna_score <- ggplot(pfc_hiv@meta.data,
                          aes(x = subgroup, y = RNA_Module1, fill = subgroup)) +
      geom_boxplot(outlier.size = 0.5) +
      theme_minimal() +
      labs(title = sprintf("RNA Module Score (Module %d)", MODULE_ID),
           x = "Subgroup", y = "Module Score") +
      theme(legend.position = "none")
    print(p_rna_score)
  }
}

if (length(rna_features_in_module) > 0) {
  rna_genes_valid <- intersect(rna_features_in_module, rownames(combined))
  if (length(rna_genes_valid) > 0) {
    ins_hiv <- AddModuleScore(ins_hiv,
                               features = list(down_genes),
                               name = "RNA_Module")
    ins_hiv$subgroup <- factor(
      ins_hiv$subgroup,
      levels = c("HIV_Group2", "HIV_Group1")
    )
    p_rna_score <- ggplot(ins_hiv@meta.data,
                          aes(x = subgroup, y = RNA_Module1, fill = subgroup)) +
      geom_boxplot(outlier.size = 0.5) +
      theme_minimal() +
      labs(title = sprintf("RNA Module Score (Module %d)", MODULE_ID),
           x = "Subgroup", y = "Module Score") +
      theme(legend.position = "none")
    print(p_rna_score)
  }
}

## Figure 5D
ins_hiv <- AddModuleScore(
  ins_hiv,
  features = list(down_genes),
  name = "RNA_Module"
)

ins_hiv$subgroup <- factor(
  ins_hiv$subgroup,
  levels = c("HIV_Group1", "HIV_Group2")
)

pfc_hiv <- AddModuleScore(
  pfc_hiv,
  features = list(down_genes),
  name = "RNA_Module"
)

pfc_hiv$subgroup <- factor(
  pfc_hiv$subgroup,
  levels = c("HIV_Group1", "HIV_Group2")
)

# Compute shared y-axis limits across both datasets
pfc_scores <- pfc_hiv@meta.data$RNA_Module1[!is.na(pfc_hiv@meta.data$subgroup) & !is.na(pfc_hiv@meta.data$RNA_Module1)]
ins_scores <- ins_hiv@meta.data$RNA_Module1[!is.na(ins_hiv@meta.data$subgroup) & !is.na(ins_hiv@meta.data$RNA_Module1)]
global_ymin <- min(c(pfc_scores, ins_scores))
global_ymax <- max(c(pfc_scores, ins_scores))
global_yrange <- global_ymax - global_ymin

plot_module_violin <- function(seurat_obj, title_text, global_ymin, global_ymax, global_yrange) {
  plot_df <- seurat_obj@meta.data
  plot_df <- plot_df[!is.na(plot_df$subgroup) & !is.na(plot_df$RNA_Module1), ]
  
  ggplot(plot_df, aes(x = subgroup, y = RNA_Module1, fill = subgroup)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black") +
    scale_fill_manual(values = c("HIV_Group2" = "#E64B35", "HIV_Group1" = "#4DBBD5")) +
    # Wilcoxon rank-sum test on per-cell module scores
    stat_compare_means(
      method = "wilcox.test",
      comparisons = list(c("HIV_Group2", "HIV_Group1")),
      label = "p.signif",
      size = 5,
      label.y = global_ymax + 0.05 * global_yrange
    ) +
    coord_cartesian(ylim = c(global_ymin, global_ymax + 0.15 * global_yrange)) +
    labs(title = title_text, x = "Subgroup", y = "Module Score") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "grey", fill = NA, linewidth = 1)
    )
}

p_pfc <- plot_module_violin(pfc_hiv, "Prefrontal Cortex", global_ymin, global_ymax, global_yrange)
p_ins <- plot_module_violin(ins_hiv, "Insular Cortex", global_ymin, global_ymax, global_yrange)

# Figure 5D
png(file.path(FIG_DIR, "fig5D.png"),width = 6.2, height = 5.5, units = "in", res = 600)
print(p_pfc + p_ins)
dev.off()
