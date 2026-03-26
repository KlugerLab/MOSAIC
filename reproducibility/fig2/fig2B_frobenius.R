# Figure 2B (left): pairwise Frobenius distances between per-sample correlation matrices, real vs shuffled null
# Input:    PFC_SEURAT_RDS
# Output:   results/fig2/panel_b_frobenius_results.RData, figures/fig2B_left.png

source("config.R")
source("R/functions.R")   # get_genes_expressed_across_samples()

library(Signac)
library(Seurat)
library(ggplot2)
library(ggsignif)
library(ArchR)
library(parallel)
set.seed(1234)

addArchRThreads(N_CORES)
addArchRGenome("hg38")

# ---- Config ----
CELL_TYPE <- "L2/3 IT"   # same population as Figure 5
MIN_CELLS <- 50          # minimum cells per sample
N_PERM    <- 10
SET_SEED  <- 42
OUT_DIR   <- file.path(RESULTS_DIR, "fig2")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Preprocessing ----
seurat_obj <- readRDS(PFC_SEURAT_RDS)
seurat_obj <- seurat_obj[, which(seurat_obj$consensus_annotation_v3 == CELL_TYPE)]

# drop samples with too few cells before feature selection
cells_per_sample <- table(seurat_obj$sample_name)
cat("Cells per sample in", CELL_TYPE, ":\n"); print(sort(cells_per_sample))
keep_samples <- names(cells_per_sample)[cells_per_sample >= MIN_CELLS]
cat("Samples kept:", length(keep_samples), "of", length(cells_per_sample),
    "(dropped:", paste(setdiff(names(cells_per_sample), keep_samples), collapse = ", "), ")\n")
seurat_obj <- seurat_obj[, seurat_obj$sample_name %in% keep_samples]

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 3000)
variable_genes <- VariableFeatures(seurat_obj)
expressed_gene_across_sample <- get_genes_expressed_across_samples(seurat_obj, "sample_name", min_expr = 10)
gene_to_keep <- variable_genes[variable_genes %in% expressed_gene_across_sample]

DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj <- FindTopFeatures(seurat_obj, min = 'q99', verbose = TRUE)
atac_to_keep <- VariableFeatures(seurat_obj)

DefaultAssay(seurat_obj) <- "RNA"
seurat_rna <- seurat_obj[gene_to_keep, ]

DefaultAssay(seurat_obj) <- "ATAC"
seurat_atac <- seurat_obj
seurat_atac[["RNA"]] <- NULL
seurat_atac <- seurat_atac[atac_to_keep, ]

seurat_rna  <- NormalizeData(seurat_rna, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_atac <- RunTFIDF(seurat_atac)

# ---- Extract data ----
rna_mat  <- seurat_rna@assays$RNA$data
atac_mat <- seurat_atac@assays$ATAC$data

shared_cells <- intersect(colnames(rna_mat), colnames(atac_mat))
rna_mat  <- rna_mat[, shared_cells]
atac_mat <- atac_mat[, shared_cells]

sample_labels <- seurat_rna@meta.data[shared_cells, "sample_name"]
sample_names  <- unique(sample_labels)
cat("Cells:", length(shared_cells), " Samples:", length(sample_names), "\n")

combined_dense <- as.matrix(t(rbind(rna_mat, atac_mat)))   # cells x features
cat("Combined matrix:", nrow(combined_dense), "cells x", ncol(combined_dense), "features\n")
rm(rna_mat, atac_mat); gc()

# ---- Helpers ----
compute_sample_cormat <- function(cell_idx, combined_dense) {
  X <- combined_dense[cell_idx, , drop = FALSE]
  cor_mat <- suppressWarnings(cor(X, method = "pearson"))
  cor_mat[is.na(cor_mat)] <- 0
  cor_mat
}

compute_pairwise_frobenius_sub <- function(cor_mats, row_idx, col_idx) {
  n <- length(cor_mats)
  sub_mats <- lapply(cor_mats, function(cm) cm[row_idx, col_idx, drop = FALSE])
  pairs <- combn(n, 2)
  vapply(seq_len(ncol(pairs)), function(k)
    sqrt(sum((sub_mats[[pairs[1, k]]] - sub_mats[[pairs[2, k]]])^2)), numeric(1))
}

# ---- Real and null cell indices (null: shuffled cell-to-sample labels) ----
n_samples <- length(sample_names)
real_cell_indices <- lapply(sample_names, function(s) which(sample_labels == s))

set.seed(SET_SEED)
null_cell_indices <- vector("list", N_PERM * n_samples)
for (pid in 1:N_PERM) {
  shuffled <- sample(sample_labels)
  for (sid in 1:n_samples) {
    null_cell_indices[[(pid - 1) * n_samples + sid]] <- which(shuffled == sample_names[sid])
  }
}
all_cell_indices <- c(real_cell_indices, null_cell_indices)
n_real <- n_samples

# ---- Correlation matrices ----
t0 <- proc.time()
all_cor_mats <- mclapply(all_cell_indices, compute_sample_cormat,
                         combined_dense = combined_dense, mc.cores = N_CORES)
cat("Correlation matrices done in", round((proc.time() - t0)[3] / 60, 1), "minutes.\n")

real_cor_mats     <- all_cor_mats[1:n_real]
null_cor_mats_all <- all_cor_mats[(n_real + 1):length(all_cor_mats)]

# ---- Distances by modality ----
feature_names <- colnames(combined_dense)
is_peak  <- grepl("^chr", feature_names)
gene_idx <- which(!is_peak)
peak_idx <- which(is_peak)
cat("Genes:", length(gene_idx), " Peaks:", length(peak_idx), "\n")

modality_pairs <- list(
  "Gene – Gene" = list(row = gene_idx, col = gene_idx),
  "Gene – Peak" = list(row = gene_idx, col = peak_idx),
  "Peak – Peak" = list(row = peak_idx, col = peak_idx)
)

results <- list()
for (mod_name in names(modality_pairs)) {
  ri <- modality_pairs[[mod_name]]$row
  ci <- modality_pairs[[mod_name]]$col
  real_d <- compute_pairwise_frobenius_sub(real_cor_mats, ri, ci)
  null_d <- unlist(lapply(1:N_PERM, function(pid) {
    idx <- ((pid - 1) * n_samples + 1):(pid * n_samples)
    compute_pairwise_frobenius_sub(null_cor_mats_all[idx], ri, ci)
  }))
  test <- wilcox.test(real_d, null_d, alternative = "greater")
  # effect size
  ratio <- median(real_d) / median(null_d)
  cat(sprintf("%s: real median %.2f, null median %.2f, ratio %.2f, P = %.2e\n",
              mod_name, median(real_d), median(null_d), ratio, test$p.value))
  results[[mod_name]] <- list(real = real_d, null = null_d, p = test$p.value, ratio = ratio)
}

# ---- Plot ----
plot_df <- do.call(rbind, lapply(names(results), function(m) rbind(
  data.frame(modality = m, group = "Real Samples",  distance = results[[m]]$real),
  data.frame(modality = m, group = "Shuffled Null", distance = results[[m]]$null))))
plot_df$group    <- factor(plot_df$group, levels = c("Real Samples", "Shuffled Null"))
plot_df$modality <- factor(plot_df$modality, levels = names(modality_pairs))

# bracket height from the upper whiskers
whisker_top <- sapply(names(modality_pairs), function(m)
  max(sapply(c("Real Samples", "Shuffled Null"), function(g)
    boxplot.stats(plot_df$distance[plot_df$modality == m & plot_df$group == g])$stats[5])))
y_max <- max(whisker_top) * 1.15

# P values typeset as 1.4 x 10^-95 (plotmath)
fmt_p <- function(p) {
  e <- floor(log10(p)); m <- signif(p / 10^e, 2)
  sprintf("italic(P) == %s %%*%% 10^%d", m, e)
}
p_labels <- sapply(names(modality_pairs), function(m) fmt_p(results[[m]]$p))

p <- ggplot(plot_df, aes(x = modality, y = distance, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, size = 0.3, fatten = 1,
               position = position_dodge(0.7)) +
  scale_fill_manual(values = c("Real Samples" = "#2166ac", "Shuffled Null" = "#d1d1d1"), name = NULL) +
  geom_signif(y_position = whisker_top * 1.05,
              xmin = c(0.825, 1.825, 2.825), xmax = c(1.175, 2.175, 3.175),
              annotations = p_labels, parse = TRUE,
              tip_length = 0.01, textsize = 1.8, size = 0.3, vjust = 0) +
  coord_cartesian(ylim = c(0, y_max)) +
  labs(x = NULL, y = "Pairwise Frobenius distance") +
  theme_classic(base_size = 7) +
  theme(text = element_text(family = "Arial"),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        axis.text = element_text(color = "black"),
        legend.position = "top",
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"))

png(file.path(FIG_DIR, "fig2B_left.png"), width = 2.5, height = 2.5, units = "in", res = 600)
print(p)
dev.off()

# ---- Save (input for fig2/fig2B_pairs.R) ----
# per-sample detection fraction of each feature
det_frac <- sapply(real_cell_indices, function(idx) colMeans(combined_dense[idx, , drop = FALSE] > 0))
colnames(det_frac) <- sample_names
save(results, real_cor_mats, sample_names, feature_names, det_frac,
     cells_per_sample, keep_samples, CELL_TYPE, MIN_CELLS,
     file = file.path(OUT_DIR, "panel_b_frobenius_results.RData"))
cat("Done.\n")
