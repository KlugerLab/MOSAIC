# Figure 2A: feature embedding benchmark on simulated data (MOSAIC, PCA, MOFA+, SIMBA; k-means ARI)
# Output:   results/fig2/fig2A_ari.csv, results/fig2/fig2A_embeddings.rds, figures/fig2A_feature_benchmark.png

source("config.R")
source("R/functions.R")
source("fig2/fig2A_functions.R")   # generate_svd_simulated_data(), run_our_method_two()
source("fig2/run_simba_in_R.R")    # run_simba_embedding(), Python set-up

library(mclust)
library(irlba)
library(MOFA2)
library(matrixStats)               # rowSums2()
library(dplyr)

out_dir <- file.path(RESULTS_DIR, "fig2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
simba_workdir <- file.path(out_dir, "simba_feature_emb")

pca_feature_embedding <- function(seurat_rna, seurat_atac,d = 50,assay1="RNA",assay2="ATAC") {
  
  # Extract normalized data for the specified cells
  rna_mat <- seurat_rna@assays[[assay1]]$data
  atac_mat <- seurat_atac@assays[[assay2]]$data
  
  # Transpose to cells x features
  rna_t <- t(rna_mat)
  atac_t <- t(atac_mat)
  
  # Scale each modality separately so neither dominates
  rna_scaled <- scale(as.matrix(rna_t), center = TRUE, scale = TRUE)
  atac_scaled <- scale(as.matrix(atac_t), center = TRUE, scale = TRUE)
  
  # Handle NaN from zero-variance features
  rna_scaled[is.nan(rna_scaled)] <- 0
  atac_scaled[is.nan(atac_scaled)] <- 0
  
  # Concatenate: cells x (genes + peaks)
  combined <- cbind(rna_scaled, atac_scaled)
  
  # SVD
  svd_res <- irlba(combined, nv = d)
  
  # Feature embedding = right singular vectors (scaled by singular values)
  embedding <- svd_res$v %*% diag(svd_res$d[1:d])
  rownames(embedding) <- c(colnames(rna_t), colnames(atac_t))
  
  k <- find_elbow_kneedle(svd_res$d)
  k <- max(k, 2)
  
  list(embedding,k)
}

mosaic_avg_emb_list <- list()
mofa_emb_list <- list()
simba_emb_list <- list()
pca_emb_list <- list()
true_label_list <- list()
ari_list <- list()

for (random_seed in 42:46){
  for (noise_level_num in seq(0.25,1,0.25)){
    
    run_name <- paste0("random seed: ",random_seed," noise level: ",noise_level_num)
    
    print(paste0("Run name: ",run_name))
    
    sim <- generate_svd_simulated_data(
      simulation_type = "DC",
      add_batch = FALSE,
      whether_rescale = FALSE,
      de_prop = 0,  # no DC features
      seed = 42,
      sample_noise_level = 0.5,
      feature_noise_level = noise_level_num,
      cell_noise_level = 0.5,
      whether_nonlinear =TRUE
    )
    
    # Get ground truth clusters
    gene_clusters <- sim$seurat_X[["originalexp"]][[]]$feature_cluster
    peak_clusters <- sim$seurat_Y[["originalexp"]][[]]$feature_cluster
    genes_sim <- rownames(sim$seurat_X)
    peaks_sim <- rownames(sim$seurat_Y)
    
    # MOSAIC
    no_shuffle_result <- run_our_method_two(sim$seurat_X,sim$seurat_Y,sample_meta="sample_id",condition_meta="condition",block_method = "cosine")
    
    emb_list <- no_shuffle_result$projected_combed_list
    emb_avg <- Reduce("+", emb_list) / length(emb_list)
    
    for (j in seq_len(ncol(emb_avg))) {
      cn <- sqrt(sum(emb_avg[, j]^2))
      if (cn > 0) emb_avg[, j] <- emb_avg[, j] / cn
    }
    
    mosaic_avg_emb_list[[run_name]] <- emb_avg
    
    
    # ---- MOFA+ ----
    combined_seurat <- sim$seurat_X
    combined_seurat[["ATAC"]] <- sim$seurat_Y[["originalexp"]]
    
    DefaultAssay(combined_seurat) <- "originalexp"
    features_rna <- rownames(combined_seurat)
    DefaultAssay(combined_seurat) <- "ATAC"
    features_atac <- rownames(combined_seurat)
    
    mofa <- create_mofa(
      combined_seurat, 
      assays = c("originalexp", "ATAC"), 
      features = list(originalexp = features_rna, ATAC = features_atac)
    )
    
    model_opts <- get_default_model_options(mofa)
    mofa <- prepare_mofa(mofa,
                         model_options = model_opts
    )
    mofa <- run_mofa(mofa, use_basilisk = TRUE)
    
    W <- get_weights(mofa)
    
    mofa_r2_values <- rowSums2(mofa@cache$variance_explained$r2_per_factor$group1)
    mofa_k <- find_elbow_kneedle(mofa_r2_values)
    mofa_k <- max(mofa_k, 2)
    cat("Adaptive MOFA+ Factors Selected:", mofa_k, "\n")
    
    W_rna <- W$originalexp[,1:mofa_k]
    W_atac <- W$ATAC[,1:mofa_k]
    
    mofa_emb <- rbind(W_rna, W_atac)
    
    mofa_emb_list[[run_name]] <- mofa_emb
    
    # ---- SIMBA ----
    simba_emb <- run_simba_embedding(combined_seurat,pbg_workers = 65,
                                     workdir = simba_workdir)
    if (dir.exists(simba_workdir)) {
      unlink(simba_workdir,
             recursive = TRUE,
             force = TRUE)
    }
    simba_emb_list[[run_name]] <- simba_emb
    
    # ---- PCA ----
    pca_result <-  pca_feature_embedding(combined_seurat, combined_seurat, d = 50, assay1= "originalexp",  assay2="ATAC")
    pca_emb <- pca_result[[1]][, 1:pca_result[[2]]]
    pca_emb_list[[run_name]] <- pca_emb
    
    # ---- Benchmark ----
    true_labels <- c(gene_clusters, peak_clusters)
    true_label_list[[run_name]] <- true_labels
    
    # k-means (10 clusters) on each embedding, ARI vs true modules
    mosaic_clusters <- kmeans(emb_avg, centers = 10, nstart = 25)$cluster
    mofa_clusters   <- kmeans(mofa_emb, centers = 10, nstart = 25)$cluster
    simba_clusters <- kmeans(simba_emb, centers = 10, nstart = 25)$cluster
    simba_clusters <- simba_clusters[names(mofa_clusters)]   # same feature order as the other methods
    pca_clusters   <- kmeans(pca_emb, centers = 10, nstart = 25)$cluster
    
    ari_list[[run_name]] <- data.frame(
      seed = random_seed,
      noise = noise_level_num,
      method = c("MOSAIC", "MOFA+", "SIMBA", "PCA"),
      ARI = c(adjustedRandIndex(true_labels, mosaic_clusters),
              adjustedRandIndex(true_labels, mofa_clusters),
              adjustedRandIndex(true_labels, simba_clusters),
              adjustedRandIndex(true_labels, pca_clusters))
    )
    print(ari_list[[run_name]])
  }
}

# ---- Save results ----
saveRDS(list(mosaic = mosaic_avg_emb_list, mofa = mofa_emb_list, simba = simba_emb_list,
             pca = pca_emb_list, true_labels = true_label_list),
        file.path(out_dir, "fig2A_embeddings.rds"))

ari_df <- do.call(rbind, ari_list)
rownames(ari_df) <- NULL
write.csv(ari_df, file.path(out_dir, "fig2A_ari.csv"), row.names = FALSE)

# ---- Plot ----
# mean ARI over replicates per noise level and method
ari_summary <- ari_df %>%
  group_by(method, noise) %>%
  summarise(ARI = mean(ARI), .groups = "drop") %>%
  mutate(method = factor(method, levels = c("MOSAIC", "MOFA+", "PCA", "SIMBA")))

method_pal <- c(
  "MOSAIC" = "#1B9E77",
  "MOFA+"  = "#7570B3",
  "PCA"    = "#666666",
  "SIMBA"  = "#E6AB02"
)

fig2A_plot <- ggplot(ari_summary, aes(x = noise, y = ARI, color = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = seq(0.25, 1, 0.25)) +
  scale_color_manual(values = method_pal) +
  labs(
    x = "Noise level",
    y = "ARI",
    color = "Method"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_line(color = "grey92", linetype = "dashed", linewidth = 0.5),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  )

png(file.path(FIG_DIR, "fig2A.png"), width = 8, height = 6, units = "in", res = 600)
print(fig2A_plot)
dev.off()
