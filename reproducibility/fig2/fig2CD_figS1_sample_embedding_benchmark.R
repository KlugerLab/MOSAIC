# Figure 2C, 2D, S1: cross-sample embedding comparability (MOSAIC vs MOFA+, SIMBA, PCA with Procrustes)
# Input:    PFC_SEURAT_RDS
# Output:   results/fig2/sample_embed_benchmark/*.rds, fig2C_silhouette.csv, figures/fig2C.png, fig2D.png, figS1.png, figures/fig2D.png, figures/figS2_PCoA.png

source("config.R")
source("R/functions.R")            # find_elbow_kneedle(), get_genes_expressed_across_samples(), eigen_block_cosine_two(), ...
source("fig2/fig2A_functions.R")   # run_our_method_two() (MOSAIC)
source("fig2/run_simba_in_R.R")    # run_simba_embedding(), Python set-up (MOFA_SIMBA_PYTHON)

library(Signac)
library(Seurat)
library(ArchR)
library(ggplot2)
library(patchwork)
library(stringr)
library(cluster)
library(reticulate)
library(MOFA2)
library(matrixStats)
library(future)
library(dplyr)
library(irlba)
library(Matrix)
library(uwot)      # umap() in plot_umap_scorch()
set.seed(1234)

addArchRThreads(threads = N_CORES)
addArchRGenome("hg38")

# Preload the CUDA libraries shipped with the Python environment (GPU MOFA+)
if (!is.na(MOFA_SIMBA_NVIDIA_DIR)) {
  py_run_string(sprintf("
import ctypes, glob, os
base = '%s'
for lib_dir in glob.glob(os.path.join(base, '*/lib')):
    for so in glob.glob(os.path.join(lib_dir, '*.so*')):
        try:
            ctypes.CDLL(so)
        except:
            pass
print('CUDA libs loaded')
", MOFA_SIMBA_NVIDIA_DIR))
}

simba_workdir <- file.path(RESULTS_DIR, "fig2", "simba_subsample_disease")

# ---- Preprocessing ----

seurat_obj <- readRDS(PFC_SEURAT_RDS)

# 6 control donors with the highest cell counts
sample_to_keep <- c("HCCRS_CTR_PFC","HCTXJ_CTR_PFC","HCTMW_CTR_PFC","HCTPK_CTR_PFC","HCTME_CTR_PFC","HCTNZ_CTR_PFC")

seurat_obj_filtered <- seurat_obj[,which(seurat_obj$sample_name %in% sample_to_keep)]

DefaultAssay(seurat_obj_filtered) <- "RNA"
seurat_obj_filtered <- NormalizeData(seurat_obj_filtered, normalization.method = "LogNormalize", scale.factor = 10000)

seurat_obj_filtered <- FindVariableFeatures(seurat_obj_filtered,nfeatures=10000)
variable_genes <- VariableFeatures(seurat_obj_filtered)
expressed_gene_across_sample <- get_genes_expressed_across_samples(seurat_obj_filtered,"sample_name",min_expr=10)

gene_to_keep <- variable_genes[variable_genes %in% expressed_gene_across_sample]

DefaultAssay(seurat_obj_filtered) <- "ATAC"
seurat_obj_filtered <- FindTopFeatures(seurat_obj_filtered, min='q99',verbose = TRUE)
atac_to_keep <- VariableFeatures(seurat_obj_filtered)

DefaultAssay(seurat_obj_filtered) <- "RNA"
seurat_rna <- seurat_obj_filtered[gene_to_keep,]

DefaultAssay(seurat_obj_filtered) <- "ATAC"
seurat_atac <- seurat_obj_filtered
seurat_atac[["RNA"]] <- NULL
seurat_atac <- seurat_atac[atac_to_keep,]

seurat_rna <- NormalizeData(seurat_rna, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_atac <- RunTFIDF(seurat_atac)

# ---- Functions ----
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

run_pca_subsamples <- function(seurat_obj,split_by="sample_name_split",d = 50,assay1="RNA",assay2="ATAC"){
  svd_emb_list <- list()
  k_list <- c()
  
  seurat_list <- SplitObject(seurat_obj, split.by = split_by)
  for (sample_to_use in names(seurat_list)){
    print(length(svd_emb_list)+1)
    seurat_subsample <- seurat_list[[sample_to_use]]
    
    DefaultAssay(seurat_subsample) <- "RNA"
    seurat_subsample <- NormalizeData(seurat_subsample)
    features_rna <- rownames(seurat_subsample)
    
    DefaultAssay(seurat_subsample) <- "ATAC"
    seurat_subsample <- RunTFIDF(seurat_subsample)
    features_atac <- rownames(seurat_subsample)
    
    pca_result <-  pca_feature_embedding(seurat_subsample, seurat_subsample, d = d, assay1= assay1,  assay2=assay2)
    
    svd_emb_list[[sample_to_use]] <- pca_result[[1]]
    k_list <- c(k_list, pca_result[[2]])
  }
  svd_emb_list <- lapply(svd_emb_list, function(x) x[, 1:max(k_list)])
  return(svd_emb_list)
}

add_split_samples_num_fixed <- function(seurat_obj, sample_col = "sample_name", 
                                        n_split = 10, num_cell = 500, seed = 123) {
  set.seed(seed)
  metadata <- seurat_obj@meta.data
  metadata$sample_name_split <- NA
  unique_samples <- unique(metadata[[sample_col]])
  
  for (sample in unique_samples) {
    sample_indices <- which(metadata[[sample_col]] == sample)
    n_cells <- length(sample_indices)
    
    # Check if enough cells
    needed <- n_split * num_cell
    if (n_cells < needed) {
      warning(sprintf("%s: only %d cells, need %d (%d splits x %d cells). Skipping.",
                      sample, n_cells, needed, n_split, num_cell))
      next
    }
    
    # shuffle and take n_split * num_cell cells
    shuffled_indices <- sample(sample_indices, size = needed)
    
    for (i in 1:n_split) {
      start <- (i - 1) * num_cell + 1
      end <- i * num_cell
      split_indices <- shuffled_indices[start:end]
      split_name <- sprintf("%s_split_%03d", sample, i)
      metadata$sample_name_split[split_indices] <- split_name
    }
    # Remaining cells stay NA
  }
  
  seurat_obj@meta.data$sample_name_split <- metadata$sample_name_split
  return(seurat_obj)
}

run_MOFA_subsamples_parallel <- function(seurat_obj, split_by = "sample_name_split", n_gpus = 8) {
  
  seurat_list <- SplitObject(seurat_obj, split.by = split_by)
  sample_names <- names(seurat_list)
  
  plan(multisession, workers = min(n_gpus, length(sample_names)))
  
  futures <- list()
  for (i in seq_along(sample_names)) {
    # Extract just this one subsample
    one_sample <- seurat_list[[sample_names[i]]]
    sample_name <- sample_names[i]
    gpu_id <- (i - 1) %% n_gpus
    
    futures[[i]] <- future({
      
      library(Seurat)
      library(MOFA2)
      library(reticulate)
      use_python(MOFA_SIMBA_PYTHON, required = TRUE)
      
      # preload CUDA libraries inside each worker
      py_run_string(sprintf("
import ctypes, glob, os
base = '%s'
for lib_dir in glob.glob(os.path.join(base, '*/lib')):
    for so in glob.glob(os.path.join(lib_dir, '*.so*')):
        try:
            ctypes.CDLL(so)
        except:
            pass
", MOFA_SIMBA_NVIDIA_DIR))
      
      DefaultAssay(one_sample) <- "RNA"
      one_sample <- NormalizeData(one_sample)
      features_rna <- rownames(one_sample)
      
      DefaultAssay(one_sample) <- "ATAC"
      one_sample <- RunTFIDF(one_sample)
      features_atac <- rownames(one_sample)
      
      mofa <- create_mofa(
        one_sample,
        assays = c("RNA", "ATAC"),
        features = list(RNA = features_rna, ATAC = features_atac)
      )
      model_opts <- get_default_model_options(mofa)
      train_opts <- get_default_training_options(mofa)
      train_opts$gpu_mode <- TRUE
      train_opts$gpu_device <- gpu_id
      mofa <- prepare_mofa(mofa,
                           model_options = model_opts,
                           training_options = train_opts)
      mofa <- run_mofa(mofa, use_basilisk = FALSE)
      
      mofa_r2_values <- rowSums2(mofa@cache$variance_explained$r2_per_factor$group1)
      mofa_k <- find_elbow_kneedle(mofa_r2_values)
      mofa_k <- max(mofa_k, 2)
      
      W <- get_weights(mofa)
      mofa_emb <- rbind(W$RNA, W$ATAC)
      
      list(emb = mofa_emb, k = mofa_k, sample = sample_name)
    }, seed = TRUE)
  }
  
  mofa_results <- lapply(futures, value)
  plan(sequential)
  
  mofa_num_eigen_list <- sapply(mofa_results, function(x) x$k)
  chosen_eigen <- median(mofa_num_eigen_list)
  
  mofa_emb_list <- lapply(mofa_results, function(x) x$emb[, 1:chosen_eigen])
  names(mofa_emb_list) <- sapply(mofa_results, function(x) x$sample)
  
  return(mofa_emb_list)
}

calculate_frobenius <- function(X, Y){
  sqrt(sum((X - Y)^2))
}

procrustes_align_generalized <- function(emb_list, max_iter = 100, tol = 1e-6) {
  # Ensure shared features in same order
  shared_features <- Reduce(intersect, lapply(emb_list, rownames))
  cat("Shared features across all embeddings:", length(shared_features), "\n")
  emb_list <- lapply(emb_list, function(m) m[shared_features, , drop = FALSE])
  
  aligned <- emb_list
  
  for (iter in seq_len(max_iter)) {
    # Compute consensus (mean of all aligned embeddings)
    consensus <- Reduce("+", aligned) / length(aligned)
    
    # Align each embedding to the consensus
    new_aligned <- lapply(aligned, function(m) {
      proc <- vegan::procrustes(consensus, m, symmetric = TRUE)
      proc$Yrot
    })
    
    # Check convergence: change in consensus
    new_consensus <- Reduce("+", new_aligned) / length(new_aligned)
    delta <- sum((new_consensus - consensus)^2)
    
    if (delta < tol) {
      cat("Converged at iteration", iter, "with delta =", delta, "\n")
      aligned <- new_aligned
      break
    }
    aligned <- new_aligned
  }
  
  if (delta >= tol) cat("Did not converge after", max_iter, "iterations. Final delta =", delta, "\n")
  
  names(aligned) <- names(emb_list)
  aligned
}

run_simba_subsamples <- function(seurat_obj,
                                 split_by = "sample_name_split",
                                 rna_assay = "RNA",
                                 atac_assay = "ATAC",
                                 n_bins = 5,
                                 n_dims = 50,
                                 pbg_workers = 65,
                                 pbg_epochs = 10,
                                 workdir_base = file.path(RESULTS_DIR, "fig2", "simba_subsample"),
                                 graph_name = "graph0",
                                 model_name = "model") {
  
  simba_emb_list <- list()
  seurat_list <- SplitObject(seurat_obj, split.by = split_by)
  
  for (sample_to_use in names(seurat_list)) {
    print(paste0("Processing sample ", length(simba_emb_list) + 1, ": ", sample_to_use))
    
    seurat_subsample <- seurat_list[[sample_to_use]]
    
    # per-sample workdir
    workdir <- file.path(workdir_base, sample_to_use)
    dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
    
    simba_emb <- run_simba_embedding(
      combined_seurat = seurat_subsample,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      n_bins = n_bins,
      n_dims = n_dims,
      pbg_workers = pbg_workers,
      pbg_epochs = pbg_epochs,
      workdir = workdir,
      graph_name = graph_name,
      model_name = model_name
    )
    
    simba_emb_list[[sample_to_use]] <- simba_emb
  }
  
  return(simba_emb_list)
}

build_distance_matrix <- function(emb_list){
  n_sample <- length(emb_list)
  similarity_matrix <- matrix(0, nrow = n_sample, ncol = n_sample)
  for (i in 1:n_sample){
    for (j in 1:n_sample){
      x <- emb_list[[i]]
      y <- emb_list[[j]]
      similarity_matrix[i,j] <- calculate_frobenius(x,y)
    }
  }
  rownames(similarity_matrix) <- names(emb_list)
  colnames(similarity_matrix) <- names(emb_list)
  return(similarity_matrix)
}

plot_umap_scorch <- function(dis_matrix,title=""){
  umap_result <-  umap(as.dist(dis_matrix), n_neighbors = 15, min_dist = 0.1)
  
  condition <- sub("_split_.*$", "", colnames(dis_matrix))
  
  sorted_donors <- sort(unique(condition))
  donor_map <- setNames(
    paste0("Donor ", seq_along(sorted_donors)),
    sorted_donors
  )
  donor_label <- donor_map[condition]
  
  disease_state <- ifelse(grepl("HIV", condition), "HIV", "CTR")
  condition_numeric <- as.integer(as.factor(condition))
  silhouette_scores <- silhouette(condition_numeric, dis_matrix)
  average_silhouette_score <- mean(silhouette_scores[, "sil_width"])
  print(paste0("sil score: ",average_silhouette_score))
  
  umap_df <- data.frame(
    UMAP_1 = umap_result[, 1],
    UMAP_2 = umap_result[, 2],
    donor = donor_label,
    condition = disease_state
  )
  
  p_umap <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = donor)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(aes(group = donor), level = 0.8, linetype = "dashed", 
                 linewidth = 0.4, show.legend = FALSE) +
    labs(title = title, color = "Donor") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10))
  
  p_umap
  
}

cal_sil_score_scorch <- function(dis_matrix){
  condition <- sub("_split_.*$", "", colnames(dis_matrix))
  disease_state <- ifelse(grepl("HIV", condition), "HIV", "CTR")
  condition_numeric <- as.integer(as.factor(condition))
  silhouette_scores <- silhouette(condition_numeric, dis_matrix)
  average_silhouette_score <- mean(silhouette_scores[, "sil_width"])
  return(average_silhouette_score)
}

plot_pcoa_scorch <- function(dis_matrix, title=""){
  
  # PCoA (classical MDS)
  pcoa_model <- cmdscale(as.dist(dis_matrix), k = 2, eig = TRUE)
  
  # variance explained per axis
  var_explained <- round(pcoa_model$eig / sum(pcoa_model$eig[pcoa_model$eig > 0]) * 100, 1)
  
  # donor labels
  condition <- sub("_split_.*$", "", colnames(dis_matrix))
  disease_state <- ifelse(grepl("HIV", condition), "HIV", "CTR")
  condition_numeric <- as.integer(as.factor(condition))
  
  
  sorted_donors <- sort(unique(condition))
  donor_map <- setNames(
    paste0("Donor ", seq_along(sorted_donors)),
    sorted_donors
  )
  donor_label <- donor_map[condition]
  
  # silhouette by donor
  silhouette_scores <- silhouette(condition_numeric, as.dist(dis_matrix))
  average_silhouette_score <- mean(silhouette_scores[, "sil_width"])
  print(paste0("sil score: ", round(average_silhouette_score, 3)))
  
  pcoa_df <- data.frame(
    PCoA_1 = pcoa_model$points[, 1],
    PCoA_2 = pcoa_model$points[, 2],
    donor = donor_label,
    condition = disease_state
  )
  
  p_pcoa <- ggplot(pcoa_df, aes(x = PCoA_1, y = PCoA_2, color = donor)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(aes(group = donor), level = 0.8, linetype = "dashed", 
                 linewidth = 0.4, show.legend = FALSE) +
    labs(title = title, 
         color = "Donor", 
         x = paste0("PCoA 1 (", var_explained[1], "%)"),
         y = paste0("PCoA 2 (", var_explained[2], "%)")) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8))
  
  return(p_pcoa)
}

quietly_run <- function(expr) {
  expr <- substitute(expr)
  invisible(
    capture.output(
      suppressMessages(
        suppressWarnings(
          eval(expr, envir = parent.frame())
        )
      )
    )
  )
}

# ---- Run benchmark ----

mofa_raw_emb_list <- list()
simba_raw_emb_list <- list()
mosaic_emb_list <- list()
pca_raw_emb_list <- list()

mofa_aligned_emb_list <- list()
simba_aligned_emb_list <- list()
pca_aligned_emb_list <- list()

# 10 pseudo-replicates per donor, cell_num cells each
for (split_seed in 42:42){
  print(split_seed)
  
  for (cell_num in c(100,200,300,400,500,600,700)){
    round_name <- paste0("Processing: Random Seed=",split_seed,", Cell Number=",cell_num)
    cat(round_name)
    seurat_rna <- add_split_samples_num_fixed(seurat_rna, sample_col = "sample_name",n_split=10,num_cell=cell_num,seed = split_seed)
    seurat_rna_subset <- seurat_rna[, !is.na(seurat_rna$sample_name_split)]
    seurat_atac_subset <- seurat_atac[,colnames(seurat_rna_subset)]
    seurat_atac_subset$sample_name_split <- seurat_rna$sample_name_split
    
    # keep features with total count >= 3 in every pseudo-replicate
    selected_gene <- get_genes_expressed_across_samples(seurat_rna_subset,"sample_name_split",min_expr=3)
    selected_peak <- get_genes_expressed_across_samples(seurat_atac_subset,"sample_name_split",min_expr=3)
    
    seurat_rna_subset <- seurat_rna_subset[selected_gene,]
    seurat_atac_subset <- seurat_atac_subset[selected_peak,]
    
    seurat_rna_subset <- NormalizeData(seurat_rna_subset, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat_atac_subset <- RunTFIDF(seurat_atac_subset)
    
    if(length(selected_gene) > 3000){
      seurat_rna_subset <- FindVariableFeatures(seurat_rna_subset,nfeatures=3000)
      variable_genes <- VariableFeatures(seurat_rna_subset)
      seurat_rna_subset <- seurat_rna_subset[variable_genes,]
    }
    
    combined_seurat <- seurat_rna_subset
    combined_seurat[["ATAC"]] <- seurat_atac_subset[["ATAC"]]
    
    cat("Start MOSAIC:")
    quietly_run({
      result_no_shuffle <- run_our_method_two(seurat_rna_subset,seurat_atac_subset,assay1 = "RNA", assay2="ATAC",
                                              sample_meta="sample_name_split",condition_meta="sample_name",block_method="cosine")
      mosaic_emb <- result_no_shuffle$projected_combed_list
    })
    cat("Finish MOSAIC!")
    
    cat("Start MOFA+:")
    quietly_run({
      mofa_emb_raw <- run_MOFA_subsamples_parallel(combined_seurat,split_by="sample_name_split")
      mofa_emb_aligned <- procrustes_align_generalized(mofa_emb_raw)
      rownames(mofa_emb_aligned) <- rownames(mofa_emb_raw)
      colnames(mofa_emb_aligned) <- colnames(mofa_emb_raw)
    })
    cat("Finish MOFA+!")
    
    cat("Start SIMBA:")
    quietly_run({
      simba_emb_raw <- run_simba_subsamples(combined_seurat,split_by="sample_name_split",
                                            workdir_base = simba_workdir)
      simba_emb_aligned <- procrustes_align_generalized(simba_emb_raw)
      rownames(simba_emb_aligned) <- rownames(simba_emb_raw)
      colnames(simba_emb_aligned) <- colnames(simba_emb_raw)
    })
    cat("Finish SIMBA!")
    
    cat("Start PCA:")
    quietly_run({
      pca_emb_raw <- run_pca_subsamples(combined_seurat,split_by="sample_name_split")
      pca_emb_aligned <- procrustes_align_generalized(pca_emb_raw)
      rownames(pca_emb_aligned) <- rownames(pca_emb_raw)
      colnames(pca_emb_aligned) <- colnames(pca_emb_raw)
    })
    cat("Finish PCA!")
    
    mofa_raw_emb_list[[round_name]] <- mofa_emb_raw
    simba_raw_emb_list[[round_name]] <- simba_emb_raw
    mosaic_emb_list[[round_name]] <- mosaic_emb
    pca_raw_emb_list[[round_name]] <- pca_emb_raw
    
    mofa_aligned_emb_list[[round_name]] <- mofa_emb_aligned
    simba_aligned_emb_list[[round_name]] <- simba_emb_aligned
    pca_aligned_emb_list[[round_name]] <- pca_emb_aligned
    
    if (dir.exists(simba_workdir)) {
      unlink(simba_workdir,
             recursive = TRUE,
             force = TRUE)
    }
  }
}

# ---- Save embeddings ----

out_dir <- file.path(RESULTS_DIR, "fig2", "sample_embed_benchmark")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

all_lists <- list(
  mofa_raw_emb = mofa_raw_emb_list,
  simba_raw_emb = simba_raw_emb_list,
  mosaic_emb = mosaic_emb_list,
  pca_raw_emb = pca_raw_emb_list,
  mofa_aligned_emb = mofa_aligned_emb_list,
  simba_aligned_emb = simba_aligned_emb_list,
  pca_aligned_emb = pca_aligned_emb_list
)

for (name in names(all_lists)) {
  saveRDS(all_lists[[name]], file = file.path(out_dir, paste0(name, ".rds")))
}

# ---- Figure 2C ----

# silhouette score per run
process_emb_list <- function(emb_list, method_name) {
  
  sil_scores <- numeric(length(emb_list))
  
  for (i in seq_along(emb_list)) {
    emb <- emb_list[[i]]
    
    # normalized distance matrix
    dist_mat <- build_distance_matrix(emb)
    dist_mat <- dist_mat / max(dist_mat)
    
    sil_scores[i] <- cal_sil_score_scorch(dist_mat)
  }
  
  df <- data.frame(
    name = names(emb_list),
    sil_score = sil_scores,
    method = method_name,
    stringsAsFactors = FALSE
  )
  
  # parse seed and cell number from the run name
  df <- df %>%
    mutate(
      seed = as.numeric(str_extract(name, "(?<=Random Seed=)\\d+")),
      cell_number = as.numeric(str_extract(name, "(?<=Cell Number=)\\d+"))
    )
  
  return(df)
}

df_mosaic <- process_emb_list(mosaic_emb_list, "MOSAIC")
df_mofa   <- process_emb_list(mofa_aligned_emb_list, "MOFA")
df_simba  <- process_emb_list(simba_aligned_emb_list, "SIMBA")
df_pca    <- process_emb_list(pca_aligned_emb_list, "PCA")

df_all <- bind_rows(df_mosaic, df_mofa, df_simba, df_pca)

# mean and SD per method and cell number
df_summary <- df_all %>%
  group_by(method, cell_number) %>%
  summarise(
    mean_sil = mean(sil_score),
    sd_sil = sd(sil_score),
    .groups = "drop"
  )

df_summary <- df_summary %>%
  dplyr::mutate(method = dplyr::if_else(method == "MOFA", "MOFA+ + Procrustes", method))%>%
  dplyr::mutate(method = dplyr::if_else(method == "SIMBA", "SIMBA + Procrustes", method))%>%
  dplyr::mutate(method = dplyr::if_else(method == "PCA", "PCA + Procrustes", method))

df_summary_filtered <- df_summary %>%
  filter(!cell_number %in% c(100, 700))

custom_pal <- c(
  "MOSAIC" = "#1B9E77",
  "MOFA+ + Procrustes" = "#7570B3",
  "PCA + Procrustes" = "#666666",
  "SIMBA + Procrustes" = "#E6AB02"
)

legend_order <- df_summary_filtered %>%
  group_by(method) %>%
  filter(cell_number == max(cell_number)) %>%
  arrange(desc(mean_sil)) %>%
  pull(method)

df_summary_filtered$method <- factor(
  df_summary_filtered$method,
  levels = legend_order
)

benchmark_plot <- ggplot(df_summary_filtered, aes(x = cell_number, y = mean_sil, color = method, group = method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_sil - sd_sil, ymax = mean_sil + sd_sil), width = 15, size = 0.8) +
  
  scale_x_continuous(breaks = unique(df_summary$cell_number)) +
  labs(
    x = "Cell Number",
    y = "Average Silhouette Score",
    color = "Method"
  ) +
  scale_color_manual(values = custom_pal) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_line(color = "grey92", linetype = "dashed", size = 0.5),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  )+scale_y_continuous(
    breaks = seq(0, 0.4, by = 0.1)
  )

png(file.path(FIG_DIR, "fig2C.png"), width = 9, height = 6, units = "in", res = 600)
print(benchmark_plot)
dev.off()

write.csv(df_all, file.path(out_dir, "fig2C_silhouette.csv"), row.names = FALSE)

# ---- Figure 2D and S1 (200 cells) ----
# element [[2]] is cell number 200
dist_mofa <- build_distance_matrix(mofa_aligned_emb_list[[2]])
dist_mofa <- dist_mofa/max(dist_mofa)
dist_mosaic <- build_distance_matrix(mosaic_emb_list[[2]])
dist_mosaic <- dist_mosaic/max(dist_mosaic)
dist_simba <- build_distance_matrix(simba_aligned_emb_list[[2]])
dist_simba <- dist_simba/max(dist_simba)
dist_pca <- build_distance_matrix(pca_aligned_emb_list[[2]])
dist_pca <- dist_pca/max(dist_pca)

# Figure 2D: UMAP of pairwise distances
set.seed(46)
umap_mofa   <- plot_umap_scorch(dist_mofa,title="MOFA+ + Procrustes")
umap_mosaic <- plot_umap_scorch(dist_mosaic,title="MOSAIC")
umap_pca    <- plot_umap_scorch(dist_pca,title="PCA + Procrustes")
umap_simba  <- plot_umap_scorch(dist_simba,title="SIMBA + Procrustes")

fig2D_plot <- ((umap_mosaic | umap_mofa | umap_simba)) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")

png(file.path(FIG_DIR, "fig2D.png"),width = 14, height = 5, units = "in", res = 600)
print(fig2D_plot)
dev.off()

# Figure S1: PCoA of pairwise distances
pcoa_mofa   <- plot_pcoa_scorch(dist_mofa,title="MOFA+ + Procrustes")
pcoa_mosaic <- plot_pcoa_scorch(dist_mosaic,title="MOSAIC")
pcoa_pca    <- plot_pcoa_scorch(dist_pca,title="PCA + Procrustes")
pcoa_simba  <- plot_pcoa_scorch(dist_simba,title="SIMBA + Procrustes")

figS1_plot <- ((pcoa_mosaic | pcoa_mofa) /(pcoa_simba | pcoa_pca)) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")

png(file.path(FIG_DIR, "figS1.png"),width = 9, height =8, units = "in", res = 600)
print(figS1_plot)
dev.off()

# average silhouette at 200 cells
sil_list <- list("MOFA+ + Procrustes" = cal_sil_score_scorch(dist_mofa),
                 "MOSAIC" = cal_sil_score_scorch(dist_mosaic),
                 "PCA + Procrustes" = cal_sil_score_scorch(dist_pca),
                 "SIMBA + Procrustes" = cal_sil_score_scorch(dist_simba))

sil_df <- data.frame(
  Method = names(sil_list),
  Silhouette = unlist(sil_list)
)
print(sil_df)
write.csv(sil_df, file.path(out_dir, "silhouette_200cells.csv"), row.names = FALSE)
