# Figure S10: four-way ablation of the two truncation stages on the Figure 2A feature-module benchmark
# Output:   results/S10/seed42_fig2a_ablation.csv, figures/figS10.png

source("config.R")
source("R/functions.R")            # eigen_block_cosine_two(), find_elbow_kneedle()
source("fig2/fig2A_functions.R")   # generate_svd_simulated_data()

suppressPackageStartupMessages({
  library(Seurat); library(RSpectra); library(mclust)
  library(dplyr); library(tidyr); library(ggplot2)
})

out_dir <- file.path(RESULTS_DIR, "S10")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

NOISE_GRID <- seq(0.25, 1, 0.25)

# ---- MOSAIC with switches for Stage 1 and Stage 2 ----
run_our_method_two_ablation <- function(rna_seurat, atac_seurat,
                                        assay1 = "originalexp",
                                        assay2 = "originalexp",
                                        sample_meta = "sample_id",
                                        condition_meta = "condition",
                                        block_method = "cosine",
                                        stage1 = TRUE, stage2 = TRUE) {
  
  rna_seurat  <- ScaleData(rna_seurat,  features = rownames(rna_seurat))
  atac_seurat <- ScaleData(atac_seurat, features = rownames(atac_seurat))
  
  rna.list  <- SplitObject(rna_seurat,  split.by = sample_meta)
  atac.list <- SplitObject(atac_seurat, split.by = sample_meta)
  
  dataset_list <- list()
  for (i in 1:length(rna.list)) {
    gm <- t(as.matrix(rna.list[[i]][[assay1]]$scale.data))
    am <- t(as.matrix(atac.list[[i]][[assay2]]$scale.data))
    dataset_list[[i]] <- list(gm, am)
  }
  
  eigen_output_cosine_list <- list(); r_cosine_list <- c()
  block_matrix_list <- list()
  for (i in 1:length(dataset_list)) {
    ms <- dataset_list[[i]]
    eo <- eigen_block_cosine_two(ms[[1]], ms[[2]])
    r_cosine_list <- c(r_cosine_list, eo[[2]])
    eigen_output_cosine_list[[i]] <- eo[c(1, 2)]
    block_matrix_list[[i]] <- eo[[3]]
  }
  names(eigen_output_cosine_list) <- names(rna.list)
  names(block_matrix_list)        <- names(rna.list)
  
  # ---- Stage 1 ----
  if (stage1) {
    sample_projection_matrix <- list()
    for (i in 1:length(eigen_output_cosine_list)) {
      ru <- max(2, r_cosine_list[i])
      Vr <- eigen_output_cosine_list[[i]][[1]]$vectors[, 1:ru, drop = FALSE]
      lv <- eigen_output_cosine_list[[i]][[1]]$values[1:ru]
      sample_projection_matrix[[i]] <- Vr %*% diag(lv, nrow = ru) %*% t(Vr)
    }
  } else {
    sample_projection_matrix <- block_matrix_list
  }
  
  sum_projection_matrix <- Reduce("+", sample_projection_matrix)
  
  # ---- Stage 2 ----
  if (stage2) {
    V_result <- RSpectra::eigs_sym(sum_projection_matrix, 50, which = "LA")
    r2 <- find_elbow_kneedle(V_result$values[1:50])
    V_matrix <- V_result$vectors[, 1:r2, drop = FALSE] %*%
      diag(V_result$values[1:r2], nrow = r2)
  } else {
    # full non-zero spectrum (dense eigendecomposition)
    ev <- eigen(sum_projection_matrix, symmetric = TRUE)
    keep <- sum(ev$values > 1e-8)
    r2 <- keep
    V_matrix <- ev$vectors[, 1:keep, drop = FALSE] %*%
      diag(ev$values[1:keep], nrow = keep)
  }
  
  projected_combed_list <- lapply(seq_along(block_matrix_list),
                                  function(i) block_matrix_list[[i]] %*% V_matrix)
  names(projected_combed_list) <- names(eigen_output_cosine_list)
  
  list(projected_combed_list = projected_combed_list,
       r1 = median(r_cosine_list), r2 = r2)
}

# ---- Figure 2A scoring: average, L2-normalise dimensions, k-means, ARI ----
score_fig2a <- function(emb_list, true_labels, km_seed = 1) {
  emb_avg <- Reduce("+", emb_list) / length(emb_list)
  for (j in seq_len(ncol(emb_avg))) {          # per-dimension L2 normalisation
    cn <- sqrt(sum(emb_avg[, j]^2))
    if (cn > 0) emb_avg[, j] <- emb_avg[, j] / cn
  }
  set.seed(km_seed)
  cl <- kmeans(emb_avg, centers = 10, nstart = 25)$cluster
  mclust::adjustedRandIndex(true_labels, cl)
}

simulate_fig2a <- function(feature_noise, seed) {
  generate_svd_simulated_data(
    simulation_type    = "DC",
    add_batch          = FALSE,
    whether_rescale    = FALSE,
    de_prop            = 0,
    seed               = seed,
    sample_noise_level = 0.5,
    feature_noise_level = feature_noise,
    cell_noise_level   = 0.5,
    whether_nonlinear  = TRUE
  )
}

# ---- Step 1: full two-stage arm only, compare with Figure 2A ----
if (TRUE) {
  repro <- do.call(rbind, lapply(NOISE_GRID, function(fn) {
    cat("reproduction, feature_noise =", fn, "\n")
    sim <- simulate_fig2a(fn, seed = 42)
    tl  <- c(sim$seurat_X[["originalexp"]][[]]$feature_cluster,
             sim$seurat_Y[["originalexp"]][[]]$feature_cluster)
    fit <- run_our_method_two_ablation(sim$seurat_X, sim$seurat_Y,
                                       stage1 = TRUE, stage2 = TRUE)
    data.frame(feature_noise = fn,
               ARI = score_fig2a(fit$projected_combed_list, tl),
               r1 = fit$r1, r2 = fit$r2)
  }))
  print(repro, row.names = FALSE)
  cat("\nCompare these four values against the MOSAIC curve in Figure 2A.\n")
}

# ---- Step 2: four arms, one shared dataset per noise level ----
if (TRUE) {
  seeds <- 42
  grid  <- expand.grid(feature_noise = NOISE_GRID, seed = seeds)
  cat(nrow(grid), "datasets x 4 arms\n")
  cat("main R pid:", Sys.getpid(), "\n")
  
  t0 <- Sys.time()
  res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(k) {
    fn <- grid$feature_noise[k]; sd <- grid$seed[k]
    sim <- simulate_fig2a(fn, seed = sd)
    tl  <- c(sim$seurat_X[["originalexp"]][[]]$feature_cluster,
             sim$seurat_Y[["originalexp"]][[]]$feature_cluster)
    
    out <- do.call(rbind, lapply(
      list(c(FALSE, FALSE), c(TRUE, FALSE), c(FALSE, TRUE), c(TRUE, TRUE)),
      function(sw) {
        fit <- run_our_method_two_ablation(sim$seurat_X, sim$seurat_Y,
                                           stage1 = sw[1], stage2 = sw[2])
        data.frame(feature_noise = fn, seed = sd,
                   arm = c("none", "stage1_only", "stage2_only", "both")[
                     1 + sw[1] + 2 * sw[2]],
                   ARI = score_fig2a(fit$projected_combed_list, tl),
                   r1 = fit$r1, r2 = fit$r2)
      }))
    
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("  %d/%d  %.1f min  eta %.1f min\n", k, nrow(grid), el,
                el / k * (nrow(grid) - k)))
    utils::flush.console()
    write.csv(out, file.path(out_dir, sprintf("seed42_fig2a_abl_cell_%03d.csv", k)), row.names = FALSE)
    out
  }))
  write.csv(res, file.path(out_dir, "seed42_fig2a_ablation.csv"), row.names = FALSE)
}

# ---- Summaries and figure ----
if (TRUE) {
  res <- read.csv(file.path(out_dir, "seed42_fig2a_ablation.csv"), stringsAsFactors = FALSE)
  res$arm <- factor(res$arm,
                    levels = c("none", "stage1_only", "stage2_only", "both"))
  
  cat("\n=== ARI by arm and noise level ===\n")
  print(as.data.frame(res %>% dplyr::group_by(feature_noise, arm) %>%
                        dplyr::summarise(ARI = round(mean(ARI), 4),
                                         se = round(sd(ARI) / sqrt(dplyr::n()), 4),
                                         r1 = round(mean(r1), 1), r2 = round(mean(r2), 1),
                                         .groups = "drop")), row.names = FALSE)
  
  w <- res %>% dplyr::select(feature_noise, seed, arm, ARI) %>%
    tidyr::pivot_wider(names_from = arm, values_from = ARI)
  
  cat("\n=== Stage 1 given Stage 2, paired by seed ===\n")
  print(as.data.frame(w %>%
                        dplyr::mutate(dd = both - stage2_only) %>%
                        dplyr::group_by(feature_noise) %>%
                        dplyr::summarise(n = dplyr::n(), delta = mean(dd),
                                         se = sd(dd) / sqrt(dplyr::n()),
                                         t = mean(dd) / (sd(dd) / sqrt(dplyr::n())),
                                         p = tryCatch(t.test(dd)$p.value, error = function(e) NA_real_),
                                         ci_lo = tryCatch(t.test(dd)$conf.int[1], error = function(e) NA_real_),
                                         ci_hi = tryCatch(t.test(dd)$conf.int[2], error = function(e) NA_real_),
                                         frac_pos = mean(dd > 0), .groups = "drop") %>%
                        dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 5)))),
        row.names = FALSE)
  
  cat("\n=== Stage 2 given nothing ===\n")
  print(as.data.frame(w %>%
                        dplyr::mutate(dd = stage2_only - none) %>%
                        dplyr::group_by(feature_noise) %>%
                        dplyr::summarise(delta = round(mean(dd), 4),
                                         se = round(sd(dd) / sqrt(dplyr::n()), 4),
                                         .groups = "drop")), row.names = FALSE)

}

res$arm <- factor(
  res$arm,
  levels = c("none", "stage1_only", "stage2_only", "both"),
  labels = c(
    "No truncation",
    "Per-sample only",
    "Post-aggregation only",
    "Full two-stage"
  )
)

p <- ggplot(res, aes(feature_noise, ARI, colour = arm)) +
  stat_summary(
    aes(linetype = arm),
    fun = mean, geom = "line",
    linewidth = 0.8
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.02
  ) +
  stat_summary(
    aes(shape = arm),
    fun = mean, geom = "point",
    size = 2.5
  ) +
  scale_linetype_manual(
    values = c("solid", "solid", "solid", "dashed")
  ) +
  scale_shape_manual(
    values = c(16, 16, 16, 17)
  ) +
  labs(
    x = expression(sigma[feature]),
    y = "ARI",
    colour = NULL,
    linetype = NULL,
    shape = NULL
  ) +
  theme_minimal()

png(file.path(FIG_DIR, "figS10.png"),width = 5.5, height = 4, units = "in", res = 600)
print(p)
dev.off()
