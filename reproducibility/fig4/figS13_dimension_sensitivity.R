# Supplemental Figure S13: sensitivity of the T-cell DC results to the embedding dimension d
# Input:    results/fig4/seurat_cd4naive.rds (or CITESEQ_H5SEURAT)
# Output:   figures/figS13.png, results/S13/sweep_metrics_full.csv
# Requires: fig4/fig4_dc_paired_test.R

source("config.R")
source("R/functions.R")   # eigen_block_cosine_two, find_elbow_kneedle, cal_peranova_per_feature_parallel_no_zero

suppressPackageStartupMessages({
  library(RSpectra)
  library(Seurat)
  library(dplyr)
})

OUT_DIR   <- file.path(RESULTS_DIR, "S13")
CACHE_RDS <- file.path(OUT_DIR, "sweep_pipeline_cache.rds")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_EIG_SWEEP <- 100   # eigenpairs extracted from the aggregate
D0_WINDOW   <- 50    # Kneedle window for d0
D_GRID      <- 2:50  # sweep range
TOPK        <- c(100, 300, 500)
TARGET_FEATURE <- "RNA_STAT5B"
NUM_PERM    <- 999   # only F_stats_list is used

diag_lines <- c()
log_diag <- function(...) {
  msg <- paste0(...)
  message(msg)
  diag_lines <<- c(diag_lines, msg)
}


## 1. find_elbow_kneedle() that also returns the distance curve

find_elbow_kneedle_curve <- function(eigenvalues) {
  ev_raw <- eigenvalues
  n_nonpos <- sum(eigenvalues <= 0)
  
  # carry forward non-positive eigenvalues
  for (i in seq_along(eigenvalues)) {
    if (eigenvalues[i] <= 0) {
      eigenvalues[i] <- eigenvalues[i - 1]
    }
  }
  
  x <- seq_along(eigenvalues)
  y <- log(eigenvalues)
  
  if (length(x) < 3) {
    return(list(elbow = length(x), distances = rep(NA_real_, length(x)),
                x_scaled = NA, y_scaled = NA, n_nonpos = n_nonpos))
  }
  
  x_scaled <- (x - min(x)) / (max(x) - min(x))
  y_scaled <- (y - min(y)) / (max(y) - min(y))
  
  distances <- perpendicular_distance(
    xp = x_scaled, yp = y_scaled,
    x1 = x_scaled[1], y1 = y_scaled[1],
    x2 = x_scaled[length(x_scaled)], y2 = y_scaled[length(y_scaled)]
  )
  
  elbow_index <- which.max(distances)
  list(
    elbow     = max(elbow_index, 2),
    distances = distances,
    x_scaled  = x_scaled,
    y_scaled  = y_scaled,
    ev_used   = eigenvalues,
    ev_raw    = ev_raw,
    n_nonpos  = n_nonpos
  )
}

# check against find_elbow_kneedle()
stopifnot(
  identical(
    find_elbow_kneedle_curve(c(10, 8, 6, 2, 1.8, 1.7, 1.6))$elbow,
    find_elbow_kneedle(c(10, 8, 6, 2, 1.8, 1.7, 1.6))
  )
)


## 2. MOSAIC pipeline returning the aggregate spectrum (n_eig eigenpairs)

run_our_method_two_sweep <- function(rna_seurat, atac_seurat,
                                     assay1 = "SCT", assay2 = "ADT",
                                     sample_meta = "sample_id",
                                     condition_meta = "condition",
                                     block_method = "cosine",
                                     n_eig = 100, d0_window = 50) {
  
  rna_seurat  <- ScaleData(rna_seurat,  features = rownames(rna_seurat))
  atac_seurat <- ScaleData(atac_seurat, features = rownames(atac_seurat))
  
  rna.list  <- SplitObject(rna_seurat,  split.by = sample_meta)
  atac.list <- SplitObject(atac_seurat, split.by = sample_meta)
  
  dataset_list <- list()
  for (i in seq_along(rna.list)) {
    gene_matrix <- t(as.matrix(rna.list[[i]][[assay1]]$scale.data))
    adt_matrix  <- t(as.matrix(atac.list[[i]][[assay2]]$scale.data))
    dataset_list[[i]] <- list(gene_matrix, adt_matrix)
  }
  
  eigen_output_cosine_list <- list()
  r_cosine_list <- c()
  block_matrix_list <- list()
  for (i in seq_along(dataset_list)) {
    message("stage 1, sample ", i)
    ms <- dataset_list[[i]]
    eo <- eigen_block_cosine_two(ms[[1]], ms[[2]])
    r_cosine_list <- c(r_cosine_list, eo[[2]])
    eigen_output_cosine_list[[i]] <- eo[c(1, 2)]
    block_matrix_list[[i]] <- eo[[3]]
  }
  r <- round(median(r_cosine_list))
  message("per-sample median r: ", r)
  
  names(eigen_output_cosine_list) <- names(rna.list)
  names(block_matrix_list)        <- names(rna.list)
  
  sample_condition_map <- unique(data.frame(
    sample_id = rna_seurat[[sample_meta]],
    condition = rna_seurat[[condition_meta]]
  ))
  colnames(sample_condition_map) <- c("sample_id", "condition")
  annotation <- data.frame(Condition = sample_condition_map$condition)
  rownames(annotation) <- sample_condition_map$sample_id
  annotation <- annotation[names(eigen_output_cosine_list), , drop = FALSE]
  
  # stage 1 reconstruction
  sample_projection_matrix <- list()
  for (i in seq_along(eigen_output_cosine_list)) {
    r_to_use <- r_cosine_list[i]
    if (r_to_use == 1) r_to_use <- 2
    V  <- eigen_output_cosine_list[[i]][[1]]$vectors[, 1:r_to_use]
    lm <- eigen_output_cosine_list[[i]][[1]]$values[1:r_to_use]
    sample_projection_matrix[[i]] <- V %*% diag(lm) %*% t(V)
  }
  
  sum_projection_matrix <- Reduce("+", sample_projection_matrix)
  
  # n_eig eigenpairs
  V_result <- RSpectra::eigs_sym(sum_projection_matrix, n_eig, which = "LA")
  
  # d0: Kneedle on the top d0_window values
  d0 <- find_elbow_kneedle(V_result$values[1:d0_window])
  message("d0 (window ", d0_window, "): ", d0)
  
  list(
    V_values          = V_result$values,
    V_vectors         = V_result$vectors,
    block_matrix_list = block_matrix_list,
    annotation        = annotation,
    r_cosine_list     = r_cosine_list,
    r_median          = r,
    d0                = d0,
    sample_names      = names(eigen_output_cosine_list)
  )
}


## 3. Run (or load) the pipeline on the Figure 4 data (CD4 naive, Day 0 and Day 7)
fig4_obj <- file.path(RESULTS_DIR, "fig4", "seurat_cd4naive.rds")
if (file.exists(fig4_obj)) {
  obj <- readRDS(fig4_obj)
  seurat_rna    <- obj$seurat_rna
  seurat_adt    <- obj$seurat_adt
  feature_names <- obj$feature_names
  rm(obj)
} else {
  source("fig4/fig4_preprocess.R")
}



if (file.exists(CACHE_RDS)) {
  log_diag("loading cached pipeline from ", CACHE_RDS)
  fit <- readRDS(CACHE_RDS)
} else {
  fit <- run_our_method_two_sweep(
    seurat_rna, seurat_adt,
    sample_meta    = "sample_id",
    condition_meta = "time",
    block_method   = "cosine",
    n_eig          = N_EIG_SWEEP,
    d0_window      = D0_WINDOW
  )
}

d0     <- fit$d0
groups <- fit$annotation[fit$sample_names, "Condition"]
n_sample <- length(fit$block_matrix_list)

log_diag("d0 = ", d0)
log_diag("per-sample stage-1 median rank r = ", fit$r_median)
log_diag("per-sample stage-1 rank range = [",
         min(fit$r_cosine_list), ", ", max(fit$r_cosine_list), "]")
log_diag("n_sample = ", n_sample)
log_diag("group sizes: ", paste(names(table(groups)), table(groups),
                                sep = "=", collapse = ", "))

n_nonpos <- sum(fit$V_values[1:N_EIG_SWEEP] <= 0)
log_diag("non-positive eigenvalues among top ", N_EIG_SWEEP, ": ", n_nonpos,
         if (n_nonpos > 0)
           "  [WARNING: Kneedle carry-forward triggered; spectrum panel will show an artificial plateau]"
         else "")


## 4. Panel A/B data: spectrum, Kneedle distance curve, cumulative squared spectral weight

kn <- find_elbow_kneedle_curve(fit$V_values[1:D0_WINDOW])
stopifnot(kn$elbow == d0)

lam2   <- fit$V_values[1:N_EIG_SWEEP]^2
cum_w  <- cumsum(lam2) / sum(lam2)

spectrum_df <- data.frame(
  index                 = 1:N_EIG_SWEEP,
  eigenvalue            = fit$V_values[1:N_EIG_SWEEP],
  log_eigenvalue        = log(pmax(fit$V_values[1:N_EIG_SWEEP], .Machine$double.eps)),
  kneedle_distance      = c(kn$distances, rep(NA_real_, N_EIG_SWEEP - D0_WINDOW)),
  cum_sq_spectral_weight = cum_w,
  is_d0                 = (1:N_EIG_SWEEP) == d0
)
write.csv(spectrum_df, file.path(OUT_DIR, "spectrum_d0.csv"), row.names = FALSE)

log_diag("cumulative squared spectral weight at d0 = ",
         sprintf("%.4f", cum_w[d0]))
for (m in c(1.25, 1.5, 2)) {
  dd <- min(round(m * d0), N_EIG_SWEEP)
  log_diag("  at ", m, "x d0 (d=", dd, ") = ", sprintf("%.4f", cum_w[dd]))
}


## 5. Search-window check

window_df <- do.call(rbind, lapply(c(25, 50, 75, 100), function(L) {
  if (L > N_EIG_SWEEP) return(NULL)
  data.frame(window_L = L, d_selected = find_elbow_kneedle(fit$V_values[1:L]))
}))
write.csv(window_df, file.path(OUT_DIR, "window_check.csv"), row.names = FALSE)
log_diag("search-window check: ",
         paste(window_df$window_L, window_df$d_selected, sep = "->", collapse = ", "))


## 6. Nested projection: column j of W does not depend on d, so E_i^(d) = E_i[, 1:d]

VL_full <- fit$V_vectors[, 1:N_EIG_SWEEP] %*% diag(fit$V_values[1:N_EIG_SWEEP])

E_full <- lapply(seq_len(n_sample), function(i) {
  message("projecting sample ", i)
  fit$block_matrix_list[[i]] %*% VL_full
})
names(E_full) <- fit$sample_names

# check nesting against a direct recomputation at d0
VL_d0    <- fit$V_vectors[, 1:d0] %*% diag(fit$V_values[1:d0])
E_direct <- fit$block_matrix_list[[1]] %*% VL_d0
stopifnot(max(abs(E_direct - E_full[[1]][, 1:d0])) < 1e-8)
log_diag("column-nesting verified (max abs deviation < 1e-8)")

n_features <- nrow(E_full[[1]])
log_diag("n_features = ", n_features)


## 7. Pseudo-F for all d from cumulative per-coordinate sums of squares

D_MAX   <- max(D_GRID)
G       <- length(unique(groups))
S       <- n_sample
grp_idx <- split(seq_len(S), as.factor(groups))

SST <- matrix(NA_real_, n_features, D_MAX)
SSW <- matrix(NA_real_, n_features, D_MAX)

for (j in seq_len(D_MAX)) {
  X  <- vapply(E_full, function(M) M[, j], numeric(n_features))  # features x samples
  SST[, j] <- rowSums((X - rowMeans(X))^2)
  w <- numeric(n_features)
  for (ix in grp_idx) {
    Xg <- X[, ix, drop = FALSE]
    w  <- w + rowSums((Xg - rowMeans(Xg))^2)
  }
  SSW[, j] <- w
}

# cumulative sums over coordinates: (SS %*% U)[f,k] = sum_{j<=k} SS[f,j]
Ucum <- upper.tri(matrix(0, D_MAX, D_MAX), diag = TRUE) * 1
cumT <- SST %*% Ucum
cumW <- SSW %*% Ucum
Fall <- ((cumT - cumW) / (G - 1)) / (cumW / (S - G))
colnames(Fall) <- paste0("d", seq_len(D_MAX))

## check against cal_peranova_per_feature_parallel_no_zero at d0
Fref <- unlist(cal_peranova_per_feature_parallel_no_zero(
  lapply(E_full, function(M) M[, 1:d0, drop = FALSE]), n_sample,
  groups = groups, dist_method = "eu", num_permutation = NUM_PERM
)$F_stats_list)
stopifnot(length(Fref) == n_features)
rel_dev <- max(abs(Fall[, d0] - Fref) / pmax(abs(Fref), 1e-12))
log_diag("fast pseudo-F vs reference routine at d0, max rel deviation: ",
         sprintf("%.3e", rel_dev))
stopifnot(rel_dev < 1e-6)

Fmat <- Fall[, D_GRID, drop = FALSE]
saveRDS(Fmat, file.path(OUT_DIR, "pseudoF_matrix.rds"))

## 8. Metrics relative to d0

D_GRID <- as.integer(sub("^d", "", colnames(Fmat)))
col_d0 <- which(D_GRID == d0)
stopifnot(length(col_d0) == 1L)
message("d0 = ", d0, " -> column ", col_d0, " (", colnames(Fmat)[col_d0], ")")

n_bad <- sum(!is.finite(Fmat))
message("non-finite pseudo-F entries: ", n_bad,
        " (", sprintf("%.3f", 100 * n_bad / length(Fmat)), "%)")
bad_rows <- which(apply(Fmat, 1, function(v) any(!is.finite(v))))
message("features with any non-finite value: ", length(bad_rows))
KEEP <- setdiff(seq_len(nrow(Fmat)), bad_rows)
message("features retained for all metrics: ", length(KEEP))

jaccard_topk <- function(a, b, K) {
  ta <- logical(length(a)); ta[order(a, decreasing = TRUE)[1:K]] <- TRUE
  tb <- logical(length(b)); tb[order(b, decreasing = TRUE)[1:K]] <- TRUE
  sum(ta & tb) / sum(ta | tb)
}
stopifnot(abs(jaccard_topk(Fmat[KEEP, col_d0], Fmat[KEEP, col_d0], 100) - 1) < 1e-12)

target_idx_full <- which(feature_names == TARGET_FEATURE)
stopifnot(length(target_idx_full) == 1L)
target_idx <- match(target_idx_full, KEEP)   # position within KEEP
stopifnot(!is.na(target_idx))

summarise_sweep <- function(M, tag) {
  M   <- M[KEEP, , drop = FALSE]
  ref <- M[, col_d0]
  nf  <- nrow(M)
  
  out <- data.frame(
    variant        = rep(tag, ncol(M)),
    d              = D_GRID,
    d_over_d0      = D_GRID / d0,
    spearman_vs_d0 = vapply(seq_len(ncol(M)),
                            function(k) stats::cor(M[, k], ref, method = "spearman"),
                            numeric(1)),
    stringsAsFactors = FALSE
  )
  for (K in TOPK) {
    out[[paste0("jaccard_top", K)]] <-
      vapply(seq_len(ncol(M)), function(k) jaccard_topk(M[, k], ref, K), numeric(1))
  }
  out$target_rank <- vapply(seq_len(ncol(M)),
                            function(k) rank(-M[, k], ties.method = "min")[target_idx],
                            numeric(1))
  out$target_pctile <- 1 - (out$target_rank - 1) / nf
  out$target_F      <- M[target_idx, ]
  out
}

sweep_df <- summarise_sweep(Fmat, "weighted")

stopifnot("spearman_vs_d0" %in% names(sweep_df),
          !all(is.na(sweep_df$spearman_vs_d0)),
          abs(sweep_df$spearman_vs_d0[sweep_df$d == d0] - 1) < 1e-12,
          abs(sweep_df$jaccard_top100[sweep_df$d == d0] - 1) < 1e-12)

write.csv(sweep_df, file.path(OUT_DIR, "sweep_metrics.csv"), row.names = FALSE)


## 9. Unweighted basis W = V (diagnostic)

if (exists("SST") && exists("SSW")) {
  D_MAX <- ncol(SST)
  lam   <- fit$V_values[1:D_MAX]
  Ucum  <- upper.tri(matrix(0, D_MAX, D_MAX), diag = TRUE) * 1
  cumTu <- sweep(SST, 2, lam^2, "/") %*% Ucum
  cumWu <- sweep(SSW, 2, lam^2, "/") %*% Ucum
  Fall_u <- ((cumTu - cumWu) / (G - 1)) / (cumWu / (S - G))
  colnames(Fall_u) <- paste0("d", seq_len(D_MAX))
  Fmat_u <- Fall_u[, D_GRID, drop = FALSE]
  sweep_u <- summarise_sweep(Fmat_u, "unweighted")
  write.csv(sweep_u, file.path(OUT_DIR, "sweep_metrics_unweighted.csv"),
            row.names = FALSE)
  saveRDS(Fmat_u, file.path(OUT_DIR, "pseudoF_matrix_unweighted.rds"))
} else {
  message("SST/SSW not in memory -- skipping unweighted diagnostic")
}


## 10. Summary

diag_lines <- c()   # reset

log_diag("=== sweep summary (weighted) ===")
log_diag("d0 = ", d0, ";  grid = ", min(D_GRID), "..", max(D_GRID),
         " (", length(D_GRID), " points)")
log_diag("features retained: ", length(KEEP), " of ", nrow(Fmat))

amb    <- 10:17     # 95%-of-peak Kneedle ambiguity interval
in_amb <- sweep_df$d %in% amb
gv <- function(col, dd) sweep_df[[col]][sweep_df$d == dd]

log_diag("Spearman at d=12 (criterion runner-up): ", sprintf("%.4f", gv("spearman_vs_d0", 12)))
log_diag("min Spearman over ambiguity interval [10,17]: ",
         sprintf("%.4f", min(sweep_df$spearman_vs_d0[in_amb])))
log_diag("min Spearman over d >  d0: ",
         sprintf("%.4f", min(sweep_df$spearman_vs_d0[sweep_df$d >  d0])))
log_diag("min Spearman over d <  d0: ",
         sprintf("%.4f", min(sweep_df$spearman_vs_d0[sweep_df$d <  d0])))

log_diag("--- low end, where degradation is expected ---")
for (dd in 2:8) {
  log_diag("  d=", dd, "  spearman=", sprintf("%.4f", gv("spearman_vs_d0", dd)),
           "  J300=", sprintf("%.4f", gv("jaccard_top300", dd)),
           "  STAT5B rank=", gv("target_rank", dd))
}

drop <- sweep_df$d[sweep_df$d < d0 & sweep_df$spearman_vs_d0 < 0.95]
log_diag("largest d below d0 with Spearman < 0.95: ",
         if (length(drop)) max(drop) else "none")
drop9 <- sweep_df$d[sweep_df$d < d0 & sweep_df$spearman_vs_d0 < 0.90]
log_diag("largest d below d0 with Spearman < 0.90: ",
         if (length(drop9)) max(drop9) else "none")

log_diag("top-300 Jaccard at d=12: ", sprintf("%.4f", gv("jaccard_top300", 12)))
log_diag("min top-300 Jaccard over [10,17]: ",
         sprintf("%.4f", min(sweep_df$jaccard_top300[in_amb])))
log_diag("STAT5B rank at d=2 / 5 / 12 / d0 / 30 / 50: ",
         paste(sapply(c(2, 5, 12, d0, 30, 50), function(dd) gv("target_rank", dd)),
               collapse = " / "))

writeLines(diag_lines, file.path(OUT_DIR, "diagnostics.txt"))
print(sweep_df[sweep_df$d %in% c(2:6, 8, 10, 12, d0, 19, 22, 30, 40, 50), ],
      row.names = FALSE)

## 11. Figure and supporting numbers

suppressPackageStartupMessages(library(ggplot2))

PANEL_DIR <- file.path(OUT_DIR, "figures")   # individual panels; FIG_DIR is the final figure
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

G_ <- length(unique(groups)); S_ <- n_sample
N_AVAIL <- ncol(E_full[[1]])          # eigenpairs extracted (100)
D_HI    <- 99                         # largest usable dimension (lambda_100 ~ 0)

## 11.1 Extend SST / SSW to d = 1 .. D_HI and rebuild the sweep table

if (ncol(SST) < D_HI) {
  message("extending SST/SSW from ", ncol(SST), " to ", D_HI, " coordinates")
  grp_idx <- split(seq_len(S_), as.factor(groups))
  add <- (ncol(SST) + 1):D_HI
  SSTa <- matrix(NA_real_, n_features, length(add))
  SSWa <- matrix(NA_real_, n_features, length(add))
  for (k in seq_along(add)) {
    j <- add[k]
    X <- vapply(E_full, function(M) M[, j], numeric(n_features))
    SSTa[, k] <- rowSums((X - rowMeans(X))^2)
    w <- numeric(n_features)
    for (ix in grp_idx) {
      Xg <- X[, ix, drop = FALSE]; w <- w + rowSums((Xg - rowMeans(Xg))^2)
    }
    SSWa[, k] <- w
  }
  SST <- cbind(SST, SSTa); SSW <- cbind(SSW, SSWa)
}

cumT_all <- t(apply(SST[, 1:D_HI, drop = FALSE], 1, cumsum))
cumW_all <- t(apply(SSW[, 1:D_HI, drop = FALSE], 1, cumsum))
Fall <- ((cumT_all - cumW_all) / (G_ - 1)) / (cumW_all / (S_ - G_))
colnames(Fall) <- paste0("d", 1:D_HI)

# check against Fmat at d0
stopifnot(max(abs(Fall[, d0] - Fmat[, col_d0])) < 1e-8)

ref     <- Fall[, d0]
rank_d0 <- rank(-ref, ties.method = "min")
top500  <- order(ref, decreasing = TRUE)[1:500]

sweep_plot <- do.call(rbind, lapply(1:D_HI, function(d) {
  Fv <- Fall[, d]
  rk <- rank(-Fv, ties.method = "min")
  row <- data.frame(
    d = d, d_over_d0 = d / d0,
    spearman_vs_d0 = stats::cor(Fv, ref, method = "spearman"),
    stringsAsFactors = FALSE)
  for (K in TOPK) row[[paste0("jaccard_top", K)]] <- jaccard_topk(Fv, ref, K)
  # largest rank shift among the top-500 features
  row$max_rank_shift_top500 <- max(abs(rk[top500] - rank_d0[top500]))
  row$target_rank <- rk[target_idx_full]
  row
}))
write.csv(sweep_plot, file.path(OUT_DIR, "sweep_metrics_full.csv"),
          row.names = FALSE)

## 11.2 Neighbourhood overlap (consensus embedding, K = 20) and its robustness

NN_NSUB <- 2000
NN_MAIN_K <- 20
NN_GRID <- sort(unique(c(1:50, seq(55, D_HI, by = 5), D_HI, d0)))
NN_DIAG_GRID <- c(1, 2, 3, 5, 8, 10, 12, d0, 20, 30, 50, D_HI)

Ebar <- Reduce("+", E_full) / length(E_full)

nn_overlap <- function(M, sub, d, K, ref_nn) {
  D <- as.matrix(dist(M[sub, 1:d, drop = FALSE])); diag(D) <- Inf
  nn <- t(apply(D, 1, function(r) order(r)[1:K]))
  if (is.null(ref_nn)) return(nn)
  mean(vapply(seq_len(nrow(nn)), function(i) sum(nn[i, ] %in% ref_nn[i, ]) / K,
              numeric(1)))
}

NN_CACHE <- file.path(OUT_DIR, "nn_overlap.csv")
if (file.exists(NN_CACHE)) {
  nn_df <- read.csv(NN_CACHE)
} else {
  set.seed(1); sub1 <- sort(sample(n_features, NN_NSUB))
  ref_nn <- nn_overlap(Ebar, sub1, d0, NN_MAIN_K, NULL)
  nn_df <- do.call(rbind, lapply(NN_GRID, function(d) {
    data.frame(d = d,
               nn_overlap = nn_overlap(Ebar, sub1, d, NN_MAIN_K, ref_nn))
  }))
  stopifnot(abs(nn_df$nn_overlap[nn_df$d == d0] - 1) < 1e-12)
  write.csv(nn_df, NN_CACHE, row.names = FALSE)
}

NN_ROB <- file.path(OUT_DIR, "nn_robustness.csv")
if (file.exists(NN_ROB)) {
  nn_rob <- read.csv(NN_ROB)
} else {
  rob <- list()
  # (a) K sensitivity, consensus embedding, seed 1
  set.seed(1); sub1 <- sort(sample(n_features, NN_NSUB))
  for (K in c(10, 20, 50)) {
    rn <- nn_overlap(Ebar, sub1, d0, K, NULL)
    for (d in NN_DIAG_GRID)
      rob[[length(rob) + 1]] <- data.frame(
        setting = paste0("consensus, K=", K), d = d,
        overlap = nn_overlap(Ebar, sub1, d, K, rn))
  }
  # (b) different random feature subsets, K = 20
  for (sd in c(2, 3)) {
    set.seed(sd); sb <- sort(sample(n_features, NN_NSUB))
    rn <- nn_overlap(Ebar, sb, d0, NN_MAIN_K, NULL)
    for (d in NN_DIAG_GRID)
      rob[[length(rob) + 1]] <- data.frame(
        setting = paste0("consensus, seed=", sd), d = d,
        overlap = nn_overlap(Ebar, sb, d, NN_MAIN_K, rn))
  }
  # (c) per-sample embeddings rather than the consensus
  set.seed(1); sub1 <- sort(sample(n_features, NN_NSUB))
  ref_ps <- lapply(E_full, function(M) nn_overlap(M, sub1, d0, NN_MAIN_K, NULL))
  for (d in NN_DIAG_GRID) {
    ov <- vapply(seq_along(E_full), function(i)
      nn_overlap(E_full[[i]], sub1, d, NN_MAIN_K, ref_ps[[i]]), numeric(1))
    rob[[length(rob) + 1]] <- data.frame(
      setting = "per-sample (mean over 16)", d = d, overlap = mean(ov))
  }
  nn_rob <- do.call(rbind, rob)
  write.csv(nn_rob, NN_ROB, row.names = FALSE)
}

## 11.3 Ambiguity interval and pairwise agreement within it

spec  <- read.csv(file.path(OUT_DIR, "spectrum_d0.csv"))
specA <- spec[spec$index <= 50, ]
kd    <- specA$kneedle_distance
peak  <- max(kd, na.rm = TRUE)
amb   <- specA$index[!is.na(kd) & kd >= 0.95 * peak]
AMB_LO <- min(amb); AMB_HI <- max(amb)
runner_up <- specA$index[!is.na(kd)][order(-kd[!is.na(kd)])][2]

pw <- expand.grid(a = AMB_LO:AMB_HI, b = AMB_LO:AMB_HI)
pw <- pw[pw$a < pw$b, ]
pw$spearman <- mapply(function(a, b)
  stats::cor(Fall[, a], Fall[, b], method = "spearman"), pw$a, pw$b)
pw$jacc100  <- mapply(function(a, b)
  jaccard_topk(Fall[, a], Fall[, b], 100), pw$a, pw$b)
write.csv(pw, file.path(OUT_DIR, "pairwise_within_ambiguity.csv"),
          row.names = FALSE)

## 11.4 Per-coordinate signal profile

bw_ratio <- (SST[, 1:D_HI] - SSW[, 1:D_HI]) / pmax(SSW[, 1:D_HI], 1e-300)
coord_df <- data.frame(
  j = 1:D_HI,
  median_BW = apply(bw_ratio, 2, median),
  share_SST = colSums(SST[, 1:D_HI]) / sum(SST[, 1:D_HI])
)
coord_df$cum_share_SST <- cumsum(coord_df$share_SST)
write.csv(coord_df, file.path(OUT_DIR, "per_coordinate_profile.csv"),
          row.names = FALSE)

## 11.5 Theme

theme_s <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
          panel.border     = element_rect(linewidth = 0.4, colour = "black", fill = NA),
          axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
          axis.text        = element_text(colour = "black"),
          legend.key       = element_blank(),
          legend.background = element_blank(),
          legend.title     = element_blank(),
          legend.key.size  = unit(8, "pt"),
          legend.text      = element_text(size = base_size - 1.5),
          legend.margin    = margin(0, 0, 0, 0),
          plot.title       = element_text(face = "bold", size = base_size + 2,
                                          hjust = 0, margin = margin(b = 3)))
}

amb_band <- annotate("rect", xmin = AMB_LO, xmax = AMB_HI,
                     ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.45)
d0_line  <- geom_vline(xintercept = d0, linetype = "dashed",
                       linewidth = 0.4, colour = "grey20")
x_lin <- scale_x_continuous(limits = c(1, 50),
                            breaks = c(1, 10, 15, 20, 30, 40, 50),
                            expand = expansion(mult = 0.015))

## Panel A: spectrum

pA <- ggplot(specA, aes(index, log_eigenvalue)) +
  amb_band + d0_line +
  geom_line(linewidth = 0.5, colour = "black") +
  geom_point(size = 0.7, colour = "black") +
  x_lin + scale_y_continuous(expand = expansion(mult = 0.06)) +
  labs(x = "dimension index j", y = expression(log~lambda[j]), title = "A") +
  theme_s() + theme(legend.position = "none")

## Panel B: Kneedle criterion curve

stopifnot(d0 == 15, runner_up == 12)
lab_df <- data.frame(
  index = c(d0, runner_up),
  y     = c(kd[d0], kd[runner_up]),
  lab   = c("selected~(italic(d)[0]==15)", "runner*'-'*up~(italic(d)==12)"),
  hj    = c(-0.12, 1.12))

pB <- ggplot(specA, aes(index, kneedle_distance)) +
  amb_band + d0_line +
  geom_line(linewidth = 0.5, colour = "black") +
  geom_point(size = 0.7, colour = "black") +
  geom_point(data = specA[specA$index == runner_up, ], size = 1.7,
             shape = 21, fill = "white", colour = "black", stroke = 0.45) +
  geom_point(data = specA[specA$index == d0, ], size = 2.0,
             shape = 21, fill = "#D55E00", colour = "black", stroke = 0.45) +
  geom_text(data = lab_df, aes(index, y, label = lab, hjust = hj),
            parse = TRUE, size = 2.6, vjust = 0, inherit.aes = FALSE) +
  x_lin + scale_y_continuous(expand = expansion(mult = c(0.12, 0.08))) +
  labs(x = "dimension index j",
       y = "normalised perpendicular distance", title = "B") +
  theme_s() + theme(legend.position = "none")

## Panel C: agreement with d0

lev <- c("pseudo-F, Spearman (all features)", "pseudo-F, Jaccard top 100",
         "pseudo-F, Jaccard top 300", "pseudo-F, Jaccard top 500",
         paste0("embedding geometry, ", NN_MAIN_K, "-NN overlap"))
longC <- rbind(
  data.frame(d = sweep_plot$d, value = sweep_plot$spearman_vs_d0, metric = lev[1]),
  data.frame(d = sweep_plot$d, value = sweep_plot$jaccard_top100,  metric = lev[2]),
  data.frame(d = sweep_plot$d, value = sweep_plot$jaccard_top300,  metric = lev[3]),
  data.frame(d = sweep_plot$d, value = sweep_plot$jaccard_top500,  metric = lev[4]),
  data.frame(d = nn_df$d,      value = nn_df$nn_overlap,           metric = lev[5]))
longC$metric <- factor(longC$metric, levels = lev)
COLS <- setNames(c("#000000", "#0072B2", "#D55E00", "#009E73", "#CC79A7"), lev)

pC <- ggplot(longC, aes(d, value, colour = metric, shape = metric)) +
  amb_band + d0_line +
  geom_line(linewidth = 0.45) + geom_point(size = 1.0) +
  scale_colour_manual(values = COLS) +
  scale_shape_manual(values = c(16, 17, 15, 18, 1)) +
  scale_x_continuous(trans = "log2", limits = c(1, D_HI),
                     breaks = c(1, 2, 4, 8, 15, 30, 50, D_HI),
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(0.25, 1.02), breaks = seq(0.25, 1, 0.25)) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  labs(x = expression("enforced dimension "*italic(d)*"  (log"[2]*" scale)"),
       y = expression("agreement with "*italic(d)[0]), title = "C") +
  theme_s() + theme(legend.position = "bottom",
                    legend.box.margin = margin(t = -4))

sat_d <- min(nn_df$d[nn_df$nn_overlap >= 0.99])


## Assemble

save_fig <- function(p, file, w, h) {
  ggsave(paste0(file, ".pdf"), p, width = w, height = h, units = "in",
         device = cairo_pdf)
  ggsave(paste0(file, ".png"), p, width = w, height = h, units = "in", dpi = 400)
}
for (nm in c("A", "B", "C"))
  save_fig(get(paste0("p", nm)), file.path(PANEL_DIR, paste0("figS13_panel", nm)),
           3.4, 2.6)

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  fig <- (pA | pB) / pC + patchwork::plot_layout(heights = c(1, 1.25))
  save_fig(fig, file.path(PANEL_DIR, "figS13"), 6.9, 5.8)
} else if (requireNamespace("cowplot", quietly = TRUE)) {
  top <- cowplot::plot_grid(pA, pB, ncol = 2, align = "h")
  fig <- cowplot::plot_grid(top, pC, ncol = 1, rel_heights = c(1, 1.25))
  save_fig(fig, file.path(PANEL_DIR, "figS13"), 6.9, 5.8)
} else if (requireNamespace("gridExtra", quietly = TRUE)) {
  fig <- gridExtra::arrangeGrob(pA, pB, pC, layout_matrix = rbind(c(1, 2), c(3, 3)))
  ggsave(file.path(PANEL_DIR, "figS13.pdf"), fig,
         width = 6.9, height = 5.8)
}

## Supplemental Figure S13
png(file.path(FIG_DIR, "figS13.png"),width = 6.9, height = 6.9, units = "in", res = 600)
print(fig)
dev.off()


## 11.6 Numbers

gv  <- function(col, dd) sweep_plot[[col]][sweep_plot$d == dd]
gnn <- function(dd) nn_df$nn_overlap[nn_df$d == dd]
ev  <- specA$eigenvalue
lr  <- log(ev[-length(ev)] / ev[-1])
d2h <- sweep_plot$d >= 2

cat("\n================ T-cell numbers ================\n")

cat("\n[1] Spectrum\n")
cat("  criterion peak at d0 =", d0, "; runner-up d =", runner_up,
    sprintf("(only %.2f%% below peak)\n", 100 * (peak - kd[runner_up]) / peak))
cat("  95%-of-peak interval = [", AMB_LO, ",", AMB_HI, "]\n")
cat("  for j >= 4, consecutive eigenvalues differ by at most",
    sprintf("%.0f%%", 100 * (exp(max(lr[4:49])) - 1)), "\n")
cat("  lambda_15 / lambda_50 =", sprintf("%.2f", ev[15] / ev[50]),
    " (a factor of <2 across 35 indices)\n")

cat("\n[2] DC ranking vs d\n")
cat("  d = 2 ..", D_HI, ": Spearman >=",
    sprintf("%.4f", min(sweep_plot$spearman_vs_d0[d2h])),
    "; top-100 Jaccard >=", sprintf("%.4f", min(sweep_plot$jaccard_top100[d2h])), "\n")
cat("  max rank shift among top-500 features, over d = 2 ..", D_HI, ":",
    max(sweep_plot$max_rank_shift_top500[d2h]), "\n")
cat("  eyeball-knee check: top-100 Jaccard at d=3 /", sprintf("%.4f", gv("jaccard_top100", 3)),
    " d=4 /", sprintf("%.4f", gv("jaccard_top100", 4)), "\n")
cat("  NEGATIVE CONTROL d=1: Spearman", sprintf("%.4f", gv("spearman_vs_d0", 1)),
    "; top-100 Jaccard", sprintf("%.4f", gv("jaccard_top100", 1)),
    "; max top-500 rank shift", gv("max_rank_shift_top500", 1), "\n")

cat("\n[3] Pairwise agreement within [", AMB_LO, ",", AMB_HI, "]\n")
cat("  min pairwise Spearman:", sprintf("%.4f", min(pw$spearman)),
    "; min pairwise top-100 Jaccard:", sprintf("%.4f", min(pw$jacc100)), "\n")

cat("\n[4] Neighbourhood overlap vs d\n")
cat("  20-NN overlap: d=1", sprintf("%.3f", gnn(1)), " d=2", sprintf("%.3f", gnn(2)),
    " d=5", sprintf("%.3f", gnn(5)), " d=8", sprintf("%.3f", gnn(8)),
    " d=10", sprintf("%.3f", gnn(10)), " d=50", sprintf("%.3f", gnn(50)), "\n")
cat("  geometry saturates at d =", min(nn_df$d[nn_df$nn_overlap >= 0.99]),
    "(overlap >= 0.99); d0 =", d0, "sits above it with margin\n")
cat("  robustness of that saturation point:\n")
for (s in unique(nn_rob$setting)) {
  z <- nn_rob[nn_rob$setting == s, ]
  cat("    ", sprintf("%-28s", s), " saturation d =",
      min(z$d[z$overlap >= 0.99]), "; overlap at d=2 =",
      sprintf("%.3f", z$overlap[z$d == 2]), "\n")
}

cat("\n[5] Per-coordinate share of SS\n")
cat("  share of total SS carried by coordinate 1 /1-2 /1-8:",
    sprintf("%.3f", coord_df$cum_share_SST[1]), "/",
    sprintf("%.3f", coord_df$cum_share_SST[2]), "/",
    sprintf("%.3f", coord_df$cum_share_SST[8]), "\n")
cat("  cumulative squared spectral weight at d=2 / d0 / 2*d0:",
    sprintf("%.4f", spec$cum_sq_spectral_weight[2]), "/",
    sprintf("%.4f", spec$cum_sq_spectral_weight[d0]), "/",
    sprintf("%.4f", spec$cum_sq_spectral_weight[2 * d0]), "\n")

cat("\n[6] STAT5B rank\n")
rr <- sweep_plot$target_rank[d2h]
cat("  rank range over d = 2 ..", D_HI, ":", min(rr), "-", max(rr),
    "of", n_features, sprintf("(%.2f%%-%.2f%% percentile)\n",
                              100 * min(rr) / n_features, 100 * max(rr) / n_features))
cat("  rank at d0:", gv("target_rank", d0), "; at d=1:", gv("target_rank", 1), "\n")

cat("\n[7] Checks\n")
cat("  column nesting exact to < 1e-8; fast pseudo-F reproduces the\n")
cat("  PERMANOVA routine at d0 to < 1e-6 relative; sweep spans d = 1 ..", D_HI,
    "\n")
