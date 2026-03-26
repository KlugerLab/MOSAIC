# Supplemental Figure S12: composition shift versus within-population analysis (simulation)
# Output:   results/S12/results_rep<k>.rds, primary_fpr_summary.csv, paired_difference_summary.csv, run_metadata.rds
# Usage:    OMP_NUM_THREADS=1 Rscript fig3/figS12_composition_simulation.R --set n_cores=65

source("config.R")

## Arms: M = all cells, composition shift only; L = same data, focal population, fixed n_sub; P = L plus within-population rewiring

suppressPackageStartupMessages({
  library(RSpectra)
  library(foreach)
  library(doParallel)
  have_vegan <- requireNamespace("vegan", quietly = TRUE)
  have_prroc <- requireNamespace("PRROC", quietly = TRUE)
})

## 1. Config
cfg <- list(
  ## generator
  n_samples        = 20L,
  cells_per_sample = 1000L,
  cell_jitter      = 200L,        # n_cells = 1000 + U{-200..200}; cond-independent
  r_latent         = 10L,
  px = 1000L, py = 800L, pz = 600L,
  n_clusters       = 10L,
  cell_center_sd   = 0.1,
  feat_center_sd   = 0.1,
  sample_noise_sd  = 0.5,
  dia_type         = "medium",
  add_batch        = FALSE,

  ## composition severity (pi of the focal population)
  focal_cluster    = 1L,
  pA_focal         = 0.10,        # condition A: uniform 1/n_clusters
  severities       = c(mild = 0.15, moderate = 0.25, strong = 0.40),
  severity_for_P   = "moderate",

  ## fixed subsample size for arms L and P
  n_sub            = 75L,

  ## DC signal (Arm P only)
  dc_prop          = 0.2,
  rescale_scope    = "population",  # "population" (default) | "mixture" (check)

  ## MOSAIC pipeline
  eig_window       = 50L,
  scale_clip       = 10,

  ## inference
  alpha_fpr        = 0.05,
  alpha_fpr_extra  = 0.01,
  alpha_bh         = 0.05,
  perm_chunk       = 4000L,

  ## replicates / execution
  n_reps           = 30L,
  seed_base        = 2100L,
  n_cores          = 1L,
  out_dir          = file.path(RESULTS_DIR, "S12")
)

dia_matrix_of <- function(type, r) switch(type,
  strong = diag(seq(10, 5,   length.out = r)),
  medium = diag(seq(5,  1,   length.out = r)),
  weak   = diag(seq(3,  0.5, length.out = r)),
  stop("unknown dia_type"))

## 2. Helpers from 0-functions.R
perpendicular_distance <- function(xp, yp, x1, y1, x2, y2) {
  A <- y1 - y2; B <- x2 - x1; C <- -A * x1 - B * y1
  abs(A * xp + B * yp + C) / sqrt(A^2 + B^2)
}

find_elbow_kneedle <- function(eigenvalues) {
  for (i in seq_along(eigenvalues))
    if (eigenvalues[i] <= 0) eigenvalues[i] <- eigenvalues[i - 1]
  x <- seq_along(eigenvalues)
  y <- log(eigenvalues)
  if (length(x) < 3) return(length(x))
  x_scaled <- (x - min(x)) / (max(x) - min(x))
  y_scaled <- (y - min(y)) / (max(y) - min(y))
  distances <- perpendicular_distance(x_scaled, y_scaled,
                                      x_scaled[1], y_scaled[1],
                                      x_scaled[length(x_scaled)],
                                      y_scaled[length(y_scaled)])
  max(which.max(distances), 2)
}

## 3. Generator

## deterministic apportionment of n cells to clusters at proportions p
allocate_counts <- function(n, p) {
  raw <- n * p
  cnt <- floor(raw)
  rem <- n - sum(cnt)
  if (rem > 0) {
    ord <- order(raw - cnt, decreasing = TRUE)
    cnt[ord[seq_len(rem)]] <- cnt[ord[seq_len(rem)]] + 1L
  }
  stopifnot(sum(cnt) == n)
  as.integer(cnt)
}

generate_composition_data <- function(seed, signal = c("none", "DC"),
                                      pB_focal, cfg) {
  signal <- match.arg(signal)
  set.seed(seed)
  S <- cfg$n_samples; r <- cfg$r_latent; C <- cfg$n_clusters
  Sx <- Sy <- Sz <- dia_matrix_of(cfg$dia_type, r)

  meta <- data.frame(sample_id = paste0("sample_", seq_len(S)),
                     condition = rep(c("A", "B"), each = S / 2),
                     stringsAsFactors = FALSE)

  probs_A <- rep((1 - cfg$pA_focal) / (C - 1), C); probs_A[cfg$focal_cluster] <- cfg$pA_focal
  probs_B <- rep((1 - pB_focal)     / (C - 1), C); probs_B[cfg$focal_cluster] <- pB_focal
  stopifnot(abs(sum(probs_A) - 1) < 1e-12, abs(sum(probs_B) - 1) < 1e-12)

  centers_U <- matrix(rnorm(C * r), nrow = C)
  U_list <- vector("list", S); cell_clusters <- vector("list", S)
  for (j in seq_len(S)) {
    n_cells <- cfg$cells_per_sample + sample(-cfg$cell_jitter:cfg$cell_jitter, 1)
    pr <- if (meta$condition[j] == "A") probs_A else probs_B
    cl <- rep.int(seq_len(C), allocate_counts(n_cells, pr))   # deterministic allocation
    U  <- centers_U[cl, , drop = FALSE] +
          matrix(rnorm(n_cells * r, sd = cfg$cell_center_sd), n_cells, r)
    U_list[[j]] <- U; cell_clusters[[j]] <- cl
  }

  make_global_V <- function(p) {
    cl <- sample.int(C, p, replace = TRUE)
    list(V = centers_U[cl, , drop = FALSE] +
           matrix(rnorm(p * r, sd = cfg$feat_center_sd), p, r), cl = cl)
  }
  gx <- make_global_V(cfg$px); gy <- make_global_V(cfg$py); gz <- make_global_V(cfg$pz)
  samp_V <- function(Vg) lapply(seq_len(S), function(i)
    Vg + matrix(rnorm(length(Vg), sd = cfg$sample_noise_sd), nrow(Vg), ncol(Vg)))
  Vx <- samp_V(gx$V); Vy <- samp_V(gy$V); Vz <- samp_V(gz$V)

  dc_x <- dc_y <- dc_z <- integer(0)
  Vx_use <- Vx; Vy_use <- Vy; Vz_use <- Vz
  if (signal == "DC") {
    constrained_shuffle <- function(idx, cl) {
      out <- idx
      for (i in seq_along(idx)) {
        pool <- idx[cl[idx] != cl[idx[i]]]
        if (length(pool) == 0) pool <- idx[idx != idx[i]]
        out[i] <- if (length(pool) == 1) pool else sample(pool, 1)
      }
      out
    }
    dc_x <- sample.int(cfg$px, round(cfg$px * cfg$dc_prop))
    dc_y <- sample.int(cfg$py, round(cfg$py * cfg$dc_prop))
    dc_z <- sample.int(cfg$pz, round(cfg$pz * cfg$dc_prop))
    sx <- constrained_shuffle(dc_x, gx$cl); sy <- constrained_shuffle(dc_y, gy$cl)
    sz <- constrained_shuffle(dc_z, gz$cl)
    for (i in seq_len(S)) if (meta$condition[i] == "B") {
      Vx_use[[i]][dc_x, ] <- Vx[[i]][sx, ]
      Vy_use[[i]][dc_y, ] <- Vy[[i]][sy, ]
      Vz_use[[i]][dc_z, ] <- Vz[[i]][sz, ]
    }
  }

  gen_mod <- function(Vu, Sm) lapply(seq_len(S), function(i) U_list[[i]] %*% Sm %*% t(Vu[[i]]))
  X_list <- gen_mod(Vx_use, Sx); Y_list <- gen_mod(Vy_use, Sy); Z_list <- gen_mod(Vz_use, Sz)

  ## pure-DC mean matching on focal-population cells (or all cells)
  if (signal == "DC") {
    ref_rows <- if (cfg$rescale_scope == "population")
      which(cell_clusters[[1]] == cfg$focal_cluster) else seq_len(nrow(U_list[[1]]))
    tgt <- function(Vg_list, Sm, dcidx)
      colMeans((U_list[[1]] %*% Sm %*% t(Vg_list[[1]]))[ref_rows, dcidx, drop = FALSE])
    tx <- tgt(Vx, Sx, dc_x); ty <- tgt(Vy, Sy, dc_y); tz <- tgt(Vz, Sz, dc_z)
    for (i in seq_len(S)) if (meta$condition[i] == "B") {
      rows_i <- if (cfg$rescale_scope == "population")
        which(cell_clusters[[i]] == cfg$focal_cluster) else seq_len(nrow(X_list[[i]]))
      shift_cols <- function(M, dcidx, t_) {
        cur <- colMeans(M[rows_i, dcidx, drop = FALSE])
        M[, dcidx] <- sweep(M[, dcidx, drop = FALSE], 2, t_ - cur, "+"); M
      }
      X_list[[i]] <- shift_cols(X_list[[i]], dc_x, tx)
      Y_list[[i]] <- shift_cols(Y_list[[i]], dc_y, ty)
      Z_list[[i]] <- shift_cols(Z_list[[i]], dc_z, tz)
    }
  }

  if (isTRUE(cfg$add_batch)) {
    nb <- 4L; meta$batch <- rep(paste0("batch_", 1:nb), each = S / nb)
    for (i in seq_len(S)) {
      b <- as.numeric(sub("batch_", "", meta$batch[i]))
      ap <- function(M, p) sweep(sweep(M, 2, rnorm(p, 1.2, 0.1), `*`), 2,
                                 rnorm(p, 10, 0.2) * (b / nb), `+`)
      X_list[[i]] <- ap(X_list[[i]], cfg$px); Y_list[[i]] <- ap(Y_list[[i]], cfg$py)
      Z_list[[i]] <- ap(Z_list[[i]], cfg$pz)
    }
  }

  dc_truth <- rep(FALSE, cfg$px + cfg$py + cfg$pz)
  if (signal == "DC")
    dc_truth[c(dc_x, cfg$px + dc_y, cfg$px + cfg$py + dc_z)] <- TRUE

  ## modality-concatenated matrices
  M_list <- vector("list", S)
  for (i in seq_len(S)) {
    M_list[[i]] <- cbind(X_list[[i]], Y_list[[i]], Z_list[[i]])
    X_list[i] <- list(NULL); Y_list[i] <- list(NULL); Z_list[i] <- list(NULL)
  }
  rm(X_list, Y_list, Z_list, Vx, Vy, Vz, Vx_use, Vy_use, Vz_use); gc(FALSE)

  list(M_list = M_list,
       meta = meta, cell_clusters = cell_clusters,
       feature_module = c(gx$cl, gy$cl, gz$cl), dc_truth = dc_truth,
       probs = list(A = probs_A, B = probs_B))
}

## 4. MOSAIC pipeline
## Pooled per-feature standardization over all cells (ScaleData analogue), streamed over samples
pooled_scale <- function(M_list, clip) {
  F_ <- ncol(M_list[[1]]); n <- sum(vapply(M_list, nrow, 1L))
  s1 <- numeric(F_)
  for (M in M_list) s1 <- s1 + colSums(M)
  mu <- s1 / n
  s2 <- numeric(F_)                                   # second pass: exact SS
  for (M in M_list) s2 <- s2 + colSums(sweep(M, 2, mu, "-")^2)
  sd_ <- sqrt(s2 / (n - 1))
  stopifnot(all(sd_ > 1e-12))                         # no constant features
  out <- lapply(M_list, function(M) {
    M <- sweep(M, 2, mu, "-"); M <- sweep(M, 2, sd_, "/")
    M[M > clip] <- clip                               # Seurat scale.max (upper)
    M
  })
  chk <- numeric(F_)
  for (M in out) chk <- chk + colSums(M)
  stopifnot(max(abs(chk / n)) < 1e-8)
  out
}

run_mosaic <- function(Msc_list, cfg) {
  S <- length(Msc_list); F_ <- ncol(Msc_list[[1]])
  k <- min(cfg$eig_window, F_ - 1L)

  ## column-normalize; zero norms set to 1, as in cosine_distance
  chat_of <- function(M) {
    nrm <- sqrt(colSums(M^2)); nrm[nrm == 0] <- 1
    sweep(M, 2, nrm, "/")
  }

  P_agg <- matrix(0, F_, F_); r_i <- integer(S)
  for (i in seq_len(S)) {
    U_i <- crossprod(chat_of(Msc_list[[i]]))  # full block cosine matrix
    e   <- RSpectra::eigs_sym(U_i, k, which = "LA")
    rm(U_i)
    ri  <- find_elbow_kneedle(e$values)       # Kneedle on all k values
    r_i[i] <- ri
    Vr  <- e$vectors[, seq_len(ri), drop = FALSE]
    P_agg <- P_agg + tcrossprod(sweep(Vr, 2, sqrt(pmax(e$values[seq_len(ri)], 0)), "*"))
    rm(e, Vr)
  }

  e2 <- RSpectra::eigs_sym(P_agg, k, which = "LA"); rm(P_agg); gc(FALSE)
  d  <- find_elbow_kneedle(e2$values[seq_len(k)])
  W  <- e2$vectors[, seq_len(d), drop = FALSE] %*% diag(e2$values[seq_len(d)])

  E <- array(NA_real_, dim = c(F_, d, S))
  for (i in seq_len(S)) {                                           # E_i = U_i W
    Ci <- chat_of(Msc_list[[i]])
    E[, , i] <- crossprod(Ci, Ci %*% W); rm(Ci)
  }

  list(E = E, d = d, r_i = r_i, W = W)
}

## 5. Exact permutation null (fixed embedding, all balanced label splits)
exact_null_test <- function(E, groups, chunk = 4000L) {
  F_ <- dim(E)[1]; d <- dim(E)[2]; S <- dim(E)[3]
  stopifnot(S %% 2L == 0L, length(groups) == S)
  gA <- groups == groups[1]
  stopifnot(sum(gA) == S / 2L)                        # balanced design required

  ## per-feature K_f = Yc_f Yc_f^T and SST_f
  Kmat <- matrix(0, F_, S * S); sst <- numeric(F_)
  for (f in seq_len(F_)) {
    Yf <- t(matrix(E[f, , ], nrow = d))               # S x d
    stopifnot(!any(rowSums(abs(Yf)) == 0))            # no all-zero rows
    Yc <- sweep(Yf, 2, colMeans(Yf), "-")
    Kmat[f, ] <- as.vector(tcrossprod(Yc)); sst[f] <- sum(Yc^2)
  }

  cSSB <- 4 / S                                       # SSB = (4/S) (g^T Yc)^2
  ssb_of <- function(g) cSSB * as.vector(Kmat %*% as.vector(tcrossprod(g)))
  ssb_obs <- ssb_of(as.numeric(gA))
  F_obs   <- (ssb_obs) / ((sst - ssb_obs) / (S - 2))

  ## one representative per complementary pair: subsets of size S/2 containing sample 1
  combs <- combn(S - 1L, S / 2L - 1L)                 # (S/2-1) x choose(S-1,S/2-1)
  nP    <- ncol(combs)
  ge    <- integer(F_)
  tol   <- 1e-9 * pmax(abs(ssb_obs), 1)
  ii    <- rep(seq_len(S), S); jj <- rep(seq_len(S), each = S)

  for (start in seq(1L, nP, by = chunk)) {
    sel <- start:min(start + chunk - 1L, nP)
    G <- matrix(0, length(sel), S); G[, 1L] <- 1
    idx <- combs[, sel, drop = FALSE] + 1L
    G[cbind(rep(seq_along(sel), each = nrow(idx)), as.vector(idx))] <- 1
    Wp   <- G[, ii, drop = FALSE] * G[, jj, drop = FALSE]           # vec(g g^T)
    SSB  <- cSSB * (Wp %*% t(Kmat))
    ge   <- ge + colSums(SSB >= rep(ssb_obs - tol, each = length(sel)))
  }

  p <- ge / nP
  list(F_obs = F_obs, p = p, p_bh = p.adjust(p, "BH"), n_perm_exact = nP)
}

assert_euF_matches_adonis2 <- function() {
  if (!have_vegan) { message("[check] vegan absent -- adonis2 cross-check SKIPPED"); return(invisible(FALSE)) }
  set.seed(1); S <- 12L; d <- 4L
  M <- matrix(rnorm(S * d), S, d); g <- rep(c("A", "B"), each = S / 2)
  E <- array(0, dim = c(1, d, S)); E[1, , ] <- t(M)
  ours <- exact_null_test(E, g)$F_obs
  ref  <- as.data.frame(vegan::adonis2(dist(M) ~ g, permutations = 2))$F[1]
  stopifnot(abs(ours - ref) < 1e-8)
  message(sprintf("[check] pseudo-F == adonis2: %.10f vs %.10f", ours, ref))
  ## enumeration size
  message(sprintf("[check] exact enumeration size for S=20: %d", choose(19, 9)))
  invisible(TRUE)
}

## 6. AUPRC
auprc <- function(scores, labels01) {
  stopifnot(length(scores) == length(labels01), any(labels01 == 1), any(labels01 == 0))
  if (have_prroc)
    return(PRROC::pr.curve(scores.class0 = scores, weights.class0 = labels01)$auc.integral)
  o <- order(scores, decreasing = TRUE)
  tp <- cumsum(labels01[o]); fp <- cumsum(1 - labels01[o])
  sum(diff(c(0, tp / sum(labels01))) * (tp / (tp + fp)))
}

## 7. Arms
subset_focal_fixed_n <- function(dat, cfg, seed_sub) {
  counts <- vapply(seq_along(dat$M_list), function(i)
    sum(dat$cell_clusters[[i]] == cfg$focal_cluster), 1L)
  if (min(counts) < cfg$n_sub)
    stop(sprintf("n_sub = %d exceeds the smallest focal count (%d)", cfg$n_sub, min(counts)))
  set.seed(seed_sub)
  lapply(seq_along(dat$M_list), function(i) {
    rows <- which(dat$cell_clusters[[i]] == cfg$focal_cluster)
    dat$M_list[[i]][sample(rows, cfg$n_sub), , drop = FALSE]   # fixed n_sub
  })
}

run_arm <- function(arm, dat, cfg, seed) {
  Msc <- if (arm == "M") pooled_scale(dat$M_list, cfg$scale_clip)
         else pooled_scale(subset_focal_fixed_n(dat, cfg, seed + 7777L), cfg$scale_clip)
  fit  <- run_mosaic(Msc, cfg)
  rm(Msc); gc(FALSE)                    # scaled cohort no longer needed
  tst  <- exact_null_test(fit$E, dat$meta$condition, cfg$perm_chunk)
  hits <- which(tst$p_bh < cfg$alpha_bh)
  out <- list(arm = arm, d = fit$d, median_ri = median(fit$r_i),
              n_perm_exact = tst$n_perm_exact,
              fpr05 = mean(tst$p < cfg$alpha_fpr),
              fpr01 = mean(tst$p < cfg$alpha_fpr_extra),
              n_bh_hits = length(hits), min_p = min(tst$p),
              p = tst$p, F_obs = tst$F_obs,
              hits_focal_module = sum(dat$feature_module[hits] == cfg$focal_cluster),
              share_focal_module = mean(dat$feature_module == cfg$focal_cluster))
  if (arm == "P") {
    out$auprc      <- auprc(tst$F_obs, as.integer(dat$dc_truth))
    out$hits_true  <- sum(dat$dc_truth[hits]); out$hits_false <- sum(!dat$dc_truth[hits])
  }
  out
}

run_replicate <- function(k, cfg) {
  seed <- cfg$seed_base + k
  res <- list(replicate = k, seed = seed, arms = list())
  for (sv in names(cfg$severities)) {
    dat <- generate_composition_data(seed, "none", cfg$severities[[sv]], cfg)
    res$arms[[paste0("M.", sv)]] <- c(run_arm("M", dat, cfg, seed), severity = sv)
    res$arms[[paste0("L.", sv)]] <- c(run_arm("L", dat, cfg, seed), severity = sv)
  }
  svP <- cfg$severity_for_P
  datP <- generate_composition_data(seed, "DC", cfg$severities[[svP]], cfg)
  res$arms[[paste0("P.", svP)]] <- c(run_arm("P", datP, cfg, seed), severity = svP)
  saveRDS(res, file.path(cfg$out_dir, sprintf("results_rep%02d.rds", k)))
  res
}

## 8. Composition-exposure diagnostic: Delta_S = sum_t (pi_t^B - pi_t^A) A_t, A_t from condition-A cells
composition_exposure_diag <- function(dat, cfg, Msc, fit, F_obs) {
  condA <- which(dat$meta$condition == "A"); C <- cfg$n_clusters
  F_ <- ncol(Msc[[1]])
  dpi <- dat$probs$B - dat$probs$A
  DeltaS <- matrix(0, F_, F_)
  for (t in seq_len(C)) {
    if (abs(dpi[t]) < 1e-12) next
    Zt <- do.call(rbind, lapply(condA, function(i)
      Msc[[i]][dat$cell_clusters[[i]] == t, , drop = FALSE]))
    DeltaS <- DeltaS + dpi[t] * (crossprod(Zt) / nrow(Zt))   # A_t = E[ZZ^T|T=t]
  }
  Q <- qr.Q(qr(fit$W))
  list(cor_full  = cor(sqrt(rowSums(DeltaS^2)), F_obs, method = "spearman"),
       cor_projW = cor(sqrt(rowSums((DeltaS %*% Q)^2)), F_obs, method = "spearman"),
       energy_in_W = norm(crossprod(Q, DeltaS) %*% Q, "F") / norm(DeltaS, "F"),
       frob = norm(DeltaS, "F"))
}

## 9. Main
main <- function(cfg) {
  dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
  assert_euF_matches_adonis2()

  t_start <- Sys.time()
  reps <- if (cfg$n_cores > 1L) {
    ## fork cluster on unix
    cl <- if (.Platform$OS.type == "unix") parallel::makeForkCluster(cfg$n_cores)
          else makeCluster(cfg$n_cores)
    registerDoParallel(cl); on.exit(stopCluster(cl), add = TRUE)
    foreach(k = seq_len(cfg$n_reps), .packages = "RSpectra",
            .export = setdiff(ls(envir = globalenv()), "cfg")) %dopar%
      run_replicate(k, cfg)
  } else lapply(seq_len(cfg$n_reps), run_replicate, cfg = cfg)
  message(sprintf("[timing] %d replicates in %.1f min",
                  cfg$n_reps, as.numeric(Sys.time() - t_start, units = "mins")))

  summ <- do.call(rbind, lapply(reps, function(r)
    do.call(rbind, lapply(names(r$arms), function(nm) {
      a <- r$arms[[nm]]
      data.frame(replicate = r$replicate, arm = substr(nm, 1, 1), severity = a$severity,
                 pB_focal = cfg$severities[[a$severity]], d = a$d, median_ri = a$median_ri,
                 fpr05 = a$fpr05, fpr01 = a$fpr01, n_bh_hits = a$n_bh_hits,
                 min_p = a$min_p, n_perm_exact = a$n_perm_exact,
                 auprc = if (is.null(a$auprc)) NA_real_ else a$auprc,
                 hits_focal_module = a$hits_focal_module,
                 stringsAsFactors = FALSE)
    }))))
  write.csv(summ, file.path(cfg$out_dir, "summary_by_run.csv"), row.names = FALSE)

  ## Mean FPR at 0.05 with t-based CI across replicates
  agg <- do.call(rbind, lapply(split(summ, list(summ$arm, summ$severity), drop = TRUE), function(s) {
    n <- nrow(s); se <- sd(s$fpr05) / sqrt(n)
    tq <- if (n > 1L) qt(0.975, df = n - 1L) else NA_real_
    data.frame(arm = s$arm[1], severity = s$severity[1], pB_focal = s$pB_focal[1],
               n_rep = n, mean_fpr05 = mean(s$fpr05), sd_fpr05 = sd(s$fpr05),
               lo = mean(s$fpr05) - tq * se, hi = mean(s$fpr05) + tq * se,
               mean_fpr01 = mean(s$fpr01), median_bh_hits = median(s$n_bh_hits),
               mean_auprc = mean(s$auprc), stringsAsFactors = FALSE)
  }))
  agg <- agg[order(agg$arm, agg$pB_focal), ]
  write.csv(agg, file.path(cfg$out_dir, "primary_fpr_summary.csv"), row.names = FALSE)
  cat("\n=== PRIMARY: false-DC rate at nominal 0.05 (mean over", cfg$n_reps,
      "independent replicates; CI across replicates) ===\n")
  print(agg, digits = 3, row.names = FALSE)

  ## Paired difference FPR_M - FPR_L at matched severity
  wide <- reshape(summ[summ$arm %in% c("M", "L"),
                       c("replicate", "severity", "arm", "fpr05")],
                  direction = "wide", idvar = c("replicate", "severity"),
                  timevar = "arm")
  names(wide) <- sub("^fpr05\\.", "", names(wide))
  wide$delta <- wide$M - wide$L
  wide$pB_focal <- unlist(cfg$severities[wide$severity])
  wide <- wide[order(wide$pB_focal, wide$replicate), ]
  write.csv(wide, file.path(cfg$out_dir, "paired_difference.csv"), row.names = FALSE)

  pair <- do.call(rbind, lapply(split(wide, wide$severity), function(s) {
    n <- nrow(s); se <- sd(s$delta) / sqrt(n)
    tq <- if (n > 1L) qt(0.975, df = n - 1L) else NA_real_
    data.frame(severity = s$severity[1], pB_focal = s$pB_focal[1], n_rep = n,
               mean_delta = mean(s$delta),
               lo = mean(s$delta) - tq * se, hi = mean(s$delta) + tq * se,
               n_positive = sum(s$delta > 0),
               p_signed_rank = if (n > 1L && length(unique(s$delta)) > 1L)
                 suppressWarnings(wilcox.test(s$delta)$p.value) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  pair <- pair[order(pair$pB_focal), ]
  write.csv(pair, file.path(cfg$out_dir, "paired_difference_summary.csv"), row.names = FALSE)
  cat("\n=== Paired: Delta_r = FPR_M - FPR_L, matched dataset ===\n")
  print(pair, digits = 3, row.names = FALSE)

  invisible(list(reps = reps, summary = summ, primary = agg))
}

## 10. Command line (--smoke, --set key=value, --norun)
args <- commandArgs(trailingOnly = TRUE)
if (any(grepl("--smoke", args))) {
  cfg$px <- 120L; cfg$py <- 100L; cfg$pz <- 80L
  cfg$cells_per_sample <- 300L; cfg$cell_jitter <- 60L
  cfg$n_sub <- 22L; cfg$n_reps <- 2L
  cfg$out_dir <- file.path(RESULTS_DIR, "S12_smoke")
  message("[smoke] reduced config active")
}
for (a in grep("^--set$", args)) {
  kv <- strsplit(args[a + 1], "=")[[1]]; v <- suppressWarnings(as.numeric(kv[2]))
  cfg[[kv[1]]] <- if (is.na(v)) kv[2] else v
}
cfg$n_samples <- as.integer(cfg$n_samples); cfg$n_reps <- as.integer(cfg$n_reps)
cfg$n_cores   <- as.integer(cfg$n_cores);   cfg$n_sub  <- as.integer(cfg$n_sub)

if (!any(grepl("--norun", args))) {
  res <- main(cfg)
  saveRDS(list(cfg = cfg, summary = res$summary, primary = res$primary,
               sessionInfo = sessionInfo()),
          file.path(cfg$out_dir, "run_metadata.rds"))
  cat("\nWrote to", normalizePath(cfg$out_dir), "\n")
}
