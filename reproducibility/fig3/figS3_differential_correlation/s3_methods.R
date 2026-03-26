# Supplemental Figure S3: methods compared (MOSAIC, no-projection ablation, DGCA, memento; sourced by the stage scripts)

suppressPackageStartupMessages({
  library(Matrix)
  library(RSpectra)
})
# Seurat is loaded only where needed
use_seurat <- function() suppressPackageStartupMessages(library(Seurat))

# MOSAIC helpers (from 0-functions.R)
cosine_distance <- function(A, B) {
    A <- t(A)
    B <- t(B)
    A_norm <- sqrt(rowSums(A^2))
    B_norm <- sqrt(rowSums(B^2))

    A_norm[A_norm == 0] <- 1
    B_norm[B_norm == 0] <- 1

    A_normalized <- A / A_norm
    B_normalized <- B / B_norm

    similarity <- tcrossprod(A_normalized, B_normalized)
    return(similarity)
}

find_elbow_kneedle <- function(eigenvalues) {
  for (i in 1:length(eigenvalues)){
    if (eigenvalues[i] <= 0){
      eigenvalues[i] <- eigenvalues[i-1]
    }
  }
  x <- 1:length(eigenvalues)
  y <- log(eigenvalues)
  if (length(x) < 3) {
    return(length(x))
  }
  x_scaled <- (x - min(x)) / (max(x) - min(x))
  y_scaled <- (y - min(y)) / (max(y) - min(y))
  distances <- perpendicular_distance(
    xp = x_scaled,
    yp = y_scaled,
    x1 = x_scaled[1],
    y1 = y_scaled[1],
    x2 = x_scaled[length(x_scaled)],
    y2 = y_scaled[length(y_scaled)]
  )
  elbow_index <- which.max(distances)
  output <- max(elbow_index,2)
  return(output)
}

perpendicular_distance <- function(xp, yp, x1, y1, x2, y2) {
  A <- y1 - y2
  B <- x2 - x1
  C <- -A * x1 - B * y1
  distance <- abs(A * xp + B * yp + C) / sqrt(A^2 + B^2)
  return(distance)
}

remove_zero_rows <- function(mat) {
  zero_rows <- apply(mat, 1, function(row) all(row == 0))
  return(list(mat[!zero_rows, , drop = FALSE],zero_rows))
}

structured_names <- function(cfg) sprintf("g%04d", seq_len(cfg$F))
background_names <- function(cfg) sprintf("bg%05d", seq_len(cfg$K_bg))
sample_names <- function(n) sprintf("s%02d", seq_len(n))

# Seurat preprocessing: NormalizeData on all genes, ScaleData on the structured genes
preprocess_seurat <- function(coh, cfg) {
  use_seurat()
  n_s <- length(coh$Y)
  sn <- sample_names(n_s)
  counts <- do.call(rbind, Map(cbind, coh$Y, coh$Ybg))                   # cells x genes
  colnames(counts) <- c(structured_names(cfg), background_names(cfg))
  cell_ids <- unlist(lapply(seq_len(n_s), function(i) sprintf("%s_c%04d", sn[i], seq_len(nrow(coh$Y[[i]])))))
  rownames(counts) <- cell_ids
  meta <- data.frame(sample_id = rep(sn, vapply(coh$Y, nrow, 1L)),
                     condition = rep(coh$conditions, vapply(coh$Y, nrow, 1L)),
                     row.names = cell_ids)
  obj <- CreateSeuratObject(counts = as(t(counts), "CsparseMatrix"), meta.data = meta,
                            min.cells = 0, min.features = 0)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- ScaleData(obj, features = structured_names(cfg), verbose = FALSE)
  obj
}

# per-sample cells x genes scaled matrices, in sample order
scaled_by_sample <- function(obj, cfg) {
  use_seurat()
  sc <- LayerData(obj, assay = "RNA", layer = "scale.data")[structured_names(cfg), , drop = FALSE]
  sid <- obj$sample_id
  lapply(sort(unique(sid)), function(s) t(as.matrix(sc[, sid == s, drop = FALSE])))
}

# per-sample F x F cosine block matrices (eigen_block_cosine_one's block_matrix)
mosaic_couplings <- function(coh, cfg) {
  lapply(scaled_by_sample(preprocess_seurat(coh, cfg), cfg), function(m) cosine_distance(m, m))
}

# MOSAIC Stage 1 and Stage 2, one modality (run_our_method_three logic)
run_mosaic_one <- function(U_list, n_eig = 50L) {
  eig <- lapply(U_list, function(U) eigs_sym(U, n_eig, which = "LA"))
  r_list <- vapply(eig, function(e) find_elbow_kneedle(e$values), 1)
  proj <- Map(function(e, r) {
    if (r == 1) r <- 2
    V <- e$vectors[, 1:r, drop = FALSE]
    V %*% diag(e$values[1:r], nrow = r) %*% t(V)
  }, eig, r_list)
  sum_proj <- Reduce(`+`, proj)
  V_result <- eigs_sym(sum_proj, n_eig, which = "LA")
  r2 <- find_elbow_kneedle(V_result$values[1:n_eig])
  W <- V_result$vectors[, 1:r2, drop = FALSE] %*% diag(V_result$values[1:r2], nrow = r2)
  list(W = W, Q = V_result$vectors[, 1:r2, drop = FALSE], d = r2, stage1_ranks = r_list,
       agg_values = V_result$values, agg_vectors = V_result$vectors,
       E = lapply(U_list, function(U) U %*% W))
}

# Two-group Euclidean pseudo-F for every feature (rows); all-zero sample rows excluded
pseudo_f_all <- function(M_list, groups) {
  S <- length(M_list); nF <- nrow(M_list[[1]])
  A <- simplify2array(M_list)                                     # F x p x S
  zero <- apply(A, c(1, 3), function(v) all(v == 0))              # F x S
  g <- as.integer(factor(groups))
  if (!any(zero)) {
    out <- numeric(nF)
    mu <- apply(A, c(1, 2), mean)
    ss_t <- rowSums(apply(A, 3, function(x) rowSums((x - mu)^2)))
    ss_w <- 0
    for (k in unique(g)) {
      Ak <- A[, , g == k, drop = FALSE]
      mk <- apply(Ak, c(1, 2), mean)
      ss_w <- ss_w + rowSums(apply(Ak, 3, function(x) rowSums((x - mk)^2)))
    }
    if (any(ss_w <= 0)) stop("zero within-group sum of squares")
    G <- length(unique(g))
    return(((ss_t - ss_w) / (G - 1)) / (ss_w / (S - G)))
  }
  vapply(seq_len(nF), function(f) {                               # rare path: per feature
    X <- t(A[f, , ])
    keep <- !zero[f, ]
    X <- X[keep, , drop = FALSE]; gg <- g[keep]
    mu <- colMeans(X); ss_t <- sum(sweep(X, 2, mu)^2)
    ss_w <- sum(vapply(unique(gg), function(k) sum(sweep(X[gg == k, , drop = FALSE], 2,
                   colMeans(X[gg == k, , drop = FALSE]))^2), 0))
    G <- length(unique(gg))
    ((ss_t - ss_w) / (G - 1)) / (ss_w / (sum(keep) - G))
  }, 0)
}

# Exact identity F_M = F_raw * R_S / R_N (asserted per feature).
projection_identity <- function(U_list, W, groups, tol = 1e-10, eps = 1e-14) {
  F_raw <- pseudo_f_all(U_list, groups)
  E_list <- lapply(U_list, function(U) U %*% W)
  F_M <- pseudo_f_all(E_list, groups)
  g <- as.integer(factor(groups)); a <- g == 1; b <- g == 2
  mean_of <- function(L) Reduce(`+`, L) / length(L)
  dhat <- mean_of(U_list[b]) - mean_of(U_list[a])
  res <- c(lapply(U_list[a], function(U) U - mean_of(U_list[a])),
           lapply(U_list[b], function(U) U - mean_of(U_list[b])))
  dW <- dhat %*% W
  num_s <- rowSums(dW^2); den_s <- rowSums(dhat^2)
  num_n <- Reduce(`+`, lapply(res, function(R) rowSums((R %*% W)^2)))
  den_n <- Reduce(`+`, lapply(res, function(R) rowSums(R^2)))
  degenerate <- (den_s < eps * den_n) | (num_n < eps)
  R_S <- ifelse(degenerate, NA, num_s / den_s)
  R_N <- ifelse(degenerate, NA, num_n / den_n)
  G <- R_S / R_N
  rel <- abs(F_M - F_raw * G) / pmax(abs(F_M), eps)
  max_err <- if (any(!degenerate)) max(rel[!degenerate]) else NA_real_
  if (!is.na(max_err) && !(max_err < tol))
    stop(sprintf("projection identity violated: max relative error %.3e", max_err))
  list(F_raw = F_raw, F_M = F_M, R_S = R_S, R_N = R_N, G = G,
       ssw_raw = den_n, ssw_W = num_n,
       max_identity_error = max_err, n_degenerate = sum(degenerate))
}

# Pairwise methods: feature score = RMS over partners of the pair statistic
rms_feature_scores <- function(gene1, gene2, stat, genes) {
  ok <- is.finite(stat)
  ss <- tapply(c(stat[ok], stat[ok])^2, c(gene1[ok], gene2[ok]), sum)
  nn <- tapply(rep(1, 2 * sum(ok)), c(gene1[ok], gene2[ok]), sum)
  out <- sqrt(ss[genes] / nn[genes]); names(out) <- genes
  out
}

# DGCA on log-normalised data, cells pooled within condition
run_dgca_scores <- function(obj, cfg) {
  use_seurat()
  suppressPackageStartupMessages(library(DGCA))
  genes <- structured_names(cfg)
  X <- as.matrix(LayerData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE])  # genes x cells
  cond <- factor(ifelse(obj$condition == 0, "A", "B"), levels = c("A", "B"))
  design <- model.matrix(~ 0 + cond); colnames(design) <- levels(cond)
  res <- ddcorAll(inputMat = X, design = design, compare = c("A", "B"), corrType = "pearson",
                  nPairs = "all", adjust = "none", nPerms = 0, classify = FALSE,
                  heatmapPlot = FALSE, getDCorAvg = FALSE, verbose = FALSE)
  if (nrow(res) != choose(length(genes), 2)) stop("DGCA did not return all pairs")
  rms_feature_scores(as.character(res$Gene1), as.character(res$Gene2), res$zScoreDiff, genes)
}

# memento in a separate Python process (counts passed as an int32 binary file)
run_memento_scores <- function(coh, cfg, random_state, python, num_cpus = cfg$n_workers,
                               num_boot = cfg$memento_boot, return_pairs = FALSE,
                               test_genes = structured_names(cfg), pairs = NULL) {
  runner <- file.path(S3_DIR, "s3_memento_runner.py")
  n_s <- length(coh$Y)
  counts <- do.call(rbind, Map(cbind, coh$Y, coh$Ybg))              # cells x genes, integer
  if (!is.integer(counts)) stop("counts must be integer")
  wd <- tempfile("memento_"); dir.create(wd)
  on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  con <- file(file.path(wd, "counts.bin"), "wb")
  writeBin(as.vector(t(counts)), con, size = 4L, endian = "little")   # row-major
  close(con)
  meta <- list(
    n_cells = nrow(counts), n_genes = ncol(counts),
    gene_names = c(structured_names(cfg), background_names(cfg)),
    sample_ids = rep(sample_names(n_s), vapply(coh$Y, nrow, 1L)),
    conditions = as.integer(rep(coh$conditions, vapply(coh$Y, nrow, 1L))),
    q = cfg$q, test_genes = test_genes, num_boot = as.integer(num_boot),
    num_cpus = as.integer(num_cpus), random_state = as.integer(random_state),
    min_cell_count = as.integer(cfg$memento_min_cell_count),
    min_perc_group = cfg$memento_min_perc_group, n_groups_expected = as.integer(n_s))
  if (!is.null(pairs)) meta$pairs <- unname(lapply(seq_len(nrow(pairs)), function(i) c(pairs[i, 1], pairs[i, 2])))
  jsonlite::write_json(meta, file.path(wd, "meta.json"), auto_unbox = TRUE, digits = NA)
  out_csv <- file.path(wd, "result.csv"); log <- file.path(wd, "log.txt")
  status <- system2(python, c(shQuote(runner), "--input", shQuote(wd), "--output", shQuote(out_csv)),
                    stdout = log, stderr = log,
                    env = c("OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1",
                            "MKL_NUM_THREADS=1", "NUMEXPR_NUM_THREADS=1"))
  if (!identical(status, 0L) || !file.exists(out_csv))
    stop("memento subprocess failed (status ", status, "):\n",
         paste(utils::tail(readLines(log), 30), collapse = "\n"))
  res <- utils::read.csv(out_csv, stringsAsFactors = FALSE)
  n_expected <- if (is.null(pairs)) choose(length(test_genes), 2) else nrow(pairs)
  if (nrow(res) != n_expected) stop("memento returned an unexpected number of pairs")
  if (!is.null(pairs)) return(res)                        # chunk: the caller merges and scores
  z <- res$corr_coef / res$corr_se
  scores <- rms_feature_scores(res$gene_1, res$gene_2, z, test_genes)
  if (return_pairs) list(scores = scores, pairs = res) else scores
}

# memento over many cohorts, one cohort per worker
run_memento_parallel <- function(tasks, cfg, workers, python, code_dir = S3_DIR) {
  cl <- parallel::makePSOCKcluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c("cfg", "code_dir", "python"), envir = environment())
  parallel::clusterEvalQ(cl, {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
               MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
    S3_DIR <- code_dir
    source(file.path(code_dir, "s3_config.R"))
    source(file.path(code_dir, "s3_generator.R"))
    source(file.path(code_dir, "s3_methods.R"))
    TRUE
  })
  parallel::parLapplyLB(cl, tasks, function(tk) {
    cfg_t <- cfg                                   # optional per-task F / cells
    if (!is.null(tk$F)) { cfg_t$F <- tk$F; cfg_t$n_m <- tk$F %/% cfg_t$r }
    if (!is.null(tk$cells)) cfg_t$n_cells <- tk$cells
    coh <- simulate_cohort(cfg_t, tk$seed, tk$kappa, tk$beta2, tk$eta)
    genes <- structured_names(cfg_t)
    if (!is.null(tk$n_test_genes)) genes <- genes[seq_len(tk$n_test_genes)]
    t0 <- proc.time()[["elapsed"]]
    out <- run_memento_scores(coh, cfg_t, tk$random_state, python = python, num_cpus = 1L,
                              num_boot = cfg_t$memento_boot, test_genes = genes, return_pairs = TRUE)
    list(task = tk, scores = out$scores, pairs = out$pairs, sec = proc.time()[["elapsed"]] - t0)
  })
}

# memento with pair chunking: each (cohort, chunk) is one single-CPU task; returns feature scores by cohort id
all_pairs <- function(genes) {
  cmb <- utils::combn(genes, 2)
  data.frame(gene_1 = cmb[1, ], gene_2 = cmb[2, ], stringsAsFactors = FALSE)
}

run_memento_chunked <- function(cohorts, cfg, workers, python, n_chunks = cfg$memento_chunks,
                                code_dir = S3_DIR, stat = c("z", "coef")) {
  # stat = "z": RMS of corr_coef / corr_se; "coef": RMS of corr_coef
  stat <- match.arg(stat)
  genes <- structured_names(cfg)
  pr <- all_pairs(genes)
  chunk_of <- if (n_chunks == 1) rep(1L, nrow(pr)) else cut(seq_len(nrow(pr)), n_chunks, labels = FALSE)
  tasks <- list()
  for (co in cohorts) for (ch in seq_len(n_chunks))
    tasks[[length(tasks) + 1]] <- c(co, list(chunk = ch, n_chunks = n_chunks))
  cl <- parallel::makePSOCKcluster(min(workers, length(tasks)))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c("cfg", "code_dir", "python", "pr", "chunk_of"), envir = environment())
  parallel::clusterEvalQ(cl, {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
               MKL_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
    S3_DIR <- code_dir
    source(file.path(code_dir, "s3_config.R"))
    source(file.path(code_dir, "s3_generator.R"))
    source(file.path(code_dir, "s3_methods.R"))
    TRUE
  })
  out <- parallel::parLapplyLB(cl, tasks, function(tk) {
    coh <- simulate_cohort(cfg, tk$seed, tk$kappa, tk$beta2, tk$eta)
    rs <- abs(seed_for(tk$seed, "memento_boot", tk$chunk)) %% .Machine$integer.max
    res <- run_memento_scores(coh, cfg, rs, python = python, num_cpus = 1L,
                              num_boot = cfg$memento_boot, pairs = pr[chunk_of == tk$chunk, , drop = FALSE])
    list(id = tk$id, chunk = tk$chunk, res = res)
  })
  ids <- vapply(out, function(x) x$id, "")
  setNames(lapply(unique(ids), function(i) {
    parts <- out[ids == i]
    parts <- parts[order(vapply(parts, function(x) x$chunk, 1))]
    if (length(parts) != n_chunks) stop("missing memento chunks for cohort ", i)
    res <- do.call(rbind, lapply(parts, function(x) x$res))
    if (!identical(res$gene_1, pr$gene_1) || !identical(res$gene_2, pr$gene_2))
      stop("memento chunks do not reassemble the full pair list for cohort ", i)
    rms_feature_scores(res$gene_1, res$gene_2,
                       if (stat == "z") res$corr_coef / res$corr_se else res$corr_coef, genes)
  }), unique(ids))
}
