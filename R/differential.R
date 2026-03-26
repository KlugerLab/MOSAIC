#' Run Differential Connectivity Test
#'
#' Tests each feature for a condition-associated shift in its MOSAIC
#' embedding. For each feature, a PERMANOVA pseudo-F statistic is computed on
#' the distances between its sample embeddings, and significance is assessed by
#' permuting condition labels across samples with the embedding held fixed.
#' The same label assignments are used for all features; each feature is
#' compared with its own null distribution.
#'
#' Without \code{blocks}, all assignments preserving the group sizes are
#' admissible; with \code{blocks}, labels are permuted within blocks. If the
#' number of admissible assignments is at most \code{max_exact}, all of them
#' are enumerated and exact P-values are returned. Otherwise \code{n_perm}
#' random assignments are used and \eqn{p = (1 + \sum_b I\{F^{(b)} \ge F\}) / (B + 1)}.
#' P-values are adjusted with the Benjamini-Hochberg procedure.
#'
#' @param mosaic_embed_list Named list of per-sample embedding matrices
#'   (features x latent dimensions), as returned by
#'   \code{\link{run_MOSAIC}()$mosaic_embed_list}. All matrices must have the
#'   same features in the same order.
#' @param groups Character or factor vector of condition labels, one per
#'   sample, in the order of \code{mosaic_embed_list}.
#' @param blocks Optional vector of block labels (e.g. donor IDs), one per
#'   sample. If given, labels are permuted only within blocks. Default
#'   \code{NULL} (all samples exchangeable).
#' @param n_perm Integer. Number of random assignments used when the
#'   admissible set is too large to enumerate. Default 999.
#' @param max_exact Numeric. Enumerate all admissible assignments if their
#'   number is at most this value. Set to 0 to always use random assignments.
#'   Default 2e5.
#' @param dist_method Distance between the sample embeddings of a feature:
#'   \code{"euclidean"} (default) or \code{"cosine"}
#'   (\eqn{1 -} cosine similarity).
#' @param seed Optional integer seed for the random assignments.
#' @param feature_chunk,perm_chunk Integers controlling memory use: number of
#'   features and of assignments processed at once.
#' @param return_null Logical. If \code{TRUE}, also return the full
#'   features x assignments matrix of permutation pseudo-F statistics and the
#'   assignment matrix. Default \code{FALSE}.
#'
#' @return A list with elements
#' \describe{
#'   \item{\code{results}}{Data frame with one row per feature: \code{feature},
#'     \code{F_stat} (observed pseudo-F), \code{R2} (\eqn{SS_B / SS_T}),
#'     \code{pvalue} and \code{padj} (Benjamini-Hochberg).}
#'   \item{\code{exact}}{Logical, whether the P-values are exact.}
#'   \item{\code{n_perm}}{Number of assignments evaluated (\eqn{|\Pi|} if
#'     exact).}
#'   \item{\code{min_pvalue}}{Smallest attainable P-value.}
#'   \item{\code{F_null}, \code{assignments}}{Only if \code{return_null = TRUE}:
#'     the permutation statistics (features x assignments) and the
#'     assignments (samples x assignments, as group labels).}
#' }
#'
#' @seealso \code{\link{run_MOSAIC}} to generate the embedding,
#'   \code{\link{compute_feature_similarity}} for per-feature sample similarity
#'   matrices used in visualization and subgroup detection.
#'
#' @examples
#' \dontrun{
#' result <- run_MOSAIC(
#'   list(RNA = seurat_rna, ADT = seurat_adt),
#'   assays = c("RNA", "ADT"),
#'   sample_meta = "sample_id",
#'   condition_meta = "time"
#' )
#'
#' # Independent samples
#' dc <- run_DC_test(result$mosaic_embed_list,
#'                   groups = result$annotation$Condition)
#'
#' # Paired design: permute time labels within donor
#' dc <- run_DC_test(result$mosaic_embed_list,
#'                   groups = result$annotation$Condition,
#'                   blocks = donor_ids)
#' head(dc$results[order(dc$results$pvalue, -dc$results$F_stat), ])
#' }
#'
#' @export
run_DC_test <- function(
    mosaic_embed_list,
    groups,
    blocks = NULL,
    n_perm = 999,
    max_exact = 2e5,
    dist_method = c("euclidean", "cosine"),
    seed = NULL,
    feature_chunk = 2000,
    perm_chunk = 2000,
    return_null = FALSE
) {
  dist_method <- match.arg(dist_method)
  .check_embed_list(mosaic_embed_list)

  S <- length(mosaic_embed_list)
  if (length(groups) != S) {
    stop("`groups` must have one label per sample (", S, ").")
  }
  groups <- as.factor(groups)
  group_codes <- as.integer(groups)
  G <- nlevels(groups)
  if (G < 2) stop("`groups` must contain at least two conditions.")
  if (S - G < 1) stop("Too few samples for the number of conditions.")
  n_g <- tabulate(group_codes, nbins = G)

  if (is.null(blocks)) {
    blocks <- rep(1L, S)
  } else if (length(blocks) != S) {
    stop("`blocks` must have one entry per sample (", S, ").")
  }
  block_idx <- split(seq_len(S), as.factor(blocks))

  # ---- admissible assignments --------------------------------------------
  n_admissible <- round(prod(vapply(block_idx, function(idx) {
    exp(lfactorial(length(idx)) - sum(lfactorial(tabulate(group_codes[idx], G))))
  }, numeric(1))))
  if (n_admissible < 2) {
    stop("The design admits only the observed assignment; no permutation test is possible.")
  }

  exact <- n_admissible <= max_exact
  if (exact) {
    assign_mat <- .enumerate_assignments(group_codes, block_idx, G)
  } else {
    if (!is.null(seed)) set.seed(seed)
    assign_mat <- .sample_assignments(group_codes, block_idx, n_perm)
  }
  P <- ncol(assign_mat)

  # ---- pseudo-F for all features and assignments -------------------------
  pairs <- which(upper.tri(matrix(0, S, S)), arr.ind = TRUE)
  ii <- pairs[, 1]
  jj <- pairs[, 2]

  n_feat <- nrow(mosaic_embed_list[[1]])
  feature_names <- rownames(mosaic_embed_list[[1]])
  if (is.null(feature_names)) feature_names <- paste0("feature_", seq_len(n_feat))

  obs_W <- .pair_weights(matrix(group_codes, ncol = 1), ii, jj, G)

  F_obs <- R2 <- rep(NA_real_, n_feat)
  n_ge <- numeric(n_feat)
  F_null <- if (return_null) matrix(NA_real_, n_feat, P) else NULL

  feat_chunks <- split(seq_len(n_feat), ceiling(seq_len(n_feat) / feature_chunk))
  perm_chunks <- split(seq_len(P), ceiling(seq_len(P) / perm_chunk))

  for (fc in feat_chunks) {
    K <- .pairwise_sq_dist(mosaic_embed_list, fc, ii, jj, dist_method)
    ss_t <- rowSums(K) / S

    stat_obs <- .pseudo_F(K, ss_t, obs_W, n_g, S, G)
    F_obs[fc] <- stat_obs$F[, 1]
    R2[fc] <- stat_obs$R2[, 1]
    tol <- 1e-9 * pmax(1, abs(F_obs[fc]))

    for (pc in perm_chunks) {
      W <- .pair_weights(assign_mat[, pc, drop = FALSE], ii, jj, G)
      F_perm <- .pseudo_F(K, ss_t, W, n_g, S, G)$F
      n_ge[fc] <- n_ge[fc] + rowSums(F_perm >= (F_obs[fc] - tol))
      if (return_null) F_null[fc, pc] <- F_perm
    }
  }

  pvalue <- if (exact) n_ge / P else (1 + n_ge) / (P + 1)
  pvalue[is.na(F_obs)] <- NA_real_
  min_p <- if (exact) .n_equivalent(assign_mat, group_codes, G) / P else 1 / (P + 1)

  results <- data.frame(
    feature = feature_names,
    F_stat = F_obs,
    R2 = R2,
    pvalue = pvalue,
    padj = stats::p.adjust(pvalue, method = "BH"),
    stringsAsFactors = FALSE
  )

  out <- list(results = results, exact = exact, n_perm = P, min_pvalue = min_p)
  if (return_null) {
    rownames(F_null) <- feature_names
    out$F_null <- F_null
    out$assignments <- matrix(levels(groups)[assign_mat], nrow = S,
                              dimnames = list(names(mosaic_embed_list), NULL))
  }
  out
}


# ---------------------------------------------------------------------------
# Internal: input checks
# ---------------------------------------------------------------------------
.check_embed_list <- function(mosaic_embed_list) {
  if (!is.list(mosaic_embed_list) || length(mosaic_embed_list) < 3) {
    stop("`mosaic_embed_list` must be a list of at least three sample matrices.")
  }
  dims <- vapply(mosaic_embed_list, function(m) dim(as.matrix(m)), integer(2))
  if (any(dims[1, ] != dims[1, 1]) || any(dims[2, ] != dims[2, 1])) {
    stop("All matrices in `mosaic_embed_list` must have the same dimensions.")
  }
  invisible(TRUE)
}


# ---------------------------------------------------------------------------
# Internal: squared distances between sample embeddings for a chunk of
# features and all sample pairs i < j. Returns |features| x |pairs|.
# ---------------------------------------------------------------------------
.pairwise_sq_dist <- function(mosaic_embed_list, feat_idx, ii, jj, dist_method) {
  X <- lapply(mosaic_embed_list, function(m) as.matrix(m)[feat_idx, , drop = FALSE])
  if (dist_method == "cosine") {
    nrm <- lapply(X, function(x) sqrt(rowSums(x^2)))
  }
  K <- matrix(0, length(feat_idx), length(ii))
  for (p in seq_along(ii)) {
    a <- X[[ii[p]]]
    b <- X[[jj[p]]]
    if (dist_method == "euclidean") {
      K[, p] <- rowSums((a - b)^2)
    } else {
      den <- nrm[[ii[p]]] * nrm[[jj[p]]]
      cs <- rowSums(a * b) / den
      cs[den == 0] <- 0
      K[, p] <- (1 - cs)^2
    }
  }
  K
}


# ---------------------------------------------------------------------------
# Internal: for each group g, the pair indicators 1{a_i = g} 1{a_j = g}
# (|pairs| x |assignments|).
# ---------------------------------------------------------------------------
.pair_weights <- function(assign_mat, ii, jj, G) {
  lapply(seq_len(G), function(g) {
    ind <- assign_mat == g
    storage.mode(ind) <- "double"
    ind[ii, , drop = FALSE] * ind[jj, , drop = FALSE]
  })
}


# ---------------------------------------------------------------------------
# Internal: PERMANOVA pseudo-F and R2 from squared distances over pairs.
# SS_T = sum_{i<j} d_ij^2 / S;  SS_W = sum_g sum_{i<j in g} d_ij^2 / n_g.
# ---------------------------------------------------------------------------
.pseudo_F <- function(K, ss_t, W, n_g, S, G) {
  ss_w <- 0
  for (g in seq_len(G)) {
    ss_w <- ss_w + (K %*% W[[g]]) / n_g[g]
  }
  ss_b <- ss_t - ss_w
  list(F = (ss_b / (G - 1)) / (ss_w / (S - G)), R2 = ss_b / ss_t)
}


# ---------------------------------------------------------------------------
# Internal: all distinct arrangements of a multiset of group codes, one
# arrangement per column.
# ---------------------------------------------------------------------------
.multiset_arrangements <- function(counts, codes = seq_along(counts)) {
  n <- sum(counts)
  if (length(counts) == 1 || n == counts[1]) {
    return(matrix(codes[1], nrow = n, ncol = 1))
  }
  sub <- .multiset_arrangements(counts[-1], codes[-1])
  if (counts[1] == 0) {
    return(sub)
  }
  cmb <- utils::combn(n, counts[1])
  out <- matrix(0L, n, ncol(cmb) * ncol(sub))
  for (k in seq_len(ncol(cmb))) {
    cols <- (k - 1) * ncol(sub) + seq_len(ncol(sub))
    out[cmb[, k], cols] <- codes[1]
    out[setdiff(seq_len(n), cmb[, k]), cols] <- sub
  }
  out
}


# ---------------------------------------------------------------------------
# Internal: enumerate all admissible assignments (labels permuted within
# blocks).
# ---------------------------------------------------------------------------
.enumerate_assignments <- function(group_codes, block_idx, G) {
  per_block <- lapply(block_idx, function(idx) {
    .multiset_arrangements(tabulate(group_codes[idx], G), seq_len(G))
  })
  grid <- as.matrix(expand.grid(lapply(per_block, function(a) seq_len(ncol(a)))))
  out <- matrix(0L, length(group_codes), nrow(grid))
  for (b in seq_along(block_idx)) {
    out[block_idx[[b]], ] <- per_block[[b]][, grid[, b], drop = FALSE]
  }
  out
}


# ---------------------------------------------------------------------------
# Internal: random admissible assignments (labels permuted within blocks).
# ---------------------------------------------------------------------------
.sample_assignments <- function(group_codes, block_idx, n_perm) {
  out <- matrix(group_codes, length(group_codes), n_perm)
  for (idx in block_idx) {
    if (length(idx) > 1) {
      out[idx, ] <- vapply(seq_len(n_perm),
                           function(b) group_codes[idx][sample.int(length(idx))],
                           integer(length(idx)))
    }
  }
  out
}


# ---------------------------------------------------------------------------
# Internal: number of enumerated assignments that necessarily reproduce the
# observed pseudo-F: the observed labeling and, for two groups, its
# complement when exchanging the labels is itself admissible.
# ---------------------------------------------------------------------------
.n_equivalent <- function(assign_mat, group_codes, G) {
  if (G != 2) return(1)
  swapped <- 3L - group_codes
  if (any(colSums(assign_mat != swapped) == 0)) 2 else 1
}
