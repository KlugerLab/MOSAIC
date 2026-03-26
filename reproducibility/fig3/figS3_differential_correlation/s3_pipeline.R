# Supplemental Figure S3: one cohort through the non-memento arms (sourced by the stage scripts)
run_cohort_nonmemento <- function(cfg, seed, kappa, beta2, eta, raw_only = FALSE, d_fixed = NULL,
                                  skip_dgca = FALSE) {
  coh <- simulate_cohort(cfg, seed, kappa, beta2, eta)
  obj <- preprocess_seurat(coh, cfg)
  U <- lapply(scaled_by_sample(obj, cfg), function(m) cosine_distance(m, m))
  labels <- coh$par$is_dc
  F_raw <- pseudo_f_all(U, coh$conditions)
  if (raw_only) return(list(scores = list(raw = F_raw), labels = labels))
  mb <- run_mosaic_one(U, cfg$n_eig)
  idn <- projection_identity(U, mb$W, coh$conditions)             # stops on violation
  F_unw <- pseudo_f_all(lapply(U, function(M) M %*% mb$Q), coh$conditions)
  F_fix <- NULL
  if (!is.null(d_fixed)) {
    Wf <- mb$agg_vectors[, 1:d_fixed, drop = FALSE] %*% diag(mb$agg_values[1:d_fixed], nrow = d_fixed)
    F_fix <- pseudo_f_all(lapply(U, function(M) M %*% Wf), coh$conditions)
  }
  dg <- if (skip_dgca) rep(NA_real_, cfg$F) else run_dgca_scores(obj, cfg)
  cond_of_sample <- coh$conditions
  list(scores = c(list(raw = idn$F_raw, mosaic = idn$F_M, dgca = unname(dg),
                       mosaic_unweighted = F_unw),
                  if (!is.null(F_fix)) list(mosaic_fixed_d = F_fix)),
       labels = labels, G = idn$G, Q = mb$Q, W = mb$W, R_N = idn$R_N,
       ssw_raw = idn$ssw_raw, ssw_W = idn$ssw_W, n_A = sum(coh$conditions == 0), n_B = sum(coh$conditions == 1),
       diag = data.frame(d = mb$d,
                         rank_median_A = median(mb$stage1_ranks[cond_of_sample == 0]),
                         rank_median_B = median(mb$stage1_ranks[cond_of_sample == 1]),
                         max_identity_error = idn$max_identity_error,
                         n_degenerate = idn$n_degenerate))
}
