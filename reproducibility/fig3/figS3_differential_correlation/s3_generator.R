# Supplemental Figure S3: single-modality RNA count generator (sourced by the stage scripts)
# log-rate: shared program + DC loading + anchor + partners + private programs + residual; counts Poisson(L * lambda), then Binomial(Z, q)

make_regular_bipartite_graph <- function(cfg, cohort_seed) {
  n_m <- cfg$n_m; m <- cfg$m
  mod <- rep(seq_len(cfg$r), each = n_m)
  MA <- which(mod == cfg$dc_modules[1]); MB <- which(mod == cfg$dc_modules[2])
  A <- matrix(0L, n_m, n_m)
  for (t in 0:(m - 1)) A[cbind(seq_len(n_m), ((seq_len(n_m) - 1 + t) %% n_m) + 1)] <- 1L
  if (!all(rowSums(A) == m) || !all(colSums(A) == m)) stop("I1: graph is not m-regular")
  perm <- with_seed(seed_for(cohort_seed, "partner_perm"), list(a = sample.int(n_m), b = sample.int(n_m)))
  idx <- which(A == 1L, arr.ind = TRUE)
  edges <- cbind(f = MA[perm$a[idx[, 1]]], g = MB[perm$b[idx[, 2]]])
  deg <- tabulate(c(edges), nbins = cfg$F)
  if (!all(deg[c(MA, MB)] == m)) stop("I1: degree check failed after permutation")
  edges
}

make_cohort_parameters <- function(cfg, cohort_seed) {
  mod <- rep(seq_len(cfg$r), each = cfg$n_m)
  is_dc <- mod %in% cfg$dc_modules
  target <- ifelse(mod == cfg$dc_modules[1], cfg$dc_modules[2],
                   ifelse(mod == cfg$dc_modules[2], cfg$dc_modules[1], NA_integer_))
  mu <- with_seed(seed_for(cohort_seed, "mu"), exp(log(cfg$mu_median) + cfg$mu_sdlog * rnorm(cfg$F)))
  list(mod = mod, is_dc = is_dc, dc = which(is_dc), target = target, mu = mu,
       edges = make_regular_bipartite_graph(cfg, cohort_seed))
}

make_sample_private_programs <- function(cfg, cohort_seed, sample) {
  genes <- with_seed(seed_for(cohort_seed, "nuisance", sample), sample.int(cfg$F))[seq_len(cfg$qn * cfg$n_m)]
  split(genes, rep(seq_len(cfg$qn), each = cfg$n_m))
}

sample_donor_amplitude <- function(cfg, cohort_seed, sample) {
  e <- with_seed(seed_for(cohort_seed, "donor", sample), rnorm(cfg$F))
  sqrt(cfg$s2) * (1 + cfg$tau * pmin(pmax(e, -3), 3))
}

build_loadings <- function(cfg, par, s_i, programs, condition, kappa, beta2, eta) {
  F <- cfg$F
  L <- matrix(0, F, cfg$r); L[cbind(seq_len(F), par$mod)] <- s_i
  P <- diag(sqrt(cfg$sigma0_sq), F)
  if (condition == 1L && beta2 > 0) {
    b <- sqrt(kappa * beta2); cc <- sqrt((1 - kappa) * beta2 / cfg$m)
    L[cbind(par$dc, par$target[par$dc])] <- L[cbind(par$dc, par$target[par$dc])] + b
    e <- par$edges                                   # same adjacency, both directions
    P[e] <- P[e] + cc
    P[e[, 2:1, drop = FALSE]] <- P[e[, 2:1, drop = FALSE]] + cc
  }
  H <- matrix(0, F, cfg$qn)
  for (k in seq_along(programs)) H[programs[[k]], k] <- sqrt(eta)
  rho2 <- cfg$V - rowSums(L^2) - rowSums(P^2) - rowSums(H^2)
  if (!all(rho2 > cfg$rho2_floor))
    stop(sprintf("I2: variance budget infeasible, min rho^2 = %.4f", min(rho2)))
  list(L = L, P = P, H = H, rho2 = rho2)
}

analytic_covariance <- function(ld) {
  tcrossprod(ld$L) + tcrossprod(ld$P) + tcrossprod(ld$H) + diag(ld$rho2)
}

population_covariance_diagnostics <- function(cfg, par, cohort_seed, kappa, beta2, eta, sample) {
  s_i <- sample_donor_amplitude(cfg, cohort_seed, sample)
  prog <- make_sample_private_programs(cfg, cohort_seed, sample)
  SA <- analytic_covariance(build_loadings(cfg, par, s_i, prog, 0L, kappa, beta2, eta))
  SB <- analytic_covariance(build_loadings(cfg, par, s_i, prog, 1L, kappa, beta2, eta))
  var_err <- max(abs(diag(SA) - cfg$V), abs(diag(SB) - cfg$V))
  outside <- !outer(par$is_dc, par$is_dc)
  closed_err <- max(abs(SB - SA)[outside])
  if (var_err >= 1e-12) stop(sprintf("I3: diag(Sigma) != V, error %.2e", var_err))
  if (closed_err >= 1e-12) stop(sprintf("I4: covariance change outside DC x DC, error %.2e", closed_err))
  c(max_variance_error = var_err, closed_set_cov_error = closed_err)
}

simulate_latent_rates <- function(cfg, par, cohort_seed, sample, condition, kappa, beta2, eta, n_cells) {
  s_i <- sample_donor_amplitude(cfg, cohort_seed, sample)
  prog <- make_sample_private_programs(cfg, cohort_seed, sample)
  ld <- build_loadings(cfg, par, s_i, prog, condition, kappa, beta2, eta)
  F <- cfg$F; n <- n_cells
  x <- with_seed(seed_for(cohort_seed, "cell_latent", sample), {
    z    <- matrix(rnorm(n * cfg$r), n, cfg$r)
    xi   <- matrix(rnorm(n * F), n, F)
    eps  <- matrix(rnorm(n * F), n, F)
    zeta <- matrix(rnorm(n * cfg$qn), n, cfg$qn)
    tcrossprod(z, ld$L) + tcrossprod(xi, ld$P) + sweep(eps, 2, sqrt(ld$rho2), `*`) + tcrossprod(zeta, ld$H)
  })
  exp(sweep(x, 2, log(par$mu) - cfg$V / 2, `+`))                    # n x F rates
}

simulate_counts <- function(cfg, cohort_seed, sample, lam) {
  n <- nrow(lam)
  lib <- with_seed(seed_for(cohort_seed, "library", sample),
                   rlnorm(n, -cfg$lib_sdlog^2 / 2, cfg$lib_sdlog))
  Y <- matrix(0L, n, ncol(lam))
  for (j in seq_len(ncol(lam))) {
    Y[, j] <- with_seed(seed_for(cohort_seed, "count", sample, j), {
      Z <- rpois(n, lib * lam[, j])
      rbinom(n, Z, cfg$q)
    })
  }
  list(Y = Y, lib = lib)
}

simulate_background_counts <- function(cfg, cohort_seed, sample, lib) {
  with_seed(seed_for(cohort_seed, "background", sample), {
    n <- length(lib)
    lam <- exp(log(cfg$bg_mu) - cfg$bg_sd^2 / 2 + cfg$bg_sd * matrix(rnorm(n * cfg$K_bg), n, cfg$K_bg))
    Z <- matrix(rpois(n * cfg$K_bg, lib * lam), n, cfg$K_bg)
    matrix(rbinom(length(Z), Z, cfg$q), n, cfg$K_bg)
  })
}

# One cohort; by default donor i is sample i and the first S/2 samples are condition A
simulate_cohort <- function(cfg, cohort_seed, kappa, beta2, eta, n_cells = cfg$n_cells,
                            donors = seq_len(cfg$S),
                            conditions = rep(0:1, c(cfg$S %/% 2, cfg$S - cfg$S %/% 2))) {
  simulate_cohort_with_par(cfg, cohort_seed, make_cohort_parameters(cfg, cohort_seed),
                           kappa, beta2, eta, n_cells, donors, conditions)
}

# CRN reference Delta: each reference donor simulated under A and B with identical draws, split-half estimates
reference_delta_crn <- function(cfg, cohort_seed, kappa, beta2) {
  ref_seed <- seed_for(cohort_seed, "crn_reference")
  nd <- cfg$ref_donors
  # the reference uses the cohort's own structure (mu, graph); only the draws come from ref_seed
  par <- make_cohort_parameters(cfg, cohort_seed)
  coh <- simulate_cohort_with_par(cfg, ref_seed, par, kappa, beta2, 0, cfg$ref_cells,
                                  donors = rep(seq_len(nd), 2), conditions = rep(0:1, each = nd))
  U <- mosaic_couplings(coh, cfg)                                # list of F x F
  diffs <- lapply(seq_len(nd), function(d) U[[nd + d]] - U[[d]])
  h <- nd %/% 2
  D1 <- Reduce(`+`, diffs[1:h]) / h
  D2 <- Reduce(`+`, diffs[(h + 1):nd]) / h
  D <- (D1 + D2) / 2
  norm2 <- rowSums(D1 * D2)                                      # unbiased ||Delta_f||^2
  list(D = D, D1 = D1, D2 = D2, norm_corr = sqrt(pmax(norm2, 0)), par = par,
       leakage = stats::median(sqrt(rowSums(D[!par$is_dc, !par$is_dc]^2))))
}

# unbiased alignment A_f = ||Q' Delta_f||^2 / ||Delta_f||^2 from split halves
reference_alignment <- function(ref, Q) {
  rowSums((ref$D1 %*% Q) * (ref$D2 %*% Q)) / rowSums(ref$D1 * ref$D2)
}

# same as simulate_cohort but with externally supplied cohort parameters
simulate_cohort_with_par <- function(cfg, seed, par, kappa, beta2, eta, n_cells, donors, conditions) {
  Y <- vector("list", length(donors)); B <- Y
  for (k in seq_along(donors)) {
    lam <- simulate_latent_rates(cfg, par, seed, donors[k], conditions[k], kappa, beta2, eta, n_cells)
    ct <- simulate_counts(cfg, seed, donors[k], lam)
    Y[[k]] <- ct$Y
    B[[k]] <- simulate_background_counts(cfg, seed, donors[k], ct$lib)
  }
  list(Y = Y, Ybg = B, conditions = conditions, par = par)
}

# pooled alignment over a gene set (ratio of summed split-half cross products)
reference_alignment_pooled <- function(ref, Q, genes) {
  num <- rowSums((ref$D1 %*% Q) * (ref$D2 %*% Q))[genes]
  den <- rowSums(ref$D1 * ref$D2)[genes]
  sum(num) / sum(den)
}

# mean over cohort seeds of the median noise-free ||Delta_f|| of DC genes
mean_norm_over_seeds <- function(cfg, seeds, kappa, b2) {
  mean(vapply(seeds, function(s) {
    r <- reference_delta_crn(cfg, s, kappa, b2); median(r$norm_corr[r$par$is_dc])
  }, 0))
}
