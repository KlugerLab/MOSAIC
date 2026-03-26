# Figure S11: controlled test of the local-to-global mechanism (nuisance strength mu from 0.03 to 2)
# Output:   results/S11/stress_v3b.csv, results/S11/stress_v3b_fixedref.csv, results/S11/supp_table_mechanism.csv, figures/figS11.png

source("config.R")
source("R/functions.R")   # find_elbow_kneedle()

suppressPackageStartupMessages({
  library(RSpectra); library(dplyr); library(tidyr); library(ggplot2)
})

out_dir <- file.path(RESULTS_DIR, "S11")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

Q_NUIS <- 3      # nuisance directions
RHO    <- 0.9    # nuisance coherence

gen_v2 <- function(theta, mu, sigma_obs, S = 20, n = 500, r = 10,
                   FX = 550, FY = 450, m_shared = 6, sv = c(5, 1),
                   q = Q_NUIS, rho = RHO, seed = 1) {
  set.seed(seed)
  F <- FX + FY
  need <- r + S * (r - m_shared) + q
  stopifnot(need <= F)
  Q   <- qr.Q(qr(matrix(rnorm(F * need), F, need)))
  V0  <- Q[, 1:r, drop = FALSE]
  W   <- Q[, (need - q + 1):need, drop = FALSE]     # shared nuisance basis
  Sig <- diag(seq(sv[1], sv[2], length.out = r))
  wts <- seq(1, 0.4, length.out = q)                # decaying, no second elbow
  
  Z <- lapply(seq_len(S), function(i) matrix(rnorm(n * r), n, r))
  E <- lapply(seq_len(S), function(i) matrix(rnorm(n * F), n, F))
  A <- lapply(seq_len(S), function(i) matrix(rnorm(n * q), n, q))
  J <- lapply(seq_len(S), function(i) qr.Q(qr(matrix(rnorm(F * q), F, q))))
  
  # Xbio: signal only (target); X0: signal + obs noise; Xn: signal + obs noise + nuisance
  mk <- function(add_obs, add_nuisance) {
    lapply(seq_len(S), function(i) {
      V <- V0
      for (k in (m_shared + 1):r)
        V[, k] <- cos(theta) * V0[, k] +
          sin(theta) * Q[, r + (i - 1) * (r - m_shared) + (k - m_shared)]
      M <- Z[[i]] %*% Sig %*% t(V)
      if (add_obs) M <- M + sigma_obs * E[[i]]
      if (add_nuisance && mu > 0) {
        Wi <- sqrt(rho) * W + sqrt(1 - rho) * J[[i]]
        M <- M + sqrt(mu) * A[[i]] %*% diag(wts, q) %*% t(Wi)
      }
      M
    })
  }
  list(Xbio = mk(FALSE, FALSE), X0 = mk(TRUE, FALSE), Xn = mk(TRUE, TRUE),
       W = W, S = S, r = r, FX = FX, FY = FY, F = F,
       theta = theta, mu = mu, sigma_obs = sigma_obs, seed = seed)
}

couple <- function(Xlist, FX, FY) {
  big <- do.call(rbind, Xlist)
  idx <- rep(seq_along(Xlist), vapply(Xlist, nrow, integer(1)))
  sc  <- function(M) { s <- apply(M, 2, sd); s[s < 1e-12] <- 1e-12
  sweep(sweep(M, 2, colMeans(M), "-"), 2, s, "/") }
  bs  <- cbind(sc(big[, 1:FX, drop = FALSE]),
               sc(big[, (FX + 1):(FX + FY), drop = FALSE]))
  lapply(seq_along(Xlist), function(i) {
    Mi <- bs[idx == i, , drop = FALSE]
    nr <- sqrt(colSums(Mi^2)); nr[nr < 1e-12] <- 1e-12
    crossprod(sweep(Mi, 2, nr, "/"))
  })
}

# ---- One configuration, on contaminated and nuisance-free data ----
run_v2 <- function(theta, mu, sigma_obs, seed, n_eig = 50, verbose = FALSE,
                   qsub_fixed = NULL, ...) {
  sim <- gen_v2(theta, mu, sigma_obs, seed = seed, ...)
  r <- sim$r; S <- sim$S
  
  Bn  <- couple(sim$Xn,   sim$FX, sim$FY)   # contaminated
  B0  <- couple(sim$X0,   sim$FX, sim$FY)   # obs noise, no nuisance
  Bb  <- couple(sim$Xbio, sim$FX, sim$FY)   # signal only
  
  # biological reference subspace from Xbio
  Qstar <- qr.Q(qr(eigs_sym(Reduce(`+`, Bb), r, which = "LA")$vectors))
  
  # nuisance subspace Q_N, re-estimated per condition unless qsub_fixed is given
  qsub <- qsub_fixed
  if (mu > 0 && is.null(qsub)) {
    D  <- Reduce(`+`, lapply(seq_len(S), function(i) Bn[[i]] - B0[[i]])) / S
    ed <- eigs_sym(D, Q_NUIS + 2, which = "LM")
    # keep positive eigenvalues only (D can be indefinite)
    pos <- which(ed$values > 0)
    if (verbose) cat("   D spectrum:", round(ed$values, 4), "\n")
    qsub <- if (length(pos) >= 1)
      ed$vectors[, pos[seq_len(min(Q_NUIS, length(pos)))], drop = FALSE] else NULL
    if (!is.null(qsub) && verbose) {
      sd_all <- apply(do.call(rbind, sim$Xn), 2, sd)
      Ww <- sim$W / pmax(sd_all, 1e-12); Ww <- qr.Q(qr(Ww))
      cat("   overlap(qsub, scaled W):",
          round(sum((t(qsub) %*% Ww)^2) / ncol(qsub), 3), "\n")
    }
  }
  
  locfit <- function(mats) {
    out <- vector("list", S); rr <- numeric(S)
    for (i in seq_len(S)) {
      eo <- eigs_sym(mats[[i]], n_eig, which = "LA")
      ri <- max(2, min(find_elbow_kneedle(eo$values), n_eig))
      rr[i] <- ri
      out[[i]] <- list(V = eo$vectors[, 1:ri, drop = FALSE],
                       lam = eo$values[1:ri])
    }
    list(bases = out, r1 = rr)
  }
  
  # returns recovery for stage2_only and both, plus the selected global rank
  eval_pair <- function(Bset) {
    L  <- locfit(Bset)
    P1 <- lapply(seq_len(S), function(i)
      L$bases[[i]]$V %*% diag(L$bases[[i]]$lam, nrow = length(L$bases[[i]]$lam)) %*%
        t(L$bases[[i]]$V))
    one <- function(stage1) {
      agg <- if (stage1) Reduce(`+`, P1) else Reduce(`+`, Bset)
      vr  <- eigs_sym(agg, n_eig, which = "LA")
      d   <- max(2, find_elbow_kneedle(vr$values[1:n_eig]))
      G   <- vr$vectors[, 1:d, drop = FALSE]
      cap <- sum((t(Qstar) %*% qr.Q(qr(G)))^2)
      list(rec = cap / r,
           # purity: captured signal per retained dimension
           purity = cap / d,
           d = d, G = G)
    }
    list(s2 = one(FALSE), full = one(TRUE), L = L)
  }
  
  En <- eval_pair(Bn)     # contaminated
  E0 <- eval_pair(B0)     # nuisance-free
  
  pi_local <- if (!is.null(qsub))
    mean(vapply(En$L$bases, function(b)
      sum((t(b$V) %*% qsub)^2) / ncol(qsub), numeric(1))) else NA_real_
  # nuisance retained in the global basis of each arm
  pi_nuis_s2   <- if (!is.null(qsub))
    sum((t(En$s2$G)   %*% qsub)^2) / ncol(qsub) else NA_real_
  pi_nuis_full <- if (!is.null(qsub))
    sum((t(En$full$G) %*% qsub)^2) / ncol(qsub) else NA_real_
  
  Lj <- vapply(seq_len(min(En$s2$d, 20)), function(j)
    median(vapply(En$L$bases, function(b)
      sum((t(b$V) %*% En$s2$G[, j])^2), numeric(1))), numeric(1))
  
  data.frame(
    theta = theta, mu = mu, sigma_obs = sigma_obs, seed = seed,
    # local rank with and without nuisance
    r1_nuis = median(En$L$r1), r1_clean = median(E0$L$r1),
    # global rank for both arms
    d2_nuis = En$s2$d, d2_clean = E0$s2$d,
    d_full_nuis = En$full$d, d_full_clean = E0$full$d,
    pi_local = pi_local,
    pi_nuis_s2 = pi_nuis_s2, pi_nuis_full = pi_nuis_full,
    pi_global = pi_nuis_s2,
    pur_s2_nuis = En$s2$purity, pur_full_nuis = En$full$purity,
    pur_s2_clean = E0$s2$purity, pur_full_clean = E0$full$purity,
    # recovery on contaminated and clean data
    rec_s2_nuis = En$s2$rec, rec_full_nuis = En$full$rec,
    rec_s2_clean = E0$s2$rec, rec_full_clean = E0$full$rec,
    L_min = min(Lj), n_weak = sum(Lj < 0.5))
}

# ---- Reference nuisance subspace at one condition ----
qn_at <- function(mu, sigma_obs, theta = 0, seed = 1, ...) {
  sim <- gen_v2(theta, mu, sigma_obs, seed = seed, ...)
  Bn  <- couple(sim$Xn, sim$FX, sim$FY)
  B0  <- couple(sim$X0, sim$FX, sim$FY)
  D   <- Reduce(`+`, lapply(seq_len(sim$S), function(i) Bn[[i]] - B0[[i]])) / sim$S
  ed  <- eigs_sym(D, Q_NUIS + 2, which = "LM")
  pos <- which(ed$values > 0)
  list(Q = ed$vectors[, pos[seq_len(min(Q_NUIS, length(pos)))], drop = FALSE],
       leading = max(abs(ed$values)))
}

# ---- Gate 1: choose sigma_obs from the mu = 0 local rank ----
if (TRUE) {
  g1 <- do.call(rbind, lapply(c(0.5, 1, 2, 3, 5), function(so) {
    rr <- vapply(1:3, function(sd) {
      sim <- gen_v2(theta = 0.6, mu = 0, sigma_obs = so, seed = sd)
      B0  <- couple(sim$X0, sim$FX, sim$FY)
      median(vapply(B0, function(U)
        max(2, find_elbow_kneedle(eigs_sym(U, 50, which = "LA")$values)),
        numeric(1)))
    }, numeric(1))
    data.frame(sigma_obs = so, r1_mean = mean(rr), r1_sd = sd(rr),
               r1_all = paste(rr, collapse = ","))
  }))
  print(g1, row.names = FALSE)
  cat("\nGATE 1: pick the sigma_obs with r1 in [9, 12] and small sd.\n",
      "If none qualifies (e.g. Kneedle at 25-40),\n",
      "this construction cannot be used.\n")
}

# ---- Main sweep: mu in {0.03, ..., 2}, theta = 0, seeds 1:10 ----
SIGMA_OBS <- 0.5      # from Gate 1
if (TRUE) {
  grid <- expand.grid(theta = c(0),
                      mu = c(0.03, 0.06, 0.12, 0.25, 0.5, 1, 2), seed = 1:10)
  t0 <- Sys.time()
  res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(k) {
    cat(sprintf("  %d/%d  %.1f min\n", k, nrow(grid),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    utils::flush.console()
    run_v2(grid$theta[k], grid$mu[k], SIGMA_OBS, grid$seed[k],
           verbose = (k <= 3))     # print D spectrum for the first few cells
  }))
  write.csv(res, file.path(out_dir, "stress_v3b.csv"), row.names = FALSE)
}

# ---- Check A: conditioning of Q_N across mu ----
if (TRUE) {
  mus <- c(0.03, 0.06, 0.12, 0.25, 0.5, 1, 2)
  QQ  <- lapply(mus, function(m) qn_at(m, SIGMA_OBS))
  cat("\n mu   leading |D| eigenvalue   overlap with mu = 2\n")
  ref <- QQ[[length(QQ)]]$Q
  for (k in seq_along(mus))
    cat(sprintf("%5.2f  %14.3e  %14.3f\n", mus[k], QQ[[k]]$leading,
                sum((t(QQ[[k]]$Q) %*% ref)^2) / ncol(ref)))
  cat("\n adjacent-mu overlaps\n")
  for (k in seq_len(length(mus) - 1))
    cat(sprintf("  %.2f vs %.2f : %.3f\n", mus[k], mus[k + 1],
                sum((t(QQ[[k]]$Q) %*% QQ[[k + 1]]$Q)^2) / Q_NUIS))
  cat("\n Adjacent overlaps near 1 with a smooth decline against the mu = 2\n",
      "reference indicate rotation of the effective direction rather than\n",
      "instability of the estimate.\n")
}

# ---- Check B: fixed Q_N reference across mu ----

if (TRUE) {
  # reference per seed, at mu = 0.06
  qn_ref <- lapply(1:10, function(sd) qn_at(0.06, SIGMA_OBS, seed = sd)$Q)
  gridB  <- expand.grid(theta = 0,
                        mu = c(0.03, 0.06, 0.12, 0.25, 0.5, 1, 2), seed = 1:10)
  resB <- do.call(rbind, lapply(seq_len(nrow(gridB)), function(k) {
    cat(sprintf("  %d/%d\n", k, nrow(gridB))); utils::flush.console()
    run_v2(gridB$theta[k], gridB$mu[k], SIGMA_OBS, gridB$seed[k],
           qsub_fixed = qn_ref[[gridB$seed[k]]])
  }))
  write.csv(resB, file.path(out_dir, "stress_v3b_fixedref.csv"), row.names = FALSE)
  
  cmp <- dplyr::bind_rows(
    read.csv(file.path(out_dir, "stress_v3b.csv")) %>% dplyr::filter(theta == 0) %>%
      dplyr::mutate(reference = "per condition"),
    resB %>% dplyr::mutate(reference = "fixed"))
  print(as.data.frame(cmp %>% dplyr::group_by(reference, mu) %>%
                        dplyr::summarise(pi_local = round(mean(pi_local), 3),
                                         pi_S2    = round(mean(pi_nuis_s2), 3),
                                         pi_Full  = round(mean(pi_nuis_full), 3),
                                         removed  = round(mean(pi_nuis_s2 - pi_nuis_full), 3),
                                         .groups = "drop")), row.names = FALSE)
}

if (TRUE) {
  res <- read.csv(file.path(out_dir, "stress_v3b.csv"), stringsAsFactors = FALSE)
  
  g <- res %>% dplyr::group_by(theta, mu) %>% dplyr::summarise(
    n = dplyr::n(),
    r1_nuis = round(mean(r1_nuis), 1), r1_clean = round(mean(r1_clean), 1),
    shift = round(mean(r1_nuis - r1_clean), 2),
    d2_nuis = round(mean(d2_nuis), 1), d2_clean = round(mean(d2_clean), 1),
    dd2 = round(mean(d2_nuis - d2_clean), 2),
    pi_local = round(mean(pi_local), 4),
    pi_nuis_s2 = round(mean(pi_nuis_s2), 4),
    .groups = "drop") %>%
    dplyr::mutate(
      gate2 = abs(shift) <= 1,          # nuisance does not move the local rank
      gate3 = pi_local < 0.1 & pi_nuis_s2 > 0.5,
      # dd2 ~ 0: nuisance displaces biology; dd2 > 0: Stage 2 raises the rank to keep nuisance
      regime = ifelse(abs(dd2) < 0.5, "displacement", "expansion"),
      mechanism = gate2 & gate3)
  
  cat("\n=== mechanism gates, with the failure mode classified ===\n")
  print(as.data.frame(g), row.names = FALSE)
  cat("\ncells satisfying the mechanism conditions:", sum(g$mechanism), "\n")
}

# ---- Denoising read-out: nuisance retained, global rank, purity ----
if (TRUE) {
  res <- read.csv(file.path(out_dir, "stress_v3b.csv"), stringsAsFactors = FALSE)
  d <- res %>% dplyr::group_by(theta, mu) %>% dplyr::summarise(
    n = dplyr::n(),
    pi_s2   = round(mean(pi_nuis_s2), 4),
    pi_full = round(mean(pi_nuis_full), 4),
    d_nuis_removed = round(mean(pi_nuis_s2 - pi_nuis_full), 4),
    se_removed = round(sd(pi_nuis_s2 - pi_nuis_full) / sqrt(dplyr::n()), 4),
    t_removed  = round(mean(pi_nuis_s2 - pi_nuis_full) /
                         (sd(pi_nuis_s2 - pi_nuis_full) / sqrt(dplyr::n())), 2),
    d_s2 = round(mean(d2_nuis), 1), d_full = round(mean(d_full_nuis), 1),
    rec_s2 = round(mean(rec_s2_nuis), 4), rec_full = round(mean(rec_full_nuis), 4),
    pur_s2 = round(mean(pur_s2_nuis), 4), pur_full = round(mean(pur_full_nuis), 4),
    .groups = "drop")
  print(as.data.frame(d), row.names = FALSE)
  
  cat("\nExpected under the mechanism:\n",
      "  pi_full << pi_s2      Stage 1 removes what Stage 2 retains\n",
      "  rec_full not < rec_s2 without costing biological recovery\n",
      "  d_full < d_s2         and without the nuisance-driven rank inflation\n",
      "  pur_full > pur_s2     leaving a cleaner retained representation\n",
      "Read these only for cells meeting the mechanism conditions above.\n")
}

# ---- Gate 4: difference-in-differences of recovery ----
if (TRUE) {
  res <- read.csv(file.path(out_dir, "stress_v3b.csv"), stringsAsFactors = FALSE)
  dd <- res %>% dplyr::mutate(
    damage_clean = rec_full_clean - rec_s2_clean,
    loss_s2      = rec_s2_clean   - rec_s2_nuis,
    loss_full    = rec_full_clean - rec_full_nuis,
    protection   = loss_s2 - loss_full,
    raw_delta    = rec_full_nuis  - rec_s2_nuis)
  
  print(as.data.frame(dd %>% dplyr::group_by(theta, mu) %>% dplyr::summarise(
    n = dplyr::n(),
    s2_clean   = round(mean(rec_s2_clean), 4),
    full_clean = round(mean(rec_full_clean), 4),
    damage     = round(mean(damage_clean), 4),
    s2_nuis    = round(mean(rec_s2_nuis), 4),
    full_nuis  = round(mean(rec_full_nuis), 4),
    raw_delta  = round(mean(raw_delta), 4),
    mean_protection = round(mean(protection), 4),
    prot_se = round(sd(protection) / sqrt(dplyr::n()), 4),
    prot_t  = round(mean(protection) / (sd(protection) / sqrt(dplyr::n())), 2),
    .groups = "drop")), row.names = FALSE)
  
  cat("\nprotection ~ 0 does not mean Stage 1 has no denoising benefit:\n",
      "nuisance retention need not cause signal-subspace loss\n",
      "(expansion regime). Denoising is read from the pi_nuis\n",
      "table above.\n")
}

# ---- Figure: (A) nuisance overlap per arm, (B) paired recovery difference; theta = 0 ----
suppressPackageStartupMessages(library(patchwork))

res <- read.csv(file.path(out_dir, "stress_v3b.csv"), stringsAsFactors = FALSE)
d0  <- res %>%
  dplyr::filter(theta == 0) %>%
  dplyr::mutate(d_rec = rec_full_nuis - rec_s2_nuis)

MU_BREAKS <- c(0.03, 0.06, 0.12, 0.25, 0.5, 1, 2)

pi_cols <- c("Local retained (Stage 1)"         = "#4D4D4D",
             "Global retained (Stage 2 only)"   = "#009E73",
             "Global retained (full two-stage)" = "#D55E00")

pA_dat <- d0 %>%
  dplyr::select(mu, seed,
                `Local retained (Stage 1)`         = pi_local,
                `Global retained (Stage 2 only)`   = pi_nuis_s2,
                `Global retained (full two-stage)` = pi_nuis_full) %>%
  tidyr::pivot_longer(-c(mu, seed), names_to = "quantity", values_to = "y") %>%
  dplyr::mutate(quantity = factor(quantity, levels = names(pi_cols)))

# log2 axis for the doubling mu sequence
base_x <- scale_x_continuous(trans = "log2", breaks = MU_BREAKS,
                             labels = MU_BREAKS)

pA <- ggplot(pA_dat, aes(mu, y, colour = quantity)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point", size = 1.6) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.08) +
  scale_colour_manual(values = pi_cols) +
  base_x +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        legend.margin = margin(b = -4),
        plot.tag = element_text(size = 11, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 8))) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Fraction of nuisance\nsubspace represented",
       colour = NULL, tag = "A")

pB <- ggplot(d0, aes(mu, d_rec)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  stat_summary(fun = mean, geom = "line", colour = "black") +
  stat_summary(fun = mean, geom = "point", colour = "black", size = 1.6) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
               width = 0.08, colour = "black") +
  coord_cartesian(ylim = c(-0.002, 0.002)) +
  base_x +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        plot.tag = element_text(size = 11, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 8))) +
  labs(x = expression(paste("Nuisance strength ", mu)),
       y = expression(atop(Delta*" biological recovery",
                           "(Full - Stage 2 only)")),
       tag = "B")

fig <- pA / pB + plot_layout(heights = c(1.25, 1))

png(file.path(FIG_DIR, "figS11.png"),width = 7, height = 5, units = "in", res = 600)
print(fig)
dev.off()

# ---- Supplementary table ----
tab <- res %>%
  dplyr::mutate(d_rec = rec_full_nuis - rec_s2_nuis) %>%
  dplyr::group_by(theta, mu) %>%
  dplyr::summarise(
    d_clean    = round(mean(d2_clean), 1),
    d_S2only   = round(mean(d2_nuis), 1),
    d_Full     = round(mean(d_full_nuis), 1),
    pi_local   = round(mean(pi_local), 3),
    pi_S2only  = round(mean(pi_nuis_s2), 3),
    pi_Full    = round(mean(pi_nuis_full), 3),
    removed    = round(mean(pi_nuis_s2 - pi_nuis_full), 3),
    removed_pc = round(100 * mean(pi_nuis_s2 - pi_nuis_full) / mean(pi_nuis_s2)),
    rec_S2only = round(mean(rec_s2_nuis), 4),
    rec_Full   = round(mean(rec_full_nuis), 4),
    .groups = "drop")
print(as.data.frame(tab), row.names = FALSE)
write.csv(tab, file.path(out_dir, "supp_table_mechanism.csv"), row.names = FALSE)

# ---- Summary values at mu = 0.25 with 95% CI ----
ci <- function(x) {
  m <- mean(x); s <- sd(x) / sqrt(length(x))
  sprintf("%.4f (95%% CI %.4f to %.4f)", m, m - 1.96 * s, m + 1.96 * s)
}
rep_cond <- res %>% dplyr::mutate(d_rec = rec_full_nuis - rec_s2_nuis) %>%
  dplyr::filter(theta == 0, mu == 0.25)
cat("\nRepresentative condition, theta = 0, mu = 0.25, n =", nrow(rep_cond), "\n")
for (nm in c("pi_local", "pi_nuis_s2", "pi_nuis_full", "rec_s2_nuis",
             "rec_full_nuis", "d_rec"))
  cat(sprintf("  %-14s %s\n", nm, ci(rep_cond[[nm]])))
