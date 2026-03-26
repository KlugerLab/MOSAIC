# Supplemental Figure S3: Stage 3, difficulty calibration and frozen configuration
# Output:   s3_frozen.rds, s3_frozen.txt, stage3/*.csv
# Usage:    Rscript s3_stage3_calibrate.R --workers 70
S3_DIR <- if (nzchar(Sys.getenv("S3_DIR"))) Sys.getenv("S3_DIR") else getwd()
for (f in c("s3_config.R", "s3_generator.R", "s3_methods.R", "s3_metrics.R", "s3_pipeline.R"))
  source(file.path(S3_DIR, f))
suppressPackageStartupMessages({ library(foreach); library(doParallel) })

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default, cast = identity) {
  i <- match(name, args); if (is.na(i)) default else cast(args[i + 1])
}
workers <- arg("--workers", S3_CFG$n_workers, as.integer)
cfg <- S3_CFG
calib_seeds <- 21:30; ref_seeds <- 21:24; check_seeds <- 25:28; align_seeds <- 21:24
fine_kappa <- c(0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 1)
window <- c(0.5, 0.7); aim <- 0.6
if ("--quick" %in% args) {                                     # smoke test only
  cfg$F <- 60L; cfg$r <- 10L; cfg$n_m <- 6L; cfg$m <- 3L; cfg$n_cells <- 100L; cfg$K_bg <- 200L
  cfg$ref_cells <- 400L; cfg$ref_donors <- 2L; cfg$n_eig <- 30L
  calib_seeds <- 21:23; align_seeds <- 21L; fine_kappa <- c(0, 0.05, 0.1, 0.3, 0.5, 1)
}
# acceptance checks stop the run (warn only with --quick)
chk_fail <- function(msg) if ("--quick" %in% args) warning(msg, call. = FALSE) else stop(msg, call. = FALSE)
registerDoParallel(cores = workers)
dir.create("stage3", showWarnings = FALSE)
# all CRN references in Stage 3 use twice the donors of the default reference
cfg_A <- cfg; cfg_A$ref_donors <- 2L * cfg$ref_donors
t0 <- proc.time()[["elapsed"]]
b2max <- beta2_max(cfg)

# A. Calibrate beta2 at kappa = 0.5 so that the unprojected-arm AUPRC is in [0.5, 0.7]
auprc_B <- function(b2) {
  a <- foreach(s = calib_seeds, .combine = c) %dopar% {
    x <- run_cohort_nonmemento(cfg, s, 0.5, b2, 0, raw_only = TRUE)
    auprc_fig3(x$scores$raw, x$labels)
  }
  mean(a)
}
# cache step A per configuration
cache_A <- "stage3/A_result.rds"
cfg_key <- digest::digest(list(cfg, calib_seeds, window, aim))
cached <- if (file.exists(cache_A)) readRDS(cache_A) else NULL
if (!is.null(cached) && identical(cached$key, cfg_key)) {
  log_A <- cached$log_A; b2_05 <- cached$b2_05
  cat(sprintf("A: reusing cached calibration (beta2(0.5) = %.4f)\n", b2_05))
} else {
log_A <- data.frame()
hi <- b2max; a_hi <- auprc_B(hi); log_A <- rbind(log_A, data.frame(beta2 = hi, auprc_B = a_hi))
if (a_hi < window[1]) stop(sprintf("arm B AUPRC %.3f < %.1f even at beta2_max: apply the knob order", a_hi, window[1]))
lo <- 0; b2_05 <- NA
if (a_hi <= window[2]) b2_05 <- hi
for (it in seq_len(14)) {
  if (!is.na(b2_05)) break
  mid <- (lo + hi) / 2; a <- auprc_B(mid)
  log_A <- rbind(log_A, data.frame(beta2 = mid, auprc_B = a))
  cat(sprintf("A: beta2 = %.4f  arm-B AUPRC = %.3f\n", mid, a))
  if (a >= window[1] && a <= window[2]) b2_05 <- mid else if (a < aim) lo <- mid else hi <- mid
}
if (is.na(b2_05)) stop("calibration did not reach the AUPRC window in 14 steps")
saveRDS(list(key = cfg_key, log_A = log_A, b2_05 = b2_05), cache_A)
}
write.csv(log_A, "stage3/A_calibration.csv", row.names = FALSE)
T_target <- mean(foreach(s = ref_seeds, .combine = c) %dopar% mean_norm_over_seeds(cfg_A, s, 0.5, b2_05))
cat(sprintf("A: beta2(0.5) = %.4f, arm-B AUPRC = %.3f, T = %.4f\n", b2_05, tail(log_A$auprc_B, 1), T_target))

# B. Capacity at kappa = 0
cap0 <- mean(foreach(s = ref_seeds, .combine = c) %dopar% mean_norm_over_seeds(cfg_A, s, 0, b2max))
if (T_target > cap0) stop(sprintf("T = %.3f exceeds kappa = 0 capacity %.3f: apply the knob order", T_target, cap0))

# C. beta2(kappa) on the fine grid: ||Delta|| at 5 beta values, quadratic fit, solve for T
fracs <- c(0.25, 0.4, 0.55, 0.7, 1)
cgrid <- expand.grid(kappa = setdiff(fine_kappa, 0.5), frac = fracs, seed = ref_seeds)
cvals <- foreach(k = seq_len(nrow(cgrid)), .combine = c) %dopar% {
  b <- sqrt(b2max) * cgrid$frac[k]
  r <- reference_delta_crn(cfg_A, cgrid$seed[k], cgrid$kappa[k], b^2)
  median(r$norm_corr[r$par$is_dc])
}
cgrid$norm <- cvals
write.csv(cgrid, "stage3/C_norm_grid.csv", row.names = FALSE)
cfit <- aggregate(norm ~ kappa + frac, data = cgrid, FUN = mean)
b2_fine <- sapply(fine_kappa, function(k) {
  if (k == 0.5) return(b2_05)
  d <- cfit[cfit$kappa == k, ]; b <- sqrt(b2max) * d$frac
  fit <- lm(d$norm ~ 0 + b + I(b^2)); a <- coef(fit)[[1]]; cc <- coef(fit)[[2]]
  # fit error relative to T
  rel <- max(abs(residuals(fit))) / T_target
  if (rel >= 0.05) chk_fail(sprintf("quadratic fit at kappa = %.2f: max |residual| / T = %.3f", k, rel))
  if (T_target > max(d$norm)) stop(sprintf("T = %.3f exceeds capacity at kappa = %.2f", T_target, k))
  bb <- if (abs(cc) < 1e-12) T_target / a else (-a + sqrt(a^2 + 4 * cc * T_target)) / (2 * cc)
  bb^2
})
names(b2_fine) <- as.character(fine_kappa)

# D. Alignment sweep (pooled A over the DC set, median over pilot seeds) and kappa grid
grid <- expand.grid(kappa = fine_kappa, seed = align_seeds)
al <- foreach(k = seq_len(nrow(grid)), .combine = rbind) %dopar% {
  kp <- grid$kappa[k]; s <- grid$seed[k]; b2 <- b2_fine[[as.character(kp)]]
  x <- run_cohort_nonmemento(cfg, s, kp, b2, 0, skip_dgca = TRUE)
  ref <- reference_delta_crn(cfg_A, s, kp, b2)
  data.frame(kappa = kp, seed = s, beta2 = b2, d = x$diag$d,
             pooled_A = reference_alignment_pooled(ref, x$Q, ref$par$is_dc),
             median_A = median(reference_alignment(ref, x$Q)[ref$par$is_dc]),
             median_norm = median(ref$norm_corr[ref$par$is_dc]))
}
write.csv(al, "stage3/D_alignment_sweep.csv", row.names = FALSE)
A_k <- tapply(al$pooled_A, al$kappa, median)                  # selection statistic: pooled A
A0 <- A_k[["0"]]; A1 <- A_k[["1"]]
targets <- A0 + (1:3) / 4 * (A1 - A0)
interior_pool <- fine_kappa[fine_kappa > 0 & fine_kappa < 1]
picked <- numeric(0)
for (a in targets) {
  pool <- setdiff(interior_pool, picked)
  picked <- c(picked, pool[which.min(abs(A_k[as.character(pool)] - a))])
}
kappa_grid <- sort(c(0, picked, 1))
if (length(unique(kappa_grid)) != 5) stop("kappa grid selection produced duplicates")
cat(sprintf("D: A(0) = %.3f, A(1) = %.3f, range %.3f%s\n", A0, A1, A1 - A0,
            if (A1 - A0 < 0.15) " (narrow; kept by rule)" else ""))
cat("D: kappa grid =", kappa_grid, " A =", round(A_k[as.character(kappa_grid)], 3), "\n")

# E. Check ||Delta|| on independent seeds, d_fixed, freeze
eg <- expand.grid(kappa = kappa_grid, seed = check_seeds)
eg$norm <- foreach(k = seq_len(nrow(eg)), .combine = c) %dopar%
  mean_norm_over_seeds(cfg_A, eg$seed[k], eg$kappa[k], b2_fine[[as.character(eg$kappa[k])]])
chk <- tapply(eg$norm, eg$kappa, mean)[as.character(kappa_grid)]
dev <- max(abs(chk - mean(chk)) / mean(chk))
cat(sprintf("E: ||Delta|| (mean over independent seeds %s) across grid: %s (max rel dev %.3f)\n",
            paste(check_seeds, collapse = ","), paste(round(chk, 3), collapse = ", "), dev))
if (dev >= 0.05) chk_fail("||Delta|| equality across kappa fails the 5% criterion")
d_fixed <- as.integer(round(median(al$d[al$kappa == 0.5])))

beta2_frozen <- b2_fine[as.character(sort(unique(c(kappa_grid, cfg$kappa_eta_panel))))]
fz <- list(cfg = cfg, T = T_target, beta2 = beta2_frozen, kappa_grid = kappa_grid,
           A_fine = A_k, A_targets = targets, d_fixed = d_fixed,
           calibration = list(beta2_05 = b2_05, arm_B_auprc = tail(log_A$auprc_B, 1), capacity_k0 = cap0),
           norm_check = setNames(chk, kappa_grid), created = format(Sys.time()))
fz$hash <- digest::digest(fz)
saveRDS(fz, "s3_frozen.rds")
writeLines(c(sprintf("frozen %s  hash %s", fz$created, fz$hash),
             sprintf("T = %.4f (arm-B AUPRC %.3f at kappa 0.5)", T_target, fz$calibration$arm_B_auprc),
             sprintf("kappa grid: %s", paste(kappa_grid, collapse = ", ")),
             sprintf("pooled A: %s", paste(round(A_k[as.character(kappa_grid)], 3), collapse = ", ")),
             sprintf("beta2: %s", paste(sprintf("%s=%.4f", names(beta2_frozen), beta2_frozen), collapse = ", ")),
             sprintf("d_fixed = %d", d_fixed)), "s3_frozen.txt")
cat(readLines("s3_frozen.txt"), sep = "\n")
cat(sprintf("\nStage 3 done in %.1f min. Configuration FROZEN.\n", (proc.time()[["elapsed"]] - t0) / 60))
