# Supplemental Figure S3: Stage 3c, predicted selectivity of the projection on pilot cohorts (not plotted)
# Input:    s3_frozen.rds
# Output:   stage3/predictions.csv, stage3/predictions.txt
# Requires: s3_stage3_calibrate.R
# Usage:    Rscript s3_stage3c_predict.R --workers 70
S3_DIR <- if (nzchar(Sys.getenv("S3_DIR"))) Sys.getenv("S3_DIR") else getwd()
for (f in c("s3_config.R", "s3_generator.R", "s3_methods.R", "s3_metrics.R", "s3_pipeline.R"))
  source(file.path(S3_DIR, f))
suppressPackageStartupMessages({ library(foreach); library(doParallel) })
args <- commandArgs(trailingOnly = TRUE)
workers <- if (is.na(match("--workers", args))) S3_CFG$n_workers else as.integer(args[match("--workers", args) + 1])

fz <- readRDS("s3_frozen.rds")
if (!identical(fz$hash, digest::digest(fz[setdiff(names(fz), "hash")]))) stop("s3_frozen.rds hash mismatch")
cfg <- fz$cfg; cfg_A <- cfg; cfg_A$ref_donors <- 2L * cfg$ref_donors
settings <- rbind(
  data.frame(setting = paste0("A_k", fz$kappa_grid), kappa = fz$kappa_grid,
             beta2 = unname(fz$beta2[as.character(fz$kappa_grid)]), eta = 0),
  data.frame(setting = paste0("B_e", cfg$etas), kappa = cfg$kappa_eta_panel,
             beta2 = unname(fz$beta2[as.character(cfg$kappa_eta_panel)]), eta = cfg$etas))
jobs <- merge(settings, data.frame(seed = 31:34), by = NULL)
registerDoParallel(cores = workers)
res <- foreach(k = seq_len(nrow(jobs)), .combine = rbind) %dopar% {
  j <- jobs[k, ]
  x <- run_cohort_nonmemento(cfg, j$seed, j$kappa, j$beta2, j$eta, skip_dgca = TRUE)
  ref <- reference_delta_crn(cfg_A, j$seed, j$kappa, j$beta2)
  dc <- x$labels
  # first-moment prediction of the multiplier, c = 1/n_A + 1/n_B
  S <- x$n_A + x$n_B; cc <- 1 / x$n_A + 1 / x$n_B
  sig_s_W <- rowSums((ref$D1 %*% x$W) * (ref$D2 %*% x$W))       # unbiased Delta W W' Delta'
  sig_s   <- rowSums(ref$D1 * ref$D2)                           # unbiased ||Delta||^2
  rs_expected <- (sig_s_W + cc * x$ssw_W / (S - 2)) / (sig_s + cc * x$ssw_raw / (S - 2))
  g_pred <- (rs_expected / x$R_N)[dc]; g_null <- x$G[!dc]
  okp <- is.finite(g_pred) & g_pred > 0; okn <- is.finite(g_null) & g_null > 0
  lp <- log(g_pred[okp]); ln <- log(g_null[okn])
  data.frame(j, pooled_A = reference_alignment_pooled(ref, x$Q, ref$par$is_dc),
             pred_selectivity = auroc(c(lp, ln), c(rep(TRUE, length(lp)), rep(FALSE, length(ln)))),
             median_log_Gpred_dc = median(lp), median_log_Ghat_null = median(ln),
             dc_dropped_frac = mean(!okp), row.names = NULL)
}
dir.create("stage3", showWarnings = FALSE)
summ <- aggregate(cbind(pooled_A, pred_selectivity, median_log_Gpred_dc, median_log_Ghat_null) ~ setting + kappa + eta,
                  data = res, FUN = median)
summ$prediction <- ifelse(summ$pred_selectivity > 0.5, "MOSAIC above unprojected profile",
                          "MOSAIC not above unprojected profile")
write.csv(res, "stage3/predictions_by_cohort.csv", row.names = FALSE)
write.csv(summ, "stage3/predictions.csv", row.names = FALSE)
writeLines(c(sprintf("predictions written %s for frozen configuration %s", format(Sys.time()), fz$hash),
             capture.output(print(summ, digits = 3))), "stage3/predictions.txt")
cat(readLines("stage3/predictions.txt"), sep = "\n")
