# Supplemental Figure S3: Stage 4, negative/positive controls and formal benchmark
# Input:    s3_frozen.rds (--mode formal)
# Output:   stage4/formal_auprc.csv, stage4/formal_mechanism.csv, stage4/formal_all.rds
# Requires: s3_stage3_calibrate.R
# Usage:    Rscript s3_stage4_run.R --mode controls --workers 70
#           Rscript s3_stage4_run.R --mode formal --workers 70
S3_DIR <- if (nzchar(Sys.getenv("S3_DIR"))) Sys.getenv("S3_DIR") else getwd()
for (f in c("s3_config.R", "s3_generator.R", "s3_methods.R", "s3_metrics.R", "s3_pipeline.R"))
  source(file.path(S3_DIR, f))
suppressPackageStartupMessages({ library(foreach); library(doParallel) })

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default, cast = identity) {
  i <- match(name, args); if (is.na(i)) default else cast(args[i + 1])
}
mode    <- arg("--mode", "formal")
workers <- arg("--workers", S3_CFG$n_workers, as.integer)
cfg <- S3_CFG
if ("--quick" %in% args) {                                    # smoke test only
  cfg$F <- 40L; cfg$r <- 10L; cfg$n_m <- 4L; cfg$m <- 3L; cfg$n_cells <- 80L; cfg$K_bg <- 200L
  cfg$memento_chunks <- 2L; cfg$n_eig <- 20L
}
stopifnot(mode %in% c("controls", "formal"))

if (mode == "controls") {
  reps <- arg("--reps", 4L, as.integer)
  seeds <- 40L + seq_len(reps)                                # pilot seeds
  settings <- rbind(data.frame(setting = "negative", panel = "control", kappa = 0.5, beta2 = 0, eta = 0),
                    data.frame(setting = "positive", panel = "control", kappa = 1, beta2 = beta2_max(cfg), eta = 0))
  d_fixed <- NULL
} else {
  fz_path <- "s3_frozen.rds"
  if (!file.exists(fz_path)) stop("s3_frozen.rds not found: run Stage 3 first")
  fz <- readRDS(fz_path)
  if (!identical(fz$hash, digest::digest(fz[setdiff(names(fz), "hash")])))
    stop("s3_frozen.rds hash mismatch: the frozen configuration was modified")
  cfg <- fz$cfg
  reps <- arg("--reps", 10L, as.integer)
  seeds <- cfg$final_seed0 + seq_len(reps) - 1L
  settings <- rbind(
    data.frame(setting = paste0("A_k", fz$kappa_grid), panel = "A", kappa = fz$kappa_grid,
               beta2 = unname(fz$beta2[as.character(fz$kappa_grid)]), eta = 0),
    data.frame(setting = paste0("B_e", cfg$etas), panel = "B", kappa = cfg$kappa_eta_panel,
               beta2 = unname(fz$beta2[as.character(cfg$kappa_eta_panel)]), eta = cfg$etas))
  d_fixed <- fz$d_fixed
}
jobs <- merge(settings, data.frame(seed = seeds), by = NULL)
jobs$id <- sprintf("%s_s%d", jobs$setting, jobs$seed)
cat(sprintf("%s: %d cohorts (%d settings x %d reps)\n", mode, nrow(jobs), nrow(settings), reps))
t0 <- proc.time()[["elapsed"]]

# Non-memento arms (fork workers)
registerDoParallel(cores = min(workers, nrow(jobs)))
nonm <- foreach(k = seq_len(nrow(jobs))) %dopar% {
  j <- jobs[k, ]
  run_cohort_nonmemento(cfg, j$seed, j$kappa, j$beta2, j$eta, d_fixed = d_fixed)
}
names(nonm) <- jobs$id
stopImplicitCluster()
cat(sprintf("non-memento arms done: %.1f s\n", proc.time()[["elapsed"]] - t0))

# memento: pair chunks on PSOCK workers
py <- setup_memento_env(cfg)
mem <- run_memento_chunked(lapply(seq_len(nrow(jobs)), function(k) as.list(jobs[k, c("id", "seed", "kappa", "beta2", "eta")])),
                           cfg, workers, python = py)
cat(sprintf("memento done: %.1f s total\n", proc.time()[["elapsed"]] - t0))

# Metrics
auprc <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(k) {
  x <- nonm[[jobs$id[k]]]; sc <- c(x$scores, list(memento = unname(mem[[jobs$id[k]]])))
  data.frame(jobs[k, c("id", "setting", "panel", "kappa", "beta2", "eta", "seed")],
             method = names(sc), auprc = vapply(sc, auprc_fig3, 0, labels = x$labels), row.names = NULL)
}))
mech <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(k) {
  x <- nonm[[jobs$id[k]]]; ok <- is.finite(x$G)
  data.frame(jobs[k, c("id", "setting", "panel", "kappa", "beta2", "eta", "seed")], x$diag,
             selectivity_auroc = auroc(log(x$G[ok]), x$labels[ok]),
             dc_percentile_shift = dc_percentile_shift(x$scores$raw, x$scores$mosaic, x$labels),
             delta_auprc = auprc_fig3(x$scores$mosaic, x$labels) - auprc_fig3(x$scores$raw, x$labels),
             row.names = NULL)
}))
dir.create("stage4", showWarnings = FALSE)
write.csv(auprc, sprintf("stage4/%s_auprc.csv", mode), row.names = FALSE)
write.csv(mech, sprintf("stage4/%s_mechanism.csv", mode), row.names = FALSE)
saveRDS(list(jobs = jobs, nonmemento = lapply(nonm, function(x) x[c("scores", "labels", "G", "diag")]),
             memento = mem, cfg = cfg), sprintf("stage4/%s_all.rds", mode))

summ <- aggregate(auprc ~ setting + method, data = auprc, FUN = mean)
print(reshape(summ, idvar = "setting", timevar = "method", direction = "wide"), digits = 3)
cat(sprintf("\ntotal wall time: %.1f min\n", (proc.time()[["elapsed"]] - t0) / 60))
if (mode == "controls") {
  main <- c("raw", "mosaic", "dgca", "memento")
  neg <- summ[summ$setting == "negative" & summ$method %in% main, ]
  pos <- summ[summ$setting == "positive" & summ$method %in% main, ]
  cat(sprintf("negative control (mean AUPRC < 0.25 for all four): %s\n", all(neg$auprc < 0.25)))
  cat(sprintf("positive control (mean AUPRC > 0.90 for all four): %s\n", all(pos$auprc > 0.90)))
}
