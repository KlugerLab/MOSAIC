# Supplemental Figure S3: Stage 4b, memento score from corr_coef only (not plotted)
# Input:    s3_frozen.rds, stage4/formal_auprc.csv
# Output:   stage4/sensitivity_memento_coef.csv, stage4/sensitivity_summary.csv
# Requires: s3_stage4_run.R --mode formal
# Usage:    Rscript s3_stage4b_memento_coef.R --workers 70
S3_DIR <- if (nzchar(Sys.getenv("S3_DIR"))) Sys.getenv("S3_DIR") else getwd()
for (f in c("s3_config.R", "s3_generator.R", "s3_methods.R", "s3_metrics.R"))
  source(file.path(S3_DIR, f))
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default, cast = identity) {
  i <- match(name, args); if (is.na(i)) default else cast(args[i + 1])
}
workers <- arg("--workers", S3_CFG$n_workers, as.integer)
boot    <- arg("--boot", 10L, as.integer)
chunks  <- arg("--chunks", 2L, as.integer)

fz <- readRDS("s3_frozen.rds")
if (!identical(fz$hash, digest::digest(fz[setdiff(names(fz), "hash")]))) stop("s3_frozen.rds hash mismatch")
cfg <- fz$cfg; cfg$memento_boot <- boot
formal <- read.csv("stage4/formal_auprc.csv")
jobs <- unique(formal[, c("id", "setting", "panel", "kappa", "beta2", "eta", "seed")])
cat(sprintf("memento coefficient-only rerun: %d cohorts, num_boot = %d, %d chunks\n", nrow(jobs), boot, chunks))

t0 <- proc.time()[["elapsed"]]
py <- setup_memento_env(cfg)
sc <- run_memento_chunked(lapply(seq_len(nrow(jobs)), function(k) as.list(jobs[k, c("id", "seed", "kappa", "beta2", "eta")])),
                          cfg, workers, python = py, n_chunks = chunks, stat = "coef")
labels <- make_cohort_parameters(cfg, jobs$seed[1])$is_dc                # DC set = fixed modules
per <- data.frame(jobs, method = "memento_coef",
                  auprc = vapply(jobs$id, function(i) auprc_fig3(unname(sc[[i]]), labels), 0), row.names = NULL)
dir.create("stage4", showWarnings = FALSE)
write.csv(per, "stage4/sensitivity_memento_coef.csv", row.names = FALSE)

all <- rbind(formal[formal$method %in% c("mosaic", "dgca", "memento"), names(per)], per)
ci <- function(x) { n <- length(x); m <- mean(x); h <- qt(0.975, n - 1) * sd(x) / sqrt(n)
                    sprintf("%.3f [%.3f, %.3f]", m, m - h, m + h) }
order_set <- unique(jobs$setting[order(jobs$panel, jobs$kappa, jobs$eta)])
summ <- do.call(rbind, lapply(order_set, function(s) {
  d <- all[all$setting == s, ]
  data.frame(setting = s, kappa = d$kappa[1], eta = d$eta[1],
             MOSAIC = ci(d$auprc[d$method == "mosaic"]),
             DGCA = ci(d$auprc[d$method == "dgca"]),
             memento_primary_coef_over_se = ci(d$auprc[d$method == "memento"]),
             memento_secondary_coef = ci(d$auprc[d$method == "memento_coef"]))
}))
write.csv(summ, "stage4/sensitivity_summary.csv", row.names = FALSE)
print(summ, right = FALSE)
cat(sprintf("\ntotal wall time: %.1f min\n", (proc.time()[["elapsed"]] - t0) / 60))
