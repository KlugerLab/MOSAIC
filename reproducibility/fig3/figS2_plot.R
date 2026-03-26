# Supplemental Figure S2: null calibration figure
# Input:    results/S2/null_calib_final/null_calib_results.rds
# Output:   figures/figS2.png
# Requires: fig3/figS2_null_calibration.R

source("config.R")

## A: type I error per cohort with mean and 95% CI; B: QQ of null p-values at --qq_sigma

suppressPackageStartupMessages({
  library(optparse); library(ggplot2); library(dplyr); library(patchwork)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--results", type = "character",
              default = file.path(RESULTS_DIR, "S2", "null_calib_final", "null_calib_results.rds")),
  make_option("--fig_dir", type = "character", default = FIG_DIR),
  make_option("--qq_sigma", type = "double", default = 0.5),
  make_option("--alpha", type = "double", default = 0.05),
  make_option("--B", type = "integer", default = 999L, help = "Number of label permutations"),
  make_option("--y_max", type = "double", default = 0.10, help = "Panel A y-axis upper limit")
)))
dir.create(opt$fig_dir, showWarnings = FALSE, recursive = TRUE)
store <- readRDS(opt$results)
col <- "#2E6F9E"

## Summaries
metrics <- do.call(rbind, lapply(store, `[[`, "metrics"))
metrics <- metrics[metrics$arm == "new", ]
ci <- function(x) { n <- length(x); if (n < 2) return(c(NA, NA))
  mean(x) + c(-1, 1) * qt(0.975, n - 1) * sd(x) / sqrt(n) }
summ <- metrics %>% group_by(sigma) %>%
  summarise(alpha_mean = mean(alpha_hat), alpha_lo = ci(alpha_hat)[1],
            alpha_hi = ci(alpha_hat)[2], .groups = "drop")
sig_lab <- function(x) format(x, nsmall = 1)
lev <- sig_lab(sort(unique(metrics$sigma)))
metrics$sigma_f <- factor(sig_lab(metrics$sigma), levels = lev)
summ$sigma_f    <- factor(sig_lab(summ$sigma), levels = lev)

## Panel A
p1 <- ggplot(summ, aes(x = sigma_f, y = alpha_mean)) +
  geom_hline(yintercept = opt$alpha, linetype = 3) +
  geom_point(data = metrics, aes(x = sigma_f, y = alpha_hat), colour = col,
             position = position_jitter(width = 0.05, height = 0, seed = 1),
             alpha = 0.35, size = 1.2) +
  geom_pointrange(aes(ymin = alpha_lo, ymax = alpha_hi), colour = col,
                  size = 0.35, linewidth = 0.6) +
  scale_y_continuous(limits = c(0, opt$y_max), expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression(sigma[sample]),
       y = sprintf("Empirical type I error (p < %.2f)", opt$alpha)) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

## Panel B
qq_keys <- names(store)[abs(sapply(store, function(x) x$het$sigma) - opt$qq_sigma) < 1e-9]
stopifnot(length(qq_keys) > 0)
qq <- do.call(rbind, lapply(qq_keys, function(k) {
  p <- store[[k]]$pvals[["new"]]
  data.frame(cohort = k, obs = -log10(sort(p)), exp = -log10(ppoints(length(p))))
}))
floor_y <- -log10(1 / (opt$B + 1))
p2 <- ggplot(qq, aes(exp, obs, group = cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  geom_hline(yintercept = floor_y, linetype = 3, colour = "grey55") +
  geom_point(size = 0.5, alpha = 0.6, colour = col) +
  labs(x = expression(Expected ~ -log[10](p)), y = expression(Observed ~ -log[10](p))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

## Combine
combo <- (p1 | p2) + plot_annotation(tag_levels = "A")

png(file.path(opt$fig_dir, "figS2.png"), width = 7, height = 4, units = "in", res = 600)
print(combo)
dev.off()

print(as.data.frame(summ), digits = 3)
