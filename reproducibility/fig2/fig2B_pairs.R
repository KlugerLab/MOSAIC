# Figure 2B (right): per-sample correlations of nine feature pairs with high cross-sample variability
# Input:    results/fig2/panel_b_frobenius_results.RData
# Output:   results/fig2/panel_b_9pairs_results.RData, figures/fig2B_right.png
# Requires: fig2/fig2B_frobenius.R

source("config.R")

library(ggplot2)
library(matrixStats)
library(ggtext)        # italic gene symbols in axis labels

# ---- Load fig2B_frobenius.R results ----
load(file.path(RESULTS_DIR, "fig2", "panel_b_frobenius_results.RData"))   # real_cor_mats, sample_names, feature_names, det_frac, ...

MIN_DET_GENE   <- 0.05   # min detection fraction in every sample
MIN_DET_PEAK   <- 0.02
CIS_WINDOW_KB  <- 500    # peak-peak pairs: same chromosome, within 500 kb

is_peak <- grepl("^chr", feature_names)
min_det <- rowMins(det_frac)
eligible <- ifelse(is_peak, min_det >= MIN_DET_PEAK, min_det >= MIN_DET_GENE)
gene_idx <- which(!is_peak & eligible)
peak_idx <- which( is_peak & eligible)
cat("Eligible genes:", length(gene_idx), " eligible peaks:", length(peak_idx), "\n")

# per-sample correlations: rows = pairs, cols = samples
pair_cors <- function(i1, i2) sapply(real_cor_mats, function(cm) cm[cbind(i1, i2)])

# rank pairs by IQR across samples
pair_table <- function(i1, i2, type) {
  C <- pair_cors(i1, i2); if (is.null(dim(C))) C <- matrix(C, nrow = 1)
  data.frame(f1_idx = i1, f2_idx = i2, f1 = feature_names[i1], f2 = feature_names[i2],
             iqr = rowIQRs(C), sd = rowSds(C), min = rowMins(C), max = rowMaxs(C), type = type)
}

# pre-filter: features whose correlation profile varies most across samples (median IQR)
profile_iqr <- function(idx) sapply(idx, function(f) {
  P <- sapply(real_cor_mats, function(cm) cm[f, c(gene_idx, peak_idx)])
  median(rowIQRs(P))
})
top_gene_idx <- gene_idx[order(-profile_iqr(gene_idx))][1:min(500, length(gene_idx))]
top_peak_idx <- peak_idx[order(-profile_iqr(peak_idx))][1:min(300, length(peak_idx))]

select_diverse <- function(df, n = 3, exclude = character(0)) {
  seen <- exclude; sel <- c()
  for (i in seq_len(nrow(df))) {
    if (!(df$f1[i] %in% seen) && !(df$f2[i] %in% seen)) { seen <- c(seen, df$f1[i], df$f2[i]); sel <- c(sel, i) }
    if (length(sel) >= n) break
  }
  df[sel, ]
}

# ---- Gene-gene ----
gg <- combn(top_gene_idx, 2)
gg_df <- pair_table(gg[1, ], gg[2, ], "Gene – Gene")
gg_df <- gg_df[order(-gg_df$iqr), ]
gg_top <- select_diverse(gg_df, 3)

# ---- Gene-peak ----
gp <- expand.grid(g = top_gene_idx[1:min(200, length(top_gene_idx))], p = top_peak_idx)
gp_df <- pair_table(gp$g, gp$p, "Gene – Peak")
gp_df <- gp_df[order(-gp_df$iqr), ]
gp_top <- select_diverse(gp_df, 3, exclude = c(gg_top$f1, gg_top$f2))

# ---- Peak-peak (cis only) ----
pp <- combn(top_peak_idx, 2)
pp_df <- pair_table(pp[1, ], pp[2, ], "Peak – Peak")
pp_df$chr1 <- sub("-.*", "", pp_df$f1); pp_df$chr2 <- sub("-.*", "", pp_df$f2)
pp_df$pos1 <- as.numeric(sub("^[^-]+-(\\d+)-.*", "\\1", pp_df$f1))
pp_df$pos2 <- as.numeric(sub("^[^-]+-(\\d+)-.*", "\\1", pp_df$f2))
pp_df$distance_kb <- abs(pp_df$pos1 - pp_df$pos2) / 1000
pp_cis <- pp_df[pp_df$chr1 == pp_df$chr2 & pp_df$distance_kb <= CIS_WINDOW_KB, ]
cat("Cis peak-peak pairs available:", nrow(pp_cis), "\n")
if (nrow(pp_cis) < 3) stop("Fewer than 3 cis peak-peak pairs; increase top_peak_idx size or lower MIN_DET_PEAK.")
pp_cis <- pp_cis[order(-pp_cis$iqr), ]
pp_top <- select_diverse(pp_cis, 3, exclude = c(gg_top$f1, gg_top$f2, gp_top$f1, gp_top$f2))

all_top <- rbind(gg_top[, c("f1","f2","type","iqr","sd","min","max")],
                 gp_top[, c("f1","f2","type","iqr","sd","min","max")],
                 pp_top[, c("f1","f2","type","iqr","sd","min","max")])
print(all_top)

# ---- Labels: gene symbols italic, peaks roman ----
hgnc_update <- c("MARCH1" = "MARCHF1")   # outdated symbols to current HGNC
fmt_feature <- function(f) {
  if (grepl("^chr", f)) {
    parts <- strsplit(f, "-")[[1]]
    return(paste0(parts[1], ":", round(as.numeric(parts[2]) / 1000), "k"))
  }
  if (f %in% names(hgnc_update)) f <- hgnc_update[[f]]
  paste0("*", f, "*")
}

plot_data <- do.call(rbind, lapply(seq_len(nrow(all_top)), function(i) {
  i1 <- match(all_top$f1[i], feature_names); i2 <- match(all_top$f2[i], feature_names)
  data.frame(pair_label = paste0(fmt_feature(all_top$f1[i]), " – ", fmt_feature(all_top$f2[i])),
             pair_type = all_top$type[i], sample = sample_names,
             correlation = sapply(real_cor_mats, function(cm) cm[i1, i2]), iqr = all_top$iqr[i])
}))
plot_data$pair_type <- factor(plot_data$pair_type, levels = c("Gene – Gene", "Gene – Peak", "Peak – Peak"))
ord <- plot_data[!duplicated(plot_data$pair_label), ]
ord <- ord[order(ord$pair_type, -ord$iqr), ]
plot_data$pair_label <- factor(plot_data$pair_label, levels = rev(ord$pair_label))

type_colors <- c("Gene – Gene" = "#4A90D9", "Gene – Peak" = "#E67E22", "Peak – Peak" = "#2ECC71")

p <- ggplot(plot_data, aes(x = correlation, y = pair_label, color = pair_type)) +
  geom_point(size = 1.8, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2.5, color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = type_colors) +
  facet_grid(pair_type ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Pearson correlation (per sample)", y = NULL) +
  theme_classic() +
  theme(axis.text.y = element_markdown(size = 7, color = "black"),
        axis.title.x = element_text(size = 10),
        strip.text.y = element_text(size = 9, face = "bold"),
        strip.background = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0.8, "lines"))

png(file.path(FIG_DIR, "fig2B_right.png"), width = 6, height = 4, units = "in", res = 600)
print(p)
dev.off()

save(all_top, gg_df, gp_df, pp_df, plot_data, file = file.path(RESULTS_DIR, "fig2", "panel_b_9pairs_results.RData"))
cat("Results saved.\n")
