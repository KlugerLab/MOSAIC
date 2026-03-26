# Supplemental Figure S4: MDS of the sample embeddings of the top DC features
# Input:    results/fig4/dc_corrected_paired.rds, results/fig4/mosaic_embedding.rds
# Output:   figures/figS4.png
# Requires: fig4/fig4_dc_paired_test.R

source("config.R")

## Top 12 features by pseudo-F, or top 10 RNA and top 2 ADT if the top 12 contain no ADT feature

library(ggplot2); library(grid)

res_dir  <- file.path(RESULTS_DIR, "fig4")
plot_dir <- FIG_DIR
dc <- readRDS(file.path(res_dir, "dc_corrected_paired.rds"))
mosaic <- readRDS(file.path(res_dir, "mosaic_embedding.rds"))     # MOSAIC embedding

n_row <- 3; n_col <- 4
page_w  <- 6.85                   # page width, inches
panel_h <- 1.8                    # height of one panel, inches
legend_h <- 0.3

day_cols <- c("Day 0" = "#56B4E9", "Day 7" = "#E69F00")
minus <- function(x) sub("^-", "\u2212", format(x, trim = TRUE, scientific = FALSE))

## 1. Feature selection
tab <- dc$DC_table[order(-dc$DC_table$F_stat), ]
top12 <- head(tab, 12)
if (any(top12$Modality == "ADT")) {
  sel <- top12
  rule <- "top 12 features by pseudo-F (includes ADT)"
} else {
  sel <- rbind(head(tab[tab$Modality == "RNA", ], 10),
               head(tab[tab$Modality == "ADT", ], 2))
  rule <- "top 10 RNA and top 2 ADT features by pseudo-F"
}
cat("selection rule:", rule, "\n")
print(sel[, c("Feature", "Modality", "F_stat", "F_rank", "exact_p", "at_p_floor")], row.names = FALSE)
stopifnot(nrow(sel) == n_row * n_col)

## 2. Sample embeddings
emb <- mosaic$projected_combed_list[dc$sample_ids]
stopifnot(length(emb) == 16L, identical(names(emb), dc$sample_ids))
for (i in seq_along(emb)) rownames(emb[[i]]) <- dc$feature_names
day <- factor(ifelse(dc$cond == "0", "Day 0", "Day 7"), levels = c("Day 0", "Day 7"))

pseudoF <- function(X, g) {
  D2 <- as.matrix(dist(X))^2; S <- nrow(X)
  SS_T <- sum(D2) / (2 * S)
  SS_W <- sum(sapply(split(seq_len(S), g), function(ix) sum(D2[ix, ix]) / (2 * length(ix))))
  (SS_T - SS_W) / (SS_W / (S - 2))
}

base_theme <- theme_classic(base_size = 8) +
  theme(axis.line = element_line(linewidth = 0.35),
        axis.ticks = element_line(linewidth = 0.35),
        plot.title = element_text(size = 8.5, hjust = 0.5, margin = margin(b = 2)),
        legend.position = "none",
        plot.margin = margin(4, 6, 4, 4))

## 3. One MDS panel per feature, as in Figure 4B
mds_panel <- function(feat, name, modality, F_ref) {
  X <- t(sapply(emb, function(m) m[feat, ]))
  stopifnot(abs(pseudoF(X, dc$cond) - F_ref) < 1e-6)       # same distances as the test
  mds <- cmdscale(dist(X), k = 2, eig = TRUE)
  pct <- 100 * mds$eig[1:2] / sum(mds$eig[mds$eig > 0])
  cd <- data.frame(m1 = mds$points[, 1], m2 = mds$points[, 2], day = day)
  if (mean(cd$m1[cd$day == "Day 7"]) < mean(cd$m1[cd$day == "Day 0"])) cd$m1 <- -cd$m1
  cen <- aggregate(cbind(m1, m2) ~ day, data = cd, FUN = mean)
  ar <- data.frame(x = cen$m1[1], y = cen$m2[1], xend = cen$m1[2], yend = cen$m2[2])
  ttl <- if (modality == "RNA") bquote(italic(.(name))) else name   # genes italic, proteins roman
  ggplot(cd, aes(m1, m2)) +
    stat_ellipse(aes(fill = day, colour = day), geom = "polygon", type = "t",
                 level = 0.68, alpha = 0.15, linetype = 2, linewidth = 0.3) +
    geom_point(aes(colour = day), size = 1.2) +
    geom_segment(data = ar, aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "grey25", linewidth = 0.5, inherit.aes = FALSE,
                 arrow = arrow(type = "closed", length = unit(0.05, "inches"))) +
    annotate("text", x = mean(range(cd$m1)), y = Inf, vjust = 1.2, size = 2.4,
             label = sprintf("italic(F)*' = '*'%.2f'", F_ref), parse = TRUE) +
    scale_colour_manual(values = day_cols) +
    scale_fill_manual(values = day_cols) +
    scale_x_continuous(labels = minus) +
    scale_y_continuous(labels = minus, expand = expansion(mult = c(0.05, 0.2))) +
    labs(title = ttl, x = sprintf("MDS 1 (%.0f%%)", pct[1]), y = sprintf("MDS 2 (%.0f%%)", pct[2])) +
    base_theme
}
panels <- lapply(seq_len(nrow(sel)), function(k)
  mds_panel(sel$Feature[k], sel$Name[k], sel$Modality[k], sel$F_stat[k]))

## 4. Shared legend
leg_plot <- panels[[1]] +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 7.5),
        legend.key.height = unit(10, "pt"), legend.key.width = unit(10, "pt"))
g <- ggplotGrob(leg_plot)
idx <- grep("guide-box", g$layout$name)
idx <- idx[!vapply(g$grobs[idx], inherits, logical(1), "zeroGrob")]
legend_grob <- g$grobs[[idx[1]]]

## 5. Assemble and save
fig_h <- legend_h + n_row * panel_h
draw_s9 <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(n_row + 1, n_col,
                                             heights = unit(c(legend_h, rep(panel_h, n_row)), "in"))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:n_col)); grid.draw(legend_grob); popViewport()
  for (k in seq_along(panels)) {
    r <- (k - 1) %/% n_col + 2; cc <- (k - 1) %% n_col + 1
    pushViewport(viewport(layout.pos.row = r, layout.pos.col = cc))
    print(panels[[k]], newpage = FALSE)
    popViewport()
  }
  popViewport()
}
draw_s9()
png(file.path(plot_dir, "figS4.png"), width = page_w, height = fig_h, units = "in", res = 600)
draw_s9(); dev.off()
pdf(file.path(plot_dir, "figS4.pdf"), width = page_w, height = fig_h)
draw_s9(); dev.off()
cat("written to", plot_dir, "| rule:", rule, "\n")
