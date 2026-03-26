# Figure 4A-D: panels and assembled figure
# Input:    results/fig4/dc_corrected_paired.rds, seurat_cd4naive.rds, mosaic_embedding.rds, coupling_rows_STAT5B.rds
# Output:   figures/fig4.png, figures/fig4{A,B,C,D}.png, results/fig4/paired_limma_rna_adt.rds
# Requires: fig4/fig4_dc_paired_test.R

source("config.R")

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(limma); library(dplyr)
  library(ggplot2); library(ggrepel); library(igraph); library(grid); library(png)
})

res_dir  <- file.path(RESULTS_DIR, "fig4")
plot_dir <- FIG_DIR
dc <- readRDS(file.path(res_dir, "dc_corrected_paired.rds"))

obj <- readRDS(file.path(res_dir, "seurat_cd4naive.rds"))
seurat_rna <- obj$seurat_rna
seurat_adt <- obj$seurat_adt
rm(obj)
mosaic <- readRDS(file.path(res_dir, "mosaic_embedding.rds"))      # MOSAIC embedding
coupling_rows <- readRDS(file.path(res_dir, "coupling_rows_STAT5B.rds"))

## Shared plotting settings
target_feat <- "RNA_STAT5B"
target_gene <- "STAT5B"
day_cols <- c("Day 0" = "#56B4E9", "Day 7" = "#E69F00")
minus <- function(x) sub("^-", "\u2212", x)
fmt_p <- function(p) if (p < 0.001) formatC(p, format = "e", digits = 1) else sprintf("%.2f", p)
base_theme <- theme_classic(base_size = 9) +
  theme(axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        legend.background = element_blank(), legend.key = element_blank(),
        plot.margin = margin(4, 6, 4, 4))


## =============================================================================
## Panel A: DC pseudo-F versus paired abundance change (limma moderated t, ~ donor + time)
## =============================================================================

meta <- data.frame(sample = dc$sample_ids,
                   donor  = factor(dc$donor),
                   time   = factor(dc$cond, levels = c("0", "7")),
                   stringsAsFactors = FALSE)
rownames(meta) <- meta$sample
stopifnot(all(table(meta$donor, meta$time) == 1L))
print(table(meta$time)); print(table(meta$donor))
design <- model.matrix(~ donor + time, data = meta)

## Donor-blocked limma fit with constant and mean-trended prior
fit_paired <- function(mat) {
  stopifnot(identical(rownames(design), colnames(mat)))   # lmFit matches by order
  list(const = topTable(eBayes(lmFit(mat, design)),               coef = "time7", number = Inf, sort.by = "none"),
       trend = topTable(eBayes(lmFit(mat, design), trend = TRUE), coef = "time7", number = Inf, sort.by = "none"))
}

## 1. RNA: sample-level aggregate of SCT counts, log2-CPM
rna_feat  <- dc$feature_names[grepl("^RNA_", dc$feature_names)]
rna_genes <- sub("^RNA_", "", rna_feat)
cnt <- GetAssayData(seurat_rna, assay = "SCT", slot = "counts")
sid <- as.character(seurat_rna$sample_id)
stopifnot(all(rna_genes %in% rownames(cnt)), all(dc$sample_ids %in% unique(sid)))
cnt <- cnt[rna_genes, , drop = FALSE]
pb <- as.matrix(sapply(dc$sample_ids, function(s) Matrix::rowSums(cnt[, sid == s, drop = FALSE])))
stopifnot(identical(colnames(pb), dc$sample_ids), identical(rownames(pb), rna_genes))
cpm    <- sweep(pb, 2, colSums(pb), "/") * 1e6
logcpm <- log2(cpm + 1)
rna_fit <- fit_paired(logcpm)

## 2. ADT: sample-level mean of CLR-normalized values
adt_feat  <- dc$feature_names[grepl("^ADT_", dc$feature_names)]
adt_names <- sub("^ADT_", "", adt_feat)
adt_mat <- GetAssayData(seurat_adt, assay = "ADT", slot = "data")
sid_adt <- as.character(seurat_adt$sample_id)
stopifnot(all(adt_names %in% rownames(adt_mat)))
pb_adt <- as.matrix(sapply(dc$sample_ids, function(s)
  Matrix::rowMeans(adt_mat[adt_names, sid_adt == s, drop = FALSE])))
stopifnot(identical(colnames(pb_adt), dc$sample_ids), identical(rownames(pb_adt), adt_names))
adt_fit <- fit_paired(pb_adt)

## 3. Constant versus trended prior
for (m in c("rna_fit", "adt_fit")) {
  f <- get(m)
  cat(sprintf("%s: Spearman(t const, t trend) = %.4f | BH<0.05 const = %d, trend = %d\n", m,
              cor(f$const$t, f$trend$t, method = "spearman"),
              sum(f$const$adj.P.Val < 0.05), sum(f$trend$adj.P.Val < 0.05)))
}
cat(sprintf("STAT5B t: const = %.3f | trend = %.3f\n",
            rna_fit$const["STAT5B", "t"], rna_fit$trend["STAT5B", "t"]))

de_all <- rbind(
  data.frame(Feature = rna_feat, de_logFC = rna_fit$trend$logFC, de_t = rna_fit$trend$t,
             de_p = rna_fit$trend$P.Value, de_padj = rna_fit$trend$adj.P.Val,
             de_t_const = rna_fit$const$t, stringsAsFactors = FALSE),
  data.frame(Feature = adt_feat, de_logFC = adt_fit$trend$logFC, de_t = adt_fit$trend$t,
             de_p = adt_fit$trend$P.Value, de_padj = adt_fit$trend$adj.P.Val,
             de_t_const = adt_fit$const$t, stringsAsFactors = FALSE))
stopifnot(setequal(de_all$Feature, dc$feature_names), !anyDuplicated(de_all$Feature))
saveRDS(list(table = de_all, logcpm = logcpm, adt = pb_adt, meta = meta),
        file.path(res_dir, "paired_limma_rna_adt.rds"))         # used by panel C

## 4. Join and summarise
pa <- dc$DC_table %>%
  dplyr::select(Feature, Name, Modality, F_stat, F_rank, exact_p, at_p_floor) %>%
  dplyr::left_join(de_all, by = "Feature")
stopifnot(nrow(pa) == length(dc$feature_names), !anyNA(pa$de_t))

n_tier <- sum(pa$at_p_floor)
stopifnot(n_tier == 113L)
cat("tier by modality:\n"); print(table(pa$Modality[pa$at_p_floor]))
cat("abundance BH < 0.05, all features:", sum(pa$de_padj < 0.05), "| in tier:",
    sum(pa$at_p_floor & pa$de_padj < 0.05), "\n")

rho_all <- cor(pa$F_stat, abs(pa$de_t), method = "spearman")
rho_rna <- with(pa[pa$Modality == "RNA", ], cor(F_stat, abs(de_t), method = "spearman"))
rho_adt <- with(pa[pa$Modality == "ADT", ], cor(F_stat, abs(de_t), method = "spearman"))
rho_const <- cor(pa$F_stat, abs(pa$de_t_const), method = "spearman")
cat(sprintf("Spearman(F, |t|): all = %.3f | RNA = %.3f | ADT = %.3f | all, constant prior = %.3f\n",
            rho_all, rho_rna, rho_adt, rho_const))

s5 <- pa[pa$Feature == "RNA_STAT5B", ]
stopifnot(nrow(s5) == 1L)
cat(sprintf("STAT5B: F = %.3f | rank %d | exact p = %.4f | t = %.2f | logFC = %.3f | p = %.3f\n",
            s5$F_stat, s5$F_rank, s5$exact_p, s5$de_t, s5$de_logFC, s5$de_p))

## 5. Panel

## Minimum-P tier versus all other features
pa$grp <- factor(ifelse(pa$at_p_floor, "tier", "other"), levels = c("other", "tier"))
pa <- pa[order(pa$grp), ]                         # tier drawn on top

col_tier <- "#D98C8C"
cols   <- c(other = "grey75", tier = col_tier)
sizes  <- c(other = 0.7, tier = 1.1)
alphas <- c(other = 0.4, tier = 0.9)
p_min_lab <- paste0("1/", round(1 / dc$p_floor))
tier_lab <- bquote("Minimum exact" ~ italic(P) ~ "=" ~ .(p_min_lab))
rho_lab <- sprintf("rho~'('*italic(F)*', |'*italic(t)*'|) = '*'%s'",
                   minus(sprintf("%.3f", rho_all)))

xlim_sym <- ceiling(max(abs(pa$de_t)) * 1.05 * 2) / 2
brk  <- seq(-floor(xlim_sym), floor(xlim_sym), by = if (xlim_sym > 4) 2 else 1)
ytop <- max(pa$F_stat) * 1.06                   # headroom for corner text

pA <- ggplot(pa, aes(x = de_t, y = F_stat)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey60") +
  geom_point(aes(colour = grp, size = grp, alpha = grp), shape = 16, stroke = 0) +
  ## STAT5B: tier colour, black outline, slightly larger
  geom_point(data = s5, shape = 21, fill = col_tier, colour = "black",
             size = 2.2, stroke = 0.6) +
  geom_text_repel(data = s5, aes(label = Name), fontface = "italic", size = 2.7,
                  nudge_x = 0.7, nudge_y = 0.9, segment.size = 0.3,
                  min.segment.length = 0, seed = 1) +
  annotate("text", x = xlim_sym, y = ytop, hjust = 1, vjust = 1, size = 2.7,
           label = rho_lab, parse = TRUE) +
  scale_colour_manual(values = cols, breaks = c("tier", "other"),
                      labels = list(tier_lab, "Other features"), name = NULL) +
  scale_size_manual(values = sizes, guide = "none") +
  scale_alpha_manual(values = alphas, guide = "none") +
  scale_x_continuous(limits = c(-xlim_sym, xlim_sym), breaks = brk, labels = minus) +
  scale_y_continuous(limits = c(0, ytop), expand = expansion(mult = c(0.02, 0))) +
  labs(x = expression("Paired abundance-change statistic (moderated " * italic(t) * ")"),
       y = expression("DC statistic (pseudo-" * italic(F) * ")")) +
  theme_classic(base_size = 9) +
  theme(legend.position = c(0.01, 1), legend.justification = c(0, 1),
        legend.background = element_blank(), legend.key = element_blank(),
        legend.key.height = unit(9, "pt"), legend.key.width = unit(8, "pt"),
        legend.text = element_text(size = 7.5), legend.text.align = 0,
        legend.margin = margin(0, 0, 0, 0),
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        plot.margin = margin(4, 6, 4, 4)) +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1)))
print(pA)

## Cell size of panel A in the assembled figure
for (ext in c("pdf", "png"))
  ggsave(file.path(plot_dir, paste0("fig4A.", ext)), pA,
         width = 3.7, height = 2.55, dpi = 600)
saveRDS(pA, file.path(res_dir, "fig4A_ggplot.rds"))
cat("written to", plot_dir, "\n")


## =============================================================================
## Panel B: classical MDS of the 16 STAT5B sample embeddings
## =============================================================================

## 1. STAT5B sample embeddings
emb <- mosaic$projected_combed_list[dc$sample_ids]
stopifnot(length(emb) == 16L, identical(names(emb), dc$sample_ids))
for (i in seq_along(emb)) rownames(emb[[i]]) <- dc$feature_names
X <- t(sapply(emb, function(m) m[target_feat, ]))          # 16 samples x d
rownames(X) <- dc$sample_ids
cat("embedding dimension:", ncol(X), "\n")

## Recompute pseudo-F and compare with the stored value
D2 <- as.matrix(dist(X))^2
S  <- nrow(X)
SS_T <- sum(D2) / (2 * S)
SS_W <- sum(sapply(split(seq_len(S), dc$cond), function(ix) sum(D2[ix, ix]) / (2 * length(ix))))
F_chk <- (SS_T - SS_W) / (SS_W / (S - 2))
F_ref <- dc$DC_table$F_stat[dc$DC_table$Feature == target_feat]
rank_s5 <- dc$DC_table$F_rank[dc$DC_table$Feature == target_feat]
cat(sprintf("pseudo-F recomputed = %.4f | stored = %.4f | rank %d\n", F_chk, F_ref, rank_s5))
stopifnot(abs(F_chk - F_ref) < 1e-6)

## 2. Classical MDS on the Euclidean distances
mds <- cmdscale(dist(X), k = 2, eig = TRUE)
ev  <- mds$eig[mds$eig > 0]
pct <- 100 * mds$eig[1:2] / sum(ev)
cat(sprintf("MDS variance explained: dim 1 = %.1f%%, dim 2 = %.1f%%\n", pct[1], pct[2]))

cd <- data.frame(sample = rownames(X),
                 m1 = mds$points[, 1], m2 = mds$points[, 2],
                 donor = dc$donor,
                 day = factor(ifelse(dc$cond == "0", "Day 0", "Day 7"),
                              levels = c("Day 0", "Day 7")))
## MDS signs are arbitrary: orient dim 1 so that Day 0 lies to the left
if (mean(cd$m1[cd$day == "Day 7"]) < mean(cd$m1[cd$day == "Day 0"])) cd$m1 <- -cd$m1

cen <- aggregate(cbind(m1, m2) ~ day, data = cd, FUN = mean)
arrow_df <- data.frame(x    = cen$m1[cen$day == "Day 0"], y    = cen$m2[cen$day == "Day 0"],
                       xend = cen$m1[cen$day == "Day 7"], yend = cen$m2[cen$day == "Day 7"])

b_lab <- sprintf("italic(F)*' = '*'%.2f'", F_ref)

## 3. Panel
pB <- ggplot(cd, aes(m1, m2)) +
  stat_ellipse(aes(fill = day, colour = day), geom = "polygon", type = "t",
               level = 0.68, alpha = 0.15, linetype = 2, linewidth = 0.35) +
  geom_point(aes(colour = day), size = 1.7) +
  geom_segment(data = arrow_df, aes(x = x, y = y, xend = xend, yend = yend),
               colour = "grey25", linewidth = 0.7, inherit.aes = FALSE,
               arrow = arrow(type = "closed", length = unit(0.07, "inches"))) +
  annotate("text", x = mean(range(cd$m1)), y = Inf, vjust = 1.3, size = 2.7,
           label = b_lab, parse = TRUE) +
  ## one legend for colour and fill
  scale_colour_manual(values = day_cols, name = NULL) +
  scale_fill_manual(values = day_cols, name = NULL) +
  scale_x_continuous(labels = minus) +
  scale_y_continuous(labels = minus, expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = sprintf("MDS 1 (%.0f%%)", pct[1]), y = sprintf("MDS 2 (%.0f%%)", pct[2])) +
  base_theme +
  theme(legend.position = "right",
        legend.text = element_text(size = 7.5),
        legend.key.height = unit(11, "pt"), legend.key.width = unit(11, "pt"),
        legend.margin = margin(0, 0, 0, -4),
        legend.box.margin = margin(0, 0, 0, 0))
print(pB)

for (ext in c("png", "pdf"))
  ggsave(file.path(plot_dir, paste0("fig4B.", ext)), pB,
         width = 3.3, height = 2.9, dpi = 600)
saveRDS(pB, file.path(res_dir, "fig4B_ggplot.rds"))
cat("written to", plot_dir, "\n")


## =============================================================================
## Panel C: STAT5B sample-level abundance (log2-CPM), paired by donor
## =============================================================================

lim <- readRDS(file.path(res_dir, "paired_limma_rna_adt.rds"))

## Data
stopifnot(target_gene %in% rownames(lim$logcpm))
bd <- data.frame(sample = colnames(lim$logcpm),
                 value  = as.numeric(lim$logcpm[target_gene, ]),
                 stringsAsFactors = FALSE)
bd$donor <- as.character(lim$meta[bd$sample, "donor"])
bd$day   <- factor(ifelse(lim$meta[bd$sample, "time"] == "0", "Day 0", "Day 7"),
                   levels = c("Day 0", "Day 7"))
stopifnot(!anyNA(bd$donor), all(table(bd$donor) == 2L))
bd$x <- as.numeric(bd$day)

st <- lim$table[lim$table$Feature == target_feat, ]
stopifnot(nrow(st) == 1L)

means <- aggregate(value ~ day, data = bd, FUN = mean)
means$x <- as.numeric(means$day) + c(-0.14, 0.14)     # beside the donor points
cat(sprintf("STAT5B: limma log2FC = %.3f | difference of means = %.3f | t = %.2f | P = %.3f | BH = %.3f\n",
            st$de_logFC, diff(means$value), st$de_t, st$de_p, st$de_padj))
stopifnot(abs(st$de_logFC - diff(means$value)) < 1e-8)

c_lab <- sprintf("atop('log'[2]*'FC = '*'%s', italic(P)*' = '*'%s')",
                 minus(sprintf("%.2f", st$de_logFC)), fmt_p(st$de_p))

## Panel
pC <- ggplot(bd, aes(x = x, y = value)) +
  geom_line(aes(group = donor), colour = "grey70", linewidth = 0.4) +
  geom_point(aes(colour = day), size = 1.7) +
  ## mean per time point
  geom_line(data = means, aes(x = x, y = value, group = 1), colour = "grey20", linewidth = 0.55) +
  geom_point(data = means, aes(x = x, y = value), shape = 23, fill = "white",
             colour = "grey20", size = 2.0, stroke = 0.5) +
  annotate("text", x = 1.5, y = Inf, vjust = 1.1, size = 2.7, label = c_lab, parse = TRUE) +
  scale_colour_manual(values = day_cols, guide = "none") +
  scale_x_continuous(breaks = 1:2, labels = c("Day 0", "Day 7"),
                     expand = expansion(add = 0.25)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(x = NULL, y = expression(italic("STAT5B") ~ "(log"[2] * " CPM)")) +
  base_theme
print(pC)

for (ext in c("png", "pdf"))
  ggsave(file.path(plot_dir, paste0("fig4C.", ext)), pC,
         width = 1.7, height = 2.6, dpi = 600)
saveRDS(pC, file.path(res_dir, "fig4C_ggplot.rds"))


## =============================================================================
## Panel D: top 20 features coupled to STAT5B at Day 0 and Day 7 (condition-averaged coupling)
## =============================================================================

top_n <- 20
bold_name <- character(0)         # no bold labels

day0_color <- "#56B4E9"
day7_color <- "#E69F00"
hero_color <- "red"
width_in <- 5.15; height_in <- 2.6  # panel D cell size (in)
lab_pad_x <- 0.50                 # label space left and right, in star radii
lab_pad_y <- 0.32                 # label space above and below, in star radii
lab_d     <- 0.15                 # gap between a node and its label
lab_cex   <- 0.7
node_cex  <- 2.2                  # partner node size (points, in character units)
hero_cex  <- 2.5                  # STAT5B node size

## Horizontal label nudges, in star radii
nudge_x <- c(ZDHHC2 = -0, NRAS = -0, BCAS4 = -0,
             GPR174 = -0, `TIM-4` = -0, SYNE2 = -0, `CD38-1` = -0)

## 1. Strongest couplings in the condition-averaged coupling matrices
## STAT5B row of each sample-specific coupling matrix, named by dc$feature_names
crow <- coupling_rows[dc$sample_ids]
stopifnot(length(crow) == 16L,
          all(vapply(crow, function(v) identical(names(v), dc$feature_names), logical(1))))

avg_coupling <- function(day) Reduce(`+`, crow[dc$sample_ids[dc$cond == day]]) / sum(dc$cond == day)
top_coupled <- function(v, f, K) {
  v <- v[names(v) != f]
  names(sort(v, decreasing = TRUE))[1:K]
}
top_neighbors_day0 <- top_coupled(avg_coupling("0"), target_feat, top_n)
top_neighbors_day7 <- top_coupled(avg_coupling("7"), target_feat, top_n)
cat("overlap between Day 0 and Day 7 lists:",
    length(intersect(top_neighbors_day0, top_neighbors_day7)), "of", top_n, "\n")
print(list(day0 = top_neighbors_day0, day7 = top_neighbors_day7))
saveRDS(list(day0 = top_neighbors_day0, day7 = top_neighbors_day7, K = top_n, space = "coupling"),
        file.path(res_dir, "stat5b_top_couplings.rds"))

## 2. Graphs
clean_name <- function(x) sub("^(RNA_|ADT_)", "", x)
modality   <- function(x) ifelse(grepl("^ADT_", x), "ADT", "RNA")
hero <- clean_name(target_feat)

build_graph <- function(neigh, partner_color) {
  nodes <- data.frame(name = c(hero, clean_name(neigh)),
                      modality = c("Hero", modality(neigh)),
                      stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(nodes$name))
  g <- graph_from_data_frame(data.frame(from = hero, to = clean_name(neigh)),
                             vertices = nodes, directed = FALSE)
  v   <- V(g)$name
  mod <- nodes$modality[match(v, nodes$name)]
  ## font: 1 roman, 3 italic (genes italic, proteins roman), +1 for bold
  font <- ifelse(mod == "ADT", 1, 3)
  font[v %in% bold_name] <- font[v %in% bold_name] + 1
  list(g = g,
       layout = layout_as_star(g, center = V(g)[v == hero]),
       color = ifelse(v == hero, hero_color, partner_color),
       shape = ifelse(mod == "ADT", "square", "circle"),
       font  = font)
}
g0 <- build_graph(top_neighbors_day0, day0_color)
g7 <- build_graph(top_neighbors_day7, day7_color)

draw <- function(x, title) {
  ## star layout scaled to a unit circle, drawn with asp = 1
  nm  <- V(x$g)$name
  isH <- nm == hero
  lay <- x$layout
  lay <- sweep(lay, 2, lay[isH, ])                      # centre on STAT5B
  lay <- lay / max(sqrt(rowSums(lay^2)))                # ring of radius 1
  lx <- 1 + lab_pad_x; ly <- 1 + lab_pad_y
  
  plot.new()
  plot.window(xlim = c(-lx, lx), ylim = c(-ly, ly), asp = 1, xaxs = "i", yaxs = "i")
  title(main = title, font.main = 1, cex.main = 1.0, line = 0.1)
  
  ## x and y units per inch (equal for a true circle)
  u <- par("usr"); p <- par("pin")
  cat(sprintf("%s: x units/in = %.4f | y units/in = %.4f\n",
              title, (u[2] - u[1]) / p[1], (u[4] - u[3]) / p[2]))
  
  segments(0, 0, lay[!isH, 1], lay[!isH, 2], col = "grey", lwd = 1)
  pch <- ifelse(x$shape == "square", 22, 21)
  points(lay[!isH, 1], lay[!isH, 2], pch = pch[!isH], bg = x$color[!isH],
         col = "black", cex = node_cex, lwd = 1)
  points(0, 0, pch = 21, bg = x$color[isH], col = "black", cex = hero_cex, lwd = 1)
  
  ang <- atan2(lay[, 2], lay[, 1])
  dx  <- ifelse(nm %in% names(nudge_x), nudge_x[match(nm, names(nudge_x))], 0)
  dx[is.na(dx)] <- 0
  for (k in which(!isH)) {
    cx <- cos(ang[k]); sy <- sin(ang[k])
    adj <- if (cx > 0.3) c(0, 0.5) else if (cx < -0.3) c(1, 0.5) else
      if (sy > 0) c(0.5, 0) else c(0.5, 1)
    text(lay[k, 1] + lab_d * cx + dx[k], lay[k, 2] + lab_d * sy,
         labels = nm[k], adj = adj, cex = lab_cex, font = x$font[k],
         family = "sans", xpd = NA)
  }
  text(0, -lab_d * 0.9, labels = nm[isH], adj = c(0.5, 1),
       cex = lab_cex, font = x$font[isH], family = "sans", xpd = NA)
}

## 3. Render at the size of the panel D cell
render <- function(file, dev) {
  if (dev == "png") png(file, width = width_in, height = height_in, units = "in", res = 600,
                        pointsize = 9, family = "sans")
  else pdf(file, width = width_in, height = height_in, pointsize = 9, family = "Helvetica")
  par(mfrow = c(1, 2), mar = c(0, 0.1, 1.1, 0.1), oma = c(1.4, 0, 0, 0))
  draw(g0, "Day 0")
  draw(g7, "Day 7")
  
  ## One-row legend
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlim = c(0, 1), ylim = c(0, 1))
  leg_lab <- c("Top 20 coupled features, Day 0", "Top 20 coupled features, Day 7",
               "RNA feature", "ADT feature")
  leg_cex <- 0.74
  tw <- strwidth(leg_lab, cex = leg_cex) + 0.02
  leg_args <- list(x = "bottom", legend = leg_lab,
                   pt.bg = c(day0_color, day7_color, "white", "white"),
                   pch = c(21, 21, 21, 22), col = "black",
                   pt.cex = 1.2, bty = "n", cex = leg_cex, horiz = TRUE,
                   x.intersp = 0.9, xpd = NA, inset = c(0, 0.005))
  ## fall back to a single text width if per-entry widths fail
  ok <- tryCatch({ do.call(legend, c(leg_args, list(text.width = tw))); TRUE },
                 error = function(e) FALSE)
  if (!ok) do.call(legend, c(leg_args, list(text.width = max(tw))))
  dev.off()
}
render(file.path(plot_dir, "fig4D.png"), "png")
render(file.path(plot_dir, "fig4D.pdf"), "pdf")
cat("written to", plot_dir, "\n")


## =============================================================================
## Assembled Figure 4: row 1 A, B; row 2 C, D (panel D placed from its PNG)
## =============================================================================

stopifnot(exists("pA"), exists("pB"), exists("pC"))

page_w <- 6.85
row_h  <- c(2.55, 2.6)
w_row1 <- c(3.70, 3.15)            # A, B (MDS)
w_row2 <- c(1.70, 5.15)            # C (abundance), D
stopifnot(abs(sum(w_row1) - page_w) < 1e-8, abs(sum(w_row2) - page_w) < 1e-8)
fig_h <- sum(row_h)

## Panels
tag_room <- theme(plot.margin = margin(14, 6, 4, 4))   # space for the panel tag
panel_A <- pA + tag_room
## MDS legend on top
panel_B <- pB + tag_room +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-4, 0, -6, 0))
panel_C <- pC + tag_room

d_file <- file.path(plot_dir, "fig4D.png")
stopifnot(file.exists(d_file))
imgD <- readPNG(d_file)
cat(sprintf("panel D image: %d x %d px (expected about %d x %d at 600 dpi)\n",
            dim(imgD)[2], dim(imgD)[1], round(w_row2[2] * 600), round(row_h[2] * 600)))
gD <- rasterGrob(imgD, interpolate = TRUE)

## Drawing
tag <- function(letter) {
  grid.text(letter, x = unit(0.04, "in"), y = unit(1, "npc") - unit(0.04, "in"),
            just = c("left", "top"), gp = gpar(fontsize = 11, fontface = "bold"))
}
place_gg <- function(p, letter, col) {
  pushViewport(viewport(layout.pos.col = col)); print(p, newpage = FALSE); tag(letter); popViewport()
}

draw_figure <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(row_h, "in"))))
  
  ## Row 1
  pushViewport(viewport(layout.pos.row = 1,
                        layout = grid.layout(1, 2, widths = unit(w_row1, "in"))))
  place_gg(panel_A, "A", 1)
  place_gg(panel_B, "B", 2)
  popViewport()
  
  ## Row 2
  pushViewport(viewport(layout.pos.row = 2,
                        layout = grid.layout(1, 2, widths = unit(w_row2, "in"))))
  place_gg(panel_C, "C", 1)
  pushViewport(viewport(layout.pos.col = 2)); grid.draw(gD); tag("D"); popViewport()
  popViewport()
  
  popViewport()
}

draw_figure()

png(file.path(plot_dir, "fig4.png"), width = page_w, height = fig_h, units = "in", res = 600)
draw_figure(); dev.off()
pdf(file.path(plot_dir, "fig4.pdf"), width = page_w, height = fig_h)
draw_figure(); dev.off()
tiff(file.path(plot_dir, "fig4.tiff"), width = page_w, height = fig_h, units = "in",
     res = 600, compression = "lzw")
draw_figure(); dev.off()
cat("written to", plot_dir, "\n")
