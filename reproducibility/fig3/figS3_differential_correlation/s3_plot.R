# Supplemental Figure S3: panels A-C
# Input:    stage4/formal_auprc.csv, s3_frozen.rds
# Output:   stage4/figS3.png
# Requires: s3_stage4_run.R --mode formal
suppressPackageStartupMessages({ library(ggplot2); library(grid) })

a  <- read.csv("stage4/formal_auprc.csv")
fz <- readRDS("s3_frozen.rds")
A_of_kappa <- fz$A_fine[as.character(fz$kappa_grid)]          # frozen pooled alignment
prevalence <- 2 / fz$cfg$r                                     # random baseline: DC = 2 of r equal modules

# ablation drawn dashed
abl <- "MOSAIC without projection (ablation)"
methods <- c(mosaic = "MOSAIC", dgca = "DGCA", memento = "memento", raw = abl)
cols    <- setNames(c("#1B9E77", "#D95F02", "#7570B3", "#7F7F7F"), methods)
ltys    <- setNames(c("solid", "solid", "solid", "dashed"), methods)
a <- a[a$method %in% names(methods), ]
a$Method <- factor(methods[a$method], levels = methods)

ci_summary <- function(d, by) {
  s <- do.call(rbind, lapply(split(d, d[by], drop = TRUE), function(x) {
    n <- nrow(x); m <- mean(x$auprc); h <- qt(0.975, n - 1) * sd(x$auprc) / sqrt(n)
    data.frame(x[1, by, drop = FALSE], mean = m, lo = m - h, hi = m + h, n = n)
  }))
  rownames(s) <- NULL; s
}

FONT <- "Helvetica"
theme_s3 <- theme_classic(base_size = 8, base_family = FONT) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.key.width = unit(1.4, "lines"), legend.text = element_text(size = 8),
        axis.text = element_text(size = 8, colour = "black"),
        axis.title = element_text(size = 8),
        plot.title = element_blank(),
        plot.margin = margin(t = 14, r = 4, b = 2, l = 4))

# A: AUPRC vs alignment
dA <- a[a$panel == "A", ]
dA$align <- unname(A_of_kappa[as.character(dA$kappa)])
sA <- ci_summary(dA, c("align", "kappa", "Method"))
labA <- sprintf("%.2f", A_of_kappa)
dodgeA <- position_dodge(width = 0.02); dodgeB <- position_dodge(width = 0.02)
pA <- ggplot(sA, aes(align, mean, colour = Method, linetype = Method, group = Method)) +
  geom_hline(yintercept = prevalence, linetype = "dotted", colour = "grey70", linewidth = 0.4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.015, position = dodgeA, linetype = "solid", linewidth = 0.3) +
  geom_line(position = dodgeA, linewidth = 0.4) + geom_point(size = 1.2, position = dodgeA) +
  scale_x_continuous(breaks = unname(A_of_kappa), labels = labA) +
  scale_colour_manual(values = cols) + scale_linetype_manual(values = ltys) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Alignment A of the condition change", y = "AUPRC") + theme_s3

# B: AUPRC vs strength of sample-private programs
dB <- a[a$panel == "B", ]
sB <- ci_summary(dB, c("eta", "Method"))
pB <- ggplot(sB, aes(eta, mean, colour = Method, linetype = Method, group = Method)) +
  geom_hline(yintercept = prevalence, linetype = "dotted", colour = "grey70", linewidth = 0.4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.015, position = dodgeB, linetype = "solid", linewidth = 0.3) +
  geom_line(position = dodgeB, linewidth = 0.4) + geom_point(size = 1.2, position = dodgeB) +
  scale_x_continuous(breaks = sort(unique(dB$eta))) +
  scale_colour_manual(values = cols) + scale_linetype_manual(values = ltys) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Strength of sample-private programs \u03b7", y = "AUPRC") + theme_s3

# C: paired difference MOSAIC minus unprojected profile
w <- reshape(dA[dA$method %in% c("mosaic", "raw"), c("align", "kappa", "seed", "method", "auprc")],
             idvar = c("align", "kappa", "seed"), timevar = "method", direction = "wide")
w$diff <- w$auprc.mosaic - w$auprc.raw
sC <- do.call(rbind, lapply(split(w, w$align), function(x) {
  n <- nrow(x); m <- mean(x$diff); h <- qt(0.975, n - 1) * sd(x$diff) / sqrt(n)
  data.frame(align = x$align[1], mean = m, lo = m - h, hi = m + h)
}))
set.seed(1)
pC <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_jitter(data = w, aes(align, diff), width = 0.006, height = 0,
              colour = "grey55", alpha = 0.7, size = 0.9) +
  geom_errorbar(data = sC, aes(align, ymin = lo, ymax = hi), width = 0.012, colour = "black", linewidth = 0.3) +
  geom_line(data = sC, aes(align, mean), colour = "black", linewidth = 0.4) +
  geom_point(data = sC, aes(align, mean), colour = "black", size = 1.2) +
  scale_x_continuous(breaks = unname(A_of_kappa), labels = labA) +
  labs(x = "Alignment A of the condition change",
       y = "\u0394AUPRC (MOSAIC \u2212 no projection)") + theme_s3

# Layout with grid; shared legend from panel A
get_legend_grob <- function(p) {
  g <- ggplotGrob(p + theme(legend.position = "bottom"))
  idx <- grep("guide-box", vapply(g$grobs, function(x) x$name, ""))
  idx <- idx[!vapply(g$grobs[idx], inherits, TRUE, what = "zeroGrob")]
  if (!length(idx)) stop("legend grob not found")
  g$grobs[[idx[1]]]
}
legend <- get_legend_grob(pA)
panels <- lapply(list(pA, pB, pC), function(p) ggplotGrob(p + theme(legend.position = "none")))

draw_figure <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 3, heights = unit(c(1, 0.11), "null"))))
  for (i in 1:3) {
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i))
    grid.draw(panels[[i]])
    grid.text(LETTERS[i], x = unit(2, "pt"), y = unit(1, "npc") - unit(1, "pt"),
              just = c("left", "top"), gp = gpar(fontsize = 12, fontface = "bold", fontfamily = FONT))
    popViewport()
  }
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:3)); grid.draw(legend); popViewport()
  popViewport()
}
dir.create("stage4", showWarnings = FALSE)
# cairo devices for the Greek letters
if (capabilities("cairo")) cairo_pdf("stage4/figS3.pdf", width = 7.2, height = 3.0) else
  pdf("stage4/figS3.pdf", width = 7.2, height = 3.0)
draw_figure(); invisible(dev.off())
png("stage4/figS3.png", width = 7.2, height = 3.0, units = "in", res = 300,
    type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
draw_figure(); invisible(dev.off())
# kappa value behind each alignment tick
cat("x-axis mapping (A <- kappa):", paste(sprintf("%.2f <- %s", A_of_kappa, names(A_of_kappa)), collapse = "; "), "\n")
cat("written stage4/figS3.pdf and .png\n")
