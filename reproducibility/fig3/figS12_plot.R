# Supplemental Figure S12: composition shift versus within-population analysis (figure)
# Input:    results/S12/results_rep*.rds, results/S12/run_metadata.rds
# Output:   figures/figS12.png
# Requires: fig3/figS12_composition_simulation.R
# Usage:    Rscript fig3/figS12_plot.R

source("config.R")

## Panels: A, false-DC rate vs composition severity; B, QQ at the strongest severity; C, Arm P AUPRC

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

args     <- commandArgs(trailingOnly = TRUE)
want_pdf <- any(args == "--pdf")
args     <- args[args != "--pdf"]
in_dir   <- if (length(args) >= 1) args[1] else file.path(RESULTS_DIR, "S12")
out_dir  <- if (length(args) >= 2) args[2] else FIG_DIR
DPI      <- 600
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## PNG device at the requested resolution, falling back if cairo is unavailable
open_png <- function(path, width, height) {
  ok <- tryCatch({
    png(path, width = width, height = height, units = "in", res = DPI, type = "cairo")
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) png(path, width = width, height = height, units = "in", res = DPI)
}

## Write one ggplot object as PNG (and optionally PDF) at a physical size
save_plot <- function(p, stem, width, height) {
  open_png(file.path(out_dir, paste0(stem, ".png")), width, height)
  print(p); dev.off()
  if (want_pdf)
    ggsave(file.path(out_dir, paste0(stem, ".pdf")), p,
           width = width, height = height)
  invisible(NULL)
}

## Appearance
COL <- c("Mixture-level" = "#C0435A", "Within-population" = "#0F7A8C")
BASE_SIZE <- 9
theme_s12 <- function() {
  theme_minimal(base_size = BASE_SIZE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
          axis.title       = element_text(size = BASE_SIZE),
          axis.text        = element_text(size = BASE_SIZE - 1, colour = "grey25"),
          legend.position  = "top",
          legend.margin    = margin(b = -4),
          legend.key.size  = unit(9, "pt"),
          plot.caption     = element_text(size = BASE_SIZE - 2.5, hjust = 0,
                                          colour = "grey30"),
          plot.title       = element_text(size = BASE_SIZE, face = "plain"),
          plot.tag         = element_text(size = BASE_SIZE + 2, face = "bold"),
          plot.tag.position = c(0, 1))
}

## Load
stopifnot(dir.exists(in_dir))
files <- list.files(in_dir, "^results_rep[0-9]+\\.rds$", full.names = TRUE)
if (!length(files)) stop("no results_rep*.rds found in ", in_dir)
reps <- lapply(files, readRDS)
meta <- if (file.exists(file.path(in_dir, "run_metadata.rds")))
  readRDS(file.path(in_dir, "run_metadata.rds")) else NULL
cfg  <- if (!is.null(meta)) meta$cfg else NULL

summ <- do.call(rbind, lapply(reps, function(r)
  do.call(rbind, lapply(names(r$arms), function(nm) {
    a <- r$arms[[nm]]
    data.frame(replicate = r$replicate, arm = substr(nm, 1, 1),
               severity = a$severity, fpr05 = a$fpr05, fpr01 = a$fpr01,
               n_bh_hits = a$n_bh_hits,
               auprc = if (is.null(a$auprc)) NA_real_ else a$auprc,
               stringsAsFactors = FALSE)
  }))))

## severity order and the abundance each label corresponds to
sev_levels <- if (!is.null(cfg)) names(cfg$severities) else unique(summ$severity)
sev_pB <- if (!is.null(cfg)) unlist(cfg$severities) else
  setNames(seq_along(sev_levels), sev_levels)
summ$severity <- factor(summ$severity, levels = sev_levels)
summ$pB       <- sev_pB[as.character(summ$severity)]
n_rep    <- length(reps)
alpha    <- if (!is.null(cfg)) cfg$alpha_fpr else 0.05
pA_focal <- if (!is.null(cfg)) cfg$pA_focal  else 0.10
dc_prev  <- if (!is.null(cfg)) cfg$dc_prop   else 0.20   # random-ranking AUPRC
lab    <- c(M = "Mixture-level", L = "Within-population")

## Panel A: false-DC rate vs severity
agg <- do.call(rbind, lapply(split(summ[summ$arm %in% c("M", "L"), ],
                                   ~ arm + severity, drop = TRUE), function(s) {
  n <- nrow(s); se <- sd(s$fpr05) / sqrt(n)
  tq <- if (n > 1L) qt(0.975, df = n - 1L) else NA_real_
  data.frame(arm = s$arm[1], severity = s$severity[1], pB = s$pB[1],
             mean_fpr05 = mean(s$fpr05),
             lo = mean(s$fpr05) - tq * se, hi = mean(s$fpr05) + tq * se)
}))
agg$Analysis <- unname(lab[agg$arm])
pts <- summ[summ$arm %in% c("M", "L"), ]; pts$Analysis <- unname(lab[pts$arm])

## y ticks: pretty breaks plus the nominal level
y_rng    <- range(c(pts$fpr05, agg$lo, agg$hi, alpha), na.rm = TRUE)
y_breaks <- pretty(y_rng, n = 5)
y_breaks <- sort(c(alpha, y_breaks[abs(y_breaks - alpha) > 0.3 * diff(y_rng) / 5]))
y_breaks <- y_breaks[y_breaks >= y_rng[1] - 1e-9]   # no tick below the plotted range

mk_A <- function(with_caption, tag = NULL) ggplot(agg, aes(pB, mean_fpr05, colour = Analysis, group = Analysis)) +
  geom_point(data = pts, aes(y = fpr05),
             position = position_jitter(width = .004, height = 0),
             size = .7, alpha = .28, stroke = 0) +
  geom_line(linewidth = .55) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .012, linewidth = .5) +
  geom_point(size = 1.7) +
  ## nominal level
  geom_hline(yintercept = alpha, linetype = "22", colour = "grey15", linewidth = .35) +
  scale_colour_manual(values = COL) +
  scale_x_continuous(breaks = sev_pB, labels = sprintf("%.2f", sev_pB)) +
  scale_y_continuous(breaks = y_breaks, labels = function(x) sprintf("%.2f", x)) +
  labs(x = sprintf("Focal-population proportion in condition B\n(condition A = %.2f)", pA_focal),
       y = bquote(Fraction~of~features~with~italic(p) < .(alpha)), colour = NULL,
       caption = if (with_caption) paste0(
         "No rewiring: conditional distributions identical across conditions.\n",
         "Dashed line, nominal level. Points, individual replicates;\nbars, 95% CI across ",
         n_rep, " replicates.") else NULL, tag = tag) +
  theme_s12()
save_plot(mk_A(TRUE), "panelA_fpr_ladder", 3.4, 3.2)

## Panel B: QQ at the strongest severity
## keep the significant tail, thin the bulk
qq_thin <- function(p, keep_tail = 500L, bulk_max = 1500L) {
  p <- sort(p); n <- length(p)
  e <- -log10(ppoints(n)); o <- -log10(p)
  idx <- seq_len(n)
  tail_idx <- idx[seq_len(min(keep_tail, n))]
  bulk_idx <- setdiff(idx, tail_idx)
  if (length(bulk_idx) > bulk_max)
    bulk_idx <- bulk_idx[round(seq(1, length(bulk_idx), length.out = bulk_max))]
  k <- sort(c(tail_idx, bulk_idx))
  data.frame(expected = e[k], observed = o[k])
}
sv_qq <- sev_levels[which.max(sev_pB)]
qq <- do.call(rbind, lapply(reps, function(r)
  do.call(rbind, lapply(names(r$arms), function(nm) {
    a <- r$arms[[nm]]
    if (substr(nm, 1, 1) == "P" || a$severity != sv_qq) return(NULL)
    d <- qq_thin(a$p)
    d$Analysis <- unname(lab[substr(nm, 1, 1)]); d$replicate <- r$replicate; d
  }))))

mk_B <- function(with_caption, tag = NULL, show_legend = TRUE, top_pad = 0, bottom_pad = 5.5) ggplot(qq, aes(expected, observed, colour = Analysis,
                          group = interaction(Analysis, replicate))) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45", linewidth = .3) +
  geom_line(linewidth = .3, alpha = .5) +
  scale_colour_manual(values = COL) +
  labs(x = expression(Expected~-log[10](italic(p))),
       y = expression(Observed~-log[10](italic(p))), colour = NULL,
       caption = if (with_caption)
         sprintf("Focal-population abundance %.2f (%s);\nexact permutation null. One line per\nreplicate (%d).",
                 sev_pB[[sv_qq]], sv_qq, n_rep) else NULL, tag = tag) +
  theme_s12() +
  theme(legend.position = if (show_legend) "top" else "none",
        plot.margin = margin(t = top_pad, r = 5.5, b = bottom_pad, l = 5.5, unit = "pt"))
save_plot(mk_B(TRUE), "panelB_qq_strong", 3.4, 3.2)

## Panel C: Arm P AUPRC
au <- summ[summ$arm == "P" & !is.na(summ$auprc), ]
mk_C <- function(with_caption, tag = NULL, top_pad = 0, bottom_pad = 5.5) ggplot(au, aes(x = "", y = auprc)) +
  geom_boxplot(width = .32, fill = "grey93", colour = "grey35",
               linewidth = .35, outlier.shape = NA) +
  geom_jitter(width = .06, height = 0, size = 1.1, alpha = .6,
              colour = COL[["Within-population"]], stroke = 0) +
  coord_cartesian(ylim = c(0.5, 1)) +
  scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
  labs(x = NULL, y = "AUPRC vs DC ground truth",
       caption = if (with_caption) paste0("Composition shift plus genuine\n",
                        "within-population rewiring,\nsame workflow.") else NULL,
       tag = tag) +
  theme_s12() +
  theme(axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(),
        plot.margin = margin(t = top_pad, r = 5.5, b = bottom_pad, l = 5.5, unit = "pt"))
save_plot(mk_C(TRUE), "panelC_auprc_arm_P", 2.1, 3.2)

## Assembled figure (grid, 1 x 3 layout)
FIG_W <- 8.0; FIG_H <- 3.0
LEGEND_PAD <- 22   # pt, legend row height in panel A
XLAB_LINE  <- 11   # pt, one line of x-axis title
draw_combined <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    1, 3, widths = unit(c(1, 1, 0.58), "null"))))
  vp <- function(col) viewport(layout.pos.row = 1, layout.pos.col = col)
  ## panel tags drawn separately at a fixed height
  print(mk_A(FALSE), vp = vp(1))
  print(mk_B(FALSE, show_legend = FALSE, top_pad = LEGEND_PAD,
             bottom_pad = 5.5 + XLAB_LINE), vp = vp(2))          # A has 2 lines, B 1
  print(mk_C(FALSE, top_pad = LEGEND_PAD,
             bottom_pad = 5.5 + 2 * XLAB_LINE + 9), vp = vp(3))  # C has none, no tick row
  for (k in 1:3) {
    pushViewport(vp(k))
    grid.text(LETTERS[k], x = unit(3, "pt"), y = unit(1, "npc") - unit(3, "pt"),
              just = c("left", "top"),
              gp = gpar(fontface = "bold", fontsize = BASE_SIZE + 3))
    popViewport()
  }
  popViewport()
}
open_png(file.path(out_dir, "figS12.png"), FIG_W, FIG_H)
draw_combined(); dev.off()
if (want_pdf) {
  pdf(file.path(out_dir, "figS12.pdf"),
      width = FIG_W, height = FIG_H)
  draw_combined(); dev.off()
}

## Numbers for the figure legend
num <- function(x, d = 4) formatC(x, format = "f", digits = d)
cat("\n================ numbers for the figure legend ============\n")
cat(sprintf("replicates: %d   features: %s   samples: %s (balanced)\n", n_rep,
            if (!is.null(cfg)) cfg$px + cfg$py + cfg$pz else "?",
            if (!is.null(cfg)) cfg$n_samples else "?"))
if (!is.null(cfg))
  cat(sprintf("exact null: %d nonredundant balanced label assignments enumerated (one per complementary pair)\n",
              choose(cfg$n_samples - 1L, cfg$n_samples / 2L - 1L)))
cat("\nFalse-DC rate at nominal", alpha, "(mean [95% CI across replicates]):\n")
for (sv in sev_levels) for (a in c("M", "L")) {
  r <- agg[agg$arm == a & agg$severity == sv, ]
  cat(sprintf("  %-18s %-9s (pB=%.2f): %s [%s, %s]\n", lab[[a]], sv, r$pB,
              num(r$mean_fpr05), num(r$lo), num(r$hi)))
}
cat("\nWithin-population arm vs nominal (one-sample t on replicate means):\n")
safe_t <- function(x, mu) {
  if (length(x) < 2L || diff(range(x)) < .Machine$double.eps) return(NA_real_)
  t.test(x, mu = mu)$p.value
}
Lm <- tapply(summ$fpr05[summ$arm == "L"], summ$replicate[summ$arm == "L"], mean)
L1 <- tapply(summ$fpr01[summ$arm == "L"], summ$replicate[summ$arm == "L"], mean)
cat(sprintf("  at 0.05: mean %s, p = %s\n", num(mean(Lm)),
            formatC(safe_t(Lm, 0.05), format = "f", digits = 3)))
cat(sprintf("  at 0.01: mean %s, p = %s\n", num(mean(L1), 5),
            formatC(safe_t(L1, 0.01), format = "f", digits = 3)))
cat("\nBH-significant false discoveries (mixture-level arm):\n")
for (sv in sev_levels) {
  x <- summ$n_bh_hits[summ$arm == "M" & summ$severity == sv]
  cat(sprintf("  %-9s: %d/%d replicates > 0, median %.0f, max %d\n",
              sv, sum(x > 0), length(x), median(x), max(x)))
}
xl <- summ$n_bh_hits[summ$arm == "L"]
cat(sprintf("  within-population arm, all severities: %d/%d runs > 0, max %d\n",
            sum(xl > 0), length(xl), max(xl)))
cat(sprintf("\nArm P AUPRC: mean %s (sd %s, range %s-%s)\n",
            num(mean(au$auprc)), num(sd(au$auprc)),
            num(min(au$auprc)), num(max(au$auprc))))

cat(sprintf("\nWrote panels and figS12 to %s (PNG at %d dpi%s)\n",
            normalizePath(out_dir), DPI, if (want_pdf) ", plus PDF" else ""))
