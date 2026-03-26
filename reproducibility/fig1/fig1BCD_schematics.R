# Figure 1B-D: schematic panels (differential connectivity, subgroup detection, outcome prediction)
# Output:   figures/fig1BCD.png

source("config.R")

library(ggplot2)
library(ggforce)
library(patchwork)
library(grid)

set.seed(42)

# ---- Color palette ----
blue    <- "#5B8BD4"
orange  <- "#E8833A"
red_cx  <- "#D94F4F"
purple  <- "#7B68AE"
gold    <- "#E6AB02"
gray_pt <- "#999999"
gray_tx <- "#888888"
dark_tx <- "#444444"
mid_tx  <- "#666666"
light_tx <- "#999999"
arrow_c <- "#AAAAAA"

# ---- Panel B: Differential Connectivity ----
# feature X at center, condition 1 partners left, condition 2 partners right

# Node positions
nodes_b <- data.frame(
  x     = c(0,    -1.4, -1.5, -1.4,   1.4,  1.5,  1.4),
  y     = c(0,     1.0,  0.0, -1.0,   1.0,  0.0, -1.0),
  label = c("Feature\nX", "Feature\nA", "Feature\nB", "Feature\nC", "Feature\nD", "Feature\nE", "Feature\nF"),
  color = c(red_cx, rep(blue, 3), rep(orange, 3)),
  size  = c(14, rep(14, 6)),
  stringsAsFactors = FALSE
)

# Edges from Gene X to partners
edges_b <- data.frame(
  x    = rep(0, 6),
  y    = rep(0, 6),
  xend = c(-1.4, -1.5, -1.4, 1.4, 1.5, 1.4),
  yend = c(1.0, 0.0, -1.0, 1.0, 0.0, -1.0),
  color = c(rep(blue, 3), rep(orange, 3)),
  stringsAsFactors = FALSE
)

# Expression bar data
bars_b <- data.frame(
  x     = c(-0.12, 0.12),
  y     = c(-1.75, -1.75),
  fill  = c(blue, orange)
)

panel_b <- ggplot() +
  # Edges
  geom_segment(data = edges_b,
               aes(x = x, y = y, xend = xend, yend = yend),
               color = edges_b$color, linewidth = 1.2) +
  # Partner nodes
  geom_point(data = nodes_b[-1, ],
             aes(x = x, y = y), color = "white",
             size = nodes_b$size[-1] + 4, shape = 16) +
  geom_point(data = nodes_b[-1, ],
             aes(x = x, y = y), color = nodes_b$color[-1],
             size = nodes_b$size[-1], shape = 16) +
  geom_text(data = nodes_b[-1, ],
            aes(x = x, y = y, label = label),
            color = "white", size = 2.5, fontface = "bold") +
  # Gene X (larger)
  geom_point(data = nodes_b[1, ],
             aes(x = x, y = y), color = "white", size = 19, shape = 16) +
  geom_point(data = nodes_b[1, ],
             aes(x = x, y = y), color = red_cx, size = 19, shape = 16) +
  geom_text(data = nodes_b[1, ],
            aes(x = x, y = y, label = label),
            color = "white", size = 3, fontface = "bold") +
  # Condition labels
  annotate("text", x = -1.45, y = 1.7, label = "Condition 1",
           color = blue, size = 3.5, fontface = "bold") +
  annotate("text", x = 1.45, y = 1.7, label = "Condition 2",
           color = orange, size = 3.5, fontface = "bold") +
  # Expression bars
  annotate("rect", xmin = -0.35, xmax = -0.05, ymin = -1.75, ymax = -1.15,
           fill = blue, color = NA) +
  annotate("rect", xmin = 0.05, xmax = 0.35, ymin = -1.75, ymax = -1.14,
           fill = orange, color = NA) +
  # n.s. bracket
  annotate("segment", x = -0.20, xend = -0.20, y = -1.11, yend = -1.06,
           color = light_tx, linewidth = 0.4) +
  annotate("segment", x = -0.20, xend = 0.20, y = -1.06, yend = -1.06,
           color = light_tx, linewidth = 0.4) +
  annotate("segment", x = 0.20, xend = 0.20, y = -1.06, yend = -1.11,
           color = light_tx, linewidth = 0.4) +
  annotate("text", x = 0, y = -0.95, label = "n.s.",
           color = light_tx, size = 2.8, fontface = "bold") +
  # Bottom annotation
  annotate("text", x = 0, y = -1.95, label = "Abundance unchanged",
           color = gray_tx, size = 2.8, fontface = "italic") +
  # Title
  annotate("text", x = 0, y = 2.15, label = "Differential Connectivity",
           color = dark_tx, size = 4, fontface = "bold") +
  coord_fixed(xlim = c(-2, 2), ylim = c(-2.1, 2.3)) +
  theme_void()

# ---- Panel C: Sample Subgroup Detection ----
n_pts <- 24

# Left: all mixed (gray)
gx <- rnorm(n_pts, -1.1, 0.45)
gy <- rnorm(n_pts, 0.4, 0.45)
df_mixed <- data.frame(x = gx, y = gy)

# Right: two separated groups
n_g <- 12
g1x <- rnorm(n_g, 1.05, 0.25)
g1y <- rnorm(n_g, 0.75, 0.25)
g2x <- rnorm(n_g, 1.05, 0.30)
g2y <- rnorm(n_g, -0.15, 0.30)

df_g1 <- data.frame(x = g1x, y = g1y, group = "Group 1")
df_g2 <- data.frame(x = g2x, y = g2y, group = "Group 2")
df_sep <- rbind(df_g1, df_g2)

panel_c <- ggplot() +
  # Left: mixed points
  geom_point(data = df_mixed, aes(x = x, y = y),
             color = gray_pt, size = 2.5, alpha = 0.6) +
  annotate("text", x = -1.1, y = -0.35, label = "?",
           color = "#AAAAAA", size = 10, fontface = "bold") +
  # Right: separated groups
  geom_point(data = df_g1, aes(x = x, y = y),
             color = purple, size = 2.5, alpha = 0.8) +
  geom_point(data = df_g2, aes(x = x, y = y),
             color = gold, size = 2.5, alpha = 0.8) +
  # Dashed ellipses
  geom_ellipse(aes(x0 = 1.05, y0 = 0.75, a = 0.50, b = 0.45, angle = 0),
               color = purple, linetype = "dashed", linewidth = 0.7) +
  geom_ellipse(aes(x0 = 1.05, y0 = -0.15, a = 0.55, b = 0.50, angle = 0),
               color = gold, linetype = "dashed", linewidth = 0.7) +
  # Group labels
  annotate("text", x = 1.915, y = 0.75, label = "Group 1",
           color = purple, size = 2.8, fontface = "bold") +
  annotate("text", x = 1.915, y = -0.15, label = "Group 2",
           color = gold, size = 2.8, fontface = "bold") +
  # Arrow
  annotate("segment", x = -0.35, xend = 0.30, y = 0.4, yend = 0.4,
           arrow = arrow(length = unit(0.08, "inches"), type = "closed"),
           color = arrow_c, linewidth = 0.8) +
  # Column labels
  annotate("text", x = -1.1, y = 1.5, label = "All Features",
           color = light_tx, size = 3.2, fontface = "bold") +
  annotate("text", x = 1.1, y = 1.5, label = "Selected Module",
           color = mid_tx, size = 3.2, fontface = "bold") +
  # Mini heatmap, left (all gray)
  annotate("rect", xmin = -1.45, xmax = -1.33, ymin = -1.15, ymax = -1.03, fill = "#C0C0C0", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = -1.31, xmax = -1.19, ymin = -1.15, ymax = -1.03, fill = "#C0C0C0", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = -1.17, xmax = -1.05, ymin = -1.15, ymax = -1.03, fill = "#C0C0C0", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = -0.93, xmax = -0.81, ymin = -1.15, ymax = -1.03, fill = "#C0C0C0", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = -1.45, xmax = -1.33, ymin = -1.27, ymax = -1.15, fill = "#C0C0C0", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = -1.31, xmax = -1.19, ymin = -1.27, ymax = -1.15, fill = "#C0C0C0", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = -1.17, xmax = -1.05, ymin = -1.27, ymax = -1.15, fill = "#C0C0C0", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = -0.93, xmax = -0.81, ymin = -1.27, ymax = -1.15, fill = "#C0C0C0", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = -1.45, xmax = -1.33, ymin = -1.39, ymax = -1.27, fill = "#C0C0C0", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = -1.31, xmax = -1.19, ymin = -1.39, ymax = -1.27, fill = "#C0C0C0", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = -1.17, xmax = -1.05, ymin = -1.39, ymax = -1.27, fill = "#C0C0C0", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = -0.93, xmax = -0.81, ymin = -1.39, ymax = -1.27, fill = "#C0C0C0", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = -1.45, xmax = -1.33, ymin = -1.51, ymax = -1.39, fill = "#C0C0C0", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = -1.31, xmax = -1.19, ymin = -1.51, ymax = -1.39, fill = "#C0C0C0", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = -1.17, xmax = -1.05, ymin = -1.51, ymax = -1.39, fill = "#C0C0C0", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = -0.93, xmax = -0.81, ymin = -1.51, ymax = -1.39, fill = "#C0C0C0", color = "white", linewidth = 0.3) +
  annotate("text", x = -1.13, y = -1.63, label = "Features \u00d7 Samples",
           color = gray_tx, size = 2.2) +
  # Mini heatmap, right (top 2 rows colored by group)
  annotate("rect", xmin = 0.75, xmax = 0.87, ymin = -1.15, ymax = -1.03, fill = purple, color = "white", linewidth = 0.3) +
  annotate("rect", xmin = 0.89, xmax = 1.01, ymin = -1.15, ymax = -1.03, fill = purple, color = "white", linewidth = 0.3) +
  annotate("rect", xmin = 1.19, xmax = 1.31, ymin = -1.15, ymax = -1.03, fill = gold, color = "white", linewidth = 0.3) +
  annotate("rect", xmin = 1.33, xmax = 1.45, ymin = -1.15, ymax = -1.03, fill = gold, color = "white", linewidth = 0.3) +
  annotate("rect", xmin = 0.75, xmax = 0.87, ymin = -1.27, ymax = -1.15, fill = purple, color = "white", linewidth = 0.3) +
  annotate("rect", xmin = 0.89, xmax = 1.01, ymin = -1.27, ymax = -1.15, fill = purple, color = "white", linewidth = 0.3) +
  annotate("rect", xmin = 1.19, xmax = 1.31, ymin = -1.27, ymax = -1.15, fill = gold, color = "white", linewidth = 0.3) +
  annotate("rect", xmin = 1.33, xmax = 1.45, ymin = -1.27, ymax = -1.15, fill = gold, color = "white", linewidth = 0.3) +
  annotate("rect", xmin = 0.75, xmax = 0.87, ymin = -1.39, ymax = -1.27, fill = "#DDDDDD", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = 0.89, xmax = 1.01, ymin = -1.39, ymax = -1.27, fill = "#DDDDDD", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = 1.19, xmax = 1.31, ymin = -1.39, ymax = -1.27, fill = "#DDDDDD", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = 1.33, xmax = 1.45, ymin = -1.39, ymax = -1.27, fill = "#DDDDDD", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = 0.75, xmax = 0.87, ymin = -1.51, ymax = -1.39, fill = "#DDDDDD", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = 0.89, xmax = 1.01, ymin = -1.51, ymax = -1.39, fill = "#DDDDDD", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = 1.19, xmax = 1.31, ymin = -1.51, ymax = -1.39, fill = "#DDDDDD", color = "white", linewidth = 0.3) +
  annotate("rect", xmin = 1.33, xmax = 1.45, ymin = -1.51, ymax = -1.39, fill = "#DDDDDD", color = "white", linewidth = 0.3) +
  annotate("text", x = 1.1, y = -1.63, label = "Module \u00d7 Samples",
           color = mid_tx, size = 2.2) +
  # Title
  annotate("text", x = 0, y = 2.15, label = "Sample Subgroup Detection",
           color = dark_tx, size = 4, fontface = "bold") +
  coord_fixed(xlim = c(-2, 2), ylim = c(-1.8, 2.3)) +
  theme_void()

# ---- Panel D: Clinical Outcome Prediction ----
n_each <- 14

# Left: 1D strip, overlapping
mod_x1 <- rnorm(n_each, -1.25, 0.30)
mod_y1 <- rnorm(n_each, 0.5,  0.08)
sev_x1 <- rnorm(n_each, -0.95, 0.30)
sev_y1 <- rnorm(n_each, 0.5,  0.08)

df_left_mod <- data.frame(x = mod_x1, y = mod_y1, group = "Moderate")
df_left_sev <- data.frame(x = sev_x1, y = sev_y1, group = "Severe")

# Right: 2D, overlap on x, separate on y (connectivity axis)
mod_x2 <- mod_x1 + 2.3
mod_y2 <- rnorm(n_each, -0.2, 0.15)
sev_x2 <- sev_x1 + 2.3
sev_y2 <- rnorm(n_each, 0.8,  0.15)

df_right_mod <- data.frame(x = mod_x2, y = mod_y2, group = "Moderate")
df_right_sev <- data.frame(x = sev_x2, y = sev_y2, group = "Severe")

shift_y <- 0.20

panel_d <- ggplot() +
  # Left: overlapping 1D strip
  geom_point(
    data = df_left_mod,
    aes(x = x, y = y - shift_y),
    color = blue, size = 2.5, alpha = 0.6
  ) +
  geom_point(
    data = df_left_sev,
    aes(x = x, y = y - shift_y),
    color = orange, size = 2.5, alpha = 0.6
  ) +
  
  # Left 1D axis
  annotate(
    "segment",
    x = -1.95, xend = -0.25,
    y = 0.15 - shift_y, yend = 0.15 - shift_y,
    arrow = arrow(length = unit(0.06, "inches"), type = "closed"),
    color = arrow_c, linewidth = 0.6
  ) +
  annotate(
    "text",
    x = -1.1, y = -0.05 - shift_y,
    label = "Abundance",
    color = gray_tx, size = 2.5
  ) +
  
  # Right: 2D scatter
  geom_point(
    data = df_right_mod,
    aes(x = x, y = y - shift_y),
    color = blue, size = 2.5, alpha = 0.7
  ) +
  geom_point(
    data = df_right_sev,
    aes(x = x, y = y - shift_y),
    color = orange, size = 2.5, alpha = 0.7
  ) +
  
  # Dashed ellipses
  geom_ellipse(
    aes(x0 = 1.05, y0 = -0.2 - shift_y,
        a = 0.50, b = 0.35, angle = 0),
    color = blue, linetype = "dashed", linewidth = 0.7
  ) +
  geom_ellipse(
    aes(x0 = 1.35, y0 = 0.8 - shift_y,
        a = 0.50, b = 0.35, angle = 0),
    color = orange, linetype = "dashed", linewidth = 0.7
  ) +
  
  # Group labels (right side)
  annotate(
    "text",
    x = 1.9, y = -0.2 - shift_y,
    label = "Moderate",
    color = blue, size = 2.8, fontface = "bold"
  ) +
  annotate(
    "text",
    x = 2.15, y = 0.8 - shift_y,
    label = "Severe",
    color = orange, size = 2.8, fontface = "bold"
  ) +
  
  # 2D axes
  annotate(
    "segment",
    x = 0.35, xend = 2.05,
    y = -0.65 - shift_y, yend = -0.65 - shift_y,
    arrow = arrow(length = unit(0.06, "inches"), type = "closed"),
    color = arrow_c, linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 1.2, y = -0.82 - shift_y,
    label = "Abundance",
    color = arrow_c, size = 2.5
  ) +
  annotate(
    "segment",
    x = 0.35, xend = 0.35,
    y = -0.65 - shift_y, yend = 1.2 - shift_y,
    arrow = arrow(length = unit(0.06, "inches"), type = "closed"),
    color = arrow_c, linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 0.15, y = 0.55 - shift_y,
    label = "Connectivity",
    color = arrow_c, size = 2.5, angle = 90
  ) +
  
  # Transition arrow
  annotate(
    "segment",
    x = -0.2, xend = 0.1,
    y = 0.5 - shift_y, yend = 0.5 - shift_y,
    arrow = arrow(length = unit(0.08, "inches"), type = "closed"),
    color = arrow_c, linewidth = 0.8
  ) +
  
  # Column labels
  annotate(
    "text",
    x = -1.1, y = 1.45 - shift_y,
    label = "Abundance\nOnly",
    color = light_tx, size = 3.2,
    fontface = "bold", lineheight = 0.9
  ) +
  annotate(
    "text",
    x = 1.1, y = 1.45 - shift_y,
    label = "Abundance +\nConnectivity",
    color = mid_tx, size = 3.2,
    fontface = "bold", lineheight = 0.9
  ) +
  
  # Title
  annotate(
    "text",
    x = 0, y = 2.15,
    label = "Clinical Outcome Prediction",
    color = dark_tx, size = 4, fontface = "bold"
  ) +
  
  coord_fixed(
    xlim = c(-2.2, 2.2),
    ylim = c(-1.5, 2.3)
  ) +
  theme_void()

# ---- Combine panels ----

combined <- panel_b + panel_c + panel_d +
  plot_layout(ncol = 3)

png(file.path(FIG_DIR, "fig1BCD.png"), width = 10.5, height = 3.5, units = "in", res = 600)
print(combined)
dev.off()

cat("All panels saved!\n")
