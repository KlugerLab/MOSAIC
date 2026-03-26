# Figure 6D, 6E: Lasso-selected predictive features of the abundance and connectivity (MOSAIC) monocyte models
# Input:    results/fig6/fig6BC_objects.rds
# Output:   figures/fig6D.png, figures/fig6E.png
# Requires: fig6/fig6BC_prediction_benchmark.R

source("config.R")

library(glmnet)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggVennDiagram)
library(grid)
library(cowplot)      # plot_grid

fig6BC_objects <- readRDS(file.path(RESULTS_DIR, "fig6", "fig6BC_objects.rds"))
X_mosaic <- fig6BC_objects$X_mosaic
X_expr   <- fig6BC_objects$X_expr
y        <- fig6BC_objects$y
cv_mos   <- fig6BC_objects$cv_mos
cv_exp   <- fig6BC_objects$cv_exp

## Lasso fits on all patients
fit_expr <- glmnet(X_expr, y, family = "binomial", alpha = 1, standardize=TRUE)
coef_expr <- coef(fit_expr, s = cv_exp$lambda.min)

fit_mos <- glmnet(X_mosaic, y, family = "binomial", alpha = 1, standardize=TRUE)
coef_mos <- coef(fit_mos, s = cv_mos$lambda.min)

## Abundance (pseudobulk) features
df_expr <- data.frame(
  Feature = rownames(coef_expr),
  Weight = as.numeric(coef_expr)
) %>%
  dplyr::filter(Feature != "(Intercept)", Weight != 0) %>%
  mutate(AbsWeight = abs(Weight)) %>%
  arrange(desc(AbsWeight)) %>%
  slice_head(n = 15) %>%
  mutate(Type = ifelse(Feature == "USP10", "Shared", "Abundance-Specific")) %>%
  mutate(RelWeight = (AbsWeight / max(AbsWeight)) * 100)

# Ascending order for plotting
df_expr$Feature <- factor(df_expr$Feature, levels = df_expr$Feature[order(df_expr$RelWeight)])

## Connectivity (MOSAIC) features
# Gene-level importance: sum of |coefficient| over embedding dimensions
mosaic_feature_name <- fig6BC_objects$mosaic_feature_name
dim_r <- fig6BC_objects$dim_r

df_mos_raw <- data.frame(
  Feature = c("(Intercept)", rep(mosaic_feature_name, dim_r)),
  Weight = as.numeric(coef_mos)
) %>%
  dplyr::filter(Feature != "(Intercept)", Weight != 0)

df_mos_agg <- df_mos_raw %>%
  mutate(Gene = gsub("_\\d+$", "", Feature)) %>% 
  group_by(Gene) %>%
  summarise(TotalWeight = sum(abs(Weight))) %>%
  arrange(desc(TotalWeight)) %>%
  slice_head(n = 15) %>%
  rename(Feature = Gene, AbsWeight = TotalWeight) %>%
  mutate(Type = ifelse(Feature == "USP10", "Shared", "Connectivity-Specific")) %>%
  mutate(RelWeight = (AbsWeight / max(AbsWeight)) * 100)

# Ascending order for plotting
df_mos_agg$Feature <- factor(df_mos_agg$Feature, levels = df_mos_agg$Feature[order(df_mos_agg$RelWeight)])

## Annotations
# Connectivity annotations (right)
anno_mos <- data.frame(
  Feature = c("KCNK6", "USP4", "SBF2", "HLA-DRB5", "USP10"), 
  Label   = c("Inflammasome Regulation", "Ubiquitin / Proteostasis", "Vesicle Trafficking", "Antigen Presentation", "Dual-Driver"),
  Color   = c("grey30", "grey30", "grey30", "grey30", "#D55E00"), 
  Font    = c("italic", "italic", "italic", "italic", "italic"),
  
  # 1 = right (inside), 0 = left (outside)
  Hjust   = c(1, 0, 0, 0, 0), 
  Nudge   = c(-2, 3, 3, 3, 3),
  
  # Inside labels sit above the line
  Vjust   = c(-0.6, 0.5, 0.5, 0.5, 0.5) 
) %>%
  left_join(df_mos_agg, by = "Feature")

# Abundance annotations (left)
anno_expr <- data.frame(
  Feature = c("UBAC2", "PRDX2", "USP10"), 
  Label   = c("Protein Quality Control", "Oxidative Stress", "Dual-Driver"),
  Color   = c("grey30", "grey30", "#D55E00"),
  Font    = c("italic", "italic", "italic"),
  
  Hjust   = c(1, 0, 0), 
  Nudge   = c(-2, 3, 3),
  
  Vjust   = c(-0.6, 0.5, 0.5)
) %>%
  left_join(df_expr, by = "Feature")

## Lollipop plot
plot_lollipop <- function(data, color_map, title_text, annotation_df = NULL) {
  
  p <- ggplot(data, aes(x = RelWeight, y = Feature)) +
    
    # Stick
    geom_segment(aes(x = 0, xend = RelWeight, y = Feature, yend = Feature), color = "grey70") +
    
    # Head
    geom_point(aes(color = Type), size = 4) +
    scale_color_manual(values = color_map) +
    
    # Annotations
    {if (!is.null(annotation_df)) 
      geom_text(data = annotation_df, 
                aes(x = RelWeight + Nudge, y = Feature, label = Label, 
                    hjust = Hjust, vjust = Vjust),
                color = annotation_df$Color, 
                fontface = annotation_df$Font,
                size = 3.5) 
    } +
    
    # Theme & Layout
    theme_bw(base_size = 14) +
    labs(title = title_text, x = "Relative Importance (%)", y = "") +
    theme(
      axis.text.y = element_text(face = "bold", color = "black"),
      panel.grid.major.y = element_line(color = "grey95"),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position = "none",
      plot.margin = margin(t=5, r=10, b=5, l=5) # Tight margins
    ) +
    # Expand x axis to fit labels
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) 
  
  return(p)
}

## Figure 6E
colors_expr <- c("Abundance-Specific" = "grey40", "Shared" = "#D55E00")
colors_mos  <- c("Connectivity-Specific" = "#E69F00", "Shared" = "#D55E00")

p_left  <- plot_lollipop(df_expr, colors_expr, "Top Abundance Drivers", anno_expr)
p_right <- plot_lollipop(df_mos_agg, colors_mos, "Top Connectivity Drivers", anno_mos)

# Top-15 features selected by both models
genes_expr_nz <- setdiff(rownames(coef_expr)[as.numeric(coef_expr) != 0], "(Intercept)")
cat("Top-15 features selected by both models (abundance list): ",
    paste(intersect(as.character(df_expr$Feature), unique(df_mos_raw$Feature)), collapse = ", "), "\n")
cat("Top-15 features selected by both models (connectivity list): ",
    paste(intersect(as.character(df_mos_agg$Feature), genes_expr_nz), collapse = ", "), "\n")

final_plot <- plot_grid(p_left, p_right, ncol = 2, align = "h")

print(final_plot)

png(file.path(FIG_DIR, "fig6E.png"),width = 6.5, height = 5, units = "in", res = 600)
print(final_plot)
dev.off()

## Figure 6D: Venn diagram of non-zero coefficients
# Expression model
genes_expr <- rownames(coef_expr)[which(as.numeric(coef_expr) != 0)]
genes_expr <- genes_expr[genes_expr != "(Intercept)"]

# MOSAIC model
genes_mos <- unique(df_mos_raw$Feature)

venn_list <- list(
  Expression = genes_expr,
  MOSAIC = genes_mos
)

# Jaccard index
intersection <- length(intersect(genes_expr, genes_mos))
union_len <- length(union(genes_expr, genes_mos))
jaccard <- round(intersection / union_len, 3)

print(paste("Jaccard Similarity Index:", jaccard))
print(paste("Shared Genes:", intersection))

# Calculate counts
n1 <- length(setdiff(venn_list[[1]], venn_list[[2]]))
n2 <- length(setdiff(venn_list[[2]], venn_list[[1]]))
n12 <- length(intersect(venn_list[[1]], venn_list[[2]]))

p_venn <- ggVennDiagram(venn_list, 
                        label_alpha = 0,
                        # Trailing spaces shift labels left
                        category.names = c("Abundance          \nFeatures          ", 
                                           "Connectivity          \nFeatures          ")) +
  scale_fill_gradient(low = "white", high = "darksalmon") +
  labs(
    title = "Lasso-Selected Predictive Features",
    subtitle = paste0("Jaccard Overlap Index: ", jaccard)
  ) +
  # Axis expansion so labels are not clipped
  scale_x_continuous(expand = expansion(mult = 0.4)) + 
  theme(plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

print(p_venn)

png(file.path(FIG_DIR, "fig6D.png"),width = 4, height = 5.7, units = "in", res = 600)
print(p_venn)
dev.off()
