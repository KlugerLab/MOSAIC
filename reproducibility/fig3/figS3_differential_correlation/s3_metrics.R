# Supplemental Figure S3: performance metrics (sourced by the stage scripts)
suppressPackageStartupMessages(library(PRROC))

# AUPRC as in the Figure 3 benchmark (Inf scores set to the largest finite score + 1)
auprc_fig3 <- function(scores, labels) {
  if (anyNA(scores)) stop("NA scores")
  if (any(is.infinite(scores))) {
    max_score <- max(scores[is.finite(scores)], na.rm = TRUE)
    scores[is.infinite(scores)] <- max_score + 1
  }
  pr.curve(scores.class0 = scores, weights.class0 = as.numeric(labels))$auc.integral
}

# AUROC of x for DC (TRUE) vs null (FALSE), Mann-Whitney form, ties counted 1/2
auroc <- function(x, is_pos) {
  r <- rank(x); n1 <- sum(is_pos); n0 <- sum(!is_pos)
  (sum(r[is_pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# mean change of the DC genes' percentile rank from F_raw to F_M (positive = moved up)
dc_percentile_shift <- function(F_raw, F_M, is_dc) {
  p_raw <- rank(F_raw) / length(F_raw); p_m <- rank(F_M) / length(F_M)
  mean(p_m[is_dc] - p_raw[is_dc])
}
