#' Find the Elbow Point Using the Kneedle Algorithm
#'
#' Automatically selects the effective dimensionality from a decreasing
#' eigenvalue sequence. The algorithm log-transforms the eigenvalues, normalizes
#' both axes to the unit interval, draws a straight line between the first and last
#' points, and returns the index of the point with the greatest perpendicular
#' distance from that line. This point corresponds to the "elbow" where the
#' eigenvalue decay transitions from steep to flat.
#'
#' MOSAIC calls this function internally at two stages:
#' (1) to choose the per-sample rank from each coupling matrix eigendecomposition,
#' and (2) to choose the global rank from the aggregated projection.
#'
#' @param eigenvalues Numeric vector of eigenvalues, typically in decreasing
#'   order. Non-positive values are automatically replaced with the preceding
#'   value to ensure a valid log transform. Must have length >= 3 for a
#'   meaningful result; if shorter, the full length is returned.
#'
#' @return An integer: the 1-based index of the elbow point. For example, a
#'   return value of 5 means the first 5 eigenvalues capture the dominant
#'   signal.
#'
#' @examples
#' # Typical eigenvalue decay with a clear elbow at position 3
#' eigenvalues <- c(10, 5, 3, 2, 1.5, 1.2, 1.1, 1.05, 1.01, 1.0)
#' find_elbow_kneedle(eigenvalues)
#' # Returns 3
#'
#' # Visualize the elbow
#' plot(eigenvalues, type = "b", pch = 19, xlab = "Index", ylab = "Eigenvalue")
#' abline(v = find_elbow_kneedle(eigenvalues), col = "red", lty = 2)
#'
#' @seealso \code{\link{plot_eigen}} for visualizing eigenvalue spectra.
#'
#' @export
find_elbow_kneedle <- function(eigenvalues) {
  # Replace non-positive values with previous value
  for (i in seq_along(eigenvalues)) {
    if (eigenvalues[i] <= 0) {
      eigenvalues[i] <- eigenvalues[i - 1]
    }
  }

  x <- seq_along(eigenvalues)
  y <- log(eigenvalues)

  if (length(x) < 3) {
    return(length(x))
  }

  # Normalize to [0, 1]
  x_scaled <- (x - min(x)) / (max(x) - min(x))
  y_scaled <- (y - min(y)) / (max(y) - min(y))

  # Perpendicular distance from line connecting first and last point
  distances <- .perpendicular_distance(
    xp = x_scaled, yp = y_scaled,
    x1 = x_scaled[1], y1 = y_scaled[1],
    x2 = x_scaled[length(x_scaled)], y2 = y_scaled[length(y_scaled)]
  )

  which.max(distances)
}


#' Calculate Empirical P-Value
#'
#' Computes an empirical p-value by comparing an observed test statistic against
#' a null distribution. This is the standard approach in MOSAIC for calibrating
#' differential connectivity significance: run the DC pipeline on real labels
#' and on shuffled labels, then compare F-statistics.
#'
#' @param test_value Numeric scalar. The observed test statistic (e.g. a
#'   PERMANOVA F-statistic from
#'   \code{\link{run_DC_test}}).
#' @param data_list Numeric vector of null distribution values (e.g. all
#'   F-statistics from a label-shuffled run).
#' @param alternative Character string specifying the alternative hypothesis:
#'   \describe{
#'     \item{\code{"greater"}}{(Default) Tests whether \code{test_value} is
#'       unusually large. P-value = proportion of null values >= test_value.}
#'     \item{\code{"less"}}{Tests whether \code{test_value} is unusually small.
#'       P-value = proportion of null values <= test_value.}
#'     \item{\code{"two.sided"}}{Tests whether \code{test_value} is extreme in
#'       either direction. P-value = 2 * one-sided p-value, capped at 1.}
#'   }
#'
#' @return A numeric p-value between 0 and 1.
#'
#' @examples
#' # Simulate a null distribution and test an observed value
#' set.seed(42)
#' null_distribution <- rnorm(10000, mean = 0, sd = 1)
#' observed <- 2.5
#'
#' # One-sided test (is the observed value unusually large?)
#' calculate_empirical_pvalue(observed, null_distribution, alternative = "greater")
#'
#' # Two-sided test
#' calculate_empirical_pvalue(observed, null_distribution, alternative = "two.sided")
#'
#' # Typical MOSAIC usage:
#' # pval <- calculate_empirical_pvalue(F_observed, F_null_vector)
#'
#' @seealso \code{\link{run_DC_test}} which
#'   produces the F-statistics to be tested.
#'
#' @export
calculate_empirical_pvalue <- function(test_value, data_list, alternative = "greater") {
  if (alternative == "greater") {
    p_value <- mean(data_list >= test_value)
  } else if (alternative == "less") {
    p_value <- mean(data_list <= test_value)
  } else if (alternative == "two.sided") {
    if (test_value > stats::median(data_list)) {
      p_value <- 2 * mean(data_list >= test_value)
    } else {
      p_value <- 2 * mean(data_list <= test_value)
    }
    p_value <- min(p_value, 1)
  } else {
    stop("alternative must be 'greater', 'less', or 'two.sided'")
  }
  return(p_value)
}


# ---------------------------------------------------------------------------
# Internal: perpendicular distance from points to a line
# ---------------------------------------------------------------------------
.perpendicular_distance <- function(xp, yp, x1, y1, x2, y2) {
  A <- y1 - y2
  B <- x2 - x1
  C <- -A * x1 - B * y1
  abs(A * xp + B * yp + C) / sqrt(A^2 + B^2)
}


# ---------------------------------------------------------------------------
# Internal: remove all-zero rows from a matrix
# ---------------------------------------------------------------------------
.remove_zero_rows <- function(mat) {
  zero_rows <- apply(mat, 1, function(row) all(row == 0))
  list(mat[!zero_rows, , drop = FALSE], zero_rows)
}


# ---------------------------------------------------------------------------
# Internal: cosine similarity between two vectors
# ---------------------------------------------------------------------------
.cosine_sim_vectors <- function(x, y) {
  nx <- sqrt(sum(x^2))
  ny <- sqrt(sum(y^2))
  if (nx == 0 || ny == 0) return(0)
  sum(x * y) / (nx * ny)
}
