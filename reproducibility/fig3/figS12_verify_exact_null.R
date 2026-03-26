# Supplemental Figure S12: check of the exact permutation null against brute-force enumeration

source("config.R")

commandArgs <- function(...) "--norun"
source("fig3/figS12_composition_simulation.R")

## Small case: all choose(S, S/2) assignments with a direct pseudo-F
set.seed(42)
S <- 10L; d <- 3L; F_ <- 7L
E <- array(rnorm(F_ * d * S), dim = c(F_, d, S))
groups <- rep(c("A", "B"), each = S / 2)

fast <- exact_null_test(E, groups, chunk = 37L)

## Direct pseudo-F
directF <- function(Y, g) {
  nA <- sum(g); nB <- length(g) - nA
  mA <- colMeans(Y[g, , drop = FALSE]); mB <- colMeans(Y[!g, , drop = FALSE])
  m  <- colMeans(Y)
  ssb <- nA * sum((mA - m)^2) + nB * sum((mB - m)^2)
  sst <- sum(sweep(Y, 2, m, "-")^2)
  ssb / ((sst - ssb) / (length(g) - 2))
}
Ys <- lapply(seq_len(F_), function(f) t(matrix(E[f, , ], nrow = d)))
gobs <- groups == "A"
F_obs_direct <- vapply(Ys, directF, 1.0, g = gobs)
cat("F_obs max abs diff:", max(abs(F_obs_direct - fast$F_obs)), "\n")

## Brute-force p over all choose(S, S/2) assignments
all_combs <- combn(S, S / 2L)
Fperm <- matrix(0, ncol(all_combs), F_)
for (b in seq_len(ncol(all_combs))) {
  g <- logical(S); g[all_combs[, b]] <- TRUE
  Fperm[b, ] <- vapply(Ys, directF, 1.0, g = g)
}
p_brute <- colMeans(Fperm >= rep(F_obs_direct - 1e-10, each = ncol(all_combs)))
cat("p max abs diff:", max(abs(p_brute - fast$p)), "\n")
cat("p_brute:", round(p_brute, 5), "\n")
cat("p_fast :", round(fast$p, 5), "\n")
cat("n_perm_exact:", fast$n_perm_exact, " choose(S-1,S/2-1):", choose(S-1, S/2-1), "\n")

## Complement invariance of pseudo-F
gc1 <- logical(S); gc1[1:5] <- TRUE
cat("complement invariance max diff:",
    max(abs(vapply(Ys, directF, 1.0, g = gc1) - vapply(Ys, directF, 1.0, g = !gc1))), "\n")

## Uniformity of exact p under a true null
set.seed(7)
E2 <- array(rnorm(200 * 4 * 10), dim = c(200, 4, 10))
p2 <- exact_null_test(E2, rep(c("A","B"), each = 5))$p
cat("null FPR@0.05:", mean(p2 < 0.05), " mean p:", round(mean(p2), 4), "\n")
