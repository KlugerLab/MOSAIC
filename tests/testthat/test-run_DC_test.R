make_embed <- function(S, n_feat = 12, d = 3, seed = 1) {
  set.seed(seed)
  emb <- lapply(seq_len(S), function(i) {
    matrix(stats::rnorm(n_feat * d), n_feat, d,
           dimnames = list(paste0("f", seq_len(n_feat)), NULL))
  })
  names(emb) <- paste0("s", seq_len(S))
  emb
}

feature_matrix <- function(emb, f) t(vapply(emb, function(e) e[f, ], numeric(ncol(emb[[1]]))))

test_that("pseudo-F matches vegan::adonis2", {
  emb <- make_embed(9)
  g <- rep(c("A", "B", "C"), c(3, 2, 4))
  res <- run_DC_test(emb, g)$results
  ref <- vapply(seq_len(nrow(res)), function(f) {
    vegan::adonis2(stats::dist(feature_matrix(emb, f)) ~ g, permutations = 0)$F[1]
  }, numeric(1))
  expect_equal(res$F_stat, ref, tolerance = 1e-10)
})

test_that("cosine pseudo-F matches adonis2 on 1 - cosine similarity", {
  emb <- make_embed(8)
  g <- rep(c("A", "B"), 4)
  res <- run_DC_test(emb, g, dist_method = "cosine")$results
  sims <- compute_feature_similarity(emb)
  ref <- vapply(sims, function(m) {
    vegan::adonis2(stats::as.dist(1 - m) ~ g, permutations = 0)$F[1]
  }, numeric(1))
  expect_equal(res$F_stat, unname(ref), tolerance = 1e-10)
})

test_that("exact P-values equal brute-force enumeration", {
  emb <- make_embed(8)
  g <- rep(c("A", "B"), c(3, 5))
  dc <- run_DC_test(emb, g)
  expect_true(dc$exact)
  expect_equal(dc$n_perm, choose(8, 3))
  cmb <- utils::combn(8, 3)
  brute <- vapply(seq_len(nrow(dc$results)), function(f) {
    D <- stats::dist(feature_matrix(emb, f))
    Fs <- apply(cmb, 2, function(ix) {
      gg <- rep("B", 8)
      gg[ix] <- "A"
      vegan::adonis2(D ~ gg, permutations = 0)$F[1]
    })
    mean(Fs >= dc$results$F_stat[f] - 1e-9)
  }, numeric(1))
  expect_equal(dc$results$pvalue, brute)
  expect_equal(dc$results$padj, stats::p.adjust(brute, "BH"))
})

test_that("paired designs permute within blocks", {
  emb <- make_embed(10)
  g <- rep(c("d0", "d7"), 5)
  donor <- rep(seq_len(5), each = 2)
  dc <- run_DC_test(emb, g, blocks = donor, return_null = TRUE)
  expect_equal(dc$n_perm, 2^5)
  expect_equal(dc$min_pvalue, 2^(1 - 5))
  expect_true(all(dc$results$pvalue >= 2^(1 - 5) - 1e-12))
  # every assignment keeps one d0 and one d7 per donor
  per_donor <- apply(dc$assignments, 2, function(a) all(tapply(a, donor, function(x) all(sort(x) == c("d0", "d7")))))
  expect_true(all(per_donor))
})

test_that("Monte Carlo P-values use (1 + count) / (B + 1)", {
  emb <- make_embed(12)
  g <- rep(c("A", "B"), 6)
  dc <- run_DC_test(emb, g, max_exact = 0, n_perm = 199, seed = 3, return_null = TRUE)
  expect_false(dc$exact)
  expect_equal(dc$n_perm, 199)
  F_obs <- dc$results$F_stat
  expected <- (1 + rowSums(dc$F_null >= F_obs - 1e-9 * pmax(1, abs(F_obs)))) / 200
  expect_equal(dc$results$pvalue, unname(expected))
  dc2 <- run_DC_test(emb, g, max_exact = 0, n_perm = 199, seed = 3)
  expect_equal(dc$results$pvalue, dc2$results$pvalue)
})

test_that("chunking does not change results", {
  emb <- make_embed(9, n_feat = 25)
  g <- rep(c("A", "B"), c(4, 5))
  a <- run_DC_test(emb, g)$results
  b <- run_DC_test(emb, g, feature_chunk = 7, perm_chunk = 13)$results
  expect_equal(a, b)
})

test_that("invalid designs are rejected", {
  emb <- make_embed(6)
  expect_error(run_DC_test(emb, rep("A", 6)))
  expect_error(run_DC_test(emb, rep(c("A", "B"), 2)))
  expect_error(run_DC_test(emb, rep(c("A", "B"), each = 3), blocks = rep(1:2, each = 3)))
})

test_that("find_elbow_kneedle is bounded below by two", {
  expect_equal(find_elbow_kneedle(c(10, 1, 0.9, 0.8, 0.7, 0.6)), 2)
  expect_equal(find_elbow_kneedle(c(3, 2)), 2)
})

test_that("calculate_empirical_pvalue is deprecated", {
  expect_warning(calculate_empirical_pvalue(1, c(0, 2)), "deprecated")
})
