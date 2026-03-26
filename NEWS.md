# MOSAIC 1.1.0

* `run_DC_test()` permutes condition labels across samples and tests each
  feature against its own null distribution. New arguments `blocks` (permute
  within blocks, e.g. donors in a paired design), `n_perm` and `max_exact`
  (exact enumeration when the number of assignments is small). Results are
  returned as a data frame with pseudo-F, R2, P-value and BH-adjusted P-value.
* `calculate_empirical_pvalue()` is deprecated; use `run_DC_test()`.
* New `compute_feature_similarity()` for per-feature sample similarity
  matrices.
* `find_elbow_kneedle()` returns at least 2.
* Analysis scripts for the paper are in `reproducibility/`.

# MOSAIC 1.0.0

* Initial release.
