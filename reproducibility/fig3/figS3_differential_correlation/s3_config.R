# Supplemental Figure S3: count simulation, shared configuration (sourced by the stage scripts)

suppressPackageStartupMessages({
  library(digest)
})

S3_CFG <- list(
  # simulation
  F           = 200L,     # structured genes
  r           = 20L,      # equal-size shared programs (modules)
  S           = 20L,      # samples; first S/2 condition A, rest condition B
  n_cells     = 300L,     # cells per sample
  m           = 4L,       # private partners per DC gene (regular bipartite graph)
  V           = 1.5,      # total log-rate variance of every structured gene
  s2          = 0.3,      # shared-program loading variance
  tau         = 0.1,      # donor amplitude sd (relative), truncated at 3 sd
  sigma0_sq   = 0.2,      # anchor (gene-specific coordinate) variance
  q           = 0.1,      # capture efficiency (binomial capture)
  qn          = 3L,       # sample-private programs per sample
  K_bg        = 500L,     # background genes (normalisation only)
  bg_mu       = 54,       # background mean transcripts per gene per cell (~90% of UMIs)
  bg_sd       = 0.1,      # background log-rate sd
  mu_median   = 10,       # structured-gene median transcripts per cell
  mu_sdlog    = 0.3,
  lib_sdlog   = 0.3,      # library-size log-normal sd (mean 1)
  dc_modules  = c(1L, 2L),
  rho2_floor  = 0,        # residual variance must stay strictly above this
  kappas      = c(0, 0.25, 0.5, 0.75, 1),
  etas        = c(0, 0.1, 0.2, 0.3),
  kappa_eta_panel = 0.5,

  # MOSAIC
  n_eig       = 50L,      # RSpectra::eigs_sym(k = 50, which = "LA")

  # CRN reference
  ref_donors  = 4L,       # per condition; must be even (split halves)
  ref_cells   = 2000L,

  # memento
  memento_commit   = "224c5c982f952cb1fd0f0de4623852dd435dbc1d",
  memento_env      = "r34_memento",
  memento_python   = "3.11",          # Python version used for the memento env
  memento_boot     = 100L,     # bootstrap replicates
  memento_chunks   = 7L,       # single-CPU tasks per cohort
  memento_min_cell_count = 20L,
  memento_min_perc_group = 0.7,

  # compute
  n_workers   = 70L,      # parallel workers
  pilot_seeds = 1:100,    # reserved for Stages 1 to 3
  final_seed0 = 1001L     # final runs start here
)
stopifnot(S3_CFG$F %% S3_CFG$r == 0, S3_CFG$m < S3_CFG$F / S3_CFG$r, S3_CFG$ref_donors %% 2 == 0)
S3_CFG$n_m <- S3_CFG$F %/% S3_CFG$r

# conservative budget capacity: every gene, sample and eta up to eta_max
beta2_max <- function(cfg, eta_max = max(cfg$etas)) {
  cfg$V - cfg$s2 * (1 + 3 * cfg$tau)^2 - cfg$sigma0_sq - eta_max - cfg$rho2_floor
}

# Seed namespaces
S3_NAMESPACES <- c("mu", "partner_perm", "donor", "nuisance", "cell_latent",
                    "library", "count", "background", "memento_boot", "crn_reference")

seed_for <- function(cohort_seed, namespace, ...) {
  stopifnot(namespace %in% S3_NAMESPACES)
  key <- paste(c("r34", cohort_seed, namespace, ...), collapse = "|")
  digest::digest2int(key)
}

# run `expr` under a fixed seed without disturbing the caller's RNG state
with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old, envir = .GlobalEnv)
  })
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  expr
}

# memento runs as a separate Python process; the interpreter is taken from S3_PYTHON
setup_memento_env <- function(cfg, python = Sys.getenv("S3_PYTHON", "")) {
  if (!nzchar(python))
    stop("S3_PYTHON is not set. Set it to a Python executable with memento installed from commit ",
         cfg$memento_commit, ", e.g.\n  export S3_PYTHON=/path/to/env/bin/python\n",
         "  (install: pip install \"git+https://github.com/yelabucsf/scrna-parameter-estimation.git@",
         cfg$memento_commit, "\")", call. = FALSE)
  if (!file.exists(python)) stop("Python executable not found (S3_PYTHON): ", python)
  verify_memento_commit(python, cfg$memento_commit)
  python
}

# Check that memento-de was installed from the pinned commit
verify_memento_commit <- function(python, commit) {
  code <- "import importlib.metadata as m; print(m.distribution('memento-de').read_text('direct_url.json') or '')"
  txt <- suppressWarnings(system2(python, c("-c", shQuote(code)), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(txt, "status")) || !length(txt) || !nzchar(paste(txt, collapse = "")))
    stop("memento-de not installed from the pinned commit in ", python, ":\n", paste(txt, collapse = "\n"))
  info <- jsonlite::fromJSON(paste(txt, collapse = "\n"))
  archive_url <- paste0("https://github.com/yelabucsf/scrna-parameter-estimation/archive/", commit, ".zip")
  ok <- identical(info$vcs_info$commit_id, commit) ||
        (!is.null(info$archive_info) && identical(info$url, archive_url))
  if (!ok) stop(sprintf("memento install (%s) does not match pinned commit %s",
                        if (is.null(info$vcs_info$commit_id)) info$url else info$vcs_info$commit_id, commit))
  invisible(commit)
}
