# Paths used by all scripts. Run scripts from the reproducibility/ directory.
# Edit the paths below or set the corresponding environment variables.

.env_or <- function(var, default) {
  val <- Sys.getenv(var, unset = "")
  if (nzchar(val)) val else default
}

# ---- Data ----
DATA_DIR <- .env_or("MOSAIC_DATA_DIR", "data")

# CITE-seq vaccination cohort (GSE164378)
CITESEQ_H5SEURAT <- file.path(DATA_DIR, "multi.h5seurat")

# COVID-19 PBMC atlas (GSE158055)
COVID_H5AD <- file.path(DATA_DIR, "Covid_19_PBMC.h5ad")

# PFC / INS multiome cohort (NeMO, UM1DA051410_spudich)
PFC_SEURAT_RDS   <- file.path(DATA_DIR, "PFC_seurat.RDS")
PFC_L23IT_RNA    <- file.path(DATA_DIR, "PFC_L23-IT_rna.RDS")
PFC_L23IT_ATAC   <- file.path(DATA_DIR, "PFC_L23-IT_atac.RDS")
PFC_COMBINED_QS2 <- file.path(DATA_DIR, "rna_atac_all_combined_PFC.qs2")
INS_COMBINED_QS2 <- file.path(DATA_DIR, "rna_atac_all_combined_INS.qs2")

# ---- Outputs ----
RESULTS_DIR <- .env_or("MOSAIC_RESULTS_DIR", "results")
FIG_DIR     <- .env_or("MOSAIC_FIG_DIR", "figures")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Python ----
# mofapy2 and simba (Figure 2)
MOFA_SIMBA_PYTHON <- .env_or("MOSAIC_MOFA_SIMBA_PYTHON", "")
# memento (Supplemental Figure S3)
MEMENTO_PYTHON <- .env_or("S3_PYTHON", "")

N_CORES <- as.integer(.env_or("MOSAIC_N_CORES", max(1L, parallel::detectCores() - 1L)))
