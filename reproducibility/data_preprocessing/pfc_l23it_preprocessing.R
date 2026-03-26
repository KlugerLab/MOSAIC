# Data preprocessing: PFC L2/3 IT neurons, paired snRNA-seq and snATAC-seq
# Input:    DATA_DIR/seurat_combined_v2.RDS, DATA_DIR/PFC_signac_integrated_v2.RDS
# Output:   PFC_L23IT_RNA, PFC_L23IT_ATAC

source("config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
})

## Input objects
SCORCH_RNA_RDS      <- file.path(DATA_DIR, "seurat_combined_v2.RDS")
SCORCH_PFC_ATAC_RDS <- file.path(DATA_DIR, "PFC_signac_integrated_v2.RDS")

atac_seurat_PFC <- readRDS(SCORCH_PFC_ATAC_RDS)
rna_seurat <- readRDS(SCORCH_RNA_RDS)

# PFC L2/3 IT neurons
rna_seurat_subset <- rna_seurat[, which(rna_seurat$region == "PFC")]
rna_seurat_subset <- rna_seurat_subset[, which(rna_seurat_subset$consensus_cell_type_coarse_l1 == "L2/3 IT")]

# rename RNA cell names to the ATAC format: <sample>_<barcode>
rna_cell_names <- colnames(rna_seurat_subset)
new_rna_cell_names <- sapply(rna_cell_names, function(name) {
  parts <- unlist(strsplit(name, split = "-"))
  paste(parts[3], paste(parts[1], parts[2], sep = "-"), sep = "_")
})
names(new_rna_cell_names) <- NULL
rna_seurat_subset <- RenameCells(rna_seurat_subset, new.names = new_rna_cell_names)

# keep cells present in both modalities
rna_seurat_subset <- rna_seurat_subset[, Cells(atac_seurat_PFC)]
atac_seurat_subset <- atac_seurat_PFC[, Cells(rna_seurat_subset)]
stopifnot(all(Cells(rna_seurat_subset) == Cells(atac_seurat_subset)))

saveRDS(rna_seurat_subset, PFC_L23IT_RNA)
saveRDS(atac_seurat_subset, PFC_L23IT_ATAC)
