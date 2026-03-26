# Figure 4, S4, S13: CITE-seq preprocessing (CD4 naive T cells, Day 0 and Day 7)
# Input:    CITESEQ_H5SEURAT
# Output:   seurat_rna, seurat_adt, feature_names, feature_names_no_preflix (in the calling environment)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(Matrix)
})

# Load the data
pbmc <- LoadH5Seurat(CITESEQ_H5SEURAT)

# subset to certain cell type and timepoint
seurat_subset <- subset(pbmc, idents = "CD4 Naive")
seurat_subset <- seurat_subset[,which(seurat_subset$time %in% c(0,7))]
rm(pbmc); invisible(gc())

# Feature selection
counts <- GetAssayData(seurat_subset, assay = "SCT", slot = "counts")
sample_ids <- unique(as.character(seurat_subset$orig.ident))
list_of_expressed_genes <- lapply(sample_ids, function(id) {
  cells_in_sample <- which(seurat_subset$orig.ident == id)
  sample_counts <- counts[, cells_in_sample, drop = FALSE]
  expressed_genes <- rownames(sample_counts)[rowSums(sample_counts) > 3]
  return(expressed_genes)
})
genes_to_keep <- Reduce(intersect, list_of_expressed_genes)

counts <- GetAssayData(seurat_subset, assay = "ADT", slot = "counts")
sample_ids <- unique(as.character(seurat_subset$orig.ident))
list_of_expressed_genes <- lapply(sample_ids, function(id) {
  cells_in_sample <- which(seurat_subset$orig.ident == id)
  sample_counts <- counts[, cells_in_sample, drop = FALSE]
  expressed_genes <- rownames(sample_counts)[rowSums(sample_counts) > 3]
  return(expressed_genes)
})
protein_to_keep <- Reduce(intersect, list_of_expressed_genes)

seurat_subset[["SCT"]] <- subset(seurat_subset[["SCT"]], features = genes_to_keep)
seurat_subset[["ADT"]] <- subset(seurat_subset[["ADT"]], features = protein_to_keep)

seurat_subset$cluster_id <- "CD4 Naive"
seurat_subset$sample_id <- seurat_subset$orig.ident
seurat_subset$group_id <- seurat_subset$time

feature_names <- c(paste0("RNA_",rownames(seurat_subset@assays$SCT$counts)),paste0("ADT_",rownames(seurat_subset@assays$ADT$counts)))
feature_names_no_preflix <- c(rownames(seurat_subset@assays$SCT$counts),rownames(seurat_subset@assays$ADT$counts))

seurat_rna <- seurat_subset
seurat_rna[["ADT"]] <- NULL
seurat_adt <- seurat_subset
DefaultAssay(seurat_adt) <- "ADT"
seurat_adt[["SCT"]] <- NULL

seurat_rna <- NormalizeData(seurat_rna)
seurat_adt <- NormalizeData(seurat_adt, assay = "ADT", normalization.method = "CLR")

rm(seurat_subset, counts, list_of_expressed_genes); invisible(gc())
