# Figure 4: MOSAIC embedding and exact within-donor DC permutation test (CITE-seq, CD4 naive, Day 0 vs Day 7)
# Input:    CITESEQ_H5SEURAT
# Output:   results/fig4/dc_corrected_paired.rds, DC_table_corrected.csv, seurat_cd4naive.rds, mosaic_embedding.rds, coupling_rows_STAT5B.rds

source("config.R")
source("R/functions.R")   # eigen_block_cosine_two, find_elbow_kneedle, plot_eigen

suppressPackageStartupMessages({
  library(Seurat); library(SeuratDisk); library(Matrix)
  library(ggplot2); library(dplyr); library(vegan)
  library(foreach); library(doParallel)
})

out_dir <- file.path(RESULTS_DIR, "fig4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## 0. Data and preprocessing
source("fig4/fig4_preprocess.R")

saveRDS(list(seurat_rna = seurat_rna, seurat_adt = seurat_adt,
             feature_names = feature_names,
             feature_names_no_preflix = feature_names_no_preflix),
        file.path(out_dir, "seurat_cd4naive.rds"))

## 1. MOSAIC embedding
run_our_method_two <- function(rna_seurat, atac_seurat, assay1 = "SCT", assay2="ADT",sample_meta="sample_id",condition_meta="condition",block_method="cosine"){
  # Global scale
  rna_seurat <- ScaleData(rna_seurat,features=rownames(rna_seurat))
  atac_seurat <- ScaleData(atac_seurat,features=rownames(atac_seurat))
  
  rna.list <- SplitObject(rna_seurat, split.by = sample_meta)
  atac.list <- SplitObject(atac_seurat, split.by = sample_meta)
  
  dataset_list <- list()
  for (i in 1:length(rna.list)){
    print(i)
    gene_seurat <- rna.list[[i]]
    gene_matrix <- gene_seurat[[assay1]]$scale.data
    gene_matrix <- t(as.matrix(gene_matrix))
    
    atac_seurat <- atac.list[[i]]
    atac_matrix <- atac_seurat[[assay2]]$scale.data
    atac_matrix <- t(as.matrix(atac_matrix))
    
    dataset_list[[i]] <- list(gene_matrix, atac_matrix)
  }
  
  eigen_output_cosine_list <- list()
  r_cosine_list <- c()
  block_matrix_list <- list()
  for (i in 1:length(dataset_list)){
    print(i)
    matrix_set <- dataset_list[[i]]
    
    if (block_method == "cosine"){
      print("cosine")
      eigen_output_cosine <- eigen_block_cosine_two(matrix_set[[1]],matrix_set[[2]])
    } else if (block_method == "gaussian"){
      print("gaussian")
      eigen_output_cosine <- eigen_block_gaussian_two(matrix_set[[1]],matrix_set[[2]])
    } else if (block_method == "emd"){
      print("emd")
      eigen_output_cosine <- eigen_block_emd_two(matrix_set[[1]],matrix_set[[2]])
    } else {
      print("block method not found")
    }
    
    r_cosine_list <- c(r_cosine_list,eigen_output_cosine[[2]])
    eigen_output_cosine_list[[i]] <- eigen_output_cosine[c(1,2)]
    block_matrix_list[[i]] <- eigen_output_cosine[[3]]
  }
  r <- round(median(r_cosine_list))
  print(paste0("r: ", r))
  names(eigen_output_cosine_list) <- names(rna.list)
  names(block_matrix_list) <- names(rna.list)
  
  sample_name <- names(eigen_output_cosine_list)
  
  
  sample_condition_map <- unique(data.frame(
    sample_id = rna_seurat[[sample_meta]],
    condition = rna_seurat[[condition_meta]]
  ))
  
  colnames(sample_condition_map) <- c("sample_id","condition")
  
  annotation <- data.frame(
    Condition = sample_condition_map$condition
  )
  rownames(annotation) <- sample_condition_map$sample_id
  
  annotation <- annotation[names(eigen_output_cosine_list), , drop = FALSE]
  
  
  normalized_r_eigenvectors <- list()
  for (i in 1:length(eigen_output_cosine_list)){
    r_to_use <- r_cosine_list[i]
    if (r_to_use == 1){
      r_to_use <- 2
    }
    if (block_method == "gaussian"){
      original_matrix <- eigen_output_cosine_list[[i]][[1]]$vectors[,2:r_to_use]
    } else {
      original_matrix <- eigen_output_cosine_list[[i]][[1]]$vectors[,1:r_to_use]
    }
    normalized_r_eigenvectors[[i]] <- original_matrix
  }
  
  
  names(normalized_r_eigenvectors) <- names(eigen_output_cosine_list)
  
  sample_projection_matrix <- list()
  for (i in 1:length(eigen_output_cosine_list)){
    r_to_use <- r_cosine_list[i]
    if (block_method == "gaussian"){
      U_U_transpose <- normalized_r_eigenvectors[[i]] %*% diag(eigen_output_cosine_list[[i]][[1]]$values[2:r_to_use]) %*% t(normalized_r_eigenvectors[[i]])
    } else {
      U_U_transpose <- normalized_r_eigenvectors[[i]] %*% diag(eigen_output_cosine_list[[i]][[1]]$values[1:r_to_use]) %*% t(normalized_r_eigenvectors[[i]])
    }
    sample_projection_matrix[[i]] <- U_U_transpose
  }
  
  sum_projection_matrix <- Reduce("+", sample_projection_matrix)
  
  V_result <- RSpectra::eigs_sym(sum_projection_matrix, 50, which = "LA")
  
  plot_eigen(V_result$values)
  
  r2 <- find_elbow_kneedle(V_result$values[1:50])
  
  print(paste0("r2: ",r2))
  
  V_matrix <- V_result$vectors[,1:r2] %*% diag(V_result$values[1:r2])
  
  projected_combed_list <- list()
  for (i in 1:length(sample_projection_matrix)){
    projected_combed_list[[i]] <- block_matrix_list[[i]] %*% V_matrix
  }
  names(projected_combed_list) <- names(eigen_output_cosine_list)
  
  return(list(projected_combed_list=projected_combed_list,eigen_output_cosine_list=eigen_output_cosine_list,annotation=annotation,block_matrix_list=block_matrix_list))
}

no_shuffle_result <- run_our_method_two(seurat_rna,seurat_adt,sample_meta="sample_id",condition_meta="time",block_method = "cosine")

## Save the embedding without the coupling matrices; keep only their STAT5B rows
saveRDS(no_shuffle_result[c("projected_combed_list", "annotation", "eigen_output_cosine_list")],
        file.path(out_dir, "mosaic_embedding.rds"))
idx_stat5b <- which(feature_names == "RNA_STAT5B")
stopifnot(length(idx_stat5b) == 1L)
rn <- rownames(no_shuffle_result$block_matrix_list[[1]])
if (!is.null(rn)) stopifnot(identical(rn, feature_names_no_preflix))
coupling_rows <- lapply(no_shuffle_result$block_matrix_list, function(M) {
  v <- M[idx_stat5b, ]; names(v) <- feature_names; v
})
saveRDS(coupling_rows, file.path(out_dir, "coupling_rows_STAT5B.rds"))
no_shuffle_result$block_matrix_list <- NULL; rm(coupling_rows); invisible(gc())

## Exact within-donor permutation test (2^8 = 256 label swaps, embedding fixed)

emb  <- no_shuffle_result$projected_combed_list
ann  <- no_shuffle_result$annotation
emb  <- emb[rownames(ann)]                      # align to annotation order
S    <- length(emb)
stopifnot(S == 16L)

n_feat <- nrow(emb[[1]])
stopifnot(length(feature_names) == n_feat)
for (i in seq_len(S)) rownames(emb[[i]]) <- feature_names

## 2. Donor and condition (sample ids like "P3_0" / "P3_7")
sample_ids <- rownames(ann)
donor <- sub("_.*$", "", sample_ids)
cond  <- as.character(ann$Condition)
print(data.frame(sample = sample_ids, donor = donor, condition = cond))
print(table(donor, cond))
stopifnot(all(table(donor, cond) == 1L))
donors  <- sort(unique(donor))
n_donor <- length(donors)
stopifnot(n_donor == 8L)

## 3. Exact within-donor permutation and pseudo-F

## 3a. All 2^8 within-donor swaps as an S x 256 group-membership matrix
swap_grid <- as.matrix(expand.grid(rep(list(c(0L, 1L)), n_donor)))
colnames(swap_grid) <- donors
B_all <- nrow(swap_grid)
base_in_A <- as.integer(cond == cond[1])                   # observed labelling
A <- matrix(NA_integer_, nrow = S, ncol = B_all)
for (b in seq_len(B_all)) {
  flip <- donor %in% donors[swap_grid[b, ] == 1L]
  A[, b] <- ifelse(flip, 1L - base_in_A, base_in_A)
}
identity_col <- which(colSums(abs(A - base_in_A)) == 0)
stopifnot(length(identity_col) == 1L)
Ac <- 1L - A
n_per_group <- S / 2L

## 3b. adonis2-consistent pseudo-F, vectorised over all 256 label assignments
sq_dist <- function(X) {                    # squared Euclidean distances between rows
  cp <- tcrossprod(X); sq <- diag(cp)
  D2 <- outer(sq, sq, "+") - 2 * cp
  D2[D2 < 0] <- 0; diag(D2) <- 0
  D2
}
pseudoF_all <- function(D2) {
  SS_T <- sum(D2) / (2 * S)
  SS_W <- (colSums(A * (D2 %*% A)) + colSums(Ac * (D2 %*% Ac))) / (2 * n_per_group)
  SS_B <- SS_T - SS_W
  (SS_B / 1) / (SS_W / (S - 2))
}

## 3c. Check against adonis2 on a few features
chk <- sapply(1:5, function(f) {
  X <- t(sapply(emb, function(m) m[f, ]))
  as.numeric(adonis2(as.dist(sqrt(sq_dist(X))) ~ cond, permutations = 9)$F[1])
})
mine <- sapply(1:5, function(f) pseudoF_all(sq_dist(t(sapply(emb, function(m) m[f, ]))))[identity_col])
print(cbind(adonis2 = chk, vectorised = mine, diff = abs(chk - mine)))
stopifnot(max(abs(chk - mine)) < 1e-8)

## 3d. Pseudo-F for all features and all 256 permutations
n_cores <- N_CORES
E <- array(NA_real_, dim = c(S, ncol(emb[[1]]), n_feat))
for (i in seq_len(S)) E[i, , ] <- t(emb[[i]])
cl <- makeCluster(n_cores); registerDoParallel(cl)
F_perm <- foreach(idx = split(seq_len(n_feat), cut(seq_len(n_feat), n_cores, labels = FALSE)),
                  .combine = rbind, .export = c("sq_dist", "pseudoF_all", "A", "Ac",
                                                "S", "n_per_group")) %dopar% {
                                                  m <- matrix(NA_real_, length(idx), ncol(A))
                                                  for (k in seq_along(idx)) m[k, ] <- pseudoF_all(sq_dist(E[, , idx[k]]))
                                                  m
                                                }
stopCluster(cl)
rownames(F_perm) <- feature_names
F_obs <- F_perm[, identity_col]

## 4. Exact p-values, floor set, ranked table

## pseudo-F is symmetric in the two labels, so the smallest p-value is 2/256
exact_p <- rowMeans(F_perm >= F_obs)
p_floor <- 2 / B_all
adj_p   <- p.adjust(exact_p, method = "BH")

floor_set <- feature_names[abs(exact_p - p_floor) < 1e-12]
cat(sprintf("exact p floor = %.7f | features at floor = %d | BH < 0.05 = %d\n",
            p_floor, length(floor_set), sum(adj_p < 0.05)))
cat("modality of the floor set:\n"); print(table(substr(floor_set, 1, 3)))

DC_table <- data.frame(
  Feature      = feature_names,
  Name         = feature_names_no_preflix,
  Modality     = substr(feature_names, 1, 3),
  F_stat       = as.numeric(F_obs),
  R2           = NA_real_,
  exact_p      = as.numeric(exact_p),
  adj_p_BH     = as.numeric(adj_p),
  at_p_floor   = feature_names %in% floor_set,
  stringsAsFactors = FALSE
)
DC_table <- DC_table[order(-DC_table$F_stat), ]
DC_table$F_rank <- seq_len(nrow(DC_table))
rownames(DC_table) <- NULL

## Where the named candidates sit
for (g in c("RNA_STAT5B", "ADT_CD38-1")) {
  if (g %in% DC_table$Feature) {
    r <- DC_table[DC_table$Feature == g, ]
    cat(sprintf("%s: F = %.3f | rank %d | exact p = %.5f | at floor: %s\n",
                g, r$F_stat, r$F_rank, r$exact_p, r$at_p_floor))
  }
}
print(head(DC_table[, c("Feature", "F_stat", "F_rank", "exact_p", "at_p_floor")], 20))

## 5. Save
write.csv(DC_table, file.path(out_dir, "DC_table_corrected.csv"), row.names = FALSE)
saveRDS(list(F_perm = F_perm, F_obs = F_obs, exact_p = exact_p, adj_p = adj_p,
             A = A, identity_col = identity_col, p_floor = p_floor,
             floor_set = floor_set, donor = donor, cond = cond,
             sample_ids = sample_ids, feature_names = feature_names,
             DC_table = DC_table),
        file.path(out_dir, "dc_corrected_paired.rds"))
cat("written to", out_dir, "\n")
