# Figure 2A, 2C-D, S1: SIMBA wrapper (run_simba_embedding) and Python set-up (MOFA_SIMBA_PYTHON)

if (!exists("MOFA_SIMBA_PYTHON")) source("config.R")
if (!nzchar(MOFA_SIMBA_PYTHON)) {
  stop("MOFA_SIMBA_PYTHON is empty. Set the environment variable ",
       "MOSAIC_MOFA_SIMBA_PYTHON (or edit config.R) to the python binary of an ",
       "environment with simba (simba-pbg), scanpy, mofapy2 and cupy.")
}

# nvidia libraries of the Python environment for GPU MOFA+; NA if absent
MOFA_SIMBA_NVIDIA_DIR <- Sys.glob(file.path(dirname(dirname(MOFA_SIMBA_PYTHON)),
                                            "lib", "python3*", "site-packages", "nvidia"))[1]
nvidia_lib_dirs <- if (!is.na(MOFA_SIMBA_NVIDIA_DIR)) {
  file.path(MOFA_SIMBA_NVIDIA_DIR, c("cublas", "cusolver", "cusparse"), "lib")
} else character(0)

Sys.setenv(
  RETICULATE_PYTHON = MOFA_SIMBA_PYTHON,
  LD_LIBRARY_PATH   = paste(
    c(nvidia_lib_dirs,
      "/usr/lib/x86_64-linux-gnu",
      "/usr/lib64",
      Sys.getenv("LD_LIBRARY_PATH")),
    collapse = ":"
  )
)

library(reticulate)
py_config()

# optional GPU check
try(py_run_string("import cupy as cp; print('CuPy:', cp.__version__); print('GPUs:', cp.cuda.runtime.getDeviceCount())"),
    silent = TRUE)

run_simba_embedding <- function(combined_seurat,
                                rna_assay = "originalexp",
                                atac_assay = "ATAC",
                                n_bins = 5,
                                n_dims = 50,
                                pbg_workers = 4,
                                pbg_epochs = 10,
                                workdir = tempdir(),
                                graph_name = "graph0",
                                model_name = "model") {
  
  rna_counts <- combined_seurat[[rna_assay]]$counts
  atac_counts <- combined_seurat[[atac_assay]]$counts
  rna_genes <- rownames(rna_counts)
  atac_peaks <- rownames(atac_counts)
  meta_df <- combined_seurat@meta.data
  
  py$rna_counts <- rna_counts
  py$atac_counts <- atac_counts
  py$rna_genes <- rna_genes
  py$atac_peaks <- atac_peaks
  py$meta_df <- meta_df
  py$workdir <- workdir
  py$n_bins <- as.integer(n_bins)
  py$n_dims <- as.integer(n_dims)
  py$pbg_workers <- as.integer(pbg_workers)
  py$pbg_epochs <- as.integer(pbg_epochs)
  py$graph_name <- graph_name
  py$model_name <- model_name
  
  py_run_string("
  


import simba as si
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix

si.settings.set_workdir(workdir)

def make_adata(counts, feature_names, metadata):
    if issparse(counts):
        X = csr_matrix(counts.T)
    else:
        X = csr_matrix(np.array(counts).T)
    adata = sc.AnnData(X=X, obs=metadata.copy())
    adata.var_names = pd.Index([str(g) for g in feature_names])
    return adata

adata_rna = make_adata(rna_counts, rna_genes, meta_df)
adata_atac = make_adata(atac_counts, atac_peaks, meta_df)

adata_rna.layers['raw'] = adata_rna.X.copy()
adata_atac.layers['raw'] = adata_atac.X.copy()

si.tl.discretize(adata_rna, n_bins=int(n_bins))
si.tl.discretize(adata_atac, n_bins=int(n_bins))

si.tl.gen_graph(list_CP=[adata_atac],
                list_CG=[adata_rna],
                copy=False,
                use_highly_variable=False,   # now uses selected features
                use_top_pcs=False,
                dirname=graph_name)

# --- Tune PBG training for speed ---
si.tl.pbg_train(auto_wd=True,
                save_wd=False,
                output=model_name)

si.load_graph_stats()
si.load_pbg_config()
dict_adata = si.read_embedding()

gene_emb_mat = dict_adata['G'].X
gene_emb_names = list(dict_adata['G'].obs_names)
peak_emb_mat = dict_adata['P'].X
peak_emb_names = list(dict_adata['P'].obs_names)
")
  
  gene_mat <- py$gene_emb_mat
  gene_names <- py$gene_emb_names
  peak_mat <- py$peak_emb_mat
  peak_names <- py$peak_emb_names
  
  colnames_emb <- paste0("dim_", seq_len(ncol(gene_mat)) - 1)
  
  gene_emb <- as.data.frame(gene_mat)
  rownames(gene_emb) <- gene_names
  colnames(gene_emb) <- colnames_emb
  
  peak_emb <- as.data.frame(peak_mat)
  rownames(peak_emb) <- peak_names
  colnames(peak_emb) <- colnames_emb
  
  rbind(gene_emb, peak_emb)
}
