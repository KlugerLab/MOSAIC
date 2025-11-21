### MOSAIC: A Spectral Framework for Integrative Phenotypic Characterization Using Population-Level Single-Cell Multi-Omics

***

**MOSAIC** (Multi-Omic Sample-wise Analysis of Inter-feature Connectivity) is a spectral framework designed to learn high-resolution **feature $\times$ sample** joint embeddings from population-scale single-cell multi-omics data.

Unlike traditional methods that focus on cell embeddings, MOSAIC explicitly models how **feature relationships** (e.g., Gene-Peak, Protein-Gene) vary across individuals. This enables the detection of regulatory network rewiring and cryptic patient subgroups.

![MOSAIC Framework](https://raw.githubusercontent.com/KlugerLab/MOSAIC/main/Figures/framework.png)

#### Key Capabilities

- **Joint Feature $\times$ Sample Embedding:** Projects features and samples into a shared latent space.
- **Differential Connectivity (DC) Analysis:** Identifies features (genes, proteins) that change their regulatory context between conditions, even without changes in abundance.
- **Unsupervised Subgroup Detection:** Discovers patient subtypes driven by specific, coherent multi-modal feature modules.
- **Scalable & Robust:** Linear complexity with respect to sample size (O(S)) and quadratic with respect to features (O(F^2)), utilizing truncated eigendecomposition. Robust to batch effects without explicit correction.

***

## Example tutorial

- The MOSAIC tutorial, available in [MOSAIC_demo.Rmd](https://github.com/KlugerLab/MOSAIC/blob/main/MOSAIC_demo.Rmd) and [MOSAIC_demo.html](https://github.com/KlugerLab/MOSAIC/blob/main/MOSAIC_demo.html).
- The functions defined MOSAIC and downstream analysis are available in [`MOSAIC_function.R`](https://github.com/KlugerLab/MOSAIC/blob/main/MOSAIC_function.R).
- Dataset used in the tutorial can be downloaded through this [link](https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat).
