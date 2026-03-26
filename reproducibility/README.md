# Analysis scripts for the MOSAIC paper

Scripts for the analyses and figures of

> Lu C, Kluger Y, Ma R. MOSAIC learns unified feature and sample joint embeddings for population-scale single-cell multi-omics.

## Usage

Set the data paths in `config.R` (or the environment variables it reads), then run scripts from this directory:

```bash
cd reproducibility
Rscript fig4/fig4_dc_paired_test.R
Rscript fig4/fig4_panels.R
```

Intermediate results are written to `results/` and figures to `figures/` (S3: `fig3/figS3_differential_correlation/stage4/`). Shared helper functions are in `R/functions.R`.

## Data

| Dataset | Figures | Source |
|---|---|---|
| CITE-seq vaccination cohort (Hao et al. 2021) | 4, S4, S13 | GEO GSE164378; h5Seurat file from https://atlas.fredhutch.org/nygc/multimodal-pbmc/ |
| COVID-19 PBMC atlas (Ren et al. 2021) | 6 | GEO GSE158055 |
| PFC and insular cortex multiome cohort (Yang et al. 2025) | 2B–D, 5, S1, S5–S9, S14 | NeMO Archive, https://scorch-portal.nemoarchive.org/, study UM1DA051410_spudich |

`data_preprocessing/pfc_l23it_preprocessing.R` extracts the PFC L2/3 IT neurons used in Figure 5. Figures 2A, 3, S2, S3, S10, S11 and S12 use simulated data.

## Scripts

Scripts are grouped by main figure; supplemental figures are in the folder of the related main figure.

| Figure | Script | Output |
|---|---|---|
| 1A | schematic, no code | |
| 1B–D | `fig1/fig1BCD_schematics.R` | `fig1BCD.png` |
| 2A | `fig2/fig2A_feature_embedding_benchmark.R` | `fig2A.png` |
| 2B | `fig2/fig2B_frobenius.R`, then `fig2/fig2B_pairs.R` | `fig2B_left.png`, `fig2B_right.png` |
| 2C, 2D | `fig2/fig2CD_figS1_sample_embedding_benchmark.R` | `fig2C.png`, `fig2D.png` |
| 3A–D | `fig3/fig3_dc_simulation_benchmark.R` | `fig3A.png` … `fig3D.png` |
| 4A–D | `fig4/fig4_dc_paired_test.R`, then `fig4/fig4_panels.R` | `fig4.png`, `fig4A.png` … `fig4D.png` |
| 5A, 5B, 5E | `fig5/fig5_figS6_S7_hiv_subgroup_discovery.R` | `fig5A.png`, `fig5B.png`, `fig5E.png` |
| 5C | `fig5/fig5C_figS9_volcano.R` | `fig5C.png` |
| 5D | `fig5/fig5D_cross_region_validation.R` (after the 5A script) | `fig5D.png` |
| 6A | `fig6/fig6A_lineage_screen.R` | `fig6A.png` |
| 6B, 6C | `fig6/fig6BC_prediction_benchmark.R` | `fig6B.png`, `fig6C.png` |
| 6D, 6E | `fig6/fig6DE_predictive_features.R` (after the 6B script) | `fig6D.png`, `fig6E.png` |
| S1 | `fig2/fig2CD_figS1_sample_embedding_benchmark.R` | `figS1.png` |
| S2 | `fig3/figS2_null_calibration.R`, then `fig3/figS2_plot.R` | `figS2.png` |
| S3 | `fig3/figS3_differential_correlation/` (see its README) | `stage4/figS3.png` |
| S4 | `fig4/figS4_top_features_mds.R` (after `fig4/fig4_dc_paired_test.R`) | `figS4.png` |
| S5A, S5C | `fig5/figS5_hiv_control_validation.R` | `figS5A.png`, `figS5C.png` |
| S5B | `fig5/figS5B_mds.R` | `figS5B.png` |
| S6, S7 | `fig5/fig5_figS6_S7_hiv_subgroup_discovery.R` | `figS6.png`, `figS7.png` |
| S8 | `fig5/figS8_pseudobulk_pca.R` (after the 5A script) | `figS8A.png`, `figS8B.png` |
| S9 | `fig5/fig5C_figS9_volcano.R` | `figS9.png` |
| S10 | `fig2/figS10_four_way_ablation.R` | `figS10.png` |
| S11 | `fig2/figS11_local_to_global.R` | `figS11.png` |
| S12 | `fig3/figS12_composition_simulation.R`, then `fig3/figS12_plot.R` | `figS12.png` |
| S13 | `fig4/figS13_dimension_sensitivity.R` (after `fig4/fig4_dc_paired_test.R`) | `figS13.png` |
| S14 | `fig5/figS14_dimension_sensitivity.R` | `figS14.png` |

Simulation scripts with command-line options list the commands used in their headers. `fig3/figS12_verify_exact_null.R` checks the exact permutation test of S12.

## Software

R 4.1.0 with Seurat 5.2.1, Signac 1.14.0, ArchR 1.0.3, SeuratDisk 0.0.0.9021, RSpectra 0.16.2, vegan 2.6-8, limma 3.50.3, DESeq2 1.34.0, edgeR 3.36.0, MAST 1.20.0, glmnet 4.1-8, pROC 1.18.5, dynamicTreeCut 1.63-1, ReactomePA 1.38.0, uwot 0.2.2, qs2 0.2.1, ggplot2 3.5.1, patchwork 1.3.0. Figure S3 also uses DGCA 1.0.3 and memento (Python, commit 224c5c9). MOFA+ and SIMBA (Figure 2) run in a Python environment set by `MOFA_SIMBA_PYTHON` in `config.R`.
