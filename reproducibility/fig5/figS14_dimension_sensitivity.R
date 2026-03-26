# Supplemental Figure S14: sensitivity of the HIV+ subgroup result to the embedding dimension d
# Input:    PFC_L23IT_RNA, PFC_L23IT_ATAC
# Output:   figures/figS14.png, results/S14/hiv_module_sensitivity.csv

## Paths and parameters
source("config.R")
FUNCTIONS_R <- "R/functions.R"
RNA_RDS     <- PFC_L23IT_RNA
ATAC_RDS    <- PFC_L23IT_ATAC
OUT_DIR     <- file.path(RESULTS_DIR, "S14")

MODULE_ID       <- 27
MIN_MODULE_SIZE <- 20
DEEP_SPLIT      <- 2
MIN_GROUP_FRAC  <- 0.25
N_EIG           <- 100
D0_WINDOW       <- 50
N_NULL_COH      <- 200
NCORE           <- N_CORES
SEED            <- 1

suppressPackageStartupMessages({
  library(Seurat); library(RSpectra); library(dynamicTreeCut)
  library(cluster); library(lsa); library(ggplot2); library(mclust)
  library(patchwork)
})
source(FUNCTIONS_R)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

HCL <- if (requireNamespace("fastcluster", quietly = TRUE))
  fastcluster::hclust else stats::hclust
if (requireNamespace("RhpcBLASctl", quietly = TRUE))
  RhpcBLASctl::blas_set_num_threads(if (NCORE > 1L) 1L else 8L)

`%ni%` <- Negate(`%in%`)
dl <- c(); lg <- function(...) { m <- paste0(...); message(m); dl <<- c(dl, m) }

## 1. Preprocessing

PREP <- file.path(OUT_DIR, "prep.rds")
if (file.exists(PREP)) {
  prep <- readRDS(PREP); lg("loaded cached preprocessing")
} else {
  rna  <- readRDS(RNA_RDS)
  atac <- readRDS(ATAC_RDS)
  DefaultAssay(atac) <- "GeneACT"
  atac[["ATAC"]] <- NULL
  
  sc   <- table(rna$sample_name)
  drop <- names(sc[sc < 10])
  rna  <- subset(rna,  subset = sample_name %ni% drop)
  atac <- subset(atac, subset = sample_name %ni% drop)
  rna  <- subset(rna,  subset = condition_new %ni% c("CTR"))
  atac <- subset(atac, subset = condition_new %ni% c("CTR"))
  
  cts <- GetAssayData(rna, assay = "RNA", slot = "counts")
  ids <- unique(as.character(rna$sample_name))
  keep_rna <- Reduce(intersect, lapply(ids, function(id) {
    m <- cts[, which(rna$sample_name == id), drop = FALSE]
    rownames(m)[rowSums(m) > 3]
  }))
  rna <- FindVariableFeatures(rna, nfeatures = 10000)
  keep_rna <- keep_rna[keep_rna %in% VariableFeatures(rna)]
  
  cts <- GetAssayData(atac, assay = "GeneACT", slot = "counts")
  ids <- unique(as.character(atac$sample_name))
  keep_atac <- Reduce(intersect, lapply(ids, function(id) {
    m <- cts[, which(atac$sample_name == id), drop = FALSE]
    rownames(m)[rowSums(m) > 5]
  }))
  atac <- FindVariableFeatures(atac, nfeatures = 3000)
  keep_atac <- keep_atac[keep_atac %in% VariableFeatures(atac)]
  
  seurat_atac <- NormalizeData(subset(atac, features = keep_atac))
  seurat_rna  <- NormalizeData(subset(rna,  features = keep_rna))
  feature_names <- c(rownames(seurat_rna@assays$RNA$counts),
                     paste0("ATAC-", rownames(seurat_atac@assays$GeneACT$counts)))
  prep <- list(rna = seurat_rna, atac = seurat_atac, feature_names = feature_names)
  saveRDS(prep, PREP)
}
feature_names <- prep$feature_names
lg("n_features = ", length(feature_names))

## 2. MOSAIC stage 1 and aggregation (N_EIG eigenpairs; d0 from values[1:D0_WINDOW])

BASIS <- file.path(OUT_DIR, "basis.rds")
if (file.exists(BASIS)) {
  bas <- readRDS(BASIS); lg("loaded cached basis")
} else {
  rna_s  <- ScaleData(prep$rna,  features = rownames(prep$rna))
  atac_s <- ScaleData(prep$atac, features = rownames(prep$atac))
  rna.list  <- SplitObject(rna_s,  split.by = "sample_name")
  atac.list <- SplitObject(atac_s, split.by = "sample_name")
  
  r_list <- c(); eig_list <- list(); block_list <- list()
  for (i in seq_along(rna.list)) {
    message("stage 1: sample ", i, "/", length(rna.list))
    gm <- t(as.matrix(rna.list[[i]][["RNA"]]$scale.data))
    am <- t(as.matrix(atac.list[[i]][["GeneACT"]]$scale.data))
    eo <- eigen_block_cosine_two(gm, am)
    r_list <- c(r_list, eo[[2]])
    eig_list[[i]]   <- eo[[1]]
    block_list[[i]] <- eo[[3]]
  }
  names(eig_list) <- names(block_list) <- names(rna.list)
  
  scm <- unique(data.frame(sample_id = rna_s[["sample_name"]],
                           condition = rna_s[["condition_new"]]))
  colnames(scm) <- c("sample_id", "condition")
  annotation <- data.frame(Condition = scm$condition,
                           row.names = scm$sample_id)[names(eig_list), , drop = FALSE]
  
  proj <- lapply(seq_along(eig_list), function(i) {
    r <- r_list[i]; if (r == 1) r <- 2
    V <- eig_list[[i]]$vectors[, 1:r]
    V %*% diag(eig_list[[i]]$values[1:r]) %*% t(V)
  })
  agg <- Reduce("+", proj); rm(proj); gc(FALSE)
  Vr  <- RSpectra::eigs_sym(agg, N_EIG, which = "LA"); rm(agg); gc(FALSE)
  
  bas <- list(values = Vr$values, vectors = Vr$vectors, r_list = r_list,
              block_matrix_list = block_list, annotation = annotation,
              sample_names = names(eig_list))
  saveRDS(bas, BASIS)
}

d0   <- find_elbow_kneedle(bas$values[1:D0_WINDOW])
samp <- bas$sample_names; S <- length(samp)
lg("d0 = ", d0, ";  n_sample = ", S,
   ";  stage-1 r range = [", min(bas$r_list), ", ", max(bas$r_list), "]")

pos  <- which(bas$values > 1e-8 * bas$values[1])
D_HI <- min(max(pos), N_EIG - 1L)
lg("largest usable dimension = ", D_HI)

VL <- bas$vectors[, 1:D_HI] %*% diag(bas$values[1:D_HI])
E_full <- lapply(seq_len(S), function(i) bas$block_matrix_list[[i]] %*% VL)
names(E_full) <- samp
F_ <- nrow(E_full[[1]]); stopifnot(F_ == length(feature_names))

VL0 <- bas$vectors[, 1:d0] %*% diag(bas$values[1:d0])
stopifnot(max(abs(bas$block_matrix_list[[1]] %*% VL0 - E_full[[1]][, 1:d0])) < 1e-8)
lg("column nesting verified (< 1e-8)")
bas$block_matrix_list <- NULL; rm(VL, VL0); gc(FALSE)

D_GRID <- sort(unique(c(1:20, 22, 25, 30, 35, 40, 50, d0)))
D_GRID <- D_GRID[D_GRID <= D_HI]
lg("D_GRID (", length(D_GRID), "): ", paste(D_GRID, collapse = ", "))

## 3. Feature stratification profiles at every d (cumulative sums)

pr <- which(upper.tri(matrix(0, S, S)), arr.ind = TRUE)
IP <- matrix(0, F_, nrow(pr)); SQ <- matrix(0, F_, S)
strat <- vector("list", length(D_GRID)); names(strat) <- as.character(D_GRID)
g <- 1L
for (k in seq_len(max(D_GRID))) {
  Ek <- vapply(E_full, function(M) M[, k], numeric(F_))
  SQ <- SQ + Ek^2
  IP <- IP + Ek[, pr[, 1]] * Ek[, pr[, 2]]
  if (g <= length(D_GRID) && k == D_GRID[g]) {
    den <- sqrt(SQ[, pr[, 1]] * SQ[, pr[, 2]]); den[den == 0] <- 1
    strat[[g]] <- IP / den
    g <- g + 1L
  }
}
rm(IP, SQ, Ek, den); gc(FALSE)

ver <- cal_peranova_per_feature_parallel_no_zero(
  lapply(E_full, function(M) M[, 1:d0, drop = FALSE]), S,
  groups = bas$annotation[samp, "Condition"], dist_method = "cosine",
  num_permutation = 9999, whether_permanova = FALSE)
ref <- do.call(rbind, lapply(ver$similarity_matrix_list,
                             function(m) m[upper.tri(m, diag = FALSE)]))
dev <- min(max(abs(ref - strat[[as.character(d0)]])),
           max(abs(ref - (1 - strat[[as.character(d0)]]))))
lg("fast stratification profiles vs reference routine at d0: max dev ",
   sprintf("%.2e", dev))
stopifnot(dev < 1e-8)
rm(ver, ref); gc(FALSE)

spearman_fast <- function(X) {
  R <- apply(X, 2, rank); s <- apply(R, 2, sd)
  bad <- !is.finite(s) | s == 0; s[bad] <- 1
  R <- sweep(sweep(R, 2, colMeans(R), "-"), 2, s, "/")
  C <- crossprod(R) / (nrow(R) - 1)
  if (any(bad)) { C[bad, ] <- 0; C[, bad] <- 0; diag(C)[bad] <- 1 }
  C
}

modules_from_gamma <- function(G) {
  DM <- 1 - G
  cutreeDynamic(dendro = HCL(as.dist(DM), method = "ward.D2"), distM = DM,
                minClusterSize = MIN_MODULE_SIZE, method = "hybrid",
                deepSplit = DEEP_SPLIT)
}

## 4. Module 27 and its 8/10 partition at d0

G0  <- spearman_fast(t(strat[[as.character(d0)]]))
cl0 <- modules_from_gamma(G0)
names(cl0) <- feature_names
lg("dynamicTreeCut at d0: ", length(unique(cl0[cl0 > 0])), " modules, ",
   sum(cl0 == 0), " unassigned")

MOD_IDX <- which(cl0 == MODULE_ID)
if (!length(MOD_IDX)) { print(table(cl0)); stop("Module ", MODULE_ID, " not found.") }
lg("Module ", MODULE_ID, ": ", length(MOD_IDX), " features")
writeLines(feature_names[MOD_IDX], file.path(OUT_DIR, "module27_features.txt"))

module_similarity <- function(idx, d) {
  X <- vapply(E_full, function(M) as.vector(M[idx, 1:d, drop = FALSE]),
              numeric(length(idx) * d))
  nrm <- sqrt(colSums(X^2)); nrm[nrm == 0] <- 1
  P <- crossprod(sweep(X, 2, nrm, "/")); dimnames(P) <- list(samp, samp); P
}
partition_of <- function(P) {                     # = find_partition_hclust()
  DM <- as.matrix(as.dist(1 - P))
  cutree(stats::hclust(dist(DM, method = "euclidean"), method = "complete"), k = 2)
}

P0   <- module_similarity(MOD_IDX, d0)
LAB0 <- partition_of(P0)
lg("Figure 5B partition at d0: sizes ",
   paste(sort(as.vector(table(LAB0))), collapse = "/"))
if (!identical(sort(as.vector(table(LAB0))), c(8L, 10L)))
  stop("Module ", MODULE_ID, " does not give the Figure 5B 8/10 split at d0.")
write.csv(data.frame(sample = samp, group = LAB0),
          file.path(OUT_DIR, "published_partition.csv"), row.names = FALSE)

MIN_GROUP <- max(3, ceiling(S * MIN_GROUP_FRAC))

## 5. Reference over all choose(18,8) bipartitions

n1   <- min(as.vector(table(LAB0)))
COMB <- utils::combn(S, n1)
MEMB <- matrix(0, S, ncol(COMB))
MEMB[cbind(as.vector(COMB), rep(seq_len(ncol(COMB)), each = n1))] <- 1
lg("enumerated ", ncol(COMB), " bipartitions of size ", n1, "/", S - n1)

sil_many <- function(DM, M) {
  nA <- colSums(M); nB <- S - nA
  s1 <- DM %*% M; s2 <- rowSums(DM) - s1
  a <- M * sweep(s1, 2, pmax(nA - 1, 1), "/") +
    (1 - M) * sweep(s2, 2, pmax(nB - 1, 1), "/")
  b <- M * sweep(s2, 2, pmax(nB, 1), "/") +
    (1 - M) * sweep(s1, 2, pmax(nA, 1), "/")
  colMeans((b - a) / pmax(a, b))
}
sil_one <- function(DM, lab) sil_many(DM, matrix(as.numeric(lab == lab[1]), ncol = 1))
stopifnot(abs(sil_one(as.matrix(as.dist(1 - P0)), LAB0) -
                mean(cluster::silhouette(LAB0, as.dist(1 - P0))[, "sil_width"])) < 1e-10)
lg("vectorised silhouette verified against cluster::silhouette")

set.seed(SEED)
NULL_SETS <- replicate(N_NULL_COH, sample.int(F_, length(MOD_IDX)), simplify = FALSE)

## 6. Sweep over d

eval_d <- function(k) {
  d <- D_GRID[k]
  G <- if (d == d0) G0 else spearman_fast(t(strat[[k]]))
  
  ## Module contrast against random sets of the same size
  contrast_of <- function(ix) {
    inn <- G[ix, ix]
    mean(inn[upper.tri(inn)]) - mean(G[ix, -ix])
  }
  delta <- contrast_of(MOD_IDX)
  nul   <- vapply(NULL_SETS, contrast_of, numeric(1))
  qd    <- quantile(nul, c(0.5, 0.95), names = FALSE)
  
  ## Module boundary under re-derivation
  cl  <- if (d == d0) cl0 else modules_from_gamma(G)
  ids <- sort(unique(cl[cl > 0]))
  if (length(ids)) {
    ov <- vapply(ids, function(m) length(intersect(which(cl == m), MOD_IDX)), numeric(1))
    sz <- vapply(ids, function(m) sum(cl == m), numeric(1))
    b  <- which.max(ov)
    containment <- ov[b] / length(MOD_IDX)
    jac <- ov[b] / (sz[b] + length(MOD_IDX) - ov[b]); best_size <- sz[b]
  } else { containment <- 0; jac <- 0; best_size <- 0L }
  rm(G); gc(FALSE)
  
  ## Donor split
  P  <- module_similarity(MOD_IDX, d)
  DM <- as.matrix(as.dist(1 - P))
  sil_pub <- sil_one(DM, LAB0)
  nulls   <- sil_many(DM, MEMB)
  qs      <- quantile(nulls, c(0.5, 0.95), names = FALSE)
  lab_d   <- partition_of(P)
  sizes   <- sort(as.vector(table(lab_d)))
  
  data.frame(
    d = d, d_over_d0 = d / d0,
    delta = delta, delta_null_med = qd[1], delta_null_q95 = qd[2],
    delta_null_max = max(nul), delta_z = (delta - mean(nul)) / sd(nul),
    max_jaccard = jac, containment = containment,
    best_module_size = best_size, n_modules = length(ids),
    sil_published = sil_pub, sil_null_med = qs[1], sil_null_q95 = qs[2],
    sil_null_max = max(nulls),
    sil_percentile = mean(nulls <= sil_pub),
    sil_rank = sum(nulls > sil_pub) + 1L,
    balanced = min(sizes) >= MIN_GROUP,
    group_sizes = paste(sizes, collapse = "/"),
    agrees_with_published = mclust::adjustedRandIndex(lab_d, LAB0) > 0.99,
    stringsAsFactors = FALSE)
}

t0  <- Sys.time()
res <- do.call(rbind, parallel::mclapply(seq_along(D_GRID), eval_d, mc.cores = NCORE))
res <- res[order(res$d), ]
lg("sweep elapsed: ", format(round(difftime(Sys.time(), t0), 1)))

stopifnot(res$max_jaccard[res$d == d0] > 0.999,
          res$agrees_with_published[res$d == d0])
write.csv(res, file.path(OUT_DIR, "hiv_module_sensitivity.csv"), row.names = FALSE)
print(res[, c("d", "delta", "delta_null_q95", "sil_published", "sil_null_max",
              "sil_rank", "agrees_with_published", "containment")],
      row.names = FALSE, digits = 3)

## 7. Spectrum and Kneedle criterion

kneedle_curve <- function(ev) {
  for (i in seq_along(ev)) if (ev[i] <= 0) ev[i] <- ev[i - 1]
  x <- seq_along(ev); y <- log(ev)
  xs <- (x - min(x)) / (max(x) - min(x)); ys <- (y - min(y)) / (max(y) - min(y))
  perpendicular_distance(xs, ys, xs[1], ys[1], xs[length(xs)], ys[length(ys)])
}
kd <- kneedle_curve(bas$values[1:D0_WINDOW])
stopifnot(max(which.max(kd), 2) == d0)
peak <- max(kd); amb <- which(kd >= 0.95 * peak)
AMB_LO <- min(amb); AMB_HI <- max(amb); runner <- order(-kd)[2]
spec <- data.frame(index = 1:D0_WINDOW,
                   log_eigenvalue = log(bas$values[1:D0_WINDOW]), kneedle = kd)
write.csv(spec, file.path(OUT_DIR, "hiv_spectrum.csv"), row.names = FALSE)

## 8. Figure

theme_s <- function(bs = 9) theme_bw(base_size = bs) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = .25, colour = "grey92"),
        panel.border = element_rect(linewidth = .4, colour = "black", fill = NA),
        axis.ticks = element_line(linewidth = .3, colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(), legend.background = element_blank(),
        legend.title = element_blank(), legend.key.size = unit(8, "pt"),
        legend.text = element_text(size = bs - 1.5),
        legend.margin = margin(0, 0, 0, 0),
        plot.title = element_text(face = "bold", size = bs + 2, hjust = 0,
                                  margin = margin(b = 3)))

band <- annotate("rect", xmin = AMB_LO, xmax = AMB_HI, ymin = -Inf, ymax = Inf,
                 fill = "grey80", alpha = .45)
vl   <- geom_vline(xintercept = d0, linetype = "dashed", linewidth = .4,
                   colour = "grey20")
xlin <- scale_x_continuous(limits = c(1, D0_WINDOW),
                           breaks = sort(unique(c(1, d0, 20, 30, 40, 50))),
                           expand = expansion(mult = .015))
xlog <- function(lo) scale_x_continuous(
  trans = "log2", limits = c(lo, max(D_GRID)),
  breaks = sort(unique(c(lo, 2, 4, 8, d0, 20, 30, 50))),
  expand = expansion(mult = .02))

## A: spectrum
pA <- ggplot(spec, aes(index, log_eigenvalue)) + band + vl +
  geom_line(linewidth = .5) + geom_point(size = .7) + xlin +
  scale_y_continuous(expand = expansion(mult = .06)) +
  labs(x = "dimension index j", y = expression(log~lambda[j]), title = "A") +
  theme_s() + theme(legend.position = "none")

## B: Kneedle criterion
labB <- data.frame(index = c(d0, runner), y = kd[c(d0, runner)],
                   lab = c(paste0("selected~(italic(d)[0]==", d0, ")"),
                           paste0("runner*'-'*up~(italic(d)==", runner, ")")),
                   hj = c(-0.2, 0.33), vj = c(0.4, 3.6))
pB <- ggplot(spec, aes(index, kneedle)) + band + vl +
  geom_line(linewidth = .5) + geom_point(size = .7) +
  geom_point(data = spec[runner, ], size = 1.7, shape = 21, fill = "white",
             colour = "black", stroke = .45) +
  geom_point(data = spec[d0, ], size = 2, shape = 21, fill = "#D55E00",
             colour = "black", stroke = .45) +
  geom_text(data = labB, aes(index, y, label = lab, hjust = hj, vjust = vj),
            parse = TRUE, size = 2.5, inherit.aes = FALSE) +
  xlin + scale_y_continuous(expand = expansion(mult = c(.18, .08))) +
  labs(x = "dimension index j", y = "normalised perpendicular distance",
       title = "B") + theme_s() + theme(legend.position = "none")

## C: donor split against all bipartitions
res$recov <- factor(ifelse(res$agrees_with_published, "recovered", "not recovered"),
                    levels = c("recovered", "not recovered"))
ylo <-min(0, min(res$sil_null_med)) - 0.02
pC <- ggplot(res, aes(d)) +
  band + vl +
  geom_ribbon(aes(ymin = sil_null_med, ymax = sil_null_q95),
              fill = "grey75", alpha = .55) +
  geom_line(aes(y = sil_null_max), colour = "grey60", linewidth = 1.1) +
  geom_line(aes(y = sil_published), colour = "grey25", linewidth = .4) +
  geom_point(aes(y = sil_published, shape = recov, fill = recov),
             size = 2.0, colour = "#D55E00", stroke = .55) +
  scale_shape_manual(values = c(21, 24), drop = FALSE) +
  scale_fill_manual(values = c("#D55E00", "white"), drop = FALSE) +
  xlog(1) +
  coord_cartesian(ylim = c(ylo, max(res$sil_published, res$sil_null_max) * 1.12)) +
  guides(shape = guide_legend(nrow = 2), fill = guide_legend(nrow = 2)) +
  labs(x = expression("enforced dimension "*italic(d)*"  (log"[2]*" scale)"),
       y = "average silhouette", title = "C") +
  theme_s() + theme(legend.position = "bottom",
                    legend.box.margin = margin(t = -4))

## D: module contrast against random 66-feature sets
resD <- res[res$d >= 2, ]
pD <- ggplot(resD, aes(d)) +
  band + vl +
  geom_ribbon(
    aes(ymin = delta_null_med, ymax = delta_null_q95),
    fill = "grey75", alpha = .55
  ) +
  geom_line(
    aes(y = delta_null_max),
    colour = "grey40", linewidth = .35,
    linetype = "22"
  ) +
  geom_line(
    aes(y = delta),
    colour = "#0072B2", linewidth = .5
  ) +
  geom_point(
    aes(y = delta),
    colour = "#0072B2", size = 1.6
  ) +
  xlog(2) +
  scale_y_continuous() +
  coord_cartesian(
    ylim = c(
      min(0, min(resD$delta_null_med, na.rm = TRUE)),
      max(resD$delta, na.rm = TRUE) * 1.15
    )
  ) +
  labs(
    x = expression("enforced dimension "*italic(d)*"  (log"[2]*" scale)"),
    y = expression("module contrast  "*Delta*"  (within " - " between)"),
    title = "D"
  ) +
  theme_s() +
  theme(legend.position = "none")

fig <- (pA | pB) / (pC | pD)

png(file.path(FIG_DIR, "figS14.png"),width = 6.9, height = 6.9, units = "in", res = 600)
print(fig)   
dev.off()

## 9. Summary

inA <- res$d >= AMB_LO & res$d <= AMB_HI
gv  <- function(col, dd) res[[col]][res$d == dd]
r2  <- res[res$d >= 2, ]

cat("\n================ HIV numbers ================\n")
cat("d0 =", d0, "; ambiguity interval [", AMB_LO, ",", AMB_HI, "]; runner-up",
    runner, sprintf("(%.2f%% below peak)\n", 100 * (peak - kd[runner]) / peak))
cat("Module", MODULE_ID, ":", length(MOD_IDX), "features; Figure 5B split 8/10\n")

cat("\n[C] DONOR SPLIT\n")
cat("  re-obtained by the pipeline's clustering rule at d =",
    paste(res$d[res$agrees_with_published], collapse = ", "),
    " [", sum(res$agrees_with_published), "of", nrow(res), "]\n")
cat("  rank among all", ncol(COMB), "bipartitions: rank 1 at",
    sum(res$sil_rank == 1), "of", nrow(res), "dimensions; worst rank",
    max(res$sil_rank), "\n")
cat("  silhouette: range", sprintf("%.3f - %.3f", min(r2$sil_published),
                                   max(r2$sil_published)), "for d >= 2\n")

cat("\n[D] MODULE CONTRAST (d >= 2)\n")
cat("  delta: range", sprintf("%.3f - %.3f", min(r2$delta), max(r2$delta)),
    "; at d0 =", sprintf("%.3f", gv("delta", d0)), "\n")
cat("  random-set 95th pct: range",
    sprintf("%.3f - %.3f", min(r2$delta_null_q95), max(r2$delta_null_q95)),
    "; random-set max: range",
    sprintf("%.3f - %.3f", min(r2$delta_null_max), max(r2$delta_null_max)), "\n")
cat("  delta exceeds random-set MAX at", sum(r2$delta > r2$delta_null_max),
    "of", nrow(r2), "dimensions\n")
cat("  z: range", sprintf("%.1f - %.1f", min(r2$delta_z), max(r2$delta_z)), "\n")

cat("\n[CSV only: module boundary under re-derivation]\n")
cat("  containment median", sprintf("%.2f", median(r2$containment)),
    "; Jaccard median", sprintf("%.2f", median(r2$max_jaccard)),
    "; best-match size", min(r2$best_module_size), "-", max(r2$best_module_size), "\n")

cat("\nNotes\n")
cat("  * Module held fixed, d varied.\n")
cat("  * Rank / percentile are stability measures, not tests;\n")
cat("    the split was defined at d0.\n")
cat("  * d = 1 excluded from D (Gamma_feat degenerates to +/-1 profiles),\n")
cat("    kept in C where P_M is built from 66-dimensional vectors.\n")

writeLines(dl, file.path(OUT_DIR, "hiv_sensitivity_diagnostics.txt"))
message("done -- ", normalizePath(FIG_DIR))
