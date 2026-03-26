"""Supplemental Figure S3: memento differential correlation (run as a subprocess by s3_methods.R)
Input:    --input DIR with counts.bin (int32, row-major cells x genes) and meta.json
Output:   --output CSV with gene_1, gene_2, corr_coef, corr_se
"""
import itertools
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import memento


def run_memento_dc(counts, gene_names, sample_ids, conditions, q, test_genes,
                   num_boot, num_cpus, random_state, min_cell_count, min_perc_group,
                   n_groups_expected, pairs=None):
    X = sp.csr_matrix(counts, dtype=np.float64)
    gene_names = [str(g) for g in gene_names]
    test_genes = [str(g) for g in test_genes]
    obs = pd.DataFrame({
        "sample": [str(s) for s in sample_ids],
        "cond": np.asarray(conditions, dtype=int),
        "capture_rate": float(q),
    })
    obs.index = [f"cell{i}" for i in range(X.shape[0])]
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=gene_names))
    assert np.isfinite(adata.X.data).all() and (adata.X.data >= 0).all()
    assert adata.var_names.is_unique
    assert not adata.obs[["sample", "cond", "capture_rate"]].isna().any().any()

    memento.setup_memento(adata, q_column="capture_rate", min_cell_count=int(min_cell_count))
    memento.create_groups(adata, label_columns=["cond", "sample"])
    memento.compute_1d_moments(adata, min_perc_group=float(min_perc_group))

    group_order = list(adata.uns["memento"]["groups"])
    if len(group_order) != int(n_groups_expected):
        raise RuntimeError(f"memento retained {len(group_order)} groups, expected {n_groups_expected}")
    kept = set(adata.var_names)
    missing = [g for g in test_genes if g not in kept]
    if missing:
        raise RuntimeError(f"memento filtered {len(missing)} tested genes, e.g. {missing[:5]}")

    # all pairs of test_genes, or an explicit subset (pair chunking across workers)
    if pairs is None:
        pairs = list(itertools.combinations(test_genes, 2))
    else:
        pairs = [(str(a), str(b)) for a, b in pairs]
        bad = [p for p in pairs if p[0] not in kept or p[1] not in kept]
        if bad:
            raise RuntimeError(f"requested pairs contain filtered genes, e.g. {bad[:3]}")
    memento.compute_2d_moments(adata, pairs)

    metadata = (adata.obs.drop_duplicates("memento_group")
                .set_index("memento_group").loc[group_order, ["cond"]])
    treatment = metadata[["cond"]].astype(float)
    covariate = pd.DataFrame({"intercept": 1.0}, index=metadata.index)
    assert treatment.index.tolist() == group_order
    assert covariate.index.tolist() == group_order

    memento.ht_2d_moments(adata, treatment=treatment, covariate=covariate,
                          num_boot=int(num_boot), num_cpus=int(num_cpus),
                          random_state=int(random_state), approx="norm",
                          resample_rep=True, backend="cpu", verbose=0)
    res = memento.get_2d_ht_result(adata)
    if len(res) != len(pairs):
        raise RuntimeError("memento returned an unexpected number of pairs")
    return res[["gene_1", "gene_2", "corr_coef", "corr_se"]].reset_index(drop=True)


def main():
    """File handoff entry point: --input DIR (counts.bin + meta.json) --output CSV."""
    import argparse
    import json
    import os
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    with open(os.path.join(a.input, "meta.json")) as fh:
        meta = json.load(fh)
    n, g = int(meta["n_cells"]), int(meta["n_genes"])
    X = np.fromfile(os.path.join(a.input, "counts.bin"), dtype="<i4")
    if X.size != n * g:
        raise RuntimeError(f"counts.bin has {X.size} values, expected {n * g}")
    X = X.reshape(n, g)
    for key in ("gene_names", "sample_ids", "conditions"):
        expected = g if key == "gene_names" else n
        if len(meta[key]) != expected:
            raise RuntimeError(f"meta {key} has length {len(meta[key])}, expected {expected}")
    res = run_memento_dc(
        counts=X, gene_names=meta["gene_names"], sample_ids=meta["sample_ids"],
        conditions=meta["conditions"], q=meta["q"], test_genes=meta["test_genes"],
        num_boot=meta["num_boot"], num_cpus=meta["num_cpus"], random_state=meta["random_state"],
        min_cell_count=meta["min_cell_count"], min_perc_group=meta["min_perc_group"],
        n_groups_expected=meta["n_groups_expected"], pairs=meta.get("pairs"))
    res.to_csv(a.output, index=False, float_format="%.17g")


if __name__ == "__main__":
    main()
