import numpy as np
import pandas as pd
import pytest
from conftest import synthetic_obs

import scbamqc as sbq

ad = pytest.importorskip("anndata")


def _adata(n_cells=300, n_genes=50, seed=0):
    from scipy.sparse import csr_matrix

    rng = np.random.default_rng(seed)
    X = rng.poisson(0.3, (n_cells, n_genes)).astype(np.float32)
    X[:20] = 0  # empty barcodes
    X[:, :3] = 0  # undetected genes
    X[20:25, 3] = 1  # gene detected in exactly 5 cells
    X[:, 3][25:] = 0
    return ad.AnnData(
        X=csr_matrix(X),
        obs=pd.DataFrame(index=[f"C{i:04d}-1" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]),
    )


def test_basic_qc():
    a = _adata()
    dense = a.X.toarray()
    out = sbq.basic_qc(a)
    cell_ok = (dense > 0).sum(1) >= 5
    gene_ok = (dense[cell_ok] > 0).sum(0) >= 5
    assert out.shape == (cell_ok.sum(), gene_ok.sum())
    kept = dense[cell_ok][:, gene_ok]
    np.testing.assert_allclose(out.obs["total_counts"], kept.sum(1))
    np.testing.assert_array_equal(out.obs["n_genes_by_counts"], (kept > 0).sum(1))
    assert out.obs["barcode_rank"].min() >= 1
    assert out.uns["scbamqc"]["basic_qc"] == {"min_genes": 5, "min_cells": 5}
    assert a.shape == (300, 50)  # input untouched


def test_add_bam_qc_keeps_adata_order():
    a = _adata()
    table = pd.DataFrame(
        {"n_reads": np.arange(100)}, index=list(a.obs_names[::3][::-1][:100])
    )
    out = sbq.add_bam_qc(a, table)
    assert list(out.obs_names) == [b for b in a.obs_names if b in set(table.index)]
    np.testing.assert_array_equal(
        out.obs["n_reads"].values, table.loc[out.obs_names, "n_reads"].values
    )


def test_add_bam_qc_no_overlap():
    a = _adata()
    table = pd.DataFrame({"n_reads": [1]}, index=["AAAA"])
    with pytest.raises(ValueError, match="barcode format"):
        sbq.add_bam_qc(a, table)


def test_qc_score_on_anndata_records_params():
    obs = synthetic_obs(n=400)
    a = ad.AnnData(X=np.zeros((len(obs), 2), dtype=np.float32), obs=obs)
    q = sbq.qc_score(a, min_reads=200)
    np.testing.assert_array_equal(a.obs["qc_score"].values, q)
    rec = a.uns["scbamqc"]["qc_score"]
    gated = (obs["n_reads"] > 200) & (obs["mean_nm"] < 1)
    assert rec["base_qual_scale"] == obs.loc[gated, "mean_base_qual"].median()
    assert rec["n_gated"] == int(np.isfinite(q).sum())


def test_base_qual_scale_median():
    df = synthetic_obs()
    gated = (df["n_reads"] > 200) & (df["mean_nm"] < 1)
    med = df.loc[gated, "mean_base_qual"].median()
    np.testing.assert_array_equal(
        sbq.qc_score(df, min_reads=200),
        sbq.qc_score(df, min_reads=200, base_qual_scale=med),
    )
    with pytest.raises(ValueError, match="median"):
        sbq.qc_score(df, min_reads=200, base_qual_scale="mean")
