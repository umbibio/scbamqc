import numpy as np
import pandas as pd
import pytest
from conftest import synthetic_obs
from scipy.stats import spearmanr

import scbamqc as sbq


def reference_bam_qc_score(df, min_reads=100):
    """Verbatim copy of `bam_qc_score` from 00_prefilter10.ipynb (cell 1)."""
    from scipy.stats import gaussian_kde

    X = np.vstack(
        [
            df["mean_clip_frac"].values,
            df["std_clip_frac"].values,
            df["mean_base_qual"].values / 10,
        ]
    )
    mask = df["n_reads"].values > min_reads
    mask &= df["mean_nm"].values < 1
    kde = gaussian_kde(X[:, mask])
    out = np.full(X.shape[1], np.nan)
    out[mask] = kde(X[:, mask])
    out = np.log10(out + 1)
    out[mask] -= out[mask].min()
    out[mask] /= out[mask].max()
    return out


def test_qc_score_matches_notebook():
    df = synthetic_obs()
    expected = reference_bam_qc_score(df, min_reads=200)
    got = sbq.qc_score(df, min_reads=200, base_qual_scale=sbq.PAPER_BASE_QUAL_SCALE)
    np.testing.assert_array_equal(np.isnan(got), np.isnan(expected))
    np.testing.assert_array_equal(got[~np.isnan(got)], expected[~np.isnan(expected)])


def test_qc_score_gate_and_range():
    df = synthetic_obs()
    q = sbq.qc_score(df, min_reads=200)
    gated = (df["n_reads"] > 200) & (df["mean_nm"] < 1)
    assert np.isfinite(q[gated]).all()
    assert np.isnan(q[~gated]).all()
    assert np.nanmin(q) == 0 and np.nanmax(q) == 1

    # gate is strict: a barcode at exactly min_reads is not scored
    df.iloc[0, df.columns.get_loc("n_reads")] = 200
    assert np.isnan(sbq.qc_score(df, min_reads=200)[0])


def test_qc_score_good_population_scores_higher():
    df = synthetic_obs()
    q = sbq.qc_score(df, min_reads=200)
    n_good = len(df) // 2
    assert np.nanmedian(q[:n_good]) > np.nanmedian(q[n_good:])


def test_base_qual_scale_preserves_ranking():
    """The KDE is affine-invariant: the scale changes values, not the order."""
    df = synthetic_obs()
    q10 = sbq.qc_score(df, min_reads=200, base_qual_scale=10.0)
    q1 = sbq.qc_score(df, min_reads=200, base_qual_scale=1.0)
    ok = np.isfinite(q10)
    assert spearmanr(q10[ok], q1[ok]).statistic == pytest.approx(1.0)
    assert not np.allclose(q10[ok], q1[ok])


def test_qc_score_too_few_barcodes():
    df = synthetic_obs(n=20)
    with pytest.raises(ValueError, match="too few"):
        sbq.qc_score(df, min_reads=10**9)


def test_select_cells_and_labels():
    df = synthetic_obs()
    df["qc_score"] = sbq.qc_score(df, min_reads=200)
    sel = sbq.select_cells(df, min_counts=100)
    labels = df["gmm_label"].values
    assert set(np.unique(labels)) <= {-1, 0, 1}
    assert (labels[~(df["qc_score"].values > 0)] == -1).all()
    # label 0 is the higher-count component
    lc = np.log10(df["total_counts"].clip(lower=1))
    assert lc[labels == 0].mean() > lc[labels == 1].mean()
    assert df["selected"].dtype == bool
    np.testing.assert_array_equal(
        sel, (labels == 0) & (df["total_counts"] > 100) & (df["qc_score"] > 0)
    )


# ---------------------------------------------------------------- AnnData helpers

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
