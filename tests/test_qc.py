import numpy as np
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
