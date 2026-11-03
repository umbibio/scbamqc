import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from conftest import synthetic_obs  # noqa: E402

import scbamqc as sbq  # noqa: E402


def _scored():
    df = synthetic_obs()
    df["qc_score"] = sbq.qc_score(df, min_reads=200)
    sbq.select_cells(df)
    return df


def test_qc_score_scatter():
    fig = sbq.pl.qc_score_scatter(_scored(), min_score=0, min_counts=100)
    assert len(fig.axes) == 5
    plt.close(fig)


def test_qc_score_hist_counts():
    df = _scored()
    fig = sbq.pl.qc_score_hist(df, threshold=0.9)
    ax = fig.axes[0]
    n_scored = int(np.isfinite(df["qc_score"]).sum())
    assert f"of {n_scored:,} scored" in ax.get_title()
    plt.close(fig)


def test_bam_qc_hist():
    df = synthetic_obs().reset_index(names="cb")
    for col in [
        "mean_mapq",
        "frac_clipped",
        "frac_clip5",
        "frac_clip3",
        "frac_multimapper",
    ]:
        df[col] = 0.1
    fig = sbq.pl.bam_qc_hist(df)
    assert len(fig.axes) == 6
    plt.close(fig)
