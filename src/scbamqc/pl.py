"""Plots of BAM-QC metrics and the quality score.

All functions return a ``matplotlib`` ``Figure`` and never call ``plt.show()``.
matplotlib is optional: ``pip install scbamqc[plot]``.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

__all__ = ["bam_qc_hist", "qc_score_hist", "qc_score_scatter"]

# Validated categorical slots (blue, orange) plus a neutral for "not kept".
_BLUE = "#2a78d6"
_ORANGE = "#eb6834"
_GRAY = "#9a9996"
_INK = "#52514e"


def _pyplot():
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:  # pragma: no cover - depends on the environment
        raise ImportError(
            "Plotting requires matplotlib. Install it with `pip install scbamqc[plot]`."
        ) from e
    return plt


def _obs(data: Any) -> pd.DataFrame:
    return data.obs if hasattr(data, "obs") else data


def _style(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(True, axis="y", color="#e4e3df", linewidth=0.6)
    ax.set_axisbelow(True)


def _id_label(barcode_tag: str) -> str:
    if barcode_tag == "CB":
        return "cell barcodes"
    if barcode_tag == "RG":
        return "cell IDs"
    return f"cell IDs ({barcode_tag})"


def bam_qc_hist(
    df: pd.DataFrame,
    barcode_tag: str = "CB",
    group_by: str = "cell",
    bins: int = 60,
):
    """Six-panel histogram of a BAM-QC table (the ``scbamqc bam`` plot)."""
    plt = _pyplot()

    row_label = "cell×contig" if group_by == "cell-contig" else "cell"
    fig, axes = plt.subplots(2, 3, figsize=(14, 9), layout="constrained")
    ax0, ax1, ax2, ax3, ax4, ax5 = axes.ravel()

    panels = [
        (ax0, "mean_mapq", f"Mean MAPQ per {row_label}", "Mean MAPQ"),
        (ax1, "mean_base_qual", f"Mean base quality per {row_label}", "Mean Phred"),
        (ax2, "frac_clipped", "Fraction of clipped reads", "Fraction clipped"),
        (ax4, "frac_multimapper", "Multimapper fraction", "Fraction NH > 1"),
        (ax5, "mean_nm", "Mean mismatches", "Mean nM"),
    ]
    for ax, col, title, xlabel in panels:
        ax.hist(df[col].dropna(), bins=bins, color=_BLUE, edgecolor="none")
        ax.set(title=title, xlabel=xlabel, ylabel="Rows")

    for col, color, label in [
        ("frac_clip5", _BLUE, "5′"),
        ("frac_clip3", _ORANGE, "3′"),
    ]:
        ax3.hist(
            df[col].dropna(),
            bins=bins,
            color=color,
            alpha=0.6,
            label=label,
            edgecolor="none",
        )
    ax3.set(title="5′ vs 3′ clipped fraction", xlabel="Fraction", ylabel="Rows")
    ax3.legend(frameon=False)

    for ax in axes.ravel():
        _style(ax)

    n_cells = df["cb"].nunique() if "cb" in df else len(df)
    fig.suptitle(
        f"BAM QC: {n_cells:,} {_id_label(barcode_tag)}, per-{row_label} distributions",
        fontsize=13,
    )
    return fig


def qc_score_scatter(
    data: Any,
    x: Sequence[str] = (
        "total_counts",
        "mean_clip_frac",
        "std_clip_frac",
        "mean_base_qual",
        "mean_nm",
    ),
    score_key: str = "qc_score",
    selected_key: str | None = "selected",
    min_score: float | None = None,
    min_counts: float | None = None,
    log_x: Sequence[str] = ("total_counts",),
    s: float = 1,
    alpha: float = 0.05,
    rasterized: bool = True,
):
    """``qc_score`` against BAM-QC and count metrics, one panel per metric.

    Barcodes in ``obs[selected_key]`` are drawn in blue, the rest in gray.
    Optional dotted lines mark ``min_score`` and ``min_counts``.
    """
    plt = _pyplot()

    obs = _obs(data)
    y = obs[score_key].values
    if selected_key is not None and selected_key in obs:
        sel = obs[selected_key].values.astype(bool)
    else:
        sel = np.ones(len(obs), dtype=bool)

    fig, axes = plt.subplots(
        1,
        len(x),
        figsize=(3.6 * len(x), 4),
        sharey=True,
        layout="constrained",
        squeeze=False,
    )
    axes = axes.ravel()
    # Draw the unselected first so the selected stay visible on top.
    order = [(~sel, _GRAY), (sel, _BLUE)]
    for ax, col in zip(axes, x):
        xv = obs[col].values.astype(float)
        if col in log_x:
            xv = np.clip(xv, 1, None)
            ax.set_xscale("log")
        for m, color in order:
            ax.scatter(
                xv[m],
                y[m],
                s=s,
                alpha=alpha,
                c=color,
                linewidths=0,
                rasterized=rasterized,
            )
        ax.set_xlabel(col)
        if min_score is not None:
            ax.axhline(min_score, color=_INK, linewidth=1, linestyle=":")
        if min_counts is not None and col == "total_counts":
            ax.axvline(min_counts, color=_INK, linewidth=1, linestyle=":")
        _style(ax)
    axes[0].set_ylabel(score_key)

    title = f"N scored = {int(np.isfinite(y).sum()):,}"
    if selected_key is not None and selected_key in obs:
        title += f", selected (blue) = {int(sel.sum()):,}"
    fig.suptitle(title)
    return fig


def qc_score_hist(
    data: Any,
    threshold: float = 0.9,
    score_key: str = "qc_score",
    selected_key: str = "selected",
    bins: int = 50,
    title: str | None = None,
    ax=None,
):
    """Stacked ``qc_score`` histogram of scored barcodes, by filter outcome.

    Gray: not selected by :func:`scbamqc.qc.select_cells`. Orange: candidate
    cell with ``qc_score <= threshold``. Blue: retained
    (``qc_score > threshold``). The dashed line marks ``threshold``.
    """
    plt = _pyplot()

    obs = _obs(data)
    y = obs[score_key].values
    scored = np.isfinite(y)
    sel = obs[selected_key].values.astype(bool) if selected_key in obs else scored
    retained = scored & sel & (y > threshold)
    candidate = scored & sel & ~retained
    background = scored & ~sel

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 3.5), layout="constrained")
    else:
        fig = ax.figure

    edges = np.linspace(0, 1, bins + 1)
    ax.hist(
        [y[background], y[candidate], y[retained]],
        bins=edges,
        stacked=True,
        color=[_GRAY, _ORANGE, _BLUE],
        label=[
            f"not selected ({background.sum():,})",
            f"candidate, ≤ {threshold:g} ({candidate.sum():,})",
            f"retained, > {threshold:g} ({retained.sum():,})",
        ],
        edgecolor="white",
        linewidth=0.5,
    )
    ax.axvline(threshold, color=_INK, linestyle="--", linewidth=1)
    ax.set(xlabel=score_key, ylabel="Barcodes", xlim=(0, 1))
    ax.legend(frameon=False, loc="upper left", fontsize=8)
    if title is None:
        title = f"{retained.sum():,} retained of {scored.sum():,} scored"
    ax.set_title(title)
    _style(ax)
    return fig
