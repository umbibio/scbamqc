"""Alignment-based quality score and cell selection for single-cell data.

The defaults reproduce the cell filter of the TgRN paper (SM1/SI1):

1. :func:`add_bam_qc` keeps the barcodes present in the BAM-QC table and joins
   its columns into ``obs``.
2. :func:`basic_qc` removes barcodes with fewer than 5 detected genes and genes
   detected in fewer than 5 barcodes, then records ``total_counts``,
   ``n_genes_by_counts`` and ``barcode_rank``.
3. :func:`qc_score` scores each gated barcode by the kernel density of its
   alignment profile.
4. :func:`select_cells` fits a two-component Gaussian mixture on
   ``(log10 total_counts, qc_score)`` and keeps the high-count component.
   The paper used ``tol=PAPER_GMM_TOL``; the default runs EM to convergence.

Functions take an ``AnnData`` (any object with ``.obs``, ``.X`` and
``obs_names``) or, where only ``obs`` columns are needed, a ``DataFrame``.
``anndata`` itself is not a dependency.
"""

from __future__ import annotations

from collections.abc import Sequence
from os import PathLike
from typing import Any, Literal

import numpy as np
import pandas as pd

__all__ = [
    "DEFAULT_FEATURES",
    "PAPER_BASE_QUAL_SCALE",
    "PAPER_GMM_TOL",
    "add_bam_qc",
    "basic_qc",
    "gmm_labels",
    "qc_score",
    "read_bam_qc",
    "select_cells",
]

#: Divisor applied to ``mean_base_qual`` in the TgRN paper (Phred / 10).
PAPER_BASE_QUAL_SCALE = 10.0

#: Convergence tolerance of the mixture fit in the TgRN paper (the
#: scikit-learn default). The fit stops before it converges.
PAPER_GMM_TOL = 1e-3

#: BAM-QC columns entering the kernel density, in order.
DEFAULT_FEATURES: tuple[str, ...] = (
    "mean_clip_frac",
    "std_clip_frac",
    "mean_base_qual",
)

_UNS_KEY = "scbamqc"


def _obs(data: Any) -> pd.DataFrame:
    return data.obs if hasattr(data, "obs") else data


def _record(data: Any, step: str, params: dict) -> None:
    uns = getattr(data, "uns", None)
    if uns is not None:
        uns.setdefault(_UNS_KEY, {})[step] = params


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


def read_bam_qc(path: str | PathLike) -> pd.DataFrame:
    """Read a per-barcode table written by ``scbamqc bam``, indexed by barcode."""
    return pd.read_csv(path, index_col=0)


def add_bam_qc(adata: Any, bam_qc: pd.DataFrame | str | PathLike) -> Any:
    """Restrict ``adata`` to barcodes in the BAM-QC table and join its columns.

    Barcodes keep their order in ``adata`` (the mixture fit in
    :func:`select_cells` depends on row order through its k-means
    initialization). The ``__no_barcode__`` row, if present, never matches.

    Parameters
    ----------
    adata
        ``AnnData`` with barcodes as ``obs_names``.
    bam_qc
        Per-barcode table (indexed by barcode) or a path to one.

    Returns
    -------
    A new ``AnnData`` with the matching barcodes and the BAM-QC columns in
    ``obs``.
    """
    if not isinstance(bam_qc, pd.DataFrame):
        bam_qc = read_bam_qc(bam_qc)

    keep = np.asarray(adata.obs_names.isin(bam_qc.index))
    if not keep.any():
        raise ValueError(
            "No barcode of adata is in the BAM-QC table. Check that both use the "
            "same barcode format (e.g. the '-1' suffix)."
        )
    overlap = adata.obs.columns.intersection(bam_qc.columns)
    if len(overlap):
        raise ValueError(f"Columns already in adata.obs: {list(overlap)}")

    out = adata[keep].copy()
    out.obs = out.obs.join(bam_qc)
    return out


def basic_qc(
    adata: Any,
    min_genes: int = 5,
    min_cells: int = 5,
    layer: str | None = None,
) -> Any:
    """Remove near-empty barcodes and rare genes and compute count metrics.

    Barcodes with fewer than ``min_genes`` detected genes are removed first,
    then genes detected in fewer than ``min_cells`` of the remaining barcodes.
    The metrics are computed after both filters, matching
    ``sclab.preprocess.qc(adata, min_counts=0)``.

    Adds to ``obs``: ``n_genes_by_counts``, ``total_counts``, ``barcode_rank``
    (rank by descending ``total_counts``, ties averaged). Adds to ``var``:
    ``n_cells_by_counts``, ``total_counts``.

    Returns
    -------
    A new, filtered ``AnnData``.
    """
    X = adata.layers[layer] if layer is not None else adata.X

    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    cell_mask = n_genes >= min_genes
    n_cells = np.asarray((X[cell_mask] > 0).sum(axis=0)).ravel()
    gene_mask = n_cells >= min_cells

    out = adata[cell_mask][:, gene_mask].copy()
    X = out.layers[layer] if layer is not None else out.X

    out.obs["n_genes_by_counts"] = np.asarray((X > 0).sum(axis=1)).ravel()
    out.obs["total_counts"] = np.asarray(X.sum(axis=1)).ravel()
    out.obs["barcode_rank"] = out.obs["total_counts"].rank(ascending=False)
    out.var["n_cells_by_counts"] = np.asarray((X > 0).sum(axis=0)).ravel()
    out.var["total_counts"] = np.asarray(X.sum(axis=0)).ravel()

    _record(out, "basic_qc", {"min_genes": min_genes, "min_cells": min_cells})
    return out


# ---------------------------------------------------------------------------
# Quality score
# ---------------------------------------------------------------------------


def qc_score(
    data: Any,
    min_reads: float = 100,
    max_mean_nm: float = 1.0,
    base_qual_scale: float | Literal["median"] = "median",
    features: Sequence[str] = DEFAULT_FEATURES,
    key_added: str | None = "qc_score",
) -> np.ndarray:
    """Kernel-density quality score of each barcode's alignment profile.

    Barcodes pass the read gate when ``n_reads > min_reads`` and
    ``mean_nm < max_mean_nm``. A Gaussian kernel density (``scipy.stats.
    gaussian_kde``, Scott's rule) is fitted on the gated barcodes' ``features``
    and evaluated at each of them. The score is the min-max scaled
    ``log10(density + 1)``, so it lies in [0, 1]. Barcodes outside the gate get
    ``NaN``.

    Parameters
    ----------
    data
        ``AnnData`` or ``DataFrame`` with ``n_reads``, ``mean_nm`` and the
        ``features`` columns.
    min_reads
        Read gate (strict ``>``). Chosen per library in the paper (50-3,000).
    max_mean_nm
        Mismatch gate (strict ``<``).
    base_qual_scale
        ``mean_base_qual`` is divided by this value before density estimation.
        ``"median"`` (default) uses the median ``mean_base_qual`` of the gated
        barcodes. To reproduce the paper, use :data:`PAPER_BASE_QUAL_SCALE`.
        The kernel density is affine-invariant, so this scale does not change
        the ranking of barcodes. It multiplies every density by the same
        constant and so changes how the ``+ 1`` in ``log10(density + 1)``
        compresses scores, which moves thresholds such as ``qc_score > 0.9``.
    features
        BAM-QC columns entering the density.
    key_added
        If ``data`` is an ``AnnData`` and this is not ``None``, the score is
        stored in ``obs[key_added]`` and the parameters in
        ``uns["scbamqc"]["qc_score"]``.

    Returns
    -------
    Array of scores aligned with the rows of ``data``.
    """
    from scipy.stats import gaussian_kde

    obs = _obs(data)
    mask = obs["n_reads"].values > min_reads
    mask &= obs["mean_nm"].values < max_mean_nm
    n = int(mask.sum())
    if n <= len(features):
        raise ValueError(
            f"Only {n} barcodes pass the read gate (n_reads > {min_reads}, "
            f"mean_nm < {max_mean_nm}); too few to estimate a density."
        )

    if isinstance(base_qual_scale, str):
        if base_qual_scale != "median":
            raise ValueError(
                f"base_qual_scale must be a number or 'median', got {base_qual_scale!r}"
            )
        base_qual_scale = float(np.median(obs["mean_base_qual"].values[mask]))
    if not base_qual_scale > 0:
        raise ValueError(f"base_qual_scale must be positive, got {base_qual_scale}")

    X = np.vstack(
        [
            obs[f].values / base_qual_scale if f == "mean_base_qual" else obs[f].values
            for f in features
        ]
    )

    kde = gaussian_kde(X[:, mask])
    out = np.full(X.shape[1], np.nan)
    out[mask] = kde(X[:, mask])

    out = np.log10(out + 1)
    out[mask] -= out[mask].min()
    out[mask] /= out[mask].max()

    if key_added is not None and hasattr(data, "obs"):
        data.obs[key_added] = out
        _record(
            data,
            "qc_score",
            {
                "min_reads": min_reads,
                "max_mean_nm": max_mean_nm,
                "base_qual_scale": base_qual_scale,
                "features": list(features),
                "n_gated": n,
            },
        )
    return out


# ---------------------------------------------------------------------------
# Mixture model and selection
# ---------------------------------------------------------------------------


def gmm_labels(
    data: Any,
    n_components: int = 2,
    random_state: int = 0,
    tol: float = 1e-10,
    max_iter: int = 100_000,
    counts_key: str = "total_counts",
    score_key: str = "qc_score",
    key_added: str | None = "gmm_label",
) -> np.ndarray:
    """Gaussian-mixture labels on ``(log10 total_counts, qc_score)``.

    The mixture is fitted on barcodes with ``qc_score > 0``. Components are
    ranked by mean ``log10`` count, the highest-count component labeled 0.
    Barcodes not in the fit are labeled -1.

    EM runs until the mean log-likelihood changes by less than ``tol``. The
    default runs to convergence, so the result does not depend on the
    initialization. Use ``tol=PAPER_GMM_TOL`` to repeat the paper, which
    stopped early.
    """
    from sklearn.mixture import GaussianMixture

    obs = _obs(data)
    counts = obs[counts_key].values
    score = obs[score_key].values

    log_counts = np.log10(np.clip(counts, 1, None))
    mask = score > 0
    X = np.column_stack([log_counts[mask], score[mask]])

    gmm = GaussianMixture(
        n_components=n_components,
        random_state=random_state,
        tol=tol,
        max_iter=max_iter,
    )
    raw_labels = gmm.fit(X).predict(X)

    mean_counts = [
        log_counts[mask][raw_labels == k].mean() for k in range(n_components)
    ]
    rank = np.argsort(np.argsort(-np.array(mean_counts)))
    labels = np.full(len(score), -1)
    labels[mask] = rank[raw_labels]

    if key_added is not None:
        obs[key_added] = labels
    return labels


def select_cells(
    data: Any,
    min_counts: float = 100,
    min_score: float = 0.0,
    n_components: int = 2,
    random_state: int = 0,
    tol: float = 1e-10,
    max_iter: int = 100_000,
    counts_key: str = "total_counts",
    score_key: str = "qc_score",
    key_added: str | None = "selected",
) -> np.ndarray:
    """Select candidate cells: high-count mixture component above both floors.

    A barcode is selected when ``total_counts > min_counts``,
    ``qc_score > min_score`` and its :func:`gmm_labels` label is 0. Also
    stores the mixture labels in ``obs["gmm_label"]``. See
    :func:`gmm_labels` for ``tol`` and ``max_iter``.

    Returns
    -------
    Boolean mask aligned with the rows of ``data``.
    """
    obs = _obs(data)
    labels = gmm_labels(
        data,
        n_components=n_components,
        random_state=random_state,
        tol=tol,
        max_iter=max_iter,
        counts_key=counts_key,
        score_key=score_key,
    )
    selected = (
        (obs[counts_key].values > min_counts)
        & (obs[score_key].values > min_score)
        & (labels == 0)
    )
    if key_added is not None:
        obs[key_added] = selected
        _record(
            data,
            "select_cells",
            {
                "min_counts": min_counts,
                "min_score": min_score,
                "n_components": n_components,
                "random_state": random_state,
                "tol": tol,
                "max_iter": max_iter,
                "n_selected": int(selected.sum()),
            },
        )
    return selected
