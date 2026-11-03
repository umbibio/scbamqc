"""Readers. Need the optional ``anndata`` and ``h5py`` packages."""

from __future__ import annotations

from os import PathLike

__all__ = ["read_10x_h5"]


def read_10x_h5(path: str | PathLike, gene_expression_only: bool = True):
    """Read a Cell Ranger v3+ feature-barcode matrix (``.h5``) as ``AnnData``.

    Equivalent to ``scanpy.read_10x_h5`` for the v3 format: ``X`` is a
    ``float32`` CSR matrix (barcodes x features), ``var_names`` are gene names
    (made unique), and ``var`` holds ``gene_ids``, ``feature_types`` and
    ``genome``.
    """
    try:
        import anndata as ad
        import h5py
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "read_10x_h5 requires anndata and h5py: `pip install anndata h5py`."
        ) from e
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    with h5py.File(path, "r") as f:
        if "matrix" not in f:
            raise ValueError(f"{path} is not a Cell Ranger v3+ matrix (no 'matrix').")
        m = f["matrix"]
        n_features, n_barcodes = m["shape"][()]
        X = csr_matrix(
            (m["data"][()].astype(np.float32), m["indices"][()], m["indptr"][()]),
            shape=(n_barcodes, n_features),
        )
        feats = m["features"]

        def _str(ds):
            return ds[()].astype(str)

        var = pd.DataFrame(
            {
                "gene_ids": _str(feats["id"]),
                "feature_types": _str(feats["feature_type"]),
                "genome": _str(feats["genome"]),
            },
            index=pd.Index(_str(feats["name"])),
        )
        obs = pd.DataFrame(index=pd.Index(_str(m["barcodes"])))

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names_make_unique()
    if gene_expression_only:
        adata = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
    return adata
