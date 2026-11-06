"""Per-barcode BAM alignment QC and alignment-based cell filtering."""

__version__ = "0.1.0"

from . import pl, qc
from .qc import (
    DEFAULT_FEATURES,
    PAPER_BASE_QUAL_SCALE,
    PAPER_GMM_TOL,
    add_bam_qc,
    basic_qc,
    gmm_labels,
    qc_score,
    read_bam_qc,
    select_cells,
)

__all__ = [
    "DEFAULT_FEATURES",
    "PAPER_BASE_QUAL_SCALE",
    "PAPER_GMM_TOL",
    "add_bam_qc",
    "basic_qc",
    "gmm_labels",
    "pl",
    "qc",
    "qc_score",
    "read_bam_qc",
    "select_cells",
]
