"""Regression against the TgRN paper, Table SI1-3 (library RH_221220.WT).

Needs <SCBAMQC_TEST_DATA>/RH-10x/RH-10x.{raw_counts.h5,bam_qc.csv.gz};
skipped otherwise.
"""

import numpy as np
import pytest
from conftest import DATA_DIR

import scbamqc as sbq

pytest.importorskip("anndata")
pytest.importorskip("h5py")

RH = DATA_DIR / "RH-10x"
pytestmark = [
    pytest.mark.paper,
    pytest.mark.skipif(not RH.exists(), reason=f"paper data not found in {RH}"),
]


@pytest.fixture(scope="module")
def rh():
    from scbamqc.io import read_10x_h5

    a = read_10x_h5(RH / "RH-10x.raw_counts.h5")
    a = sbq.add_bam_qc(a, RH / "RH-10x.bam_qc.csv.gz")
    a = sbq.basic_qc(a)
    sbq.qc_score(a, min_reads=150, base_qual_scale=sbq.PAPER_BASE_QUAL_SCALE)
    sbq.select_cells(a, min_counts=100, min_score=0)
    return a


def test_rh_gated(rh):
    assert np.isfinite(rh.obs["qc_score"]).sum() == 20_272


def test_rh_candidates(rh):
    assert rh.obs["selected"].sum() == 16_120
    assert rh.obs.loc[rh.obs["selected"], "qc_score"].min() == pytest.approx(
        0.680, abs=5e-4
    )


def test_rh_thresholds(rh):
    q = rh.obs.loc[rh.obs["selected"], "qc_score"]
    assert (q > 0.8).sum() == 13_461
    assert (q > 0.9).sum() == 8_609
    assert (q > 0.95).sum() == 4_536
