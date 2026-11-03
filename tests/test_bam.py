import numpy as np
import pandas as pd
import pytest
from conftest import make_bam

from scbamqc.bam import (
    NO_BARCODE,
    WelfordAcc,
    collect_metrics,
    count_tags,
    filter_contigs,
    parse_soft_clips,
    read_length_from_cigar,
    to_dataframe,
)

pytest.importorskip("pysam")


def test_cigar_helpers():
    assert parse_soft_clips("5S90M5S") == (5, 5)
    assert parse_soft_clips("10S90M") == (10, 0)
    assert parse_soft_clips("90M") == (0, 0)
    assert parse_soft_clips("*") == (0, 0)
    assert read_length_from_cigar("5S40M1000N50M2I3S") == 100


def test_filter_contigs():
    c = ["chr1", "chr2", "chrX", "scaffold_1"]
    assert filter_contigs(c) == c
    assert filter_contigs(c, contig_list=["chrX", "chr1"]) == ["chr1", "chrX"]
    assert filter_contigs(c, contig_prefixes=["scaf"], contig_list=["chr2"]) == [
        "chr2",
        "scaffold_1",
    ]


def test_welford_merge_matches_numpy():
    rng = np.random.default_rng(1)
    x = rng.normal(size=1000)
    a, b = WelfordAcc(), WelfordAcc()
    for v in x[:300]:
        a.update(v)
    for v in x[300:]:
        b.update(v)
    m = a.merge(b)
    assert m.n == 1000
    assert m.mean == pytest.approx(x.mean())
    assert m.std == pytest.approx(x.std())  # population SD
    a.update(float("nan"))
    assert a.n == 300


def _reads():
    reads = []
    # cell A: 4 reads on chr1, 2 on chr2; one soft-clipped, one multimapper
    for i in range(4):
        reads.append(
            dict(
                contig="chr1",
                pos=100 + i,
                cigar="100M",
                qual=30,
                tags={"CB": "AAAA-1", "nM": i % 2, "NM": 9, "NH": 1},
            )
        )
    reads.append(
        dict(
            contig="chr2",
            pos=10,
            cigar="10S90M",
            qual=20,
            mapq=3,
            tags={"CB": "AAAA-1", "nM": 2, "NM": 9, "NH": 2},
        )
    )
    reads.append(
        dict(
            contig="chr2",
            pos=20,
            cigar="100M",
            qual=20,
            tags={"CB": "AAAA-1", "nM": 0, "NM": 9, "NH": 1},
        )
    )
    # cell B: 2 reads on chr1
    for i in range(2):
        reads.append(
            dict(
                contig="chr1",
                pos=500 + i,
                cigar="95M5S",
                qual=40,
                tags={"CB": "CCCC-1", "nM": 1, "NM": 9, "NH": 1},
            )
        )
    # no barcode
    reads.append(dict(contig="chr1", pos=900, cigar="100M", tags={"nM": 0}))
    return reads


@pytest.fixture
def bam(tmp_path):
    return make_bam(tmp_path / "test.bam", _reads())


def test_collect_metrics_per_cell(bam):
    stats = collect_metrics(str(bam), workers=2)
    df = to_dataframe(stats).set_index("cb")
    assert list(df.index) == ["AAAA-1", "CCCC-1"]

    a = df.loc["AAAA-1"]
    assert a["n_reads"] == 6
    assert a["frac_clipped"] == round(1 / 6, 4)
    assert a["frac_clip5"] == round(1 / 6, 4)
    assert a["frac_clip3"] == 0
    assert a["frac_multimapper"] == round(1 / 6, 4)
    assert a["mean_nm"] == round(np.mean([0, 1, 0, 1, 2, 0]), 3)
    assert a["mean_base_qual"] == round(np.mean([30] * 4 + [20] * 2), 3)
    assert a["mean_mapq"] == round(np.mean([255] * 5 + [3]), 3)
    assert a["mean_clip_frac"] == round(0.1 / 6, 4)

    b = df.loc["CCCC-1"]
    assert b["n_reads"] == 2
    assert b["frac_clip3"] == 1.0
    assert b["mean_clip_frac"] == 0.05
    assert b["std_clip_frac"] == 0.0


def test_collect_metrics_options(bam):
    stats = collect_metrics(str(bam), workers=1)
    df = to_dataframe(stats, include_no_cb=True)
    assert NO_BARCODE in set(df["cb"])

    by_contig = to_dataframe(stats, group_by="cell-contig")
    assert set(zip(by_contig["cb"], by_contig["contig"])) == {
        ("AAAA-1", "chr1"),
        ("AAAA-1", "chr2"),
        ("CCCC-1", "chr1"),
    }

    only2 = to_dataframe(collect_metrics(str(bam), workers=1, contig_list=["chr2"]))
    assert list(only2["cb"]) == ["AAAA-1"] and only2["n_reads"].item() == 2

    nm = to_dataframe(collect_metrics(str(bam), workers=1, nm_tag="NM"))
    assert (nm["mean_nm"] == 9).all()

    with pytest.raises(ValueError, match="No contigs remain"):
        collect_metrics(str(bam), workers=1, contig_list=["chrZ"])


def test_count_tags(bam):
    n, counts = count_tags(str(bam))
    assert n == 9
    assert counts["CB"] == 8 and counts["nM"] == 9 and counts["NM"] == 8


def test_cli_bam_writes_outputs(bam, tmp_path, monkeypatch):
    from scbamqc.cli import app

    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit) as e:
        app(["bam", str(bam), "-j", "1", "--out-prefix", "out", "--no-progress"])
    assert e.value.code in (0, None)
    df = pd.read_csv(tmp_path / "out_per_cell.csv.gz")
    assert list(df["cb"]) == ["AAAA-1", "CCCC-1"]
    assert not (tmp_path / "out_plots.png").exists()  # plots only with --plots


def test_cli_plots(bam, tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    from scbamqc.cli import app

    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        app(["bam", str(bam), "-j", "1", "--out-prefix", "o", "--plots"])
    assert (tmp_path / "o_plots.png").exists()
    assert (tmp_path / "o_per_cell.csv.gz").exists()


def test_cli_warns_on_missing_tag(bam, tmp_path, monkeypatch, capsys):
    from scbamqc.cli import app

    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        app(
            [
                "bam",
                str(bam),
                "-j",
                "1",
                "--barcode-tag",
                "XC",
                "--no-csv",
            ]
        )
    out = capsys.readouterr().out
    assert "XC" in out and "CB" in out
