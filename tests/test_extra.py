"""Tests of the opt-in extra metrics (``--extra-metrics``)."""

import gzip
from io import StringIO

import numpy as np
import pandas as pd
import pytest
from conftest import (
    GENOMIC,
    GOLDEN_DIR,
    GOLDEN_RUNS,
    GTF_CONTIGS,
    TINY_GTF,
    TOXODB_GTF,
    TSO,
    cellranger_reads,
    clip_reads,
    golden_fixtures,
    strand_clip_reads,
    gtf_reads,
    make_bam,
    revcomp,
)

from scbamqc.annot import (
    ANTISENSE,
    EXONIC,
    INTERGENIC,
    INTRONIC,
    ContigAnnot,
    load_gtf,
)
from scbamqc.bam import (
    EXTRA_COLUMNS,
    NO_BARCODE,
    CellStats,
    ExtraStats,
    UmiCounter,
    collect_metrics,
    scan_clips,
    to_dataframe,
)

pytest.importorskip("pysam")

DEFAULT_COLUMNS = list(CellStats().to_dict("x"))


def _run_cli(bam, tmp_path, monkeypatch, *args) -> str:
    from scbamqc.cli import app

    monkeypatch.chdir(tmp_path)
    cmd = ["bam", str(bam), "-j", "2", "--out-prefix", "o", "--no-progress"]
    with pytest.raises(SystemExit) as e:
        app([*cmd, "--no-check-tags", *args])
    assert e.value.code in (0, None)
    name = (
        "o_per_cell_extra.csv.gz" if "--extra-metrics" in args else "o_per_cell.csv.gz"
    )
    return gzip.open(tmp_path / name, "rt").read()


def test_cli_output_names(cr_bam, tmp_path, monkeypatch):
    """--extra-metrics writes *_per_cell_extra.csv.gz; the default name is kept."""
    _run_cli(cr_bam, tmp_path, monkeypatch)
    assert sorted(p.name for p in tmp_path.glob("o_*")) == ["o_per_cell.csv.gz"]
    _run_cli(cr_bam, tmp_path, monkeypatch, "--extra-metrics")
    assert sorted(p.name for p in tmp_path.glob("o_*")) == [
        "o_per_cell.csv.gz",
        "o_per_cell_extra.csv.gz",
    ]


# ---------------------------------------------------------------------------
# Default output
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("run", list(GOLDEN_RUNS))
@pytest.mark.parametrize("fixture", list(golden_fixtures()))
def test_default_output_unchanged(fixture, run, tmp_path, monkeypatch):
    """Without --extra-metrics, the CSV text equals the output of 0.4.0.dev0."""
    reads, contigs = golden_fixtures()[fixture]
    bam = make_bam(tmp_path / "t.bam", reads, contigs=contigs)
    text = _run_cli(bam, tmp_path, monkeypatch, *GOLDEN_RUNS[run])
    assert text == (GOLDEN_DIR / f"{fixture}__{run}.csv").read_text()


@pytest.mark.parametrize("run", list(GOLDEN_RUNS))
@pytest.mark.parametrize("fixture", list(golden_fixtures()))
def test_extra_metrics_only_append_columns(fixture, run, tmp_path, monkeypatch):
    """With --extra-metrics, the default columns keep their order and values."""
    reads, contigs = golden_fixtures()[fixture]
    bam = make_bam(tmp_path / "t.bam", reads, contigs=contigs)
    text = _run_cli(bam, tmp_path, monkeypatch, *GOLDEN_RUNS[run], "--extra-metrics")
    golden = pd.read_csv(GOLDEN_DIR / f"{fixture}__{run}.csv")
    extra = pd.read_csv(StringIO(text))
    n = golden.shape[1]
    assert list(extra.columns[n:]) == list(EXTRA_COLUMNS)
    pd.testing.assert_frame_equal(extra.iloc[:, :n], golden)


def test_default_dataframe_has_no_extra_columns(tmp_path):
    bam = make_bam(tmp_path / "t.bam", cellranger_reads())
    stats = collect_metrics(str(bam), workers=1)
    assert all(cs.extra is None for cs in stats.values())
    assert list(to_dataframe(stats).columns) == DEFAULT_COLUMNS
    assert list(to_dataframe({}).columns) == DEFAULT_COLUMNS


# ---------------------------------------------------------------------------
# Cell Ranger tags: regions, UMIs, pa/ts, MAPQ and nM distributions
# ---------------------------------------------------------------------------


@pytest.fixture
def cr_bam(tmp_path):
    return make_bam(tmp_path / "cr.bam", cellranger_reads())


@pytest.mark.parametrize("workers", [1, 2])
def test_tags_path(cr_bam, workers):
    stats = collect_metrics(str(cr_bam), workers=workers, extra_metrics=True)
    df = to_dataframe(stats).set_index("cb")
    a = df.loc["AAAA-1"]
    assert a["n_reads"] == 7
    assert a["region_source"] == "tags"
    # 6 reads have RE: E, E, N, I, antisense (AN only), E (AN + GX)
    assert a["frac_exonic"] == round(3 / 6, 4)
    assert a["frac_intronic"] == round(1 / 6, 4)
    assert a["frac_intergenic"] == round(1 / 6, 4)
    assert a["frac_antisense"] == round(1 / 6, 4)
    # chr1: U1 x2, U2, U3; chr2: U1, U4 (the secondary U4 read is not counted)
    assert a["n_umi"] == 5
    assert a["mean_reads_per_umi"] == round(6 / 5, 3)
    assert a["frac_singleton_umi"] == round(4 / 5, 4)
    assert a["n_umi_counted"] == 3  # xf & 8
    assert a["polya_tso_source"] == "tags"
    assert a["frac_polya"] == round(1 / 7, 4)
    assert a["frac_tso"] == round(1 / 7, 4)
    assert a["frac_mapq_unique"] == round(5 / 7, 4)
    assert a["frac_nm0"] == round(4 / 7, 4)
    assert a["frac_nm_ge2"] == round(2 / 7, 4)

    c = df.loc["CCCC-1"]
    assert c["frac_exonic"] == 0.5 and c["frac_intergenic"] == 0.5
    assert c["n_umi"] == 0  # the BAM has UMIs, this barcode has none
    assert c["n_umi_counted"] == 1
    assert np.isnan(c["mean_reads_per_umi"])
    assert np.isnan(c["frac_singleton_umi"])
    assert c["frac_polya"] == 0 and c["frac_tso"] == 0


def test_umi_cross_contig(cr_bam):
    """A UMI on two contigs counts once per contig, and the sums merge."""
    stats = collect_metrics(str(cr_bam), workers=2, extra_metrics=True)
    by_contig = to_dataframe(stats, group_by="cell-contig")
    a = by_contig[by_contig["cb"] == "AAAA-1"].set_index("contig")
    assert a.loc["chr1", "n_umi"] == 3
    assert a.loc["chr1", "mean_reads_per_umi"] == round(4 / 3, 3)
    assert a.loc["chr1", "frac_singleton_umi"] == round(2 / 3, 4)
    assert a.loc["chr2", "n_umi"] == 2
    assert a.loc["chr2", "frac_singleton_umi"] == 1.0
    assert a.loc["chr1", "n_umi_counted"] == 2
    assert a.loc["chr2", "n_umi_counted"] == 1
    # contig columns stay before the extra columns
    assert list(by_contig.columns) == [*DEFAULT_COLUMNS, "contig", *EXTRA_COLUMNS]

    # U1 is on chr1 and chr2: 5 UMIs per barcode, not 4
    per_cell = to_dataframe(stats).set_index("cb")
    assert per_cell.loc["AAAA-1", "n_umi"] == 5


def test_extra_stats_merge():
    a = ExtraStats(region_source="gtf", region=[1, 2, 3, 4], n_umi=2, n_nm0=1)
    a.add_umis(2, 4, 1)
    b = ExtraStats(region_source="gtf", region=[1, 0, 0, 1], n_tso=1, n_xf=2)
    b.add_umis(1, 1, 1)
    m = a.merge(b)
    assert m.region == [2, 2, 3, 5]
    assert (m.n_umi, m.n_umi_reads, m.n_umi_singleton) == (5, 5, 2)
    assert m.n_nm0 == 1 and m.n_tso == 1 and m.n_xf == 2
    cs = CellStats(n_reads=3, extra=a).merge(CellStats(n_reads=1, extra=b))
    assert cs.extra.region == [2, 2, 3, 5]
    assert CellStats().merge(CellStats()).extra is None


def test_mapq_unique_threshold(cr_bam):
    stats = collect_metrics(str(cr_bam), workers=1, extra_metrics=True, mapq_unique=3)
    a = to_dataframe(stats).set_index("cb").loc["AAAA-1"]
    assert a["frac_mapq_unique"] == round(6 / 7, 4)


def test_no_barcode_row_has_extra_metrics(cr_bam):
    stats = collect_metrics(str(cr_bam), workers=1, extra_metrics=True)
    df = to_dataframe(stats, include_no_cb=True).set_index("cb")
    assert df.loc[NO_BARCODE, "n_umi"] == 1
    assert df.loc[NO_BARCODE, "frac_exonic"] == 1.0


def test_tags_win_over_gtf(cr_bam):
    stats = collect_metrics(
        str(cr_bam), workers=1, extra_metrics=True, gtf=str(TINY_GTF)
    )
    df = to_dataframe(stats)
    assert (df["region_source"] == "tags").all()


# ---------------------------------------------------------------------------
# GTF regions
# ---------------------------------------------------------------------------


def test_load_gtf_tiny():
    annot = load_gtf(TINY_GTF)
    assert set(annot) == {"chr1", "chr2"}
    chr1 = annot["chr1"]
    assert chr1.exons["+"] == ([1000, 1800, 2500], [1500, 1900, 3000])
    assert chr1.genes["+"] == ([1000], [3000])
    assert chr1.exons["-"] == ([5000], [6000])
    # G5 has strand ".": on both strands
    assert annot["chr2"].genes["-"] == ([1000, 5000], [2000, 5500])
    assert annot["chr2"].genes["+"] == ([5000], [5500])


def test_load_gtf_without_gene_rows(tmp_path):
    """Without gene and transcript rows, each gene_id spans its exons."""
    gtf = tmp_path / "exons.gtf.gz"
    lines = [
        'chr1\tx\texon\t101\t200\t.\t+\t.\tgene_id "A"; transcript_id "A1";',
        'chr1\tx\texon\t401\t500\t.\t+\t.\tgene_id "A"; transcript_id "A1";',
        'chr1\tx\texon\t901\t950\t.\t.\t.\tgene_id "B"; transcript_id "B1";',
    ]
    with gzip.open(gtf, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    chr1 = load_gtf(gtf)["chr1"]
    assert chr1.genes["+"] == ([100, 900], [500, 950])
    assert chr1.genes["-"] == ([900], [950])  # strand "." is on both strands
    assert chr1.classify([(250, 300)], "+") == INTRONIC
    assert chr1.classify([(250, 300)], "-") == ANTISENSE
    assert chr1.classify([(910, 920)], "-") == EXONIC


def test_classify_rules():
    ca = load_gtf(TINY_GTF)["chr1"]
    assert ca.classify([(1100, 1200)], "+") == EXONIC
    assert ca.classify([(1100, 1200)], "-") == ANTISENSE
    assert ca.classify([(2000, 2100)], "+") == INTRONIC
    assert ca.classify([(4000, 4100)], "+") == INTERGENIC
    assert ca.classify([(1450, 1500), (2500, 2550)], "+") == EXONIC
    assert ca.classify([(1450, 1500)], "+") == EXONIC
    assert ca.classify([(1450, 1550)], "+") == EXONIC  # exactly 50 %
    assert ca.classify([(1449, 1550)], "+") == EXONIC
    assert ca.classify([(1451, 1550)], "+") == INTRONIC
    assert ContigAnnot().classify([(0, 10)], "+") == INTERGENIC
    # exon of G4 inside the intron of G1
    assert ca.classify([(1810, 1860)], "+") == EXONIC
    # unstranded: both strands are sense
    assert ca.classify([(1100, 1200)], None) == EXONIC
    assert ca.classify([(2000, 2100)], None) == INTRONIC
    assert ca.classify([(5100, 5200)], None) == EXONIC
    assert ca.classify([(4000, 4100)], None) == INTERGENIC


@pytest.fixture
def gtf_bam(tmp_path):
    return make_bam(tmp_path / "gtf.bam", gtf_reads(), contigs=GTF_CONTIGS)


def _expected_class(cb: str, flip: bool = False) -> str:
    cls = cb.split("-")[0]
    if flip:
        cls = {"exonic": "antisense", "antisense": "exonic"}.get(cls, cls)
    return cls


@pytest.mark.parametrize("workers", [1, 3])
def test_gtf_path_forward(gtf_bam, workers):
    stats = collect_metrics(
        str(gtf_bam), workers=workers, extra_metrics=True, gtf=str(TINY_GTF)
    )
    df = to_dataframe(stats).set_index("cb")
    assert (df["region_source"] == "gtf").all()
    for cb, row in df.iterrows():
        cls = _expected_class(cb)
        assert row[f"frac_{cls}"] == 1.0, (cb, row.filter(like="frac_"))
    # no UMI tag, no xf tag, no pa/ts tags
    assert df["n_umi"].isna().all()
    assert df["frac_singleton_umi"].isna().all()
    assert df["n_umi_counted"].isna().all()
    assert (df["polya_tso_source"] == "seq").all()
    assert (df["frac_polya"] == 0).all()


def test_gtf_path_reverse_strand(gtf_bam):
    stats = collect_metrics(
        str(gtf_bam),
        workers=2,
        extra_metrics=True,
        gtf=str(TINY_GTF),
        strand="reverse",
    )
    df = to_dataframe(stats).set_index("cb")
    # Reverse library: sense and antisense swap. Intronic reads become
    # antisense (the gene is on the other strand), so check the clear cases.
    for cb in ["exonic-1", "antisense-1", "exonic-2", "antisense-3", "exonic-3"]:
        cls = _expected_class(cb, flip=True)
        assert df.loc[cb, f"frac_{cls}"] == 1.0, cb
    assert df.loc["exonic-7", "frac_antisense"] == 1.0
    assert df.loc["intronic-1", "frac_antisense"] == 1.0
    assert df.loc["antisense-2", "frac_intronic"] == 1.0
    assert df.loc["intergenic-4", "frac_intergenic"] == 1.0
    # strand "." gene: sense on both strands
    assert df.loc["exonic-9", "frac_exonic"] == 1.0
    assert df.loc["exonic-10", "frac_exonic"] == 1.0


def test_gtf_path_unstranded(gtf_bam):
    stats = collect_metrics(
        str(gtf_bam),
        workers=2,
        extra_metrics=True,
        gtf=str(TINY_GTF),
        strand="unstranded",
    )
    df = to_dataframe(stats).set_index("cb")
    assert (df["frac_antisense"] == 0).all()
    unstranded = {
        "antisense-1": "exonic",
        "antisense-2": "intronic",
        "antisense-3": "exonic",
        "antisense-4": "intronic",
    }
    for cb, row in df.iterrows():
        cls = unstranded.get(cb, _expected_class(cb))
        assert row[f"frac_{cls}"] == 1.0, (cb, row.filter(like="frac_"))


def test_gtf_path_read2_of_pair(tmp_path):
    """Read 2 of a pair has the opposite strand of the transcript."""
    reads = [
        dict(contig="chr1", pos=1100, cigar="100M", flag=0x1 | 0x80, reverse=True,
             tags={"CB": "r2-rev"}),
        dict(contig="chr1", pos=1100, cigar="100M", flag=0x1 | 0x40, reverse=False,
             tags={"CB": "r1-fwd"}),
        dict(contig="chr1", pos=1100, cigar="100M", flag=0x1 | 0x80, reverse=False,
             tags={"CB": "r2-fwd"}),
        dict(contig="chr1", pos=1100, cigar="100M", flag=0x1 | 0x40, reverse=True,
             tags={"CB": "r1-rev"}),
    ]  # fmt: skip
    bam = make_bam(tmp_path / "pe.bam", reads)
    stats = collect_metrics(str(bam), workers=1, extra_metrics=True, gtf=TINY_GTF)
    df = to_dataframe(stats).set_index("cb")
    assert df.loc["r2-rev", "frac_exonic"] == 1.0
    assert df.loc["r1-fwd", "frac_exonic"] == 1.0
    assert df.loc["r2-fwd", "frac_antisense"] == 1.0
    assert df.loc["r1-rev", "frac_antisense"] == 1.0


def test_strand_must_be_valid(gtf_bam):
    with pytest.raises(ValueError, match="strand"):
        collect_metrics(str(gtf_bam), workers=1, extra_metrics=True, strand="both")


def test_no_region_source_gives_nan(gtf_bam):
    stats = collect_metrics(str(gtf_bam), workers=1, extra_metrics=True)
    df = to_dataframe(stats)
    assert (df["region_source"] == "none").all()
    for col in ["frac_exonic", "frac_intronic", "frac_intergenic", "frac_antisense"]:
        assert df[col].isna().all()
    for col in ["n_umi", "mean_reads_per_umi", "frac_singleton_umi"]:
        assert df[col].isna().all()
    assert (df["frac_nm0"] == 1.0).all()
    assert (df["frac_nm_ge2"] == 0.0).all()


def test_no_nm_tag_gives_nan(tmp_path):
    reads = [dict(contig="chr1", pos=10, cigar="50M", tags={"CB": "A"})]
    bam = make_bam(tmp_path / "nonm.bam", reads)
    df = to_dataframe(collect_metrics(str(bam), workers=1, extra_metrics=True))
    assert df["frac_nm0"].isna().all() and df["frac_nm_ge2"].isna().all()
    assert df["frac_mapq_unique"].item() == 1.0


def test_empty_table_has_extra_columns(gtf_bam):
    stats = collect_metrics(str(gtf_bam), workers=1, extra_metrics=True)
    only_nocb = {k: v for k, v in stats.items() if k[0] == "exonic-1"}
    only_nocb = {(NO_BARCODE, k[1]): v for k, v in only_nocb.items()}
    df = to_dataframe(only_nocb)
    assert df.empty
    assert list(df.columns) == [*DEFAULT_COLUMNS, *EXTRA_COLUMNS]


@pytest.mark.skipif(not TOXODB_GTF.exists(), reason=f"{TOXODB_GTF} not found")
def test_load_gtf_toxodb():
    annot = load_gtf(TOXODB_GTF)
    chroms = [c for c in annot if c.startswith("TGME49_chr")]
    assert len(chroms) == 14
    # 436 annotated contigs: 14 chromosomes and 422 KE* contigs
    assert len(annot) == 436
    assert sum(c.startswith("KE") for c in annot) == 422
    chr1a = annot["TGME49_chrIa"]
    # TGME49_500051 (+): exon 6510-6758; TGME49_287280 (-): exons 10836-10992,
    # 11304-11310 (1-based)
    assert chr1a.classify([(6600, 6650)], "+") == EXONIC
    assert chr1a.classify([(6600, 6650)], "-") == ANTISENSE
    assert chr1a.classify([(11100, 11150)], "-") == INTRONIC
    assert chr1a.classify([(11100, 11150)], "+") == ANTISENSE
    assert chr1a.classify([(8000, 8050)], "+") == INTERGENIC
    # KE138841: TGME49_302060 (-), 1-454, one exon
    assert annot["KE138841"].classify([(100, 150)], "-") == EXONIC


# ---------------------------------------------------------------------------
# PolyA and TSO
# ---------------------------------------------------------------------------


def test_scan_clips_unit():
    g = GENOMIC
    # plus=True: transcript on +; plus=False: on -; None: both
    assert scan_clips(g + "A" * 10, 0, 10, True) == (True, False)
    assert scan_clips(g + "A" * 9, 0, 9, True) == (False, False)
    assert scan_clips("T" * 10 + g, 10, 0, False) == (True, False)
    assert scan_clips("A" * 10 + g, 10, 0, False) == (False, False)
    assert scan_clips(TSO + g, 30, 0, True) == (False, True)
    assert scan_clips(TSO[:12] + g, 12, 0, True) == (False, True)
    assert scan_clips(TSO[-12:] + g, 12, 0, True) == (False, True)
    assert scan_clips(TSO[-11:] + g, 11, 0, True) == (False, False)
    assert scan_clips(g + revcomp(TSO), 0, 30, False) == (False, True)
    assert scan_clips(g + TSO, 0, 30, False) == (False, False)
    assert scan_clips(g, 0, 0, True) == (False, False)
    assert scan_clips("T" * 10 + g + "A" * 10, 10, 10, True) == (True, False)
    assert scan_clips("T" * 10 + g, 10, 0, None) == (True, False)
    assert scan_clips(g + "A" * 10, 0, 10, None) == (True, False)
    assert scan_clips(g + revcomp(TSO), 0, 30, None) == (False, True)
    assert scan_clips(TSO + g, 30, 0, None) == (False, True)
    assert scan_clips(g + TSO, 0, 30, None) == (False, False)


def test_polya_tso_sequence_path(tmp_path):
    bam = make_bam(tmp_path / "clip.bam", clip_reads())
    stats = collect_metrics(str(bam), workers=1, extra_metrics=True)
    df = to_dataframe(stats).set_index("cb")
    assert (df["polya_tso_source"] == "seq").all()
    assert {"polyA-hardclip", "tso-hardclip"} <= set(df.index)
    for cb, row in df.iterrows():
        assert row["frac_polya"] == (1.0 if cb.startswith("polyA") else 0.0), cb
        assert row["frac_tso"] == (1.0 if cb.startswith("tso") else 0.0), cb


@pytest.mark.parametrize(
    "strand, polya, tso",
    [
        # forward alignments; transcript strand + (forward) or - (reverse)
        ("forward", {"trailA", "r2-leadT"}, {"leadTSO"}),
        ("reverse", {"leadT"}, {"trailTSOrc"}),
        ("unstranded", {"leadT", "trailA", "r2-leadT"}, {"leadTSO", "trailTSOrc"}),
    ],
)
def test_polya_tso_follow_library_strand(tmp_path, strand, polya, tso):
    bam = make_bam(tmp_path / "s.bam", strand_clip_reads())
    stats = collect_metrics(str(bam), workers=1, extra_metrics=True, strand=strand)
    df = to_dataframe(stats).set_index("cb")
    assert set(df.index[df["frac_polya"] == 1.0]) == polya
    assert set(df.index[df["frac_tso"] == 1.0]) == tso
    assert set(df.index[df["frac_polya"] == 0.0]) == set(df.index) - polya


def test_polya_tso_tags_path_ignores_sequence(tmp_path):
    """With pa/ts tags in the BAM, the clipped sequence is not scanned."""
    reads = clip_reads() + [
        dict(contig="chr2", pos=10, cigar="50M", seq=GENOMIC,
             tags={"CB": "tagged", "pa": 15, "ts": 30}),
        dict(contig="chr2", pos=20, cigar="50M", seq=GENOMIC,
             tags={"CB": "tagged", "pa": 0}),
    ]  # fmt: skip
    bam = make_bam(tmp_path / "tagged.bam", reads)
    stats = collect_metrics(str(bam), workers=2, extra_metrics=True)
    df = to_dataframe(stats).set_index("cb")
    assert (df["polya_tso_source"] == "tags").all()
    assert df.loc["tagged", "frac_polya"] == 0.5
    assert df.loc["tagged", "frac_tso"] == 0.5
    assert (df.drop(index="tagged")[["frac_polya", "frac_tso"]] == 0).all().all()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def test_cli_extra_metrics_gtf(gtf_bam, tmp_path, monkeypatch, capsys):
    text = _run_cli(
        gtf_bam,
        tmp_path,
        monkeypatch,
        "--extra-metrics",
        "--gtf",
        str(TINY_GTF),
        "--strand",
        "forward",
        "--umi-tag",
        "UB",
        "--mapq-unique",
        "255",
    )
    out = " ".join(capsys.readouterr().out.split())  # undo the line wrap
    assert "regions: GTF" in out
    assert "2/3 contigs in the GTF" in out  # chr3 is not in the GTF
    assert "No read has the UMI tag" in out
    header = text.splitlines()[0].split(",")
    assert header == [*DEFAULT_COLUMNS, *EXTRA_COLUMNS]


def test_cli_warns_on_extra_options_without_flag(cr_bam, tmp_path, monkeypatch, capsys):
    text = _run_cli(cr_bam, tmp_path, monkeypatch, "--umi-tag", "UR")
    assert "have no effect" in " ".join(capsys.readouterr().out.split())
    assert text.splitlines()[0].split(",") == DEFAULT_COLUMNS


def test_cli_missing_gtf(cr_bam, tmp_path, monkeypatch):
    from scbamqc.cli import app

    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit) as e:
        app(["bam", str(cr_bam), "--extra-metrics", "--gtf", "missing.gtf"])
    assert e.value.code == 1


def test_umi_counter_packed_equals_strings():
    rng = np.random.default_rng(7)
    bases = np.array(list("ACGT"))
    umis = ["".join(rng.choice(bases, 12)) for _ in range(300)]
    umis += ["ACGTNACGTACG", "ACGTACGTACGTAC", "ACG", "", "acgtacgtacgt"]
    barcodes = [f"BC{i}" for i in range(40)] + [NO_BARCODE]
    packed, strings = UmiCounter(), UmiCounter(packed=False)
    for _ in range(5000):
        cb = barcodes[rng.integers(len(barcodes))]
        umi = umis[rng.integers(len(umis))]
        packed.add(cb, umi)
        strings.add(cb, umi)
    assert len(packed._keys) > 0 and len(strings._keys) == 0
    assert packed.per_barcode() == strings.per_barcode()
    # same UMI bases under two barcodes are two keys
    c = UmiCounter()
    for cb, umi in [("A", "AAAA"), ("B", "AAAA"), ("A", "AAAA"), ("A", "AAAN")]:
        c.add(cb, umi)
    assert c.per_barcode() == {"A": (2, 3, 1), "B": (1, 1, 1)}


def test_source_detection_uses_selected_contigs(tmp_path):
    reads = gtf_reads()
    for r in reads:
        if r["contig"] == "chr2":
            r["tags"]["RE"] = "E"
    bam = make_bam(tmp_path / "re2.bam", reads, contigs=GTF_CONTIGS)
    only1 = collect_metrics(
        str(bam), workers=1, extra_metrics=True, contig_list=["chr1"]
    )
    only2 = collect_metrics(
        str(bam), workers=1, extra_metrics=True, contig_list=["chr2"]
    )
    assert (to_dataframe(only1)["region_source"] == "none").all()
    assert (to_dataframe(only2)["region_source"] == "tags").all()


def test_umi_tag_integer(tmp_path):
    """A UMI tag with integer values counts like strings."""
    reads = [
        dict(contig="chr1", pos=10 + i, cigar="50M", tags={"CB": "A", "UB": u})
        for i, u in enumerate([7, 7, 8, 9])
    ]
    bam = make_bam(tmp_path / "intumi.bam", reads)
    df = to_dataframe(collect_metrics(str(bam), workers=1, extra_metrics=True))
    row = df.set_index("cb").loc["A"]
    assert row["n_umi"] == 3
    assert row["mean_reads_per_umi"] == round(4 / 3, 3)
    assert row["frac_singleton_umi"] == round(2 / 3, 4)
    c = UmiCounter()
    for u in [7, "7", "ACGT"]:
        c.add("A", u)
    assert c.per_barcode() == {"A": (2, 3, 1)}


def test_region_source_forced_gtf_on_tagged_bam(cr_bam, gtf_bam):
    # cr_bam has RE tags; tiny.gtf has G1 (+) exon [1000, 1500) on chr1 only,
    # so the cr_bam reads (chr1 100-700, chr2) are intergenic with the GTF.
    auto = to_dataframe(
        collect_metrics(str(cr_bam), workers=1, extra_metrics=True, gtf=TINY_GTF)
    )
    assert (auto["region_source"] == "tags").all()
    forced = to_dataframe(
        collect_metrics(
            str(cr_bam),
            workers=1,
            extra_metrics=True,
            gtf=TINY_GTF,
            region_source="gtf",
        )
    ).set_index("cb")
    assert (forced["region_source"] == "gtf").all()
    assert (forced["frac_intergenic"] == 1.0).all()
    # the same reads as gtf_bam give the same result as auto with a GTF
    kw = dict(workers=1, extra_metrics=True, gtf=TINY_GTF)
    a = to_dataframe(collect_metrics(str(gtf_bam), **kw))
    g = to_dataframe(collect_metrics(str(gtf_bam), region_source="gtf", **kw))
    pd.testing.assert_frame_equal(a, g)


def test_region_source_tags_on_gtf_bam(cr_bam, gtf_bam):
    with pytest.raises(ValueError, match="no RE tag"):
        collect_metrics(
            str(gtf_bam), workers=1, extra_metrics=True, region_source="tags"
        )
    df = to_dataframe(
        collect_metrics(
            str(cr_bam), workers=1, extra_metrics=True, region_source="tags"
        )
    )
    assert (df["region_source"] == "tags").all()


def test_region_source_gtf_needs_gtf(cr_bam, tmp_path, monkeypatch, capsys):
    with pytest.raises(ValueError, match="needs a GTF"):
        collect_metrics(str(cr_bam), workers=1, extra_metrics=True, region_source="gtf")
    with pytest.raises(ValueError, match="region_source must be"):
        collect_metrics(str(cr_bam), workers=1, extra_metrics=True, region_source="x")

    from scbamqc.cli import app

    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit) as e:
        app(["bam", str(cr_bam), "--extra-metrics", "--region-source", "gtf"])
    assert e.value.code == 1
    assert "needs --gtf" in capsys.readouterr().err


def test_cli_region_source_gtf(cr_bam, tmp_path, monkeypatch, capsys):
    text = _run_cli(
        cr_bam,
        tmp_path,
        monkeypatch,
        "--extra-metrics",
        "--gtf",
        str(TINY_GTF),
        "--region-source",
        "gtf",
    )
    out = " ".join(capsys.readouterr().out.split())  # undo the line wrap
    assert "the RE tags are not used" in out
    df = pd.read_csv(StringIO(text))
    assert (df["region_source"] == "gtf").all()
