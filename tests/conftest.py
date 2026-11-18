import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("SCBAMQC_TEST_DATA", REPO.parent.parent / "data"))


def make_bam(
    path: Path, reads: list[dict], contigs=(("chr1", 10_000), ("chr2", 10_000))
):
    """Write and index a small coordinate-sorted BAM.

    Each read dict: contig, pos, cigar, mapq, qual (int), tags (dict), and
    optional seq (str, default all ``A``), reverse (bool), flag (int, OR-ed
    into the SAM flag, e.g. 256 for a secondary alignment).
    """
    pysam = pytest.importorskip("pysam")
    from scbamqc.bam import read_length_from_cigar

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": n, "LN": ln} for n, ln in contigs],
    }
    names = [n for n, _ in contigs]
    unsorted = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for i, r in enumerate(reads):
            a = pysam.AlignedSegment(out.header)
            a.query_name = f"r{i}"
            a.reference_id = names.index(r["contig"])
            a.reference_start = r["pos"]
            a.cigarstring = r["cigar"]
            qlen = read_length_from_cigar(r["cigar"])
            seq = r.get("seq", "A" * qlen)
            assert len(seq) == qlen, (r, len(seq), qlen)
            a.flag |= r.get("flag", 0)
            a.is_reverse = r.get("reverse", False)
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array(
                chr(r.get("qual", 30) + 33) * qlen
            )
            a.mapping_quality = r.get("mapq", 255)
            a.set_tags(list(r.get("tags", {}).items()))
            out.write(a)
    pysam.sort("-o", str(path), str(unsorted))
    pysam.index(str(path))
    unsorted.unlink()
    return path


def synthetic_obs(n=3000, seed=0) -> pd.DataFrame:
    """BAM-QC-like table: a tight 'good' population plus a diffuse background."""
    rng = np.random.default_rng(seed)
    n_good = n // 2
    good = pd.DataFrame(
        {
            "n_reads": rng.integers(500, 5000, n_good),
            "mean_nm": rng.uniform(0, 0.8, n_good),
            "mean_clip_frac": rng.normal(0.03, 0.01, n_good).clip(0),
            "std_clip_frac": rng.normal(0.08, 0.02, n_good).clip(0),
            "mean_base_qual": rng.normal(38, 0.5, n_good),
            "total_counts": rng.lognormal(7, 0.5, n_good),
        }
    )
    m = n - n_good
    bad = pd.DataFrame(
        {
            "n_reads": rng.integers(1, 3000, m),
            "mean_nm": rng.uniform(0, 2, m),
            "mean_clip_frac": rng.uniform(0, 0.5, m),
            "std_clip_frac": rng.uniform(0, 0.3, m),
            "mean_base_qual": rng.uniform(25, 39, m),
            "total_counts": rng.lognormal(4, 1, m),
        }
    )
    df = pd.concat([good, bad], ignore_index=True)
    df.index = [f"BC{i:05d}-1" for i in range(n)]
    return df


# ---------------------------------------------------------------------------
# Fixtures for the extra metrics
# ---------------------------------------------------------------------------

TESTS_DATA = Path(__file__).resolve().parent / "data"
TINY_GTF = TESTS_DATA / "tiny.gtf"
GOLDEN_DIR = TESTS_DATA / "golden"
TOXODB_GTF = Path(
    os.environ.get("SCBAMQC_TOXODB_GTF", DATA_DIR / "ToxoDB-68_TgondiiME49.gtf")
)

TSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG"
# 50 nt that contain no run of A/T and no TSO 12-mer
GENOMIC = "ACGTCCGTAGCATGCGACTGCAGTCGATCGGCTAGCTGACGTGCACGTCA"


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def cellranger_reads() -> list[dict]:
    """Cell Ranger-like reads: RE/AN/GX/TX, UB, pa/ts tags.

    Cell AAAA-1 (7 alignments, one of them secondary):
      chr1  E  UB=U1  nM=0                 sense (GX, TX)
      chr1  E  UB=U1  nM=0
      chr1  N  UB=U2  nM=1  pa=12
      chr1  I  UB=U3  nM=2  ts=5  MAPQ 3
      chr2  E  UB=U1  nM=3  AN, no GX/TX  -> antisense
      chr2  E  UB=U4  nM=0  AN and GX     -> exonic (sense wins)
      chr2  -  UB=U4  nM=0  secondary, MAPQ 0, no RE
    xf: 25 (bit 8 set) on 3 alignments of AAAA-1 and 1 of CCCC-1.
    Cell CCCC-1: 2 reads on chr1, no UB.
    One read without CB.
    """
    gx = {"GX": "G1", "TX": "T1,+100,100M"}
    a = "AAAA-1"
    return [
        dict(contig="chr1", pos=100, cigar="100M",
             tags={"CB": a, "UB": "U1", "nM": 0, "NH": 1, "RE": "E", "xf": 25,
                   **gx}),
        dict(contig="chr1", pos=110, cigar="100M",
             tags={"CB": a, "UB": "U1", "nM": 0, "NH": 1, "RE": "E", "xf": 17,
                   **gx}),
        dict(contig="chr1", pos=120, cigar="90M10S",
             tags={"CB": a, "UB": "U2", "nM": 1, "NH": 1, "RE": "N", "pa": 12,
                   "xf": 25}),
        dict(contig="chr1", pos=130, cigar="10S90M", mapq=3,
             tags={"CB": a, "UB": "U3", "nM": 2, "NH": 2, "RE": "I", "ts": 5,
                   "xf": 0}),
        dict(contig="chr2", pos=100, cigar="100M",
             tags={"CB": a, "UB": "U1", "nM": 3, "NH": 1, "RE": "E",
                   "AN": "T9,-100,100M", "xf": 0}),
        dict(contig="chr2", pos=200, cigar="100M",
             tags={"CB": a, "UB": "U4", "nM": 0, "NH": 1, "RE": "E",
                   "AN": "T9,-200,100M", "xf": 25, **gx}),
        dict(contig="chr2", pos=300, cigar="100M", mapq=0, flag=256,
             tags={"CB": a, "UB": "U4", "nM": 0, "NH": 2}),
        dict(contig="chr1", pos=500, cigar="100M", qual=35,
             tags={"CB": "CCCC-1", "nM": 0, "NH": 1, "RE": "E", "xf": 25, **gx}),
        dict(contig="chr1", pos=600, cigar="100M", qual=35,
             tags={"CB": "CCCC-1", "nM": 1, "NH": 1, "RE": "I"}),
        dict(contig="chr1", pos=900, cigar="100M",
             tags={"UB": "U9", "nM": 0, "RE": "E", **gx}),
    ]  # fmt: skip


# tiny.gtf (1-based, closed):
#   chr1  G1 +  1001-3000  exons 1001-1500, 2501-3000
#   chr1  G4 +  1801-1900  exon  1801-1900 (inside the G1 intron)
#   chr1  G2 -  5001-6000  exon  5001-6000
#   chr2  G3 -  1001-2000  exons 1001-1200, 1801-2000
#   chr2  G5 .  5001-5500  exon  5001-5500 (no strand: both strands)
# chr3 is not in the GTF.
GTF_CASES = [
    # (barcode = expected class + id, contig, pos (0-based), cigar, reverse)
    ("exonic-1", "chr1", 1100, "100M", False),
    ("antisense-1", "chr1", 1100, "100M", True),
    ("intronic-1", "chr1", 2000, "100M", False),
    ("antisense-2", "chr1", 2000, "100M", True),
    ("intergenic-1", "chr1", 4000, "100M", False),
    ("intergenic-2", "chr1", 4000, "100M", True),
    ("exonic-2", "chr1", 5100, "100M", True),
    ("antisense-3", "chr1", 5100, "100M", False),
    ("exonic-3", "chr2", 1100, "100M", True),
    ("intronic-2", "chr2", 1500, "100M", True),
    ("antisense-4", "chr2", 1500, "100M", False),
    # spliced: blocks [1450, 1500) and [2500, 2550), both in G1 exons
    ("exonic-4", "chr1", 1450, "50M1000N50M", False),
    # 40 of 100 bases in an exon, 60 in the intron
    ("intronic-3", "chr1", 1460, "100M", False),
    # 60 of 100 bases in an exon, 40 in the intron
    ("exonic-5", "chr1", 1440, "100M", False),
    # soft clips do not count: 80 aligned bases, all in the exon
    ("exonic-6", "chr1", 1420, "80M20S", False),
    # 40 bases in the exon, 60 after the gene end: not exonic, not intronic
    ("intergenic-3", "chr1", 2960, "100M", False),
    ("intergenic-4", "chr3", 100, "100M", False),
    # exon of G4 inside the intron of G1 (same strand)
    ("exonic-7", "chr1", 1810, "50M", False),
    # deletion and splice gap: blocks [1450, 1470), [1475, 1500), [2500, 2550)
    ("exonic-8", "chr1", 1450, "20M5D25M1000N50M", False),
    # blocks [1550, 1600) intron and [2600, 2640) exon: 40 of 90 bases exonic.
    # The N gap covers the G4 exon and the G1 exon start; they do not count.
    ("intronic-4", "chr1", 1550, "50M1000N40M", False),
    # gene with strand "." is sense on both strands
    ("exonic-9", "chr2", 5100, "100M", False),
    ("exonic-10", "chr2", 5100, "100M", True),
]


def gtf_reads() -> list[dict]:
    """STAR-like reads (no RE, no UB, no pa/ts) for the GTF region tests."""
    from scbamqc.bam import read_length_from_cigar

    return [
        dict(
            contig=c,
            pos=p,
            cigar=cg,
            reverse=rev,
            seq=("ACGT" * 50)[: read_length_from_cigar(cg)],
            tags={"CB": cb, "nM": 0},
        )
        for cb, c, p, cg, rev in GTF_CASES
    ]


GTF_CONTIGS = (("chr1", 10_000), ("chr2", 10_000), ("chr3", 10_000))


def clip_reads() -> list[dict]:
    """Reads for the polyA/TSO soft-clip scan (no pa/ts tags).

    Barcode = expected result: polyA-*, tso-* or none-*.
    """
    g = GENOMIC
    return [
        # read 3' clip with A10+ (forward: trailing clip)
        dict(contig="chr1", pos=100, cigar="50M20S", seq=g + "A" * 15 + "CCCCC",
             tags={"CB": "polyA-fwd"}),
        # reverse alignment: read 3' end is the leading clip, shown as T
        dict(contig="chr1", pos=200, cigar="20S50M", seq="GGGGG" + "T" * 15 + g,
             reverse=True, tags={"CB": "polyA-rev"}),
        # A run in the read 5' clip is not polyA
        dict(contig="chr1", pos=300, cigar="15S50M", seq="A" * 15 + g,
             tags={"CB": "none-polyA-5p"}),
        # A run of 9 is too short
        dict(contig="chr1", pos=400, cigar="50M9S", seq=g + "A" * 9,
             tags={"CB": "none-polyA-short"}),
        # full TSO in the read 5' clip
        dict(contig="chr1", pos=500, cigar="30S50M", seq=TSO + g,
             tags={"CB": "tso-fwd"}),
        # reverse alignment: read 5' end is the trailing clip, reverse complement
        dict(contig="chr1", pos=600, cigar="50M30S", seq=g + revcomp(TSO),
             reverse=True, tags={"CB": "tso-rev"}),
        # 14-nt TSO suffix only
        dict(contig="chr1", pos=700, cigar="14S50M", seq=TSO[-14:] + g,
             tags={"CB": "tso-suffix"}),
        # TSO with one mismatch in the middle
        dict(contig="chr1", pos=800, cigar="30S50M",
             seq=TSO[:15] + ("C" if TSO[15] != "C" else "G") + TSO[16:] + g,
             tags={"CB": "tso-mismatch"}),
        # TSO in the read 3' clip is not counted
        dict(contig="chr1", pos=900, cigar="50M30S", seq=g + TSO,
             tags={"CB": "none-tso-3p"}),
        # 11-nt TSO suffix is too short
        dict(contig="chr1", pos=1000, cigar="11S50M", seq=TSO[-11:] + g,
             tags={"CB": "none-tso-short"}),
        dict(contig="chr1", pos=1100, cigar="50M", seq=g, tags={"CB": "none-noclip"}),
        # hard clip outside the soft clip
        dict(contig="chr1", pos=1200, cigar="50M20S5H", seq=g + "A" * 15 + "CCCCC",
             tags={"CB": "polyA-hardclip"}),
        dict(contig="chr1", pos=1300, cigar="5H30S50M", seq=TSO + g,
             tags={"CB": "tso-hardclip"}),
        dict(contig="chr2", pos=100, cigar="50M", seq=g, tags={"CB": "none-chr2"}),
    ]  # fmt: skip


def strand_clip_reads() -> list[dict]:
    """Forward alignments for the polyA/TSO strand tests (no pa/ts tags)."""
    g = GENOMIC
    return [
        dict(contig="chr1", pos=100, cigar="20S50M", seq="CCCCC" + "T" * 15 + g,
             tags={"CB": "leadT"}),
        dict(contig="chr1", pos=200, cigar="50M20S", seq=g + "A" * 15 + "CCCCC",
             tags={"CB": "trailA"}),
        dict(contig="chr1", pos=300, cigar="50M30S", seq=g + revcomp(TSO),
             tags={"CB": "trailTSOrc"}),
        dict(contig="chr1", pos=400, cigar="30S50M", seq=TSO + g,
             tags={"CB": "leadTSO"}),
        # read 2 of a pair: opposite strand of read 1
        dict(contig="chr1", pos=500, cigar="20S50M", seq="CCCCC" + "T" * 15 + g,
             flag=0x1 | 0x80, tags={"CB": "r2-leadT"}),
    ]  # fmt: skip


def golden_fixtures() -> dict[str, tuple[list[dict], tuple]]:
    """BAM fixtures with a golden default output in ``tests/data/golden``.

    The golden files were written by scbamqc 0.4.0.dev0 before the extra
    metrics were added (``tests/data/golden/make_golden.py``).
    """
    from test_bam import _reads

    default_contigs = (("chr1", 10_000), ("chr2", 10_000))
    return {
        "basic": (_reads(), default_contigs),
        "cellranger": (cellranger_reads(), default_contigs),
        "gtf": (gtf_reads(), GTF_CONTIGS),
        "clip": (clip_reads(), default_contigs),
    }


GOLDEN_RUNS = {
    # name: extra CLI arguments
    "cell": [],
    "cell-contig": ["--group-by", "cell-contig"],
    "cell-nocb": ["--include-no-cb"],
    "nm-NM": ["--nm-tag", "NM"],
    "barcode-UB": ["--barcode-tag", "UB", "--include-no-cb"],
    "contigs-chr2": ["--contigs", "chr2"],
    "prefix-chr1": ["--contig-prefix", "chr1", "--group-by", "cell-contig"],
}
