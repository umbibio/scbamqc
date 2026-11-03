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

    Each read dict: contig, pos, cigar, mapq, qual (int), tags (dict).
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
            a.query_sequence = "A" * qlen
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
