"""Streaming per-barcode BAM QC.

Works with Cell Ranger BAMs (``CB`` tag) and STARsolo SmartSeq BAMs (``RG`` tag).
Workers stream reads one contig at a time and accumulate per-(cell, contig)
statistics with Welford's online algorithm, so no per-read rows are stored.
Peak memory is O(unique barcodes x contigs touched), not O(reads).

``pysam`` is an optional dependency (``pip install scbamqc[bam]``) and is only
imported when a BAM file is read.

``collect_metrics(..., extra_metrics=True)`` adds count-based metrics
(:class:`ExtraStats`): read regions, UMI counts, polyA and TSO artefacts and
per-read MAPQ and mismatch distributions. Without it, no extra per-read work is
done and the output is the same as ``bam_qc_summary.py`` v0.2.1.
"""

from __future__ import annotations

import multiprocessing as mp
import re
from array import array
import threading
from collections import defaultdict
from dataclasses import dataclass, field
from os import PathLike
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd

from .annot import ANTISENSE, EXONIC, INTERGENIC, INTRONIC, ContigAnnot, load_gtf

if TYPE_CHECKING:
    from rich.console import Console
    from rich.progress import Progress

__all__ = [
    "EXTRA_COLUMNS",
    "NO_BARCODE",
    "TSO_10X",
    "CellStats",
    "ExtraConfig",
    "ExtraStats",
    "WelfordAcc",
    "collapse_to_cell",
    "collect_metrics",
    "count_tags",
    "filter_contigs",
    "parse_soft_clips",
    "UmiCounter",
    "read_length_from_cigar",
    "scan_clips",
    "to_dataframe",
]

GroupBy = Literal["cell", "cell-contig"]
StatsKey = tuple[str, str]  # (cell_id, contig)

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
NO_BARCODE = "__no_barcode__"
_REPORT_EVERY = 50_000  # a worker reports progress every N processed reads
_TAG_SAMPLE = 20_000  # mapped reads sampled to detect the tags of a BAM

#: 10x Genomics template switch oligo (Cell Ranger ``ts`` tag documentation).
TSO_10X = "AAGCAGTGGTATCAACGCAGAGTACATGGG"
POLYA_MIN_RUN = 10  # A bases in a row in the 3' soft clip
TSO_MIN_MATCH = 12  # TSO bases in a row in the 5' soft clip


def _import_pysam():
    try:
        import pysam
    except ImportError as e:  # pragma: no cover - depends on the environment
        raise ImportError(
            "Reading BAM files requires pysam. Install it with "
            "`pip install scbamqc[bam]`."
        ) from e
    return pysam


# ---------------------------------------------------------------------------
# CIGAR helpers
# ---------------------------------------------------------------------------


def parse_soft_clips(cigar_str: str) -> tuple[int, int]:
    """Return ``(clip_5prime, clip_3prime)`` soft-clipped bases of a CIGAR string."""
    ops = CIGAR_RE.findall(cigar_str)
    if not ops:
        return 0, 0
    clip5 = int(ops[0][0]) if ops[0][1] == "S" else 0
    clip3 = int(ops[-1][0]) if ops[-1][1] == "S" else 0
    return clip5, clip3


def read_length_from_cigar(cigar_str: str) -> int:
    """Total query-consuming length (M + I + S + = + X)."""
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar_str) if op in "MISX=")


# ---------------------------------------------------------------------------
# PolyA / TSO soft-clip scan
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _kmer_regex(seq: str, k: int) -> re.Pattern:
    kmers = sorted({seq[i : i + k] for i in range(len(seq) - k + 1)})
    return re.compile("|".join(kmers))


_TSO_FWD = _kmer_regex(TSO_10X, TSO_MIN_MATCH)
_TSO_REV = _kmer_regex(_revcomp(TSO_10X), TSO_MIN_MATCH)
_POLYA = "A" * POLYA_MIN_RUN
_POLYT = "T" * POLYA_MIN_RUN


def scan_clips(
    seq: str, clip_left: int, clip_right: int, plus: bool | None
) -> tuple[bool, bool]:
    """Find polyA and TSO sequence in the soft clips of one alignment.

    ``seq`` is the BAM ``SEQ`` (reference orientation); ``clip_left`` and
    ``clip_right`` are the leading and trailing soft clips. ``plus`` is the
    genome strand of the transcript that the read comes from (``True``: ``+``),
    as for the GTF regions (library strand and read-2 flip applied). ``None``
    (unstranded library) checks both strands.

    The transcript 3' end has the polyA tail and the transcript 5' end has the
    TSO. In reference orientation:

    - transcript on ``+``: polyA is a run of ``POLYA_MIN_RUN`` (10) ``A`` in
      the trailing clip; TSO is ``TSO_MIN_MATCH`` (12) or more consecutive
      bases of ``TSO_10X`` in the leading clip.
    - transcript on ``-``: polyA is a run of 10 ``T`` in the leading clip; TSO
      is 12 or more consecutive bases of the reverse complement of ``TSO_10X``
      in the trailing clip.

    Any 12-mer of the TSO matches, so a TSO prefix, a TSO suffix or a full TSO
    with one mismatch is found.

    Returns
    -------
    ``(has_polya, has_tso)``
    """
    n = len(seq)
    left = seq[:clip_left] if clip_left else ""
    right = seq[n - clip_right :] if clip_right else ""
    polya = tso = False
    if plus is None or plus:
        polya = len(right) >= POLYA_MIN_RUN and _POLYA in right
        tso = len(left) >= TSO_MIN_MATCH and _TSO_FWD.search(left) is not None
    if plus is None or not plus:
        polya = polya or (len(left) >= POLYA_MIN_RUN and _POLYT in left)
        tso = tso or (
            len(right) >= TSO_MIN_MATCH and _TSO_REV.search(right) is not None
        )
    return polya, tso


def _scan_soft_clips(read, cigar: str, clip5: int, clip3: int) -> tuple[int, int]:
    """Leading and trailing soft clips, also when a hard clip is outside them.

    ``parse_soft_clips`` (default metrics) only looks at the first and last
    CIGAR operation, so ``5H10S80M`` gives no clip there. The scan uses the
    soft clip next to the hard clip.
    """
    if "H" not in cigar:
        return clip5, clip3
    ops = read.cigartuples
    i, j = 0, len(ops) - 1
    if ops[i][0] == 5:  # BAM_CHARD_CLIP
        i += 1
    if ops[j][0] == 5:
        j -= 1
    left = ops[i][1] if i <= j and ops[i][0] == 4 else 0  # BAM_CSOFT_CLIP
    right = ops[j][1] if j > i and ops[j][0] == 4 else 0
    return left, right


# ---------------------------------------------------------------------------
# UMI counting
# ---------------------------------------------------------------------------

_UMI_DIGITS = str.maketrans("ACGT", "0123")


class UmiCounter:
    """Reads per (barcode, UMI) on one contig, in a compact form.

    Barcodes get small integer ids. A UMI of ``A``/``C``/``G``/``T`` only is
    packed into 2 bits per base. Each read appends the 64-bit key
    ``barcode_id << (2 * L) | umi_bits`` to an ``array('Q')`` (8 bytes per
    read). ``L`` is the length of the first packed UMI (at most 28). A UMI with
    another letter (e.g. ``N``) or another length goes into a small dict with
    the string key ``"<barcode_id>\t<umi>"``. ``packed=False`` uses the dict
    for all UMIs (for tests). :meth:`per_barcode` sorts the keys in place,
    counts the runs and frees the keys, so call it once. Measured peak: about
    25 bytes per read at 1 read per UMI, 12 at 5 reads per UMI.
    """

    __slots__ = (
        "_bc_ids",
        "_bcs",
        "_keys",
        "_other",
        "_umi_len",
        "_max_bid",
        "_packed",
    )

    def __init__(self, packed: bool = True) -> None:
        self._bc_ids: dict = {}
        self._bcs: list = []
        self._keys = array("Q")
        self._other: dict[str, int] = {}
        self._umi_len = 0
        self._max_bid = 0
        self._packed = packed

    def add(self, cb, umi) -> None:
        if not isinstance(umi, str):
            umi = str(umi)  # e.g. an integer tag
        bid = self._bc_ids.get(cb)
        if bid is None:
            bid = self._bc_ids[cb] = len(self._bcs)
            self._bcs.append(cb)
        n = len(umi)
        if (
            self._packed
            and (n == self._umi_len or (self._umi_len == 0 and 0 < n <= 28))
            and umi.isalpha()
        ):
            digits = umi.translate(_UMI_DIGITS)
            if digits.isdigit():  # only A, C, G, T
                if self._umi_len == 0:
                    self._umi_len = n
                    self._max_bid = (1 << (64 - 2 * n)) - 1
                if bid <= self._max_bid:
                    self._keys.append((bid << (2 * n)) | int(digits, 4))
                    return
        key = f"{bid}\t{umi}"
        other = self._other
        other[key] = other.get(key, 0) + 1

    def per_barcode(self) -> dict:
        """``{barcode: (n_umi, n_reads, n_singleton_umi)}``."""
        nb = len(self._bcs)
        n_umi = np.zeros(nb, dtype=np.int64)
        n_reads = np.zeros(nb, dtype=np.int64)
        n_single = np.zeros(nb, dtype=np.int64)
        r = len(self._keys)
        if r:
            keys = np.frombuffer(self._keys, dtype=np.uint64)  # a view, no copy
            keys.sort()  # in place
            first = np.empty(r, dtype=bool)  # first key of each run
            first[0] = True
            np.not_equal(keys[1:], keys[:-1], out=first[1:])
            starts = np.flatnonzero(first)
            del first
            bids = keys[starts]
            bids >>= np.uint64(2 * self._umi_len)
            bids = bids.view(np.int64)
            del keys
            self._keys = array("Q")  # free the read keys before the counts
            counts = np.empty_like(starts)
            np.subtract(starts[1:], starts[:-1], out=counts[:-1])
            counts[-1] = r - starts[-1]
            del starts
            n_umi += np.bincount(bids, minlength=nb)
            n_reads += np.bincount(bids, weights=counts, minlength=nb).astype(np.int64)
            n_single += np.bincount(bids[counts == 1], minlength=nb)
        for key, v in self._other.items():
            bid = int(key.split("\t", 1)[0])
            n_umi[bid] += 1
            n_reads[bid] += v
            n_single[bid] += v == 1
        return {
            self._bcs[b]: (int(n_umi[b]), int(n_reads[b]), int(n_single[b]))
            for b in np.flatnonzero(n_umi)
        }


# ---------------------------------------------------------------------------
# Contig filtering
# ---------------------------------------------------------------------------


def filter_contigs(
    contigs: list[str],
    contig_list: list[str] | None = None,
    contig_prefixes: list[str] | None = None,
) -> list[str]:
    """Keep contigs in ``contig_list`` or starting with any of ``contig_prefixes``.

    The two filters are combined as a union. With neither given, ``contigs`` is
    returned unchanged. The input order is preserved.
    """
    if not contig_list and not contig_prefixes:
        return contigs

    keep: set[str] = set()
    if contig_list:
        allowed = set(contig_list)
        keep.update(c for c in contigs if c in allowed)
    if contig_prefixes:
        prefix_tup = tuple(contig_prefixes)
        keep.update(c for c in contigs if c.startswith(prefix_tup))

    return [c for c in contigs if c in keep]


# ---------------------------------------------------------------------------
# Accumulators
# ---------------------------------------------------------------------------


@dataclass
class WelfordAcc:
    """Online mean and population variance (Welford; merge by Chan et al.)."""

    n: int = 0
    mean: float = 0.0
    M2: float = 0.0

    def update(self, x: float) -> None:
        if not np.isfinite(x):
            return
        self.n += 1
        delta = x - self.mean
        self.mean += delta / self.n
        self.M2 += delta * (x - self.mean)

    @property
    def variance(self) -> float:
        return self.M2 / self.n if self.n > 1 else 0.0

    @property
    def std(self) -> float:
        return self.variance**0.5

    def merge(self, other: WelfordAcc) -> WelfordAcc:
        """Return an accumulator combining ``self`` and ``other``."""
        if self.n == 0:
            return other
        if other.n == 0:
            return self
        c = WelfordAcc()
        c.n = self.n + other.n
        delta = other.mean - self.mean
        c.mean = (self.n * self.mean + other.n * other.mean) / c.n
        c.M2 = self.M2 + other.M2 + delta**2 * self.n * other.n / c.n
        return c


RegionSource = Literal["tags", "gtf", "none"]
Strand = Literal["forward", "reverse", "unstranded"]
ArtefactSource = Literal["tags", "seq"]

#: Columns that ``extra_metrics=True`` appends, in order.
EXTRA_COLUMNS = (
    "region_source",
    "frac_exonic",
    "frac_intronic",
    "frac_intergenic",
    "frac_antisense",
    "n_umi",
    "mean_reads_per_umi",
    "frac_singleton_umi",
    "n_umi_counted",
    "polya_tso_source",
    "frac_polya",
    "frac_tso",
    "frac_mapq_unique",
    "frac_nm0",
    "frac_nm_ge2",
)
UMI_COLUMNS = ("n_umi", "mean_reads_per_umi", "frac_singleton_umi")


def _frac(num: int, den: int, digits: int = 4) -> float:
    return round(num / den, digits) if den else float("nan")


@dataclass
class ExtraConfig:
    """Settings of the extra metrics, shared by all workers.

    ``region_source`` and ``artefact_source`` are decided once for the BAM by
    :func:`collect_metrics`.
    """

    region_source: RegionSource = "none"
    artefact_source: ArtefactSource = "seq"
    strand: Strand = "forward"
    umi_tag: str = "UB"
    mapq_unique: int = 255


@dataclass
class ExtraStats:
    """Counters of the extra metrics for one (cell, contig) pair.

    All fields are counts, so two accumulators merge by addition.

    ``region`` counts reads by class: exonic, intronic, intergenic, antisense
    (the ``scbamqc.annot`` constants). The classes do not overlap; reads with
    no class (``RE`` tag missing) are not counted.

    UMI counts are reduced from (barcode, UMI) pairs at the end of each contig
    worker. A UMI with reads on two contigs therefore counts as two UMIs.

    ``n_umi_counted`` counts alignments with the Cell Ranger ``xf`` tag bit 8
    (the read that represents a counted molecule); ``n_xf`` counts alignments
    with an ``xf`` tag.
    """

    region_source: str = "none"
    polya_tso_source: str = "seq"
    region: list[int] = field(default_factory=lambda: [0, 0, 0, 0])
    n_umi: int = 0
    n_umi_reads: int = 0
    n_umi_singleton: int = 0
    n_xf: int = 0
    n_umi_counted: int = 0
    n_artefact_checked: int = 0
    n_polya: int = 0
    n_tso: int = 0
    n_mapq_unique: int = 0
    n_nm: int = 0
    n_nm0: int = 0
    n_nm_ge2: int = 0

    def merge(self, other: ExtraStats) -> ExtraStats:
        return ExtraStats(
            region_source=self.region_source,
            polya_tso_source=self.polya_tso_source,
            region=[a + b for a, b in zip(self.region, other.region)],
            n_umi=self.n_umi + other.n_umi,
            n_umi_reads=self.n_umi_reads + other.n_umi_reads,
            n_umi_singleton=self.n_umi_singleton + other.n_umi_singleton,
            n_xf=self.n_xf + other.n_xf,
            n_umi_counted=self.n_umi_counted + other.n_umi_counted,
            n_artefact_checked=self.n_artefact_checked + other.n_artefact_checked,
            n_polya=self.n_polya + other.n_polya,
            n_tso=self.n_tso + other.n_tso,
            n_mapq_unique=self.n_mapq_unique + other.n_mapq_unique,
            n_nm=self.n_nm + other.n_nm,
            n_nm0=self.n_nm0 + other.n_nm0,
            n_nm_ge2=self.n_nm_ge2 + other.n_nm_ge2,
        )

    def add_umis(self, n_umi: int, n_reads: int, n_singleton: int) -> None:
        """Add the UMI counts of one barcode on one contig."""
        self.n_umi += n_umi
        self.n_umi_reads += n_reads
        self.n_umi_singleton += n_singleton

    def to_dict(self, n_reads: int) -> dict:
        n_region = sum(self.region) if self.region_source != "none" else 0
        ex, intr, inter, anti = self.region
        return {
            "region_source": self.region_source,
            "frac_exonic": _frac(ex, n_region),
            "frac_intronic": _frac(intr, n_region),
            "frac_intergenic": _frac(inter, n_region),
            "frac_antisense": _frac(anti, n_region),
            "n_umi": self.n_umi,
            "mean_reads_per_umi": _frac(self.n_umi_reads, self.n_umi, 3),
            "frac_singleton_umi": _frac(self.n_umi_singleton, self.n_umi),
            "n_umi_counted": self.n_umi_counted,
            "polya_tso_source": self.polya_tso_source,
            "frac_polya": _frac(self.n_polya, self.n_artefact_checked),
            "frac_tso": _frac(self.n_tso, self.n_artefact_checked),
            "frac_mapq_unique": _frac(self.n_mapq_unique, n_reads),
            "frac_nm0": _frac(self.n_nm0, self.n_nm),
            "frac_nm_ge2": _frac(self.n_nm_ge2, self.n_nm),
        }


@dataclass
class CellStats:
    """Running statistics for one (cell, contig) pair.

    ``extra`` is ``None`` unless the extra metrics are on.
    """

    n_reads: int = 0
    n_clipped: int = 0
    n_clip5: int = 0
    n_clip3: int = 0
    n_multimapper: int = 0
    mapq: WelfordAcc = field(default_factory=WelfordAcc)
    base_qual: WelfordAcc = field(default_factory=WelfordAcc)
    clip_frac: WelfordAcc = field(default_factory=WelfordAcc)
    nm: WelfordAcc = field(default_factory=WelfordAcc)
    extra: ExtraStats | None = None

    def update(
        self,
        mapq_val: int,
        qual_val: float,
        clip5: int,
        clip3: int,
        clip_frac_val: float,
        nh_val: float,
        nm_val: float,
    ) -> None:
        self.n_reads += 1
        if clip5 + clip3 > 0:
            self.n_clipped += 1
        if clip5 > 0:
            self.n_clip5 += 1
        if clip3 > 0:
            self.n_clip3 += 1
        if np.isfinite(nh_val) and nh_val > 1:
            self.n_multimapper += 1
        self.mapq.update(mapq_val)
        self.base_qual.update(qual_val)
        self.clip_frac.update(clip_frac_val)
        self.nm.update(nm_val)

    def merge(self, other: CellStats) -> CellStats:
        return CellStats(
            n_reads=self.n_reads + other.n_reads,
            n_clipped=self.n_clipped + other.n_clipped,
            n_clip5=self.n_clip5 + other.n_clip5,
            n_clip3=self.n_clip3 + other.n_clip3,
            n_multimapper=self.n_multimapper + other.n_multimapper,
            mapq=self.mapq.merge(other.mapq),
            base_qual=self.base_qual.merge(other.base_qual),
            clip_frac=self.clip_frac.merge(other.clip_frac),
            nm=self.nm.merge(other.nm),
            extra=_merge_extra(self.extra, other.extra),
        )

    def to_dict(self, cb: str, contig: str | None = None) -> dict:
        nr = self.n_reads
        d = {
            "cb": cb,
            "n_reads": nr,
            "mean_mapq": round(self.mapq.mean, 3),
            "std_mapq": round(self.mapq.std, 3),
            "mean_base_qual": round(self.base_qual.mean, 3),
            "std_base_qual": round(self.base_qual.std, 3),
            "mean_clip_frac": round(self.clip_frac.mean, 4),
            "std_clip_frac": round(self.clip_frac.std, 4),
            "frac_clipped": round(self.n_clipped / nr, 4) if nr else 0.0,
            "frac_clip5": round(self.n_clip5 / nr, 4) if nr else 0.0,
            "frac_clip3": round(self.n_clip3 / nr, 4) if nr else 0.0,
            "frac_multimapper": round(self.n_multimapper / nr, 4) if nr else 0.0,
            "mean_nm": round(self.nm.mean, 3),
        }
        if contig is not None:
            d["contig"] = contig
        if self.extra is not None:
            d.update(self.extra.to_dict(nr))
        return d


def _merge_extra(a: ExtraStats | None, b: ExtraStats | None) -> ExtraStats | None:
    if a is None:
        return b
    if b is None:
        return a
    return a.merge(b)


# ---------------------------------------------------------------------------
# Per-contig worker
# ---------------------------------------------------------------------------

# Module-level references so worker processes can reach the shared objects
_shared_counter = None
_shared_lock = None
_shared_barcode_tag: str = "CB"
_shared_nm_tag: str = "nM"
_shared_extra: ExtraConfig | None = None


def _worker_init(
    counter,
    lock,
    barcode_tag: str = "CB",
    nm_tag: str = "nM",
    extra: ExtraConfig | None = None,
) -> None:
    global _shared_counter, _shared_lock, _shared_barcode_tag, _shared_nm_tag
    global _shared_extra
    _shared_counter = counter
    _shared_lock = lock
    _shared_barcode_tag = barcode_tag
    _shared_nm_tag = nm_tag
    _shared_extra = extra


def _region_from_tags(read) -> int | None:
    """Region class from Cell Ranger tags, or ``None`` without an ``RE`` tag.

    Assumed Cell Ranger semantics (BAM tag documentation, Cell Ranger 7+):

    - ``RE:A``: region of the alignment, ``E`` exonic, ``N`` intronic,
      ``I`` intergenic. ``RE`` does not depend on the strand.
    - ``TX:Z`` / ``GX:Z``: transcripts / genes on the same strand as the read
      (sense) that are compatible with the alignment.
    - ``AN:Z``: like ``TX``, but for transcripts (or, for intronic reads, genes)
      on the opposite strand.

    A read is antisense when it has ``AN`` and neither ``GX`` nor ``TX``.
    Otherwise ``RE`` gives the class. So the four classes do not overlap.
    """
    if not read.has_tag("RE"):
        return None
    if read.has_tag("AN") and not (read.has_tag("GX") or read.has_tag("TX")):
        return ANTISENSE
    re_val = read.get_tag("RE")
    if re_val == "E":
        return EXONIC
    if re_val == "N":
        return INTRONIC
    if re_val == "I":
        return INTERGENIC
    return None


def _transcript_plus(read, strand: str) -> bool | None:
    """Genome strand of the transcript that the read comes from.

    ``"forward"``: the read (a single read such as 10x 3' R2, or read 1 of a
    pair) is on the transcript strand. Read 2 of a pair is on the opposite
    strand. ``"reverse"`` flips the result. ``"unstranded"`` gives ``None``.
    """
    if strand == "unstranded":
        return None
    plus = not read.is_reverse
    if read.is_paired and read.is_read2:
        plus = not plus
    if strand == "reverse":
        plus = not plus
    return plus


def _update_extra(
    ex: ExtraStats,
    read,
    cb: str,
    cfg: ExtraConfig,
    contig_annot: ContigAnnot | None,
    cigar: str,
    clip5: int,
    clip3: int,
    nm_val: float,
    umis: UmiCounter,
) -> None:
    """Add one alignment to the extra-metric counters (``extra_metrics=True``)."""
    if read.mapping_quality >= cfg.mapq_unique:
        ex.n_mapq_unique += 1
    if nm_val == nm_val:  # not NaN: the read has the mismatch tag
        ex.n_nm += 1
        if nm_val == 0:
            ex.n_nm0 += 1
        elif nm_val >= 2:
            ex.n_nm_ge2 += 1

    # Read region
    if cfg.region_source == "tags":
        cls = _region_from_tags(read)
        if cls is not None:
            ex.region[cls] += 1
    elif cfg.region_source == "gtf":
        if contig_annot is None:
            ex.region[INTERGENIC] += 1  # contig not in the GTF
        else:
            plus = _transcript_plus(read, cfg.strand)
            sense = None if plus is None else ("+" if plus else "-")
            ex.region[contig_annot.classify(read.get_blocks(), sense)] += 1

    # UMI: primary alignments only, so that a multimapper counts once
    if (
        not read.is_secondary
        and not read.is_supplementary
        and read.has_tag(cfg.umi_tag)
    ):
        umis.add(cb, read.get_tag(cfg.umi_tag))
    # Cell Ranger xf:i bit 8 (value 8): this read represents a molecule that
    # is counted in the feature-barcode matrix.
    if read.has_tag("xf"):
        ex.n_xf += 1
        if read.get_tag("xf") & 8:
            ex.n_umi_counted += 1

    # PolyA and TSO
    if cfg.artefact_source == "tags":
        # Cell Ranger: pa:i = polyA bases trimmed from the 3' end of R2,
        # ts:i = TSO bases trimmed from the 5' end of R2. A read without the
        # tag had no bases trimmed.
        ex.n_artefact_checked += 1
        if read.has_tag("pa") and read.get_tag("pa") > 0:
            ex.n_polya += 1
        if read.has_tag("ts") and read.get_tag("ts") > 0:
            ex.n_tso += 1
    else:
        seq = read.query_sequence
        if seq is not None:  # secondary alignments can have SEQ "*"
            ex.n_artefact_checked += 1
            left, right = _scan_soft_clips(read, cigar, clip5, clip3)
            if left >= POLYA_MIN_RUN or right >= POLYA_MIN_RUN:
                polya, tso = scan_clips(
                    seq, left, right, _transcript_plus(read, cfg.strand)
                )
                ex.n_polya += polya
                ex.n_tso += tso


def _worker(
    args: tuple[str, str] | tuple[str, str, ContigAnnot | None],
) -> dict[StatsKey, CellStats]:
    """Stream all reads of one contig, keyed by (cell_id, contig).

    ``args`` is ``(bam_path, contig)``, plus the contig annotation when the
    extra metrics use a GTF.
    """
    pysam = _import_pysam()
    bam_path, contig = args[0], args[1]
    contig_annot = args[2] if len(args) > 2 else None
    cfg = _shared_extra
    umi_counter: UmiCounter | None = None
    if cfg is None:
        stats: dict[StatsKey, CellStats] = defaultdict(CellStats)
    else:
        src, art = cfg.region_source, cfg.artefact_source
        stats = defaultdict(
            lambda: CellStats(extra=ExtraStats(region_source=src, polya_tso_source=art))
        )
        # (barcode, UMI) -> number of reads, for this contig only
        umi_counter = UmiCounter()
    local_count = 0
    tag = _shared_barcode_tag
    nm_tag = _shared_nm_tag
    counter, lock = _shared_counter, _shared_lock

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(contig):
            if read.is_unmapped or read.cigarstring is None:
                continue

            cb: str = read.get_tag(tag) if read.has_tag(tag) else NO_BARCODE  # type: ignore[assignment]

            cigar = read.cigarstring
            clip5, clip3 = parse_soft_clips(cigar)
            rlen = read_length_from_cigar(cigar)
            clip_frac = (clip5 + clip3) / rlen if rlen > 0 else 0.0

            quals = read.query_qualities
            mean_qual = (
                float(np.mean(quals))
                if quals is not None and len(quals) > 0
                else float("nan")
            )

            nh_val = float(read.get_tag("NH")) if read.has_tag("NH") else float("nan")
            nm_val = (
                float(read.get_tag(nm_tag)) if read.has_tag(nm_tag) else float("nan")
            )

            cs = stats[(cb, contig)]
            cs.update(
                mapq_val=read.mapping_quality,
                qual_val=mean_qual,
                clip5=clip5,
                clip3=clip3,
                clip_frac_val=clip_frac,
                nh_val=nh_val,
                nm_val=nm_val,
            )
            if umi_counter is not None and cfg is not None and cs.extra is not None:
                _update_extra(
                    cs.extra,
                    read,
                    cb,
                    cfg,
                    contig_annot,
                    cigar,
                    clip5,
                    clip3,
                    nm_val,
                    umi_counter,
                )

            local_count += 1
            if local_count % _REPORT_EVERY == 0 and counter is not None:
                if lock is not None:
                    with lock:
                        counter.value += _REPORT_EVERY

    remainder = local_count % _REPORT_EVERY
    if remainder and counter is not None and lock is not None:
        with lock:
            counter.value += remainder

    if umi_counter is not None:
        # Reduce (barcode, UMI) counts to per-barcode sums before returning,
        # so that only O(barcodes) data leaves the worker.
        for cb, counts in umi_counter.per_barcode().items():
            ex = stats[(cb, contig)].extra
            if ex is not None:
                ex.add_umis(*counts)
        del umi_counter

    return dict(stats)


# ---------------------------------------------------------------------------
# Progress display
# ---------------------------------------------------------------------------


def _make_progress(console: Console) -> Progress:
    from rich.progress import (
        BarColumn,
        MofNCompleteColumn,
        Progress,
        SpinnerColumn,
        TaskProgressColumn,
        TextColumn,
        TimeElapsedColumn,
        TimeRemainingColumn,
    )

    return Progress(
        SpinnerColumn(spinner_name="dots", style="bold cyan"),
        TextColumn("[bold]{task.description}"),
        BarColumn(
            bar_width=36,
            style="cyan",
            complete_style="bold cyan",
            finished_style="bold green",
        ),
        TaskProgressColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
        transient=False,
    )


def _poll_progress(counter, lock, progress, task_id, stop_evt, interval=0.25) -> None:
    """Background thread: poll the shared counter and advance the progress bar."""
    prev = 0
    while not stop_evt.is_set():
        with lock:
            cur = counter.value
        if cur != prev:
            progress.advance(task_id, cur - prev)
            prev = cur
        stop_evt.wait(interval)
    with lock:
        cur = counter.value
    if cur != prev:
        progress.advance(task_id, cur - prev)


# ---------------------------------------------------------------------------
# Parallel collector
# ---------------------------------------------------------------------------


def collect_metrics(
    bam_path: str,
    workers: int = 4,
    barcode_tag: str = "CB",
    nm_tag: str = "nM",
    contig_list: list[str] | None = None,
    contig_prefixes: list[str] | None = None,
    show_progress: bool = False,
    console: Console | None = None,
    extra_metrics: bool = False,
    gtf: str | PathLike | None = None,
    strand: Strand = "forward",
    umi_tag: str = "UB",
    mapq_unique: int = 255,
    region_source: Literal["auto", "tags", "gtf"] = "auto",
) -> dict[StatsKey, CellStats]:
    """Shard an indexed BAM by contig and collect per-(cell, contig) metrics.

    Parameters
    ----------
    bam_path
        Path to a coordinate-sorted, indexed BAM file.
    workers
        Number of worker processes.
    barcode_tag
        BAM tag holding the cell identity (``CB`` for Cell Ranger, ``RG`` for
        STARsolo SmartSeq). Reads without it are counted under ``NO_BARCODE``.
    nm_tag
        BAM tag holding the mismatch count, averaged into ``mean_nm``. ``nM``
        (STAR, Cell Ranger; mismatches per read pair) by default, ``NM`` for
        aligners that only write the standard edit distance. Reads without it
        do not contribute to ``mean_nm`` (``NaN`` when no read has it).
    contig_list, contig_prefixes
        Restrict to these contigs (union of both filters).
    show_progress
        Show a live progress bar (requires ``console``).
    console
        ``rich`` console for status messages. ``None`` runs silently.
    extra_metrics
        Also count the metrics of :class:`ExtraStats` (columns
        ``EXTRA_COLUMNS``, appended after the default columns). Off by
        default; then the output is the same as ``bam_qc_summary.py`` v0.2.1.
    gtf
        GTF annotation for the read regions. Used only when the BAM has no
        Cell Ranger ``RE`` tag. Without ``RE`` tags and without a GTF, the
        region fractions are ``NaN``.
    strand
        Library strand, for the GTF regions and the polyA/TSO sequence scan.
        ``"forward"``: the read (read 1 of a pair) is on the transcript strand,
        as 10x 3' R2; read 2 of a pair is on the opposite strand.
        ``"reverse"``: the opposite. ``"unstranded"`` (e.g. Smart-seq2): both
        strands are sense, so no read is antisense.
    umi_tag
        BAM tag with the corrected UMI (``UB`` for Cell Ranger and STARsolo).
        Without it on any read, the UMI columns are ``NaN``.
    mapq_unique
        Reads with ``MAPQ >= mapq_unique`` count as unique in
        ``frac_mapq_unique``. STAR and Cell Ranger give 255 to unique reads.
    region_source
        Source of the read regions. ``"auto"``: the ``RE`` tags if found, else
        the GTF if given, else none (``NaN``). ``"tags"``: the ``RE``/``AN``
        tags. ``"gtf"``: the GTF, also when the BAM has ``RE`` tags (to compare
        the two sources).

    The ``RE`` tag (regions) and the ``pa``/``ts`` tags (polyA/TSO) are
    searched in the first 20,000 mapped reads of the selected contigs. If
    found, the tags are used for the whole BAM. Otherwise the GTF (regions)
    and the soft-clipped sequence (polyA/TSO) are used.

    Raises
    ------
    ValueError
        If the BAM index is empty or missing, or no contig remains after
        filtering, or ``strand`` or ``region_source`` is not valid, or
        ``region_source="tags"`` and the sample has no ``RE`` tag, or
        ``region_source="gtf"`` without ``gtf``.
    """
    pysam = _import_pysam()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        index_stats = bam.get_index_statistics()
        contigs = [s.contig for s in index_stats if s.total > 0]
        read_counts = {s.contig: s.total for s in index_stats if s.total > 0}

    if not contigs:
        raise ValueError(
            f"BAM index of {bam_path} is empty or missing. Run `samtools index` first."
        )

    contigs_before = len(contigs)
    contigs = filter_contigs(
        contigs, contig_list=contig_list, contig_prefixes=contig_prefixes
    )
    if not contigs:
        raise ValueError(
            "No contigs remain after filtering. Check the contig list and prefixes."
        )

    total_reads = sum(read_counts[c] for c in contigs)
    tasks: list[tuple] = [(bam_path, c) for c in contigs]
    counter = mp.Value("q", 0)  # signed 64-bit int
    lock = mp.Lock()

    if console is not None:
        filter_note = ""
        if contigs_before != len(contigs):
            filter_note = (
                f"  [dim](filtered {contigs_before} → {len(contigs)} contigs)[/dim]"
            )
        console.print(
            f"[bold cyan]scbamqc[/bold cyan]  "
            f"[dim]{len(tasks)} contigs[/dim]  "
            f"[dim]{workers} workers[/dim]  "
            f"[dim]~{total_reads:,} reads[/dim]  "
            f"[dim]tag={barcode_tag}[/dim]  "
            f"[dim]nm-tag={nm_tag}[/dim]"
            f"{filter_note}"
        )

    extra_cfg: ExtraConfig | None = None
    if extra_metrics:
        extra_cfg, annot = _extra_config(
            bam_path,
            gtf,
            strand,
            umi_tag,
            mapq_unique,
            contigs,
            console,
            region_source,
        )
        if annot is not None:
            tasks = [(bam_path, c, annot.get(c)) for c in contigs]

    def _run_pool() -> list[dict[StatsKey, CellStats]]:
        with mp.Pool(
            processes=workers,
            initializer=_worker_init,
            initargs=(counter, lock, barcode_tag, nm_tag, extra_cfg),
        ) as pool:
            return pool.map(_worker, tasks)

    if show_progress and console is not None:
        progress = _make_progress(console)
        task_id = progress.add_task("processing reads", total=total_reads)
        stop_evt = threading.Event()
        poll_thread = threading.Thread(
            target=_poll_progress,
            args=(counter, lock, progress, task_id, stop_evt),
            daemon=True,
        )
        with progress:
            poll_thread.start()
            try:
                shard_results = _run_pool()
            finally:
                stop_evt.set()
                poll_thread.join()
    else:
        shard_results = _run_pool()

    merged: dict[StatsKey, CellStats] = {}
    for shard in shard_results:
        for key, cs in shard.items():
            merged[key] = merged[key].merge(cs) if key in merged else cs

    if console is not None:
        n_cells = len({cb for cb, _ in merged if cb != NO_BARCODE})
        n_contigs = len({ctg for _, ctg in merged})
        n_reads = sum(cs.n_reads for cs in merged.values())
        console.print(
            f"[bold green]✓[/bold green] "
            f"[bold]{n_reads:,}[/bold] mapped reads  "
            f"[bold]{n_cells:,}[/bold] cell IDs  "
            f"[bold]{n_contigs:,}[/bold] contigs"
        )
        if extra_cfg is not None and not _has_umis(merged):
            console.print(
                f"[yellow]Note:[/yellow] No read has the UMI tag "
                f"[bold]{umi_tag}[/bold]. The UMI columns are NaN."
            )
    return merged


def _extra_config(
    bam_path: str,
    gtf: str | PathLike | None,
    strand: Strand,
    umi_tag: str,
    mapq_unique: int,
    contigs: list[str],
    console: Console | None,
    source: str = "auto",
) -> tuple[ExtraConfig, dict[str, ContigAnnot] | None]:
    """Choose the region and polyA/TSO sources and load the GTF if needed."""
    if strand not in ("forward", "reverse", "unstranded"):
        raise ValueError(
            f"strand must be 'forward', 'reverse' or 'unstranded', got {strand!r}"
        )
    if source not in ("auto", "tags", "gtf"):
        raise ValueError(
            f"region_source must be 'auto', 'tags' or 'gtf', got {source!r}"
        )
    if source == "gtf" and gtf is None:
        raise ValueError("region_source 'gtf' needs a GTF file (--gtf).")
    n_sampled, tag_counts = count_tags(bam_path, n_reads=_TAG_SAMPLE, contigs=contigs)
    has_re = bool(tag_counts.get("RE"))
    if source == "tags" and not has_re:
        raise ValueError(
            f"region_source 'tags': no RE tag in the first {n_sampled:,} mapped "
            "reads. Use --gtf with --region-source auto or gtf."
        )
    annot = None
    if has_re and source != "gtf":
        region_source = "tags"
        region_note = "RE/AN/GX/TX tags"
        if gtf is not None:
            region_note += " (the GTF is not used)"
    elif gtf is not None:
        region_source = "gtf"
        annot = load_gtf(gtf)
        n_match = sum(c in annot for c in contigs)
        region_note = (
            f"GTF, strand={strand}, {n_match}/{len(contigs)} contigs in the GTF"
        )
        if has_re:
            region_note += " (the RE tags are not used)"
        if n_match == 0 and console is not None:
            console.print(
                "[yellow]Warning:[/yellow] No BAM contig is in the GTF. "
                "All reads are intergenic. Check the contig names."
            )
    else:
        region_source = "none"
        region_note = "none (no RE tag, no GTF): NaN"
    has_art_tags = bool(tag_counts.get("pa") or tag_counts.get("ts"))
    artefact_source = "tags" if has_art_tags else "seq"
    cfg = ExtraConfig(
        region_source=region_source,
        artefact_source=artefact_source,
        strand=strand,
        umi_tag=umi_tag,
        mapq_unique=mapq_unique,
    )
    if console is not None:
        art_note = "pa/ts tags" if has_art_tags else f"soft clips, strand={strand}"
        console.print(
            f"[bold cyan]extra metrics[/bold cyan]  "
            f"[dim]regions: {region_note}[/dim]  "
            f"[dim]polyA/TSO: {art_note}[/dim]  "
            f"[dim]umi-tag={umi_tag}[/dim]  "
            f"[dim]mapq-unique>={mapq_unique}[/dim]"
        )
    return cfg, annot


def _has_umis(stats: dict[StatsKey, CellStats]) -> bool:
    return any(cs.extra is not None and cs.extra.n_umi > 0 for cs in stats.values())


# ---------------------------------------------------------------------------
# Tag inspection
# ---------------------------------------------------------------------------


def count_tags(
    bam_path: str, n_reads: int = 100_000, contigs: list[str] | None = None
) -> tuple[int, dict[str, int]]:
    """Count how many of the first ``n_reads`` mapped reads carry each BAM tag.

    Use it to find the cell-identity and mismatch tags of a BAM from an
    unfamiliar pipeline. Reads are taken in file order, so a coordinate-sorted
    BAM is sampled from its first contigs. With ``contigs`` (needs an index),
    only reads of these contigs are sampled, in the given order.

    Returns
    -------
    ``(n_sampled, {tag: n_reads_with_tag})``, tags sorted by frequency.
    """
    pysam = _import_pysam()
    counts: dict[str, int] = defaultdict(int)
    n = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        if contigs is None:
            reads = bam.fetch(until_eof=True)
        else:
            reads = (r for c in contigs for r in bam.fetch(c))
        for read in reads:
            if read.is_unmapped:
                continue
            for tag, _ in read.get_tags():
                counts[tag] += 1
            n += 1
            if n >= n_reads:
                break
    return n, dict(sorted(counts.items(), key=lambda kv: -kv[1]))


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------


def collapse_to_cell(stats: dict[StatsKey, CellStats]) -> dict[str, CellStats]:
    """Merge all (cell, contig) entries into one ``CellStats`` per cell."""
    by_cell: dict[str, CellStats] = {}
    for (cb, _contig), cs in stats.items():
        by_cell[cb] = by_cell[cb].merge(cs) if cb in by_cell else cs
    return by_cell


def to_dataframe(
    stats: dict[StatsKey, CellStats],
    include_no_cb: bool = False,
    group_by: GroupBy = "cell",
) -> pd.DataFrame:
    """Build the per-barcode table from the output of :func:`collect_metrics`.

    Parameters
    ----------
    stats
        Mapping of (cell, contig) to ``CellStats``.
    include_no_cb
        Keep the row of reads that lacked the barcode tag.
    group_by
        ``"cell"`` (one row per cell, contigs collapsed) or ``"cell-contig"``
        (one row per cell x contig).
    """
    if group_by == "cell":
        rows = [
            cs.to_dict(cb)
            for cb, cs in collapse_to_cell(stats).items()
            if include_no_cb or cb != NO_BARCODE
        ]
    elif group_by == "cell-contig":
        rows = [
            cs.to_dict(cb, contig)
            for (cb, contig), cs in stats.items()
            if include_no_cb or cb != NO_BARCODE
        ]
    else:
        raise ValueError(f"group_by must be 'cell' or 'cell-contig', got {group_by!r}")

    has_extra = any(cs.extra is not None for cs in stats.values())
    if not rows:
        empty = CellStats(extra=ExtraStats() if has_extra else None)
        cols = list(empty.to_dict("", "" if group_by == "cell-contig" else None))
        return pd.DataFrame(columns=cols)
    df = (
        pd.DataFrame(rows)
        .sort_values(["cb", "n_reads"], ascending=[True, False])
        .reset_index(drop=True)
    )
    if has_extra and not _has_umis(stats):
        # no UMI tag in the BAM (e.g. Smart-seq2): n_umi is unknown, not 0
        df[list(UMI_COLUMNS)] = np.nan
    if has_extra and not any(
        cs.extra is not None and cs.extra.n_xf > 0 for cs in stats.values()
    ):
        df["n_umi_counted"] = np.nan  # no Cell Ranger xf tag
    return df
