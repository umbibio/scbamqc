"""Streaming per-barcode BAM QC.

Works with Cell Ranger BAMs (``CB`` tag) and STARsolo SmartSeq BAMs (``RG`` tag).
Workers stream reads one contig at a time and accumulate per-(cell, contig)
statistics with Welford's online algorithm, so no per-read rows are stored.
Peak memory is O(unique barcodes x contigs touched), not O(reads).

``pysam`` is an optional dependency (``pip install scbamqc[bam]``) and is only
imported when a BAM file is read.
"""

from __future__ import annotations

import multiprocessing as mp
import re
import threading
from collections import defaultdict
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from rich.console import Console
    from rich.progress import Progress

__all__ = [
    "NO_BARCODE",
    "CellStats",
    "WelfordAcc",
    "collapse_to_cell",
    "collect_metrics",
    "count_tags",
    "filter_contigs",
    "parse_soft_clips",
    "read_length_from_cigar",
    "to_dataframe",
]

GroupBy = Literal["cell", "cell-contig"]
StatsKey = tuple[str, str]  # (cell_id, contig)

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
NO_BARCODE = "__no_barcode__"
_REPORT_EVERY = 50_000  # a worker reports progress every N processed reads


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


@dataclass
class CellStats:
    """Running statistics for one (cell, contig) pair."""

    n_reads: int = 0
    n_clipped: int = 0
    n_clip5: int = 0
    n_clip3: int = 0
    n_multimapper: int = 0
    mapq: WelfordAcc = field(default_factory=WelfordAcc)
    base_qual: WelfordAcc = field(default_factory=WelfordAcc)
    clip_frac: WelfordAcc = field(default_factory=WelfordAcc)
    nm: WelfordAcc = field(default_factory=WelfordAcc)

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
        return d


# ---------------------------------------------------------------------------
# Per-contig worker
# ---------------------------------------------------------------------------

# Module-level references so worker processes can reach the shared objects
_shared_counter = None
_shared_lock = None
_shared_barcode_tag: str = "CB"
_shared_nm_tag: str = "nM"


def _worker_init(counter, lock, barcode_tag: str = "CB", nm_tag: str = "nM") -> None:
    global _shared_counter, _shared_lock, _shared_barcode_tag, _shared_nm_tag
    _shared_counter = counter
    _shared_lock = lock
    _shared_barcode_tag = barcode_tag
    _shared_nm_tag = nm_tag


def _worker(args: tuple[str, str]) -> dict[StatsKey, CellStats]:
    """Stream all reads of one contig, keyed by (cell_id, contig)."""
    pysam = _import_pysam()
    bam_path, contig = args
    stats: dict[StatsKey, CellStats] = defaultdict(CellStats)
    local_count = 0
    tag = _shared_barcode_tag
    nm_tag = _shared_nm_tag

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(contig):
            if read.is_unmapped or read.cigarstring is None:
                continue

            cb = read.get_tag(tag) if read.has_tag(tag) else NO_BARCODE

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

            stats[(cb, contig)].update(
                mapq_val=read.mapping_quality,
                qual_val=mean_qual,
                clip5=clip5,
                clip3=clip3,
                clip_frac_val=clip_frac,
                nh_val=nh_val,
                nm_val=nm_val,
            )

            local_count += 1
            if local_count % _REPORT_EVERY == 0 and _shared_counter is not None:
                with _shared_lock:
                    _shared_counter.value += _REPORT_EVERY

    remainder = local_count % _REPORT_EVERY
    if remainder and _shared_counter is not None:
        with _shared_lock:
            _shared_counter.value += remainder

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

    Raises
    ------
    ValueError
        If the BAM index is empty or missing, or no contig remains after
        filtering.
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
    tasks = [(bam_path, c) for c in contigs]
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

    def _run_pool() -> list[dict[StatsKey, CellStats]]:
        with mp.Pool(
            processes=workers,
            initializer=_worker_init,
            initargs=(counter, lock, barcode_tag, nm_tag),
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
    return merged


# ---------------------------------------------------------------------------
# Tag inspection
# ---------------------------------------------------------------------------


def count_tags(bam_path: str, n_reads: int = 100_000) -> tuple[int, dict[str, int]]:
    """Count how many of the first ``n_reads`` mapped reads carry each BAM tag.

    Use it to find the cell-identity and mismatch tags of a BAM from an
    unfamiliar pipeline. Reads are taken in file order, so a coordinate-sorted
    BAM is sampled from its first contigs.

    Returns
    -------
    ``(n_sampled, {tag: n_reads_with_tag})``, tags sorted by frequency.
    """
    pysam = _import_pysam()
    counts: dict[str, int] = defaultdict(int)
    n = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
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

    if not rows:
        cols = list(CellStats().to_dict("", "" if group_by == "cell-contig" else None))
        return pd.DataFrame(columns=cols)
    return (
        pd.DataFrame(rows)
        .sort_values(["cb", "n_reads"], ascending=[True, False])
        .reset_index(drop=True)
    )
