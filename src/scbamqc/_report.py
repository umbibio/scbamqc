"""Terminal summary of a BAM-QC run (``rich``)."""

from __future__ import annotations

import pandas as pd
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from .bam import NO_BARCODE, CellStats, StatsKey
from .pl import _id_label


def _pct(part: int, total: int) -> str:
    return f"{100 * part / total:.1f}%" if total else "n/a"


def print_summary(
    console: Console,
    df: pd.DataFrame,
    stats: dict[StatsKey, CellStats],
    barcode_tag: str = "CB",
    group_by: str = "cell",
) -> None:
    total_reads = sum(cs.n_reads for cs in stats.values())
    no_cb_reads = sum(cs.n_reads for (cb, _), cs in stats.items() if cb == NO_BARCODE)
    cb_reads = total_reads - no_cb_reads

    hdr = Table.grid(padding=(0, 2))
    hdr.add_column(style="bold cyan")
    hdr.add_column()
    hdr.add_row("total mapped reads", f"{total_reads:>14,}")
    hdr.add_row(
        f"reads with {barcode_tag} tag",
        f"{cb_reads:>14,}  ({_pct(cb_reads, total_reads)})",
    )
    hdr.add_row(
        f"reads without {barcode_tag} tag",
        f"{no_cb_reads:>14,}  ({_pct(no_cb_reads, total_reads)})",
    )
    hdr.add_row(f"unique {_id_label(barcode_tag)}", f"{df['cb'].nunique():>14,}")
    if group_by == "cell-contig":
        hdr.add_row("contigs with reads", f"{df['contig'].nunique():>14,}")
        hdr.add_row("rows (cell × contig)", f"{len(df):>14,}")
    else:
        hdr.add_row("rows (per cell)", f"{len(df):>14,}")

    console.print(
        Panel(
            hdr,
            title="[bold]BAM QC summary[/bold]",
            border_style="cyan",
            padding=(1, 2),
        )
    )
    if df.empty:
        return

    def _metric_table(title: str, col: str, digits: int) -> Table:
        t = Table(
            title=title,
            title_style="bold",
            header_style="bold cyan",
            border_style="dim",
            padding=(0, 1),
        )
        t.add_column("metric", style="cyan", no_wrap=True)
        t.add_column("value", justify="right")
        for k, v in df[col].describe().round(digits).items():
            t.add_row(str(k), str(v))
        return t

    console.print()
    console.print(
        _metric_table("MAPQ (per-cell means)", "mean_mapq", 2),
        _metric_table("Base quality (per-cell means)", "mean_base_qual", 2),
        _metric_table("Clip fraction (per-cell means)", "mean_clip_frac", 4),
    )
    console.print()
