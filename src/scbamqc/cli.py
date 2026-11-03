"""Command-line interface: ``scbamqc bam`` and ``scbamqc tags``."""

from __future__ import annotations

import multiprocessing as mp
import sys
from pathlib import Path
from typing import Annotated, Literal

import cyclopts
from rich.console import Console
from rich.table import Table

from . import __version__

console = Console()
err_console = Console(stderr=True)

app = cyclopts.App(
    name="scbamqc",
    help="Alignment QC for single-cell BAM files.",
    version=__version__,
)

# Tags that commonly carry cell identity or mismatches, by pipeline.
_BARCODE_TAGS = ("CB", "CR", "XC", "RG", "BC", "DB")
_NM_TAGS = ("nM", "NM")
_TAG_SAMPLE = 20_000


def _fail(msg: str) -> None:
    err_console.print(f"[bold red]Error:[/bold red] {msg}")
    sys.exit(1)


def _check_tags(bam: Path, barcode_tag: str, nm_tag: str) -> None:
    """Warn when the chosen tags are missing from most sampled reads."""
    from .bam import count_tags

    n, counts = count_tags(str(bam), n_reads=_TAG_SAMPLE)
    if n == 0:
        return
    for tag, kind, candidates, flag in [
        (barcode_tag, "barcode", _BARCODE_TAGS, "--barcode-tag"),
        (nm_tag, "mismatch", _NM_TAGS, "--nm-tag"),
    ]:
        frac = counts.get(tag, 0) / n
        if frac < 0.5:
            present = [f"{t} ({counts[t] / n:.0%})" for t in candidates if t in counts]
            hint = f" Tags found: {', '.join(present)}." if present else ""
            console.print(
                f"[yellow]Warning:[/yellow] Only {frac:.0%} of the first {n:,} "
                f"mapped reads have the {kind} tag [bold]{tag}[/bold].{hint} "
                f"Set the tag with [bold]{flag}[/bold]. "
                "To see all tags, run [bold]scbamqc tags[/bold]."
            )


@app.command
def bam(
    bam: Annotated[
        Path,
        cyclopts.Parameter(help="BAM file, sorted by coordinate and indexed."),
    ],
    *,
    workers: Annotated[
        int,
        cyclopts.Parameter(
            name=["-j", "--workers"],
            help="Number of worker processes. 0 uses half of the CPU cores.",
        ),
    ] = 0,
    barcode_tag: Annotated[
        str,
        cyclopts.Parameter(
            name="--barcode-tag",
            help="BAM tag with the cell barcode. Cell Ranger: CB. STARsolo "
            "SmartSeq: RG. For other pipelines, run `scbamqc tags`.",
        ),
    ] = "CB",
    nm_tag: Annotated[
        str,
        cyclopts.Parameter(
            name="--nm-tag",
            help="BAM tag with the mismatch count, for mean_nm. STAR and Cell "
            "Ranger: nM. Other aligners: NM.",
        ),
    ] = "nM",
    group_by: Annotated[
        Literal["cell", "cell-contig"],
        cyclopts.Parameter(
            name="--group-by",
            help="cell: one row for each cell. "
            "cell-contig: one row for each cell and contig.",
        ),
    ] = "cell",
    contigs: Annotated[
        list[str] | None,
        cyclopts.Parameter(
            name="--contigs",
            help="Use only these contigs, for example: chr1 chr2. "
            "Adds to --contig-prefix.",
            negative_iterable="",
        ),
    ] = None,
    contig_prefix: Annotated[
        list[str] | None,
        cyclopts.Parameter(
            name="--contig-prefix",
            help="Use only contigs with names that start with these prefixes. "
            "Adds to --contigs.",
            negative_iterable="",
        ),
    ] = None,
    out_prefix: Annotated[
        str,
        cyclopts.Parameter(name="--out-prefix", help="Prefix of the output files."),
    ] = "bam_qc",
    include_no_cb: Annotated[
        bool,
        cyclopts.Parameter(
            name="--include-no-cb",
            help="Add one row for the reads without the barcode tag.",
            negative_bool="",
        ),
    ] = False,
    show_progress: Annotated[
        bool,
        cyclopts.Parameter(
            name="--progress",
            help="Show the progress bar. Off when the output is not a terminal.",
        ),
    ] = True,
    make_plots: Annotated[
        bool,
        cyclopts.Parameter(
            name="--plots",
            help="Also write histograms (PNG). Needs matplotlib.",
            negative_bool="",
        ),
    ] = False,
    write_csv: Annotated[
        bool,
        cyclopts.Parameter(name="--csv", help="Write the table (CSV)."),
    ] = True,
    check_tags: Annotated[
        bool,
        cyclopts.Parameter(
            name="--check-tags",
            help="Check the tags on the first reads. Warn if they are missing.",
        ),
    ] = True,
):
    """Calculate alignment statistics for each cell barcode.

    Writes OUT_PREFIX_per_cell.csv.gz, and OUT_PREFIX_plots.png with --plots.
    Memory use depends on the number of barcodes, not on the number of reads.
    """
    from ._report import print_summary
    from .bam import collect_metrics, to_dataframe

    if not bam.exists():
        _fail(f"File not found: {bam}")

    if not any(p.exists() for p in (bam.with_suffix(".bai"), Path(f"{bam}.bai"))):
        console.print(
            f"[yellow]Warning:[/yellow] No .bai index for [cyan]{bam.name}[/cyan]. "
            "Run [bold]samtools index[/bold]."
        )

    if make_plots:
        try:
            import matplotlib

            matplotlib.use("Agg")
        except ImportError:
            _fail("--plots needs matplotlib. Run: pip install scbamqc[plot]")

    if check_tags:
        _check_tags(bam, barcode_tag, nm_tag)

    n_workers = workers if workers > 0 else max(1, mp.cpu_count() // 2)
    try:
        stats = collect_metrics(
            str(bam),
            workers=n_workers,
            barcode_tag=barcode_tag,
            nm_tag=nm_tag,
            contig_list=contigs,
            contig_prefixes=contig_prefix,
            show_progress=show_progress and console.is_terminal,
            console=console,
        )
    except ValueError as e:
        _fail(str(e))

    df = to_dataframe(stats, include_no_cb=include_no_cb, group_by=group_by)
    print_summary(console, df, stats, barcode_tag=barcode_tag, group_by=group_by)

    if make_plots and not df.empty:
        import matplotlib.pyplot as plt

        from .pl import bam_qc_hist

        fig = bam_qc_hist(df, barcode_tag=barcode_tag, group_by=group_by)
        plot_path = f"{out_prefix}_plots.png"
        fig.savefig(plot_path, dpi=150)
        plt.close(fig)
        console.print(f"[bold green]✓[/bold green] plots → [cyan]{plot_path}[/cyan]")

    if write_csv:
        csv_path = f"{out_prefix}_per_cell.csv.gz"
        df.to_csv(csv_path, index=False, compression="gzip")
        console.print(
            f"[bold green]✓[/bold green] per-cell CSV ({len(df):,} rows) "
            f"→ [cyan]{csv_path}[/cyan]"
        )


@app.command
def tags(
    bam: Annotated[Path, cyclopts.Parameter(help="Path to a BAM file.")],
    *,
    n_reads: Annotated[
        int,
        cyclopts.Parameter(
            name=["-n", "--n-reads"], help="Number of mapped reads to examine."
        ),
    ] = 100_000,
):
    """List the BAM tags on the first mapped reads.

    Use this to find --barcode-tag and --nm-tag for a new pipeline.
    """
    from .bam import count_tags

    if not bam.exists():
        _fail(f"File not found: {bam}")
    n, counts = count_tags(str(bam), n_reads=n_reads)
    if n == 0:
        _fail("The BAM file has no mapped reads.")

    t = Table(title=f"Tags on the first {n:,} mapped reads", header_style="bold cyan")
    t.add_column("tag", style="bold")
    t.add_column("reads", justify="right")
    t.add_column("fraction", justify="right")
    t.add_column("note", style="dim")
    notes = {
        **{k: "possible cell barcode" for k in _BARCODE_TAGS},
        "nM": "mismatches in the read pair (STAR)",
        "NM": "edit distance to the reference",
        "NH": "number of alignments of the read",
    }
    for tag, c in counts.items():
        t.add_row(tag, f"{c:,}", f"{c / n:.1%}", notes.get(tag, ""))
    console.print(t)


def main() -> None:
    app()


if __name__ == "__main__":
    main()
