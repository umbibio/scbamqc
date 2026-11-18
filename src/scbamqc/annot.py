"""Minimal GTF loader and read-region classifier for ``scbamqc bam --extra-metrics``.

Only the feature columns are used: ``seqname``, ``feature``, ``start``,
``end``, ``strand`` and the ``gene_id`` attribute (only when a GTF has no
``gene`` or ``transcript`` rows). Gene type attributes (``gene_biotype`` in
Ensembl, ``gene_ebi_biotype`` in VEuPathDB) are not read, so every annotated
gene counts.

Coordinates are converted from GTF (1-based, closed) to 0-based half-open
intervals, the same as pysam alignment blocks. For each (contig, strand) the
exons and the gene spans are merged into sorted, non-overlapping intervals.
Gene spans also include the exons, so an exon is always inside a gene span.
Features with strand ``.`` are put on both strands.

Classification of one read (aligned blocks, from ``AlignedSegment.get_blocks``,
so spliced ``N`` gaps and soft clips do not count; ``sense`` is the strand of
the transcript that the read comes from):

1. ``INTERGENIC`` or ``ANTISENSE`` if no aligned base is in a gene span on the
   sense strand: ``ANTISENSE`` if an aligned base is in a gene span on the
   opposite strand, ``INTERGENIC`` otherwise.
2. ``EXONIC`` if at least 50 % of the aligned bases are in sense exons.
3. ``INTRONIC`` if an aligned base is in a sense gene span but not in a sense
   exon.
4. ``INTERGENIC`` otherwise (less than 50 % exonic, rest outside the gene).

Rules 2 to 4 follow the Cell Ranger ``RE`` definition (exonic if >= 50 % of
the read is in an exon, intronic if not exonic and it intersects an intron).

For an unstranded library, both strands are sense: the read gets the first of
exonic, intronic, intergenic that one of the two strands gives.
"""

from __future__ import annotations

import gzip
import re
from bisect import bisect_right
from collections import defaultdict
from dataclasses import dataclass, field
from os import PathLike

import numpy as np

__all__ = [
    "ANTISENSE",
    "EXONIC",
    "INTERGENIC",
    "INTRONIC",
    "ContigAnnot",
    "load_gtf",
]

EXONIC, INTRONIC, INTERGENIC, ANTISENSE = 0, 1, 2, 3

Intervals = tuple[list[int], list[int]]  # (starts, ends), sorted, non-overlapping

_GENE_ID_RE = re.compile(r'gene_id\s+"?([^";]+)"?')


def _merge(starts: list[int], ends: list[int]) -> Intervals:
    """Sort and merge half-open intervals (touching intervals are merged)."""
    if not starts:
        return [], []
    s = np.asarray(starts, dtype=np.int64)
    e = np.asarray(ends, dtype=np.int64)
    order = np.argsort(s, kind="stable")
    s, e = s[order], e[order]
    out_s: list[int] = [int(s[0])]
    out_e: list[int] = [int(e[0])]
    for a, b in zip(s[1:].tolist(), e[1:].tolist()):
        if a <= out_e[-1]:
            if b > out_e[-1]:
                out_e[-1] = b
        else:
            out_s.append(a)
            out_e.append(b)
    return out_s, out_e


def _overlap(iv: Intervals, blocks: list[tuple[int, int]]) -> int:
    """Number of bases of ``blocks`` inside the merged intervals ``iv``."""
    starts, ends = iv
    n = len(starts)
    if n == 0:
        return 0
    tot = 0
    for s, e in blocks:
        i = bisect_right(ends, s)  # first interval that ends after s
        while i < n and starts[i] < e:
            tot += min(e, ends[i]) - max(s, starts[i])
            i += 1
    return tot


@dataclass
class ContigAnnot:
    """Merged exons and gene spans of one contig, by strand (``+``, ``-``)."""

    exons: dict[str, Intervals] = field(default_factory=dict)
    genes: dict[str, Intervals] = field(default_factory=dict)

    def classify(self, blocks: list[tuple[int, int]], sense: str | None) -> int:
        """Return ``EXONIC``, ``INTRONIC``, ``INTERGENIC`` or ``ANTISENSE``.

        ``sense`` is ``"+"`` or ``"-"``: the strand of the transcript that the
        read comes from (after the library strand is applied). ``None``
        (unstranded library): both strands are sense. The result is the first
        of exonic, intronic, intergenic that one of the strands gives, and is
        never ``ANTISENSE``.
        """
        if sense is None:
            best = min(self.classify(blocks, "+"), self.classify(blocks, "-"))
            return INTERGENIC if best == ANTISENSE else best
        empty: Intervals = ([], [])
        g_ov = _overlap(self.genes.get(sense, empty), blocks)
        if g_ov == 0:
            anti = "-" if sense == "+" else "+"
            if _overlap(self.genes.get(anti, empty), blocks) > 0:
                return ANTISENSE
            return INTERGENIC
        total = sum(e - s for s, e in blocks)
        ex_ov = _overlap(self.exons.get(sense, empty), blocks)
        if ex_ov > 0 and 2 * ex_ov >= total:
            return EXONIC
        if g_ov > ex_ov:
            return INTRONIC
        return INTERGENIC


def _open_text(path: str | PathLike):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p)


def load_gtf(path: str | PathLike) -> dict[str, ContigAnnot]:
    """Load the exons and gene spans of a GTF (or GFF-like) file.

    Gene spans come from ``gene`` rows. Without ``gene`` rows, they come from
    ``transcript`` rows. Without both, each ``gene_id`` spans its exons.
    Files ending in ``.gz`` are read with gzip.

    Returns
    -------
    ``{contig: ContigAnnot}``. Contigs with no exon and no gene are absent.
    """
    # (contig, strand) -> ([starts], [ends])
    exons: dict[tuple[str, str], tuple[list[int], list[int]]] = defaultdict(
        lambda: ([], [])
    )
    genes: dict[tuple[str, str], tuple[list[int], list[int]]] = defaultdict(
        lambda: ([], [])
    )
    transcripts: dict[tuple[str, str], tuple[list[int], list[int]]] = defaultdict(
        lambda: ([], [])
    )
    # (contig, strand, gene_id) -> [min start, max end], from exons
    exon_genes: dict[tuple[str, str, str], list[int]] = {}

    with _open_text(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            feature = cols[2]
            if feature == "exon":
                target = exons
            elif feature == "gene":
                target = genes
            elif feature == "transcript" or feature == "mRNA":
                target = transcripts
            else:
                continue
            contig = cols[0]
            start = int(cols[3]) - 1
            end = int(cols[4])
            strands = ("+", "-") if cols[6] not in ("+", "-") else (cols[6],)
            for strand in strands:
                ss, ee = target[(contig, strand)]
                ss.append(start)
                ee.append(end)
                if feature == "exon" and len(cols) > 8:
                    m = _GENE_ID_RE.search(cols[8])
                    if m:
                        key = (contig, strand, m.group(1))
                        span = exon_genes.get(key)
                        if span is None:
                            exon_genes[key] = [start, end]
                        else:
                            span[0] = min(span[0], start)
                            span[1] = max(span[1], end)

    if not genes:
        if transcripts:
            genes = transcripts
        else:
            for (contig, strand, _gid), (s, e) in exon_genes.items():
                ss, ee = genes[(contig, strand)]
                ss.append(s)
                ee.append(e)

    out: dict[str, ContigAnnot] = {}
    for (contig, strand), (es, ee) in exons.items():
        out.setdefault(contig, ContigAnnot()).exons[strand] = _merge(es, ee)
    for (contig, strand), (gs, ge) in genes.items():
        # gene spans include the exons, so that exon overlap <= gene overlap
        es, ee = exons.get((contig, strand), ([], []))
        out.setdefault(contig, ContigAnnot()).genes[strand] = _merge(gs + es, ge + ee)
    for (contig, strand), (es, ee) in exons.items():
        ca = out[contig]
        if strand not in ca.genes:
            ca.genes[strand] = _merge(es, ee)
    return out
