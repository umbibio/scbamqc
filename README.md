# scbamqc

Alignment QC for single-cell BAM files, and cell selection from that QC.

- `scbamqc bam` writes one row of alignment statistics for each cell barcode.
- `scbamqc.qc` calculates a quality score for each barcode and selects cells.
- `scbamqc.pl` plots the statistics and the score.

## Install

```bash
pip install scbamqc
```

Add extras for more functions:

| Extra | Adds | Use it for |
|---|---|---|
| `bam` | pysam | `scbamqc bam` and `scbamqc tags` (Linux, macOS) |
| `anndata` | anndata, h5py | `scbamqc.io.read_10x_h5` |
| `plot` | matplotlib | `scbamqc.pl` and `scbamqc bam --plots` |
| `all` | all of the above | |

For example: `pip install "scbamqc[bam,plot]"`.

## Command line

The BAM file must be sorted by coordinate and indexed (`samtools index`).

```bash
# Cell Ranger
scbamqc bam possorted_genome_bam.bam -j 16 --out-prefix S1

# STARsolo SmartSeq: one read group for each cell
scbamqc bam Aligned.sortedByCoord.out.bam --barcode-tag RG
```

The tool writes `S1_per_cell.csv.gz`, with one row for each barcode. To also
write histograms to `S1_plots.png`, add `--plots`.

For all options, run `scbamqc bam --help`.

### Other pipelines

Different pipelines write the barcode and the mismatch count in different BAM
tags. Find the tags before you run `bam`:

```bash
scbamqc tags sample.bam
scbamqc bam sample.bam --barcode-tag CB --nm-tag NM
```

`bam` shows a warning if most reads do not have the tags you set.

> [!WARNING]
> Without a mismatch tag, `mean_nm` is empty (`NaN`). Then no barcode gets a
> quality score.

## Python

```python
import scbamqc as sbq
from scbamqc.io import read_10x_h5

adata = read_10x_h5("raw_feature_bc_matrix.h5")
adata = sbq.add_bam_qc(adata, "S1_per_cell.csv.gz")
adata = sbq.basic_qc(adata)
sbq.qc_score(adata, min_reads=150)
sbq.select_cells(adata)

sbq.pl.qc_score_scatter(adata, min_score=0, min_counts=100)
sbq.pl.qc_score_hist(adata, threshold=0.9)

cells = adata[adata.obs["selected"] & (adata.obs["qc_score"] > 0.9)].copy()
```

| Step | Result |
|---|---|
| `add_bam_qc` | Keeps the barcodes in the BAM table. Adds the table columns to `obs`. |
| `basic_qc` | Removes barcodes with fewer than 5 genes and genes in fewer than 5 barcodes. Adds `total_counts`. |
| `qc_score` | Adds `obs["qc_score"]`, from 0 to 1. |
| `select_cells` | Adds `obs["gmm_label"]` and `obs["selected"]`. |

### The quality score

1. A barcode gets a score only if `n_reads > min_reads` and `mean_nm < 1`.
   Other barcodes get `NaN`.
2. A kernel density is fitted to three values of these barcodes:
   `mean_clip_frac`, `std_clip_frac` and `mean_base_qual` divided by its
   median.
3. The score is `log10(density + 1)`, scaled from 0 to 1.

A high score shows that the barcode aligns like most other barcodes.

Set `min_reads` for each library. The paper used values from 50 to 3,000.

`base_qual_scale` sets the divisor of `mean_base_qual`. The default is the
median of the scored barcodes. The divisor does not change the order of the
barcodes. It can change the score values only when the densities are near 1.

To use the same settings as the TgRN paper:

```python
sbq.qc_score(adata, min_reads=150, base_qual_scale=sbq.PAPER_BASE_QUAL_SCALE)
sbq.select_cells(adata, tol=sbq.PAPER_GMM_TOL)
```

### Cell selection

`select_cells` fits a Gaussian mixture with two components to
`log10(total_counts)` and `qc_score`. A barcode is selected if:

- it is in the component with the higher counts,
- `total_counts > min_counts` (default 100),
- `qc_score > min_score` (default 0).

The mixture fit runs until it converges. The paper stopped the fit early
(`tol=PAPER_GMM_TOL`). On the paper sample, the converged fit selects 15,536
cells and the paper fit selects 16,120. All cells of the converged fit are in
the paper set.

## R

The R package in [`r/`](r/README.md) does the same steps for Seurat objects.
It reads the table that `scbamqc bam` writes.

```r
remotes::install_github("umbibio/scbamqc", subdir = "r")
```

## Development

```bash
uv sync
uv run pytest
```

The tests in `tests/test_paper.py` compare the results with the paper. They
need the data in `../../data`, or in the folder set by `SCBAMQC_TEST_DATA`.
Without the data, pytest skips them.

## License

BSD 3-Clause.
