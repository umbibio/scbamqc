# scbamqc for R

Quality score and cell selection for Seurat objects and data frames.

The R package does not read BAM files. Use the `scbamqc bam` command of the
Python package to make the per-barcode table (see the
[main README](../README.md)).

## Install

```r
# install.packages("remotes")
remotes::install_github("umbibio/scbamqc", subdir = "r")
```

The package needs only `mclust`. Install `SeuratObject` to use Seurat objects
and `ggplot2` to make plots.

## Use

```r
library(scbamqc)

obj <- add_bam_qc(obj, "S1_per_cell.csv.gz")  # Seurat object with raw counts
obj <- basic_qc(obj)
obj <- qc_score(obj, min_reads = 150)
obj <- select_cells(obj)

qc_score_scatter(obj, min_score = 0, min_counts = 100)
qc_score_hist(obj, threshold = 0.9)

cells <- subset(obj, cells = colnames(obj)[obj$selected & obj$qc_score > 0.9])
```

Each function returns the object with new metadata columns. `qc_score()`,
`select_cells()` and the plots also accept a data frame.

## Differences from the Python package

| Step | Result |
|---|---|
| `qc_score()` | Same scores as Python (difference less than 1e-12). |
| `select_cells()` | Almost the same cells as Python. Both fits run until they converge. |

On the paper sample (RH_221220.WT), R selects 15,533 cells and Python 15,536.
The R package cannot repeat the paper fit, which stopped early and selected
16,120 cells. To repeat the paper, use the Python package with
`tol=PAPER_GMM_TOL`.
