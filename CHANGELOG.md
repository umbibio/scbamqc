# Changelog

All notable changes to this project are in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- `scbamqc bam --extra-metrics`: adds columns for read regions, UMIs, polyA and TSO reads, unique MAPQ and mismatch counts. Without the option, the output does not change. With the option, the output file is `PREFIX_per_cell_extra.csv.gz`.
- Read regions come from the Cell Ranger tags `RE`, `AN`, `GX` and `TX`. Without these tags, they come from the GTF file of `--gtf`.
- Options for `--extra-metrics`: `--gtf`, `--strand` (`forward`, `reverse` or `unstranded`), `--umi-tag` and `--mapq-unique`.
- `--region-source` option (`auto`, `tags` or `gtf`): selects the source of the read regions. Use `gtf` to compare the GTF regions with the Cell Ranger `RE` tags.
- `n_umi_counted` column: the reads that Cell Ranger counts as UMIs (`xf` tag, bit 8).
- `scbamqc.annot.load_gtf()`: reads the exons and genes of a GTF file. It does not need gene type attributes.
- `collect_metrics()` options: `extra_metrics`, `gtf`, `strand`, `umi_tag`, `mapq_unique`, `region_source`.
- `count_tags()` option `contigs`: sample only the reads of these contigs.

### Removed
- Support for Python 3.11. scbamqc needs Python 3.12 or later.

## [0.3.0] - 2026-11-17
### Changed
- The version numbers continue from the `bam_qc_summary.py` script. The paper uses version 0.2.1 of that script. The package goes from 0.1.0 to 0.3.0, so that no package version has the same number as the script.
- The BAM metrics are the same as in `bam_qc_summary.py` v0.2.1. For the cell selection of the paper, use `base_qual_scale=PAPER_BASE_QUAL_SCALE` and `tol=PAPER_GMM_TOL`.
- The R package has the same version number as the Python package. `bump-my-version bump release` changes the version in `r/DESCRIPTION`.

## [0.1.0] - 2026-11-06
### Added
- `scbamqc bam`: the `bam_qc_summary.py` v0.2.1 script as a command. The options and output files are the same.
- `--nm-tag` option for aligners that write `NM` and not `nM`. The default is `nM`.
- `scbamqc tags`: lists the BAM tags on the first mapped reads.
- `scbamqc bam` shows a warning if most reads do not have the barcode tag or the mismatch tag. Use `--no-check-tags` to skip the check.
- `scbamqc.qc`: `add_bam_qc()`, `basic_qc()`, `qc_score()`, `gmm_labels()`, `select_cells()`. With `base_qual_scale=PAPER_BASE_QUAL_SCALE` and `tol=PAPER_GMM_TOL`, they give the same cells as the TgRN paper.
- `qc_score()` options: `base_qual_scale`, `max_mean_nm`, `features`. The default `base_qual_scale` is the median `mean_base_qual` of the scored barcodes. For the paper value, use `PAPER_BASE_QUAL_SCALE` (10).
- `select_cells()` and `gmm_labels()` options: `tol`, `max_iter`. The default mixture fit runs until it converges. For the paper fit, which stops early, use `tol=PAPER_GMM_TOL`.
- Optional extras: `bam` (pysam), `anndata` (anndata, h5py), `plot` (matplotlib), `all`.
- `scbamqc.pl`: `bam_qc_hist()`, `qc_score_scatter()`, `qc_score_hist()`.
- `scbamqc.io.read_10x_h5()`. Scanpy is not necessary.
- Support for Python 3.11 to 3.14.
- R package in `r/` for Seurat objects and data frames: `add_bam_qc()`, `basic_qc()`, `qc_score()`, `select_cells()` and ggplot2 plots.

### Changed
- `scbamqc bam` writes plots only with `--plots`. matplotlib is optional.

### Fixed
- Version 0.2.1 wrote the plots only with `--no-plots`. Now `--plots` writes them.
- `scbamqc bam` does not fail when no read has the barcode tag.
