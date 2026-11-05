#' Divisor of `mean_base_qual` in the TgRN paper (Phred / 10)
#'
#' Use it as `base_qual_scale` in [qc_score()] to repeat the paper's scores.
#' @export
PAPER_BASE_QUAL_SCALE <- 10

#' BAM-QC columns that enter the kernel density
#' @export
DEFAULT_FEATURES <- c("mean_clip_frac", "std_clip_frac", "mean_base_qual")

# ---------------------------------------------------------------------------
# Helpers for Seurat objects and data frames
# ---------------------------------------------------------------------------

is_seurat <- function(x) inherits(x, "Seurat")

check_input <- function(x) {
  if (is_seurat(x)) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("Install the 'SeuratObject' package to use Seurat objects.", call. = FALSE)
    }
  } else if (!is.data.frame(x)) {
    stop("`object` must be a Seurat object or a data frame.", call. = FALSE)
  }
  invisible(x)
}

get_meta <- function(x) if (is_seurat(x)) x[[]] else x

set_meta <- function(x, name, value) {
  x[[name]] <- value
  x
}

need_cols <- function(meta, cols) {
  miss <- setdiff(cols, colnames(meta))
  if (length(miss)) {
    stop("Missing columns: ", paste(miss, collapse = ", "),
         ". Run add_bam_qc() and basic_qc() first.", call. = FALSE)
  }
}

record <- function(x, step, params) {
  if (is_seurat(x)) {
    rec <- SeuratObject::Misc(x, slot = "scbamqc")
    if (is.null(rec)) rec <- list()
    rec[[step]] <- params
    x <- suppressWarnings(SeuratObject::`Misc<-`(x, slot = "scbamqc", value = rec))
  } else {
    rec <- attr(x, "scbamqc")
    if (is.null(rec)) rec <- list()
    rec[[step]] <- params
    attr(x, "scbamqc") <- rec
  }
  x
}

with_seed <- function(seed, expr) {
  has <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (has) old <- get(".Random.seed", envir = globalenv())
  on.exit(
    if (has) assign(".Random.seed", old, envir = globalenv())
    else rm(".Random.seed", envir = globalenv())
  )
  set.seed(seed)
  expr
}

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

#' Read a BAM-QC table
#'
#' Reads the per-barcode table that `scbamqc bam` writes
#' (`<prefix>_per_cell.csv.gz`).
#'
#' @param path Path to the table.
#' @return A data frame with one row for each barcode. The row names are the
#'   barcodes.
#' @export
read_bam_qc <- function(path) {
  utils::read.csv(path, row.names = 1, check.names = FALSE)
}

#' Add the BAM-QC table to the cell metadata
#'
#' Keeps only the barcodes that are in the BAM-QC table, in their order in
#' `object`. Adds the table columns to the metadata.
#'
#' @param object A Seurat object, or a data frame with barcodes as row names.
#' @param bam_qc A data frame from [read_bam_qc()], or a path to the table.
#' @return `object` with the matching barcodes and the BAM-QC columns.
#' @export
add_bam_qc <- function(object, bam_qc) {
  check_input(object)
  if (!is.data.frame(bam_qc)) bam_qc <- read_bam_qc(bam_qc)

  cells <- if (is_seurat(object)) colnames(object) else rownames(object)
  keep <- cells[cells %in% rownames(bam_qc)]
  if (!length(keep)) {
    stop("No barcode of `object` is in the BAM-QC table. ",
         "Make sure that both use the same barcode format (for example the '-1' suffix).",
         call. = FALSE)
  }
  overlap <- intersect(colnames(get_meta(object)), colnames(bam_qc))
  if (length(overlap)) {
    stop("Columns already in the metadata: ", paste(overlap, collapse = ", "),
         call. = FALSE)
  }

  add <- bam_qc[keep, , drop = FALSE]
  if (is_seurat(object)) {
    object <- subset(object, cells = keep)
    SeuratObject::AddMetaData(object, add)
  } else {
    cbind(object[keep, , drop = FALSE], add)
  }
}

#' Remove empty barcodes and rare genes, and add count metrics
#'
#' Removes barcodes with fewer than `min_genes` detected genes. Then removes
#' genes detected in fewer than `min_cells` of the remaining barcodes. The
#' metrics are calculated after both steps.
#'
#' Adds to the cell metadata: `n_genes_by_counts`, `total_counts` and
#' `barcode_rank` (rank by decreasing `total_counts`).
#'
#' @param object A Seurat object.
#' @param min_genes Minimum number of detected genes for a barcode.
#' @param min_cells Minimum number of barcodes for a gene.
#' @param assay Assay to use. The default is the active assay.
#' @param layer Layer with the raw counts.
#' @return The filtered Seurat object.
#' @export
basic_qc <- function(object, min_genes = 5, min_cells = 5, assay = NULL,
                     layer = "counts") {
  if (!is_seurat(object)) stop("`object` must be a Seurat object.", call. = FALSE)
  check_input(object)
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)

  counts <- SeuratObject::LayerData(object, assay = assay, layer = layer)
  n_genes <- Matrix::colSums(counts > 0)
  cells <- colnames(counts)[n_genes >= min_genes]
  n_cells <- Matrix::rowSums(counts[, cells, drop = FALSE] > 0)
  genes <- rownames(counts)[n_cells >= min_cells]

  object <- subset(object, cells = cells, features = genes)
  counts <- SeuratObject::LayerData(object, assay = assay, layer = layer)
  total <- Matrix::colSums(counts)
  object <- set_meta(object, "n_genes_by_counts", Matrix::colSums(counts > 0))
  object <- set_meta(object, "total_counts", total)
  object <- set_meta(object, "barcode_rank", rank(-total, ties.method = "average"))

  record(object, "basic_qc", list(min_genes = min_genes, min_cells = min_cells))
}

# ---------------------------------------------------------------------------
# Quality score
# ---------------------------------------------------------------------------

#' Quality score from the alignment profile
#'
#' A barcode gets a score only if `n_reads > min_reads` and
#' `mean_nm < max_mean_nm`. A Gaussian kernel density (Scott's rule, as
#' `scipy.stats.gaussian_kde`) is fitted to the `features` of these barcodes.
#' The score is `log10(density + 1)`, scaled from 0 to 1. Other barcodes get
#' `NA`.
#'
#' @param object A Seurat object or a data frame with the BAM-QC columns.
#' @param min_reads Minimum number of reads. Set it for each library.
#' @param max_mean_nm Maximum mean number of mismatches.
#' @param base_qual_scale Divisor of `mean_base_qual`. `"median"` uses the
#'   median of the scored barcodes. Use [PAPER_BASE_QUAL_SCALE] to repeat the
#'   paper. The divisor does not change the order of the barcodes.
#' @param features Columns that enter the density.
#' @param key_added Name of the new metadata column.
#' @return `object` with the score in `key_added`.
#' @export
qc_score <- function(object, min_reads = 100, max_mean_nm = 1,
                     base_qual_scale = "median", features = DEFAULT_FEATURES,
                     key_added = "qc_score") {
  check_input(object)
  meta <- get_meta(object)
  need_cols(meta, c("n_reads", "mean_nm", features))

  mask <- meta$n_reads > min_reads & meta$mean_nm < max_mean_nm
  mask[is.na(mask)] <- FALSE
  n <- sum(mask)
  if (n <= length(features)) {
    stop(sprintf(paste0("Only %d barcodes pass the read gate (n_reads > %s, ",
                        "mean_nm < %s): too few to estimate a density."),
                 n, min_reads, max_mean_nm), call. = FALSE)
  }

  if (identical(base_qual_scale, "median")) {
    base_qual_scale <- stats::median(meta$mean_base_qual[mask])
  } else if (!is.numeric(base_qual_scale) || length(base_qual_scale) != 1 ||
             !(base_qual_scale > 0)) {
    stop("`base_qual_scale` must be a positive number or \"median\".", call. = FALSE)
  }

  X <- vapply(features, function(f) {
    v <- as.numeric(meta[[f]][mask])
    if (f == "mean_base_qual") v / base_qual_scale else v
  }, numeric(n))
  X <- matrix(X, nrow = n)

  v <- log10(kde_at_points(X) + 1)
  v <- v - min(v)
  score <- rep(NA_real_, nrow(meta))
  score[mask] <- v / max(v)

  object <- set_meta(object, key_added, score)
  record(object, "qc_score", list(
    min_reads = min_reads, max_mean_nm = max_mean_nm,
    base_qual_scale = base_qual_scale, features = features, n_gated = n
  ))
}

# ---------------------------------------------------------------------------
# Mixture model and selection
# ---------------------------------------------------------------------------

#' Gaussian mixture labels on counts and score
#'
#' Fits a Gaussian mixture ([mclust::Mclust()], full covariance) to
#' `log10(total_counts)` and `qc_score` of the barcodes with `qc_score > 0`.
#' The component with the highest mean count gets label 0. Barcodes not in
#' the fit get label -1.
#'
#' EM runs until the log-likelihood changes by less than `tol`. With a small
#' `tol`, the fit does not depend on `seed`.
#'
#' @inheritParams qc_score
#' @param n_components Number of mixture components.
#' @param seed Random seed for the initialization. The global random state
#'   does not change.
#' @param tol Convergence tolerance of EM (relative change of log-likelihood).
#' @param counts_key Metadata column with the total counts.
#' @param score_key Metadata column with the quality score.
#' @return `object` with the labels in `key_added`.
#' @export
gmm_labels <- function(object, n_components = 2, seed = 0, tol = 1e-10,
                       counts_key = "total_counts", score_key = "qc_score",
                       key_added = "gmm_label") {
  check_input(object)
  meta <- get_meta(object)
  need_cols(meta, c(counts_key, score_key))

  log_counts <- log10(pmax(meta[[counts_key]], 1))
  score <- meta[[score_key]]
  mask <- !is.na(score) & score > 0
  X <- cbind(log_counts[mask], score[mask])

  fit <- with_seed(seed, mclust::Mclust(
    X, G = n_components, modelNames = "VVV", verbose = FALSE,
    control = mclust::emControl(tol = c(tol, sqrt(.Machine$double.eps)))
  ))
  if (is.null(fit)) stop("The Gaussian mixture fit failed.", call. = FALSE)

  cl <- fit$classification
  mean_counts <- tapply(log_counts[mask], factor(cl, levels = seq_len(n_components)), mean)
  rk <- rank(-mean_counts, ties.method = "first") - 1L
  labels <- rep(-1L, nrow(meta))
  labels[mask] <- as.integer(rk[cl])

  set_meta(object, key_added, labels)
}

#' Select candidate cells
#'
#' A barcode is selected if it is in the mixture component with the highest
#' counts (see [gmm_labels()]), `total_counts > min_counts` and
#' `qc_score > min_score`. Also adds the mixture labels as `gmm_label`.
#'
#' @inheritParams gmm_labels
#' @param min_counts Minimum total counts.
#' @param min_score Minimum quality score.
#' @return `object` with a logical column `key_added`.
#' @export
select_cells <- function(object, min_counts = 100, min_score = 0,
                         n_components = 2, seed = 0, tol = 1e-10,
                         counts_key = "total_counts", score_key = "qc_score",
                         key_added = "selected") {
  object <- gmm_labels(object, n_components = n_components, seed = seed,
                       tol = tol, counts_key = counts_key, score_key = score_key)
  meta <- get_meta(object)
  score <- meta[[score_key]]
  selected <- meta[[counts_key]] > min_counts & !is.na(score) &
    score > min_score & meta$gmm_label == 0
  object <- set_meta(object, key_added, selected)
  record(object, "select_cells", list(
    min_counts = min_counts, min_score = min_score, n_components = n_components,
    seed = seed, tol = tol, n_selected = sum(selected)
  ))
}
