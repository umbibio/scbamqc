skip_if_not_installed("SeuratObject")
skip_if_not_installed("Matrix")

make_seurat <- function(n_cells = 300, n_genes = 50, seed = 0) {
  set.seed(seed)
  X <- matrix(rpois(n_cells * n_genes, 0.3), n_genes, n_cells)
  X[, 1:20] <- 0 # empty barcodes
  X[1:3, ] <- 0 # undetected genes
  X[4, ] <- 0
  X[4, 21:25] <- 1 # gene detected in exactly 5 cells
  dimnames(X) <- list(paste0("g", seq_len(n_genes)), sprintf("C%04d-1", seq_len(n_cells)))
  suppressWarnings(SeuratObject::CreateSeuratObject(Matrix::Matrix(X, sparse = TRUE)))
}

test_that("basic_qc filters barcodes, then genes, and adds metrics", {
  obj <- make_seurat()
  X <- as.matrix(SeuratObject::LayerData(obj, layer = "counts"))
  cell_ok <- colSums(X > 0) >= 5
  gene_ok <- rowSums(X[, cell_ok] > 0) >= 5
  kept <- X[gene_ok, cell_ok]

  out <- basic_qc(obj)
  expect_equal(dim(out), c(sum(gene_ok), sum(cell_ok)))
  expect_equal(unname(out$total_counts), unname(colSums(kept)))
  expect_equal(unname(out$n_genes_by_counts), unname(colSums(kept > 0)))
  expect_gte(min(out$barcode_rank), 1)
  expect_equal(SeuratObject::Misc(out, "scbamqc")$basic_qc,
               list(min_genes = 5, min_cells = 5))
  expect_equal(dim(obj), c(50L, 300L)) # input unchanged
})

test_that("the full workflow runs on a Seurat object", {
  obj <- make_seurat(n_cells = 1500)
  tab <- read.csv(test_path("fixtures", "synthetic_python.csv.gz"), row.names = 1)
  tab <- tab[, c("n_reads", "mean_nm", "mean_clip_frac", "std_clip_frac", "mean_base_qual")]
  rownames(tab) <- colnames(obj)
  tab <- tab[-(1:10), ] # barcodes without BAM-QC rows are removed

  obj <- add_bam_qc(obj, tab)
  expect_equal(ncol(obj), 1490)
  obj <- basic_qc(obj)
  obj <- qc_score(obj, min_reads = 200)
  obj <- select_cells(obj, min_counts = 10)

  meta <- obj[[]]
  expect_true(all(c("qc_score", "gmm_label", "selected") %in% colnames(meta)))
  expect_type(meta$selected, "logical")
  expect_named(SeuratObject::Misc(obj, "scbamqc"),
               c("basic_qc", "qc_score", "select_cells"))
})
