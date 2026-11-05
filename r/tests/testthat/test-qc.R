py <- function() {
  read.csv(test_path("fixtures", "synthetic_python.csv.gz"), row.names = 1)
}

test_that("qc_score equals the Python package (scipy gaussian_kde)", {
  df <- py()
  paper <- qc_score(df, min_reads = 200, base_qual_scale = PAPER_BASE_QUAL_SCALE)
  expect_equal(is.na(paper$qc_score), is.na(df$py_q_paper))
  expect_equal(paper$qc_score, df$py_q_paper, tolerance = 1e-10)

  med <- qc_score(df, min_reads = 200)
  expect_equal(med$qc_score, df$py_q_median, tolerance = 1e-10)
})

test_that("qc_score gate is strict and scores are from 0 to 1", {
  df <- py()
  q <- qc_score(df, min_reads = 200)$qc_score
  gated <- df$n_reads > 200 & df$mean_nm < 1
  expect_true(all(!is.na(q[gated])))
  expect_true(all(is.na(q[!gated])))
  expect_equal(range(q, na.rm = TRUE), c(0, 1))

  df$n_reads[1] <- 200
  expect_true(is.na(qc_score(df, min_reads = 200)$qc_score[1]))
})

test_that("base_qual_scale changes values, not the order", {
  df <- py()
  q10 <- qc_score(df, min_reads = 200, base_qual_scale = 10)$qc_score
  q1 <- qc_score(df, min_reads = 200, base_qual_scale = 1)$qc_score
  ok <- !is.na(q10)
  expect_equal(stats::cor(q10[ok], q1[ok], method = "spearman"), 1)
  expect_false(isTRUE(all.equal(q10[ok], q1[ok])))
})

test_that("qc_score records the median divisor and checks its input", {
  df <- py()
  out <- qc_score(df, min_reads = 200)
  gated <- df$n_reads > 200 & df$mean_nm < 1
  rec <- attr(out, "scbamqc")$qc_score
  expect_equal(rec$base_qual_scale, median(df$mean_base_qual[gated]))
  expect_equal(rec$n_gated, sum(gated))
  expect_error(qc_score(df, base_qual_scale = "mean"), "median")
  expect_error(qc_score(df, min_reads = 1e9), "too few")
  expect_error(qc_score(df[, -1]), "Missing columns: n_reads")
})

test_that("select_cells agrees with a converged scikit-learn fit", {
  df <- py()
  df$qc_score <- df$py_q_median
  out <- select_cells(df, min_counts = 100)
  expect_type(out$selected, "logical")
  expect_gt(mean(out$selected == df$py_selected_converged), 0.99)

  lab <- out$gmm_label
  expect_setequal(unique(lab), c(-1L, 0L, 1L))
  expect_true(all(lab[is.na(df$qc_score) | df$qc_score <= 0] == -1L))
  lc <- log10(pmax(df$total_counts, 1))
  expect_gt(mean(lc[lab == 0]), mean(lc[lab == 1]))
  expect_equal(out$selected,
               lab == 0 & df$total_counts > 100 & !is.na(df$qc_score) & df$qc_score > 0)
})

test_that("select_cells does not depend on the seed or change the global RNG", {
  df <- py()
  df$qc_score <- df$py_q_median
  set.seed(42)
  before <- .Random.seed
  a <- select_cells(df, seed = 0)$selected
  expect_identical(.Random.seed, before)
  b <- select_cells(df, seed = 7)$selected
  expect_identical(a, b)
})

test_that("add_bam_qc keeps the order of the object", {
  df <- py()
  obj <- data.frame(row.names = sprintf("C%04d-1", 1:300), x = 1:300)
  tab <- df[1:100, 1:3]
  rownames(tab) <- rev(rownames(obj)[seq(1, 300, by = 3)])
  out <- add_bam_qc(obj, tab)
  expect_equal(rownames(out), rownames(obj)[rownames(obj) %in% rownames(tab)])
  expect_equal(out$n_reads, tab[rownames(out), "n_reads"])
  expect_error(add_bam_qc(obj, data.frame(n_reads = 1, row.names = "AAAA")),
               "barcode format")
  expect_error(add_bam_qc(out, tab), "already in the metadata")
})

test_that("read_bam_qc reads the table of scbamqc bam", {
  path <- tempfile(fileext = ".csv.gz")
  on.exit(unlink(path))
  con <- gzfile(path, "w")
  writeLines(c("cb,n_reads,mean_nm", "AAAC-1,10,0.5", "AAAG-1,20,0.1"), con)
  close(con)
  tab <- read_bam_qc(path)
  expect_equal(rownames(tab), c("AAAC-1", "AAAG-1"))
  expect_equal(tab$n_reads, c(10L, 20L))
})
