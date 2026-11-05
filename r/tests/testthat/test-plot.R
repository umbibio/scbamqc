skip_if_not_installed("ggplot2")

scored <- function() {
  df <- read.csv(test_path("fixtures", "synthetic_python.csv.gz"), row.names = 1)
  df$qc_score <- df$py_q_median
  df$selected <- df$py_selected_converged
  df
}

test_that("qc_score_scatter builds one panel for each metric", {
  p <- qc_score_scatter(scored(), x = c("total_counts", "mean_nm"),
                        min_score = 0, min_counts = 100)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_equal(nrow(b$layout$layout), 2)
  expect_match(p$labels$title, "selected \\(blue\\)")
})

test_that("qc_score_hist counts each group", {
  df <- scored()
  p <- qc_score_hist(df, threshold = 0.9)
  expect_s3_class(p, "ggplot")
  ok <- !is.na(df$qc_score)
  n_ret <- sum(ok & df$selected & df$qc_score > 0.9)
  expect_equal(p$labels$title,
               sprintf("%s retained of %s scored", format(n_ret, big.mark = ","),
                       format(sum(ok), big.mark = ",")))
  b <- ggplot2::ggplot_build(p)
  expect_equal(sum(b$data[[1]]$count), sum(ok))
})
