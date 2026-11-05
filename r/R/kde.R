# Gaussian kernel density, equal to scipy.stats.gaussian_kde with its default
# bandwidth (Scott's rule, full covariance), evaluated at the fitted points.
#
# X: numeric matrix, one row per point. Returns the density at each row.
# Points are processed in blocks of `chunk` rows to limit memory use.
kde_at_points <- function(X, chunk = 1000L) {
  n <- nrow(X)
  d <- ncol(X)
  # Kernel covariance: sample covariance (n - 1) times Scott's factor squared.
  R <- chol(stats::cov(X) * n^(-2 / (d + 4)))
  # Whiten: Euclidean distances of Y are Mahalanobis distances under the kernel.
  Y <- sweep(X, 2, colMeans(X)) %*% backsolve(R, diag(d))
  dens <- numeric(n)
  for (start in seq(1L, n, by = chunk)) {
    i <- start:min(start + chunk - 1L, n)
    # Sum the squared differences directly. The shorter |a|^2 + |b|^2 - 2ab
    # leaves rounding residue, which changes the scores of far outliers.
    d2 <- 0
    for (k in seq_len(d)) d2 <- d2 + outer(Y[i, k], Y[, k], "-")^2
    dens[i] <- rowSums(exp(-0.5 * d2))
  }
  dens / (n * (2 * pi)^(d / 2) * prod(diag(R)))
}
