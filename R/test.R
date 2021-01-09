score_test <- function(y, X, Z, Beta, sigma, lambda, lam_idx, test_idx){
  components <- log_lik(y = y, X = X, Z = Z, Beta = Beta, sigma = sigma,
                        lambda = lambda, lam_idx = lam_idx, diffs = 2)
  theta <- c(Beta, sigma, lambda)
  m <- length(theta)
  stopifnot(is.atomic(test_idx), is.null(dim(test_idx)), all(test_idx %in% 1:m),
            length(test_idx) == length(unique(test_idx)))
  s <- components[[2]][test_idx]
  I <- components[[3]][test_idx, test_idx]
  test_stat <- sum(s * qr.solve(I, s))
  df <- length(test_idx)
  p_val <- pchisq(q = test_stat, df = df, lower.tail = FALSE)
  c("chi_sq" = test_stat, "df" = df, "p_val" = p_val)
}