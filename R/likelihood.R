log_lik <- function(y, X, Z, Beta, sigma, lambda, lam_idx, diffs){
  stopifnot(is.atomic(y), is.null(dim(y)))
  n <- length(y)
  stopifnot(is.matrix(X), nrow(X) == n)
  p <- ncol(X)
  stopifnot(is.matrix(Z), nrow(Z) == n)
  q <- ncol(Z)
  if(q >= n) stop("n <= q not yet supported")
  stopifnot(is.atomic(Beta), is.null(dim(Beta)), length(Beta) == p)
  stopifnot(is.atomic(sigma), length(sigma) == 1)
  if(abs(sigma) < 1e-12){
    warning("sigma is approximately zero, computations may be unstable")
  }
  stopifnot(is.atomic(lambda), is.null(dim(lambda)))
  d <- length(lambda)
  stopifnot(is.atomic(lam_idx), length(lam_idx) == q, all(sort(unique(lam_idx))
  == 1:d))
  stopifnot(diffs %in% 0:2)

  log_lik_rcpp(y, X, Z, Beta, sigma, lambda, lam_idx, diffs)
}

score <- function(y, X, Z, Beta, sigma, lambda, lam_idx)
{
  log_lik(y = y, X = X, Z = Z, Beta = Beta, sigma = sigma, lambda = lambda,
    lam_idx = lam_idx, diffs = 1)[[2]] * c(rep(1, length(Beta)), sigma, lambda)
}

fish_inf <- function(y, X, Z, Beta, sigma, lambda, lam_idx)
{
  out <- log_lik(y = y, X = X, Z = Z, Beta = Beta, sigma = sigma, lambda =
    lambda, lam_idx = lam_idx, diffs = 2)
  scale_vec <- c(rep(1, length(Beta)), sigma, lambda)
  scale_vec * out[[3]] * scale_vec
}
