loglik <- function(y, X, Z, Beta, sigma, lambda, id, diffs)
{
  stopifnot(is.atomic(y), is.null(dim(y)))
  n <- length(y)
  stopifnot(is.matrix(X), nrow(X) == n)
  p <- ncol(X)
  stopifnot(is.matrix(Z), nrow(Z) == n)
  q <- ncol(Z)
  stopifnot(is.atomic(Beta), length(Beta) == p)
  stopifnot(is.atomic(sigma), length(sigma) == 1, sigma > 0)
  stopifnot(is.atomic(lambda))
  d <- length(lambda)
  stopifnot(is.atomic(id), length(id) == q, sort(unique(id)) == 1:d)
  stopifnot(diffs %in% 0:2)
  
  loglik_rcpp(y, X, Z, Beta, sigma, lambda, id, diffs)
}