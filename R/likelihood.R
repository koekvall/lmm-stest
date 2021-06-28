#' Compute log-likelihood, scaled score, and scaled information
#'
#' @description{
#' Computed quantities assume a linear mixed model with independent random
#' random effects, parameterized in terms of regression coefficients and
#' standard deviations of the random effects and the error term (see details).
#' 
#' The order of parameters in the score function is as in c(Beta, sigma,
#' lambda). 
#' 
#' The elements of the score corresponding to scale parameters (sigma,
#' lambda) are scaled by (1 / scale parameter) (see details).
#' }
#' @param y A vector of observed responses.
#' @param X A matrix of predictors whose i:th row corresponds to the i:th
#'   element in y.
#' @param Z A design matrix for the random effects whose i:th row corresponds
#'   to the i:th element in y.
#' @param Beta A vector of regression coefficients of length ncol(X).
#' @param sigma The standard deviation of the error term.
#' @param lambda A vector of scale parameters (standard deviations) for the
#'   random effects.
#' @param lam_idx A vector of length ncol(Z) whose j:th element indicates
#'   which element of lambda scales the j:th random effect.
#' @param diffs An integer indicating whether to compute only the
#'   log-likelihood (0), the log-likelihood and scaled score (1), or
#'   the log-likelihood, scaled score, and scaled information (2).
#'
#' @return A list with log-likelihood ("loglik"), scaled score ("scaled_score")
#'   and scaled information ("scaled_inf").
#'
#' @details{
#'   The linear mixed model assumed is y = X \%*\% Beta + Z \%*\% u + e,
#'   where the vector of random effects u is from a multivariate normal
#'   distribution with mean zero and diagonal covariance matrix. The j:th
#'   diagonal element of that covariance matrix -- the variance of the j:th
#'   random effect -- is equal to lambda[lam_idx[j]]^2.
#'  
#'   The elements of e are independent draws from a normal distribution
#'   with mean 0 and variance sigma^2.
#'   
#'   The scaling of the score is done for scale parameters only. For example,
#'   the scaled score for sigma is the usual score times (1 / sigma).
#'   The regular score can be obtained by loglik(..)[[2]] * c(rep(1,
#'   length(Beta)), sigma, lambda). This calculation is implemented in the
#'   non-exported function 'score'. Similarly, the regular expected Fisher
#'   information can be obtained as diag(c(rep(1, length(Beta)), sigma, lambda))
#'   \%*\% loglik(...)[[3]] \%*\% diag(c(rep(1, length(Beta)), sigma, lambda)),
#'   and this calculation is implemented in the non-exported function
#'   'fish_inf'.
#'   
#'   When one or more scale parameter are zero, the corresponding elements
#'   of the usual score function are identically zero regardless of the data,
#'   but the scaled score function, which can be defined as a limit, is still
#'   useful. In particular, it has a positive definite covariance matrix.
#'   }
#'
#' @export
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
  diag(scale_vec, length(scale_vec)) %*% out[[3]] %*% diag(scale_vec,
  length(scale_vec))
}


