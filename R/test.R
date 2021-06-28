#' Compute the (modified) score test statistic
#'
#' @description{
#'   The statistic is based on a modified score function that, unlike the usual
#'   score function, has a full rank covariance matrix on the whole parameter
#'   space. See ?log_lik for details on the parameterization.
#'   
#'   The test is reliable at and near the boundary of the
#'   parameter space where some scale parameters (sigma, lambda) are near zero.
#' }
#'
#' @param y A vector of observed responses.
#' @param X A matrix of predictors whose i:th row corresponds to the i:th
#'   element in y.
#' @param Z A design matrix for the random effects whose i:th row corresponds
#'   to the i:th element in y.
#' @param Beta A vector of regression coefficients of length ncol(X).
#' @param sigma The standard deviation of the error term.
#' @param lambda A vector of scale parameters (standard deviations) of the
#'   random effects.
#' @param lam_idx A vector of length ncol(Z) whose j:th element indicates
#'   which element of lambda scales the j:th random effect.
#' @param test_idx A vector of integers indicating for which elements of
#'   theta = c(Beta, sigma, lambda) the test statistic is to be computed.
#' @param fix_idx A vector of integeres indicating which elements of
#'   theta are treated as fixed and known
#' @param efficient If TRUE, use efficient Fisher information (Schur complement)
#'  for tested parameters. 
#'
#' @return A vector with test statistic ("chi_sq"), degrees of freedom ("df"),
#' and p-value ("p_val").
#'
#' @details{
#'   The linear mixed model assumed is y = X \%*\% Beta + Z \%*\% u + e,
#'   where the vector of random effects u is from a multivariate normal
#'   distribution with mean zero and diagonal covariance matrix. The j:th
#'   diagonal element of that covariance matrix -- the variance of the j:th
#'   random effect -- is equal to lambda(lam_idx[j])^2.
#'  
#'   The elements of e are independent draws from a normal distribution
#'   with mean 0 and variance sigma^2.
#'   
#'   The test statistic is computed for theta[test_idx], where theta = c(Beta,
#'   sigma, lambda). Elements of theta not restricted under the null hypothesis
#'   should ideally be evaluated at their maximum likelihood estimates under the
#'   null hypothesis.
#'   
#'   When no standard deviation parameters are zero under the
#'   null hypothesis, the test statistic is the usual score test statistic
#'   standardized by expected Fisher information.
#' }
#' @useDynLib lmmstest, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @export
score_test <- function(y, X, Z, Beta, sigma, lambda, lam_idx, test_idx,
                       fix_idx = NULL, efficient = TRUE){
  components <- log_lik(y = y, X = X, Z = Z, Beta = Beta, sigma = sigma,
                        lambda = lambda, lam_idx = lam_idx, diffs = 2)
  theta <- c(Beta, sigma, lambda)
  m <- length(theta)
  stopifnot(is.atomic(test_idx), is.null(dim(test_idx)), all(test_idx %in% 1:m),
            length(test_idx) == length(unique(test_idx)), all(fix_idx %in% 1:m),
            length(fix_idx) == length(unique(fix_idx)),
            all(!(fix_idx %in% test_idx)))
  no_idx <- (1:m)[-unique(c(test_idx, fix_idx))]
  s <- components[[2]][test_idx]
  I <- components[[3]]
  I_block <- I[test_idx, test_idx]
  if(efficient & (length(no_idx) > 0)){
    I_block <- I_block -
      I[test_idx, no_idx] %*% solve(I[no_idx, no_idx], I[no_idx, test_idx])
  }

  test_stat <- sum(s * qr.solve(I_block, s))
  df <- length(test_idx)
  p_val <- stats::pchisq(q = test_stat, df = df, lower.tail = FALSE)
  c("chi_sq" = test_stat, "df" = df, "p_val" = p_val)
}
