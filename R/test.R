#' Compute (modified) score test statistic in a linear mixed model
#'
#' The statistic is based on a modified score function that, unlike the usual
#' score function, has a full rank covariance matrix on the whole parameter
#' space. In particular, the test is reliable at and near the boundary of the
#' parameter space where some scale parameters (sigma, lambda) are near zero.
#'
#' @param y a vector of observed responses
#' @param X a matrix of predictors whose i:th row corresponds to the i:th
#' element in y
#' @param Z a design matrix for the random effects whose i:th row corresponds
#' to the i:th element in y
#' @param Beta a vector of regression coefficients of length ncol(X)
#' @param sigma the standard deviation of the error term
#' @param lambda a vector of scale parameters (standard deviations) of the
#' random effects
#' @param lam_idx a vector of length ncol(Z) whose j:th element indicates
#' which element of lambda scales the j:th random effect
#' @param test_idx a vector of integers indicating for which elements of
#' theta = c(Beta, sigma, lambda) the test statistic is to be computed
#'
#' @return a vector with test-statistic ("chi_sq"), degrees of freedom ("df"),
#' and p-value ("p_val")
#'
#' @details The modified score replaces the first partial derivative of the
#'   log-likelihood with respect to a standard deviation parameter by the second
#'   partial derivative at points where the former is almost surely equal to
#'   zero, which happens if and only if the standard deviation parameter is
#'   equal to zero. When no standard deviation parameters are zero under the
#'   null hypothesis, the test statistic is the usual score test statistic
#'   standardized by expected Fisher information.
#' @useDynLib lmmstest, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @export
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
  p_val <- stats::pchisq(q = test_stat, df = df, lower.tail = FALSE)
  c("chi_sq" = test_stat, "df" = df, "p_val" = p_val)
}
