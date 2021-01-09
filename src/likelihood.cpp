#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List log_lik_rcpp(arma::vec y, const arma::mat X,
  const arma::mat Z, const arma::vec beta, const double sigma,
  const arma::vec lambda, const arma::uvec lam_idx, const uint diffs)
{
  // Raw Y not used
  y -= X * beta;

  // Define constants
  const uint n = y.n_elem;
  const uint p = X.n_cols;
  const uint q = lam_idx.n_elem;
  const uint dim = std::min(n, q);
  const uint d = lambda.n_elem;
  arma::vec Lam(q);
  for(size_t jj = 0; jj < q; jj ++){
    Lam(jj) = lambda(lam_idx(jj) - 1); // Indexing differ in R and C
  }
  double val = 0;
  arma::vec g(p + d + 1, arma::fill::zeros);
  arma::mat I(p + d + 1, p + d + 1, arma::fill::zeros);

  // Main SVD computation
  arma::mat U(n, dim, arma::fill::zeros);
  arma::mat V(dim, dim, arma::fill::zeros);
  arma::vec D(dim, arma::fill::zeros);
  arma::mat F = Z * arma::diagmat(arma::abs(Lam));
  arma::svd_econ(U, D, V, F);
  // Start computing quantities
  if(n > q){
    // Common quantities
    arma::vec v = 1.0 / (arma::pow(D, -2.0) + std::pow(sigma, -2.0));
    arma::vec w;
    arma::mat SIZ;
    arma::mat ZtSIZ;
    arma::mat ZtSIy;
    if(diffs >= 1){
      w = U * (v % (U.t() * y));
      SIZ = std::pow(sigma, -2.0) * Z;
      SIZ -= std::pow(sigma, -4.0) * (U * arma::diagmat(v) * U.t() * Z);
      ZtSIZ = Z.t() * SIZ;
      ZtSIy = SIZ.t() * y;
    }

    // LOG-LIKELIHOOD
    // log determinant term
    val = -0.5 * (n - q) * std::log(sigma * sigma) -
          0.5 * arma::accu(arma::log(D % D + sigma * sigma));
    // quadratic form using Woodbury identity
    val += -0.5 * std::pow(sigma, -2.0) * arma::accu(arma::square(y));
    val += 0.5 * std::pow(sigma, -4.0) * arma::accu(arma::square(arma::sqrt(v) %
          (U.t() * y )));
    val -= 0.5 * n * std::log(2.0 * arma::datum::pi);
    if(diffs >= 1){
      // SCORE FUNCTION
      // beta score
      g.subvec(0, p - 1) = std::pow(sigma, -2.0) * X.t() * y;
      g.subvec(0, p - 1) -= std::pow(sigma, -4.0) * (X.t() * w);

      // sigma scaled (by 1 / sigma) score
      g(p) = - ((n - q) * std::pow(sigma, -2.0) +
            arma::accu(1.0 / (D % D + std::pow(sigma, 2))));
      g(p) += arma::accu(arma::square(std::pow(sigma, -2.0) * y -
              std::pow(sigma, -4.0) * w));

      // lambda scaled (by 1 / lambda) score
      for(size_t jj = 0; jj < d; jj++){
        arma::uvec idx = arma::find(lam_idx == jj + 1);
        g(p + 1 + jj) = -arma::trace(ZtSIZ.submat(idx, idx));
        g(p + 1 + jj) += arma::accu(arma::square(ZtSIy.elem(idx)));
      }
    }

    // Information
    if(diffs >= 2){
      // Information for beta
      I.submat(0, 0, p - 1, p - 1) = std::pow(sigma, -2.0) * X.t() * X;
      I.submat(0, 0, p - 1, p - 1) -= std::pow(sigma, -4.0) * X.t() * U *
                                       arma::diagmat(v) * U.t() * X;
      // Information for sigma
      I(p, p) = 2.0 * (n - q) * std::pow(sigma, -4.0);
      I(p, p) += 2.0 * arma::accu(arma::square(1.0 / (D % D +
                std::pow(sigma, 2.0))));

      // Information for lambda
      for(size_t jj = 0; jj < d; jj++){
        arma::uvec idx_j = arma::find(lam_idx == jj + 1);
        // Cross with sigma
        I(p + 1 + jj, p) = 2.0 * arma::accu(arma::square(SIZ.cols(idx_j)));
        I(p, p + 1 + jj) = I(p + 1 + jj, p);
        // Within lambda
        for(size_t kk = 0; kk <= jj; kk ++){
          arma::uvec idx_k = arma::find(lam_idx == kk + 1);
          I(p + 1 + jj, p + 1 + kk) =
            2.0 * arma::accu(arma::square(ZtSIZ.submat(idx_j, idx_k)));
          I(p + 1 + kk, p + 1 + jj) = I(p + 1 + jj, p + 1 + kk);
        }
      }
    }
  } else{
    // Not supported yet
  }
  // Force symmetry
  I = 0.5 * (I + I.t());
  return Rcpp::List::create(Rcpp::Named("loglik") = val,
                           Rcpp::Named("scaled_score") = g,
                           Rcpp::Named("scaled_inf") = I);
}
