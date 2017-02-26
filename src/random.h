#include <RcppArmadillo.h>

//' Multivariate Gaussian random generator in C++
//'
//' \code{rmvnorm} generates random samples from a multivariate Gaussian distribution with mean \eqn{\mu} and variance-covariance matrix \eqn{\Sigma} assumed positive semi-definite.
//'
//' @param n Number of samples to be drawn from the distribution (default: 1).
//' @param mean Mean vector of the distribution (default: 0).
//' @param squareRootSigma Square-root of the variance-covariance matrix of the distribution (default = 1).
//' @return A matrix of size \eqn{n\times} \code{length(mean)} containing \eqn{n} randomly sampled vectors of size \code{length(mean)}.
//' @keywords internal
// [[Rcpp::export]]
arma::mat
rmvnorm(const unsigned int n,
        const arma::colvec& mean,
        const arma::mat& squareRootSigma);
