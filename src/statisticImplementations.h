// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Positive eigenvalues of the sample covariance matrix
//'
//' \code{eigenvalues} computes the positive eigenvalues of the sample covariance matrix using the \code{RcppArmadillo} interface for fast eigen-decomposition via divide and conquer method.
//'
//' @param x Dataframe or matrix containing the data collected from 1st population.
//' @param y Dataframe or matrix containing the data collected from 2nd population (default: \code{NULL}).
//' @return A numeric vector with the positive eigenvalues of the sample covariance matrix or sample pooled covariance matrix when a second dataset is provided.
//' @keywords internal
// [[Rcpp::export]]
arma::colvec
  eigenvalues(const arma::mat& x,
              const arma::mat& y);

//' Pseudoinverse of the sample covariance matrix
//'
//' \code{pseudoinverse} computes the Moore-Penrose inverse of the sample covariance matrix using the \code{RcppArmadillo} interface for fast pseudoinverse calculations via divide and conquer method.
//'
//' @param x Dataframe or matrix containing the data collected from 1st population.
//' @param y Dataframe or matrix containing the data collected from 2nd population (default: \code{NULL}).
//' @param tolerance Positive value \eqn{\tau} such that positive eigenvalues of the sample covariance matrix that are below \eqn{\tau} tr\eqn{\Sigma} are discarded for inverse computation (default: 0).
//' @return A matrix containing the pseudoinverse of the sample covariance matrix.
//' @keywords internal
// [[Rcpp::export]]
arma::mat
  pseudoinverse(const arma::mat& x,
                const arma::mat& y,
                const double tolerance = 0.0);

//' Compute Hotelling's \eqn{T^2} statistic
//'
//' \code{stat_hotelling_impl} computes the generalized Hotelling's \eqn{T^2} statistic using the \code{RcppArmadillo} interface for fast pseudoinverse calculations via divide and conquer method.
//'
//' @inheritParams statistics
//' @param use_correction Flag to turn on and off a statistic adjustment factor that makes it Fisher-distributed under the normality assumption (default: \code{FALSE}).
//' @param tolerance Positive value \eqn{\tau} such that positive eigenvalues of the sample covariance matrix that are below \eqn{\tau} tr\eqn{\Sigma} are discarded for inverse computation (default: 0).
//' @return A named numeric vector of size 1 storing the Hotelling's \eqn{T^2} statistic defined as \deqn{\int_D}.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector
  stat_hotelling_impl(const arma::mat& x,
                      const arma::mat& y,
                      const arma::colvec& mu,
                      const bool paired = false,
                      const double step_size = 0.0,
                      const bool use_correction = false,
                      const double tolerance = 0.0);

//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector
  stat_L1_impl(const arma::mat& x,
               const arma::mat& y,
               const arma::colvec& mu,
               const bool paired = false,
               const double step_size = 0.0);

//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector
  stat_L2_impl(const arma::mat& x,
               const arma::mat& y,
               const arma::colvec& mu,
               const bool paired = false,
               const double step_size = 0.0);

//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector
  stat_Linf_impl(const arma::mat& x,
                 const arma::mat& y,
                 const arma::colvec& mu,
                 const bool paired = false,
                 const double step_size = 0.0);

//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector
  stat_L1_std_impl(const arma::mat& x,
                   const arma::mat& y,
                   const arma::colvec& mu,
                   const bool paired = false,
                   const double step_size = 0.0);

//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector
  stat_L2_std_impl(const arma::mat& x,
                   const arma::mat& y,
                   const arma::colvec& mu,
                   const bool paired = false,
                   const double step_size = 0.0);

//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector
  stat_Linf_std_impl(const arma::mat& x,
                     const arma::mat& y,
                     const arma::colvec& mu,
                     const bool paired = false,
                     const double step_size = 0.0);

//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector
  stat_all_impl(const arma::mat& x,
                const arma::mat& y,
                const arma::colvec& mu,
                const bool paired = false,
                const double step_size = 0.0);
