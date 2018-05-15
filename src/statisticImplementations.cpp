#include "statisticImplementations.h"
#include "statisticClass.h"

arma::colvec
  eigenvalues(const arma::mat& x,
              const arma::mat& y)
  {
    fdahotelling::Statistic stat;
    stat.SetInput1(x);

    if (y.n_cols != 0)
      stat.SetInput2(y);

    stat.SetTolerance(0.0);

    return stat.GetEigenValues();
  }

arma::mat
  pseudoinverse(const arma::mat& x,
                const arma::mat& y,
                const double tolerance)
  {
    fdahotelling::Statistic stat;
    stat.SetInput1(x);

    if (y.n_cols != 0)
      stat.SetInput2(y);

    stat.SetTolerance(tolerance);

    return stat.GetPseudoInverse();
  }

Rcpp::NumericVector
  stat_hotelling_impl(const arma::mat& x,
                      const arma::mat& y,
                      const arma::colvec& mu,
                      const bool paired,
                      const double step_size,
                      const bool use_correction,
                      const double tolerance)
  {
    fdahotelling::Statistic stat;
    stat.SetInput1(x);
    stat.SetUsePairedTest(paired);
    stat.SetStepSize(step_size);
    stat.SetUseCorrection(use_correction);
    stat.SetTolerance(tolerance);

    if (y.n_cols != 0)
      stat.SetInput2(y);

    stat.SetDelta0(mu);

    return stat.GetValue(fdahotelling::Statistic::Hotelling);
  }

Rcpp::NumericVector
  stat_L1_impl(const arma::mat& x,
               const arma::mat& y,
               const arma::colvec& mu,
               const bool paired,
               const double step_size)
  {
    fdahotelling::Statistic stat;
    stat.SetInput1(x);
    stat.SetUsePairedTest(paired);
    stat.SetStepSize(step_size);

    if (y.n_cols != 0)
      stat.SetInput2(y);

    stat.SetDelta0(mu);

    return stat.GetValue(fdahotelling::Statistic::L1);
  }

Rcpp::NumericVector
  stat_L2_impl(const arma::mat& x,
               const arma::mat& y,
               const arma::colvec& mu,
               const bool paired,
               const double step_size)
  {
    fdahotelling::Statistic stat;
    stat.SetInput1(x);
    stat.SetUsePairedTest(paired);
    stat.SetStepSize(step_size);

    if (y.n_cols != 0)
      stat.SetInput2(y);

    stat.SetDelta0(mu);

    return stat.GetValue(fdahotelling::Statistic::L2);
  }

Rcpp::NumericVector
  stat_Linf_impl(const arma::mat& x,
                 const arma::mat& y,
                 const arma::colvec& mu,
                 const bool paired,
                 const double step_size)
  {
    fdahotelling::Statistic stat;
    stat.SetInput1(x);
    stat.SetUsePairedTest(paired);
    stat.SetStepSize(step_size);

    if (y.n_cols != 0)
      stat.SetInput2(y);

    stat.SetDelta0(mu);

    return stat.GetValue(fdahotelling::Statistic::Linf);
  }

Rcpp::NumericVector
  stat_L1_std_impl(const arma::mat& x,
                   const arma::mat& y,
                   const arma::colvec& mu,
                   const bool paired,
                   const double step_size)
  {
    fdahotelling::Statistic stat;
    stat.SetInput1(x);
    stat.SetUsePairedTest(paired);
    stat.SetStepSize(step_size);

    if (y.n_cols != 0)
      stat.SetInput2(y);

    stat.SetDelta0(mu);

    return stat.GetValue(fdahotelling::Statistic::StdL1);
  }

Rcpp::NumericVector
  stat_L2_std_impl(const arma::mat& x,
                   const arma::mat& y,
                   const arma::colvec& mu,
                   const bool paired,
                   const double step_size)
  {
    fdahotelling::Statistic stat;
    stat.SetInput1(x);
    stat.SetUsePairedTest(paired);
    stat.SetStepSize(step_size);

    if (y.n_cols != 0)
      stat.SetInput2(y);

    stat.SetDelta0(mu);

    return stat.GetValue(fdahotelling::Statistic::StdL2);
  }

Rcpp::NumericVector
  stat_Linf_std_impl(const arma::mat& x,
                     const arma::mat& y,
                     const arma::colvec& mu,
                     const bool paired,
                     const double step_size)
  {
    fdahotelling::Statistic stat;
    stat.SetInput1(x);
    stat.SetUsePairedTest(paired);
    stat.SetStepSize(step_size);

    if (y.n_cols != 0)
      stat.SetInput2(y);

    stat.SetDelta0(mu);

    return stat.GetValue(fdahotelling::Statistic::StdLinf);
  }

Rcpp::NumericVector
  stat_all_impl(const arma::mat& x,
                const arma::mat& y,
                const arma::colvec& mu,
                const bool paired,
                const double step_size)
  {
    fdahotelling::Statistic stat;
    stat.SetInput1(x);
    stat.SetUsePairedTest(paired);
    stat.SetStepSize(step_size);
    stat.SetUseCorrection(false);
    stat.SetTolerance(0);

    if (y.n_cols != 0)
      stat.SetInput2(y);

    stat.SetDelta0(mu);

    return stat.GetValue(fdahotelling::Statistic::All);
  }
