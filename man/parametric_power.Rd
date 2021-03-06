% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulations.R
\name{parametric_power}
\alias{parametric_power}
\title{Power calculation for the Hotelling's parametric test}
\usage{
parametric_power(mu1 = 0, Sigma = diag(length(mu1)), n1 = 10L,
  n2 = NULL, MC = 1000L, alpha = 0.05, paired = FALSE, step_size = 0)
}
\arguments{
\item{mu1}{True mean function or difference between mean functions
(default: 0).}

\item{Sigma}{True covariance matrix \eqn{\Sigma} (default: 1).}

\item{n1}{Sample size of data pertaining to the 1st population (default: 10).}

\item{n2}{Sample size of data pertaining to the 2nd population (default:
\code{NULL}).}

\item{MC}{Number of Monte-Carlo runs to estimate statistical power (default:
1000).}

\item{alpha}{Significance level (default: 0.05).}

\item{paired}{Is the input data paired? (default: \code{FALSE}).}

\item{step_size}{The step size used to perform integral approximation via the
method of rectangles (default: \code{0}). When set to \code{0}, it assumes
that we are dealing with multivariate data rather than functional data and
thus no integration is necessary.}
}
\value{
An estimate of the statistical power of the test.
}
\description{
\code{parametric_power} computes an estimate of the statistical power of
Hotelling's T^2 parametric test on the mean function (or on the difference
between the mean functions) using Monte-Carlo simulations.
}
\details{
This function computes the statistical power of Hotelling's T^2 parametric
test based on a specific generative model. It is assumed that data is
generated from a Gaussian distribution. The user-defined inputs are
\itemize{
\item \code{mu1}: a numeric vector containing the actual mean function of the
distribution (or difference between mean functions) evaluated in a pointwise
fashion on a uniform grid.
\item \code{Sigma}: a numeric matrix containing the actual covariance kernel
of the distribution (or pooled covariance kernel) evaluated in a pointwise
fashion on a uniform grid.
}
}
\examples{
# Set the sample sizes for a two-sample test
n1 <- 10
n2 <- 10
# Set the dimensionality for curve approximation
p <- 100
# Set the actual covariance kernel
Sigma <- diag(1, p)
# The following lines of code computes the actual significance level of the
# test.
mu1 <- rep(0, p)
parametric_power(mu1, Sigma, n1, n2, MC = 1000)
# We can use more complex covariance matrices and compute an estimate of the
# statistical power of the test for a given non-null mean difference.
mu1 <- rep(4, p)
s <- seq(-1, 1, length.out = 100)
Sigma <- outer(s, s, function(t, s) exp(-abs(t - s)))
parametric_power(mu1, Sigma, n1, n2, MC = 1000)
}
\seealso{
The underlying statistical test is described in details in Secchi,
 P., Stamm, A., & Vantini, S. (2013). \emph{Inference for the mean of large
 \eqn{p} small \eqn{n} data: A finite-sample high-dimensional generalization
 of Hotelling theorem}. Electronic Journal of Statistics, 7, pp. 2005-2031.
 doi:10.1214/13-EJS833, available online at
 \url{http://projecteuclid.org/euclid.ejs/1375708877}.
}
