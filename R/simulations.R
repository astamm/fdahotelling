#' Power calculation for the Hotelling's parametric test
#'
#' \code{parametric_power} computes an estimate of the statistical power of
#' Hotelling's T^2 parametric test on the mean function (or on the difference
#' between the mean functions) using Monte-Carlo simulations.
#'
#' This function computes the statistical power of Hotelling's T^2 parametric
#' test based on a specific generative model. It is assumed that data is
#' generated from a Gaussian distribution. The user-defined inputs are
#' \itemize{
#' \item \code{mu1}: a numeric vector containing the actual mean function of the
#' distribution (or difference between mean functions) evaluated in a pointwise
#' fashion on a uniform grid.
#' \item \code{Sigma}: a numeric matrix containing the actual covariance kernel
#' of the distribution (or pooled covariance kernel) evaluated in a pointwise
#' fashion on a uniform grid.
#' }
#'
#' @inheritParams statistics
#' @param mu1 True mean function or difference between mean functions
#'  (default: 0).
#' @param Sigma True covariance matrix \eqn{\Sigma} (default: 1).
#' @param n1 Sample size of data pertaining to the 1st population (default: 10).
#' @param n2 Sample size of data pertaining to the 2nd population (default:
#' \code{NULL}).
#' @param MC Number of Monte-Carlo runs to estimate statistical power (default:
#'  1000).
#' @param alpha Significance level (default: 0.05).
#' @return An estimate of the statistical power of the test.
#' @seealso The underlying statistical test is described in details in Secchi,
#'  P., Stamm, A., & Vantini, S. (2013). \emph{Inference for the mean of large
#'  \eqn{p} small \eqn{n} data: A finite-sample high-dimensional generalization
#'  of Hotelling theorem}. Electronic Journal of Statistics, 7, pp. 2005-2031.
#'  doi:10.1214/13-EJS833, available online at
#'  \url{http://projecteuclid.org/euclid.ejs/1375708877}.
#' @examples
#' # Set the sample sizes for a two-sample test
#' n1 <- 10
#' n2 <- 10
#' # Set the dimensionality for curve approximation
#' p <- 100
#' # Set the actual covariance kernel
#' Sigma <- diag(1, p)
#' # The following lines of code computes the actual significance level of the
#' # test.
#' mu1 <- rep(0, p)
#' parametric_power(mu1, Sigma, n1, n2, MC = 1000)
#' # We can use more complex covariance matrices and compute an estimate of the
#' # statistical power of the test for a given non-null mean difference.
#' mu1 <- rep(4, p)
#' s <- seq(-1, 1, length.out = 100)
#' Sigma <- outer(s, s, function(t, s) exp(-abs(t - s)))
#' parametric_power(mu1, Sigma, n1, n2, MC = 1000)
#' @export
parametric_power <- function(mu1 = 0, Sigma = diag(length(mu1)), n1 = 10L,
                             n2 = NULL, MC = 1000L, alpha = 0.05,
                             paired = FALSE, step_size = 0) {

  # Checks on dimensionality and sample sizes
  p <- dim(Sigma)[1]
  if (p < 1)
    stop("Problem must have at least dimension 1.")
  if (dim(Sigma)[2] != p)
    stop("Input variance-covariance matrix should be square.")
  if (length(mu1) != p)
    stop("Input mu vector should match covariance dimension.")

  oneSampleTest <- is.null(n2)

  if (paired && !oneSampleTest) {
    if (n2 != n1)
      stop("Sample sizes of x and y datasets must match in a paired scenario.")
    oneSampleTest <- TRUE
  }

  if (oneSampleTest && n1 < 2) # one-sample test
    stop("One-sample inference requires at least two samples.")
  if (!oneSampleTest && (n1 < 1 || n2 < 1 || n1 + n2 < 3)) # two-sample test
    stop("Two-sample inference requires at least one observation per sample and
         at least three observations in total.")

  # Put mu in matrix form
  mu0 <- rep(0, p)
  deltaMatrix <- matrix(mu1, nrow = n1, ncol = p, byrow = T)

  # Square root of covariance matrix
  squareRootSigma <- matrix_power(Sigma, 0.5)

  if (oneSampleTest) n2 <- 1
  q_fisher <- stats::qf(1 - alpha,
                        df1 = min(n1 + n2 - 2, p),
                        df2 = abs(n1 + n2 - 2 - p) + 1)
  df2 <- matrix(nrow = 0, ncol = 0)

  if (oneSampleTest) n2 <- 0
  idx1 <- seq_len(n1)
  idx2 <- seq_len(n2) + n1

  pos <- 0
  for (i in 1:MC) {
    df <- rmvnorm(n = n1 + n2, mean = mu0, squareRootSigma = squareRootSigma)
    df1 <- df[idx1, , drop = FALSE] + deltaMatrix
    if (!oneSampleTest)
      df2 <- df[idx2, , drop = FALSE]
    FS <- stat_hotelling_impl(x = df1, y = df2, mu = 0, step_size = step_size,
                             use_correction = TRUE)
    pos <- pos + (FS > q_fisher)
  }
  pos / MC
}

#' Power calculation for the permutation test
#'
#' \code{permutation_power} computes an estimate of the statistical power of
#' permutation-based test on the mean function (or on the difference between
#' mean functions) using Monte-Carlo simulations.
#'
#' This function computes the statistical power of permutation-based tests using
#' a set of user-specified statistics. The power calculation relies on a
#' specific generative model. It is assumed that data is generated from a
#' Gaussian distribution. The user-defined inputs are
#' \itemize{
#' \item \code{mu1}: a numeric vector containing the actual mean function of the
#' distribution (or difference between mean functions) evaluated in a pointwise
#' fashion on a uniform grid.
#' \item \code{Sigma}: a numeric matrix containing the actual covariance kernel
#' of the distribution (or pooled covariance kernel) evaluated in a pointwise
#' fashion on a uniform grid.
#' }
#'
#' @inheritParams parametric_power
#' @inheritParams tests
#' @param mc.cores Number of cores to run the estimation on (default: 1).
#' @return An estimate of the statistical power of the test.
#' @seealso The underlying statistical test is described in details in the
#'  technical report by Pini, A., Stamm, A., & Vantini, S. (2015).
#'  \emph{Hotelling \eqn{T^2} in functional Hilbert spaces}, available online
#'  at \url{https://mox.polimi.it/publication-results/?id=524&tipo=add_qmox}.
#' @importFrom purrr set_names
#' @examples
#' # Set the sample sizes for a two-sample test
#' n1 <- 10
#' n2 <- 10
#' # Set the dimensionality for curve approximation
#' p <- 100
#' # Set the actual covariance kernel
#' Sigma <- diag(1, p)
#' #----------------------------------------------------------------------------
#' # The following lines of code computes the actual significance level of the
#' # test and runs in about 1 minute on a single core. Computation time is
#' # linear in the number of cores and does not depend on the statistic or the
#' # number of statistics used in the testing procedure.
#' mu1 <- rep(0, p)
#' permutation_power(mu1, Sigma, n1, n2, MC = 100, B = 250)
#' #----------------------------------------------------------------------------
#' # Usually, it is recommended to set MC = 1000 and B = 1000, which would take
#' # about 40 minutes. However, this function is parallelized and can be run on
#' # multiple cores by setting the optional argument mc.cores to a suitable
#' # integer. On a computer with 4 physical cores, the following lines of code
#' # runs in about 10 minutes:
#' \dontrun{
#' mu1 <- rep(0, p)
#' permutation_power(mu1, Sigma, n1, n2, MC = 1000, B = 1000, mc.cores = 4)
#' }
#' #----------------------------------------------------------------------------
#' # We can use more complex covariance matrices and compute an estimate of the
#' # statistical power of the test for a given non-null mean difference.
#' \dontrun{
#' mu1 <- rep(4, p)
#' s <- seq(-1, 1, length.out = 100)
#' Sigma <- outer(s, s, function(t, s) exp(-abs(t - s)))
#' permutation_power(mu1, Sigma, n1, n2, MC = 1000, B = 1000, mc.cores = 4)
#' }
#' @export
permutation_power <- function(mu1 = 0, Sigma = diag(length(mu1)), n1 = 10L,
                              n2 = NULL, MC = 1000L, alpha = 0.05,
                              paired = FALSE, step_size = 0, B = 1000L,
                              statistic = "Hotelling", mc.cores = 1L) {
  # Checks on dimensionality and sample sizes
  p <- dim(Sigma)[1]
  if (p < 1)
    stop("Problem must have at least dimension 1.")
  if (dim(Sigma)[2] != p)
    stop("Input variance-covariance matrix should be square.")
  if (length(mu1) != p)
    stop("Input mu vector should match covariance dimension.")

  oneSampleTest <- is.null(n2)

  if (paired && !oneSampleTest) {
    if (n2 != n1)
      stop("Sample sizes of x and y datasets must match in a paired scenario.")
    oneSampleTest <- TRUE
  }

  if (oneSampleTest && n1 < 2) # one-sample test
    stop("One-sample inference requires at least two samples.")
  if (!oneSampleTest && (n1 < 1 || n2 < 1 || n1 + n2 < 3)) # two-sample test
    stop("Two-sample inference requires at least one observation per sample and
         at least three observations in total.")

  # Warning pval resolution
  if (oneSampleTest)
    n_comb <- 2^n1
  else
    n_comb <- choose(n1 + n2, n1)
  p_min <- 1 / min(B, n_comb)
  writeLines(paste(" - P-value resolution:", p_min))
  if (n_comb < B)
    writeLines(" - Using exact p-value calculations.")
  else
    writeLines(" - Using approximate p-value calculations by Monte-Carlo
               simulations.")
  if (alpha < p_min)
    warning("Significance level cannot be reached given current sample sizes.")

  mu0 <- rep(0, p)
  deltaMatrix <- matrix(mu1, nrow = n1, ncol = p, byrow = T)
  squareRootSigma <- matrix_power(Sigma, 0.5)
  df2 <- matrix(nrow = 0, ncol = 0)
  if (oneSampleTest) n2 <- 0
  idx1 <- seq_len(n1)
  idx2 <- seq_len(n2) + n1

  # Define internal function to perform test with H0: mu = 0
  single_test <- function(i, n1, n2, mu0, squareRootSigma, deltaMatrix, B,
                          statistic, step_size) {
    print(paste("Iteration", i))
    df <- rmvnorm(n = n1 + n2, mean = mu0, squareRootSigma = squareRootSigma)
    df1 <- df[idx1, , drop = FALSE] + deltaMatrix
    if (!oneSampleTest)
      df2 <- df[idx2, , drop = FALSE]
    test <- permutation_test(x = df1, y = df2, mu = mu0, step_size = step_size,
                             B = B, statistic = statistic, verbose = FALSE,
                             skip_check = TRUE)
    test$pValue %>% set_names(test$statName)
  }

  system.time({
    pvalList <- NULL
    if (mc.cores == 1L)
      pvalList <- lapply(seq_len(MC), single_test, n1, n2, mu0, squareRootSigma,
                         deltaMatrix, B, statistic, step_size)
    else {
      cl <- parallel::makeCluster(mc.cores)
      tryCatch({
        pvalList <- parallel::parLapply(cl, seq_len(MC), single_test, n1, n2,
                                        mu0, squareRootSigma, deltaMatrix, B,
                                        statistic, step_size)
      }, finally = {
        parallel::stopCluster(cl)
      })
    }
  }) %>% print()

  pvalList %>%
    purrr::transpose() %>%
    purrr::map(`<`, alpha) %>%
    purrr::map_dbl(mean)
}
