#' Statistics for functional data
#'
#' Compute various statistics for functional data using numerical integration
#' when needed. A total of 7 statistics are available in the package and the
#' user can even define its own statistic for use within the statistical test
#' functions.
#'
#' @seealso A detailed description of the individual statistics for functional
#'   data provided by the \pkg{fdahotelling} package can be found in the
#'   vignette \emph{Available statistics for functional data}.
#' @param x Dataframe or matrix containing the data collected from 1st
#'   population.
#' @param y Dataframe or matrix containing the data collected from 2nd
#'   population (default: \code{NULL}).
#' @param mu True mean value (or difference between mean values) under the null
#'   hypothesis (default = 0).
#' @param paired Is the input data paired? (default: \code{FALSE}).
#' @param step_size The step size used to perform integral approximation via the
#'   method of rectangles (default: 0). When set to 0 (default), it assumes that
#'   we are dealing with multivariate data rather than functional data and thus
#'   no integration is necessary.
#' @return All \code{get_*} functions return a named numeric vector of size 1
#'   storing the value of the corresponding statistic except \code{get_all} that
#'   outputs a vector of size 7 containing the values of all 7 statistics.
#' @name get-statistic
NULL

#' @rdname get-statistic
#' @export
get_hotelling <- function(x, y = NULL, mu = NULL,
                          paired = FALSE, step_size = 0) {
  l <- check_arguments(x = x, y = y, mu = mu, step_size = step_size)
  get_hotelling_impl(x = l$x, y = l$y, mu = l$mu,
                     paired = paired, step_size = step_size)
}

#' @rdname get-statistic
#' @export
get_L1 <- function(x, y = NULL, mu = NULL,
                   paired = FALSE, step_size = 0) {
  l <- check_arguments(x = x, y = y, mu = mu, step_size = step_size)
  get_L1_impl(x = l$x, y = l$y, mu = l$mu,
              paired = paired, step_size = step_size)
}

#' @rdname get-statistic
#' @export
get_L1_std <- function(x, y = NULL, mu = NULL,
                       paired = FALSE, step_size = 0) {
  l <- check_arguments(x = x, y = y, mu = mu, step_size = step_size)
  get_L1_std_impl(x = l$x, y = l$y, mu = l$mu,
                  paired = paired, step_size = step_size)
}

#' @rdname get-statistic
#' @export
get_L2 <- function(x, y = NULL, mu = NULL,
                   paired = FALSE, step_size = 0) {
  l <- check_arguments(x = x, y = y, mu = mu, step_size = step_size)
  get_L2_impl(x = l$x, y = l$y, mu = l$mu,
              paired = paired, step_size = step_size)
}

#' @rdname get-statistic
#' @export
get_L2_std <- function(x, y = NULL, mu = NULL,
                       paired = FALSE, step_size = 0) {
  l <- check_arguments(x = x, y = y, mu = mu, step_size = step_size)
  get_L2_std_impl(x = l$x, y = l$y, mu = l$mu,
                  paired = paired, step_size = step_size)
}

#' @rdname get-statistic
#' @export
get_Linf <- function(x, y = NULL, mu = NULL,
                     paired = FALSE, step_size = 0) {
  l <- check_arguments(x = x, y = y, mu = mu, step_size = step_size)
  get_Linf_impl(x = l$x, y = l$y, mu = l$mu,
                paired = paired, step_size = step_size)
}

#' @rdname get-statistic
#' @export
get_Linf_std <- function(x, y = NULL, mu = NULL,
                         paired = FALSE, step_size = 0) {
  l <- check_arguments(x = x, y = y, mu = mu, step_size = step_size)
  get_Linf_std_impl(x = l$x, y = l$y, mu = l$mu,
                    paired = paired, step_size = step_size)
}

#' @rdname get-statistic
#' @export
get_all <- function(x, y = NULL, mu = NULL,
                    paired = FALSE, step_size = 0) {
  l <- check_arguments(x = x, y = y, mu = mu, step_size = step_size)
  get_all_impl(x = l$x, y = l$y, mu = l$mu,
               paired = paired, step_size = step_size)
}
