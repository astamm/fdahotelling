#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`

#' Power of symmetric positive definite matrix
#'
#' \code{matrix_power} computes the power of a symmetric positive definite
#' matrix using spectral decomposition via the \code{eigen} function of the R
#' \code{base} package in its symmetric version.
#'
#' @param matrix Symmetric positive definite matrix to compute the power from.
#' @param power Desired power.
#' @return A symmetric positive definite matrix with same eigenvectors as input
#'   matrix and eigenvalues being the powered eigenvalues of the input matrix.
#' @keywords internal
matrix_power <- function(matrix, power) {
  eigs <- eigen(matrix, symmetric = T)
  isPD <- TRUE
  for (l in rev(eigs$values)) {
    if (l <= .Machine$double.eps) {
      isPD <- FALSE
      break;
    }
  }
  if (!isPD)
    stop("Input matrix should be positive definite.")
  return(eigs$vectors %*% diag(eigs$values^power) %*% t(eigs$vectors))
}

#' Random index generator for data splitting
#'
#' \code{random_group} randomly splits a vector of indices (defined by the
#' number of lines of a dataset) into subgroups.
#'
#' @param n Total number of indices to be splitted into smaller groups.
#' @param probs A named vector of probabilities defining the names and
#' proportions of indices defining each subgroup.
#' @return A character vector of size \code{n} in which each entry contains the
#' name of one of the subgroups defined in the \code{probs} vector.
#' @keywords internal
random_group <- function(n, probs) {
  probs <- probs / sum(probs)
  g <- findInterval(seq(0, 1, length = n), c(0, cumsum(probs)),
                    rightmost.closed = TRUE)
  names(probs)[sample(g)]
}

split_matrix <- function(x, indices) {
  nn <- levels(as.factor(indices))
  res <- list(x[indices == nn[1], , drop = FALSE],
              x[indices == nn[2], , drop = FALSE])
  names(res) <- nn
  res
}

#' Data partitioning
#'
#' \code{partition} generates random partitions of the input dataset into
#' smaller groups.
#'
#' @param df Input dataset to be partitioned.
#' @param n Number of random partitions.
#' @param probs A named vector of probabilities defining the names and
#' relative sizes of each subsequent subgroup.
#' @return A character vector of size \code{n} in which each entry contains the
#' name of one of the subgroups defined in the \code{probs} vector.
#' @keywords internal
partition <- function(df, n, probs) {
  if (is.null(dim(n))) {
    return(n %>%
             replicate(
               expr = split_matrix(df, random_group(nrow(df), probs)),
               simplify = FALSE
             ) %>%
             purrr::transpose() %>%
             dplyr::as_data_frame()
    )
  }
  dplyr::data_frame(combId = seq_len(dim(n)[2])) %>%
    dplyr::mutate(D1 = purrr::map(combId, ~ df[ n[, .x], ]),
                  D2 = purrr::map(combId, ~ df[-n[, .x], ])) %>%
    dplyr::select(-combId) %>%
    purrr::set_names(names(probs))
}

check_arguments <- function(x, y = NULL, mu = NULL, step_size = 0) {

  if (inherits(x, "fd")) {
    if (step_size == 0)
      stop("You must specify the step size when using objects of class fd.")
    if (!requireNamespace("fda", quietly = TRUE))
      stop("The package fda is required for using this function with inputs of
           type fd. Please install it.", call. = FALSE)
  }


  if (inherits(x, "matrix")) {

    if ((!is.null(y) && !inherits(y, "matrix")) ||
        (!is.null(mu) && !inherits(mu, "numeric")))
      stop("The argument x is a matrix. Hence, y should be a matrix too and mu
           should be a vector.")

    if (is.null(y))
      y <- matrix(nrow = 0, ncol = 0)

    if (is.null(mu))
      mu <- rep(0, dim(x)[2])

  } else if (inherits(x, "data.frame")) {

    if ((!is.null(y) && !inherits(y, "data.frame")) ||
        (!is.null(mu) && !inherits(mu, "numeric")))
      stop("The argument x is a data frame. Hence, y should be a data frame too
           and mu should be a vector.")

    x <- as.matrix(x)

    if (is.null(y))
      y <- matrix(nrow = 0, ncol = 0)
    else
      y <- as.matrix(y)

    if (is.null(mu))
      mu <- rep(0, dim(x)[2])

  } else if (inherits(x, "fd")) {
    if ((!is.null(y) && !inherits(y, "fd")) ||
        (!is.null(mu) && !inherits(mu, "fd")))
      stop("The argument x is an object of class fd. Hence, both y and mu should
           also be objects of class fd.")

    rangeval <- x$basis$rangeval

    if (!is.null(y)) {
      rangeval[1] <- max(rangeval[1], y$basis$rangeval[1])
      rangeval[2] <- min(rangeval[2], y$basis$rangeval[2])
    }

    if (!is.null(mu)) {
      rangeval[1] <- max(rangeval[1], mu$basis$rangeval[1])
      rangeval[2] <- min(rangeval[2], mu$basis$rangeval[2])
    }

    if (rangeval[1] > rangeval[2])
      stop("Common range of abscissa values between inputs is empty.")

    abscissa <- seq(rangeval[1], rangeval[2], by = step_size * (rangeval[2] - rangeval[1]))

    x <- t(fda::eval.fd(fdobj = x, evalarg = abscissa))

    if (is.null(y))
      y <- matrix(nrow = 0, ncol = 0)
    else
      y <- t(fda::eval.fd(fdobj = y, evalarg = abscissa))

    if (is.null(mu)) {
        mu <- fda::fd(mu)
        mu$basis$rangeval <- rangeval
    }

    mu <- t(fda::eval.fd(fdobj = mu, evalarg = abscissa))

  } else {
    stop("Unrecognised type of input data. Please provide the data either in
         matrix format, or data frame format or fd format.")
  }

  list(x = x, y = y, mu = mu)
}
