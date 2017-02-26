get_hotelling.fd <- function(x, y, mu, paired = FALSE, grid_step = 1) {
  if (!inherits(x, "fd") || !inherits(y, "fd") || !inherits(mu, "fd"))
    stop("Input must be a functional data object of class fd as defined in the
         fda package.")

  rangeval <- x$basis$rangeval
  if (any(y$basis$rangeval != rangeval))
    stop("Inputs have not the same range.")
  if (any(mu$basis$rangeval != rangeval))
    stop("The assumed mean function does not have the same range as the input
         data.")

  dx <- (rangeval[2] - rangeval[1]) * grid_step
  abscissa <- seq(rangeval[1], rangeval[2], by = dx)

  xx <- t(fda::eval.fd(fdobj = x, evalarg = abscissa))
  yy <- t(fda::eval.fd(fdobj = y, evalarg = abscissa))
  mmu <- t(fda::eval.fd(fdobj = mu, evalarg = abscissa))
  get_hotelling(x = xx, y = yy, mu = mmu, paired = paired)
}

test <- function(x, y, mu, paired = FALSE, grid_step = 1) {
  if (!inherits(x, "matrix")) {

  }
}
